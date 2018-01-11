import datetime
import logging
import os
import re
from collections import Counter

import numpy as np
import pandas as pd
import pysam

from cnv.Targets.TargetCollection import TargetCollection


class ReadGroups(object):
    """ This is a class to hold all the ReadGroup (RG) tags in a BAM header file.  We use it
     to verify that if the file follows the typical GATK conventions
     (e.g. https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)
      we can grab the SM, LB and PU tags, and verify that only one LB tag is present. """
    __slots__ = ("sample", "library", "flowcell", "rgs")

    def __init__(self, bam_file):
        """
        Create a new set of readgroups, formed from a list of readgroup dictionaries returned by
        pysam.AlignmentFile.header.get('RG', [])

        :param bam_file: A pysam AlignmentFile
        """
        self.rgs = bam_file.header.get('RG', [])
        self._set_item_or_raise_error("sample", "SM", bam_file.filename, bam_file.filename)
        self._set_item_or_raise_error("library", "LB", bam_file.filename)
        self._set_item_or_raise_error("flowcell", "PU", bam_file.filename)

    def _set_item_or_raise_error(self, attr, tag, bam_file_name, default = "NA"):
        """ There might be several read groups, but they should all have the same SM and LB tags, this method verifies this.
        and sets the attribute equal to this shared value.  In the case of the PU tag, which should hold the flow cell
        information, we clip off the lane information.  e.g. HM2NMAFXX-L004 -> HM2NMAFXX-L004"""
        tag_vals = [x.get(tag) for x in self.rgs if x.has_key(tag)]
        if tag == "PU":
            # Clip lane number at end of e.g. HM2NMAFXX-L004 to get flow cell info
            tag_vals = [re.sub("\-L[0-9]{3}", "", k) for k in tag_vals]
        tag_vals = set(tag_vals)
        if not tag_vals:
            logging.info('File: {} missing RG info for tag  {}.  Using value {}  instead'.format(bam_file_name,
                                                                                                 tag,
                                                                                                 default))
            setattr(self, attr, default)
        elif len(tag_vals) == 1:
            setattr(self, attr, tag_vals.pop())
        elif tag == 'SM':
            raise IOError("File {} had two or more SM tags.  The CNV package expects only one sample per BAM.".format(bam_file_name))
        else:
            logging.warning("File {} had two or more {} tags. Recording all library-preps/flow-cells.".format(bam_file_name, tag))
            setattr(self, attr, '|'.join(tag_vals))

class WrappedBAM(object):
    """ A read-only BAM file with some special methods to get header information useful when creating
     CoverageMatrices """
    def __init__(self, bam_name):
        """
        Opens a new BAM for reading
        :param bam_name: File name
        """
        if not os.path.exists(bam_name):
            raise IOError("File: " + bam_name + " does not exist")
        self.fname = bam_name
        self.file = pysam.AlignmentFile(self.fname, 'rb')
        if not self.file.has_index():
            raise IOError('{} is missing an index'.format(bam_name))
        self.read_groups = ReadGroups(self.file)

    def get_date_modified(self):
        """ Returns the last modified time of the file """
        mtime = os.path.getmtime(self.fname)
        return datetime.datetime.fromtimestamp(mtime)

    def get_bwa_version(self):
        """Look for the BWA version in the headers"""
        return next((PG.get('VN') for PG in self.file.header.get('PG', []) if PG.get('ID') == 'bwa'), None)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()


class CoverageMatrix(object):
    # Collection of filters applied to reads
    default_checks = [
        (lambda read, insert, max_insert: read.is_unmapped, 'unmapped'),
        (lambda read, insert, max_insert: read.mapping_quality != 60, 'MAPQ_below_60'),
        (lambda read, insert, max_insert: read.is_duplicate, 'PCR_duplicate'),
        (lambda read, insert, max_insert: read.mate_is_unmapped, 'mate_is_unmapped'),
        (lambda read, insert, max_insert: not read.is_proper_pair, 'not_proper_pair'),
        (lambda read, insert, max_insert: read.is_reverse == read.mate_is_reverse, 'tandem_pair'),
        (lambda read, insert, max_insert: insert <= 0, 'negative_insert_length'),
        # To ensure that no read pairs overlap multiple targets, skip all reads with
        # insert length greater than the distance used to merge intervals
        (lambda read, insert, max_insert: insert >= max_insert, 'insert_length_greater_than_merge_distance'),
        (lambda read, insert, max_insert: min(read.reference_start, read.next_reference_start) + insert < read.reference_end, 'pair_end_less_than_reference_end')
    ]

    def __init__(self, unwanted_filters=None):
        self.logger = logging.getLogger(__name__)
        self.list_of_checks = self.filter_list_of_checks(unwanted_filters) if unwanted_filters else self.default_checks


    def filter_list_of_checks(self, unwanted_filters):
        """ Remove any unwanted filters from the list of checks to perform on each read """

        # Ensure that the provided unwanted_filters are real filters that can be removed
        if not isinstance(unwanted_filters, (list, tuple)):
            raise TypeError('argument unwanted_filters must be a list or tuple, the following is invalid: {}'.format(unwanted_filters))
        check_names = [check_name for check, check_name in self.default_checks]
        for unwanted_filter in unwanted_filters:
            if unwanted_filter in check_names:
                self.logger.info('Removing {} from list_of_checks to perform on each read'.format(unwanted_filter))
            else:
                raise RuntimeError('{} is not a valid check_name to remove from the list_of_checks'.format(unwanted_filter))

        # Filter the list of checks to remove the unwanted filters
        return filter(lambda x: x[1] not in unwanted_filters, self.default_checks)

    def passes_checks(self, read, insert_length, min_interval_separation, skipped_counts):
        """ Only count reads that pass the necessary quality checks, and keep counts of those that don't """
        for check, check_name in self.list_of_checks:
            if check(read, insert_length, min_interval_separation):
                if skipped_counts is not None:
                    skipped_counts.update([check_name])
                return False
        return True

    def get_subject_coverage(self, bamfile, targets, skipped_counts=None):
        """ Get vector of coverage counts for any given bamfile across any provided target regions """
        coverage_vector = []
        assert isinstance(targets, TargetCollection)
        for target in targets:
            read_pairs = {} # Counting dictionary
            # Scan through all reads in each target, and count the number of unique read pairs
            for read in bamfile.file.fetch(reference=target.chrom, start=target.start, end=target.end):
                insert_length = read.template_length
                # Multiply insert_length by -1 if the read is reverse
                if read.is_reverse:
                    insert_length *= -1
                # Only count reads that pass the necessary quality checks, and keep counts of those that don't
                if self.passes_checks(read, insert_length, targets.min_dist, skipped_counts):
                    # Keep track of each read pair, and count coverage at the end in order to only count each read pair once
                    pair_start = min(read.reference_start, read.next_reference_start)
                    tpl = (read.query_name, pair_start, insert_length)
                    read_pairs[tpl] = read_pairs.get(tpl, 0) + 1

            duplicate_read_pairs = {key: value for key, value in read_pairs.items() if value > 2}
            if duplicate_read_pairs:
                self.logger.warning('For {}, the following read_pairs appeared more than twice within {}: {}'.format(
                    bamfile.filename, target.label, duplicate_read_pairs))

            # Count the number of unique read pairs as the amount of coverage for any target
            target_coverage = len(read_pairs)
            coverage_vector.append(target_coverage)

        return coverage_vector

    def create_coverage_matrix(self, bamfiles_fofn, targets):
        """  Create coverage matrix with exons as columns, samples as rows, and amount of coverage in each exon as the values,
        plus extra columns for identifying info for each sample.

        :param bamfiles_fofn: Either a list of files names or the name of one file containing a BAM file name on each line
        :param targets: A TargetCollection object
        :return: a pandas data frame with the coverage data
        """

        # Initiate matrix headers
        headers = ['sample', 'library', 'flow_cell_id', 'bwa_version', 'date_modified']
        headers += [target.label for target in targets]

        skipped_counts = Counter()

        # Load list of BAM files
        if isinstance(bamfiles_fofn, list):
            bamfile_paths = bamfiles_fofn
        else:
            # If the bamfile paths are not already a list, create the list from the provided fofn (file of file names)
            with open(bamfiles_fofn) as f:
                bamfile_paths = [bamfile_path.strip() for bamfile_path in f.readlines()]
        file_count = len(bamfile_paths)

        # Iterate over all the provided bamfile paths and create the coverage_matrix
        logging.info('\nCreating coverage_matrix with {} files'.format(file_count))
        coverage_matrix = []
        for bamfile_path in bamfile_paths:
            logging.info('Getting coverage for {}'.format(bamfile_path))
            with WrappedBAM(bamfile_path) as bamfile:
                # collect meta-data
                bam_info = [bamfile.read_groups.sample,
                            bamfile.read_groups.library,
                            bamfile.read_groups.flowcell,
                            bamfile.get_bwa_version(),
                            bamfile.get_date_modified()]

                # Get subject coverage vector
                subj_coverage_vector = self.get_subject_coverage(bamfile, targets, skipped_counts=skipped_counts)
                if subj_coverage_vector.count(0) * 2 > len(targets):
                    self.logger.warning('{} is missing coverage for more than half of its targets'.format(bamfile.read_groups.sample))

                bam_info += subj_coverage_vector

                if len(bam_info) != len(headers):
                    raise RuntimeError('Unequal number of columns ({}) vs headers ({})'.format(len(bam_info), len(headers)))
                coverage_matrix.append(bam_info)

        coverage_df = pd.DataFrame(coverage_matrix, columns=headers)

        # Add a column for sum of all baseline counts, if any baseline targets exist
        # Merged baselines must also contain 'Baseline'
        baseline_columns = [column for column in coverage_df.columns if 'Baseline' in column]
        if baseline_columns:
            coverage_df['BaselineSum'] = np.sum(coverage_df[baseline_columns], axis=1)

        # Log counts of skipped reads
        for key, count in skipped_counts.iteritems():
            self.logger.info('{} reads were skipped from: {}'.format(count, key))
        return coverage_df
