import os
import re
import numpy as np
import pandas as pd
import pysam
from genepeeks.common import utilities as util


class CoverageMatrix(object):
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

    def __init__(self, unwanted_filters=None, min_interval_separation=629):
        super(CoverageMatrix, self).__init__()
        self.logger = util.create_logging()
        self.list_of_checks = self.filter_list_of_checks(unwanted_filters) if unwanted_filters else self.default_checks
        self.min_interval_separation = min_interval_separation

    def filter_list_of_checks(self, unwanted_filters):
        """ Remove any unwanted filters from the list of checks to perform on each read """

        # Ensure that the provided unwanted_filters are real filters that can be removed
        if not isinstance(unwanted_filters, (list, tuple)):
            util.stop_err('unwanted_filters must be a list or tuple, the following is invalid: {}'.format(unwanted_filters))
        check_names = [check_name for check, check_name in self.default_checks]
        for unwanted_filter in unwanted_filters:
            if unwanted_filter in check_names:
                self.logger.info('Removing {} from list_of_checks to perform on each read'.format(unwanted_filter))
            else:
                util.stop_err('{} is not a valid check_name to remove from the list_of_checks'.format(unwanted_filter))

        # Filter the list of checks to remove the unwanted filters
        return filter(lambda x: x[1] not in unwanted_filters, self.default_checks)

    def passes_checks(self, read, insert_length, skipped_counts):
        """ Only count reads that pass the necessary quality checks, and keep counts of those that don't """
        for check, check_name in self.list_of_checks:
            if check(read, insert_length, self.min_interval_separation):
                if skipped_counts is not None:
                    util.add_to_dict(skipped_counts, check_name)
                return False
        return True

    def get_unique_panel_intervals(self, print_counts=True, panel_chrom='chrX'):
        """ Get the intervals that are unique to each panel """
        panel_intervals = {
            'TSID': {'file': os.path.join(os.path.realpath(os.path.dirname(__file__)), 'inputs', 'TruSight_Inherited_Disease_Manifest_A.bed')},
            'TSO': {'file': os.path.join(os.path.realpath(os.path.dirname(__file__)), 'inputs', 'TruSight-One-BED-May-2014.txt')}
        }

        for name, intrv_info in panel_intervals.items():
            df = pd.read_csv(intrv_info['file'], delimiter='\t', header=None, names=('chrom', 'start', 'end', 'id'))
            filtered_df = df[df['chrom'] == panel_chrom].sort_values(by='start')
            intrv_info['intrv_list'] = map(dict, dict(filtered_df.T).values())

        # For each panel, get the intervals that are in that panel but missing from the other panel
        unique_panel_intervals = {
            'TSID': util.interval_diff(panel_intervals['TSO']['intrv_list'], panel_intervals['TSID']['intrv_list'], extend_by=self.min_interval_separation),
            'TSO': util.interval_diff(panel_intervals['TSID']['intrv_list'], panel_intervals['TSO']['intrv_list'], extend_by=self.min_interval_separation)
        }

        # Merge each panel's unique intervals and add a label to each interval
        for panel, unique_intervals in unique_panel_intervals.items():
            if print_counts:
                self.logger.info('{} only: {} intervals over {} bp'.format(panel, len(unique_intervals), util.add_intervals(unique_intervals)))
            unique_intervals_merged = util.merge_intervals(unique_intervals, min_dist=self.min_interval_separation, print_counts=print_counts)
            for i, intrv in enumerate(unique_intervals_merged):
                intrv['chrom'] = panel_chrom
                intrv['label'] = 'unique_{}_Target_{}'.format(panel, i + 1)
            unique_panel_intervals[panel] = unique_intervals_merged
        return unique_panel_intervals

    @staticmethod
    def get_sample_info(RG):
        """ Gather identifying info for each sample """
        subject = RG.get('SM')

        # Get the sex for GP subjects as the first letter in the subject ID
        sex = subject[0] if subject and subject[0] in ('M', 'F') else None

        specimen = RG.get('LB')
        if specimen:
            # Get the baits for GP subjects, which indicate if the sample was sequenced under TruSight One (T), Inherited Disease (P), or mixed (M)
            baits_results = re.findall(r'([MTP])\d$', specimen)
            baits = baits_results[0] if baits_results else None
        else:
            baits = None

        # Get the flow_cell_id, removing the lane if provided
        flow_cell_id = RG.get('PU')
        if flow_cell_id and flow_cell_id[-5:-2] == '-L0':
            flow_cell_id, __ = flow_cell_id.rsplit('-', 1)

        sample_info = [subject, specimen, sex, baits, flow_cell_id]
        return sample_info

    def get_subject_info(self, bamfile):
        """ Gather identifying info for each subject, checking for differences between samples """
        subject_info = []
        for RG in bamfile.header.get('RG', []):
            sample_info = self.get_sample_info(RG)
            if not subject_info:
                subject_info = sample_info
            elif subject_info != sample_info:
                # Check if a subject has multiple different values for any of the headers
                for i, (subject_value, sample_value) in enumerate(zip(subject_info, sample_info)):
                    if subject_value != sample_value:
                        if subject_value is None:
                            # If the field is empty in the subject_info, use the value from the sample_value
                            subject_info[i] = sample_value
                        elif sample_value is not None and sample_value not in subject_value.split('|'):
                            # If the field has data in both subject_info and sample_info, combine the values together
                            subject_info[i] += '|{}'.format(sample_value)
        return subject_info

    def get_unique_panel_reads(self, bamfile_path, unique_panel_intervals, subject_id):
        """ Count the reads that fall in intervals anywhere in the X chromosome that are unique to each panel """
        unique_panel_reads = {}
        for panel, unique_intervals in unique_panel_intervals.iteritems():
            aligned_bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
            coverage_vector = self.get_subject_coverage(aligned_bamfile, unique_intervals, subject_id)
            unique_panel_reads[panel] = sum(coverage_vector)
        return unique_panel_reads

    def get_subject_coverage(self, bamfile, targets, subject_id, skipped_counts=None):
        """ Get vector of coverage counts for any given bamfile across any provided target regions """

        coverage_vector = []
        for target in targets:
            read_pairs = {}
            # Scan through all reads in each target, and count the number of unique read pairs
            reference = re.sub(r'^chr', '', target['chrom'])
            for read in bamfile.fetch(reference=reference, start=target['start'], end=target['end']):
                insert_length = read.template_length
                # Multiply insert_length by -1 if the read is reverse
                if read.is_reverse:
                    insert_length *= -1

                # Only count reads that pass the necessary quality checks, and keep counts of those that don't
                if self.passes_checks(read, insert_length, skipped_counts):
                    # Keep track of each read pair, and count coverage at the end in order to only count each read pair once
                    pair_start = min(read.reference_start, read.next_reference_start)
                    # util.add_to_dict(read_pairs, read.query_name, nested_key=(pair_start, insert_length))
                    util.add_to_dict(read_pairs, (read.query_name, pair_start, insert_length))

            duplicate_read_pairs = {key: value for key, value in read_pairs.items() if value > 2}
            if duplicate_read_pairs:
                self.logger.warning('For {}, the following read_pairs appeared more than twice within {}: {}'.format(
                    subject_id, target['label'], duplicate_read_pairs))

            # Count the number of unique read pairs as the amount of coverage for any target
            target_coverage = len(read_pairs)
            coverage_vector.append(target_coverage)

        return coverage_vector

    def create_coverage_matrix(self, bamfiles_fofn, targets):
        """ Create coverage matrix with exons as columns, samples as rows, and amount of coverage in each exon as the values,
        plus extra columns for identifying info for each sample """

        unique_panel_intervals = self.get_unique_panel_intervals()

        # Initiate matrix headers
        base_headers = ['subject', 'specimen', 'sex', 'baits', 'flow_cell_id']
        extra_headers = ['bwa_version', 'date_modified', 'TSID_only', 'TSO_only']
        full_headers = base_headers + extra_headers + [target['label'] for target in targets]

        skipped_counts = {}
        coverage_matrix = []

        if isinstance(bamfiles_fofn, list):
            bamfile_paths = bamfiles_fofn
        else:
            # If the bamfile paths are not already a list, create the list from the provided fofn (file of file names)
            with open(bamfiles_fofn) as f:
                bamfile_paths = [bamfile_path.strip() for bamfile_path in f.readlines()]
        file_count = len(bamfile_paths)

        starting_message = '\nCreating coverage_matrix with {} subjects'.format(file_count)
        timing_fields = util.initiate_timer(message=starting_message, add_counts=True, logger=self.logger,
                                            total_counts=file_count, count_steps=15 if file_count > 100 else None)

        # Iterate over all the provided bamfile paths and create the coverage_matrix
        for bamfile_path in bamfile_paths:
            if not os.path.exists(bamfile_path):
                util.stop_err('The bamfile path {} does not exist'.format(bamfile_path))

            bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
            if not bamfile.has_index():
                util.stop_err('{} is missing an index'.format(bamfile_path))

            # Get identifying info for each subject
            subject_info = self.get_subject_info(bamfile)
            if not subject_info:
                subject_info = [None] * len(base_headers)
            if not subject_info[0]:
                bamfile_name = os.path.basename(bamfile_path)
                subject_info[0] = bamfile_name
                self.logger.warning('Missing subject_id for {} and using bamfile name instead'.format(bamfile_name))
            subject_id = subject_info[0]
            if len(subject_info) != len(base_headers):
                util.stop_err('Unequal number of values in subject_info vs base_headers: {} vs {}'.format(base_headers, subject_info))

            date_modified = os.path.getmtime(bamfile_path)
            bwa_version = next((PG.get('VN') for PG in bamfile.header.get('PG', []) if PG.get('ID') == 'bwa'), None)

            # Get subject coverage vector
            subj_coverage_vector = self.get_subject_coverage(bamfile, targets, subject_id, skipped_counts=skipped_counts)
            if subj_coverage_vector.count(0) * 2 > len(targets):
                self.logger.warning('{} is missing coverage for more than half of its targets'.format(subject_id))

            # Get counts of reads in unique panel regions
            unique_panel_reads = self.get_unique_panel_reads(bamfile_path, unique_panel_intervals, subject_id)
            if unique_panel_reads['TSID'] == unique_panel_reads['TSO'] == 0:
                self.logger.warning('{} does not have any coverage that is unique to either TSO or TSID'.format(subject_id))

            # Create subject row with all needed info for the subject, and add to the coverage_matrix
            extra_data = [bwa_version, date_modified, unique_panel_reads['TSID'], unique_panel_reads['TSO']]
            full_subj_vector = subject_info + extra_data + subj_coverage_vector
            if len(full_subj_vector) != len(full_headers):
                util.stop_err('Unequal number of columns ({}) vs headers ({})'.format(len(full_subj_vector), len(full_headers)))
            coverage_matrix.append(full_subj_vector)

            util.get_timing(timing_fields, display_counts=True)

        coverage_df = pd.DataFrame(coverage_matrix, columns=full_headers)

        # Add a column with the ratio of inherited disease only reads compared to Trusight One only reads
        coverage_df['TSID_ratio'] = coverage_df['TSID_only'] / (coverage_df['TSID_only'] + coverage_df['TSO_only'])

        # Add a column for sum of all baseline counts, if any baseline targets exist
        baseline_columns = [column for column in coverage_df.columns if column.startswith('Baseline')]
        if baseline_columns:
            coverage_df['BaselineSum'] = np.sum(coverage_df[baseline_columns], axis=1)

        # Log counts of skipped reads
        for key, count in skipped_counts.items():
            self.logger.info('{} reads were skipped from: {}'.format(count, key))
        return coverage_df
