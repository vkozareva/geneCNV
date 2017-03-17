import os
import pysam
import pandas as pd
from mando import command, main

from genepeeks.common import utilities as util
import utilities as cnv_util

""" Only count reads that pass all of the following checks """
list_of_checks = [
    (lambda read, insert, max_insert: read.is_unmapped, 'unmapped'),
    (lambda read, insert, max_insert: read.mapping_quality != 60, 'MAPQ below 60'),
    (lambda read, insert, max_insert: read.is_duplicate, 'PCR_duplicate'),
    (lambda read, insert, max_insert: read.mate_is_unmapped, 'mate_is_unmapped'),
    (lambda read, insert, max_insert: not read.is_proper_pair, 'not a proper pair'),
    (lambda read, insert, max_insert: read.is_reverse == read.mate_is_reverse, 'tandem_pair'),
    (lambda read, insert, max_insert: insert <= 0, 'negative insert_length'),
    # To ensure that no read pairs overlap multiple targets, skip all reads with
    # insert length greater than the distance used to merge intervals
    (lambda read, insert, max_insert: insert >= max_insert, 'insert_length greater than interval merge distance'),
    (lambda read, insert, max_insert: min(read.reference_start, read.next_reference_start) + insert < read.reference_end, 'pair_end is less than reference_end'),
]


def passes_checks(read, insert_length, max_insert_length, skipped_counts):
    """ Only count reads that pass the necessary quality checks, and keep counts of those that don't """
    for check, check_name in list_of_checks:
        if check(read, insert_length, max_insert_length):
            if skipped_counts is not None:
                util.add_to_dict(skipped_counts, check_name)
            return False
    return True


class CoverageMatrix(object):
    """docstring for CoverageMatrix

    :param min_interval_separation: Any two intervals that are closer than this distance will be merged together,
        and any read pairs with insert lengths greater than this distance will be skipped. The default value of 629
        was derived to be one less than the separation between intervals for Exon 69 and Exon 70 of DMD.

    """
    def __init__(self, min_interval_separation=629):
        super(CoverageMatrix, self).__init__()
        self.logger = util.create_logging()
        self.min_interval_separation = min_interval_separation

    def get_unique_panel_intervals(self):
        """ Get the intervals that are unique to each panel """
        panel_intervals = {
            'TSID': {'file': os.path.join(os.path.realpath(os.path.dirname(__file__)), 'inputs', 'TruSight_Inherited_Disease_Manifest_A.bed')},
            'TSO': {'file': os.path.join(os.path.realpath(os.path.dirname(__file__)), 'inputs', 'TruSight-One-BED-May-2014.txt')}
        }

        for name, intrv_info in panel_intervals.items():
            df = pd.read_csv(intrv_info['file'], delimiter='\t', header=None, names=('chrom', 'start', 'end', 'id'))
            X_df = df[df['chrom'] == 'chrX'].sort_values(by='start')
            intrv_info['intrv_list'] = map(dict, dict(X_df.T).values())

        # For each panel, get the intervals that are in that panel but missing from the other panel
        self.unique_panel_intervals = {
            'TSID': util.interval_diff(panel_intervals['TSO']['intrv_list'], panel_intervals['TSID']['intrv_list'], extend_by=self.min_interval_separation),
            'TSO': util.interval_diff(panel_intervals['TSID']['intrv_list'], panel_intervals['TSO']['intrv_list'], extend_by=self.min_interval_separation)
        }

        # Merge each panel's unique intervals and add a label to each interval
        for panel, unique_intervals in self.unique_panel_intervals.items():
            unique_intervals_merged = util.merge_intervals(unique_intervals, min_dist=self.min_interval_separation, print_counts=True)
            for i, intrv in enumerate(unique_intervals_merged):
                intrv['label'] = 'unique_{}_Target_{}'.format(panel, i + 1)
            self.logger.info('{} only: {} intervals over {} bp'.format(panel, len(unique_intervals_merged), util.add_intervals(unique_intervals_merged)))
            self.unique_panel_intervals[panel] = unique_intervals_merged
        assert util.interval_intersect(*self.unique_panel_intervals.values()) == []

    def get_sample_info(self, RG, bwa_version, date_modified):
        """ Gather identifying info for each sample """
        try:
            # normal RG['ID'] format: FCLR-GP01-2121_1-M1-1_HGGF5AFXX-L004
            subject, specimen_sample, flow_cell_lane = RG['ID'].split('_')
        except ValueError:
            # older RG['ID'] format: FPWB-0000-429L_1-P1-1
            subject, specimen_sample = RG['ID'].split('_')
            flow_cell_id = None
        else:
            flow_cell_id, __ = flow_cell_lane.rsplit('-', 1)

        # simulated CNV subjects have one of these suffixes in this field
        if 'del' in RG['SM'] or 'dup' in RG['SM']:
            subject = RG['SM']

        gender = subject[0]
        if specimen_sample.startswith(('ACGT', 'Omega')):
            lab, specimen_num, baits, sample = specimen_sample.split('.')
            specimen_num = '{}_{}'.format(lab, specimen_num)
        else:
            specimen_num, baits, sample = specimen_sample.split('-')
        specimen = '{}_{}'.format(subject, specimen_num)
        sample = '{}_{}'.format(subject, specimen_sample)
        full_id = RG['ID']

        sample_info = [full_id, subject, specimen, sample, gender, baits[0], flow_cell_id, bwa_version, date_modified]
        return sample_info

    def get_unique_panel_reads(self, bamfile_path):
        """ Count the reads that fall in intervals anywhere in the X chromosome that are unique to each panel """
        unique_panel_reads = {}
        for panel, unique_intervals in self.unique_panel_intervals.iteritems():
            aligned_bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
            coverage_vector = self.get_subject_coverage(aligned_bamfile, unique_intervals, max_insert_length=self.min_interval_separation, allow_zero_coverage=True)
            unique_panel_reads[panel] = sum(coverage_vector)
        return unique_panel_reads

    @staticmethod
    def get_subject_coverage(bamfile, targets, skipped_counts=None, max_insert_length=629, allow_zero_coverage=False):
        """ Get vector of coverage counts for any given bamfile across any provided target regions """

        coverage_vector = []
        for target in targets:
            read_pairs = {}
            # Scan through all reads in each target, and count the number of unique read pairs
            for read in bamfile.fetch('X', start=target['start'], end=target['end']):
                insert_length = read.template_length
                # Multiply insert_length by -1 if the read is reverse
                if read.is_reverse:
                    insert_length *= -1

                # Only count reads that pass the necessary quality checks, and keep counts of those that don't
                if passes_checks(read, insert_length, max_insert_length, skipped_counts):
                    # Keep track of each read pair, and count coverage at the end in order to only count each read pair once
                    pair_start = min(read.reference_start, read.next_reference_start)
                    util.add_to_dict(read_pairs, (read.query_name, pair_start, insert_length))

            duplicate_read_pairs = {key: value for key, value in read_pairs.items() if value > 2}
            if duplicate_read_pairs:
                util.stop_err('The following read_pairs appeared more than twice within {}: {}'.format(target['label'], duplicate_read_pairs))

            # Count the number of unique read pairs as the amount of coverage for any target
            target_coverage = len(read_pairs)

            # Unless allow_zero_coverage is set to true, do not include subjects that have zero coverage for any of the targets
            if target_coverage or allow_zero_coverage:
                coverage_vector.append(target_coverage)
            else:
                return None

        return coverage_vector

    def create_coverage_matrix(self, bamfiles_fofn, targets):
        """ Create coverage matrix with exons as columns, samples as rows, and amount of coverage in each exon as the values,
        plus extra columns for identifying info for each sample """

        self.get_unique_panel_intervals()

        # Initiate matrix headers
        base_headers = [
            'id', 'subject', 'specimen', 'sample', 'gender', 'baits', 'flow_cell_id',
            'bwa_version', 'date_modified', 'TSID_only', 'TSO_only']
        full_headers = base_headers + [target['label'] for target in targets]

        skipped_counts = {}
        coverage_matrix = []

        # Load the bamfile paths from the provided fofn (file of file names), or use example.bam
        with open(bamfiles_fofn) as f:
            bamfile_paths = [bamfile_path.strip() for bamfile_path in f.readlines()]
        file_count = len(bamfile_paths)

        starting_message = '\nCreating coverage_matrix with {} subjects'.format(file_count)
        timing_fields = util.initiate_timer(message=starting_message, add_counts=True, logger=self.logger,
                                            total_counts=file_count, count_steps=3 if file_count > 100 else None)

        # Iterate over all the provided bamfile paths and create the coverage_matrix
        for bamfile_path in bamfile_paths:
            if not os.path.exists(bamfile_path):
                util.stop_err('The bamfile path {} does not exist'.format(bamfile_path))
            bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
            if bamfile.has_index():
                date_modified = os.path.getmtime(bamfile_path)

                # Get identifying sample info
                bwa_version = next(PG['VN'] for PG in bamfile.header['PG'] if PG.get('ID') == 'bwa')
                sample_info = self.get_sample_info(bamfile.header['RG'][0], bwa_version, date_modified)

                # Get subject coverage info
                subj_coverage_vector = self.get_subject_coverage(bamfile, targets, skipped_counts=skipped_counts)
                if subj_coverage_vector is not None:
                    unique_panel_reads = self.get_unique_panel_reads(bamfile_path)
                    full_subj_vector = sample_info + [unique_panel_reads['TSID'], unique_panel_reads['TSO']] + subj_coverage_vector
                    if len(full_subj_vector) != len(full_headers):
                        util.stop_err('Unequal number of columns ({}) vs headers ({})'.format(len(full_subj_vector), len(full_headers)))

                    coverage_matrix.append(full_subj_vector)
                else:
                    self.logger.warning('{} does not have coverage for every target, and is thus not being included in the coverage_matrix'.format(sample_info[1]))
            else:
                self.logger.warning('{} is missing an index'.format(bamfile_path))

            util.get_timing(timing_fields, display_counts=True)

        coverage_df = pd.DataFrame(coverage_matrix, columns=full_headers)

        # Add a column with the ratio of inherited disease only reads compared to Trusight One only reads
        coverage_df['TSID_ratio'] = coverage_df['TSID_only'] / (coverage_df['TSID_only'] + coverage_df['TSO_only'])

        # Log counts of skipped reads
        for key, count in skipped_counts.items():
            self.logger.info('{} reads were skipped from: {}'.format(count, key))
        return coverage_df


@command('run-matrix')
def run_matrix(bamfiles_fofn, to_csv=False, wanted_gene='DMD', min_dist=629):
    """ Create coverage_matrix from given bam directory.

    :param bamfiles_fofn: File containing the paths to all bedfiles to be included in the coverage_matrix
    :param to_csv: Boolean of whether to create a csv file of the matrix
    :param wanted_gene: Name of the gene for where to get targets from
    :param min_dist: Any two intervals that are closer than this distance will be merged together,
        and any read pairs with insert lengths greater than this distance will be skipped

    """
    targets = cnv_util.combine_panel_intervals(wanted_gene=wanted_gene, min_dist=min_dist)

    matrix_instance = CoverageMatrix(min_interval_separation=min_dist)
    coverage_matrix_df = matrix_instance.create_coverage_matrix(bamfiles_fofn, targets)
    if to_csv:
        outfile_name = '{}_coverage_matrix.csv'.format(wanted_gene)
        coverage_matrix_df.to_csv("../exon_data/{}".format(outfile_name))
        print 'Finished creating {}'.format(outfile_name)

if __name__ == "__main__":
    main()
