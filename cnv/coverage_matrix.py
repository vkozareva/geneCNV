import os
import pysam
import pandas as pd
from mando import command, main
from collections import Counter

from genepeeks.common import utilities as util
import utilities as cnv_util


class coverageMatrix(object):
    """docstring for coverageMatrix"""
    def __init__(self, min_interval_separation=629):
        super(coverageMatrix, self).__init__()
        self.logger = util.create_logging()
        self.min_interval_separation = min_interval_separation

    def get_unique_panel_intervals(self):
        """ Get the intervals that are unique to each panel """
        panel_intervals = {
            'TSID': {'file': os.path.join('..', 'inputs', 'TruSight_Inherited_Disease_Manifest_A.bed')},
            'TSO': {'file': os.path.join('..', 'inputs', 'TruSight-One-BED-May-2014.txt')}
        }

        for name, intrv_info in panel_intervals.items():
            df = pd.read_csv(intrv_info['file'], delimiter='\t', header=None, names=('chrom', 'start', 'end', 'id'))
            X_df = df[df['chrom'] == 'chrX'].sort_values(by='start')
            intrv_info['intrv_list'] = map(dict, dict(X_df.T).values())

        self.unique_panel_intervals = {
            'TSID': util.interval_diff(panel_intervals['TSO']['intrv_list'], panel_intervals['TSID']['intrv_list'], extend_by=self.min_interval_separation),
            'TSO': util.interval_diff(panel_intervals['TSID']['intrv_list'], panel_intervals['TSO']['intrv_list'], extend_by=self.min_interval_separation)
        }
        for panel, unique_intervals in self.unique_panel_intervals.items():
            self.logger.info('{} only: {} intervals over {} bp'.format(panel, len(unique_intervals), util.add_intervals(unique_intervals)))

    def filter_bamfiles(self, file_name, files, subj_name_filter):
        """ Filter out unwanted bam files. Return True if bamfile should be used, otherwise return False """
        if not file_name.endswith('.bam'):
            return False

        # The following subject does not have legit data
        if 'FPWB-0001-0309' in file_name:
            return False

        if subj_name_filter is not None:
            if isinstance(subj_name_filter, list):
                subject = os.path.splitext(file_name)[0]
                if subject not in subj_name_filter:
                    return False
            elif subj_name_filter not in file_name:
                return False

        if '{}.bai'.format(file_name) not in files:
            self.logger.info('{} is missing an index file'.format(file_name))
            return False
        return True

    def get_sample_info(self, RG, bwa_version, date_modified, root=None):
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
        if root and root.endswith('re86'):
            full_id += '_re86'
            is_re86 = True
        else:
            is_re86 = False

        sample_info = [full_id, subject, specimen, sample, gender, baits[0], flow_cell_id, bwa_version, date_modified, is_re86]
        return sample_info

    def get_unique_panel_reads(self, bamfile_path):
        """ Count the reads that fall in intervals anywhere in the X chromosome that are unique to each panel """
        unique_panel_reads = {}
        for panel, unique_intervals in self.unique_panel_intervals.iteritems():
            aligned_bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
            unique_intervals_merged = util.merge_intervals(unique_intervals, min_dist=self.min_interval_separation)
            coverage_vector = self.get_subject_coverage(aligned_bamfile, unique_intervals_merged, max_insert_length=self.min_interval_separation)
            unique_panel_reads[panel] = sum(coverage_vector)
        return unique_panel_reads

    @staticmethod
    def get_subject_coverage(bamfile, targets, skipped_counts=None, max_insert_length=629):
        """ Get vector of coverage counts for any given bamfile across any provided target regions """

        coverage_vector = []
        for i, target in enumerate(targets):
            read_pairs = {}
            for read in bamfile.fetch('X', start=target['start'], end=target['end']):
                # Skip reads that don't meet necessary quality checks and keep counts
                if read.is_unmapped:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'unmapped')
                    continue

                if read.mapping_quality != 60:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'MAPQ below 60')
                    continue

                if read.is_duplicate:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'PCR_duplicate')
                    continue

                if read.mate_is_unmapped:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'mate_is_unmapped')
                    continue

                if not read.is_proper_pair:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'not a proper pair')
                    continue

                insert_length = read.template_length
                if read.is_reverse == read.mate_is_reverse:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'tandem_pair')
                    continue
                elif read.is_reverse:
                    insert_length *= -1

                if insert_length <= 0:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'negative_insert_length')
                    continue

                # To ensure that no read pairs overlap multiple targets, skip all reads with
                # insert length greater than the distance used to merge intervals
                if insert_length >= max_insert_length:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'insert_length greater than {}'.format(max_insert_length))
                    continue

                # Determine the start and end of each read pair
                pair_start = min(read.reference_start, read.next_reference_start)
                pair_end = pair_start + insert_length
                if pair_end < read.reference_end:
                    if skipped_counts is not None:
                        util.add_to_dict(skipped_counts, 'pair_end is less than reference_end')
                    continue

                # Keep track of each read pair, and count coverage at the end in order to only count each read pair once
                util.add_to_dict(read_pairs, (pair_start, pair_end))

            if skipped_counts is not None:
                read_pair_counts = dict(Counter(read_pairs.values()))
                duplicate_read_pairs = {key: value for key, value in read_pair_counts.items() if key not in [1, 2]}
                if duplicate_read_pairs:
                    util.add_to_dict(skipped_counts, 'duplicate_read_pairs', sum(duplicate_read_pairs.values()))

            # Count each read pair once towards the target coverage
            target_coverage = len(read_pairs)
            coverage_vector.append(target_coverage)

        return coverage_vector

    def create_coverage_matrix(self, intervals, interval_labels, bam_dir=None, subj_name_filter=None):
        """ Create coverage matrix with exons as columns, samples as rows, and amount of coverage in each exon as the values,
        plus extra columns for identifying info for each sample """

        if len(intervals) != len(interval_labels):
            util.stop_err('Unequal number of intervals ({}) vs interval labels ({})'.format(len(intervals), len(interval_labels)))

        self.get_unique_panel_intervals()

        # Initiate matrix headers
        base_headers = [
            'id', 'subject', 'specimen', 'sample', 'gender', 'baits', 'flow_cell_id',
            'bwa_version', 'date_modified', 'is_rerun', 'TSID_only', 'TSO_only']
        full_headers = base_headers + interval_labels
        skipped_counts = {}
        coverage_matrix = []
        if bam_dir is None:
            bam_dir = '/mnt/vep/subjects'

        # Count the number of bamfiles that will be used in order to use timing_fields
        file_count = 0
        for root, dirs, files in os.walk(bam_dir):
            for file_name in files:
                use_bamfile = self.filter_bamfiles(file_name, files, subj_name_filter)
                if use_bamfile:
                    file_count += 1
        starting_message = '\nCreating coverage_matrix with {} subjects'.format(file_count)
        timing_fields = util.initiate_timer(message=starting_message, add_counts=True, logger=self.logger,
                                            total_counts=file_count, count_steps=3 if file_count > 100 else None)

        # Iterate over all bamfiles in the directory and create the coverage_matrix
        for root, dirs, files in os.walk(bam_dir):
            for file_name in files:
                use_bamfile = self.filter_bamfiles(file_name, files, subj_name_filter)
                if use_bamfile:
                    bamfile_path = os.path.join(root, file_name)
                    bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
                    date_modified = os.path.getmtime(bamfile_path)

                    # Get identifying sample info
                    bwa_version = next(PG['VN'] for PG in bamfile.header['PG'] if PG.get('ID') == 'bwa')
                    sample_info = self.get_sample_info(bamfile.header['RG'][0], bwa_version, date_modified, root=root)

                    # Get subject coverage info
                    subj_coverage_vector = self.get_subject_coverage(bamfile, intervals, skipped_counts=skipped_counts)
                    unique_panel_reads = self.get_unique_panel_reads(bamfile_path)
                    full_subj_vector = sample_info + [unique_panel_reads['TSID'], unique_panel_reads['TSO']] + subj_coverage_vector
                    if len(full_subj_vector) != len(full_headers):
                        util.stop_err('Unequal number of columns ({}) vs headers ({})'.format(len(full_subj_vector), len(full_headers)))

                    coverage_matrix.append(full_subj_vector)

                    util.get_timing(timing_fields, display_counts=True)

        coverage_df = pd.DataFrame(coverage_matrix, columns=full_headers)

        # Add a column with the ratio of inherited disease only reads compared to Trusight One only reads
        coverage_df['TSID_ratio'] = coverage_df['TSID_only'] / (coverage_df['TSID_only'] + coverage_df['TSO_only'])

        # Log counts of skipped reads
        for key, count in skipped_counts.items():
            self.logger.info('{} reads were skipped due to {}'.format(count, key))
        return coverage_df


@command('run-matrix')
def run_matrix(bam_dir='../../library_files/inputs/bam_files', subj_filter=None, to_csv=False, wanted_gene='DMD', min_dist=629):
    """ Create coverage_matrix from given bam directory. Use subj_filter to only include certain bamfiles, and use to_csv to create a csv file of the matrix """
    targets, target_labels = cnv_util.combine_panel_intervals(wanted_gene=wanted_gene, min_dist=min_dist)

    # sample subj_name_filter: 'FRMR-00AW-8645' or 'RMR'
    subj_name_filter = subj_filter.split(',') if subj_filter and ',' in subj_filter else subj_filter

    matrix_instance = coverageMatrix(min_interval_separation=min_dist)
    coverage_matrix_df = matrix_instance.create_coverage_matrix(targets, target_labels, bam_dir=bam_dir, subj_name_filter=subj_name_filter)
    if to_csv:
        outfile_name = '{}_coverage_matrix{}.csv'.format(wanted_gene, '_' + subj_filter if subj_filter else '')
        coverage_matrix_df.to_csv("../exon_data/{}".format(outfile_name))
        print 'Finished creating {}'.format(outfile_name)

if __name__ == "__main__":
    main()
