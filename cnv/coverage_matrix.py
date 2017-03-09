from genepeeks.common import utilities as util
import DMD_utilities as DMD_util
import pysam
import os
import pandas as pd
from mando import command, main


class coverageMatrix(object):
    """docstring for coverageMatrix"""
    def __init__(self):
        super(coverageMatrix, self).__init__()
        self.logger = util.create_logging()

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
            'TSID': util.interval_diff(panel_intervals['TSO']['intrv_list'], panel_intervals['TSID']['intrv_list'], extend_by=300),
            'TSO': util.interval_diff(panel_intervals['TSID']['intrv_list'], panel_intervals['TSO']['intrv_list'], extend_by=300)
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
            flow_cell_id = lane = None
        else:
            flow_cell_id, lane = flow_cell_lane.rsplit('-', 1)

        # simulated CNV subjects have one of these suffixes in this field
        if 'del' in RG['SM'] or 'dup' in RG['SM']:
            subject = RG['SM']

        gender = subject[0]
        if specimen_sample.startswith(('ACGT', 'Omega')):
            lab, specimen_num, sequencer, sample = specimen_sample.split('.')
            specimen_num = '{}_{}'.format(lab, specimen_num)
        else:
            specimen_num, sequencer, sample = specimen_sample.split('-')
        specimen = '{}_{}'.format(subject, specimen_num)
        sample = '{}_{}'.format(subject, specimen_sample)
        full_id = RG['ID']
        if root and root.endswith('re86'):
            full_id += '_re86'
            is_re86 = True
        else:
            is_re86 = False

        sample_info = [full_id, subject, specimen, sample, gender, sequencer, flow_cell_id, lane, bwa_version, date_modified, is_re86]
        return sample_info

    def find_unique_panel_reads(self, subject_coverages, bamfile_path):
        """ Count the reads that fall in intervals anywhere in the X chromosome that are unique to each panel """

        for panel, unique_intervals in self.unique_panel_intervals.iteritems():
            aligned_bamfile = pysam.AlignmentFile(bamfile_path, 'rb')
            for i, target in enumerate(unique_intervals):
                for read in aligned_bamfile.fetch('X', start=target['start'] - 1, end=target['end'] + 1):
                    if not read.is_unmapped and read.mapping_quality == 60:
                        # If a read is in multiple targets, only count it once by waiting until the last target it is in to count it
                        if i != len(unique_intervals) - 1 and read.reference_end >= unique_intervals[i + 1]['start']:
                            pass
                        else:
                            subject_coverages[read.get_tag('RG')][len(self.base_headers) - (2 if panel == 'TSID' else 1)] += 1

    def get_subject_coverage_matrix(self, bamfile_path, intervals, skipped_counts, root=None):
        """ Create matrix of exon coverage for any given subject """

        date_modified = os.path.getmtime(bamfile_path)

        bamfile = pysam.AlignmentFile(bamfile_path, 'rb')

        # Gather identifying info for each sample
        subject_coverages = {}
        bwa_version = next(PG['VN'] for PG in bamfile.header['PG'] if PG.get('ID') == 'bwa')
        for RG in bamfile.header['RG']:
            sample_info = self.get_sample_info(RG, bwa_version, date_modified, root=root)

            # Initialize each row with identifying info for the sample plus each exon's coverage of 0.
            # Also create 2 extra exons at end for coding regions of first and last exon
            initialized_row = sample_info + [0] * (len(intervals) + 2)
            if len(initialized_row) != len(self.full_headers):
                util.stop_err('Unequal number of columns ({}) vs headers ({})'.format(len(initialized_row), len(self.full_headers)))

            subject_coverages[RG['ID']] = initialized_row

        self.find_unique_panel_reads(subject_coverages, bamfile_path)

        # Get coverage data for each sample within each exon
        for i, target in enumerate(intervals):
            for read in bamfile.fetch('X', start=target['start'] - 1, end=target['end'] + 1):
                if not read.is_unmapped:
                    if read.mapping_quality == 60:
                        # Check if the read falls in multiple targets, and skip if it does
                        in_multiple_targets = (
                            (i != 0 and read.reference_start <= intervals[i - 1]['end']) or
                            (i != len(intervals) - 1 and read.reference_end >= intervals[i + 1]['start'])
                        )
                        if in_multiple_targets:
                            util.add_to_dict(skipped_counts, 'in_multiple_targets')
                        else:
                            subject_coverages[read.get_tag('RG')][i + len(self.base_headers)] += 1
                    else:
                        util.add_to_dict(skipped_counts, 'MAPQ below 60')

        return subject_coverages

    def create_coverage_matrix(self, intervals, interval_labels, bam_dir=None, subj_name_filter=None):
        """ Create coverage matrix with exons as columns, samples as rows, and amount of coverage in each exon as the values,
        plus extra columns for identifying info for each sample """

        if len(intervals) != len(interval_labels):
            util.stop_err('Unequal number of intervals ({}) vs interval labels ({})'.format(len(intervals), len(interval_labels)))

        self.get_unique_panel_intervals()

        # Initiate matrix headers
        self.base_headers = [
            'id', 'subject', 'specimen', 'sample', 'gender', 'sequencer', 'flow_cell_id',
            'lane', 'bwa_version', 'date_modified', 'is_rerun', 'TSID_only', 'TSO_only']
        self.full_headers = self.base_headers + interval_labels
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
                    subject_coverages = self.get_subject_coverage_matrix(bamfile_path, intervals, skipped_counts, root=root)
                    coverage_matrix += subject_coverages.values()

                    util.get_timing(timing_fields, display_counts=True)

        coverage_df = pd.DataFrame(coverage_matrix, columns=self.full_headers)

        # Add a column with the ratio of inherited disease only reads compared to Trusight One only reads
        coverage_df['TSID_ratio'] = coverage_df['TSID_only'] / (coverage_df['TSID_only'] + coverage_df['TSO_only'])

        # Log counts of skipped reads
        for key, count in skipped_counts.items():
            self.logger.info('{} reads were skipped due to {}'.format(count, key))
        return coverage_df


@command('run-matrix')
def run_matrix(bam_dir='../../library_files/inputs/bam_files', subj_filter=None, to_csv=False, wanted_gene='DMD'):
    """ Create coverage_matrix from given bam directory. Use subj_filter to only include certain bamfiles, and use to_csv to create a csv file of the matrix """
    exons_merged, exon_labels = DMD_util.get_merged_exons(wanted_gene=wanted_gene)

    # sample subj_name_filter: 'FRMR-00AW-8645' or 'RMR'
    subj_name_filter = subj_filter.split(',') if subj_filter and ',' in subj_filter else subj_filter
    matrix_instance = coverageMatrix()
    coverage_matrix_df = matrix_instance.create_coverage_matrix(exons_merged, exon_labels, bam_dir=bam_dir, subj_name_filter=subj_name_filter)
    if to_csv:
        outfile_name = '{}_coverage_matrix{}.csv'.format(wanted_gene, '_' + subj_filter if subj_filter else '')
        coverage_matrix_df.to_csv("../exon_data/{}".format(outfile_name))
        print 'Finished creating {}'.format(outfile_name)

if __name__ == "__main__":
    main()
