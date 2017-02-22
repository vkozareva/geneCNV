import pysam
import numpy as np
import random
#import matplotlib.pyplot as plt
import pandas as pd
from mando import command, main

from genepeeks.common import utilities as util

class CNVDetector(object):
    def __init__(self, param_inputs):
        for key, value in param_inputs.items():
            setattr(self, key, value)

    def detect_CNV(self, test_df):
        coverage_df_f = self.coverage_df[self.coverage_df.gender == 'F']
        coverage_df_M1 = self.coverage_df[self.coverage_df.sequencer == 'M1']

        fem_norm = self.reshape_df(coverage_df_f, include_stats=True, subject_droplist=self.subject_list)
        M1_norm = self.reshape_df(coverage_df_M1, include_stats=True, subject_droplist=self.subject_list)

        merged = M1_norm.merge(test_df, on='Exon')

        for subject in self.subject_list:
            merged['{}_z_score'.format(subject)] = (pd.to_numeric(merged[subject]) - pd.to_numeric(merged.Mean)) / pd.to_numeric(merged.SD)
            for index, score in enumerate(merged['{}_z_score'.format(subject)]):
                compare = 'greater' if score > 0 else 'less'
                mutation = 'duplication' if score > 0 else 'deletion'
                if 2 < abs(score) < 3:
                    print 'Subject {} has {} than normal coverage for exon {}-- investigate further, score {}'.format(subject,
                                                                                            compare, merged.Exon[index], score)
                elif abs(score) > 3:
                    print 'Subject {} has much {} than normal coverage for exon {}-- indicates {}, score {}'.format(subject, compare,
                                                                                        merged.Exon[index], mutation, score)
        if self.output_csv:
            merged.to_csv('cnv_detection_{}_through_{}.csv'.format(self.subject_list[0], self.subject_list[-1]))

    def get_exon_data(self):
        DMD_ensembl = util.Mongo.get_collection_data('gene', wanted_db='prod', query={'_id': 'DMD'}, find_one=True, single_field='ensembl')
        DMD_exons = util.get_nested_value(DMD_ensembl, ('is_primary', 'transcripts', 'is_primary', 'exons'))
        self.DMD_exons_merged = util.merge_intervals(DMD_exons, min_dist=200, include_index=True)
        self.exon_labels = ['Ex' + exon['index'] for exon in self.DMD_exons_merged]

    @staticmethod
    def reshape_df(df, include_stats=False, subject_droplist=None):

        df_grouped = df.groupby(['subject']).sum()
        df_norm = df_grouped.div(df_grouped.sum(axis=1), axis=0)
        df_norm = df_norm.transpose().reset_index()
        df_norm.rename(columns={'index': 'Exon'}, inplace=True)

        if subject_droplist:
            for subject in subject_droplist:
                df_norm.drop(subject, axis=1, inplace=True)
        if include_stats:
            df_norm['Mean'] = df_norm.mean(axis=1)
            df_norm['SD'] = df_norm.std(axis=1)
        return df_norm

    def generate_test_frame(self):
        base_headers = ['subject', 'specimen', 'sample', 'gender', 'sequencer', 'flow_cell_id', 'lane']
        coverage_matrix = [base_headers + self.exon_labels]
        for file in self.test_bamfiles:
            subject_coverages = {}
            bamfile = pysam.AlignmentFile(file, 'rb')

            for RG in bamfile.header['RG']:
                subject, specimen_sample, flow_cell_lane = RG['ID'].split('_')
                gender = subject[0]
                specimen_num, sequencer, sample = specimen_sample.split('-')
                specimen = '{}_{}'.format(subject, specimen_num)
                sample = '{}_{}'.format(subject, specimen_sample)
                flow_cell_id, lane = flow_cell_lane.rsplit('-', 1)
                row = [subject, specimen, sample, gender, sequencer, flow_cell_id, lane]
                subject_coverages[RG['ID']] = row + [0] * len(self.DMD_exons_merged)

            for read in bamfile.fetch('X', start=31137345, end=33229636):
                if not read.is_unmapped and read.mapping_quality == 60:
                    # Find what exon each read falls in, and increase that exon's coverage by 1
                    interval_info = util.in_interval(read.reference_start, self.DMD_exons_merged, get_interval=True)
                    if not interval_info[0]:
                        # If the start of the read is not in an exon, check the end of the read
                        interval_info = util.in_interval(read.reference_end, self.DMD_exons_merged, get_interval=True)

                    if interval_info[0]:
                        exon_num = interval_info[1]
                        subject_coverages[read.get_tag('RG')][exon_num + len(base_headers)] += 1

            coverage_matrix += subject_coverages.values()
        coverage_matrix = np.array(coverage_matrix)

        test_subjects = pd.DataFrame(coverage_matrix[1:], columns=coverage_matrix[0])
        test_subjects_clean = test_subjects.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        self.subject_list = list(test_subjects_clean.subject.unique())

        test_norm = self.reshape_df(test_subjects_clean)
        return test_norm

@command('call_CNV')
def call_CNV(test_bamfile_paths, control_subjects_csv_path='exon_data/coverage_matrix.csv', output_csv=False):
    '''Calls large deletions and duplications in input bamfiles based on exon coverage calculations. '''

    coverage_df = pd.read_csv('exon_data/coverage_matrix.csv', header=1, index_col=0)
    coverage_df.index.name = None

    if output_csv:
        output_dir = 'cnv_outputs'
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

    Parameters = {
        'test_bamfiles': test_bamfile_paths.split('|'),
        'coverage_df': coverage_df,
        'output_csv': output_csv,
    }

    detector = CNVDetector(Parameters)
    detector.get_exon_data()
    test_df = detector.generate_test_frame()
    print type(test_df)
    detector.detect_CNV(test_df)

if __name__ == '__main__':
    main()

