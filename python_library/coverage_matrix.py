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

    def get_subject_coverage_matrix(self, bamfile_path, exons_merged, skipped_counts, root=None):
        """ Create matrix of exon coverage for any given subject """

        date_modified = os.path.getmtime(bamfile_path)

        bamfile = pysam.AlignmentFile(bamfile_path, "rb")

        # Gather identifying info for each sample
        subject_coverages = {}
        bwa_version = next(PG['VN'] for PG in bamfile.header['PG'] if PG.get('ID') == 'bwa')
        for RG in bamfile.header['RG']:
            sample_info = self.get_sample_info(RG, bwa_version, date_modified, root=root)
            if len(sample_info) != len(self.base_headers):
                util.stop_err('Unequal number of sample info fields vs base headers: {}'.format(zip(self.base_headers, sample_info)))

            # Initialize each row with identifying info for the sample plus each exon's coverage of 0.
            # Also create 2 extra exons at end for coding regions of first and last exon
            subject_coverages[RG['ID']] = sample_info + [0] * (len(exons_merged) + 2)

        # Get coverage data for each sample within each exon
        for read in bamfile.fetch('X', start=exons_merged[0]['start'], end=exons_merged[-1]['end']):
            if not read.is_unmapped:
                if read.mapping_quality == 60:
                    # Find what exon each read falls in, and increase that exon's coverage by 1
                    exon_num = DMD_util.get_exon_num(read.reference_start, read.reference_end, exons_merged, skipped_counts)
                    if exon_num is not None:
                        subject_coverages[read.get_tag('RG')][exon_num + len(self.base_headers)] += 1

                        # For first and last exon, also check if the read falls in the coding region
                        if exon_num == 0:
                            if read.reference_end >= exons_merged[exon_num]['coding_start']:
                                subject_coverages[read.get_tag('RG')][-2] += 1
                        elif exon_num == len(exons_merged) - 1:
                            if read.reference_start <= exons_merged[exon_num]['coding_end']:
                                subject_coverages[read.get_tag('RG')][-1] += 1
                else:
                    util.add_to_dict(skipped_counts, 'MAPQ below 60')

        return subject_coverages

    def create_coverage_matrix(self, DMD_exons_merged, exon_labels, bam_dir='/mnt/vep/subjects', subj_name_filter=None):
        """ Create coverage matrix with exons as columns, samples as rows, and amount of coverage in each exon as the values,
        plus extra columns for identifying info for each sample """

        self.base_headers = ['id', 'subject', 'specimen', 'sample', 'gender', 'sequencer', 'flow_cell_id', 'lane', 'bwa_version', 'date_modified', 'is_rerun']
        full_headers = self.base_headers + exon_labels + ['Ex1_coding', 'Ex79_coding']
        subject_count = 0
        skipped_counts = {}
        coverage_matrix = []
        for root, dirs, files in os.walk(bam_dir):
            for file_name in files:
                if file_name.endswith('.bam'):
                    # The following subject does not have legit data
                    if 'FPWB-0001-0309' in file_name:
                        continue
                    if subj_name_filter is not None:
                        if isinstance(subj_name_filter, list):
                            subject = os.path.splitext(file_name)[0]
                            if subject not in subj_name_filter:
                                continue
                        elif subj_name_filter not in file_name:
                            continue
                    if '{}.bai'.format(file_name) not in files:
                        self.logger.info('{} is missing an index file'.format(file_name))
                        continue

                    bamfile_path = os.path.join(root, file_name)
                    subject_coverages = self.get_subject_coverage_matrix(bamfile_path, DMD_exons_merged, skipped_counts, root=root)
                    coverage_matrix += subject_coverages.values()
                    subject_count += 1
                    if subject_count % 20 == 0:
                        self.logger.info('Finished parsing {} subjects'.format(subject_count))

        coverage_matrix_df = pd.DataFrame(coverage_matrix, columns=full_headers)

        # Log counts of skipped reads
        self.logger.info('Finished parsing all {} subjects'.format(subject_count))
        for key, count in skipped_counts.items():
            self.logger.info('{} reads were skipped due to {}'.format(count, key))
        return coverage_matrix_df


@command('run-matrix')
def run_matrix(bam_dir='../../library_files/inputs/bam_files', subj_filter=None, to_csv=False):
    """ Create coverage_matrix from given bam directory. Use subj_filter to only include certain bamfiles, and use to_csv to create a csv file of the matrix """
    DMD_exons_merged, exon_labels = DMD_util.get_DMD_exons_merged()

    # sample subj_name_filter: 'FRMR-00AW-8645' or 'RMR'
    subj_name_filter = subj_filter.split(',') if subj_filter and ',' in subj_filter else subj_filter
    matrix_instance = coverageMatrix()
    coverage_matrix_df = matrix_instance.create_coverage_matrix(DMD_exons_merged, exon_labels, bam_dir=bam_dir, subj_name_filter=subj_name_filter)
    if to_csv:
        outfile_name = 'coverage_matrix{}.csv'.format('_' + subj_filter if subj_filter else '')
        coverage_matrix_df.to_csv("../exon_data/{}".format(outfile_name))
        print 'Finished creating {}'.format(outfile_name)

if __name__ == "__main__":
    main()
