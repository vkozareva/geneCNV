import unittest
import pysam
from cnv import coverage_matrix as cm, utilities as cnv_util
from genepeeks.common import utilities as util
from test_resources import EXAMPLE_BAM_PATH


class CoverageMatrixTests(unittest.TestCase):
    min_dist = 629
    targets = cnv_util.combine_panel_intervals(min_dist=min_dist)
    matrix_instance = cm.CoverageMatrix(min_interval_separation=min_dist)

    def test_unwanted_filters(self):
        original_number_of_checks = len(self.matrix_instance.list_of_checks)
        list_of_checks = self.matrix_instance.filter_list_of_checks(['PCR_duplicate'])
        self.assertEqual(original_number_of_checks, len(list_of_checks) + 1)

    def test_targets(self):
        for i, target in enumerate(self.targets):
            self.assertIn('label', target, 'Every interval must have a "label" key')
            self.assertGreater(target['end'], target['start'], 'The {} interval has its start after its end'.format(target['label']))
            if i != 0:
                self.assertGreater(target['start'], self.targets[i - 1]['end'] + self.min_dist,
                                   'The {} and {} intervals are less than {} base pairs apart'.format(i - 1, i, self.min_dist))

    def test_unique_panel_intervals(self):
        unique_panel_intervals = self.matrix_instance.get_unique_panel_intervals(print_counts=False)
        self.assertEqual(util.interval_intersect(*unique_panel_intervals.values()), [],
                         'The unique_panel_intervals are supposed to be unique to each panel, but they overlap with each other')

        unique_panel_reads = self.matrix_instance.get_unique_panel_reads(EXAMPLE_BAM_PATH, unique_panel_intervals)
        self.assertNotEqual(sum(unique_panel_reads.values()), 0, 'The example bam has no coverage unique to either panel')
        # Following lines are commented out for now because example.bam is only for DMD, but unique panel intervals
        # looks at coverage across the entire X chromosome
        # for panel, unique_panel_coverage in unique_panel_reads.items():
        #     if not unique_panel_coverage:
        #         self.assertNotEqual(unique_panel_coverage, 0, 'The example bam has no coverage unique to {}'.format(panel))

    def test_subject_coverage(self):
        aligned_bamfile = pysam.AlignmentFile(EXAMPLE_BAM_PATH, 'rb')
        coverage_vector = self.matrix_instance.get_subject_coverage(aligned_bamfile, self.targets)
        self.assertEqual(len(coverage_vector), len(self.targets),
                         'There are {} targets but the coverage_vector has length {}'.format(len(self.targets), len(coverage_vector)))
        self.assertEqual(coverage_vector.count(0), 0, 'The example bamfile has 0 coverage in one of the targets')
