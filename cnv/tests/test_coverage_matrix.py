from collections import Counter
import unittest
import pysam

from cnv import inputs
from cnv.Targets.TargetCollection import DEFAULT_MERGE_DISTANCE
from cnv.Targets.TargetCollection import TargetCollection
from cnv import coverage_matrix as cm, utilities as cnv_util
from cnv.coverage_matrix import WrappedBAM
from test_resources import EXAMPLE_BAM_PATH


class CoverageMatrixTests(unittest.TestCase):
    targets = TargetCollection.load_from_txt_file(inputs.get_dmd_exons(),)
    matrix_instance = cm.CoverageMatrix()

    def test_unwanted_filters(self):
        original_number_of_checks = len(self.matrix_instance.list_of_checks)
        chk = 'PCR_duplicate'
        list_of_checks = self.matrix_instance.filter_list_of_checks([chk])
        self.assertEqual(original_number_of_checks, len(list_of_checks) + 1)
        self.assertEqual(self.matrix_instance.list_of_checks.count(chk), 0)



    def test_targets(self):
        for i, target in enumerate(self.targets):
            self.assertTrue(hasattr(target, 'label'))
            self.assertTrue(hasattr(target, 'chrom'))
            self.assertGreater(target.end, target.start, "Target has start after end")
            if i != 0:
                self.assertGreater(target.start, self.targets[i - 1].end + DEFAULT_MERGE_DISTANCE,
                                   'The {} and {} intervals are less than {} base pairs apart'.format(i - 1, i, DEFAULT_MERGE_DISTANCE))

    def test_subject_coverage(self):
        aligned_bamfile = WrappedBAM(EXAMPLE_BAM_PATH)
        skipped_counts = Counter()
        coverage_vector = self.matrix_instance.get_subject_coverage(aligned_bamfile, self.targets, skipped_counts)
        self.assertEqual(len(coverage_vector), len(self.targets),
                         'There are {} targets but the coverage_vector has length {}'.format(len(self.targets), len(coverage_vector)))
        self.assertEqual(coverage_vector.count(0), 0, 'The example bamfile has 0 coverage in one of the targets')
