import os
import unittest

from cnv import inputs

_data_direc = os.path.join(os.path.dirname(__file__), "data")
EXAMPLE_BAM_PATH = os.path.join(_data_direc, 'example.bam')
MATLAB_HLN_PATH = os.path.join(_data_direc, 'hln_test.mat')
TEST_HLN_PARAMS = os.path.join(_data_direc, 'test_hln_params.p')

class ResourceTest(unittest.TestCase):

    def test_tso_bed_exists(self):
        tsobed = inputs.get_true_sight_one_bed()
        self.assertTrue(os.path.exists(tsobed))


    def test_tso_id_bed_exists(self):
        tsidbed = inputs.get_true_sight_inherited_disease_bed()
        self.assertTrue(os.path.exists(tsidbed))

    def test_dmd_exons_exist(self):
        dmdexons = inputs.get_dmd_exons()
        self.assertTrue(os.path.exists(dmdexons))

if __name__ == '__main__':
    unittest.main()