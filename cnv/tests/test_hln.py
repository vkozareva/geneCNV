import unittest
import numpy as np
import scipy.io as sio
from test_resources import *

from cnv.LogisticNormal import hln_EM

class LogisticNormalTest(unittest.TestCase):
    def test_FileExists(self):
        testFile = MATLAB_HLN_PATH
        self.assertTrue(os.path.exists(testFile))

    def test_ResultsMatch(self):
        testFile = sio.loadmat(MATLAB_HLN_PATH)

        testData = testFile['Y']
        ml_mu, ml_cov = hln_EM(testData)
        # no assertAlmostEqual for lists
        for i in xrange(len(ml_mu)):
            self.assertAlmostEqual(ml_mu[i], testFile['mu2'][i])
        for i in xrange(len(ml_cov.flatten())):
            self.assertAlmostEqual(ml_cov.flatten()[i], testFile['cov2'].flatten()[i])

if __name__ == '__main__':
    unittest.main()