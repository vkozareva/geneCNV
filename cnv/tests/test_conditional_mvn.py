import unittest
import numpy as np

from cnv.MCMC.TargetJointDistribution import TargetJointDistribution

class ConditionalMVNTest(unittest.TestCase):
    def test_conditional_mvn(self):
        test_mu = np.array([175., 71.]).reshape((-1, 1))
        test_cov = np.array([[550., 40.], [40., 8.]])
        test_input = np.array([0, 70])
        test_index = 0

        true_results = [170., 350.]

        test_mu_bar, test_cov_bar = TargetJointDistribution.get_conditional_mvn(test_mu, test_cov, test_index, test_input)
        self.assertListEqual([test_mu_bar, test_cov_bar], true_results)

if __name__ == '__main__':
    unittest.main()