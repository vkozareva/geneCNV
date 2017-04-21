"""
Test code for method of MCMC Sampling to detect DMD CNVs

Proposed production DMD CNV flow:
    Training
        Compute a set of intensity model hyperparameters from all female subjects run with desired TSO/TSID ratios.
    Per subject
        Create the coverage matrix for the subject.
        Run MCMC Sampling on that data using the intensity model to estimate exon copy numbers.

Test plan:
    Training (not done)
        Create a set of H1 random hyperparameters.
        Draw n synthetic subjects with r reads each from the intensity model.
        Train H2 hyperparameters on the synthetic subject draws.
        Show that as r increases H2 converges on H1.

    MCMC sampling
        Random test (done, w/o intensity model)
            Compute a set of random copy numbers.
            Draw synthetic reads based on the copy numbers and intensity model.
            Today the intensity model is uniform but in the future it will be a real model.
            Run MCMC sampling on the reads to identify most probable copy numbers.
            Sampling results should be idential to the random copy numbers.

        Test with actual subject data (not done)
            Commit a set of intensity model hyperparameters.
            Extract the coverage matrices from a set of test subjects and commit those.
            Run evaluate_sample on the test subject data asserting that the correct copy numbers are found.

        Geweke test (not done)
"""

import sys, os, unittest, logging
import numpy as np
import cPickle
from cnv.MCMC.IntensitiesDistribution import IntensitiesDistribution
from cnv.MCMC.CopyNumberDistribution import CopyNumberDistribution
from cnv.MCMC.PloidyModel import PloidyModel
from cnv.hln_parameters import HLN_Parameters

from test_resources import *


class TestPloidyModel(unittest.TestCase):
    def shortDescription(self):
        """Turn off the nosetests "feature" of using the first line of the doc string as the test name."""
        return None

    def test_01_random_data(self):
        """Generate a random copy number vector with uniform distribution over support.
        Generate random intensity vector from prior (using test parameters).
        Use multinomial to generate counts that reflect that distribution.
        MCMC sample the posterior of the prior + data generated from the prior.
        Compare the converged posterior probabilities with the random copy numbers.
        The copy number with the highest probability for each exon should be the initial random copy number."""

        n_draws = 40000
        cnv_support = [1e-10, 1, 2, 3]

        # loading the params as a dict here because for pickled user class instances
        # cPickle requires the class to be on the same directory level as when instance was pickled
        test_hln_params = cPickle.load(open(TEST_HLN_PARAMS, 'rb'))
        n_targets = len(test_hln_params['targets'])
        test_params = HLN_Parameters(test_hln_params['targets'], test_hln_params['mu'], test_hln_params['covariance'])

        copy_numbers = CopyNumberDistribution(n_targets, support=cnv_support).sample_prior(n_targets)
        logging.info('Initial target copy numbers: {}'.format(copy_numbers))
        intensities = IntensitiesDistribution(test_params.mu, test_params.covariance).sample()
        p_vector = np.multiply(copy_numbers, np.exp(intensities))
        p_vector /= float(np.sum(p_vector))

        ploidy_instance = PloidyModel(cnv_support, test_params, data=np.random.multinomial(n_draws, p_vector))
        ploidy_instance.RunMCMC()
        copy_posteriors = ploidy_instance.ReportMCMCData()
        self.assertEqual([ploidy_instance.cnv_support[np.where(mcmc_target_result==max(mcmc_target_result))[0][0]]
                          for mcmc_target_result in copy_posteriors],
                         list(copy_numbers))


if __name__ == '__main__':
    unittest.main()
