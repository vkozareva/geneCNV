"""
Test code for method of Gibbs Sampling to detect DMD CNVs

Proposed production DMD CNV flow:
    Training
        Compute a set of intensity model hyperparameters from all female subjects run with desired TSO/TSID ratios.
    Per subject
        Create the coverage matrix for the subject.
        Run Gibbs Sampling on that data using the intensity model to estimate exon copy numbers.

Test plan:
    Training (not done)
        Create a set of H1 random hyperparameters.
        Draw n synthetic subjects with r reads each from the intensity model.
        Train H2 hyperparameters on the synthetic subject draws.
        Show that as r increases H2 converges on H1.

    Gibbs sampling
        Random test (done, w/o intensity model)
            Compute a set of random copy numbers.
            Draw synthetic reads based on the copy numbers and intensity model.
            Today the intensity model is uniform but in the future it will be a real model.
            Run Gibbs sampling on the reads to identify most probable copy numbers.
            Sampling results should be idential to the random copy numbers.
        
        Test with actual subject data (not done)
            Commit a set of intensity model hyperparameters.
            Extract the coverage matrices from a set of test subjects and commit those.
            Run evaluate_sample on the test subject data asserting that the correct copy numbers are found.

        Geweke test (not done)
"""

import sys, os, unittest, logging
import numpy as np
import pandas as pd
from genepeeks.common import utilities as util
from cnv.Gibbs.IntensitiesDistribution import IntensitiesDistribution
from cnv.Gibbs.PloidyModel import PloidyModel


class TestPloidyModel(unittest.TestCase):
    logger = util.create_logging(to_console=False, to_file=True,
                                 console_fmt='%(levelname)s: %(message)s',
                                 file_path='test_gibbs_sampling.log', file_mode='w')

    def shortDescription(self):
        """Turn off the nosetests "feature" of using the first line of the doc string as the test name."""
        return None

    def test_01_random_data(self):
        """Generate a random copy number vector with uniform distribution over support.
        Compute a fixed intensity vector which is uniform.
        Use multinomial to generate counts that reflect that distribution.
        Gibbs sample the posterior of the prior + data generated from the prior.
        Compare the converged posterior probabilities with the random copy numbers.
        The copy number with the highest probability for each exon should be the initial random copy number."""

        n_targets = 78
        n_draws = 20000
        cnv_support = [1, 2, 3]
        copy_numbers = np.random.choice(cnv_support, n_targets)
        self.logger.info('copy_numbers: {}'.format(copy_numbers))
        intensities = np.ones(n_targets)/n_targets
        p_vector = copy_numbers * intensities
        p_vector /= float(np.sum(p_vector))

        ploidy = PloidyModel(cnv_support, data=np.random.multinomial(n_draws, p_vector), intensities=intensities, logger=self.logger)
        ploidy.RunGibbsSampler()
        gibbs_data_results, gibbs_df = ploidy.ReportGibbsData()
        self.assertEqual([ploidy.cnv_support[np.where(gibbs_target_result==max(gibbs_target_result))[0][0]]
                          for gibbs_target_result in gibbs_data_results],
                         list(copy_numbers))


if __name__ == '__main__':
    unittest.main()
