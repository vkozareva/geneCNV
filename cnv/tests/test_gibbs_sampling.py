"""
Test code for method of Gibbs Sampling to detect DMD CNVs

Proposed production DMD CNV flow based on the current sampling code:
    Initialization
        Compute the coverage matrix of exon intensities from female RMA subjects.
        Compute a set of priors (X_probs) from those intensities.
    Per subject
        Create the coverage matrix for the subject.
        Run Gibbs Sampling on that data using the priors to estimate exon copy numbers.

DMD test plan:
    create_coverage_matrix()
        Create 10 .bam files that each include a single exon read of a single type for two exons.
        Run create_coverage_matrix() on each file and assert that exon counts are all zero, except for the included exon.
        Run create_coverage_matrix() on the aggregation of the test bam files and assert that the counts sum from above.

    generate_gibbs_df()
        Compute a set of X_probs from normal subjects.
        Extract the coverage matrices from a set of test subjects and commit those.
        Run generate_gibbs_df() on a set of normal and test subjects, asserting that the correct CNVs are found.
        Until we have test subject data, use simulated mutations.

Open Issues
    - What data files need to be committed to github?
    - Or could we just do some kind of rsync from dropbox to get the large bam files?
    - Need some threshold-based code for determining copy number from the posterior distribution.
"""

import sys, os, unittest, logging
import numpy as np
import pandas as pd
from genepeeks.common import utilities as util
from cnv.Gibbs.IntensitiesDistribution import IntensitiesDistribution
from cnv.Gibbs.CopyNumberDistribution import CopyNumberDistribution


class TestPloidyModel(unittest.TestCase):
    cnv_support = [1,2,3]
    file_dir = os.path.dirname(__file__)
    exon_data_dir = os.path.join(file_dir, '..', '..', 'exon_data75')
    # Do we want to use logger here?
    logger = util.create_logging(to_file=True,
                                 console_fmt='%(levelname)s: %(message)s',
                                 file_path='test_gibbs_sampling.log', file_mode='w')

    def shortDescription(self):
        """Turn off the nosetests "feature" of using the first line of the doc string as the test name."""
        return None

    def test_01_random_data(self):
        """Generate a random copy number vector with uniform distribution over support.
        Compute a fixed intensity vector which is uniform.
        Use multinomial to generate the counts for 20000 draws that reflect that distribution.
        Gibbs Sample the posterior of the prior + data generated from the prior.
        Compare the converged posterior probabilities with the random copy numbers.
        They should match, i.e. the copy number with the highest probability for each exon should be identical to the random copy number."""

        n_targets = 78
        copy_numbers = np.random.choice(TestPloidyModel.cnv_support, n_targets)
        self.logger.info('copy_numbers: {}'.format(copy_numbers))
        intensities = IntensitiesDistribution(np.ones(n_targets)/n_targets)
        p_vector = copy_numbers * intensities.intensities
        p_vector /= float(np.sum(p_vector))
        X_priors = np.random.multinomial(20000, p_vector)
        
        ploidy = CopyNumberDistribution(TestPloidyModel.cnv_support, data=X_priors, logger=self.logger)

        n_iterations = 1100
        gibbs_cnv_data = np.zeros((n_targets, n_iterations))
        likelihoods = np.zeros(n_iterations)

        # Gibbs Sampling of the posterior likelihood.
        for i in xrange(n_iterations):
            intensities.sample(ploidy)
            likelihood, cnv_probs = ploidy.sample(intensities)
            for exon in range(n_targets):
                gibbs_cnv_data[exon, i] = ploidy.cnv[exon]
            likelihoods[i] = likelihood
            if i % (n_iterations / 20) == 0:
                self.logger.info('After {} iterations:\ncnv: {}\nlikelihood: {}\ncnv_probs: {}'.format(
                    i + 1, ploidy.cnv, likelihood, cnv_probs))

        # Get proportions using burn-in of 1000 iterations.
        gibbs_data_results = np.zeros((n_targets, len(TestPloidyModel.cnv_support)))
        for target_i in range(n_targets):
            # Exclude samples before burn in and then take only every 100th sample to reduce autocorrelation.
            gibbs_slice = gibbs_cnv_data[target_i][1000:][::100]
            gibbs_data_results[target_i] = np.bincount(gibbs_slice.astype(np.int64), 
                                                       minlength=len(TestPloidyModel.cnv_support) + 1)[1:]
        gibbs_data_results /= float(len(gibbs_slice))

        gibbs_df = pd.DataFrame(gibbs_data_results, columns=['copy_{}'.format(cnv) for cnv in TestPloidyModel.cnv_support])
        gibbs_df['Exon'] = ['exon_{}'.format(i) for i in range(1, 79)]
        self.logger.info('Gibbs Results:\n{}'.format(gibbs_df))

        for target_i, gibbs_target_result in enumerate(gibbs_data_results):
            self.assertEqual(TestPloidyModel.cnv_support[list(gibbs_target_result).index(max(gibbs_target_result))],
                             copy_numbers[target_i],
                             'Gibbs results mismatch: {} versus {}'.format(gibbs_target_result, copy_numbers[target_i]))


if __name__ == '__main__':
    unittest.main()
