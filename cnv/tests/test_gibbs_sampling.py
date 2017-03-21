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
from cnv.Gibbs.PloidyModel import PloidyModel


class TestPloidyModel(unittest.TestCase):
    file_dir = os.path.dirname(__file__)
    exon_data_dir = os.path.join(file_dir, '..', '..', 'exon_data75')
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
        cnv_support = [1, 2, 3]
        copy_numbers = np.random.choice(cnv_support, n_targets)
        self.logger.info('copy_numbers: {}'.format(copy_numbers))
        intensities = np.ones(n_targets)/n_targets
        p_vector = copy_numbers * intensities
        p_vector /= float(np.sum(p_vector))

        ploidy = PloidyModel(cnv_support, data=np.random.multinomial(20000, p_vector), intensities=intensities, logger=self.logger)
        ploidy.RunGibbsSampler()
        gibbs_data_results, gibbs_df = ploidy.ReportGibbsData()
        self.assertEqual([ploidy.cnv_support[list(gibbs_target_result).index(max(gibbs_target_result))]
                          for gibbs_target_result in gibbs_data_results],
                         list(copy_numbers))


if __name__ == '__main__':
    unittest.main()
