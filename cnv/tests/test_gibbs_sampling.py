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

import sys, os, unittest
from cnv.Gibbs.PloidyModel import PloidyModel


class TestPlodyModel(unittest.TestCase):
    cnv_support = [1,2,3]
    file_dir = os.path.dirname(__file__)
    exon_data_dir = os.path.join(file_dir, '..', '..', 'exon_data')

    def shortDescription(self):
        """Turn off the nosetests "feature" of using the first line of the doc string as the test name."""
        return None

    def test_01_random_data(self):
        """Generate a random copy number vector with uniform distribution over support.
        Compute a fixed intensity vector which is uniform.
        Use multinomial to generate the counts for 20000 draws that reflect that distribution.
        Use Gibbs' method to sample the posterior of the prior + data generated from the prior.
        Compare the converged posterior probabilities with the random copy numbers.
        They should match, i.e. the copy number with the highest probability for each exon should be identical to the random copy number."""

        n_targets = 78
        copy_numbers = np.random.choice(TestPloidyModel.ancnv_support, n_targets)
        intensity = np.ones(len(copy_numbers))/len(copy_numbers)
        p_vector = copy_numbers * intensity
        p_vector /= float(np.sum(p_vector))
        X_priors = np.random.multinomial(20000, p_vector)
        
        ploidy = PloidyModel(TestPlodyModel.cnv_support, TestPlodyModel.X_priors)
        ploidy.RunGibbsSampler()
        self.gibbs_cnv_data, self.gibbs_X, gibbs_data_results, self.likelihoods = ploidy.OutputGibbsData(None)


if __name__ == '__main__':
    unittest.main()
