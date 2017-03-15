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
from mando import command, main

from Gibbs.PloidyModel import PloidyModel


class TestPlodyModel(unittest.TestCase):
    cnv_support = [1,2,3]

    def shortDescription(self):
        """Turn off the nosetests "feature" of using the first line of the doc string as the test name."""
        return None

    def test_01_random_data(self):
        ploidy = PloidyModel(TestPlodyModel.cnv_support, TestPlodyModel.X_probs38)
        ploidy.RunGibbsSampler()
        self.gibbs_cnv_data, self.gibbs_X, gibbs_data_results, self.likelihoods = ploidy.OutputGibbsData(None)

    @classmethod
    def compute_X_probs():
        # get full dataset and also subsets based on gender, sequencer, etc
        coverage_df = pd.read_csv('../exon_data/coverage_matrix.csv', header=0, index_col=0)
        coverage_df.is_rerun.values == False
        # remove rerun data and coding region counts
        coverage_df = coverage_df[coverage_df.is_rerun == False]
        coverage_df.drop(['is_rerun'], axis=1, inplace=True)
        coverage_df.index.name = None
        coverage_df.date_modified = pd.to_datetime(coverage_df.date_modified, unit='s')
        coverage_df['date'] = coverage_df.date_modified.dt.date
        coverage_df_RMA = coverage_df[coverage_df.subject.str.contains('FRMR')]
    
        # checking dates of RMA samples -- all RMA samples seem to have been run within a day of each other in June 2016
        # in fact, this is just the last time the bams were modified, the samples were run earlier in the year, see below
        RMA_dates = coverage_df_RMA.date.unique()
    
        # Use RMA samples for initial intensity vector--note that all RMA individuals used the M1 mixin panel.
        # This is only the females in RMA.
        RMA_subset = coverage_df_RMA[coverage_df_RMA.subject.isin(subjects38)]
        gibbs_columns = ['subject'] + [column for column in rma_subset.columns if 'Ex' in column]
        rma38 = rma_subset[gibbs_columns]
        rma38_norm = reshape_df(rma38, include_stats=True)
        cls.X_probs38 = np.array(rma38_norm.Mean)


@command('extract_test_coverage_matrices')
def save_test_subject_coverage_matrices():
    test_subjects_num = cov.coverageMatrix().create_coverage_matrix(DMD_exons_merged, exon_labels, bam_dir='../bams/test_subjects')
    
def load_test_subject_coverage_matrices():
    return 


if __name__ == "__main__":
    main()
