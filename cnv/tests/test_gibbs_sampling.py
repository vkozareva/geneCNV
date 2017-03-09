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

from python_library import DMD_utilities as DMD_util


@command('extract_test_coverage_matrices')
def save_test_subject_coverage_matrices():
    test_subjects_num = cov.coverageMatrix().create_coverage_matrix(DMD_exons_merged, exon_labels, bam_dir='../bams/test_subjects')
    
def load_test_subject_coverage_matrices():
    return 


if __name__ == "__main__":
    main()
