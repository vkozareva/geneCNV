"""The main Gibbs Sampling Model """

import numpy as np
import pandas as pd


from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution
from cnv.Targets.TargetCollection import TargetCollection

class PloidyModel(object):
    """This is the full statistical model and class that runs the Gibbs Sampling. It is responsible for taking a
    parameter set, a BAM file, and a TargetCollection and running a Gibbs sampler to determine the
    posterior probability of different ploidy states."""

    def __init__(self, cnv_support, cnv=None, targets=None, n_targets=78, parameterFileName=None, bamFileName=None,
                 data=None, intensities=None, logger=None):
        """Initialize the data model with its input arguments.
        Load the parameters and calculate the coverage at each interval in the BAM file."""
        self.logger = logger

        # Fake code:
        isinstance(targets, TargetCollection)

        # For now either targets or n_targets (optionally) may be specified.
        if targets is not None:
            self.targets = targets
            self.n_targets = len(self.targets)
        else:
            self.targets = ['Exon_{}'.format(i) for i in range(1, n_targets + 1)]
            self.n_targets = n_targets

        # For bamFileName or data may be specified.
        if bamFileName is not None:
            self.data = targets.getData(bamFileName)
        else:
            self.data = data

        self.cnv_support = cnv_support
        self.intensities = IntensitiesDistribution(parameterFileName=parameterFileName, intensities=intensities)
        self.ploidy = CopyNumberDistribution(self.targets, data=data, support=self.cnv_support, copies=cnv)

    def RunGibbsSampler(self, n_iterations=10000):
        """Gibbs sampling of the posterior likelihood"""
        self.gibbs_cnv_data = np.zeros((self.n_targets, n_iterations))
        self.gibbs_X = np.zeros((n_iterations, self.n_targets))
        self.likelihoods = np.zeros(n_iterations)

        for i in xrange(n_iterations):
            self.intensities.sample(self.ploidy)
            likelihood, cnv_probs = self.ploidy.sample(self.intensities)
            self.gibbs_cnv_data[:,i] = self.ploidy.copies
            self.likelihoods[i] = likelihood
            self.gibbs_X[i] = self.intensities.intensities

            # Log some convergence info at decile intervals.
            if (i + 1) % (n_iterations / 10) == 0:
                self.logger.info('After {} iterations:\ncnv: {}\nlikelihood: {}\ncnv_probs: {}'.format(
                    i + 1, self.ploidy.copies, likelihood, cnv_probs))

    def ReportGibbsData(self, out_file_name=None, burn_in=1000):
        """Output a results file and a PDF of the posterior probabilities plot.
        Report on the posterior distribution obtained by the sampling procedure"""
        # Get proportions using burn-in of 1000 iterations.
        gibbs_data_results = np.zeros((self.n_targets, len(self.cnv_support)))

        for target_i in xrange(self.n_targets):
            # Exclude samples before burn in and then take only every 100th sample to reduce autocorrelation.
            gibbs_slice = self.gibbs_cnv_data[target_i][burn_in:][::100]
            gibbs_data_results[target_i] = np.bincount(gibbs_slice.astype(np.int64),
                                                       minlength=len(self.cnv_support) + 1)[1:]
        gibbs_data_results /= float(len(gibbs_slice))

        # Return df for easier visualization.
        # TODO: This df is kind of lame and perhaps should be eliminated.
        gibbs_df = pd.DataFrame(gibbs_data_results, columns=['copy_{}'.format(cnv) for cnv in self.cnv_support])
        gibbs_df['Exon'] = self.targets

        if out_file_name is not None:
            # TODO: Write the output file.
            pass
        else:
            # Log the results
            self.logger.info('Gibbs Results:\n{}'.format(gibbs_df))

        return gibbs_data_results, gibbs_df
