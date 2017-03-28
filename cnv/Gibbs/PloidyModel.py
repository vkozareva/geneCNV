"""The main Gibbs Sampling Model """

import logging
import numpy as np
import pandas as pd


from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution
from TargetJointDistribution import TargetJointDistribution
from cnv.Targets.TargetCollection import TargetCollection

class PloidyModel(object):
    """This is the full statistical model and class that runs the Gibbs Sampling. It is responsible for taking a
    parameter set, a BAM file, and a TargetCollection and running a Gibbs sampler to determine the
    posterior probability of different ploidy states."""

    def __init__(self, cnv_support, hln_parameters, data=None, cnv=None, parameterFileName=None, bamFileName=None,
                intensities=None, include_covar=False):
        """Initialize the data model with its input arguments.
        Load the parameters and calculate the coverage at each interval in the BAM file."""
        # Fake code:
        # isinstance(targets, TargetCollection)

        # For now either targets or n_targets (optionally) may be specified.
        # if targets is not None:
        #     self.targets = targets
        #     self.n_targets = len(self.targets)
        # else:
        #     self.targets = ['Exon_{}'.format(i) for i in range(1, n_targets + 1)]
        #     self.n_targets = n_targets

        self.mu = hln_parameters.mu
        self.covariance = hln_parameters.covariance
        self.targets = hln_parameters.targets
        self.n_targets = len(self.targets)

        # For bamFileName or data may be specified.
        if bamFileName is not None:
            self.data = targets.getData(bamFileName)
        else:
            self.data = data

        self.cnv_support = cnv_support
        # self.intensities = IntensitiesDistribution(parameterFileName=parameterFileName, intensities=intensities)
        # self.ploidy = CopyNumberDistribution(self.targets, data=data, support=self.cnv_support, copies=cnv)

        # initialize values
        self.intensities = IntensitiesDistribution(self.mu, self.covariance).sample()
        self.ploidy = CopyNumberDistribution(self.n_targets, support=self.cnv_support).sample()

        self.joint_target = TargetJointDistribution(self.targets, self.mu, self.covariance, self.data, self.support,
                                                    include_covar=include_covar)

    def RunMCMC(self, n_iterations=10000):
        """Metropolis Hastings sampling of the posterior likelihood"""
        self.mcmc_copy_data = np.zeros((self.n_targets, n_iterations))
        self.mcmc_intens = np.zeros((n_iterations, self.n_targets))
        self.likelihoods = np.zeros(n_iterations)
        self.acceptance = np.zeros((n_iterations, self.n_targets))

        for i in xrange(n_iterations):
            for target_i in xrange(self.n_targets):
                self.ploidy[target_i], self.intensities[target_i], self.acceptance[i, target_i] = self.joint_target.sample(self.ploidy,
                                                                                                    self.intensities, target_i)
                self.mcmc_copy_data[target_i, i] = self.ploidy[target_i]
                self.mcmc_intens[i, target_i] = self.intensities[target_i]


            self.likelihoods[i] = self.joint_target.log_joint_likelihood(self.ploidy, self.intensities)

            # Log some convergence info at decile intervals.
            if (i + 1) % (n_iterations / 10) == 0:
                logging.info('After {} iterations:\ncnv: {}\nlikelihood: {}\n'.format(
                    i + 1, self.ploidy, self.likelihoods[i]))

        # Log acceptance ratio at end
        logging.info('Acceptance ratio: '.format(np.mean(self.acceptance)))


    def ReportMCMCData(self, out_file_name=None, burn_in=1000, autocor_slice=100):
        """Output a results file and a PDF of the posterior probabilities plot.
        Report on the posterior distribution obtained by the sampling procedure"""
        # Get proportions using burn-in of 1000 iterations.
        mcmc_posteriors = np.zeros((self.n_targets, len(self.cnv_support)))

        for target_i in xrange(self.n_targets):
            # Exclude samples before burn in and then take only every 100th sample to reduce autocorrelation.
            mcmc_slice = self.mcmc_posteriors[target_i][burn_in:][::autocor_slice]
            mcmc_posteriors[target_i] = np.bincount(mcmc_slice.astype(np.int64),
                                                       minlength=len(self.cnv_support) + 1)[1:]
        mcmc_posteriors /= float(len(mcmc_slice))

        # Return df for easier visualization.
        # TODO: This df is kind of lame and perhaps should be eliminated.
        mcmc_df = pd.DataFrame(mcmc_posteriors, columns=['copy_{}'.format(cnv) for cnv in self.cnv_support])
        mcmc_df['Target'] = self.targets

        # if out_file_name is not None:
        #     # TODO: Write the output file.
        #     pass
        # else:
        #     # Log the results
        #     logging.info('MCMC Results:\n{}'.format(mcmc_df))

        return mcmc_posteriors, mcmc_df
