import logging
import numpy as np

from cnv.Targets.TargetCollection import TargetCollection

class CopyNumberDistribution(object):
    """Describe the distribution of ploidy across a set of targets.
    Random data is substituted if optional subject data is omitted.
    Includes methods for Gibbs Sampling of the posterior likelihood."""
    def __init__(self, targets, data=None, cnv=None, support=[1, 2, 3], sim_reads=3e4, logger=logging.getLogger()):
        """Initialize the data to hold past and current samples of ploidy levels for all values in targets."""
        isinstance(targets, TargetCollection)
        self.support = support
        # Initialization of target copy numbers
        if cnv is not None:
            self.cnv = cnv
        else:
            # Generate initial guess for exon copy numbers using uniform prior distribution.
            self.cnv = np.random.choice(self.support, size=len(data))
        logger.info('Initial target guess: {}'.format(self.cnv))
        # For testing only.
        if data is not None:
            self.data = data
        else:
            self.data = np.random.multinomial(sim_reads, self.intensities.intensities)

    def sample(self, intensities):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state."""
        # Sample all cnv values.
        cnv_probs = np.zeros((len(intensities.intensities), len(self.support)))
        for exon in range(len(intensities.intensities)):
            cnv_log_probs = np.zeros(len(self.support))
            for value in self.support:
                self.cnv[exon] = value
                # get new normed probabilities given test value and priors for exon intensities
                log_likelihood = self.log_p(intensities.intensities)
                cnv_log_probs[value - 1] = log_likelihood
            cnv_log_probs = cnv_log_probs - np.max(cnv_log_probs)
            cnv_probs[exon] = np.exp(cnv_log_probs)
            cnv_probs[exon] = cnv_probs[exon] / np.sum(cnv_probs[exon])
            new_cnv = np.random.choice(self.support, p = cnv_probs[exon])
            self.cnv[exon] = new_cnv

        log_probs = np.log(np.multiply(self.cnv, intensities.intensities) / np.sum(np.multiply(self.cnv, intensities.intensities)))
        likelihood = np.sum(np.multiply(log_probs, self.data))

        return likelihood, cnv_probs

    def log_p(self, X_probs):
        normed_probs = np.multiply(self.cnv, X_probs) / np.sum(np.multiply(self.cnv, X_probs)) 
        log_likelihood =  np.sum(np.multiply(np.log(normed_probs), self.data))
        return log_likelihood
