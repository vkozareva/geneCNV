import logging
import numpy as np

from cnv.Targets.TargetCollection import TargetCollection

class CopyNumberDistribution(object):
    """Describe the distribution of ploidy across a set of targets.
    Random data is substituted if optional subject data is omitted.
    Includes methods for Gibbs Sampling of the posterior likelihood."""
    def __init__(self, targets, data=None, cnv=None, support=None, sim_reads=3e4, logger=logging.getLogger()):
        """Initialize the data to hold past and current samples of ploidy levels for all values in targets."""
        isinstance(targets, TargetCollection)
        if support is not None:
            self.support = support
        else:
            self.support = [1, 2, 3]
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
            p_vector = self.cnv * self.intensities.intensities
            p_vector /= np.sum(p_vector)
            self.data = np.random.multinomial(sim_reads, p_vector)

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
                cnv_log_probs[value - 1] = self.log_likelihood(intensities.intensities)
            cnv_log_probs = cnv_log_probs - np.max(cnv_log_probs)
            cnv_probs[exon] = np.exp(cnv_log_probs)
            cnv_probs[exon] = cnv_probs[exon] / np.sum(cnv_probs[exon])
            new_cnv = np.random.choice(self.support, p = cnv_probs[exon])
            self.cnv[exon] = new_cnv

        return self.log_likelihood(intensities.intensities), cnv_probs

    def log_likelihood(self, X_probs):
        normed_probs = self.cnv * X_probs
        normed_probs /= np.sum(normed_probs)
        return np.sum(np.log(normed_probs) * self.data)
