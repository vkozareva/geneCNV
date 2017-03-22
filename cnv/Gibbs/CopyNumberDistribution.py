import logging
import numpy as np

from cnv.Targets.TargetCollection import TargetCollection

class CopyNumberDistribution(object):
    """Describe the distribution of ploidy across a set of targets.
    Random data is substituted if optional subject data is omitted.
    Includes methods for Gibbs Sampling of the posterior likelihood."""
    def __init__(self, targets, data=None, support=None, copies=None, sim_reads=20000, logger=logging.getLogger()):
        """Initialize the data to hold past and current samples of ploidy levels for all values in targets."""
        isinstance(targets, TargetCollection)
        self.data = data
        if support is not None:
            self.support = support
        else:
            self.support = [1, 2, 3]
        self.copies = copies

    def sample(self, intensities):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state."""
        # Sample all copy number values.
        cnv_probs = np.zeros((len(intensities.intensities), len(self.support)))
        for exon in range(len(intensities.intensities)):
            cnv_log_probs = np.zeros(len(self.support))
            for value in self.support:
                self.copies[exon] = value
                # get new normed probabilities given test value and priors for exon intensities
                cnv_log_probs[value - 1] = self.log_likelihood(intensities.intensities)
            cnv_log_probs = cnv_log_probs - np.max(cnv_log_probs)
            cnv_probs[exon] = np.exp(cnv_log_probs)
            cnv_probs[exon] = cnv_probs[exon] / np.sum(cnv_probs[exon])
            new_copies = np.random.choice(self.support, p = cnv_probs[exon])
            self.copies[exon] = new_copies

        return self.log_likelihood(intensities.intensities), cnv_probs

    def log_likelihood(self, X_probs):
        normed_probs = self.copies * X_probs
        normed_probs /= np.sum(normed_probs)
        return np.sum(np.log(normed_probs) * self.data)
