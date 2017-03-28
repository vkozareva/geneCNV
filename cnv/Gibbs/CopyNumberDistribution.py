import numpy as np
from cnv.Targets.TargetCollection import TargetCollection

class CopyNumberDistribution(object):
    """Describe the distribution of ploidy across a set of targets.
    Random data is substituted if optional subject data is omitted.
    Includes methods for sampling of the prior."""
    def __init__(self, n_targets, support=None):
        """Initialize the data to hold past and current samples of ploidy levels for all values in targets."""
        if support is not None:
            self.support = support
        else:
            self.support = [1, 2, 3]
        self.n_targets = n_targets

    def sample_prior(self):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state."""

        copy_sample = np.random.choice(self.support, size=self.n_targets)

        return copy_sample