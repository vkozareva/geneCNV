import numpy as np

class CopyNumberDistribution(object):
    """Describe the distribution of ploidy across a set of targets.
    Random data is substituted if optional subject data is omitted.
    Includes methods for sampling of the prior."""
    def __init__(self, n_targets, support=None):
        self.support = [1, 2, 3] if support is None else support
        self.n_targets = n_targets

    def sample_prior(self):
        """Sample a new ploidy state from a uniform prior, given support"""
        copy_sample = np.random.choice(self.support, size=self.n_targets)
        return copy_sample