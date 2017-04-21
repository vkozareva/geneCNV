import numpy as np

class CopyNumberDistribution(object):
    """Describe the distribution of ploidy across a set of targets.
    Random data is substituted if optional subject data is omitted.
    Includes methods for sampling of the prior."""
    def __init__(self, n_targets, support=None):
        self.support = [1, 2, 3] if support is None else support
        self.n_targets = n_targets

    def sample_prior(self, first_baseline_i):
        """Sample a new ploidy state from a uniform prior, given support"""
        # Set all baseline targets to ploidy 2
        copy_sample = 2. * np.ones(self.n_targets)
        copy_sample[:first_baseline_i] = np.random.choice(self.support, size=first_baseline_i)

        return copy_sample