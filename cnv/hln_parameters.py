class HLN_Parameters(object):
    """Container for the subject testing model parameters
    This currently includes the target intervals, and the model hyperparameters."""
    def __init__(self, targets, mu, covariance):
        self.targets = targets
        self.mu = mu
        self.covariance = covariance
