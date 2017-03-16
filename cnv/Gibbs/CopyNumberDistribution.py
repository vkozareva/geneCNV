import numpy as np

from cnv.Targets.TargetCollection import TargetCollection
from IntensitiesDistribution import IntensitiesDistribution

class CopyNumberDistribution(object):
    """A class which describes the distribution of ploidy across the different exons
    given a support vector for CNVs, an intensity vector, and subject data (optional).
    Includes methods for Gibbs Sampling."""
    def __init__(self, cnv_support, X_priors, scale, data=None, cnv=None, sim_reads=3e4, exon_labels=None):
        self.cnv_support = cnv_support
        # Initialization of cnv counts.
        if cnv is not None:
            self.cnv = cnv
        else:
            # generate initial guess for exon copy numbers using uniform prior distribution
            self.cnv = np.random.choice(self.cnv_support, size=len(X_priors)) 
        # TODO: initialize intensity vector.
        print self.cnv
        self.scale = scale
        self.exon_labels = exon_labels

        # we're using a scaling factor for the X_priors (small ratios give too much variability in the dirichlet)
        # TODO: Consider pulling this out into PloidyModel:
        self.X_priors = X_priors
        intensities = IntensitiesDistribution(None, scale)
        X_vect = intensities.sample(self.X_priors, data)
        # consider returning this vector along with initial cnv vector
        normed_probs_first = np.multiply(self.cnv, X_vect) / np.sum(np.multiply(self.cnv, X_vect))
        # for testing only
        if data is not None:
            self.data = data
        else:
            self.data = np.random.multinomial(sim_reads, normed_probs_first)

    def sample(self, intensities):
        """Given a current set of intensities, and the current ploidy state maintained in this class,
        sample a new ploidy state."""
        # sample all cnv values
        cnv = np.copy(self.cnv)
        X_vect = np.random.dirichlet(self.X_priors * self.scale)
        for exon in range(len(self.X_priors)):
            test = np.zeros(len(self.cnv_support))
            for value in self.cnv_support:
                cnv[exon] = value
                # get new normed probabilities given test value and priors for exon intensities
                log_likelihood = self.joint_log_p(cnv, X_vect)
                test[value - 1] = log_likelihood
            test = test - np.max(test)
            sample_probs = np.exp(test)
            sample_probs = sample_probs / np.sum(sample_probs)
            new_cnv = np.random.choice(self.cnv_support, p = sample_probs)
            cnv[exon] = new_cnv
        # sample new X_vect after updating with data
        X_vect = intensities.sample(self.X_priors, self.scale, self.data)
        gibbs_X[i] = X_vect
        
        log_probs = np.log(np.multiply(cnv, X_vect) / np.sum(np.multiply(cnv, X_vect)))
        likelihoods = np.sum(np.multiply(log_probs, self.data))

        return cnv, X_vect, likelihoods

    def joint_log_p(self, cnv, X_vect):
        log_probs_norm = np.sum(self.data) * -1 * np.log(np.sum(np.multiply(cnv, X_vect)))
        log_likelihood = (log_probs_norm +
                          np.sum(np.multiply(np.log(cnv), self.data) +
                                 np.multiply((self.data + (self.X_priors * self.scale) - 1), np.log(X_vect))))
