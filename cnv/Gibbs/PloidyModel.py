"""The main Gibbs Sampling Model """
from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution
from cnv.Targets.TargetCollection import TargetCollection

class PloidyModel(object):
    """This is the full statistical model and class that runs the Gibbs Sampling. It is responsible for taking a
    parameter set, a BAM file, and a TargetCollection and running a Gibbs sampler to determine the
    posterior probability of different ploidy states."""

    def __init__(self, cnv_support, X_priors, bamFileName, targets, parameterFileName):
        """Initialize the data model with its input arguments, we should load the parameters and
        calculate the coverage at each interval in the BAM file"""
        isinstance(targets, TargetCollection)
        # Fake code below
        data = TargetCollection.getData(bamFileName)
        self.intensities = IntensitiesDistribution(parameterFileName)

        self.cnv_support = cnv_support
        self.X_priors = X_priors
        self.ploidy = CopyNumberDistribution(cnv_support, X_priors, targets, data)

    def RunGibbsSampler(self, numIterations=10000):
        """ Run a Gibbs sampler for several iterations """
        self.gibbs_cnv_data = np.zeros((len(X_priors), iterations))
        self.gibbs_X = np.zeros((iterations, len(X_priors)))
        self.likelihoods = np.zeros(iterations)

        for i in xrange(0, numIterations):
            if (i+1) % (iterations / 20) == 0:
                print 'Finished {} iterations'.format(i)
            self.intensities.sample()
            cnv, X_vect, likelihoods = self.ploidy.sample(self.intensities)
            for exon in range(len(self.X_priors)):
                self.gibbs_cnv_data[exon, i] = cnv[exon]
            self.gibbs_X[i] = X_vect
            self.likelihoods[i] = likelihoods

    def OutputGibbsData(self, out_file_name, burn_in=1000, df_wanted=True):
        """Output a results file and a PDF of the posterior probabilities plot.
        Report on the posterior distribution obtained by the sampling procedure"""
        # get proportions using burn-in of 1000 iterations 
        gibbs_data_results = np.zeros((len(self.X_priors), len(self.cnv_support)))
        for index in range(len(self.X_priors)):
            # exclude samples before burn in and then take only every 100th sample to reduce autocorrelation
            gibbs_slice = self.gibbs_cnv_data[index][burn_in:][::100]
            gibbs_data_results[index] = np.bincount(gibbs_slice.astype(np.int64), 
                                                    minlength=len(self.cnv_support)+1)[1:]
        gibbs_data_results = gibbs_data_results / float(len(gibbs_slice))
    
        # return df for easier visualization
        if df_wanted:
            gibbs_df = pd.DataFrame(gibbs_data_results, columns =['copy_{}'.format(cnv) for cnv in self.cnv_support])
            gibbs_df['Exon'] = self.exon_labels
            return self.gibbs_cnv_data, self.gibbs_X, gibbs_data_results, self.likelihoods, gibbs_df

        return self.gibbs_cnv_data, self.gibbs_X, gibbs_data_results, self.likelihoods
