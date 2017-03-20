"""The main Gibbs Sampling Model """

import os
import numpy as np
import pandas as pd


from IntensitiesDistribution import IntensitiesDistribution
from CopyNumberDistribution import CopyNumberDistribution
from cnv.Targets.TargetCollection import TargetCollection

class PloidyModel(object):
    """This is the full statistical model and class that runs the Gibbs Sampling. It is responsible for taking a
    parameter set, a BAM file, and a TargetCollection and running a Gibbs sampler to determine the
    posterior probability of different ploidy states."""

    def __init__(self, bamFileName, targets, parameterFileName):
        """Initialize the data model with its input arguments.
        Load the parameters and calculate the coverage at each interval in the BAM file."""
        # Fake code below
        isinstance(targets, TargetCollection)
        data = targets.getData(bamFileName)
        self.intensities = IntensitiesDistribution(parameterFileName)
        self.cnv_support = [1, 2, 3]
        self.ploidy = CopyNumberDistribution(targets, data)

    @staticmethod
    def reshape_coverage_df(df, include_stats=False, groupby='subject', subject_droplist=None, df_counts_wanted=False):
        """Reshape data frames so that exons are observations (rows) and subjects are variables (columns).
        Make sure datetimes have been converted to datetime objects before using."""
        df_grouped = df.groupby([groupby]).sum()
        df_norm = df_grouped.div(df_grouped.sum(axis=1), axis=0)
        df_norm = df_norm.transpose().reset_index()
        df_norm.rename(columns={'index': 'Exon'}, inplace=True)

        if subject_droplist:
            for subject in subject_droplist:
                df_norm.drop(subject, axis=1, inplace=True)
        if include_stats:
            df_norm['Mean'] = df_norm.mean(axis=1)
            df_norm['SD'] = df_norm.std(axis=1)
        if df_counts_wanted:
            return df_norm, df_grouped
        else:
            return df_norm

    @staticmethod
    def compute_X_probs(exon_data_dir):
        #TODO: Remove any data loading code from here

        # get full dataset and also subsets based on gender, sequencer, etc
        coverage_df = pd.read_csv(os.path.join(exon_data_dir, 'coverage_matrix.csv'), header=0, index_col=0)
        coverage_df.is_rerun.values = False
        # remove rerun data and coding region counts
        coverage_df = coverage_df[coverage_df.is_rerun == False]
        coverage_df.drop(['is_rerun'], axis=1, inplace=True)
        coverage_df.index.name = None
        coverage_df.date_modified = pd.to_datetime(coverage_df.date_modified, unit='s')
        coverage_df['date'] = coverage_df.date_modified.dt.date
        coverage_df_RMA = coverage_df[coverage_df.subject.str.contains('FRMR')]

        # checking dates of RMA samples -- all RMA samples seem to have been run within a day of each other in June 2016
        # in fact, this is just the last time the bams were modified, the samples were run earlier in the year, see below
        RMA_dates = coverage_df_RMA.date.unique()

        # Use RMA samples for initial intensity vector--note that all RMA individuals used the M1 mixin panel.
        # This is only the females in RMA.
        by_flow = coverage_df_RMA.groupby(['subject', 'flow_cell_id']).sum()
        by_flow['tsid_ratio'] = by_flow.TSID_only / (by_flow.TSID_only + by_flow.TSO_only)
        by_flow.reset_index(inplace=True)
        subjects38 = by_flow[(by_flow.tsid_ratio < 0.39) & (by_flow.tsid_ratio > 0.37)]['subject']
        RMA_subset = coverage_df_RMA[coverage_df_RMA.subject.isin(subjects38)]
        gibbs_columns = ['subject'] + [column for column in RMA_subset.columns if 'Ex' in column]
        rma38 = RMA_subset[gibbs_columns]
        rma38_norm = PloidyModel.reshape_coverage_df(rma38, include_stats=True)
        return np.array(rma38_norm.Mean)

    def RunGibbsSampler(self, cnv=None, numIterations=10000):
        self.gibbs_cnv_data = np.zeros((len(self.intensities.intensities), numIterations))
        self.gibbs_X = np.zeros((numIterations, len(self.intensities.intensities)))
        self.likelihoods = np.zeros(numIterations)

        for i in xrange(0, numIterations):
            if (i+1) % (numIterations / 20) == 0:
                print 'Finished {} iterations'.format(i)
            self.intensities.sample(self.ploidy)
            likelihood = self.ploidy.sample(self.intensities)
            for exon in range(len(self.X_priors)):
                self.gibbs_cnv_data[exon, i] = cnv[exon]
            self.gibbs_X[i] = self.intensities.intensities
            self.likelihoods[i] = likelihood

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
