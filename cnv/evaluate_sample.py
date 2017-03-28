from mando import main
from mando import command
import logging
import matplotlib
matplotlib.use(‘Agg’)
import matplotlib.pyplot as plt
from coverage_matrix import CoverageMatrix
from Gibbs.PloidyModel import PloidyModel
#pylint: disable=unused-argument

def visualize_mcmc(mcmc_df, title, outputFile):
    f, ax = plt.subplots(1, figsize=(20,10))
    bar_width = 1

    bar_l = [i for i in range(len(df['copy_1']))]
    tick_pos = [i + bar_width for i in bar_l]

    # loop through these later
    ax.bar(bar_l, df.copy_1, label='1 Copy', alpha=0.9, color='#019600', width=bar_width, edgecolor='white')

    ax.bar(bar_l, df.copy_2, bottom=df.copy_1, label='2 Copy', alpha=0.9, color='#3C5F5A', width=bar_width,
           edgecolor='white')

    ax.bar(bar_l, df.copy_3, bottom=[i+j for i,j in zip(df.copy_1, df.copy_2)], label='3 Copy',
           alpha=0.9, color='#219AD8', width=bar_width, edgecolor='white')

    plt.xticks(tick_pos, df['Target'])
    ax.set_ylabel("Probabilities")
    ax.set_xlabel("Targets")

    plt.xlim([min(tick_pos)-2*bar_width, max(tick_pos)+bar_width])
    plt.ylim(-0.1, 1.1)

    plt.setp(plt.gca().get_xticklabels(), rotation=60, horizontalalignment='right')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
               fancybox=True, shadow=True, ncol=3)
    plt.title(title)
    plt.savefig(outputFile)


@command('evaluate-sample')
def evaluate_sample(subjectBamfilePath, parametersFile, outputFile, include_covar=False):
    """Test for copy number variation in a given sample

    :param subjectID: The ID of the subject
    :param parametersFile: The files with the model parameters
    :param outputFile: Output file name.

    """
    logging.info("Running evaluate samples")


    # Read the parameters file.
    with open(parametersFile, 'w') as f:
        hln_parameters = pickle.load(f)

    targets = hln_parameters.targets
    # Parse subject bamfile
    subject_df = CoverageMatrix().create_coverage_matrix([subjectBamfilePath], targets)
    subject_id = subject_df['subject'][0]
    target_columns = [col for col in subject_df.columns if ('Target' in col) or ('Ex' in col)]
    subject_data = subject_df[target_columns].values.astype('float')

    cnv_support = [1, 2, 3]

    ploidy_model = PloidyModel(cnv_support, hln_parameters, data=subject_data, include_covar=include_covar)
    ploidy_model.RunMCMC()
    mcmc_posteriors, mcmc_df = ploidy_model.ReportMCMCData()

    visualize_mcmc(mcmc_df, 'Copy Number Posteriors for Subject {}'.format(subject_id), outputFile)




@command('train-model')
def train_model(subjectIDFile, outputFile):
    """Train a model that detects copy number variation.

    :param subjectIDFile: File with list of sample names
    :param outputFile: Output file name.

    """
    logging.info("Running sample training")


if __name__ ==  '__main__':
    main()
