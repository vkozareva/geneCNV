from mando import main
from mando import command
import logging
#pylint: disable=unused-argument

@command('evaluate-sample')
def evaluate_sample(subjectID, parametersFile, outputFile):
    """Test for copy number variation in a given sample

    :param subjectID: The ID of the subject
    :param parametersFile: The files with the model parameters
    :param outputFile: Output file name.

    """
    logging.info("Running evaluate samples")

@command('train-model')
def train_model(subjectIDFile, outputFile):
    """Train a model that detects copy number variation.

    :param subjectIDFile: File with list of sample names
    :param outputFile: Output file name.

    """
    logging.info("Running sample training")


if __name__ ==  '__main__':
    main()
