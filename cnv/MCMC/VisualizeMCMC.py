import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class VisualizeMCMC(object):
    """A class which visualizes the results of MCMC sampling"""
    def __init__(self, cnv_support, targets, copy_posteriors, intensity_posteriors=None):
        self.cnv_support = cnv_support
        # first three are historically used colors -- can be changed to match rest of color map easily
        self.colors =  ['tomato', '#019600', '#3C5F5A', '#219AD8'] + [(0, i / 20.0, 0, 1) for i in range(len(cnv_support)-4)]
        self.targets = targets
        self.copy_posteriors = copy_posteriors
        self.intensity_posteriors = intensity_posteriors  # may want to visualize intensity posteriors in some way later

    def visualize_copy_numbers(self, title, outputFile=None):
        """Generate stacked bar chart for copy number posteriors, given an array of posteriors from MCMC sampling."""
        f, ax = plt.subplots(1, figsize=(20,10))
        bar_width = 1

        bar_l = [i for i in range(len(self.copy_posteriors))]
        tick_pos = [i + bar_width for i in bar_l]

        # loop through and add bar to stack for each copy number in support
        for i in range(len(self.cnv_support)):
            args = {'label': '{} Copy'.format(self.cnv_support[i]),
                    'alpha': 0.9,
                    'color': self.colors[i],
                    'width': bar_width,
                    'edgecolor': 'white'}

            # build on other bars if not first
            if i != 0:
                args['bottom'] = np.sum(self.copy_posteriors[:, :i], axis=1)
            ax.bar(bar_l, self.copy_posteriors[:,i], **args)

        plt.xticks(tick_pos, self.targets)
        ax.set_ylabel('Probabilities')
        ax.set_xlabel('Targets')

        plt.xlim([min(tick_pos)-2*bar_width, max(tick_pos)+bar_width])
        plt.ylim(-0.02, 1.02)

        plt.setp(plt.gca().get_xticklabels(), rotation=60, horizontalalignment='right')
        plt.legend(loc='upper center', bbox_to_anchor=(1.05, 0.8),
                   fancybox=True, shadow=True, ncol=1)
        plt.title(title)
        if outputFile:
            plt.savefig(outputFile)
        else:
            plt.show()