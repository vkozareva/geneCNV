# Install geneCNV on top of genepeeks-science.
FROM genepeeks-science

# Required by pysam.
USER root
RUN apt-get -y install liblzma-dev

COPY . geneCNV/
RUN chown -R genepeeks.genepeeks geneCNV
USER genepeeks
WORKDIR geneCNV
# Host's pyenv version isn't present; remove it.
RUN rm -f .python-version

RUN python setup.py install
RUN pip install nose
RUN ./runtests.sh

WORKDIR ..
RUN rm -rf geneCNV

# Copy model training data including coverage matrix and hyperparameters.
# For now just assume they were previously built in the training directory.
COPY training/training_bam.fofn \
     training/training_coverage_matrix_wanted.csv \
     training/training_targets_wanted.pickle \
     training/training_wanted_parameters.pickle \
     training/DMD_targets.bed \
     training/DMD_coverage_matrix.csv \
     training/DMD_targets.pickle \
     training/DMD_parameters.pickle \
     training/DMD_with_baseline_targets.bed \
     training/DMD_with_baseline_coverage_matrix.csv \
     training/DMD_with_baseline_targets.pickle \
     training/DMD_with_baseline_parameters.pickle \
     training/
USER root
RUN chown -R genepeeks.genepeeks training
USER genepeeks
WORKDIR training

ENTRYPOINT ["/home/genepeeks/.pyenv/shims/cnv", "evaluate-sample"]
