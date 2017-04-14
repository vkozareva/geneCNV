# Install dmd on top of genepeeks-science.
FROM genepeeks-science

# Required by pysam.
USER root
RUN apt-get -y install liblzma-dev

COPY . dmd/
RUN chown -R genepeeks.genepeeks dmd
USER genepeeks
WORKDIR dmd
# Host's pyenv version isn't present; remove it.
RUN rm -f .python-version

RUN python setup.py install
RUN pip install nose
RUN ./runtests.sh

WORKDIR ..
RUN rm -rf dmd
