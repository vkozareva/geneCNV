# Source Image
FROM ubuntu:16.04

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        python \
        python-pip \
        zlib1g-dev

# apt clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Update pip and setuptools
RUN pip install -U pip setuptools

# Pre-install requirements.
# Not sure why but this avoids lengthy compilation.
COPY . geneCNV/
WORKDIR geneCNV
RUN pip install -U -r requirements.txt
RUN python setup.py install
WORKDIR ..
RUN rm -rf geneCNV

# Define default command
CMD ["genecnv"]

# File Author / Maintainer
MAINTAINER Carlos Borroto <carlos@genepeeks.com>
