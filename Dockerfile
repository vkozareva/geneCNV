# Source Image
FROM ubuntu:16.04

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        build-essential \
        gfortran \
        libatlas-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblapack-dev \
        liblzma-dev \
        libopenblas-dev \
        libssl-dev \
        python \
        python-dev \
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

# Set scipy building environment variables ATLAS, BLAS and LAPACK
ENV ATLAS=/usr/lib/libatlas.so.3
ENV BLAS=/usr/lib/libblas.so.3
ENV LAPACK=/usr/lib/liblapack.so.3

# Set libblas to use openblas alternative, recommended
RUN update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libblas.so.3

# Install geneCNV
RUN python setup.py install

# Define default command
CMD ["genecnv"]

# File Author / Maintainer
MAINTAINER Carlos Borroto <carlos@genepeeks.com>
