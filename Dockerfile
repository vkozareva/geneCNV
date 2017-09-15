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

# Update pip and setuptools
RUN pip install -U pip
RUN pip install -U setuptools

# Install requirements that are failing to install with setup.py
RUN pip install numpy==1.11.1

# Install requirements that are failing to install with pip and setup.py
RUN apt-get install -y  \
        python-matplotlib

# apt clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set scipy building environment variables ATLAS, BLAS and LAPACK
ENV ATLAS=/usr/lib/libatlas.so.3
ENV BLAS=/usr/lib/libblas.so.3
ENV LAPACK=/usr/lib/liblapack.so.3

# Set libblas to use openblas alternative, recommended
RUN update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libblas.so.3

# Install geneCNV
COPY . /geneCNV
RUN cd /geneCNV && \
    python setup.py install

# Define default command
CMD ["genecnv"]

# File Author / Maintainer
MAINTAINER Carlos Borroto <carlos@genepeeks.com>
