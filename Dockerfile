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
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libssl-dev \
        python \
        python-matplotlib \
        python-nose \
        python-numpy \
        python-pandas \
        python-pip \
        python-scipy \
        zlib1g-dev

# apt clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install geneCNV
COPY . /geneCNV
RUN cd /geneCNV && \
    python setup.py install


# Define default command
CMD ["cnv"]

# File Author / Maintainer
MAINTAINER Carlos Borroto <carlos@genepeeks.com>
