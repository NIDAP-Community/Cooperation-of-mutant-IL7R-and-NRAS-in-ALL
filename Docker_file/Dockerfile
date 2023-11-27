# Built Docker Image with Ubuntu 20.04 System
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Set system overall permission in the image
RUN chmod -R ugo+rx /root \
  && chmod -R ugo+rx /opt \
  && chmod -R ugo+rx /tmp \
  && umask u=rwx,g=rwx,o=rx	\
  && echo "TMPDIR=/mnt" > /root/.Renviron
  
RUN apt-get update && apt-get install -y software-properties-common

# Install System Libary and Packages
RUN apt-get update \
    && add-apt-repository --enable-source --yes "ppa:marutter/rrutter4.0" \
    && add-apt-repository --enable-source --yes "ppa:c2d4u.team/c2d4u4.0+" \
    && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libxml2-dev \
        libssl-dev \
        libpng-dev \
        libhdf5-dev \
        libquadmath0 \
        libtiff5-dev \
        libjpeg-dev \
        libfreetype6-dev \
        libgfortran5 \
        libgmp-dev \
        libmpc-dev \
        libopenblas0-pthread \
        libgeos-dev \
        libcairo2 \
        libcairo2-dev \
        libpango1.0-dev \
        libfftw3-dev \
        locales \
        vim-tiny \
        wget \
        ca-certificates \
        fonts-texgyre \
        curl \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Set System Encoding
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh \
    && bash Miniconda3-py37_4.8.2-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py37_4.8.2-Linux-x86_64.sh

# Add Miniconda Path to System
ENV PATH="/root/miniconda3/bin:${PATH}"

# Config Conda and install mamba
RUN conda config --append channels bioconda \
    && conda config --append channels conda-forge \
    && conda install -c conda-forge mamba

# Initiate Conda Bash Terminal
RUN conda init bash

RUN conda install -c conda-forge conda-build
# Create local conda channel
RUN mkdir /local_channel
RUN mkdir /local_channel/linux-64
COPY *.tar.bz2 /local_channel/linux-64/
RUN conda index /local_channel
RUN conda config --add channels file:///local_channel

# Build Conda Environment from environment.yml
COPY environment.yml .
RUN apt-get update \
    && apt-get install -y libarchive-dev
RUN mamba env update -f environment.yml

# Ask conda to activate the built environment when starting
ARG built_env_name="$(grep -E '^name:' environment.yml | sed 's/name: *//')"
RUN echo "source activate $built_env_name" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

CMD ["bash"]
