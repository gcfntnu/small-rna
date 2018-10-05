FROM phusion/baseimage:latest

#CMD ["/sbin/my_init"]

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y \
    sed \
    wget \
    unzip \
    bzip2 
    
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

RUN apt-get -qq update

RUN conda update conda
RUN conda config --add channels defaults
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda update -y --all
RUN conda install -y bowtie \
    cutadapt \
    perl-compress-raw-zlib \
    perl-config-simple \
    perl-mime-base64 \
    perl-parallel-forkmanager \
    perl-xsloader \
    r-reshape2 \
    r-ggplot2 \
    r-latticeExtra \
    r-scales


ENV QUICKMIRSEQ_URL "https://downloads.sourceforge.net/project/quickmirseq/QuickMIRSeq.tgz"
RUN wget --quiet $QUICKMIRSEQ_URL -O- | tar xzv -C /opt
WORKDIR /opt/QuickMIRSeq
RUN chmod +x -R .
RUN cp merge-tRNA.pl merge_tRNA.pl
RUN sed -e 's/ftp:/http:/' -e 's/release-83/release-92/' database/mouse/download.config > database/mouse/download.config
RUN sed -e 's/ftp:/http:/' -e 's/release-83/release-92/' database/human/download.config > database/human/download.config
RUN sed -e 's/ftp:/http:/' -e 's/release-83/release-92/' database/rat/download.config > database/rat/download.config

RUN grep -rl '^#!/usr/bin/perl' ./ | xargs sed -i '1s:#!/usr/bin/perl:#!/usr/bin/env perl:'

ENV QuickMIRSeq /opt/QuickMIRSeq
ENV PATH $QuickMIRSeq:$PATH

