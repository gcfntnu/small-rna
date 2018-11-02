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

ENV TINI_VERSION v0.18.0
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
    bbmap \
    bowtie2 \
    bzip2 \
    atropos \
    fastp \
    samtools \
    fastx_toolkit \
    sortmerna \
    trimmomatic\
    ea-utils \
    fastqc 

