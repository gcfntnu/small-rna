FROM ubuntu:18.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

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

RUN conda install -y python=2.7 pip wheel bowtie samtools cutadapt viennarna pigz matplotlib pandas biopython numpy=1.15.1 scipy scikit-learn reportlab forgi

WORKDIR /opt
ADD https://api.github.com/repos/flatberg/miRge/git/refs/heads/master version.json
RUN git clone https://github.com/flatberg/miRge.git
RUN cd miRge && \
    rm ./dist/* && \
    python setup.py bdist_wheel && \
    pip install ./dist/mirge-2.0-py2-none-any.whl &&\
    rm -rf /opt/miRge