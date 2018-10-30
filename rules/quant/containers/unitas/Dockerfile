FROM phusion/baseimage:latest

CMD ["/sbin/my_init"]

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y \
    wget \
    unzip \
    bzip2 \
    perl \
    libarchive-extract-perl \
    cpanminus \
    build-essential \
    gcc-multilib \
    dos2unix
    
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN cpanm install LWP::Simple \
    File::Path \
    Archive::Zip \
    File::Copy \
    File::Basename \
    Getopt::Long \
    Archive::Extract \
    Scalar::Util
 
ENV UNITAS_VERSION 1.6.1
RUN wget --quiet http://www.smallrnagroup.uni-mainz.de/software/unitas_${UNITAS_VERSION}.zip
RUN unzip unitas_${UNITAS_VERSION}.zip -d /opt
RUN rm unitas_${UNITAS_VERSION}.zip
WORKDIR /opt/unitas
RUN dos2unix unitas_${UNITAS_VERSION}.pl
RUN chmod +755 unitas_${UNITAS_VERSION}.pl
RUN ln -s unitas_${UNITAS_VERSION}.pl unitas.pl
ENV PATH /opt/unitas/:$PATH

