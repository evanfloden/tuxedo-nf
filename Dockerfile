FROM debian:jessie

MAINTAINER evan.floden@crg.eu

RUN apt-get update && apt-get install apt-utils build-essential vim wget unzip \
                                      libz-dev libncurses5-dev libncursesw5-dev \
                                      python --yes --force-yes

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O /opt/samtools-1.3.1.tar.bz2 && \
    cd opt && \
    tar xvfj /opt/samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    make && \
    make prefix=opt/samtools-1.3.1 install && \ 
    rm /opt/samtools-1.3.1.tar.bz2
ENV PATH opt/samtools-1.3.1:$PATH

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-source.zip -O /opt/hisat2-2.0.4-source.zip && \
    cd /opt && \
    unzip hisat2-2.0.4-source.zip && \
    cd hisat2-2.0.4 && \
    make && \ 
    rm /opt/hisat2-2.0.4-source.zip
ENV PATH opt/hisat2-2.0.4:$PATH

RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.tar.gz -O /opt/stringtie-1.3.0.tar.gz && \
    cd /opt && \
    tar xvfz stringtie-1.3.0.tar.gz && \
    cd stringtie-1.3.0 && \
    make release && \
    rm /opt/stringtie-1.3.0.tar.gz
ENV PATH opt/stringtie-1.3.0:$PATH

RUN echo "deb http://cran.rstudio.com/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list &&\
 apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 &&\
 apt-get update --fix-missing && \
 apt-get -y install r-base

RUN R -e 'source("http://bioconductor.org/biocLite.R"); library(BiocInstaller); biocLite("ballgown");'


