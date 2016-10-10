FROM debian:jessie

MAINTAINER evan.floden@crg.eu

RUN apt-get update && apt-get install apt-utils build-essential vim wget unzip \
                                      libz-dev libncurses5-dev libncursesw5-dev \
                                      libmagic-dev libhdf5-dev git openjdk-7-jdk \
                                      libcurl4-openssl-dev libxml2-dev libssl-dev \
                                      python --yes --force-yes

RUN mkdir /opt/ncbi && \
    cd /opt/ncbi && \
    git clone https://github.com/ncbi/ngs.git && \
    cd ngs/ngs-sdk && \
    ./configure && \
    make && make install && \
    cd /opt/ncbi && \
    git clone https://github.com/ncbi/ncbi-vdb && \
    cd ncbi-vdb && \
    ./configure && \
    make && make install

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-source.zip -O /opt/hisat2-2.0.4-source.zip && \
    cd /opt && \
    unzip hisat2-2.0.4-source.zip && \
    cd hisat2-2.0.4 && \
    make USE_SRA=1 NCBI_NGS_DIR=/usr/local/ngs/ngs-sdk/ NCBI_VDB_DIR=/usr/local/ncbi/ncbi-vdb/ && \
    rm /opt/hisat2-2.0.4-source.zip
ENV PATH /opt/hisat2-2.0.4:$PATH

RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.0.tar.gz -O /opt/stringtie-1.3.0.tar.gz && \
    cd /opt && \
    tar xvfz stringtie-1.3.0.tar.gz && \
    cd stringtie-1.3.0 && \
    make release && \
    rm /opt/stringtie-1.3.0.tar.gz
ENV PATH /opt/stringtie-1.3.0:$PATH

RUN echo "deb http://cran.rstudio.com/bin/linux/debian jessie-cran3/" >>  /etc/apt/sources.list &&\
 apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 &&\
 apt-get update --fix-missing && \
 apt-get -y install r-base

RUN R -e 'source("http://bioconductor.org/biocLite.R"); library(BiocInstaller); biocLite("ballgown"); install.packages("devtools", repos="https://cloud.r-project.org"); devtools::install_github("alyssafrazee/RSkittleBrewer"); biocLite("genefilter"); install.packages("dplyr", repos="https://cloud.r-project.org")'




