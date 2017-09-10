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

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O /opt/samtools-1.3.1.tar.bz2 && \
    cd opt && \
    tar xvfj /opt/samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    make && \
    make prefix=opt/samtools-1.3.1 install && \ 
    rm /opt/samtools-1.3.1.tar.bz2
ENV PATH /opt/samtools-1.3.1:$PATH


RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.0/sratoolkit.2.8.0-ubuntu64.tar.gz -O /opt/ncbi/sratoolkit.2.8.0-ubuntu64.tar.gz && \
    cd /opt/ncbi && tar zxf /opt/ncbi/sratoolkit.2.8.0-ubuntu64.tar.gz && \
    rm /opt/ncbi/sratoolkit.2.8.0-ubuntu64.tar.gz
ENV PATH /opt/ncbi/sratoolkit.2.8.0-ubuntu64/bin:$PATH

RUN groupadd -g 1001 sra_user
RUN useradd -d /home/sra_user -u 1001 -g 1001 -m -s /bin/bash sra_user
USER sra_user
ENV HOME /home/sra_user
RUN \
    cd /home/sra_user/ && \
    wget http://download.asperasoft.com/download/sw/connect/3.6.2/aspera-connect-3.6.2.117442-linux-64.tar.gz && \
    tar -xzvf aspera-connect-3.6.2.117442-linux-64.tar.gz && \
    sh aspera-connect-3.6.2.117442-linux-64.sh && \
    rm aspera-connect-3.6.2.117442-linux-64.tar.gz
USER root
ENV HOME /root/
ENV PATH /home/sra_user/.aspera/connect/bin/:$PATH
ENV ASCPKEY /home/sra_user/.aspera/connect/etc/asperaweb_id_dsa.openssh

RUN mkdir /ncbi_site && mkdir /ncbi_local && \
    vdb-config --root -s /repository/user/main/public/root=/ncbi_local && \
    echo '/repository/site/main/public/apps/file/volumes/flat = "files"' >> /root/.ncbi/user-settings.mkfg && \
    echo '/repository/site/main/public/apps/nakmer/volumes/nakmerFlat = "nannot"' >> /root/.ncbi/user-settings.mkfg && \
    echo '/repository/site/main/public/apps/nannot/volumes/nannotFlat = "nannot"' >> /root/.ncbi/user-settings.mkfg && \
    echo '/repository/site/main/public/apps/refseq/volumes/refseq = "refseq"' >> /root/.ncbi/user-settings.mkfg && \
    echo '/repository/site/main/public/apps/sra/volumes/sraFlat = "sra"' >> /root/.ncbi/user-settings.mkfg && \
    echo '/repository/site/main/public/apps/wgs/volumes/wgsFlat = "wgs"' >> /root/.ncbi/user-settings.mkfg && \
    vdb-config --root -s /repository/site/main/public/root=/ncbi_site

RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.9.8.tar.gz -O /opt/gffcompare-0.9.8.tar.gz && \
    cd /opt/ && \
    tar zxf /opt/gffcompare-0.9.8.tar.gz && \
    cd gffcompare-0.9.8 && \
    make && \
    chmod +x gffcompare
ENV PATH /opt/gffcompare-0.9.8:$PATH

RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    chmod 755 /FastQC/fastqc  && \
    ln -s /FastQC/fastqc bin/fastqc && \
    rm fastqc_v0.11.5.zip

