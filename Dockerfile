FROM ubuntu:18.04

LABEL This Dockerfile is for HAMRLNC. It is maintained by Harry Li <harrli02@sas.upenn.edu> & Chosen Obih <chosenobih@arizona.edu>
ENV DEBIAN_FRONTEND=noninteractive

USER root

RUN apt-get update && apt-get install -y g++ \
		build-essential \
	   	make \
		git \
		libcurl4 \
		libcurl4-openssl-dev \
		libssl-dev \
		libncurses5-dev \
		libsodium-dev \
		libmariadb-dev \
		libbz2-dev \
		liblzma-dev \
		libssl-dev \
		zlib1g-dev \
		libxml2-dev \
		libfontconfig1-dev \
		libharfbuzz-dev \
		libfribidi-dev \
		libfreetype6-dev \
		libpng-dev \
		libtiff5-dev \
		libjpeg-dev \
		libpq-dev \
		openssl \
		default-jdk \
		lbzip2 \
		unzip \
		bzip2 \
		tzdata \ 
		perl \
		wget \
		bcftools \
		curl

RUN ldconfig
RUN apt-get install -y locales && locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

# Install additional system dependencies for R and R packages
RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository universe
RUN apt-get update && apt-get install -y \
#    software-properties-common \
    dirmngr \
	pkgconf \
    gnupg

# Manually add the GPG key for the CRAN repository
RUN apt-key adv --fetch-keys https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' #18.04
#RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/'   #22.04

#RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
#RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'


# Install R
#RUN apt-get update && apt-get install -y r-base
RUN apt-get update && apt-get install -y --no-install-recommends \
r-base

#r-base=4.1.2-1ubuntu2 \
#r-base-dev=4.1.2-1ubuntu2 \
#r-recommended=4.1.2-1ubuntu2

RUN echo "r-base hold" | dpkg --set-selections


# Install BiocManager in R
RUN R -e "options(repos = list(CRAN = 'http://cran.rstudio.com')); install.packages('BiocManager')"

# Install Biostrings package using BiocManager and log the output
RUN R -e "BiocManager::install('Biostrings', ask=FALSE)"
# Install additional R packages
RUN R -e "install.packages(c('dplyr', 'RPostgreSQL', 'httr', 'openssl', 'splitstackshape', 'getopt'))"
# additional R packages
RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages(c('reshape2', 'janitor', 'ggplot2', 'readr', 'tidyr', 'pheatmap'))"

# Downlaod and install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_23.11.0-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda
    #conda config --add channels r

# Conda packages
RUN conda install -y -c conda-forge libstdcxx-ng && \
	conda install cutadapt==4.9 -c bioconda -y && \
 	conda install bedops==2.4.41 -c bioconda -y && \
    conda install bedtools==2.31.1 -c bioconda -y && \
	conda install trim-galore==0.6.10 -c bioconda -y && \
	conda install star==2.7.10a -c bioconda -y && \
	conda install gffread==0.12.7 -y && \
 	conda install gffcompare==0.12.6 -c bioconda -y && \
	conda install subread==2.0.1 -c bioconda -y && \
	conda install stringtie==2.1.5 -c bioconda -y && \
	conda install bioawk==1.0 -c bioconda -y && \
 	conda install python -y && \
	conda install numpy==1.21.5 -y && \
	conda install pandas==1.3.5 -y && \
 	conda install pysam==0.22.1 -y && \
	conda install last==1454-0 -c bioconda -y && \
	conda install diamond==0.9.10 -c bioconda -y && \
	conda install transdecoder==5.5.0 -c bioconda -y && \
	conda install matplotlib-base -c conda-forge -y && \
	conda install fastp==0.23.4 -c bioconda -y && \
	conda install seqkit==2.8.2 -c bioconda -y

#Install biopython
RUN pip3 install biopython

# Required files
WORKDIR /
ENV BINPATH /usr/bin

# cpan
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm

# rFAM database
RUN apt-get install infernal
RUN mkdir Rfam && cd Rfam
RUN wget http://ftp.ebi.ac.uk/pub/databases/Rfam/14.10/Rfam.cm.gz && \
	gunzip Rfam.cm.gz && \
	wget http://ftp.ebi.ac.uk/pub/databases/Rfam/14.10/Rfam.clanin && \
	cmpress Rfam.cm && \
	cd ..

# CPC2
RUN wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz && \
	gzip -dc v1.0.1.tar.gz | tar xf - && \
	cd CPC2_standalone-1.0.1 && \
 	cd libs/libsvm && \
  	gzip -dc libsvm-3.18.tar.gz | tar xf - && \
   	cd libsvm-3.18 && \
    make clean && make 
WORKDIR /

## HAMR (python 3 compatible)
RUN git clone https://github.com/harrlol/HAMR.git
RUN chmod +x /HAMR/hamr.py
RUN cp /HAMR/hamr.py $BINPATH

# GATK (4.3.0.0)
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip && \
	unzip gatk-4.3.0.0.zip && rm gatk-4.3.0.0.zip
ENV PATH="/gatk-4.3.0.0:${PATH}"

# request
RUN apt-get update
RUN apt-get install -y libgdal-dev
RUN apt-get install python3-venv -y
RUN pip3 install requests

# jq
RUN apt-get install -y jq

# panther
RUN git clone https://github.com/pantherdb/pantherapi-pyclient.git && \
    cd pantherapi-pyclient && \
    python3 -m venv env && \
    . env/bin/activate && \
    pip3 install -r requirements.txt

WORKDIR /

RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz && \
 	tar -vxzf sratoolkit.tar.gz

RUN apt-get install bc -y

#HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
	tar -xvjf htslib-1.17.tar.bz2 && \
	cd htslib-1.17 && \
	./configure --prefix=/usr/bin && \
	make && \
	make install
WORKDIR /


ADD /scripts/*.py /scripts/
ADD /scripts/*.R /scripts/
RUN chmod +x /scripts/*.R && cp -r /scripts/ $BINPATH
ENV scripts /scripts
ADD util /util/
RUN chmod +x /util/*.pl && cp -r /util/ $BINPATH
ENV util /util

# Setting paths to all the softwares
ENV PATH /usr/bin/:$PATH
ENV PATH /HAMR/hamr.py:$PATH
ENV PATH /HAMR/:$PATH
ENV PATH  /sratoolkit.3.1.1-ubuntu64/bin:$PATH

# HAMRLNC wrapper script
ADD HAMRLNC.sh $BINPATH
RUN chmod +x $BINPATH/HAMRLNC.sh

ENTRYPOINT ["HAMRLNC.sh"]
