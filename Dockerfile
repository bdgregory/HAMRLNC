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
		libmariadb-client-lgpl-dev \
		libbz2-dev \
		liblzma-dev \
		libssl-dev \
		zlib1g-dev \
		libcurl4-openssl-dev \ 
		openssl \
		default-jdk \
		lbzip2 \
		unzip \
		bzip2 \
		# python3 \
		# python3-pip \
		# python-matplotlib \
		# python-numpy \
  		# python-pandas \
	 	# python-pysam \
		tzdata \ 
		perl \
		wget \
		bcftools \
		# python2.7 \
		# python2.7-dev \
		curl

RUN ldconfig
RUN apt-get install -y locales && locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

# Install system dependencies for R and R packages
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    gnupg \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libpq-dev

# Manually add the GPG key for the CRAN repository
RUN apt-key adv --fetch-keys https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'

# Install R
RUN apt-get update && apt-get install -y r-base

# Install BiocManager in R
RUN R -e "options(repos = list(CRAN = 'http://cran.rstudio.com')); install.packages('BiocManager')"

# Install Biostrings package using BiocManager and log the output
RUN R -e "BiocManager::install('Biostrings', ask=FALSE)" \
    && R -e "packageVersion('Biostrings')"
# Install additional R packages
RUN R -e "install.packages(c('dplyr', 'RPostgreSQL', 'httr', 'openssl', 'splitstackshape', 'getopt'))"

# Downlaod and install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_23.11.0-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels r

# Conda packages
RUN conda install cutadapt -c bioconda -y && \
	# conda install hisat2==2.2.1 -c bioconda -y && \
	# conda install bowtie2==2.2.5 -c bioconda -y && \
 	conda install bedops==2.4.41 -c bioconda -y && \
    	conda install bedtools==2.31.1 -c bioconda -y && \
    	conda install htslib==1.16 -c bioconda -y && \
    	conda install sra-tools==3.0.10 -c bioconda -y && \
	conda install trim-galore==0.6.10 -c bioconda -y && \
	conda install bedtools==2.31.0 -c bioconda -y && \
	conda install samtools==1.16.1 -c bioconda -y && \
	conda install star==2.7.10a -c bioconda -y && \
	conda install gffread -y && \
 	conda install gffcompare==0.12.6 -c bioconda -y && \
	conda install subread==2.0.1 -c bioconda -y && \
	conda install stringtie==2.1.5 -c bioconda -y && \
	conda install bioawk==1.0 -c bioconda -y && \
 	conda install python -y && \
	conda install numpy -y && \
	conda install pandas -y && \
 	conda install pysam -y && \
	conda install last==1454-0 -c bioconda -y && \
	conda install diamond==0.9.10 -c bioconda -y && \
	conda install transdecoder==5.5.0 -c bioconda -y && \
	conda install matplotlib-base -c conda-forge -y

#Install biopython
RUN pip3 install biopython

# Required files
WORKDIR /
ENV BINPATH /usr/bin

## evolinc-part-I
#RUN git clone https://github.com/chosenobih/Evolinc-I.git
#RUN cp -R /Evolinc-I /evolinc_docker
#WORKDIR /evolinc_docker

# Cufflinks
# RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
#     tar -zxvf cufflinks-2.2.1.Linux_x86_64.tar.gz && \
#     rm -rf cufflinks-2.2.1.Linux_x86_64.tar.gz
# ENV PATH="/cufflinks-2.2.1.Linux_x86_64:${PATH}"

# cpan
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm

# Remove the existing symbolic link (if it exists), create a symbolic link to make 'python' refer to 'python3',
# And make Python 3 the default Python 
# RUN rm /usr/bin/python && \
#     ln -sf /usr/bin/python3 /usr/bin/python && \
#     apt-get update && apt-get install -y python3 && \
#     update-alternatives --install /usr/bin/python python /usr/bin/python3 1

# Uniprot database
#ADD https://github.com/iPlantCollaborativeOpenSource/docker-builds/releases/download/evolinc-I/uniprot_sprot.dmnd.gz /evolinc_docker/
#RUN gzip -d /evolinc_docker/uniprot_sprot.dmnd.gz && \
#	chmod +r /evolinc_docker/uniprot_sprot.dmnd

# rFAM database
#ADD https://de.cyverse.org/dl/d/12EF1A2F-B9FC-456D-8CD9-9F87197CACF2/rFAM_sequences.fasta /evolinc_docker/
RUN apt-get install infernal
RUN mkdir Rfam && cd Rfam
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz && \
	gunzip Rfam.cm.gz && \
	wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin && \
	cmpress Rfam.cm && \
	cd ..

RUN R -e "install.packages('tidyr')"

# CPC2
RUN wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz && \
	gzip -dc v1.0.1.tar.gz | tar xf - && \
	cd CPC2_standalone-1.0.1 && \
 	cd libs/libsvm && \
  	gzip -dc libsvm-3.18.tar.gz | tar xf - && \
   	cd libsvm-3.18 && \
    	make clean && make 

# evolinc-part-I wrapper script
#RUN chmod +x /evolinc_docker/evolinc-part-I.sh && cp /evolinc_docker/evolinc-part-I.sh $BINPATH

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

# RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
# 	tar -vxzf sratoolkit.tar.gz

RUN apt-get install bc -y

# fix to pysam not found
# RUN pip3 wheel pysam && pip3 install pysam*.whl

# additional R packages
RUN R -e "install.packages(c('reshape2', 'janitor', 'ggplot2', 'readr', 'pheatmap'))"


ADD /scripts/*.R /scripts/
RUN chmod +x /scripts/*.R && cp -r /scripts/ $BINPATH
ENV scripts /scripts
ADD util /util/
RUN chmod +x /util/*.pl && cp -r /util/ $BINPATH
ENV util /util

# Setting paths to all the softwares
ENV PATH /evolinc_docker/cufflinks-2.2.1.Linux_x86_64/:$PATH
#ENV PATH /evolinc_docker/bin/:$PATH
#ENV PATH /evolinc_docker/CPC2-beta/bin/:$PATH
#ENV PATH /evolinc_docker/:$PATH
ENV PATH /usr/bin/:$PATH
ENV PATH /HAMR/hamr.py:$PATH
ENV PATH /HAMR/:$PATH
# ENV PATH /sratoolkit.3.0.10-ubuntu64/bin:$PATH

# HAMRLNC wrapper script
ADD HAMRLNC.sh $BINPATH
RUN chmod +x $BINPATH/HAMRLNC.sh

ENTRYPOINT ["HAMRLNC.sh"]
