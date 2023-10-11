FROM ubuntu:18.04

LABEL This Dockerfile is for HAMRLINC. It is maintained by Harry Li <harrli02@sas.upenn.edu> & Chosen Obih <chosenobih@arizona.edu>
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
		python2.7 \
		python2.7-dev \
		python-pip \
		python-matplotlib \
		python-numpy \
       	python-pandas \
		tzdata \ 
		perl \
		wget \
		bcftools \
		curl

RUN ldconfig
RUN apt-get install -y locales && locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

#Install biopython
RUN pip2 install biopython==1.76

# Downlaod and install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Conda channels
RUN conda config --add channels conda-forge && \
	conda config --add channels bioconda

# Conda packages
RUN conda install cutadapt==1.9.1 -c bioconda -y && \
	conda install sra-tools==3.0.5 -c bioconda -y && \
	conda install trim-galore==0.6.10 -c bioconda -y && \
	conda install bedtools==2.31.0 -c bioconda -y && \
	conda install samtools==1.17 -c bioconda -y && \
	conda install biopython==1.76 -c anaconda -y && \
	conda install star==2.7.10a -c bioconda -y && \
	conda install gffread==0.12.1 -c bioconda -y && \
	conda install subread==2.0.1 -c bioconda -y && \
	conda install cufflinks==2.2.1 -c bioconda -y && \
	conda install stringtie==2.1.5 -c bioconda -y && \
	conda install bowtie2==2.2.5 -c bioconda -y && \
	conda install tophat==2.1.1 -c bioconda -y && \
	conda install numpy -y && \
	conda install pandas -y && \
	conda install last==1454-0 -c bioconda -y && \
	conda install diamond==0.9.10 -c bioconda && \
	conda install matplotlib-base -c conda-forge -y && \
	conda install python -y

# Required files
WORKDIR /

## evolinc-part-I
RUN git clone https://github.com/chosenobih/Evolinc-I.git
RUN cp -R /Evolinc-I /evolinc_docker
ENV BINPATH /usr/bin
WORKDIR /evolinc_docker

# Transdecoder
RUN wget -O- https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz | tar xzvf -

# Bedops tool
RUN wget -O- https://github.com/bedops/bedops/releases/download/v2.4.16/bedops_linux_x86_64-v2.4.16.tar.bz2 | tar jxvf -

# cpan
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm URI/Escape.pm


# R libraries
RUN apt-get update && apt-get upgrade -y && \
	apt-get -y install ca-certificates software-properties-common gnupg2 gnupg1 gnupg && \
	#apt-key adv --keyserver hkp://pgp.mit.edu --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" && \
	apt-get install -y r-base && \
	Rscript -e 'install.packages("openssl", dependencies = TRUE,  repos="http://cran.rstudio.com/")' && \
	Rscript -e 'install.packages("splitstackshape", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
	Rscript -e 'install.packages("dplyr", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
	Rscript -e 'install.packages("tidyr", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("data.table", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("ggplot2", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("janitor", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("pheatmap", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("readr", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("reshape2", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("stringr", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
    Rscript -e 'install.packages("viridislite", dependencies = TRUE, repos="http://cran.rstudio.com/");' && \
	Rscript -e 'install.packages("getopt", dependencies = TRUE, repos="http://cran.rstudio.com/");'

# Uniprot database
ADD https://github.com/iPlantCollaborativeOpenSource/docker-builds/releases/download/evolinc-I/uniprot_sprot.dmnd.gz /evolinc_docker/
RUN gzip -d /evolinc_docker/uniprot_sprot.dmnd.gz && \
	chmod +r /evolinc_docker/uniprot_sprot.dmnd

# rFAM database
ADD https://de.cyverse.org/dl/d/12EF1A2F-B9FC-456D-8CD9-9F87197CACF2/rFAM_sequences.fasta /evolinc_docker/

# CPC2
WORKDIR /evolinc_docker/CPC2-beta/libs/libsvm/
RUN tar xvf libsvm-3.22.tar.gz
WORKDIR libsvm-3.22
RUN make clean && make
WORKDIR /

# evolinc-part-I wrapper script
RUN chmod +x /evolinc_docker/evolinc-part-I.sh && cp /evolinc_docker/evolinc-part-I.sh $BINPATH

## HAMR (now working under Python 3)
ENV BINPATH /usr/bin
RUN git clone https://github.com/harrlol/HAMR && \
	chmod +x /HAMR/hamr.py && cp /HAMR/hamr.py $BINPATH && \
	cp -R /HAMR/models /usr/bin/hamr_models
ENV hamr_model=/usr/bin/hamr_models

# GATK (4.3.0.0)
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip && \
	unzip gatk-4.3.0.0.zip && rm gatk-4.3.0.0.zip
ENV PATH="/gatk-4.3.0.0:${PATH}"

# pip3
RUN apt-get update
RUN apt-get install -y libgdal-dev
RUN apt-get install -y python3-pip
RUN pip3 install --upgrade pip
RUN apt-get install python3-venv -y

# panther
RUN git clone https://github.com/pantherdb/pantherapi-pyclient.git && \
    cd pantherapi-pyclient && \
    python3 -m venv env && \
    . env/bin/activate && \
    pip3 install -r requirements.txt

#HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
	tar -xvjf htslib-1.17.tar.bz2 && \
	cd htslib-1.17 && \
	./configure --prefix=/usr/bin && \
	make && \
	make install
WORKDIR /

RUN apt-get install bc -y

ADD /scripts/*.R /scripts/
RUN chmod +x /scripts/*.R && cp -r /scripts/ $BINPATH
ENV scripts /scripts
ADD util /util/
RUN chmod +x /util/*.pl && cp -r /util/ $BINPATH
ENV util /util

# Setting paths to all the softwares
ENV PATH /evolinc_docker/TransDecoder-2.0.1/:$PATH
ENV PATH /evolinc_docker/ncbi-blast-2.4.0+/bin/:$PATH
ENV PATH /evolinc_docker/bin/:$PATH
ENV PATH /evolinc_docker/CPC2-beta/bin/:$PATH
ENV PATH /evolinc_docker/:$PATH
ENV PATH /usr/bin/:$PATH
ENV PATH /HAMR/hamr.py:$PATH
ENV PATH /HAMR/:$PATH

# Caching the sra data
RUN vdb-config --root -s /repository/user/cache-disabled="true"

# HAMRLINC wrapper script
ADD HAMRLINC.sh $BINPATH
RUN chmod +x $BINPATH/HAMRLINC.sh

ENTRYPOINT ["HAMRLINC.sh"]
