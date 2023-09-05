FROM ubuntu:18.04

LABEL This Dockerfile is for hamrbox. It is maintained by Harry Li <harrli02@sas.upenn.edu> & Chosen Obih <chosenobih@arizona.edu>
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
		python3.6-dev \
        python3-pip \
        apt-utils \
		tzdata \
		perl \
		wget \
		bcftools \
		curl

RUN ldconfig
RUN apt-get install -y locales && locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

# Downlaod and install conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-py27_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    #wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_23.5.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH /opt/conda/bin:$PATH

# Conda channels
RUN conda config --add channels conda-forge && \
	conda config --add channels bioconda

# Conda packages
RUN conda install cutadapt==1.9.1 -y && \
	conda install sra-tools==3.0.5 -y && \
	conda install trim-galore==0.6.10 -y && \
	conda install STAR==2.7.10b -y && \
	conda install bedtools==2.31.0 -y && \
	conda install samtools==1.17 -y && \
	conda install python -y && \
	conda install tophat -y && \
	conda install bowtie2 -y

# Required files
WORKDIR /

# R libraries
RUN apt-get -y install ca-certificates software-properties-common gnupg2 gnupg1 gnupg && \
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

## HAMR (now working under Python 3)
ENV BINPATH /usr/bin
RUN git clone https://github.com/harrlol/HAMR && \
	chmod +x /HAMR/hamr.py && cp /HAMR/hamr.py $BINPATH && \
	cp -R /HAMR/models /usr/bin/hamr_models
ENV HAMR_MODELS_PATH=/usr/bin/hamr_models

# GATK (4.3.0.0)
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip && \
	unzip gatk-4.3.0.0.zip && rm gatk-4.3.0.0.zip
ENV PATH="/gatk-4.3.0.0:${PATH}"

# pip3
RUN apt-get update
RUN apt-get install -y libgdal-dev
RUN pip3 install --upgrade pip
RUN apt-get install python3-venv -y

# panther
RUN git clone https://github.com/pantherdb/pantherapi-pyclient.git && \
    cd pantherapi-pyclient && \
    python3 -m venv env && \
    . env/bin/activate && \
    pip3 install -r requirements.txt

# Setting paths to all the softwares
ENV PATH /usr/bin/:$PATH
ENV PATH /HAMR/hamr.py:$PATH
ENV PATH /HAMR/:$PATH

# Caching the sra data
RUN vdb-config --root -s /repository/user/cache-disabled="true"

# pamlinc wrapper script
ADD hamrbox.sh $BINPATH
RUN chmod +x $BINPATH/hamrbox.sh

ENTRYPOINT ["hamrbox.sh"]
