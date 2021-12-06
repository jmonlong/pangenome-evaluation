FROM ubuntu:18.04

MAINTAINER jmonlong@ucsc.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true
ARG THREADS=4

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        wget \
        curl \
        less \
        python3 \
        python3-pip \
        python3-setuptools \
        python \
        python-dev \
        python-pip \
        python-setuptools \
        python3-dev \
        gcc \ 
        bcftools \
        samtools \
        tabix \
        tzdata \
        make \
        pigz \
        gawk \
        graphviz \
        imagemagick \
        gfortran-8 \
        bzip2 \
        git \
        sudo \
        pkg-config \
        libxml2-dev libssl-dev libmariadbclient-dev libcurl4-openssl-dev \ 
        apt-transport-https software-properties-common dirmngr gpg-agent \ 
        && rm -rf /var/lib/apt/lists/*

ENV TZ=America/Los_Angeles

WORKDIR /build

## GNU time
RUN wget https://ftp.gnu.org/gnu/time/time-1.9.tar.gz && \
        tar -xzvf time-1.9.tar.gz && \
        cd time-1.9 && \
        ./configure && \
        make && \
        make install

## AWS cli and modules for Snakemake
RUN pip3 install --upgrade pip

RUN pip3 install --no-cache-dir requests awscli snakemake boto3 pandas numpy

## R and pandoc
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
        && apt-get update \
        && apt-get install -y r-base r-base-dev

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz')"

RUN R -e "install.packages(c('dplyr', 'knitr', 'ggplot2', 'rmarkdown', 'DT', 'ggrepel', 'tidyr', 'BiocManager', 'plotly', 'DiagrammeR', 'visNetwork', 'shiny'))"

RUN R -e "BiocManager::install(c('VariantAnnotation'))"

RUN apt-get update \
    && wget https://github.com/jgm/pandoc/releases/download/2.10.1/pandoc-2.10.1-1-amd64.deb \
    && dpkg -i pandoc-2.10.1-1-amd64.deb \
    && apt-get install -f \
    && rm pandoc-2.10.1-1-amd64.deb

## seqtk
RUN git clone https://github.com/lh3/seqtk.git && \
        cd seqtk && \
        make

ENV PATH /build/seqtk/:$PATH

## bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
        tar -xjf bcftools-1.10.2.tar.bz2 && \
        cd bcftools-1.10.2 && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -rf bcftools-1.10.2 bcftools-1.10.2.tar.bz2

## Singularity
RUN apt-get update && \
	apt-get install -y build-essential \
	libseccomp-dev pkg-config squashfs-tools cryptsetup

RUN wget --quiet https://golang.org/dl/go1.15.5.linux-amd64.tar.gz && \
	tar -C /usr/local -xzf go1.15.5.linux-amd64.tar.gz

ENV GOPATH=/build/go
ENV PATH=$PATH:/usr/local/go/bin

RUN curl -sSfL https://raw.githubusercontent.com/golangci/golangci-lint/master/install.sh | sh -s -- -b $(go env GOPATH)/bin

WORKDIR /build/go/src/github.com/hpcng

RUN git clone https://github.com/hpcng/singularity.git && \
	cd singularity && \
	git checkout v3.8.0 && \
	./mconfig && \
	cd ./builddir && \
	make && \
	sudo make install

## add screen and nano for interacive mode
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	screen nano \
        && rm -rf /var/lib/apt/lists/*

ENV PS1 '${debian_chroot:+($debian_chroot)}\u \W\$ '

## launcher script for ksub in the home directory
WORKDIR /home
ADD run.sh /home
