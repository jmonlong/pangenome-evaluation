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

## minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 > /build/minimap2-2.17_x64-linux.tar.bz2 \
        && tar -jxvf /build/minimap2-2.17_x64-linux.tar.bz2 \
        && rm /build/minimap2-2.17_x64-linux.tar.bz2

ENV PATH=/build/minimap2-2.17_x64-linux/:$PATH

## seqtk
RUN git clone https://github.com/lh3/seqtk.git && \
        cd seqtk && \
        make

ENV PATH /build/seqtk/:$PATH

## minigraph
RUN wget --no-check-certificate -O v0.10.tar.gz https://github.com/lh3/minigraph/archive/v0.10.tar.gz && \
        tar -xzvf v0.10.tar.gz && \
        cd minigraph-0.10 && \
        make

ENV PATH /build/minigraph-0.10/:$PATH

## bcftools
RUN wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
        tar -xjf bcftools-1.10.2.tar.bz2 && \
        cd bcftools-1.10.2 && \
        ./configure && \
        make && \
        make install && \
        cd .. && \
        rm -rf bcftools-1.10.2 bcftools-1.10.2.tar.bz2

###### vg
##

# fetch the desired git revision of vg
ARG vg_git_revision=9e1b22c4fb506aa6206018238ab43d44aa29af02

RUN git clone https://github.com/vgteam/vg.git /vg

WORKDIR /vg

RUN git fetch --tags origin && \
        git checkout "$vg_git_revision" && \
        git submodule update --init --recursive

# If we're trying to build from a non-recursively-cloned repo, go get the
# submodules.
RUN bash -c "[[ -e deps/sdsl-lite/CMakeLists.txt ]] || git submodule update --init --recursive"

RUN make get-deps && \
        . ./source_me.sh && \
        env && \
        make include/vg_git_version.hpp && \
        make -j $((THREADS < $(nproc) ? THREADS : $(nproc))) && \
        make static && \
        strip -d bin/vg

ENV PATH /vg/bin:$PATH

RUN rm /vg/bin/bgzip /vg/bin/tabix
##
######

WORKDIR /build

## fpa
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs -o rust.sh && \
        sh rust.sh -y --no-modify-path

ENV PATH /root/.cargo/bin:$PATH

RUN git clone https://github.com/natir/fpa.git && \
        cd fpa && \
        git checkout v0.5.1 && \
        cargo build && \
        cargo install --path .

## seqwish
RUN apt-get update && \
        apt-get -qqy install zlib1g zlib1g-dev libomp-dev && \
        apt-get -qqy install build-essential software-properties-common && \
        add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
        apt-get update > /dev/null && \
        apt-get -qqy install gcc-snapshot && \
        apt-get update > /dev/null && \
        apt-get -qqy install gcc-8 g++-8 && \
        update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 80 --slave /usr/bin/g++ g++ /usr/bin/g++-8 && \
        apt-get -qqy install cmake git

ARG seqwish_git_revision=39d3aa32252c4739e2a5c5abb19cc581223eab8f

RUN git clone --recursive https://github.com/ekg/seqwish.git && \
        cd seqwish && \
        git checkout "$seqwish_git_revision" && \
        cmake -H. -Bbuild && \
        cmake --build build -- 

RUN apt-get -qy autoremove

ENV PATH /build/seqwish/bin:$PATH

## odgi
ARG odgi_git_revision=1029a5a69a3c0f68e07ad0afca9aa54bb5bdc33f

RUN git clone --recursive https://github.com/vgteam/odgi.git && \
        cd odgi && \
        git fetch --tags origin && \
        git checkout "$odgi_git_revision" && \
        git submodule update --init --recursive && \
        cmake -H. -Bbuild && \
        cmake --build build -- -j $THREADS

ENV PATH /build/odgi/bin:$PATH

## smoothxg
ARG smoothxg_git_revision=ed8321be39a8430b50e1a0d3a275e711da60bf40

RUN git clone --recursive https://github.com/ekg/smoothxg.git && \
        cd smoothxg && \
        git fetch --tags origin && \
        git checkout "$smoothxg_git_revision" && \
        git submodule update --init --recursive && \
        cmake -H. -Bbuild && \
        cmake --build build -- -j $THREADS

ENV PATH /build/smoothxg/bin:$PATH

## GraphAligner
RUN wget -O Miniconda-latest-Linux.sh https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
        bash Miniconda-latest-Linux.sh -b -p /miniconda

ENV PATH /miniconda/bin:$PATH

SHELL ["/bin/bash", "-c"] 

RUN conda init bash && source ~/.bashrc && conda update -n base -c defaults conda

RUN git clone --recursive https://github.com/maickrau/GraphAligner && \
        cd GraphAligner && \
        git fetch --tags origin && \
        git checkout "v1.0.12" && \
        git submodule update --init --recursive && \
        conda env create -f CondaEnvironment.yml && \
        source activate GraphAligner && \
        make bin/GraphAligner

ENV PATH /build/GraphAligner/bin:$PATH

## edyeet
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        libgsl-dev \
        && rm -rf /var/lib/apt/lists/*

ARG edyeet_git_revision=8f827e824016757842d33034c9364598976cede0

RUN git clone --recursive https://github.com/ekg/edyeet.git && \
        cd edyeet && \
        git fetch --tags origin && \
        git checkout "$edyeet_git_revision" && \
        ./bootstrap.sh && \
        ./configure && \
        make && \
        make install


ENV VG_COMMIT $vg_git_revision
ENV SEQWISH_COMMIT $seqwish_git_revision
ENV SMOOTHXG_COMMIT $smoothxg_git_revision
ENV EDYEET_COMMIT $edyeet_git_revision

## SVanalyzer
# RUN conda config --add channels bioconda && conda install svanalyzer && conda update --all

RUN apt-get update \
        && apt-get install -y --no-install-recommends \
        singularity-container \
        && rm -rf /var/lib/apt/lists/*

## launcher script for ksub in the home directory
WORKDIR /home
ADD run.sh /home
