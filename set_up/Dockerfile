FROM ubuntu:16.04

# Note: Installing SVPV in a Docker container will essentially prohibit GUI use

MAINTAINER Jacob Munro j.munro@victorchang.edu.au

# install prerequisites
RUN apt-get update -y && apt-get install -y \
    cmake \
    git \
    graphicsmagick-imagemagick-compat \
    libncurses-dev \
    python \
    python-numpy \
    python-tk \
    zlib1g-dev \
    && apt-get install --no-install-recommends -y \
    asciidoc \
    libxml2-utils \
    xmlto \
    r-base-core \
    && apt-get clean

# install htslib, samtools and bcftools
RUN cd /usr/local \
    && git clone https://github.com/samtools/htslib.git \
    && cd /usr/local/htslib \
    && make \
    && make lib-static \
    && make install
RUN cd /usr/local \
    && git clone https://github.com/samtools/bcftools.git \
    && cd /usr/local/bcftools \
    && make \
    && make docs \
    && make install
RUN cd /usr/local \
    && git clone https://github.com/samtools/samtools.git \
    && cd /usr/local/samtools \
    && make \
    && make install

# install SVPV
RUN cd /usr/local \
    && git clone https://github.com/VCCRI/SVPV.git

# add user
RUN groupadd -r -g 1000 svpv_usr && useradd -r -g svpv_usr -u 1000 -m svpv_usr
USER svpv_usr

# start bash
CMD ["/bin/bash"]
