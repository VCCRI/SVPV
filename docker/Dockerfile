FROM ubuntu:16.04

MAINTAINER Jacob Munro j.munro@victorchang.edu.au
# python-dev
RUN apt-get update && apt-get install -y \
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
    && git clone https://github.com/VCCRI/SVPV.git \

# add user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

# start bash
CMD ["/bin/bash"]

