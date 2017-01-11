FROM ubuntu:16.04

MAINTAINER Jacob Munro j.munro@victorchang.edu.au

RUN apt-get update && apt-get install -y \
    cmake \
    git \
    python \
    python-dev \
    python-numpy \
    python-tk \
    graphicsmagick-imagemagick-compat \
    zlib1g-dev \
    && apt-get clean

# install htslib, samtools and bcftools
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install
RUN cd /opt \
    && git clone https://github.com/samtools/bcftools.git \
    && cd /opt/bcftools \
    && make \
    && make docs \
    && make install
RUN cd /opt \
    && git clone https://github.com/samtools/samtools.git \
    && cd /opt/samtools \
    && make \
    && make install

# install SVPV
RUN cd ~
    && git clone https://github.com/VCCRI/SVPV.git \

# add user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

# start bash
CMD ["/bin/bash"]
