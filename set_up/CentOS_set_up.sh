#!/usr/bin/env bash

###############################################
# install bcftools/samtools/SVPV dependencies #
# tested on fresh CentOS 7 image              #
###############################################

yum update -y && yum install -y \
gcc \
cmake \
git \
ImageMagick \
ncurses-devel \
python \
numpy \
tkinter \
zlib-devel \
epel-release \
asciidoc \
libxml2-devel \
xmlto

yum install -y R && yum clean all

cd /usr/local \
&& git clone https://github.com/samtools/htslib.git \
&& cd /usr/local/htslib \
&& make \
&& make lib-static \
&& make install

cd /usr/local \
&& git clone https://github.com/samtools/bcftools.git \
&& cd /usr/local/bcftools \
&& make \
&& make docs \
&& make install

cd /usr/local \
&& git clone https://github.com/samtools/samtools.git \
&& cd /usr/local/samtools \
&& make \
&& make install

