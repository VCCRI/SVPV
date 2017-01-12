#!/usr/bin/env bash
apt-get update && apt-get install -y \
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

