FROM ubuntu:18.04
MAINTAINER Roberto Semeraro <robe.semeraro@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.6 python3-pip python3-setuptools python3-dev

#Git
RUN apt-get install -y git

#Samtools
RUN apt-get install -y samtools && \
    apt-get install -y bcftools && \
    apt-get install -y tabix

#XomeBlender
WORKDIR /root

RUN git clone --single-branch -b 3.1 https://github.com/rsemeraro/XomeBlender.git && \
    cd XomeBlender && \
    git checkout 3.1 &&  \
    python3 setup.py install
