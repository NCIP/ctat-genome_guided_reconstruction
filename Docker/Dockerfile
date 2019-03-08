FROM ubuntu:18.04
MAINTAINER bhaas@broadinstitute.org


RUN apt-get update && \
    apt-get install -y python git && \
	apt-get clean

ENV BIN /usr/local/bin
ENV SRC /usr/local/src

WORKDIR $SRC

RUN git clone https://github.com/NCIP/ctat-genome_guided_reconstruction.git && \
    mv ctat-genome_guided_reconstruction/* $BIN

RUN apt-get install -y wget

## STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.0e.tar.gz && \
    tar xvf 2.7.0e.tar.gz && \
    mv STAR-2.7.0e/bin/Linux_x86_64_static/STAR $BIN

## Stringtie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.5.Linux_x86_64.tar.gz && \
    tar xvf stringtie-1.3.5.Linux_x86_64.tar.gz && \
    mv stringtie-1.3.5.Linux_x86_64/stringtie $BIN


