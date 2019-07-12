FROM ubuntu:14.04

LABEL Description="Provide Ubuntu 14.04 build environment Anaquin project"

ENV DEBIAN_FRONTEND noninteractive

RUN \
        apt-get update --assume-yes && \
        apt-get install --assume-yes build-essential software-properties-common && \
        add-apt-repository ppa:ubuntu-toolchain-r/test -y && \
        apt-get update --assume-yes && \
        apt-get install --assume-yes --no-install-recommends \
                gcc-7 \
                g++-7 \
                gcc \
                g++ \
                zlib1g-dev \
                autoconf \
                automake \
                make \
                cmake3 && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/*

