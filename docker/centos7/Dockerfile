FROM centos:7

RUN set -x && \
    yum --assumeyes install epel-release && \
    yum --assumeyes install centos-release-scl && \
    yum --assumeyes install \
    devtoolset-7-toolchain \
    gcc-c++ \
    gcc \
    zlib-devel \
    autoconf \
    automake \
    cmake3 \
    make \
    rpm-build \
    rpm && \
    yum --assumeyes clean all
