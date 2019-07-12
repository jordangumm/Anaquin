FROM centos:6

RUN set -x && \
    yum --assumeyes install epel-release && \
    yum --assumeyes install centos-release-scl && \
    yum --assumeyes install \
    devtoolset-7-toolchain \
    gcc-c++ \
    gcc \
    zlib-devel \
    autoconf268 \
    automake \
    cmake3 \
    make \
    rpm-build \
    rpm && \
    yum --assumeyes clean all

