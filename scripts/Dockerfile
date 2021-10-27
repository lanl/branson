## Ref: https://github.com/lanl/Draco/regression/Dockerfile

FROM kinetictheory/draco-travis-ci

# Use ubuntu if building from scratch
#FROM ubuntu:latest

MAINTAINER KineticTheory "https://github.com/KineticTheory"

# See draco/.travis-run-tests.sh for some instructions.

## Environment needed by 'docker build' ----------------------------------------

## for apt to be noninteractive
## https://stackoverflow.com/questions/8671308/non-interactive-method-for-dpkg-reconfigure-tzdata
## https://spack.readthedocs.io/en/latest/workflows.html?highlight=docker
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true
ENV SPACK_ROOT=/vendors/spack
ENV DRACO_TPL="cmake@3.11.4 random123@1.09 openmpi@3.1.1 netlib-lapack@3.8.0 metis@5.1.0"
ENV FORCE_UNSAFE_CONFIGURE=1
ENV DISTRO=bionic
ENV CLANG_FORMAT_VER=6.0
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1

## Update packages on the raw Ubuntu image -------------------------------------------

RUN sed -i '/DISTRO/d' /etc/apt/sources.list
# try to eliminate warning about "mesg: ttyname failed: Inappropriate ioctl for device"
RUN sed -i 's/mesg/tty -s \&\& mesg/' /root/.profile
RUN cat /root/.profile

## preesed tzdata, update package index, upgrade packages and install needed software
RUN echo "tzdata tzdata/Areas select US" > /tmp/preseed.txt; \
    echo "tzdata tzdata/Zones/US select Mountain" >> /tmp/preseed.txt; \
    debconf-set-selections /tmp/preseed.txt && \
    apt-get update && \
    apt-get install -y tzdata

## Basic admin tools
RUN apt-get -y install apt-utils autoconf python software-properties-common

## Basic developer tools
RUN apt-get -y install build-essential ca-certificates coreutils ccache curl doxygen environment-modules gcc-7 g++-7 gfortran-7 git grace graphviz python-pip tar tcl tk unzip vim wget
# RUN apg-get upgrade
RUN if ! test -f /etc/profile.d/modules.sh; then \
      echo "source /usr/share/modules/init/bash" > /etc/profile.d/modules.sh; \
    fi

## LLVM tools like clang-format
## Note: we can't use variables in the add-apt-repository commmand as this
##       creates an invalid entry in /etc/apt/sources.list
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - ; \
    add-apt-repository 'deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-6.0 main' ; \
    apt-get update; \
    apt-get install -y clang-format-${CLANG_FORMAT_VER}; \
    export PATH=$PATH:/usr/bin

## code cov plugin...
RUN python -m pip install --upgrade pip
RUN python -m pip install codecov

## SPACK -----------------------------------------------------------------------------

# install/setup spack
RUN mkdir -p $SPACK_ROOT
RUN curl -s -L https://api.github.com/repos/spack/spack/tarball | tar xzC $SPACK_ROOT --strip 1
# note: if you wish to change default settings, add files alongside
#       the Dockerfile with your desired settings. Then uncomment this line
#COPY packages.yaml modules.yaml $SPACK_ROOT/etc/spack/
RUN if ! test -f /etc/profile.d/spack.sh; then \
      echo "source $SPACK_ROOT/share/spack/setup-env.sh" > /etc/profile.d/spack.sh; \
    fi

## Provide some TPLs
RUN export PATH=$SPACK_ROOT/bin:$PATH && spack install ${DRACO_TPL} && spack clean -a

# image run hook: the -l will make sure /etc/profile.d/*.sh environments are loaded
CMD /bin/bash -l
