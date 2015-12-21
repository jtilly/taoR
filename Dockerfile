## start with the latest Ubuntu image
FROM ubuntu:latest

## Maintainer
MAINTAINER "Jan Tilly" jantilly@gmail.com

## Set working directory
WORKDIR /root

## Set path to petsc
ENV PETSC_DIR /usr/local/lib/python2.7/dist-packages/petsc

## Remain current
RUN apt-get update -qq \
    && apt-get dist-upgrade -y

## Add compilers, make, and valgrind
RUN apt-get install -y --no-install-recommends \
        g++ \
        gcc \
        gfortran \
        make \
        valgrind

## Add lapack / blas 
RUN apt-get install -y --no-install-recommends \
        libblas3gf \
        libblas-dev  \
        liblapack3gf \
        liblapack-dev

## Install python-dev and python-pip
RUN apt-get install -y --no-install-recommends \
        python-dev \
        python-pip

## Install petsc
RUN pip install petsc --allow-external petsc
    
## Install R
RUN apt-get install -y --no-install-recommends apt-transport-https \
    && echo "deb https://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list \ 
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 \
    && apt-get update -qq \
    && apt-get install -y --no-install-recommends r-base-dev \
    && apt-get install -y --no-install-recommends r-cran-rcpp

## Clone taoR and install it
RUN cd ~ \
    && git clone https://github.com/jtilly/taoR \
    && R CMD build taoR \
    && export PKG_FILE_NAME=$(ls -1t *.tar.gz | head -n 1) \
    && R CMD INSTALL $PKG_FILE_NAME

## Clean-up
RUN rm -rf /tmp/* && rm -rf /var/lib/apt/lists/*
