FROM ubuntu:trusty
MAINTAINER Matei David <matei.david.at.oicr.on.ca>

RUN apt-get update && \
    apt-get install -y \
        software-properties-common && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y \
        build-essential \
        g++-4.9 \
        libhdf5-dev \
        libboost1.55-dev \
        libboost-python1.55-dev \
        python2.7-minimal \
        python-setuptools \
        python-virtualenv

ENV TZ=${TZ}
RUN ln -snf /usr/share/zoneinfo/${TZ} /etc/localtime && echo ${TZ} > /etc/timezone

RUN useradd --create-home --uid ${USER_ID} ${USER_NAME}
USER ${USER_NAME}

VOLUME /data
WORKDIR /data
