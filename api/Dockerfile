# Use ubuntu as a parent image - 18.04 has LTS
############################################################
######################## SUMO ##############################
############################################################
FROM ubuntu:18.04

# MAINTAINER David Schemm <david.schemm@tum.de>, based on docker-image 
# by Bo Gao <bogao@dcs.warwick.ac.uk> and Jakob Smretschnig <jakob.smretschnig@tum.de>
# LABEL Description="Dockerised SUMO (Simulation of Urban MObility)"

ENV SUMO_HOME /sumo

# # Install system dependencies
RUN apt-get -qq update && DEBIAN_FRONTEND=noninteractive apt-get -qq install \
    cmake g++ libxerces-c-dev libfox-1.6-dev libgdal-dev libproj-dev libgl2ps-dev swig git \
    ffmpeg \
    wget \
    sudo \
    xorg 
    # python3-dev \
    # python3-pip \
    # python3-rtree

RUN apt-get -qq install python3.7 curl python3-distutils && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3.7 get-pip.py

WORKDIR /
# # Download and build sumo
RUN git clone --recursive https://github.com/eclipse/sumo --verbose --progress && mkdir sumo/build/cmake-build
WORKDIR /sumo/build/cmake-build
RUN cmake ../.. && make -j8

############################################################
######################## API ###############################
############################################################
WORKDIR /
RUN apt-get -qq update && DEBIAN_FRONTEND=noninteractive apt-get -qq install python3.7-dev python3-rtree

RUN cp /usr/bin/python3.7 /usr/bin/python
# Setting up poetry
COPY poetry.lock /
COPY pyproject.toml .
RUN pip3 install poetry && \
    poetry config settings.virtualenvs.create false && \
    poetry install -v

COPY . /

EXPOSE 8000
EXPOSE 4040

CMD gunicorn app.main:app --timeout 500 --workers 5 -k uvicorn.workers.UvicornWorker -b 0.0.0.0

