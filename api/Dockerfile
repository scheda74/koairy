FROM ubuntu:18.04
# FROM python:3.7-slim

LABEL Description="Dockerised SUMO (Simulation of Urban MObility)"
ENV SUMO_HOME /sumo

# Install system dependencies
RUN apt-get -qq update \ 
    && DEBIAN_FRONTEND=noninteractive apt-get -qq install \
    cmake python g++ libxerces-c-dev libfox-1.6-dev libgdal-dev libproj-dev libgl2ps-dev swig git \
    ffmpeg python-dev \
    wget \
    python2.7 \
    sudo \
    xorg \
    python-pip \
    python-rtree

# Download and build sumo
RUN git clone --recursive https://github.com/eclipse/sumo --verbose --progress
RUN mkdir sumo/build/cmake-build
WORKDIR sumo/build/cmake-build
RUN cmake ../..
RUN make -j8

ADD . /app
WORKDIR /app

RUN pip install -r requirements.txt 
ENV QUART_APP emviz:app
EXPOSE 5000
# CMD quart run --host=0.0.0.0 --port=80
ENTRYPOINT [ "python", "emviz.py" ]