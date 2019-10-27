# FROM ubuntu:18.04
# # FROM python:3.7-slim

# LABEL Description="Dockerised SUMO (Simulation of Urban MObility)"
# ENV SUMO_HOME /sumo

# # Install system dependencies
# RUN apt-get -qq update \ 
#     && DEBIAN_FRONTEND=noninteractive apt-get -qq install \
#     cmake python g++ libxerces-c-dev libfox-1.6-dev libgdal-dev libproj-dev libgl2ps-dev swig git \
#     ffmpeg python-dev \
#     wget \
#     python2.7 \
#     sudo \
#     xorg \
#     python-pip \
#     python-rtree

# # Download and build sumo
# RUN git clone --recursive https://github.com/eclipse/sumo --verbose --progress
# RUN mkdir sumo/build/cmake-build
# WORKDIR sumo/build/cmake-build
# RUN cmake ../..
# RUN make -j8

# ADD . /app
# WORKDIR /app

# RUN pip install -r requirements.txt 
# ENV QUART_APP emviz:app
# EXPOSE 5000
# # CMD quart run --host=0.0.0.0 --port=80
# ENTRYPOINT [ "python", "start_server.py" ]


# dockerfile

# Use ubuntu as a parent image - 18.04 has LTS
############################################################
######################## SUMO ##############################
############################################################
FROM ubuntu:18.04

# # MAINTAINER Jakob Smretschnig <jakob.smretschnig@tum.de>, based on docker-image by Bo Gao <bogao@dcs.warwick.ac.uk> 
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
WORKDIR sumo/build/cmake-build
RUN cmake ../.. && make -j8

# # Install python dependencies & Copy the RESTful service
# COPY requirements.txt /
# RUN mkdir -p /sumo 
# && pip install -r /requirements.txt

# COPY ./app /app
# FROM python:3.7
############################################################
######################## API ###############################
############################################################
WORKDIR /
RUN apt-get -qq update && DEBIAN_FRONTEND=noninteractive apt-get -qq install python3.7-dev python3-rtree
# RUN pip3 --version && python3.7 --version
# RUN echo "alias poetry='python3.7 $HOME/.poetry/bin/poetry'" >> ~/.bashrc
RUN cp /usr/bin/python3.7 /usr/bin/python
# Setting up poetry
COPY poetry.lock /
COPY pyproject.toml .
RUN pip3 install poetry && \
    poetry config settings.virtualenvs.create false && \
    poetry install -v


COPY . /

EXPOSE 8000
# EXPOSE 80
EXPOSE 4040
# alembic upgrade head && \
CMD gunicorn app.main:app --timeout 500 --workers 10 -k uvicorn.workers.UvicornWorker -b 0.0.0.0
# CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80"]




# for shared volume
# WORKDIR /ACPS

