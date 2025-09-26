FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

ENV PYTHONEXE=python3
ENV HDF5_ROOT_DIR=/usr/lib/x86_64-linux-gnu/hdf5/openmpi


RUN apt-get update \
    && apt-get install -y --no-install-recommends \    
        git \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-wheel \
        python3-dev \
        build-essential \
        gcc \
        make \
        mpich \
        libboost-all-dev \
        libhdf5-mpich-dev \
        libhdf5-openmpi-dev \
        ffmpeg \
    && rm -rf /var/lib/apt/lists/*


RUN pip3 install -U --no-cache-dir \
    ipykernel==6.29.5 \
    ipython \
    h5py==3.10.0 \
    jupyter==1.0.0 \
    jupyterlab==4.1.4 \
    matplotlib==3.8.3 \
    numpy==1.26.4 \
    scipy==1.11.3 \ 
    pint \
    sphinx \
    ipympl==0.9.3 \
    ipywidgets==8.1.3
    


