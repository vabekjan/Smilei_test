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
        bash-completion \
        # Qt/X runtime for Spyder
        # libgl1 \
        # libx11-6 libxext6 libxrender1 \
        # libxcb1 libxi6 libsm6 \
        # libxkbcommon0 libxkbcommon-x11-0 \
        # libdbus-1-3 xauth \
        libgl1 \
        libx11-6 libxext6 libxrender1 \
        libx11-xcb1 libxcb1 libxi6 libsm6 libice6 libxtst6 \
        libxkbcommon0 libxkbcommon-x11-0 libdbus-1-3 xauth \
        libglib2.0-0 libfontconfig1 libfreetype6 \
        libxcb-render0 libxcb-render-util0 \
        libxcb-icccm4 libxcb-image0 libxcb-keysyms1 \
        libxcb-randr0 libxcb-shape0 libxcb-xfixes0 \
        libxcb-xinerama0 libxcb-xkb1 libxcb-cursor0 libxcb-xinput0 \
    && rm -rf /var/lib/apt/lists/*

RUN echo 'if [ -f /etc/bash_completion ]; then . /etc/bash_completion; fi' >> /root/.bashrc


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
    ipywidgets==8.1.3 \
    jupyterlab-h5web==12.5.0 \
    spyder==6.0.0 \
    spyder-kernels==3.0.0 \
    PyQt6
    
# 3. Make JupyterLab terminals start Bash (login shell) instead of sh
RUN jupyter server --generate-config && \
    printf "\nc.ServerApp.terminado_settings = {'shell_command': ['/bin/bash','-l']}\n" \
    >> /root/.jupyter/jupyter_server_config.py

ENV SHELL=/bin/bash
ENV QT_X11_NO_MITSHM=1
ENV QT_QPA_PLATFORM=xcb