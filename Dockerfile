FROM debian:buster-slim

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get -y update \
    && apt-get -y --no-install-recommends install \
    # Install Python
    python3=3.7.3-1 \
    python3-pip=18.1-5 \
    # Install build packages
    gcc=4:8.3.0-1 \
    g++=4:8.3.0-1 \
    make=4.2.1-1.2 \
    python3.7-dev=3.7.3-2+deb10u1 \
    swig=3.0.12-2 \
    libhdf5-dev=1.10.4+repack-10 \
    libboost-dev=1.67.0.1 \
    # Install Python build packages
    libfftw3-dev=3.3.8-2 \
    libgsl-dev=2.5+dfsg-6 \
    && rm -rf /var/lib/apt/lists/* \
    # Install common Python packages
    && pip3 install --no-cache-dir setuptools==41.0.1 \
    && pip3 install --no-cache-dir \
        cycler==0.10.0 \
        cython==0.29.2 \
        h5py==2.10.0 \
        kiwisolver==1.2.0 \
        matplotlib==3.1.2 \
        numpy==1.18.1 \
        psd==1.5 \
        pyfftw==0.12.0 \
        pyparsing==2.4.7 \
        python-dateutil==2.8.1 \
        pyyaml==3.13 \
        scipy==1.4.1 \        
        six==1.14.0

RUN mkdir /codes
COPY LISANode /codes/LISANode

RUN mkdir /codes/LDC
COPY ldc /codes/LDC/ldc
COPY data_generation /codes/LDC/data_generation
COPY setup.py requirements.txt /codes/LDC/
WORKDIR /codes/LDC
RUN pip3 install -r requirements.txt
RUN python3 setup.py install
WORKDIR /codes/LDC/ldc/lisa/orbits/lib
RUN make


WORKDIR /codes/LISANode
RUN pip3 install -e .


