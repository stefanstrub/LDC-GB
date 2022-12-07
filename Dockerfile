FROM gitlab-registry.in2p3.fr/lisa/lisa-cde:0.3.0
USER root
RUN apt-get -y update \
&& apt-get -y --no-install-recommends install wget liblapacke-dev nano
RUN pip3 install --no-cache-dir mistune==0.8.4 sphinx sphinx_rtd_theme mkdocs-material m2r2 tqdm

RUN mkdir /codes
COPY FastEMRIWaveforms /codes/FastEMRIWaveforms
RUN rm -r /codes/FastEMRIWaveforms/docs
COPY constants /codes/constants

WORKDIR /codes/constants
RUN python3 setup.py install install_headers clean --all

RUN mkdir /codes/LDC
COPY README.md /codes/LDC/README.md
COPY ldc /codes/LDC/ldc
COPY test /codes/LDC/test
COPY doc /codes/LDC/doc
COPY data_generation /codes/LDC/data_generation
COPY setup.py requirements.txt requirements.jupyter.txt /codes/LDC/
WORKDIR /codes/LDC
RUN pip3 install -r requirements.txt
RUN python3 setup.py install --with-fastGB --with-imrphenomD --with-fastAK 
RUN python3 setup.py clean --all
#WORKDIR /codes/LDC/ldc/lisa/orbits/lib
#RUN make

WORKDIR /codes/FastEMRIWaveforms
RUN sed -i 's/"hdf5_hl"/"hdf5_hl", "lapacke"/' setup.py
RUN LDSHARED='h5c++ -shlib -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions  -Wl,-z,relro -g -fwrapv -O2    ' CC='h5c++ -shlib' CXX='h5c++ -shlib' python3 setup.py install clean --all
WORKDIR /codes/LDC/test
RUN python3 install_emri_files.py

USER lisauser
