FROM gitlab-registry.in2p3.fr/lisa/lisa-cde:0.1.0
RUN pip3 install --no-cache-dir snakemake==5.7.4

RUN mkdir /codes
COPY LISANode /codes/LISANode

RUN mkdir /codes/LDC
COPY ldc /codes/LDC/ldc
COPY data_generation /codes/LDC/data_generation
COPY setup.py requirements.txt /codes/LDC/
WORKDIR /codes/LDC
RUN pip3 install -r requirements.txt
RUN python3 setup.py build_liborbits 
RUN python3 setup.py install
#WORKDIR /codes/LDC/ldc/lisa/orbits/lib
#RUN make


WORKDIR /codes/LISANode
RUN pip3 install -e .


