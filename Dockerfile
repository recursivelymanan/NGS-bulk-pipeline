FROM ubuntu:22.04

RUN apt-get update && \
    apt-get install -y wget bzip2

# Setup miniconda and add to PATH
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -u -p /miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /miniconda3/bin/conda init bash

ENV PATH=/miniconda3/bin/:$PATH

# Create conda env with necessary packages
RUN conda create -n nf bioconda::fastqc bioconda::star bioconda::multiqc -y

CMD ["/bin/bash"]