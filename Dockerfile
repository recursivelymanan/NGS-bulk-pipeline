FROM ubuntu:22.04

# Install jdk amd nextflow
RUN apt-get update && \
    apt-get install -y curl wget default-jdk && \ 
    curl -s https://get.nextflow.io | bash && \
    chmod +x nextflow && \
    mv nextflow /usr/bin

# Install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -u -p /miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    /miniconda3/bin/conda init bash

ENV PATH=/miniconda3/bin/:$PATH

# Create conda env with necessary packages
RUN conda create -n nf bioconda::fastqc bioconda::star bioconda::multiqc -y && \
    conda init && \
    conda activate nf

# Bring nextflow workflow files and test data into container
COPY workflow.nf nextflow.config /app/
COPY data/ /app/data/

WORKDIR /app

# Run the nextflow executable when the container starts. Parameters defined at runtime
ENTRYPOINT ["nextflow"]