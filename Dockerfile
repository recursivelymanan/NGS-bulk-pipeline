FROM continuumio/miniconda3

# Install jdk and nextflow
RUN apt-get update && \
    apt-get install -y curl wget default-jdk unzip && \ 
    curl -s https://get.nextflow.io | bash && \
    chmod +x nextflow && \
    mv nextflow /usr/bin

# Create conda env with necessary packages
RUN conda create -n nf bioconda::fastqc bioconda::star bioconda::multiqc bioconda::subread conda-forge::ncbi-datasets-cli -y

# Bring nextflow workflow files and test data into container
COPY workflow.nf nextflow.config /app/
COPY data/ /app/data/

WORKDIR /app

# Run the nextflow executable when the container starts. Parameters defined at runtime
ENTRYPOINT ["conda"]

CMD ["run", "-n", "nf", "--live-stream", "nextflow", "run", "workflow.nf", "--input_dir", "data"]
