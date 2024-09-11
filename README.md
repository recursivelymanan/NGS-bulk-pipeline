# Bulk RNA-seq Processing Pipeline in Nextflow
Nextflow pipeline meant to streamline the RNAseq QC/processing workflow, including quality control, genome indexing, alignment and quantification. Written entirely in Nextflow and developed to be run directly with nextflow or through Docker. 

## Installing and running the workflow
You can run this workflow directly through Nextflow or using the provided Dockerfile.
### Running with Nextflow
1. [Install nextflow](https://www.nextflow.io/docs/latest/install.html)
2. Using the package manager of your choice, create an environment with the dependencies listed in the provided `env.yml`. Using miniconda, this would look like: `conda env create -f environment.yml`
3. Clone this repo to your local machine using `git clone https://github.com/recursivelymanan/NGS-bulk-pipeline`.
4. Navigate into the repo directory and run the workflow using `nextflow run workflow.nf --input_dir <Input directory>` and any additional options.

### Running with Docker
1. [Install docker](https://docs.docker.com/engine/install/)
2. Clone this repo to your local machine using `git clone https://github.com/recursivelymanan/NGS-bulk-pipeline`.
3. Build the Docker image using `docker build -t <image name> .`
4. Run the Docker image using `docker run -v /path/to/local/output:/app/workflow_output <image name>`

> **_NOTE:_**  The -v option in the docker run command mounts a local directory to the docker container so that output from nextflow workflow is not lost when the process ends and the container is closed. See [here](https://docs.docker.com/engine/storage/bind-mounts/) for more info. 

> **_NOTE:_**  If you need to modify the Docker runtime command, navigate into the `Dockerfile` and change the options in the last `CMD` line. By default, when running with Docker, all inputs (reads, genome files) are expected to be in an input directory named `data`, and all output will be directed to a folder named `workflow_output` inside of the container.

## Workflow options 
