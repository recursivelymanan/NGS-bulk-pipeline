# Bulk RNA-seq Processing Pipeline in Nextflow

Nextflow pipeline meant to streamline the RNAseq QC/processing workflow, including quality control, genome indexing, alignment and quantification. Written entirely in Nextflow and developed to be run directly with nextflow or through Docker.

## Installing and running the workflow

You can run this workflow directly through Nextflow or using the provided Dockerfile.

### Running with Nextflow

1. [Install nextflow](https://www.nextflow.io/docs/latest/install.html)
2. Using the package manager of your choice, create an environment with the dependencies listed in the provided `env.yml`. Using miniconda, this would look like: `conda env create -f env.yml`
3. Clone this repo to your local machine using `git clone https://github.com/recursivelymanan/NGS-bulk-pipeline`.
4. Navigate into the repo directory and run the workflow using `nextflow run workflow.nf --inputDir <input directory>` and any additional options.

### Running with Docker

1. [Install docker](https://docs.docker.com/engine/install/)
2. Clone this repo to your local machine using `git clone https://github.com/recursivelymanan/NGS-bulk-pipeline`.
3. Build the Docker image using `docker build -t <image name> .`
4. Run the Docker image using `docker run -v /path/to/local/output:/app/workflow_output <image name>`

> **_NOTE:_** The -v option in the docker run command mounts a local directory to the docker container so that output from nextflow workflow is not lost when the process ends and the container is closed. See [here](https://docs.docker.com/engine/storage/bind-mounts/) for more info.

> **_NOTE:_** If you need to modify the Docker runtime command, navigate into the `Dockerfile` and change the options in the last `CMD` line. By default, when running with Docker, all inputs (reads, genome files) are expected to be in an input directory named `data`, and all output will be directed to a folder named `workflow_output` inside of the container.

## Workflow options

#### Synopsis

`nextflow run workflow.nf --inputDir <input dir> [options]`

#### Options

| Option                     | Default           | Description                                                                                                                                   |
| -------------------------- | ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `--help`                   | N/A               | Display the help message.                                                                                                                     |
| `--outputDir <output dir>` | `workflow_output` | Specify the output directory.                                                                                                                 |
| `--downloadReferenceFiles` | N/A               | Downloaded reference files from NCBI. If this option is not selected, required genome files are assumed to be located in the input directory. |
| `--aligner <aligner>`      | `hisat2`          | Pick which program to use for alignment, currently only HISAT2 is supported.                                                                  |
| `--paired`                 | N/A               | Denotes that experiment is paired-end. If not selected, it is assumed that the experiment is single-end reads.                                |

#### Custom Arguments

Each main program used in the workflow has an associated parameter with which you can add custom arguments. Please refer to the documentation for each program to ensure arguments are being passed properly.

| Program       | Parameter |
| ------------- | --------- |
| FastQC        | --fqc     |
| MultiQC       | --mqc     |
| HISAT2        | --hs2     |
| featureCounts | --fc      |

> **_NOTE:_** Do not pass custom arguments related to the format of reads (single vs. paired end), as these options are taken care of based on the --paired parameter.

##### Example:

Using MultiQC as an example, if you wanted to change the title of the MultiQC report you could do so with the MultiQC --title parameter as follows:

`nextflow run workflow.nf --inputDir data --mqc "--title CustomTitle"`

## Planned Additions

- Functionality for use with the STAR aligner
- A `nextflow-aws.config` file to enable easier use with AWS
- ~~Ability to pass on more arguments into the various tools to allow for more customizable workflows~~ Added 9-12-24
