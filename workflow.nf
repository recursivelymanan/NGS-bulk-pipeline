params.output_dir = "./workflow_output"
params.input_dir = null
params.help = null
params.download_reference_files = false

// Show help message
if (params.help == true) {
    print("To use the quality control pipeline, use the --input_dir parameter to indicate the path to the directory containing your .fastq files. Optionally, set an output directory using --output_dir")
}

// Throw error if no input_folder is provided
if (params.input_dir == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --input_dir")
}

println("Saving all workflow output to $params.output_dir")
println("Starting workflow...\n")

/*
Process for running FastQC quality control on given input files. Outputs HTML QC report to folder qc_output.
*/
process fastQC {
    publishDir(
        path: "$params.output_dir/qc_reports",
        mode: "copy"
    )

    input:
    path fqs

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    """
    fastqc $fqs
    """
}

/*
Process for combining previously generated HTML QC reports into one QC summary using MultiQC. Requires that the fastqc process created a folder named qc_output
*/
process multiQC {
    publishDir (
        path: "$params.output_dir/qc_reports",
        mode: "copy"
    )

    input:
    path fq_html
    path fq_zip

    output:
    path "multiqc_report.html"

    """
    multiqc .
    """
}

/*
Retrieve genomic sequence .fna and gene transfer format .gtf file for the human genome for downstream processes. 
*/
process retrieveFilesHuman {
    output:
    path "human_genome/ncbi_dataset/data/GCF_000001405.40/*.fna", emit: genome
    path "human_genome/ncbi_dataset/data/GCF_000001405.40/*.gtf", emit: gtf

    """
    datasets download genome accession GCF_000001405.40 --dehydrated --include genome,gtf --filename human_GRCh38_dataset.zip
    unzip human_GRCh38_dataset.zip -d human_genome
    datasets rehydrate --directory human_genome/
    """
}

/*
Rename the genome files retrieved by retrieveFilesHuman() in order to remove the long directory chain resulting from the ncbi download. 
*/
process renameGenomeFiles {
    publishDir (
        path: "$params.output_dir/ref",
        mode: "copy"
    )

    input:
    path genome
    path gtf

    output:
    path "*.fna", emit: genome
    path "*.gtf", emit: gtf

    """
    cp $genome genome.fna
    cp $gtf genome.gtf
    """
}

/*
Prepare the genome index directory to prepare for mapping, using STAR with --runMode genomeGenerate.
*/
process alignmentSetup {
    publishDir(
        path: "$params.output_dir/STAR/genome",
        mode: "copy"
    )

    input:
    path genome_fasta
    path genome_gtf

    output:
    path "STAR", emit: index

    """
    STAR --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles $genome_fasta
    --sjdbGTFfile $genome_gtf
    --sjdbOverhang 99
    """
}

/*
Map reads using STAR.
*/
process alignment {
    publishDir(
        path: "$params.output_dir/STAR",
        mode: "copy"
    )

    input:
    path paired_fqs
    path genome_index

    output:
    path "*.bam", emit: bam

    """
    STAR \
    --genomeDir $genome_index \
    --readFilesIn $paired_fqs \
    --outSAMtype BAM
    """
}

workflow onlyQC {
    fqs = channel.from(file("$params.input_dir/*.fastq"))
    
    fastQC(fqs)
    multiQC(fastQC.out.html.collect(), fastQC.out.zip.collect())
}

workflow skipQC {
    paired_fqs = channel.fromFilePairs(file("$params.input_dir/*{1,2}.fastq"))

    if (params.download_reference_files) {
        println("Downloading required genome files from NCBI...")
        retrieveFilesHuman()
        renameGenomeFiles(retrieveFilesHuman.out.genome.collect(), retrieveFilesHuman.out.gtf.collect())
        genome_fasta = retrieveFilesHuman.out.genome.collect()
        genome_gtf = retrieveFilesHuman.out.gtf.collect()
    }
    else {
        println("Using provided files...")
        genome_fasta = channel.from(file("$params.input_dir/*.fna"))
        genome_gtf = channel.from(file("$params.input_dir/*.gtf"))
    }
    alignmentSetup(genome_fasta, genome_gtf)
    //alignment(paired_fqs, alignment_setup.out.index.collect())
}

workflow {
    onlyQC()
    skipQC()
}
