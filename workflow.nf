params.output_dir = "./workflow_output"
params.input_dir = null
params.help = null

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

process alignment_setup {
    publishDir(
        path: "$params.output_dir/STAR/genome"
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

process alignment {
    publishDir(
        path: "$params.output_dir/STAR"
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
    
    """
}

workflow onlyQC {
    fqs = channel.from(file("$params.input_dir/*.fastq"))
    
    fastQC(fqs)
    multiQC(fastQC.out.html.collect(), fastQC.out.zip.collect())
}

workflow skipQC {
    paired_fqs = channel.fromFilePairs(file("$params.input_dir/*{1,2}.fastq"))
    genome_fasta = channel.from(file("$params.input_dir/*.fna"))
    genome_gtf = channel.from(file("$params.input_dir/*.gtf"))

    alignment_setup(genome_fasta, genome_gtf, input_dir)
    alignment(paired_fqs, alignment_setup.out.index.collect())
}

workflow {
    onlyQC()
    skipQC()
}
