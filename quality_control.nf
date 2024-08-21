params.output_dir = "./workflow_output"
params.input_dir = null
params.help = null
params.genome_fasta = null
params.genome_gtf = null

// Show help message
if (params.help == true) {
    print("To use the quality control pipeline, use the --input_dir parameter to indicate the path to the directory containing your .fastq files. Optionally, set an output directory using --output_dir")
}

// Throw error if no input_folder is provided
if (params.input_dir == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --input_dir")
}

/*
Process for running FastQC quality control on given input files. Outputs HTML QC report to folder qc_output.
*/
process fastQC {
    input:
    val fqs
    path out_dir

    output:
    stdout

    """
    mkdir -p $out_dir
    mkdir -p $out_dir/qc_reports
    fastqc -o $out_dir/qc_reports $fqs
    """
}

/*
Process for combining previously generated HTML QC reports into one QC summary using MultiQC. Requires that the fastqc process created a folder named qc_output
*/
process multiQC {
    input:
    val ready
    path out_dir

    """
    multiqc --outdir $out_dir/qc_reports $out_dir/qc_reports
    """
}

process alignment_setup {
    input:
    path genome_fasta
    path genome_gtf
    path input_dir
    path out_dir

    output:
    path genome_index
    stdout

    """
    mkdir -p $PWD/$out_dir
    mkdir -p $PWD/$out_dir/STAR
    mkdir -p $PWD/$out_dir/STAR/genome
    STAR --runMode genomeGenerate \
    --genomeDir $PWD/$out_dir/STAR/genome \
    --genomeFastaFiles $PWD/$input_dir/$genome_fasta
    --sjdbGTFfile $PWD/$input_dir/$genome_gtf
    --sjdbOverhang 99
    """
}

process alignment {
    input:
    path fqs

    """
    
    """
}

workflow onlyQC {
    println("========== Running QC ==========\n")

    out_dir = file("$params.output_dir")
    fqs = channel.from(file("$params.input_dir/*.fastq"))
    
    qc_ready = fastQC(fqs, out_dir) 
    qc_ready.view { print it }
    multiQC(qc_ready.collect(), out_dir)
}

workflow skipQC {
    println("========== Running Alignment & Quantification ==========\n")

    input_dir = file("$params.input_dir")
    out_dir = file("$params.output_dir")
    fqs = channel.from(file("./$params.input_dir/*.fastq"))
    genome_fasta = channel.from(file("./$params.input_dir/*.fna"))
    genome_gtf = channel.from(file("./$params.input_dir/*.gtf"))

    genome_index = alignment_setup(genome_fasta, genome_gtf, input_dir, out_dir)
}

workflow {
    onlyQC()
    skipQC()
}
