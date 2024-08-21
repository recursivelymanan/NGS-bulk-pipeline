params.output_dir = "workflow_output"
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

/*
Process for running FastQC quality control on given input files. Outputs HTML QC report to folder qc_output.
*/
process fastQC {
    input:
    path fqs
    path out_dir

    output:
    stdout

    """
    mkdir -p $PWD/$out_dir
    mkdir -p $PWD/$out_dir/qc_reports
    fastqc -o $PWD/$out_dir/qc_reports $fqs
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
    multiqc --outdir $PWD/$out_dir/qc_reports $PWD/$out_dir/qc_reports
    """
}

process alignment {
    input:
    path fqs
    path index

    output:
    path alignments

    """
    
    """
}

process dummy {
    output:
    stdout

    """
    echo successfully_waited
    """
}

workflow onlyQC {
    println("========== QC Only workflow selected, running QC ==========\n")

    out_dir = file("$params.output_dir")
    fqs = channel.from(file("./$params.input_dir/*"))
    
    qc_ready = fastQC(fqs, out_dir) 
    qc_ready.view { print it }
    multiQC(qc_ready.collect(), out_dir)
}

workflow skipQC {
    
}

workflow {
    out_dir = file("$params.output_dir")
    fqs = channel.from(file("./$params.input_dir/*"))
    
    // Quality Control 
    qc_ready = fastQC(fqs, out_dir) 
    qc_ready.view { print it }
    multiQC(qc_ready.collect(), out_dir)

    // Alignment


    // Quantification

}
