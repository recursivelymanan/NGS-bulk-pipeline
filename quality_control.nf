params.output_dir = "qc_output"
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
process fastqc {
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
process multiqc {
    input:
    val ready
    path out_dir

    """
    multiqc --outdir $PWD/$out_dir $PWD/$out_dir/qc_reports
    """
}

workflow  {      
    out_dir = file("$params.output_dir")
    fqs = channel.from(file("./$params.input_dir/*"))
    ready = fastqc(fqs, out_dir) 
    ready.view { print it }
    multiqc(ready.collect(), out_dir)
}