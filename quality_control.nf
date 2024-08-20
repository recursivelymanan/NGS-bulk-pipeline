params.input_folder = null
params.help = null

// Show help message
if (params.help == true) {
    print("To use the quality control pipeline, use the --input_folder parameter to indicate the path to the folder containing your .fastq files.")
}

// Throw error if no input_folder is provided
if (params.input_folder == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --input_folder")
}

/*
Process for running FastQC quality control on given input files. Outputs HTML QC report to folder qc_output.
*/
process fastqc {
    input:
    path fqs

    output:
    stdout

    """
    mkdir -p $PWD/qc_output
    mkdir -p $PWD/qc_output/qc_reports
    fastqc -o $PWD/qc_output/qc_reports $fqs
    """
}

/*
Process for combining previously generated HTML QC reports into one QC summary using MultiQC. Requires that the fastqc process created a folder named qc_output
*/
process multiqc {
    input:
    val ready

    """
    multiqc --outdir $PWD/qc_output $PWD/qc_output/qc_reports
    """
}

workflow  {       
    fqs = channel.from(file("./$params.input_folder/*"))
    ready = fastqc(fqs) 
    ready.view { print it }
    multiqc(ready.collect())
}