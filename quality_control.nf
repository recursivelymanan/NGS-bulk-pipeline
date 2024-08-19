params.input_folder = null
params.help = null

if (params.help == true) {
    print("To use the quality control pipeline, use the --input_folder parameter to indicate the path to the folder containing your .fastq files.")
}

if (params.input_folder == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --input_folder")
}

process fastqc {
    input:
    path fqs

    output:
    stdout

    """
    mkdir -p $PWD/qc_output
    fastqc -o $PWD/qc_output $fqs
    """
}

workflow  {       
    fqs = channel.from(file("./$params.input_folder/*"))
    fastqc(fqs) | view { it }
}