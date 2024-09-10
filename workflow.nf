
// Show help message
if (params.help == true) {
    print("To use the quality control pipeline, use the --input_dir parameter to indicate the path to the directory containing your .fastq files. Optionally, set an output directory using --output_dir")
}

// Throw error if no input_folder is provided
if (params.input_dir == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --input_dir")
}

// Throw error if aligner is not hisat2 or star
if (params.aligner != "star" && params.aligner != "hisat2") {
    throw new RuntimeException("Invalid aligner. --aligner must be set to 'star' or 'hisat2'")
}

println("Reading input from directory: $params.input_dir")
println("Saving all workflow output to directory: $params.output_dir\n")

if (params.download_reference_files) {
    println("--download_reference_files option was selected, necessary files will be downloaded.\n")
}
else {
    println("--download_reference_files option was not selected, using provided genome files in $params.input_dir\n")
}

println("STARTING WORKFLOW...\n")



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
Retrieve genomic sequence .fna and gene transfer format .gtf file for the human genome for downstream processes. (Used for STAR workflow)
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
Rename the genome files retrieved by retrieveFilesHuman() in order to remove the long directory chain resulting from the ncbi download. (Used for STAR workflow)
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
Retrieve only the gene transfer format .gtf file for the human genome for downstream processes. (Used for HISAT2 workflow)
*/
process retrieveGTFhuman {
    output:
    path "human_genome/ncbi_dataset/data/GCF_000001405.40/*.gtf", emit: gtf

    """
    datasets download genome accession GCF_000001405.40 --dehydrated --include gtf --filename human_GRCh38_dataset.zip
    unzip human_GRCh38_dataset.zip -d human_genome
    datasets rehydrate --directory human_genome/
    """
}

/*
Rename the GTF file retrieved by retrieveGTFhuman() in order to remove the long directory chain resulting from the ncbi download. (Used for HISAT2 workflow)
*/
process renameGTFfile {
    publishDir (
        path: "$params.output_dir/ref",
        mode: "copy"
    )

    input:
    path gtf

    output:
    path "*.gtf", emit: gtf

    """
    cp $gtf genome.gtf
    """
}

/*
Prepare the genome index directory to prepare for mapping with HISAT2. Genome index is not generated locally. 
*/
process alignmentSetupHISAT {
    publishDir(
        path: "$params.output_dir/HISAT2/genome",
        mode: "copy"
    )

    output:
    path "grch38/*.ht2*", emit: indices

    """
    curl https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz --output hisat2_index.tar.gz
    tar -xvzf hisat2_index.tar.gz
    """
}

/*
Map reads using HISAT2.
*/
process alignmentHISAT {
    publishDir(
        path: "$params.output_dir/HISAT2",
        mode: "copy"
    )

    input:
    path indices
    tuple val(sampleID), path(read1), path(read2)

    output:
    path "*.sam", emit: sam

    """
    hisat2 -x /app/$params.output_dir/HISAT2/genome/grch38/genome -1 $read1 -2 $read2 -S ${sampleID}.sam
    """
}

/*
Convert SAM files to BAM format.
*/
process convertToBAM {
    publishDir(
	path: "$params.output_dir/HISAT2",
	mode: "copy"
    )

    input:
    path sam

    output:
    path "*.bam", emit: bam

    script:
    def name = ${sam}.getBaseName()
    """
    samtools view -bSh $sam > ${name}.bam
    """
}

/*
Prepare the genome index directory to prepare for mapping, using STAR with --runMode genomeGenerate.
*/
process alignmentSetupSTAR {
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
process alignmentSTAR {
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

/*
Quantify expression with featureCounts
*/
process quantify {
    publishDir(
        path: "$params.output_dir/featureCounts",
        mode: "copy"
    )

    input:
    path gtf
    path bams

    output:
    path "*.txt", emit: counts

    script:
    def out_name = 
    """
    featureCounts -a $gtf -o 
    """

}

/*
Workflow for running QC with FastQC and MultiQC. 
*/
workflow QC {
    fqs = channel.from(file("$params.input_dir/*.fastq"))
    
    fastQC(fqs)
    multiQC(fastQC.out.html.collect(), fastQC.out.zip.collect())
}

/*
Workflow for post-QC processing with indexing, alignment and quantification. 
*/
workflow processing {
    paired_fqs = channel.fromFilePairs("$params.input_dir/*{1,2}.fastq", flat: true)
    
    // Processing workflow using HISAT2 
    if (params.aligner == "hisat2") {
        if (params.download_reference_files) {
            retrieveGTFhuman()
            genome_gtf = retrieveGTFhuman.out.gtf.collect()
        }
        else {
            genome_gtf = channel.from(file("$params.input_dir/*.gtf"))
        }
        alignmentSetupHISAT()
        alignmentHISAT(alignmentSetupHISAT.out.indices.collect(), paired_fqs)
	    convertToBAM(alignmentHISAT.out.sam)
        quantify(genome_gtf, convertToBAM.out.bam)
    }

    // Processing workflow using STAR
    else {
        if (params.download_reference_files) {
            retrieveFilesHuman()
            renameGenomeFiles(retrieveFilesHuman.out.genome.collect(), retrieveFilesHuman.out.gtf.collect())
            genome_fasta = retrieveFilesHuman.out.genome.collect()
            genome_gtf = retrieveFilesHuman.out.gtf.collect()
        }
        else {
            genome_fasta = channel.from(file("$params.input_dir/*.fna"))
            genome_gtf = channel.from(file("$params.input_dir/*.gtf"))
        }
    
        alignmentSetupSTAR(genome_fasta, genome_gtf)
        alignmentSTAR(paired_fqs, alignment_setup.out.index)
        quantify(genome_gtf, alignmentSTAR.out.bam)
    }
}

/*
Main workflow. 
*/
workflow {
    QC()
    processing()
}
