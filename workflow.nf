
// Show help message
if (params.help == true) {
    print("To use the quality control pipeline, use the --inputDir parameter to indicate the path to the directory containing your .fastq files. Optionally, set an output directory using --outputDir")
}

// Throw error if no input_folder is provided
if (params.inputDir == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --inputDir")
}

// Throw error if aligner is not hisat2 or star
if (params.aligner != "hisat2" && params.aligner != "star") {
    throw new RuntimeException("Invalid aligner. --aligner must be set to 'star' or 'hisat2'")
}

println("Reading input from directory: ${params.inputDir}")
println("Saving all workflow output to directory: ${params.outputDir}\n")

if (params.downloadReferenceFiles) {
    println("--downloadReferenceFiles option was selected, necessary files will be downloaded to ${params.genomeDir}\n")

}
else {
    println("--downloadReferenceFiles option was not selected, using provided genome files in ${params.genomeDir}\n")
}

println("STARTING WORKFLOW...\n")

/*
Process for running FastQC quality control on given input files. Outputs HTML QC report to folder qc_output.
*/
process fastQC {
    publishDir(
        path: "${params.outputDir}/qc_reports",
        mode: "copy"
    )

    input:
    path fqs
    val args

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    """
    fastqc $args $fqs
    """
}

/*
Process for combining previously generated HTML QC reports into one QC summary using MultiQC. Requires that the fastqc process created a folder named qc_output
*/
process multiQC {
    publishDir (
        path: "${params.outputDir}/qc_reports",
        mode: "copy"
    )

    input:
    path fq_html
    path fq_zip
    val args

    output:
    path "*.html"

    """
    multiqc $args .
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
Prepare the genome index directory to prepare for mapping with HISAT2. Genome index is not generated locally. 
*/
process alignmentSetupHISAT {
    publishDir(
        path: "${params.outputDir}/HISAT2/genome",
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
Map reads using HISAT2. Process can handle single or paired-end reads, as determined by the sampleID value, which is set to null in the case of single reads.
*/
process alignmentHISAT {
    publishDir(
        path: "${params.outputDir}/HISAT2",
        mode: "copy"
    )

    input:
    val paired
    path indices
    tuple val(sampleID), path(read1), path(read2)
    val args

    output:
    path "*.sam", emit: sam
    
    script:
    if (paired) {
        """
        hisat2 -x genome -1 $read1 -2 $read2 -S ${sampleID}.sam $args
        """
    }
    else {
        """
        hisat2 -x genome -U $read1 -S ${sampleID}.sam $args
        """
    }

}

process alignmentSetupSTAR {
    publishDir(
        path: "${params.outputDir}/STAR/genome_files",
        mode: "copy"
    )

    input:
    path genome
    path gtf
    val args

    output:
    path "STAR/*", emit: genomeFiles

    """
    STAR --runMode genomeGenerate \
    --genomeDir ./STAR/ \
    --genomeFastaFiles $genome \
    --sjdbGTFfile $gtf $args
    """
}

process alignmentSTAR {
    publishDir(
        path: "${params.outputDir}/STAR",
        mode: "copy"
    )

    input:
    val paired
    path genomeFiles
    tuple val(sampleID), path(read1), path(read2)
    val args

    output:
    path "*.{sam,bam}", emit: sam

    script:
    if (paired) {
        """
        STAR --genomeDir $genomeFiles \
        --readFilesIn $read1 $read2 \
        --outFileNamePrefix ./ $args
        """
    }
    else {
        """
        STAR --genomeDir $genomeFiles \
        --readFilesIn $read1 \
        --outFileNamePrefix ./ $args
        """        
    }
    
}

/*
Convert SAM files to BAM format.
*/
process convertToBAM {
    publishDir(
	path: "${params.outputDir}/HISAT2",
	mode: "copy"
    )

    input:
    path sam

    output:
    path "*.bam", emit: bam

    script:
    def fileType = sam.getExtension()
    def name = sam.getBaseName()
    if (fileType == "sam") {
        """"
        samtools view -bSh $sam > ${name}.bam
        """
    }
    else {
        """
        echo Conversion not necessary as input file is already .bam
        """
    }
}

/*
Quantify expression with featureCounts
*/
process quantify {
    publishDir(
        path: "${params.outputDir}/featureCounts",
        mode: "copy"
    )

    input:
    val paired
    path gtf
    path bam
    val args

    output:
    path "*.txt"

    script:
    def out_name = bam.getBaseName()

    if (paired) {
        """
        featureCounts -a $gtf -o ${out_name}.txt -p $args $bam 
        """
    }
    else {
        """
        featureCounts -a $gtf -o ${out_name}.txt $args $bam
        """        
    }
}

/*
Workflow for running QC with FastQC and MultiQC. 
*/
workflow QC {
    fqs = channel.fromPath("${params.inputDir}/*.fastq")
    
    fastQC(fqs, params.fqc)
    multiQC(fastQC.out.html.collect(), fastQC.out.zip.collect(), params.mqc)
}

/*
Workflow for post-QC processing with indexing, alignment and quantification. 
*/
workflow processing {
    if (params.paired) {
        fqs = channel.fromFilePairs("${params.inputDir}/*{1,2}.fastq", flat: true)
    }
    else {
        fqs = channel.fromPath("${params.inputDir}/*.fastq")
            .map {
                val1 ->
                    tuple(val1.getName(), val1, file("assets/NO_FILE"))
            }
    }
    
    // Processing workflow using HISAT2 
    if (params.aligner == "hisat2") {
        if (params.downloadReferenceFiles) {
            retrieveGTFhuman()
            genome_gtf = retrieveGTFhuman.out.gtf.collect()
        }
        else {
            genome_gtf = channel.fromPath("${params.genomeDir}/*.gtf")
        }
        alignmentSetupHISAT()
        alignmentHISAT(params.paired, alignmentSetupHISAT.out.indices.collect(), fqs, params.hs2)
        convertToBAM(alignmentHISAT.out.sam)
        quantify(params.paired, genome_gtf, convertToBAM.out.bam, params.fc)
    }
    else if (params.aligner == "star") {
        if (params.downloadReferenceFiles) {
            retrieveFilesHuman()
            genome_gtf = retrieveFilesHuman.out.gtf
            genome_fa = retrieveFilesHuman.out.genome
        }
        else {
            genome_gtf = channel.fromPath("${params.genomeDir}/*.gtf")
            genome_fa = channel.fromPath("${params.genomeDir}/*.fa")
        }

        alignmentSetupSTAR(genome_fa, genome_gtf, params.ss)
        alignmentSTAR(params.paired, alignmentSetupSTAR.out.genomeFiles, fqs, params.star)
        convertToBAM(alignmentSTAR.out.sam)
        quantify(params.paired, genome_gtf, convertToBAM.out.bam, params.fc)
    }
}

/*
Main workflow. 
*/
workflow {
    QC()
    processing()
}
