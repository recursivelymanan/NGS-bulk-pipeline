
// Show help message
if (params.help == true) {
    print("To use the quality control pipeline, use the --inputDir parameter to indicate the path to the directory containing your .fastq files. Optionally, set an output directory using --outputDir")
}

// Throw error if no input_folder is provided
if (params.inputDir == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --inputDir")
}

// Throw error if aligner is not hisat2 or star
if (params.aligner != "star" && params.aligner != "hisat2") {
    throw new RuntimeException("Invalid aligner. --aligner must be set to 'star' or 'hisat2'")
}

println("Reading input from directory: ${params.inputDir}")
println("Saving all workflow output to directory: ${params.outputDir}\n")

if (params.downloadReferenceFiles) {
    println("--downloadReferenceFiles option was selected, necessary files will be downloaded.\n")
}
else {
    println("--downloadReferenceFiles option was not selected, using provided genome files in ${params.inputDir}\n")
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
        path: "${params.outputDir}/qc_reports",
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
        path: "${params.outputDir}/ref",
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
        path: "${params.outputDir}/ref",
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
    path indices
    tuple val(sampleID), path(read1), path(read2)

    output:
    path "*.sam", emit: sam
    
    script:
    if (read2 == null) {
        """
        hisat2 -x genome -U $read1 -S ${sampleID}.sam
        """
    }
    else {
        """
        hisat2 -x genome -1 $read1 -2 $read2 -S ${sampleID}.sam
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
    def name = sam.getBaseName()
    """
    samtools view -bSh $sam > ${name}.bam
    """
}

/*
Prepare the genome index directory to prepare for mapping, using STAR with --runMode genomeGenerate.
*/
process alignmentSetupSTAR {
    publishDir(
        path: "${params.outputDir}/STAR/genome",
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
        path: "${params.outputDir}/STAR",
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
        path: "${params.outputDir}/featureCounts",
        mode: "copy"
    )

    input:
    val paired
    path gtf
    path bam

    output:
    path "*.txt"

    script:
    def out_name = bam.getBaseName()

    if (paired) {
        """
        featureCounts -a $gtf -o ${out_name}.txt -p $bam
        """
    }
    else {
        """
        featureCounts -a $gtf -o ${out_name}.txt $bam
        """        
    }
}

/*
Workflow for running QC with FastQC and MultiQC. 
*/
workflow QC {
    fqs = channel.from(file("${params.inputDir}/*.fastq"))
    
    fastQC(fqs)
    multiQC(fastQC.out.html.collect(), fastQC.out.zip.collect())
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
                    tuple(val1.getName(), val1, null)
            }
    }
    
    // Processing workflow using HISAT2 
    if (params.aligner == "hisat2") {
        if (params.downloadReferenceFiles) {
            retrieveGTFhuman()
            genome_gtf = retrieveGTFhuman.out.gtf.collect()
        }
        else {
            genome_gtf = channel.from(file("${params.inputDir}/*.gtf"))
        }
        alignmentSetupHISAT()
        alignmentHISAT(alignmentSetupHISAT.out.indices.collect(), fqs)
        convertToBAM(alignmentHISAT.out.sam)
        quantify(params.paired, genome_gtf, convertToBAM.out.bam)
    }

    // Processing workflow using STAR
    else {
        if (params.downloadReferenceFiles) {
            retrieveFilesHuman()
            renameGenomeFiles(retrieveFilesHuman.out.genome.collect(), retrieveFilesHuman.out.gtf.collect())
            genome_fasta = retrieveFilesHuman.out.genome.collect()
            genome_gtf = retrieveFilesHuman.out.gtf.collect()
        }
        else {
            genome_fasta = channel.from(file("${params.inputDir}/*.fna"))
            genome_gtf = channel.from(file("${params.inputDir}/*.gtf"))
        }
    
        alignmentSetupSTAR(genome_fasta, genome_gtf)
        alignmentSTAR(fqs, alignment_setup.out.index)
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

workflow testWF {
    channel.fromPath("${params.inputDir}/*.fastq")
        .map {
            val1 ->
                tuple(val1.getName(), val1, null)
        }
        .view()
}
