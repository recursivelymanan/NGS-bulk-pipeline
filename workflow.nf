def currentDate = new Date().format("MM-dd-yyyy-HH:mm:ss")

println("Saving command used to execute this run in workflow_logs/runcmd_${currentDate}.txt")

"mkdir -p ${params.logDir}".execute().waitFor()
new File("${params.logDir}/runcmd_${currentDate}.txt").withWriter { writer ->
        writer.write(workflow.commandLine)
    }

// Show help message
if (params.help == true) {
    print("To use the quality control pipeline, use the --inputDir parameter to indicate the path to the directory containing your .fastq files. Optionally, set an output directory using --outputDir")
}

// Throw error if no input_folder is provided
if (params.inputDir == null) {  
    throw new RuntimeException("Please provide the path to the folder containing .fastq files by using --inputDir")
}

// Throw error if aligner is not hisat2
if (params.aligner != "hisat2") {
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
    path "*.txt", emit: counts

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
Merge count tables from featureCounts for input into DESeq2
*/
process mergeCounts {
    publishDir(
        path: "${params.outputDir}/featureCounts",
        mode: "copy"
    )

    input:
    path py_script
    path countfiles

    output:
    path "*.txt", emit: counts

    """
    python $py_script $countfiles
    """
}

/*
Differential expression analysis with DESeq2
*/
process diffExpr {
    publishDir(
        path: "${params.outputDir}/DESeq2",
        mode: "copy"
    )

    input:
    path counts
    path exp_design

    output:
    path "*.csv", emit: results

    script:
    if (params.comparisons != null) {
        """
        Rscript ./scripts/difexp_analysis.R $counts $exp_design $params.comparisons
        """
    }
    else {
        """
        Rscript ./scripts/difexp_analysis.R $counts $exp_design
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
    main:
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
            genome_gtf = channel.from(file("${params.inputDir}/*.gtf"))
        }
        alignmentSetupHISAT()
        alignmentHISAT(params.paired, alignmentSetupHISAT.out.indices.collect(), fqs, params.hs2)
        convertToBAM(alignmentHISAT.out.sam)
        quantify(params.paired, genome_gtf, convertToBAM.out.bam, params.fc)
        mergeCounts(channel.from(file("./scripts/merge_counts.py")), quantify.out.counts)
    }

    emit:
    mergeCounts.out.counts
}

workflow diffExp {
    take: 
    counts

    main:

    diffExpr(counts, params.expDesign)
}

/*
Workflow without differential expression analysis. 
*/
workflow noDiffExp {
    QC()
    processing()
}

/*
Main workflow. 
*/
workflow {
    println("By default, DESeq2 analysis is included in the workflow. If you want to run DESeq2, make sure you have a properly set up meta data table. Refer to the README for more instructions.\n")

    QC()
    processing()
    diffExp(processing.out)
}
