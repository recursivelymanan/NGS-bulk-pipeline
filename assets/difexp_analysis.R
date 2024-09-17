library ("DESeq2")

# Handle parameters from commandline
args <- commandArgs(trailingOnly = TRUE)

custom_comparisons == FALSE
if (length(args) >= 2) {
    cts <- args[1]
    col_data <- args[2]
    comparisons <- unique(col_data$condition)

    if (length(args) == 3) {
        custom_comparisons == TRUE
        comparisons <- strsplit(args[3], split = ",")
    }
} else {
    stop("Incorrect number of parameters passed.")
}

# Load in count and metadata tables
cts <- read.table(cts, sep = "\t", header = TRUE)
col_data <- read.table(col_data, sep = "\t", header = TRUE)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design = ~ condition)
                        
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

dds <- DeSeq(dds)

# Get normalized counts and save file
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = "normalized_counts.csv")

# Get pairwise comparisons and save files
for (x in comparisons) {
    for (y in comparisons) {
        if (x != y) {
            res <- results(dds, contrast = c("condition", x, y))
            write.csv(res,
                      file = paste0("./deseq2_results_",
                                    y,
                                    "-vs-",
                                    x,
                                    ".csv"))
        }
    }
}