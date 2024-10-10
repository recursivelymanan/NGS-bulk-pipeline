library ("DESeq2")

# Handle parameters from commandline
args <- commandArgs(trailingOnly = TRUE)

custom_comparisons <- FALSE
if (length(args) >= 2) {
    # Load in count and metadata tables
    cts <- read.csv(args[1], sep = "\t", header = TRUE, row.names = 1)
    col_data <- read.csv(args[2], header = TRUE)

    if (length(args) == 3) {
        custom_comparisons <- TRUE
        comparisons <- strsplit(args[3], split = ",")
    }
    else {
        comparisons <- unique(col_data$condition)
    }
} else {
    stop("Incorrect number of parameters passed.")
}

# Load in count and metadata tables
#cts <- read.table(cts, sep = "\t", header = TRUE)
#col_data <- read.table(col_data, sep = "\t", header = TRUE)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = col_data,
                              design = ~ condition)
                        
keep <- rowSums(counts(dds) >= 10) >= 1
dds <- dds[keep, ]

dds <- DESeq(dds)

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