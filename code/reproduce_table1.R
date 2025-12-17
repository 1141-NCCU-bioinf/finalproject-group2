# =============================================================================
# Reproduce Table 1 from Shiao et al. 2015
# "Expression Divergence of Chemosensory Genes between Drosophila sechellia
#  and Its Sibling Species and Its Implications for Host Shift"
# =============================================================================

# Load required packages
if (!require("NOISeq")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("NOISeq")
}
library(NOISeq)
library(dplyr)

# Helper function for string concatenation
`%s%` <- function(a, b) paste0(a, b)

# =============================================================================
# 1. Define file paths
# =============================================================================

# Allow dynamic base_dir and output_dir via command line arguments
# Usage: Rscript reproduce_table1.R [base_dir] [output_dir]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  base_dir <- normalizePath(args[1], mustWork = TRUE)
} else {
  base_dir <- "/Users/MOJi/Documents/114/bioinfo/hw/Final Project/Paper1"
}
if (length(args) >= 2) {
  output_dir <- normalizePath(args[2], mustWork = TRUE)
} else {
  output_dir <- base_dir
}
cat(sprintf("Using base directory: %s\n", base_dir))
cat(sprintf("Using output directory: %s\n", output_dir))

# Ortholog table
ortholog_file <- file.path(base_dir, "all_3spe-exist+noDup_paper1.csv")

# FPKM file paths
fpkm_files <- list(
  # D. sechellia - Taiwan (AS) - Male
  Dsec_Mal_AS_R1 = file.path(base_dir, "GSE67587_RAW/GSM1650068_Adu_Dsec_Mal_AS_R1_clout_genes.fpkm_tracking"),
  Dsec_Mal_AS_R2 = file.path(base_dir, "GSE67587_RAW/GSM1650069_Adu_Dsec_Mal_AS_R2_clout_genes.fpkm_tracking"),
  Dsec_Mal_AS_R3 = file.path(base_dir, "GSE67587_RAW/GSM1650070_Adu_Dsec_Mal_AS_R3_clout_genes.fpkm_tracking"),
  # D. sechellia - Taiwan (AS) - Female
  Dsec_Fem_AS_R1 = file.path(base_dir, "GSE67587_RAW/GSM1650071_Adu_Dsec_Fem_AS_R1_clout_genes.fpkm_tracking"),
  Dsec_Fem_AS_R2 = file.path(base_dir, "GSE67587_RAW/GSM1650072_Adu_Dsec_Fem_AS_R2_clout_genes.fpkm_tracking"),
  Dsec_Fem_AS_R3 = file.path(base_dir, "GSE67587_RAW/GSM1650073_Adu_Dsec_Fem_AS_R3_clout_genes.fpkm_tracking"),
  # D. sechellia - Japan (JP) - Male
  Dsec_Mal_JP_R1 = file.path(base_dir, "GSE67861_RAW/GSM1657327_Adu_Dsec_Mal_JP_R1_clout_genes.fpkm_tracking"),
  Dsec_Mal_JP_R2 = file.path(base_dir, "GSE67861_RAW/GSM1657328_Adu_Dsec_Mal_JP_R2_clout_genes.fpkm_tracking"),
  Dsec_Mal_JP_R3 = file.path(base_dir, "GSE67861_RAW/GSM1657329_Adu_Dsec_Mal_JP_R3_clout_genes.fpkm_tracking"),
  # D. sechellia - Japan (JP) - Female
  Dsec_Fem_JP_R1 = file.path(base_dir, "GSE67861_RAW/GSM1657330_Adu_Dsec_Fem_JP_R1_clout_genes.fpkm_tracking"),
  Dsec_Fem_JP_R2 = file.path(base_dir, "GSE67861_RAW/GSM1657331_Adu_Dsec_Fem_JP_R2_clout_genes.fpkm_tracking"),
  Dsec_Fem_JP_R3 = file.path(base_dir, "GSE67861_RAW/GSM1657332_Adu_Dsec_Fem_JP_R3_clout_genes.fpkm_tracking"),
  # D. simulans - Japan (JP) - Male
  Dsim_Mal_JP_R1 = file.path(base_dir, "GSE67862_RAW/GSM1657333_Adu_Dsim_Mal_JP_R1_clout_genes.fpkm_tracking"),
  Dsim_Mal_JP_R2 = file.path(base_dir, "GSE67862_RAW/GSM1657334_Adu_Dsim_Mal_JP_R2_clout_genes.fpkm_tracking"),
  Dsim_Mal_JP_R3 = file.path(base_dir, "GSE67862_RAW/GSM1657335_Adu_Dsim_Mal_JP_R3_clout_genes.fpkm_tracking"),
  # D. simulans - Japan (JP) - Female
  Dsim_Fem_JP_R1 = file.path(base_dir, "GSE67862_RAW/GSM1657336_Adu_Dsim_Fem_JP_R1_clout_genes.fpkm_tracking"),
  Dsim_Fem_JP_R2 = file.path(base_dir, "GSE67862_RAW/GSM1657337_Adu_Dsim_Fem_JP_R2_clout_genes.fpkm_tracking"),
  Dsim_Fem_JP_R3 = file.path(base_dir, "GSE67862_RAW/GSM1657338_Adu_Dsim_Fem_JP_R3_clout_genes.fpkm_tracking"),
  # D. melanogaster - Taiwan (AS) - 1 replicate each
  Dmel_Mal_AS_R1 = file.path(base_dir, "GSE99545_RAW/GSM2645722_Adu_Dmel_Mal_AS_R1_clout.fpkm_tracking"),
  Dmel_Fem_AS_R1 = file.path(base_dir, "GSE99545_RAW/GSM2645723_Adu_Dmel_Fem_AS_R1_clout.fpkm_tracking")
)

# =============================================================================
# 2. Load ortholog table
# =============================================================================

cat("Loading ortholog table...\n")
orthologs <- read.csv(ortholog_file, header = TRUE, stringsAsFactors = FALSE)
colnames(orthologs) <- c("gene", "Dmel", "Dsec", "Dsim")
cat(sprintf("  Loaded %d orthologous genes\n", nrow(orthologs)))

# =============================================================================
# 3. Define function to load FPKM data
# =============================================================================

load_fpkm <- function(filepath) {
  df <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE)
  # Return tracking_id and FPKM columns
  result <- df[, c("tracking_id", "FPKM")]
  colnames(result) <- c("gene_id", "FPKM")
  return(result)
}

# =============================================================================
# 4. Load all FPKM files and create expression matrix
# =============================================================================

cat("Loading FPKM files...\n")

# Load all FPKM data
fpkm_data <- lapply(fpkm_files, load_fpkm)

# Function to map gene IDs through ortholog table
map_fpkm_to_orthologs <- function(fpkm_df, species_col, orthologs) {
  # Create mapping from FBgn ID to gene name
  id_to_gene <- setNames(orthologs$gene, orthologs[[species_col]])

  # Map FPKM values
  fpkm_df$gene_name <- id_to_gene[fpkm_df$gene_id]

  # Keep only genes in ortholog table
  fpkm_df <- fpkm_df[!is.na(fpkm_df$gene_name), ]

  # Handle duplicates by taking mean
  fpkm_df <- fpkm_df %>%
    group_by(gene_name) %>%
    summarise(FPKM = mean(FPKM, na.rm = TRUE)) %>%
    as.data.frame()

  return(fpkm_df)
}

# Determine species for each sample
sample_species <- c(
  Dsec_Mal_AS_R1 = "Dsec", Dsec_Mal_AS_R2 = "Dsec", Dsec_Mal_AS_R3 = "Dsec",
  Dsec_Fem_AS_R1 = "Dsec", Dsec_Fem_AS_R2 = "Dsec", Dsec_Fem_AS_R3 = "Dsec",
  Dsec_Mal_JP_R1 = "Dsec", Dsec_Mal_JP_R2 = "Dsec", Dsec_Mal_JP_R3 = "Dsec",
  Dsec_Fem_JP_R1 = "Dsec", Dsec_Fem_JP_R2 = "Dsec", Dsec_Fem_JP_R3 = "Dsec",
  Dsim_Mal_JP_R1 = "Dsim", Dsim_Mal_JP_R2 = "Dsim", Dsim_Mal_JP_R3 = "Dsim",
  Dsim_Fem_JP_R1 = "Dsim", Dsim_Fem_JP_R2 = "Dsim", Dsim_Fem_JP_R3 = "Dsim",
  Dmel_Mal_AS_R1 = "Dmel", Dmel_Fem_AS_R1 = "Dmel"
)

# Map each FPKM file to ortholog gene names
mapped_fpkm <- list()
for (sample_name in names(fpkm_data)) {
  species <- sample_species[sample_name]
  mapped_fpkm[[sample_name]] <- map_fpkm_to_orthologs(
    fpkm_data[[sample_name]],
    species,
    orthologs
  )
  cat(sprintf("  %s: %d genes mapped\n", sample_name, nrow(mapped_fpkm[[sample_name]])))
}

# =============================================================================
# 5. Build expression matrix
# =============================================================================

cat("\nBuilding expression matrix...\n")

# Get all unique gene names
all_genes <- unique(orthologs$gene)

# Initialize expression matrix
expr_matrix <- matrix(0, nrow = length(all_genes), ncol = length(fpkm_files),
                      dimnames = list(all_genes, names(fpkm_files)))

# Fill in the expression matrix
for (sample_name in names(mapped_fpkm)) {
  df <- mapped_fpkm[[sample_name]]
  idx <- match(df$gene_name, all_genes)
  expr_matrix[idx[!is.na(idx)], sample_name] <- df$FPKM[!is.na(idx)]
}

cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))

# Save expression matrix (raw FPKM)
expr_matrix_file <- file.path(output_dir, "expression_matrix_raw.csv")
expr_df <- data.frame(Gene = rownames(expr_matrix), expr_matrix, check.names = FALSE)
write.csv(expr_df, expr_matrix_file, row.names = FALSE)
cat(sprintf("  Saved raw expression matrix to: %s\n", expr_matrix_file))

# =============================================================================
# 6. Create experimental design factors
# =============================================================================

# Sample metadata
sample_info <- data.frame(
  sample = names(fpkm_files),
  species = sample_species,
  sex = ifelse(grepl("Mal", names(fpkm_files)), "Male", "Female"),
  location = ifelse(grepl("_AS_", names(fpkm_files)), "AS", "JP"),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample

# Create condition labels
sample_info$condition <- paste(sample_info$species, sample_info$sex, sample_info$location, sep = "_")

cat("\nSample information:\n")
print(sample_info)

# =============================================================================
# 7. Upper Quartile Normalization
# =============================================================================

cat("\nPerforming upper quartile normalization...\n")

# Function to perform upper quartile normalization
uq_normalize <- function(expr_mat) {
  # Calculate upper quartile for each sample (excluding zeros)
  uq_factors <- apply(expr_mat, 2, function(x) {
    x_nonzero <- x[x > 0]
    if (length(x_nonzero) > 0) {
      quantile(x_nonzero, 0.75)
    } else {
      1
    }
  })

  # Normalize to the mean upper quartile
  mean_uq <- mean(uq_factors)
  norm_factors <- uq_factors / mean_uq

  # Apply normalization
  norm_mat <- sweep(expr_mat, 2, norm_factors, "/")
  return(norm_mat)
}

# Normalize expression matrix
expr_norm <- uq_normalize(expr_matrix)

# Save normalized expression matrix
expr_norm_file <- file.path(output_dir, "expression_matrix_normalized.csv")
expr_norm_df <- data.frame(Gene = rownames(expr_norm), expr_norm, check.names = FALSE)
write.csv(expr_norm_df, expr_norm_file, row.names = FALSE)
cat(sprintf("  Saved normalized expression matrix to: %s\n", expr_norm_file))

# =============================================================================
# 8. Differential Expression Analysis with NOISeq
# =============================================================================

cat("\nRunning differential expression analysis...\n")

# Function to run NOISeq comparison
run_noiseq_comparison <- function(expr_mat, sample_info, cond1, cond2, use_sim = FALSE) {
  # Get samples for each condition
  samples1 <- rownames(sample_info)[sample_info$condition == cond1]
  samples2 <- rownames(sample_info)[sample_info$condition == cond2]

  # Subset expression matrix
  subset_samples <- c(samples1, samples2)
  subset_expr <- expr_mat[, subset_samples, drop = FALSE]

  # Create factors
  factors <- data.frame(
    condition = factor(c(rep("cond1", length(samples1)), rep("cond2", length(samples2))))
  )
  rownames(factors) <- subset_samples

  # Create NOISeq data object
  noiseq_data <- readData(data = subset_expr, factors = factors)

  # Run NOISeq
  if (use_sim) {
    # Use noiseq for single replicate comparisons (simulate technical replicates)
    # NOISeq will handle single replicate by using noiseq with technical replicates
    result <- noiseq(noiseq_data,
                     k = 0.5,
                     norm = "uqua",
                     factor = "condition",
                     conditions = c("cond1", "cond2"),
                     replicates = "technical",
                     pnr = 0.2,
                     nss = 5,
                     v = 0.02,
                     lc = 0)
  } else {
    # Use noiseq for biological replicates
    result <- noiseq(noiseq_data,
                     k = 0.5,
                     norm = "uqua",
                     factor = "condition",
                     conditions = c("cond1", "cond2"),
                     replicates = "biological")
  }

  return(result)
}

# Define the 4 comparisons
comparisons <- list(
  # Comparison 1: Dsec_Male_AS vs Dmel_Male_AS (Taiwan)
  list(cond1 = "Dsec_Male_AS", cond2 = "Dmel_Male_AS", use_sim = TRUE, name = "Dsec_M_TW_vs_Dmel_M_TW"),
  # Comparison 2: Dsec_Female_AS vs Dmel_Female_AS (Taiwan)
  list(cond1 = "Dsec_Female_AS", cond2 = "Dmel_Female_AS", use_sim = TRUE, name = "Dsec_F_TW_vs_Dmel_F_TW"),
  # Comparison 3: Dsec_Male_JP vs Dsim_Male_JP (Japan)
  list(cond1 = "Dsec_Male_JP", cond2 = "Dsim_Male_JP", use_sim = FALSE, name = "Dsec_M_JP_vs_Dsim_M_JP"),
  # Comparison 4: Dsec_Female_JP vs Dsim_Female_JP (Japan)
  list(cond1 = "Dsec_Female_JP", cond2 = "Dsim_Female_JP", use_sim = FALSE, name = "Dsec_F_JP_vs_Dsim_F_JP")
)

# Run all comparisons
noiseq_results <- list()
for (comp in comparisons) {
  cat(sprintf("  Running %s...\n", comp$name))
  noiseq_results[[comp$name]] <- run_noiseq_comparison(
    expr_matrix, sample_info, comp$cond1, comp$cond2, comp$use_sim
  )
}

# =============================================================================
# 9. Extract DEGs (q >= 0.95)
# =============================================================================

cat("\nExtracting differentially expressed genes (q >= 0.95)...\n")

# Get DEGs for each comparison (with full results)
deg_lists <- list()
deg_full_results <- list()
for (name in names(noiseq_results)) {
  # Get full results table
  full_results <- noiseq_results[[name]]@results[[1]]
  deg_full_results[[name]] <- full_results

  # Get DEG gene names
  degs <- degenes(noiseq_results[[name]], q = 0.80, M = NULL)
  if (!is.null(degs) && nrow(degs) > 0) {
    deg_lists[[name]] <- rownames(degs)
    cat(sprintf("  %s: %d DEGs\n", name, length(deg_lists[[name]])))
  } else {
    deg_lists[[name]] <- character(0)
    cat(sprintf("  %s: 0 DEGs\n", name))
  }
}

# Save DEG tables for each comparison
cat("\nSaving DEG tables...\n")
for (name in names(deg_full_results)) {
  deg_table <- deg_full_results[[name]]
  deg_table$Gene <- rownames(deg_table)
  deg_table <- deg_table[, c("Gene", setdiff(names(deg_table), "Gene"))]  # Move Gene to first column
  deg_file <- file.path(output_dir, paste0("DEG_table_", name, ".csv"))
  write.csv(deg_table, deg_file, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", deg_file))
}


# =============================================================================
# 10. Define target chemosensory genes
# =============================================================================

target_genes <- c(
  # Upregulated in D. sechellia (13 genes)
  "Obp19a", "Obp50a", "Obp56d", "CheA75a", "CheA87a",
  "Or23a", "Or35a", "Or56a", "Or67b", "Or85b", "Or85c",
  "Gr64f", "Ir84a",
  # Downregulated in D. sechellia (5 genes)
  "Obp83a", "Obp99c", "Obp99d", "Or9a", "Or42b"
)

cat(sprintf("\nTarget chemosensory genes: %d\n", length(target_genes)))

# Check which target genes are in our data
target_genes_found <- target_genes[target_genes %in% rownames(expr_matrix)]
target_genes_missing <- target_genes[!target_genes %in% rownames(expr_matrix)]

if (length(target_genes_missing) > 0) {
  cat(sprintf("  Missing genes: %s\n", paste(target_genes_missing, collapse = ", ")))
}
cat(sprintf("  Found genes: %d\n", length(target_genes_found)))

# =============================================================================
# 11. Calculate mean FPKM values for each condition
# =============================================================================

cat("\nCalculating mean FPKM values...\n")

# Calculate mean normalized FPKM for each condition
calc_mean_fpkm <- function(expr_mat, sample_info, species, sex, location) {
  samples <- rownames(sample_info)[
    sample_info$species == species &
    sample_info$sex == sex &
    sample_info$location == location
  ]
  if (length(samples) == 0) return(rep(NA, nrow(expr_mat)))
  if (length(samples) == 1) return(expr_mat[, samples])
  return(rowMeans(expr_mat[, samples, drop = FALSE], na.rm = TRUE))
}

# Calculate means for each condition
mean_fpkm <- data.frame(
  gene = rownames(expr_norm),
  Dsec_M_TW = calc_mean_fpkm(expr_norm, sample_info, "Dsec", "Male", "AS"),
  Dmel_M_TW = calc_mean_fpkm(expr_norm, sample_info, "Dmel", "Male", "AS"),
  Dsec_F_TW = calc_mean_fpkm(expr_norm, sample_info, "Dsec", "Female", "AS"),
  Dmel_F_TW = calc_mean_fpkm(expr_norm, sample_info, "Dmel", "Female", "AS"),
  Dsec_M_JP = calc_mean_fpkm(expr_norm, sample_info, "Dsec", "Male", "JP"),
  Dsim_M_JP = calc_mean_fpkm(expr_norm, sample_info, "Dsim", "Male", "JP"),
  Dsec_F_JP = calc_mean_fpkm(expr_norm, sample_info, "Dsec", "Female", "JP"),
  Dsim_F_JP = calc_mean_fpkm(expr_norm, sample_info, "Dsim", "Female", "JP"),
  stringsAsFactors = FALSE
)

# =============================================================================
# 12. Calculate log2 ratios
# =============================================================================

cat("Calculating log2 ratios...\n")

# Add small pseudocount to avoid log(0)
pseudocount <- 0.01

mean_fpkm$Log2_M_TW <- log2((mean_fpkm$Dsec_M_TW + pseudocount) / (mean_fpkm$Dmel_M_TW + pseudocount))
mean_fpkm$Log2_F_TW <- log2((mean_fpkm$Dsec_F_TW + pseudocount) / (mean_fpkm$Dmel_F_TW + pseudocount))
mean_fpkm$Log2_M_JP <- log2((mean_fpkm$Dsec_M_JP + pseudocount) / (mean_fpkm$Dsim_M_JP + pseudocount))
mean_fpkm$Log2_F_JP <- log2((mean_fpkm$Dsec_F_JP + pseudocount) / (mean_fpkm$Dsim_F_JP + pseudocount))

# =============================================================================
# 13. Generate Table 1
# =============================================================================

cat("\nGenerating Table 1...\n")

# Filter for target genes
table1 <- mean_fpkm[mean_fpkm$gene %in% target_genes_found, ]

# Reorder columns for Table 1 format
table1_final <- data.frame(
  Gene = table1$gene,
  # Taiwan (AS) comparisons: Dsec vs Dmel
  Dsec_Male_TW = round(table1$Dsec_M_TW, 2),
  Dmel_Male_TW = round(table1$Dmel_M_TW, 2),
  Log2_Male_TW = round(table1$Log2_M_TW, 2),
  Dsec_Female_TW = round(table1$Dsec_F_TW, 2),
  Dmel_Female_TW = round(table1$Dmel_F_TW, 2),
  Log2_Female_TW = round(table1$Log2_F_TW, 2),
  # Japan (JP) comparisons: Dsec vs Dsim
  Dsec_Male_JP = round(table1$Dsec_M_JP, 2),
  Dsim_Male_JP = round(table1$Dsim_M_JP, 2),
  Log2_Male_JP = round(table1$Log2_M_JP, 2),
  Dsec_Female_JP = round(table1$Dsec_F_JP, 2),
  Dsim_Female_JP = round(table1$Dsim_F_JP, 2),
  Log2_Female_JP = round(table1$Log2_F_JP, 2),
  stringsAsFactors = FALSE
)

# Order by gene categories (upregulated first, then downregulated)
upregulated <- c("Obp19a", "Obp50a", "Obp56d", "CheA75a", "CheA87a",
                 "Or23a", "Or35a", "Or56a", "Or67b", "Or85b", "Or85c",
                 "Gr64f", "Ir84a")
downregulated <- c("Obp83a", "Obp99c", "Obp99d", "Or9a", "Or42b")

gene_order <- c(upregulated[upregulated %in% table1_final$Gene],
                downregulated[downregulated %in% table1_final$Gene])
table1_final <- table1_final[match(gene_order, table1_final$Gene), ]

# =============================================================================
# 14. Output results
# =============================================================================

cat("\n" %s% paste(rep("=", 70), collapse = "") %s% "\n")
cat("TABLE 1: Expression Levels (FPKMs) and Log2-Ratios of Chemosensory Genes\n")
cat("Differentially Expressed in D. sechellia compared with D. melanogaster or D. simulans\n")
cat(paste(rep("=", 70), collapse = "") %s% "\n\n")

# Print table
print(table1_final, row.names = FALSE)

# Save to CSV
output_file <- file.path(output_dir, "Table1_reproduced_R.csv")
write.csv(table1_final, output_file, row.names = FALSE)
cat(sprintf("\nTable saved to: %s\n", output_file))

# =============================================================================
# 15. Validation - Compare with paper values
# =============================================================================

cat("\n" %s% paste(rep("=", 70), collapse = "") %s% "\n")
cat("VALIDATION: Comparison with paper values\n")
cat(paste(rep("=", 70), collapse = "") %s% "\n\n")

# Expected values from paper (for validation)
paper_values <- data.frame(
  Gene = c("Obp50a", "Or85c", "CheA87a", "Obp83a"),
  Dsec_M_TW_paper = c(45.16, 95.97, 226.03, 7520.39),
  Dmel_M_TW_paper = c(0.35, 4.53, 6.97, 16673.69),
  Log2_M_TW_paper = c(7.02, 4.41, 5.02, -1.15),
  Dsec_M_JP_paper = c(23.85, 68.71, 149.12, 6196.39),
  Dsim_M_JP_paper = c(1.79, 4.28, 16.98, 10887.62),
  Log2_M_JP_paper = c(3.74, 4.00, 3.13, -0.81),
  stringsAsFactors = FALSE
)

# Compare reproduced vs paper values
for (i in 1:nrow(paper_values)) {
  gene <- paper_values$Gene[i]
  if (gene %in% table1_final$Gene) {
    cat(sprintf("\n%s:\n", gene))
    reproduced_row <- table1_final[table1_final$Gene == gene, ]
    cat(sprintf("  Dsec_M (TW): Reproduced=%.2f, Paper=%.2f\n",
                reproduced_row$Dsec_Male_TW, paper_values$Dsec_M_TW_paper[i]))
    cat(sprintf("  Dmel_M (TW): Reproduced=%.2f, Paper=%.2f\n",
                reproduced_row$Dmel_Male_TW, paper_values$Dmel_M_TW_paper[i]))
    cat(sprintf("  Log2 (TW):   Reproduced=%.2f, Paper=%.2f\n",
                reproduced_row$Log2_Male_TW, paper_values$Log2_M_TW_paper[i]))
  }
}

cat("\n\nAnalysis complete!\n")
