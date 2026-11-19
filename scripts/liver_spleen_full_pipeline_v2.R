# 可用 THESIS_DATA 和 LNP_RESULTS 指定输入和输出目录。
data_root <- path.expand(Sys.getenv("THESIS_DATA", "~/Downloads/毕业论文数据"))
result_root <- path.expand(Sys.getenv("LNP_RESULTS", "LNP_results"))
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
# 仅在图形显示层使用 eLNP。
plot_label <- function(x) gsub("(?<![A-Za-z])LNP", "eLNP", as.character(x), perl = TRUE)

#!/usr/bin/env Rscript

# Shared configuration and small utility functions.
# Edit paths and global cutoffs here. Figure-specific settings stay in each script.

######## 00_config_and_helpers ###########
options(stringsAsFactors = FALSE)
set.seed(20260820)

input_dir <- "/nfsdata/disk11/yanwankun_seq/plot/20250118/Liver_Spleen_Bulk"
count_file <- file.path(input_dir, "Liver_Spleen_Bulk_RawCountMatrix3.csv")
tpm_file <- file.path(input_dir, "All_Group_TPM.csv")
kegg_gmt_file <- file.path(input_dir, "kegg_mouse_gageData.gmt")
output_root <- file.path(result_root, "Organ_Bulk")

organ_order <- c("Liver", "Spleen")
treatment_order <- c("Control", "LNP")
contrast_treatments <- "LNP"

deg_fdr_cutoff <- 0.05
deg_pvalue_cutoff <- 0.05
deg_abs_logFC_cutoff <- 1
go_overlap_pvalue_cutoff <- 0.05

group_colors <- c(
  Control = "#6B6B6B",
  LNP = "#FF7F0E"
)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "edgeR", "limma", "ggplot2", "dplyr", "tidyr", "readr", "tibble",
  "scales", "ComplexHeatmap", "circlize", "clusterProfiler",
  "AnnotationDbi", "org.Mm.eg.db"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- as.data.frame(x)
  if (!ncol(x)) writeLines("", path) else readr::write_csv(x, path, na = "")
}

save_gg <- function(plot, stem, width, height, dpi = 350) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
}

read_expression_matrix <- function(path, matrix_type = c("count", "tpm")) {
  matrix_type <- match.arg(matrix_type)
  dat <- readr::read_csv(
    path, show_col_types = FALSE, name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character()), na = character()
  ) %>% as.data.frame()
  gene <- dat[[1]]
  values <- dat[-1]
  numeric_values <- lapply(values, function(x) suppressWarnings(as.numeric(x)))
  numeric_count <- vapply(numeric_values, function(x) sum(!is.na(x)), integer(1))
  schema_rows <- which(rowSums(!is.na(as.data.frame(numeric_values))) == 0)
  if (length(schema_rows)) {
    warning("Removing ", length(schema_rows), " schema/type row(s) from ", basename(path))
    gene <- gene[-schema_rows]
    numeric_values <- lapply(numeric_values, function(x) x[-schema_rows])
  }
  mat <- as.matrix(as.data.frame(numeric_values, check.names = FALSE))
  rownames(mat) <- make.unique(as.character(gene))
  storage.mode(mat) <- "numeric"
  if (matrix_type == "count") mat <- round(mat)
  mat
}

parse_samples <- function(samples) {
  bind_rows(lapply(samples, function(x) {
    part <- strsplit(x, "-", fixed = TRUE)[[1]]
    if (length(part) < 3) stop("Cannot parse sample name: ", x)
    treatment <- case_when(
      tolower(part[1]) == "control" ~ "Control",
      tolower(part[1]) == "lnp" ~ "LNP",
      TRUE ~ NA_character_
    )
    organ <- case_when(
      tolower(part[2]) == "liver" ~ "Liver",
      tolower(part[2]) == "spleen" ~ "Spleen",
      TRUE ~ NA_character_
    )
    replicate <- suppressWarnings(as.integer(part[3]))
    data.frame(Sample = x, Treatment = treatment, Organ = organ,
               Replicate = replicate, Batch = factor(replicate))
  }))
}

symbol_to_entrez <- function(symbols) {
  symbols <- unique(na.omit(as.character(symbols)))
  if (!length(symbols)) return(data.frame(SYMBOL = character(), ENTREZID = character()))
  suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db, keys = symbols, keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )) %>%
    as.data.frame() %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(SYMBOL, .keep_all = TRUE)
}

read_offline_kegg_gmt <- function(path = kegg_gmt_file) {
  if (!file.exists(path)) stop("Offline KEGG GMT not found: ", path)
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 3]
  labels <- vapply(fields, `[[`, character(1), 1)
  ids <- sub(" .*", "", labels)
  names <- sub("^[^ ]+ +", "", labels)
  term2gene <- bind_rows(lapply(seq_along(fields), function(i) {
    data.frame(term = ids[i], gene = unique(fields[[i]][-(1:2)]))
  }))
  term2name <- data.frame(term = ids, name = names) %>% distinct(term, .keep_all = TRUE)
  list(term2gene = term2gene, term2name = term2name)
}

load_prepared <- function() {
  prep_dir <- file.path(output_root, "00_prepared_data")
  metadata <- readr::read_csv(file.path(prep_dir, "sample_metadata.csv"), show_col_types = FALSE)
  metadata$Treatment <- factor(metadata$Treatment, levels = treatment_order)
  metadata$Organ <- factor(metadata$Organ, levels = organ_order)
  metadata$Batch <- factor(metadata$Batch)
  tpm <- readRDS(file.path(prep_dir, "TPM_matrix.rds"))
  log_tpm <- readRDS(file.path(prep_dir, "log2TPM_matrix.rds"))
  log_tpm_batch_removed <- readRDS(file.path(prep_dir, "log2TPM_batch_removed.rds"))
  list(metadata = metadata, tpm = tpm, log_tpm = log_tpm,
       log_tpm_batch_removed = log_tpm_batch_removed)
}

load_deg_table <- function(organ, treatment) {
  readr::read_csv(
    file.path(output_root, "01_differential_expression",
              paste0(organ, "_", treatment, "_vs_Control_all_genes.csv.gz")),
    show_col_types = FALSE
  )
}

load_deg_union <- function(organ, basis = c("rawP", "FDR")) {
  basis <- match.arg(basis)
  suffix <- if (basis == "rawP") "rawP_exploratory" else "FDR"
  path <- file.path(output_root, "01_differential_expression",
                    paste0(organ, "_LNP_DEGs_", suffix, ".csv"))
  readr::read_csv(path, show_col_types = FALSE)$gene
}

######## 01_prepare_tpm_and_metadata ######## 

#!/usr/bin/env Rscript

# Run once before the figure scripts.
# Outputs raw TPM, log2(TPM+1), batch-removed log2TPM and sample metadata.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

out_dir <- file.path(output_root, "00_prepared_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 使用归档的原始计数与 TPM，按样本名关联实验记录的批次。
archive <- path.expand(Sys.getenv("ORGAN_PREPARED", "~/Downloads/Liver_Spleen_Modular_Analysis/00_prepared_data"))
counts <- readRDS(file.path(archive, "raw_count_matrix.rds"))
tpm <- readRDS(file.path(archive, "TPM_matrix.rds"))
metadata <- read.csv(file.path(archive, "sample_metadata.csv"), stringsAsFactors = FALSE)
metadata <- subset(metadata, Treatment %in% treatment_order & Organ %in% organ_order)
samples <- metadata$Sample
stopifnot(!anyDuplicated(samples), all(samples %in% colnames(counts)), all(samples %in% colnames(tpm)),
          all(table(factor(metadata$Organ, organ_order), factor(metadata$Treatment, treatment_order)) >= 2))
counts <- counts[, samples, drop = FALSE]
tpm <- tpm[, samples, drop = FALSE]
metadata$Treatment <- factor(metadata$Treatment, treatment_order)
metadata$Organ <- factor(metadata$Organ, organ_order)
metadata$Batch <- factor(metadata$Batch)

matrix_qc <- data.frame(
  Metric = c("Count_genes", "TPM_genes", "Shared_genes",
             "Count_NA_cells", "TPM_NA_cells"),
  Value = c(nrow(counts), nrow(tpm),
            length(intersect(rownames(counts), rownames(tpm))),
            sum(!is.finite(counts)), sum(!is.finite(tpm)))
)
write_csv_safe(matrix_qc, file.path(out_dir, "count_TPM_gene_QC.csv"))

keep <- rowSums(tpm) > 0
tpm <- tpm[keep, , drop = FALSE]
log_tpm <- log2(tpm + 1)
condition <- interaction(metadata$Organ, metadata$Treatment, drop = TRUE)
display_design <- model.matrix(~ 0 + condition)
log_tpm_batch_removed <- limma::removeBatchEffect(
  log_tpm, batch = metadata$Batch, design = display_design
)

saveRDS(counts, file.path(out_dir, "raw_count_matrix.rds"))
saveRDS(tpm, file.path(out_dir, "TPM_matrix.rds"))
saveRDS(log_tpm, file.path(out_dir, "log2TPM_matrix.rds"))
saveRDS(log_tpm_batch_removed, file.path(out_dir, "log2TPM_batch_removed.rds"))
write_csv_safe(metadata, file.path(out_dir, "sample_metadata.csv"))
write_csv_safe(metadata %>% count(Organ, Treatment, Batch, name = "n"),
               file.path(out_dir, "sample_group_batch_counts.csv"))

message("Prepared data saved to: ", out_dir)

######## 02_differential_expression ######## 

#!/usr/bin/env Rscript

# Differential expression only.
# RAW COUNTS: edgeR filterByExpr + TMM + limma-voom; Batch is included in model.
# Saves complete statistics, strict FDR lists, exploratory raw-P lists and unions.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

prep <- load_prepared()
metadata <- prep$metadata
counts <- readRDS(file.path(output_root, "00_prepared_data", "raw_count_matrix.rds"))
out_dir <- file.path(output_root, "01_differential_expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# A count matrix may contain isolated blank/non-numeric cells even after a
# schema/type row was removed. Such rows cannot be used by edgeR. Do not turn
# unknown values into zero; record and remove the affected genes instead.
bad_count_rows <- which(rowSums(!is.finite(counts)) > 0)
bad_count_report <- if (length(bad_count_rows)) {
  data.frame(
    gene = rownames(counts)[bad_count_rows],
    MissingCellCount = rowSums(!is.finite(counts))[bad_count_rows],
    stringsAsFactors = FALSE
  )
} else {
  data.frame(gene = character(), MissingCellCount = integer())
}
write_csv_safe(bad_count_report, file.path(out_dir, "count_genes_removed_for_missing_values.csv"))
if (length(bad_count_rows)) {
  warning("Removing ", length(bad_count_rows),
          " count gene row(s) containing NA/non-finite values; see count_genes_removed_for_missing_values.csv")
  counts <- counts[-bad_count_rows, , drop = FALSE]
}
if (any(counts < 0)) stop("Negative values detected in raw counts.")
if (anyNA(counts) || any(!is.finite(counts))) stop("Count matrix still contains invalid values.")

all_results <- list()
summary_rows <- list()

for (organ_name in organ_order) {
  meta <- metadata %>% filter(Organ == organ_name) %>%
    arrange(match(Treatment, treatment_order), Replicate)
  y <- edgeR::DGEList(counts = counts[, meta$Sample, drop = FALSE])
  design <- model.matrix(~ Treatment + Batch, data = meta)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit, robust = TRUE)
  
  write_csv_safe(data.frame(gene = rownames(counts), PassedFilter = rownames(counts) %in% rownames(y)),
                 file.path(out_dir, paste0(organ_name, "_filterByExpr_status.csv")))
  write_csv_safe(rownames_to_column(as.data.frame(design), "Sample"),
                 file.path(out_dir, paste0(organ_name, "_design_matrix.csv")))
  saveRDS(v$E, file.path(out_dir, paste0(organ_name, "_voom_logCPM.rds")))
  
  organ_tables <- list()
  for (trt in contrast_treatments) {
    coef_name <- paste0("Treatment", trt)
    tab <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
      rownames_to_column("gene") %>%
      mutate(
        Organ = organ_name, Treatment = trt,
        significant = case_when(
          adj.P.Val < deg_fdr_cutoff & logFC >= deg_abs_logFC_cutoff ~ "Up",
          adj.P.Val < deg_fdr_cutoff & logFC <= -deg_abs_logFC_cutoff ~ "Down",
          TRUE ~ "Not"
        ),
        exploratory_significant = case_when(
          P.Value < deg_pvalue_cutoff & logFC >= deg_abs_logFC_cutoff ~ "Up",
          P.Value < deg_pvalue_cutoff & logFC <= -deg_abs_logFC_cutoff ~ "Down",
          TRUE ~ "Not"
        )
      )
    organ_tables[[trt]] <- tab
    all_results[[paste(organ_name, trt, sep = "_")]] <- tab
    
    write_csv_safe(tab, file.path(out_dir, paste0(organ_name, "_", trt, "_vs_Control_all_genes.csv.gz")))
    write_csv_safe(filter(tab, significant != "Not"),
                   file.path(out_dir, paste0(organ_name, "_", trt, "_vs_Control_DEG_FDR.csv")))
    write_csv_safe(filter(tab, exploratory_significant != "Not"),
                   file.path(out_dir, paste0(organ_name, "_", trt, "_vs_Control_DEG_rawP_exploratory.csv")))
    
    summary_rows[[paste(organ_name, trt)]] <- data.frame(
      Organ = organ_name, Treatment = trt, GenesTested = nrow(tab),
      Up_FDR = sum(tab$significant == "Up"), Down_FDR = sum(tab$significant == "Down"),
      Up_rawP = sum(tab$exploratory_significant == "Up"),
      Down_rawP = sum(tab$exploratory_significant == "Down")
    )
  }
  
  union_fdr <- unique(unlist(lapply(organ_tables, function(x) x$gene[x$significant != "Not"])))
  union_rawp <- unique(unlist(lapply(organ_tables, function(x) x$gene[x$exploratory_significant != "Not"])))
  write_csv_safe(data.frame(gene = union_fdr),
                 file.path(out_dir, paste0(organ_name, "_LNP_DEGs_FDR.csv")))
  write_csv_safe(data.frame(gene = union_rawp),
                 file.path(out_dir, paste0(organ_name, "_LNP_DEGs_rawP_exploratory.csv")))
}

write_csv_safe(bind_rows(summary_rows), file.path(out_dir, "DEG_count_summary.csv"))
saveRDS(all_results, file.path(out_dir, "all_DEG_tables.rds"))
message("Differential-expression results saved to: ", out_dir)


######## 03_figure_PCA ######## 

#!/usr/bin/env Rscript

# Figure 1: PCA based on batch-removed log2(TPM + 1), all expressed genes.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

dat <- load_prepared()
metadata <- dat$metadata
mat <- dat$log_tpm_batch_removed
out_dir <- file.path(output_root, "02_PCA")

pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
plot_df <- metadata %>% mutate(
  PC1 = pca$x[Sample, 1], PC2 = pca$x[Sample, 2],
  Shape = Organ, Label = Sample
)
write_csv_safe(plot_df, file.path(out_dir, "PCA_sample_coordinates.csv"))
write_csv_safe(data.frame(PC = paste0("PC", seq_along(variance)), VariancePercent = variance),
               file.path(out_dir, "PCA_variance_explained.csv"))

p <- ggplot(plot_df, aes(PC1, PC2, color = Treatment, shape = Organ, label = plot_label(Label))) +
  geom_point(size = 4) +
  geom_text(nudge_y = 1.5, show.legend = FALSE, size = 3.5) +
  scale_color_manual(values = group_colors, labels=plot_label) +
  scale_fill_manual(values = group_colors, labels=plot_label) +
  scale_shape_manual(values = c(Liver = 16, Spleen = 17)) +
  labs(x = sprintf("PC1 (%.1f%%)", variance[1]),
       y = sprintf("PC2 (%.1f%%)", variance[2])) +
  theme_classic(base_size = 13)

save_gg(p, file.path(out_dir, "Figure_01_PCA"), 9, 6.5)

######## 04_figure_sample_similarity ######## 
#!/usr/bin/env Rscript

# Figure 2: sample similarity based on TPM.
# For each organ, use the UNION of three treatment-vs-Control exploratory DEGs.
# TPM -> log2(TPM + 1) -> remove batch while preserving Organ x Treatment
# -> select DEG union -> Pearson correlation.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

dat <- load_prepared()
metadata <- dat$metadata

# Prepared in 01_prepare_tpm_and_metadata.txt. This matrix has already undergone
# exactly one log2(TPM + 1) transformation and visualization-only batch removal.
mat_all <- dat$log_tpm_batch_removed
out_dir <- file.path(output_root, "03_sample_similarity")

# Organ-specific color scales. These affect display only, not correlations.
color_parameters <- data.frame(
  Organ = c("Liver", "Spleen"),
  Midpoint = c(0.90, 0.94),
  LowerReference = c(0.80, 0.85)
)
write_csv_safe(
  color_parameters,
  file.path(out_dir, "sample_similarity_color_parameters.csv")
)

for (organ_name in organ_order) {
  meta <- metadata %>%
    filter(Organ == organ_name) %>%
    arrange(match(Treatment, treatment_order), Replicate)
  
  genes <- intersect(
    load_deg_union(organ_name, "rawP"),
    rownames(mat_all)
  )
  if (length(genes) < 2) {
    stop("Fewer than two TPM-matched DEG-union genes: ", organ_name)
  }
  
  # mat_all is already batch-removed log2(TPM + 1). Do not log again.
  mat <- mat_all[genes, meta$Sample, drop = FALSE]
  if (any(!is.finite(mat))) {
    stop(organ_name, " batch-removed log2TPM matrix contains NA/Inf values")
  }
  
  cor_mat <- cor(mat, method = "pearson", use = "everything")
  
  write_csv_safe(
    data.frame(gene = genes),
    file.path(out_dir, paste0(organ_name, "_TPM_genes_used.csv"))
  )
  write_csv_safe(
    rownames_to_column(as.data.frame(mat), "gene"),
    file.path(out_dir, paste0(organ_name, "_batch_removed_log2TPM_DEG_matrix.csv.gz"))
  )
  write_csv_safe(
    rownames_to_column(as.data.frame(cor_mat), "Sample"),
    file.path(out_dir, paste0(organ_name, "_TPM_sample_correlation.csv"))
  )
  
  long <- as.data.frame(as.table(cor_mat))
  colnames(long) <- c("SampleY", "SampleX", "Correlation")
  long$SampleX <- factor(long$SampleX, levels = meta$Sample)
  long$SampleY <- factor(long$SampleY, levels = rev(meta$Sample))
  
  organ_color <- color_parameters[color_parameters$Organ == organ_name, ]
  color_midpoint <- organ_color$Midpoint
  color_lower_limit <- min(organ_color$LowerReference, min(cor_mat))
  
  p <- ggplot(long, aes(SampleX, SampleY, fill = Correlation)) + scale_x_discrete(labels=plot_label) + scale_y_discrete(labels=plot_label) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2B55A1", mid = "white", high = "#E53B2F",
      midpoint = color_midpoint,
      limits = c(color_lower_limit, 1),
      oob = scales::squish
    ) +
    coord_fixed() +
    labs(title = organ_name, x = NULL, y = NULL, fill = "Correlation") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      axis.text.x = element_text(angle = 40, hjust = 1),
      panel.grid = element_blank()
    )
  
  save_gg(
    p,
    file.path(out_dir, paste0("Figure_02_", organ_name, "_TPM_similarity")),
    7, 6
  )
}

######## 05_figure_DEG_counts ###########

#!/usr/bin/env Rscript

# Figure 3: up/down DEG counts. Default plot uses exploratory raw P < 0.05.
# Strict FDR counts remain in DEG_count_summary.csv.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

out_dir <- file.path(output_root, "04_DEG_counts")
summary <- readr::read_csv(file.path(output_root, "01_differential_expression", "DEG_count_summary.csv"),
                           show_col_types = FALSE)
plot_df <- summary %>%
  select(Organ, Treatment, Up = Up_rawP, Down = Down_rawP) %>%
  pivot_longer(c(Up, Down), names_to = "Direction", values_to = "Count") %>%
  mutate(PlotCount = ifelse(Direction == "Down", -Count, Count),
         Treatment = factor(Treatment, levels = contrast_treatments),
         Fill = paste(Treatment, Direction, sep = "."))
write_csv_safe(plot_df, file.path(out_dir, "DEG_count_plot_data.csv"))

p <- ggplot(plot_df, aes(Treatment, PlotCount, fill = Fill)) + scale_x_discrete(labels=plot_label) +
  geom_col(width = 0.68) +
  geom_text(aes(label = Count),
            vjust = ifelse(plot_df$PlotCount >= 0, -0.35, 1.35), size = 3.8) +
  facet_wrap(~Organ, nrow = 1, scales = "free_y") +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.14))) +
  scale_fill_manual(values = c(
    "LNP.Up" = "#FF7F0E", "LNP.Down" = "#F8C68F"
  ), labels=plot_label) +
  labs(x = NULL, y = "Differentially expressed genes", fill = NULL) +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold", size = 16))

save_gg(p, file.path(out_dir, "Figure_03_DEG_counts"), 8, 5.6)

######## 06_figure_DEG_heatmap ###########

#!/usr/bin/env Rscript

# Figure 4: DEG expression heatmap.
# Genes: three-treatment exploratory DEG union; values: row Z-score of
# batch-removed log2(TPM + 1).

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize)})

dat <- load_prepared()
metadata <- dat$metadata
mat_all <- dat$log_tpm_batch_removed
out_dir <- file.path(output_root, "05_DEG_heatmap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

row_zscore <- function(x) {
  z <- t(apply(x, 1, function(v) if (sd(v) == 0) rep(0, length(v)) else (v - mean(v)) / sd(v)))
  rownames(z) <- rownames(x); colnames(z) <- colnames(x)
  pmax(pmin(z, 2), -2)
}

for (organ_name in organ_order) {
  meta <- metadata %>% filter(Organ == organ_name) %>%
    arrange(match(Treatment, treatment_order), Replicate)
  genes <- intersect(load_deg_union(organ_name, "rawP"), rownames(mat_all))
  z <- row_zscore(mat_all[genes, meta$Sample, drop = FALSE])
  write_csv_safe(rownames_to_column(as.data.frame(z), "gene"),
                 file.path(out_dir, paste0(organ_name, "_heatmap_Zscore_matrix.csv")))
  ha <- HeatmapAnnotation(
    Treatment = plot_label(meta$Treatment),
    col = list(Treatment = setNames(group_colors,plot_label(names(group_colors)))), show_annotation_name = FALSE
  )
  ht <- Heatmap(
    z, name = "Z-score", top_annotation = ha,
    col = circlize::colorRamp2(c(-2, 0, 2), c("#3B5BA9", "white", "#E84A3A")),
    show_row_names = FALSE, show_column_names = TRUE, column_labels=plot_label(colnames(z)),
    cluster_columns = FALSE, column_order = meta$Sample,
    column_names_rot = 40, border = TRUE
  )
  pdf(file.path(out_dir, paste0("Figure_04_", organ_name, "_DEG_heatmap.pdf")), 7, 10)
  draw(ht); dev.off()
  png(file.path(out_dir, paste0("Figure_04_", organ_name, "_DEG_heatmap.png")),
      width = 2450, height = 3500, res = 350)
  draw(ht); dev.off()
}

######## 07_figure_GO_BP_up ###########

#!/usr/bin/env Rscript

# Figure 5: GO Biological Process for UP-regulated exploratory DEGs.
# Overlap threshold is raw enrichment pvalue < 0.05 (NOT adjusted P value).
# Edit priority_terms to control the ten displayed pathways.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

direction <- "Up"
out_dir <- file.path(output_root, "06_GO_BP_up")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

priority_terms <- c(
  "regulation of response to cytokine stimulus",
  "regulation of apoptotic signaling pathway",
  "positive regulation of leukocyte mediated immunity",
  "positive regulation of inflammatory response",
  "inflammatory response to antigenic stimulus",
  "chronic inflammatory response",
  "cellular response to oxidative stress",
  "cellular response to biotic stimulus", "autophagy",
  "acute inflammatory response"
)

go_results <- list()
for (organ_name in organ_order) {
  for (trt in contrast_treatments) {
    tab <- load_deg_table(organ_name, trt)
    genes <- tab$gene[tab$exploratory_significant == direction]
    mapped <- symbol_to_entrez(genes)
    universe <- symbol_to_entrez(tab$gene)$ENTREZID
    eg <- if (nrow(mapped) >= 3) suppressMessages(clusterProfiler::enrichGO(
      gene = unique(mapped$ENTREZID), universe = unique(universe),
      OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",
      pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
    )) else NULL
    result <- if (is.null(eg)) data.frame() else as.data.frame(eg)
    if (nrow(result)) {
      result <- result %>% mutate(Organ = organ_name, Treatment = trt, Direction = direction)
    }
    go_results[[paste(organ_name, trt, sep = "_")]] <- result
    write_csv_safe(result, file.path(out_dir, paste0(organ_name, "_", trt, "_GO_BP_all.csv")))
  }
}
go_all <- bind_rows(go_results)
write_csv_safe(go_all, file.path(out_dir, "GO_BP_up_all_comparisons.csv"))

if (nrow(go_all)) {
available <- go_all %>% group_by(ID, Description) %>%
  summarise(best_p = min(pvalue, na.rm = TRUE), .groups = "drop") %>% arrange(best_p)
selected_names <- intersect(priority_terms, available$Description)
if (length(selected_names) < 10) {
  selected_names <- c(selected_names,
                      head(setdiff(available$Description, selected_names), 10 - length(selected_names)))
}
selected_names <- head(selected_names, 10)
plot_data <- expand_grid(Description = selected_names, Organ = organ_order,
                         Treatment = contrast_treatments) %>%
  left_join(go_all %>% select(Description, Organ, Treatment, pvalue),
            by = c("Description", "Organ", "Treatment")) %>%
  group_by(Description, Organ, Treatment) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(score = ifelse(is.na(pvalue), 0, pmin(-log10(pmax(pvalue, 1e-300)), 30)),
         Description = factor(Description, levels = rev(selected_names)),
         Treatment = factor(Treatment, levels = contrast_treatments))
write_csv_safe(plot_data, file.path(out_dir, "GO_BP_up_selected_plot_data.csv"))

p_heat <- ggplot(plot_data, aes(Treatment, Description, fill = score)) + scale_x_discrete(labels=plot_label) +
  geom_tile(color = "white", linewidth = 0.8) + facet_wrap(~Organ, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#A7191F") +
  labs(x = NULL, y = NULL, fill = "-log10(P)") + theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 15))
save_gg(p_heat, file.path(out_dir, "Figure_05_GO_BP_up_heatmap"), 9, 6.5)
}

######## 08_figure_GO_BP_down ###########

#!/usr/bin/env Rscript

# Figure 6: GO Biological Process for DOWN-regulated exploratory DEGs.
# Overlap threshold is raw enrichment pvalue < 0.05 (NOT adjusted P value).
# Edit priority_terms to control the ten displayed pathways.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

direction <- "Down"
out_dir <- file.path(output_root, "07_GO_BP_down")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

priority_terms <- c(
  "regulation of lipid storage", "regulation of anoikis",
  "positive regulation of macroautophagy",
  "positive regulation of fatty acid beta-oxidation",
  "positive regulation of autophagy of mitochondrion",
  "oxidoreductase activity, acting on the CH-OH group of donors",
  "NK T cell differentiation", "NF-kappaB binding",
  "negative regulation of NF-kappaB transcription factor activity",
  "intrinsic apoptotic signaling pathway"
)

go_results <- list()
for (organ_name in organ_order) {
  for (trt in contrast_treatments) {
    tab <- load_deg_table(organ_name, trt)
    genes <- tab$gene[tab$exploratory_significant == direction]
    mapped <- symbol_to_entrez(genes)
    universe <- symbol_to_entrez(tab$gene)$ENTREZID
    eg <- if (nrow(mapped) >= 3) suppressMessages(clusterProfiler::enrichGO(
      gene = unique(mapped$ENTREZID), universe = unique(universe),
      OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",
      pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE
    )) else NULL
    result <- if (is.null(eg)) data.frame() else as.data.frame(eg)
    if (nrow(result)) {
      result <- result %>% mutate(Organ = organ_name, Treatment = trt, Direction = direction)
    }
    go_results[[paste(organ_name, trt, sep = "_")]] <- result
    write_csv_safe(result, file.path(out_dir, paste0(organ_name, "_", trt, "_GO_BP_all.csv")))
  }
}
go_all <- bind_rows(go_results)
write_csv_safe(go_all, file.path(out_dir, "GO_BP_down_all_comparisons.csv"))

if (nrow(go_all)) {
available <- go_all %>% group_by(ID, Description) %>%
  summarise(best_p = min(pvalue, na.rm = TRUE), .groups = "drop") %>% arrange(best_p)
selected_names <- intersect(priority_terms, available$Description)
if (length(selected_names) < 10) selected_names <- c(
  selected_names, head(setdiff(available$Description, selected_names), 10 - length(selected_names))
)
selected_names <- head(selected_names, 10)
plot_data <- expand_grid(Description = selected_names, Organ = organ_order,
                         Treatment = contrast_treatments) %>%
  left_join(go_all %>% select(Description, Organ, Treatment, pvalue),
            by = c("Description", "Organ", "Treatment")) %>%
  group_by(Description, Organ, Treatment) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(score = ifelse(is.na(pvalue), 0, pmin(-log10(pmax(pvalue, 1e-300)), 30)),
         Description = factor(Description, levels = rev(selected_names)),
         Treatment = factor(Treatment, levels = contrast_treatments))
write_csv_safe(plot_data, file.path(out_dir, "GO_BP_down_selected_plot_data.csv"))

p_heat <- ggplot(plot_data, aes(Treatment, Description, fill = score)) + scale_x_discrete(labels=plot_label) +
  geom_tile(color = "white", linewidth = 0.8) + facet_wrap(~Organ, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#244A7C") +
  labs(x = NULL, y = NULL, fill = "-log10(P)") + theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 15))
save_gg(p_heat, file.path(out_dir, "Figure_06_GO_BP_down_heatmap"), 9, 6.5)
}

######## 10_figure_GSEA ###########

#!/usr/bin/env Rscript

# Figure 8: offline KEGG GSEA, ranked by log2FC from all tested genes.
# Only pathways with NES > 0 are eligible for the selected display.
# Edit priority_terms to change displayed pathways.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

out_dir <- file.path(output_root, "09_GSEA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
kegg_gmt_file <- Sys.getenv("MOUSE_KEGG_GMT", kegg_gmt_file)
if (file.exists(kegg_gmt_file)) {
kegg <- read_offline_kegg_gmt()

priority_terms <- c(
  "Hepatocellular carcinoma", "Efferocytosis", "IL-17 signaling pathway",
  "Cytokine-cytokine receptor interaction", "Natural killer cell mediated cytotoxicity",
  "Apoptosis", "Inflammatory bowel disease", "Phagosome",
  "TNF signaling pathway", "Chemokine signaling pathway", "Lysosome"
)

results <- list()
for (organ_name in organ_order) for (trt in contrast_treatments) {
  tab <- load_deg_table(organ_name, trt)
  mapped <- symbol_to_entrez(tab$gene)
  ranked <- tab %>% inner_join(mapped, by = c("gene" = "SYMBOL")) %>%
    filter(is.finite(logFC)) %>% group_by(ENTREZID) %>%
    slice_max(abs(logFC), n = 1, with_ties = FALSE) %>% ungroup() %>% arrange(desc(logFC))
  gene_list <- ranked$logFC; names(gene_list) <- ranked$ENTREZID
  set.seed(20260820)
  gsea <- suppressMessages(clusterProfiler::GSEA(
    geneList = gene_list, TERM2GENE = kegg$term2gene, TERM2NAME = kegg$term2name,
    minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
    pAdjustMethod = "BH", verbose = FALSE, seed = TRUE
  ))
  out <- as.data.frame(gsea)
  if (nrow(out)) out <- out %>% mutate(
    Organ = organ_name, Treatment = trt,
    CoreGeneCount = lengths(strsplit(core_enrichment, "/", fixed = TRUE)),
    CoreGeneRatio = CoreGeneCount / setSize
  )
  results[[paste(organ_name, trt, sep = "_")]] <- out
  write_csv_safe(out, file.path(out_dir, paste0(organ_name, "_", trt, "_GSEA_all.csv")))
}

all_gsea <- bind_rows(results)
positive <- all_gsea %>% filter(is.finite(NES), NES > 0)
write_csv_safe(all_gsea, file.path(out_dir, "GSEA_all_comparisons.csv"))
write_csv_safe(positive, file.path(out_dir, "GSEA_NES_positive_only.csv"))

available <- positive %>% group_by(Description) %>%
  summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>% arrange(best_padj)
selected <- intersect(priority_terms, available$Description)
if (length(selected) < 10) selected <- c(
  selected, head(setdiff(available$Description, selected), 10 - length(selected))
)
selected <- head(selected, 10)
plot_data <- positive %>% filter(Description %in% selected) %>% mutate(
  Description = factor(Description, levels = rev(selected)),
  Treatment = factor(Treatment, levels = contrast_treatments),
  Organ = factor(Organ, levels = organ_order),
  logP = pmin(-log10(pmax(pvalue, 1e-300)), 30)
)
write_csv_safe(plot_data, file.path(out_dir, "GSEA_selected_positive_plot_data.csv"))

p <- ggplot(plot_data, aes(NES, Description)) +
  geom_segment(aes(x = 0, xend = NES, yend = Description), color = "#D4D4D4") +
  geom_point(aes(shape = Treatment, size = CoreGeneRatio, fill = logP), color = "black") +
  facet_wrap(~Organ, nrow = 1) +
  scale_shape_manual(values = c(LNP = 22),labels=plot_label) +
  scale_fill_gradient(low = "#FFE79A", high = "#D7191C") +
  scale_size_continuous(range = c(2.5, 9)) +
  labs(x = "NES", y = NULL, fill = "-log10(P)", size = "Core gene ratio") +
  theme_classic(base_size = 12) + theme(strip.text = element_text(face = "bold", size = 17))
save_gg(p, file.path(out_dir, "Figure_08_GSEA_NES_positive"), 11, 7)


} else message("KEGG GSEA 需要 MOUSE_KEGG_GMT 指向小鼠 KEGG GMT 文件。")

######## 11_figure_key_gene_TPM_heatmap ###########

#!/usr/bin/env Rscript

# Figure 9: selected single-gene heatmap.
# Uses the original TPM matrix directly: log2(mean TPM of replicates + 1).
# No batch removal and no row Z-score.
# Edit key_gene_sets to change genes/categories.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
#source(file.path(script_dir, "00_config_and_helpers.txt"))

out_dir <- file.path(output_root, "10_key_gene_TPM_heatmap")
dat <- load_prepared()

key_gene_sets <- list(
  Apoptosis = c("Aif1", "Ifitm3", "Socs3"),
  Senescence = c("Cd68", "Lcn2", "Saa1", "Trp53"),
  Inflammation = c("Ccl6", "Cxcl1", "Cxcl9", "Hp", "Orm2", "Pf4",
                   "S100a6", "S100a8", "S100a9"),
  Tumorigenesis = c("Bst2", "Lrg1", "Ly6a","Ctnnb1")
)
annotation <- bind_rows(lapply(names(key_gene_sets), function(category) {
  data.frame(gene = key_gene_sets[[category]], Category = category)
}))
genes <- intersect(annotation$gene, rownames(dat$tpm))
write_csv_safe(data.frame(gene = setdiff(annotation$gene, genes)),
               file.path(out_dir, "missing_key_genes.csv"))

plot_data <- bind_rows(lapply(organ_order, function(organ_name) {
  meta <- dat$metadata %>% filter(Organ == organ_name)
  bind_rows(lapply(treatment_order, function(trt) {
    mean_tpm <- rowMeans(dat$tpm[genes, meta$Sample[meta$Treatment == trt], drop = FALSE])
    data.frame(gene = names(mean_tpm), Organ = organ_name, Treatment = trt,
               MeanTPM = mean_tpm, Log2TPM = log2(mean_tpm + 1))
  }))
})) %>% left_join(annotation, by = "gene") %>% mutate(
  Treatment = factor(Treatment, levels = treatment_order),
  Organ = factor(Organ, levels = organ_order),
  Category = factor(Category, levels = names(key_gene_sets)),
  gene = factor(gene, levels = annotation$gene)
)
write_csv_safe(plot_data, file.path(out_dir, "key_gene_TPM_group_means_and_plot_data.csv"))

p <- ggplot(plot_data, aes(gene, Treatment, fill = Log2TPM)) + scale_y_discrete(labels=plot_label) +
  geom_tile(color = "white", linewidth = 0.8) +
  facet_grid(Organ ~ Category, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = c("#5B4A9B", "white", "#E63B22"),
                       values = scales::rescale(c(0, 5, 10)),
                       limits = c(0, 10), oob = scales::squish) +
  labs(x = NULL, y = NULL, fill = "log2(mean TPM + 1)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        strip.text = element_text(face = "bold", size = 13),
        panel.spacing = grid::unit(0.35, "lines"))
save_gg(p, file.path(out_dir, "Figure_09_key_gene_TPM_heatmap"), 14, 6.5)




