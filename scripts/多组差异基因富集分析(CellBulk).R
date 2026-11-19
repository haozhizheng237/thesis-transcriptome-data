# 可用 THESIS_DATA 和 LNP_RESULTS 指定输入和输出目录。
data_root <- path.expand(Sys.getenv("THESIS_DATA", "~/Downloads/毕业论文数据"))
result_root <- path.expand(Sys.getenv("LNP_RESULTS", "LNP_results"))
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
# 仅在图形显示层使用 eLNP。
plot_label <- function(x) gsub("(?<![A-Za-z])LNP", "eLNP", as.character(x), perl = TRUE)

# 按样本名选列，保留原始计数，用批次参与建模。
fit_counts <- function(tab, meta, contrasts, out, quantile_voom = FALSE) {
  stopifnot(!anyDuplicated(meta$sample), all(meta$sample %in% names(tab)))
  gene <- tab[[intersect(c("gene", "gene_id"), names(tab))[1]]]
  counts <- as.matrix(tab[, meta$sample, drop = FALSE])
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(colSums(counts) > 0),
            !anyNA(gene), all(nzchar(gene)))
  counts <- rowsum(counts, gene, reorder = FALSE)
  meta$batch <- factor(meta$batch)
  design <- if (nlevels(meta$batch) > 1) model.matrix(~ 0 + group + batch, meta) else model.matrix(~ 0 + group, meta)
  colnames(design) <- sub("^group", "", colnames(design))
  stopifnot(qr(design)$rank == ncol(design), nrow(design) > ncol(design))
  y <- edgeR::DGEList(counts)
  y <- y[edgeR::filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  v <- limma::voom(y, design, normalize.method = if (quantile_voom) "quantile" else "none")
  model <- limma::lmFit(v, design)
  cm <- limma::makeContrasts(contrasts = unname(contrasts), levels = design)
  colnames(cm) <- names(contrasts)
  model <- limma::eBayes(limma::contrasts.fit(model, cm), robust = !quantile_voom)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(meta, file.path(out, "samples.csv"), row.names = FALSE)
  expression <- v$E
  if (nlevels(meta$batch) > 1)
    expression <- limma::removeBatchEffect(expression, batch = meta$batch, design = model.matrix(~ group, meta))
  pcs <- prcomp(t(expression))$x[, 1:2, drop = FALSE]
  write.csv(cbind(meta, pcs), file.path(out, "PCA.csv"), row.names = FALSE)
  write.csv(cor(expression), file.path(out, "sample_correlations.csv"))
  pdf(file.path(out, "PCA.pdf"))
  plot(pcs, col = as.integer(meta$group), pch = 19)
  legend("topright", plot_label(levels(meta$group)), col = seq_along(levels(meta$group)), pch = 19)
  dev.off()
  result <- lapply(names(contrasts), function(key) {
    t <- limma::topTable(model, coef = key, number = Inf)
    t$gene <- rownames(t)
    t$significant <- if (quantile_voom) {
      ifelse(t$P.Value < 0.05 & abs(t$logFC) >= 1.5, ifelse(t$logFC > 0,"Up","Down"), "NS")
    } else {
      ifelse(t$adj.P.Val < 0.05 & abs(t$logFC) > 1, ifelse(t$logFC > 0,"Up","Down"), "NS")
    }
    write.csv(t, file.path(out, paste0(key, ".csv")), row.names = FALSE)
    t
  })
  list(tables = setNames(result, names(contrasts)), expression = expression)
}
input_dir <- file.path(data_root, "Cell_Bulk")
meta_all <- read.csv(file.path(input_dir, "three_cellline_sample_metadata.csv"), check.names = FALSE)
files <- c(AML12 = "eLNP_AML12_raw_counts.csv", `TCMK-1` = "eLNP_TCMK_1_raw_counts.csv", `293T` = "eLNP_293t_raw_counts.csv")
for (model in names(files)) {
  m <- subset(meta_all, cell_line == model & treatment %in% c("Control","LNP"))
  m$sample <- m$count_column
  m$group <- factor(m$treatment, c("Control","LNP"))
  stopifnot(all(table(m$group) >= 2))
  fit_counts(read.csv(file.path(input_dir, files[[model]]), check.names = FALSE), m,
             c(LNP_vs_Control = "LNP-Control"), file.path(result_root, "Cell_Bulk", model), TRUE)
}
###################差异基因富集分析#####################

#!/usr/bin/env Rscript

# Shared functions for CellBulk GO BP scripts 07 and 08.
# Run 07/08, not this file directly.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
})

input_dir <- file.path(result_root, "Cell_Bulk")
treatments <- "LNP"
overlap_pvalue_cutoff <- 0.05
group_colors <- c(LNP = "#FF7F0E")

cell_config <- data.frame(
  CellLine = c("TCMK-1", "AML12", "293T"),
  Species = c("mouse", "mouse", "human"),
  DEGPattern = c(
    "TCMK-1/%s_vs_Control.csv",
    "AML12/%s_vs_Control.csv",
    "293T/%s_vs_Control.csv"
  ),
  stringsAsFactors = FALSE
)

write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(x), path, na = "")
}

save_gg <- function(p, stem, width, height, dpi = 350) {
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height)
  ggplot2::ggsave(
    paste0(stem, ".png"), p, width = width, height = height,
    dpi = dpi, bg = "white"
  )
}

find_column_index <- function(dat, candidates, required_name, path) {
  nm <- colnames(dat)
  hit <- which(tolower(nm) %in% tolower(candidates))
  if (!length(hit)) stop("Cannot find ", required_name, " in ", path)
  hit[1]
}

read_deg_table <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  dat <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal") %>%
    as.data.frame(check.names = FALSE)
  
  sig_col <- find_column_index(
    dat, c("significant"), "the 'significant' column", path
  )
  
  # These DEG files store gene symbols in the first unnamed column.
  if (is.na(colnames(dat)[1]) || colnames(dat)[1] == "") {
    gene_col <- 1L
  } else {
    candidates <- c(
      "gene", "symbol", "gene_symbol", "genesymbol", "X", "X1", "...1"
    )
    hit <- which(tolower(colnames(dat)) %in% tolower(candidates))
    gene_col <- if (length(hit)) hit[1] else 1L
  }
  
  gene <- trimws(as.character(dat[[gene_col]]))
  significant <- trimws(as.character(dat[[sig_col]]))
  keep <- !is.na(gene) & gene != ""
  data.frame(gene = gene[keep], significant = significant[keep])
}

map_to_entrez <- function(genes, species) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != ""]
  if (!length(genes)) return(character())
  
  is_ensembl <- grepl("^ENS[A-Z]*G[0-9]+(?:\\.[0-9]+)?$", toupper(genes))
  from_type <- if (mean(is_ensembl) >= 0.5) "ENSEMBL" else "SYMBOL"
  if (from_type == "ENSEMBL") genes <- sub("\\.[0-9]+$", "", genes)
  orgdb <- if (species == "mouse") {
    org.Mm.eg.db::org.Mm.eg.db
  } else {
    org.Hs.eg.db::org.Hs.eg.db
  }
  
  mapped <- tryCatch(
    suppressMessages(clusterProfiler::bitr(
      genes, fromType = from_type, toType = "ENTREZID", OrgDb = orgdb
    )),
    error = function(e) {
      warning("Gene-ID mapping failed: ", conditionMessage(e))
      data.frame()
    }
  )
  if (is.null(mapped) || !nrow(mapped)) character() else unique(mapped$ENTREZID)
}

empty_go_result <- function() {
  data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), RichFactor = numeric(), FoldEnrichment = numeric(),
    zScore = numeric(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    CellLine = character(), Treatment = character(), Direction = character(),
    Ontology = character(), stringsAsFactors = FALSE
  )
}

run_one_enrichment <- function(tab, direction, species, cell_line, treatment) {
  selected <- tab$gene[tolower(tab$significant) == tolower(direction)]
  selected_entrez <- map_to_entrez(selected, species)
  universe_entrez <- map_to_entrez(tab$gene, species)
  
  if (length(selected_entrez) < 3 || length(universe_entrez) < 10) {
    warning(cell_line, " ", treatment, " ", direction,
            ": too few mapped genes for GO BP enrichment")
    return(empty_go_result())
  }
  
  orgdb <- if (species == "mouse") {
    org.Mm.eg.db::org.Mm.eg.db
  } else {
    org.Hs.eg.db::org.Hs.eg.db
  }
  eg <- suppressMessages(clusterProfiler::enrichGO(
    gene = selected_entrez,
    universe = universe_entrez,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ))
  result <- if (is.null(eg)) data.frame() else as.data.frame(eg)
  if (!nrow(result)) return(empty_go_result())
  
  result %>% mutate(
    CellLine = cell_line,
    Treatment = treatment,
    Direction = direction,
    Ontology = "BP"
  )
}

run_cellline_go_bp <- function(
    direction, requested_bp_terms, excluded_non_bp_terms = character(),
    output_subdir, figure_number, heat_color) {
  
  out_dir <- file.path(input_dir, output_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (length(excluded_non_bp_terms)) {
    write_csv_safe(
      data.frame(
        RequestedTerm = excluded_non_bp_terms,
        Reason = "Not GO Biological Process; replaced using TCMK-1 BP ranking"
      ),
      file.path(out_dir, "excluded_non_BP_requested_terms.csv")
    )
  }
  
  go_results <- list()
  mapping_summary <- list()
  
  for (i in seq_len(nrow(cell_config))) {
    cfg <- cell_config[i, ]
    for (trt in treatments) {
      path <- file.path(input_dir, sprintf(cfg$DEGPattern, trt))
      tab <- read_deg_table(path)
      key <- paste(cfg$CellLine, trt, sep = "_")
      direction_genes <- tab$gene[tolower(tab$significant) == tolower(direction)]
      mapped_direction <- map_to_entrez(direction_genes, cfg$Species)
      mapped_universe <- map_to_entrez(tab$gene, cfg$Species)
      
      write_csv_safe(
        data.frame(gene = unique(direction_genes)),
        file.path(out_dir, paste0(cfg$CellLine, "_", trt, "_", direction, "_DEGs.csv"))
      )
      mapping_summary[[key]] <- data.frame(
        CellLine = cfg$CellLine, Treatment = trt, Direction = direction,
        Input_DEGs = length(unique(direction_genes)),
        Mapped_DEG_Entrez = length(mapped_direction),
        Universe_genes = length(unique(tab$gene)),
        Mapped_universe_Entrez = length(mapped_universe)
      )
      
      result <- run_one_enrichment(
        tab, direction, cfg$Species, cfg$CellLine, trt
      )
      go_results[[key]] <- result
      write_csv_safe(
        result,
        file.path(out_dir, paste0(cfg$CellLine, "_", trt, "_GO_BP_all.csv"))
      )
    }
  }
  
  go_all <- dplyr::bind_rows(go_results)
  write_csv_safe(go_all, file.path(out_dir, paste0("GO_BP_", tolower(direction), "_all_comparisons.csv")))
  write_csv_safe(dplyr::bind_rows(mapping_summary), file.path(out_dir, "gene_mapping_summary.csv"))
  
  # Fill missing slots only from TCMK-1 GO BP terms, ranked by the smallest raw
  # p value observed among SiO2/LNP/EXO. This is also used to replace explicitly
  # excluded non-BP terms supplied by the user.
  tcmk_ranked <- go_all %>%
    dplyr::filter(CellLine == "TCMK-1", is.finite(pvalue)) %>%
    dplyr::group_by(ID, Description) %>%
    dplyr::summarise(best_pvalue = min(pvalue), .groups = "drop") %>%
    dplyr::arrange(best_pvalue, Description)
  
  n_needed <- 10L - length(requested_bp_terms)
  replacement_terms <- character()
  if (n_needed > 0) {
    replacement_terms <- tcmk_ranked %>%
      dplyr::filter(!Description %in% requested_bp_terms) %>%
      dplyr::slice_head(n = n_needed) %>%
      dplyr::pull(Description)
    if (length(replacement_terms) < n_needed) {
      message("Available GO BP terms retained: ", length(replacement_terms))
    }
  }
  selected_terms <- c(requested_bp_terms, replacement_terms)
  selected_info <- data.frame(
    DisplayOrder = seq_along(selected_terms),
    Description = selected_terms,
    SelectionSource = c(
      rep("User-specified GO BP", length(requested_bp_terms)),
      rep("TCMK-1 BP replacement ranked by raw pvalue", length(replacement_terms))
    )
  ) %>%
    dplyr::left_join(
      tcmk_ranked %>% dplyr::select(Description, TCMK1_best_raw_pvalue = best_pvalue),
      by = "Description"
    )
  write_csv_safe(selected_info, file.path(out_dir, "selected_GO_BP_pathways.csv"))
  
  best_per_cell <- go_all %>%
    dplyr::filter(Description %in% selected_terms, is.finite(pvalue)) %>%
    dplyr::group_by(Description, CellLine, Treatment) %>%
    dplyr::slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(Description, CellLine, Treatment, ID, pvalue, p.adjust,
                  GeneRatio, Count, geneID, Ontology)
  
  plot_data <- tidyr::expand_grid(
    Description = selected_terms,
    CellLine = cell_config$CellLine,
    Treatment = treatments
  ) %>%
    dplyr::left_join(best_per_cell, by = c("Description", "CellLine", "Treatment")) %>%
    dplyr::mutate(
      score = ifelse(is.na(pvalue), 0, pmin(-log10(pmax(pvalue, 1e-300)), 30)),
      rawP_lt_0.05 = !is.na(pvalue) & pvalue < overlap_pvalue_cutoff,
      Description = factor(Description, levels = rev(selected_terms)),
      CellLine = factor(CellLine, levels = cell_config$CellLine),
      Treatment = factor(Treatment, levels = treatments)
    )
  write_csv_safe(plot_data, file.path(out_dir, paste0("GO_BP_", tolower(direction), "_selected_plot_data.csv")))
  
  p_heat <- ggplot(plot_data, aes(Treatment, Description, fill = score)) + scale_x_discrete(labels=plot_label) +
    geom_tile(color = "white", linewidth = 0.8) +
    facet_wrap(~CellLine, nrow = 1) +
    scale_fill_gradient(low = "white", high = heat_color) +
    labs(x = NULL, y = NULL, fill = "-log10(P)") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 15)
    )
  save_gg(
    p_heat,
    file.path(out_dir, paste0("Figure_", figure_number, "_CellLines_GO_BP_", tolower(direction), "_heatmap")),
    12, 7
  )
  
  message("Analysis complete: ", out_dir)
}

#### Up DEGs部分
#!/usr/bin/env Rscript

# Cell-line GO BP analysis for significant UP-regulated DEGs.
# GO-term overlap threshold: raw enrichment pvalue < 0.05.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) {
  #dirname(normalizePath(sub("^--file=", "", script_arg[1])))
#} else {
  #getwd()
#}
#source(file.path(script_dir, "00_CellLines_GO_BP_common.txt"))

requested_bp_terms <- c(
  "regulation of response to cytokine stimulus",
  "regulation of cytokine-mediated signaling pathway",
  "regulation of apoptotic signaling pathway",
  "positive regulation of leukocyte mediated immunity",
  "positive regulation of inflammatory response",
  "positive regulation of immune effector process",
  "inflammatory response to antigenic stimulus",
  "cytokine production involved in inflammatory response",
  "cellular response to oxidative stress",
  "acute inflammatory response"
)

run_cellline_go_bp(
  direction = "Up",
  requested_bp_terms = requested_bp_terms,
  excluded_non_bp_terms = character(),
  output_subdir = "07_CellLines_GO_BP_up",
  figure_number = "07",
  heat_color = "#A7191F"
)

##### Down DEGs 部分
#!/usr/bin/env Rscript

# Cell-line GO BP analysis for significant DOWN-regulated DEGs.
# GO-term overlap threshold: raw enrichment pvalue < 0.05.
# Two requested MF terms are excluded and replaced by the top-ranked TCMK-1
# GO BP terms according to the smallest raw pvalue across three treatments.

#script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#script_dir <- if (length(script_arg)) {
  #dirname(normalizePath(sub("^--file=", "", script_arg[1])))
#} else {
  #getwd()
#}
#source(file.path(script_dir, "00_CellLines_GO_BP_common.txt"))

requested_bp_terms <- c(
  "regulation of response to extracellular stimulus",
  "positive regulation of monocyte chemotaxis",
  "positive regulation of macrophage migration",
  "intrinsic apoptotic signaling pathway",
  "interleukin-12 production",
  "cytokine-mediated signaling pathway",
  "cellular response to extracellular stimulus",
  "apoptotic cell clearance"
)

excluded_non_bp_terms <- c(
  "histone deacetylase activity",
  "deacetylase activity"
)

run_cellline_go_bp(
  direction = "Down",
  requested_bp_terms = requested_bp_terms,
  excluded_non_bp_terms = excluded_non_bp_terms,
  output_subdir = "08_CellLines_GO_BP_down",
  figure_number = "08",
  heat_color = "#244A7C"
)





contrast_treatments <- "LNP"
read_deg_table <- function(path) {
  if (!file.exists(path)) stop("DEG file not found: ", path)
  dat <- readr::read_csv(path, show_col_types = FALSE, name_repair = "minimal") %>%
    as.data.frame(check.names = FALSE)
  nm <- colnames(dat)
  sig_hit <- which(tolower(nm) == "significant")
  logfc_hit <- which(tolower(nm) %in% c("logfc", "log2fc", "log2foldchange"))
  if (!length(sig_hit)) stop("No significant column in ", path)
  if (!length(logfc_hit)) stop("No logFC column in ", path)
  if (is.na(nm[1]) || nm[1] == "") {
    gene_col <- 1L
  } else {
    gene_hit <- which(tolower(nm) %in% c(
      "gene", "symbol", "gene_symbol", "genesymbol", "x", "x1", "...1"
    ))
    gene_col <- if (length(gene_hit)) gene_hit[1] else 1L
  }
  out <- data.frame(
    gene = trimws(as.character(dat[[gene_col]])),
    logFC = suppressWarnings(as.numeric(dat[[logfc_hit[1]]])),
    significant = trimws(as.character(dat[[sig_hit[1]]]))
  )
  out[!is.na(out$gene) & out$gene != "", , drop = FALSE]
}

gene_key <- function(x) {
  x <- trimws(as.character(x))
  key <- toupper(x)
  ens <- grepl("^ENS[A-Z]*G[0-9]+\\.[0-9]+$", key)
  key[ens] <- sub("\\.[0-9]+$", "", key[ens])
  key
}

match_to_expression <- function(genes, expression_genes) {
  exp_keys <- gene_key(expression_genes)
  keep <- !duplicated(exp_keys) & exp_keys != "" & !is.na(exp_keys)
  lookup <- setNames(expression_genes[keep], exp_keys[keep])
  mapped <- unname(lookup[gene_key(genes)])
  unique(mapped[!is.na(mapped) & mapped != ""])
}

map_gene_table <- function(genes, species) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != ""]
  if (!length(genes)) return(data.frame(gene = character(), ENTREZID = character()))
  ens <- grepl("^ENS[A-Z]*G[0-9]+(?:\\.[0-9]+)?$", toupper(genes))
  from_type <- if (mean(ens) >= 0.5) "ENSEMBL" else "SYMBOL"
  keys <- if (from_type == "ENSEMBL") sub("\\.[0-9]+$", "", genes) else genes
  orgdb <- if (species == "mouse") org.Mm.eg.db::org.Mm.eg.db else org.Hs.eg.db::org.Hs.eg.db
  mapped <- tryCatch(
    suppressMessages(clusterProfiler::bitr(
      keys, fromType = from_type, toType = "ENTREZID", OrgDb = orgdb
    )), error = function(e) data.frame()
  )
  if (!nrow(mapped)) return(data.frame(gene = character(), ENTREZID = character()))
  names(mapped)[1] <- "map_key"
  source <- data.frame(gene = genes, map_key = keys)
  source %>% inner_join(mapped, by = "map_key") %>%
    dplyr::select(gene, ENTREZID) %>% distinct()
}
read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 3]
  labels <- vapply(fields, `[[`, character(1), 1)
  ids <- sub(" .*", "", labels)
  desc <- sub("^[^ ]+ +", "", labels)
  list(
    term2gene = bind_rows(lapply(seq_along(fields), function(i) {
      data.frame(term = ids[i], gene = unique(fields[[i]][-(1:2)]))
    })),
    term2name = data.frame(term = ids, name = desc) %>% distinct(term, .keep_all = TRUE)
  )
}

load_kegg_collection <- function(species) {
  file_candidates <- if (species == "mouse") c(
    file.path(input_dir, "kegg_mouse_KEGGREST.gmt"),
    file.path(input_dir, "kegg_mouse_gageData.gmt"),
    file.path(dirname(input_dir), "Liver_Spleen_Bulk", "kegg_mouse_gageData.gmt")
  ) else c(
    file.path(input_dir, "kegg_human_KEGGREST.gmt"),
    file.path(input_dir, "kegg_human_gageData.gmt"),
    file.path(dirname(input_dir), "Liver_Spleen_Bulk", "kegg_human_gageData.gmt")
  )
  exists_and_nonempty <- file.exists(file_candidates) &
    !is.na(file.info(file_candidates)$size) & file.info(file_candidates)$size > 0
  existing <- file_candidates[exists_and_nonempty]
  if (length(existing)) return(read_gmt(existing[1]))
  
  if (requireNamespace("gageData", quietly = TRUE)) {
    dataset <- if (species == "mouse") "kegg.sets.mm" else "kegg.sets.hs"
    env <- new.env(parent = emptyenv())
    utils::data(list = dataset, package = "gageData", envir = env)
    sets <- get(dataset, envir = env)
    labels <- names(sets)
    ids <- sub(" .*", "", labels)
    desc <- sub("^[^ ]+ +", "", labels)
    collection <- list(
      term2gene = bind_rows(lapply(seq_along(sets), function(i) {
        data.frame(term = ids[i], gene = unique(as.character(sets[[i]])))
      })),
      term2name = data.frame(term = ids, name = desc) %>%
        distinct(term, .keep_all = TRUE)
    )
    # Persist the package data as a local GMT, matching the Liver/Spleen setup.
    cache_path <- file.path(
      input_dir,
      if (species == "mouse") "kegg_mouse_gageData.gmt" else "kegg_human_gageData.gmt"
    )
    gmt_lines <- vapply(seq_along(sets), function(j) {
      paste(c(labels[j], "gageData", unique(as.character(sets[[j]]))), collapse = "\t")
    }, character(1))
    writeLines(gmt_lines, cache_path)
    message("Cached local gageData KEGG GMT: ", cache_path)
    return(collection)
  }
  
  # For human data, stop immediately instead of performing slow OrgDb scans or
  # waiting for a blocked KEGG REST connection. Run the supplied one-time GMT
  # preparation script first.
  if (species == "human") {
    warning(
      "Human KEGG GMT is absent and gageData is not installed. Run ",
      "with a human GMT for human KEGG GSEA."
    )
    return(NULL)
  }
  
  # Offline fallback 1: MSigDB KEGG collections, if msigdbr is installed.
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    message("Using installed msigdbr as the offline KEGG source for ", species)
    msig <- tryCatch({
      call_args <- list(species = if (species == "mouse") "Mus musculus" else "Homo sapiens")
      if ("collection" %in% names(formals(msigdbr::msigdbr))) {
        call_args$collection <- "C2"
      } else {
        call_args$category <- "C2"
      }
      do.call(msigdbr::msigdbr, call_args)
    }, error = function(e) data.frame())
    if (nrow(msig)) {
      subcollection_col <- intersect(c("gs_subcollection", "gs_subcat"), colnames(msig))
      if (length(subcollection_col)) {
        msig <- msig[grepl("KEGG", msig[[subcollection_col[1]]], ignore.case = TRUE), , drop = FALSE]
      } else {
        msig <- msig[grepl("KEGG", msig$gs_name, ignore.case = TRUE), , drop = FALSE]
      }
      entrez_col <- intersect(c("ncbi_gene", "entrez_gene"), colnames(msig))
      if (nrow(msig) && length(entrez_col)) {
        term2gene <- data.frame(
          term = msig$gs_name,
          gene = as.character(msig[[entrez_col[1]]]),
          stringsAsFactors = FALSE
        ) %>% filter(!is.na(gene), gene != "") %>% distinct()
        term2name <- data.frame(
          term = unique(msig$gs_name),
          name = tools::toTitleCase(tolower(gsub("^KEGG(_MEDICUS|_LEGACY)?_", "", unique(msig$gs_name)))),
          stringsAsFactors = FALSE
        )
        if (nrow(term2gene)) return(list(term2gene = term2gene, term2name = term2name))
      }
    }
  }
  
  # Offline fallback 2: legacy KEGG PATH mappings bundled in some OrgDb
  # releases. When pathway descriptions are unavailable, retain the pathway ID;
  # a dictionary supplies names for the priority pathways used by script 10.
  orgdb <- if (species == "mouse") org.Mm.eg.db::org.Mm.eg.db else org.Hs.eg.db::org.Hs.eg.db
  if ("PATH" %in% AnnotationDbi::columns(orgdb)) {
    message("Using OrgDb PATH mappings as the offline KEGG source for ", species)
    path_map <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        orgdb,
        keys = AnnotationDbi::keys(orgdb, keytype = "ENTREZID"),
        keytype = "ENTREZID", columns = "PATH"
      )), error = function(e) data.frame()
    )
    if (nrow(path_map) && all(c("ENTREZID", "PATH") %in% colnames(path_map))) {
      path_map <- path_map[!is.na(path_map$PATH), , drop = FALSE]
      organism <- if (species == "mouse") "mmu" else "hsa"
      term2gene <- data.frame(
        term = paste0(organism, path_map$PATH),
        gene = as.character(path_map$ENTREZID), stringsAsFactors = FALSE
      ) %>% distinct()
      known_names <- c(
        `05225` = "Hepatocellular carcinoma",
        `04148` = "Efferocytosis",
        `04657` = "IL-17 signaling pathway",
        `04060` = "Cytokine-cytokine receptor interaction",
        `04650` = "Natural killer cell mediated cytotoxicity",
        `04210` = "Apoptosis",
        `05321` = "Inflammatory bowel disease",
        `04145` = "Phagosome",
        `04668` = "TNF signaling pathway",
        `04062` = "Chemokine signaling pathway",
        `04142` = "Lysosome"
      )
      path_ids <- unique(as.character(path_map$PATH))
      path_names <- unname(known_names[path_ids])
      path_names[is.na(path_names)] <- paste("KEGG pathway", path_ids[is.na(path_names)])
      term2name <- data.frame(
        term = paste0(organism, path_ids), name = path_names,
        stringsAsFactors = FALSE
      )
      if (nrow(term2gene)) return(list(term2gene = term2gene, term2name = term2name))
    }
  }
  
  warning("Missing offline KEGG collection for ", species)
  NULL
}

out_dir <- file.path(input_dir, "10_CellLines_GSEA_logFC_ranked_NES_positive")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

priority_terms <- unique(c(
  "Hepatocellular carcinoma",
  "Efferocytosis",
  "IL-17 signaling pathway",
  "Cytokine-cytokine receptor interaction",
  "Natural killer cell mediated cytotoxicity",
  "Apoptosis",
  "Inflammatory bowel disease",
  "TNF signaling pathway",
  "Chemokine signaling pathway",
  "Lysosome",
  "Efferocytosis",
  "TNF signaling pathway",
  "Chemokine signaling pathway",
  "B cell receptor signaling pathway",
  "Cytokine-cytokine receptor interaction",
  "IL-17 signaling pathway",
  "Viral protein interaction with cytokine and cytokine receptor",
  "NF-kappa B signaling pathway",
  "Chemical carcinogenesis-DNA adducts",
  "Inflammatory mediator regulation of TRP Channels"
))

empty_gsea <- function() {
  data.frame(
    ID = character(), Description = character(), setSize = integer(),
    enrichmentScore = numeric(), NES = numeric(), pvalue = numeric(),
    p.adjust = numeric(), qvalue = numeric(), rank = integer(),
    leading_edge = character(), core_enrichment = character(),
    CellLine = character(), Treatment = character(), Species = character(),
    CoreGeneCount = integer(), CoreGeneRatio = numeric()
  )
}

results <- list()
for (i in seq_len(nrow(cell_config))) {
  cfg <- cell_config[i, ]
  kegg <- load_kegg_collection(cfg$Species)
  
  if (is.null(kegg)) {
    warning(cfg$CellLine, ": GSEA skipped because no offline KEGG collection is available")
    for (trt in contrast_treatments) {
      key <- paste(cfg$CellLine, trt, sep = "_")
      results[[key]] <- empty_gsea()
      write_csv_safe(
        data.frame(
          CellLine = cfg$CellLine, Treatment = trt,
          Status = "Skipped: install gageData or provide a human KEGG GMT"
        ),
        file.path(out_dir, paste0(cfg$CellLine, "_", trt, "_GSEA_status.csv"))
      )
    }
    next
  }
  
  for (trt in contrast_treatments) {
    message("GSEA: ", cfg$CellLine, " ", trt, " vs Control")
    path <- file.path(input_dir, sprintf(cfg$DEGPattern, trt))
    tab <- read_deg_table(path)
    mapped <- map_gene_table(tab$gene, cfg$Species)
    ranked <- tab %>%
      inner_join(mapped, by = "gene") %>%
      filter(is.finite(logFC), !is.na(ENTREZID)) %>%
      group_by(ENTREZID) %>%
      slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(logFC))
    
    if (nrow(ranked) < 100) {
      stop(cfg$CellLine, " ", trt, ": fewer than 100 mapped genes for GSEA")
    }
    gene_list <- ranked$logFC
    names(gene_list) <- ranked$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)
    write_csv_safe(
      ranked %>% dplyr::select(gene, ENTREZID, logFC),
      file.path(out_dir, paste0(cfg$CellLine, "_", trt, "_GSEA_logFC_ranked_gene_list.csv"))
    )
    
    set.seed(20260822)
    gsea <- suppressMessages(clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = kegg$term2gene,
      TERM2NAME = kegg$term2name,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE
    ))
    result <- if (is.null(gsea)) data.frame() else as.data.frame(gsea)
    if (nrow(result)) {
      result <- result %>% mutate(
        CellLine = cfg$CellLine,
        Treatment = trt,
        Species = cfg$Species,
        CoreGeneCount = ifelse(
          is.na(core_enrichment) | core_enrichment == "", 0L,
          lengths(strsplit(core_enrichment, "/", fixed = TRUE))
        ),
        CoreGeneRatio = CoreGeneCount / setSize
      )
    } else {
      result <- empty_gsea()
    }
    key <- paste(cfg$CellLine, trt, sep = "_")
    results[[key]] <- result
    write_csv_safe(
      result,
      file.path(out_dir, paste0(cfg$CellLine, "_", trt, "_GSEA_all.csv"))
    )
  }
}

all_gsea <- bind_rows(results)
positive <- all_gsea %>% filter(is.finite(NES), NES > 0)
write_csv_safe(all_gsea, file.path(out_dir, "CellLines_GSEA_all_comparisons.csv"))
write_csv_safe(positive, file.path(out_dir, "CellLines_GSEA_NES_positive_only.csv"))
if (nrow(positive)) {

available <- positive %>%
  group_by(Description) %>%
  summarise(best_padjust = min(p.adjust, na.rm = TRUE),
            best_pvalue = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_padjust, best_pvalue)
selected <- intersect(priority_terms, available$Description)
if (length(selected) < 10) {
  selected <- c(
    selected,
    head(setdiff(available$Description, selected), 10 - length(selected))
  )
}
selected <- head(selected, 10)
write_csv_safe(
  data.frame(
    DisplayOrder = seq_along(selected),
    Description = selected,
    SelectionSource = ifelse(selected %in% priority_terms,
                             "Priority pathway", "Best available NES-positive pathway")
  ) %>% left_join(available, by = "Description"),
  file.path(out_dir, "GSEA_selected_NES_positive_pathways.csv")
)

plot_data <- positive %>%
  filter(Description %in% selected) %>%
  mutate(
    Description = factor(Description, levels = rev(selected)),
    Treatment = factor(Treatment, levels = contrast_treatments),
    CellLine = factor(CellLine, levels = cell_config$CellLine),
    logP = pmin(-log10(pmax(pvalue, 1e-300)), 30)
  )
write_csv_safe(plot_data, file.path(out_dir, "GSEA_selected_NES_positive_plot_data.csv"))

p <- ggplot(plot_data, aes(NES, Description)) +
  geom_segment(aes(x = 0, xend = NES, yend = Description), color = "#D4D4D4") +
  geom_point(
    aes(shape = Treatment, size = CoreGeneRatio, fill = logP),
    color = "black"
  ) +
  facet_wrap(~CellLine, nrow = 1) +
  scale_shape_manual(values = c(LNP = 22),labels=plot_label) +
  scale_fill_gradient(low = "#FFE79A", high = "#D7191C") +
  scale_size_continuous(range = c(2.5, 9)) +
  labs(x = "NES", y = NULL, fill = "-log10(P)", size = "Core gene ratio") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 16))
save_gg(p, file.path(out_dir, "Figure_10_CellLines_GSEA_NES_positive"), 14, 7)

message("Analysis complete: ", out_dir)


}

library(ggplot2)
plot_assay <- function(values, label) {
  df <- stack(values)
  df$ind <- factor(df$ind, names(values))
  p <- ggplot(df, aes(ind, values)) + scale_x_discrete(labels=plot_label) + stat_summary(fun = mean, geom = "col", fill = "grey80") +
    stat_summary(fun.data = function(x) data.frame(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x)), geom="errorbar", width=.2) +
    geom_jitter(width=.08) + labs(x=NULL,y="Cell viability (%)",title=label) + theme_classic()
  ggsave(file.path(result_root, "Cell_Bulk", paste0(label,".pdf")), p)
}
pbs_vals <- c(101.08, 119.31, 92.57, 93.78, 93.26)
lnp_1e8    <- c(98.35, 97.74, 107.85, 100.43)
lnp_1e9    <- c(94.03, 88.67, 91.21, 93.69)
lnp_1e10   <- c(74.53, 86.73, 72.84, 88.30)
lnp_1e11   <- c(59.69, 55.98, 58.93, 53.556)
lnp_1e12   <- c(43.78, 43.17, 46.78)
control_0h <- c(96.6, 96.6, 92.57, 92.57, 104.14)
lnp_12h <- c(72.74, 82.51, 81.65, 72.45, 75.61)
lnp_24h <- c(74.37, 70.37, 75.99, 77.69, 84.59)
lnp_36h <- c(53.23, 55.20, 50.60, 57.50, 56.52)
lnp_48h <- c(21.55, 31.28, 25.14, 28.89, 22.15)
plot_assay(list(PBS=pbs_vals, LNP_1e8=lnp_1e8, LNP_1e9=lnp_1e9, LNP_1e10=lnp_1e10, LNP_1e11=lnp_1e11, LNP_1e12=lnp_1e12), "TCMK-1 concentration")
plot_assay(list(Control_0h=control_0h, LNP_12h=lnp_12h, LNP_24h=lnp_24h, LNP_36h=lnp_36h, LNP_48h=lnp_48h), "TCMK-1 time")

key_gene_sets <- list(
  Apoptosis = c(
    "Dap", "Ddit3", "Dffa", "Gadd45a", "Gadd45b",
    "Gadd45g", "Hmox1", "Pmaip1", "Rhob", "Sqstm1"
  ),
  `Cellular Senescence` = c("Cdkn1a", "Cdkn2a", "Hmga1", "Mapk3"),
  `Inflammatory Response` = c(
    "Ccl2", "Cd68", "Cxcl1", "Cxcl10", "Cxcl5", "Pvr", "Slc11a2"
  ),
  Cancer = c("Atf3", "Lcn2", "Nupr1", "Plscr1", "Rhoc", "Tgif1")
)

# 原 TPM 矩阵可通过 CELL_TPM_DIR 指定；按列名选择两个分组。
tpm_dir <- Sys.getenv("CELL_TPM_DIR",file.path(data_root,"Cell_Bulk","gene_exp_matrix"))
tpm_names <- c(AML12="All_Group_TPM_AML12.csv", `TCMK-1`="All_Group_TPM_TCMK_1.csv", `293T`="All_Group_TPM_293T.csv")
for (model in names(tpm_names)) {
  path <- file.path(tpm_dir,tpm_names[[model]])
  if (!file.exists(path)) {message("TPM heatmap skipped; missing ",path); next}
  expr <- read.csv(path,check.names=FALSE)
  key <- intersect(c("gene","gene_id","symbol"),names(expr))
  genes <- expr[[if (length(key)) key[1] else 1L]]
  requested <- unique(unlist(key_gene_sets,use.names=FALSE))
  i <- match(toupper(requested),toupper(genes))
  hit <- !is.na(i)
  sample_names <- grep("^(control|lnp)-[0-9]+$",names(expr),value=TRUE,ignore.case=TRUE)
  groups <- ifelse(grepl("^lnp-",sample_names,ignore.case=TRUE),"LNP","Control")
  stopifnot(all(c("Control","LNP") %in% groups))
  if (sum(hit)<2) next
  mat <- as.matrix(expr[i[hit],sample_names,drop=FALSE])
  storage.mode(mat) <- "numeric"
  stopifnot(all(is.finite(mat)),all(mat>=0))
  means <- sapply(c("Control","LNP"),function(g) log2(1+rowMeans(mat[,groups==g,drop=FALSE])))
  rownames(means) <- requested[hit]
  write.csv(means,file.path(input_dir,model,"key_gene_log2_mean_TPM.csv"))
  pdf(file.path(input_dir,model,"key_gene_TPM_heatmap.pdf"),width=6,height=8)
  heatmap(means,scale="none",Colv=NA,labCol=plot_label(colnames(means)))
  dev.off()
}
