# 可用 THESIS_DATA 和 LNP_RESULTS 指定输入和输出目录。
data_root <- path.expand(Sys.getenv("THESIS_DATA", "~/Downloads/毕业论文数据"))
result_root <- path.expand(Sys.getenv("LNP_RESULTS", "LNP_results"))
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
# 仅在图形显示层使用 eLNP。
plot_label <- function(x) gsub("(?<![A-Za-z])LNP", "eLNP", as.character(x), perl = TRUE)

selected_groups <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("control", "ctrl", "ctr"), "Control",
         ifelse(x %in% c("lnp", "elnp"), "LNP", NA_character_))
}
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(1234)

configs <- list(
  Spleen = list(
    rds = Sys.getenv("SPLEEN_CELLCHAT_RDS", file.path(result_root, "Spleen_CellChat", "Spleen_cellchat_LNP_Control.rds")),
    out = file.path(result_root, "Spleen_CellChat_pathways"),
    groups = c("Control", "LNP"),
    cell_types = c("B cell", "Neutrophil", "T cell", "NK cell", "Macrophage", "Monocyte", "pDC"),
    labels = c("B", "Neut", "T", "NK", "Macro", "Mono", "pDC"),
    selected = c("ANNEXIN", "GALECTIN", "MIF", "COMPLEMENT", "MHC-II", "SELPLG"),
    representative = c("ANNEXIN", "MIF", "COMPLEMENT")
  ),
  Liver = list(
    rds = Sys.getenv("LIVER_CELLCHAT_RDS", file.path(result_root, "Liver_CellChat", "Liver_cellchat_LNP_Control.rds")),
    out = file.path(result_root, "Liver_CellChat_pathways"),
    groups = c("Control", "LNP"),
    cell_types = c("T cell", "B cell", "Endothelial cell", "Neutrophil", "Hepatocyte", "NK cell", "Kupffer cell", "HSC"),
    labels = c("T", "B", "Endo", "Neut", "Hep", "NK", "Kupffer", "HSC"),
    selected = c("ANNEXIN", "COLLAGEN", "SAA", "COMPLEMENT", "GALECTIN", "CXCL", "TNF"),
    representative = c("ANNEXIN", "COLLAGEN", "TNF")
  )
)

group_colours <- c(
  Control = "#8C8C8C",
  LNP = "#F28E2B"
)

category_levels <- c(
  "Immune activation and recruitment",
  "Inflammation and death signaling",
  "ECM remodeling and angiogenesis",
  "Metabolic stress and growth",
  "Other"
)

category_colours <- c(
  "Immune activation and recruitment" = "#6A51A3",
  "Inflammation and death signaling" = "#CB181D",
  "ECM remodeling and angiogenesis" = "#238B45",
  "Metabolic stress and growth" = "#F39C34",
  "Other" = "#C7C7C7"
)

classify_pathway <- function(pathway) {
  x <- stringr::str_to_upper(pathway)
  dplyr::case_when(
    stringr::str_detect(x, "CXCL|CCL|COMPLEMENT|ICAM|VCAM|SELPLG|SELL|MHC|CD80|CD86|CD40|CD45|CD22|CD48|BST2") ~ category_levels[[1]],
    stringr::str_detect(x, "TNF|IL[0-9]|IFN|MIF|CSF|GALECTIN|SPP1|FASLG|TRAIL|ANNEXIN|SAA|OSM|PAR") ~ category_levels[[2]],
    stringr::str_detect(x, "VEGF|ANGPT|PDGF|TGF|COLLAGEN|LAMININ|FN1|THBS|NOTCH|JAM|PECAM|CDH5|BMP|SEMA|VTN|HSPG|EPHB") ~ category_levels[[3]],
    stringr::str_detect(x, "HGF|FGF|EGF|IGF|WNT|VISFATIN|GAS|PROS|KIT") ~ category_levels[[4]],
    TRUE ~ category_levels[[5]]
  )
}

align_matrix <- function(mat, cell_types) {
  ans <- matrix(0, length(cell_types), length(cell_types), dimnames = list(cell_types, cell_types))
  rr <- intersect(rownames(mat), cell_types)
  cc <- intersect(colnames(mat), cell_types)
  ans[rr, cc] <- mat[rr, cc, drop = FALSE]
  ans
}

matrix_long <- function(mat, value_name = "value") {
  x <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(x) <- c("source", "target", value_name)
  x
}

pathway_matrix <- function(cellchat, pathway, cell_types) {
  pathways <- dimnames(cellchat@netP$prob)[[3]]
  if (!pathway %in% pathways) {
    return(matrix(0, length(cell_types), length(cell_types), dimnames = list(cell_types, cell_types)))
  }
  align_matrix(cellchat@netP$prob[, , pathway, drop = TRUE], cell_types)
}

circle_layout <- function(cell_types, labels) {
  theta <- seq(pi / 2, pi / 2 - 2 * pi, length.out = length(cell_types) + 1)[-length(cell_types) - 1]
  tibble(
    Cell_type = cell_types,
    label = labels,
    x = cos(theta),
    y = sin(theta)
  )
}

save_plot <- function(plot, stem, width, height, dpi = 320) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = "pdf", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
}

make_information_flow <- function(cellchat_list, cfg, tissue) {
  all_pathways <- sort(unique(unlist(lapply(cellchat_list, function(x) dimnames(x@netP$prob)[[3]]))))
  raw <- purrr::map_dfr(cfg$groups, function(group) {
    cc <- cellchat_list[[group]]
    paths <- dimnames(cc@netP$prob)[[3]]
    flow <- apply(cc@netP$prob, 3, sum, na.rm = TRUE)
    tibble(Group = group, Pathway = paths, information_flow = as.numeric(flow))
  }) |>
    tidyr::complete(Group = cfg$groups, Pathway = all_pathways, fill = list(information_flow = 0)) |>
    mutate(Group = factor(Group, levels = cfg$groups))

  selected <- intersect(cfg$selected, all_pathways)
  if (length(selected) < length(cfg$selected)) {
    warning(tissue, ": absent selected pathways: ", paste(setdiff(cfg$selected, selected), collapse = ", "))
  }
  heat <- raw |>
    filter(Pathway %in% selected) |>
    group_by(Pathway) |>
    mutate(relative_flow = if (max(information_flow) > 0) information_flow / max(information_flow) else 0) |>
    ungroup() |>
    mutate(
      Pathway = factor(Pathway, levels = rev(selected)),
      label = if_else(information_flow >= 0.01, formatC(information_flow, format = "f", digits = 2), "<0.01")
    )

  p <- ggplot(heat, aes(Group, Pathway, fill = relative_flow)) + scale_x_discrete(labels=plot_label) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = label), size = 3.2, fontface = "bold") +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1), name = "Relative flow") +
    labs(title = paste0(tissue, " selected signaling pathways"), subtitle = "CellChat information flow; colour normalized within each pathway", x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"), axis.text.x = element_text(face = "bold"), legend.position = "right")

  list(data = raw, selected = selected, plot = p)
}

make_pathway_network <- function(cellchat_list, cfg, tissue, pathway) {
  layout <- circle_layout(cfg$cell_types, cfg$labels)
  mats <- lapply(cfg$groups, function(g) pathway_matrix(cellchat_list[[g]], pathway, cfg$cell_types))
  names(mats) <- cfg$groups
  all_positive <- unlist(lapply(mats, function(x) x[x > 0]))
  cutoff <- if (length(all_positive) > 28) as.numeric(quantile(all_positive, 0.45, na.rm = TRUE)) else 0
  max_edge <- max(c(all_positive, 1e-8), na.rm = TRUE)
  hub_values <- unlist(lapply(mats, function(x) rowSums(x) + colSums(x)))
  max_hub <- max(c(hub_values, 1e-8), na.rm = TRUE)

  plots <- lapply(seq_along(cfg$groups), function(i) {
    group <- cfg$groups[[i]]
    show_legend <- i == 1
    mat <- mats[[group]]
    nodes <- tibble(
      Cell_type = cfg$cell_types,
      hub = rowSums(mat) + colSums(mat),
      present = cfg$cell_types %in% as.character(unique(cellchat_list[[group]]@idents))
    ) |>
      left_join(layout, by = "Cell_type") |>
      arrange(match(Cell_type, cfg$cell_types))
    edges <- matrix_long(mat, "prob") |>
      filter(prob > 0, prob >= cutoff) |>
      transmute(from = source, to = target, prob)
    graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE, node_key = "Cell_type")

    p <- ggraph(graph, layout = "manual", x = x, y = y) +
      geom_edge_link(aes(width = prob), colour = group_colours[[group]], alpha = 0.72,
                     arrow = grid::arrow(type = "closed", length = grid::unit(2.1, "mm")),
                     start_cap = circle(2.7, "mm"), end_cap = circle(3.7, "mm"), lineend = "round",
                     show.legend = show_legend) +
      geom_edge_loop(aes(width = prob), colour = group_colours[[group]], alpha = 0.58,
                     arrow = grid::arrow(type = "closed", length = grid::unit(1.8, "mm")),
                     start_cap = circle(2.7, "mm"), end_cap = circle(3.7, "mm"), show.legend = FALSE) +
      geom_node_point(aes(size = hub, alpha = present), shape = 21, fill = "#F7FBFF", colour = "#263238", stroke = 0.6,
                      show.legend = show_legend) +
      geom_node_text(aes(label = label), size = 3.2, fontface = "bold") +
      scale_edge_width_continuous(range = c(0.25, 3.2), limits = c(0, max_edge), name = "Communication probability") +
      scale_size_continuous(range = c(4.5, 12), limits = c(0, max_hub), name = "Pathway hub score") +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.15), guide = "none") +
      coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.18, 1.18), clip = "off") +
      labs(title = plot_label(group)) +
      theme_graph(base_family = "sans") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            legend.position = "bottom", legend.box = "vertical",
            legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))
    p
  })

  patchwork::wrap_plots(plots, nrow = 1, guides = "collect") +
    plot_annotation(
      title = paste0(tissue, " - ", pathway, " signaling network"),
      subtitle = "Common node and edge scales across groups; weakest edges are filtered for readability",
      theme = theme(plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 9))
    ) & theme(legend.position = "bottom")
}

make_lr_bubble <- function(cellchat_list, cfg, selected, tissue) {
  all_lr <- purrr::map_dfr(cfg$groups, function(group) {
    d <- CellChat::subsetCommunication(cellchat_list[[group]], signaling = selected)
    if (nrow(d) == 0) return(tibble())
    as_tibble(d) |>
      mutate(Group = group)
  }) |>
    mutate(
      Group = factor(Group, levels = cfg$groups),
      communication = paste0(interaction_name_2, " | ", source, " -> ", target),
      minus_log10_p = -log10(pmax(pval, 1e-6))
    )

  plot_data <- all_lr |>
    group_by(Group, pathway_name) |>
    slice_max(prob, n = 3, with_ties = FALSE) |>
    ungroup() |>
    mutate(communication = stringr::str_wrap(communication, 42))

  p <- ggplot(plot_data, aes(Group, communication)) + scale_x_discrete(labels=plot_label) +
    geom_point(aes(size = prob, colour = minus_log10_p), alpha = 0.85) +
    facet_wrap(~pathway_name, scales = "free_y", ncol = 2) +
    scale_colour_gradient(low = "#9ECAE1", high = "#CB181D", name = "-log10(p-value)") +
    scale_size_continuous(range = c(2, 8), name = "Probability") +
    labs(title = paste0(tissue, " representative ligand-receptor communications"), subtitle = "Top three interactions per group and pathway; full data are provided in CSV", x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 14), strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.major = element_line(linewidth = 0.2))

  list(data = all_lr, plot = p)
}

make_overall_summary <- function(cellchat_list, cfg, tissue) {
  layout <- circle_layout(cfg$cell_types, cfg$labels)
  rows <- list()
  nodes_all <- list()
  burden <- list()

  for (group in cfg$groups) {
    cc <- cellchat_list[[group]]
    weight <- align_matrix(cc@net$weight, cfg$cell_types)
    count <- align_matrix(cc@net$count, cfg$cell_types)
    nodes_all[[group]] <- tibble(
      Group = group, Cell_type = cfg$cell_types,
      hub = rowSums(weight) + colSums(weight),
      present = cfg$cell_types %in% as.character(unique(cc@idents))
    )
    burden[[group]] <- tibble(
      Group = group,
      `Total interaction number` = sum(count, na.rm = TRUE),
      `Total interaction strength` = sum(weight, na.rm = TRUE)
    )

    probs <- cc@netP$prob
    paths <- dimnames(probs)[[3]]
    cats <- classify_pathway(paths)
    category_mats <- lapply(category_levels, function(cat) {
      idx <- which(cats == cat)
      if (length(idx) == 0) return(matrix(0, nrow(weight), ncol(weight), dimnames = dimnames(weight)))
      apply(probs[, , idx, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
    })
    category_mats <- lapply(category_mats, align_matrix, cell_types = cfg$cell_types)
    arr <- simplify2array(category_mats)
    dominant_idx <- apply(arr, c(1, 2), which.max)
    dominant_strength <- apply(arr, c(1, 2), max)
    total_category <- apply(arr, c(1, 2), sum)
    dominant_prop <- ifelse(total_category > 0, dominant_strength / total_category, 0)
    dominant_cat <- matrix(category_levels[dominant_idx], nrow = length(cfg$cell_types), dimnames = dimnames(weight))

    rows[[group]] <- matrix_long(weight, "strength") |>
      left_join(matrix_long(dominant_prop, "dominant_proportion"), by = c("source", "target")) |>
      left_join(matrix_long(dominant_cat, "dominant_category"), by = c("source", "target")) |>
      mutate(Group = group)
  }

  edges <- bind_rows(rows)
  nodes <- bind_rows(nodes_all)
  burden_wide <- bind_rows(burden) |>
    mutate(Group = factor(Group, levels = cfg$groups))
  burden_long <- burden_wide |>
    pivot_longer(-Group, names_to = "Metric", values_to = "Value")
  positive <- edges |>
    filter(strength > 0, source != target) |>
    pull(strength)
  cutoff <- if (length(positive)) as.numeric(quantile(positive, 0.50, na.rm = TRUE)) else 0
  max_edge <- max(c(edges$strength, 1e-8), na.rm = TRUE)
  max_hub <- max(c(nodes$hub, 1e-8), na.rm = TRUE)

  burden_plot <- ggplot(burden_long, aes(Group, Value, fill = Group)) + scale_x_discrete(labels=plot_label) +
    geom_col(width = 0.68, colour = "#333333", linewidth = 0.25) +
    geom_text(aes(label = if_else(Metric == "Total interaction number", formatC(Value, format = "f", digits = 0), formatC(Value, format = "f", digits = 1))), vjust = -0.45, size = 3.1, fontface = "bold") +
    facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = group_colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(title = paste0(tissue, " communication burden"), x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

  network_plots <- lapply(seq_along(cfg$groups), function(i) {
    group <- cfg$groups[[i]]
    show_legend <- i == 1
    node_dat <- nodes |>
      filter(Group == group) |>
      left_join(layout, by = "Cell_type") |>
      arrange(match(Cell_type, cfg$cell_types))
    edge_dat <- edges |>
      filter(Group == group, source != target, strength > 0, strength >= cutoff) |>
      transmute(from = source, to = target, strength, dominant_category = factor(dominant_category, levels = category_levels), dominant_proportion)
    graph <- tidygraph::tbl_graph(nodes = node_dat, edges = edge_dat, directed = TRUE, node_key = "Cell_type")
    ggraph(graph, layout = "manual", x = x, y = y) +
      geom_edge_link(aes(width = strength, colour = dominant_category, alpha = dominant_proportion),
                     arrow = grid::arrow(type = "closed", length = grid::unit(2, "mm")),
                     start_cap = circle(2.8, "mm"), end_cap = circle(3.8, "mm"), lineend = "round",
                     show.legend = show_legend) +
      scale_edge_colour_manual(
        values = category_colours, drop = FALSE, name = "Dominant signaling category",
        labels = c("Immune recruitment", "Inflammation/death", "ECM/angiogenesis", "Metabolic/growth", "Other"),
        guide = guide_legend(nrow = 2, byrow = TRUE)
      ) +
      scale_edge_width_continuous(
        range = c(0.25, 3.3), limits = c(0, max_edge), name = "Interaction strength",
        guide = guide_legend(nrow = 1)
      ) +
      scale_edge_alpha_continuous(
        range = c(0.22, 0.95), limits = c(0, 1), name = "Dominant fraction",
        guide = guide_legend(nrow = 1)
      ) +
      geom_node_point(aes(size = hub, alpha = present), shape = 21, fill = "#E8F1F8", colour = "#263238", stroke = 0.6,
                      show.legend = show_legend) +
      geom_node_text(aes(label = label), size = 3.1, fontface = "bold") +
      scale_size_continuous(range = c(4.5, 12), limits = c(0, max_hub), name = "Hub score") +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.15), guide = "none") +
      coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.18, 1.18), clip = "off") +
      labs(title = plot_label(group)) +
      theme_graph(base_family = "sans") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            legend.position = "bottom", legend.box = "vertical",
            legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))
  })

  networks <- wrap_plots(network_plots, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))
  final <- burden_plot / networks +
    plot_layout(heights = c(0.9, 3.1), guides = "collect") +
    plot_annotation(
      title = paste0(tissue, " communication network summary"),
      subtitle = "Node size: hub score; edge width: total strength; edge colour: dominant functional category; strongest 50% of non-self edges shown",
      theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 9))
    ) & theme(legend.position = "bottom")

  list(plot = final, burden = burden_wide, nodes = nodes, edges = edges)
}

run_tissue <- function(cfg, tissue) {
  message("Reading ", tissue, " CellChat objects...")
  cellchat_list <- readRDS(cfg$rds)
  aliases <- selected_groups(names(cellchat_list))
  cellchat_list <- cellchat_list[!is.na(aliases)]
  names(cellchat_list) <- na.omit(aliases)
  stopifnot(!anyDuplicated(names(cellchat_list)))
  missing_groups <- setdiff(cfg$groups, names(cellchat_list))
  if (length(missing_groups)) stop(tissue, ": missing groups: ", paste(missing_groups, collapse = ", "))
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)

  info <- make_information_flow(cellchat_list, cfg, tissue)
  write.csv(info$data, file.path(cfg$out, "Pathway_information_flow.csv"), row.names = FALSE)
  selected_table <- tibble(
    Tissue = tissue,
    Pathway = info$selected,
    Category = classify_pathway(info$selected),
    Representative_panel = info$selected %in% cfg$representative
  )
  write.csv(selected_table, file.path(cfg$out, "Selected_pathways.csv"), row.names = FALSE)
  save_plot(info$plot, file.path(cfg$out, "Figure_A_pathway_information_flow_heatmap"), 8.4, max(4.5, 0.55 * length(info$selected) + 1.8))

  pathway_plots <- list()
  for (pathway in info$selected) {
    message(tissue, ": drawing ", pathway)
    p <- make_pathway_network(cellchat_list, cfg, tissue, pathway)
    pathway_plots[[pathway]] <- p
    save_plot(p, file.path(cfg$out, paste0("Figure_B_", gsub("[^A-Za-z0-9_-]", "_", pathway), "_LNP_Control_network")), 15.5, 5.2)
  }

  lr <- make_lr_bubble(cellchat_list, cfg, info$selected, tissue)
  write.csv(lr$data, file.path(cfg$out, "Ligand_receptor_bubble_data.csv"), row.names = FALSE)
  save_plot(lr$plot, file.path(cfg$out, "Figure_E_ligand_receptor_bubble"), 13, max(10, 2.3 * ceiling(length(info$selected) / 2)))

  overall <- make_overall_summary(cellchat_list, cfg, tissue)
  write.csv(overall$burden, file.path(cfg$out, paste0(tissue, "_overall_communication_burden.csv")), row.names = FALSE)
  write.csv(overall$nodes, file.path(cfg$out, paste0(tissue, "_overall_node_metrics.csv")), row.names = FALSE)
  write.csv(overall$edges, file.path(cfg$out, paste0(tissue, "_overall_edge_metrics.csv")), row.names = FALSE)
  save_plot(overall$plot, file.path(cfg$out, "Figure_F_overall_network_summary"), 16, 10.5)

  reps <- intersect(cfg$representative, names(pathway_plots))
  preview_parts <- c(list(info$plot), pathway_plots[reps], list(lr$plot, overall$plot))
  preview <- wrap_plots(preview_parts, ncol = 1, heights = c(1.2, rep(1.15, length(reps)), 2.2, 2.5)) +
    plot_annotation(
      title = paste0(tissue, " CellChat figure panel preview"),
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 18))
    )
  save_plot(preview, file.path(cfg$out, "Figure_complete_panel_preview"), 17, 30, dpi = 260)

  message(tissue, " outputs: ", normalizePath(cfg$out))
  invisible(list(info = info, pathway_plots = pathway_plots, lr = lr, overall = overall))
}

results <- purrr::imap(configs, run_tissue)
message("All CellChat pathway panels completed.")
