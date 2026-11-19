# 可用 THESIS_DATA 和 LNP_RESULTS 指定输入和输出目录。
data_root <- path.expand(Sys.getenv("THESIS_DATA", "~/Downloads/毕业论文数据"))
result_root <- path.expand(Sys.getenv("LNP_RESULTS", "LNP_results"))
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
# 仅在图形显示层使用 eLNP。
plot_label <- function(x) gsub("(?<![A-Za-z])LNP", "eLNP", as.character(x), perl = TRUE)

########## scRNA data cell-cell communication whole network ########## 
#
# Spleen CellChat v2: LNP and control communication network summary
# Required Seurat object: Spleen_final
# Metadata columns: Group, Cell_type, orig.ident
# 

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
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
future::plan(strategy = "multisession", workers = 4)

# 
# 1. Parameters
# 

Group_levels <- c("Control", "LNP")

Cell_type_levels <- c(
  "B cell",
  "Neutrophil",
  "T cell",
  "NK cell",
  "Macrophage",
  "Monocyte",
  "pDC"
)

min_cell <- 3
population_size <- FALSE
mean_method <- "triMean"
trim <- 0.1
nboot <- 100
edge_quantile <- 0

output_dir <- file.path(result_root, "Spleen_CellChat")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

selected_groups <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("control", "ctrl", "ctr"), "Control",
         ifelse(x %in% c("lnp", "elnp"), "LNP", NA_character_))
}
Spleen_final <- readRDS(file.path(data_root, "scRNA_data", "Spleen_scRNA_data.rds"))
groups <- selected_groups(Spleen_final$Group)
Spleen_final <- Spleen_final[, !is.na(groups)]
Spleen_final$Group <- groups[!is.na(groups)]
stopifnot(all(Group_levels %in% Spleen_final$Group))
stopifnot(inherits(Spleen_final, "Seurat"))
stopifnot(all(c("Group", "Cell_type", "orig.ident") %in% colnames(Spleen_final[[]])))

Spleen_final$Group <- factor(
  as.character(Spleen_final$Group),
  levels = Group_levels
)

Spleen_final$Cell_type <- factor(
  trimws(as.character(Spleen_final$Cell_type)),
  levels = Cell_type_levels
)

if (anyNA(Spleen_final$Group)) {
  stop("Group contains values not listed in Group_levels.")
}

if (anyNA(Spleen_final$Cell_type)) {
  stop("Cell_type contains values not listed in Cell_type_levels.")
}

print(table(Spleen_final$Group, Spleen_final$Cell_type))

# 
# 2. CellChat database check
# 

data("CellChatDB.mouse", package = "CellChat")
CellChatDB.use <- CellChatDB.mouse

DB_genes <- CellChat::extractGene(CellChatDB.use)
overlap_genes <- intersect(rownames(Spleen_final), DB_genes)

message(
  "CellChatDB interactions = ", nrow(CellChatDB.use$interaction),
  "; DB genes = ", length(DB_genes),
  "; matched genes = ", length(overlap_genes)
)

if (length(overlap_genes) < 3) {
  stop("Fewer than three genes overlap with CellChatDB.mouse.")
}

# 
# 3. Run CellChat by group
# 

run_cellchat_group <- function(
    seu,
    group_name,
    group_col = "Group",
    cell_type_col = "Cell_type",
    assay = "RNA",
    species = "mouse",
    min_cells = 3,
    population_size = FALSE,
    mean_method = c("triMean", "truncatedMean"),
    trim = 0.1,
    nboot = 100,
    seed = 1234) {
  
  species <- match.arg(species)
  mean_method <- match.arg(mean_method)
  set.seed(seed)
  
  meta <- seu[[]]
  cells_use <- rownames(meta)[
    as.character(meta[[group_col]]) == group_name &
      !is.na(meta[[cell_type_col]])
  ]
  
  if (length(cells_use) == 0) {
    stop("No cells found for group: ", group_name)
  }
  
  seu_group <- subset(seu, cells = cells_use)
  cell_number <- table(seu_group[[]][[cell_type_col]])
  keep_types <- names(cell_number)[cell_number >= min_cells]
  
  if (length(keep_types) < 2) {
    stop(
      group_name,
      " has fewer than two cell types with at least ",
      min_cells,
      " cells."
    )
  }
  
  keep_cells <- rownames(seu_group[[]])[
    as.character(seu_group[[]][[cell_type_col]]) %in% keep_types
  ]
  
  seu_group <- subset(seu_group, cells = keep_cells)
  seu_group[[cell_type_col]] <- factor(
    as.character(seu_group[[]][[cell_type_col]]),
    levels = keep_types
  )
  
  message(
    "Running ", group_name,
    "; cells = ", ncol(seu_group),
    "; cell types = ", paste(keep_types, collapse = ", ")
  )
  
  cellchat <- CellChat::createCellChat(
    object = seu_group,
    group.by = cell_type_col,
    assay = assay
  )
  
  if (species == "mouse") {
    data("CellChatDB.mouse", package = "CellChat", envir = environment())
    database_use <- CellChatDB.mouse
  } else {
    data("CellChatDB.human", package = "CellChat", envir = environment())
    database_use <- CellChatDB.human
  }
  
  # Do not run subsetDB() without explicit arguments.
  cellchat@DB <- database_use
  cellchat <- CellChat::subsetData(cellchat)
  
  if (nrow(cellchat@data.signaling) < 3) {
    stop(group_name, ": fewer than three signaling genes were retained.")
  }
  
  message(
    group_name,
    ": data.signaling = ", nrow(cellchat@data.signaling),
    " genes x ", ncol(cellchat@data.signaling), " cells"
  )
  
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  
  if (mean_method == "triMean") {
    cellchat <- CellChat::computeCommunProb(
      cellchat,
      type = "triMean",
      population.size = population_size,
      nboot = nboot
    )
  } else {
    cellchat <- CellChat::computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = trim,
      population.size = population_size,
      nboot = nboot
    )
  }
  
  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = min_cells
  )
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat
}

cellchat_list <- setNames(
  lapply(Group_levels, function(group_name) {
    run_cellchat_group(
      seu = Spleen_final,
      group_name = group_name,
      group_col = "Group",
      cell_type_col = "Cell_type",
      assay = "RNA",
      species = "mouse",
      min_cells = min_cell,
      population_size = population_size,
      mean_method = mean_method,
      trim = trim,
      nboot = nboot,
      seed = 1234
    )
  }),
  Group_levels
)

saveRDS(
  cellchat_list,
  file = file.path(output_dir, "Spleen_cellchat_LNP_Control.rds")
)

liver <- readRDS(file.path(data_root,"scRNA_data","Liver_scRNA_data.rds"))
groups <- selected_groups(liver$Group)
liver <- liver[,!is.na(groups)]
liver$Group <- groups[!is.na(groups)]
liver_chats <- setNames(lapply(Group_levels,function(g)
  run_cellchat_group(seu=liver,group_name=g,group_col="Group",cell_type_col="Cell_type",
                     min_cells=min_cell,population_size=population_size,mean_method=mean_method,
                     trim=trim,nboot=nboot,seed=1234)),Group_levels)
liver_dir <- file.path(result_root,"Liver_CellChat")
dir.create(liver_dir,recursive=TRUE,showWarnings=FALSE)
saveRDS(liver_chats,file.path(liver_dir,"Liver_cellchat_LNP_Control.rds"))
rm(liver,liver_chats)

# 
# 4. Signaling categories
# 

category_levels <- c(
  "Immune activation and recruitment",
  "Inflammation and death signaling",
  "ECM remodeling and angiogenesis",
  "Metabolic stress and senescence",
  "Other"
)

category_colours <- c(
  "Immune activation and recruitment" = "#6A51A3",    # 深紫（免疫）
  "Inflammation and death signaling" = "#CB181D",     # 亮红（炎症死亡）
  "ECM remodeling and angiogenesis" = "#238B45",      # 深绿（修复/重塑）
  "Metabolic stress and senescence" = "#F39C34",      # 橙黄（应激）
  "Other" = "#D3D3D3"                               # 浅灰（忽略）
)

classify_signaling_pathway <- function(pathway) {
  pathway_upper <- stringr::str_to_upper(pathway)
  
  dplyr::case_when(
    stringr::str_detect(
      pathway_upper,
      "CXCL|CCL|COMPLEMENT|ICAM|VCAM|SELECTIN|SELPLG|SELL|MHC|CD80|CD86|TCR|BCR|CTLA4|PD-L1|PD-1|CD40"
    ) ~ "Immune activation and recruitment",
    stringr::str_detect(
      pathway_upper,
      "TNF|IL[0-9]|INTERFERON|IFN|MIF|CSF|GALECTIN|SPP1|FASLG|TRAIL|TNFRSF|CASPASE|BAX|BCL|RIPK|MLKL"
    ) ~ "Inflammation and death signaling",
    stringr::str_detect(
      pathway_upper,
      paste(
        "VEGF|VEGFR|ANGPT|PDGF|TGF|COLLAGEN|LAMININ|FN1|THBS|
        PERIOSTIN|NOTCH|JAM|PECAM|CDH5|BMP|SEMA|TIE|S1P|EDN1|CTGF|PLGF|
        INTEGRIN|ITG|CD44|SEMA3",
        sep = "|"
      )
    ) ~ "ECM remodeling and angiogenesis",
    stringr::str_detect(
      pathway_upper,
      paste(
        "HGF|FGF|EGF|IGF|INSULIN|GLUCAGON|ADIPONECTIN|VISFATIN|RESISTIN|WNT|PTN|MK|GAS|PROS|APELIN|P16|P21|P53|CDKN|TERT",
        sep = "|"
      )
    ) ~ "Metabolic stress and senescence",
    TRUE ~ "Other"
  )
}

pathway_annotation <- purrr::imap_dfr(
  cellchat_list,
  function(cellchat, group_name) {
    tibble::tibble(
      Group = group_name,
      pathway = cellchat@netP$pathways
    )
  }
) |>
  dplyr::distinct(pathway) |>
  dplyr::mutate(
    category = classify_signaling_pathway(pathway),
    category = factor(category, levels = category_levels)
  ) |>
  dplyr::arrange(category, pathway)

utils::write.csv(
  pathway_annotation,
  file.path(output_dir, "CellChat_pathway_category_annotation.csv"),
  row.names = FALSE
)

# 
# 5. Extract network summaries
# 

align_cellchat_matrix <- function(mat, cell_types) {
  output <- matrix(
    0,
    nrow = length(cell_types),
    ncol = length(cell_types),
    dimnames = list(cell_types, cell_types)
  )
  
  common_rows <- intersect(rownames(mat), cell_types)
  common_cols <- intersect(colnames(mat), cell_types)
  output[common_rows, common_cols] <- mat[
    common_rows,
    common_cols,
    drop = FALSE
  ]
  output
}

matrix_to_long <- function(mat, value_name) {
  output <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(output) <- c("source", "target", value_name)
  output |>
    dplyr::mutate(
      source = as.character(source),
      target = as.character(target)
    )
}

extract_cellchat_summary <- function(
    cellchat,
    group_name,
    all_cell_types,
    category_levels) {
  
  count_mat <- align_cellchat_matrix(cellchat@net$count, all_cell_types)
  weight_mat <- align_cellchat_matrix(cellchat@net$weight, all_cell_types)
  
  node_table <- tibble::tibble(
    Cell_type = all_cell_types,
    outgoing = rowSums(weight_mat),
    incoming = colSums(weight_mat)
  ) |>
    dplyr::mutate(
      hub_score = outgoing + incoming,
      Group = group_name,
      present = Cell_type %in% as.character(unique(cellchat@idents))
    )
  
  edge_table <- matrix_to_long(weight_mat, "interaction_strength") |>
    dplyr::left_join(
      matrix_to_long(count_mat, "interaction_number"),
      by = c("source", "target")
    )
  
  pathway_prob <- cellchat@netP$prob
  if (length(dim(pathway_prob)) != 3) {
    stop(group_name, ": netP$prob is not a three-dimensional array.")
  }
  
  pathway_names <- dimnames(pathway_prob)[[3]]
  pathway_key <- tibble::tibble(
    pathway = pathway_names,
    category = classify_signaling_pathway(pathway_names)
  )
  
  # ---------- 按功能类别汇总通路通讯 ----------
  category_edge_table <- purrr::map_dfr(
    category_levels,
    function(category_name) {
      pathway_index <- which(pathway_key$category == category_name)
      
      if (length(pathway_index) == 0) {
        category_mat <- matrix(
          0,
          nrow = dim(pathway_prob)[1],
          ncol = dim(pathway_prob)[2],
          dimnames = dimnames(pathway_prob)[1:2]
        )
      } else {
        category_mat <- apply(
          pathway_prob[, , pathway_index, drop = FALSE],
          MARGIN = c(1, 2),
          FUN = sum,
          na.rm = TRUE
        )
      }
      
      category_mat <- align_cellchat_matrix(category_mat, all_cell_types)
      matrix_to_long(category_mat, "category_strength") |>
        dplyr::mutate(category = category_name)
    }
  )
  
  # ---------- 计算每个细胞对的总通讯强度 ----------
  category_edge_table <- category_edge_table |>
    dplyr::group_by(source, target) |>
    dplyr::mutate(
      total_strength = sum(category_strength, na.rm = TRUE)
    ) |>
    dplyr::ungroup()
  
  # ---------- 确定优势类别及其可信度 ----------
  dominant_category <- category_edge_table |>
    dplyr::group_by(source, target) |>
    dplyr::arrange(
      dplyr::desc(category_strength),
      factor(category, levels = category_levels),  # 极少数平局时的兜底
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      # 优势类别占总通讯的比例
      dominant_proportion = dplyr::if_else(
        total_strength > 0,
        category_strength / total_strength,
        0
      ),
      # 添加可信度分层（便于后续过滤）
      dominance_confidence = dplyr::case_when(
        total_strength == 0 ~ "No communication",
        dominant_proportion >= 0.7 ~ "High (dominant)",
        dominant_proportion >= 0.5 ~ "Medium (major)",
        TRUE ~ "Low (mixed)"
      )
    ) |>
    dplyr::select(
      source, target,
      dominant_category = category,
      dominant_category_strength = category_strength,
      total_strength,
      dominant_proportion,
      dominance_confidence
    )
  
  # ---------- 合并到边表 ----------
  edge_table <- edge_table |>
    dplyr::left_join(dominant_category, by = c("source", "target")) |>
    dplyr::mutate(
      # 无通讯的边标记为 "Other"（不再用"Basal homeostatic"）
      dominant_category = dplyr::if_else(
        is.na(dominant_category) | is.na(total_strength) | total_strength == 0,
        "Other",
        dominant_category
      ),
      dominant_category = factor(dominant_category, levels = category_levels),
      Group = group_name
    )
  
  list(
    nodes = node_table,
    edges = edge_table,
    category_edges = category_edge_table |>
      dplyr::mutate(Group = group_name),
    count_matrix = count_mat,
    weight_matrix = weight_mat
  )
}

network_summary <- purrr::imap(
  cellchat_list,
  function(cellchat, group_name) {
    extract_cellchat_summary(
      cellchat = cellchat,
      group_name = group_name,
      all_cell_types = Cell_type_levels,
      category_levels = category_levels
    )
  }
)

nodes_all <- purrr::map_dfr(network_summary, "nodes") |>
  dplyr::mutate(Group = factor(Group, levels = Group_levels))

edges_all <- purrr::map_dfr(network_summary, "edges") |>
  dplyr::mutate(Group = factor(Group, levels = Group_levels))

category_edges_all <- purrr::map_dfr(network_summary, "category_edges") |>
  dplyr::mutate(Group = factor(Group, levels = Group_levels))

utils::write.csv(
  nodes_all,
  file.path(output_dir, "Spleen_network_node_metrics.csv"),
  row.names = FALSE
)

utils::write.csv(
  edges_all,
  file.path(output_dir, "Spleen_network_edge_metrics.csv"),
  row.names = FALSE
)

# 
# 6. Communication burden
# 

burden_table <- edges_all |>
  dplyr::group_by(Group) |>
  dplyr::summarise(
    total_interaction_number = sum(interaction_number, na.rm = TRUE),
    total_interaction_strength = sum(interaction_strength, na.rm = TRUE),
    .groups = "drop"
  )

burden_long <- burden_table |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("total_"),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      total_interaction_number = "Total interaction number",
      total_interaction_strength = "Total interaction strength"
    )
  )

Group_colours <- c(
  "Control" = "grey70",
  "LNP" = "#FF7F0E"
)

y_limits <- burden_long |>
  dplyr::group_by(metric) |>
  dplyr::summarise(
    y_max = max(value, na.rm = TRUE) * 1.25  # 顶部留出25%空间
  )

# 创建辅助列用于设置y轴上限
burden_long <- burden_long |>
  dplyr::left_join(y_limits, by = "metric")

burden_plot <- ggplot2::ggplot(
  burden_long,
  ggplot2::aes(x = Group, y = value, fill = Group)
) +
  ggplot2::geom_col(width = 0.68, colour = "black", linewidth = 0.3) +
  
  # ====== 分面时使用自定义的y轴上限 
  ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
  
  # ====== 使用geom_blank()强制每个分面使用不同的y轴上限 
  ggplot2::geom_blank(
    data = burden_long,
    ggplot2::aes(y = y_max)
  ) +
  
  # ====== 或者使用更简洁的方法：直接设置expand参数 
  # 注意：如果上面用了geom_blank，下面这行可以省略
  # 如果不用geom_blank，用这行：
  # ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.3))) +
  
  ggplot2::scale_x_discrete(labels=plot_label) + ggplot2::scale_fill_manual(values = Group_colours,labels=plot_label) +
  
  # ====== 数值标签：使用vjust和nudge_y调整位置
  ggplot2::geom_text(
    ggplot2::aes(
      label = round(value, 1),
      y = value  # 默认就在柱子顶端
    ),
    vjust = -0.8,        # 向上偏移（负值=向上）
    size = 3.5,
    fontface = "bold"
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = NULL,
    title = "Spleen communication burden"
  ) +
  
  ggplot2::theme_classic(base_size = 11) +
  
  ggplot2::theme(
    legend.position = "none",
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    # ====== 增加分面板与柱子的间距 
    strip.placement = "outside"
  )


# 
# 7. Network plots
# 

# 节点布局（固定坐标，保证跨组可比）
node_layout <- tibble::tribble(
  ~Cell_type,          ~x,    ~y,    ~label,
  "Macrophage",        0.00,  0.00,  "Macro",
  "Monocyte",        -0.62,  0.38,  "Mono",
  "pDC",              0.62,  0.38,  "pDC",
  "Neutrophil",       -1.00,  0.92,  "Neut",
  "NK cell",          -0.92, -0.48,  "NK",
  "T cell",            0.92, -0.48,  "T",
  "B cell",            1.00,  0.92,  "B"
)

# 计算边强度阈值（只显示前40%强边）
positive_edges <- edges_all |>
  dplyr::filter(interaction_strength > 0, source != target) |>
  dplyr::pull(interaction_strength)

edge_cutoff <- if (length(positive_edges) > 0) {
  as.numeric(stats::quantile(
    positive_edges,
    probs = edge_quantile,  # 0.60
    na.rm = TRUE
  ))
} else {
  0
}

# 全局最大尺度（保证所有子图使用统一标尺）
max_hub <- max(nodes_all$hub_score, na.rm = TRUE)
max_edge <- max(edges_all$interaction_strength, na.rm = TRUE)
if (!is.finite(max_hub) || max_hub <= 0) max_hub <- 1
if (!is.finite(max_edge) || max_edge <= 0) max_edge <- 1

# 核心绘图函数
plot_Spleen_network <- function(
    group_name,
    nodes_all,
    edges_all,
    node_layout,
    edge_cutoff,
    category_colours,
    max_hub,
    max_edge) {
  
  # 准备节点数据
  node_dat <- nodes_all |>
    dplyr::filter(Group == group_name) |>
    dplyr::left_join(node_layout, by = "Cell_type") |>
    dplyr::arrange(match(Cell_type, node_layout$Cell_type))
  
  # 准备边数据（过滤：非自环、强度≥阈值）
  edge_dat <- edges_all |>
    dplyr::filter(
      Group == group_name,
      source != target,
      interaction_strength >= edge_cutoff,
      interaction_strength > 0
    )
  
  # 构建图对象
  graph_object <- tidygraph::tbl_graph(
    nodes = node_dat,
    edges = edge_dat |>
      dplyr::transmute(
        from = source,
        to = target,
        interaction_strength,
        dominant_category,
        dominance_confidence  # 新增：可信度信息
      ),
    directed = TRUE,
    node_key = "Cell_type"
  )
  
  # 组名显示（SiO₂带下标2）
  display_name <- plot_label(group_name)
  
  # 绘图
  p <- ggraph::ggraph(graph_object, layout = "manual", x = x, y = y) +
    
    # ---- 边图层 ----
  ggraph::geom_edge_link(
    ggplot2::aes(
      width = interaction_strength,
      colour = dominant_category,
      edge_alpha = dominance_confidence
    ),
    lineend = "round",
    arrow = grid::arrow(
      type = "closed",
      length = grid::unit(2.5, "mm")
    ),
    start_cap = ggraph::circle(2.5, "mm"),
    end_cap = ggraph::circle(3.5, "mm")
  ) +
    
    # ---- 边透明度（按可信度分层） ----
  ggraph::scale_edge_alpha_manual(
    values = c(
      "High (dominant)" = 1.00,
      "Medium (major)" = 0.5,
      "Low (mixed)" = 0.25,
      "No communication" = 0
    ),
    limits = c(
      "High (dominant)",
      "Medium (major)",
      "Low (mixed)"
    ),
    breaks = c(
      "High (dominant)",
      "Medium (major)",
      "Low (mixed)"
    ),
    drop = FALSE,
    name = "Dominance confidence"
  )  +
    
    # ---- 节点图层 ----
  ggraph::geom_node_point(
    ggplot2::aes(size = hub_score, alpha = present),
    shape = 21,
    fill = "#E8F1F8",
    colour = "#263238",
    stroke = 0.6
  ) +
    
    # ---- 节点标签 ----
  ggraph::geom_node_text(
    ggplot2::aes(label = label),
    size = 3.4,
    fontface = "bold",
    repel = FALSE
  ) +
    
    # ---- 在节点下方显示hub_score数值（方便精确读取） ----
  ggraph::geom_node_text(
    ggplot2::aes(
      label = ifelse(present, round(hub_score, 1), ""),
      x = x,
      y = y - 0.05
    ),
    size = 3,
    colour = "#666666",
    fontface = "plain"
  ) +
    
    # ---- 图例与标度 ----
  ggraph::scale_edge_colour_manual(
    values = category_colours,
    drop = FALSE,
    name = "Dominant signaling"
  ) +
    
    ggraph::scale_edge_width_continuous(
      range = c(0.3, 4.0),
      limits = c(0, max_edge),
      name = "Interaction strength"
    ) +
    
    ggplot2::scale_size_continuous(
      range = c(3.5, 13),
      limits = c(0, max_hub),
      name = "Hub score"
    ) +
    
    ggplot2::scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.15),
      guide = "none"
    ) +
    
    # ---- 坐标系与主题 ----
  ggplot2::coord_equal(
    xlim = c(-1.25, 1.25),
    ylim = c(-1.05, 1.15),
    clip = "off"
  ) +
    
    ggplot2::labs(title = display_name) +
    
    ggraph::theme_graph(base_family = "sans") +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = ggplot2::margin(t = -5, b = 2)
    )
  
  return(p)
}

# 生成两组网络图
network_plots <- lapply(Group_levels, function(group_name) {
  plot_Spleen_network(
    group_name = group_name,
    nodes_all = nodes_all,
    edges_all = edges_all,
    node_layout = node_layout,
    edge_cutoff = edge_cutoff,
    category_colours = category_colours,
    max_hub = max_hub,
    max_edge = max_edge
  )
})

# 使用patchwork横向排列（共享图例）
network_row <- patchwork::wrap_plots(
  network_plots,
  nrow = 1,
  guides = "collect"
) &
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

# 组合最终面板：上方为通讯负荷柱状图，下方为两个网络图
final_panel <- (burden_plot / network_row) +
  patchwork::plot_layout(
    heights = c(0.9, 3.2),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "Spleen communication network summary",
    subtitle = paste0(
      "Node size: total hub score (value shown below node); ",
      "edge width: communication strength; ",
      "edge colour: dominant signaling category; ",
      "edge transparency: dominance confidence"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0)
    )
  )

ggplot2::ggsave(
  file.path(output_dir, "Figure4_Spleen_network_summary.pdf"),
  plot = final_panel,
  width = 16,
  height = 9,
  units = "in",
  device = "pdf",
  limitsize = FALSE
)

ggplot2::ggsave(
  file.path(output_dir, "Figure4_Spleen_network_summary.png"),
  plot = final_panel,
  width = 16,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

# 
# 8. Summary tables
# 

main_hub_table <- nodes_all |>
  dplyr::group_by(Group) |>
  dplyr::slice_max(order_by = hub_score, n = 3, with_ties = FALSE) |>
  dplyr::arrange(Group, dplyr::desc(hub_score), .by_group = TRUE) |>
  dplyr::mutate(hub_rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::select(
    Group,
    hub_rank,
    Cell_type,
    outgoing,
    incoming,
    hub_score
  )

dominant_program_table <- category_edges_all |>
  dplyr::group_by(Group, category) |>
  dplyr::summarise(
    category_information_flow = sum(category_strength, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(Group) |>
  dplyr::mutate(
    total_information_flow = sum(category_information_flow),
    fraction = dplyr::if_else(
      total_information_flow > 0,
      category_information_flow / total_information_flow,
      0
    )
  ) |>
  dplyr::arrange(Group, dplyr::desc(fraction)) |>
  dplyr::ungroup()

main_program_table <- dominant_program_table |>
  dplyr::group_by(Group) |>
  dplyr::slice_max(
    order_by = category_information_flow,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

Group_summary_table <- burden_table |>
  dplyr::left_join(
    main_hub_table |>
      dplyr::filter(hub_rank == 1) |>
      dplyr::select(
        Group,
        main_hub = Cell_type,
        main_hub_score = hub_score
      ),
    by = "Group"
  ) |>
  dplyr::left_join(
    main_program_table |>
      dplyr::select(
        Group,
        dominant_program = category,
        dominant_program_fraction = fraction
      ),
    by = "Group"
  )

utils::write.csv(
  burden_table,
  file.path(output_dir, "Figure4_Spleen_communication_burden.csv"),
  row.names = FALSE
)

utils::write.csv(
  main_hub_table,
  file.path(output_dir, "Figure4_Spleen_main_hubs.csv"),
  row.names = FALSE
)

utils::write.csv(
  dominant_program_table,
  file.path(output_dir, "Figure4_Spleen_dominant_programs.csv"),
  row.names = FALSE
)

utils::write.csv(
  Group_summary_table,
  file.path(output_dir, "Figure4_Spleen_network_group_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    cellchat_list = cellchat_list,
    network_summary = network_summary,
    nodes_all = nodes_all,
    edges_all = edges_all,
    category_edges_all = category_edges_all,
    burden_table = burden_table,
    main_hub_table = main_hub_table,
    dominant_program_table = dominant_program_table,
    Group_summary_table = Group_summary_table,
    final_panel = final_panel
  ),
  file = file.path(output_dir, "Spleen_CellChat_complete_results.rds")
)

message("Analysis completed. Results saved in: ", normalizePath(output_dir))
