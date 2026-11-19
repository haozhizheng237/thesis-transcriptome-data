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
for (model in c("AML12","Liver")) {
  folder <- if (model == "AML12") "Cell_Bulk" else "Organ_Bulk"
  tab <- read.csv(file.path(data_root, folder, paste0("UA_rescue_",model,"_raw_counts.csv")), check.names = FALSE)
  # 批次按原实验记录对应到保留样本，允许输入列重新排序。
  labels <- if (model == "AML12") c("Control-1","Control-2","Control-3","LNP-1","LNP-2","LNP-3","LNP-UA-1","LNP-UA-2") else
    c("Control-1","Control-2","Control-3","LNP-1","LNP-2","LNP-UA-1","LNP-UA-2")
  batch <- if (model == "AML12") c(1,3,3,2,3,4,2,4) else c(1,2,3,3,4,3,4)
  group <- factor(gsub("-", "_", sub("-[0-9]+$", "", labels)), c("Control","LNP","LNP_UA"))
  meta <- data.frame(sample=paste(model, labels, sep="-"), group=group, batch=batch)
  out <- file.path(result_root, "UA_rescue", model)
  fit <- fit_counts(tab, meta, c(LNP_vs_Control="LNP-Control", LNP_UA_vs_LNP="LNP_UA-LNP"), out)
  main <- fit$tables$LNP_vs_Control
  ua <- fit$tables$LNP_UA_vs_LNP
  paired <- merge(main, ua, by="gene", suffixes=c("_LNP","_UA"))
  reversed <- subset(paired, significant_LNP != "NS" & significant_UA != "NS" & logFC_LNP * logFC_UA < 0)
  write.csv(reversed, file.path(out,"opposite_direction_DEGs.csv"), row.names=FALSE)
  if (nrow(reversed) >= 2) {
    mat <- fit$expression[reversed$gene,,drop=FALSE]
    means <- sapply(levels(group), function(g) rowMeans(mat[, group == g, drop=FALSE]))
    pdf(file.path(out,"UA_rescue_heatmap.pdf"))
    heatmap(means, labCol=plot_label(colnames(means)), scale="row", Colv=NA, labRow=NA)
    dev.off()
  }
  for (direction in c("Up","Down")) {
    genes <- reversed$gene[reversed$significant_LNP == direction]
    if (length(genes) >= 3) {
      go <- clusterProfiler::enrichGO(genes, OrgDb=org.Mm.eg.db::org.Mm.eg.db, keyType="SYMBOL",
                 ont="BP", universe=main$gene, pvalueCutoff=1, qvalueCutoff=1)
      write.csv(as.data.frame(go),file.path(out,paste0("rescue_GO_",direction,".csv")),row.names=FALSE)
    }
  }
}
# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 1. 数据
data <- data.frame(
  Control = c(0.905474, 0.996219, 0.936588, 1.113227, 0.974182, 0.956574924, 0.968807339, 1.040061162, 0.934556575),
  LNP = c(0.819266055, 0.780122324, 0.762996942, 0.772782875, NA, NA, NA, NA, NA),
  `LNP+1 nmol Ursolic acid` = c(0.882568807, 0.817431193, 0.788073394, NA, NA, NA, NA, NA, NA),
  `LNP+10 nmol Ursolic acid` = c(0.863914373, 0.83883792, 0.785321101, NA, NA, NA, NA, NA, NA),
  `LNP+50 nmol Ursolic acid` = c(0.946411, 0.881109, 0.908556, NA, NA, NA, NA, NA, NA),
  `LNP+100 nmol Ursolic acid` = c(0.817677285, 0.776744, 0.7815619, NA, NA, NA, NA, NA, NA),
  check.names = FALSE  # 保持原始列名
)

# 2. 转换为长格式（明确使用tidyr::pivot_longer）
data_long <- data %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Treatment",
    values_to = "Value"
  ) %>%
  tidyr::drop_na()

out <- file.path(result_root,"UA_rescue")
reference <- na.omit(data$LNP)
doses <- grep("^LNP[+]", names(data), value=TRUE)
dose_stats <- data.frame(dose=doses, p_value=vapply(doses, function(d) t.test(na.omit(data[[d]]),reference)$p.value, numeric(1)))
dose_stats$FDR <- p.adjust(dose_stats$p_value, "BH")
write.csv(dose_stats,file.path(out,"UA_dose_comparisons.csv"),row.names=FALSE)
p <- ggplot(data_long,aes(factor(Treatment,levels=names(data)),Value)) + geom_boxplot() + geom_jitter(width=.08) + theme_classic() + scale_x_discrete(labels=plot_label) + labs(x=NULL)
ggsave(file.path(out,"UA_dose_response.pdf"),p,width=10,height=5)
data <- data.frame(
  LNP = c(5.212, 6.156, 5.481, 5.117, 6.731),
  LNP_UA = c(6.374, 5.672, 5.723, 5.856, 5.063)
)
delivery <- t.test(data$LNP_UA, data$LNP, paired=TRUE)
write.csv(data,file.path(out,"delivery_paired_measurements.csv"),row.names=FALSE)
capture.output(delivery,file=file.path(out,"delivery_paired_t_test.txt"))
pdf(file.path(out,"delivery_efficiency.pdf"))
matplot(1:2, t(as.matrix(data)),type="b",xaxt="n",xlab="",ylab="Delivery efficiency")
axis(1,1:2,plot_label(c("LNP","LNP+UA")))
dev.off()
