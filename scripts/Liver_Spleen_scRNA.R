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
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
set.seed(1234)
for (organ in c("Liver","Spleen")) {
  obj <- readRDS(file.path(data_root,"scRNA_data",paste0(organ,"_scRNA_data.rds")))
  groups <- selected_groups(obj$Group)
  obj <- obj[, !is.na(groups)]
  obj$Group <- factor(groups[!is.na(groups)],c("Control","LNP"))
  stopifnot(all(table(obj$Group)>0), all(c("Cell_type","orig.ident") %in% names(obj[[]])))
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]],"Assay5")) obj <- SeuratObject::JoinLayers(obj,assay="RNA")
  obj <- NormalizeData(obj, verbose=FALSE)
  out <- file.path(result_root,"scRNA",organ)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  meta <- obj[[]]
  cells <- as.data.frame(table(meta$orig.ident,meta$Cell_type),stringsAsFactors=FALSE)
  names(cells) <- c("sample","cell_type","cells")
  cells <- merge(cells, unique(data.frame(sample=meta$orig.ident,group=meta$Group)),by="sample")
  cells$fraction <- cells$cells / ave(cells$cells,cells$sample,FUN=sum)
  write.csv(cells,file.path(out,"sample_cell_composition.csv"),row.names=FALSE)
  p <- ggplot(cells,aes(sample,fraction,fill=cell_type)) + geom_col() + facet_grid(~group,labeller=labeller(group=plot_label),scales="free_x",space="free_x") + theme_bw()
  ggsave(file.path(out,"sample_cell_composition.pdf"),p)
  obj$figure_group <- plot_label(obj$Group)
  if ("umap" %in% Reductions(obj)) ggsave(file.path(out,"annotated_UMAP.pdf"),
      DimPlot(obj,reduction="umap",group.by="Cell_type",split.by="figure_group"),width=12,height=6)
  for (ct in unique(na.omit(as.character(meta$Cell_type)))) {
    sub <- obj[,rownames(meta)[which(meta$Cell_type == ct)]]
    if (any(table(factor(sub$Group,levels=c("Control","LNP"))) < 3)) next
    tab <- FindMarkers(sub,ident.1="LNP",ident.2="Control",group.by="Group",min.pct=.25)
    tab$gene <- rownames(tab)
    stem <- file.path(out,paste0(gsub("[^[:alnum:]_-]","_",ct),"_LNP_vs_Control"))
    tab$regulation <- ifelse(tab$p_val < .05 & abs(tab$avg_log2FC)>1,ifelse(tab$avg_log2FC>0,"Up","Down"),"NS")
    write.csv(tab,paste0(stem,".csv"),row.names=FALSE)
    genes <- tab$gene[tab$regulation != "NS"]
    if (length(genes)>1) {
      expr <- AverageExpression(sub,assays="RNA",group.by="Group",features=genes,verbose=FALSE)$RNA
      mat <- as.matrix(log2(expr+1))
      mat <- mat[apply(mat,1,sd)>0,,drop=FALSE]
      if (nrow(mat)>1) {
        pdf(paste0(stem,"_heatmap.pdf")); heatmap(mat,labCol=plot_label(colnames(mat)),scale="row",Colv=NA,labRow=NA); dev.off()
      }
    }
    for (direction in c("Up","Down")) {
      genes <- tab$gene[tab$regulation==direction]
      if (length(genes)<3) next
      go <- enrichGO(gene=genes,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP",pvalueCutoff=1,qvalueCutoff=1)
      write.csv(as.data.frame(go),paste0(stem,"_GO_",direction,".csv"),row.names=FALSE)
    }
  }
  if (organ=="Spleen" && requireNamespace("escape",quietly=TRUE) && any(obj$Cell_type=="Neutrophil")) {
    neut <- obj[,colnames(obj)[which(obj$Cell_type=="Neutrophil")]]
    gene_sets <- new.env()
    utils::data("escape.gene.sets",package="escape",envir=gene_sets)
    if (exists("escape.gene.sets",gene_sets,inherits=FALSE)) {
      mouse_sets <- lapply(gene_sets$escape.gene.sets,stringr::str_to_title)
      scores <- escape::escape.matrix(neut,gene.sets=mouse_sets,method="AUCell",groups=1000,min.size=3)
      saveRDS(scores,file.path(out,"neutrophil_AUCell_scores.rds"))
    }
  }
  if (requireNamespace("CommPath",quietly=TRUE)) {
    for (g in levels(obj$Group)) {
      group_obj <- obj[,colnames(obj)[obj$Group==g]]
      comm <- CommPath::createCommPath(expr.mat=as.matrix(GetAssayData(group_obj,assay="RNA",slot="data")),
                                      cell.info=as.character(group_obj$Cell_type),species="mmusculus")
      comm <- CommPath::findLRmarker(comm,method="wilcox.test")
      comm <- CommPath::findLRpairs(comm,logFC.thre=0,p.thre=.01)
      saveRDS(comm,file.path(out,paste0("CommPath_",g,".rds")))
    }
  }
}
