# 可用 THESIS_DATA 和 LNP_RESULTS 指定输入和输出目录。
data_root <- path.expand(Sys.getenv("THESIS_DATA", "~/Downloads/毕业论文数据"))
result_root <- path.expand(Sys.getenv("LNP_RESULTS", "LNP_results"))
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
# 仅在图形显示层使用 eLNP。
plot_label <- function(x) gsub("(?<![A-Za-z])LNP", "eLNP", as.character(x), perl = TRUE)

library(ggplot2)
out <- file.path(result_root,"Organ_enrichment")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
radiance <- data.frame(
  time = c(0,12,24,36,48,60,72,84,96),
  LNP_1 = c(2000000000, 20900000000, 32900000000, 33200000000, 22900000000, 12300000000, 9340000000, 7650000000, 5230000000),
  LNP_2 = c(1890000000, 24100000000, 38400000000, 36600000000, 26300000000, 21900000000, 8890000000, 6860000000, 5110000000),
  LNP_3 = c(1930000000, 18600000000, 34900000000, 38900000000, 22800000000, 9340000000, 6340000000, 6770000000, 5440000000)
)
values <- as.matrix(radiance[,paste0("LNP_",1:3)])
curve <- data.frame(time = c(0,4,8,12,24,36,48,60,72), mean=rowMeans(values), sd=apply(values,1,sd))
write.csv(curve,file.path(out,"LNP_radiance.csv"),row.names=FALSE)
p <- ggplot(curve,aes(time,mean)) + geom_line(colour="#FF7F0E") + geom_point() +
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=1) + labs(x="Time (h)",y="Total radiance") + theme_classic()
ggsave(file.path(out,"LNP_radiance.pdf"),p)
for (organ in c("Liver","Spleen")) {
  tab <- read.csv(file.path(data_root,"Organ_Bulk",paste0("eLNP_",organ,"_DEG.csv")), check.names=FALSE)
  stopifnot(all(c("gene","logFC","P.Value") %in% names(tab)))
  for (direction in c("Up","Down")) {
    selected <- tab$P.Value < 0.05 & (if (direction=="Up") tab$logFC >= 1 else tab$logFC <= -1)
    genes <- unique(na.omit(tab$gene[which(selected)]))
    stem <- file.path(out,paste(organ,"LNP_vs_Control",direction,sep="_"))
    write.csv(tab[which(selected),], paste0(stem,"_genes.csv"), row.names=FALSE)
    if (length(genes)<3) next
    go <- clusterProfiler::enrichGO(gene=genes, universe=unique(tab$gene), OrgDb=org.Mm.eg.db::org.Mm.eg.db,
                keyType="SYMBOL",ont="BP",pvalueCutoff=1,qvalueCutoff=1,readable=TRUE)
    write.csv(as.data.frame(go),paste0(stem,"_GO_BP.csv"),row.names=FALSE)
    if (nrow(as.data.frame(go))) ggsave(paste0(stem,"_GO_BP.pdf"),enrichplot::dotplot(go,showCategory=10))
  }
}
