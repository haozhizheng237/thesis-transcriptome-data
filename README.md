# 转录组数据

本仓库收录毕业论文中 AML12、293T、TCMK-1 三种细胞系以及肝脏、脾脏的转录组数据、LNP 分析脚本和第三章主要富集结果配套表。

主分析以 Control 为参照，比较 **LNP vs Control**。UA 挽救分析比较 **LNP+UA vs LNP**。图中的 **eLNP** 对应统计分析中的 **LNP** 组。

## 数据与结果

| 位置 | 内容 |
| --- | --- |
| [Cell_Bulk](Cell_Bulk/) | 三种细胞系的计数矩阵、差异表达结果和样本信息，以及 AML12 的 UA 挽救实验计数矩阵 |
| [Organ_Bulk](Organ_Bulk/) | 肝脏、脾脏的计数矩阵和差异表达结果，以及肝脏 UA 挽救实验计数矩阵 |
| [scRNA_data](scRNA_data/) | 肝脏、脾脏各一个单细胞 RDS 对象，使用 Git LFS 存储 |
| [第三章主要富集结果配套数据表](tables/第三章主要富集结果配套数据表.xlsx) | 六张工作表：批量 GO、趋势簇 KEGG、批量 KEGG GSEA、跨模型 GO GSEA、单细胞 GO、肝脏 UA 逆转 |

## 分析脚本

| 脚本 | 用途 |
| --- | --- |
| [细胞系 bulk 分析](scripts/多组差异基因富集分析%28CellBulk%29.R) | 按细胞系比较 LNP 与 Control，分析差异表达和通路富集 |
| [肝脾富集分析](scripts/多组差异基因富集分析%28LiverSpleen%29.R) | 从肝、脾差异表达表提取上下调基因，进行 GO 富集，并绘制 LNP 辐射强度时间曲线 |
| [肝脾 bulk 流程](scripts/liver_spleen_full_pipeline_v2.R) | 从计数、TPM 和批次信息开展差异分析、样本相似度分析、富集分析和绘图 |
| [肝脾单细胞分析](scripts/Liver_Spleen_scRNA.R) | 统计样本的细胞类型组成，按细胞类型比较 LNP 与 Control |
| [CellChat 通讯网络](scripts/scRNA_cellcell_communication_whole_network.R) | 按组构建 CellChat 对象并汇总通讯网络 |
| [CellChat 通路面板](scripts/CellChat_pathway_panels_Spleen_Liver.R) | 读取上一步生成的肝、脾 CellChat 对象，绘制通路和配体–受体通讯图 |
| [UA 挽救实验](scripts/Wet_experiments.R) | 分析 AML12 和肝脏的 UA 挽救实验，筛选两个比较中显著且变化方向相反的基因 |

## 运行示例

安装 Git LFS、R 和相应脚本调用的 R 包后，在仓库根目录运行以下命令。示例读取肝脏和脾脏的差异表达表，结果写入工作目录下的 `LNP_results/`。

```bash
git lfs pull
export THESIS_DATA="$PWD"
export LNP_RESULTS="$PWD/LNP_results"
Rscript 'scripts/多组差异基因富集分析(LiverSpleen).R'
```

肝脾 bulk 完整流程还需要通过 `ORGAN_PREPARED` 指定含 `raw_count_matrix.rds`、`TPM_matrix.rds` 和 `sample_metadata.csv` 的预处理目录。该流程的 KEGG GSEA 使用 `MOUSE_KEGG_GMT` 指定小鼠基因集文件。细胞系关键基因热图可通过 `CELL_TPM_DIR` 指定细胞系 TPM 矩阵目录。这些额外输入需另行提供。
