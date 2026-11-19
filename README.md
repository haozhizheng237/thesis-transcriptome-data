# 毕业论文转录组数据

本仓库收录纳米颗粒研究中三种细胞系及肝脏、脾脏的转录组数据。文件保留原始命名。

| 目录 | 内容 |
| --- | --- |
| Cell_Bulk | AML12、293T、TCMK-1 的计数矩阵、差异表达结果及细胞系样本信息，包括 AML12 的 UA 干预计数矩阵 |
| Organ_Bulk | 肝脏和脾脏的计数矩阵、差异表达结果，以及肝脏 UA 干预计数矩阵 |
| scRNA_data | 肝脏和脾脏的两个单细胞 RDS 对象，通过 Git LFS 存储 |

安装 Git LFS 后，在克隆的仓库目录运行 `git lfs pull` 下载 RDS 文件。可在 R 中使用 `readRDS("scRNA_data/Liver_scRNA_data.rds")` 读取肝脏对象。
