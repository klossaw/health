
# Header1
## Header2

```{r}
knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo=FALSE) 
pkgs <- c("fs", "stringr", "ggpubr", "ggthemes", 
          "glue", "tidyverse", "dplyr", "DESeq2", "kableExtra", 
          "readxl", "org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler", 
          "ComplexHeatmap", "png", "knitr", "jhuanglabGO", "ggforce", 
          "FactoMineR", "factoextra", "RColorBrewer", "optparse")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- params$project
dataset <- params$dataset
species <- params$species
case_group <- params$case
ctrl_group <- params$control
batch <- params$batch
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/rnaseq")
workdir %>% fs::dir_create() %>% setwd()

```


```{r}
## load tables
cnt_dat_fn <- glue("~/projects/{project}/analysis/{dataset}/{species}/rnaseq/tables/{dataset}_{species}_counts.csv")
count_dat <- read_csv(cnt_dat_fn)
fpkm_fn <- glue("~/projects/{project}/analysis/{dataset}/{species}/rnaseq/tables/{dataset}_{species}.csv")
fpkm <- read_csv(fpkm_fn)

## Specify control and experimental
sampleinfo_fn <- glue("~/projects/{project}/docs/{dataset}/sampleinfo_{dataset}.xlsx")
sampleinfo <- readxl::read_excel(sampleinfo_fn, sheet = "rnaseq")
case <- sampleinfo %>% dplyr::filter(grepl(case_group, .$sample_id), grepl(batch, .$sample_id)) %>% pull(sample_id)
control <- sampleinfo %>% dplyr::filter(grepl(ctrl_group, .$sample_id), grepl(batch, .$sample_id)) %>% pull(sample_id)
contrast <- glue("{batch}_{case_group}_vs_{ctrl_group}")
# specify which output directory
output_dir <- glue::glue("~/projects/{project}/analysis/{dataset}/{species}/rnaseq/{contrast}")
output_dir %>% fs::dir_create() %>% setwd()

```

## 基因表达的整体情况

```{r}

pca_plot_rna(fpkm = fpkm, case = case, control = control, output_dir = output_dir)
jhuanglabGO::varied_genes_htmap(case, control, fpkm, top_n = 1000, output_dir)
# # when using pdftools to convert the pdf to png/jpeg, the studio sever got corrupted with no idea
pdftools::pdf_convert(pdf = glue("{output_dir}/qc/pca.pdf"), filenames = "pca.png", dpi = 300)
fs::file_move("pca.png", new_path = glue("{output_dir}/qc/pca.png"))
pdftools::pdf_convert(pdf = glue("{output_dir}/qc/top1000_htmap.pdf"), filenames = "top1000_htmap.png", dpi = 300)
fs::file_move("top1000_htmap.png", new_path = glue("{output_dir}/qc/top1000_htmap.png"))

```

主成分分析 (PCA)常用于评估组间差异及组内样本重复情况，PCA采用线性代数的计算方法，对数以万计的基因变量进行降维及主成分提取。我们对所有样本的基因表达值进行PCA分析，如下图所示。

```{r.width='80%'} 
knitr::include_graphics(glue("{output_dir}/qc/pca.png"))
```

结果文件：[mettl3_cre_vs_nc/qc/pca.png](mettl3_cre_vs_nc/qc/pca.png)

为进一步确认组间基因表达模式的区别，本流程绘制了基于前1000个高变基因的聚类热图，如下图所示。图中每一行代表一个高变基因，每一列代表一个检测样本。顶部的颜色代表不同的实验条件，结合横向的层级聚类可观察不同样本或不同实验条件对于基因表达模式的影响。

```{r.width=6, out.height='100%', out.width='80%'}
knitr::include_graphics(glue("{output_dir}/qc/top1000_htmap.png"))
```

结果文件：[mettl3_cre_vs_nc/qc/top1000_htmap.png](mettl3_cre_vs_nc/qc/top1000_htmap.png)

## 差异基因

```{r out.width='80%'}
res <- jhuanglabGO::get_diff_genes(case, control, count_dat, species = species, output_dir)
# volcano_plot_rna(res, output_dir = output_dir)
pdftools::pdf_convert(pdf = glue("{output_dir}/diff/volcano_plot.pdf"), filenames = glue("volcano_plot.png"), dpi = 300)
fs::file_move("volcano_plot.png", new_path = glue("{output_dir}/diff/volcano_plot.png"))

```

基因表达定量完成后，需要对其表达数据进行统计学分析，筛选样本在不同状态下表达水平显著差异的基因。
通常以在两组样品中的表达量差异达到两倍以上 (即｜log2FoldChange｜\> 1)，且多重检验矫正的padj \< 0.05为阈值。本次测序共发现`r res %>% dplyr::filter(abs(log2FoldChange) > 1, padj < 0.05) %>% nrow()`个差异基因，上升的基因数目为`r res %>% dplyr::filter(log2FoldChange > 1, padj < 0.05) %>% nrow()`，下降的基因数目为`r res %>% dplyr::filter(log2FoldChange < -1, padj < 0.05) %>% nrow()`。差异最显著的10个基因如下表所示，其余结果见[mettl3_cre_vs_nc/diff/diff_genes.xlsx](mettl3_cre_vs_nc/diff/diff_genes.xlsx)。

```{r}
res <- readxl::read_excel(glue("{output_dir}/diff/diff_genes.xlsx"))
res %>% head(n = 10) %>% dplyr::select(-11) %>% 
  knitr::kable(caption = "表1：差异基因top 10", row.names = NA, col.names = colnames(.)) %>%
  column_spec(7, width = "22cm")

```

结果文件：[mettl3_cre_vs_nc/diff/diff_genes.xlsx](mettl3_cre_vs_nc/diff/diff_genes.xlsx)

为直观展示每个比较组合的差异基因分布情况，本次分析绘制了差异基因火山图，如下图所示。图中横坐标表示基因在处理和对照两组中的表达倍数变化(log2FoldChange)，纵坐标表示基因在处理和对照两组中表达差异的显著性水平 (-log10padj)。蓝色表示在对照组显著上升的基因，红色代表在处理组显著上升的基因。

```{r}
knitr::include_graphics(glue("{output_dir}/diff/volcano_plot.png"))

```

结果文件：[mettl3_cre_vs_nc/diff/volcano_plot.png](mettl3_cre_vs_nc/diff/volcano_plot.png)

## 基于GO数据库的差异基因富集分析

```{r}
# jhuanglabGO::get_GO_enrichment(diff = res, output_dir = output_dir, species = species)
tds <- c("all", "down", "up")
ctgs <- c("BP", "CC", "MF")
pdf_fn <- expand_grid(tds, ctgs) %>% mutate(pdf_fn = paste0(tds, "_", ctgs, "_barplot")) %>% pull(pdf_fn)
idx <- file.exists(paste0(output_dir, "/go/", pdf_fn, ".pdf"))
pdf_fn <- pdf_fn[idx]
parallel::mclapply(pdf_fn, \(fn){
  pdftools::pdf_convert(pdf = glue("{output_dir}/go/{fn}.pdf"), 
                filenames = glue("{fn}.png"), dpi = 300)
  fs::file_move(glue("{fn}.png"), new_path = glue("{output_dir}/go/"))
}, mc.cores = length(pdf_fn))

```

GO (Gene Ontology)数据库，收集的是对各种物种基因功能进行限定和描述的标准词汇 (term)，是国际标准化的基因功能描述分类系统。根据基因产物的相关生物学过程 (Biological Process)、细胞组分 (Cellular Component)以及分子功能 (Molecular Function)三个大类分别给予定义。

本流程对显著变化基因、显著下降基因与显著上升基因分别进行了基于GO数据库的富集分析，并绘制了柱状图。基于显著变化基因的GO富集分析结果如下图所示。横坐标代表富集到的基因的数目，纵轴为富集到的GO条目，颜色代表多重校正后的p值 (p.adjust)。

```{r}
knitr::include_graphics(glue("{output_dir}/go/all_BP_barplot.png"))

```

结果文件目录：[mettl3_cre_vs_nc/go](mettl3_cre_vs_nc/go)

## 基于KEGG数据库的差异基因富集分析

KEGG (Kyoto Encyclopedia of Genes and Genomes)是一个综合数据库，整合了基因组信息、化学信息和生化系统功能信息，目前包含了16个子数据库。比如，KEGG PATHWAY数据库包含了图解的细胞代谢、膜转运、信号传导等通路信息； KEGG GENES数据库、KEGG GENOME数据库则包含了部分或者完整序列的基因/基因组信息；KEGG Orthology (KO)是KEGG直系同源数据库，将各个KEGG注释系统联系在一起，将分子网络和基因组信息联系起来，根据直系同源关系，实现跨物种的基因组或转录组的功能注释。

```{r}
# enrich_kegg <- jhuanglabGO::get_kegg_enrichment(diff = res, output_dir = output_dir, species = species)

```

结果文件目录：[mettl3_cre_vs_nc/kegg](mettl3_cre_vs_nc/kegg)

## 基因集富集分析 (GSEA)

```{r}
# gsea_dir <- glue::glue("{output_dir}/gsea")
# jhuanglabGSEA_all(dat = fpkm, case_samples = case,
#                   control_samples = control, out_dir = gsea_dir, gmt_version = "v7.5")

```

GSEA (Gene Set Enrichment Analysis)，即基因集富集分析。基本思想是使用预定义的基因集 (通常来自功能注释或先前实验的结果)，将基因按照在两类样本中的差异表达程度排序，然后检验预先设定的基因集合是否在这个排序表的顶端或者底端富集。基因集合富集分析检测基因集合而不是单个基因的表达变化，因此可以包含这些细微的表达变化，预期得到更为理想的结果。

本流程的GSEA使用了博德研究所 (Broad Institute) 开发的GSEAclient软件，结合了分子特征数据库 (The Molecular Signatures Database, MSigDB) 的9个基因集合，从多个角度描述了不同实验条件对细胞类型、免疫状态、基因模式表达的影响。下面两张表格分别表示了一些在处理组或对照组中富集的基因集合。

```{r}
gsea_id <- list.files(glue("{output_dir}/gsea/c5")) %>% stringr::str_sub(start = nchar("jhuanglabGSEA.Gsea.") + 1)
# enriched in Case group
case_fn <- glue("{output_dir}/gsea/c5/jhuanglabGSEA.Gsea.{gsea_id}/gsea_report_for_Case_{gsea_id}.tsv")
gsea_res_case <- readr::read_tsv(case_fn) %>% dplyr::select(c(-2, -3, -12))
gsea_top10_case <- gsea_res_case %>% head(n = 10) %>% 
  knitr::kable(., col.names = colnames(.), align = rep('c', ncol(.)), 
               caption = "表2: 处理组富集的部分基因集合 (c5)")
gsea_top10_case
# case_fn <- glue("{output_dir}/gsea/c5/jhuanglabGSEA.Gsea.{gsea_id}/gsea_report_for_Case_{gsea_id}.html")
# htmltools::includeHTML(path = case_fn)
# encirhed in Control group
# control_fn <- glue("{output_dir}/gsea/c5/jhuanglabGSEA.Gsea.{gsea_id}/gsea_report_for_Control_{gsea_id}.html")
# htmltools::includeHTML(path = control_fn)
control_fn <- glue("{output_dir}/gsea/c5/jhuanglabGSEA.Gsea.{gsea_id}/gsea_report_for_Control_{gsea_id}.tsv")
gsea_res_control <- readr::read_tsv(control_fn) %>% dplyr::select(c(-2, -3, -12))
gsea_top10_control <- gsea_res_control %>% head(n = 10) %>% 
  knitr::kable(., col.names = colnames(.), align = rep('c', ncol(.)), 
               caption = "表2: 对照组富集的部分基因集合 (c5)")
gsea_top10_control

```

结果文件目录：[mettl3_cre_vs_nc/gsea](mettl3_cre_vs_nc/gsea)
