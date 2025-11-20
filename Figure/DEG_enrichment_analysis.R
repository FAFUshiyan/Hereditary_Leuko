library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(enrichplot)
library(ggplot2)


file_up   <- "genes_up.txt"    # 上调基因
file_down <- "genes_down.txt"  # 下调基因
#file_bg   <- "genes_bg.txt"    # 所有 filter 之后的背景基因
genes_up_raw   <- unique(scan(file_up,   what = "character", quiet = TRUE))
genes_down_raw <- unique(scan(file_down, what = "character", quiet = TRUE))
cat("原始上调基因数:",   length(genes_up_raw),   "\n")
cat("原始下调基因数:",   length(genes_down_raw), "\n")
#cat("原始背景基因数:",   length(genes_bg_raw),   "\n")



if (length(genes_up_raw) > 0) {
  ego_bp_up <- enrichGO(
    gene          = genes_up_raw,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  go_bp_up_df <- as.data.frame(ego_bp_up)
  write.table(
    go_bp_up_df,
    file      = "GO_BP_up_from_list_all.tsv",
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  cat("已输出 GO BP 上调富集结果：GO_BP_up_from_list_all.tsv\n")
} else {
  cat("上调 ENTREZ 基因集为空，跳过 GO BP 上调富集。\n")
}

if (length(genes_down_raw) > 0) {
  ego_bp_down <- enrichGO(
    gene          = genes_down_raw,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  go_bp_down_df <- as.data.frame(ego_bp_down)
  write.table(
    go_bp_down_df,
    file      = "GO_BP_down_from_list_all.tsv",
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  cat("已输出 GO BP 下调富集结果：GO_BP_down_from_list_all.tsv\n")
} else {
  cat("下调 ENTREZ 基因集为空，跳过 GO BP 下调富集。\n")
}


genes_up <- bitr(
  genes_up_raw,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

genes_up_entrez <- unique(na.omit(genes_up$ENTREZID))

if (length(genes_up_entrez) > 0) {
  react_up <- enrichPathway(
    gene          = genes_up_entrez,
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 1,
    readable      = TRUE
  )

  react_up_df <- as.data.frame(react_up)
  write.table(
    react_up_df,
    file      = "Reactome_up_from_list_all.tsv",
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  cat("已输出 Reactome 上调富集结果：Reactome_up_from_list_all.tsv\n")
} else {
  cat("上调 ENTREZ 基因集为空，跳过 Reactome 上调富集。\n")
}

genes_down <- bitr(
  genes_down_raw,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

genes_down_entrez <- unique(na.omit(genes_down$ENTREZID))

## 5.2 下调 Reactome
if (length(genes_down_entrez) > 0) {
  react_down <- enrichPathway(
    gene          = genes_down_entrez,
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 1,
    readable      = TRUE
  )

  react_down_df <- as.data.frame(react_down)
  write.table(
    react_down_df,
    file      = "Reactome_down_from_list_all.tsv",
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  cat("已输出 Reactome 下调富集结果：Reactome_down_from_list_all.tsv\n")
} else {
  cat("下调 ENTREZ 基因集为空，跳过 Reactome 下调富集。\n")
}

cat("=== 基于 2 个 gene list 的 GO + Reactome 富集分析完成 ===\n")


###############################
## 1. 读入 DESeq2 结果
###############################

# 输入文件：DESeq2_case_vs_control_sex_age_lfcShrink.csv
res_df <- read.csv(
  "DESeq2_case_vs_control_with_sex_age_lfcShrink.csv",
  row.names   = 1,
  check.names = FALSE
)

# 确认关键列存在
stopifnot("log2FoldChange" %in% colnames(res_df))

cat("DESeq2 结果行数(基因数):", nrow(res_df), "\n")

###############################
## 2. ENSEMBL ID 处理 & ID 映射
###############################

# 如果 ENSEMBL ID 带 .13 这种版本号，先去掉
res_df2 <- res_df %>%
  rownames_to_column("GENE_ID") %>%
  mutate(GENE_ID = sub("\\.\\d+$", "", GENE_ID))

FROM_TYPE <- "ENSEMBL"  # 你的基因 ID 类型

# ENSEMBL -> ENTREZID + SYMBOL
id_map <- bitr(
  res_df2$GENE_ID,
  fromType = FROM_TYPE,
  toType   = c("ENTREZID", "SYMBOL"),
  OrgDb    = org.Hs.eg.db
)

# 改名为统一的 GENE_ID
id_map <- id_map %>%
  rename(GENE_ID = !!FROM_TYPE) %>%
  distinct(GENE_ID, .keep_all = TRUE)

cat("成功映射到 ENTREZID 的 GENE_ID 数:", nrow(id_map), "\n")

# 合并回 DESeq2 结果，只保留有 ENTREZID 的行
res_annot <- res_df2 %>%
  left_join(id_map, by = "GENE_ID") %>%
  filter(!is.na(ENTREZID), !is.na(log2FoldChange))

cat("用于 GSEA 的行数(有ENTREZID+log2FC):", nrow(res_annot), "\n")

###############################
## 3. 构造 GSEA 所需 geneList
###############################
# 要求: 每个 ENTREZID 只出现一次 & 按数值递减排序

gsea_df <- res_annot %>%
  group_by(ENTREZID) %>%
  summarise(
    log2FoldChange = mean(log2FoldChange, na.rm = TRUE)
    # 若想更激进，可改为:
    # log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))]
  ) %>%
  ungroup()

cat("按 ENTREZID 合并后的基因数:", nrow(gsea_df), "\n")

geneList <- gsea_df$log2FoldChange
names(geneList) <- gsea_df$ENTREZID

# 去 NA
geneList <- geneList[!is.na(geneList)]

# 按数值从大到小排序 (GSEA 要求)
geneList <- sort(geneList, decreasing = TRUE)

# 检查是否有重复名字（不允许）
stopifnot(length(geneList) == length(unique(names(geneList))))
cat("GSEA geneList 长度:", length(geneList), "\n")

###############################
## 4. 运行 GSEA（GO BP）
###############################

gsea_go <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",       # Biological Process
  minGSSize    = 10,
  maxGSSize    = 500,
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  verbose      = TRUE
)

gsea_go_df <- as.data.frame(gsea_go)

# 输出完整结果表
write.table(
  gsea_go_df,
  file      = "GSEA_GO_BP_from_DESeq2_all_terms.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

cat("已输出 GSEA GO BP 结果：GSEA_GO_BP_from_DESeq2_all_terms.tsv\n")

###############################
## 5. 可选：画一张 top 通路 GSEA 曲线图 (PDF)
###############################

if (nrow(gsea_go_df) > 0) {
  # 取 NES 绝对值最大的那条通路
  top_id <- gsea_go_df$ID[which.max(abs(gsea_go_df$NES))]
  top_desc <- gsea_go_df$Description[which.max(abs(gsea_go_df$NES))]

  pdf("GSEA_GO_BP_topPath_gseaplot.pdf", width = 7, height = 5)
  print(
    gseaplot2(
      gsea_go,
      geneSetID = top_id,
      title     = top_desc
    )
  )
  dev.off()
  cat("已输出 top 通路的 GSEA 曲线图：GSEA_GO_BP_topPath_gseaplot.pdf\n")
} else {
  cat("GSEA 未检测到显著通路，未生成 GSEA 曲线图。\n")
}

cat("=== GSEA 分析完成 ===\n")

