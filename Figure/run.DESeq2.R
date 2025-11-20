#!/usr/bin/env Rscript

################################
## 0. 解析命令行参数
################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("用法:\n")
  cat("  Rscript run_DESeq2_PLAAT4.R <count_mat.tsv> <sample_info.tsv> <out_prefix>\n\n")
  cat("示例:\n")
  cat("  Rscript run_DESeq2_PLAAT4.R ../PLAAT4.37samples.RNA-seq.count.tsv sample_info.tsv PLAAT4_DESeq2\n")
  quit(save = "no", status = 1)
}

count_file  <- args[1]
sample_file <- args[2]
out_prefix  <- args[3]

cat("计数矩阵文件: ", count_file,  "\n")
cat("样本信息文件: ", sample_file, "\n")
cat("输出前缀:     ", out_prefix,  "\n")

################################
## 1. 加载 R 包
################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(tidyverse)
})

################################
## 2. 读入计数矩阵和样本信息
################################

# 2.1 读 counts
count_mat <- read.table(
  file        = count_file,
  header      = TRUE,
  row.names   = 1,
  sep         = "\t",
  check.names = FALSE
)

cat("原始计数矩阵维度: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples\n")

# 2.2 读 sample_info
coldata <- read.table(
  file             = sample_file,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

if (!("sample_id" %in% colnames(coldata))) {
  stop("样本信息文件必须包含一列名为 'sample_id' 的列，与计数矩阵列名对应。")
}

################################
## 3. 对齐样本：自动取交集并排序
################################

# 3.1 样本交集
samples_counts <- colnames(count_mat)
samples_meta   <- coldata$sample_id

common_samples <- intersect(samples_counts, samples_meta)

cat("计数矩阵中的样本数: ", length(samples_counts), "\n")
cat("样本信息中的样本数: ", length(samples_meta),   "\n")
cat("两者交集样本数:     ", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("计数矩阵列名与样本信息中的 sample_id 完全没有交集，请检查文件。")
}

# 3.2 只保留交集样本，并按照计数矩阵的顺序对齐
count_mat <- count_mat[, common_samples, drop = FALSE]

coldata_filtered <- coldata[coldata$sample_id %in% common_samples, , drop = FALSE]
# 按照 count_mat 列顺序重排 coldata
coldata_filtered <- coldata_filtered[match(common_samples, coldata_filtered$sample_id), , drop = FALSE]

# 设置行名为 sample_id
rownames(coldata_filtered) <- coldata_filtered$sample_id

# 再次检查顺序是否完全一致
if (!all(colnames(count_mat) == rownames(coldata_filtered))) {
  stop("对齐后 colnames(count_mat) 与 rownames(coldata_filtered) 仍不一致，请检查。")
}

cat("对齐后计数矩阵维度: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples\n")

################################
## 4. 处理分组和协变量
################################

coldata <- coldata_filtered

# condition
if (!("condition" %in% colnames(coldata))) {
  stop("样本信息中必须包含 'condition' 列（如 case/control）。")
}
coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = "control")  # 假定 control 为对照

# sex
if (!("sex" %in% colnames(coldata))) {
  stop("样本信息中必须包含 'sex' 列。")
}
coldata$sex <- factor(coldata$sex)
coldata$sex <- relevel(coldata$sex, ref = "F")  # 以 F 为参考

# age
if (!("age" %in% colnames(coldata))) {
  stop("样本信息中必须包含 'age' 列。")
}
coldata$age <- as.numeric(coldata$age)

################################
## 5. 构建 DESeqDataSet
################################

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ sex + age + condition
)

################################
## 6. 过滤低表达基因
################################

# 按你之前的规则：至少在 min_group_n 个样本中 counts >= 10
# 总样本数
n_samples <- ncol(dds)

# 要求至少 20% 的样本中 counts >= 10
prop_threshold   <- 0.20
min_samples_n    <- ceiling(prop_threshold * n_samples)

cat("总样本数: ", n_samples, "\n")
cat("过滤条件: 至少 ", min_samples_n,
    " 个样本 counts >= 10 (约等于 ", prop_threshold*100, "% 样本)\n")

keep <- rowSums(counts(dds) >= 10) >= min_samples_n
dds  <- dds[keep, ]

cat("过滤后保留基因数: ", nrow(dds), "\n")

################################
## 7. 运行 DESeq2
################################

dds <- DESeq(dds)

################################
## 8. 提取 condition 的差异分析结果
################################

# 自动找到 condition 的系数名
coef_names <- resultsNames(dds)
cat("DESeq2 系数名称:\n")
print(coef_names)

coef_condition <- grep("^condition_", coef_names, value = TRUE)

if (length(coef_condition) == 0) {
  stop("没有找到以 'condition_' 开头的系数名，请检查 design 或 condition 设置。")
}
if (length(coef_condition) > 1) {
  warning("找到多个以 'condition_' 开头的系数，默认使用第一个：", coef_condition[1])
  coef_condition <- coef_condition[1]
}

cat("用于比较的系数名: ", coef_condition, "\n")

res <- results(
  dds,
  name  = coef_condition,
  alpha = 0.05
)

################################
## 9. log2FC 收缩
################################

res_shrunk <- lfcShrink(
  dds,
  coef = coef_condition,
  type = "apeglm"
)

res_shrunk_ordered <- res_shrunk[order(res_shrunk$padj), ]
res_df <- as.data.frame(res_shrunk_ordered)

# ---- 这里改成 TSV 输出 ----
out_tsv_all <- paste0(out_prefix, ".DESeq2_lfcShrink.tsv")
write.table(
  res_df,
  file      = out_tsv_all,
  sep       = "\t",
  quote     = FALSE,
  row.names = TRUE
)

cat("已输出 DESeq2 全基因结果 (TSV): ", out_tsv_all, "\n")

################################
## 10. 质量控制图：PCA & MA
################################

vsd <- vst(dds, blind = FALSE)

# PCA
pca_pdf <- paste0(out_prefix, ".PCA_vst_condition_sex.pdf")
pdf(pca_pdf, width = 6, height = 5)
print(plotPCA(vsd, intgroup = c("condition", "sex")))
dev.off()
cat("已输出 PCA 图: ", pca_pdf, "\n")

# MA
ma_pdf <- paste0(out_prefix, ".MA_case_vs_control.pdf")
pdf(ma_pdf, width = 6, height = 5)
plotMA(res_shrunk, ylim = c(-5, 5))
dev.off()
cat("已输出 MA 图: ", ma_pdf, "\n")

################################
## 11. 提取显著差异基因子集
################################

sig_res <- res_df %>%
  dplyr::filter(!is.na(padj),
                padj < 0.05,
                abs(log2FoldChange) > 1)

out_tsv_sig <- paste0(out_prefix, ".DESeq2_sig_genes_padj0.05_log2FC1.tsv")
write.table(
  sig_res,
  file      = out_tsv_sig,
  sep       = "\t",
  quote     = FALSE,
  row.names = TRUE
)

cat("已输出显著差异基因列表 (TSV): ", out_tsv_sig, "\n")
cat("=== DESeq2 分析完成 ===\n")
