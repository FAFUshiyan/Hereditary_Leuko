#!/usr/bin/env Rscript

#############################################
## 0. 解析命令行参数
#############################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("用法:\n")
  cat("  Rscript plot_GO_selected.R <GO_enrich.tsv> <GO_terms.txt> [out_prefix]\n\n")
  cat("参数说明:\n")
  cat("  <GO_enrich.tsv> : clusterProfiler 等产生的 GO 富集结果 TSV 文件\n")
  cat("  <GO_terms.txt>  : 想展示的 GO term 列表文件, 每行一个 Description\n")
  cat("  [out_prefix]    : (可选) 输出文件前缀, 默认为 'GO_selected'\n")
  quit(save = "no", status = 1)
}

go_file   <- args[1]
term_file <- args[2]
out_prefix <- ifelse(length(args) >= 3, args[3], "GO_selected")

cat("GO 富集结果文件: ", go_file,   "\n")
cat("GO term 列表文件:", term_file, "\n")
cat("输出前缀:        ", out_prefix, "\n\n")

#############################################
## 1. 加载 R 包
#############################################

suppressPackageStartupMessages({
  library(tidyverse)
})

#############################################
## 2. 读入数据
#############################################

# 2.1 读入 GO 富集结果
go_df <- read.table(
  go_file,
  header = TRUE,
  sep    = "\t",
  quote  = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2.2 读入想展示的 Description 列表
keep_terms <- readLines(term_file)
keep_terms <- keep_terms[keep_terms != ""]  # 去掉空行

cat("从 term 列表中读取到 ", length(keep_terms), " 个 GO 描述.\n")

#############################################
## 3. 计算 GeneRatio 数值列
#############################################

# GeneRatio 通常是形如 "41/517" 的字符串, 这里把它转成数值
parse_gene_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(v) {
    v <- trimws(v)
    if (length(v) == 2) {
      as.numeric(v[1]) / as.numeric(v[2])
    } else {
      # 如果本身已经是数字, 或格式异常, 就直接当作数值
      as.numeric(v[1])
    }
  })
}

if (!"GeneRatio" %in% colnames(go_df)) {
  stop("输入的 GO 结果文件中没有 'GeneRatio' 列, 请检查文件格式.\n")
}

if (!"p.adjust" %in% colnames(go_df)) {
  stop("输入的 GO 结果文件中没有 'p.adjust' 列, 请检查文件格式.\n")
}

if (!"Count" %in% colnames(go_df)) {
  stop("输入的 GO 结果文件中没有 'Count' 列, 请检查文件格式.\n")
}

go_df <- go_df %>%
  mutate(
    GeneRatio_num   = parse_gene_ratio(GeneRatio),
    minus_log10_FDR = -log10(p.adjust)
  )

#############################################
## 4. 挑选需要展示的 term, 并按 GeneRatio 排序
#############################################

if (!"Description" %in% colnames(go_df)) {
  stop("输入的 GO 结果文件中没有 'Description' 列, 请检查文件格式.\n")
}

go_sel <- go_df %>%
  filter(Description %in% keep_terms)

cat("在富集结果中找到了 ", nrow(go_sel), " 个匹配的 GO term.\n")
missing_terms <- setdiff(keep_terms, go_sel$Description)
if (length(missing_terms) > 0) {
  cat("警告: 有 ", length(missing_terms), " 个 term 在富集结果中没有找到:\n")
  cat("  ", paste(missing_terms, collapse = "; "), "\n\n")
}

# 按 GeneRatio 排序 (从小到大), 然后用这个顺序做 factor,
# 图上 y 轴从下到上就是 GeneRatio 由小到大 (上面的是 GeneRatio 最大的)
go_sel <- go_sel %>%
  arrange(GeneRatio_num) %>%
  mutate(Description = factor(Description, levels = Description))

#############################################
## 5. 画图: barplot + dotplot (PDF 输出)
#############################################

library(ggplot2)

## 5.1 Barplot
pdf(paste0(out_prefix, ".GO_barplot.pdf"), width = 8, height = 6)

p_bar <- ggplot(go_sel,
                aes(x = GeneRatio_num,
                    y = Description,
                    fill = p.adjust)) +
  geom_col(width = 0.7) +
  scale_fill_gradient(
    name  = "Adjusted\nP value",
    low   = "#d73027",
    high  = "#4575b4"
  ) +
  labs(
    title = "Selected GO Enrichment (barplot)",
    x     = "GeneRatio (k / n)",
    y     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )

print(p_bar)
dev.off()

## 5.2 Dotplot
pdf(paste0(out_prefix, ".GO_dotplot.pdf"), width = 8, height = 6)

p_dot <- ggplot(go_sel,
                aes(x = GeneRatio_num,
                    y = Description,
                    size = Count,
                    color = p.adjust)) +
  geom_point() +
  scale_color_gradient(
    name  = "Adjusted\nP value",
    low   = "#d73027",
    high  = "#4575b4"
  ) +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Selected GO Enrichment (dotplot)",
    x     = "GeneRatio (k / n)",
    y     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.text.y  = element_text(size = 9)
  )

print(p_dot)
dev.off()

cat("完成: 已输出 PDF:\n")
cat("  ", paste0(out_prefix, ".GO_barplot.pdf"), "\n")
cat("  ", paste0(out_prefix, ".GO_dotplot.pdf"), "\n")
