library(tidyverse)
library(limma)

## 1) 读入 Olink NPX 表
dat <- read.table(
  "~/Downloads/PLAAT.olink.NPX.tsv",      # 换成你的文件名
  header = TRUE,
  sep    = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NULL", "NaN", "", "NA")   # 关键：把这些都当成 NA
)

## 看一下结构，确认分组列和蛋白列
str(dat)
head(dat)

## 2) 选出蛋白列：去掉分组和样本 ID 列
protein_cols <- setdiff(colnames(dat), c("Assay","treatment"))

## 3) 把蛋白列全部强制转为 numeric
dat <- dat %>%
  mutate(
    across(
      all_of(protein_cols),
      ~ as.numeric(.x)
    )
  )

## 再检查一下每列的类型是否都是 numeric
sapply(dat[, protein_cols], class)

## 4) 构建表达矩阵：行=蛋白，列=样本
expr <- t(as.matrix(dat[, protein_cols]))
colnames(expr) <- dat$sample
rownames(expr) <- protein_cols

## 再确认一下 expr 的类型
is.numeric(expr)   # 应该返回 TRUE
mode(expr)         # 应该是 "numeric"

## 5) 分组信息（保证对照在前，病例在后）
group <- factor(dat$treatment, levels = c("HC", "Patient"))

design <- model.matrix(~ 0 + group)
colnames(design) <- c("HC", "Patient")

## 6) 用 limma 做差异分析（NA 是允许的） ----------------------------
fit  <- lmFit(expr, design)
cont <- makeContrasts(Patient_vs_HC = Patient - HC, levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

res <- topTable(
  fit2,
  coef          = "Patient_vs_HC",
  number        = Inf,
  adjust.method = "BH"
)

res$Assay <- rownames(res)
res <- res %>% relocate(Assay)

write.table(
  res,
  file      = "~/Downloads/olink_DE_limma.cleaned.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

