## 如果还没装这些包，先运行：
## install.packages(c("ggplot2", "ggsignif"))

library(ggplot2)
library(ggsignif)

## 1. 把你在问题里贴的表保存为一个 TSV 文件，例如 cytokine.tsv
##    注意用 Tab 分隔，并保留表头这一行。

dat <- read.table(
  "~/Downloads/factor.tsv",          # 修改成你的文件名
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  na.strings = c("NULL")   # 把字符串 "NULL" 当成 NA
)

## 2. 构建分组因子（保证 HC 在左，Patient 在右）
dat$treatment <- factor(dat$treatment, levels = c("HC", "Patient"))

## 3. 选择要画的指标，这里以 CSF1 为例 --------------------------
##    想画其他指标（如 "IFN-γ"、"IL17C"）只要把下面的 CSF1 换掉即可

marker <- "IL17C"

## 4. 计算两组的 P 值（这里用 Wilcoxon 检验；如需 t 检验改成 t.test）
f <- as.formula(paste(marker, "~ treatment"))
wilcox_res <- wilcox.test(f, data = dat)
p_val   <- wilcox_res$p.value
p_label <- paste0("P = ", signif(p_val, 3))

## 5. 设定 y 轴最大值，方便放 P 值标注
y_max <- max(dat[[marker]], na.rm = TRUE)

## 6. 作图对象 -------------------------------------------------------
p <- ggplot(dat, aes(x = treatment, y = .data[[marker]], fill = treatment)) +
  # 小提琴
  geom_violin(
    width = 0.9,
    alpha = 0.3,
    color = NA,
    trim  = FALSE
  ) +
  # 箱线
  geom_boxplot(
    width = 0.25,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  # 单个样本点
  geom_jitter(
    width = 0.08,
    size  = 1.5,
    alpha = 0.7
  ) +
  # P 值括号
  geom_signif(
    comparisons = list(c("HC", "Patient")),
    annotations = p_label,
    y_position  = y_max * 1.07,
    tip_length  = 0.02,
    textsize    = 3.5
  ) +
  scale_fill_manual(values = c("HC" = "#E64B35FF", "Patient" = "#4DBBD5FF")) +
  labs(
    x = "",
    y = paste0(marker, " (A.U.)")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 11),
    axis.text.y     = element_text(size = 11)
  )

## 7. 输出 PDF -------------------------------------------------------
pdf("~/Downloads/IL17C_HC_vs_Patient_violin.pdf", width = 3.2, height = 3.5)
print(p)
dev.off()
