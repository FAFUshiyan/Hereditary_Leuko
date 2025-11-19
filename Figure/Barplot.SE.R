# ===================================================================
# 1. 设置和加载库
# ===================================================================
library(tidyverse)
library(scales)
library(RColorBrewer) # 用于生成颜色

# ===================================================================
# 2. !! 用户需要修改的部分 !!
# ===================================================================
# 请将路径替换为您的实际文件路径
all_data_file <- "/Users/shiyan/Library/CloudStorage/OneDrive-Personal/2025工作目录/工作汇报/PLAAT4/Figure-project/metabilism_Figure/mice_lipid_reanalysis_Oct_29/堆积图/WM.mice.data.tsv"          # 亚类分类文件
classification_file   <- "/Users/shiyan/Downloads/Barplot/Mice/WM.CAR.class.tsv"       # 包含所有脂质数据的TSV文件
output_dir          <- "~/Downloads/Barplot/Mice/WM"      # 图片保存的文件夹名称

# ===================================================================
# 3. 解析亚类分类文件
# ===================================================================
# 创建输出文件夹
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 读取所有行，并移除空行
lines <- readLines(classification_file)
lines <- lines[lines != ""]

# 提取每行开头的亚类前缀 (如 "BMP", "CAR", "CE")
prefixes <- sub("\\s.*", "", lines)

# 根据前缀将所有脂质ID分组到一个命名列表中
# 列表的名称是亚类 (e.g., "BMP"), 列表的内容是该亚类的所有脂质ID
lipid_subclasses <- split(lines, prefixes)

cat("成功解析出以下脂质亚类:\n")
print(names(lipid_subclasses))

# ===================================================================
# 4. 读入总数据文件
# ===================================================================
# 读入数据 (保留带空格/括号的列名)
df_all <- read.table(all_data_file, sep = "\t", header = TRUE, check.names = FALSE)
cat(sprintf("\n成功读入总数据文件，共 %d 行, %d 列。\n", nrow(df_all), ncol(df_all)))

# ===================================================================
# 5. 循环绘图
# ===================================================================
cat("\n开始循环绘制每个亚类的图表...\n")

# 遍历解析出的每一个亚类
for (subclass_name in names(lipid_subclasses)) {
  
  # 获取当前亚类的所有脂质ID
  current_lipid_ids <- lipid_subclasses[[subclass_name]]
  
  # 检查数据文件中是否存在这些脂质列，只保留存在的
  lipids_in_data <- intersect(current_lipid_ids, colnames(df_all))
  
  if (length(lipids_in_data) == 0) {
    cat(sprintf("-> 跳过亚类 '%s'，因为在数据文件中找不到对应的脂质列。\n", subclass_name))
    next # 如果一个都没有，就跳到下一个循环
  }
  
  cat(sprintf("-> 正在处理亚类: %s (%d 个脂质)\n", subclass_name, length(lipids_in_data)))
  
  # ---- A. 数据准备和筛选 ----
  # 从总数据中只选择当前亚类相关的列，以及样本信息列
  df_subset <- df_all %>%
    select(Species, treatment, all_of(lipids_in_data))
  
  # ---- B. 执行您之前的绘图逻辑 ----
  # 宽转长
  long <- df_subset %>%
    pivot_longer(
      cols = all_of(lipids_in_data),
      names_to = "lipid", values_to = "abundance"
    )
  
  # 组内统计
  stats <- long %>%
    group_by(treatment, lipid) %>%
    summarise(
      mean = mean(abundance, na.rm = TRUE),
      sd   = sd(abundance,   na.rm = TRUE),
      n    = n_distinct(Species),
      se   = sd / sqrt(n),
      .groups = "drop"
    )
  
  # 以 control 的均值从大到小排序物种；固定组别顺序
  lipid_order <- stats %>% filter(treatment == "control") %>%
    arrange(desc(mean)) %>% pull(lipid)
  
  stats <- stats %>%
    mutate(lipid = factor(lipid, levels = lipid_order),
           treatment = forcats::fct_relevel(treatment, "control", "case")) %>%
    arrange(treatment, lipid) %>%
    group_by(treatment) %>%
    mutate(xstart = cumsum(dplyr::lag(mean, default = 0)),
           xend   = xstart + mean,
           err    = se,
           xmin_err = xend - err,
           xmax_err = xend + err) %>%
    ungroup()
  
  # ---- C. 动态生成颜色 ----
  # 为当前亚类的脂质动态生成足够数量的颜色
  num_lipids <- length(lipids_in_data)
  if (num_lipids <= 12) {
    # 如果脂质少，可以用预设的高质量调色板
    pal <- brewer.pal(max(3, num_lipids), "Paired")[1:num_lipids]
  } else {
    # 如果脂质多，就生成渐变色
    pal <- colorRampPalette(brewer.pal(9, "Set1"))(num_lipids)
  }
  # 确保颜色名称与lipid因子水平对应
  names(pal) <- lipid_order
  
  # ---- D. 绘图 ----
  p <- ggplot(stats, aes(y = treatment)) +
    geom_col(aes(x = mean, fill = lipid), width = 0.62,
             position = position_stack(reverse = TRUE)) +
    geom_errorbarh(aes(xmin = xmin_err, xmax = xmax_err),
                   height = 0.02, linewidth = 0.4, colour = "grey51") +
    scale_fill_manual(values = pal, name = NULL, breaks = lipid_order) +
    labs(x = "Mean abundance", y = NULL, title = sprintf("Lipid Subclass: %s", subclass_name)) +
    theme_minimal(base_size = 12)+
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0),
      axis.line = element_blank(),      # 取消轴线，避免与边框叠加变粗
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1) # 只保留统一的外框
    )

  
  # ---- E. 保存图像 ----
  output_filename <- file.path(output_dir, sprintf("plot_%s.pdf", subclass_name))
  ggsave(output_filename, p, width = 10, height = 4 + num_lipids * 0.05)
}

cat("\n所有图表绘制完成！\n")
