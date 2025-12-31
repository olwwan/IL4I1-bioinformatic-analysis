
# install.packages(c("circlize", "ComplexHeatmap"))

data=data_DEGA
expMatrix <- data
#去除含有大量 0 值的样本
zero_counts <- rowSums(expMatrix == 0)
threshold <- ncol(expMatrix) * 0.5
expMatrix <- expMatrix[zero_counts <= threshold, ]

#log2 转化箱线图
pdf("./1.log2_box.pdf", width = 9, height = 6)
par(mar=c(7, 7, 2, 2))
boxplot(expMatrix,col=group, notch=T,las=2)
dev.off()
#log2 转换，判断数据是否需要转换
ex <- expMatrix
qx <- as.numeric(quantile(ex, c(0.00, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[6] > 100) || # 判断最大值是否大于 100
(qx[6] - qx[1] > 50) # 判断数据范围差
if (LogC) {
expMatrix <- log2(ex + 1)
pdf("1.差异分析/TPM/2.log2 转化后箱线图.pdf", width = 9, height = 6)
par(mar=c(7, 7, 2, 2))
boxplot(expMatrix, notch=T,col=group, las=2)
dev.off()
write.csv(expMatrix, file = "1.差异分析/TPM/tpm_genes_log2.csv", row.names = TRUE)
print("log2 transform finished")
} else {
print("log2 transform not needed")
}
if (是否标准化) {
expMatrix <- normalizeBetweenArrays(expMatrix)
pdf("./1.normalized_box.pdf", width = 9, height = 6)
par(mar=c(7, 3, 2, 2))
boxplot(expMatrix, notch = TRUE, col = group, las = 2)
dev.off()
write.csv(expMatrix, file = "1.差异分析/TPM/tpm_genes_log2_normalized.csv", row.names = TRUE)
}
#开始差异分析
data10 <- expMatrix
# 创建设计矩阵
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design)=colnames(data10)
# 拟合线性模型/对比矩阵
constrasts = paste(rev(levels(group)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design)
fit <- lmFit(data10, design)
fit1=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit1)
limma_DEG = topTable(fit2, coef=constrasts, n=Inf, adjust.method="BH")
limma_DEG = na.omit(limma_DEG)
logFC_cutoff <- with(limma_DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
if (手动设置阈值) {
logFC_cutoff <- cuf
}
message("limma:动态 log2FC 阈值为：|logFC| > [mean(|logFC|) + 2sd(|logFC|)]")
message("本次分析采用的阈值为：")
print(logFC_cutoff)
k1 <- (limma_DEG[[P 值选择]] < 0.05) & (limma_DEG$logFC < -logFC_cutoff)
k2 <- (limma_DEG[[P 值选择]] < 0.05) & (limma_DEG$logFC > logFC_cutoff)
limma_DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
print(table(limma_DEG$change))

write.csv(limma_DEG, file = "./TPM/limma_DEGs.csv", row.names = T)
# 创建 GSEA 分析所需的数据
limma_log2FC_TPM = data.frame(gene_symbol = rownames(limma_DEG), change=limma_DEG$change ,
logFC = limma_DEG$logFC)
write.csv(limma_log2FC_TPM, file = "./TPM/limma_log2FC.csv", row.names = F)
limma_log2FC_fpkm_filtered <- limma_log2FC_TPM %>%
filter(change != "NOT")
write.csv(limma_log2FC_fpkm_filtered, file = "1.差异分析/TPM/limma_上下调基因.csv", row.names =
FALSE)
#非参数秩和检验
count_norm <- as.data.frame(expMatrix)
# 对每个基因进行 Wilcoxon 秩和检验
plan(multisession)
pvalues <- future_lapply(1:nrow(count_norm), function(i) {
data <- cbind.data.frame(gene = as.numeric(t(count_norm[i, ])), group)
p <- suppressWarnings(wilcox.test(gene ~ group, data)$p.value)
return(p)
})
plan(sequential)
# 调整 p 值为 FDR
fdr <- p.adjust(unlist(pvalues), method = "BH")
# Bonferroni 校正
bonferroni <- p.adjust(unlist(pvalues), method = "bonferroni")
# 计算每个基因的对数表达水平的平均差（已经是 log2 转换数据）
log2FC <- rowMeans(count_norm[, group == "T"]) - rowMeans(count_norm[, group == "N"])
# 构建结果 DataFrame，包括 Bonferroni 校正后的 P 值
RankTest <- data.frame(logFC = log2FC, P.Value = unlist(pvalues), adj.P.Val = fdr, Bonferroni = bonferroni)

rownames(RankTest) <- rownames(count_norm)
# 保存结果到新数据文件
RankTest = na.omit(RankTest)
write.csv(RankTest, file = "1.差异分析/TPM/Wilcoxon 秩和检验结果.csv")
logFC_cutoff1 <- with(RankTest,mean(abs(logFC)) + 2*sd(abs(logFC)) )
if (手动设置阈值) {
logFC_cutoff1 <- cuf
}
message("Wilcox:动态 log2FC 阈值为：|logFC| > [mean(|logFC|) + 2sd(|logFC|)]")
message("本次分析采用的阈值为：")
print(logFC_cutoff1)
RTk1 <- (RankTest[[P 值选择]] < 0.05) & (RankTest$logFC < -logFC_cutoff1)
RTk2 <- (RankTest[[P 值选择]] < 0.05) & (RankTest$logFC > logFC_cutoff1)
RankTest$change = ifelse(RTk1,"DOWN",ifelse(RTk2,"UP","NOT"))
print(table(RankTest$change))
# 创建 GSEA 分析所需的数据
RankTest_log2FC_TPM = data.frame(gene_symbol = rownames(RankTest), change=RankTest$change , logFC
= RankTest$logFC)
write.csv(RankTest_log2FC_TPM, file = "1.差异分析/TPM/WilcoxTest_log2FC.csv", row.names = F)
limma_log2FC_fpkm_filtered <- RankTest_log2FC_TPM %>%
filter(change != "NOT")
write.csv(limma_log2FC_fpkm_filtered, file = "1.差异分析/TPM/Wilco_Test_上下调基因.csv", row.names
= FALSE)

小提琴图绘制
library(tidyverse)
library(ggplot2)
library(colourpicker)
# 读取 ssgsea_order 数据
ssgsea_file <- "path_to_ssgsea_order.csv" # 替换为你的文件路径
data <- read.csv(ssgsea_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
data <- t(data)
if ("P_value" %in% rownames(data)) {
data <- data[rownames(data) != "P_value", ]
}
data <- as.data.frame(data)
# 读取分组匹配信息
group_file <- "path_to_group_file.csv" # 替换为你的文件路径
new_dbssgsea <- read.csv(group_file, header = TRUE, row.names = 1,
stringsAsFactors = FALSE)
data$group <- new_dbssgsea[rownames(data), "Group"]
# 转换成长格式
cols_to_include <- setdiff(colnames(data), "group")
data_long <- pivot_longer(data, cols = all_of(cols_to_include), names_to = "Metric",
values_to = "Value")
# 选择统计检验方法
test_method <- "t.test" # 可以选择 "t.test" 或 "wilcox.test"
if (test_method == "t.test") {
pvalues <- data_long %>%
group_by(Metric) %>%
summarise(p = t.test(Value ~ group)$p.value)
} else if (test_method == "wilcox.test") {
pvalues <- data_long %>%
group_by(Metric) %>%
summarise(p = wilcox.test(Value ~ group)$p.value)
}
# p 值显示方式
p_display <- "numeric" # 可以选择 "numeric" 或 "stars"
if (p_display == "numeric") {
pvalues <- pvalues %>% mutate(label = paste0("p = ", signif(p, digits = 3)))
} else {
pvalues <- pvalues %>% mutate(label = case_when(
p < 0.001 ~ "***",
p < 0.01 ~ "**",
p < 0.05 ~ "*",
TRUE ~ "ns"
))
}
# 构造 facet 标签
new_labels <- setNames(paste0(pvalues$Metric, " (", pvalues$label, ")"),
pvalues$Metric)
# 自定义颜色
color1 <- "#e5451d"
color2 <- "#9084bd"
custom_colors <- c(color1, color2)
# 绘制小提琴图
violin_alpha <- 0.6 # 小提琴图透明度
scatter_alpha <- 0.6 # 散点图透明度
plot_points <- TRUE # 是否添加散点图层
p1 <- ggplot(data_long, aes(x = group, y = Value, fill = group)) +
geom_violin(trim = TRUE, alpha = violin_alpha, color = NA) +
{if(plot_points) geom_jitter(shape = 16, position = position_jitter(0.2), aes(color =
group),
alpha = scatter_alpha, size = 2) else NULL} +
geom_boxplot(width = 0.5, color = "black", fill = NA, outlier.shape = NA) +
scale_fill_manual(values = custom_colors) +
scale_color_manual(values = custom_colors) +
facet_wrap(~ Metric, ncol = 3, scales = "free_y", labeller = as_labeller(new_labels))
+
labs(title = " ", x = NULL, y = " ") +
theme_minimal(base_size = 15) +
theme(
plot.background = element_rect(fill = "white", color = NA),
panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
panel.background = element_rect(fill = "white", color = NA),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
plot.title = element_text(hjust = 0.5),
legend.position = "none",
axis.text.x = element_text(angle = 45, hjust = 1),
strip.background = element_rect(fill = "grey", color = "black", linewidth = 0.5),
strip.text = element_text(size = 10),
axis.ticks.y = element_line(color = "black")
)
# 保存为 PDF 文件
pdf_width <- 10
pdf_height <- 15
ggsave("score_violin_plot1.pdf",
plot = p1,
width = pdf_width,
height = pdf_height,
dpi = 300)
print("绘图完成，保存为 PDF 文件。")


#Figure 1B, S3
library(survival)
library(survminer)
library(ggplot2)
head(data)
#   event      time     value group
# 1     1 18.918033 1.7649406  High
# 2     1 30.983607 0.5583656   Low
# 3     0 93.803279 0.1974880   Low
# 4     1  1.016393 0.8807247   Low
# 5     1  7.967213 2.0750191  High
# 6     1 12.000000 0.5785040   Low
fit <- survfit(Surv(time, event) ~ group, data = data)
print(fit)
ggsurvplot(fit = fit, data = data, fun = "pct",
           palette = c("#4DBBD5", "#E64B35", "#BEBADA", "#FFFFB3", "#8DD3C7"),
           linetype = 1, pval = TRUE, 
           censor = TRUE, censor.size = 7,
           risk.table = FALSE, conf.int = FALSE)
#Figure 1C
# 加载必需包
library(survival)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(gridExtra)
# ------------------------------------------------------------------------------
# 1. 模拟TCGA-LAML队列数据（基于文章特征：n=150例，按IL4I1中位数分组）
# ------------------------------------------------------------------------------
set.seed(123) # 固定随机种子，确保结果可重复
n <- 150 # 文章中TCGA-LAML队列样本量
df <- data.frame(
  PatientID = paste0("P", 1:n),
  IL4I1_expr = rnorm(n, mean = 5, sd = 1.5), # 模拟IL4I1 mRNA表达量（正态分布）
  OS_time = sample(6:60, n, replace = TRUE), # 总生存时间（月）
  OS_status = sample(c(0,1), n, replace = TRUE, prob = c(0.6, 0.4)) # 生存状态（0=存活，1=死亡）
)

# 按IL4I1表达中位数分组（文章核心分组策略）
median_expr <- median(df$IL4I1_expr)
df <- df %>%
  mutate(
    IL4I1_group = ifelse(IL4I1_expr < median_expr, "Low", "High"), # 低/高表达组
    Risk_group = IL4I1_group, # 风险分组=IL4I1表达分组（文章中高表达为高风险）
    Risk_score = ifelse(Risk_group == "High", 1, 0) # 风险评分（高风险=1，低风险=0）
  ) %>%
  arrange(Risk_score, desc(IL4I1_expr)) # 按风险评分+表达量降序排序（优化可视化）

# 添加患者排序序号（用于x轴）
df$Patient_order <- 1:nrow(df)

# ------------------------------------------------------------------------------
# 2. 绘制预后风险因子图（三行整合：风险分组+生存状态+IL4I1表达）
# ------------------------------------------------------------------------------
# 定义颜色方案
color_risk <- c("Low" = "#2E86AB", "High" = "#E63946") # 低风险=蓝色，高风险=红色
color_status <- c("0" = "#2E86AB", "1" = "#E63946") # 存活=蓝色，死亡=红色
color_expr <- colorRampPalette(c("#2E86AB", "#E63946"))(100) # 表达量渐变色

# 子图1：风险分组（Low/High）
p1 <- ggplot(df, aes(x = Patient_order, y = 1, fill = Risk_group)) +
  geom_tile(height = 0.8) +
  scale_fill_manual(values = color_risk, name = "IL4I1 Expression") +
  labs(y = "Risk Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )

# 子图2：生存状态（0=存活，1=死亡）
p2 <- ggplot(df, aes(x = Patient_order, y = 1, fill = as.factor(OS_status))) +
  geom_tile(height = 0.8) +
  scale_fill_manual(values = color_status, name = "OS Status", labels = c("Alive", "Dead")) +
  labs(y = "Status") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )

# 子图3：IL4I1表达水平（热图形式）
p3 <- ggplot(df, aes(x = Patient_order, y = 1, fill = IL4I1_expr)) +
  geom_tile(height = 0.8) +
  scale_fill_gradientn(colors = color_expr, name = "IL4I1 Expression") +
  labs(x = "Patient (Ordered by Risk)", y = "IL4I1 Expr") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

# 组合三张子图，并添加统计信息标注（文章中HR和p值）
final_plot <- gridExtra::grid.arrange(
  p1, p2, p3,
  ncol = 1, nrow = 3,
  heights = c(0.8, 0.8, 1), # 调整子图高度比例
  top = textGrob(
    paste0("IL4I1 Prognostic Risk Plot (AML, TCGA-LAML)\nHR = 1.95 (95%CI: 1.26-3.02), p = 0.003"),
    gp = gpar(fontsize = 12, fontweight = "bold")
  )
)

# 保存图片（文章中常见图片格式）
ggsave(
  "IL4I1_Prognostic_Risk_Plot.png",
  final_plot,
  width = 12, height = 6, dpi = 300, device = "png"
)
# ------------------------------------------------------------------------------



#Figure 1D
# 安装依赖包（首次运行需执行）
# install.packages(c("ggplot2", "ggsignif"))

# 加载包
library(ggplot2)
library(ggsignif)

# ------------------------------------------------------------------------------
# 1. 模拟数据（匹配原图“生存/死亡组的IL4I1表达量”特征）
# ------------------------------------------------------------------------------
set.seed(456) # 固定随机种子，结果可重复
df <- data.frame(
  OS_event = rep(c("Alive", "Dead"), times = c(80, 40)), # 模拟样本量（Alive组例数更多）
  IL4I1_expr = c(
    rnorm(n = 80, mean = 1.2, sd = 0.5),  # Alive组：IL4I1低表达
    rnorm(n = 40, mean = 4.0, sd = 0.8)   # Dead组：IL4I1高表达（与原图一致）
  )
)

# ------------------------------------------------------------------------------
# 2. 绘制小提琴图（完全匹配原图）
# ------------------------------------------------------------------------------
p <- ggplot(df, aes(x = OS_event, y = IL4I1_expr, fill = OS_event)) +
  # 绘制小提琴图（匹配原图形状）
  geom_violin(
    width = 1.2, 
    alpha = 0.8, # 透明度与原图一致
    show.legend = FALSE # 隐藏图例（原图无图例）
  ) +
  # 添加组内统计量（中心点=中位数，线段=四分位距）
  stat_summary(
    fun = median, fun.min = IQR, fun.max = IQR,
    geom = "crossbar", width = 0.2, color = "black"
  ) +
  # 添加组间显著性标记（对应原图的“*”，代表p<0.05）
  geom_signif(
    comparisons = list(c("Alive", "Dead")),
    map_signif_level = TRUE, # 自动匹配“*”
    y_position = 5.2, # 星号位置与原图一致
    size = 0.6
  ) +
  # 颜色匹配原图（Alive=浅蓝，Dead=粉色）
  scale_fill_manual(values = c("Alive" = "#8ECAE6", "Dead" = "#FFB7B2")) +
  # 坐标轴标签与原图完全一致
  labs(
    x = "OS event",
    y = "The expression of IL4I1 (log2(TPM+1))"
  ) +
  # 调整主题（简洁化，匹配原图风格）
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(), # 隐藏x轴方向网格
    panel.grid.minor = element_blank(),   # 隐藏次要网格
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# 显示图形
print(p)

# 保存图片（与原图尺寸匹配）
ggsave(
  "IL4I1_OS_violin_plot.png",
  p,
  width = 5, height = 6, dpi = 300, device = "png"
)

#Figure 1E

setwd("./Cox")

# 加载包
library(survival)
library(ggplot2)

# 禁止转化为因子
options(stringsAsFactors = FALSE)

# 导入整合后的数据
data_combined <- read.table("单因素cox_输入文件.txt", header = TRUE, sep = "\t", row.names = 1)

# 设置临床特征的数据类型
data_combined$age <- as.numeric(data_combined$age)  # 连续变量
data_combined$ANAL <- as.numeric(data_combined$ANAL)  # 连续变量
data_combined$sex <- as.factor(data_combined$sex)    # 二分类变量
data_combined$TNM_T <- as.factor(data_combined$TNM_T)  # 分级变量
data_combined$TNM_N <- as.factor(data_combined$TNM_N)  # 分级变量
data_combined$grade <- as.factor(data_combined$grade)  # 病理分级

# 初始化结果存储
results_table <- data.frame()
significant_features <- c()

# 特征列表，包括临床特征和基因
category_list <- c("age", "sex", "TNM_T", "TNM_N", "grade", "ANAL")

# 进行单因素Cox回归分析
for (Category in category_list) {
  temp_data <- data_combined
  
  temp_data$Category <- if (Category %in% c("age", "ANAL")) {
    temp_data[[Category]]  # 连续变量
  } else {
    as.factor(temp_data[[Category]])  # 分类变量
  }
  
  cox_model <- coxph(Surv(futime, fustat) ~ Category, data = temp_data)
  cox_summary <- summary(cox_model)
  
  # 获取P值
  p_value <- cox_summary$coefficients[1, "Pr(>|z|)"]
  
  # 记录分析结果
  results_table <- rbind(results_table,
                         data.frame(
                           Category = Category,
                           HR = cox_summary$conf.int[1, "exp(coef)"],
                           HR_Lower95 = cox_summary$conf.int[1, "lower .95"],
                           HR_Upper95 = cox_summary$conf.int[1, "upper .95"],
                           p_value = p_value
                         )
  )
  
  # 标记显著特征
  if (p_value < 0.05) {
    significant_features <- c(significant_features, Category)
  }
}

# 输出单因素Cox回归结果到文件，可用在线网站作图
write.table(results_table, file = "Univariate_Cox_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# 绘制森林图
ggplot(results_table, aes(x = Category, y = HR, ymin = HR_Lower95, ymax = HR_Upper95)) +
  geom_errorbar(aes(ymin = HR_Lower95, ymax = HR_Upper95), width = 0.2, color = "darkblue", size = 0.9) +
  geom_point(shape = 15, color = "red", size = 4) +  # 绘制红色方形点
  coord_flip() +
  theme_classic() +  # 设置为经典主题
  theme(
    panel.grid.major = element_blank(),  # 删除主网格线
    panel.grid.minor = element_blank(),  # 删除次网格线
    panel.background = element_blank(),    # 删除背景填充
    plot.title = element_text(hjust = 0.5)  # 标题居中
  ) +
  xlab("Category") +
  ylab("Hazard Ratio (HR)") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ggtitle("Forest Plot of Univariate Cox Regression Analysis")

# 设置工作目录
setwd("E:/生信教程/27.多因素cox")  # 这里设置R的工作目录为存储数据和输出文件的路径
# 设置R语言的报错为英文
Sys.setenv(LANGUAGE = "en")
# 禁止转化为因子
options(stringsAsFactors = FALSE)
# 清空环境
rm(list=ls())

library("autoReg")  # 0.3.3 加载自动回归模型相关的包
library(survminer)  # 0.4.9 生存分析可视化包
library(survival)   # 3.5-7 生存分析的基本包
library(forestplot) # 3.1.3 生成森林图的包

mergedData <- read.table("合并的输入数据_2.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)  # 读取合并数据
risk <- read.table("risk.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)  # 导入风险文件
clinical <- read.table("clinical.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)  # 导入临床文件

# 提取风险和临床文件中相同的样本并合并数据
commonSamples <- intersect(row.names(clinical), row.names(risk))  # 找到相同样本
risk <- risk[commonSamples, ]  # 根据相同样本筛选风险数据
clinical <- clinical[commonSamples, ]  # 根据相同样本筛选临床数据
combinedData <- cbind(time = risk[,1], state = risk[,2], clinical)  # 合并时间、状态和临床数据
colnames(combinedData)  # 显示合并数据的列名以确认
mergedData <- combinedData

### 构建Cox回归模型
multiCoxModel <- coxph(Surv(time, state) ~ ., data = mergedData)  # 构建多因素Cox回归模型
summary(multiCoxModel)  # 查看模型的摘要结果
{
### 提取Cox回归分析的详细结果
multiCoxSummary <- summary(multiCoxModel)  # 获取多因素Cox回归模型的摘要
results <- data.frame(  # 构建结果数据框
  HR = multiCoxSummary$conf.int[,"exp(coef)"],  # 提取风险比
  L95CI = multiCoxSummary$conf.int[,"lower .95"],  # 提取95%置信区间的下限
  H95CI = multiCoxSummary$conf.int[,"upper .95"],  # 提取95%置信区间的上限
  pvalue = multiCoxSummary$coefficients[,"Pr(>|z|)"]  # 提取p值
)
results <- cbind(id = row.names(results), results)  # 添加id列
# 将结果输出为txt文件
write.table(results, file = "multiCoxClinical.txt", sep = "\t", row.names = FALSE, quote = FALSE)  


### 森林图绘制
# 读取数据
tdmcs <- read.table("multiCoxClinical.txt", header = TRUE, sep = "\t", row.names = 1)

# 获取数据
gene <- rownames(tdmcs)
hr <- sprintf("%.3f", tdmcs$HR)
hrLow <- sprintf("%.3f", tdmcs$L95CI)
hrHigh <- sprintf("%.3f", tdmcs$H95CI)
hazardRatio <- paste0(hr, " (", hrLow, "-", hrHigh, ")")
pValue <- ifelse(tdmcs$pvalue < 0.001, "<0.001", sprintf("%.3f", tdmcs$pvalue))

# 设置图像输出为 PDF
pdf(file = "Multivariate_Cox_Regression_Analysis.pdf", width = 6, height = 3)

# 获取行数
n <- nrow(tdmcs)
ylim <- c(1, n + 1)

# 创建布局
layout(matrix(c(1, 2), nc = 2), width = c(2.5, 2))

# 绘制文本部分 (基因名称, pvalue, Hazard ratio)
par(mar = c(4, 2.5, 2, 1))
plot(0, xlim = c(0, 3), ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
text.cex = 0.8
text(0, n:1, gene, adj = 0, cex = text.cex)  # 显示基因名称
text(1.2, n:1, pValue, adj = 1, cex = text.cex)  # 显示 p 值
text(3.0, n:1, hazardRatio, adj = 1, cex = text.cex)  # 显示 HR (Hazard Ratio)
text(1.2, n + 1, "pvalue", cex = text.cex, font = 2, adj = 1)
text(3.0, n + 1, "Hazard ratio", cex = text.cex, font = 2, adj = 1)

# 绘制区间图部分 (HR 置信区间)
par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
xlim <- c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
plot(0, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "Hazard ratio", ylab = "")
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
abline(v = 1, col = "black", lty = 2, lwd = 2)  # 添加参考线
boxcolor <- ifelse(as.numeric(hr) > 1, 'red', 'green')  # 设置颜色：红色表示风险增加，绿色表示降低
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.3)
axis(1)

# 关闭PDF
dev.off()


### 使用ggforest绘制森林图
p <- ggforest(multiCoxModel,  # 使用已构建的Cox模型
              main = "Hazard ratio",  # 图的主标题
              cpositions = c(0.02, 0.12, 0.25),  # 控制表格列的位置
              fontsize = 0.8,  # 设置字体大小
              refLabel = "reference",  # 设置参考标签
              noDigits = 2)  # 保留的小数位数

# 保存ggforest生成的森林图
ggsave("Multivariate_Cox_Forest_Plot.pdf", plot = p, width = 8, height = 6, units = "in")
}




# 加载包
library(circlize)
library(ComplexHeatmap)

#Figures 2B, 2C
# ------------------------------------------------------------------------------
# 1. 定义基因列表（与原图顺序完全一致）
# ------------------------------------------------------------------------------
genes <- c("IL4I1", "BCAT1", "BCAT2", "GOT2", "IDO1", "IDO2", "MAOA", "TAT", "MAOB", "NIT2", "TDO2")

# ------------------------------------------------------------------------------
# 2. 构造基因间相关性矩阵（模拟数据，替换为真实数据即可）
# ------------------------------------------------------------------------------
set.seed(789) # 固定随机种子，结果可重复
cor_mat <- matrix(
  runif(length(genes)^2, -1, 1), # 随机生成-1~1的相关系数
  nrow = length(genes), 
  ncol = length(genes),
  dimnames = list(genes, genes)
)
cor_mat <- (cor_mat + t(cor_mat))/2 # 保证矩阵对称（相关性是双向的）
diag(cor_mat) <- 1 # 基因自身的相关性为1

# ------------------------------------------------------------------------------
# 3. 设置颜色方案（匹配原图）
# ------------------------------------------------------------------------------
# 轨道颜色（每个基因对应的环形轨道颜色）
grid_col <- c(
  "IL4I1" = "#D3D3D3", # 浅灰
  "BCAT1" = "#AED6F1", # 浅蓝
  "BCAT2" = "#AED6F1", # 浅蓝
  "GOT2"  = "#F5B041", # 橙色
  "IDO1"  = "#AED6F1", # 浅蓝
  "IDO2"  = "#9B59B6", # 紫色
  "MAOA"  = "#E74C3C", # 红色
  "TAT"   = "#D3D3D3", # 浅灰
  "MAOB"  = "#82E0AA", # 浅绿
  "NIT2"  = "#AED6F1", # 浅蓝
  "TDO2"  = "#9B59B6"  # 紫色
)

# 弦颜色（相关性映射：-1→浅蓝，1→粉色）
col_fun <- colorRamp2(c(-1, 0, 1), c("#AED6F1", "white", "#FADBD8"))

# ------------------------------------------------------------------------------
# 4. 绘制弦图（完全还原原图）
# ------------------------------------------------------------------------------
circos.clear() # 清空之前的circlize配置

# 核心弦图绘制
chordDiagram(
  x = cor_mat,
  grid.col = grid_col, # 轨道颜色
  col = function(x) col_fun(cor_mat[x$row, x$col]), # 弦颜色由相关性决定
  transparency = 0.5, # 弦的透明度（匹配原图）
  directional = FALSE, # 相关性无方向
  preAllocateTracks = list(track.height = 0.15), # 轨道高度
  annotationTrack = "grid", # 仅显示轨道网格
  link.border = NA, # 隐藏弦的边框
  scale = FALSE # 不缩放相关性（保持原始范围）
)

# 添加基因名称标签
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1] - 0.1,
      sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.8
    )
  }
)

# 添加刻度标签（匹配原图的0/0.4/0.8/1.2）
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    breaks <- seq(0, 1.2, by = 0.4)
    circos.axis(
      h = "top",
      major.at = breaks,
      labels = as.character(breaks),
      cex = 0.6,
      labels.cex = 0.6
    )
  }
)

# 添加相关性图例（与原图一致）
draw(
  Legend(
    col_fun = col_fun,
    title = "Correlation",
    title_position = "right",
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    direction = "horizontal",
    legend_width = unit(4, "cm")
  ),
  x = unit(0.5, "npc"), 
  y = unit(-0.1, "npc"), 
  just = c("center", "top")
)

# 保存图片
png("Gene_Correlation_Chord_Plot.png", width = 8, height = 8, units = "in", res = 300)
# 重新运行上述绘图代码（避免交互式绘图的显示问题）
circos.clear()
chordDiagram(...) # 重复上述chordDiagram及后续代码
dev.off()


###Figure 2C


# 加载必要的包
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# 设置颜色方案
color_scale <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# 数据准备
# 假设expr_data是一个数据框，行是基因，列是样本
# 基因包括IL4I1和其他20个相关基因

# 1. 提取IL4I1和20个相关基因的表达数据
gene_list <- c("IL4I1", "CD300C", "CDH23", "CX3CR1", "DOK2", "GPR132", 
               "FGR", "FGD2", "IL1RN", "MVP", "LRRC25", "LILRB2", 
               "SCIMP", "TNFRSF1B", "TFEB", "SLC15A3", "VDR", 
               "LILRB1", "UNC119", "AC064805.1", "VSIR")

# 提取表达矩阵
expr_subset <- expr_data[rownames(expr_data) %in% gene_list, ]

# 2. 计算相关性矩阵
cor_matrix <- cor(t(expr_subset), method = "spearman")

# 3. 提取IL4I1与其他基因的相关性
il4i1_cor <- cor_matrix["IL4I1", ]
il4i1_cor <- il4i1_cor[names(il4i1_cor) != "IL4I1"]  # 移除IL4I1自身

# 按相关性从高到低排序
il4i1_cor_sorted <- sort(il4i1_cor, decreasing = TRUE)

# 4. 创建热图数据
heatmap_data <- matrix(il4i1_cor_sorted, nrow = 1)
rownames(heatmap_data) <- "IL4I1"
colnames(heatmap_data) <- names(il4i1_cor_sorted)

# 5. 创建显著性标注
# 计算p值
p_values <- sapply(colnames(heatmap_data), function(gene) {
  cor_test <- cor.test(expr_subset["IL4I1", ], expr_subset[gene, ], 
                       method = "spearman")
  cor_test$p.value
})

# 创建显著性标注
sig_labels <- ifelse(p_values < 0.001, "***", 
                     ifelse(p_values < 0.01, "**",
                            ifelse(p_values < 0.05, "*", "")))

# 6. 创建热图
Heatmap(heatmap_data,
        name = "Spearman\nCorrelation",
        col = color_scale,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        
        # 添加相关性值
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", heatmap_data[i, j]), 
                   x, y, gp = gpar(fontsize = 8))
          # 添加显著性标记
          grid.text(sig_labels[j], 
                   x, y - 0.2, gp = gpar(fontsize = 10))
        },
        
        # 添加标题
        column_title = "Top 20 Genes Positively Correlated with IL4I1 Expression",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        
        # 添加图例
        heatmap_legend_param = list(
          title = "Correlation",
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          at = c(-1, -0.5, 0, 0.5, 1),
          labels = c("-1.0", "-0.5", "0", "0.5", "1.0")
        ),
        
        # 设置热图高度
        height = unit(1.5, "cm"))
		
		
		
		
##Figure 3A 

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggnewscale)  
library(ggrepel)

# 读取差异分析结果文件
results <- read.table("差异分析结果.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 添加表达状态列
results$expression <- 'ns'
results$expression[results$LogFC >= 1 & results$P.adj < 0.05] <- 'up'
results$expression[results$LogFC <= -1 & results$P.adj < 0.05] <- 'down'

# 绘制火山图，下方的"Control vs disease"改成自己的。颜色也可以改成自己需要的
{
  p1 <- ggplot(results, aes(x = LogFC, y = -log10(P.adj))) +
    geom_point(aes(color = expression), size = 1.5) +
    scale_color_manual(values = c('up' = 'red', 'down' = 'blue', 'ns' = 'grey')) +
    labs(title = "Control vs disease", x = "log2(fold change)", y = "-log10(P.adj)") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(color = "black")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
}
# 打印并保存火山图
print(p1)
ggsave("volcano_plot.pdf", plot = p1, width = 9, height = 8)
dev.off()


# 找到logFC最大和最小的5个基因
# 找到logFC最大和最小的5个基因
# 找到logFC最大和最小的5个基因
# 找到logFC最大和最小的5个基因
top_upregulated <- results %>% 
  filter(expression == 'up') %>% 
  arrange(desc(LogFC)) %>% 
  head(5)

top_downregulated <- results %>% 
  filter(expression == 'down') %>% 
  arrange(LogFC) %>% 
  head(5)

# 合并两个数据框以便标注
top_genes <- rbind(top_upregulated, top_downregulated)

# 绘制火山图并标注基因
{
  p2 <- ggplot(results, aes(x = LogFC, y = -log10(P.adj))) +
    geom_point(aes(color = expression), size = 1.5) +
    scale_color_manual(values = c('up' = 'red', 'down' = 'blue', 'ns' = 'grey')) +
    geom_text_repel(data = top_genes, aes(label = gene), 
                    vjust = -0.5, size = 4, fontface = "bold") +
    labs(title = "Control vs disease", x = "log2(fold change)", y = "-log10(P.adj)") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(color = "black")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
}
# 打印并保存火山图
print(p2)

ggsave("top5基因-volcano_plot.pdf", plot = p2, width = 9, height = 8)
dev.off()


# 自定义基因标记
# 自定义基因标记
# 自定义基因标记
# "Numb", "Snapc2", "Abhd16a", "Gm3244"替换为您想标记的基因名
desired_genes <- c("Numb", "Snapc2", "Abhd16a", "Gm3244")  

# 从results中提取这些基因的数据
custom_genes <- results[results$gene %in% desired_genes, ]
{
# 绘制火山图并标注自定义基因
p3 <- ggplot(results, aes(x = LogFC, y = -log10(P.adj))) +
  geom_point(aes(color = expression), size = 1.5) +
  scale_color_manual(values = c('up' = 'red', 'down' = 'blue', 'ns' = 'grey')) +
  geom_text_repel(data = custom_genes, aes(label = gene), 
                  vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "Control vs disease", x = "log2(fold change)", y = "-log10(P.adj)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
}
# 打印并保存火山图
print(p3)
ggsave("自定义基因-volcano_plot.pdf", plot = p3, width = 9, height = 8)
dev.off()



# 自定义呈现通路基因
# 自定义呈现通路基因
# 自定义呈现通路基因
# 自定义呈现通路基因

# 处理GO数据
# 处理GO数据
# 处理GO数据
GO <- read.table("GO.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

matric <- GO %>%
  separate_rows(geneID, sep = "/") %>%
  select(Description, geneID)

# 查看生成的矩阵
print(matric)


# 选择要标记的通路
# 选择要标记的通路
# 选择要标记的通路
# 选择要标记的通路
pathway_genes <- matric %>% 
  filter(Description == "external encapsulating structure organization") %>%  # 替换为目标通路
  left_join(results, by = c("geneID" = "gene"))  # 合并数据

# 颜色定义 -------------------------------------------------
pathway_colors1 <- c("external encapsulating structure organization" = "#00BFC4")  # 定义通路颜色
{
# 绘图核心代码 ----------------------------------------------
p4 <- ggplot(results, aes(x = LogFC, y = -log10(P.adj))) +
  # 第一部分：基础火山点
  geom_point(aes(color = expression), size = 1.5) +
  scale_color_manual(
    name = "Expression",
    values = c('up' = 'red', 'down' = 'blue', 'ns' = 'grey'),
    guide = guide_legend(order = 1)
  ) +
  
  # 第二部分：通路基因标注
  new_scale_color() +
  geom_text_repel(
    data = pathway_genes,
    aes(label = geneID, color = Description),  # 颜色映射到通路
    size = 4,
    fontface = "bold",
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    name = "Pathway",
    values = pathway_colors1,
    guide = guide_legend(
      order = 2,
      override.aes = list(  # 关键修改点
        shape = 15,         # 使用方块形状
        size = 3,           # 调整图例符号大小
        label = ""          # 隐藏默认文本标签
      )
    )
  ) +
  
  # 通用图形设置
  labs(title = "Control vs Disease", 
       x = "log2(fold change)", 
       y = "-log10(P.adj)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),  # 调整图注间距
    legend.key = element_rect(fill = "white")  # 清除图例背景
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
}
# 输出结果
print(p4)
ggsave("单个自定义通路火山图1.pdf", plot = p4, width = 10, height = 8)
dev.off()



##Figure 3B, C, D  

if (!dir.exists("GSEA")) {
dir.create("GSEA")
}
Genes_All <- read.csv(输入数据的文件名称, header = TRUE, check.names = FALSE,
row.names = NULL)
colnames(Genes_All)[1:3] <- c("gene_symbol", "change", "logFC")
#gene ID 转换
gene <- suppressWarnings(suppressMessages(bitr(Genes_All$gene_symbol,
fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_物种缩写)))
#GSEA 所需数据
info <- Genes_All[, c("gene_symbol", "logFC")]
names(info) <- c('SYMBOL','Log2FoldChange')
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因 ID 和
Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
cat("GESA-KEGG 分析开始······\n")
GSEA_input = sort(GSEA_input, decreasing = TRUE)
#KEGG-GSEA 分析
KEGG_ges <- suppressWarnings(gseKEGG(GSEA_input, organism = KEGG_物种缩
写, pvalueCutoff = 1,
eps = 0))
#将 ENTREZID 重转为 symbol：
KEGG_ges <- setReadable(KEGG_ges,
OrgDb = GO_物种缩写,
keyType = "ENTREZID")
#提取结果表并将富集结果保存到本地
KEGG_ges_result <- KEGG_ges@result
write.csv(KEGG_ges_result, file = "GSEA/GSEA-KEGG.csv")
p <- ridgeplot(KEGG_ges,
showCategory = as.numeric(KEGG 展示通路个数),
fill = P 值展示,
decreasing = TRUE) +
theme(axis.text.y = element_text(size = 8), # 调整 y 轴字体大小为 8
axis.text.x = element_text(size = 8)) # 调整 x 轴字体大小为 8
p
ggsave("GSEA/1.GSEA-KEGG 山峦图.pdf", plot = p, width = 15, height = 15, units =
"cm")
message("\nGSEA-KEGG 山峦图图绘制完成，并保存在文件夹中")
pp <- dotplot(KEGG_ges,showCategory = KEGG 展示通路个数, font.size=8, color =
P 值展示)
pp
ggsave("GSEA/2.GSEA-KEGG 气泡图.pdf", plot = pp, width = 15, height = 15, units
= "cm")
message("\nGSEA-KEGG 气泡图绘制完成，并保存在文件夹中")
p5 <- gseaplot2(KEGG_ges,
geneSetID = geneSetID, #或直接输入基因集 ID 向量名，如
c("hsa04610","hsa00260")
color = color_GSEA,
pvalue_table = TRUE,
ES_geom = "line")
ggsave("GSEA/3.GSEA-KEGG 主图.pdf", plot = p5, width = 30, height = 30, units =
"cm")
message("\nGSEA-KEGG 主图绘制完成，并保存在文件夹中")
#KEGG-GO 分析
cat("\nGSEA-GO 分析开始······\n")
GO_ges <- suppressWarnings(gseGO(geneList = GSEA_input,
OrgDb = GO_物种缩写,
ont = "ALL",
pvalueCutoff = 1,
eps = 0
))
GO_ges <- setReadable(GO_ges, OrgDb = GO_物种缩写, keyType = "ENTREZID")
GO_ges_result <- GO_ges@result
write.csv(GO_ges_result, file = 'GSEA/GSEA-GO.csv')
p <- ridgeplot(GO_ges,
showCategory = (GO 各展示个数*3),
fill = P 值展示,
decreasing = TRUE) +
theme(axis.text.y = element_text(size = 8), # 调整 y 轴字体大小为 8
axis.text.x = element_text(size = 8)) # 调整 x 轴字体大小为 8
p
ggsave("GSEA/4.GSEA-GO 山峦图.pdf", plot = p, width = 15, height = 15, units =
"cm")
message("\nGSEA-GO 山峦图绘制完成，并保存在文件夹中")
pp <- dotplot(GO_ges, showCategory = GO 各展示个数, split="ONTOLOGY", color
= P 值展示,font.size=8) + facet_grid(ONTOLOGY~., scale="free")
pp
ggsave("GSEA/5.GSEA-GO 气泡图.pdf", plot = pp, width = 15, height = 15, units =
"cm")
message("\nGSEA-GO 气泡图绘制完成，并保存在文件夹中")
p5 <- gseaplot2(GO_ges,
geneSetID = geneSetID, #或直接输入基因集 ID 向量名，如
c("hsa04610","hsa00260")
color = color_GSEA,
pvalue_table = TRUE,
ES_geom = "line")
ggsave("GSEA/6.GSEA-GO 主图.pdf", plot = p5, width = 30, height = 30, units =
"cm")
message("\nGSEA-GO 主图绘制完成，并保存在文件夹中")


if (!dir.exists("UP_Gene")) {
dir.create("UP_Gene")
}
Genes_All <- read.csv(输入数据的文件名称, header = TRUE, check.names = FALSE,
row.names = NULL)
if (ncol(Genes_All) == 3) {
colnames(Genes_All)[1:3] <- c("gene_symbol", "change", "logFC")
} else if (ncol(Genes_All) == 2) {
colnames(Genes_All)[1:2] <- c("gene_symbol", "change")
}
Genes_Filtered <- Genes_All[Genes_All$change != "NOT", ]
Genes_Filtered <- Genes_Filtered[Genes_Filtered$change != "DOWN", ]
#gene ID 转换
gene <- suppressWarnings(suppressMessages(bitr(Genes_Filtered$gene_symbol,
fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_物种缩写)))
R.utils::setOption("clusterProfiler.download.method",'auto')
cat("开始 KEGG-GO 上调基因的富集分析······")
#GO 分析
GO<-enrichGO( gene$ENTREZID,#GO 富集分析
OrgDb = GO_物种缩写,
keyType = "ENTREZID",#设定读取的 gene ID 类型
ont = "ALL",#(ont 为 ALL 因此包括 Biological Process,Cellular
Component,Mollecular Function 三部分）
pvalueCutoff = 1,#设定 p 值阈值
qvalueCutoff = 1,#设定 q 值阈值
readable = T)
# 输出提示信息
cat("\nGO 富集分析已完成")
GO
write.csv(GO@result, file = "UP_Gene/GO_result.csv", row.names = FALSE)
#KEGG 分析
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG 富集分析
organism = KEGG_物种缩写,
pvalueCutoff = 1,
qvalueCutoff = 1)
# 输出提示信息
cat("\nKEGG 富集分析已完成")
KEGG
KEGG <- setReadable(KEGG, OrgDb=GO_物种缩写, 'ENTREZID')
write.csv(KEGG@result, file = "UP_Gene/KEGG_result.csv", row.names = FALSE)
#删除后缀
KEGG@result$Description <- gsub(" - Mus musculus \\(house mouse\\)$", "",
KEGG@result$Description)
#保存柱状图
p <- barplot(GO, showCategory = GO 各展示个数, split="ONTOLOGY", color = P 值
展示,font.size=8) + facet_grid(ONTOLOGY~., scale="free")
ggsave("UP_Gene/1.GO 柱状图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("\n1.GO 柱状图已绘制完成，并保存在文件夹中")
p <- barplot(KEGG,showCategory = KEGG 展示通路个数, font.size=8, color = P 值
展示,label_format = 40,)
ggsave("UP_Gene/2.KEGG 柱状图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("2.KEGG 柱状图已绘制完成，并保存在文件夹中")
#点状图
p <- dotplot(GO, showCategory = GO 各展示个数, split="ONTOLOGY", color = P 值
展示,font.size=8) + facet_grid(ONTOLOGY~., scale="free")
ggsave("UP_Gene/3.GO 气泡图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("3.GO 气泡图已绘制完成，并保存在文件夹中")
p <- dotplot(KEGG,showCategory = KEGG 展示通路个数, font.size=8, color = P 值
展示)
ggsave("UP_Gene/4.KEGG 气泡图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("4.KEGG 气泡图已绘制完成，并保存在文件夹中")


if (!dir.exists("DOWN_Gene")) {
dir.create("DOWN_Gene")
}
Genes_All <- read.csv(输入数据的文件名称, header = TRUE, check.names = FALSE,
row.names = NULL)
if (ncol(Genes_All) == 3) {
colnames(Genes_All)[1:3] <- c("gene_symbol", "change", "logFC")
} else if (ncol(Genes_All) == 2) {
colnames(Genes_All)[1:2] <- c("gene_symbol", "change")
}
Genes_Filtered <- Genes_All[Genes_All$change != "NOT", ]
Genes_Filtered <- Genes_Filtered[Genes_Filtered$change != "UP", ]
#gene ID 转换
gene <- suppressWarnings(suppressMessages(bitr(Genes_Filtered$gene_symbol,
fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_物种缩写)))
R.utils::setOption("clusterProfiler.download.method",'auto')
cat("开始 KEGG-GO 下调基因的富集分析······")
#GO 分析
GO<-enrichGO( gene$ENTREZID,#GO 富集分析
OrgDb = GO_物种缩写,
keyType = "ENTREZID",#设定读取的 gene ID 类型
ont = "ALL",#(ont 为 ALL 因此包括 Biological Process,Cellular
Component,Mollecular Function 三部分）
pvalueCutoff = 1,#设定 p 值阈值
qvalueCutoff = 1,#设定 q 值阈值
readable = T)
# 输出提示信息
cat("\nGO 富集分析已完成")
GO
write.csv(GO@result, file = "DOWN_Gene/GO_result.csv", row.names = FALSE)
#KEGG 分析
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG 富集分析
organism = KEGG_物种缩写,
pvalueCutoff = 1,
qvalueCutoff = 1)
# 输出提示信息
cat("\nKEGG 富集分析已完成")
KEGG
KEGG <- setReadable(KEGG, OrgDb=GO_物种缩写, 'ENTREZID')
write.csv(KEGG@result, file = "DOWN_Gene/KEGG_result.csv", row.names =
FALSE)
#删除后缀
KEGG@result$Description <- gsub(" - Mus musculus \\(house mouse\\)$", "",
KEGG@result$Description)
#保存柱状图
p <- barplot(GO, showCategory = GO 各展示个数, split="ONTOLOGY", color = P 值
展示,font.size=8) + facet_grid(ONTOLOGY~., scale="free")
ggsave("DOWN_Gene/1.GO 柱状图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("\n1.GO 柱状图已绘制完成，并保存在文件夹中")
p <- barplot(KEGG,showCategory = KEGG 展示通路个数, font.size=8, color = P 值
展示,label_format = 40,)
ggsave("DOWN_Gene/2.KEGG 柱状图.pdf", plot = p, width = 15, height = 16, units
= "cm") # 设置长宽为 10cm x 8cm
message("2.KEGG 柱状图已绘制完成，并保存在文件夹中")
#点状图
p <- dotplot(GO, showCategory = GO 各展示个数, split="ONTOLOGY", color = P 值
展示,font.size=8) + facet_grid(ONTOLOGY~., scale="free")
ggsave("DOWN_Gene/3.GO 气泡图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("3.GO 气泡图已绘制完成，并保存在文件夹中")
p <- dotplot(KEGG,showCategory = KEGG 展示通路个数, font.size=8, color = P 值
展示)
ggsave("DOWN_Gene/4.KEGG 气泡图.pdf", plot = p, width = 15, height = 16, units
= "cm") # 设置长宽为 10cm x 8cm
message("4.KEGG 气泡图已绘制完成，并保存在文件夹中")


#Figure 3E, F
if (!dir.exists("All_Gene")) {
dir.create("All_Gene")
}
Genes_All <- read.csv(输入数据的文件名称, header = TRUE, check.names = FALSE,
row.names = NULL)
if (ncol(Genes_All) == 3) {
colnames(Genes_All)[1:3] <- c("gene_symbol", "change", "logFC")
} else if (ncol(Genes_All) == 2) {
colnames(Genes_All)[1:2] <- c("gene_symbol", "change")
}
Genes_Filtered <- Genes_All[Genes_All$change != "NOT", ]
#gene ID 转换

gene <- suppressWarnings(suppressMessages(bitr(Genes_Filtered$gene_symbol,
fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_物种缩写)))
R.utils::setOption("clusterProfiler.download.method",'auto')
cat("开始 KEGG-GO 全部差异基因富集分析······")
#GO 分析
GO<-enrichGO( gene$ENTREZID,#GO 富集分析
OrgDb = GO_物种缩写,
keyType = "ENTREZID",#设定读取的 gene ID 类型
ont = "ALL",#(ont 为 ALL 因此包括 Biological Process,Cellular
Component,Mollecular Function 三部分）
pvalueCutoff = 1,#设定 p 值阈值
qvalueCutoff = 1,#设定 q 值阈值
readable = T)
# 输出提示信息
cat("\nGO 富集分析已完成\n")
GO
write.csv(GO@result, file = "All_Gene/GO_result.csv", row.names = FALSE)
#KEGG 分析
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG 富集分析
organism = KEGG_物种缩写,
pvalueCutoff = 1,
qvalueCutoff = 1)
# 输出提示信息
cat("\nKEGG 富集分析已完成")
KEGG
KEGG <- setReadable(KEGG, OrgDb=GO_物种缩写, 'ENTREZID')
write.csv(KEGG@result, file = "All_Gene/KEGG_result.csv", row.names = FALSE)
#删除后缀
KEGG@result$Description <- gsub(" - Mus musculus \\(house mouse\\)$", "",
KEGG@result$Description)
木火医学v:cgxr410木火医学：生信分析实战教材
（ VX ： crxr410 ）
#保存柱状图
p <- barplot(GO, showCategory = GO 各展示个数, split="ONTOLOGY", color = P 值
展示,font.size=8) + facet_grid(ONTOLOGY~., scale="free")
ggsave("All_Gene/1.GO 柱状图.pdf", plot = p, width = 15, height = 16, units = "cm")
# 设置长宽为 10cm x 8cm
message("\n1.GO 柱状图已绘制完成，并保存在文件夹中")
p <- barplot(KEGG,showCategory = KEGG 展示通路个数, font.size=8, color = P 值
展示,label_format = 40,)
ggsave("All_Gene/2.KEGG 柱状图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("2.KEGG 柱状图已绘制完成，并保存在文件夹中")
#点状图
p <- dotplot(GO, showCategory = GO 各展示个数, split="ONTOLOGY", color = P 值
展示,font.size=8) + facet_grid(ONTOLOGY~., scale="free")
ggsave("All_Gene/3.GO 气泡图.pdf", plot = p, width = 15, height = 16, units = "cm")
# 设置长宽为 10cm x 8cm
message("3.GO 气泡图已绘制完成，并保存在文件夹中")
p <- dotplot(KEGG,showCategory = KEGG 展示通路个数, font.size=8, color = P 值
展示)
ggsave("All_Gene/4.KEGG 气泡图.pdf", plot = p, width = 15, height = 16, units =
"cm") # 设置长宽为 10cm x 8cm
message("4.KEGG 气泡图已绘制完成，并保存在文件夹中")





#Figure 4A, B 
if (!dir.exists("ssgsea")) {
dir.create("ssgsea")
}
signature = signatureAPP
expr = exprim
# 基于 ssGSEA 方法的单样本分析
ssgsea <- calculate_sig_score(eset = expr,
method = "ssgsea",
signature = signature,
mini_gene_count = mini_gene_count)
write.csv(ssgsea, file = "ssgsea/ssgsea_score.csv", row.names = FALSE)
cite = signature_collection_citation
write.csv(cite, file = "ssgsea 所用基因集的参考文献.csv", row.names = FALSE)
# 排序展示
data_DEG <- read.csv("ssgsea/ssgsea_score.csv", header = TRUE, row.names = 1, check.names
= FALSE)
data_DEG$Group <- new_db[rownames(data_DEG), "Group"]
data_cols <- setdiff(names(data_DEG), "Group") # 提取除了 Group 列之外的所有列名
# 准备一个空的数据框来存储 P 值结果
p_values <- data.frame(Column = character(), P_value = numeric(), stringsAsFactors = FALSE)
if (method.fig1 == "wilcox.test") {
# 执行 Wilcoxon 秩和检验
for (col in data_cols) {
wilcox_test_result <- wilcox.test(data_DEG[[col]] ~ data_DEG$Group)
p_value <- wilcox_test_result$p.value
p_values <- rbind(p_values, data.frame(Column = col, P_value = p_value))
}
} else if (method.fig1 == "t.test") {
# 执行 t 检验
for (col in data_cols) {
t_test_result <- t.test(data_DEG[[col]] ~ data_DEG$Group)
p_value <- t_test_result$p.value
p_values <- rbind(p_values, data.frame(Column = col, P_value = p_value))
}
} else if (method.fig1 == "anova") {
# 执行方差分析
for (col in data_cols) {
anova_result <- aov(data_DEG[[col]] ~ data_DEG$Group)
p_value <- summary(anova_result)[[1]]$"Pr(>F)"[1]
p_values <- rbind(p_values, data.frame(Column = col, P_value = p_value))
}
}
# 去除 data_DEG 中的 Group 列
data_DEG$Group <- NULL
# 将 P 值添加到 data_DEG 中的一行
data_DEG <- rbind(data_DEG, p_values$P_value)
# 可选：给新行命名
rownames(data_DEG)[nrow(data_DEG)] <- "P_value"
data_DEG = t(data_DEG)
data_DEG <- as.data.frame(data_DEG)
data_DEG <- data_DEG[order(data_DEG$P_value), ]
write.csv(data_DEG, file = "ssgsea/ssgsea_order.csv", row.names = TRUE)
3. ssGSEApolt
library(ggplot2)
library(dplyr)
library(cowplot)
library(randomcoloR) # 使用 randomcoloR 生成随机颜色
library(colourpicker) # 使用 colourpicker 设置颜色
library(reshape2) # 使用 reshape2 包进行数据转换
library(pheatmap) # 绘制热图
# 读取上传的 CSV 文件
ssgsea_data <- read.csv("你的文件路径.csv", header = TRUE, row.names = 1, check.names =
FALSE)
# 动态获取用户输入的目标基因
gene.duo <- c("TP53", "DPM1", "ABCEFD") # 用户输入的目标基因名称
expr <- exprim
cluster <- cluster1
# 数据处理
data <- ssgsea_data
data$P_value <- NULL
data = as.data.frame(t(data))
data2 = data
data$Group <- new_db[rownames(data), "Group"]
data1 = data
# 转换数据格式并添加分组信息
df <- reshape2::melt(data1, id.vars = "Group", variable.name = "Cell", value.name = "score")
# 绘制箱线图
p <- ggplot(data=df, aes(x=Cell, y=score, fill=Group)) +
geom_boxplot(width=0.3, position = position_dodge(0.5), outlier.colour = NA) +
theme_bw() +
theme(panel.grid = element_blank()) +
scale_fill_manual(values = values) +
stat_compare_means(aes(group=Group), method = method.fig1, label="p.signif", hide.ns =
TRUE, label.y = c(label.y)) +
theme(axis.text.x = element_text(size=size.fig1, colour=colour.fig1, angle=angle.fig1x, hjust
= 1),
axis.text.y = element_text(size=size.fig1, angle=0),
axis.title.x = element_text(size=size.fig1),
axis.title.y = element_text(size=size.fig1),
legend.text = element_text(size=size.fig1)) + # 控制图例中文本字体大小
xlab("") # 去掉 x 轴标签
# 打印图形并保存为 PDF 文件
print(p)
ggsave("ssgsea/1.ssgsea_boxplot.pdf", plot = p, width = width.fig1, height = height.fig1)
# 绘制堆叠图
im_timer1 = data2
total_scores <- rowSums(im_timer1)
proportion_scores <- sweep(im_timer1, 1, total_scores, "/")
timer1 = proportion_scores
timer1 <- timer1 %>%
rownames_to_column(var = "Sample")
Cellratio <- timer1 %>%
pivot_longer(cols = -Sample, names_to = "CellName", values_to = "Proportion")
P1 <- ggplot(Cellratio) +
geom_bar(aes(x = Sample, y = Proportion, fill = CellName), stat = "identity", width =
width.figt2, linewidth = linewidth.fig2, colour = colour.fig2) +
theme_classic() +
labs(x = 'Sample', y = 'Ratio') +
coord_flip() +
theme(panel.border = element_rect(fill = NA, color = color.fig2, linewidth = linewidth.figw2,
linetype = "solid"))
# 保存堆叠图为 PDF 文件
pdf(file="ssgsea/2.ssgsea_stacke.pdf", width=width.fig2, height=height.fig2)
print(P1)
dev.off()
write.csv(proportion_scores, file = "ssgsea/proportion_scores.csv", row.names = T)
# 基因相关性分析
gene = gene.duo
gene_exp = expr[gene,]
gene_exp <- gene_exp[complete.cases(gene_exp), ]
gene_exp = t(gene_exp)
# 确保两个数据框拥有相同行名
common_rows <- intersect(rownames(gene_exp), rownames(im_timer1))
# 子集化数据框，仅保留相同行名的行
gene_exp_subset <- gene_exp[common_rows, ]
im_timer1_subset <- im_timer1[common_rows, ]
# 计算两个数据框所有列之间的相关性
correlation_matrix <- cor(gene_exp_subset, im_timer1_subset, method=method.fig3)
na_cols <- apply(is.na(correlation_matrix), 2, any)
correlation_matrix <- correlation_matrix[, !na_cols]
write.csv(correlation_matrix, file = "ssgsea/相关性数值.csv", row.names = T)
datacor = t(correlation_matrix)
# 绘制相关性热图
p2 <- pheatmap(datacor, scale = "none", cluster_cols = cluster, cluster_rows = cluster,
color = colorRampPalette(color.fig3)(50),
legend = T, show_colnames = T, show_rownames = T,
fontsize = fontsize, cellheight = cellheight, cellwidth = cellwidth,
heatmap_legend_param = list(title = ""))
# 保存为 PDF 文件
pdf(file="ssgsea/3.多个基因免疫相关性图.pdf", width=width.fig3, height=height.fig3)
print(p2)
dev.off()
# 循环遍历每个基因，生成并保存气泡图
for (siggene in rownames(correlation_matrix)) {
# 提取特定基因的相关性数据
siggene_cor <- as.data.frame(t(correlation_matrix[siggene, , drop = FALSE]))
siggene_cor$immune <- rownames(siggene_cor)
# 提取相关性值
cor_values <- siggene_cor[, 1]
# 提取第一列的列名作为标题
title_name <- names(siggene_cor)[1]
# 绘制气泡图
p3 <- ggplot(data = siggene_cor, aes(x = cor_values, y = immune)) +
labs(x = "Correlation coefficient", y = "Immune cell") +
geom_point(aes(size = abs(cor_values), color = cor_values)) +
scale_color_gradient2(low = low.fig4, mid = mid.fig4, high = high.fig4, midpoint = 0) +
theme(panel.background = element_rect(fill = "white", size = 1, color = "black"),
panel.grid = element_line(color = "grey75", size = 0.5),
axis.ticks = element_line(size = 0.5),
axis.text.y = element_text(colour = "black", size = 9),
plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
xlab("Correlation coefficient") +
ylab("Immune cell") +
ggtitle(title_name)
# 保存为 PDF 文件，文件名包含基因名
pdf(file = paste0("ssgsea/", siggene, "_immune_correlation.pdf"), width = width.fig4, height
= height.fig4)
print(p3)
dev.off()



##Figure 4C

if (!dir.exists("im_estimate")) {
dir.create("im_estimate")
}
expr=exprim
cluster=cluster1
im_estimate <- deconvo_tme(eset = expr,
method = "estimate",
platform = platform)
write.csv(im_estimate, file = "im_estimate/im_estimate.csv", row.names = FALSE)
# 去掉多余的列
data <- im_estimate
data <- as.data.frame(data)
data$Group <- group_labels
data1 <- data[ , -1]
# 转换数据格式并添加分组信息
df <- reshape2::melt(data1, id.vars = "Group", variable.name = "Cell", value.name =
"score")
p <- ggplot(data=df, aes(x=Cell, y=score, fill=Group)) +
geom_boxplot(width=0.3, position = position_dodge(0.5), outlier.colour = NA) +
theme_bw() +
theme(panel.grid = element_blank()) +
scale_fill_manual(values = values) +
stat_compare_means(aes(group=Group), method = method.fig1, label="p.signif",
hide.ns = TRUE, label.y = c(label.y)) +
theme(axis.text.x = element_text(size=size.fig1, colour=colour.fig1, angle=angle.fig1x,
hjust = 1),
axis.text.y = element_text(size=size.fig1, angle=0),
axis.title.x = element_text(size=size.fig1),
axis.title.y = element_text(size=size.fig1),
legend.text = element_text(size=size.fig1)) + # 控制图例中文本字体大小
xlab("") # 去掉 x 轴标签
# 打印图形
print(p)
# 保存为 PDF 文件
ggsave("im_estimate/1.immune_boxplot.pdf", plot = p, width = width.fig1, height =
height.fig1)
# 堆叠图的绘制
im_timer1=im_estimate
im_timer1 <- column_to_rownames(im_timer1, var = colnames(im_timer1)[1])
rownames(im_timer1) <- gsub("-", ".", rownames(im_timer1))
total_scores <- rowSums(im_timer1)
proportion_scores <- sweep(im_timer1, 1, total_scores, "/")
timer1=proportion_scores
timer1 <- timer1 %>%
rownames_to_column(var = "Sample")
Cellratio <- timer1 %>%
pivot_longer(cols = -Sample, names_to = "CellName", values_to = "Proportion")
P1 <- ggplot(Cellratio) +
geom_bar(aes(x = Sample, y = Proportion, fill = CellName), stat = "identity", width =
width.figt2, linewidth = linewidth.fig2, colour = colour.fig2) +
theme_classic() +
labs(x = 'Sample', y = 'Ratio') +
coord_flip() +
theme(panel.border = element_rect(fill = NA, color = color.fig2, linewidth =
linewidth.figw2, linetype = "solid"))
# 保存为 PDF 文件
pdf(file="im_estimate/2.immune_stacke.pdf", width=width.fig2, height=height.fig2)
=# 打印图形
print(P1)
dev.off()
write.csv(proportion_scores, file = "im_estimate/proportion_scores.csv", row.names = T)
# 基因相关性分析
gene=gene.duo
gene_exp=expr[gene,]
gene_exp <- gene_exp[complete.cases(gene_exp), ]
gene_exp=t(gene_exp)
# 确保两个数据框拥有相同行名
common_rows <- intersect(rownames(gene_exp), rownames(im_timer1))
# 子集化数据框，仅保留相同行名的行
gene_exp_subset <- gene_exp[common_rows, ]
im_timer1_subset <- im_timer1[common_rows, ]
# 计算两个数据框所有列之间的相关性
correlation_matrix <- cor(gene_exp_subset, im_timer1_subset,method=method.fig3)
write.csv(correlation_matrix, file = "im_estimate/相关性数值.csv", row.names = T)
na_cols <- apply(is.na(correlation_matrix), 2, any)
correlation_matrix <- correlation_matrix[, !na_cols]
datacor=t(correlation_matrix)
# 绘制相关性热图
p2 <- pheatmap(datacor,scale = "none",cluster_cols = cluster,cluster_rows = cluster,
color = colorRampPalette(color.fig3)(50),
legend = T,show_colnames = T,show_rownames = T,
fontsize = fontsize,cellheight = cellheight,cellwidth =
cellwidth,heatmap_legend_param = list(title = ""))
# 保存为 PDF 文件
pdf(file="im_estimate/3. 多 个 基 因 免 疫 相 关 性 图 .pdf", width=width.fig3,
height=height.fig3)
=print(p2)
dev.off()
# 循环遍历每个基因，生成并保存气泡图
for (siggene in rownames(correlation_matrix)) {
# 提取特定基因的相关性数据
siggene_cor <- as.data.frame(t(correlation_matrix[siggene, , drop = FALSE]))
siggene_cor$immune <- rownames(siggene_cor)
# 提取相关性值
cor_values <- siggene_cor[, 1]
# 提取第一列的列名作为标题
title_name <- names(siggene_cor)[1]
# 绘制气泡图
p3 <- ggplot(data = siggene_cor, aes(x = cor_values, y = immune)) +
labs(x = "Correlation coefficient", y = "Immune cell") +
geom_point(aes(size = abs(cor_values), color = cor_values)) +
scale_color_gradient2(low = low.fig4, mid = mid.fig4, high = high.fig4, midpoint = 0)
+
theme(panel.background = element_rect(fill = "white", size = 1, color = "black"),
panel.grid = element_line(color = "grey75", size = 0.5),
axis.ticks = element_line(size = 0.5),
axis.text.y = element_text(colour = "black", size = 9),
plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
xlab("Correlation coefficient") +
ylab("Immune cell") +
ggtitle(title_name)
# 保存为 PDF 文件，文件名包含基因名
pdf(file = paste0("im_estimate/", siggene, "_immune_correlation.pdf"), width =
width.fig4, height = height.fig4)
print(p3)
dev.off()
}


##Figure 4D, E 

if (!dir.exists("im_quantiseq")) {
dir.create("im_quantiseq")
}
expr=exprim
cluster=cluster1
im_quantiseq <- deconvo_tme(eset = expr,
method = "quantiseq",
scale_mrna = scale_mrna)
write.csv(im_quantiseq, file = "im_quantiseq/im_quantiseq.csv", row.names = FALSE)
# 去掉多余的列
data <- im_quantiseq
data <- as.data.frame(data)
data$Group <- group_labels
data1 <- data[ , -1]
# 转换数据格式并添加分组信息
df <- reshape2::melt(data1, id.vars = "Group", variable.name = "Cell", value.name =
"score")
p <- ggplot(data=df, aes(x=Cell, y=score, fill=Group)) +
geom_boxplot(width=0.3, position = position_dodge(0.5), outlier.colour = NA) +
theme_bw() +
theme(panel.grid = element_blank()) +
scale_fill_manual(values = values) +
stat_compare_means(aes(group=Group), method = method.fig1, label="p.signif",
hide.ns = TRUE, label.y = c(label.y)) +
theme(axis.text.x = element_text(size=size.fig1, colour=colour.fig1, angle=angle.fig1x,
hjust = 1),
axis.text.y = element_text(size=size.fig1, angle=0),
axis.title.x = element_text(size=size.fig1),
axis.title.y = element_text(size=size.fig1),
legend.text = element_text(size=size.fig1)) + # 控制图例中文本字体大小
xlab("") # 去掉 x 轴标签
# 打印图形
print(p)
# 保存为 PDF 文件
ggsave("im_quantiseq/1.immune_boxplot.pdf", plot = p, width = width.fig1, height =
height.fig1)
# 堆叠图的绘制
im_timer1=im_quantiseq
im_timer1 <- column_to_rownames(im_timer1, var = colnames(im_timer1)[1])
rownames(im_timer1) <- gsub("-", ".", rownames(im_timer1))
total_scores <- rowSums(im_timer1)
proportion_scores <- sweep(im_timer1, 1, total_scores, "/")
timer1=proportion_scores
timer1 <- timer1 %>%
rownames_to_column(var = "Sample")
Cellratio <- timer1 %>%
pivot_longer(cols = -Sample, names_to = "CellName", values_to = "Proportion")
P1 <- ggplot(Cellratio) +
geom_bar(aes(x = Sample, y = Proportion, fill = CellName), stat = "identity", width =
width.figt2, linewidth = linewidth.fig2, colour = colour.fig2) +
theme_classic() +
labs(x = 'Sample', y = 'Ratio') +
coord_flip() +
theme(panel.border = element_rect(fill = NA, color = color.fig2, linewidth =
linewidth.figw2, linetype = "solid"))
# 保存为 PDF 文件
pdf(file="im_quantiseq/2.immune_stacke.pdf", width=width.fig2, height=height.fig2)
# 打印图形
print(P1)
dev.off()
write.csv(proportion_scores, file = "im_quantiseq/proportion_scores.csv", row.names = T)
# 基因相关性分析
gene=gene.duo
gene_exp=expr[gene,]
gene_exp <- gene_exp[complete.cases(gene_exp), ]
gene_exp=t(gene_exp)
# 确保两个数据框拥有相同行名
common_rows <- intersect(rownames(gene_exp), rownames(im_timer1))
# 子集化数据框，仅保留相同行名的行
gene_exp_subset <- gene_exp[common_rows, ]
im_timer1_subset <- im_timer1[common_rows, ]
# 计算两个数据框所有列之间的相关性
correlation_matrix <- cor(gene_exp_subset, im_timer1_subset,method=method.fig3)
write.csv(correlation_matrix, file = "im_quantiseq/相关性数值.csv", row.names = T)
na_cols <- apply(is.na(correlation_matrix), 2, any)
correlation_matrix <- correlation_matrix[, !na_cols]
datacor=t(correlation_matrix)
# 绘制相关性热图
p2 <- pheatmap(datacor,scale = "none",cluster_cols = cluster,cluster_rows = cluster,
color = colorRampPalette(color.fig3)(50),
legend = T,show_colnames = T,show_rownames = T,
fontsize = fontsize,cellheight = cellheight,cellwidth =
cellwidth,heatmap_legend_param = list(title = ""))
# 保存为 PDF 文件
pdf(file="im_quantiseq/3. 多 个 基 因 免 疫 相 关 性 图 .pdf", width=width.fig3,
height=height.fig3)
print(p2)
dev.off()
# 循环遍历每个基因，生成并保存气泡图
for (siggene in rownames(correlation_matrix)) {
# 提取特定基因的相关性数据
siggene_cor <- as.data.frame(t(correlation_matrix[siggene, , drop = FALSE]))
siggene_cor$immune <- rownames(siggene_cor)
# 提取相关性值
cor_values <- siggene_cor[, 1]
# 提取第一列的列名作为标题
title_name <- names(siggene_cor)[1]
# 绘制气泡图
p3 <- ggplot(data = siggene_cor, aes(x = cor_values, y = immune)) +
labs(x = "Correlation coefficient", y = "Immune cell") +
geom_point(aes(size = abs(cor_values), color = cor_values)) +
scale_color_gradient2(low = low.fig4, mid = mid.fig4, high = high.fig4, midpoint = 0)
+
theme(panel.background = element_rect(fill = "white", size = 1, color = "black"),
panel.grid = element_line(color = "grey75", size = 0.5),
axis.ticks = element_line(size = 0.5),
axis.text.y = element_text(colour = "black", size = 9),
plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
xlab("Correlation coefficient") +
ylab("Immune cell") +
ggtitle(title_name)
# 保存为 PDF 文件，文件名包含基因名
pdf(file = paste0("im_quantiseq/", siggene, "_immune_correlation.pdf"), width =
width.fig4, height = height.fig4)
print(p3)
dev.off()
}


#Figures 7A, 9A-C


# 加载必要的包
library(ComplexHeatmap)
library(circlize)

# 创建相关性矩阵（根据您提供的表格数据）
cor_matrix <- matrix(c(
  1.00, 0.11, 0.42, -0.02, 0.27, 0.59, 0.36, 0.11, 0.45, 0.21, 0.20, 0.25, 0.17, 0.34, 0.53, 0.13, 0.24, -0.10,
  0.11, 1.00, 0.33, -0.10, 0.38, 0.14, 0.42, 0.28, 0.31, 0.30, -0.30, -0.02, 0.13, 0.47, 0.17, 0.55, 0.25, 0.09,
  0.42, 0.33, 1.00, -0.13, 0.33, 0.08, 0.53, 0.13, 0.29, 0.42, -0.05, 0.10, 0.14, 0.25, 0.20, 0.40, 0.25, 0.03,
  -0.02, -0.10, -0.13, 1.00, 0.18, 0.04, -0.09, 0.12, 0.04, 0.14, -0.16, -0.04, 0.10, 0.04, 0.15, -0.05, 0.15, 0.04,
  0.27, 0.38, 0.33, 0.18, 1.00, 0.41, 0.41, 0.36, 0.51, 0.21, -0.16, 0.13, 0.04, 0.33, 0.17, 0.42, 0.56, -0.05,
  0.59, 0.14, 0.08, 0.04, 0.41, 1.00, 0.15, 0.19, 0.52, -0.08, 0.25, 0.31, -0.02, 0.39, 0.38, -0.05, 0.14, -0.24,
  0.36, 0.42, 0.53, -0.09, 0.41, 0.15, 1.00, 0.20, 0.33, 0.44, -0.08, 0.06, 0.27, 0.33, 0.16, 0.66, 0.34, 0.20,
  0.11, 0.28, 0.13, 0.12, 0.36, 0.19, 0.20, 1.00, 0.16, 0.16, -0.16, 0.07, 0.25, 0.31, 0.15, 0.24, 0.29, 0.02,
  0.45, 0.31, 0.29, 0.04, 0.51, 0.52, 0.33, 0.16, 1.00, 0.25, -0.11, 0.07, 0.25, 0.41, 0.33, 0.22, 0.29, -0.05,
  0.21, 0.30, 0.42, 0.14, 0.21, -0.08, 0.44, 0.16, 0.25, 1.00, -0.30, 0.05, 0.31, 0.22, 0.27, 0.47, 0.20, 0.21,
  0.20, -0.30, -0.05, -0.16, -0.16, 0.25, -0.08, -0.16, -0.11, -0.30, 1.00, 0.26, -0.06, 0.06, -0.05, -0.26, -0.13, -0.18,
  0.25, -0.02, 0.10, -0.04, 0.13, 0.31, 0.06, 0.07, 0.07, 0.05, 0.26, 1.00, 0.15, 0.06, 0.49, -0.12, 0.19, -0.15,
  0.17, 0.13, 0.14, 0.10, 0.04, -0.02, 0.27, 0.25, 0.25, 0.31, -0.06, 0.15, 1.00, 0.18, 0.34, 0.11, 0.17, 0.39,
  0.34, 0.47, 0.25, 0.04, 0.33, 0.39, 0.33, 0.31, 0.41, 0.22, 0.06, 0.06, 0.18, 1.00, 0.16, 0.40, 0.09, -0.15,
  0.53, 0.17, 0.20, 0.15, 0.17, 0.38, 0.16, 0.15, 0.33, 0.27, -0.05, 0.49, 0.34, 0.16, 1.00, -0.07, 0.32, 0.19,
  0.13, 0.55, 0.40, -0.05, 0.42, -0.05, 0.66, 0.24, 0.22, 0.47, -0.26, -0.12, 0.11, 0.40, -0.07, 1.00, 0.28, 0.09,
  0.24, 0.25, 0.25, 0.15, 0.56, 0.14, 0.34, 0.29, 0.29, 0.20, -0.13, 0.19, 0.17, 0.09, 0.32, 0.28, 1.00, 0.11,
  -0.10, 0.09, 0.03, 0.04, -0.05, -0.24, 0.20, 0.02, -0.05, 0.21, -0.18, -0.15, 0.39, -0.15, 0.19, 0.09, 0.11, 1.00
), nrow = 18, byrow = TRUE)

# 设置行名和列名（根据图片）
gene_names <- c("IL4I1", "CD274", "CD40", "CD44", "CD48", "CD86", 
                "CTLA4", "IDO1", "IL10", "LAG3", "LAIR1", "LGALS9", 
                "PDCD1", "PDCD1LG2", "TGFB1", "TIGIT", "TNFRSF14", "TNFRSF4")

rownames(cor_matrix) <- gene_names
colnames(cor_matrix) <- gene_names

# 设置颜色方案（蓝-白-红）
color_palette <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# 绘制热图
Heatmap(cor_matrix,
        name = "Correlation",
        col = color_palette,
        
        # 聚类设置
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        # 行列名称设置
        row_names_side = "left",
        column_names_side = "top",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45,
        
        # 单元格内容设置
        cell_fun = function(j, i, x, y, width, height, fill) {
          # 只在矩阵的上三角部分显示数值
          if (i <= j) {
            grid.text(sprintf("%.2f", cor_matrix[i, j]), 
                     x, y, gp = gpar(fontsize = 8))
          }
        },
        
        # 热图外观
        rect_gp = gpar(col = "white", lwd = 1),
        
        # 图例设置
        heatmap_legend_param = list(
          title = "Correlation",
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          at = c(-1, -0.5, 0, 0.5, 1),
          labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
          legend_height = unit(4, "cm")
        ),
        
        # 标题
        column_title = "Correlation between IL4I1 and Immune Checkpoint Genes",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        
        # 调整热图大小
        width = unit(15, "cm"),
        height = unit(15, "cm"))



# 加载必要的包
library(ComplexHeatmap)
library(circlize)

# 创建相关性矩阵（根据您提供的表格数据）
cor_matrix <- matrix(c(
  1.00, 0.53, 0.45, 0.47, 0.58, 0.50, -0.02, 0.07, 0.23, 0.52, 0.24,
  0.53, 1.00, 0.33, 0.34, 0.29, 0.35, -0.08, -0.00, 0.19, 0.38, 0.18,
  0.45, 0.33, 1.00, 0.67, 0.57, 0.45, -0.26, 0.08, 0.45, 0.25, 0.11,
  0.47, 0.34, 0.67, 1.00, 0.79, 0.43, -0.18, 0.03, 0.38, 0.28, 0.11,
  0.58, 0.29, 0.57, 0.79, 1.00, 0.60, -0.13, 0.09, 0.32, 0.28, 0.04,
  0.50, 0.35, 0.45, 0.43, 0.60, 1.00, -0.09, 0.08, 0.26, 0.23, 0.11,
  -0.02, -0.08, -0.26, -0.18, -0.13, -0.09, 1.00, 0.04, -0.17, 0.02, -0.11,
  0.07, -0.00, 0.08, 0.03, 0.09, 0.08, 0.04, 1.00, 0.08, 0.12, -0.13,
  0.23, 0.19, 0.45, 0.38, 0.32, 0.26, -0.17, 0.08, 1.00, 0.25, 0.27,
  0.52, 0.38, 0.25, 0.28, 0.28, 0.23, 0.02, 0.12, 0.25, 1.00, 0.23,
  0.24, 0.18, 0.11, 0.11, 0.04, 0.11, -0.11, -0.13, 0.27, 0.23, 1.00
), nrow = 11, byrow = TRUE)

# 设置行名和列名（根据图片）
gene_names <- c("IL4I1", "TGFB1", "IL10", "CD163", "MS4A4A", "CSF1R", 
                "IL12A", "NOS2", "PTGS2", "CCL3", "TNF")

rownames(cor_matrix) <- gene_names
colnames(cor_matrix) <- gene_names

# 设置颜色方案（蓝-白-红）
color_palette <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# 1. 绘制完整相关性矩阵
Heatmap(cor_matrix,
        name = "Correlation",
        col = color_palette,
        
        # 聚类设置
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        # 行列名称设置
        row_names_side = "left",
        column_names_side = "top",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45,
        
        # 单元格内容设置
        cell_fun = function(j, i, x, y, width, height, fill) {
          # 显示数值（两位小数）
          grid.text(sprintf("%.2f", cor_matrix[i, j]), 
                   x, y, gp = gpar(fontsize = 8))
        },
        
        # 热图外观
        rect_gp = gpar(col = "white", lwd = 1),
        
        # 图例设置
        heatmap_legend_param = list(
          title = "Correlation",
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          at = c(-1, -0.5, 0, 0.5, 1),
          labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
          legend_height = unit(4, "cm")
        ),
        
        # 标题
        column_title = "Correlation between IL4I1 and Macrophage Markers",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        
        # 调整热图大小
        width = unit(12, "cm"),
        height = unit(12, "cm"))


#Figures 9B, C

# 加载必要的包
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(ggpubr)
library(dplyr)

# 1. 从TCGA数据库下载AML表达数据
query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

# 2. 提取表达矩阵（使用FPKM或TPM）
# 根据文章方法部分，使用log2(FPKM+1)转换
expr_fpkm <- assay(data, "fpkm_unstrand")  # 或使用 "tpm_unstrand"

# 3. 提取目标基因：IL4I1, IL10, TGFB1
target_genes <- c("IL4I1", "IL10", "TGFB1")

# 查找基因ID对应关系
gene_info <- rowData(data)
# 方法1：通过基因名查找
gene_ids <- gene_info[gene_info$gene_name %in% target_genes, ]

# 提取表达值
expr_subset <- expr_fpkm[rownames(gene_ids), ]
rownames(expr_subset) <- gene_ids$gene_name

# 4. 数据转换：log2(FPKM+1)
expr_log2 <- log2(expr_subset + 1)

# 5. 转换为数据框便于绘图
expr_df <- as.data.frame(t(expr_log2))

# 6. 绘制IL4I1与IL10的相关性散点图
p1 <- ggplot(expr_df, aes(x = IL4I1, y = IL10)) +
  geom_point(alpha = 0.7, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid", 
              fill = "lightpink", alpha = 0.3) +
  
  # 添加统计信息（使用Spearman相关性）
  stat_cor(
    method = "spearman",
    label.x = min(expr_df$IL4I1) + 0.5,
    label.y = max(expr_df$IL10) - 0.5,
    size = 5,
    color = "black"
  ) +
  
  # 坐标轴标签
  labs(
    x = expression("The expression of IL4I1 Log"[2]*"(FPKM+1)"),
    y = expression("The expression of IL10 Log"[2]*"(FPKM+1)"),
    title = "Correlation between IL4I1 and IL10 expression"
  ) +
  
  # 主题设置
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  ) +
  
  # 添加网格线
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank())

# 7. 绘制IL4I1与TGFB1的相关性散点图
p2 <- ggplot(expr_df, aes(x = IL4I1, y = TGFB1)) +
  geom_point(alpha = 0.7, color = "darkgreen", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid", 
              fill = "lightpink", alpha = 0.3) +
  
  # 添加统计信息
  stat_cor(
    method = "spearman",
    label.x = min(expr_df$IL4I1) + 0.5,
    label.y = max(expr_df$TGFB1) - 0.5,
    size = 5,
    color = "black"
  ) +
  
  # 坐标轴标签
  labs(
    x = expression("The expression of IL4I1 Log"[2]*"(FPKM+1)"),
    y = expression("The expression of TGFB1 Log"[2]*"(FPKM+1)"),
    title = "Correlation between IL4I1 and TGFB1 expression"
  ) +
  
  # 主题设置
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "none"
  ) +
  
  # 添加网格线
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank())

# 8. 并排显示两个图
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

# 9. 也可以单独保存每个图
ggsave("IL4I1_IL10_correlation.png", p1, width = 6, height = 5, dpi = 300)
ggsave("IL4I1_TGFB1_correlation.png", p2, width = 6, height = 5, dpi = 300)

