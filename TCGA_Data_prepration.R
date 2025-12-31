
# 设置工作目录，设置当前工作目录为TCGA下载与整理的文件夹
setwd("E:/tcga") #

# 加载必要的包
library(jsonlite)
library(tidyverse)
library(dplyr)

# 读取JSON文件中的元数据
meta_data <- fromJSON("metadata.cart.2025-11-24.json") # 读取JSON格式的元数据文件

# 构建样本信息数据框
# 使用map2_df函数结合元数据中的关联实体信息与文件名，创建一个样本信息的数据框
samples_df <- map2_df(meta_data$associated_entities, seq_along(meta_data$associated_entities),
                      ~ tibble(sample_id = .x[,1], file_name = meta_data$file_name[.y]))

# 获取计数数据文件列表
count_file_paths <- list.files('gdc_download_20251124_091719.405130/', pattern = '*.tsv', recursive = TRUE) # 获取当前目录下所有tsv文件的路径

# 提取纯文件名
file_names_only <- sapply(strsplit(count_file_paths, split='/'), function(x) x[2]) # 提取文件名，不包含路径信息

# 创建空的表达矩阵
expr_matrix <- data.frame() # 初始化空的数据框，用于存储表达数据

# 循环读取每个文件，并生成TPM类型文件
for (file_index in seq_along(count_file_paths)) {
  full_path <- paste0('gdc_download_20251124_091719.405130/', count_file_paths[file_index]) # 拼接完整的文件路径
  count_data <- read.delim(full_path, fill = TRUE, header = FALSE, row.names = 1) # 读取计数数据，使用行名作为列名
  colnames(count_data) <- count_data[2,] # 设置列名，跳过第一行
  count_data <- count_data[-(1:6),] # 删除前六行，这些行包含 header 信息
  
  # 提取感兴趣的TPM数据列
  tpm_data <- count_data[6]
  
  # 设置TPM数据的列名
  sample_id <- samples_df$sample_id[which(samples_df$file_name == file_names_only[file_index])][[1]]
  colnames(tpm_data) <- sample_id
  
  if (nrow(expr_matrix) == 0) {
    # 如果expr_matrix还是空的，直接赋值
    expr_matrix <- tpm_data
  } else {
    # 检查行名是否一致，确保数据的一致性
    if (!all(rownames(expr_matrix) == rownames(tpm_data))) {
      stop("行名不一致，无法合并数据")
    }
    expr_matrix <- cbind(expr_matrix, tpm_data)
  }
}


# 循环读取每个文件，并生成count类型文件
for (file_index in seq_along(count_file_paths)) {
  full_path <- paste0('gdc_download_20251124_091719.405130/', count_file_paths[file_index]) # 拼接完整的文件路径
  count_data <- read.delim(full_path, fill = TRUE, header = FALSE, row.names = 1) # 读取计数数据，使用行名作为列名
  colnames(count_data) <- count_data[2,] # 设置列名，跳过第一行
  count_data <- count_data[-(1:6),] # 删除前六行，这些行包含 header 信息
  
  # 提取感兴趣的TPM数据列
  tpm_data <- count_data[3]
  
  # 设置TPM数据的列名
  sample_id <- samples_df$sample_id[which(samples_df$file_name == file_names_only[file_index])][[1]]
  colnames(tpm_data) <- sample_id
  
  if (nrow(expr_matrix) == 0) {
    # 如果expr_matrix还是空的，直接赋值
    expr_matrix <- tpm_data
  } else {
    # 检查行名是否一致，确保数据的一致性
    if (!all(rownames(expr_matrix) == rownames(tpm_data))) {
      stop("行名不一致，无法合并数据")
    }
    expr_matrix <- cbind(expr_matrix, tpm_data)
  }
}
  # 循环读取每个文件，并生成fpkm类型文件
for (file_index in seq_along(count_file_paths)) {
  full_path <- paste0('gdc_download_20251124_091719.405130/', count_file_paths[file_index]) # 拼接完整的文件路径
  count_data <- read.delim(full_path, fill = TRUE, header = FALSE, row.names = 1) # 读取计数数据，使用行名作为列名
  colnames(count_data) <- count_data[2,] # 设置列名，跳过第一行
  count_data <- count_data[-(1:6),] # 删除前六行，这些行包含 header 信息
  
  # 提取感兴趣的TPM数据列
  tpm_data <- count_data[7]
  
  # 设置TPM数据的列名
  sample_id <- samples_df$sample_id[which(samples_df$file_name == file_names_only[file_index])][[1]]
  colnames(tpm_data) <- sample_id
  
  if (nrow(expr_matrix) == 0) {
    # 如果expr_matrix还是空的，直接赋值
    expr_matrix <- tpm_data
  } else {
    # 检查行名是否一致，确保数据的一致性
    if (!all(rownames(expr_matrix) == rownames(tpm_data))) {
      stop("行名不一致，无法合并数据")
    }
    expr_matrix <- cbind(expr_matrix, tpm_data)
  }
} 
# 读取基因信息
gene_info <- as.matrix(read.delim(paste0('gdc_download_20251124_091719.405130/', count_file_paths[1]), fill = TRUE, header = FALSE, row.names = 1))
genes <- gene_info[-(1:6), 1] # 获取基因名称
gene_types <- gene_info[-(1:6), 2] # 获取基因类型
    
# 添加基因符号和类型
expr_matrix <- cbind(gene_type = gene_types, gene_symbol = genes, expr_matrix) # 将基因类型和符号添加到表达矩阵中
    
# 聚合数据，保留最大表达量
Bexpr_matrix <- aggregate(. ~ gene_symbol, data = expr_matrix, max) # 对每个基因符号的最大表达量进行聚合
# 整合
expr_matrix <- semi_join(expr_matrix, Bexpr_matrix, by = "gene_symbol") # 合并原始表达矩阵和聚合后的表达矩阵
expr_matrix <- Bexpr_matrix # 将聚合后的矩阵赋值给expr_matrix

    
# 将gene_symbol列设为行名,并转化为导出格式
rownames(expr_matrix) <- expr_matrix[, "gene_symbol"] # 将基因符号列设为行名
expr_matrix <- expr_matrix[, -1] # 移除基因符号列
expr_matrix <- data.frame(ID = rownames(expr_matrix), expr_matrix) # 创建ID列，并将表达数据转化为data.frame格式
colnames(expr_matrix) <- gsub('[.]', '-', colnames(expr_matrix)) # 将列名中的点替换为连字符，以符合导出文件的要求
expr_matrix <- subset(x = expr_matrix, gene_type == "protein_coding") # 筛选出蛋白质编码基因的表达数据
expr_matrix <- expr_matrix[, -1] # 移除多余的列
expr_matrix <- expr_matrix[, -1] # 再次移除多余的列（可能是由于上一步筛选导致的额外列）
# 导出最终表格
write.table(expr_matrix, 'TCGAAAAAA_LUSC_TPM.txt', sep="\t", quote=FALSE, col.names = NA) # 将最终的表达矩阵导出为tsv格式文件

（2）GTEx数据下载后处理
# 1. 加载数据（假设是GTEx表达矩阵）
# gtex_data <- read.csv("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct", sep = "\t", skip = 2)

# 2. 简化示例数据
gtex_data <- data.frame(
    gene_id = c("ENSG00000223972.5", "ENSG00000227232.5", "ENSG00000243485.5"),
    Heart = c(10.5, 8.2, 15.3),
    Liver = c(9.1, 7.8, 14.2),
    Brain = c(11.2, 9.1, 16.5)
)

# 3. 提取基本基因ID（去掉版本号）
gtex_data$gene_id_base <- sub("\\.[0-9]+$", "", gtex_data$gene_id)

# 4. 使用org.Hs.eg.db进行注释
library(org.Hs.eg.db)
library(AnnotationDbi)

# 创建注释映射
symbols <- mapIds(org.Hs.eg.db,
                  keys = gtex_data$gene_id_base,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

gene_names <- mapIds(org.Hs.eg.db,
                     keys = gtex_data$gene_id_base,
                     column = "GENENAME",
                     keytype = "ENSEMBL",
                     multiVals = "first")

entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gtex_data$gene_id_base,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# 5. 添加到数据框
gtex_data$symbol <- symbols
gtex_data$gene_name <- gene_names
gtex_data$entrez_id <- entrez_ids

# 6. 查看结果
head(gtex_data)

# 7. 可选：移除没有注释的基因
gtex_clean <- gtex_data[!is.na(gtex_data$symbol), ]

# 8. 可选：保存结果
# write.csv(gtex_clean, "gtex_annotated.csv", row.names = FALSE)


（3）按照基因名对数据进行分组
geneGroup_data_file = 'maxmean__unique.csv' #表达矩阵文件
group_by = 'IL4I1' #按哪个基因分组
if (!dir.exists("中位数分组结果")) {
dir.create("中位数分组结果")
}
data_Select_group <- read.csv(geneGroup_data_file, header = TRUE,
row.names = 1, check.names = FALSE)
data = data_Select_group
prompt <- "按中位数进行分组，开始"
print(prompt)
median_value <- median(as.numeric(data[group_by, ]), na.rm = TRUE)
cat(" 该基因中位数为：")
print(median_value)
new_colnames <- ifelse(data[group_by, ] > median_value, paste(colnames(data),
"_high", sep = ""), paste(colnames(data), "_low", sep = ""))
colnames(data) <- new_colnames
print("分组后的列名")
print(colnames(data))
data_gp_byMedian <- data[, order(-grepl("_high$", names(data)),
grepl("_low$", names(data)))]
print("对列名进行排序")
print(colnames(data_gp_byMedian))
write.csv(data_gp_byMedian, file = "中位数分组结果/new_group_by_median.csv",
row.names = T)
prompt <- "结束，分组后的文件已下载，文件名为 new_group_by_median"
cat("此处无需写作交代\n")

（4）TCGA 分组排序
TCGA_group_file = 'lusc_mRNA_count.csv' #表达矩阵文件的名称
if (!dir.exists("TCGA 分组排序结果")) {
dir.create("TCGA 分组排序结果")
}
data_group <- read.csv(TCGA_group_file, header = TRUE, row.names = 1,
check.names = FALSE)
col_names <- colnames(data_group)
sort_function <- function(col_name) {
num <- as.numeric(substr(col_name, 14, 15))
return(num < 10)
}
sorted_col_names <- c(col_names[sapply(col_names, sort_function)],
col_names[!sapply(col_names, sort_function)])
data_group1 <- data_group[, sorted_col_names]
write.csv(data_group1, file = "TCGA 分组排序结果/after_group_TCGA.csv",
row.names = T)
cat("此处无需写作交代\n")

（12）去除低表达基因
data_DEGA <- read.csv(file, header = TRUE, row.names = 1)
#去除低质量、低表达基因
data10=data_DEGA
data=data_DEGA
data10 = as.data.frame(data10)
data10=data[apply(cpm(data),1,sum)>zuidiz,]
data10 <- round(as.matrix(data10))
write.csv(data10, "去除低表达基因（Count）.csv", row.names = TRUE)
gene_means <- rowMeans(data_DEGA)
filtered_data <- data_DEGA[gene_means > threshold, ]
write.csv(filtered_data, "去除低表达基因（非Count）.csv", row.names =
TRUE)
#去除0值的比例
data0123=data_DEGA
expMatrix <- data0123
#去除含有大量0值的样本
zero_counts <- rowSums(expMatrix == 0)
threshold <- ncol(expMatrix) * bili0
expMatrix <- expMatrix[zero_counts <= threshold, ]
write.csv(expMatrix, "去除高0值比例.csv", row.names = TRUE)



