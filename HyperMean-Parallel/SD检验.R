SD_test <- function(sam1, sam2) {
  # 记录函数开始执行时间
  start_time <- Sys.time()
  
  # 计算样本大小和维度
  n1 <- nrow(sam1)  # 第一组样本的观测数
  n2 <- nrow(sam2)  # 第二组样本的观测数
  num_features <- ncol(sam1)  # 特征维度
  total_n <- n1 + n2 - 2  # 合并自由度
  
  # 计算标准化因子
  scaling_factor <- (n1 * n2) / (n1 + n2)
  
  # 计算组间均值差异
  mean_sam1 <- colMeans(sam1)  # 第一组样本的均值向量
  mean_sam2 <- colMeans(sam2)  # 第二组样本的均值向量
  mean_diff <- mean_sam1 - mean_sam2  # 组间均值差异向量
  
  # 计算合并协方差矩阵
  pooled_cov <- ((n1 - 1) * cov(sam1) + (n2 - 1) * cov(sam2)) / total_n
  diag_vars <- diag(pooled_cov)  # 提取方差项（对角线元素）
  
  # 计算相关矩阵和其平方的迹
  cor_matrix <- cov2cor(pooled_cov)  # 将协方差矩阵转换为相关矩阵
  trace_cor_sq <- sum(cor_matrix^2)  # 相关矩阵平方的迹
  
  # 计算检验统计量
  normalized_diff <- sum(mean_diff^2 / diag_vars)  # 方差归一化的均值差异
  numerator <- scaling_factor * normalized_diff - (total_n * num_features) / (total_n - 2)
  
  adjustment_term <- 1 + trace_cor_sq / (num_features^1.5)  # Srivastava调整系数
  denominator <- sqrt(2 * (trace_cor_sq - (num_features^2) / total_n) * adjustment_term)
  
  test_statistic <- numerator / denominator  
  
  # 计算单侧p值
  p_value <- pnorm(-test_statistic)
  
  # 计算总耗时
  computation_time <- Sys.time() - start_time
  
  # 返回结果
  return(list(
    p_value = p_value,
    computation_time = computation_time
  ))
}
SD_test(sam1,sam2)
