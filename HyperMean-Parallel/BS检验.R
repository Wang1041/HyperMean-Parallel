BS_test<- function(sam1, sam2) {
  # Bai & Saranadasa (1996) 高维两样本均值检验
  # 输入：
  #   sam1, sam2 - 样本矩阵（行=观测值，列=变量）

  # 记录函数开始执行时间（用于性能监控）
  start_time <- Sys.time()
  
  # 计算样本大小
  n1 <- nrow(sam1)  # 第一组样本的观测数
  n2 <- nrow(sam2)  # 第二组样本的观测数
  total_n <- n1 + n2 - 2  # 合并自由度(用于协方差估计)
  
  # 计算样本均值差异
  mean_sam1 <- colMeans(sam1)  # 第一组样本的均值向量
  mean_sam2 <- colMeans(sam2)  # 第二组样本的均值向量
  mean_difference <- mean_sam1 - mean_sam2  # 组间均值差异向量
  
  # 检验统计量计算
  scaling_factor <- (n1 * n2) / (n1 + n2)  # 标准化的缩放因子
  squared_norm_diff <- sum(mean_difference^2)  # 均值差异向量的L2范数平方
  
  # 合并协方差矩阵估计 (类似ANOVA的合并方差)
  cov <- ((n1 - 1) * cov(sam1) + (n2 - 1) * cov(sam2)) / total_n
  
  # 计算协方差矩阵的迹及其平方
  trace_cov <- sum(diag(cov))  # 协方差矩阵的迹
  trace_sq_cov <- sum(cov^2) - (trace_cov^2 / ncol(sam1))  # 迹的平方修正项
  
  # 构造Bai检验统计量
  numerator <- scaling_factor * squared_norm_diff - trace_cov  # 检验统计量分子
  denominator <- sqrt(2 * (total_n + 1) * total_n / ((total_n + 2) * (total_n - 1)) * trace_sq_cov)  # 分母
  bai_statistic <- numerator / denominator  # 最终检验统计量
  
  # 计算单侧p值
  p_value <- 1 - pnorm(bai_statistic)
  
  # 返回结果列表
  return(list(
    p_value = p_value,  # 检验的p值
    computation_time = Sys.time() - start_time  # 计算耗时
  ))
}



