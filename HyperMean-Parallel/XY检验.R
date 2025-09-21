XY_test <- function(sam1, sam2, N = 1000) {
  # Xue & Yao (2020) 高维两样本均值检验
  # 使用自助法（Bootstrap）计算p值
  # 输入：
  #   sam1, sam2 - 样本矩阵（行=观测值，列=变量）
  #   N - Bootstrap重抽样次数，默认1000次
  
  start <- Sys.time()
  
  # 获取样本维度
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  p <- ncol(sam1)  # 变量维度
  
  # ========== 中心化处理 ==========
  # 将每个样本减去其列均值（中心化）
  t1 <- sam1 - matrix(colMeans(sam1), byrow = TRUE, nrow = n1, ncol = p)
  t2 <- sam2 - matrix(colMeans(sam2), byrow = TRUE, nrow = n2, ncol = p)
  
  # ========== Bootstrap重抽样 ==========
  T_n.e <- numeric(N)  # 存储Bootstrap统计量
  
  for (i in 1:N) {
    # 生成标准正态随机变量
    E1 <- rnorm(n1)  # 样本1的随机权重
    E2 <- rnorm(n2)  # 样本2的随机权重
    
    # 计算Bootstrap统计量：最大标准化差异
    T_n.e[i] <- max(abs(colSums(E1 * t1) / sqrt(n1) - colSums(E2 * t2) / sqrt(n2)))
  }
  
  # ========== 计算实际统计量 ==========
  # 实际样本的标准化最大均值差异
  T_XY <- max(abs(colMeans(sam1) - colMeans(sam2))) * sqrt(n1)
  
  # ========== 计算p值 ==========
  # Bootstrap p值：实际统计量超过Bootstrap统计量的比例
  P_value <- sum(T_XY > T_n.e) / N
  
  # ========== 返回结果 ==========
  out <- list()
  out$p值 <- P_value  # Bootstrap p值
  out$time <- Sys.time() - start  # 计算时间
  return(out)
}
XY_test(sam1,sam2)

