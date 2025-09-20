CQ_test<- function(sam1, sam2) {
  # Chen et al.(2010)的高维两样本检验
  start_time <- Sys.time()
  
  # 基本样本信息
  n1 <- dim(sam1)[1]
  n2 <- dim(sam2)[1]
  p <- dim(sam1)[2]
  
  # 均值向量
  sam1_mean <- colSums(sam1)
  sam2_mean <- colSums(sam2)
  
  # U统计量估计协方差矩阵的迹
  sigma1_squared_estimate <- 0
  sigma1_sigma2_estimate <- 0
  
  for (j in 1:n1) {
    sam1_j <- sam1[j, ]
    
    # tr(Σ1²)估计：使用双重中心化
    for (k in 1:n1) {
      if (j != k) {
        sam1_k <- sam1[k, ]
        temp_mean1 <- (sam1_mean - sam1_j - sam1_k) / (n1 - 2)
        sigma1_squared_estimate <- sigma1_squared_estimate + 
          sum(sam1_j * (sam1_k - temp_mean1)) * 
          sum(sam1_k * (sam1_j - temp_mean1))
      }
    }
    
    # tr(Σ1Σ2)估计：交叉样本计算
    for (k in 1:n2) {
      sam2_k <- sam2[k, ]
      sigma1_sigma2_estimate <- sigma1_sigma2_estimate +
        sum(sam1_j * (sam2_k - (sam2_mean - sam2_k)/(n2 - 1))) *
        sum(sam2_k * (sam1_j - (sam1_mean - sam1_j)/(n1 - 1)))
    }
  }
  
  # 标准化迹估计
  trace_sigma1_squared <- sigma1_squared_estimate / (n1 * (n1 - 1))
  trace_sigma1_sigma2 <- sigma1_sigma2_estimate / (n1 * n2)
  
  # tr(Σ2²)估计
  sigma2_squared_estimate <- 0
  for (j in 1:n2) {
    for (k in 1:n2) {
      if (j != k) {
        sam2_j <- sam2[j, ]
        sam2_k <- sam2[k, ]
        temp_mean2 <- (sam2_mean - sam2_j - sam2_k) / (n2 - 2)
        sigma2_squared_estimate <- sigma2_squared_estimate +
          sum(sam2_j * (sam2_k - temp_mean2)) * 
          sum(sam2_k * (sam2_j - temp_mean2))
      }
    }
  }
  trace_sigma2_squared <- sigma2_squared_estimate / (n2 * (n2 - 1))
  
  # 检验统计量的方差估计
  deno <- 2/(n1*(n1-1))*trace_sigma1_squared + 
    2/(n2*(n2-1))*trace_sigma2_squared + 
    4/(n1*n2)*trace_sigma1_sigma2
  
  # 快速计算平方和（避免存储大矩阵）
  sam1_sum <- sum(sam1 %*% t(sam1)) - sum(diag(sam1 %*% t(sam1)))
  sam2_sum <- sum(sam2 %*% t(sam2)) - sum(diag(sam2 %*% t(sam2)))
  cross_sum <- sum(sam1 %*% t(sam2))
  
  # 检验统计量
  nume <- sam1_sum/(n1*(n1-1)) + sam2_sum/(n2*(n2-1)) - 2*cross_sum/(n1*n2)
  
  # 正态近似计算p值
  list(p_value = pnorm(-nume/sqrt(deno)),
       computation_time = Sys.time() - start_time)
}
CQ_test(sam1,sam2)
