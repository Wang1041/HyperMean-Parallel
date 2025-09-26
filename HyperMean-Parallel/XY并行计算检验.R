XY_parallel_test  <- function(sam1, sam2, N = 1000,num_cores=availableCores() - 1) {
  # Xue et al.(2020) 并行高维两样本均值检验
  # 使用并行Bootstrap方法计算p值
  # 输入：
  #   sam1, sam2 - 样本矩阵（行=观测值，列=变量）
  #   N - Bootstrap重抽样次数，默认1000次
  # 设置并行运算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  start <- Sys.time()
  
  # 获取样本维度
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  p <- ncol(sam1)
  
  # ========== 中心化处理 ==========
  # 转置方式实现列中心化（更高效）
  t1 <- t(t(sam1) - colMeans(sam1))  # 样本1中心化
  t2 <- t(t(sam2) - colMeans(sam2))  # 样本2中心化
  
  # ========== 并行Bootstrap设置 ==========
  block_size <- 10  # 每个并行任务处理10次Bootstrap重抽样
  num_blocks <- ceiling(N / block_size)  # 总任务块数
  
  # ========== 并行Bootstrap计算 ==========
  result_s_level <- foreach(b = 1:num_blocks, .combine = rbind) %dopar% {
    result_block <- NULL
    start_idx <- (b - 1) * block_size + 1
    end_idx <- min(b * block_size, N)
    
    # 处理当前任务块内的所有Bootstrap重抽样
    for (i in start_idx:end_idx) {
      # 生成标准正态随机权重
      E1 <- rnorm(n1)  # 样本1的随机权重
      E2 <- rnorm(n2)  # 样本2的随机权重
      
      # 计算Bootstrap统计量：最大标准化差异
      T_n.e <- max(abs(colSums(E1 * t1) / sqrt(n1) - colSums(E2 * t2) / sqrt(n2)))
      result_block <- rbind(result_block, T_n.e)
    }
    result_block  # 返回当前任务块的结果
  }
  
  # ========== 计算实际统计量 ==========
  T_XY <- max(abs(colMeans(sam1) - colMeans(sam2))) * sqrt(n1)
  
  # ========== 计算Bootstrap p值 ==========
  P_value <- sum(T_XY > result_s_level) / N
  
  # ========== 返回结果 ==========
  out <- list()
  out$用时 <- Sys.time() - start  # 计算用时
  out$P值 <- P_value            # Bootstrap p值
  # 关闭并行环境
  stopCluster(cl)
  return(out)
}



