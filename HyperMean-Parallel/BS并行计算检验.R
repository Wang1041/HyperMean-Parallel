BS_parallel_test <- function(sam1, sam2, num_cores=availableCores() - 1) {
  # 设置并行运算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # 记录函数并行开始执行的时间
  start_time <- Sys.time()
  
  # ======================== 数据准备与分块 ========================
  # 获取样本维度信息
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  num_features <- ncol(sam1)
  degrees_of_freedom <- n1 + n2 - 2
  scaling_factor <- (n1 * n2) / (n1 + n2)
  
  # 确定数据分块数量
  num_groups <- num_cores
  
  # 计算每块的特征数和余数
  features_per_group <- floor(num_features / num_groups)
  remainder_features <- num_features %% num_groups
  
  # 数据分块函数：将高维特征分成多个子集
  split_features <- function(data) {
    lapply(1:num_groups, function(group_idx) {
      # 计算当前组的特征索引范围
      if (group_idx <= remainder_features) {
        start_idx <- (group_idx - 1) * (features_per_group + 1) + 1
        end_idx <- group_idx * (features_per_group + 1)
      } else {
        start_idx <- (group_idx - 1) * features_per_group + remainder_features + 1
        end_idx <- min(group_idx * features_per_group + remainder_features, num_features)
      }
      data[, start_idx:end_idx, drop = FALSE]
    })
  }
  
  # 将两个样本集按特征分块
  sam1_groups <- split_features(sam1)
  sam2_groups <- split_features(sam2)
  
  # ======================== 并行计算核心统计量 ========================
  
  # 并行计算各组统计量
  group_results <- foreach(
    group_idx = 1:num_groups,
    .combine = rbind
  ) %dopar% {
    # 获取当前特征组的子样本
    sam1_sub <- sam1_groups[[group_idx]]
    sam2_sub <- sam2_groups[[group_idx]]
    
    # 计算组间均值差异
    mean_sam1_sub <- colMeans(sam1_sub)
    mean_sam2_sub <- colMeans(sam2_sub)
    mean_diff <- mean_sam1_sub - mean_sam2_sub
    squared_diff_norm <- sum(mean_diff^2)
    
    # 计算合并方差估计
    # 样本1的离差平方和
    dev_sq_sam1 <- sum((sam1_sub - 
                          matrix(mean_sam1_sub, nrow(sam1_sub), ncol(sam1_sub), byrow = TRUE))^2)
    # 样本2的离差平方和
    dev_sq_sam2 <- sum((sam2_sub - 
                          matrix(mean_sam2_sub, nrow(sam2_sub), ncol(sam2_sub), byrow = TRUE))^2)
    # 合并方差估计
    pooled_variance <- (dev_sq_sam1 + dev_sq_sam2) / (n1 + n2 - 2)
    
    # 计算协方差交叉项
    cov_cross_term <- 0
    for (j in 1:num_groups) {
      # 计算当前组与其他组的协方差
      cov_matrix <- ((n1 - 1) * cov(sam1_sub, sam1_groups[[j]]) + 
                       (n2 - 1) * cov(sam2_sub, sam2_groups[[j]])) / degrees_of_freedom
      cov_cross_term <- cov_cross_term + sum(cov_matrix^2)
    }
    
    # 返回当前组的计算结果
    c(squared_diff_norm, pooled_variance, cov_cross_term)
  }
  

  # ======================== 计算最终统计量和p值 ========================
  # 汇总所有分块结果
  total_results <- colSums(group_results)
  
  # 计算Bai检验统计量
  numerator <- scaling_factor * total_results[1] - total_results[2]
  trace_adjustment <- total_results[3] - (total_results[2]^2 / num_features)
  denominator <- sqrt(2 * degrees_of_freedom * (degrees_of_freedom + 1) / 
                        ((degrees_of_freedom - 1) * (degrees_of_freedom + 2)) * 
                        trace_adjustment)
  bai_statistic <- numerator / denominator
  
  # 计算单侧p值
  p_value <- pnorm(-bai_statistic)
  
  # 计算过程总耗时
  computation_time <- Sys.time() - start_time
  
  # 返回结果
  return(list(
    p_value = p_value,
    time = computation_time
  ))
  
  # 关闭并行环境
  stopCluster(cl)
}


