SD_parallel_test <- function(sam1, sam2, num_cores = availableCores() - 1) {
  # 设置并行计算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  start_time <- Sys.time()
  
  # ======================== 数据准备与分块 ========================
  # 获取样本维度信息
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  num_features <- ncol(sam1)
  degrees_of_freedom <- n1 + n2 - 2
  scaling_factor <- (n1 * n2) / (n1 + n2)
  
  # 确定数据分块数量（等于CPU核心数）
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
  # 并行计算每个任务的结果
  group_results <- foreach(
    task_id = 1:num_groups ,
    .combine = rbind,
    .packages = "foreach"
  ) %dopar% {
    # 确定当前任务对应的特征组
    group_idx <- ((task_id - 1) %% num_groups) + 1
    
    # 获取当前特征组的子样本
    sam1_sub <- sam1_groups[[group_idx]]
    sam2_sub <- sam2_groups[[group_idx]]
    num_sub_features <- ncol(sam1_sub)
    
    # 计算组间均值差异
    mean_sam1_sub <- colMeans(sam1_sub)
    mean_sam2_sub <- colMeans(sam2_sub)
    mean_diff <- mean_sam1_sub - mean_sam2_sub
    
    # 计算合并方差估计（更精确的方法）
    # 样本1的离差平方和
    dev_sq_sam1 <- colSums((sweep(sam1_sub, 2, mean_sam1_sub, "-"))^2)
    # 样本2的离差平方和
    dev_sq_sam2 <- colSums((sweep(sam2_sub, 2, mean_sam2_sub, "-"))^2)
    # 合并方差估计
    pooled_variance <- (dev_sq_sam1 + dev_sq_sam2) / degrees_of_freedom
    
    # 计算方差归一化统计量
    normalized_diff <- sum(mean_diff^2 / pooled_variance)
    numerator_part <- scaling_factor * normalized_diff - 
                     (degrees_of_freedom * num_sub_features) / (degrees_of_freedom - 2)
    
    # 计算标准差倒数（用于相关矩阵计算）
    std_inv <- 1 / sqrt(pooled_variance)
    
    # 计算相关矩阵平方迹的部分
    cor_sq_sum <- 0
    for (j in 1:num_groups) {
      # 计算当前组与其他组的协方差
      cov_matrix <- ((n1 - 1) * cov(sam1_sub, sam1_groups[[j]]) + 
                   (n2 - 1) * cov(sam2_sub, sam2_groups[[j]])) / degrees_of_freedom
      
      # 应用标准化变换
      if (j == group_idx) {
        # 同一组时使用当前组的std_inv
        scaled_cov <- cov_matrix * outer(std_inv, std_inv)
      } else {
        # 获取j组的标准化因子
        std_inv_j <- 1 / sqrt(
          (colSums((sweep(sam1_groups[[j]], 2, colMeans(sam1_groups[[j]]), "-"))^2) +
          colSums((sweep(sam2_groups[[j]], 2, colMeans(sam2_groups[[j]]), "-"))^2)) / degrees_of_freedom
        )
        scaled_cov <- cov_matrix * outer(std_inv, std_inv_j)
      }
      
      # 累加平方和
      cor_sq_sum <- cor_sq_sum + sum(scaled_cov^2)
    }
    
    # 返回当前组的计算结果
    c(numerator_part, cor_sq_sum)
  }
  
  # ======================== 汇总结果并计算p值 ========================
    rep_results <- colSums(group_results)
    # 计算检验统计量
    numerator <- rep_results[1]
    trace_cor_sq <- rep_results[2]
    
    # 计算Srivastava调整因子
    adjustment_factor <- 1 + trace_cor_sq / (num_features^1.5)
    
    # 计算分母
    denominator <- sqrt(2 * (trace_cor_sq - (num_features^2) / degrees_of_freedom) * adjustment_factor)
    
    # 计算检验统计量
    test_statistic <- numerator / denominator
    
    # 计算单侧p值
    p_value <-pnorm(-test_statistic)
  
  # ======================== 整理结果并关闭集群 ========================
  computation_time <- Sys.time() - start_time
  return(list(
    p_value = p_value,  
    computation_time = computation_time
  ))
     # 关闭并行环境
  stopCluster(cl)
}
