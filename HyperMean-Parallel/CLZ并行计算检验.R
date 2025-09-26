CLZ_parallel_test <- function(sam1, sam2, eq.cov = F, eta = 0.05, num_cores=availableCores() - 1) {
  # Chen et al.(2014)的并行高维两样本检验
  # eq.cov: 是否假设协方差矩阵相等
  # eta: 阈值参数，默认0.05
  # num_cores: 并行计算使用的核心数
  
  # 设置并行运算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  start <- Sys.time()
  
  # 获取样本维度
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  p <- ncol(sam1)
  
  # ========== 1. 数据并行分块处理 ==========
  # 确定最优分块数（能被p整除的最大核心数）
  num_class <- num_cores
  while(p %% num_class != 0) {
    num_class <- num_class - 1
  }
  cols_per_group <- p / num_class  # 每块的变量数
  
  # 数据分块函数
  split_data <- function(sam) {
    lapply(seq_len(num_class), function(i) {
      start_col <- (i - 1) * cols_per_group + 1
      end_col <- i * cols_per_group
      sam[, start_col:end_col, drop = FALSE]  # 保持矩阵格式
    })
  }
  
  # 并行分块
  sam1_list <- split_data(sam1)
  sam2_list <- split_data(sam2)
  
  # 极值分布参数
  a_f <- sqrt(2*log(log(p)))
  b_f <- 2*log(log(p)) + 0.5*log(log(log(p))) - 0.5*log(4*pi/(1 - 0.05)^2)
  
  if(eq.cov) {  # 等协方差情况
    # ========== 2. 并行计算各数据块 ==========
    result_s_level <- foreach(i = 1:num_class, .combine = rbind) %dopar% {
      sam1_part <- sam1_list[[i]]
      sam2_part <- sam2_list[[i]]
      
      # 合并方差估计
      diag.cov <- ((colSums(sam1_part^2) - n1*colMeans(sam1_part)^2 + 
                      (colSums(sam2_part^2) - n2*colMeans(sam2_part)^2))) / (n1 + n2 - 2)
      diag.cov[diag.cov <= 1e-10] <- 1e-10  # 数值稳定
                   
      # 计算标准化统计量
      T_orig <- (colMeans(sam1_part) - colMeans(sam2_part))^2 / 
                ((1/n1 + 1/n2)*diag.cov*2*log(p))
                   
      # 阈值筛选
      s_level <- T_orig[T_orig > 0 & T_orig <= (1 - eta)]
                   
      # 对齐输出长度（便于合并结果）
      max_len <- ncol(sam1_part)
      list(T_orig = T_orig,
            s_level = c(s_level, rep(0, max_len - length(s_level))))
    }
    
  } else {  # 不等协方差情况
    result_s_level <- foreach(i = 1:num_class, .combine = rbind) %dopar% {
      sam1_part <- sam1_list[[i]]
      sam2_part <- sam2_list[[i]]
      
      # 分别估计方差
      diag1 <- (colSums(sam1_part^2) - n1*colMeans(sam1_part)^2) / (n1 - 1)
      diag2 <- (colSums(sam2_part^2) - n2*colMeans(sam2_part)^2) / (n2 - 1)
      diag1[diag1 <= 1e-10] <- 1e-10
      diag2[diag2 <= 1e-10] <- 1e-10
      
      # 计算标准化统计量
      T_orig <- (colMeans(sam1_part) - colMeans(sam2_part))^2 / 
        ((diag1/n1 + diag2/n2)*2*log(p))
      
      # 阈值筛选
      s_level <- T_orig[T_orig > 0 & T_orig <= (1 - eta)]
      
      max_len <- ncol(sam1_part)
      list(T_orig = T_orig,
           s_level = c(s_level, rep(0, max_len - length(s_level))))
    }
  }
  
  # ========== 3. 合并并行结果 ==========
  # 提取并合并所有块的T_orig和s_level
  T_orig_all <- unlist(lapply(result_s_level[,1], function(x) x))
  s_level_all <- unlist(lapply(result_s_level[,2], function(x) x[x > 0]))  # 过滤填充的0
  
  # ========== 4. 计算最终统计量 ==========
  s_m <- matrix(s_level_all*2*log(p), length(s_level_all), p, byrow = FALSE)
  T_m <- matrix(T_orig_all*2*log(p)-1, length(s_level_all), p, byrow = TRUE)
  T_m[T_m + 1 < s_m] <- 0  # 阈值筛选
  
  thr <- rowSums(T_m)
  mean_thr <- 2*sqrt(s_level_all*log(p)/pi)*p^(1-s_level_all)
  sd_thr <- sqrt(2/sqrt(2*pi)*p^(1-s_level_all)*((2*s_level_all*log(p))^1.5 + 
                                                   (2*s_level_all*log(p))^0.5) + 4*p - 4*p*pnorm(sqrt(2*s_level_all*log(p))))
  
  max_threshold <- max((thr - mean_thr)/sd_thr)*a_f - b_f
  pval <- 1 - exp(-exp(-as.numeric(max_threshold)))
  
  # ========== 5. 返回结果 ==========
  list(P值 = pval,
       time = Sys.time() - start)
}

