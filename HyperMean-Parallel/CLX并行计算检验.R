CLX_parallel_test <- function(sam1, sam2, num_cores=availableCores() - 1) {
  # Cai et al.(2014)的并行高维两样本均值检验
  # 输入：
  #   sam1, sam2 - 样本数据矩阵（行是观测，列是变量）
  #   num_cores - 使用的CPU核心数
  # 输出：包含P值和计算时间的列表
  
  
  # 设置并行运算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  start <- Sys.time()  # 记录开始时间
  
  # 1. 获取样本维度信息
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  p <- ncol(sam1)  # 变量维度
  
  # 2. 数据并行分块处理
  num_class <- num_cores  # 分块数等于核心数
  cols_per_group <- floor(p / num_class)  # 每块基本列数
  remainder <- p %% num_class  # 余数列数
  
  # 数据分块函数（处理不能整除的情况）
  split_data <- function(sam) {
    lapply(seq_len(num_class), function(i) {
      if(i <= remainder) {
        # 前remainder块每块多分1列
        start_col <- (i - 1) * (cols_per_group + 1) + 1
        end_col <- i * (cols_per_group + 1)
      } else {
        # 剩余块按基本列数分配
        start_col <- (i - 1) * cols_per_group + remainder + 1
        end_col <- min(i * cols_per_group + remainder, p)
      }
      sam[, start_col:end_col, drop = FALSE]  # 保持矩阵结构
    })
  }
  
  # 数据分块
  sam1_list <- split_data(sam1)
  sam2_list <- split_data(sam2)
  
  # 3. 并行计算各数据块的统计量
  P_cai_zhi <- foreach(i = 1:num_class, .combine = rbind) %dopar% {
    sam1_part <- sam1_list[[i]]
    sam2_part <- sam2_list[[i]]
    
    # 计算当前分块的统计量
    x1_bar <- colMeans(sam1_part)
    x2_bar <- colMeans(sam2_part)
    x_cha <- x1_bar - x2_bar
    nume <- x_cha^2
    deno <- (colSums(sam1_part^2)/n1 - x1_bar^2)/(n1-1) + 
      (colSums(sam2_part^2)/n2 - x2_bar^2)/(n2-1)
    deno[deno <= 1e-10] <- 1e-10  # 数值稳定性处理
    
    max(nume/deno)  # 返回当前分块的最大值
  }
  
  # 4. 合并结果并计算最终统计量
  T_CLX_result <- max(P_cai_zhi)  # 所有分块中的最大值
  M <- T_CLX_result - 2*log(p) + log(log(p)) + log(pi)  # 极值分布调整
  alpha <- 1 - exp(-exp(-M/2))  # p值计算
  
  # 5. 返回结果
  out <- list()
  out$P值 <- alpha  # 检验p值
  out$time <- Sys.time() - start  # 计算耗时
  # 关闭并行环境
  stopCluster(cl)
  
  return(out)
}
CLX_parallel_test (sam1,sam2,3)
