CQ_parallel_test <- function(sam1, sam2, num_cores) {
  # 设置并行运算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # 初始化计时和并行计算环境
  start <- Sys.time()
  
  # 获取样本维度信息
  n1 <- dim(sam1)[1]  # 样本1的观测数
  n2 <- dim(sam2)[1]  # 样本2的观测数
  p <- dim(sam1)[2]   # 变量维度
  
  # 计算列总和
  x1_bar <- colSums(sam1)
  x2_bar <- colSums(sam2)
  
  # 并行计算6个核心统计量
  results <- foreach(i = 1:6, .combine = 'rbind') %dopar% {
    if (i == 1) {
      # 1. 计算tr(Σ1²)的U统计量
      Sigma1_square <- 0
      for (j in 1:n1) {
        sam1j <- sam1[j, ] 
        for(k in 1:n1) {
          if (j != k) {
            sam1k <- sam1[k, ]
            tempmean1 <- (x1_bar - sam1j - sam1k) / (n1 - 2)
            Sigma1_square <- Sigma1_square + 
              sum(sam1j * (sam1k - tempmean1)) * 
              sum(sam1k * (sam1j - tempmean1))
          }
        }
      }
      return(c(Sigma1_square / (n1 * (n1 - 1)), NA, NA, NA, NA, NA))
      
    } else if(i == 2) {
      # 2. 计算tr(Σ2²)
      Sigma2_square <- 0
      for (j in 1:n2) {
        for(k in 1:n2) {
          if (j != k) {
            sam2j <- sam2[j, ]; sam2k <- sam2[k, ]
            tempmean2 <- (x2_bar - sam2j - sam2k) / (n2 - 2)
            Sigma2_square <- Sigma2_square + 
              sum(sam2j * (sam2k - tempmean2)) * 
              sum(sam2k * (sam2j - tempmean2))
          }
        }
      }
      return(c(NA, Sigma2_square / (n2 * (n2 - 1)), NA, NA, NA, NA))
      
    } else if(i == 3) {
      # 3. 计算tr(Σ1Σ2)的交叉项
      Sigma1_2 <- 0
      for (j in 1:n1) {
        sam1j <- sam1[j, ] 
        for (k in 1:n2) {
          sam2k <- sam2[k, ]
          Sigma1_2 <- Sigma1_2 +
            sum(sam1j * (sam2k - (x2_bar - sam2k) / (n2 - 1))) *
            sum(sam2k * (sam1j - (x1_bar - sam1j) / (n1 - 1)))
        }
      }
      return(c(NA, NA, Sigma1_2 / (n1 * n2), NA, NA, NA))
      
    } else if(i == 4) {
      # 4. 计算样本1的内积和（排除对角线）
      T1 <- sam1 %*% t(sam1)
      return(c(NA, NA, NA, sum(T1) - sum(diag(T1)), NA, NA))
      
    } else if(i == 5) {
      # 5. 计算样本2的内积和（排除对角线）
      T2 <- sam2 %*% t(sam2)
      return(c(NA, NA, NA, NA, sum(T2) - sum(diag(T2)), NA))
      
    } else if(i == 6) {
      # 6. 计算样本间交叉内积和
      return(c(NA, NA, NA, NA, NA, sum(sam1 %*% t(sam2))))
    }
  }
  
  # 组装最终结果
  trSigma1_square <- results[1, 1]
  trSigma2_square <- results[2, 2]
  tr_Sigma1_Sigma2 <- results[3, 3]
  
  # 计算检验统计量
  nume <- results[4, 4] / (n1 * (n1 - 1)) + 
    results[5, 5] / (n2 * (n2 - 1)) - 
    2 * results[6, 6] / (n1 * n2)
  
  deno <- 2 / (n1 * (n1 - 1)) * trSigma1_square + 
    2 / (n2 * (n2 - 1)) * trSigma2_square + 
    4 / (n1 * n2) * tr_Sigma1_Sigma2
  # 返回结果列表
  list(
    p值 = pnorm(-nume / sqrt(deno)),  # 正态近似p值
    time = Sys.time() - start         # 总计算时间
  )
    # 关闭并行环境
  stopCluster(cl)
}

