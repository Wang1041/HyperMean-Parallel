HXWP_parallel_test  <- function(sam1, sam2, num_cores=availableCores() - 1) {
  # He et al.(2021) 并行高维两样本检验
  # 使用并行计算加速统计量计算
  # 输入：sam1, sam2 - 样本矩阵，num_cores - 并行核心数
  # 设置并行运算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  start <- Sys.time()

  # 辅助函数：计算U统计量（高效算法）
  tmpU <- function(tmpV1, tmpV2, tmpV3, tmpV4, tmpV5, tmpV6, p) {
    npow <- 6  # 有限范数个数
    tmpRes <- matrix(0, npow, p)
    
    # 使用累积量计算中心矩
    tmpRes[1,] <- tmpV1
    tmpRes[2,] <- tmpV1^2 - tmpV2
    tmpRes[3,] <- tmpV1^3 - 3 * tmpV2 * tmpV1 + 2 * tmpV3
    tmpRes[4,] <- tmpV1^4 - 6 * tmpV2 * tmpV1^2 + 3 * tmpV2^2 + 8 * tmpV3 * tmpV1 - 6 * tmpV4
    tmpRes[5,] <- tmpV1^5 - 10 * tmpV2 * tmpV1^3 + 15 * tmpV2^2 * tmpV1 + 20 * tmpV3 * tmpV1^2 - 
      20 * tmpV3 * tmpV2 - 30 * tmpV4 * tmpV1 + 24 * tmpV5
    tmpRes[6,] <- tmpV1^6 - 15 * tmpV2 * tmpV1^4 + 40 * tmpV3 * tmpV1^3 + 45 * tmpV1^2 * tmpV2^2 - 
      90 * tmpV1^2 * tmpV4 - 120 * tmpV1 * tmpV2 * tmpV3 + 144 * tmpV1 * tmpV5 - 
      15 * tmpV2^3 + 90 * tmpV2 * tmpV4 + 40 * tmpV3^2 - 120 * tmpV6
    
    return(tmpRes)
  }
  
  # 辅助函数：计算排列数 P(n, a) = n!/(n-a)!
  permnum <- function(n, a) {
    factorial(n) / factorial(n - a)
  }
  
  # 获取样本基本信息
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  p <- ncol(sam1)
  npow <- 6  # 有限范数个数
  
  # 计算合并协方差矩阵
  n <- n1 + n2 - 2
  Sn <- ((n1 - 1) * cov(sam1) + (n2 - 1) * cov(sam2)) / n
  
  # ========== 并行计算排列数 ==========
  pxy <- foreach(i = 1:npow, .combine = rbind) %dopar% {
    px <- permnum(n1, i)
    py <- permnum(n2, i)
    c(px, py)  # 返回两个样本的排列数
  }
  
  # ========== 数据分块策略 ==========
  # 将变量维度p均匀分配到各个核心
  base_size <- p %/% num_cores  # 每个核心的基本块大小
  remainder <- p %% num_cores   # 剩余变量数
  
  idx <- rep(base_size, num_cores)
  if(remainder != 0) {
    idx[1:remainder] <- idx[1:remainder] + 1  # 前remainder个核心多分1个变量
  }
  
  # 计算累积索引用于数据分割
  cumulative_idx <- cumsum(c(0, idx))
  
  # ========== 并行计算样本1的统计量 ==========
  xresult_UV <- foreach(i = 1:num_cores, .combine = c) %dopar% {
    # 提取当前核心负责的数据块
    crossXtmpd <- sam1[, (cumulative_idx[i] + 1):cumulative_idx[i + 1], drop = FALSE]
    
    # 计算1-6阶矩
    tmpV1xd <- colSums(crossXtmpd)
    tmpV2xd <- colSums(crossXtmpd^2)
    tmpV3xd <- colSums(crossXtmpd^3)
    tmpV4xd <- colSums(crossXtmpd^4)
    tmpV5xd <- colSums(crossXtmpd^5)
    tmpV6xd <- colSums(crossXtmpd^6)
    
    # 计算U统计量并标准化
    tmpResXd <- tmpU(tmpV1xd, tmpV2xd, tmpV3xd, tmpV4xd, tmpV5xd, tmpV6xd, idx[i]) / pxy[, 1]
    list(tmpResXd)  # 返回当前块的结果
  }
  
  # ========== 并行计算样本2的统计量 ==========
  yresult_UV <- foreach(i = 1:num_cores, .combine = c) %dopar% {
    # 提取当前核心负责的数据块
    crossYtmpd <- sam2[, (cumulative_idx[i] + 1):cumulative_idx[i + 1], drop = FALSE]
    
    # 计算1-6阶矩
    tmpV1yd <- colSums(crossYtmpd)
    tmpV2yd <- colSums(crossYtmpd^2)
    tmpV3yd <- colSums(crossYtmpd^3)
    tmpV4yd <- colSums(crossYtmpd^4)
    tmpV5yd <- colSums(crossYtmpd^5)
    tmpV6yd <- colSums(crossYtmpd^6)
    
    # 计算U统计量并标准化
    tmpResYd <- tmpU(tmpV1yd, tmpV2yd, tmpV3yd, tmpV4yd, tmpV5yd, tmpV6yd, idx[i]) / pxy[, 2]
    list(tmpResYd)  # 返回当前块的结果
  }
  
  # ========== 合并并行计算结果 ==========
  # 组合两个样本的统计量矩阵
  result_UV <- rbind(1, do.call(cbind, xresult_UV), 1, do.call(cbind, yresult_UV))
  
  # ========== 并行计算检验统计量和方差 ==========
  tmpSPU <- foreach(i = 1:npow, .combine = rbind) %dopar% {
    tmpResbothd <- 0
    # 组合两个样本的统计量
    for (k in 0:i) {
      Cak <- choose(i, k)  # 二项式系数
      tmpResbothd <- tmpResbothd + (-1)^(i - k) * Cak * 
        (result_UV[k + 1, ] * result_UV[i - k + 8, ])
    }
    
    # 计算统计量和标准差
    tmpSPU_value <- sum(tmpResbothd)
    sigma <- sqrt(factorial(i) * ((n1 + n2) / (n1 * n2))^i * sum(Sn^i))
    
    c(tmpSPU_value, sigma)  # 返回统计量和标准差
  }
  
  # ========== 计算有限范数的p值 ==========
  Ts <- tmpSPU[, 1] / tmpSPU[, 2]  # 标准化统计量
  pval.f <- 2 * (1 - pnorm(abs(Ts)))  # 双边检验p值
  
  # ========== 计算极值统计量 (Cai et al. 2014) ==========
  x1_bar <- colMeans(sam1)
  x2_bar <- colMeans(sam2)
  diff <- x1_bar - x2_bar
  nume <- diff^2
  
  # 方差估计
  if(n1 == n2) {
    deno <- diag(Sn) / n1 * 2
  } else {
    deno <- (colMeans(sam1^2) - x1_bar^2) / (n1 - 1) + 
      (colMeans(sam2^2) - x2_bar^2) / (n2 - 1)
  }
  deno[deno <= 1e-10] <- 1e-10  # 数值稳定性
  
  T_CLX <- max(nume / deno)  # 极值统计量
  M <- T_CLX - 2 * log(p) + log(log(p)) + log(pi)
  pval_inf <- 1 - exp(-exp(-M / 2))  # 极值分布p值
  
  # ========== 合并p值 ==========
  pval.f <- c(pval.f, pval_inf)  # 7个p值（6有限+1极值）
  pval.min <- min(pval.f)
  pval <- 1 - (1 - pval.min)^7  # Bonferroni调整
  
  # ========== 返回结果 ==========
  out <- list()
  out$pval.f <- pval.f  # 各范数单独p值
  out$p值 <- pval       # 最终合并p值
  out$time <- Sys.time() - start  # 计算时间
  # 关闭并行环境
  stopCluster(cl)
  return(out)
}


