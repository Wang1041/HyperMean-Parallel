HLWP_test <- function(sam1, sam2) {
  # He & Xu & Wu & Pan (2021) 高维两样本检验方法
  # 结合有限范数统计量与极值统计量
  # 输入：sam1, sam2 - 样本矩阵（行=观测值，列=变量）
  
  start <- Sys.time()
  
  # 辅助函数：使用高效算法计算U统计量
  tmpU <- function(tmpV1, tmpV2, tmpV3, tmpV4, tmpV5, tmpV6) {
    npow <- 6  # 要计算的有限范数个数
    tmpRes <- matrix(0, npow, p)
    
    # 使用累积量高效计算中心矩
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
  
  # 辅助函数：排列数 P(n, a) = n!/(n-a)!
  permnum <- function(n, a) {
    factorial(n) / factorial(n - a)
  }
  
  # 初始化参数
  npow <- 6  # 有限范数个数（γ = 1,...,6）
  n1 <- nrow(sam1)
  n2 <- nrow(sam2)
  p <- ncol(sam1)
  n <- n1 + n2 - 2  # 合并协方差的自由度
  
  # 计算合并协方差矩阵
  Sn <- ((n1 - 1) * cov(sam1) + (n2 - 1) * cov(sam2)) / n
  
  # ========== 计算样本1的U统计量 ==========
  tmpV1xd <- colSums(sam1)       # 一阶矩
  tmpV2xd <- colSums(sam1^2)     # 二阶矩
  tmpV3xd <- colSums(sam1^3)     # 三阶矩
  tmpV4xd <- colSums(sam1^4)     # 四阶矩
  tmpV5xd <- colSums(sam1^5)     # 五阶矩
  tmpV6xd <- colSums(sam1^6)     # 六阶矩
  tmpResXd <- tmpU(tmpV1xd, tmpV2xd, tmpV3xd, tmpV4xd, tmpV5xd, tmpV6xd)
  
  # ========== 计算样本2的U统计量 ==========
  tmpV1yd <- colSums(sam2)
  tmpV2yd <- colSums(sam2^2)
  tmpV3yd <- colSums(sam2^3)
  tmpV4yd <- colSums(sam2^4)
  tmpV5yd <- colSums(sam2^5)
  tmpV6yd <- colSums(sam2^6)
  tmpResYd <- tmpU(tmpV1yd, tmpV2yd, tmpV3yd, tmpV4yd, tmpV5yd, tmpV6yd)
  
  # ========== 标准化并合并统计量 ==========
  # 创建第一行为1的矩阵以便高效计算
  tmpResXPd <- matrix(1, nrow = npow + 1, ncol = p)
  tmpResYPd <- matrix(1, nrow = npow + 1, ncol = p)
  
  # 用排列数进行标准化
  for (i in 2:(npow + 1)) {
    tmpResXPd[i, ] <- tmpResXd[i - 1, ] / permnum(n1, i - 1)
    tmpResYPd[i, ] <- tmpResYd[i - 1, ] / permnum(n2, i - 1)
  }
  
  # 合并两个样本的统计量
  tmpResbothd <- matrix(0, nrow = npow, ncol = p)
  for (i in 1:npow) {
    a <- i
    for (k in 0:a) {
      Cak <- choose(a, k)  # 二项式系数
      tmpResbothd[i, ] <- tmpResbothd[i, ] + 
        (-1)^(a - k) * Cak * (tmpResXPd[k + 1, ] * tmpResYPd[a - k + 1, ])
    }
  }
  
  # ========== 计算检验统计量 ==========
  tmpSPU <- rowSums(tmpResbothd)  # 对每个范数求变量和
  
  # 计算每个范数的渐近方差
  sigma <- numeric(6)
  for(a in 1:6) {
    sigma[a] <- factorial(a) * ((n1 + n2) / (n1 * n2))^a * sum(Sn^a)
  }
  
  # 计算有限范数的p值
  Ts <- tmpSPU / sqrt(sigma)
  pval.f <- 2 * (1 - pnorm(abs(Ts)))  # 双边检验
  
  # ========== 计算极值统计量 (Cai et al. 2014方法) ==========
  x1_bar <- colMeans(sam1)
  x2_bar <- colMeans(sam2)
  diff <- x1_bar - x2_bar
  nume <- diff^2
  
  # 方差估计
  if(n1 == n2) {
    deno <- diag(Sn) / n1 * 2  # 等样本量时的方差估计
  } else {
    deno <- (colMeans(sam1^2) - x1_bar^2) / (n1 - 1) +  # 样本1方差
      (colMeans(sam2^2) - x2_bar^2) / (n2 - 1)    # 样本2方差
  }
  deno[deno <= 1e-10] <- 1e-10  # 数值稳定性处理
  
  T_CLX <- max(nume / deno)  # 极值统计量
  M <- T_CLX - 2 * log(p) + log(log(p)) + log(pi)  # 极值分布调整
  pval_inf <- 1 - exp(-exp(-M / 2))  # 基于Gumbel分布的p值
  
  # ========== 合并p值 ==========
  pval.f <- c(pval.f, pval_inf)  # 共7个p值（6个有限范数 + 1个极值）
  pval.min <- min(pval.f)         # 取最小p值
  pval <- 1 - (1 - pval.min)^7   # Bonferroni类型调整
  
  # ========== 返回结果 ==========
  out <- list()
  out$pval.f <- pval.f  # 各范数的单独p值
  out$p值 <- pval       # 最终合并的p值
  out$time <- Sys.time() - start  # 计算耗时
  return(out)
}
HLWP_test(sam1,sam2)

