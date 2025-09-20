CLZ_test <- function(sam1, sam2, eq.cov = F, eta = 0.05) {
  # Chen et al.(2014)的高维两样本检验
  # eq.cov: 是否假设协方差矩阵相等
  # eta: 阈值参数，默认0.05
  
  if(eq.cov) {  # 等协方差情况
    start <- Sys.time()
    n1 <- dim(sam1)[1]; n2 <- dim(sam2)[1]; p <- dim(sam1)[2]
    
    # 极值分布参数
    a_f <- sqrt(2*log(log(p)))
    b_f <- 2*log(log(p)) + 0.5*log(log(log(p))) - 0.5*log(4*pi/(1 - 0.05)^2)
    
    # 合并方差估计
    diag.cov <- ((colSums(sam1^2)-n1*colMeans(sam1)^2) + 
                   (colSums(sam2^2)-n2*colMeans(sam2)^2))/(n1 + n2 - 2)
    diag.cov[diag.cov <= 10^(-10)] <- 10^(-10)  # 防止数值问题
    
    # 标准化检验统计量
    T_orig <- (colMeans(sam1) - colMeans(sam2))^2/((1/n1 + 1/n2)*diag.cov*2*log(p))
    
    # 阈值筛选过程
    s_level <- T_orig[sign(T_orig) > 0 & T_orig <= (1 - eta)] 
    s_m <- matrix(s_level*2*log(p), length(s_level), p, byrow = F)
    T_m <- matrix(T_orig*2*log(p)-1, length(s_level), p, byrow = T)
    T_m[T_m + 1 < s_m ] <- 0  # 阈值筛选
    
    # 极值统计量计算
    thr <- rowSums(T_m) 
    mean_thr <- 2*sqrt(s_level*log(p)/pi)*p^(1-s_level)
    sd_thr <- sqrt(2/sqrt(2*pi)*p^(1-s_level)*((2*s_level*log(p))^1.5+(2*s_level*log(p))^0.5) +
                     4*p - 4*p*pnorm(sqrt(2*s_level*log(p))))
    max_threshold <- max((thr-mean_thr)/sd_thr)*a_f - b_f
    
    # p值计算（极值分布）
    pval <- 1 - exp(-exp(-as.numeric(max_threshold)))
    
    out <- list(pval = pval, time = Sys.time()-start)
    return(out)
    
  } else {  # 不等协方差情况
    start <- Sys.time()
    n1 <- dim(sam1)[1]; n2 <- dim(sam2)[1]; p <- dim(sam1)[2]
    
    # 极值分布参数（同上）
    a_f <- sqrt(2*log(log(p)))
    b_f <- 2*log(log(p)) + 0.5*log(log(log(p))) - 0.5*log(4*pi/(1 - 0.05)^2)
    
    # 分别估计方差
    diag1 <- (colSums(sam1^2)-n1*colMeans(sam1)^2)/(n1-1)
    diag2 <- (colSums(sam2^2)-n2*colMeans(sam2)^2)/(n2-1)
    diag1[diag1 <= 10^(-10)] <- 10^(-10)
    diag2[diag2 <= 10^(-10)] <- 10^(-10)
    
    # 标准化检验统计量
    T_orig <- (colMeans(sam1) - colMeans(sam2))^2/((diag1/n1 + diag2/n2)*2*log(p))
    
    # 以下步骤与等协方差情况相同
    s_level <- T_orig[sign(T_orig) > 0 & T_orig <= (1 - eta)] 
    s_m <- matrix(s_level*2*log(p), length(s_level), p, byrow = F)
    T_m <- matrix(T_orig*2*log(p)-1, length(s_level), p, byrow = T)
    T_m[T_m + 1 < s_m ] <- 0
    thr <- rowSums(T_m) 
    mean_thr <- 2*sqrt(s_level*log(p)/pi)*p^(1-s_level)
    sd_thr <- sqrt(2/sqrt(2*pi)*p^(1-s_level)*((2*s_level*log(p))^1.5+(2*s_level*log(p))^0.5) +
                     4*p - 4*p*pnorm(sqrt(2*s_level*log(p))))
    max_threshold <- max((thr-mean_thr)/sd_thr)*a_f - b_f
    
    pval <- 1 - exp(-exp(-as.numeric(max_threshold)))
    
    out <- list(pval = pval, time = Sys.time()-start)
    return(out)
  }
}

CLZ_test(sam1,sam2)
