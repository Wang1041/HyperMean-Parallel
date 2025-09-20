XLWP_parallel_else <- function(sam1, sam2, eq.cov = TRUE, G_greek = c(1:6, Inf), 
                                 num_cores=availableCores() - 1, cv.fold = 5, norm = "F", 
                                 cov.same, cov1, cov2, bandwidth, bandwidth1, bandwidth2) {
  # Xu et al.(2016) 并行高维两样本检验
  # 使用并行计算加速SPU检验和带宽选择
  # 输入参数：
  #   sam1, sam2: 样本矩阵（行=观测，列=变量）
  #   eq.cov: 是否假设协方差矩阵相等
  #   G_greek: γ值集合，默认1-6和无穷大
  #   num_cores: 并行核心数
  #   cv.fold: 交叉验证折数
  #   norm: 矩阵范数类型
  #   cov.same/cov1/cov2: 预指定的协方差矩阵
  #   bandwidth/bandwidth1/bandwidth2: 带宽参数
  # 设置并行运算环境
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  start <- Sys.time()
  
  # 并行带宽选择函数
  best.band_parallel <- function(sam, bandwidth, cv.fold, norm) {
    p <- ncol(sam)
    n <- nrow(sam)
    fold.size <- round(n / cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    matrix.id <- abs(row(diag(p)) - col(diag(p)))  # 距离矩阵
    
    # 并行交叉验证
    diff.norm <- foreach(i = 1:cv.fold, .combine = rbind) %dopar% {
      if(i == cv.fold) {
        temp.idx <- sam.idx[((i - 1) * fold.size + 1):n]
      } else {
        temp.idx <- sam.idx[((i - 1) * fold.size + 1):(i * fold.size)]
      }
      
      sam.train <- sam[-temp.idx, ]
      sam.test <- sam[temp.idx, ]
      sam.train.cov <- cov(sam.train)
      sam.test.cov <- cov(sam.test)
      
      diff <- numeric(n.bandwidth)
      for(j in 1:n.bandwidth) {
        sam.train.cov.band <- sam.train.cov
        sam.train.cov.band[matrix.id > bandwidth[j]] <- 0  # 带宽截断
        diff[j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
      }
      diff  # 返回当前fold的差异
    }
    
    # 选择最优带宽
    diff.norm <- colMeans(diff.norm)
    best <- which.min(diff.norm)
    return(bandwidth[best[1]])
  }
  
  # 并行等协方差检验函数
  Xu_same_parallel <- function(sam1, sam2, cov.same, G_greek = c(1:6, Inf)) {
    n1 <- nrow(sam1); n2 <- nrow(sam2); p <- ncol(sam1)
    x1_bar <- colMeans(sam1)
    x2_bar <- colMeans(sam2)
    diff <- x1_bar - x2_bar
    
    # 确保包含无穷大γ值
    if(!is.element(Inf, G_greek)) {
      G_greek <- c(G_greek, Inf)
    }
    G_greek0 <- G_greek[G_greek != Inf]
    
    # 划分奇偶γ值集合
    odd.ga <- G_greek0[G_greek0 %% 2 == 1]  # 奇数γ
    even.ga <- G_greek0[G_greek0 %% 2 == 0] # 偶数γ
    n.odd.ga <- length(odd.ga)
    n.even.ga <- length(even.ga)
    
    # γ=1时的方差
    var1 <- (1/n1 + 1/n2) * sum(cov.same)
    
    # 并行计算SPU统计量的均值（偶数γ）
    result_P0 <- foreach(gamma = 2*(1:6), .combine = rbind) %dopar% {
      gamma.half <- gamma / 2
      u_i <- 0
      for(d in 0:gamma.half) {
        u_i <- u_i + 1 / (factorial(d) * factorial(gamma.half - d) * n1^d * n2^(gamma.half - d))
      }
      u <- (diag(cov.same))^gamma.half * factorial(gamma) * u_i / 2^gamma.half 
      return(u)
    }
    
    # 并行生成Greek_A索引矩阵（方差计算）
    Greek_A <- foreach(i = 2:6, .combine = rbind) %dopar% {
      result_Greek_A <- matrix(NA, 0, 6)
      t <- s <- i
      for(c3 in 0:min(c(t, s))) {
        for(d3 in 0:min(c(t - c3, s - c3))) {
          if(c3 + d3 > 0) {
            for(c1 in 0:floor((t - c3 - d3) / 2)) {
              for(c2 in 0:floor((s - c3 - d3) / 2)) {
                if((t - 2*c1 - c3 - d3) %% 2 == 0 & (s - 2*c2 - c3 - d3) %% 2 == 0) {
                  d1 <- (t - 2*c1 - c3 - d3) / 2
                  d2 <- (s - 2*c2 - c3 - d3) / 2
                  result_Greek_A <- rbind(result_Greek_A, c(c1, c2, c3, d1, d2, d3))
                }
              }
            }
          }
        }
      }
      return(result_Greek_A)
    }
    
    # 并行生成Greek_Ast索引矩阵（协方差计算）
    Greek_Ast <- foreach(t = 2:6, .combine = rbind) %dopar% {
      result_Greek_A <- matrix(NA, 0, 6)
      for(s in 1:(t-1)) {
        for(c3 in 0:s) {
          for(d3 in 0:(s - c3)) {
            if(c3 + d3 > 0) {
              for(c1 in 0:floor((s - c3 - d3) / 2)) {
                for(c2 in 0:floor((t - c3 - d3) / 2)) {
                  if((t - 2*c2 - c3 - d3) %% 2 == 0 & (s - 2*c1 - c3 - d3) %% 2 == 0) {
                    d1 <- (s - 2*c1 - c3 - d3) / 2
                    d2 <- (t - 2*c2 - c3 - d3) / 2
                    result_Greek_A <- rbind(result_Greek_A, c(c1, c2, c3, d1, d2, d3))
                  }
                }
              }
            }
          }
        }
      }
      return(result_Greek_A)
    }
    
    # 准备矩阵计算
    diag_cov <- diag(cov.same)
    p <- length(diag_cov)
    mat1 <- matrix(rep(diag_cov, p), p, p, byrow = FALSE)
    mat2 <- matrix(rep(diag_cov, p), p, p, byrow = TRUE)
    diag(cov.same) <- 0
    
    # 并行计算方差第三项P3
    result_P3 <- foreach(i = 1:nrow(Greek_A), .combine = rbind) %dopar% {
      c1 <- Greek_A[i, 1]; c2 <- Greek_A[i, 2]; c3 <- Greek_A[i, 3]
      d1 <- Greek_A[i, 4]; d2 <- Greek_A[i, 5]; d3 <- Greek_A[i, 6]
      
      mat <- mat1^(c1 + d1) * mat2^(c2 + d2) * cov.same^(c3 + d3)
      nume <- sum(mat)
      deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * factorial(c1) * factorial(c2) *
                 factorial(d1) * factorial(d2) * factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
      nume / deno
    }
    
    # 方差计算（各γ值）
    P3 <- c(4 * sum(result_P3[1:3]), 36 * sum(result_P3[4:15]),
            576 * sum(result_P3[16:32]), 14400 * sum(result_P3[33:72]),
            518400 * sum(result_P3[73:126]))
    
    u_sum <- rowSums(result_P0)
    P1 <- u_sum[-1]
    P2 <- c(sum(result_P0[1, ]^2), 0, sum(result_P0[2, ]^2), 0, sum(result_P0[3, ]^2))
    var <- c(var1, P1 - P2 + P3)
    
    # 计算检验统计量
    XLWP <- (colSums(sapply(1:6, function(i) diff^i)) - c(0, u_sum[1], 0, u_sum[2], 0, u_sum[3])) / sqrt(var)
    T_O <- max(abs(XLWP[c(1, 3, 5)]))  # 奇数γ的最大值
    T_E <- max(XLWP[c(2, 4, 6)])       # 偶数γ的最大值
    
    # 并行计算协方差第三项
    cov_result_P3 <- foreach(i = 1:nrow(Greek_Ast), .combine = rbind) %dopar% {
      c1 <- Greek_Ast[i, 1]; c2 <- Greek_Ast[i, 2]; c3 <- Greek_Ast[i, 3]
      d1 <- Greek_Ast[i, 4]; d2 <- Greek_Ast[i, 5]; d3 <- Greek_Ast[i, 6]
      
      mat <- mat1^(c1 + d1) * mat2^(c2 + d2) * cov.same^(c3 + d3)
      nume <- sum(mat)
      deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * factorial(c1) * factorial(c2) *
                 factorial(d1) * factorial(d2) * factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
      nume / deno
    }
    
    # 协方差计算
    cov_P3 <- c(6 * sum(cov_result_P3[1:4]), 48 * sum(cov_result_P3[5:10]),
                120 * sum(cov_result_P3[11:16]), 720 * sum(cov_result_P3[17:36]),
                1440 * sum(cov_result_P3[37:45]), 17280 * sum(cov_result_P3[46:73]))
    
    cov_P1 <- c(u_sum[2], u_sum[3], u_sum[3], u_sum[4], u_sum[4], u_sum[5])
    cov_P2 <- c(0, sum(result_P0[1, ] * result_P0[2, ]), 0, 0, 
                sum(result_P0[1, ] * result_P0[3, ]), sum(result_P0[2, ] * result_P0[3, ]))
    cov <- cov_P1 - cov_P2 + cov_P3
    
    # 计算相关矩阵
    cov_e <- cov[c(2, 5, 6)] / sqrt(c(var[2] * var[4], var[2] * var[6], var[4] * var[6]))
    cov_o <- cov[c(1, 3, 4)] / sqrt(c(var[1] * var[3], var[1] * var[5], var[5] * var[3]))
    
    R_E <- matrix(c(1, cov_e[1], cov_e[2], cov_e[1], 1, cov_e[3], cov_e[2], cov_e[3], 1), 3, 3)
    R_O <- matrix(c(1, cov_o[1], cov_o[2], cov_o[1], 1, cov_o[3], cov_o[2], cov_o[3], 1), 3, 3)
    
    # 计算p值（多元正态分布）
    pval_O <- 1 - pmvnorm(lower = -rep(T_O, 3), upper = rep(T_O, 3), mean = rep(0, 3), sigma = R_O)
    pval_E <- 1 - pmvnorm(lower = rep(-Inf, 3), upper = rep(T_E, 3), mean = rep(0, 3), sigma = R_E)
    
    # 计算γ=∞的情况（Cai et al. 2014）
    nume <- diff^2
    if(n1 == n2) {
      deno <-  diag_cov  / n1 * 2
    } else {
      deno <- (colMeans(sam1^2) - x1_bar^2) / (n1 - 1) + (colMeans(sam2^2) - x2_bar^2) / (n2 - 1)
    }
    deno[deno <= 1e-10] <- 1e-10
    
    T_CLX <- max(nume / deno)
    M <- T_CLX - 2 * log(p) + log(log(p)) + log(pi)
    pval_inf <- 1 - exp(-exp(-M / 2))
    
    # 最终p值（Bonferroni调整）
    pval.min <- min(c(pval_O, pval_E, pval_inf))
    pval <- 1 - (1 - pval.min)^3
    
    return(list(P_value = pval))
  }
  
  # 不等协方差情况下的并行检验函数
  Xu_diff_parallel <- function(sam1, sam2, cov1, cov2, G_greek = c(1:6, Inf)) {
    n1 <- nrow(sam1); n2 <- nrow(sam2); p <- ncol(sam1)
    x1_bar <- colMeans(sam1)
    x2_bar <- colMeans(sam2)
    diff <- x1_bar - x2_bar
    
    # 确保包含无穷大γ值
    if(!is.element(Inf, G_greek)) {
      G_greek <- c(G_greek, Inf)
    }
    G_greek0 <- G_greek[G_greek != Inf]
    
    # 划分奇偶γ值集合
    odd.ga <- G_greek0[G_greek0 %% 2 == 1]  # 奇数γ
    even.ga <- G_greek0[G_greek0 %% 2 == 0] # 偶数γ
    n.odd.ga <- length(odd.ga)
    n.even.ga <- length(even.ga)
    
    # γ=1时的方差（不等协方差版本）
    var1 <- sum(cov1)/n1 + sum(cov2)/n2
    
    # 并行计算SPU统计量的均值（偶数γ）- 不等协方差版本
    result_P0 <- foreach(gamma = 2*(1:6), .combine = rbind) %dopar% {
      gamma.half <- gamma / 2
      u_i <- 0
      for(d in 0:gamma.half) {
        u_i <- u_i + 1 / (factorial(d) * factorial(gamma.half - d) * n1^d * n2^(gamma.half - d))
      }
      # 不等协方差：分别使用cov1和cov2的对角线元素
      u <- (diag(cov1)^d) * (diag(cov2)^(gamma.half - d)) * factorial(gamma) * u_i / 2^gamma.half 
      return(u)
    }
    
    # 并行生成Greek_A索引矩阵（方差计算）- 与等协方差相同结构
    Greek_A <- foreach(i = 2:6, .combine = rbind) %dopar% {
      result_Greek_A <- matrix(NA, 0, 6)
      t <- s <- i
      for(c3 in 0:min(c(t, s))) {
        for(d3 in 0:min(c(t - c3, s - c3))) {
          if(c3 + d3 > 0) {
            for(c1 in 0:floor((t - c3 - d3) / 2)) {
              for(c2 in 0:floor((s - c3 - d3) / 2)) {
                if((t - 2*c1 - c3 - d3) %% 2 == 0 & (s - 2*c2 - c3 - d3) %% 2 == 0) {
                  d1 <- (t - 2*c1 - c3 - d3) / 2
                  d2 <- (s - 2*c2 - c3 - d3) / 2
                  result_Greek_A <- rbind(result_Greek_A, c(c1, c2, c3, d1, d2, d3))
                }
              }
            }
          }
        }
      }
      return(result_Greek_A)
    }
    
    # 并行生成Greek_Ast索引矩阵（协方差计算）
    Greek_Ast <- foreach(t = 2:6, .combine = rbind) %dopar% {
      result_Greek_A <- matrix(NA, 0, 6)
      for(s in 1:(t-1)) {
        for(c3 in 0:s) {
          for(d3 in 0:(s - c3)) {
            if(c3 + d3 > 0) {
              for(c1 in 0:floor((s - c3 - d3) / 2)) {
                for(c2 in 0:floor((t - c3 - d3) / 2)) {
                  if((t - 2*c2 - c3 - d3) %% 2 == 0 & (s - 2*c1 - c3 - d3) %% 2 == 0) {
                    d1 <- (s - 2*c1 - c3 - d3) / 2
                    d2 <- (t - 2*c2 - c3 - d3) / 2
                    result_Greek_A <- rbind(result_Greek_A, c(c1, c2, c3, d1, d2, d3))
                  }
                }
              }
            }
          }
        }
      }
      return(result_Greek_A)
    }
    
    # 准备矩阵计算 - 不等协方差版本
    diag_cov1 <- diag(cov1)
    diag_cov2 <- diag(cov2)
    p <- length(diag_cov1)
    
    # 创建四个矩阵分别对应两个协方差矩阵的行列组合
    mat1_cov1 <- matrix(rep(diag_cov1, p), p, p, byrow = FALSE)
    mat2_cov1 <- matrix(rep(diag_cov1, p), p, p, byrow = TRUE)
    mat1_cov2 <- matrix(rep(diag_cov2, p), p, p, byrow = FALSE)
    mat2_cov2 <- matrix(rep(diag_cov2, p), p, p, byrow = TRUE)
    
    diag(cov1) <- diag(cov2) <- 0  # 剔除对角线
    
    # 并行计算方差第三项P3 - 不等协方差版本
    result_P3 <- foreach(i = 1:nrow(Greek_A), .combine = rbind) %dopar% {
      c1 <- Greek_A[i, 1]; c2 <- Greek_A[i, 2]; c3 <- Greek_A[i, 3]
      d1 <- Greek_A[i, 4]; d2 <- Greek_A[i, 5]; d3 <- Greek_A[i, 6]
      
      # 不等协方差：使用两个协方差矩阵的组合
      mat <- mat1_cov1^c1 * mat1_cov2^d1 * mat2_cov1^c2 * mat2_cov2^d2 * cov1^c3 * cov2^d3
      nume <- sum(mat)
      deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * factorial(c1) * factorial(c2) *
                 factorial(d1) * factorial(d2) * factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
      nume / deno
    }
    
    # 方差计算（各γ值）
    P3 <- c(4 * sum(result_P3[1:3]), 36 * sum(result_P3[4:15]),
            576 * sum(result_P3[16:32]), 14400 * sum(result_P3[33:72]),
            518400 * sum(result_P3[73:126]))
    
    u_sum <- rowSums(result_P0)
    P1 <- u_sum[-1]
    P2 <- c(sum(result_P0[1, ]^2), 0, sum(result_P0[2, ]^2), 0, sum(result_P0[3, ]^2))
    var <- c(var1, P1 - P2 + P3)
    
    # 计算检验统计量
    XLWP <- (colSums(sapply(1:6, function(i) diff^i)) - c(0, u_sum[1], 0, u_sum[2], 0, u_sum[3])) / sqrt(var)
    T_O <- max(abs(XLWP[c(1, 3, 5)]))  # 奇数γ的最大值
    T_E <- max(XLWP[c(2, 4, 6)])       # 偶数γ的最大值
    
    # 并行计算协方差第三项 - 不等协方差版本
    cov_result_P3 <- foreach(i = 1:nrow(Greek_Ast), .combine = rbind) %dopar% {
      c1 <- Greek_Ast[i, 1]; c2 <- Greek_Ast[i, 2]; c3 <- Greek_Ast[i, 3]
      d1 <- Greek_Ast[i, 4]; d2 <- Greek_Ast[i, 5]; d3 <- Greek_Ast[i, 6]
      
      # 不等协方差：使用两个协方差矩阵的组合
      mat <- mat1_cov1^c1 * mat1_cov2^d1 * mat2_cov1^c2 * mat2_cov2^d2 * cov1^c3 * cov2^d3
      nume <- sum(mat)
      deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * factorial(c1) * factorial(c2) *
                 factorial(d1) * factorial(d2) * factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
      nume / deno
    }
    
    # 协方差计算
    cov_P3 <- c(6 * sum(cov_result_P3[1:4]), 48 * sum(cov_result_P3[5:10]),
                120 * sum(cov_result_P3[11:16]), 720 * sum(cov_result_P3[17:36]),
                1440 * sum(cov_result_P3[37:45]), 17280 * sum(cov_result_P3[46:73]))
    
    cov_P1 <- c(u_sum[2], u_sum[3], u_sum[3], u_sum[4], u_sum[4], u_sum[5])
    cov_P2 <- c(0, sum(result_P0[1, ] * result_P0[2, ]), 0, 0, 
                sum(result_P0[1, ] * result_P0[3, ]), sum(result_P0[2, ] * result_P0[3, ]))
    cov <- cov_P1 - cov_P2 + cov_P3
    
    # 计算相关矩阵
    cov_e <- cov[c(2, 5, 6)] / sqrt(c(var[2] * var[4], var[2] * var[6], var[4] * var[6]))
    cov_o <- cov[c(1, 3, 4)] / sqrt(c(var[1] * var[3], var[1] * var[5], var[5] * var[3]))
    
    R_E <- matrix(c(1, cov_e[1], cov_e[2], cov_e[1], 1, cov_e[3], cov_e[2], cov_e[3], 1), 3, 3)
    R_O <- matrix(c(1, cov_o[1], cov_o[2], cov_o[1], 1, cov_o[3], cov_o[2], cov_o[3], 1), 3, 3)
    
    # 计算p值（多元正态分布）
    pval_O <- 1 - pmvnorm(lower = -rep(T_O, 3), upper = rep(T_O, 3), mean = rep(0, 3), sigma = R_O)
    pval_E <- 1 - pmvnorm(lower = rep(-Inf, 3), upper = rep(T_E, 3), mean = rep(0, 3), sigma = R_E)
    
    # 计算γ=∞的情况（Cai et al. 2014）- 不等协方差版本
    nume <- diff^2
    deno <- diag_cov1/n1 + diag_cov2/n2  # 不等协方差的方差估计
    deno[deno <= 1e-10] <- 1e-10
    
    T_CLX <- max(nume / deno)
    M <- T_CLX - 2 * log(p) + log(log(p)) + log(pi)
    pval_inf <- 1 - exp(-exp(-M / 2))
    
    # 最终p值（Bonferroni调整）
    pval.min <- min(c(pval_O, pval_E, pval_inf))
    pval <- 1 - (1 - pval.min)^3
    
    return(list(P_value = pval))
  }
  
  # 主程序逻辑
  if(eq.cov) {
    # 等协方差情况
    n1 <- nrow(sam1); n2 <- nrow(sam2); p <- ncol(sam1)
    sam.cov <- ((n1 - 1) * cov(sam1) + (n2 - 1) * cov(sam2)) / (n1 + n2 - 2)
    
    # 带宽参数处理
    if(missing(bandwidth)) {
      bandwidth <- seq(from = 0, to = p, by = floor(p / 50))
    }
    bandwidth[bandwidth < 0] <- 0
    bandwidth <- floor(bandwidth)
    
    # 协方差矩阵估计
    if(missing(cov.same)) {
      if(length(bandwidth) > 1) {
        centered_data <- rbind(
          sam1 - matrix(colMeans(sam1), byrow = TRUE, nrow = n1, ncol = p),
          sam2 - matrix(colMeans(sam2), byrow = TRUE, nrow = n2, ncol = p)
        )
        optim.bandwidth <- best.band_parallel(centered_data, bandwidth, cv.fold, norm)
      } else {
        optim.bandwidth <- bandwidth
      }
      
      if(optim.bandwidth > 0) {
        cov.same <- sam.cov
        cov.same[abs(row(cov.same) - col(cov.same)) > optim.bandwidth] <- 0
      } else {
        cov.same <- diag(diag(sam.cov))
      }
    }
    
    # 执行等协方差检验
    out <- Xu_same_parallel(sam1, sam2, cov.same, G_greek)
    out$time <- Sys.time() - start
    return(out)
    
  } else {
    # 不等协方差情况
    sam.cov1 <- cov(sam1)
    sam.cov2 <- cov(sam2)
    p <- ncol(sam1)
    
    # 带宽参数处理（两个样本分别处理）
    if(missing(bandwidth1)) {
      bandwidth1 <- seq(from = 0, to = p, by = floor(p / 50))
    }
    bandwidth1[bandwidth1 < 0] <- 0
    bandwidth1 <- floor(bandwidth1)
    
    if(missing(bandwidth2)) {
      bandwidth2 <- seq(from = 0, to = p, by = floor(p / 50))
    }
    bandwidth2[bandwidth2 < 0] <- 0
    bandwidth2 <- floor(bandwidth2)
    
    # 协方差矩阵估计（样本1）
    if(missing(cov1)) {
      if(length(bandwidth1) > 1) {
        optim.bandwidth1 <- best.band_parallel(sam1, bandwidth1, cv.fold, norm)
      } else {
        optim.bandwidth1 <- bandwidth1
      }
      
      if(optim.bandwidth1 > 0) {
        cov1 <- sam.cov1
        cov1[abs(row(cov1) - col(cov1)) > optim.bandwidth1] <- 0
      } else {
        cov1 <- diag(diag(sam.cov1))
      }
    }
    
    # 协方差矩阵估计（样本2）
    if(missing(cov2)) {
      if(length(bandwidth2) > 1) {
        optim.bandwidth2 <- best.band_parallel(sam2, bandwidth2, cv.fold, norm)
      } else {
        optim.bandwidth2 <- bandwidth2
      }
      
      if(optim.bandwidth2 > 0) {
        cov2 <- sam.cov2
        cov2[abs(row(cov2) - col(cov2)) > optim.bandwidth2] <- 0
      } else {
        cov2 <- diag(diag(sam.cov2))
      }
    }
    
    # 执行不等协方差检验
    out <- Xu_diff_parallel(sam1, sam2, cov1, cov2, G_greek)
    out$time <- Sys.time() - start
    return(out)
  }
  # 关闭并行环境
  stopCluster(cl)
}
XLWP_parallel_else(sam1,sam2,T)
XLWP_parallel_else(sam1,sam2,F)
