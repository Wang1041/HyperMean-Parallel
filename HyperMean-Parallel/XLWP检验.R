XLWP_test <- function(sam1, sam2, eq.cov = TRUE, cov.same, cov1, cov2, 
                   G_greek = c(1:6, Inf), bandwidth, bandwidth1, bandwidth2, 
                   cv.fold = 5, norm = "F") {
  # Xu et al.(2016) 高维两样本均值检验
  # 结合SPU检验统计量与极值统计量
  # 输入参数：
  #   sam1, sam2: 样本矩阵（行=观测，列=变量）
  #   eq.cov: 是否假设协方差矩阵相等
  #   cov.same/cov1/cov2: 预指定的协方差矩阵
  #   G_greek: γ值集合，默认1-6和无穷大
  #   bandwidth/bandwidth1/bandwidth2: 带宽参数
  #   cv.fold: 交叉验证折数
  #   norm: 矩阵范数类型
  
  start <- Sys.time()
  
  # 辅助函数：最优带宽选择
  best.band <- function(sam, bandwidth, cv.fold, norm) {
    p <- ncol(sam)
    n <- nrow(sam)
    fold.size <- round(n / cv.fold)
    sam.idx <- sample(1:n, size = n, replace = FALSE)
    n.bandwidth <- length(bandwidth)
    diff.norm <- matrix(0, cv.fold, n.bandwidth)
    matrix.id <- abs(row(diag(p)) - col(diag(p)))  # 距离矩阵
    
    for(i in 1:cv.fold) {
      # 划分训练集和测试集
      if(i == cv.fold) {
        temp.idx <- sam.idx[((i - 1) * fold.size + 1):n]
      } else {
        temp.idx <- sam.idx[((i - 1) * fold.size + 1):(i * fold.size)]
      }
      
      sam.train <- sam[-temp.idx, ]
      sam.test <- sam[temp.idx, ]
      sam.train.cov <- cov(sam.train)
      sam.test.cov <- cov(sam.test)
      
      # 评估不同带宽
      for(j in 1:n.bandwidth) {
        sam.train.cov.band <- sam.train.cov
        sam.train.cov.band[matrix.id > bandwidth[j]] <- 0  # 带宽截断
        diff.norm[i, j] <- norm(sam.train.cov.band - sam.test.cov, type = norm)
      }
    }
    
    # 选择最优带宽
    diff.norm <- colMeans(diff.norm)
    best <- which.min(diff.norm)
    return(bandwidth[best[1]])
  }
  
  if(eq.cov) {
    # ========== 等协方差情况 ==========
    n1 <- nrow(sam1)
    n2 <- nrow(sam2)
    p <- ncol(sam1)
    
    # 计算合并协方差矩阵
    sam.cov <- ((n1 - 1) * cov(sam1) + (n2 - 1) * cov(sam2)) / (n1 + n2 - 2)
    
    # 带宽参数处理
    if(missing(bandwidth)) {
      bandwidth <- seq(from = 0, to = p, by = floor(p / 50))
    }
    bandwidth[bandwidth < 0] <- 0
    bandwidth <- floor(bandwidth)
    
    # 协方差矩阵估计（带带宽选择）
    if(missing(cov.same)) {
      output.opt.bw <- TRUE
      if(length(bandwidth) > 1) {
        # 中心化数据用于带宽选择
        centered_data <- rbind(
          sam1 - matrix(colMeans(sam1), byrow = TRUE, nrow = n1, ncol = p),
          sam2 - matrix(colMeans(sam2), byrow = TRUE, nrow = n2, ncol = p)
        )
        optim.bandwidth <- best.band(centered_data, bandwidth, cv.fold, norm)
      } else {
        optim.bandwidth <- bandwidth
      }
      
      # 根据最优带宽处理协方差矩阵
      if(optim.bandwidth > 0) {
        cov.same <- sam.cov
        cov.same[abs(row(cov.same) - col(cov.same)) > optim.bandwidth] <- 0
      } else {
        cov.same <- diag(diag(sam.cov))  # 只保留对角线
      }
    } else {
      output.opt.bw <- FALSE
    }
    
    # 内部函数：等协方差情况下的检验
    xu2016_samecov <- function(sam1, sam2, cov.same, G_greek = c(1:6, Inf)) {
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
      
      # 计算SPU统计量的均值
      SPU_u <- function(n1, n2, gamma, cov.same) {
        if(gamma %% 2 == 0) {
          gamma.half <- gamma / 2
          u_i <- 0
          for(d in 0:gamma.half) {
            u_i <- u_i + 1 / (factorial(d) * factorial(gamma.half - d) * n1^d * n2^(gamma.half - d))
          }
          u <- (diag(cov.same))^gamma.half * factorial(gamma) * u_i / 2^gamma.half 
        } else {
          u <- 0  # 奇数γ的均值为0
        }
        return(u)
      }
      
      # 计算方差第三项
      var_other <- function(n1, n2, gamma, cov.same) {
        diag_cov <- diag(cov.same)
        p <- length(diag_cov)
        mat1 <- matrix(rep(diag_cov, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diag_cov, p), p, p, byrow = TRUE)
        diag(cov.same) <- 0  # 剔除对角线
        
        P3 <- 0
        for(c3 in 0:gamma) {
          for(d3 in 0:(gamma - c3)) {
            if(c3 + d3 > 0) {
              for(c1 in 0:floor((gamma - c3 - d3) / 2)) {
                for(c2 in 0:floor((gamma - c3 - d3) / 2)) {
                  if((gamma - 2*c1 - c3 - d3) %% 2 == 0 & (gamma - 2*c2 - c3 - d3) %% 2 == 0) {
                    d1 <- (gamma - 2*c1 - c3 - d3) / 2
                    d2 <- (gamma - 2*c2 - c3 - d3) / 2
                    mat <- mat1^(c1 + d1) * mat2^(c2 + d2) * cov.same^(c3 + d3)
                    nume <- (factorial(gamma))^2 * sum(mat)
                    deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * 
                               factorial(c1) * factorial(c2) * factorial(d1) * factorial(d2) *
                               factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
                    P3 <- P3 + nume / deno
                  }
                }
              }
            }
          }
        }
        return(P3)
      }
      
      # 计算SPU统计量的方差
      SPU_var <- function(n1, n2, gamma, cov.same) {
        p <- nrow(cov.same)
        if(gamma == 1) {
          var <- (1/n1 + 1/n2) * sum(cov.same)
        } else {
          P1 <- sum(SPU_u(n1, n2, 2*gamma, cov.same))
          P2 <- sum((SPU_u(n1, n2, gamma, cov.same))^2)
          P3 <- var_other(n1, n2, gamma, cov.same)
          var <- P1 - P2 + P3
        }
        return(var)
      }
      
      # 计算奇偶统计量T_O和T_E
      T_O <- T_E <- 0
      Var.ss <- numeric(length(G_greek0))
      
      for(i in G_greek0) {
        Var_i <- SPU_var(n1, n2, i, cov.same)
        Var.ss[i] <- Var_i
        XLWP <- (sum(diff^i) - sum(SPU_u(n1, n2, i, cov.same))) / sqrt(Var_i)
        
        if(i %% 2 == 1) {
          if(abs(XLWP) >= T_O) T_O <- abs(XLWP)
        } else {
          if(XLWP >= T_E) T_E <- XLWP
        }
      }
      
      # 计算不同γ值间的相关系数
      cov_other <- function(n1, n2, s, t, cov.same) {
        diag_cov <- diag(cov.same)
        p <- length(diag_cov)
        mat1 <- matrix(rep(diag_cov, p), p, p, byrow = FALSE)
        mat2 <- matrix(rep(diag_cov, p), p, p, byrow = TRUE)
        diag(cov.same) <- 0
        
        P3 <- 0
        for(c3 in 0:min(c(t, s))) {
          for(d3 in 0:min(c(t - c3, s - c3))) {
            if(c3 + d3 > 0) {
              for(c1 in 0:floor((t - c3 - d3) / 2)) {
                for(c2 in 0:floor((s - c3 - d3) / 2)) {
                  if((t - 2*c1 - c3 - d3) %% 2 == 0 & (s - 2*c2 - c3 - d3) %% 2 == 0) {
                    d1 <- (t - 2*c1 - c3 - d3) / 2
                    d2 <- (s - 2*c2 - c3 - d3) / 2
                    mat <- mat1^(c1 + d1) * mat2^(c2 + d2) * cov.same^(c3 + d3)
                    nume <- factorial(s) * factorial(t) * sum(mat)
                    deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * 
                               factorial(c1) * factorial(c2) * factorial(d1) * factorial(d2) *
                               factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
                    P3 <- P3 + nume / deno
                  }
                }
              }
            }
          }
        }
        return(P3)
      }
      
      SPU_cov <- function(n1, n2, s, t, cov.same) {
        p <- nrow(cov.same)
        P1 <- sum(SPU_u(n1, n2, s + t, cov.same))
        P2 <- sum(SPU_u(n1, n2, s, cov.same) * SPU_u(n1, n2, t, cov.same))
        P3 <- cov_other(n1, n2, s, t, cov.same)
        return(P1 - P2 + P3)
      }
      
      # 计算相关矩阵R_O和R_E
      R_O <- matrix(0, n.odd.ga, n.odd.ga)
      R_E <- matrix(0, n.even.ga, n.even.ga)
      
      for(s in odd.ga) {
        for(t in odd.ga) {
          if(s < t) {
            io <- which(odd.ga == s)
            jo <- which(odd.ga == t)
            i <- which(G_greek0 == s)
            j <- which(G_greek0 == t)
            R_O[io, jo] <- SPU_cov(n1, n2, s, t, cov.same) / sqrt(Var.ss[i] * Var.ss[j])
          }
        }
      }
      
      for(s in even.ga) {
        for(t in even.ga) {
          if(s < t) {
            ie <- which(even.ga == s)
            je <- which(even.ga == t)
            i <- which(G_greek0 == s)
            j <- which(G_greek0 == t)
            R_E[ie, je] <- SPU_cov(n1, n2, s, t, cov.same) / sqrt(Var.ss[i] * Var.ss[j])
          }
        }
      }
      
      # 对称化相关矩阵
      R_O <- R_O + t(R_O)
      diag(R_O) <- 1
      R_E <- R_E + t(R_E)
      diag(R_E) <- 1
      
      # 计算奇偶集合的p值（多元正态分布）
      pval_O <- 1 - pmvnorm(lower = -rep(T_O, n.odd.ga), upper = rep(T_O, n.odd.ga),
                            mean = rep(0, n.odd.ga), sigma = R_O)
      
      pval_E <- 1 - pmvnorm(lower = rep(-Inf, n.even.ga), upper = rep(T_E, n.even.ga),
                            mean = rep(0, n.even.ga), sigma = R_E)
      
      # 计算γ=∞的情况（Cai et al. 2014方法）
      nume <- diff^2
      if(n1 == n2) {
        deno <- diag(cov.same) / n1 * 2
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
    
    # 执行等协方差检验
    out <- xu2016_samecov(sam1, sam2, cov.same, G_greek)
    out$time <- Sys.time() - start
    return(out)
    
  } else {
    # ========== 不等协方差情况 ==========
    # （代码结构与等协方差情况类似，主要区别在于协方差矩阵的处理）
    # 分别处理两个样本的协方差矩阵
    
    sam.cov1 <- cov(sam1)
    sam.cov2 <- cov(sam2)
    p <- ncol(sam1)
    
    # 带宽参数处理（样本1）
    if(missing(bandwidth1)) {
      bandwidth1 <- seq(from = 0, to = p, by = floor(p / 50))
    }
    bandwidth1[bandwidth1 < 0] <- 0
    bandwidth1 <- floor(bandwidth1)
    
    # 带宽参数处理（样本2）
    if(missing(bandwidth2)) {
      bandwidth2 <- seq(from = 0, to = p, by = floor(p / 50))
    }
    bandwidth2[bandwidth2 < 0] <- 0
    bandwidth2 <- floor(bandwidth2)
    
    # 协方差矩阵估计（样本1）
    if(missing(cov1)) {
      output.opt.bw1 <- TRUE
      if(length(bandwidth1) > 1) {
        optim.bandwidth1 <- best.band(sam1, bandwidth1, cv.fold, norm)
      } else {
        optim.bandwidth1 <- bandwidth1
      }
      
      if(optim.bandwidth1 > 0) {
        cov1 <- sam.cov1
        cov1[abs(row(cov1) - col(cov1)) > optim.bandwidth1] <- 0
      } else {
        cov1 <- diag(diag(sam.cov1))
      }
    } else {
      output.opt.bw1 <- FALSE
    }
    
    # 协方差矩阵估计（样本2）
    if(missing(cov2)) {
      output.opt.bw2 <- TRUE
      if(length(bandwidth2) > 1) {
        optim.bandwidth2 <- best.band(sam2, bandwidth2, cv.fold, norm)
      } else {
        optim.bandwidth2 <- bandwidth2
      }
      
      if(optim.bandwidth2 > 0) {
        cov2 <- sam.cov2
        cov2[abs(row(cov2) - col(cov2)) > optim.bandwidth2] <- 0
      } else {
        cov2 <- diag(diag(sam.cov2))
      }
    } else {
      output.opt.bw2 <- FALSE
    }
    
    # 内部函数：不等协方差情况下的检验
    xu2016_diffcov <- function(sam1, sam2, cov1, cov2, G_greek = c(1:6, Inf)) {
      n1 <- nrow(sam1); n2 <- nrow(sam2); p <- ncol(sam1)
      x1_bar <- colMeans(sam1)
      x2_bar <- colMeans(sam2)
      diff <- x1_bar - x2_bar
      
      # 确保包含无穷大γ值
      if(!is.element(Inf, G_greek)) {
        G_greek <- c(G_greek, Inf)
      }
      G_greek0 <- G_greek[G_greek != Inf]
      
      # 划分奇偶集合
      odd.ga <- G_greek0[G_greek0 %% 2 == 1]
      even.ga <- G_greek0[G_greek0 %% 2 == 0]
      n.odd.ga <- length(odd.ga)
      n.even.ga <- length(even.ga)
      
      # 计算统计量均值 - 不等协方差版本
      SPU_u_diff <- function(n1, n2, gamma, cov1, cov2) {
        if(gamma %% 2 == 0) {
          gamma.half <- gamma / 2
          u_i <- 0
          for(d in 0:gamma.half) {
            u_i <- u_i + 1 / (factorial(d) * factorial(gamma.half - d) * n1^d * n2^(gamma.half - d))
          }
          # 不等协方差：分别使用两个协方差矩阵的对角线元素
          u <- (diag(cov1)^d) * (diag(cov2)^(gamma.half - d)) * factorial(gamma) * u_i / 2^gamma.half 
        } else {
          u <- 0  # 奇数γ的均值为0
        }
        return(u)
      }
      
      # 方差的第三块计算 - 不等协方差版本
      var_diff_other <- function(n1, n2, gamma, cov1, cov2) {
        t <- gamma
        s <- gamma
        
        # 提取对角元素并剔除对角元素
        diag_cov1 <- diag(cov1)
        diag_cov2 <- diag(cov2)
        p <- length(diag_cov1)
        
        # 创建四个矩阵分别对应两个协方差矩阵的行列组合
        mat1_cov1 <- matrix(rep(diag_cov1, p), p, p, byrow = FALSE)
        mat2_cov1 <- matrix(rep(diag_cov1, p), p, p, byrow = TRUE)
        mat1_cov2 <- matrix(rep(diag_cov2, p), p, p, byrow = FALSE)
        mat2_cov2 <- matrix(rep(diag_cov2, p), p, p, byrow = TRUE)
        
        diag(cov1) <- diag(cov2) <- 0  # 剔除对角线
        
        P3 <- 0
        for(c3 in 0:min(c(t, s))) {
          for(d3 in 0:min(c(t - c3, s - c3))) {
            if(c3 + d3 > 0) {
              for(c1 in 0:floor((t - c3 - d3) / 2)) {
                for(c2 in 0:floor((s - c3 - d3) / 2)) {
                  if((t - 2*c1 - c3 - d3) %% 2 == 0 & (s - 2*c2 - c3 - d3) %% 2 == 0) {
                    d1 <- (t - 2*c1 - c3 - d3) / 2
                    d2 <- (s - 2*c2 - c3 - d3) / 2
                    
                    # 不等协方差：使用两个协方差矩阵的组合
                    mat <- mat1_cov1^c1 * mat1_cov2^d1 * mat2_cov1^c2 * mat2_cov2^d2 * cov1^c3 * cov2^d3
                    nume <- (factorial(gamma))^2 * sum(mat)
                    deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * factorial(c1) * factorial(c2) *
                             factorial(d1) * factorial(d2) * factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
                    P3 <- P3 + nume / deno
                  }
                }
              }
            }
          }
        }
        return(P3)
      }
      
      # 统计量的方差计算 - 不等协方差版本
      SPU_var_diff <- function(n1, n2, gamma, cov1, cov2) {
        p <- nrow(cov1)
        if(gamma == 1) {
          var <- sum(cov1)/n1 + sum(cov2)/n2  # 不等协方差的方差计算
        } else {
          P1 <- sum(SPU_u_diff(n1, n2, 2*gamma, cov1, cov2))
          P2 <- sum((SPU_u_diff(n1, n2, gamma, cov1, cov2))^2)
          P3 <- var_diff_other(n1, n2, gamma, cov1, cov2)
          var <- P1 - P2 + P3
        }
        return(var)
      }
      
      # 计算To和Te的值
      T_O <- T_E <- 0
      Var.ss <- numeric(length(G_greek0))
      
      for(i in G_greek0) {
        Var_i <- SPU_var_diff(n1, n2, i, cov1, cov2)
        Var.ss[i] <- Var_i
        XLWP <- (sum(diff^i) - sum(SPU_u_diff(n1, n2, i, cov1, cov2))) / sqrt(Var_i)
        
        if(i %% 2 == 1) {
          if(abs(XLWP) >= T_O) T_O <- abs(XLWP)
        } else {
          if(XLWP >= T_E) T_E <- XLWP
        }
      }
      
      # 计算不同gamma值之间的相关系数 - 不等协方差版本
      cov_diff_other <- function(n1, n2, s, t, cov1, cov2) {
        # 提取对角元素并剔除对角元素
        diag_cov1 <- diag(cov1)
        diag_cov2 <- diag(cov2)
        p <- length(diag_cov1)
        
        # 创建四个矩阵分别对应两个协方差矩阵的行列组合
        mat1_cov1 <- matrix(rep(diag_cov1, p), p, p, byrow = FALSE)
        mat2_cov1 <- matrix(rep(diag_cov1, p), p, p, byrow = TRUE)
        mat1_cov2 <- matrix(rep(diag_cov2, p), p, p, byrow = FALSE)
        mat2_cov2 <- matrix(rep(diag_cov2, p), p, p, byrow = TRUE)
        
        diag(cov1) <- diag(cov2) <- 0  # 剔除对角线
        
        P3 <- 0
        for(c3 in 0:min(c(t, s))) {
          for(d3 in 0:min(c(t - c3, s - c3))) {
            if(c3 + d3 > 0) {
              for(c1 in 0:floor((t - c3 - d3) / 2)) {
                for(c2 in 0:floor((s - c3 - d3) / 2)) {
                  if((t - 2*c1 - c3 - d3) %% 2 == 0 & (s - 2*c2 - c3 - d3) %% 2 == 0) {
                    d1 <- (t - 2*c1 - c3 - d3) / 2
                    d2 <- (s - 2*c2 - c3 - d3) / 2
                    
                    # 不等协方差：使用两个协方差矩阵的组合
                    mat <- mat1_cov1^c1 * mat1_cov2^d1 * mat2_cov1^c2 * mat2_cov2^d2 * cov1^c3 * cov2^d3
                    nume <- factorial(s) * factorial(t) * sum(mat)
                    deno <- (n1^(c1 + c2 + c3) * n2^(d1 + d2 + d3) * factorial(c1) * factorial(c2) *
                             factorial(d1) * factorial(d2) * factorial(c3) * factorial(d3) * 2^(c1 + c2 + d1 + d2))
                    P3 <- P3 + nume / deno
                  }
                }
              }
            }
          }
        }
        return(P3)
      }
      
      # 统计量的协方差计算 - 不等协方差版本
      SPU_cov_diff <- function(n1, n2, s, t, cov1, cov2) {
        p <- nrow(cov1)
        P1 <- sum(SPU_u_diff(n1, n2, s + t, cov1, cov2))
        P2 <- sum(SPU_u_diff(n1, n2, s, cov1, cov2) * SPU_u_diff(n1, n2, t, cov1, cov2))
        P3 <- cov_diff_other(n1, n2, s, t, cov1, cov2)
        cov <- P1 - P2 + P3
        return(cov)
      }
      
      # 计算RO、RE
      R_E <- matrix(0, n.even.ga, n.even.ga)
      R_O <- matrix(0, n.odd.ga, n.odd.ga)
      
      for(s in odd.ga) {
        for(t in odd.ga) {
          if(s < t) {
            io <- which(odd.ga == s)
            jo <- which(odd.ga == t)
            i <- which(G_greek0 == s)
            j <- which(G_greek0 == t)
            R_O[io, jo] <- SPU_cov_diff(n1, n2, s, t, cov1, cov2) / sqrt(Var.ss[i] * Var.ss[j])
          }
        }
      }
      
      for(s in even.ga) {
        for(t in even.ga) {
          if(s < t) {
            ie <- which(even.ga == s)
            je <- which(even.ga == t)
            i <- which(G_greek0 == s)
            j <- which(G_greek0 == t)
            R_E[ie, je] <- SPU_cov_diff(n1, n2, s, t, cov1, cov2) / sqrt(Var.ss[i] * Var.ss[j])
          }
        }
      }
      
      # 对称化相关矩阵
      R_O <- R_O + t(R_O)
      diag(R_O) <- 1
      R_E <- R_E + t(R_E)
      diag(R_E) <- 1
      
      # 计算奇偶集合的P值
      pval_O <- 1 - pmvnorm(lower = -rep(T_O, n.odd.ga), upper = rep(T_O, n.odd.ga),
                            mean = rep(0, n.odd.ga), sigma = R_O)
      
      pval_E <- 1 - pmvnorm(lower = rep(-Inf, n.even.ga), upper = rep(T_E, n.even.ga),
                            mean = rep(0, n.even.ga), sigma = R_E)
      
      # 计算gamma为∞的情况（Cai et al. 2014）- 不等协方差版本
      nume <- diff^2
      deno <- diag(cov1)/n1 + diag(cov2)/n2  # 不等协方差的方差估计
      deno[deno <= 1e-10] <- 1e-10
      
      T_CLX <- max(nume / deno)
      M <- T_CLX - 2 * log(p) + log(log(p)) + log(pi)
      pval_inf <- 1 - exp(-exp(-M / 2))
      
      # 计算最终的p值
      pval.min <- min(c(pval_O, pval_E, pval_inf))
      pval <- 1 - (1 - pval.min)^3
      
      return(list(P_value = pval))
    }
}


