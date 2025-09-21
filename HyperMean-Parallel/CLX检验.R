CLX_test <- function(sam1, sam2) {
  # Cai & Liu & Xia (2014)的高维两样本均值检验
  # 输入：sam1和sam2为两个样本数据矩阵（行是观测，列是变量）
  # 输出：包含p值和计算时间的列表
  
  start <- Sys.time()  # 记录开始时间
  
  # 1. 计算样本均值
  x1_bar <- colMeans(sam1)  # 样本1的均值向量
  x2_bar <- colMeans(sam2)  # 样本2的均值向量
  x_cha <- x1_bar - x2_bar  # 均值差异向量
  
  # 2. 计算检验统计量的分子和分母
  nume <- x_cha^2  # 分子：均值差的平方
  n1 <- nrow(sam1); n2 <- nrow(sam2)  # 样本大小
  p <- ncol(sam1)  # 变量维度
  
  # 分母：方差估计
  deno <- (colSums(sam1^2)/n1 - x1_bar^2)/(n1-1) +  # 样本1方差
    (colSums(sam2^2)/n2 - x2_bar^2)/(n2-1)     # 样本2方差
  deno[deno <= 1e-10] <- 1e-10  # 防止数值下溢
  
  # 3. 计算极值统计量
  T_CLX <- max(nume/deno)  # 最大标准化差异
  M <- T_CLX - 2*log(p) + log(log(p)) + log(pi)  # 极值分布调整
  
  # 4. 计算p值（基于极值分布）
  alpha <- 1 - exp(-exp(-M/2))
  
  # 5. 返回结果
  out <- list()
  out$p值 <- alpha  # 检验p值
  out$time <- Sys.time() - start  # 计算耗时
  return(out)
}
CLX_test(sam1,sam2)

