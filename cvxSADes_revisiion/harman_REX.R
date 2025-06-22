library(OptimalDesign)

# 1. 建立 candidate grid
N1 <- 11  # or any number
grid <- expand.grid(replicate(7, 
                              seq(-1, 1, length.out = N1),
                              simplify = FALSE))
names(grid) <- paste0("x", 1:7)

# 2. 計算 Fisher 資訊的 Fx
theta0 <- c(-0.4926,-0.6280,-0.3283,0.4378,0.5283,
            0.6120,-0.6837,-0.2061)
Fx <- Fx_glm(~ x1 + x2 + x3 + x4+ x5 +x6 + x7, 
             theta0 = theta0,
             glm.model = "bin-logit",
             lower = rep(-1, 7),
             upper = rep(+1, 7),
             n.levels = rep(N1, 7))


# 3. 用 od_REX 計算 approximate D‑optimal 設計
res_REX <- od_REX(Fx, crit="D", alg.AA="REX", eff=0.999999, track=TRUE)

# 4. 抽出最終的支撐點及其權重
design_REX <- data.frame(grid[res_REX$supp, ], weight = res_REX$w.supp) |> arrange(x1, x2, x3, x4, x5, x6, x7)


## From our paper

(design_mine <- data.frame(
  x1 = c(rep(-1, 13), rep(1, 16)),
  x2 = c(rep(-1, 5), rep(1, 8), rep(-1, 9), rep(1, 7)),
  x3 = c(rep(-1, 3), rep(1, 2), rep(-1, 4), rep(1, 4), rep(-1, 3),
         rep(1, 6), rep(-1, 2), rep(1, 5)),
  x4 = c(-1, rep(1, 2), rep(-1, 4), rep(1, 2), -1, rep(1, 3), -1, rep(1, 2), rep(-1,3), rep(1, 3), -1, 1, -1, rep(1, 4)),
  x5 = c(rep(-1, 2),  rep(1, 3), -1, 1, -1, 1, rep(-1, 2), rep(1,2),rep(-1, 2), 1, rep(-1, 2), 1, -1, rep(1,2), rep(-1,6), 1),
  x6 = c(rep(1, 2), rep(-1, 2), 1, rep(-1,2), 1, rep(-1, 2), 1, -1, 1, rep(-1, 4), 1, -1, 1, -1,  1, rep(-1, 4), rep(1, 2), -1),
  x7 = c(1, -1,  1, rep(-1, 4), 1, -1, 1, -1, rep(1, 6), rep(-1, 2), 1, rep(-1, 4), rep(1, 2), -1, rep(1, 2)),
  weight = c(0.0627, 0.0732, 0.0487, 0.0499, 0.0460, 0.0088, 0.0561,
             0.0212, 0.0226, 0.0840, 0.0306, 0.0023, 0.0730, 0.0135,
             0.0217, 0.0415, 0.0409, 0.0375, 0.0073, 0.0255, 0.0489,
             0.0100, 0.0491, 0.0404, 0.0042, 0.0058, 0.0420, 0.0033, 0.0295))) |> 
  arrange(x1, x2, x3, x4, x5, x6, x7)


calc_info_matrix <- function(design_df, cols = 1:7, intercept = TRUE, digits = 4) {
  
  X <- as.matrix(design_df[, cols])  # extract design matrix
  if (intercept) {
    X <- cbind(1, X)  # add intercept column
  }
  w <- design_df$weight/sum(design_df$weight)  # extract weights
  
  M <- matrix(0, ncol = ncol(X), nrow = ncol(X))
  for (i in seq_len(nrow(X))) {
    g <- X[i, , drop = FALSE]
    M <- M + w[i] * t(g) %*% g
  }
  
  return(round(M, digits))
}

# Assuming design_REX is your data frame
design_REX$weight <- round(design_REX$weight,4) # normalize weights

FIM_REX <- calc_info_matrix(design_REX)
FIM_mine <- calc_info_matrix(design_mine)
det(solve(FIM_REX))^(1/8)
det(solve(FIM_mine))^(1/8)



