library(OptimalDesign)

# Step 1: Define number of levels per dimension (reduce if too slow)
N1 <- 4  # e.g., 3 levels per variable ⇒ 3^7 = 2187 candidate points

D <- 7
# Step 2: Set parameter vector theta0 (logistic regression coefficients)
theta0 <- c(-0.4926,-0.6280,-0.3283,0.4378,0.5283,
            0.6120,-0.6837,-0.2061)

# Step 3: Construct candidate design matrix for binomial-logit GLM
Fx <- Fx_glm(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, 
             theta0 = theta0,
             glm.model = "bin-logit",
             lower = rep(-1, D),
             upper = rep(+1, D),
             n.levels = rep(2, D))

# Step 4: Compute D-optimal exact design using KL exchange algorithm
res <- od_KL(Fx, N = 30, crit = "D", bin = TRUE, t.max = 30)

# Step 5: Summarize and evaluate
cat("Support indices:\n")
print(res$supp)

cat("\nDesign weights (replicates at each support):\n")
print(res$w.supp)

cat("\nDesign efficiency:\n")
print(res$eff.best)

cat("\nD-optimal criterion value (log-determinant of Fisher info):\n")
print(res$Phi.best)

# Step 6: Visual inspection of information matrix
M <- res$M.best
cat("\nInformation matrix M:\n")
print(M)

# (Optional) Check criterion again with optcrit
crit_val <- optcrit(Fx, res$w.best, crit = "D")
cat("\nValidated D-optimal criterion:\n")
print(crit_val)


Fx[res$supp, ]

## DO this manually 

# Step 1: Build the exact grid
grid_points <- expand.grid(rep(list(c(-1, 0, 1)), 7))
colnames(grid_points) <- paste0("x", 1:7)
n <- nrow(grid_points)

# Step 2: Create design matrix (including intercept)
X <- model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, data = grid_points)

# Step 3: Logistic function
theta0 <- c(-0.4926, -0.6280, -0.3283, 0.4378, 0.5283,
            0.6120, -0.6837, -0.2061)

eta <- X %*% theta0
mu <- 1 / (1 + exp(-eta))
W <- as.vector(mu * (1 - mu))  # weight: variance function

# Step 4: Multiply rows of X by sqrt(W)
X_weighted <- X * sqrt(W)

# Step 5: This is your Fx
Fx <- X_weighted
res <- od_KL(Fx, N = 21, crit = "D", bin = TRUE, t.max = 30)

my_res <- X[res$supp, ] |> as.matrix()

my_res_df <- as.data.frame(my_res)
my_res_df$n <- res$w.supp
rownames(my_res_df) <- NULL
my_res_df
