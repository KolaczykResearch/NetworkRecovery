# libraries ----
library(igraph)
library(parallel)

# functions ----
source("/R/functions.R")

# parameters ----
n <- c(50, 100, 200)            # |V|
m <- c(5, 10)                   # number of networks
x.grid <- seq(.3, 2.3, .2)      # sqrt(delta) * log(n)
reps <- 10                      # number of replicates

grid <- expand.grid(n, m, x.grid, stringsAsFactors = FALSE)
colnames(grid) <- c("n", "m", "x.grid")

res <- matrix(nrow = nrow(grid), ncol = 9)
colnames(res) <- c("n", "m", "p", "beta", "alpha", "x"
                   , "Frobenius", "Accuracy", "Recovery")

for (g in seq_len(nrow(grid))) {
  n <- grid[g, "n"]
  m <- grid[g, "m"]
  x <- grid[g, "x.grid"]
  
  p <- log(n)^2 / n
  beta <- (x / log(n))^2
  alpha <- beta * p / (1-p) # edge unbiased
  
  start <- Sys.time()
  vals <- simplify2array(
    mclapply(seq_len(reps), mc.cores = detectCores(), FUN = function (r) {
      set.seed(r)
      
      # sample parent A ~ ER(n, p)
      A <- erdos.renyi.game(n, p, type = "gnp")
      
      E <- length(E(A))
      alpha <- beta * E / (n*(n-1)/2 - E)
      
      # m noisy copies of A
      A_i <- rcorrER(m, A, alpha, beta)
      
      # scramble edges
      pi_star <- matrix(replicate(m-1, sample(n)), nrow = m-1, byrow = TRUE)
      for (i in 2:m) {
        A_i[[i]] <- A_i[[i]][pi_star[i-1, ], pi_star[i-1, ]] 
      }
      
      A_i <- lapply(A_i, graph_from_adjacency_matrix, mode = "undirected")
      
      # match ----
      match <- match.multi(A_i, parallel = TRUE, cleanup = TRUE)
      A_match <- match$A_match
      pi_hat <- match$pi_hat
      
      # average ----
      A_avg <- Reduce("+", A_match) / m
      A_hat <- as.matrix(A_avg)
      A_hat[A_hat >= .5] <- 1
      A_hat[A_hat < .5] <- 0
      
      # evaluate ----
      Frobenius <- sum((get.adjacency(A) - A_avg)^2) / (n*(n-1))
      Accuracy <- sum(get.adjacency(A) == A_hat) / n^2
      Recovery <- sum(pi_hat == pi_star) / (n*(m-1))

      return(c(Frobenius, Accuracy, Recovery))
    }))
  
  if (is.null(dim(vals))) {
    message("Investigate ", g)
    vals <- simplify2array(vals[which(sapply(vals, is.numeric))])
  } else {
    message(g, " (", round(difftime(Sys.time(), start, units = "secs")), " secs) : "
            , " m = ", m
            , " , n = ", n
            , " , x = ", x
            , " , Frobenius = ", round(median(vals[1, ]), 2)
            , " , Accuracy = ", round(median(vals[2, ]), 2)
            , " , Recovery = ", round(median(vals[3, ]), 2))
  }
  
  res[g, "n"] <- n  
  res[g, "m"] <- m
  res[g, "p"] <- p
  res[g, "alpha"] <- alpha
  res[g, "beta"] <- beta
  res[g, "x"] <- x
  res[g, "Frobenius"] <- median(vals[1, ])
  res[g, "Accuracy"] <- median(vals[2, ])
  res[g, "Recovery"] <- median(vals[3, ])
}

# plot ----
library(reshape2)
library(ggplot2)

df <- melt(as.data.frame(res), id.vars = c("n", "m", "p", "alpha", "beta", "x"))
df[, c(1:6, 8)] <- lapply(df[, c(1:6, 8)], as.numeric)
names(df) <- c("n", "m", "p", "alpha", "beta", "x", "Measure", "Value")
df$Measure <- ordered(df$Measure, levels = c("Frobenius", "Accuracy", "Recovery"))

ggplot(df, aes(x = x, y = Value, group = factor(n), color = factor(n))) +
  geom_point() +
  geom_line() +
  xlab(expression(sqrt(beta)~log(n))) + ylab("") + labs(color = "n") +
  facet_grid(~m~Measure, scales = "free_y") +
  scale_color_viridis_d(option = "C") +
  theme_bw()
