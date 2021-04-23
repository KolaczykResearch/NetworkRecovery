# Libraries ----
library(igraph)
library(parallel)

# Functions ----
source("/functions.R")

# Simulation ----
n <- c(50, 100, 250, 500, 1000) # |V|
x.grid <- c(.9, 1.1)            # sqrt(delta) * log(n)
reps <- 1000                    # number of replicates

grid <- expand.grid(n, x.grid, stringsAsFactors = FALSE)
colnames(grid) <- c("n", "x.grid")

res <- matrix(nrow = nrow(grid), ncol = 5)
colnames(res) <- c("n", "x", "all", "match", "proportion")

set.seed(022321)
for (g in seq_len(nrow(grid))) {
  
  # Parameters ----
  n <- grid[g, "n"]
  x <- grid[g, "x.grid"]
  p <- log(n)^2 / n
  beta <- (x / log(n))^2
  
  start <- Sys.time()
  vals <- simplify2array(
    mclapply(seq_len(reps), mc.cores = detectCores(), FUN = function (r) {
      set.seed(r)
      
      # sample parent A ~ ER(n, p)
      A <- erdos.renyi.game(n, p, type = "gnp")
      
      E <- length(E(A))
      E_c <- n*(n-1)/2 - E
      
      log.lik <- E*log(p) + E_c*log(1-p)
      
      alpha <- beta * E / E_c # edge unbiased
      
      # noisy copies of A
      A_i <- rcorrER(m = 2, A, alpha, beta)
      A_i <- lapply(A_i, graph_from_adjacency_matrix, mode = "undirected")
      
      # scramble nodes of A^(i)
      pi_star <- sample(n)
      A_i[[2]] <- graph_from_adjacency_matrix(get.adjacency(A_i[[2]])[pi_star, pi_star]
                                              , mode = "undirected")
      
      # match A^(1) to A^(2)
      pi_hat <- match(A_i[[1]], A_i[[2]], verbose = FALSE)$pi_hat
      
      c(log.lik, identical(pi_star[pi_hat], 1:n))
      
    })) # end replicates
  
  message(g, " (", round(difftime(Sys.time(), start, units = "secs")), " secs) : "
          , " n = ", n
          , " , x = ", x)
  
  rownames(vals) <- c("log.lik", "match")
  
  res[g, "n"] <- n
  res[g, "x"] <- x
  res[g, "all"] <- mean(vals["log.lik", ])
  res[g, "match"] <- mean(vals["log.lik", vals["match", ] == 1])
  res[g, "proportion"] <- mean(vals["match", ])
  
} # end grid

mean(res[, 4] >= res[, 3], na.rm = TRUE)