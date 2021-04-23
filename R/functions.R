rcorrER <- function(m = 1, G, alpha, beta) {
  # random samples from correlated ER graph with noise
  #
  # m     : number of samples
  # G     : true underlying graph
  # alpha : Type I rate
  # beta  : Type II rate
  
  require(Matrix)
  
  if (is.igraph(G)) G <- get.adjacency(G)
  
  n <- nrow(G)
  true.edges <- which(G == 1)
  
  G_i <- replicate(m, expr = { # draw m samples
    
    # Add noisy edges wp alpha
    g <- rsparsematrix(n, n, nnz = rbinom(1, n^2, alpha)
                       , rand.x = function(n) return(1))
    
    # Add true edges from G but thin wp beta
    g[true.edges] <- sample(0:1, length(true.edges)
                            , replace = TRUE
                            , prob = c(beta, 1 - beta))
    
    # Make simple and undirected
    diag(g) <- 0
    g <- forceSymmetric(g, uplo = "L")
    
    return(g)
  })
  
  if (m == 1) G_i <- G_i[[1]]
  
  return(G_i)
}

degree.profile <- function(A) {
  # Compute degree profiles
  #
  # A : graph A
  
  require(igraph)
  
  n <- length(V(A))
  
  # compute degree for i = [n]
  # eq (17)
  a_i <- degree(A)
  
  # compute closed neighborhood, N_A[i], for i = [n]
  N_A <- neighborhood(A)
  
  return(sapply(1:n, function(i) a_i[setdiff(N_A[[i]], i)]))
  
}

match <- function(A, B, nc = 100, verbose = FALSE) {
  # Implementation of Algorithm 1 (p13)
  #
  # A, B    : graphs to be matched
  # nc      : number of clean up steps
  
  require(parallel)
  require(transport)
  
  n <- length(V(A))
  
  u_i <- degree.profile(A)
  v_i <- degree.profile(B)
  
  # matrix of 1-Wasserstein distance between degree profiles
  # SEE: remarks on middle of p52
  Z <- matrix(unlist(lapply(u_i
                            , function(i) lapply(v_i
                                                 , function(k)
                                                   wasserstein1d(i, k)
                            )), use.names = FALSE)
              , ncol = n, byrow = TRUE)
  
  # top n indices
  S <- arrayInd(which(Z %in% sort(Z)[1:n]), dim(Z), useNames = FALSE)
  
  # check if perfect match found
  if (identical(sort(S[, 1]), 1:n) && identical(sort(S[, 2]), 1:n)) {
    
    return(list(pi_hat = S[order(S[, 2]), 1], Z = Z[S]))
    
  } else {
    if (verbose) print("error -- solving assignment problem given Z")
    
    # return approximation to assignment problem given Z
    require(clue)
    
    lp <- as.numeric(solve_LSAP(Z))
    
    require(rlapjv) # faster implementation than clue
    
    # lp <- lapjv(Z)
    
    # clean up procedure in Algorithm 5 (p75)
    A <- get.adjacency(A)
    B <- get.adjacency(B)
    for (t in seq_len(nc)) {
      lp_t <- lapmod(A %*% diag(nrow = n)[lp, ] %*% B, maximize = TRUE)
      
      if (identical(lp, lp_t)) {
        break
      } else {
        lp <- lp_t
      }
    }
    
    lp.Z <- diag(nrow = n)[lp, ] * Z
    
    return(list(pi_hat = lp
                , Z = colSums(lp.Z)))
    
  }
  
}

match.multi <- function(A, parallel = TRUE, cleanup = TRUE, nc = 100) {
  # Implementation of Algorithm 1 for m > 2
  #
  # A    : list of graphs to be matched
  
  m <- length(A)
  
  # sequential pairwise match via DP+
  if (parallel) {
    require(parallel)
    
    pi <- mclapply(mc.cores = detectCores(), seq_len(m-1), function(i) {
      match(A[[i]], A[[i+1]])
    })
  } else {
    pi <- lapply(seq_len(m-1), function(i) match(A[[i]], A[[i+1]]))
  }
  
  # compose
  pi_hat <- t(sapply(pi, "[[", "pi_hat"))
  for (i in 2:(m-1)) {
    pi_hat[i, ] <- pi_hat[i, ][pi_hat[i-1, ]] # composition of permutations
  }
  
  # orient to A_1
  A <- A_match <- lapply(A, get.adjacency)
  for (i in 2:m) {
    A_match[[i]] <- A_match[[i]][pi_hat[i-1, ], pi_hat[i-1, ]]
  }
  
  # cleanup
  if (cleanup) {
    require(rlapjv)
    
    ind.1 <- 1
    for (t in seq_len(nc)) {
      ind.2 <- sample(seq_len(m)[-ind.1], 1)
      ind <- sort(c(ind.1, ind.2))
      for (tt in seq_len(nc)) {
        lp_t <- lapmod(A_match[[ind[1]]] %*% diag(nrow = n)[pi_hat[ind[2]-1, ], ] %*% A[[ind[2]]]
                       , maximize = TRUE)
        
        if (identical(pi_hat[ind[2]-1, ], lp_t)) {
          break
        } else {
          pi_hat[ind[2]-1, ] <- lp_t
        }
      }
      A_match[[ind[2]]] <- A[[ind[2]]][pi_hat[ind[2]-1, ], pi_hat[ind[2]-1, ]]
      ind.1 <- ind.2
    }
    
  }
  
  return(list(A_match = A_match, pi_hat = t(apply(pi_hat, 1, order))))
}

compose <- function(pi) {
  # Compose permutations
  #
  # pi   : list of permutations
  
  Reduce(`%*%`, pi)
}

acc <- function(pi_hat, pi_star) {
  # compute accuracy as fraction of correctly matched pairs
  #
  # pi_hat   : estimated permutation
  # pi_star  : true permutation
  
  
  mean(pi_hat == pi_star)
}