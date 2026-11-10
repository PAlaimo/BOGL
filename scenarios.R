
# Packages ----------------------------------------------------------------
library(mvtnorm)


# Utility functions -------------------------------------------------------


make_sigma <- function(p, rho) {
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  Sigma
}


# Scenario structure -------------------------------------------------------

build_random_groups <- function(p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  G  <- p / 2
  ng <- sample(1:5, G, replace = TRUE)
  
  pg <- vector("list", G)
  for (g in 1:G) {
    pg[[g]] <- sort(sample(1:p, ng[g]))
  }
  
  covered <- sort(unique(unlist(pg)))
  missing <- setdiff(1:p, covered)
  
  if (length(missing) > 0) {
    for (j in missing) {
      g_add <- sample(1:G, 1)
      pg[[g_add]] <- sort(unique(c(pg[[g_add]], j)))
    }
  }
  
  pg
}


build_interaction_design <- function(n, rho) {
  # fixed scenario with p_base = 4
  p_base <- 4
  
  Sigma <- make_sigma(p_base, rho)
  X0    <- rmvnorm(n, sigma = Sigma)
  colnames(X0) <- c("x1", "x2", "x3", "x4")
  
  X2 <- X0^2
  colnames(X2) <- c("x1^2", "x2^2", "x3^2", "x4^2")
  
  X1_cross <- cbind(X0[, 1] * X0[, -1])
  colnames(X1_cross) <- c("x1x2", "x1x3", "x1x4")
  
  X2_cross <- cbind(X0[, 2] * X0[, -c(1:2)])
  colnames(X2_cross) <- c("x2x3", "x2x4")
  
  X3_cross <- cbind(X0[, 3] * X0[, -c(1:3)])
  colnames(X3_cross) <- c("x3x4")
  
  X <- cbind(X0, X2, X1_cross, X2_cross, X3_cross)
  
  pg <- vector("list", 14)
  
  # singleton groups for nonlinear/interaction terms
  for (g in 1:10) {
    pg[[g]] <- g + 4
  }
  
  # Hierarchical groups on principal effect
  pg[[11]] <- c(1, 5, 9, 10, 11)
  pg[[12]] <- c(2, 6, 9, 12, 13)
  pg[[13]] <- c(3, 7, 10, 12, 14)
  pg[[14]] <- c(4, 8, 11, 13, 14)
  
  list(
    X = X,
    Groups = pg,
    p = ncol(X)
  )
}

build_nested_groups <- function(p) {
  pg <- vector("list", p)
  for (g in 1:p) {
    pg[[g]] <- g:p
  }
  pg
}


# Random generating functions --------------------------------------------


make_beta <- function(p, Groups, null_groups = c(1, 2), seed = NULL) {
  
  # Setting the seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Penalized coefficients
  beta <- rnorm(p)
  # Set to 0
  null_idx <- sort(unique(unlist(Groups[null_groups])))
  beta[null_idx] <- 0

  return(beta)
}


make_response <- function(X, beta, sigma = 1, seed = NULL) {
  
  # Setting the seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # evaluating signal noise
  eps <- rnorm(nrow(X), sd = sigma)
  y   <- as.numeric(X %*% beta) + eps
  
  list(
    y = y,
    eps = eps,
    signal_sd = sigma
  )
}

# One simulation replicate ------------------------------------------------------------


generate_one_scenario <- function(
    scenario = c("random"),
    rho = 0,
    n = 200,
    p0 = 1,
    int0 = T,
    p_random = 10,
    p_bilevel = 16,
    p_nested = 10,
    seed = NULL
) {
  
  # Setting the seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  if (scenario == "random") {
    # Size of penalized coefficients
    p <- p_random
    # Group constructions
    Groups <- build_random_groups(p)
    # Correlation of covariates of penalized coefficients
    Sigma <- make_sigma(p, rho)
    # Simulated covariates of penalized coefficients
    X <- rmvnorm(n, sigma = Sigma)
    
    # beta generation
    beta <- make_beta(
      p = p,
      Groups = Groups,
      null_groups = c(1, 2)
    )
    
  }  else if (scenario == "interactions") {
    tmp <- build_interaction_design(n = n, rho = rho)
    X <- tmp$X
    Groups <- tmp$Groups
    p <- tmp$p
    beta <- make_beta(
      p = p,
      Groups = Groups,
      null_groups = c(1, 2)
    )
    
  } else if (scenario == "nested") {
    p <- p_nested
    Sigma <- make_sigma(p, rho)
    X <- rmvnorm(n, sigma = Sigma)
    Groups <- build_nested_groups(p)
    beta <- make_beta(
      p = p,
      Groups = Groups,
      null_groups = 3:10
    )
    
  }
  
  # Design matrix of unpenalized coefficients
  X0 <- as.matrix(rep(1, n))
  if(p0>1) X0 <- cbind(X0, rmvnorm(n, diag(p0-1)))
  
  # Unpenalized coefficients
  beta0 <- rnorm(p0)
  if(int0) beta0[1] <- 0
  

  # Inactive penalized coefficients
  null_idx   <- if(scenario != "nested") sort(unique(unlist(Groups[c(1, 2)]))) else sort(unique(unlist(Groups[-c(1, 2)])))

  # Active penalized coefficients
  active_idx <- c((1:p)[-null_idx])
  
  list(
    scenario = scenario,
    rho = rho,
    n = n,
    p0 = p0,
    X0 = X0,
    beta0 = beta0,
    p = p,
    X = X,
    beta = beta,
    Groups = Groups,
    null_groups = c(1, 2),
    null_idx = null_idx,
    active_idx = active_idx
  )
}

# Full simulation design  ------------------------------------------------------------

generate_simulation_design <- function(
    scenario = c("random"),
    rho = c(0),
    n_rep = 100,
    n = 200,
    p0 = 1,
    int0 = T,
    p_random = 10,
    p_bilevel = 16,
    p_nested = 10,
    noise_sd = 1,
    seed = 123
) {
  
  # Storing objects
  out <- list()

  # Generating data
  datnoy <- generate_one_scenario(
    scenario = scenario,
    rho = rho,
    n = n,
    p0 = 1,
    int0 = T,
    p_random = p_random,
    p_bilevel = p_bilevel,
    p_nested = p_nested,
    seed = seed
  )
  datnoy$noise_sd <- noise_sd
  datnoy$sc <- paste0(scenario, "_rho", rho)
  
  out[[1]] <- datnoy
  
  # Full design matrix
  Xtilde <- cbind(datnoy$X0, datnoy$X)
  # Full vector of coefficients
  betatilde <- c(datnoy$beta0, datnoy$beta)
  
  for (rep in 1:n_rep) {
     
    # Generation of y
    y_out <- make_response(X = Xtilde, beta = betatilde, sigma = noise_sd, seed = seed+rep)$y
    
    
    out[[rep+1]] <- y_out
  }

  out
}
