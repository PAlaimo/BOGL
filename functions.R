
# Packages ----------------------------------------------------------------

library(mvtnorm)
library(statmod)



# Utility functions -------------------------------------------------------


log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

rinvgamma <- function(n, shape, scale = 1) {
  1 / rgamma(n, shape = shape, rate = scale)
}



# Prior tau2 --------------------------------------------------------------


pitau2 <- function(tau2, lambda2, p, B, m_g) {
  
  # g_j(tau2): diagonal element of D(tau2), i.e. prior precision
  gj <- sapply(1:p, function(j) sum(m_g[B[[j]]] * tau2[B[[j]]]))
  
  if (any(gj <= 0) || any(tau2 <= 0)) {
    return(0)
  }
  
  out <- prod(gj^(-1/2)) *
    prod(tau2^(-3/2)) *
    exp(-0.5 * lambda2 * sum(1 / tau2))
  
  return(out)
}



# MCMC for Bayesian Overlap Group Lasso S&S -----------------------------------------------------------------

sampler_BOGL_SpikeSlab <- function(
    y,
    X0 = NULL,
    X,
    Groups,
    weighted = T,
    d0 = 1e-6,   # prior precision for non-penalized coefficients
    alpha_sigma2 = 1,
    eta_sigma2 = 1,
    pi_g = 0.5,
    beta0_init = rep(0, ncol(X0)),
    beta_init = rep(0, ncol(X)),
    sigma2_init = NULL,
    tau2_init = 1,
    gamma_init = rep(1, length(Groups)),
    lambda2_update = c("Gibbs", "fixed"),
    lambda2_fixed = NULL,   # richiesto se lambda2_update = "fixed"
    r = 1,            # shape prior on lambda2
    delta = 1,        # rate prior on lambda2
    lambda2_init = 1,
    nMCMC = 5000,
    burn = 0.2,
    thin = 1,
    verbose = TRUE){
  
  
  #### INITIALIZATION ####
  # Set lambda2 update
  lambda2_update <- match.arg(lambda2_update)
  if (lambda2_update == "fixed" && is.null(lambda2_fixed)) {
    stop("Con lambda2_update = 'fixed' devi fornire lambda2_fixed.")
  }
  
  # Store and convert data inputs
  # Outcome
  y <- as.numeric(y)
  n <- length(y)
  
  #  Storing and checking data
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(X0)) {
    X0 <- matrix(1, nrow = n, ncol = 1)
  } else {
    X0 <- as.matrix(X0)
  }
  p0 <- ncol(X0)
  X <- as.matrix(X)
  p <- ncol(X)
  d0vec <- rep(d0, p0)

  #  Grouping structure
  G  <- length(Groups)
  B <- vector("list", p)
  for (j in 1:p) {
    B[[j]] <- which(sapply(Groups, function(gr) j %in% gr))
    if (length(B[[j]]) == 0) {
      stop(sprintf("Coefficient %d does not belong to any group.", j))
    }
  }
  # Weighting penalty by group size
  m_g <- unlist(lapply(Groups, length))
  w   <- if (weighted) m_g else rep(1, G)
  
  # Log-marginal likelihod of y for gamma update
  log_marginal_y <- function(A, sigma2_val, omj_val){
    
    # Restricting to active set
    XtildeA  <- as.matrix(Xtilde[, c(rep(T, p0), A)])
    dtildeA  <- c(d0vec, omj_val[A])
    
    # Variance matrix
    if(nrow(Xtilde) >= ncol(Xtilde)){
      VA <- diag(dtildeA, nrow=length(dtildeA)) + crossprod(XtildeA)
      RA <- chol(VA)
      ldetA <- 2*sum(log(diag(RA)))
      quadA <- sum(y^2) - sum(forwardsolve(t(RA), crossprod(XtildeA, y))^2)
    } else {
      VA <- diag(1, n) + crossprod(t(XtildeA)/sqrt(dtildeA))
      RA <- chol(VA)
      ldetA <- 2*sum(log(diag(RA)))
      quadA <- sum(forwardsolve(t(RA), y)^2)
    }
    
    # Return log of prop to marginal of y
    log_marginal <- - 0.5 * (ldetA +  quadA / sigma2_val)
    
    return(as.numeric(log_marginal))
  }
  
  
  #  Chain storing objects
  gamma_post   <- matrix(NA, nrow = nMCMC, ncol = G)
  xi_post      <- matrix(NA, nrow = nMCMC, ncol = p)
  beta0_post   <- matrix(NA, nrow = nMCMC, ncol = p0)
  b_post       <- matrix(NA, nrow = nMCMC, ncol = p)
  beta_post    <- matrix(NA, nrow = nMCMC, ncol = p)
  tau2_post    <- matrix(NA, nrow = nMCMC, ncol = G)
  sigma2_post  <- rep(NA, nMCMC)
  lambda2_post <- rep(NA, nMCMC)
  active_size  <- rep(NA, nMCMC)
  loglik_post <- matrix(NA, nrow = nMCMC, ncol = n)   # log-lik puntuale
  
  # Prior probabilities per group: all equal if one provided
  if(length(pi_g)==1) pi_g <- rep(pi_g, G)
  
  #  Starting values
  gamma_post[1, ] <- gamma_init
  xi_post[1, ] <- sapply(1:p, function(j) as.integer(all(gamma_post[1, ][B[[j]]] == 1)))
  beta0_post[1, ] <- beta0_init
  b_post[1, ]  <- beta_init
  beta_post[1, ] <- b_post[1, ]*xi_post[1, ]
  sigma2_post[1] <- if (is.null(sigma2_init)) stats::var(y) else sigma2_init
  tau2_post[1, ] <- rep(tau2_init, G)
  lambda2_post[1] <- if (lambda2_update == "fixed") lambda2_fixed else lambda2_init

  # Active size
  active_size[1] <- sum(xi_post[1, ] == 1)
  
  # Evaluating loglikelihood at starting value
  mu_init <- as.numeric(X0 %*% beta0_post[1, ] + X %*% beta_post[1, ])
  loglik_post[1, ] <- dnorm(y, mean = mu_init, sd = sqrt(sigma2_post[1]), log = TRUE)
  
  # Pre-evaluating useful quantities
  Xtilde  <- cbind(X0, X)
  XtildeTXtilde <- t(Xtilde) %*% Xtilde
  XtildeTy <- t(Xtilde) %*% y
  
  #### LOOPING ####
  for (iter in 2:nMCMC){
    
    if (verbose && iter %% 250 == 0) {
      message("Iteration: ", iter, "/", nMCMC)
    }
    
    # Current values
    gamma_curr   <- gamma_post[iter - 1, ]
    xi_curr <- sapply(1:p, function(j) as.integer(all(gamma_curr[B[[j]]] == 1)))
    A_curr <- (xi_curr==1)
    beta0_curr  <- beta0_post[iter - 1, ]
    b_curr <- b_post[iter - 1, ]
    beta_curr  <- b_curr*xi_curr
    sigma2_curr  <- sigma2_post[iter - 1]
    tau2_curr    <- tau2_post[iter - 1, ]
    omj_curr <- sapply(1:p, function(j) sum(w[B[[j]]] * tau2_curr[B[[j]]]))
    lambda2_curr <- lambda2_post[iter - 1]
    if(lambda2_curr==0){
      omj_curr <- rep(d0, p)
    }
    
    #### GAMMA AND BETA UPDATE ####
    # GAMMA UPDATE
    # Looping over indicators
    # Sets to compare
    # Evalauting log-marginals
    log_m_curr <- log_marginal_y(A_curr, sigma2_curr, omj_curr)
    for (g in 1:G) {
      
      # Sets to compare
      if(gamma_curr[g]==1){
        gamma1 <- gamma_curr
        log_m1 <- log_m_curr
        
        gamma0 <- gamma_curr
        gamma0[g] <- 0
        xi0 <- sapply(1:p, function(j) as.integer(all(gamma0[B[[j]]] == 1)))
        A0 <- (xi0==1)
      
        # Evaluating missing log-marginal
        log_m0 <- log_marginal_y(A0, sigma2_curr, omj_curr)
        
      }else{
        gamma0 <- gamma_curr
        log_m0 <- log_m_curr
        
        gamma1 <- gamma_curr
        gamma1[g] <- 1
        xi1 <- sapply(1:p, function(j) as.integer(all(gamma1[B[[j]]] == 1)))
        A1 <- (xi1==1)
      
        # Evaluating missing log-marginal
        log_m1 <- log_marginal_y(A1, sigma2_curr, omj_curr)
      }
      
      # Evaluating ratio components
      log_num <- log(pi_g[g]) + log_m1
      log_den <- log(exp(log(pi_g[g]) + log_m1) + exp(log(1 - pi_g[g]) + log_m0))
      
      # Posterior probability of inclusion
      post_prob_1 <- exp(log_num - log_den)
      post_prob_1 <- min(max(post_prob_1, 0), 1)
      
      # sample gamma_g
      gamma_curr[g] <- rbinom(1, size = 1, prob = post_prob_1)
      
      if(gamma_curr[g]==1){
        log_m_curr <- log_m1
      } else {
        log_m_curr <- log_m0
      }
    }
    # Update the chains
    gamma_post[iter, ] <- gamma_curr
    xi_curr <- sapply(1:p, function(j) as.integer(all(gamma_curr[B[[j]]] == 1)))
    xi_post[iter, ] <- xi_curr
    A_curr <- (xi_curr==1)
    n_A_curr <- sum(xi_curr)
    active_size[iter] <- n_A_curr
    
    # BETA UPDATE
    # Active set
    Xtilde_A_curr  <- as.matrix(Xtilde[, c(rep(T, p0), A_curr)])
    dtilde_A_curr  <- c(d0vec, omj_curr[A_curr])
    Dtilde_A_curr <- diag(dtilde_A_curr, nrow=length(dtilde_A_curr))
    Q_A <- XtildeTXtilde[c(rep(T, p0), A_curr), c(rep(T, p0), A_curr), drop=F] + Dtilde_A_curr
    sd_na <- sqrt(sigma2_curr/omj_curr[!A_curr])

    # Sample active coefficients
    RQ_A  <- chol(Q_A)
    mu_A_curr <- solve(Q_A, XtildeTy[c(rep(T, p0), A_curr), , drop=F])
    btilde_A_curr <- mu_A_curr + sqrt(sigma2_curr) * backsolve(RQ_A, rnorm(n = p0+n_A_curr))
    
    # Sample inactive coefficients
    nn_A_curr <- (p-n_A_curr)
    bn_A_curr <- rnorm(nn_A_curr, 0, sd_na)
    
    # Update the chains
    beta0_curr <- btilde_A_curr[1:p0]
    b_curr[A_curr] <- btilde_A_curr[-(1:p0)]
    b_curr[!A_curr] <- bn_A_curr
    beta_curr  <- b_curr*xi_curr
    
    beta0_post[iter, ] <- beta0_curr
    b_post[iter, ]  <- b_curr
    beta_post[iter, ]  <- beta_curr
    
    
    #### SIGMA2 UPDATE ####
    btilde_A_curr <- c(beta0_curr, b_curr[A_curr])
    mu_A_curr <- as.numeric(Xtilde_A_curr %*% btilde_A_curr)
    resid <- y - mu_A_curr
    shape_sigma2 <- alpha_sigma2 + (n + p0 + p) / 2
    dtilde_curr <- c(d0vec, omj_curr)
    btilde_curr <- c(beta0_curr, b_curr)
    rate_sigma2  <- eta_sigma2 + 0.5 * ( sum(resid^2) + sum(dtilde_curr * btilde_curr^2) )
    # Sampling
    sigma2_curr <- 1 / rgamma(n = 1, shape = shape_sigma2, rate = rate_sigma2)
    
    # Storing in the chain
    sigma2_post[iter] <- sigma2_curr
    
    ####  TAU2 UPDATE   ####
    for (g in 1:G) {
      
      # Penalized betas of current group
      bg <- b_curr[Groups[[g]]]
      norm2_bg <- sum(bg^2)
      
      if (norm2_bg > 1e-8){
        # Evaluate mean
        mu_g <- exp(0.5 * (log(lambda2_curr) + log(sigma2_curr) - log(w[g]) - log(norm2_bg)))
        
        # Sampling
        tau2_curr[g]  <- rinvgauss(n = 1, mean = mu_g, shape = lambda2_curr)
      }else{
        # limite ||beta_g|| -> 0 : InvGamma(1/2, scale = lambda2/2)
        tau2_curr[g]  <- 1/rgamma(n = 1, shape = 1/2, rate = lambda2_curr / 2)
      }
      
    }
    # Storing in the chain
    tau2_post[iter, ] <- tau2_curr
    

    ####  LAMBDA2 UPDATE   #### 
    # Gibbs update
    if(lambda2_update=="Gibbs"){
      
      lambda2_post[iter] <- rgamma(1, shape =  r + (p + G) / 2, rate = delta + 1 / 2 * sum(1 / tau2_curr))
      
    }else{
      
      # Fixed update
      lambda2_post[iter] <- lambda2_fixed
    }
    
    ####  Evaluating and storing Log-likelihood   #### 
    loglik_post[iter, ] <- dnorm(y, mean = mu_A_curr, sd = sqrt(sigma2_curr), log = TRUE)
    
  }
  
  #### FINALIZING CHAIN   #### 
  # Burning and thinning the chain
  burnin <- floor(burn * nMCMC)
  thin <- seq(from=burnin, to=nMCMC, by=thin)
  keep <- seq.int(from = burnin + 1, to = nMCMC)
  
  # Final chains
  gamma_post   <- gamma_post[keep, , drop = FALSE]
  xi_post      <- xi_post[keep, , drop = FALSE]
  active_post  <- active_size[keep]
  beta0_post   <- beta0_post[keep, , drop = FALSE]
  b_post       <- b_post[keep, , drop = FALSE]
  beta_post    <- beta_post[keep, , drop = FALSE]
  sigma2_post  <- sigma2_post[keep]
  tau2_post    <- tau2_post[keep, , drop = FALSE]
  lambda2_post <- lambda2_post[keep]
  
  # Calculate WAIC
  loglik_post <- loglik_post[keep, , drop = FALSE]

  lppd   <- sum(apply(loglik_post, 2, log_mean_exp))   # sum_i log E[p(y_i|theta)]
  p_waic <- sum(apply(loglik_post, 2, var))            # sum_i Var[log p(y_i|theta)]
  waic   <- -2 * (lppd - p_waic)
  
  # Evaluating some final summaries
  beta0_summary <- rbind(
    mean = apply(beta0_post, 2, mean),
    sd   = apply(beta0_post, 2, sd),
    q025 = apply(beta0_post, 2, quantile, probs = 0.025),
    q975 = apply(beta0_post, 2, quantile, probs = 0.975)
  )
  
  beta_summary <- rbind(
    mean = apply(beta_post, 2, mean),
    sd   = apply(beta_post, 2, sd),
    q025 = apply(beta_post, 2, quantile, probs = 0.025),
    q975 = apply(beta_post, 2, quantile, probs = 0.975)
  )
  
  res <- list()
  res$gamma                 <- gamma_post
  res$xi                    <- xi_post
  res$active_size           <- active_post
  res$beta0                 <- beta0_post
  res$b                     <- b_post
  res$beta                  <- beta_post
  res$sigma2                <- sigma2_post
  res$tau2                  <- tau2_post
  res$lambda2               <- lambda2_post
  res$beta0_summary         <- beta0_summary
  res$beta_summary          <- beta_summary
  res$group_membership      <- Groups
  res$coefficient_to_groups <- B
  res$loglik                <- loglik_post
  res$lppd                  <- lppd
  res$p_waic                <- p_waic
  res$waic                  <- waic
  
  return(res)
}

# MCMC for Bayesian Overlap Group Lasso -----------------------------------------------------------------

sampler_BOGL <- function(
    y,
    X0 = NULL,
    X,
    Groups,          
    weighted = TRUE, # penalty lambda * sum_g sqrt(m_g) ||beta_g||
    d0 = 1e-6,       # prior precision for non-penalized coefficients
    alpha_sigma2 = 1,
    eta_sigma2 = 1,
    beta0_init = rep(0, ncol(X0)),
    beta_init = rep(0, ncol(X)),
    sigma2_init = NULL,
    tau2_init = 1,
    lambda2_update = c("Gibbs", "fixed"),
    lambda2_fixed = NULL, # required only if lambda2_update=="fixed"
    r = 1,                # shape prior on lambda2, required only if lambda2_update!="fixed"
    delta = 1,            # rate prior on lambda2, required only if lambda2_update!="fixed"
    lambda2_init = 1,     # required only if lambda2_update!="fixed"
    nMCMC = 5000,
    burn = 0.2,
    thin = 1,
    verbose = TRUE){
  
 
  #### INITIALIZATION   #### 
  # Choose the lambda2 updating strategy
  lambda2_update <- match.arg(lambda2_update)
  if (lambda2_update == "fixed" && is.null(lambda2_fixed)) {
    stop("Con lambda2_update = 'fixed' devi fornire lambda2_fixed.")
  }
  
  #  Storing and checking data
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(X0)) {
    X0 <- matrix(1, nrow = n, ncol = 1)
  } else {
    X0 <- as.matrix(X0)
  }
  p0 <- ncol(X0)
  X <- as.matrix(X)
  p <- ncol(X)
  d0vec <- rep(d0, p0)
  
  #  Grouping structure
  G  <- length(Groups)
  B <- vector("list", p)
  for (j in 1:p) {
    B[[j]] <- which(sapply(Groups, function(gr) j %in% gr))
    if (length(B[[j]]) == 0) {
      stop(sprintf("Coefficient %d does not belong to any group.", j))
    }
  }
  # Weighting penalty by group size
  m_g <- unlist(lapply(Groups, length))
  w   <- if (weighted) m_g else rep(1, G)
  
  #  Chain storing objects
  beta0_post   <- matrix(NA, nrow = nMCMC, ncol = p0)
  beta_post    <- matrix(NA, nrow = nMCMC, ncol = p)
  sigma2_post  <- rep(NA, nMCMC)
  tau2_post    <- matrix(NA, nrow = nMCMC, ncol = G)
  lambda2_post <- rep(NA, nMCMC)
  loglik_post  <- matrix(NA, nrow = nMCMC, ncol = n)       # Pointwise loglikelihood
  
  #  Starting values
  beta0_post[1, ] <- beta0_init
  beta_post[1, ]  <- beta_init
  sigma2_post[1]  <- if (is.null(sigma2_init)) stats::var(y) else sigma2_init
  tau2_post[1, ]  <- rep(tau2_init, G)
  lambda2_post[1] <- if (lambda2_update == "fixed") lambda2_fixed else lambda2_init
  
  # Evaluating loglikelihood at starting value
  mu_init <- as.numeric(X0 %*% beta0_post[1, ] + X %*% beta_post[1, ])
  loglik_post[1, ] <- dnorm(y, mean = mu_init, sd = sqrt(sigma2_post[1]), log = TRUE)
  
  # Pre-evaluating useful quantities
  Xtilde  <- cbind(X0, X)
  XtildeTXtilde <- t(Xtilde) %*% Xtilde
  XtildeTy <- t(Xtilde) %*% y
  
  #### LOOPING   #### 
  for (iter in 2:nMCMC){
    
    if (verbose && iter %% 250 == 0) {
      message("Iteration: ", iter, "/", nMCMC)
    }
    
    # Current values
    beta0_curr  <- beta0_post[iter - 1, ]
    beta_curr  <- beta_post[iter - 1, ]
    sigma2_curr  <- sigma2_post[iter - 1]
    tau2_curr    <- tau2_post[iter - 1, ]
    omj_curr  <- sapply(1:p, function(j) sum(w[B[[j]]] * tau2_curr[B[[j]]]))
    lambda2_curr <- lambda2_post[iter - 1]
    if(lambda2_curr==0){
      omj_curr <- rep(d0, p)
    }
    
    ####  BETA UPDATE   #### 
    dtilde_curr   <- c(d0vec, omj_curr)
    Dtilde_curr   <- diag(dtilde_curr)
    # Conditional precision matrix
    Q_beta <- XtildeTXtilde + Dtilde_curr           
    RQ_beta <- chol(Q_beta)                         
    mu_beta <- solve(Q_beta, XtildeTy)              # Conditional mean
    
    # Sampling
    betatilde_curr <- mu_beta + sqrt(sigma2_curr) * backsolve(RQ_beta, rnorm(length(mu_beta)))
    
    # Splitting unpenalized and penalized
    beta0_curr <- betatilde_curr[1:p0]
    beta_curr  <- betatilde_curr[-c(1:p0)]
    # Storing in the chain
    beta0_post[iter, ] <- beta0_curr
    beta_post[iter, ]  <- beta_curr
    
    ####  SIGMA2 UPDATE   #### 
    mu_iter <- as.numeric(Xtilde %*% betatilde_curr)
    resid  <- y - mu_iter          # Evaluate residuals
    shape_sigma2 <- alpha_sigma2 + (n + p0 + p) / 2
    rate_sigma2  <- eta_sigma2 + 0.5 * sum(resid^2) + 0.5 * sum(dtilde_curr * betatilde_curr^2)
    # Sampling
    sigma2_curr <- 1 / rgamma(n = 1, shape = shape_sigma2, rate  = rate_sigma2)
    
    # Storing in the chain
    sigma2_post[iter] <- sigma2_curr

    ####  TAU2 UPDATE   #### 
    # Looping over dimensions
    for (g in 1:G) {
      
      # Penalized betas of current group
      betag <- beta_curr[Groups[[g]]]
      norm2_betag <- sum(betag^2)
      
      if (norm2_betag > 1e-8){
        # Evaluate mean
        mu_g <- exp(0.5 * (log(lambda2_curr) + log(sigma2_curr) - log(w[g]) - log(norm2_betag)))

        # Sampling
        tau2_curr[g]  <- rinvgauss(n = 1, mean = mu_g, shape = lambda2_curr)
      }else{
        # limite ||beta_g|| -> 0 : InvGamma(1/2, scale = lambda2/2)
        tau2_curr[g]  <- 1/rgamma(n = 1, shape = 1/2, rate = lambda2_curr / 2)
      }
    }
    
    tau2_post[iter, ] <- tau2_curr
    
    ####  LAMBDA2 UPDATE   #### 
    # Gibbs update
    if(lambda2_update=="Gibbs"){
      
      lambda2_post[iter] <- rgamma(1, shape =  r + (p+G) / 2, rate = delta + 1 / 2 * sum(1 / tau2_curr))
      
    
    } else {
      
      # Fixed update
      lambda2_post[iter] <- lambda2_fixed
    }
    
    
    ####  Evaluating and storing Log-likelihood   #### 
    loglik_post[iter, ] <- dnorm(y, mean = mu_iter, sd = sqrt(sigma2_curr), log = TRUE)
    
  }
  
  #### FINALIZING CHAIN   #### 
  # Burning and thinning the chain
  burnin <- floor(burn * nMCMC)
  thin <- seq(from=burnin, to=nMCMC, by=thin)
  keep <- seq.int(from = burnin + 1, to = nMCMC)
  
  # Final chains
  beta0_post   <- beta0_post[keep, , drop = FALSE]
  beta_post    <- beta_post[keep, , drop = FALSE]
  sigma2_post  <- sigma2_post[keep]
  tau2_post    <- tau2_post[keep, , drop = FALSE]
  lambda2_post <- lambda2_post[keep]
  
  # Calculate WAIC
  loglik_post <- loglik_post[keep, , drop = FALSE]
  
  lppd   <- sum(apply(loglik_post, 2, log_mean_exp))   # sum_i log E[p(y_i|theta)]
  p_waic <- sum(apply(loglik_post, 2, var))            # sum_i Var[log p(y_i|theta)]
  waic   <- -2 * (lppd - p_waic)
  
  
  # Evaluating some final summaries
  beta0_summary <- rbind(
    mean = apply(beta0_post, 2, mean),
    sd   = apply(beta0_post, 2, sd),
    q025 = apply(beta0_post, 2, quantile, probs = 0.025),
    q975 = apply(beta0_post, 2, quantile, probs = 0.975)
  )
  
  beta_summary <- rbind(
    mean = apply(beta_post, 2, mean),
    sd   = apply(beta_post, 2, sd),
    q025 = apply(beta_post, 2, quantile, probs = 0.025),
    q975 = apply(beta_post, 2, quantile, probs = 0.975)
  )
  
  # Creating output
  res <- list()
  res$beta0                 <- beta0_post
  res$beta                  <- beta_post
  res$sigma2                <- sigma2_post
  res$tau2                  <- tau2_post
  res$lambda2               <- lambda2_post
  res$beta0_summary         <- beta0_summary
  res$beta_summary          <- beta_summary
  res$group_membership      <- Groups
  res$coefficient_to_groups <- B
  res$weights               <- w
  res$lambda2_update        <- lambda2_update
  res$loglik                <- loglik_post
  res$lppd                  <- lppd
  res$p_waic                <- p_waic
  res$waic                  <- waic
  
  # Returning the output
  return(res)
}


# Wrapper -----------------------------------------------------------------


run_BOGL <- function(
    y,
    X0 = NULL,
    X,  
    Groups,
    weighted = TRUE,
    d0 = 1e-6,        # prior variance on unpenalized coefficients
    alpha_sigma2 = 1, # shape prior on sigma2
    eta_sigma2 = 1, # rate prior on sigma2
    beta0_init = rep(0, ncol(X0)),
    beta_init = rep(0, ncol(X)),
    sigma2_init = NULL,
    tau2_init = 1,
    SpikeandSlab = FALSE,
    pi_g = 0.5,             # only if SpikeandSlab = TRUE
    lambda2_update = c("Gibbs", "fixed"),
    lambda2_fixed = NULL, # required only if lambda2_update != "fixed"
    r = 1,            # shape prior on lambda2, required only if lambda2_update != "fixed"
    delta = 1,        # rate prior on lambda2, required only if lambda2_update != "fixed"
    lambda2_init = 1,
    nMCMC = 5000,
    burn = 0.2,
    thin = 1,
    verbose = TRUE){
  
  
  #### Choose lambda2 updating strategy ####
  lambda2_update <- match.arg(lambda2_update)
  # Check if a lambda2 is given when fixed is the strategy
  if (lambda2_update == "fixed" && is.null(lambda2_fixed)) {
    stop("Con lambda2_update = 'fixed' devi fornire lambda2_fixed.")
  }
  
  # Spike and Slab model
  if (SpikeandSlab == TRUE) {
    
    res <- sampler_BOGL_SpikeSlab(
      y             = y,
      X0            = X0,
      X             = X,
      Groups        = Groups,
      weighted      = weighted,
      d0            = d0,
      alpha_sigma2  = alpha_sigma2,
      eta_sigma2    = eta_sigma2,
      pi_g          = pi_g,
      beta0_init    = rep(0, ncol(X0)),
      beta_init     = rep(0, ncol(X)),
      tau2_init     = tau2_init,
      sigma2_init   = sigma2_init,
      lambda2_update = lambda2_update,
      lambda2_fixed = lambda2_fixed,
      r             = r,
      delta         = delta,
      lambda2_init  = lambda2_init,
      nMCMC         = nMCMC,
      burn          = burn,
      thin          = thin,
      verbose       = verbose
    )
    
  } else {
    
    res <- sampler_BOGL(
      y             = y,
      X0            = X0,
      X             = X,
      Groups        = Groups,
      weighted      = weighted,
      d0            = d0,
      alpha_sigma2  = alpha_sigma2,
      eta_sigma2    = eta_sigma2,
      beta0_init    = rep(0, ncol(X0)),
      beta_init     = rep(0, ncol(X)),
      tau2_init     = tau2_init,
      sigma2_init   = sigma2_init,
      lambda2_update = lambda2_update,
      lambda2_fixed = lambda2_fixed,
      r             = r,
      delta         = delta,
      lambda2_init  = lambda2_init,
      nMCMC         = nMCMC,
      burn          = burn,
      thin          = thin,
      verbose       = verbose
    )
  }
  
  res$SpikeandSlab <- SpikeandSlab
  
  return(res)
}
