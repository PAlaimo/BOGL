
# Packages ----------------------------------------------------------------

library(foreach)
library(doSNOW)
library(parallel)


# Scenario setting --------------------------------------------------------

scenarios <- c("random", "nested", "interactions")
rhos <- c(0, 0.8)
ex <- expand.grid(scenarios,rhos)
nr <- nrow(ex)

# Training set
prop_tr <- 0.7

# Simulation settings
n_rep <- 100
n <- 300

# Simulation and fit ------------------------------------------------------


# Drawing the seeds
seeds <- sample(1:1000, size = nr)
# save(seeds, file = paste0("seeds_scenario_", cv, ".RData"))

# Parallelization
cl <- makeCluster(nr, type = "SOCK")
registerDoSNOW(cl)


res <- foreach(cv = 1:nr, .packages = c("mvtnorm","mniw","cluster","glmnet"), 
               .export = c("ex" ,"n_rep", "n", "seeds"))%dopar%{
  # Load functions
  source("scenarios.R")
  source("functions.R")
  source("metrics.R")
  

  # Storing object
  out <- list()
  
  # Simulate data
  # COMMENT THIS IF DATA HAVE BEEN GENERATED AND MUST BE LOADED
  data <- generate_simulation_design(
    scenario    = ex[cv,1],
    rho         = ex[cv,2],
    n_rep       = n_rep,
    n           = n,
    p0          = 1,
    int0        = T,
    p_random    = 25,
    p_bilevel   = 16,   # multiplo di 8
    p_nested    = 10,
    noise_sd    = 1,
    seed        = seeds[cv]
  )
  # DECOMMENT THIS IF DATA HAVE BEEN GENERATED AND MUST BE LOADED
  # load("simdata/simdat.RData")
  # data <- simdat[[cv]]

  for(j in 1:n_rep){
    
    # To test
    # j <- 1
    
    # Simulated data
    y <- data[[j+1]]
    n <- length(y)
    X0 <- data[[1]]$X0
    p0 <- ncol(X0)
    X <- data[[1]]$X
    p <- ncol(X)
    
    # Setting and storing seeds, scenarios, repetition
    set.seed(seeds[cv] + j)
    save(seeds, cv, j, file = paste0("res_rep_", j, "_scenario_", cv, ".RData"))
    
    # Train-test split 70-30
    id_tr <- sample(1:n, floor(prop_tr * n))
    n_tr <- length(id_tr)
    id_te <- (1:n)[-id_tr]
    n_te <- n-n_tr
    X0_tr <- as.matrix(X0[id_tr, ])
    X_tr <- X[id_tr, ]
    y_tr <- y[id_tr]
    X0_te <- as.matrix(X0[id_te, ])
    X_te <- X[id_te, ]
    y_te <- y[id_te]
    
    # Final storing objects
    out_j <- list()
    
    out_j$y <- y
    out_j$id_tr <- id_tr
    out_j$id_te <- id_te

      
    # Storing true underlying parameters
    act   <- data[[1]]$active_idx
    inact <- data[[1]]$null_idx
    beta0_true <- data[[1]]$beta0
    beta_true <- data[[1]]$beta
    sigma2_true <- data[[1]]$noise_sd^2
    pg <- data[[1]]$Groups
    
    # Groups for standard lasso
    pg_lasso <- as.list(1:p)

    # Bi-level group creation
    pg_bilevel <- c(pg, pg_lasso)

    # Fit ---------------------------------------------------------------------

    # lambda hyperparameters
    r       <- 1
    delta   <- 1
    
    fit_reg <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg,
      weighted       = TRUE,
      SpikeandSlab   = FALSE,
      lambda2_update = "fixed",
      lambda2_fixed  = 0, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
    
    fit_blasso <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg_lasso,
      weighted       = TRUE,
      SpikeandSlab   = FALSE,
      lambda2_update = "Gibbs",
      r              = r,
      delta          = delta, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
  
    
    fit_bogl <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg,
      weighted       = TRUE,
      SpikeandSlab   = FALSE,
      lambda2_update = "Gibbs",
      lambda2_fixed  = NULL,
      r              = r,
      delta          = delta, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
    
    fit_bogl_bilevel <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg_bilevel,
      weighted       = TRUE,
      SpikeandSlab   = FALSE,
      lambda2_update = "Gibbs",
      lambda2_fixed  = NULL,
      r              = r,
      delta          = delta, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
    
    fit_ss <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg_lasso,
      weighted       = TRUE,
      SpikeandSlab   = TRUE,
      lambda2_update = "fixed",
      lambda2_fixed  = 0, 
      r              = r,
      delta          = delta, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
    
    fit_blasso_ss <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg_lasso,
      weighted       = TRUE,
      SpikeandSlab   = TRUE,
      lambda2_update = "Gibbs",
      lambda2_fixed  = NULL,
      r              = r,
      delta          = delta, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
    
    fit_bogl_ss <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg,
      weighted       = TRUE,
      SpikeandSlab   = TRUE,
      lambda2_update = "Gibbs",
      lambda2_fixed  = NULL,
      r              = r,
      delta          = delta, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
    
    fit_bogl_bilevel_ss <- run_BOGL(
      y              = y_tr,
      X0             = X0_tr,
      X              = X_tr,
      Groups         = pg_bilevel,
      weighted       = TRUE,
      SpikeandSlab   = TRUE,
      lambda2_update = "Gibbs",
      lambda2_fixed  = NULL,
      r              = r,
      delta          = delta, 
      nMCMC          = 10000,
      burn           = 0.2,
      verbose        = TRUE
    )
    
    # reg
    out_j$beta_reg <- apply(fit_reg$beta, 2, mean)
    out_j$waic_reg      <- fit_reg$waic
    out_j$rmse_beta_reg <- mean(apply(fit_reg$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_reg <- mean(apply(fit_reg$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_reg <- mean(apply(fit_reg$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_reg <- mean(apply(fit_reg$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    
    # blasso
    out_j$beta_blasso <- apply(fit_blasso$beta, 2, mean)
    out_j$waic_blasso      <- fit_blasso$waic
    out_j$rmse_beta_blasso <- mean(apply(fit_blasso$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_blasso <- mean(apply(fit_blasso$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_blasso <- mean(apply(fit_blasso$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_blasso <- mean(apply(fit_blasso$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    
    # bogl
    out_j$beta_bogl <- apply(fit_bogl$beta, 2, mean)
    out_j$waic_bogl      <- fit_bogl$waic
    out_j$rmse_beta_bogl <- mean(apply(fit_bogl$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_bogl <- mean(apply(fit_bogl$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_bogl <- mean(apply(fit_bogl$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_bogl <- mean(apply(fit_bogl$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    
    # bogl-bilevel
    out_j$beta_bogl_bilevel <- apply(fit_bogl_bilevel$beta, 2, mean)
    out_j$waic_bogl_bilevel      <- fit_bogl_bilevel$waic
    out_j$rmse_beta_bogl_bilevel <- mean(apply(fit_bogl_bilevel$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_bogl_bilevel <- mean(apply(fit_bogl_bilevel$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_bogl_bilevel <- mean(apply(fit_bogl_bilevel$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_bogl_bilevel <- mean(apply(fit_bogl_bilevel$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    
    # S&S
    out_j$beta_mean_xi1_ss <- apply(fit_ss$beta, 2, function(z) mean(z[z!=0]))
    out_j$xi_mean_ss <- apply(fit_ss$xi, 2, mean)
    out_j$waic_ss      <- fit_ss$waic
    out_j$rmse_beta_ss <- mean(apply(fit_ss$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_ss <- mean(apply(fit_ss$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_ss <- mean(apply(fit_ss$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_ss <- mean(apply(fit_ss$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    out_j$prop_act_ss <- mean(apply(fit_ss$xi, 1, function(x) mean(x[act])))
    out_j$prop_inact_ss <- mean(apply(fit_ss$xi, 1, function(x) mean(x[inact])))
    
    # blasso S&S
    out_j$beta_mean_xi1_blasso_ss <- apply(fit_blasso_ss$beta, 2, function(z) mean(z[z!=0]))
    out_j$xi_mean_blasso_ss <- apply(fit_blasso_ss$xi, 2, mean)
    out_j$waic_blasso_ss      <- fit_blasso_ss$waic
    out_j$rmse_beta_blasso_ss <- mean(apply(fit_blasso_ss$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_blasso_ss <- mean(apply(fit_blasso_ss$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_blasso_ss <- mean(apply(fit_blasso_ss$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_blasso_ss <- mean(apply(fit_blasso_ss$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    out_j$prop_act_blasso_ss <- mean(apply(fit_blasso_ss$xi, 1, function(x) mean(x[act])))
    out_j$prop_inact_blasso_ss <- mean(apply(fit_blasso_ss$xi, 1, function(x) mean(x[inact])))
    
    # bogl S&S
    out_j$beta_mean_xi1_bogl_ss <- apply(fit_bogl_ss$beta, 2, function(z) mean(z[z!=0]))
    out_j$xi_mean_bogl_ss <- apply(fit_bogl_ss$xi, 2, mean)
    out_j$waic_bogl_ss      <- fit_bogl_ss$waic
    out_j$rmse_beta_bogl_ss <- mean(apply(fit_bogl_ss$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_bogl_ss <- mean(apply(fit_bogl_ss$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_bogl_ss <- mean(apply(fit_bogl_ss$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_bogl_ss <- mean(apply(fit_bogl_ss$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    out_j$prop_act_bogl_ss <- mean(apply(fit_bogl_ss$xi, 1, function(x) mean(x[act])))
    out_j$prop_inact_bogl_ss <- mean(apply(fit_bogl_ss$xi, 1, function(x) mean(x[inact])))
    
    # bogl-bilevel S&S
    out_j$beta_mean_xi1_bogl_bilevel_ss <- apply(fit_bogl_bilevel_ss$beta, 2, function(z) mean(z[z!=0]))
    out_j$xi_mean_bogl_bilevel_ss <- apply(fit_bogl_bilevel_ss$xi, 2, mean)
    out_j$waic_bogl_bilevel_ss      <- fit_bogl_bilevel_ss$waic
    out_j$rmse_beta_bogl_bilevel_ss <- mean(apply(fit_bogl_bilevel_ss$beta, 1, function(z) mean((beta_true - z)^2)))
    out_j$rmse_beta_active_bogl_bilevel_ss <- mean(apply(fit_bogl_bilevel_ss$beta, 1, function(z) mean((beta_true[act] - z[act])^2)))
    out_j$rmse_beta_inactive_bogl_bilevel_ss <- mean(apply(fit_bogl_bilevel_ss$beta, 1, function(z) mean((beta_true[inact] - z[inact])^2)))
    out_j$mse_pred_bogl_bilevel_ss <- mean(apply(fit_bogl_bilevel_ss$beta, 1, function(z) mean((y_te - X_te %*% z)^2)))
    out_j$prop_act_bogl_bilevel_ss <- mean(apply(fit_bogl_bilevel_ss$xi, 1, function(x) mean(x[act])))
    out_j$prop_inact_bogl_bilevel_ss <- mean(apply(fit_bogl_bilevel_ss$xi, 1, function(x) mean(x[inact])))
    
    out_j$selection_metrics_ss           <- var_sel_metrics(fit_ss, data[[1]]$beta, null_groups = data[[1]]$null_groups)
    out_j$selection_metrics_bogl_ss           <- var_sel_metrics(fit_bogl_ss, data[[1]]$beta, null_groups = data[[1]]$null_groups)
    out_j$selection_metrics_blasso_ss         <- var_sel_metrics(fit_blasso_ss, data[[1]]$beta)
    out_j$selection_metrics_bogl_bilevel_ss  <- var_sel_metrics(fit_bogl_bilevel_ss, data[[1]]$beta, null_groups = data[[1]]$null_groups)
    
    
    out_j$rank <- rank(c(fit_reg$waic, fit_blasso$waic, fit_bogl$waic, fit_bogl_bilevel$waic,
                         fit_ss$waic, fit_blasso_ss$waic, fit_bogl_ss$waic, fit_bogl_bilevel_ss$waic))
    
    out[[j]] <- out_j
  }
  
  out[[n_rep+1]] <- data
  
  return(out)
}


res$structure <- ex
  
stopCluster(cl)
save(res,file="results.RData")
