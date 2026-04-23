# =============================================================================
# Ensemble Learning and GFR Framework: Complete Reproducible Code
# Mili, Argoubi & Mili (2025)
# R version 4.1.2
# =============================================================================

# --- 1. PACKAGE SETUP --------------------------------------------------------

required_packages <- c(
  "mclust",        # 5.4.9  — GMM estimation
  "glmnet",        # 4.1.4  — Ridge regularization
  "rpart",         # 4.1.15 — CART
  "randomForest",  # 4.7.1  — Random forests
  "xgboost",       # 1.6.0.1— XGBoost
  "MASS",          # 7.3.55 — GLM support
  "Matrix",        # 1.4.0  — Sparse matrices
  "ggplot2",       # 3.3.5  — Visualization
  "dplyr"          # 1.0.8  — Data manipulation
)

installed <- rownames(installed.packages())
to_install <- required_packages[!required_packages %in% installed]
if (length(to_install) > 0) install.packages(to_install)
lapply(required_packages, library, character.only = TRUE)

set.seed(42)

# =============================================================================
# 2. CORE GFR FUNCTIONS
# =============================================================================

# --- 2.1 Fit GMM to habitat availability -------------------------------------

fit_gmm <- function(availability_data, max_k = 6) {
  # Fits Gaussian mixture model to landscape availability distribution
  # Returns optimal GMM selected by BIC
  # availability_data: matrix [n_cells x n_habitat_vars]
  
  bic_vals <- mclust::mclustBIC(availability_data, G = 1:max_k)
  mod      <- mclust::Mclust(availability_data, x = bic_vals)
  
  list(
    weights = mod$parameters$pro,
    means   = mod$parameters$mean,       # [n_vars x K]
    vars    = mod$parameters$variance$sigmasq,  # diagonal variances
    K       = mod$G,
    model   = mod
  )
}

# --- 2.2 Polynomial GFR moments ----------------------------------------------

compute_poly_moments <- function(gmm, max_order = 12) {
  # Computes moments E[X^m]_b from GMM for polynomial GFR
  # Returns list of moments per variable per order
  
  K       <- gmm$K
  n_vars  <- length(gmm$means[, 1])
  moments <- list()
  
  for (j in seq_len(n_vars)) {
    moments[[j]] <- numeric(max_order)
    for (m in seq_len(max_order)) {
      # Gaussian moment: sum_k w_k * E[X^m | component k]
      # Uses closed-form raw moments of Gaussian
      raw_moment <- 0
      for (k in seq_len(K)) {
        mu_k  <- gmm$means[j, k]
        sig_k <- sqrt(gmm$vars[k])
        # Raw moment of N(mu, sigma^2): sum_{i=0}^{floor(m/2)} C(m,2i)*(2i-1)!!*sigma^(2i)*mu^(m-2i)
        raw_moment <- raw_moment + gmm$weights[k] *
          sum(sapply(0:floor(m / 2), function(i) {
            choose(m, 2 * i) *
              prod(seq(1, 2 * i - 1, by = 2)) *
              sig_k^(2 * i) *
              mu_k^(m - 2 * i)
          }))
      }
      moments[[j]][m] <- raw_moment
    }
  }
  moments
}

# --- 2.3 RBF convolution integral I_{j,m,b} ----------------------------------

compute_rbf_integral <- function(gmm, centers, bandwidths) {
  # Closed-form convolution of RBF kernels with GMM availability
  # centers:    [n_vars x M_j] matrix of RBF centers
  # bandwidths: [n_vars x M_j] matrix of RBF bandwidths
  # Returns:    [n_vars x M_j] matrix of integral values
  
  K      <- gmm$K
  n_vars <- nrow(centers)
  M_j    <- ncol(centers)
  I      <- matrix(0, nrow = n_vars, ncol = M_j)
  
  for (j in seq_len(n_vars)) {
    for (m in seq_len(M_j)) {
      xi_jm  <- centers[j, m]
      sig_jm <- bandwidths[j, m]
      val    <- 0
      for (k in seq_len(K)) {
        mu_jk  <- gmm$means[j, k]
        sig_jk <- sqrt(gmm$vars[k])
        denom  <- sqrt(sig_jk^2 + sig_jm^2)
        val    <- val + gmm$weights[k] * (sig_jm / denom) *
          exp(-(mu_jk - xi_jm)^2 / (2 * (sig_jk^2 + sig_jm^2)))
      }
      I[j, m] <- val
    }
  }
  I
}

# --- 2.4 Build design matrix -------------------------------------------------

build_design_matrix <- function(local_x, gmm_list, type = "poly",
                                max_order = 5, centers = NULL, bw = NULL) {
  # Constructs design matrix X for GLM fitting
  # local_x:   [n_obs x n_vars] local habitat values
  # gmm_list:  list of GMMs, one per landscape block
  # type:      "poly" or "rbf"
  
  n_obs   <- nrow(local_x)
  n_vars  <- ncol(local_x)
  n_land  <- length(gmm_list)
  
  if (type == "poly") {
    # Feature: local_x_j * E[X_j^m]_b — landscape moment weighted local covariate
    X_list <- lapply(seq_len(n_land), function(b) {
      moments <- compute_poly_moments(gmm_list[[b]], max_order)
      do.call(cbind, lapply(seq_len(n_vars), function(j) {
        outer(local_x[, j], moments[[j]][seq_len(max_order)])
      }))
    })
  } else {
    # Feature: local RBF evaluated at x_j, scaled by integral I_{j,m,b}
    if (is.null(centers)) stop("Centers required for RBF design matrix.")
    X_list <- lapply(seq_len(n_land), function(b) {
      I_jmb <- compute_rbf_integral(gmm_list[[b]], centers, bw)
      do.call(cbind, lapply(seq_len(n_vars), function(j) {
        rbf_vals <- sapply(seq_len(ncol(centers)), function(m) {
          exp(-(local_x[, j] - centers[j, m])^2 / (2 * bw[j, m]^2))
        })
        sweep(rbf_vals, 2, I_jmb[j, ], "*")
      }))
    })
  }
  
  # Stack all landscape-specific matrices
  do.call(rbind, X_list)
}

# --- 2.5 Place RBF centers and bandwidths ------------------------------------

setup_rbf <- function(habitat_data, n_basis = 5) {
  # Places RBF centers at quantiles; bandwidths = max adjacent spacing
  n_vars  <- ncol(habitat_data)
  probs   <- seq(0, 1, length.out = n_basis + 2)[c(-1, -(n_basis + 2))]
  centers <- matrix(0, nrow = n_vars, ncol = n_basis)
  bw      <- matrix(0, nrow = n_vars, ncol = n_basis)
  
  for (j in seq_len(n_vars)) {
    q           <- quantile(habitat_data[, j], probs = probs)
    centers[j, ] <- q
    spacings     <- diff(c(min(habitat_data[, j]), q, max(habitat_data[, j])))
    bw[j, ]      <- rep(max(spacings), n_basis)
  }
  list(centers = centers, bandwidths = bw)
}

# =============================================================================
# 3. MODEL FITTING FUNCTIONS
# =============================================================================

# --- 3.1 Standard GLM (baseline) ---------------------------------------------

fit_glm <- function(y, X, family = "poisson") {
  df  <- data.frame(y = y, X)
  mod <- glm(y ~ ., data = df, family = family)
  mod
}

# --- 3.2 GFR / RBF-GFR (unregularized) --------------------------------------

fit_gfr <- function(y, X_design, family = "poisson") {
  # Fits GFR model via IRLS (standard GLM on GFR design matrix)
  df  <- data.frame(y = y, X_design)
  mod <- glm(y ~ ., data = df, family = family)
  mod
}

# --- 3.3 Regularized GFR (Ridge) ---------------------------------------------

fit_gfr_ridge <- function(y, X_design, family = "poisson",
                          lambda_grid = 10^seq(-4, 2, length.out = 100)) {
  # Selects lambda by BIC using effective df from singular values
  
  fit_bic <- function(lam) {
    mod   <- glmnet::glmnet(X_design, y, family = family,
                            alpha = 0, lambda = lam)
    # Effective df via singular values
    svd_X <- svd(X_design)$d
    p_eff <- sum(svd_X^2 / (svd_X^2 + lam))
    ll    <- -deviance(mod) / 2   # approximate
    bic   <- -2 * ll + log(length(y)) * p_eff
    bic
  }
  
  bic_vals   <- sapply(lambda_grid, fit_bic)
  best_lam   <- lambda_grid[which.min(bic_vals)]
  final_mod  <- glmnet::glmnet(X_design, y, family = family,
                               alpha = 0, lambda = best_lam)
  list(model = final_mod, lambda = best_lam)
}

# --- 3.4 CART + GFR ----------------------------------------------------------

fit_cart_gfr <- function(y, X_local, X_gfr, family = "poisson") {
  # Grows CART on local covariates; fits GFR in each leaf
  
  df_tree <- data.frame(y = y, X_local)
  tree    <- rpart::rpart(
    y ~ ., data = df_tree,
    method  = ifelse(family == "binomial", "class", "poisson"),
    control = rpart::rpart.control(cp = 0, minsplit = 20)
  )
  
  # Prune by 1-SE rule
  cp_table  <- tree$cptable
  min_idx   <- which.min(cp_table[, "xerror"])
  threshold <- cp_table[min_idx, "xerror"] + cp_table[min_idx, "xstd"]
  best_cp   <- cp_table[cp_table[, "xerror"] <= threshold, "CP"][1]
  tree      <- rpart::prune(tree, cp = best_cp)
  
  # Assign leaf IDs
  leaf_ids  <- rpart::predict(tree, df_tree, type = "matrix")[, 1]
  leaves    <- unique(leaf_ids)
  
  # Fit GFR model within each leaf
  leaf_models <- lapply(leaves, function(l) {
    idx <- which(leaf_ids == l)
    if (length(idx) < 10) return(NULL)
    fit_gfr(y[idx], X_gfr[idx, , drop = FALSE], family)
  })
  names(leaf_models) <- as.character(leaves)
  
  list(tree = tree, leaf_models = leaf_models, leaves = leaves)
}

# --- 3.5 Random Forest + GFR -------------------------------------------------

fit_rf_gfr <- function(y, X_local, X_gfr, family = "poisson",
                       n_trees = 500) {
  # RF partitions data; GFR fitted in each leaf node
  
  is_class <- family == "binomial"
  y_rf     <- if (is_class) as.factor(y) else y
  mtry_val <- if (is_class) floor(ncol(X_local) / 3) else
    floor(sqrt(ncol(X_local)))
  
  rf <- randomForest::randomForest(
    x     = X_local,
    y     = y_rf,
    ntree = n_trees,
    mtry  = mtry_val,
    keep.inbag = TRUE
  )
  
  # Extract terminal node IDs for each observation
  node_matrix <- attr(randomForest::predict(rf, X_local, nodes = TRUE),
                      "nodes")  # [n_obs x n_trees]
  
  # For each tree, fit GFR per leaf
  tree_models <- vector("list", n_trees)
  for (t in seq_len(n_trees)) {
    nodes_t     <- node_matrix[, t]
    unique_nodes <- unique(nodes_t)
    node_fits   <- lapply(unique_nodes, function(nd) {
      idx <- which(nodes_t == nd)
      if (length(idx) < 5) return(NULL)
      tryCatch(
        fit_gfr(y[idx], X_gfr[idx, , drop = FALSE], family),
        error = function(e) NULL
      )
    })
    names(node_fits)  <- as.character(unique_nodes)
    tree_models[[t]]  <- list(nodes = nodes_t, fits = node_fits)
  }
  
  list(rf = rf, tree_models = tree_models)
}

# --- 3.6 XGBoost + GFR -------------------------------------------------------

fit_xgb_gfr <- function(y, X_local, X_gfr, family = "poisson",
                        iter_grid = c(2, 5, 10, 15, 20, 40, 80,
                                      100, 200, 300, 400, 500),
                        n_folds = 5) {
  # Nested CV to select optimal iterations; GFR fitted in leaf nodes
  
  obj     <- if (family == "binomial") "binary:logistic" else "count:poisson"
  dtrain  <- xgboost::xgb.DMatrix(data = as.matrix(X_local), label = y)
  
  # Nested CV for iteration selection
  fold_ids  <- cut(seq_along(y), breaks = n_folds, labels = FALSE)
  cv_scores <- matrix(NA, nrow = n_folds, ncol = length(iter_grid))
  
  for (f in seq_len(n_folds)) {
    idx_val   <- which(fold_ids == f)
    idx_train <- setdiff(seq_along(y), idx_val)
    # Split train into 80/20 tuning split
    n_tune    <- floor(0.8 * length(idx_train))
    idx_tune  <- idx_train[seq_len(n_tune)]
    idx_vali  <- idx_train[(n_tune + 1):length(idx_train)]
    
    d_tune  <- xgboost::xgb.DMatrix(as.matrix(X_local[idx_tune, ]), label = y[idx_tune])
    d_vali  <- xgboost::xgb.DMatrix(as.matrix(X_local[idx_vali, ]), label = y[idx_vali])
    
    for (i_iter in seq_along(iter_grid)) {
      mod_tmp <- xgboost::xgboost(
        data      = d_tune,
        nrounds   = iter_grid[i_iter],
        objective = obj,
        verbose   = 0
      )
      preds <- predict(mod_tmp, d_vali)
      r2    <- 1 - sum((y[idx_vali] - preds)^2) /
        sum((y[idx_vali] - mean(y[idx_vali]))^2)
      cv_scores[f, i_iter] <- r2
    }
  }
  
  best_iter <- iter_grid[which.max(apply(cv_scores, 2, median))]
  
  # Fit final XGBoost with best iterations
  xgb_final <- xgboost::xgboost(
    data      = dtrain,
    nrounds   = best_iter,
    objective = obj,
    verbose   = 0
  )
  
  # Extract leaf assignments and fit GFR per leaf
  leaf_preds  <- predict(xgb_final, dtrain, predleaf = TRUE)  # [n_obs x n_trees]
  n_trees_xgb <- ncol(leaf_preds)
  
  tree_models_xgb <- vector("list", n_trees_xgb)
  for (t in seq_len(n_trees_xgb)) {
    nodes_t      <- leaf_preds[, t]
    unique_nodes <- unique(nodes_t)
    node_fits    <- lapply(unique_nodes, function(nd) {
      idx <- which(nodes_t == nd)
      if (length(idx) < 5) return(NULL)
      tryCatch(
        fit_gfr(y[idx], X_gfr[idx, , drop = FALSE], family),
        error = function(e) NULL
      )
    })
    names(node_fits)       <- as.character(unique_nodes)
    tree_models_xgb[[t]]  <- list(nodes = nodes_t, fits = node_fits)
  }
  
  list(xgb = xgb_final, tree_models = tree_models_xgb, best_iter = best_iter)
}

# =============================================================================
# 4. PREDICTION FUNCTIONS
# =============================================================================

predict_glm_gfr <- function(model, X_new, type = "response") {
  predict(model, newdata = data.frame(X_new), type = type)
}

predict_ridge_gfr <- function(model_obj, X_new) {
  as.numeric(predict(model_obj$model, newx = X_new, type = "response"))
}

predict_rf_gfr <- function(model_obj, X_local_new, X_gfr_new) {
  rf          <- model_obj$rf
  tree_models <- model_obj$tree_models
  n_trees     <- length(tree_models)
  n_obs       <- nrow(X_local_new)
  
  pred_matrix <- matrix(NA, nrow = n_obs, ncol = n_trees)
  node_mat    <- attr(
    randomForest::predict(rf, X_local_new, nodes = TRUE), "nodes"
  )
  
  for (t in seq_len(n_trees)) {
    nodes_t <- node_mat[, t]
    for (i in seq_len(n_obs)) {
      nd      <- as.character(nodes_t[i])
      fit_nd  <- tree_models[[t]]$fits[[nd]]
      if (!is.null(fit_nd)) {
        pred_matrix[i, t] <- predict(
          fit_nd,
          newdata = data.frame(t(X_gfr_new[i, ])),
          type = "response"
        )
      }
    }
  }
  rowMeans(pred_matrix, na.rm = TRUE)
}

predict_xgb_gfr <- function(model_obj, X_local_new, X_gfr_new) {
  xgb         <- model_obj$xgb
  tree_models <- model_obj$tree_models
  n_trees     <- length(tree_models)
  n_obs       <- nrow(X_local_new)
  
  dtest       <- xgboost::xgb.DMatrix(as.matrix(X_local_new))
  leaf_preds  <- predict(xgb, dtest, predleaf = TRUE)
  
  pred_matrix <- matrix(NA, nrow = n_obs, ncol = n_trees)
  for (t in seq_len(n_trees)) {
    nodes_t <- leaf_preds[, t]
    for (i in seq_len(n_obs)) {
      nd     <- as.character(nodes_t[i])
      fit_nd <- tree_models[[t]]$fits[[nd]]
      if (!is.null(fit_nd)) {
        pred_matrix[i, t] <- predict(
          fit_nd,
          newdata = data.frame(t(X_gfr_new[i, ])),
          type = "response"
        )
      }
    }
  }
  rowMeans(pred_matrix, na.rm = TRUE)
}

# =============================================================================
# 5. PERFORMANCE METRICS
# =============================================================================

compute_r2 <- function(y_obs, y_pred) {
  ss_res <- sum((y_obs - y_pred)^2)
  ss_tot <- sum((y_obs - mean(y_obs))^2)
  1 - ss_res / ss_tot
}

compute_r2_dev <- function(y_obs, y_pred) {
  # Deviance-based R² for count data (Cameron & Windmeijer, 1996)
  y_pred  <- pmax(y_pred, 1e-10)
  y_bar   <- mean(y_obs)
  
  dev_mod  <- sum(y_obs * log(y_obs / y_pred) - (y_obs - y_pred),
                  na.rm = TRUE)
  dev_null <- sum(y_obs * log(y_obs / y_bar), na.rm = TRUE)
  1 - dev_mod / dev_null
}

compute_mad <- function(r2_vec) {
  median(abs(r2_vec - median(r2_vec))) * 1.4826
}

# =============================================================================
# 6. BLOCK CROSS-VALIDATION FRAMEWORK
# =============================================================================

run_block_cv <- function(data, block_var, y_var, habitat_vars,
                         family = "poisson", n_basis = 5,
                         max_poly_order = 12, n_trees_rf = 500) {
  # Main benchmark function
  # data:         data.frame with all variables
  # block_var:    column name identifying landscape/colony/pack blocks
  # y_var:        response variable name
  # habitat_vars: character vector of habitat covariate names
  
  blocks    <- unique(data[[block_var]])
  n_blocks  <- length(blocks)
  
  # Storage for results
  models    <- c("GLM", "GFR", "RBF-GFR", "Reg-GFR", "Reg-RBF-GFR",
                 "GFR-CART", "RBF-GFR-CART", "GFR-RF", "RBF-GFR-RF",
                 "GFR-XGBoost", "RBF-GFR-XGBoost")
  results   <- matrix(NA, nrow = n_blocks, ncol = length(models))
  colnames(results) <- models
  
  for (b_idx in seq_len(n_blocks)) {
    test_block  <- blocks[b_idx]
    train_data  <- data[data[[block_var]] != test_block, ]
    test_data   <- data[data[[block_var]] == test_block, ]
    
    y_train     <- train_data[[y_var]]
    y_test      <- test_data[[y_var]]
    X_train     <- as.matrix(train_data[, habitat_vars])
    X_test      <- as.matrix(test_data[, habitat_vars])
    
    # Fit GMMs per training block
    train_blocks <- unique(train_data[[block_var]])
    gmm_list     <- lapply(train_blocks, function(bl) {
      block_data <- train_data[train_data[[block_var]] == bl, habitat_vars]
      fit_gmm(as.matrix(block_data))
    })
    
    # Select polynomial order by BIC
    best_poly <- 1
    best_bic  <- Inf
    for (ord in 1:max_poly_order) {
      tryCatch({
        X_des <- build_design_matrix(X_train, gmm_list, type = "poly",
                                     max_order = ord)
        df    <- data.frame(y = y_train, X_des)
        mod   <- glm(y ~ ., data = df, family = family)
        bic   <- BIC(mod)
        if (bic < best_bic) {
          best_bic  <- bic
          best_poly <- ord
        }
      }, error = function(e) NULL)
    }
    
    # Select RBF basis count by BIC
    best_basis <- 1
    best_bic   <- Inf
    rbf_setup  <- NULL
    for (nb in 1:12) {
      tryCatch({
        rs    <- setup_rbf(X_train, n_basis = nb)
        X_des <- build_design_matrix(X_train, gmm_list, type = "rbf",
                                     centers = rs$centers, bw = rs$bandwidths)
        df    <- data.frame(y = y_train, X_des)
        mod   <- glm(y ~ ., data = df, family = family)
        bic   <- BIC(mod)
        if (bic < best_bic) {
          best_bic   <- bic
          best_basis <- nb
          rbf_setup  <- rs
        }
      }, error = function(e) NULL)
    }
    if (is.null(rbf_setup)) rbf_setup <- setup_rbf(X_train, n_basis = 1)
    
    # Build final design matrices
    X_poly_train <- build_design_matrix(X_train, gmm_list, "poly",
                                        max_order = best_poly)
    X_rbf_train  <- build_design_matrix(X_train, gmm_list, "rbf",
                                        centers = rbf_setup$centers,
                                        bw      = rbf_setup$bandwidths)
    
    # Test GMM (single block for test landscape)
    gmm_test     <- list(fit_gmm(X_test))
    X_poly_test  <- build_design_matrix(X_test, gmm_test, "poly",
                                        max_order = best_poly)
    X_rbf_test   <- build_design_matrix(X_test, gmm_test, "rbf",
                                        centers = rbf_setup$centers,
                                        bw      = rbf_setup$bandwidths)
    
    # ----- Fit all 11 models -------------------------------------------------
    
    # 1. GLM
    tryCatch({
      m_glm <- fit_glm(y_train, X_train, family)
      p     <- predict(m_glm, newdata = data.frame(X_test), type = "response")
      results[b_idx, "GLM"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 2. GFR (polynomial, unregularized)
    tryCatch({
      m_gfr <- fit_gfr(y_train, X_poly_train, family)
      p     <- predict(m_gfr, newdata = data.frame(X_poly_test), type = "response")
      results[b_idx, "GFR"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 3. RBF-GFR (unregularized)
    tryCatch({
      m_rbf <- fit_gfr(y_train, X_rbf_train, family)
      p     <- predict(m_rbf, newdata = data.frame(X_rbf_test), type = "response")
      results[b_idx, "RBF-GFR"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 4. Reg-GFR (ridge, polynomial)
    tryCatch({
      m_reg_gfr <- fit_gfr_ridge(y_train, X_poly_train, family)
      p         <- predict_ridge_gfr(m_reg_gfr, X_poly_test)
      results[b_idx, "Reg-GFR"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 5. Reg-RBF-GFR (ridge, RBF)
    tryCatch({
      m_reg_rbf <- fit_gfr_ridge(y_train, X_rbf_train, family)
      p         <- predict_ridge_gfr(m_reg_rbf, X_rbf_test)
      results[b_idx, "Reg-RBF-GFR"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 6. GFR-CART
    tryCatch({
      m_cart_gfr <- fit_cart_gfr(y_train, X_train, X_poly_train, family)
      leaf_test  <- predict(m_cart_gfr$tree,
                            data.frame(X_test), type = "matrix")[, 1]
      preds <- sapply(seq_along(y_test), function(i) {
        nd  <- as.character(leaf_test[i])
        fit <- m_cart_gfr$leaf_models[[nd]]
        if (is.null(fit)) return(NA)
        predict(fit, newdata = data.frame(t(X_poly_test[i, ])),
                type = "response")
      })
      results[b_idx, "GFR-CART"] <- compute_r2(y_test, preds)
    }, error = function(e) NULL)
    
    # 7. RBF-GFR-CART
    tryCatch({
      m_cart_rbf <- fit_cart_gfr(y_train, X_train, X_rbf_train, family)
      leaf_test  <- predict(m_cart_rbf$tree,
                            data.frame(X_test), type = "matrix")[, 1]
      preds <- sapply(seq_along(y_test), function(i) {
        nd  <- as.character(leaf_test[i])
        fit <- m_cart_rbf$leaf_models[[nd]]
        if (is.null(fit)) return(NA)
        predict(fit, newdata = data.frame(t(X_rbf_test[i, ])),
                type = "response")
      })
      results[b_idx, "RBF-GFR-CART"] <- compute_r2(y_test, preds)
    }, error = function(e) NULL)
    
    # 8. GFR-RF
    tryCatch({
      m_rf_gfr <- fit_rf_gfr(y_train, X_train, X_poly_train, family, n_trees_rf)
      p        <- predict_rf_gfr(m_rf_gfr, X_test, X_poly_test)
      results[b_idx, "GFR-RF"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 9. RBF-GFR-RF
    tryCatch({
      m_rf_rbf <- fit_rf_gfr(y_train, X_train, X_rbf_train, family, n_trees_rf)
      p        <- predict_rf_gfr(m_rf_rbf, X_test, X_rbf_test)
      results[b_idx, "RBF-GFR-RF"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 10. GFR-XGBoost
    tryCatch({
      m_xgb_gfr <- fit_xgb_gfr(y_train, X_train, X_poly_train, family)
      p         <- predict_xgb_gfr(m_xgb_gfr, X_test, X_poly_test)
      results[b_idx, "GFR-XGBoost"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    # 11. RBF-GFR-XGBoost
    tryCatch({
      m_xgb_rbf <- fit_xgb_gfr(y_train, X_train, X_rbf_train, family)
      p         <- predict_xgb_gfr(m_xgb_rbf, X_test, X_rbf_test)
      results[b_idx, "RBF-GFR-XGBoost"] <- compute_r2(y_test, p)
    }, error = function(e) NULL)
    
    cat(sprintf("Block %d / %d complete\n", b_idx, n_blocks))
  }
  
  results
}

# =============================================================================
# 7. STATISTICAL COMPARISON (Demsar, 2006)
# =============================================================================

run_friedman_test <- function(results_matrix) {
  # Friedman test across models (columns) over blocks (rows)
  # results_matrix: [n_blocks x n_models] R² values
  
  # Convert to ranks within each block (row)
  ranked <- t(apply(results_matrix, 1, rank, na.last = "keep"))
  
  n  <- nrow(ranked)    # number of datasets/blocks
  k  <- ncol(ranked)    # number of models
  
  R_j    <- colMeans(ranked, na.rm = TRUE)
  chi_F  <- (12 * n) / (k * (k + 1)) *
    (sum(R_j^2) - k * (k + 1)^2 / 4)
  
  p_val  <- pchisq(chi_F, df = k - 1, lower.tail = FALSE)
  
  cat(sprintf("Friedman test: chi²_F = %.3f, df = %d, p = %.4f\n",
              chi_F, k - 1, p_val))
  
  # Wilcoxon post-hoc for top-3 pairs
  top3    <- names(sort(colMeans(results_matrix, na.rm = TRUE),
                        decreasing = TRUE)[1:3])
  cat("\nWilcoxon signed-rank post-hoc (top 3 models):\n")
  pairs   <- combn(top3, 2)
  for (i in seq_len(ncol(pairs))) {
    a   <- pairs[1, i]
    b   <- pairs[2, i]
    wt  <- wilcox.test(results_matrix[, a], results_matrix[, b],
                       paired = TRUE, exact = FALSE)
    cat(sprintf("  %s vs %s: W = %.1f, p = %.4f\n",
                a, b, wt$statistic, wt$p.value))
  }
  
  list(chi_F = chi_F, p_val = p_val, mean_ranks = R_j)
}

# =============================================================================
# 8. SUMMARY TABLES AND VISUALIZATION
# =============================================================================

summarize_results <- function(results_matrix) {
  med_r2 <- apply(results_matrix, 2, median, na.rm = TRUE)
  mad_r2 <- apply(results_matrix, 2, compute_mad)
  
  summary_df <- data.frame(
    Model    = names(med_r2),
    Median_R2 = round(med_r2, 3),
    MAD       = round(mad_r2, 3)
  )
  summary_df <- summary_df[order(-summary_df$Median_R2), ]
  rownames(summary_df) <- NULL
  summary_df
}

plot_performance <- function(results_list, dataset_names) {
  # results_list: named list of results matrices, one per dataset
  
  plot_data <- do.call(rbind, lapply(seq_along(results_list), function(i) {
    mat <- results_list[[i]]
    data.frame(
      Model   = rep(colnames(mat), each = nrow(mat)),
      R2      = as.vector(mat),
      Dataset = dataset_names[i]
    )
  }))
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(Model, R2, median),
                                          y = R2, fill = Dataset)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Out-of-sample R² by model and dataset",
      x     = NULL,
      y     = expression(Out-of-sample~R^2)
    ) +
    ggplot2::theme_bw(base_size = 12)
}

# =============================================================================
# 9. EXAMPLE USAGE (Wolf telemetry dataset)
# =============================================================================

# Replace with actual data path
# wolf_data <- read.csv("data/wolf_telemetry.csv")
#
# wolf_results <- run_block_cv(
#   data         = wolf_data,
#   block_var    = "pack_id",
#   y_var        = "used",
#   habitat_vars = c("dist_human", "dist_edge", "slope",
#                    "burnt", "alpine", "shrub", "rock", "herbaceous"),
#   family       = "binomial",
#   n_basis      = 5,
#   n_trees_rf   = 500
# )
#
# summarize_results(wolf_results)
# run_friedman_test(wolf_results)

# =============================================================================
# Session info for reproducibility
# =============================================================================

sessionInfo()
