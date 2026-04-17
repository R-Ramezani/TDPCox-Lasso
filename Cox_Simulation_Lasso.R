# ============================================================
# Time-dependent Cox Simulation + LASSO Path (C++) + LASSO (glmnet)
# Authors: Rohollah Ramezani, Mohammad Reza Rabiei, Mohammad Arashi
# Description:
# This script performs a Monte Carlo simulation study for 
# time-dependent Cox proportional hazards models using:
# 1) Custom C++ implementation of LASSO path (oplasso)
# 2) Standard LASSO via glmnet
#
# Outputs:
# - Performance metrics (FP, FN, Precision, Recall, F1)
# - Estimated coefficients (mean & sd)
# - Computational time
# - Results saved incrementally into an Excel file
# ============================================================


## ============================================================
## 0) Install & Load Required Packages
## ============================================================
required_pkgs <- c("Rcpp", "RcppArmadillo", "survival", "glmnet", "openxlsx")

for(pkg in required_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}


## ============================================================
## 1) Load C++ Implementation
## ============================================================
# Path to the C++ file containing core algorithms
file_path <- "fit_cox_lasso.cpp"

if(!file.exists(file_path)) stop("C++ file not found!")

Rcpp::sourceCpp(file_path)


## ============================================================
## 2) Utility Function: Append Results to Excel
## ============================================================
# This function appends a new row of results to an Excel sheet.
# If the file does not exist, it creates a new one.
save_append_excel_one_sheet <- function(df_row, file, sheet = "Results"){
  
  if(!file.exists(file)){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df_row)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    
  } else {
    wb <- openxlsx::loadWorkbook(file)
    
    if(!(sheet %in% names(wb))) {
      openxlsx::addWorksheet(wb, sheet)
    }
    
    old <- openxlsx::read.xlsx(file, sheet = sheet)
    
    # Align columns before binding
    all_cols <- union(names(old), names(df_row))
    old[setdiff(all_cols, names(old))] <- NA
    df_row[setdiff(all_cols, names(df_row))] <- NA
    
    out <- rbind(old[, all_cols], df_row[, all_cols])
    
    openxlsx::writeData(wb, sheet, out)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
}


# ============================================================
# 3) Simulation Grid Setup
# ============================================================
n_vec     <- c(200)
ratio_vec <- c(0.5,1,1.5)
rho_vec   <- c(0,0.25,0.5)
pc_vec    <- c(0.4,0.2,0)
fpenalty  <- 0   # LASSO penalty type

param_grid <- expand.grid(
  n = n_vec, 
  ratio = ratio_vec, 
  rho = rho_vec, 
  pc = pc_vec, 
  fpen = fpenalty
)

iterate <- 2
m <- 3
K_FOLDS <- 3

out_all <- "Results_Simulation.xlsx"
sheet_all <- "Results"


# ============================================================
# 4) Main Simulation Loop
# ============================================================
for(g in seq_len(nrow(param_grid))){
  
  n <- param_grid$n[g]
  rho <- param_grid$rho[g]
  pc <- param_grid$pc[g]
  ratio <- param_grid$ratio[g]
  fpen <- param_grid$fpen[g]
  
  p <- round(ratio * n)
  p1 <- 6              # Number of true non-zero coefficients
  p2 <- p - p1         # Number of zero coefficients
  
  cat(sprintf("\n--- Grid %d/%d: n=%d, p=%d ---\n", g, nrow(param_grid), n, p))
  
  # True coefficients
  beta_true <- c(1.2, 0.8, 0.5, -0.5, -0.8, -1.2, rep(0, p2))
  true_active <- which(beta_true != 0)
  
  # Storage variables
  beta_store_oplasso <- matrix(NA_real_, iterate, p)
  beta_store_lasso   <- matrix(NA_real_, iterate, p)
  
  FP_oplasso <- FN_oplasso <- rep(NA, iterate)
  FP_lasso   <- FN_lasso   <- rep(NA, iterate)
  
  time_oplasso_total <- time_lasso_total <- rep(NA, iterate)
  
  lambda_oplasso_vec     <- rep(NA, iterate)
  lambda_lasso_1se_vec   <- rep(NA, iterate)
  
  # Accuracy metrics
  PREC_oplasso <- REC_oplasso <- F1_oplasso <- rep(NA, iterate)
  PREC_lasso   <- REC_lasso   <- F1_lasso   <- rep(NA, iterate)
  
  
  for (iter in seq_len(iterate)) {
    
    print(iter)
    
    # --------------------------------------------------------
    # 1) Data Generation
    # --------------------------------------------------------
    sim <- generate_cox_simulation_data_cpp(n, m, p, rho, pc, beta_true, 0.1)
    
    X <- sim$X
    y_start <- sim$y_start
    y_stop  <- sim$y_stop
    y_event <- sim$y_event
    
    y_surv <- survival::Surv(y_start, y_stop, y_event)
    
    
    # --------------------------------------------------------
    # 2) oplasso (C++ LASSO Path + BIC Selection)
    # --------------------------------------------------------
    t_oplasso <- system.time({
      
      res_oplasso <- tryCatch({
        
        score0 <- cox_score_at_zero_cpp(X, y_start, y_stop, y_event)
        l_max  <- max(abs(score0), na.rm = TRUE)
        
        l_seq  <- exp(seq(log(l_max), log(l_max*0.05), length.out = 30))
        
        b_path <- fit_cox_path_cpp(X, y_start, y_stop, y_event, l_seq, fpen)
        
        bic_res <- fast_bic_selection(b_path, X, y_start, y_stop, y_event, n)
        
        metrics <- calculate_metrics_cpp(bic_res$best_beta, true_active)
        
        list(
          beta = bic_res$best_beta,
          lambda = l_seq[bic_res$best_idx],
          fp = metrics$fp,
          fn = metrics$fn,
          precision = metrics$precision,
          recall = metrics$recall,
          f1 = metrics$f1
        )
        
      }, error = function(e) NULL)
    })
    
    
    if(!is.null(res_oplasso)){
      beta_store_oplasso[iter,] <- res_oplasso$beta
      lambda_oplasso_vec[iter]  <- res_oplasso$lambda
      
      FP_oplasso[iter] <- res_oplasso$fp
      FN_oplasso[iter] <- res_oplasso$fn
      
      PREC_oplasso[iter] <- res_oplasso$precision
      REC_oplasso[iter]  <- res_oplasso$recall
      F1_oplasso[iter]   <- res_oplasso$f1
      
      time_oplasso_total[iter] <- t_oplasso[["elapsed"]]
    }
    
    
    # --------------------------------------------------------
    # 3) Standard LASSO (glmnet)
    # --------------------------------------------------------
    t_lasso <- system.time({
      
      res_lasso <- tryCatch({
        
        cv_fit <- glmnet::cv.glmnet(
          X, y_surv, 
          family = "cox",
          nfolds = K_FOLDS,
          standardize = FALSE
        )
        
        b_lasso <- as.numeric(coef(cv_fit, s = "lambda.1se"))
        
        m_lasso <- calculate_metrics_cpp(b_lasso, true_active)
        
        list(
          beta = b_lasso,
          l_1se = cv_fit$lambda.1se,
          fp = m_lasso$fp,
          fn = m_lasso$fn,
          precision = m_lasso$precision,
          recall = m_lasso$recall,
          f1 = m_lasso$f1
        )
        
      }, error = function(e) NULL)
    })
    
    
    if(!is.null(res_lasso)){
      beta_store_lasso[iter,] <- res_lasso$beta
      lambda_lasso_1se_vec[iter] <- res_lasso$l_1se
      
      FP_lasso[iter] <- res_lasso$fp
      FN_lasso[iter] <- res_lasso$fn
      
      PREC_lasso[iter] <- res_lasso$precision
      REC_lasso[iter]  <- res_lasso$recall
      F1_lasso[iter]   <- res_lasso$f1
      
      time_lasso_total[iter] <- t_lasso[["elapsed"]]
    }
  }
  
  
  # ============================================================
  # 5) Summary Function
  # ============================================================
  prepare_summary <- function(method, beta_mat, fp, fn, time, l_vec,
                              prec, rec, f1) {
    
    df_metrics <- data.frame(
      method = method, 
      n = n, 
      p1 = p1,   
      p2 = p2,
      p = p, 
      rho = rho, 
      pc = pc,
      
      lambda_median = median(l_vec, na.rm=TRUE),
      
      mean_FP = mean(fp, na.rm=TRUE), 
      sd_FP   = sd(fp, na.rm=TRUE),
      
      mean_FN = mean(fn, na.rm=TRUE), 
      sd_FN   = sd(fn, na.rm=TRUE),
      
      mean_precision = mean(prec, na.rm=TRUE),
      mean_recall    = mean(rec, na.rm=TRUE),
      mean_F1        = mean(f1, na.rm=TRUE),
      
      sd_precision   = sd(prec, na.rm=TRUE),
      sd_recall      = sd(rec, na.rm=TRUE),
      sd_F1          = sd(f1, na.rm=TRUE),
      
      median_time_solver_path = median(time, na.rm=TRUE),
      mean_time_solver_path   = mean(time, na.rm=TRUE),
      
      median_time = median(time, na.rm=TRUE),
      mean_time   = mean(time, na.rm=TRUE),
      sd_time     = sd(time, na.rm=TRUE)
    )
    
    # Coefficient summaries
    b_mean <- colMeans(beta_mat, na.rm=TRUE)
    b_sd   <- apply(beta_mat, 2, sd, na.rm=TRUE)
    
    combined_betas <- matrix(NA, nrow = 1, ncol = 2 * p)
    combined_betas[1, seq(1, 2 * p, by = 2)] <- b_mean
    combined_betas[1, seq(2, 2 * p, by = 2)] <- b_sd
    
    beta_names <- character(2 * p)
    beta_names[seq(1, 2 * p, by = 2)] <- paste0("beta", 1:p, "_mean")
    beta_names[seq(2, 2 * p, by = 2)] <- paste0("beta", 1:p, "_sd")
    
    beta_df <- as.data.frame(combined_betas)
    names(beta_df) <- beta_names
    
    return(cbind(df_metrics, beta_df))
  }
  
  
  # ------------------------------------------------------------
  # Save Results
  # ------------------------------------------------------------
  res_oplasso <- prepare_summary(
    "oplasso", beta_store_oplasso,
    FP_oplasso, FN_oplasso, time_oplasso_total,
    lambda_oplasso_vec,
    PREC_oplasso, REC_oplasso, F1_oplasso
  )
  
  res_lasso <- prepare_summary(
    "LASSO", beta_store_lasso,
    FP_lasso, FN_lasso, time_lasso_total,
    lambda_lasso_1se_vec,
    PREC_lasso, REC_lasso, F1_lasso
  )
  
  final_row <- rbind(res_oplasso, res_lasso)
  
  save_append_excel_one_sheet(final_row, out_all, sheet = sheet_all)
  
  rm(res_oplasso, res_lasso, final_row)
  gc()
}

cat("\n✅ Simulation completed successfully.\n")