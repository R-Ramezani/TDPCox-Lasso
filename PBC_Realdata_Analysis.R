# ============================================================
# Real Data Analysis: PBC Dataset (Time-dependent Cox Model)
# Authors: Rohollah Ramezani, Mohammad Reza Rabiei, Mohammad Arashi
#
# Description:
# This script applies:
# - AACD (C++ Cox LASSO path)
# - Standard LASSO (glmnet)
# on the PBC (Primary Biliary Cirrhosis) dataset
#
# Outputs:
# - Coefficient estimates + standard errors
# - C-index
# - Time-dependent AUC (1, 3, 5 years)
# - Execution time
# - Results saved to Excel file
# ============================================================


############################################################
# 0) Global Options
############################################################
set.seed(123)
options(stringsAsFactors = FALSE)


############################################################
# 1) Load Required Packages
############################################################
pkgs <- c(
  "survival",
  "glmnet",
  "timeROC",
  "splines",
  "dplyr",
  "openxlsx"
)

for(p in pkgs){
  if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}


############################################################
# 2) Load C++ Solver
############################################################
# Make sure this file exists in the working directory
Rcpp::sourceCpp("fit_cox_Lasso.cpp")


############################################################
# 3) Prepare Time-dependent PBC Data
############################################################
data(pbc)
data(pbcseq)

# Baseline dataset
base <- subset(pbc, select = c(id:sex, stage))

# Create survival structure
pbc2 <- tmerge(base, base, id = id,
               death = event(time, status))

# Add time-dependent covariates
pbc2 <- tmerge(
  pbc2, pbcseq, id = id,
  ascites   = tdc(day, ascites),
  bili      = tdc(day, bili),
  albumin   = tdc(day, albumin),
  protime   = tdc(day, protime),
  alk.phos  = tdc(day, alk.phos),
  ast       = tdc(day, ast),
  chol      = tdc(day, chol),
  platelet  = tdc(day, platelet),
  edema     = tdc(day, edema),
  hepato    = tdc(day, hepato),
  spiders   = tdc(day, spiders)
)

# Remove missing values and invalid intervals
pbc2 <- na.omit(pbc2)
pbc2 <- pbc2[pbc2$tstart < pbc2$tstop, ]


############################################################
# 4) Construct Model Matrix
############################################################
X <- model.matrix(
  ~ log(bili) + log(protime) +
    log(alk.phos + 1) + log(ast + 1) +
    albumin + chol +
    platelet +
    stage + ascites + edema + hepato + spiders + sex,
  data = pbc2
)[,-1]

# Standardization
X <- scale(X)


############################################################
# 5) Survival Variables
############################################################
start <- pbc2$tstart
stop  <- pbc2$tstop
event <- as.numeric(pbc2$death == 2)


############################################################
# 6) AACD Model (C++ LASSO Path)
############################################################
time_aacd <- system.time({
  
  fpenalty <- 0
  
  score0 <- cox_score_at_zero_cpp(X, start, stop, event)
  
  lmax <- max(abs(score0))
  
  lseq <- exp(seq(
    log(lmax),
    log(lmax * 0.001),
    length.out = 100
  ))
  
  bpath <- fit_cox_path1_cpp(
    X, start, stop, event,
    lseq,
    fpenalty
  )
  
  # Model selection via EBIC/BIC
  bic_res <- fast_bic_selection(
    bpath, X, start, stop, event,
    nrow(X)
  )
  
  beta_aacd <- bic_res$best_beta
})

# Risk score
risk_aacd <- as.vector(X %*% beta_aacd)


############################################################
# 7) Standard LASSO (glmnet)
############################################################
y <- Surv(start, stop, event)

time_lasso <- system.time({
  
  cvfit <- cv.glmnet(
    X, y,
    family = "cox",
    alpha = 1,
    nfolds = 5,
    standardize = FALSE
  )
  
  beta_lasso <- as.vector(coef(cvfit, s = "lambda.min"))
})

risk_lasso <- as.vector(X %*% beta_lasso)


############################################################
# 8) Concordance Index (C-index)
############################################################
cindex_aacd  <- concordance(Surv(stop, event) ~ I(-risk_aacd))$concordance
cindex_lasso <- concordance(Surv(stop, event) ~ I(-risk_lasso))$concordance


############################################################
# 9) Time-dependent AUC
############################################################
times <- c(365, 1095, 1825)

roc_aacd <- timeROC(
  T = stop,
  delta = event,
  marker = risk_aacd,
  cause = 1,
  times = times
)

roc_lasso <- timeROC(
  T = stop,
  delta = event,
  marker = risk_lasso,
  cause = 1,
  times = times
)


############################################################
# 10) Coefficient Tables
############################################################
varnames <- colnames(X)

coef_aacd <- data.frame(
  Variable = varnames,
  Coefficient = beta_aacd
)

coef_lasso <- data.frame(
  Variable = varnames,
  Coefficient = beta_lasso
)


############################################################
# 11) Standard Errors via Cox Refitting
############################################################
selected_aacd  <- which(beta_aacd != 0)
selected_lasso <- which(beta_lasso != 0)

cox_aacd <- coxph(
  Surv(start, stop, event) ~ X[, selected_aacd]
)

cox_lasso <- coxph(
  Surv(start, stop, event) ~ X[, selected_lasso]
)

se_aacd  <- sqrt(diag(vcov(cox_aacd)))
se_lasso <- sqrt(diag(vcov(cox_lasso)))

coef_aacd$SE <- NA
coef_aacd$SE[selected_aacd] <- se_aacd

coef_lasso$SE <- NA
coef_lasso$SE[selected_lasso] <- se_lasso


############################################################
# 12) Final Performance Table
############################################################
perf <- data.frame(
  Method = c("AACD", "LASSO"),
  
  Runtime_seconds = c(
    time_aacd[3],
    time_lasso[3]
  ),
  
  Cindex = c(
    cindex_aacd,
    cindex_lasso
  ),
  
  AUC1 = c(
    roc_aacd$AUC[1],
    roc_lasso$AUC[1]
  ),
  
  AUC3 = c(
    roc_aacd$AUC[2],
    roc_lasso$AUC[2]
  ),
  
  AUC5 = c(
    roc_aacd$AUC[3],
    roc_lasso$AUC[3]
  )
)

print(perf)


############################################################
# 13) Save Results to Excel
############################################################
wb <- createWorkbook()

addWorksheet(wb, "Performance")
addWorksheet(wb, "AACD_coefficients")
addWorksheet(wb, "LASSO_coefficients")

writeData(wb, "Performance", perf)
writeData(wb, "AACD_coefficients", coef_aacd)
writeData(wb, "LASSO_coefficients", coef_lasso)

saveWorkbook(
  wb,
  "Results_PBC_DataSet.xlsx",
  overwrite = TRUE
)