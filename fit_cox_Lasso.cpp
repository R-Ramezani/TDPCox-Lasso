// ============================================================
// Time-dependent Cox Model - LASSO Path (C++ Implementation)
// Authors:  Rohollah Ramezani, Mohammad Reza Rabiei, Mohammad Arashi
//
// Description:
// This file implements core algorithms for:
// - Cox partial likelihood derivatives
// - LASSO path estimation via coordinate descent
// - BIC/EBIC model selection
// - Simulation data generation
// - Performance metrics calculation
//
// Used via Rcpp in simulation.R
// ============================================================

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;


// ============================================================
// Soft Thresholding Operator (for LASSO)
// ============================================================
inline double soft_threshold_cpp(double z, double gamma) {
  if (z > gamma) return z - gamma;
  else if (z < -gamma) return z + gamma;
  else return 0.0;
}


// ============================================================
// Internal Function: Cox Score & Information (Fast Computation)
// ============================================================
void get_cox_derivs_fast_internal(
    const arma::vec& beta, const arma::mat& X, const arma::vec& start,
    const arma::vec& stop, const arma::vec& event, const arma::vec& uniq_times,
    arma::vec& U, arma::vec& I, double eps = 1e-9) {
  
  U.zeros();
  I.zeros();
  
  arma::vec eta = X * beta;
  double max_eta = eta.max();
  arma::vec exp_eta = exp(eta - max_eta);
  
  for (uword i = 0; i < uniq_times.n_elem; i++) {
    double t = uniq_times[i];
    
    arma::uvec risk_mask = find((start < t) && (stop >= t));
    arma::uvec ev_mask   = find((stop == t) && (event == 1));
    
    if (risk_mask.n_elem == 0 || ev_mask.n_elem == 0) continue;
    
    arma::vec w = exp_eta.elem(risk_mask);
    double s0 = sum(w);
    if (s0 <= 1e-12) continue;
    
    arma::mat X_risk = X.rows(risk_mask);
    arma::vec s1 = X_risk.t() * w;
    
    arma::mat Xw = X_risk.each_col() % sqrt(w);
    arma::vec s2 = sum(square(Xw), 0).t();
    
    arma::vec E = s1 / s0;
    arma::vec V = s2 / s0 - square(E);
    
    V.clean(0.0);
    
    double m_events = (double)ev_mask.n_elem;
    arma::vec sum_X_event = sum(X.rows(ev_mask), 0).t();
    
    U += sum_X_event - m_events * E;
    I += m_events * V;
  }
  
  // Avoid division by zero
  I.replace(0.0, eps); 
}


// ============================================================
// Score Function at Beta = 0
// ============================================================
// [[Rcpp::export]]
arma::vec cox_score_at_zero_cpp(
    const arma::mat& X, 
    const arma::vec& start, 
    const arma::vec& stop, 
    const arma::vec& event) {
  
  arma::mat Xs = X;
  arma::vec starts = start;
  arma::vec stops  = stop;
  arma::vec events = event;
  
  // Sort by stop time
  arma::uvec ord = arma::sort_index(stops);
  stops  = stops.elem(ord);
  starts = starts.elem(ord);
  events = events.elem(ord);
  Xs     = Xs.rows(ord);
  
  arma::uvec event_idx = find(events == 1);
  arma::vec uniq_times = unique(stops.elem(event_idx));
  
  arma::vec score = zeros<vec>(Xs.n_cols);
  
  for (uword i = 0; i < uniq_times.n_elem; i++) {
    double t = uniq_times[i];
    
    arma::uvec risk_mask = find((starts < t) && (stops >= t));
    arma::uvec ev_mask   = find((stops == t) && (events == 1));
    
    if (risk_mask.n_elem == 0 || ev_mask.n_elem == 0) continue;
    
    double s0 = (double)risk_mask.n_elem;
    arma::vec s1 = sum(Xs.rows(risk_mask), 0).t();
    
    arma::vec E = s1 / s0;
    
    score += sum(Xs.rows(ev_mask), 0).t() - 
             (double)ev_mask.n_elem * E;
  }
  
  return score;
}


// ============================================================
// Cox LASSO Path via Coordinate Descent
// ============================================================
// [[Rcpp::export]]
arma::mat fit_cox_path1_cpp(
    const arma::mat& X, const arma::vec& start, const arma::vec& stop,
    const arma::vec& event, const arma::vec& lambda_seq, int fpenalty,
    double a = 3.7, double kappa = 1.0, int max_lla_iter = 5,
    double tol_lla = 1e-4, int max_cd_iter = 100, double tol_cd = 1e-5) {

  int p = X.n_cols;
  int n_lambda = lambda_seq.n_elem;

  // Sort data by stop time
  arma::uvec ord = arma::sort_index(stop);
  arma::mat  Xs  = X.rows(ord);
  arma::vec  sta = start.elem(ord);
  arma::vec  sto = stop.elem(ord);
  arma::vec  ev  = event.elem(ord);

  arma::uvec event_idx = arma::find(ev == 1);
  arma::vec  uniq_times = arma::unique(sto.elem(event_idx));

  arma::mat beta_path(p, n_lambda, arma::fill::zeros);
  arma::vec beta_current(p, arma::fill::zeros);

  for (int i = 0; i < n_lambda; i++) {
    
    double lambda = lambda_seq[i];

    for (int lla_it = 0; lla_it < max_lla_iter; lla_it++) {
      
      arma::vec beta_old_lla = beta_current;

      // Penalty weights (LASSO)
      arma::vec weights;
      if (fpenalty == 0)
        weights = arma::vec(p, arma::fill::ones) * lambda;

      arma::vec eta = Xs * beta_current;

      for (int cd_it = 0; cd_it < max_cd_iter; cd_it++) {
        
        arma::vec beta_cd_old = beta_current;

        for (int j = 0; j < p; j++) {
          
          double Uj = 0.0, Ij = 0.0;
          
          for (arma::uword ti = 0; ti < uniq_times.n_elem; ti++) {
            
            double t = uniq_times[ti];
            
            arma::uvec risk = arma::find((sta < t) % (sto >= t));
            if (risk.n_elem == 0) continue;
            
            arma::vec exp_eta_risk = arma::exp(eta.elem(risk));
            double S0 = arma::sum(exp_eta_risk);
            if (S0 < 1e-10) continue;
            
            arma::vec xj_risk = Xs(risk, arma::uvec({(arma::uword)j}));
            
            double S1j = arma::dot(xj_risk, exp_eta_risk) / S0;
            double S2j = arma::dot(xj_risk % xj_risk, exp_eta_risk) / S0;
            
            arma::uvec ev_t = arma::find((sto == t) % (ev == 1));
            int d_t = ev_t.n_elem;
            
            arma::vec xj_evt = Xs(ev_t, arma::uvec({(arma::uword)j}));
            double sum_xj_evt = arma::sum(xj_evt);
            
            Uj += sum_xj_evt - d_t * S1j;
            Ij += d_t * (S2j - S1j * S1j);
          }
          
          double den = std::max(Ij, 1e-7);
          double zj  = beta_current[j] + Uj / den;
          
          double beta_new = soft_threshold_cpp(zj, weights[j] / den);
          
          double delta = beta_new - beta_current[j];
          
          if (std::abs(delta) > 1e-10) {
            eta += delta * Xs.col(j);
          }
          
          beta_current[j] = beta_new;
        }

        double change = arma::max(arma::abs(beta_current - beta_cd_old)) /
                        (arma::max(arma::abs(beta_cd_old)) + 1e-5);
        
        if (change < tol_cd) break;
      }

      double lla_change = arma::max(arma::abs(beta_current - beta_old_lla)) /
                          (arma::max(arma::abs(beta_old_lla)) + 1e-5);
      
      if (lla_change < tol_lla) break;
    }

    beta_path.col(i) = beta_current;
  }

  return beta_path;
}


// ============================================================
// Partial Log-Likelihood
// ============================================================
// [[Rcpp::export]]
double partial_loglik_cpp(
    const arma::vec& beta, const arma::mat& X, 
    const arma::vec& start, const arma::vec& stop, 
    const arma::vec& event) {
  
  if (arma::sum(event) == 0) return 0.0;

  arma::vec eta = X * beta;
  arma::vec exp_eta = arma::exp(eta);

  arma::uvec event_indices = arma::find(event == 1);
  arma::vec stop_event = stop.elem(event_indices);
  arma::vec uniq_times = arma::unique(stop_event);
  
  double loglik = 0.0;
  
  for (uword i = 0; i < uniq_times.n_elem; i++) {
    
    double t = uniq_times[i];
    
    arma::uvec in_risk = arma::find((start < t) && (stop >= t));
    arma::uvec ev_at_t = arma::find((stop == t) && (event == 1));
    
    if (in_risk.n_elem > 0 && ev_at_t.n_elem > 0) {
      
      double s0 = arma::sum(exp_eta.elem(in_risk));
      
      if (s0 > 1e-12) {
        double sum_eta_event = arma::sum(eta.elem(ev_at_t));
        
        loglik += sum_eta_event - 
                  (double)ev_at_t.n_elem * std::log(s0);
      }
    }
  }
  
  return loglik;
}


// ============================================================
// Simulation Data Generator (Time-dependent Cox)
// ============================================================
// [[Rcpp::export]]
Rcpp::List generate_cox_simulation_data_cpp(
    int n, int m, int p, double rho, 
    double target_censoring_rate, 
    arma::vec beta_true, 
    double baseline_hazard = 0.1) {

  mat Sigma(p, p);

  for(int i = 0; i < p; ++i) {
    for(int j = 0; j < p; ++j) {
      Sigma(i, j) = std::pow(rho, std::abs(i - j));
    }
  }

  mat R = chol(Sigma);
  mat X_subj = randn<mat>(n, p) * R;
  
  vec linpred = X_subj * beta_true;
  vec T_event = -log(randu<vec>(n)) / (baseline_hazard * exp(linpred));
  
  vec C(n);

  if(target_censoring_rate <= 0) {
    C.fill(T_event.max() * 2.0);
  } else {
    double avg_hazard = mean(baseline_hazard * exp(linpred));
    double censor_lambda =
      (target_censoring_rate / (1.0 - target_censoring_rate)) * avg_hazard;
    
    C = -log(randu<vec>(n)) / censor_lambda;
  }

  // Use Rcpp::pmin for numerical stability
  Rcpp::NumericVector T_vec = Rcpp::wrap(T_event);
  Rcpp::NumericVector C_vec = Rcpp::wrap(C);
  Rcpp::NumericVector time_rcpp = Rcpp::pmin(T_vec, C_vec);

  arma::vec time = Rcpp::as<arma::vec>(time_rcpp);
  vec event = conv_to<vec>::from(T_event <= C);

  int Nlong = n * m;

  mat X_long(Nlong, p);
  vec y_start(Nlong), y_stop(Nlong), y_event(Nlong), id(Nlong);

  for(int i = 0; i < n; ++i) {
    
    double ti = std::max(time[i], 1e-10);
    
    for(int j = 0; j < m; ++j) {
      
      int idx = i * m + j;
      
      id[idx] = i + 1;
      y_start[idx] = (ti / m) * j;
      y_stop[idx]  = (ti / m) * (j + 1);
      
      y_event[idx] = (j == (m - 1) && event[i] == 1) ? 1.0 : 0.0;
      
      X_long.row(idx) = X_subj.row(i);
    }
  }

  return Rcpp::List::create(
    Rcpp::_["X"] = X_long,
    Rcpp::_["y_start"] = y_start,
    Rcpp::_["y_stop"] = y_stop,
    Rcpp::_["y_event"] = y_event
  );
}


// ============================================================
// Fast BIC / EBIC Model Selection
// ============================================================
// [[Rcpp::export]]
Rcpp::List fast_bic_selection(
    const arma::mat& beta_path, const arma::mat& X, 
    const arma::vec& start, const arma::vec& stop, 
    const arma::vec& event, double n_samples,
    double gamma = 0.5) {

  int n_lambda = beta_path.n_cols;
  int p = X.n_cols;

  vec bic_vec(n_lambda);
  double d = arma::sum(event);

  for(int j = 0; j < n_lambda; ++j) {
    
    vec current_beta = beta_path.col(j);
    
    double loglik = partial_loglik_cpp(current_beta, X, start, stop, event);
    
    double df = arma::sum(
      arma::abs(current_beta) / (arma::abs(current_beta) + 1e-2)
    );

    bic_vec[j] = -2.0 * loglik +
                 std::log(d) * df +
                 2.0 * gamma * std::log(p) * df;
  }

  uword best_idx = bic_vec.index_min();

  return Rcpp::List::create(
    Rcpp::_["best_beta"] = beta_path.col(best_idx),
    Rcpp::_["best_idx"]  = best_idx + 1,
    Rcpp::_["bic_vec"]   = bic_vec
  );
}


// ============================================================
// Performance Metrics (FP, FN, Precision, Recall, F1)
// ============================================================
// [[Rcpp::export]]
Rcpp::List calculate_metrics_cpp(
    arma::vec beta_est, arma::uvec true_active_idx) {

  double thr = 1e-3;

  arma::uvec est_active = arma::find(arma::abs(beta_est) > thr);

  int fp = 0;
  for(uword i=0; i<est_active.n_elem; ++i) {
    bool found = false;
    for(uword j=0; j<true_active_idx.n_elem; ++j) {
      if(est_active[i] == (true_active_idx[j]-1)) {
        found = true; break;
      }
    }
    if(!found) fp++;
  }

  int fn = 0;
  for(uword i=0; i<true_active_idx.n_elem; ++i) {
    bool found = false;
    for(uword j=0; j<est_active.n_elem; ++j) {
      if((true_active_idx[i]-1) == est_active[j]) {
        found = true; break;
      }
    }
    if(!found) fn++;
  }

  int tp = true_active_idx.n_elem - fn;

  double precision = (tp + fp) > 0 ? (double)tp / (tp + fp) : 0.0;
  double recall    = (tp + fn) > 0 ? (double)tp / (tp + fn) : 0.0;

  double f1 = (precision + recall) > 0 ?
    2.0 * precision * recall / (precision + recall) : 0.0;

  return Rcpp::List::create(
    Rcpp::_["fp"] = fp,
    Rcpp::_["fn"] = fn,
    Rcpp::_["tp"] = tp,
    Rcpp::_["precision"] = precision,
    Rcpp::_["recall"] = recall,
    Rcpp::_["f1"] = f1
  );
}