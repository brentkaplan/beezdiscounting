/// @file ChoiceDiscounting.h
/// Trial-level SS-vs-LL choice family (Family 2) via TMB (beezdiscounting).
/// =============================================================================
/// Binomial/logit GLMM of choice (1 = chose larger-later). mode switch:
///   mode 0 = STRUCTURAL: eta = beta0*has_intercept
///                              + gamma * ((ll/ss) * D(k, delay) - 1)
///            with k = exp(X*beta_k + sigma_u * u(subj)),  gamma = exp(log_gamma).
///            D: eqn_type 0 = mazur 1/(1+k*x); 1 = exponential exp(-k*x).
///            Scale-invariant relative comparison; beta0=0 => P=0.5 at indiff.
///            RE: non-centered random intercept on log k, u(subj) ~ N(0,1).
///   mode 1 = DESCRIPTIVE (Young 2018): eta = Z*theta + Zre * (L * b(subj)),
///            no discount function. Correlated per-subject random SLOPES on the
///            Young predictors (log(ll/ss), log(delay+1)). Sigma = diag(sd) R
///            diag(sd), rho = tanh(cor_re), L = chol(Sigma); b(subj) ~ N(0, I_q).
///            n_re = q (2 = random slopes; 0 = pooled fixed-effect logistic).
/// =============================================================================
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type ChoiceDiscounting(objective_function<Type>* obj) {
  DATA_INTEGER(mode);            // 0 structural, 1 descriptive
  DATA_INTEGER(eqn_type);        // 0 mazur, 1 exponential (structural)
  DATA_INTEGER(has_intercept);   // structural beta0 on/off
  DATA_VECTOR(choice);           // 0/1 (1 = LL), length n_obs
  DATA_VECTOR(ss_amount);
  DATA_VECTOR(ll_amount);
  DATA_VECTOR(delay);
  DATA_IVECTOR(subject_id);      // 0-indexed
  DATA_MATRIX(X);                // n_obs x ncol(X) design for log k (structural)
  DATA_MATRIX(Z);                // n_obs x p_z fixed design (descriptive)
  DATA_MATRIX(Zre);              // n_obs x q random-slope design (descriptive)
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_subjects);
  DATA_INTEGER(n_re);            // q (descriptive RE dim; 0 = pooled)

  PARAMETER_VECTOR(beta_k);      // fixed effects on log k (structural)
  PARAMETER(log_sigma_u);        // log SD of the random intercept on log k
  PARAMETER(log_gamma);          // log choice sensitivity (structural)
  PARAMETER(beta0);              // choice-bias intercept (structural)
  PARAMETER_MATRIX(u);           // n_subjects x 1 standardized RE (structural)
  PARAMETER_VECTOR(theta);       // descriptive fixed sensitivities
  PARAMETER_VECTOR(log_sd_re);   // descriptive log SDs (len q)
  PARAMETER_VECTOR(cor_re);      // descriptive unconstrained correlation(s)
  PARAMETER_MATRIX(b);           // n_subjects x q standardized RE (descriptive)

  Type nll = Type(0.0);

  if (mode == 0) {
    Type sigma_u = exp(log_sigma_u);
    Type gamma = exp(log_gamma);
    for (int i = 0; i < n_subjects; i++) {
      nll -= dnorm(u(i, 0), Type(0.0), Type(1.0), true);
    }
    for (int i = 0; i < n_obs; i++) {
      int subj = subject_id(i);
      vector<Type> x_i = X.row(i);
      Type log_k_i = (x_i * beta_k).sum() + sigma_u * u(subj, 0);
      Type k_i = exp(log_k_i);
      Type D;
      if (eqn_type == 0) {
        D = Type(1.0) / (Type(1.0) + k_i * delay(i));   // mazur
      } else if (eqn_type == 1) {
        D = exp(-k_i * delay(i));                        // exponential
      } else {
        error("ChoiceDiscounting: unknown eqn_type (expected 0=mazur, 1=exponential).");
      }
      Type eta = gamma * ((ll_amount(i) / ss_amount(i)) * D - Type(1.0));
      if (has_intercept == 1) eta += beta0;
      nll -= dbinom_robust(choice(i), Type(1.0), eta, true);
    }
    ADREPORT(beta_k);
    ADREPORT(gamma);
    ADREPORT(sigma_u);
    return nll;
  } else if (mode == 1) {
    // Build the q x q covariance + lower Cholesky from log-SDs + correlation.
    // q == 2 this slice (Young's two predictors); general-q deferred.
    matrix<Type> L(n_re, n_re);
    matrix<Type> re(n_subjects, (n_re > 0 ? n_re : 1));
    if (n_re > 0) {
      vector<Type> sd_re = exp(log_sd_re);
      matrix<Type> Sigma(n_re, n_re);
      if (n_re == 2) {
        Type rho = tanh(cor_re(0));
        Sigma(0, 0) = sd_re(0) * sd_re(0);
        Sigma(1, 1) = sd_re(1) * sd_re(1);
        Sigma(0, 1) = rho * sd_re(0) * sd_re(1);
        Sigma(1, 0) = Sigma(0, 1);
        ADREPORT(rho);
      } else if (n_re == 1) {
        Sigma(0, 0) = sd_re(0) * sd_re(0);
      } else {
        error("ChoiceDiscounting: descriptive mode supports n_re in {0,1,2} this slice.");
      }
      L = Sigma.llt().matrixL();
      for (int i = 0; i < n_subjects; i++) {
        vector<Type> b_i = b.row(i);
        vector<Type> re_i = L * b_i;
        for (int j = 0; j < n_re; j++) re(i, j) = re_i(j);
        for (int j = 0; j < n_re; j++) {
          nll -= dnorm(b(i, j), Type(0.0), Type(1.0), true);  // N(0,I) prior
        }
      }
      ADREPORT(sd_re);
      ADREPORT(Sigma);
    }
    for (int i = 0; i < n_obs; i++) {
      int subj = subject_id(i);
      vector<Type> z_i = Z.row(i);
      Type eta = (z_i * theta).sum();
      if (n_re > 0) {
        vector<Type> zre_i = Zre.row(i);
        for (int j = 0; j < n_re; j++) eta += zre_i(j) * re(subj, j);
      }
      nll -= dbinom_robust(choice(i), Type(1.0), eta, true);
    }
    return nll;
  } else {
    error("ChoiceDiscounting: unknown mode (expected 0=structural, 1=descriptive).");
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
