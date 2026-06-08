/// @file ChoiceDiscounting.h
/// Trial-level SS-vs-LL choice family (Family 2) via TMB (beezdiscounting).
/// =============================================================================
/// Binomial/logit GLMM of choice (1 = chose larger-later). mode switch:
///   mode 0 = STRUCTURAL: eta = beta0*has_intercept
///                              + gamma * ((ll/ss) * D(k, delay) - 1)
///            with k = exp(X*beta_k + sigma_u * u(subj)),  gamma = exp(log_gamma).
///            D: eqn_type 0 = mazur 1/(1+k*x); 1 = exponential exp(-k*x).
///            Scale-invariant relative comparison; beta0=0 => P=0.5 at indiff.
///   mode 1 = DESCRIPTIVE: not implemented in this slice (errors).
/// Non-centered random intercept on log k: u(subj) ~ N(0,1).
/// =============================================================================
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type ChoiceDiscounting(objective_function<Type>* obj) {
  DATA_INTEGER(mode);            // 0 structural, 1 descriptive (not impl.)
  DATA_INTEGER(eqn_type);        // 0 mazur, 1 exponential (structural)
  DATA_INTEGER(has_intercept);   // structural beta0 on/off
  DATA_VECTOR(choice);           // 0/1 (1 = LL), length n_obs
  DATA_VECTOR(ss_amount);
  DATA_VECTOR(ll_amount);
  DATA_VECTOR(delay);
  DATA_IVECTOR(subject_id);      // 0-indexed
  DATA_MATRIX(X);                // n_obs x ncol(X) design for log k
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_subjects);

  PARAMETER_VECTOR(beta_k);      // fixed effects on log k
  PARAMETER(log_sigma_u);        // log SD of the random intercept on log k
  PARAMETER(log_gamma);          // log choice sensitivity (inverse temperature)
  PARAMETER(beta0);              // choice-bias intercept (used iff has_intercept)
  PARAMETER_MATRIX(u);           // n_subjects x 1 standardized random intercepts

  if (mode != 0) {
    error("ChoiceDiscounting: only mode 0 (structural) is implemented.");
  }

  Type sigma_u = exp(log_sigma_u);
  Type gamma = exp(log_gamma);

  Type nll = Type(0.0);
  for (int i = 0; i < n_subjects; i++) {
    nll -= dnorm(u(i, 0), Type(0.0), Type(1.0), true);
  }

  for (int i = 0; i < n_obs; i++) {
    int subj = subject_id(i);
    vector<Type> x_i = X.row(i);
    Type log_k_i = (x_i * beta_k).sum() + sigma_u * u(subj, 0);
    Type k_i = exp(log_k_i);

    // Discount value (SS immediate => V_SS = ss_amount).
    Type D;
    if (eqn_type == 0) {
      D = Type(1.0) / (Type(1.0) + k_i * delay(i));   // mazur
    } else if (eqn_type == 1) {
      D = exp(-k_i * delay(i));                        // exponential
    } else {
      error("ChoiceDiscounting: unknown eqn_type (expected 0=mazur, 1=exponential).");
    }
    // Scale-invariant relative comparison.
    // ss_amount > 0 is guaranteed by .dd_validate_choice (no zero-division guard here).
    Type eta = gamma * ((ll_amount(i) / ss_amount(i)) * D - Type(1.0));
    if (has_intercept == 1) {
      eta += beta0;
    }

    // Binomial/logit log-likelihood, numerically robust (logit-parameterized).
    nll -= dbinom_robust(choice(i), Type(1.0), eta, true);
  }

  ADREPORT(beta_k);
  ADREPORT(gamma);
  ADREPORT(sigma_u);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
