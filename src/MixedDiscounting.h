/// @file MixedDiscounting.h
/// IP-family Mixed-Effects Discounting via TMB (beezdiscounting 0.4.0)
/// =============================================================================
///
/// MODEL STRUCTURE:
/// One-parameter discounting (Mazur hyperbolic or exponential) for indifference
/// proportions y in [0,1], with a random intercept on log k and between-subject
/// fixed effects via design matrix X. Two observation families dispatch at
/// runtime: scale-location-truncated beta (SLT-beta) and Gaussian.
///
///   eqn_type: 0 = mazur        mu = 1 / (1 + k*x)
///             1 = exponential  mu = exp(-k*x)
///   family:   0 = sltb         y ~ SLTBeta(mu, phi)      aux = phi
///             1 = gaussian     y ~ N(mu, sigma_e)        aux = sigma_e
///
/// LINEAR PREDICTOR (log k scale, single random intercept, non-centered):
///   log_k_i = (X.row(i) * beta_k).sum() + sigma_u * u(subj, 0)
///   k_i     = exp(log_k_i)
///   sigma_u = exp(log_sigma_u),  u(i,0) ~ N(0,1)
///
/// AUXILIARY: log_aux is a single generic scalar.
///   family == 0 (sltb):     phi     = exp(log_aux)
///   family == 1 (gaussian): sigma_e = exp(log_aux)
///
/// SLT-BETA DENSITY (constants s = 1.0000001, l = 1e-8; verified in
/// dev/sltb-verification/verify_sltb.R):
///   a = mu*phi, b = (1-mu)*phi
///   log f = lgamma(a+b) - lgamma(a) - lgamma(b)
///         + (a-1)*log(y/s + l) + (b-1)*log(1 - (y/s + l))
///         - log(s)
///         - log( pbeta(1/s + l, a, b) - pbeta(l, a, b) )
/// The truncation normalizer Z = pbeta(1/s+l,a,b) - pbeta(l,a,b) is load-bearing.
///
/// MU is guarded to [1e-6, 1-1e-6] via CppAD::CondExp (exponential underflows to
/// 0 at long delays). subject_id is 0-indexed from R.
/// =============================================================================

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
Type MixedDiscounting(objective_function<Type>* obj) {
  // =========================================================================
  // DATA
  // =========================================================================
  DATA_VECTOR(y);               // indifference proportions in [0,1], length n_obs
  DATA_VECTOR(x);               // delays, length n_obs
  DATA_IVECTOR(subject_id);     // 0-indexed subject index per observation
  DATA_MATRIX(X);               // n_obs x ncol(X) fixed-effect design for log k
  DATA_INTEGER(eqn_type);       // 0 = mazur, 1 = exponential
  DATA_INTEGER(family);         // 0 = sltb, 1 = gaussian
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_subjects);

  // =========================================================================
  // PARAMETERS
  // =========================================================================
  PARAMETER_VECTOR(beta_k);     // length ncol(X): fixed effects on log k
  PARAMETER(log_sigma_u);       // log SD of random intercept on log k
  PARAMETER(log_aux);           // log phi (sltb) OR log sigma_e (gaussian)
  PARAMETER_MATRIX(u);          // n_subjects x 1 standardized random intercepts

  // =========================================================================
  // TRANSFORM SCALARS
  // =========================================================================
  Type sigma_u = exp(log_sigma_u);
  Type aux = exp(log_aux);      // phi (family 0) or sigma_e (family 1)

  // SLT-beta fixed constants (verified reference values).
  Type s = Type(1.0000001);
  Type l = Type(1e-8);

  // mu guard bounds.
  Type mu_lo = Type(1e-6);
  Type mu_hi = Type(1.0) - Type(1e-6);

  Type log2pi = log(Type(2.0) * M_PI);

  // =========================================================================
  // NEGATIVE LOG-LIKELIHOOD
  // =========================================================================
  Type nll = Type(0.0);

  // Standard-normal prior on the standardized random intercepts (non-centered;
  // sigma_u enters multiplicatively in the linear predictor below).
  for (int i = 0; i < n_subjects; i++) {
    nll -= dnorm(u(i, 0), Type(0.0), Type(1.0), true);
  }

  // Data likelihood.
  for (int i = 0; i < n_obs; i++) {
    int subj = subject_id(i);

    // Linear predictor on log k: fixed effects + random intercept.
    vector<Type> x_i = X.row(i);
    Type log_k_i = (x_i * beta_k).sum() + sigma_u * u(subj, 0);
    Type k_i = exp(log_k_i);

    // Mean via the discounting function (identity link).
    Type mu_raw;
    if (eqn_type == 0) {
      // Mazur hyperbolic.
      mu_raw = Type(1.0) / (Type(1.0) + k_i * x(i));
    } else {
      // Exponential.
      mu_raw = exp(-k_i * x(i));
    }

    // Guard mu to [mu_lo, mu_hi] (exponential underflow at long delays).
    // Never branch on a Type-valued condition with `if`; use CppAD::CondExp.
    Type mu = CppAD::CondExpLt(mu_raw, mu_lo, mu_lo, mu_raw);
    mu = CppAD::CondExpGt(mu, mu_hi, mu_hi, mu);

    if (family == 0) {
      // ------------------------------------------------------------------
      // SLT-beta density. aux = phi.
      // ------------------------------------------------------------------
      Type phi = aux;
      Type a = mu * phi;
      Type b = (Type(1.0) - mu) * phi;

      Type yt = y(i) / s + l;                  // scaled-location transform
      Type Z = pbeta(Type(1.0) / s + l, a, b)  // truncation normalizer
               - pbeta(l, a, b);

      Type logf = lgamma(a + b) - lgamma(a) - lgamma(b)
                  + (a - Type(1.0)) * log(yt)
                  + (b - Type(1.0)) * log(Type(1.0) - yt)
                  - log(s)
                  - log(Z);
      nll -= logf;

    } else {
      // ------------------------------------------------------------------
      // Gaussian density. aux = sigma_e.
      // ------------------------------------------------------------------
      Type sigma_e = aux;
      Type resid = (y(i) - mu) / sigma_e;
      nll -= -log_aux - Type(0.5) * log2pi - Type(0.5) * resid * resid;
    }
  }

  // =========================================================================
  // ADREPORT
  // =========================================================================
  ADREPORT(beta_k);     // fixed effects on log k
  ADREPORT(sigma_u);    // SD of random intercept on log k
  ADREPORT(aux);        // phi (family 0) or sigma_e (family 1)

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
