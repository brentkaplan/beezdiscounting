# canonical dd stancode snapshot is stable

    Code
      cat(brms::make_stancode(spec$formula, data = d, family = spec$family))
    Output
      // generated with brms 2.22.0
      functions {
      }
      data {
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        int<lower=1> K_logk;  // number of population-level effects
        matrix[N, K_logk] X_logk;  // population-level design matrix
        // covariates for non-linear functions
        vector[N] C_1;
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_logk_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
      }
      parameters {
        vector[K_logk] b_logk;  // regression coefficients
        real<lower=0> phi;  // precision parameter
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        array[M_1] vector[N_1] z_1;  // standardized group-level effects
      }
      transformed parameters {
        vector[N_1] r_1_logk_1;  // actual group-level effects
        real lprior = 0;  // prior contributions to the log posterior
        r_1_logk_1 = (sd_1[1] * (z_1[1]));
        lprior += gamma_lpdf(phi | 0.01, 0.01);
        lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
          - 1 * student_t_lccdf(0 | 3, 0, 2.5);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] nlp_logk = rep_vector(0.0, N);
          // initialize non-linear predictor term
          vector[N] mu;
          nlp_logk += X_logk * b_logk;
          for (n in 1:N) {
            // add more terms to the linear predictor
            nlp_logk[n] += r_1_logk_1[J_1[n]] * Z_1_logk_1[n];
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            mu[n] = (1 / (10 ^ 6) + (1 - 2 / (10 ^ 6)) * (1 / (1 + exp(nlp_logk[n]) * C_1[n])));
          }
          target += beta_lpdf(Y | mu * phi, (1 - mu) * phi);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);
      }
      generated quantities {
      }

