////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 0, Group Level, Mixture
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept, and observation level location ID.
// No boat-level predictors.
// Based on M0, but uses a mixture distribution for the possibility of outliers.
// Time-stamp: <2017-04-27 15:45:31 (slane)>
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

data{
  // Data to be supplied to sampler
  /* Number of (fully observed) observations */
  int<lower=1> N;
  /* Number of (censored) observations */
  int<lower=1> nCens;
  /* Categorical predictors, entered as matrices of indicators */
  /* Location of measurement, hull as base case */
  int<lower=1> numLoc;
  matrix[N, numLoc - 1] locID;
  matrix[nCens, numLoc - 1] locIDCens;
  /* Boat random effect */
  int<lower=1> numBoat;
  int<lower=1,upper=numBoat> boatID[N];
  int<lower=1,upper=numBoat> boatIDCens[nCens];
  /* Observed data */
  real<lower=1.5> Y[N];
  /* Truncated data (brute force, all equal */
  real<upper=min(Y)> U;
}

parameters{
  // Parameters for the model
  /* Intercept */
  real mu1;
  real mu2;
  /* Betas for categorical indicators */
  vector[numLoc - 1] betaLoc1;
  vector[numLoc - 1] betaLoc2;
  /* Alphas for modelled random effect */
  vector[numBoat] alphaBoat1;
  vector[numBoat] alphaBoat2;
  /* Errors for categorical predictors */
  real<lower=0> sigma_alphaBoat1;
  real<lower=0> sigma_alphaBoat2;
  /* Error */
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  /* Mixing proportions */
  simplex[2] theta;
}

transformed parameters{
  // Make it easier for some sampling statements (not necessary)
  /* Regression for observed data */
  vector[N] muHat1;
  vector[N] muHat2;
  /* Regression for censored data */
  vector[nCens] muHatCens1;
  vector[nCens] muHatCens2;
  for(i in 1:N){
    muHat1[i] = mu1 + locID[i] * betaLoc1 + alphaBoat1[boatID[i]];
    muHat2[i] = mu2 + locID[i] * betaLoc2 + alphaBoat2[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens1[j] = mu1 + locIDCens[j] * betaLoc1 + alphaBoat1[boatIDCens[j]];
    muHatCens2[j] = mu2 + locIDCens[j] * betaLoc2 + alphaBoat2[boatIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Temp variable for mixing and log-proportions */
  vector[2] temp_p;
  vector[2] log_theta;
  log_theta = log(theta);
  /* Priors for categorical indicators */
  betaLoc1 ~ student_t(3, 0, 1);
  betaLoc2 ~ student_t(3, 0, 1);
  /* Priors for modelled random effect */
  mu1 ~ normal(0, 5);
  sigma_alphaBoat1 ~ cauchy(0, 2.5);
  alphaBoat1 ~ cauchy(0, sigma_alphaBoat1);
  mu2 ~ normal(0, 5);
  sigma_alphaBoat2 ~ cauchy(0, 2.5);
  alphaBoat2 ~ cauchy(0, sigma_alphaBoat2);
  /* Prior for observation (model) error */
  sigma1 ~ cauchy(0, 2.5);
  sigma2 ~ cauchy(0, 5);
  /* Observed log-likelihood */
  for(i in 1:N){
    temp_p[1] = log_theta[1] + lognormal_lpdf(Y[i] | muHat1[i], sigma1);
    temp_p[2] = log_theta[2] + lognormal_lpdf(Y[i] | muHat2[i], sigma2);
    target += log_sum_exp(temp_p);
  }
  /* Censored log-likelihood */
  for(i in 1:nCens){
    temp_p[1] = log_theta[1] + lognormal_lcdf(U | muHat1[i], sigma1);
    temp_p[2] = log_theta[2] + lognormal_lcdf(U | muHat2[i], sigma2);
    target += log_sum_exp(temp_p);
  }
}

generated quantities{
  /* Need the log likelihood for leave one out prediction errors */
  vector[N + nCens] log_lik;
  {
    real linPred1;
    real linPred2;
    vector[2] log_theta;
    log_theta = log(theta);
    for(i in 1:N){
      linPred1 = mu1 + locID[i] * betaLoc1 + alphaBoat1[boatID[i]];
      linPred2 = mu2 + locID[i] * betaLoc2 + alphaBoat2[boatID[i]];
      log_lik[i] = sum(log_theta) + lognormal_lpdf(Y[i] | linPred1, sigma1) +
	lognormal_lpdf(Y[i] | linPred2, sigma2);
    }
    for(j in 1:nCens){
      linPred1 = mu1 + locIDCens[j] * betaLoc1 + alphaBoat1[boatIDCens[j]];
      linPred2 = mu2 + locIDCens[j] * betaLoc2 + alphaBoat2[boatIDCens[j]];
      log_lik[N + j] = sum(log_theta) + lognormal_lcdf(U | linPred1, sigma1) +
	lognormal_lcdf(U | linPred2, sigma2);
    }
  }
}
