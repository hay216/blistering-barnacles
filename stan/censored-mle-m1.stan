////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 1
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept. Edited to include predictions. Removed some
// predictors which don't seem to do anything.
// Time-stamp: <2017-04-21 13:39:29 (slane)>
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

data{
  // Data to be supplied to sampler
  /* Number of (fully observed) observations */
  int<lower=1> N;
  /* Number of (censored) observations */
  int<lower=1> nCens;
  /* Numeric/ordinal predictors */
  real days2[N];
  real days2Cens[nCens];
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
  real mu;
  /* Betas for continuous */
  real betaDays2;
  /* Betas for categorical indicators */
  vector[numLoc - 1] betaLoc;
  /* Alphas for modelled random effect */
  vector[numBoat] alphaBoat;
  /* Errors for categorical predictors */
  real<lower=0> sigma_alphaBoat;
  /* Censored data */
  vector<lower=0,upper=U>[nCens] yCens;
  /* Error */
  real<lower=0> sigma;
}

transformed parameters{
  // Make it easier for some sampling statements (not necessary)
    /* Regression for observed data */
  vector[N] muHat;
  /* Regression for censored data */
  vector[nCens] muHatCens;
  for(i in 1:N){
    muHat[i] = mu + betaDays2 * days2[i] + locID[i] * betaLoc + alphaBoat[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = mu + betaDays2 * days2Cens[j] + locIDCens[j] * betaLoc + alphaBoat[boatIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Priors for regression parameters */
  mu ~ normal(0, 5);
  betaDays2 ~ student_t(3, 0, 1);
  /* Priors for categorical indicators */
  betaLoc ~ student_t(3, 0, 1);
  /* Priors for modelled random effect */
  sigma_alphaBoat ~ cauchy(0, 2.5);
  alphaBoat ~ cauchy(0, sigma_alphaBoat);
  /* Prior for observation (model) error */
  sigma ~ cauchy(0, 2.5);
  /* Observed log-likelihood */
  Y ~ lognormal(muHat, sigma);
  /* Censored log-likelihood */
  yCens ~ lognormal(muHatCens, sigma);
}

generated quantities{
  /* Need the log likelihood for leave one out prediction errors */
  vector[N + nCens] log_lik;
  {
    real linPred;
    for(i in 1:N){
      linPred = mu + betaDays2 * days2[i] + locID[i] * betaLoc + alphaBoat[boatID[i]];
      log_lik[i] = lognormal_lpdf(Y[i] | linPred, sigma);
    }
    for(j in 1:nCens){
      linPred = mu + betaDays2 * days2Cens[j] + locIDCens[j] * betaLoc + alphaBoat[boatIDCens[j]];
      log_lik[N + j] = lognormal_lpdf(yCens[j] | linPred, sigma);
    }
  }
}
