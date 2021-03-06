////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE 2 (Boat level)
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept. Edited to include predictions.
// Time-stamp: <2017-03-10 14:47:00 (slane)>
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

data{
  // Data to be supplied to sampler
  /* Number of (fully observed) observations */
  int<lower=1> N;
  /* Number of (censored) observations */
  int<lower=1> nCens;
  /* Numeric/ordinal predictors */
  real days1[N];
  real days1Cens[nCens];
  real days2[N];
  real days2Cens[nCens];
  real midTrips[N];
  real midTripsCens[nCens];
  /* Categorical predictors */
  int<lower=1> numPaint;
  int<lower=1,upper=numPaint> paintType[N];
  int<lower=1,upper=numPaint> paintTypeCens[nCens];
  int<lower=1> numLoc;
  int<lower=1,upper=numLoc> locID[N];
  int<lower=1,upper=numLoc> locIDCens[nCens];
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
  real beta1;
  real beta2;
  real beta3;
  /* Alphas for categorical */
  vector[numPaint] alpha1;
  vector[numLoc] alpha2;
  vector[numBoat] alphaBoat;
  /* Errors for categorical predictors (i.e. hierarchical now) */
  real<lower=0> sigma_alpha1;
  real<lower=0> sigma_alpha2;
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
    muHat[i] = mu + beta1 * days1[i] + beta2 * days2[i] + beta3 * midTrips[i] + alpha1[paintType[i]] + alpha2[locID[i]] + alphaBoat[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = mu + beta1 * days1Cens[j] + beta2 * days2Cens[j] + beta3 * midTripsCens[j] + alpha1[paintTypeCens[j]] + alpha2[locIDCens[j]] + alphaBoat[boatIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Priors for regression parameters */
  mu ~ normal(0, 10);
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  beta3 ~ normal(0, 5);
  /* Priors for categoricals */
  sigma_alpha1 ~ cauchy(0, 2.5);
  alpha1 ~ normal(0, sigma_alpha1);
  sigma_alpha2 ~ cauchy(0, 2.5);
  alpha2 ~ cauchy(0, sigma_alpha2);
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
  
}
