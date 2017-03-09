////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Time-stamp: <2017-03-09 15:49:12 (slane)>
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
  /* Observed data */
  real<lower=1.5> Y[N];
  /* Truncated data (brute force, all equal */
  real<upper=1.5> U[nCens];
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
  real alpha1[numPaint];
  real alpha2[numLoc];
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
    muHat[i] = mu + beta1 * days1[i] + beta2 * days2[i] + beta3 * midTrips[i] + alpha1[paintType[i]] + alpha2[locID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = mu + beta1 * days1Cens[j] + beta2 * days2Cens[j] + beta3 * midTripsCens[j] + alpha1[paintTypeCens[j]] + alpha2[locIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Priors for categoricals */
  alpha1 ~ cauchy(0, 2.5);
  alpha2 ~ cauchy(0, 2.5);
  /* Leave as is for the moment (calculate MLE) */
  /* Observed log-likelihood */
  target += lognormal_lpdf(Y | muHat, sigma);
  /* Truncated log-likelihood */
  target += lognormal_lcdf(U | muHatCens, sigma);
}

generated quantities{
  // Statements for predictive outputs, e.g. new data
}
