////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Title: Censored MLE, Model 0
// Author: Steve Lane
// Date: Wednesday, 08 March 2017
// Synopsis: Sampling statements to fit a regression with censored outcome data.
// Includes boat-level intercept. Edited to include predictions. All predictors
// entered.
// Time-stamp: <2017-04-20 15:33:43 (slane)>
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
  real hullSA[N];
  real hullSACens[nCens];
  /* Categorical predictors, entered as matrices of indicators */
  /* Paint type, ablative as base case */
  int<lower=1> numPaint;
  matrix[N, numPaint - 1] paintType;
  matrix[nCens, numPaint - 1] paintTypeCens;
  /* Location of measurement, hull as base case */
  int<lower=1> numLoc;
  matrix[N, numLoc - 1] locID;
  matrix[nCens, numLoc - 1] locIDCens;
  /* Boat type, yacht as base case */
  int<lower=1> numType;
  matrix[N, numType - 1] boatType;
  matrix[nCens, numType - 1] boatTypeCens;
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
  real betaDays1;
  real betaDays2;
  real betaMidTrips;
  real betaHullSA;
  /* Betas for categorical indicators */
  vector[numPaint - 1] betaPaint;
  vector[numLoc - 1] betaLoc;
  vector[numType - 1] betaType;
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
    muHat[i] = mu + betaDays1 * days1[i] + betaDays2 * days2[i] + betaMidTrips * midTrips[i] + betaHullSA * hullSA[i] + locID[i] * betaLoc + paintType[i] * betaPaint + boatType[i] * betaType + alphaBoat[boatID[i]];
  }
  for(j in 1:nCens){
    muHatCens[j] = mu + betaDays1 * days1Cens[j] + betaDays2 * days2Cens[j] + betaMidTrips * midTripsCens[j] + betaHullSA * hullSACens[j] + locIDCens[j] * betaLoc + paintTypeCens[j] * betaPaint + boatTypeCens[j] * betaType + alphaBoat[boatIDCens[j]];
  }
}

model{
  // Model sampling statements
  /* Priors for regression parameters */
  mu ~ normal(0, 5);
  betaDays1 ~ student_t(3, 0, 1);
  betaDays2 ~ student_t(3, 0, 1);
  betaMidTrips ~ student_t(3, 0, 1);
  /* Priors for categorical indicators */
  betaPaint ~ student_t(3, 0, 1);
  betaLoc ~ student_t(3, 0, 1);
  betaType ~ student_t(3, 0, 1);
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
