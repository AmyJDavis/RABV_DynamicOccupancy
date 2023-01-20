# RABV_DynamicOccupancy

Code for fitting a custom coded Bayesian hierarchical dynamic occupancy model for rabies virus in raccoons. This code is presented as a companion to the publication "Raccoon rabies control and elimination in the northeastern U.S. and southern Québec, Canada" by Amy J. Davis, Marianne Gagnier, Ariane Massé, Kathleen M. Nelson, Jordona D. Kirby, Ryan Wallace, Xiaoyue Ma, Christine Fehlner-Gardiner, Richard B. Chipman, and Amy T. Gilbert. 

## How to use this code 

The Bayesian hierarchical model that is provided here is custom coded and requires data to be processed for input in this model. The model requires 14 inputs which are described below. This model is based on a 10km x 10km grid overlayed on a study area.  The number of total samples (n) and the number of positive samples (y) should be summarized by grid cell (also referred to as site), season, and surveillance method. The covariates associated with the colonization rates need to be provided in a design matrix, which includes an intercept, and can be spatial and/or temporal in nature. The Markov Chain Monte Carlo (MCMC) iterations are set by the user (n.mcmc) and the tuning and start parameters need to be provided.  

![alt text](https://github.com/AmyJDavis/RABV_DynamicOccupancy/blob/main/ModelFlow1.jpg?raw=true)
Figure 1. Diagram showing the model input elements (data and model inputs) and the model outputs.  The notation is described in the lists below.


  1. y = array of number of rabies positive animals by site, season, and method
  2. n = array of number of animals samples by site, season, and method
  3. Xlocs = data fram of Grid ID values and associated latitude and longitude centroids
  4. catnames = vector of names for the surveillance categories
  5. xgamnames = vector names of covariates on colonization parameter (gamma)
  6. xgam = matrix of covariates for colonization, including at minimum an intercept, columns should match xgamnames
  7. gam.tune = tuning parameter for gamma, can be a single value or a vector of the same length as xgamnames
  8. GWtQU = TRUE/FALSE value as to whether the neighbor effect should be considered (Queen's neighbor)
  9. Gtimesince = TRUE/FALSE value as to whether the time since a cell was occupied should be included as a covariate
  10. a.psi = prior value for alpha for psi1 (initial occupancy)
  11. b.psi = prior value for beta for psi1 (initial occupancy)
  12. a.e = prior value for alpha for epsilon (extinction probability)
  13. b.e = prior value for beta for epsilon (extinction probability) 
  14. n.mcmc = number of MCMC iterations

The model will result in a list of posterior values for:

  1. z = The estimated presence or absence of RABV
  2. psi = Occupancy probabilities
  3. p = Detection probabilities
  4. eps = Local extinction probabilities
  5. gammas = Local colonization probabilities
  6. betagamma = Covariate estimates
