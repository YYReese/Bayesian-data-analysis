require(stochvol)
require(fBasics)
data("exrates")
print.data.frame(exrates[1:6,])
cat(paste("Data from ", min(exrates$date)," until ",max(exrates$date)))

library(rstan)
options(mc.cores = parallel::detectCores())
model_string_2 <-
  "data{
int<lower=0> n_days;
int<lower=0> n_observed;
int<lower=0> n_pred;
int<lower=0> d; #number of different exchange rates considered, 2 in this question 
matrix[d,n_observed] y_obs;
array[n_days] int obs_or_not;
}

parameters
{
  matrix[d,n_days+n_pred-n_observed+1] y_missing; //y_missing[1] corresponds to y_0
  matrix[d,n_days+n_pred+2] h; //h[1] corresponds to h_0 in the formula
  matrix[d,d] phi;
  vector[d] beta0;
  matrix<lower=-1,upper=1>[d,d] beta1;
  
  corr_matrix[2*d] Omega;        // correlation matrix of Sigma
  vector<lower=0>[2*d] Sqrt_Sigma_eig;    
  // square root of eigenvalues of covariance matrix, with positivity constraint
}

transformed parameters
{
  matrix[2*d,2*d] Sigma=quad_form_diag(Omega, Sqrt_Sigma_eig);
  matrix[d,d] Sigma11=Sigma[1:d,1:d];
  matrix[d,d] Sigma12=Sigma[1:d,(d+1):(2*d)];
  matrix[d,d] Sigma21=Sigma[(d+1):(2*d),1:d];
  matrix[d,d] Sigma22=Sigma[(d+1):(2*d),(d+1):(2*d)];
  matrix[d,d] invSigma22;
  matrix[d,d] M1;
  matrix[d,d] M2;
  matrix[d,d] cholSigma22;
  matrix[d,d] cholM2;

  invSigma22=inverse(Sigma[(d+1):(2*d),(d+1):(2*d)]);
  M1=Sigma12 * invSigma22;
  M2=Sigma11-Sigma12 * invSigma22 * Sigma21;
  cholSigma22=cholesky_decompose(Sigma22);
  cholM2=cholesky_decompose(M2);

  matrix[d,n_days+n_pred+1] eta;
  eta[1:d,1]=rep_vector(1.0,d);
  eta[1:d,2:(n_days+n_pred+1)]=h[1:d,3:(n_days+n_pred+2)]-phi*h[1:d,2:(n_days+n_pred+1)];
  
  matrix[d,n_days+n_pred+1] y;
  //y[1] corresponds to y_0 in the formula
  y[1:d,1]=y_missing[1:d,1];
  { 
    int ind_obs=1;
    int ind_missing=2;
  
    for(i in 1:n_days)
    {
      if(obs_or_not[i]==1)
      {
        y[1:d,i+1]=y_obs[1:d,ind_obs];
        ind_obs=ind_obs+1;
      }
      else
      {
        y[1:d,i+1]=y_missing[1:d,ind_missing];
        ind_missing=ind_missing+1;
      }
    }
    for(i in (n_days+1):(n_days+n_pred))
    {
        y[1:d,i+1]=y_missing[1:d,ind_missing];
        ind_missing=ind_missing+1;
    }
  }
}

model{
  Omega ~ lkj_corr(1);
  Sqrt_Sigma_eig ~ cauchy(0, 0.02);
  
  beta0~normal(0,1);
  h[1:d,1]~normal(0,1);

  for(i in 1:d)
  {
    for(j in 1:d)
    {
      phi[i,j]~uniform(-1,1);
      beta1[i,j]~uniform(-1,1);
    }
  }

  y[1,1]~normal(0.99,0.01);
  y[1,2]~normal(1.61,0.02);

  for(t in 1:(n_days+n_pred))
  {
      h[1:d,t+1]~multi_normal(phi*h[1:d,t],Sigma22);
  }
  for(t in 1:(n_days+n_pred))
  {
    y[1:d,t+1]~multi_normal_cholesky(beta0+beta1 * y[1:d,t]+exp(h[1:d,t+1]/2).* (M1 * eta[1:d,t+1]),
    diag_matrix(exp(h[1:d,t+1]/2))*cholM2);
  }

}

generated quantities
{
  matrix[d,n_observed] y_obs_rep;
  {
    int ind_obs=1;
    
    for(t in 1:n_days)
    {
      if(obs_or_not[t]==1)
      {
        y_obs_rep[1:d,ind_obs]= multi_normal_cholesky_rng(beta0+beta1 * y[1:d,t]+exp(h[1:d,t+1]/2) .* (M1 * eta[1:d,t+1]),
        diag_matrix(exp(h[1:d,t+1]/2))*cholM2);
        ind_obs=ind_obs+1;
      }
    }

  }
}
"


#finding dates up to 2000-04-02
n_days=as.numeric(as.Date("2000-04-02")-as.Date("2000-01-03"))+1;
n_pred=3;
y_obs=t(cbind(exrates$USD[1:n_observed],exrates$GBP[1:n_observed]));
d=2;
obs_or_not=rep(0,n_days)
for (it in 1:n_observed)
{
  date_it=exrates$date[it]-as.Date("2000-01-03")+1
  obs_or_not[date_it]=1
}

data=list(n_days=n_days,n_observed=n_observed, y_obs=y_obs, obs_or_not=obs_or_not,n_pred=n_pred,d=d)

fname="model_1e.stan"
cat(model_string_2,file=fname,append=FALSE)
# list with data and hyperparameters

#passing the model string to STAN
res2<- stan(file = fname, data = data, 
             # Below are optional arguments
             iter = 6000,
             #iter is the number of iterations, including the burn-in
             #the burn-in period is set to iter/2 by default, it can be set to
             chains = 8,cores = parallel::detectCores(),refresh=0);
print(res2,pars=c("phi", "beta0", "beta1", "Sigma"))
