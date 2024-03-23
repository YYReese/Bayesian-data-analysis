require(stochvol)
require(fBasics)
data("exrates")
print.data.frame(exrates[1:6,])
cat(paste("Data from ", min(exrates$date)," until ",max(exrates$date)))

library(rstan)
options(mc.cores = parallel::detectCores())

#model in STAN language
model_string_1 <-
  "data{
int<lower=0> n_days;
int<lower=0> n_observed;
int<lower=0> n_pred;
vector[n_observed] y_obs;
array[n_days] int obs_or_not;
}

parameters
{
  vector[n_days+n_pred-n_observed+1] y_missing; //y_missing[1] corresponds to y_0
  vector[n_days+n_pred+2] h_std;
  real<lower=-1,upper=1> rho;
  real<lower=-1,upper=1> phi;
  real<lower=0> sigma2;
  real beta0;
  real <lower=-1,upper=1> beta1;
  real mu;
  real <lower=0> nu; #degrees of freedom
}

transformed parameters
{
  real sigma=sqrt(sigma2);
  real sqrt1mphi2=sqrt(1-phi^2);

  vector[n_days+n_pred+1] h;
  
  h=mu+sigma/sqrt1mphi2*h_std;
  
  vector[n_days+n_pred+1] eta;
  for(t in 1:(n_days+n_pred))
  {
     eta[t+1]=(h_std[t+2]-phi*h_std[t+1])/sqrt1mphi2;
  }
  
  vector[n_days+n_pred+1] y;
  //y[1] corresponds to y_0 in the formula
  y[1]=y_missing[1];
  
  {
  
    int ind_obs=1;
    int ind_missing=2;
  
    for(i in 1:n_days)
    {
      if(obs_or_not[i]==1)
      {
        y[i+1]=y_obs[ind_obs];
        ind_obs=ind_obs+1;
      }
      else
      {
        y[i+1]=y_missing[ind_missing];
        ind_missing=ind_missing+1;
      }
    }
    for(i in (n_days+1):(n_days+n_pred))
    {
        y[i+1]=y_missing[ind_missing];
        ind_missing=ind_missing+1;
    }
    
  }
}



model{
  sigma2~inv_gamma(1e-4,1e-4);
  phi~uniform(-1,1);
  rho~uniform(-1,1);
  mu~normal(-8,1);
  beta0~normal(0,1);
  beta1~uniform(-1,1);
  nu~gamma(2,0.1);

  h_std[1]~normal(0, 1);
  y[1]~normal(0.99,0.01);

  for(t in 1:(n_days+n_pred+1))
  {
     h_std[t+1]~normal(phi*h_std[t],sqrt1mphi2);
  }
  
  for(t in 1:(n_days+n_pred))
  {
     y[t+1]~student_t(nu,beta0+beta1*y[t]+exp(h[t+1]/2)*rho*eta[t+1],exp(h[t+1]/2));
  }

}

generated quantities
{
  vector[n_observed] y_obs_rep;
  
  {
    int ind_obs=1;
 
    for(t in 1:n_days)
    {
      if(obs_or_not[t]==1)
      {
      y_obs_rep[ind_obs]=student_t_rng(nu,beta0+beta1*y[t]+
      exp(h[t+1]/2)*rho*eta[t+1],exp(h[t+1]/2));
      ind_obs=ind_obs+1;
      }
    }
  
  }
}
"


#data

#finding dates up to 2000-04-02
n_observed=max(which(exrates$date<=as.Date("2000-04-02")));
n_days=as.numeric(as.Date("2000-04-02")-as.Date("2000-01-03"))+1
n_pred=3
y_obs=exrates$USD[1:n_observed];
obs_or_not=rep(0,n_days);
for (it in 1:n_observed)
{
  date_it=exrates$date[it]-as.Date("2000-01-03")+1;
  obs_or_not[date_it]=1;
}


data=list(n_days=n_days,n_pred=n_pred,n_observed=n_observed,
          y_obs=y_obs, obs_or_not=obs_or_not);

fname="model_1.stan";
cat(model_string_1,file=fname,append=FALSE);
# list with data and hyperparameters

#passing the model string to STAN
res1<- stan(file = fname, data = data, 
             # Below are optional arguments
             iter = 200000, warmup = 40000,
             #iter is the number of iterations, including the burn-in
             #the burn-in period is set to iter/2 by default, it can be set to
             #something else using the warmup parameter
             chains = 8,cores = parallel::detectCores(),thin=10,refresh=0);
print(res1,pars=c("mu", "rho", "phi", "beta0", "beta1", "sigma", "nu"))



  


y_obs_rep=extract(res1)$y_obs_rep;
y_obs_rep_median=apply(y_obs_rep,1,median)

hist(y_obs_rep_median,col="gray40",
     main="Posterior predictive check for median)", breaks=15)
abline(v=median(y_obs),col="red",lwd=2)

y_obs_rep_skewness=apply(y_obs_rep,1,skewness)
hist(y_obs_rep_skewness,col="gray40",
     main="Posterior predictive check for skewness)", breaks=80,xlim=c(0,1))
abline(v=skewness(y_obs),col="red",lwd=2)

mean_square_diff<-function(v){l=length(v); return(mean((v[2:l]-v[1:(l-1)])^2)); }

y_obs_rep_msqdiff=apply(y_obs_rep,1,mean_square_diff)

hist(y_obs_rep_msqdiff[y_obs_rep_msqdiff<3e-4],col="gray40",
     main="Posterior predictive check for m.s.d.)", breaks=15)
abline(v=mean_square_diff(y_obs),col="red",lwd=2)
