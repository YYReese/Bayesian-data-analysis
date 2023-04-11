

Data block
```{r}
n <- nrow(car.data)
body.dum <- dummy(car.data.ctr$veh_body)
agecat.dum <- dummy(car.data.ctr$agecat)
gen.dum <- dummy(car.data.ctr$gender)
area.dum <- dummy(car.data.ctr$area)
car.logit.data <- list(n=n, body=body.dum[,-1],
                       agecat=agecat.dum[,-1],
                       gender=gen.dum[,-1],
                       area=area.dum[,-1],
                       veh_value=car.data.ctr$veh_value,
                       exposure=car.data.ctr$exposure,
                       veh_age=car.data.ctr$veh_age)

car.logit.data2 <- list(n=n, body.convt=body.dum[,2],
                        body.coupe=body.dum[,3],
                        body.hback=body.dum[,4],
                        body.hdtop=body.dum[,5],
                        body.mcara=body.dum[,6],
                        body.mibus=body.dum[,7],
                        body.panvn=body.dum[,8],
                        body.rdstr=body.dum[,9],
                        body.sedan=body.dum[,10],
                        body.stnwg=body.dum[,11],
                        body.truck=body.dum[,12],
                        body.ute=body.dum[,13],
                        agecat2=agecat.dum[,2],
                        agecat3=agecat.dum[,3],
                        agecat4=agecat.dum[,4],
                        agecat5=agecat.dum[,5],
                        agecat6=agecat.dum[,6],
                        gender=gen.dum[,1],
                        area.B=area.dum[,2],
                        area.C=area.dum[,3],
                        area.D=area.dum[,4],
                        area.E=area.dum[,5],
                        area.F=area.dum[,6],
                        veh_value=car.data.ctr$veh_value,
                        exposure=car.data.ctr$exposure,
                        veh_age=car.data.ctr$veh_age)
```

```{r}
car.logit.model2 <- "model {
 beta.mu.0   <- 0
 beta.tau.0  <- 1/(log(15)^2)
 beta.mu <- 0
 beta.tau <- 1/(log(10)/2)^2
 # prior
 beta0   ~ dnorm(beta.mu.0,beta.tau.0)
 beta.veh_value  ~ dnorm(beta.mu,beta.tau)
 beta.exposure ~ dnorm(beta.mu,beta.tau)
 beta.veh_age    ~ dnorm(beta.mu,beta.tau)
 
 beta.gender.m ~ dnorm(beta.mu,beta.tau)
 
 beta.body.convt ~ dnorm(beta.mu,beta.tau)
 beta.body.coupe ~ dnorm(beta.mu,beta.tau)
 beta.body.hback ~ dnorm(beta.mu,beta.tau)
 beta.body.hdtop ~ dnorm(beta.mu,beta.tau)
 beta.body.mcara ~ dnorm(beta.mu,beta.tau)
 beta.body.mibus ~ dnorm(beta.mu,beta.tau)
 beta.body.panvn ~ dnorm(beta.mu,beta.tau)
 beta.body.rdstr ~ dnorm(beta.mu,beta.tau)
 beta.body.sedan ~ dnorm(beta.mu,beta.tau)
 beta.body.stnwg ~ dnorm(beta.mu,beta.tau)
 beta.body.truck ~ dnorm(beta.mu,beta.tau)
 beta.body.ute ~ dnorm(beta.mu,beta.tau)
 
 beta.area.B ~ dnorm(beta.mu,beta.tau)
 beta.area.C ~ dnorm(beta.mu,beta.tau)
 beta.area.D ~ dnorm(beta.mu,beta.tau)
 beta.area.E ~ dnorm(beta.mu,beta.tau)
 beta.area.F ~ dnorm(beta.mu,beta.tau)
 
 beta.agecat.2 ~ dnorm(beta.mu,beta.tau)
 beta.agecat.3 ~ dnorm(beta.mu,beta.tau)
 beta.agecat.4 ~ dnorm(beta.mu,beta.tau)
 beta.agecat.5 ~ dnorm(beta.mu,beta.tau)
 beta.agecat.6 ~ dnorm(beta.mu,beta.tau)
 
 #Likelihood
 for(i in 1:n) {
   logit(mu[i])  <- beta0+beta.veh_value*veh_value[i] +
               beta.exposure*exposure[i] +beta.veh_age*veh_age[i]
               +beta.gender.m*gender[i]
               +beta.body.convt*body.convt[i]
               +beta.body.coupe *body.coupe[i]
               +beta.body.hback *body.hback[i]
               +beta.body.hdtop *body.hdtop[i]
               +beta.body.mcara *body.mcara[i]
               +beta.body.mibus *body.mibus[i]
               +beta.body.panvn *body.panvn[i]
               +beta.body.rdstr *body.rdstr[i]
               +beta.body.sedan *body.sedan[i]
               +beta.body.stnwg *body.stnwg[i]
               +beta.body.truck *body.truck[i]
               +beta.body.ute *body.ute[i]
               +beta.area.B *area.B[i]
               +beta.area.C *area.C[i]
               +beta.area.D *area.D[i]
               +beta.area.E *area.E[i]
               +beta.area.F *area.F[i]
               +beta.agecat.2 *agecat2[i]
               +beta.agecat.3 *agecat3[i]
               +beta.agecat.4 *agecat4[i]
               +beta.agecat.5 *agecat5[i]
               +beta.agecat.6 *agecat6[i]

   clm[i] ~ dbern(mu[i])
   clm.rep[i] ~ dbern(mu[i])
 } 
}"
```
```{r}
# Run JAGS to the completion of the "adaption" stage 
results.car.logit2 <- jags.model(file=textConnection(car.logit.model2), 
                                 data=car.logit.data2,
                                 n.chains=3)

# Burn-in of 10000 iterations
update(results.car.logit2, n.iter=10000)

# Longer run for making inferences, assuming chains have converged
results.car.logit2 <- coda.samples(results.car.logit2,
                                   variable.names=c("clm.rep","beta0",
                                                    "beta.veh_value",
                                                    "beta.exposure" ,"beta.veh_age",
                                                    "beta.gender.m",
                                                    "beta.body.convt",
                                                    "beta.body.coupe",
                                                    "beta.body.hback"  ,
                                                    "beta.body.hdtop" ,
                                                    "beta.body.mcara" ,
                                                    "beta.body.mibus",
                                                    "beta.body.panvn",
                                                    "beta.body.rdstr"  ,
                                                    "beta.body.sedan"  ,
                                                    "beta.body.stnwg"  ,
                                                    "beta.body.truck",
                                                    "beta.body.ute"  ,
                                                    "beta.area.B"  ,
                                                    "beta.area.C" ,
                                                    "beta.area.D" ,
                                                    "beta.area.E" ,
                                                    "beta.area.F"  ,
                                                    "beta.agecat.2",
                                                    "beta.agecat.3",
                                                    "beta.agecat.4"  ,
                                                    "beta.agecat.5"  ,
                                                    "beta.agecat.6"),
                                   n.iter=3000)

# Summary 

```
Model string
```{r}
car.logit.model <- "model {
 beta.mu.0   <- 0
 beta.tau.0  <- 1/(log(15)^2)
 beta.mu <- 0
 beta.tau <- 1/(log(5)/2)^2
 # prior
 beta0       ~ dnorm(beta.mu.0,beta.tau.0)
 beta.veh_value   ~ dnorm(beta.mu,beta.tau)
 beta.exposure ~ dnorm(beta.mu,beta.tau)
 beta.veh_age    ~ dnorm(beta.mu,beta.tau)
 beta.gender ~ dnorm(beta.mu,beta.tau)
 for (i in 1:length(body[1,])){
     beta.body[i] ~ dnorm(beta.mu,beta.tau)
 }
 for (i in 1:length(agecat[1,])){
     beta.agecat[i] ~ dnorm(beta.mu,beta.tau)
 }
 for (i in 1:length(area[1,])){
     beta.area[i] ~ dnorm(beta.mu,beta.tau)
 }
 
 #Likelihood
 for(i in 1:n) {
   logit(mu[i])  <- beta0+beta.veh_value*veh_value[i] +
               beta.exposure*exposure[i] + beta.veh_age*veh_age[i] +
               beta.gender*gender[i] +
               sum(beta.body%*%body[i,]) +
               sum(beta.agecat%*%agecat[i,]) +
               sum(beta.area%*%area[i,])

   clm[i] ~ dbern(mu[i])
 } 
}"
```

```{r}
# Run JAGS to the completion of the "adaption" stage 
results.car.logit <- jags.model(file=textConnection(car.logit.model), 
                                data=car.logit.data,
                                n.chains=3)

# Burn-in of 10000 iterations
update(results.car.logit, n.iter=1000)

# Longer run for making inferences, assuming chains have converged
results.car.logit <- coda.samples(results.car.logit,
                                  variable.names=c("beta0","beta.veh_value",
                                                   "beta.exposure","beta.veh_age",
                                                   "beta.gender","beta.body",
                                                   "beta.agecat","beta.area"), 
                                  n.iter=1000)

# Summary 
summary(results.car.logit)
```

```{r}

#Interpretation of parameter values
#Joining all the chains in one data.frame
results.car.logit.output <- 
  do.call(rbind.data.frame, results.car.logit)
#interpretation
cat("E[beta_LightMed] is",
    mean(results.car.logit.output$`beta.agecat[1]`,3),"\n")

cat("E[beta_LightMed] is",
    mean(results.car.logit.output$`beta.agecat[2]`,3),"\n")

cat("E[beta_LightMed] is",
    mean(results.car.logit.output$`beta.agecat[3]`,3))

#... similarly for the other colours
```
```{r}
#interpretation
ilogit=function(x){ 1/(1+exp(-x))} #Inverse logistic function
cat("E[ilogit(beta0)] is",round(mean(ilogit(results.car.logit.output$beta0)),3),"\n")

```

```{r}
cat("E[ilogit(beta0+beta1)]/E[ilogit(beta0)] is",
    round(mean(ilogit(results.beetles.output$beta0
                      +results.beetles.output$beta1))/
            mean(ilogit(results.beetles.output$beta0)),3),"\n")
```


```{r}

#Create new rows in the dataset for these widths,
#for all 4 colours, with the response Satellites set to NA
newdata <- car.inla.data %>%
  mutate(clm=rep(NA,len=nrow(car.inla.data)))
#Join the new rows with the original dataframe,
#and include the centered width covariate
newdata=rbind(newdata,car.inla.data)

#whenever we use a model with link function that is not the identity
#and we want to compute the posterior marginal and mean of mu_i
#and not eta_i,we need to set link=1 in control.predictor
m3.I <- inla(formula = clm ~ veh_value+exposure+veh_body+
               veh_age+gender+area+agecat, family="binomial",Ntrials=1,
             control.family=list(link="logit"),
             data=newdata, control.fixed=prior.beta,
             control.predictor = list(compute = TRUE,link=1),
             control.compute=list(config=TRUE))

```