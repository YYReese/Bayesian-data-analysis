## Introduction
In this project, we study a dataset about currency exchange rates. The
exrates dataset of the stochvol package contains the daily average
exchange rates of 24 currencies versus the EUR, from 2000-01-03 until
2012-04-04. We are going to fit a various stochastic volatility
models on this dataset

## Model 1
Consider the following leveraged Stochastic Volatility
(SV) model.

$$\begin{aligned} y_t&=\beta_0+\beta_1 y_{t-1}+\exp(h_t/2)\epsilon_t \quad \text{for}\quad 1\le t\le T,\\ h_{t+1}&=\mu+\phi(h_t-\mu)+\sigma \eta_t\quad \text{for} \quad 0\le t\le T, \quad h_0\sim N(\mu, \sigma^2/(1-\phi^2)),\(\epsilon_t,\eta_t)&\sim N\left(0, \Sigma_{\rho}\right)\quad \text{ for } \quad \Sigma_{\rho}=\left(\begin{matrix} 1 & \rho\\ \rho & 1\end{matrix}\right). \end{aligned}$$

Here $t$ is the time index, $y_t$ are the observations (such
as daily USD/EUR rate), $h_t$ are the log-variance process,
$\epsilon_t$ is the observation noise, and $\eta_t$ is the
log-variance process noise (which are correlated, but independent for
different values of $t$). The hyperparameters are
$\beta_0, \beta_1, \mu, \phi, \sigma, \rho$.

For stability, it is necessary to have $\phi\in (-1,1)$, and by
the definition of correlation matrices, we have $\rho\in [-1,1]$.

First, we implement the model in JAGS. Here the unobserved components of
$y$ will be included as NA in the data, and JAGS will automatically
create stochastic nodes for them.

We start by explaining our prior choices. We use
Gamma($10^{-4},10^{-4}$) prior for $\tau=1/\sigma^2$, as this has mean
$1$ and variance $10^4$, making it quite uninformative. We use
uniform(-1,1) priors for parameters $\rho$ and $\phi$, as they are
constrained in the interval $(-1,1)$ according to the question
statement. We also use uniform $(-1,1)$ prior for $\beta_1$, since
values larger than 1 would be unstable (i.e. $y_t$ would grow
exponentially in $t$), and we do not believe that this is expected of
the major currencies we are modelling here. For $\beta_0$, we expect
this parameter to be on the order of magnitude 1 for $USD/EUR$, so we
set a $N(0,1)$ prior. We also need to put a prior distribution on $y_0$,
since this is part of the model for the first observation $y_1$. Based
on the last few days of USD/EUR exchange rate in December 1999 available
on <https://www.federalreserve.gov/releases/h10/20000103>, the rate is
around $0.99$ with low variability, so we select a prior
$y_0\sim N(0.99, 0.01^2)$.

In JAGS, multivariate nodes cannot be partially observed. Hence we do
not use multivariate normal distributions in our implementation.
Instead, we need to compute the conditional distribution of $\epsilon_t$
given $\eta_t$. Since these random variables are jointly normal with
covariance matrix
$\Sigma_{\rho}=\left(\begin{matrix}1 & \rho\\ \rho & 1\end{matrix}\right)$,
based on the formula on <https://statproofbook.github.io/P/mvn-cond>,
the conditional distribution of
$\epsilon_t|\eta_t\sim N(\rho \eta_t, 1-\rho^2)$. This is used in our
definition of $\epsilon_t$ in the model code.


## Model 2




