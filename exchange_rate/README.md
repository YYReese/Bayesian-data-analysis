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


## Model 2

we are going to look use a multivariate stochastic
volatility model with leverage to study the USD/EUR and GBP/EUR exchange
rates jointly. The model is described as follows,

$$\begin{aligned}\boldsymbol{y}_t&=\boldsymbol{\beta}_0+\boldsymbol{\beta}_1 \boldsymbol{y}_{t-1}+\exp(h_t/2)\boldsymbol{\epsilon}_t \quad \text{for}\quad 1\le t\le T,\\ 
\boldsymbol{h}_{t+1}&=\boldsymbol{\phi}(\boldsymbol{h}_t)+\boldsymbol{\eta}_t\quad \text{for} \quad 0\le t\le T, \quad h_0\sim N(
\mu, I),\\ (\epsilon_t,\eta_t)&\sim N\left(0, \Sigma\right).\end{aligned}$$

Here I denotes the 2 x 2 identity matrix,
$\boldsymbol{y}_t, \boldsymbol{\beta}_0, \boldsymbol{h}_t, \boldsymbol{\eta}_t, \boldsymbol{\epsilon}_t$
are 2 dimensional vectors, $\boldsymbol{\beta}_1$ and
$\boldsymbol{\phi}$ are 2 x 2 matrices, $\boldsymbol{\Sigma}$ is a
4 x 4 covariance matrix. At each time step $t$, the two components
of $y_t$ will be used to model the USD/EUR and GBP/EUR exchange
rates, respectively.





