# BayesGPfit

**An R package for Bayesian Gaussian process regression on regular grid points based on modified exponential sqaured kernel**

* **Install and load the package**
```{r}
devtools::install_github("kangjian2016/BayesGPfit")
library(BayesGPfit)
```


* **Simulate curve on d-dimensional Euclidean space**
```{r}
library(lattice)
set.seed(1224)
dat = list()
dat$x = GP.generate.grids(d=2,num_grids = 100)
curve = GP.simulate.curve.fast(dat$x,a=0.01,b=0.5,poly_degree=20L)
GP.plot.curve(curve,main="Simulated Curve")
```

* **Bayesian model fitting based on two methods**
```{r}
dat$f = curve$f + rnorm(length(curve$f),sd=1)
fast_fit = GP.fast.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
reg_fit = GP.Bayes.fit(dat$f,dat$x,a=0.01,b=0.5,poly_degree=20L,progress_bar = TRUE)
mse = c(reg = mean((reg_fit$f - curve$f)^2),
       fast = mean((fast_fit$f - curve$f)^2))
print(mse)
```
