---
title:  "Bayes R Packages"
date:   2020-10-26 20:35:00
categories: bayesian rstan rstanarm rethinking
toc: true
---

# Introduction

In Part 1, I skipped over model diagnostics so I’m going to address that
in this post. I’ll demonstrate *bayesplot* functions here, but
*shinystan* is a nice, interactive alternative.

# Setup Environment

Recall that I’m using the following model for the *mtcars* data set:

\[mpg = a + b*disp + \epsilon\]

``` r
rm(list=ls())

library(tidyverse)
library(rstanarm)
library(rethinking)
library(bayesplot)
library(ggplot2)
library(rstan)
#library(tidybayes)

#knitr::opts_chunk$set(out.width = "50%")
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(message=FALSE)
knitr::opts_chunk$set(warning=FALSE)

options("scipen" = 1, "digits" = 4)

library(datasets)
data(mtcars)
```

# Package *rstanarm*

I’ll demonstrate diagnostic plots first using the model generated from
*stan\_glm* with default priors.

``` r
mdl1 <- stan_glm(mpg ~ disp, data = mtcars)
```

Once the model has been fit, we can use either *as.matrix* or *as.array*
to extract the posterior draws. The key difference is that *as.array*
keeps the chains separate.

``` r
post <- as.array(mdl1)
str(post)
```

    ##  num [1:1000, 1:4, 1:3] 29.4 29.5 28.9 30.4 28.6 ...
    ##  - attr(*, "dimnames")=List of 3
    ##   ..$ iterations: NULL
    ##   ..$ chains    : chr [1:4] "chain:1" "chain:2" "chain:3" "chain:4"
    ##   ..$ parameters: chr [1:3] "(Intercept)" "disp" "sigma"

Note that the default is four chains; we can specify the number of
chains with the *chains=* argument in *stan\_glm*.

## Trace Plots with *bayesplot*

The *bayesplot* package provides the function *mcmc\_trace* which plots
the MCMC draws.

``` r
mcmc_trace(post, pars=c("disp", "sigma"))
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part2_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />
There are three things we are looking for in the trace plot of each
chain:

1.  *Good mixing* - In other words, the chain is rapidly changing values
    across the full region versus getting “stuck” near a particular
    value and slowly changing.

2.  *Stationarity* - The mean of the chain is relatively stable.

3.  *Convergence* - All of the chains spend most of the time around the
    same high-probability values.

The trace plots above look good. However, sometimes it can be hard to
tell when there are multiple chains overlaid on the same plot, so two
alternatives are shown below.

## Trace Plots with *ggplot2*

One alternative is to manually plot each chain separately. Here’s one
way to do it with *ggplot2*.

``` r
library(gridExtra)

pars <- c("disp", "sigma")

plts <- list()
for (par in pars)
{
  df <- as.data.frame(post[,,par]) %>%
    mutate(iteration = row_number()) %>%
    pivot_longer(cols=c(-"iteration"), values_to="value", names_to="chain")

  plts[[par]] <- df %>%
    ggplot() +
    geom_line(mapping=aes(x=iteration, y=value), color="blue") +
    facet_wrap(~chain, ncol=1) +
    labs(title=par)
}

grid.arrange(grobs=plts, nrow=1)
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part2_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

## Trace Rank Plots with *rethinking*

Another alternative is to use the *trankplot* function from the
*rethinking* package. This function plots a trace rank plot which is the
distribution of the ranked samples.

``` r
# Note that trankplot requires the $stanfit parameter from the stan_glm object
trankplot(mdl1$stanfit, pars=c("disp", "sigma"))
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part2_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

## Effective Sample Size

Notice that the *trankplot* function also displays the effective sample
size, *n\_eff*. Since MCMC samples are usually correlated, *n\_eff* is
often less than the number of samples. There is no hard and fast rule
for what is an acceptable number for *n\_eff*. McElreath’s guidance is
it depends on what you are trying to estimate. If you are interested
mostly in the posterior mean, then *n\_eff* = 200 can be enough. But if
you are interested in the tails of the distribution and it’s highly
skewed then you’ll need *n\_eff* to be much larger. There are two
parameters, *iter* and *warmup*, which you can adjust in *stan\_glm* (or
*map2stan* or *stan* itself) if a larger *n\_eff* is needed.

There are a couple other functions besides *trankplot* which will return
*n\_eff* for a *stan\_glm* object.

The *summary* function display *n\_eff* (and a lot of other information)
for the *stan\_glm* object.

``` r
summary(mdl1)
```

    ## 
    ## Model Info:
    ##  function:     stan_glm
    ##  family:       gaussian [identity]
    ##  formula:      mpg ~ disp
    ##  algorithm:    sampling
    ##  sample:       4000 (posterior sample size)
    ##  priors:       see help('prior_summary')
    ##  observations: 32
    ##  predictors:   2
    ## 
    ## Estimates:
    ##               mean   sd   10%   50%   90%
    ## (Intercept) 29.6    1.3 27.9  29.6  31.2 
    ## disp         0.0    0.0  0.0   0.0   0.0 
    ## sigma        3.4    0.5  2.8   3.3   4.0 
    ## 
    ## Fit Diagnostics:
    ##            mean   sd   10%   50%   90%
    ## mean_PPD 20.1    0.9 19.0  20.1  21.2 
    ## 
    ## The mean_ppd is the sample average posterior predictive distribution of the outcome variable (for details see help('summary.stanreg')).
    ## 
    ## MCMC diagnostics
    ##               mcse Rhat n_eff
    ## (Intercept)   0.0  1.0  3480 
    ## disp          0.0  1.0  3328 
    ## sigma         0.0  1.0  3104 
    ## mean_PPD      0.0  1.0  4048 
    ## log-posterior 0.0  1.0  1286 
    ## 
    ## For each parameter, mcse is Monte Carlo standard error, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence Rhat=1).

The *precis* function from the *rethinking* package is another.

``` r
precis(mdl1$stanfit)
```

    ##                    mean      sd      5.5%     94.5% n_eff  Rhat4
    ## (Intercept)    29.55287 1.31308  27.46931  31.60126  3480 1.0000
    ## disp           -0.04098 0.00509  -0.04891  -0.03301  3328 0.9999
    ## sigma           3.37280 0.46295   2.72785   4.18390  3104 0.9995
    ## mean_PPD       20.11456 0.86421  18.77320  21.49677  4048 0.9996
    ## log-posterior -89.51834 1.34816 -91.99463 -88.09986  1286 1.0008

# Package *rstan*

We can do the same diagnostic checks for a model fit with *stan*.

``` r
mdlstan <-
  'data{
      int<lower=1> N;
      real mpg[N];
      real disp[N];
  }
  parameters{
      real a;
      real<lower=-0.1,upper=0> b;
      real<lower=0> sigma;
  }
  model{
      vector[N] mu;
      sigma ~ exponential( 0.17 );
      // b ~ uniform( -0.1 , 0 );
      a ~ normal( 20 , 10 );
      for ( i in 1:N ) {
          mu[i] = a + b * (disp[i] - 230.7);
      }
      mpg ~ normal( mu , sigma );
  }
  generated quantities{
      vector[N] mu;
      for ( i in 1:N ) {
          mu[i] = a + b * (disp[i] - 230.7);
      }
  }'
  
# Allow parallel processing if supported
options(mc.cores = parallel::detectCores())
# Save the compiled Stan program to disk
rstan_options(auto_write = TRUE)

# Note that we must reformat the data as a list
mtcars_dat <- list(N = nrow(mtcars), mpg = mtcars$mpg, disp=mtcars$disp)
mdl2 <- stan(model_code = mdlstan, data=mtcars_dat)
```

The functions for display trace and trace rank plots work the same way
with the *stan* object. The only difference is the parameter names. The
coefficient for *disp* is explicitly named *b* in the *stan* model.

``` r
mcmc_trace(mdl2, pars=c("b", "sigma"))
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part2_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

``` r
trankplot(mdl2, pars=c("b", "sigma"))
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part2_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

# Package *rethinking*

Here’s how to do the same diagnostic checks for a model fit with the
*rethinking* package.

``` r
# Define model

f <- alist(
  mpg ~ dnorm(mu, sigma),
  mu ~ a + b * (disp - 230.7),
  a ~ dnorm(20, 10),
  b ~ dunif(-0.1, 0),
  sigma ~ dexp(0.17)
)

# Fit model - specify chains = 4
mdl3 <- map2stan(f, mtcars, chain=4)
```

Note that again the coefficient for *disp* is explicitly named in the
*rethinking* model. Also note the subtle difference for accessing the
stanfit object (@ versus $). This is because *map2stan* returns an S4
object wheras *stan\_glm* returns an S3 object.

``` r
traceplot(mdl3@stanfit, pars=c("b", "sigma"))
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part2_files/figure-gfm/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

``` r
trankplot(mdl3, pars=c("b", "sigma"))
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part2_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

# Package *rstan*
