---
title:  "Bayes R Packages"
date:   2020-10-08 19:54:00
categories: bayesian rstan rstanarm rethinking
toc: true
---

# Introduction

The first course I took on Bayesian methods focused mostly on theory,
and since the course was only one semester there wasn’t time to learn
about some of the modern software packages that are commonly used for
Bayesian analysis. This blog post serves as a gentle introduction to
some of these R packages.

Since most of the people in my stats program were only familiar with R
and SAS (and maybe a little Python), I think the following is an easy
way to work up to *rstan* which has a more C-like syntax:

1.  rstanarm
      - Pro: Functions are syntactically very similar to frequentist
        functions with which users are already familiar.
    
      - Pro: Default priors are generally appropriate so the user isn’t
        required to specify priors.
    
      - Con: The user isn’t required to specify priors (i.e., caveat
        emptor).
2.  rethinking
      - Pro: Uses the R formula syntax with which users are already
        familiar.
    
      - Pro: The user is required to specify all priors (i.e., no
        shortcuts).
    
      - Pro: You can get the rstan model out of the rethinking model, so
        this is a nice bridge between R and stan syntax.
    
      - Con: None that I’ve found yet, other than it’s built on top of
        rstan so some folks might prefer to just go right to the source.
3.  rstan
      - Pro: It’s the R interface to stan which is the Bayesian MCMC
        software that runs on multiple platforms and supports multiple
        languages.
    
      - Con: If you aren’t familiar with C, Java or C++ then it’s a
        completely new syntax to learn on top of the Bayesian concepts

# Approach

Fitting a model with any of these packages should include the following
steps:

1.  Fit the model
2.  Examine the prior predictive distribution
3.  *Examine diagnostics*
4.  Examine posterior distribution
5.  Examine the posterior predictive distribution

In this post, I’m only going to illustrate steps 1, 2, 4 and 5. I’ll
save step 3, *Examine Diagnostics* for another time since this post is
already quite lengthy. Plus, many of the same functions and analysis
tools work with the models generated from *rstanarm*,
*rethinking* or *rstan* so it will be more efficient to discuss once
rather than repeating multiple times here.

# Setup Environment

Some basic R environment setup

``` r
rm(list=ls())

library(tidyverse)
library(rstanarm)
library(rethinking)
library(bayesplot)

knitr::opts_chunk$set(out.width = "50%")
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(message=FALSE)
knitr::opts_chunk$set(warning=FALSE)

options("scipen" = 1, "digits" = 4)
```

# Define the Model

I’ll use the mtcars dataset in the examples. To keep things simple, I’m
just going to model *mpg* with a single predictor *disp*.

\[mpg = a + b*disp + \epsilon\]

``` r
library(datasets)
data(mtcars)
head(mtcars)
```

    ##                    mpg cyl disp  hp drat    wt  qsec vs am gear carb
    ## Mazda RX4         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
    ## Mazda RX4 Wag     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
    ## Datsun 710        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
    ## Hornet 4 Drive    21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
    ## Hornet Sportabout 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
    ## Valiant           18.1   6  225 105 2.76 3.460 20.22  1  0    3    1

``` r
mtcars %>%
  ggplot(aes(x=disp, y=mpg)) +
  geom_point(aes(color=factor(cyl))) +
  stat_smooth(method="lm")
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-3-1.png" width="50%" style="display: block; margin: auto;" />

Note that a linear model isn’t a great fit to the data–I’ll deal with
that in the next blog post. This post is focused on the mechanics of how
to use each package.

One more thing, let’s calculate the mean and standard deviation of both
*mpg* and *disp*. We’ll need this info later.

``` r
mu <- mtcars %>% select(mpg, disp) %>% colMeans()
sigma <- mtcars %>% select(mpg, disp) %>% apply(2,sd)

knitr::kable(cbind(mu, sigma), col.names = c("Mean", "Std Dev"))
```

|      |   Mean | Std Dev |
| :--- | -----: | ------: |
| mpg  |  20.09 |   6.027 |
| disp | 230.72 | 123.939 |

# Package rstanarm

We will use the *stan\_glm* function from the *rstanarm* package for the
linear model. As you’ll see, the syntax is very similar to *lm*.

## Default Priors

Let’s start with the default priors. When using the default priors,
*stan\_glm* automatically standardizes the parameters so we don’t need
to do that beforehand.

### Fit Model

``` r
mdl1 <- stan_glm(mpg ~ disp, data = mtcars, cores=parallel::detectCores())
```

### Prior Predictive Distribution

We can check whether or not the defaults priors seem reasonable with the
prior predictive distribution. The *prior\_summary* function shows the
default priors for the model as well as the adjusted priors after
automatic scaling. See
<http://mc-stan.org/rstanarm/articles/priors.html> if you are interested
in the details about the default and adjusted priors.

``` r
prior_summary(mdl1)
```

    ## Priors for model 'mdl1' 
    ## ------
    ## Intercept (after predictors centered)
    ##   Specified prior:
    ##     ~ normal(location = 20, scale = 2.5)
    ##   Adjusted prior:
    ##     ~ normal(location = 20, scale = 15)
    ## 
    ## Coefficients
    ##   Specified prior:
    ##     ~ normal(location = 0, scale = 2.5)
    ##   Adjusted prior:
    ##     ~ normal(location = 0, scale = 0.12)
    ## 
    ## Auxiliary (sigma)
    ##   Specified prior:
    ##     ~ exponential(rate = 1)
    ##   Adjusted prior:
    ##     ~ exponential(rate = 0.17)
    ## ------
    ## See help('prior_summary.stanreg') for more details

``` r
# Plot prior predictive distribution
N <- 100

prior_samples <- data.frame(a = rnorm(N, 20, 15),
                            b = rnorm(N, 0, 0.12))

D <- seq(min(mtcars$disp), max(mtcars$disp), length.out = N)

res <- as.data.frame(apply(prior_samples, 1, function(x) x[1] + x[2] * (D-230.7))) %>%
  mutate(disp = D) %>%
  pivot_longer(cols=c(-"disp"), names_to="iter") 

res %>%
  ggplot() +
  geom_line(aes(x=disp, y=value, group=iter), alpha=0.2) +
  labs(x="disp", y="prior predictive mpg")
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-7-1.png" width="50%" style="display: block; margin: auto;" />

Two observations from this plot stand out: 1) negative mpg is
unrealistic and 2) increasing mpg as displacement increases also seems
unlikely in the real-world. Later on I’ll choose a more informative
prior that incorporates this additional knowledge. But the adjusted
default priors aren’t totally unreasonable so I’ll proceed with the
analysis.

### Posterior Distribution

The Bayesian posterior point estimates for *a* and *b* are shown below.

``` r
coef(mdl1)
```

    ## (Intercept)        disp 
    ##    29.61562    -0.04135

The 89% credible intervals for all *a*, *b* and *sigma* are shown below.

``` r
knitr::kable(posterior_interval(mdl1, prob=0.89))
```

|             |     5.5% |    94.5% |
| :---------- | -------: | -------: |
| (Intercept) |  27.5252 |  31.6364 |
| disp        | \-0.0491 | \-0.0336 |
| sigma       |   2.7228 |   4.1633 |

### Posterior Predictive Distribution

Finally, let’s check the posterior predictive distribution.

``` r
newdata <- data.frame(disp=seq(min(mtcars$disp), max(mtcars$disp)))

y_rep <- as.data.frame(t(posterior_linpred(mdl1, newdata=newdata, draws=20))) %>%
  cbind(newdata) %>%
  pivot_longer(cols=starts_with("V"), names_to="grp", values_to="mpg")

y_rep %>%
  ggplot(aes(x=disp, y=mpg)) +
  geom_line(aes(group=grp), alpha=0.2) +
  geom_point(data = mtcars, aes(color=factor(cyl))) 
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-10-1.png" width="50%" style="display: block; margin: auto;" />

Given our assumption of a linear model (which we already know isn’t
really appropriate for this data), the resulting posterior predictive
distribution is consistent with the observed data.

## User-Specified Priors

This time I’ll specify priors instead of using the defaults. But first,
I’ll standardize both *mpg* and *disp*.

``` r
# Standardize
df <- data.frame(mtcars %>% select(mpg, disp) %>% scale())
df['cyl'] = mtcars$cyl

mdl2 <- stan_glm(mpg ~ disp, data = df,
                 prior = normal(0,1/sqrt(2)), # prior for slope
                 prior_intercept = normal(0,1/sqrt(2)), # prior for intercept
                 cores=parallel::detectCores())
```

### Prior Predictive Distribution

Again, let’s do a sanity check with the prior predictive distribution.

``` r
prior_summary(mdl2)
```

    ## Priors for model 'mdl2' 
    ## ------
    ## Intercept (after predictors centered)
    ##  ~ normal(location = 0, scale = 0.71)
    ## 
    ## Coefficients
    ##  ~ normal(location = 0, scale = 0.71)
    ## 
    ## Auxiliary (sigma)
    ##  ~ exponential(rate = 1)
    ## ------
    ## See help('prior_summary.stanreg') for more details

``` r
# Plot prior predictive distribution
N <- 100

prior_samples <- data.frame(a = rnorm(N, 0, 1/sqrt(2)),
                            b = rnorm(N, 0, 1/sqrt(2)))

D <- seq(min(df$disp), max(df$disp), length.out = N)

res <- as.data.frame(apply(prior_samples, 1, function(x) x[1] + x[2] * (D))) %>%
  mutate(disp = D) %>%
  pivot_longer(cols=c(-"disp"), names_to="iter") 

res %>%
  ggplot() +
  geom_line(aes(x=disp, y=value, group=iter), alpha=0.2) +
  labs(x="disp", y="prior predictive mpg")
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-13-1.png" width="50%" style="display: block; margin: auto;" />

Remember, I standardized *mpg* & *disp* so that’s why the scales are
different in this plot. Also, the negative values of *mpg* aren’t
necessarily unrealistic after standardization. However, the unrealistic situations 
where *mpg* increases as *disp* increases are still present. 
This is because the prior I chose for *b* is still symmetric
about 0. In a later example, I’ll choose yet another prior for *b* that
is even further refined based on real-world information.

### Posterior Distribution

Let’s check the posterior estimates:

``` r
mdl2
```

    ## stan_glm
    ##  family:       gaussian [identity]
    ##  formula:      mpg ~ disp
    ##  observations: 32
    ##  predictors:   2
    ## ------
    ##             Median MAD_SD
    ## (Intercept)  0.0    0.1  
    ## disp        -0.8    0.1  
    ## 
    ## Auxiliary parameter(s):
    ##       Median MAD_SD
    ## sigma 0.5    0.1   
    ## 
    ## ------
    ## * For help interpreting the printed output see ?print.stanreg
    ## * For info on the priors used see ?prior_summary.stanreg

And the 89% posterior credible intervals:

``` r
posterior_interval(mdl2, prob=0.89)
```

    ##                5.5%   94.5%
    ## (Intercept) -0.1580  0.1517
    ## disp        -0.9959 -0.6654
    ## sigma        0.4496  0.6862

The above are standardized, so let’s convert back to the orginal scale
and compare to the results using the defaults priors.

``` r
a_prime_mdl2 <- mu['mpg'] + sigma['mpg']*coef(mdl2)[1] - coef(mdl2)[2] * sigma['mpg'] * mu['disp'] / sigma['disp']
b_prime_mdl2 <- coef(mdl2)[2]*sigma['mpg'] / sigma['disp']

knitr::kable(cbind(coef(mdl1), c(a_prime_mdl2, b_prime_mdl2)), 
             col.names = c("Default", "User-Specified"))
```

|             |  Default | User-Specified |
| :---------- | -------: | -------------: |
| (Intercept) |  29.6156 |        29.4235 |
| disp        | \-0.0414 |       \-0.0404 |

Voila\! The results are very similar as expected.

### Posterior Predictive Distribution

Finally, let’s check the posterior predictive distribution using the
*posterior\_linepred* function.

``` r
newdata <- data.frame(disp=seq(min(df$disp), max(df$disp)))

y_rep <- as.data.frame(t(posterior_linpred(mdl2, newdata=newdata, draws=20))) %>%
  cbind(newdata) %>%
  pivot_longer(cols=starts_with("V"), names_to="grp", values_to="mpg")

y_rep %>%
  ggplot(aes(x=disp, y=mpg)) +
  geom_line(aes(group=grp), alpha=0.2) +
  geom_point(data = df, aes(color=factor(cyl))) 
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-17-1.png" width="50%" style="display: block; margin: auto;" />

And again, the results are consistent with our assumptions and
expectations.

# Package rethinking

## Original data

Again, I’ll start with the original data. First, define the model as
shown below. Note that the *rethinking* package requires you to define
all priors–there are no defaults. I’ll use the same priors for *a* and
*sigma* as *rstanarm’s* adjusted default priors, but now I’ll use a
uniform\[-0.1, 0\] prior for *b*.

``` r
# Define model

f <- alist(
  mpg ~ dnorm(mu, sigma),
  mu ~ a + b * (disp - 230.7),
  a ~ dnorm(20, 10),
  b ~ dunif(-0.1, 0),
  sigma ~ dexp(0.17)
)
```

``` r
# Fit model
mdl3 <- map2stan(f,mtcars)
```

### Prior Predictive Distribution

You’ll see the effect of my choice of priors in the prior predictive
distribution plot below.

``` r
# Plot prior predictive distribution
N <- 100

prior_samples <- data.frame(a = rnorm(N, 25, 15),
                            b = runif(N, -0.1, 0))

D <- seq(min(mtcars$disp), max(mtcars$disp), length.out = N)

res <- as.data.frame(apply(prior_samples, 1, function(x) x[1] + x[2] * (D))) %>%
  mutate(disp = D) %>%
  pivot_longer(cols=c(-"disp"), names_to="iter") 

res %>%
  ggplot() +
  geom_line(aes(x=disp, y=value, group=iter), alpha=0.2) +
  labs(x="disp", y="prior predictive mpg")
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-20-1.png" width="50%" style="display: block; margin: auto;" />

Note that now all of the slopes are non-positive. This reflects my prior
belief that *mpg* cannot increase as *disp* increases.

### Posterior Distribution

The *precis* function from the *rethinking* package gives us the point
estimate, credible intervals and some additional information.

``` r
precis(mdl3)
```

    ##           mean       sd     5.5%    94.5% n_eff Rhat4
    ## a     20.08232 0.591827 19.17169 21.04652 886.7 1.001
    ## b     -0.04105 0.004896 -0.04896 -0.03339 758.8 1.007
    ## sigma  3.35813 0.440448  2.72045  4.11441 748.3 1.000

### Posterior Predictive Distribution

Finally, we check the posterior predictive distribution using the
*extract.samples* function.

``` r
N <- 20
ppd <- as.data.frame(extract.samples(mdl3, N)) %>%
  mutate(x_lwr = c(rep(min(mtcars$disp),N)),
         x_upr = c(rep(max(mtcars$disp), N)),
         grp = 1:N) %>%
  pivot_longer(cols=starts_with("x_"), names_to="x", values_to="disp") %>%
  mutate(y = a + b * (disp - 230.7))


ggplot(data=mtcars, mapping=aes(x=disp, y=mpg)) +
  geom_point(aes(color=factor(cyl))) +
  geom_line(data=ppd, mapping=aes(x=disp, y=y, group=grp), color="black", alpha=0.2)
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-22-1.png" width="50%" style="display: block; margin: auto;" />

## Standardized data

Now, we’ll standardize the \(mpg\) and \(disp\) and then define the
model as follows

``` r
# Standardize 
df <- as.data.frame(mtcars %>% select(mpg, disp) %>% scale())
df['cyl'] <- mtcars$cyl

# Define model
f <- alist(
  mpg ~ dnorm(mu, sigma),
  mu ~ a + b * disp,
  a ~ dnorm(0,0.7), # map2stan doesn't like 1/sqrt(2) here so I use 0.7 instead
  b ~ dnorm(0,0.7),
  sigma ~ dexp(1)
)

# Fit model
mdl4 <- map2stan(f,df)
```

### Prior Predictive Distribution

Same as [*rstanarm* prior predictive distribution.](#prior-predictive-distribution-1)

### Posterior Distribution

``` r
precis(mdl4)
```

    ##            mean      sd    5.5%   94.5%  n_eff  Rhat4
    ## a      0.005393 0.09108 -0.1359  0.1529 1150.8 0.9996
    ## b     -0.829030 0.10559 -0.9938 -0.6524  814.6 1.0056
    ## sigma  0.555785 0.07629  0.4468  0.6818  929.2 0.9990

Again, let’s convert back to the original scale for comparison.

``` r
# Let's convert back to original scale for comparison
a_prime_mdl4 <- mu['mpg'] + sigma['mpg']*coef(mdl4)['a'] - coef(mdl4)['b'] * sigma['mpg'] * mu['disp'] / sigma['disp']
b_prime_mdl4 <- coef(mdl4)['b']*sigma['mpg'] / sigma['disp']

knitr::kable(cbind(coef(mdl3)[1:2], c(a_prime_mdl4, b_prime_mdl4)), 
             col.names = c("Default", "User-Specified"))
```

|   | Default | User-Specified |
| :- | ------: | -------------: |
| a |  20.082 |        29.4246 |
| b | \-0.041 |       \-0.0403 |

# Package rstan

First, he’s a handy [cheat
sheet](http://www.sumsar.net/files/posts/2017-bayesian-tutorial-exercises/stan_cheat_sheet2.12.pdf)
on *stan’s* syntax.

## Original Data

Remember way at the top when I said you can get the *rstan* model from
the *rethinking* model? Here it is:

``` r
#Below is output from stancode(mdl2)

(mdlstan <- stancode(mdl3))
```

    ## //2020-10-08 20:12:29
    ## data{
    ##     int<lower=1> N;
    ##     real mpg[N];
    ##     real disp[N];
    ## }
    ## parameters{
    ##     real a;
    ##     real<lower=-0.1,upper=0> b;
    ##     real<lower=0> sigma;
    ## }
    ## model{
    ##     vector[N] mu;
    ##     sigma ~ exponential( 0.17 );
    ##     // b ~ uniform( -0.1 , 0 );
    ##     a ~ normal( 20 , 10 );
    ##     for ( i in 1:N ) {
    ##         mu[i] = a + b * (disp[i] - 230.7);
    ##     }
    ##     mpg ~ normal( mu , sigma );
    ## }
    ## generated quantities{
    ##     vector[N] mu;
    ##     for ( i in 1:N ) {
    ##         mu[i] = a + b * (disp[i] - 230.7);
    ##     }
    ## }

    ## [1] "//2020-10-08 20:12:29\ndata{\n    int<lower=1> N;\n    real mpg[N];\n    real disp[N];\n}\nparameters{\n    real a;\n    real<lower=-0.1,upper=0> b;\n    real<lower=0> sigma;\n}\nmodel{\n    vector[N] mu;\n    sigma ~ exponential( 0.17 );\n    // b ~ uniform( -0.1 , 0 );\n    a ~ normal( 20 , 10 );\n    for ( i in 1:N ) {\n        mu[i] = a + b * (disp[i] - 230.7);\n    }\n    mpg ~ normal( mu , sigma );\n}\ngenerated quantities{\n    vector[N] mu;\n    for ( i in 1:N ) {\n        mu[i] = a + b * (disp[i] - 230.7);\n    }\n}\n\n"

Now let’s fit this model with *rstan*.

``` r
# Allow parallel processing if supported
options(mc.cores = parallel::detectCores())
# Save the compiled Stan program to disk
rstan_options(auto_write = TRUE)

# Note that we must reformat the data as a list
mtcars_dat <- list(N = nrow(mtcars), mpg = mtcars$mpg, disp=mtcars$disp)
fit <- stan(model_code = mdlstan, data=mtcars_dat)
```

### Prior Predictive Distribution

Same as [*rethinking* prior predictive distribution.](#prior-predictive-distribution-2)

### Posterior Distribution

``` r
print(fit, pars=c("a", "b", "sigma"), probs=c(0.055, 0.945))
```

    ## Inference for Stan model: ea0216f179c9cc087250ab0ac3760980.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##        mean se_mean   sd  5.5% 94.5% n_eff Rhat
    ## a     20.10    0.01 0.58 19.17 21.03  3327    1
    ## b     -0.04    0.00 0.00 -0.05 -0.03  3199    1
    ## sigma  3.36    0.01 0.46  2.70  4.16  3223    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Oct 08 20:22:49 2020.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

``` r
plot(fit, pars=c("a", "b", "sigma"), ci_level=0.89)
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-29-1.png" width="50%" style="display: block; margin: auto;" />

### Posterior Predictive Distribution

``` r
# extract method keeps the chains separate, 
# as.data.frame merges them

N <- 20
ppd <- as.data.frame(fit, pars=c("a", "b")) %>% # Returns 1000 rows
  head(N) %>%  # keep first N rows
  mutate(x_lwr = c(rep(min(mtcars$disp),N)),
         x_upr = c(rep(max(mtcars$disp), N)),
         grp = 1:N) %>%
  pivot_longer(cols=starts_with("x_"), names_to="x", values_to="disp") %>%
  mutate(y = a + b * (disp - 230.7))


ggplot(data=mtcars, mapping=aes(x=disp, y=mpg)) +
  geom_point(aes(color=factor(cyl))) +
  geom_line(data=ppd, mapping=aes(x=disp, y=y, group=grp), color="black", alpha=0.2)
```

<img src="{{site.baseurl}}/images/Bayes_Pkgs_Part1_files/figure-gfm/unnamed-chunk-30-1.png" width="50%" style="display: block; margin: auto;" />

# Summary

Let’s compare the posterior point estimates for the three packages.
Recall that I used different priors for *b* in some of the examples, so
I’m summarizing the results into two tables to make sure I’m doing
apples to apples comparisons.

The first table shows the posterior estimates for a Normal prior for
*b*, and the results are very close for all three packages as expected.

|             | rstanarm (default prior) | rstanarm (user specified) | rethinking |
| :---------- | -----------------------: | ------------------------: | ---------: |
| (Intercept) |                  29.6156 |                   29.4235 |    29.4246 |
| disp        |                 \-0.0414 |                  \-0.0404 |   \-0.0403 |


The next table shows the posterior estimates when the prior for *b* was
U\[-0.1, 0\]. As expected, the results are very similar for the two
packages. Note: If you are wondering why I didn’t demonstrate *rstanarm* with this
prior, it’s because *stan\_glm* has a limited set of user-specified
prior distributions it supports.


|   | rethinking |    rstan |
| :- | ---------: | -------: |
| a |     20.082 |  20.1012 |
| b |    \-0.041 | \-0.0412 |


Finally, by comparing the two tables, we get a sense of how the choice
of prior affects the posterior. It only really affects the posterior
estimate of the intercept which is generally not of interest. The
posterior estimate for *b* is essentially the same for both priors which
tells us that observed data is driving the estimate, not the priors.
Thus, *assuming a linear model is appropriate for this data*, I would be
fairly confident that this posterior estimate for *b* is not sensitive
to my choice of prior.
