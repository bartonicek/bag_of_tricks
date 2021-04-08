Handy Algorithms
================
Adam Bartonicek
(last updated: 2021-04-08)

-   [Determinant by row-operations](#determinant-by-row-operations)
-   [Inverse transform random
    sampling](#inverse-transform-random-sampling)
-   [Pseudo-random number generation (linear congruential
    generator)](#pseudo-random-number-generation-linear-congruential-generator)

## Determinant by row-operations

Find the determinant of a square matrix using row-operations. Based on
the following properties of determinants:

1.  Multiplying row *i* by *c* mutliplies the determinant by *c*
2.  Subtracting *c* times row *i* from row *j* does not affect the
    determinant

``` r
det_row <- function(mat) {
  
  nr <- nrow(mat)
  fact <- numeric(nr - 1)
  
  for (i in seq_len(nr - 1)) {
    
    # Normalize current row (1 in first column) & save the multiplication factor
    fact[i] <- mat[i, i]
    mat[i, ] <- mat[i, ] / mat[i, i]
    
    # Subtract multiples of the current row from all succeeding rows
    # to get a column of zeros below
    for (j in (i + 1):nr)
    
      mat[j, ] <- mat[j, ] - mat[j, i] * mat[i, ] 
  
  }
  
  # Take the product of the multiplication factors 
  # & the last value in the matrix (unnormalized)
  prod(fact) * mat[nr, nr]

}

# Example

set.seed(12345)
mat <- matrix(sample(1:9), ncol = 3)

det_row(mat)
```

    ## [1] -95

``` r
det(mat)
```

    ## [1] -95

## Inverse transform random sampling

Generate random samples for an arbitrary distribution by using inverse
cumulative density function (or numerical approximation of it). This is
useful when we have a probability density function (PDF) but no function
to draw samples from it (e.g. like `rnorm()`, `rpois()`, etc…). If we
have cumulative density function (CDF), we can use it directly,
otherwise we can also make a CDF via R’s numerical integration.

(from Ben Lambert’s [Student Guide to Bayesian Statistics problem
sets](https://study.sagepub.com/lambert))

``` r
# Probability density function (normal)
pdf <- function(x) (1 / sqrt(2 * pi)) * exp(-1/2 * (x)^2)

# ...suppose we didn't have rnorm() but had the pdf() above

# Use R's numerical integration to make a CDF
cdf <- function(x) {
  integrate(pdf, 0, x)[[1]]
}

# Make inverse CDF by approximating a vector of cumulative densities
# numerically (via linear approximation)
cumdens <- sapply(seq(0.01, 5, 0.01), cdf)
icdf <- approxfun(cumdens, seq(0.01, 5, 0.01))

# Draw random samples & transform them via ICDF (and also randomly flip sign
# so that we get both positive & negative values)
samples <- runif(10000, 0, 0.5)
x <- icdf(samples) * ifelse(runif(1000, 0, 2) > 1, 1, -1)

par(mfrow = c(1, 2))

plot(seq(-3, 3, 0.1), pdf(seq(-3, 3, 0.1)), type = 'l',
     xlab = 'x', ylab = 'Density', main = 'Probability density function')
hist(x, main = 'Random samples')
box()
```

![](handy_algorithms_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Pseudo-random number generation (linear congruential generator)

Generate “random” samples using the following function:

*s*<sub>*t*</sub> = (*a* ⋅ *s*<sub>*t* − 1</sub> + *b*) mod *M*

where *a*, *b*, and *M* are integers and *s*<sub>*t*</sub> is the
current sample.

(from Ben Lambert’s [Student Guide to Bayesian Statistics problem
sets](https://study.sagepub.com/lambert))

``` r
s0 <- 1
a <- 1597
b <- 51749
M <- 244944

nsamples <- 300

s <- numeric(nsamples)
s[1] <- s0

for (i in 2:nsamples) {
  s[i] = (a * s[i - 1] + b) %% M
}

plot(1:nsamples, s, type = 'l',
     xlab = 'Sample', ylab = 'Value')
```

![](handy_algorithms_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
