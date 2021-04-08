Handy Algorithms
================
Adam Bartonicek
(last updated: 2021-04-08)

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
