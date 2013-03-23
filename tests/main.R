library(CondReg)
data(condreg_data)

#### kmax <- 3
## Covariance estimation
crcov3 <- condreg(X,3)

## Inspect output
sigma.hat3 <- round(crcov3$S,6)    ## estimate of sigma matrix
omega.hat3 <- round(crcov3$invS,6) ## estimate of inverse of sigma matrix

## Known solutions
sigma.hat3.true <- c( 0.860712,  0.033816, -0.005833,  0.004839,  0.032028,
                      0.033816,  0.696350,  0.336740,  0.012817,  0.025672,
                     -0.005833,  0.336740,  0.667493, -0.009572, -0.031694,
                      0.004839,  0.012817, -0.009572,  1.005860, -0.006578,
                      0.032028,  0.025672, -0.031694, -0.006578,  1.011117)
omega.hat3.true <- c( 1.166516, -0.078653,  0.048222, -0.004370, -0.033471,
                     -0.078653,  1.912133, -0.969452, -0.033715, -0.076664,
                      0.048222, -0.969452,  1.992161,  0.031641,  0.085737,
                     -0.004370, -0.033715,  0.031641,  0.994981,  0.008459,
                     -0.033471, -0.076664,  0.085737,  0.008459,  0.994754)
sigma.hat3.true <- matrix(sigma.hat3.true, ncol=ncol(X))
omega.hat3.true <- matrix(omega.hat3.true, ncol=ncol(X))

## Compare computed vs known solutions
if (!all.equal(sigma.hat3.true,sigma.hat3)){ stop("Computed sigma.hat3 is wrong") }
if (!all.equal(omega.hat3.true,omega.hat3)){ stop("Computed omega.hat3 is wrong") }


#### kmax <- 5
## Covariance estimation
crcov5 <- condreg(X,5)

## Inspect output
sigma.hat5 <- round(crcov5$S,6)    ## estimate of sigma matrix
omega.hat5 <- round(crcov5$invS,6) ## estimate of inverse of sigma matrix

## Known solutions
sigma.hat5.true <- c( 0.862413,  0.050283,  0.003136,  0.009209,  0.029314,
                      0.050283,  0.737123,  0.470320,  0.044699,  0.013306,
                      0.003136,  0.470320,  0.696367,  0.017572, -0.051488,
                      0.009209,  0.044699,  0.017572,  1.016169, -0.012329,
                      0.029314,  0.013306, -0.051488, -0.012329,  1.013818)
omega.hat5.true <- c( 1.167868, -0.130368,  0.080884, -0.006587, -0.028030,
                     -0.130368,  2.417887, -1.638750, -0.078200, -0.112141,
                      0.080884, -1.638750,  2.552753,  0.029020,  0.149167,
                     -0.006587, -0.078200,  0.029020,  0.987265,  0.014697,
                     -0.028030, -0.112141,  0.149167,  0.014697,  0.996407)
sigma.hat5.true <- matrix(sigma.hat5.true, ncol=ncol(X))
omega.hat5.true <- matrix(omega.hat5.true, ncol=ncol(X))

## Compare computed vs known solutions
if (!all.equal(sigma.hat5.true,sigma.hat5)){ stop("Computed sigma.hat5 is wrong") }
if (!all.equal(omega.hat5.true,omega.hat5)){ stop("Computed omega.hat5 is wrong") }
