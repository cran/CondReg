#' @name X
#' @title Sample data matrix X for CondReg examples
#' @description n-by-p (200-by-5) matrix of multivariate normal data
#' @docType data
NULL

#' @name sigma
#' @title True covariance matrix of X
#' @description p-by-p (5-by-5) covariance matrix
#' @docType data
NULL

#' Return a vector of grid of penalties for cross-validation
#' @param gridmax maximum value in penalty grid
#' @param numpts number of points in penalty grid
#' @return vector of penalties between 1 and approximately
#'   \code{gridmax} with logarithmic spacing
#' @examples
#' gmax <- 20 ## maximum value for the grid of points
#' npts <- 10 ## number of grid points returned
#' gridpts <- kgrid(gmax,npts)
#' @export
kgrid <- function( gridmax, numpts ){
  gridmax <- round(exp(gridmax-1))
  x <- seq(1, gridmax, length.out=numpts)
  log(1/x) + log(max(x)) + 1
}

#' Compute the best condition number regularized based 
#' based on cross-validation selected penalty parameter
#' @param X n-by-p matrix of data
#' @param k vector of penalties for cross-validation
#' @param ... parameters for \code{select_kmax}
#' @return list of condition number regularized covariance matrix S
#' and its inverse invS
#' @examples
#' ## True covariance matrix
#' \dontrun{sigma <- diag(5)}
#' \dontrun{sigma[3,2] <- sigma[2,3] <- 0.8}
#'
#' ## Generate normal random samples
#' \dontrun{library(MASS)}
#' \dontrun{X <- mvrnorm(200,rep(0,5),sigma)}
#'
#' ## load pre-generated data
#' data(condreg_data)
#'
#' ## Covariance estimation
#' gridpts <- kgrid(50,100)           ## generate grid of penalties to search over
#' crcov <- select_condreg(X,gridpts) ## automatically selects penalty parameter
#'
#' ## Inspect output
#' str(crcov)              ## returned object
#' sigma.hat <- crcov$S    ## estimate of sigma matrix
#' omega.hat <- crcov$invS ## estimate of inverse of sigma matrix
#' @export
select_condreg <- function(X, k, ...){

  n <- nrow(X)
  p <- ncol(X)

  kmax <- select_kmax(X, k, ...)
  
  S <- (t(X) %*% X)/n
  svdS <- svd(S, nu=p)
  
  QLQ <- list(Q=svdS$u, L=svdS$d)
  soln <- condreg(QLQ, kmax$kmax)

  list(S=soln$S, invS=soln$invS, kmax=kmax$kmax)
}

#' Compute shrinkage of eigenvalues for condreg
#' @param L vector of eigenvalues
#' @param k vector of penalties
#' @return list of vector of shrinked eigenvalues \code{Lbar},
#'   optimal u value \code{uopt} and interval indicator \code{intv}.
ml_solver <- function(L, k){

  p <- length(L)
  g <- length(k)
  
  Lbar <- matrix(0, g, p)
  uopt <- rep(0, g)
  intv <- rep(0, g)

  L[L<.Machine$double.eps] <- .Machine$double.eps

  degenindx <- (k > (L[1]/L[p]))
  if (sum(degenindx)>0) {
    Lbar[which(degenindx),] <- matrix(L, sum(degenindx), p, byrow=TRUE)
    uopt[degenindx] <- pmax(1/k[degenindx]/L[p], 1/L[1])
    intv[degenindx] <- TRUE
  }

  if ( any(!degenindx) ){
    kmax1 <- k[!degenindx]
    path <- path_solver(L)

    tmp <- approx(path$k,1/path$u,kmax1)
    uopt <- 1/tmp$y
    
    lambda <- pmin(matrix(rep(kmax1*uopt,p), ncol=p),
                   pmax(matrix(rep(uopt,p), ncol=p),
                        matrix(rep(1/L,length(kmax1)), ncol=p, byrow=TRUE)))

    
    d_shrunk <- 1/lambda
    
    ## cat('mat1:',dim(matrix(rep(kmax1*uopt,p), ncol=p)),'\n')
    ## cat('mat2:',dim(matrix(rep(uopt,p), ncol=p)),'\n')
    ## cat('mat3:',dim(matrix(rep(1/L,g), ncol=p, byrow=TRUE)),'\n')

    ## cat('dim(d_shrunk)',length(lambda),'\n')
    ## cat('degenindx',length(degenindx),'\n')
    ## cat('sum(degenindx)',sum(degenindx),'\n')
    ## cat('sum(!degenindx)',sum(!degenindx),'\n')

    Lbar[!degenindx,] <- d_shrunk
    uopt[!degenindx] <- uopt
    intv[!degenindx] <- FALSE
  }
  
  list(Lbar=Lbar, uopt=uopt, intv=intv)
}

#' Compute optimal u of Lemma 1 in JRSSB paper
#' @param L vector of eigenvalues
path_solver <- function(L) {

  p <- length(L)

  idxzero <- L<.Machine$double.eps
  numzero <- sum(idxzero)
  L[idxzero] <- .Machine$double.eps

  beta <- p

  slope_denom <- sum(L[beta:p])

  u_cur <- 1/L[1]
  v_cur <- 1/L[beta]

  alpha <- 1
  slope_num <- sum(L[1:alpha])

  u <- u_cur
  v <- v_cur
  kmax <- v_cur/u_cur

  isDone <- FALSE

  while (!isDone){

    h_bottom <- 1/L[beta-1]
    v_right <- 1/L[alpha+1]

    u_new <- (slope_num*u_cur + slope_denom*v_cur)/(slope_num+slope_denom)
    v_new <- u_new
    if (u_new < v_right & v_new > h_bottom) {
      isDone <- TRUE
    } else {
      v_new <- h_bottom;
      u_new <- u_cur - slope_denom*(v_new - v_cur)/slope_num;
      if ( u_new < v_right ){	## u_new inside the grid
        beta <- beta - 1;
        slope_denom <- slope_denom + L[beta];
      } else { ## u_new outside the grid. compute for the vertical line u=1/L(alpha+1).
        u_new <- v_right;
        v_new <- v_cur - slope_num*(u_new - u_cur)/slope_denom;
        if ( v_new > h_bottom ) { ## v_new inside the grid
          alpha <- alpha + 1;
          slope_num <- slope_num + L[alpha];
        } else { ## intersection at the bottom right corner of the grid
          u_new <- v_right;
          v_new <- h_bottom;
          alpha <- alpha + 1;
          beta <- beta - 1;
          slope_num <- slope_num + L[alpha];
          slope_denom <- slope_denom + L[beta];
        }
      }
    }
    new_kmax <- v_new/u_new

    if ( any(abs(new_kmax - kmax) > .Machine$double.eps) ) {
      u = c(u_new, u);
      v = c(v_new, v);
      kmax = c(new_kmax, kmax);
    }
    u_cur = u_new;
    v_cur = v_new;
  }

  list(k=kmax, u=u, v=v)
  
}

#' Selection of penalty parameter based on cross-validation
#' @param X n-by-p data matrix
#' @param k vector of penalties for cross-validation
#' @param fold number of folds for cross-validation
#' @export
select_kmax <- function(X, k, fold=min(nrow(X),10)){

  n <- nrow(X)
  p <- ncol(X)
  g <- length(k)

  sections <- cut(1:n,breaks=fold,labels=c(1:fold))
  condmax <- 1

  negloglikelihood <- matrix(0, fold, g);
  
  for (i in c(1:fold)){
    tsindx <- which(sections==i)
    trindx <- which(sections!=i)
    
    Xtr <- X[trindx,,drop=FALSE]
    Xts <- X[tsindx,,drop=FALSE]

    ntr <- nrow(Xtr)
    nts <- nrow(Xts)

    Str <- (t(Xtr)%*%Xtr)/ntr

    soln <- crbulk( Str, k )
    ytest <- Xts %*% soln$Q
    a <- array(ytest^2, dim=c(nts,p,g))
    a <- aperm(a, c(2,1,3))
    b <- array(1/soln$Lbar, dim=c(g,p,nts))
    b <- aperm(b, c(2,3,1))
    ztest <- a * b
    negloglikelihood[i,] <-
      rowSums(aperm(ztest,dim=c(3,1,2)))/nts + rowSums(log(soln$Lbar))

    L <- rep(0,p)
    L[1:min(ntr,p)] <- soln$L[1:min(ntr,p)]
    condmax <- max(condmax, L[1]/L[min(ntr,p)])

  }

  nL <- colSums(negloglikelihood)
  minind <- floor(median(which(nL==min(nL))))
  kmaxopt <- min(k[minind], condmax)

  list(kmax=kmaxopt, negL=nL)
  
}

#' Compute the condition number with given penalty parameter
#' @param data_in data in one of two formats. 1) If a list, elements
#' should be an orthogonal matrix Q and a vector of eigenvalues L from
#' the spectral decomposition of sample covariance matrix of the
#' data. 2) If a n-by-p matrix, each of the n rows are the p variate
#' random observations.
#' @param kmax scalar regularization parameter
#' @return list of condition number regularized covariance matrix s
#' and its inverse invS. 
#' @examples
#' ## True covariance matrix
#' \dontrun{sigma <- diag(5)}
#' \dontrun{sigma[3,2] <- sigma[2,3] <- 0.8}
#'
#' ## Generate normal random samples
#' \dontrun{library(MASS)}
#' \dontrun{X <- mvrnorm(200,rep(0,5),sigma)}
#'
#' ## load pre-generated data
#' data(condreg_data)
#'
#' ## Covariance estimation
#' crcov <- condreg(X,3)
#'
#' ## Inspect output
#' str(crcov)              ## returned object
#' sigma.hat <- crcov$S    ## estimate of sigma matrix
#' omega.hat <- crcov$invS ## estimate of inverse of sigma matrix
#' @export
condreg <- function(data_in, kmax){
  
  if (class(data_in)!='list'){
    n <- nrow(data_in)
    p <- ncol(data_in)
    
    S <- (t(data_in) %*% data_in)/n
    svdS <- svd(S, nu=p)
    
    data_in <- list(Q=svdS$u, L=svdS$d) ## overwrite with spectral decomposition
  }
  
  sol <- ml_solver(data_in$L, kmax)

  Lbar <- as.numeric(sol$Lbar)

  S <- data_in$Q %*% diag(Lbar) %*% t(data_in$Q)
  invS <- data_in$Q %*% diag(1/Lbar) %*% t(data_in$Q)

  if (!sol$intv) u <- sol$uopt
  else           u <- mean(sol$uopt)

  list(S=S, invS=invS)
  
}

#' Computes multiple solutions
#' @param S sample covariance matrix
#' @param k vector of regularization parameters
#' @return list of orthogonal matrix Q, shrinked eigenvalues Lbar
#'   (shrinkage depending on penalty parameters) and
#'   sample eigenvalues L
crbulk <- function(S, k){

  n <- nrow(S)
  p <- ncol(S)
  g <- length(k)
  
  svdS <- svd(S)
  soln <- ml_solver( svdS$d, k )

  ## size of 'soln' components
  ## soln$eigvals: g-by-n
  ## soln$uopt   : g
  ## soln$isInt  : g

  ## if n<p there are some 0 eigenvalues
  if (n<p) {
    soln$Lbar <- cbind(soln$Lbar, matrix(0, g, max(p-n,0)))
  }
  
  list(Q=svdS$u, Lbar=soln$Lbar, L=svdS$d)  
}
