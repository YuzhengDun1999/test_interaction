#' @importFrom CompQuadForm davies
#' @importFrom glmnet cv.glmnet
# Change for October 4, 2019
# Use tryCatch instead of try in ridge.select.linear when calling solve.

# Jan 29, 2014
# v12: change from "Davies" to "davies" and "Liu" to "liu"
# Version 11: Jan 27, 2014
# identical to V10, except BT functions removed
# Verstion 10: Oct 1, 2013
# identical to Version 9, only SKAT_davies() function modified so that .qfc() calls from iSKAT
# Version 9: Sept 25, 2013
# identical to Version 8, only SKAT_davies() function modified so that .qfc() calls from SKAT
# this is to resolve namespace problems.
# Version 8: Jan 28, 2013
# change to GxEscore.linear.GCV() such that varhat can accomodate larger p
# Version 7: Dec 20, 2012
# change to ridge.select.linear() and ridge.linear()
# such that an intercept is now included but Y is not centered
# Version 6: Dec 17, 2012
# the only change is to add a try() function in ridge.select.linear()
# such that the GCV will only select a model that converges/ matrix is invertible
# Version 5: March 15, 2012
# Changes in V5:
# GxEscore() renamed to GxEscore.linear.GCV() for consistency compared to logistic model
# chooseridge() renamed to chooseridge.linear()
# ridge() renamed to ridge.linear()
# ridge.select() renamed to ridge.select.linear()
# GxEscore.linear.GCV() implements Davies method from v5 onwards
# GxEscore.linear.GCV have changes made to speed up computation, adds 3 arugments related to scale and weights
# chooseridge.linear() modified slightly, add 3 other arguments center.Z, scale.Z, and weights.Z
# ridge.linear and ridge.select.linear() modified to add 3 arguments center.Z, scale.Z, and weights.Z, matrix computation made faster
# v5 adds 3 helper functions to get davies p-value which are called by GxEscore.linear.GCV()
# SKAT_davies()
# Get_Lambda() # modified to include only.values=T in eigen()
# Get_PValue_GESAT()
# All 3 helper functions were copied unmodified from GxE-scoretest-logistic-snpset-v19.R except Get_Lambda() #modified to include only.values=T in eigen()

# Version 4 sets the upper limit of lambda to be 9*floor(sqrt(n)), where n=sample size, modification made in GxEscore()
# Version 3 sets the upper limit of lambda to be 3*floor(sqrt(n)), where n=sample size, modification made in GxEscore() and chooseridge()
# Version 2 was modified so that GxEscore() calculates min p-value as well, changed the input to function
# Version 1 adapted from GxE-scoretest-v9.R
# converts it into a function
# chooseridge() was modified such that lambda=0 is not an option
# p-value calculated without perturbations

#----------------------------------------------------------------------------------------------------------
# Note: these are the functions for linear (v5) and logistic regression (v20) respectively:
# GxEscore.logistic.GCV()                            : main function
# ridge.logistic()                                   : fit final null model
# chooseridge.logistic() and ridge.select.logistic() : select ridge parameter
# Burdentests.GE.logistic                            : Burden tests

# GxEscore.linear.GCV()                              : main function
# ridge.linear()                                     : fit final null model
# chooseridge.linear() and ridge.select.linear()     : select ridge parameter
# Burdentests.GE.linear                              : Burden tests

# The analogous functions in each take in exactly the same arguments in linear and logistic codes
# ridge.logistic() returns slightly different values from ridge.linear()
# other analogous functions return the same values in both linear and logistic codes
#----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
# START: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------
GxEscore.linear.GCV <- function(Y, Xtilde, Z, V, ridge.penalty.factor=rep(1, ncol(Z)), lasso.select=F, lasso.criterion="lambda.min", lasso.ols=F, ols=F, type="davies",
                                lower=NULL, upper=NULL, nintervals=NULL, plotGCV=F, plotfile=NA, scale.Z=T, weights.Z=NULL, weights.V=NULL){

  # Y (n x 1 matrix):  continuous outcome variable
  # Xtilde (n x qtilde matrix): are the variables adjusted for (not penalized)
  # Do NOT include intercept as part of Xtilde (Y is centered in ridge, so no intercept needed)
  # Z (n x p matrix):  genetic covariates that are adjusted for (penalized)
  # V (n x p matrix): GxE terms that we are testing
  # ridge.penalty.factor: Separate ridge penalty factors can be applied to each coefficient. NaN means excluding the variables
  # lasso.criterion: How to select best lasso coefficients. "lambda.1se" or "lambda.min"
  # lasso.ols: refit OLS after using lasso select variables.
  # n = no. of samples
  # lower cannot be zero
  # to use a fixed lambda, set nintervals=1
  # NB: if nintervals=1, upper is used and lower is ignored

  if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  if(nrow(V)!= nrow(Y)) stop("dimensions of V and Y don't match")
  if(type!="davies"&&type!="liu") stop("type has to be either davies or liu")
  if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")
  if(length(ridge.penalty.factor)!= ncol(Z)) stop("dimensions of penalty factors and Z don't match")
  
  n <- drop(nrow(Y))

  #---------------------------------------------------------------------------
  # fit ridge regression model under the null
  # Note that all variables should always be centered as otherwise scaling by scale() will be incorrect
  #---------------------------------------------------------------------------
  if (ols | lasso.ols){
    upper = 0
    nintervals = 1
  }
  if(nintervals>1){
    if(is.null(weights.Z)==F & scale.Z==F){
      Z <- t(t(Z) * (weights.Z))
    }
    Z <-  scale(Z, center=T, scale=scale.Z)
    if (lasso.select){
      lasso.fit = cv.glmnet(x = cbind(Xtilde, Z), y = Y, family = c("gaussian"),
                            alpha = 1, penalty.factor = c(rep(0, ncol(Xtilde)), rep(1, ncol(Z))))
      beta_lasso = c(coef(lasso.fit, s = lasso.criterion)[-1])[(ncol(Xtilde) + 1):(ncol(Xtilde) + ncol(Z))]
      ridge.penalty.factor = rep(1, ncol(Z))
      ridge.penalty.factor[which(beta_lasso != 0)] = 0.000001
    }
    lambdahat <- chooseridge.linear(Y, Xtilde, Z, ridge.penalty.factor, lambdastart=lower, lambdaend=upper,
                                    intervals=nintervals, plot=plotGCV, file=plotfile,
                                    center.Z=F, scale.Z=F, weights.Z=NULL)
    ridgemodel <- ridge.linear(Y, Xtilde, Z, ridge.penalty.factor, lambda = lambdahat, center.Z=F, scale.Z=F, weights.Z=NULL)
    Yhat <- ridgemodel$Yhat
  }else{
    if (lasso.ols) {
      #library(glmnet)
      lasso.fit = cv.glmnet(x = cbind(Xtilde, Z), y = Y, family = c("gaussian"),
                            alpha = 1, penalty.factor = c(rep(0, ncol(Xtilde)), rep(1, ncol(Z))))
      beta_lasso = c(coef(lasso.fit, s = lasso.criterion)[-1])[(ncol(Xtilde) + 1):(ncol(Xtilde) + ncol(Z))]
      ridge.penalty.factor = rep(1, ncol(Z))
      ridge.penalty.factor[which(beta_lasso == 0)] = NaN
      upper = 0.000001
    }
    lambdahat <- upper
    ridgemodel <- ridge.linear(Y, Xtilde, Z, ridge.penalty.factor, lambda = lambdahat, center.Z=T, scale.Z=scale.Z, weights.Z=weights.Z, lasso.ols=lasso.ols)
    Yhat <- ridgemodel$Yhat
  }
  #---------------------------------------------------------------------------
  # Score statistic
  #---------------------------------------------------------------------------

  if(is.null(weights.V)==F){
    V <- t(t(V) * (weights.V))
  }

  # Q <- t(Y-Yhat) %*% V %*% t(V) %*% (Y-Yhat)
  #Q1 <- t(Y-Yhat) %*% V
  #Q <- Q1 %*% t(Q1)


  #varhat <- var(Y-Yhat)                             # Change in GxE-scoretest-v8.R
  #df1 <- sum(ridgemodel$W * t(ridgemodel$invW))      # Change in GxE-scoretest-v8.R
  #varhat <- var(Y-Yhat) * (n-1) / (n - df1)          # Change in GxE-scoretest-v8.R
  s2 = sum((Y - Yhat)^2)
  
  D0 = diag(n)
  P0 = D0 - ridgemodel$W %*% ridgemodel$invW #### I-projection matrix
  P0_square = P0 %*% P0
  K = V %*% t(V)
  PKP = P0 %*% K %*% P0
  q = as.numeric(t(Y) %*% PKP %*% Y / s2)
  A = PKP - q * P0_square

  #M1 <- t(V) - t(V) %*% ridgemodel$W %*% ridgemodel$invW     # W = ridgemodel$W = all variables under the null with appropriate scaling, centering, invW = inverse %*% t(W)
  #M2 <- M1 %*% t(M1)
  #M3 <- drop(varhat)*M2

  #---------------------------------------------------------------------------
  # p-value from non-central chi-square approximation
  #---------------------------------------------------------------------------
  if(type=="liu"){
    pvalue <- Liu.pval(0, A, Yhat)
    Is_converge <- 1
  }else{

    #---------------------------------------------------------------------------
    # p-value from davies
    #---------------------------------------------------------------------------
    daviesout <- Get_PValue_GESAT(A, Yhat)
    pvalue <- daviesout$p.value
    Is_converge <- daviesout$is_converge

    if(Is_converge<=0){
      pvalue <- Liu.pval(0, A, Yhat)
    }

  }

  return(list(pvalue=pvalue, Is_converge=Is_converge, lambda=drop(lambdahat)))
  #return(list(pvalue=pvalue, beta_lasso = beta_lasso, OLS = ridgemodel$thetahat))
}


ridge.linear <- function(Y, Xtilde, Z, ridge.penalty.factor, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL, lasso.ols=F){
  # modified in v7, intercept is included, Y is not centered
  # modified in v5 to return W, invW instead of H and GCV, matrix computations made faster,
  # add 3 arguments center.Z=T, scale.Z=T, weights.Z=NULL
  # computes the ridge estimator for each value of the ridge parameter, lambda
  # Xtilde is a n*qtilde matrix,
  # where the qtilde covariates do not have penalty imposed
  # Z is a n*p matrix where a penalty is imposed on the p covariates
  # Y is the n*1 matrix of outcomes
  # returns a (qtilde+p)*1 matrix for thetahat, n*1 matrix for yhat, (qtilde+p)*(qtilde+p) matrix for Hat matrix
  # when lambda=0, Yhat is the same as that from predict(lm(Y~Xtilde+Z))
  # NB: if scale.Z=T, regardless of what weights.Z is, this corresponds to beta(MAF; 0.5, 0.5) weight
  # to use weights specified as weights.Z, set scale.Z=F

  #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  #if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")

  n <- nrow(Y)
  if (any(is.nan(ridge.penalty.factor))) {
    Z = Z[, -(which(is.nan(ridge.penalty.factor))), drop = FALSE]
    ridge.penalty.factor = ridge.penalty.factor[-which(is.nan(ridge.penalty.factor))]
  }
  qtilde <- ncol(Xtilde) + 1 		 # +1 is for intercept, new in v7
  p <- ncol(Z)
  if(p==1){
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda * drop(ridge.penalty.factor))))
  }else{
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda, p) * ridge.penalty.factor)))
  }
  
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  W <- cbind(rep(1,n), scale(Xtilde), scale(Z, center=center.Z, scale=scale.Z)) # doesn't matter if Xtilde is scaled or not

  transW <- t(W)
  invW <- solve(transW %*% W + temp, transW) #invW <- solve(t(W) %*% W + temp, t(W))
  thetahat <- invW %*% Y

  Yhat <- W %*% thetahat
  
  return(list(lambda=lambda, thetahat=thetahat, Yhat = Yhat, W=W, invW=invW))
}


ridge.select.linear <- function(Y, Xtilde, Z, ridge.penalty.factor, lambda=0, center.Z=T, scale.Z=T, weights.Z=NULL){
  # modified in v7 to include intercept, but not center Y
  # modified in v6 to add a try() function to matrix inversion
  # modified in v5 to add 3 arguments center.Z=T scale.Z=T, weights.Z=NULL and matrix computations made faster
  # function is identical to ridge.linear(), except it only returns GCV and lambda
  # computes the ridge estimator for each value of the ridge parameter, lambda
  # Xtilde is a n*qtilde matrix,
  # where the qtilde covariates do not have penalty imposed
  # Z is a n*p matrix where a penalty is imposed on the p covariates
  # Y is the n*1 matrix of outcomes
  # when lambda=0, Yhat is the same as that from predict(lm(Y~Xtilde+Z))
  # NB: if scale.Z=T, regardless of what weights.Z is, this corresponds to beta(MAF; 0.5, 0.5) weight
  # to use weights specified as weights.Z, set scale.Z=F

  #if(nrow(Xtilde)!= nrow(Y)) stop("dimensions of Xtilde and Y don't match")
  #if(nrow(Z)!= nrow(Y)) stop("dimensions of Z and Y don't match")
  #if(scale.Z==T & is.null(weights.Z)==F) print("Warning: since scale.Z=T, weights.Z are ignored! To use weights as weights.Z, set scale.Z=F")

  n <- nrow(Y)
  qtilde <- ncol(Xtilde) + 1			 # +1 is for intercept, new in v7
  p <- ncol(Z)
  if(p==1){
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), matrix(lambda * drop(ridge.penalty.factor))))
  }else{
    temp <- rbind(matrix(0, nrow=qtilde,ncol=qtilde+p), cbind(matrix(0, nrow=p, ncol=qtilde), diag(rep(lambda, p) * ridge.penalty.factor)))
  }
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  W <- cbind(rep(1,n), scale(Xtilde), scale(Z, center=center.Z, scale=scale.Z)) # doesn't matter if Xtilde is scaled or not
  #Ys <- Y - mean(Y)

  #inverse <- solve(t(W) %*% W + temp)
  #invW <- inverse %*% t(W)
  transW <- t(W)

  # change in v6 to add a try function
  GCV_effective.df <- tryCatch({
    invW <- solve(transW %*% W + temp, transW)
    thetahat <- invW %*% Y
    Yhat <- W %*% thetahat
    equivH <- invW %*% W         					# tr(H) = tr(equivH)
    effective.df <- mtrace(equivH)-1  			 	# -1 for intercept
    GCV <- sum((Y-Yhat)^2)/(n*(1-effective.df/n)^2)
    return(list(GCV = GCV, effective.df = effective.df))
  },
  error = function(cnd) {
    return(list(GCV = NA, effective.df = NA))

  })

  return(c(list(lambda = lambda), GCV_effective.df))
}


chooseridge.linear <- function(Y, Xtilde, Z, ridge.penalty.factor, lambdastart=NULL, lambdaend=NULL, intervals=5, plot=F, file=NA, center.Z=T, scale.Z=T, weights.Z=NULL){
  # modified in v5
  # need lambdastart, lambdaend >=0, lambdastart <= lambdaend
  n <- nrow(Y)
  if(is.null(weights.Z)==F & scale.Z==F){
    Z <- t(t(Z) * (weights.Z))
  }
  sd_y <- sqrt(var(Y) * (n - 1) / n)[1, 1]
  W <- cbind(scale(Xtilde), scale(Z, center=center.Z, scale=scale.Z)) # doesn't matter if Xtilde is scaled or not
  
  if (is.null(lambdastart) == FALSE && is.null(lambdaend) == FALSE){
    lambda <- c(exp(seq(log(lambdastart),log(lambdaend),length=intervals)))
    lambda_glmnet = lambda * sd_y / n
    fit_glmnet <- cv.glmnet(W, Y, alpha=0, lambda = lambda_glmnet, penalty.factor = c(rep(0, ncol(Xtilde)), ridge.penalty.factor), standardize = F)
  } else{
    fit_glmnet <- cv.glmnet(W, Y, alpha=0, penalty.factor = c(rep(0, ncol(Xtilde)), ridge.penalty.factor), standardize = F)
  }
  lambdafinal <- fit_glmnet$lambda.min * n / sd_y

  return(lambdafinal)
}

#--------------------------------------------------------------------------------------------------------------
# END: Functions specific to linear
#
#--------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# START: common functions shared by both linear and logistic models
# NB: Functions common to both have same names in both scripts
# NB: Common functions are identical in both GxE-scoretest-logistic-snpset-v19.R and GxE-scoretest-snpset-v5.R
#-------------------------------------------------------------------------------------------------------
Get_PValue_GESAT <- function(A, Yhat, acc = 1e-6, lim=1e6){
  
  ee = eigen(A, symmetric = T)
  idx0 = which(abs(ee$values) >= 1e-10)
  lambda0 = ee$value[idx0]
  U = ee$vectors[, idx0]
  mu_norm = t(U) %*% Yhat
  noncentralparam = mu_norm ^ 2
  
  out <- davies(0, lambda0, delta = noncentralparam[,1], acc = acc, lim = lim)

  p.val <- out$Qq
  is_converge <- 1
  # check p-value
  if(p.val > 1 || p.val < 0 | out$ifault > 0){
    p.val = NA
    is_converge <- -2
  }

  # check convergence
  if(length(lambda0) == 1){
    p.val <-  NA
    is_converge <- -1
  } else if(out$ifault != 0){
    p.val <- NA
    is_converge <- 0
  }


  return(list(p.value=p.val, is_converge=is_converge))

}

Liu.pval = function(Q, A, Yhat){
  AA <- A %*% A
  AAA = A %*% AA
  AAAA = A %*% AAA
  kappa1 <-  mtrace(A) + mtrace(t(Yhat) %*% A %*% Yhat)	          # = tr(M3)
  kappa2 <-  2 * (mtrace(AA) + mtrace(t(Yhat) %*% AA %*% Yhat))           # = 2 tr(M3 M3) = 2 tr(M4)
  kappa3 <-  8 * (mtrace(AAA) + mtrace(t(Yhat) %*% AAA %*% Yhat))      # = 8 tr(M3 M3 M3) = 8 tr(M3 M4) = 8 sum(M3 * M4')
  kappa4 <-  48 * (mtrace(AAAA) + mtrace(t(Yhat) %*% AAAA %*% Yhat))     # = 48 tr(M3 M3 M3 M3) = 48 tr(M4 M4) = 48 sum(M4 * M4')
  
  approx2 <- noncentralapproxdirect(kappa2, kappa3, kappa4)
  Q.Norm2 <-((Q - kappa1)/sqrt(kappa2))*approx2$sigmaX +  approx2$muX
  pvalue <- pchisq(Q.Norm2, df=approx2$df, ncp = approx2$ncp,  lower.tail=F)
  return(pvalue)
}



qqplots <- function(pvalues, header = "Quantile-Quantile Plot of P-values", filename = "qqplot.png"){
  temp <- sort(pvalues)
  N <- length(pvalues)
  observed.quantiles <- -log10(temp)
  expected.quantiles <- -log10((1:N)/(N+1))
  png(file = filename, width = 800, height = 800)
  par(mar=c(4.5,4.5,3,1.5)+0.0)
  plot(expected.quantiles, observed.quantiles, col = "red2", pch = 19,
       xlab = "-log(Expected P-value)", ylab = "-log(Observed P-value)", main = header,
       cex.main = 1, family = "serif")
  abline(0,1)
  dev.off()
}


qqplotspdf <- function(pvalues, header = "Quantile-Quantile Plot of P-values", filename = "qqplot.pdf"){
  temp <- sort(pvalues)
  N <- length(pvalues)
  observed.quantiles <- -log10(temp)
  expected.quantiles <- -log10((1:N)/(N+1))
  pdf(file = filename, width = 6, height = 6)
  par(mar=c(4.5,4.5,3,1.5)+0.0)
  plot(expected.quantiles, observed.quantiles, col = "red2", pch = 19,
       xlab = "-log(Expected P-value)", ylab = "-log(Observed P-value)", main = header,
       cex.main = 1, family = "serif")
  abline(0,1)
  dev.off()
}


qqplotspdfNA <- function(inputpvalues, header = "Quantile-Quantile Plot of P-values", filename = "qqplot.pdf"){
  # similar to qqplotspdf except it first removes NA values
  pvalues <- inputpvalues[is.na(inputpvalues)==F]
  temp <- sort(pvalues)
  N <- length(pvalues)
  observed.quantiles <- -log10(temp)
  expected.quantiles <- -log10((1:N)/(N+1))
  pdf(file = filename, width = 6, height = 6)
  par(mar=c(4.5,4.5,3,1.5)+0.0)
  plot(expected.quantiles, observed.quantiles, col = "red2", pch = 19,
       xlab = "-log(Expected P-value)", ylab = "-log(Observed P-value)", main = header,
       cex.main = 1, family = "serif")
  abline(0,1)
  dev.off()
}


skewness <-  function(x){
  return((mean((x-mean(x))^3))/(sd(x)^3))
}


kurtosis <- function(x){
  return((mean((x-mean(x))^4))/(sd(x)^4))
}


noncentralapproxdirect <- function(k2, k3, k4){
  # k2, k3, k4 are the 2nd, 3rd and 4th cumulants

  s1 <- k3/((k2^1.5)*sqrt(8))
  s2 <- k4/(k2*k2*12)

  if(s1^2>s2){
    a <- 1/(s1-sqrt(s1^2-s2))
    delta <- s1*a^3-a^2
    l <- a^2-2*delta
  }else{
    a <- 1/s1
    delta <- 0
    l <- 1/s1^2
  }
  return(list(df=l, ncp=delta, muX=l+delta, sigmaX=sqrt(2)*a))
}


mtrace <- function(X){
  return(sum(diag(X)))
}


mafcall <- function(x){
  return(min(sum(x, na.rm=T), 2*sum(is.na(x)==F)-sum(x, na.rm=T))/(2*sum(is.na(x)==F)))
}


checkpolymorphic <- function(x){
  return(length(unique(x))>1)
}
#-------------------------------------------------------------------------------------------------------
# END: common functions shared by both linear and logistic models
#
#-------------------------------------------------------------------------------------------------------
