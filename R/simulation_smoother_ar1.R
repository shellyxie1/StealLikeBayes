#' @title Samples draws from a multivariate normal distribution using the 
#' precision sampler by 
#' 
#' @description Samples random numbers from an \eqn{N}-variate normal distribution
#' specified by the precision matrix \eqn{P} and location vector \eqn{L} as in:
#' \deqn{N(P^{-1}L, P^{-1})}
#' where the precision matrix \eqn{P} is bi-diagonal with the diagonal elements 
#' given in the vector argument \code{precision_diag} and the off-diagonal element 
#' is given in the scalar argument \code{precision_offdiag}, and the location 
#' vector \eqn{L} is provided in the vector argument \code{location}.
#' 
#' @details This function is based on C++ code from the R package \pkg{stochvol}
#' by Hosszejni D., Kastner G. (2025) and Kastner G. (2016).
#' 
#' @param precision_diag an \eqn{N}-vector with the diagonal elements of the 
#' precision matrix \eqn{P}
#' @param precision_offdiag a numeric scalar with the off-diagonal element of 
#' the precision matrix \eqn{P}
#' @param location an \eqn{N}-vector with the location parameter \eqn{L}
#' @return  an \eqn{N}-vector with random draws from the multivariate normal 
#' distribution
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' 
#' Chan J.C.C., Jeliazkov I. (2009). Efficient simulation and integrated 
#' likelihood estimation in state space models. International Journal of 
#' Mathematical Modelling and Numerical Optimisation, 
#' 1(1/2), <doi:10.1504/IJMMNO.2009.030090>.
#' 
#' Kastner G. (2016). Dealing with Stochastic Volatility in Time Series Using 
#' the R Package stochvol. Journal of Statistical Software, 69(5), 1–30. 
#' <doi:10.18637/jss.v069.i05>.
#' 
#' Hosszejni D., Kastner G. (2025). stochvol: Efficient Bayesian Inference for 
#' Stochastic Volatility (SV) Models. R package version 3.2.8, 
#' <doi:10.32614/CRAN.package.stochvol>
#' 
#' @export
precision_sampler_ar1 <- function(precision_diag, precision_offdiag, location) {
  
  stopifnot(
    "The argument precision_diag must be a numeric vector with real numbers." =
    is.numeric(precision_diag) & all(!is.na(precision_diag)) 
  )
  stopifnot(
    "The argument precision_offdiag must be a real number." =
    is.numeric(precision_offdiag) & length(precision_offdiag) == 1
  )
  stopifnot(
    "The argument location must be a numeric vector with real numbers." =
    is.numeric(location) & all(!is.na(location)) 
  )
  stopifnot(
    "Arguments precision_diag and location must be of the same length." =
    length(location) == length(precision_diag)  
  )

  out = .Call(`_StealLikeBayes_precision_sampler_ar1`, precision_diag, precision_offdiag, location)
  
  return(out)
}

