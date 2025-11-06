#ifndef _SIMULATION_SMOOTHER_AR1_
#define _SIMULATION_SMOOTHER_AR1_

#include <RcppArmadillo.h>

Rcpp::List cholesky_tridiagonal(
    const arma::vec&    omega_diag,
    const double&       omega_offdiag
);


arma::vec forward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector
);


arma::vec backward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp
);


arma::vec precision_sampler_ar1(
    const arma::vec&     precision_diag,
    const double&        precision_offdiag,
    const arma::vec&     location
);

#endif