// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// col_normal_em
List col_normal_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP rselection);
RcppExport SEXP oda_col_normal_em(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP, SEXP rselectionSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rselection(rselectionSEXP );
        List __result = col_normal_em(ryo, rxo, rxa, rd, rlam, rpriorprob, rselection);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// col_normal_em_soft
List col_normal_em_soft(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob);
RcppExport SEXP oda_col_normal_em_soft(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        List __result = col_normal_em_soft(ryo, rxo, rxa, rlam, rpriorprob);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// col_normal_gibbs
List col_normal_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter);
RcppExport SEXP oda_col_normal_gibbs(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP, SEXP rburninSEXP, SEXP rniterSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rburnin(rburninSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rniter(rniterSEXP );
        List __result = col_normal_gibbs(ryo, rxo, rxa, rlam, rpriorprob, rburnin, rniter);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// col_normal_var
List col_normal_var(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob);
RcppExport SEXP oda_col_normal_var(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        List __result = col_normal_var(ryo, rxo, rxa, rlam, rpriorprob);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// col_mixture_gibbs
List col_mixture_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rpriorprob, SEXP rburnin, SEXP rniter, SEXP ralpha);
RcppExport SEXP oda_col_mixture_gibbs(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rpriorprobSEXP, SEXP rburninSEXP, SEXP rniterSEXP, SEXP ralphaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rburnin(rburninSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rniter(rniterSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ralpha(ralphaSEXP );
        List __result = col_mixture_gibbs(ryo, rxo, rxa, rpriorprob, rburnin, rniter, ralpha);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// lasso_em
List lasso_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, SEXP rlasso);
RcppExport SEXP oda_lasso_em(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlassoSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rlasso(rlassoSEXP );
        List __result = lasso_em(ryo, rxo, rxa, rd, rlasso);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mixture_var
List mixture_var(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP ralpha, SEXP rrho_min, SEXP ranneal_steps);
RcppExport SEXP oda_mixture_var(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP, SEXP ralphaSEXP, SEXP rrho_minSEXP, SEXP ranneal_stepsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ralpha(ralphaSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rrho_min(rrho_minSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ranneal_steps(ranneal_stepsSEXP );
        List __result = mixture_var(ryo, rxo, rxa, rd, rlam, rpriorprob, ralpha, rrho_min, ranneal_steps);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// normal_em
List normal_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob);
RcppExport SEXP oda_normal_em(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        List __result = normal_em(ryo, rxo, rxa, rd, rlam, rpriorprob);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// normal_gibbs
List normal_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter);
RcppExport SEXP oda_normal_gibbs(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP, SEXP rburninSEXP, SEXP rniterSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rburnin(rburninSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rniter(rniterSEXP );
        List __result = normal_gibbs(ryo, rxo, rxa, rd, rlam, rpriorprob, rburnin, rniter);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// normal_var
List normal_var(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP rrho_min, SEXP ranneal_steps);
RcppExport SEXP oda_normal_var(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP, SEXP rrho_minSEXP, SEXP ranneal_stepsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rrho_min(rrho_minSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ranneal_steps(ranneal_stepsSEXP );
        List __result = normal_var(ryo, rxo, rxa, rd, rlam, rpriorprob, rrho_min, ranneal_steps);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// oda_gdp_em
List oda_gdp_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rpriorprob, SEXP ralpha, SEXP reta);
RcppExport SEXP oda_oda_gdp_em(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rpriorprobSEXP, SEXP ralphaSEXP, SEXP retaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ralpha(ralphaSEXP );
        Rcpp::traits::input_parameter< SEXP >::type reta(retaSEXP );
        List __result = oda_gdp_em(ryo, rxo, rxa, rd, rpriorprob, ralpha, reta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// or_gdp_em
List or_gdp_em(NumericVector ryo, NumericMatrix rxo, SEXP ralpha, SEXP reta);
RcppExport SEXP oda_or_gdp_em(SEXP ryoSEXP, SEXP rxoSEXP, SEXP ralphaSEXP, SEXP retaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ralpha(ralphaSEXP );
        Rcpp::traits::input_parameter< SEXP >::type reta(retaSEXP );
        List __result = or_gdp_em(ryo, rxo, ralpha, reta);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// or_lasso_em
List or_lasso_em(NumericVector ryo, NumericMatrix rxo, SEXP rlasso);
RcppExport SEXP oda_or_lasso_em(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rlassoSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rlasso(rlassoSEXP );
        List __result = or_lasso_em(ryo, rxo, rlasso);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// probit_normal_em
List probit_normal_em(NumericVector rbit, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob);
RcppExport SEXP oda_probit_normal_em(SEXP rbitSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type rbit(rbitSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        List __result = probit_normal_em(rbit, rxo, rxa, rd, rlam, rpriorprob);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// probit_normal_gibbs
List probit_normal_gibbs(NumericVector rbit, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter);
RcppExport SEXP oda_probit_normal_gibbs(SEXP rbitSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rlamSEXP, SEXP rpriorprobSEXP, SEXP rburninSEXP, SEXP rniterSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type rbit(rbitSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rlam(rlamSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rburnin(rburninSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rniter(rniterSEXP );
        List __result = probit_normal_gibbs(rbit, rxo, rxa, rd, rlam, rpriorprob, rburnin, rniter);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP oda_rcpparma_hello_world() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        arma::mat __result = rcpparma_hello_world();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP oda_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP );
        arma::mat __result = rcpparma_outerproduct(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP oda_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP );
        double __result = rcpparma_innerproduct(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP oda_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP );
        Rcpp::List __result = rcpparma_bothproducts(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// t_em
List t_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rpriorprob, SEXP ralpha);
RcppExport SEXP oda_t_em(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rpriorprobSEXP, SEXP ralphaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ralpha(ralphaSEXP );
        List __result = t_em(ryo, rxo, rxa, rd, rpriorprob, ralpha);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// mixture_gibbs
List mixture_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rpriorprob, SEXP rburnin, SEXP rniter, SEXP ralpha);
RcppExport SEXP oda_mixture_gibbs(SEXP ryoSEXP, SEXP rxoSEXP, SEXP rxaSEXP, SEXP rdSEXP, SEXP rpriorprobSEXP, SEXP rburninSEXP, SEXP rniterSEXP, SEXP ralphaSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ryo(ryoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxo(rxoSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type rxa(rxaSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rd(rdSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rpriorprob(rpriorprobSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rburnin(rburninSEXP );
        Rcpp::traits::input_parameter< SEXP >::type rniter(rniterSEXP );
        Rcpp::traits::input_parameter< SEXP >::type ralpha(ralphaSEXP );
        List __result = mixture_gibbs(ryo, rxo, rxa, rd, rpriorprob, rburnin, rniter, ralpha);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
