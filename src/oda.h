#ifndef ODA_H
#define ODA_H

#include <Rmath.h>
#include <R.h>
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>

//LAPACK Variables//
extern char uplo;
extern char transN;
extern char transT;
extern char unit_tri;
extern int info;
extern int nrhs;
extern int inc;
extern double inputscale1;
extern double inputscale0;
extern double unity;
extern double nunity;

extern "C" {
	void dpotrf_(char *  UPLO,int * N,double * A,int * LDA,int * INFO );
	void dpotrs_(char * UPLO, int * N, int * NRHS, double * A, int * LDA, double * B, int * LDB, int *  INFO);
	void dtrmv_(char * UPLO, char * TRANS, char * DIAG, int * N,double *  A,int * LDA,double * X,int * INCX);
	void dtrsv_(char * UPLO, char * TRANS, char * DIAG, int * N, double * A, int * LDA, double * X, int * INCX);
	void dgemv_(char * TRANS, int * M, int * N, double * ALPHA, double * A, int * LDA, double * X, int * INCX, double * BETA, double * Y, int * INCY);
	void daxpy_(int * N, double * DA, double * DX, int * INCX, double * DY, int * INCY);
	double ddot_(int * N, double * DX, int * INCX, double * DY, int * INCY);
}



void submatrices_uncollapsed(bool gamma_diff, const std::vector<double> B, std::vector<double> &xog, std::vector<double> &xag,  std::vector<double> &lamg, std::vector<double> &Bg, const std::vector<int> &gamma,  const std::vector<double> &lam, const std::vector<double> &xo, const std::vector<double> &xa,  int p_gamma, double &b, int p, int no, int na);

void submatrices_collapsed(bool gamma_diff, std::vector<double> &mu, std::vector<double> &xag, std::vector<double> &xogxog_Lamg, std::vector<double> &xogyo, std::vector<double> &lamg, std::vector<double> &Bg, const std::vector<int> &gamma, const std::vector<double> &xoyo, const std::vector<double> &lam, const std::vector<double> &xa, const std::vector<double> &xoxo,  int p_gamma, double &b, const double yoyo, int p, int na);

void bernoulli_probabilities(std::vector<double> &prob, std::vector<double> &odds, std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &xoyo, const std::vector<double> &xaya, const std::vector<double> &priorprob, const std::vector<double> &lam, const double phi);

void betabinomial_probabilities(std::vector<double> &prob, std::vector<double> &odds, std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &xoyo, const std::vector<double> &xaya, const double theta, const std::vector<double> &lam, const double phi);

void uniform_probabilities(std::vector<double> &prob, std::vector<double> &odds, std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &xoyo, const std::vector<double> &xaya,  const std::vector<double> &lam, const double phi);

void draw_collapsed_xaya(std::vector<double> &xaya, std::vector<double> &xa, std::vector<double> &xag, std::vector<double> &mu, double phi, std::vector<double> &Z, std::vector<double> &xogxog_Lamg, int na, int p, int p_gamma);

void draw_uncollapsed_xaya(std::vector<double> &xaya, std::vector<double> &xa, std::vector<double> &xag, std::vector<double> Bg,  double phi, std::vector<double> &Z, int na, int p, int p_gamma);

double draw_collapsed_phi( double &b, int p_gamma, const int no, const double yoyo, std::vector<double> &xogyo, std::vector<double> &Bg);

double draw_uncollapsed_phi(int p_gamma, int no, const std::vector<double> &yo, std::vector<double> &xog, std::vector<double> &Bg, const std::vector<double> &lamg, const double yoyo);

void draw_gamma(std::vector<int> &gamma, int &p_gamma, const std::vector<double> prob);

void draw_beta(const std::vector<int> &gamma, std::vector<double> &B, const std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &lam, double phi);

void draw_lambda_t( std::vector<double> &lam, const std::vector<int> &gamma, const double alpha, const std::vector<double> &B, const double phi);

bool gamma_change(const int * gamma_mcmc, int t, int p);

void rao_blackwell(double * B_rb, double * prob_rb, const std::vector<double> B, const std::vector<double> prob, int burnin, int niter);
#endif
