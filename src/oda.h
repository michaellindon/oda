#ifndef ODA_H
#define ODA_H

#include <Rmath.h>
#include <R.h>
#include <vector>
#include <numeric>

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

extern "C" {
	void dpotrf_(char *  UPLO,int * N,double * A,int * LDA,int * INFO );
	void dpotrs_(char * UPLO, int * N, int * NRHS, double * A, int * LDA, double * B, int * LDB, int *  INFO);
	void dtrmv_(char * UPLO, char * TRANS, char * DIAG, int * N,double *  A,int * LDA,double * X,int * INCX);
	void dtrsv_(char * UPLO, char * TRANS, char * DIAG, int * N, double * A, int * LDA, double * X, int * INCX);
	void dgemv_(char * TRANS, int * M, int * N, double * ALPHA, double * A, int * LDA, double * X, int * INCX, double * BETA, double * Y, int * INCY);
	void daxpy_(int * N, double * DA, double * DX, int * INCX, double * DY, int * INCY);
	double ddot_(int * N, double * DX, int * INCX, double * DY, int * INCY);
}


void submatrices_collapsed(std::vector<double> &mu, std::vector<double> &xag, std::vector<double> &xogxog_Lamg, std::vector<double> &xogyo, std::vector<double> &lamg, std::vector<double> &Bg, const std::vector<int> &gamma, const std::vector<double> &xoyo, const std::vector<double> &lam, const std::vector<double> &xa, const std::vector<double> &xoxo,  int p_gamma, double &b, const double yoyo, int p, int na);

void fixed_probabilities(std::vector<double> &prob, std::vector<double> &odds, std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &xoyo, const std::vector<double> &xaya, const std::vector<double> &priorodds, const std::vector<double> &ldl, const std::vector<double> &dli, double phi);

void draw_collapsed_xaya(std::vector<double> &xaya, std::vector<double> &xa, std::vector<double> &xag, std::vector<double> &mu, double phi, std::vector<double> &Z, std::vector<double> &xogxog_Lamg, int na, int p, int p_gamma);

void draw_gamma(std::vector<int> &gamma, const std::vector<double> prob);

bool gamma_change(const std::vector<int> &gamma_mcmc, int t, int p);

#endif
