#include <vector>
#include <algorithm>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


//LAPACK Variables//
char uplo = 'U';
char transN = 'N';
char transT = 'T';
char unit_tri = 'N';
int info;
int nrhs=1;
int inc=1;
double inputscale1=1.0;
double inputscale0=0.0;
double unity=1.0;


extern "C" {
	void dpotrf_(char *  UPLO,int * N,double * A,int * LDA,int * INFO );
	void dpotrs_(char * UPLO, int * N, int * NRHS, double * A, int * LDA, double * B, int * LDB, int *  INFO);
	void dtrmv_(char * UPLO, char * TRANS, char * DIAG, int * N,double *  A,int * LDA,double * X,int * INCX);
	void dtrsv_(char * UPLO, char * TRANS, char * DIAG, int * N, double * A, int * LDA, double * X, int * INCX);
	void dgemv_(char * TRANS, int * M, int * N, double * ALPHA, double * A, int * LDA, double * X, int * INCX, double * BETA, double * Y, int * INCY);
	void daxpy_(int * N, double * DA, double * DX, int * INCX, double * DY, int * INCY);
	double ddot_(int * N, double * DX, int * INCX, double * DY, int * INCY);
}

inline void fixed_probabilities(vector<double> &prob, vector<double> &odds, vector<double> &Bols, const vector<double> &d, const vector<double> &xoyo, const vector<double> &xaya, const vector<double> &priorodds, const vector<double> &ldl, const vector<double> &dli, double phi){
	for(size_t i=0; i!=prob.size(); ++i)
	{
		Bols[i]=(1/d[i])*(xoyo[i]+xaya[i]);
		odds[i]=priorodds[i]*ldl[i]*exp(0.5*phi*dli[i]*d[i]*d[i]*Bols[i]*Bols[i]);
		prob[i]=odds[i]/(1+odds[i]);
		if(prob[i]!=prob[i]) prob[i]=1;	 //Catch NaN exponential
	}
};

inline void draw_collapsed_xaya(vector<double> &xaya, vector<double> &xa, vector<double> &xag, vector<double> &mu, double phi, vector<double> &Z, vector<double> &xogxog_Lamg, int na, int p, int p_gamma){

	double sd=sqrt(1/phi);
	Z.resize(na);
	for(vector<double>::iterator it=Z.begin(); it!=Z.end(); ++it) *it=R::rnorm(0,1);
	if(p_gamma!=0){
		xaya=mu;
		daxpy_(&p, &sd, &*Z.begin(), &inc, &*xaya.begin(), &inc);
		Z.resize(p_gamma);
		for(vector<double>::iterator it=Z.begin(); it!=Z.end(); ++it) *it=R::rnorm(0,1);
		//Computes R^{-1}Z where xogxog_Lamg^{-1}=R^{-1}R^{-T}
		dtrsv_(&uplo, &transN, &unit_tri, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &*Z.begin(), &inc);
		dgemv_(&transN , &na, &p_gamma, &sd, &*xag.begin(), &na, &*Z.begin(), &inc, &inputscale1, &*xaya.begin(), &inc);
		dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
	}else{
		for(size_t i=0; i!=xaya.size(); ++i) xaya[i]=sd*Z[i];
		dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
	}
};

inline void draw_gamma(vector<int> &gamma, vector<double> prob, vector<double> &U){
	for(vector<double>::iterator it=U.begin(); it!=U.end(); ++it) *it=R::runif(0,1);
	for(size_t i = 0; i!=prob.size() ; ++i)
	{
		if(U[i]<prob[i]){
			gamma[i]=1;
		}else{
			gamma[i]=0;
		}
	}
}

inline bool gamma_change(const vector<int> &gamma_mcmc, int t, int p){

	for(int i=0; i<p; ++i)	if(gamma_mcmc[t*p+i]!=gamma_mcmc[(t-1)*p+i]) return(true);
	return(false);
}

// [[Rcpp::export]]
List col_normal_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericMatrix rxoxo, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter){

	//Dimensions//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();

	//Read RData Into Matrix Classes//
	vector<double> xo(rxo.begin(), rxo.begin()+no*p);
	vector<double> xa(rxa.begin(), rxa.begin()+na*p);
	vector<double> yo(ryo.begin(), ryo.end()); 
	double yobar=std::accumulate(yo.begin(),yo.end(),0.0)/yo.size();
	for(vector<double>::iterator it=yo.begin(); it!=yo.end(); ++it) *it-=yobar;
	vector<double> lam(rlam.begin(),rlam.end());
	vector<double> d(rd.begin(),rd.end());
	vector<double> xoxo(rxoxo.begin(),rxoxo.end());
	vector<double> priorprob(rpriorprob.begin(),rpriorprob.end());
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);

	//Create Matrices//
	vector<double> xoyo(p);
	dgemv_(&transT , &no, &p, &unity, &*xo.begin(), &no, &*yo.begin(), &inc, &inputscale0, &*xoyo.begin(), &inc);
	vector<double> xogyo; xogyo.reserve(p);
	vector<double> xogxog_Lamg; xogxog_Lamg.reserve(p*p);
	vector<double> xag; xag.reserve(na*p);
	vector<double> lamg; lamg.reserve(p); //vector instead of diagonal pxp matrix

	//Phi Variables//
	double a=(double)0.5*(no-1);
	double b=1;
	double phi=1;
	double yoyo=ddot_(&no, &*yo.begin(), &inc, &*yo.begin(), &inc);

	//Ya Variables//
	vector<double> mu(na);
	vector<double> Z; Z.reserve(na);
	vector<double> xaya(p);

	//Beta Variables//
	vector<double> Bg; Bg.reserve(p);

	//Gamma Variables//
	vector<int> gamma(p);
	vector<double> U(p);
	vector<double> Bols(p);
	vector<double> prob(p);
	vector<double> odds(p);
	vector<double> priorodds(p);
	vector<double> ldl(p);
	vector<double> dli(p);
	for(size_t i=0; i!=priorodds.size(); ++i){
		priorodds[i]=priorprob[i]/(1-priorprob[i]);
		ldl[i]=sqrt(lam[i]/(d[i]+lam[i]));
		dli[i]=1/(d[i]+lam[i]);
	}
	bool gamma_diff=true;
	int p_gamma;

	//Allocate Space for MCMC Draws//
	vector<double> prob_mcmc(p*niter);
	vector<int>  gamma_mcmc(p*niter);
	vector<double> phi_mcmc(niter);
	std::copy(prob.begin(),prob.end(),prob_mcmc.begin());
	std::copy(gamma.begin(),gamma.end(),gamma_mcmc.begin());
	phi_mcmc[0]=phi;

	//Begin Gibbs Sampling Algorithm//
	for (int t = 1; t < niter; ++t)
	{
		//If gamma[t]!=Gamma[t-1] Form Submatrices//
		if(gamma_diff)
		{
			p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);
			if(p_gamma!=0){
				xag.resize(0);
				xogxog_Lamg.resize(0);
				xogyo.resize(0);
				lamg.resize(0);
				Bg.resize(0);
				for(int i=0; i<p; ++i)
				{
					if(gamma[i]==1)
					{
						xogyo.push_back(xoyo[i]);
						Bg.push_back(xoyo[i]);
						lamg.push_back(lam[i]);
						for(int j=0; j<na; ++j) xag.push_back(xa[i*na+j]);
						for(int j=0; j<p; ++j) if(gamma[j]==1) xogxog_Lamg.push_back(xoxo[i*p+j]);
					}
				}
				for(int i=0; i<p_gamma; ++i) xogxog_Lamg[i*p_gamma+i]+=lamg[i];
				//Positive Definite Cholesky Factorization//
				dpotrf_(&uplo, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &info); //xerbla handles info error
				//xogxog_Lamg now stores R upper triangular where xogxog_Lamg=R'R
				//Triangular Positive Definite Solve via Cholesky
				dpotrs_(&uplo,  &p_gamma, &nrhs, &*xogxog_Lamg.begin(), &p_gamma, &*Bg.begin(),  &p_gamma, &info);

				//b=0.5*(yoyo-std::inner_product(xogyo.begin(),xogyo.end(),Bg.begin(),0));
				b=0.5*(yoyo-ddot_(&p_gamma, &*xogyo.begin(), &inc, &*Bg.begin(), &inc));
				dgemv_(&transN , &na, &p_gamma, &unity, &*xag.begin(), &na, &*Bg.begin(), &inc, &inputscale0, &*mu.begin(), &inc);
			}else{
				b=0.5*yoyo;
			}
		}

		//Draw Phi//
		phi=R::rgamma(a,(1/b)); //rgamma uses scale

		//Draw Ya//
		draw_collapsed_xaya(xaya,xa,xag,mu,phi,Z,xogxog_Lamg,na,p,p_gamma);

		//Draw Gamma//
		fixed_probabilities(prob, odds, Bols, d, xoyo, xaya, priorodds, ldl, dli, phi);
		draw_gamma(gamma,prob,U);

		//Store Draws//
		std::copy(gamma.begin(),gamma.end(),gamma_mcmc.begin()+p*t);
		std::copy(prob.begin(),prob.end(),prob_mcmc.begin()+p*t);
		phi_mcmc[t]=phi;

		//Has Gamma Changed?//
		gamma_diff=gamma_change(gamma_mcmc,t,p);

	}


	//Return Data to R//
	return Rcpp::List::create(
			Rcpp::Named("phi_mcmc") = phi_mcmc,
			Rcpp::Named("prob_mcmc") = prob_mcmc,
			Rcpp::Named("gamma_mcmc") = gamma_mcmc
			) ;
}


//Notes: Lam is O(pxp) space complexity, terrible, just use lam O(p) instead.
//Better to pass xoxo and d instead of xoxo & xaxa and computing d inside the code, as this requires O(pxp) instead of O(p) [xa is needed anyway] and xa.t()*xa takes forever
//dtritri otherwise need to do dpotrs and then dpotrf again
