#include <numeric>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


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
	double ddot_(int * N, double * DX, int * INCX, double * DY, int * INCY);
}

inline void fixed_probabilities(Col<double> &prob, Col<double> &odds, Col<double> &Bols, const Col<double> &d, const Col<double> &xoyo, const Col<double> &xaya, const Col<double> &priorodds, const Col<double> &ldl, const Col<double> &dli, double phi){
	for(int i=0; i<prob.n_elem; ++i)
	{
		Bols(i)=(1/d(i))*(xoyo(i)+xaya(i));
		odds(i)=priorodds(i)*ldl(i)*exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i));
		prob(i)=odds(i)/(1+odds(i));
		if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN exponential
	}
};

inline void draw_collapsed_xaya(Col<double> &xaya, Mat<double> &xa, Mat<double> &xag, Col<double> &mu, double phi,Col<double> &Z, Mat<double> &xogxog_Lamg, int na, int p_gamma){

	double sd=sqrt(1/phi);
	Z.set_size(na);
	for(Col<double>::iterator it=Z.begin(); it!=Z.end(); ++it) *it=R::rnorm(0,1);
	if(p_gamma!=0){
		xaya=mu+sqrt(1/phi)*Z;
		Z.set_size(p_gamma);
		for(Col<double>::iterator it=Z.begin(); it!=Z.end(); ++it) *it=R::rnorm(0,1);
		//Computes R^{-1}Z where xogxog_Lamg^{-1}=R^{-1}R^{-T}
		dtrsv_(&uplo, &transN, &unit_tri, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &*Z.begin(), &inc);
		dgemv_(&transN , &na, &p_gamma, &sd, &*xag.begin(), &na, &*Z.begin(), &inc, &inputscale1, &*xaya.begin(), &inc);
		dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
	}else{
		xaya=sqrt(1/phi)*Z;
		dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
	}
};

inline void draw_gamma(Col<uword> &gamma, Col<double> prob, Col<double> &U){
	for(Col<double>::iterator it=U.begin(); it!=U.end(); ++it) *it=R::runif(0,1);
	for (int i = 0; i < prob.n_elem ; ++i)
	{
		if(U(i)<prob(i)){
			gamma(i)=1;
		}else{
			gamma(i)=0;
		}
	}
}

inline bool gamma_change(const Mat<uword> &gamma_mcmc, int t){

	if(sum(abs(gamma_mcmc.col(t)-gamma_mcmc.col(t-1)))==0)
	{
		return(false);
	}else
	{
		return(true);
	}
}

// [[Rcpp::export]]
List col_normal_gibbs(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob, SEXP rburnin, SEXP rniter){

	//Dimensions//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();

	//Read RData Into Matrix Classes//
	Mat<double> xo(rxo.begin(), no, p, false);
	Mat<double> xa(rxa.begin(), na, p, false);
	Col<double> yo(ryo.begin(), ryo.size(), false); yo-=mean(yo);
	Col<double> lam(rlam.begin(),rlam.size(), false);
	Col<double> priorprob(rpriorprob.begin(),rpriorprob.size(), false);
	int niter=Rcpp::as<int >(rniter);
	int burnin=Rcpp::as<int >(rburnin);

	//Create Matrices//
	Col<double> xoyo=xo.t()*yo;
	Col<double> xogyo=xoyo;
	Mat<double> xoxo=xo.t()*xo;
	Mat<double> xaxa=xa.t()*xa;
	Mat<double> xogxog_Lamg(p,p);
	Col<double> d(p); for(int i=0; i<p; ++i) d(i)=xoxo(i,i)+xaxa(i,i);
	Mat<double> xag=xa;
	Col<double> lamg=lam; //vector instead of diagonal pxp matrix

	//Phi Variables//
	double a=(double)0.5*(no-1);
	double b=1;
	double phi=1;
	double yoyo=ddot_(&no, &*yo.begin(), &inc, &*yo.begin(), &inc);

	//Ya Variables//
	Col<double> ya(na);
	Col<double> mu(na);
	Col<double> Z(na);
	Col<double> xaya(p);

	//Beta Variables//
	Col<double> Bg(p,fill::zeros);

	//Gamma Variables//
	Col<uword> gamma(p,fill::zeros);
	Col<uword> inc_indices(p,fill::ones);
	Col<double> U(p);
	Col<double> Bols(p,fill::zeros);
	Col<double> prob(p,fill::ones);
	Col<double> odds(p);
	Col<double> priorodds=priorprob/(1-priorprob);
	Col<double> ldl=sqrt(lam/(d+lam));
	Col<double> dli=1/(d+lam);
	bool gamma_diff=true;
	int p_gamma;

	//Allocate Space for MCMC Draws//
	Mat<double> ya_mcmc(na,niter);
	Mat<double> prob_mcmc(p,niter);
	Mat<uword>  gamma_mcmc(p,niter);
	Col<double> phi_mcmc(niter);
	ya_mcmc.col(0)=ya;
	phi_mcmc(0)=phi;
	gamma_mcmc.col(0)=gamma;
	prob_mcmc.col(0)=prob;


	//Begin Gibbs Sampling Algorithm//
	for (int t = 1; t < niter; ++t)
	{
		//If gamma[t]!=Gamma[t-1] Form Submatrices//
		if(gamma_diff)
		{
			p_gamma=sum(gamma);
			if(p_gamma!=0){
				inc_indices=find(gamma);
				xag=xa.cols(inc_indices);
				xogyo=xoyo.elem(inc_indices);
				lamg=lam.elem(inc_indices);
				xogxog_Lamg=xoxo.submat(inc_indices,inc_indices);
				for(int d=0; d<p_gamma; ++d)	xogxog_Lamg(d,d)+=lamg(d); 
				//Positive Definite Cholesky Factorization//
				dpotrf_(&uplo, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &info); //xerbla handles info error
				//xogxog_Lamg now stores R upper triangular where xogxog_Lamg=R'R
				Bg=xogyo; 
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
		draw_collapsed_xaya(xaya,xa,xag,mu,phi,Z,xogxog_Lamg,na,p_gamma);

		//Draw Gamma//
		fixed_probabilities(prob, odds, Bols, d, xoyo, xaya, priorodds, ldl, dli, phi);
		draw_gamma(gamma,prob,U);

		//Store Draws//
		gamma_mcmc.col(t)=gamma;
		prob_mcmc.col(t)=prob;
		ya_mcmc.col(t)=ya;
		phi_mcmc(t)=phi;

		//Has Gamma Changed?//
		gamma_diff=gamma_change(gamma_mcmc,t);

	}


	//Return Data to R//
	return Rcpp::List::create(
			Rcpp::Named("phi_mcmc") = phi_mcmc,
			Rcpp::Named("prob_mcmc") = prob_mcmc,
			Rcpp::Named("prob") = mean(prob_mcmc.cols(burnin-1,niter-1),1),
			Rcpp::Named("gamma_mcmc") = gamma_mcmc,
			Rcpp::Named("ya_mcmc") = ya_mcmc
			) ;
}


//Notes: Lam is O(pxp) space complexity, terrible, just use lam O(p) instead.
//Better to pass xoxo and d instead of xoxo & xaxa and computing d inside the code, as this requires O(pxp) instead of O(p) [xa is needed anyway] and xa.t()*xa takes forever
//dtritri otherwise need to do dpotrs and then dpotrf again
