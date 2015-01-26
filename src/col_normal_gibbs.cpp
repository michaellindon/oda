#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

extern "C" {
	void dpotrf_(char *  UPLO,int * N,double * A,int * LDA,int * INFO );
	void dpotrs_(char * UPLO, int * N, int * NRHS, double * A, int * LDA, double * B, int * LDB, int *  INFO );
	void dtrmv_(char * UPLO, char * TRANS, char * DIAG, int * N,double *  A,int * LDA,double * X,int * INCX);
	void dtrsv_(char * UPLO, char * TRANS, char * DIAG, int * N, double * A, int * LDA, double * X, int * INCX );
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
	double yoyo=dot(yo,yo);

	//Ya Variables//
	Col<double> ya(na,fill::zeros);
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

	//Allocate Space for MCMC Draws//
	Mat<double> ya_mcmc(na,niter,fill::zeros);
	Mat<double> prob_mcmc(p,niter,fill::zeros);
	Mat<uword>  gamma_mcmc(p,niter,fill::ones);
	Col<double> phi_mcmc(niter,fill::ones);
	ya_mcmc.col(0)=ya;
	phi_mcmc(0)=phi;
	gamma_mcmc.col(0)=gamma;
	prob_mcmc.col(0)=prob;

	//LAPACK Variables//
	char uplo = 'U';
	char transN = 'N';
	char transT = 'T';
	char diag = 'N';
	int info;
	int nrhs=1;
	int incx=1;
	int p_gamma;

	//Begin Gibbs Sampling Algorithm//
	for (int t = 1; t < niter; ++t)
	{
		//Form Submatrices//
		p_gamma=sum(gamma);
		if(gamma_diff)
		{
			if(p_gamma!=0){
				inc_indices=find(gamma);
				xag=xa.cols(inc_indices);
				xogyo=xoyo.elem(inc_indices);
				lamg=lam.elem(inc_indices);
				xogxog_Lamg=xoxo.submat(inc_indices,inc_indices);
				xogxog_Lamg.diag()+=lamg; 
				//Positive Definite Cholesky Factorization//
				dpotrf_(&uplo, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &info); //xerbla handles info error
				//xogxog_Lamg now stores R upper triangular where xogxog_Lamg=R'R
				Bg=xogyo; 
				//Triangular Positive Definite Solve via Cholesky
				dpotrs_(&uplo,  &p_gamma, &nrhs, &*xogxog_Lamg.begin(), &p_gamma, &*Bg.begin(),  &p_gamma, &info);

				b=0.5*(yoyo-dot(xogyo, Bg)); 
				mu=xag*Bg;
			}else{
				b=0.5*yoyo;
			}
		}

		//Draw Phi//
		phi=R::rgamma(a,(1/b)); //rgamma uses scale

		//Draw Ya//
		Z.set_size(na);
		Z.randn();
		if(p_gamma!=0){
			ya=mu+sqrt(1/phi)*Z;
			Z.set_size(p_gamma);
			Z.randn();
			//Computes R^{-1}Z where xogxog_Lamg^{-1}=R^{-1}R^{-T}
			dtrsv_(&uplo, &transN, &diag, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &*Z.begin(), &incx);
			ya+=sqrt(1/phi)*xag*Z;
			xaya=ya;
			dtrmv_(&uplo, &transT, &diag, &na, &*xa.begin(), &na, &*xaya.begin(), &incx);
		}else{
			ya=sqrt(1/phi)*Z;
			xaya=ya;
			dtrmv_(&uplo, &transT, &diag, &na, &*xa.begin(), &na, &*xaya.begin(), &incx);
		}

		//Draw Gamma//
		Bols=(1/d)%(xoyo+xaya);
		odds=priorodds%ldl%trunc_exp(0.5*phi*dli%d%d%Bols%Bols);
		prob=odds/(1+odds);
		U.randu();
		for (int i = 0; i < p; ++i)
		{
			if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN
			if(U(i)<prob(i)){
				gamma(i)=1;
			}else{
				gamma(i)=0;
			}
		}

		//Store Draws//
		gamma_mcmc.col(t)=gamma;
		prob_mcmc.col(t)=prob;
		ya_mcmc.col(t)=ya;
		phi_mcmc(t)=phi;


		if(sum(abs(gamma_mcmc.col(t)-gamma_mcmc.col(t-1)))==0)
		{
			gamma_diff=false;
		}else
		{
			gamma_diff=true;
		}
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
