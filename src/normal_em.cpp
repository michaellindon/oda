#include <chrono>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double log_posterior_density(int no,const Col<double>& lam, const Col<uword>& gamma,const Col<double>& priorodds, double a, double b,int p);

// [[Rcpp::export]]
List normal_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	int t=1;
	double a=(double)0.5*(no+na-1);
	double b=1;
	double deltaP,deltaB;
	double phi;
	Mat<double> xag;
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> D(p,p);
	Mat<double> Lam(p,p);
	Mat<double> Lamg(p,p);
	Mat<double> E(na,na);
	Mat<double> xog;
	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Mat<double> P1(no,no);
	Mat<double> prob_trace(p,10000,fill::zeros);
	Mat<double> B_trace(p,10000,fill::zeros);
	Mat<uword>  gamma_trace(p,10000,fill::ones);
	Col<double> lpd_trace(10000,fill::zeros);
	Col<double> a_trace(10000);
	Col<double> b_trace(10000);
	Col<double> one(no,fill::ones);
	Col<double> mu_ya(na,fill::zeros);
	Col<double> ya(na);
	Col<double> d(p);
	Col<double> yoc(no);
	Col<double> Bols(p);
	Col<double> B(p,fill::zeros);
	Col<double> Bg;
	Col<double> xoyo(p);
	Col<double> xamu_ya(p);
	Col<double> prob(p,fill::ones);
	Col<double> priorodds(p);
	Col<double> odds(p,fill::ones);
	Col<double> ldl(p);
	Col<double> dli(p);
	Col<uword> gamma(p,fill::zeros);
	Col<uword> inc_indices(p,fill::ones);
	Col<uword> top_model;

	//Copy RData Into Matrix Classes//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);

	//Create Matrices//
	xoyo=xo.t()*yo;
	xoxo=xo.t()*xo;
	xaxa=xa.t()*xa;
	D=xaxa+xoxo;
	d=D.diag();
	P1=one*(one.t()*one).i()*one.t();
	yoc=(Ino-P1)*yo;
	prob=priorprob;
	Lam=diagmat(lam);

	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);

	//Randomize Initial Gammas//
//	for (int i = 0; i < p; ++i) if(R::runif(0,1)>0.5) gamma(i)=1;
	B=solve(xoxo+2*no*Ina,xoyo); //Initialize at Ridge
	mu_ya=xa*B;
	a=(double)0.5*(no-1);
	b=(double)0.5*dot(yoc-xo*B,yoc-xo*B);
	Bols=(1/d)%(xoyo+xa.t()*mu_ya);
	odds=priorodds%ldl%trunc_exp(0.5*(a/b)*dli%d%d%Bols%Bols);
	prob=odds/(1+odds);
	for(int i=0; i<p; i++){
		if(R::runif(0,1)<prob(i)){
			gamma(i)=1;
		}else{
			gamma(i)=0;
		}
	}

	//Run Gibbs Sampler//
	gamma_trace.col(0)=gamma;
	prob_trace.col(0)=prob;
	auto start = std::chrono::steady_clock::now();
	do{
		//Form Submatrices
		inc_indices=find(gamma);
		Bg=B.elem(inc_indices);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);

		//E-Step//
		//Phi//
		b=(double)0.5*dot(yoc-xog*Bg,yoc-xog*Bg)+0.5*dot(Bg,Lamg*Bg);
		a=(double)0.5*(no+sum(gamma)-1);
		phi=a/b;
    lpd_trace(t-1)= log_posterior_density( no,lam, gamma, priorodds, a,  b, p);

		//Ya//
		mu_ya=xag*Bg;
		xamu_ya=xa.t()*mu_ya;

		//Gamma//
		for (int i = 0; i < p; ++i)
		{
			Bols(i)=(1/d(i))*(xoyo(i)+xamu_ya(i));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i));
			prob(i)=odds(i)/(1+odds(i));
			//if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN

			//M-Step//
			//Choose Median Probability Model//
			if(prob(i)*sqrt( ((d(i)+lam(i))*(phi))/(2*M_PI) )>(1-prob(i))){
				gamma(i)=1;
				B(i)=dli(i)*d(i)*Bols(i);
			}else{
				gamma(i)=0;
				B(i)=0;
			}
		}



		//Store Values//
		a_trace(t)=a;
		b_trace(t)=b;
		gamma_trace.col(t)=gamma;
		B_trace.col(t)=B;
		prob_trace.col(t)=prob;


		deltaP=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		deltaB=dot(B_trace.col(t)-B_trace.col(t-1),B_trace.col(t)-B_trace.col(t-1));
		t=t+1;
	} while(deltaP>0.00000001 || deltaB>0.00000000001);

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed=end-start;

	std::cout <<  elapsed.count() << " sec - Total Runtime" << std::endl;
	std::cout <<  elapsed.count()/(t-1) << " sec - Per Iteration (avg)" << std::endl;




	gamma_trace.resize(p,t);
	prob_trace.resize(p,t);
	B_trace.resize(p,t);
	lpd_trace.resize(t-1);
	a_trace.resize(t);
	b_trace.resize(t);

	top_model=find(gamma)+1;

	cout << "Top Model Predictors" << endl;
	cout << top_model << endl;

	return Rcpp::List::create(
			Rcpp::Named("top_model")=top_model,
			Rcpp::Named("a") = a,
			Rcpp::Named("a_trace") = a_trace,
			Rcpp::Named("b") = b,
			Rcpp::Named("b_trace") = b_trace,
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("B") = B,
			Rcpp::Named("B_trace") = B_trace,
			Rcpp::Named("gamma") = gamma,
			Rcpp::Named("gamma_trace") = gamma_trace,
			Rcpp::Named("log_posterior_density_trace") = lpd_trace
			) ;

}




double log_posterior_density(int no,const Col<double>& lam, const Col<uword>& gamma, const Col<double>& priorodds, double a, double b,int p){

	double lpd;
	Col<uword> one(p,fill::ones);

	lpd=lgamma(a)-a*log(b)+0.5*sum(gamma%log(lam))-0.5*(no+sum(gamma)-1)*log(2*M_PI)+sum(gamma%log(priorodds));
  return(lpd);

}
