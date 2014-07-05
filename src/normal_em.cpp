#include <chrono>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double log_posterior_density(int no,const Col<double>& lam, const Col<uword>& gamma,const Col<double>& priorodds, double a, double b,int p);

// [[Rcpp::export]]
List normal_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob){


	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	arma::mat xo(rxo.begin(), no, p, false);
	//arma::mat xoxo(rxoxo.begin(), p, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec d(rd.begin(),rd.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);
  

	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Col<double> xoyo=xo.t()*yo;
//  Col<double> B=solve(xoxo+no*Ina,xoyo); //Initialize at Ridge
// Col<double> B=(1/no)*xoyo-(1/(no*no))*xo.t()*solve(Ino+(1/no)*xo*xo.t(),xo*xoyo);
  Col<double> B(p,fill::zeros);
  Col<double> Bols=B;
	Col<double> mu_ya=xa*B;
  Mat<double> xat=xa.t();
	Col<double> xamu_ya=xat*mu_ya;
  
	Mat<double> Lam=diagmat(lam);
	Col<double> prob=priorprob;
	Col<double> priorodds=priorprob/(1-priorprob);
	Col<double> odds=priorodds;
	Col<double> ldl=sqrt(lam/(d+lam));
	Col<double> dli=1/(d+lam);

	Col<uword> gamma(p,fill::zeros);
	Col<uword> inc_indices(p,fill::ones);
	//Randomize Initial Gammas//
	for (int i = 0; i < p; ++i) if(R::runif(0,1)>0.5) gamma(i)=1;

	//Create Submatrices
	Mat<double> xog;
	Mat<double> xag;
	Mat<double> Lamg(p,p);
	Col<double> Bg;

	//Create Trace Matrices
	Mat<double> prob_trace(p,10000);
	Mat<double> B_trace(p,10000);
	Mat<uword>  gamma_trace(p,10000);
	Col<double> lpd_trace(10000);
	Col<double> a_trace(10000);
	Col<double> b_trace(10000);

	double a;
	double b;
	double deltaP;
	double deltaB;
	double phi;

	//Run Gibbs Sampler//
	gamma_trace.col(0)=gamma;
	prob_trace.col(0)=prob;
	int t=1;
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
		b=(double)0.5*dot(yo-xog*Bg,yo-xog*Bg)+0.5*dot(Bg,Lamg*Bg);
		a=(double)0.5*(no+sum(gamma)-1);
		phi=a/b;
		lpd_trace(t-1)= log_posterior_density( no,lam, gamma, priorodds, a,  b, p);

		//Ya//
		mu_ya=xag*Bg;
		xamu_ya=xat*mu_ya;

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
	} while(deltaP>0.0001 || deltaB>0.0001);

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

	Col<uword> top_model=find(gamma)+1;

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
