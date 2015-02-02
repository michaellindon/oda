#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//double log_posterior_density(int no,const Col<double>& lam, const Col<uword>& gamma,const Col<double>& priorodds, double a, double b,int p);

// [[Rcpp::export]]
List oda_gdp_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rpriorprob, SEXP ralpha, SEXP reta){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();

	//Create Data//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec d(rd.begin(),rd.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);

	//Pre-Processing//
	Col<double> lam(p,fill::ones);
	Col<double> xoyo=xo.t()*yo;
	Col<double> B(p,fill::ones);
	Col<double> Babs=B;
	Col<double> Bols=B;
	Col<double> ya(na);
	Mat<double> xat=xa.t();
	Col<double> xaya(p);
	Col<double> prob=priorprob;
	Col<double> priorodds=priorprob/(1-priorprob);
	Col<double> odds=priorodds;
	Col<double> ldl=sqrt(lam/(d+lam));
	Col<double> dli=1/(d+lam);
	Col<double> explam=lam;
	double a;
	double b;
	double deltaP;
	double deltaB;
	double phi=1;
	double eta=Rcpp::as<double >(reta);
	double alpha=Rcpp::as<double >(ralpha);

	//Create Submatrices
	Col<uword> gamma(p,fill::ones);
	Col<uword> inc_indices(p,fill::ones);
	Mat<double> xog=xo;
	Mat<double> xag=xa;
	Col<double> Bg=B;
	Col<double> lamg=lam;

	//Create Trace Matrices
	Mat<double> prob_trace(p,10000);
	Mat<double> B_trace(p,10000);
	Mat<uword>  gamma_trace(p,10000);
	Col<double> lpd_trace(10000);
	Col<double> phi_trace(10000);

	//Run EM Algorithm//
	gamma_trace.col(0)=gamma;
	prob_trace.col(0)=prob;
	int t=1;
	do{



		//	lpd_trace(t-1)= log_posterior_density( no,lam, gamma, priorodds, a,  b, p);

		//Expectation//
		ya=xag*Bg;
		Babs=abs(B);
		lam=(alpha+1)/(eta*sqrt(phi)*Babs + phi*Babs%Babs);
		for(int i=0; i<p; ++i){
			explam(i)=lam(i);
			dli(i)=1/(d(i)+lam(i));
			ldl(i)=sqrt(1/(d(i)/lam(i)+1));
		}


		//M-Step (Beta,Gamma)//
		xaya=xat*ya;
		for (int i = 0; i < p; ++i)
		{

			Bols(i)=(1/d(i))*(xoyo(i)+xaya(i));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i));
			prob(i)=odds(i)/(1+odds(i));
			//if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN

			//Choose Median Probability Model//
			if(sqrt(phi*(d(i)+1)/(2*M_PI))*prob(i)>(1-prob(i))){
				gamma(i)=1;
				B(i)=dli(i)*d(i)*Bols(i);
			}else{
				gamma(i)=0;
				B(i)=0;
			}
		}

		//Form Submatrices
		inc_indices=find(gamma);
		Bg=B.elem(inc_indices);
		lamg=lam.elem(inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);


		b=(double)0.5*dot(yo-xog*Bg,yo-xog*Bg)+0.5*dot(Bg,lamg%Bg);
		a=(double)0.5*(no+sum(gamma)-3);
		phi=a/b;

		//Store Values//
		phi_trace(t)=phi;
		gamma_trace.col(t)=gamma;
		B_trace.col(t)=B;
		prob_trace.col(t)=prob;

		deltaP=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		deltaB=dot(B_trace.col(t)-B_trace.col(t-1),B_trace.col(t)-B_trace.col(t-1));
		t=t+1;
	} while(deltaP>0.0000001 || deltaB>0.00000000000001);


	//Resize Trace Matrices//
	gamma_trace.resize(p,t);
	prob_trace.resize(p,t);
	B_trace.resize(p,t);
	lpd_trace.resize(t-1);
	phi_trace.resize(t);

	//Report Top Model//
	Col<uword> top_model=find(gamma)+1;
	cout << "Top Model Predictors" << endl;
	cout << top_model << endl;

	return Rcpp::List::create(
			Rcpp::Named("top_model")=top_model,
			Rcpp::Named("phi") = phi,
			Rcpp::Named("phi_trace") = phi_trace,
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("B") = B,
			Rcpp::Named("B_trace") = B_trace,
			Rcpp::Named("gamma") = gamma,
			Rcpp::Named("gamma_trace") = gamma_trace,
			Rcpp::Named("lpd_trace") = lpd_trace
			) ;

}



/*
   double log_posterior_density(int no,const Col<double>& lam, const Col<uword>& gamma, const Col<double>& priorodds, double a, double b,int p){

   double lpd;
   Col<uword> one(p,fill::ones);

   lpd=lgamma(a)-a*log(b)+0.5*sum(gamma%log(lam))-0.5*(no+sum(gamma)-1)*log(2*M_PI)+sum(gamma%log(priorodds));
   return(lpd);

   }*/
