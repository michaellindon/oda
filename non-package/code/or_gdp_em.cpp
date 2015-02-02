#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double or_gdp_log_posterior_density(int no, int p, double alpha, double eta, const Col<double>& yo, const Mat<double>& xo, const Col<double>& B,double phi);

// [[Rcpp::export]]
List or_gdp_em(NumericVector ryo, NumericMatrix rxo, SEXP ralpha, SEXP reta){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	double eta=Rcpp::as<double >(reta);
	double alpha=Rcpp::as<double >(ralpha);

	//Create Data//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);

	//Pre-Processing//
	Col<double> xoyo=xo.t()*yo;
	Col<double> B=xoyo/no;
	Col<double> Babs=abs(B);
	Mat<double> xoxo=xo.t()*xo;
	Mat<double> D=eye(p,p);
	Mat<double> Ip=eye(p,p);
	double yoyo=dot(yo,yo);
	double deltaB;
	double deltaphi;
	double phi=no/dot(yo-xo*B,yo-xo*B);
	double lp;

	//Create Trace Matrices
	Mat<double> B_trace(p,20000);
	Col<double> phi_trace(20000);
	Col<double> lpd_trace(20000);

	//Run EM Algorithm//
	cout << "Beginning EM Algorithm" << endl;
	int t=0;
	B_trace.col(t)=B;
	phi_trace(t)=phi;
	lpd_trace(t)=or_gdp_log_posterior_density(no,p,alpha,eta,yo,xo,B,phi);
	do{
		t=t+1;


		Babs=abs(B);
		D=diagmat(sqrt(((eta+sqrt(phi)*Babs)%(sqrt(phi)*Babs))/(alpha+1)));
		B=D*solve(D*xoxo*D+Ip,D*xoyo);

		phi=(no+p-3)/(yoyo-dot(xoyo,B));

		//Store Values//
		B_trace.col(t)=B;
		phi_trace(t)=phi;
		lpd_trace(t)=or_gdp_log_posterior_density(no,p,alpha,eta,yo,xo,B,phi);

		deltaB=dot(B_trace.col(t)-B_trace.col(t-1),B_trace.col(t)-B_trace.col(t-1));
		deltaphi=phi_trace(t)-phi_trace(t-1);
	} while((deltaB>0.00001 || deltaphi>0.00001) && t<19999);
	cout << "EM Algorithm Converged in " << t << " Iterations" << endl;

	//Resize Trace Matrices//
	B_trace.resize(p,t);
	phi_trace.resize(t);
	lpd_trace.resize(t);

	return Rcpp::List::create(
			Rcpp::Named("B") = B,
			Rcpp::Named("B_trace") = B_trace,
			Rcpp::Named("phi") = phi,
			Rcpp::Named("phi_trace") = phi_trace,
			Rcpp::Named("lpd_trace") = lpd_trace
			) ;

}



double or_gdp_log_posterior_density(int no, int p, double alpha, double eta, const Col<double>& yo, const Mat<double>& xo, const Col<double>& B,double phi){

	double lpd;
	double xi=eta/(sqrt(phi)*alpha);
	lpd=(double)0.5*((double)no-1)*log(phi/(2*M_PI))-p*log(2*xi)-(alpha+1)*sum(log(1+abs(B)/(alpha*xi)))-0.5*phi*dot(yo-xo*B,yo-xo*B)-log(phi);
	return(lpd);

}
