#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double log_posterior_density(int no, int p, double lasso, const Col<double>& yo, const Mat<double>& xo, const Col<double>& B,double phi);

// [[Rcpp::export]]
List lasso_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, SEXP rlasso){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	double lasso=Rcpp::as<double >(rlasso);

	//Create Data//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec d(rd.begin(),rd.size(), false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);

	//Pre-Processing//
	Col<double> xoyo=xo.t()*yo;
	Col<double> B(p,fill::ones);
	Col<double> Babs=B;
	Col<double> ya=xa*B;
	Mat<double> xat=xa.t();
	Col<double> xaya(p);
	Col<double> xcyc(p);
	double yoyo=dot(yo,yo);
	double deltaB;
	double deltaphi;
	double phi=1;
	double lp;

	//Create Submatrices
	Col<uword> inc_indices(p,fill::ones);
	Mat<double> xag(na,p);
	Col<double> Bg(p,fill::zeros);


	//Create Trace Matrices
	Mat<double> B_trace(p,10000);
	Col<double> phi_trace(10000);
	Col<double> lpd_trace(10000);

	//Run EM Algorithm//
	int t=1;
	B_trace.col(0)=B;
	phi_trace(0)=phi;
	do{
		//Form Submatrices
		inc_indices=find(B);
		Bg=B.elem(inc_indices);
		xag=xa.cols(inc_indices);

		ya=xag*Bg;

		xcyc=xoyo+xat*ya;
		lp=sqrt(lasso/phi);
		Babs=abs(B);
		for(int i=0; i<p; ++i){
			B(i)=Babs(i)/(d(i)*Babs(i)+lp)*xcyc(i);
		}
		phi=(no+na+p-3)/(yoyo+dot(ya,ya)-dot(xcyc,B));


		//Store Values//
		B_trace.col(t)=B;
		phi_trace(t)=phi;
		lpd_trace(t-1)=log_posterior_density(no,p,lasso,yo,xo,B,phi);

		deltaB=dot(B_trace.col(t)-B_trace.col(t-1),B_trace.col(t)-B_trace.col(t-1));
		deltaphi=phi_trace(t)-phi_trace(t-1);
		t=t+1;
	} while((deltaB>0.001 || deltaphi>0.001) && t<10000);


	//Resize Trace Matrices//
	B_trace.resize(p,t);
	phi_trace.resize(t);
	lpd_trace.resize(t-1);

	//Report Top Model//
	Col<uword> top_model=find(B)+1;
	cout << "Top Model Predictors" << endl;
	cout << top_model << endl;

	return Rcpp::List::create(
			Rcpp::Named("top_model")=top_model,
			Rcpp::Named("B") = B,
			Rcpp::Named("B_trace") = B_trace,
			Rcpp::Named("phi") = phi,
			Rcpp::Named("phi_trace") = phi_trace,
			Rcpp::Named("lpd_trace") = lpd_trace
			) ;

}




double log_posterior_density(int no, int p, double lasso, const Col<double>& yo, const Mat<double>& xo, const Col<double>& B,double phi){

	double lpd;
	lpd=0.5*(no-1)*log(phi/(2*M_PI))-0.5*phi*dot(yo-xo*B,yo-xo*B)+0.5*p*log(phi*lasso)-sqrt(phi*lasso)*sum(abs(B))-log(phi);
	return(lpd);

}
