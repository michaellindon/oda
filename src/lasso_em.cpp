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
	//	Col<double> B(p,fill::randn);
	Col<double> B=xoyo/no;
	Col<double> Babs=abs(B);
	Col<double> ya(na);
	Mat<double> xat=xa.t();
	Col<double> xaya(p);
	Col<double> xcyc(p);
	double yoyo=dot(yo,yo);
	double deltaB;
	double deltaphi;
	double phi=no/dot(yo-xo*B,yo-xo*B);
	//double phi=no/yoyo;
	double lp;

	//Create Submatrices
	Col<uword> inc_indices(p,fill::ones);
	Mat<double> xag(na,p);
	Col<double> Bg(p,fill::zeros);

	//Create Trace Matrices
	Mat<double> B_trace(p,20000);
	Col<double> phi_trace(20000);
	Col<double> lpd_trace(20000);

	//Run EM Algorithm//
	cout << "Beginning EM Algorithm" << endl;
	int t=0;
	B_trace.col(t)=B;
	phi_trace(t)=phi;
	lpd_trace(t)=log_posterior_density(no,p,lasso,yo,xo,B,phi);
	do{
		t=t+1;

		/*//Form Submatrices
		  inc_indices=find(B);
		  Bg=B.elem(inc_indices);
		  xag=xa.cols(inc_indices);*/

		//E Step
		ya=xa*B;
		xcyc=xoyo+xat*ya;
		lp=sqrt(lasso/phi);

		Babs=abs(B);
		for(int i=0; i<p; ++i){
			B(i)=(Babs(i)/(Babs(i)*d(i)+lp))*xcyc(i);
		}



		phi=(no+na+p-3)/((na/phi)+dot(yo,yo)+dot(ya,ya)-dot(xcyc,B));

		//Store Values//
		B_trace.col(t)=B;
		phi_trace(t)=phi;
		lpd_trace(t)=log_posterior_density(no,p,lasso,yo,xo,B,phi);

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




double log_posterior_density(int no, int p, double lasso, const Col<double>& yo, const Mat<double>& xo, const Col<double>& B,double phi){

	double lpd;
	lpd=(double)0.5*((double)no-1)*log(phi/(2*M_PI))-0.5*phi*dot(yo-xo*B,yo-xo*B)+0.5*(double)p*log(phi*lasso)-sqrt(phi*lasso)*sum(abs(B))-log(phi);
	return(lpd);

}
