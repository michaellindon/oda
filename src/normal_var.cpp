#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List normal_var(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	int t=1;
	double a=(no-1)/2;
	double b;
	double phi;
	double varyo;
	double delta;
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> xoxog(p,p);
	Mat<double> D(p,p);
	Mat<double> Q(p,p);
	Mat<double> Lam(p,p);
	Mat<double> xcxcLami(p,p);
	Mat<double> H(na,na,fill::eye);
	Mat<double> Ino=eye(no,no);
	Mat<double> P=eye(p,p);
	Mat<double> Pi=eye(p,p);
	Mat<double> Ina=eye(na,na);
	Mat<double> P1(no,no);
	Mat<double> Px(no,no);
	Mat<double> mu_trace(na,1000,fill::zeros);
	//Mat<double> H_trace(na*na,1000,fill::zeros);
	Mat<double> prob_trace(p,1000,fill::zeros);
	Col<double> b_trace(1000);
	Col<double> one(no,fill::ones);
	Col<double> mu(na,fill::zeros);
	Col<double> d(p);
	Col<double> Bols(p);
	Col<double> xoyo(p);
	Col<double> prob(p,fill::ones);
	Col<double> priorodds(p);
	Col<double> odds(p);
	Col<double> ldl(p);
	Col<double> dli(p);

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
	Lam=diagmat(lam);
	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);
	xcxcLami=diagmat(1/(d+lam));
	P1=one*(one.t()*one).i()*one.t();
	Px=xo*(xoxo).i()*xo.t();
	varyo=dot(yo,(Ino-P1)*yo);

	//Randomize Initial Probabilities//
	for (int i = 0; i < p; ++i) prob(i)=R::runif(0,1);
	P.diag()=prob;
	Pi.diag()=1/prob;

	//Run Variational//
	prob_trace.col(0)=prob;
	mu_trace.col(0)=mu;
	//H_trace.col(0)=vectorise(H);
	b_trace(0)=b;
	do{
    cout << t << endl;
		//Phi Step//
		Q=(D+Lam)*Pi;
		b=0.5*(varyo-dot(xoyo,solve(Q-xaxa,xoyo)));
		phi=((double)a)/b;

		//Ya Step//
		mu=solve(Ina-xa*P*xcxcLami*xa.t(),xa*P*xcxcLami*xoyo);
		H=Ina-xa*P*xcxcLami*xa.t();

		//Probability Step//
		for (int i = 0; i < p; i++)
		{
			Bols(i)=(1/d(i))*(xoyo(i)+dot(xa.col(i),mu));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*phi*dli(i)*(d(i)*d(i)*Bols(i)*Bols(i)+0.5*dli(i)*dot(xa.col(i),solve(H,xa.col(i))) ));
			prob(i)=odds(i)/(1+odds(i));
		}
		P.diag()=prob;
		Pi.diag()=1/prob;

		//Store Values//
		prob_trace.col(t)=prob;
		mu_trace.col(t)=mu;
	//	H_trace.col(t)=vectorise(H);
		b_trace(t)=b;

		delta=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
    cout << delta << endl;
		t=t+1;
	}while (delta>0.001);

	prob_trace.resize(p,t);
	mu_trace.resize(p,t);
//	H_trace.resize(na*na,t);
	b_trace.resize(t);

	return Rcpp::List::create(
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("mu_trace") = mu_trace,
			Rcpp::Named("mu") = mu,
		//	Rcpp::Named("H_trace") = H_trace,
			Rcpp::Named("H") = H,
			Rcpp::Named("b_trace") = b_trace,
			Rcpp::Named("b") = b
			) ;

}
