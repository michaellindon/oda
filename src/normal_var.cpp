#include <google/profiler.h>
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
	double a;
	double b;
	double phi;
	double delta;
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> D(p,p);
	Mat<double> Lam(p,p);
	Mat<double> xc(no+na,p);
	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Mat<double> Inc=eye(no+na,no+na);
	Mat<double> P=eye(p,p);
	Mat<double> E_B=eye(p,p);
	Mat<double> P1(no,no);
	Mat<double> mu_ya_trace(na,10000,fill::zeros);
	Mat<double> mu_B_trace(p,10000,fill::zeros);
	Mat<double> prob_trace(p,100000,fill::zeros);
	Col<double> b_trace(10000);
	Col<double> one(no,fill::ones);
	Col<double> mu_ya(na,fill::zeros);
	Col<double> mu_B(p,fill::zeros);
	Col<double> d(p);
	Col<double> Bols(p);
	Col<double> yoc(no);
	Col<double> xoyo(p);
	Col<double> prob(p,fill::ones);
	Col<double> priorodds(p);
	Col<double> odds(p);
	Col<double> ldl(p);
	Col<double> dli(p);

//  ProfilerStart("profile.out");
	//Copy RData Into Matrix Classes//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);

	//Create Matrices//
	xc=join_cols(xo,xa);
	xoyo=xo.t()*yo;
	xoxo=xo.t()*xo;
	xaxa=xa.t()*xa;
	D=xaxa+xoxo;
	d=D.diag();
	D=diagmat(d);
	Lam=diagmat(lam);
	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);
	E_B=diagmat(dli);
	mu_B=solve(xo.t()*xo,xo.t()*yo);
	mu_ya=xa*mu_B;
	P1=one*(one.t()*one).i()*one.t();
	yoc=(Ino-P1)*yo;
  phi=1;

	//Randomize Initial Probabilities//
	//for (int i = 0; i < p; ++i) prob(i)=R::rbeta(1,10); // prob(i)=R::runif(0,1);
	prob=rbeta(p,0.5,0.5);
  P=diagmat(prob);

	//Run Variational//
	prob_trace.col(0)=prob;
	mu_ya_trace.col(0)=mu_ya;
	b_trace(0)=b;
	do{
		//Phi Maximization Step//
		//b=0.5*dot(yoc-xo*P*mu_B,yoc-xo*P*mu_B)+0.5*dot(mu_B,Lam*P*mu_B)+0.5*trace(D*P*E_B);
    b=0.5*dot(yoc-xo*P*mu_B,yoc-xo*P*mu_B)+0.5*dot(mu_B,Lam*P*mu_B)+0.5*dot(mu_B,(P-P*P)*D*mu_B)+0.5*sum(prob)/phi;
		a=0.5*(no+sum(prob)-1);
		phi=((double)a)/b;

		//Ya Maximization Step//
		mu_ya=xa*P*mu_B;

		//Probability Step//
		Bols=(1/d)%(xoyo+xa.t()*mu_ya);
		odds=priorodds%ldl%trunc_exp(0.5*phi*dli%d%d%Bols%Bols);
		prob=odds/(1+odds);
		P=diagmat(prob);

		mu_B=dli%d%Bols;
		E_B=diagmat(dli)/phi;

		//Store Values//
		prob_trace.col(t)=prob;
		mu_ya_trace.col(t)=mu_ya;
		b_trace(t)=b;
		delta=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		t=t+1;
	}while (delta>0.001);

	prob_trace.resize(p,t);
	mu_ya_trace.resize(p,t);
	b_trace.resize(t);

//ProfilerStop();
	return Rcpp::List::create(
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("mu_ya_trace") = mu_ya_trace,
			Rcpp::Named("mu_ya") = mu_ya,
			Rcpp::Named("mu_B") = mu_B,
			Rcpp::Named("b_trace") = b_trace,
			Rcpp::Named("b") = b,
			Rcpp::Named("a") = a
			) ;

}
