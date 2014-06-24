//#include <google/profiler.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double lower_bound(double a,double b,const Col<double>& mu_ya,const Col<double>& prob,const Mat<double>& P,const Col<double>& lam, const Col<double>& yo,const Mat<double>& xc,const Mat<double>& xo,const Mat<double>& xa, const Col<double>& d, const Mat<double>& D,int no,int na, int p,const Col<double>& priorprob,const Mat<double>& Aaai,const Col<double>& ldl,const Col<double>& dli);

// [[Rcpp::export]]
List col_normal_var(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob){
  
	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	int t=1;
	double a=(double)0.5*(no-1);
	double b;
	double phi;
	double delta;
	double lb=0;
	Mat<double> Aooa(no,no);
	Mat<double> A(no+na,no+na);
	Mat<double> Aoo(no,no);
	Mat<double> Aaa(na,na);
	Mat<double> Aoa(no,na);
	Mat<double> Aao(na,no);
	Mat<double> Aaai(na,na);
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> xcxcLami(p,p);
	Mat<double> D(p,p);
	Mat<double> Lam(p,p);
	Mat<double> xc(no+na,p);
	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Mat<double> Inc=eye(no+na,no+na);
	Mat<double> P=eye(p,p);
	Mat<double> P1(no,no);
	Mat<double> mu_trace(na,1000,fill::zeros);
	Mat<double> prob_trace(p,1000,fill::zeros);
	Col<double> lb_trace(1000);
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
	xcxcLami=diagmat(1/(d+lam));
	Lam=diagmat(lam);
	priorodds=priorprob/(1-priorprob);
	ldl=sqrt(lam/(d+lam));
	dli=1/(d+lam);
	P1=one*(one.t()*one).i()*one.t();

	//Randomize Initial Values//
	for (int i = 0; i < p; ++i) prob(i)=R::runif(0,1);
	P=diagmat(prob);
	A=Inc-xc*P*xcxcLami*xc.t();
	Aoo=A.submat(0,0,no-1,no-1)-P1;
	Aaa=A.submat(no,no,no+na-1,no+na-1);
	Aoa=A.submat(0,no,no-1,no+na-1);
	Aao=A.submat(no,0,no+na-1,no-1);
	Aaai=inv_sympd(Aaa);
	Aooa=Aoo-Aoa*Aaai*Aao;
	b=0.5*(dot(yo,Aooa*yo));
	phi=((double)a)/b;
	mu=-Aaai*Aao*yo;

	//Run Variational//
	prob_trace.col(0)=prob;
	mu_trace.col(0)=mu;
	b_trace(0)=b;
	do{


		//Probability Step//
		for (int i = 0; i < p; i++)
		{
			Bols(i)=(1/d(i))*(xoyo(i)+dot(xa.col(i),mu));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i)+0.5*dli(i)*dot(xa.col(i),Aaai*xa.col(i)) );
			prob(i)=odds(i)/(1+odds(i));

		}
		P=diagmat(prob);

 cout <<  lower_bound( a, b, mu, prob, P,lam,yo,xc,xo,xa ,d, D, no, na, p, priorprob,Aaai,ldl,dli) << endl;

		//Phi Step//
		A=Inc-xc*P*xcxcLami*xc.t();
		Aoo=A.submat(0,0,no-1,no-1)-P1;
		Aaa=A.submat(no,no,no+na-1,no+na-1);
		Aoa=A.submat(0,no,no-1,no+na-1);
		Aao=A.submat(no,0,no+na-1,no-1);
		Aaai=inv_sympd(Aaa);
		Aooa=Aoo-Aoa*Aaai*Aao;
		b=0.5*(dot(yo,Aooa*yo));
		phi=((double)a)/b;

		//Ya Step//
		mu=-Aaai*Aao*yo;


  cout << lower_bound( a, b, mu, prob, P,lam,yo,xc,xo,xa ,d, D, no, na, p, priorprob,Aaai,ldl,dli) << endl;
		
		//Compute Lower Bound//
		/*lb=-a*log(2*M_PI)+lgamma(a)-a*log(b)-0.5*log(det(Aaa))+0.5*sum(prob%log(ldl));
		for (int i = 0; i < p; ++i)
		{
			if(prob(i)!=1 && prob(i)!=0){
				lb-=(prob(i)*log(prob(i))+(1-prob(i))*log(1-prob(i)));
			}
		}*/


		//Store Values//
		prob_trace.col(t)=prob;
		mu_trace.col(t)=mu;
		b_trace(t)=b;
		lb_trace(t)=lb;

		delta=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		t=t+1;
	}while (delta>0.001);

	prob_trace.resize(p,t);
	mu_trace.resize(p,t);
	b_trace.resize(t);
  lb_trace.resize(2*(t-1));

//ProfilerStop();
	return Rcpp::List::create(
			Rcpp::Named("lb_trace") = lb_trace,
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("mu_trace") = mu_trace,
			Rcpp::Named("mu") = mu,
			Rcpp::Named("E") = Aaai,
			Rcpp::Named("b_trace") = b_trace,
			Rcpp::Named("b") = b
			) ;

}





double lower_bound(double a,double b,const Col<double>& mu_ya,const Col<double>& prob,const Mat<double>& P,const Col<double>& lam, const Col<double>& yo,const Mat<double>& xc,const Mat<double>& xo,const Mat<double>& xa, const Col<double>& d, const Mat<double>& D,int no,int na, int p,const Col<double>& priorprob, const Mat<double>& Aaai,const Col<double>& ldl,const Col<double>& dli){

	double L_yc=0;
	double L_g=0;
	double L_phi=0;
	double entropy_g=0;
	double entropy_yaphi=0;

	L_yc=0.5*(no+na-1)*(Rf_digamma(a)-log(b)-log(2*M_PI))-0.5*(a/b)*(dot(yo,yo)+dot(mu_ya,mu_ya))-0.5*trace(Aaai);
	for (int i = 0; i < p; i++) {
		L_yc+=0.5*prob(i)*log(ldl(i))+0.5*prob(i)*dli(i)*dot(xa.col(i),Aaai*xa.col(i))+0.5*(a/b)*prob(i)*dli(i)*(dot(xo.col(i),yo)+dot(xa.col(i),mu_ya))*(dot(xo.col(i),yo)+dot(xa.col(i),mu_ya));
		L_g+=prob(i)*log(priorprob(i))+(1-prob(i))*log(1-priorprob(i));
		if(prob(i)!=0 && prob(i)!=1) entropy_g-=(prob(i)*log(prob(i))+(1-prob(i))*log(1-prob(i)));
	}
	L_phi=-(Rf_digamma(a)-log(b));
	entropy_yaphi=-(0.5*(no+na-1)-1)*(Rf_digamma(a)-log(b))+0.5*log(det(Aaai))+0.5*na*log(2*M_PI)-a*log(b)+lgamma(a)+a;
	return(L_yc+L_g+L_phi+entropy_g+entropy_yaphi);
}
