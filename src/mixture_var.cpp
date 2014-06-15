#include <chrono>
//#include <google/profiler.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double lower_bound(double a,double b,const Col<double>& mu_B,const Col<double>& sigma2_B,const Col<double>& mu_ya,const Col<double>& prob,const Mat<double>& P,const Col<double>& lam, const Col<double>& yo,const Mat<double>& xc,const Mat<double>& xo,const Mat<double>& xa, const Col<double>& d, const Mat<double>& D,int no,int na, int p,const Col<double>& priorprob,const Col<double>& a_lam, const Col<double>& b_lam, double alpha);

// [[Rcpp::export]]
List mixture_var(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rlam, NumericVector rpriorprob, SEXP ralpha){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	int t=1;
	double a;
	double b;
	double phi;
	double delta;
	double alpha=Rcpp::as<double >(ralpha);
	Mat<double> xaxa(p,p);
	Mat<double> xoxo(p,p);
	Mat<double> D(p,p);
	Mat<double> Lam(p,p);
	Mat<double> xc(no+na,p);
	Mat<double> Ino=eye(no,no);
	Mat<double> Ina=eye(na,na);
	Mat<double> Inc=eye(no+na,no+na);
	Mat<double> P=eye(p,p);
	Mat<double> P1(no,no);
	Mat<double> mu_ya_trace(na,10000,fill::zeros);
	Mat<double> mu_B_trace(p,10000,fill::zeros);
	Mat<double> prob_trace(p,10000,fill::zeros);
	Col<double> b_trace(10000);
	Col<double> lb_trace(10000);
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
	Col<double> sigma2_B(p);
	Col<double> a_lam(p,fill::ones);
	Col<double> b_lam(p,fill::ones);
  Col<uword> top_model;

	//	ProfilerStart("profile.out");
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
	sigma2_B=dli;
	mu_B=solve(xoxo+Lam,xoyo); //Initialize at Ridge
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
    auto start = std::chrono::steady_clock::now();
	do{
		//Phi Maximization Step//
		b=0.5*dot(yoc-xo*P*mu_B,yoc-xo*P*mu_B)+0.5*dot(mu_B,Lam*P*mu_B)+0.5*dot(mu_B,(P-P*P)*D*mu_B)+0.5*sum(prob%sigma2_B%(lam+d));
		a=0.5*(no+sum(prob)-1);
		phi=((double)a)/b;

		//Ya Maximization Step//
		mu_ya=xa*P*mu_B;

		lb_trace[3*t-3]= lower_bound( a, b, mu_B, sigma2_B, mu_ya, prob, P,lam,yo,xc,xo,xa ,d, D, no, na, p, priorprob,a_lam,b_lam,alpha);

		//Probability Step//
		Bols=(1/d)%(xoyo+xa.t()*mu_ya);
		odds=priorodds%ldl%trunc_exp(0.5*phi*dli%d%d%Bols%Bols);
		prob=odds/(1+odds);
		P=diagmat(prob);

		mu_B=dli%d%Bols;
		sigma2_B=(dli)/phi;

		lb_trace[3*t-2]= lower_bound( a, b, mu_B, sigma2_B, mu_ya, prob, P,lam,yo,xc,xo,xa ,d, D, no, na, p, priorprob,a_lam,b_lam,alpha);

		for (int i = 0; i < p; ++i)
		{
			a_lam(i)=0.5*(alpha+prob(i));
			b_lam(i)=0.5*(alpha+phi*prob(i)*(mu_B(i)*mu_B(i)+sigma2_B(i)));
			lam(i)=(alpha+prob(i))/(alpha+phi*prob(i)*(mu_B(i)*mu_B(i)+sigma2_B(i)));
		}
		Lam=diagmat(lam);
		ldl=sqrt(lam/(d+lam));
		dli=1/(d+lam);

		lb_trace[3*t-1]= lower_bound( a, b, mu_B, sigma2_B, mu_ya, prob, P,lam,yo,xc,xo,xa ,d, D, no, na, p, priorprob,a_lam,b_lam,alpha);

		//Store Values//
		prob_trace.col(t)=prob;
		mu_ya_trace.col(t)=mu_ya;
		b_trace(t)=b;

		delta=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		t=t+1;
	}while (delta>0.001*p);

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed=end-start;

  std::cout <<  elapsed.count() << " sec - Total Runtime" << std::endl;
	std::cout <<  elapsed.count()/(t-1) << " sec - Per Iteration (avg)" << std::endl;

	prob_trace.resize(p,t);
	mu_ya_trace.resize(p,t);
	b_trace.resize(t);
	lb_trace.resize(3*(t-1));

top_model=find(prob>0.5)+1;

cout << "Top Model Predictors" << endl;
cout << top_model << endl;

	//	ProfilerStop();
	return Rcpp::List::create(
  		Rcpp::Named("top_model") = top_model,
Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("mu_ya_trace") = mu_ya_trace,
			Rcpp::Named("mu_ya") = mu_ya,
			Rcpp::Named("mu_B") = mu_B,
                	Rcpp::Named("sigma2_B") = sigma2_B,
			Rcpp::Named("b_trace") = b_trace,
			Rcpp::Named("b") = b,
			Rcpp::Named("lam") = lam,
			Rcpp::Named("lb_trace") = lb_trace,
			Rcpp::Named("a") = a
			) ;

}



double lower_bound(double a,double b,const Col<double>& mu_B,const Col<double>& sigma2_B,const Col<double>& mu_ya,const Col<double>& prob,const Mat<double>& P,const Col<double>& lam, const Col<double>& yo,const Mat<double>& xc,const Mat<double>& xo,const Mat<double>& xa, const Col<double>& d, const Mat<double>& D,int no,int na, int p,const Col<double>& priorprob,const Col<double>& a_lam, const Col<double>& b_lam, double alpha){

	double L_yc=0;
	double L_B=0;
	double L_g=0;
	double L_phi=0;
	double L_lam=0;
	double entropy_Bg=0;
	double entropy_lam=0;
	double entropy_yaphi=0;

	L_yc=0.5*(no+na-1)*(Rf_digamma(a)-log(b)-log(2*M_PI))-(0.5*na)-0.5*(a/b)*(dot(yo,yo)+dot(mu_ya,mu_ya)+dot(mu_B,P*D*mu_B)-2*dot(yo,xo*P*mu_B)-2*dot(mu_ya,xa*P*mu_B)+sum(d%prob%sigma2_B));
	for (int i = 0; i < p; i++) {
		L_B+=0.5*prob(i)*(Rf_digamma(a)-log(b)+Rf_digamma(a_lam(i))-log(b_lam(i))-log(2*M_PI))-0.5*(a/b)*prob(i)*(a_lam(i)/b_lam(i))*(mu_B(i)*mu_B(i)+sigma2_B(i));
		L_g+=prob(i)*log(priorprob(i))+(1-prob(i))*log(1-priorprob(i));
		L_lam+=0.5*alpha*log(0.5*alpha)-lgamma(0.5*alpha)+(0.5*alpha-1)*(Rf_digamma(a_lam(i))-log(b_lam(i)))-0.5*alpha*(a_lam(i)/b_lam(i));
		entropy_Bg+=(0.5*prob(i)*log(2*M_PI*sigma2_B(i))+0.5*prob(i));
		entropy_lam+=-a_lam(i)*log(b_lam(i))+lgamma(a_lam(i))-(a_lam(i)-1)*(Rf_digamma(a_lam(i))-log(b_lam(i)))+a_lam(i);
		if(prob(i)!=0 && prob(i)!=1) entropy_Bg-=(prob(i)*log(prob(i))+(1-prob(i))*log(1-prob(i)));
	}
	L_phi=-(Rf_digamma(a)-log(b));
	entropy_yaphi=-0.5*na*(Rf_digamma(a)-log(b)-log(2*M_PI))+0.5*na-a*log(b)+lgamma(a)-(a-1)*(Rf_digamma(a)-log(b))+a;
	return(L_yc+L_B+L_g+L_phi+L_lam+entropy_lam+entropy_Bg+entropy_yaphi);
}
