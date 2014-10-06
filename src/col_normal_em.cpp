#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double log_posterior_density_colcn(int no,const Mat<double>& Lamg, const Mat<double>& xoxog, const Col<uword>& gamma,const Col<double>& priorodds, double b);

// [[Rcpp::export]]
List col_normal_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP rselection){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxa.nrow();
	int selection=Rcpp::as<int >(rselection);


	//Create Data//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xa(rxa.begin(), na, p, false);
	arma::colvec d(rd.begin(),rd.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);

	//Pre-Processing//
	Col<double> xoyo=xo.t()*yo;
	Mat<double> xoxo=xo.t()*xo;
	Col<double> B(p);
	Col<double> Bols=B;
	Col<double> mu_ya(na);
	Mat<double> xat=xa.t();
	Col<double> xamu_ya(p);
	Mat<double> Lam=diagmat(lam);
	Col<double> prob=priorprob;
	Col<double> priorodds=priorprob/(1-priorprob);
	Col<double> odds=priorodds;
	Col<double> ldl=sqrt(lam/(d+lam));
	Col<double> dli=1/(d+lam);

	double a=(double)0.5*(no-1);
	double b;
	double delta;
	double yoyo=dot(yo,yo);

	//Create Submatrices
	Col<uword> gamma(p,fill::zeros);
	Col<uword> inc_indices(p);
	Mat<double> xog(no,p);
	Mat<double> xoxog(p,p);
	Mat<double> xag(na,p);
	Mat<double> Lamg(p,p);
	Mat<double> E(na,na);

	//Create Matrices//
	Mat<double> prob_trace(p,1000,fill::zeros);
	Mat<uword>  gamma_trace(p,1000,fill::ones);
	Col<double> lpd_trace(1000);

	//Create Initial Gammas//
	//Forward Selection//
	if(selection==0) gamma.fill(0);
	
	//Backward Selection
	if(selection==1) gamma.fill(1);

	//Randomize Selection//
	if(selection==2){
	double init_prob=R::runif(0,1);
	for(int i=0; i<p; ++i){
		if(R::runif(0,1)<init_prob) gamma(i)=1;
	}
	}


	//Run EM Algorithm//
	cout << "Beginning EM Algorithm" << endl;
	gamma_trace.col(0)=gamma;
	prob_trace.col(0)=prob;
	int t=1;
	do{
		//Form Submatrices
		inc_indices=find(gamma);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xag=xa.cols(inc_indices);
		xog=xo.cols(inc_indices);
		xoxog=xoxo.submat(inc_indices,inc_indices);

		B=solve(xoxog+Lamg,xog.t()*yo);

		//E-Step//
		//Phi//
		b=0.5*(yoyo - dot(B,(xoxog+Lamg)*B));
		lpd_trace(t-1)=log_posterior_density_colcn( no, Lamg, xoxog, gamma, priorodds,b);

		//Ya//
		mu_ya=xag*B;
		//	E=Ina+xag*(xoxog+Lamg).i()*xag.t(); Identity Matrix requires a lot of memory to store
		E=xag*(xoxog+Lamg).i()*xag.t();
		for(int i=0; i<p; ++i) E(i,i)+=1;

		xamu_ya=xat*mu_ya;
		//Gamma//
		for (int i = 0; i < p; ++i)
		{
			Bols(i)=(1/d(i))*(xoyo(i)+xamu_ya(i));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*(a/b)*dli(i)*d(i)*d(i)*Bols(i)*Bols(i)+0.5*dli(i)*dot(xa.col(i),E*xa.col(i)));
			prob(i)=odds(i)/(1+odds(i));
			//if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN

			//M-Step//
			//Choose Median Probability Model//
			if(0.5<prob(i)){
				gamma(i)=1;
			}else{
				gamma(i)=0;
			}
		}

		//Store Values//
		gamma_trace.col(t)=gamma;
		prob_trace.col(t)=prob;

		delta=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		t=t+1;
	} while(delta>0.0001);
	cout << "EM Algorithm Complete" << endl;

	//Report Top Model//
	Col<uword> top_model=find(gamma)+1;
	cout << "Top Model Predictors" << endl;
	cout << top_model << endl;

	//Resize Trace Matrices
	gamma_trace.resize(p,t);
	prob_trace.resize(p,t);
	lpd_trace.resize(t-1);

	return Rcpp::List::create(
			Rcpp::Named("top_model")=top_model,
			Rcpp::Named("prob") = prob,
			Rcpp::Named("prob_trace") = prob_trace,
			Rcpp::Named("gamma") = gamma,
			Rcpp::Named("gamma_trace") = gamma_trace,
			Rcpp::Named("lpd_trace") = lpd_trace
			) ;

}

double log_posterior_density_colcn(int no,const Mat<double>& Lamg, const Mat<double>& xoxog, const Col<uword>& gamma,const Col<double>& priorodds, double b){
     return(0.5*log(det(Lamg))-0.5*log(det(Lamg+xoxog))-0.5*(no-1)*log(b)+sum(gamma%log(priorodds)));
}
