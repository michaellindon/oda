#include <stack>
#include <vector>
#include <map>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double log_posterior_density(int no,const Col<double>& lam, const Col<uword>& gamma,const Col<double>& priorodds, double a, double b,int p);

// [[Rcpp::export]]
List normal_em(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxaxa, NumericVector rd, NumericVector rlam, NumericVector rpriorprob, SEXP rselection, SEXP rmaxiter){

	std::map<std::vector<int>, bool> visitedModels;
	std::stack<std::vector<double> > stackB;
	std::stack<std::vector<int> > stackgamma;
	std::stack<double> stackphi;


	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxaxa.nrow();
	int selection=Rcpp::as<int >(rselection);
	int maxiter=Rcpp::as<int >(rmaxiter);

	//Create Data//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xaxa(rxaxa.begin(), na, p, false);
	arma::colvec d(rd.begin(),rd.size(), false);
	arma::colvec lam(rlam.begin(),rlam.size(), false);
	arma::colvec priorprob(rpriorprob.begin(),rpriorprob.size(), false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);

	//Pre-Processing//
	Col<double> xoyo=xo.t()*yo;
	Col<double> B(p,fill::zeros);
	Col<double> Bols=B;
	Col<double> mu_ya(na);
	Col<double> xaya(p);
	Mat<double> Lam=diagmat(lam);
	Col<double> prob=priorprob;
	Col<double> priorodds=priorprob/(1-priorprob);
	Col<double> odds=priorodds;
	Col<double> ldl=sqrt(lam/(d+lam));
	Col<double> dli=1/(d+lam);
	double a;
	double b;
	double deltaP;
	double deltaB;
	double phi;

	//Create Submatrices
	Col<uword> gamma(p,fill::zeros);
	Col<uword> inc_indices(p,fill::ones);
	Mat<double> xog(no,p);
	Mat<double> xaxag(na,p);
	Mat<double> Lamg(p,p);
	Col<double> Bg(p,fill::zeros);

	//Randomize Initial Gammas//
	for (int i = 0; i < p; ++i) if(R::runif(0,1)>0.5) gamma(i)=1;

	//Create Trace Matrices
	Mat<double> prob_trace(p,maxiter);
	Mat<double> B_trace(p,maxiter);
	Mat<uword>  gamma_trace(p,maxiter);
	Col<double> lpd_trace(maxiter);
	Col<double> a_trace(maxiter);
	Col<double> b_trace(maxiter);

	//Create Initial Gammas//
	//Forward Selection//
	if(selection==0) gamma.fill(0);

	//Backward Selection
	if(selection==1) gamma.fill(1);

	//Randomize Selection//
	if(selection==2){
		double init_prob=R::runif(0,1);
		for(int i=0; i<p; ++i){
			if(R::runif(0,1)<0.5){
				gamma(i)=1; //Note gamma is initialized at zero
				B(i)=R::rnorm(0,1);
			}
		}
	}

	//Run EM Algorithm//
	gamma_trace.col(0)=gamma;
	prob_trace.col(0)=prob;
	B_trace.col(0)=B;
	int t=1;
	do{
		//Form Submatrices
		inc_indices=find(gamma);
		Bg=B.elem(inc_indices);
		Lamg=Lam.submat(inc_indices,inc_indices);
		xaxag=xaxa.cols(inc_indices);
		xog=xo.cols(inc_indices);

		//E-Step Phi//
		b=(double)0.5*dot(yo-xog*Bg,yo-xog*Bg)+0.5*dot(Bg,Lamg*Bg);
		a=(double)0.5*(no+sum(gamma)-1);
		phi=a/b;
		lpd_trace(t-1)= log_posterior_density( no,lam, gamma, priorodds, a,  b, p);

		//E-Step Ya//
		xaya=xaxag*Bg;

		//M-Step (Beta,Gamma)//
		for (int i = 0; i < p; ++i)
		{
			Bols(i)=(1/d(i))*(xoyo(i)+xaya(i));
			odds(i)=priorodds(i)*ldl(i)*trunc_exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i));
			prob(i)=odds(i)/(1+odds(i));
			//if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN

			//Choose Median Probability Model//
			if(priorprob(i)*sqrt(lam(i)/(2*M_PI))*trunc_exp(0.5*Rf_digamma(a)-0.5*log(b))*trunc_exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i))>(1-priorprob(i))){
				gamma(i)=1;
				B(i)=dli(i)*d(i)*Bols(i);
			}else{
				gamma(i)=0;
				B(i)=0;
			}
		}

		//Activates when there is a transition
		if(sum(abs(gamma-gamma_trace.col(t-1)))!=0)
		{
			std::vector<int> gammastl(gamma.begin(),gamma.end());
			//Activate when model has not been visited before
			if(visitedModels.find(gammastl)==visitedModels.end()){
				Col<double> Bprev=B_trace.col(t-1);
				Col<uword> gammaprev=gamma_trace.col(t-1);
				std::vector<int> gammaprevstl(gammaprev.begin(),gammaprev.end());
				std::vector<double> Bprevstl(Bprev.begin(),Bprev.end());
				stackB.push(Bprevstl);
				stackgamma.push(gammaprevstl);
				visitedModels[gammastl]=true;

			}else{ //If gamma changes but model has already been visited, revert to old gamma
			gamma=gamma_trace.col(t-1);
			B=B%gamma;
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
	} while(deltaB>0.000000001 && t<maxiter);

	std::cout << visitedModels.size() << std::endl;

	//Resize Trace Matrices//
	gamma_trace.resize(p,t);
	prob_trace.resize(p,t);
	B_trace.resize(p,t);
	lpd_trace.resize(t-1);
	a_trace.resize(t);
	b_trace.resize(t);

	//Report Top Model//
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
			Rcpp::Named("lpd_trace") = lpd_trace
			) ;

}




double log_posterior_density(int no,const Col<double>& lam, const Col<uword>& gamma, const Col<double>& priorodds, double a, double b,int p){

	double lpd;
	Col<uword> one(p,fill::ones);

	lpd=lgamma(a)-a*log(b)+0.5*sum(gamma%log(lam))-0.5*(no+sum(gamma)-1)*log(2*M_PI)+sum(gamma%log(priorodds));
	return(lpd);

}

