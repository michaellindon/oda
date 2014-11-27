#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double log_posterior_density_t_bb(const Col<double>& Bg, const Col<uword>& gamma, double phi,const Col<double>& yo,const Mat<double>& xog,double alpha, double beta, int no, double nu,int p);

// [[Rcpp::export]]
List t_em_bb(NumericVector ryo, NumericMatrix rxo, NumericMatrix rxaxa, NumericVector rd, SEXP ralpha, SEXP rbeta, SEXP rnu, SEXP rselection, SEXP rmaxiter){

	//Define Variables//
	int p=rxo.ncol();
	int no=rxo.nrow();
	int na=rxaxa.nrow();
	int selection=Rcpp::as<int >(rselection);
	int maxiter=Rcpp::as<int >(rmaxiter);
	double alpha=Rcpp::as<double >(ralpha);
	double beta=Rcpp::as<double >(rbeta);

	//Create Data//
	arma::mat xo(rxo.begin(), no, p, false);
	arma::mat xaxa(rxaxa.begin(), na, p, false);
	arma::colvec d(rd.begin(),rd.size(), false);
	arma::colvec yo(ryo.begin(), ryo.size(), false);
	yo-=mean(yo);

	//Pre-Processing//
	Col<double> lam(p,fill::ones);
	Col<double> xoyo=xo.t()*yo;
	Col<double> B(p,fill::randn);
	Col<double> Bols=B;
	Col<double> ya(na);
	Col<double> xaya(p);
	Col<double> prob(p,fill::zeros);;
	double priorodds=0.5;
	Col<double> odds(p);
	Col<double> ldl=sqrt(lam/(d+lam));
	Col<double> dli=1/(d+lam);
	Col<double> explam=lam;
	double a;
	double b;
	double deltaP;
	double deltaB;
	uword max_index;
	double phi=(no+p-3)/dot(yo-xo*B,yo-xo*B);
	double nu=Rcpp::as<double >(rnu);

	//Create Submatrices
	Col<uword> gamma(p,fill::zeros);
	Col<uword> inc_indices(p,fill::ones);
	Mat<double> xog=xo;
	Mat<double> xaxag=xaxa;
	Col<double> Bg=B;
	Col<double> lamg=lam;

	//Create Trace Matrices
	Mat<double> prob_trace(p,maxiter);
	Mat<double> B_trace(p,maxiter);
	Mat<uword>  gamma_trace(p,maxiter);
	Col<double> lpd_trace(maxiter);
	Col<double> phi_trace(maxiter);

	//Create Initial Gammas//
	//Forward Selection//
	if(selection==0) {
		gamma.fill(0);
		B.fill(0.0);
		Bg.fill(0.0);
	}

	//Backward Selection
	if(selection==1) gamma.fill(1);

	//Randomize Selection//
	if(selection==2){
		double init_prob=R::runif(0,1);
		for(int i=0; i<p; ++i){
			if(R::runif(0,1)<1/p) gamma(i)=1; //Note gamma initialized at zero
			//if(R::runif(0,1)<(double)(4.0/(double)p)) gamma(i)=1; //Note gamma initialized at zero
		}
	}

	//Run EM Algorithm//
	gamma_trace.col(0)=gamma;
	prob_trace.col(0)=prob;
	B_trace.col(0)=B;
	phi_trace(0)=phi;
	int t=1;
	do{

		lpd_trace(t-1)=log_posterior_density_t_bb(Bg, gamma, phi,yo,xog, alpha, beta, no, nu,p);

		//Expectation//
		lam=(nu+gamma)/(nu+phi*gamma%B%B);
		for(int i=0; i<p; ++i){
			explam(i)=trunc_exp(Rf_digamma(nu+gamma(i)))/(nu+gamma(i)*phi*B(i)*B(i));
			dli(i)=1/(d(i)+lam(i));
			ldl(i)=sqrt(explam(i)/(d(i)+explam(i)));
		}

		//CM-Step (Beta,Gamma)|phi//
		xaya=xaxag*Bg;
		priorodds=trunc_exp(Rf_digamma(alpha+sum(gamma))-Rf_digamma(p-sum(gamma)+beta));
		for (int i = 0; i < p; ++i)
		{

			Bols(i)=(1/d(i))*(xoyo(i)+xaya(i));
			odds(i)=priorodds*ldl(i)*trunc_exp(0.5*phi*dli(i)*d(i)*d(i)*Bols(i)*Bols(i));
			prob(i)=odds(i)/(1+odds(i));
			if(prob(i)!=prob(i)) prob(i)=1;	 //Catch NaN

			//Choose Median Probability Model//
			if(sqrt(phi*(d(i)+explam(i))/(2*M_PI))*prob(i)>(1-prob(i))){
				gamma(i)=1;
				B(i)=dli(i)*d(i)*Bols(i);
			}else{
				gamma(i)=0;
				B(i)=0;
			}
		}

		//Avoid Null Model Entrapment
	/*	if(t<=3 && sum(gamma)==0){
			prob.max(max_index);
			gamma(max_index)=1;
			B(max_index)=dli(max_index)*d(max_index)*Bols(max_index);
		}*/

		//Form Submatrices
		inc_indices=find(gamma);
		Bg=B.elem(inc_indices);
		lamg=lam.elem(inc_indices);
		xaxag=xaxa.cols(inc_indices);
		xog=xo.cols(inc_indices);

		//CM-Step Phi|(Beta,Gamma)//
		b=(double)0.5*dot(yo-xog*Bg,yo-xog*Bg)+0.5*dot(Bg,lamg%Bg);
		a=(double)0.5*(no+sum(gamma)-3);
		phi=a/b;

		//Store Values//
		phi_trace(t)=phi;
		gamma_trace.col(t)=gamma;
		B_trace.col(t)=B;
		prob_trace.col(t)=prob;

		//deltaP=dot(prob_trace.col(t)-prob_trace.col(t-1),prob_trace.col(t)-prob_trace.col(t-1));
		deltaB=dot(B_trace.col(t)-B_trace.col(t-1),B_trace.col(t)-B_trace.col(t-1));
		t=t+1;
	} while(deltaB>0.00000000000000000001 && t<maxiter);


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




double log_posterior_density_t_bb(const Col<double>& Bg, const Col<uword>& gamma, double phi,const Col<double>& yo,const Mat<double>& xog,double alpha, double beta, int no, double nu, int p){

	double lpd;

	lpd=0.5*(no+sum(gamma)-3)*log(phi)-0.5*phi*dot(yo-xog*Bg,yo-xog*Bg)-0.5*(nu+1)*sum(log(1+phi*Bg%Bg/nu))+Rf_lgammafn(alpha+beta)-Rf_lgammafn(alpha)-Rf_lgammafn(beta)+Rf_lgammafn(sum(gamma)+alpha)+Rf_lgammafn(p-sum(gamma)+beta)-Rf_lgammafn(p+alpha+beta);
	return(lpd);

}
