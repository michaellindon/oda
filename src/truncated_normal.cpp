#include "oda.h"

extern "C" void truncated_normal(double * result, double *mu, double *mu_minus, double *var){
	GetRNGstate();
	if(*mu_minus<=*mu)
	{
		double proposal;
		do{
			proposal=*mu+sqrt(*var)*Rf_rnorm(0,1);
		}while(proposal<*mu_minus);
	//	return proposal;
	*result=proposal;
	}else{
		*mu_minus=(*mu_minus-*mu)/sqrt(*var);
		double rate=0.5*(*mu_minus+sqrt(*mu_minus * *mu_minus+4));
		double proposal,u,prob;
		do{
			proposal=*mu_minus+Rf_rexp(rate);
			prob=exp(-0.5*(proposal-rate)*(proposal-rate));
			u=Rf_runif(0,1);
		}while(u>prob);
	//	return(*mu+sqrt(*var)*proposal);
		*result=*mu+sqrt(*var)*proposal;
	}
	PutRNGstate();
}

//References
//[1] Christian P. Robert. Simulation of truncated normal variables.
//    Statistics and Computing, 5(2):121–125, June 1995.
