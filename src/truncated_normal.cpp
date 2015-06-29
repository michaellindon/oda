#include "oda.h"

extern "C" double truncated_normal(double mu, double mu_minus, double var){
	if(mu_minus<=mu)
	{
		double proposal;
		do{
			proposal=mu+sqrt(var)*Rf_rnorm(0,1);
		}while(proposal<mu_minus);
		return proposal;
	}else{
		double rate=0.5*(mu_minus+sqrt(mu_minus*mu_minus+4));
		double proposal,u,prob;
		do{
			proposal=mu_minus+Rf_rexp(rate);
			prob=exp(-0.5*(proposal-rate)*(proposal-rate));
			u=Rf_runif(0,1);
		}while(u>prob);
		return(mu+sqrt(var)*proposal);
	}
}

//References
//[1] Christian P. Robert. Simulation of truncated normal variables.
//    Statistics and Computing, 5(2):121–125, June 1995.
