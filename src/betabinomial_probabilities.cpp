#include "oda.h"


void betabinomial_probabilities(std::vector<double> &prob, std::vector<double> &odds, std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &xoyo, const std::vector<double> &xaya, const double theta, const std::vector<double> &lam, const double phi){
	for(size_t i=0; i!=prob.size(); ++i)
	{
		Bols[i]=(1/d[i])*(xoyo[i]+xaya[i]);
		odds[i]=(theta/(1-theta))*sqrt(lam[i]/(lam[i]+d[i]))*exp(0.5*phi*((d[i]*d[i])/(d[i]+lam[i]))*Bols[i]*Bols[i]);
		prob[i]=odds[i]/(1+odds[i]);
		if(prob[i]!=prob[i]) prob[i]=1;	 //Catch NaN exponential
	}
};
