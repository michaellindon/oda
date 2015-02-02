#include "oda.h"

void fixed_probabilities(std::vector<double> &prob, std::vector<double> &odds, std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &xoyo, const std::vector<double> &xaya, const std::vector<double> &priorodds, const std::vector<double> &ldl, const std::vector<double> &dli, double phi){
	for(size_t i=0; i!=prob.size(); ++i)
	{
		Bols[i]=(1/d[i])*(xoyo[i]+xaya[i]);
		odds[i]=priorodds[i]*ldl[i]*exp(0.5*phi*dli[i]*d[i]*d[i]*Bols[i]*Bols[i]);
		prob[i]=odds[i]/(1+odds[i]);
		if(prob[i]!=prob[i]) prob[i]=1;	 //Catch NaN exponential
	}
};
