#include "oda.h"

void draw_lambda_t( std::vector<double> &lam, const std::vector<int> &gamma, const double alpha, const std::vector<double> &B, const double phi)
{
	for(size_t i=0; i!=lam.size(); ++i)
	{
		if(gamma[i]==0)
		{
			lam[i]=Rf_rgamma(0.5*alpha,2/alpha); //rgamma uses scale
		}else
		{
			lam[i]=Rf_rgamma(0.5*(alpha+1),2/(alpha+phi*B[i]*B[i])); //rgamma uses scale
		}
	}
}
