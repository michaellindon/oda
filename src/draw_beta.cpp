#include "oda.h"

void draw_beta(const std::vector<int> &gamma, std::vector<double> &B, const std::vector<double> &Bols, const std::vector<double> &d, const std::vector<double> &lam, double phi)
{
	for(size_t i=0; i!=B.size(); ++i)
	{
		if(gamma[i]==0)
		{
			B[i]=0;
		}else
		{
			B[i]=(d[i]/(d[i]+lam[i]))*Bols[i]+sqrt(1/(phi*(d[i]+lam[i])))*Rf_rnorm(0,1);
		}
	}
}
