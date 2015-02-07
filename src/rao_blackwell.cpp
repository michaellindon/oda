#include "oda.h"

void rao_blackwell(double * B_rb, double * prob_rb, const std::vector<double> B, const std::vector<double> prob, int burnin, int niter)
{
	for(size_t i=0; i!=B.size(); ++i)
	{
		prob_rb[i]+=prob[i]/(niter-burnin);
		B_rb[i]+=B[i]*prob[i]/(niter-burnin);
	}
}
