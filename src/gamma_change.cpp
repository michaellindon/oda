#include "oda.h"

bool gamma_change(const std::vector<int> &gamma_mcmc, int t, int p){

	for(int i=0; i<p; ++i)	if(gamma_mcmc[t*p+i]!=gamma_mcmc[(t-1)*p+i]) return(true);
	return(false);
}
