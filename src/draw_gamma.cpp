#include "oda.h"

void draw_gamma(std::vector<int> &gamma, int &p_gamma, const std::vector<double> prob){
	for(size_t i = 0; i!=prob.size() ; ++i)
	{
		if(Rf_runif(0,1)<prob[i]){
			gamma[i]=1;
		}else{
			gamma[i]=0;
		}
	}
	p_gamma=std::accumulate(gamma.begin(),gamma.end(),0);
}
