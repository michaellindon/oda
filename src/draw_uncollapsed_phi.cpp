#include "oda.h"

//Draw Gamma for Uncollapsed Gibbs. Note Rf_rgamma uses scale //

double draw_uncollapsed_phi(int p_gamma, int no, const std::vector<double> &yo, std::vector<double> &xog, std::vector<double> &Bg, const std::vector<double> &lamg, const double yoyo)
{
	if(p_gamma){
		double a=0.5*(no+p_gamma-1);
		std::vector<double> residual=yo;
		dgemv_(&transN , &no, &p_gamma, &nunity, &*xog.begin(), &no, &*Bg.begin(), &inc, &inputscale1, &*residual.begin(), &inc);
		double b=0.5*ddot_(&no, &*residual.begin(), &inc, &*residual.begin(), &inc);
		for(size_t i=0; i!=Bg.size(); ++i) b+=0.5*Bg[i]*Bg[i]*lamg[i]; 
		return(Rf_rgamma(a,(1/b))); 
	}else{
		double a=0.5*(no-1);
		double b=0.5*yoyo;
		return(Rf_rgamma(a,(1/b)));
	}


}
