#include "oda.h"

double draw_collapsed_phi(double &b, int p_gamma, const int no, const double yoyo, std::vector<double> &xogyo, std::vector<double> &Bg)
{
double	a=0.5*((double)no-1.0);
if(p_gamma){
	b=0.5*(yoyo-ddot_(&p_gamma, &*xogyo.begin(), &inc, &*Bg.begin(), &inc));
}else{
	b=0.5*yoyo;
}
return(Rf_rgamma(a,(1.0/b))); //rgamma uses scale
}
