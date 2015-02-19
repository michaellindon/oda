#include "oda.h"

void draw_xoyo(const std::vector<double> &Z, std::vector<double> &xoyo, std::vector<double> &xo, std::vector<double> &xog, std::vector<double> Bg,  double phi, int no, int p, int p_gamma)
{
	std::vector<double> U(no);
	for(std::vector<double>::iterator it=U.begin(); it!=U.end(); ++it) *it=Rf_runif(0,1);
	std::vector<double> yo(no,0.0);

	//Set Mean Values//
	if(p_gamma!=0)	dgemv_(&transN , &no, &p_gamma, &unity, &*xog.begin(), &no, &*Bg.begin(), &inc, &inputscale0, &*yo.begin(), &inc);

	for (size_t i = 0; i < Z.size(); ++i)
	{
		if(Z[i]==1){
			double Fa=pnorm(-yo[i],0,1,1,0);
			yo[i]+=qnorm(Fa+U[i]*(1-Fa),0,1,1,0);
		}else{
			yo[i]+=qnorm(U[i]*pnorm(-yo[i],0,1,1,0),0,1,1,0);
		}

	}

	//Multiply by Xo'//
	dgemv_( &transT, &no, &p, &unity, &*xo.begin(), &no, &*yo.begin(), &inc, &inputscale0, &*xoyo.begin(), &inc);
}
