#include "oda.h"

void draw_collapsed_xaya(std::vector<double> &xaya, std::vector<double> &xa, std::vector<double> &xag, std::vector<double> &mu, double phi, std::vector<double> &Z, std::vector<double> &xogxog_Lamg, int na, int p, int p_gamma){

	double sd=sqrt(1/phi);
	Z.resize(na);
	for(std::vector<double>::iterator it=Z.begin(); it!=Z.end(); ++it) *it=Rf_rnorm(0,1);
	if(p_gamma!=0){
		xaya=mu;
		daxpy_(&p, &sd, &*Z.begin(), &inc, &*xaya.begin(), &inc);
		Z.resize(p_gamma);
		for(std::vector<double>::iterator it=Z.begin(); it!=Z.end(); ++it) *it=Rf_rnorm(0,1);
		//Computes R^{-1}Z where xogxog_Lamg^{-1}=R^{-1}R^{-T}
		dtrsv_(&uplo, &transN, &unit_tri, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &*Z.begin(), &inc);
		dgemv_(&transN , &na, &p_gamma, &sd, &*xag.begin(), &na, &*Z.begin(), &inc, &inputscale1, &*xaya.begin(), &inc);
		dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
	}else{
		for(size_t i=0; i!=xaya.size(); ++i) xaya[i]=sd*Z[i];
		dtrmv_(&uplo, &transT, &unit_tri, &na, &*xa.begin(), &na, &*xaya.begin(), &inc);
	}
};
