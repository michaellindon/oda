#include "oda.h" 

void submatrices_collapsed(const bool gamma_diff, std::vector<double> &mu, std::vector<double> &xag, std::vector<double> &xogxog_Lamg, std::vector<double> &xogyo, std::vector<double> &lamg, std::vector<double> &Bg, const std::vector<int> &gamma, const std::vector<double> &xoyo, const std::vector<double> &lam, const std::vector<double> &xa, const std::vector<double> &xoxo,  int p_gamma, double &b, const double yoyo, int p, int na){

	if(gamma_diff)
	{
		xag.resize(0);
		xogyo.resize(0);
		for(int i=0; i<p; ++i)
		{
			if(gamma[i]==1)
			{
				xogyo.push_back(xoyo[i]);
				for(int j=0; j<=i; ++j) xag.push_back(xa[j+(i+1)*i/2]); 
				for(int j=i+1; j!=p; ++j) xag.push_back(0);
			}
		}
	}


	xogxog_Lamg.resize(0);
	lamg.resize(0);
	Bg.resize(0);
	for(int i=0; i<p; ++i)
	{
		if(gamma[i]==1)
		{
			Bg.push_back(xoyo[i]);
			lamg.push_back(lam[i]);
			for(int j=0; j<p; ++j) if(gamma[j]==1) xogxog_Lamg.push_back(xoxo[i*p+j]);
		}
	}

	for(int i=0; i<p_gamma; ++i) xogxog_Lamg[i*p_gamma+i]+=lamg[i];
	dpotrf_( &uplo, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &info); //xerbla handles info error
	dpotrs_( &uplo, &p_gamma, &nrhs, &*xogxog_Lamg.begin(), &p_gamma, &*Bg.begin(), &p_gamma, &info);
	b=0.5*(yoyo-ddot_(&p_gamma, &*xogyo.begin(), &inc, &*Bg.begin(), &inc));
	dgemv_( &transN , &na, &p_gamma, &unity, &*xag.begin(), &na, &*Bg.begin(), &inc, &inputscale0, &*mu.begin(), &inc);
}
