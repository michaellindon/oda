#include "oda.h" 

void submatrices_collapsed(std::vector<double> &mu, std::vector<double> &xag, std::vector<double> &xogxog_Lamg, std::vector<double> &xogyo, std::vector<double> &lamg, std::vector<double> &Bg, const std::vector<int> &gamma, const std::vector<double> &xoyo, const std::vector<double> &lam, const std::vector<double> &xa, const std::vector<double> &xoxo,  int p_gamma, double &b, const double yoyo, int p, int na){
	xag.resize(0);
	xogxog_Lamg.resize(0);
	xogyo.resize(0);
	lamg.resize(0);
	Bg.resize(0);
	for(int i=0; i<p; ++i)
	{
		if(gamma[i]==1)
		{
			xogyo.push_back(xoyo[i]);
			Bg.push_back(xoyo[i]);
			lamg.push_back(lam[i]);
			for(int j=0; j<na; ++j) xag.push_back(xa[i*na+j]);
			for(int j=0; j<p; ++j) if(gamma[j]==1) xogxog_Lamg.push_back(xoxo[i*p+j]);
		}
	}
	for(int i=0; i<p_gamma; ++i) xogxog_Lamg[i*p_gamma+i]+=lamg[i];
	//Positive Definite Cholesky Factorization//
	dpotrf_( &uplo, &p_gamma, &*xogxog_Lamg.begin(), &p_gamma, &info); //xerbla handles info error
	//xogxog_Lamg now stores R upper triangular where xogxog_Lamg=R'R
	//Triangular Positive Definite Solve via Cholesky
	dpotrs_( &uplo, &p_gamma, &nrhs, &*xogxog_Lamg.begin(), &p_gamma, &*Bg.begin(), &p_gamma, &info);

	//b=0.5*(yoyo-std::inner_product(xogyo.begin(),xogyo.end(),Bg.begin(),0));
	b=0.5*(yoyo-ddot_(&p_gamma, &*xogyo.begin(), &inc, &*Bg.begin(), &inc));
	dgemv_( &transN , &na, &p_gamma, &unity, &*xag.begin(), &na, &*Bg.begin(), &inc, &inputscale0, &*mu.begin(), &inc);
}
