#include "oda.h" 

void submatrices_uncollapsed(bool gamma_diff, const std::vector<double> B, std::vector<double> &xog, std::vector<double> &xag,  std::vector<double> &lamg, std::vector<double> &Bg, const std::vector<int> &gamma,  const std::vector<double> &lam, const std::vector<double> &xo, const std::vector<double> &xa,  int p_gamma, int p, int no, int na){

	if(gamma_diff){
		xog.resize(0);
		xag.resize(0);
		for(int i=0; i<p; ++i)
		{
			if(gamma[i]==1)
			{
				for(int j=0; j<no; ++j) xog.push_back(xo[i*no+j]);
				for(int j=0; j<=i; ++j) xag.push_back(xa[j+(i+1)*i/2]); 
				for(int j=i+1; j!=p; ++j) xag.push_back(0);
			}
		}
	}

	lamg.resize(0);
	Bg.resize(0);
	for(int i=0; i<p; ++i)
	{
		if(gamma[i]==1)
		{
			Bg.push_back(B[i]);
			lamg.push_back(lam[i]);
		}
	}
}
