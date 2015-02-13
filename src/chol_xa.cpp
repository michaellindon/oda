#include "oda.h"

void chol_xa(std::vector<double> &xa, std::vector<double> &xaxa, std::vector<double> &xoxo, std::vector<double> &d, int p)
{
	char eigenvec='N';
	char eigenrange='I';
	double eigenprec=0.001;
	int no_eigenval=0;
	double eigenvalues[1];
	double * vl=NULL;
	double * vu=NULL;
	double * eigenvectors=NULL;
	int iwork[5*p];
	int ifail[p];
	int lwork=8*p;
	double work[lwork];
	dsyevx_(&eigenvec,&eigenrange,&uplo,&p, &*xaxa.begin(), &p,vl,vu,&p,&p,&eigenprec,&no_eigenval,&eigenvalues[0], eigenvectors, &p, &work[0], &lwork, &iwork[0], &ifail[0], &info);
	std::copy(xoxo.begin(),xoxo.end(),xaxa.begin());
	std::transform(xaxa.begin(), xaxa.end(), xaxa.begin(), bind2nd(std::multiplies<double>(),-1.0));	
	for(int i=0; i!=p; ++i)
	{
		d[i]=eigenvalues[0]+0.01;
		xaxa[i*p+i]+=d[i];
	}
	dpotrf_( &uplo, &p, &*xaxa.begin(), &p, &info); //xerbla handles info error
	dtrttp_(&uplo, &p, &*xaxa.begin(), &p, &*xa.begin(), &info);
}


