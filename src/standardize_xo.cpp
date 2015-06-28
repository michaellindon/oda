#include "oda.h"

void standardize_xo(std::vector<double> &xo, int no, int p)
{
	for(int c=0; c!=p; ++c)
	{
		double xobar=std::accumulate(xo.begin()+c*no,xo.begin()+(c+1)*no,0.0)/no;
		std::transform(xo.begin()+c*no, xo.begin()+(c+1)*no, xo.begin()+c*no, bind2nd(std::minus<double>(),xobar));	
		double xoscale=ddot_(&no, &*(xo.begin()+c*no), &inc, &*(xo.begin()+c*no), &inc);
		xoscale=sqrt(no/xoscale);
		std::transform(xo.begin()+c*no,xo.begin()+(c+1)*no,xo.begin()+c*no,bind2nd(std::multiplies<double>(),xoscale));	
	}
	//Possibly return list of scalings//
}
