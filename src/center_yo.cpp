#include "oda.h"

double center_yo(std::vector<double> &yo)
{
	double yobar=std::accumulate(yo.begin(),yo.end(),0.0)/yo.size();
	std::transform(yo.begin(),yo.end(),yo.begin(),bind2nd(std::minus<double>(),yobar));	
	return(yobar);
}
