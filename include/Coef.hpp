#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "Functor.hpp"
#include "Jump.h" 
// only used to store some constants and functions
// "Jump.h" isn't a proper name

typedef vector< vector< double > > stdMat;

// Originally struct intercoefficient InterCoef_cpp()
// this part is just copy and paste since there's no elegant way to do this...
// test this function in `test_coef.cpp`
stdVec interCoef(Functor& f, double x, double interval) {
	stdVec ret(4);
	
	ret[0] = f(x) + ((11*f(x) - 18*f(x+interval) + 9*f(x+2*interval) - 
		2*f(x+3*interval))*x)/(6.*interval) + ((f(x) - (5*f(x+interval))/2. + 2*f(x+2*interval) - 
		f(x+3*interval)/2.)*pow(x,2))/pow(interval,2) + ((f(x) - 3*f(x+interval) + 
		3*f(x+2*interval) - f(x+3*interval))*pow(x,3))/(6.*pow(interval,3));

	ret[1] = (-11*f(x) + 18*f(x+interval) - 9*f(x+2*interval) +
	 2*f(x+3*interval))/(6.*interval) + ((-2*f(x) + 5*f(x+interval) - 
	 	4*f(x+2*interval) + f(x+3*interval))*x)/pow(interval,2) + ((-f(x) + 3*f(x+interval) - 
	 	3*f(x+2*interval) + f(x+3*interval))*pow(x,2))/(2.*pow(interval,3));

	ret[2] = (f(x) - (5*f(x+interval))/2. + 2*f(x+2*interval) - 
		f(x+3*interval)/2.)/pow(interval,2) + ((f(x) - 3*f(x+interval) + 3*f(x+2*interval) - 
			f(x+3*interval))*x)/(2.*pow(interval,3));
	
	ret[3] = (-f(x) + 3*f(x+interval) - 3*f(x+2*interval) + f(x+3*interval))/(6.*pow(interval,3));
	return ret;
}


// use stdVec to store coef, I think it's better than using a struct...



/*
	Main Dish!
*/

/*
	Arguments:
		1. marginal transition expansion coef \alpha_{k, i}^{(l)}, vector of length 4*nbasis
		2. marginal moment expansion coef \beta_{k, i}^{(l)}, square matrix of size (4*nbasis)^2
		calculate 1 & 2 using
		3. (i-1) th y0
		4. i th y
		5. par = {mu, kappa, theta, xi, rho};
*/
void BasisCoef(stdVec& alpha, //stdMat& beta, 
	double y0, double y, const stdVec& par) {
	// determine the interpolation nodes, from the stationary gamma distribution
	// same as method in InitialMoments.hpp
	double stLogMean = log(par[2]);
	double stLogSD = log(sqrt( par[2]*par[3]*par[3] / (2*par[1]) ));
	// create the (nbasis+1) nodes we want to interpolate
	stdVec interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	for(int i = 0; i < 22; ++i) 
		interNodes[i] = exp(stLogMean - interNodes[i]*stLogSD) / (par[3]*par[3]/(2*par[1]));
	/*
		prepare for the integral
	*/
	// conditional mean and sd
	/*
	double mu_c = cond_v_mean(par, v0);
	double var_c = cond_v_var(par, v0);
	double sd_c = sqrt(cond_v_var(par, v0));
	NormalPdf pdf_cond(mu_c, var_c);
	Normalcdf cdf_cond(mu_c, var_c);
	stdVec p_c(nbasis+1);
	stdVec phi_c(nbasis+1);
	transform(interNodes.begin(), interNodes.end(), p_c.begin(), (pdf_cond());
	transform(interNodes.begin(), interNodes.end(), p_c.begin(), (cdf_cond());
	*/
	
	int i, j;
	// first calculate alpha easily 
	MarTranDensity p1_y(y, y0); // the marginal density functional object
	double stepSize;
	stdVec tmp(4);
	auto p_alpha = alpha.begin();
	for(i = 0; i < nbasis; ++i) {
		stepSize = (interNodes[i+1] - interNodes[i])/3.0;
		tmp = interCoef(p1_y, interNodes[i], stepSize);
		copy(tmp.begin(), tmp.end(), p_alpha + 4*i);
	}
	
	/* 
		calculate beta matrix (row first):
		- first traverse the k nodes, to get the B_k at each interval
		- then for each B_k, expand it to each j nodes
		Implementation:
			1. first initialize the four integral functional object I,
				which should be f(x_i) = p1_y * I_{1, 2, 3, 4},
				which should be calculated recursively.
			2. use transform to get the values of 4 Intgeral at each interval
				as well as p1_y
	*/
}





