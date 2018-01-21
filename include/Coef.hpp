#pragma once
#include <cmath>
#include <algorithm>

using namespace std;
#include "Functor.hpp"

// Originally struct intercoefficient InterCoef_cpp()
// this part is just copy and paste since there's no elegant way to do this...
// test this function in `test_coef.cpp`
stdVec interCoef(const Functor& f, double x, double interval) {
	// use stdVec to store coef, I think it's better than using a struct...
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


/*
	Main Dish!
		Arguments:
		1. marginal transition expansion coef \alpha_{k, i}^{(l)}, vector of length 4*nbasis
		2. marginal moment expansion coef \beta_{k, i}^{(l)}, square matrix of size (4*nbasis)^2
		calculate 1 & 2 using
		3. (i-1) th y0
		4. i th y
		5. parameters
*/
void BasisCoef(stdVec& alpha, stdMat& beta, double y, double y0, const stdVec& par) {
	// determine the interpolation nodes, from the stationary gamma distribution
	// same as method in InitialMoments.hpp
	double stLogMean = log(par[2]);
	double stLogSD = log(sqrt( par[2]*par[3]*par[3] / (2*par[1]) ));
	// create the (nbasis+1) nodes we want to interpolate
	stdVec interNodes = {-0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35,
		 -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.3, 0.4};
	int i;
	for(int i = 0; i < 22; ++i) 
		interNodes[i] = exp(stLogMean - interNodes[i]*stLogSD);
	// first calculate alpha easily 
	MarTranDensity p1_y(y, y0, par); // the marginal density functional object, i.e. prestfun
	double stepSize;
	stdVec tmp(4);
	for(i = 0; i < nbasis; ++i) {
		stepSize = (interNodes[i+1] - interNodes[i])/3.0;
		tmp = interCoef(p1_y, interNodes[i], stepSize);
		alpha[i] = tmp[0];
		alpha[i+nbasis] = tmp[1];
		alpha[i+2*nbasis] = tmp[2];
		alpha[i+3*nbasis] = tmp[3];
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
	
	int k, j;
	// traverse v first then v0, row first
	for(k = 0; k < nbasis; ++k) {
		// initialize the four functor
		Integral0 I0(y, y0, interNodes[k], interNodes[k+1], par);
		Integral1 I1(y, y0, interNodes[k], interNodes[k+1], par);
		Integral2 I2(y, y0, interNodes[k], interNodes[k+1], par);
		Integral3 I3(y, y0, interNodes[k], interNodes[k+1], par);
		
		for(j = 0; j < nbasis; ++j) {
			if(abs(k-j) > 3)
				continue;
			stepSize = (interNodes[j+1] - interNodes[j])/3.0;
			
			tmp = interCoef(I0, interNodes[j], stepSize);
			beta[k][j] = tmp[0];          beta[k][j+nbasis] = tmp[1];
			beta[k][j+2*nbasis] = tmp[2]; beta[k][j+3*nbasis] = tmp[3];
			//beta[4*k+0][4*j+0] = tmp[0]; beta[4*k+0][4*j+1] = tmp[1];
			//beta[4*k+0][4*j+2] = tmp[2]; beta[4*k+0][4*j+3] = tmp[3];

			tmp = interCoef(I1, interNodes[j], stepSize);
			beta[k+nbasis][j] = tmp[0];          beta[k+nbasis][j+nbasis] = tmp[1];
			beta[k+nbasis][j+2*nbasis] = tmp[2]; beta[k+nbasis][j+3*nbasis] = tmp[3];
			//beta[4*k+1][4*j+0] = tmp[0]; beta[4*k+1][4*j+1] = tmp[1];
			//beta[4*k+1][4*j+2] = tmp[2]; beta[4*k+1][4*j+3] = tmp[3];
			
			tmp = interCoef(I2, interNodes[j], stepSize);
			beta[k+2*nbasis][j] = tmp[0];          beta[k+2*nbasis][j+nbasis] = tmp[1];
			beta[k+2*nbasis][j+2*nbasis] = tmp[2]; beta[k+2*nbasis][j+3*nbasis] = tmp[3];
			//beta[4*k+2][4*j+0] = tmp[0]; beta[4*k+2][4*j+1] = tmp[1];
			//beta[4*k+2][4*j+2] = tmp[2]; beta[4*k+2][4*j+3] = tmp[3];
			
			tmp = interCoef(I3, interNodes[j], stepSize);
			beta[k+3*nbasis][j] = tmp[0];          beta[k+3*nbasis][j+nbasis] = tmp[1];
			beta[k+3*nbasis][j+2*nbasis] = tmp[2]; beta[k+3*nbasis][j+3*nbasis] = tmp[3];
			//beta[4*k+3][4*j+0] = tmp[0]; beta[4*k+3][4*j+1] = tmp[1];
			//beta[4*k+3][4*j+2] = tmp[2]; beta[4*k+3][4*j+3] = tmp[3];
		}
	}
}





