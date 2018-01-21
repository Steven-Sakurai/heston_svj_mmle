#pragma once

/*
	This file is for controlling global parameters:
		- Heston parameters `par`
		- jump parameters
		......
*/

#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

#define Delta 1.0/252.0
#define nbasis 21

extern double lambda = 0;
extern double mu_s = -0.3;
extern double var_s = 1.0;

#include <vector>
#include <cmath>
using namespace std;

typedef vector< double > stdVec;
typedef vector< vector< double > > stdMat;

// extern stdVec par = {0.05, 2.0, 0.2, 0.25, -0.8};

/*
	Some support function to save coding time
*/

/*
	The full one-step transition density parameters
*/
inline double tran_y_mean(const stdVec& par, double y0, double v0) {
	return (y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s);
}

inline double tran_y_var(const stdVec& par, double v0) {
	return (v0*Delta + lambda*Delta*(var_s + pow(mu_s, 2)));
}

inline double tran_v_mean(const stdVec& par, double v0) {
	return (v0 + par[1]*(par[2] - v0)*Delta);
}

inline double tran_v_var(const stdVec& par, double v0) {
	return (pow(par[3], 2)*v0*Delta);
}

/*
	The one-step transition density decomposition: p(y, v) = p(y) * p(v|y) 
*/
double cond_v_var(const stdVec& par, double v0)
{
	return (pow(par[3], 2)*v0*Delta * (1 - pow(par[4], 2)));
}

double cond_v_mean(const stdVec& par, double y, double y0, double v0)
{
	return ((v0 + par[1]*(par[2] - v0)*Delta) + par[3]*par[4]*(y - (y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s)));
}
