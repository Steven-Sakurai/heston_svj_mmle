#pragma once

#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

#define Delta 1.0/252.0  // one trading day, interval between obs
#define nbasis 21

extern const double lambda = 0;
extern const double mu_s = -0.3;
extern const double var_s = 1.0;

#include <vector>

typedef vector< double > stdVec;

stdVec par = {0.05, 2.0, 0.2, 0.25, -0.8};

/*
	Some support function to save coding time
*/

/*
	The full one-step transition density parameters
*/
inline double tran_y_mean(const stdVec& par, double y0, double v0) {
	return y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s;
}

inline double tran_y_var(const stdVec& par, double v0) {
	return v0*Delta + lambda*Delta*(var_s + pow(mu_s, 2));
}

inline double tran_v_mean(const stdVec& par, double v0) {
	return v0 + par[1]*(par[2] - v0)*Delta;
}

inline double tran_v_var(const stdVec& par, double v0) {
	return pow(par[3], 2)*v0*Delta;
}

inline double tran_cov(const stdVec& par, double v0) {
	return par[4]*par[3]*v0*Delta; 
}

inline double tran_rho(const stdVec& par, double v0) {
	return tran_cov(par, v0) / sqrt( tran_v_var(par, v0)*tran_y_var(par, v0) );
}

/*
	The one-step transition density decomposition: p(y, v) = p(y) * p(v|y) 
*/
inline double cond_v_var(const stdVec& par, double v0)
{
	return tran_v_var(par, v0) * (1 - pow(tran_rho(par, v0), 2));
}

inline double cond_v_mean(const stdVec& par, double y, double y0, double v0)
{
	return tran_v_mean(par, v0) + sqrt( tran_v_var(par, v0)/tran_y_var(par, v0) ) * 
		tran_rho(par, v0) * (y - tran_y_mean(par, y0, v0));
}
