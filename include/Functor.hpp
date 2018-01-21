#pragma once

#include <cmath>
#include "Jump.hpp"

// scalar 1d function
class Functor {
public:
    Functor() {}
    virtual ~Functor() {}
    virtual double operator() (double x) const = 0;
};

// taking two arguments
class Functor2 {
public:
    Functor2() {}
    virtual ~Functor2() {}
    virtual double operator() (double x, double y) const = 0;
};

class NormalPdf: public Functor{
public:
    double mean;
    double variance;
    
    NormalPdf() {}
    explicit NormalPdf(double mu, double var):
        mean(mu), variance(var) {}
    ~NormalPdf() {}
    
	double operator() (double x) const { 
        return exp(-pow(x-mean, 2)/(2*variance)) / sqrt(2*Pi*variance);
    }
};

class NormalCdf: public Functor{
public:
	double mean;
	double variance;
	
    NormalCdf() {}
    explicit NormalCdf(double mu, double var):
        mean(mu), variance(var) {}
    ~NormalCdf() {}
    
	double operator() (double x) const { 
        double z = (x - mean)/sqrt(variance);
		return (1 - 0.5*erfc(z/sqrt(2)));
    }
};

class MarTranDensity: public Functor {
public:
	double y, y0;
	stdVec par;	
	
	MarTranDensity() {}
	explicit MarTranDensity(double y_, double y0_, const stdVec& par_):
		y(y_), y0(y0_), par(par_) {}
	~MarTranDensity() {}
	
    double operator() (double v0) const { 
		double a = pow(y - tran_y_mean(par, y0, v0), 2);
		double b = tran_y_var(par, v0);
		return (exp(-a/(2*b)) / sqrt(2*Pi*b));
    }
};

class CondPdf: public Functor2{
public:
	double y, y0;
	stdVec par;
	
	CondPdf() {}
	explicit CondPdf(double y_, double y0_, const stdVec& par_):
		y(y_), y0(y0_), par(par_) {}
	~CondPdf() {}
	
	double operator() (double v, double v0) const { 
		double mu = v0 + par[1]*(par[2] - v0)*Delta + par[3]*par[4]*(y - (y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s));
		double var = pow(par[3], 2) *v0*Delta*(1 - pow(par[4], 2));
		return exp(-pow(v-mu, 2)/(2*var)) / sqrt(2*Pi*var);
	}
};

class CondCdf: public Functor2{
public:
	double y, y0;
	stdVec par;
	
	CondCdf() {}
	explicit CondCdf(double y_, double y0_, const stdVec& par_):
		y(y_), y0(y0_), par(par_) {}
	~CondCdf() {}
	
	double operator() (double v, double v0) const { 
		double mu = v0 + par[1]*(par[2] - v0)*Delta + par[3]*par[4]*(y - (y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s));
		double var = pow(par[3], 2) *v0*Delta*(1 - pow(par[4], 2));
		double x = (v-mu)/sqrt(var);
		return (1 - 0.5*erfc(x/sqrt(2)));
	}
};

class Integral0: public Functor{
public:
	double y, y0, v1, v2;
	stdVec par;
	
	Integral0() {}
	explicit Integral0(double y_, double y0_, double v1_, double v2_, const stdVec& par_):
		y(y_), y0(y0_), v1(v1_), v2(v2_), par(par_) {}
	~Integral0() {}
	
	double operator() (double v0) const { 
		MarTranDensity prePy(y, y0, par);
		CondCdf cdf_c(y, y0, par);
		return prePy(v0) * (cdf_c(v2, v0) - cdf_c(v1, v0));
	}
};

class Integral1: public Functor{
public:
	double y, y0, v1, v2;
	stdVec par;
	
	Integral1() {}
	explicit Integral1(double y_, double y0_, double v1_, double v2_, const stdVec& par_):
		y(y_), y0(y0_), v1(v1_), v2(v2_), par(par_) {}
	~Integral1() {}
	
	double operator() (double v0) const { 
		MarTranDensity prePy(y, y0, par);
		CondPdf pdf_c(y, y0, par);
		Integral0 I0(y, y0, v1, v2, par);
		double mu_c = (v0 + par[1]*(par[2] - v0)*Delta) + par[3]*par[4]*(y - (y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s));
		double var_c = (pow(par[3], 2)*v0*Delta * (1 - pow(par[4], 2)));
		return I0(v0)*mu_c  - prePy(v0)*var_c*(pdf_c(v2, v0) - pdf_c(v1, v0));
	}
};

class Integral2: public Functor{
public:	
	double y, y0, v1, v2;
	stdVec par;
	
	Integral2() {}
	explicit Integral2(double y_, double y0_, double v1_, double v2_, const stdVec& par_):
		y(y_), y0(y0_), v1(v1_), v2(v2_), par(par_) {}
	~Integral2() {}
	
	double operator() (double v0) const { 
		MarTranDensity prePy(y, y0, par);
		CondPdf pdf_c(y, y0, par);
		Integral0 I0(y, y0, v1, v2, par);
		Integral1 I1(y, y0, v1, v2, par);
		double mu_c = (v0 + par[1]*(par[2] - v0)*Delta) + par[3]*par[4]*(y - (y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s));
		double var_c = (pow(par[3], 2)*v0*Delta * (1 - pow(par[4], 2)));
		return I0(v0)*var_c + I1(v0)*mu_c - prePy(v0)*var_c*(pdf_c(v2, v0)*v2 - pdf_c(v1, v0)*v1);
	}
}; 

class Integral3: public Functor{
public:
	double y, y0, v1, v2;
	stdVec par;
	
	Integral3() {}
	explicit Integral3(double y_, double y0_, double v1_, double v2_, const stdVec& par_):
		y(y_), y0(y0_), v1(v1_), v2(v2_), par(par_) {}
	~Integral3() {}
	
	double operator() (double v0) const { 
		MarTranDensity prePy(y, y0, par);
		CondPdf pdf_c(y, y0, par);
		Integral1 I1(y, y0, v1, v2, par);
		Integral2 I2(y, y0, v1, v2, par);
		double mu_c = (v0 + par[1]*(par[2] - v0)*Delta) + par[3]*par[4]*(y - (y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s));
		double var_c = (pow(par[3], 2)*v0*Delta * (1 - pow(par[4], 2)));
		return 2*I1(v0)*var_c + I2(v0)*mu_c - prePy(v0)*var_c*(pdf_c(v2, v0)*v2*v2 - pdf_c(v1, v0)*v1*v1);
	}
}; 
