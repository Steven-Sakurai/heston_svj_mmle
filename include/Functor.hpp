#pragma once

#include <cmath>

#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

#include <vector>
typedef std::vector< double > stdVec;

#include "Jump.h"

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
	
	MarTranDensity() {}
	explicit MarTranDensity(double y_, double y0_):
		y(y_), y0(y0_) {}
	~MarTranDensity() {}
	
    double operator() (double v0) const { 
		double a = pow(-y + y0 + (par[0]-0.5*v0)*Delta + lambda*Delta*mu_s, 2);
		double b = v0*Delta + lambda*Delta*(var_s + pow(mu_s, 2));
		return (exp(-a/(2*b)) / sqrt(2*Pi*b));
    }
};


class CondPdf: public Functor2{
public:
	double y, y0;
	
	CondPdf() {}
	explicit CondPdf(double y_, double y0_):
		y(y_), y0(y0_) {}
	~CondPdf() {}
	
	double operator() (double v, double v0) const { 
		
	}
};

class CondCdf: public Functor2{
public:
	double y, y0;
	
	CondCdf() {}
	explicit CondCdf(double y_, double y0_):
		y(y_), y0(y0_) {}
	~CondCdf() {}
	
	double operator() (double x) const { 
	}
};
