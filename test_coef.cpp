#include "Coef.hpp"
#include <iostream>

class Square: public Functor {
public:
    Square() {}
    ~Square() {}
    double operator() (double x) const { 
        return x*x;
    }
};

int main() {
    NormalPdf my_fun(5.0, 1.0);
    cout << my_fun(5.0) << endl;
    Square my_square;
    stdVec coef = interCoef(my_square, -2, 1);
    for(int i = 0; i < 4; ++i)
        cout << coef[i] << endl;
    return 0;
}
