#include "Coef.hpp"

#include <iostream>
#include <iomanip>

int main() {
    stdVec alpha(84);
    double y0 = 4.6052, y1 = 4.6195;
    stdVec par = {0.05, 2.0, 0.2, 0.25, -0.8};
    BasisCoef(alpha, y0, y1, par);
    cout << setprecision(16) << setw(16) << endl;
    for(int i = 0; i < 84; ++i) {
        cout << alpha[i] << '\t';
        if((i+1) % 4 == 0)
            cout << endl;
    }
    cout << alpha.size() << endl;
}
