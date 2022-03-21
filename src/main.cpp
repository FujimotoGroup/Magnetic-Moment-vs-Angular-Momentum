#include "parameters.hpp"

int main(){
    initialize();

    int n = 50;
    std::string file_name = "./dat/T.dat";
    std::ofstream file(file_name);
    for (int i=-n; i<=n; i++) {
        double kx = 5e-2*g0*double(i)/double(n);
        vectorReal k = {kx,0e0,0e0};
        matrixComplex H_T = set_T(k);

        vectorReal E_T = diagonalize_N(H_T);

        file << std::scientific << kx << ", ";
        for(int i=0; i<bandsT; i++) {
            file << std::scientific << E_T[i];
            if (i != bandsT-1) file << ", ";
        }
        file << std::endl;

    };

}
