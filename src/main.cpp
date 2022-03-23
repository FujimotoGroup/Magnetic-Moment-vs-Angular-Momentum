#include "parameters.hpp"

int main(){
    initialize();

// output kx-dependence @ T {{{
    int n = 50;
    std::string file_name = "./dat/EigenValue_T_"+std::to_string(bandsT)+"bands.dat";
    std::ofstream file(file_name);
    for (int i=-n; i<=n; i++) {
        double kx = cutoff*double(i)/double(2*n);
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
/// }}}
// output kx-dependence @ Ls {{{
    for(int valley=0; valley<valleys; valley++) {
        std::string file_name = "./dat/EigenValue_L"+std::to_string(valley+1)+"_"+std::to_string(bandsL)+"bands.dat";
        std::ofstream file(file_name);
        for (int i=-n; i<=n; i++) {
            double kx = cutoff*double(i)/double(2*n);
            vectorReal k = {kx,0e0,0e0};
            matrixComplex H_L = set_L(valley, k);

            vectorReal E_L = diagonalize_N(H_L);

            file << std::scientific << kx << ", ";
            for(int i=0; i<bandsL; i++) {
                file << std::scientific << E_L[i];
                if (i != bandsL-1) file << ", ";
            }
            file << std::endl;

        };
    };
/// }}}

    

}
