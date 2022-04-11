#include "parameters.hpp"

matrixComplex product(matrixComplex A, matrixComplex B) { // {{{
    const int N = A.size();
    matrixComplex C(N, vectorComplex(N, 0e0));
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            for(int k=0; k<N; k++) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}; // }}}
