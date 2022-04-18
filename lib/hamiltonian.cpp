#include "parameters.hpp"

matrixComplex set_T(double k[3]) { // {{{
    matrixComplex value(bandsT, vectorComplex(bandsT, 0e0));
    if ( (std::abs(k[0]) > cutoff) | (std::abs(k[1]) > cutoff) | (std::abs(k[2]) > cutoff) ) {
        for(int i=0; i<bandsT; i++) {
            value[i][i] = -10e0;
        }
        return value;
    }
    double ene = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    ene = ene * (charge*hbar*hbar/angstrom/angstrom/mass*5e-1);
    ene = 0e0;
    for(int i=0; i<bandsT; i++) {
        value[i][i] = ET[i] + ene;
        for(int j=i+1; j<bandsT; j++) {
            for(int axis=0; axis<space_dim; axis++) {
                value[i][j] += hbar/angstrom*v0*k[axis]*vT[axis][i][j];
                value[j][i] += hbar/angstrom*v0*k[axis]*vT[axis][j][i];
            }
        }
    }
    return value;
}; // }}}

matrixComplex set_L(int valley, double k[3]) { // {{{
    matrixComplex value(bandsL, vectorComplex(bandsL, 0e0));
    double ene = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    ene = ene * (charge*hbar*hbar/angstrom/angstrom/mass*5e-1);
    ene = 0e0;
    for(int i=0; i<bandsL; i++) {
        value[i][i] = EL[valley][i] + ene;
        for(int j=i+1; j<bandsL; j++) {
            for(int axis=0; axis<space_dim; axis++) {
                value[i][j] += hbar/angstrom*v0*k[axis]*vL[valley][axis][i][j];
                value[j][i] += hbar/angstrom*v0*k[axis]*vL[valley][axis][j][i];
            }
        }
    }
//    double m = (std::abs(EL[valley][2] - EL[valley][0]))*5e-1;
//    double vF = 1e6;
//    double hk= hbar/angstrom*v0;
//    value[0][0] = m;
//    value[0][1] = 0e0;
//    value[0][2] = zi*hk*k[2];
//    value[0][3] = zi*hk*k[0] + hk*k[1];
//    value[1][0] = 0e0;
//    value[1][1] = m;
//    value[1][2] =-zi*hk*k[0] - hk*k[1];
//    value[1][3] =-zi*hk*k[2];
//    value[2][0] =-zi*hk*k[2];
//    value[2][1] =-zi*hk*k[0] - hk*k[1];
//    value[2][2] = -m;
//    value[2][3] = 0e0;
//    value[3][0] =-zi*hk*k[0] + hk*k[1];
//    value[3][1] = zi*hk*k[2];
//    value[3][2] = 0e0;
//    value[3][3] = -m;
    return value;
}; // }}}

vectorReal diagonalize_N(matrixComplex A) { // {{{
    const int N = A.size();
    Complex H[N][N];
    Complex U[N][N];

    for(int i=0; i<N; i++ ){
        for(int j=0; j<N; j++ ){
            H[i][j] = A[i][j];
            U[i][j] = A[j][i];
        }
    }

    double E[N];

    int info;
    Complex cwork[4*N];
    double  rwork[4*N];

    zheev_( 'V', 'L', N, (Complex**)U, N, E, cwork, 4*N, rwork, info, 1, 1 );

    for(int i=0; i<N; i++ ){
        for(int j=i+1; j<N; j++ ){
            Complex temp;
            temp = U[i][j];
            U[i][j] = U[j][i];
            U[j][i] = temp;
        }
    }

    vectorReal res(N);
    for(int i=0; i<N; i++ ){
        res[i] = E[i];
    }

    return res;
}; // }}}
