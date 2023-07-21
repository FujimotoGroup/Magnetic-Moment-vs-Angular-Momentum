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

matrixComplex set_T(chemical_potential mu, double k[3]) { // {{{
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
        value[i][i] = ET[i] + ene - mu;
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
//    const tensor2Complex sigma = { { {0e0, 1e0}, {1e0, 0e0} },
//                                   { {0e0,- zi}, { zi, 0e0} },
//                                   { {1e0, 0e0}, {0e0,-1e0} } };
//    for(int l=0; l<space_dim; l++) {
//        for(int i=0; i<2; i++) {
//            for(int j=0; j<2; j++) {
//                value[i][2+j] +=  zi*hk*k[l]*sigma[l][i][j];
//                value[2+i][j] += -zi*hk*k[l]*sigma[l][i][j];
//            }
//        }
//    }
//    value[0][0] =  m;
//    value[1][1] =  m;
//    value[2][2] = -m;
//    value[3][3] = -m;
    return value;
}; // }}}

matrixComplex set_L(int valley, chemical_potential mu, double k[3]) { // {{{
    matrixComplex value(bandsL, vectorComplex(bandsL, 0e0));
    double ene = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];
    ene = ene * (charge*hbar*hbar/angstrom/angstrom/mass*5e-1);
    ene = 0e0;
    for(int i=0; i<bandsL; i++) {
        value[i][i] = EL[valley][i] + ene - mu;
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
//    const tensor2Complex sigma = { { {0e0, 1e0}, {1e0, 0e0} },
//                                   { {0e0,- zi}, { zi, 0e0} },
//                                   { {1e0, 0e0}, {0e0,-1e0} } };
//    for(int l=0; l<space_dim; l++) {
//        for(int i=0; i<2; i++) {
//            for(int j=0; j<2; j++) {
//                value[i][2+j] +=  zi*hk*k[l]*sigma[l][i][j];
//                value[2+i][j] += -zi*hk*k[l]*sigma[l][i][j];
//            }
//        }
//    }
//    value[0][0] =  m;
//    value[1][1] =  m;
//    value[2][2] = -m;
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

    zheev_( 'N', 'L', N, (Complex**)U, N, E, cwork, 4*N, rwork, info, 1, 1 );

    vectorReal res(N);
    for(int i=0; i<N; i++ ){
        res[i] = E[i];
    }

    return res;
}; // }}}

diag_set diagonalize_V(matrixComplex A) { // {{{
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

    diag_set eigen;
    eigen.values.resize(N);
    eigen.vectors.resize(N);
    for(int i=0; i<N; i++){
        eigen.vectors[i].resize(N);
    }

    for(int i=0; i<N; i++){
        eigen.values[i] = E[i];

        eigen.vectors[i][i] = U[i][i];
        for(int j=i+1; j<N; j++){
            eigen.vectors[i][j] = U[j][i];
            eigen.vectors[j][i] = U[i][j];
        }
    }

    return eigen;
}; // }}}
