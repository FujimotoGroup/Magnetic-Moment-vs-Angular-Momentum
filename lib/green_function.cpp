#include "parameters.hpp"

matrixComplex get_inverse(matrixComplex A) { // {{{
    const int N = A.size();

    int lda = N;
    int lwork = 16*N;
    Complex AT[N*N];
    int ipiv[N];
    Complex work[lwork];

    int info;

    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            AT[N*i+j] = A[i][j];
        }
    }

    zgetrf_(N, N, AT, lda, ipiv, info);
    if (info != 0) {
        std::cout << info << std::endl;
    }

    zgetri_(N, AT, lda, ipiv, work, lwork, info);
    if (info != 0) {
        std::cout << info << std::endl;
    }

    matrixComplex A_inv(N, vectorComplex(N, 0e0));

    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A_inv[i][j] = AT[N*i+j];
        }
    }

    return A_inv;
}; // }}}

Green_function get_green_function_T(Complex ene, kpoint k) { // {{{
    Green_function G(bandsT, vectorComplex(bandsT, 0e0));

    G = set_T(k.vec);
    for(int i=0; i<bandsT; i++) {
        G[i][i] = ene - G[i][i];
        for(int j=i+1; j<bandsT; j++) {
            G[i][j] = - G[i][j];
            G[j][i] = - G[j][i];
        }
    }

//    matrixComplex Identity(bandsT, vectorComplex(bandsT, 0e0));
//    matrixComplex inv = get_inverse(G);
//    for(int i=0; i<bandsT; i++) {
//        for(int j=0; j<bandsT; j++) {
//            std::cout << inv[i][j];
//        }
//        std::cout << std::endl;
//    }
//    for(int i=0; i<bandsT; i++) {
//        for(int j=0; j<bandsT; j++) {
//            for(int l=0; l<bandsT; l++) {
//                Identity[i][j] += inv[i][l]*G[l][j];
//            }
//            std::cout << std::setprecision(5) << Identity[i][j] << ", ";
//        }
//        std::cout << std::endl;
//    }

    G = get_inverse(G);
    return G;
}; // }}}

Green_function get_green_function_L(Complex ene, int valley, kpoint k) { // {{{
    Green_function G(bandsL, vectorComplex(bandsL, 0e0));

    G = set_L(valley, k.vec);
    for(int i=0; i<bandsL; i++) {
        G[i][i] = ene - G[i][i];
        for(int j=i+1; j<bandsL; j++) {
            G[i][j] = - G[i][j];
            G[j][i] = - G[j][i];
        }
    }

    G = get_inverse(G);
    return G;
}; // }}}

Green_function get_full_green_function_T(Energy ene, kpoint k, Self_energy se) { // {{{
    Green_function G(bandsT, vectorComplex(bandsT, 0e0));

    G = set_T(k.vec);
    for(int i=0; i<bandsT; i++) {
        G[i][i] = ene - G[i][i] - se[i][i];
        for(int j=i+1; j<bandsT; j++) {
            G[i][j] = - G[i][j] - se[i][j];
            G[j][i] = - G[j][i] - se[j][i];
        }
    }

    G = get_inverse(G);
    return G;
}; // }}}

Green_function get_full_green_function_L(Energy ene, int valley, kpoint k, Self_energy se) { // {{{
    Green_function G(bandsL, vectorComplex(bandsL, 0e0));

    G = set_L(valley, k.vec);
    for(int i=0; i<bandsL; i++) {
        G[i][i] = ene - G[i][i] - se[i][i];
        for(int j=i+1; j<bandsL; j++) {
            G[i][j] = - G[i][j] - se[i][j];
            G[j][i] = - G[j][i] - se[j][i];
        }
    }

    G = get_inverse(G);
    return G;
}; // }}}
