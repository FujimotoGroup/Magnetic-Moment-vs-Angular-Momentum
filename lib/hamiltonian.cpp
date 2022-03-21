#include "parameters.hpp"

matrixComplex set_T(vectorReal k) {
    matrixComplex value(bandsT, vectorComplex(bandsT, 0e0));
    for(int i=0; i<bandsT; i++) {
        double ene = std::inner_product(k.begin(), k.end(), k.begin(), 0e0);
        ene = ene * (charge*hbar*hbar/angstrom/angstrom/mass);
        value[i][i] = ET[i] + ene;
        for(int j=i+1; j<bandsT; j++) {
            for(int axis=0; axis<space_dim; axis++) {
                value[i][j] += hbar/angstrom*k[axis]*vT[axis][i][j];
                value[j][i] += hbar/angstrom*k[axis]*vT[axis][j][i];
            }
        }
    }
    return value;
};


extern "C" {
  void zheev_ ( const char& JOBZ, const char& UPLO,
		const int& N, Complex** A, const int& LDA,
		double* W, Complex* WORK, const int& LWORK,
		double* RWORK,
		int& INFO, int JOBZlen, int UPLOlen );
};

template <int N> int zheev( Complex H[N][N], Complex U[N][N], double E[N] )
{
  int i, j, info;




  return info;
}

vectorReal diagonalize_N(matrixComplex A) {
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
};
