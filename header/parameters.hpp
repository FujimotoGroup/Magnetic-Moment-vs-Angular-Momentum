#ifndef   PARAMETERS_HPP
#define   PARAMETERS_HPP

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <numeric>

using Complex = std::complex<double>;
using vectorReal    = std::vector<double>;
using vectorComplex = std::vector<std::complex<double>>;
using matrixComplex = std::vector<vectorComplex>;

extern const double angstrom;
extern const double hbar;
extern const double mass;
extern const double charge;

extern const int bands;
extern const int bandsT;
extern const int bandsL;
extern const int lowest_band_T;
extern const int lowest_band_L;

extern const double a;
extern const double c;
extern const double g0;

extern const vectorReal b1;
extern const vectorReal b2;
extern const vectorReal b3;

extern const int space_dim;
extern const int spin_dim;
extern const double eps_phys;

extern const double cutoff;

extern const std::string axises[];
extern const int valleys;

extern std::vector<matrixComplex> vT;
extern std::vector<std::vector<matrixComplex>> vL;
extern vectorReal ET;
extern std::vector<vectorReal> EL;

extern matrixComplex Hamiltonian;
void initialize();
matrixComplex set_T(vectorReal k);
matrixComplex set_L(int valley, vectorReal k);
vectorReal diagonalize_N(matrixComplex A);

extern "C" {
    #include "gts.h"
}

extern int band_index;
extern int valley_index;
void dispersionT(gdouble ** f, GtsCartesianGrid g, guint k, gpointer data);
void dispersionL(gdouble ** f, GtsCartesianGrid g, guint k, gpointer data);

extern "C" {
    void zheev_ ( const char& JOBZ, const char& UPLO,
        const int& N, Complex** A, const int& LDA,
        double* W, Complex* WORK, const int& LWORK,
        double* RWORK,
        int& INFO, int JOBZlen, int UPLOlen );
};

int test();

#endif // PARAMETERS_HPP
