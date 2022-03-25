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
#include <thread>

using Complex = std::complex<double>;
using vectorReal    = std::vector<double>;
using vectorComplex = std::vector<std::complex<double>>;
using matrixComplex = std::vector<vectorComplex>;
using chemical_potential = double;

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
extern const int mesh;

extern const std::string axises[];
extern const int valleys;

extern std::vector<matrixComplex> vT;
extern std::vector<std::vector<matrixComplex>> vL;
extern vectorReal ET;
extern std::vector<vectorReal> EL;

extern matrixComplex Hamiltonian;
void initialize();
matrixComplex set_T(double k[3]);
matrixComplex set_L(int valley, double k[3]);
vectorReal diagonalize_N(matrixComplex A);

extern "C" {
    void zheev_ ( const char& JOBZ, const char& UPLO,
        const int& N, Complex** A, const int& LDA,
        double* W, Complex* WORK, const int& LWORK,
        double* RWORK,
        int& INFO, int JOBZlen, int UPLOlen );
};

struct vector3 {
    double vec[3];
};

using kpoint = vector3;
using velocity= vector3;

struct fermi_surface {
    double e;
    std::vector<kpoint> kset;
    std::vector<velocity> vset;
};

struct band {
    int index;
    chemical_potential mu_min, mu_max, dmu;
    int mu_mesh;
    std::vector<fermi_surface> fs;
};

struct sys {
    std::vector<band> bands;
};

fermi_surface get_fermi_suraceT(int band_index, chemical_potential mu);
int fermi_surface_write(fermi_surface fs, std::string filename);
velocity get_velocity(int band_index, chemical_potential mu, kpoint k);

void get_band_T(band& b, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh);

sys get_T();
void sys_write(sys s);


#endif // PARAMETERS_HPP
