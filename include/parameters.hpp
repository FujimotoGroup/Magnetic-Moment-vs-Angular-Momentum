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
#include <mutex>

extern std::mutex mtx;

using Complex = std::complex<double>;
using vectorReal    = std::vector<double>;
using vectorComplex = std::vector<std::complex<double>>;
using matrixComplex = std::vector<vectorComplex>;
using tensor2Complex = std::vector<matrixComplex>;
using tensor3Complex = std::vector<tensor2Complex>;
using matrixReal = std::vector<vectorReal>;
using chemical_potential = double;

extern const double angstrom;
extern const double hbar;
extern const double mass;
extern const double charge;

extern const double v0;

extern const double eps;
extern const double pi;
extern const Complex zi;

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

extern vectorReal kT;
extern vectorReal kL[3];

extern const int space_dim;
extern const int spin_dim;
extern const double eps_phys;
extern const double eps_num;

extern const double cutoff;
extern double dk[3];
extern const int k_mesh;
extern const int k_mesh_more;
extern const int mu_mesh;

extern const std::string axises[];
extern const int valleys;

extern std::vector<matrixComplex> vT;
extern std::vector<std::vector<matrixComplex>> vL;
extern vectorReal ET;
extern std::vector<vectorReal> EL;

extern std::vector<matrixComplex> mu_s_T;
extern std::vector<std::vector<matrixComplex>> mu_s_L;

extern std::vector<std::vector<matrixComplex>> v_s_T;
extern std::vector<std::vector<std::vector<matrixComplex>>> v_s_L;

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
    void zgetrf_ ( const int& N, const int& M, Complex* A, const int& LDA, int* IPIV, int& INFO);
    void zgetri_ ( const int& N, Complex* A, const int& LDA, int* IPIV, Complex* WORK, const int& LWORK, int& INFO);
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

struct face {
    int face[3];
    double dS;
    double center[3];
    double normal[3];
};

struct triangles {
    double ene;
    std::vector<vector3> vertexes;
    std::vector<vector3> normals;
    std::vector<face> faces;
};

triangles get_triangles(fermi_surface fs);
fermi_surface get_fermi_surace_T(int band_index, chemical_potential mu);
velocity get_velocity_T(int band_index, chemical_potential mu, kpoint k);
fermi_surface get_fermi_surace_L(int valley, int band_index, chemical_potential mu);
velocity get_velocity_L(int band_index, chemical_potential mu, kpoint k);
int fermi_surface_write(fermi_surface fs, std::string filename);
int triangles_write_T(triangles tri, std::string filename);
int triangles_write_L(triangles tri, std::string filename, int valley);
triangles get_triangles_T(int band_index, chemical_potential mu);
triangles get_triangles_L(int valley, int band_index, chemical_potential mu);

void init(double& value);
void init(Complex& value);
void init(vectorComplex& value);
void init(matrixComplex& value);

double        add(double value, double a);
Complex       add(Complex value, Complex a);
vectorComplex add(vectorComplex value, vectorComplex a);
matrixComplex add(matrixComplex value, matrixComplex a);

double        times(double value, double a);
Complex       times(Complex value, double a);
vectorComplex times(vectorComplex value, double a);
matrixComplex times(matrixComplex value, double a);

matrixComplex product(matrixComplex A, matrixComplex B);
Complex tr(matrixComplex A);

template<class Fn, class N> void integrate_triangles_T(Fn fn, N& res, triangles tri, int band_index, chemical_potential mu) { // {{{
    init(res);

    int size = tri.faces.size();

    double factor = 1e0 / (2e0*pi)*(2e0*pi)*(2e0*pi);

    for(int i=0; i<size; i++) {
        kpoint center = {tri.faces[i].center[0], tri.faces[i].center[1], tri.faces[i].center[2]};
        double norm = 0e0;
        for(int axis=0; axis<space_dim; axis++) {
            norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
        }
        norm = std::sqrt(norm);
        double dS = tri.faces[i].dS / norm * factor;
        N c = times(fn(band_index, mu, center), dS);
        res = add(res, c);
    }
}; // }}}

template<class Fn, class N> void integrate_triangles_L(Fn fn, N& res, triangles tri, int valley, int band_index, chemical_potential mu) { // {{{
    init(res);

    int size = tri.faces.size();

    double factor = 1e0 / (2e0*pi)*(2e0*pi)*(2e0*pi);

    for(int i=0; i<size; i++) {
        kpoint center = {tri.faces[i].center[0], tri.faces[i].center[1], tri.faces[i].center[2]};
        double norm = 0e0;
        for(int axis=0; axis<space_dim; axis++) {
            norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
        }
        norm = std::sqrt(norm);
        double dS = tri.faces[i].dS / norm * factor;
        N c = times(fn(valley, band_index, mu, center), dS);
        res = add(res, c);
    }
}; // }}}

double get_DOS_T(triangles tri, int band_index, chemical_potential mu);
double get_DOS_L(triangles tri, int valley, int band_index, chemical_potential mu);

struct band {
    int index;
    chemical_potential mu_min, mu_max, dmu;
    int mu_mesh;
    std::vector<triangles> tri;
    std::vector<double> dos;
};

struct sys_T {
    std::vector<band> bands;
};

struct sys_L {
    std::vector<std::vector<band>> bands;
};

void get_band_T(band& b, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh);
sys_T get_T();
void sys_T_write(sys_T s);

void get_band_L(band& b, int valley, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh);
sys_L get_L();
void sys_L_write(sys_L s);

band set_band_T(int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh);
band set_band_L(int valley, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh);

using Green_function = matrixComplex;
Green_function get_green_function_T(Complex ene, kpoint k);
Green_function get_green_function_L(Complex ene, int valley, kpoint k);

using SHC = std::vector<matrixReal>;
SHC get_SHC_T(band b);
SHC get_SHC_L(band b, int valley);

using Conductivity = matrixComplex;
Conductivity get_conductivity_T(band b);
Conductivity get_conductivity_L(band b, int valley);

#endif // PARAMETERS_HPP
