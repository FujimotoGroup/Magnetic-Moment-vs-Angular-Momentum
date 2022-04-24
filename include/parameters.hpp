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
#include <algorithm>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

extern std::mutex mtx;
extern int thread_num;

using Complex = std::complex<double>;
using vectorReal    = std::vector<double>;
using vectorComplex = std::vector<std::complex<double>>;
using matrixComplex = std::vector<vectorComplex>;
using tensor2Complex = std::vector<matrixComplex>;
using tensor3Complex = std::vector<tensor2Complex>;
using matrixReal = std::vector<vectorReal>;
using chemical_potential = double;
using Energy = double;

extern const double angstrom;
extern const double hbar;
extern const double mass;
extern const double charge;
extern const double muB;

extern const double v0;

extern const double eps;
extern const double pi;
extern const Complex zi;

extern const int bands;
extern const int bandsT;
extern const int bandsL;
extern const int lowest_band_T;
extern const int lowest_band_L;

extern double mu_cutoff_T;
extern double mu_cutoff_L;
extern int mu_cutoff_mesh_T;
extern int mu_cutoff_mesh_L;

extern vectorReal band_edge_T;
extern vectorReal band_edge_T_sign;
extern matrixReal band_edge_L;
extern matrixReal band_edge_L_sign;

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
void set_isotropic();
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

using Green_function = matrixComplex;
using SHC = tensor2Complex;
using Conductivity = matrixComplex;

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

void init(double& value, double& res);
void init(Complex& value, Complex& res);
void init(vectorComplex& value, vectorComplex& res);
void init(matrixComplex& value, matrixComplex& res);
void init(tensor2Complex& value, tensor2Complex& res);

double         add(double value, double a);
Complex        add(Complex value, Complex a);
vectorComplex  add(vectorComplex value, vectorComplex a);
matrixComplex  add(matrixComplex value, matrixComplex a);
tensor2Complex add(tensor2Complex value, tensor2Complex a);

vectorComplex  minus(vectorComplex value, vectorComplex a);
matrixComplex  minus(matrixComplex value, matrixComplex a);
tensor2Complex minus(tensor2Complex value, tensor2Complex a);

double         times(double value, double a);
Complex        times(Complex value, double a);
vectorComplex  times(vectorComplex value, double a);
matrixComplex  times(matrixComplex value, double a);
tensor2Complex times(tensor2Complex value, double a);

matrixComplex product(matrixComplex A, matrixComplex B);
Complex tr(matrixComplex A);

template<class Fn, class N> void integrate_triangles_T(Fn fn, N& res, triangles tri, int band_index, chemical_potential mu) { // {{{
    init(res, res);

    std::vector<std::thread> threads;
    threads.resize(thread_num);
    std::vector<N> part;
    part.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        init(part[i_thread], res);
        auto func = [](int i_thread, triangles tri, Fn& fn, N& part, int band_index, chemical_potential mu) {
            double coef = 1e0 / (2e0*pi)*(2e0*pi)*(2e0*pi);
            int size = tri.faces.size();
            for(int i=i_thread; i<size; i=i+thread_num) {
                kpoint center = {tri.faces[i].center[0], tri.faces[i].center[1], tri.faces[i].center[2]};
                double norm = 0e0;
                for(int axis=0; axis<space_dim; axis++) {
                    norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
                }
                norm = std::sqrt(norm);
                double dS = tri.faces[i].dS / norm * coef;
                N c = times(fn(band_index, mu, center), dS);
                part = add(part, c);
            }
        };
        threads[i_thread] = std::thread(func, i_thread, tri, std::ref(fn), std::ref(part[i_thread]), band_index, mu);
    }

    for(auto& thread : threads){
        thread.join();
    }

    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        res = add(res, part[i_thread]);
    }
}; // }}}

template<class Fn, class N> void integrate_triangles_L(Fn fn, N& res, triangles tri, int valley, int band_index, chemical_potential mu) { // {{{
    init(res, res);

    std::vector<std::thread> threads;
    threads.resize(thread_num);
    std::vector<N> part;
    part.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        init(part[i_thread], res);
        auto func = [](int i_thread, triangles tri, Fn& fn, N& part, int valley, int band_index, chemical_potential mu) {
            int size = tri.faces.size();
            for(int i=i_thread; i<size; i=i+thread_num) {
                kpoint center = {tri.faces[i].center[0], tri.faces[i].center[1], tri.faces[i].center[2]};
                double norm = 0e0;
                for(int axis=0; axis<space_dim; axis++) {
                    norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
                }
                norm = std::sqrt(norm);
                double dS = tri.faces[i].dS / norm;
                N c = times(fn(valley, band_index, mu, center), dS);
                part = add(part, c);
            }
        };
        threads[i_thread] = std::thread(func, i_thread, tri, std::ref(fn), std::ref(part[i_thread]), valley, band_index, mu);
    }

    for(auto& thread : threads){
        thread.join();
    }

    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        res = add(res, part[i_thread]);
    }
    double coef = 1e0 / ((2e0*pi)*(2e0*pi)*(2e0*pi));
    res = times(res, coef);
}; // }}}

double get_E_T(int band_index, kpoint k);
double get_E_L(int valley, int band_index, kpoint k);

double get_DOS_T(triangles tri, int band_index, chemical_potential mu);
double get_DOS_L(triangles tri, int valley, int band_index, chemical_potential mu);

struct band {
    int index;
    std::vector<chemical_potential> ene;
    int mesh;
    std::vector<triangles> tri;
    std::vector<double> dos;
};

band set_band_T(int band_index, Energy ene_min, Energy ene_max, int ene_mesh);
band set_band_2n_T(int band_index, Energy ene_center, Energy delta, int n);
band set_band_L(int valley, int band_index, Energy ene_min, Energy ene_max, int ene_mesh);
band set_band_2n_L(int valley, int band_index, Energy ene_center, Energy delta, int n, double power);

void init_band_L(band& b, int valley, int band_index);
void add_fs_L(band& b, int valley, chemical_potential mu);

band combine_band(band b1, band b2);
band combine_band_2n(band b_global, band b_local);

void write_res(Conductivity sigma, chemical_potential mu,  std::string filename);
void write_res(SHC sigma, chemical_potential mu,  std::string filename);

template<class Fn, class N> void integrate_band_L(Fn fn, N& res, band b, int valley, chemical_potential mu) { // {{{
//    std::string filename = "sigma.csv";
//    std::string filename1 = "sigma1.csv";
    N sigma;
    Energy dmu;
    int i_mu = 0;
        init(sigma, res);
        integrate_triangles_L(fn, sigma, b.tri[i_mu], valley, b.index, mu);
        dmu = (b.ene[i_mu+1] - b.ene[i_mu])*5e-1;
//        write_res(sigma, b.ene[i_mu]-mu, filename1);
        sigma = times(sigma, dmu);
//        write_res(sigma, b.ene[i_mu]-mu, filename);
        res = add(res, sigma);
    for(i_mu=1; i_mu<b.mesh-1; i_mu++) {
        init(sigma, res);
        integrate_triangles_L(fn, sigma, b.tri[i_mu], valley, b.index, mu);
        dmu = (b.ene[i_mu+1] - b.ene[i_mu-1])*5e-1;
//        write_res(sigma, b.ene[i_mu]-mu, filename1);
        sigma = times(sigma, dmu);
//        write_res(sigma, b.ene[i_mu]-mu, filename);
        res = add(res, sigma);
    }
    i_mu = b.mesh-1;
        init(sigma, res);
        integrate_triangles_L(fn, sigma, b.tri[i_mu], valley, b.index, mu);
        dmu = (b.ene[i_mu] - b.ene[i_mu-1])*5e-1;
//        write_res(sigma, b.ene[i_mu]-mu, filename1);
        sigma = times(sigma, dmu);
//        write_res(sigma, b.ene[i_mu]-mu, filename);
        res = add(res, sigma);
}; // }}}

Green_function get_green_function_T(Complex ene, kpoint k);
Green_function get_green_function_L(Complex ene, int valley, kpoint k);

SHC get_SHC_T1(band b);
SHC get_SHC_T2(band b);
SHC get_SHC_L1(band b, Energy epsilon, chemical_potential mu, int valley);
SHC get_SHC_L2(band b, Energy epsilon, chemical_potential mu, int valley);

Conductivity get_conductivity_T(band b);
Conductivity get_conductivity_L(band b, Energy epsilon, chemical_potential mu, int valley);

void set_response_L(chemical_potential ene_min, chemical_potential ene_max, int ene_mesh, int valley, int band_index);
#endif // PARAMETERS_HPP
