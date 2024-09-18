#include "parameters.hpp"
#include <algorithm>

int thread_num = std::thread::hardware_concurrency();
//int thread_num = 1;

const double angstrom = 1e-10; // [m]
const double hbar     = 6.582119569e-16; // [eV s]
const double mass     = 9.10938356e-31; // [kg]
const double charge   = 1.60217662e-19; // [C]
const double muB      = 9.2740100783e-24; // [J/T]
const double muBeV    = 5.7883818060e-5; // [eV/T]

const double v0       = 1e6; // [m/s]

const double pi = 3.141592653589793e0; // pi
const double eps = 1e-12;
const Complex zi(0e0, 1e0);

const int H_dim = 16;
const int bands = H_dim;

//const int bandsT = 6; const int lowest_band_T = 7; // minimal
//const int bandsT = 8; const int lowest_band_T = 5; // not good
//const int bandsT = 8; const int lowest_band_T = 7; // better
//const int bandsT = 10; const int lowest_band_T = 7; // more better
const int bandsT = 12; const int lowest_band_T = 5; // almost full
//const int bandsT = 16; const int lowest_band_T = 1; // full

vectorReal band_edge_T;
vectorReal band_edge_T_sign0 = { 1e0, 1e0, 1e0, 1e0,-1e0,-1e0,-1e0,-1e0,-1e0,-1e0, 1e0, 1e0, -1e0,-1e0, 1e0, 1e0};
vectorReal band_edge_T_sign;

double mu_cutoff_T;
double mu_cutoff_L;
int mu_cutoff_mesh_T;
int mu_cutoff_mesh_L;

const int bandsL = 4; const int lowest_band_L = 9;
//const int bandsL = 12; const int lowest_band_L = 5;

matrixReal band_edge_L;
matrixReal band_edge_L_sign0 = { {-1e0,-1e0, 1e0, 1e0,-1e0,-1e0,-1e0,-1e0,-1e0,-1e0, 1e0, 1e0, -1e0,-1e0, 1e0, 1e0},
                                 {-1e0,-1e0, 1e0, 1e0,-1e0,-1e0,-1e0,-1e0,-1e0,-1e0, 1e0, 1e0, -1e0,-1e0, 1e0, 1e0},
                                 {-1e0,-1e0, 1e0, 1e0,-1e0,-1e0,-1e0,-1e0,-1e0,-1e0, 1e0, 1e0, -1e0,-1e0, 1e0, 1e0} };
matrixReal band_edge_L_sign;

const int space_dim = 3;
const int spin_dim = 3;
const double eps_phys = 1e-8;
const double eps_num = 1e-4;

const double a =  4.5332e0; // angstrom
const double c = 11.7967e0; // angstrom
const double g0 = 1.3861e0; // angstrom^-1

const double cutoff = 2e-1*g0;
double dk[3];
const int mu_mesh_T = 160;
const int mu_mesh_L = 160;
const int fermi_surface_mesh_lim_T = 3000;
const int fermi_surface_mesh_lim_L = 3000;
const double damping_T = 1.0e-4;
const double damping_L = 1.6e-4;

const vectorReal b1 = {-g0    ,-std::sqrt(3e0)*g0/3e0       , (a/c)*g0 };
const vectorReal b2 = { g0    ,-std::sqrt(3e0)*g0/3e0       , (a/c)*g0 };
const vectorReal b3 = { 0e0   , 2e0*std::sqrt(3e0)*g0/3e0   , (a/c)*g0 };

vectorReal kT;
vectorReal kL[3];

const std::string axises[] = {"x", "y", "z"};
const int valleys = 3;

const double sigma_imp = 10e0; // angstrom
//const vectorReal damping = {1e-5, 2e-5, 3e-5, 5e-5, 7e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4}; // eV
const vectorReal damping = {1e-5, 5e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4}; // eV
double nu_F;
double nu_F_T;
double nu_F_L[3];

std::vector<matrixComplex> vT;
std::vector<std::vector<matrixComplex>> vL;
vectorReal ET(bandsT);
std::vector<vectorReal> EL(valleys, vectorReal(bandsL));

std::vector<matrixComplex> mu_s_T;
std::vector<std::vector<matrixComplex>> mu_s_L;

std::vector<std::vector<matrixComplex>> v_s_T;
std::vector<std::vector<std::vector<matrixComplex>>> v_s_L;

std::vector<matrixComplex> sigma_T;
std::vector<std::vector<matrixComplex>> sigma_L;

std::vector<std::vector<matrixComplex>> v_sigma_T;
std::vector<std::vector<std::vector<matrixComplex>>> v_sigma_L;

matrixComplex impurityV1_T, impurityV2_T;
std::vector<matrixComplex> impurityV1_L, impurityV2_L;

std::mutex mtx;

void initialize() {
    using namespace std;

// set kT, kL[3] {{{
    kT.resize(3);
    kT = vectorReal(space_dim, 0e0);
    for (int valley=0; valley<valleys; valley++) {
        kL[valley].resize(3);
        kL[valley] = vectorReal(space_dim, 0e0);
    };
    for (int i=0; i<3; i++) {
        kT[i] = 5e-1* ( b1[i] + b2[i] + b3[i] );
        kL[0][i] = 5e-1 * b1[i];
        kL[1][i] = 5e-1 * b2[i];
        kL[2][i] = 5e-1 * b3[i];
    };
/// }}}
    // load data by Dr. Izaki {{{
    { // T point {{{
    int l = lowest_band_T - 1;
    int n = lowest_band_T + bandsT - 1;
//  velocity {{{
    vT.resize(space_dim);
    for (int axis=0; axis<space_dim; axis++) {
        matrixComplex v;
        string file_name = "./lib/izaki_Bi_SPH_data/v"+axises[axis]+"_T.dat";
//        cout << file_name << endl;
        ifstream file(file_name);
        string line;
        while ( getline(file,line) ) {
            istringstream stream(line);
            Complex c;
            vectorComplex row;
            while ( stream >> c ) {
                row.push_back(c);
            }
            v.push_back(row);
        }
        vT[axis].resize(bandsT);
        for(int i=l; i<n; i++) {
            vT[axis][i-l].resize(bandsT);
            for(int j=l; j<n; j++) {
                vT[axis][i-l][j-l] = v[i][j] / v0;
            }
        }

    }
// }}}
//  spin {{{
    sigma_T.resize(spin_dim);
    for (int axis=0; axis<spin_dim; axis++) {
        matrixComplex v;
        string file_name = "./lib/izaki_Bi_SPH_sigma_data/sig_"+axises[axis]+"_T.dat";
//        cout << file_name << endl;
        ifstream file(file_name);
        string line;
        while ( getline(file,line) ) {
            istringstream stream(line);
            Complex c;
            vectorComplex row;
            while ( stream >> c ) {
                row.push_back(c);
            }
            v.push_back(row);
        }
        sigma_T[axis].resize(bandsT);
        for(int i=l; i<n; i++) {
            sigma_T[axis][i-l].resize(bandsT);
            for(int j=l; j<n; j++) {
                sigma_T[axis][i-l][j-l] = v[i][j];
            }
        }
    }
// }}}
// EigenEnergy {{{
    vector<double> e(bands);
    string file_name = "./lib/izaki_Bi_SPH_data/Liu_allen_E_T.dat";
//    cout << file_name << endl;
    ifstream file(file_name);
    string line;
    getline(file,line);
    int i = 0;
    while ( getline(file,line) ) {
        istringstream stream(line);
        stream >> e[i];
        i++;
    }
    ET.assign(&e[l], &e[n]);
//    for (int i=0; i<ET.size(); i++) {
//        cout << ET[i] << endl;
//    }
// }}}
// magnetic moment {{{
    mu_s_T.resize(space_dim);
    for (int spin=0; spin<spin_dim; spin++) {
        matrixComplex v(bands, vectorComplex(bands, 0e0));
        string file_name = "./lib/izaki_Bi_SPH_data/mu_s_"+axises[spin]+"_T.dat";
//        cout << file_name << endl;
        ifstream file(file_name);
        string line;
        while ( getline(file,line) ) {
            istringstream stream(line);
            int i1, i2;
            Complex c;
            stream >> i1 >> i2 >> c;
//            std::cout << i1 << ", " << i2 << ", " << c << std::endl;
            v[i1-1][i2-1] = c;
        }
        mu_s_T[spin].resize(bandsT);
        for(int i=l; i<n; i++) {
            mu_s_T[spin][i-l].resize(bandsT);
            for(int j=l; j<n; j++) {
                mu_s_T[spin][i-l][j-l] = v[i][j] / muB * muBeV;
            }
        }
    }
// }}}
    { // impurity potential 1 {{{
    matrixComplex v(bands, vectorComplex(bands, 0e0));
    string file_name = "./lib/izaki_Bi_SPH_I8/check_I81_T.dat";
//    cout << file_name << endl;
    ifstream file(file_name);
    string line;
    while ( getline(file,line) ) {
        istringstream stream(line);
        int i1, i2;
        Complex c;
        stream >> i1 >> i2 >> c;
//        std::cout << i1 << ", " << i2 << ", " << c << std::endl;
        v[i1-1][i2-1] = c;
    }
    impurityV1_T.resize(bandsT);
    for(int i=l; i<n; i++) {
        impurityV1_T[i-l].resize(bandsT);
        for(int j=l; j<n; j++) {
            impurityV1_T[i-l][j-l] = v[i][j];
        }
    }
    } // }}}
    { // impurity potential 2 {{{
    matrixComplex v(bands, vectorComplex(bands, 0e0));
    string file_name = "./lib/izaki_Bi_SPH_I8/check_I82_T.dat";
//    cout << file_name << endl;
    ifstream file(file_name);
    string line;
    while ( getline(file,line) ) {
        istringstream stream(line);
        int i1, i2;
        Complex c;
        stream >> i1 >> i2 >> c;
//        std::cout << i1 << ", " << i2 << ", " << c << std::endl;
        v[i1-1][i2-1] = c;
    }
    impurityV2_T.resize(bandsT);
    for(int i=l; i<n; i++) {
        impurityV2_T[i-l].resize(bandsT);
        for(int j=l; j<n; j++) {
            impurityV2_T[i-l][j-l] = v[i][j];
        }
    }
    } // }}}
    } // }}}
    { // L points {{{
    int l = lowest_band_L - 1;
    int n = lowest_band_L + bandsL - 1;
    vL.resize(valleys);
// velocity {{{
    mu_s_L.resize(valleys);
    for (int valley=0; valley<valleys; valley++) {
        vL[valley].resize(space_dim);
        for (int axis=0; axis<space_dim; axis++) {
            matrixComplex v;
            string file_name = "./lib/izaki_Bi_SPH_data/v"+axises[axis]+"_L-"+to_string(valley+1)+".dat";
            ifstream file(file_name);
            string line;
            while ( getline(file,line) ) {
                istringstream stream(line);
                Complex c;
                vectorComplex row;
                while ( stream >> c ) {
                    row.push_back(c);
                }
                v.push_back(row);
            }
            vL[valley][axis].resize(bandsL);
            for(int i=l; i<n; i++) {
                vL[valley][axis][i-l].resize(bandsL);
                for(int j=l; j<n; j++) {
                    vL[valley][axis][i-l][j-l] = v[i][j] / v0;
                }
            }
        }
// }}}
// spin {{{
    sigma_L.resize(valleys);
    for (int valley=0; valley<valleys; valley++) {
        sigma_L[valley].resize(spin_dim);
        for (int axis=0; axis<spin_dim; axis++) {
            matrixComplex v;
            string file_name = "./lib/izaki_Bi_SPH_sigma_data/sig_"+axises[axis]+"_L-"+to_string(valley+1)+".dat";
            ifstream file(file_name);
            string line;
            while ( getline(file,line) ) {
                istringstream stream(line);
                complex<double> c;
                vector<complex<double>> row;
                while ( stream >> c ) {
                    row.push_back(c);
                }
                v.push_back(row);
            }
            sigma_L[valley][axis].resize(bandsL);
            for(int i=l; i<n; i++) {
                sigma_L[valley][axis][i-l].resize(bandsL);
                for(int j=l; j<n; j++) {
                    sigma_L[valley][axis][i-l][j-l] = v[i][j];
                }
            }
        }
    }
// }}}
// EigenEnergy {{{
        vector<double> e(bands);
        string file_name = "./lib/izaki_Bi_SPH_data/Liu_allen_E_L-"+to_string(valley+1)+".dat";
        ifstream file(file_name);
        string line;
        getline(file,line);
        int i = 0;
        while ( getline(file,line) ) {
            istringstream stream(line);
            stream >> e[i];
            i++;
        }
        EL[valley].assign(&e[lowest_band_L-1], &e[lowest_band_L+bandsL-1]);
// }}}
// magnetic moment {{{
        mu_s_L[valley].resize(space_dim);
        for (int axis=0; axis<space_dim; axis++) {
            matrixComplex v(bands, vectorComplex(bands, 0e0));
            string file_name = "./lib/izaki_Bi_SPH_data/mu_s_"+axises[axis]+"_L-"+to_string(valley+1)+".dat";
//            cout << file_name << endl;
            ifstream file(file_name);
            string line;
            while ( getline(file,line) ) {
                istringstream stream(line);
                int i1, i2;
                Complex c;
                stream >> i1 >> i2 >> c;
//                std::cout << i1 << ", " << i2 << ", " << c << std::endl;
                v[i1-1][i2-1] = c;
            }
            mu_s_L[valley][axis].resize(bandsL);
            for(int i=l; i<n; i++) {
                mu_s_L[valley][axis][i-l].resize(bandsL);
                for(int j=l; j<n; j++) {
                    mu_s_L[valley][axis][i-l][j-l] = v[i][j] / muB * muBeV;
                }
            }

        }
    }
// }}}
    impurityV1_L.resize(valleys);
    { // impurity potential 1 {{{
    for (int valley=0; valley<valleys; valley++) {
        matrixComplex v(bands, vectorComplex(bands, 0e0));
        string file_name = "./lib/izaki_Bi_SPH_I8/check_I81_L-"+to_string(valley+1)+".dat";
//        cout << file_name << endl;
        ifstream file(file_name);
        string line;
        while ( getline(file,line) ) {
            istringstream stream(line);
            int i1, i2;
            Complex c;
            stream >> i1 >> i2 >> c;
//            std::cout << i1 << ", " << i2 << ", " << c << std::endl;
            v[i1-1][i2-1] = c;
        }
        impurityV1_L[valley].resize(bandsL);
        for(int i=l; i<n; i++) {
            impurityV1_L[valley][i-l].resize(bandsL);
            for(int j=l; j<n; j++) {
                impurityV1_L[valley][i-l][j-l] = v[i][j];
            }
        }
    }
    } // }}}
    impurityV2_L.resize(valleys);
    { // impurity potential 2 {{{
    for (int valley=0; valley<valleys; valley++) {
        matrixComplex v(bands, vectorComplex(bands, 0e0));
        string file_name = "./lib/izaki_Bi_SPH_I8/check_I82_L-"+to_string(valley+1)+".dat";
//        cout << file_name << endl;
        ifstream file(file_name);
        string line;
        while ( getline(file,line) ) {
            istringstream stream(line);
            int i1, i2;
            Complex c;
            stream >> i1 >> i2 >> c;
//            std::cout << i1 << ", " << i2 << ", " << c << std::endl;
            v[i1-1][i2-1] = c;
        }
        impurityV2_L[valley].resize(bandsL);
        for(int i=l; i<n; i++) {
            impurityV2_L[valley][i-l].resize(bandsL);
            for(int j=l; j<n; j++) {
                impurityV2_L[valley][i-l][j-l] = v[i][j];
            }
        }
    }
    } // }}}
    } // }}}
    // }}}
    // v_s_T[axis][spin] {{{
    v_s_T.resize(space_dim);
    for(int axis=0; axis<space_dim; axis++) {
        v_s_T[axis].resize(bandsT);
        for(int spin=0; spin<spin_dim; spin++) {
            v_s_T[axis][spin].resize(bandsT);
            for(int i=0; i<bandsT; i++) {
                v_s_T[axis][spin][i].resize(bandsT);
            }

            for(int i=0; i<bandsT; i++) {
                for(int j=0; j<bandsT; j++) {
                    Complex c = 0e0;
                    for(int k=0; k<bandsT; k++) {
                        c += (vT[axis][i][k]*mu_s_T[spin][k][j] + mu_s_T[spin][i][k]*vT[axis][k][j])*5e-1;
                    }
                    v_s_T[axis][spin][i][j] = c;
                }
            }
        }
    }
    // }}}
    // v_s_L[valley][axis][spin] {{{
    v_s_L.resize(valleys);
    for(int valley=0; valley<valleys; valley++) {
        v_s_L[valley].resize(space_dim);
        for(int axis=0; axis<space_dim; axis++) {
            v_s_L[valley][axis].resize(bandsL);
            for(int spin=0; spin<spin_dim; spin++) {
                v_s_L[valley][axis][spin].resize(bandsL);
                for(int i=0; i<bandsL; i++) {
                    v_s_L[valley][axis][spin][i].resize(bandsL);
                }

                for(int i=0; i<bandsL; i++) {
                    for(int j=0; j<bandsL; j++) {
                        Complex c = 0e0;
                        for(int k=0; k<bandsL; k++) {
                            c += (vL[valley][axis][i][k]*mu_s_L[valley][spin][k][j] + mu_s_L[valley][spin][i][k]*vL[valley][axis][k][j])*5e-1;
                        }
                        v_s_L[valley][axis][spin][i][j] = c;
                    }
                }
            }
        }
    }
    // }}}
    // v_sigma_T[axis][spin] {{{
    v_sigma_T.resize(space_dim);
    for(int axis=0; axis<space_dim; axis++) {
        v_sigma_T[axis].resize(bandsT);
        for(int spin=0; spin<spin_dim; spin++) {
            v_sigma_T[axis][spin].resize(bandsT);
            for(int i=0; i<bandsT; i++) {
                v_sigma_T[axis][spin][i].resize(bandsT);
            }

            for(int i=0; i<bandsT; i++) {
                for(int j=0; j<bandsT; j++) {
                    Complex c = 0e0;
                    for(int k=0; k<bandsT; k++) {
                        c += (vT[axis][i][k]*sigma_T[spin][k][j] + sigma_T[spin][i][k]*vT[axis][k][j])*5e-1;
                    }
                    v_sigma_T[axis][spin][i][j] = c;
                }
            }
        }
    }
    // }}}
    // v_sigma_L[valley][axis][spin] {{{
    v_sigma_L.resize(valleys);
    for(int valley=0; valley<valleys; valley++) {
        v_sigma_L[valley].resize(space_dim);
        for(int axis=0; axis<space_dim; axis++) {
            v_sigma_L[valley][axis].resize(bandsL);
            for(int spin=0; spin<spin_dim; spin++) {
                v_sigma_L[valley][axis][spin].resize(bandsL);
                for(int i=0; i<bandsL; i++) {
                    v_sigma_L[valley][axis][spin][i].resize(bandsL);
                }

                for(int i=0; i<bandsL; i++) {
                    for(int j=0; j<bandsL; j++) {
                        Complex c = 0e0;
                        for(int k=0; k<bandsL; k++) {
                            c += (vL[valley][axis][i][k]*sigma_L[valley][spin][k][j] + sigma_L[valley][spin][i][k]*vL[valley][axis][k][j])*5e-1;
                        }
                        v_sigma_L[valley][axis][spin][i][j] = c;
                    }
                }
            }
        }
    }
    // }}}
    // set band_edge T {{{
    band_edge_T.resize(bandsT);
    band_edge_T_sign.resize(bandsT);
    int l = lowest_band_T - 1;
    int n = lowest_band_T + bandsT - 1;
    for(int i=l; i<n; i++) {
        band_edge_T[i-l] = ET[i-l];
        band_edge_T_sign[i-l] = band_edge_T_sign0[i];
    }
    // }}}
    // set band_edge L {{{
    band_edge_L.resize(valleys);
    band_edge_L_sign.resize(valleys);
    l = lowest_band_L - 1;
    n = lowest_band_L + bandsL - 1;
    for(int valley=0; valley<valleys; valley++) {
        band_edge_L[valley].resize(bandsL);
        band_edge_L_sign[valley].resize(bandsL);
        for(int i=l; i<n; i++) {
            band_edge_L[valley][i-l] = EL[valley][i-l];
            band_edge_L_sign[valley][i-l] = band_edge_L_sign0[valley][i];
        }
    }
    // }}}
//load data DOS @ Fermi level {{{
    string file_name = "./lib/data/dos_mu0e0.dat";
//    cout << file_name << endl;
    ifstream file(file_name);
    string line;
    getline(file,line);
    istringstream stream(line);
    double r;
    vectorReal dos;
    while ( stream >> r ) {
        dos.push_back(r);
    }
    nu_F_T    = dos[0];
    nu_F_L[0] = dos[1];
    nu_F_L[1] = dos[2];
    nu_F_L[2] = dos[3];
//    std::cout << "DOS: " << nu_F_T << ", " << nu_F_L[0] << ", " << nu_F_L[1] << ", " << nu_F_L[2] << std::endl;

    nu_F = 0e0;
    for (auto nu : dos) nu_F += nu;
// }}}
};

void get_DOS_at_Fermi_level() { // {{{
    int band_index = 4;
    chemical_potential mu = 0e0;
    triangles tri = get_triangles_T(band_index, mu);
    nu_F_T = get_DOS_T(tri, band_index, mu);
    nu_F = nu_F_T;
    std::cout << std::scientific << std::setprecision(15) << "DOS of T = " << nu_F_T << std::endl;

    for(int valley=0; valley<valleys; valley++) {
        int band_index = 6; // for 12bands
        triangles tri = get_triangles_L(valley, band_index, mu);
        nu_F_L[valley] = get_DOS_L(tri, valley, band_index, mu);
        std::cout <<  std::scientific << std::setprecision(15) << "DOS of L" << std::to_string(valley) << " = " << nu_F_L[valley] << std::endl;

        nu_F = nu_F + nu_F_L[valley];
    }

    std::string filename = "lib/data/dos_mu0e0.dat";
    std::ofstream of(filename);
    of << std::scientific << std::setprecision(15) << nu_F_T;
    for(int valley=0; valley<valleys; valley++) {
        of <<  std::scientific << std::setprecision(15) << "  " << nu_F_L[valley];
    }
    of << std::endl;
}; // }}}

void show_impurity_potentials() { // {{{
    for(int valley=0; valley<valleys; valley++) {
        for (auto v: impurityV1_L[valley]) {
            for (auto c: v) {
                std::cout << std::right << std::fixed << std::setprecision(6)
                          << std::setw(10) << c.real()
                          << ", ";
            }
            std::cout << std::left << std::setw(5) << "" << std::left << std::setw(5) << "|";
            for (auto c: v) {
                std::cout << std::right << std::fixed << std::setprecision(6)
                          << std::setw(10) << c.imag()
                          << ", ";
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;

        for (auto v: impurityV2_L[valley]) {
            for (auto c: v) {
                std::cout << std::right << std::fixed << std::setprecision(6)
                          << std::setw(10) << c.real()
                          << ", ";
            }
            std::cout << std::left << std::setw(5) << "" << std::left << std::setw(5) << "|";
            for (auto c: v) {
                std::cout << std::right << std::fixed << std::setprecision(6)
                          << std::setw(10) << c.imag()
                          << ", ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}; // }}}

void set_isotropic() { // {{{
    const matrixComplex sigma_x = { {0e0, 1e0}, {1e0, 0e0} };
    const matrixComplex sigma_y = { {0e0,- zi}, { zi, 0e0} };
    const matrixComplex sigma_z = { {1e0, 0e0}, {0e0,-1e0} };
    tensor2Complex v(3, matrixComplex(4, vectorComplex(4, 0e0)));

    tensor2Complex mu_s(3, matrixComplex(4, vectorComplex(4, 0e0)));
    tensor2Complex sigma(3, matrixComplex(4, vectorComplex(4, 0e0)));

    for(int k=0; k<2; k++) {
        for(int l=0; l<2; l++) {
            v[0][k][2+l] =  zi*sigma_x[k][l];
            v[0][2+k][l] = -zi*sigma_x[k][l];

            v[1][k][2+l] =  zi*sigma_y[k][l];
            v[1][2+k][l] = -zi*sigma_y[k][l];

            v[2][k][2+l] =  zi*sigma_z[k][l];
            v[2][2+k][l] = -zi*sigma_z[k][l];

            mu_s[0][k][l]     =   sigma_x[k][l];
            mu_s[0][2+k][2+l] = - sigma_x[k][l];

            mu_s[1][k][l]     =   sigma_y[k][l];
            mu_s[1][2+k][2+l] = - sigma_y[k][l];

            mu_s[2][k][l]     =   sigma_z[k][l];
            mu_s[2][2+k][2+l] = - sigma_z[k][l];
      }
  }
    for (int valley=0; valley<valleys; valley++) {
        for (int axis=0; axis<space_dim; axis++) {
            for(int i=0; i<4; i++) {
                for(int j=0; j<4; j++) {
                    vL[valley][axis][i][j] = v[axis][i][j];
                }
            }
        }
        for (int spin=0; spin<spin_dim; spin++) {
            for(int i=0; i<4; i++) {
                for(int j=0; j<4; j++) {
                    mu_s_L[valley][spin][i][j] = mu_s[spin][i][j];
                }
            }
        }
        for (int axis=0; axis<space_dim; axis++) {
            for (int spin=0; spin<spin_dim; spin++) {
                for(int i=0; i<4; i++) {
                    for(int j=0; j<4; j++) {
                        Complex c = 0e0;
                        for(int k=0; k<4; k++) {
                            c += (vL[valley][axis][i][k]*mu_s_L[valley][spin][k][j] + mu_s_L[valley][spin][i][k]*vL[valley][axis][k][j])*5e-1;
                        }
                        v_s_L[valley][axis][spin][i][j] = c;

                        c = 0e0;
                        for(int k=0; k<4; k++) {
                            c += (vL[valley][axis][i][k]*sigma_L[valley][spin][k][j] + sigma_L[valley][spin][i][k]*vL[valley][axis][k][j])*5e-1;
                        }
                        v_sigma_L[valley][axis][spin][i][j] = c;
                    }
                }
            }
        }
    }
} // }}}

void init(double& value, double& res) { // {{{
    value = 0e0;
}; // }}}

void init(Complex& value, Complex& res) { // {{{
    value = 0e0;
}; // }}}

void init(vectorComplex& value, vectorComplex& res) { // {{{
    value.resize(res.size());
    for(int i=0; i<res.size(); i++) {
        value[i] = 0e0;
    }
}; // }}}

void init(matrixComplex& value, matrixComplex& res) { // {{{
    value.resize(res.size());
    for(int i=0; i<res.size(); i++) {
        value[i].resize(res.at(0).size());
        for(int j=0; j<res.at(0).size(); j++) {
            value[i][j] = 0e0;
        }
    }
}; // }}}

void init(tensor2Complex& value, tensor2Complex& res) { // {{{
    value.resize(res.size());
    for(int i=0; i<res.size(); i++) {
        value[i].resize(res.at(0).size());
        for(int j=0; j<res.at(0).size(); j++) {
            value[i][j].resize(res.at(0).at(0).size());
            for(int k=0; k<res.at(0).at(0).size(); k++) {
                value[i][j][k] = 0e0;
            }
        }
    }
}; // }}}

double add(double value, double a) { // {{{
    return value+a;
}; // }}}

Complex add(Complex value, Complex a) { // {{{
    return value+a;
}; // }}}

vectorComplex add(vectorComplex value, vectorComplex a) { // {{{
    vectorComplex res(value.size(), 0e0);
    for(int i=0; i<value.size(); i++) {
        res[i] = value[i]+a[i];
    }
    return res;
}; // }}}

matrixComplex add(matrixComplex value, matrixComplex a) { // {{{
    matrixComplex res(value.size(), vectorComplex(value.at(0).size(), 0e0));
    for(int i=0; i<value.size(); i++) {
        for(int j=0; j<value.size(); j++) {
            res[i][j] = value[i][j] + a[i][j];
        }
    }
    return res;
}; // }}}

tensor2Complex add(tensor2Complex value, tensor2Complex a) { // {{{
    tensor2Complex res(value.size(), matrixComplex(value.at(0).size(), vectorComplex(value.at(0).at(0).size(), 0e0)));
    for(int i=0; i<value.size(); i++) {
        for(int j=0; j<value.at(0).size(); j++) {
            for(int k=0; k<value.at(0).at(0).size(); k++) {
                res[i][j][k] = value[i][j][k] + a[i][j][k];
            }
        }
    }
    return res;
}; // }}}

vectorComplex minus(vectorComplex value, vectorComplex a) { // {{{
    vectorComplex res(value.size(), 0e0);
    for(int i=0; i<value.size(); i++) {
        res[i] = value[i]-a[i];
    }
    return res;
}; // }}}

matrixComplex minus(matrixComplex value, matrixComplex a) { // {{{
    matrixComplex res(value.size(), vectorComplex(value.at(0).size(), 0e0));
    for(int i=0; i<value.size(); i++) {
        for(int j=0; j<value.size(); j++) {
            res[i][j] = value[i][j] - a[i][j];
        }
    }
    return res;
}; // }}}

tensor2Complex minus(tensor2Complex value, tensor2Complex a) { // {{{
    tensor2Complex res(value.size(), matrixComplex(value.at(0).size(), vectorComplex(value.at(0).at(0).size(), 0e0)));
    for(int i=0; i<value.size(); i++) {
        for(int j=0; j<value.at(0).size(); j++) {
            for(int k=0; k<value.at(0).at(0).size(); k++) {
                res[i][j][k] = value[i][j][k] - a[i][j][k];
            }
        }
    }
    return res;
}; // }}}

double times(double value, double a) { // {{{
    return value*a;
}; // }}}

Complex times(Complex value, double a) { // {{{
    return value*a;
}; // }}}

vectorComplex times(vectorComplex value, double a) { // {{{
    vectorComplex res(value.size(), 0e0);
    for(int i=0; i<value.size(); i++) {
        res[i] = value[i]*a;
    }
    return res;
}; // }}}

matrixComplex times(matrixComplex value, double a) { // {{{
    matrixComplex res(value.size(), vectorComplex(value.at(0).size(), 0e0));
    for(int i=0; i<value.size(); i++) {
        for(int j=0; j<value.at(0).size(); j++) {
            res[i][j] = value[i][j] * a;
        }
    }
    return res;
}; // }}}

tensor2Complex times(tensor2Complex value, double a) { // {{{
    tensor2Complex res(value.size(), matrixComplex(value.at(0).size(), vectorComplex(value.at(0).at(0).size(), 0e0)));
    for(int i=0; i<value.size(); i++) {
        for(int j=0; j<value.at(0).size(); j++) {
            for(int k=0; k<value.at(0).at(0).size(); k++) {
                res[i][j][k] = value[i][j][k] * a;
            }
        }
    }
    return res;
}; // }}}

matrixComplex product(matrixComplex A, matrixComplex B) { // {{{
    const int N = A.size();
    matrixComplex C(N, vectorComplex(N, 0e0));
    Complex da;
    for(int i=0; i<N; i++) {
        for(int k=0; k<N; k++) {
            da = A[i][k];
            for(int j=0; j<N; j++) {
                C[i][j] = C[i][j] + da*B[k][j];
            }
        }
    }
    return C;
}; // }}}

Complex tr(matrixComplex A) { // {{{
    Complex res = 0e0;
    for(int i=0; i<A.size(); i++) {
        res += A[i][i];
    }
    return res;
//    const int n = A.size();
//
//    vectorReal diag;
//    diag.resize(n+1);
//    for(int i=0; i<n; i++) {
//        diag[i] = A[i][i].real();
//    }
//    diag[n] = 0e0;
//    std::sort(diag.begin(), diag.end());
//    auto zero = std::find(diag.begin(), diag.end(), 0e0);
//    const int index = std::distance(diag.begin(), zero);
//    const int min = std::min(index, n-index);
//
//    if (min != index) {
//        auto start = std::find(diag.begin(), zero, diag[min]);
//        std::sort(start, zero, std::greater<int>{});
//    }
//    for(int i=min; i<n-min+1; i++) {
//        res += diag[i];
//    }
//    for(int i=0; i<min; i++) {
//        res += (diag[i] + diag[n-i]);
//    }
//
//    return res;
}; // }}}

