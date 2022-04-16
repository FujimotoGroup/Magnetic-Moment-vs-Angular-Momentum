#include "parameters.hpp"
#include <algorithm>

int thread_num = std::thread::hardware_concurrency();

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
const int bandsT = 8; const int lowest_band_T = 7; // better
//const int bandsT = 10; const int lowest_band_T = 7; // more better
//const int bandsT = 12; const int lowest_band_T = 5; // almost full
//const int bandsT = 16; const int lowest_band_T = 1; // full

const int bandsL = 4; const int lowest_band_L = 9;
//const int bandsL = 16; const int lowest_band_L = 1;

const int space_dim = 3;
const int spin_dim = 3;
const double eps_phys = 1e-8;
const double eps_num = 1e-4;

const double a =  4.5332e0; // angstrom
const double c = 11.7967e0; // angstrom
const double g0 = 1.3861e0; // angstrom^-1

const double cutoff = 2e-1*g0;
double dk[3];
const int k_mesh = 100;
const int k_mesh_more = 70;
const int mu_mesh = 60;
//const int k_mesh = 30; const int k_mesh_more = 50; const int mu_mesh = 32;

const vectorReal b1 = {-g0    ,-std::sqrt(3e0)*g0/3e0       , (a/c)*g0 };
const vectorReal b2 = { g0    ,-std::sqrt(3e0)*g0/3e0       , (a/c)*g0 };
const vectorReal b3 = { 0e0   , 2e0*std::sqrt(3e0)*g0/3e0   , (a/c)*g0 };

vectorReal kT;
vectorReal kL[3];

const std::string axises[] = {"x", "y", "z"};
const int valleys = 3;

std::vector<matrixComplex> vT;
std::vector<std::vector<matrixComplex>> vL;
vectorReal ET(bandsT);
std::vector<vectorReal> EL(valleys, vectorReal(bandsL));

std::vector<matrixComplex> mu_s_T;
std::vector<std::vector<matrixComplex>> mu_s_L;

std::vector<std::vector<matrixComplex>> v_s_T;
std::vector<std::vector<std::vector<matrixComplex>>> v_s_L;

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
    // T point {{{
    int l = lowest_band_T - 1;
    int n = lowest_band_T + bandsT - 1;
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
    // L points {{{
    l = lowest_band_L - 1;
    n = lowest_band_L + bandsL - 1;
    vL.resize(valleys);
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
                complex<double> c;
                vector<complex<double>> row;
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
                v_s_L[valley][axis][spin].resize(bandsT);
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
};

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
    const int n = A.size();

    vectorReal diag;
    diag.resize(n+1);
    for(int i=0; i<n; i++) {
        diag[i] = A[i][i].real();
    }
    diag[n] = 0e0;
    std::sort(diag.begin(), diag.end());
    auto zero = std::find(diag.begin(), diag.end(), 0e0);
    const int index = std::distance(diag.begin(), zero);
    const int min = std::min(index, n-index);

    if (min != index) {
        auto start = std::find(diag.begin(), zero, diag[min]);
        std::sort(start, zero, std::greater<int>{});
    }
    for(int i=min+1; i<n-min+1; i++) {
        res += diag[i];
    }
    for(int i=0; i<min; i++) {
        res += diag[i] + diag[n-i];
    }

    return res;
}; // }}}

