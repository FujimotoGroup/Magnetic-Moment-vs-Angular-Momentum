#include "parameters.hpp"

const double angstrom = 1e-10;
const double hbar     = 6.582119569e-16;
const double mass     = 9.10938356e-31;
const double charge   = 1.60217662e-19;

const int H_dim = 16;
const int bands = H_dim;
const int bandsT = 8; const int lowest_band_T = 5;
const int bandsL = 4; const int lowest_band_L = 9;

const int space_dim = 3;
const int spin_dim = 3;
const double eps_phys = 1e-8;

const double a =  4.5332e0; // angstrom
const double c = 11.7967e0; // angstrom
const double g0 = 1.3861e0; // angstrom^-1

const vectorReal b1 = {-g0    ,-std::sqrt(3e0)*g0/3e0       , (a/c)*g0 };
const vectorReal b2 = { g0    ,-std::sqrt(3e0)*g0/3e0       , (a/c)*g0 };
const vectorReal b3 = { 0e0   , 2e0*std::sqrt(3e0)*g0/3e0   , (a/c)*g0 };

const std::string axises[] = {"x", "y", "z"};
const int valleys = 3;

std::vector<matrixComplex> vT(space_dim);
std::vector<std::vector<matrixComplex>> vL;
vectorReal ET(bandsT);
std::vector<vectorReal> EL(valleys, vectorReal(bandsL));

void initialize() {
    using namespace std;

    // load data by Dr. Izaki {{{
    // T point {{{
    int l = lowest_band_T - 1;
    int n = lowest_band_T + bandsT - 1;
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
                vT[axis][i-l][j-l] = v[i][j];
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
    // }}}
    // L points {{{
    for (int valley=0; valley<valleys; valley++) {
        for (int axis=0; axis<space_dim; axis++) {
            string file_name = "./lib/izaki_Bi_SPH_data/v"+axises[axis]+"_L"+to_string(valley+1)+".dat";
//            cout << file_name << endl;
            ifstream file(file_name);
            string line;
            while ( getline(file,line) ) {
                istringstream stream(line);
                complex<double> c;
                vector<complex<double>> row;
                while ( stream >> c ) {
                    row.push_back(c);
                }
                vL[valley][axis].push_back(row);
            }
        }

        vector<double> e(bands);
        string file_name = "./lib/izaki_Bi_SPH_data/Liu_allen_E_L-"+to_string(valley+1)+".dat";
//        cout << file_name << endl;
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
//        for (int i=0; i<EL[valley].size(); i++) {
//            cout << EL[valley][i] << endl;
//        }
        }
    // }}}
    // }}}
}
