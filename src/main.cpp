#include "parameters.hpp"

int main(){
    initialize();

//  output config.ini {{{
    {
    std::string file_name = "./config.ini";
    std::ofstream file(file_name);
        file << std::scientific << "[physics]" << std::endl;
        file << std::scientific << "a  = " << a << std::endl;
        file << std::scientific << "c  = " << c << std::endl;
        file << std::scientific << "g0 = " << g0 << std::endl;
        file << std::scientific << "b1 = [" << b1[0] << ", " << b1[1] << ", " << b1[2] << "]" << std::endl;
        file << std::scientific << "b2 = [" << b2[0] << ", " << b2[1] << ", " << b2[2] << "]" << std::endl;
        file << std::scientific << "b3 = [" << b3[0] << ", " << b3[1] << ", " << b3[2] << "]" << std::endl;
        file << std::endl;
        file << std::scientific << "kT = [" << kT[0] << ", " << kT[1] << ", " << kT[2] << "]" << std::endl;
        file << std::scientific << "kL1 = [" << kL[0][0] << ", " << kL[0][1] << ", " << kL[0][2] << "]" << std::endl;
        file << std::scientific << "kL2 = [" << kL[1][0] << ", " << kL[1][1] << ", " << kL[1][2] << "]" << std::endl;
        file << std::scientific << "kL3 = [" << kL[2][0] << ", " << kL[2][1] << ", " << kL[2][2] << "]" << std::endl;
        file << std::endl;
        file << std::scientific << "bands  = " <<  bands  << std::endl;
        file << std::scientific << "bandsT = " <<  bandsT << std::endl;
        file << std::scientific << "bandsL = " <<  bandsL << std::endl;
        file << std::scientific << "lowest_band_T = " << lowest_band_T << std::endl;
        file << std::scientific << "lowest_band_L = " << lowest_band_L << std::endl;
        file << std::scientific << "[numeric]" << std::endl;
    }
// }}}
// output kx-dependence @ T {{{
    int n = 50;
    std::string file_name = "./dat/EigenValue_T_"+std::to_string(bandsT)+"bands.dat";
    std::ofstream file(file_name);
    for (int i=-n; i<=n; i++) {
        double kx = cutoff*double(i)/double(n);
        double k[3] = {0e0,0e0,kx};
        matrixComplex H_T = set_T(k);

        vectorReal E_T = diagonalize_N(H_T);

        file << std::scientific << kx << ", ";
        for(int i=0; i<bandsT; i++) {
            file << std::scientific << E_T[i];
            if (i != bandsT-1) file << ", ";
        }
        file << std::endl;

    };
/// }}}
// output kx-dependence @ Ls {{{
    for(int valley=0; valley<valleys; valley++) {
        std::string file_name = "./dat/EigenValue_L"+std::to_string(valley+1)+"_"+std::to_string(bandsL)+"bands.dat";
        std::ofstream file(file_name);
        for (int i=-n; i<=n; i++) {
            double kx = cutoff*double(i)/double(n);
            double k[3] = {kx,0e0,0e0};
            matrixComplex H_L = set_L(valley, k);

            vectorReal E_L = diagonalize_N(H_L);

            file << std::scientific << kx << ", ";
            for(int i=0; i<bandsL; i++) {
                file << std::scientific << E_L[i];
                if (i != bandsL-1) file << ", ";
            }
            file << std::endl;

        };
    };
/// }}}

    int valley = 0;
    int band_index = 2;
    chemical_potential mu = 0e0;
    triangles tri = get_triangles_L(valley, band_index, mu);
    double dos = get_DOS_L(tri, valley, band_index, mu);
    std::string dos_file = "dos"+std::to_string(band_index)+".csv";
    std::ofstream ofs(dos_file, std::ios::app);
    ofs << std::scientific << tri.faces.size() << ", " << dos << std::endl;
    std::string filename = "tri_k"+std::to_string(k_mesh)+".csv";
    triangles_write(tri, filename);

//    chemical_potential mu_max = double(EL[valley][band_index])+5e-2;
//    chemical_potential mu_min = double(EL[valley][band_index])-1e-2;
//    band bt;
//    get_band_L(bt, valley, band_index, mu_min, mu_max, mu_mesh);

//    sys_T T = get_T();
//    sys_T_write(T);

//    sys_L L = get_L();
//    sys_L_write(L);

    return 0;
}
