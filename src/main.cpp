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
    std::string file_name_x = "./dat/EigenValue_T_x_"+std::to_string(bandsT)+"bands.dat"; // x {{{
    std::ofstream file_x(file_name_x);
    for (int i=-n; i<=n; i++) {
        double kx = cutoff*double(i)/double(n);
        double k[3] = {kx,0e0,0e0};
        matrixComplex H_T = set_T(k);

        vectorReal E_T = diagonalize_N(H_T);

        file_x << std::scientific << kx << ", ";
        for(int i=0; i<bandsT; i++) {
            file_x << std::scientific << E_T[i];
            if (i != bandsT-1) file_x << ", ";
        }
        file_x << std::endl;

    }; // }}}
    std::string file_name_y = "./dat/EigenValue_T_y_"+std::to_string(bandsT)+"bands.dat"; // y {{{
    std::ofstream file_y(file_name_y);
    for (int i=-n; i<=n; i++) {
        double kx = cutoff*double(i)/double(n);
        double k[3] = {0e0,kx,0e0};
        matrixComplex H_T = set_T(k);

        vectorReal E_T = diagonalize_N(H_T);

        file_y << std::scientific << kx << ", ";
        for(int i=0; i<bandsT; i++) {
            file_y << std::scientific << E_T[i];
            if (i != bandsT-1) file_y << ", ";
        }
        file_y << std::endl;

    }; // }}}
    std::string file_name_z = "./dat/EigenValue_T_z_"+std::to_string(bandsT)+"bands.dat"; // z {{{
    std::ofstream file_z(file_name_z);
    for (int i=-n; i<=n; i++) {
        double kx = cutoff*double(i)/double(n);
        double k[3] = {0e0,0e0,kx};
        matrixComplex H_T = set_T(k);

        vectorReal E_T = diagonalize_N(H_T);

        file_z << std::scientific << kx << ", ";
        for(int i=0; i<bandsT; i++) {
            file_z << std::scientific << E_T[i];
            if (i != bandsT-1) file_z << ", ";
        }
        file_z << std::endl;

    }; // }}}
/// }}}
// output kx-dependence @ Ls {{{
    for(int valley=0; valley<valleys; valley++) {
        std::string file_name_x = "./dat/EigenValue_L"+std::to_string(valley+1)+"_x_"+std::to_string(bandsL)+"bands.dat"; // x {{{
        std::ofstream file_x(file_name_x);
        for (int i=-n; i<=n; i++) {
            double kx = cutoff*double(i)/double(n);
            double k[3] = {kx,0e0,0e0};
            matrixComplex H_L = set_L(valley, k);

            vectorReal E_L = diagonalize_N(H_L);

            file_x << std::scientific << kx << ", ";
            for(int i=0; i<bandsL; i++) {
                file_x << std::scientific << E_L[i];
                if (i != bandsL-1) file_x << ", ";
            }
            file_x << std::endl;

        }; // }}}
        std::string file_name_y = "./dat/EigenValue_L"+std::to_string(valley+1)+"_y_"+std::to_string(bandsL)+"bands.dat"; // y {{{
        std::ofstream file_y(file_name_y);
        for (int i=-n; i<=n; i++) {
            double kx = cutoff*double(i)/double(n);
            double k[3] = {0e0,kx,0e0};
            matrixComplex H_L = set_L(valley, k);

            vectorReal E_L = diagonalize_N(H_L);

            file_y << std::scientific << kx << ", ";
            for(int i=0; i<bandsL; i++) {
                file_y << std::scientific << E_L[i];
                if (i != bandsL-1) file_y << ", ";
            }
            file_y << std::endl;

        }; // }}}
        std::string file_name_z = "./dat/EigenValue_L"+std::to_string(valley+1)+"_z_"+std::to_string(bandsL)+"bands.dat"; // z {{{
        std::ofstream file_z(file_name_z);
        for (int i=-n; i<=n; i++) {
            double kx = cutoff*double(i)/double(n);
            double k[3] = {0e0,0e0,kx};
            matrixComplex H_L = set_L(valley, k);

            vectorReal E_L = diagonalize_N(H_L);

            file_z << std::scientific << kx << ", ";
            for(int i=0; i<bandsL; i++) {
                file_z << std::scientific << E_L[i];
                if (i != bandsL-1) file_z << ", ";
            }
            file_z << std::endl;

        }; // }}}
    };
// }}}

//    int valley = 0;
//    int band_index = 2;
//    chemical_potential mu = 0.1e-2;
//    triangles tri = get_triangles_T(band_index, mu);
//    double dos = get_DOS_T(tri, band_index, mu);
//    std::cout << std::scientific << tri.faces.size() << ", " << dos << std::endl;
////    std::string dos_file = "dos"+std::to_string(band_index)+".csv";
////    std::ofstream ofs(dos_file, std::ios::app);
////    ofs << std::scientific << tri.faces.size() << ", " << dos << std::endl;
////    std::string filename = "tri_k"+std::to_string(k_mesh)+".csv";
////    triangles_write(tri, filename);

    int band_index = 2;
    chemical_potential mu_max = double(ET[band_index])+1e-2;
    chemical_potential mu_min = double(ET[band_index])-1e-1;
    band bL;
    bL = set_band_T(band_index, mu_min, mu_max, mu_mesh);

//    int band_index = 2;
////    double delta = (EL[valley][2] - EL[valley][0])*5e-1;
////    chemical_potential mu_min = delta - 1e-2;
////    chemical_potential mu_max = delta + 8e-1;
//    chemical_potential mu_max = double(EL[valley][band_index])+4e-1;
//    chemical_potential mu_min = double(EL[valley][band_index])-1e-4;
//    band bL;
//    bL = set_band_L(valley, band_index, mu_min, mu_max, mu_mesh);

//    sys_T T = get_T();
//    sys_T_write(T);

//    sys_L L = get_L();
//    sys_L_write(L);

    return 0;
}
