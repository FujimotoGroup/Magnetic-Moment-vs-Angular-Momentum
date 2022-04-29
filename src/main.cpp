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
        file << std::scientific << "b1x = " << b1[0]  << std::endl;
        file << std::scientific << "b1y = " << b1[1]  << std::endl;
        file << std::scientific << "b1z = " << b1[2]  << std::endl;
        file << std::scientific << "b2x = " << b2[0]  << std::endl;
        file << std::scientific << "b2y = " << b2[1]  << std::endl;
        file << std::scientific << "b2z = " << b2[2]  << std::endl;
        file << std::scientific << "b3x = " << b3[0]  << std::endl;
        file << std::scientific << "b3y = " << b3[1]  << std::endl;
        file << std::scientific << "b3z = " << b3[2]  << std::endl;
        file << std::endl;
        file << std::scientific << "kTx = " << kT[0] << std::endl;
        file << std::scientific << "kTy = " << kT[1] << std::endl;
        file << std::scientific << "kTz = " << kT[2] << std::endl;
        file << std::scientific << "kL1x = " << kL[0][0] << std::endl;
        file << std::scientific << "kL1y = " << kL[0][1] << std::endl;
        file << std::scientific << "kL1z = " << kL[0][2] << std::endl;
        file << std::scientific << "kL2x = " << kL[1][0] << std::endl;
        file << std::scientific << "kL2y = " << kL[1][1] << std::endl;
        file << std::scientific << "kL2z = " << kL[1][2] << std::endl;
        file << std::scientific << "kL3x = " << kL[2][0] << std::endl;
        file << std::scientific << "kL3y = " << kL[2][1] << std::endl;
        file << std::scientific << "kL3z = " << kL[2][2] << std::endl;
        file << std::endl;
        file << std::scientific << "bands  = " <<  bands  << std::endl;
        file << std::scientific << "bandsT = " <<  bandsT << std::endl;
        file << std::scientific << "bandsL = " <<  bandsL << std::endl;
        file << std::scientific << "lowest_band_T = " << lowest_band_T << std::endl;
        file << std::scientific << "lowest_band_L = " << lowest_band_L << std::endl;
        file << std::endl;
        for(int i=0; i<bandsT; i++)
        file << std::scientific << "ET" << i << " = " <<  ET[i]  << std::endl;
        file << std::endl;
        for(int i=0; i<bandsL; i++)
        file << std::scientific << "EL" << i << " = " <<  EL[0][i]  << std::endl;
        file << std::endl;
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
//// get triangles (vertex and face) {{{
//    int band_index = 2;
//    chemical_potential mu = 0e0;
//    triangles tri = get_triangles_T(band_index, mu);
//    std::string name = "./dat/triangle_T";
//    triangles_write_T(tri, name);
//    for(int valley=0; valley<valleys; valley++) {
//        int band_index = 2;
//        chemical_potential mu = 0e0;
//        triangles tri = get_triangles_L(valley, band_index, mu);
//        std::string name = "./dat/triangle_L"+std::to_string(valley+1);
//        triangles_write_L(tri, name, valley);
//    }
//// }}}

//// isotropic Dirac (required modification in lib/hamiltonian.cpp) {{{
//    set_isotropic();
//    int valley = 0;
//    int band_index = 0;
//    double delta = (EL[valley][2] - EL[valley][0])*5e-1;
//    mu_cutoff_L = - 1e-1;
//    mu_cutoff_mesh_L = 100;
//    band_edge_L[valley][0] = -delta; band_edge_L[valley][1] = -delta;
//    band_edge_L[valley][2] =  delta; band_edge_L[valley][3] =  delta;
//    chemical_potential mu_min = - delta - 5e-2;
//    chemical_potential mu_max = 0e0;
//    set_response_L(mu_min, mu_max, mu_mesh, valley, band_index);
//
//    band_index = 2;
//    mu_cutoff_L = 1e-1;
//    mu_cutoff_mesh_L = 100;
//    mu_min = 0e0;
//    mu_max = delta + 5e-2;
//    set_response_L(mu_min, mu_max, mu_mesh, valley, band_index);
//// }}}

//    int band_index = 2;
//    chemical_potential mu = 0e0;
//    int valley = 0;
//    chemical_potential mu = 0e0;
//    chemical_potential mu = -1.4e-1;
//    triangles tri = get_triangles_T(band_index, mu);
//    double dos = get_DOS_T(tri, band_index, mu);
//    std::string name = "./triangle_T";
//    triangles_write_T(tri, name);
//    std::cout << std::setprecision(15) << tri.faces.size() << ", " << dos << std::endl;
//
//    triangles tri = get_triangles_L(valley, band_index, mu);
//    double dos = get_DOS_L(tri, valley, band_index, mu);
//    std::string name = "./triangle_L1";
//    triangles_write_L(tri, name, valley);
//    std::cout << std::setprecision(15) << tri.faces.size() << ", " << dos << std::endl;
////    std::string dos_file = "dos"+std::to_string(band_index)+".csv";
////    std::ofstream ofs(dos_file, std::ios::app);
////    ofs << std::scientific << tri.faces.size() << ", " << dos << std::endl;
////    std::string filename = "tri_k"+std::to_string(k_mesh)+".csv";
////    triangles_write(tri, filename);
//
//    int band_index = 5;
//    chemical_potential mu_max = double(ET[band_index])-5e-5;
////    chemical_potential mu_min = double(ET[band_index])-2.5e-1;
////    chemical_potential mu_max = -0.5e-1;
//    chemical_potential mu_min = -1.4e-1;
//    band bT;
//    bT = set_band_T(band_index, mu_min, mu_max, mu_mesh);
//    Conductivity sigma_T = get_conductivity_T(bT);
//    SHC SHC1_T = get_SHC_T1(bT);
//    SHC SHC2_T = get_SHC_T2(bT);

    int valley = 0;
    int band_index = 6;
//    chemical_potential mu = 0e0;
////    for(int valley=0; valley<valleys; valley++) {
        mu_cutoff_L = 8e-1;
        mu_cutoff_mesh_L = 100;
        chemical_potential mu_max = double(EL[valley][band_index])+4e-1;
        chemical_potential mu_min = double(EL[valley][band_index])-1e-2;
        set_response_L(mu_min, mu_max, mu_mesh, valley, band_index);
////    }

//    sys_T T = get_T();
//    sys_T_write(T);

//    sys_L L = get_L();
//    sys_L_write(L);
//

    return 0;
}
