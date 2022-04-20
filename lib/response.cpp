#include "parameters.hpp"
#include <filesystem>

int i_epsilon_min = 4;
int i_epsilon_max = 5;

void set_output_directory(std::string dir) { // {{{
    namespace fs = std::filesystem;
    std::string directory_name("dat/"+dir);
    fs::create_directory(directory_name)?
         std::cout << "created directory - " << directory_name << std::endl :
         std::cout << "create_directory() failed" << std::endl;
}; // }}}

// T {{{
Conductivity get_conductivity_T(band b) { // {{{
    Conductivity response;
    response.resize(space_dim); // external index;
    for(int external=0; external<space_dim; external++) {
        response[external].resize(space_dim); // response index;
        for(int axis=0; axis<space_dim; axis++) {
            response[external][axis] = 0e0;
        }
    }

    std::string dir = "T"+std::to_string(bandsT)+"bands";
    set_output_directory(dir);

    for(int i=i_epsilon_min; i<i_epsilon_max; i++) {
        double epsilon = 1e0;
        for(int j=0; j<i; j++) {
            epsilon *= 1e-1;
        }
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        auto fn = [=](int band_index, chemical_potential mu, kpoint k) {
            Green_function GR = get_green_function_T(mu + epsilon*zi, k);
            Green_function GA = get_green_function_T(mu - epsilon*zi, k);

            matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));

            for(int external=0; external<space_dim; external++) {
                matrixComplex C1 = product(vT[external], GR);
                for(int axis=0; axis<space_dim; axis++) {
                    matrixComplex C2 = product(vT[axis], GA);

                    res[external][axis] = tr(product(C1, C2));
                }
            }

            return res;
        };

        std::string filename = "dat/"+dir+"/conductivity_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofs(filename);
        ofs << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                ofs << ", " << axises[external]+axises[axis];
            }
        }
        for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
            matrixComplex sigma(space_dim, vectorComplex(space_dim, 0e0));
            chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
            integrate_triangles_T(fn, sigma, b.tri[i_mu], b.index, mu);
            ofs << std::scientific << mu;
            for( auto s : sigma ) {
                for( auto v : s ) {
                    ofs << std::scientific << ", " << v.real();
                }
            }
            ofs << std::endl;
        }
    }

    return response;
}; // }}}

SHC get_SHC_T1(band b) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    std::string dir = "T"+std::to_string(bandsT)+"bands";
    set_output_directory(dir);

    for(int i=i_epsilon_min; i<i_epsilon_max; i++) {
        double epsilon = 1e0;
        for(int j=0; j<i; j++) {
            epsilon *= 1e-1;
        }
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        auto fn = [](int band_index, chemical_potential mu, kpoint k) {
            SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
            Green_function GR = get_green_function_T(mu + eps_num*zi, k);
            Green_function GA = get_green_function_T(mu - eps_num*zi, k);
            for(int external=0; external<space_dim; external++) {
//                matrixComplex vGR    = product(vT[external], GR);
                matrixComplex vGA    = product(vT[external], GA);
                for(int axis=0; axis<space_dim; axis++) {
                    for(int spin=0; spin<spin_dim; spin++) {
                        matrixComplex vsGR   = product(v_s_T[axis][spin], GR);
//                        matrixComplex vsGA   = product(v_s_T[axis][spin], GA);

//                        res[external][axis][spin] = tr(product(vsGR, vGA)) - 5e-1*(tr(product(vsGR, vGR)) + tr(product(vsGA, vGA)));
                        res[external][axis][spin] = tr(product(vsGR, vGA));
                    }
                }
            }
            return res;
        };

        std::string filename = "dat/"+dir+"/spin_Hall_conductivity1_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofs(filename);
        ofs << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofs << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofs << std::endl;
        for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
            SHC sigma(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
            chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
            integrate_triangles_T(fn, sigma, b.tri[i_mu], b.index, mu);
            ofs << std::scientific << mu;
            for( auto e : sigma ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofs << std::scientific << ", " << s.real();
                    }
                }
            }
            ofs << std::endl;
        }
    }

    return response;
}; // }}}

SHC get_SHC_T2(band b) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    std::string dir = "T"+std::to_string(bandsT)+"bands";
    set_output_directory(dir);

    for(int i=i_epsilon_min; i<i_epsilon_max; i++) {
        double epsilon = 1e0;
        for(int j=0; j<i; j++) {
            epsilon *= 1e-1;
        }
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        auto fn = [](int band_index, chemical_potential mu, kpoint k) {
            SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
            Green_function GR = get_green_function_T(mu + eps_num*zi, k);
            Green_function GA = get_green_function_T(mu - eps_num*zi, k);

            for(int external=0; external<space_dim; external++) {
                matrixComplex vGR = product(vT[external], GR);
                matrixComplex vGA = product(vT[external], GA);
                matrixComplex GRv = product(GR, vT[external]);
                matrixComplex GAv = product(GA, vT[external]);
                matrixComplex QR  = minus(GRv, vGR);
                matrixComplex QA  = minus(GAv, vGA);
                matrixComplex vR  = product(GR, product(QR, GR));
                matrixComplex vA  = product(GA, product(QA, GA));
                for(int axis=0; axis<space_dim; axis++) {
                    for(int spin=0; spin<spin_dim; spin++) {
                        matrixComplex R = product(v_s_T[axis][spin], vR);
                        matrixComplex A = product(v_s_T[axis][spin], vA);

                        res[external][axis][spin] = tr(R);
//                        res[external][axis][spin] = tr(R) - tr(A);
                    }
                }
            }
            return res;
        };

        std::string filename = "dat/"+dir+"/spin_Hall_conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofs(filename);
        ofs << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofs << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofs << std::endl;
        for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
            SHC sigma(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
            chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
            integrate_triangles_T(fn, sigma, b.tri[i_mu], b.index, mu);
            ofs << std::scientific << mu;
            for( auto e : sigma ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofs << std::scientific << ", " << s.real();
                    }
                }
            }
            ofs << std::endl;
        }
    }

    return response;
}; // }}}
// }}}

// L {{{
//
//Conductivity get_conductivity_L(band b, int valley) { // {{{
//    Conductivity response(space_dim, vectorComplex(space_dim, 0e0));
//
//    std::string dir = "L"+std::to_string(valley+1)+"_"+std::to_string(bandsL)+"bands";
//    set_output_directory(dir);
//
//    for(int i=i_epsilon_min; i<i_epsilon_max; i++) {
//        double epsilon = 1e0;
//        for(int j=0; j<i; j++) {
//            epsilon *= 1e-1;
//        }
//        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;
//
//        auto fn = [=](int valley, int band_index, chemical_potential mu, kpoint k) {
//            Green_function GR = get_green_function_L(mu + epsilon*zi, valley, k);
//            Green_function GA = get_green_function_L(mu - epsilon*zi, valley, k);
//
//            matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));
//
//            for(int external=0; external<space_dim; external++) {
//                matrixComplex C1 = product(vL[valley][external], GR);
//                for(int axis=0; axis<space_dim; axis++) {
//                    matrixComplex C2 = product(vL[valley][axis], GA);
//
//                    res[external][axis] = tr(product(C1, C2));
//                }
//            }
//
//            return res;
//        };
//
//        std::string filename = "dat/"+dir+"/conductivity_eps"+std::to_string(epsilon)+".csv";
//        std::ofstream ofs(filename);
//        ofs << "mu";
//        for(int external=0; external<space_dim; external++) {
//            for(int axis=0; axis<space_dim; axis++) {
//                ofs << ", " << axises[external]+axises[axis];
//            }
//        }
//        filename = "dat/"+dir+"/conductivity_imag_eps"+std::to_string(epsilon)+".csv";
//        std::ofstream ofss(filename);
//        ofss << "mu";
//        for(int external=0; external<space_dim; external++) {
//            for(int axis=0; axis<space_dim; axis++) {
//                ofss << ", " << axises[external]+axises[axis];
//            }
//        }
//        for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
//            matrixComplex sigma(space_dim, vectorComplex(space_dim, 0e0));
//            chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
//            integrate_triangles_L(fn, sigma, b.tri[i_mu], valley, b.index, mu);
//            ofs << std::scientific << mu;
//            for( auto s : sigma ) {
//                for( auto v : s ) {
//                    ofs << std::scientific << ", " << v.real();
//                }
//            }
//            ofs << std::endl;
//            ofss << std::scientific << mu;
//            for( auto s : sigma ) {
//                for( auto v : s ) {
//                    ofss << std::scientific << ", " << v.imag();
//                }
//            }
//            ofss << std::endl;
//        }
//    }
//
//    return response;
//}; // }}}
//
void set_response_L(chemical_potential ene_min, chemical_potential ene_max, int ene_mesh, int valley, int band_index) { // {{{
    std::string dir = "L"+std::to_string(valley+1)+"_"+std::to_string(bandsL)+"bands";
    set_output_directory(dir);

    for(int i=i_epsilon_min; i<i_epsilon_max; i++) {
        double epsilon = 1e0;
        for(int j=0; j<i; j++) {
            epsilon *= 1e-1;
        }
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        std::string filename;
// init dos file {{{
        filename = "dat/"+dir+"/dos_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofdos(filename);
// }}}
// init conductivity file {{{
        filename = "dat/"+dir+"/conductivity_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofsigma(filename);
        ofsigma << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                ofsigma << ", " << axises[external]+axises[axis];
            }
        }
        ofsigma << std::endl;
// }}}
// init spin Hall conductivity 1 file {{{
        filename = "dat/"+dir+"/spin_Hall_conductivity1_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc1(filename);
        ofshc1 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc1 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc1 << std::endl;
// }}}

        chemical_potential d_ene = (ene_max - ene_min) / double(ene_mesh-1);
//        for(int i_ene=0; i_ene<ene_mesh; i_ene++) {
        for(int i_ene=0; i_ene<1; i_ene++) {
            const chemical_potential mu = ene_min + d_ene*double(i_ene);

            chemical_potential e_min = mu - epsilon*5e0;
            chemical_potential e_max = mu + epsilon*5e0;
            int e_mesh = 51;

            band bL;
            bL = set_band_L(valley, band_index, e_min, e_max, e_mesh);

// dos output {{{
            for(int i_e=0; i_e<e_mesh; i_e++) {
                ofdos << std::scientific << bL.tri[i_e].ene << ", " << bL.dos[i_e] << std::endl;
            }
// }}}

            Conductivity sigma(space_dim, vectorComplex(space_dim, 0e0));
            sigma = get_conductivity_L(bL, mu, valley);
// conductivity output {{{
            ofsigma << std::scientific << mu;
            for( auto s : sigma ) {
                for( auto v : s ) {
                    ofsigma << std::scientific << ", " << v.real();
                }
            }
            ofsigma << std::endl;
// }}}

            SHC SHC1(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
            SHC1 = get_SHC_L1(bL, mu, valley);
// SHC1 output {{{
            ofshc1 << std::scientific << mu;
            for( auto e : SHC1 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc1 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc1 << std::endl;
// }}}
        }
    }
}; // }}}

Conductivity get_conductivity_L(band b, chemical_potential mu, int valley) { // {{{
    Conductivity response(space_dim, vectorComplex(space_dim, 0e0));

    double epsilon = b.dmu;

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_L(e + epsilon*zi, valley, k);
        Green_function GA = get_green_function_L(e - epsilon*zi, valley, k);

        matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));

        for(int external=0; external<space_dim; external++) {
//            matrixComplex C1 = product(vL[valley][external], GR);
            matrixComplex C1 = GR;
            for(int axis=0; axis<space_dim; axis++) {
//                matrixComplex C2 = product(vL[valley][axis], GA);
                matrixComplex C2 = GA;

                res[external][axis] = tr(product(C1, C2));
            }
        }

        return res;
    };

    std::string filename = "conductivity.csv";
    std::ofstream ofsigma(filename, std::ios::app);
    for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
        matrixComplex sigma(space_dim, vectorComplex(space_dim, 0e0));
        integrate_triangles_L(fn, sigma, b.tri[i_mu], valley, b.index, mu);
        ofsigma << b.tri[i_mu].ene - mu << ", " << sigma[0][0].real() << ", " << sigma[0][0].imag() << std::endl;
        sigma = times(sigma, b.dmu);
        response = add(response, sigma);
    }

    return response;
}; // }}}

SHC get_SHC_L1(band b, chemical_potential mu, int valley) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    double epsilon = b.dmu;

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_L(e + eps_num*zi, valley, k);
        Green_function GA = get_green_function_L(e - eps_num*zi, valley, k);
        for(int external=0; external<space_dim; external++) {
//            matrixComplex vGR    = product(vL[valley][external], GR);
            matrixComplex vGA    = product(vL[valley][external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    matrixComplex vsGR   = product(v_s_L[valley][axis][spin], GR);
//                    matrixComplex vsGA   = product(v_s_L[valley][axis][spin], GA);

//                    res[external][axis][spin] = tr(product(vsGR, vGA)) - 5e-1*(tr(product(vsGR, vGR)) + tr(product(vsGA, vGA)));
                    res[external][axis][spin] = tr(product(vsGR, vGA));
                }
            }
        }
        return res;
    };

    for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
        SHC sigma(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        integrate_triangles_L(fn, sigma, b.tri[i_mu], valley, b.index, mu);
        sigma = times(sigma, b.dmu);
        response = add(response, sigma);
    }

    return response;
}; // }}}
// }}}
