#include "parameters.hpp"
#include <filesystem>

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

    for(int i=1; i<7; i++) {
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
                for(int axis=0; axis<space_dim; axis++) {
                    matrixComplex C1 = product(vT[external], GR);
                    matrixComplex C2 = product(vT[axis], GA);

                    res[external][axis] = tr(product(C1, C2)) * epsilon;
                }
            }

            return res;
        };

        std::string filename = "dat/"+dir+"/conductivity_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofs(filename);
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

    for(int i=3; i<4; i++) {
        double epsilon = 1e0;
        for(int j=0; j<i; j++) {
            epsilon *= 1e-1;
        }
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        auto fn = [](int band_index, chemical_potential mu, kpoint k) {
            SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
            for(int external=0; external<space_dim; external++) {
                for(int axis=0; axis<space_dim; axis++) {
                    for(int spin=0; spin<spin_dim; spin++) {
                        Green_function GR = get_green_function_T(mu + eps_num*zi, k);
                        Green_function GA = get_green_function_T(mu - eps_num*zi, k);

                        matrixComplex vsGR   = product(v_s_T[axis][spin], GR);
//                        matrixComplex vsGA   = product(v_s_T[axis][spin], GA);
//                        matrixComplex vGR    = product(vT[external], GR);
                        matrixComplex vGA    = product(vT[external], GA);

//                        res[external][axis][spin] = tr(product(vsGR, vGA)) - 5e-1*(tr(product(vsGR, vGR)) + tr(product(vsGA, vGA)));
                        res[external][axis][spin] = tr(product(vsGR, vGA));
                    }
                }
            }
            return res;
        };

        std::string filename = "dat/"+dir+"/spin_Hall_conductivity1_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofs(filename);
        ofs << "# mu";
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

    for(int i=3; i<4; i++) {
        double epsilon = 1e0;
        for(int j=0; j<i; j++) {
            epsilon *= 1e-1;
        }
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        auto fn = [](int band_index, chemical_potential mu, kpoint k) {
            SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
            for(int external=0; external<space_dim; external++) {
                for(int axis=0; axis<space_dim; axis++) {
                    for(int spin=0; spin<spin_dim; spin++) {
                        Green_function GR = get_green_function_T(mu + eps_num*zi, k);
                        Green_function GA = get_green_function_T(mu - eps_num*zi, k);

                        matrixComplex GRGR = product(GR, GR);
                        matrixComplex GAGA = product(GA, GA);

                        matrixComplex vsGR   = product(v_s_T[axis][spin], GR);
                        matrixComplex vsGRGR = product(v_s_T[axis][spin], GRGR);
                        matrixComplex vsGA   = product(v_s_T[axis][spin], GA);
                        matrixComplex vsGAGA = product(v_s_T[axis][spin], GAGA);
                        matrixComplex vGR    = product(vT[external], GR);
                        matrixComplex vGRGR  = product(vT[external], GRGR);
                        matrixComplex vGA    = product(vT[external], GA);
                        matrixComplex vGAGA  = product(vT[external], GAGA);

                        Complex c = tr(product(vsGRGR, vGR)) - tr(product(vsGR, vGRGR)) - tr(product(vsGAGA, vGA)) + tr(product(vsGA, vGAGA));
                        res[external][axis][spin] = - c;
                    }
                }
            }
            return res;
        };

        std::string filename = "dat/"+dir+"/spin_Hall_conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofs(filename);
        ofs << "# mu";
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
Conductivity get_conductivity_L(band b, int valley) { // {{{
    Conductivity response;
    response.resize(space_dim); // external index;
    for(int external=0; external<space_dim; external++) {
        response[external].resize(space_dim); // response index;
        for(int axis=0; axis<space_dim; axis++) {
            response[external][axis] = 0e0;
        }
    }

    std::string dir = "L"+std::to_string(valley)+"_"+std::to_string(bandsL)+"bands";
    set_output_directory(dir);

    std::vector<std::thread> threads;
    threads.resize(thread_num);

    for(int i=1; i<7; i++) {
        double epsilon = 1e0;
        for(int j=0; j<i; j++) {
            epsilon *= 1e-1;
        }
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        auto fn = [=](int band_index, chemical_potential mu, int valley, kpoint k) {
            Green_function GR = get_green_function_L(mu + epsilon*zi, valley, k);
            Green_function GA = get_green_function_L(mu - epsilon*zi, valley, k);

            matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));

            for(int external=0; external<space_dim; external++) {
                for(int axis=0; axis<space_dim; axis++) {
                    matrixComplex C1 = product(vL[valley][external], GR);
                    matrixComplex C2 = product(vL[valley][axis], GA);

                    res[external][axis] = tr(product(C1, C2)) * epsilon;
                }
            }

            return res;
        };

        std::string filename = "dat/"+dir+"/conductivity_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofs(filename);
        for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
            matrixComplex sigma(space_dim, vectorComplex(space_dim, 0e0));
            chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
            integrate_triangles_L(fn, sigma, b.tri[i_mu], valley, b.index, mu);
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

SHC get_SHC_L(band b, int valley) { // {{{
    SHC response;
    response.resize(space_dim); // external index;
    for(int external=0; external<space_dim; external++) {
        response[external].resize(space_dim);
        for(int axis=0; axis<space_dim; axis++) {
            response[external][axis].resize(spin_dim);
            for(int spin=0; spin<spin_dim; spin++) {
                response[external][axis][spin] = 0e0;
            }
        }
    }

    auto fn = [](int valley, int band_index, chemical_potential mu, kpoint k) {
        int external = 1; // x
        int axis = 1;     // y
        int spin = 2;     // z
        Green_function GR = get_green_function_L(mu + eps_num*zi, valley, k);
        Green_function GA = get_green_function_L(mu - eps_num*zi, valley, k);

        Complex res = 0e0;

        matrixComplex C1 = product(vL[valley][external], GR);
        matrixComplex C2 = product(vL[valley][axis], GA);

        res = tr(product(C1, C2));

        return res;
    };

    for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
        Complex dos = 0e0;
        chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
        integrate_triangles_L(fn, dos, b.tri[i_mu], valley, b.index, mu);
    }

    return response;
}; // }}}
// }}}
