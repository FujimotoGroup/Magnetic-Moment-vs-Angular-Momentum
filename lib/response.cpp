#include "parameters.hpp"

void set_output_directory(std::string dir) { // {{{
    namespace fs = std::filesystem;
    std::string directory_name("dat/"+dir);
    fs::create_directories(directory_name)?
         std::cout << "created directory - " << directory_name << std::endl :
         std::cout << "create_directory() failed" << std::endl;
}; // }}}

// T {{{
void set_response_T(chemical_potential ene_min, chemical_potential ene_max, int ene_mesh, int band_index) { // {{{
    std::string dir = "T_"+std::to_string(bandsT)+"bands/band_index"+std::to_string(band_index);
    set_output_directory(dir);

    for(int i=5; i<6; i++) {
        double epsilon = double(i)*1e-4;
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        chemical_potential d_ene = (ene_max - ene_min) / double(ene_mesh-1);

        band b;
        band b_main = set_band_T(band_index, ene_min, ene_max, ene_mesh);
        if (mu_cutoff_T < ene_min) {
            band b_cut = set_band_T(band_index, mu_cutoff_T, ene_min-d_ene, mu_cutoff_mesh_T);
            b = combine_band(b_cut, b_main);
        } else if (ene_max < mu_cutoff_T) {
            band b_cut = set_band_T(band_index, ene_max+d_ene, mu_cutoff_T, mu_cutoff_mesh_T);
            b = combine_band(b_main, b_cut);
        } else {
            std::cerr << "mu_cutoff_T shoudl be smaller than ene_min or larger than ene_max" << std::endl;
            exit(0);
        }
        if ( (ene_min < band_edge_T[band_index]) & (band_edge_T[band_index] < ene_max) ) {
            int e_mesh = 15;
            double e_cut = 1.9e0*epsilon;
            double power = 7e-1;
            band b_edge = set_band_2n_T(band_index, band_edge_T[band_index], e_cut, e_mesh, power);
            b_main = combine_band_2n(b_main, b_edge);

        }

        std::vector<Conductivity> sigma_e;
        std::vector<SHC> sigma_m1, sigma_m2;
        std::vector<SHC> sigma_a1, sigma_a2;
        std::vector<SHC> sigma_r1, sigma_r2;
        sigma_e.resize(b_main.mesh);
        sigma_m1.resize(b_main.mesh);
        sigma_m2.resize(b_main.mesh);
        sigma_a1.resize(b_main.mesh);
        sigma_a2.resize(b_main.mesh);
        sigma_r1.resize(b_main.mesh);
        sigma_r2.resize(b_main.mesh);

        std::string filename;
// init dos file {{{
        filename = "dat/"+dir+"/dos_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofdos(filename);
// }}}
// init conductivity file {{{
        filename = "dat/"+dir+"/sigma_e_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofsigma1(filename);
        ofsigma1 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                ofsigma1 << ", " << axises[external]+axises[axis];
            }
        }
        ofsigma1 << std::endl;

        filename = "dat/"+dir+"/sigma_e_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofsigma2(filename);
        ofsigma2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                ofsigma2 << ", " << axises[external]+axises[axis];
            }
        }
        ofsigma2 << std::endl;
// }}}
// init spin magnetic Hall conductivity 1 file {{{
        filename = "dat/"+dir+"/sigma_m1_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc11(filename);
        ofshc11 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc11 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc11 << std::endl;

        filename = "dat/"+dir+"/sigma_m1_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc12(filename);
        ofshc12 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc12 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc12 << std::endl;
// }}}
// init spin magnetic Hall conductivity 2 file {{{
        filename = "dat/"+dir+"/sigma_m2_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc21(filename);
        ofshc21 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc21 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc21 << std::endl;

        filename = "dat/"+dir+"/sigma_m2_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc22(filename);
        ofshc22 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc22 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc22 << std::endl;
// }}}
// init spin angular Hall conductivity 1 file {{{
        filename = "dat/"+dir+"/sigma_a1_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc11(filename);
        of_angular_shc11 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc11 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc11 << std::endl;

        filename = "dat/"+dir+"/sigma_a1_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc12(filename);
        of_angular_shc12 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc12 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc12 << std::endl;
// }}}
// init spin angular Hall conductivity 2 file {{{
        filename = "dat/"+dir+"/sigma_a2_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc21(filename);
        of_angular_shc21 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc21 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc21 << std::endl;

        filename = "dat/"+dir+"/sigma_a2_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc22(filename);
        of_angular_shc22 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc22 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc22 << std::endl;
// }}}

        int i_ene;
        for(i_ene=0; i_ene<b_main.mesh; i_ene++) {
            const chemical_potential mu = b_main.ene[i_ene];
            std::cout << "start: i_ene = " << i_ene << ": mu = " << mu << std::endl;

// dos output {{{
            ofdos << std::scientific << mu << ", " << b_main.dos[i_ene] / (angstrom*angstrom*angstrom) << std::endl;
// }}}

            int e_mesh = 47;
            double e_cut = 59e0*epsilon;
            double power = 9e-1;

            band bT = set_band_2n_T(band_index, mu, e_cut, e_mesh, power);

            band b_sum = combine_band_2n(b, bT);

            sigma_e[i_ene] = get_conductivity_T(b_sum, epsilon, mu);
// conductivity output {{{
            ofsigma1 << std::scientific << mu;
            ofsigma2 << std::scientific << mu;
            for( auto s : sigma_e[i_ene] ) {
                for( auto v : s ) {
                    ofsigma1 << std::scientific << ", " << v.real();
                    ofsigma2 << std::scientific << ", " << v.imag();
                }
            }
            ofsigma1 << std::endl;
            ofsigma2 << std::endl;
// }}}

            sigma_m1[i_ene] = get_SHC_T1(b_sum, epsilon, mu);
// magnetic SHC1 output {{{
            ofshc11 << std::scientific << mu;
            ofshc12 << std::scientific << mu;
            for( auto e : sigma_m1[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc11 << std::scientific << ", " << s.real();
                        ofshc12 << std::scientific << ", " << s.imag();
                    }
                }
            }
            ofshc11 << std::endl;
            ofshc12 << std::endl;
// }}}

            sigma_m2[i_ene] = get_SHC_T2(b_sum, epsilon, mu);
// magnetic SHC2 output {{{
            ofshc21 << std::scientific << mu;
            ofshc22 << std::scientific << mu;
            for( auto e : sigma_m2[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc21 << std::scientific << ", " << s.real();
                        ofshc22 << std::scientific << ", " << s.imag();
                    }
                }
            }
            ofshc21 << std::endl;
            ofshc22 << std::endl;
// }}}

            sigma_a1[i_ene] = get_angular_SHC_T1(b_sum, epsilon, mu);
// angular SHC1 output {{{
            of_angular_shc11 << std::scientific << mu;
            of_angular_shc12 << std::scientific << mu;
            for( auto e : sigma_a1[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc11 << std::scientific << ", " << s.real();
                        of_angular_shc12 << std::scientific << ", " << s.imag();
                    }
                }
            }
            of_angular_shc11 << std::endl;
            of_angular_shc12 << std::endl;
// }}}

            sigma_a2[i_ene] = get_angular_SHC_T2(b_sum, epsilon, mu);
// angular SHC2 output {{{
            of_angular_shc21 << std::scientific << mu;
            of_angular_shc22 << std::scientific << mu;
            for( auto e : sigma_a2[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc21 << std::scientific << ", " << s.real();
                        of_angular_shc22 << std::scientific << ", " << s.imag();
                    }
                }
            }
            of_angular_shc21 << std::endl;
            of_angular_shc22 << std::endl;
// }}}
        }

        Conductivity Sigma(space_dim, vectorComplex(space_dim, 0e0));
        SHC magnetic_SHC1(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        SHC magnetic_SHC2(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        SHC angular_SHC1(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        SHC angular_SHC2(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

// mu dependence @ T = 0 {{{
        dir = dir+"/mu-dependence";
        set_output_directory(dir);
// init Sigma file {{{
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
// init magnetic SHC1 file {{{
        filename = "dat/"+dir+"/spin-magnetic-conductivity1_eps"+std::to_string(epsilon)+".csv";
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
// init magnetic SHC2 file {{{
        filename = "dat/"+dir+"/spin-magnetic-conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc2(filename);
        ofshc2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc2 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc2 << std::endl;
// }}}
// init angular SHC1 file {{{
        filename = "dat/"+dir+"/spin-angular-conductivity1_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc1(filename);
        of_angular_shc1 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc1 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc1 << std::endl;
// }}}
// init angular SHC2 file {{{
        filename = "dat/"+dir+"/spin-angular-conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc2(filename);
        of_angular_shc2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc2 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc2 << std::endl;
// }}}
// init real SHC1 file {{{
        filename = "dat/"+dir+"/real-spin-angular-conductivity1_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_real_shc1(filename);
        of_real_shc1 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_real_shc1 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_real_shc1 << std::endl;
// }}}
// init real SHC2 file {{{
        filename = "dat/"+dir+"/real-spin-angular-conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_real_shc2(filename);
        of_real_shc2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_real_shc2 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_real_shc2 << std::endl;
// }}}

        chemical_potential mu;
        for(int i_ene=0; i_ene<b_main.mesh; i_ene++) {
            mu = b_main.ene[i_ene];
            Sigma = sigma_e[i_ene];
// Sigma output {{{
            ofsigma << std::scientific << mu;
            for( auto s : Sigma ) {
                for( auto v : s ) {
                    ofsigma << std::scientific << ", " << v.real();
                }
            }
            ofsigma << std::endl;
// }}}
            magnetic_SHC1 = sigma_m1[i_ene];
// magnetic SHC1 output {{{
            ofshc1 << std::scientific << mu;
            for( auto e : magnetic_SHC1 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc1 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc1 << std::endl;
// }}}
            angular_SHC1 = sigma_a1[i_ene];
// angular SHC1 output {{{
            of_angular_shc1 << std::scientific << mu;
            for( auto e : angular_SHC1 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc1 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc1 << std::endl;
// }}}
        }

        double de;
        SHC magnetic_SHC2_mu;
        i_ene = 0;
            mu = b_main.ene[i_ene];
// magnetic SHC2 output {{{
            ofshc2 << std::scientific << mu;
            for( auto e : magnetic_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc2 << std::endl;
// }}}
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene])*5e-1;
            magnetic_SHC2_mu = times(sigma_m2[i_ene], de);
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
        for(int i_ene=1; i_ene<b_main.mesh-1; i_ene++) {
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene-1])*2.5e-1;
            magnetic_SHC2_mu = times(sigma_m2[i_ene], de);
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
// magnetic SHC2 output {{{
            ofshc2 << std::scientific << mu;
            for( auto e : magnetic_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc2 << std::endl;
// }}}
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
        }
        i_ene = b_main.mesh-1;
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene] - b_main.ene[i_ene-1])*5e-1;
            magnetic_SHC2_mu = times(sigma_m2[i_ene], de);
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
// magnetic SHC2 output {{{
            ofshc2 << std::scientific << mu;
            for( auto e : magnetic_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc2 << std::endl;
// }}}

        SHC angular_SHC2_mu;
        i_ene = 0;
            mu = b_main.ene[i_ene];
// angular SHC2 output {{{
            of_angular_shc2 << std::scientific << mu;
            for( auto e : angular_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc2 << std::endl;
// }}}
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene])*5e-1;
            angular_SHC2_mu = times(sigma_a2[i_ene], de);
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
        for(int i_ene=1; i_ene<b_main.mesh-1; i_ene++) {
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene-1])*2.5e-1;
            angular_SHC2_mu = times(sigma_a2[i_ene], de);
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
// angular SHC2 output {{{
            of_angular_shc2 << std::scientific << mu;
            for( auto e : angular_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc2 << std::endl;
// }}}
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
        }
        i_ene = b_main.mesh-1;
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene] - b_main.ene[i_ene-1])*5e-1;
            angular_SHC2_mu = times(sigma_a2[i_ene], de);
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
// angular SHC2 output {{{
            of_angular_shc2 << std::scientific << mu;
            for( auto e : angular_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc2 << std::endl;
// }}}
    }
// }}}
}; // }}}

Conductivity get_conductivity_T(band b, Energy epsilon, chemical_potential mu) { // {{{
    Conductivity response(space_dim, vectorComplex(space_dim, 0e0));

    auto fn = [=](int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_T(e + epsilon*zi, k);
        Green_function GA = get_green_function_T(e - epsilon*zi, k);

        matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));

        for(int external=0; external<space_dim; external++) {
            matrixComplex vGRe = product(vT[external], GR);
            matrixComplex vGAe = product(vT[external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                matrixComplex vGRa = product(vT[axis], GR);
                matrixComplex vGAa = product(vT[axis], GA);

                res[external][axis] = tr(product(vGRa, vGAe)) - 5e-1*(tr(product(vGRa, vGRe)) + tr(product(vGAa, vGAe)));
            }
        }

        return res;
    };

    integrate_band_T(fn, response, b, mu);

    double coef = (charge*v0)*(charge*v0) * hbar / (4e0*pi) / (angstrom*angstrom*angstrom) / charge;
    response = times(response, coef);

    return response;
}; // }}}

void set_conductivity_damping_dependence_at_Fermi_level_T(int band_index) { // {{{
    std::vector<Conductivity> sigma(damping.size());

    for (int i=0; i<damping.size(); i++) {
        double damping_constant = damping[i];
        std::cout << "damping constant = " << damping_constant << "[" << i << "]" <<std::endl;
        Energy epsilon = damping_constant;
        chemical_potential mu = 0e0;
        int e_mesh = 47;
        double e_cut = 59e0*epsilon;
        double power = 9e-1;

        band bT = set_band_2n_T(band_index, mu, e_cut, e_mesh, power);

        nu_F_T = bT.dos[e_mesh];
        double coef = damping_constant / nu_F_T;
        Self_energy se = get_self_energy_born_T(bT, 0e0, mu, epsilon, coef);
        se = add(product(impurityV1_T, product(se, impurityV1_T)), product(impurityV2_T, product(se, impurityV2_T)));
        sigma[i] = get_conductivity_with_self_energy_T(bT, 0e0, mu, se);
    }
// init Sigma file {{{
    std::string dir = "T"+std::to_string(bandsT)+"bands/damping-dependence";
    set_output_directory(dir);
    std::string filename = "dat/"+dir+"/conductivity_T.csv";
    std::ofstream ofsigma(filename);
    ofsigma << "damping";
    for(int external=0; external<space_dim; external++) {
        for(int axis=0; axis<space_dim; axis++) {
            ofsigma << ", " << axises[external]+axises[axis];
        }
    }
    ofsigma << std::endl;

    for (int i=0; i<damping.size(); i++) {
        ofsigma << std::scientific << damping[i];
        Conductivity s = sigma[i];
        for( auto v : s ) {
            for( auto r : v ) {
                ofsigma << std::scientific << ", " << r.real();
            }
        }
        ofsigma << std::endl;
    }
// }}}
}; // }}}

Conductivity get_conductivity_with_self_energy_T(band b, Energy epsilon, chemical_potential mu, Self_energy se) { // {{{
    Conductivity response(space_dim, vectorComplex(space_dim, 0e0));

    Self_energy SigmaR = se;
    Self_energy SigmaA(bandsT, vectorComplex(bandsT));
    for(int i=0; i<bandsT; i++) {
        SigmaA[i][i] = std::conj(SigmaR[i][i]);
        for(int j=i+1; j<bandsT; j++) {
            SigmaA[i][j] = std::conj(SigmaR[i][j]);
            SigmaA[j][i] = std::conj(SigmaR[j][i]);
        }
    }

    auto fn = [=](int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_full_green_function_T(e + eps_phys*zi, k, SigmaR);
        Green_function GA = get_full_green_function_T(e - eps_phys*zi, k, SigmaA);

        matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));

        for(int external=0; external<space_dim; external++) {
            matrixComplex vGRe = product(vT[external], GR);
            matrixComplex vGAe = product(vT[external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                matrixComplex vGRa = product(vT[axis], GR);
                matrixComplex vGAa = product(vT[axis], GA);

                res[external][axis] = tr(product(vGRa, vGAe)) - 5e-1*(tr(product(vGRa, vGRe)) + tr(product(vGAa, vGAe)));
            }
        }

        return res;
    };

    integrate_band_T(fn, response, b, mu);

    double coef = (charge*v0)*(charge*v0) * hbar / (4e0*pi) / (angstrom*angstrom*angstrom) / charge;
    response = times(response, coef);

    return response;
}; // }}}

SHC get_SHC_T1(band b, Energy epsilon, chemical_potential mu) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int band_index, chemical_potential e, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_T(e + epsilon*zi, k);
        Green_function GA = get_green_function_T(e - epsilon*zi, k);
        for(int external=0; external<space_dim; external++) {
//            matrixComplex vGR    = product(vT[external], GR);
            matrixComplex vGA    = product(vT[external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    matrixComplex vsGR   = product(v_s_T[axis][spin], GR);
//                    matrixComplex vsGA   = product(v_s_T[axis][spin], GA);

//                    res[external][axis][spin] = tr(product(vsGR, vGA)) - 5e-1*(tr(product(vsGR, vGR)) + tr(product(vsGA, vGA)));
                    res[external][axis][spin] = tr(product(vsGR, vGA));
                }
            }
        }
        return res;
    };

    integrate_band_T(fn, response, b, mu);

    double coef = - charge * v0*v0 * hbar / (2e0*pi) / (angstrom*angstrom*angstrom) / muBeV;
    response = times(response, coef);

    return response;
}; // }}}

SHC get_SHC_T2(band b, Energy epsilon, chemical_potential mu) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int band_index, chemical_potential mu, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_T(mu + epsilon*zi, k);
        Green_function GA = get_green_function_T(mu - epsilon*zi, k);

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

                    res[external][axis][spin] = tr(R) - tr(A);
                }
            }
        }
        return res;
    };

    integrate_band_T(fn, response, b, mu);

    double coef = - charge * v0*v0 * hbar / (4e0*pi) / (angstrom*angstrom*angstrom) / muBeV;
    response = times(response, coef);

    return response;
}; // }}}

SHC get_angular_SHC_T1(band b, Energy epsilon, chemical_potential mu) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int band_index, chemical_potential e, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_T(e + epsilon*zi, k);
        Green_function GA = get_green_function_T(e - epsilon*zi, k);
        for(int external=0; external<space_dim; external++) {
//            matrixComplex vGR    = product(vT[external], GR);
            matrixComplex vGA    = product(vT[external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    matrixComplex vsGR   = product(v_sigma_T[axis][spin], GR);
//                    matrixComplex vsGA   = product(v_sigma_T[axis][spin], GA);

//                    res[external][axis][spin] = tr(product(vsGR, vGA)) - 5e-1*(tr(product(vsGR, vGR)) + tr(product(vsGA, vGA)));
                    res[external][axis][spin] = tr(product(vsGR, vGA));
                }
            }
        }
        return res;
    };

    integrate_band_T(fn, response, b, mu);

    double coef = - charge * v0*v0 * hbar / (2e0*pi) / (angstrom*angstrom*angstrom);
    response = times(response, coef);

    return response;
}; // }}}

SHC get_angular_SHC_T2(band b, Energy epsilon, chemical_potential mu) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int band_index, chemical_potential mu, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_T(mu + epsilon*zi, k);
        Green_function GA = get_green_function_T(mu - epsilon*zi, k);

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
                    matrixComplex R = product(v_sigma_T[axis][spin], vR);
                    matrixComplex A = product(v_sigma_T[axis][spin], vA);

                    res[external][axis][spin] = tr(R) - tr(A);
                }
            }
        }
        return res;
    };

    integrate_band_T(fn, response, b, mu);

    double coef = - charge * v0*v0 * hbar / (4e0*pi) / (angstrom*angstrom*angstrom);
    response = times(response, coef);

    return response;
}; // }}}
// }}}

// L {{{
void set_response_L(chemical_potential ene_min, chemical_potential ene_max, int ene_mesh, int valley, int band_index) { // {{{
    std::string dir = "L"+std::to_string(valley+1)+"_"+std::to_string(bandsL)+"bands/band_index"+std::to_string(band_index);
//    std::string dir = "L"+std::to_string(bandsL)+"bands-isotropic/band_index"+std::to_string(band_index);
    set_output_directory(dir);

    for(int i=5; i<6; i++) {
        double epsilon = double(i)*1e-4;
        std::cout << std::scientific << "epsilon = " << epsilon << std::endl;

        chemical_potential d_ene = (ene_max - ene_min) / double(ene_mesh-1);

        band b;
        band b_main = set_band_L(valley, band_index, ene_min, ene_max, ene_mesh);
        if (mu_cutoff_L < ene_min) {
            band b_cut = set_band_L(valley, band_index, mu_cutoff_L, ene_min-d_ene, mu_cutoff_mesh_L);
            b = combine_band(b_cut, b_main);
        } else if (ene_max < mu_cutoff_L) {
            band b_cut = set_band_L(valley, band_index, ene_max+d_ene, mu_cutoff_L, mu_cutoff_mesh_L);
            b = combine_band(b_main, b_cut);
        } else {
            std::cerr << "mu_cutoff_L shoudl be smaller than ene_min or larger than ene_max" << std::endl;
            exit(0);
        }
        if ( (ene_min < band_edge_L[valley][band_index]) & (band_edge_L[valley][band_index] < ene_max) ) {
            int e_mesh = 15;
            double e_cut = 1.9e0*epsilon;
            double power = 7e-1;
            band b_edge = set_band_2n_L(valley, band_index, band_edge_L[valley][band_index], e_cut, e_mesh, power);
            b_main = combine_band_2n(b_main, b_edge);

        }

        std::vector<Conductivity> sigma_e;
        std::vector<SHC> sigma_m1, sigma_m2;
        std::vector<SHC> sigma_a1, sigma_a2;
        std::vector<SHC> sigma_r1, sigma_r2;
        sigma_e.resize(b_main.mesh);
        sigma_m1.resize(b_main.mesh);
        sigma_m2.resize(b_main.mesh);
        sigma_a1.resize(b_main.mesh);
        sigma_a2.resize(b_main.mesh);
        sigma_r1.resize(b_main.mesh);
        sigma_r2.resize(b_main.mesh);

        std::string filename;
// init dos file {{{
        filename = "dat/"+dir+"/dos_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofdos(filename);
// }}}
// init conductivity file {{{
        filename = "dat/"+dir+"/sigma_e_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofsigma1(filename);
        ofsigma1 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                ofsigma1 << ", " << axises[external]+axises[axis];
            }
        }
        ofsigma1 << std::endl;

        filename = "dat/"+dir+"/sigma_e_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofsigma2(filename);
        ofsigma2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                ofsigma2 << ", " << axises[external]+axises[axis];
            }
        }
        ofsigma2 << std::endl;
// }}}
// init spin magnetic Hall conductivity 1 file {{{
        filename = "dat/"+dir+"/sigma_m1_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc11(filename);
        ofshc11 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc11 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc11 << std::endl;

        filename = "dat/"+dir+"/sigma_m1_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc12(filename);
        ofshc12 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc12 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc12 << std::endl;
// }}}
// init spin magnetic Hall conductivity 2 file {{{
        filename = "dat/"+dir+"/sigma_m2_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc21(filename);
        ofshc21 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc21 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc21 << std::endl;

        filename = "dat/"+dir+"/sigma_m2_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc22(filename);
        ofshc22 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc22 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc22 << std::endl;
// }}}
// init spin angular Hall conductivity 1 file {{{
        filename = "dat/"+dir+"/sigma_a1_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc11(filename);
        of_angular_shc11 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc11 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc11 << std::endl;

        filename = "dat/"+dir+"/sigma_a1_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc12(filename);
        of_angular_shc12 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc12 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc12 << std::endl;
// }}}
// init spin angular Hall conductivity 2 file {{{
        filename = "dat/"+dir+"/sigma_a2_real_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc21(filename);
        of_angular_shc21 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc21 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc21 << std::endl;

        filename = "dat/"+dir+"/sigma_a2_imag_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc22(filename);
        of_angular_shc22 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc22 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc22 << std::endl;
// }}}

        int i_ene;
        for(i_ene=0; i_ene<b_main.mesh; i_ene++) {
            const chemical_potential mu = b_main.ene[i_ene];
            std::cout << "start: i_ene = " << i_ene << ": mu = " << mu << std::endl;

// dos output {{{
            ofdos << std::scientific << mu << ", " << b_main.dos[i_ene] / (angstrom*angstrom*angstrom) << std::endl;
// }}}

            int e_mesh = 47;
            double e_cut = 59e0*epsilon;
            double power = 9e-1;
            band bL = set_band_2n_L(valley, band_index, mu, e_cut, e_mesh, power);

            band b_sum = combine_band_2n(b, bL);

            sigma_e[i_ene] = get_conductivity_L(b_sum, epsilon, mu, valley);
// conductivity output {{{
            ofsigma1 << std::scientific << mu;
            ofsigma2 << std::scientific << mu;
            for( auto s : sigma_e[i_ene] ) {
                for( auto v : s ) {
                    ofsigma1 << std::scientific << ", " << v.real();
                    ofsigma2 << std::scientific << ", " << v.imag();
                }
            }
            ofsigma1 << std::endl;
            ofsigma2 << std::endl;
// }}}

            sigma_m1[i_ene] = get_SHC_L1(b_sum, epsilon, mu, valley);
// magnetic SHC1 output {{{
            ofshc11 << std::scientific << mu;
            ofshc12 << std::scientific << mu;
            for( auto e : sigma_m1[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc11 << std::scientific << ", " << s.real();
                        ofshc12 << std::scientific << ", " << s.imag();
                    }
                }
            }
            ofshc11 << std::endl;
            ofshc12 << std::endl;
// }}}

            sigma_m2[i_ene] = get_SHC_L2(b_sum, epsilon, mu, valley);
// magnetic SHC2 output {{{
            ofshc21 << std::scientific << mu;
            ofshc22 << std::scientific << mu;
            for( auto e : sigma_m2[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc21 << std::scientific << ", " << s.real();
                        ofshc22 << std::scientific << ", " << s.imag();
                    }
                }
            }
            ofshc21 << std::endl;
            ofshc22 << std::endl;
// }}}

            sigma_a1[i_ene] = get_angular_SHC_L1(b_sum, epsilon, mu, valley);
// angular SHC1 output {{{
            of_angular_shc11 << std::scientific << mu;
            of_angular_shc12 << std::scientific << mu;
            for( auto e : sigma_a1[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc11 << std::scientific << ", " << s.real();
                        of_angular_shc12 << std::scientific << ", " << s.imag();
                    }
                }
            }
            of_angular_shc11 << std::endl;
            of_angular_shc12 << std::endl;
// }}}

            sigma_a2[i_ene] = get_angular_SHC_L2(b_sum, epsilon, mu, valley);
// angular SHC2 output {{{
            of_angular_shc21 << std::scientific << mu;
            of_angular_shc22 << std::scientific << mu;
            for( auto e : sigma_a2[i_ene] ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc21 << std::scientific << ", " << s.real();
                        of_angular_shc22 << std::scientific << ", " << s.imag();
                    }
                }
            }
            of_angular_shc21 << std::endl;
            of_angular_shc22 << std::endl;
// }}}
        }

        Conductivity Sigma(space_dim, vectorComplex(space_dim, 0e0));
        SHC magnetic_SHC1(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        SHC magnetic_SHC2(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        SHC angular_SHC1(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        SHC angular_SHC2(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

// mu dependence @ T = 0 {{{
        dir = dir+"/mu-dependence";
        set_output_directory(dir);
// init Sigma file {{{
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
// init magnetic SHC1 file {{{
        filename = "dat/"+dir+"/spin-magnetic-conductivity1_eps"+std::to_string(epsilon)+".csv";
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
// init magnetic SHC2 file {{{
        filename = "dat/"+dir+"/spin-magnetic-conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream ofshc2(filename);
        ofshc2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    ofshc2 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        ofshc2 << std::endl;
// }}}
// init angular SHC1 file {{{
        filename = "dat/"+dir+"/spin-angular-conductivity1_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc1(filename);
        of_angular_shc1 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc1 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc1 << std::endl;
// }}}
// init angular SHC2 file {{{
        filename = "dat/"+dir+"/spin-angular-conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_angular_shc2(filename);
        of_angular_shc2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_angular_shc2 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_angular_shc2 << std::endl;
// }}}
// init real SHC1 file {{{
        filename = "dat/"+dir+"/real-spin-angular-conductivity1_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_real_shc1(filename);
        of_real_shc1 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_real_shc1 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_real_shc1 << std::endl;
// }}}
// init real SHC2 file {{{
        filename = "dat/"+dir+"/real-spin-angular-conductivity2_eps"+std::to_string(epsilon)+".csv";
        std::ofstream of_real_shc2(filename);
        of_real_shc2 << "mu";
        for(int external=0; external<space_dim; external++) {
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    of_real_shc2 << ", " << axises[external]+axises[axis]+axises[spin];
                }
            }
        }
        of_real_shc2 << std::endl;
// }}}

        chemical_potential mu;
        for(int i_ene=0; i_ene<b_main.mesh; i_ene++) {
            mu = b_main.ene[i_ene];
            Sigma = sigma_e[i_ene];
// Sigma output {{{
            ofsigma << std::scientific << mu;
            for( auto s : Sigma ) {
                for( auto v : s ) {
                    ofsigma << std::scientific << ", " << v.real();
                }
            }
            ofsigma << std::endl;
// }}}
            magnetic_SHC1 = sigma_m1[i_ene];
// magnetic SHC1 output {{{
            ofshc1 << std::scientific << mu;
            for( auto e : magnetic_SHC1 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc1 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc1 << std::endl;
// }}}
            angular_SHC1 = sigma_a1[i_ene];
// angular SHC1 output {{{
            of_angular_shc1 << std::scientific << mu;
            for( auto e : angular_SHC1 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc1 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc1 << std::endl;
// }}}
        }

        double de;
        SHC magnetic_SHC2_mu;
        i_ene = 0;
            mu = b_main.ene[i_ene];
// magnetic SHC2 output {{{
            ofshc2 << std::scientific << mu;
            for( auto e : magnetic_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc2 << std::endl;
// }}}
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene])*5e-1;
            magnetic_SHC2_mu = times(sigma_m2[i_ene], de);
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
        for(int i_ene=1; i_ene<b_main.mesh-1; i_ene++) {
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene-1])*2.5e-1;
            magnetic_SHC2_mu = times(sigma_m2[i_ene], de);
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
// magnetic SHC2 output {{{
            ofshc2 << std::scientific << mu;
            for( auto e : magnetic_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc2 << std::endl;
// }}}
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
        }
        i_ene = b_main.mesh-1;
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene] - b_main.ene[i_ene-1])*5e-1;
            magnetic_SHC2_mu = times(sigma_m2[i_ene], de);
            magnetic_SHC2 = add(magnetic_SHC2, magnetic_SHC2_mu);
// magnetic SHC2 output {{{
            ofshc2 << std::scientific << mu;
            for( auto e : magnetic_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        ofshc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            ofshc2 << std::endl;
// }}}

        SHC angular_SHC2_mu;
        i_ene = 0;
            mu = b_main.ene[i_ene];
// angular SHC2 output {{{
            of_angular_shc2 << std::scientific << mu;
            for( auto e : angular_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc2 << std::endl;
// }}}
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene])*5e-1;
            angular_SHC2_mu = times(sigma_a2[i_ene], de);
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
        for(int i_ene=1; i_ene<b_main.mesh-1; i_ene++) {
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene+1] - b_main.ene[i_ene-1])*2.5e-1;
            angular_SHC2_mu = times(sigma_a2[i_ene], de);
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
// angular SHC2 output {{{
            of_angular_shc2 << std::scientific << mu;
            for( auto e : angular_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc2 << std::endl;
// }}}
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
        }
        i_ene = b_main.mesh-1;
            mu = b_main.ene[i_ene];
            de = (b_main.ene[i_ene] - b_main.ene[i_ene-1])*5e-1;
            angular_SHC2_mu = times(sigma_a2[i_ene], de);
            angular_SHC2 = add(angular_SHC2, angular_SHC2_mu);
// angular SHC2 output {{{
            of_angular_shc2 << std::scientific << mu;
            for( auto e : angular_SHC2 ) {
                for( auto a : e ) {
                    for( auto s : a ) {
                        of_angular_shc2 << std::scientific << ", " << s.real();
                    }
                }
            }
            of_angular_shc2 << std::endl;
// }}}
    }
// }}}
}; // }}}

Conductivity get_conductivity_L(band b, Energy epsilon, chemical_potential mu, int valley) { // {{{
    Conductivity response(space_dim, vectorComplex(space_dim, 0e0));

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_L(e + epsilon*zi, valley, k);
        Green_function GA = get_green_function_L(e - epsilon*zi, valley, k);

        matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));

        for(int external=0; external<space_dim; external++) {
            matrixComplex vGRe = product(vL[valley][external], GR);
            matrixComplex vGAe = product(vL[valley][external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                matrixComplex vGRa = product(vL[valley][axis], GR);
                matrixComplex vGAa = product(vL[valley][axis], GA);

                res[external][axis] = tr(product(vGRa, vGAe)) - 5e-1*(tr(product(vGRa, vGRe)) + tr(product(vGAa, vGAe)));
            }
        }

        return res;
    };

    integrate_band_L(fn, response, b, valley, mu);

    double coef = (charge*v0)*(charge*v0) * hbar / (4e0*pi) / (angstrom*angstrom*angstrom) / charge;
    response = times(response, coef);

    return response;
}; // }}}

void set_conductivity_damping_dependence_at_Fermi_level_L(int valley, int band_index) { // {{{
    vectorReal damping = {1e-5, 2e-5, 3e-5, 5e-5, 7e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4};

    std::vector<Conductivity> sigma(damping.size());
    for (int i=0; i<damping.size(); i++) {
        double damping_constant = damping[i];
        Energy epsilon = damping_constant;
        chemical_potential mu = 0e0;
        int e_mesh = 47;
        double e_cut = 59e0*epsilon;
        double power = 9e-1;

        band bL = set_band_2n_L(valley, band_index, mu, e_cut, e_mesh, power);

        nu_F_L[valley] = bL.dos[e_mesh];
        double coef = damping_constant / nu_F_L[valley];
        Self_energy se = get_self_energy_born_L(bL, 0e0, valley, mu, epsilon, coef);
        se = add(product(impurityV1_L[valley], product(se, impurityV1_L[valley])), product(impurityV2_L[valley], product(se, impurityV2_L[valley])));
        sigma[i] = get_conductivity_with_self_energy_L(bL, 0e0, mu, valley, se);
    }
// init Sigma file {{{
    std::string dir = "L"+std::to_string(bandsL)+"bands/damping-dependence";
    set_output_directory(dir);
    std::string filename = "dat/"+dir+"/conductivity_L"+std::to_string(valley+1)+".csv";
    std::ofstream ofsigma(filename);
    ofsigma << "damping";
    for(int external=0; external<space_dim; external++) {
        for(int axis=0; axis<space_dim; axis++) {
            ofsigma << ", " << axises[external]+axises[axis];
        }
    }
    ofsigma << std::endl;

    for (int i=0; i<damping.size(); i++) {
        ofsigma << std::scientific << damping[i];
        Conductivity s = sigma[i];
        for( auto v : s ) {
            for( auto r : v ) {
                ofsigma << std::scientific << ", " << r.real();
            }
        }
        ofsigma << std::endl;
    }
// }}}
}; // }}}

Conductivity get_conductivity_with_self_energy_L(band b, Energy epsilon, chemical_potential mu, int valley, Self_energy se) { // {{{
    Conductivity response(space_dim, vectorComplex(space_dim, 0e0));

    Self_energy SigmaR = se;
    Self_energy SigmaA(bandsL, vectorComplex(bandsL));
    for(int i=0; i<bandsL; i++) {
        SigmaA[i][i] = std::conj(SigmaR[i][i]);
        for(int j=i+1; j<bandsL; j++) {
            SigmaA[i][j] = std::conj(SigmaR[i][j]);
            SigmaA[j][i] = std::conj(SigmaR[j][i]);
        }
    }

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_full_green_function_L(e + eps_phys*zi, valley, k, SigmaR);
        Green_function GA = get_full_green_function_L(e - eps_phys*zi, valley, k, SigmaA);

        matrixComplex res(space_dim, vectorComplex(space_dim, 0e0));

        for(int external=0; external<space_dim; external++) {
            matrixComplex vGRe = product(vL[valley][external], GR);
            matrixComplex vGAe = product(vL[valley][external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                matrixComplex vGRa = product(vL[valley][axis], GR);
                matrixComplex vGAa = product(vL[valley][axis], GA);

                res[external][axis] = tr(product(vGRa, vGAe)) - 5e-1*(tr(product(vGRa, vGRe)) + tr(product(vGAa, vGAe)));
            }
        }

        return res;
    };

    integrate_band_L(fn, response, b, valley, mu);

    double coef = (charge*v0)*(charge*v0) * hbar / (4e0*pi) / (angstrom*angstrom*angstrom) / charge;
    response = times(response, coef);

    return response;
}; // }}}

SHC get_SHC_L1(band b, Energy epsilon, chemical_potential mu, int valley) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_L(e + epsilon*zi, valley, k);
        Green_function GA = get_green_function_L(e - epsilon*zi, valley, k);
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

    integrate_band_L(fn, response, b, valley, mu);

    double coef = - charge * v0*v0 * hbar / (2e0*pi) / (angstrom*angstrom*angstrom) / muBeV;
    response = times(response, coef);

    return response;
}; // }}}

SHC get_SHC_L2(band b, Energy epsilon, chemical_potential mu, int valley) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int valley, int band_index, chemical_potential mu, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_L(mu + epsilon*zi, valley, k);
        Green_function GA = get_green_function_L(mu - epsilon*zi, valley, k);

        for(int external=0; external<space_dim; external++) {
            matrixComplex vGR = product(vL[valley][external], GR);
            matrixComplex vGA = product(vL[valley][external], GA);
            matrixComplex GRv = product(GR, vL[valley][external]);
            matrixComplex GAv = product(GA, vL[valley][external]);
            matrixComplex QR  = minus(GRv, vGR);
            matrixComplex QA  = minus(GAv, vGA);
            matrixComplex vR  = product(GR, product(QR, GR));
            matrixComplex vA  = product(GA, product(QA, GA));
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    matrixComplex R = product(v_s_L[valley][axis][spin], vR);
                    matrixComplex A = product(v_s_L[valley][axis][spin], vA);

                    res[external][axis][spin] = tr(R) - tr(A);
                }
            }
        }
        return res;
    };

    integrate_band_L(fn, response, b, valley, mu);

    double coef = - charge * v0*v0 * hbar / (4e0*pi) / (angstrom*angstrom*angstrom) / muBeV;
    response = times(response, coef);

    return response;
}; // }}}

SHC get_angular_SHC_L1(band b, Energy epsilon, chemical_potential mu, int valley) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_L(e + epsilon*zi, valley, k);
        Green_function GA = get_green_function_L(e - epsilon*zi, valley, k);
        for(int external=0; external<space_dim; external++) {
//            matrixComplex vGR    = product(vL[valley][external], GR);
            matrixComplex vGA    = product(vL[valley][external], GA);
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    matrixComplex vsGR   = product(v_sigma_L[valley][axis][spin], GR);
//                    matrixComplex vsGA   = product(v_sigma_L[valley][axis][spin], GA);

//                    res[external][axis][spin] = tr(product(vsGR, vGA)) - 5e-1*(tr(product(vsGR, vGR)) + tr(product(vsGA, vGA)));
                    res[external][axis][spin] = tr(product(vsGR, vGA));
                }
            }
        }
        return res;
    };

    integrate_band_L(fn, response, b, valley, mu);

    double coef = - charge * v0*v0 * hbar / (2e0*pi) / (angstrom*angstrom*angstrom);
    response = times(response, coef);

    return response;
}; // }}}

SHC get_angular_SHC_L2(band b, Energy epsilon, chemical_potential mu, int valley) { // {{{
    SHC response(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));

    auto fn = [=](int valley, int band_index, chemical_potential mu, kpoint k) {
        SHC res(space_dim, matrixComplex(space_dim, vectorComplex(spin_dim, 0e0)));
        Green_function GR = get_green_function_L(mu + epsilon*zi, valley, k);
        Green_function GA = get_green_function_L(mu - epsilon*zi, valley, k);

        for(int external=0; external<space_dim; external++) {
            matrixComplex vGR = product(vL[valley][external], GR);
            matrixComplex vGA = product(vL[valley][external], GA);
            matrixComplex GRv = product(GR, vL[valley][external]);
            matrixComplex GAv = product(GA, vL[valley][external]);
            matrixComplex QR  = minus(GRv, vGR);
            matrixComplex QA  = minus(GAv, vGA);
            matrixComplex vR  = product(GR, product(QR, GR));
            matrixComplex vA  = product(GA, product(QA, GA));
            for(int axis=0; axis<space_dim; axis++) {
                for(int spin=0; spin<spin_dim; spin++) {
                    matrixComplex R = product(v_sigma_L[valley][axis][spin], vR);
                    matrixComplex A = product(v_sigma_L[valley][axis][spin], vA);

                    res[external][axis][spin] = tr(R) - tr(A);
                }
            }
        }
        return res;
    };

    integrate_band_L(fn, response, b, valley, mu);

    double coef = - charge * v0*v0 * hbar / (4e0*pi) / (angstrom*angstrom*angstrom);
    response = times(response, coef);

    return response;
}; // }}}
// }}}

void write_res(Conductivity sigma, chemical_potential mu, std::string filename) { // {{{
     std::ofstream ofs(filename, std::ios::app);
     ofs << std::scientific << mu;
     for( auto s : sigma ) {
         for( auto v : s ) {
             ofs << std::scientific << ", " << v.real();
         }
     }
     ofs << std::endl;

}; // }}}

void write_res(SHC sigma, chemical_potential mu, std::string filename) { // {{{
     std::ofstream ofs(filename, std::ios::app);
     ofs << std::scientific << mu;
     for( auto e : sigma ) {
         for( auto a : e ) {
             for( auto s : a ) {
                 ofs << std::scientific << ", " << s.real();
             }
         }
     }
     ofs << std::endl;
}; // }}}
