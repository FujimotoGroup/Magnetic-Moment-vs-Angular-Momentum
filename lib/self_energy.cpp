#include "parameters.hpp"

Self_energy get_self_energy_born_T(band b, Energy ene, chemical_potential mu, Energy epsilon) { // {{{
    Self_energy self_ene(bandsT, vectorComplex(bandsT, 0e0));

    auto fn = [=](int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_T(e + epsilon*zi, k);

        Self_energy se(bandsT, vectorComplex(bandsT, 0e0));

        for(int i=0; i<bandsT; i++) {
            se[i][i] = zi*GR[i][i].imag();
            for(int j=i+1; j<bandsT; j++) {
                se[i][j] = zi * GR[i][j].imag();
                se[j][i] = zi * GR[j][i].imag();
            }
        }

        return se;
    };

    integrate_band_T(fn, self_ene, b, ene + mu);

    return self_ene;

} // }}}

Self_energy get_self_energy_born_L(band b, Energy ene, int valley, chemical_potential mu, Energy epsilon) { // {{{
    Self_energy self_ene(bandsL, vectorComplex(bandsL, 0e0));

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_L(e + epsilon*zi, valley, k);

        Self_energy se(bandsL, vectorComplex(bandsL, 0e0));

        for(int i=0; i<bandsL; i++) {
            se[i][i] = zi * GR[i][i].imag();
            for(int j=i+1; j<bandsL; j++) {
                se[i][j] = zi * GR[i][j].imag();
                se[j][i] = zi * GR[j][i].imag();
            }
//            se[i][i] = GR[i][i];
//            for(int j=i+1; j<bandsL; j++) {
//                se[i][j] = GR[i][j];
//                se[j][i] = GR[j][i];
//            }
        }

        return se;
    };

    integrate_band_L(fn, self_ene, b, valley, ene + mu);

    return self_ene;

} // }}}

vectorReal get_lifetime_T(int band_index, Self_energy se, triangles tri) { // {{{
    int size = tri.vertexes.size();

    vectorReal lifetime(size, 0e0);

    for (int i=0; i<size; i++) {
        kpoint k = tri.vertexes[i];
        matrixComplex H = set_T(k.vec);
        diag_set eigen = diagonalize_V(H);

        vectorComplex U = eigen.vectors[band_index];
        Complex tau = 0e0;
        for (int j=0; j<U.size(); j++) {
            for (int l=0; l<U.size(); l++) {
                tau = tau + std::conj(U[j])*se[j][l]*U[l];
            }
        }

        lifetime[i] = - hbar / (2e0*tau.imag());
    }

    return lifetime;
} // }}}

vectorReal get_lifetime_L(int band_index, Self_energy se, int valley, triangles tri) { // {{{
    int size = tri.vertexes.size();

    vectorReal lifetime(size, 0e0);

    for (int i=0; i<size; i++) {
        kpoint k = tri.vertexes[i];
        matrixComplex H = set_L(valley, k.vec);
        diag_set eigen = diagonalize_V(H);

        vectorComplex U = eigen.vectors[band_index];
        Complex tau = 0e0;
        for (int j=0; j<U.size(); j++) {
            for (int l=0; l<U.size(); l++) {
                tau += std::conj(U[j])*se[j][l]*U[l];
            }
        }

        lifetime[i] = - hbar / (2e0*tau.imag());
    }

    return lifetime;
} // }}}

Self_energy get_self_energy_born_k_T(band b, Energy ene, kpoint kp, chemical_potential mu, Energy epsilon) { // {{{
    Self_energy self_ene(bandsT, vectorComplex(bandsT, 0e0));

    auto fn = [=](int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_T(e + epsilon*zi, k);

        Self_energy se(bandsT, vectorComplex(bandsT, 0e0));

        double coeff = 0e0;
        for (int axis=0; axis<space_dim; axis++) {
            double c = kp.vec[axis]-k.vec[axis];
            coeff += c*c;
        }
        coeff = std::exp(- coeff * sigma_imp*sigma_imp);

        for(int i=0; i<bandsL; i++) {
            se[i][i] = zi*coeff*GR[i][i].imag();
            for(int j=i+1; j<bandsL; j++) {
                se[i][j] = zi * coeff * GR[i][j].imag();
                se[j][i] = zi * coeff * GR[j][i].imag();
            }
        }

        return se;
    };

    integrate_band_T(fn, self_ene, b, ene + mu);

    return self_ene;

} // }}}

Self_energy get_self_energy_born_k_L(band b, Energy ene, int valley, kpoint kp, chemical_potential mu, Energy epsilon) { // {{{
    Self_energy self_ene(bandsL, vectorComplex(bandsL, 0e0));

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_L(e + epsilon*zi, valley, k);

        Self_energy se(bandsL, vectorComplex(bandsL, 0e0));

        double coeff = 0e0;
        for (int axis=0; axis<space_dim; axis++) {
            double c = kp.vec[axis]-k.vec[axis];
            coeff += c*c;
        }
        coeff = std::exp(- coeff * sigma_imp*sigma_imp);

        for(int i=0; i<bandsL; i++) {
            se[i][i] = zi*coeff*GR[i][i].imag();
            for(int j=i+1; j<bandsL; j++) {
                se[i][j] = zi * coeff * GR[i][j].imag();
                se[j][i] = zi * coeff * GR[j][i].imag();
            }
        }

        return se;
    };

    integrate_band_L(fn, self_ene, b, valley, ene + mu);

    return self_ene;

} // }}}

vectorReal get_lifetime_Gaussian_T(band b, chemical_potential mu, triangles tri, Energy epsilon) { // {{{ int size = tri.vertexes.size();
    int size = tri.vertexes.size();

    vectorReal lifetime(size, 0e0);

    for (int i=0; i<size; i++) {
        kpoint k = tri.vertexes[i];
        matrixComplex H = set_T(k.vec);
        diag_set eigen = diagonalize_V(H);

        Self_energy se = get_self_energy_born_k_T(b, 0e0, k, mu, epsilon);

        vectorComplex U = eigen.vectors[b.index];
        Complex gamma = 0e0;
        for (int j=0; j<U.size(); j++) {
            for (int l=0; l<U.size(); l++) {
                gamma += std::conj(U[j])*se[j][l]*U[l];
            }
        }

        lifetime[i] = - hbar / (2e0*gamma.imag());
    }

    return lifetime;
} // }}}

vectorReal get_lifetime_Gaussian_L(band b, int valley, chemical_potential mu, triangles tri, Energy epsilon) { // {{{
    int size = tri.vertexes.size();

    vectorReal lifetime(size, 0e0);

    std::vector<std::vector<int>> indexes(thread_num);
    std::vector<std::thread> threads;
    threads.resize(thread_num);
    std::vector<vectorReal> part;
    part.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        for (int i=i_thread; i<size; i=i+thread_num) {
            indexes[i_thread].push_back(i);
        }
    }

    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        auto func =[&](band b, int i_thread, int valley) {
            for(int i=0; i<indexes[i_thread].size(); i++) {
                mtx.lock();
//                if (i%10 == 0) std::cout << i << std::endl;
                std::cout << i << std::endl;
                mtx.unlock();
                kpoint k = tri.vertexes[indexes[i_thread][i]];
                matrixComplex H = set_L(valley, k.vec);
                diag_set eigen = diagonalize_V(H);

                Self_energy se = get_self_energy_born_k_L(b, 0e0, valley, k, mu, epsilon);
                se = add(product(impurityV1_L[valley], product(se, impurityV1_L[valley])), product(impurityV2_L[valley], product(se, impurityV2_L[valley])));

                vectorComplex U = eigen.vectors[b.index];
                Complex gamma = 0e0;
                for (int j=0; j<U.size(); j++) {
                    for (int l=0; l<U.size(); l++) {
                        gamma += std::conj(U[j])*se[j][l]*U[l];
                    }
                }
//                double tau = - hbar / (2e0*gamma.imag());
                double tau = - gamma.imag();

                part[i_thread].push_back(tau);
            }
        };
        threads[i_thread] = std::thread(func, b, i_thread, valley);
    }

    for(auto& thread : threads){
        thread.join();
    }

    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        int j = 0;
        for (auto i : indexes[i_thread]) {
            lifetime[i] = part[i_thread][j];
            j++;
        }
    }

    return lifetime;
} // }}}

void write_res(Self_energy sigma, chemical_potential mu, std::string filename) { // {{{
     std::ofstream ofs(filename, std::ios::app);
     ofs << std::scientific << mu;
     for( auto s : sigma ) {
         for( auto v : s ) {
             ofs << std::scientific << ", " << v.imag();
         }
     }
     ofs << std::endl;

}; // }}}
