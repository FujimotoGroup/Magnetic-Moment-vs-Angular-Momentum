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
            se[i][i] = GR[i][i].imag();
            for(int j=i+1; j<bandsL; j++) {
                se[i][j] = GR[i][j].imag();
                se[j][i] = GR[j][i].imag();
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

        lifetime[i] = - tau.imag();
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
                tau = tau + std::conj(U[j])*se[j][l]*U[l];
            }
        }

        lifetime[i] = - tau.imag();
    }

    return lifetime;
} // }}}

Self_energy get_self_energy_born_k_T(band b, Energy ene, kpoint kp, chemical_potential mu, Energy epsilon) { // {{{
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

Self_energy get_self_energy_born_k_L(band b, Energy ene, int valley, chemical_potential mu, Energy epsilon) { // {{{
    Self_energy self_ene(bandsL, vectorComplex(bandsL, 0e0));

    auto fn = [=](int valley, int band_index, chemical_potential e, kpoint k) {
        Green_function GR = get_green_function_L(e + epsilon*zi, valley, k);

        Self_energy se(bandsL, vectorComplex(bandsL, 0e0));

        for(int i=0; i<bandsL; i++) {
            se[i][i] = zi*GR[i][i].imag();
            for(int j=i+1; j<bandsL; j++) {
                se[i][j] = zi * GR[i][j].imag();
                se[j][i] = zi * GR[j][i].imag();
            }
        }

        return se;
    };

    integrate_band_L(fn, self_ene, b, valley, ene + mu);

    return self_ene;

} // }}}

