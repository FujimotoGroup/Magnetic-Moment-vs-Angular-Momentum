#include "parameters.hpp"

void get_band_T(band& b, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh) { // {{{
    b.index   = band_index;
    b.mu_min  = mu_min;
    b.mu_max  = mu_max;
    b.dmu     = (mu_max - mu_min) / double(mu_mesh-1);
    b.mu_mesh = mu_mesh;
    b.fs.resize(mu_mesh);
    for(int i_mu=0; i_mu<mu_mesh; i_mu++) {
        chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
        b.fs[i_mu] = get_fermi_surace_T(band_index, mu);
    }
}; // }}}

void get_band_L(band& b, int valley, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh) { // {{{
    b.index   = band_index;
    b.mu_min  = mu_min;
    b.mu_max  = mu_max;
    b.dmu     = (mu_max - mu_min) / double(mu_mesh-1);
    b.mu_mesh = mu_mesh;
    for(int i_mu=0; i_mu<mu_mesh; i_mu++) {
        chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
        b.fs.push_back(get_fermi_surace_L(valley, band_index, mu));
    }
}; // }}}
