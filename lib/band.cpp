#include "parameters.hpp"

void get_band_T(band& b, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh) { // {{{
    b.index   = band_index;
    b.mu_min  = mu_min;
    b.mu_max  = mu_max;
    b.dmu     = (mu_max - mu_min) / double(mu_mesh-1);
    b.mu_mesh = mu_mesh;
    b.tri.resize(mu_mesh);
    b.dos.resize(mu_mesh);
    mtx.lock();
    std::cout << "T point; band#" << band_index << " search start" << std::endl;
    mtx.unlock();
    for(int i_mu=0; i_mu<mu_mesh; i_mu++) {
        chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
        b.tri[i_mu] = get_triangles_T(band_index, mu);
        b.dos[i_mu] = get_DOS_T(b.tri[i_mu], band_index, mu);
    }
    std::cout << "T point; band#" << band_index << " search end" << std::endl;
}; // }}}

void get_band_L(band& b, int valley, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh) { // {{{
    b.index   = band_index;
    b.mu_min  = mu_min;
    b.mu_max  = mu_max;
    b.dmu     = (mu_max - mu_min) / double(mu_mesh-1);
    b.mu_mesh = mu_mesh;
    b.tri.resize(mu_mesh);
    mtx.lock();
    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search start" << std::endl;
    mtx.unlock();
    for(int i_mu=0; i_mu<mu_mesh; i_mu++) {
        chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
        b.tri[i_mu] = get_triangles_L(valley, band_index, mu);
    }
    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search end" << std::endl;
}; // }}}
