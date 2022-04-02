#include "parameters.hpp"

int thread_num = std::thread::hardware_concurrency();

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
    b.dos.resize(mu_mesh);
    mtx.lock();
    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search start" << std::endl;
    mtx.unlock();
    std::string dos = "dat/L_dos"+std::to_string(band_index)+".csv";
    std::ofstream ofs(dos);
    for(int i_mu=0; i_mu<mu_mesh; i_mu++) {
        chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
        std::cout << "mu = " << mu << std::endl;
        b.tri[i_mu] = get_triangles_L(valley, band_index, mu);
        b.dos[i_mu] = get_DOS_L(b.tri[i_mu], valley, band_index, mu);
        ofs << std::scientific << b.tri[i_mu].ene << ", " << b.dos[i_mu] << std::endl;
    }
    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search end" << std::endl;
}; // }}}

band set_band_L(int valley, int band_index, chemical_potential mu_min, chemical_potential mu_max, int mu_mesh) { // {{{
    band b;
    b.index   = band_index;
    b.mu_min  = mu_min;
    b.mu_max  = mu_max;
    b.dmu     = (mu_max - mu_min) / double(mu_mesh-1);
    b.mu_mesh = mu_mesh;
    b.tri.resize(mu_mesh);
    b.dos.resize(mu_mesh);
    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search start" << std::endl;
    std::string dos = "dat/L_dos"+std::to_string(band_index)+".csv";
    std::ofstream ofs(dos);

    std::vector<std::thread> threads;
    threads.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        auto func =[](int i_thread, band& b, int valley, int band_index) {
            for(int i_mu=i_thread; i_mu<b.mu_mesh; i_mu=i_mu+thread_num) {
                chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
                mtx.lock();
                std::cout << "mu = " << mu << std::endl;
                mtx.unlock();
                b.tri[i_mu] = get_triangles_L(valley, band_index, mu);
                b.dos[i_mu] = get_DOS_L(b.tri[i_mu], valley, band_index, mu);
            }
        };
        threads[i_thread] = std::thread(func, i_thread, std::ref(b), valley, band_index);
    }
    for(auto& thread : threads){
        thread.join();
    }

    for(int i_mu=0; i_mu<b.mu_mesh; i_mu++) {
        chemical_potential mu = b.mu_min + b.dmu*double(i_mu);
        ofs << std::scientific << b.tri[i_mu].ene << ", " << b.dos[i_mu] << std::endl;
    }

    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search end" << std::endl;

    return b;
}; // }}}
