#include "parameters.hpp"

band set_band_T(int band_index, Energy ene_min, Energy ene_max, int ene_mesh) { // {{{
    band b;
    b.index = band_index;
    b.mesh  = ene_mesh;
    b.ene.resize(b.mesh);
    b.tri.resize(b.mesh);
    b.dos.resize(b.mesh);
    double dmu = (ene_max - ene_min) / double(b.mesh-1);
    for(int i=0; i<b.mesh; i++) {
        b.ene[i] = ene_min + dmu*double(i);
    }

    std::cout << "T point; band#" << band_index << " search start" << std::endl;

    std::vector<std::thread> threads;
    threads.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        auto func =[](int i_thread, band& b, int band_index) {
            for(int i_mu=i_thread; i_mu<b.mesh; i_mu=i_mu+thread_num) {
                mtx.lock();
                std::cout << "mu = " << b.ene[i_mu] << std::endl;
                mtx.unlock();
                b.tri[i_mu] = get_triangles_T(band_index, b.ene[i_mu]);
                b.dos[i_mu] = get_DOS_T(b.tri[i_mu], band_index, b.ene[i_mu]);
            }
        };
        threads[i_thread] = std::thread(func, i_thread, std::ref(b), band_index);
    }
    for(auto& thread : threads){
        thread.join();
    }

    std::cout << "T point; band#" << band_index << " search end" << std::endl;

    return b;
}; // }}}

band set_band_2n_T(int band_index, Energy ene_center, Energy delta, int n) { // {{{
    band b;
    b.index = band_index;
    b.mesh  = 2*n+1;
    b.ene.resize(b.mesh);
    b.tri.resize(b.mesh);
    b.dos.resize(b.mesh);

    b.ene[n] = ene_center;
    for(int i=0; i<n; i++) {
        double dn = delta;
        for(int j = 0; j<i; j++) dn *= 5e-1;
        b.ene[i] = ene_center - dn;

        dn = delta;
        for(int j = 0; j<n-i-1; j++) dn *= 5e-1;
        b.ene[n+i+1] = ene_center + dn;
    }

    std::cout << "T point; band#" << band_index << " search start" << std::endl;

    std::vector<std::thread> threads;
    threads.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        auto func =[](int i_thread, band& b, int band_index) {
            for(int i_mu=i_thread; i_mu<b.mesh; i_mu=i_mu+thread_num) {
                mtx.lock();
                std::cout << "mu = " << b.ene[i_mu] << std::endl;
                mtx.unlock();
                b.tri[i_mu] = get_triangles_T(band_index, b.ene[i_mu]);
                b.dos[i_mu] = get_DOS_T(b.tri[i_mu], band_index, b.ene[i_mu]);
            }
        };
        threads[i_thread] = std::thread(func, i_thread, std::ref(b), band_index);
    }
    for(auto& thread : threads){
        thread.join();
    }

    std::cout << "T point; band#" << band_index << " search end" << std::endl;

    return b;
}; // }}}

band set_band_L(int valley, int band_index, Energy ene_min, Energy ene_max, int ene_mesh) { // {{{
    band b;
    b.index = band_index;
    b.mesh  = ene_mesh;
    b.ene.resize(b.mesh);
    b.tri.resize(b.mesh);
    b.dos.resize(b.mesh);
    double dmu = (ene_max - ene_min) / double(b.mesh-1);
    for(int i=0; i<b.mesh; i++) {
        b.ene[i] = ene_min + dmu*double(i);
    }

    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search start" << std::endl;

    std::vector<std::thread> threads;
    threads.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        auto func =[](int i_thread, band& b, int valley, int band_index) {
            for(int i_mu=i_thread; i_mu<b.mesh; i_mu=i_mu+thread_num) {
                mtx.lock();
                std::cout << "mu = " << b.ene[i_mu] << std::endl;
                mtx.unlock();
                b.tri[i_mu] = get_triangles_L(valley, band_index, b.ene[i_mu]);
                b.dos[i_mu] = get_DOS_L(b.tri[i_mu], valley, band_index, b.ene[i_mu]);
            }
        };
        threads[i_thread] = std::thread(func, i_thread, std::ref(b), valley, band_index);
    }
    for(auto& thread : threads){
        thread.join();
    }

    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search end" << std::endl;

    return b;
}; // }}}

band set_band_2n_L(int valley, int band_index, Energy ene_center, Energy delta, int n) { // {{{
    band b;
    b.index = band_index;
    b.mesh  = 2*n+1;
    b.ene.resize(b.mesh);
    b.tri.resize(b.mesh);
    b.dos.resize(b.mesh);

    double p = 7e-1;
    b.ene[n] = ene_center;
    for(int i=0; i<n; i++) {
        double dn = delta;
        for(int j = 0; j<i; ++j) dn *= p;
        b.ene[i] = ene_center - dn;

        dn = delta;
        for(int j = 0; j<n-i-1; j++) dn *= p;
        b.ene[n+i+1] = ene_center + dn;
    }

    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search start" << std::endl;

    std::vector<std::thread> threads;
    threads.resize(thread_num);
    for (int i_thread=0; i_thread<thread_num; i_thread++) {
        auto func =[](int i_thread, band& b, int valley, int band_index) {
            for(int i_mu=i_thread; i_mu<b.mesh; i_mu=i_mu+thread_num) {
                mtx.lock();
//                std::cout << "mu = " << b.ene[i_mu] << std::endl;
                mtx.unlock();
                b.tri[i_mu] = get_triangles_L(valley, band_index, b.ene[i_mu]);
                b.dos[i_mu] = get_DOS_L(b.tri[i_mu], valley, band_index, b.ene[i_mu]);
            }
        };
        threads[i_thread] = std::thread(func, i_thread, std::ref(b), valley, band_index);
    }
    for(auto& thread : threads){
        thread.join();
    }

    std::cout << "L point; valley#" << valley << ", band#" << band_index << " search end" << std::endl;

    return b;
}; // }}}

void init_band_L(band& b, int valley, int band_index) { // {{{
    b.index = band_index;
    b.mesh  = 0;
}; // }}}

void add_fs_L(band& b, int valley, chemical_potential mu) { // {{{
    triangles tri = get_triangles_L(valley, b.index, mu);
    double dos = get_DOS_L(tri, valley, b.index, mu);

    b.mesh += 1;
    b.ene.push_back(mu);
    b.tri.push_back(tri);
    b.dos.push_back(dos);
    std::cout << "hoge" << std::endl;
}; // }}}
