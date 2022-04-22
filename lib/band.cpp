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

band set_band_2n_L(int valley, int band_index, Energy ene_center, Energy delta, int n, double power) { // {{{
    band b;
    b.index = band_index;
    b.mesh  = 2*n+1;
    b.ene.resize(b.mesh);
    b.tri.resize(b.mesh);
    b.dos.resize(b.mesh);

    b.ene[n] = ene_center;
    for(int i=0; i<n; i++) {
        double dn = delta;
        for(int j = 0; j<i; ++j) dn *= power;
        dn = dn * double(n-i)/double(n);
        b.ene[i] = ene_center - dn;

        dn = delta;
        for(int j = 0; j<n-i-1; j++) dn *= power;
        dn = dn * double(i+1)/double(n);
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
}; // }}}

band combine_band(band b1, band b2) { // {{{
    band b;
    b.index = b1.index;
    b.mesh  = b1.mesh + b2.mesh;

    for(int i=0; i<b1.mesh; i++) {
        b.ene.push_back(b1.ene[i]);
        b.tri.push_back(b1.tri[i]);
        b.dos.push_back(b1.dos[i]);
    }
    for(int i=0; i<b2.mesh; i++) {
        b.ene.push_back(b2.ene[i]);
        b.tri.push_back(b2.tri[i]);
        b.dos.push_back(b2.dos[i]);
    }

    return b;
}; // }}}

band combine_band_2n(band b_global, band b_local) { // {{{
    band b;
    b.index = b_global.index;

    int index = 0;
    for(int i=0; i<b_global.mesh; i++) {
        if( b_global.ene[i] > b_local.ene[0] ) {
            index = i;
            break;
        }
        b.ene.push_back(b_global.ene[i]);
        b.tri.push_back(b_global.tri[i]);
        b.dos.push_back(b_global.dos[i]);
    }
    for(int i=0; i<b_local.mesh; i++) {
        b.ene.push_back(b_local.ene[i]);
        b.tri.push_back(b_local.tri[i]);
        b.dos.push_back(b_local.dos[i]);
    }
    for(int i=index+1; i<b_global.mesh; i++) {
        if( b_global.ene[i] < b_local.ene[b_local.mesh-1] ) continue;
        b.ene.push_back(b_global.ene[i]);
        b.tri.push_back(b_global.tri[i]);
        b.dos.push_back(b_global.dos[i]);
    }

    b.mesh  = b.tri.size();
    return b;
}; // }}}
