#include "parameters.hpp"

int thread_num = std::thread::hardware_concurrency();

// set T {{{
sys_T get_T() { // {{{
    sys_T s;
    chemical_potential mu_min[bandsT], mu_max[bandsT];
    for(int i=0; i<bandsT; i++) {
        mu_max[i] = double(ET[i])+5e-2;
        mu_min[i] = double(ET[i])-5e-2;
    }

    std::vector<std::thread> threads;
    std::vector<band> b;
    b.resize(bandsT);
    for(int i=0; i<bandsT; i++) {
        threads.emplace_back(std::thread(get_band_T, std::ref(b[i]), i, mu_min[i], mu_max[i], mu_mesh));
    }


    for(auto& thread : threads){
        thread.join();
    }

    s.bands = b;

    return s;
}; // }}}

void sys_T_write(sys_T s) { // {{{
    for(int band_index=0; band_index<bandsT; band_index++) {
        for(int i_mu=0; i_mu<s.bands[band_index].mu_mesh; i_mu++) {
            std::string filename = "dat/T_band"+std::to_string(band_index)+"k"+std::to_string(i_mu)+".csv";
            fermi_surface_write(s.bands[band_index].fs[i_mu], filename);
        }
    }
}; // }}}
// }}}

// set L {{{
sys_L get_L() { // {{{

    sys_L s;
    chemical_potential mu_min[bandsL], mu_max[bandsL];
    const int mu_mesh = 10;
    s.bands.resize(valleys);
    for(int valley=0; valley<valleys; valley++) {
        for(int i=0; i<bandsL; i++) {
            mu_max[i] = double(EL[valley][i])+1e-1;
            mu_min[i] = double(EL[valley][i])-1e-1;
        }

        std::vector<std::thread> threads;
        std::vector<band> b;
        b.resize(bandsL);

        for(int i=0; i<bandsL; i++) {
            threads.emplace_back(std::thread(get_band_L, std::ref(b[i]), valley, i, mu_min[i], mu_max[i], mu_mesh));
        }

        for(auto& thread : threads){
            thread.join();
        }
        s.bands[valley] = b;
    };

    return s;
}; // }}}

void sys_L_write(sys_L s) { // {{{
    for(int valley=0; valley<valleys; valley++) {
        for(int band_index=0; band_index<bandsL; band_index++) {
            for(int i_mu=0; i_mu<s.bands[valley][band_index].mu_mesh; i_mu++) {
                std::string filename = "dat/L"+std::to_string(valley)+"_band"+std::to_string(band_index)+"k"+std::to_string(i_mu)+".csv";
                fermi_surface_write(s.bands[valley][band_index].fs[i_mu], filename);
            }
        }
    }
}; // }}}
// }}}
