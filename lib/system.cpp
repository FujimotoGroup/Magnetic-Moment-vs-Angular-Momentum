#include "parameters.hpp"
#include <functional>

int thread_num = std::thread::hardware_concurrency();

sys get_T() {
    std::vector<std::thread> threads;

    sys s;
    chemical_potential mu_min[bandsT], mu_max[bandsT];
    int mu_mesh = 10;
    for(int i=0; i<bandsT; i++) {
        mu_max[i] = double(ET[i])+1e-1;
        mu_min[i] = double(ET[i])-1e-1;
    }

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
}

void sys_write(sys s) {
    for(int band_index=0; band_index<bandsT; band_index++) {
        for(int i_mu=0; i_mu<s.bands[band_index].mu_mesh; i_mu++) {
            std::cout << i_mu << std::endl;
            std::string filename = "dat/band"+std::to_string(band_index)+"k"+std::to_string(i_mu)+".csv";
            fermi_surface_write(s.bands[band_index].fs[i_mu], filename);
        }
    }
}
