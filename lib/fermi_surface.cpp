#include "parameters.hpp"

// for T point {{{
double get_E_T(int band_index, kpoint k) { // {{{
    matrixComplex H_T = set_T(k.vec);
    vectorReal E_T = diagonalize_N(H_T);
    return E_T[band_index];
}; // }}}

kpoint bisec_T(int band_index, chemical_potential mu, double ene, kpoint k, int axis, double r) { // {{{
    bool flag = false;
    kpoint p = k;
    double e;
    int loop_max = 50;
    double x1 = p.vec[axis];
    double x2 = r;

    for(int i=0; i<loop_max; i++) {
        p.vec[axis] = (x1 + x2)*0.5e0;
        e = get_E_T(band_index, p) - mu;
        if ( std::abs(e) < eps_phys ) {
            flag = true;
            break;
        }

        if (e * ene > 0e0) {
            x1 = p.vec[axis];
        } else {
            x2 = p.vec[axis];
        }
    }
    if (flag == true) {
        return p;
    } else {
        return {0e0,0e0,0e0};
    }
}; // }}}

velocity get_velocity_T(int band_index, chemical_potential mu, kpoint k) { // {{{
    velocity v = {0e0, 0e0, 0e0};
    for(int axis=0; axis<space_dim; axis++) {
        for(int pm=-1; pm <1; pm=pm+2) {
            double p = double(pm);
            kpoint kp = k;
            kp.vec[axis] += eps_phys*p;
            double ene = get_E_T(band_index, kp) - mu;
            v.vec[axis] += ene*p;
        }
        v.vec[axis] = v.vec[axis] / (2e0*eps_phys);
    }
    return v;
}; // }}}

fermi_surface get_fermi_surace_T(int band_index, chemical_potential mu) { // {{{
    int n = k_mesh;

    fermi_surface fs;
    fs.e = mu;

    double z, x[2*n], y[2*n];
    double ene[2*n][2*n];
    double dn = 1e0/double(n);
    for (int i_z = 0; i_z < 2*n; i_z++) {
        z = cutoff * double(i_z-n)*dn;
        for (int i_y = 0; i_y < 2*n; i_y++) {
            y[i_y] = cutoff * double(i_y-n)*dn;
            for (int i_x = 0; i_x < 2*n; i_x++) {
                x[i_x] = cutoff * double(i_x-n)*dn;
                kpoint p = {x[i_x], y[i_y], z};
                ene[i_y][i_x] = get_E_T(band_index, p) - mu;
            }
        }

        for (int i_y = 0; i_y < 2*n-1; i_y++) {
            for (int i_x = 0; i_x < 2*n-1; i_x++) {
                if (ene[i_y][i_x+1] * ene[i_y][i_x] < 0e0) {
                    int axis = 0;
                    kpoint p = {x[i_x], y[i_y], z};
                    kpoint k = bisec_T(band_index, mu, ene[i_y][i_x], p, axis, x[i_x+1]);
                    fs.kset.push_back(k);
                    fs.vset.push_back(get_velocity_T(band_index, mu, k));
                }

                if (ene[i_y+1][i_x] * ene[i_y][i_x] < 0e0) {
                    int axis = 1;
                    kpoint p = {x[i_x], y[i_y], z};
                    kpoint k = bisec_T(band_index, mu, ene[i_y][i_x], p, axis, y[i_y+1]);
                    fs.kset.push_back(k);
                    fs.vset.push_back(get_velocity_T(band_index, mu, k));
                }
            }
        }

    }

    return fs;
}; // }}}
// }}}

// for L points {{{
double get_E_L(int valley, int band_index, kpoint k) { // {{{
    matrixComplex H_L = set_L(valley, k.vec);
    vectorReal E_L = diagonalize_N(H_L);
    return E_L[band_index];
}; // }}}

kpoint bisec_L(int valley, int band_index, chemical_potential mu, double ene, kpoint k, int axis, double r) { // {{{
    bool flag = false;
    kpoint p = k;
    double e;
    int loop_max = 50;
    double x1 = p.vec[axis];
    double x2 = r;

    for(int i=0; i<loop_max; i++) {
        p.vec[axis] = (x1 + x2)*0.5e0;
        e = get_E_L(valley, band_index, p) - mu;
        if ( std::abs(e) < eps_phys ) {
            flag = true;
            break;
        }

        if (e * ene > 0e0) {
            x1 = p.vec[axis];
        } else {
            x2 = p.vec[axis];
        }
    }
    if (flag == true) {
        return p;
    } else {
        return {0e0,0e0,0e0};
    }
}; // }}}

velocity get_velocity_L(int valley, int band_index, chemical_potential mu, kpoint k) { // {{{
    velocity v = {0e0, 0e0, 0e0};
    for(int axis=0; axis<space_dim; axis++) {
        for(int pm=-1; pm <1; pm=pm+2) {
            double p = double(pm);
            kpoint kp = k;
            kp.vec[axis] += eps_phys*p;
            double ene = get_E_L(valley, band_index, kp) - mu;
            v.vec[axis] += ene*p;
        }
        v.vec[axis] = v.vec[axis] / (2e0*eps_phys);
    }
    return v;
}; // }}}

fermi_surface get_fermi_surace_L(int valley, int band_index, chemical_potential mu) { // {{{
    int n = k_mesh;

    fermi_surface fs;
    fs.e = mu;

    double z, x[2*n], y[2*n];
    double ene[2*n][2*n];
    double dn = 1e0/double(n);
    for (int i_z = 0; i_z < 2*n; i_z++) {
        z = cutoff * double(i_z-n)*dn;
        for (int i_y = 0; i_y < 2*n; i_y++) {
            y[i_y] = cutoff * double(i_y-n)*dn;
            for (int i_x = 0; i_x < 2*n; i_x++) {
                x[i_x] = cutoff * double(i_x-n)*dn;
                kpoint p = {x[i_x], y[i_y], z};
                ene[i_y][i_x] = get_E_L(valley, band_index, p) - mu;
            }
        }

        for (int i_y = 0; i_y < 2*n-1; i_y++) {
            for (int i_x = 0; i_x < 2*n-1; i_x++) {
                if (ene[i_y][i_x+1] * ene[i_y][i_x] < 0e0) {
                    int axis = 0;
                    kpoint p = {x[i_x], y[i_y], z};
                    kpoint k = bisec_L(valley, band_index, mu, ene[i_y][i_x], p, axis, x[i_x+1]);
                    fs.kset.push_back(k);
                    fs.vset.push_back(get_velocity_L(valley, band_index, mu, k));
                }

                if (ene[i_y+1][i_x] * ene[i_y][i_x] < 0e0) {
                    int axis = 1;
                    kpoint p = {x[i_x], y[i_y], z};
                    kpoint k = bisec_L(valley, band_index, mu, ene[i_y][i_x], p, axis, y[i_y+1]);
                    fs.kset.push_back(k);
                    fs.vset.push_back(get_velocity_L(valley, band_index, mu, k));
                }
            }
        }

    }


    return fs;
}; // }}}
// }}}

int fermi_surface_write(fermi_surface fs, std::string filename) { // {{{
    int size = fs.kset.size();
    if (size > 0) {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cout << filename << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            ofs << std::scientific
                << fs.kset[i].vec[0] << ", " << fs.kset[i].vec[1] << ", " << fs.kset[i].vec[2] << ", "
                << fs.vset[i].vec[0] << ", " << fs.vset[i].vec[1] << ", " << fs.vset[i].vec[2] << ", "
                << fs.e
                << std::endl;
        }
    }

    return 0;
}; // }}}
