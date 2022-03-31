#include "parameters.hpp"
#include <random>

std::vector<std::string> split(const std::string& input, char delimiter) { // {{{
    std::istringstream stream(input);

    std::string field;
    std::vector<std::string> result;
    while (std::getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}; // }}}

triangles get_triangles(fermi_surface fs) { // {{{
    triangles tri;

    int size = fs.kset.size();
    if (size > 0) {
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<> dst(0, 100000);
        int id = dst(engine);

        std::cout << id << std::endl;

        std::string input = "./fs"+std::to_string(id)+".csv";
        std::string output = "./triangles"+std::to_string(id)+".obj";
        fermi_surface_write(fs, input);
        std::string exec = "python ./py/get_triangles.py "+input+" "+output;
        system(exec.c_str());

        std::ifstream ifs(output);
        if (ifs.fail()) {
            std::cerr << "Failed to open " << output << std::endl;
        }
        std::string buf;
        while (std::getline(ifs, buf)) {
            std::vector<std::string> data;
            data = split(buf, ' ');
            if (data[0] == "#") {
                continue;
            } else if (data[0] == "v") {
                vector3 v = {std::stod(data[1]), std::stod(data[2]), std::stod(data[3])};
                tri.vertexes.push_back(v);
            } else if (data[0] == "vn") {
                vector3 vn = {std::stod(data[1]), std::stod(data[2]), std::stod(data[3])};
                tri.normals.push_back(vn);
            } else if (data[0] == "f") {
                std::vector<std::string> f1, f2, f3;
                f1 = split(data[1], '/');
                f2 = split(data[2], '/');
                f3 = split(data[3], '/');
                face f  = {std::stoi(f1[0]), std::stoi(f2[0]), std::stoi(f3[0])};
                tri.faces.push_back(f);
            }
        }

        if (tri.vertexes.size() != tri.normals.size())
            std::cerr << "Failed to get triangles; vertexes# " << tri.vertexes.size() << ", normals# " << tri.normals.size() << std::endl;

        exec = "rm "+input+" "+output;
        system(exec.c_str());
    }

    return tri;
}; // }}}

// for T point {{{
double get_E_T(int band_index, kpoint k) { // {{{
    matrixComplex H_T = set_T(k.vec);
    vectorReal E_T = diagonalize_N(H_T);
    return E_T[band_index];
}; // }}}

kpoint bisec_T2points(int band_index, chemical_potential mu, double ene, kpoint k1, kpoint k2) { // {{{
    bool flag = false;
    kpoint p, p1, p2;
    double e;
    int loop_max = 50;
    p1 = k1;
    p2 = k2;

    for(int i=0; i<loop_max; i++) {
        for (int axis=0; axis<space_dim; axis++) {
            p.vec[axis] = (p1.vec[axis] + p2.vec[axis])*0.5e0;
        }
        e = get_E_T(band_index, p) - mu;
        if ( std::abs(e) < eps_phys ) {
            flag = true;
            break;
        }

        if (e * ene > 0e0) {
            p1 = p;
        } else {
            p2 = p;
        }
    }
    if (flag == true) {
        return p;
    } else {
        std::cerr << "Error, bisec_T2points did not succeeded." << std::endl;
        return {0e0,0e0,0e0};
    }
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

triangles get_triangles_T(int band_index, chemical_potential mu) { // {{{
    fermi_surface fs = get_fermi_surace_T(band_index, mu);

    triangles tri = get_triangles(fs);

    if (tri.vertexes.size() == 0) {
        return tri;
    }

    int size = tri.vertexes.size();
    double mean = 0e0;
    face f = tri.faces[1];
    for(int axis=0; axis<space_dim; axis++) {
        mean += (tri.vertexes[f.face[1]].vec[axis]-tri.vertexes[f.face[0]].vec[axis])*(tri.vertexes[f.face[1]].vec[axis]-tri.vertexes[f.face[0]].vec[axis]);
    }
    mean = std::sqrt(mean);

    for (int i=0; i<size; i++) {
        double ene = get_E_T(band_index, tri.vertexes[i]) - mu;
        if (std::abs(ene) > eps_phys) {
//            std::cout << std::scientific << "#" << i << ", 誤差が大きい: e = " << ene << std::endl;

            velocity v = get_velocity_T(band_index, mu, tri.vertexes[i]);
            double norm = 0e0;
            for(int axis=0; axis<space_dim; axis++) {
                norm += v.vec[axis]*v.vec[axis];
            }
            norm = std::sqrt(norm);

            kpoint k;
            double pm = ene / std::abs(ene);
            for(int axis=0; axis<space_dim; axis++) {
                k.vec[axis] = tri.vertexes[i].vec[axis] - pm * v.vec[axis]/norm*mean;
            }

            double e = get_E_T(band_index, k);
            if ( e*ene < 0e0 ) {
                k = bisec_T2points(band_index, mu, e, k, tri.vertexes[i]);
                e = get_E_T(band_index, k);
                tri.vertexes[i] = k;
                tri.normals[i] = get_velocity_T(band_index, mu, k);
//                std::cout << std::scientific << "#" << i << ", 再構成OK: e = " << e << std::endl;
            }
        }

    }

    return tri;
}; // }}}

// }}}

// for L points {{{
double get_E_L(int valley, int band_index, kpoint k) { // {{{
    matrixComplex H_L = set_L(valley, k.vec);
    vectorReal E_L = diagonalize_N(H_L);
    return E_L[band_index];
}; // }}}

kpoint bisec_L2points(int valley, int band_index, chemical_potential mu, double ene, kpoint k1, kpoint k2) { // {{{
    bool flag = false;
    kpoint p, p1, p2;
    double e;
    int loop_max = 50;
    p1 = k1;
    p2 = k2;

    for(int i=0; i<loop_max; i++) {
        for (int axis=0; axis<space_dim; axis++) {
            p.vec[axis] = (p1.vec[axis] + p2.vec[axis])*0.5e0;
        }
        e = get_E_L(valley, band_index, p) - mu;
        if ( std::abs(e) < eps_phys ) {
            flag = true;
            break;
        }

        if (e * ene > 0e0) {
            p1 = p;
        } else {
            p2 = p;
        }
    }
    if (flag == true) {
        return p;
    } else {
        std::cerr << "Error, bisec_L2points did not succeeded." << std::endl;
        return {0e0,0e0,0e0};
    }
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

triangles get_triangles_L(int valley, int band_index, chemical_potential mu) { // {{{
    fermi_surface fs = get_fermi_surace_L(valley, band_index, mu);

    triangles tri = get_triangles(fs);

    int size = tri.vertexes.size();
    double mean = 0e0;
    face f = tri.faces[1];
    for(int axis=0; axis<space_dim; axis++) {
        mean += (tri.vertexes[f.face[1]].vec[axis]-tri.vertexes[f.face[0]].vec[axis])*(tri.vertexes[f.face[1]].vec[axis]-tri.vertexes[f.face[0]].vec[axis]);
    }
    mean = std::sqrt(mean);

    for (int i=0; i<size; i++) {
        double ene = get_E_L(valley, band_index, tri.vertexes[i]) - mu;
        if (std::abs(ene) > eps_phys) {
//            std::cout << std::scientific << "#" << i << ", 誤差が大きい: e = " << ene << std::endl;

            velocity v = get_velocity_L(valley, band_index, mu, tri.vertexes[i]);
            double norm = 0e0;
            for(int axis=0; axis<space_dim; axis++) {
                norm += v.vec[axis]*v.vec[axis];
            }
            norm = std::sqrt(norm);

            kpoint k;
            double pm = ene / std::abs(ene);
            for(int axis=0; axis<space_dim; axis++) {
                k.vec[axis] = tri.vertexes[i].vec[axis] - pm * v.vec[axis]/norm*mean;
            }

            double e = get_E_T(band_index, k);
            if ( e*ene < 0e0 ) {
                k = bisec_L2points(valley, band_index, mu, e, k, tri.vertexes[i]);
                e = get_E_L(valley, band_index, k);
                tri.vertexes[i] = k;
                tri.normals[i] = get_velocity_L(valley, band_index, mu, k);
//                std::cout << std::scientific << "#" << i << ", 再構成OK: e = " << e << std::endl;
            }
        }

    }

    return tri;
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

int triangles_write(triangles tri, std::string filename) { // {{{
    int size = tri.vertexes.size();
    if (size > 0) {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cout << filename << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            ofs << std::scientific
                << tri.vertexes[i].vec[0] << ", " << tri.vertexes[i].vec[1] << ", " << tri.vertexes[i].vec[2] << ", "
                << tri.normals[i].vec[0] << ", " << tri.normals[i].vec[1] << ", " << tri.normals[i].vec[2] << ", "
                << std::endl;
        }
    }

    return 0;
}; // }}}
