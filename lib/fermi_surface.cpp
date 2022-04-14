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

double NRsqrt(double a) { // {{{
    return std::sqrt(a);

//    double x=1e0;
//    int count = 0;
//
//    if (a < 0e0) {
//        std::cerr << "NRsqrt is not valid for negative double" << std::endl;
//    } else if (a == 0e0) {
//        return a;
//    }
//
//    do {
//        count++;
//        x = ( x + a/x ) * 5e-1;
//        if (count > 100000) {
//            std::cerr << "sqrt function does not converge." << std::endl;
//            break;
//        }
//    } while (std::abs(a - x*x) > 1e-13);
//
//    return x;
} // }}}

triangles get_triangles(fermi_surface fs) { // {{{
    triangles tri;
    tri.ene = fs.e;

    int size = fs.kset.size();
    if (size > 0) {
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_int_distribution<> dst(0, 100000);
        int id = dst(engine);

//        std::cout << id << std::endl;

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

                face f;
                f.face[0] = std::stoi(f1[0])-1;
                f.face[1] = std::stoi(f2[0])-1;
                f.face[2] = std::stoi(f3[0])-1;
                tri.faces.push_back(f);
            }
        }

        if (tri.vertexes.size() == 0)
            std::cerr << "Failed to get triangles; vertexes# " << tri.vertexes.size() << ", normals# " << tri.normals.size() << std::endl;

        exec = "rm "+input+" "+output;
        system(exec.c_str());
    }

    return tri;
}; // }}}

double get_length(kpoint k1, kpoint k2) { // {{{
    double len;
    for(int axis=0; axis<space_dim; axis++) {
        len = (k1.vec[axis] - k2.vec[axis])*(k1.vec[axis] - k2.vec[axis]);
    }
    len = NRsqrt(len);
    return len;
} // }}}

kpoint get_center(kpoint k1, kpoint k2, kpoint k3) { // {{{
    kpoint center;
    for(int axis=0; axis<space_dim; axis++) {
        center.vec[axis] = (k1.vec[axis] + k2.vec[axis] + k3.vec[axis])/3e0;
    }
    return center;
} // }}}

kpoint get_normal(kpoint n1, kpoint n2, kpoint n3) { // {{{
    kpoint normal;
    for(int axis=0; axis<space_dim; axis++) {
        normal.vec[axis] = (n1.vec[axis] + n2.vec[axis] + n3.vec[axis])/3e0;
    }
    return normal;
} // }}}

double get_dS(kpoint k1, kpoint k2, kpoint k3) { // {{{
    double p = 0e0, q = 0e0, r = 0e0;
    for(int axis=0; axis<space_dim; axis++) {
        p += (k1.vec[axis] - k2.vec[axis])*(k1.vec[axis] - k2.vec[axis]);
        q += (k2.vec[axis] - k3.vec[axis])*(k2.vec[axis] - k3.vec[axis]);
        r += (k3.vec[axis] - k1.vec[axis])*(k3.vec[axis] - k1.vec[axis]);
    }
    p = NRsqrt(p);
    q = NRsqrt(q);
    r = NRsqrt(r);
    double s = (p + q + r)*5e-1;
    return NRsqrt( s*(s-p)*(s-q)*(s-r) );
} // }}}

double get_min_length(kpoint k1, kpoint k2, kpoint k3) { // {{{
    double p = 0e0, q = 0e0, r = 0e0;
    for(int axis=0; axis<space_dim; axis++) {
        p += (k1.vec[axis] - k2.vec[axis])*(k1.vec[axis] - k2.vec[axis]);
        q += (k2.vec[axis] - k3.vec[axis])*(k2.vec[axis] - k3.vec[axis]);
        r += (k3.vec[axis] - k1.vec[axis])*(k3.vec[axis] - k1.vec[axis]);
    }
    double min = std::min(std::min(p, q), r);
    return NRsqrt(min);

//    p = NRsqrt(p);
//    q = NRsqrt(q);
//    r = NRsqrt(r);
//    return std::min(std::min(p, q), r);
} // }}}

// for T point {{{
double get_E_T(int band_index, kpoint k) { // {{{
    matrixComplex H_T = set_T(k.vec);
    vectorReal E_T = diagonalize_N(H_T);
    return E_T[band_index];
}; // }}}

kpoint bisec_T2points(int band_index, chemical_potential mu, double ene1, double ene2, kpoint k1, kpoint k2) { // {{{
    bool flag = false;
    kpoint p, p1, p2;
    double e;
    int loop_max = 50;
    p1 = k1;
    p2 = k2;
//    std::cout << "p1 = " << p1.vec[0] << ", " << p1.vec[1] << ", " << p1.vec[2] << std::endl;
//    std::cout << "p2 = " << p2.vec[0] << ", " << p2.vec[1] << ", " << p2.vec[2] << std::endl;

    for(int i=0; i<loop_max; i++) {
        for (int axis=0; axis<space_dim; axis++) {
            p.vec[axis] = (p1.vec[axis] + p2.vec[axis])*0.5e0;
        }
        e = get_E_T(band_index, p) - mu;

//        double diff = 0e0;
//        for(int axis=0; axis<space_dim; axis++)
//            diff += (p1.vec[axis] - p2.vec[axis])*(p1.vec[axis] - p2.vec[axis]);
//        diff = std::sqrt(diff);
//        std::cout << "#" << i << ": diff = " << diff << ", e = " << e <<std::endl;

        if ( std::abs(e) < eps_phys ) {
            flag = true;
            break;
        }

        if (e * ene1 > 0e0) {
            p1 = p;
        } else if (e * ene2 > 0e0)  {
            p2 = p;
        } else {
            std::cerr << "Error!!!!" << std::endl;
        }
    }
    if (flag == true) {
        return p;
    } else {
        std::cout << "Error, bisec_T2points did not succeeded." << std::endl;
        std::cout << "k1 = " << k1.vec[0] << ", " << k1.vec[1] << ", " << k1.vec[2] << ", e = " << ene1 << std::endl;
        std::cout << "k2 = " << k2.vec[0] << ", " << k2.vec[1] << ", " << k2.vec[2] << ", e = " << ene2 << std::endl;
        return {0e0, 0e0, 0e0};
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
        for(int i=0; i<2; i++) {
            double p = double(2*i-1);
            kpoint kp = k;
            kp.vec[axis] += eps_phys*p;
            double ene = get_E_T(band_index, kp);
            v.vec[axis] += ene*p;
        }
        v.vec[axis] = v.vec[axis] / (2e0*eps_phys);
    }
    return v;
}; // }}}

fermi_surface get_fermi_surface_more_T(fermi_surface fs, int band_index, chemical_potential mu) { // {{{
    int n = k_mesh_more;

    fermi_surface new_fs;
    new_fs.e = mu;

    double x_max = 0e0, x_min = 0e0;
    double y_max = 0e0, y_min = 0e0;
    double z_max = 0e0, z_min = 0e0;
    int size = fs.kset.size();
    for(int i=0; i<size; i++) {
        x_max = std::max(x_max, fs.kset[i].vec[0]);
        x_min = std::min(x_min, fs.kset[i].vec[0]);
        y_max = std::max(y_max, fs.kset[i].vec[1]);
        y_min = std::min(y_min, fs.kset[i].vec[1]);
        z_max = std::max(z_max, fs.kset[i].vec[2]);
        z_min = std::min(z_min, fs.kset[i].vec[2]);
    }

    vectorReal len = vectorReal(3, 0e0);
    len[0] = x_max - x_min;
    len[1] = y_max - y_min;
    len[2] = z_max - z_min;

    for(int axis=0; axis<space_dim; axis++) {
        if (len[axis] == 0e0) len[axis] = 2e0*dk[axis];
    }

    kpoint center;
    center.vec[0] = (x_max + x_min)* 5e-1;
    center.vec[1] = (y_max + y_min)* 5e-1;
    center.vec[2] = (z_max + z_min)* 5e-1;

//    std::cout << "center = " << center.vec[0] << ", " << center.vec[1] << ", " << center.vec[2] << std::endl;
//    std::cout << "length = " << len[0] << ", " << len[1] << ", " << len[2] << std::endl;

    double z, x[2*n], y[2*n];
    double ene[2*n][2*n];
    double dn = 1.3e0/double(n);
    for (int i_z = 0; i_z < 2*n; i_z++) {
        z = center.vec[2] + double(i_z-n)*dn*len[2];
        for (int i_y = 0; i_y < 2*n; i_y++) {
            y[i_y] = center.vec[1] + double(i_y-n)*dn*len[1];
            for (int i_x = 0; i_x < 2*n; i_x++) {
                x[i_x] = center.vec[0] + double(i_x-n)*dn*len[0];
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
                    new_fs.kset.push_back(k);
                    new_fs.vset.push_back(get_velocity_T(band_index, mu, k));
                }

                if (ene[i_y+1][i_x] * ene[i_y][i_x] < 0e0) {
                    int axis = 1;
                    kpoint p = {x[i_x], y[i_y], z};
                    kpoint k = bisec_T(band_index, mu, ene[i_y][i_x], p, axis, y[i_y+1]);
                    new_fs.kset.push_back(k);
                    new_fs.vset.push_back(get_velocity_T(band_index, mu, k));
                }
            }
        }

    }

//    std::cout << "new_fs#" << new_fs.kset.size() << std::endl;
    return new_fs;
}; // }}}

fermi_surface get_fermi_surace_T(int band_index, chemical_potential mu) { // {{{
    int n = k_mesh;

    fermi_surface fs;
    fs.e = mu;

    double z, x[2*n], y[2*n];
    double ene[2*n][2*n];
    double dn = 1e0/double(n);
    for(int axis=0; axis<space_dim; axis++) {
        dk[axis] = cutoff*dn;
    }

    for (int i_z = 0; i_z < 2*n; i_z++) {
        z = double(i_z-n)*dk[2];
        for (int i_y = 0; i_y < 2*n; i_y++) {
            y[i_y] = double(i_y-n)*dk[1];
            for (int i_x = 0; i_x < 2*n; i_x++) {
                x[i_x] = double(i_x-n)*dk[0];
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

    int size_fs = fs.kset.size();
    if (size_fs < 5) {
        fs.kset.resize(0);
        return fs;
    }

    if ((size_fs > 5) && (size_fs < 1000)) {
        do {
            std::cout << "re-constructing Fermi surface because of poor mesh number." << std::endl;
            fs = get_fermi_surface_more_T(fs, band_index, mu);
            size_fs = fs.kset.size();
            std::cout << "fs#" << size_fs << std::endl;
        } while (fs.kset.size() < 400);
        size_fs = fs.kset.size();
    }

    std::cout << "final fs#" << size_fs << std::endl;

    return fs;
}; // }}}

triangles get_triangles_T(int band_index, chemical_potential mu) { // {{{
    fermi_surface fs = get_fermi_surace_T(band_index, mu);
    int size_fs = fs.kset.size();

    triangles tri;
    if (size_fs == 0) {
        tri.ene = fs.e;
        return tri;
    }

    tri = get_triangles(fs);

    int size;
    size = tri.vertexes.size();
    if (tri.normals.size() == 0) tri.normals.resize(size);

    for (int i=0; i<size; i++) {
        tri.normals[i] = get_velocity_T(band_index, mu, tri.vertexes[i]);
    }

    size = tri.faces.size();
    for (int i=0; i<size; i++) {
        int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
        vector3 k1 = tri.vertexes[v[0]];
        vector3 k2 = tri.vertexes[v[1]];
        vector3 k3 = tri.vertexes[v[2]];
        vector3 n1 = tri.normals[v[0]];
        vector3 n2 = tri.normals[v[1]];
        vector3 n3 = tri.normals[v[2]];
        for(int axis=0; axis<space_dim; axis++) {
            tri.faces[i].center[axis] = k1.vec[axis];
            tri.faces[i].normal[axis] = n1.vec[axis];
        }
        tri.faces[i].dS = get_dS(k1, k2, k3);
    }

    return tri;

//    for (int i=0; i<size; i++) {
//        double ene = get_E_T(band_index, tri.vertexes[i]) - mu;
//        if (std::abs(ene) > eps_phys) {
////            std::cout << std::scientific << "#" << i << ", 誤差が大きい: e = " << ene << std::endl;
//
//            velocity v = get_velocity_T(band_index, mu, tri.vertexes[i]);
//            double norm = 0e0;
//            for(int axis=0; axis<space_dim; axis++) {
//                norm += v.vec[axis]*v.vec[axis];
//            }
//            norm = NRsqrt(norm);
//
//            double mean;
//            int f[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
//            vector3 k1 = tri.vertexes[f[0]];
//            vector3 k2 = tri.vertexes[f[1]];
//            vector3 k3 = tri.vertexes[f[2]];
//            double p = 0e0, q = 0e0, r = 0e0;
//            for(int axis=0; axis<space_dim; axis++) {
//                p += (k1.vec[axis] - k2.vec[axis])*(k1.vec[axis] - k2.vec[axis]);
//                q += (k2.vec[axis] - k3.vec[axis])*(k2.vec[axis] - k3.vec[axis]);
//                r += (k3.vec[axis] - k1.vec[axis])*(k3.vec[axis] - k1.vec[axis]);
//            }
//            mean = std::min(NRsqrt(p), NRsqrt(q));
//            mean = std::min(mean, NRsqrt(r)) * 5e-2;
//
//            kpoint k;
//            double pm = ene / std::abs(ene);
//            for(int axis=0; axis<space_dim; axis++) {
//                k.vec[axis] = tri.vertexes[i].vec[axis] - pm * v.vec[axis]/norm*mean*2e0;
//            }
//            double ek = get_E_T(band_index, k) - mu;
//            if ( ek*ene < 0e0 ) {
////                std::cout << std::scientific << "ek = " << ek << ", ene = " << ene << std::endl;
//                kpoint q = bisec_T2points(band_index, mu, ek, ene, k, tri.vertexes[i]);
//                ek = get_E_T(band_index, q) - mu;
////                std::cout << std::scientific << "#" << i << ", 再構成OK: ek = " << ek << std::endl;
//                tri.vertexes[i] = q;
//
//            }
//        }
//
//        tri.normals[i] = get_velocity_T(band_index, mu, tri.vertexes[i]);
//
//    }
//
//    size = tri.faces.size();
//    for (int i=0; i<size; i++) {
//        int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
//        vector3 k1 = tri.vertexes[v[0]];
//        vector3 k2 = tri.vertexes[v[1]];
//        vector3 k3 = tri.vertexes[v[2]];
//        vector3 n1 = tri.normals[v[0]];
//        vector3 n2 = tri.normals[v[1]];
//        vector3 n3 = tri.normals[v[2]];
//
//        double p = 0e0, q = 0e0, r = 0e0;
//        for(int axis=0; axis<space_dim; axis++) {
//            tri.faces[i].center[axis] = (k1.vec[axis] + k2.vec[axis] + k3.vec[axis])/3e0;
////            tri.faces[i].normal[axis] = (n1.vec[axis] + n2.vec[axis] + n3.vec[axis])/3e0;
//            tri.faces[i].normal[axis] = n1.vec[axis];
//            p += (k1.vec[axis] - k2.vec[axis])*(k1.vec[axis] - k2.vec[axis]);
//            q += (k2.vec[axis] - k3.vec[axis])*(k2.vec[axis] - k3.vec[axis]);
//            r += (k3.vec[axis] - k1.vec[axis])*(k3.vec[axis] - k1.vec[axis]);
//        }
//        p = NRsqrt(p);
//        q = NRsqrt(q);
//        r = NRsqrt(r);
//        double s = (p + q + r)*5e-1;
//        tri.faces[i].dS = NRsqrt( s*(s-p)*(s-q)*(s-r) );
//
////        double n = 0e0;
////        for(int axis=0; axis<space_dim; axis++) {
////            k1.vec[axis] = k1.vec[axis] - k3.vec[axis];
////            k2.vec[axis] = k2.vec[axis] - k3.vec[axis];
////        }
////        k3.vec[0] = k1.vec[1]*k2.vec[2] - k1.vec[2]*k2.vec[1];
////        k3.vec[1] = k1.vec[2]*k2.vec[0] - k1.vec[0]*k2.vec[2];
////        k3.vec[2] = k1.vec[0]*k2.vec[1] - k1.vec[1]*k2.vec[0];
////        for(int axis=0; axis<space_dim; axis++) {
////            n += k3.vec[axis]*k3.vec[axis];
////        }
////        n = NRsqrt(n)*5e-1;
//    }
//
//
//    return tri;
}; // }}}

double get_DOS_T(triangles tri, int band_index, chemical_potential mu) { // {{{
    double dos = 0e0;
    auto fn = [](int band_index, chemical_potential mu, kpoint k) { return 1e0; };
    integrate_triangles_T(fn, dos, tri, band_index, mu);

    return dos;
}; // }}}

int triangles_write_T(triangles tri, std::string name) { // {{{
    int size = tri.vertexes.size();
    std::string vertex_name = name+"_vertex.csv";
    if (size > 0) {
        std::ofstream ofs(vertex_name);
        if (!ofs) {
            std::cout << vertex_name << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            ofs << std::scientific
                << kT[0]+tri.vertexes[i].vec[0] << ", "
                << kT[1]+tri.vertexes[i].vec[1] << ", "
                << kT[2]+tri.vertexes[i].vec[2] << ", "
                << tri.normals[i].vec[0] << ", " << tri.normals[i].vec[1] << ", " << tri.normals[i].vec[2] << ", "
                << std::endl;
        }
    }

    size = tri.faces.size();
    std::string face_name = name+"_face.csv";
    if (size > 0) {
        std::ofstream ofs(face_name);
        if (!ofs) {
            std::cout << face_name << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {

            int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};

            double norm = 0e0;
            for(int axis=0; axis<space_dim; axis++) {
                norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
            }
            norm = NRsqrt(norm);

            ofs << std::scientific
                << kT[0]+tri.vertexes[v[0]].vec[0] << ", "
                << kT[1]+tri.vertexes[v[0]].vec[1] << ", "
                << kT[2]+tri.vertexes[v[0]].vec[2]
                << std::endl;
            ofs << std::scientific
                << kT[0]+tri.vertexes[v[1]].vec[0] << ", "
                << kT[1]+tri.vertexes[v[1]].vec[1] << ", "
                << kT[2]+tri.vertexes[v[1]].vec[2]
                << std::endl;
            ofs << std::scientific
                << kT[0]+tri.vertexes[v[2]].vec[0] << ", "
                << kT[1]+tri.vertexes[v[2]].vec[1] << ", "
                << kT[2]+tri.vertexes[v[2]].vec[2]
                << std::endl;
            ofs << std::scientific
                << kT[0]+tri.vertexes[v[0]].vec[0] << ", "
                << kT[1]+tri.vertexes[v[0]].vec[1] << ", "
                << kT[2]+tri.vertexes[v[0]].vec[2]
                << std::endl;
            ofs << ",," << std::endl;

        }
    }

    return 0;
}; // }}}
// }}}

// for L points {{{
double get_E_L(int valley, int band_index, kpoint k) { // {{{
    matrixComplex H_L = set_L(valley, k.vec);
    vectorReal E_L = diagonalize_N(H_L);
    return E_L[band_index];
}; // }}}

kpoint bisec_L2points(int valley, int band_index, chemical_potential mu, double ene1, double ene2, kpoint k1, kpoint k2) { // {{{
    bool flag = false;
    kpoint p, p1, p2;
    double e;
    int loop_max = 50;
    p1 = k1;
    p2 = k2;
//    std::cout << "p1 = " << p1.vec[0] << ", " << p1.vec[1] << ", " << p1.vec[2] << std::endl;
//    std::cout << "p2 = " << p2.vec[0] << ", " << p2.vec[1] << ", " << p2.vec[2] << std::endl;

    for(int i=0; i<loop_max; i++) {
        for (int axis=0; axis<space_dim; axis++) {
            p.vec[axis] = (p1.vec[axis] + p2.vec[axis])*0.5e0;
        }
        e = get_E_L(valley, band_index, p) - mu;

//        double diff = 0e0;
//        for(int axis=0; axis<space_dim; axis++)
//            diff += (p1.vec[axis] - p2.vec[axis])*(p1.vec[axis] - p2.vec[axis]);
//        diff = std::sqrt(diff);
//        std::cout << "#" << i << ": diff = " << diff << ", e = " << e <<std::endl;

        if ( std::abs(e) < eps_phys ) {
            flag = true;
            break;
        }

        if (e * ene1 > 0e0) {
            p1 = p;
        } else if (e * ene2 > 0e0)  {
            p2 = p;
        } else {
            std::cerr << "Error!!!!" << std::endl;
        }
    }
    if (flag == true) {
        return p;
    } else {
        std::cout << "Error, bisec_T2points did not succeeded." << std::endl;
        std::cout << "k1 = " << k1.vec[0] << ", " << k1.vec[1] << ", " << k1.vec[2] << ", e = " << ene1 << std::endl;
        std::cout << "k2 = " << k2.vec[0] << ", " << k2.vec[1] << ", " << k2.vec[2] << ", e = " << ene2 << std::endl;
        return {0e0, 0e0, 0e0};
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
        for(int i=0; i<2; i++) {
            double p = double(2*i-1);
            kpoint kp = k;
            kp.vec[axis] += eps_phys*p;
            double ene = get_E_L(valley, band_index, kp);
            v.vec[axis] += ene*p;
        }
        v.vec[axis] = v.vec[axis] / (2e0*eps_phys);
    }
    return v;
}; // }}}

fermi_surface get_fermi_surface_more_L(fermi_surface fs, int valley, int band_index, chemical_potential mu) { // {{{
    const int n = k_mesh_more;

    fermi_surface new_fs;
    new_fs.e = mu;

    double x_max = 0e0, x_min = 0e0;
    double y_max = 0e0, y_min = 0e0;
    double z_max = 0e0, z_min = 0e0;
    int size = fs.kset.size();
//    std::cout << "reconstructing fermi surface: " << size << std::endl;
    for(int i=0; i<size; i++) {
        x_max = std::max(x_max, fs.kset[i].vec[0]);
        x_min = std::min(x_min, fs.kset[i].vec[0]);
        y_max = std::max(y_max, fs.kset[i].vec[1]);
        y_min = std::min(y_min, fs.kset[i].vec[1]);
        z_max = std::max(z_max, fs.kset[i].vec[2]);
        z_min = std::min(z_min, fs.kset[i].vec[2]);
    }

    vectorReal len = vectorReal(3, 0e0);
    len[0] = x_max - x_min;
    len[1] = y_max - y_min;
    len[2] = z_max - z_min;

    for(int axis=0; axis<space_dim; axis++) {
        if (len[axis] == 0e0) len[axis] = 2e0*dk[axis];
    }

    kpoint center;
    center.vec[0] = (x_max + x_min)* 5e-1;
    center.vec[1] = (y_max + y_min)* 5e-1;
    center.vec[2] = (z_max + z_min)* 5e-1;

//    std::cout << "center = " << center.vec[0] << ", " << center.vec[1] << ", " << center.vec[2] << std::endl;
//    std::cout << "length = " << len[0] << ", " << len[1] << ", " << len[2] << std::endl;

    double z, x[2*n], y[2*n];
    double ene[2*n][2*n];
    double dn = 1.5e0/double(n);
    for (int i_z = 0; i_z < 2*n; i_z++) {
        z = center.vec[2] + double(i_z-n)*dn*len[2];
//        std::cout << "z = " << z << std::endl;
        for (int i_y = 0; i_y < 2*n; i_y++) {
            y[i_y] = center.vec[1] + double(i_y-n)*dn*len[1];
//            std::cout << "y = " << y[i_y] << std::endl;
            for (int i_x = 0; i_x < 2*n; i_x++) {
                x[i_x] = center.vec[0] + double(i_x-n)*dn*len[0];
//                std::cout << "x = " << x[i_x] << std::endl;
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
                    new_fs.kset.push_back(k);
                    new_fs.vset.push_back(get_velocity_L(valley, band_index, mu, k));
                }

                if (ene[i_y+1][i_x] * ene[i_y][i_x] < 0e0) {
                    int axis = 1;
                    kpoint p = {x[i_x], y[i_y], z};
                    kpoint k = bisec_L(valley, band_index, mu, ene[i_y][i_x], p, axis, y[i_y+1]);
                    new_fs.kset.push_back(k);
                    new_fs.vset.push_back(get_velocity_L(valley, band_index, mu, k));
                }
            }
        }

    }

//    std::cout << "new_fs#" << new_fs.kset.size() << std::endl;
    return new_fs;
}; // }}}

fermi_surface get_fermi_surace_L(int valley, int band_index, chemical_potential mu) { // {{{
    int n = k_mesh;

    fermi_surface fs;
    fs.e = mu;

    double z, x[2*n], y[2*n];
    double ene[2*n][2*n];
    double dn = 1e0/double(n);
    for(int axis=0; axis<space_dim; axis++) {
        dk[axis] = cutoff*dn;
    }

    for (int i_z = 0; i_z < 2*n; i_z++) {
        z = double(i_z-n)*dk[2];
        for (int i_y = 0; i_y < 2*n; i_y++) {
            y[i_y] = double(i_y-n)*dk[1];
            for (int i_x = 0; i_x < 2*n; i_x++) {
                x[i_x] = double(i_x-n)*dk[0];
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

    int size_fs = fs.kset.size();
    if (size_fs < 5) {
        fs.kset.resize(0);
        return fs;
    }

    if ((size_fs > 5) && (size_fs < 1000)) {
        std::cout << "re-constructing Fermi surface because of poor mesh number." << std::endl;
        fs = get_fermi_surface_more_L(fs, valley, band_index, mu);
        size_fs = fs.kset.size();
    }

    std::cout << "final fs#" << size_fs << std::endl;

    return fs;
}; // }}}

void divide_triangle(triangles& tri, int valley, int band_index, chemical_potential mu, int i) { // {{{
    double dS = tri.faces[i].dS;
//    std::cout << i << ", dS = " << dS << std::endl;

    int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
    vector3 k1 = tri.vertexes[v[0]];
    vector3 k2 = tri.vertexes[v[1]];
    vector3 k3 = tri.vertexes[v[2]];
    vector3 n1 = tri.normals[v[0]];
    vector3 n2 = tri.normals[v[1]];
    vector3 n3 = tri.normals[v[2]];

//    std::cout << "e1 = " << get_E_L(valley, band_index, k1) - mu << ", "
//              << "e2 = " << get_E_L(valley, band_index, k2) - mu << ", "
//              << "e3 = " << get_E_L(valley, band_index, k3) - mu << ", "
//              << std::endl;

    kpoint center = get_center(k1, k2, k3);
    double ene = get_E_L(valley, band_index, center) - mu;
//    vector3 c1, c2, normal;
//    for(int axis=0; axis<space_dim; axis++) {
//        c1.vec[axis] = k1.vec[axis] - k3.vec[axis];
//        c2.vec[axis] = k2.vec[axis] - k3.vec[axis];
//    }
//    normal.vec[0] = c1.vec[1]*c2.vec[2] - c1.vec[2]*c2.vec[1];
//    normal.vec[1] = c1.vec[2]*c2.vec[0] - c1.vec[0]*c2.vec[2];
//    normal.vec[2] = c1.vec[0]*c2.vec[1] - c1.vec[1]*c2.vec[0];

    velocity vel = get_velocity_L(valley, band_index, mu, center);

    double norm = 0e0;
    for(int axis=0; axis<space_dim; axis++) {
        norm += vel.vec[axis]*vel.vec[axis];
//        norm += normal.vec[axis]*normal.vec[axis];
    }
    norm = NRsqrt(norm);

    kpoint k;
    double pm = ene / std::abs(ene);
    double min = get_min_length(k1, k2, k3);
    double ek;
    for(int i=1; i<50; i++) {
        for(int axis=0; axis<space_dim; axis++) {
            k.vec[axis] = center.vec[axis] - pm * vel.vec[axis]/norm*min*5e-1*double(i);
//            k.vec[axis] = center.vec[axis] - pm * normal.vec[axis]/norm*min*5e-1*double(i);
        }
        ek = get_E_L(valley, band_index, k) - mu;
        if (ek*ene < 0e0) break;
    }

    if ( ek*ene < 0e0 ) {
//        std::cout << std::scientific << "ek = " << ek << ", ene = " << ene << std::endl;
        kpoint q = bisec_L2points(valley, band_index, mu, ek, ene, k, center);
//        ek = get_E_L(valley, band_index, q) - mu;
//        std::cout << std::scientific << "#" << i << ", 再構成OK: ek = " << ek << std::endl;
        tri.vertexes.push_back(q);
        int v_new = tri.vertexes.size()-1;
        velocity n = get_velocity_L(valley, band_index, mu, q);

        int v[4] = {v_new, tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
//        std::cout << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << std::endl;
        kpoint k[4] = {q, k1, k2, k3};
        int j1, j2, j3;
        double r = get_length(q, k1);
        double s = get_length(q, k2);
        double t = get_length(q, k3);
//        std::cout << "hoge = " << r << ", " << s << ", " << t << std::endl;
        if (r > s) {
            if (r > t) { // r > t, s
                j1 = 1;
                j2 = 2;
                j3 = 3;
            } else { // t > r > s
                j1 = 3;
                j2 = 1;
                j3 = 2;
            }
        } else {
            if (s > t) { // s > r, t
                j1 = 2;
                j2 = 3;
                j3 = 1;
            } else { // t > s > r
                j1 = 3;
                j2 = 1;
                j3 = 2;
            }
        }
        tri.faces[i].face[j1-1] = v[0];
        tri.faces[i].dS = get_dS(k[0], k[j2], k[j3]);
        face f1, f2;
        f1.face[0] = v[0];
        f1.face[1] = v[j1];
        f1.face[2] = v[j2];
        f1.dS = get_dS(k[0], k[j1], k[j2]);
        f2.face[0] = v[0];
        f2.face[1] = v[j3];
        f2.face[2] = v[j1];
        f2.dS = get_dS(k[0], k[j3], k[j1]);

//        kpoint center1  = get_center(q, k1, k2);
//        kpoint center2  = get_center(q, k1, k3);
//        kpoint center3  = get_center(q, k2, k3);
//        vector3 normal1 = get_normal(n, n1, n2);
//        vector3 normal2 = get_normal(n, n1, n3);
//        vector3 normal3 = get_normal(n, n2, n3);

        kpoint center1  = q;
        kpoint center2  = k1;
        kpoint center3  = k2;
        vector3 normal1 = n;
        vector3 normal2 = n1;
        vector3 normal3 = n2;
        for(int axis=0; axis<space_dim; axis++) {
            f1.center[axis] = center1.vec[axis];
            f1.normal[axis] = normal1.vec[axis];
            f2.center[axis] = center2.vec[axis];
            f2.normal[axis] = normal2.vec[axis];
            tri.faces[i].center[axis] = center3.vec[axis];
            tri.faces[i].normal[axis] = normal3.vec[axis];
        }

        tri.faces.push_back(f1);
        tri.faces.push_back(f2);

//        std::cout << i << ", dS1 = " << tri.faces[i].dS <<
//                          ", dS2 = " << f1.dS <<
//                          ", dS3 = " << f2.dS << std::endl;
//

    }
} // }}}

triangles get_triangles_L(int valley, int band_index, chemical_potential mu) { // {{{
    fermi_surface fs = get_fermi_surace_L(valley, band_index, mu);
    int size_fs = fs.kset.size();

    triangles tri;
    if (size_fs == 0) {
        tri.ene = fs.e;
        return tri;
    }

    tri = get_triangles(fs);


    int size;
    size = tri.vertexes.size();
    if (tri.normals.size() == 0) tri.normals.resize(size);

    for (int i=0; i<size; i++) {
        tri.normals[i] = get_velocity_L(valley, band_index, mu, tri.vertexes[i]);
    }

    size = tri.faces.size();
    for (int i=0; i<size; i++) {
        int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
        vector3 k1 = tri.vertexes[v[0]];
        vector3 k2 = tri.vertexes[v[1]];
        vector3 k3 = tri.vertexes[v[2]];
        vector3 n1 = tri.normals[v[0]];
        vector3 n2 = tri.normals[v[1]];
        vector3 n3 = tri.normals[v[2]];
        for(int axis=0; axis<space_dim; axis++) {
            tri.faces[i].center[axis] = k1.vec[axis];
            tri.faces[i].normal[axis] = n1.vec[axis];
        }
        tri.faces[i].dS = get_dS(k1, k2, k3);
    }

    return tri;

//    for (int i=0; i<size; i++) {
//        double ene = get_E_L(valley, band_index, tri.vertexes[i]) - mu;
//        if (std::abs(ene) > eps_phys) {
////            std::cout << std::scientific << "#" << i << ", 誤差が大きい: e = " << ene << std::endl;
//
//            velocity v = get_velocity_L(valley, band_index, mu, tri.vertexes[i]);
//            double norm = 0e0;
//            for(int axis=0; axis<space_dim; axis++) {
//                norm += v.vec[axis]*v.vec[axis];
//            }
//            norm = NRsqrt(norm);
//
//            double min;
//            int f[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
//            vector3 k1 = tri.vertexes[f[0]];
//            vector3 k2 = tri.vertexes[f[1]];
//            vector3 k3 = tri.vertexes[f[2]];
//            min = get_min_length(k1, k2, k3) * 5e-2;
//
//            kpoint k;
//            double ek;
//            double pm = ene / std::abs(ene);
//            for(int i=1; i<50; i++) {
//                for(int axis=0; axis<space_dim; axis++) {
//                    k.vec[axis] = tri.vertexes[i].vec[axis] - pm * v.vec[axis]/norm*min*double(i);
//                }
//                ek = get_E_L(valley, band_index, k) - mu;
//                if (ek*ene < 0e0) break;
//            }
//            if ( ek*ene < 0e0 ) {
////                std::cout << std::scientific << "ek = " << ek << ", ene = " << ene << std::endl;
//                kpoint q = bisec_L2points(valley, band_index, mu, ek, ene, k, tri.vertexes[i]);
//                ek = get_E_L(valley, band_index, q) - mu;
////                std::cout << std::scientific << "#" << i << ", 再構成OK: ek = " << ek << std::endl;
//                tri.vertexes[i] = q;
//
//            }
//        }
//
//        tri.normals[i] = get_velocity_L(valley, band_index, mu, tri.vertexes[i]);
//    }
//
//    size = tri.faces.size();
//    for (int i=0; i<size; i++) {
//        int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
//        vector3 k1 = tri.vertexes[v[0]];
//        vector3 k2 = tri.vertexes[v[1]];
//        vector3 k3 = tri.vertexes[v[2]];
//        vector3 n1 = tri.normals[v[0]];
//        vector3 n2 = tri.normals[v[1]];
//        vector3 n3 = tri.normals[v[2]];
//
////        kpoint center = get_center(k1, k2, k3);
//        kpoint center = k1;
////        vector3 normal = get_normal(n1, n2, n3);
//        vector3 normal = n1;
//        for(int axis=0; axis<space_dim; axis++) {
//            tri.faces[i].center[axis] = center.vec[axis];
//            tri.faces[i].normal[axis] = normal.vec[axis];
//        }
//        tri.faces[i].dS = get_dS(k1, k2, k3);
//
////        double n = 0e0;
////        for(int axis=0; axis<space_dim; axis++) {
////            k1.vec[axis] = k1.vec[axis] - k3.vec[axis];
////            k2.vec[axis] = k2.vec[axis] - k3.vec[axis];
////        }
////        k3.vec[0] = k1.vec[1]*k2.vec[2] - k1.vec[2]*k2.vec[1];
////        k3.vec[1] = k1.vec[2]*k2.vec[0] - k1.vec[0]*k2.vec[2];
////        k3.vec[2] = k1.vec[0]*k2.vec[1] - k1.vec[1]*k2.vec[0];
////        for(int axis=0; axis<space_dim; axis++) {
////            n += k3.vec[axis]*k3.vec[axis];
////        }
////        n = NRsqrt(n)*5e-1;
//
//    }
//
////    int count = 0;
////    for (int i=0;i<2;i++) {
//////    do  {
////        count = 0;
////        int size = tri.faces.size();
////        for (int i=0; i<size; i++) {
////            if (tri.faces[i].dS > 1e-5) {
//////                std::cout << i << ", dS = " << tri.faces[i].dS << std::endl;
////                count++;
////                divide_triangle(tri, valley, band_index, mu, i);
//////                std::cout << i << ", " << tri.faces[i].dS << ", " <<tri.vertexes.size() << std::endl;
////            }
//////            if (count == 10) exit(0);
////        }
////        std::cout << "count = " << count << std::endl;
//////    } while (count > 0);
////    };
//////    std::string check_name = "check.csv";
//////    std::ofstream ofs(check_name);
//////    size = tri.faces.size();
//////    for (int i=0; i<size; i++) {
//////        if (tri.faces[i].dS > 1e-4) {
//////        ofs << std::scientific
//////            << tri.vertexes[tri.faces[i].face[0]].vec[0] << ", "
//////            << tri.vertexes[tri.faces[i].face[0]].vec[1] << ", "
//////            << tri.vertexes[tri.faces[i].face[0]].vec[2]
//////            << std::endl;
//////        ofs << std::scientific
//////            << tri.vertexes[tri.faces[i].face[1]].vec[0] << ", "
//////            << tri.vertexes[tri.faces[i].face[1]].vec[1] << ", "
//////            << tri.vertexes[tri.faces[i].face[1]].vec[2]
//////            << std::endl;
//////        ofs << std::scientific
//////            << tri.vertexes[tri.faces[i].face[2]].vec[0] << ", "
//////            << tri.vertexes[tri.faces[i].face[2]].vec[1] << ", "
//////            << tri.vertexes[tri.faces[i].face[2]].vec[2]
//////            << std::endl;
//////        ofs << std::scientific
//////            << tri.vertexes[tri.faces[i].face[0]].vec[0] << ", "
//////            << tri.vertexes[tri.faces[i].face[0]].vec[1] << ", "
//////            << tri.vertexes[tri.faces[i].face[0]].vec[2]
//////            << std::endl;
//////        ofs << ",," << std::endl;
//////        }
//////    }

    return tri;
}; // }}}

double get_DOS_L(triangles tri, int valley, int band_index, chemical_potential mu) { // {{{
    double dos = 0e0;
    auto fn = [](int valley, int band_index, chemical_potential mu, kpoint k) { return 1e0; };
    integrate_triangles_L(fn, dos, tri, valley, band_index, mu);

    return dos;
}; // }}}

int triangles_write_L(triangles tri, std::string name, int valley) { // {{{
    int size = tri.vertexes.size();
    std::string vertex_name = name+"_vertex.csv";
    if (size > 0) {
        std::ofstream ofs(vertex_name);
        if (!ofs) {
            std::cout << vertex_name << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            ofs << std::scientific
                << kL[valley][0]+tri.vertexes[i].vec[0] << ", "
                << kL[valley][1]+tri.vertexes[i].vec[1] << ", "
                << kL[valley][2]+tri.vertexes[i].vec[2] << ", "
                << tri.normals[i].vec[0] << ", " << tri.normals[i].vec[1] << ", " << tri.normals[i].vec[2] << ", "
                << std::endl;
        }
    }

    size = tri.faces.size();
    std::string face_name = name+"_face.csv";
    if (size > 0) {
        std::ofstream ofs(face_name);
        if (!ofs) {
            std::cout << face_name << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};

//            double norm = 0e0;
//            for(int axis=0; axis<space_dim; axis++) {
//                norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
//            }
//            norm = NRsqrt(norm);

            ofs << std::scientific
                << kL[valley][0]+tri.vertexes[v[0]].vec[0] << ", "
                << kL[valley][1]+tri.vertexes[v[0]].vec[1] << ", "
                << kL[valley][2]+tri.vertexes[v[0]].vec[2]
                << std::endl;
            ofs << std::scientific
                << kL[valley][0]+tri.vertexes[v[1]].vec[0] << ", "
                << kL[valley][1]+tri.vertexes[v[1]].vec[1] << ", "
                << kL[valley][2]+tri.vertexes[v[1]].vec[2]
                << std::endl;
            ofs << std::scientific
                << kL[valley][0]+tri.vertexes[v[2]].vec[0] << ", "
                << kL[valley][1]+tri.vertexes[v[2]].vec[1] << ", "
                << kL[valley][2]+tri.vertexes[v[2]].vec[2]
                << std::endl;
            ofs << std::scientific
                << kL[valley][0]+tri.vertexes[v[0]].vec[0] << ", "
                << kL[valley][1]+tri.vertexes[v[0]].vec[1] << ", "
                << kL[valley][2]+tri.vertexes[v[0]].vec[2]
                << std::endl;
            ofs << ",," << std::endl;

        }

    }
    std::string check_name = name+"_check.csv";
    std::ofstream oft(check_name);
    for (int i=0; i<size; i++) {
        oft << std::scientific << i << ", " << tri.faces[i].dS << std::endl;
    }

    return 0;
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
    int size = tri.faces.size();
    if (size > 0) {
        std::ofstream ofs(filename);
        if (!ofs) {
            std::cout << filename << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};

            double norm = 0e0;
            for(int axis=0; axis<space_dim; axis++) {
                norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
            }
            norm = NRsqrt(norm);

            ofs << std::scientific
                << tri.vertexes[v[0]].vec[0] << ", " << tri.vertexes[v[0]].vec[1] << ", " << tri.vertexes[v[0]].vec[2]
                << std::endl;
            ofs << std::scientific
                << tri.vertexes[v[1]].vec[0] << ", " << tri.vertexes[v[1]].vec[1] << ", " << tri.vertexes[v[1]].vec[2]
                << std::endl;
            ofs << std::scientific
                << tri.vertexes[v[2]].vec[0] << ", " << tri.vertexes[v[2]].vec[1] << ", " << tri.vertexes[v[2]].vec[2]
                << std::endl;
            ofs << std::scientific
                << tri.vertexes[v[0]].vec[0] << ", " << tri.vertexes[v[0]].vec[1] << ", " << tri.vertexes[v[0]].vec[2]
                << std::endl;
            ofs << std::endl;
        }
    }

    return 0;
}; // }}}
