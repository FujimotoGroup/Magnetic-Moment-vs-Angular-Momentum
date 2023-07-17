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

triangles set_triangles(chemical_potential mu, Surface_mesh mesh) { // {{{
    triangles tri;
    tri.ene = mu;

    for(Surface_mesh::Vertex_index vd : mesh.vertices()) {
        vector3 v = {mesh.point(vd).x(), mesh.point(vd).y(), mesh.point(vd).z()};
        tri.vertexes.push_back(v);
    }

    for(Surface_mesh::Face_index fd : mesh.faces()) {
        face f;
        int i = 0;
        for(Surface_mesh::Vertex_index vd : vertices_around_face(mesh.halfedge(fd), mesh)) {
            f.face[i] = vd;
            i++;
        }
        tri.faces.push_back(f);
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

Surface_mesh get_triangles_cgal_T(int band_index, chemical_potential mu, double bounce) { // {{{
    double c = cutoff;
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    auto dispersionT = [&](Point_3 p) {
        const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
        kpoint k;
        k.vec[0] = double(p.x());
        k.vec[1] = double(p.y());
        k.vec[2] = double(p.z());
        const FT e = FT(get_E_T(band_index, k)) - mu;
//        std::cout << p.x() << ", " << p.y() << ", " << p.z() << ", " << e << std::endl;
        return e;
    };
    Surface_3 surface(dispersionT,             // pointer to function
                      Sphere_3(CGAL::ORIGIN, c), 5e-9);  // bounding sphere

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                       bounce,  // radius bound
                                                       bounce); // distance bound

    // meshing surface
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    Surface_mesh sm;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
    for(Surface_mesh::Face_index fd : sm.faces()) {
        for(Surface_mesh::Vertex_index vd : vertices_around_face(sm.halfedge(fd), sm)) {
        kpoint k = {sm.point(vd).x(), sm.point(vd).y(), sm.point(vd).z()};
        double e = get_E_T(band_index, k) - mu;
        if (std::abs(e) > 1e-8)
            std::cout << vd << ", " << e  << std::endl;
        }
    }

//    std::ofstream out("sphere.off");
//    out << sm << std::endl;
//    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";

    return sm;
} // }}}

triangles get_triangles_T(int band_index, chemical_potential mu) { // {{{
    triangles tri;

    double sign = band_edge_T_sign[band_index];
    if(sign*mu <= sign*band_edge_T[band_index]) {
        return tri;
    }

    Surface_mesh mesh;
    double c = 1e-1;
    int size = 0;
    do {
        c *= 5e-1;
        mesh = get_triangles_cgal_T(band_index, mu, c);
        size = mesh.number_of_vertices();
    } while (size < fermi_surface_mesh_lim_T);
    size = mesh.number_of_faces();
    std::cout << "mu = " << mu << "; final fs face# " << size << " and vertex# " << mesh.number_of_vertices() << std::endl;

    tri = set_triangles(mu, mesh);

    size = tri.vertexes.size();
    tri.normals.resize(size);
    for (int i=0; i<size; i++) {
        tri.normals[i] = get_velocity_T(band_index, mu, tri.vertexes[i]);
    }

    size = tri.faces.size();
    for (int i=0; i<size; i++) {
        int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
        vector3 k1 = tri.vertexes[v[0]];
        vector3 k2 = tri.vertexes[v[1]];
        vector3 k3 = tri.vertexes[v[2]];
        vector3 center;
        vector3 normal;
        for(int axis=0; axis<space_dim; axis++) {
            center.vec[axis] = k1.vec[axis];
        }
        normal = get_velocity_T(band_index, mu, center);
        for(int axis=0; axis<space_dim; axis++) {
            tri.faces[i].center[axis] = center.vec[axis];
            tri.faces[i].normal[axis] = normal.vec[axis];
        }
        tri.faces[i].dS = get_dS(k1, k2, k3);
    }

    return tri;
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

//            double norm = 0e0;
//            for(int axis=0; axis<space_dim; axis++) {
//                norm += tri.faces[i].normal[axis] * tri.faces[i].normal[axis];
//            }
//            norm = NRsqrt(norm);

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
    std::string check_name = name+"_check.csv";
    std::ofstream oft(check_name);
    for (int i=0; i<size; i++) {
        oft << std::scientific << i << ", " << tri.faces[i].dS << std::endl;
    }

    return 0;
}; // }}}

int triangles_write_T(triangles tri, std::string name, vectorReal values) { // {{{
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
                << values[i]
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

    size = tri.faces.size();
    std::string tri_name = name+"_tri.csv";
    if (size > 0) {
        std::ofstream ofs(tri_name);
        if (!ofs) {
            std::cout << tri_name << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};

            ofs << v[0] << ", " << v[1] << ", " << v[2] << std::endl;
        }

    }
//    std::string check_name = name+"_check.csv";
//    std::ofstream oft(check_name);
//    for (int i=0; i<size; i++) {
//        oft << std::scientific << i << ", " << tri.faces[i].dS << std::endl;
//    }

    return 0;
}; // }}}
// }}}

// for L points {{{
double get_E_L(int valley, int band_index, kpoint k) { // {{{
    matrixComplex H_L = set_L(valley, k.vec);
    vectorReal E_L = diagonalize_N(H_L);
    return E_L[band_index];
}; // }}}

double get_E_L(int valley, int band_index, chemical_potential mu, kpoint k) { // {{{
    matrixComplex H_L = set_L(valley, mu, k.vec);
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
    double epsilon = 1e-8;
    for(int axis=0; axis<space_dim; axis++) {
        for(int i=0; i<2; i++) {
            double p = std::pow(-1e0, i);
            kpoint kp = k;
            kp.vec[axis] += epsilon*p;
            double ene = get_E_L(valley, band_index, mu, kp);
            v.vec[axis] += ene*p;
        }
        v.vec[axis] = v.vec[axis] / (2e0*epsilon);
    }

//    velocity v_old = {0e0, 0e0, 0e0};
//    for (int j=7; j<=8; j++) {
//        v = {0e0, 0e0, 0e0};
//        double epsilon = std::pow(1e-1,j);
//        for(int axis=0; axis<space_dim; axis++) {
//            std::cout << epsilon << ": ";
//            for(int i=0; i<2; i++) {
//                double p = std::pow(-1e0, i);
//                kpoint kp = k;
//                kp.vec[axis] += epsilon*p;
//                double ene = get_E_L(valley, band_index, mu, kp);
//                std::cout << ene*p << ", ";
//                v.vec[axis] += ene*p;
//            }
//            std::cout << v.vec[axis] << std::endl;
//            v.vec[axis] = v.vec[axis] / (2e0*epsilon);
//        }
//
//        if (j >= 8) {
//            for(int axis=0; axis<space_dim; axis++) {
//                double c = std::abs(v_old.vec[axis] - v.vec[axis]);
//                if (c > 1e-9) {
//                    std::cout << "error: " << epsilon << ", " << axis << ", " << c << ", v = " << v.vec[axis] << ", v_old = " << v_old.vec[axis] << std::endl;
//                }
//            }
//        }
//
//        v_old = v;
//    }
    return v;
}; // }}}

Surface_mesh get_triangles_cgal_L(int valley, int band_index, chemical_potential mu, double bounce) { // {{{
    double c = cutoff;
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    auto dispersionL = [&](Point_3 p) {
        const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
        kpoint k;
        k.vec[0] = double(p.x());
        k.vec[1] = double(p.y());
        k.vec[2] = double(p.z());
        const FT e = FT(get_E_L(valley, band_index, k)) - mu;
//        std::cout << p.x() << ", " << p.y() << ", " << p.z() << ", " << e << std::endl;
        return e;
    };
    Surface_3 surface(dispersionL,             // pointer to function
                      Sphere_3(CGAL::ORIGIN, c), 5e-9);  // bounding sphere

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
                                                       bounce,  // radius bound
                                                       bounce); // distance bound

    // meshing surface
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    Surface_mesh sm;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
    for(Surface_mesh::Face_index fd : sm.faces()) {
        for(Surface_mesh::Vertex_index vd : vertices_around_face(sm.halfedge(fd), sm)) {
        kpoint k = {sm.point(vd).x(), sm.point(vd).y(), sm.point(vd).z()};
        double e = get_E_L(valley, band_index, k) - mu;
        if (std::abs(e) > 1e-8)
            std::cout << vd << ", " << e << std::endl;
        }
    }

//    std::ofstream out("sphere.off");
//    out << sm << std::endl;
//    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";

    return sm;
} // }}}

triangles get_triangles_L(int valley, int band_index, chemical_potential mu) { // {{{
    triangles tri;

    double sign = band_edge_L_sign[valley][band_index];
    if(sign*mu <= sign*band_edge_L[valley][band_index]) {
        return tri;
    }

    Surface_mesh mesh;
    double c = 1e-1;
    int size = 0;
    do {
        c *= 5e-1;
        mesh = get_triangles_cgal_L(valley, band_index, mu, c);
        size = mesh.number_of_vertices();
    } while (size < fermi_surface_mesh_lim_L);
    size = mesh.number_of_faces();
    std::cout << "ene = " << mu << "; final fs face# " << size << " and vertex# " << mesh.number_of_vertices() << ", c = " << c << std::endl;

    tri = set_triangles(mu, mesh);

    size = tri.vertexes.size();
    tri.normals.resize(size);
    for (int i=0; i<size; i++) {
//        tri.normals[i] = {0e0, 0e0, 0e0};
        tri.normals[i] = get_velocity_L(valley, band_index, mu, tri.vertexes[i]);

    }

    size = tri.faces.size();
    for (int i=0; i<size; i++) {
        int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};
        vector3 k1 = tri.vertexes[v[0]];
        vector3 k2 = tri.vertexes[v[1]];
        vector3 k3 = tri.vertexes[v[2]];
        vector3 center;
        vector3 normal;
        for(int axis=0; axis<space_dim; axis++) {
            center.vec[axis] = k1.vec[axis];
        }
        normal = get_velocity_L(valley, band_index, mu, center);
        for(int axis=0; axis<space_dim; axis++) {
            tri.faces[i].center[axis] = center.vec[axis];
            tri.faces[i].normal[axis] = normal.vec[axis];
        }
        tri.faces[i].dS = get_dS(k1, k2, k3);

    }

//    std::string name = "./triangles_L"+std::to_string(valley+1)+"-mu"+std::to_string(mu);
//    triangles_write_L(tri, name, valley);

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

int triangles_write_L(triangles tri, std::string name, int valley, vectorReal values) { // {{{
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
                << values[i]
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

    size = tri.faces.size();
    std::string tri_name = name+"_tri.csv";
    if (size > 0) {
        std::ofstream ofs(tri_name);
        if (!ofs) {
            std::cout << tri_name << " cannot be opened." << std::endl;
            return 1;
        }

        for (int i=0; i<size; i++) {
            int v[3] = {tri.faces[i].face[0], tri.faces[i].face[1], tri.faces[i].face[2]};

            ofs << v[0] << ", " << v[1] << ", " << v[2] << std::endl;
        }

    }
//    std::string check_name = name+"_check.csv";
//    std::ofstream oft(check_name);
//    for (int i=0; i<size; i++) {
//        oft << std::scientific << i << ", " << tri.faces[i].dS << std::endl;
//    }

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

int fermi_surface_write(fermi_surface fs, std::string filename, vectorReal value) { // {{{
    int size = fs.kset.size();
    if (size != value.size()) {
        std::cout << value.size() << " does not match " << size << std::endl;
    }

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
                << ", " << value[i]
                << std::endl;
        }
    }

    return 0;
}; // }}}
