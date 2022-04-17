#include "parameters.hpp"
//

FT dispersionT(Point_3 p) {
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    kpoint k;
    k.vec[0] = double(p.x());
    k.vec[1] = double(p.y());
    k.vec[2] = double(p.z());
    const FT e = FT(get_E_T(band_index, k)) - mu;
//    std::cout << p.x() << ", " << p.y() << ", " << p.z() << ", " << e << std::endl;
    return e;
}

Surface_mesh get_triangles_cgal_T(CGAL::Surface_mesh_default_criteria_3<Tr> criteria) {
    double c = cutoff;
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    Surface_3 surface(dispersionT,             // pointer to function
                      Sphere_3(CGAL::ORIGIN, c), 1e-8);  // bounding sphere
    // meshing surface
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    Surface_mesh sm;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

    std::ofstream out("sphere.off");
    out << sm << std::endl;
    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";

    return sm;
}
