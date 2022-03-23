#include "parameters.hpp"

static void write_vertex (GtsVertex * v, guint * nv) {
  printf ("  GtsVertex * v%u = gts_vertex_new (gts_vertex_class (), %g, %g, %g);\n", *nv, GTS_POINT (v)->x, GTS_POINT (v)->y, GTS_POINT (v)->z);
  GTS_OBJECT (v)->reserved = GUINT_TO_POINTER ((*nv)++);
};

static void write_face (GtsTriangle * t, guint * nf)
{
    printf ("  GtsFace * f%u = gts_face_new (gts_face_class (),\n"
    "                                           %g, %g, %g);\n",
    (*nf)++,
    GTS_POINT (GTS_OBJECT (t->e1->segment.v1))->x,
    GTS_POINT (GTS_OBJECT (t->e2->segment.v1))->x,
    GTS_POINT (GTS_OBJECT (t->e3->segment.v1))->x);
};

int main(){
    initialize();

// output kx-dependence @ T {{{
    int n = 50;
    std::string file_name = "./dat/EigenValue_T_"+std::to_string(bandsT)+"bands.dat";
    std::ofstream file(file_name);
    for (int i=-n; i<=n; i++) {
        double kx = cutoff*double(i)/double(2*n);
        vectorReal k = {kx,0e0,0e0};
        matrixComplex H_T = set_T(k);

        vectorReal E_T = diagonalize_N(H_T);

        file << std::scientific << kx << ", ";
        for(int i=0; i<bandsT; i++) {
            file << std::scientific << E_T[i];
            if (i != bandsT-1) file << ", ";
        }
        file << std::endl;

    };
/// }}}
// output kx-dependence @ Ls {{{
    for(int valley=0; valley<valleys; valley++) {
        std::string file_name = "./dat/EigenValue_L"+std::to_string(valley+1)+"_"+std::to_string(bandsL)+"bands.dat";
        std::ofstream file(file_name);
        for (int i=-n; i<=n; i++) {
            double kx = cutoff*double(i)/double(2*n);
            vectorReal k = {kx,0e0,0e0};
            matrixComplex H_L = set_L(valley, k);

            vectorReal E_L = diagonalize_N(H_L);

            file << std::scientific << kx << ", ";
            for(int i=0; i<bandsL; i++) {
                file << std::scientific << E_L[i];
                if (i != bandsL-1) file << ", ";
            }
            file << std::endl;

        };
    };
/// }}}

    GtsCartesianGrid g;
    g.nx = 10;
    g.ny = 10;
    g.nz = 10;
    /* interval is [-10:10][-10:10][-10:10] */
    g.x = -cutoff; g.dx = 2e0*cutoff/(gdouble) (g.nx - 1);
    g.y = -cutoff; g.dy = 2e0*cutoff/(gdouble) (g.ny - 1);
    g.z = -cutoff; g.dz = 2e0*cutoff/(gdouble) (g.nz - 1);

    GtsSurface * surface;
    gdouble iso;
    GtsIsoCartesianFunc func = dispersionT;

    iso = 0e0;

    surface = gts_surface_new (gts_surface_class (),
                               gts_face_class (),
                               gts_edge_class (),
                               gts_vertex_class ());
    band_index = 5;
    gts_isosurface_cartesian (surface, g, func, NULL, iso);
    std::cout << surface->faces << std::endl;

//    gts_surface_print_stats (surface, stderr);

    guint nv = 0, nf = 0;
//    gts_surface_foreach_vertex (surface, (GtsFunc) write_vertex, &nv);
    gts_surface_foreach_face (surface, (GtsFunc) write_face, &nf);

    return 0;
}
