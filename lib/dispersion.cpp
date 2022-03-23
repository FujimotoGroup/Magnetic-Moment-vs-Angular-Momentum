#include "parameters.hpp"

int band_index = 0;
int valley_index = 0;

void dispersionT(gdouble ** f, GtsCartesianGrid g, guint k, gpointer data) {
    initialize();

    gdouble x, y, z = g.z;
    guint i, j;

    for (i = 0, x = g.x; i < g.nx; i++, x += g.dx) {
        for (j = 0, y = g.y; j < g.ny; j++, y += g.dy) {
            vectorReal momentum = {x, y, z};
            matrixComplex H_T = set_T(momentum);
            vectorReal E_T = diagonalize_N(H_T);
            f[i][j] = E_T[band_index];
//            f[i][j] = x*x + y*y + z*z;
        }
    }
}

void dispersionL(gdouble ** f, GtsCartesianGrid g, guint k, gpointer data) {
    initialize();

    gdouble x, y, z = g.z;
    guint i, j;

    for (i = 0, x = g.x; i < g.nx; i++, x += g.dx) {
        for (j = 0, y = g.y; j < g.ny; j++, y += g.dy) {
            vectorReal momentum = {x, y, z};
            matrixComplex H_L = set_L(valley_index, momentum);
            vectorReal E_L = diagonalize_N(H_L);
            f[i][j] = E_L[band_index];
//            f[i][j] = x*x + y*y + z*z;
        }
    }
}

