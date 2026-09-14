#include <stdio.h>

#include "structure.h"
#include "prototype.h"

double lai_mass(struct Pmas *mass, struct Pchar *pchar) { return mass->lai; }
double f_gpp(struct Grid *grid, struct Loct *loct, struct Pchar *pchar,
             struct Pmas *mass) { return 0.0; }
double f_rfm(struct Grid *grid, struct Loct *loct, struct Pchar *pchar,
             struct Pmas *mass) { return 0.0; }
double f_rcm(struct Grid *grid, struct Loct *loct, struct Pchar *pchar,
             struct Pmas *mass) { return 0.0; }
double f_rrm(struct Grid *grid, struct Loct *loct, struct Pchar *pchar,
             struct Pmas *mass) { return 0.0; }
double f_rfg(struct Grid *grid, struct Pchar *pchar, struct Pflx *flux) { return 0.0; }
double f_rcg(struct Grid *grid, struct Pchar *pchar, struct Pflx *flux) { return 0.0; }
double f_rrg(struct Grid *grid, struct Pchar *pchar, struct Pflx *flux) { return 0.0; }
void f_allocation(struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux) {
    flux->tpf = flux->tpc = flux->tpr = flux->tpg = flux->tpp = 0.0;
}
void reallocation_survival(struct Grid *grid, struct Pchar *pchar,
                           struct Pmas *mass, struct Pflx *flux) {}
void f_n_alloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar,
               struct Pmas *mass, struct Pflx *flux) {}
void f_n_realloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar,
                 struct Pmas *mass, struct Pflx *flux) {}
void f_n_abandon_salvage(struct Grid *grid, struct Loct *loct,
                         struct Pchar *pchar, struct Pmas *mass,
                         struct Pflx *flux) {}

int main(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    struct Pchar pchar = {0};
    struct Pmas mass = {0};
    struct Pflx flux = {0};

    if (scanf("%lf %lf %lf", &mass.fol, &mass.stm, &mass.rot) != 3) return 2;
    if (scanf("%lf %lf %lf %lf", &pchar.lf, &pchar.lc, &pchar.lr,
              &pchar.dcd) != 4) return 2;
    if (scanf("%ld %ld %lf %lf %lf %ld %ld",
              &pchar.season, &pchar.crop_stage, &mass.lai,
              &loct.pot_total_h, &loct.tmp_2m, &pchar.phototype,
              &grid.veg_type) != 7) return 2;

    plant_process(&grid, &loct, &flux, &pchar, &mass);
    printf("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
           mass.fol, mass.stm, mass.rot, flux.lf, flux.lc, flux.lr);
    return 0;
}