#include <stdio.h>
#include <string.h>

#include "structure.h"
#include "prototype.h"

/* Effective PM inputs isolate hydro_balance.c from hydro_flows.c. */
static double pm_values[7];

double pm_incep(struct Grid *grid, struct Loct *loct, long layer) {
    (void)grid;
    (void)loct;
    return pm_values[layer - 1];
}

double pm_evap(struct Grid *grid, struct Loct *loct) {
    (void)grid;
    (void)loct;
    return pm_values[3];
}

double pm_transp(struct Grid *grid, struct Loct *loct, long layer) {
    (void)grid;
    (void)loct;
    return pm_values[3 + layer];
}

int main(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    struct Echar echar = {0};
    struct Mass mass = {0};

    if (scanf(
            "%lf %lf %lf "
            "%lf %lf %lf "
            "%lf %lf %lf %lf %lf %lf %lf "
            "%lf %lf %lf "
            "%lf %lf %lf %lf %lf "
            "%31s",
            &mass.snwa, &mass.sw30, &mass.sww,
            &loct.prate_sfc, &loct.tmp_2m, &loct.tmp200_soil,
            &pm_values[0], &pm_values[1], &pm_values[2], &pm_values[3],
            &pm_values[4], &pm_values[5], &pm_values[6],
            &grid.fieldcap30, &grid.fieldcap, &grid.hyd_cond,
            &mass.tree.lai, &mass.c3.lai, &mass.c4.lai,
            &loct.funder_c3, &loct.funder_c4,
            grid.site_id) != 22) {
        return 2;
    }

    echar.tree.gc = 1.0;
    echar.c3.gc = 1.0;
    echar.c4.gc = 1.0;
    f_hydrology(&grid, &loct, &echar, &mass);

    printf(
        "%.17g,%.17g,%.17g,"
        "%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g\n",
        mass.snwa, mass.sw30, mass.sww,
        loct.snp, loct.thaw,
        loct.incep_tree, loct.incep_c3, loct.incep_c4, loct.incep,
        loct.ro1, loct.evpr,
        loct.trnsp_tree, loct.trnsp_c3, loct.trnsp_c4, loct.trnsp,
        loct.ro2, loct.aet, loct.pet
    );
    return 0;
}
