#include <math.h>
#include <stdio.h>
#include <string.h>

#include "structure.h"
#include "prototype.h"

static int run_pm(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    double evaporation, incep_tree, incep_c3, incep_c4;
    double transp_tree, transp_c3, transp_c4;

    if (scanf(
            "%lf %lf %lf %lf %lf %lf "
            "%lf %lf %lf %lf "
            "%lf %lf %lf "
            "%lf %lf",
            &loct.tmp_2m, &loct.air_prsr, &loct.vp, &loct.vpd,
            &loct.wnd_10m, &loct.daylen[0],
            &loct.rn_tree, &loct.rn_c3, &loct.rn_c4, &loct.rn_ground,
            &loct.gc_tree, &loct.gc_c3, &loct.gc_c4,
            &loct.soilwtr_l, &grid.fieldcap30) != 15) {
        return 2;
    }
    loct.doy = 0;
    loct.vps = f_vap_pre_sat(&grid, &loct);
    loct.slope_vps = f_slope_vps(&grid, &loct);
    loct.r_aero = f_r_aero(&grid, &loct);
    loct.air_dns = f_airdens(&grid, &loct);
    evaporation = pm_evap(&grid, &loct);
    incep_tree = pm_incep(&grid, &loct, 1);
    incep_c3 = pm_incep(&grid, &loct, 2);
    incep_c4 = pm_incep(&grid, &loct, 3);
    transp_tree = pm_transp(&grid, &loct, 1);
    transp_c3 = pm_transp(&grid, &loct, 2);
    transp_c4 = pm_transp(&grid, &loct, 3);

    printf(
        "%.17g,%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g\n",
        loct.vps, loct.slope_vps, loct.air_dns, loct.r_aero,
        evaporation,
        incep_tree, incep_c3, incep_c4,
        transp_tree, transp_c3, transp_c4, loct.rsst_soil
    );
    return 0;
}

static int run_radiation(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    struct Echar echar = {0};
    struct Mass mass = {0};
    double eK0_tree, eK0_c3, eK0_c4;
    double aaa, bbb, ccc;
    long hour;

    if (scanf(
            "%lf %lf %lf %lf %lf "
            "%lf %lf %lf "
            "%lf %lf %lf "
            "%lf %lf %lf %lf "
            "%lf %lf %lf",
            &loct.tmp_2m, &loct.tmp_sfc, &loct.vp, &loct.tcdc_clm,
            &loct.dswrf_sfc,
            &mass.tree.lai, &mass.c3.lai, &mass.c4.lai,
            &eK0_tree, &eK0_c3, &eK0_c4,
            &echar.tree.eK, &echar.c3.eK, &echar.c4.eK,
            &echar.tree.albedo,
            &echar.c3.albedo, &echar.c4.albedo, &echar.soil.albedo) != 18) {
        return 2;
    }
    if (scanf("%lf %lf", &loct.funder_c3, &loct.funder_c4) != 2) {
        return 2;
    }
    aaa = 1.0 - exp(-eK0_tree * mass.tree.lai);
    bbb = 1.0 - exp(-eK0_c3 * mass.c3.lai);
    ccc = 1.0 - exp(-eK0_c4 * mass.c4.lai);
    loct.fcover_tree = aaa;
    loct.fcover_c3 = (1.0 - aaa) * loct.funder_c3 * bbb;
    loct.fcover_c4 = (1.0 - aaa) * loct.funder_c4 * ccc;
    loct.fcover_ground = 1.0 - loct.fcover_tree - loct.fcover_c3 - loct.fcover_c4;
    for (hour = 0; hour < DSTEP; ++hour) {
        loct.sfcrad_h[hour] = loct.dswrf_sfc;
    }
    f_net_rad(&grid, &loct, &echar, &mass);
    printf(
        "%.17g,%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g\n",
        loct.rn_tree, loct.rn_c3, loct.rn_c4, loct.rn_ground, loct.rn_eco,
        loct.fcover_tree, loct.fcover_c3, loct.fcover_c4, loct.fcover_ground,
        loct.rn_short_tree, loct.rn_short_c3, loct.rn_short_c4, loct.rn_short_ground,
        loct.rn_long_tree, loct.rn_long_c3, loct.rn_long_c4, loct.rn_long_ground
    );
    return 0;
}

static int run_conductance(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    struct Pchar pchar = {0};
    struct Pmas mass = {0};
    if (scanf(
            "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &pchar.psat, &pchar.eK, &pchar.lue, &pchar.ppfd_t,
            &mass.lai, &loct.aCO2, &pchar.cmpcd, &loct.vpd,
            &pchar.gs_b0, &pchar.gs_b1, &pchar.gs_b2) != 11) {
        return 2;
    }
    printf("%.17g\n", f_canopy_cond(&grid, &loct, &pchar, &mass));
    return 0;
}

static int run_lai(void) {
    struct Pchar pchar = {0};
    struct Pmas mass = {0};
    if (scanf("%lf %lf", &mass.fol, &pchar.sla) != 2) return 2;
    printf("%.17g\n", lai_mass(&mass, &pchar));
    return 0;
}

static int run_extinction(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    struct Pchar pchar = {0};
    double solar_height;
    if (scanf("%lf %lf", &pchar.eK0, &solar_height) != 2) return 2;
    loct.doy = 0;
    loct.hour = 0;
    loct.solhgt_h[0][0] = solar_height;
    printf("%.17g\n", irr_attn(&grid, &loct, &pchar));
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    if (strcmp(argv[1], "pm") == 0) return run_pm();
    if (strcmp(argv[1], "radiation") == 0) return run_radiation();
    if (strcmp(argv[1], "conductance") == 0) return run_conductance();
    if (strcmp(argv[1], "lai") == 0) return run_lai();
    if (strcmp(argv[1], "extinction") == 0) return run_extinction();
    return 2;
}
