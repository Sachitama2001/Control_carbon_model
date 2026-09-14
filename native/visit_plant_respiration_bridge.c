#include <stdio.h>

#include "structure.h"
#include "prototype.h"

int main(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    struct Pchar pchar = {0};
    struct Pmas mass = {0};
    struct Pflx flux = {0};
    long positive_epp;

    if (scanf("%lf %lf %lf", &mass.fol, &mass.stm, &mass.rot) != 3) return 2;
    if (scanf("%lf %lf %lf", &pchar.rgf, &pchar.rgc, &pchar.rgr) != 3) return 2;
    if (scanf("%lf %lf %lf %lf %lf", &pchar.rmf0, &pchar.rmc_s,
              &pchar.rmc_h, &pchar.rmr_s, &pchar.rmr_h) != 5) return 2;
    if (scanf("%lf %lf %lf %lf %lf", &pchar.qTf0, &pchar.qTc0,
              &pchar.qTr0, &pchar.f_sz_s, &pchar.f_sz_r) != 5) return 2;
    if (scanf("%lf %lf", &loct.tmp_sfc, &loct.tmp10_soil) != 2) return 2;
    if (scanf("%lf %lf %lf %ld", &flux.tpf, &flux.tpc, &flux.tpr,
              &positive_epp) != 4) return 2;

    f_q10_ar(&loct, &pchar);
    f_spcfc_resp(&pchar, &mass);
    flux.rfm = f_rfm(&grid, &loct, &pchar, &mass);
    flux.rcm = f_rcm(&grid, &loct, &pchar, &mass);
    flux.rrm = f_rrm(&grid, &loct, &pchar, &mass);
    if (positive_epp) {
        flux.rfg = f_rfg(&grid, &pchar, &flux);
        flux.rcg = f_rcg(&grid, &pchar, &flux);
        flux.rrg = f_rrg(&grid, &pchar, &flux);
    }

    printf("%.17g,%.17g,%.17g,", pchar.qTf, pchar.qTc, pchar.qTr);
    printf("%.17g,%.17g,%.17g,", pchar.rmf, pchar.rmc, pchar.rmr);
    printf("%.17g,%.17g,%.17g,", flux.rfm, flux.rcm, flux.rrm);
    printf("%.17g,%.17g,%.17g,", flux.rfg, flux.rcg, flux.rrg);
    printf("%.17g,%.17g,%.17g,%.17g\n",
           mass.stm_sp, mass.rot_fn, mass.stm_ht, mass.rot_tp);
    return 0;
}