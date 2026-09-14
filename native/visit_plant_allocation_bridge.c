#include <stdio.h>

#include "structure.h"
#include "prototype.h"

int main(void) {
    struct Pchar pchar = {0};
    struct Pmas mass = {0};
    struct Pflx flux = {0};

    if (scanf("%lf %lf %lf", &pchar.alloc_ass, &pchar.alloc_abg,
              &pchar.sla) != 3) return 2;
    if (scanf("%lf %lf %ld %ld", &mass.lai, &pchar.opt_lai,
              &pchar.season, &pchar.crop_stage) != 4) return 2;
    if (scanf("%lf %lf %lf %lf %lf", &flux.epp, &flux.gpp,
              &flux.rfm, &flux.rcm, &flux.rrm) != 5) return 2;

    f_allocation(&pchar, &mass, &flux);
    printf("%.17g,%.17g,%.17g,%.17g,", pchar.malloc_f, pchar.malloc_c,
           pchar.malloc_r, pchar.malloc_g);
    printf("%.17g,%.17g,%.17g,%.17g,%.17g\n",
           flux.tpf, flux.tpc, flux.tpr, flux.tpg, flux.tpp);
    return 0;
}