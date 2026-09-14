#include <stdio.h>
#include <stdlib.h>

#include "structure.h"
#include "prototype.h"

short SA_PARA = 0;
short SA_PARA_EN = 0;
double SA_PARA_VAR = 0.0;

void f_n_mineralz(struct Grid *grid, struct Loct *loct, struct Smas *mass,
                   struct Sflx *flux) {}
void f_nh3_volatilization(struct Grid *grid, struct Loct *loct,
                          struct Smas *mass, struct Sflx *flux) {}
void f_n_leaching(struct Grid *grid, struct Loct *loct, struct Smas *mass,
                  struct Sflx *flux) {}
void f_n_immoblz(struct Grid *grid, struct Loct *loct, struct Schar *schar,
                 struct Smas *mass, struct Sflx *flux) {}
void f_n_mcrb_abdn(struct Grid *grid, struct Loct *loct, struct Schar *schar,
                   struct Smas *mass, struct Sflx *flux) {}

static int read_double(double *value) {
    return scanf("%lf", value) == 1;
}

int main(void) {
    struct Grid grid = {0};
    struct Loct loct = {0};
    struct Schar schar = {0};
    struct Smas mass = {0};
    struct Sflx flux = {0};
    long n_steps;
    long step;

    double *states[] = {
        &mass.ltr_tf, &mass.ltr_tc, &mass.ltr_tr,
        &mass.ltr_gf, &mass.ltr_gc, &mass.ltr_gr,
        &mass.msl_a, &mass.msl_i, &mass.msl_p,
    };
    double *parameters[] = {
        &schar.sr_lf, &schar.sr_lc, &schar.sr_lr,
        &schar.sr_ha, &schar.sr_hi, &schar.sr_hp,
        &schar.f_co2_lf, &schar.f_co2_lc, &schar.f_co2_lr,
        &schar.f_hm_a, &schar.f_hm_i, &schar.f_hm_p,
        &schar.kml, &schar.kmh, &schar.kmsl, &schar.kmsh,
    };
    double *environment[] = {
        &loct.tmp10_soil, &loct.tmp200_soil,
        &loct.soilwtr_l, &loct.soilwtr_h,
        &loct.soilappr_l, &loct.soilappr_w,
        &grid.fieldcap30, &grid.fieldcap,
    };
    double *inputs[] = {
        &flux.li_tf, &flux.li_tc, &flux.li_tr,
        &flux.li_gf, &flux.li_gc, &flux.li_gr,
    };

    for (size_t index = 0; index < 9; ++index) {
        if (!read_double(states[index])) return 2;
    }
    for (size_t index = 0; index < 16; ++index) {
        if (!read_double(parameters[index])) return 2;
    }
    if (scanf("%ld", &n_steps) != 1 || n_steps < 0) return 2;

    for (step = 0; step < n_steps; ++step) {
        for (size_t index = 0; index < 6; ++index) {
            if (!read_double(inputs[index])) return 2;
        }
        for (size_t index = 0; index < 8; ++index) {
            if (!read_double(environment[index])) return 2;
        }
        f_cycle_soil(&grid, &loct, &schar, &mass, &flux);
        printf("%ld", step);
        for (size_t index = 0; index < 9; ++index) {
            printf(",%.17g", *states[index]);
        }
        printf(",%.17g,%.17g,%.17g\n",
               flux.hr, schar.f_tm_l, schar.f_tm_h);
    }
    return 0;
}