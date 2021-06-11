#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "xdrfile_xtc.h"
#include "analyze.h"
#include "cblas.h"

void get_msd(const rvec *curr, const rvec *prev, const float l,
             rvec delta_r, rvec delta_alpha, rvec delta_phi) {
  rvec x, y, z;
  rvec x_prev, y_prev, z_prev;

  get_eckart(curr[0], curr[1], curr[2], l, x, y, z);
  get_eckart(prev[0], prev[1], prev[2], l, x_prev, y_prev, z_prev);
  get_delta(curr[0], prev[0], 0.0, x, y, z, x_prev, y_prev, z_prev, l,
            delta_r, delta_alpha, delta_phi);
}

void rxcrossyz(const rvec x, const rvec y, rvec z) {
  z[0] = x[1] * y[2] - x[2] * y[1];
  z[1] = x[2] * y[0] - x[0] * y[2];
  z[2] = x[0] * y[1] - x[1] * y[0];
}

void rxpyz(const rvec x, const rvec y, rvec z) {
  int i;

  for (i = 0; i < DIM; i++) {
      z[i] = x[i] + y[i];
  }
}

void print_rvec(const rvec x, const char *name) {
  printf("%s = [ %f %f %f ]\n", name, x[0], x[1], x[2]);
}

void get_eckart(const rvec ow, const rvec h1, const rvec h2, const float l,
                rvec x, rvec y, rvec z) {
  rvec oh1, oh2;
  int i;

  for (i = 0; i < DIM; i++) {
    oh1[i] = pbc_corr(h1[i] - ow[i], l);
    oh2[i] = pbc_corr(h2[i] - ow[i], l);
    z[i] = oh1[i] + oh2[i];
    x[i] = oh1[i] - oh2[i];
  }
  cblas_sscal(DIM, 1.0 / cblas_snrm2(DIM, x, INC), x, INC);
  rxcrossyz(z, x, y);
  cblas_sscal(DIM, 1.0 / cblas_snrm2(DIM, y, INC), y, INC);
  rxcrossyz(x, y, z);
  cblas_sscal(DIM, 1.0 / cblas_snrm2(DIM, z, INC), z, INC);
}

/* TODO: pass-by-reference better? */
float pbc_corr(const float x, const float l) {
  return fmod(x + 1.5 * l, l) - 0.5 * l;
}

void get_delta(const rvec ow, const rvec ow_prev, const float delta_z,
               const rvec x, const rvec y, const rvec z,
               const rvec x_prev, const rvec y_prev, const rvec z_prev,
               const float l,
               rvec delta_r, rvec delta_alpha, rvec delta_phi) {
  int i;

  for (i = 0; i < DIM; i++) {
    delta_r[i] = pbc_corr(ow[i] + delta_z * z[i]
                          - ow_prev[i] + delta_z * z_prev[i], l);
  }
  delta_alpha[0] = 0.5 * (cblas_sdot(DIM, delta_r, INC, x, INC)
                          + cblas_sdot(DIM, delta_r, INC, x_prev, INC));
  delta_alpha[1] = 0.5 * (cblas_sdot(DIM, delta_r, INC, y, INC)
                          + cblas_sdot(DIM, delta_r, INC, y_prev, INC));
  delta_alpha[2] = 0.5 * (cblas_sdot(DIM, delta_r, INC, z, INC)
                          + cblas_sdot(DIM, delta_r, INC, z_prev, INC));
  delta_phi[0] = 0.5 * (cblas_sdot(DIM, y, INC, z_prev, INC)
                        - cblas_sdot(DIM, z, INC, y_prev, INC));
  delta_phi[1] = 0.5 * (cblas_sdot(DIM, z, INC, x_prev, INC)
                        - cblas_sdot(DIM, x, INC, z_prev, INC));
  delta_phi[2] = 0.5 * (cblas_sdot(DIM, x, INC, y_prev, INC)
                        - cblas_sdot(DIM, y, INC, x_prev, INC));
}
