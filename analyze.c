#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "xdrfile_xtc.h"
#include "analyze.h"
#include "cblas.h"

int get_msd(const rvec *curr, const rvec *prev, const float l) {
  rvec x, y, z;

  get_eckart(curr[0], curr[1], curr[2], l, x, y, z);
  printf("\n");
  print_rvec(x, "x");
  print_rvec(y, "y");
  print_rvec(z, "z");

  return EXIT_SUCCESS;
}

void rxcrossyz(const rvec x, const rvec y, rvec z) {
  z[0] = x[1] * y[2] - x[2] * y[1];
  z[1] = x[2] * y[0] - x[0] * y[2];
  z[2] = x[0] * y[1] - x[1] * y[0];
}

void print_rvec(const rvec x, const char *name) {
  printf("%s = [ %f %f %f ]\n", name, x[0], x[1], x[2]);
}

void get_eckart(const rvec o, const rvec h1, const rvec h2, const float l,
                rvec x, rvec y, rvec z) {
  rvec oh1, oh2;
  int i;

  for (i = 0; i < DIM; i++) {
    oh1[i] = fmod(h1[i] - o[i] + 1.5 * l, l) - 0.5 * l;
    oh2[i] = fmod(h2[i] - o[i] + 1.5 * l, l) - 0.5 * l;
    z[i] = oh1[i] + oh2[i];
    x[i] = oh1[i] - oh2[i];
  }
  cblas_sscal(DIM, 1.0 / cblas_snrm2(DIM, x, INC), x, INC);
  rxcrossyz(z, x, y);
  cblas_sscal(DIM, 1.0 / cblas_snrm2(DIM, y, INC), y, INC);
  rxcrossyz(x, y, z);
  cblas_sscal(DIM, 1.0 / cblas_snrm2(DIM, z, INC), z, INC);
}
