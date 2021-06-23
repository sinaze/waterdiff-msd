#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "xdrfile_xtc.h"
#include "analyze.h"
#include "cblas.h"

void tau_avrg(const rvec *r_msd_tau, const rvec *alpha_msd_tau,
              const rvec *phi_msd_tau, const int nframes,
              long int *n_tau, float *r_msd, rvec *alpha_msd,
              rvec *phi_msd) {
  int tau, i, j;
  rvec delta_msd;

  for (tau = 0; tau < nframes; tau++) {
    n_tau[tau] += nframes - tau;
    for (i = 0; i < nframes - tau; i++) {
      rxmyz(r_msd_tau[i+tau], r_msd_tau[i], delta_msd);
      r_msd[tau] += cblas_sdot(DIM, delta_msd, INC, delta_msd, INC);
      rxmyz(alpha_msd_tau[i+tau], alpha_msd_tau[i], delta_msd);
      for (j = 0; j < DIM; j++) {
        alpha_msd[tau][j] += pow(delta_msd[j], 2);
      }
      rxmyz(phi_msd_tau[i+tau], phi_msd_tau[i], delta_msd);
      for (j = 0; j < DIM; j++) {
        phi_msd[tau][j] += pow(delta_msd[j], 2);
      }
    }
  }
}

void log_tau_avrg(const rvec *r_msd_tau, const rvec *alpha_msd_tau,
                  const rvec *phi_msd_tau, const int nframes,
                  const float a, long int *n_tau,
                  float *r_msd, rvec *alpha_msd, rvec *phi_msd) {
  int tau, i, j, ntau, t, t_prev;
  rvec delta_msd;

  ntau = logspace(nframes, a);
  for (tau = 0; tau < ntau; tau++) {
    t = ilogspace(t, a);
    if (t > t_prev) {
      n_tau[tau] += nframes - t;
      for (i = 0; i < nframes - t; i++) {
        rxmyz(r_msd_tau[i+t], r_msd_tau[i], delta_msd);
        r_msd[t] += cblas_sdot(DIM, delta_msd, INC, delta_msd, INC);
        rxmyz(alpha_msd_tau[i+t], alpha_msd_tau[i], delta_msd);
        for (j = 0; j < DIM; j++) {
          alpha_msd[tau][j] += pow(delta_msd[j], 2);
        }
        rxmyz(phi_msd_tau[i+t], phi_msd_tau[i], delta_msd);
        for (j = 0; j < DIM; j++) {
          phi_msd[tau][j] += pow(delta_msd[j], 2);
        }
      }
    t_prev = t;
    }
  }
}

int logspace(const int x, const float a) {
  return (int) ceil(log10(x) / a);
}

int ilogspace(const int x, const float a) {
  return (int) floor(pow(10, a * x));
}

void get_msd(const int mol, const rvec *curr, const rvec *prev,
             const float l, const float l_prev, const float delta_z,
             rvec delta_r, rvec delta_alpha, rvec delta_phi) {
  rvec x, y, z;
  rvec x_prev, y_prev, z_prev;
  int i_mol;

  i_mol = mol * DIM;
  get_eckart(curr[i_mol], curr[i_mol + 1], curr[i_mol + 2], l, x, y, z);
  get_eckart(prev[i_mol], prev[i_mol + 1], prev[i_mol + 2], l_prev,
             x_prev, y_prev, z_prev);
  get_delta(curr[i_mol], prev[i_mol], delta_z, x, y, z, x_prev, y_prev, z_prev, l,
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

void rxppyz(const rvec x, const rvec y, rvec z) {
  int i;

  for (i = 0; i < DIM; i++) {
      z[i] += x[i] + y[i];
  }
}

void rxppy(const rvec x, rvec y) {
  int i;

  for (i = 0; i < DIM; i++) {
      y[i] += x[i];
  }
}

void rxmyz(const rvec x, const rvec y, rvec z) {
  int i;

  for (i = 0; i < DIM; i++) {
      z[i] = x[i] - y[i];
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
