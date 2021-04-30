#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xdrfile/xdrfile_xtc.h"

#define rPI 3.1416

// prototypes
void read_params(int argc, char *argv[]);
void alloc_mem(void);
void free_mem(void);
void analyze_MSDs(int mol_ind);
int map(int nb);
int invmap(int nb);
double vect_norm(double *a);
double vect_scalar_prod(double *a, double *b);
int vect_cross_prod(double *a, double *b, double *c);

// global variables
char fname[300], fname_out[300];
int nb_atoms, nb_mols, nb_frames;
int delta_frame=1, nb_mols_fr=47;
int nb_tau;
long int *N;
float *DATA, *Ltab;
double L, timestep, delta_z, e10=0.02;
double *traj, *traj_loc, *traj_th;
double *loc_frame;
double *MSD, *MSD_loc1, *MSD_loc2, *MSD_loc3;
double *MSD_th1, *MSD_th2, *MSD_th3;
rvec *X;
matrix box;
XDRFILE *traj_xtc;
FILE *fp_out;


/****************************************************************************/


int main(int argc, char **argv) {
  int loop, index, i, j, k, mol_offset;
  int status, step;
  float time, prec;
  char fname_end[300];

  // read in args: fname, delta_z, no. of frames
  read_params(argc, argv);
  // store no. of atoms in nb_atoms
  read_xtc_natoms(fname, &nb_atoms);
  // one HOH-molecule has 3 atoms
  nb_mols = nb_atoms / 3;
  printf("# %d atoms corresponding to %d molecules\n", nb_atoms, nb_mols);

  // allocate memory
  alloc_mem();

  // loop over batch of molecules and calculate MSD
  for (mol_offset = 0; mol_offset + nb_mols_fr <= nb_mols; mol_offset += nb_mols_fr) {
    printf("# Reading positions of molecules %d-%d", mol_offset, mol_offset + nb_mols_fr - 1);
    fflush(stdout);

    // allocate memory for positions
    X = calloc(nb_atoms, sizeof(X[0]));
    // open trajectory file
    traj_xtc = xdrfile_open(fname, "r");
    // read first frame
    status = read_xtc(traj_xtc, nb_atoms, &step, &time, box, X, &prec);
    // printf("# Box: %8.3f %8.3f %8.3f\n", box[0][0], box[1][1], box[2][2]);
    // printf("# Precision: %8.3f\n", prec);

    // loop over time-frames for batch of nb_mols_fr molecules
    loop = 0;
    index = 0;
    do {
      // store timestep from first frame
      if (loop == 1) {
        timestep = time;
      }
      // abort when no. of frames reached
      if (loop >= nb_frames) {
        break;
      }
      // read in every delta_frame frame
      if (loop % delta_frame == 0) {
        // tabulate box length
        Ltab[index] = box[0][0];
        // i = H2O molecule #1 .. #nb_mols_fr
        for (i = 0; i < nb_mols_fr; i++) {
          // j = O, H1, H2
          for (j = 0; j < 3; j++) {
            // k = x, y, z
            for (k = 0; k < 3; k++) {
              // flattened array of positions for every time-frame (index)
              DATA[3 * (3*index*nb_mols_fr + 3*i + j) + k] = X[3 * (mol_offset + i) + j][k];
            }
          }
        }
        index++;
      }
      loop++;
    } while(!read_xtc(traj_xtc, nb_atoms, &step, &time, box, X, &prec));
    // close trajectory and free memory
    xdrfile_close(traj_xtc);
    free(X);

    printf(" - analyzing");
    fflush(stdout);
    // i = no. of molecule in current batch
    for (i = 0; i < nb_mols_fr; i++) {
      analyze_MSDs(i);
    }
    printf(" - done\n");
    fflush(stdout);
  }

  printf("# %d frames read, %d frames analyzed, timestep = %3.15g ps\n", loop, index, timestep);
  sprintf(fname_end, "_dz%2.2f_%dfr_%2.3fps.txt", delta_z, nb_frames, timestep);
  sprintf(fname_out, "trans_MSDs%s", fname_end);
  printf("# Printing %s ", fname_out);

  // write results to file
  fp_out = fopen(fname_out, "w");
  for (i = 0; i < nb_tau; i++) {
    if (N[i] != 0) {
      fprintf(fp_out, "%3.3g %3.8g %3.8g %3.8g %3.8g\n", invmap(i) * timestep, MSD[i] / N[i], MSD_loc1[i] / N[i], MSD_loc2[i] / N[i], MSD_loc3[i] / N[i]);
    }
  }
  fclose(fp_out);
  printf("- done\n");
  fflush(stdout);

  sprintf(fname_out, "rot_MSDs%s", fname_end);
  printf("# Printing %s ", fname_out);
  fp_out = fopen(fname_out, "w");
  for (i = 0; i < nb_tau; i++) {
    if (N[i] != 0) {
      fprintf(fp_out, "%3.3g %3.8g %3.8g %3.8g\n", invmap(i) * timestep, MSD_th1[i] / N[i], MSD_th2[i] / N[i], MSD_th3[i] / N[i]);
    }
  }
  fclose(fp_out);
  printf("- done\n");
  fflush(stdout);

  // free memory and exit
  free_mem();
  return 0;
}


/****************************************************************************/


void read_params(int argc, char* argv[]) {
  switch (argc) {
    case 4:
      sscanf(argv[3], "%d", &nb_frames);
      sscanf(argv[2], "%lf", &delta_z);
      sscanf(argv[1], "%s", fname);
      break;
    default: {
      fprintf(stderr, "%s\nUsage: filename, delta_z, no. of frames\n", argv[0]);
      exit(0);
    }
  }
}


void alloc_mem(void) {
  int tmp, i;

  nb_tau = map(nb_frames);

  // tabulated box-length L for every time-frame
  Ltab = calloc(nb_frames, sizeof(double));

  // flattened array of positions for every time-frame:
  // no. of frames * no. of molecules curr. batch * 3 atoms per mol. * 3 dim.
  tmp = 3 * nb_frames * nb_mols_fr;
  DATA = calloc(3 * tmp, sizeof(float));

  traj = calloc(tmp, sizeof(double));
  traj_loc = calloc(tmp, sizeof(double));
  traj_th = calloc(tmp, sizeof(double));
  loc_frame = calloc(3 * tmp, sizeof(double));

  MSD = calloc(nb_tau, sizeof(double));
  MSD_loc1 = calloc(nb_tau, sizeof(double));
  MSD_loc2 = calloc(nb_tau, sizeof(double));
  MSD_loc3 = calloc(nb_tau, sizeof(double));
  MSD_th1 = calloc(nb_tau, sizeof(double));
  MSD_th2 = calloc(nb_tau, sizeof(double));
  MSD_th3 = calloc(nb_tau, sizeof(double));

  N = calloc(nb_tau, sizeof(long int));

  printf("# Memory for %d frames allocated\n", nb_frames);
}


void free_mem(void) {
  free(Ltab);
  free(DATA);
  free(traj);
  free(traj_loc);
  free(traj_th);
  free(loc_frame);
  free(MSD);
  free(MSD_loc1);
  free(MSD_loc2);
  free(MSD_loc3);
  free(MSD_th1);
  free(MSD_th2);
  free(MSD_th3);
  free(N);
  printf("# Memory free\n");
}


void analyze_MSDs(int mol_ind) {
  int i, j, O_ind, T, dt, lastdt=-1;
  double dtmp;
  double OH1[3], OH2[3], X_h2o[3], Y_h2o[3], Z_h2o[3], delta_R[3];
  double X_h2o_last[3], Y_h2o_last[3], Z_h2o_last[3];

  // initialize
  for (j = 0; j < 3; j++) {
    traj[j] = 0;
    traj_loc[j] = 0;
    traj_th[j] = 0;
  }

  // loop over time-frames
  for (i = 0; i < nb_frames; i++) {
    // box-length in current time-frame
    L = Ltab[i];
    // index of O atom #mol_ind for all time-frames i in DATA / 3
    O_ind = 3 * (nb_mols_fr * i + mol_ind);
    // Eckart frame unit vectors Z, X
    // j = x, y, z
    for (j = 0; j < 3; j++) {
      // OH vector corrected for PBC if distance > L/2
      OH1[j] = (double) fmod(DATA[3 * (O_ind + 1) + j] - DATA[3*O_ind + j] + 1.5*L, L) - 0.5*L;
      OH2[j] = (double) fmod(DATA[3 * (O_ind + 2) + j] - DATA[3*O_ind + j] + 1.5*L, L) - 0.5*L;
      Z_h2o[j] = OH1[j] + OH2[j];
      X_h2o[j] = OH1[j] - OH2[j];
    }
    // normalize X unit vector
    dtmp = vect_norm(X_h2o);
    for (j = 0; j < 3; j++) {
      X_h2o[j] = X_h2o[j] / dtmp;
    }
    // Y = Z x X
    vect_cross_prod(Z_h2o, X_h2o, Y_h2o);
    // normalize Y unit vector
    dtmp = vect_norm(Y_h2o);
    for (j = 0; j < 3; j++) {
      Y_h2o[j] = Y_h2o[j] / dtmp;
    }
    // recalculate Z unit vector and normalize
    vect_cross_prod(X_h2o, Y_h2o, Z_h2o);
    dtmp = vect_norm(Z_h2o);
    for (j = 0; j < 3; j++) {
      Z_h2o[j] = Z_h2o[j] / dtmp;
    }
    // store Eckart unit vectors in flat array
    for (j = 0; j < 3; j++) {
      loc_frame[9*i + j] = X_h2o[j];
      loc_frame[9*i + 3 + j] = Y_h2o[j];
      loc_frame[9*i + 6 + j] = Z_h2o[j];
    }

    // calculate MSD Deltas
    if (i > 0) {
      // lab-frame
      for (j = 0; j < 3; j++) {
        // correct PBC
        delta_R[j] = (double) fmod(DATA[3*O_ind + j] + delta_z * Z_h2o[j] - DATA[3 * (O_ind - 3*nb_mols_fr) + j] - delta_z * Z_h2o_last[j] + 1.5*L, L) - 0.5*L;
        traj[3*i + j] = traj[3 * (i - 1) + j] + delta_R[j];
      }
      // Eckart frame tranlations Delta alpha
      traj_loc[3*i] = traj_loc[3 * (i - 1)] + 0.5 * (vect_scalar_prod(delta_R, X_h2o_last) + vect_scalar_prod(delta_R, X_h2o));
      traj_loc[3*i + 1] = traj_loc[3 * (i - 1) + 1] + 0.5 * (vect_scalar_prod(delta_R, Y_h2o_last) + vect_scalar_prod(delta_R, Y_h2o));
      traj_loc[3*i + 2] = traj_loc[3 * (i - 1) + 2] + 0.5 * (vect_scalar_prod(delta_R, Z_h2o_last) + vect_scalar_prod(delta_R, Z_h2o));
      // Eckart frame rotations Delta phi
      traj_th[3*i] = traj_th[3 * (i - 1)] + 0.5 * (vect_scalar_prod(Y_h2o, Z_h2o_last) - vect_scalar_prod(Z_h2o, Y_h2o_last));
      traj_th[3*i + 1] = traj_th[3 * (i - 1) + 1] + 0.5 * (vect_scalar_prod(Z_h2o, X_h2o_last) - vect_scalar_prod(X_h2o, Z_h2o_last));
      traj_th[3*i + 2] = traj_th[3 * (i - 1) + 2] + 0.5 * (vect_scalar_prod(X_h2o, Y_h2o_last) - vect_scalar_prod(Y_h2o, X_h2o_last));
    }
    // store last frame's unit vectors (also if i == 0)
    for (j = 0; j < 3; j++) {
      X_h2o_last[j] = X_h2o[j];
      Y_h2o_last[j] = Y_h2o[j];
      Z_h2o_last[j] = Z_h2o[j];
    }
  }

  for (T = 0; T < nb_tau; T++) {
    dt = invmap(T);
    if (dt > lastdt) {
      N[T] += nb_frames - dt;
      for (i = 0; i < nb_frames - dt; i++) {
        dtmp = 0.0;
        for (j = 0; j < 3; j++) {
          dtmp += pow(traj[3 * (i + dt) + j] - traj[3*i +j], 2);
        }
        MSD[T] += dtmp;
        MSD_loc1[T] += pow(traj_loc[3 * (i + dt)] - traj_loc[3*i], 2);
        MSD_loc2[T] += pow(traj_loc[3 * (i + dt) + 1] - traj_loc[3*i + 1], 2);
        MSD_loc3[T] += pow(traj_loc[3 * (i + dt) + 2] - traj_loc[3*i + 2], 2);
        MSD_th1[T] += pow(traj_th[3 * (i + dt)] - traj_th[3*i], 2);
        MSD_th2[T] += pow(traj_th[3 * (i + dt) + 1] - traj_th[3*i + 1], 2);
        MSD_th3[T] += pow(traj_th[3 * (i + dt) + 2] - traj_th[3*i + 2], 2);
      }
    }
    lastdt = dt;
  }
}


int map(int nb) {
  return (int) ceil(log10(nb) / e10);
  // return nb;
}


int invmap(int nb) {
  return (int) floor(pow(10, e10 * nb));
  // return nb;
}


double vect_norm(double *a) {
  int k;
  double dtmp = 0;
  for (k = 0; k < 3; k++) {
    dtmp += pow(a[k], 2);
  }
  return sqrt(dtmp);
}


double vect_scalar_prod(double *a, double *b) {
  int k;
  double dtmp = 0;
  for (k = 0; k < 3; k++) {
    dtmp += a[k] * b[k];
  }
  return dtmp;
}


int vect_cross_prod(double *a, double *b, double *c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}
