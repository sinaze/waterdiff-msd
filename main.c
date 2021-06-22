/* main.c */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <sys/types.h>

#include "xdrfile_xtc.h"
#include "xtc_seek.h"
#include "tools.h"
#include "analyze.h"

extern int errno;
extern char *optarg;
extern int opterr, optind;


int main(int argc, char *argv[]) {
    int opt;
    options_t options = { 0.0, 0, "\0" , 1, 0, "msd_out.txt" };

    int natoms, nmol;
    int mol;
    int molbatch = 53;

    XDRFILE *xd;
    int step;
    float time, prec;
    matrix box;
    rvec *x;

    int n_frames;
    int est_nframes;
    int64_t **offsets;

    int nframes;

    int i_frame = 0;
    float delta_t;
    rvec *curr;
    rvec *prev;

    rvec delta_r;
    rvec delta_alpha;
    rvec delta_phi;

    rvec *r_msd_tau;
    rvec *alpha_msd_tau;
    rvec *phi_msd_tau;

    long int *n_tau;
    float *r_msd;
    rvec *alpha_msd;
    rvec *phi_msd;

    FILE *fp;

    float l, l_prev;

    int i, i_mol;

    opterr = 0;

    while ((opt = getopt(argc, argv, OPTSTR)) != EOF)
       switch(opt) {
           case 'f':
              strcpy(options.fname, optarg);
              break;

           case 'z':
              options.delta_z = strtof(optarg, NULL);
              break;

           case 'n':
              options.max_frames = (uint64_t) strtoul(optarg, NULL, 10);
              n_frames = options.max_frames;
              break;

           case 's':
              options.stride = atoi(optarg);
              break;

           case 'i':
              options.i_mol = atoi(optarg);
              break;
           case 'o':
              strcpy(options.oname, optarg);
              break;

           case 'h':
           default:
              usage(basename(argv[0]), opt);
              /* NOTREACHED */
              break;
       }

    if (options.fname[0] == '\0') {
      perror(ERR_NOFNAME);
      exit(EXIT_FAILURE);
      /* NOTREACHED */
    }

    if (do_the_needful(&options) != EXIT_SUCCESS) {
       perror(ERR_DO_THE_NEEDFUL);
       exit(EXIT_FAILURE);
       /* NOTREACHED */
    }

    printf("-------- optarg debug --------\n");
    printf("delta_z    = %f\n", options.delta_z);
    printf("max_frames = %llu\n", options.max_frames);
    printf("fname      = %s\n", options.fname);
    printf("stride     = %d\n", options.stride);
    printf("-------- gubed gratpo --------\n\n");

    read_xtc_natoms(options.fname, &natoms);
    nmol = natoms / 3;

    // read_xtc_n_frames(options.fname, &n_frames, &est_nframes, offsets);

    printf("# %d atoms corresponding to %d molecules\n", natoms, nmol);

    if (!(xd = xdrfile_open(options.fname, "r"))) {
      perror(ERR_FOPEN_INPUT);
      exit(EXIT_FAILURE);
      /* NOTREACHED */
    }
    x = calloc(natoms, sizeof(x[0]));
    if ((read_xtc(xd, natoms, &step, &time, box, x, &prec))) {
      perror(ERR_FOPEN_INPUT);
      exit(EXIT_FAILURE);
      /* NOTREACHED */
    }

    printf("-------- xdr debug --------\n");
    printf("n_frames = %d\n", n_frames);
    printf("natoms = %d\n", natoms);
    printf("step   = %d\n", step);
    printf("time   = %f\n", time);
    printf("box    = [ %f, %f, %f ]\n", box[0][0], box[1][1], box[2][2]);
    printf("prec   = %f\n", prec);
    printf("-------- gubed rdx --------\n\n");

    curr = calloc(natoms, sizeof(curr[0]));
    prev = calloc(natoms, sizeof(prev[0]));

    nframes = (options.max_frames < n_frames) ? options.max_frames : n_frames;
    r_msd_tau = calloc(nframes, sizeof(r_msd_tau[0]));
    for (i = 0; i < DIM; i++) {
      r_msd_tau[0][i] = 0;
    }
    alpha_msd_tau = calloc(nframes, sizeof(alpha_msd_tau[0]));
    for (i = 0; i < DIM; i++) {
      alpha_msd_tau[0][i] = 0;
    }
    phi_msd_tau = calloc(nframes, sizeof(phi_msd_tau[0]));
    for (i = 0; i < DIM; i++) {
      phi_msd_tau[0][i] = 0;
    }

    do {
      if (i_frame == 1)
        delta_t = time;
      if (options.max_frames > 0 && i_frame >= options.max_frames) {
        printf("\nMaximum number of frames reached.\n");
        break;
      }
      if (i_frame % options.stride == 0) {
        printf("\rReading frame %d of %d, t = %.2f ps ", i_frame+1, n_frames, time);
        fflush(stdout);
        if (i_frame == 0) {
          memcpy(curr, x, natoms*sizeof(x[0]));
          l = box[0][0];
        }
        else if (i_frame > 0) {
          i_mol = options.i_mol;
          memcpy(prev, curr, natoms*sizeof(x[0]));
          memcpy(curr, x, natoms*sizeof(x[0]));
          l_prev = l;
          l = box[0][0];
          get_msd(i_mol, curr, prev, l, l_prev, options.delta_z,
                  delta_r, delta_alpha, delta_phi);
          rxpyz(r_msd_tau[i_frame-1], delta_r, r_msd_tau[i_frame]);
          rxpyz(alpha_msd_tau[i_frame-1], delta_alpha, alpha_msd_tau[i_frame]);
          rxpyz(phi_msd_tau[i_frame-1], delta_phi, phi_msd_tau[i_frame]);
        }
      }
      i_frame++;
    } while(!read_xtc(xd, natoms, &step, &time, box, x, &prec));

    free(x);
    free(curr);
    free(prev);

    printf("\nTau averaging ...\n");
    fflush(stdout);
    n_tau = calloc(nframes, sizeof(long int));
    r_msd = calloc(nframes, sizeof(float));
    alpha_msd = calloc(nframes, sizeof(alpha_msd[0]));
    phi_msd = calloc(nframes, sizeof(phi_msd[0]));
    tau_avrg(r_msd_tau, alpha_msd_tau, phi_msd_tau, nframes,
             n_tau, r_msd, alpha_msd, phi_msd);
    free(r_msd_tau);
    free(alpha_msd_tau);
    free(phi_msd_tau);

    printf("Writing to file %s ...\t", options.oname);
    fflush(stdout);
    fp = fopen(options.oname, "w");
    for (i = 0; i < nframes; i++) {
      fprintf(fp, "%3.3g %3.8g %3.8g %3.8g %3.8g %3.8g %3.8g %3.8g\n",
              i * delta_t, r_msd[i] / n_tau[i], alpha_msd[i][0] / n_tau[i],
              alpha_msd[i][1] / n_tau[i], alpha_msd[i][2] / n_tau[i],
              phi_msd[i][0] / n_tau[i], phi_msd[i][1] / n_tau[i],
              phi_msd[i][2] / n_tau[i]);
    }
    fclose(fp);
    free(n_tau);
    free(r_msd);
    free(alpha_msd);
    free(phi_msd);

    printf("done\n");
    fflush(stdout);

    return EXIT_SUCCESS;
}
