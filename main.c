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
    options_t options = { 0.0, 0, "\0" , 1 };

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

    int i_frame = 0;
    int t_frame;
    rvec *curr;
    rvec *prev;

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
              break;

           case 's':
              options.stride = atoi(optarg);
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

    read_xtc_n_frames(options.fname, &n_frames, &est_nframes, offsets);

    printf("# %d atoms corresponding to %d molecules\n", natoms, nmol);

    // for (mol = 0; mol+molbatch <= nmol; mol+=molbatch) {
    //   printf("# Reading positions of molecules %d-%d\n", mol, mol+molbatch-1);
    //   fflush(stdout);
    // }

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

    do {
      if (i_frame == 1)
        t_frame = time;
      if (options.max_frames > 0 && i_frame >= options.max_frames) {
        printf("\nMaximum number of frames reached.\n");
        break;
      }
      if (i_frame % options.stride == 0) {
        printf("\rReading frame %d of %d, t = %.2f ps ", i_frame+1, n_frames, time);
        fflush(stdout);

        if (i_frame == 0) {
          memcpy(curr, x, natoms*DIM);
          // printf("\nx[1] = [ %f %f %f ]\n", x[1][0], x[1][1], x[1][2]);
          // printf("curr[1] = [ %f %f %f ]\n", curr[1][0], curr[1][1], curr[1][2]);
        }
        else if (i_frame > 0) {
          memcpy(prev, curr, natoms*DIM);
          memcpy(curr, x, natoms*DIM);
          // printf("\nx[1] = [ %f %f %f ]\n", x[1][0], x[1][1], x[1][2]);
          // printf("curr[1] = [ %f %f %f ]\n", curr[1][0], curr[1][1], curr[1][2]);
          // printf("prev[1] = [ %f %f %f ]\n", prev[1][0], prev[1][1], prev[1][2]);
        }
        get_msd(curr, prev);

      }
      i_frame++;
    } while(!read_xtc(xd, natoms, &step, &time, box, x, &prec));

    return EXIT_SUCCESS;
}
