/* main.c */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <libgen.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <sys/types.h>

#include "xdrfile_xtc.h"

#define OPTSTR "f:z:n:h"
#define USAGE_FMT  "%s [-f trajectory] [-z delta_z] [-n number of frames] [-h]\n"
#define ERR_FOPEN_INPUT  "fopen(input, r)"
#define ERR_FOPEN_OUTPUT "fopen(output, w)"
#define ERR_DO_THE_NEEDFUL "do_the_needful blew up"
#define DEFAULT_PROGNAME "locmsd"

extern int errno;
extern char *optarg;
extern int opterr, optind;

typedef struct {
  // FILE         *input;
  float         delta_z;
  uint64_t      nframes;
  char          fname[300];
  int           batch;
} options_t;

int dumb_global_variable = -11;

void usage(char *progname, int opt);
int  do_the_needful(options_t *options);

int main(int argc, char *argv[]) {
    int opt;
    options_t options = { 0, 0, "" , 0};

    int natoms, nmol;
    int mol;
    int molbatch = 0;

    opterr = 0;

    while ((opt = getopt(argc, argv, OPTSTR)) != EOF)
       switch(opt) {
           case 'f':
              sscanf(optarg, "%s", options.fname);
              // if (!(options.input = fopen(optarg, "r")) ){
              //    perror(ERR_FOPEN_INPUT);
              //    exit(EXIT_FAILURE);
              //    /* NOTREACHED */
              // }
              break;

           case 'z':
              options.delta_z = strtof(optarg, NULL);
              break;

           case 'n':
              options.nframes = (uint64_t) strtoul(optarg, NULL, 10);
              break;

           case 'b':
              options.batch += 1;

           case 'h':
           default:
              usage(basename(argv[0]), opt);
              /* NOTREACHED */
              break;
       }

    if (do_the_needful(&options) != EXIT_SUCCESS) {
       perror(ERR_DO_THE_NEEDFUL);
       exit(EXIT_FAILURE);
       /* NOTREACHED */
    }

    printf("-------- optarg debug --------\n");
    printf("delta_z = %f\n", options.delta_z);
    printf("nframes = %llu\n", options.nframes);
    printf("fname   = %s\n", options.fname);
    printf("-------- gubed gratpo --------\n\n");

    read_xtc_natoms(options.fname, &natoms);
    nmol = natoms / 3;

    printf("# %d atoms corresponding to %d molecules\n", natoms, nmol);

    for (mol = 0; mol + molbatch <= nmol; mol+=molbatch) {
      printf("# Reading positions of molecules %d-%d\n", mol, mol+molbatch-1);
      fflush(stdout);
    }

    return EXIT_SUCCESS;
}

void usage(char *progname, int opt) {
   fprintf(stderr, USAGE_FMT, progname?progname:DEFAULT_PROGNAME);
   exit(EXIT_FAILURE);
   /* NOTREACHED */
}

int do_the_needful(options_t *options) {

   if (!options) {
     errno = EINVAL;
     return EXIT_FAILURE;
   }

   // if (!options->input || !options->output) {
   //   errno = ENOENT;
   //   return EXIT_FAILURE;
   // }

   /* XXX do needful stuff */

   return EXIT_SUCCESS;
}
