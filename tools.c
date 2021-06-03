#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "tools.h"

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
