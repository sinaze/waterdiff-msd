#define OPTSTR "f:z:n:s:i:o:h"
#define USAGE_FMT  "%s -f trajectory [-z delta_z] [-n number of frames] [-s stride] [-i molecule no.] [-o output] [-h]\n"
#define ERR_FOPEN_INPUT  "fopen(input, r)"
#define ERR_FOPEN_OUTPUT "fopen(output, w)"
#define ERR_NOFNAME "-f is mandatory!\n"
#define ERR_DO_THE_NEEDFUL "do_the_needful blew up"
#define DEFAULT_PROGNAME "locmsd"

typedef struct {
  float         delta_z;
  long long     max_frames;
  char          fname[300];
  int           stride;
  int           i_mol;
  char          oname[300];
} options_t;

void usage(char *progname, int opt);

int  do_the_needful(options_t *options);
