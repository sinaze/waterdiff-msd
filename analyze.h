#define INC 1

int get_msd(const rvec *curr, const rvec *prev, const float l);

void rxcrossyz(const rvec x, const rvec y, rvec z);

void print_rvec(const rvec x, const char *name);

void get_eckart(const rvec o, const rvec h1, const rvec h2, const float l,
                rvec x, rvec y, rvec z);
