#define INC 1

void get_msd(const rvec *curr, const rvec *prev, const float l,
             rvec delta_r, rvec delta_alpha, rvec delta_phi);

void rxcrossyz(const rvec x, const rvec y, rvec z);

void rxpyz(const rvec x, const rvec y, rvec z);

void print_rvec(const rvec x, const char *name);

void get_eckart(const rvec o, const rvec h1, const rvec h2, const float l,
                rvec x, rvec y, rvec z);

float pbc_corr(const float x, const float l);

void get_delta(const rvec ow, const rvec ow_prev, const float delta_z,
               const rvec x, const rvec y, const rvec z,
               const rvec x_prev, const rvec y_prev, const rvec z_prev,
               const float l,
               rvec delta_r, rvec delta_alpha, rvec delta_phi);
