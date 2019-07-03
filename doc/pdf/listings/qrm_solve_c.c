void sqrm_solve_c(qrm_mat_c, transp,
                  b, x, nrhs);
struct sqrm_spmat_type_c *qrm_mat_c;
const char transp;
float *b;
float *x;
const int nrhs;
