void sqrm_residual_norm_c(qrm_spmat_c, 
                          b, x, 
                          nrhs, nrm);
struct sqrm_spmat_type_c *qrm_spmat_c;
float *b;
float *x;
const int nrhs;
float *nrm;
