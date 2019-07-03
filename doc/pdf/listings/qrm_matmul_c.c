void sqrm_matmul_c(qrm_mat_c, transp,
                   alpha, x, 
                   beta, y, 
                   nrhs);
struct sqrm_spmat_type_c *qrm_mat_c,;
const char transp;
const float alpha;
float *x; 
const float beta;
float *y;
const int nrhs;
