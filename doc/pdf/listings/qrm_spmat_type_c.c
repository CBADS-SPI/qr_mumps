struct sqrm_spmat_type_c{
  /* Row and column indices */
  int *irn, *jcn;
  /* Numerical values */
  float *val;
  /* Number of rows, columns 
     and nonzeroes */
  int m, n, nz;
  /* A pointer to an array 
     containing a column permutation 
     provided by the user */
   int *cperm_in;
  /* Integer control parameters */
  int icntl[20];
  /* Collected statistics */
  long int gstats[10];
};
