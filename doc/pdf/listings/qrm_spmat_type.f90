type sqrm_spmat_type
   ! Row and column indices
   integer, pointer :: irn(:), jcn(:)
   ! Numerical values
   real, pointer :: val(:)
   ! Number of rows, columns
   ! and nonzeroes
   integer :: m, n, nz
   ! A pointer to an array 
   ! containing a column permutation 
   ! provided by the user
   integer, pointer :: cperm_in(:)
   ! Integer control parameters
   integer :: icntl(20)
   ! Collected statistics
   integer(kind=8) :: gstats(10)
end type sqrm_spmat_type
