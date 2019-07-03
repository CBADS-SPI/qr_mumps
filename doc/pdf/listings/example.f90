program sqrm_test_small
  use sqrm_mod
  implicit none

  type(sqrm_spmat_type)  :: qrm_mat
  integer                :: ierr, nargs, i, nrhs
  real, allocatable      :: b(:), x(:), r(:)
  real                   :: rnrm, onrm

  call qrm_init()
  
  ! initialize the control data structure. 
  call qrm_spmat_init(qrm_mat)

  ! allocate arrays for the input matrix
  call qrm_alloc(qrm_mat%irn, 13)
  call qrm_alloc(qrm_mat%jcn, 13)
  call qrm_alloc(qrm_mat%val, 13)
  ! initialize the input matrix
  qrm_mat%jcn = (/1,1,1,2,2,3,3,3,3,4,4,5,5/)
  qrm_mat%irn = (/2,3,6,1,6,2,4,5,7,2,3,2,4/)
  qrm_mat%val = (/0.7,0.6,0.4,0.1,0.1,0.3,0.6,0.7,0.2,0.5,0.2,0.1,0.6/)
  qrm_mat%m   = 7
  qrm_mat%n   = 5
  qrm_mat%nz  = 13

  write(*,'("Starting Analysis")')
  call qrm_analyse(qrm_mat)

  write(*,'("Starting Factorization")')
  call qrm_factorize(qrm_mat)

  call qrm_alloc(b, qrm_mat%m)
  call qrm_alloc(r, qrm_mat%m)
  call qrm_alloc(x, qrm_mat%n)
  
  b = 1.e0
  ! as by is changed when applying Q', we save a copy in r for later use
  r = b
  call qrm_apply(qrm_mat, 't', b)
  call qrm_solve(qrm_mat, 'n', b, x)

  ! compute the residual
  call qrm_residual_norm(qrm_mat, r, x, rnrm)
  call qrm_residual_orth(qrm_mat, r, onrm)   
  write(*,'("||r||/||A||    = ",e10.2)')rnrm
  write(*,'("||A^tr||/||r|| = ",e10.2)')onrm

  call qrm_dealloc(b)
  call qrm_dealloc(r)
  call qrm_dealloc(x)
  call qrm_spmat_destroy(qrm_mat, all=.true.)
  write(*,'("  Nonzeroes in R           : ",i20)')qrm_mat%gstats(qrm_nnz_r_)
  write(*,'("  Total flops at facto     : ",i20)')qrm_mat%gstats(qrm_e_facto_flops_)
  write(*,'("  Global memory peak       : ",f9.3," MB")') &
       &real(qrm_max_mem,kind(1.d0))/1024.d0/1024.d0

  call qrm_finalize()
  
  stop
end program sqrm_test_small

