  type(sqrm_spmat_type)  :: qrm_mat
  real(kind(1.e0))       :: b(:,:)
  real(kind(1.e0))       :: x(:,:)

  type(sqrm_rhs_type)    :: x_rhs, b_rhs
  type(qrm_dscr_type)    :: qrm_dscr

  call qrm_init()

  !init the matrix data structure
  call qrm_spmat_init(qrm_mat)

  ! fill up the matrix and b
  ! ...

  ! init the rhs data structures
  call qrm_rhs_init(b_rhs, b)
  call qrm_rhs_init(x_rhs, x)
  
  ! init the descriptor data structure
  call qrm_dscr_init(qrm_dscr)

  ! submit analysis, facto, apply and solve operations
  call qrm_analyse_async(qrm_dscr, qrm_mat, 'n')
  call qrm_factorize_async(qrm_dscr, qrm_mat, 'n')

  call qrm_apply_async(qrm_dscr, qrm_mat, 't', b_rhs, err)
  call qrm_solve_async(qrm_dscr, qrm_mat, 'n', b_rhs, x_rhs, err)

  ! wait for their completion
  call qrm_barrier(qrm_dscr)

  ! cleanup
  call qrm_dscr_destroy(qrm_dscr)
  call qrm_rhs_destroy(b_rhs)
  call qrm_rhs_destroy(x_rhs)
  call qrm_spmat_destroy(qrm_mat)
