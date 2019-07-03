interface qrm_residual_norm

   subroutine sqrm_residual_norm1d(qrm_mat, b, x, nrm, info)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: b(:)
     real                  :: x(:)
     real                  :: nrm
     integer, optional     :: info
   end subroutine sqrm_residual_norm1d

   subroutine sqrm_residual_norm2d(qrm_mat, b, x, nrm, info)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: b(:,:)
     real                  :: x(:,:)
     real                  :: nrm
     integer, optional     :: info
   end subroutine sqrm_residual_norm2d

end interface qrm_residual_norm
