interface qrm_least_squares
   
   subroutine sqrm_least_squares1d(qrm_mat, b, x, info)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: b(:)
     real                  :: x(:)
     integer, optional     :: info
   end subroutine sqrm_least_squares1d

   subroutine sqrm_least_squares2d(qrm_mat, b, x, info)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: b(:,:)
     real                  :: x(:,:)
     integer, optional     :: info
   end subroutine sqrm_least_squares2d

end interface qrm_least_squares
