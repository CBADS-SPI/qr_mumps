interface qrm_matmul

   subroutine sqrm_matmul1d(qrm_mat, transp, alpha, x, beta, y)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: y(:)
     real                  :: x(:)
     real                  :: alpha, beta
     character             :: transp
   end subroutine sqrm_matmul1d
   
   subroutine sqrm_matmul2d(qrm_mat, transp, alpha, x, beta, y)
     type(sqrm_spmat_type) :: qrm_mat
     real                  :: y(:,:)
     real                  :: x(:,:)
     real                  :: alpha, beta
     character             :: transp
   end subroutine sqrm_matmul2d
   
end interface qrm_matmul
