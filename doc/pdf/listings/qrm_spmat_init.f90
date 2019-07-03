interface qrm_spmat_init

   subroutine sqrm_spmat_init(qrm_spmat, info)
     type(sqrm_spmat_type) :: qrm_mat
     integer, optional     :: info
   end subroutine sqrm_spmat_init

end interface qrm_spmat_init
