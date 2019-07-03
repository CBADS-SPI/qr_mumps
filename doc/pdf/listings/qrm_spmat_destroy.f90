interface qrm_spmat_destroy

   subroutine sqrm_spmat_destroy(qrm_mat, all, info)
     type(sqrm_spmat_type) :: qrm_mat
     logical, optional     :: all
     integer, optional     :: info
   end subroutine sqrm_spmat_destroy

end interface qrm_spmat_destroy
