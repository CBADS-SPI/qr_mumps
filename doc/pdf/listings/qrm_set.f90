interface qrm_set

   subroutine sqrm_pseti(qrm_mat, string, ival, info)
     type(sqrm_spmat_type) :: qrm_mat
     character(len=*)      :: string
     integer               :: ival
     integer, optional     :: info
   end subroutine sqrm_pseti

   subroutine qrm_gseti(string, ival, info)
     character(len=*)      :: string
     integer               :: ival
     integer, optional     :: info
   end subroutine qrm_gseti

end interface qrm_set
