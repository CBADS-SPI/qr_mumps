interface qrm_get

   subroutine sqrm_pgeti(qrm_mat, string, ival, info)
     type(sqrm_spmat_type) :: qrm_mat
     character(len=*)      :: string
     integer               :: ival
     integer, optional     :: info
   end subroutine sqrm_pgeti
   
   subroutine sqrm_pgetii(qrm_mat, string, ival, info)
     type(sqrm_spmat_type) :: qrm_mat
     character(len=*)      :: string
     integer(kind=8)       :: ival
     integer, optional     :: info
   end subroutine sqrm_pgetii
   
   subroutine qrm_ggeti(string, ival, info)
     character(len=*)      :: string
     integer               :: ival
     integer, optional     :: info
   end subroutine qrm_ggeti

   subroutine qrm_ggetii(string, ival, info)
     character(len=*)      :: string
     integer(kind=8)       :: ival
     integer, optional     :: info
   end subroutine qrm_ggetii

end interface qrm_get
