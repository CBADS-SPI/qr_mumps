interface qrm_factorize

   subroutine sqrm_factorize(qrm_mat, transp, info)
     type(sqrm_spmat_type):: qrm_mat
     character, optional  :: transp
     integer, optional    :: info
   end subroutine sqrm_factorize

end interface qrm_factorize
