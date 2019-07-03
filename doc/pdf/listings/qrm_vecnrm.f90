interface qrm_vecnrm

   subroutine sqrm_vecnrm1d(vec, n, ntype, nrm)
     real              :: vec(:)
     integer           :: n
     character         :: ntype
     real              :: nrm
     integer, optional :: info
   end subroutine sqrm_vecnrm1d

   subroutine sqrm_vecnrm2d(vec, n, ntype, nrm)
     real              :: vec(:,:)
     integer           :: n
     character         :: ntype
     real              :: nrm(:)
     integer, optional :: info
   end subroutine sqrm_vecnrm2d

end interface qrm_vecnrm
