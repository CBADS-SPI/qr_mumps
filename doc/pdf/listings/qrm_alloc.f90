interface qrm_alloc

   subroutine qrm_aalloc_s(a, m, info)
     real(kind(1.e0)), allocatable  :: a(:)
     integer                        :: m
     integer, optional              :: info
   end subroutine qrm_aalloc_s
   
   subroutine qrm_aalloc_2s(a, m, n, info)
     real(kind(1.e0)), allocatable  :: a(:,:)
     integer                        :: m, n
     integer, optional              :: info
   end subroutine qrm_aalloc_2s

   subroutine qrm_palloc_s(a, m, info)
     real(kind(1.e0)), pointer      :: a(:)
     integer                        :: m
     integer, optional              :: info
   end subroutine qrm_palloc_s
   
   subroutine qrm_palloc_2s(a, m, n, info)
     real(kind(1.e0)), pointer      :: a(:,:)
     integer                        :: m, n
     integer, optional              :: info
   end subroutine qrm_palloc_2s

end interface qrm_alloc


interface qrm_dealloc
   
   subroutine qrm_adealloc_s(a, info)
     real(kind(1.e0)), allocatable  :: a(:)
     integer, optional              :: info
   end subroutine qrm_adealloc_s
   
   subroutine qrm_adealloc_2s(a, info)
     real(kind(1.e0)), allocatable  :: a(:,:)
     integer, optional              :: info
   end subroutine qrm_adealloc_2s
   
   subroutine qrm_pdealloc_s(a, info)
     real(kind(1.e0)), pointer      :: a(:)
     integer, optional              :: info
   end subroutine qrm_pdealloc_s
   
   subroutine qrm_pdealloc_2s(a, info)
     real(kind(1.e0)), pointer      :: a(:,:)
     integer, optional              :: info
   end subroutine qrm_pdealloc_2s
   
end interface qrm_dealloc
