subroutine tmp_prs
  use tables
  include 'real8.com'
! -------------------------------------------------------
  do it=1,itmax
!
     call rho_tmp
!
     drho = (prs-prs0)/dpdrho
     rho = min(max(rho - drho,rho/2),2*rho)
     if(abs(prs-prs0).lt.diuk*prs0) then
        keyerr=0
        return
     endif
  enddo
  write(*,*) ' tmp_prs diverged '
  keyerr=12
! --------------------------------------------
end subroutine tmp_prs
