subroutine rho_prs
  use tables
  include 'real8.com'
!   -------------------------------------------------------
  do 100 it=1,itmax
!
     call rho_tmp
!
     dt = (prs-prs0)/dpdtmp
     tmp = min(max(tmp - dt,0.9*tmp),1.1*tmp)
!          print *,' it,dt,t ',it,dt,tmp
     if(abs(prs-prs0).lt.diuk*prs0) then
        keyerr=0
        return
     endif
100 enddo
  keyerr=12
  write(*,*) "rho_prs did not converge"
  stop
!     --------------------------------------------
end subroutine rho_prs
