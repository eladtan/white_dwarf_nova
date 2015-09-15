subroutine enr_prs
  use tables
  include 'real8.com'
!   -------------------------------------------------------
  do 100 it=1,itmax
!
     call rho_tmp
!
     det=dpdrho*dedtmp-dpdtmp*dedrho
     if(abs(det).lt.1.d-12) stop' zero det in eospe'
     dp=prs-prs0
     de=enr-enr0
     if(abs(dp).lt.diuk*prs0.and.abs(de).lt.diuk*abs(enr0)) then
        keyerr=0
        return
     endif
     drho = -(dedtmp*dp-dpdtmp*de)/det
     dtmp = -(dpdrho*de-dedrho*dp)/det
!
     rhon=rho+drho
     tmpn=tmp+dtmp
     delta=max(abs(rhon-rho)/rho,abs(tmpn-tmp)/tmp)+1.d-3
     fac=min(1.d0,0.5/delta)
     rho=rho+fac*(rhon-rho)
     tmp=tmp+fac*(tmpn-tmp)
!
100 end do
  write(*,*) ' enr_prs diverged '
  keyerr=12
!     --------------------------------------------
end subroutine enr_prs


