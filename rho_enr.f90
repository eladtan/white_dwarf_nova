subroutine rho_enr
  use tables
  include 'real8.com'
! -------------------------------------------------------
  t0=tmp0
  knisa=0
  keyerr=0
10 continue
  knisa=knisa+1
  tmp=t0
  tmp_cold=tmp_min*2.d0
  tmpmin=tmp_cold
  tmpmax=tmp_max
!
  do 100 it=1,itmax
!
     call rho_tmp
     if(knisa.gt.1) write(*,11) it,tmp,enr,enr0,(enr0-enr)/enr0
11   format('rho_enr-it,tmp,enr,enr0,de_rel ',i5,6es12.4)
!
     if(abs(enr-enr0).lt.diuk*abs(enr0)) then
        keyerr=0
        return
     endif
!
     dt = -(enr-enr0)/dedtmp
     if(dedtmp.le.0.d0) then
        print*,' neg dedtmp - rho,tmp,enr,dedtmp ' &
             ,rho,tmp,enr,dedtmp
        if(knisa.lt.2) go to 10
        return
     endif
     if(knisa.gt.1) write(*,12) dt,tmpmin,tmpmax
12   format(' dt,tmpmin,tmpmax ',10es18.10)
     if(dt.gt.0.d0) tmpmin=tmp
     if(dt.le.0.d0) tmpmax=tmp
     tmpn=tmp+dt
     tmpn=min(tmpn,tmpmax)
     tmpn=max(tmpn,tmpmin)
     adt=abs(tmpn-tmp)
     if(it.eq.1) go to 80
     if(adt.lt.(tmpmax-tmpmin)*0.6d0) go to 80
     tmpn=0.5d0*(tmpmin+tmpmax)
80   continue
     tmp=max(min(tmpn,1.2*tmp),0.8*tmp)
!
     if(tmp.lt.1.001*tmp_cold) then
        tmp=tmp_cold
        call rho_tmp
        keyerr=-1
        return
     endif
!
100 enddo
!
  write(*,*) ' rho_enr diverged - rho,tmp,e0,en,p ' &
       ,  rho,tmp,enr0,enr,prs
  if(knisa.lt.2) go to 10
!!!  stop ' RHO_ENR DIVERGES '   !!!   DEBUG
!     --------------------------------------------
end subroutine rho_enr
