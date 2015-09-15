subroutine nse_equilib (indxeos,rho,tmp0,enr,x,a,matters  &
     ,aa,zz,de_nse,tmpsaf_nse,key_done)
!
  use global
  include 'real8.com'
  include 'network.com'
  dimension x(matters),a(matters)
  dimension y(matters),y0(matters),y1(matters),y2(matters)
  real(8),save :: eps_converge=1.d-4
  integer,save :: it_max=50
! ----------------------------------------------------
  key_done=0 
  keypr=0
  iter=0
  enr0=enr
  tmp=tmp0
  de_nse=0.d0
  enr_min=0.1*enr
  enr_max=10*enr
  do m=1,matters
     y0(m)=x(m)/a(m)
  enddo
  tnse_min=0.d0
  tnse_max=1.d12
!
  keyte=0
  phasex=0.d0
  stop 'implementation revoked'
!  call eq_of_state (0,'nse_equilib0'            &
!       ,indxeos,keyte,aa,zz,rho,enr,tmp,prs     &
!       ,x,phasex                                &
!       ,ss,dpdro,dpde,dedro,dedt,soundsqr,keyerr)
! ============ iterations of final tmp =================
10 continue
  iter=iter+1
  if(iter.gt.30) keypr=1  !!!!!!!!!!!
!
  phasex=0.d0
  stop 'implementation revoked'
!  call eq_of_state (0,'nse_equilib1'            &
!       ,indxeos,keyte,aa,zz,rho,enr,tmp,prs   &
!       ,x,phasex                              &
!       ,ss,dpdro,dpde,dedro,dedt,soundsqr,keyerr)
!
  tmp1=max(tmp,3.d9)  !!! ad-hoc
  dtmp=1.d-2*tmp1
  call nse (rho,tmp1,a,y1)
  tmp2=tmp1+dtmp
  call nse (rho,tmp2,a,y2)
  if(keypr.ne.0) then
     do m=1,-matters
        print 1234,iter,m,a(m),y1(m),y2(m)
     enddo
1234 format(' NSE_equilib : iter,m,a,y1,y2 ',2i5,10es12.4)
  end if
! -------------------------------------------------------------------
  de_nse1=sum((y1(1:matters)-y0(1:matters))*excess(1:matters)) * emev_nuc
  de_nse2=sum((y2(1:matters)-y0(1:matters))*excess(1:matters)) * emev_nuc
! -------------------------------------------------------------------
  de_nse=de_nse1
  de_nse_dtmp=(de_nse2-de_nse1)/dtmp
  de_nse_denr=de_nse_dtmp/dedt
!
  fq=de_nse-(enr-enr0)
  if(abs(fq).gt.eps_converge*enr) then
     dfde=de_nse_denr-1.d0
     de=-fq/dfde
     if(de.gt.0.d0) enr_min=max(enr_min,enr)
     if(de.le.0.d0) enr_max=min(enr_max,enr)
     enrn=enr+de
     if(abs(de).gt.0.8*(enr_max-enr_min)) then
        enrn=(enr_min+enr_max)/2
     endif
     if(keypr.ne.0) print 71,iter,rho,tmp,fq/enr,enr0,enr,enrn,de,enr_min,enr_max
71   format(' iter,rho,tmp,fq/enr0,enr0,enr,enrn,de,enr_min,enr_max ',i5,10es16.8)
     enr=min(max(0.5d0*enr,enrn),1.5d0*enr)
     if(iter.gt.it_max) stop 'iter gt it_max in nse_equilib'
     go to 10
  endif
! -------------    convergence  -----------------------------
  key_done=1
  do m=1,matters
     x(m)=y1(m)*a(m)
  enddo
  sumx=sum(x(1:matters))
  x(1:matters)=x(1:matters)/sumx * xnorm
! ----------------------------------------------------
end subroutine nse_equilib
