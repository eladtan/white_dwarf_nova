subroutine nse (ro,t,a,y)
!
  use global
  include 'real8.com'
  include 'network.com'
!
  dimension a(nmat),y(nmat)
  dimension x(nmat),yn(nmat),b(nmat)
! -----------------------------------------------------
  knisa=0
10 continue
  knisa=knisa+1
  iter=0
  y(1)=1.d0/a(1)
!xxx  y(1)=1.d-10   !!!
  t9=t*1.d-9
  t932=t9**1.5
  ymin=1.d-100
  ymax=1.d0
!
  qfac=emev_nuc/(boltz*avo)*1.d-9
  cc=9.867425d9/ro*t932
  qval=(3.d0*excess(1)-excess(2))*qfac
  ratio=(a(1)**3/a(2))**1.5
  b(2)=cc**2*ratio*exp(qval/t9)
  bfac=b(2)*ro**2
  do m=3,nmat
     qval=(excess(1)+excess(m-1)-excess(m))*qfac
     ratio=(a(1)*a(m-1)/a(m))**1.5
     b(m)=cc*ratio*exp(qval/t9)
     bfac=b(m)*ro
  enddo
! -------------------------------------------------------------
20 continue
  iter=iter+1
  if(iter.gt.300) go to 77
!
  sumx=y(1)*a(1)
  dsum=a(1)
  y(2)=y(1)**3/b(2)
  dymdy1=3.d0*y(2)/y(1)
  sumx=sumx+y(2)*a(2)
  dsum=dsum+dymdy1*a(2)
  do m=3,nmat
     y(m)=y(1)*y(m-1)/b(m)
     dymdy1=(y(m-1)+y(1)*dymdy1)/b(m)
     sumx=sumx+y(m)*a(m)
     dsum=dsum+dymdy1*a(m)
  enddo
!
!xxx  yn1=y(1)-(sumx-1.d0)/dsum
  yn1=y(1)-(sumx-xnorm)/dsum
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  if(sumx.gt.xnorm) then
     ymax=y(1)
  else
     ymin=y(1)
  endif
  xlg=(sumx/xnorm)**0.2d0
  yn1=y(1)/xlg
  if(abs(yn1-y(1)).gt.(ymax-ymin)*0.8) then
     yn1=(ymax+ymin)/2
  endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  
  error=abs(sumx-xnorm)
  if(knisa.gt.1) then
     if(iter.eq.1) then
        do m=2,nmat
           print*,' m,b(m) ',m,b(m)
        end do
     end if
     print 22,iter,error,ro,t9,y(1),yn1
22   format(' NSE : iter,error,ro,t9,y(1),yn1 ',i6,10es12.4)
  endif
!!!xxx  y(1)=min(max(yn1,y(1)/2),2*y(1))
  y(1)=yn1
  if(error.gt.1.d-8) go to 20
! -----------------    convergence  ---------------------------
  y(2)=y(1)**3/b(2)
  do m=3,nmat
     y(m)=y(1)*y(m-1)/b(m)
  enddo
!
  return
! -------------------------------------------------------------
!                     iteration diverge
77 continue
  if(knisa.lt.2) go to 10
  stop' nse diverges'
! --------------------------------------------------------------
end subroutine nse
