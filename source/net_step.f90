subroutine net_step (rho,tmp,x,a,z,dedtmp,matters,dtime,qreac,delx,key_done)
!
  use global
  include 'real8.com'
  include 'network.com'
!      parameter (nd=matters)
  dimension x(matters),a(matters),xn(matters)
  dimension rr(matters+1,matters+1),rhs(matters)
  dimension y(matters),yn(matters),y0(matters),y00(matters)
  dimension dy(matters),dydt(matters),dyidyj(matters,matters+1)
  dimension dx(matters)
  data xneg/-1.d-3/
! -----------------------------------------------------------------
  qreac=0.d0
  t9=min(tmp/1.d9,20.d0)
  if(t9.lt.0.05) then
     key_done=1
     return
  end if
! --------------------------------------------------------------------
  do  m=1,matters
     y00(m)=x(m)/a(m)
     y(m)=y00(m)
  enddo
  tim=0.d0
  dt=dtime
  mat=matters
  explic=0.5d0
  if(tmp.gt.2.d9) explic=0.d0
  emplic=1.d0-explic
  nsteps=0
!  ----------------------  begin dt step  ----------------------------
7 continue
  nsteps=nsteps+1
  nreps=0
  keypr=0
  dt=min(dt,dtime-tim)
  if(dtime-(tim+dt).lt.0.3*dt) dt=dtime-tim
  y0(1:matters)=y(1:matters)
! ---------------------- repeat time step ----------------------------
9 continue
  nreps=nreps+1
  iter=0
!xxxxxxxxxxxxxx
  if(nreps.gt.30.or.dt.lt.1.d-6*dtime) then
     print*,' BURN_STEP Divereges : nsteps,nreps,dtime,dt ',nsteps,nreps,dtime,dt
     print 323,rho,tmp
323  format(' BURN_STEP Divereges : rho,tmp ',10es12.4)
     print*,'network_master, nreps,sumx,del ',nreps,sumx,del
!xxx      write(*,323) iter,( y(m)*a(m),m=1,matters)
!xxx      write(*,434) iter,(dy(m)*a(m),m=1,matters)
324  format(' iter, x ',i3,100es12.4)
434  format(' iter,dx ',i3,100es12.4)
!xxx     stop ' BURN_STEP Divereges '
     key_done=0
     return
  endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
! rrrrrrrrrrrrrrrrrrrrrrrrr  compute rates  rrrrrrrrrrrrrrrrrrrrrr
!                       sigmav and screening
  aa=1.d0/sum(y(1:matters))
  zz=aa/2
  call sigmav_rates (matters,rho,tmp,aa,zz)
!                           iterations
  del=0.d0
10 continue
  iter=iter+1
  prev_del=del
  do m=1,matters
     yn(m)=explic*y0(m)+emplic*y(m)
  enddo
!
  call alpha_rates (matters,rho,tmp,yn,dydt,dyidyj)
! -------------------------------------------------------------
!     build the linear system for Newton iteration
  rr(:,:)=0.d0
  do m=1,matters
     rhs(m)=y0(m)-y(m)+dt*dydt(m)
     rr(m,m)=1.d0
  enddo
  xdt=emplic*dt
  do m1=1,matters 
!xxx     print 37,m1,yn(m1),dydt(m1),rhs(m1)
37   format(' m,yn,dydt,rhs ',i4,10es12.4)
     do m2=1,matters
        rr(m1,m2)=rr(m1,m2)+dyidyj(m1,m2)*xdt
     enddo
  enddo
!  -----------------------------------------------------------
!                           solve the linear system
  key=3
  call linsys_burn (rr,rhs,dy,key,matters,matters+1)
!
  sumx_old=0.d0
  sumx=0.d0
  xneg=0.d0
  dx(1:matters)=dy(1:matters)*a(1:matters)
  sumdx=sqrt(sum(dx(1:matters)**2))
  facdy=min(1.d0,0.1d0/(sumdx+1.d-6))
  dy(1:matters)=dy(1:matters)*facdy
!
  do m=1,matters
     ynm=y(m)+dy(m)
     xneg=min(xneg,ynm*a(m))
     yn(m)=max(y(m)+dy(m),0.d0)
     xn(m)=yn(m)*a(m)
     sumx_old=sumx_old+x(m)
     sumx=sumx+xn(m)
  enddo
!
  sumx=sumx/xnorm
!
  xn(1:matters)=xn(1:matters)/sumx
  yn(1:matters)=yn(1:matters)/sumx
  del=maxval(abs(xn(1:matters)-x(1:matters)))
  x(1:matters)=xn(1:matters)
  y(1:matters)=yn(1:matters)
! ---------------------------------------------------------
  if(nreps.gt.10.or.iter.gt.20) print 47,nreps,iter,tmp,dtime,tim,dt,del,xneg
47 format(' net-step : nreps,iter,tmp,dtime,tim,dt,del,xneg ',2i6,10es12.4)
! ----------------------------------------------------------
  if(abs(xneg).gt.del/2) then
     explic=0.d0
     emplic=1.d0
  end if
! ---------------------------------------------------------
  if(del.lt.1.d-6) go to 400
  if(del.lt.prev_del/2) go to 10
  if_rep=0
  if(iter.gt.20) if_rep=1
  if(iter.gt.10.and.del.gt.0.1) if_rep=1
  if(if_rep.ne.0) then
     dt=dt*min(0.02d0/del,0.5d0)
     y(1:matters)=y0(1:matters)
     go to 9
  endif
  if(iter.gt.100) then
     stop' net_step diverges '
     call flush(6)
  end if
  go to 10
!     -------------------      convergence     -----------------------------
400 continue
  tim=tim+dt
  dt=dt*min(1.2d0,0.1d0/(0.05d0+del))
  if(tim.lt.0.999d0*dtime) go to 7
!     -------------------------------------------------------
500 continue
  qreac=0.d0
  delx=0.d0
  do m=1,matters
     x(m)=y(m)*a(m)
     delx=max(delx,abs(x(m)-y00(m)*a(m)))
     qreac=qreac+(y(m)-y00(m))*excess(m)
  enddo
  qreac=qreac*emev_nuc/dtime
  key_done=1
! --------------------------------------------------------
end subroutine net_step
! ==============================================================================
