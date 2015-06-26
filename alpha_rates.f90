subroutine alpha_rates (matters,rho,tmp,y,dydt,dyidyj)
!
  include 'real8.com'
  include 'network.com'
  dimension y(matters),dydt(matters),dyidyj(matters,matters+1)
!     ----------------------------------------------------------------
!   y(i)=x(i)/a(i)  =  number abundances
!   dydt(i)= rate(rho,tmp,y(1,,,matters)
! rrrrrrrrrrrrrrrrrrrrrrrrr  compute rates  rrrrrrrrrrrrrrrrrrrrr
!         print*,' alphanet-matters,rho,tmp ',matters,rho,tmp
!
  dydt(:)=0.d0
  dyidyj(:,:)=0.d0
  do 100 nr=1,numreac
     rate=rates(nr)
!
     if(npart(nr).eq.1) then
!     one particle reaction ( beta decay)
        m=inpt(nr,1)
        ratey=y(m)*rate
        dydt(m)=dydt(m)-ratey
        dyidyj(m,m)=dyidyj(m,m)+rate
        do j=1,3
           ij=iout(nr,j)
           if(ij.gt.0) then
              dydt(ij)=dydt(ij)+ratey
              dyidyj(ij,m)=dyidyj(ij,m)-rate
           endif
        enddo
!     two particeles reactions
     else if (npart(nr).eq.2) then
!
        rtro=rate*rho
        m=inpt(nr,1)
        n=inpt(nr,2)
        ratem=rtro*y(m)
        raten=rtro*y(n)
        ratemn=rtro*y(m)*y(n)
        if(m.eq.n) ratemn=ratemn/2
        dydt(m)=dydt(m)-ratemn
        dydt(n)=dydt(n)-ratemn
        dyidyj(m,m)=dyidyj(m,m)+raten
        dyidyj(n,n)=dyidyj(n,n)+ratem
        dyidyj(m,n)=dyidyj(m,n)+ratem
        dyidyj(n,m)=dyidyj(n,m)+raten
        do j=1,3
           ij=iout(nr,j)
           if(ij.gt.0) then
              dydt(ij)=dydt(ij)+ratemn
              dyidyj(ij,m)=dyidyj(ij,m)-raten
              dyidyj(ij,n)=dyidyj(ij,n)-ratem
           endif
        enddo
!     three particles   reactions    (3*alfa-c only)
     else
        if(npart(nr).ne.3) stop ' npart err in rates '
        rtro=rate*rho**2/6
        rtroy=rtro*y(1)**3
        dydt(1)=dydt(1)-3*rtroy
        dydt(2)=dydt(2)+rtroy
        dyidyj(1,1)=dyidyj(1,1)+rtro*9.d0*y(1)**2
        dyidyj(2,1)=dyidyj(2,1)-rtro*3.d0*y(1)**2
     endif
!
100 enddo
!xxxxxxxxxxxxxxxxxxx
  keypr=0
  if(keypr.gt.0) then
     open(13,file='alpha.out',status='unknown',form='formatted')
     do i=1,matters
        write(13,1313) i,(dyidyj(i,j),j=1,matters)
1313    format(' i=',i3,13es9.1)
     enddo
     stop 'keypr in alpha_rates'
  endif
!xxxxxxxxxxxxxxxxxxx
! --------------------------------------------------------
end subroutine alpha_rates
