subroutine sigmav_rates (matters,rho,tmp,aa,zz)
!
  include 'real8.com'
  include 'network.com'
!
  rates(:)=0.d0
  drate_dtmp(:)=0.d0
!      drate(:)=0.d0
  e_screen=1.d0
!
  t9=min(tmp/1.d9,20.d0)
  if(t9.lt.0.05) return
  t913=t9**(1./3.)
  t923=t913**2
  t932=t9**1.5
  t953=t913**2*t9
  t9log=log(t9)
!
!
  do  nr=1,numreac
!
     call sigmav(t9,t913,t953,t9log,nr,rates(nr),drate_dtmp(nr))
!               print*,' nr,rate ',nr,rates(nr)   !!
!
     e_screen=1.d0
     if (if_screen.eq.2) then
        if (npart(nr).eq.2) then
           nin=2
           z1=zmater(inpt(nr,1))
           z2=zmater(inpt(nr,2))
           call screening (if_screen,e_screen     &
                ,rho,tmp,aa,zz,nin,z1,z2)
        else if(npart(nr).eq.3) then
           nin=3
           z1=zmater(1)
           z2=zmater(1)
           call screening (if_screen,e_screen     &
                ,rho,tmp,aa,zz,nin,z1,z2)
        endif
     endif
     rates(nr)=rates(nr)*e_screen
     drate_dtmp(nr)=drate_dtmp(nr)*e_screen
!
  enddo
!
end subroutine sigmav_rates
