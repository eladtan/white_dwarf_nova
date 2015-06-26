subroutine sigmav (t9,t913,t953,t9log,nr,rate,drdtmp)
!                             by theilemann , 1989
  include 'real8.com'
  include 'network.com'
!  ------------------------------------------------------------
  rate=1.d-100
  drdt=0.d0
  rate=crate(1,nr)+crate(2,nr)/t9+crate(3,nr)/t913+crate(4,nr)*t913 &
       +crate(5,nr)*t9+crate(6,nr)*t953+crate(7,nr)*t9log
  rate=max(rate,-200.d0)
  if(rate.gt.600.d0) write(*,*) ' nr,t9,r ',nr,t9,rate
  rate=min(rate,600.d0)
  rate=exp(rate)
  drdtmp=rate*(-crate(2,nr)/t9**2-crate(3,nr)/3.d0/t913/t9 &
       +crate(4,nr)/3.d0*t913/t9+crate(5,nr) &
       +5.d0/3.d0*crate(6,nr)*t913**2+crate(7,nr)/t9) /1.d9
!
end subroutine sigmav
