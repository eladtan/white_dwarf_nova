subroutine linsys_burn (a,b,x,key,neq,mdim)
  include 'real8.com'
  dimension a(mdim,mdim),b(mdim),x(mdim)
!  ---------------------------------------------------
  if(neq.le.1) then
     x(1)=b(1)/a(1,1)
     return
  endif
!  -------------------------------------------------
  if(key.eq.2) go to 100
!                            Gauss elimination of a
  do  i=1,neq-1
     do  j=i+1,neq
        p=a(j,i)/a(i,i)
        if(abs(p).gt.1.d-12) then
           a(j,i)=p
           do k=i+1,neq
              a(j,k)=a(j,k)-p*a(i,k)
           end do
        endif
     end do
  end do
  if (key.eq.1) return
!     ----------------------------------------------
100 continue
!                           Gauss elimination of b
  do i=1,neq-1
     do j=i+1,neq
        b(j)=b(j)-a(j,i)*b(i)
     end do
  end do
! ---------------------  back substitution
  x(neq)=b(neq)/a(neq,neq)
  do ii=2,neq
     i=neq+1-ii
     s=b(i)
     do j=i+1,neq
        s=s-a(i,j)*x(j)
     end do
     x(i)=s/a(i,i)
  end do
!--------------------------------------------------------------
end subroutine linsys_burn
