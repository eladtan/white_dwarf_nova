module global
  real(8),save :: xnorm
end module global
! =================================================================
subroutine burn_step (indxeos,rho,enr,tmp,x,amol,zmol,dedtmp  &
     ,matters,dtime,qreac,nse,tmpsaf_nse,key_done,screen_type)
!
  use global
  include 'real8.com'
  include 'network.com'
  dimension x(matters),xn(matters),xnse(matters)
  character(*) :: screen_type
!  ------------------------------------------------------------
  xnorm=sum(x(1:matters))
!  ------------------------------------------------------------
  key_done=0
  if(nse.ne.0.and.tmp.gt.tmpsaf_nse) then
     ense=enr
     tmpx_nse=tmp
     xnse(:)=x(:)
!
     call nse_equilib (indxeos,rho,tmpx_nse,ense,xnse,amater,matters  &
          ,amol,zmol,de_nse,tmpsaf_nse,key_done)
     if(key_done.eq.0) then
        print*,' nse_equlib fails : tmp,tmpx_nse,tmpsaf_nse ',tmp,tmpx_nse,tmpsaf_nse
        go to 100
     end if
     x(:)=xnse(:)
     qreac=(ense-enr)/dtime
     key_done=2
     return
  endif
!  ------------  rate integration ----------------------------
100 continue
  xn(:)=x(:)
  call net_step (rho,tmp,xn,amater,zmater,dedtmp,matters,dtime,qreac,delx,key_done)
!
  if(key_done.ne.0) then
     x(:)=xn(:)
  end if
end subroutine burn_step
! ====================================================================================
