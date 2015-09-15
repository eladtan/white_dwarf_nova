subroutine eos_fermi (keyeos,im_gas,im_photons,im_Coulomb   &
     ,rho_in,enr_in,tmp_in,prs_in,ent_in,abar,zbar          &
     ,chem_pot,dpdro,dpde,dedro,dedt,sounds,keyerror)
!
  use tables
  !  include 'real8.com'
  integer, intent(in) :: keyeos
  integer, intent(in) :: im_gas
  integer, intent(in) :: im_photons
  integer, intent(in) :: im_Coulomb
  real*8, intent(inout) :: rho_in
  real*8, intent(inout) :: enr_in
  real*8, intent(inout) :: tmp_in
  real*8, intent(inout) :: prs_in
  real*8, intent(inout) :: ent_in
  real*8, intent(in) :: abar
  real*8, intent(in) :: zbar
  real*8, intent(out) :: chem_pot
  real*8, intent(out) :: dpdro
  real*8, intent(out) :: dpde
  real*8, intent(out) :: dedro
  real*8, intent(out) :: dedt
  real*8, intent(out) :: sounds
  integer, intent(out) :: keyerror
  ! ------------------------------------------------------------------
  rho0=rho_in
  enr0=enr_in
  tmp0=tmp_in
  prs0=prs_in
  rho =rho_in
  enr =enr_in
  tmp =tmp_in
  anum=abar
  znum=zbar
  inc_gas=im_gas
  inc_photons=im_photons
  key_Coulomb=im_Coulomb
! ---------------------------------------------------------------
  keyerror=0
  keyerr=0
  select case(keyeos)
  case(0)
     call rho_enr
  case(1)
     call rho_tmp
  case(2)
     call rho_prs
  case(4)
     call tmp_prs
  case(5)
     call enr_prs
  case default
     print*,' wrong keyeos in eos_tabuler - keyeos=',keyeos
     stop ' STOP wrong keyeos in eos_tabuler '
  end select
  keyerror=keyerr
! ---------------------------------------------
  rho_in=rho
  tmp_in=tmp
  enr_in=enr
  prs_in=prs
  ent_in=entropy
  chem_pot=xfermi
!
  dedtmp=max(dedtmp,1.d-10)   !!!!!!!!!!!
  dedt=dedtmp
  dpde = dpdtmp /dedtmp
  dpdro= dpdrho
  dedro= dedrho
!         dpdroe - derivative for fixed e
  dpdroe= dpdrho-dpde*dedrho
  if(keyeos.eq.0) dpdro=dpdroe
!
  sounds=max(dpdroe+dpde*prs/rho**2,1.d-30)
  sounds=sqrt(sounds)
!  --------------------------------------------------------
end subroutine eos_fermi
