module tables
  real(8),allocatable,save :: xtab(:),ytab(:)          &
       ,etab(:,:),ptab(:,:),efer(:,:)
  real(8),save :: dxtab,dytab,diuk=1.d-9
  real(8),save :: rho0,tmp0,enr0,prs0,ent0             & 
       ,rho,tmp,enr,prs,entropy,xfermi                 &
       ,anum,znum,dedrho,dedtmp,dpdrho,dpdtmp          &
       ,tmp_min,tmp_max
!
  integer,save :: init=0,ninterp=4,itmax=100,itab,jtab &
       ,inc_gas=1,inc_photons=1,keyerr,key_Coulomb=1
  real(8),save :: cinterp(4,5)
! =========================================================================
contains
  subroutine rd_tables (tab_file)
! --------------------------------------------
    include 'real8.com'
    character*80 header,label
    character*(*) tab_file
! --------------------------------------------------------------------------
!                define interpolation constants
!      
    nint=ninterp
    if(nint.eq.4) then
!       3trd order- 4 coeeficients
       cinterp(1,1:nint)= (/0.d0,-0.5d0, 1.0d0,-0.5d0/)
       cinterp(2,1:nint)= (/1.d0,0.d0  ,-2.5d0, 1.5d0/)
       cinterp(3,1:nint)= (/0.d0,0.5d0 ,2.d0  ,-1.5d0/)
       cinterp(4,1:nint)= (/0.d0,0.d0  ,-0.5d0, 0.5d0/)
    else if (nint.eq.5) then
!       4th order- 5 coeeficients
       cinterp(1,1:nint)= (/0.d0,-0.5d0, 2.0d0,-2.5d0, 1.d0/)
       cinterp(2,1:nint)= (/1.d0,0.d0  ,-3.5d0, 3.5d0,-1.d0/)
       cinterp(3,1:nint)= (/0.d0,0.5d0 ,1.d0  ,0.5d0 ,-1.d0/)
       cinterp(4,1:nint)= (/0.d0,0.d0  ,0.5d0 ,-1.5d0, 1.d0/)
    endif
! --------------------------------------------------------------------------
    nt=77
    ibin=0        !!!
    if(ibin.eq.0) then
!     ------------------------ Coded File --------------------------------
!       print*,' open coded table file in rd_tables -',tab_file
       open (nt,file=tab_file,status='old')
       read(nt,*) header
!       write (* ,111) header
111    format(1x,a80)
       read(nt,*) itab,jtab,dxtab,dytab
!       write(*,124) itab,jtab,dxtab,dytab
124    format(' itab,jtab,dxtab,dytab ',2i4,1p2e12.4)
!
       allocate(xtab(itab),ytab(jtab)  &
            ,etab(itab,jtab),ptab(itab,jtab),efer(itab,jtab))
!         
!       read(nt,*) xtab(1:itab)
!       read(nt,*) ytab(1:jtab)
!       read(nt,*) etab(1:itab,1:jtab)
!       read(nt,*) ptab(1:itab,1:jtab)
!       read(nt,*) efer(1:itab,1:jtab)
       read (nt,133) (xtab(i),i=1,itab)
       read (nt,133) (ytab(j),j=1,jtab)
       read (nt,133) ((etab(i,j),i=1,itab),j=1,jtab)
       read (nt,133) ((ptab(i,j),i=1,itab),j=1,jtab)
       read (nt,133) ((efer(i,j),i=1,itab),j=1,jtab)
133    format(2x,8es15.7)
    else
!    --------------------- Binary File ------------------------
!       print*,' open binary table file in rd_tables -',tab_file
       open (nt,file=tab_file,status='old',form='unformatted')
!       print*,' open succeed'
       read (nt) header
       read (nt) itab,jtab,dxtab,dytab
!       print*,' itab,jtab ',itab,jtab
!
       allocate(xtab(itab),ytab(jtab)  &
            ,etab(itab,jtab),ptab(itab,jtab),efer(itab,jtab))
!         
       read (nt) xtab   !!! (1:itab)
       read (nt) ytab   !!! (1:jtab)
       read (nt) etab   !!! (1:itab,1:jtab)
       read (nt) ptab   !!! (1:itab,1:jtab)
       read (nt) efer   !!! (1:itab,1:jtab)
    endif
    close (nt)
! -----------------------------------------------------------------
    tmp_min=exp(ytab(1))
    tmp_max=exp(ytab(jtab))
!    write(*,144) tmp_min,tmp_max
144 format(' EOS- tmp_min,tmp_max ',2es12.4)
  end subroutine rd_tables
end module tables
! ========================================================================
subroutine init_tabular (tab_file)
  use tables
  include 'real8.com'
  !character*(*) tab_file
  character(80) :: tab_file
  call rd_tables (tab_file)
end subroutine init_tabular
! ========================================================================
