      parameter (maxreac=100,maxmat=15)
!      ---------------------------------------------------------------
      parameter(shlish=1.d0/3.d0,half=0.5d0,pai=3.141592736d0)
      parameter(pai2=2.d0*pai,pai4=2.d0*pai2,pai8=pai*8.d0)
      parameter(avo=6.022045d23,boltz=1.380662d-16,charge=4.80286d-10)
      parameter(clite=2.997925d10,emev=1.60219d-6,fmelct=9.10953d-28)
      parameter(planck=6.62618d-27,tmev=emev/boltz,rgas=avo*boltz)
      parameter(elect=fmelct*clite**2,compton=planck/(fmelct*clite))
      parameter(emev_nuc=9.648846d17)
!     ------------------------------------------------------------------
      save /network/
      common/network/numreac,nreac,nmat,if_screen       &
      ,npart(maxreac),inpt(maxreac,3),iout(maxreac,3)
      save /worknet/
      common/worknet/crate(7,maxreac)                   &
      ,amater(maxmat),zmater(maxmat),excess(maxmat)     &
      ,rates(maxreac),drate_dtmp(maxreac),tmp_nse
      save /cnet/
      common/cnet/xmater(maxmat)
      character*5 xmater
! ----------------------------------------------------------------------
