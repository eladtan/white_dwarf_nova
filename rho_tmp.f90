subroutine rho_tmp
!
! ---------------------------------------------------------
  use tables
  include 'real8.com'
  parameter (arad=7.564d-15,elmass=9.108d-28, clight=2.99793d10 )
  parameter (planck=6.625d-27,avogadro=6.025d23,boltz=1.380662d-16)
  parameter (third=1.d0/3.d0, xmc=elmass*clight,gascon=avogadro*boltz)
  parameter (xmc2=xmc*clight, echarge=4.803d-10)
  parameter (f43=4.d0/3.d0,pi=3.14159265d0,pai43=f43*pi)
!      parameter (coulomb=-0.9d0*pai43**third*avogadro**f43*echarge**2)
  parameter (coul=-0.9d0*echarge**2)
!   -------------------------------------------------------
  dimension alfat(4),alfas(4),dalfat(4),dalfas(4)
  dimension ix(4),jx(4)
  dimension sp(6),dsp(6),tp(6),dtp(6)
! -------------------------------------------------------
  entropy_const=2.5d0+1.5d0*log(2*pi/avogadro)-log(avogadro) &
       -3*log(planck)
  coulomb=coul*pai43**third*avogadro**f43
  ye=znum/anum
!
  etabl=0.d0
  detdrho=0.d0
  detdtmp=0.d0
  ptabl=0.d0
  dptdrho=0.d0
  dptdtmp=0.d0
  stabl=0.d0
!
  tions=1.d4
  tl=log(tmp)
!   if(tl.lt.ytab(2)) go to 100         not good
!   if(tmp.lt.tions) ye=ye*(tl-ytab(2))/(tions-ytab(2))
  roe=rho*ye
!  ---------------  Ideal Gas for electrons ---------------------
!
  if(log(roe).le.xtab(2)) then
     ptabl   = gascon * roe * tmp
     etabl   =  1.5d0*ptabl
     dptdtmp=  ptabl / tmp
     dptdrho=  ptabl / rho
     detdtmp=  etabl / tmp
     detdrho=  etabl / rho
     go to 100
  endif
! ---------------- Compute table's e and p ------------------------
  eps=1.d-5
  x=min(max(log(roe),xtab(2)+eps*dxtab),xtab(itab-1))
  y=min(max(log(tmp),ytab(2)+eps*dytab),ytab(jtab-1))

!
  ii=min(int(floor((x-xtab(1))/dxtab))+2,itab-1)
  jj=min(int(floor((y-ytab(1))/dytab))+2,jtab-1)
   ii=min(max(ii,3),itab-1)
  jj=min(max(jj,3),jtab-1)
!
  s=(x-(xtab(1)+(ii-2)*dxtab))/dxtab
  t=(y-(ytab(1)+(jj-2)*dytab))/dytab
!xxx  s=(x-xtab(ii-1))/dxtab
!xxx  t=(y-ytab(jj-1))/dytab
  s=max(0.d0,min(s,1.d0))
  t=max(0.d0,min(t,1.d0))
  sp(:)=0.d0
  dsp(:)=0.d0
  tp(:)=0.d0
  dtp(:)=0.d0
  do j=1,ninterp
     if(s.gt.0) then
        sp(j)=s**(j-1)
        dsp(j)=(j-1)*sp(j)/s
     endif
     if(t.gt.0) then
        tp(j)=t**(j-1)
        dtp(j)=(j-1)*tp(j)/t
     endif
  enddo
  sp(1)=1.d0
  tp(1)=1.d0
  dsp(2)=1.d0
  dtp(2)=1.d0
!
  do i=1,4
     alfat(i)=sum(cinterp(i,1:ninterp)*tp(1:ninterp))
     alfas(i)=sum(cinterp(i,1:ninterp)*sp(1:ninterp))
     dalfat(i)=sum(cinterp(i,1:ninterp)*dtp(1:ninterp))
     dalfas(i)=sum(cinterp(i,1:ninterp)*dsp(1:ninterp))
  enddo
!                          e,deds,dedt
  ix(1)=ii-2
  ix(2)=ii-1
  ix(3)=ii
  ix(4)=ii+1
!     ----------------
  jx(1)=jj-2
  jx(2)=jj-1
  jx(3)=jj
  jx(4)=jj+1
!     ----------------
!     find e deds dedt
!     ----------------
  ex=0.d0
  dexds=0.d0
  dexdt=0.d0
  px=0.d0
  dpxds=0.d0
  dpxdt=0.d0
  efx=0.d0
  do i=1,4
     do j=1,4
! -- e
        ex=ex+alfas(i)*alfat(j)*etab(ix(i),jx(j))
        dexds=dexds+dalfas(i)*alfat(j)*etab(ix(i),jx(j))
        dexdt=dexdt+dalfat(j)*alfas(i)*etab(ix(i),jx(j))
! -- p
        px=px+alfas(i)*alfat(j)*ptab(ix(i),jx(j))
        dpxds=dpxds+dalfas(i)*alfat(j)*ptab(ix(i),jx(j))
        dpxdt=dpxdt+dalfat(j)*alfas(i)*ptab(ix(i),jx(j))
! -- ef        
        efx=efx+alfas(i)*alfat(j)*efer(ix(i),jx(j))
     enddo
  enddo
!     ---------------------------------------------------------
  etabl=exp(ex)
  detdrho=etabl/rho*dexds/dxtab
  detdtmp=etabl/tmp*dexdt/dytab
!     ---------------------------------------------------------
  ptabl=px*etabl
  dptdrho=px*detdrho+etabl/rho*dpxds/dxtab
  dptdtmp=px*detdtmp+etabl/tmp*dpxdt/dytab
!
  stabl= (etabl+ptabl)/tmp - gascon*efx*roe
!!!  stabl=max(stabl,0.d0)
  xfermi=efx
! ==================================================================
100 continue
!                                       RADIATION
!
  if(inc_photons.eq.0) then
     prad = 0.d0
  else
     prad = arad / 3 * tmp**4
  endif
  erad = 3*prad
  srad = (erad+prad)/tmp
!                               Ideal Gas for Ions
  if(inc_gas.eq.0) then
     pgas = 0.d0
     sgas = 0.d0
  else
     pgas = gascon * rho * tmp / anum
!  sgas=entropy_const+1.5d0*log(boltz*tmp)-log(rho)+1.5d0*log(anum)
     sgas=entropy_const+1.5d0*log(boltz*tmp)-log(rho)+2.5d0*log(anum)
     sgas=sgas*gascon/anum * rho
  endif
  egas   =  1.5d0*pgas
!
  dpgasdt=  pgas / tmp
  dpgasdd=  pgas / rho
  degasdt=  egas / tmp
  degasdd=  egas / rho
!                                         Coulomb Correction
  rho43=(rho/anum)**f43
  ecol=0.d0
  if(inc_gas.ne.0.and.rho.gt.1.2d+16) ecol=coulomb*znum**2*rho43  !!
  pcol=ecol/3.d0
  decoldd=ecol*f43/rho
  dpcoldd=decoldd/3.d0
! ====================== SCREENING - Debye corrections =================
  e_screen=0.d0
  p_screen=0.d0
  s_screen=0.d0
  dp_screendt=0.d0
  dp_screendd=0.d0
  de_screendt=0.d0
  de_screendd=0.d0
!
  if(key_Coulomb.eq.1) then
     xion=avogadro*rho/anum
     bol_tmp=boltz*tmp
     gam_plasma=(znum*echarge)**2/bol_tmp  * (pai43*xion)**third
     gp=gam_plasma
!xxx     print*,' gam_plasma= ',gam_plasma
     dgdt=-gam_plasma/tmp
     dgdd=+gam_plasma*third/rho
     if(gam_plasma.gt.1.d0) then    !!!  gam_plasma > 1 
        gp4=gp**(0.25d0)
!-----------------------------
        a=-0.9d0
        b=0.97d0
        c=0.22d0
        d=-0.86d0
!-----------------------------
        gpol=(a*gp+b*gp4+c/gp4+d)
        e_screen=xion*bol_tmp * gpol
        de_screendt=e_screen/tmp + xion*bol_tmp * (a+b*gp4/4/gp-c/gp4/4/gp)*dgdt
        de_screendd=e_screen/rho + xion*bol_tmp * (a+b*gp4/4/gp-c/gp4/4/gp)*dgdd
     else                            !!!  gam_plasma < 1 
        a=0.29d0
        b=-0.1d0
        gpol=(a*gam_plasma**(1.5d0)+b*gam_plasma**2)
        e_screen=-3.d0*pgas * gpol
        de_screendt=-3.d0*dpgasdt*gpol - 3.d0*pgas*(a*1.5d0*gp**(0.5d0)+b*2*gp)*dgdt
        de_screendd=-3.d0*dpgasdd*gpol - 3.d0*pgas*(a*1.5d0*gp**(0.5d0)+b*2*gp)*dgdd 
     end if
!xxxxxxxxx
!xxx     print*,' gp,e_screen ',gp,e_screen
!xxxxxxxxxxxxxx
     p_screen=e_screen/3
     dp_screendt=de_screendt/3
     dp_screendd=de_screendd/3
     s_screen=0.d0   !!!   ad-hoc
  end if
! ===============================   TOTALS ==========================================
  prs    = prad + pgas + pcol + ptabl + p_screen
  enr    = ( erad + egas + ecol + etabl + e_screen)/rho
  entropy  = (stabl   + sgas + srad + s_screen )/rho / gascon
  dpdtmp   = 4*prad / tmp + dpgasdt   + dptdtmp + dp_screendt
  dpdrho   = dpgasdd + dpcoldd + dptdrho + dp_screendd
  dedtmp   = (4*erad / tmp   + detdtmp + degasdt + de_screendt)/rho
  dedrho   = (-enr + detdrho + degasdd + decoldd + de_screendd)/rho
!    -------------------------------------------------------------
!        print*,' i,j,x,y ',i,j,x,y
!        print*,' prad,pgas,ptabl ',prad,pgas,ptabl
!        print*,' erad,egas,etabl ',erad,egas,etabl
!     ----------------------------------------------------
end subroutine rho_tmp
! =========================================================================
