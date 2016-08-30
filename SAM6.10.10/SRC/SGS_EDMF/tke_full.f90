
subroutine tke_full

!	this subroutine solves the TKE equation

use vars
use sgs
use params
use microphysics, only : qp
implicit none

real def2(nx,ny,nzm)
real grd,betdz,Ck,smix,Pr,Cee,pblh,thetav_c,thetav_b,thetav_k
real buoy_sgs,a_prod_sh,a_prod_bu,a_diss
real lstarn, lstarp, bbb, omn, omp
real qsatt,dqsat
integer i,j,k,kc,kb
real thetavs, sgs_thv_flux, tketau, l23, tke_thvflx, wstar
real dtkedtsum, dtkedtmin
real, parameter :: xkar=0.4
real wthl, wqt

real tabs_interface, qp_interface, qtot_interface, qsat_check, qctot

call t_startf('tke_full')

Ck=0.5
Cee= 0.16 * 2.5
Pr=1. 
pblh=1000.

thvflx = 0.

if(RUN3D) then
  call shear_prod3D(def2)
else
  call shear_prod2D(def2)
endif

do k=1,nzm
  kb=k-1
  kc=k+1

  grd=dz*adz(k)

  if(k.eq.1) then
    kb=1
    kc=2
  end if
  if(k.eq.nzm) then
    kb=nzm-1
    kc=nzm
  end if

  do j=1,ny
  do i=1,nx

  tke(i,j,k)=max(0.,tke(i,j,k))
  !N**2 based on virtual potential temperature
  thetav_c = (1.+epsv*qv(i,j,kc))*tabs(i,j,kc)*(pres0/pres(kc))**(rgas/cp)
  thetav_k = (1.+epsv*qv(i,j,k))*tabs(i,j,k)*(pres0/pres(k))**(rgas/cp)
  thetav_b = (1.+epsv*qv(i,j,kb))*tabs(i,j,kb)*(pres0/pres(kb))**(rgas/cp)
  buoy_sgs=ggr/thetav_k * (thetav_c-thetav_b)/ (z(kc)-z(kb))

  !thetavs
  thetavs = (1.+epsv*qv(i,j,1))*tabs(i,j,1)*(pres0/pres(1))**(rgas/cp)
  sgs_thv_flux = (1.+epsv*qv(i,j,1))*fluxbt(i,j) + epsv*tabs(i,j,1)*(pres0/pres(1))**(rgas/cp)*fluxbq(i,j)
  tketau= max(0.5 * pblh /  ((ggr/thetavs*sgs_thv_flux*pblh)**(1./3.)),0.0)
  l23 = (tketau*sqrt(tke(i,j,k)))**(-1)
  l23 = l23**(-1)
  smix=  l23 + (xkar*z(k)-l23)*exp(-z(k)/100.)
  tk(i,j,k) = Ck*smix*sqrt(tke(i,j,k))


  ! buoyancy flux
  wthl = (0.5*(twsb3(i,j,k) + twsb3(i,j,k+1)))* (pres0/pres(k))**(rgas/cp) 
  if (qp(i,j,k).gt.0.0) then
    wthl = wthl +  ((fac_cond*qpl(i,j,k)+fac_sub*qpi(i,j,k))/qp(i,j,k) * 0.5*(mkwsb3(i,j,k,2) + mkwsb3(i,j,k+1,2)))  &
           * (pres0/pres(k))**(rgas/cp) 
  end if
  wqt  = 0.5*(mkwsb3(i,j,k,1) + mkwsb3(i,j,k+1,1))

  ! get buoyancy flux from PDF scheme
  !thvflx(i,j,k) = cthl(i,j,k) * wthl + cqt(i,j,k) * wqt
  ! assume fully unsaturated for now
  thvflx(i,j,k) = (1.+epsv*qv(i,j,k)) *  wthl + epsv * thetav_k/(1.+epsv*qv(i,j,k)) * wqt 
 
  a_prod_sh=(tk(i,j,k)+0.001)*def2(i,j,k)
  a_prod_bu= ggr/thetav_k * thvflx(i,j,k)
  a_diss=Cee / smix*tke(i,j,k)**1.5 
  dtkedtsum = a_prod_sh+a_prod_bu-a_diss
  dtkedtmin = -tke(i,j,k)/dtn
  if (dtkedtsum.lt.dtkedtmin) then
    dtkedtsum = dtkedtmin / dtkedtsum
    a_prod_bu = dtkedtsum * a_prod_bu
    a_prod_sh = dtkedtsum * a_prod_sh
    a_diss    = dtkedtsum * a_diss
  end if
  tke(i,j,k)=tke(i,j,k)+dtn*(a_prod_sh+a_prod_bu-a_diss)
  wstar=max(0.,(ggr/thetavs*sgs_thv_flux*pblh)**(1./3.))

  tkelediss(k) = 0.
  tkesbdiss(k) = 0.
  tkesbshear(k)= 0.
  tkesbbuoy(k) = 0.

  tkh(i,j,k)=Pr*tk(i,j,k)

  tkelediss(k) = tkelediss(k) - a_prod_sh
  tkesbdiss(k) = tkesbdiss(k) + a_diss
  tkesbshear(k)= tkesbshear(k)+ a_prod_sh
  tkesbbuoy(k) = tkesbbuoy(k) + a_prod_bu

  end do ! i
  end do ! j

  tkelediss(k) = tkelediss(k)/float(nx*ny)


end do ! k

call t_stopf('tke_full')

end


