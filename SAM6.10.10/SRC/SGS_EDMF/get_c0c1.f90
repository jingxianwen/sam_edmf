
subroutine get_c0c1(cthl,cqt)

!output:
! c0,c1     ... coefficients needed for buoyancy flux


use vars
use params
use grid
use micro_params
use microphysics
use sgs


implicit none

real,dimension(nx,ny,nz),intent(out) :: cthl, cqt

!local variables
integer i,j,k, kb, kc, n,ks
real an, omn, thetali
real :: lambdaf, alphaf, qsl, totheta, tl, qsw, qsi



an = 1./(tbgmax-tbgmin) 

do i=1,nx
do j=1,ny

do k=1,nzm
 
 
  totheta=(pres(k)/1000.)**(rgas/cp)

  !thetali
  thetali = (t(i,j,k)+fac_cond*qpl(i,j,k)+fac_sub*qpi(i,j,k)-gamaz(k))  /  totheta 
  ! get liquid/ice water temperature Tl
  tl     = totheta * thetali

  ! get saturation mixing ratio at Tl
  omn = max(0.,min(1.,(tl-tbgmin)*an))
  qsw=qsatw(tl,pres(k))
  qsi=qsati(tl,pres(k))
  qsl =  omn * qsw + (1.-omn) * qsi  
  ! get dqsdT
  alphaf = (omn * qsw * lcond/(tl**2*rv) + (1.-omn) * qsi * lsub/(tl**2*rv))
  lambdaf = (1. + alphaf * (omn*fac_cond + (1.-omn) * fac_sub) )**(-1.)
  alphaf =alphaf * totheta

  ! Bechtold (1995) Eqn. (B5)
  cthl(i,j,k) = cfrac_pdf(i,j,k) * (1.+epsv*q(i,j,k)-(fac_cond/totheta - (1.+epsv) * &
         tabs(i,j,k)/totheta)*lambdaf*alphaf) + (1.-cfrac_pdf(i,j,k)) * (1.+epsv*q(i,j,k))
  cqt(i,j,k)  = cfrac_pdf(i,j,k) * (epsv*tabs(i,j,k)/totheta+(fac_cond/totheta - (1.+epsv) *tabs(i,j,k)/totheta)*lambdaf ) &
           + (1.-cfrac_pdf(i,j,k)) * epsv * tabs(i,j,k)/totheta

end do


end do ! j
end do ! i

 
end subroutine get_c0c1

