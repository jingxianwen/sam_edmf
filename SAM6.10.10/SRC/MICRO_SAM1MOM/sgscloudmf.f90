
subroutine sgscloudmf()


!output:
! cfrac_pdf ... stratiform cloud cover 
! qn        ... domain mean non-precip cloud condensate
! qv        ... domain mean watervapor
! qcl       ... domain mean liquid condensate
! qci       ... domain mean frozen condensate
! tabs      ... domain mean absolute temperature
! c0,c1     ... coefficients needed for buoyancy flux


use vars
use params
use grid
use micro_params
use microphysics, only : q,qp,qn


implicit none

!local variables
integer i,j,k, kb, kc, n,ks
real dtabs, an, omn
real, dimension(nzm) :: qte
real  qne, tabse, ds, q1
integer,parameter :: niter=10
real :: lambdaf, alphaf, qsl, totheta, tl, qsw, qsi, tabsnoql, dqndt, domndt, frac_mf2
real,dimension (2) :: tabs1, tabs2

real, parameter :: q_crit=1.6
real, parameter :: lcld = 150.

real :: cab,ckk,leps
real,dimension(nzm) :: qtqt, thlthl, qtthl, thetali, thetaligrad, qtgrad, sigmas

an = 1./(tbgmax-tbgmin) 

!IF (tke(i,j,k).gt.1.d-3.and.dosgscloud) then
IF (dosgscloud) then


do i=1,nx
do j=1,ny

do k=1,nzm
  totheta=(pres(k)/pres0)**(rgas/cp)
  !total water
  qte(k) = q(i,j,k)
  !thetali
  thetali(k) = (t(i,j,k)+fac_cond*qpl(i,j,k)+fac_sub*qpi(i,j,k)-gamaz(k))  /  totheta 
end do

do k=1,nzm


  
! compute vertical gradients of qt and thetali
  kb=k-1
  kc=k+1
  if(k.eq.1) then
   kb=1
   kc=2
  elseif (k.eq.nzm) then
   kb=nzm-1
   kc=nzm
  end if

  thetaligrad(k) = (thetali(kc)-thetali(kb))/(z(kc)-z(kb))
  qtgrad(k)      = (qte(kc)-qte(kb))/(z(kc)-z(kb))

! use balance between dissipation and production
! to diagnose (co-)variances, then apply limiters
  qtqt(k)  = lcld**2*qtgrad(k)**2. 
  thlthl(k)= lcld**2*thetaligrad(k)**2.
  qtthl(k) = lcld**2*qtgrad(k)*thetaligrad(k)

  ! limit (co-)variances
  ! Kay's limiters
  !thlthl(k)=max(min(thlthl(k),10.),0.01)
  !qtqt(k)=max(min(qtqt(k),1.e-5),1.e-8)
  !qtthl(k)=sign(max(min(abs(qtthl(k)),1.e-2),1.e-5),qtthl(k))
  ! max/min taken from Heinze's (2005, JAMES) LES
  thlthl(k)=max(min(thlthl(k),3.),0.001)
  qtqt(k)=max(min(qtqt(k),2.e-6),1.e-8)
  qtthl(k)=max(min(qtthl(k),0.02d-3),-4.d-3)

 ! get liquid/ice water temperature Tl
 totheta=(pres(k)/pres0)**(rgas/cp)
 tl     = totheta * thetali(k)
 tabsnoql = tl ! save tabs assuming no condensate (which equals Tl in that case)
 tabse=tabsnoql! set initial guess for temperature to temp without condensate

 ! get saturation mixing ratio at Tl
 omn = max(0.,min(1.,(tl-tbgmin)*an))
 qsw=qsatw(tl,pres(k))
 qsi=qsati(tl,pres(k))
 qsl =  omn * qsw + (1.-omn) * qsi  
 ! get dqsdT
 alphaf = (omn * qsw * lcond/(tl**2*rv) + (1.-omn) * qsi * lsub/(tl**2*rv))
 lambdaf = (1. + alphaf * (omn*fac_cond + (1.-omn) * fac_sub) )**(-1.)
 alphaf =alphaf * totheta

 ! mean saturation deficit
 ds=lambdaf *(qte(k)-qsl)
 ! and its std dev
 sigmas(k) = lambdaf * (max(0.0,qtqt(k) + alphaf**2 * thlthl(k)**2. - 2. * alphaf * qtthl(k)))**(0.5)
 if (dozerosigma) sigmas(k) = 0.0

 ! The scheme is actually diagnostic, but we find that the initial guess can be slightly improved
 ! (ie the Taylor approximation to the mean sat. deficit can be improved)
 ! That way, for sigma=0, the scheme gives exactly the same result as cloud() 
 ! only difference, here we don't consider the enthalpy change due to a repartitioning of qpl and qpi, 
 ! ie thetali is defined based on non-falling condensate only
 n=0
 dtabs=100.
 do while (abs(dtabs).gt.0.01.and.n.le.niter) 

 if (n.gt.0) then ! ds is known for first guess, but is recomputed otherwise based on new tabs
   tabse  = tabse+dtabs
   omn = max(0.,min(1.,(tabse-tbgmin)*an))
   qsw=qsatw(tabse,pres(k))
   qsi=qsati(tabse,pres(k))
   qsl = omn * qsw + (1.-omn) * qsi
   ds=qte(k)-qsl
 end if
 
 if (tabse.gt.tbgmax.or.tabse.lt.tbgmin) then  
    domndt = 0.
 else
    domndt  = 1./an
 end if

 IF (sigmas(k).le.0.0) then
   ! homogeneous grid box
   if (donoenvcloud) then
   cfrac_pdf(i,j,k) = 0.
   qne = 0.
   else
   cfrac_pdf(i,j,k) =   ABS ( (SIGN(1.0,ds)+1.0)*0.5 )
   qne = cfrac_pdf(i,j,k) * ds
   end if
   q1=-999.
   dqndt =  - (omn  * dtqsatw(tabse,pres(k)) + (1.-omn) * dtqsati(tabse,pres(k)))
 ELSE
   q1= ds/sigmas(k)
   cfrac_pdf(i,j,k) = MIN ( 1.0, MAX ( 0.0, &
                                        0.5 * (1.0+q1/q_crit) ) )
   IF ( q1 .le. - q_crit ) THEN
      qne = 0.0
      dqndt = 0.
   ELSEIF ( q1 .ge. q_crit ) THEN
      qne = sigmas(k) * q1 
      dqndt = - (omn  * dtqsatw(tabse,pres(k)) + (1.-omn) * dtqsati(tabse,pres(k)))
   ELSE
      qne = sigmas(k) * (q1+q_crit) * (q1+q_crit) / (2.*(q_crit+q_crit))
      dqndt = - (q1+q_crit)/(2.*q_crit) * (omn  * dtqsatw(tabse,pres(k)) + (1.-omn) * dtqsati(tabse,pres(k)))
   ENDIF
 END IF
 
 
 ! compute missing tabs increment (considering both errors in tabs, qn, and partitioning into ql/qi)
 ! the derivative of qn with respect to tabs is based on the PDF parameterization (instead of on dqs/dtabs as in cloud())
 dtabs= (tabsnoql - tabse + (omn*fac_cond+(1.-omn)*fac_sub) * qne) / (1.- (omn*fac_cond+(1.-omn)+fac_sub) * dqndt + qne * fac_fus * domndt )

 n=n+1

 end do ! end of iterative loop

 frac_mf2 = 0.5*(frac_mf(i,j,k)+frac_mf(i,j,k+1))
 qcl(i,j,k) = (1.-frac_mf2) * omn * qne + 0.5 * frac_mf2 * (qcsgs_mf(i,j,k)+qcsgs_mf(i,j,k+1))
 qci(i,j,k) = (1.-frac_mf2) * (1.-omn) * qne + 0.5 * frac_mf2 * (qisgs_mf(i,j,k)+qisgs_mf(i,j,k+1)) 
 cfrac_pdf(k) = min(frac_mf2,0.5*(cfrac_mf(i,j,k+1)+cfrac_mf(i,j,k))) + (1.-frac_mf2) * cfrac_pdf(i,j,k)

 tabs(i,j,k)  = t(i,j,k)-gamaz(k) + fac_cond*(qpl(i,j,k)+qcl(i,j,k)) + fac_sub*(qpi(i,j,k)+qci(i,j,k)) 

 qn(i,j,k) = max(qcl(i,j,k) + qci(i,j,k),0.)
 qv(i,j,k) = max(0.,q(i,j,k) - qn(i,j,k))
 qp(i,j,k) = max(0.,qp(i,j,k))


end do ! k
end do ! j
end do ! i

else

 call cloud()

end if ! dosgscloud


end subroutine sgscloudmf

