! WL this routine provides the fluxes sumM(nz), sumMrt(nz), sumMu(nz), sumMv(nz), sumMt(nz), and sumMtke(nz)
! WL that result from the multiplume model provided by Kay Suselj and adapted by WL
subroutine edmf()

use vars
use grid
use params
implicit none

! ==============================================================================
! Multiplume stochastic EDMF model (some comments by Kay)
!  - does vertical  integration for the single horizontal grid-point (loop over horizontal points must be outside of this routine) 
!  - nup plumes initialized at the surface, integrated until the level where their vertical velocity becomes zero
!  - surface conditions - following Cheinet 2003
!  - plume integration - Suselj et al., 2013 and Suselj et al., 2014
!  - entrainment rate - stochastic process, following Poisson distribution
!
! output needed for tendencies (on faces)
! sumM = sum (a_i x w_i) [ms-1]
! sumMt = sum (a_i x w_i x t_i) [Kms-1]
! sumMrt = sum (a_i x w_i x rt_i)  [ms-1]
! sumMu = sum (a_i x w_i x u_i)     = 0. [ms-1]^2
! sumMv = sum (a_i x w_i x v_i)     = 0. [ms-1]^2
! sumMtke = sum (a_i x w_i x tke_i) = 0. [ms-1]^3
!
! Input needed
! zw ... heights of the half levels
! p,u,v,thl,thv,qt (kts:kte) - profiles of grid-box mean variables 
! ust [ms-1],wthl [Kms-1],wqt [ms-1] - surface fluxes (ustar, sensible and latent heat)
! pblh - boundary layer height
! ==============================================================================
! local variables:
 ! entrainment variables     
      REAl,DIMENSION(1:nzm,1:nup)    :: ENTf
      INTEGER,DIMENSION(1:nzm,1:nup) :: ENTi

      INTEGER :: K,I
      REAL :: wthv,wqt,wthl,qstar,thstar,sigmaW,sigmaQT,sigmaTH,sigmaTHV,zs, &
           pwmax,wmin,wmax,wlv,wtv
      REAL :: QTn,Tn,THVn,QCLn,QCIn,Un,Vn,Wn2,EntEXP,EntW, hlp, acrit, Wa

! w parameters
      REAL,PARAMETER :: &
        ! witek
        !&Wa=2.,& 
        !&Wb= 0.0 , & 
        !&Wc= 1.0 
        ! de roode
     !   &Wa=0.5,& 
     Wb= 0.0,&  
     Wc= 0.9
        ! suselj
        !&Wa=2./3.,& 
        !&Wb= 0.002 , & 
        !&Wc= 1.5 

! entrainment parameters
      REAL,PARAMETER :: &
        & L0=100.,&
        & ENT0=0.1,&
        & deldz=2.e-3/1000. ! increase in detrainment above cloud base

! termination fractional criterion: exit, if ai < facrit * ai0
      REAL,PARAMETER :: &
        & facrit=1.e-4


!initialize plume properties
 UPM=0.
 UPW=0.
 UPT=0.
 UPTABS=0.
 UPTHV=0.
 UPQT=0.
 UPQCL=0.
 UPQCI=0.
 UPA=0.
 UPCF=0.
 UPU=0.
 UPV=0.
 ENT=0.
 DET=0.
 BUOY=0.
 UPTHD=0.
 qcsgs_mf=0.
 qisgs_mf=0.
 qtsgs_mf=0.
 cfrac_mf=0.
 tabs_mf=0.
 tke_mf=0.
 frac_mf=0.


 
 ! surface fluxes
 ! sensible heat flux (K m s-1)
 wthl = fluxbt(i,j,)
 ! latent heat flux (m s-1)
 wqt  = fluxbq(i,j) 
 ! virtual pot. temp flux
 wthv = (1.+epsv*qv(i,j,1))*fluxbt(i,j) + epsv*tabs(i,j,1)*(pres0/pres(1))**(rgas/cp)*fluxbq(i,j)

 ! quit in case of a non-positive buoyancy flux
 if (wthv.le.0.0) return
 
 Wa = 0.4* Wc + 0.3
 ! get entrainment rate
  !do i=1,nup
  !  do k=1,nzm
  !    ENTf(k,i)=(zi(k+1)-zi(k))/L0
  !   enddo
  ! enddo

 ! get random Poisson number
   !call Poisson(1,nup,1,nzm,ENTf,ENTi)

 ! Calculate entrainment rate: Ent=Ent0/dz*P(dz/L0) for each layer             
  !  do i=1,nup
  !  do k=1,nzm
  !    if (fixedeps) then
  !      ENT(k,i)= eps0  
  !    else
  !      ENT(k,i)=real(ENTi(k,i))*Ent0/(zi(k+1)-zi(k))
  !    end if
  !  enddo
  !  enddo



! set initial conditions for updrafts
    zs=50.
    pwmax=3.5

! see Lenschow et al. (1980), JAS
    wstar=max(0.,(ggr/thetav(1)*wthv*pblh)**(1./3.))
    qstar=wqt/wstar
    thstar=wthl/wstar 
    sigmaW=1.34*wstar*(zs/pblh)**(1./3.)*(1.-0.8*zs/pblh)
    sigmaQT=1.34*qstar*(zs/pblh)**(-1./3.)
    sigmaTH=1.34*thstar*(zs/pblh)**(-1./3.)
 
! now get sigmaTHV from linearization of pot. temp. and using a corr. coeff between qv and theta of 0.75 (see Sobjan 1991)
    sigmaTHV=sqrt(sigmaTH**2+epsv**2*theta(1)**2*sigmaQT**2 &
        + 2.*epsv*theta(1)*0.75*sigmaTH*sigmaQT)

    wmin=sigmaW*pwmin
    wmax=sigmaW*pwmax
   
    if (dosingleplume) then
       ! following Soares (2004) and Witek (2011)

       UPA(1,1) = 0.1
       UPW(1,1) = 0.0
       UPM(1,1) = 0.0
       UPU(1,1)=u(1)
       UPV(1,1)=v(1)
       ! beta=0.3
       ! note that tke here is fully explicit only (at step n) if dosequential=.false., otherwise
       ! the tendency from buoyancy production, shear production, and dissipation have been added.
       ! This, if dosequential=.false., tke could be 0 and we simply add 0.01  to avoid division by zero
       UPQT(1,1)=qt(1)+beta*wqt/(sqrt(0.2)*wstar) !(sqrt(tke(1)) + 0.001 )
       UPTHV(1,1)=thetav(1)+beta*wthv/ (sqrt(0.2)*wstar) !(sqrt(tke(1)) + 0.001 )
       UPTABS(1,1)=UPTHV(1,1)/(1.+epsv*UPQT(1,1)) * (pres(1)/p00)**(rgas/cp) 
       UPQCL(1,1)=qcl(1)
       UPQCI(1,1)=qci(1)
       UPT(1,1)= (p00/pres(1))**(rgas/cp) *&
       (UPTABS(1,1) - fac_cond*(qcl(1)) - fac_sub*(qci(1)))    
       UPCF(1,1) = 0.0
       UPTHD(1,1) = UPTHV(1,1)/(1.+epsv*(UPQT(1,1)-UPQCL(1,1)-UPQCI(1,1))) - thetav(1)/(1.+epsv*qv(1)-qn(1))
       frac_mf(1) = UPA(1,1)
       qcsgs_mf(1) = qcsgs_mf(1) + UPA(1,i)*UPQCL(1,i)
       qisgs_mf(1) = qisgs_mf(1) + UPA(1,i)*UPQCI(1,i)
       qtsgs_mf(1) = qtsgs_mf(1) + UPA(1,i)*UPQT (1,i)
       tabs_mf(1)  = tabs_mf(1)  + UPA(1,I)*UPTABS(1,i)
       tke_mf(1)   = tke_mf(1)   + UPA(1,I)* 3./4. * UPW(1,i)**2

    else
      ! following Cheinet

      DO i=1,nup

         wlv=wmin+(wmax-wmin)/nup*(i-1)
         wtv=wmin+(wmax-wmin)/nup*i

         UPA(1,I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))
         !UPW(1,I)=0.5*(wlv+wtv)
         UPW(1,I)=  sigmaW/(UPA(1,I)*sqrt(2.*pi)) * (exp(-(wlv**2)/(2.*sigmaW**2)) - exp(-(wtv**2)/(2.*sigmaW**2))  ) 
 
         UPM(1,I) = UPA(1,I) * UPW(1,I)

         UPU(1,I)=u(1)
         UPV(1,I)=v(1)

         ! specific humidity needed (will be convert back in the end)
         UPQT(1,I)=qt(1)+0.32*UPW(1,I)*sigmaQT/sigmaW
         ! according to cheinet the 0.58 is for thetav, hence thetav is initialized (instead of theta)
         UPTHV(1,I)=thetav(1)+0.58*UPW(1,I)*sigmaTHV/sigmaW
         UPTABS(1,I)=UPTHV(1,I)/(1.+epsv*UPQT(1,I)) * (pres(1)/p00)**(rgas/cp) 
         UPQCL(1,I)=qcl(1)
         UPQCI(1,I)=qci(1)
         UPT(1,I)= (p00/pres(1))**(rgas/cp) *&
         (UPTABS(1,I) - fac_cond*(qcl(1)) - fac_sub*(qci(1)))    
         UPCF(1,I) = 0.0
         UPTHD(1,I) = UPTHV(1,I)/(1.+epsv*(UPQT(1,I)-UPQCL(1,I)-UPQCI(1,I))) - thetav(1)/(1.+epsv*qv(1)-qn(1))
         frac_mf(1) = frac_mf(1)+UPA(1,I)
         qcsgs_mf(1) = qcsgs_mf(1) + UPA(1,i)*UPQCL(1,i)
         qisgs_mf(1) = qisgs_mf(1) + UPA(1,i)*UPQCI(1,i)
         qtsgs_mf(1) = qtsgs_mf(1) + UPA(1,i)*UPQT (1,i)
         tabs_mf(1)  = tabs_mf(1)  + UPA(1,I)*UPTABS(1,i)
         tke_mf(1)   = tke_mf(1)   + UPA(1,I)* 3./4. * UPW(1,i)**2
      ENDDO

    end if

    qcsgs_mf(1)=qcsgs_mf(1)/frac_mf(1)
    qisgs_mf(1)=qisgs_mf(1)/frac_mf(1)
    qtsgs_mf(1)=qtsgs_mf(1)/frac_mf(1)
    tabs_mf(1)=tabs_mf(1)/frac_mf(1)

    

        !    write(*,*) ' z=',zi(1),' w=',upw(1,1),' thv=',upthv(1,1),' thl=',upthl(1,1),' qt=',upqt(1,1)*1000.,' qc=',upqc(1,1) *1000.
        !    write(*,*) 'qv=',qv(1)/(1.+qv(1))*1000. 
! updraft integration     
    DO i=1,nup
       DO k=2,nzm
 
          if (fixedeps) then
            ENT(k-1,i) = eps0
          elseif (neggerseps) then
            !ENT(k-1,i) = 2.*nuneggers*wstar/(UPW(k-1,i)+1.0e-6)/pblh 
            !ENT(k-1,i) = 2.*tauneggers*wstar/(UPW(k-1,i)+1.0e-6)/pblh 
            ENT(k-1,i) = 1./(tauneggers*UPW(k-1,i)+1.0e-6) 
          elseif (witekeps) then
            ENT(k-1,i) = 0.7/(smix(k-1)+1.0e-06)
          elseif (randomeps) then
            ! not implemented yet
          end if

          EntExp=exp(-ENT(k-1,i)*(zi(k)-zi(k-1)))

          QTn=qt(k-1)*(1.-EntExp)+UPQT(k-1,i)*EntExp
          Tn=(p00/pres(k-1))**(rgas/cp) * (t(k-1)/cp-ggr/cp*z(k-1)+fac_cond*qpl(k-1)+fac_sub*qpi(k-1))*(1.-EntExp)+UPT(k-1,i)*EntExp
          !Un=u(k)*(1.-EntExp)+UPU(k-1,i)*EntExp
          !Vn=v(k)*(1.-EntExp)+UPV(k-1,i)*EntExp

          ! compute mass flux
          ! linear increase in detrainment above cloud base
          if (k.gt.2.and.(UPQCL(k-1,i)+UPQCI(k-1,i).gt.0.0 .or. DET(k-2,i).gt. del0)) then
             DET(k-1,i) = DET(k-2,i) + deldz * (z(k-1) - z(k-2))
          ! constant detrainment below cloud base
          else
             DET(k-1,i) = del0
          end if
          EntExp=exp((ENT(k-1,i)-DET(k-1,i))*(zi(k)-zi(k-1)))
          if (.not.fixedfa) then
            UPM(k,i) = UPM(k-1,i) * EntExp
          end if


          ! all-or-nothing condensation scheme
            call condensation_edmf(QTn,Tn,presi(k),zi(k),THVn,QCLn,QCIn)
            if (donoplumesat) then
              QCLn=0.0
              QCIn=0.0
              THVn = Tn * (1.+epsv*QTn)
            end if

!          end if
       
          ! based on density potential temperature (without qp effect since homogeneous across cell)
          BUOY(k-1,i)=ggr*(0.5*(THVn+UPTHV(k-1,i))/thetar(k-1)-1.)

          !EntW=exp(-2.*(Wb+Wc*ENT(k-1,i))*(zi(k)-zi(k-1)))
          !Wn2=UPW(k-1,i)**2*EntW + (1.-EntW)*Wa*BUOY(k-1,i)/(Wb+Wc*ENT(k-1,i))

          ! new (old standard) approach that allows entrainment to force w to zero even if B>0:
          Wn2=UPW(k-1,i)**2 + 2. * (- (Wb+Wc*ENT(k-1,i)) * UPW(k-1,i)**2 &
              +  Wa * BUOY(k-1,i) ) * (zi(k)-zi(k-1))

 
          IF (Wn2 >0.) THEN
             UPW(k,i)=sqrt(Wn2) 
             if (.not.fixedfa) then 
                UPA(k,i) = UPM(k,i) / UPW(k,i) 
             else
                UPA(k,i) = UPA(k-1,i)
                UPM(k,i) = UPW(k,i) * UPA(k,i)  
             end if
             IF (fixedfa.or.(UPA(k,i).ge.facrit*UPA(1,i))) THEN
               UPTHV(k,i)=THVn
               UPTABS  (k,i)=THVn/(1.+epsv*(QTn-QCLn-QCIn)-QCLn-QCIn)*(presi(k)/p00)**(rgas/cp)
               UPT(k,i)=Tn
               UPQT(k,i)=QTn
               UPQCL(k,i)=QCLn
               UPQCI(k,i)=QCIn
             !  UPU(k,i)=Un
             !  UPV(k,i)=Vn
             ELSE
                UPW(k,i)= 0.
                UPA(k,i)= 0.
                EXIT
             END IF
             UPTHD(k-1,i)=0.5*(UPTHV(k,i)/(1.+epsv*(UPQT(k,i)-UPQCL(k,i)-UPQCI(k,i))-UPQCL(k,i)-UPQCI(k,i))+&
     UPTHV(k-1,i)/(1.+epsv*(UPQT(k-1,i)-UPQCL(k-1,i)-UPQCI(k,i))-UPQCL(k-1,i)-UPQCI(k-1,i)))-theta(k-1)
          ELSE
            EXIT
          END IF 
       ENDDO
    ENDDO

! computing variables needed for tendency calculation

! mass flux is zero at surface and top
    sumM(1)       = 0.0
    sumMt(1)      = 0.0
    sumMrt(1)     = 0.0
    sumMrp(1)     = 0.0
    sumMu(1)      = 0.0
    sumMv(1)      = 0.0
    sumMtke(1)    = 0.0
    sumMthv(1)    = 0.0
    sumDEF2(1)     = 0.0
    sumM(nz)       = 0.0
    sumMt(nz)      = 0.0
    sumMrt(nz)     = 0.0
    sumMrp(nz)     = 0.0
    sumMu(nz)      = 0.0
    sumMv(nz)      = 0.0
    sumMtke(nz)    = 0.0
    sumMthv(nz)    = 0.0
    sumDEF2(nz)     = 0.

    sumMu      = 0.0
    sumMv      = 0.0
    sumMtke    = 0.0
    sumMrp     = 0.0

    DO k=2,nzm
      DO i=1,nup
        sumM(k)      =sumM(k)      +UPA(K,I)*UPW(K,I)
        if (dotlflux) then
          ! flux in terms of liquid/ice pot temperature (needs to include precip contribution since its considered in env)
          sumMt(k)     =sumMt(k)     +UPA(K,i)*(UPT(K,I) - (p00/presi(k))**(rgas/cp)*&
           (fac_cond*0.5*(qpl(k)+qpl(k-1)) +fac_sub*0.5*(qpi(k)+qpi(k-1))  ))*UPW(K,I)
        else
          ! flux in terms of liquid/frozen water static energy (convert to h)
          sumMt(k)     =sumMt(k)     +UPA(K,i)*(cp*(presi(k)/p00)**(rgas/cp)*UPT(K,I)-&
          lcond*0.5*(qpl(k)+qpl(k-1))-lsub*0.5*(qpi(k)+qpi(k-1))+ggr*zi(k))*UPW(K,I)
        end if
        sumMrt(k)    =sumMrt(k)    +UPA(K,i)*UPQT(K,I)*UPW(K,I) 
        sumMthv(k)   =sumMthv(k)   +UPA(K,i)*UPW(K,I) * (UPTHV(k,i) - 0.5 *(thetar(k-1)+thetar(k)))
        sumDEF2(k)    = sumDEF2(k) + 3. * UPA(K,i) * ((UPW(K+1,I)-UPW(K-1,I))/(zi(k+1)-zi(k-1)))**2
        !sumMu(k)  =sumMu(k)+UPA(K,i)*UPW(K,I)*UPU(K,I)
        !sumMv(k)  =sumMv(k)+UPA(K,i)*UPW(K,I)*UPV(K,I)
        qcsgs_mf(k) = qcsgs_mf(k) + UPA(K,i)*UPQCL(k,i)
        qisgs_mf(k) = qisgs_mf(k) + UPA(K,i)*UPQCI(k,i)
        qtsgs_mf(k) = qtsgs_mf(k) + UPA(K,i)*UPQT (k,i)
        tabs_mf(k)  = tabs_mf(k)  + UPA(K,I)*UPTABS(k,i)
        tke_mf(k)   = tke_mf(k)   + UPA(K,I)* 3./4. * UPW(k,i)**2
        frac_mf(k)  = frac_mf(k)  + UPA(K,i)
        if (UPQCL(k,i)+UPQCI(k,i).gt.0.0) cfrac_mf(k) = cfrac_mf(k) + UPA(K,i)
      ENDDO
      if (frac_mf(k).gt.0.) then
         qcsgs_mf(k) = qcsgs_mf(k) / frac_mf(k)
         qisgs_mf(k) = qisgs_mf(k) / frac_mf(k)
         qtsgs_mf(k) = qtsgs_mf(k) / frac_mf(k)
         tabs_mf(k)  = tabs_mf(k)  / frac_mf(k)
      end if
    ENDDO
    

end subroutine edmf

subroutine condensation_edmf(QT,THLI,P,zlev,THV,QC,QI)
!
! zero or one condensation for edmf: calculates THV and QC, and QI
!
use params
use micro_params
use vars, only : qsatw, qsati

implicit none
real,intent(in) :: QT,THLI,P, zlev
real,intent(out):: THV,QC,QI

integer :: niter,i
real :: diff,t,qs,qnold, an, bn, qn, om

an = 1./(tbgmax-tbgmin)
bn = tbgmin * an

! number of iterations
niter=100
! minimum difference
diff=1.e-7

QC=0.
QI=0.
QN=0.

!print*, '+++++++++++++++++++++++++'
!print*, '+++++++++++++qc=0++++++++'
!print*, '+++++++++++++++++++++++++'
do i=1,NITER
! get tabs
T = THLI* (P/p00)**(rgas/cp) +fac_cond*QC+fac_sub*QI
! WL get saturation mixing ratio
if (T.ge.tbgmax) then
  QS=qsatw(T,P)
  om=1.
elseif (T.le.tbgmin) then
  QS=qsati(T,P)
  om=0.
else
  om = an*T-bn
  QS=om*qsatw(T,P)+(1.-om)*qsati(T,P)
end if
QNOLD=QN
QN=max(0.5*QN+0.5*(QT-QS),0.)
QC= om * QN
QI= (1.-om) * QN
!write(*,'(5A)') ' tabs','  ','   qs','  ','   qc'
!write(*,'(F5.1,2x,F5.2,2x,F5.2)'),T,qs*1000.,qc*1000.
if (abs(QN-QNOLD)<Diff) exit
!print*, '+++++++++++++++++++++++++'
enddo

T = THLI* (P/p00)**(rgas/cp) +fac_cond*QC+fac_sub*QI
if (T.ge.tbgmax) then
  QS=qsatw(T,P)
elseif (T.le.tbgmin) then
  QS=qsati(T,P)
else
  om = an*T-bn
  QS=om*qsatw(T,P)+(1.-om)*qsati(T,P)
end if
QN=max(QT-QS,0.)
QC= om * QN
QI= (1.-om) * QN

THV = (THLI +(fac_cond*QC+fac_sub*QI)*(p00/P)**(rgas/cp)) * (1.+epsv*(QT-QN)-QN)

end subroutine condensation_edmf

