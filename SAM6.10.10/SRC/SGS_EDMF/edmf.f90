! WL this routine provides the fluxes sumM(nz), sumMrt(nz), sumMu(nz), sumMv(nz), sumMt(nz), and sumMtke(nz)
! WL that result from the multiplume model provided by Kay Suselj and adapted by WL
subroutine edmf

use vars
use grid
use params
use microphysics, only : q, qn, qp
use sgs
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

      REAL, DIMENSION(1:nz,1:nup) :: UPM, UPW, UPT, UPTABS, UPTHV, UPQT, UPQCL,&
                                     UPQCI,UPA, UPU, UPV  
      REAL, DIMENSION(1:nz,1:nup) :: ENT, UPCF, BUOY

      INTEGER :: K,N,i,j, ic, jc
      REAL :: wthv,wqt,wthl,qstar,thstar,sigmaW,sigmaQT,sigmaTH,sigmaTHV,zs, &
           pwmax,wmin,wmax,wlv,wtv,thetav1,theta
      REAL :: QTn,Tn,THVn,QCLn,QCIn,Un,Vn,Wn2,EntEXP,EntW, hlp, acrit, Wa, thetavenv

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
 BUOY=0.
 qcsgs_mf=0.
 qisgs_mf=0.
 cfrac_mf=0.
 frac_mf=0.
 cfrac_mf1D=0.
 frac_mf1D=0.

 Wa = 0.4* Wc + 0.3
! set initial conditions for updrafts
 zs=50.
 pwmax=3.5

 do i=1,nx
 do j=1,ny
 
 ! surface fluxes
 ! sensible heat flux (K m s-1)
 wthl = fluxbt(i,j)
 ! latent heat flux (m s-1)
 wqt  = fluxbq(i,j) 
 ! virtual pot. temp flux
 wthv = (1.+epsv*qv(i,j,1))*fluxbt(i,j) + epsv*tabs(i,j,1)*(1000./pres(1))**(rgas/cp)*fluxbq(i,j)

 thetav1 = (1.+epsv*qv(i,j,1)-(qn(i,j,1) +qp(i,j,1)))*tabs(i,j,1)*(1000./pres(1))**(rgas/cp)
 theta   = thetav1/(1.+epsv*qv(i,j,1)-(qn(i,j,1) +qp(i,j,1)))

 wstar(i,j)=max(0.,(ggr/thetav1*wthv*pblh(i,j))**(1./3.))

 ! quit in case of a non-positive buoyancy flux
 if (wthv.le.0.0) cycle
 
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


! see Lenschow et al. (1980), JAS
    qstar=wqt/wstar(i,j)
    thstar=wthl/wstar(i,j) 
    sigmaW=1.34*wstar(i,j)*(zs/pblh(i,j))**(1./3.)*(1.-0.8*zs/pblh(i,j))
    sigmaQT=1.34*qstar*(zs/pblh(i,j))**(-1./3.)
    sigmaTH=1.34*thstar*(zs/pblh(i,j))**(-1./3.)
 
! now get sigmaTHV from linearization of pot. temp. and using a corr. coeff between qv and theta of 0.75 (see Sobjan 1991)
    sigmaTHV=sqrt(sigmaTH**2+epsv**2*theta**2*sigmaQT**2 &
        + 2.*epsv*theta*0.75*sigmaTH*sigmaQT)

    wmin=sigmaW*pwmin
    wmax=sigmaW*pwmax
   
    if (dosingleplume) then
       ! following Soares (2004) and Witek (2011)

       UPA(1,1) = 0.1
       UPW(1,1) = 0.0
       UPM(1,1) = 0.0
       ic=i+1
       UPU(1,1)=0.5*(u(i,j,1)+u(ic,j,1))
       if (RUN3D) then
         jc=j+YES3D
         UPV(1,1)=0.5*(v(i,j,1)+v(i,jc,1))
       else
         UPV(1,1)=v(i,j,1)
       end if
       ! beta=0.3
       ! note that tke here is fully explicit only (at step n) if dosequential=.false., otherwise
       ! the tendency from buoyancy production, shear production, and dissipation have been added.
       ! This, if dosequential=.false., tke could be 0 and we simply add 0.01  to avoid division by zero
       UPQT(1,1)=q(i,j,1)+beta*wqt/(sqrt(0.2)*wstar(i,j)) !(sqrt(tke(1)) + 0.001 )
       UPTHV(1,1)=thetav1+beta*wthv/ (sqrt(0.2)*wstar(i,j)) !(sqrt(tke(1)) + 0.001 )
       UPTABS(1,1)=UPTHV(1,1)/(1.+epsv*UPQT(1,1)) * (pres(1)/1000.)**(rgas/cp) 
       UPQCL(1,1)=qcl(i,j,1) 
       UPQCI(1,1)=qci(i,j,1) 
       UPT(1,1)= (1000./pres(1))**(rgas/cp) *&
       (UPTABS(1,1) - fac_cond*(qcl(i,j,1)) - fac_sub*(qci(i,j,1)))    
       UPCF(1,1) = 0.0
       frac_mf(i,j,1) = UPA(1,1)
       qcsgs_mf(i,j,1) =  UPA(1,1)*UPQCL(1,1)
       qisgs_mf(i,j,1) =  UPA(1,1)*UPQCI(1,1)

    else
      ! following Cheinet

      DO N=1,nup

         wlv=wmin+(wmax-wmin)/nup*(N-1)
         wtv=wmin+(wmax-wmin)/nup*N

         UPA(1,N)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))
         !UPW(1,N)=0.5*(wlv+wtv)
         UPW(1,N)=  sigmaW/(UPA(1,N)*sqrt(2.*pi)) * (exp(-(wlv**2)/(2.*sigmaW**2)) - exp(-(wtv**2)/(2.*sigmaW**2))  ) 
 
         UPM(1,N) = UPA(1,N) * UPW(1,N)

         ic=i+1
         UPU(1,N)=0.5*(u(i,j,1)+u(ic,j,1))
         if (RUN3D) then
           jc=j+YES3D
           UPV(1,N)=0.5*(v(i,j,1)+v(i,jc,1))
         else
           UPV(1,N)=v(i,j,1)
         end if

         ! specific humidity needed (will be convert back in the end)
         UPQT(1,N)=q(i,j,1)+0.32*UPW(1,N)*sigmaQT/sigmaW
         ! according to cheinet the 0.58 is for thetav, hence thetav is initialized (instead of theta)
         UPTHV(1,N)=thetav1+0.58*UPW(1,N)*sigmaTHV/sigmaW
         UPTABS(1,N)=UPTHV(1,N)/(1.+epsv*UPQT(1,N)) * (pres(1)/1000.)**(rgas/cp) 
         UPQCL(1,N)=qcl(i,j,1) 
         UPQCI(1,N)=qci(i,j,1)
         UPT(1,N)= (1000./pres(1))**(rgas/cp) *&
         (UPTABS(1,N) - fac_cond*(qcl(i,j,1)) - fac_sub*(qci(i,j,1)))    
         UPCF(1,N) = 0.0
         frac_mf(i,j,1) = frac_mf(i,j,1)+UPA(1,N)
         qcsgs_mf(i,j,1) = qcsgs_mf(i,j,1) + UPA(1,N)*UPQCL(1,N)
         qisgs_mf(i,j,1) = qisgs_mf(i,j,1) + UPA(1,N)*UPQCI(1,N)
      ENDDO

    end if

    qcsgs_mf(i,j,1)=qcsgs_mf(i,j,1)/frac_mf(i,j,1)
    qisgs_mf(i,j,1)=qisgs_mf(i,j,1)/frac_mf(i,j,1)
    frac_mf1D(1)  = frac_mf1D(1) + frac_mf(i,j,1)
    cfrac_mf1D(1) = cfrac_mf1D(1) + cfrac_mf(i,j,1)

    

        !    write(*,*) ' z=',zi(1),' w=',upw(1,1),' thv=',upthv(1,1),' thl=',upthl(1,1),' qt=',upqt(1,1)*1000.,' qc=',upqc(1,1) *1000.
        !    write(*,*) 'qv=',qv(1)/(1.+qv(1))*1000. 
! updraft integration     
    DO N=1,nup
       DO k=2,nzm
 
          if (fixedeps) then
            ENT(k-1,N) = eps0
          else
            ENT(k-1,N) = 1./(tauneggers*UPW(k-1,N)+1.0e-6) 
          end if

          EntExp=exp(-ENT(k-1,N)*(zi(k)-zi(k-1)))

          QTn=q(i,j,k-1)*(1.-EntExp)+UPQT(k-1,N)*EntExp
          Tn=(1000./pres(k-1))**(rgas/cp) * (t(i,j,k-1)-ggr/cp*z(k-1)+fac_cond*qpl(i,j,k-1) + fac_sub*qpi(i,j,k-1))*(1.-EntExp)+UPT(k-1,N)*EntExp
          !Un=0.5*(u(i,j,k-1)+u(i+1,j,k-1))*(1.-EntExp)+UPU(k-1,i)*EntExp
          !if (RUN3D) then
          !Vn=0.5*(v(i,j,k-1)+v(i,j+YES3D,k-1))*(1.-EntExp)+UPV(k-1,i)*EntExp
          !else
          !Vn=v(i,j,k-1)*(1.-EntExp)+UPV(k-1,i)*EntExp
          !end if
          


          ! all-or-nothing condensation scheme
            call condensation_edmf(QTn,Tn,presi(k),zi(k),THVn,QCLn,QCIn)
            if (donoplumesat) then
              QCLn=0.0
              THVn = Tn * (1.+epsv*QTn)
            end if

!          end if
       
          ! based on density potential temperature (without qp effect since homogeneous across cell)
          thetavenv = (1.+epsv*qv(i,j,k-1)-qn(i,j,k-1))*tabs(i,j,k-1)*(1000./pres(k-1))**(rgas/cp)
          BUOY(k-1,N)=ggr*(0.5*(THVn+UPTHV(k-1,N))/thetavenv-1.)

          !EntW=exp(-2.*(Wb+Wc*ENT(k-1,i))*(zi(k)-zi(k-1)))
          !Wn2=UPW(k-1,i)**2*EntW + (1.-EntW)*Wa*BUOY(k-1,i)/(Wb+Wc*ENT(k-1,i))

          ! new (old standard) approach that allows entrainment to force w to zero even if B>0:
          Wn2=UPW(k-1,N)**2 + 2. * (- (Wb+Wc*ENT(k-1,N)) * UPW(k-1,N)**2 &
              +  Wa * BUOY(k-1,N) ) * (zi(k)-zi(k-1))

 
          IF (Wn2 .gt.0.) THEN
             UPW(k,N)=sqrt(Wn2) 
             UPA(k,N) = UPA(k-1,N)
             UPM(k,N) = UPW(k,N) * UPA(k,N)  
             UPTHV(k,N)=THVn
             UPTABS  (k,N)=THVn/(1.+epsv*(QTn-QCLn-QCIn)-QCLn-QCIn)*(presi(k)/1000.)**(rgas/cp)
             UPT(k,N)=Tn
             UPQT(k,N)=QTn
             UPQCL(k,N)=QCLn
             UPQCI(k,N)=QCIn
           !  UPU(k,i)=Un
           !  UPV(k,i)=Vn
          ELSE
            EXIT
          END IF 
       ENDDO
    ENDDO

! computing variables needed for tendency calculation

! mass flux is zero at surface and top

    DO k=2,nzm
      DO N=1,nup
        sgs_field_sumM(i,j,k,1)=sgs_field_sumM(i,j,k,1) + UPA(K,N)*UPW(K,N)
        sgs_field_sumM(i,j,k,5)=sgs_field_sumM(i,j,k,5) + UPA(K,N)*((presi(k)/1000.)**(rgas/cp)*UPT(K,N)+ggr/cp*zi(k)&
            -fac_cond*qpl(i,j,k)-fac_sub*qpi(i,j,k))*UPW(K,N)
        sgs_field_sumM(i,j,k,6)=sgs_field_sumM(i,j,k,6) + UPA(K,N)*UPQT(K,N)*UPW(K,N)
        !sumMu(k)  =sumMu(k)+UPA(K,i)*UPW(K,I)*UPU(K,I)
        !sumMv(k)  =sumMv(k)+UPA(K,i)*UPW(K,I)*UPV(K,I)
        qcsgs_mf(i,j,k) = qcsgs_mf(i,j,k) + UPA(K,N)*UPQCL(k,N)
        qisgs_mf(i,j,k) = qisgs_mf(i,j,k) + UPA(K,N)*UPQCI(k,N)
        frac_mf(i,j,k)  = frac_mf(i,j,k)  + UPA(K,N)
        if (UPQCL(k,N)+UPQCI(k,N).gt.0.0) cfrac_mf(i,j,k) = cfrac_mf(i,j,k) + UPA(K,N)
      ENDDO
      if (frac_mf(i,j,k).gt.0.) then
         qcsgs_mf(i,j,k) = qcsgs_mf(i,j,k) / frac_mf(i,j,k)
         qisgs_mf(i,j,k) = qisgs_mf(i,j,k) / frac_mf(i,j,k)
      end if
      frac_mf1D(k)  = frac_mf1D(k) + frac_mf(i,j,k)
      cfrac_mf1D(k) = cfrac_mf1D(k) + cfrac_mf(i,j,k)
    ENDDO


  end do
  end do
    

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
T = THLI* (P/1000.)**(rgas/cp) +fac_cond*QC+fac_sub*QI
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

T = THLI* (P/1000.)**(rgas/cp) +fac_cond*QC+fac_sub*QI
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

THV = (THLI +(fac_cond*QC+fac_sub*QI)*(1000./P)**(rgas/cp)) * (1.+epsv*(QT-QN)-QN)

end subroutine condensation_edmf

