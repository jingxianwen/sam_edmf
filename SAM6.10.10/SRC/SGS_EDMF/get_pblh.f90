subroutine get_pblh(i,j,thetav,thvflux,pblh)

use vars
use params
use grid
use sgs, only : pblhfluxmin, pblhthvgrad, tke

implicit none

!input
integer, intent(in) :: i,j
real, dimension(nzm),intent(in) :: thetav, thvflux
!output
real :: pblh


    !---------------------------------------------------------------
    !             NOTES ON THE PBLH FORMULATION
    !
    !The 1.5-theta-increase method defines PBL heights as the level at 
    !which the potential temperature first exceeds the minimum potential 
    !temperature within the boundary layer by 1.5 K. When applied to 
    !observed temperatures, this method has been shown to produce PBL-
    !height estimates that are unbiased relative to profiler-based 
    !estimates (Nielsen-Gammon et al. 2008). However, their study did not
    !include LLJs. Banta and Pichugina (2008) show that a TKE-based 
    !threshold is a good estimate of the PBL height in LLJs. Therefore,
    !a hybrid definition is implemented that uses both methods, weighting
    !the TKE-method more during stable conditions (PBLH < 400 m).
    !A variable tke threshold (TKEeps) is used since no hard-wired
    !value could be found to work best in all conditions.
    !---------------------------------------------------------------
! local variables
    REAL ::  PBLH_TKE,qtke,qtkem1,wt,maxqke,TKEeps,minthv, X1, X2, X3, Y1, Y2, Y3, aa, bb, dzavg
    REAL :: delt_thv   !delta theta-v; dependent on land/sea point
    REAL,dimension(nzm) :: dthvdz
    REAL, PARAMETER :: sbl_lim  = 200. !typical scale of stable BL (m).
    REAL, PARAMETER :: sbl_damp = 400. !transition length for blending (m).
    INTEGER :: K,kthv,ktke

    IF (pblhfluxmin) then ! compute cpbl top as level with flux minimum


    pblh = zi(2)
    !FIND MIN Flux
    k = 1
    kthv = 1
    minthv = 9.E9
    DO WHILE (k.le.nzm)
       ! use flux from previous step
       qtke  = thvflux(k) ! 
       IF (minthv > qtke) then
           minthv = qtke
           kthv = k
       ENDIF
       k = k+1
    ENDDO 
    pblh=z(kthv)

    ELSEIF (pblhthvgrad) THEN

    k = 1
    kthv = 1
    maxqke = -1000000000000000.
    dthvdz=0.0
    DO WHILE (k .LE. nzm-2)
       if (k.eq.1) then
         dthvdz(k) = (thetav(k+1) - thetav(k))/(z(k+1)-z(k))
       elseif (k.eq.nzm-1.or.k.eq.2) then
         dthvdz(k) = (thetav(k+1) - thetav(k-1))/(z(k+1)-z(k-1))
       else 
         dzavg = 0.5*(0.25*(z(k+2)-z(k-2)) + 0.5*(z(k+1)-z(k-1))) 
         dthvdz(k) = (-thetav(k+2) +8.* thetav(k+1) - 8.* thetav(k-1) + thetav(k-2))/(12.*dzavg)
       end if
       IF (maxqke < dthvdz(k)) then
           maxqke = dthvdz(k)
           kthv = k
       ENDIF
       k = k+1
    ENDDO
    if (kthv.eq.nzm-2) kthv=nzm-3
 

    !compute quadratic fit to thetav/dz using three points
    !then let height of vertex be the PBL height to allow for PBL height to fall in between levels
    if (kthv .eq. 1) then
       pblh = z(1)
    else

       Y1 = dthvdz(kthv-1)
       Y2 = dthvdz(kthv)
       Y3 = dthvdz(kthv+1)
       X1=z(kthv-1)
       X2=z(kthv)
       X3=z(kthv+1)
       aa=((Y2-Y1)*(X1-X3) + (Y3-Y1)*(X2-X1))/((X1-X3)*(X2**2-X1**2) + (X2-X1)*(X3**2-X1**2))
       bb=((Y2 - Y1) - aa*(X2**2 - X1**2)) / (X2-X1)


      if (aa.eq.0.) then
        ! linear profile of gradient with bb.ne.0.0
        if (bb.gt.0.0) then
          pblh = z(kthv+1)
        elseif(bb.lt.0.0) then
          pblh = z(kthv-1)
        else
          !constant gradient
          pblh = z(kthv)
        end if
      else
        pblh =  - bb/(2.*aa)
      end if

    end if


    ELSE 


    !FIND MAX TKE AND MIN THETAV IN THE LOWEST 500 M
    k = 1
    kthv = 1
    ktke = 1
    maxqke = 0.
    minthv = 9.E9
    DO WHILE (zi(k+1) .LE. 500.)
       qtke  =MAX(tke(i,j,k),0.)   ! maximum QKE
       IF (maxqke < qtke) then
           maxqke = qtke
           ktke = k
       ENDIF
       IF (minthv > thetav(k)) then
           minthv = thetav(k)
           kthv = k
       ENDIF
       k = k+1
    ENDDO
    !TKEeps = maxtke/20. = maxqke/40.
    TKEeps = maxqke/40. 
    TKEeps = MAX(TKEeps,0.025)

    !FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
    IF(ocean) THEN                                            
        ! WATER
        delt_thv = 0.75
    ELSE         
        ! LAND     
        delt_thv = 1.5  
    ENDIF

    pblh=0.
    k = kthv+1
    DO WHILE (pblh .EQ. 0.) 
       IF (thetav(k) .GE. (minthv + delt_thv))THEN
          pblh = z(k) - adzw(k)*dz* &
             & MIN(1.,(thetav(k)-(minthv + delt_thv))/ &
             & MAX(thetav(k)-thetav(k-1),1.E-6))
       ENDIF
       k = k+1
       IF (k .EQ. nzm-1) pblh = zi(2) !EXIT SAFEGUARD
    ENDDO
    !print*,"IN GET_PBLH:",thsfc,zi

    !FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
    !THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
    !THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE 
    !WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.

    PBLH_TKE=0.
    k = 2
    DO WHILE (PBLH_TKE .EQ. 0.) 
       !QKE CAN BE NEGATIVE (IF CKmod == 0)... MAKE TKE NON-NEGATIVE.
       qtke  =MAX(tke(i,j,k),0.)      ! maximum TKE
       qtkem1=MAX(tke(i,j,k-1),0.)
       IF (qtke .LE. TKEeps) THEN
           PBLH_TKE = z(k) - adzw(k)*dz* &
             & MIN((TKEeps-qtke)/MAX(qtkem1-qtke, 1.E-6), 1.0)
           !IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
           PBLH_TKE = MAX(PBLH_TKE,zi(2))
           !print *,"PBLH_TKE:",i,j,PBLH_TKE, Qke1D(k)/2., zw1D(kts+1)
       ENDIF
       k = k+1
       IF (k .EQ. nzm-1) PBLH_TKE = zi(2) !EXIT SAFEGUARD
    ENDDO

    !With TKE advection turned on, the TKE-based PBLH can be very large 
    !in grid points with convective precipitation (> 8 km!),
    !so an artificial limit is imposed to not let PBLH_TKE exceed 4km.
    !This has no impact on 98-99% of the domain, but is the simplest patch
    !that adequately addresses these extremely large PBLHs.
    PBLH_TKE = MIN(PBLH_TKE,4000.)

    !BLEND THE TWO PBLH TYPES HERE: 
    wt=.5*TANH((pblh - sbl_lim)/sbl_damp) + .5
    pblh=PBLH_TKE*(1.-wt) + pblh*wt

    END IF  


end subroutine get_pblh
