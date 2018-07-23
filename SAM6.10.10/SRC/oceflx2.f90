
subroutine oceflx2(pmidm1, ubot, vbot, tbot, qbot, thbot, zbot, ts, &
                  shf, lhf, taux, tauy, ssq)

use params, only: salt_factor
implicit none	

!
! Compute ocean to atmosphere surface fluxes of sensible, latent heat
! and stress components:
!
! Assume:
!   1) Neutral 10m drag coeff: 
!         cdn = .0027/U10N + .000142 + .0000764 U10N
!   2) Neutral 10m stanton number: 
!         ctn = .0327 sqrt(cdn), unstable
!         ctn = .0180 sqrt(cdn), stable
!   3) Neutral 10m dalton number:  
!         cen = .0346 sqrt(cdn)
!
! Note:
!   1) here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   2) wind speeds should all be above a minimum speed (umin)
!
!--------------------------Code History---------------------------------
!
! Original:      Bill Large/M.Vertenstein, Sep. 1995 for CCM3.5
! Standardized:  L. Buja,     Feb 1996
! Reviewed:      B. Briegleb, March 1996
! Adopted for LES by Marat Khairoutdinov, July 1998
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
!                         
      real pmidm1  ! Bottom level pressure, mb
      real ubot    ! Bottom level u wind
      real vbot    ! Bottom level v wind
      real tbot    ! Bottom level temperature
      real qbot    ! Bottom level specific humidity
      real thbot   ! Bottom level potential temperature ---> SEE NOTE BELOW.
      real zbot    ! Bottom level height above surface
      real ts      ! Surface temperature
!
!bloss(2015-01-09): !!NOTE!! that liquid-ice static energy at the surface
!   is input into this routine in SAM rather than potential temperature.
!   As written, the potential temperature would only work if the surface
!   pressure were used as the reference pressure in place of the usual 
!   1000 hPa.
!                         
! Output arguments        
!                         
      real shf     ! Initial sensible heat flux (Km/s)
      real lhf     ! Initial latent heat flux (m/s)
      real taux    ! X surface stress (N/m2)
      real tauy    ! Y surface stress (N/m2)
      real ssq     ! Surface saturation specific humidity
!
!---------------------------Local variables-----------------------------
!
      real tref    	  ! 2m reference temperature 
      real ustar          ! ustar
      real tstar          ! tstar
      real qstar          ! qstar
      real u10n           ! neutral 10 m wind speed over ocean
      real vmag           ! Surface wind magnitude
      real thvbot         ! Bottom lev virtual potential temp
      real delt           ! potential T difference (K)
      real delq           ! specific humidity difference (kg/kg)
      real rdn            ! sqrt of neutral exchange coeff (momentum)
      real rhn            ! sqrt of neutral exchange coeff (heat)
      real ren            ! sqrt of neutral exchange coeff (tracers)          
      real rd             ! sqrt of exchange coeff (momentum)
      real rh             ! sqrt of exchange coeff (heat)
      real re             ! sqrt of exchange coeff (tracers)
      real hol            ! Ref hgt (10m) / monin-obukhov length
      real xsq            ! Temporary variable
      real xqq            ! Temporary variable
      real alz            ! ln(zbot/z10)
      real cp             ! Specific heat of moist air
      real tau     ! Reference height stress
      real psimh          ! Stability funct at ref lev (momentum)
      real psixh          ! Stability funct at ref lev (heat & tracers) 
      real stable         ! Stability factor
      real rbot    	  ! Density at bottom model level
      real bn             ! exchange coef funct for interpolation
      real bh             ! exchange coef funct for interpolation
      real fac            ! interpolation factor
      real ln0            ! log factor for interpolation
      real ln3            ! log factor for interpolation
      real ztref          ! reference height for air temperature
      real gravit	  ! =9.81 m/s2
      real ltheat  	  ! Latent heat for given srf conditions
      
      real xkar     ! Von Karman constant
      real zref     ! 10m reference height
      real umin     ! Minimum wind speed at bottom level

      parameter  (xkar = 0.4 , &
            zref =   10.0      , &
            umin =    1.0, &
            ztref=2.0 )


      real hlp



!---------------------------------------------------------------
! Set up necessary variables
!---------------------------------------------------------------
         cp     = 1005.
	 gravit = 9.81
	 ltheat = 2.51e6

         rbot = pmidm1*100.  / (287.*tbot )
         delt   = thbot-gravit/cp*zbot  - ts


         taux  = 0.0 
         tauy  = 0.0 

        
         hlp = 50./cp/rbot/2.

         shf  = hlp * 2. - hlp * (delt+2.)
         lhf  = 0.0

end subroutine oceflx2

 
