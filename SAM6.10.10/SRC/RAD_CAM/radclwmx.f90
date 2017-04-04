!bloss #include <misc.h>
!bloss #include <params.h>
subroutine radclwmx(lchnk   ,ncol    ,                            &
                    lwupcgs ,tnm     ,qnm     ,o3vmr   , &
                    pmid    ,pint    ,pmln    ,piln    ,          &
                             n2o     ,ch4     ,cfc11   ,cfc12   , &
                    cld     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                    flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                    flut    ,flutc   , &
                    aer_mass,fnl     ,fsul, fsdl    ,ful     , fdl)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute longwave radiation heating rates and boundary fluxes
! 
! Method: 
! Uses broad band absorptivity/emissivity method to compute clear sky;
! assumes randomly overlapped clouds with variable cloud emissivity to
! include effects of clouds.
!
! Computes clear sky absorptivity/emissivity at lower frequency (in
! general) than the model radiation frequency; uses previously computed
! and stored values for efficiency
!
! Note: This subroutine contains vertical indexing which proceeds
!       from bottom to top rather than the top to bottom indexing
!       used in the rest of the model.
! 
! Author: B. Collins
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r4 => shr_kind_r4
   use ppgrid
   use radae, only: nbands, radems, radabs, radtpl, abstot_3d, absnxt_3d, emstot_3d
   use volcrad
   use quicksort, only: quick_sort

   implicit none

   integer pverp2,pverp3,pverp4
   parameter (pverp2=pver+2,pverp3=pver+3,pverp4=pver+4)

   real(r4) cldmin
   parameter (cldmin = 1.0e-30)
!------------------------------Commons----------------------------------
!bloss #include <comctl.h>
!-----------------------------------------------------------------------
!bloss #include <crdcon.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   integer, intent(in) :: nmxrgn(pcols)         ! Number of maximally overlapped regions
   real(r4), intent(in) :: pmxrgn(pcols,pverp)  ! Maximum values of pmid for each
                                                !    maximally overlapped region.
                                                !    0->pmxrgn(i,1) is range of pmid for
                                                !    1st region, pmxrgn(i,1)->pmxrgn(i,2) for
                                                !    2nd region, etc

   real(r4), intent(in) :: lwupcgs(pcols)       ! Longwave up flux in CGS units
!
! Input arguments which are only passed to other routines
!
   real(r4), intent(in) :: tnm(pcols,pver)      ! Level temperature
   real(r4), intent(in) :: qnm(pcols,pver)      ! Level moisture field
   real(r4), intent(in) :: o3vmr(pcols,pver)    ! ozone volume mixing ratio
   real(r4), intent(in) :: pmid(pcols,pver)     ! Level pressure
   real(r4), intent(in) :: pint(pcols,pverp)    ! Model interface pressure
   real(r4), intent(in) :: pmln(pcols,pver)     ! Ln(pmid)
   real(r4), intent(in) :: piln(pcols,pverp)    ! Ln(pint)
   real(r4), intent(in) :: n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
   real(r4), intent(in) :: ch4(pcols,pver)      ! methane mass mixing ratio
   real(r4), intent(in) :: cfc11(pcols,pver)    ! cfc11 mass mixing ratio
   real(r4), intent(in) :: cfc12(pcols,pver)    ! cfc12 mass mixing ratio
   real(r4), intent(in) :: cld(pcols,pver)      ! Cloud cover
   real(r4), intent(in) :: emis(pcols,pver)     ! Cloud emissivity
   real(r4), intent(in) :: aer_mass(pcols,pver) ! STRAER mass in layer

!
! Output arguments
!
   real(r4), intent(out) :: qrl(pcols,pver)      ! Longwave heating rate
   real(r4), intent(out) :: flns(pcols)          ! Surface cooling flux
   real(r4), intent(out) :: flnt(pcols)          ! Net outgoing flux
   real(r4), intent(out) :: flut(pcols)          ! Upward flux at top of model
   real(r4), intent(out) :: flnsc(pcols)         ! Clear sky surface cooing
   real(r4), intent(out) :: flntc(pcols)         ! Net clear sky outgoing flux
   real(r4), intent(out) :: flutc(pcols)         ! Upward clear-sky flux at top of model
   real(r4), intent(out) :: flwds(pcols)         ! Down longwave flux at surface
   real(r4),intent(out)  :: fsul(pcols,pverp)    ! clear sky net flux at interfaces
   real(r4),intent(out)  :: fsdl(pcols,pverp)    ! clear sky net flux at interfaces
   real(r4),intent(out)  :: fnl(pcols,pverp)     ! net flux at interfaces
   real(r4),intent(out)  :: ful(pcols,pverp)     ! upward flux at interfaces
   real(r4),intent(out)  :: fdl(pcols,pverp)     ! downward flux at interfaces
!
!---------------------------Local variables-----------------------------
!
   integer i                 ! Longitude index
   integer ilon              ! Longitude index
   integer ii                ! Longitude index
   integer iimx              ! Longitude index (max overlap)
   integer k                 ! Level index
   integer k1                ! Level index
   integer k2                ! Level index
   integer k3                ! Level index
   integer km                ! Level index
   integer km1               ! Level index
   integer km3               ! Level index
   integer km4               ! Level index
   integer irgn              ! Index for max-overlap regions
   integer l                 ! Index for clouds to overlap
   integer l1                ! Index for clouds to overlap
   integer n                 ! Counter

!
   real(r4) :: plco2(pcols,pverp)   ! Path length co2
   real(r4) :: plh2o(pcols,pverp)   ! Path length h2o
   real(r4) tmp(pcols)           ! Temporary workspace
   real(r4) tmp2(pcols)          ! Temporary workspace
   real(r4) tmp3(0:pverp)        ! Temporary workspace
   real(r4) tmp4                 ! Temporary workspace
   real(r4) tfdl                 ! Temporary workspace
   real(r4) tful                 ! Temporary workspace
   real(r4) absbt(pcols)         ! Downward emission at model top
   real(r4) plol(pcols,pverp)    ! O3 pressure wghted path length
   real(r4) plos(pcols,pverp)    ! O3 path length
   real(r4) aer_mpp(pcols,pverp) ! STRAER path above kth interface level
   real(r4) co2em(pcols,pverp)   ! Layer co2 normalized planck funct. derivative
   real(r4) co2eml(pcols,pver)   ! Interface co2 normalized planck funct. deriv.
   real(r4) delt(pcols)          ! Diff t**4 mid layer to top interface
   real(r4) delt1(pcols)         ! Diff t**4 lower intrfc to mid layer
   real(r4) bk1(pcols)           ! Absrptvty for vertical quadrature
   real(r4) bk2(pcols)           ! Absrptvty for vertical quadrature
   real(r4) cldp(pcols,pverp)    ! Cloud cover with extra layer
   real(r4) fcnl(pcols,pverp)    ! Net Clear sky longwv flux
   real(r4) fclb4(pcols,-1:pver)    ! Sig t**4 for cld bottom interfc
   real(r4) fclt4(pcols,0:pver)    ! Sig t**4 for cloud top interfc
   real(r4) s(pcols,pverp,pverp) ! Flx integral sum
   real(r4) tplnka(pcols,pverp)  ! Planck fnctn temperature
   real(r4) s2c(pcols,pverp)     ! H2o cont amount
   real(r4) tcg(pcols,pverp)     ! H2o-mass-wgted temp. (Curtis-Godson approx.)
   real(r4) w(pcols,pverp)       ! H2o path
   real(r4) tplnke(pcols)        ! Planck fnctn temperature
   real(r4) h2otr(pcols,pverp)   ! H2o trnmsn for o3 overlap
   real(r4) co2t(pcols,pverp)    ! Prs wghted temperature path
   real(r4) tint(pcols,pverp)    ! Interface temperature
   real(r4) tint4(pcols,pverp)   ! Interface temperature**4
   real(r4) tlayr(pcols,pverp)   ! Level temperature
   real(r4) tlayr4(pcols,pverp)  ! Level temperature**4
   real(r4) plh2ob(nbands,pcols,pverp)! Pressure weighted h2o path with 
                                      !    Hulst-Curtis-Godson temp. factor 
                                      !    for H2O bands 
   real(r4) wb(nbands,pcols,pverp)    ! H2o path length with 
                                      !    Hulst-Curtis-Godson temp. factor 
                                      !    for H2O bands 

   real(r4) cld0                 ! previous cloud amt (for max overlap)
   real(r4) cld1                 ! next cloud amt (for max overlap)
   real(r4) emx(0:pverp)         ! Emissivity factors (max overlap)
   real(r4) emx0                 ! Emissivity factors for BCs (max overlap)
   real(r4) trans                ! 1 - emis
   real(r4) asort(pver)          ! 1 - cloud amounts to be sorted for max ovrlp.
   real(r4) atmp                 ! Temporary storage for sort when nxs = 2
   real(r4) maxcld(pcols)        ! Maximum cloud at any layer

   integer indx(pcols)       ! index vector of gathered array values
   integer indxmx(pcols+1,pverp)! index vector of gathered array values
!   integer indxmx(pcols,pverp)! index vector of gathered array values
!    (max overlap)
   integer nrgn(pcols)       ! Number of max overlap regions at longitude
   integer npts              ! number of values satisfying some criterion
   integer ncolmx(pverp)     ! number of columns with clds in region
   integer kx1(pcols,pverp)  ! Level index for top of max-overlap region
   integer kx2(pcols,0:pverp)! Level index for bottom of max-overlap region
   integer kxs(0:pverp,pcols,pverp)! Level indices for cld layers sorted by cld()
!    in descending order
   integer nxs(pcols,pverp)  ! Number of cloudy layers between kx1 and kx2
   integer nxsk              ! Number of cloudy layers between (kx1/kx2)&k
   integer ksort(0:pverp)    ! Level indices of cloud amounts to be sorted
!    for max ovrlp. calculation
   integer ktmp              ! Temporary storage for sort when nxs = 2

  real(r4) aer_trn_ttl(pcols,pverp,pverp,bnd_nbr_LW) ! [fraction] Total
!                               ! transmission between interfaces k1 and k2  
!
! Pointer variables to 3d structures
!
   real(r4), pointer :: abstot(:,:,:)
   real(r4), pointer :: absnxt(:,:,:)
   real(r4), pointer :: emstot(:,:)

!
! Trace gas variables
!
   real(r4) ucfc11(pcols,pverp)  ! CFC11 path length
   real(r4) ucfc12(pcols,pverp)  ! CFC12 path length
   real(r4) un2o0(pcols,pverp)   ! N2O path length
   real(r4) un2o1(pcols,pverp)   ! N2O path length (hot band)
   real(r4) uch4(pcols,pverp)    ! CH4 path length
   real(r4) uco211(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r4) uco212(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r4) uco213(pcols,pverp)  ! CO2 9.4 micron band path length
   real(r4) uco221(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r4) uco222(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r4) uco223(pcols,pverp)  ! CO2 10.4 micron band path length
   real(r4) bn2o0(pcols,pverp)   ! pressure factor for n2o
   real(r4) bn2o1(pcols,pverp)   ! pressure factor for n2o
   real(r4) bch4(pcols,pverp)    ! pressure factor for ch4
   real(r4) uptype(pcols,pverp)  ! p-type continuum path length
   real(r4) abplnk1(14,pcols,pverp)  ! non-nearest layer Plack factor
   real(r4) abplnk2(14,pcols,pverp)  ! nearest layer factor
!
!
!-----------------------------------------------------------------------
!
!
! Set pointer variables
!
   abstot => abstot_3d(:,:,:,lchnk)
   absnxt => absnxt_3d(:,:,:,lchnk)
   emstot => emstot_3d(:,:,lchnk)
!
! accumulate mass path from top of atmosphere
!
  call aer_pth(aer_mass, aer_mpp, ncol)

!
! Calculate some temperatures needed to derive absorptivity and
! emissivity, as well as some h2o path lengths
!
   call radtpl(lchnk   ,ncol    ,                            &
               tnm     ,lwupcgs ,qnm     ,pint    ,plco2   ,plh2o   , &
               tplnka  ,s2c     ,tcg     ,w       ,tplnke  , &
               tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    , &
               piln    ,plh2ob  ,wb      )


   if (doabsems) then
!
! Compute ozone path lengths at frequency of a/e calculation.
!
      call radoz2(lchnk, ncol, o3vmr   ,pint    ,plol    ,plos, ntoplw    )
!
! Compute trace gas path lengths
!
      call trcpth(lchnk   ,ncol    ,                            &
                  tnm     ,pint    ,cfc11   ,cfc12   ,n2o     , &
                  ch4     ,qnm     ,ucfc11  ,ucfc12  ,un2o0   , &
                  un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                  uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
                  bch4    ,uptype  )

!     Compute transmission through STRAER absorption continuum
      call aer_trn(aer_mpp, aer_trn_ttl)

!
!
! Compute total emissivity:
!
      call radems(lchnk   ,ncol    ,                            &
                  s2c     ,tcg     ,w       ,tplnke  ,plh2o   , &
                  pint    ,plco2   ,tint    ,tint4   ,tlayr   , &
                  tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  , &
                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  , &
                  bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  , &
                  co2t    ,h2otr   ,abplnk1 ,abplnk2 ,emstot  , &
                  plh2ob  ,wb      , &
                  aer_trn_ttl)
!
! Compute total absorptivity:
!
      call radabs(lchnk   ,ncol    ,                            &
                  pmid    ,pint    ,co2em   ,co2eml  ,tplnka  , &
                  s2c     ,tcg     ,w       ,h2otr   ,plco2   , &
                  plh2o   ,co2t    ,tint    ,tlayr   ,plol    , &
                  plos    ,pmln    ,piln    ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  , &
                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  , &
                  bn2o0   ,bn2o1   ,bch4    ,abplnk1 ,abplnk2 , &
                  abstot  ,absnxt  ,plh2ob  ,wb      , &
                  aer_mpp ,aer_trn_ttl)

   end if
!
! Compute sums used in integrals (all longitude points)
!
! Definition of bk1 & bk2 depends on finite differencing.  for
! trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
! layers only.
!
! delt=t**4 in layer above current sigma level km.
! delt1=t**4 in layer below current sigma level km.
!
   do i=1,ncol
      delt(i) = tint4(i,pver) - tlayr4(i,pverp)
      delt1(i) = tlayr4(i,pverp) - tint4(i,pverp)
      s(i,pverp,pverp) = stebol*(delt1(i)*absnxt(i,pver,1) + delt (i)*absnxt(i,pver,4))
      s(i,pver,pverp)  = stebol*(delt (i)*absnxt(i,pver,2) + delt1(i)*absnxt(i,pver,3))
   end do
   do k=ntoplw,pver-1
      do i=1,ncol
         bk2(i) = (abstot(i,k,pver) + abstot(i,k,pverp))*0.5
         bk1(i) = bk2(i)
         s(i,k,pverp) = stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
      end do
   end do
!
! All k, km>1
!
   do km=pver,ntoplw+1,-1
      do i=1,ncol
         delt(i)  = tint4(i,km-1) - tlayr4(i,km)
         delt1(i) = tlayr4(i,km) - tint4(i,km)
      end do
!CSD$ PARALLEL DO PRIVATE( bk1, bk2 )
      do k=pverp,ntoplw,-1
         if (k == km) then
            do i=1,ncol
               bk2(i) = absnxt(i,km-1,4)
               bk1(i) = absnxt(i,km-1,1)
            end do
         else if (k == km-1) then
            do i=1,ncol
               bk2(i) = absnxt(i,km-1,2)
               bk1(i) = absnxt(i,km-1,3)
            end do
         else
            do i=1,ncol
               bk2(i) = (abstot(i,k,km-1) + abstot(i,k,km))*0.5
               bk1(i) = bk2(i)
            end do
         end if
         do i=1,ncol
            s(i,k,km) = s(i,k,km+1) + stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
         end do
      end do
!CSD$ END PARALLEL DO
   end do
!
! Computation of clear sky fluxes always set first level of fsul
!
   do i=1,ncol
      fsul(i,pverp) = lwupcgs(i)
   end do
!
! Downward clear sky fluxes store intermediate quantities in down flux
! Initialize fluxes to clear sky values.
!
   do i=1,ncol
      tmp(i) = fsul(i,pverp) - stebol*tint4(i,pverp)
      fsul(i,ntoplw) = fsul(i,pverp) - abstot(i,ntoplw,pverp)*tmp(i) + s(i,ntoplw,ntoplw+1)
      fsdl(i,ntoplw) = stebol*(tplnke(i)**4)*emstot(i,ntoplw)
   end do
!
! fsdl(i,pverp) assumes isothermal layer
!
   do k=ntoplw+1,pver
      do i=1,ncol
         fsul(i,k) = fsul(i,pverp) - abstot(i,k,pverp)*tmp(i) + s(i,k,k+1)
         fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k) - (s(i,k,ntoplw+1) - s(i,k,k+1))
      end do
   end do
!
! Store the downward emission from level 1 = total gas emission * sigma
! t**4.  fsdl does not yet include all terms
!
   do i=1,ncol
      absbt(i) = stebol*(tplnke(i)**4)*emstot(i,pverp)
      fsdl(i,pverp) = absbt(i) - s(i,pverp,ntoplw+1)
   end do

   do k = 1,pverp
      do i = 1,ncol
         fcnl(i,k) = fsul(i,k) - fsdl(i,k)
      end do
   end do
!
!----------------------------------------------------------------------
! Modifications for clouds -- max/random overlap assumption
!
! The column is divided into sets of adjacent layers, called regions,
!   in which the clouds are maximally overlapped.  The clouds are
!   randomly overlapped between different regions.  The number of
!   regions in a column is set by nmxrgn, and the range of pressures
!   included in each region is set by pmxrgn.  The max/random overlap
!   can be written in terms of the solutions of random overlap with
!   cloud amounts = 1.  The random overlap assumption is equivalent to
!   setting the flux boundary conditions (BCs) at the edges of each region
!   equal to the mean all-sky flux at those boundaries.  Since the
!   emissivity array for propogating BCs is only computed for the
!   TOA BC, the flux BCs elsewhere in the atmosphere have to be formulated
!   in terms of solutions to the random overlap equations.  This is done
!   by writing the flux BCs as the sum of a clear-sky flux and emission
!   from a cloud outside the region weighted by an emissivity.  This
!   emissivity is determined from the location of the cloud and the
!   flux BC.
!
! Copy cloud amounts to buffer with extra layer (needed for overlap logic)
!
   cldp(:ncol,ntoplw:pver) = cld(:ncol,ntoplw:pver)
   cldp(:ncol,pverp) = 0.0
!
!
! Select only those locations where there are no clouds
!    (maximum cloud fraction <= 1.e-3 treated as clear)
!    Set all-sky fluxes to clear-sky values.
!
   maxcld(1:ncol) = maxval(cldp(1:ncol,ntoplw:pver),dim=2)

   npts = 0
   do i=1,ncol
      if (maxcld(i) < cldmin) then
         npts = npts + 1
         indx(npts) = i
      end if
   end do

   do ii = 1, npts
      i = indx(ii)
      do k = ntoplw, pverp
         fdl(i,k) = fsdl(i,k)
         ful(i,k) = fsul(i,k)
      end do
   end do
!
! Select only those locations where there are clouds
!
   npts = 0
   do i=1,ncol
      if (maxcld(i) >= cldmin) then
         npts = npts + 1
         indx(npts) = i
      end if
   end do

!
! Initialize all-sky fluxes. fdl(i,1) & ful(i,pverp) are boundary conditions
!
   do ii = 1, npts
      i = indx(ii)
      fdl(i,ntoplw) = fsdl(i,ntoplw)
      fdl(i,pverp)  = 0.0
      ful(i,ntoplw) = 0.0
      ful(i,pverp)  = fsul(i,pverp)
      do k = ntoplw+1, pver
         fdl(i,k) = 0.0
         ful(i,k) = 0.0
      end do
!
! Initialize Planck emission from layer boundaries
!
      do k = ntoplw, pver
         fclt4(i,k-1) = stebol*tint4(i,k)
         fclb4(i,k-1) = stebol*tint4(i,k+1)
      enddo
      fclb4(i,ntoplw-2) =  stebol*tint4(i,ntoplw)
      fclt4(i,pver)     = stebol*tint4(i,pverp)
!
! Initialize indices for layers to be max-overlapped
!
      do irgn = 0, nmxrgn(i)
         kx2(i,irgn) = ntoplw-1
      end do
      nrgn(i) = 0
   end do

!----------------------------------------------------------------------
! INDEX CALCULATIONS FOR MAX OVERLAP

!CSD$ PARALLEL DO PRIVATE( ii, ilon, irgn, n, k1, k2, k, indxmx, ncolmx, iimx, i, ksort ) &
!CSD$ PRIVATE( asort, ktmp, atmp, km1, km4, k3, emx0, nxsk, emx, cld0 ) &
!CSD$ PRIVATE( tmp4, l1, tmp3, tfdl, l, cld1, trans, km3, tful )
   do ii = 1, npts
      ilon = indx(ii)

!
! Outermost loop over regions (sets of adjacent layers) to be max overlapped
!
      do irgn = 1, nmxrgn(ilon)
!
! Calculate min/max layer indices inside region.
!
         n = 0
         if (kx2(ilon,irgn-1) < pver) then
            nrgn(ilon) = irgn
            k1 = kx2(ilon,irgn-1)+1
            kx1(ilon,irgn) = k1
            kx2(ilon,irgn) = 0
!cdir novector
            do k2 = pver, k1, -1
               if (pmid(ilon,k2) <= pmxrgn(ilon,irgn)) then
                  kx2(ilon,irgn) = k2
                  exit
               end if
            end do
!
! Identify columns with clouds in the given region.
!
!cdir novector
            do k = k1, k2
               if (cldp(ilon,k) >= cldmin) then
                  n = n+1
                  indxmx(n,irgn) = ilon
                  exit
               endif
            end do
         endif
         ncolmx(irgn) = n
!
! Dummy value for handling clear-sky regions
!
         indxmx(ncolmx(irgn)+1,irgn) = ncol+1
!
! Outer loop over columns with clouds in the max-overlap region
!
         do iimx = 1, ncolmx(irgn)
            i = indxmx(iimx,irgn)
!
! Sort cloud areas and corresponding level indices.
!
            n = 0
!cdir novector
            do k = kx1(i,irgn),kx2(i,irgn)
               if (cldp(i,k) >= cldmin) then
                  n = n+1
                  ksort(n) = k
!
! We need indices for clouds in order of largest to smallest, so
!    sort 1-cld in ascending order
!
                  asort(n) = 1.0-cldp(i,k)
               end if
            end do
            nxs(i,irgn) = n

!
! If nxs(i,irgn) eq 1, no need to sort.
! If nxs(i,irgn) eq 2, sort by swapping if necessary
! If nxs(i,irgn) ge 3, sort using local sort routine
!
            if (nxs(i,irgn) == 2) then
               if (asort(2) < asort(1)) then
                  ktmp = ksort(1)
                  ksort(1) = ksort(2)
                  ksort(2) = ktmp

                  atmp = asort(1)
                  asort(1) = asort(2)
                  asort(2) = atmp
               endif
            else if (nxs(i,irgn) >= 3) then
               call quick_sort(asort(1:nxs(i,irgn)),ksort(1:nxs(i,irgn)))
            endif

!cdir novector
            do l = 1, nxs(i,irgn)
               kxs(l,i,irgn) = ksort(l)
            end do
!
! End loop over longitude i for fluxes
!
         end do
!
! End loop over regions irgn for max-overlap
!
      end do
!
!----------------------------------------------------------------------
! DOWNWARD FLUXES:
! Outermost loop over regions (sets of adjacent layers) to be max overlapped
!
      do irgn = 1, nmxrgn(ilon)
!
! Compute clear-sky fluxes for regions without clouds
!
         iimx = 1
         if (ilon < indxmx(iimx,irgn) .and. irgn <= nrgn(ilon)) then
!
! Calculate emissivity so that downward flux at upper boundary of region
!    can be cast in form of solution for downward flux from cloud above
!    that boundary.  Then solutions for fluxes at other levels take form of
!    random overlap expressions.  Try to locate "cloud" as close as possible
!    to TOA such that the "cloud" pseudo-emissivity is between 0 and 1.
!
            k1 = kx1(ilon,irgn)
!cdir novector
            do km1 = ntoplw-2, k1-2
               km4 = km1+3
               k2 = k1
               k3 = k2+1
               tmp(ilon) = s(ilon,k2,min(k3,pverp))*min(1,pverp2-k3)
               emx0 = (fdl(ilon,k1)-fsdl(ilon,k1))/ &
                      ((fclb4(ilon,km1)-s(ilon,k2,km4)+tmp(ilon))- fsdl(ilon,k1))
               if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
            end do
            km1 = min(km1,k1-2)
!cdir novector
            do k2 = kx1(ilon,irgn)+1, kx2(ilon,irgn)+1
               k3 = k2+1
               tmp(ilon) = s(ilon,k2,min(k3,pverp))*min(1,pverp2-k3)
               fdl(ilon,k2) = (1.0-emx0)*fsdl(ilon,k2) + &
                               emx0*(fclb4(ilon,km1)-s(ilon,k2,km4)+tmp(ilon))
            end do
         else if (ilon==indxmx(iimx,irgn) .and. iimx<=ncolmx(irgn)) then
            iimx = iimx+1
         end if
!
! Outer loop over columns with clouds in the max-overlap region
!
         do iimx = 1, ncolmx(irgn)
            i = indxmx(iimx,irgn)

!
! Calculate emissivity so that downward flux at upper boundary of region
!    can be cast in form of solution for downward flux from cloud above that
!    boundary.  Then solutions for fluxes at other levels take form of
!    random overlap expressions.  Try to locate "cloud" as close as possible
!    to TOA such that the "cloud" pseudo-emissivity is between 0 and 1.
!
            k1 = kx1(i,irgn)
!cdir novector
            do km1 = ntoplw-2,k1-2
               km4 = km1+3
               k2 = k1
               k3 = k2 + 1
               tmp(i) = s(i,k2,min(k3,pverp))*min(1,pverp2-k3)
               tmp2(i) = s(i,k2,min(km4,pverp))*min(1,pverp2-km4)
               emx0 = (fdl(i,k1)-fsdl(i,k1))/((fclb4(i,km1)-tmp2(i)+tmp(i))-fsdl(i,k1))
               if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
            end do
            km1 = min(km1,k1-2)
            ksort(0) = km1 + 1
!
! Loop to calculate fluxes at level k
!
            nxsk = 0
            do k = kx1(i,irgn), kx2(i,irgn)
!
! Identify clouds (largest to smallest area) between kx1 and k
!    Since nxsk will increase with increasing k up to nxs(i,irgn), once
!    nxsk == nxs(i,irgn) then use the list constructed for previous k
!
               if (nxsk < nxs(i,irgn)) then
                  nxsk = 0
!cdir novector
                  do l = 1, nxs(i,irgn)
                     k1 = kxs(l,i,irgn)
                     if (k >= k1) then
                        nxsk = nxsk + 1
                        ksort(nxsk) = k1
                     endif
                  end do
               endif
!
! Dummy value of index to insure computation of cloud amt is valid for l=nxsk+1
!
               ksort(nxsk+1) = pverp
!
! Initialize iterated emissivity factors
!
!cdir novector
               do l = 1, nxsk
                  emx(l) = emis(i,ksort(l))
               end do
!
! Initialize iterated emissivity factor for bnd. condition at upper interface
!
               emx(0) = emx0
!
! Initialize previous cloud amounts
!
               cld0 = 1.0
!
! Indices for flux calculations
!
               k2 = k+1
               k3 = k2+1
               tmp4 = s(i,k2,min(k3,pverp))*min(1,pverp2-k3)
!
! Special case nxsk == 0
!
               if ( nxsk == 0 ) then
                  fdl(i,k2) = fdl(i,k2)+fsdl(i,k2)
                  if ( emx0 /= 0.0_r4 ) then
                     km1 = ksort(0)-1
                     km4 = km1+3
                     fdl(i,k2) = fdl(i,k2)+emx0* &
                                 (fclb4(i,km1)-(s(i,k2,min(km4,pverp))*min(1,pverp2-km4))+tmp4-fsdl(i,k2))
                  end if
                  cycle
               end if                   ! nxsk == 0

!
! Loop over number of cloud levels inside region (biggest to smallest cld area)
!
!cdir novector
                do l1 = 0, nxsk
                   km1 = ksort(l1)-1
                   km4 = km1+3
                   tmp3(l1) = fclb4(i,km1)-(s(i,k2,min(km4,pverp))*min(1,pverp2-km4))+tmp4-fsdl(i,k2)
               end do

               tfdl = 0.0_r4

               do l = 1, nxsk+1
!
! Calculate downward fluxes
!
                  cld1 = cldp(i,ksort(l))*min(1,nxsk+1-l)
                  if (cld0 /= cld1) then
                     tfdl = tfdl+(cld0-cld1)*fsdl(i,k2)
!cdir novector
                     do l1 = 0, l - 1
                        tfdl = tfdl+(cld0-cld1)*emx(l1)*tmp3(l1)
                     end do
                  endif
                  cld0 = cld1
!
! Multiply emissivity factors by current cloud transmissivity
!
                  if (l <= nxsk) then
                     k1 = ksort(l)
                     trans = 1.0-emis(i,k1)
!
! Ideally the upper bound on l1 would be l-1, but the sort routine
!    scrambles the order of layers with identical cloud amounts
!
!cdir novector
                     do l1 = 0, nxsk
                        if (ksort(l1) < k1) then
                           emx(l1) = emx(l1)*trans
                        endif
                     end do
                  end if
!
! End loop over number l of cloud levels
!
               end do
               fdl(i,k2) = tfdl
!
! End loop over level k for fluxes
!
            end do
!
! End loop over longitude i for fluxes
!
         end do
!
! End loop over regions irgn for max-overlap
!
      end do

!
!----------------------------------------------------------------------
! UPWARD FLUXES:
! Outermost loop over regions (sets of adjacent layers) to be max overlapped
!
      do irgn = nmxrgn(ilon), 1, -1
!
! Compute clear-sky fluxes for regions without clouds
!
         iimx = 1
         if (ilon < indxmx(iimx,irgn) .and. irgn <= nrgn(ilon)) then
!
! Calculate emissivity so that upward flux at lower boundary of region
!    can be cast in form of solution for upward flux from cloud below that
!    boundary.  Then solutions for fluxes at other levels take form of
!    random overlap expressions.  Try to locate "cloud" as close as possible
!    to surface such that the "cloud" pseudo-emissivity is between 0 and 1.
! Include allowance for surface emissivity (both numerator and denominator
!    equal 1)
!
            k1 = kx2(ilon,irgn)+1
            if (k1 < pverp) then
!cdir novector
               do km1 = pver-1,kx2(ilon,irgn),-1
                  km3 = km1+2
                  k2 = k1
                  k3 = k2+1
                  tmp(ilon) = s(ilon,k2,min(km3,pverp))* min(1,pverp2-km3)
                  emx0 = (ful(ilon,k1)-fsul(ilon,k1))/ &
                         ((fclt4(ilon,km1)+s(ilon,k2,k3)-tmp(ilon))- fsul(ilon,k1))
                  if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
               end do
               km1 = max(km1,kx2(ilon,irgn))
            else
               km1 = k1-1
               km3 = km1+2
               emx0 = 1.0
            endif

!cdir novector
            do k2 = kx1(ilon,irgn), kx2(ilon,irgn)
               k3 = k2+1
!
! If km3 == pver+2, one of the s integrals = 0 (integration limits both = p_s)
!
               tmp(ilon) = s(ilon,k2,min(km3,pverp))* min(1,pverp2-km3)
               ful(ilon,k2) =(1.0-emx0)*fsul(ilon,k2) + emx0* &
                             (fclt4(ilon,km1)+s(ilon,k2,k3)-tmp(ilon))
            end do
         else if (ilon==indxmx(iimx,irgn) .and. iimx<=ncolmx(irgn)) then
            iimx = iimx+1
         end if
!
! Outer loop over columns with clouds in the max-overlap region
!
         do iimx = 1, ncolmx(irgn)
            i = indxmx(iimx,irgn)

!
! Calculate emissivity so that upward flux at lower boundary of region
!    can be cast in form of solution for upward flux from cloud at that
!    boundary.  Then solutions for fluxes at other levels take form of
!    random overlap expressions.  Try to locate "cloud" as close as possible
!    to surface such that the "cloud" pseudo-emissivity is between 0 and 1.
! Include allowance for surface emissivity (both numerator and denominator
!    equal 1)
!
            k1 = kx2(i,irgn)+1
            if (k1 < pverp) then
!cdir novector
               do km1 = pver-1,kx2(i,irgn),-1
                  km3 = km1+2
                  k2 = k1
                  k3 = k2+1
                  tmp(i) = s(i,k2,min(km3,pverp))*min(1,pverp2-km3)
                  emx0 = (ful(i,k1)-fsul(i,k1))/((fclt4(i,km1)+s(i,k2,k3)-tmp(i))-fsul(i,k1))
                  if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
               end do
               km1 = max(km1,kx2(i,irgn))
            else
               emx0 = 1.0
               km1 = k1-1
            endif
            ksort(0) = km1 + 1

!
! Loop to calculate fluxes at level k
!
            nxsk = 0
            do k = kx2(i,irgn), kx1(i,irgn), -1
!
! Identify clouds (largest to smallest area) between k and kx2
!    Since nxsk will increase with decreasing k up to nxs(i,irgn), once
!    nxsk == nxs(i,irgn) then use the list constructed for previous k
!
               if (nxsk < nxs(i,irgn)) then
                  nxsk = 0
!cdir novector
                  do l = 1, nxs(i,irgn)
                     k1 = kxs(l,i,irgn)
                     if (k <= k1) then
                        nxsk = nxsk + 1
                        ksort(nxsk) = k1
                     endif
                  end do
               endif
!
! Dummy value of index to insure computation of cloud amt is valid for l=nxsk+1
!
               ksort(nxsk+1) = pverp
!
! Initialize iterated emissivity factors
!
!cdir novector
               do l = 1, nxsk
                  emx(l) = emis(i,ksort(l))
               end do
!
! Initialize iterated emissivity factor for bnd. condition at lower interface
!
               emx(0) = emx0
!
! Initialize previous cloud amounts
!
               cld0 = 1.0
!
! Indices for flux calculations
!
               k2 = k
               k3 = k2+1
!
! Special case nxsk == 0
!
               if ( nxsk == 0 ) then
                  ful(i,k2) = ful(i,k2)+fsul(i,k2)
                  if ( emx0 /= 0.0_r4 ) then
                     km1 = ksort(0)-1
                     km3 = km1+2
!
! If km3 == pver+2, one of the s integrals = 0 (integration limits both = p_s)
!
                     ful(i,k2) = ful(i,k2)+emx0* &
                                (fclt4(i,km1)+s(i,k2,k3)-(s(i,k2,min(km3,pverp))*min(1,pverp2-km3))-fsul(i,k2))

                  end if
                  cycle
               end if
!
! Loop over number of cloud levels inside region (biggest to smallest cld area)
!
!cdir novector
               do l1 = 0, nxsk
                  km1 = ksort(l1)-1
                  km3 = km1+2
!
! If km3 == pver+2, one of the s integrals = 0 (integration limits both = p_s)
!
                  tmp3(l1) = fclt4(i,km1)+s(i,k2,k3)-(s(i,k2,min(km3,pverp))*min(1,pverp2-km3))- fsul(i,k2)
               end do

               tful = 0.0_r4

               do l = 1, nxsk+1
!
! Calculate upward fluxes
!
                  cld1 = cldp(i,ksort(l))*min(1,nxsk+1-l)
                  if (cld0 /= cld1) then
                     tful = tful+(cld0-cld1)*fsul(i,k2)
!cdir novector
                     do l1 = 0, l - 1
                        tful = tful+(cld0-cld1)*emx(l1)*tmp3(l1)
                     end do
                  endif
                  cld0 = cld1
!
! Multiply emissivity factors by current cloud transmissivity
!
                  if (l <= nxsk) then
                     k1 = ksort(l)
                     trans = 1.0-emis(i,k1)
!
! Ideally the upper bound on l1 would be l-1, but the sort routine
!    scrambles the order of layers with identical cloud amounts
!
!cdir novector
                     do l1 = 0, nxsk
                        if (ksort(l1) > k1) then
                           emx(l1) = emx(l1)*trans
                        endif
                     end do
                  end if
!
! End loop over number l of cloud levels
!
               end do
               ful(i,k2) = tful
!
! End loop over level k for fluxes
!
            end do
!
! End loop over longitude i for fluxes
!
         end do
!
! End loop over regions irgn for max-overlap
!
      end do
!
! End outermost longitude loop
!
   end do
!CSD$ END PARALLEL DO
!
! End cloud modification loops
!
!----------------------------------------------------------------------
! All longitudes: store history tape quantities
!
   do i=1,ncol
      flwds(i) = fdl (i,pverp )
      flns(i)  = ful (i,pverp ) - fdl (i,pverp )
      flnsc(i) = fsul(i,pverp ) - fsdl(i,pverp )
      flnt(i)  = ful (i,ntoplw) - fdl (i,ntoplw)
      flntc(i) = fsul(i,ntoplw) - fsdl(i,ntoplw)
      flut(i)  = ful (i,ntoplw)
      flutc(i) = fsul(i,ntoplw)
   end do

   do k = 1,pverp
      do i = 1,ncol
         fnl(i,k) = ful(i,k) - fdl(i,k)
      end do
   end do
!
! Computation of longwave heating (J/kg/s)
!
   do k=ntoplw,pver
      do i=1,ncol
         qrl(i,k) = (ful(i,k) - fdl(i,k) - ful(i,k+1) + fdl(i,k+1))* &
                     1.E-4*gravit/((pint(i,k) - pint(i,k+1)))
      end do
   end do
! Return 0 above solution domain
   if ( ntoplw > 1 )then
      qrl(:ncol,:ntoplw-1) = 0.
   end if
!
   return
end subroutine radclwmx

