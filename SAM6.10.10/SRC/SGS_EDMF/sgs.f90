module sgs

! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
! Marat Khairoutdinov, 2012

use grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s 
use params, only: dosgs,doedmf
use microphysics, only : nmicro_fields
use tracers, only : ntracers
implicit none

!----------------------------------------------------------------------
! Required definitions:

!!! prognostic scalar (need to be advected arround the grid):

! sum of betap and betam has to be one
real, parameter :: betap = 1., betam = 0.

integer, parameter :: nsgs_fields = 1   ! total number of prognostic sgs vars

real sgs_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nsgs_fields)

!!! sgs diagnostic variables that need to exchange boundary information (via MPI):

integer, parameter :: nsgs_fields_diag = 2   ! total number of diagnostic sgs vars

! diagnostic fields' boundaries:
integer, parameter :: dimx1_d=0, dimx2_d=nxp1, dimy1_d=1-YES3D, dimy2_d=nyp1

real sgs_field_diag(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm, nsgs_fields_diag)
! mass fluxes from multiplume model: sumM, sumMu, sumMv, sumMtke, sumMt, sumMmicro(:), sumMtr(:)
real sgs_field_sumM(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nz, 5+nmicro_fields+ntracers)
! note that an additional array has to be used since vertical size is different; boundary exchange has been added to task_boundary and periodic

! these fluxes hold sum of ED and MF parts of fluxes
real twsb3 (nx,ny,nz)                         ! sgs vertical flux of h/cp
real mkwsb3(nx,ny,nz,1:nmicro_fields)         ! sgs vertical flux of qx
real uwsb3(nx,ny,nz)                          ! sgs vertical flux of u
real vwsb3(nx,ny,nz)                          ! sgs vertical flux of v
real tkewsb3(nx,ny,nz)                        ! sgs vertical flux of tke
real thvflx(nx,ny,nzm)                        ! sgs buoyancy flux in TKE equation
real thvflx1D(nzm)                        ! sgs buoyancy flux in TKE equation

! MF things
real twsb3_mf (nx,ny,nz)                         ! sgs vertical flux of h/cp
real mkwsb3_mf(nx,ny,nz)                         ! sgs vertical flux of qt
real cfrac_mf1D(nz)
real frac_mf1D(nz)

logical:: advect_sgs = .true. ! advect prognostics
logical, parameter:: do_sgsdiag_bound = .true.  ! exchange boundaries for diagnostics fields

! SGS fields that output by default (if =1).
integer, parameter :: flag_sgs3Dout(nsgs_fields) = (/0/)
integer, parameter :: flag_sgsdiag3Dout(nsgs_fields_diag) = (/0,0/)

real fluxbsgs (nx, ny, 1:nsgs_fields) ! surface fluxes 
real fluxtsgs (nx, ny, 1:nsgs_fields) ! top boundary fluxes 

!!! these arrays may be needed for output statistics:

real sgswle(nz,1:nsgs_fields)  ! resolved vertical flux
real sgswsb(nz,1:nsgs_fields)  ! SGS vertical flux
real sgsadv(nz,1:nsgs_fields)  ! tendency due to vertical advection
real sgslsadv(nz,1:nsgs_fields)  ! tendency due to large-scale vertical advection
real sgsdiff(nz,1:nsgs_fields)  ! tendency due to vertical diffusion


!------------------------------------------------------------------
! internal (optional) definitions:

! make aliases for prognostic variables:

real tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
equivalence (tke(dimx1_s,dimy1_s,1),sgs_field(dimx1_s,dimy1_s,1,1))

! make aliases for diagnostic variables:

real tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
real tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
equivalence (tk(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,1))
equivalence (tkh(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,2))


real grdf_x(nzm)! grid factor for eddy diffusion in x
real grdf_y(nzm)! grid factor for eddy diffusion in y
real grdf_z(nzm)! grid factor for eddy diffusion in z

logical :: dofixedtau   ! if true, tau=600 sec
                          ! and assuming tke(z=0) = 0
logical :: fixedeps
logical :: donoplumesat
logical :: dosingleplume
logical :: doedmfpw
logical :: domffxdsflx
real :: tauneggers
real :: ctketau
real :: beta
real :: pwmin
integer :: nup
real ::eps0
real :: alphaqt
real :: alphathv
real :: Wc

real :: wstar(nx,ny)

! Local diagnostics:

real tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz)

CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm (namelist) file

subroutine sgs_setparm()

  use grid, only: case, caseid, masterproc
  use params, only: pblhfluxmin,dopblh
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder

  !======================================================================
  NAMELIST /SGS_TKE/ &
       dofixedtau,ctketau,fixedeps,tauneggers,&
       dosingleplume,beta,donoplumesat,pwmin,nup,eps0,doedmfpw,domffxdsflx,&
       alphaqt,alphathv, Wc

  NAMELIST /BNCUIODSBJCB/ place_holder

  dofixedtau = .false.
  if (dofixedtau) then
    ctketau = 700.
  else
    ctketau = 0.5
  end if
  fixedeps=.false.
  donoplumesat=.false.
  tauneggers=500. 
  dosingleplume=.false.
  doedmfpw=.false.
  domffxdsflx=.false.
  beta=0.3
  pwmin=1.4
  nup = 40
  eps0=1.0e-3
  alphaqt=0.32
  alphathv=0.58
  Wc=0.5

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  read (55,SGS_TKE,IOSTAT=ios)


  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in SGS_TKE namelist'
        call task_abort()
     end if
  end if

  if (dopblh.and.pblhfluxmin .and..not.dosgs) then
      if (masterproc) write(*,*) '****** ERROR: bad specification in SGS_TKE namelist for dopblh'
        call task_abort()
  end if
  if (doedmf.and..not.dopblh) then
      if (masterproc) write(*,*) '****** ERROR: bad specification in SGS_TKE namelist for dopblh&doedmf'
        call task_abort()
  end if

  close(55)
   ! write namelist values out to file for documentation
   if(masterproc) then
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.nml', form='formatted', position='append')
      write (unit=55,nml=SGS_TKE,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
   end if

  ! END UW ADDITION
  !======================================================================

end subroutine sgs_setparm

!----------------------------------------------------------------------
!!! Initialize sgs:


subroutine sgs_init()

  use grid, only: nrestart, dx, dy, dz, adz, masterproc
  use params, only: LES
  use vars, only : qcsgs_mf, qisgs_mf, cfrac_mf, frac_mf, twsbmf
  integer k

  if(nrestart.eq.0) then

     sgs_field = 0.
     sgs_field_sumM = 0.
     sgs_field_diag = 0.

     fluxbsgs = 0.
     fluxtsgs = 0.

     twsb3 = 0.
     twsb3_mf = 0.
     mkwsb3 = 0.
     twsbmf = 0.
     mkwsb3_mf = 0.
     uwsb3 = 0. 
     vwsb3 = 0. 
     tkewsb3 = 0.       

     qcsgs_mf=0.
     qisgs_mf=0.
     cfrac_mf=0.
     frac_mf=0.
     cfrac_mf1D=0.
     frac_mf1D=0.
     thvflx1D=0.

  end if

  if(masterproc) then
     write(*,*) 'EDMF scheme'
  end if

  if(LES) then
    do k=1,nzm
       grdf_x(k) = dx**2/(adz(k)*dz)**2
       grdf_y(k) = dy**2/(adz(k)*dz)**2
       grdf_z(k) = 1.
    end do
  else
    do k=1,nzm
       grdf_x(k) = min(16.,dx**2/(adz(k)*dz)**2)
       grdf_y(k) = min(16.,dy**2/(adz(k)*dz)**2)
       grdf_z(k) = 1.
    end do
  end if

  sgswle = 0.
  sgswsb = 0.
  sgsadv = 0.
  sgsdiff = 0.
  sgslsadv = 0.


end subroutine sgs_init

!----------------------------------------------------------------------
!!! make some initial noise in sgs:
!
subroutine setperturb_sgs(ptype)

use vars, only: q0, z
use grid, only: masterproc
integer, intent(in) :: ptype
integer i,j,k

select case (ptype)


  case(-1)

    if (masterproc) print*, 'Initial TKE=0.000'
    tke =0.001

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(k.le.4) then
            tke(i,j,k)=0.04*(5-k)
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.6.e-3) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do

  case(2)

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.0.5e-3) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.150.) then
            tke(i,j,k)=0.15*(1.-z(k)/150.)
         endif
       end do
      end do
     end do


  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.3000.) then
            tke(i,j,k)=1.-z(k)/3000.
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.6.e-3) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case default

end select

end subroutine setperturb_sgs

!----------------------------------------------------------------------
!!! Estimate Courant number limit for SGS
!

subroutine kurant_sgs(cfl)

use grid, only: dt, dx, dy, dz, adz, adzw
implicit none

real, intent(out) :: cfl

integer k
real tkhmax(nz)

cfl=0.
if (betap.ne.1.) then
! Only needed if not fully implicit
! Also, only vertical

do k = 1,nzm
 tkhmax(k) = maxval(tkh(1:nx,1:ny,k))
end do

do k=1,nzm
  cfl = max(cfl,        &
     !0.5*tkhmax(k)*grdf_z(k)*dt/(dz*adzw(k))**2, &
     !0.5*tkhmax(k)*grdf_x(k)*dt/dx**2, &
     !YES3D*0.5*tkhmax(k)*grdf_y(k)*dt/dy**2)
     2.0*tkhmax(k)*grdf_z(k)*dt/(dz*adzw(k))**2) 
end do
end if

end subroutine kurant_sgs


!----------------------------------------------------------------------
!!! compute sgs diffusion of momentum:
!
subroutine sgs_mom()

   call diffuse_mom()

end subroutine sgs_mom

!----------------------------------------------------------------------
!!! compute sgs diffusion of scalars:
!
subroutine sgs_scalars()

  use vars
  use microphysics
  use tracers
  use params, only: dotracers,dosurface,doedmf
  use grid, only: nx, ny, nz, adz, dz 
  implicit none

    real dummy(nz)
    real dummy3(nx,ny,nz), dummy3b(nx,ny,nz)
    real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
    integer i,j,k

      dummy3 = sgs_field_sumM(1:nx,1:ny,1:nz,5)

      call diffuse_scalar_edmf(t,fluxbt,fluxtt,dummy3,tdiff,twsb,twsbmf,twsb3, &
                           t2lediff,t2lediss,twlediff,.true.,doedmf,twsb3_mf)
    
      if (.not.dosurface) ustar=0.
      if(advect_sgs) then
          do i=1,nx
          do j=1,ny
            ! compute a downward flux of tke assuming tke=0 at sfc and K_sfc=K_1
            tkewsb3(i,j,1) = - 2.*tk(i,j,1)/adz(1)/dz * (tke(i,j,1) - 3.75*ustar(i,j)**2 - 0.2*wstar(i,j)**2)  
          end do
          end do
         fluxbtmp(1:nx,1:ny) = tkewsb3(1:nx,1:ny,1)
         dummy3 = sgs_field_sumM(1:nx,1:ny,1:nz,4)
         call diffuse_scalar(tke,fluxbtmp,fzero,dummy3,sgsdiff(1:nz,1),sgswsb,tkewsb3, &
                                    dummy,dummy,dummy,.false.,.false.)
      end if


!
!    diffusion of microphysics prognostics:
!
      call micro_flux()

      total_water_evap = total_water_evap - total_water()

      do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k)
           dummy3 =  sgs_field_sumM(1:nx,1:ny,1:nz,5+k)
           if (flag_precip(k).eq.0) then
           call diffuse_scalar_edmf(micro_field(:,:,:,k),fluxbtmp,fluxttmp,dummy3, &
                mkdiff(:,k),mkwsb(:,k),mkwsbmf,mkwsb3(:,:,:,k), dummy,dummy,dummy,.false.,doedmf,mkwsb3_mf)
           else
           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp,dummy3, &
                mkdiff(:,k),mkwsb(:,k),mkwsb3(:,:,:,k),dummy,dummy,dummy,.false.,.false.)
           end if
       end if
      end do

      total_water_evap = total_water_evap + total_water()

 ! diffusion of tracers:

      if(dotracers) then

        call tracers_flux()

        do k = 1,ntracers

          fluxbtmp = fluxbtr(:,:,k)
          fluxttmp = fluxttr(:,:,k)
          dummy3b = sgs_field_sumM(1:nx,1:ny,1:nz,5+nmicro_fields+k)
          call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp,dummy3b, &
               trdiff(:,k),trwsb(:,k),dummy3, &
               dummy,dummy,dummy,.false.,.false.)
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)

        end do

      end if



end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc()

   use grid, only: nstep,dt,icycle,dompi,nzm,ncycle
   use params, only: dosmoke, dotracers, dosgs, dopblh
   use vars
   use microphysics
   use tracers



   integer :: i,j,k
   real :: pwmax,coef1


     pw_xyinst=0.
     do k=1,nzm
        coef1 = rho(k)*dz*adz(k)
        do i=1,nx
        do j=1,ny
         pw_xyinst(i,j) = pw_xyinst(i,j)+qv(i,j,k)*coef1 ! not weighted by dtn/dt
                                                         !since used in edmf each cycle
        end do
        end do
     end do

     
      ! get maximum of instantaneous PW distribution
     pw_globalmax=0.0
     if (doedmfpw.and.doedmf.and.icycle.eq.1) then

       pw_xyinst=0.
       do k=1,nzm
          coef1 = rho(k)*dz*adz(k)
          do i=1,nx
          do j=1,ny
           pw_xyinst(i,j) = pw_xyinst(i,j)+qv(i,j,k)*coef1 ! not weighted by dtn/dt
                                                           !since used in edmf
                                                         !each cycle
          end do
          end do
       end do
       pwmax=maxval(pw_xyinst)
       if (dompi) then
         call task_max_real(pwmax,pw_globalmax,1)
       else
         pw_globalmax=pwmax
       end if
 
     elseif (icycle.eq.ncycle .and. .not.doedmf) then
    
       pw_xyinst=0.
       do k=1,nzm
          coef1 = rho(k)*dz*adz(k)
          do i=1,nx
          do j=1,ny
           pw_xyinst(i,j) = pw_xyinst(i,j)+qv(i,j,k)*coef1 ! not weighted by dtn/dt
                                                           !since used in edmf
                                                         !each cycle
          end do
          end do
       end do
     

     end if



     if (dopblh) call get_pblh()


     sgs_field_sumM = 0.
     if (dosgs.and.doedmf) call edmf()

!    SGS TKE equation:

     if(dosgs) call tke_full()

! add exchange of sumM fields which isn't done in task_boundaries 

 !if (dompi) then
!
!   do i = 1,5+nmicro_fields+ntracers
!       call task_exchange(sgs_field_sumM(:,:,:,i),dimx1_d,dimx2_d,dimy1_d,dimy2_d,nz, &
!                                                             1+dimx1_d,dimx2_d-nx,YES3D+dimy1_d,1-YES3D+dimy2_d-ny,&
!                                                             4+nsgs_fields+nsgs_fields_diag+nmicro_fields+ntracers+i)
!   end do
! else
!
!   do i = 1,5+ntracers+nmicro_fields
!    if(dosgs.and.do_sgsdiag_bound) &
!     call bound_exchange(sgs_field_sumM(:,:,:,i),dimx1_d,dimx2_d,dimy1_d,dimy2_d,nz, &
!                                                           1+dimx1_d,dimx2_d-nx,YES3D+dimy1_d,1-YES3D+dimy2_d-ny,&
!                                                           4+nsgs_fields+nsgs_fields_diag+nmicro_fields+ntracers+i)
!   end do
!
! end if


end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine sgs_diagnose()
! None 

end subroutine sgs_diagnose


!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!
subroutine sgs_statistics()
  
  use vars
  use hbuffer, only: hbuf_put, hbuf_avg_put
  use params, only : lcond

  real tmp(2), factor_xy 
  real tkz(nzm), tkhz(nzm)
  integer i,j,k,n
  character(LEN=6) :: statname  !bloss: for conditional averages

  call t_startf ('sgs_statistics')

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
    tkz(k) = 0.
    tkhz(k) = 0.
    do j=1,ny
    do i=1,nx
      tkz(k)=tkz(k)+tk(i,j,k)
      tkhz(k)=tkhz(k)+tkh(i,j,k)
    end do
    end do
  end do

  call hbuf_avg_put('TKES',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)

  call hbuf_put('TK',tkz,factor_xy)
  call hbuf_put('TKH',tkhz,factor_xy)

!---------------------------------------------------------
! SGS TKE Budget:

         call hbuf_put('ADVTRS',sgsdiff(:,1),factor_xy/dtn)
         call hbuf_put('BUOYAS',tkesbbuoy,factor_xy)
         call hbuf_put('SHEARS',tkesbshear,factor_xy)
         call hbuf_put('DISSIPS',tkesbdiss,factor_xy)
         call hbuf_put('TKEWTHV',thvflx1D,factor_xy)

!---------------------------------------------------------
! MF Massflux properties:

  call hbuf_put('FRAC_MF',frac_mf1D,factor_xy)
  call hbuf_put('CFRAC_MF',cfrac_mf1D,factor_xy)


  call t_stopf ('sgs_statistics')

end subroutine sgs_statistics

!----------------------------------------------------------------------
! called when stepout() called

subroutine sgs_print()

 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)

end subroutine sgs_print

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics 
!
subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
use vars
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,sgscount


   count = count + 1
   sgscount = sgscount + 1
   namelist(count) = 'CFRAC_MF'
   deflist(count) = 'Cloud fraction from plumes'
   unitlist(count) = '1'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   sgscount = sgscount + 1
   namelist(count) = 'FRAC_MF'
   deflist(count) = 'Plume fractional area'
   unitlist(count) = '1'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   sgscount = sgscount + 1
   namelist(count) = 'TKEWTHV'
   deflist(count) = 'Buoyancy flux in tke equation'
   unitlist(count) = 'K m/s'
   status(count) = 1
   average_type(count) = 0

end subroutine sgs_hbuf_init


end module sgs



