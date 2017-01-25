subroutine diffuse_scalar_edmf(f,fluxb,fluxt,sumMs, &
                          fdiff,flux,fluxmf,flux3,f2lediff,f2lediss,fwlediff,doit,massflux,flux3mf)

use grid
use vars, only: rho, rhow
use sgs, only: tkh, sgs_field_sumM, betap, betam
implicit none

! input:	
real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real sumMs(1:nx, 1:ny,nz)		! MF flux of scalar
real flux(nz),fluxmf(nz)
real flux3(nx,ny,nz),flux3mf(nx,ny,nz)
real fdiff(nz)
real f2lediff(nzm)
real f2lediss(nzm)
real fwlediff(nzm)
real var(nzm),sumMs1D(nz),tkh1D(nzm)
real,dimension(nzm) :: a, b, c, d
logical doit, massflux, ispresent
! Local
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real f0(nzm),df0(nzm),factor_xy
real r2dx,r2dy,r2dx0,r2dy0,r2dz
integer i,j,k,kb,kc,jb,jc

call t_startf ('diffuse_scalars')


if(dostatis) then
	
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = f(i,j,k)
     end do
    end do
  end do

endif

flux = 0.
flux3 = 0.
fluxmf = 0.
flux3mf = 0.
do i=1,nx
  do j=1,ny
    var = f(i,j,1:nzm)
    sumMs1D = sumMs(i,j,1:nz)
    tkh1D   = tkh(i,j,1:nzm)
    call get_abcd(i,j,betap,betam,var,sumMs1D,tkh1D,a,b,c,d, massflux, fluxb(i,j))
    call tridiag(a,b,c,d)
    flux(1) = flux(1) + fluxb(i,j) * rhow(1) / dz
    flux3(i,j,1) = fluxb(i,j)
    flux3(i,j,nz) = 0.0
    flux(nz)= 0.
    flux3(i,j,2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tkh(i,j,1:nzm-1) + tkh(i,j,2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(f(i,j,2:nzm)-f(i,j,1:nzm-1)) ) 
    if (massflux) then
      flux3mf(i,j,2:nzm) = (sumMs(i,j,2:nzm) - (betap*d(2:nzm) + betam*f(i,j,2:nzm)) * sgs_field_sumM(i,j,2:nzm,1) )
      fluxmf(2:nzm) = fluxmf(2:nzm) + rhow(2:nzm) * flux3mf(i,j,2:nzm)/dz
    end if
    flux(2:nzm) = flux(2:nzm) + rhow(2:nzm) * flux3(i,j,2:nzm)/dz
    f(i,j,1:nzm) = d(1:nzm)
  end do
end do


if(dostatis) then
	
  do k=1,nzm
    fdiff(k)=0.
    do j=1,ny
     do i=1,nx
      fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
     end do
    end do
  end do

endif

if(dostatis.and.doit) then
	
  call stat_varscalar(f,df,f0,df0,f2lediff)
  call stat_sw2(f,df,fwlediff)

  factor_xy=1./float(nx*ny)
  r2dx0=1./(2.*dx)
  r2dy0=1./(2.*dy)
  do k=1,nzm
    f2lediss(k)=0.
    kc=min(nzm,k+1)
    kb=max(1,k-1)
    r2dz=2./((kc-kb)*(adzw(k+1)+adzw(k))*dz)
    r2dx=r2dx0*sqrt((kc-kb)*dx*r2dz) ! grid anisotropy correction
    r2dy=r2dy0*sqrt((kc-kb)*dx*r2dz)
    f2lediss(k)=0.
    do j=1,ny
     jc=j+YES3D
     jb=j-YES3D
     do i=1,nx
      f2lediss(k)=f2lediss(k)-tkh(i,j,k)*( &
                       ((f(i+1,j,k)-f(i-1,j,k))*r2dx)**2+ &
                       ((f(i,jc,k)-f(i,jb,k))*r2dy)**2+ &
                       ((f(i,j,kc)-f0(kc)-f(i,j,kb)+f0(kb))*r2dz)**2 )
     end do
    end do
    f2lediss(k)=f2lediss(k)*2.*factor_xy
  end do

endif

call t_stopf ('diffuse_scalars')

end subroutine diffuse_scalar_edmf 
