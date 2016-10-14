
subroutine diffuse_mom2D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk, grdf_x, grdf_z, betap, betam, uwsb3, vwsb3, sgs_field_sumM
use params, only: docolumn, dowallx
implicit none

real rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
real dxz,dzx

integer i,j,k,ic,ib,kc,kcu
real tkx, tkz, rhoi, iadzw, iadz
real fu(0:nx,1,nz),fv(0:nx,1,nz),fw(0:nx,1,nz)
real var(nzm), sumMs1D(nz), tk1D(nzm)

real,dimension(nzm) :: a, b, c, d
logical :: massflux

rdx2=1./dx/dx
rdx25=0.25*rdx2

dxz=dx/dz

j=1


if(.not.docolumn) then

if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
         v(0,j,k) = v(1,j,k)
         w(0,j,k) = w(1,j,k)
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
         v(nx+1,j,k) = v(nx,j,k)
         w(nx+1,j,k) = w(nx,j,k)
    end do
  end if

end if

end if

uwsb3 = 0.
uwsb = 0.
massflux=.false.
do i=1,nx
  ib=i-1
    var = u(i,j,1:nzm)
    sumMs1D = 0.5*(sgs_field_sumM(i,j,1:nz,2)+sgs_field_sumM(ib,j,1:nz,2) )
    tk1D = 0.5*(tk(i,j,1:nzm)+tk(ib,j,1:nzm))
    call get_abcd(i,j,betap,betam,var,sumMs1D,tk1D,a,b,c,d, massflux, fluxbu(i,j))
    call tridiag(a,b,c,d)
    uwsb(1) = uwsb(1) + fluxbu(i,j) * rhow(1) / dz
    uwsb3(i,j,1) = fluxbu(i,j)
    uwsb3(i,j,nz) = 0.0
    uwsb(nz)= 0.
    uwsb3(i,j,2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.25*(tk(i,j,1:nzm-1) + tk(i,j,2:nzm)+tk(ib,j,1:nzm-1) + tk(ib,j,2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(u(i,j,2:nzm)-u(i,j,1:nzm-1)) ) &
                           +(0.5*(sgs_field_sumM(i,j,2:nzm,2)+sgs_field_sumM(ib,j,2:nzm,2)) - (betap*d(2:nzm) + betam*u(i,j,2:nzm)) * sgs_field_sumM(i,j,2:nzm,1) )
    uwsb(2:nzm) = uwsb(2:nzm) + rhow(2:nzm) * uwsb3(i,j,2:nzm) / dz
    dudtdiff(i,j,:,na) = (d-u(i,j,:))/dtn
    dudt(i,j,:,na) = dudt(i,j,:,na) + dudtdiff(i,j,:,na)
end do

vwsb3 = 0.
vwsb = 0.
massflux=.false.
do i=1,nx
    var = v(i,j,1:nzm)
    sumMs1D = sgs_field_sumM(i,j,1:nz,3)
    tk1D = tk(i,j,1:nzm)
    call get_abcd(i,j,betap,betam,var,sumMs1D,tk1D,a,b,c,d, massflux, fluxbv(i,j))
    call tridiag(a,b,c,d)
    vwsb(1) = vwsb(1) + fluxbv(i,j) * rhow(1) / dz
    vwsb3(i,j,1) = fluxbv(i,j)
    vwsb3(i,j,nz) = 0.0
    vwsb(nz)= 0.
    vwsb3(i,j,2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.5*(tk(i,j,1:nzm-1) + tk(i,j,2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(v(i,j,2:nzm)-v(i,j,1:nzm-1)) ) &
                           +(sgs_field_sumM(i,j,2:nzm,3) - (betap*d(2:nzm) + betam*v(i,j,2:nzm)) * sgs_field_sumM(i,j,2:nzm,1) )
    vwsb(2:nzm) = vwsb(2:nzm) + rhow(2:nzm) * vwsb3(i,j,2:nzm) / dz
    dvdtdiff(i,j,:,na) = (d-v(i,j,:))/dtn
    dvdt(i,j,:,na) = dvdt(i,j,:,na) + dvdtdiff(i,j,:,na)
end do


end subroutine diffuse_mom2D


