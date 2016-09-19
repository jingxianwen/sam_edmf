
subroutine diffuse_mom3D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk, grdf_x, grdf_y, grdf_z, betap, betam
use params, only: docolumn, dowallx, dowally
implicit none

real rdx2,rdy2,rdz2,rdz,rdx25,rdy25
real rdx21,rdy21,rdx251,rdy251,rdz25
real dxy,dxz,dyx,dyz,dzx,dzy

integer i,j,k,ic,ib,jb,jc,kc,kcu
real tkx, tky, tkz, rhoi, iadzw, iadz
real fu(0:nx,0:ny,nz),fv(0:nx,0:ny,nz),fw(0:nx,0:ny,nz)

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)

rdx25=0.25*rdx2
rdy25=0.25*rdy2

dxy=dx/dy
dxz=dx/dz
dyx=dy/dx
dyz=dy/dz

!-----------------------------------------
if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
     do j=1,ny
         v(0,j,k) = v(1,j,k)
         w(0,j,k) = w(1,j,k)
     end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
     do j=1,ny
         v(nx+1,j,k) = v(nx,j,k)
         w(nx+1,j,k) = w(nx,j,k)
     end do
    end do
  end if

end if

if(dowally) then

  if(rank.lt.nsubdomains_x) then
    do k=1,nzm
       do i=1,nx
         u(i,1-YES3D,k) = u(i,1,k)
         w(i,1-YES3D,k) = w(i,1,k)
       end do
    end do
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    do k=1,nzm
       do i=1,nx
         u(i,ny+YES3D,k) = u(i,ny,k)
         w(i,ny+YES3D,k) = w(i,ny,k)
       end do
    end do
  end if

end if

uwsb3 = 0.
uwsb = 0.
do i=1,nx
  ib=i-1
  do j=1,ny
    call get_abcd(i,j,betap,betam,u(i,j,:),0.5*(sgs_field_sumM(i,j,:,2)+sgs_field_sumM(ib,j,:,2) ),0.5*(tk(i,j,:)+tk(ib,j,:)),a,b,c,d, massflux, fluxbu(i,j))
    call tridiag(a,b,c,d)
    uwsb(1) = uwsb(1) + fluxbu(i,j)
    uwsb3(i,j,1) = fluxbu(i,j)
    uwsb3(i,j,nz) = 0.0
    uwsb(nz)= 0.
    uwsb3(i,j,2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.25*(tk(i,j,1:nzm-1) + tk(i,j,2:nzm)+tk(ib,j,1:nzm-1) + tk(ib,j,2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(u(i,j,2:nzm)-u(i,j,1:nzm-1)) ) &
                           +(0.5*(sgs_field_sumM(i,j,2:nzm,2)+sgs_field_sumM(ib,j,2:nzm),2) - (betap*d(2:nzm) + betam*u(i,j,2:nzm)) * sgs_field_sumM(i,j,2:nzm,1) )
    uwsb(2:nzm) = uwsb(2:nzm) + rhow(2:nzm) * uwsb3(i,j,2:nzm)
    dudt(i,j,:,na) = dudt(i,j,:,na) + (d-u(i,j,:))/dtn
  end do
end do

vwsb3 = 0.
vwsb = 0.
do i=1,nx
  do j=1,ny
  jb=j-1
    call get_abcd(i,j,betap,betam,v(i,j,:),0.5*(sgs_field_sumM(i,j,:,3)+sgs_field_sumM(i,jb,:,3) ),0.5*(tk(i,j,:)+tk(i,jb,:)),a,b,c,d, massflux, fluxbv(i,j))
    call tridiag(a,b,c,d)
    vwsb(1) = vwsb(1) + fluxbv(i,j)
    vwsb3(i,j,1) = fluxbv(i,j)
    vwsb3(i,j,nz) = 0.0
    vwsb(nz)= 0.
    vwsb3(i,j,2:nzm) =  (-1.) /adzw(2:nzm)/dz *         0.25*(tk(i,j,1:nzm-1) + tk(i,j,2:nzm)+tk(i,jb,1:nzm-1) + tk(i,jb,2:nzm)) *                  &
                           (betap* (d(2:nzm)- d(1:nzm-1))+betam*(v(i,j,2:nzm)-v(i,j,1:nzm-1)) ) &
                           +(0.5*(sgs_field_sumM(i,j,2:nzm,3)+sgs_field_sumM(i,jb,2:nzm),3) - (betap*d(2:nzm) + betam*v(i,j,2:nzm)) * sgs_field_sumM(i,j,2:nzm,1) )
    vwsb(2:nzm) = vwsb(2:nzm) + rhow(2:nzm) * vwsb3(i,j,2:nzm)
    dvdt(i,j,:,na) = dvdt(i,j,:,na) + (d-v(i,j,:))/dtn
  end do
end do

end subroutine diffuse_mom3D
