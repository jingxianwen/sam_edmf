subroutine write_sgsfields3D(nf)

use vars
use sgs
use params

implicit none

integer,intent(inout) :: nf
integer :: k,j,i
character *80 long_name
character *8 name
character *10 units
real(4) tmp(nx,ny,nzm)


nf=nf+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=twsb3(i,j,k)
    end do
   end do
  end do
  name='TLFLUXED'
  long_name='Liquid water static energy (SGS ED)'
  units='K m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

if (doedmf) then
nf=nf+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)= 0.0
    end do
   end do
  end do
  name='TLFLUXMF'
  long_name='Liquid water static energy (SGS MF)'
  units='K m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

nf=nf+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mkwsb3(i,j,k,1)
    end do
   end do
  end do
  name='QTFLUXED'
  long_name='Non-precipitating water flux (SGS ED)'
  units='g/g m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

if (doedmf) then
nf=nf+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=0.0
    end do
   end do
  end do
  name='QTFLUXMF'
  long_name='Non-precipitating water flux (SGS MF)'
  units='g/g m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nf=nf+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=0.0
    end do
   end do
  end do
  name='CFRACMF'
  long_name='Fractional cloud cover by plumes'
  units='1'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

end if


end subroutine

