
subroutine shear_prod1D_3Dgrid(def2)
	
use vars
use params, only : docolumn
implicit none
	
real def2(nx,ny,nzm), tmpu(nzm), tmpv(nzm)
	
real rdx0,rdx,rdx_up,rdx_dn
real rdy0,rdy,rdy_up,rdy_dn
real rdz,rdzw_up,rdzw_dn
integer i,j,k,ib,ic,jb,jc,kb,kc

if (docolumn) then
  tmpu=0.
  tmpv=0.
else
  tmpu=u0
  tmpv=v0
end if


rdx0=1./dx 
rdy0=1./dy

do k=2,nzm-1  

 kb=k-1
 kc=k+1
 rdz = 1./(dz*adz(k))
 rdzw_up = 1./(dz*adzw(kc))
 rdzw_dn = 1./(dz*adzw(k))
 rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
 rdy=rdy0 * sqrt(dy*rdz) 
 rdx_up=rdx0 * sqrt(dx*rdzw_up) 
 rdy_up=rdy0 * sqrt(dy*rdzw_up) 
 rdx_dn=rdx0 * sqrt(dx*rdzw_dn) 
 rdy_dn=rdy0 * sqrt(dy*rdzw_dn) 

 do j=1,ny
   jb=j-YES3D
   jc=j+YES3D
   do i=1,nx
     ib=i-1
     ic=i+1
	 
      def2(i,j,k)= &
        + 0.25 * ( &
          ( (u(ic,j,kc)-tmpu(kc)-u(ic,j, k)+tmpu(k))*rdzw_up)**2+ &
          ( (u(i ,j,kc)-tmpu(kc)-u(i ,j, k)+tmpu(k))*rdzw_up)**2+ &
          ( (u(ic,j,k )-tmpu(k)-u(ic,j,kb)+tmpu(kb))*rdzw_dn)**2+ &
          ( (u(i ,j,k )-tmpu(k)-u(i ,j,kb)+tmpu(kb))*rdzw_dn)**2 )
      def2(i,j,k)=def2(i,j,k) &	
        + 0.25 * ( & 
          ( (v(i,jc,kc)-tmpv(kc)-v(i,jc, k)+tmpv(k))*rdzw_up)**2+ &
          ( (v(i,j ,kc)-tmpv(kc)-v(i,j , k)+tmpv(k))*rdzw_up)**2+ &
          ( (v(i,jc,k )-tmpv(k)-v(i,jc,kb)+tmpv(kb))*rdzw_dn)**2+ &
          ( (v(i,j ,k )-tmpv(k)-v(i,j ,kb)+tmpv(kb))*rdzw_dn)**2 ) 
    end do
 end do
end do ! k


k=1
kc=k+1

rdz = 1./(dz*adz(k))
rdzw_up = 1./(dz*adzw(kc))
rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
rdy=rdy0 * sqrt(dy*rdz) 
rdx_up=rdx0 * sqrt(dx*rdzw_up) 
rdy_up=rdy0 * sqrt(dy*rdzw_up) 
	 
do j=1,ny
  jb=j-YES3D
  jc=j+YES3D
  do i=1,nx
     ib=i-1
     ic=i+1
	 	
      def2(i,j,k)= &
	 + 0.5 * ( &
          ( (v(i,jc,kc)-tmpv(kc)-v(i,jc, k)+tmpv(k))*rdzw_up)**2+ &
          ( (v(i,j ,kc)-tmpv(kc)-v(i,j , k)+tmpv(k))*rdzw_up)**2) &
	 + 0.5 * ( &
          ( (u(ic,j,kc)-tmpu(kc)-u(ic,j, k)+tmpu(k))*rdzw_up)**2+ &
          ( (u(i ,j,kc)-tmpu(kc)-u(i ,j, k)+tmpu(k))*rdzw_up)**2)

   end do 
end do
	 
	
k=nzm
kc=k+1
kb=k-1

rdz = 1./(dz*adz(k))
rdzw_dn = 1./(dz*adzw(k))
rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
rdy=rdy0 * sqrt(dy*rdz) 
rdx_dn=rdx0 * sqrt(dx*rdzw_dn) 
rdy_dn=rdy0 * sqrt(dy*rdzw_dn) 

do j=1,ny
  jb=j-1*YES3D
  jc=j+1*YES3D
  do i=1,nx
      ib=i-1
      ic=i+1
      def2(i,j,k)= & 
       + 0.5 * ( &
           ( (v(i,jc,k )-tmpv(k)-v(i,jc,kb)+tmpv(kb))*rdzw_dn)**2+ &
           ( (v(i,j ,k )-tmpv(k)-v(i,j ,kb)+tmpv(kb))*rdzw_dn)**2) &
 	+ 0.5 * ( &
           ( (u(ic,j,k )-tmpu(k)-u(ic,j,kb)+tmpu(kb))*rdzw_dn)**2+ &
           ( (u(i ,j,k )-tmpu(k)-u(i ,j,kb)+tmpu(kb))*rdzw_dn)**2)
  end do 
end do
	
end

