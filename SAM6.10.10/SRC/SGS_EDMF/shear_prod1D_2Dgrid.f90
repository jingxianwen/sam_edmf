
subroutine shear_prod1D_2Dgrid(def2)
	
use vars
use params,only :  docolumn
implicit none
	
real def2(nx,ny,nzm), tmpu(nzm), tmpv(nzm)
	
real rdx0,rdx,rdx_up,rdx_dn
real rdz,rdzw_up,rdzw_dn
integer i,j,k,ib,ic,kb,kc

rdx0=1./dx 
j=1

if (docolumn) then
  tmpu=0.
  tmpv=0.
else
  tmpu=u0
  tmpv=v0
end if


do k=2,nzm-1  

  kb=k-1
  kc=k+1
  rdz = 1./(dz*adz(k))
  rdzw_up = 1./(dz*adzw(kc))
  rdzw_dn = 1./(dz*adzw(k))
  rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
  rdx_up=rdx0 * sqrt(dx*rdzw_up) 
  rdx_dn=rdx0 * sqrt(dx*rdzw_dn) 

  do i=1,nx
    ib=i-1
    ic=i+1
	 
      def2(i,j,k)= &
    + 0.25 * ( &
            ( (u(ic,j,kc)-tmpu(kc)-u(ic,j, k)+tmpu(k))*rdzw_up)**2+ &
            ( (u(i ,j,kc)-tmpu(kc)-u(i ,j, k)+tmpu(k))*rdzw_up)**2+ &
            ( (u(ic,j,k )-tmpu(k)-u(ic,j,kb)+tmpu(kb))*rdzw_dn)**2+ &
            ( (u(i ,j,k )-tmpu(k)-u(i ,j,kb)+tmpu(kb))*rdzw_dn)**2)+ &
      0.5 *  (&
            ( (v(i,j ,kc)-tmpv(kc)-v(i,j , k)+tmpv(k))*rdzw_up )**2 + &
            ( (v(i,j ,k )-tmpv(k)-v(i,j ,kb)+tmpv(kb))*rdzw_dn )**2 )

  end do
end do ! k


k=1
kc=k+1

rdz = 1./(dz*adz(k))
rdzw_up = 1./(dz*adzw(kc))
rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
rdx_up=rdx0 * sqrt(dx*rdzw_up) 
	 
do i=1,nx
 ib=i-1
 ic=i+1
	 	
      def2(i,j,k)=( (v(i,j ,kc)-tmpv(kc)-v(i,j,k)+tmpv(k))*rdzw_up )**2 &
   + 0.5 * ( &
           ( (u(ic,j,kc)-tmpu(kc)-u(ic,j, k)+tmpu(k))*rdzw_up)**2+ &
           ( (u(i ,j,kc)-tmpu(kc)-u(i ,j, k)+tmpu(k))*rdzw_up)**2 )
end do 
	 
k=nzm
kc=k+1
kb=k-1

rdz = 1./(dz*adz(k))
rdzw_dn = 1./(dz*adzw(k))
rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
rdx_dn=rdx0 * sqrt(dx*rdzw_dn) 


do i=1,nx
 ib=i-1
 ic=i+1

      def2(i,j,k)=  ( (v(i,j ,k )-tmpv(k)-v(i,j ,kb)+tmpv(kb))*rdzw_dn )**2 &
   + 0.5 * ( &
         ( (u(ic,j,k )-tmpu(k)-u(ic,j,kb)+tmpu(kb))*rdzw_dn)**2+ &
         ( (u(i ,j,k )-tmpu(k)-u(i ,j,kb)+tmpu(kb))*rdzw_dn)**2)

end do 
	
end

