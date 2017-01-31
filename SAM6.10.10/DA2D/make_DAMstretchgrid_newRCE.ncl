load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; five layer structure originally used for cold pool runs
n1=15
n2=35
n3=42
n4=2
n5=2
dz1=50.
dz3=600.
dz5=600.
nz=n1+n2+n3+n4+n5
dz = new(nz,float)
z=new(nz+1,float)
z(0)=0.0
alpha=1.0

dz2  = 0.5*(dz1+dz3)
dz4  = 0.5*(dz3+dz5)
a=(1.+int2flt(n2))/2.
b=(1.+int2flt(n4))/2.

do i =0,nz-1
  k=i+1

  if (k.le.n1) then
    dz(i)=dz1
  end if
  if (k.ge.(n1+n2+1).and.k.le.(n1+n2+n3)) then
    dz(i)=dz3
  end if
  if (k.ge.(n1+n2+n3+1)) then
    dz(i)=dz5
  end if
  ;if ( k.gt.n1 .and. k.le.(n1+n2)  ) then
  ;   dz(i) = dzm + (dzm-dz1) * (2.0*(k-n1-1)/(int2flt(n2)-1.) - 1.  )^3.
  ;end if
  if ( k.gt.n1 .and. k.le.(n1+n2)  ) then
     dz(i) = dz2 + (dz1-dz2)/tanh(2.*alpha) * tanh(2.0*alpha/(1.-a) * (int2flt(k-n1)-a ))
  end if
  if ( k.gt.(n1+n2+n3) .and. k.le.(n1+n2+n3+n4)  ) then
     dz(i) = dz4 + (dz3-dz4)/tanh(2.*alpha) * tanh(2.0*alpha/(1.-b) * (int2flt(k-(n1+n2+n3) )-b ))
  end if
 
    z(i+1) = z(i) + dz(i)
end do


;************************************************
 wks   = gsn_open_wks ("eps","dz_100mrce")                ; open workstation

 res                  = True                     ; plot mods desired
 res@gsnDraw = False
 res@gsnFrame = False
 res@trXMaxF = 620.;max(dz)+100.
 res@trXMinF = 0.;max(dz)+100.
 res@trYMaxF = max(z)
 res@trYMinF = 000.;max(z)+100.
 res@xyLineThicknessF = 3.
 res@xyMarkLineMode = "MarkLines"
 res@xyMarker = 16

 res@tiXAxisString = "!F8!d!F!z (m)"
 res@tiYAxisString = "Height (m)"


 plot=new(2,graphic)
 plot(0)  = gsn_csm_xy (wks,dz,0.5*(z(0:nz-1)+z(1:nz) ) ,res) ; create plot
 res@trXMinF = 1.;max(dz)+100.
 res@trXMaxF = int2flt(nz);max(dz)+100.
 res@tiXAxisString = "Level"
 plot(1)  = gsn_csm_xy (wks,ispan(1,nz,1),0.5*(z(0:nz-1)+z(1:nz) ) ,res) ; create plot
; plot  = gsn_csm_xy (wks,ispan(1,nz,1),0.5*(z(0:nz-1)+z(1:nz) ) ,res) ; create plot

 gsn_panel(wks,plot,(/1,2/),False)

 mout=new((/nz,1/),float)
 mout(:,0) = 0.5*(z(0:nz-1)+z(1:nz))
 opt=True
 opt@fout="grd_da2d"
 write_matrix(mout,"f11.5",opt)

