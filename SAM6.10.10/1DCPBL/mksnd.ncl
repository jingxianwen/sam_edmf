
rd=287.
cp=1004.

z=fspan(10.,4990.,250)
nzm=dimsizes(z)
p=z
theta=z
qv=z

zinvbottom  = 1350.
theta0      = 300.
qv0         = 5.e-03
dthetadzfree= 2.e-03
dqvdzpbl    = - 3.7e-07
dqvdzfree   = - 9.4e-07

do k = 0,nzm-1
   if (z(k).lt.zinvbottom) then
      theta(k) = theta0
      qv(k) = qv0+ z(k) * dqvdzpbl
   else
      theta(k) = theta0+ (z(k)-zinvbottom) * dthetadzfree
      qv(k) = qv0+ zinvbottom * dqvdzpbl+ (z(k)-zinvbottom) * dqvdzfree
   end if
end do

qv = where(qv.lt.0.0,0.0,qv)
rv=(qv/(1.-qv))*1000.
delete(qv)

u=theta
v=theta
u=0.
v=0.
p=-999.

day0=0.0
psfc=1000.


 alist = [/z, p, theta, rv, u, v/]

 header = (/"z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]", \
            "0., "+sprinti("%4.0i",nzm)+" "+sprintf("%4.1f",psfc) +" day,levels,pres0"/)

 hlist = [/header/]
 write_table("snd", "w", hlist, "%s")
 write_table("snd", "a", alist, "%8.2f%6.0f%6.2f%5.2f%6.2f%6.2f")

 delete(header) 
 delete(hlist)
 header = (/"1000000., "+sprinti("%4.0i",nzm)+" "+sprintf("%4.1f",psfc) +" day,levels,pres0"/)
 hlist = [/header/]
 write_table("snd", "a", hlist, "%s")
 write_table("snd", "a", alist, "%8.2f%6.0f%6.2f%5.2f%6.2f%6.2f")



