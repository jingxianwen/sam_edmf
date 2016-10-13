
rd=287.
cp=1004.

z=fspan(10.,4990.,250)
nzm=dimsizes(z)
p=z
theta=z
qv=z
u=theta
v=theta
u=7.
v=-5.5

zinvbottom  = 840.
theta0      = 289.
qv0         = 9.0e-03
dthetadzfree= 2.e-03
dqvdzpbl    = - 3.7e-07
dqvdzfree   = - 9.4e-07

do k = 0,nzm-1
   if (z(k).le.zinvbottom) then
      theta(k) = theta0
      qv(k) = qv0
   else
      theta(k) = 297.5 + (z(k)-zinvbottom)^(1./3.)
      qv(k) = 1.5e-03
   end if
end do

qv = where(qv.lt.0.0,0.0,qv)
rv=(qv/(1.-qv))*1000.
delete(qv)

p=-999.

day0=0.0
psfc=1017.8


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



