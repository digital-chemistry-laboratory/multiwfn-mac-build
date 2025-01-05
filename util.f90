module util
implicit real*8 (a-h,o-z)

interface sort
	module procedure sortr8
	module procedure sorti4
end interface
interface invarr
	module procedure invarrr8
	module procedure invarri4
end interface
!!------------------- Root and weight of Hermite polynomial
real*8 Rhm(10,10),Whm(10,10)
data Rhm(1,1) /  0D0                      /
data Rhm(2,1) / -0.70710678118654752440D+00 /
data Rhm(2,2) /  0.70710678118654752440D+00 /
data Rhm(3,1) / -1.22474487139158904910D+00 /
data Rhm(3,2) /  0D0                      /
data Rhm(3,3) /  1.22474487139158904910D+00 /
data Rhm(4,1) / -1.65068012388578455588D+00 /
data Rhm(4,2) / -0.52464762327529031788D+00 /
data Rhm(4,3) /  0.52464762327529031788D+00 /
data Rhm(4,4) /  1.65068012388578455588D+00 /
data Rhm(5,1) / -2.02018287045608563293D+00 /
data Rhm(5,2) / -0.95857246461381850711D+00 /
data Rhm(5,3) /  0D0                      /
data Rhm(5,4) /  0.95857246461381850711D+00 /
data Rhm(5,5) /  2.02018287045608563293D+00 /
data Rhm(6,1) / -2.35060497367449222283D+00 /
data Rhm(6,2) / -1.33584907401369694971D+00 /
data Rhm(6,3) / -0.43607741192761650868D+00 /
data Rhm(6,4) /  0.43607741192761650868D+00 /
data Rhm(6,5) /  1.33584907401369694971D+00 /
data Rhm(6,6) /  2.35060497367449222283D+00 /
data Rhm(7,1) / -2.65196135683523349245D+00 /
data Rhm(7,2) / -1.67355162876747144503D+00 /
data Rhm(7,3) / -0.81628788285896466304D+00 /
data Rhm(7,4) /  0D0                      /
data Rhm(7,5) /  0.81628788285896466304D+00 /
data Rhm(7,6) /  1.67355162876747144503D+00 /
data Rhm(7,7) /  2.65196135683523349245D+00 /
data Rhm(8,1) / -2.93063742025724401922D+00 /
data Rhm(8,2) / -1.98165675669584292585D+00 /
data Rhm(8,3) / -1.15719371244678019472D+00 /
data Rhm(8,4) / -0.38118699020732211685D+00 /
data Rhm(8,5) /  0.38118699020732211685D+00 /
data Rhm(8,6) /  1.15719371244678019472D+00 /
data Rhm(8,7) /  1.98165675669584292585D+00 /
data Rhm(8,8) /  2.93063742025724401922D+00 /
data Rhm(9,1) / -3.19099320178152760723D+00 /
data Rhm(9,2) / -2.26658058453184311180D+00 /
data Rhm(9,3) / -1.46855328921666793167D+00 /
data Rhm(9,4) / -0.72355101875283757332D+00 /
data Rhm(9,5) /  0D0                      /
data Rhm(9,6) /  0.72355101875283757332D+00 /
data Rhm(9,7) /  1.46855328921666793167D+00 /
data Rhm(9,8) /  2.26658058453184311180D+00 /
data Rhm(9,9) /  3.19099320178152760723D+00 /
data Rhm(10,1) /  -3.43615911883773760333D+00 /
data Rhm(10,2) /  -2.53273167423278979641D+00 /
data Rhm(10,3) /  -1.75668364929988177345D+00 /
data Rhm(10,4) /  -1.03661082978951365418D+00 /
data Rhm(10,5) /  -0.34290132722370460879D+00 /
data Rhm(10,6) /   0.34290132722370460879D+00 /
data Rhm(10,7) /   1.03661082978951365418D+00 /
data Rhm(10,8) /   1.75668364929988177345D+00 /
data Rhm(10,9) /   2.53273167423278979641D+00 /
data Rhm(10,10) /  3.43615911883773760333D+00 /
data Whm(1,1) / 1.77245385090551602730D+00 / ! SQRT(PI)
data Whm(2,1) / 8.86226925452758013649D-01 /
data Whm(2,2) / 8.86226925452758013649D-01 /
data Whm(3,1) / 2.95408975150919337883D-01 /
data Whm(3,2) / 1.18163590060367735153D+00 /
data Whm(3,3) / 2.95408975150919337883D-01 /
data Whm(4,1) / 8.13128354472451771430D-02 /
data Whm(4,2) / 8.04914090005512836506D-01 /
data Whm(4,3) / 8.04914090005512836506D-01 /
data Whm(4,4) / 8.13128354472451771430D-02 /
data Whm(5,1) / 1.99532420590459132077D-02 /
data Whm(5,2) / 3.93619323152241159828D-01 /
data Whm(5,3) / 9.45308720482941881226D-01 /
data Whm(5,4) / 3.93619323152241159828D-01 /
data Whm(5,5) / 1.99532420590459132077D-02 /
data Whm(6,1) / 4.53000990550884564086D-03 /
data Whm(6,2) / 1.57067320322856643916D-01 /
data Whm(6,3) / 7.24629595224392524092D-01 /
data Whm(6,4) / 7.24629595224392524092D-01 /
data Whm(6,5) / 1.57067320322856643916D-01 /
data Whm(6,6) / 4.53000990550884564086D-03 /
data Whm(7,1) / 9.71781245099519154149D-04 /
data Whm(7,2) / 5.45155828191270305922D-02 /
data Whm(7,3) / 4.25607252610127800520D-01 /
data Whm(7,4) / 8.10264617556807326765D-01 /
data Whm(7,5) / 4.25607252610127800520D-01 /
data Whm(7,6) / 5.45155828191270305922D-02 /
data Whm(7,7) / 9.71781245099519154149D-04 /
data Whm(8,1) / 1.99604072211367619206D-04 /
data Whm(8,2) / 1.70779830074134754562D-02 /
data Whm(8,3) / 2.07802325814891879543D-01 /
data Whm(8,4) / 6.61147012558241291030D-01 /
data Whm(8,5) / 6.61147012558241291030D-01 /
data Whm(8,6) / 2.07802325814891879543D-01 /
data Whm(8,7) / 1.70779830074134754562D-02 /
data Whm(8,8) / 1.99604072211367619206D-04 /
data Whm(9,1) / 3.96069772632643819046D-05 /
data Whm(9,2) / 4.94362427553694721722D-03 /
data Whm(9,3) / 8.84745273943765732880D-02 /
data Whm(9,4) / 4.32651559002555750200D-01 /
data Whm(9,5) / 7.20235215606050957124D-01 /
data Whm(9,6) / 4.32651559002555750200D-01 /
data Whm(9,7) / 8.84745273943765732880D-02 /
data Whm(9,8) / 4.94362427553694721722D-03 /
data Whm(9,9) / 3.96069772632643819046D-05 /
data Whm(10,1) /  7.64043285523262062916D-06 /
data Whm(10,2) /  1.34364574678123269220D-03 /
data Whm(10,3) /  3.38743944554810631362D-02 /
data Whm(10,4) /  2.40138611082314686417D-01 /
data Whm(10,5) /  6.10862633735325798784D-01 /
data Whm(10,6) /  6.10862633735325798784D-01 /
data Whm(10,7) /  2.40138611082314686417D-01 /
data Whm(10,8) /  3.38743944554810631362D-02 /
data Whm(10,9) /  1.34364574678123269220D-03 /
data Whm(10,10) / 7.64043285523262062916D-06 /

contains
!Content sequences:
!!Geometry operation
!!Array, Vector
!!String process
!!Matrix calculation
!!Misc

!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Geometry operation !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!

!!---------- Get angle (degree) between two vectors in Degree
real*8 function vecang(vec1x,vec1y,vec1z,vec2x,vec2y,vec2z)
real*8 vec1x,vec1y,vec1z,vec2x,vec2y,vec2z
pi=3.141592653589793D0
rnorm1=dsqrt(vec1x**2+vec1y**2+vec1z**2)
rnorm2=dsqrt(vec2x**2+vec2y**2+vec2z**2)
costheta=(vec1x*vec2x+vec1y*vec2y+vec1z*vec2z)/rnorm1/rnorm2
if (costheta>1D0) costheta=1
vecang=acos(costheta)/pi*180
end function


!!---------- Get distance of point 0 to plane (defined by point 1,2,3) 
real*8 function potpledis(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0)
real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz
call pointprjple(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz)
potpledis=dsqrt((x0-prjx)**2+(y0-prjy)**2+(z0-prjz)**2)
end function


!!---------- Project a point (x0,y0,z0) to a plane defined by three points x/y/z-1/2/3, prjx/y/z are XYZ of the resulting point
subroutine pointprjple(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz)
real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,prjx,prjy,prjz,A,B,C,D,t
call pointABCD(x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D)
! (x0-x)/A=(y0-y)/B=(z0-z)/C=t, t is unknown ---> x=x0-t*A y=y0-t*B z=z0-t*C , substitute into Ax+By+Cz+D=0 solve t
t=(D+A*x0+B*y0+C*z0)/(A**2+B**2+C**2)
prjx=x0-t*A
prjy=y0-t*B
prjz=z0-t*C
end subroutine


!!---------- Project a point (x0,y0,z0) to a plane defined by A,B,C,D, prjx/y/z are XYZ of the resulting point
subroutine pointprjpleABCD(A,B,C,D,x0,y0,z0,prjx,prjy,prjz)
real*8 x0,y0,z0,prjx,prjy,prjz,A,B,C,D,t
t=(D+A*x0+B*y0+C*z0)/(A**2+B**2+C**2)
prjx=x0-t*A
prjy=y0-t*B
prjz=z0-t*C
end subroutine


!!---------- Get distance of a point to the plane defined by point A*x+B*y+C*z+D=0
!See https://en.wikipedia.org/wiki/Plane_(geometry)
!imode=0: Return the absolute distance  =1: Consider sign
subroutine pointABCDdis(x,y,z,pleA,pleB,pleC,pleD,dist,imode)
real*8 x,y,z,pleA,pleB,pleC,pleD,dist
integer imode
if (imode==0) dist=abs(pleA*x+pleB*y+pleC*z+pleD)/dsqrt(pleA**2+pleB**2+pleC**2)
if (imode==1) dist=(pleA*x+pleB*y+pleC*z+pleD)/dsqrt(pleA**2+pleB**2+pleC**2)
end subroutine


!!-------- Input three points, get ABCD of Ax+By+Cz+D=0
subroutine pointABCD(x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D)
real*8 v1x,v1y,v1z,v2x,v2y,v2z,x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D
v1x=x2-x1
v1y=y2-y1
v1z=z2-z1
v2x=x3-x1
v2y=y3-y1
v2z=z3-z1
! Solve determinant (Vector multiply) to get the normal vector (A,B,C):
!  i   j   k   //unit vector
! v1x v1y v1z
! v2x v2y v2z
A=v1y*v2z-v1z*v2y
B=-(v1x*v2z-v1z*v2x)
C=v1x*v2y-v1y*v2x
D=A*(-x1)+B*(-y1)+C*(-z1)
end subroutine


!!-------- Input three points 0,1,2, get the vertical projection point of 0 to the line linking 1-2
subroutine pointprjline(x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz)
real*8 x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz,v12x,v12y,v12z,t
v12x=x2-x1
v12y=y2-y1
v12z=z2-z1
!Since prjx=x1+t*v12x prjy=y1+t*v12y prjz=z1+t*v12z
!So v12x*(x0-prjx)+v12y*(y0-prjy)+v12z*(z0-prjz)=0
!v12x*(x0-x1)-v12x*v12x*t + v12y*(y0-y1)-v12y*v12y*t + v12z*(z0-z1)-v12z*v12z*t =0
t=( v12x*(x0-x1)+v12y*(y0-y1)+v12z*(z0-z1) )/(v12x**2+v12y**2+v12z**2)
prjx=x1+t*v12x
prjy=y1+t*v12y
prjz=z1+t*v12z
end subroutine


!!---------- Get distance of point 0 to the line 1-2
real*8 function potlinedis(x0,y0,z0,x1,y1,z1,x2,y2,z2)
real*8 x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz
call pointprjline(x0,y0,z0,x1,y1,z1,x2,y2,z2,prjx,prjy,prjz)
potlinedis=dsqrt((x0-prjx)**2+(y0-prjy)**2+(z0-prjz)**2)
end function


!!--------- Input two points, return their distance in Bohr
!To calculate distance between two atoms, "atomdist" should be used instead
real*8 function xyz2dist(x1,y1,z1,x2,y2,z2)
real*8 x1,y1,z1,x2,y2,z2
xyz2dist=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
end function


!!--------- Input three points, return angle between 1-2 and 2-3 (in degree)
!To calculate angle between three atoms, "atomang" should be used instead
real*8 function xyz2angle(x1,y1,z1,x2,y2,z2,x3,y3,z3)
real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3
real*8 :: pi=3.141592653589793D0
vec1x=x1-x2
vec1y=y1-y2
vec1z=z1-z2
vec2x=x3-x2
vec2y=y3-y2
vec2z=z3-z2
dotprod=vec1x*vec2x+vec1y*vec2y+vec1z*vec2z
rnormv1=dsqrt( vec1x**2+vec1y**2+vec1z**2 )
rnormv2=dsqrt( vec2x**2+vec2y**2+vec2z**2 )
tmp=dotprod/(rnormv1*rnormv2)
if (tmp>1D0) tmp=1D0
if (tmp<-1D0) tmp=-1D0
xyz2angle=acos(tmp)/pi*180
end function


!!--------- Input four points, return dihedral (in degree)
!The value is in line with internal coordinate dihedral sign rule, namely [-180,180]
!To calculate dihedral between four atoms, "atomdih" should be used instead
real*8 function xyz2dih_sign(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,phi,tmp
real*8 :: pi=3.141592653589793D0
v12x=x2-x1
v12y=y2-y1
v12z=z2-z1
v23x=x3-x2
v23y=y3-y2
v23z=z3-z2
v34x=x4-x3
v34y=y4-y3
v34z=z4-z3
call vecprod(v12x,v12y,v12z,v23x,v23y,v23z,p1x,p1y,p1z)
call vecprod(v23x,v23y,v23z,v34x,v34y,v34z,p2x,p2y,p2z)
tmp=(p1x*p2x+p1y*p2y+p1z*p2z)/(sqrt(p1x*p1x+p1y*p1y+p1z*p1z)*sqrt(p2x*p2x+p2y*p2y+p2z*p2z))
if (tmp>1D0) tmp=1D0
if (tmp<-1D0) tmp=-1D0 
phi=acos(tmp)
xyz2dih_sign=phi/pi*180
call vecprod(p1x,p1y,p1z,p2x,p2y,p2z,tx,ty,tz)
tmp=tx*v23x+ty*v23y+tz*v23z
if (tmp<0) xyz2dih_sign=-xyz2dih_sign
end function



!!-------- Input index of an atom (iatm), return the atom closest to it (iclose). PBC is taken into account
subroutine nearest_atom(iatm,iclose)
use defvar
integer iatm,iclose
distmin=1E10
do itmp=1,ncenter
    if (itmp==iatm) cycle
    tmpdist=atomdist(iatm,itmp,1)
    if (tmpdist<distmin) then
        distmin=tmpdist
        iclose=itmp
    end if
end do
end subroutine



!!-------- Return distance between two atoms in Bohr, consider PBC if iallowPBC=1
real*8 function atomdist(iatm,jatm,iallowPBC)
use defvar
integer iatm,jatm,iallowPBC
real*8 xyzA(3),xyzB(3),xyzB2(3)
if (ifPBC==0.or.iallowPBC==0) then
    atomdist=dsqrt((a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2)
else !Consider PBC, finding nearest jatm that connected to iatm
    xyzA(1)=a(iatm)%x
    xyzA(2)=a(iatm)%y
    xyzA(3)=a(iatm)%z
    xyzB(1)=a(jatm)%x
    xyzB(2)=a(jatm)%y
    xyzB(3)=a(jatm)%z
    call nearest_mirror(xyzA,xyzB,xyzB2)
    atomdist=dsqrt(sum((xyzA-xyzB2)**2))
end if
end function
!!-------- Return distance between two atoms in Angstrom
real*8 function atomdistA(iatm,jatm,iallowPBC)
use defvar
integer iatm,jatm,iallowPBC
atomdistA=atomdist(iatm,jatm,iallowPBC)*b2a
end function


!!-------- Return angle (Degree) between three atoms
real*8 function atomang(iatm,jatm,katm,iallowPBC)
use defvar
integer iatm,jatm,katm,iallowPBC
real*8 xyzA(3),xyzB(3),xyzC(3),xyzA2(3),xyzC2(3)
if (ifPBC==0.or.iallowPBC==0) then
    atomang=xyz2angle(a(iatm)%x,a(iatm)%y,a(iatm)%z,a(jatm)%x,a(jatm)%y,a(jatm)%z,a(katm)%x,a(katm)%y,a(katm)%z)
else !Consider PBC, finding nearest iatm and katm that connected to jatm
    xyzA(1)=a(iatm)%x
    xyzA(2)=a(iatm)%y
    xyzA(3)=a(iatm)%z
    xyzB(1)=a(jatm)%x
    xyzB(2)=a(jatm)%y
    xyzB(3)=a(jatm)%z
    xyzC(1)=a(katm)%x
    xyzC(2)=a(katm)%y
    xyzC(3)=a(katm)%z
    call nearest_mirror(xyzB,xyzA,xyzA2)
    call nearest_mirror(xyzB,xyzC,xyzC2)
    atomang=xyz2angle(xyzA2(1),xyzA2(2),xyzA2(3),xyzB(1),xyzB(2),xyzB(3),xyzC2(1),xyzC2(2),xyzC2(3))
end if
end function


!!-------- Return dihedral between four atoms
real*8 function atomdih(iatm,jatm,katm,latm,iallowPBC)
use defvar
integer iatm,jatm,katm,latm,iallowPBC
real*8 xyzA(3),xyzB(3),xyzC(3),xyzD(3),xyzA2(3),xyzC2(3),xyzD2(3)
if (ifPBC==0.or.iallowPBC==0) then
	atomdih=xyz2dih_sign(a(iatm)%x,a(iatm)%y,a(iatm)%z,a(jatm)%x,a(jatm)%y,a(jatm)%z,a(katm)%x,a(katm)%y,a(katm)%z,a(latm)%x,a(latm)%y,a(latm)%z)
else !Consider PBC, finding nearest iatm, katm and latm that connected to jatm
    xyzA(1)=a(iatm)%x
    xyzA(2)=a(iatm)%y
    xyzA(3)=a(iatm)%z
    xyzB(1)=a(jatm)%x
    xyzB(2)=a(jatm)%y
    xyzB(3)=a(jatm)%z
    xyzC(1)=a(katm)%x
    xyzC(2)=a(katm)%y
    xyzC(3)=a(katm)%z
    xyzD(1)=a(latm)%x
    xyzD(2)=a(latm)%y
    xyzD(3)=a(latm)%z
    call nearest_mirror(xyzB,xyzA,xyzA2)
    call nearest_mirror(xyzB,xyzC,xyzC2)
    call nearest_mirror(xyzB,xyzD,xyzD2)
    atomdih=xyz2dih_sign(xyzA2(1),xyzA2(2),xyzA2(3),xyzB(1),xyzB(2),xyzB(3),xyzC2(1),xyzC2(2),xyzC2(3),xyzD2(1),xyzD2(2),xyzD2(3))
end if
end function


!!--------------- Get area of a triangle, need input coordinates of three points
function gettriangarea(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz)
implicit real*8 (a-h,o-z)
real*8 gettriangarea,pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz
! a---b             va=pb-pa vb=pc-pa
! |
! V
! c
va1=pbx-pax
va2=pby-pay
va3=pbz-paz
vb1=pcx-pax
vb2=pcy-pay
vb3=pcz-paz
call vecprod(va1,va2,va3,vb1,vb2,vb3,vc1,vc2,vc3)  !vc=va×vb=|va||vb|sinθ*i  where i is unit vector perpendicular to va and vb
absvc=dsqrt(vc1**2+vc2**2+vc3**2)
gettriangarea=0.5D0*absvc
end function


!!--------------- Get volume of a tetrahedron, need input coordinates of four points
function gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz)
implicit real*8 (a-h,o-z)
real*8 pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz
! real*8 volmat(4,4)
! volmat(:,1)=(/ pax,pbx,pcx,pdx /)
! volmat(:,2)=(/ pay,pby,pcy,pdy /)
! volmat(:,3)=(/ paz,pbz,pcz,pdz /)
! volmat(:,4)=1D0
! gettetravol=abs(detmat(volmat))/6D0
! call showmatgau(volmat)
!vol=abs( (a-d)·((b-d)×(c-d)) )/6,  see http://en.wikipedia.org/wiki/Tetrahedron
vec1x=pax-pdx
vec1y=pay-pdy
vec1z=paz-pdz
vec2x=pbx-pdx
vec2y=pby-pdy
vec2z=pbz-pdz
vec3x=pcx-pdx
vec3y=pcy-pdy
vec3z=pcz-pdz
call vecprod(vec2x,vec2y,vec2z,vec3x,vec3y,vec3z,vec2x3x,vec2x3y,vec2x3z)
gettetravol=abs(vec1x*vec2x3x+vec1y*vec2x3y+vec1z*vec2x3z)/6D0
end function


!!----- Generate rotation matrix using the algorithm "Efficiently Building a Matrix to Rotate One Vector to Another (1999)"
!vec1 will point towards to vec2, return rotation matrix "mat"
subroutine rotmat_vec1_vec2(vec1,vec2,mat)
implicit real*8 (a-h,o-z)
real*8 vec1(3),vec2(3),mat(3,3)
real*8 vecv(3),vecu(3),vecc(3),vecx(3)

test=abs(sum(vec1*vec2))

if (test<0.99D0) then !The two inputted vectors are not nearly orthogonal with each other
	call vecprod(vec1(1),vec1(2),vec1(3),vec2(1),vec2(2),vec2(3),vecv(1),vecv(2),vecv(3))
	vecu=vecv/dsqrt(sum(vecv**2))
	valc=sum(vec1*vec2)
	valh=(1-valc)/(1-valc**2)
	vx=vecv(1)
	vy=vecv(2)
	vz=vecv(3)
	mat(1,1)=valc+valh*vx**2
	mat(1,2)=valh*vx*vy-vz
	mat(1,3)=valh*vx*vz+vy
	mat(2,1)=valh*vx*vy+vz
	mat(2,2)=valc+valh*vy**2
	mat(2,3)=valh*vy*vz-vx
	mat(3,1)=valh*vx*vz-vy
	mat(3,2)=valh*vy*vz+vx
	mat(3,3)=valc+valh*vz**2
    
else
	!Find vecx, which is the Cartesian axis most nearly orthgonal with vec1
	ixdir=1
    absmin=abs(vec1(1))
	do i=2,3
		if (abs(vec1(i))<absmin) then
			ixdir=i
            absmin=abs(vec1(i))
        end if
    end do
    vecx=0
    if (ixdir==1) vecx(1)=1
    if (ixdir==2) vecx(2)=1
    if (ixdir==3) vecx(3)=1
	vecu=vecx-vec1
	vecv=vecx-vec2
    vecusqr=sum(vecu**2)
    vecvsqr=sum(vecv**2)
    do i=1,3
		do j=1,3
			deltmp=0
            if (i==j) deltmp=1
			mat(i,j)= deltmp - 2*vecu(i)*vecu(j)/vecusqr - 2*vecv(i)*vecv(j)/vecvsqr + 4*sum(vecu*vecv)/vecusqr/vecvsqr*vecv(i)*vecu(j)
        end do
    end do
end if
end subroutine


!!---------- Get rotation matrix by rotating specific angle around Z,Y,X axis
!See "General rotations" of https://en.wikipedia.org/wiki/Rotation_matrix
!a,b,c are rotation angle around Z,Y,X axis in degree, mat is rotation matrix
subroutine get_rotmat_aroundXYZ(ain,bin,cin,mat)
implicit real*8 (a-h,o-z)
real*8 a,b,c,ain,bin,cin,mat(3,3)
real*8 :: pi=3.141592653589793D0
a=ain/180*pi
b=bin/180*pi
c=cin/180*pi
sina=sin(a)
sinb=sin(b)
sinc=sin(c)
cosa=cos(a)
cosb=cos(b)
cosc=cos(c)
mat(1,1)=cosa*cosb
mat(1,2)=cosa*sinb*sinc-sina*cosc
mat(1,3)=cosa*sinb*cosc+sina*sinc
mat(2,1)=sina*cosb
mat(2,2)=sina*sinb*sinc+cosa*cosc
mat(2,3)=sina*sinb*cosc-cosa*sinc
mat(3,1)=-sinb
mat(3,2)=cosb*sinc
mat(3,3)=cosb*cosc
end subroutine


!!---------- Get rotation matrix by rotating specific angle around a vector
!See "Rotation matrix from axis and angle" of https://en.wikipedia.org/wiki/Rotation_matrix
!rotdeg is rotation angle around the vector (vec) in degree, mat is rotation matrix
subroutine get_rotmat_aroundvec(vec,rotdeg,mat)
implicit real*8 (a-h,o-z)
real*8 vec(3),rotdeg,mat(3,3)
real*8 :: pi=3.141592653589793D0
cosv=cos(rotdeg/180*pi)
sinv=sin(rotdeg/180*pi)
vec=vec/dsqrt(sum(vec**2)) !Normalize
vx=vec(1)
vy=vec(2)
vz=vec(3)
mat(1,1)=cosv+vx**2*(1-cosv)
mat(1,2)=vx*vy*(1-cosv)-vz*sinv
mat(1,3)=vx*vz*(1-cosv)+vy*sinv
mat(2,1)=vy*vx*(1-cosv)+vz*sinv
mat(2,2)=cosv+vy**2*(1-cosv)
mat(2,3)=vy*vz*(1-cosv)-vx*sinv
mat(3,1)=vz*vx*(1-cosv)-vy*sinv
mat(3,2)=vz*vy*(1-cosv)-vx*sinv
mat(3,3)=cosv+vz**2*(1-cosv)
end subroutine


!!---------- Input a set of atoms to fit a best plane. Return A, B, C, D of plane equation A*x+B*y+C*z+D=0
!"atmarr" with size of "natm" is the array containing atom indices in the ring. rmsfit measures RMS fitting error in Bohr
!  Based on joriki's answer: https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
!  "Subtract out the centroid, form a 3×N matrix X out of the resulting coordinates and calculate its singular value decomposition. &
!The normal vector of the best-fitting plane is the left singular vector corresponding to the least singular value"
subroutine ptsfitplane(atmarr,natm,planeA,planeB,planeC,planeD,rmsfit)
use defvar
implicit real*8 (a-h,o-z)
integer natm,atmarr(natm)
real*8 planeA,planeB,planeC,planeD,rmsfit
real*8 tmpmat(3,natm),singval(3),matU(3,3),matV(natm,natm)

cenx=sum(a(atmarr(:))%x)/natm
ceny=sum(a(atmarr(:))%y)/natm
cenz=sum(a(atmarr(:))%z)/natm
tmpmat(1,:)=a(atmarr(:))%x-cenx
tmpmat(2,:)=a(atmarr(:))%y-ceny
tmpmat(3,:)=a(atmarr(:))%z-cenz

call SVDmat(1,tmpmat,matU,matV,singval,info)
ileast=minloc(singval,1)

planeA=matU(1,ileast)
planeB=matU(2,ileast)
planeC=matU(3,ileast)
planeD=-(planeA*cenx+planeB*ceny+planeC*cenz)

!Check fitting quality using A*x+B*y+C*z+D=0
accum=0
do iatm=1,natm
    tmp=planeA*a(atmarr(iatm))%x + planeB*a(atmarr(iatm))%y + planeC*a(atmarr(iatm))%z + planeD 
    accum=accum+abs(tmp)
end do
rmsfit=dsqrt(accum/natm)
end subroutine


!!--------- Input an angle in radian, return degree
real*8 function rad2ang(rad)
use defvar
real*8 rad
rad2ang=rad/pi*180D0
end function


!!--------- Input an angle in degree, return radian
real*8 function ang2rad(ang)
use defvar
real*8 ang
ang2rad=ang/180D0*pi
end function


!!--------- Get projection of a 3*3 matrix along a given vector. For example, component of magnetic shielding tensor
!The vector will be automatically normalized here
real*8 function prjmat(mat,vec)
real*8 mat(3,3),vec(3),uvec(3)
uvec=vec/dsqrt(sum(vec**2))
prjmat=0
do idir=1,3
    do jdir=1,3
        prjmat=prjmat+uvec(idir)*uvec(jdir)*mat(idir,jdir)
    end do
end do
end function




!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! Array, Vector !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!

!!--------- Invert array, from istart to iend
!If "list" is presented, also exchange the index during invert. list should have same size as array
!Real*8 version
subroutine invarrr8(array,list,list2)
real*8 array(:)
integer,optional :: list(:),list2(:)
len=size(array)
ilist=0
if (present(list)) ilist=1
ilist2=0
if (present(list2)) ilist2=1
do i=1,int(len/2D0)
	j=len-i+1
	tmp=array(i)
	array(i)=array(j)
	array(j)=tmp
	if (ilist==1) then
		itemp=list(i)
		list(i)=list(j)
		list(j)=itemp
	end if
	if (ilist2==1) then
		itemp=list2(i)
		list2(i)=list2(j)
		list2(j)=itemp
	end if
end do
end subroutine
!Integer 4 version
subroutine invarri4(array,list,list2)
integer array(:)
integer,optional :: list(:),list2(:)
len=size(array)
ilist=0
if (present(list)) ilist=1
ilist2=0
if (present(list2)) ilist2=1
do i=1,int(len/2D0)
	j=len-i+1
	itmp=array(i)
	array(i)=array(j)
	array(j)=itmp
	if (ilist==1) then
		itemp=list(i)
		list(i)=list(j)
		list(j)=itemp
	end if
	if (ilist2==1) then
		itemp=list2(i)
		list2(i)=list2(j)
		list2(j)=itemp
	end if
end do
end subroutine


!-------- Sort value from small to big by Bubble method
!inmode =abs: sort by absolute value, =val: sort by value. Default is by value
!If "list" is presented, also exchange the index in the list during sorting. list should have same size as array
!If want to sort from big to small, then you should use invarrr8 or invarri4 to invert the array
!BTW: To sort in descending order, one modify the code as
!do i=1,N-1
!	do j=1,N-i
!		if (array(j)<array(j+1)) then
!		...
!       end if
!   end do
!end do
!
!Real*8 version
subroutine sortr8(array,inmode,list,list2)
real*8 array(:)
character,optional :: inmode*3
integer,optional :: list(:),list2(:)
N=size(array)
mode=1
if (present(inmode)) then
	if (inmode=="abs") mode=2
end if
ilist=0
if (present(list)) ilist=1
ilist2=0
if (present(list2)) ilist2=1
if (mode==1) then
	do i=1,N-1
		do j=i+1,N
			if (array(i)>array(j)) then
				temp=array(i)
				array(i)=array(j)
				array(j)=temp
				if (ilist==1) then
					itemp=list(i)
					list(i)=list(j)
					list(j)=itemp
				end if
				if (ilist2==1) then
					itemp=list2(i)
					list2(i)=list2(j)
					list2(j)=itemp
				end if
			end if
		end do
	end do
else if (mode==2) then
	do i=1,N-1
		do j=i+1,N
			if (abs(array(i))>abs(array(j))) then
				temp=array(i)
				array(i)=array(j)
				array(j)=temp
				if (ilist==1) then
					itemp=list(i)
					list(i)=list(j)
					list(j)=itemp
				end if
				if (ilist2==1) then
					itemp=list2(i)
					list2(i)=list2(j)
					list2(j)=itemp
				end if
			end if
		end do
	end do
end if
end subroutine
!Integer 4 version
subroutine sorti4(array,inmode,list,list2)
integer array(:)
character,optional :: inmode*3
integer,optional :: list(:),list2(:)
N=size(array)
mode=1
if (present(inmode)) then
	if (inmode=="abs") mode=2
end if
ilist=0
if (present(list)) ilist=1
ilist2=0
if (present(list2)) ilist2=1
if (mode==1) then
	do i=1,N-1
		do j=i+1,N
			if (array(i)>array(j)) then
				itemp=array(i)
				array(i)=array(j)
				array(j)=itemp
				if (ilist==1) then
					itemp=list(i)
					list(i)=list(j)
					list(j)=itemp
				end if
				if (ilist2==1) then
					itemp=list2(i)
					list2(i)=list2(j)
					list2(j)=itemp
				end if
			end if
		end do
	end do
else if (mode==2) then
	do i=1,N-1
		do j=i+1,N
			if (abs(array(i))>abs(array(j))) then
				itemp=array(i)
				array(i)=array(j)
				array(j)=itemp
				if (ilist==1) then
					itemp=list(i)
					list(i)=list(j)
					list(j)=itemp
				end if
				if (ilist2==1) then
					itemp=list2(i)
					list2(i)=list2(j)
					list2(j)=itemp
				end if
			end if
		end do
	end do
end if
end subroutine


!-------- Sort an index array, namely array(1:nslot,1:nmember). For example, array(:,1)=1,3,9  array(:,2)=2,6,4  etc.
!Sort from left to right (slot 1 to slot "nslot"); when processing slot i, arrays are reordered so that the index of this slot is from small to large, but order of slot i-1 remains unchanged
!For example, original list is
!7,5,6
!3,8,2
!4,6,9
!7,5,1
!4,2,3
!After sorting the list will be
!3,8,2
!4,2,3
!4,6,9
!7,5,1
!7,5,6
subroutine sortidxlist(array,nslot,nmember)
integer nslot,nmember
integer array(nslot,nmember),arrtmp(nslot)
do islot=1,nslot
	do imember=1,nmember-1
cycjmem:do jmember=imember+1,nmember-1
			if (array(islot,imember)>array(islot,jmember)) then
                !If changing order of islot will reverse order of islot-1, then do not change order
                if (islot>=2) then
                    do itest=1,islot-1
					    if (array(itest,imember)<array(itest,jmember)) cycle cycjmem
                    end do
                end if
				arrtmp(:)=array(:,imember)
				array(:,imember)=array(:,jmember)
				array(:,jmember)=arrtmp(:)
            end if
        end do cycjmem
    end do
end do
end subroutine


!!--------- Evaluate standard deviation of array elements
real*8 function stddevarray(array)
real*8 array(:),avg
avg=sum(array)/size(array)
stddevarray=dsqrt(sum((array-avg)**2)/size(array))
end function


!!--------- Evaluate covariant of two array elements
real*8 function covarray(array1,array2)
real*8 array1(:),array2(:),avg1,avg2
avg1=sum(array1)/size(array1)
avg2=sum(array2)/size(array2)
covarray=sum((array1-avg1)*(array2-avg2))/size(array1)
end function


!--- Cross product (vector product) of two vectors, return a new vector (x,y,z)
subroutine vecprod(x1,y1,z1,x2,y2,z2,x,y,z)
real*8 x1,y1,z1,x2,y2,z2,x,y,z
! |i  j  k |
! |x1 y1 z1|
! |x2 y2 z2|
x=  y1*z2-z1*y2
y=-(x1*z2-z1*x2)
z=  x1*y2-y1*x2
end subroutine


!--- The same as vecprod, but input and output are vectors
subroutine vecprodv(vec1,vec2,newvec)
real*8 x1,y1,z1,x2,y2,z2
real*8 vec1(3),vec2(3),newvec(3)
x1=vec1(1)
y1=vec1(2)
z1=vec1(3)
x2=vec2(1)
y2=vec2(2)
z2=vec2(3)
newvec(1)=  y1*z2-z1*y2
newvec(2)=-(x1*z2-z1*x2)
newvec(3)=  x1*y2-y1*x2
end subroutine


!--- Cross product of two vectors, return a new vector (x,y,z)
!Unlike vecprodv, this take care of special case when the two vectors are collinear with each other
subroutine crossprod(vec1,vec2,newvec)
implicit real*8 (a-h,o-z)
real*8 vec1(3),vec2(3),newvec(3),vect(3)
call test_vec_parallel(vec1,vec2,istat)
if (istat/=0) then !Collinear
	vect=(/ 1D0,-1D0,1D0 /) !An arbitrary vector
    call test_vec_parallel(vec1,vect,istat)
    if (istat/=0) then !Collinear
		vect=(/ -1D0,1D0,1D0 /) !Another arbitrary vector
		call vecprodv(vec1,vect,newvec)
	else !Not collinear
		call vecprodv(vec1,vect,newvec)
	end if
else !Not collinear
	call vecprodv(vec1,vec2,newvec)
end if
end subroutine


!-------- Test if two vectors are parallel or antiparallel
!istat=0: Not parallel or antiparallel
!istat=1: Nearly parallel
!istat=2: Nearly antiparallel
!ref.: https://stackoverflow.com/questions/7572640/how-do-i-know-if-two-vectors-are-near-parallel
subroutine test_vec_parallel(vec1,vec2,istat)
real*8 vec1(3),vec2(3),norm1,norm2,crit
integer istat
norm1=dsqrt(sum(vec1**2))
norm2=dsqrt(sum(vec2**2))
cosang=sum(vec1*vec2)/norm1/norm2
crit=1D-5
if (cosang>1D0-crit) then
	istat=1
else if (cosang<-1D0+crit) then
	istat=2
else
	istat=0
end if
end subroutine


!--- Calculate area defined by two vectors, namely |v1 x v2|
subroutine vec2area(vec1,vec2,area)
real*8 vec1(3),vec2(3),vecnew(3),area
call vecprod(vec1(1),vec1(2),vec1(3),vec2(1),vec2(2),vec2(3),vecnew(1),vecnew(2),vecnew(3))
area=dsqrt(vecnew(1)**2+vecnew(2)**2+vecnew(3)**2)
end subroutine


!--- Generate full arrangement (all permutation) array
!ncol is the number of elements, nrow should have size of (ncol)!
!Example: nrow=3!=6, ncol=3, the returned arr will be:
! 1           3           2
! 2           1           3
! 2           3           1
! 3           1           2
! 3           2           1
! 1           2           3
subroutine fullarrange(arr,nrow,ncol)
integer nrow,ncol,arr(nrow,ncol),seq(ncol)
seq=(/ (i,i=1,ncol) /)
arr(1,:)=seq !The first array will be 1,2,3,4...
do icyc=2,nrow
	do i=ncol-1,1,-1
		if (seq(i)<seq(i+1)) exit
	end do
	do j=ncol,1,-1
		if (seq(j)>seq(i)) exit
	end do
	itmp=seq(i)
	seq(i)=seq(j)
	seq(j)=itmp
	call invarr(seq(i+1:ncol))
	arr(icyc,:)=seq
end do
end subroutine







!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! String process !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!

!-------- Parse inputted integer string (e.g. 3,4,5,6,7-25,32-35,55,88) to array
!nelement will be return, which is the number of indices in the string
!if "array" is present, the terms will be written into it, else if "array" is not present, only count the number of terms as "nelement"
!array must be large enough to contain the elements
subroutine str2arr(inpstr,nelement,array)
character(len=*) inpstr
character c80tmp*80
integer nelement
integer,optional :: array(:)
if (present(array)) array=0
icommaold=0
nelement=0
!Loop the string to find comma
do ipos=1,len_trim(inpstr)
	!Deal with the string between current comma(or terminal) and last comma(or beginning)
	if (inpstr(ipos:ipos)==','.or.ipos==len_trim(inpstr)) then
		icommanew=ipos
		if (ipos==len_trim(inpstr)) icommanew=ipos+1
		read(inpstr(icommaold+1:icommanew-1),*) c80tmp !The sub-string to be parsed
		if (index(c80tmp,'-')==0) then !The sub-string corresponds to a single number
			nelement=nelement+1
			if (present(array)) read(c80tmp,*) array(nelement)
		else
			iheng=index(c80tmp,'-')
			read(c80tmp(1:iheng-1),*) ilow
			read(c80tmp(iheng+1:len_trim(c80tmp)),*) ihigh
			if (present(array)) then
				do itmp=ilow,ihigh
					nelement=nelement+1
					array(nelement)=itmp
				end do
			else
				nelement=nelement+ihigh-ilow+1
			end if
		end if
		icommaold=icommanew
	end if
end do
end subroutine


!-------- Convert an array containing integers to a compact string (str), no "-" will be used
!"str" must have large enough characters
subroutine arr2str(array,str)
character(len=*) str
character c80tmp*80
integer :: array(:)
nele=size(array)
ipos=1
str=" "
do iele=1,nele
    write(c80tmp,*) array(iele)
    nlen=len_trim(adjustl(c80tmp))
    if (iele<nele) then
        str(ipos:ipos+nlen)=trim(adjustl(c80tmp))//','
        ipos=ipos+nlen+1
    else
        str(ipos:ipos+nlen-1)=trim(adjustl(c80tmp))
    end if
end do
end subroutine

!!-------- The same as arr2str, but may contain "-" symbol
!All value in the array must be >=0, negative is not supported, after all, the negative sign must conflict with - symbol
!Algorithm: Assume the maximal index is 8, then arrtmp(0:8) is allocated; if an index is in the array, the corresponding slot will be 1, else 0
!  we scan from 0 to 8, 
subroutine arr2str_2(array,str)
character(len=*) str
character c80tmp*80,c2000tmp*2000
integer :: array(:)
integer,allocatable :: arrtmp(:)
ntmp=maxval(array)
allocate(arrtmp(0:ntmp))
arrtmp=0
do itmp=1,size(array)
    arrtmp(array(itmp))=1
end do
c2000tmp=" "
do itmp=0,ntmp
    if (arrtmp(itmp)==1) then
        write(c80tmp,*) itmp
        if (c2000tmp/=" ") then !String has already recorded at least one index
            if (arrtmp(itmp-1)==1) then !The current value is contiguous to the last one, may use - to connect to last one
                if (itmp<ntmp) then !This is not the last index in the inputted array
                    if (arrtmp(itmp+1)==1) cycle !If the next value also occurs in the inputted array, the current value will be outputted to the string
                end if
                ifmod=1
                if (itmp-2>=0) then
					if (arrtmp(itmp-2)==0) ifmod=0
                end if
                if (ifmod==1) then
					ntmp2=len_trim(c2000tmp) !Modify last comma as -
					c2000tmp(ntmp2:ntmp2)='-'
                end if
            end if
        end if
        c2000tmp=trim(c2000tmp)//trim(adjustl(c80tmp))
        if (any(arrtmp(itmp+1:)==1)) write(c2000tmp(len_trim(c2000tmp)+1:),"(a)") ','
    end if
end do
str=c2000tmp
end subroutine


!---------Input path name, e.g. C:\ltwd\MIO.wfn, output file name, e.g. MIO
subroutine path2filename(pathnamein,filenameout)
character(len=*) pathnamein,filenameout
do i=len_trim(pathnamein),1,-1
	if (pathnamein(i:i)=='/'.or.pathnamein(i:i)=='\') exit 
end do
filenameout=pathnamein(i+1:)
i=index(filenameout,'.')
if (i/=0) filenameout(i:)=" "
end subroutine


!-------- Add name of input file as prefix of given string
subroutine addprefix(inname)
use defvar
character c200tmp*200
character(len=*) inname
call path2filename(filename,c200tmp)
inname=trim(c200tmp)//'_'//trim(inname)
end subroutine


!!--------- Convert a character to lower case
subroutine uc2lc(inc)
character inc
itmp=ichar(inc)
if (itmp>=65.and.itmp<=90) itmp=itmp+32
inc=char(itmp)
end subroutine


!!--------- Convert a character to upper case
subroutine lc2uc(inc)
character inc
itmp=ichar(inc)
if (itmp>=97.and.itmp<=122) itmp=itmp-32
inc=char(itmp)
end subroutine


!!--------- Convert a string to lower case
subroutine struc2lc(str)
character(len=*) str
do i=1,len_trim(str)
	call uc2lc(str(i:i))
end do
end subroutine


!!--------- Convert a string to upper case
subroutine strlc2uc(str)
character(len=*) str
do i=1,len_trim(str)
	call lc2uc(str(i:i))
end do
end subroutine


!!--------- Return the position of a given character in a string, itime is the number of times that it occurs
!e.g. strcharpos(sfisi1123,'i',1) returns 3, strcharpos(sfisi1123,'i',2) returns 5. If not found, return 0
function strcharpos(str,char,itime)
character(len=*) str
character char
integer itime
strcharpos=0
icount=0
do i=1,len_trim(str)
    if (str(i:i)==char) icount=icount+1
    if (icount==itime) then
        strcharpos=i
        return
    end if
end do
end function


!!--------- Return the number of a character in a string
!e.g. strcharnum(sfisi1123,'i') returns 2
function strcharnum(str,char)
character(len=*) str
character char
integer itime
strcharnum=0
do i=1,len_trim(str)
    if (str(i:i)==char) strcharnum=strcharnum+1
end do
end function


!!-------- Read float data after the last specific sign (can be multiple characters) from a string
subroutine readaftersign_float(str,sign,val)
character(len=*) str,sign
real*8 val
itmp=index(trim(str),sign,back=.true.)
read(str(itmp+len(sign):),*) val
end subroutine
!!-------- Read integer data after the last specific sign from a string
subroutine readaftersign_int(str,sign,val)
character(len=*) str,sign
integer val
itmp=index(trim(str),sign,back=.true.)
read(str(itmp+len(sign):),*) val
end subroutine

!!-------- Read float data from "ifileid" file
!Locate to the line containing "label" first, then read data after "sign" to "val". If the label was not found, "val" keeps original value
subroutine readoption_float(ifileid,label,sign,val)
integer ifileid
character(len=*) label,sign
character str*200
real*8 val
call loclabel(ifileid,label,ifound)
if (ifound==1) then
    read(ifileid,"(a)") str
    itmp=index(trim(str),sign,back=.true.)
    read(str(itmp+len(sign):),*) val
end if
end subroutine
!!----- Same as above, but for loading integer value
subroutine readoption_int(ifileid,label,sign,val)
integer ifileid
character(len=*) label,sign
character str*200
integer val
call loclabel(ifileid,label,ifound)
if (ifound==1) then
    read(ifileid,"(a)") str
    itmp=index(trim(str),sign,back=.true.)
    read(str(itmp+len(sign):),*) val
end if
end subroutine
!!----- Same as above, but for loading a float vector
subroutine readoption_vec_float(ifileid,label,sign,vec)
integer ifileid
character(len=*) label,sign
character str*200
real*8 vec(:)
call loclabel(ifileid,label,ifound)
if (ifound==1) then
    read(ifileid,"(a)") str
    itmp=index(trim(str),sign) !Load information after the first given "sign"
    read(str(itmp+len(sign):),*) vec
end if
end subroutine
!!----- Same as above, but for loading a integer vector
subroutine readoption_vec_int(ifileid,label,sign,vec)
integer ifileid
character(len=*) label,sign
character str*200
integer vec(:)
call loclabel(ifileid,label,ifound)
if (ifound==1) then
    read(ifileid,"(a)") str
    itmp=index(trim(str),sign) !Load information after the first given "sign"
    read(str(itmp+len(sign):),*) vec
end if
end subroutine


!!-------- Locate to the line containing given "label" in "ifileid" file, return the string ready for loading data via read(str,*)
!For example, call get_option_str(ifileid,"graph2Dsize=",str), will locate to the following line
!   graph2Dsize= 1500,1200 // Width and height...
!Then "str" will be 1500,1200
!If str=" ", that means the line was not found
subroutine get_option_str(ifileid,label,str)
integer ifileid
character(len=*) label,str
character c200tmp*200
str=" "
call loclabel(ifileid,label,ifound)
if (ifound==1) then
	read(ifileid,"(a)") c200tmp
	ibeg=index(c200tmp,'=')
    iend=index(c200tmp,'//')
    if (iend==0) iend=len_trim(c200tmp)+1
    str=adjustl(c200tmp(ibeg+1:iend-1))
end if
end subroutine


!!------- Determine how many data in the given string. Maximally test 1000
subroutine numdatastr(str,ndata)
character(len=*) str
character(len=100) arr(1000)
integer ndata
do i=1,1000
	read(str,*,iostat=ierror) arr(1:i)
	if (ierror/=0) then
        ndata=i-1
        return
    end if
end do
end subroutine


!!------ Remove all parentheses and the data within them from a string, namely setting them as space
!For example, "lutian (saint) is weida?" will be "lutian         is weida?" 
subroutine remove_parentheses(str)
character(len=*) str
do while(.true.)
    if (index(str,'(')==0) exit
    itmp=index(str,'(')
    jtmp=index(str,')')
    str(itmp:jtmp)=" "
end do
end subroutine


!!------- Remove the first and the last " or ' symbol, because directly dragging file into the window will result in " or ' symbol, which is unrecognized
subroutine remove_quotemark(str)
character(len=*) str
ltmp=len_trim(str)
if (str(1:1)=='"'.or.str(1:1)=="'") str(1:1)=" "
if (str(ltmp:ltmp)=='"'.or.str(ltmp:ltmp)=="'") str(ltmp:ltmp)=" "
str=adjustl(str)
end subroutine


!!------- Remove replace a character from a string with " " if any
subroutine remove_char(str,char)
character(len=*) str
character char
itmp=index(str,char)
if (itmp/=0) str(itmp:itmp)=" "
end subroutine


!!------- Show content of "nline" lines of fileid on screen, and the position of reading is unchanged
subroutine showlines(ifileid,nline)
character c200tmp*200
do iline=1,nline
	read(ifileid,"(a)") c200tmp
    write(*,"(a)") trim(c200tmp)
end do
do iline=1,nline
	backspace(ifileid)
end do
end subroutine





!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Matrix calculation !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!

!!----- Make the matrix to upper trigonal matrix
subroutine ratio_upper(mat)
real*8 :: mat(:,:),m,st
real*8,allocatable :: temp(:),s(:),divided(:)
integer a,i,j,t
a=size(mat,1)
allocate(temp(a))
allocate(s(a))
allocate(divided(a))
do i=1,a
	s(i)=maxval(abs(mat(i,1:a)))
end do
do i=1,a-1
	divided(i:a)=mat(i:a,i)/s(i:a)
	t=maxloc(abs(divided(i:a)),dim=1)
	temp(:)=mat(i,:)
	mat(i,:)=mat(i+t-1,:)
	mat(i+t-1,:)=temp(:)
	st=s(i)
	s(i)=s(i+t-1)
	s(i+t-1)=st
	do j=i+1,a
		m=mat(j,i)/mat(i,i)
		mat(j,i:a)=mat(j,i:a)-mat(i,i:a)*m
	end do
end do
deallocate(temp,s,divided)
end subroutine


!!----- Get value of determinant of a matrix
real*8 function detmat(mat)
real*8 mat(:,:)
real*8,allocatable :: mattmp(:,:)
isizemat=size(mat,1)
detmat=1D0
NOTlowertri=0
NOTuppertri=0

outter1: do i=1,isizemat !Check if already is lower-trigonal matrix
	do j=i+1,isizemat
		if (abs(mat(i,j))>1D-12) then
			NOTlowertri=1 !There are at least one big value at upper trigonal part, hence not lower trigonal matrix
			exit outter1
		end if
	end do
end do outter1
outter2: do i=1,isizemat !Check if already is upper-trigonal matrix
	do j=1,i-1
		if (abs(mat(i,j))>1D-12) then
			NOTuppertri=1 !There are at least one big value at lower trigonal part, hence not upper trigonal matrix
			exit outter2
		end if
	end do
end do outter2

if (NOTlowertri==0.or.NOTuppertri==0) then !Is lower or upper trigonal matrix, don't need to convert to trigonal matrix
	do i=1,isizemat
		detmat=detmat*mat(i,i)
	end do
else !Not upper or lower trigonal matrix
	allocate(mattmp(isizemat,isizemat))
	mattmp=mat
	call ratio_upper(mattmp)
	detmat=1D0
	do i=1,isizemat
		detmat=detmat*mattmp(i,i)
	end do
end if
end function



!!-------- Get trace of a matrix
real*8 function mattrace(mat)
real*8 mat(:,:)
mattrace=0
do i=1,size(mat,1)
	mattrace=mattrace+mat(i,i)
end do
end function


!!--- Use Jacobi method to diagonalize matrix, simple, but much slower than diagsymat and diaggemat if the matrix is large
subroutine diagmat(mat,S,eigval,inmaxcyc,inthres)
! mat: input and will be diagonalized matrix, S:eigenvector matrix(columns correspond to vectors), eigval:eigenvalue vector
! inmaxcyc: max cycle, inthres: expected threshold
implicit real*8 (a-h,o-z)
integer,optional :: inmaxcyc
real*8,optional :: inthres
real*8 thres,mat(:,:),S(:,:),eigval(:)
real*8,allocatable :: R(:,:)
n=size(mat,1)
allocate(R(n,n))
maxcyc=200
thres=1D-9
if (present(inmaxcyc)) maxcyc=inmaxcyc
if (present(inthres)) thres=inthres
S=0
do i=1,n
	S(i,i)=1D0
end do
do k=1,maxcyc+1
	R=0
	do i=1,n
		R(i,i)=1D0
	end do
	i=1
	j=2
	do ii=1,n
		do jj=ii+1,n
			if (abs(mat(ii,jj))>abs(mat(i,j))) then
				i=ii
				j=jj
			end if
		end do
	end do
	if (abs(mat(i,j))<thres) exit
	if (k==maxcyc+1) write(*,*) "Note: Matrix diagonalization exceed max cycle before convergence"
	phi=atan(2*mat(i,j)/(mat(i,i)-mat(j,j)))/2D0
	R(i,i)=cos(phi)
	R(j,j)=R(i,i)
	R(i,j)=-sin(phi)
	R(j,i)=-R(i,j)
	mat=matmul(matmul(transpose(R),mat),R)
	S=matmul(S,R)
end do
do i=1,n
	eigval(i)=mat(i,i)
end do
end subroutine


!!------------ Diagonalize a symmetry matrix 
!Repack the extremely complex "DSYEV" routine in lapack to terse form
!if istat/=0, means error occurs
subroutine diagsymat(mat,eigvecmat,eigvalarr,istat)
integer istat
real*8 mat(:,:),eigvecmat(:,:),eigvalarr(:)
real*8,allocatable :: lworkvec(:)
isize=size(mat,1)
allocate(lworkvec(3*isize-1))
call DSYEV('V','U',isize,mat,isize,eigvalarr,lworkvec,3*isize-1,istat)
eigvecmat=mat
mat=0D0
forall (i=1:isize) mat(i,i)=eigvalarr(i)
end subroutine


!!------------ Diagonalize a general matrix 
!Repack the extremal complex "DGEEV" routine in lapack to terse form
!eigvecmat is right eigenvector matrix
!eigvalarr is real part of eigenvalue, imaginary parts are discarded
!if istat/=0, means error appears
subroutine diaggemat(mat,eigvecmat,eigvalarr,istat)
integer istat,lwork
real*8 mat(:,:),eigvecmat(:,:),eigvalarr(:),tmpmat(1,1)
real*8,allocatable :: lworkvec(:),eigvalimgarr(:)
isize=size(mat,1)
lwork=8*isize !4*isize is enough, but for better performance we use larger size
allocate(lworkvec(lwork),eigvalimgarr(isize))
call DGEEV('N','V',isize,mat,isize,eigvalarr,eigvalimgarr,tmpmat,1,eigvecmat,isize,lworkvec,lwork,istat)
mat=0D0
forall (i=1:isize) mat(i,i)=eigvalarr(i)
end subroutine


!!------------ Singular value decomposition (SVD) for general matrix: A = U * SIGMA * transpose(V)
!Repack the extremal complex "DGESVD" and "DGESDD" routine in lapack to terse form
!imethod=1: Standard algorithm (DGESVD) =2: Divide-and-conquer algorithm (DGESDD), faster but need more workspace and not as stable as DGESVD
!The inputted A will not be rewritten. singval is diagonal terms of SIGMA, sorted in descending order
!If A is a m*n matrix, then inputted U should be m*m and V should be n*n
!If info=0, that means the calculation is successfully finished
!Note that the returned matV is exactly V, rather than tranpose(V)
subroutine SVDmat(imethod,matA,matU,matV,singval,info)
real*8 matA(:,:),matU(:,:),matV(:,:),singval(:)
real*8,allocatable :: matAbk(:,:),workarr(:)
integer info,imethod
integer,allocatable :: workarri(:)
mdim=size(matA,1)
ndim=size(matA,2)
allocate(matAbk(mdim,ndim))
matAbk=matA
allocate(workarr(1))
!Obtain optimal lwork size
if (imethod==1) then
	call dgesvd('A','A',mdim,ndim,matA,mdim,singval,matU,mdim,matV,ndim,workarr,-1,info) !Standard algorithm
else
	call dgesdd('A',mdim,ndim,matA,mdim,singval,matU,mdim,matV,ndim,workarr,-1,workarri,info) !Divide-and-conquer algorithm, faster but need work workspace
end if
lwork=nint(workarr(1))
deallocate(workarr)
allocate(workarr(lwork))
!Below routines return transpose(V) and destory A. Note that dgesvd returns transpose(V), while current routine returns V, so final we do matV=transpose(matV)
if (imethod==1) then
	call dgesvd('A','A',mdim,ndim,matA,mdim,singval,matU,mdim,matV,ndim,workarr,lwork,info) !Standard algorithm
else
	allocate(workarri(8*min(mdim,ndim))) !Only needed by dgesdd
	call dgesdd('A',mdim,ndim,matA,mdim,singval,matU,mdim,matV,ndim,workarr,lwork,workarri,info) !Divide-and-conquer algorithm, faster but need work workspace
end if
matA=matAbk
matV=transpose(matV)
end subroutine


!!----------- Pseudo inversion (Moore–Penrose inverse) for an inputted matrix (mat) and returns matinv
!If mat has dimension of n*m, then matinv has dimension of m*n
!Realized based on singular value decomposition (SVD), ref.: https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse
!The matrix is assumed to be real
subroutine pseudoinverse(mat,matinv)
implicit real*8 (a-h,o-z)
real*8 mat(:,:),matinv(:,:)
real*8,allocatable :: singval(:),matU(:,:),matV(:,:),matSiginv(:,:)
n1=size(mat,1)
n2=size(mat,2)
nmin=min(n1,n2)
allocate(singval(nmin),matU(n1,n1),matV(n2,n2),matSiginv(n2,n1))
call SVDmat(1,mat,matU,matV,singval,info) !info=0: Successful
if (idebug==1.and.info/=0) write(*,*) "Warning: Singular value decomposition (SVD) is failed!"
matSiginv=0
do i=1,nmin
    if (abs(singval(i))>1D-10) matSiginv(i,i)=1/singval(i)
end do
matinv=matmul(matmul(matV,matSiginv),transpose(matU))
end subroutine


!!-------Return matrix multiplication of two double float matrices. This is a wrapper of DGEMM routine in BLAS so that invoking could be easier
!In fact, I found at -O2 level with /Qopt-matmul option, the matmul() is even much faster than this, therefore this function may be useless ,
!However, when MKL is linked, due to parallellization, this routine is much faster than matmul() with /Qopt-matmul
! nArow is the number of rows of matrix A, namely size(matA,1). If tranBin=1, it should be size(matA,2)
! nBcol is the number of columns of matrix B, namely size(matB,2). If tranBin=1, it should be size(matB,1)
! The returned matrix is (nArow,nBcol) dimension
! tranAin (tranBin) is optional, if it is 1 rather than 0, then matA (matB) will be transposed before doing the multiplication
!This statement is wrong: If the inputted matrix is not square, the transpose should be done using transpose() prior to input
function matmul_blas(matA,matB,nArow,nBcol,tranAin,tranBin)
real*8 matA(:,:),matB(:,:),matmul_blas(nArow,nBcol)
integer nArow,nBcol
integer,optional :: tranAin,tranBin
character tranA,tranB
nAcol=size(matA,2)
tranA='N'
lda=nArow
tranB='N'
ldb=nAcol
if (present(tranAin)) then
	if (tranAin==1) then
		tranA='T'
		lda=nAcol
	end if
end if
if (present(tranBin)) then
	if (tranBin==1) then
		tranB='T'
		ldb=nBcol
	end if
end if
call dgemm(tranA,tranB,nArow,nBcol,nAcol,1D0,matA,lda,matB,ldb,0D0,matmul_blas,nArow)
end function


!!--------------- A function to output inverted matrix, inputted matrix will not be affected. Essentially is a warpper of KROUT
function invmat(mat,N)
integer N,ierr
real*8 :: mat(N,N),invmat(N,N),tmpvec(N)
invmat=mat
call KROUT(0,N,0,invmat,N,tmpvec,N,ierr)
end function
!!--------------- A routine to invert matrix. The inputted matrix will be taken placed by inverted matrix. Essentially is a warpper of KROUT
subroutine invmatsub(mat,N)
integer N,ierr
real*8 :: mat(N,N),tmpvec(N)
call KROUT(0,N,0,mat,N,tmpvec,N,ierr)
end subroutine
!Taken and adapted from krout.f, which can be downloaded at http://jblevins.org/mirror/amiller/
!-----------------------------------------------------------------------
!  CROUT PROCEDURE FOR INVERTING MATRICES AND SOLVING EQUATIONS
!-----------------------------------------------------------------------
!  A IS A MATRIX OF ORDER N WHERE N IS GREATER THAN OR EQUAL TO 1.
!  IF MO = 0 THEN THE INVERSE OF A IS COMPUTED AND STORED IN A.
!  IF MO IS NOT 0 THEN THE INVERSE IS NOT COMPUTED.

!  IF M IS GREATER THAN 0 THEN B IS A MATRIX HAVING N ROWS AND M COLUMNS.
!  IN THIS CASE AX = B IS SOLVED AND THE SOLUTION X IS STORED IN B.
!  IF M=0 THEN THERE ARE NO EQUATIONS TO BE SOLVED.
!  N.B. B is passed as a VECTOR not as a matrix.

!  KA = THE LENGTH OF THE COLUMNS OF THE ARRAY A
!  KB = THE LENGTH OF THE COLUMNS OF THE ARRAY B (IF M > 0)

!  IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. WHEN
!  THE ROUTINE TERMINATES IERR HAS ONE OF THE FOLLOWING VALUES ...
!     IERR =  0   THE REQUESTED TASK WAS PERFORMED.
!     IERR = -1   EITHER N, KA, OR KB IS INCORRECT.
!     IERR =  K   THE K-TH PIVOT ELEMENT IS 0.
!-----------------------------------------------------------------------
! Adapted from the routine KROUT in the NSWC Math. Library by Alan Miller
! Latest revision - 3 August 1998
subroutine KROUT(MO, N, M, A, KA, B, KB, IERR)
implicit none
integer, intent(in)                            :: MO
integer, intent(in)                            :: N
integer, intent(in)                            :: M
real*8, intent(in out), dimension(:,:) :: A     ! a(ka,n)
integer, intent(in)                            :: KA
real*8, intent(in out), dimension(:)   :: B
integer, intent(in)                            :: KB
integer, intent(out)                           :: IERR
integer, allocatable, dimension(:)         :: INDX
real*8, allocatable, dimension(:)  :: TEMP
integer        :: I, J, JP1, K, KJ, KM1, KP1, L, LJ, MAXB, NJ, NMJ, NMK, NM1, ONEJ
real*8 :: D, DSUM, P, T
real*8, parameter :: ZERO = 0D0, ONE = 1D0
if (N < 1 .or. KA < N) then
  IERR = -1
  return
end if
if (M > 0 .and. KB < N) then
  IERR = -1
  return
end if
IERR = 0
if (N < 2) then
  D = A(1,1)
  if (D == ZERO) then
    IERR = N
    return
  end if
  if (MO == 0) A(1,1) = ONE / D
  if (M <= 0) return
  MAXB = KB*M
  do KJ = 1,MAXB,KB
    B(KJ) = B(KJ)/D
  end do
  return
end if
if (MO == 0) then
  allocate( INDX(N-1), TEMP(N) )
end if
NM1 = N - 1
do K = 1,NM1
  KP1 = K + 1
  P = abs(A(K,K))
  L = K
  do I = KP1,N
    T = abs(A(I,K))
    if (P >= T) cycle
    P = T
    L = I
  end do
  if (P == ZERO) then
    IERR = K
    return
  end if
  P = A(L,K)
  if (MO == 0) then
    INDX(K) = L
  end if
  if (K /= L) then
    do J = 1,N
      T = A(K,J)
      A(K,J) = A(L,J)
      A(L,J) = T
    end do
    if (M > 0) then
      KJ = K
      LJ = L
      do J = 1,M
        T = B(KJ)
        B(KJ) = B(LJ)
        B(LJ) = T
        KJ = KJ + KB
        LJ = LJ + KB
      end do
    end if
  end if
  if (K <= 1) then
    do J = KP1,N
      A(K,J) = A(K,J)/P
    end do
  else
    do J = KP1,N
      DSUM = A(K,J) - dot_product( A(K,1:KM1), A(1:KM1,J) )
      A(K,J) = DSUM / P
    end do
  end if
  do I = KP1,N
    DSUM = A(I,KP1) - dot_product( A(I,1:K), A(1:K,KP1) )
    A(I,KP1) = DSUM
  end do
  KM1 = K
end do
if (A(N,N) == ZERO) then
  IERR = N
  return
end if
if (M > 0) then
  MAXB = KB*M
  do ONEJ = 1,MAXB,KB
    KJ = ONEJ
    B(KJ) = B(KJ)/A(1,1)
    do K = 2,N
      KJ = KJ + 1
      DSUM = B(KJ)
      KM1 = K - 1
      LJ = ONEJ
      do L = 1,KM1
        DSUM = DSUM - A(K,L)*B(LJ)
        LJ = LJ + 1
      end do
      B(KJ) = DSUM / A(K,K)
    end do
  end do
  do NJ = N,MAXB,KB
    KJ = NJ
    do NMK = 1,NM1
      K = N - NMK
      LJ = KJ
      KJ = KJ - 1
      DSUM = B(KJ)
      KP1 = K + 1
      do L = KP1,N
        DSUM = DSUM - A(K,L)*B(LJ)
        LJ = LJ + 1
      end do
      B(KJ) = DSUM
    end do
  end do
end if
if (MO /= 0) return
do J = 1,NM1
  A(J,J) = ONE / A(J,J)
  JP1 = J + 1
  do I = JP1,N
    DSUM = dot_product( A(I,J:I-1), A(J:I-1,J) )
    A(I,J) = -DSUM / A(I,I)
  end do
end do
A(N,N) = ONE / A(N,N)
do NMK = 1,NM1
  K = N - NMK
  KP1 = K + 1
  do J = KP1,N
    TEMP(J) = A(K,J)
    A(K,J) = ZERO
  end do
  do J = 1,N
    DSUM = A(K,J) - dot_product( TEMP(KP1:N), A(KP1:N,J) )
    A(K,J) = DSUM
  end do
end do
do NMJ = 1,NM1
  J = N - NMJ
  K = INDX(J)
  if (J == K) cycle
  do I = 1,N
    T = A(I,J)
    A(I,J) = A(I,K)
    A(I,K) = T
  end do
end do
if (MO == 0) deallocate( INDX, TEMP )
end subroutine


!------- Calculate how much is a square matrix deviates from identity matrix
!error=∑[i,j]abs( abs(mat(i,j))-δ(i,j) )
real*8 function identmaterr(mat)
implicit real*8 (a-h,o-z)
real*8 mat(:,:)
nsize=size(mat,1)
identmaterr=0D0
do i=1,nsize
	do j=1,nsize
		if (i==j) then
			identmaterr=identmaterr+abs(abs(mat(i,j))-1D0)
		else
			identmaterr=identmaterr+abs(mat(i,j))
		end if
	end do
end do
end function


!------- Return maximal deviation of given square matrix from identity matrix
!errdiag is maximum absolute deviation of diagonal elements from 1, idiag is its index
!errndiag is maximum absolute deviation of non-diagonal elements from 0, indiag,jndiag are its row and column indices
subroutine identmatmaxerr(mat,errdiag,idiag,errndiag,indiag,jndiag)
implicit real*8 (a-h,o-z)
real*8 mat(:,:)
nsize=size(mat,1)
errdiag=0D0
idiag=1
errndiag=0D0
indiag=1
jndiag=1
do i=1,nsize
	do j=1,nsize
		if (i==j) then
			absdev=abs(mat(i,j)-1)
            if (absdev>errdiag) then
				idiag=i
                errdiag=absdev
            end if
		else
			absdev=abs(mat(i,j))
            if (absdev>errndiag) then
				indiag=i
				jndiag=j
                errndiag=absdev
            end if
		end if
	end do
end do
end subroutine


!----- Convert a square matrix to an array. imode=1/2/3: Full matrix; Lower half matrix; Upper half matrix
!For mode=1,2, "arr" should be nsize*(nsize+1)/2
subroutine mat2arr(mat,arr,imode)
implicit real*8 (a-h,o-z)
real*8 mat(:,:),arr(:)
nsize=size(mat,1)
itmp=0
if (imode==1) then !Full matrix
	do i=1,nsize
		do j=1,nsize
			itmp=itmp+1
			arr(itmp)=mat(i,j)
		end do
	end do
else if (imode==2) then !Lower half matrix
	do i=1,nsize
		do j=1,i
			itmp=itmp+1
			arr(itmp)=mat(i,j)
		end do
	end do
else !Upper half matrix
	do i=1,nsize
		do j=i,nsize
			itmp=itmp+1
			arr(itmp)=mat(i,j)
		end do
	end do
end if
end subroutine



!----- Calculate outer product of two arrays "arr1" and "arr2" with size n1 and n2 to yield a new matrix "mat"(n1,n2)
!The inputted two arrays are considered as column arrays, the "arr2" will be transposed
subroutine vecextprod(mat,arr1,arr2,n1,n2)
implicit real*8 (a-h,o-z)
integer n1,n2
real*8 mat(n1,n2)
!$OMP PARALLEL DO SHARED(mat) PRIVATE(i,j) schedule(auto) NUM_THREADS(nthreads)
do i=1,n1
    do j=1,n2
        mat(i,j)=arr1(i)*arr2(j)
    end do
end do
!$OMP END PARALLEL DO
end subroutine


!----- Calculate matrix product of "mat1"(na*np) and "mat2"(np,nb) to yield "mat"(na,nb), using OpenMP
!imode=1: Normal product
!imode=2: matmul(mat1,transpose(mat2))
subroutine matprod(imode,mat,mat1,mat2)
real*8 mat(:,:),mat1(:,:),mat2(:,:)
integer imode
na=size(mat,1)
nb=size(mat,2)
np=size(mat1,2)
mat=0
if (imode==1) then
    !$OMP PARALLEL DO SHARED(mat) PRIVATE(ia,ib,ip) schedule(dynamic) NUM_THREADS(nthreads)
    do ia=1,na
        do ib=1,nb
            do ip=1,np
                mat(ia,ib)=mat(ia,ib)+mat1(ia,ip)*mat2(ip,ib)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
else if (imode==2) then
    !$OMP PARALLEL DO SHARED(mat) PRIVATE(ia,ib,ip) schedule(dynamic) NUM_THREADS(nthreads)
    do ia=1,na
        do ib=1,nb
            do ip=1,np
                mat(ia,ib)=mat(ia,ib)+mat1(ia,ip)*mat2(ib,ip)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
end if
end subroutine






!===============================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Misc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================================================!


!!------- Returns a unique second value. The difference between two times of invoking this routine is the consumed wall clock time
subroutine walltime(inow)
character nowdate*80,nowtime*10
integer :: inow,i1,i2,yyyy,mm,dd

call date_and_time(nowdate,nowtime)
!i1 converts yyyymmdd to a unique integer of second resolution
!see https://community.intel.com/t5/Intel-Fortran-Compiler/Convert-some-date-to-seconds-since-1970/m-p/805010
read(nowdate(1:4),*) yyyy
read(nowdate(5:6),*) mm
read(nowdate(7:8),*) dd
i1=dd-32075+1461*(yyyy+4800+(mm-14)/12)/4 + &
367*(mm-2-((mm-14)/12)*12)/12- &
3*((yyyy + 4900 + (mm - 14)/12)/100)/4

!i2 records seconds passed today
read(nowtime(1:2),*) inowhour
read(nowtime(3:4),*) inowminute
read(nowtime(5:6),*) inowsecond
i2=inowhour*3600+inowminute*60+inowsecond

inow=i1+i2
end subroutine



!!----- Find the position of specific value in cube
subroutine findvalincub(cubfile,value,i,j,k)
real*8 cubfile(:,:,:),value
integer i,j,k
do ii=1,size(cubfile,1)
	do jj=1,size(cubfile,2)
		do kk=1,size(cubfile,3)
			if (cubfile(ii,jj,kk)==value) then
				i=ii
				j=jj
				k=kk
			end if
		end do
	end do
end do
end subroutine



!!------ Print matrix in a format similar to Gaussian. The number of columns is 5 (unadjustable)
!mat: The matrix to be printed
!  Below are optional
!Label: The title information. If the content is empty, title will not be printed
!insemi: If 1, print as lower trigonal matrix. Default is 0 (full matrix)
!form: Format to print, total width must be 14 characters. Default is 1PE14.5E3
!fildid: Output destination, 6 corresponds to outputting to screen
!usern1 and usern2: Dimensions of the matrix, default or -1 means determining automatically
!inncol: seems controls spacing between number labels of each frame
!formindex: Format of index. Default is "i8,6x". Note that if you manually set inncol, you should also set this to broaden or narrow index spacing
subroutine showmatgau(mat,label,insemi,form,fileid,usern1,usern2,inncol,formindex)
implicit real*8 (a-h,o-z)
real*8 :: mat(:,:)
character(*),optional :: label,form,formindex
integer,optional :: insemi,fileid,usern1,usern2,inncol
integer semi,ides,ncol
semi=0
ides=6
ncol=5
n1=size(mat,1)
n2=size(mat,2)
if (present(usern1)) then
	if (usern1/=-1) n1=usern1
end if
if (present(usern2)) then
	if (usern2/=-1) n2=usern2
end if
if (present(insemi)) semi=insemi
if (present(fileid)) ides=fileid
if (present(inncol)) ncol=inncol
if (present(label)) then
	if (label/='') then
		nspc=int((79-len(label))/2D0)
		write(ides,"(a)",advance='no') " "
		do i=1,nspc-2
			write(ides,"(a)",advance='no') '*'
		end do
		write(ides,"(a)",advance='no') ' '//label//' '
		do i=1,nspc-2
			write(ides,"(a)",advance='no') '*'
		end do
		write(ides,*)
    end if
end if
nf=ceiling(n2/float(ncol))
do i=1,nf !How many frame
	ns=(i-1)*5+1 !Absolute starting index to read
	if (i/=nf) then
        ne=i*ncol !Absolute ending index to read
	else
        ne=n2
    end if
	!Write basis number in separate line
	write(ides,"(6x)",advance='no')
	do j=ns,ne
		if (present(formindex)) then
			write(ides,'('//formindex//')',advance='no') j
		else
			write(ides,"(i8,6x)",advance='no') j
		end if
	end do
	write(ides,*)
	!Write content in each regular line
	do k=1,n1
		if (k<ns.and.semi==1) cycle !The lines have been outputted are skipped
		write(ides,"(i6)",advance='no') k
		do j=ns,ne
			if (semi==1.and.k<j) cycle !Upper trigonal element were passed
			if (present(form)) then
				write(ides,'('//form//')',advance='no') mat(k,j)			
			else
				write(ides,"(1PE14.5E3)",advance='no') mat(k,j)
			end if
		end do
		write(ides,*) !Change to next line
	end do
end do
end subroutine



!!------- Read matrix in Gaussian output file, such as that printed by IOp(3/33=1) (or may be matrix printed by e.g. NBO)
subroutine readmatgau(fileid,mat,insemi,inform,inskipcol,inncol,innspace,iostat)
!e.g. Lower trigonal matrix: readmatgau(10,tmpmat,1,"D14.6",7,5)   full matrix: readmatgau(10,tmpmat,0,"?",7,5)
!mat: The matrix to return
!Below are optional options
!insemi: 1 means the matrix is a lower trigonal matrix. Default is full matrix
!inform: Format used to read data, default is D14.6. e.g. "f8.4 " (must be 5 character!!!). inform="?" means use free format to read, recommended!
!inskipcol: Number of column (marker) to be skipped at the beginning of each row, default is 7
!inncol: Number of columns in each line, default is 5
!innspace: Number of useless lines in between each frame, default is 1
!iostat: Returned reading status. /=0 means error occurs

!Before use, loclabel should be used to move reading position to title line of the matrix，namely move to "*** Overlap ***"
! *** Overlap ***
!                 1             2             3             4             5
!       1  0.100000D+01
!       2  0.236704D+00  0.100000D+01
!       ...
!    Nbasis ...
implicit real*8 (a-h,o-z)
real*8 :: mat(:,:)
character(len=*),optional :: inform
integer,optional :: inskipcol,inncol,insemi,innspace,iostat
character form*7,c80tmp*79,c200tmp*200
integer fileid, skipcol,ncol,nspace,semi
form="(D14.6)" !Suitable for Sbas,Kinene,Potene,Hcore by IOp(3/33=1) ,Fockmat,Densmat by IOp(5/33=3)
skipcol=7
ncol=5
nspace=1
semi=0
if (present(iostat)) iostat=0
if (present(inform)) form='('//inform//')'
if (present(inskipcol)) skipcol=inskipcol
if (present(inncol)) ncol=inncol
if (present(innspace)) nspace=innspace
if (present(insemi)) semi=insemi
n1=size(mat,1)
n2=size(mat,2)
nf=ceiling(n2/float(ncol))
mat=0D0
read(fileid,*) !Skip title line
do i=1,nf !Number of frames
	ns=(i-1)*ncol+1 !Absolute starting index to read
	if (i/=nf) then !Absolute ending index to read
        ne=i*ncol
	else
        ne=n2
    end if
	do itmp=1,nspace !Skip rows when reading each new frame
		read(fileid,*)
	end do
	do k=1,n1 !Scan rows in each frame
 		!read(fileid,"(a)") c80tmp
 		!write(*,"(a)") c80tmp
 		!backspace(fileid)
        if (semi==1.and.ns>k) cycle
        if (form/="(?)") then !Use fixed format to read
		    do itmp=1,skipcol !Skip marker columns in each row
			    read(fileid,"(1x)",advance='no')
		    end do
		    do j=ns,ne !Scan elements in each row
			    if (semi==1.and.j>k) cycle
			    read(fileid,form,advance='no',iostat=ierror) mat(k,j)
                if (ierror/=0) then
                    if (present(iostat)) iostat=ierror
                    return
                end if
 			    !write(*,*) i,k,j,mat(k,j)
		    end do
		    read(fileid,*)
        else !Use free format to read, assume the line is up to 200 characters
            read(fileid,"(a)") c200tmp
            if (semi==0) then
                read(c200tmp(1+skipcol:),*,iostat=ierror) mat(k,ns:ne)
            else
                read(c200tmp(1+skipcol:),*,iostat=ierror) mat(k,ns:min(k,ne))
            end if
            if (ierror/=0) then
                if (present(iostat)) iostat=ierror
                return
            end if
        end if
	end do
end do
if (semi==1) then !For lower trigonal matrix, it is assumed to be symmetric. Now fill another part
	mat=mat+transpose(mat)
	do i=1,n1
		mat(i,i)=mat(i,i)/2D0
	end do
end if
end subroutine


!!--------------------- Determine how many lines in the fileid
!If imode==1, space line will be regarded as the sign of end of file. If imode==2, will count actual number of lines in the file
integer function totlinenum(fileid,imode)
integer fileid,ierror,imode
character c80*80
totlinenum=0
rewind(fileid)
do while(.true.)
	read(fileid,"(a)",iostat=ierror) c80
	if (imode==1) then
		if (ierror/=0.or.c80==" ") exit
	else if (imode==2) then
		if (ierror/=0) exit
	end if
	totlinenum=totlinenum+1
end do
rewind(fileid)
end function


!!-------- Locate the line where the label first appears in fileid
!Return ifound=1 if found the label, else return 0
!Default is rewind, if irewind=0 then will not rewind
!If the current line just has the label, calling this subroutine will do nothing
!maxline define the maximum number of lines that will be searched, default is search the whole file
subroutine loclabel(fileid,label,ifound,irewind,maxline)
integer fileid,ierror
integer,optional :: ifound,irewind,maxline
character c200tmp*200
character(len=*) label
if (.not.present(irewind)) then
	rewind(fileid)
else
	if (irewind==1) rewind(fileid)
end if
if (.not.present(maxline)) then
	do while(.true.)
		read(fileid,"(a)",iostat=ierror) c200tmp
		if (ierror/=0) exit
		if (index(c200tmp,label)/=0) then
			backspace(fileid)
			if (present(ifound)) ifound=1 !Found result
			return
		end if
	end do
else
	do iline=1,maxline
		read(fileid,"(a)",iostat=ierror) c200tmp
		if (ierror/=0) exit
		if (index(c200tmp,label)/=0) then
			backspace(fileid)
			if (present(ifound)) ifound=1 !Found result
			return
		end if
	end do
end if
if (present(ifound)) ifound=0
end subroutine


!-------- Locate to the final label, and meantime returns the number of matches. Based on "loclabel"
subroutine loclabelfinal(fileid,label,nfound)
integer fileid,nfound,ifound
character(len=*) label
nfound=0
rewind(fileid)
do while(.true.)
    call loclabel(fileid,label,ifound,0)
    if (ifound==0) then
        exit
    else
        nfound=nfound+1
        read(fileid,*)
    end if
end do
rewind(fileid)
do ifound=1,nfound
    call loclabel(fileid,label,inouse,0)
    if (ifound<nfound) read(fileid,*)
end do
end subroutine


!!-------- Skip specific number of lines in specific fileid
subroutine skiplines(id,nskip)
integer id,nskip
do i=1,nskip
	read(id,*)
end do
end subroutine


!!---------------- Calculate factorial
integer function ft(i)
integer i
ft=i
if (i==0) ft=1
do j=i-1,1,-1
	ft=ft*j
end do
end function


!!---- Calculate gamma(Lval+1/2), see http://en.wikipedia.org/wiki/Gamma_function
real*8 function gamma_ps(n)
use defvar
integer n
gamma_ps=ft(2*n)*dsqrt(pi)/4**n/ft(n)
end function


!!-------- Get all combinations of any ncomb elements of array, which length is ntot
!outarray(A,B) is output array, A is the length and must be ntot!/ncomb!/(ntot-ncomb)!, B is generated array, should be equal to ncomb
subroutine combarray(array,ntot,ncomb,outarray)
integer array(ntot),ntot,ncomb,idxarr(ncomb),outarray(:,:)
ipos=ncomb !Current position in the array
forall (i=1:ncomb) idxarr(i)=i !Used to record index
ioutput=1
ncount=0
do while(ipos>0)
	if (ioutput==1) then
		ncount=ncount+1
		outarray(ncount,:)=array(idxarr(:))
	end if
	ioutput=0
	idxarr(ipos)=idxarr(ipos)+1
	if (idxarr(ipos)>ntot) then
		ipos=ipos-1 !Go back to last position
		cycle
	end if
	if (ipos<ncomb) then
		ipos=ipos+1
		idxarr(ipos)=idxarr(ipos-1)
		cycle
	end if
	if (ipos==ncomb) ioutput=1
end do
end subroutine


!!--------- Convert XY scatter data to density distribution
!xarr and yarr records the points. nlen is array length
!mat is the outputted matrix, matnx and matny are its number of element in X and Y
!x/ymin, x/ymax are lower and upper limit of "mat", the X and Y ranges contain nvalx and nvaly data
!e.g. n=5
!   |   1   |   2   |   3   |   4   |   5    |
! xmin-------------------------------------xmax
subroutine xypt2densmat(xarr,yarr,nlen,mat,nvalx,nvaly,xmin,xmax,ymin,ymax)
integer nlen,nvalx,nvaly
real*8 xarr(nlen),yarr(nlen),mat(nvalx,nvaly),xmin,xmax,ymin,ymax
mat=0D0
spcx=(xmax-xmin)/nvalx
spcy=(ymax-ymin)/nvaly
!If enable parallel, program often prompts memory is not enough, I don't know how to solve this
! !$OMP PARALLEL DO SHARED(mat) PRIVATE(ix,iy,ipt,xlow,xhigh,ylow,yhigh) schedule(dynamic) NUM_THREADS(nthreads)
do ix=1,nvalx
	xlow=xmin+(ix-1)*spcx
	xhigh=xmin+ix*spcx
	do iy=1,nvaly
		ylow=ymin+(iy-1)*spcy
		yhigh=ymin+iy*spcy
		do ipt=1,nlen
			if (xarr(ipt)>xlow.and.xarr(ipt)<=xhigh.and.yarr(ipt)>ylow.and.yarr(ipt)<=yhigh) then
				mat(ix,iy)=mat(ix,iy)+1D0
			end if
		end do
	end do
end do
! !$OMP END PARALLEL DO
end subroutine



!--------- Determine the present file is output file of which code
!The file must has been opended as "ifileid"
!iprog: 1=Outputted by Gaussian, 2=Outputted by ORCA, 3=Outputted by GAMESS-US, 4=Outputted by Firefly, 5=CP2K, 6=xTB, 7=BDF, 0=Undetermined
!info (optional): 1 means output the type of this file
subroutine outputprog(ifileid,iprog,info)
integer ifileid,iprog
integer,optional :: info
call loclabel(ifileid,"Gaussian, Inc",ifound,maxline=500)
if (ifound==0) call loclabel(ifileid,"Entering Gaussian System",ifound,maxline=200)
if (ifound==1) then
    iprog=1
    if (present(info)) write(*,*) "Note: This file is recognized as a Gaussian output file"
    return
end if
call loclabel(ifileid,"O   R   C   A",ifound,maxline=500)
if (ifound==1) then
    iprog=2
    if (present(info)) write(*,*) "Note: This file is recognized as an ORCA output file"
    return
end if
call loclabel(ifileid,"GAMESS VERSION =",ifound,maxline=500)
if (ifound==1) then
    iprog=3
    if (present(info)) write(*,*) "Note: This file is recognized as a GAMESS-US output file"
    return
end if
call loclabel(ifileid,"Firefly Project",ifound,maxline=500)
if (ifound==1) then
    iprog=4
    if (present(info)) write(*,*) "Note: This file is recognized as a Firefly output file"
    return
end if
call loclabel(10,"CP2K|",ifound,maxline=500)
if (ifound==1) then
    iprog=5
    if (present(info)) write(*,*) "Note: This file is recognized as a CP2K output file"
    return
end if
call loclabel(10,"x T B",ifound,maxline=500)
if (ifound==1) then
    iprog=6
    if (present(info)) write(*,*) "Note: This file is recognized as a xTB output file"
    return
end if
call loclabel(10,"BDF",ifound,maxline=500)
if (ifound==1) then
    iprog=7
    if (present(info)) write(*,*) "Note: This file is recognized as a BDF output file"
    return
end if
iprog=0
if (present(info)) write(*,*) "Note: This file is treated as a plain text file"
end subroutine

!--------- Determine the present file is input file of which code when the answer cannot be determined from file extension
!The file must have been opened as "ifileid"
!iprog=0: Cannot be identified  =1: CP2K  =2: ORCA  =3: Quantum ESPRESSO
subroutine inputprog(ifileid,iprog)
integer ifileid,iprog
iprog=0
call loclabel(ifileid,"&COORD ",ifound,maxline=10000) !Various information may present before &COORD, so maxline should be large
if (ifound==1) then
    iprog=1
    return
end if
call loclabel(ifileid,"* xyz ",ifound,maxline=500)
if (ifound==1) then
    iprog=2
    return
end if
call loclabel(ifileid,"&control ",ifound,maxline=1000)
if (ifound==0) call loclabel(ifileid,"&CONTROL ",ifound,maxline=1000)
if (ifound==1) then
    iprog=3
    return
end if
end subroutine


!Determine if the real space function currently to be studied (ifuncsel) involve ESP, and thus should call doinitlibreta first
!If should do, ifdoESP=.true., else ifdoESP=.false.
logical function ifdoESP(ifuncsel)
use defvar
integer ifuncsel
ifdoESP=.false.
if (ifuncsel==12) then
    ifdoESP=.true.
else if (ifuncsel==100) then
    if (iuserfunc==8) ifdoESP=.true.
    if (iuserfunc==14) ifdoESP=.true.
    if (iuserfunc==34) ifdoESP=.true.
    if (iuserfunc==39) ifdoESP=.true.
    if (iuserfunc>=60.and.iuserfunc<=68) ifdoESP=.true.
    if (iuserfunc==101) ifdoESP=.true.
    if (iuserfunc==102) ifdoESP=.true.
    if (iuserfunc==103) ifdoESP=.true.
    if (iuserfunc==999) ifdoESP=.true.
    if (iuserfunc==110) ifdoESP=.true.
    if (iuserfunc==111) ifdoESP=.true.
end if
end function



!showorbinfo1, showorbinfo2, showorbinfo3 print more information than showorbinfo, &
!and are mainly invoked by GUI subroutines, but can also be invoked by other subroutines

!----- Show all orbitals
subroutine showorbinfo1(id)
use defvar
integer,intent (in) :: id
character(len=3) :: orbtype(0:2)=(/ "A+B"," A "," B " /)
character :: symstr*6
symstr=" "
naorb=count(MOtype==1)
write(*,*) "Orbital list:"
do i=1,nmo
	if (allocated(MOsym)) symstr='('//MOsym(i)//')'
	if (wfntype==0.or.wfntype==2.or.wfntype==3) then
		write(*,"(' Orb:',i6,' Ene(au/eV):',f13.6,f13.4,' Occ:',f9.6,' Type:',a,1x,a)") &
		i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
	else
		if (MOtype(i)==1) then
			write(*,"(i10,5x,' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
			i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
		else
			write(*,"(i6,' (',i6,')',' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
			i,i-naorb,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
		end if
	end if
end do
if (any(MOtype==2)) write(*,"(a)") " Note: For beta orbitals, &
&the index in the parenthese shown above is the index counted from the first beta orbital"
end subroutine

!----- Show orbitals up to LUMO+10, works for wfntype==0,1,2
subroutine showorbinfo2(id)
use defvar
integer,intent (in) :: id
character(len=3) :: orbtype(0:2)=(/ "A+B"," A "," B " /)
character :: symstr*6
symstr=" "
naorb=count(MOtype==1)
if (wfntype==0.or.wfntype==2) then
	write(*,*) "Orbital list:"
	do nmoend=1,nmo
		if (MOocc(nmoend)==0D0) exit
	end do
	nmoend=nmoend+10
	if (nmoend>nmo) nmoend=nmo
	do i=1,nmoend
		if (allocated(MOsym)) symstr='('//MOsym(i)//')'
		write(*,"(' Orb:',i6,' Ene(au/eV):',f13.6,f13.4,' Occ:',f9.6,' Type:',a,1x,a)") &
		i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
	end do
else if (wfntype==1) then
	do iLUMOa=1,nmo
		if (MOocc(iLUMOa)==0) exit
	end do
	do iLUMOb=nmo,1,-1
		if (MOocc(iLUMOb)==1) exit
	end do
	iLUMOb=iLUMOb+1
	do ibeta=1,nmo
		if (MOtype(ibeta)==2) exit
	end do
	iaend=iLUMOa+10
	if (iaend>=ibeta) iaend=ibeta-1
	ibend=iLUMOb+10
	if (ibend>nmo) ibend=nmo
	write(*,*) "Alpha orbital list:"
	do i=1,iaend
		if (allocated(MOsym)) symstr='('//MOsym(i)//')'
		write(*,"(i10,5x,' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
		i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
	end do
	write(*,*) "Beta orbital list:"
	do i=ibeta,ibend
		if (allocated(MOsym)) symstr='('//MOsym(i)//')'
		write(*,"(i6,' (',i6,')',' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
		i,i-naorb,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
	end do
	write(*,"(a)") " Note: For beta orbitals, &
	&the index in the parenthese shown above is the index counted from the first beta orbital"
end if
end subroutine

!----- Show all occupied orbitals
subroutine showorbinfo3(id)
use defvar
integer,intent (in) :: id
character(len=3) :: orbtype(0:2)=(/ "A+B"," A "," B " /)
character symstr*6
symstr=" "
naorb=count(MOtype==1)
write(*,*) "Orbital list:"
do i=1,nmo
	if (MOocc(i)==0) cycle
	if (allocated(MOsym)) symstr='('//MOsym(i)//')'
	if (wfntype==0.or.wfntype==2.or.wfntype==3) then
		write(*,"(' Orb:',i6,' Ene(au/eV):',f13.6,f13.4,' Occ:',f9.6,' Type:',a,1x,a)") &
		i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
	else
		if (MOtype(i)==1) then
			write(*,"(i10,5x,' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
			i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
		else
			write(*,"(i6,' (',i6,')',' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
			i,i-naorb,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
		end if
	end if
end do
if (any(MOtype==2)) write(*,"(a)") " Note: For beta orbitals, &
&the index in the parenthese shown above is the index counted from the first beta orbital"
end subroutine



!!--------- Calculate switching function based on error function at x
!Decrease from x= 0 to inf
!xhalf is the position where switching function is 0.5
!xscale is the scale factor of x value. The larger the smoother. When it is 0.25, the switching function will decay from 1 to 0 approximately within [xhalf-0.5,xhalf+0.5]
!Approximate sharpness relationship: xscale=0.3 is basically iter=3 of Becke, xscale=0.5 is basically iter=1 of Becke
real*8 function switch_erf(x,xhalf,xscale)
implicit real*8 (a-h,o-z)
real*8 x,xhalf,xscale
switch_erf=1-(0.5D0*erf((x-xhalf)/xscale)+0.5D0)
end function

!!--------- Calculate switching function based on Becke's function at x
!Decrease from x= xhalf-1 (value=1) to x= xhalf+1 (value=0)
!xhalf is the position where switching function is 0.5
!nBeckeiter is number of iterations (>=1), the larger the sharper the variation. If =1, variation will be linear
real*8 function switch_Becke(x,xhalf,nBeckeiter)
implicit real*8 (a-h,o-z)
real*8 x,xhalf
integer nBeckeiter
tmps=x-xhalf
if (tmps>1) then
	switch_Becke=0
else if (tmps<-1) then
	switch_Becke=1
else
	do iter=1,nBeckeiter
		tmps=1.5D0*tmps-0.5D0*tmps**3
	end do
	switch_Becke=0.5D0*(1-tmps)
end if
end function

!!--------- Calculate switching function based on (unnormalized) Gaussian function at x
!Decrease from 1 to 0 as x goes from 0 to inf
!gauFWHM is FWHM of Gaussian function
real*8 function switch_Gauss(x,gauFWHM)
real*8 x,gauFWHM,parmc
parmc=gauFWHM/2.35482D0 !Parameter c of Gaussian
switch_Gauss=exp(-x**2/(2*parmc**2))
end function


!Test code for examining various switching function from 0.0~8.0
!open(10,file="switch_erf.txt",status="replace")
!open(11,file="switch_Becke.txt",status="replace")
!open(12,file="switch_Gauss.txt",status="replace")
!xhalf=vdwr(6)/b2a !Radius of carbon
!do i=1,800
!	x=i*0.01D0
!	write(10,"(f8.3,4f12.6)") x,switch_erf(x,xhalf,0.25D0),switch_erf(x,xhalf,0.5D0),switch_erf(x,xhalf,1D0),switch_erf(x,xhalf,2D0)
!	write(11,"(f8.3,4f12.6)") x,switch_Becke(x,xhalf,3),switch_Becke(x,xhalf,2),switch_Becke(x,xhalf,1),switch_Becke(x,xhalf,0)
!	write(12,"(f8.3,3f12.6)") x,switch_Gauss(x,2*xhalf),switch_Gauss(x,1.5*xhalf),switch_Gauss(x,1.0*xhalf)
!end do
!close(10)
!close(11)
!close(12)
end module



!================= OUTSIDE MODULE =================



!!--------- Lagrange interpolation in 1D, produce interpolated value, 1st and 2nd derivatives
!NOTE: 4 adjacent data points will be used to interpolate
!ptpos and ptval are the position and value of the input array, npt is the number of its element. ptpos must vary from small to large
!r is the point to be studied, the resultant val, der1, der2 are its value, 1st and 2nd derivatives
!itype=1: only calculate value    =2: also calculate 1st-derv.    =3: also calculate 2nd-derv.
subroutine lagintpol(ptpos,ptval,npt,r,val,der1,der2,itype)
implicit real*8 (a-h,o-z)
integer npt,itype
real*8 ptpos(npt),ptval(npt),r,val,der1,der2
if (r<=ptpos(1)) then !Out of boundary
	der1=(ptval(2)-ptval(1))/(ptpos(2)-ptpos(1))
	der2=0D0
	val=ptval(1)-(ptpos(1)-r)*der1 !Use linear interpolation
	return
else if (r>=ptpos(npt)) then
	val=0D0 !Because this function is mainly used to interpolate radial atomic density, at long distance the value must be zero
	der1=0D0
	der2=0D0
! 	val=ptval(npt)
! 	der1=(ptval(npt)-ptval(npt-1))/(ptpos(npt)-ptpos(npt-1))
 	return
end if
do i=1,npt !Determine which four data points will be used to interpolation
	if (ptpos(i)>r) exit
end do
! if (i==npt+1) i=npt !i==npt+1 means i exceeded the last point
iend=i+1
istart=i-2
if (istart<1) then
	istart=istart+1
	iend=iend+1
else if (iend>npt) then
	iend=iend-1
	istart=istart-1
end if
!Calculate interpolated value
val=0D0
do m=istart,iend
	poly=1D0
	do j=istart,iend
		if (j==m) cycle
		poly=poly* (r-ptpos(j))/(ptpos(m)-ptpos(j))
	end do
	val=val+ptval(m)*poly
end do
!Calculate interpolated 1st-derv.
if (itype<2) return
der1=0D0
do m=istart,iend
	suml=0D0
	do l=istart,iend
		if (l==m) cycle
		poly=1D0
		do j=istart,iend
			if (j==m.or.j==l) cycle
			poly=poly* (r-ptpos(j))/(ptpos(m)-ptpos(j))
		end do
		suml=suml+poly/(ptpos(m)-ptpos(l))
	end do
	der1=der1+ptval(m)*suml
end do
!Calculate interpolated 2nd-derv.
if (itype<3) return
der2=0D0
do m=istart,iend
	suml=0D0
	do l=istart,iend
		if (l==m) cycle
		sumn=0D0
		do n=istart,iend
			if (n==l.or.n==m) cycle
			poly=1D0
			do j=istart,iend
				if (j==l.or.j==m.or.j==n) cycle
				poly=poly* (r-ptpos(j))/(ptpos(m)-ptpos(j))
			end do
			sumn=sumn+poly/(ptpos(m)-ptpos(n))
		end do
		suml=suml+sumn/(ptpos(m)-ptpos(l))
	end do
	der2=der2+ptval(m)*suml
end do
end subroutine




!!--------- Linear interpolation in 1D
!r: position
!val: interpolated value
!ptpos: must vary from small to large
subroutine linintpol(ptpos,ptval,npt,r,val)
implicit real*8 (a-h,o-z)
real*8 ptpos(npt),ptval(npt),r,val
integer highpt

if (r<=ptpos(1)) then
	val=ptval(1)
    return
else if (r>=ptpos(npt)) then
	val=ptval(npt)
    return
end if
do ipt=1,npt
	if (ptpos(ipt)<r) then
		lowpt=ipt
    else
		highpt=ipt
        exit
    end if
end do
der=(ptval(highpt)-ptval(lowpt))/(ptpos(highpt)-ptpos(lowpt))
val=ptval(lowpt)+der*(r-ptpos(lowpt))
end subroutine




!--------- Show progress bar
subroutine showprog(inow,nall)
integer inow,nall
integer :: itmp=0
character c80tmp*80
iprog=int(dfloat(inow)/nall*50)
c80tmp=' Progress: ['
c80tmp(13:62)=repeat('#',iprog)
c80tmp(13+iprog:62)=repeat('-',50-iprog)
c80tmp(63:63)=']'
write(c80tmp(64:),"(f8.1,' %')") dfloat(inow)/nall*100
itmp=itmp+1
if (itmp==1) c80tmp(79:79)='-'
if (itmp==2) c80tmp(79:79)='\'
if (itmp==3) c80tmp(79:79)='|'
if (itmp==4) then
	c80tmp(79:79)='/'
	itmp=0
end if
write(*,"(2a$)") trim(c80tmp),char(13)
if (inow>=nall) write(*,*)
end subroutine


!--------- Run system command by inputting command string
subroutine runcommand(cmd)
use defvar
character(len=*) cmd
!Windows removes double quotation at the two sides of inputted string, therefore I add additional double quotation to protect those in the string
write(*,"(a)") " Running: "//trim(cmd)
if (isys==1) then
    call system(""""//cmd//"""")
else
    call system(cmd)
end if
end subroutine



!!-------- Get an integer argument from command line. e.g. call getarg_int("-nt",nthreads,ifound)
!argname: The label for which the value after it should be extracted
!argval: Returned value
!ifound=1/0: Found / Not found the argument
subroutine getarg_int(argname,argval,ifound)
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200
character(len=*) argname
integer ifound,argval
narg=command_argument_count()
iarg=1
ifound=0
do while(iarg<=narg)
    call get_command_argument(iarg,c200tmp)
    if (c200tmp==argname) then
        iarg=iarg+1
        call get_command_argument(iarg,c80tmp)
        read(c80tmp,*) argval
        ifound=1
        exit
    end if
    iarg=iarg+1
end do
end subroutine

!!--------- Same as above, but return float
subroutine getarg_float(argname,argval,ifound)
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200
character(len=*) argname
integer ifound
real*8 argval
narg=command_argument_count()
iarg=1
ifound=0
do while(iarg<=narg)
    call get_command_argument(iarg,c200tmp)
    if (c200tmp==argname) then
        iarg=iarg+1
        call get_command_argument(iarg,c80tmp)
        read(c80tmp,*) argval
        ifound=1
        exit
    end if
    iarg=iarg+1
end do
end subroutine

!!--------- Same as above, but return string
subroutine getarg_str(argname,argstr,ifound)
implicit real*8 (a-h,o-z)
character c80tmp,c200tmp*200
character(len=*) argname,argstr
integer ifound
narg=command_argument_count()
iarg=1 !The input file name must be the first argument
ifound=0
do while(iarg<=narg)
    call get_command_argument(iarg,c200tmp)
    if (c200tmp==argname) then
        iarg=iarg+1
        call get_command_argument(iarg,argstr)
        ifound=1
        exit
    end if
    iarg=iarg+1
end do
end subroutine

!!--------- Test if an argument is existed
subroutine testarg(argname,ifound)
implicit real*8 (a-h,o-z)
character c80tmp,c200tmp*200
character(len=*) argname
integer ifound
narg=command_argument_count()
iarg=1 !The input file name must be the first argument
ifound=0
do while(iarg<=narg)
    call get_command_argument(iarg,c200tmp)
    if (c200tmp==argname) then
        ifound=1
        exit
    end if
    iarg=iarg+1
end do
end subroutine



!!--------- Test if a directory is existed
subroutine inquire_dir(dirname,diralive)
use defvar
character(len=*) dirname
logical diralive
if (isys==1) then !Windows
	call execute_command_line("mkdir "//trim(dirname)//" >NUL 2>&1",exitstat=ierror)
	if (ierror==0) then !Directory can be created, so the directory doesn't exist before
		diralive=.false.
		!Delete temporarily created directory
		call execute_command_line("rmdir /Q /S "//trim(dirname)//" >NUL 2>&1")
	else
		diralive=.true.
	end if
else
    call execute_command_line("[ -d "//dirname//" ]",exitstat=ierror)
    if (ierror==0) then
		diralive=.true.
    else
		diralive=.false.
    end if
end if
end subroutine
