!!-------- Show cell information
subroutine showcellinfo
use defvar
implicit real*8 (a-h,o-z)
character c200tmp*200
real*8 rcellv1(3),rcellv2(3),rcellv3(3)

if (ifPBC>0) then
	write(*,*)
    write(*,*) "Cell information"
    write(*,*) "In Bohr:"
    write(*,"(' Cell vector 1,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") cellv1,dsqrt(sum(cellv1**2))
    if (ifPBC>1) then
        write(*,"(' Cell vector 2,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") cellv2,dsqrt(sum(cellv2**2))
        if (ifPBC>2) then
            write(*,"(' Cell vector 3,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") cellv3,dsqrt(sum(cellv3**2))
        end if
    end if
    write(*,*) "In Angstrom:"
    write(*,"(' Cell vector 1,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") cellv1*b2a,dsqrt(sum(cellv1**2))*b2a
    if (ifPBC>1) then
        write(*,"(' Cell vector 2,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") cellv2*b2a,dsqrt(sum(cellv2**2))*b2a
        if (ifPBC>2) then
            write(*,"(' Cell vector 3,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") cellv3*b2a,dsqrt(sum(cellv3**2))*b2a
        end if
    end if
    if (ifPBC==3) then
        call getrcellvec(rcellv1,rcellv2,rcellv3)
        write(*,*) "Reciprocal vectors in Angstrom^-1:"
        write(*,"(' Vector 1,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") rcellv1/b2a,dsqrt(sum(rcellv1**2))/b2a
        write(*,"(' Vector 2,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") rcellv2/b2a,dsqrt(sum(rcellv2**2))/b2a
        write(*,"(' Vector 3,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Norm:',f11.5)") rcellv3/b2a,dsqrt(sum(rcellv3**2))/b2a
    end if
    if (ifPBC==3) then
        call getcellabc(asize,bsize,csize,alpha,beta,gamma)
        write(*,"(' Cell angles:  Alpha=',f9.4,'  Beta=',f9.4,'  Gamma=',f9.4,' degree')") alpha,beta,gamma
        call calc_cellvol(cellvol)
        call calc_vecarea(cellv1*b2a,cellv2*b2a,area)
        write(*,"(' Area between a and b:',f12.4,' Angstrom^2')") area
        call calc_vecarea(cellv2*b2a,cellv3*b2a,area)
        write(*,"(' Area between b and c:',f12.4,' Angstrom^2')") area
        call calc_vecarea(cellv1*b2a,cellv3*b2a,area)
        write(*,"(' Area between a and c:',f12.4,' Angstrom^2')") area
        write(*,"(' Cell volume:',f16.4,' Bohr^3    (',f16.4,' Angstrom^3 )')") cellvol,cellvol*b2a**3
        dens=sum(atmwei(a%index))*amu2kg*1000  /(cellvol*b2a**3/1D24)
        write(*,"(' Density:',f10.5,' g/cm^3    (',f10.3,' kg/m^3 )')") dens,dens*1000
        write(*,*)
        write(c200tmp,"(a,3f10.5,3f8.3,a)") "pbc set {",asize,bsize,csize,alpha,beta,gamma," } -all"
        write(*,"(a,/,1x,a)") " Command of showing box in VMD program: (then run ""pbc box"")",trim(c200tmp)
        write(*,*)
        do idx=1,8
            call cellvertex(idx,vertx,verty,vertz)
            write(*,"(' Vertex',i2,' of the cell:',3f10.4,' Bohr')") idx,vertx,verty,vertz
        end do
    end if
end if
end subroutine

    
    
!!!------ Initialize some PBC information and settings. After using this, periodic code will be fully compatible with any case
subroutine init_PBC
use defvar
use util
ifdoPBCx=ifdoPBCx_in
ifdoPBCy=ifdoPBCy_in
ifdoPBCz=ifdoPBCz_in
PBCnx=PBCnx_in
PBCny=PBCny_in
PBCnz=PBCnz_in
!Force some directions do not consider periodicity
if (ifPBC==0) then !Isolated system
    ifdoPBCx=0
    ifdoPBCy=0
    ifdoPBCz=0
else if (ifPBC==1) then !1D system
    ifdoPBCy=0
    ifdoPBCz=0
    !Set temporary cellv2 and cellv3 orthogonal to v1, otherwise fractional coordinate cannot be generated from Cartesian coordinate
    call vecprod(cellv1(1),cellv1(2),cellv1(3),5D0,7D0,-4D0,cellv2(1),cellv2(2),cellv2(3)) !(5,7,-4) is an arbitrarily selected vector
    call vecprod(cellv1(1),cellv1(2),cellv1(3),cellv2(1),cellv2(2),cellv2(3),cellv3(1),cellv3(2),cellv3(3))
    cellv2=cellv2/dsqrt(sum(cellv2**2))
    cellv3=cellv3/dsqrt(sum(cellv3**2))
else if (ifPBC==2) then !2D system
    ifdoPBCz=0
    !Set temporary cellv3 orthogonal to v1 and v2, otherwise fractional coordinate cannot be generated from Cartesian coordinate
    call vecprod(cellv1(1),cellv1(2),cellv1(3),cellv2(1),cellv2(2),cellv2(3),cellv3(1),cellv3(2),cellv3(3))
    cellv3=cellv3/dsqrt(sum(cellv3**2))
end if
!When PBC is not considered in a direction, PBCnx/y/z should be set to 0 to avoid consider neighbouring cells in corresponding direction
if (ifdoPBCx==0) PBCnx=0
if (ifdoPBCy==0) PBCny=0
if (ifdoPBCz==0) PBCnz=0
end subroutine



!!!------ Convert fractional coordinates to Cartesian coordinates
subroutine fract2Cart(fract,Cart)
use defvar
real*8 fract(3),Cart(3),Amat(3,3),Fcoord(3,1),rcoord(3,1)
Fcoord(:,1)=fract(:)
Amat(:,1)=cellv1(:)
Amat(:,2)=cellv2(:)
Amat(:,3)=cellv3(:)
rcoord=matmul(Amat,Fcoord) !New xyz coordinate
Cart(:)=rcoord(:,1)
end subroutine



!!!------ Convert Cartesian coordinates to fractional coordinates
subroutine Cart2fract(Cart,fract)
use defvar
use util
real*8 Cart(3),fract(3),Amat(3,3),Bmat(3,3),Fcoord(3,1),rcoord(3,1)
rcoord(:,1)=Cart(:)
Amat(:,1)=cellv1(:)
Amat(:,2)=cellv2(:)
Amat(:,3)=cellv3(:)
Bmat=invmat(Amat,3)
Fcoord=matmul(Bmat,rcoord)
fract(:)=Fcoord(:,1)
end subroutine


    
!!!------ If the inputted Cartesian coordinate is out of cell, move it to proper position in cell according to PBC
!See Ian Bush's reply in https://scicomp.stackexchange.com/questions/20165/periodic-boundary-conditions-for-triclinic-box for algorithm
subroutine move_to_cell(xyzin,xyzout)
use defvar
use util
real*8 xyzin(3),xyzout(3),Amat(3,3),Bmat(3,3),Fcoord(3,1),rcoord(3,1)
rcoord(:,1)=xyzin(:)
Amat(:,1)=cellv1(:)
Amat(:,2)=cellv2(:)
Amat(:,3)=cellv3(:)
Bmat=invmat(Amat,3)
!Get fractional coordinate
Fcoord=matmul(Bmat,rcoord)
!Move to central cell
if (ifdoPBCx==1) Fcoord(1,1)=Fcoord(1,1)-floor(Fcoord(1,1))
if (ifdoPBCy==1) Fcoord(2,1)=Fcoord(2,1)-floor(Fcoord(2,1))
if (ifdoPBCz==1) Fcoord(3,1)=Fcoord(3,1)-floor(Fcoord(3,1))
!New xyz coordinte
rcoord=matmul(Amat,Fcoord)
xyzout(:)=rcoord(:,1)
end subroutine


!!!------- The same as move_to_cell, but input x,y,z component, then they will be updated
subroutine move_to_cell_scalar(xin,yin,zin)
real*8 xin,yin,zin,xyzin(3),xyzout(3)
xyzin(1)=xin
xyzin(2)=yin
xyzin(3)=zin
call move_to_cell(xyzin,xyzout)
xin=xyzout(1)
yin=xyzout(2)
zin=xyzout(3)
end subroutine



!!!------ Input X,Y,Z of a point, get index of the cell containing this point
subroutine getpointcell(x,y,z,icell,jcell,kcell)
use defvar
use util
real*8 x,y,z,vec(3),Amat(3,3),Bmat(3,3),Fcoord(3,1),rcoord(3,1)
integer icell,jcell,kcell
if (ifPBC==0) then
    icell=0
    jcell=0
    kcell=0
    return
end if
vec(1)=x
vec(2)=y
vec(3)=z
rcoord(:,1)=vec(:)
Amat(:,1)=cellv1(:)
Amat(:,2)=cellv2(:)
Amat(:,3)=cellv3(:)
Bmat=invmat(Amat,3)
Fcoord=matmul(Bmat,rcoord) !Get fractional coordinate
icell=0;jcell=0;kcell=0
if (ifdoPBCx==1) icell=floor(Fcoord(1,1))
if (ifdoPBCy==1) jcell=floor(Fcoord(2,1))
if (ifdoPBCz==1) kcell=floor(Fcoord(3,1))
end subroutine



!!!-------- Input xyz of points A and B, return coordinate of B closest to point A according to PBC (xyzB2)
subroutine nearest_mirror(xyzA,xyzB,xyzB2)
use defvar
real*8 xyzA(3),xyzB(3),xyztmp(3),xyzB2(3),tvec(3),dist2tmp
dist2min=999999D0
do icell=-ifdoPBCx,ifdoPBCx
    do jcell=-ifdoPBCy,ifdoPBCy
        do kcell=-ifdoPBCz,ifdoPBCz
            call tvec_PBC(icell,jcell,kcell,tvec)
            xyztmp=xyzB+tvec
            dist2tmp=sum((xyzA-xyztmp)**2)
            if (dist2tmp<dist2min) then
                dist2min=dist2tmp
                xyzB2=xyztmp
            end if
        end do
    end do
end do
end subroutine



!!!-------- Input indices of atoms A and B, return their nearest distance, and coordinate of B closest to A
!For isolated system, the B closest to A is still the inputted B
subroutine nearest_atmdistxyz(iatm,jatm,dist,atmx,atmy,atmz)
use defvar
integer iatm,jatm
real*8 dist,atmx,atmy,atmz
real*8 xyzA(3),xyzB(3),xyzB2(3)
xyzA(1)=a(iatm)%x
xyzA(2)=a(iatm)%y
xyzA(3)=a(iatm)%z
xyzB(1)=a(jatm)%x
xyzB(2)=a(jatm)%y
xyzB(3)=a(jatm)%z
if (ifPBC==0) then
    atmx=xyzB(1)
    atmy=xyzB(2)
    atmz=xyzB(3)
    dist=dsqrt(sum((xyzA-xyzB)**2))
else
    call nearest_mirror(xyzA,xyzB,xyzB2)
    atmx=xyzB2(1)
    atmy=xyzB2(2)
    atmz=xyzB2(3)
    dist=dsqrt(sum((xyzA-xyzB2)**2))
end if
end subroutine



!!!-------- Input xyz of points A and B, return their nearest distance according to PBC
subroutine nearest_dist(xyzA,xyzB,dist)
real*8 xyzA(3),xyzB(3),xyzB2(3),dist
call nearest_mirror(xyzA,xyzB,xyzB2)
dist=dsqrt(sum((xyzA-xyzB2)**2))
end subroutine
!!!-------- The same as above, but return vector A->B
subroutine nearest_distvec(xyzA,xyzB,distvec)
real*8 xyzA(3),xyzB(3),xyzB2(3),distvec(3)
call nearest_mirror(xyzA,xyzB,xyzB2)
distvec(:)=xyzB2(:)-xyzA(:)
end subroutine



!!!----------- Input number of translation times in various directions, return translation vector
subroutine tvec_PBC(icell,jcell,kcell,tvec)
use defvar
integer icell,jcell,kcell
real*8 tvec(3)
tvec(1)=icell*cellv1(1)+jcell*cellv2(1)+kcell*cellv3(1)
tvec(2)=icell*cellv1(2)+jcell*cellv2(2)+kcell*cellv3(2)
tvec(3)=icell*cellv1(3)+jcell*cellv2(3)+kcell*cellv3(3)
end subroutine



!!!----- Return coordinate of specific vertex of the cell
!idx is the index of following map
!    5-------8
!   /|      /|
!  6-+-----7 |
!  | |     | |
!  | 1-----+-4
!  |/      |/
!  2-------3
!
!   Z
!   |
!   0---Y    
!  / 
! X
subroutine cellvertex(idx,vertx,verty,vertz)
use defvar
integer idx
real*8 vertx,verty,vertz,v1tmp(3),v2tmp(3),v3tmp(3)
if (idx==1) then
    icell=0;jcell=0;kcell=0
else if (idx==2) then
    icell=1;jcell=0;kcell=0
else if (idx==3) then
    icell=1;jcell=1;kcell=0
else if (idx==4) then
    icell=0;jcell=1;kcell=0
else if (idx==5) then
    icell=0;jcell=0;kcell=1
else if (idx==6) then
    icell=1;jcell=0;kcell=1
else if (idx==7) then
    icell=1;jcell=1;kcell=1
else if (idx==8) then
    icell=0;jcell=1;kcell=1
end if
v1tmp=0
v2tmp=0
v3tmp=0
if (ifdoPBCx==1) v1tmp=cellv1
if (ifdoPBCy==1) v2tmp=cellv2
if (ifdoPBCz==1) v3tmp=cellv3
vertx=icell*v1tmp(1)+jcell*v2tmp(1)+kcell*v3tmp(1)
verty=icell*v1tmp(2)+jcell*v2tmp(2)+kcell*v3tmp(2)
vertz=icell*v1tmp(3)+jcell*v2tmp(3)+kcell*v3tmp(3)
end subroutine

!!----- Wrapper of cellvertex, input two vertex indices return two sets of coordinates
subroutine cellvertex2(idx,jdx,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
integer idx,jdx
real*8 vert1x,vert1y,vert1z,vert2x,vert2y,vert2z
call cellvertex(idx,vert1x,vert1y,vert1z)
call cellvertex(jdx,vert2x,vert2y,vert2z)
end subroutine



!!!------- Return max x,y,z of all cell vertices
subroutine cellmaxxyz(xmax,ymax,zmax)
use defvar
implicit real*8 (a-h,o-z)
real*8 xmax,ymax,zmax
xmax=-999
ymax=-999
zmax=-999
do idx=1,8
    call cellvertex(idx,vertx,verty,vertz)
    if (vertx>xmax) xmax=vertx
    if (verty>ymax) ymax=verty
    if (vertz>zmax) zmax=vertz
end do
end subroutine


!!!------- Return mix x,y,z of all cell vertices
subroutine cellminxyz(xmin,ymin,zmin)
use defvar
implicit real*8 (a-h,o-z)
real*8 xmin,ymin,zmin
xmin=999
ymin=999
zmin=999
do idx=1,8
    call cellvertex(idx,vertx,verty,vertz)
    if (vertx<xmin) xmin=vertx
    if (verty<ymin) ymin=verty
    if (vertz<zmin) zmin=vertz
end do
end subroutine


!!!------- Return size of x,y,z of current box
subroutine cellxyzsize(xsize,ysize,zsize)
implicit real*8 (a-h,o-z)
real*8 xsize,ysize,zsize
call cellminxyz(xmin,ymin,zmin)
call cellmaxxyz(xmax,ymax,zmax)
xsize=xmax-xmin
ysize=ymax-ymin
zsize=zmax-zmin
end subroutine



!!!----- Return cell volume according to cell information
subroutine calc_cellvol(cellvol)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 cellvol,vec(3),mat(3,3)
!V=|a.(bxc)|
call vecprod(cellv2(1),cellv2(2),cellv2(3),cellv3(1),cellv3(2),cellv3(3),vec(1),vec(2),vec(3))
cellvol=abs(sum(cellv1(:)*vec(:)))
!Another way to solve
!mat(:,1)=cellv1(:)
!mat(:,2)=cellv2(:)
!mat(:,3)=cellv3(:)
!cellvol=abs(detmat(mat))
end subroutine



!!---------- Input two vectors, return area     
subroutine calc_vecarea(vec1,vec2,area)
implicit real*8 (a-h,o-z)
real*8 vec1(3),vec2(3),area
vec2norm=dsqrt(sum(vec2**2))
tmp=sum(vec1*(vec2/vec2norm)) !Length of vec1 along vector 2
tmp2=dsqrt(sum(vec1**2)-tmp**2) !Length of vec1 perpendicular to vector 2
area=tmp2*vec2norm
end subroutine



!!!------ Calculate cell vectors cellv1/2/3 based on inputted a,b,c (in Bohr) and alpha,beta,gamma (in degree), reference: http://gisaxs.com/index.php/Unit_cell
subroutine abc2cellv(asize,bsize,csize,alpha,beta,gamma)
use defvar
implicit real*8 (a-h,o-z)
real*8 asize,bsize,csize,alpha,beta,gamma
alpha=alpha/180*pi
beta=beta/180*pi
gamma=gamma/180*pi
cellv1(1)=asize
cellv1(2)=0
cellv1(3)=0
cellv2(1)=bsize*cos(gamma)
cellv2(2)=bsize*sin(gamma)
cellv2(3)=0
cellv3(1)=csize*cos(beta)
cellv3(2)=csize*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
cellv3(3)=csize*dsqrt( 1 - cos(beta)**2 - ((cos(alpha)-cos(beta)*cos(gamma))/sin(gamma))**2 )
end subroutine



!!!------ Return a,b,c (in Angstrom!) and alpha,beta,gamma based on current cell information
subroutine getcellabc(asize,bsize,csize,alpha,beta,gamma)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 asize,bsize,csize,alpha,beta,gamma
asize=dsqrt(sum(cellv1**2))*b2a
bsize=dsqrt(sum(cellv2**2))*b2a
csize=dsqrt(sum(cellv3**2))*b2a
alpha=vecang(cellv2(1),cellv2(2),cellv2(3),cellv3(1),cellv3(2),cellv3(3))
beta =vecang(cellv1(1),cellv1(2),cellv1(3),cellv3(1),cellv3(2),cellv3(3))
gamma=vecang(cellv1(1),cellv1(2),cellv1(3),cellv2(1),cellv2(2),cellv2(3))
end subroutine



!!!------- Get reciprocal vectors
subroutine getrcellvec(rcellv1,rcellv2,rcellv3)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 rcellv1(3),rcellv2(3),rcellv3(3)
call calc_cellvol(cellvol)
call vecprod(cellv2(1),cellv2(2),cellv2(3),cellv3(1),cellv3(2),cellv3(3),rcellv1(1),rcellv1(2),rcellv1(3))
call vecprod(cellv1(1),cellv1(2),cellv1(3),cellv3(1),cellv3(2),cellv3(3),rcellv2(1),rcellv2(2),rcellv2(3))
call vecprod(cellv1(1),cellv1(2),cellv1(3),cellv2(1),cellv2(2),cellv2(3),rcellv3(1),rcellv3(2),rcellv3(3))
rcellv1=rcellv1*2*pi/cellvol
rcellv2=rcellv2*2*pi/cellvol
rcellv3=rcellv3*2*pi/cellvol
end subroutine



!!!--------- Make molecule truncated by box whole, called by subroutine geom_operation
subroutine makemolwhole
use defvar
use util
implicit real*8 (a-h,o-z)
integer fraglist(ncenter,ncenter) !fraglist(1:frgnatm(ifrag),ifrag) is indices of atoms in fragment ifrag
integer fraglist_old(ncenter,ncenter)
integer frgnatm(ncenter)
integer nfrag !Total number of fragments
integer iffrag(ncenter),intarr(ncenter)
real*8 tvec(3)
write(*,*)
write(*,*) "Please wait..."
!Construct fragments
intarr=1
iatm=1
nfrag=0
frgnatm(:)=0
do while(.true.)
    call getfragatoms(iatm,iffrag)
    nfrag=nfrag+1
    do jatm=1,ncenter
        if (iffrag(jatm)==1) then
            frgnatm(nfrag)=frgnatm(nfrag)+1
            fraglist(frgnatm(nfrag),nfrag)=jatm
            intarr(jatm)=0
        end if
    end do
    if (all(intarr==0)) exit
    do iatm=1,ncenter
        if (intarr(iatm)==1) exit
    end do
end do

!Sort fragments from small to large number of members
forall(i=1:nfrag) intarr(i)=i
call sorti4(frgnatm(1:nfrag),list=intarr)
fraglist_old=fraglist
do ifrag=1,nfrag
    fraglist(:,ifrag)=fraglist_old(:,intarr(ifrag))
end do
!do ifrag=1,nfrag
!    write(*,"(' Fragment',i5,':')") ifrag
!    write(*,"(15i5)") fraglist(1:frgnatm(ifrag),ifrag)
!end do

!Move fragments if new bond can be detected to other atoms
ncombine=0
combinefrag: do while(.true.)
    do ifrag=1,nfrag
        if (frgnatm(ifrag)==0) cycle !This fragment has been combined to another
        !write(*,"(' Testing fragment',i5)") ifrag
        imove=0
        testfrag: do icell=-1,1
            do jcell=-1,1
                do kcell=-1,1
                    if (icell==0.and.jcell==0.and.kcell==0) cycle
                    call tvec_PBC(icell,jcell,kcell,tvec)
                    do idx=1,frgnatm(ifrag)
                        iatm=fraglist(idx,ifrag)
                        do jfrag=1,nfrag
                            if (jfrag==ifrag) cycle
                            if (frgnatm(jfrag)==0) cycle !This fragment has been combined to another
                            do jdx=1,frgnatm(jfrag)
                                jatm=fraglist(jdx,jfrag)
                                dist=dsqrt((a(iatm)%x+tvec(1)-a(jatm)%x)**2+(a(iatm)%y+tvec(2)-a(jatm)%y)**2+(a(iatm)%z+tvec(3)-a(jatm)%z)**2)
                                if (dist < bondcrit*(covr(a(iatm)%index)+covr(a(jatm)%index))) then
                                    imove=1
                                    exit testfrag
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do testfrag
        if (imove==1) then
            !write(*,"(' Moving fragment',i5,' by X/Y/Z:',3f10.4,' Bohr')") ifrag,tvec
            a(fraglist(1:frgnatm(ifrag),ifrag))%x=a(fraglist(1:frgnatm(ifrag),ifrag))%x+tvec(1)
            a(fraglist(1:frgnatm(ifrag),ifrag))%y=a(fraglist(1:frgnatm(ifrag),ifrag))%y+tvec(2)
            a(fraglist(1:frgnatm(ifrag),ifrag))%z=a(fraglist(1:frgnatm(ifrag),ifrag))%z+tvec(3)
            !write(*,"(' Combining fragment',i5,' and',i5)") ifrag,jfrag
            fraglist( frgnatm(jfrag)+1:frgnatm(jfrag)+frgnatm(ifrag) , jfrag ) = fraglist(1:frgnatm(ifrag),ifrag)
            frgnatm(jfrag)=frgnatm(jfrag)+frgnatm(ifrag)
            frgnatm(ifrag)=0
            ncombine=ncombine+1
            cycle combinefrag
        end if
    end do
    if (ifrag==nfrag+1) exit
end do combinefrag
if (allocated(connmat)) deallocate(connmat) !getfragatoms generated connmat, while after changing order, the relationship is no longer valid
write(*,"(a,i7)") " Done! Number of combinations performed:",ncombine
end subroutine




!!!--------- Construct global array a_tmp, which contains real atoms and replicated boundary atoms
!ncenter_tmp is returned variable containing number of all atoms including the replicated boundary ones
subroutine construct_atmp_withbound(ncenter_tmp)
use defvar
implicit real*8 (a-h,o-z)
real*8 Cart(3),fract(3),ftest(3)
integer ncenter_tmp

call getcellabc(asize,bsize,csize,alpha,beta,gamma)
devthres=1D-2 !If distance to cell boundary is less than 0.01 Bohr, then it will be regarded as boundary atom
fdeva=devthres/asize
fdevb=devthres/bsize
fdevc=devthres/csize

!itime=1: Test how many atoms to be added, =2: Actual constract a_tmp
do itime=1,2
    nadd=0
    do iatm=1,ncenter
        Cart(1)=a(iatm)%x
        Cart(2)=a(iatm)%y
        Cart(3)=a(iatm)%z
        call Cart2fract(Cart,fract)
        xdiff=asize*abs(fract(1)-nint(fract(1)))
        ydiff=bsize*abs(fract(2)-nint(fract(2)))
        zdiff=csize*abs(fract(3)-nint(fract(3)))
        if (xdiff<devthres.or.ydiff<devthres.or.zdiff<devthres) then !Existing boundary atom is at least very close to one of cell walls
            !Try to replicate the existing boundary atom to all possible neighbouring mirror sites
            do ix=-1,1
                do iy=-1,1
                    do iz=-1,1
                        if (ix==0.and.iy==0.and.iz==0) cycle
                        ftest(1)=fract(1)+ix
                        ftest(2)=fract(2)+iy
                        ftest(3)=fract(3)+iz
                        !Mirror boundary must be within or quasi within current cell
                        if (ftest(1)>-fdeva.and.ftest(1)<1+fdeva.and.ftest(2)>-fdevb.and.ftest(2)<1+fdevb.and.ftest(3)>-fdevc.and.ftest(3)<1+fdevc) then
                            call fract2Cart(ftest,Cart)
                            !Check if the mirror boundary atom is too close to existing atoms
                            iadd=1
                            do jatm=1,ncenter
                                if (jatm==iatm) cycle
                                dist=dsqrt( (Cart(1)-a(jatm)%x)**2 + (Cart(2)-a(jatm)%y)**2 + (Cart(3)-a(jatm)%z)**2 )
                                if (dist<devthres) then !Too close, skip it
                                    iadd=0
                                    exit
                                end if
                            end do
                            if (iadd==1) then
                                nadd=nadd+1
                                if (itime==2) then
                                    a_tmp(ncenter+nadd)=a(iatm)
                                    a_tmp(ncenter+nadd)%x=Cart(1)
                                    a_tmp(ncenter+nadd)%y=Cart(2)
                                    a_tmp(ncenter+nadd)%z=Cart(3)
                                    a_tmp_idx(ncenter+nadd)=iatm
                                end if
                            end if
                        end if
                    end do
                end do
            end do
        end if
    end do
    if (itime==1) then
        ncenter_tmp=ncenter+nadd
        if (allocated(a_tmp)) deallocate(a_tmp)
        if (allocated(a_tmp_idx)) deallocate(a_tmp_idx)
        allocate(a_tmp(ncenter_tmp),a_tmp_idx(ncenter_tmp))
        a_tmp(1:ncenter)=a(1:ncenter)
        forall(i=1:ncenter) a_tmp_idx(i)=i
    end if
end do
end subroutine



!!!------- Check if a point is close enough to any cell plane
!pos: XYZ coordinate
!thres: Threshold of determination in Bohr
!ionplane=1: On plane, =0: Not
subroutine check_on_plane(pos,thres,ionplane)
use defvar
implicit real*8 (a-h,o-z)
real*8 pos(3),fract(3),thres
integer ionplane

call getcellabc(asize,bsize,csize,alpha,beta,gamma)
call Cart2fract(pos,fract)
xdiff=asize*abs(fract(1)-nint(fract(1)))
ydiff=bsize*abs(fract(2)-nint(fract(2)))
zdiff=csize*abs(fract(3)-nint(fract(3)))
ionplane=0
if (xdiff<thres.or.ydiff<thres.or.zdiff<thres) ionplane=1
end subroutine



!!!-------- Back up current PBC information to memory
subroutine savePBCinfo
use defvar
ifPBC_bk=ifPBC
cellv1_bk=cellv1
cellv2_bk=cellv2
cellv3_bk=cellv3
end subroutine



!!!-------- Restore backed up PBC information
subroutine loadPBCinfo
use defvar
ifPBC=ifPBC_bk
cellv1=cellv1_bk
cellv2=cellv2_bk
cellv3=cellv3_bk
call init_PBC
end subroutine



!!!-------- Convert index of a grid outside the cell into intracell index
!Note that grid index starts from 1, so 0 in X direction should be treated as nx
subroutine PBCgrididx(ix,iy,iz)
use defvar
integer ix,iy,iz
ix=ix-floor(float(ix-1)/nx)*nx
iy=iy-floor(float(iy-1)/ny)*ny
iz=iz-floor(float(iz-1)/nz)*nz
end subroutine
!!!-------- Wrap index of a grid in direction 1 into intracell index
subroutine PBCgrididx1(ix)
use defvar
integer ix
ix=ix-floor(float(ix-1)/nx)*nx
end subroutine
!!!-------- Wrap index of a grid in direction 2 into intracell index
subroutine PBCgrididx2(iy)
use defvar
integer iy
iy=iy-floor(float(iy-1)/ny)*ny
end subroutine
!!!-------- Wrap index of a grid in direction 3 into intracell index
subroutine PBCgrididx3(iz)
use defvar
integer iz
iz=iz-floor(float(iz-1)/nz)*nz
end subroutine



!!!-------- Load cell information from [Cell].txt in current folder if available
!The format of this file should be identical to content below [Cell] field of .molden file
subroutine loadcellinfo_txt
use defvar
implicit real*8 (a-h,o-z)
character selectyn,c80tmp*80

ifileid=11
inquire(file="[Cell].txt",exist=alive)
if (alive) then
    write(*,"(a)") " Cell information is not defined in input file, but [Cell].txt is found in current folder, do you want to load cell information from it? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
        open(ifileid,file="[Cell].txt",status="old")
        read(ifileid,"(a)") c80tmp
        read(c80tmp,*,iostat=ierror) alen,blen,clen,anga,angb,angc
        if (ierror==0) then !Load as cell parameters
            call abc2cellv(alen/b2a,blen/b2a,clen/b2a,anga,angb,angc)
            ifPBC=3
        else !Load as cell vectors
            read(c80tmp,*) cellv1
            ifPBC=1
            read(ifileid,*,iostat=ierror) cellv2
            if (ierror/=0) then
                cellv2=0
            else
                ifPBC=ifPBC+1
                read(ifileid,*,iostat=ierror) cellv3
                if (ierror/=0) then
                    cellv3=0
                else
                    ifPBC=ifPBC+1
                end if
            end if
            cellv1=cellv1/b2a
            cellv2=cellv2/b2a
            cellv3=cellv3/b2a
        end if
        close(ifileid)
    end if
end if
end subroutine




!!!-------- Calculate spacing between cell planes
!h,k,l: Miller index
!d: Returned spacing in Angstrom
!Ref: Triclinic form in https://en.wikipedia.org/wiki/Crystal_structure
subroutine cellplane_spacing(h,k,l,d)
use defvar
use util
implicit real*8 (a-h,o-z)
integer h,k,l
real*8 d

call getcellabc(asize,bsize,csize,alpha,beta,gamma)
arad=ang2rad(alpha)
brad=ang2rad(beta)
grad=ang2rad(gamma)
sina=sin(arad)
sinb=sin(brad)
sing=sin(grad)
cosa=cos(arad)
cosb=cos(brad)
cosg=cos(grad)

tmp1=h**2/asize**2*sina**2 + k**2/bsize**2*sinb**2 + l**2/csize**2*sing**2
tmp2=2*k*l/bsize/csize*(cosb*cosg-cosa) + 2*h*l/asize/csize*(cosg*cosa-cosb) + 2*h*k/asize/bsize*(cosa*cosb-cosg)
tmp3=1-cosa**2-cosb**2-cosg**2+2*cosa*cosb*cosg

d=(tmp1+tmp2)/tmp3
d=1D0/dsqrt(d)
end subroutine