module plot
use defvar
use dislin
use util
implicit real*8 (a-h,o-z)

contains


!!!------------- Draw molecular structure and show isosurface of orbitals or grid data
!Also draw topology analysis results, surface extrema, basin or domain space etc.
subroutine drawmol
use topo
use surfvertex
use basinintmod
implicit real*8 (a-h,o-z)
integer i,j,idxtmp,iret,screenx,screeny,ipath,ipathp1,ipt,icp1,icp2,ipathtype,ipathmidpt,isurf,interval
real*8 abslenx,absleny,abslenz,plotlenx,plotleny,plotlenz !absolute and real 3D coordinate
real*8 plot2abs,xplotcoor,yplotcoor,absx,absy,absz,dist,textheighmod,extsiz,augplotlen
real*8 trianglex(3),triangley(3),trianglez(3)
real*8 arrayx(nx),arrayy(ny),arrayz(nz)
character ctemp*5,c80tmp*80

!Set viewpoint
!Note that due to limitation of DISLIN, it is not possible to view molecule in all viewpoints. The YVU should be limited to between -90 and 90, else the viewpoint will suddently flip
XVUold=XVU
YVUold=YVU
if (YVU==90) YVU=89.999D0 !Temporarily modify YVU, otherwise when YVU equals to 90 or -90, the perspective will jump suddenly
if (YVU==-90) YVU=-89.999D0

!Determine x/y/zlow, x/y/zhigh and plotlenx/y/z, which are lower, upper positions and length of the axis
if (ifPBC>0) then !For PBC case, make the axes large enough to show all atoms
    call cellminxyz(xmin,ymin,zmin)
    xlow=min(xmin-3,minval(a%x)-3)
    ylow=min(ymin-3,minval(a%y)-3)
    zlow=min(zmin-3,minval(a%z)-3)
    call cellmaxxyz(xmax,ymax,zmax)
    xhigh=max(xmax+3,maxval(a%x)+3)
    yhigh=max(ymax+3,maxval(a%y)+3)
    zhigh=max(zmax+3,maxval(a%z)+3)
    plotlenx=xhigh-xlow
    plotleny=yhigh-ylow
    plotlenz=zhigh-zlow
else if (ifiletype==7.or.ifiletype==8) then !For cub/grd/cube/vti file, make the axis able to show entire grid data
    call gridminxyz(xmin,ymin,zmin)
    call gridmaxxyz(xmax,ymax,zmax)
	xlow=xmin-1
	ylow=ymin-1
	zlow=zmin-1
	xhigh=xmax+1
	yhigh=ymax+1
	zhigh=zmax+1
    plotlenx=xhigh-xlow
    plotleny=yhigh-ylow
    plotlenz=zhigh-zlow
else !Other cases, determine displayed spatial scope by boundary atoms of the system
    !Determine axis range
    augplotlen=8D0 !Common augmentation of plotting axis length (the sum of both direction)
    if (GUI_mode==4.and.ifunctopo==25) augplotlen=24D0 !Minima of vdW potential are often very far from atoms, use larger augmentation
    if (GUI_mode==5) augplotlen=12D0 !Because surface extreme points laying on vdW surface, the scatter region is large, so use larger augmentation
    if (GUI_mode==7) augplotlen=12D0 !Using GUI to set box, the case is complicated, so use larger augmentation
    if ((idrawisosur==1.and.aug3D>4).or.GUI_mode==6) augplotlen=aug3D*2.2D0 !Shouldn't be 2.0 as expected, otherwise sometimes there will be a band occuring at boundary

    !Below comment and code are redundant, ignore them...
    !Commonly, augment distances in all side are idential. However, for very long chain system, the molecule cannot be displayed unless ZVU is first set to
    !fairly large value, namely extremely zooming out. To circumvent this problem, we first calculate geometric average in X,Y,Z (Davg)
    !if ratio between molecular size in any direction and Davg is by far lower than a threshold, then we make all axis lengths identical to the longest one,
    !in this case the molecule must be able to be shown
    xmolsize=maxval(a%x)-minval(a%x)
    ymolsize=maxval(a%y)-minval(a%y)
    zmolsize=maxval(a%z)-minval(a%z)
    !avgD=(xmolsize*ymolsize*zmolsize)**(1D0/3D0)
    !ratiox=xmolsize/avgD
    !ratioy=ymolsize/avgD
    !ratioz=zmolsize/avgD
    !thresratio=0.2D0
    !write(*,*) xmolsize, ymolsize, zmolsize,augplotlen
    !if (xmolsize<thresratio .or. ymolsize<thresratio .or. zmolsize<thresratio) then
    !    rlong=max(max(xmolsize,ymolsize),zmolsize)
    !    augplotlenx=augplotlen+(rlong-xmolsize)/2
    !    augplotleny=augplotlen+(rlong-ymolsize)/2
    !    augplotlenz=augplotlen+(rlong-zmolsize)/2
    !else
        augplotlenx=augplotlen
        augplotleny=augplotlen
        augplotlenz=augplotlen
    !end if

	xlow=minval(a%x)-augplotlenx/2
	ylow=minval(a%y)-augplotleny/2
	zlow=minval(a%z)-augplotlenz/2
	plotlenx=xmolsize+augplotlenx
	plotleny=ymolsize+augplotleny
	plotlenz=zmolsize+augplotlenz
    xhigh=xlow+plotlenx
    yhigh=ylow+plotleny
    zhigh=zlow+plotlenz
end if

!Initialize DISLIN
abslenx=2D0 !Absolute length in DISLIN
absleny=abslenx*plotleny/plotlenx
abslenz=abslenx*plotlenz/plotlenx
plot2abs=abslenx/plotlenx !The relationship between molecular coordinate and absolute coordinate
if (isavepic==0) then
	call METAFL('CONS')
	if (GUI_mode==6) call METAFL('CONS') !Namely showing basin, using opengl by default to accelerate displaying, however when savepic, if still use opengl, things cannot be properly shown
	CALL setxid(idisgraph,'WIDGET')
else if (isavepic==1) then
    if (iorbvis==0) then
        call setfil("dislin."//trim(graphformat))
    else
        write(c80tmp,"(i6.6,'.',a)") iorbvis,trim(graphformat)
        call setfil(trim(c80tmp))
    end if
	CALL setxid(0,'NONE')
	call METAFL(graphformat)
	call winsiz(graph3Dwidth,graph3Dheight) !Actual image size is set by this routine, namely 770*770
end if
CALL SCRMOD('revers')
CALL PAGE(3000,3000)
CALL IMGFMT("RGB") !If don't call this routine, the saved picture use palette color mode by default, the color is not smooth
CALL DISINI
CALL ORIGIN(ORIGIN_3D_X,ORIGIN_3D_Y) !Set origin of the map
if (iorthoview==0) then
	CALL PROJ3D("PERSPECTIVE")
	call VANG3D(VANG3DANG)
else if (iorthoview==1) then
	CALL PROJ3D("ORTHO")
end if
CALL VFOC3D(XFOC,YFOC,ZFOC,"ABS")
CALL VUP3D(camrotang)
call ERRMOD("ALL","OFF")
! call linmod("ON","SMOOTH") !It seems that Anti-aliased doesn't take effect
CALL LABDIG(1,"X")
CALL LABDIG(1,"Y")
CALL LABDIG(1,"Z")
if (ishowaxis==0) call NOGRAF
CALL VIEW3D(XVU,YVU,ZVU,"ANGLE")
if (iorthoview==1) call vscl3d(XFAC)
CALL erase

!Set font
if (isavepic==0.or.graphformat=="pdf ") then
	CALL HWFONT
else if (isavepic==1) then
	if (ttfontfile=="none") then
		CALL HELVES
    else
		CALL TTFONT(ttfontfile)
    end if
    CALL SHDCHA
end if
if (isavepic==0.and.isys==1) then
	call height(60)
	CALL HNAME(60)
else
	call height(50)
	CALL HNAME(50)
end if

!Set axis
CALL NAME('X-axis (Bohr)','X')
CALL NAME('Y-axis (Bohr)','Y')
CALL NAME('Z-axis (Bohr)','Z')
CALL NAMDIS(100,'XYZ')
call axis3D(abslenx,absleny,abslenz)
CALL AXSPOS(400,2600)
CALL AX3LEN(2600,2600,2600)
call labl3D("horizontal")
CALL LIGHT('on')
CALL SHDMOD('SMOOTH','SURFACE') !By default, each surface triangle has the same color and thus doesn't look smooth
nsclx=8
spcx=ceiling(plotlenx/nsclx) !Expected step size in X
spcy=ceiling(plotleny/ (plotleny/plotlenx*nsclx) )
spcz=ceiling(plotlenz/ (plotlenz/plotlenx*nsclx) /1.5D0 )
shiftx=mod(xlow,spcx) !Shift the initial point of the axis to make a label just occur at 0.0
shifty=mod(ylow,spcy)
shiftz=mod(zlow,spcz)
! CALL FLAB3D !If use this, the starting label of Y and Z axis will be plotted, however this may suppress the starting label of X axis
CALL GRAF3D(xlow,xhigh,xlow-shiftx,spcx, ylow,yhigh,ylow-shifty,spcy, zlow,zhigh,zlow-shiftz,spcz)
if (ishowaxis==1) CALL GRID3D (2,2,'bottom')
CALL ZBFINI(IRET) !Enable Z-buffer to determine visibility
call litpos(1,XVU,YVU,ZVU,'ANGLE')
call litmod(1,'on') !Dislin default light 1 is on, and off for others

!Construct global array a_tmp, for non-PBC system it is the same as a, while for PBC case it contains real atoms and replicated boundary atoms
if (idrawmol==1.or.ishowatmlab==1) then
	if (ishowboundaryatom==1.and.ifPBC>0) then
		call construct_atmp_withbound(ncenter_tmp)
	else
		if (allocated(a_tmp)) deallocate(a_tmp)
		ncenter_tmp=ncenter
		allocate(a_tmp(ncenter),a_tmp_idx(ncenter))
		a_tmp=a
        forall(i=1:ncenter) a_tmp_idx(i)=i
	end if
end if

if (idrawmol==1) then
	do ilight=2,8
		call litmod(ilight,'off')
	end do
	!Draw atoms
	do i=1,ncenter_tmp
        if (ishowhydrogen==0.and.a_tmp(i)%index==1) cycle
		CALL MATOP3(atm3Dclr(a_tmp(i)%index,1),atm3Dclr(a_tmp(i)%index,2),atm3Dclr(a_tmp(i)%index,3),'diffuse')
        if (allocated(highlightatomlist)) then
            if (any(highlightatomlist==i)) CALL MATOP3(0.3D0,1D0,1D0,'diffuse')
        end if
		if (a_tmp(i)%index/=0) then !Normal atoms
			if (isavepic==1) then
                CALL SPHE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,vdwr(a_tmp(i)%index)/4*ratioatmsphere,40,40)
            else !Resolution of 30,30 is nearly perfect. For large system lower the resolution for faster view
                if (ncenter_tmp<100) then
                    CALL SPHE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,vdwr(a_tmp(i)%index)/4*ratioatmsphere,30,30)
                else if (ncenter_tmp<300) then
                    CALL SPHE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,vdwr(a_tmp(i)%index)/4*ratioatmsphere,20,20)
                else
                    CALL SPHE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,vdwr(a_tmp(i)%index)/4*ratioatmsphere,10,10)
                end if
            end if
		else !Bq atoms. The size keeps fixed
			CALL SPHE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,vdwr(a_tmp(i)%index)/4,20,20)
		end if
	end do
	!Draw bonds. If connectivity is available, then do not automatically determine bonding
	CALL MATOP3(bondclrR,bondclrG,bondclrB,'diffuse')
    if (allocated(connmat)) then !Note that when using defined connectivity, bonds for boundary atoms will not be plotted
	    do i=1,ncenter
		    do j=i+1,ncenter
                if (ishowhydrogen==0.and.(a(i)%index==1.or.a(j)%index==1)) cycle
			    if (a(i)%index==0.or.a(j)%index==0) cycle !Never make Bq atoms bonding
			    if (connmat(i,j)/=0) then
                    if (ncenter<100.or.isavepic==1) then
                        CALL TUBE3D(a(i)%x,a(i)%y,a(i)%z,a(j)%x,a(j)%y,a(j)%z,bondradius,20,20) !Resolution of 20,20 is visually perfect
                    else if (ncenter<300) then
                        CALL TUBE3D(a(i)%x,a(i)%y,a(i)%z,a(j)%x,a(j)%y,a(j)%z,bondradius,15,15)
                    else
                        CALL TUBE3D(a(i)%x,a(i)%y,a(i)%z,a(j)%x,a(j)%y,a(j)%z,bondradius,8,8)
                    end if
                end if
		    end do
	    end do
    else
	    do i=1,ncenter_tmp
		    do j=i+1,ncenter_tmp
                if (ishowhydrogen==0.and.(a_tmp(i)%index==1.or.a_tmp(j)%index==1)) cycle
			    if (a_tmp(i)%index==0.or.a_tmp(j)%index==0) cycle !Never make bonding
			    dist=dsqrt( (a_tmp(i)%x-a_tmp(j)%x)**2+(a_tmp(i)%y-a_tmp(j)%y)**2+(a_tmp(i)%z-a_tmp(j)%z)**2 )
			    !If the distance between two atoms exceeds 15% of sum of their covalent radii, they will be seemed to be not bonded
			    if (dist<( covr(a_tmp(i)%index)+covr(a_tmp(j)%index) )*bondcrit) then
                    if (ncenter_tmp<100.or.isavepic==1) then
                        CALL TUBE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,a_tmp(j)%x,a_tmp(j)%y,a_tmp(j)%z,bondradius,20,20) !Resolution of 20,20 is visually perfect
                    else if (ncenter_tmp<300) then
                        CALL TUBE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,a_tmp(j)%x,a_tmp(j)%y,a_tmp(j)%z,bondradius,15,15)
                    else
                        CALL TUBE3D(a_tmp(i)%x,a_tmp(i)%y,a_tmp(i)%z,a_tmp(j)%x,a_tmp(j)%y,a_tmp(j)%z,bondradius,8,8)
                    end if
                end if
		    end do
	    end do
    end if
end if

!Draw critical points
if (ishow3n3==1.or.ishow3n1==1.or.ishow3p1==1.or.ishow3p3==1) then
	numcp_tmp=numcp
	CPpos_tmp=CPpos
    CPtype_tmp=CPtype
    forall(i=1:numCP) CP_tmp_idx(i)=i
	if (ishowboundarytopo==1.and.ifPBC>0) call construct_CPtmp_withbound(numcp_tmp)
end if
if (ishow3n3==1) then
	CALL MATOP3(CP3n3RGB(1),CP3n3RGB(2),CP3n3RGB(3),'diffuse') !Purple
	do i=1,numcp_tmp
		if (CPtype_tmp(i)==1) CALL SPHE3D(CPpos_tmp(1,i),CPpos_tmp(2,i),CPpos_tmp(3,i),0.15D0*ratioCPsphere,20,20)
	end do
end if
if (ishow3n1==1) then
	CALL MATOP3(CP3n1RGB(1),CP3n1RGB(2),CP3n1RGB(3),'diffuse') !Orange
	do i=1,numcp_tmp
		if (CPtype_tmp(i)==2) CALL SPHE3D(CPpos_tmp(1,i),CPpos_tmp(2,i),CPpos_tmp(3,i),0.1D0*ratioCPsphere,20,20)
	end do
end if
if (ishow3p1==1) then
	CALL MATOP3(CP3p1RGB(1),CP3p1RGB(2),CP3p1RGB(3),'diffuse') !Yellow
	do i=1,numcp_tmp
		if (CPtype_tmp(i)==3) CALL SPHE3D(CPpos_tmp(1,i),CPpos_tmp(2,i),CPpos_tmp(3,i),0.1D0*ratioCPsphere,20,20)
	end do
end if
if (ishow3p3==1) then
	CALL MATOP3(CP3p3RGB(1),CP3p3RGB(2),CP3p3RGB(3),'diffuse') !Green
	do i=1,numcp_tmp
		if (CPtype_tmp(i)==4) CALL SPHE3D(CPpos_tmp(1,i),CPpos_tmp(2,i),CPpos_tmp(3,i),0.1D0*ratioCPsphere,20,20)
	end do
end if

!Show trace of movement to attractor in basin generating process. If you want to enable this, move declaration of "ntrjgrid" and "trjgrid" from "generatebasin" to basinintmod
! CALL MATOP3(0.1D0, 0.5D0, 0.85D0, 'diffuse')
! do igrd=2,ntrjgrid
! 	attx=orgx+(trjgrid(igrd,1)-1)*dx
! 	atty=orgy+(trjgrid(igrd,2)-1)*dy
! 	attz=orgz+(trjgrid(igrd,3)-1)*dz
! 	CALL SPHE3D(attx,atty,attz,0.05D0,20,20)
! 	call tube3D(attx,atty,attz,orgx+(trjgrid(igrd-1,1)-1)*dx,orgy+(trjgrid(igrd-1,2)-1)*dy,orgz+(trjgrid(igrd-1,3)-1)*dx,0.015D0,30,30)
! end do

!For basin integration analysis
if (numatt>0.and.ishowatt==1) then
	!Draw attractors
	do iatt=1,numatt
		if (attval(iatt)>=0) then
			CALL MATOP3(0.65D0, 0.9D0, 0.45D0, 'diffuse')
		else
			CALL MATOP3(0.35D0, 0.7D0, 0.9D0, 'diffuse')
		end if
		CALL SPHE3D(attxyz(1,iatt),attxyz(2,iatt),attxyz(3,iatt),attsphsize,20,20)
	end do
	!Draw basin
	if (idrawbasinidx/=-10) then !-10 means don't draw basins
		if (idrawbasinidx<=0) then !Positive basin or the basins consist of boundary or unassigned grids
			CALL MATOP3(0D0, 0.8D0, 0D0, 'diffuse') !Green
		else if (realattval(idrawbasinidx)>=0) then
			CALL MATOP3(0D0, 0.8D0, 0D0, 'diffuse') !Green
		else !Negative basins
			CALL MATOP3(0.3D0, 0.6D0, 1D0, 'diffuse') !Light blue
		end if
        if (basinsphsize==0) then !Sphere size has not been initialized
			call calc_dvol(dvol)
            !Convert grid volume to effective radius: dvol=3/4*pi*r^3, so r=(dvol*4/3/pi)**(1D0/3D0)
            basinsphsize=(dvol*4/3/pi)**(1D0/3D0)*1.2D0 !1.2 is used to slightly enhance to fully fill the grid
			!basinsphsize=dsqrt(dx**2+dy**2+dz**2)/2D0 !Not suitable for non-orthogonal grid
        end if
		do iz=2,nz-1
			do iy=2,ny-1
				do ix=2,nx-1
					if (gridbas(ix,iy,iz)==idrawbasinidx) then
						call getgridxyz(ix,iy,iz,tmpx,tmpy,tmpz)
                        if (ishowbasinmethod==1) then !Entire basin
						    if (interbasgrid(ix,iy,iz)) CALL SPHE3D(tmpx,tmpy,tmpz,basinsphsize,4,4) !Show sphere of interbasin grids
						    if (idrawinternalbasin==1) then !Draw quad at box boundary to shield internal region of basin
							    if (ix==2.or.ix==nx-1.or.iy==2.or.iy==ny-1.or.iz==2.or.iz==nz-1) CALL QUAD3D(tmpx,tmpy,tmpz,dx,dy,dz)
						    end if
                        else if (ishowbasinmethod==2) then !Basin within vdW surface
                            if (rhocub(ix,iy,iz)<0.001D0) cycle
                            if (interbasgrid(ix,iy,iz)) then !Interbasin grid
                                CALL SPHE3D(tmpx,tmpy,tmpz,basinsphsize,4,4)
                            else !Internal grid in the basin
                                if (idrawinternalbasin==1) CALL SPHE3D(tmpx,tmpy,tmpz,basinsphsize,4,4)
                            end if
                        end if
					end if
				end do
			end do
		end do
	end if
end if

!Draw domain defined by isosurface as grids
if (idrawdomain==1.and.idrawdomainidx/=0) then
	CALL MATOP3(0D0, 0.8D0, 0D0, 'diffuse') !Green
	tmpsphrad=dsqrt(dx**2+dy**2+dz**2)/2D0
	do igrd=1,domainsize(idrawdomainidx)
		idx=domaingrid(igrd,idrawdomainidx)
		CALL SPHE3D(gridxyz(1,idx),gridxyz(2,idx),gridxyz(3,idx),tmpsphrad,4,4)
	end do
end if

!Draw topology paths
if (idrawpath==1) then
	numpath_tmp=numpath
	pathnumpt_tmp=pathnumpt
    topopath_tmp=topopath
    forall(i=1:numpath) path_tmp_idx(i)=i
	if (ishowboundarytopo==1.and.ifPBC>0) call construct_pathtmp_withbound(numpath_tmp) !Generate paths at cell boundary
    ntuberes=5
    if (numpath_tmp>200) ntuberes=3 !For showing faster, use lower resolution
	do ipath=1,numpath_tmp
		call path_cp(path_tmp_idx(ipath),icp1,icp2,ipathtype)
		!ipathtype=0: other   =1: (3,-1)->(3,-3) =2: (3,+1)->(3,+3) =3: (3,-1)<-->(3,+1)
		if (ipathtype==0) call MATOP3(0.8D0, 0.8D0, 0.8D0, 'diffuse')
		if (ipathtype==1) call MATOP3(0.9D0, 0.7D0, 0.0D0, 'diffuse')
		if (ipathtype==2) call MATOP3(0.2D0, 0.7D0, 0.2D0, 'diffuse')
		if (ipathtype==3) call MATOP3(0.8D0, 0.8D0, 0.3D0, 'diffuse')
		do ipt=2,pathnumpt_tmp(ipath)
            xlast=topopath_tmp(1,ipt-1,ipath)
            ylast=topopath_tmp(2,ipt-1,ipath)
            zlast=topopath_tmp(3,ipt-1,ipath)
            xthis=topopath_tmp(1,ipt,ipath)
            ythis=topopath_tmp(2,ipt,ipath)
            zthis=topopath_tmp(3,ipt,ipath)
            dist=dsqrt((xlast-xthis)**2+(ylast-ythis)**2+(zlast-zthis)**2)
            !If distance between two path points is very large, then the two points must cross box, so do not draw
            if (dist<pathstepsize*2) CALL TUBE3D(xlast,ylast,zlast,xthis,ythis,zthis,0.03D0,ntuberes,ntuberes)
		end do
	end do
end if

!From now on, all objects plotted below use ambient of (1,1,1), which makes objects brighter than default (0.2,0.2,0.2)
call MATOP3(1D0,1D0,1D0,'ambient')

!Draw interbasin surfaces
if (idrawbassurf==1.and.numbassurf>0) then
	if (isurfstyle==1) then !A bunch of paths to represent the interbasin surface
		call MATOP3(0.7D0, 0.7D0, 0.8D0, 'diffuse')
		do isurf=1,numbassurf
			do ipath=1,nsurfpathpercp
				do ipt=2,nsurfpt
					CALL TUBE3D(bassurpath(1,ipt-1,ipath,isurf),bassurpath(2,ipt-1,ipath,isurf),bassurpath(3,ipt-1,ipath,isurf)&
					,bassurpath(1,ipt,ipath,isurf),bassurpath(2,ipt,ipath,isurf),bassurpath(3,ipt,ipath,isurf),0.005D0,5,5)
				end do
			end do
		end do
	else if (isurfstyle==2) then !Use triangles to consist the interbasin surfaces
		call MATOP3(IBSclrR, IBSclrG, IBSclrB, 'diffuse')
		interval=2
		do isurf=1,numbassurf
			do ipath=1,nsurfpathpercp
				ipathp1=ipath+1
				if (ipath==nsurfpathpercp) ipathp1=1
				do ipt=1,nsurfpt-interval,interval
					!Use two triangle to comprise a tile of face. Define the two triangle in clockwise manner first
					trianglex(:)=(/ bassurpath(1,ipt,ipath,isurf),bassurpath(1,ipt+interval,ipath,isurf),bassurpath(1,ipt+interval,ipathp1,isurf) /)
					triangley(:)=(/ bassurpath(2,ipt,ipath,isurf),bassurpath(2,ipt+interval,ipath,isurf),bassurpath(2,ipt+interval,ipathp1,isurf) /)
					trianglez(:)=(/ bassurpath(3,ipt,ipath,isurf),bassurpath(3,ipt+interval,ipath,isurf),bassurpath(3,ipt+interval,ipathp1,isurf) /)
					call TRIA3D(trianglex,triangley,trianglez)
					trianglex(:)=(/ bassurpath(1,ipt,ipath,isurf),bassurpath(1,ipt+interval,ipathp1,isurf),bassurpath(1,ipt,ipathp1,isurf) /)
					triangley(:)=(/ bassurpath(2,ipt,ipath,isurf),bassurpath(2,ipt+interval,ipathp1,isurf),bassurpath(2,ipt,ipathp1,isurf) /)
					trianglez(:)=(/ bassurpath(3,ipt,ipath,isurf),bassurpath(3,ipt+interval,ipathp1,isurf),bassurpath(3,ipt,ipathp1,isurf) /)
					call TRIA3D(trianglex,triangley,trianglez)
					!Define the two triangle in anti-clockwise manner. We must draw the surfaces two times from different view,
					!otherwise the surface will be invisible from certain viewpoint
					trianglex(:)=(/ bassurpath(1,ipt,ipath,isurf),bassurpath(1,ipt+interval,ipathp1,isurf),bassurpath(1,ipt+interval,ipath,isurf) /)
					triangley(:)=(/ bassurpath(2,ipt,ipath,isurf),bassurpath(2,ipt+interval,ipathp1,isurf),bassurpath(2,ipt+interval,ipath,isurf) /)
					trianglez(:)=(/ bassurpath(3,ipt,ipath,isurf),bassurpath(3,ipt+interval,ipathp1,isurf),bassurpath(3,ipt+interval,ipath,isurf) /)
					call TRIA3D(trianglex,triangley,trianglez)
					trianglex(:)=(/ bassurpath(1,ipt,ipath,isurf),bassurpath(1,ipt,ipathp1,isurf),bassurpath(1,ipt+interval,ipathp1,isurf) /)
					triangley(:)=(/ bassurpath(2,ipt,ipath,isurf),bassurpath(2,ipt,ipathp1,isurf),bassurpath(2,ipt+interval,ipathp1,isurf) /)
					trianglez(:)=(/ bassurpath(3,ipt,ipath,isurf),bassurpath(3,ipt,ipathp1,isurf),bassurpath(3,ipt+interval,ipathp1,isurf) /)
					call TRIA3D(trianglex,triangley,trianglez)
				end do
			end do
		end do
	else if (isurfstyle==3) then !Use a lot of cylinders between adjacent paths to portray interbasin surfaces
		call MATOP3(IBSclrR, IBSclrG, IBSclrB, 'diffuse')
		do isurf=1,numbassurf
			do ipath=1,nsurfpathpercp
				ipathp1=ipath+1
				if (ipath==nsurfpathpercp) ipathp1=1
				do ipt=1,nsurfpt
					CALL TUBE3D(bassurpath(1,ipt,ipath,isurf),bassurpath(2,ipt,ipath,isurf),bassurpath(3,ipt,ipath,isurf)&
					,bassurpath(1,ipt,ipathp1,isurf),bassurpath(2,ipt,ipathp1,isurf),bassurpath(3,ipt,ipathp1,isurf),0.03D0,3,3)
				end do
			end do
		end do
	end if
end if

!Draw spheres corresponding to local minimum and maximum derived from molecular surface analysis
if (GUI_mode==5) then
	if (ishowlocminpos==1) then
		CALL MATOP3(0D0, 0D0, 1D0, 'diffuse') !Blue
		do i=1,nsurlocmin
			idxtmp=surlocminidx(i)
			if (idxtmp==0) cycle !The extreme has already been discarded
			CALL SPHE3D(survtx(idxtmp)%x,survtx(idxtmp)%y,survtx(idxtmp)%z,0.15D0,20,20)
		end do
	end if
	if (ishowlocmaxpos==1) then
		CALL MATOP3(1D0, 0D0, 0D0, 'diffuse') !Red
		do i=1,nsurlocmax
			idxtmp=surlocmaxidx(i)
			if (idxtmp==0) cycle !The extreme has already been discarded
			CALL SPHE3D(survtx(idxtmp)%x,survtx(idxtmp)%y,survtx(idxtmp)%z,0.15D0,20,20)
		end do
	end if
end if

!When isosur1style==5(transparent), isosur2style must be also 5. For transparent style, TPRINI/TPRFIN conflict with Z-buffer, so here we need to recall ZBFFIN earlier
!When one of isosur1style and isosur2style is unequal to 5, then another must not be 5. Overall, we ensure that the circumstance that only one isosurface is transparent will not occured
if (isosur1style==5) CALL ZBFFIN

if (idrawisosur==1) then
	!Set lighting parameter for showing both isosurface 1 and 2
 	call litpos(1,XVU,YVU,ZVU,'ANGLE')
	call litpos(2,XVU+90,YVU+90,ZVU,'ANGLE')
	call litpos(3,XVU+90,YVU-90,ZVU,'ANGLE')
 	call litpos(4,XVU-90,YVU+90,ZVU,'ANGLE')
 	call litpos(5,XVU-90,YVU-90,ZVU,'ANGLE')
!  	call litpos(6,XVU+50,YVU,ZVU,'ANGLE')
!  	call litpos(7,XVU,YVU+50,ZVU,'ANGLE')
!  	call litpos(8,XVU-50,YVU,ZVU,'ANGLE')
!  	call litpos(8,-2*abslenx,2*absleny,2.5*abslenz,'ANGLE')
	if (ienablelight1==1) call litmod(1,'on')
	if (ienablelight2==1) call litmod(2,'on')
	if (ienablelight3==1) call litmod(3,'on')
	if (ienablelight4==1) call litmod(4,'on')
	if (ienablelight5==1) call litmod(5,'on')
	if (ienablelight1==0) call litmod(1,'off')
	if (ienablelight2==0) call litmod(2,'off')
	if (ienablelight3==0) call litmod(3,'off')
	if (ienablelight4==0) call litmod(4,'off')
	if (ienablelight5==0) call litmod(5,'off')

	!Draw isosurface for cubmat
	do ix=1,nx
		arrayx(ix)=orgx+(ix-1)*gridv1(1)
	end do
	do iy=1,ny
		arrayy(iy)=orgy+(iy-1)*gridv2(2)
	end do
	do iz=1,nz
		arrayz(iz)=orgz+(iz-1)*gridv3(3)
	end do
	if (isosur1style==5) CALL TPRVAL(opacitycub1)
	nplottime=1
	if (isosurshowboth==1) nplottime=2 !Show both positive and negative regions
	do iplottime=1,nplottime
		CALL MATOP3(clrRcub1same,clrGcub1same,clrBcub1same,'diffuse') !Set color for solid isosurface 1 with the same sign of set isovalue
		call setrgb(clrRcub1samemeshpt,clrGcub1samemeshpt,clrBcub1samemeshpt)
        if (GUI_mode==1) then
			sur_valuenow=sur_value_orb
        else
			sur_valuenow=sur_value
        end if
		if (iplottime==2) then
			CALL MATOP3(clrRcub1oppo,clrGcub1oppo,clrBcub1oppo,'diffuse') !Set color for solid isosurface 1 with the opposite sign of set isovalue
			call setrgb(clrRcub1oppomeshpt,clrGcub1oppomeshpt,clrBcub1oppomeshpt)
			if (GUI_mode==1) then
				sur_valuenow=-sur_value_orb
			else
				sur_valuenow=-sur_value
			end if
		end if
		if (isosur1style==2) then !Plotted as mesh rather than solid face
            call surmsh("LINES")
		else if (isosur1style==-2) then !Plotted as mesh rather than solid face
            if (iplottime==2) call surmsh("LINES") !Use mesh tyle only for negative part
        end if
		if (isosur1style==3) call surmsh("POINTS") !Plotted as points rather than solid face
		if (isosur1style==4) then !Face+lines
            call surmsh("ON")
		else if (isosur1style==-4) then !Face+lines only for negative part
            if (iplottime==2) call surmsh("ON")
        end if
		if (isosur1style==5) CALL TPRINI !Must be called individually for same and opposite survalue, else same part will completely overlay opposite part
		call suriso(arrayx,nx,arrayy,ny,arrayz,nz,cubmat,sur_valuenow)
		if (isosur1style==5) CALL TPRFIN
		call surmsh("OFF") !If don't set this, then other things (atoms, bonds) will be drawn as LINES too
	end do
	
	!Draw isosurface for cubmattmp at the same time
	if (isosur2style==5) CALL TPRVAL(opacitycub2)
	if (isosursec==1.and.allocated(cubmattmp)) then
		nplottime=1
		if (isosurshowboth==1) nplottime=2 !Show both positive and negative region
		do iplottime=1,nplottime
			CALL MATOP3(clrRcub2same,clrGcub2same,clrBcub2same,'diffuse') !Set color for solid isosurface 2 with the same sign of set isovalue
			call setrgb(clrRcub2samemeshpt,clrGcub2samemeshpt,clrBcub2samemeshpt)
			if (GUI_mode==1) then
				sur_valuenow=sur_value_orb
			else
				sur_valuenow=sur_value
			end if
			if (iplottime==2) then
				CALL MATOP3(clrRcub2oppo,clrGcub2oppo,clrBcub2oppo,'diffuse') !Set color for solid isosurface 2 with the opposite sign of set isovalue
				call setrgb(clrRcub2oppomeshpt,clrGcub2oppomeshpt,clrBcub2oppomeshpt)
				if (GUI_mode==1) then
					sur_valuenow=-sur_value_orb
				else
					sur_valuenow=-sur_value
				end if
			end if
			if (isosur2style==2) call surmsh("LINES") !Plotted as mesh rather than solid face
			if (isosur2style==3) call surmsh("POINTS") !Plotted as points rather than solid face
			if (isosur2style==4) call surmsh("ON") !face+lines
			if (isosur1style==5) CALL TPRINI
			call suriso(arrayx,nx,arrayy,ny,arrayz,nz,cubmattmp,sur_valuenow)
			if (isosur1style==5) CALL TPRFIN
			call surmsh("OFF")
		end do
	end if
end if

!Draw a 3D rectangle box to show spatial range of present grid data or grid data to be calculated
tubethk=0.07D0
if (ishowdatarange==1) then
	CALL MATOP3(0D0,0D0,0.8D0,'diffuse')
    call gridvertex2(1,2,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(2,3,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(3,4,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(1,4,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(1,5,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(2,6,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(3,7,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(4,8,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(5,6,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(6,7,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(7,8,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call gridvertex2(5,8,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
end if

!Show cell frame
if (ishowcell==1) then
	!Along a
	CALL MATOP3(0.7D0,0.0D0,0.0D0,'diffuse')
    call cellvertex2(1,2,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(3,4,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(5,6,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(7,8,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    !Along b
	CALL MATOP3(0.0D0,0.7D0,0.0D0,'diffuse')
    call cellvertex2(2,3,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(1,4,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(6,7,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(5,8,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    !Along c
	CALL MATOP3(0.0D0,0.0D0,0.7D0,'diffuse')
    call cellvertex2(1,5,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(2,6,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(3,7,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
    call cellvertex2(4,8,vert1x,vert1y,vert1z,vert2x,vert2y,vert2z)
	call tube3D(vert1x,vert1y,vert1z,vert2x,vert2y,vert2z,tubethk,30,30)
end if

if (isosur1style/=5) CALL ZBFFIN !Ending of Z-buffer

!Draw label of atom name and index of atoms/CPs/paths/surface extremes/real attractors in 3D plot
if ((ishowatmlab==1.or.ishowCPlab==1.or.ishowpathlab==1.or.ishowlocminlab==1.or.ishowlocmaxlab==1.or.(ishowattlab==1.and.numatt>0)).and.textheigh>0) then
    textheighmod=textheigh+155*plot2abs**1.3D0 !Change text size according to molecule size
	call height(int(textheighmod))
	if (ishowatmlab==1) then
		call setRGB(atmlabclrR,atmlabclrG,atmlabclrB)
		do i=1,ncenter_tmp
            if (ishowhydrogen==0.and.a_tmp(i)%index==1) cycle
			write(ctemp,"(i5)") a_tmp_idx(i)
			absx=(a_tmp(i)%x-(xhigh+xlow)/2) * plot2abs !Find atomic absolute coordinate in the axis
			absy=(a_tmp(i)%y-(yhigh+ylow)/2) * plot2abs
 			absz=(a_tmp(i)%z-(zhigh+zlow)/2) * plot2abs
			call abs3pt(absx,absy,absz,xplotcoor,yplotcoor) !Convert atomic absolute coordinate in the axis to screen coordinate (pixel)
			screeny=nint(yplotcoor-textheighmod/1.8D0)
			if (iatmlabtype3D==1) then !Show element
				screenx=nint(xplotcoor-textheighmod/2D0)
				if (a_tmp(i)%name(2:2)/=" ") screenx=nint(xplotcoor-textheighmod/1.3D0)
				call messag(trim(a_tmp(i)%name),screenx,screeny)
			else if (iatmlabtype3D==2) then !Show index all atoms (including Bq)
				screenx=nint(xplotcoor-textheighmod/2D0)
				if (i>=10) screenx=nint(xplotcoor-textheighmod/1.7D0)
				call messag(ADJUSTL(ctemp),screenx,screeny)
			else if (iatmlabtype3D==4.or.iatmlabtype3D==5) then !Show index	only for Bq atoms
				if (a_tmp(i)%index/=0) cycle
				screenx=nint(xplotcoor-textheighmod/2D0)
				if (iatmlabtype3D==4) then
					if (i>=10) screenx=nint(xplotcoor-textheighmod/1.7D0)
				else
					itmp=count(a_tmp(1:i)%index==0)
					write(ctemp,"(i5)") itmp
					if (itmp>=10) screenx=nint(xplotcoor-textheighmod/1.7D0)
				end if
				call messag(ADJUSTL(ctemp),screenx,screeny)
			else !Show element and index
                if (iatmlabtype3D==6.and.a_tmp(i)%index==0) cycle !Do not show Bq label
				if (i<10) then !Slightly move center of text so that they can at center of atom
					screenx=nint(xplotcoor-textheighmod/1.1D0)
				else !Move in X more, since the index has two or more digitals
					screenx=nint(xplotcoor-textheighmod/0.8D0)
				end if
				call messag(trim(a_tmp(i)%name)//ADJUSTL(ctemp),screenx,screeny)
			end if
		end do
	end if
	if (ishowCPlab==1) then !Draw label for critical points
		CALL SETRGB(CPlabclrR,CPlabclrG,CPlabclrB)
		do i=1,numcp_tmp
			if (CPtype_tmp(i)==1.and.ishow3n3==0) cycle
			if (CPtype_tmp(i)==2.and.ishow3n1==0) cycle
			if (CPtype_tmp(i)==3.and.ishow3p1==0) cycle
			if (CPtype_tmp(i)==4.and.ishow3p3==0) cycle
            if (lab_oneCP>0.and.lab_oneCP/=i) cycle
			write(ctemp,"(i5)") CP_tmp_idx(i)
			absx=(CPpos_tmp(1,i)-(xhigh+xlow)/2) * plot2abs !Find atomic absolute coordinate
			absy=(CPpos_tmp(2,i)-(yhigh+ylow)/2) * plot2abs
 			absz=(CPpos_tmp(3,i)-(zhigh+zlow)/2) * plot2abs
			call abs3pt(absx,absy,absz,xplotcoor,yplotcoor) !Convert atomic absolute coordinate to screen coordinate(pixel)
			screeny=nint(yplotcoor-textheighmod/1.8)
            if (i<10) then
				screenx=nint(xplotcoor-textheighmod/2.6)
            else if (i<100) then
				screenx=nint(xplotcoor-textheighmod/1.5)
            else
				screenx=nint(xplotcoor-textheighmod)
			end if
			call messag(trim(ADJUSTL(ctemp)),screenx,screeny)
		end do
	end if
	if (ishowpathlab==1) then
		call color("RED")
		do ipath=1,numpath_tmp
			write(ctemp,"(i5)") path_tmp_idx(ipath)
			ipathmidpt=nint(pathnumpt_tmp(ipath)/2D0)
			absx=(topopath_tmp(1,ipathmidpt,ipath)-(xhigh+xlow)/2) * plot2abs !Find atomic absolute coordinate
			absy=(topopath_tmp(2,ipathmidpt,ipath)-(yhigh+ylow)/2) * plot2abs
 			absz=(topopath_tmp(3,ipathmidpt,ipath)-(zhigh+zlow)/2) * plot2abs
			call abs3pt(absx,absy,absz,xplotcoor,yplotcoor) !Convert atomic absolute coordinate to screen coordinate(pixel)
			screenx=nint(xplotcoor-textheighmod/2.6D0)
			screeny=nint(yplotcoor-textheighmod/1.8D0)
			call messag(trim(ADJUSTL(ctemp)),screenx,screeny)
		end do
	end if
	if (ishowlocminlab==1) then
		call color("MAGENTA")
		do i=1,nsurlocmin
			idxtmp=surlocminidx(i)
			if (idxtmp==0) cycle !The extreme has already been discarded
			write(ctemp,"(i5)") i
			absx=(survtx(idxtmp)%x-(xhigh+xlow)/2) * plot2abs !Find absolute coordinate
			absy=(survtx(idxtmp)%y-(yhigh+ylow)/2) * plot2abs
 			absz=(survtx(idxtmp)%z-(zhigh+zlow)/2) * plot2abs
 			call abs3pt(absx,absy,absz,xplotcoor,yplotcoor) !Convert atomic absolute coordinate to screen coordinate(pixel)
			screenx=nint(xplotcoor-textheighmod/2.6D0)
			screeny=nint(yplotcoor-textheighmod/1.8D0)
			call messag(trim(ADJUSTL(ctemp)),screenx,screeny)
		end do
	end if
	if (ishowlocmaxlab==1) then
		call color("GREEN")
		do i=1,nsurlocmax
			idxtmp=surlocmaxidx(i)
			if (idxtmp==0) cycle !The extreme has already been discarded
			write(ctemp,"(i5)") i
			absx=(survtx(idxtmp)%x-(xhigh+xlow)/2) * plot2abs !Find absolute coordinate
			absy=(survtx(idxtmp)%y-(yhigh+ylow)/2) * plot2abs
 			absz=(survtx(idxtmp)%z-(zhigh+zlow)/2) * plot2abs
 			call abs3pt(absx,absy,absz,xplotcoor,yplotcoor) !Convert atomic absolute coordinate to screen coordinate(pixel)
			screenx=nint(xplotcoor-textheighmod/2.6D0)
			screeny=nint(yplotcoor-textheighmod/1.8D0)
			call messag(trim(ADJUSTL(ctemp)),screenx,screeny)
		end do
	end if
	if (ishowattlab==1.and.numatt>0) then
		CALL SETRGB(0.70D0,0.3D0,0.9D0)
		do iatt=1,numatt
			irealatt=attconv(iatt)
			write(ctemp,"(i5)") irealatt
			absx=(attxyz(1,iatt)-(xhigh+xlow)/2) * plot2abs !Find atomic absolute coordinate
			absy=(attxyz(2,iatt)-(yhigh+ylow)/2) * plot2abs
 			absz=(attxyz(3,iatt)-(zhigh+zlow)/2) * plot2abs
			call abs3pt(absx,absy,absz,xplotcoor,yplotcoor) !Convert atomic absolute coordinate to screen coordinate(pixel)
			screenx=nint(xplotcoor-textheighmod/2.6)
			screeny=nint(yplotcoor-textheighmod/1.8)
			call messag(trim(ADJUSTL(ctemp)),screenx,screeny)
		end do
	end if
end if
CALL DISFIN
XVU=XVUold
YVU=YVUold
if (allocated(a_tmp)) deallocate(a_tmp,a_tmp_idx)
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!------------------------- Draw property in a line
!if atomr1==atomr2, then don't plot two points to highlight nuclear positions
!Input coordinate must be in Bohr
subroutine drawcurve(curvex,curvey,N,curvexmin,curvexmax,steplabx,curveymin,curveymax,steplaby,status,atomr1,atomr2,axisnamey)
implicit real*8 (a-h,o-z)
real*8 curvex(N),curvey(N),atompointx(2),atompointy(2),curvexmin,curvexmax,curveymin,curveymax,steplabx,steplaby
real*8,optional :: atomr1,atomr2
real*8 vertlinex(2),vertliney(2)
integer N
character status*4
character(len=*),optional :: axisnamey
!Set conversion factor, determined by global variable ilenunit1D
scll=1D0 !Default, namely Bohr as unit of X-axis
if (ilenunit1D==2) scll=b2a !Conver X-axis unit to Angstrom

call SCRMOD('REVERSE')
CALL PAGE(2970,2100)
CALL IMGFMT("RGB")
CALL setxid(0,'NONE') !If we don't set this, after we draw a graph embedded in GUI(e.g. relif map), curve map will not be shown 
if (status=="show") then
	call METAFL('xwin')
	call window(200,100,900,600)
else if (status=="save") then
	call METAFL(graphformat)
	call winsiz(graph1Dwidth,graph1Dheight)
end if
CALL DISINI
! call LINMOD ('ON','SMOOTH') !It seems that Anti-aliased doesn't take effect
if (status=="show".and.isys==1) then
	call WINTIT("Curve graph, click right mouse button to continue...")
	call height(45)
	CALL HNAME(45)
else
	call height(40) !The text shown in graphic file is strangely larger than window, so slight decrease it
end if
CALL NAMDIS(40,'X')
CALL NAMDIS(50,'Y')
if (status=="show".or.graphformat=="pdf ") then
	CALL HWFONT
else if (status=="save") then
	if (ttfontfile=="none") then
		CALL HELVES
    else
		CALL TTFONT(ttfontfile)
    end if
    CALL SHDCHA
end if
call center
nysize=nint(2300*curvexyratio)
call AXSLEN(2300,nysize)
if (nysize>1800) call AXSLEN(nint(2300/(nysize/1800D0)),1800)
if (ilog10y==0) then
    call AXSSCL('lin','Y')
    CALL LABDIG(numdigliney,"Y")
else if (ilog10y==1) then
    call AXSSCL('log','Y')
    call labels('log','Y')
    CALL LABDIG(-1,"Y") !Do not show digit in Y
end if
if (ilenunit1D==1) CALL NAME('Position (Bohr)','X')
if (ilenunit1D==2) CALL NAME('Position (Angstrom)','X')
if (present(axisnamey)) then
	if (axisnamey==" ") then
		CALL NAME('Value (a.u.)','Y')
    else
		CALL NAME(trim(axisnamey),'Y')
    end if
else
	CALL NAME('Value (a.u.)','Y')
end if
CALL LABDIG(numdiglinex,"X")
CALL TICPOS("REVERS","XYZ")
call height(curve_axistextsize)
CALL HNAME(curve_axisnamesize)
call ERRMOD("ALL","OFF")
shifty=mod(curveymin,steplaby)
if (ilog10y==0) then
	CALL GRAF(curvexmin*scll,curvexmax*scll,curvexmin*scll,steplabx, curveymin,curveymax,curveymin-shifty,steplaby)
else if (ilog10y==1) then
	CALL GRAF(curvexmin*scll,curvexmax*scll,curvexmin*scll,steplabx, curveymin,curveymax,curveymin-shifty,1D0)
end if
CALL SETRGB(0.6D0,0.6D0,0.6D0)
CALL DASH
CALL XAXGIT !Draw the dashed line of Y=0

if (icurve_vertlinex==1) then
	vertlinex=curve_vertlinex !Draw a vertical line throughing out the whole graph to help locating special position
	if (ilog10y==0) then
		vertliney(1)=curveymin
		vertliney(2)=curveymax
	else if (ilog10y==1) then
		vertliney(1)=10**(curveymin)
		vertliney(2)=10**(curveymax)
	end if
	CALL LINWID(icurvethick-2)
	CALL CURVE(vertlinex,vertliney,2)
end if

call setcolor(iclrcurve)
CALL SOLID
CALL LINWID(icurvethick)
CALL CURVE(curvex*scll,curvey,N)
CALL LINWID(1)

if (present(atomr1)) then !Draw position of the two atom selected
	if (atomr1/=atomr2) then
		atompointx(1)=atomr1*scll !the position of atom shown in curve graph
		atompointx(2)=atomr2*scll
		atompointy=curveymin
		CALL SETRGB(0.9D0,0D0,0D0)
		CALL INCMRK(-1)
		CALL MARKER(21)
		CALL HSYMBL(25)
		CALL CURVE(atompointx,atompointy,2)
    end if
end if
CALL DISFIN
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!------------ Draw scatter graph for two functions
!This routine can also simultaneously plot two set of scatter data, as long as arrayx2/y2 are defined, they will be plotted as red points
!iratio=1: graph is 4:3  =2: graph is 1:1
subroutine drawscatter(arrayx,arrayy,ntot,xmin,xmax,ymin,ymax,iratio,stringx,stringy,arrayx2,arrayy2,ntot2)
real*8 arrayx(:),arrayy(:),xmin,xmax,ymin,ymax
integer ntot,iratio
character(len=*),optional :: stringx,stringy
integer,optional :: ntot2
real*8,optional :: arrayx2(:),arrayy2(:)
call SCRMOD('REVERSE')
if (iratio==1) then
	CALL PAGE(3000,2250) !4:3
	if (isavepic==0) then
		call METAFL('xwin')
		call window(100,100,800,600)
	else if (isavepic==1) then
		call METAFL(graphformat)
		call winsiz(graph1Dwidth,nint(graph1Dwidth/4D0*3D0)) !Ensure 4:3
	end if
else if (iratio==2) then
	CALL PAGE(3000,3000) !1:1
	if (isavepic==0) then
		call METAFL('xwin')
		call window(100,100,700,700)
	else if (isavepic==1) then
		call METAFL(graphformat)
		call winsiz(graph1Dwidth,graph1Dwidth) !Ensure 1:1
	end if
end if
CALL setxid(0,'NONE')
CALL DISINI
if (isavepic==0) call WINTIT("Scatter graph between two functions, click right mouse button to continue")
call ERRMOD("ALL","OFF")
if (isavepic==0.or.graphformat=="pdf ") then
	call hwfont
else if (isavepic==1) then
	if (ttfontfile=="none") then
		CALL HELVES
	else
		CALL TTFONT(ttfontfile)
	end if
    CALL SHDCHA
end if
call center
if (iratio==1) then
	call AXSLEN(2400,1800) !4:3
else if (iratio==2) then
	call AXSLEN(2400,2400) !1:1
end if
CALL HNAME(45)
CALL height(40)
call TEXMOD("ON")
if (present(stringx)) then
	CALL NAME(stringx,'X')
else
	CALL NAME('Function 1 value','X')
end if
if (present(stringy)) then
	CALL NAME(stringy,'Y')
else
	CALL NAME('Function 2 value','Y')
end if
CALL LABDIG(numdiglinex,"X")
CALL LABDIG(numdigliney,"Y")
CALL TICPOS("REVERS","XYZ")
CALL GRAF(xmin,xmax,xmin,(xmax-xmin)/10,ymin,ymax,ymin,(ymax-ymin)/10)
CALL INCMRK(-1)
CALL MARKER(21)
CALL HSYMBL(symbolsize)
call curve(arrayx,arrayy,ntot)
if (present(arrayx2)) then
	call color("RED")
	call curve(arrayx2,arrayy2,ntot2)
	call color("WHITE") !Restore to default (black)
end if
CALL DISFIN
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!----------Draw a matrix as color-mapped graph
!mat is the 2D matrix to be plotted, the dimension is numx,numy
!ninterpo determines the number of interpolation
!x/ymin,x/ymax are the lower and upper limit of X and Y axes
!zmin and zmax are the lower and upper limit of color bar (Z axis)
!stepx/y/z are stepsizes between labels of X/Y/Z axes
!If nlabdig=-1, the labels of X and Y axes will be integer, this is suitable for natively representation of a matrix. If =n, there will be n digitals
!textheight is size of label, optional, default is 45
!nlabdigz is the number of digial of labels in Z axis, optional, default is 4
subroutine drawmatcolor(mat,numx,numy,xmin,xmax,ymin,ymax,zmin,zmax,stepx,stepy,ninterpo,nlabdig,textheight,nlabdigz,stepzin)
implicit real*8 (a-h,o-z)
integer numx,numy,ninterpo,nlabdig
real*8 mat(numx,numy),xmin,xmax,ymin,ymax,zmin,zmax,stepx,stepy
integer :: lengthx=2300
integer,optional :: textheight,nlabdigz
real*8,optional :: stepzin
!Truncate the values larger than and lower than color scale, so that these regions will not be shown as white and black, respectively
if (iclrtrans/=0) then
	where (mat>zmax) mat=zmax-1D-10   !Augment by a minimal value to avoid numerical noise
	where (mat<zmin) mat=zmin+1D-10
end if
call SCRMOD('REVERSE')
CALL setxid(0,'NONE')
CALL PAGE(3200,2700)
if (isavepic==0) then
	call METAFL('xwin')
	call window(200,100,900,770)
else if (isavepic==1) then
	call METAFL(graphformat)
	call winsiz(graph2Dwidth,graph2Dheight) !Because color bar need to draw, so width is bigger than height
	CALL IMGFMT('RGB')
end if
call DISINI
call ERRMOD("ALL","OFF") !If don't set this, when atom label in contour map is out of page range, DISLIN annoys users
if (isavepic==0.or.graphformat=="pdf ") then
	call HWFONT
else if (isavepic==1) then
	if (ttfontfile=="none") then
		CALL HELVES
    else
		CALL TTFONT(ttfontfile)
    end if
	CALL SHDCHA
end if
CALL LABDIG(nlabdig,"X") !-1 means integer label
CALL LABDIG(nlabdig,"Y")
if (present(nlabdigz)) then
	CALL LABDIG(nlabdigz,"Z")
else
	CALL LABDIG(4,"Z")
end if
if (present(textheight)) then
	call height(textheight)
else
	call height(45)
end if
if (present(stepzin)) then
	stepz=stepzin
else
	stepz=(zmax-zmin)/10D0
end if
call ticks(1,"XYZ")
call WINTIT("Colored matrix")
call center
call AUTRES(numx,numy)
call AX3LEN(lengthx,nint(lengthx*dfloat(numy)/numx),nint(lengthx*dfloat(numy)/numx))
call setcolortable(iclrtrans)
if (ninterpo==1) then !Don't interpolate
	call sursze(1D0,dfloat(numx),1D0,dfloat(numx)) !Manually set center position of starting and ending grids to ensure boundary grids have the same size as internal grids
	call GRAF3(xmin-0.5D0,xmax+0.5D0,xmin,stepx,ymin-0.5D0,ymax+0.5D0,ymin,stepy,zmin,zmax,zmin,stepz)
else
	call GRAF3(xmin,xmax,xmin,stepx,ymin,ymax,ymin,stepy,zmin,zmax,zmin,stepz)
end if
call CRVMAT(mat,numx,numy,ninterpo,ninterpo)
call DISFIN
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!------------------------- Draw property on a plane (or output critical points and paths, when iplaneoutall=1)
!init and end is coordinate range in each dimension
!Input length unit must be in Bohr, but you can use ilenunit2D to change the displayed length unit (=1/2 denote Bohr/Angstrom)
subroutine drawplane(init1inp,end1inp,init2inp,end2inp,init3,end3)
use topo
use util
use functions
implicit real*8 (a-h,o-z)
real*8 init1inp,end1inp,init2inp,end2inp,init1,end1,init2,end2,init3,end3
real*8 xcoord(ngridnum1),ycoord(ngridnum2),gradd1tmp(ngridnum1,ngridnum2),gradd2tmp(ngridnum1,ngridnum2)
real*8 dx2D,dy2D,pix2usr,n1,n2,n1_2,n2_2
real*8 planetrunc(ngridnum1,ngridnum2),planetrunc2(ngridnum1,ngridnum2) !Store truncated planemat
integer lengthx !length of x axis
integer i,j
integer :: inplane(ncenter+1) !If =1, then the label is close enough to the plotting plane
character atmlabtext*20,tmpstr*20
!For plotting extrema on a contour
integer,parameter :: maxctrpt=10000,maxctrline=300
real*8 ctrptx(maxctrpt),ctrpty(maxctrpt) !X,Y coordinate of points on all contour lines of the specific isovalue
real*8 ctrptval(maxctrpt),ctrptval_tmp(maxctrpt),ctrptval_org(maxctrpt) !Another function on the contour line points
integer ctrptnum(maxctrline) !Number of points on each contour line
integer ifopenctr(maxctrline)

scll=1D0
if (ilenunit2D==2) scll=b2a
init1=init1inp*scll !We use init1inp as an intermediate variable rather than directly use init1, because orgx2D will be passed as init1, we don't want to alternate orgx2D
init2=init2inp*scll
end1=end1inp*scll
end2=end2inp*scll
disshowlabel=disshowlabel*scll
disshowCP=disshowCP*scll
disshowpath=disshowpath*scll
if (iplaneoutall==1) goto 10 !Don't plot anything, but only output coordinate of critical points and paths, and then directly return
if (ilenunit2D==2) call convgridlenunit(1) !Convert plane parameters to Angstrom. At the end of the routine, they will be converted back

lengthx=2300
call SCRMOD('REVERSE')
if (isavepic==0) then
	if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
		CALL setxid(idisgraph, 'WIDGET')
		call METAFL('CONS')
		CALL PAGE(3000,3000)
	else
		CALL setxid(0,'NONE')
		call METAFL('xwin')
		call window(200,100,913,770) !Same ratio as 3200:2700
		CALL PAGE(3200,2700)
	end if
else if (isavepic==1) then
	CALL setxid(0,'NONE')
	call METAFL(graphformat)
	if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
		call winsiz(graph3Dwidth,graph3Dheight)
	else
		call winsiz(graph2Dwidth,graph2Dheight) !Because color bar need to draw, so width is bigger than height
	end if
	if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
		CALL PAGE(3000,3000)
	else
		CALL PAGE(3200,2700)
	end if
	CALL IMGFMT('RGB')
end if
call DISINI
call ERRMOD("ALL","OFF") !If don't set this, when atom label in contour map is out of page range, DISLIN annoys users
CALL VIEW3D(XVU,YVU,ZVU,"ANGLE")
CALL erase
if (isavepic==0.or.graphformat=="pdf ") then
	CALL HWFONT
else if (isavepic==1) then
	if (ttfontfile=="none") then
		CALL HELVES
    else
		CALL TTFONT(ttfontfile)
    end if
	CALL SHDCHA
end if
if (itickreverse==1) CALL TICPOS("REVERSE","XYZ")
if (numdigx>0) then
    CALL LABDIG(numdigx,"X")
else
    CALL LABDIG(-1,"X") !Integer label, no decimal place
end if
if (numdigy>0) then
    CALL LABDIG(numdigy,"Y")
else
    CALL LABDIG(-1,"Y") !Integer label, no decimal place
end if
if (numdigz>0) then
    CALL LABDIG(numdigz,"Z")
else
    CALL LABDIG(-1,"Z") !Integer label, no decimal place
end if
if (isavepic==0.and.isys==1) then
	call height(plane_axistextsize)
	CALL HNAME(plane_axisnamesize)
else !The text shown in graphic file is strangely larger than window, so slight decrease it
	call height(nint(plane_axistextsize*0.81))
	CALL HNAME(nint(plane_axisnamesize*0.9))
end if
call ticks(iticks,"XYZ")
if (ilenunit2D==1) CALL NAME('Length unit: Bohr', 'X')
if (ilenunit2D==2) CALL NAME('Length unit: Angstrom', 'X')
dx2D=(end1-init1)/(ngridnum1-1)
dy2D=(end2-init2)/(ngridnum2-1)
do i=1,ngridnum1
	xcoord(i)=init1+(i-1)*dx2D
end do
do i=1,ngridnum2
	ycoord(i)=init2+(i-1)*dy2D
end do
shiftx=mod(init1,planestpx) !Shift of axis origin, so that there is a label just at exactly midpoint
shifty=mod(init2,planestpy)
shiftz=0
! shiftz=mod(init3,planestpz)

!1:Color-filled map, 2:Contour line map, 6:Gradient lines map, 7:Vector field map
if (idrawtype==1.or.idrawtype==2.or.idrawtype==6.or.idrawtype==7) then
	if (isavepic==0) then
		if (idrawtype==1) call WINTIT("Color-filled map, click right mouse button to continue")
		if (idrawtype==2) call WINTIT("Contour line map, click right mouse button to continue")
		if (idrawtype==6) call WINTIT("Gradient line map, click right mouse button to continue")
		if (idrawtype==7) call WINTIT("Contour lines with gradient lines, click right mouse button to continue")
		if (idrawtype==8) call WINTIT("Contour lines + vector field, click right mouse button to continue")
	end if
	!Length of Y is determined according to X by proportion, if length of Y axis is too large, then scale length x (and thus length y, automatically)
	if ( lengthx*(end2-init2)/(end1-init1)>2400 ) lengthx=2400/(end2-init2)*(end1-init1)
	lengthy=int(lengthx*(end2-init2)/(end1-init1))
    
    if ((idrawtype==2.or.idrawtype==6.or.idrawtype==7).and.ifillctrline==1.and.ishowclrfill_bar==1) then
        !Filling color for contour lines and meantime showing color bar, we need to proper set position of axis &
        !to leave enough region for the bar. "call center" should not used otherwise AXSPOS will not work
	    CALL AXSPOS((3200-(lengthx+300))/2,lengthy+(2700-(lengthy+150))/2)
    else
	    call center
    end if
    
    !Set color table and generate truncated plane data, they will be used in color-filled map or filling colors for contour lines
    call setcolortable(iclrtrans) !This routine must be invoked prior to GRAF
    if (iclrtrans/=0) then !Truncate the values larger than and lower than color scale, so that these regions will not be shown as white and black, respectively
		planetrunc=planemat
		where (planetrunc>end3) planetrunc=end3-1D-10   !Augment by a minimal value to avoid numerical noise
		where (planetrunc<init3) planetrunc=init3+1D-10
    else
        planetrunc=planemat
    end if
    
	if (idrawtype==1) then !Draw color color-filled map
		call AUTRES(ngridnum1,ngridnum2)
		call AX3LEN(lengthx,lengthy,lengthy) !The length of color bar is identical to Y axis
		call GRAF3(init1,end1,init1-shiftx,planestpx, init2,end2,init2-shifty,planestpy, init3,end3,init3-shiftz,planestpz)
        call CRVMAT(planetrunc,ngridnum1,ngridnum2,fillcoloritpx,fillcoloritpy)
	else !Set axis for other kinds of maps
		call axslen(lengthx,lengthy)
		call GRAF(init1,end1,init1-shiftx,planestpx, init2,end2,init2-shifty,planestpy)
	end if

! 	CALL LINWID(10) !Draw an arrow 
! 	call RLvec(0D0,0.216790D0,0D0,0.216790D0-0.229901D0*3,1421)
! 	call RLvec(1.424912D0,-0.867160D0,1.424912D0+0.17136000D0*3,-0.867160D0-0.14455800D0*3,1421)
! 	CALL LINWID(1)

	if (idrawcontour==1) then !Draw contour lines, may be also used for gradient lines and vector field map
        if (ifillctrline==1) then !Filling colors for contour lines
            call GETPOS(iPSX,iPSY) !Get position of actual user coordinate
            call ZSCALE(init3,end3) !Set lower and upper limit of color bar
            if (ishowclrfill_bar==1) call ZAXIS(init3,end3,init3,planestpz,lengthy," ",0,0,iPSX+lengthx+50,iPSY) !Show color bar
            CALL CONSHD(xcoord,ngridnum1,ycoord,ngridnum2,planetrunc,ctrval,ncontour) !Filling contour lines
        end if
		if (ilabel_on_contour==1) then !Enable showing isovalue on contour lines
			call height(ictrlabsize)
			call labels('float','contur')
			call labdig(numdigctr,'contur')
		end if
		do i=1,ncontour !Use different type of line to draw contour line
			if (ctrval(i)>=0) then
				CALL LINWID(iwidthposctr)
				call setcolor(iclrindctrpos)
				nsizestyle=2
				if (ctrposstyle(2)==0) nsizestyle=1
				call myline(ctrposstyle,nsizestyle) !Set line style
				call contur(xcoord,ngridnum1,ycoord,ngridnum2,planemat,ctrval(i))
				if (iorbsel2/=0) call contur(xcoord,ngridnum1,ycoord,ngridnum2,planemattmp,ctrval(i)) !Also plot another orbital if it was specified
			else
				CALL LINWID(iwidthnegctr)
				call setcolor(iclrindctrneg)
				nsizestyle=2
				if (ctrnegstyle(2)==0) nsizestyle=1
				call myline(ctrnegstyle,nsizestyle) !Set line style
				call contur(xcoord,ngridnum1,ycoord,ngridnum2,planemat,ctrval(i))
				if (iorbsel2/=0) call contur(xcoord,ngridnum1,ycoord,ngridnum2,planemattmp,ctrval(i))
			end if
		end do
		if (allocated(boldlinelist)) then !Bold some contour lines
			CALL LINWID(10)
			do i=1,size(boldlinelist)
				if (ctrval(boldlinelist(i))>=0D0) then
					call solid
				else
					call dashm
                end if
				call contur(xcoord,ngridnum1,ycoord,ngridnum2,planemat,ctrval(boldlinelist(i)))		
			end do
			CALL LINWID(1) !Restore to default
		end if
	end if

	if (idrawtype==6) then !Draw gradient line map
		gradd1tmp=gradd1 !Create temporary gradd1,d2 arrays
		gradd2tmp=gradd2
		if (iinvgradvec==1) then
			gradd1tmp=-gradd1tmp
			gradd2tmp=-gradd2tmp
        end if
		call solid
		call setcolor(iclrindgradline)
		if (igrad_arrow==1) then
			call stmmod('on','arrow')
		else
			call stmmod('off','arrow')
        end if
        call stmmod(stream_intmethod,'integration') !Integration method of stream line
		call stmval(gradplotstep,'step')
		call stmval(gradplotdis,'distance') !Controls number of lines, smaller value means more lines
		call stmval(gradplottest,'test')
		call stmval(0.3D0,'arrows') !Set interval distance of arrows in the same line
		call LINWID(iwidthgradline)
		call stream(gradd1tmp,gradd2tmp,ngridnum1,ngridnum2,xcoord,ycoord,(/ 0D0 /),(/ 0D0 /),0)
		CALL LINWID(1) !Restore to default
        
	else if (idrawtype==7) then !Draw vector field map
		gradd1tmp=gradd1 !Create temporary gradd1,d2 arrays
		gradd2tmp=gradd2
		do i=1,ngridnum1
			do j=1,ngridnum2
				rnorm=dsqrt(gradd1tmp(i,j)**2+gradd2tmp(i,j)**2)
				if (rnorm>cutgradvec) then
					gradd1tmp(i,j)=gradd1tmp(i,j)/rnorm*cutgradvec
					gradd2tmp(i,j)=gradd2tmp(i,j)/rnorm*cutgradvec
				end if
			end do
		end do
		if (iinvgradvec==1) then
			gradd1tmp=-gradd1tmp
			gradd2tmp=-gradd2tmp
        end if
		if (icolorvecfield==1) then
			call vecclr(-2)
		else
			call vecclr(vecclrind)
		end if
		call vecmat(gradd1tmp,gradd2tmp,ngridnum1,ngridnum2,xcoord,ycoord,1501)
	end if

    !Drawing of color-filled, contour, vector maps... has been finished, below we draw more widgets on the map

    !Draw vdW contour line (electron density=0.001)
    if (idrawplanevdwctr==1) then
	    call setcolor(ivdwclrindctr) !This routine must be invoked before LINWID, else no effect
	    if (ivdwctrlabsize/=0) then  !Enable showing isovalue on contour lines
		    call DISALF
		    call height(ivdwctrlabsize)
		    call labels('float','contur')
		    call labdig(3,'contur')
	    else
		    call labels('NONE','contur')		
	    end if
	    if (vdwctrstyle(2)==0) nsizestyle=1
	    call myline(vdwctrstyle,nsizestyle)
	    CALL LINWID(iwidthvdwctr)
	    call contur(xcoord,ngridnum1,ycoord,ngridnum2,planemattmp,0.001D0)
	    CALL LINWID(1) !Restore to default
	    call color("WHITE") !Restore to default
    end if

    !Draw all kinds of topology information on color-filled/contour/gradient/vector field graph
    10	continue
	!Draw topology paths on color-filled/contour/gradient/vector field graph. If iplaneoutall==1, then output the data points
	if (numpath>0.and.imarkpath==1) then
		if (iplaneoutall==0) then
			call setcolor(iclrpath)
			call HSYMBL(sizemarkpath)
		else if (iplaneoutall==1) then
			open(10,file="planepath.txt",status="replace")
		end if
		do ipath=1,numpath
			if ( plesel==1.and.any(abs(topopath(3,1:pathnumpt(ipath),ipath)*scll-orgz2D) > disshowpath) ) cycle
			if ( plesel==2.and.any(abs(topopath(2,1:pathnumpt(ipath),ipath)*scll-orgy2D) > disshowpath) ) cycle
			if ( plesel==3.and.any(abs(topopath(1,1:pathnumpt(ipath),ipath)*scll-orgx2D) > disshowpath) ) cycle
			if ( plesel==4.or.plesel==5.or.plesel==6.or.plesel==7 ) then
				ioutplane=0
				do ipt=1,pathnumpt(ipath)
					if (potpledis(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,&
					topopath(1,ipt,ipath)*scll,topopath(2,ipt,ipath)*scll,topopath(3,ipt,ipath)*scll)>disshowpath) then
						ioutplane=1
						exit
					end if
				end do
				if (ioutplane==1) cycle
			end if
			do ipt=1,pathnumpt(ipath)
				posmarkx=topopath(1,ipt,ipath)*scll
				posmarky=topopath(2,ipt,ipath)*scll
				posmarkz=topopath(3,ipt,ipath)*scll
				if (plesel==1) then
					if (posmarkx<init1.or.posmarkx>end1.or.posmarky<init2.or.posmarky>end2) cycle !To avoid path out of plotting range
					if (iplaneoutall==0) call rlsymb(21,posmarkx,posmarky)
					if (iplaneoutall==1) write(10,"(2f12.6)") posmarkx*b2a,posmarky*b2a
				else if (plesel==2) then
					if (posmarkx<init1.or.posmarkx>end1.or.posmarkz<init2.or.posmarkz>end2) cycle !To avoid path out of plotting range
					if (iplaneoutall==0) call rlsymb(21,posmarkx,posmarkz)
					if (iplaneoutall==1) write(10,"(2f12.6)") posmarkx*b2a,posmarkz*b2a
				else if (plesel==3) then
					if (posmarky<init1.or.posmarky>end1.or.posmarkz<init2.or.posmarkz>end2) cycle !To avoid path out of plotting range
					if (iplaneoutall==0) call rlsymb(21,posmarky,posmarkz)
					if (iplaneoutall==1) write(10,"(2f12.6)") posmarky*b2a,posmarkz*b2a
				else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
					call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,posmarkx,posmarky,posmarkz,prjx,prjy,prjz)
					if (abs(v1x*v2y-v2x*v1y)>1D-8) then
						det2_2=v1x*v2y-v2x*v1y
						n1=( (prjx-orgx2D)*v2y-v2x*(prjy-orgy2D) )/det2_2
						n2=( v1x*(prjy-orgy2D)-(prjx-orgx2D)*v1y )/det2_2
					else if (abs(v1x*v2z-v2x*v1z)>1D-8) then
						det2_2=v1x*v2z-v2x*v1z
						n1=( (prjx-orgx2D)*v2z-v2x*(prjz-orgz2D) )/det2_2
						n2=( v1x*(prjz-orgz2D)-(prjx-orgx2D)*v1z )/det2_2
					else if (abs(v1y*v2z-v2y*v1z)>1D-8) then
						det2_2=v1y*v2z-v2y*v1z
						n1=( (prjy-orgy2D)*v2z-v2y*(prjz-orgz2D) )/det2_2
						n2=( v1y*(prjz-orgz2D)-(prjy-orgy2D)*v1z )/det2_2
					end if
					tmp1=n1*d1; tmp2=n2*d2
					if (tmp1<init1.or.tmp1>end1.or.tmp2<init2.or.tmp2>end2) cycle !To avoid path out of plotting range
					if (iplaneoutall==0) call rlsymb(21,tmp1,tmp2)
					if (iplaneoutall==1) write(10,"(2f12.6)") n1*d1*b2a,n2*d2*b2a
				end if
			end do
			if (iplaneoutall==1) write(10,*) !Leave a blank line between each paths
		end do
		if (iplaneoutall==1) close(10)
	end if
	
	!Draw interbasin paths on color-filled/contour/gradient/vector field graph. If iplaneoutall==1, then output the data points
	if (nple3n1path>0.and.idrawintbasple==1) then
		if (iplaneoutall==0) then
			call setcolor(iclr3n1path)
			call HSYMBL(sizemark3n1path)
		else if (iplaneoutall==1) then
			open(10,file="planeinterbasin.txt",status="replace")
		end if
		do icp=1,nple3n1path
			do idir=1,2
				do ipt=1,n3n1plept
					posmarkx=ple3n1path(1,ipt,idir,icp)*scll
					posmarky=ple3n1path(2,ipt,idir,icp)*scll
					posmarkz=ple3n1path(3,ipt,idir,icp)*scll
					if (plesel==1) then
						if (abs(posmarkz-orgz2D) > disshowlabel) cycle
						if (posmarkx<init1.or.posmarkx>end1.or.posmarky<init2.or.posmarky>end2) cycle !To avoid path out of plotting range
						if (iplaneoutall==0) call rlsymb(21,posmarkx,posmarky)
						if (iplaneoutall==1) write(10,"(2f12.6)") posmarkx*b2a,posmarky*b2a
					else if (plesel==2) then
						if (abs(posmarky-orgy2D) > disshowlabel) cycle
						if (posmarkx<init1.or.posmarkx>end1.or.posmarkz<init2.or.posmarkz>end2) cycle !To avoid path out of plotting range
						if (iplaneoutall==0) call rlsymb(21,posmarkx,posmarkz)
						if (iplaneoutall==1) write(10,"(2f12.6)") posmarkx*b2a,posmarkz*b2a
					else if (plesel==3) then
						if (abs(posmarkx-orgx2D) > disshowlabel) cycle
						if (posmarky<init1.or.posmarky>end1.or.posmarkz<init2.or.posmarkz>end2) cycle !To avoid path out of plotting range
						if (iplaneoutall==0) call rlsymb(21,posmarky,posmarkz)
						if (iplaneoutall==1) write(10,"(2f12.6)") posmarky*b2a,posmarkz*b2a
					else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
						call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,posmarkx,posmarky,posmarkz,prjx,prjy,prjz)
						if ( (posmarkx-prjx)**2+(posmarky-prjy)**2+(posmarkz-prjz)**2 > disshowlabel**2) cycle
						if (abs(v1x*v2y-v2x*v1y)>1D-8) then
							det2_2=v1x*v2y-v2x*v1y
							n1=( (prjx-orgx2D)*v2y-v2x*(prjy-orgy2D) )/det2_2
							n2=( v1x*(prjy-orgy2D)-(prjx-orgx2D)*v1y )/det2_2
						else if (abs(v1x*v2z-v2x*v1z)>1D-8) then
							det2_2=v1x*v2z-v2x*v1z
							n1=( (prjx-orgx2D)*v2z-v2x*(prjz-orgz2D) )/det2_2
							n2=( v1x*(prjz-orgz2D)-(prjx-orgx2D)*v1z )/det2_2
						else if (abs(v1y*v2z-v2y*v1z)>1D-8) then
							det2_2=v1y*v2z-v2y*v1z
							n1=( (prjy-orgy2D)*v2z-v2y*(prjz-orgz2D) )/det2_2
							n2=( v1y*(prjz-orgz2D)-(prjy-orgy2D)*v1z )/det2_2
						end if
						tmp1=n1*d1; tmp2=n2*d2
						if (tmp1<init1.or.tmp1>end1.or.tmp2<init2.or.tmp2>end2) cycle !To avoid path out of plotting range
						if (iplaneoutall==0) then
							call rlsymb(21,tmp1,tmp2)
						else if (iplaneoutall==1) then
							write(10,"(2f12.6)") n1*d1*b2a,n2*d2*b2a
                        end if
					end if
				end do
				if (iplaneoutall==1) write(10,*) !Leave a blank line between each paths
			end do
		end do
		if (iplaneoutall==1) close(10)
	end if
	
	!Draw critical points on color-filled/contour/gradient/vector field graph. If iplaneoutall==1, then output the data points
	if (numcp>0) then
		if (iplaneoutall==0) then
			call HSYMBL(sizemarkcp)
		else if (iplaneoutall==1) then
			open(10,file="planeCP.txt",status="replace")
		end if
		do icp=1,numcp
			if (CPtype(icp)==1.and.imark3n3==0) cycle
			if (CPtype(icp)==2.and.imark3n1==0) cycle
			if (CPtype(icp)==3.and.imark3p1==0) cycle
			if (CPtype(icp)==4.and.imark3p3==0) cycle
			if (iplaneoutall==0) then
				if (CPtype(icp)==1) call setrgb(CP3n3RGB_2D(1),CP3n3RGB_2D(2),CP3n3RGB_2D(3))
				if (CPtype(icp)==2) call setrgb(CP3n1RGB_2D(1),CP3n1RGB_2D(2),CP3n1RGB_2D(3))
				if (CPtype(icp)==3) call setrgb(CP3p1RGB_2D(1),CP3p1RGB_2D(2),CP3p1RGB_2D(3))
				if (CPtype(icp)==4) call setrgb(CP3p3RGB_2D(1),CP3p3RGB_2D(2),CP3p3RGB_2D(3))
			end if
			posmarkx=CPpos(1,icp)*scll
			posmarky=CPpos(2,icp)*scll
			posmarkz=CPpos(3,icp)*scll
			if (plesel==1.and.abs(posmarkz-orgz2D)<disshowCP) then
				if (posmarkx<init1.or.posmarkx>end1.or.posmarky<init2.or.posmarky>end2) cycle !To avoid path out of plotting range
				if (iplaneoutall==0) call rlsymb(21,posmarkx,posmarky)
				if (iplaneoutall==1) write(10,"(2f12.6,i4)") posmarkx*b2a,posmarky*b2a,CPtype(icp)
			else if (plesel==2.and.abs(posmarky-orgy2D)<disshowCP) then
				if (posmarkx<init1.or.posmarkx>end1.or.posmarkz<init2.or.posmarkz>end2) cycle !To avoid path out of plotting range
				if (iplaneoutall==0) call rlsymb(21,posmarkx,posmarkz)
				if (iplaneoutall==1) write(10,"(2f12.6,i4)") posmarkx*b2a,posmarkz*b2a,CPtype(icp)
			else if (plesel==3.and.abs(posmarkx-orgx2D)<disshowCP) then
				if (posmarky<init1.or.posmarky>end1.or.posmarkz<init2.or.posmarkz>end2) cycle !To avoid path out of plotting range
				if (iplaneoutall==0) call rlsymb(21,posmarky,posmarkz)
				if (iplaneoutall==1) write(10,"(2f12.6,i4)") posmarky*b2a,posmarky*b2a,CPtype(icp)
			else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
				if (potpledis(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,posmarkx,posmarky,posmarkz)<disshowCP) then
					call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,posmarkx,posmarky,posmarkz,prjx,prjy,prjz)
					if (abs(v1x*v2y-v2x*v1y)>1D-8) then
						det2_2=v1x*v2y-v2x*v1y
						n1=( (prjx-orgx2D)*v2y-v2x*(prjy-orgy2D) )/det2_2
						n2=( v1x*(prjy-orgy2D)-(prjx-orgx2D)*v1y )/det2_2
					else if (abs(v1x*v2z-v2x*v1z)>1D-8) then
						det2_2=v1x*v2z-v2x*v1z
						n1=( (prjx-orgx2D)*v2z-v2x*(prjz-orgz2D) )/det2_2
						n2=( v1x*(prjz-orgz2D)-(prjx-orgx2D)*v1z )/det2_2
					else if (abs(v1y*v2z-v2y*v1z)>1D-8) then
						det2_2=v1y*v2z-v2y*v1z
						n1=( (prjy-orgy2D)*v2z-v2y*(prjz-orgz2D) )/det2_2
						n2=( v1y*(prjz-orgz2D)-(prjy-orgy2D)*v1z )/det2_2
					end if
					tmp1=n1*d1; tmp2=n2*d2
					if (tmp1<init1.or.tmp1>end1.or.tmp2<init2.or.tmp2>end2) cycle !To avoid path out of plotting range
					if (iplaneoutall==0) then
						call rlsymb(21,tmp1,tmp2)
					else if (iplaneoutall==1) then
						write(10,"(2f12.6,i4)") n1*d1*b2a,n2*d2*b2a,CPtype(icp)
                    end if
				end if
			end if
		end do
		if (iplaneoutall==1) close(10)
	end if
	if (iplaneoutall==1) return
    
    !Draw extrema of a function on a contour line of specific value of present function
    if (iextrema_on_contour==1) then
		call conpts(xcoord,ngridnum1,ycoord,ngridnum2,planetrunc,ctrval_2Dextrema,ctrptx,ctrpty,maxctrpt,ctrptnum,maxctrline,nctrlines)
        if (nctrlines==0) then
			write(*,"(a,1PE16.8,a)") " Warning: No contour line of isovalue of",ctrval_2Dextrema," can be generated!"
        else
			write(*,"(i4,a)") nctrlines," contour line(s) are generated"
			!Show coordinate of contour line points on the map
	   !     istart=0
	   !     do ictr=1,nctrlines
				!write(*,"(/,a,i5)") " Contour line",ictr
	   !         do ipt=istart+1,istart+ctrptnum(ictr)
					!write(*,"(i6,2f18.12)") ipt,ctrptx(ipt),ctrpty(ipt)
	   !         end do
	   !         istart=istart+ctrptnum(ictr)
	   !     end do
			!Determine if the contour lines are open
			istart=0
			do ictr=1,nctrlines
				dist=dsqrt((ctrptx(istart+1)-ctrptx(istart+ctrptnum(ictr)))**2 + (ctrpty(istart+1)-ctrpty(istart+ctrptnum(ictr)))**2)
				if (dist>dsqrt(dx2D**2+dy2D**2)) then
					ifopenctr(ictr)=1
					write(*,"(' Contour',i5,' is open')") ictr
				else
					ifopenctr(ictr)=0
					write(*,"(' Contour',i5,' is closed')") ictr
				end if
				istart=istart+ctrptnum(ictr)
			end do
			write(*,*) "Calculating function value on the contour lines..."
			istart=0
			do ictr=1,nctrlines
				!write(*,"(/,a,i5)") " Contour line",ictr
				!$OMP PARALLEL DO SHARED(ctrptval) PRIVATE(ipt,fract1,fract2,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads)
				do ipt=istart+1,istart+ctrptnum(ictr)
					!Transform 2D coordinate in the map to real space 3D coordinate, then calculate function value
					fract1=(ctrptx(ipt)-init1)/dx2D
					fract2=(ctrpty(ipt)-init2)/dy2D
					tmpx=orgx2D+fract1*v1x+fract2*v2x
					tmpy=orgy2D+fract1*v1y+fract2*v2y
					tmpz=orgz2D+fract1*v1z+fract2*v2z
					ctrptval(ipt)=calcfuncall(ifunc_2Dextrema,tmpx,tmpy,tmpz)
					!write(*,"(i6,3f16.10,1PE18.10)") ipt,tmpx,tmpy,tmpz,ctrptval(ipt)
					!write(10,"(3f16.10)") tmpx*b2a,tmpy*b2a,tmpz*b2a
				end do
				!$OMP END PARALLEL DO
				istart=istart+ctrptnum(ictr)
			end do
			ctrptval_org=ctrptval
			!Hint: If there is still artificial extrema, set ntime=3 may be helpful. If some extrema are missing, try to only consider nearest neighours (ignore considering secondary ones)
			write(*,*) "Detecting extrema on the contour lines..."
			istart=0
			do ictr=1,nctrlines
				write(*,"(/,a,i5)") " Contour line",ictr
                if (ifopenctr(ictr)==0) then
					ibeg=istart+1
					iend=istart+ctrptnum(ictr)
                else
					ibeg=istart+3
					iend=istart+ctrptnum(ictr)-2
                end if
				nmaxi=0
				nmini=0
				!Denoising by taking average between neighours, so that minor numerical fluctuations (causing artifical extrema) can be removed
				ntime=2 !I found denoising twice is fully adequate
				do itime=1,ntime
					ctrptval_tmp=ctrptval
					do ipt=ibeg,iend
						iprev=ipt-1
						if (iprev<ibeg) iprev=iend
						iprev2=iprev-1
						if (iprev2<ibeg) iprev2=iend
						inext=ipt+1
						if (inext>iend) inext=ibeg
						inext2=inext+1
						if (inext2>iend) inext2=ibeg
						!ctrptval(ipt)=(ctrptval_tmp(ipt)+ctrptval_tmp(iprev)+ctrptval_tmp(inext))/3
						ctrptval(ipt)=(ctrptval_tmp(ipt)+ctrptval_tmp(iprev)+ctrptval_tmp(iprev2)+ctrptval_tmp(inext)+ctrptval_tmp(inext2))/5
					end do
				end do
				!Check effect of denoising: Twice, once, no denoising
	   !         open(10,file="checksmooth.txt",status="replace")
				!do ipt=ibeg,iend
				!	write(10,"(i4,3f16.10)") ipt,ctrptval(ipt),ctrptval_tmp(ipt),ctrptval_org(ipt)
	   !         end do
	   !         close(10)
				do ipt=ibeg,iend
					!Consider both nearest and second nearest neighbours is important even after denoising data, otherwise there may be many meaningless extrema
					iprev=ipt-1
					if (iprev<ibeg) iprev=iend
					iprev2=iprev-1
					if (iprev2<ibeg) iprev2=iend
					inext=ipt+1
					if (inext>iend) inext=ibeg
					inext2=inext+1
					if (inext2>iend) inext2=ibeg
					if (ctrptval(ipt)>ctrptval(iprev).and.ctrptval(ipt)>ctrptval(inext).and.ctrptval(ipt)>ctrptval(iprev2).and.ctrptval(ipt)>ctrptval(inext2)) then
						write(*,"(a,2f12.6,'  Value:',f16.10)") ' Maximum found at',ctrptx(ipt),ctrpty(ipt),ctrptval_org(ipt)
						call color('red')
						call rlsymb(21,ctrptx(ipt),ctrpty(ipt))
						nmaxi=nmaxi+1
					else if (ctrptval(ipt)<ctrptval(iprev).and.ctrptval(ipt)<ctrptval(inext).and.ctrptval(ipt)<ctrptval(iprev2).and.ctrptval(ipt)<ctrptval(inext2)) then
						write(*,"(a,2f12.6,'  Value:',f16.10)") ' Minimum found at',ctrptx(ipt),ctrpty(ipt),ctrptval_org(ipt)
						call color('blue')
						call rlsymb(21,ctrptx(ipt),ctrpty(ipt))
						nmini=nmini+1
					end if
				end do
				write(*,"(' Total maxima:',i6)") nmaxi
				write(*,"(' Total minima:',i6)") nmini
				istart=istart+ctrptnum(ictr)
			end do
			write(*,*)
			write(*,"(a)") " On the map, red and blue spheres correspond to maxima and minima on the contour line, respectively"
			call color("WHITE") !Restore to default
        end if
    end if

    !Construct "inplane" list, if =1, the atom or reference point is close enough to the plotting plane
    !The "inplane" will be used for plotting atomic labels and bonds later
    inplane=0
    nallpoints=ncenter
    if (imarkrefpos==1) nallpoints=ncenter+1
    do ipt=1,nallpoints
	    if (ipt<=ncenter) then
		    posmarkx=a(ipt)%x*scll; posmarky=a(ipt)%y*scll; posmarkz=a(ipt)%z*scll
	    else if (ipt==ncenter+1) then !Plot reference point of correlation hole/factor, source function...
		    posmarkx=refx*scll; posmarky=refy*scll; posmarkz=refz*scll
	    end if
	    if (plesel==1) then
		    if (abs(posmarkz-orgz2D)<disshowlabel) inplane(ipt)=1
	    else if (plesel==2) then
		    if (abs(posmarky-orgy2D)<disshowlabel) inplane(ipt)=1
	    else if (plesel==3) then
		    if (abs(posmarkx-orgx2D)<disshowlabel) inplane(ipt)=1
	    else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
		    call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,posmarkx,posmarky,posmarkz,prjx,prjy,prjz)
		    if ((prjx-posmarkx)**2+(prjy-posmarky)**2+(prjz-posmarkz)**2<disshowlabel**2) inplane(ipt)=1
	    end if
    end do

    !Plot bonds as lines only for the atoms whose labels are shown on the plane
    if (ibond_on_plane==1) then
	    call setcolor(iclrindbndlab); call solid; call LINWID(bondthick2D)
	    do ipt=1,ncenter
		    posmarkx=a(ipt)%x*scll
		    posmarky=a(ipt)%y*scll
		    posmarkz=a(ipt)%z*scll
		    if (inplane(ipt)==0) cycle
		    if (plesel<=3) then !Cartesian plane
			    do iatm=ipt+1,ncenter
				    if (inplane(iatm)==0) cycle
				    if (atomdist(ipt,iatm,0) < ( covr(a(ipt)%index)+covr(a(iatm)%index) )*bondcrit) then
					    if (plesel==1) call rline(posmarkx,posmarky,a(iatm)%x*scll,a(iatm)%y*scll)
					    if (plesel==2) call rline(posmarkx,posmarkz,a(iatm)%x*scll,a(iatm)%z*scll)
					    if (plesel==3) call rline(posmarky,posmarkz,a(iatm)%y*scll,a(iatm)%z*scll)
				    end if
			    end do
		    else if (plesel<=7) then !Plane defined by three atoms/points
			    call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,posmarkx,posmarky,posmarkz,prjx,prjy,prjz)
			    !Get position of ipt in the plotting coordinate. Comment can be found in similar part below
			    if (abs(v1x*v2y-v2x*v1y)>1D-8) then
				    det2_2=v1x*v2y-v2x*v1y
				    n1=( (prjx-orgx2D)*v2y-v2x*(prjy-orgy2D) )/det2_2
				    n2=( v1x*(prjy-orgy2D)-(prjx-orgx2D)*v1y )/det2_2
			    else if (abs(v1x*v2z-v2x*v1z)>1D-8) then
				    det2_2=v1x*v2z-v2x*v1z
				    n1=( (prjx-orgx2D)*v2z-v2x*(prjz-orgz2D) )/det2_2
				    n2=( v1x*(prjz-orgz2D)-(prjx-orgx2D)*v1z )/det2_2
			    else if (abs(v1y*v2z-v2y*v1z)>1D-8) then
				    det2_2=v1y*v2z-v2y*v1z
				    n1=( (prjy-orgy2D)*v2z-v2y*(prjz-orgz2D) )/det2_2
				    n2=( v1y*(prjz-orgz2D)-(prjy-orgy2D)*v1z )/det2_2
			    end if
			    if ( n1>0.and.n2>0.and.n1<(ngridnum1-1).and.n2<(ngridnum2-1) ) then !The ipt is within the scope of drawing range
				    do iatm=ipt+1,ncenter
					    if (inplane(iatm)==0) cycle
					    if (atomdist(ipt,iatm,0) < ( covr(a(ipt)%index)+covr(a(iatm)%index) )*bondcrit) then
						    !Get position of iatm in the plotting coordinate
						    call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,a(iatm)%x*scll,a(iatm)%y*scll,a(iatm)%z*scll,prjx,prjy,prjz)
						    if (abs(v1x*v2y-v2x*v1y)>1D-8) then
							    det2_2=v1x*v2y-v2x*v1y
							    n1_2=( (prjx-orgx2D)*v2y-v2x*(prjy-orgy2D) )/det2_2
							    n2_2=( v1x*(prjy-orgy2D)-(prjx-orgx2D)*v1y )/det2_2
						    else if (abs(v1x*v2z-v2x*v1z)>1D-8) then
							    det2_2=v1x*v2z-v2x*v1z
							    n1_2=( (prjx-orgx2D)*v2z-v2x*(prjz-orgz2D) )/det2_2
							    n2_2=( v1x*(prjz-orgz2D)-(prjx-orgx2D)*v1z )/det2_2
						    else if (abs(v1y*v2z-v2y*v1z)>1D-8) then
							    det2_2=v1y*v2z-v2y*v1z
							    n1_2=( (prjy-orgy2D)*v2z-v2y*(prjz-orgz2D) )/det2_2
							    n2_2=( v1y*(prjz-orgz2D)-(prjy-orgy2D)*v1z )/det2_2
						    end if
						    if ( n1_2>0.and.n2_2>0.and.n1_2<(ngridnum1-1).and.n2_2<(ngridnum2-1) ) then
							    call rline(n1*d1,n2*d2,n1_2*d1,n2_2*d2)
						    end if
					    end if
				    end do
			    end if
		    end if
	    end do
	    CALL LINWID(1) !Restore to default
    end if

    !Draw atomic labels or reference point on graph
    pix2usr=(end1-init1)/lengthx  !Convert actual pixel to user coordinate
    if (iatom_on_plane==1) then
	    call height(pleatmlabsize)
	    movetext=pleatmlabsize/2 !Slight movement on the labels to make the center of label just appear at expected position
	    do ipt=1,nallpoints
 		    call SERIF
		    CALL SHDCHA
		    if (ipt<=ncenter) then !Plot atomic label on plane
			    posmarkx=a(ipt)%x*scll; posmarky=a(ipt)%y*scll; posmarkz=a(ipt)%z*scll
	    	    if (iatmlabtype==1) then !Only plot element for atomic label
	    	        atmlabtext=a(ipt)%name
	    	    else if (iatmlabtype==2.or.iatmlabtype==3) then !Plot index
	    	        write(atmlabtext,"(i6)") ipt
	    	        atmlabtext=adjustl(atmlabtext)
				    tmpstr=atmlabtext
	    	        if (iatmlabtype==3) write(atmlabtext,"(a,a)") trim(a(ipt)%name),trim(tmpstr) !Plot both element and index
	            end if
	            call setcolor(iclrindatmlab)
		    else !Plot reference point of correlation hole/factor, source function...
			    posmarkx=refx*scll; posmarky=refy*scll; posmarkz=refz*scll
			    call color('blue')
			    call HSYMBL(pleatmlabsize+20)
		    end if
		    if (plesel==1) then !XY plane
			    if (inplane(ipt)==1) then !Close enough to the plane
				    if (ipt<=ncenter) call rlmess(trim(atmlabtext),posmarkx-movetext*pix2usr,posmarky+movetext*pix2usr)
				    if (ipt==ncenter+1) call rlsymb(8,posmarkx,posmarky)
			    else if (iatom_on_plane_far==1) then !Far away from the plane
				    call DUPLX
				    call rlmess(trim(atmlabtext),posmarkx-movetext*pix2usr,posmarky+movetext*pix2usr)
			    end if
		    else if (plesel==2) then !XZ plane
			    if (inplane(ipt)==1) then
				    if (ipt<=ncenter) call rlmess(trim(atmlabtext),posmarkx-movetext*pix2usr,posmarkz+movetext*pix2usr)
				    if (ipt==ncenter+1) call rlsymb(8,posmarkx,posmarkz)
			    else if (iatom_on_plane_far==1) then
				    call DUPLX
				    call rlmess(trim(atmlabtext),posmarkx-movetext*pix2usr,posmarkz+movetext*pix2usr)
			    end if
		    else if (plesel==3) then !YZ plane
			    if (inplane(ipt)==1) then
				    if (ipt<=ncenter) call rlmess(trim(atmlabtext),posmarky-movetext*pix2usr,posmarkz+movetext*pix2usr)
				    if (ipt==ncenter+1) call rlsymb(8,posmarky,posmarkz)
			    else if (iatom_on_plane_far==1) then
				    call DUPLX
				    call rlmess(trim(atmlabtext),posmarky-movetext*pix2usr,posmarkz+movetext*pix2usr)
			    end if
		    else if (plesel==4.or.plesel==5.or.plesel==6.or.plesel==7) then
			    call pointprjple(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,posmarkx,posmarky,posmarkz,prjx,prjy,prjz)
			    !prjx-orgx2D=n1*v1x+n2*v2x, n1*d1 equals x in user coordinate. prjy as well.  prjx,y,z is projected position from posmark point
			    !prjy-orgy2D=n1*v1y+n2*v2y
			    !prjz-orgz2D=n1*v1z+n2*v2z
			    !Use Kramer rule to solve this linear equation to get n1 and n2
			    !We can use any two conditions, the precondition is det2_2 is not almost zero
			    if (abs(v1x*v2y-v2x*v1y)>1D-8) then
				    det2_2=v1x*v2y-v2x*v1y
				    n1=( (prjx-orgx2D)*v2y-v2x*(prjy-orgy2D) )/det2_2
				    n2=( v1x*(prjy-orgy2D)-(prjx-orgx2D)*v1y )/det2_2
			    else if (abs(v1x*v2z-v2x*v1z)>1D-8) then
				    det2_2=v1x*v2z-v2x*v1z
				    n1=( (prjx-orgx2D)*v2z-v2x*(prjz-orgz2D) )/det2_2
				    n2=( v1x*(prjz-orgz2D)-(prjx-orgx2D)*v1z )/det2_2
			    else if (abs(v1y*v2z-v2y*v1z)>1D-8) then
				    det2_2=v1y*v2z-v2y*v1z
				    n1=( (prjy-orgy2D)*v2z-v2y*(prjz-orgz2D) )/det2_2
				    n2=( v1y*(prjz-orgz2D)-(prjy-orgy2D)*v1z )/det2_2
			    end if
			    if (n1>0.and.n2>0.and.n1<(ngridnum1-1).and.n2<(ngridnum2-1)) then !Avoid atom label out of range
				    if (inplane(ipt)==1) then
					    if (ipt<=ncenter) call rlmess(trim(atmlabtext),n1*d1-movetext*pix2usr,n2*d2+movetext*pix2usr)
					    if (ipt==ncenter+1) call rlsymb(8,n1*d1,n2*d2) ! Note: rlsymb automatically avoids plotting symbol out of axis range
				    else if (iatom_on_plane_far==1) then
					    call DUPLX
					    call rlmess(trim(atmlabtext),n1*d1-movetext*pix2usr,n2*d2+movetext*pix2usr)
				    end if
			    end if
		    end if
	    end do
	    call color("WHITE") !Restore to default
    end if

!Relief map & shaded relief map with/without projection
else if (idrawtype==3.or.idrawtype==4.or.idrawtype==5) then
    call setcolortable(iclrtrans) !This routine must be invoked prior to GRAF
    CALL AXSPOS(100,2800) !Make position of coordinate proper
	CALL ORIGIN(ORIGIN_3D_X,ORIGIN_3D_Y) !Set origin of the map
	planetrunc=planemat
	planetrunc2=planemat
	!Now truncate the value in planemat to uplimit of Z-scale of relief map and save to planetrunc, else the color scale will range from
	!minimal to maximum, so the color transition is not completely between lower and upper limit of relief map, and effect is not good
	if (inucespplot==1) then !Now cut value, nuclear attraction potential is too big so treat it separately
		where (planetrunc>50) planetrunc=50
	else
! 		if (any(planetrunc>3).or.any(planetrunc<-3)) write(*,"(a)") " Note: The function values are temporarily truncated &
! 		at +/- 3.0 to garanteee that data range is completely within range of Z-axis"
		where (planetrunc>3)
			planetrunc=3
		elsewhere (planetrunc<-3)
			planetrunc=-3
		end where
	end if
	CALL AXSLEN(3000,3000) !Should not smaller than xxx of page(xxx,xxx), else project map couldn't show completely
    
	!!! Set axis
	call axis3D(2D0,(end2-init2)/(end1-init1)*2D0,2D0)
	if (idrawtype==5) then !Employ large negative part in Z to avoid relief map overlay the projected map
		if (inucespplot==1.and.ifiletype/=4) then !Nuclear ESP, since it is not from atomic charge (ifiletype==4), the value will be huge, use large Z uplimit
			CALL GRAF3D(init1,end1,init1-shiftx,planestpx, init2,end2,init2-shifty,planestpy, -50D0,50D0,-50D0,5D0)
		else !Common cases
			CALL GRAF3D(init1,end1,init1-shiftx,planestpx, init2,end2,init2-shifty,planestpy, -5D0,3D0,-5D0,1D0)
		end if
	else if (idrawtype==3.or.idrawtype==4) then
		if (inucespplot==1.and.ifiletype/=4) then !Nuclear ESP
			CALL GRAF3D(init1,end1,init1-shiftx,planestpx, init2,end2,init2-shifty,planestpy, -3D0,50D0,-3D0,5D0)
		else !Common cases
			CALL GRAF3D(init1,end1,init1-shiftx,planestpx, init2,end2,init2-shifty,planestpy, -3D0,3D0,-3D0,1D0)
		end if
	end if
	!!! Plot
	call light('on')
	if (idrawtype==3) then
! 		call surmat(planetrunc,ngridnum1,ngridnum2,1,1) !Works well for most functions, however not good for Laplacian
		call surmsh("only") !Draw grids on shaded surface
		CALL SHDMOD('SMOOTH','SURFACE')
		call surshd(xcoord,ngridnum1,ycoord,ngridnum2,planetrunc) !Draw shaded relief map
		call surmsh("OFF") !Recover to default, otherwise when drawing molecule structure later, black lines will present on object surfaces
	else if (idrawtype==4.or.idrawtype==5) then
		if (idrawtype==5) then !Draw projection
            !Note that it is not possible to use "call rlmess" to plot atomic label on the projected map, it only works for 2D map
            planetrunc2=planemat
            if (iclrtrans/=0) then !Truncate the values larger than and lower than color scale, so that these regions will not be shown as white and black, respectively
			    where (planetrunc2>end3) planetrunc2=end3-1D-10   !Augment by a minimal value to avoid numerical noise
			    where (planetrunc2<init3) planetrunc2=init3+1D-10
		    end if
			CALL GRFINI(-1D0,-(end2-init2)/(end1-init1),-1D0, 1D0,-(end2-init2)/(end1-init1),-1D0, 1D0,(end2-init2)/(end1-init1),-1D0)
			call SETGRF('none','none','none','none') !Do not show axis ticks around the plane map at bottom, since it has been show in 3D axis
			call AUTRES(ngridnum1,ngridnum2)
			call VKXBAR(170)
	        call height(67) !Slightly increase text size of color bar
			call GRAF3(init1,end1,init1-shiftx,planestpx, init2,end2,init2-shifty,planestpy, init3,end3,init3-shiftz,planestpz)
			call CRVMAT(planetrunc2,ngridnum1,ngridnum2,fillcoloritpx,fillcoloritpy)
			CALL GRFFIN
		end if
		call litmod(1,'on')
		call litmod(3,'on')
		call litpos(1,XVU,YVU,ZVU,'ANGLE')
		call surmsh(drawsurmesh) !Draw grids on shaded surface
		CALL SHDMOD('SMOOTH','SURFACE')
		call zscale(surcolorzmin,surcolorzmax) !For shaded relief map, set different color to shad defined by user
		call surshd(xcoord,ngridnum1,ycoord,ngridnum2,planetrunc)
		call surmsh("OFF") !Recover to default, otherwise when drawing molecule structure later, black lines will present on object surfaces
	end if
end if

call DISFIN

!Convert to original length unit
disshowlabel=disshowlabel/scll
disshowCP=disshowCP/scll
disshowpath=disshowpath/scll
if (ilenunit2D==2) call convgridlenunit(2)
end subroutine


end module



!---- Interface of selecting color by index
subroutine selcolor(clrind)
integer clrind
write(*,*) "1  = Red        2  = Green"
write(*,*) "3  = Blue       4  = White"
write(*,*) "5  = Black      6  = Gray"
write(*,*) "7  = Cyan       8  = Yellow"
write(*,*) "9  = Orange     10 = Magenta"
write(*,*) "11 = Crimson    12 = Dark green"
write(*,*) "13 = Purple     14 = Brown"
write(*,*) "15 = Dark blue  16 = Pink"
read(*,*) clrind
end subroutine


!---- Set color used by DISLIN routine by index
subroutine setcolor(clrind)
integer clrind
if (clrind==1) call color('RED')
if (clrind==2) call color('GREEN')
if (clrind==3) call color('BLUE')
if (clrind==4) call color('BLACK') !Due to current mode is REVERSE, so WHITE=BLACK
if (clrind==5) call color('WHITE')
if (clrind==6) call setRGB(0.65D0,0.65D0,0.7D0) !Gray
if (clrind==7) call color('CYAN')
if (clrind==8) call color('YELLOW')
if (clrind==9) call color('ORANGE')
if (clrind==10) call color('MAGENTA')
if (clrind==11) call setRGB(0.7D0,0D0,0D0) !Crimson
if (clrind==12) call setRGB(0D0,0.7D0,0D0) !Dark green
if (clrind==13) call setRGB(0.4D0,0D0,0.84D0) !Purple
if (clrind==14) call setRGB(0.7D0,0.5D0,0.4D0) !Brown
if (clrind==15) call setRGB(0D0,0D0,0.5D0) !Dark blue
if (clrind==16) call setRGB(1D0,0.5D0,1D0) !Pink
end subroutine


!------ Translate color index (defined by colorname) to R,G,B components
subroutine clridx2RGB(iclrind,Rcomp,Gcomp,Bcomp)
integer iclrind
real*8 Rcomp,Gcomp,Bcomp
Rcomp=0;Gcomp=0;Bcomp=0	!White
if (iclrind==1) then !Red
	Rcomp=1
else if (iclrind==2) then !Green
	Gcomp=1
else if (iclrind==3) then !Blue
	Bcomp=1
else if (iclrind==4) then !BLACK
	Rcomp=1;Gcomp=1;Bcomp=1
else if (iclrind==5) then !WHITE
	continue
else if (iclrind==6) then !Gray
	Rcomp=0.65D0;Gcomp=0.65D0;Bcomp=0.7D0
else if (iclrind==7) then !CYAN
	Gcomp=1D0;Bcomp=1D0
else if (iclrind==8) then !YELLOW
	Rcomp=1D0;Gcomp=1D0
else if (iclrind==9) then !ORANGE
	Rcomp=1D0;Gcomp=0.5D0
else if (iclrind==10) then !MAGENTA
	Rcomp=1D0;Bcomp=1D0
else if (iclrind==11) then !Crimson
	Rcomp=0.7D0
else if (iclrind==12) then !Dark green
	Gcomp=0.7D0
else if (iclrind==13) then !Purple
	Rcomp=0.4D0;Bcomp=0.84D0
else if (iclrind==14) then !Brown
	Rcomp=0.7D0;Gcomp=0.5D0;Bcomp=0.4D0
else if (iclrind==15) then !Dark blue
	Bcomp=0.5D0
else if (iclrind==16) then !Pink
	Rcomp=1D0;Gcomp=0.5D0;Bcomp=1D0
end if
end subroutine


subroutine setgraphformat
use defvar
use dislin
write(*,*) "Input index to select a format"
write(*,*) "Note: 1~4 are pixel formats, 5~9 are vector formats"
write(*,*) "1 png"
write(*,*) "2 gif"
write(*,*) "3 tiff"
write(*,*) "4 bmp"
write(*,*) "5 ps"
write(*,*) "6 eps"
write(*,*) "7 pdf"
write(*,*) "8 wmf"
write(*,*) "9 svg"
read(*,*) itmp
if (itmp==1) graphformat="png"
if (itmp==2) graphformat="gif"
if (itmp==3) graphformat="tiff"
if (itmp==4) graphformat="bmp"
if (itmp==5) graphformat="ps"
if (itmp==6) graphformat="eps"
if (itmp==7) graphformat="pdf"
if (itmp==8) graphformat="wmf"
if (itmp==9) graphformat="svg"
call setfil("dislin."//trim(graphformat)) !The file name of saved image file may have been modified, recover to default one
end subroutine


!!----------- Select color table
subroutine selcolortable
use defvar
write(*,*) "Select a color transition method"
write(*,*) "0  Rainbow with black/white for values exceeding lower/higher color limit"
write(*,*) "1  Rainbow               2 Reversed rainbow"
write(*,*) "3  Rainbow starting from white"
write(*,*) "4  Spectrum (Pink-Blue-Green-Red)  5 Reversed Spectrum"
write(*,*) "6  Grey (Black-White)    7  Reversed Grey"
write(*,*) "8  Blue-White-Red        9  Red-White-Blue"
write(*,*) "10 Blue-Green-Red        11 Red-Green-Blue"
write(*,*) "12 White-Dark red        13 Black-Orange-Yellow"
write(*,*) "14 White-Dark green      15 Black-Green"
write(*,*) "16 White-Dark blue       17 Black-Blue-Cyan"
write(*,*) "18 Viridis               19 Yellow-Orange-Black"
read(*,*) iclrtrans
end subroutine

!!----------- Set color table of Dislin, will affect color transition of color-filled map and heat map
subroutine setcolortable(isel) !"isel" corresponds to global variable "iclrtrans"
use dislin
implicit real*8 (a-h,o-z)
integer isel
integer,parameter :: nlevel=255
real*8 Rarr(0:nlevel),Garr(0:nlevel),Barr(0:nlevel),tmparr(0:nlevel)
Rarr=0
Garr=0
Barr=0
if (isel==0.or.isel==1) then !Rainbow. Note that for isel>0, in the plot.f90, the plane data will be automatically truncated before invoking plotting routines
    CALL SETVLT("RAIN")
else if (isel==2) then !Reversed rainbow
    CALL SETVLT("RRAIN")
else if (isel==4) then !Spectrum
    CALL SETVLT("SPEC")
else if (isel==5) then !Reversed spectrum
    CALL SETVLT("RSPEC")
else if (isel==6) then !Black-White
    !CALL SETVLT("GREY") !This will also make label color become grey
    do i=0,nlevel
        Rarr(i)=dfloat(i)/nlevel
        Garr(i)=dfloat(i)/nlevel
        Barr(i)=dfloat(i)/nlevel
    end do
else if (isel==7) then !White-Black
    !CALL SETVLT("RGREY")
    do i=0,nlevel
        Rarr(i)=dfloat(nlevel-i)/nlevel
        Garr(i)=dfloat(nlevel-i)/nlevel
        Barr(i)=dfloat(nlevel-i)/nlevel
    end do
else if (isel==3) then !White-Rainbow
    !Below setting makes appear of blue to early, leading to visually detectable color levels during white-blue transition
    !!0: White
    !!51: Blue(0,0,1)
    !!102: Cyan(0,1,1)
    !!153: Green(0,1,0)
    !!204: Yellow(1,1,0)
    !!255: Red(1,0,0)
    !Rarr(0)=1
    !Garr(0)=1
    !Barr(0)=1
    !do i=1,51 !White(1,1,1)-Blue(0,0,1)
    !    Rarr(i)=dfloat(51-i)/51
    !    Garr(i)=dfloat(51-i)/51
    !    Barr(i)=1
    !end do
    !j=0
    !do i=52,102 !Blue(0,0,1)-Cyan(0,1,1)
    !    j=j+1
    !    Rarr(i)=0
    !    Garr(i)=dfloat(j)/51
    !    Barr(i)=1
    !end do
    !j=0
    !do i=103,153 !Cyan(0,1,1)-Green(0,1,0)
    !    j=j+1
    !    Rarr(i)=0
    !    Garr(i)=1
    !    Barr(i)=dfloat(51-j)/51
    !end do
    !j=0
    !do i=154,204 !Green(0,1,0)-Yellow(1,1,0)
    !    j=j+1
    !    Rarr(i)=dfloat(j)/51
    !    Garr(i)=1
    !    Barr(i)=0
    !end do
    !j=0
    !do i=205,255 !Yellow(1,1,0)-Red(1,0,0)
    !    j=j+1
    !    Rarr(i)=1
    !    Garr(i)=dfloat(51-j)/51
    !    Barr(i)=0
    !end do
    !0: White
    !91: Blue(0,0,1)
    !132: Cyan(0,1,1)
    !173: Green(0,1,0)
    !214: Yellow(1,1,0)
    !255: Red(1,0,0)
    Rarr(0)=1
    Garr(0)=1
    Barr(0)=1
    do i=1,91 !White(1,1,1)-Blue(0,0,1)
        Rarr(i)=dfloat(91-i)/91
        Garr(i)=dfloat(91-i)/91
        Barr(i)=1
    end do
    j=0
    do i=92,132 !Blue(0,0,1)-Cyan(0,1,1)
        j=j+1
        Rarr(i)=0
        Garr(i)=dfloat(j)/41
        Barr(i)=1
    end do
    j=0
    do i=133,173 !Cyan(0,1,1)-Green(0,1,0)
        j=j+1
        Rarr(i)=0
        Garr(i)=1
        Barr(i)=dfloat(41-j)/41
    end do
    j=0
    do i=174,214 !Green(0,1,0)-Yellow(1,1,0)
        j=j+1
        Rarr(i)=dfloat(j)/41
        Garr(i)=1
        Barr(i)=0
    end do
    j=0
    do i=215,255 !Yellow(1,1,0)-Red(1,0,0)
        j=j+1
        Rarr(i)=1
        Garr(i)=dfloat(41-j)/41
        Barr(i)=0
    end do
else if (isel==8) then !Blue-White-Red
    !First half, Blue-White
    j=0
    do i=0,127
        j=j+1
        Rarr(i)=dfloat(2*j)/256
        Garr(i)=dfloat(2*j)/256
        Barr(i)=1
    end do
    !Second half, White-Red
    j=0
    do i=128,255
        j=j+1
        Rarr(i)=1
        Garr(i)=dfloat(256-2*j)/256
        Barr(i)=dfloat(256-2*j)/256
    end do
else if (isel==9) then !Red-White-Blue
    !First half, Red-White
    j=0
    do i=0,127
        j=j+1
        Rarr(i)=1
        Garr(i)=dfloat(2*j)/256
        Barr(i)=dfloat(2*j)/256
    end do
    !Second half, White-Blue
    j=0
    do i=128,255
        j=j+1
        Rarr(i)=dfloat(256-2*j)/256
        Garr(i)=dfloat(256-2*j)/256
        Barr(i)=1
    end do
else if (isel==10) then !Blue-Green-Red
    !Blue(0,0,1)-Green(0,1,0)
    j=0
    do i=0,127
        j=j+1
        Rarr(i)=0
        Garr(i)=dfloat(2*j)/256
        Barr(i)=dfloat(256-2*j)/256
    end do
    !Green(0,1,0)-Red(1,0,0)
    j=0
    do i=128,255
        j=j+1
        Rarr(i)=dfloat(2*j)/256
        Garr(i)=dfloat(256-2*j)/256
        Barr(i)=0
    end do
else if (isel==11) then !Red-Green-Blue
    !Red(1,0,0)-Green(0,1,0)
    j=0
    do i=0,127
        j=j+1
        Rarr(i)=dfloat(256-2*j)/256
        Garr(i)=dfloat(2*j)/256
        Barr(i)=0
    end do
    !Green(0,1,0)-Blue(0,0,1)
    j=0
    do i=128,255
        j=j+1
        Rarr(i)=0
        Garr(i)=dfloat(256-2*j)/256
        Barr(i)=dfloat(2*j)/256
    end do
else if (isel==12) then !White-Dark red
    do i=0,nlevel
        Rarr(i)=1-dfloat(i)/nlevel*0.3
        Garr(i)=dfloat(nlevel-i)/nlevel
        Barr(i)=dfloat(nlevel-i)/nlevel
    end do
else if (isel==13.or.isel==19) then !Black-Orange-Yellow (mimic black-hole map) and its invert
    do i=0,170
        Rarr(i)=dfloat(i)/170
    end do
    Rarr(171:255)=1
    if (isel==19) then
		tmparr=Rarr
        do i=0,255
			Rarr(i)=tmparr(255-i)
        end do
    end if
    do i=70,nlevel
        Garr(i)=dfloat(i-70)/185
    end do
    if (isel==19) then
		tmparr=Garr
        do i=0,255
			Garr(i)=tmparr(255-i)
        end do
    end if
    do i=155,nlevel
        Barr(i)=dfloat(i-155)/100*0.6D0
    end do
    if (isel==19) then
		tmparr=Barr
        do i=0,255
			Barr(i)=tmparr(255-i)
        end do
    end if
else if (isel==14) then !White-Dark green
    do i=0,nlevel
        Rarr(i)=dfloat(nlevel-i)/nlevel
        Garr(i)=1-dfloat(i)/nlevel*0.3
        Barr(i)=dfloat(nlevel-i)/nlevel
    end do
else if (isel==15) then !Black-Green
    do i=0,nlevel
        Garr(i)=dfloat(i)/nlevel
    end do
else if (isel==16) then !White-Dark blue
    do i=0,nlevel
        Rarr(i)=dfloat(nlevel-i)/nlevel
        Garr(i)=dfloat(nlevel-i)/nlevel
        Barr(i)=1-dfloat(i)/nlevel*0.2
    end do
else if (isel==17) then !Black-Blue-Cyan
    do i=0,200
        Barr(i)=dfloat(i)/200
    end do
    Barr(201:255)=1
    do i=100,nlevel
        Garr(i)=dfloat(i-100)/155
    end do
else if (isel==18) then !Viridis, come from https://raw.githubusercontent.com/Gnuplotting/gnuplot-palettes/master/viridis.pal
	call defineViridis(Rarr,Garr,Barr)
end if
if ((isel>=6.and.isel<=19).or.isel==3) CALL MYVLT(Rarr,Garr,Barr,nlevel)
end subroutine



!!-------- For a color transition (icolor), input a percent (0~1), return red, green and blue components (0~1)
!icolor=1: Viridis
subroutine percent2RGB(icolor,percent,Rval,Gval,Bval)
integer icolor
real*8 percent,Rval,Gval,Bval,Rarr(0:255),Garr(0:255),Barr(0:255)
if (icolor==1) then
	call defineViridis(Rarr,Garr,Barr)
end if
Rval=Rarr(nint(percent*255))
Gval=Garr(nint(percent*255))
Bval=Barr(nint(percent*255))
end subroutine



!!---------- Define RGB arrays of Viridis color transition
subroutine defineViridis(Rarr,Garr,Barr)
real*8 Rarr(0:255),Garr(0:255),Barr(0:255)
Rarr(0  )=0.267004D0;Garr(0  )=0.004874D0;Barr(0  )=0.329415D0
Rarr(1  )=0.268510D0;Garr(1  )=0.009605D0;Barr(1  )=0.335427D0
Rarr(2  )=0.269944D0;Garr(2  )=0.014625D0;Barr(2  )=0.341379D0
Rarr(3  )=0.271305D0;Garr(3  )=0.019942D0;Barr(3  )=0.347269D0
Rarr(4  )=0.272594D0;Garr(4  )=0.025563D0;Barr(4  )=0.353093D0
Rarr(5  )=0.273809D0;Garr(5  )=0.031497D0;Barr(5  )=0.358853D0
Rarr(6  )=0.274952D0;Garr(6  )=0.037752D0;Barr(6  )=0.364543D0
Rarr(7  )=0.276022D0;Garr(7  )=0.044167D0;Barr(7  )=0.370164D0
Rarr(8  )=0.277018D0;Garr(8  )=0.050344D0;Barr(8  )=0.375715D0
Rarr(9  )=0.277941D0;Garr(9  )=0.056324D0;Barr(9  )=0.381191D0
Rarr(10 )=0.278791D0;Garr(10 )=0.062145D0;Barr(10 )=0.386592D0
Rarr(11 )=0.279566D0;Garr(11 )=0.067836D0;Barr(11 )=0.391917D0
Rarr(12 )=0.280267D0;Garr(12 )=0.073417D0;Barr(12 )=0.397163D0
Rarr(13 )=0.280894D0;Garr(13 )=0.078907D0;Barr(13 )=0.402329D0
Rarr(14 )=0.281446D0;Garr(14 )=0.084320D0;Barr(14 )=0.407414D0
Rarr(15 )=0.281924D0;Garr(15 )=0.089666D0;Barr(15 )=0.412415D0
Rarr(16 )=0.282327D0;Garr(16 )=0.094955D0;Barr(16 )=0.417331D0
Rarr(17 )=0.282656D0;Garr(17 )=0.100196D0;Barr(17 )=0.422160D0
Rarr(18 )=0.282910D0;Garr(18 )=0.105393D0;Barr(18 )=0.426902D0
Rarr(19 )=0.283091D0;Garr(19 )=0.110553D0;Barr(19 )=0.431554D0
Rarr(20 )=0.283197D0;Garr(20 )=0.115680D0;Barr(20 )=0.436115D0
Rarr(21 )=0.283229D0;Garr(21 )=0.120777D0;Barr(21 )=0.440584D0
Rarr(22 )=0.283187D0;Garr(22 )=0.125848D0;Barr(22 )=0.444960D0
Rarr(23 )=0.283072D0;Garr(23 )=0.130895D0;Barr(23 )=0.449241D0
Rarr(24 )=0.282884D0;Garr(24 )=0.135920D0;Barr(24 )=0.453427D0
Rarr(25 )=0.282623D0;Garr(25 )=0.140926D0;Barr(25 )=0.457517D0
Rarr(26 )=0.282290D0;Garr(26 )=0.145912D0;Barr(26 )=0.461510D0
Rarr(27 )=0.281887D0;Garr(27 )=0.150881D0;Barr(27 )=0.465405D0
Rarr(28 )=0.281412D0;Garr(28 )=0.155834D0;Barr(28 )=0.469201D0
Rarr(29 )=0.280868D0;Garr(29 )=0.160771D0;Barr(29 )=0.472899D0
Rarr(30 )=0.280255D0;Garr(30 )=0.165693D0;Barr(30 )=0.476498D0
Rarr(31 )=0.279574D0;Garr(31 )=0.170599D0;Barr(31 )=0.479997D0
Rarr(32 )=0.278826D0;Garr(32 )=0.175490D0;Barr(32 )=0.483397D0
Rarr(33 )=0.278012D0;Garr(33 )=0.180367D0;Barr(33 )=0.486697D0
Rarr(34 )=0.277134D0;Garr(34 )=0.185228D0;Barr(34 )=0.489898D0
Rarr(35 )=0.276194D0;Garr(35 )=0.190074D0;Barr(35 )=0.493001D0
Rarr(36 )=0.275191D0;Garr(36 )=0.194905D0;Barr(36 )=0.496005D0
Rarr(37 )=0.274128D0;Garr(37 )=0.199721D0;Barr(37 )=0.498911D0
Rarr(38 )=0.273006D0;Garr(38 )=0.204520D0;Barr(38 )=0.501721D0
Rarr(39 )=0.271828D0;Garr(39 )=0.209303D0;Barr(39 )=0.504434D0
Rarr(40 )=0.270595D0;Garr(40 )=0.214069D0;Barr(40 )=0.507052D0
Rarr(41 )=0.269308D0;Garr(41 )=0.218818D0;Barr(41 )=0.509577D0
Rarr(42 )=0.267968D0;Garr(42 )=0.223549D0;Barr(42 )=0.512008D0
Rarr(43 )=0.266580D0;Garr(43 )=0.228262D0;Barr(43 )=0.514349D0
Rarr(44 )=0.265145D0;Garr(44 )=0.232956D0;Barr(44 )=0.516599D0
Rarr(45 )=0.263663D0;Garr(45 )=0.237631D0;Barr(45 )=0.518762D0
Rarr(46 )=0.262138D0;Garr(46 )=0.242286D0;Barr(46 )=0.520837D0
Rarr(47 )=0.260571D0;Garr(47 )=0.246922D0;Barr(47 )=0.522828D0
Rarr(48 )=0.258965D0;Garr(48 )=0.251537D0;Barr(48 )=0.524736D0
Rarr(49 )=0.257322D0;Garr(49 )=0.256130D0;Barr(49 )=0.526563D0
Rarr(50 )=0.255645D0;Garr(50 )=0.260703D0;Barr(50 )=0.528312D0
Rarr(51 )=0.253935D0;Garr(51 )=0.265254D0;Barr(51 )=0.529983D0
Rarr(52 )=0.252194D0;Garr(52 )=0.269783D0;Barr(52 )=0.531579D0
Rarr(53 )=0.250425D0;Garr(53 )=0.274290D0;Barr(53 )=0.533103D0
Rarr(54 )=0.248629D0;Garr(54 )=0.278775D0;Barr(54 )=0.534556D0
Rarr(55 )=0.246811D0;Garr(55 )=0.283237D0;Barr(55 )=0.535941D0
Rarr(56 )=0.244972D0;Garr(56 )=0.287675D0;Barr(56 )=0.537260D0
Rarr(57 )=0.243113D0;Garr(57 )=0.292092D0;Barr(57 )=0.538516D0
Rarr(58 )=0.241237D0;Garr(58 )=0.296485D0;Barr(58 )=0.539709D0
Rarr(59 )=0.239346D0;Garr(59 )=0.300855D0;Barr(59 )=0.540844D0
Rarr(60 )=0.237441D0;Garr(60 )=0.305202D0;Barr(60 )=0.541921D0
Rarr(61 )=0.235526D0;Garr(61 )=0.309527D0;Barr(61 )=0.542944D0
Rarr(62 )=0.233603D0;Garr(62 )=0.313828D0;Barr(62 )=0.543914D0
Rarr(63 )=0.231674D0;Garr(63 )=0.318106D0;Barr(63 )=0.544834D0
Rarr(64 )=0.229739D0;Garr(64 )=0.322361D0;Barr(64 )=0.545706D0
Rarr(65 )=0.227802D0;Garr(65 )=0.326594D0;Barr(65 )=0.546532D0
Rarr(66 )=0.225863D0;Garr(66 )=0.330805D0;Barr(66 )=0.547314D0
Rarr(67 )=0.223925D0;Garr(67 )=0.334994D0;Barr(67 )=0.548053D0
Rarr(68 )=0.221989D0;Garr(68 )=0.339161D0;Barr(68 )=0.548752D0
Rarr(69 )=0.220057D0;Garr(69 )=0.343307D0;Barr(69 )=0.549413D0
Rarr(70 )=0.218130D0;Garr(70 )=0.347432D0;Barr(70 )=0.550038D0
Rarr(71 )=0.216210D0;Garr(71 )=0.351535D0;Barr(71 )=0.550627D0
Rarr(72 )=0.214298D0;Garr(72 )=0.355619D0;Barr(72 )=0.551184D0
Rarr(73 )=0.212395D0;Garr(73 )=0.359683D0;Barr(73 )=0.551710D0
Rarr(74 )=0.210503D0;Garr(74 )=0.363727D0;Barr(74 )=0.552206D0
Rarr(75 )=0.208623D0;Garr(75 )=0.367752D0;Barr(75 )=0.552675D0
Rarr(76 )=0.206756D0;Garr(76 )=0.371758D0;Barr(76 )=0.553117D0
Rarr(77 )=0.204903D0;Garr(77 )=0.375746D0;Barr(77 )=0.553533D0
Rarr(78 )=0.203063D0;Garr(78 )=0.379716D0;Barr(78 )=0.553925D0
Rarr(79 )=0.201239D0;Garr(79 )=0.383670D0;Barr(79 )=0.554294D0
Rarr(80 )=0.199430D0;Garr(80 )=0.387607D0;Barr(80 )=0.554642D0
Rarr(81 )=0.197636D0;Garr(81 )=0.391528D0;Barr(81 )=0.554969D0
Rarr(82 )=0.195860D0;Garr(82 )=0.395433D0;Barr(82 )=0.555276D0
Rarr(83 )=0.194100D0;Garr(83 )=0.399323D0;Barr(83 )=0.555565D0
Rarr(84 )=0.192357D0;Garr(84 )=0.403199D0;Barr(84 )=0.555836D0
Rarr(85 )=0.190631D0;Garr(85 )=0.407061D0;Barr(85 )=0.556089D0
Rarr(86 )=0.188923D0;Garr(86 )=0.410910D0;Barr(86 )=0.556326D0
Rarr(87 )=0.187231D0;Garr(87 )=0.414746D0;Barr(87 )=0.556547D0
Rarr(88 )=0.185556D0;Garr(88 )=0.418570D0;Barr(88 )=0.556753D0
Rarr(89 )=0.183898D0;Garr(89 )=0.422383D0;Barr(89 )=0.556944D0
Rarr(90 )=0.182256D0;Garr(90 )=0.426184D0;Barr(90 )=0.557120D0
Rarr(91 )=0.180629D0;Garr(91 )=0.429975D0;Barr(91 )=0.557282D0
Rarr(92 )=0.179019D0;Garr(92 )=0.433756D0;Barr(92 )=0.557430D0
Rarr(93 )=0.177423D0;Garr(93 )=0.437527D0;Barr(93 )=0.557565D0
Rarr(94 )=0.175841D0;Garr(94 )=0.441290D0;Barr(94 )=0.557685D0
Rarr(95 )=0.174274D0;Garr(95 )=0.445044D0;Barr(95 )=0.557792D0
Rarr(96 )=0.172719D0;Garr(96 )=0.448791D0;Barr(96 )=0.557885D0
Rarr(97 )=0.171176D0;Garr(97 )=0.452530D0;Barr(97 )=0.557965D0
Rarr(98 )=0.169646D0;Garr(98 )=0.456262D0;Barr(98 )=0.558030D0
Rarr(99 )=0.168126D0;Garr(99 )=0.459988D0;Barr(99 )=0.558082D0
Rarr(100)=0.166617D0;Garr(100)=0.463708D0;Barr(100)=0.558119D0
Rarr(101)=0.165117D0;Garr(101)=0.467423D0;Barr(101)=0.558141D0
Rarr(102)=0.163625D0;Garr(102)=0.471133D0;Barr(102)=0.558148D0
Rarr(103)=0.162142D0;Garr(103)=0.474838D0;Barr(103)=0.558140D0
Rarr(104)=0.160665D0;Garr(104)=0.478540D0;Barr(104)=0.558115D0
Rarr(105)=0.159194D0;Garr(105)=0.482237D0;Barr(105)=0.558073D0
Rarr(106)=0.157729D0;Garr(106)=0.485932D0;Barr(106)=0.558013D0
Rarr(107)=0.156270D0;Garr(107)=0.489624D0;Barr(107)=0.557936D0
Rarr(108)=0.154815D0;Garr(108)=0.493313D0;Barr(108)=0.557840D0
Rarr(109)=0.153364D0;Garr(109)=0.497000D0;Barr(109)=0.557724D0
Rarr(110)=0.151918D0;Garr(110)=0.500685D0;Barr(110)=0.557587D0
Rarr(111)=0.150476D0;Garr(111)=0.504369D0;Barr(111)=0.557430D0
Rarr(112)=0.149039D0;Garr(112)=0.508051D0;Barr(112)=0.557250D0
Rarr(113)=0.147607D0;Garr(113)=0.511733D0;Barr(113)=0.557049D0
Rarr(114)=0.146180D0;Garr(114)=0.515413D0;Barr(114)=0.556823D0
Rarr(115)=0.144759D0;Garr(115)=0.519093D0;Barr(115)=0.556572D0
Rarr(116)=0.143343D0;Garr(116)=0.522773D0;Barr(116)=0.556295D0
Rarr(117)=0.141935D0;Garr(117)=0.526453D0;Barr(117)=0.555991D0
Rarr(118)=0.140536D0;Garr(118)=0.530132D0;Barr(118)=0.555659D0
Rarr(119)=0.139147D0;Garr(119)=0.533812D0;Barr(119)=0.555298D0
Rarr(120)=0.137770D0;Garr(120)=0.537492D0;Barr(120)=0.554906D0
Rarr(121)=0.136408D0;Garr(121)=0.541173D0;Barr(121)=0.554483D0
Rarr(122)=0.135066D0;Garr(122)=0.544853D0;Barr(122)=0.554029D0
Rarr(123)=0.133743D0;Garr(123)=0.548535D0;Barr(123)=0.553541D0
Rarr(124)=0.132444D0;Garr(124)=0.552216D0;Barr(124)=0.553018D0
Rarr(125)=0.131172D0;Garr(125)=0.555899D0;Barr(125)=0.552459D0
Rarr(126)=0.129933D0;Garr(126)=0.559582D0;Barr(126)=0.551864D0
Rarr(127)=0.128729D0;Garr(127)=0.563265D0;Barr(127)=0.551229D0
Rarr(128)=0.127568D0;Garr(128)=0.566949D0;Barr(128)=0.550556D0
Rarr(129)=0.126453D0;Garr(129)=0.570633D0;Barr(129)=0.549841D0
Rarr(130)=0.125394D0;Garr(130)=0.574318D0;Barr(130)=0.549086D0
Rarr(131)=0.124395D0;Garr(131)=0.578002D0;Barr(131)=0.548287D0
Rarr(132)=0.123463D0;Garr(132)=0.581687D0;Barr(132)=0.547445D0
Rarr(133)=0.122606D0;Garr(133)=0.585371D0;Barr(133)=0.546557D0
Rarr(134)=0.121831D0;Garr(134)=0.589055D0;Barr(134)=0.545623D0
Rarr(135)=0.121148D0;Garr(135)=0.592739D0;Barr(135)=0.544641D0
Rarr(136)=0.120565D0;Garr(136)=0.596422D0;Barr(136)=0.543611D0
Rarr(137)=0.120092D0;Garr(137)=0.600104D0;Barr(137)=0.542530D0
Rarr(138)=0.119738D0;Garr(138)=0.603785D0;Barr(138)=0.541400D0
Rarr(139)=0.119512D0;Garr(139)=0.607464D0;Barr(139)=0.540218D0
Rarr(140)=0.119423D0;Garr(140)=0.611141D0;Barr(140)=0.538982D0
Rarr(141)=0.119483D0;Garr(141)=0.614817D0;Barr(141)=0.537692D0
Rarr(142)=0.119699D0;Garr(142)=0.618490D0;Barr(142)=0.536347D0
Rarr(143)=0.120081D0;Garr(143)=0.622161D0;Barr(143)=0.534946D0
Rarr(144)=0.120638D0;Garr(144)=0.625828D0;Barr(144)=0.533488D0
Rarr(145)=0.121380D0;Garr(145)=0.629492D0;Barr(145)=0.531973D0
Rarr(146)=0.122312D0;Garr(146)=0.633153D0;Barr(146)=0.530398D0
Rarr(147)=0.123444D0;Garr(147)=0.636809D0;Barr(147)=0.528763D0
Rarr(148)=0.124780D0;Garr(148)=0.640461D0;Barr(148)=0.527068D0
Rarr(149)=0.126326D0;Garr(149)=0.644107D0;Barr(149)=0.525311D0
Rarr(150)=0.128087D0;Garr(150)=0.647749D0;Barr(150)=0.523491D0
Rarr(151)=0.130067D0;Garr(151)=0.651384D0;Barr(151)=0.521608D0
Rarr(152)=0.132268D0;Garr(152)=0.655014D0;Barr(152)=0.519661D0
Rarr(153)=0.134692D0;Garr(153)=0.658636D0;Barr(153)=0.517649D0
Rarr(154)=0.137339D0;Garr(154)=0.662252D0;Barr(154)=0.515571D0
Rarr(155)=0.140210D0;Garr(155)=0.665859D0;Barr(155)=0.513427D0
Rarr(156)=0.143303D0;Garr(156)=0.669459D0;Barr(156)=0.511215D0
Rarr(157)=0.146616D0;Garr(157)=0.673050D0;Barr(157)=0.508936D0
Rarr(158)=0.150148D0;Garr(158)=0.676631D0;Barr(158)=0.506589D0
Rarr(159)=0.153894D0;Garr(159)=0.680203D0;Barr(159)=0.504172D0
Rarr(160)=0.157851D0;Garr(160)=0.683765D0;Barr(160)=0.501686D0
Rarr(161)=0.162016D0;Garr(161)=0.687316D0;Barr(161)=0.499129D0
Rarr(162)=0.166383D0;Garr(162)=0.690856D0;Barr(162)=0.496502D0
Rarr(163)=0.170948D0;Garr(163)=0.694384D0;Barr(163)=0.493803D0
Rarr(164)=0.175707D0;Garr(164)=0.697900D0;Barr(164)=0.491033D0
Rarr(165)=0.180653D0;Garr(165)=0.701402D0;Barr(165)=0.488189D0
Rarr(166)=0.185783D0;Garr(166)=0.704891D0;Barr(166)=0.485273D0
Rarr(167)=0.191090D0;Garr(167)=0.708366D0;Barr(167)=0.482284D0
Rarr(168)=0.196571D0;Garr(168)=0.711827D0;Barr(168)=0.479221D0
Rarr(169)=0.202219D0;Garr(169)=0.715272D0;Barr(169)=0.476084D0
Rarr(170)=0.208030D0;Garr(170)=0.718701D0;Barr(170)=0.472873D0
Rarr(171)=0.214000D0;Garr(171)=0.722114D0;Barr(171)=0.469588D0
Rarr(172)=0.220124D0;Garr(172)=0.725509D0;Barr(172)=0.466226D0
Rarr(173)=0.226397D0;Garr(173)=0.728888D0;Barr(173)=0.462789D0
Rarr(174)=0.232815D0;Garr(174)=0.732247D0;Barr(174)=0.459277D0
Rarr(175)=0.239374D0;Garr(175)=0.735588D0;Barr(175)=0.455688D0
Rarr(176)=0.246070D0;Garr(176)=0.738910D0;Barr(176)=0.452024D0
Rarr(177)=0.252899D0;Garr(177)=0.742211D0;Barr(177)=0.448284D0
Rarr(178)=0.259857D0;Garr(178)=0.745492D0;Barr(178)=0.444467D0
Rarr(179)=0.266941D0;Garr(179)=0.748751D0;Barr(179)=0.440573D0
Rarr(180)=0.274149D0;Garr(180)=0.751988D0;Barr(180)=0.436601D0
Rarr(181)=0.281477D0;Garr(181)=0.755203D0;Barr(181)=0.432552D0
Rarr(182)=0.288921D0;Garr(182)=0.758394D0;Barr(182)=0.428426D0
Rarr(183)=0.296479D0;Garr(183)=0.761561D0;Barr(183)=0.424223D0
Rarr(184)=0.304148D0;Garr(184)=0.764704D0;Barr(184)=0.419943D0
Rarr(185)=0.311925D0;Garr(185)=0.767822D0;Barr(185)=0.415586D0
Rarr(186)=0.319809D0;Garr(186)=0.770914D0;Barr(186)=0.411152D0
Rarr(187)=0.327796D0;Garr(187)=0.773980D0;Barr(187)=0.406640D0
Rarr(188)=0.335885D0;Garr(188)=0.777018D0;Barr(188)=0.402049D0
Rarr(189)=0.344074D0;Garr(189)=0.780029D0;Barr(189)=0.397381D0
Rarr(190)=0.352360D0;Garr(190)=0.783011D0;Barr(190)=0.392636D0
Rarr(191)=0.360741D0;Garr(191)=0.785964D0;Barr(191)=0.387814D0
Rarr(192)=0.369214D0;Garr(192)=0.788888D0;Barr(192)=0.382914D0
Rarr(193)=0.377779D0;Garr(193)=0.791781D0;Barr(193)=0.377939D0
Rarr(194)=0.386433D0;Garr(194)=0.794644D0;Barr(194)=0.372886D0
Rarr(195)=0.395174D0;Garr(195)=0.797475D0;Barr(195)=0.367757D0
Rarr(196)=0.404001D0;Garr(196)=0.800275D0;Barr(196)=0.362552D0
Rarr(197)=0.412913D0;Garr(197)=0.803041D0;Barr(197)=0.357269D0
Rarr(198)=0.421908D0;Garr(198)=0.805774D0;Barr(198)=0.351910D0
Rarr(199)=0.430983D0;Garr(199)=0.808473D0;Barr(199)=0.346476D0
Rarr(200)=0.440137D0;Garr(200)=0.811138D0;Barr(200)=0.340967D0
Rarr(201)=0.449368D0;Garr(201)=0.813768D0;Barr(201)=0.335384D0
Rarr(202)=0.458674D0;Garr(202)=0.816363D0;Barr(202)=0.329727D0
Rarr(203)=0.468053D0;Garr(203)=0.818921D0;Barr(203)=0.323998D0
Rarr(204)=0.477504D0;Garr(204)=0.821444D0;Barr(204)=0.318195D0
Rarr(205)=0.487026D0;Garr(205)=0.823929D0;Barr(205)=0.312321D0
Rarr(206)=0.496615D0;Garr(206)=0.826376D0;Barr(206)=0.306377D0
Rarr(207)=0.506271D0;Garr(207)=0.828786D0;Barr(207)=0.300362D0
Rarr(208)=0.515992D0;Garr(208)=0.831158D0;Barr(208)=0.294279D0
Rarr(209)=0.525776D0;Garr(209)=0.833491D0;Barr(209)=0.288127D0
Rarr(210)=0.535621D0;Garr(210)=0.835785D0;Barr(210)=0.281908D0
Rarr(211)=0.545524D0;Garr(211)=0.838039D0;Barr(211)=0.275626D0
Rarr(212)=0.555484D0;Garr(212)=0.840254D0;Barr(212)=0.269281D0
Rarr(213)=0.565498D0;Garr(213)=0.842430D0;Barr(213)=0.262877D0
Rarr(214)=0.575563D0;Garr(214)=0.844566D0;Barr(214)=0.256415D0
Rarr(215)=0.585678D0;Garr(215)=0.846661D0;Barr(215)=0.249897D0
Rarr(216)=0.595839D0;Garr(216)=0.848717D0;Barr(216)=0.243329D0
Rarr(217)=0.606045D0;Garr(217)=0.850733D0;Barr(217)=0.236712D0
Rarr(218)=0.616293D0;Garr(218)=0.852709D0;Barr(218)=0.230052D0
Rarr(219)=0.626579D0;Garr(219)=0.854645D0;Barr(219)=0.223353D0
Rarr(220)=0.636902D0;Garr(220)=0.856542D0;Barr(220)=0.216620D0
Rarr(221)=0.647257D0;Garr(221)=0.858400D0;Barr(221)=0.209861D0
Rarr(222)=0.657642D0;Garr(222)=0.860219D0;Barr(222)=0.203082D0
Rarr(223)=0.668054D0;Garr(223)=0.861999D0;Barr(223)=0.196293D0
Rarr(224)=0.678489D0;Garr(224)=0.863742D0;Barr(224)=0.189503D0
Rarr(225)=0.688944D0;Garr(225)=0.865448D0;Barr(225)=0.182725D0
Rarr(226)=0.699415D0;Garr(226)=0.867117D0;Barr(226)=0.175971D0
Rarr(227)=0.709898D0;Garr(227)=0.868751D0;Barr(227)=0.169257D0
Rarr(228)=0.720391D0;Garr(228)=0.870350D0;Barr(228)=0.162603D0
Rarr(229)=0.730889D0;Garr(229)=0.871916D0;Barr(229)=0.156029D0
Rarr(230)=0.741388D0;Garr(230)=0.873449D0;Barr(230)=0.149561D0
Rarr(231)=0.751884D0;Garr(231)=0.874951D0;Barr(231)=0.143228D0
Rarr(232)=0.762373D0;Garr(232)=0.876424D0;Barr(232)=0.137064D0
Rarr(233)=0.772852D0;Garr(233)=0.877868D0;Barr(233)=0.131109D0
Rarr(234)=0.783315D0;Garr(234)=0.879285D0;Barr(234)=0.125405D0
Rarr(235)=0.793760D0;Garr(235)=0.880678D0;Barr(235)=0.120005D0
Rarr(236)=0.804182D0;Garr(236)=0.882046D0;Barr(236)=0.114965D0
Rarr(237)=0.814576D0;Garr(237)=0.883393D0;Barr(237)=0.110347D0
Rarr(238)=0.824940D0;Garr(238)=0.884720D0;Barr(238)=0.106217D0
Rarr(239)=0.835270D0;Garr(239)=0.886029D0;Barr(239)=0.102646D0
Rarr(240)=0.845561D0;Garr(240)=0.887322D0;Barr(240)=0.099702D0
Rarr(241)=0.855810D0;Garr(241)=0.888601D0;Barr(241)=0.097452D0
Rarr(242)=0.866013D0;Garr(242)=0.889868D0;Barr(242)=0.095953D0
Rarr(243)=0.876168D0;Garr(243)=0.891125D0;Barr(243)=0.095250D0
Rarr(244)=0.886271D0;Garr(244)=0.892374D0;Barr(244)=0.095374D0
Rarr(245)=0.896320D0;Garr(245)=0.893616D0;Barr(245)=0.096335D0
Rarr(246)=0.906311D0;Garr(246)=0.894855D0;Barr(246)=0.098125D0
Rarr(247)=0.916242D0;Garr(247)=0.896091D0;Barr(247)=0.100717D0
Rarr(248)=0.926106D0;Garr(248)=0.897330D0;Barr(248)=0.104071D0
Rarr(249)=0.935904D0;Garr(249)=0.898570D0;Barr(249)=0.108131D0
Rarr(250)=0.945636D0;Garr(250)=0.899815D0;Barr(250)=0.112838D0
Rarr(251)=0.955300D0;Garr(251)=0.901065D0;Barr(251)=0.118128D0
Rarr(252)=0.964894D0;Garr(252)=0.902323D0;Barr(252)=0.123941D0
Rarr(253)=0.974417D0;Garr(253)=0.903590D0;Barr(253)=0.130215D0
Rarr(254)=0.983868D0;Garr(254)=0.904867D0;Barr(254)=0.136897D0
Rarr(255)=0.993248D0;Garr(255)=0.906157D0;Barr(255)=0.143936D0
end subroutine