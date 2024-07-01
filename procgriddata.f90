module procgridmod 
use deftype

type pleequa
real*8 a,b,c,d    !plane pleequa ax+by+cz+d=0
end type

contains !some trivial math functions, mainly invoked by doproject

!--------- distance between two 3D points
real*8 function dist2p(pta,ptb) 
type(content) pta,ptb
dist2p=dsqrt((pta%x-ptb%x)**2+(pta%y-ptb%y)**2+(pta%z-ptb%z)**2)
end function
!--------- Return vector from two content type points
function vecfrom2p(veca,vecb)
type(content) veca,vecb,vecfrom2p
vecfrom2p%x=vecb%x-veca%x
vecfrom2p%y=vecb%y-veca%y
vecfrom2p%z=vecb%z-veca%z
vecfrom2p%value=0D0
end function
!--------- Vector multiplication
function vecmul(a,b)
type(content) a,b,vecmul
vecmul%x=a%y*b%z-a%z*b%y
vecmul%y=-(a%x*b%z-a%z*b%x)
vecmul%z=a%x*b%y-a%y*b%x
end function
!--------- project a point to a given plane
function protoplane_pos(point,plane,fa)
type(content) protoplane_pos,point,fa !fa is normal vector
type(pleequa) plane
real*8 scale !!scale of vector of trans point from a1 to a2
! write(*,*) "1",point
! write(*,*) "2",plane
! write(*,*) "3",fa
scale=(plane%a*point%x+plane%b*point%y+plane%c*point%z+plane%d)/&
      (plane%a*fa%x+plane%b*fa%y+plane%c*fa%z)
protoplane_pos%x=point%x-scale*fa%x
protoplane_pos%y=point%y-scale*fa%y
protoplane_pos%z=point%z-scale*fa%z
protoplane_pos%value=point%value
end function
!--------- Value of determinant of 2x2 matrix
function det2_2(vala,valb,valc,vald)
real*8 vala,valb,valc,vald,det2_2
det2_2=vala*vald-valb*valc
end function
!--------- cos value of two vectors
real*8 function cos_2v(veca,vecb)  
type(content) veca,vecb
cos_2v=abs(veca%x*vecb%x+veca%y*vecb%y+veca%z*vecb%z)/dsqrt(veca%x**2+veca%y**2+veca%z**2)/dsqrt(vecb%x**2+vecb%y**2+vecb%z**2)
end function

end module

!!-------- Process grid data. This module are mainly transplanted from GsGrid
!We use Bohr in interal processes, but convert to Angstrom when input and output
subroutine procgriddata
use deftype
use defvar
use util
use GUI
implicit real*8 (a-h,o-z)
real*8,allocatable :: avgdata(:,:)
integer,allocatable :: atmlist(:),atmlist2(:)
logical,allocatable :: cub_do(:,:,:)
character gridfile2*200,gridfilenew*200,atmidxfile*200,atmidxfile2*200,c200tmp*200,tmpchar,c2000tmp*2000,selectyn
type(content) useratom(3),maxv,minv

if (.not.allocated(cubmat)) then
	write(*,"(a)") " Error: Grid data has not been loaded or generated! If you want to load a grid data file now, &
	&input its path, e.g. C:\nico.cub, else input 0 to return"
	do while(.true.)
		read(*,"(a)") c200tmp
		if (c200tmp(1:1)=='0') then
			return
		else
			inquire(file=c200tmp,exist=alive)
			if (alive) then
				inamelen=len_trim(c200tmp)
				!Only load grid data, do not perturb other variables
				if (c200tmp(inamelen-2:inamelen)=="grd") then
					call readgrd(c200tmp,1,1)
				else if (c200tmp(inamelen-2:inamelen)=="cub".or.c200tmp(inamelen-3:inamelen)=="cube") then
					call readcube(c200tmp,1,1)
				else
					write(*,*) "Error: Unknown file type, input again"
					cycle
				end if
				exit
			else
				write(*,*) "Cannot find the file, input again"
			end if
		end if
	end do
end if

do while(.true.)
	write(*,*)
	write(*,*) "                 ============= Process grid data =============="
	write(*,*) "-2 Visualize isosurface of present grid data"
	write(*,*) "-1 Return to main menu"
	write(*,*) "0 Export present grid data to Gaussian-type cube file (.cub)"
	write(*,*) "1 Output all data points with value and coordinate to output.txt"
	write(*,*) "2 Output data points in a XY plane by specifying Z to output.txt"
	write(*,*) "3 Output data points in a YZ plane by specifying X to output.txt"
	write(*,*) "4 Output data points in a XZ plane by specifying Y to output.txt"
	write(*,*) "5 Output average data of XY planes in a range of Z to output.txt"
	write(*,*) "6 Output average data of YZ planes in a range of X to output.txt"
	write(*,*) "7 Output average data of XZ planes in a range of Y to output.txt"
	write(*,*) "8 Output data points in a plane defined by three atom indices to output.txt"
	write(*,*) "9 Output data points in a plane defined by three points to output.txt"
	write(*,*) "10 Output data points in specified value range to output.txt"
	write(*,*) "11 Grid data calculation"
	write(*,"(a)") " 12 Map values of a cube file to specified isosurface of present grid data"
	write(*,*) "13 Set value of the grid points that far away from / close to some atoms"
	write(*,*) "14 Set value of the grid points outside overlap region of two fragments"
	write(*,*) "15 If data value is within certain range, set it to a specified value"
	write(*,*) "16 Scale data range of present grid data"
	write(*,*) "17 Show statistic data of grid points in specific spatial and value ranges"
	write(*,*) "18 Plot (local) integral curve or plane-averaged in X/Y/Z direction"
	read(*,*) isel
    
    if (isel>=2.and.isel<=7.and.ifgridortho()/=1) then
        write(*,*) "Error: This function is only available for orthogonal grid!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        cycle
    end if
	
	if (isel==-2) then
		call drawisosurgui(1)
	else if (isel==-1) then
		exit
	else if (isel==0) then
        call outcube_wrapper
	else if (isel==1) then
		open(10,file="output.txt",status="replace")
		write(*,*) "Outputting data, please wait..."
		ii=0
		do k=1,nz
			do j=1,ny
				do i=1,nx
					call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					write(10,"(3f10.5,f22.15)") tmpx*b2a,tmpy*b2a,tmpz*b2a,cubmat(i,j,k)
				end do
			end do
			call showprog(k,nz)
		end do
		write(*,*)
		write(*,*) "The data have been saved to output.txt in current folder"
		write(*,"(a)") " The first three columns correspond to X,Y,Z, unit is Angstrom, the last column is data value"
		close(10)
	else if (isel==2) then
		write(*,*) "Input Z (in Angstrom) to define a XY plane, e.g. 0.52"
		read(*,*) posZ
		posZ=posZ/b2a
		rmindist=abs(orgz-posZ)
		k=1
		do icycz=2,nz
			disttmp=abs((orgz+(icycz-1)*dz)-posZ)
			if (disttmp<rmindist) then 
				rmindist=disttmp
				k=icycz
			end if
		end do
		write(*,"(a,f14.6,a)") " The X-Y plane closest to your input is Z=",(orgz+(k-1)*dz)*b2a," Angstrom"
		open(10,file="output.txt",status="replace")
		do i=1,nx
			do j=1,ny
				write(10,"(3f11.6,f22.15)") (orgx+(i-1)*dx)*b2a,(orgy+(j-1)*dy)*b2a,(orgz+(k-1)*dz)*b2a,cubmat(i,j,k)
			end do
		end do
		close(10)
		write(*,*) "The data have been exported to output.txt in current folder"
		write(*,*) "The first three columns correspond to X,Y,Z, unit is Angstrom"
		
	else if (isel==3) then
		write(*,*) "Input X (in Angstrom) to define a YZ plane, e.g. 0.52"
		read(*,*) posX
		posX=posX/b2a
		rmindist=abs(orgx-posX)
		i=1
		do icycx=2,nx
			disttmp=abs((orgx+(icycx-1)*dx)-posX)
			if (disttmp<rmindist) then 
				rmindist=disttmp
				i=icycx
			end if
		end do
		write(*,"(a,f14.6,a)") " The YZ plane closest to your input is X=",(orgx+(i-1)*dx)*b2a," Angstrom"
		open(10,file="output.txt",status="replace")
		do j=1,ny
			do k=1,nz
				write(10,"(3f11.6,f22.15)") (orgx+(i-1)*dx)*b2a,(orgy+(j-1)*dy)*b2a,(orgz+(k-1)*dz)*b2a,cubmat(i,j,k)
			end do
		end do
		close(10)
		write(*,*) "The data have been exported to output.txt in current folder"
		write(*,*) "The first three columns correspond to X,Y,Z, unit is Angstrom"
		
	else if (isel==4) then
		write(*,*) "Input Y (in Angstrom) to define a XZ plane, e.g. 0.52"
		read(*,*) posY
		posY=posY/b2a
		rmindist=abs(orgy-posY)
		j=1
		do icycy=2,ny
			disttmp=abs((orgy+(icycy-1)*dy)-posY)
			if (disttmp<rmindist) then 
				rmindist=disttmp
				j=icycy
			end if
		end do
		write(*,"(a,f14.6,a)") " The XZ plane closest to your input is Y=",(orgy+(j-1)*dy)*b2a," Angstrom"
		open(10,file="output.txt",status="replace")
		do i=1,nx
			do k=1,nz
				write(10,"(3f11.6,f22.15)") (orgx+(i-1)*dx)*b2a,(orgy+(j-1)*dy)*b2a,(orgz+(k-1)*dz)*b2a,cubmat(i,j,k)
			end do
		end do
		close(10)
		write(*,*) "The data have been exported to output.txt in current folder"
		write(*,*) "The column 1/2/3/4 correspond to X,Y,Z,value respectively unit is Angstrom"
		
	else if (isel==5) then
		if (allocated(avgdata)) deallocate(avgdata)
		allocate(avgdata(nx,ny))
		avgdata=0D0
		write(*,*) "Input the range of Z (Angstrom) for XY planes, e.g. 3.0 12.4"
		read(*,*) rangelow,rangehigh
		rangelow=rangelow/b2a
		rangehigh=rangehigh/b2a
		nlayers=0
		do icyck=1,nz
			if ((orgz+(icyck-1)*dz)>rangelow.and.(orgz+(icyck-1)*dz)<rangehigh) then
				avgdata=avgdata+cubmat(:,:,icyck)
				nlayers=nlayers+1
			end if
		end do
		avgdata=avgdata/nlayers
		write(*,"(' There are ',i8,' layers within the range')") nlayers
		open(10,file="output.txt",status="replace")
		do i=1,nx
			do j=1,ny
				write(10,"(2f11.6,f22.15)") (orgx+(i-1)*dx)*b2a,(orgy+(j-1)*dy)*b2a,avgdata(i,j)
			end do
		end do
		close(10)
		write(*,*) "The data have been exported to output.txt in current folder"
		write(*,*) "Column 1,2,3 correspond to X,Y,value respectively"
		
	else if (isel==6) then
		if (allocated(avgdata)) deallocate(avgdata)
		allocate(avgdata(ny,nz))
		avgdata=0D0
		write(*,*) "Input the range of X (Angstrom) for YZ planes, e.g. 3.0 12.4"
		read(*,*) rangelow,rangehigh
		rangelow=rangelow/b2a
		rangehigh=rangehigh/b2a
		nlayers=0
		do icyci=1,nx
			if ((orgx+(icyci-1)*dx)>rangelow.and.(orgx+(icyci-1)*dx)<rangehigh) then
				avgdata=avgdata+cubmat(icyci,:,:)
				nlayers=nlayers+1
			end if
		end do
		avgdata=avgdata/nlayers
		write(*,"(' There are ',i8,' layers within the range')") nlayers
		open(10,file="output.txt",status="replace")
		do j=1,ny
			do k=1,nz
				write(10,"(2f11.6,f22.15)") (orgy+(j-1)*dy)*b2a,(orgz+(k-1)*dz)*b2a,avgdata(j,k)
			end do
		end do
		close(10)
		write(*,*) "The data have been exported to output.txt in current folder"
		write(*,*) "Column 1,2,3 correspond to Y,Z,value respectively"
				
	else if (isel==7) then
		if (allocated(avgdata)) deallocate(avgdata)
		allocate(avgdata(nx,nz))
		avgdata=0D0
		write(*,*) "Input the range of Y (Angstrom) for XZ planes, e.g. 3.0 12.4"
		read(*,*) rangelow,rangehigh
		rangelow=rangelow/b2a
		rangehigh=rangehigh/b2a
		nlayers=0
		do icycj=1,ny
			if ((orgy+(icycj-1)*dy)>=rangelow.and.(orgy+(icycj-1)*dy)<=rangehigh) then
				avgdata=avgdata+cubmat(:,icycj,:)
				nlayers=nlayers+1
			end if
		end do
		avgdata=avgdata/nlayers
		write(*,"(' There are ',i8,' layers within the range')") nlayers
		open(10,file="output.txt",status="replace")
		do i=1,nx
			do k=1,nz
				write(10,"(2f11.6,f22.15)") (orgx+(i-1)*dx)*b2a,(orgz+(k-1)*dz)*b2a,avgdata(i,k)
			end do
		end do
		close(10)
		write(*,*) "The data have been exported to output.txt in current folder"
		write(*,*) "Column 1,2,3 correspond to X,Z,value respectively"
		
	else if (isel==8) then
		do while(.true.)
			write(*,*) "Please input indices of three atoms to define a plane, e.g. 4,8,3"
			read(*,*) iatm1,iatm2,iatm3
			if (iatm1<=ncenter.and.iatm1>=1.and.iatm2<=ncenter.and.iatm2>=1.and.iatm3<=ncenter.and.iatm3>=1) exit
			write(*,*) "The atom indices are invalid, input again"
		end do
		useratom(1)%x=a(iatm1)%x
		useratom(1)%y=a(iatm1)%y
		useratom(1)%z=a(iatm1)%z
		useratom(2)%x=a(iatm2)%x
		useratom(2)%y=a(iatm2)%y
		useratom(2)%z=a(iatm2)%z
		useratom(3)%x=a(iatm3)%x
		useratom(3)%y=a(iatm3)%y
		useratom(3)%z=a(iatm3)%z
		write(*,"(' The coordinate (Angstrom) of the three atoms you inputted:')")
		write(*,"(' Atom 1  x, y, z:',3f12.6)") useratom(1)%x*b2a,useratom(1)%y*b2a,useratom(1)%z*b2a
		write(*,"(' Atom 2  x, y, z:',3f12.6)") useratom(2)%x*b2a,useratom(2)%y*b2a,useratom(2)%z*b2a
		write(*,"(' Atom 3  x, y, z:',3f12.6)") useratom(3)%x*b2a,useratom(3)%y*b2a,useratom(3)%z*b2a
		call doproject(useratom)
	else if (isel==9) then
	
		write(*,*) "Please input three points to define a plane (in Angstrom)"
		do ipt=1,3
			write(*,"('Now input x,y,z of point ',i5,'  e.g. 3.0,1.2,4.5')") ipt
			read(*,*) useratom(ipt)%x,useratom(ipt)%y,useratom(ipt)%z
		end do
		useratom%x=useratom%x/b2a !Convert to Bohr
		useratom%y=useratom%y/b2a
		useratom%z=useratom%z/b2a
		call doproject(useratom)
	else if (isel==10) then
		write(*,*) "Input the value range of data points you want to output, e.g. 0.5,0.9"
		write(*,"(a)") " Note: If the two values are identical, points with value within deviation of 3% will be outputted"
		read(*,*) rangelow,rangehigh
		if (rangelow==rangehigh) then
			rangelow=rangelow-0.03D0*abs(rangelow)
			rangehigh=rangehigh+0.03D0*abs(rangehigh)
		end if
		write(*,"(' Points with value within ',1PE13.5, ' and ',1PE13.5,' will be outputted')") rangelow,rangehigh
		write(*,*) "Outputting, please wait..."

		open(10,file="output.txt",status="replace")
		do k=1,nz
			do j=1,ny
				do i=1,nx
					if (cubmat(i,j,k)>=rangelow.and.cubmat(i,j,k)<=rangehigh) then
                        call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					    write(10,"(3f11.6,f22.15)") tmpx*b2a,tmpy*b2a,tmpz*b2a,cubmat(i,j,k)
                    end if
				end do
			end do
		end do
		close(10)
		write(*,*) "The data have been exported to output.txt in current folder"
		write(*,*) "The column 1/2/3/4 correspond to X,Y,Z,value respectively, unit is Angstrom"
	
	else if (isel==11) then
		write(*,*) "               --------------- Grid data calculation ------------"
		write(*,*) "0 Return"
		write(*,*) "1 Add a constant                       e.g. A+0.1=C"
		write(*,*) "2 Add a grid file                      i.e. A+B=C"
		write(*,*) "3 Subtract a constant                  e.g. A-0.1=C"
		write(*,*) "4 Subtract a grid file                 i.e. A-B=C"
		write(*,*) "5 Multiplied by a constant             e.g. A*0.1=C"
		write(*,*) "6 Multiplied by a grid file            i.e. A*B=C"
		write(*,*) "7 Divided by a constant                e.g. A/5.2=C"
		write(*,*) "8 Divided by a grid file               i.e. A/B=C"
		write(*,*) "9 Exponentiation                       e.g. A^1.3=C"
		write(*,*) "10 Square sum with a grid file         i.e. A^2+B^2=C"
		write(*,*) "11 Square subtract with a grid file    i.e. A^2-B^2=C"
		write(*,*) "12 Get average with a grid file        i.e. (A+B)/2=C"
		write(*,*) "13 Get absolute value                  i.e. |A|=C"
		write(*,*) "14 Get exponential value with base 10  i.e. 10^A=C"
		write(*,*) "15 Get logarithm with base 10          i.e. log10(A)=C"
		write(*,*) "16 Get natural exponential value       i.e. e^A=C"
		write(*,*) "17 Get natural logarithm               i.e. ln(A)=C"
		write(*,*) "18 Add a grid file multiplied by a value  i.e. A+0.4*B=C"
		write(*,*) "19 The same as 6 but with weighting function min(|A|,|B|)/max(|A|,|B|)"
		write(*,*) "20 Multiplied by a coordinate variable"
		write(*,*) "21 Get minimal value with another function    i.e. min(A,B)"
		write(*,*) "22 Get min(|A|,|B|)"
		read(*,*) isel2
		
		if (isel2==0) then
			continue
		else if (isel2==1.or.isel2==3.or.isel2==5.or.isel2==7.or.isel2==9) then
			write(*,*) "Input the value for the calculation, e.g. 2.3"
			read(*,*) calconstant
			if (isel2==1) cubmat=cubmat+calconstant
			if (isel2==3) cubmat=cubmat-calconstant
			if (isel2==5) cubmat=cubmat*calconstant
			if (isel2==7) cubmat=cubmat/calconstant
			if (isel2==9) cubmat=cubmat**calconstant
		else if (isel2==2.or.isel2==4.or.isel2==6.or.isel2==8.or.isel2==10.or.isel2==11.or.isel2==12.or.isel2==18.or.isel2==19.or.isel2==21.or.isel2==22) then
			do while(.true.)
				write(*,*) "Input another file containing grid data (.cub, .grd or CHGCAR/CHG)"
                write(*,*) "e.g. C:\mius\Umi.cub"
				read(*,"(a)") gridfile2
				inquire(file=gridfile2,exist=alive)
				if (alive) exit
				write(*,*) "File not found, input again"
				write(*,*)
			end do
			inamelen=len_trim(gridfile2)
			if (gridfile2(inamelen-2:inamelen)=="cub".or.gridfile2(inamelen-3:inamelen)=="cube") then
				call readcubetmp(gridfile2,1,inconsis)
			else if (gridfile2(inamelen-2:inamelen)=="grd") then
				call readgrdtmp(gridfile2,inconsis)
			else if (index(gridfile2,"CHG")/=0) then
				call readVASPgrdtmp(gridfile2,inconsis)
			end if
			if (inconsis==1) cycle
			if (isel2==2) cubmat=cubmat+cubmattmp
			if (isel2==4) cubmat=cubmat-cubmattmp
			if (isel2==6) cubmat=cubmat*cubmattmp
			if (isel2==8) cubmat=cubmat/cubmattmp
			if (isel2==10) cubmat=cubmat**2+cubmattmp**2
			if (isel2==11) cubmat=cubmat**2-cubmattmp**2
			if (isel2==12) cubmat=(cubmat+cubmattmp)/2D0
			if (isel2==19) cubmat=cubmat*cubmattmp*min(abs(cubmat),abs(cubmattmp)) / max(abs(cubmat),abs(cubmattmp))
			if (isel2==21) cubmat=min(cubmat,cubmattmp)
			if (isel2==22) cubmat=min(abs(cubmat),abs(cubmattmp))
			
			if (isel2==18) then
				write(*,*) "Input the value to be multiplied to the just loaded cube file"
				read(*,*) tmpval
				cubmat=cubmat+tmpval*cubmattmp
			end if
			deallocate(cubmattmp)
		else if (isel2==13) then
			cubmat=abs(cubmat)
		else if (isel2==14) then
			cubmat=10**cubmat
		else if (isel2==15) then
			cubmat=log10(cubmat)
		else if (isel2==16) then
			cubmat=exp(cubmat)
		else if (isel2==17) then
			cubmat=log(cubmat)
		else if (isel2==20) then
			write(*,*) "Multiplied by which variable? Input ""X"" or ""Y"" or ""Z"""
			read(*,*) tmpchar
		    do k=1,nz
			    do j=1,ny
				    do i=1,nx
                        call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					    if (tmpchar=='x'.or.tmpchar=='X') then
                            cubmat(i,j,k)=cubmat(i,j,k)*tmpx
                        else if (tmpchar=='y'.or.tmpchar=='Y') then
                            cubmat(i,j,k)=cubmat(i,j,k)*tmpy
                        else if (tmpchar=='z'.or.tmpchar=='Z') then
                            cubmat(i,j,k)=cubmat(i,j,k)*tmpz
                        end if
				    end do
			    end do
		    end do
		end if
		if (isel2/=0) write(*,*) "Done, the grid data has been updated"
	
	else if (isel==12) then
		write(*,*) "Input a value to define a isosurface of present grid data, e.g. 0.001"
		read(*,*) value_ref
		write(*,*) "Input allowed deviation (%)   e.g. 4"
		read(*,*) deviation
		rangehigh=value_ref+abs(value_ref)*0.01D0*deviation
		rangelow=value_ref-abs(value_ref)*0.01D0*deviation
		write(*,"(' Value between ',1PE13.5, ' and ',1PE13.5,' are regarded as isosurface points')") rangelow,rangehigh

		do while(.true.)
			write(*,*) "Map which grid file to the isosurface?  e.g. C:\t.cub"
			read(*,"(a)") gridfile2
			inquire(file=gridfile2,exist=alive)
			if (alive) exit
			write(*,*) "File not found, input again"
			write(*,*)
		end do
		call readcubetmp(gridfile2,1,inconsis)
		open(10,file="output.txt",status="replace")
        do k=1,nz
			do j=1,ny
		        do i=1,nx
					if (cubmat(i,j,k)>=rangelow.and.cubmat(i,j,k)<=rangehigh) then
                        call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					    write(10,"(3f11.6,f22.15)") tmpx*b2a,tmpy*b2a,tmpz*b2a,cubmattmp(i,j,k)
                    end if
				end do
			end do
		end do
		close(10)
		write(*,"(a)") "The mapped data on the isosurface have been exported to output.txt in current folder"
		write(*,*) "The column 1/2/3/4 correspond to X,Y,Z,value respectively, unit is Angstrom"
		deallocate(cubmattmp)
		
	else if (isel==13) then
		write(*,*) "Input the value for scaling vdW radius, e.g. 1.3"
		write(*,"(a)") " Note: Inputting a positive(negative) value means the grid points outside(inside) the distance cutoff will be set to a given value"
		read(*,*) vdwscale
		write(*,*) "Set to what value? e.g. 100"
		read(*,*) setval
		write(*,*) "Select mode of inputting atomic indices"
		write(*,*) "1 Load atomic indices from a plain text file"
		write(*,*) "2 Inputting atomic indices by hand"
		read(*,*) imodinput
		if (imodinput==1) then
			write(*,*) "Input the filename that recorded atom indices, e.g. C:\niconiconi\atmlist.txt"
			do while(.true.)
				read(*,"(a)") atmidxfile
				inquire(file=atmidxfile,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			open(10,file=atmidxfile,status="old")
			read(10,*) numatmidx
			allocate(atmlist(numatmidx))
			read(10,*) atmlist
			close(10)
		else if (imodinput==2) then
			write(*,*) "Input atomic indices, e.g. 3,5,7-20,25"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,numatmidx)
			allocate(atmlist(numatmidx))
			call str2arr(c2000tmp,numatmidx,atmlist)
		end if

		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,xnow,ynow,znow)
					if (vdwscale>0) then
						do ind=1,numatmidx
							iatm=atmlist(ind)
							dist2=(xnow-a(iatm)%x)**2+(ynow-a(iatm)%y)**2+(znow-a(iatm)%z)**2
							if (dist2< (vdwscale*vdwr(a(iatm)%index))**2 ) exit !If this point is near any atom in the list, it will be unchanged
							if (ind==numatmidx) cubmat(i,j,k)=setval !Final cycle, suggests that this point is outside any atom in the list, set its value
						end do
					else if (vdwscale<0) then
						do ind=1,numatmidx
							iatm=atmlist(ind)
							dist2=(xnow-a(iatm)%x)**2+(ynow-a(iatm)%y)**2+(znow-a(iatm)%z)**2
							if (dist2< (vdwscale*vdwr(a(iatm)%index))**2 ) then
								cubmat(i,j,k)=setval
								exit
							end if
						end do
					end if
				end do
			end do
		end do
		if (vdwscale>0) write(*,"(a,f10.6,a,E16.8)") " Done! The value of all points ouside",abs(vdwscale)," times vdW radius of specified atoms have been set to",setval
		if (vdwscale<0) write(*,"(a,f10.6,a,E16.8)") " Done! The value of all points inside",abs(vdwscale)," times vdW radius of specified atoms have been set to",setval
		deallocate(atmlist)
		
	else if (isel==14) then
		write(*,*) "Input the value for scaling vdW radius  e.g. 1.3"
		read(*,*) vdwscale
		write(*,*) "Set to what value?"
		read(*,*) setval
		write(*,*) "Select mode of inputting atomic indices"
		write(*,*) "1 Load atomic indices from a plain text file"
		write(*,*) "2 Inputting atomic indices by hand"
		read(*,*) imodinput
		if (imodinput==1) then
			write(*,*) "Input path of the file recording atom indices of fragment 1"
			write(*,*) "e.g. C:\nicomaki\frag1.txt"
			do while(.true.)
				read(*,"(a)") atmidxfile
				inquire(file=atmidxfile,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			write(*,*) "Input path of the file recording atom indices of fragment 2"
			write(*,*) "e.g. C:\nicomaki\frag2.txt"
			do while(.true.)
				read(*,"(a)") atmidxfile2
				inquire(file=atmidxfile2,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do

			open(10,file=atmidxfile,status="old")
			read(10,*) numatmidx
			allocate(atmlist(numatmidx))
			read(10,*) atmlist
			close(10)
			open(10,file=atmidxfile2,status="old")
			read(10,*) numatmidx2
			allocate(atmlist2(numatmidx2))
			read(10,*) atmlist2
			close(10)
		else if (imodinput==2) then
			write(*,*) "Input atomic indices for fragment 1, e.g. 3,5,7-20,25"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,numatmidx)
			allocate(atmlist(numatmidx))
			call str2arr(c2000tmp,numatmidx,atmlist)
			write(*,*) "Input atomic indices for fragment 2, e.g. 3,5,7-20,25"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,numatmidx2)
			allocate(atmlist2(numatmidx2))
			call str2arr(c2000tmp,numatmidx2,atmlist2)
		end if

		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,xnow,ynow,znow)
					iwithoutfrag1=0
					do ind=1,numatmidx
						iatm=atmlist(ind)
						dist2=(xnow-a(iatm)%x)**2+(ynow-a(iatm)%y)**2+(znow-a(iatm)%z)**2
						if (dist2< (vdwscale*vdwr(a(iatm)%index))**2 ) exit
						if (ind==numatmidx) iwithoutfrag1=1
					end do
					if (iwithoutfrag1==1) then !This point is without fragment 1, hence mustn't be within overlap of frag.1 and 2, so set its value now
						cubmat(i,j,k)=setval
					else !This point is within fragment 1, check if it is also within fragment 2, if yes then lzeave it unchanged, else set its value
						do ind=1,numatmidx2
							iatm=atmlist2(ind)
							dist2=(xnow-a(iatm)%x)**2+(ynow-a(iatm)%y)**2+(znow-a(iatm)%z)**2
							if (dist2< (vdwscale*vdwr(a(iatm)%index))**2 ) exit
							if (ind==numatmidx2) cubmat(i,j,k)=setval !Final cycle, suggest this point is simutaneously without fragment 1 and 2
						end do				
					end if
				end do
			end do
		end do
		write(*,"(a,1PE16.8)") " Done! The value of all points outside overlap region of fragment 1 and 2 have been set to",setval
		deallocate(atmlist,atmlist2)
		
	else if (isel==15) then
		write(*,*) "Input lower limit and upper limit  e.g. 3.9,4.2"
		read(*,*) flow,fhigh
		write(*,"(' If the value >= ',E15.8,' and <=',E15.8,', set it to what value?')") flow,fhigh
		read(*,*) calconstant
		where (cubmat>=flow.and.cubmat<=fhigh) cubmat=calconstant
		
	else if (isel==16) then !Scale certain range of original data to specified range
		write(*,*) "Input lower and upper limit of original data range, e.g. 0,1.5"
		write(*,"(a)") " Note: If input 0,0 then minimum and maximum value of present grid data will be used"
		read(*,*) orglow,orgup
		if (orglow==0D0.and.orgup==0D0) then
			orglow=minval(cubmat)
			orgup=maxval(cubmat)
		end if
		write(*,*) "Input lower and upper limit of new data range, e.g. 0,255"
		read(*,*) rnewlow,rnewup
		where (cubmat>orgup) cubmat=orgup
		where (cubmat<orglow) cubmat=orglow
		cubmat=cubmat-orglow
		ratiofac=(rnewup-rnewlow)/(orgup-orglow)
		cubmat=cubmat*ratiofac
		cubmat=cubmat+rnewlow
		write(*,*) "Done!"
		
	else if (isel==17) then
		call gridmaxxyz(rhighx,rhighy,rhighz)
		call gridminxyz(rlowx,rlowy,rlowz)
		rlowv=minval(cubmat)
		rhighv=maxval(cubmat)
		write(*,*) "1 Obtain statistic data for all grid points"
		write(*,*) "2 Obtain statistic data for grid points in specific spatial and value ranges"
		read(*,*) iselrange
        if (iselrange==1) then
			iregion=0
		else if (iselrange==2) then
			write(*,"(a)") " Input the lower and upper limits of value, e.g. -2,3.25"
            write(*,*) "If you do not want to set this constraint, press ENTER button directly"
			read(*,"(a)") c200tmp
			if (index(c200tmp,'a')==0.and.c200tmp/=" ") read(c200tmp,*) rlowv,rhighv
            write(*,*) "How to define the spatial region?"
            write(*,*) "1 Rectangular region"
            write(*,*) "2 Cylindrical region"
            write(*,*) "3 Spherical region"
            read(*,*) iregion
            if (iregion==1) then
				write(*,"(a)") " Input the lower and upper limits of X coordinate (in Angstrom), e.g. -40,33.5"
				write(*,*) "If you do not want to set constraint, input ""a"""
				read(*,"(a)") c200tmp
				if (index(c200tmp,'a')/=0) then
					read(c200tmp,*) rlowx,rhighx
					rlowx=rlowx/b2a
					rhighx=rhighx/b2a
				end if
				write(*,"(a)") " Input the lower and upper limits of Y coordinate (in Angstrom)"
				write(*,*) "If you do not want to set this constraint, input ""a"""
				read(*,"(a)") c200tmp
				if (index(c200tmp,'a')/=0) then
					read(c200tmp,*) rlowy,rhighy
					rlowy=rlowy/b2a
					rhighy=rhighy/b2a
				end if
				write(*,"(a)") " Input the lower and upper limits of Z coordinate (in Angstrom), e.g. -40,33.5"
				write(*,*) "If you do not want to set this constraint, input ""a"""
				read(*,"(a)") c200tmp
				if (index(c200tmp,'a')/=0) then
					read(c200tmp,*) rlowz,rhighz
					rlowz=rlowz/b2a
					rhighz=rhighz/b2a
				end if
            else if (iregion==2) then
				write(*,*) "Input X,Y,Z of the first end point of cylinder (in Angstrom), e.g. 0,1,-1.5"
                read(*,*) cylx1,cyly1,cylz1
				write(*,*) "Input X,Y,Z of the second end point of cylinder (in Angstrom), e.g. 0,1,1.5"
                read(*,*) cylx2,cyly2,cylz2
                write(*,*) "Input radius of the cylinder (in Angstrom), e.g. 2.1"
                read(*,*) cylrad
                cylx1=cylx1/b2a
                cyly1=cyly1/b2a
                cylz1=cylz1/b2a
                cylx2=cylx2/b2a
                cyly2=cyly2/b2a
                cylz2=cylz2/b2a
                cylrad=cylrad/b2a
            else if (iregion==3) then
				write(*,*) "Input X,Y,Z coordinate of sphere center (in Angstrom), e.g. 1.0,5.4,-0.1"
                read(*,*) sphcenx,sphceny,sphcenz
                write(*,*) "Input radius of the sphere (in Angstrom), e.g. 3.5"
                read(*,*) sphrad
                sphcenx=sphcenx/b2a
                sphceny=sphceny/b2a
                sphcenz=sphcenz/b2a
                sphrad=sphrad/b2a
            end if
		end if
		
        if (iregion==0) then
			write(*,*) "The geometry range for statistics"
			write(*,"(' Lower and upper limit of X:',2f14.8,' Bohr')") rlowx,rhighx
			write(*,"(' Lower and upper limit of Y:',2f14.8,' Bohr')") rlowy,rhighy
			write(*,"(' Lower and upper limit of Z:',2f14.8,' Bohr')") rlowz,rhighz
		end if
        
		maxv%value=cubmat(1,1,1)
		maxv%x=orgx
		maxv%y=orgy
		maxv%z=orgz
		minv%value=cubmat(1,1,1)
		minv%x=orgx
		minv%y=orgy
		minv%z=orgz
		sumuppos=0D0
		sumupneg=0D0
		cenxpos=0D0
		cenypos=0D0
		cenzpos=0D0
		cenxneg=0D0
		cenyneg=0D0
		cenzneg=0D0
		igoodpointpos=0
		igoodpointneg=0
		sumupsqrtot=0
        allocate(cub_do(nx,ny,nz))
        cub_do=.false.
		do k=1,nz
			do j=1,ny
				do i=1,nx
					valtmp=cubmat(i,j,k)
					if (valtmp<rlowv.or.valtmp>rhighv) cycle
                    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                    if (iregion==1) then !Rectangle
						if (tmpx<rlowx.or.tmpx>rhighx) cycle
						if (tmpy<rlowy.or.tmpy>rhighy) cycle
						if (tmpz<rlowz.or.tmpz>rhighz) cycle
                    else if (iregion==2) then !Cylinder
						rdist=potlinedis(tmpx,tmpy,tmpz,cylx1,cyly1,cylz1,cylx2,cyly2,cylz2)
                        if (rdist>cylrad) cycle
                        cenx=(cylx1+cylx2)/2
                        ceny=(cyly1+cyly2)/2
                        cenz=(cylz1+cylz2)/2
                        dirx=cylx1-cylx2
                        diry=cyly1-cyly2
                        dirz=cylz1-cylz2
                        !Plane equation:  dirx*(x-cenx)+diry*(y-ceny)+dirz*(z-cenz)=0
                        !namely, dirx*x + diry*y + dirz*z -dirx*cenx -diry*ceny -dirz*cenz = 0
                        pleD=-dirx*cenx-diry*ceny-dirz*cenz 
                        call pointABCDdis(tmpx,tmpy,tmpz,dirx,diry,dirz,pleD,dist,0)
                        cyllength=dsqrt(dirx**2+diry**2+dirz**2)
                        if (dist>cyllength/2) cycle
                    else if (iregion==3) then !Sphere
						rdist=dsqrt((tmpx-sphcenx)**2+(tmpy-sphceny)**2+(tmpz-sphcenz)**2)
                        if (rdist>sphrad) cycle
                    end if
                    cub_do(i,j,k)=.true.
					if (valtmp>0) then
						sumuppos=sumuppos+valtmp
						cenxpos=cenxpos+tmpx*valtmp
						cenypos=cenypos+tmpy*valtmp
						cenzpos=cenzpos+tmpz*valtmp
						igoodpointpos=igoodpointpos+1
					else if (valtmp<0) then
						sumupneg=sumupneg+valtmp
						cenxneg=cenxneg+tmpx*valtmp
						cenyneg=cenyneg+tmpy*valtmp
						cenzneg=cenzneg+tmpz*valtmp
						igoodpointneg=igoodpointneg+1
					end if
					if (valtmp>maxv%value) then
						maxv%value=cubmat(i,j,k)
						call getgridxyz(i,j,k,maxv%x,maxv%y,maxv%z)
					end if
					if (valtmp<minv%value) then
						minv%value=cubmat(i,j,k)
						call getgridxyz(i,j,k,minv%x,minv%y,minv%z)
					end if
					sumupsqrtot=sumupsqrtot+valtmp**2
				end do
			end do
		end do
		sumuptot=sumuppos+sumupneg
		cenxtot=(cenxpos+cenxneg)/sumuptot
		cenytot=(cenypos+cenyneg)/sumuptot
		cenztot=(cenzpos+cenzneg)/sumuptot
		cenxpos=cenxpos/sumuppos
		cenypos=cenypos/sumuppos
		cenzpos=cenzpos/sumuppos
		cenxneg=cenxneg/sumupneg
		cenyneg=cenyneg/sumupneg
		cenzneg=cenzneg/sumupneg
		numpt=nx*ny*nz
        call calc_dvol(dvol)
		avgtot=sumuptot/numpt
		stddev=0
        
		!Calculate standard deviation
		do k=1,nz
			do j=1,ny
				do i=1,nx
                    if (cub_do(i,j,k)) then
						stddev=stddev+(cubmat(i,j,k)-avgtot)**2
                    end if
				end do
			end do
		end do
		stddev=dsqrt(stddev/numpt)
        
        write(*,*)
		write(*,"(' The minimum value:',E16.8,' at',3f12.6,' Bohr')") minv%value,minv%x,minv%y,minv%z
		write(*,"(' The maximum value:',E16.8,' at',3f12.6,' Bohr')") maxv%value,maxv%x,maxv%y,maxv%z
		write(*,"(' Differential element:',f15.10,' Bohr^3')") dvol
		write(*,"(' Average value:',E16.8)") avgtot
		write(*,"(' Root mean square (RMS):',E16.8)") dsqrt(sumupsqrtot/numpt)
		write(*,"(' Standard deviation:',E16.8)") stddev
		write(*,*)
		write(*,"(' Volume of positive value space:',f30.10,' Bohr^3')") igoodpointpos*dvol
		write(*,"(' Volume of negative value space:',f30.10,' Bohr^3')") igoodpointneg*dvol
		write(*,"(' Volume of all space:           ',f30.10,' Bohr^3')") (igoodpointpos+igoodpointneg)*dvol
		write(*,*)
		write(*,"(' Summing up positive values:',f30.10)") sumuppos
		write(*,"(' Summing up negative values:',f30.10)") sumupneg
		write(*,"(' Summing up all values:     ',f30.10)") sumuptot
		write(*,*)
		write(*,"(' Integral of positive data:',f30.10)") sumuppos*dvol
		write(*,"(' Integral of negative data:',f30.10)") sumupneg*dvol
		write(*,"(' Integral of all data:     ',f30.10)") sumuptot*dvol
		write(*,*)
		write(*,"(' X,Y,Z of barycenter (in Bohr)')")
		write(*,"(' Positive part:',3f20.8)") cenxpos,cenypos,cenzpos
		write(*,"(' Negative part:',3f20.8)") cenxneg,cenyneg,cenzneg
        
		!If positive and negative cancel each other exactly, namely sumuptot is about zero, total barycenter will be infinitely large
		if (abs(sumuptot)>0.001D0) write(*,"(' Total:        ',3f20.8)") cenxtot,cenytot,cenztot
        
        if (iselrange==2) then
			write(*,*)
			write(*,"(a)") " Do you want to export the grids actually involved in statistics to grid.xyz in current folder for visual check? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y'.or.selectyn=='Y') then
				open(10,file="grid.xyz",status="replace")
                write(10,*) count(cub_do.eqv..true.)
                write(10,*)
				do k=1,nz
					do j=1,ny
						do i=1,nx
							if (cub_do(i,j,k)) then
								call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
								write(10,"(a4,3f12.6)") "X",tmpx*b2a,tmpy*b2a,tmpz*b2a
							end if
						end do
					end do
				end do
                close(10)
                write(*,"(a)") " Done! grid.xyz has been exported in current folder. The coordinates are in Angstrom. You may use VMD program to visualze it"
            end if
        end if
		deallocate(cub_do)
        
	else if (isel==18) then !Integral curve
        call drawintcurve
	end if
end do
end subroutine



!!----------- Project grid data to a specified plane, selectly project to XY plane, then output to plain text file
subroutine doproject(useratom)
use defvar
use procgridmod
implicit real*8 (a-h,o-z)
type(content) vec1,vec2,vec3,useratom(3),pttemp,point1,point2,cubmatpt
type(content),allocatable :: planedata(:)
type(pleequa) userplane
integer flipmethod
vec1=vecfrom2p(useratom(1),useratom(2))
vec2=vecfrom2p(useratom(1),useratom(3))
vec3=vecmul(vec1,vec2)  !normal vector of userplane!!!
userplane%a=vec3%x  !ax+by+cz+d=0
userplane%b=vec3%y
userplane%c=vec3%z
userplane%d=-vec3%x*useratom(1)%x-vec3%y*useratom(1)%y-vec3%z*useratom(1)%z
cosang=cos_2v(content(0,0,1,0),vec3) !cos between userplane and XY plane
ang=2*180*acos(cosang)/(2D0*pi)
write(*,"(' The angle between your defined plane and XY plane is ',f6.2,' degrees')") ang

if (userplane%a==0.and.userplane%b==0) then
	write(*,"(' Warning: You defined plane is parallel to XY plane, use function 2 instead! Now exit...')")
	return
else if (userplane%b==0.and.userplane%c==0) then
	write(*,"(' Warning: You defined plane is parallel to YZ plane, use function 3 instead! Now exit...')")
	return
else if (userplane%a==0.and.userplane%c==0) then
	write(*,"(' Warning: You defined plane is parallel to XZ plane, use function 4 instead! Now exit...')")
	return
end if

distoler=dsqrt(dx**2+dy**2+dz**2)/4D0  !Default tolerance distance between grid to userplane
write(*,*)
write(*,"(' Input tolerance distance for projecting grids to your defined plane, unit is Angstrom. (Note: You can input 0 to use recommended value: ',f9.4,' Angstrom)')") distoler*b2a
read(*,*) temp
temp=temp/b2a
if (temp/=0) distoler=temp

allocate(planedata(nx*ny*nz))
npt=0
do k=1,nz
	do j=1,ny
		do i=1,nx
			cubmatpt%value=cubmat(i,j,k)
            call getgridxyz(i,j,k,cubmatpt%x,cubmatpt%y,cubmatpt%z)
			pttemp=protoplane_pos(cubmatpt,userplane,vec3)  !Project original data a(i,j,k) to userplane
			if (dist2p(pttemp,cubmatpt)<distoler) then
				npt=npt+1
				planedata(npt)=pttemp
			end if
		end do
	end do
end do
write(*,"(i10,' points are presented in the plane you defined')") npt

write(*,*) "If projecting data to XY plane? 1=yes 2=no"
read(*,*) itest

if (itest==1) then
	point1%x=1
	point1%y=(-userplane%d-userplane%a*point1%x)/userplane%b  !find arbitrary two points(x=1 and x=2) in userplane and XY plane(Z=0) crossing
	point1%z=0
	point2%x=2
	point2%y=(-userplane%d-userplane%a*point2%x)/userplane%b
	point2%z=0
	vec1=vecfrom2p(point1,point2)  !crossing vector of XY plane and userplane
	detvalue=det2_2(vec1%x, vec1%y, userplane%a, userplane%b)
	scalex=abs(vec3%x)/dsqrt(vec3%x**2+vec3%y**2)  !The ratio of x component to norm of vec3. scalex/y have direction, namely positive or negative sign
	scaley=abs(vec3%y)/dsqrt(vec3%x**2+vec3%y**2)
	
	if (ang>80) then
		flipmethod=1
	else
		flipmethod=0
	end if

	do i=1,npt
!Use Cramer law to find a point (point1) in crossing vector (namely vec1) of userplane and XY plane
!point1 fulfill two condition (of course, its z=0, since it is laying on XY plane):
!(1) The vector of point1 to planedata(n) is vertical to crossing vector (2)the point is in userplane. so we can use cramer law to solve linear equation to get point1
		point1%x=det2_2(vec1%x*planedata(i)%x+vec1%y*planedata(i)%y, vec1%y, -userplane%d, userplane%b)/detvalue
		point1%y=det2_2(vec1%x, vec1%x*planedata(i)%x+vec1%y*planedata(i)%y, userplane%a, -userplane%d)/detvalue
		temp2=dist2p(planedata(i),point1) - dist2p(planedata(i),point1)*cosang !temp2 is distance between point1 and original grid point i

		if (flipmethod==0) then
			if (planedata(i)%x>point1%x) then       ! used to determine adding direction
				planedata(i)%x=planedata(i)%x+scalex*temp2
			else
				planedata(i)%x=planedata(i)%x-scalex*temp2
			end if
			if (planedata(i)%y>point1%y) then
				planedata(i)%y=planedata(i)%y+scaley*temp2
			else
				planedata(i)%y=planedata(i)%y-scaley*temp2
			end if
		else !Angle near 90 degree
			if (userplane%a*userplane%b>0) then
				if (planedata(i)%z>0) then
					planedata(i)%x=point1%x+dist2p(planedata(i),point1)*scalex
					planedata(i)%y=point1%y+dist2p(planedata(i),point1)*scaley
				else
					planedata(i)%x=point1%x-dist2p(planedata(i),point1)*scalex
					planedata(i)%y=point1%y-dist2p(planedata(i),point1)*scaley
				end if
			else
				if (planedata(i)%z>0) then
					planedata(i)%x=point1%x-dist2p(planedata(i),point1)*scalex
					planedata(i)%y=point1%y+dist2p(planedata(i),point1)*scaley
				else
					planedata(i)%x=point1%x+dist2p(planedata(i),point1)*scalex
					planedata(i)%y=point1%y-dist2p(planedata(i),point1)*scaley
				end if
			end if
		end if
		planedata(i)%z=0
	end do
end if

open(10,file="output.txt",status="replace")
do i=1,npt
	write(10,"(3f11.6,f22.15)") planedata(i)%x*b2a,planedata(i)%y*b2a,planedata(i)%z*b2a,planedata(i)%value
end do
close(10)
write(*,*) "Completed!"
write(*,"(' In the output.txt, column 1/2/3/4 correspond to x/y/z/value respectively, unit is Angstrom')")
deallocate(planedata)
end subroutine


!!------- Integrate grid data along a direction to obtain integral curve
subroutine drawintcurve
use defvar
use util
use GUI
implicit real*8 (a-h,o-z)
character c200tmp*200,c80tmp*80
real*8,allocatable :: intcurve(:),locintcurve(:),pleavgcurve(:),curvepos(:)

if (gridv1(2)==0.and.gridv1(3)==0) then !Can integrate along X axis
    continue
else if (gridv2(1)==0.and.gridv2(3)==0) then !Can integrate along Y axis
    continue
else if (gridv3(1)==0.and.gridv3(2)==0) then !Can integrate along Z axis
    continue
else
    write(*,"(a)") "Error: This function cannot be used for present grid. The grid must satisfy at least one below conditions:"
    write(*,*) "1 Grid vector 1 is parallel along X axis and the other two are in YZ plane"
    write(*,*) "2 Grid vector 2 is parallel along Y axis and the other two are in XZ plane"
    write(*,*) "3 Grid vector 3 is parallel along Z axis and the other two are in XY plane"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
    return
end if

write(*,*) "Integrating in which direction? Input ""X"" or ""Y"" or ""Z"""
read(*,*) c200tmp
if (index(c200tmp,'X')/=0.or.index(c200tmp,'x')/=0) idir=1
if (index(c200tmp,'Y')/=0.or.index(c200tmp,'y')/=0) idir=2
if (index(c200tmp,'Z')/=0.or.index(c200tmp,'z')/=0) idir=3
if (idir==1) then
    if (gridv1(2)==0.and.gridv1(3)==0.and.sum(gridv1*gridv2)==0.and.sum(gridv1*gridv3)==0) then
        continue
    else
        write(*,"(a)") " Error: To integrate along X axis, grid vector 1 must be parallel &
        &to X axis and the other two should be in YZ plane. The present grid does not meet this condition"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
        return
    end if
else if (idir==2) then
    if (gridv2(1)==0.and.gridv2(3)==0.and.sum(gridv2*gridv1)==0.and.sum(gridv2*gridv3)==0) then
        continue
    else
        write(*,"(a)") " Error: To integrate along Y axis, grid vector 2 must be parallel &
        &to Y axis and the other two should be in XZ plane. The present grid does not meet this condition"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
        return
    end if
else if (idir==3) then
    if (gridv3(1)==0.and.gridv3(2)==0.and.sum(gridv3*gridv1)==0.and.sum(gridv3*gridv2)==0) then
        continue
    else
        write(*,"(a)") " Error: To integrate along Z axis, grid vector 3 must be parallel &
        &to Z axis and the other two should be in XY plane. The present grid does not meet this condition"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
        return
    end if
end if

write(*,*) "Input lower and upper limits of the position (in Angstrom), e.g. -30,51.5"
write(*,*) "Press ENTER button directly or input ""a"" can choose the entire range"
read(*,"(a)") c200tmp
if (c200tmp==" ".or.index(c200tmp,"a")/=0) then
	continue
else
	read(c200tmp,*) rlow,rhigh
	rlow=rlow/b2a
	rhigh=rhigh/b2a
end if
write(*,*) "Calculating data..."
tmpintval=0
if (idir==1) then !Direction 1
	allocate(intcurve(nx),locintcurve(nx),curvepos(nx),pleavgcurve(nx))
	ncurpt=nx
    v1len=dsqrt(sum(gridv1**2))
	if (c200tmp==" ".or.index(c200tmp,"a")/=0) then
		rlow=orgx
		rhigh=orgx+(nx-1)*v1len
	end if
	do ix=1,nx
		curvepos(ix)=orgx+(ix-1)*v1len
		if (curvepos(ix)<rlow.or.curvepos(ix)>rhigh) cycle
        call vec2area(gridv2(:),gridv3(:),area)
		locintcurve(ix)=sum(cubmat(ix,:,:))*area
        pleavgcurve(ix)=sum(cubmat(ix,:,:))/(ny*nz)
		tmpintval=tmpintval+locintcurve(ix)*v1len
		intcurve(ix)=tmpintval
	end do
else if (idir==2) then !Direction 2
	allocate(intcurve(ny),locintcurve(ny),curvepos(ny),pleavgcurve(ny))
	ncurpt=ny
    v2len=dsqrt(sum(gridv2**2))
	if (c200tmp==" ".or.index(c200tmp,"a")/=0) then
		rlow=orgy
		rhigh=orgy+(ny-1)*v2len
	end if
	do iy=1,ny
		curvepos(iy)=orgy+(iy-1)*v2len
		if (curvepos(iy)<rlow.or.curvepos(iy)>rhigh) cycle
        call vec2area(gridv1(:),gridv3(:),area)
		locintcurve(iy)=sum(cubmat(:,iy,:))*area
        pleavgcurve(iy)=sum(cubmat(:,iy,:))/(nx*nz)
		tmpintval=tmpintval+locintcurve(iy)*v2len
		intcurve(iy)=tmpintval
	end do
else if (idir==3) then !Direction 3
	allocate(intcurve(nz),locintcurve(nz),curvepos(nz),pleavgcurve(nz))
	ncurpt=nz
    v3len=dsqrt(sum(gridv3**2))
	if (c200tmp==" ".or.index(c200tmp,"a")/=0) then
		rlow=orgz
		rhigh=orgz+(nz-1)*v3len
	end if
	do iz=1,nz
		curvepos(iz)=orgz+(iz-1)*v3len
		if (curvepos(iz)<rlow.or.curvepos(iz)>rhigh) cycle
        call vec2area(gridv1(:),gridv2(:),area)
		locintcurve(iz)=sum(cubmat(:,:,iz))*area
        pleavgcurve(iz)=sum(cubmat(:,:,iz))/(nx*ny)
		tmpintval=tmpintval+locintcurve(iz)*v3len
		intcurve(iz)=tmpintval
	end do
end if

write(*,"(/,a)") " Options 1~3 are mainly used for previewing purpose. To gain better graphical effect, &
&you can use options 7~9 to export curve data and plot the curves using third-part tool, such as gnuplot"
write(*,*)
write(*,"(' Minimum and maximum of plane-averaged curve:',2E16.8)") minval(pleavgcurve),maxval(pleavgcurve)
write(*,"(' Minimum and maximum of local integral curve:',2E16.8)") minval(locintcurve),maxval(locintcurve)
write(*,"(' Minimum and maximum of integral curve:      ',2E16.8)") minval(intcurve),maxval(intcurve)
rlowlim=1
ruplim=1

do while(.true.)
	write(*,*)
	write(*,*) "-2 Set format of saving graphical file, current: "//graphformat
	if (ilenunit1D==1) write(*,*) "-1 Change length unit of the graph from Bohr to Angstrom"
	if (ilenunit1D==2) write(*,*) "-1 Change length unit of the graph from Angstrom to Bohr"
	write(*,*) "0 Return"
	write(*,*) "1 Plot graph of integral curve"
	write(*,*) "2 Plot graph of local integral curve"
	write(*,*) "3 Plot graph of plane-averaged curve"
	write(*,*) "4 Save graph of integral curve to current folder"
	write(*,*) "5 Save graph of local integral curve to current folder"
	write(*,*) "6 Save graph of plane-averaged curve to current folder"
	write(*,*) "7 Export data of integral curve to intcurve.txt in current folder"
	write(*,*) "8 Export data of local integral curve to locintcurve.txt in current folder"
	write(*,*) "9 Export data of plane-averaged curve to pleavgcurve.txt in current folder"
	if (rlowlim==ruplim) then
		write(*,*) "10 Set lower and upper limits of Y-axis, current: Auto"
    else
		write(*,"(a,2E14.6)") " 10 Set lower and upper limits of Y-axis, current:",rlowlim,ruplim
    end if
    write(*,*) "11 Calculate value a given position (via linear interpolation)"
	read(*,*) isel2
    if (isel2==-2) then
        call setgraphformat
	else if (isel2==-1) then
		if (ilenunit1D==1) then
			ilenunit1D=2
		else if (ilenunit1D==2) then
			ilenunit1D=1
		end if
	else if (isel2==0) then
		exit
	else if (isel2==1.or.isel2==4) then !Integral curve
        if (rlowlim==ruplim) then
			disminmax=maxval(intcurve)-minval(intcurve)
			ylow=minval(intcurve)-0.1D0*disminmax
			yhigh=maxval(intcurve)+0.1D0*disminmax
        else
			ylow=rlowlim
            yhigh=ruplim
        end if
		stplabx=(rhigh-rlow)/10
		stplaby=(yhigh-ylow)/10
        call showcurveminmax(ncurpt,curvepos,intcurve,ilenunit1D)
		if (isel2==1) then
			call drawcurve(curvepos,intcurve,ncurpt,rlow,rhigh,stplabx,ylow,yhigh,stplaby,"show")
		else if (isel2==4) then
			call drawcurve(curvepos,intcurve,ncurpt,rlow,rhigh,stplabx,ylow,yhigh,stplaby,"save")
            write(*,"(a)") " Done! Image file of integral curve has been saved to current folder with DISLIN prefix"
        end if
	else if (isel2==2.or.isel2==5) then !Local integral curve
        if (rlowlim==ruplim) then
			disminmax=maxval(locintcurve)-minval(locintcurve)
			ylow=minval(locintcurve)-0.1D0*disminmax
			yhigh=maxval(locintcurve)+0.1D0*disminmax
        else
			ylow=rlowlim
            yhigh=ruplim
        end if
		stplabx=(rhigh-rlow)/10
		stplaby=(yhigh-ylow)/10
        call showcurveminmax(ncurpt,curvepos,locintcurve,ilenunit1D)
		if (isel2==2) then
			call drawcurve(curvepos,locintcurve,ncurpt,rlow,rhigh,stplabx,ylow,yhigh,stplaby,"show")
		else if (isel2==5) then
			call drawcurve(curvepos,locintcurve,ncurpt,rlow,rhigh,stplabx,ylow,yhigh,stplaby,"save")
            write(*,"(a)") " Done! Image file of local integral curve has been saved to current folder with DISLIN prefix"
        end if
	else if (isel2==3.or.isel2==6) then !Plane-averaged curve
        if (rlowlim==ruplim) then
			disminmax=maxval(pleavgcurve)-minval(pleavgcurve)
			ylow=minval(pleavgcurve)-0.1D0*disminmax
			yhigh=maxval(pleavgcurve)+0.1D0*disminmax
        else
			ylow=rlowlim
            yhigh=ruplim
        end if
		stplabx=(rhigh-rlow)/10
		stplaby=(yhigh-ylow)/10
        call showcurveminmax(ncurpt,curvepos,pleavgcurve,ilenunit1D)
		if (isel2==3) then
			call drawcurve(curvepos,pleavgcurve,ncurpt,rlow,rhigh,stplabx,ylow,yhigh,stplaby,"show")
		else if (isel2==6) then
			call drawcurve(curvepos,pleavgcurve,ncurpt,rlow,rhigh,stplabx,ylow,yhigh,stplaby,"save")
            write(*,"(a)") " Done! Image file of plane-averaged curve has been saved to current folder with DISLIN prefix"
        end if
	else if (isel2==7) then
		open(10,file="intcurve.txt",status="replace")
		do ipt=1,ncurpt
			if (curvepos(ipt)<rlow.or.curvepos(ipt)>rhigh) cycle
			write(10,"(2f14.8,f20.10)") curvepos(ipt),curvepos(ipt)*b2a,intcurve(ipt)
		end do
		close(10)
		write(*,"(a)") " Done! 1,2,3 columns correspond to the coordinate (in Bohr and in Angstrom) in the direction you selected and the integral, respectively"
	else if (isel2==8) then
		open(10,file="locintcurve.txt",status="replace")
		do ipt=1,ncurpt
			if (curvepos(ipt)<rlow.or.curvepos(ipt)>rhigh) cycle
			write(10,"(2f14.8,f20.10)") curvepos(ipt),curvepos(ipt)*b2a,locintcurve(ipt)
		end do
		close(10)
		write(*,"(a)") " Done! 1,2,3 columns correspond to the coordinate (in Bohr and in Angstrom) in the direction you selected and the local integral, respectively"
	else if (isel2==9) then
		open(10,file="pleavgcurve.txt",status="replace")
		do ipt=1,ncurpt
			if (curvepos(ipt)<rlow.or.curvepos(ipt)>rhigh) cycle
			write(10,"(2f14.8,f20.10)") curvepos(ipt),curvepos(ipt)*b2a,pleavgcurve(ipt)
		end do
		close(10)
		write(*,"(a)") " Done! 1,2,3 columns correspond to the coordinate (in Bohr and in Angstrom) in the direction you selected and the plane-averaged value, respectively"						
	else if (isel2==10) then
		write(*,*) "Input lower and upper limits of Y-axis, e.g. -2.3,16.8"
        write(*,*) "If pressing ENTER button directly, the limits will be automatically determined"
        read(*,"(a)") c80tmp
        if (c80tmp==" ") then
			rlowlim=1
            ruplim=1
        else
	        read(c80tmp,*) rlowlim,ruplim
        end if
    else if (isel2==11) then
		write(*,*) "Calculate value a given point for which curve?"
        write(*,*) "1 Integral curve"
        write(*,*) "2 Local integral curve"
        write(*,*) "3 Plane-averaged curve"
        read(*,*) itype
        if (ilenunit1D==1) write(*,*) "Input position in Bohr, e.g. 2.3"
        if (ilenunit1D==2) write(*,*) "Input position in Angstrom, e.g. 2.3"
        read(*,*) posinp
        if (ilenunit1D==2) posinp=posinp/b2a
        if (posinp<minval(curvepos)) then
			write(*,*) "Error: The inputted position is smaller than lower limit!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            cycle
        else if (posinp>maxval(curvepos)) then
			write(*,*) "Error: The inputted position is larger than upper limit!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            cycle
        end if
        if (itype==1) call linintpol(curvepos,intcurve,ncurpt,posinp,val)
        if (itype==2) call linintpol(curvepos,locintcurve,ncurpt,posinp,val)
        if (itype==3) call linintpol(curvepos,pleavgcurve,ncurpt,posinp,val)
        write(*,"(' The value is ',1PE18.8)") val
	end if
end do
end subroutine