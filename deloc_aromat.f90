!---------------------------------------------------------------------
!-------- Interface of various delocalization and aromaticity analyses
!---------------------------------------------------------------------
subroutine deloc_aromat
implicit real*8 (a-h,o-z)
do while(.true.)
	write(*,*)
	write(*,*) "   ================ Delocalization and aromaticity analyses ==============="
    write(*,*) " 0 Return to main menu"
    write(*,*) " 1 Multi-center index (MCI), also known as multicenter bond order (MCBO)"
    write(*,*) "-1 Multi-center index (MCI) in NAO basis"
    write(*,*) " 2 AV1245 index"
    write(*,*) " 3 Iso-chemical shielding surface (ICSS)"
    write(*,*) " 4 NICS_ZZ for non-planar or tilted system"
    write(*,*) " 5 ELF-sigma/pi and LOL-sigma/pi"
    write(*,*) " 6 Harmonic oscillator measure of aromaticity (HOMA) and Bird indices"
    write(*,*) " 7 Shannon aromaticity index"
    write(*,*) " 8 Para-delocalization index (PDI)"
    write(*,*) " 9 Aromatic fluctuation index (FLU) and FLU-pi"
    write(*,*) "10 Para linear response index (PLR)"
    write(*,*) "11 Information-theoretic (ITA) aromaticity index"
    write(*,*) "12 Properties of ring critical point"
    read(*,*) isel
    if (isel==0) then
		return
    else if (isel==1) then
        call ask_Sbas_PBC
		call multicenter(2)
	else if (isel==-1) then
		call multicenterNAO
	else if (isel==2) then
		call AV1245
	else if (isel==3) then
		call ICSS
	else if (isel==4) then
		call NICS_ZZ
	else if (isel==5) then
		write(*,"(a)") " Multiwfn is able to analyze ELF-sigma/pi and LOL-sigma/pi in different ways, such as plotting plane and isosurface maps, &
        performing topology analysis to obtain bifurcation value, etc., which rely on different main functions. Please check &
        Section 4.4.9, 4.5.3 and 4.100.22 of manunal for example of realizing these analyses."
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else if (isel==6) then
		call HOMA_Bird
	else if (isel==7) then
		write(*,"(a)") " To realize this analysis, you should use topology analysis module (main function 2), see Section 3.14.6 of manual &
        for introduction and Section 4.2.1 for practical example."
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else if (isel==8) then
		write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.6 of manual &
        for introduction and Section 4.15.2 for practical example."
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else if (isel==9) then
		write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.7 of manual &
        for introduction and Section 4.15.2 for practical example."
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else if (isel==10) then
		write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.9 of manual &
        for introduction and Section 4.15.2 for practical example."
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else if (isel==11) then
		write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.11 of manual for detail."
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else if (isel==12) then
		write(*,"(a)") " To realize this analysis, you should use topology analysis module (main function 2), see Section 3.14.6 of manual &
        for introduction and Section 4.2.1 for practical example."
        write(*,*) "Press ENTER button to continue"
        read(*,*)
    end if
end do
end subroutine




!------------------------------------------------------------
!---------- Iso-chemical shielding surfaces (ICSS) ----------
!------------------------------------------------------------
subroutine ICSS
use defvar
use util
use GUI
implicit real*8 (a-h,o-z)
character c200tmp*200,gauinpfile*200,gauoutfile*200,selectyn,suffix*4
character,allocatable :: gauinpcontent(:)*79

write(*,*) "Citation of ICSS: J. Chem. Soc. Perkin Trans. 2, 2001, 1893"
write(*,*) "Citation of ICSS_XX/YY/ZZ and ICSS_ani:"
write(*,*) "Carbon, 165, 468 (2020) DOI: 10.1016/j.carbon.2020.04.099"

!Set grid for calculating NICS
aug3D=12
call setgrid(0,itmp)
numbqper=NICSnptlim-ncenter
write(*,"(' The number of Bq per batch:',i10)") numbqper
write(*,"(' The number of center per file (NICSnptlim in settings.ini):',i10)") NICSnptlim
npttot=nx*ny*nz
nfile=ceiling(dfloat(npttot)/numbqper)
!Generate Gaussian input file
write(*,*)
write(*,*) "If skip generating Gaussian input file of NMR task? (y/n)"
read(*,*) selectyn
if (selectyn=='n'.or.selectyn=='N') then
    if (ifiletype==12) then
        write(*,"(a)") " Note: The file loaded when Multiwfn boots up is regarded as template Gaussian input file"
        gauinpfile=filename
    else
	    write(*,*) "Input the path of template Gaussian input file, e.g. C:\ltwd.gjf"
	    do while(.true.)
		    read(*,"(a)") gauinpfile
		    inquire(file=gauinpfile,exist=alive)
		    if (alive) exit
		    write(*,*) "Error: Cannot find corresponding files, input again"
	    end do
    end if
	open(10,file=gauinpfile,status="old")
	numgauline=totlinenum(10,2)
	allocate(gauinpcontent(numgauline))
	numblank=0
	iendcoord=numgauline !Which line is the last line recording coordinates
	do i=1,numgauline
		read(10,"(a)") gauinpcontent(i)
		if (gauinpcontent(i)==" ") then
			numblank=numblank+1
			if (numblank==3) iendcoord=i-1
		end if
		if (index(gauinpcontent(i),'#')/=0) then
			gauinpcontent(i)=trim(gauinpcontent(i))//" NMR"
			if (index(gauinpcontent(i),'conn')==0) gauinpcontent(i)=trim(gauinpcontent(i))//" geom=connectivity"
		end if
	end do
	close(10)

	gauinpfile="NICS"
	do ifile=1,nfile
		write(c200tmp,"(a,i4.4,a)") trim(gauinpfile),ifile,".gjf"
		open(10,file=c200tmp,status="replace")
		write(*,"(a,a,a)") " Outputting ",trim(c200tmp)," to current folder..."
		!Move keywords and coordinate from template file into NMR input file
		do i=1,iendcoord
			if (ifile>1.and.index(gauinpcontent(i),'#')/=0) then
				write(10,"(a)") trim(gauinpcontent(i))//" guess=read"
			else
				write(10,"(a)") trim(gauinpcontent(i))
			end if
		end do
		!Write Bq information
		itmp=0
		do i=1,nx
			do j=1,ny
				do k=1,nz
					itmp=itmp+1
                    call getgridxyz(i,j,k,rnowx,rnowy,rnowz)
					if (itmp<=(ifile-1)*numbqper) cycle
					if (itmp>ifile*numbqper) exit
					write(10,"('Bq ',3f12.6)") rnowx*b2a,rnowy*b2a,rnowz*b2a
				end do
			end do
		end do
		write(10,*)
		!Write connectivity explicitly
		if (ifile/=nfile) then
			do i=1,NICSnptlim
				write(10,"(i7)") i
			end do
		else if (ifile==nfile) then
			do i=1,npttot-(ifile-1)*numbqper+ncenter
				write(10,"(i7)") i
			end do
		end if
		!Write remaining part
		do i=iendcoord+1,numgauline
			write(10,"(a)") gauinpcontent(i)
		end do
		close(10)
	end do
end if
write(*,*) "Now please run these input files by Gaussian"
write(*,*)
!Load NICS from Gaussian output file
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
write(*,*) "Input the path of Gaussian output file of NMR task"
write(*,"(a)") " Assume that you input ""C:\ltwd\NICS"", then C:\ltwd\NICS0001.out, C:\ltwd\NICS0002.out, &
C:\ltwd\NICS0003.out... will be loaded (.log suffix is also allowed)"
do while(.true.)
	read(*,"(a)") gauoutfile
	suffix=".out"
	inquire(file=trim(gauoutfile)//"0001"//suffix,exist=alive)
	if (alive) exit
	if (.not.alive) then
		suffix=".log"
		inquire(file=trim(gauoutfile)//"0001"//suffix,exist=alive)
	end if
	if (alive) exit
	write(*,"(a)") " Error: Unable to find either "//trim(gauoutfile)//"0001.out or "//trim(gauoutfile)//"0001.log"
	write(*,*) "Please input the path again"
end do
100 write(*,*) "Load which term?"
write(*,*) "1: Isotropic  2: Anisotropy  3: XX component  4: YY component  5: ZZ component"
read(*,*) iload
do ifile=1,nfile
	write(c200tmp,"(a,i4.4,a)") trim(gauoutfile),ifile,suffix
	inquire(file=c200tmp,exist=alive)
	if (.not.alive) then
		write(*,"(' Error: Unable to find ',a)") trim(c200tmp)
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		return
	end if
	write(*,"(' Loading ',a,'...')") trim(c200tmp)
	open(10,file=c200tmp,status="old")
	call loclabel(10,"GIAO Magnetic shielding tensor",ifound)
    if (ifound==0) then
		write(*,*) "Error: Unable to find ""GIAO Magnetic shielding tensor"" field!"
        write(*,*) "Please check this file to make sure that the task has normally finished"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        return
    end if
	read(10,*)
	!Detect format. The NMR output format changes since G09 D.01 to leave more space for atomic index
	read(10,"(a80)") c200tmp
	backspace(10)
	iformat=1
	if (c200tmp(25:25)=='=') iformat=2 !Since G09 D.01
	
	do i=1,ncenter !Skip atom's result
		read(10,*)
		read(10,*)
		read(10,*)
		read(10,*)
		read(10,*)
	end do
	itmp=0
	iloadthis=0
	do i=1,nx
		do j=1,ny
			do k=1,nz
				itmp=itmp+1
				if (itmp<=(ifile-1)*numbqper) cycle
				if (itmp>ifile*numbqper) exit
				if (iload==1.or.iload==2) then
					iloadthis=iloadthis+1
                    if (iformat==1) then
					    if (iload==1) read(10,"(22x,f10.4)",iostat=ierror) cubmat(i,j,k)
					    if (iload==2) read(10,"(48x,f10.4)",iostat=ierror) cubmat(i,j,k)
					else if (iformat==2) then
					    if (iload==1) read(10,"(26x,f10.4)",iostat=ierror) cubmat(i,j,k)
					    if (iload==2) read(10,"(52x,f10.4)",iostat=ierror) cubmat(i,j,k)
					end if
					if (ierror/=0) then
						write(*,"(' Error: Unable to load the',i7,'th Bq in this file!')") iloadthis
						write(*,"(' This Bq should correspond to the',i7,'th center in this file')") ncenter+iloadthis
						write(*,*) "Please double check your grid setting and ""NICSnptlim"" in settings.ini"
                        write(*,*) "Press ENTER button to exit"
						read(*,*)
						return
					end if
					read(10,*)
					read(10,*)
					read(10,*)
					read(10,*)
				else
					read(10,*)
					if (iload==3) then
						read(10,"(8x,f10.4)") cubmat(i,j,k)
						read(10,*)
						read(10,*)
					else if (iload==4) then
						read(10,*)
						read(10,"(24x,f10.4)") cubmat(i,j,k)
						read(10,*)
					else if (iload==5) then
						read(10,*)
						read(10,*)
						read(10,"(42x,f10.4)") cubmat(i,j,k)
					end if
					read(10,*)
				end if
			end do
		end do
	end do
	close(10)
end do
write(*,*) "Loading finished!"
do while(.true.)
	write(*,*)
	write(*,*) "-1 Load another ICSS form"
	write(*,*) "0 Return to main menu"
	if (iload==1) write(*,*) "1 Visualize iso-chemical shielding surface"
	if (iload==2) write(*,*) "1 Visualize aniso-chemical shielding surface"
	if (iload==3) write(*,*) "1 Visualize XX-chemical shielding surface"
	if (iload==4) write(*,*) "1 Visualize YY-chemical shielding surface"
	if (iload==5) write(*,*) "1 Visualize ZZ-chemical shielding surface"
	if (iload==1) write(*,*) "2 Export the grid data to ICSS.cub current folder"
	if (iload==2) write(*,*) "2 Export the grid data to AICSS.cub current folder"
	if (iload==3) write(*,*) "2 Export the grid data to ICSSXX.cub current folder"
	if (iload==4) write(*,*) "2 Export the grid data to ICSSYY.cub current folder"
	if (iload==5) write(*,*) "2 Export the grid data to ICSSZZ.cub current folder"
	read(*,*) isel
	if (isel==-1) then
		goto 100
	else if (isel==0) then
		exit
	else if (isel==1) then
		call drawisosurgui(1)
	else if (isel==2) then
		if (iload==1) open(10,file="ICSS.cub",status="replace")
		if (iload==2) open(10,file="AICSS.cub",status="replace")
		if (iload==3) open(10,file="ICSSXX.cub",status="replace")
		if (iload==4) open(10,file="ICSSYY.cub",status="replace")
		if (iload==5) open(10,file="ICSSZZ.cub",status="replace")
		call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,"(a)") " The cube file has been exported to current folder"
	end if
end do
end subroutine




!---------------------------------------------
!------- HOMA / Bird aromaticity index -------
!---------------------------------------------
subroutine HOMA_Bird
use defvar
use util
implicit real*8 (a-h,o-z)
integer :: numHOMAatm=6,HOMAatm(1000),numBirdatm=5,Birdatm(1000)
real*8 :: HOMArefbond(nelesupp,nelesupp)=-1D0,HOMAsigma(nelesupp,nelesupp)=0D0
real*8 :: Birda(nelesupp,nelesupp)=-1D0,Birdb(nelesupp,nelesupp)=-1D0
real*8 :: BirdVref(100)=-1
real*8,allocatable :: BirdN(:)
character C200inp*200

if (all(HOMArefbond==-1D0)) then !If =-1, means missing reference value
	HOMArefbond(5,6)=1.4235D0
	HOMArefbond(6,5)=HOMArefbond(5,6)
	HOMArefbond(5,7)=1.402D0
	HOMArefbond(7,5)=HOMArefbond(5,7)
	HOMArefbond(6,6)=1.388D0
	HOMArefbond(6,7)=1.334D0
	HOMArefbond(7,6)=HOMArefbond(6,7)
	HOMArefbond(6,8)=1.265D0
	HOMArefbond(8,6)=HOMArefbond(6,8)
	HOMArefbond(6,15)=1.698D0
	HOMArefbond(15,6)=HOMArefbond(6,15)
	HOMArefbond(6,16)=1.677D0
	HOMArefbond(16,6)=HOMArefbond(6,16)
	HOMArefbond(7,7)=1.309D0
	HOMArefbond(7,8)=1.248D0
	HOMArefbond(8,7)=HOMArefbond(7,8)
	HOMAsigma(5,6)=104.507D0
	HOMAsigma(6,5)=HOMAsigma(5,6)
	HOMAsigma(5,7)=72.03D0
	HOMAsigma(7,5)=HOMAsigma(5,7)
	HOMAsigma(6,6)=257.7D0
	HOMAsigma(6,7)=93.52D0
	HOMAsigma(7,6)=HOMAsigma(6,7)
	HOMAsigma(6,8)=157.38D0
	HOMAsigma(8,6)=HOMAsigma(6,8)
	HOMAsigma(6,15)=118.91D0
	HOMAsigma(15,6)=HOMAsigma(6,15)
	HOMAsigma(6,16)=94.09D0
	HOMAsigma(16,6)=HOMAsigma(6,16)
	HOMAsigma(7,7)=130.33D0
	HOMAsigma(7,8)=57.21D0
	HOMAsigma(8,7)=HOMAsigma(7,8)
end if
if (all(Birda==-1D0)) then !If =-1, means missing reference value
	Birda(5,7)=7.72D0
	Birda(7,5)=Birda(5,7)
	Birdb(5,7)=2.16D0
	Birdb(7,5)=Birdb(5,7)
	Birda(6,6)=6.80D0
	Birdb(6,6)=1.71D0
	Birda(6,7)=6.48D0
	Birda(7,6)=Birda(6,7)
	Birdb(6,7)=2.00D0
	Birdb(7,6)=Birdb(6,7)
	Birda(6,8)=5.75D0
	Birda(8,6)=Birda(6,8)
	Birdb(6,8)=1.85D0
	Birdb(8,6)=Birdb(6,8)
	Birda(6,16)=11.9D0
	Birda(16,6)=Birda(6,16)
	Birdb(6,16)=2.59D0
	Birdb(16,6)=Birdb(6,16)
	Birda(7,8)=4.98D0
	Birda(8,7)=Birda(7,8)
	Birdb(7,8)=1.41D0
	Birdb(8,7)=Birdb(7,8)
	Birda(7,7)=5.28D0
	Birdb(7,7)=1.41D0
end if
if (all(BirdVref==-1D0)) then
	BirdVref(5)=35D0
	BirdVref(6)=33.2D0	
end if

do while(.true.)
	write(*,*)
	write(*,*) "           ------------- HOMA / Bird aromaticity index -------------"
	write(*,*) "-1 Return"
	write(*,*) "0 Start calculation for HOMA!"
	write(*,*) "1 Adjust parameters for HOMA calculation"
	write(*,*) "2 Start calculation for Bird aromaticity index!"
	write(*,*) "3 Adjust a and b parameters for Bird aromaticity index calculation"
	write(*,*) "4 Adjust reference V parameter for Bird aromaticity index calculation"
	read(*,*) isel
	if (isel==-1) then
		exit
	else if (isel==0) then
		write(*,*) "Current reference bond length (Angstrom) and sigma parameters:"
		do iref=1,nelesupp
			do jref=iref,nelesupp
				if (HOMArefbond(iref,jref)/=-1) write(*,"(' ',a,a,a,a,2f12.4)") ind2name(iref),'-',ind2name(jref),':',HOMArefbond(iref,jref),HOMAsigma(iref,jref)
			end do
		end do
		write(*,*)
		do while(.true.)
			write(*,*) "Input indices of the atoms involved in the ring, e.g. 1,5,6,7,8,12"
			write(*,*) "(Input q can return)"
			read(*,"(a)") c200inp
			if (c200inp=='q'.or.c200inp=='Q') exit
			call str2arr(c200inp,numHOMAatm,HOMAatm)
			HOMAval=1D0
            write(*,*)
			write(*,*) "        Atom pair         Contribution  Bond length(Angstrom)"
			do iidx=1,numHOMAatm
				jidx=iidx+1
				if (iidx==numHOMAatm) jidx=1
				iatm=HOMAatm(iidx) !Actual atom index in present system
				jatm=HOMAatm(jidx)
				iatmeleidx=a(iatm)%index !Index in periodic table
				jatmeleidx=a(jatm)%index
				refbondlen=HOMArefbond(iatmeleidx,jatmeleidx)
				refsigma=HOMAsigma(iatmeleidx,jatmeleidx)
				if (refbondlen==-1D0) then
					write(*,"(' Error: Missing reference parameter for ',a,'-',a)") ind2name(iatmeleidx),ind2name(jatmeleidx)
					exit
				end if
				paircontri=-refsigma/numHOMAatm*(refbondlen-distmat(iatm,jatm)*b2a)**2
				write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),paircontri,distmat(iatm,jatm)*b2a
				HOMAval=HOMAval+paircontri
				if (iidx==numHOMAatm) write(*,"(a,f12.6)") " HOMA value is",HOMAval
			end do
			write(*,*)
		end do
	else if (isel==1) then
		write(*,*) "Current reference bond length (in Angstrom) and sigma paramters:"
		do iref=1,nelesupp
			do jref=iref,nelesupp
				if (HOMArefbond(iref,jref)/=-1) write(*,"(' ',a,a,a,a,2f12.5)") ind2name(iref),'- ',ind2name(jref),':',HOMArefbond(iref,jref),HOMAsigma(iref,jref)
			end do
		end do
		do while(.true.)
			write(*,*) "Input two element indices and new bond length and sigma parameter"
			write(*,"(a)") " e.g. 6,7,1.334,93.52 means set reference bond length and sigma for C-N to 1.334 Angstrom and 93.52 respectively"
			write(*,*) "(Input q can return)"
			read(*,"(a)") C200inp
			if (C200inp(1:1)=='q'.or.C200inp(1:1)=='Q') exit
			read(C200inp,*) itmp,jtmp,refbondtmp,sigmatmp
			HOMArefbond(itmp,jtmp)=refbondtmp
			HOMArefbond(jtmp,itmp)=refbondtmp
			HOMAsigma(itmp,jtmp)=sigmatmp
			HOMAsigma(jtmp,itmp)=sigmatmp
			write(*,*) "Done!"
		end do
	
	else if (isel==2) then
	    write(*,*) "Note: Default a and b parameters are taken from Tetrahedron, 57, 5715 (2001)"
		write(*,*) "Current a and b parameters:"
		do iref=1,nelesupp
			do jref=iref,nelesupp
				if (Birda(iref,jref)/=-1) write(*,"(' ',a,a,a,a,2f12.5)") ind2name(iref),'- ',ind2name(jref),':',Birda(iref,jref),Birdb(iref,jref)
			end do
		end do
		write(*,*)
		write(*,*) "Current reference V parameters:"	
		do itmp=1,size(BirdVref)
			if (BirdVref(itmp)==-1) cycle
			write(*,"(i4,' centers, value:',f10.4)") itmp,BirdVref(itmp)
		end do
		write(*,*)
		do while(.true.)
			write(*,*) "Input indices of the atoms involved, e.g. 1,5,6,7,8,12"
			write(*,*) "(Input q can return)"
			read(*,"(a)") c200inp
			if (c200inp=='q'.or.c200inp=='Q') exit
			call str2arr(c200inp,numBirdatm,Birdatm)
			if (BirdVref(numBirdatm)==-1) then
				write(*,*) "Error: Missing reference V parameter for this number of centers!"
				exit
			end if
			allocate(BirdN(numBirdatm))
			write(*,*) "        Atom pair             N term   Bond length(Angstrom)"
			do iidx=1,numBirdatm
				jidx=iidx+1
				if (iidx==numBirdatm) jidx=1
				iatm=Birdatm(iidx) !Actual atom index in present system
				jatm=Birdatm(jidx)
				iatmeleidx=a(iatm)%index !Index in periodic table
				jatmeleidx=a(jatm)%index
				Birdanow=Birda(iatmeleidx,jatmeleidx)
				Birdbnow=Birdb(iatmeleidx,jatmeleidx)
				if (Birdanow==-1D0) then
					write(*,"(' Error: Missing a and b parameters for ',a,'-',a)") ind2name(iatmeleidx),ind2name(jatmeleidx)
					exit
				end if
				BirdN(iidx)=Birdanow/(distmat(iatm,jatm)*b2a)**2-Birdbnow
				write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),BirdN(iidx),distmat(iatm,jatm)*b2a
				if (iidx==numBirdatm) then
					avgBirdN=sum(BirdN)/numBirdatm
					BirdVnow=dsqrt(sum((BirdN(:)-avgBirdN)**2)/numBirdatm)*100/avgBirdN
					Birdval=100*(1-BirdVnow/BirdVref(numBirdatm))
					write(*,"(a,f12.6)") " Bird aromaticity index is",Birdval
				end if
			end do
			deallocate(BirdN)
			write(*,*)
		end do
	else if (isel==3) then
		write(*,*) "Current a and b parameters:"
		do iref=1,nelesupp
			do jref=iref,nelesupp
				if (Birda(iref,jref)/=-1) write(*,"(' ',a,a,a,a,2f12.5)") ind2name(iref),'- ',ind2name(jref),':',Birda(iref,jref),Birdb(iref,jref)
			end do
		end do
		do while(.true.)
			write(*,*) "Input two element indices and a and b parameters"
			write(*,"(a)") " e.g. 6,7,6.941,2.205 means set a and b parameters for C-N to 6.941 and 2.205 respectively"
			write(*,*) "(Input q can return)"
			read(*,"(a)") C200inp
			if (C200inp(1:1)=='q'.or.C200inp(1:1)=='Q') exit
			read(C200inp,*) itmp,jtmp,Birdatmp,Birdbtmp
			Birda(itmp,jtmp)=Birdatmp
			Birda(jtmp,itmp)=Birda(itmp,jtmp)
			Birdb(itmp,jtmp)=Birdbtmp
			Birdb(jtmp,itmp)=Birdb(itmp,jtmp)
			write(*,*) "Done!"
		end do
	else if (isel==4) then
		write(*,*) "Current reference V parameters:"	
		do itmp=1,size(BirdVref)
			if (BirdVref(itmp)==-1) cycle
			write(*,"(i4,' centers, value:',f10.4)") itmp,BirdVref(itmp)
		end do
		write(*,*) "Set the parameter for how many centers? e.g. 5"
		read(*,*) nBirdcentmp
		write(*,*) "Set to which value? e.g. 33.4"
		read(*,*) BirdVref(nBirdcentmp)
		write(*,*) "Done!"
	end if
end do
end subroutine




!!-----------------------------------------
!!------------ AV1245 index ---------------
!!-----------------------------------------
subroutine AV1245
use defvar
use util
use NAOmod
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer,allocatable :: atmarr(:),atmarrorg(:)
integer cenind(2000),minidx(4)
real*8,allocatable :: PSmat(:,:),PSmatA(:,:),PSmatB(:,:)

iopsh=0
if (allocated(CObasa)) then !Calculate AV1245 in original basis
    write(*,*) "Calculating PS matrix, please wait..."
    if (allocated(Palpha)) then !Open shell
        iopsh=1
        allocate(PSmatA(nbasis,nbasis),PSmatB(nbasis,nbasis))
        PSmatA=matmul_blas(Palpha,Sbas,nbasis,nbasis,0,0)
        PSmatB=matmul_blas(Pbeta,Sbas,nbasis,nbasis,0,0)
    else
        allocate(PSmat(nbasis,nbasis))
        PSmat=matmul_blas(Ptot,Sbas,nbasis,nbasis,0,0)
    end if
    ifNAO=0
else !Load NAO and DMNAO information
    write(*,"(a)") " Basis information is not presented, therefore trying to load natural atomic orbital (NAO) information from input file"
    open(10,file=filename,status="old")
    call checkNPA(ifound);if (ifound==0) return
    call loadNAOinfo
    write(*,*) "Loading NAO information finished!"
    call checkDMNAO(ifound);if (ifound==0) return
    call loadDMNAO
    close(10)
    write(*,*) "Loading density matrix in NAO basis finished!"
    write(*,*) "The AV1245 will be calculated based on NAOs"
    if (iopshNAO==0) then
        allocate(PSmat(numNAO,numNAO))
        PSmat=DMNAO
    else if (iopshNAO==1) then !Open shell
        iopsh=1
        allocate(PSmatA(numNAO,numNAO),PSmatB(numNAO,numNAO))
        PSmatA=DMNAOa
        PSmatB=DMNAOb
    end if
    nbasis=numNAO
    ifNAO=1
    !Move information from NAO variables to common variables, so that multi-center bond order routines could be used
    if (allocated(basstart)) deallocate(basstart,basend)
    allocate(basstart(ncenter),basend(ncenter))
    basstart=NAOinit
    basend=NAOend
end if
iMCBOtype_old=iMCBOtype
iMCBOtype=2

do while(.true.)
    write(*,*)
    write(*,*) "                      --------- AV1245 and AVmin ---------"
    write(*,*) "Input index of the atoms in the order of connectivity, e.g. 2,3,7,18,19,20"
    write(*,*) "To exit, input ""q"""
    !When NAO information is loaded form NBO output file, geometry information is not available and cannot generate connectivity
    if (ifNAO==0) write(*,"(a)") " Hint: If input ""d"" and press ENTER button, then you can input the indices in arbitrary order because the actual order &
    will be automatically guessed, however in this case any atom should not connect to more than two atoms in the ring"
    read(*,"(a)") c2000tmp
    
    if (index(c2000tmp,'q')/=0) then
        exit
    else if (index(c2000tmp,'d')/=0) then
        if (.not.allocated(connmat)) call genconnmat(1,0) !Generate connectivity matrix
        write(*,*)
        write(*,*) "Input index of the atoms, the order is arbitrary"
        write(*,*) "For example: 1,3-4,6-8,10-14"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,natm)
        allocate(atmarr(natm),atmarrorg(natm))
        call str2arr(c2000tmp,natm,atmarrorg)
        !Reorganize the atmarrorg to correct sequence as atmarr according to connectivity
        atmarr=0
        atmarr(1)=atmarrorg(1)
        inow=atmarr(1) !Current atom
        atmarrorg(1)=0 !This atom has been picked out, so set to zero
        do idx=2,natm
            do jdx=1,natm
                if (atmarrorg(jdx)==0) cycle
                jatm=atmarrorg(jdx)
                if (connmat(inow,jatm)/=0) then
                    inow=jatm
                    atmarr(idx)=inow
                    atmarrorg(jdx)=0
                    exit
                end if
            end do
            if (jdx==natm+1) then
                write(*,"(' Failed to determine connectivity of atom',i6)") inow
                exit
            end if
        end do
        deallocate(atmarrorg)
        if (any(atmarr<=0)) then
            write(*,"(a)") " Unfortunately, the order was not successfully recognized, you should manually input &
            the atom indices according to connectivity"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            deallocate(atmarr)
            cycle
        else
            write(*,*) "The order of the atoms in the ring has been successfully identified"
            write(*,*)
        end if
    else
        call str2arr(c2000tmp,natm)
        allocate(atmarr(natm))
        call str2arr(c2000tmp,natm,atmarr)
    end if
    
    write(*,"(' Number of selected atoms:',i6)") natm
    write(*,*) "Atomic sequence:"
    write(*,"(12i6)") atmarr
    write(*,*)
    totval=0
    ipos=1
    AVmin=1D10
    do while(.true.)
        cenind(1)=atmarr(ipos)
        if (ipos+1>natm) then
            cenind(2)=atmarr(ipos+1-natm)
        else
            cenind(2)=atmarr(ipos+1)
        end if
        if (ipos+3>natm) then
            cenind(3)=atmarr(ipos+3-natm)
        else
            cenind(3)=atmarr(ipos+3)
        end if
        if (ipos+4>natm) then
            cenind(4)=atmarr(ipos+4-natm)
        else
            cenind(4)=atmarr(ipos+4)
        end if
        if (iopsh==0) then
            call calcmultibndord(4,cenind,PSmat,nbasis,tmpval)
        else
            call calcmultibndord(4,cenind,PSmatA,nbasis,tmpvalA)
            call calcmultibndord(4,cenind,PSmatB,nbasis,tmpvalB)
            tmpval=8*(tmpvalA+tmpvalB) !8=2^(n-1)
        end if
        tmpval=tmpval/3 !Convert 4c-MCI to 4c-ESI according to Eq.10 of AV1245 paper
        write(*,"(' 4-center electron sharing index of',4i6,':',f14.8)") cenind(1:4),tmpval
        if (abs(tmpval)<AVmin) then
            AVmin=abs(tmpval)
            minidx(:)=cenind(1:4)
        end if
        totval=totval+tmpval
        if (ipos==natm) exit
        ipos=ipos+1
    end do
    
    totval=totval/natm
    write(*,"(/,a,f14.8)") " AV1245 times 1000 for the selected atoms is",totval*1000
    write(*,"(a,f12.6,' (',4i5,')')") " AVmin times 1000 for the selected atoms is ",AVmin*1000,minidx(:)
    !write(*,"(a,f14.8)") " AV1245 times 1000 for the selected atoms is",totval*1000*0.635 !mimic data of AV1245 paper
    deallocate(atmarr)
end do

iMCBOtype=iMCBOtype_old
end subroutine




!---------------------------------------
!----------- Multi-center index --------
!---------------------------------------
!ientry denotes how this module is entered. =2 Normal MCBO =-3 MCBO under Lowdin symmetrized basis
subroutine multicenter(ientry)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: PSmat(:,:),PSmata(:,:),PSmatb(:,:),PSmattot(:,:)
real*8 maxbndord
integer ientry
integer cenind(2000) !Can maximally deal with 2000 centers
integer maxcenind(2000) !Used to store the combination that has the maximum bond order during automatic search
character c2000tmp*2000

allocate(PSmat(nbasis,nbasis),PSmattot(nbasis,nbasis))
ntime=1 !Closed-shell
if (ientry==2) then
    write(*,*) "Calculating PS matrix, please wait..."
    PSmattot=matmul_blas(Ptot,Sbas,nbasis,nbasis,0,0)
else
    PSmattot=Ptot
end if
if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open-shell
    allocate(PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis))
	ntime=3
    if (ientry==2) then
        PSmata=matmul_blas(Palpha,Sbas,nbasis,nbasis,0,0)
        PSmatb=matmul_blas(Pbeta,Sbas,nbasis,nbasis,0,0)
    else
        PSmata=Palpha
        PSmatb=Pbeta
    end if
end if

if (iMCBOtype==0) then
    if (ientry==2) write(*,"(a)") " Note: Because the PS matrix may be not symmetric, the result may be dependent of &
    inputting order of atomic indices. Settings the iMCBOtype in settings.ini to 1 is recommended, in this case the results corresponding &
    to forward and reverse inputting orders will be automatically averaged"
else if (iMCBOtype==1) then
    write(*,"(a)") " Note: Since iMCBOtype in settings.ini has been set to 1, the multi-center bond order will be &
    reported after averaging the results corresponding to forward and reverse inputting directions"
else if (iMCBOtype==2) then
    write(*,"(a)") " Note: Since iMCBOtype in settings.ini has been set to 2, the multi-center bond order will be &
    reported by taking all possible permutations into account"
end if

do while(.true.)
	write(*,*)
	write(*,*) "Input atom indices, e.g. 3,4,7,8,10   (number of centers is arbitrary)"
    write(*,"(a)") " Note: The input order must be in consistency with atomic connectivity. You can also input e.g. 4-10 if the indices are contiguous"
	write(*,*) "Input -3/-4/-5/-6 will search all possible three/four/five/six-center bonds"
	write(*,*) "Input 0 can return to upper level menu"
	read(*,"(a)") c2000tmp

	if (c2000tmp(1:1)=='0') then
		Return
	else if (c2000tmp(1:1)/='-') then
		call str2arr(c2000tmp,nbndcen,cenind)
		
		do itime=1,ntime
			if (wfntype==0.or.wfntype==3) then
				PSmat=PSmattot
			else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
				if (itime==1) PSmat=PSmattot
				if (itime==2) PSmat=PSmata
				if (itime==3) PSmat=PSmatb
			end if
			call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum) !accum is pristine result without any factor
			if (itime==1) bndordmix=accum
			if (itime==2) bndordalpha=accum*2**(nbndcen-1)
			if (itime==3) bndordbeta=accum*2**(nbndcen-1)
		end do
		if (wfntype==0.or.wfntype==3) then
            if (accum>1D-8) then
			    write(*,"(a,f16.10)") " The multicenter bond order:",accum
            else
			    write(*,"(a,1PE20.10)") " The multicenter bond order:",accum
            end if
			!Normalized multicenter bond order, see Electronic Aromaticity Index for Large Rings DOI: 10.1039/C6CP00636A
			!When it is negative, first obtain **(1/n) using its absolute value, then multiply it by -1
			write(*,"(a,f16.10)") " The normalized multicenter bond order:",accum/abs(accum) * (abs(accum)**(1D0/nbndcen))
			
		else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
			totbndorder=bndordalpha+bndordbeta
            if (totbndorder>1D-5) then
			    write(*,"(a,f13.7)") " The multicenter bond order from alpha density matrix:",bndordalpha
			    write(*,"(a,f13.7)") " The multicenter bond order from beta density matrix: ",bndordbeta
			    write(*,"(a,f13.7)") " The sum of multicenter bond order from alpha and beta parts:    ",totbndorder
            else
			    write(*,"(a,1PE20.10)") " The multicenter bond order from alpha density matrix:",bndordalpha
			    write(*,"(a,1PE20.10)") " The multicenter bond order from beta density matrix: ",bndordbeta
			    write(*,"(a,1PE20.10)") " The sum of multicenter bond order from alpha and beta parts:    ",totbndorder
            end if
			write(*,"(a,f13.7)") " Above result in normalized form:",totbndorder/abs(totbndorder) * (abs(totbndorder)**(1D0/nbndcen))
            if (bndordmix>1D-5) then
			    write(*,"(a,f13.7)") " The multicenter bond order from mixed alpha&beta density matrix:",bndordmix
            else
			    write(*,"(a,1PE13.6)") " The multicenter bond order from mixed alpha&beta density matrix:",bndordmix
            end if
			write(*,"(a,f13.7)") " Above result in normalized form:",bndordmix/abs(bndordmix) * (abs(bndordmix)**(1D0/nbndcen))
		end if
		
	else if (c2000tmp(1:1)=='-') then !Automatic search
		read(c2000tmp,*) nbndcen
		nbndcen=abs(nbndcen)
		PSmat=PSmattot
		!Search all combinations. Owing to simplicity and efficiency consideration, for open-shell system, compulsory to use mixed alpha&beta density matrix
		if (wfntype==1.or.wfntype==2.or.wfntype==4) write(*,"(a)") "Note: The bond order considered here comes from mixed alpha&beta density matrix"
		write(*,*)
		write(*,*) "Input magnitude threshold for printing bond orders, e.g. 0.03"
		read(*,*) thres
		
		nfound=0
		maxbndord=0D0
		if (nbndcen/=3) write(*,*) "Note: The search may be not exhaustive. Please wait..."
		if (nbndcen==3) then
			!All combinations
			do iatm=1,ncenter
				do jatm=iatm+1,ncenter
					do katm=jatm+1,ncenter
						cenind(1)=iatm
						cenind(2)=jatm
						cenind(3)=katm
						!Clockwise and anticlockwise
						do iseq=1,2
							if (iseq==2) call invarr(cenind(1:nbndcen))
							call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
							if (abs(accum)>=thres) then
								tmp=accum/abs(accum) * (abs(accum)**(1D0/nbndcen))
								write(*,"(3i6,'  Result:',f10.6,'  Normalized:',f10.5)") cenind(1:nbndcen),accum,tmp
								nfound=nfound+1
							end if
							if (abs(accum)>maxbndord) then
								maxbndord=accum
								maxcenind=cenind
							end if
						end do
					end do
				end do
			end do
		else if (nbndcen==4) then
			!$OMP PARALLEL DO private(iatm,jatm,katm,latm,cenind,iseq,accum,tmp) shared(nfound) schedule(dynamic) NUM_THREADS(nthreads)
			do iatm=1,ncenter
				do jatm=iatm+1,ncenter
					do katm=jatm+1,ncenter
						do latm=katm+1,ncenter
							cenind(1)=iatm
							cenind(2)=jatm
							cenind(3)=katm
							cenind(4)=latm
							do iseq=1,2
								if (iseq==2) call invarr(cenind(1:nbndcen))
								call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
								if (abs(accum)>=thres) then
									tmp=accum/abs(accum) * (abs(accum)**(1D0/nbndcen))
									write(*,"(4i6,'  Result:',f10.6,'  Normalized:',f10.5)") cenind(1:nbndcen),accum,tmp
									nfound=nfound+1
								end if
								if (abs(accum)>maxbndord) then
									maxbndord=accum
									maxcenind=cenind
								end if
							end do
						end do
					end do
				end do
			end do
			!$OMP end parallel do
		else if (nbndcen==5) then
		 	do iatm=1,ncenter
				!$OMP PARALLEL DO private(jatm,katm,latm,matm,cenind,iseq,accum,tmp) shared(nfound) schedule(dynamic) NUM_THREADS(nthreads)
				do jatm=iatm+1,ncenter
					do katm=jatm+1,ncenter
						do latm=katm+1,ncenter
							do matm=latm+1,ncenter
								cenind(1)=iatm
								cenind(2)=jatm
								cenind(3)=katm
								cenind(4)=latm
								cenind(5)=matm
								do iseq=1,2
									if (iseq==2) call invarr(cenind(1:nbndcen))
									call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
									if (abs(accum)>=thres) then
										tmp=accum/abs(accum) * (abs(accum)**(1D0/nbndcen))
										write(*,"(5i6,'  Result:',f10.6,'  Normalized:',f10.5)") cenind(1:nbndcen),accum,tmp
										nfound=nfound+1
									end if
									if (abs(accum)>maxbndord) then
										maxbndord=accum
										maxcenind=cenind
									end if
								end do
							end do
						end do
					end do
				end do
				!$OMP end parallel do
			end do
		else if (nbndcen==6) then
			do iatm=1,ncenter
				!$OMP PARALLEL DO private(jatm,katm,latm,matm,cenind,iseq,accum,tmp) shared(nfound) schedule(dynamic) NUM_THREADS(nthreads)
				do jatm=iatm+1,ncenter
					do katm=jatm+1,ncenter
						do latm=katm+1,ncenter
							do matm=latm+1,ncenter
								do natm=matm+1,ncenter
									cenind(1)=iatm
									cenind(2)=jatm
									cenind(3)=katm
									cenind(4)=latm
									cenind(5)=matm
									cenind(6)=natm
									do iseq=1,2
										if (iseq==2) call invarr(cenind(1:nbndcen))
										call calcmultibndord(nbndcen,cenind,PSmat,nbasis,accum)
										if (abs(accum)>=thres) then
											tmp=accum/abs(accum) * (abs(accum)**(1D0/nbndcen))
											write(*,"(6i6,'  Result:',f10.6,'  Normalized:',f10.5)") cenind(1:nbndcen),accum,tmp
											nfound=nfound+1
										end if
										if (abs(accum)>maxbndord) then
											maxbndord=accum
											maxcenind=cenind
										end if
									end do
								end do
							end do
						end do
					end do
				end do
				!$OMP end parallel do
			end do
		end if
		if (nfound==0) then
			write(*,*) "No multi-center bonds above criteria were found"
			cycle
		end if
		write(*,*)
		write(*,*) "The maximum bond order is"
		tmp=maxbndord/abs(maxbndord) * (abs(maxbndord)**(1D0/nbndcen))
		if (nbndcen==3) write(*,"(3i6,'  Result:',f10.6,'  Normalized:',f10.5)") maxcenind(1:nbndcen),maxbndord,tmp
		if (nbndcen==4) write(*,"(4i6,'  Result:',f10.6,'  Normalized:',f10.5)") maxcenind(1:nbndcen),maxbndord,tmp
		if (nbndcen==5) write(*,"(5i6,'  Result:',f10.6,'  Normalized:',f10.5)") maxcenind(1:nbndcen),maxbndord,tmp
		if (nbndcen==6) write(*,"(6i6,'  Result:',f10.6,'  Normalized:',f10.5)") maxcenind(1:nbndcen),maxbndord,tmp
	end if
end do
end subroutine




!-------------------------------------------------
!-------- Multi-center index in NAO basis --------
!-------------------------------------------------
subroutine multicenterNAO
use defvar
use util
use NAOmod
implicit real*8 (a-h,o-z)
integer cenind(2000) !Can maximally deal with 2000 centers
character :: c2000tmp*2000

!Load NAO and DMNAO information
open(10,file=filename,status="old")
call checkNPA(ifound);if (ifound==0) return
call loadNAOinfo
call checkDMNAO(ifound);if (ifound==0) return
call loadDMNAO
close(10)

!Move information from NAO variables to common variables, so that multi-center bond order routines could be used
ncenter=ncenter_NAO
if (allocated(basstart)) deallocate(basstart,basend)
allocate(basstart(ncenter),basend(ncenter))
basstart=NAOinit
basend=NAOend

if (iMCBOtype==2) then
    write(*,"(a)") " Note: Since iMCBOtype in settings.ini has been set to 2, the multi-center bond order will be &
    reported by taking all possible permutations into account"
end if

do while(.true.)
	write(*,*)
	write(*,*) "Input atom indices, e.g. 3,4,7,8,10   (number of centers is arbitrary)"
    write(*,"(a)") " Note: The input order must be in consistency with atomic connectivity. You can also input e.g. 4-10 if the indices are contiguous"
	write(*,*) "Input 0 can exit"
	read(*,"(a)") c2000tmp
	if (c2000tmp(1:1)=='0') then
		deallocate(basstart,basend)
		return
	else
		call str2arr(c2000tmp,nbndcen,cenind)

		if (.not.allocated(DMNAOb)) then !Closed shell
			call calcmultibndord(nbndcen,cenind,DMNAO,numNAO,bndord)
			if (nbndcen==2) then
				write(*,"(a,f16.10)") " The Wiberg bond order:",bndord
			else
				write(*,"(a,f16.10)") " The multicenter bond order:",bndord
				write(*,"(a,f16.10)") " The normalized multicenter bond order:",bndord/abs(bndord) * (abs(bndord)**(1D0/nbndcen))
			end if
		else !Open shell
            write(*,*) "Calculating based on alpha density matrix..."
			call calcmultibndord(nbndcen,cenind,DMNAOa,numNAO,bndordalpha)
            write(*,*) "Calculating based on beta density matrix..."
			call calcmultibndord(nbndcen,cenind,DMNAOb,numNAO,bndordbeta)
            write(*,*) "Calculating based on mixed alpha&beta density matrix..."
			call calcmultibndord(nbndcen,cenind,DMNAO,numNAO,bndordmix)
			bndordalpha=bndordalpha*2**(nbndcen-1)
			bndordbeta=bndordbeta*2**(nbndcen-1)
			totbndorder=bndordalpha+bndordbeta
			if (nbndcen==2) then
				write(*,"(a,f16.10)") " The bond order from alpha density matrix:",bndordalpha
				write(*,"(a,f16.10)") " The bond order from beta density matrix: ",bndordbeta
				write(*,"(a,f16.10)") " The sum of above two terms:",bndordalpha+bndordbeta
				write(*,"(a,f16.10)") " The bond order from mixed alpha&beta density matrix: ",bndordmix
			else
				write(*,"(a,f13.7)") " The multicenter bond order from alpha density matrix:",bndordalpha
				write(*,"(a,f13.7)") " The multicenter bond order from beta density matrix: ",bndordbeta
				write(*,"(a,f13.7)") " The sum of multicenter bond order from alpha and beta parts:    ",totbndorder
				write(*,"(a,f13.7)") " Above result in normalized form:",totbndorder/abs(totbndorder) * (abs(totbndorder)**(1D0/nbndcen))
				write(*,"(a,f13.7)") " The multicenter bond order from mixed alpha&beta density matrix:",bndordmix
				write(*,"(a,f13.7)") " Above result in normalized form:",bndordmix/abs(bndordmix) * (abs(bndordmix)**(1D0/nbndcen))
			end if
		end if
	end if
end do
end subroutine



!------ General routine directly calculates two- or multi-center bond order without complex things
!This is a wrapper of subroutine "calcmultibndord_do" for returning different definitions of multi-center bond order
!Shared by subroutine multicenter, multicenterNAO, AV1245 and others
!MCBOtype=0, return the MCBO in usual manner
!MCBOtype=1, return the averaged result of positive order and reverse order of inputted atoms
!MCBOtype=2, return the MCBO calculated as Eq. 9 of the AV1245 paper, namely taking all permutation into account
!Note that for open-shell cases, the returned result should then be multiplied by a proper factor
!  Input variables:
!nbndcen: Actual number of atoms to be calculated
!PSmat: Commonly constructed as e.g. matmul(Ptot,Sbas)
!cenind: Atomic indices, must be size of 2000
!matdim: Commonly is nbasis
subroutine calcmultibndord(nbndcen,cenind,PSmat,matdim,result)
use defvar
use util
implicit real*8 (a-h,o-z)
integer nbndcen,cenind(2000),cenindtmp(2000),matdim
real*8 PSmat(matdim,matdim),result
integer,allocatable :: allperm(:,:)

if (iMCBOtype==0) then
    call calcmultibndord_do(nbndcen,cenind,PSmat,matdim,result)
else if (iMCBOtype==1) then
    call calcmultibndord_do(nbndcen,cenind,PSmat,matdim,result1)
    do i=1,nbndcen !Reverse order
        cenindtmp(nbndcen-i+1)=cenind(i)
    end do
    call calcmultibndord_do(nbndcen,cenindtmp,PSmat,matdim,result2)
    result=(result1+result2)/2
else if (iMCBOtype==2) then
    nperm=ft(nbndcen)
    allocate(allperm(nperm,nbndcen))
    call fullarrange(allperm,nperm,nbndcen) !Generate all possible permutation sequence
    result=0
    !$OMP PARALLEL DO SHARED(result) PRIVATE(iperm,cenindtmp,resulttmp) schedule(dynamic) NUM_THREADS(nthreads)
    do iperm=1,nperm
        cenindtmp(1:nbndcen)=cenind(allperm(iperm,:))
        call calcmultibndord_do(nbndcen,cenindtmp,PSmat,matdim,resulttmp)
        !$OMP CRITICAL
        result=result+resulttmp
        !$OMP END CRITICAL
    end do
    !$OMP END PARALLEL DO
    result=result/(2*nbndcen)
end if
end subroutine



!!------- Actual working horse for multi-center index calculation
subroutine calcmultibndord_do(nbndcen,cenind,PSmat,matdim,result)
use defvar
implicit real*8 (a-h,o-z)
integer nbndcen,cenind(2000),matdim !Can maximally deal with 2000 centers
real*8 PSmat(matdim,matdim),result
real*8,allocatable :: mat1(:,:),mat2(:,:)

result=0
if (nbndcen==2) then !Special case, only two atoms
	do ib=basstart(cenind(2)),basend(cenind(2))
		do ia=basstart(cenind(1)),basend(cenind(1))
	        result=result+PSmat(ia,ib)*PSmat(ib,ia)
		end do
	end do
    return
end if

allocate(mat1(matdim,matdim),mat2(matdim,matdim))

!jatm is the atom index to be looped and summed up in each time of contraction
!mat1 and mat2 are two matrices used alternately to perform tensor contraction
iatm=nbndcen-1
katm=1
kbeg=basstart(cenind(katm))
kend=basend(cenind(katm))
do icontract=1,nbndcen-2
    !write(*,"(' Doing contraction',i3)") icontract
    ibeg=basstart(cenind(iatm))
    iend=basend(cenind(iatm))
    jatm=iatm+1
    jbeg=basstart(cenind(jatm))
    jend=basend(cenind(jatm))
    if (icontract==1) then !Initialize
        icalc=1
        mat2=PSmat
    end if
    
    if (icalc==1) then
        do ibas=ibeg,iend
            do kbas=kbeg,kend
                mat1(ibas,kbas)=sum(PSmat(ibas,jbeg:jend)*mat2(jbeg:jend,kbas))
            end do
        end do
        icalc=2
    else if (icalc==2) then
        do ibas=ibeg,iend
            do kbas=kbeg,kend
                mat2(ibas,kbas)=sum(PSmat(ibas,jbeg:jend)*mat1(jbeg:jend,kbas))
            end do
        end do
        icalc=1
    end if
    iatm=iatm-1
end do

result=0
do ibas=basstart(cenind(1)),basend(cenind(1))
    jbeg=basstart(cenind(2))
    jend=basend(cenind(2))
    if (icalc==1) then
        result=result+sum(PSmat(ibas,jbeg:jend)*mat2(jbeg:jend,ibas))
    else if (icalc==2) then
        result=result+sum(PSmat(ibas,jbeg:jend)*mat1(jbeg:jend,ibas))
    end if
end do
end subroutine

!Extremely SLOW, foolish, lengthy code! In addition, compiling this code is very slow under ifort with O1 or O2
!!---- The actual working horse for multi-center bond order calculation
!subroutine calcmultibndord_do(nbndcen,cenind,PSmat,matdim,result)
!use defvar
!implicit real*8 (a-h,o-z)
!integer nbndcen,cenind(12),matdim
!real*8 PSmat(matdim,matdim),result
!
!result=0D0
!if (nbndcen==2) then
!	do ib=basstart(cenind(2)),basend(cenind(2))
!		do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ia)
!		end do
!	end do
!else if (nbndcen==3) then
!	do ic=basstart(cenind(3)),basend(cenind(3))
!		do ib=basstart(cenind(2)),basend(cenind(2))
!			do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,ia)
!			end do
!		end do
!	end do
!else if (nbndcen==4) then
!	do id=basstart(cenind(4)),basend(cenind(4))
!		do ic=basstart(cenind(3)),basend(cenind(3))
!			do ib=basstart(cenind(2)),basend(cenind(2))
!				do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ia)
!				end do
!			end do
!		end do
!	end do
!else if (nbndcen==5) then
!	do ie=basstart(cenind(5)),basend(cenind(5))
!		do id=basstart(cenind(4)),basend(cenind(4))
!			do ic=basstart(cenind(3)),basend(cenind(3))
!				do ib=basstart(cenind(2)),basend(cenind(2))
!					do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,ia)
!					end do
!				end do
!			end do
!		end do
!	end do
!else if (nbndcen==6) then
!	do i_f=basstart(cenind(6)),basend(cenind(6))
!		do ie=basstart(cenind(5)),basend(cenind(5))
!			do id=basstart(cenind(4)),basend(cenind(4))
!				do ic=basstart(cenind(3)),basend(cenind(3))
!					do ib=basstart(cenind(2)),basend(cenind(2))
!						do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,ia)
!						end do
!					end do
!				end do
!			end do
!		end do
!	end do
!else if (nbndcen==7) then
!	itmp=0
!	ntot=basend(cenind(7))-basstart(cenind(7))+1
!	do i_g=basstart(cenind(7)),basend(cenind(7))
!		call showprog(itmp,ntot)
!		do i_f=basstart(cenind(6)),basend(cenind(6))
!			do ie=basstart(cenind(5)),basend(cenind(5))
!				do id=basstart(cenind(4)),basend(cenind(4))
!					do ic=basstart(cenind(3)),basend(cenind(3))
!						do ib=basstart(cenind(2)),basend(cenind(2))
!							do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,ia)
!							end do
!						end do
!					end do
!				end do
!			end do
!		end do
!		itmp=itmp+1
!	end do
!else if (nbndcen==8) then
!	itmp=0
!	ntot=basend(cenind(8))-basstart(cenind(8))+1
!	do i_h=basstart(cenind(8)),basend(cenind(8))
!		call showprog(itmp,ntot)
!		do i_g=basstart(cenind(7)),basend(cenind(7))
!			do i_f=basstart(cenind(6)),basend(cenind(6))
!				do ie=basstart(cenind(5)),basend(cenind(5))
!					do id=basstart(cenind(4)),basend(cenind(4))
!						do ic=basstart(cenind(3)),basend(cenind(3))
!							do ib=basstart(cenind(2)),basend(cenind(2))
!								do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,ia)
!								end do
!							end do
!						end do
!					end do
!				end do
!			end do
!		end do
!		itmp=itmp+1
!	end do
!else if (nbndcen==9) then
!	itmp=0
!	ntot=basend(cenind(9))-basstart(cenind(9))+1
!	do i_i=basstart(cenind(9)),basend(cenind(9))
!		call showprog(itmp,ntot)
!		do i_h=basstart(cenind(8)),basend(cenind(8))
!			do i_g=basstart(cenind(7)),basend(cenind(7))
!				do i_f=basstart(cenind(6)),basend(cenind(6))
!					do ie=basstart(cenind(5)),basend(cenind(5))
!						do id=basstart(cenind(4)),basend(cenind(4))
!							do ic=basstart(cenind(3)),basend(cenind(3))
!								do ib=basstart(cenind(2)),basend(cenind(2))
!									do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,ia)
!									end do
!								end do
!							end do
!						end do
!					end do
!				end do
!			end do
!		end do
!		itmp=itmp+1
!	end do
!else if (nbndcen==10) then
!	itmp=0
!	ntot=basend(cenind(10))-basstart(cenind(10))+1
!	do i_j=basstart(cenind(10)),basend(cenind(10))
!		call showprog(itmp,ntot)
!		do i_i=basstart(cenind(9)),basend(cenind(9))
!			do i_h=basstart(cenind(8)),basend(cenind(8))
!				do i_g=basstart(cenind(7)),basend(cenind(7))
!					do i_f=basstart(cenind(6)),basend(cenind(6))
!						do ie=basstart(cenind(5)),basend(cenind(5))
!							do id=basstart(cenind(4)),basend(cenind(4))
!								do ic=basstart(cenind(3)),basend(cenind(3))
!									do ib=basstart(cenind(2)),basend(cenind(2))
!										do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,i_j)*PSmat(i_j,ia)
!										end do
!									end do
!								end do
!							end do
!						end do
!					end do
!				end do
!			end do
!		end do
!		itmp=itmp+1
!	end do
!else if (nbndcen==11) then
!	itmp=0
!	ntot=( basend(cenind(11))-basstart(cenind(11))+1 ) * ( basend(cenind(10))-basstart(cenind(10))+1 )
!	do i_k=basstart(cenind(11)),basend(cenind(11))
!		do i_j=basstart(cenind(10)),basend(cenind(10))
!			call showprog(itmp,ntot)
!			do i_i=basstart(cenind(9)),basend(cenind(9))
!				do i_h=basstart(cenind(8)),basend(cenind(8))
!					do i_g=basstart(cenind(7)),basend(cenind(7))
!						do i_f=basstart(cenind(6)),basend(cenind(6))
!							do ie=basstart(cenind(5)),basend(cenind(5))
!								do id=basstart(cenind(4)),basend(cenind(4))
!									do ic=basstart(cenind(3)),basend(cenind(3))
!										do ib=basstart(cenind(2)),basend(cenind(2))
!											do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,i_j)*PSmat(i_j,i_k)*PSmat(i_k,ia)
!											end do
!										end do
!									end do
!								end do
!							end do
!						end do
!					end do
!				end do
!			end do
!			itmp=itmp+1
!		end do
!	end do
!else if (nbndcen==12) then
!	itmp=0
!	ntot=( basend(cenind(12))-basstart(cenind(12))+1 ) * ( basend(cenind(11))-basstart(cenind(11))+1 ) * ( basend(cenind(10))-basstart(cenind(10))+1 )
!	do i_l=basstart(cenind(12)),basend(cenind(12))
!		do i_k=basstart(cenind(11)),basend(cenind(11))
!			do i_j=basstart(cenind(10)),basend(cenind(10))
!				call showprog(itmp,ntot)
!				do i_i=basstart(cenind(9)),basend(cenind(9))
!					do i_h=basstart(cenind(8)),basend(cenind(8))
!						do i_g=basstart(cenind(7)),basend(cenind(7))
!							do i_f=basstart(cenind(6)),basend(cenind(6))
!								do ie=basstart(cenind(5)),basend(cenind(5))
!									do id=basstart(cenind(4)),basend(cenind(4))
!										do ic=basstart(cenind(3)),basend(cenind(3))
!											do ib=basstart(cenind(2)),basend(cenind(2))
!												do ia=basstart(cenind(1)),basend(cenind(1))
!	result=result+PSmat(ia,ib)*PSmat(ib,ic)*PSmat(ic,id)*PSmat(id,ie)*PSmat(ie,i_f)*PSmat(i_f,i_g)*PSmat(i_g,i_h)*PSmat(i_h,i_i)*PSmat(i_i,i_j)*PSmat(i_j,i_k)*PSmat(i_k,i_l)*PSmat(i_l,ia)
!												end do
!											end do
!										end do
!									end do
!								end do
!							end do
!						end do
!					end do
!				end do
!				itmp=itmp+1
!			end do
!		end do
!	end do
!end if
!if (nbndcen>=7) then
!    call showprog(ntot,ntot)
!    write(*,*)
!end if
!end subroutine




!------------------------------------------------------------------------------------------------
!--------- Generate Bq position for calculating NICS_ZZ, and obtain NICS_ZZ for non-planar system
!------------------------------------------------------------------------------------------------
subroutine NICS_ZZ
use util
use defvar
implicit real*8 (a-h,o-z)
character c1000*1000,c2000tmp*2000,c200*200
integer,allocatable :: tmparr(:)
real*8 hesstmp(3,3)
write(*,*) "Input center coordinate of the ring (in Angstrom), e.g. 2.0,2.4,1.1"
write(*,*) "Note 1: If use Bohr as unit, the first letter should be ""b"", e.g. b3.0,3.8,2.2"
write(*,"(a)") " Note 2: If press ENTER button directly, the geometric center of the atoms inputted in the next stage will be employed as the ring center"
read(*,"(a)") c1000
if (c1000(1:1)=='b') then
	read(c1000(2:),*) tmpx,tmpy,tmpz
else if (c1000==" ") then
    continue
else
	read(c1000,*) tmpx,tmpy,tmpz
	tmpx=tmpx/b2a
	tmpy=tmpy/b2a
	tmpz=tmpz/b2a
end if

write(*,*) "Input indices of at least three atoms in the ring to fit a plane"
write(*,*) "For example: 2,3,7-10,15"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,ntmp)
allocate(tmparr(ntmp))
call str2arr(c2000tmp,ntmp,tmparr)
call ptsfitplane(tmparr,ntmp,xnor,ynor,znor,rnouse,rmsfit)
write(*,"(' RMS error of the plane fitting:',f12.6,' Angstrom')") rmsfit*b2a
facnorm=sqrt(xnor**2+ynor**2+znor**2)
xnor=xnor/facnorm !Normalize normal vector, then (xnor,ynor,znor) is the unit vector normal to the fitted plane
ynor=ynor/facnorm
znor=znor/facnorm
write(*,"(' The unit normal vector is',3f14.8)") xnor,ynor,znor

if (c1000==" ") then
    tmpx=sum(a(tmparr)%x)/ntmp
    tmpy=sum(a(tmparr)%y)/ntmp
    tmpz=sum(a(tmparr)%z)/ntmp
    write(*,"(' The ring center is',3f12.5,' Angstrom')") tmpx*b2a,tmpy*b2a,tmpz*b2a
end if
write(*,*)
write(*,"(a)") " The X,Y,Z coordinate of the points below and above 1 Angstrom of the plane from the point you defined, respectively:"
write(*,"(3f16.10,' Angstrom')") (tmpx-xnor/b2a)*b2a,(tmpy-ynor/b2a)*b2a,(tmpz-znor/b2a)*b2a
write(*,"(3f16.10,' Angstrom')") (tmpx+xnor/b2a)*b2a,(tmpy+ynor/b2a)*b2a,(tmpz+znor/b2a)*b2a
write(*,*)
write(*,"(a)") " Now input magnetic shielding tensor outputted by your ab-initio program, then Multiwfn will calculate the shielding value in the direction perpendicular to the plane. &
You can also input ""q"" to return"
! If you are a Gaussian user, you can also directly copy NMR shielding tensor (such as below) from Gaussian output file to present window
! (Note that this is not always a symmetric matrix)
!    XX=     6.2246   YX=   -54.8190   ZX=    94.7322
!    XY=    -3.1233   YY=    51.8304   ZY=     4.0216
!    XZ=    96.6247   YZ=    38.6474   ZZ=    94.4768
write(*,*) "Input XX,YX,ZX component of the tensor, e.g. 12.5150,-2.3289,3.7151"
read(*,"(a)") c1000
if (index(c1000,'q')/=0) then
	return
else if (index(c1000,"XX")==0) then
	read(c1000,*) hesstmp(1,1),hesstmp(1,2),hesstmp(1,3)
	write(*,*) "Input XY,YY,ZY of the tensor, e.g. -2.3289,10.0169,-2.2232"
	read(*,*) hesstmp(2,1),hesstmp(2,2),hesstmp(2,3)
	write(*,*) "Input XZ,YZ,ZZ of the tensor, e.g. 3.7151,-2.2232,12.1699"
	read(*,*) hesstmp(3,1),hesstmp(3,2),hesstmp(3,3)
else
	read(c1000,*) c200,hesstmp(1,1),c200,hesstmp(1,2),c200,hesstmp(1,3)
	read(*,"(a)") c1000	
	read(c1000,*) c200,hesstmp(2,1),c200,hesstmp(2,2),c200,hesstmp(2,3)
	read(*,"(a)") c1000			
	read(c1000,*) c200,hesstmp(3,1),c200,hesstmp(3,2),c200,hesstmp(3,3)
end if
write(*,*) "The magnetic shielding tensor you inputted is:"
write(*,"(' XX=',f12.6,'   YX=',f12.6,'   ZX=',f12.6)") hesstmp(1,:)
write(*,"(' XY=',f12.6,'   YY=',f12.6,'   ZY=',f12.6)") hesstmp(2,:)
write(*,"(' XZ=',f12.6,'   YZ=',f12.6,'   ZZ=',f12.6)") hesstmp(3,:)
shieldperpen=xnor*xnor*hesstmp(1,1)+xnor*ynor*hesstmp(1,2)+xnor*znor*hesstmp(1,3)+&
			 ynor*xnor*hesstmp(2,1)+ynor*ynor*hesstmp(2,2)+ynor*znor*hesstmp(2,3)+&
			 znor*xnor*hesstmp(3,1)+znor*ynor*hesstmp(3,2)+znor*znor*hesstmp(3,3)
write(*,"(' The shielding value normal to the plane is',f20.10)") shieldperpen
write(*,"(' The NICS_ZZ value is thus',f20.10,/)") -shieldperpen
end subroutine