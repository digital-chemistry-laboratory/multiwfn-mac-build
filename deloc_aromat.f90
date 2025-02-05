!---------------------------------------------------------------------
!-------- Interface of various delocalization and aromaticity analyses
!---------------------------------------------------------------------
subroutine deloc_aromat
implicit real*8 (a-h,o-z)
character c80tmp*80
do while(.true.)
	write(*,*)
	write(*,*) "   ================ Delocalization and aromaticity analyses ==============="
    write(*,*) " 0 Return to main menu"
    write(*,*) " 1 Multicenter bond order"
    write(*,*) "-1 Multicenter bond order in NAO basis"
    write(*,*) " 2 AV1245 index"
    write(*,*) " 3 Iso-chemical shielding surface (ICSS)"
    write(*,*) " 4 NICS_ZZ for non-planar or tilted system"
    write(*,*) " 5 ELF-sigma/pi and LOL-sigma/pi"
    write(*,*) " 6 Harmonic oscillator measure of aromaticity (HOMA) and Bird indices"
    write(*,*) " 6a HOMAc (reparameterized HOMA)   6b HOMER (HOMA for excited states)"
    write(*,*) " 7 Shannon aromaticity index"
    write(*,*) " 8 Para-delocalization index (PDI)"
    write(*,*) " 9 Aromatic fluctuation index (FLU) and FLU-pi"
    write(*,*) "10 Para linear response index (PLR)"
    write(*,*) "11 Information-theoretic (ITA) aromaticity index"
    write(*,*) "12 Properties of ring critical point"
    write(*,*) "13 NICS-1D scan curve map, integral NICS (INICS) and FiPC-NICS"
    write(*,*) "14 NICS-2D scan plane map"
    read(*,"(a)") c80tmp
    if (c80tmp=="6a") then
		call HOMA_reparm(1)
    else if (c80tmp=="6b") then
		call HOMA_reparm(2)
    else
		read(c80tmp,*) isel
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
			&performing topology analysis to obtain bifurcation value, etc., which rely on different main functions. Please check &
			&Section 4.4.9, 4.5.3 and 4.100.22 of manunal for example of realizing these analyses."
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		else if (isel==6) then
			call HOMA_Bird
		else if (isel==7) then
			write(*,"(a)") " To realize this analysis, you should use topology analysis module (main function 2), see Section 3.14.6 of manual &
			&for introduction and Section 4.2.1 for practical example."
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		else if (isel==8) then
			write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.6 of manual &
			&for introduction and Section 4.15.2 for practical example."
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		else if (isel==9) then
			write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.7 of manual &
			&for introduction and Section 4.15.2 for practical example."
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		else if (isel==10) then
			write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.9 of manual &
			&for introduction and Section 4.15.2 for practical example."
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		else if (isel==11) then
			write(*,"(a)") " To realize this analysis, you should use fuzzy analysis module (main function 15), see Section 3.18.11 of manual for detail."
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		else if (isel==12) then
			write(*,"(a)") " To realize this analysis, you should use topology analysis module (main function 2), see Section 3.14.6 of manual &
			&for introduction and Section 4.2.1 for practical example."
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		else if (isel==13) then
			call NICS_1D
		else if (isel==14) then
			call study2dim(1,0,0)
		end if
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
			if (index(gauinpcontent(i),'nmr')==0.and.index(gauinpcontent(i),'NMR')==0) gauinpcontent(i)=trim(gauinpcontent(i))//" NMR"
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
&C:\ltwd\NICS0003.out... will be loaded (.log suffix is also allowed)"
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
	call menutitle("HOMA / Bird aromaticity index",10,1)
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
		write(*,"(a)") " Current reference bond length (Angstrom) and sigma (Angstrom^-2) parameters:"
		do iref=1,nelesupp
			do jref=iref,nelesupp
				if (HOMArefbond(iref,jref)/=-1) write(*,"(' ',a,a,a,a,2f12.4)") ind2name(iref),'-',ind2name(jref),':',HOMArefbond(iref,jref),HOMAsigma(iref,jref)
			end do
		end do
		write(*,*)
		do while(.true.)
			write(*,"(a)") " Input indices of the atoms according to bonding relationship in the ring, e.g. 1,5,6,7,8,12"
			write(*,*) "(Input q can return)"
			read(*,"(a)") c200inp
			if (c200inp=='q') exit
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
				paircontri=-refsigma/numHOMAatm*(refbondlen-atomdist(iatm,jatm,1)*b2a)**2
				write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),paircontri,atomdist(iatm,jatm,1)*b2a
				HOMAval=HOMAval+paircontri
				if (iidx==numHOMAatm) write(*,"(a,f12.6)") " HOMA value is",HOMAval
			end do
			write(*,*)
		end do
	else if (isel==1) then
		write(*,"(a)") " Current reference bond length (Angstrom) and sigma (Angstrom^-2) parameters:"
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
			if (C200inp(1:1)=='q') exit
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
			write(*,"(a)") " Input indices of the atoms according to bonding relationship in the ring, e.g. 1,5,6,7,8,12"
			write(*,*) "(Input q can return)"
			read(*,"(a)") c200inp
			if (c200inp=='q') exit
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
				iatmeleidx=a(iatm)%index !Element index in periodic table
				jatmeleidx=a(jatm)%index
				Birdanow=Birda(iatmeleidx,jatmeleidx)
				Birdbnow=Birdb(iatmeleidx,jatmeleidx)
				if (Birdanow==-1D0) then
					write(*,"(' Error: Missing a and b parameters for ',a,'-',a)") ind2name(iatmeleidx),ind2name(jatmeleidx)
					exit
				end if
				BirdN(iidx)=Birdanow/(atomdist(iatm,jatm,1)*b2a)**2-Birdbnow
				write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),BirdN(iidx),atomdist(iatm,jatm,1)*b2a
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
			if (C200inp(1:1)=='q') exit
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




!------------------------------------------------------------------------
!------- Reparameterized HOMA (HOMAc and HOMER) aromaticity index -------
!------------------------------------------------------------------------
!itype=1: HOMAc   =2: HOMER
subroutine HOMA_reparm(itype)
use defvar
use util
implicit real*8 (a-h,o-z)
integer :: numHOMAatm=6,HOMAatm(1000)
real*8 :: refbond(nelesupp,nelesupp)=-1D0,sigma(nelesupp,nelesupp)=0D0,HOMAval
character C200tmp*200,label*5

!Table 1 of HOMAc original paper
if (itype==1) then !HOMAc
	refbond(6,6)=1.392D0
	refbond(6,7)=1.333D0
	refbond(7,6)=refbond(6,7)
	refbond(6,8)=1.315D0
	refbond(8,6)=refbond(6,8)
	refbond(7,7)=1.318D0
	sigma(6,6)=153.37D0
	sigma(6,7)=111.83D0
	sigma(7,6)=sigma(6,7)
	sigma(6,8)=335.16D0
	sigma(8,6)=sigma(6,8)
	sigma(7,7)=98.99D0
    label="HOMAc"
	call menutitle("HOMAc (HOMA based on computational bonds)",10,1)
else if (itype==2) then !HOMER
	refbond(6,6)=1.437D0
	refbond(6,7)=1.390D0
	refbond(7,6)=refbond(6,7)
	refbond(6,8)=1.379D0
	refbond(8,6)=refbond(6,8)
	refbond(7,7)=1.375D0
	sigma(6,6)=950.74D0
	sigma(6,7)=506.43D0
	sigma(7,6)=sigma(6,7)
	sigma(6,8)=164.96D0
	sigma(8,6)=sigma(6,8)
	sigma(7,7)=187.36D0
    label="HOMER"
	call menutitle("HOMER (Harmonic Oscillator Model of Excited-state aRomaticity)",10,1)
end if

write(*,"(a)") " Reference bond length (Angstrom) and sigma parameters (Angstrom^-2) defined by "//label//":"
do iref=1,nelesupp
	do jref=iref,nelesupp
		if (refbond(iref,jref)/=-1) write(*,"(' ',a,a,a,a,2f12.4)") ind2name(iref),'-',ind2name(jref),':',refbond(iref,jref),sigma(iref,jref)
	end do
end do
do while(.true.)
	write(*,*)
	write(*,"(a)") " Input indices of the atoms according to bonding relationship in the ring, e.g. 1,5,6,7,8,12"
	write(*,*) "(Input q can exit)"
	read(*,"(a)") C200tmp
	if (C200tmp=='q') exit
	call str2arr(C200tmp,numHOMAatm,HOMAatm)
	HOMAval=1D0
    write(*,*)
	write(*,*) "        Atom pair         Contribution  Bond length(Angstrom)"
	do iidx=1,numHOMAatm
		jidx=iidx+1
		if (iidx==numHOMAatm) jidx=1
		iatm=HOMAatm(iidx) !Actual atom index in present system
		jatm=HOMAatm(jidx)
		iatmeleidx=a(iatm)%index !Element index in periodic table
		jatmeleidx=a(jatm)%index
		refbondlen=refbond(iatmeleidx,jatmeleidx)
		refsigma=sigma(iatmeleidx,jatmeleidx)
		if (refbondlen==-1D0) then
			write(*,"(' Error: Missing reference parameter for ',a,'-',a)") ind2name(iatmeleidx),ind2name(jatmeleidx)
			exit
		end if
		paircontri=-refsigma/numHOMAatm*(refbondlen-atomdist(iatm,jatm,1)*b2a)**2
		write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),paircontri,atomdist(iatm,jatm,1)*b2a
		HOMAval=HOMAval+paircontri
		if (iidx==numHOMAatm) write(*,"(1x,a,f12.6)") label//" value is",HOMAval
	end do
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
call ask_Sbas_PBC

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
    &will be automatically guessed, however in this case any atom should not connect to more than two atoms in the ring"
    read(*,"(a)") c2000tmp
    
    if (index(c2000tmp,'q')/=0) then
        exit
    else if (index(c2000tmp,'d')/=0) then
        if (.not.allocated(connmat)) call genconnmat(1,1) !Generate connectivity matrix
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
            &the atom indices according to connectivity"
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
    &inputting order of atomic indices. Settings the iMCBOtype in settings.ini to 1 is recommended, in this case the results corresponding &
    &to forward and reverse inputting orders will be automatically averaged"
else if (iMCBOtype==1) then
    write(*,"(a)") " Note: Since iMCBOtype in settings.ini has been set to 1, the multi-center bond order will be &
    &reported after averaging the results corresponding to forward and reverse inputting directions"
else if (iMCBOtype==2) then
    write(*,"(a)") " Note: Since iMCBOtype in settings.ini has been set to 2, the multi-center bond order will be &
    &reported by taking all possible permutations into account"
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
			    write(*,"(a,1PE17.10)") " The sum of multicenter bond order from alpha and beta parts:    ",totbndorder
            end if
			write(*,"(a,f13.7)") " Above result in normalized form:",totbndorder/abs(totbndorder) * (abs(totbndorder)**(1D0/nbndcen))
            if (bndordmix>1D-5) then
			    write(*,"(a,f13.7)") " The multicenter bond order from mixed alpha&beta density matrix:",bndordmix
            else
			    write(*,"(a,1PE16.7)") " The multicenter bond order from mixed alpha&beta density matrix:",bndordmix
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
    &reported by taking all possible permutations into account"
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
&You can also input ""q"" to return"
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
!This part can also use "function prjmat" to evaluate
shieldperpen=xnor*xnor*hesstmp(1,1)+xnor*ynor*hesstmp(1,2)+xnor*znor*hesstmp(1,3)+&
			 ynor*xnor*hesstmp(2,1)+ynor*ynor*hesstmp(2,2)+ynor*znor*hesstmp(2,3)+&
			 znor*xnor*hesstmp(3,1)+znor*ynor*hesstmp(3,2)+znor*znor*hesstmp(3,3)
write(*,"(' The shielding value normal to the plane is',f20.10)") shieldperpen
write(*,"(' The NICS_ZZ value is thus',f20.10,/)") -shieldperpen
end subroutine




!-------------------------------------------------------
!--------- NICS-1D scan curve map and integral ---------
!-------------------------------------------------------
subroutine NICS_1D
use defvar
use util
use plot
implicit real*8 (a-h,o-z)
real*8 endpt1(3),endpt2(3),normalvec(3),uvec(3),cenpos(3),tmpvec(3)
character c80tmp*80,c2000tmp*2000,graphformat_old*4
integer,allocatable :: tmparr(:)
integer npt !Number of ghost atoms
real*8,allocatable :: ptxyz(:,:) !XYZ of scanning points, ptpos(3,npt)
real*8,allocatable :: ptpos(:) !Relative position of scanning points. In the case of inputting two points, relative to the first one. Case of relative to plane, relative to plane center
real*8,allocatable :: pttens(:,:,:) !pttens(1:3,1:3,ipt) is 3*3 NICS tensor of ipt (negative of magnetic shielding tensor)
real*8,allocatable :: ptNICS(:) !Specific component value of NICS of scanning points
real*8,allocatable :: ptNICSiso(:),ptNICSaniso(:) !Isotropic and anisotropic NICS of scanning points
real*8,allocatable :: compin(:),compout(:) !In-plane and out-of-plane components of NICS

write(*,*)
call menutitle("NICS-1D scan and integral",10,1)
write(*,*) "0 Return"
write(*,*) "Choose the way of defining the two end points of the line for scanning"
write(*,*) "1 Directly input Cartesian coordinates of two end points"
write(*,"(a)") " 2 The two end points are above and below the center of a plane fitted for specific atoms, and the line perpendicularly passes through their center"
read(*,*) iway

if (iway==0) then
    return
else if (iway==1) then
    write(*,*) "Input X,Y,Z of the first end point in Angstrom, e.g. 3.2,0.1,-9.5"
    read(*,*) endpt1(:)
    endpt1(:)=endpt1(:)/b2a
    write(*,*) "Input X,Y,Z of the second end point in Angstrom, e.g. 3.2,0.1,-9.5"
    read(*,*) endpt2(:)
    endpt2(:)=endpt2(:)/b2a
    tmpvec=endpt2-endpt1
    uvec=tmpvec/dsqrt(sum(tmpvec**2))
else
    write(*,*) "Input the atoms constituting the plane, e.g. 1-3,9,10"
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,ntmp)
    allocate(tmparr(ntmp))
    call str2arr(c2000tmp,ntmp,tmparr)
    normalvec=0
    call ptsfitplane(tmparr,ntmp,normalvec(1),normalvec(2),normalvec(3),rnouse,rmsfit)
    uvec=normalvec/dsqrt(sum(normalvec**2))
    write(*,"(' Unit normal vector of fitted plane:',3f10.5)") uvec(:)
    write(*,*)
    write(*,*) "Input X,Y,Z of ring center in Angstrom, e.g. 3.2,1.3,-2.1"
    write(*,*) "If press ENTER button directly, geometric center will be used"
    read(*,"(a)") c2000tmp
    if (c2000tmp==" ") then
        cenpos(1)=sum(a(tmparr(:))%x)/ntmp
        cenpos(2)=sum(a(tmparr(:))%y)/ntmp
        cenpos(3)=sum(a(tmparr(:))%z)/ntmp
        write(*,"(' Center:',3f10.5,' Angstrom')") cenpos*b2a
    else
        read(c2000tmp,*) cenpos(:)
    end if
    write(*,"(/,a)") " Input distance between the point below the plane and the plane center in Angstrom"
    write(*,*) "For example, 9.8"
    write(*,*) "If press ENTER button directly, 10 Angstrom will be used, which is usually sufficient"
    read(*,"(a)") c80tmp
    if (c80tmp==" ") then
		distbelow=10
    else
	    read(c80tmp,*) distbelow
    end if
    distbelow=distbelow/b2a
    endpt1(:)=cenpos(:)-distbelow*uvec(:)
    write(*,"(a)") " Input distance between the point above the plane and the plane center in Angstrom"
    write(*,*) "For example, 9.8"
    write(*,*) "If press ENTER button directly, 10 Angstrom will be used, which is usually sufficient"
    read(*,"(a)") c80tmp
    if (c80tmp==" ") then
		distabove=10
    else
	    read(c80tmp,*) distabove
    end if
    distabove=distabove/b2a
    endpt2(:)=cenpos(:)+distabove*uvec(:)
    write(*,"(' Starting point:',3f10.5,' Angstrom')") endpt1(:)*b2a
    write(*,"(' Ending point:  ',3f10.5,' Angstrom')") endpt2(:)*b2a
end if

write(*,*)
write(*,*) "How many points evenly distributing in the scanning line? e.g. 50"
npttmp=nint(dsqrt(sum((endpt2-endpt1)**2))*b2a/0.1)
write(*,"(a,i6,a)") " If press ENTER button directly,",npttmp," points will be used, which is fine enough"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
	npt=npttmp
else
	read(c80tmp,*) npt
end if

allocate(ptxyz(3,npt),ptpos(npt),compin(npt),compout(npt))
tmpvec(:)=endpt2(:)-endpt1(:)
dist=dsqrt(sum(tmpvec(:)**2))
spacing=dist/(npt-1)
tmpvec=tmpvec/dsqrt(sum(tmpvec**2))
do ipt=1,npt
    ptxyz(:,ipt)=endpt1(:)+(ipt-1)*spacing*tmpvec(:)
    if (iway==1) then
        ptpos(ipt)=(ipt-1)*spacing
    else if (iway==2) then
        ptpos(ipt)=(ipt-1)*spacing-distbelow
    end if
    !write(*,"(i5,3f12.6)") ipt,ptxyz(:,ipt)*b2a
end do

do while(.true.)
    write(*,*)
    write(*,*) "0 Exit"
    write(*,*) "1 Generate Gaussian input file for NICS-1D scanning"
    write(*,*) "2 Load Gaussian output file of NICS-1D scanning"
    read(*,*) isel

    if (isel==0) then
        return
    else if (isel==1) then
        write(*,*) "Input the path of Gaussian template file of performing NMR task"
        write(*,*) "e.g. D:\Aqours\Mari\shiny.gjf"
        write(*,"(a)") " Note: In this file, the coordinate part should be recorded as [geometry], which will be automatically replaced with current geometry"
        do while(.true.)
	        read(*,"(a)") c2000tmp
	        inquire(file=c2000tmp,exist=alive)
	        if (alive) exit
	        write(*,*) "Cannot find the file, input again!"
        end do

        open(10,file=c2000tmp,status="old")
        open(11,file="NICS_1D.gjf",status="replace")
        do while (.true.)
            read(10,"(a)",iostat=ierror) c2000tmp
            if (ierror/=0) exit
            if (index(c2000tmp,"#")/=0) then
                if (index(c2000tmp,"geom=conn")==0) then
                    write(11,"(a)") trim(c2000tmp)//" geom=connectivity"
                else
                    write(11,"(a)") trim(c2000tmp)
                end if
            else if (index(c2000tmp,"[geometry]")/=0) then
                do iatm=1,ncenter
                    write(11,"(a,3f12.6)") ind2name(a(iatm)%index),a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
                end do
                do ipt=1,npt
                    write(11,"('Bq',3f12.6)") ptxyz(:,ipt)*b2a
                end do
                write(11,*)
                do itmp=1,ncenter+npt
                    write(11,"(i6)") itmp
                end do
            else
                write(11,"(a)") trim(c2000tmp)
            end if
        end do
        write(11,*)
        write(11,*)
        close(10)
        close(11)
        write(*,"(a)") " NICS_1D.gjf has been generated in current folder! Please check it, and then run it by Gaussian manually"
    else if (isel==2) then
        write(*,*) "Input path of Gaussian output file of NICS-1D scanning task"
        write(*,*) "e.g. D:\Aqours\Mari\shiny.out"
        do while(.true.)
	        read(*,"(a)") c2000tmp
	        inquire(file=c2000tmp,exist=alive)
	        if (alive) exit
	        write(*,*) "Cannot find the file, input again!"
        end do
        exit
    end if
end do

!Load shielding tensor of ghost atoms
open(10,file=c2000tmp,status="old")
call loclabel(10,"Isotropic =",ifound)
if (ifound==0) then
	close(10)
    write(*,"(a)") " Error: Unable to find magnetic shielding tensor in this file! Please check keywords. Press ENTER button to return"
    read(*,*)
    return
end if
allocate(pttens(3,3,npt),ptNICSiso(npt),ptNICSaniso(npt))
call loclabel(10,"Bq",ifound,0)
do ipt=1,npt
    write(*,"(' Loading Bq',i8)") ipt
	read(10,"(a)") c80tmp
    read(c80tmp(26:),*) ptNICSiso(ipt)
    read(c80tmp(52:),*) ptNICSaniso(ipt)
    read(10,"(6x,f11.4,6x,f11.4,6x,f11.4)") pttens(1,1,ipt),pttens(1,2,ipt),pttens(1,3,ipt)
    read(10,"(6x,f11.4,6x,f11.4,6x,f11.4)") pttens(2,1,ipt),pttens(2,2,ipt),pttens(2,3,ipt)
    read(10,"(6x,f11.4,6x,f11.4,6x,f11.4)") pttens(3,1,ipt),pttens(3,2,ipt),pttens(3,3,ipt)
    read(10,*)
end do
close(10)
pttens=-pttens
ptNICSiso=-ptNICSiso
write(*,*) "Loading finished!"

allocate(ptNICS(npt))
ptNICS=0
do ipt=1,npt !By default, take the component along the scanning line
    ptNICS(ipt)=prjmat(pttens(:,:,ipt),uvec(:))
end do
icomp=4

do while(.true.)
    write(*,*)
    call menutitle("Post-processing menu",10,1)
    write(*,*) "-3 Invert direction of X-axis"
    write(*,*) "-2 Multiply NICS data by a factor"
    if (icomp==-1) write(*,*) "-1 Select component of NICS, current: Anisotropy"
    if (icomp==0) write(*,*) "-1 Select component of NICS, current: Isotropy"
    if (icomp==1) write(*,*) "-1 Select component of NICS, current: XX component"
    if (icomp==2) write(*,*) "-1 Select component of NICS, current: YY component"
    if (icomp==3) write(*,*) "-1 Select component of NICS, current: ZZ component"
    if (icomp==4) write(*,*) "-1 Select component of NICS, current: Along the scanning line"
    if (icomp==5) write(*,*) "-1 Select component of NICS, current: Along an inputted vector"
    write(*,*) "0 Exit to main menu"
    write(*,*) "1 Plot NICS curve along the line"
    write(*,*) "2 Save image file of NICS curve along the line"
    write(*,*) "3 Export NICS curve data along the line"
    write(*,*) "4 Calculate integral of NICS along the line"
    write(*,*) "5 Find minima and maxima along NICS curve"
    write(*,*) "6 Calculate FiPC-NICS"
    read(*,*) isel
    
    if (isel==0) then
        return
        
    else if (isel==-3) then
		do ipt=1,floor(npt/2D0)
			tmpval=ptpos(ipt)
            ptpos(ipt)=ptpos(npt+1-ipt)
            ptpos(npt+1-ipt)=tmpval
			tmpvec(:)=ptxyz(:,ipt)
            ptxyz(:,ipt)=ptxyz(:,npt+1-ipt)
            ptxyz(:,npt+1-ipt)=tmpvec(:)
        end do
        write(*,*) "Done!"
        
    else if (isel==-2) then
        write(*,*) "Input the value to be multiplied to NICS data, e.g. 2.5"
        write(*,*) "Obviously, if you input -1, then the sign of NICS will be inverted"
        read(*,*) tmpval
        ptNICS=ptNICS*tmpval
        pttens=pttens*tmpval
        write(*,*) "Done!"
        
    else if (isel==-1) then
        write(*,*) "Please choose the NICS component:"
        write(*,*) "-1 Anisotropy"
        write(*,*) " 0 Isotropy"
        write(*,*) " 1 XX component"
        write(*,*) " 2 YY component"
        write(*,*) " 3 ZZ component"
        write(*,*) " 4 Along the scanning line"
        write(*,*) " 5 Along an inputted vector"
        read(*,*) icomp
        if (icomp==-1) then
            ptNICS(:)=ptNICSaniso(:)
        else if (icomp==0) then
            ptNICS(:)=ptNICSiso(:)
        else if (icomp==1) then
            ptNICS(:)=pttens(1,1,:)
        else if (icomp==2) then
            ptNICS(:)=pttens(2,2,:)
        else if (icomp==3) then
            ptNICS(:)=pttens(3,3,:)
        else if (icomp==4) then
            ptNICS=0
            do ipt=1,npt
                ptNICS(ipt)=prjmat(pttens(:,:,ipt),uvec(:))
            end do
        else if (icomp==5) then
            write(*,*) "Input the vector, e.g. 3.5,0,-1.2"
            read(*,*) tmpvec(:)
            do ipt=1,npt
                ptNICS(ipt)=prjmat(pttens(:,:,ipt),tmpvec(:))
            end do
        end if
        write(*,*) "Done!"
        
    else if (isel==1.or.isel==2.or.isel==4) then
        rintval=sum(ptNICS(:))*spacing*b2a
        write(*,"(' Integral of NICS component:',f10.2,' ppm*Angstrom')") rintval
        if (isel==1.or.isel==2) then
            curvexmin=minval(ptpos(:))
            curvexmax=maxval(ptpos(:))
            tmpval=(maxval(ptNICS(:))-minval(ptNICS(:)))/10D0
            curveymin=floor(minval(ptNICS(:))-tmpval)
            curveymax=ceiling(maxval(ptNICS(:))+tmpval)
            stepx=nint((curvexmax-curvexmin)*b2a/10D0)
            stepy=nint((curveymax-curveymin)/10D0)
            ilenunit1D=2
            numdiglinex_old=numdiglinex
            numdigliney_old=numdigliney
            numdiglinex=1
            numdigliney=1
            if (isel==1) then
                call drawcurve(ptpos(:),ptNICS(:),npt,curvexmin,curvexmax,stepx,curveymin,curveymax,stepy,"show",0D0,0D0,"NICS (ppm)")
            else
                graphformat_old=graphformat
                graphformat="pdf "
                call setfil("NICS_1D.pdf") !The file name of saved image file may have been modified, recover to default one
                call drawcurve(ptpos(:),ptNICS(:),npt,curvexmin,curvexmax,stepx,curveymin,curveymax,stepy,"save",0D0,0D0,"NICS (ppm)")
                write(*,*) "Graphical file has been saved to NICS_1D.pdf in current folder"
                call setfil("dislin."//trim(graphformat)) !Recover to default one
            end if
            numdiglinex=numdiglinex_old
            numdigliney=numdigliney_old
            graphformat=graphformat_old
        end if
        
    else if (isel==3) then
        write(*,*) "Input path of the file to export, e.g. D:\Genjitsu\no\Yohane.txt"
        write(*,*) "If press ENTER button directly, data will be exported to NICS_1D.txt in current folder"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") c2000tmp="NICS_1D.txt"
        open(10,file=c2000tmp,status="replace")
        do ipt=1,npt
            write(10,"(i6,4f12.6,f10.4)") ipt,ptxyz(:,ipt)*b2a,ptpos(ipt)*b2a,ptNICS(ipt)
        end do
        close(10)
        write(*,*) "Exporting finished!"
        write(*,*) "Column 1: Point index"
        write(*,*) "Column 2/3/4: X/Y/Z position of the point in Angstrom"
        write(*,*) "Column 5: Distance of the point relative to starting point in Angstrom"
        write(*,*) "Column 6: NICS data at the point"
        
    else if (isel==5) then
		call showcurveminmax(npt,ptpos(:),ptNICS(:),2)
        
    else if (isel==6) then !FiPC-NICS
        if (     abs(uvec(2))<1D-3.and.abs(uvec(3))<1D-3) then
			idir=1 !X scan
            write(*,*) "The scanning direction is found to be X"
        else if (abs(uvec(1))<1D-3.and.abs(uvec(3))<1D-3) then
			idir=2 !Y scan
            write(*,*) "The scanning direction is found to be Y"
        else if (abs(uvec(1))<1D-3.and.abs(uvec(2))<1D-3) then
			idir=3 !Z scan
            write(*,*) "The scanning direction is found to be Z"
        else
			write(*,"(a)") " Error: This function is only available when scanning direction is X or Y or Z, however, the current scanning direction does not belong to any of them"
            write(*,*) "Press ENTER button to return"
			read(*,*)
            cycle
        end if
        open(10,file="FiPC-NICS.txt",status="replace")
		do ipt=1,npt
            if (idir==1) then
				compin(ipt)=(pttens(2,2,ipt)+pttens(3,3,ipt))/3D0
				compout(ipt)=pttens(1,1,ipt)/3D0
            else if (idir==2) then
				compin(ipt)=(pttens(1,1,ipt)+pttens(3,3,ipt))/3D0
				compout(ipt)=pttens(2,2,ipt)/3D0
            else if (idir==3) then
				compin(ipt)=(pttens(1,1,ipt)+pttens(2,2,ipt))/3D0
				compout(ipt)=pttens(3,3,ipt)/3D0
            end if
			write(10,"(i5,f10.3,2f12.6)") ipt,ptpos(ipt)*b2a,compin(ipt),compout(ipt)
		end do
		close(10)
        write(*,*)
        write(*,*) "FiPC-NICS.txt has been generated in current folder, meaning of each column:"
        write(*,*) "Column 1: Point index"
        write(*,*) "Column 2: Scanning distance in Angstrom"
        write(*,*) "Column 3: In-plane component of NICS"
        write(*,*) "Column 4: Out-of-plane component of NICS"
        !Determine position for evaluating FiPC-NICS
        do ipt=2,npt
			if (compin(ipt-1)*compin(ipt)<0) exit
        end do
        if (ipt==npt+1) then
			write(*,*)
			write(*,*) "Error: Unable to determine the position for evaluating FiPC-NICS!"
        else
			tmpval = compout(ipt-1) + (compout(ipt)-compout(ipt-1))*compin(ipt-1)/(compin(ipt-1)-compin(ipt))
            tmpr = ptpos(ipt-1) + (ptpos(ipt)-ptpos(ipt-1))*compin(ipt-1)/(compin(ipt-1)-compin(ipt))
			write(*,"(/,' FiPC-NICS is',f12.6,' ppm, at',f8.3,' Angstrom')") tmpval,tmpr*b2a
        end if
    end if
end do

end subroutine