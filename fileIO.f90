!!--------- Use suitable routine to read in the file, infomode=0/1 show/don't show info
subroutine readinfile(thisfilename,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) thisfilename
character fchname*200,moldenname*200,c80tmp*20,filenameonly*200
integer infomode,inamelen

call path2filename(thisfilename,filenameonly) !Remove folder part

iresinfo=0 !First assume residue information is not available from input file (but will set to 1 in the subroutine of loading pdb/pqr/gro)

inamelen=len_trim(thisfilename)
if (infomode==0) write(*,*) "Please wait..."
if (thisfilename(inamelen-2:inamelen)=="fch".or.thisfilename(inamelen-3:inamelen)=="fchk"&
.or.thisfilename(inamelen-3:inamelen)=="FCH".or.thisfilename(inamelen-3:inamelen)=="FChk"&
.or.thisfilename(inamelen-3:inamelen)=="FCHK") then
	call readfch(thisfilename,infomode)
else if (thisfilename(inamelen-3:inamelen)=="mwfn") then
    call readmwfn(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="wfn".or.thisfilename(inamelen-2:inamelen)=="WFN") then
	call readwfn(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="wfx".or.thisfilename(inamelen-2:inamelen)=="WFX") then
	call readwfx(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="chg".or.thisfilename(inamelen-2:inamelen)=="CHG") then
	call readchg(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="pqr".or.thisfilename(inamelen-2:inamelen)=="PQR") then
	call readpdb_pqr(thisfilename,infomode,1)
else if (thisfilename(inamelen-2:inamelen)=="pdb".or.thisfilename(inamelen-2:inamelen)=="PDB") then
	call readpdb_pqr(thisfilename,infomode,0)
else if (thisfilename(inamelen-2:inamelen)=="xyz".or.thisfilename(inamelen-2:inamelen)=="XYZ") then
	call readxyz(thisfilename,infomode,1)
else if (thisfilename(inamelen-1:inamelen)=="31") then
	call read31(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="grd") then
	call readgrd(thisfilename,infomode,0)
else if (thisfilename(inamelen-2:inamelen)=="vti") then
	call readvti(thisfilename,infomode,0)
else if (thisfilename(inamelen-1:inamelen)=="dx") then
	call readdx(thisfilename,infomode,0)
else if (thisfilename(inamelen-2:inamelen)=="gro") then
	call readgro(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="cub".or.thisfilename(inamelen-3:inamelen)=="cube") then
	call readcube(thisfilename,infomode,0)
else if (thisfilename(inamelen-2:inamelen)=="gms") then
	call readgms(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="mol".or.thisfilename(inamelen-2:inamelen)=="sdf") then
	call readmol(thisfilename,infomode)
else if (thisfilename(inamelen-3:inamelen)=="mol2") then
	call readmol2(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="gjf".or.thisfilename(inamelen-2:inamelen)=="com") then
	call readgjf(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="mop".or.thisfilename(inamelen-2:inamelen)=="MOP") then
	call readmop(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="cif") then
	call readcif(thisfilename,infomode)
else if (index(thisfilename,".molden")/=0.or.index(thisfilename,".molden.input")/=0.or.index(thisfilename,"molden.inp")/=0.or.filenameonly=="MOLDEN") then
	!ORCA uses .molden.input as suffix, Dalton uses molden.inp
	call readmolden(thisfilename,infomode)
else if (thisfilename(inamelen-2:inamelen)=="chk") then
	if (formchkpath==" ") then
		write(*,"(a)") " Error: .chk file is not directly supported. You should use formchk to convert it to fch/fchk before loading. &
		&Alternatively, you can set ""formchkpath"" parameter in settings.ini to actual path, so that Multiwfn can directly open .chk file"
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		stop
	end if
	inquire(file=formchkpath,exist=alive)
	if (.not.alive) then
		write(*,"(a)") " Note: Albeit ""formchkpath"" parameter in settings.ini has been defined, &
		&the formchk executable file cannot be located, therefore the .chk file cannot be directly opened by Multiwfn"
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		stop
	else
        fchname=thisfilename
		call chk2fch(fchname)
		inquire(file=fchname,exist=alive)
		if (alive) then
            !Check file content, conversion failure may also result in file, but the content is empty
            open(10,file=fchname,status="old")
            read(10,"(a)",iostat=ierror) c80tmp
            close(10)
            if (c80tmp==" ".or.ierror/=0) then
			    write(*,*) "Error: formchk conversion failed!"
			    write(*,*) "Press ENTER button to exit"
			    read(*,*)
			    stop
            end if
			call readfch(fchname,infomode)
			call delfile(fchname)
		else
			write(*,*) "Error: formchk conversion failed!"
			write(*,*) "Press ENTER button to exit"
			read(*,*)
			stop
		end if
	end if
else if (thisfilename(inamelen-2:inamelen)=="gbw") then
	if (orca_2mklpath==" ") then
		write(*,"(a)") " Error: .gbw file is not directly supported. You should use orca_2mkl to convert it to .molden before loading. &
		&Alternatively, you can set ""orca_2mklpath"" parameter in settings.ini to actual path, so that Multiwfn can directly open .gbw file"
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		stop
	end if
	inquire(file=orca_2mklpath,exist=alive)
	if (.not.alive) then
		write(*,"(a)") " Note: Albeit ""orca_2mklpath"" parameter in settings.ini has been defined, &
		&the orca_2mkl executable file cannot be located, therefore the .gbw file cannot be directly opened by Multiwfn"
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		stop
	else
        moldenname=thisfilename
		call gbw2molden(moldenname)
		inquire(file=moldenname,exist=alive)
		if (alive) then
			call readmolden(moldenname,infomode)
			call delfile(moldenname)
		else
			write(*,*) "Error: orca_2mkl conversion failed!"
			write(*,*) "Press ENTER button to exit"
			read(*,*)
			stop
		end if
	end if
else if (index(thisfilename,"POSCAR")/=0.and.index(thisfilename,".in")==0) then !If file name contains POSCAR while do not have .in or .inp, load as VASP POSCAR
    write(*,"(a)") " Because the input file name contains ""POSCAR"", this file will be loaded as the POSCAR file of VASP"
	call readPOSCAR(thisfilename,infomode)
else if (index(filenameonly,"CHGCAR")/=0.or.index(filenameonly,"CHG")/=0.or.index(filenameonly,"ELFCAR")/=0.or.index(filenameonly,"LOCPOT")/=0) then
    call readVASPgrd(thisfilename,infomode)
else !Plain text file
	ifiletype=0
    !Try to load as output file of a code
    open(10,file=thisfilename,status="old")
    call outputprog(10,iprog)
    close(10)
    if (iprog==1) then
        call loadGauout(10,thisfilename)
    else if (iprog==2) then
        call loadORCAout(10,thisfilename)
    end if
    !Try to load as input file of a code
    if (iprog==0) then
        open(10,file=thisfilename,status="old")
        call inputprog(10,iprog)
        close(10)
        if (iprog==1) then !Read as CP2K input file
            if (thisfilename(inamelen-2:inamelen)=="inp".or.thisfilename(inamelen-6:inamelen)=="restart") call readcp2k(thisfilename)
        else if (iprog==2) then !Read as ORCA input file
            call readORCAinp(thisfilename,infomode)
        else if (iprog==3) then !Read as Quantum ESPRESSO input file
            call readQEinp(thisfilename)
        end if
    end if
    !Try to load as Turbomole input file ($coord in the first line)
    open(10,file=thisfilename,status="old")
    read(10,"(a)") c80tmp
    close(10)
    if (index(c80tmp,"$coord")/=0) then
        write(*,"(a)") " Loading the file as coordinate file of Turbomole since the first line is ""$coord"""
        call readTurbomole(thisfilename,infomode)
    end if
end if

if (iresinfo==0.and.allocated(a)) then !Residue information was not loaded, initialize empty content
    a%resid=1
    a%resname=" "
end if

if (allocated(b)) then
    !Determine how to supply EDF information for the file containing GTF information when pseudopotential basis set is involved
    !For .wfn file, EDF has already been added in due case in subroutine "readwfn", and EDF may has been added in subroutine readwfx if EDF is available in the file and readEDF==1
	ifEDF=0
	do iatm=1,ncenter
		if (a(iatm)%charge==0) cycle !Bq atom, ignore
        if (a(iatm)%index/=a(iatm)%charge) then
			ifEDF=1
            exit
        end if
    end do
	if (ifEDF==1.and.nEDFprims==0) then
		if (isupplyEDF==0) then !Do nothing
			continue
		else if (isupplyEDF==1) then !Supply EDF from .wfx file
			call readEDFatmwfx
		else if (isupplyEDF==2) then !Supply EDF from bulit-in library
			call readEDFlib(infomode)
		end if
	end if
end if

call genfragatm !Generate fragatm array containing all atoms

call init_PBC !Initialize some PBC information and settings
if (ifPBC/=0) then
	ishowcell=1 !Show cell
    iorthoview=1 !Use orthographic view
end if

call init_func !Some user-defined real space functions need special initialization for an system

!Current Multiwfn only initializes LIBRETA when ESP is to be calculated, so below codes are commented
!if (if_initlibreta==0) then
!    write(*,*) "Initializing LIBRETA library for ESP evaluation, please wait..."
!    call initlibreta()
!    if_initlibreta=1 !Global variable
!end if
end subroutine


!!------ Assume ifileid is a Gaussian output file, when iloadGaugeom>0, then load the final geometry as well number of electrons
!iloadGaugeom=1: First try to load input orientation, if it is not found, load standard orientation instead
!iloadGaugeom=2: Load standard orientation
subroutine loadGauout(ifileid,gaufilename)
use defvar
use util
implicit real*8 (a-h,o-z)
character gaufilename*200,c80tmp*80,locstr*40
open(ifileid,file=gaufilename,status="old")
if (iloadGaugeom>0) then
    write(*,*) "Trying to load geometry from this file..."
    do itime=1,2
        if (iloadGaugeom==2.and.itime==1) cycle
        if (itime==1) then 
            locstr="Input orientation:"
        else if (itime==2) then 
            locstr="Standard orientation:"
        end if
        nskip=0
        do while(.true.) !Find the final geometry. nskip is the number of labels to be skipped
            call loclabel(ifileid,locstr,ifound,0)
            if (ifound==0) exit
            nskip=nskip+1
            read(ifileid,*)
        end do
        if (nskip>0) then !Found at least once
            call loclabel(ifileid,locstr,ifound)
            call skiplines(ifileid,5)
            ncenter=0
            do while(.true.)
                read(ifileid,"(a)") c80tmp
                if (index(c80tmp,"----")/=0) exit
                ncenter=ncenter+1
            end do
            if (allocated(a)) deallocate(a)
            allocate(a(ncenter))
            rewind(ifileid)
            do iload=1,nskip
                call loclabel(ifileid,locstr,ifound,0)
                read(ifileid,*)
            end do
            call skiplines(ifileid,4)
            do iatm=1,ncenter
                read(ifileid,*) inouse,a(iatm)%index,inouse,a(iatm)%x,a(iatm)%y,a(iatm)%z
                a(iatm)%name=ind2name(a(iatm)%index)
            end do
            a%charge=a%index !Incompatible with ECP case
            a%x=a%x/b2a
            a%y=a%y/b2a
            a%z=a%z/b2a
            if (itime==1) write(*,*) "Geometry (final, input orientation) has been loaded from this file"
            if (itime==2) write(*,*) "Geometry (final, standard orientation) has been loaded from this file"
            exit
        else
            if (itime==1) then
                write(*,"(a)") " Note: Unable to load geometry in input orientation, trying to load geometry in standard orientation instead"
                rewind(ifileid)
            else if (itime==2) then
                write(*,*) "Warning: Failed to load geometry in standard orientation from this file"
            end if
        end if
    end do
    !Load electrons
    call loclabel(ifileid," alpha electrons")
    read(ifileid,"(a)") c80tmp
    read(c80tmp,*) naelec
    read(c80tmp(24:),*) nbelec
    nelec=naelec+nbelec
    write(*,"(a,3i8)") " Number of alpha/beta/total electrons:",nint(naelec),nint(nbelec),nint(nelec)
    !Load net charge and spin multiplicity, so that when save to gjf they will be identical to output file
    call loclabel(ifileid," Multiplicity =",ifound)
    if (ifound==1) then
        read(ifileid,"(a)") c80tmp
        ieqs=index(c80tmp,'=')
        read(c80tmp(ieqs+1:),*) loadcharge
        ieqs=index(c80tmp,'=',back=.true.)
        read(c80tmp(ieqs+1:),*) loadmulti
    end if
else
    write(*,"(a)") " Note: This file seems to be a Gaussian output file. If you set ""iloadGaugeom"" in settings.ini to 1, &
    then Multiwfn will load geometry (final, input orientation) as well as number of electrons from this file"
end if
close(ifileid)
end subroutine



!!------ Assume ifileid is a ORCA output file, when iloadORCAgeom=1, then load the final geometry as well number of electrons
subroutine loadORCAout(ifileid,ORCAfilename)
use defvar
use util
implicit real*8 (a-h,o-z)
character ORCAfilename*200,c80tmp*80
character :: locstr*40="CARTESIAN COORDINATES (ANGSTROEM)"

open(ifileid,file=ORCAfilename,status="old")
if (iloadORCAgeom>0) then
    write(*,*) "Trying to load geometry from this file..."
    nfound=0
    do while(.true.) !Find the final geometry. nskip is the number of labels to be skipped
        call loclabel(ifileid,locstr,ifound,0)
        if (ifound==0) exit
        nfound=nfound+1
        read(ifileid,*)
    end do
    if (nfound>0) then !Found at least once
        call loclabel(ifileid,locstr,ifound)
        call skiplines(ifileid,2)
        ncenter=0
        do while(.true.)
            read(ifileid,"(a)") c80tmp
            if (c80tmp==" ") exit
            ncenter=ncenter+1
        end do
        if (allocated(a)) deallocate(a)
        allocate(a(ncenter))
        rewind(ifileid)
        do iload=1,nfound-1
            call loclabel(ifileid,locstr,ifound,0)
            read(ifileid,*)
        end do
        call loclabel(ifileid,locstr,ifound,0)
        call skiplines(ifileid,2)
        do i=1,ncenter
            read(ifileid,*) a(i)%name,a(i)%x,a(i)%y,a(i)%z
			call elename2idx(a(i)%name,a(i)%index)
        end do
        a%charge=a%index !Incompatible with ECP case
        a%x=a%x/b2a
        a%y=a%y/b2a
        a%z=a%z/b2a
        write(*,*) "Geometry (final) has been loaded from this file"
    else
        write(*,"(a)") " Warning: Failed to load geometry from this file, """//trim(locstr)//""" field cannot be located"
    end if
else
    write(*,"(a)") " Note: This file seems to be an ORCA output file. If you set ""iloadORCAgeom"" in settings.ini to 1, &
    then Multiwfn will load (final) geometry as well as number of electrons from this file"
end if
close(ifileid)
end subroutine



!!-------- Convert inputted chk file to fch/fchk via formchk, return the path of generated (may be failed) .fch/fchk file
subroutine chk2fch(fchname)
use defvar
implicit real*8 (a-h,o-z)
character(len=200) fchname,c500tmp,fchnamenew
itmp=index(fchname,'.',back=.true.)
if (isys==1) then
    fchnamenew=fchname(:itmp)//"fch"
	c500tmp=trim(formchkpath)//' "'//trim(fchname)//'" "'//trim(fchnamenew)//'" > NUL'
else
    fchnamenew=fchname(:itmp)//"fchk"
	c500tmp=trim(formchkpath)//' "'//trim(fchname)//'" "'//trim(fchnamenew)//'" > /dev/null'
end if
write(*,"(a)") " Running: "//trim(c500tmp)
call system(trim(c500tmp))
fchname=fchnamenew
end subroutine



!!-------- Convert gbw to molden via orca_2mkl, return the path of generated (may be failed) .molden file
subroutine gbw2molden(gbwname)
use defvar
implicit real*8 (a-h,o-z)
character(len=200) gbwname,c200tmp
inamelen=len_trim(gbwname)
ico=index(gbwname,'.',back=.true.)
if (isys==1) then
	c200tmp=trim(orca_2mklpath)//' "'//trim(gbwname(:ico-1))//'" -molden > NUL'
else
	c200tmp=trim(orca_2mklpath)//' "'//trim(gbwname(:ico-1))//'" -molden > /dev/null'
end if
write(*,"(a)") " Running: "//trim(c200tmp)
call system(trim(c200tmp))
gbwname=trim(gbwname(:ico))//"molden.input"
end subroutine



!!-----------------------------------------------------------------
!!------------------------- Read Gaussian formatted check file
!Some temporary arrays, including shelltype, shell2atom, shellcon, primexp, concoeff will be stored to corresponding global arrays after some revisions
!infomode=0 means output info, =1 silent
subroutine readfch(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character selectyn,c80*80,fchtitle*79,levelstr*80
integer,allocatable :: shelltype(:),shell2atom(:),shellcon(:) !Degree of shell contraction
integer,allocatable :: tmparrint(:)
real*8,allocatable :: primexp(:),concoeff(:),SPconcoeff(:),amocoeff(:,:),bmocoeff(:,:),tmpmat(:,:),tmparr(:)
integer :: s2f(-5:5,21)=0 !Give shell type & orbital index to get functype
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
!For backing up spherical basis functions
integer,allocatable :: shelltype5D(:),MOtype5D(:)
real*8,allocatable :: CObasa5D(:,:),CObasb5D(:,:),MOocc5D(:),MOene5D(:),CO5D(:,:)
integer,allocatable :: shelltype6D(:) !Temporarily store Cartesian shell types containing SP shell
real*8,external :: normgau
ifiletype=1
imodwfn=0
s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
s2f(-1,1:4)=(/ 1,2,3,4 /)
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /) !Note: The sequence of f functions in Multiwfn is not identical to .fch, so convert here. While spdgh are identical
s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)
call gensphcartab(1,conv5d6d,conv7f10f,conv9g15g,conv11h21h)

open(10,file=name,status="old")
read(10,"(a)") fchtitle
if (infomode==0) write(*,*) "Loading various information of the wavefunction"
isaveNO=0
isaveNBOocc=0
isaveNBOene=0
if (index(fchtitle,'saveNBOocc')/=0.or.index(fchtitle,'SaveNBOocc')/=0) isaveNBOocc=1
if (index(fchtitle,'saveNBOene')/=0.or.index(fchtitle,'SaveNBOene')/=0) isaveNBOene=1
if (index(fchtitle,'saveNO')/=0.or.index(fchtitle,'SaveNO')/=0) isaveNO=1
if ((isaveNBOocc==1.or.isaveNBOene==1).and.infomode==0) write(*,*) "The file contains NBO information"
if (isaveNO==1.and.infomode==0) write(*,*) "The file contains natural orbitals information"
read(10,"(a)") levelstr
if (index(levelstr,"  ZDO")/=0) then
    write(*,"(/,a)") " Warning: It seems that this fch file comes from semi-empirical calculation, the wavefunction generated using this kind of method &
    is not supported by Multiwfn, all wavefunction analysis results will be meaningless!"
    write(*,"(a)") " Hint: If you want to use Multiwfn to perform wavefunction analysis &
    at a level as cheap as semi-empirical methods, you may consider to use xtb code to conduct GFN-xTB calculation, the .molden file generated by xtb is supported by Multiwfn"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
end if
if (levelstr(11:11)=="R") wfntype=0
if (levelstr(11:11)=="U") wfntype=1
if (levelstr(11:12)=="RO") wfntype=2

call loclabel(10,'Number of electrons')
read(10,"(49x,f12.0)") nelec
read(10,"(49x,f12.0)") naelec
read(10,"(49x,f12.0)") nbelec

!Because closed-shell calculation of OPBE, O3LYP, etc. also has RO prefix, so need to check spin multiplicity. &
!RO must correspond to naelec>nbelec
if (naelec==nbelec.and.wfntype==2) wfntype=0 !Should change to restricted form

if (naelec/=nbelec.and.wfntype==0) wfntype=1 !This is often redundant, but considering that sometimes U is not properly recognized, this maybe useful
if (index(levelstr,"CASSCF")/=0.and.naelec/=nbelec.and.isaveNO/=1) then !Suitable for CASSCF calculation of spin multiplicity >1
    wfntype=2
    write(*,"(a)") " This is a CASSCF wavefunction with spin multiplicity >1, &
    and pseudo-canonical orbitals are recorded, the current wavefunction is treated as ROHF"
end if
if (isaveNBOocc==1.or.isaveNBOene==1.or.isaveNO==1) then
	if (wfntype==0) wfntype=3
	if (wfntype==1) wfntype=4
end if

call loclabel(10,'Number of basis functions')
read(10,"(49x,i12)") nbasis
nindbasis=nbasis
call loclabel(10,'Number of independent functions',ifound) !G09
if (ifound==0) call loclabel(10,'Number of independant functions',ifound) !G03
if (ifound==1) read(10,"(49x,i12)") nindbasis !Number of linear independant functions
call loclabel(10,'Number of translation vectors',ifound,maxline=1000) !G09
if (ifound==1) read(10,"(49x,i12)") ifPBC
virialratio=2D0
call loclabel(10,'Virial Ratio',ifound)
if (ifound==1) read(10,"(49x,1PE22.15)") virialratio
totenergy=0D0
call loclabel(10,'Total Energy',ifound) !if no this entry, loclabel return ifound=0, else =1
if (ifound==1) read(10,"(49x,1PE22.15)") totenergy
call loclabel(10,'Atomic numbers')
read(10,"(49x,i12)") ncenter
allocate(tmparrint(ncenter))
read(10,"(6i12)") tmparrint
if (any(tmparrint==0).and.infomode==0) then
	write(*,*)
	write(*,*) "Ghost atoms (Bq) are found in this file. Do you want to load them? (y/n)"
	write(*,"(a)") " Note: If all Bq atoms have corresponding basis functions, then they can be safely loaded. &
	&However, if some of them do not have basis functions, in general they should not be loaded, otherwise problems or crashes may occur &
	&when performing analyses based on wavefunction"
	read(*,"(a)") selectyn
	if (selectyn=='n'.or.selectyn=='N') ncenter=count(tmparrint/=0)
end if
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
a%index=tmparrint(1:ncenter)
a%name=ind2name(a%index)
call loclabel(10,'Nuclear charges') !If ECP was used, nuclear charge /= atomic number
read(10,*)
read(10,"(5(1PE16.8))") (a(i)%charge,i=1,ncenter)
call loclabel(10,'Current cartesian coordinates')
read(10,*)
read(10,"(5(1PE16.8))") (a(i)%x,a(i)%y,a(i)%z,i=1,ncenter)
if (ifPBC>0) then
    call loclabel(10,'Translation vectors')
    read(10,*)
    if (ifPBC==1) read(10,*) cellv1(:)
    if (ifPBC==2) read(10,*) cellv1(:),cellv2(:)
    if (ifPBC==3) read(10,*) cellv1(:),cellv2(:),cellv3(:)
end if
call loclabel(10,'Shell types')
read(10,"(49x,i12)") nshell
allocate(shelltype(nshell))
read(10,"(6i12)") (shelltype(i),i=1,nshell)

!Note that Multiwfn allows Cartesian and spherical harmonic basis functions mixed together. If any basis function is spherical harmonic type, then isphergau=1.
!Only the spherical harmonic ones will be treated specially
if (infomode==0) write(*,"(' The highest angular moment basis functions is ',a)") shtype2name(maxval(abs(shelltype))) 
isphergau=0
if (any(shelltype<=-2)) isphergau=1
if (any(abs(shelltype)>5).and.infomode==0) then
	write(*,"(a)") " Error: GTFs with angular moment higher than h are found, while Multiwfn currently only supports up to h. Press ENTER button to exit"
	read(*,*)
	stop
end if

if (infomode==0) write(*,*) "Loading basis set definition..."
call loclabel(10,'Number of primitives per shell')
read(10,*)
allocate(shellcon(nshell))
read(10,"(6i12)") (shellcon(i),i=1,nshell)
call loclabel(10,'Shell to atom map')
read(10,*)
allocate(shell2atom(nshell))
read(10,"(6i12)") (shell2atom(i),i=1,nshell)
call loclabel(10,'Primitive exponents')
read(10,"(49x,i12)") nprimshell
allocate(primexp(nprimshell))
read(10,"(5(1PE16.8))") (primexp(i),i=1,nprimshell)
call loclabel(10,'Contraction coefficients')
read(10,*)
allocate(concoeff(nprimshell))
read(10,"(5(1PE16.8))") (concoeff(i),i=1,nprimshell)
read(10,"(a)") c80
if (index(c80,"P(S=P) Contraction coefficients")/=0) then
	backspace(10)
	read(10,*)
	allocate(SPconcoeff(nprimshell))
	read(10,"(5(1PE16.8))") (SPconcoeff(i),i=1,nprimshell)
end if

if (infomode==0) write(*,*) "Loading orbitals..."
!Note: Some basis maybe removed by linear dependence checking, hence the number of orbitals in .fch may less than nbasis(always equals to nmo in Multiwfn)
!Hence when reading information involving the number of orbitals in .fch, use nindbasis instead nmo
!The expansion coefficients, energies in those undefined orbitals are all set to zero
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Restricted/restricted open-shell, saveNO/NBO
	nmo=nbasis
	allocate(MOene(nmo))
	allocate(MOocc(nmo))
	allocate(MOtype(nmo))
	allocate(amocoeff(nmo,nbasis))
	MOtype=0
	MOocc=0D0
	MOene=0D0
	amocoeff=0D0
	if (wfntype==0.or.wfntype==3) then
		MOocc(1:nint(nelec/2))=2D0
	else if (wfntype==2) then
		MOocc(1:nint(nbelec))=2D0  !Alpha electrons is always more than beta counterpart
		MOocc(nint(nbelec+1):nint(naelec))=1D0
		MOtype(nint(nbelec+1):nint(naelec))=1
	end if
	call loclabel(10,'Alpha Orbital Energies',ifound) 
	if (ifound==0) call loclabel(10,'orbital energies',ifound) !PSI4 since 1.2 uses "orbital energies" for closed-shell cases
    if (ifound==0) then
        write(*,"(a)") " Error: Unable to find orbital information from your input file! If the file was produced by ONIOM calculation of Gaussian, &
        you should use special way to convert chk to fch file, see http://sobereva.com/wfnbbs/viewtopic.php?pid=2182"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        stop
    end if
	read(10,*)
	read(10,"(5(1PE16.8))") (MOene(i),i=1,nindbasis)
	call loclabel(10,'Alpha MO coefficients',ifound)
	if (ifound==0) call loclabel(10,'MO coefficients') !PSI4 since 1.2 uses "MO coefficients" for closed-shell cases
	read(10,*)
	read(10,"(5(1PE16.8))") ((amocoeff(imo,ibasis),ibasis=1,nbasis),imo=1,nindbasis)
else if (wfntype==1.or.wfntype==4) then !Unrestricted single-determinant or multiconfiguration wavefunction
	!Note that PSI4 since 1.2 uses lower case, this has been taken into account
	nmo=2*nbasis
	allocate(MOene(nmo))
	allocate(MOocc(nmo))
	allocate(MOtype(nmo))
	allocate(amocoeff(nbasis,nbasis))
	allocate(bmocoeff(nbasis,nbasis))
	MOocc=0D0
	MOocc(1:nint(naelec))=1D0
	MOocc(nbasis+1:nbasis+nint(nbelec))=1D0
	MOtype(1:nbasis)=1
	MOtype(nbasis+1:nmo)=2
	MOene=0D0
	amocoeff=0D0
	bmocoeff=0D0
	call loclabel(10,'Alpha Orbital Energies',ifound)
	if (ifound==0) call loclabel(10,'alpha orbital energies') !For PSI4
	read(10,*)
	read(10,"(5(1PE16.8))") (MOene(i),i=1,nindbasis)
	call loclabel(10,'Beta Orbital Energies',ifound)
	if (ifound==0) call loclabel(10,'beta orbital energies',ifound) !For PSI4
	if (ifound==0) then
		write(*,*) "Error: Beta orbital information was not found but expected!"
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		stop
	end if
	read(10,*)
	read(10,"(5(1PE16.8))") (MOene(i),i=nbasis+1,nbasis+nindbasis)

	call loclabel(10,'lpha MO coefficients') !Compatible with PSI4
	read(10,*)
	read(10,"(5(1PE16.8))") ((amocoeff(imo,ibasis),ibasis=1,nbasis),imo=1,nindbasis)
	call loclabel(10,'eta MO coefficients') !Compatible with PSI4
	read(10,*)
	read(10,"(5(1PE16.8))") ((bmocoeff(imo,ibasis),ibasis=1,nbasis),imo=1,nindbasis)
end if

if (isaveNBOocc==1.or.isaveNO==1) then
	MOocc=MOene
	MOene=0D0
end if
if (isaveNBOene==1) MOocc=0D0 !For saveNBO, the automatically determined occupation number is meaningless
where (MOocc==1000) MOocc=0D0 !When saveNBO is used, the latest several occupation/energy of NBO are 1000, we modify them to zero
where (MOene==1000) MOene=0D0

close(10)

!!!!!! Reading have finished, now generate basis information

!Backup spherical basis information (some of them may be Cartesian ones) with 5D suffix (of course, may be actually 7f, 9g, 11h...),
!convert them to Cartesian type temporarily, at final stage recover them back
if (isphergau==1) then
	allocate(shelltype5D(nshell))
	shelltype5D=shelltype
	where (shelltype<=-2) shelltype=-shelltype !Convert to Cartesian type
	nbasis5D=nbasis
	nbasis=0
	do i=1,nshell
		nbasis=nbasis+shtype2nbas(shelltype(i))
	end do
end if
allocate(shelltype6D(nshell)) !Back up shell information of Cartesian basis
shelltype6D=shelltype
nbasisCar=nbasis

!Allocate space for arrays
nprims=0
do i=1,nshell
	nprims=nprims+shtype2nbas(shelltype(i))*shellcon(i)
end do
allocate(b(nprims),CO(nmo,nprims),basshell(nbasis),bascen(nbasis),bastype(nbasis),primstart(nbasis),&
primend(nbasis),primconnorm(nprims),basstart(ncenter),basend(ncenter))
!Fill CObasa and CObasb
if (isphergau==0) then
	allocate(CObasa(nbasis,nbasis))
	CObasa=transpose(amocoeff)
	if (wfntype==1.or.wfntype==4) then
		allocate(CObasb(nbasis,nbasis))
		CObasb=transpose(bmocoeff)
	end if
else if (isphergau==1) then
	allocate(CObasa(nbasis,nbasis),CObasa5D(nbasis5D,nbasis5D))
	CObasa5D=transpose(amocoeff)
	CObasa=0D0
	if (wfntype==1.or.wfntype==4) then
		allocate(CObasb(nbasis,nbasis),CObasb5D(nbasis5D,nbasis5D))
		CObasb5D=transpose(bmocoeff)
		CObasb=0D0
	end if
	!Map 5D coefficient to 6D coefficient
	ipos5D=1
	ipos6D=1
	do ish=1,nshell
		ishtyp5D=shelltype5D(ish)
		ishtyp6D=shelltype(ish)
		numshorb5D=shtype2nbas(ishtyp5D)
		numshorb6D=shtype2nbas(ishtyp6D)
		if (ishtyp5D>=-1) then !S or P or SP or other Cartesian shells
			CObasa(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasa5D(ipos5D:ipos5D+numshorb5D-1,:)
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasb5D(ipos5D:ipos5D+numshorb5D-1,:)			
		else if (ishtyp5D==-2) then
			!5D->6D
			CObasa(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasb5D(ipos5D:ipos5D+4,:))
		else if (ishtyp5D==-3) then
			!7F->10F
			CObasa(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasb5D(ipos5D:ipos5D+6,:))
		else if (ishtyp5D==-4) then
			!9G->15G
			CObasa(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasb5D(ipos5D:ipos5D+8,:))
		else if (ishtyp5D==-5) then
			!11H->21H
			CObasa(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasb5D(ipos5D:ipos5D+10,:))
		end if
		ipos5D=ipos5D+numshorb5D
		ipos6D=ipos6D+numshorb6D
	end do
end if

if (infomode==0) write(*,*) "Converting basis function information to GTF information..."
!Distribute exponent, functype to every GTF and generate CO(:,:) from amocoeff/bmocoeff
!Fill: b,basshell,bascen,bastype,co,primstart,primend,primconnorm
k=1 !Current position of GTF
iexp=1
ibasis=1 !Current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !cycle each shell
	b(k:k+shellcon(i)*shtype2nbas(shelltype(i))-1)%center=shell2atom(i)
	basshell(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=i !!! Set basis attributed to which shell
	bascen(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=shell2atom(i) !!! Set basis attributed to which center
	do j=1,shtype2nbas(shelltype(i)) !cycle each basis(orbital) in each shell
		b(k:k+shellcon(i)-1)%type=s2f(shelltype(i),j)
		bastype(ibasis)=s2f(shelltype(i),j) !!! set basis type
		primstart(ibasis)=k !!! From where the GTFs attributed to ibasis'th basis
		primend(ibasis)=k+shellcon(i)-1 !!! To where the GTFs attributed to ibasis'th basis
        !write(*,*) i,j,ibasis,primstart(ibasis),primend(ibasis)
		do l=1,shellcon(i) !cycle each GTF in each basis in each shell
			b(k)%exp=primexp(iexp+l-1)
			tnormgau=normgau(b(k)%type,b(k)%exp)  !!!Normalization coefficient of GTFs
			if (ifchprog==2) tnormgau=1D0 !In the .fch file of old version of Q-chem (may be older than 5.0), normalization factor of GTF has already been multiplied into contraction coefficient of GTFs
			temp=concoeff(iexp+l-1)  !!!Contraction coefficient of GTFs
			if (shelltype(i)==-1.and.j/=1) temp=SPconcoeff(iexp+l-1)
			primconnorm(k)=temp*tnormgau !Combines contraction and normalization coefficient
			do imo=1,nmo
				if (wfntype==0.or.wfntype==2.or.wfntype==3) then !R or RO
					CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
				else if (wfntype==1.or.wfntype==4) then !U
					if (isphergau==1) then
						if (imo<=nbasis5D) CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
						if (imo>nbasis5D) CO(imo,k)=CObasb(ibasis,imo-nbasis5D)*temp*tnormgau
					else
						if (imo<=nbasis) CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
						if (imo>nbasis) CO(imo,k)=CObasb(ibasis,imo-nbasis)*temp*tnormgau
					end if
				end if
			end do
			k=k+1
		end do
		ibasis=ibasis+1
	end do
	iexp=iexp+shellcon(i)
end do

!Convert information from Cartesian basis to spherical basis
if (isphergau==1) then
	if (iloadasCart==1) then !For special purpose, keep Cartesian basis functions, e.g. convert spherical .fch/.molden to .47 file
		!Currently nbasis and dimension of all matrix (except for CO) correspond to full Cartesian case, but nmo &
		!and related arrays as well as CO still correspond to spherical harmonic case and thus need to be "expanded", the MO energies are artifically set to 100
		allocate(MOocc5D(nmo),MOene5D(nmo),MOtype5D(nmo),CO5D(nmo,nprims))
		MOocc5D=MOocc
		MOene5D=MOene
		MOtype5D=MOtype
		CO5D=CO
		deallocate(MOocc,MOene,MOtype,CO)
		if (wfntype==0.or.wfntype==2.or.wfntype==3) nmo=nbasis !R, RO
		if (wfntype==1.or.wfntype==4) nmo=2*nbasis !U
		allocate(MOocc(nmo),MOene(nmo),MOtype(nmo),CO(nmo,nprims))
		MOocc=0
		MOene=100
		CO=0
		if (wfntype==0.or.wfntype==2.or.wfntype==3) then !R, RO
			MOtype=0
			MOocc(1:nbasis5D)=MOocc5D
			MOene(1:nbasis5D)=MOene5D
			MOtype(1:nbasis5D)=MOtype5D
			CO(1:nbasis5D,:)=CO5D
		else !U
			MOtype(:nbasis)=1
			MOtype(nbasis+1:)=2
			MOocc(:nbasis5D)=MOocc5D(:nbasis5D)
			MOocc(nbasis+1:nbasis+nbasis5D)=MOocc5D(nbasis5D+1:)
			MOene(:nbasis5D)=MOene5D(:nbasis5D)
			MOene(nbasis+1:nbasis+nbasis5D)=MOene5D(nbasis5D+1:)
			MOtype(:nbasis5D)=MOtype5D(:nbasis5D)
			MOtype(nbasis+1:nbasis+nbasis5D)=MOtype5D(nbasis5D+1:)
			CO(:nbasis5D,:)=CO5D(:nbasis5D,:)
			CO(nbasis+1:nbasis+nbasis5D,:)=CO5D(nbasis5D+1:,:)
		end if
		isphergau=0
		
	else !Common case, transform to spherical functions
		if (infomode==0) write(*,*) "Back converting basis function information from Cartesian to spherical type..."
		!Recover spherical Gaussian basis function information
		nbasis=nbasis5D
		shelltype=shelltype5D
		ibasis=1
		do i=1,nshell
			basshell(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=i
			bascen(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=shell2atom(i)
			do j=1,shtype2nbas(shelltype(i))
				bastype(ibasis)=s2f(shelltype(i),j)
				ibasis=ibasis+1
			end do
		end do
		deallocate(CObasa)
		allocate(CObasa(nbasis,nbasis))
		CObasa=CObasa5D
		if (wfntype==1.or.wfntype==4) then
			deallocate(CObasb)
			allocate(CObasb(nbasis,nbasis))
			CObasb=CObasb5D
		end if
	end if
end if

!Split SP shell to S and P shells
noldshell=nshell
noldprimshell=nprimshell
ibasis=1
do i=1,nshell !Count how many basis shells and primitive shells after splitting SP as S and P, and meantime update basshell
	if (shelltype(i)==-1) then
		nshell=nshell+1
		nprimshell=nprimshell+shellcon(i)
		basshell(ibasis+1:nbasis)=basshell(ibasis+1:nbasis)+1 !The shell index of the basis function after current one should be augmented by 1, since current shell is splitted
	end if
	ibasis=ibasis+shtype2nbas(shelltype(i))
end do
allocate(shtype(nshell),shtypeCar(nshell),shcen(nshell),shcon(nshell),primshexp(nprimshell),primshcoeff(nprimshell)) !Global array
jsh=1 !New basis shell index
iprsh=1 !Old primitive shell index
jprsh=1 !New primitive shell index
do ish=1,noldshell !Finally determine global arrays shtype,shcen,shcon,primshexp,primshcoeff shell arrays, in which SP shells are not presented
	ncon=shellcon(ish)
	if (shelltype(ish)/=-1) then !Normal shell
		shtype(jsh)=shelltype(ish)
        shtypeCar(jsh)=shelltype6D(ish)
		shcen(jsh)=shell2atom(ish)
		shcon(jsh)=ncon
		primshexp(jprsh:jprsh+ncon-1)=primexp(iprsh:iprsh+ncon-1)
		primshcoeff(jprsh:jprsh+ncon-1)=concoeff(iprsh:iprsh+ncon-1)
		jsh=jsh+1
		jprsh=jprsh+ncon
	else !SP shell
		shtype(jsh)=0 !S
		shtype(jsh+1)=1 !P
		shtypeCar(jsh)=0 !S
		shtypeCar(jsh+1)=1 !P
		shcen(jsh:jsh+1)=shell2atom(ish)
		shcon(jsh:jsh+1)=ncon
		primshexp(jprsh:jprsh+ncon-1)=primexp(iprsh:iprsh+ncon-1)
		primshexp(jprsh+ncon:jprsh+2*ncon-1)=primexp(iprsh:iprsh+ncon-1)
		primshcoeff(jprsh:jprsh+ncon-1)=concoeff(iprsh:iprsh+ncon-1)
		primshcoeff(jprsh+ncon:jprsh+2*ncon-1)=SPconcoeff(iprsh:iprsh+ncon-1)
		jsh=jsh+2
		jprsh=jprsh+2*ncon
	end if
	iprsh=iprsh+ncon
end do

!Generate basstart and basend
call bascen2basstart_end

!Generate one-particle matrix in basis functions
if (igenP==1) then
	if (infomode==0) then
		if (isaveNO==0) write(*,*) "Generating density matrix based on SCF orbitals..."
		if (isaveNO==1) write(*,*) "Generating density matrix based on natural orbitals..."
	end if
	call genP
end if

if (ifPBC==0) then !For PBC case, Sbas will be generated when need it because it may be relatively expensive
    if (infomode==0) write(*,*) "Generating overlap matrix..."
    call genSbas_curr
    devtmp=abs(sum(Sbas*Ptot)-nint(nelec))
    if (devtmp>0.01D0) then
	    write(*,"(' Deviation of Tr(S*P) to the total number of electrons:',f12.6)") devtmp
	    write(*,"(' Tr(S*P) is',f12.6)") sum(Sbas*Ptot)
	    write(*,"(' Total number of electrons is',i10)") nint(nelec)
	    write(*,"(/,a)") " Warning: The loaded wavefunction is incorrect! That means this file is problematic"
	    write(*,"(a)") " If you really want to proceed, the result may be incorrect"
	    !write(*,"(a)") " If you really want to proceed, press ENTER button now, but the result may be incorrect"
	    !read(*,*)
    end if
end if

!Output summary of present wavefunction
if (infomode==0) then
	write(*,*)
	write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
	write(*,"(' Net charge:',f12.5,'      Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
	write(*,"(' Atoms:',i7,',  Basis functions:',i7,',  GTFs:',i7)") ncenter,nbasis,nprims
	write(*,"(' Total energy:',f21.12,' Hartree,   Virial ratio:',f12.8)") totenergy,virialratio
	if (wfntype==0) then
		write(*,"(' This is a restricted single-determinant wavefunction')")
		write(*,"(' Orbitals from 1 to',i6,' are occupied')") nint(nelec/2)
	else if (wfntype==1) then
		write(*,"(' This is an unrestricted single-determinant wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are occupied')") 1,nbasis,1,nint(naelec)
		write(*,"(' Orbitals from ',i6,' to',i6,' are beta,  from',i6,' to',i6,' are occupied')") nbasis+1,nmo,nbasis+1,nbasis+nint(nbelec)
	else if (wfntype==2) then
		write(*,"(' This is a restricted open-shell wavefunction')")
		write(*,"(' Orbitals from',i6,' to',i6,' are doubly occupied')") 1,nint(nbelec)
		write(*,"(' Orbitals from',i6,' to',i6,' are singly occupied')") nint(nbelec)+1,nint(naelec)
	else if (wfntype==3) then
		write(*,"(' This is a restricted multiconfiguration wavefunction')")
	else if (wfntype==4) then
		write(*,"(' This is an unrestricted multiconfiguration wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are beta')") 1,nbasis,nbasis+1,nmo
	end if
	write(*,"(' Title line of this file: ',a)") trim(fchtitle)
end if

call getHOMOidx !Find out index of HOMO, will be used in some cases
end subroutine


!!--------------------------------------------------------------------
!!-------------- Read Gaussian formatted check file for AdNDP analysis
!CObasa is not read from .fch file, which has already written by AdNDP module, which contains AO expansion coefficients of candidate or saved orbitals
subroutine readfchadndp(fchfilename,iusespin,orbocc,adndpCObas,numorb)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) fchfilename
character c80*80
integer iusespin,numorb
real*8 orbocc(numorb),adndpCObas(nbasis,numorb)
integer,allocatable :: shelltype(:),shell2atom(:),shellcon(:) !Degree of shell contraction
real*8,allocatable :: primexp(:),concoeff(:),SPconcoeff(:)
integer :: s2f(-5:5,21)=0 !Give shell type & orbital index to get functype
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
 !For backing up spherical basis functions
integer,allocatable :: shelltype5D(:)
real*8,allocatable :: CObasa5D(:,:)
real*8,external :: normgau
nmo=numorb
s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
s2f(-1,1:4)=(/ 1,2,3,4 /)
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /)
s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)
call gensphcartab(1,conv5d6d,conv7f10f,conv9g15g,conv11h21h)

open(10,file=fchfilename,status="old")

call loclabel(10,'Nuclear charges') !If ECP was used, nuclear charge /= atomic number
read(10,*)
read(10,"(5(1PE16.8))") (a(i)%charge,i=1,ncenter)
call loclabel(10,'Current cartesian coordinates')
read(10,*)
read(10,"(5(1PE16.8))") (a(i)%x,a(i)%y,a(i)%z,i=1,ncenter)

call loclabel(10,'Shell types')
read(10,"(49x,i12)") nshell
allocate(shelltype(nshell))
read(10,"(6i12)") (shelltype(i),i=1,nshell)
isphergau=0
if (any(shelltype<=-2)) isphergau=1
call loclabel(10,'Number of primitives per shell')
read(10,*)
allocate(shellcon(nshell))
read(10,"(6i12)") (shellcon(i),i=1,nshell)
call loclabel(10,'Shell to atom map')
read(10,*)
allocate(shell2atom(nshell))
read(10,"(6i12)") (shell2atom(i),i=1,nshell)
call loclabel(10,'Primitive exponents')
read(10,"(49x,i12)") nprimshell
allocate(primexp(nprimshell))
read(10,"(5(1PE16.8))") (primexp(i),i=1,nprimshell)
call loclabel(10,'Contraction coefficients')
read(10,*)
allocate(concoeff(nprimshell))
read(10,"(5(1PE16.8))") (concoeff(i),i=1,nprimshell)
read(10,"(a)") c80
if (index(c80,"P(S=P) Contraction coefficients")/=0) then
	backspace(10)
	read(10,*)
	allocate(SPconcoeff(nprimshell))
	read(10,"(5(1PE16.8))") (SPconcoeff(i),i=1,nprimshell)
end if

if (allocated(MOene)) deallocate(MOene)
if (allocated(MOocc)) deallocate(MOocc)
if (allocated(MOtype)) deallocate(MOtype)
allocate(MOene(nmo),MOocc(nmo),MOtype(nmo))
MOene=0D0
MOocc=orbocc
MOtype=iusespin

nbasis_org=nbasis
if (isphergau==1) then
	allocate(shelltype5D(nshell))
	shelltype5D=shelltype
	where (shelltype<=-2) shelltype=-shelltype !Convert to Cartesian type
	nbasis5D=nbasis
	nbasis=0
	do i=1,nshell
		nbasis=nbasis+shtype2nbas(shelltype(i))
	end do
end if
if (allocated(shtypeCar)) deallocate(shtypeCar)
allocate(shtypeCar(nbasis)) !Store shell information of Cartesian basis into global array, which may be used later
shtypeCar=shelltype
nbasisCar=nbasis

!Allocate space for arrays
nprims=0
do i=1,nshell
	nprims=nprims+shtype2nbas(shelltype(i))*shellcon(i)
end do
if (.not.allocated(b)) allocate(b(nprims))
if (allocated(CO)) deallocate(CO)
allocate(CO(nmo,nprims))

if (isphergau==0) then
	if (allocated(CObasa)) deallocate(CObasa)
	allocate(CObasa(nbasis,nmo))
	CObasa=adndpCObas
else if (isphergau==1) then
	if (allocated(CObasa)) deallocate(CObasa)
	allocate(CObasa(nbasis,nmo),CObasa5D(nbasis5D,nmo))
	CObasa5D=adndpCObas
	CObasa=0D0
	
	!Map 5D coefficient to 6D coefficient
	ipos5D=1
	ipos6D=1
	do ish=1,nshell
		ishtyp5D=shelltype5D(ish)
		ishtyp6D=shelltype(ish)
		numshorb5D=shtype2nbas(ishtyp5D)
		numshorb6D=shtype2nbas(ishtyp6D)
		if (ishtyp5D==0.or.ishtyp5D==1.or.ishtyp5D==-1) then !S or P or SP
			CObasa(ipos6D:ipos6D+numshorb6D-1,1:nmo)=CObasa5D(ipos5D:ipos5D+numshorb5D-1,1:nmo)
		else if (ishtyp5D==-2) then
			!5D->6D
			CObasa(ipos6D:ipos6D+5,1:nmo)=matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,1:nmo))
		else if (ishtyp5D==-3) then
			!7F->10F
			CObasa(ipos6D:ipos6D+9,1:nmo)=matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,1:nmo))
		else if (ishtyp5D==-4) then
			!9G->15G
			CObasa(ipos6D:ipos6D+14,1:nmo)=matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,1:nmo))
		else if (ishtyp5D==-5) then
			!11H->21H
			CObasa(ipos6D:ipos6D+20,1:nmo)=matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,1:nmo))
		end if
		ipos5D=ipos5D+numshorb5D
		ipos6D=ipos6D+numshorb6D
	end do
end if

!Distribute exponent, functype to every GTF and generate CO(:,:) from CObasa
k=1 !current position of GTF
iexp=1
ibasis=1 !current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !cycle each shell
	b(k:k+shellcon(i)*shtype2nbas(shelltype(i))-1)%center=shell2atom(i)
	do j=1,shtype2nbas(shelltype(i)) !cycle each basis(orbital) in each shell
		b(k:k+shellcon(i)-1)%type=s2f(shelltype(i),j)
		do l=1,shellcon(i) !cycle each GTF in each basis in each shell
			b(k)%exp=primexp(iexp+l-1)
			tnormgau=normgau(b(k)%type,b(k)%exp)
			temp=concoeff(iexp+l-1)
			if (shelltype(i)==-1.and.j/=1) temp=SPconcoeff(iexp+l-1)
			do imo=1,nmo
				CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
			end do
			k=k+1
		end do
		ibasis=ibasis+1
	end do
	iexp=iexp+shellcon(i)
end do

close(10)

nbasis=nbasis_org
end subroutine



!!------------------------------------------------------------------------
!!------------- Read density matrix from .fch file
!fileid is the file id that may be used, fchname is path of fch file
!For closed-shell wavefunction, Ptot will be updated; For open-shell, Palpha and Pbeta will also be updated.
subroutine readfchdensmat(fileid,fchname)
use defvar
use util
real*8,allocatable :: Pspin(:,:)
integer fileid
character(len=*) fchname
open(fileid,file=fchname,status="old")
call loclabel(10,"Total SCF Density")
read(fileid,*)
read(fileid,*) ((Ptot(i,j),j=1,i),i=1,nbasis)
if (allocated(Pbeta)) then !Must be an open shell system
    call loclabel(fileid,"Spin SCF Density",ifound)
    if (ifound==1) then
        allocate(Pspin(nbasis,nbasis))
        read(fileid,*)
        read(fileid,*) ((Pspin(i,j),j=1,i),i=1,nbasis)
        Palpha=(Ptot+Pspin)/2
        Pbeta=(Ptot-Pspin)/2
    else
        write(*,*) "Error: Unable to find ""Spin SCF Density"", Palpha and Pbeta are not updated"
        return
    end if
end if
do i=1,nbasis
	do j=i+1,nbasis
		Ptot(i,j)=Ptot(j,i)
    end do
end do
if (allocated(Pbeta)) then
    do i=1,nbasis
	    do j=i+1,nbasis
		    Palpha(i,j)=Palpha(j,i)
		    Pbeta(i,j)=Pbeta(j,i)
        end do
    end do
end if
close(fileid)
end subroutine



!!------------------------------------------------------------------------
!!------------- Read .chg file that only contain atomic charge information
subroutine readchg(name,infomode) !infomode=0 means output info, =1 silent
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
integer i
real*8 dipx,dipy,dipz,tdip
ifiletype=4
open(10,file=name,status="old")
ncenter=totlinenum(10,1)
allocate(a(ncenter))
do i=1,ncenter
	read(10,*) a(i)%name,a(i)%x,a(i)%y,a(i)%z,a(i)%charge
    call elename2idx(a(i)%name,a(i)%index)
end do
close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
dipx=sum(a(:)%x*a(:)%charge)
dipy=sum(a(:)%y*a(:)%charge)
dipz=sum(a(:)%z*a(:)%charge)
if (infomode==0) then
	write(*,"(' Totally',i8,' atoms')") ncenter
	write(*,"(' Sum of all atomic charges:',f12.6)") sum(a%charge)
	write(*,"(' Component of electric dipole moment:')")
	write(*,"(' X=',f12.6,' a.u.  (',f12.6,' Debye )')") dipx,dipx*au2debye
	write(*,"(' Y=',f12.6,' a.u.  (',f12.6,' Debye )')") dipy,dipy*au2debye
	write(*,"(' Z=',f12.6,' a.u.  (',f12.6,' Debye )')") dipz,dipz*au2debye
	tdip=dsqrt(dipx**2+dipy**2+dipz**2)
	write(*,"(' Total electric dipole moment:',f12.6,' a.u.  (',f12.6,' Debye )')") tdip,tdip*au2debye
end if
end subroutine



!!----------------------------------------------------------------
!!-------------------- Read .pdb or .pqr file --------------------
! infomode=0: Output summary, =1: do not
! ipqr=0: read as pdb, =1: read as pqr
subroutine readpdb_pqr(name,infomode,ipqr)
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode,ipqr
character(len=*) name
character test*6,tmpname*4,tmpname_up*4,element*3
character c80tmp*80
ifiletype=5
open(10,file=name,status="old")
ncenter=0
do while(.true.)
	read(10,"(a6)",iostat=ierror) test
	if (ierror/=0) exit
	if (test(1:3)=="END") exit
	if (test=="HETATM".OR.test=="ATOM  ") ncenter=ncenter+1
end do
rewind(10)
allocate(a(ncenter))
i=0
do while(.true.)
	read(10,"(a6)",iostat=ierror) test
	if (ierror/=0) exit
	if (test(1:3)=="END") exit
	if (test=="HETATM".or.test=="ATOM  ") then
		backspace(10)
		i=i+1
		read(10,"(12x,a4,1x,a3,2x,i4,4x,3f8.3,22x,a2)") tmpname,a(i)%resname,a(i)%resid,a(i)%x,a(i)%y,a(i)%z,element
		if (ipqr==1) then !Use free-format to load atomic charge
			backspace(10)
			read(10,"(a)") c80tmp
			read(c80tmp(55:76),*) a(i)%charge
		end if
		element=adjustl(element)
		ifound=0
		!Use "element" term to determine actual element
		do j=0,nelesupp
			if (ind2name_up(j)==element(1:2).or.ind2name(j)==element(1:2)) then
				a(i)%index=j
				ifound=1
				exit
			end if
		end do
		!"Element" term is missing, use atomic name to determine element
		if (ifound==0) then
			tmpname_up=adjustl(tmpname)
			call strlc2uc(tmpname_up)
            ifound=0
			do j=1,nelesupp !Recognize e.g. CL1, Li
                if (ind2name_up(j)==tmpname_up(1:2)) then
                    a(i)%index=j
					ifound=1
					exit
                end if
            end do
            if (ifound==0) then
			    do j=1,nelesupp
				    if (ind2name_up(j)==tmpname_up(1:1)//' ' .or.  & !Recognize e.g. C, N11, Li
				    ((ichar(tmpname_up(1:1))<=57).and.ind2name_up(j)==tmpname_up(2:2)//' ')) then !Recognize such as 1H5*
					    a(i)%index=j
					    ifound=1
					    exit
				    end if
			    end do
            end if
		end if
		if (ifound==0) then
			write(*,"(3a)") " Warning: Found unknown element """,tmpname,""" , assume it is carbon"
			a(i)%index=12
		end if
		if (ipqr==0) a(i)%charge=a(i)%index
		a(i)%name=ind2name(a(i)%index)
	end if
end do
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a

iresinfo=1 !Residue information is available in this file

!Try to load cell information
call loclabel(10,"CRYST1  ",ifound)
if (ifound==1) then
    read(10,"(6x,3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
    asize=asize/b2a
    bsize=bsize/b2a
    csize=csize/b2a
    call abc2cellv(asize,bsize,csize,alpha,beta,gamma)
    ifPBC=3
end if
close(10)

if (infomode==0) then
	write(*,"(' Totally',i8,' atoms')") ncenter
	if (ipqr==1) write(*,"(' Sum of atomic charges:',f12.6)") sum(a%charge)
end if

call guessnelec

end subroutine



!!--------------------------------------------------------
!!-------------------- Read .xyz file --------------------
! infomode=0: Output summary, =1: do not
! iopen=0 means don't open and close file in this routine, used to continuously read trjectory. =1 means do these
!Rule of detecting atom element indices: If there is a .pdb file with same path, use element from it to guess index if element is not blank, &
!otherwise guess index from atom name in xyz
subroutine readxyz(name,infomode,iopen) 
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode,iopen
character(len=*) name
character titleline*200,pdbpath*200,teststr*6,c200tmp*200

ifiletype=5
if (iopen==1) open(10,file=name,status="old")
ncenter=0
read(10,*) ncenter
read(10,"(a)") titleline
itmp=index(titleline,'spinmulti')
if (itmp/=0) then
    read(titleline(itmp+9:),*) loadmulti
    write(*,"(' Spin multiplicity:',i4)") loadmulti
end if
itmp=index(titleline,'charge ')
if (itmp/=0) then
    read(titleline(itmp+6:),*,iostat=ierror) loadcharge
    if (ierror/=0) then
        loadcharge=-99 !Default
    else
        write(*,"(' Net charge:',i4)") loadcharge
    end if
end if

!Load PBC information. If Multiwfn loaded a periodic system, then the exported xyz file will carry this kind of information
itmp=index(titleline,'Tv_1:')
if (itmp/=0) then
    read(titleline(itmp+5:),*) cellv1(:)
    cellv1=cellv1/b2a
    ifPBC=1
end if
itmp=index(titleline,'Tv_2:')
if (itmp/=0) then
    read(titleline(itmp+5:),*) cellv2(:)
    cellv2=cellv2/b2a
    ifPBC=2
end if
itmp=index(titleline,'Tv_3:')
if (itmp/=0) then
    read(titleline(itmp+5:),*) cellv3(:)
    cellv3=cellv3/b2a
    ifPBC=3
end if
itmp=index(titleline,'attice=') !Extended xyz format use Lattice= to record cell vectors
if (itmp/=0) then
	c200tmp=trim(titleline)
	itmp=index(c200tmp,'"')
    c200tmp(:itmp)=" "
	itmp=index(c200tmp,'"')
    c200tmp(itmp:)=" "
    read(c200tmp,*) cellv1(:),cellv2(:),cellv3(:)
    cellv1=cellv1/b2a
    cellv2=cellv2/b2a
    cellv3=cellv3/b2a
    ifPBC=3
end if

if (allocated(a)) deallocate(a)
allocate(a(ncenter))
ipos=index(name,'.',back=.true.)
pdbpath=name(:ipos-1)//".pdb"
inquire(file=pdbpath,exist=alive)
if (alive) then
    if (infomode==0) write(*,"(a)") " Note: Found "//trim(pdbpath)//", element names in this file will be used instead if given"
    open(11,file=pdbpath,status="old")
    do while(.true.) !Find the first atom section
	    read(11,"(a6)") teststr
	    if (teststr=="HETATM".OR.teststr=="ATOM  ") then
            backspace(11)
            exit
        end if
    end do
end if

do i=1,ncenter
	read(10,*) a(i)%name,a(i)%x,a(i)%y,a(i)%z
    if (alive) then
        do while(.true.) !Read next atom from pdb (the line must have HETATM or ATOM)
            read(11,"(a)") c200tmp
            if (index(c200tmp,"HETATM")==0.and.index(c200tmp,"ATOM")==0) cycle
            read(c200tmp,"(76x,a2)") teststr(1:2) !Load element from pdb file
            if (teststr(1:2)/=" ") a(i)%name=adjustl(teststr(1:2)) !If pdb file doesn't provide element name, still use name in .xyz
            exit
        end do
    end if
    call elename2idx(a(i)%name,idx)
    if (idx/=0) then
        a(i)%index=idx
	else !Only use the first letter of atom name to try to assign. For example OW, HW may be in .xyz file
		do j=1,nelesupp
			if ( ind2name(j)(2:2)==" ".and.a(i)%name(1:1)==ind2name(j)(1:1) ) then
				a(i)%index=j
				exit
			end if
		end do
	    if (j==nelesupp+1) then
		    write(*,"(/,a)") " Warning: Found unknown element """//a(i)%name//""", assume it is carbon"
		    a(i)%index=12
	    end if
	end if
end do
if (iopen==1) close(10)
if (alive) close(11)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
a%name=ind2name(a%index)
if (infomode==0) then
	write(*,"(' Title line: ',a)") trim(titleline)
	write(*,"(' Totally',i8,' atoms')") ncenter
end if
call guessnelec
end subroutine


!!---------- Read a frame of xyz trajectory form a given ifileid
!Only coordinates will be loaded, other information will not be loaded as "readxyz"
!Commonly readxyz should be used prior to this to assign other information
!This routine will not open/close file, and will not allocate/deallocate "a"
subroutine readxyztrj(ifileid) 
use defvar
implicit real*8 (a-h,o-z)
integer ifileid
character c80tmp*80
read(ifileid,*)
read(ifileid,*)
do i=1,ncenter
	read(ifileid,*) c80tmp,a(i)%x,a(i)%y,a(i)%z
end do
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
end subroutine


!!--------- Get number of frame in .xyz file. The xyz file must have been opened with id=ifileid
subroutine xyzfile_nframe(ifileid,nframe)
use defvar
implicit real*8 (a-h,o-z)
integer ifileid,nframe
character c80tmp*80
nframe=0
do while(.true.)
    read(ifileid,*,iostat=ierror) c80tmp
    if (ierror/=0.or.c80tmp==" ") then
        return
    else
        nframe=nframe+1
    end if
    read(ifileid,*)
    do i=1,ncenter
	    read(ifileid,*) c80tmp,tmpx,tmpy,tmpz
    end do
end do
end subroutine



!!---------------------------------------------------------------------
!!------------ Load Turobmole coordinate file ($coord at the first line)
! infomode=0: Output summary, =1: do not
subroutine readTurbomole(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character c80tmp*80

ifiletype=5 !The same as xyz
open(10,file=name,status="old")

call loclabel(10,"$coord")
read(10,*)
ncenter=0
do while(.true.)
	read(10,*,iostat=ierror) c80tmp
    if (c80tmp==" ".or.index(c80tmp,'$')/=0.or.ierror/=0) exit
    ncenter=ncenter+1
end do
if (infomode==0) write(*,"(' Totally',i8,' atoms')") ncenter

if (allocated(a)) deallocate(a)
allocate(a(ncenter))

call loclabel(10,"$coord")
read(10,*)
do i=1,ncenter
	read(10,*) a(i)%x,a(i)%y,a(i)%z,a(i)%name
    call elename2idx(a(i)%name,a(i)%index)
end do

call loclabel(10,"$periodic",ifound)
if (ifound==1) then
	read(10,*) c80tmp,ifPBC
	call loclabel(10,"$lattice",ifound)
	if (ifound==1) then
		read(10,*)
        cellv1=0
        cellv2=0
        cellv3=0
        if (ifPBC==1) then !On X axis
			read(10,*) cellv1(1)
        else if (ifPBC==2) then !On xy plane
			read(10,*) cellv1(1:2)
			read(10,*) cellv2(1:2)
        else if (ifPBC==3) then
			read(10,*) cellv1
			read(10,*) cellv2
			read(10,*) cellv3
        end if
    else !Turbomole can also use $cell to record cell information by specifying cell lengths and angles, but this is not supported here
		write(*,*) "Error: Unable to find ""$lattice"", now program exit"
        read(*,*)
        stop
	end if
end if

close(10)

a%charge=a%index
a%name=ind2name(a%index)
call guessnelec

end subroutine



!!---------------------------------------------------
!!-------------------- Read .gro --------------------
! infomode=0: Output summary, =1: Do not
subroutine readgro(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character resname*5,resname_up*5,tmpname*5,tmpname_up*5,c200tmp*200
ifiletype=15
open(10,file=name,status="old")
read(10,*)
read(10,*) ncenter
allocate(a(ncenter))
do iatm=1,ncenter
    read(10,"(i5,2a5,i5,3f8.3)") a(iatm)%resid,resname,tmpname,inouse,a(iatm)%x,a(iatm)%y,a(iatm)%z
	tmpname_up=adjustl(tmpname)
    resname_up=adjustl(resname)
	call strlc2uc(tmpname_up)
	call strlc2uc(resname_up)
    ifound=0
    if (tmpname_up==resname_up) then !If resname and atomname are identical, for example NA, CL, CA, then they should be ion, directly recognize
        do i=1,nelesupp
            if (ind2name_up(i)==tmpname_up) then
			    a(iatm)%index=i
			    ifound=1
			    exit
		    end if
        end do
    end if
    if (ifound==0) then !Guess element according to atomname
	    do i=1,nelesupp
		    if (ind2name_up(i)==tmpname_up(1:2)) then !Find matched element name
			    a(iatm)%index=i
			    ifound=1
			    exit
		    end if
	    end do
        if (ifound==0) then
			do i=1,nelesupp
				if (ind2name_up(i)==tmpname_up(1:1)//' '.or. & !Recognize e.g. C, N11, Li
				((ichar(tmpname_up(1:1))<=57).and.ind2name_up(i)==tmpname_up(2:2)//' ')) then !Recognize such as 1H5*
					a(iatm)%index=i
					ifound=1
					exit
				end if
			end do
        end if
    end if
    if (ifound==0) then
	    write(*,"(a)") " Warning: Element of """//trim(tmpname)//""" cannot be recognized, assume it is carbon"
	    a(iatm)%index=12
    end if
    a(iatm)%name=ind2name(a(iatm)%index)
    a(iatm)%resname=resname(1:4)
end do
read(10,"(a)") c200tmp
read(c200tmp,*,iostat=ierror) cellv1(1),cellv2(2),cellv3(3),cellv1(2),cellv1(3),cellv2(1),cellv2(3),cellv3(1),cellv3(2)
if (ierror/=0) then
    cellv1=0;cellv2=0;cellv3=0
    read(c200tmp,*) cellv1(1),cellv2(2),cellv3(3)
end if
cellv1=cellv1*10/b2a !to Angstrom, then to Bohr
cellv2=cellv2*10/b2a
cellv3=cellv3*10/b2a
ifPBC=3
close(10)

a%x=a%x*10/b2a
a%y=a%y*10/b2a
a%z=a%z*10/b2a
a%charge=a%index
iresinfo=1 !Residue information is available in this file
if (infomode==0) write(*,"(' Totally',i8,' atoms')") ncenter
call guessnelec
end subroutine



!!--------------------------------------------------------
!!-------------------- Read .gjf file --------------------
! Only support Cartesian coordinate currently
! infomode=0: Output summary, =1: Do not
subroutine readgjf(name,infomode) 
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character c40tmp*40,c200tmp*200
ifiletype=12
open(10,file=name,status="old")
iskip=0 !How many lines before charge and multiplicity line
iblank=0
do while(.true.)
	read(10,"(a)") c200tmp
	iskip=iskip+1
	if (c200tmp==" ") iblank=iblank+1
	if (iblank==2) exit
end do
read(10,*) loadcharge,loadmulti

ncenter=0
do while(.true.)
	read(10,"(a)",iostat=ierror) c200tmp
	if (c200tmp==" ".or.index(c200tmp,"Tv")/=0.or.index(c200tmp,"TV")/=0.or.index(c200tmp,"tv")/=0.or.ierror/=0) exit
	ncenter=ncenter+1
end do

if (allocated(a)) deallocate(a)
allocate(a(ncenter))
rewind(10)
do i=1,iskip
	read(10,*)
end do
read(10,*)

do i=1,ncenter
    read(10,"(a)") c200tmp
    ibeg=index(c200tmp,'(')
    iend=index(c200tmp,')')
    if (iend>ibeg) c200tmp(ibeg:iend)=" " !Purify the line to clean (...) information
	read(c200tmp,*,iostat=ierror) a(i)%name,c40tmp
	if (ierror/=0) then
		write(*,"(a)") " Error: Unable to load atom information. The input file is too complicated &
        or Z-matrix is used, Multiwfn does not support these cases"
		write(*,*) "Press ENTER button to exit program"
		read(*,*)
		stop
	end if
    if (index(c40tmp,'.')/=0) then !Didn't use 0 and -1 after atom name to indicate freeze status
	    read(c200tmp,*) a(i)%name,a(i)%x,a(i)%y,a(i)%z
    else
        read(c200tmp,*) a(i)%name,inouse,a(i)%x,a(i)%y,a(i)%z
    end if
	if (iachar(a(i)%name(1:1))>=48.and.iachar(a(i)%name(1:1))<=57) then !Use element index
		read(a(i)%name,*) a(i)%index
	else
		call lc2uc(a(i)%name(1:1)) !Convert to upper case
		call uc2lc(a(i)%name(2:2)) !Convert to lower case
        if (a(i)%name(1:2)=='X ') then !Dummy atom will be finally recognized as Bq
            a(i)%index=0
        else
			call elename2idx(a(i)%name,a(i)%index)
        end if
	end if
end do
!Try to load translation vector (Tv)
read(10,"(a)",iostat=ierror) c200tmp
if (ierror==0.and.(index(c200tmp,"Tv")/=0.or.index(c200tmp,"TV")/=0.or.index(c200tmp,"tv")/=0)) then
    ifPBC=1
    read(c200tmp,*) c40tmp,cellv1
    cellv1=cellv1/b2a
    read(10,"(a)",iostat=ierror) c200tmp
    if (ierror==0.and.(index(c200tmp,"Tv")/=0.or.index(c200tmp,"TV")/=0.or.index(c200tmp,"tv")/=0)) then
        ifPBC=ifPBC+1
        read(c200tmp,*) c40tmp,cellv2
        cellv2=cellv2/b2a
        read(10,"(a)",iostat=ierror) c200tmp
        if (ierror==0.and.(index(c200tmp,"Tv")/=0.or.index(c200tmp,"TV")/=0.or.index(c200tmp,"tv")/=0)) then
            ifPBC=ifPBC+1
            read(c200tmp,*) c40tmp,cellv3
            cellv3=cellv3/b2a
        end if
    end if
end if
close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
a%name=ind2name(a%index)
nelec=sum(a%index)-loadcharge
naelec=(nint(nelec)+loadmulti-1)/2
nbelec=nelec-naelec
if (infomode==0) then
	write(*,"(' Totally',i8,' atoms')") ncenter
	write(*,"(' The number of alpha and beta electrons:',2i8)") nint(naelec),nint(nbelec)
end if
end subroutine



!!--------------------------------------------------------
!!-------------------- Read ORCA .inp file ---------------
! Only support Cartesian coordinate currently
! infomode=0: Output summary, =1: Do not
subroutine readORCAinp(name,infomode) 
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character c40tmp*40,c200tmp*200
ifiletype=12
open(10,file=name,status="old")
call loclabel(10,"* xyz ",ifound)
read(10,*) c40tmp,c40tmp,loadcharge,loadmulti

!Test how many centers
ncenter=0
do while(.true.)
	read(10,"(a)",iostat=ierror) c200tmp
	if (index(c200tmp,'*')/=0.or.index(c200tmp,'end')/=0) exit
	if (index(c200tmp,':')/=0) cycle !Ignore ghost atom
	ncenter=ncenter+1
end do
if (allocated(a)) deallocate(a)
allocate(a(ncenter))

!Load atom information
call loclabel(10,"* xyz ",ifound)
read(10,*)
do i=1,ncenter
    read(10,"(a)") c200tmp
    ibeg=index(c200tmp,'(')
    iend=index(c200tmp,')')
    if (iend>ibeg) c200tmp(ibeg:iend)=" " !Purify the line to clean (...) information
	if (index(c200tmp,':')/=0) cycle !Ignore ghost atom
    read(c200tmp,*,iostat=ierror) a(i)%name,a(i)%x,a(i)%y,a(i)%z
	if (ierror/=0) then
		write(*,"(a)") " Error: Unable to successfully load atom information. The input file may be too complicated"
		write(*,*) "Press ENTER button to exit program"
		read(*,*)
		stop
	end if
    if (a(i)%name(1:2)=='DA') then !Dummy atom will be finally recognized as Bq
        a(i)%index=0
        cycle
    end if
    call elename2idx(a(i)%name,a(i)%index)
end do
close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
a%name=ind2name(a%index)
nelec=sum(a%index)-loadcharge
naelec=(nint(nelec)+loadmulti-1)/2
nbelec=nelec-naelec
if (infomode==0) then
	write(*,"(' Totally',i8,' atoms')") ncenter
	write(*,"(' Number of alpha and beta electrons:',2i8)") nint(naelec),nint(nbelec)
end if
end subroutine



!!--------------------------------------------------------
!!------------------- Read MOPAC .mop file ---------------
! Only support Cartesian coordinate currently
! infomode=0: Output summary, =1: Do not
subroutine readmop(name,infomode) 
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character c40tmp*40,c200tmp*200

ifiletype=12
open(10,file=name,status="old")
read(10,"(a)") c200tmp

!Load net charge
loadcharge=0
i1=index(c200tmp,"charge")
if (i1/=0) then
	read(c200tmp(i1+7:),*) loadcharge
else
	i2=index(c200tmp,"CHARGE")
    if (i2/=0) read(c200tmp(i2+7:),*) loadcharge
end if

!Load spin multiplicity
loadmulti=1
if (index(c200tmp,"doublet")/=0.or.index(c200tmp,"DOUBLET")/=0) loadmulti=2
if (index(c200tmp,"triplet")/=0.or.index(c200tmp,"TRIPLET")/=0) loadmulti=3
if (index(c200tmp,"quartet")/=0.or.index(c200tmp,"QUARTET")/=0) loadmulti=4
if (index(c200tmp,"quintet")/=0.or.index(c200tmp,"QUINTET")/=0) loadmulti=5
if (index(c200tmp,"sextet")/=0 .or.index(c200tmp,"SEXTET")/=0 ) loadmulti=6

!Test how many centers
call skiplines(10,2)
ncenter=0
do while(.true.)
	read(10,"(a)",iostat=ierror) c200tmp
	if (c200tmp==" ".or.ierror/=0) exit
	ncenter=ncenter+1
end do
if (allocated(a)) deallocate(a)
allocate(a(ncenter))

!Load atom information
rewind(10)
call skiplines(10,3)
do i=1,ncenter
    read(10,"(a)") c200tmp
    read(c200tmp,*,iostat=ierror) a(i)%name,a(i)%x,i1,a(i)%y,i2,a(i)%z,i3 !The line may contain integer for controlling optimization
    if (ierror/=0) read(c200tmp,*,iostat=ierror) a(i)%name,a(i)%x,a(i)%y,a(i)%z
	if (ierror/=0) then
		write(*,"(a)") " Error: Unable to successfully load atom information. The input file may be too non-standard"
		write(*,*) "Press ENTER button to exit program"
		read(*,*)
		stop
	end if
    call elename2idx(a(i)%name,a(i)%index)
end do
close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
a%name=ind2name(a%index)
nelec=sum(a%index)-loadcharge
naelec=(nint(nelec)+loadmulti-1)/2
nbelec=nelec-naelec
if (infomode==0) then
	write(*,"(' Number of alpha and beta electrons:',2i8)") nint(naelec),nint(nbelec)
end if
end subroutine



!!------------------------------------------------------------------------------
!!-------------------- Read CP2K input file or restart file --------------------
subroutine readcp2k(name) 
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character c40tmp*40,c200tmp*200
real*8 vec1(3),vec2(3)

ifiletype=16
open(10,file=name,status="old")

call loclabel(10,"&SUBSYS",ifound) !Will avoid to locate to &CELL in &PRINT
if (ifound==0) then !This is the simplest .inp file, only &COORD and &CELL fields
	call loclabel(10,"&CELL",ifound)
else
	call loclabel(10,"&CELL ",ifound,0) !Must be "&CELL " rather than "&CELL", the latter may locate to e.g. &CELL_OPT
end if

do while(.true.)
    read(10,"(a)") c200tmp
    if (index(c200tmp,'&CELL_REF')/=0) then !Skip &CELL_REF part
        do while(.true.)
            read(10,"(a)") c200tmp
            if (index(c200tmp,'&END')/=0) then
                read(10,*) c200tmp
                exit
            end if
        end do
    end if
    if (index(c200tmp,'&END')/=0) exit
	itmp=index(c200tmp,'[') !There may be [angstrom], remove it
	jtmp=index(c200tmp,']')
    if (itmp/=0) c200tmp(itmp:jtmp)=" "
    if (index(trim(c200tmp),'A ')/=0) read(c200tmp,*) c40tmp,cellv1
    if (index(trim(c200tmp),'B ')/=0) read(c200tmp,*) c40tmp,cellv2
    if (index(trim(c200tmp),'C ')/=0.and.index(c200tmp,'PERIODIC')==0) read(c200tmp,*) c40tmp,cellv3
    if (index(c200tmp,'ABC')/=0) then
        read(c200tmp,*) c40tmp,asize,bsize,csize
        alpha=90;beta=90;gamma=90
        call loclabel(10,"ALPHA_BETA_GAMMA",ifound)
        if (ifound==1) read(10,*) c40tmp,alpha,beta,gamma
        call abc2cellv(asize,bsize,csize,alpha,beta,gamma)
        exit
    end if
end do
cellv1=cellv1/b2a
cellv2=cellv2/b2a
cellv3=cellv3/b2a
ifPBC=3

call loclabel(10,"&SUBSYS",ifound)
if (ifound==1) then !Load &COORD within &SUBSYS
    call loclabel(10,"&COORD",ifound,0)
else !Locate to &COORD (may be a native CP2K input file, only containing geometry and cell information)
    call loclabel(10,"&COORD",ifound,1)
end if
if (ifound==0) then
    write(*,*) "Error: Unable to find &COORD field, this file cannot be loaded!"
    write(*,*) "Press ENTER button to exit"
    read(*,*)
    stop
end if
read(10,*)
ncenter=0
do while(.true.)
    read(10,"(a)") c200tmp
    if (index(c200tmp,'&END')/=0) exit
    read(c200tmp,*,iostat=ierror) c40tmp,tmpx,tmpy,tmpz
    if (ierror/=0) cycle !This line may be e.g. "UNIT angstrom" or "SCALED  F"
    ncenter=ncenter+1
end do
if (allocated(a)) deallocate(a)
allocate(a(ncenter))

call loclabel(10,"&SUBSYS",ifound)
if (ifound==1) then !Load &COORD within &SUBSYS
    call loclabel(10,"&COORD",ifound,0)
else !Locate to &COORD (may be a native CP2K input file, only containing geometry and cell information)
    call loclabel(10,"&COORD",ifound,1)
end if
read(10,*)
iscale=0
do i=1,ncenter
    read(10,"(a)") c200tmp
    if (index(c200tmp,"scaled")/=0.or.index(c200tmp,"SCALED")/=0) then
        if (index(c200tmp,"f")==0.and.index(c200tmp,"F")==0) then !Not set to false
            iscale=1
            read(10,"(a)") c200tmp
        end if
    end if
	read(c200tmp,*,iostat=ierror) a(i)%name,a(i)%x,a(i)%y,a(i)%z
    itmp=index(a(i)%name,'_') !User may use e.g. Fe_1 to indicate special atom kind
    if (itmp/=0) a(i)%name(itmp:)=" "
    do ichar=1,len_trim(a(i)%name) !Remove possible number, e.g. Fe3 to Fe
		if (iachar(a(i)%name(ichar:ichar))>=48.and.iachar(a(i)%name(ichar:ichar))<=57) a(i)%name(ichar:ichar)=" "
    end do
    call elename2idx(a(i)%name,a(i)%index)
end do
do while(.true.) !User may put SCALED after coordinates
    read(10,"(a)") c200tmp
    if (index(c200tmp,'&')/=0) exit
    if (index(c200tmp,"scaled")/=0.or.index(c200tmp,"SCALED")/=0) then
        if (index(c200tmp,"f")==0.and.index(c200tmp,"F")==0) then !Not set to false
            iscale=1
            exit
        end if
    end if
end do
close(10)

if (iscale==1) then !Fractional coordinate to Cartesian coordinate
    do iatm=1,ncenter
        vec1(1)=a(iatm)%x
        vec1(2)=a(iatm)%y
        vec1(3)=a(iatm)%z
        call fract2Cart(vec1,vec2)
        a(iatm)%x=vec2(1)
        a(iatm)%y=vec2(2)
        a(iatm)%z=vec2(3)
    end do
else !Directly loaded Cartesian coordinate in Angstrom, convert to Bohr
    a(:)%x=a(:)%x/b2a
    a(:)%y=a(:)%y/b2a
    a(:)%z=a(:)%z/b2a
end if
    
a%charge=a%index
a%name=ind2name(a%index)
call guessnelec
end subroutine



!!--------------------------------------------------------------------------
!!-------------------- Read Quantum ESPRESSO input file --------------------
! Only support ibrav=0
subroutine readQEinp(name)
use defvar
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200
character(len=*) name

ifiletype=19
open(10,file=name,status="old")

call loclabel(10,"&SYSTEM",ifound)
if (ifound==0) call loclabel(10,"&system",ifound)
if (ifound==0) then
	write(*,*) "Error: &system field cannot be found! Press ENTER button to exit"
    read(*,*)
    stop
end if

!Load &SYSTEM
alat=0
read(10,*)
do while(.true.)
	read(10,"(a)") c200tmp
	if (index(c200tmp,'/')/=0) exit
    
    itmp=index(c200tmp,'=')
    read(c200tmp(:itmp-1),*) c80tmp
	if (c80tmp=="ibrav") then
		call readaftersign_int(c200tmp,'=',itmp)
		if (itmp/=0) then
			write(*,*) "Error: Only ibrav=0 is supported! Press ENTER button to exit"
			read(*,*)
			stop
		end if
	else if (c80tmp=="nat") then
		call readaftersign_int(c200tmp,'=',ncenter)
		if (allocated(a)) deallocate(a)
		allocate(a(ncenter))
	else if (c80tmp=="celldm(1)") then
		call readaftersign_float(c200tmp,'=',alat)
	else if (c80tmp=="a".or.c80tmp=="A") then
		call readaftersign_float(c200tmp,'=',alat)
        alat=alat/b2a
    end if
end do
if (alat/=0) write(*,"(' alat is',f12.6,' Bohr')") alat

!Load CELL_PARAMETERS
call loclabel(10,"CELL_PARAMETERS",ifound)
if (ifound==0) call loclabel(10,"cell_parameters",ifound)
if (ifound==0) then
	write(*,*) "Error: CELL_PARAMETERS card cannot be found! Press ENTER button to exit"
    read(*,*)
    stop
else
	read(10,"(a)") c200tmp
	read(10,*) cellv1(:)
	read(10,*) cellv2(:)
	read(10,*) cellv3(:)
	if (index(c200tmp,'angstrom')/=0.or.index(c200tmp,'Angstrom')/=0) then
		cellv1=cellv1/b2a
		cellv2=cellv2/b2a
		cellv3=cellv3/b2a
	else if (index(c200tmp,'alat')/=0) then
		cellv1=cellv1*alat
		cellv2=cellv2*alat
		cellv3=cellv3*alat
	end if
	ifPBC=3
end if

!Load ATOMIC_POSITIONS
call loclabel(10,"ATOMIC_POSITIONS",ifound)
if (ifound==0) call loclabel(10,"atomic_positions",ifound)
if (ifound==0) then
	write(*,*) "Error: ATOMIC_POSITIONS card cannot be found! Press ENTER button to exit"
    read(*,*)
    stop
else
	read(10,"(a)") c200tmp
	do iatm=1,ncenter
		read(10,*) c80tmp,a(iatm)%x,a(iatm)%y,a(iatm)%z
        !User may use customized name, such as Fe2. So we only keep the first two letter
        a(iatm)%name=c80tmp(1:2)
		call lc2uc(a(iatm)%name(1:1))
		call uc2lc(a(iatm)%name(2:2))
		call elename2idx(a(iatm)%name,a(iatm)%index)
	end do
	if (index(c200tmp,'angstrom')/=0.or.index(c200tmp,'Angstrom')/=0) then
		a%x=a%x/b2a
		a%y=a%y/b2a
		a%z=a%z/b2a
	else if (index(c200tmp,'alat')/=0) then
		a%x=a%x*alat
		a%y=a%y*alat
		a%z=a%z*alat
	end if
end if

close(10)

a%charge=a%index
call guessnelec
end subroutine



!!----------------------------------------
!!------------- Read POSCAR of VASP
!See https://www.vasp.at/wiki/index.php/POSCAR for format description
subroutine readPOSCAR(name,infomode) !infomode=0 means output info, =1 silent
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode,anum(200)
character(len=*) name
character titleline*200,c200tmp*200,atype(200)*2,ctest
real*8 fract(3),Cart(3)

ifiletype=18
open(10,file=name,status="old")
read(10,"(a)") titleline
read(10,*) sclfac

read(10,*) cellv1
read(10,*) cellv2
read(10,*) cellv3
cellv1=cellv1*sclfac/b2a
cellv2=cellv2*sclfac/b2a
cellv3=cellv3*sclfac/b2a
ifPBC=3

!Test how many types
read(10,"(a)") c200tmp
read(c200tmp,*,iostat=ierror) i
if (ierror/=0) then !Element name is explicitly given
    backspace(10)
	call numdatastr(c200tmp,ntype)
    read(10,*) atype(1:ntype) !Read element names
    read(10,*) anum(1:ntype) !Read number of atoms of each type
else !Element name is not explicitly given
    call numdatastr(c200tmp,ntype)
    read(c200tmp,*) anum(1:ntype)
    write(*,*) "Element name is not explicitly given in POSCAR, please manually input"
    do itype=1,ntype
        write(*,"(' Input element name for type',i3, ' ( Number of atoms is',i5,' )  e.g. Fe')") itype,anum(itype)
        read(*,*) atype(itype)
    end do
end if

read(10,*) ctest
if (ctest=="s".or.ctest=="S") then !The read line is selective dynamics, read next line for test coordinate type
    read(10,*) ctest
end if
icoord=0 !Assume to be "direct", namely fractional coordinate
if (ctest=='k'.or.ctest=='K'.or.ctest=='c'.or.ctest=='C') icoord=1 !Cartesian coordinate

ncenter=sum(anum(1:ntype))
if (allocated(a)) deallocate(a)
allocate(a(ncenter))

ncenter=0
do itype=1,ntype
    nthis=anum(itype)
    do iatm=ncenter+1,ncenter+nthis
        read(10,*) xtmp,ytmp,ztmp
        if (icoord==0) then !Fractional coordinate
            fract(1)=xtmp
            fract(2)=ytmp
            fract(3)=ztmp
            call fract2Cart(fract,Cart)
            a(iatm)%x=Cart(1)
            a(iatm)%y=Cart(2)
            a(iatm)%z=Cart(3)
        else if (icoord==1) then !Cartesian coordinate
            a(iatm)%x=xtmp/b2a*sclfac
            a(iatm)%y=ytmp/b2a*sclfac
            a(iatm)%z=ztmp/b2a*sclfac
        end if
    end do
    !Assign element index
    call elename2idx(atype(itype),idx)
    a(ncenter+1:ncenter+nthis)%index=idx
    ncenter=ncenter+nthis
end do

close(10)

a%charge=a%index
a%name=ind2name(a%index)
call guessnelec
if (infomode==0) then
    write(*,*) "Title line of this file:"
    write(*,"(a)") trim(titleline)
	write(*,"(' Totally',i8,' atoms')") ncenter
end if
end subroutine



!!----------------------------------------------------------------------------------
!!------------- Read grid data of VASP, CHGCAR/CHG/ELFCAR/LOCPOT files are supported
!CHGCAR/CHG (electron density multipled by cell volume) is outputted by VASP by default
!ELFCAR (ELF) can be generated by using LELF=.TRUE.
!LOCPOT (recorded in eV) can be generated by using LVTOT=.TRUE. If LVHAR=.TRUE., V_XC is not included, and the grid data &
!  corresponds to Hartree potential, namely negative of ESP. If LVHAR=.FALSE., the grid data corresponds to external potential in KS operator
!In spin polarized case (ISPIN=2), there are two sets of data:
!  CHGCAR/CHG records total density and spin density, ELFCAR records ELF for alpha and beta, LOCPOT records potential for alpha and beta
!See https://www.vasp.at/wiki/index.php/CHGCAR for format description
!Also see https://www.vasp.at/wiki/index.php/ELFCAR and https://www.vasp.at/wiki/index.php/LOCPOT
subroutine readVASPgrd(name,infomode) !infomode=0 means output info, =1 silent
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character c200tmp*200,c80tmp*80,selectyn
real*8,allocatable :: tmpeffchg(:)
character,allocatable :: tmpelestr(:)*2

!The first part of CHGCAR/CHG is identical to POSCAR, load atom and cell information from this part
call readPOSCAR(name,0)

ifiletype=7 !Same as cube

open(10,file=name,status="old")
do while(.true.)
    read(10,"(a)") c200tmp
	if (index(c200tmp,"Nval")/=0) then
		call numdatastr(c200tmp(6:),ndata)
		nread=ndata/2
		allocate(tmpeffchg(nread),tmpelestr(nread))
		read(c200tmp(6:),*) (tmpelestr(i),tmpeffchg(i),i=1,nread)
		if (nread>0) then
			write(*,*) "Loaded effective nuclear charges"
			do i=1,nread
				write(*,"(1x,a,f12.6)") tmpelestr(i),tmpeffchg(i)
				call lc2uc(tmpelestr(i)(1:1))
				call uc2lc(tmpelestr(i)(2:2))
				do iatm=1,ncenter
					if (a(iatm)%name==tmpelestr(i)) a(iatm)%charge=tmpeffchg(i)
				end do
			end do
		end if
	else if (c200tmp==" ") then !Find the first blank line, the following data are grid data
		exit
    end if
end do

read(10,*) nx,ny,nz
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
write(*,*) "Loading grid data..."
read(10,*) (((cubmat(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)

!Check if there is second grid data
write(c80tmp,"(3i5)") nx,ny,nz
call loclabel(10,trim(c80tmp),ifound,0)
if (ifound==1) then
    write(*,*)
    if (index(name,"CHG")/=0) then
        write(*,"(a)") " This file also contains spin density, do you want to load it instead of the loaded data (total electron density)? (y/n)"
    else if (index(name,"ELFCAR")/=0) then
        write(*,"(a)") " This file also contains ELF of beta electrons, do you want to load it instead of the loaded data (ELF of alpha electrons)? (y/n)"
    else if (index(name,"LOCPOT")/=0) then
        write(*,"(a)") " This file also contains potential for beta electrons, do you want to load it instead of potential for alpha electrons? (y/n)"
    end if
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') read(10,*) (((cubmat(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
end if

if (index(name,"CHG")/=0) then
    write(*,"(a)") " The CHGCAR/CHG file originally records rho*V_cell, do you want to convert it to electron density (rho, in a.u.)? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
        call calc_cellvol(cellvol)
        cubmat=cubmat/cellvol
    end if
else if (index(name,"LOCPOT")/=0) then
    write(*,"(/,a)") " The LOCPOT file originally records potential felt by single electron, &
    do you want to convert it to electrostatic potential (namely taking negative of that)? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') cubmat=-cubmat
    write(*,"(/,a)") " The grid data in LOCPOT file was originally recorded in eV, do you want to convert it to a.u. to follow common convention? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') cubmat=cubmat/au2eV
end if
close(10)

orgx=0
orgy=0
orgz=0
gridv1=cellv1/nx
gridv2=cellv2/ny
gridv3=cellv3/nz
dx=gridv1(1)
dy=gridv2(2)
dz=gridv3(3)
call getgridend !Generate endx,endy,endz
if (infomode==0) then
    call showgridinfo(1) !Calculate grid information
    call showgridinfo(2) !Calculate statistical information
end if
end subroutine


!!---------------------------------------------------------------------
!!----- The same as readVASPgrd, but only load grid data of CHGCAR/CHG, and the number of grid will be checked so that they are in agreement with cubmat
!inconsis is returned value, 1 means the grid setting of this cube file is inconsistent with that of cubmat
subroutine readVASPgrdtmp(name,inconsis)
use defvar
use util
implicit real*8 (a-h,o-z)
integer inconsis
character(len=*) name
character c200tmp*200,c80tmp*80,selectyn

!Direct locate to the first blank line and load grid data
open(10,file=name,status="old")
do while(.true.)
    read(10,"(a)") c200tmp
    if (c200tmp==" ") exit
end do
read(10,*) nxtmp,nytmp,nztmp

inconsis=0
if (nxtmp/=nx.or.nytmp/=ny.or.nztmp/=nz) then
	write(*,"(a)") " Error: The number of grid of this file is inconsistent with that of the grid data stored in memory!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
	inconsis=1
    close(10)
	return
end if

if (allocated(cubmattmp)) deallocate(cubmattmp)
allocate(cubmattmp(nx,ny,nz))
write(*,*) "Loading grid data..."
read(10,*) (((cubmattmp(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)

!Check if there is second grid data
write(c80tmp,"(3i5)") nx,ny,nz
call loclabel(10,trim(c80tmp),ifound,0)
if (ifound==1) then
    write(*,*)
    write(*,"(a)") " This file also contains spin density, do you want to load it instead of loaded data (total density)? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') read(10,*) (((cubmattmp(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
end if

if (index(name,"CHG")/=0) then
    write(*,"(a)") " The CHGCAR/CHG file originally records rho*V_cell, do you want to convert it to rho? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
        call calc_cellvol(cellvol)
        cubmattmp=cubmattmp/cellvol
    end if
end if
close(10)

end subroutine




!!------------------- Read MDL .mol file (V2000) -------------------
! sdf can also be load, as it is a wrapper of sdf with additional information
! Format description: https://en.wikipedia.org/wiki/Chemical_table_file
! infomode=0: Output summary, =1: do not
subroutine readmol(name,infomode) 
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name

ifiletype=11
open(10,file=name,status="old")
!Molden file containing vibrational modes exported by CP2K also has .mol extension. If it is the case, regarded as plain text file
call loclabel(10,"[Molden Format]",ifound,maxline=100)
if (ifound==1) then
	ifiletype=0
    close(10)
    return
end if
rewind(10)
read(10,"(a)") titleline
read(10,*)
read(10,*)
read(10,"(2i3)") ncenter,nbond
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
do i=1,ncenter
	read(10,*) a(i)%x,a(i)%y,a(i)%z,a(i)%name
    call elename2idx(a(i)%name,a(i)%index)
end do
if (allocated(connmat)) deallocate(connmat)
allocate(connmat(ncenter,ncenter))
connmat=0
do ibond=1,nbond
	read(10,"(3i3)") i,j,ntmp
	connmat(i,j)=ntmp
	connmat(j,i)=ntmp
end do
close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
if (infomode==0) write(*,"(' Totally',i8,' atoms')") ncenter
call guessnelec
end subroutine




!!------------------- Read .mol2 file -------------------
! infomode=0: Output summary, =1: do not
subroutine readmol2(name,infomode) 
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character c80tmp*80
ifiletype=13
open(10,file=name,status="old")
call loclabel(10,"@<TRIPOS>MOLECULE")
read(10,*)
read(10,"(a)") titleline
read(10,*) ncenter,nbond

if (allocated(a)) deallocate(a)
allocate(a(ncenter))
call loclabel(10,"@<TRIPOS>ATOM")
read(10,*)
do i=1,ncenter
	read(10,*) inouse,c80tmp,a(i)%x,a(i)%y,a(i)%z,a(i)%name
	ico=index(a(i)%name,'.')
	if (ico/=0) a(i)%name(ico:)=" "
    call elename2idx(a(i)%name,a(i)%index)
end do

if (allocated(connmat)) deallocate(connmat)
allocate(connmat(ncenter,ncenter))
call loclabel(10,"@<TRIPOS>BOND")
read(10,*)
connmat=0
!gview records aromatic carbon as type of 4 in .mol, but record it as Ar
!When loading such kind of atom, Multiwfn still recognize it as type 4
do ibond=1,nbond
	read(10,*) inouse,i,j,c80tmp
	if (c80tmp=="ar".or.c80tmp=="Ar") then
		ntmp=4
	else if (c80tmp=="am") then
		ntmp=1
	else if (c80tmp=="un".or.c80tmp=="nc".or.c80tmp=="du") then
		ntmp=0
	else
		read(c80tmp,*) ntmp
	end if
	connmat(i,j)=ntmp
	connmat(j,i)=ntmp
end do

call loclabel(10,"@<TRIPOS>CRYSIN",ifound)
if (ifound==1) then
    ifPBC=3
    read(10,*)
    read(10,*) asize,bsize,csize,alpha,beta,gamma
    asize=asize/b2a
    bsize=bsize/b2a
    csize=csize/b2a
    call abc2cellv(asize,bsize,csize,alpha,beta,gamma)
end if

close(10)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
a%charge=a%index
if (infomode==0) write(*,"(' Totally',i8,' atoms')") ncenter
call guessnelec
end subroutine




!------------- Only read connectivity from .mol or .mol2 file
subroutine readmolconn(name) 
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character c80tmp*80

if (allocated(connmat)) deallocate(connmat)
allocate(connmat(ncenter,ncenter))
connmat=0
    
open(10,file=name,status="old")
if (index(name,".mol2")/=0) then !Read mol2
	call loclabel(10,"@<TRIPOS>MOLECULE")
    read(10,*)
    read(10,*)
    read(10,*) ntmp,nbond
    if (ntmp/=ncenter) then
		write(*,"(a)") " Warning: The number of atoms in the .mol2 file is not the same as present system!"
		!write(*,"(a)") " Press ENTER button to continue (Multiwfn may crash)"
  !      read(10,*)
  !      pause
    end if
	call loclabel(10,"@<TRIPOS>BOND")
	read(10,*)
	do ibond=1,nbond
		read(10,*) inouse,i,j,c80tmp
		if (c80tmp=="ar".or.c80tmp=="Ar") then
			ntmp=4
		else if (c80tmp=="am") then
			ntmp=1
		else if (c80tmp=="un".or.c80tmp=="nc".or.c80tmp=="du") then
			ntmp=0
		else
			read(c80tmp,*) ntmp
		end if
		connmat(i,j)=ntmp
		connmat(j,i)=ntmp
	end do
else !Read mol
	read(10,*)
	read(10,*)
	read(10,*)
	read(10,"(2i3)") ntmp,nbond
    if (ntmp/=ncenter) then
		write(*,"(a)") " Warning: The number of atoms in the .mol2 file is not the same as present system!"
		!write(*,"(a)") " Press ENTER button to continue (Multiwfn may crash)"
  !      read(10,*)
  !      pause
    end if
	do i=1,ntmp
		read(10,*)
	end do
	do ibond=1,nbond
		read(10,"(3i3)") i,j,ntmp
		connmat(i,j)=ntmp
		connmat(j,i)=ntmp
	end do
end if
close(10)
end subroutine




!!-----------------------------------------------------------------
!! Read .31 and one of .32 to .40 file generated by NBO program
!! After loading, only GTF information is yielded, while basis functions are discarded
subroutine read31(name,infomode) !mode=0 means output wfn property,=1 not
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character :: name2*200=" ",chartemp*80,tmpc2*2
integer infomode
integer bastype2func(500) !Convert basis type in NBO (NBO5.0 Program Manual P103) to function type in .wfn
integer,allocatable :: shellcon(:),shell2atom(:),shellnumbas(:),shell2prmshl(:),bastype31(:),shellnumbas5D(:)
integer bastype_p(3) !The NBO plot file outputted by ORCA+NBO use special order of p shell, therefore record actual order and use later
real*8,allocatable :: orbcoeff(:,:),prmshlexp(:),cs(:),cp(:),cd(:),cf(:),cg(:),orbcoeff5D(:,:)
real*8 :: conv5d6d(6,5)=0D0,conv7f10f(10,7)=0D0,conv9g15g(15,9)=0D0
real*8,external :: normgau
ifiletype=6
bastype2func(1)=1 !s
bastype2func(101)=2 !x
bastype2func(102)=3 !y
bastype2func(103)=4 !z
bastype2func(201)=5 !xx
bastype2func(202)=8 !xy
bastype2func(203)=9 !xz
bastype2func(204)=6 !yy
bastype2func(205)=10 !yz
bastype2func(206)=7 !zz
bastype2func(301)=11 !xxx
bastype2func(302)=14 !xxy
bastype2func(303)=15 !xxz
bastype2func(304)=17 !xyy
bastype2func(305)=20 !xyz
bastype2func(306)=18 !xzz
bastype2func(307)=12 !yyy
bastype2func(308)=16 !yyz
bastype2func(309)=19 !yzz
bastype2func(310)=13 !zzz
!Below g sequence comes from line 47384 in NBO_5 src
bastype2func(401)=35 !XXXX
bastype2func(402)=34 !XXXY
bastype2func(403)=33 !XXXZ
bastype2func(404)=32 !XXYY
bastype2func(405)=31 !XXYZ
bastype2func(406)=30 !XXZZ
bastype2func(407)=29 !XYYY
bastype2func(408)=28 !XYYZ
bastype2func(409)=27 !XYZZ
bastype2func(410)=26 !XZZZ
bastype2func(411)=25 !YYYY
bastype2func(412)=24 !YYYZ
bastype2func(413)=23 !YYZZ
bastype2func(414)=22 !YZZZ
bastype2func(415)=21 !ZZZZ
!Conversion relationship used below can be found in http://sobereva.com/97
!Used to convert coefficient matrix from 5D to 6D
!5D sequence in .31: 255   252   253   254   251
!namely -0.5*XX-0.5*YY+ZZ, XZ, YZ, √3/2*(XX-YY), XY
!to 6D: XX,XY,XZ,YY,YZ,ZZ, namely 201~206, the indexes are consecutive
conv5d6d(1,1)=-0.5D0
conv5d6d(4,1)=-0.5D0
conv5d6d(6,1)=1D0
conv5d6d(3,2)=1D0
conv5d6d(5,3)=1D0
conv5d6d(1,4)=sqrt(3D0)/2D0
conv5d6d(4,4)=-sqrt(3D0)/2D0
conv5d6d(2,5)=1D0
!Used to convert coefficient matrix from 7F to 10F, Standard f set
!7F sequence in .31: 351   352   353   354   355   356   357
!to 10F: XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ, namely 301~310, the indexes are consecutive
conv7f10f(3,1)=-0.721962098225322D0
conv7f10f(8,1)=-0.721962098225322D0
conv7f10f(10,1)=0.481308065483548D0
conv7f10f(1,2)=-0.281160203343101D0
conv7f10f(4,2)=-0.281160203343101D0
conv7f10f(6,2)=1.1246408133724D0
conv7f10f(2,3)=-0.281160203343101D0
conv7f10f(7,3)=-0.281160203343101D0
conv7f10f(9,3)=1.1246408133724D0
conv7f10f(3,4)=sqrt(3D0)/2D0
conv7f10f(8,4)=-sqrt(3D0)/2D0
conv7f10f(5,5)=1D0
conv7f10f(1,6)=0.369693511996758D0
conv7f10f(4,6)=-1.10908053599027D0
conv7f10f(2,7)=1.10908053599027D0
conv7f10f(7,7)=-0.369693511996758D0

write(*,*) "Input filename with suffix ranging from .32 to .40 (e.g. ltwd.35)"
write(*,*) "Note: 32=PNAO 33=NAO 34=PNHO 35=NHO 36=PNBO 37=NBO 38=PNLMO 39=NLMO 40=MO"
do while(.true.)
	read(*,"(a)") name2
	if (name2(1:2)=="32".or.name2(1:2)=="33".or.name2(1:2)=="34".or.name2(1:2)=="35"&
	.or.name2(1:2)=="36".or.name2(1:2)=="37".or.name2(1:2)=="38".or.name2(1:2)=="39".or.name2(1:2)=="40") then
		tmpc2=name2(1:2)
		name2(1:len(name))=name
		name2(len_trim(name2)-1:len_trim(name2))=tmpc2
	end if
	inquire(file=name2,exist=alive)
	if (alive.eqv..true.) exit
	write(*,*) "File not found, input again"
end do

itmplen=len_trim(name2)
if (name2(itmplen-1:itmplen)=="32") write(*,*) "Loading .32 file(PNAO)"
if (name2(itmplen-1:itmplen)=="33") write(*,*) "Loading .33 file(NAO)"
if (name2(itmplen-1:itmplen)=="34") write(*,*) "Loading .34 file(PNHO)"
if (name2(itmplen-1:itmplen)=="35") write(*,*) "Loading .35 file(NHO)"
if (name2(itmplen-1:itmplen)=="36") write(*,*) "Loading .36 file(PNBO)"
if (name2(itmplen-1:itmplen)=="37") write(*,*) "Loading .37 file(NBO)"
if (name2(itmplen-1:itmplen)=="38") write(*,*) "Loading .38 file(PNLMO)"
if (name2(itmplen-1:itmplen)=="39") write(*,*) "Loading .39 file(NLMO)"
if (name2(itmplen-1:itmplen)=="40") write(*,*) "Loading .40 file(MO)"
open(10,file=name,status="old")
read(10,*)
read(10,*)
read(10,*)
read(10,*) ncenter,nshell,nprimshell
allocate(a(ncenter),shellcon(nshell),shell2atom(nshell),shellnumbas(nshell),shell2prmshl(nshell))
allocate(prmshlexp(nprimshell),cs(nprimshell),cp(nprimshell),cd(nprimshell),cf(nprimshell),cg(nprimshell))
allocate(bastype31(nshell*15)) !We don't know how many basis before read the file, so use the maximum value(up to g function)
read(10,*)
do i=1,ncenter
	read(10,*,iostat=ierror) a(i)%index,a(i)%x,a(i)%y,a(i)%z
    if (ierror/=0) then
        write(*,"(a,i5,a)") " Error: Coordinate of atom",i," cannot be normally loaded!"
        write(*,"(a)") " Please open .31 file by text editor and check if coordinate of this atom (and may be other ones) is ******. &
        If yes, you should manually replace it with actual coordinate"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        stop
    end if
	a(i)%charge=a(i)%index !.31 doesn't record charge and index separately, so we have to make them indentical
end do
a%name=ind2name(a%index)
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
read(10,*)
j=1
do i=1,nshell
	read(10,*) shell2atom(i),shellnumbas(i),shell2prmshl(i),shellcon(i)
	read(10,*) bastype31(j:j+shellnumbas(i)-1)
    if (shellnumbas(i)==3) bastype_p(1:3)=bastype31(j:j+shellnumbas(i)-1) !Record actual order of p shell
	j=j+shellnumbas(i)
end do
isphergau=0
if (any(bastype31==251).or.any(bastype31==351).or.any(bastype31==451)) isphergau=1
!The conversion relationship between Cartesian and pure functions used by NBO is different to mainstream quantum chemistry code, the relationship for f is documented
!in the NBO manual, however there is no way to find out that for g-type. So if g-type is involved, Cartesian type must be used.
if (any(bastype31==451)) then
	write(*,"(a)") " Error: Multiwfn does not support spherical harmonic Gaussian functions with g or higher angular moment in NBO plot files. &
	&If you used Gaussian to generate them, you should add ""6d 10f"" keywords and regenerate these files"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	stop
end if
read(10,*)
read(10,*) prmshlexp
read(10,*)
read(10,*) cs
read(10,*)
read(10,*) cp
read(10,*)
read(10,*) cd
read(10,*)
read(10,*) cf
read(10,"(a)",iostat=ierror) chartemp
if (ierror==0) then
	backspace(10)
	read(10,*) cg
end if
close(10)

totenergy=0D0
virialratio=0D0
nbasis=sum(shellnumbas) !Calculate the number of basis
write(*,"(' Expected number of basis functions:',i10)") nbasis
open(10,file=name2,status="old")
read(10,*)
read(10,*)
read(10,*)
read(10,"(a)") chartemp
naelec=0D0
nbelec=0D0
!Note: When diffuse functions are used, although the number of AO in NBO plot files are always identical to nbasis, &
!the number of NAOs and thus the resulting NBOs etc. may be smaller than nbasis (and can also be different to Nbsuse in Gaussian output file)
!"nload" is used to test how many frames (orbitals) can be normally loaded
if (chartemp(1:11)==" ALPHA SPIN".or.chartemp(1:11)==" alpha spin") then
	wfntype=4
	nmo=2*nbasis
	allocate(orbcoeff(nbasis,nmo),MOocc(nmo),MOene(nmo),MOtype(nmo))
	MOtype(1:nbasis)=1 !alpha
	MOtype(nbasis+1:nmo)=2 !beta
	orbcoeff=0
	MOocc=0
	nload=0
	do iorb=1,nbasis !Note that the occupation section may be erronesouly loaded as orbital coefficient here when linear dependency problem occur
		read(10,"(1x,5f15.9)",iostat=ierror) orbcoeff(1:nbasis,iorb)
		if (ierror/=0) exit
		nload=nload+1
	end do
	if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") then !Need to load occupation
		if (nload/=nbasis) then
			write(*,"(/,a)") " Warning: The number of orbitals is smaller than basis functions! This is often because diffuse functions are used. &
			&Please now input the actual number of NAOs (you can easily find it from output of NBO program), e.g. 374"
            write(*,"(' If you really do not know how to input, try to input',i6)") nload
			read(*,*) nNAOs
			call loclabel(10,chartemp(1:11))
			read(10,*)
			do iorb=1,nNAOs
				read(10,*) orbcoeff(1:nbasis,iorb)
			end do
			orbcoeff(1:nbasis,nNAOs+1:nbasis)=0
		else
			nNAOs=nbasis
		end if
        !write(*,*) nNAOs,nload,nbasis
		read(10,*) MOocc(1:nNAOs)
	end if
	if (chartemp(1:11)==" ALPHA SPIN") then
		call loclabel(10," BETA  SPIN",ifound)
	else
		call loclabel(10," beta  spin",ifound)
	end if
	if (ifound==1) then
		read(10,*)
		do iorb=nbasis+1,nbasis+nNAOs
			read(10,*) orbcoeff(1:nbasis,iorb)
		end do
		if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") read(10,*) MOocc(nbasis+1:nbasis+nNAOs)
	else
		write(*,*) "Warning: Beta-spin information was not found in this file!"
		write(*,*)
	end if
	naelec=sum(MOocc(1:nbasis))
	nbelec=sum(MOocc(nbasis+1:nmo))
	nelec=naelec+nbelec
else !Close shell system
	wfntype=3
	nmo=nbasis
	allocate(orbcoeff(nbasis,nmo),MOocc(nmo),MOene(nmo),MOtype(nmo))
	MOtype=0
	orbcoeff=0
	MOocc=0
	nload=0
	backspace(10)
	do iorb=1,nmo
		read(10,"(1x,5f15.9)",iostat=ierror) orbcoeff(1:nbasis,iorb)
        !if (ierror/=0) write(*,"(5f15.9)") orbcoeff(1:nbasis,iorb)
		if (ierror/=0) exit
		nload=nload+1
	end do
    !write(*,*) nload,nbasis
	if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") then
		if (nload/=nbasis) then
			write(*,"(/,a)") " Warning: The number of orbitals is smaller than basis functions! This is often because diffuse functions are used. &
			&Please now input the actual number of NAOs (you can easily find it from output of NBO program), e.g. 374"
            write(*,"(' If you really do not know how to input, try to input',i6)") nload
			read(*,*) nNAOs
            rewind(10)
            read(10,*)
            read(10,*)
            read(10,*)
			do iorb=1,nNAOs
				read(10,*) orbcoeff(1:nbasis,iorb)
			end do
			orbcoeff(1:nbasis,nNAOs+1:nbasis)=0
		else
			nNAOs=nbasis
		end if
		read(10,*) MOocc(1:nNAOs)
	end if
	nelec=sum(MOocc)
	naelec=nelec/2
	nbelec=naelec
end if
close(10)
MOene=0D0

!Temporarily convert spherical harmonic Gaussian functions' information to Cartesian type
if (isphergau==1) then
	!Calculate how many Cartesian basis after conversion
	nbasis5D=nbasis
	nbasis=0
	do i=1,nshell
		if (shellnumbas(i)==5) then !D
			nbasis=nbasis+6
		else if (shellnumbas(i)==7) then !F
			nbasis=nbasis+10
		else if (shellnumbas(i)==9) then !G
			nbasis=nbasis+15
		else !S,P,SP, or Cartesian shells
			nbasis=nbasis+shellnumbas(i)
		end if
	end do
	allocate(shellnumbas5D(nshell))
	shellnumbas5D=shellnumbas !Backup
	where(shellnumbas==5) shellnumbas=6 !Convert number of orbitals in each shell from 5D to Cartesian type
	where(shellnumbas==7) shellnumbas=10
	where(shellnumbas==9) shellnumbas=15

	allocate(orbcoeff5D(nbasis5D,nmo))
	orbcoeff5D=orbcoeff !Backup
	deallocate(orbcoeff)
	allocate(orbcoeff(nbasis,nmo)) !Enlarge size from spherical type to Cartesian type
	orbcoeff=0D0

	deallocate(bastype31)
	allocate(bastype31(nbasis)) !Enlarge size from spherical type to Cartesian type

	!Generate Cartesian .31 basis type, the indexes are consecutive, in line with conv5d6d and conv7f10f
	i=1
	do ish=1,nshell
		if (shellnumbas(ish)==1) then !s
			bastype31(i)=1
		else if (shellnumbas(ish)==3) then !p
            bastype31(i:i+2)=bastype_p(:)
		else if (shellnumbas(ish)==4) then !sp
			bastype31(i)=1
			bastype31(i+1)=101
			bastype31(i+2)=102
			bastype31(i+3)=103
		else if (shellnumbas(ish)==6) then !d
			do j=1,6
				bastype31(i+j-1)=200+j
			end do
		else if (shellnumbas(ish)==10) then !f
			do j=1,10
				bastype31(i+j-1)=300+j
			end do
		else if (shellnumbas(ish)==15) then !g
			do j=1,15
				bastype31(i+j-1)=400+j
			end do
		end if
		i=i+shellnumbas(ish)
	end do
	
	!Map 5D coefficient to 6D coefficient
	ipos5D=1
	ipos6D=1
	do ish=1,nshell
		n5D=shellnumbas5D(ish)
		n6D=shellnumbas(ish)
		if (n5D==1.or.n5D==3.or.n5D==4) then !S or P or SP
			if (wfntype==3) orbcoeff(ipos6D:ipos6D+n6D-1,1:nbasis5D)=orbcoeff5D(ipos5D:ipos5D+n5D-1,:)
			if (wfntype==4) orbcoeff(ipos6D:ipos6D+n6D-1,1:2*nbasis5D)=orbcoeff5D(ipos5D:ipos5D+n5D-1,:)
		else if (n5D==5) then
			!5D->6D
			if (wfntype==3) orbcoeff(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,orbcoeff5D(ipos5D:ipos5D+4,:))
			if (wfntype==4) orbcoeff(ipos6D:ipos6D+5,1:2*nbasis5D)=matmul(conv5d6d,orbcoeff5D(ipos5D:ipos5D+4,:))
		else if (n5D==7) then
			!7F->10F
			if (wfntype==3) orbcoeff(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,orbcoeff5D(ipos5D:ipos5D+6,:))
			if (wfntype==4) orbcoeff(ipos6D:ipos6D+9,1:2*nbasis5D)=matmul(conv7f10f,orbcoeff5D(ipos5D:ipos5D+6,:))
		else if (n5D==9) then
			!9G->15G
			if (wfntype==3) orbcoeff(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,orbcoeff5D(ipos5D:ipos5D+8,:))
			if (wfntype==4) orbcoeff(ipos6D:ipos6D+14,1:2*nbasis5D)=matmul(conv9g15g,orbcoeff5D(ipos5D:ipos5D+8,:))
		end if
		ipos5D=ipos5D+n5D
		ipos6D=ipos6D+n6D
	end do
end if

nprims=0
do i=1,nshell
	nprims=nprims+shellcon(i)*shellnumbas(i)
end do
allocate(b(nprims),CO(nmo,nprims))

iGTF=1 !current GTF index
ibasis=1
do i=1,nshell !cycle each shell
	b(iGTF:iGTF+shellcon(i)*shellnumbas(i)-1)%center=shell2atom(i)
	do j=1,shellnumbas(i) !cycle each basis function in each shell
		b(iGTF:iGTF+shellcon(i)-1)%type=bastype2func(bastype31(ibasis))
		do k=1,shellcon(i) !cycle each GTF in each basis in each shell
			iprmshlpos=shell2prmshl(i)+k-1
			b(iGTF)%exp=prmshlexp(iprmshlpos)
			if (bastype31(ibasis)==1) then !s
				contract=cs(iprmshlpos)
			else if (bastype31(ibasis)<=200) then !p
				contract=cp(iprmshlpos)
			else if (bastype31(ibasis)<=300) then !d
			!Contract coefficient in .31 contains normalization coefficient, however for d type, the
			!normalization coefficient is for XX,YY,ZZ, for XY,XZ,YZ, we need refresh normalization coefficient
				contract=cd(iprmshlpos)
				if (bastype31(ibasis)==202.or.bastype31(ibasis)==203.or.bastype31(ibasis)==205) then
					valnorm31=normgau(5,prmshlexp(iprmshlpos)) !Normalization coefficient for XX,YY,ZZ are identical
					valnormnew=normgau(8,prmshlexp(iprmshlpos)) !Normalization coefficient for XY,XZ,YZ are identical
					contract=contract/valnorm31*valnormnew
				end if
			else if (bastype31(ibasis)<=400) then !f
				contract=cf(iprmshlpos)
				!For f shell, in .31 normalization coefficient is for XXX,YYY,ZZZ, now refresh
				if (bastype31(ibasis)/=301.and.bastype31(ibasis)/=307.and.bastype31(ibasis)/=310) then !not XXX,YYY,ZZZ
					valnorm31=normgau(11,prmshlexp(iprmshlpos))
					if (bastype31(ibasis)==302.or.bastype31(ibasis)==303.or.bastype31(ibasis)==304&
					.or.bastype31(ibasis)==306.or.bastype31(ibasis)==308.or.bastype31(ibasis)==309) then
						valnormnew=normgau(14,prmshlexp(iprmshlpos)) !XXY,XXZ,XYY,XZZ,YYZ,YZZ are identical
					else if (bastype31(ibasis)==305) then 
						valnormnew=normgau(20,prmshlexp(iprmshlpos)) !XYZ
					end if
					contract=contract/valnorm31*valnormnew
				end if
			else if (bastype31(ibasis)<=500) then !g
				contract=cg(iprmshlpos)
				!For g shell, in .31 normalization coefficient is for XXXX,YYYY,ZZZZ, now refresh
				!Note: I haven't verified that, since .37 outputted by NBO3.1 in Gaussian didn't contains g information
				nt=bastype31(ibasis) !now type
				if (nt/=401.and.nt/=411.and.nt/=415) then !not XXXX,YYYY,ZZZZ (4,0)
					valnorm31=normgau(21,prmshlexp(iprmshlpos))
					if (nt==402.or.nt==403.or.nt==407.or.nt==410.or.nt==412.or.nt==414) then
						valnormnew=normgau(22,prmshlexp(iprmshlpos)) !XXXY,XXXZ,XYYY,XZZZ,YYYZ,YZZZ (3,1)
					else if (nt==404.or.nt==406.or.nt==413) then
						valnormnew=normgau(23,prmshlexp(iprmshlpos)) !XXYY,XXZZ,YYZZ (2,2)
					else if (nt==405.or.nt==408.or.nt==409) then
						valnormnew=normgau(27,prmshlexp(iprmshlpos)) !XYZZ,XYYZ,XXYZ (2,1,1)
					end if
					contract=contract/valnorm31*valnormnew
				end if
			end if
			CO(:,iGTF)=orbcoeff(ibasis,:)*contract
			iGTF=iGTF+1
		end do
		ibasis=ibasis+1
	end do
end do
nbasisCar=nbasis
if (isphergau==1) nbasis=nbasis5D

if (infomode==0) then
	if (name2(itmplen-1:itmplen)=="37".or.name2(itmplen-1:itmplen)=="39") write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
	write(*,"(' Atoms:',i6,',  Basis functions:',i6,',  Orbitals:',i6,',  GTFs:',i6)") ncenter,nbasis,nmo,nprims
	if (wfntype==3) write(*,*) " This is a closed-shell system"
	if (wfntype==4) then
		write(*,*) " This is an open-shell system"
		write(*,"(' Orbitals from 1 to',i6,' are alpha, from',i6,' to',i6,' are beta')") nbasis,nbasis+1,nmo
	end if
	write(*,*)
end if
end subroutine




!!-----------------------------------------------------------------
!!---------------- Read Gaussian cube file and store in cubmat
!infomode=0 means output cube details, =1 do not, =2 also do not print loading process
!ionlygrid=1 means only read grid data, but do not perturb atom information, =0 means do all
subroutine readcube(cubname,infomode,ionlygrid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) cubname
character titleline1*79,titleline2*79,selectyn,c200tmp*200,c80tmp*80
integer infomode,ionlygrid
integer,allocatable :: mo_serial(:)
real*8,allocatable :: temp_readdata(:),tmpreadcub(:,:,:),tmpeffchg(:)
character,allocatable :: tmpelestr(:)*2

if (ionlygrid==0) ifiletype=7
open(10,file=cubname,status="old")
read(10,"(a)") titleline1
read(10,"(a)") titleline2
read(10,*) ncentertmp,orgx,orgy,orgz
read(10,*) nx,gridv1
read(10,*) ny,gridv2
read(10,*) nz,gridv3
if (ionlygrid==0) ncenter=ncentertmp
dx=gridv1(1)
dy=gridv2(2)
dz=gridv3(3)
call getgridend !Generate endx,endy,endz

mo_number=0
if (ncenter<0) then
	mo_number=1 !This cube file contains at least one MO data
	ncenter=abs(ncenter)
end if
if (ionlygrid==0) then
	if (allocated(a)) deallocate(a)
	allocate(a(ncenter))
end if
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))

if (infomode==0) then
    write(*,*)
	write(*,*) "Title line of this file:"
	write(*,"(a)") trim(titleline1)
	write(*,"(a)") trim(titleline2)
	write(*,"(/,' Total number of atoms:',i8)") ncenter
    call showgridinfo(1)
end if
if (index(titleline1,"box2cell")/=0) then
	write(*,*)
	write(*,*) "Note: Box of grid data has been set to cell!"
    call grid2cellinfo
end if

if (ionlygrid==0) then
	do i=1,ncenter
		read(10,*) a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z !%value is its charge, if ECP was used, it not equal to atomindex
	end do
	a%name=ind2name(a%index)
else if (ionlygrid==1) then
	do i=1,ncentertmp
		read(10,*) !Do not read atomic information, simply skip
	end do
end if
if (infomode<2) write(*,*)

if (mo_number==1) then
	read(10,"(i5)",advance="no") mo_number !Get actual number of MO
	if (mo_number>1) then
		allocate(mo_serial(mo_number))
		allocate(temp_readdata(nz*mo_number))
		read(10,*) mo_serial
		write(*,"(' There are ',i6,' MOs in this cube file, the indices are: ')") mo_number
		do i=1,mo_number
			write(*,"(' Index ',i6,', corresponds to MO',i6)") i,mo_serial(i)
		end do
		write(*,*) "Which MO do you want to load? Input the index, e.g. 2"
		do while(.true.)
			read(*,*) mo_select
			if (mo_select>0.and.mo_select<=mo_number) exit
			write(*,*) "Error: Invalid input, input again"
		end do
	else
		read(10,*) !Only one MO, pass the MO serial line
	end if
end if

if (infomode<2) write(*,*) "Loading grid data, please wait..."
!Note that data in cube file is recorded in reverse order as Fortran array
if (mo_number==0.or.mo_number==1) then !Commonly case, below code has the best compatibility
	allocate(tmpreadcub(nz,ny,nx))
	read(10,*) tmpreadcub(:,:,:)
	do i=1,nx
		do j=1,ny
			do k=1,nz
				cubmat(i,j,k)=tmpreadcub(k,j,i)
			end do
		end do
	end do
	deallocate(tmpreadcub)
else !Load specified of many orbitals
	do i=1,nx
		do j=1,ny
			read(10,*) temp_readdata
			cubmat(i,j,:)=temp_readdata(mo_select:size(temp_readdata):mo_number)
		end do
        if (infomode<2) call showprog(i,nx)
	end do
end if

if (all(a%charge==0)) then
	write(*,"(/,a)") " Warning: Effective nuclear charges recorded in this file are all zero. If this file was produced by CP2K, it is a bug"
	rewind(10)
    read(10,"(a)") c200tmp
    read(c200tmp,*,iostat=ierror) c80tmp,itmp !Can load a string and number, user may have specified element and effective nuclear charge
    if (ierror==0) then
		write(*,"(/,a)") " Do you want to let Multiwfn load effective nuclear charge of each element from the first line of this file? (y/n)"
		write(*,"(a)") " Note: If you choose ""y"", the effctive nuclear charges should be provided in the first line such as ""B 3 N 5 Cl 7"""
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			call numdatastr(c200tmp,ndata)
			nread=ndata/2
			allocate(tmpeffchg(nread),tmpelestr(nread))
			read(c200tmp,*) (tmpelestr(i),tmpeffchg(i),i=1,nread)
			if (nread>0) then
				write(*,*) "Loaded effctive nuclear charges"
				do i=1,nread
					write(*,"(1x,a,f12.6)") tmpelestr(i),tmpeffchg(i)
					call lc2uc(tmpelestr(i)(1:1))
					call uc2lc(tmpelestr(i)(2:2))
					do iatm=1,ncenter
						if (a(iatm)%name==tmpelestr(i)) a(iatm)%charge=tmpeffchg(i)
					end do
				end do
			end if
        end if
    end if
end if

close(10)

if (infomode==0) call showgridinfo(2) !Calculate statistical information
call guessnelec
end subroutine



!!-----------------------------------------------------------------
!!---- Read Gaussian cube file and save to cubmattmp, this is a simple version of readcube, can only be invoked after cubmat has been loaded
!don't read atomic information, don't modify loaded grid infomation such as nx/y/z,orgx/y/z...
!and don't output statistic information, don't specify coordinate for grid points...
!infomode=1: Print loading information, =2: do not print anything
!inconsis is returned value, 1 means the grid setting of this cube file is inconsistent with that of cubmat
subroutine readcubetmp(cubname,infomode,inconsis)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) cubname
integer inconsis
integer,allocatable :: mo_serial(:)
real*8,allocatable :: temp_readdata(:)
real*8 gridvectmp1(3),gridvectmp2(3),gridvectmp3(3)

open(10,file=cubname,status="old")
read(10,*)
read(10,*)
read(10,*) ncentertmp
read(10,*) nxtmp,gridvectmp1
read(10,*) nytmp,gridvectmp2
read(10,*) nztmp,gridvectmp3
inconsis=0
if (nxtmp/=nx.or.nytmp/=ny.or.nztmp/=nz.or.&
maxval(abs(gridv1-gridvectmp1))>0.005D0.or.maxval(abs(gridv2-gridvectmp2))>0.005D0.or.maxval(abs(gridv3-gridvectmp3))>0.005D0) then
	write(*,"(a)") " Error: The grid setting of this cube file is inconsistent with that of the grid data stored in memory!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
	inconsis=1
	return
end if

mo_number=0
if (ncentertmp<0) then
	mo_number=1 !This cube file contains at least one MO data
	ncentertmp=abs(ncentertmp)
end if
if (allocated(cubmattmp)) deallocate(cubmattmp)
allocate(cubmattmp(nx,ny,nz)) !nx,ny,nz is identical to cubmat that already loaded into memory

do i=1,ncentertmp
	read(10,*) !Skip a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z !%value is its charge, if ECP was used, it not equal to atomindex
end do

if (mo_number==1) then
	read(10,"(i5)",advance="no") mo_number !Get actual number of MO
	if (mo_number>1) then
		allocate(mo_serial(mo_number))
		allocate(temp_readdata(nz*mo_number))
		read(10,*) mo_serial
		write(*,"(' There are ',i6,' MOs in this grid file, the serial number are: ')") mo_number
		do i=1,mo_number
			write(*,"(' Number ',i6,' : MO= ',i5)") i,mo_serial(i)
		end do
		write(*,*) "Which MO do you want to load? Input the serial number"
		do while(.true.)
			read(*,*) mo_select
			if (mo_select>0.and.mo_select<=mo_number) exit
			write(*,*) "Invalid input, input again"
		end do
	else
		read(10,*) !Only one MO, pass the MO serial line
	end if
end if

if (infomode<2) write(*,*)
if (infomode<2) write(*,*) "Loading grid data, please wait..."
!Load data
ii=0
do i=1,nx
	do j=1,ny
		if (mo_number==0.or.mo_number==1) then
			read(10,*) cubmattmp(i,j,:)
		else !Load the specified MO from vast of MOs
			read(10,*) temp_readdata
			cubmattmp(i,j,:)=temp_readdata(mo_select:size(temp_readdata):mo_number)
		end if
	end do
	call showprog(i,nx)
end do
close(10)
end subroutine




!!-----------------------------------------------------------------
!!----------- Load .dx file (may be generated by volmap plugin in VMD)
!!infomode=0 means output dx details, =1 not
!ionlygrid=1 means only read grid data, but do not perturb any other variables, =0 means do all
subroutine readdx(dxname,infomode,ionlygrid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) dxname
character c80tmp*80
integer infomode,ionlygrid
type(content) maxv,minv

open(10,file=dxname,status="old")
call loclabel(10,"counts")
read(10,"(a)") c80tmp
itmp=index(c80tmp,"counts")
read(c80tmp(itmp+6:),*) nx,ny,nz
call loclabel(10,"origin")
read(10,*) c80tmp,orgx,orgy,orgz
orgx=orgx/b2a
orgy=orgy/b2a
orgz=orgz/b2a
call loclabel(10,"delta")
read(10,*) c80tmp,gridv1(:)
read(10,*) c80tmp,gridv2(:)
read(10,*) c80tmp,gridv3(:)
gridv1=gridv1/b2a
gridv2=gridv2/b2a
gridv3=gridv3/b2a
dx=gridv1(1)
dy=gridv2(2)
dz=gridv3(3)
call getgridend !Generate endx,endy,endz

if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
if (infomode==0) call showgridinfo(1)
write(*,*)
write(*,*) "Loading grid data, please wait..."
call loclabel(10,"object 3")
read(10,*)
read(10,*) (((cubmat(i,j,k),k=1,nz),j=1,ny),i=1,nx)
write(*,*) "Done!"
ifiletype=8
close(10)

if (infomode==0) call showgridinfo(2) !Show statistical information

end subroutine




!!-----------------------------------------------------------------
!!----------- Load Dmol3 .grd file, see format description in Material Studio help file
!infomode=0 means output grd details, =1 not
!ionlygrid=1 means only read grid data, but do not perturb any other variables, =0 means do all
subroutine readgrd(grdname,infomode,ionlygrid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) grdname
character titleline*79
integer infomode,ionlygrid
type(content) maxv,minv
if (ionlygrid==0) ifiletype=8
open(10,file=grdname,status="old")
read(10,"(a)") titleline
read(10,*)
read(10,*) flenx,fleny,flenz,anga,angb,angc !Notice that the length unit is in Angstrom in .grd file
read(10,*) nx,ny,nz !Here nx,ny,nz are total spaces (between neighbour grid point) in each direction
read(10,*) ifast,ixback,ixforw,iyback,iyforw,izback,izforw
if (ifast/=1) then !ifast=1 means x varies fastest
	write(*,*) "Error: The first integer in the fifth line must be 1!"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	stop
end if
call abc2cellv(flenx/b2a,fleny/b2a,flenz/b2a,anga,angb,angc)
gridv1=cellv1/nx
gridv2=cellv2/ny
gridv3=cellv3/nz
dx=gridv1(1)
dy=gridv2(2)
dz=gridv3(3)
nx=nx+1 !Convert the number of spacings to the number of points
ny=ny+1
nz=nz+1
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
orgx=gridv1(1)*ixback+gridv2(1)*ixback+gridv3(1)*ixback
orgy=gridv1(2)*iyback+gridv2(2)*iyback+gridv3(2)*iyback
orgz=gridv1(3)*izback+gridv2(3)*izback+gridv3(3)*izback
call getgridend !Generate endx,endy,endz

if (infomode==0) then
	write(*,"(' Title line of this file: ',a)") trim(titleline)
	write(*,*)
	call showgridinfo(1)
end if
write(*,*)
write(*,*) "Loading grid data, please wait..."
ii=0
do k=1,nz   !a(x,y,z)
	do j=1,ny
		do i=1,nx
			read(10,*) cubmat(i,j,k)
		end do
	end do
	call showprog(k,nz)
end do
write(*,*) "Done!"
close(10)

if (infomode==0) call showgridinfo(2) !Show statistical information
end subroutine


!!-----------------------------------------------------------------
!!----------- Load Dmol3 .grd file and save to cubmattmp
!inconsis is returned value, 1 means the grid setting of this cube file is inconsistent with that of cubmat
subroutine readgrdtmp(grdname,inconsis)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) grdname
open(10,file=grdname,status="old")
read(10,*)
read(10,*)
read(10,*)
read(10,*) nxtmp,nytmp,nztmp
read(10,*)
if (nxtmp+1/=nx.or.nytmp+1/=ny.or.nztmp+1/=nz) then
	write(*,"(a)") " Error: The grid setting of this grd file is inconsistent with that of the grid data stored in memory!"
	inconsis=1
	return
end if

if (allocated(cubmattmp)) deallocate(cubmattmp)
allocate(cubmattmp(nx,ny,nz))
write(*,*)
write(*,*) "Loading data, please wait..."
do k=1,nz
	do j=1,ny
		do i=1,nx
			read(10,*) cubmattmp(i,j,k)
		end do
	end do
	call showprog(k,nz)
end do
write(*,*) "Grid data loading completed!"
close(10)
end subroutine




!!---------------------------------------------------------------------
!!----------- Load vti file, which can be yielded by ParaView and GIMIC
!infomode=0 means output grd details, =1 not
!ionlygrid=1 means only read grid data, but do not perturb any other variables, =0 means do all
subroutine readvti(vtiname,infomode,ionlygrid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) vtiname
character c1000tmp*1000
integer infomode,ionlygrid
type(content) maxv,minv
if (ionlygrid==0) ifiletype=8
open(10,file=vtiname,status="old")

call loclabel(10,"<ImageData")
call readxmlline(10,c1000tmp)
i=strcharpos(c1000tmp,'"',1)
j=strcharpos(c1000tmp,'"',2)
read(c1000tmp(i+1:j-1),*) ibeg,iend,jbeg,jend,kbeg,kend
nx=iend-ibeg+1
ny=jend-jbeg+1
nz=kend-kbeg+1
i=strcharpos(c1000tmp,'"',3)
j=strcharpos(c1000tmp,'"',4)
read(c1000tmp(i+1:j-1),*) orgx,orgy,orgz
i=strcharpos(c1000tmp,'"',5)
j=strcharpos(c1000tmp,'"',6)
read(c1000tmp(i+1:j-1),*) dx,dy,dz
gridv1=0;gridv2=0;gridv3=0
gridv1(1)=dx;gridv2(2)=dy;gridv3(3)=dz
call getgridend !Generate endx,endy,endz

if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))

if (infomode==0) call showgridinfo(1)
write(*,*)
write(*,*) "Loading grid data, please wait..."
call loclabel(10,"<DataArray")
call readxmlline(10,c1000tmp)
if (index(c1000tmp,"scalars")==0) then
    write(*,*) "Error: Only vti file recording scalar data is supported!"
    write(*,*) "Press ENTER button to exit program"
    read(*,*)
    stop
end if
read(10,*) (((cubmat(i,j,k),i=1,nx),j=1,ny),k=1,nz)
write(*,*) "Done!"
close(10)

!Perform statistic
maxv%value=cubmat(1,1,1)
maxv%x=orgx
maxv%y=orgy
maxv%z=orgz

if (infomode==0) call showgridinfo(2) !Show statistical information
end subroutine


!Return content enclosed by "<" in current line and the next ">" in a file of XML format. At most read 1000 characters
!Before invoking this routine, current position must be the line containing <
subroutine readxmlline(fileid,str)
character str*1000,c1000tmp2*1000
integer fileid,itmp,jtmp
read(fileid,"(a)") str
if (index(str,'>')==0) then
    do while(.true.)
        read(10,"(a)") c1000tmp2
        str=trim(str)//trim(c1000tmp2)
        if (index(c1000tmp2,'>')/=0) exit
    end do
end if
itmp=index(str,'<')
jtmp=index(str,'>')
str=str(itmp+1:jtmp-1)
end subroutine




!!-----------------------------------------------------------------
!!------------- Read .wfn file, infomode=0 means output wavefunction property and related information, =1 not
subroutine readwfn(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
CHARACTER(LEN=*) name
character wfntitle*80,lastline*80,c80tmp*80,c80tmp2*80
real*8,allocatable :: tmpCO(:,:),tmpMOocc(:),tmpMOene(:)
integer,allocatable :: tmpMOtype(:),orbidx(:)
integer i,j,infomode
!Original .wfn format doesn't support g, however the .wfn outputted by Multiwfn, Molden2AIM and G09 since B.01 formally supports g
!Below is the g sequence used in Molden2AIM, .wfx, .molden and the .wfn outputted by Multiwfn and G09 since B.01
! 21 XXXX 22 YYYY 23 ZZZZ 24 XXXY 25 XXXZ
! 26 XYYY 27 YYYZ 28 XZZZ 29 YZZZ 30 XXYY
! 31 XXZZ 32 YYZZ 33 XXYZ 34 XYYZ 35 XYZZ
!Below is the g sequence internally used in Multiwfn, identical to .fch
! 21 ZZZZ 22 YZZZ 23 YYZZ 24 YYYZ 25 YYYY
! 26 XZZZ 27 XYZZ 28 XYYZ 29 XYYY 30 XXZZ
! 31 XXYZ 32 XXYY 33 XXXZ 34 XXXY 35 XXXX
! convGseq is used to convert g used in .wfn to the internal sequence of Multiwfn 
! I assume that h sequence is identical to .wfx and Multiwfn
integer :: convGseq(35)=(/ (0,i=1,20), 35,25,21,34,33, 29,24,26,22,32, 30,23,31,28,27 /)
ifiletype=2
imodwfn=0
open(10,file=name,status="old")
read(10,"(a)") wfntitle
!read(10,"(a8,i15,13x,i7,11x,i9)") c80tmp,nmo,nprims,ncenter
read(10,*) c80tmp,nmo,c80tmp,c80tmp,nprims,c80tmp,ncenter !Use free format for compatible with .wfn produced by Q-Chem
ibasmode=1
if (index(c80tmp,"SLATER")/=0) ibasmode=2
if (ibasmode==2) then
	write(*,"(a)") " Error: Multiwfn does not support the wfn file recording Slater type orbitals! Press ENTER button to exit"
	read(*,*)
	stop
end if
if (infomode==0.and.index(wfntitle,"Run Type")/=0) then
	write(*,"(a)") " Warning: It seems that this .wfn file was generated by ORCA. Notice that the .wfn file generated by ORCA is often non-standard, &
	&and usually makes Multiwfn crash. Using .molden file as input file instead is recommended"
end if

allocate(a(ncenter))
allocate(b(nprims),CO(nmo,nprims))
allocate(MOocc(nmo),MOene(nmo),MOtype(nmo),orbidx(nmo))

do i=1,ncenter
	read(10,"(a24,3f12.8,10x,f5.1)") c80tmp,a(i)%x,a(i)%y,a(i)%z,a(i)%charge
	read(c80tmp,*) a(i)%name
    call elename2idx(a(i)%name,a(i)%index)
	if (a(i)%index==0) then
        if (infomode==0) then
			write(*,"(/,3a,i5,'!')") " Warning: Found unknown element ",a(i)%name," with atom index of",i
			write(*,*) "This atom now is recognized as Bq (ghost atom)"
        end if
		a(i)%index=0
        a(i)%name=ind2name(0)
	end if
end do
read(10,"(20x,20i3)") (b(i)%center,i=1,nprims)
read(10,"(20x,20i3)") (b(i)%type,i=1,nprims)
do i=1,nprims	
	if (b(i)%type>=21.and.b(i)%type<=35) b(i)%type=convGseq(b(i)%type)
end do

!Read exponents
read(10,"(10x,5E14.7)") (b(i)%exp,i=1,nprims)

!From Gaussian09 B.01, if ECP is used, additional CENTER, TYPE, EXPONENTS field present to represent EDF(electron density functions) likewise .wfx file
!However such wfn file is foolish, the number of GTFs used to represent EDF is not explicitly recorded, so we must use trick to guess the number
!The coefficient of EDF is not written, how we use these EDF information? Obviously impossible! So we just skip these mad fields
read(10,"(a)") c80tmp
if (c80tmp(1:6)=="CENTRE") then
	call loclabel(10,"MO",ifound,0)
! 	ncenline=1
! 	do while(.true.) !Move to the first line of TYPE field, and count we passed how many rows (the number of lines of CENTER field)
! 		read(10,"(a)") c80tmp
! 		if (c80tmp(1:6)/="CENTRE") exit
! 		ncenline=ncenline+1
! 	end do
! 	backspace(10)
! 	backspace(10) !Return to the last line of CENTER field
! 	read(10,"(a)") c80tmp
! 	nEDFprims=(len_trim(c80tmp)-20)/3+(ncenline-1)*20
! 	do itmp=1,ncenline !Move to the first line of CENTER field, make ready to load EDF information
! 		backspace(10)
! 	end do
! 	allocate(b_EDF(nEDFprims),CO_EDF(nEDFprims))
! 	read(10,"(20x,20i3)") (b_EDF(i)%center,i=1,nEDFprims)
! 	read(10,"(20x,20i3)") (b_EDF(i)%type,i=1,nEDFprims)
! 	read(10,"(10x,5E14.7)") (b_EDF(i)%exp,i=1,nEDFprims)
else
	backspace(10)
end if

!Read orbitals
do i=1,nmo
	read(10,"(a)") c80tmp
    ipos=index(c80tmp,"MO")
    read(c80tmp(ipos+2:),*) orbidx(i)
    ipos=index(c80tmp,'=')
    read(c80tmp(ipos+1:),*) MOocc(i)
    ipos=index(c80tmp,'=',back=.true.)
    read(c80tmp(ipos+1:),*) MOene(i)
	read(10,"(5E16.8)") CO(i,:)
end do
read(10,*)
!Use free format to read in energy and virial to ensure compatibility
read(10,"(a)") lastline
iequalsign1=0
iequalsign2=0
do i=1,80
	if (lastline(i:i)=='=') then
		iequalsign1=i
		exit
	end if
end do
do i=80,1,-1
	if (lastline(i:i)=='=') then
		iequalsign2=i
		exit
	end if
end do
totenergy=0
virialratio=2
if (iequalsign1/=0) read(lastline(iequalsign1+1:),*) totenergy
if (iequalsign1==0) write(*,*) "Warning: Unable to find system energy in this file!"
if (iequalsign2/=0) read(lastline(iequalsign2+1:),*) virialratio
if (iequalsign2==0) write(*,*) "Warning: Unable to find Virial in this file!"
call loclabel(10,"$MOSPIN",ifoundmospin,0) !Also read spin-type from $MOSPIN $END field outputted by Molden2AIM since v2.0.5, if detected
if (ifoundmospin==1) then !Have defined spin-type explicitly, don't reset spin-type by guessing later
	if (infomode==0) write(*,"(a)") " Note: Found $MOSPIN field, orbital spin types are directly loaded rather than automatically determined"
	read(10,*)
	read(10,*) MOtype
	where (MOtype==3) MOtype=0
end if

!Load cell information if any. May be 1/2/3-dimension
call loclabel(10,"[Cell]",ifound)
if (ifound==0) call loclabel(10,"[cell]",ifound)
if (ifound==0) call loclabel(10,"[CELL]",ifound)
if (ifound==1) then
    read(10,*)
    read(10,*) cellv1
    ifPBC=1
    read(10,*,iostat=ierror) cellv2
    if (ierror/=0) then
        cellv2=0
    else
        ifPBC=ifPBC+1
        read(10,*,iostat=ierror) cellv3
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

close(10)

!Determine MOtype of all orbitals and wfntype from occupancy
!write(*,*) MOocc
!pause
if (sum(MOocc)==2*nmo.and.all(int(MOocc)==MOocc)) then
	wfntype=0 !This is restricted wavefunction
	MOtype=0
else if (sum(MOocc)==nmo.and.all(int(MOocc)==MOocc)) then
	wfntype=1 !This is unrestricted wavefunction
	if (ifoundmospin==0) then
		MOtype=1 !Set all MO is alpha
        if (nmo>1) then
		    do i=2,nmo !If nmo=1, i will be set to 2, and no errors will appear
			    if ( MOene(i)<MOene(i-1) .or. orbidx(i)>orbidx(i-1)+1 ) exit !If energy or index is not contiguous w.r.t. last orbital, beta orbital must be encountered
		    end do
		    MOtype(i:nmo)=2 !Beta
        end if
	end if
else if (any(MOocc/=int(MOocc))) then
	if (nint(maxval(MOocc))==2) then !Maximum occupation close to 2, so considered as restricted multiconfiguration wavefunction
		wfntype=3
		MOtype=0
	else
		wfntype=4 !This is unrestricted multiconfiguration wavefunction
		if (ifoundmospin==0) then
			MOtype=1
            if (nmo>1) then
			    do i=2,nmo
				    if (MOocc(i)>MOocc(i-1) .or. orbidx(i)>orbidx(i-1)+1 ) exit !If occupation or index is not contiguous w.r.t. last orbital, beta orbital must be encountered
			    end do
		        MOtype(i:nmo)=2 !Beta
            end if
		end if
	end if
else
	wfntype=2 !This is RO wavefunction
	MOtype=0
	do i=1,nmo
		if (MOocc(i)==1) MOtype(i)=1 !alpha
	end do
end if

!Count electrons
call updatenelec

!Sort orbitals so that the orbitals with same spin-type are contiguous, because the wfn file outputted by Molden2AIM is not always in line with this convention
if (ifoundmospin==1.and.(wfntype==1.or.wfntype==4)) then
	allocate(tmpCO(nmo,nprims),tmpMOtype(nmo),tmpMOocc(nmo),tmpMOene(nmo))
	ipos=0
	do itime=1,2
		do imo=1,nmo
			if ((itime==1.and.MOtype(imo)==1).or.(itime==2.and.MOtype(imo)==2)) then !Move alpha orbitals to tmp arrays
				ipos=ipos+1
				tmpCO(ipos,:)=CO(imo,:)
				tmpMOocc(ipos)=MOocc(imo)
				tmpMOene(ipos)=MOene(imo)
				tmpMOtype(ipos)=MOtype(imo)
			end if
		end do
	end do
	CO=tmpCO
	MOocc=tmpMOocc
	MOene=tmpMOene
	MOtype=tmpMOtype
	if (infomode==0) write(*,*) "Note: Sequence of orbitals has been sorted according to spin type"
	deallocate(tmpCO,tmpMOocc,tmpMOene,tmpMOtype)
end if

!Introduce EDF information
if (any(a%index/=nint(a%charge))) then
	if (isupplyEDF==0) then !Do nothing
		continue
	else if (isupplyEDF==2) then !Supply EDF from bulit-in library
		call readEDFlib(infomode)
	end if
end if

!Summary
if (infomode==0) then
	write(*,*)
	write(*,"(' Total energy:',f21.12,' Hartree,   Virial ratio:',f12.8)") totenergy,virialratio
	write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
	write(*,"(' Net charge:',f12.5,'    Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
	write(*,"(' The number of orbitals:',i6,',  Atoms:',i7,',  GTFs:',i7)") nmo,ncenter,nprims
	if (wfntype==0) write(*,"(' This is a restricted closed-shell single-determinant wavefunction')")
	if (wfntype==1) write(*,"(' This is an unrestricted single-determinant wavefunction')")
	if (wfntype==2) write(*,"(' This is a restricted open-shell wavefunction')")
	if (wfntype==3) write(*,"(' This is a restricted multiconfiguration wavefunction')")
	if (wfntype==4) write(*,"(' This is an unrestricted multiconfiguration wavefunction')")
	if (wfntype==1.or.wfntype==4) then
		do i=1,nmo
			if (MOtype(i)==2) exit
		end do
		if (any(MOtype==2)) then
            write(*,"(' Orbitals from 1 to',i6,' are alpha type, from',i6,' to',i6,' are beta type')") i-1,i,nmo
        else
            write(*,"(' All orbitals are alpha type')")
        end if
	end if
	write(*,"(' Title line of this file: ',a)") trim(wfntitle)
end if
end subroutine



!!-----------------------------------------------------------------
!!-------- Read .wfx. mode=0 means output related information, =1 Do not output
subroutine readwfx(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
CHARACTER(LEN=*) name
character spintype*20
integer infomode
!Below is the g sequence used in Molden2AIM, .wfx, .molden and the .wfn outputted by Multiwfn and G09 since B.01
! 21 XXXX 22 YYYY 23 ZZZZ 24 XXXY 25 XXXZ
! 26 XYYY 27 YYYZ 28 XZZZ 29 YZZZ 30 XXYY
! 31 XXZZ 32 YYZZ 33 XXYZ 34 XYYZ 35 XYZZ
!Below is the g sequence internally used in Multiwfn, identical to .fch
! 21 ZZZZ 22 YZZZ 23 YYZZ 24 YYYZ 25 YYYY
! 26 XZZZ 27 XYZZ 28 XYYZ 29 XYYY 30 XXZZ
! 31 XXYZ 32 XXYY 33 XXXZ 34 XXXY 35 XXXX
! convGseq is used to convert g used in .wfx to the internal sequence of Multiwfn
! PS: spdfh sequence in .wfx is identical to Multiwfn
integer :: convGseq(35)=(/ (0,i=1,20), 35,25,21,34,33, 29,24,26,22,32, 30,23,31,28,27 /)
ifiletype=3
imodwfn=0
open(10,file=name,status="old")
call loclabel(10,"<Number of Nuclei>")
read(10,*)
read(10,*) ncenter
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
call loclabel(10,"<Number of Primitives>")
read(10,*)
read(10,*) nprims
allocate(b(nprims))
call loclabel(10,"<Number of Occupied Molecular Orbitals>")
read(10,*)
read(10,*) nmo
allocate(MOocc(nmo),MOene(nmo),MOtype(nmo),CO(nmo,nprims))
call loclabel(10,"<Nuclear Names>")
read(10,*)
read(10,*) a%name
do i=1,ncenter !Multiwfn does not allow number included in atom name
	if (iachar(a(i)%name(2:2))<=57.and.iachar(a(i)%name(2:2))>=48) a(i)%name(2:2)=' '
end do
call loclabel(10,"<Atomic Numbers>")
read(10,*)
read(10,*) a%index
call loclabel(10,"<Nuclear Charges>")
read(10,*)
read(10,*) a%charge
call loclabel(10,"<Nuclear Cartesian Coordinates>")
read(10,*)
do i=1,ncenter
    read(10,*) a(i)%x,a(i)%y,a(i)%z
end do
call loclabel(10,"<Number of Electrons>")
read(10,*)
read(10,*) nelec
call loclabel(10,"<Number of Alpha Electrons>")
read(10,*)
read(10,*) naelec
call loclabel(10,"<Number of Beta Electrons>")
read(10,*)
read(10,*) nbelec
call loclabel(10,"<Primitive Centers>")
read(10,*)
read(10,*) b%center
call loclabel(10,"<Primitive Types>")
read(10,*)
read(10,*) b%type
!The g sequence in .wfx is not identical to Multiwfn, convert them here
do i=1,nprims	
	if (b(i)%type>=21.and.b(i)%type<=35) b(i)%type=convGseq(b(i)%type)
end do
call loclabel(10,"<Primitive Exponents>")
read(10,*)
read(10,*) b%exp
!------ Process EDF information
call loclabel(10,"<Number of EDF Primitives>",ifound)
if (ifound==1.and.readEDF==1) then
	write(*,*) "Loading electron density functions (EDF) field in this file..."
	read(10,*)
	read(10,*) nEDFprims
	allocate(b_EDF(nEDFprims),CO_EDF(nEDFprims))
	call loclabel(10,"<EDF Primitive Centers>")
	read(10,*)
	read(10,*) b_EDF%center
	call loclabel(10,"<EDF Primitive Types>")
	read(10,*)
	read(10,*) b_EDF%type !We assume all the type index is 1 (S type)
	if (maxval(b_EDF%type)>1) then
		write(*,*) "Error: All GTFs of electron density function must be S type! Press ENTER button to exit"
		read(*,*)
		stop
	end if
	call loclabel(10,"<EDF Primitive Exponents>")
	read(10,*)
	read(10,*) b_EDF%exp
	call loclabel(10,"<EDF Primitive Coefficients>")
	read(10,*)
	read(10,*) CO_EDF
	call loclabel(10,"<Number of Core Electrons>")
	read(10,*)
	read(10,*) nEDFelec
    allocate(nEDFelecatm(ncenter))
    nEDFelecatm(:)=a(:)%index-a(:)%charge
	if (infomode==0) write(*,"(a,i6,a)") " Note: EDF information represents",nEDFelec," inner-core electrons"
end if
call loclabel(10,"<Molecular Orbital Occupation Numbers>")
read(10,*)
read(10,*) MOocc
call loclabel(10,"<Molecular Orbital Energies>")
read(10,*)
read(10,*) MOene
call loclabel(10,"<Molecular Orbital Spin Types>")
read(10,*)
do i=1,nmo
	read(10,"(a20)") spintype
	if (adjustl(spintype)=="Alpha and Beta") MOtype(i)=0 !adjustl is needed, because the wfx outputted by ORCA is non-standard
	if (adjustl(spintype)=="Alpha") MOtype(i)=1
	if (adjustl(spintype)=="Beta") MOtype(i)=2
end do
call loclabel(10,"<Molecular Orbital Primitive Coefficients>")
read(10,*)
do i=1,nmo
	read(10,*)
	read(10,*)
	read(10,*)
	read(10,*) CO(i,:)
end do
call loclabel(10,"<Energy = T + Vne + Vee + Vnn>",ifound,0) !Don't rewind, otherwise Multiwfn will scan lots of information of MO field
read(10,*)
read(10,*) totenergy
call loclabel(10,"<Virial Ratio (-V/T)>",ifound,0)
read(10,*)
read(10,*) virialratio
close(10)

!Angular moment of GTF should be no higher than h
if (any(b%type>56)) then
	write(*,"(' Warning: Angular moment of one or more GTFs exceeds h, Multiwfn is unable to deal with this case! Its/their contributions will be discarded')")
	write(*,*) "Press ENTER button to continue"
	read(*,*)
	do iGTF=1,nprims
		if (b(iGTF)%type>56) then
			b(iGTF)%type=1 !Assume it is S type
			CO(:,iGTF)=0D0
		end if
	end do
end if

if ( all(MOocc==nint(MOocc)) ) then
	wfntype=2
	if (nmo==nint(nelec)) wfntype=1
	if (nmo==nint(nelec)/2) wfntype=0
else !post-HF
	if (naelec==nbelec) wfntype=3
	if (naelec/=nbelec) wfntype=4
end if

if (infomode==0) then
	write(*,*)
	write(*,"(' Total energy:',f21.12,' Hartree,   Virial ratio:',f12.8)") totenergy,virialratio
	write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
	write(*,"(' Number of orbital:',i6,',  Atoms:',i7,',  GTFs:',i7)") nmo,ncenter,nprims
	if (wfntype==0) write(*,"(' This is restricted closed-shell single-determinant wavefunction')")
	if (wfntype==1) write(*,"(' This is unrestricted single-determinant wavefunction')")
	if (wfntype==2) write(*,"(' This is restricted open-shell wavefunction')")
	if (wfntype==3) write(*,"(' This is restricted multiconfiguration wavefunction')")
	if (wfntype==4) write(*,"(' This is unrestricted multiconfiguration wavefunction')")
	if (wfntype==1.or.wfntype==4) then
		do i=1,nmo
			if (MOtype(i)==2) exit
		end do
		write(*,"(' Orbitals from 1 to',i6,' are alpha type, from',i6,' to',i6,' are beta type')") i-1,i,nmo
	end if
	write(*,*)
end if
end subroutine



!!------- Load EDF information from external atomic .wfx files
subroutine readEDFatmwfx
use defvar
use util
implicit real*8 (a-h,o-z)
character elewfxfilename(110)*200,c200tmp*200
real*8,allocatable :: EDFCOtmp(:),EDFexptmp(:)
integer,allocatable :: EDFtypetmp(:)
integer atmsel(nelesupp,ncenter),natmsel(nelesupp)
iwfxtime=1
nEDFprims=0
nEDFelec=0
allocate(nEDFelecatm(ncenter))
nEDFelecatm=0
do while(.true.)
	write(*,*) "Load the inner-core density (EDF information) for which element? e.g. Fe"
	write(*,*) "You can also input atomic indices, e.g. 5,8-10,31 means selecting 5,8,9,10,31"
	write(*,*) "Note: If finished, input ""q"""
	read(*,"(a)") c200tmp
	itmp=ichar(c200tmp(1:1))
	if (c200tmp=='q') then
		exit
	else if (itmp>=48.and.itmp<=57) then !Inputted is atomic indices
		call str2arr(c200tmp,natmsel(iwfxtime),atmsel(iwfxtime,:))
	else !Inputted is element name
		call lc2uc(c200tmp(1:1)) !Make the first/second character in upper/lower case
		call uc2lc(c200tmp(2:2))
		natmsel(iwfxtime)=0
		do iatm=1,ncenter
			if (a(iatm)%name==c200tmp(1:2)) then
				natmsel(iwfxtime)=natmsel(iwfxtime)+1
				atmsel(iwfxtime,natmsel(iwfxtime))=iatm
			end if
		end do
	end if
	if (natmsel(iwfxtime)==0) then
		write(*,*) "No atoms are selected, input again"
		write(*,*)
	else if (natmsel(iwfxtime)>0) then
		write(*,"(' The number of atoms selected is',i7,',  including:')") natmsel(iwfxtime)
		write(*,"(12i6)") atmsel(iwfxtime,1:natmsel(iwfxtime))
		write(*,*)
		write(*,*) "Load EDF information from which file? e.g. C:\ltwd\Fe_lanl2.wfx"
		do while(.true.)
			read(*,*) elewfxfilename(iwfxtime)
			inquire(file=elewfxfilename(iwfxtime),exist=alive)
			if (alive) exit
			write(*,*) "Cannot find this file, input again"
		end do
		open(10,file=elewfxfilename(iwfxtime),status="old") !Count how many EDF GTFs in this file
		call loclabel(10,"<Number of EDF Primitives>",ifound)
		if (ifound==0) then
			write(*,*) "Error: Unable to find EDF information from this file!"
			cycle
		end if
		read(10,*)
		read(10,*) nEDFtmp
		call loclabel(10,"<Number of Core Electrons>")
		read(10,*)
		read(10,*) nEDFelectmp
		close(10)
		write(*,"(' The number of EDF primitives in this file is',i5,/)") nEDFtmp
		nEDFprims=nEDFprims+nEDFtmp*natmsel(iwfxtime)
		nEDFelec=nEDFelec+nEDFelectmp*natmsel(iwfxtime)
        nEDFelecatm(atmsel(iwfxtime,1:natmsel(iwfxtime)))=nEDFelectmp
		iwfxtime=iwfxtime+1
	end if
end do
nwfxtime=iwfxtime-1
write(*,"(' The total number of EDF primitives is',i7)") nEDFprims
write(*,"(' The total number of inner-core electrons represented by EDF is',i8)") nEDFelec
allocate(b_EDF(nEDFprims),CO_EDF(nEDFprims))
ipos=1
do iwfxtime=1,nwfxtime
	open(10,file=elewfxfilename(iwfxtime),status="old")
	call loclabel(10,"<Number of EDF Primitives>",ifound)
	read(10,*)
	read(10,*) nEDFtmp
	allocate(EDFCOtmp(nEDFtmp),EDFexptmp(nEDFtmp),EDFtypetmp(nEDFtmp))
	call loclabel(10,"<EDF Primitive Types>")
	read(10,*)
	read(10,*) EDFtypetmp
	call loclabel(10,"<EDF Primitive Exponents>")
	read(10,*)
	read(10,*) EDFexptmp
	call loclabel(10,"<EDF Primitive Coefficients>")
	read(10,*)
	read(10,*) EDFCOtmp
	do iatm=1,natmsel(iwfxtime)
		b_EDF(ipos:ipos+nEDFtmp-1)%type=EDFtypetmp
		b_EDF(ipos:ipos+nEDFtmp-1)%exp=EDFexptmp
		b_EDF(ipos:ipos+nEDFtmp-1)%center=atmsel(iwfxtime,iatm)
		CO_EDF(ipos:ipos+nEDFtmp-1)=EDFCOtmp
		ipos=ipos+nEDFtmp
	end do
	deallocate(EDFCOtmp,EDFexptmp,EDFtypetmp)
	close(10)
end do
write(*,*) "The EDF information has been loaded"
end subroutine


!!------ Load EDF information from EDFlib provided by Zork
!See http://bbs.keinsci.com/forum.php?mod=viewthread&tid=5354 for description
!infomode=0/1 show/don't show info
subroutine readEDFlib(infomode)
use defvar
implicit real*8 (a-h,o-z)
real*8 EDFcoeff(100),EDFexp(100)
if (infomode==0) then
    write(*,*) "Loading electron density functions (EDF) information from built-in EDF library"
    write(*,*) "(The library is available at https://github.com/zorkzou/Molden2AIM)"
end if
nEDFprims=0
nEDFelec=0
allocate(nEDFelecatm(ncenter))
nEDFelecatm=0
!First time, find total number of EDF GTFs so that b_EDF and CO_EDF can be allocated
do iatm=1,ncenter
	natmcore=a(iatm)%index-nint(a(iatm)%charge)
	if (natmcore==0) cycle !Full electron
	if (a(iatm)%index==0.or.a(iatm)%charge==0) cycle !Skip Bq
	call EDFLIB(a(iatm)%index,natmcore,nfun,EDFexp,EDFcoeff)
	if (infomode==0) write(*,"(1x,a,'(',i5,')      Core electrons:',i3,'     EDF primitive GTFs:',i3)") a(iatm)%name,iatm,natmcore,nfun
	if (nfun==0) then
		if (infomode==0) write(*,*) "Warning: Unable to find proper EDF information for this atom!"
    else
	    nEDFelecatm(iatm)=natmcore
    end if
	nEDFprims=nEDFprims+nfun
end do
nEDFelec=sum(nEDFelecatm(:))
if (infomode==0) then
    write(*,"(' The number of total inner-core electrons:',i6)") nEDFelec
    write(*,"(' The number of total EDF primitive GTFs:',i6)") nEDFprims
end if

allocate(b_EDF(nEDFprims),CO_EDF(nEDFprims))
ifun=0
b_EDF%type=0
do iatm=1,ncenter
	natmcore=a(iatm)%index-nint(a(iatm)%charge)
	if (natmcore==0) cycle
	call EDFLIB(a(iatm)%index,natmcore,nfun,EDFexp,EDFcoeff)
	if (nfun==0) cycle !Didn't find corresponding EDF information
	b_EDF(ifun+1:ifun+nfun)%exp=EDFexp(1:nfun)
	b_EDF(ifun+1:ifun+nfun)%center=iatm
	CO_EDF(ifun+1:ifun+nfun)=EDFcoeff(1:nfun)
	ifun=ifun+nfun
end do
if (infomode==0) write(*,*) "Loading EDF library finished!"
end subroutine



!!-----------------------------------------------------------------
!!--------- Read Molden input file, get coordinate, basis function and GTF information
!Known issue:
!CFour sometimes fail (e.g. benzene)
subroutine readmolden(name,infomode) !infomode=0 means output info, =1 silent
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character c80*80,c80tmp*80,c200tmp*200,symtmp*10
integer,allocatable :: shelltype(:),shellcon(:),shell2atom(:) !The definition of shelltype is identical to .fch
integer :: s2f(-5:5,21)=0 !Give shell type & orbital index to get functype
real*8,allocatable :: primexp(:),concoeff(:)
real*8,allocatable :: amocoeff(:,:),bmocoeff(:,:)
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
!For backing up spherical basis functions
integer,allocatable :: shelltype5D(:),MOtype5D(:)
character(len=4),allocatable :: MOsym5D(:)
real*8,allocatable :: CObasa5D(:,:),CObasb5D(:,:),MOene5D(:),MOocc5D(:),CO5D(:,:)
real*8,external :: normgau
ifiletype=9
imodwfn=0
s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
s2f(-1,1:4)=(/ 1,2,3,4 /)
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
!---------- The sequence of f functions in Multiwfn (=wfn=wfx) is not identical to Molden, so convert here
!11  12  13  14  15  16  17  18  19  20  !Multiwfn sequence
!XXX YYY ZZZ XXY XXZ YYZ XYY XZZ YZZ XYZ
!xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz !Molden sequence
s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /)
!---------- The sequence of g functions in Multiwfn (=fch) is not identical to Molden, so convert here
! 21   22   23   24   25   26   27   28   29   30   31   32   33   34   35  !Multiwfn sequence
!ZZZZ YZZZ YYZZ YYYZ YYYY XZZZ XYZZ XYYZ XYYY XXZZ XXYZ XXYY XXXZ XXXY XXXX
!xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy !Molden sequence
s2f(4,1:15)=(/ 35,25,21,34,33,29,24,26,22,32,30,23,31,28,27 /)
!---------- The sequence of h functions in Multiwfn (=fch=wfx) is not identical to Molden, so convert here
!Note that h angular moment is not formally supported by .molden format!!!
! 36    37    38    39    40    41    42    43    44    45    46  !Multiwfn sequence
!ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ XYYYY
!xxxxx yyyyy zzzzz xxxxy xxxxz xyyyy xzzzz yyyyz yzzzz xxxyy xxxzz !Molden sequence
! 47    48    49    50    51    52    53    54    55    56  !Multiwfn sequence
!XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
!xxyyy xxzzz yyyzz yyzzz xxxyz xyyyz xyzzz xxyyz xxyzz xyyzz !Molden sequence
!s2f(5,1:21)=(/ 56,41,36,55,54,46,42,40,37,53,51,50,47,39,38,52,45,43,49,48,44 /)
!I forgot how the "Molden sequence" was previously determined by me... So, do not alter the sequence, in this case the molden generated by Multiwfn can be normally loaded
!Respective of if reordering the sequence, Multiwfn is unable to correctly load wavefunction from .molden containing h generated by ORCA (at least for 4.2.1)
forall(i=1:21) s2f(5,i)=i+35

call gensphcartab(2,conv5d6d,conv7f10f,conv9g15g,conv11h21h)

open(10,file=name,status="old")

if (infomode==0) write(*,*) "Loading various information of the wavefunction"

!!!!! Load cell information if any. May be 1/2/3-dimension
call loclabel(10,"[Cell]",ifound,maxline=50)
if (ifound==0) call loclabel(10,"[cell]",ifound,maxline=50)
if (ifound==0) call loclabel(10,"[CELL]",ifound,maxline=50)
if (ifound==1) then
    read(10,*)
    read(10,"(a)") c80tmp
    read(c80tmp,*,iostat=ierror) alen,blen,clen,anga,angb,angc
    if (ierror==0) then !Load as cell parameters
        call abc2cellv(alen/b2a,blen/b2a,clen/b2a,anga,angb,angc)
        ifPBC=3
    else !Load as cell vectors
		ipos=index(c80tmp,'A')
        if (ipos/=0) c80tmp(ipos:ipos)=' ' !The line may be "A    20.73430000     0.00000000     0.00000000", remove 'A'
        read(c80tmp,*) cellv1
        ifPBC=1
		read(10,"(a)") c80tmp
		ipos=index(c80tmp,'B')
        if (ipos/=0) c80tmp(ipos:ipos)=' ' !The line may be "B     0.00000000    20.73430000     0.00000000", remove 'B'
        read(c80tmp,*,iostat=ierror) cellv2
        if (ierror/=0) then
            cellv2=0
        else
            ifPBC=ifPBC+1
			read(10,"(a)") c80tmp
			ipos=index(c80tmp,'C')
			if (ipos/=0) c80tmp(ipos:ipos)=' ' !The line may be "C     0.00000000     0.00000000    20.73430000", remove 'C'
            read(c80tmp,*,iostat=ierror) cellv3
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
else
	call loadcellinfo_txt !Load cell information from [Cell].txt in current folder if available
end if

!!!!! Load atom information
call loclabel(10,"[Atoms]",ifound)
if (ifound==0) call loclabel(10,"[ATOMS]",ifound)
ilenunit=1 !Default length unit is a.u. =2 means Angstrom
read(10,"(a)") c80
if (index(c80,"Angs")/=0) ilenunit=2
ncenter=0
do while(.true.) !Count the number of atoms
	read(10,"(a)") c80
	if (c80==" ".or.index(c80,"[")/=0) then !Have passed atom field
		exit
	end if
	ncenter=ncenter+1
end do
if (allocated(a)) deallocate(a)
allocate(a(ncenter))
call loclabel(10,"[Atoms]",ifound) !Return to [Atoms]
if (ifound==0) call loclabel(10,"[ATOMS]",ifound)
read(10,"(a)") c80
!NOTICE: molden input file has a severe drawback, namely atomic charge is not explicitly recorded, this will be problematic when ECP is used
!In Multiwfn, atomic index is determined according to atomic name, while "atom number" column is read as atomic charge. Therefore, if you already know
!current file have used ECP, then you can open the file and manually change the atomic number to atomic charge.
do iatm=1,ncenter
	read(10,"(a)") c200tmp
	read(c200tmp,*,iostat=ierror) c80,nouse,a(iatm)%charge,a(iatm)%x,a(iatm)%y,a(iatm)%z
    if (ierror/=0) then
        read(c200tmp,*) c80,nouse,c80tmp,a(iatm)%x,a(iatm)%y,a(iatm)%z
    end if
	call lc2uc(c80(1:1)) !Convert to upper case
	call uc2lc(c80(2:2)) !Convert to lower case
	do i=0,nelesupp
		if (iachar(c80(2:2))>=48.and.iachar(c80(2:2))<=57) c80(2:2)=' ' !Dalton adds number after element name
		if (c80(1:2)==ind2name(i)) then
			a(iatm)%index=i
			a(iatm)%name=ind2name(i)
			exit
		end if
	    if (i==nelesupp) then
            if (infomode==0) then
		        write(*,"(/,3a,i5,'!')") " Warning: Found unknown element ",trim(c80(1:2))," with atom index of",iatm
		        write(*,*) "This atom now is recognized as Bq (ghost atom)"
            end if
		    a(iatm)%index=0
            a(iatm)%name=ind2name(0)
	    end if
	end do
end do
if (ilenunit==2) then !Angstrom->a.u.
	a%x=a%x/b2a
	a%y=a%y/b2a
	a%z=a%z/b2a
end if

!If user manually specified number of valence electrons of specific element, use it
call loclabel(10,"[Nval]",ifound,maxline=500)
if (ifound==0) call loclabel(10,"[nval]",ifound)
if (ifound==1) then
    read(10,*)
    do while(.true.)
        read(10,"(a)",iostat=ierror) c80
        if (ierror/=0.or.c80==" ".or.index(c80,'[')/=0) exit
        c80=adjustl(c80)
	    call lc2uc(c80(1:1)) !Convert to upper case
	    call uc2lc(c80(2:2)) !Convert to lower case
	    do iele=0,nelesupp
            if (c80(1:2)==ind2name(iele)) exit
        end do
        if (iele<=nelesupp) then
            read(c80,*) c80tmp,ichg
            if (infomode==0) write(*,"(' Note: Nuclear charge of ',a,' has been set to',i4)") ind2name(iele),ichg
            do iatm=1,ncenter
                if (a(iatm)%index==iele) a(iatm)%charge=ichg
            end do
        end if
    end do
end if

!Detect the Molden input file is produced by which program. Special treatment is needed for ORCA, xtb and CFOUR
!If iorca/ixtb/icfour=1, that means this file was generated by orca/xtb/CFOUR
rewind(10)
iorca=0
ixtb=0
icfour=0
do while(.true.)
	read(10,"(a)") c80
	if (index(c80,"[GTO]")/=0) exit
	call struc2lc(c80)
	if (index(c80,"orca")/=0) then
		iorca=1
		exit
	else if (index(c80,"cfour")/=0) then
		icfour=1
		exit
	end if
end do
if (iorca==1.and.infomode==0) write(*,"(a)") " This file is recognized to be generated by ORCA because there is ""orca"" word in title line. &
&Special treatments are applied..."
if (icfour==1.and.infomode==0) write(*,"(a)") " This file is recognized to be generated by CFour because there is ""cfour"" word in title line. &
&Special treatments are applied..."
if (iorca==0.and.icfour==0) then !Test if it may be generated by xtb using the first three lines
	rewind(10)
	read(10,"(a)") c80
	if (c80=="[Molden Format]") then
		read(10,"(a)") c80
		if (c80=="[Title]") then
			read(10,"(a)") c80
			if (c80=="[Atoms] AU") then
				ixtb=1
				if (ixtb==1.and.infomode==0) write(*,*) "This file is found to be generated by xtb! Special treatment is applied..."
			end if
		end if
	end if
end if

!For ORCA >= 6.0, if there is [Pseudo] after [Atoms] field, load it
if (iorca==1) then
	call loclabel(10,"[Atoms]")
    call skiplines(10,ncenter+1)
    read(10,*) c80tmp
    if (index(c80tmp,"[Pseudo]")/=0) then
		do while(.true.)
			read(10,"(a)",iostat=ierror) c80
			if (ierror/=0.or.c80==" ".or.index(c80,'[')/=0) exit
            read(c80,*) c80tmp,idx,ichg
            a(idx)%charge=ichg
            if (infomode==0) write(*,"(' Note: Nuclear charge of atom',i6,' has been set to',i4)") idx,ichg
		end do
    end if
end if

!!!!! Load basis set, and build up GTF information
if (infomode==0) write(*,*) "Loading basis set definition..."
call loclabel(10,"[GTO]",ifound)
if (ifound==0) then
	write(*,*) "Error: Unable to locate [GTO] field! Press ENTER button to exit"
	write(*,*) "Note: [STO] is currently not supported"
	read(*,*)
	stop
end if

!First time, we count how many shells are there to allocate proper size of allocatable arrays
nshell=0
nprimshell=0
read(10,*)
do while(.true.)
	read(10,*) iatm
	do while(.true.)
		read(10,"(a)") c200tmp
		if (c200tmp==" ") exit !This atom does not carry any basis function
		read(c200tmp,*) c80,ncon
		nshell=nshell+1
		iaddnprmsh=0
		do ish=1,ncon
			read(10,*) tmpv1,tmpv2
			if (tmpv2/=0D0) iaddnprmsh=iaddnprmsh+1 !Many GTF shells outputted by Molpro have zero concontraction coefficient, they shouldn't be read in
		end do
		nprimshell=nprimshell+iaddnprmsh
		if (index(c80,"sp")/=0.or.index(c80,"SP")/=0) then !sp shell will be separated as s and p
			nshell=nshell+1
			nprimshell=nprimshell+iaddnprmsh
		end if
		read(10,"(a)") c80
		if (c80==" ") exit
		backspace(10)
	end do
	read(10,"(a)") c80
	if (c80==" ".or.index(c80,"[")/=0) exit !Finished reading [GTO] field
	backspace(10)
end do

!Second time, read basis set information actually
allocate(shelltype(nshell),shellcon(nshell),shell2atom(nshell))
allocate(primexp(nprimshell),concoeff(nprimshell))
call loclabel(10,"[GTO]",ifound)
shellcon=0
ishell=0
iprimshell=0
read(10,*) 
do while(.true.)
	read(10,*) iatm
	do while(.true.)
		read(10,"(a)") c200tmp
		if (c200tmp==" ") exit !This atom does not carry any basis function
		read(c200tmp,*) c80,ncon
		ishell=ishell+1
		shell2atom(ishell)=iatm
		!Determine shell type of basis function, here we first assume they are all Cartesian type
		if (index(c80,"sp")/=0.or.index(c80,"SP")/=0) then
			shelltype(ishell)=-1
		else if (index(c80,"s")/=0.or.index(c80,"S")/=0) then
			shelltype(ishell)=0
		else if (index(c80,"p")/=0.or.index(c80,"P")/=0) then
			shelltype(ishell)=1
		else if (index(c80,"d")/=0.or.index(c80,"D")/=0) then
			shelltype(ishell)=2
		else if (index(c80,"f")/=0.or.index(c80,"F")/=0) then
			shelltype(ishell)=3
		else if (index(c80,"g")/=0.or.index(c80,"G")/=0) then
			shelltype(ishell)=4
		else if (index(c80,"h")/=0.or.index(c80,"H")/=0) then
			shelltype(ishell)=5
		else if (index(c80,"i")/=0.or.index(c80,"I")/=0) then
			write(*,"(a)") " Error: Multiwfn supports angular moment of basis function up to h, but this file contains basis functions with angular moment of i, &
            Multiwfn cannot deal with this case. Please use smaller basis set."
            write(*,*) "Press ENTER button to exit"
            read(*,*)
            stop
		end if
		iprimshellold=iprimshell
		do ish=1,ncon !Read exponents and contraction coefficients. For SP, here load the S one
			read(10,*) exptmp,concoefftmp
			if (concoefftmp==0D0) cycle !The shell with zero contraction coefficients will be ripped out
			iprimshell=iprimshell+1
			shellcon(ishell)=shellcon(ishell)+1
			primexp(iprimshell)=exptmp
			if (iorca==1.or.ixtb==1) then !ORCA/xtb doesn't present SP shell in Molden input file, so don't worry about -1
				!The normalization factor of spherical harmonic GTFs are weirdly multiplied into contraction coefficients,&
				!so here normalization factor is eliminated out from contraction coefficients to meet common convention
				rnorm=renormgau_ORCA(primexp(iprimshell),shelltype(ishell))
				concoefftmp=concoefftmp/rnorm
			end if
			concoeff(iprimshell)=concoefftmp
! 			write(*,"(2i4,2f18.10)") iatm,shelltype(ishell),primexp(iprimshell),concoefftmp
		end do
		nprmshadd=iprimshell-iprimshellold
		!For ORCA (or may be also other programs), d,f,g basis functions are not properly normalized, &
		!e.g. for ORCA at current stage d are normalized to 3, so renormalization (i.e. modifying contraction coefficients) is required and thus performed here
		!The LCAO coefficients always correspond to normalized basis functions (but someones normalized to -1, will be fixed at later stage)
		if (shelltype(ishell)/=-1) call renormmoldengau(nprmshadd,shelltype(ishell),&
		primexp(iprimshell-nprmshadd+1:iprimshell),concoeff(iprimshell-nprmshadd+1:iprimshell))
		
		if (shelltype(ishell)==-1) then !Separate SP shell as S and P shells
			shelltype(ishell)=0 !s
			call renormmoldengau(nprmshadd,shelltype(ishell),primexp(iprimshell-nprmshadd+1:iprimshell),concoeff(iprimshell-nprmshadd+1:iprimshell))
			ishell=ishell+1
			shelltype(ishell)=1 !p
			shellcon(ishell)=shellcon(ishell-1)
			shell2atom(ishell)=shell2atom(ishell-1)
			primexp(iprimshell+1:iprimshell+nprmshadd)=primexp(iprimshellold+1:iprimshell)
			do itmp=1,ncon !Backspace and load P contract coefficient
				backspace(10)
			end do
			do itmp=1,ncon
				read(10,*) exptmp,rnouse,concoefftmp
				if (concoefftmp==0D0) cycle
				iprimshell=iprimshell+1
				concoeff(iprimshell)=concoefftmp
			end do
			call renormmoldengau(nprmshadd,shelltype(ishell),primexp(iprimshell-nprmshadd+1:iprimshell),concoeff(iprimshell-nprmshadd+1:iprimshell))
		end if
		read(10,"(a)") c80
		if (c80==" ") exit
		backspace(10)
	end do
	read(10,"(a)") c80
	if (c80==" ".or.index(c80,"[")/=0) exit !Finished reading [GTO] field
	backspace(10)
end do

!Determine if the basis functions are Cartesian or spherical harmonic type. Admixture Cartesian and spherical type are permitted
isphergau=0 !Default is Cartesian type
i5D=0
i7F=0
i10Flabel=0
i9G=0
i11H=0

if (iorca==1) then !ORCA only uses spherical harmonic type and thus need not to be tested
    i5D=1;i7F=1;i9G=1;i11H=1
else
    imaxL=maxval(shelltype)
    if (infomode==0) write(*,"(' The highest angular moment basis functions is ',a)") shtype2name(imaxL) 
    if (imaxL>=2) then !Containing d or higher basis functions
	    rewind(10)
	    do while(.true.)
		    read(10,"(a)") c80
		    if (index(c80,'[5D')/=0.or.index(c80,'[5d')/=0) i5D=1
		    if (index(c80,'[7F')/=0.or.index(c80,'[7f')/=0) i7F=1
		    if (index(c80,'10F')/=0) i10Flabel=1
		    if (index(c80,'9G')/=0.or.index(c80,'9g')/=0) i9G=1
		    if (index(c80,'11H')/=0.or.index(c80,'11h')/=0) i11H=1
		    if (index(c80,'[MO]')/=0) exit
	    end do
	    if (i5D==1) then !By default, using 5D also implies 7F is used, unless 10F is explicitly specified
		    i7F=1
		    if (i10Flabel==1) i7F=0
	    end if
        if (imaxL==4.and.i7F==1.and.i9G==0) then
            write(*,"(a)") " Note: f functions are harmonic spherical but [9G] label was not found. All g functions are assumed to be harmonic spherical type"
            i9G=1
        end if
        if (imaxL==5.and.i9G==1.and.i11H==0) then
            write(*,"(a)") " Note: g functions are harmonic spherical but [11H] label was not found. All h functions are assumed to be harmonic spherical type"
            i11H=1
        end if
    end if
end if

if (i5D==1.or.i7F==1.or.i9G==1.or.i11H==1) isphergau=1
if (i5D==1) where(shelltype==2) shelltype=-2
if (i7F==1) where(shelltype==3) shelltype=-3
if (i9G==1) where(shelltype==4) shelltype=-4
if (i11H==1) where(shelltype==5) shelltype=-5
if (infomode==0) then
	if (isphergau==0) then
		write(*,*) "All basis functions are Cartesian type"
	else if (isphergau==1) then
		if (i5D==1.and.any(abs(shelltype)==2)) write(*,*) "All D basis functions are spherical harmonic type"
		if (i7F==1.and.any(abs(shelltype)==3)) write(*,*) "All F basis functions are spherical harmonic type"
		if (i9G==1.and.any(abs(shelltype)==4)) write(*,*) "All G basis functions are spherical harmonic type"
		if (i11H==1.and.any(abs(shelltype)==5)) write(*,*) "All H basis functions are spherical harmonic type"
	end if
end if
nbasis=0
do ishell=1,nshell
	nbasis=nbasis+shtype2nbas(shelltype(ishell))
end do

!!!!! Load orbital information. The sequence: Alpha(high occ / low ene) -> Alpha(low occ / high ene) -> Beta(high occ / low ene) -> Beta(low occ / high ene)
!Close shell orbitals are formally marked as "Alpha" spin. For singly occupied orbitals of ROHF, the spin marker are also Alpha
if (infomode==0) write(*,*) "Loading orbitals..."
!call walltime(iwalltime1)
nmo=nbasis
! Here I don't use call loclabel(10,"Beta",ibeta) to check if there are Beta orbitals, because for very large file, this will be quite expensive
! I assume that when the first orbital has occupation number <1.05, then the wavefunction must be unrestricted
ibeta=0
call loclabel(10,"[MO]")
do while(.true.)
	read(10,"(a)") c80
	if (index(c80,"OCCUP=")/=0.or.index(c80,"Occup=")/=0) then
		read(c80,*) c80tmp,occtmp
		if (occtmp<1.05D0) ibeta=1
		exit
	end if
end do
!Allocate size for arrays
if (ibeta==0) then
	nmo=nbasis
	allocate(amocoeff(nmo,nbasis),MOocc(nmo),MOene(nmo),MOtype(nmo),MOsym(nmo))
	amocoeff=0D0
else if (ibeta==1) then
	nmo=2*nbasis
	allocate(amocoeff(nbasis,nbasis),bmocoeff(nbasis,nbasis),MOocc(nmo),MOene(nmo),MOtype(nmo),MOsym(nmo))
	amocoeff=0D0
	bmocoeff=0D0
end if
MOsym=" "
MOocc=0D0
MOene=0D0
!Start to load orbitals
call loclabel(10,"[MO]")
read(10,*)
iMOa=0
iMOb=0
iloadorb=0
do while(.true.)
	iloadorb=iloadorb+1
	read(10,"(a)") c80 !Test if it is "Sym=", some programs do not output this field
	isym=index(c80,"Sym")
	if (isym==0) isym=index(c80,"SYM")
	if (isym/=0) then
		ieq=index(c80,'=')
		read(c80(ieq+1:),*) symtmp
		!Remove digitals before the IRREP, e.g. 23B1 should be changed to B1
		do jtmp=1,len(symtmp)
			if (ichar(symtmp(jtmp:jtmp))<48.or.ichar(symtmp(jtmp:jtmp))>57) exit !Find the first position of non-digital
		end do
		symtmp=symtmp(jtmp:len(symtmp))
	else
		symtmp="?"
		backspace(10)
	end if
	read(10,"(a)") c80 !Read orbital energy. Since there may be no spacing between = and energy (Bagel program), position of = should be tested
	ieq=index(c80,'=')
	read(c80(ieq+1:),*) enetmp
 	!write(*,*) iloadorb,enetmp,nbasis,nmo  !<<------ If encountering problem when loading MOs, using this to locate the problematic MO
    !PS: When loading is failed, it is suggested to search *** in the molden file
	
	read(10,"(a)") c80 !Read orbital spin
	ispintmp=1 !Alpha
	if (index(c80,"Beta")/=0) ispintmp=2 !Beta
	read(10,*) c80,occtmp !Read orbital occupation number
	if (ispintmp==1) then
		iMOa=iMOa+1
		MOocc(iMOa)=occtmp
		MOene(iMOa)=enetmp
        call corrsymlab(MOsym(iMOa),symtmp) !Purify symmetry label to meet standard
		do while(.true.) !This loading way is slow but compatible with Dalton and .molden by all other programs
			read(10,"(a)",iostat=ierror) c80
			if (ierror/=0) then !May reached end of file
				backspace(10)
				exit
			end if
			read(c80,*,iostat=ierror) itmp,tmpval
			if (ierror/=0) then !New orbital
				backspace(10)
				exit
			end if
			amocoeff(iMOa,itmp)=tmpval
		end do
	else
		iMOb=iMOb+1
		MOocc(nbasis+iMOb)=occtmp
		MOene(nbasis+iMOb)=enetmp
        call corrsymlab(MOsym(nbasis+iMOb),symtmp) !Purify symmetry label to meet standard
		do while(.true.)
			read(10,"(a)",iostat=ierror) c80
			if (ierror/=0) then !May read to end of file
				backspace(10)
				exit
			end if
			read(c80,*,iostat=ierror) itmp,tmpval
			if (ierror/=0) then !New orbital
				backspace(10)
				exit
			end if
			bmocoeff(iMOb,itmp)=tmpval
		end do
	end if
	read(10,"(a)",iostat=ierror) c80 !Test if the end of [MO] field has reached
	if (ierror/=0.or.c80==" ".or.index(c80,'[')/=0) exit
	backspace(10)
end do
!call walltime(iwalltime2)
!write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

!Fix orbital coefficients for ORCA. ORCA is rather rather frantic:
!the F(+3,-3) and G(+3,-3,+4,-4) and H(+3,-3,+4,-4) in ORCA are normalized to -1 rather than 1
!therefore the sign of their coefficients in all orbitals must be inverted! Hey ORCA, why did you use so strange convention????? Totally non-understandable!
if (iorca==1) then
	ibasis=0
	do ishell=1,nshell
		if (shelltype(ishell)==-3) then !f
			amocoeff(:,ibasis+6:ibasis+7)=-amocoeff(:,ibasis+6:ibasis+7)
			if (ibeta==1) bmocoeff(:,ibasis+6:ibasis+7)=-bmocoeff(:,ibasis+6:ibasis+7)
		else if (shelltype(ishell)==-4) then !g
			amocoeff(:,ibasis+6:ibasis+9)=-amocoeff(:,ibasis+6:ibasis+9)
			if (ibeta==1) bmocoeff(:,ibasis+6:ibasis+9)=-bmocoeff(:,ibasis+6:ibasis+9)
		else if (shelltype(ishell)==-5) then !h
            amocoeff(:,ibasis+6:ibasis+9)=-amocoeff(:,ibasis+6:ibasis+9)
		end if
		ibasis=ibasis+shtype2nbas(shelltype(ishell))
	end do
end if
!Fix orbital coefficients for CFour according to Molden2aim. CFour only use Cartesian type basis function
!Notice that at current stage the GTO recording sequence has not been reordered, which is still identical to Molden sequence
!Also notice that even we do this, the result is still incorrect for e.g. benzene. But if we don't do this, we can't even obtain correct result for test case of Molden2aim
if (icfour==1) then
	ibasis=0
	do ishell=1,nshell
		if (shelltype(ishell)==2) then
			!xx, yy, zz, xy, xz, yz
			amocoeff(:,ibasis+1:ibasis+3)=amocoeff(:,ibasis+1:ibasis+3)*sqrt(3D0) !d(xx,yy,zz)*sqrt(3)
			if (ibeta==1) bmocoeff(:,ibasis+1:ibasis+3)=bmocoeff(:,ibasis+1:ibasis+3)*sqrt(3D0)
		else if (shelltype(ishell)==3) then
			!xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz !Molden sequence
			amocoeff(:,ibasis+1:ibasis+3)=amocoeff(:,ibasis+1:ibasis+3)*sqrt(15D0) !f(xxx,yyy,zzz)*sqrt(15)
			amocoeff(:,ibasis+4:ibasis+9)=amocoeff(:,ibasis+4:ibasis+9)*sqrt(3D0) !f(xyy,xzz,yxx,yzz,zxx,zyy)*sqrt(3)
			if (ibeta==1) then
				bmocoeff(:,ibasis+1:ibasis+3)=bmocoeff(:,ibasis+1:ibasis+3)*sqrt(15D0)
				bmocoeff(:,ibasis+4:ibasis+9)=bmocoeff(:,ibasis+4:ibasis+9)*sqrt(3D0)
			end if
		else if (shelltype(ishell)==4) then
			!xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy !Molden sequence
			amocoeff(:,ibasis+1:ibasis+3)=amocoeff(:,ibasis+1:ibasis+3)*sqrt(105D0) !g(x4,y4,z4)*sqrt(105)
			amocoeff(:,ibasis+4:ibasis+9)=amocoeff(:,ibasis+4:ibasis+9)*sqrt(15D0) !g(x3y,x3z,y3x,y3z,z3x,z3y)*sqrt(15)
			amocoeff(:,ibasis+10:ibasis+12)=amocoeff(:,ibasis+10:ibasis+12)*3D0 !g(x2y2,x2z2,y2z2)*3.0
			amocoeff(:,ibasis+13:ibasis+15)=amocoeff(:,ibasis+13:ibasis+15)*sqrt(3D0) !g(x2yz,y2xz,z2xy)*sqrt(3)
			if (ibeta==1) then
				bmocoeff(:,ibasis+1:ibasis+3)=bmocoeff(:,ibasis+1:ibasis+3)*sqrt(105D0)
				bmocoeff(:,ibasis+4:ibasis+9)=bmocoeff(:,ibasis+4:ibasis+9)*sqrt(15D0)
				bmocoeff(:,ibasis+10:ibasis+12)=bmocoeff(:,ibasis+10:ibasis+12)*3D0
				bmocoeff(:,ibasis+13:ibasis+15)=bmocoeff(:,ibasis+13:ibasis+15)*sqrt(3D0)
			end if
		end if
		ibasis=ibasis+shtype2nbas(shelltype(ishell))
	end do	
end if
!xtb ignore core orbitals (including subvalence ones), therefore the nuclear charge of atoms must be substracted by number of core electrons
if (ixtb==1) then
	do i=1,ncenter
		if (a(i)%index>2.and.a(i)%index<=10) a(i)%charge=a(i)%charge-2
		if (a(i)%index>=11.and.a(i)%index<=18) a(i)%charge=a(i)%charge-10 !Na~Ar
		if (a(i)%index>=19.and.a(i)%index<=29) a(i)%charge=a(i)%charge-18 !K~Cu
		if (a(i)%index>=30.and.a(i)%index<=36) a(i)%charge=a(i)%charge-28 !Zn~Kr
		if (a(i)%index>=37.and.a(i)%index<=47) a(i)%charge=a(i)%charge-36 !Rb~Ag
		if (a(i)%index>=48.and.a(i)%index<=54) a(i)%charge=a(i)%charge-46 !Cd~Xe
		if (a(i)%index==55.and.a(i)%index==56) a(i)%charge=a(i)%charge-54 !Cs~Ba
		if (a(i)%index>=57.and.a(i)%index<=71) a(i)%charge=3 !La~Lu. I have confirmed this by practical calculation
		if (a(i)%index>=72.and.a(i)%index<=79) a(i)%charge=a(i)%charge-68 !Hf~Au
		if (a(i)%index>=80.and.a(i)%index<=86) a(i)%charge=a(i)%charge-78 !Hg~Rn
        !GFN-xTB only supports up to 86
	end do
end if

!Determine wavefunction type
if (ibeta==0) then
	MOtype=0 !Close shell orbital
	wfntype=0 !RHF
	if (any(MOocc/=nint(MOocc))) then
		wfntype=3 !R-post-HF
	else if (any(MOocc==1D0)) then
		wfntype=2 !ROHF
		do imo=1,nmo
			if (MOocc(imo)==1D0) MOtype(imo)=1
		end do
	end if
	if (infomode==0) write(*,"( ' The actual number of orbitals read:',i10)") iMOa
else if (ibeta==1) then
	wfntype=1 !UHF
	if (any(MOocc/=nint(MOocc))) wfntype=4 !U-post-HF
	MOtype(1:nbasis)=1
	MOtype(nbasis+1:nmo)=2
	if (infomode==0) write(*,"( ' The actual number of Alpha/Beta orbitals read:',i10,'  /',i10)") iMOa,iMOb
end if
call updatenelec !Cound the number of electrons

close(10)

!!!!!! All reading have finished, now generate basis information
!Below codes are adapted from readfch

!Backup spherical Gaussian basis information with 5D suffix (of course, may be 7f, 9g... in fact), &
!convert them to Cartesian type temporarily, at final stage recover them back
if (isphergau==1) then
	allocate(shelltype5D(nshell))
	shelltype5D=shelltype
	where (shelltype<=-2) shelltype=-shelltype !Convert to Cartesian type
	nbasis5D=nbasis
	nbasis=0
	do i=1,nshell
		nbasis=nbasis+shtype2nbas(shelltype(i))
	end do
end if
allocate(shtypeCar(nbasis)) !Store shell information of Cartesian basis into global array, which may be used later
shtypeCar=shelltype
nbasisCar=nbasis

!Allocate space for arrays
nprims=0
do i=1,nshell
	nprims=nprims+shtype2nbas(shelltype(i))*shellcon(i)
end do
allocate(b(nprims),CO(nmo,nprims),basshell(nbasis),bascen(nbasis),bastype(nbasis),primstart(nbasis),&
primend(nbasis),primconnorm(nprims),basstart(ncenter),basend(ncenter))

!Fill CObasa and CObasb
if (isphergau==0) then
	allocate(CObasa(nbasis,nbasis))
	CObasa=transpose(amocoeff)
	if (wfntype==1.or.wfntype==4) then
		allocate(CObasb(nbasis,nbasis))
		CObasb=transpose(bmocoeff)
	end if
else if (isphergau==1) then !Since we have artifically converted spherical shells to Cartesian shells, here the orbital coefficients are also correspondingly converted
	allocate(CObasa(nbasis,nbasis),CObasa5D(nbasis5D,nbasis5D))
	CObasa5D=transpose(amocoeff)
	CObasa=0D0
	if (wfntype==1.or.wfntype==4) then
		allocate(CObasb(nbasis,nbasis),CObasb5D(nbasis5D,nbasis5D))
		CObasb5D=transpose(bmocoeff)
		CObasb=0D0
	end if
	!Map 5D coefficient to 6D coefficient. Since the number of spherical basis functions is more than Cartesian ones, 
	!therefore CObasa (6D) will have some orbitals with vacant coefficients, only orbitals (1~nbasis5D) are filled
	ipos5D=1
	ipos6D=1
	do ish=1,nshell
		ishtyp5D=shelltype5D(ish)
		ishtyp6D=shelltype(ish)
		numshorb5D=shtype2nbas(ishtyp5D)
		numshorb6D=shtype2nbas(ishtyp6D)
		if (ishtyp5D>=-1) then !S or P or SP or other Cartesian shells, directly copy
			CObasa(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasa5D(ipos5D:ipos5D+numshorb5D-1,:)
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasb5D(ipos5D:ipos5D+numshorb5D-1,:)			
		else if (ishtyp5D==-2) then
			!5D->6D
			CObasa(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasb5D(ipos5D:ipos5D+4,:))
		else if (ishtyp5D==-3) then
			!7F->10F
			CObasa(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasb5D(ipos5D:ipos5D+6,:))
		else if (ishtyp5D==-4) then
			!9G->15G
			CObasa(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasb5D(ipos5D:ipos5D+8,:))
		else if (ishtyp5D==-5) then
			!11H->21H
			CObasa(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,:))
			if (wfntype==1.or.wfntype==4) CObasb(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasb5D(ipos5D:ipos5D+10,:))
		end if
		ipos5D=ipos5D+numshorb5D
		ipos6D=ipos6D+numshorb6D
	end do
end if

if (infomode==0) write(*,*) "Converting basis function information to GTF information..."
!Distribute exponent, functype to every GTF and generate CO(:,:) from amocoeff/bmocoeff
!Fill: b,basshell,bascen,bastype,co,primstart,primend,primconnorm
k=1 !current position of GTF
iexp=1
ibasis=1 !current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !cycle each basis shell
	b(k:k+shellcon(i)*shtype2nbas(shelltype(i))-1)%center=shell2atom(i)
	basshell(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=i !!! set basis attributed to which shell
	bascen(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=shell2atom(i) !!! set basis attributed to which center
	do j=1,shtype2nbas(shelltype(i)) !cycle each basis function in each basis shell
		b(k:k+shellcon(i)-1)%type=s2f(shelltype(i),j)
		bastype(ibasis)=s2f(shelltype(i),j) !!! set basis type
		primstart(ibasis)=k !!! From where the GTFs attributed to ibasis'th basis
		primend(ibasis)=k+shellcon(i)-1 !!! To where the GTFs attributed to ibasis'th basis
		do l=1,shellcon(i) !cycle each GTF in each basis function
			b(k)%exp=primexp(iexp+l-1)
			tnormgau=normgau(b(k)%type,b(k)%exp)  !!!Normalization coefficient of Cartesian GTFs
			temp=concoeff(iexp+l-1)  !!!Contraction coefficient of GTFs
			primconnorm(k)=temp*tnormgau !Combines contraction and normalization coefficient
			do imo=1,nmo
				if (wfntype==0.or.wfntype==2.or.wfntype==3) then
					CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
! 					if (imo==4) write(*,"(4f18.10)") CO(imo,k),CObasa(ibasis,imo),temp,tnormgau
				else if (wfntype==1.or.wfntype==4) then
					if (isphergau==1) then
						if (imo<=nbasis5D) CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
						if (imo>nbasis5D) CO(imo,k)=CObasb(ibasis,imo-nbasis5D)*temp*tnormgau
					else
						if (imo<=nbasis) CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
						if (imo>nbasis) CO(imo,k)=CObasb(ibasis,imo-nbasis)*temp*tnormgau
					end if
				end if
			end do
			k=k+1
		end do
		ibasis=ibasis+1
	end do
	iexp=iexp+shellcon(i)
end do

!Check normalizaiton of basis functions
! do i=1,size(sbas,1)
! 	write(*,"(i10,f12.6)") i,sbas(i,i)
! end do

if (isphergau==1) then
	if (iloadasCart==1) then !For special purpose, keep Cartesian basis functions, e.g. convert spherical .fch/.molden to .47 file
		!Currently nbasis and dimension of all matrix (except for CO) correspond to full Cartesian case, but nmo &
		!and related arrays as well as CO still correspond to spherical harmonic case and thus need to be "expanded", the MO energies are artifically set to 100
		allocate(MOocc5D(nmo),MOene5D(nmo),MOtype5D(nmo),MOsym5D(nmo),CO5D(nmo,nprims))
		MOocc5D=MOocc
		MOene5D=MOene
		MOtype5D=MOtype
		MOsym5D=MOsym
		CO5D=CO
		deallocate(MOocc,MOene,MOtype,MOsym,CO)
		if (wfntype==0.or.wfntype==2.or.wfntype==3) nmo=nbasis !R, RO
		if (wfntype==1.or.wfntype==4) nmo=2*nbasis !U
		allocate(MOocc(nmo),MOene(nmo),MOtype(nmo),MOsym(nmo),CO(nmo,nprims))
		MOocc=0
		MOene=100
		MOsym=" "
		CO=0
		if (wfntype==0.or.wfntype==2.or.wfntype==3) then !R, RO
			MOtype=0
			MOocc(1:nbasis5D)=MOocc5D
			MOene(1:nbasis5D)=MOene5D
			MOtype(1:nbasis5D)=MOtype5D
			MOsym(1:nbasis5D)=MOsym5D
			CO(1:nbasis5D,:)=CO5D
		else !U
			MOtype(:nbasis)=1
			MOtype(nbasis+1:)=2
			MOocc(:nbasis5D)=MOocc5D(:nbasis5D)
			MOocc(nbasis+1:nbasis+nbasis5D)=MOocc5D(nbasis5D+1:)
			MOene(:nbasis5D)=MOene5D(:nbasis5D)
			MOene(nbasis+1:nbasis+nbasis5D)=MOene5D(nbasis5D+1:)
			MOtype(:nbasis5D)=MOtype5D(:nbasis5D)
			MOtype(nbasis+1:nbasis+nbasis5D)=MOtype5D(nbasis5D+1:)
			MOsym(:nbasis5D)=MOsym5D(:nbasis5D)
			MOsym(nbasis+1:nbasis+nbasis5D)=MOsym5D(nbasis5D+1:)
			CO(:nbasis5D,:)=CO5D(:nbasis5D,:)
			CO(nbasis+1:nbasis+nbasis5D,:)=CO5D(nbasis5D+1:,:)
		end if
		isphergau=0
		
	else !Commonly case, transform to spherical harmonic functions
		if (infomode==0) write(*,*) "Back converting basis function information from Cartesian to spherical type..."
		!Recover spherical Gaussian basis function information
		nbasis=nbasis5D
		shelltype=shelltype5D
		ibasis=1
		do i=1,nshell
			basshell(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=i
			bascen(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=shell2atom(i)
			do j=1,shtype2nbas(shelltype(i))
				bastype(ibasis)=s2f(shelltype(i),j)
				ibasis=ibasis+1
			end do
		end do
		deallocate(CObasa)
		allocate(CObasa(nbasis,nbasis))
		CObasa=CObasa5D
		if (wfntype==1.or.wfntype==4) then
			deallocate(CObasb)
			allocate(CObasb(nbasis,nbasis))
			CObasb=CObasb5D
		end if
	end if
end if

!Move local shell arrays to the global ones, they will be used in other functions
allocate(shtype(nshell),shcen(nshell),shcon(nshell),primshexp(nprimshell),primshcoeff(nprimshell))
shtype=shelltype
shcen=shell2atom
shcon=shellcon
primshexp=primexp
primshcoeff=concoeff

!Generate basstart and basend
call bascen2basstart_end

!Generate one-particle density matrix for basis functions
if (igenP==1) then
	if (infomode==0) write(*,*) "Generating density matrix..."
	call genP
end if
if (ifPBC==0) then !For PBC case, Sbas will be generated when need it because it may be relatively expensive
    if (infomode==0) write(*,*) "Generating overlap matrix..."
    call genSbas_curr
end if

!Check wavefunction sanity for non-PBC case
if (infomode==0.and.iorca==0.and.ifPBC==0) then !For ORCA with angular moment >f, warning has already been shown before
	devtmp=abs(sum(Sbas*Ptot)-nint(nelec))
	if (devtmp>0.01D0) then
		write(*,"( ' Deviation of Tr(S*P) to the total number of electrons:',f12.6)") devtmp
		write(*,"(/,a)") " Warning: The wavefunction loaded is problematic! That means this Molden input file cannot be well supported by Multiwfn"
        write(*,"(a)") " If the input file was generated by CP2K, please carefully check Section 2.9.2.1 of Multiwfn manual on how to properly modify this file. &
        For other cases, perhaps the molden file is non-standard, try to use Molden2AIM program to standardize this file before loading it into Multiwfn, see Section 5.1 of Multiwfn manual for detail"
		!write(*,"(a)") " If you really want to proceed, press ENTER button, but notice that the result will not be correct"
		!read(*,*)
	end if
end if

!Output summary of present wavefunction
if (infomode==0) then
	write(*,*)
	write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
	write(*,"(' Net charge:',f12.5,'      Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
	write(*,"(' Atoms:',i7,',  Basis functions:',i7,',  GTFs:',i7)") ncenter,nbasis,nprims
	if (wfntype==0) then
		write(*,"(' This is a restricted single-determinant wavefunction')")
		write(*,"(' Orbitals from 1 to',i6,' are occupied')") nint(nelec/2)
	else if (wfntype==1) then
		write(*,"(' This is an unrestricted single-determinant wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are occupied')") 1,nbasis,1,nint(naelec)
		write(*,"(' Orbitals from ',i6,' to',i6,' are beta,  from',i6,' to',i6,' are occupied')") nbasis+1,nmo,nbasis+1,nbasis+nint(nbelec)
	else if (wfntype==2) then
		write(*,"(' This is a restricted open-shell wavefunction')")
		write(*,"(' Orbitals from',i6,' to',i6,' are doubly occupied')") 1,nint(nbelec)
		write(*,"(' Orbitals from',i6,' to',i6,' are singly occupied')") nint(nbelec)+1,nint(naelec)
	else if (wfntype==3) then
		write(*,"(' This is a restricted multiconfiguration wavefunction')")
	else if (wfntype==4) then
		write(*,"(' This is an unrestricted multiconfiguration wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are beta')") 1,nbasis,nbasis+1,nmo
	end if
end if

call getHOMOidx !Find out index of HOMO, will be used in some cases

tmpnet=sum(a%charge)-nelec
if (infomode==0.and.tmpnet>18.and.any(a%index>18)) then
    write(*,"(/,a)") " !! Warning! Warning! Warning! Warning! Warning! Warning! Warning! Warning !!"
    write(*,"(a,i5,a)") " The net charge of this system is quite large (",nint(tmpnet),")! Probably ECP is employed while you &
    forgot to modify the atomic indices in [atoms] field of the .molden file to actual nuclear charges, in this case some &
    analysis results will be problematic! Please check ""Molden"" part of Section 2.5 of Multiwfn manual to understand why and how to modify the file"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
end if

end subroutine


!--- Purify symmetry label to meet Multiwfn standard, used by readmolden. e.g. 61a1g --> A1g
subroutine corrsymlab(symnew,symorg)
use util
implicit real*8 (a-h,o-z)
character(len=*) symorg,symnew
symorg=adjustl(symorg)
do i=1,len_trim(symorg)
    itmp=ichar(symorg(i:i))
    if ((itmp>=65.and.itmp<=90).or.(itmp>=97.and.itmp<=122)) exit
end do
symnew=" "
symnew=trim(symorg(i:))
itmp=ichar(symnew(1:1))
if (itmp>=97.and.itmp<=122) call lc2uc(symnew(1:1)) !To upper case
do i=2,len_trim(symnew)
    itmp=ichar(symnew(i:i))
    if (itmp>=65.and.itmp<=90) call uc2lc(symnew(i:i)) !To lower case
end do
end subroutine




!---------- Read GAMESS-US and Firefly output file to get GTF and basis information, the suffix must be "gms"
!GAMESS-US always print LCAO cofficients as Cartesian basis functions, while the number of MOs correspond to spherical harmonic functions (if used)
!infomode=0 means output info, =1 silent
subroutine readgms(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
character c80*80,chartmp
integer :: s2f(0:5,21)=0 !Give shell type & orbital index to get functype
real*8,allocatable :: amocoeff(:,:),bmocoeff(:,:)
integer,allocatable :: zcore(:) !The number of core electrons replaced with ECP
real*8,external :: normgau

ifiletype=10
imodwfn=0
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
!---------- The sequence of f functions in Multiwfn (=wfn=wfx) is not identical to GAMESS, so convert here
!11  12  13  14  15  16  17  18  19  20  !Multiwfn sequence
!XXX YYY ZZZ XXY XXZ YYZ XYY XZZ YZZ XYZ
!xxx yyy zzz xxy xxz xyy yyz zzx zzy xyz !GAMESS sequence
s2f(3,1:10)=(/ 11,12,13,14,15,17,16,18,19,20 /)
!---------- The sequence of g functions in Multiwfn (=fch) is not identical to GAMESS(=molden), so convert here
! 21   22   23   24   25   26   27   28   29   30   31   32   33   34   35  !Multiwfn sequence
!ZZZZ YZZZ YYZZ YYYZ YYYY XZZZ XYZZ XYYZ XYYY XXZZ XXYZ XXYY XXXZ XXXY XXXX
!xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy !GAMESS sequence
s2f(4,1:15)=(/ 35,25,21,34,33,29,24,26,22,32,30,23,31,28,27 /)
!---------- The sequence of h functions in Multiwfn (=fch=wfx) is not identical to GAMESS, so convert here
! 36    37    38    39    40    41    42    43    44    45    46  !Multiwfn sequence
!ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ XYYYY
!xxxxx yyyyy zzzzz xxxxy xxxxz xyyyy yyyyz xzzzz yzzzz xxxyy xxxzz !GAMESS sequence
! 47    48    49    50    51    52    53    54    55    56  !Multiwfn sequence
!XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
!xxyyy yyyzz xxzzz yyzzz xxxyz xyyyz xyzzz xxyyz xxyzz xyyzz !GAMESS sequence
s2f(5,1:21)=(/ 56,41,36,55,54,46,40,42,37,53,51, 50,39,47,38,52,45,43,49,48,44 /)

open(10,file=name,status="old")
call loclabel(10,"Firefly Project",ifirefly,maxline=100)
if (ifirefly==1) write(*,*) "Note: This file will be recognized as Firefly output"
if (infomode==0) write(*,*) "Loading various information of the wavefunction"

!!!!! Sanity check
call loclabel(10,"RUNTYP=OPTIMIZE",iopt,maxline=10000)
if (iopt==1) then
	!Although GAMESS-US output final geometry and corresponding wavefunction (labelled by "MOLECULAR ORBITALS"), 
	!the number of orbitals is much smaller than expected, therefore unable to provide enough information
	write(*,"(a)") " Warning: This is an optimization task, only wavefunction corresponding to initial geometry will be loaded"
	write(*,*) "Press ENTER button to continue"
	read(*,*)
end if
!If point group is not C1, only basis set definition of symmetry unique atoms is printed
call loclabel(10,"THE POINT GROUP OF THE MOLECULE IS C1",iC1,maxline=1000)
if (iC1==0) then
    write(*,"(a)") " Error: The point group is not C1, Multiwfn does not support this case"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
    stop
end if

!!!!! Load atom information
call loclabel(10,"TOTAL NUMBER OF ATOMS",ifound)
if (ifirefly==0) read(10,"(47x,i5)") ncenter
if (ifirefly==1) read(10,"(38x,i5)") ncenter
allocate(a(ncenter))
call loclabel(10,"ATOM      ATOMIC",ifound)
read(10,*)
read(10,*)
do iatm=1,ncenter
	read(10,*) c80,a(iatm)%index,a(iatm)%x,a(iatm)%y,a(iatm)%z
end do
a%name=ind2name(a%index)
a%charge=a%index
rewind(10)

!!!!! Load basis set, and build up GTF information
!First time, we count how many shells are there to allocate proper size of allocatable arrays
if (ifirefly==0) call loclabel(10,"SHELL TYPE  PRIMITIVE",ifound)
if (ifirefly==1) call loclabel(10,"SHELL TYPE PRIM",ifound)
if (ifound==0) then
	write(*,"(' Error: Unable to find basis set definition section!')")
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	stop
else
	if (infomode==0) write(*,*) "Loading basis set definition..."
end if
read(10,*)
read(10,*)
iaddshell=0
nshell=0
nprimshell=0
nshell_tmp=0 !Number of basis shells without splitting SP shell
do iatm=1,ncenter
	read(10,*) !Atom name
	read(10,*) !Space line
	do while(.true.)
		read(10,"(a)") c80
		if (c80==" ") then !Finished loading last shell
			nshell=nshell+1
            nshell_tmp=nshell_tmp+1
			if (iaddshell==1) nshell=nshell+1 !Last shell is L
			iaddshell=0
			read(10,"(a)") c80;backspace(10)
			if (iachar(c80(2:2))>=65.or.c80(2:2)=='*') exit !The second column of next line is a letter (element name (* for Bq) or string), indicating that reading present atom has finished
			if (c80==" ") exit !The next line is blank (may occur in the case of Firefly), indicating the whole field of basis set definition has finished
		else
            if (iachar(c80(2:2))>=65.or.c80(2:2)=='*') then !The second column of next line is a letter (element name (* for Bq) or string), indicating that reading present atom has finished
                backspace(10)
                exit
            end if
			nprimshell=nprimshell+1
			if (index(c80,"L")/=0) then
				nprimshell=nprimshell+1 !sp shell will be separated as s and p
				iaddshell=1
			end if
		end if
	end do
    !write(*,*) iatm,nprimshell,nshell,nshell_tmp
end do
write(*,"(' Number of primitive shells (L is splitted as S and P):',i8)") nprimshell
write(*,"(' Number of basis function shells (L is splitted as S and P):',i8)") nshell
write(*,"(' Number of basis function shells (without splitting L):',i8)") nshell_tmp

!Second time, read basis set information actually
allocate(shtype(nshell),shcon(nshell),shcen(nshell))
allocate(primshexp(nprimshell),primshcoeff(nprimshell))
if (ifirefly==0) call loclabel(10,"SHELL TYPE  PRIMITIVE",ifound)
if (ifirefly==1) call loclabel(10,"SHELL TYPE PRIM",ifound)
read(10,*)
read(10,*)
shcon=0
ishell=1
iprimshell=0
do iatm=1,ncenter
    !write(*,"(' Loading atom',i6,' (',a,')')") iatm,ind2name(a(iatm)%index)
	read(10,"(a)") c80 !Atom name
	read(10,"(a)") c80
	do while(.true.)
		read(10,"(a)") c80
		if (c80==" ") then !Finished loading last basis shell
			if (shtype(ishell)==-1) then !Separate SP shell as S and P shells
				shtype(ishell)=0 !s
				ishell=ishell+1
				shtype(ishell)=1 !p
				shcon(ishell)=shcon(ishell-1)
				shcen(ishell)=shcen(ishell-1)
				do itmp=1,shcon(ishell)+1 !Backspace and load P contraction coefficient
					backspace(10)
				end do
				do itmp=1,shcon(ishell)
					iprimshell=iprimshell+1
					if (ifirefly==0) then
						read(10,*) inouse,chartmp,inouse,primshexp(iprimshell),rnouse,primshcoeff(iprimshell)
					else
						read(10,"(a)") c80
						read(c80,*) inouse,chartmp,inouse,primshexp(iprimshell)
						read(c80(69:78),*) primshcoeff(iprimshell)
					end if
				end do
				read(10,*)
			end if
			ishell=ishell+1
			read(10,"(a)") c80;backspace(10)
			if (iachar(c80(2:2))>=65.or.c80(2:2)=='*') exit !The second column of next line is a letter (element name (* for Bq) or string), indicating that reading present atom has finished
			if (c80==" ") exit !The next line is blank (may occur in the case of Firefly), indicating the whole field of basis set definition has finished
		else !Loading next primitive shell in current basis function shell
            !In very special case, there are two Bq successively occur at the end, below is an example. &
            !In this case Bq doesn't have any line defining primitive shell, so we need to exit according to the second character
            !...
            ! 84   L     251             1.1695961    0.399512826089    0.607683718598
            ! 84   L     252             0.3803890    0.700115468880    0.391957393099
            !
            ! *         
            !
            ! *         
            !
            ! TOTAL NUMBER OF BASIS SET SHELLS             =   84
            !...
            if (iachar(c80(2:2))>=65.or.c80(2:2)=='*') then !The second column of next line is a letter (element name (* for Bq) or string), indicating that reading present atom has finished
                backspace(10)
                exit
            end if
			iprimshell=iprimshell+1
			shcon(ishell)=shcon(ishell)+1
			if (ifirefly==0) then !GAMESS-US
				read(c80,*) inouse,chartmp,inouse,primshexp(iprimshell),primshcoeff(iprimshell)
			else
				read(c80,*) inouse,chartmp,inouse,primshexp(iprimshell)
				read(c80(43:52),*) primshcoeff(iprimshell)
! 				write(*,*) iprimshell,chartmp,primshcoeff(iprimshell)
			end if
			shcen(ishell)=iatm
			if (chartmp=="L") then
				shtype(ishell)=-1
			else if (chartmp=="S") then
				shtype(ishell)=0
			else if (chartmp=="P") then
				shtype(ishell)=1
			else if (chartmp=="D") then
				shtype(ishell)=2
			else if (chartmp=="F") then
				shtype(ishell)=3
			else if (chartmp=="G") then
				shtype(ishell)=4
			else if (chartmp=="H") then
				shtype(ishell)=5
			end if
		end if
	end do
end do
allocate(shtypeCar(nbasis)) !Store shell information of Cartesian basis into global array, which may be used later
shtypeCar=shtype

nbasis=0
do ishell=1,nshell
	nbasis=nbasis+shtype2nbas(shtype(ishell))
end do
nbasisCar=nbasis

if (infomode==0) write(*,*) "Loading orbitals..."

!  "nmoactual" below is the actual number of MOs (for each spin in unrestricted case), less or equal to nbasis
!  In GAMESS-US, the actual number of orbitals is nbasis - A - B, where
!A is "NUMBER OF SPHERICAL CONTAMINANTS DROPPED", namely the difference betweem number of spherical and cartesian basis functions
!B is "NUMBER OF LINEARLY DEPENDENT MOS DROPPED", namely the number of removed linearly dependent basis functions
!  In Firefly, the actual number of orbitals is nbasis - "TOTAL NUMBER OF CONTAMINANTS DROPPED"
if (ifirefly==0) then
	call loclabel(10,"TOTAL NUMBER OF MOS IN VARIATION SPACE=",ispher)
	if (ispher==0) then !Cartesian functions
		nmoactual=nbasis
	else
		read(10,"(40x,i8)") nmoactual
	end if
else
	nremove=0
	call loclabel(10,"TOTAL NUMBER OF CONTAMINANTS DROPPED",ifound)
    read(10,"(a)") c80
	if (ifound==1) call readaftersign_int(c80,'=',nremove)
    nmoactual=nbasis-nremove
end if
call loclabel(10,"----- BETA SET -----",ibeta)
if (ibeta==0) then !Only one set of orbitals
	nmo=nbasis
	allocate(amocoeff(nbasis,nmoactual),MOocc(nmo),MOene(nmo),MOtype(nmo),MOsym(nmo))
	MOocc=0
	MOene=100
	MOsym="?"
	call loclabel(10,"          EIGENVECTORS",ifound)
    if (ifound==0) then
        write(*,"(a)") " Error: Unable to find orbital coefficients field! The NPRINT parameter must be set to default value"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        stop
    end if
	call readgmsLCAO(10,ifirefly,nbasis,nmoactual,amocoeff,MOene(1:nmoactual),MOsym(1:nmoactual))
else
	nmo=2*nbasis
	allocate(amocoeff(nbasis,nmoactual),bmocoeff(nbasis,nmoactual),MOocc(nmo),MOene(nmo),MOtype(nmo),MOsym(nmo))
	MOocc=0
	MOene=100
	MOsym="?"
	call loclabel(10,"          EIGENVECTORS",ifound)
    if (ifound==0) then
        write(*,"(a)") " Error: Unable to find orbital coefficients field! The NPRINT parameter must be set to default value"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        stop
    end if
	call readgmsLCAO(10,ifirefly,nbasis,nmoactual,amocoeff,MOene(1:nmoactual),MOsym(1:nmoactual))
	call loclabel(10,"          EIGENVECTORS",ifound,0) !Don't rewind
	call readgmsLCAO(10,ifirefly,nbasis,nmoactual,bmocoeff,MOene(nbasis+1:nbasis+nmoactual),MOsym(nbasis+1:nbasis+nmoactual))
end if

!Determine the number of electrons and load the number of core electrons replaced with ECP
call loclabel(10,"NUMBER OF OCCUPIED ORBITALS (ALPHA)",ifound)
if (ifirefly==0) read(10,"(47x,f5.0)") naelec
if (ifirefly==1) read(10,"(38x,f5.0)") naelec
call loclabel(10,"NUMBER OF OCCUPIED ORBITALS (BETA )",ifound)
if (ifirefly==0) read(10,"(47x,f5.0)") nbelec
if (ifirefly==1) read(10,"(38x,f5.0)") nbelec
call loclabel(10,"THE ECP RUN REMOVES",iECP,maxline=10000) !ECP is used
if (iECP==1) then
	allocate(zcore(ncenter))
	zcore=0
	rewind(10)
	do while(.true.)
		call loclabel(10,"WITH ZCORE",ifound,0,maxline=10000)
		if (ifound==1) then
			read(10,"(a)") c80
			itmp=index(c80,"ATOM")
			read(c80(itmp+4:),*) iatm
			itmp=index(c80,"ZCORE")
			read(c80(itmp+5:),*) zcore(iatm)
		else
			exit
		end if
	end do
	rewind(10)
	do while(.true.)
		call loclabel(10,"ARE THE SAME AS ATOM",ifound,0,maxline=10000)
		if (ifound==1) then
			read(10,"(a)") c80
			itmp=index(c80,"ATOM")
			read(c80(itmp+4:),*) iatm
			itmp=index(c80,"ATOM",back=.true.)
			read(c80(itmp+4:),*) jatm
			zcore(iatm)=zcore(jatm)
		else
			exit
		end if
	end do
	a%charge=a%charge-zcore !Change "charge" from atom indices to actual number of electrons in ECP calculation
	naelec=naelec-sum(zcore)/2
	nbelec=nbelec-sum(zcore)/2
end if
nelec=naelec+nbelec

!Determine wavefunction type
if (ibeta==0) then
	if (any(MOocc/=nint(MOocc))) then
! 		wfntype=3 !R-post-HF
	else if (naelec/=nbelec) then
		wfntype=2 !ROHF
		MOocc(1:nint(nbelec))=2D0
		MOocc(nint(nbelec+1):nint(naelec))=1D0
		MOtype=0 !Close shell orbital
		MOtype(nint(nbelec+1):nint(naelec))=1
	else !RHF
		wfntype=0
		MOocc(1:nint(naelec))=2D0
		MOtype=0 !Close shell orbital		
	end if
else if (ibeta==1) then
	wfntype=1 !UHF
! 	if (any(MOocc/=nint(MOocc))) wfntype=4 !U-post-HF
	MOtype(1:nbasis)=1
	MOtype(nbasis+1:nmo)=2
	MOocc(1:nint(naelec))=1D0
	MOocc(nbasis+1:nbasis+nint(nbelec))=1D0
end if

close(10)

!!!!!! All reading have finished, now generate basis information
 
!Allocate space for arrays
nprims=0
do i=1,nshell
	nprims=nprims+shtype2nbas(shtype(i))*shcon(i)
end do
allocate(b(nprims),CO(nmo,nprims),basshell(nbasis),bascen(nbasis),bastype(nbasis),primstart(nbasis),&
primend(nbasis),primconnorm(nprims),basstart(ncenter),basend(ncenter))

!Fill CObasa and CObasb, the gap spaces due to difference between Cartesian and spherical harmonic functions are filled by zero
allocate(CObasa(nbasis,nbasis))
CObasa=0
CObasa(:,1:nmoactual)=amocoeff
if (wfntype==1.or.wfntype==4) then
	allocate(CObasb(nbasis,nbasis))
	CObasb=0
	CObasb(:,1:nmoactual)=bmocoeff
end if

if (infomode==0) write(*,*) "Converting basis function information to GTF information..."
!Distribute exponent, functype to every GTF and generate CO(:,:) from amocoeff/bmocoeff
!Fill: b,basshell,bascen,bastype,co,primstart,primend,primconnorm
k=1 !current position of GTF
iexp=1
ibasis=1 !current position of basis
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !cycle each basis shell
	b(k:k+shcon(i)*shtype2nbas(shtype(i))-1)%center=shcen(i)
	basshell(ibasis:ibasis+shtype2nbas(shtype(i))-1)=i !!! set basis attributed to which shell
	bascen(ibasis:ibasis+shtype2nbas(shtype(i))-1)=shcen(i) !!! set basis attributed to which center
	do j=1,shtype2nbas(shtype(i)) !cycle each basis function in each basis shell
		b(k:k+shcon(i)-1)%type=s2f(shtype(i),j)
		bastype(ibasis)=s2f(shtype(i),j) !!! set basis type
		primstart(ibasis)=k !!! From where the GTFs attributed to ibasis'th basis
		primend(ibasis)=k+shcon(i)-1 !!! To where the GTFs attributed to ibasis'th basis
		do l=1,shcon(i) !cycle each GTF in each basis function
			b(k)%exp=primshexp(iexp+l-1)
			tnormgau=normgau(b(k)%type,b(k)%exp)  !!!Normalization coefficient of Cartesian GTFs
			temp=primshcoeff(iexp+l-1)  !!!Contraction coefficient of GTFs
			primconnorm(k)=temp*tnormgau !Combines contraction and normalization coefficient
			do imo=1,nmo
				if (wfntype==0.or.wfntype==2.or.wfntype==3) then
					CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
				else if (wfntype==1.or.wfntype==4) then
					if (imo<=nbasis) CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
					if (imo>nbasis) CO(imo,k)=CObasb(ibasis,imo-nbasis)*temp*tnormgau
				end if
			end do
			k=k+1
		end do
		ibasis=ibasis+1
	end do
	iexp=iexp+shcon(i)
end do

!Generate basstart and basend
call bascen2basstart_end

!Generate one-particle matrices for basis functions
if (igenP==1) then
	if (infomode==0) write(*,*) "Generating density matrix..."
	call genP
end if
if (infomode==0) write(*,*) "Generating overlap matrix..."
call genSbas_curr

!Check wavefunction sanity
if (infomode==0) then
	devtmp=abs(sum(Sbas*Ptot)-nint(nelec))
	if (devtmp>0.01D0) then
		write(*,"( ' Deviation of Tr(S*P) to the total number of electrons:',f12.6)") devtmp
		write(*,"(/,a)") " Warning: The wavefunction loaded is problematic! That means this Molden input file cannot be well supported by Multiwfn"
	end if
end if
 
!Output summary of present wavefunction
if (infomode==0) then
	write(*,*)
	write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
	write(*,"(' Net charge:',f12.5,'      Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
	write(*,"(' Atoms:',i7,',  Basis functions:',i7,',  GTFs:',i7)") ncenter,nbasis,nprims
	if (wfntype==0) then
		write(*,"(' This is a restricted single-determinant wavefunction')")
		write(*,"(' Orbitals from 1 to',i6,' are occupied')") nint(nelec/2)
	else if (wfntype==1) then
		write(*,"(' This is an unrestricted single-determinant wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are occupied')") 1,nbasis,1,nint(naelec)
		write(*,"(' Orbitals from ',i6,' to',i6,' are beta,  from',i6,' to',i6,' are occupied')") nbasis+1,nmo,nbasis+1,nbasis+nint(nbelec)
	else if (wfntype==2) then
		write(*,"(' This is a restricted open-shell wavefunction')")
		write(*,"(' Orbitals from',i6,' to',i6,' are doubly occupied')") 1,nint(nbelec)
		write(*,"(' Orbitals from',i6,' to',i6,' are singly occupied')") nint(nbelec)+1,nint(naelec)
	else if (wfntype==3) then
		write(*,"(' This is a restricted multiconfiguration wavefunction')")
	else if (wfntype==4) then
		write(*,"(' This is an unrestricted multiconfiguration wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are beta')") 1,nbasis,nbasis+1,nmo
	end if
end if

call getHOMOidx !Find out index of HOMO, will be used in some cases
end subroutine

!----- A routine used to read GAMESS-US LCAO matrix. One should first move pointer to the line containing "EIGENVECTORS"
subroutine readgmsLCAO(fileid,ifirefly,n1,n2,mat,ene,sym)
implicit real*8 (a-h,o-z)
integer fileid,ifirefly,n1,n2
real*8 :: mat(n1,n2),ene(n2)
character(len=4) :: sym(n2)
character c80tmp*80
read(fileid,*)
read(fileid,*)
ncol=5
nt=ceiling(n2/float(ncol))
do i=1,nt !Number of frames
! 	write(*,*) i !Show frame
	read(fileid,*)
	ns=(i-1)*ncol+1
	if (i/=nt) ne=(i-1)*ncol+ncol
	if (i==nt) ne=n2
	read(fileid,*)
    
 	!read(fileid,"(a)") c80tmp !Used to debug which frame has problem
 	!write(*,*) trim(c80tmp)
  !  backspace(fileid)
    
	read(fileid,"(15x)",advance="no")
    !write(*,*) ns,ne,n1 !Debug
	do j=ns,ne
		read(fileid,"(f11.4)",advance="no") ene(j)
	end do
	read(fileid,*)
	read(fileid,*) sym(ns:ne)
	do k=1,n1 !Scan rows in each frame
		read(fileid,"(15x)",advance='no') !Skip marker columns in each row
		do j=ns,ne !Scan elements in each row
 			!write(*,*) i,nt,k,j
			read(fileid,"(f11.6)",advance='no') mat(k,j)
		end do
		read(fileid,*)
	end do
end do
end subroutine









!=======================================================================
!=======================================================================
!!!!!!!!!!!!!!!! Below routines are used to output files !!!!!!!!!!!!!!!
!=======================================================================
!=======================================================================




!!---------- Interface of outputting pdb file
subroutine outpdb_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for outputting pdb file, e.g. C:\ltwd.pdb"
write(*,"(a)") " If press ENTER button directly, the system will be exported to "//trim(c200tmp)//".pdb in current folder"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".pdb"
call outpdb(outname,10)
end subroutine
!!---------- Output current coordinate to pdb file
subroutine outpdb(outpdbname,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outpdbname
character prefixname*6,tmpresname*3
integer ifileid,resid

open(ifileid,file=outpdbname,status="replace")
write(ifileid,"('REMARK   Generated by Multiwfn, Totally',i10,' atoms')") ncenter
if (ifPBC>0) then
    call getcellabc(asize,bsize,csize,alpha,beta,gamma)
    write(ifileid,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
end if
do i=1,ncenter
    call getprefixname(i,prefixname)
    tmpresname="MOL"
    if (a(i)%resname/=" ") tmpresname=a(i)%resname
	write(ifileid,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
	prefixname,i,' '//ind2name_up(a(i)%index)//' ',tmpresname, 'A',a(i)%resid,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,1D0,0D0,adjustr(ind2name_up(a(i)%index))
end do
write(ifileid,"('END')")
close(ifileid)
write(*,*) "Exporting pdb file finished!"
end subroutine


!!------- The same as subroutine outpdb, but also output image atoms at cell boundary if any. Residue index corresponds to index of unique atom
subroutine outpdb_PBC(outpdbname,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outpdbname
character prefixname*6,tmpresname*3
integer ifileid

call construct_atmp_withbound(ncenter_tmp)
open(ifileid,file=outpdbname,status="replace")
write(ifileid,"('REMARK   Generated by Multiwfn, Totally',i10,' atoms')") ncenter_tmp
call getcellabc(asize,bsize,csize,alpha,beta,gamma)
write(ifileid,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
do i=1,ncenter_tmp
    call getprefixname(a_tmp_idx(i),prefixname)
    tmpresname="MOL"
    if (a_tmp(i)%resname/=" ") tmpresname=a_tmp(i)%resname
	write(ifileid,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
	prefixname,i,' '//ind2name_up(a_tmp(i)%index)//' ',tmpresname, 'A',a_tmp_idx(i),a_tmp(i)%x*b2a,a_tmp(i)%y*b2a,a_tmp(i)%z*b2a,1D0,0D0,adjustr(ind2name_up(a_tmp(i)%index))
end do
write(ifileid,"('END')")
close(ifileid)
write(*,*) "Exporting pdb file finished!"
if (allocated(a_tmp)) deallocate(a_tmp,a_tmp_idx)
end subroutine


!!---------- Interface of outputting pqr file
subroutine outpqr_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for outputting pqr file, e.g. C:\ltwd.pqr"
write(*,"(a)") " If press ENTER button directly,the system will be exported to "//trim(c200tmp)//".pqr in current folder"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".pqr"
call outpqr(outname,10)
write(*,*) "Exporting pqr file finished!"
write(*,"(a)") " This file contains atomic charges that originally recorded in your .chg file. The radius column corresponds to Bondi vdW radii"
end subroutine
!!---------- Output current coordinate and atomic charges to pqr file
subroutine outpqr(outpqrname,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outpqrname
character resname*4,prefixname*6
integer i,ifileid
open(ifileid,file=outpqrname,status="replace")
write(ifileid,"('REMARK   Generated by Multiwfn, Totally',i10,' atoms')") ncenter
if (ifPBC>0) then
    call getcellabc(asize,bsize,csize,alpha,beta,gamma)
    write(ifileid,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
end if
do i=1,ncenter
    call getprefixname(i,prefixname)
	write(ifileid,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
	prefixname,i,' '//ind2name_up(a(i)%index)//' ',a(i)%resname, 'A',a(i)%resid,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,a(i)%charge,vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
end do
write(ifileid,"('END')")
close(ifileid)
end subroutine

!!-------- Determine "ATOM" or "HETATM" of an atom
subroutine getprefixname(iatm,prefixname)
use defvar
character prefixname*6,tmpname*3
integer iatm
prefixname="HETATM"
tmpname=adjustl(a(iatm)%resname)
if (tmpname=="ASP".or.tmpname=="GLU".or.tmpname=="HIS".or.tmpname=="LYS".or.tmpname=="ARG".or.&
    tmpname=="ASN".or.tmpname=="GLN".or.tmpname=="SER".or.tmpname=="THR".or.tmpname=="CYS".or.&
    tmpname=="VAL".or.tmpname=="ALA".or.tmpname=="LEU".or.tmpname=="PHE".or.tmpname=="TYR".or.&
    tmpname=="MET".or.tmpname=="ILE".or.tmpname=="PRO".or.tmpname=="GLY".or.tmpname=="TRP") then
    prefixname="ATOM  "
end if
end subroutine


!!---------- Interface of outputting xyz file
subroutine outxyz_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for outputting xyz file, e.g. C:\ltwd.xyz"
write(*,"(a)") " If press ENTER button directly, the system will be exported to "//trim(c200tmp)//".xyz in current folder"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".xyz"
call outxyz(outname,10)
write(*,*) "Exporting xyz file finished!"
end subroutine
!!---------- Output current coordinate to xyz file
subroutine outxyz(outxyzname,ifileid)
use defvar
character(len=*) outxyzname
integer i,ifileid
open(ifileid,file=outxyzname,status="replace")
write(ifileid,"(i6)") ncenter
if (ifPBC==0) then
    write(ifileid,*) "Generated by Multiwfn (http://sobereva.com/multiwfn)"
else !Write cell translation vectors, e.g. Tv_1:   9.901260   0.000000   0.000000 Tv_2:  -4.879808   8.533788   0.000000 Tv_3:   0.000000   0.000000  10.000000
    if (ifPBC==1) write(ifileid,"(' Tv_1:',3f11.6)") cellv1*b2a
    if (ifPBC==2) write(ifileid,"(' Tv_1:',3f11.6,' Tv_2:',3f11.6)") cellv1*b2a,cellv2*b2a
    if (ifPBC==3) write(ifileid,"(' Tv_1:',3f11.6,' Tv_2:',3f11.6,' Tv_3:',3f11.6)") cellv1*b2a,cellv2*b2a,cellv3*b2a
end if
do i=1,ncenter
	write(ifileid,"(a,3f16.8)") ind2name(a(i)%index),a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
close(ifileid)
end subroutine


!!---------- Output current coordinate to cml file
!iunit=0: In Angstrom, iunit=1: In Bohr
subroutine outcml(outcmlname,ifileid,iunit)
use defvar
character(len=*) outcmlname
character tmpstr1*80,tmpstr2*80,tmpstr3*80,tmpstr4*80
integer i,ifileid,iunit
if (.not.allocated(connmat)) call genconnmat(1,0) !Generate connectivity matrix

open(ifileid,file=outcmlname,status="replace")
write(ifileid,"(a)") '<molecule>'
write(ifileid,"(a)") ' <atomArray>'
fac=1
if (iunit==0) fac=b2a !Convert to Angstrom

do i=1,ncenter
    write(tmpstr1,*) i
    write(tmpstr2,"(f12.6)") a(i)%x*fac
    write(tmpstr3,"(f12.6)") a(i)%y*fac
    write(tmpstr4,"(f12.6)") a(i)%z*fac
    write(ifileid,"(a)") &
    '  <atom id="a'//trim(adjustl(tmpstr1))//'" elementType="'//trim(adjustl(a(i)%name))//'" x3="'//trim(adjustl(tmpstr2))&
    //'" y3="'//trim(adjustl(tmpstr3))//'" z3="'//trim(adjustl(tmpstr4))//'"/>'
end do
write(ifileid,"(a)") ' </atomArray>'
write(ifileid,"(a)") ' <bondArray>'
do i=1,ncenter
    do j=i+1,ncenter
        if (connmat(i,j)>0) then
            write(tmpstr1,*) i
            write(tmpstr2,*) j
            write(tmpstr3,*) connmat(i,j)
            write(ifileid,"(a)") '  <bond atomRefs2="a'//trim(adjustl(tmpstr1))//' a'//&
            trim(adjustl(tmpstr2))//'" order="'//trim(adjustl(tmpstr3))//'"/>'
        end if
    end do
end do
write(ifileid,"(a)") ' </bondArray>'
write(ifileid,"(a)") '</molecule>'

close(ifileid)
write(*,*) "Exporting cml file finished!"
end subroutine


!!---------- Output current coordinate to chg file
subroutine outchg(outchgname,ifileid)
use defvar
character(len=*) outchgname
integer i,ifileid
open(ifileid,file=outchgname,status="replace")
do i=1,ncenter
	write(ifileid,"(a,4f16.8)") ind2name(a(i)%index),a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,a(i)%charge
end do
close(ifileid)
write(*,*) "Exporting chg file finished!"
end subroutine



!!---------- Interface of outputting Gaussian input file
subroutine outgjf_wrapper
use util
use defvar
character(len=200) outname,c200tmp
character c2000tmp*2000
integer :: icoordtype=1 !1=Cartesian   2=Z-matrix with variable names   -2=Z-matrix directly with geometry parameters
integer,allocatable :: tmparr(:)

icustom=0
call path2filename(filename,c200tmp)
do while(.true.)
    write(*,*)
    write(*,*) "Input path for generating Gaussian input file, e.g. C:\ltwd.gjf"
    write(*,"(a)") " If press ENTER button directly, "//trim(c200tmp)//".gjf will be generated in current folder"
    if (ifPBC==0) then
        if (icoordtype==1) then
            write(*,"(a)") " Note: Cartesian coordinate will be outputted. If you want to output Z-matrix instead, input ""zmat"""
        else if (icoordtype==2) then
            write(*,"(a)") " Note: Z-matrix will be outputted with variable names. If you want to output Cartesian coordinates instead, input ""cart""; &
            if you want to output geometry parameters without variable names, input ""zmat2"""
        else if (icoordtype==-2) then
            write(*,"(a)") " Note: Z-matrix will be outputted without variable names. If you want to output Cartesian coordinates instead, input ""cart""; &
            if you want to output geometry parameters with variable names, input ""zmat1"""
        end if
    end if
    write(*,"(a)") " Hint 1: If template.gjf is presented in current folder, &
    then it will be used as template file, the line containing [geometry] or [GEOMETRY] will be replaced with the present coordinate, &
    [name] will be replaced with name of the new input file (without suffix), net charge and spin multiplicity will correspond to present system"
    write(*,"(a)") " Hint 2: If inputting ""custom"", you can specify charge, spin multiplicity and indices of the atoms for outputting"
    read(*,"(a)") outname
    if (outname=="zmat".or.outname=="zmat1") then
        icoordtype=2
    else if (outname=="zmat2") then
        icoordtype=-2
    else if (outname=="cart") then
        icoordtype=1
    else if (outname=="custom") then !Temporarily replace charge, spin multiplicit and atoms to the selected ones
		icustom=1
        loadcharge_tmp=loadcharge !Backup
        loadmulti_tmp=loadmulti
        ncenter_tmp=ncenter
		allocate(a_tmp(ncenter))
        a_tmp=a
        write(*,*) "Input net charge and spin multiplicity, e.g. 1 3"
        read(*,*) loadcharge,loadmulti
        write(*,*) "Input index of the atoms to output in the .gjf file, e.g. 2,5,7-10"
        read(*,"(a)") c2000tmp
		call str2arr(c2000tmp,ncenter)
        deallocate(a)
        allocate(a(ncenter),tmparr(ncenter))
		call str2arr(c2000tmp,ncenter,tmparr)
        do idx=1,ncenter
			a(idx)=a_tmp(tmparr(idx))
        end do
    else if (outname=="Bq") then !The same as custom, but inputted atoms will be treated as ghost atom. Mainly used in sobEDA shell script
		if (icoordtype/=1) then
			write(*,*) "Error: Only Cartesian coordinates can be used in this case"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            cycle
        end if
		icustom=2
        loadcharge_tmp=loadcharge !Backup
        loadmulti_tmp=loadmulti
        write(*,*) "Input net charge and spin multiplicity, e.g. 1 3"
        read(*,*) loadcharge,loadmulti
        write(*,*) "Input index of real atoms (all other atoms will be treated as ghost atoms)"
        write(*,*) "For example, 2,5,7-10,22"
        read(*,"(a)") c2000tmp
		call str2arr(c2000tmp,ntmp)
        allocate(tmparr(ntmp))
		call str2arr(c2000tmp,ntmp,tmparr)
        a%charge=-a%charge
        do idx=1,ntmp !Use negative charge to indicate ghost atoms
			iatm=tmparr(idx)
			a(iatm)%charge=abs(a(iatm)%charge)
        end do
    else if (outname==" ") then
        outname=trim(c200tmp)//".gjf"
        exit
    else
        exit
    end if
end do

call outgjf(outname,10,icoordtype)

if (icustom==1) then
    loadcharge=loadcharge_tmp
    loadmulti=loadmulti_tmp
    ncenter=ncenter_tmp
    deallocate(a)
    allocate(a(ncenter))
    a=a_tmp
	deallocate(a_tmp)
else if (icustom==2) then
	a(:)%charge=abs(a(:)%charge)
end if
end subroutine

!!---------- Output current coordinates and cell information to Gaussian input file
!1=Cartesian   2=Z-matrix with variable names   -2=Z-matrix directly with geometry parameters
subroutine outgjf(outgjfname,ifileid,icoordtype)
use defvar
use util
character(len=*) outgjfname
character selectyn,c10tmp*10,tmpname*200,c200tmp*200
integer icoordtype,Zmat(ncenter,3)
character(len=6) bondstr(ncenter),anglestr(ncenter),dihstr(ncenter)
real*8 bondval(ncenter),angleval(ncenter),dihval(ncenter)

selectyn='n'
if (allocated(CObasa).and.wfntype<=2) then
    write(*,"(a)") " Do you also want to write current wavefunction as initial guess into the .gjf file? (y/n)"
    read(*,*) selectyn
end if

if (abs(icoordtype)==2) then
    call genZmat(Zmat,ierror)
    if (ierror==1) then
        write(*,"(a)") " Error: Generation of internal coordinate was failed! May be this system cannot be properly represented by internal coordinate"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
end if

call path2filename(outgjfname,tmpname)
open(ifileid,file=outgjfname,status="replace")

!Determine net charge and spin multiplicity
if (loadcharge==-99) then !Not loaded from input file
    netcharge=nint(sum(a%charge)-nelec)
    if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
else
    netcharge=loadcharge
end if
if (loadmulti==-99) then !Not loaded from input file
    multi=nint(naelec-nbelec)+1
else
    multi=loadmulti
end if

!Write information before atomic coordinates
inquire(file="template.gjf",exist=alive)
if (alive) then !Write information in template.gjf to current gjf file except for coordinate part
    write(*,"(a)") " Note: template.gjf was found in current folder, it will be used as template file to generate new .gjf file"
    open(ifileid+1,file="template.gjf",status="old")
    nspace=0
    do while(.true.)
        read(ifileid+1,"(a)",iostat=ierror) c200tmp
        if (c200tmp==" ".and.nspace<2) nspace=nspace+1
        if (nspace==2) then !Write charge and spin multiplicity of present system
			if (loadmulti/=-99.or.allocated(b)) then !Charge and spin multiplicity are loaded or specified manually, or wavefunction information is available
				read(ifileid+1,*) !Skip charge and spin multiplicity line
				write(ifileid,"(/,2i3)") netcharge,multi
            else !Use charge and spin multiplicity in template.gjf
                read(ifileid+1,"(a)") c200tmp
				write(ifileid,"(/,a)") trim(c200tmp)
            end if
            nspace=nspace+1
            cycle
        end if
		if (index(c200tmp,"[geometry]")==0.and.index(c200tmp,"[GEOMETRY]")==0) then
			itmp=index(c200tmp,"[name]")
			if (itmp/=0) then
				call path2filename(outgjfname,tmpname)
				c200tmp=c200tmp(1:itmp-1)//trim(tmpname)//trim(c200tmp(itmp+6:))
			end if
			write(ifileid,"(a)") trim(c200tmp)
		else
			exit !Encountered [geometry] or [GEOMETRY]
		end if
    end do
else !Common case
    write(ifileid,"(a)") "%chk="//trim(tmpname)//".chk"
    if (selectyn=='y'.or.selectyn=='Y') then
        write(ifileid,"(a,/,/,a,/)") "#P B3LYP/6-31G* guess=cards","Generated by Multiwfn"
    else
        write(ifileid,"(a,/,/,a,/)") "#P B3LYP/6-31G*","Generated by Multiwfn"
    end if
    write(ifileid,"(2i3)") netcharge,multi
end if

!Write atomic coordinates
if (icoordtype==1) then !Cartesian
    do i=1,ncenter
	    if (a(i)%charge>=0) then
			write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
        else !Ghost atoms have negative charge temporarily
			write(ifileid,"(a,1x,3f14.8)") trim(a(i)%name)//"-Bq",a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
        end if
    end do
else if (icoordtype==2) then !Z-matrix with variable names
    nbond=0
    nangle=0
    ndih=0
    write(ifileid,"(a)") a(1)%name
    do iatm=2,ncenter
        write(ifileid,"(a)",advance="no") a(iatm)%name
        nbond=nbond+1
        write(c10tmp,"(i5)") nbond
        bondstr(nbond)="B"//trim(adjustl(c10tmp))
        i1=zmat(iatm,1)
        bondval(nbond)=atomdist(iatm,i1,0)*b2a
        write(ifileid,"(i5,3x,a)",advance="no") zmat(iatm,1),bondstr(nbond)
        if (iatm>=3) then
            nangle=nangle+1
            write(c10tmp,"(i5)") nangle
            anglestr(nangle)="A"//trim(adjustl(c10tmp))
            i2=zmat(iatm,2)
            angleval(nangle)=atomang(iatm,i1,i2,0)
            write(ifileid,"(3x,i5,3x,a)",advance="no") zmat(iatm,2),anglestr(nangle)
            if (iatm>=4) then
                ndih=ndih+1
                write(c10tmp,"(i5)") ndih
                dihstr(ndih)="D"//trim(adjustl(c10tmp))
                i3=zmat(iatm,3)
                dihval(ndih)=atomdih(iatm,i1,i2,i3,0)
                write(ifileid,"(3x,i5,3x,a)",advance="no") zmat(iatm,3),trim(dihstr(ndih))
            end if
        end if
        write(ifileid,*)
    end do
    write(ifileid,*)
    do ibond=1,nbond
        write(ifileid,"(a,f15.8)") bondstr(ibond),bondval(ibond)
    end do
    do iangle=1,nangle
        write(ifileid,"(a,f15.8)") anglestr(iangle),angleval(iangle)
    end do
    do idih=1,ndih
        write(ifileid,"(a,f15.8)") dihstr(idih),dihval(idih)
    end do
else if (icoordtype==-2) then !Z-matrix directly with geometry parameters
    write(ifileid,"(a)") a(1)%name
    do iatm=2,ncenter
        write(ifileid,"(a)",advance="no") a(iatm)%name
        i1=zmat(iatm,1)
        write(ifileid,"(i5,f12.6)",advance="no") zmat(iatm,1),atomdist(iatm,i1,0)*b2a
        if (iatm>=3) then
            i2=zmat(iatm,2)
            tmpval=atomang(iatm,i1,i2,0)
            write(ifileid,"(3x,i5,f10.4)",advance="no") zmat(iatm,2),tmpval
            if (iatm>=4) then
                i3=zmat(iatm,3)
                tmpval=atomdih(iatm,i1,i2,i3,0)
                write(ifileid,"(3x,i5,f10.4)",advance="no") zmat(iatm,3),tmpval
            end if
        end if
        write(ifileid,*)
    end do
end if

!Write translation vectors
if (ifPBC>=1) then
	write(ifileid,"('Tv',3f12.6)") cellv1*b2a
	if (ifPBC>=2) then
		write(ifileid,"('Tv',3f12.6)") cellv2*b2a
		if (ifPBC==3) write(ifileid,"('Tv',3f12.6)") cellv3*b2a
    end if
end if

!Write initial guess
if (selectyn=='y') then
    write(ifileid,"(/,'(5E16.9)',/,'-1')")
    do i=1,nbasis
	    if (wfntype==0.or.wfntype==2) write(ifileid,"('! Orbital:',i6,' Occ:',f10.6)") i,MOocc(i)
	    if (wfntype==1) write(ifileid,"('! Alpha orbital:',i6,' Occ:',f10.6)") i,MOocc(i)
	    write(ifileid,"(5E16.9)") (CObasa(j,i),j=1,nbasis)
    end do
    if (wfntype==1) then
	    write(ifileid,"('-1')")
	    do i=1,nbasis
		    write(ifileid,"('! Beta orbital:',i6,' Occ:',f10.6)") i,MOocc(nbasis+i)
		    write(ifileid,"(5E16.9)") (CObasb(j,i),j=1,nbasis)
	    end do
    end if
    write(ifileid,"('0',/)")
end if

!Write rest part of template.gjf into present gjf file
if (alive) then
    do while(.true.)
        read(ifileid+1,"(a)",iostat=ierror) c200tmp
        if (ierror/=0) then
            exit
        else
            itmp=index(c200tmp,"[name]")
            if (itmp/=0) then
                call path2filename(outgjfname,tmpname)
                c200tmp=c200tmp(1:itmp-1)//trim(tmpname)//trim(c200tmp(itmp+6:))
            end if
            write(ifileid,"(a)") trim(c200tmp)
        end if
    end do
    close(ifileid+1)
else
    write(ifileid,*)  !Two blank lines at the end of the file
end if

close(ifileid)
write(*,"(a)") " Exporting Gaussian input file finished! It corresponds to single point task at B3LYP/6-31G* level"
if (selectyn=='y') write(*,"(a)") " Note that you must specify the basis set to the one &
&originally used to yield the wavefunction, otherwise Gaussian calculation must be failed"
end subroutine


!!---------- Output current coordinate to GAMESS-US input file
!The level exactly corresponds to the B3LYP-D3(BJ)/6-31G* of Gaussian
subroutine outGAMESSinp(outname,ifileid)
use defvar
character(len=*) outname
character SCFTYPE*3,selectyn
ioutguess=0
if (allocated(CObasa)) then
	write(*,*) "If write initial guess information? (y/n)"
	read(*,*) selectyn
	if (selectyn=='y'.or.selectyn=='Y') ioutguess=1
end if
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
mult=nint(naelec-nbelec)+1
iopsh=0
if (mult/=1.or.wfntype==1.or.wfntype==4) iopsh=1
SCFTYPE="RHF"
if (iopsh==1) SCFTYPE="UHF"
open(ifileid,file=outname,status="replace")
write(ifileid,"(a,i1,a,i1,a)") " $CONTRL SCFTYP="//SCFTYPE//" MULT=",mult," ICHARG=",netcharge," RUNTYP=ENERGY"
write(ifileid,"(a)") " DFTTYP=B3LYPV1R ISPHER=0 MAXIT=60 NPRINT=-5 $END"
write(ifileid,"(a)") " $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=0 DIFFSP=.F. DIFFS=.F. $END"
write(ifileid,"(a)") " $SYSTEM MWORDS=200 $END"
!  $lmoeda matom(1)=3,4 mcharg(1)=0,0 mmult(1)=1,1 $end
write(ifileid,"(a)") " $SCF DIRSCF=.T. $END"
write(ifileid,"(a)") " $DFT IDCVER=4 $END"
if (ioutguess==1) write(ifileid,"(a,i5,a)") " $GUESS GUESS=MOREAD NORB=",nbasis," $END"
write(ifileid,"(a)") " $DATA"
write(ifileid,"(a)") "Generated by Multiwfn"
write(ifileid,"(a)") "C1"
do i=1,ncenter
	write(ifileid,"(a,1x,f5.1,1x,3f14.8)") a(i)%name,dfloat(a(i)%index),a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,"(a)") " $END"
! Write $VEC information
if (ioutguess==1) then
	write(ifileid,*)
	write(ifileid,"(a)") " $VEC"
	ntime=ceiling(nbasis/5D0)
	idxmo=1
	do imo=1,nbasis
		do itime=1,ntime
			if (itime<ntime) then
				write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasa((itime-1)*5+1:itime*5,imo)
			else
				write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasa((itime-1)*5+1:nbasis,imo)
			end if
		end do
		idxmo=idxmo+1
		if (idxmo==100) idxmo=0
	end do
	if (iopsh==1) then
		idxmo=1
		do imo=1,nbasis
			do itime=1,ntime
				if (itime<ntime) then
					write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasb((itime-1)*5+1:itime*5,imo)
				else
					write(ifileid,"(i2,i3,5(1PE15.8))") idxmo,itime,CObasb((itime-1)*5+1:nbasis,imo)
				end if
			end do
			idxmo=idxmo+1
			if (idxmo==100) idxmo=0
		end do
	end if
	write(ifileid,"(a)") " $END"
end if
close(ifileid)
write(*,"(a)") " Exporting GAMESS-US input file finished! It corresponds to single point task at B3LYP-D3(BJ)/6-31G* level"
end subroutine


!!---------- Interface of outputting ORCA input file
subroutine outORCAinp_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,"(/,a)") " Note: If you benefit from this function when creating ORCA input file, please cite Multiwfn original paper, thanks!"
write(*,*)
write(*,*) "Input path for generating ORCA input file, e.g. C:\ltwd.inp"
write(*,"(a)") " If press ENTER button directly, will be exported to "//trim(c200tmp)//".inp"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".inp"
call outORCAinp(outname,10)
end subroutine
!!---------- Output current coordinate to ORCA input file
subroutine outORCAinp(outname,ifileid)
use defvar
use util
character(len=*) outname
character keyword*200,c80tmp*80,c200tmp*200,solvname*30,c2000tmp*2000,grd4str*20,grd5str*20
character(len=50) atmbasname(ncenter)
integer :: itask=1,idiffuse=0,frozenatm(ncenter),CPfrag1(ncenter),CPfrag2(ncenter),ioptHonly=0,idxarr(ncenter),ifixHbond=0
integer :: maxcore=1000
integer :: iORCAver=1 !ORCA version, =0: 4.x, =1: 5.0
real*8 efield(3)
nprocs=nthreads
isolv=0
nfrozenatm=0
efield=0
atmbasname=" "

do while(.true.)
    write(*,*)
    write(*,*) "-100 Use template input file provided by user to generate new input file"
    if (iORCAver==0) write(*,*) "-11 Choose ORCA version compatibility, current: ORCA 4.x"
    if (iORCAver==1) write(*,*) "-11 Choose ORCA version compatibility, current: >ORCA 5.0"
    write(*,"(a,i4,a,i6,a)") " -10 Set computational resources, core:",nprocs," memory/core:",maxcore," MB"
    write(*,"(a,i5,a)") " -5 Set mixed basis set, current: Defined for",count(atmbasname/=" ")," atoms"
    write(*,*) "-4 Other settings"
    if (itask==2.or.itask==4.or.itask==5) then
        if (ioptHonly==0) then
            write(*,"(a,', current:',i6' atoms')") " -3 Set atom freeze in optimization",nfrozenatm
        else if (ioptHonly==1) then
            write(*,"(a)") " -3 Set atom freeze in optimization, current: Only optimize hydrogens"
        end if
    end if
    if (idiffuse==1) write(*,*) "-2 Toggle adding diffuse functions, current: Yes"
    if (idiffuse==0) write(*,*) "-2 Toggle adding diffuse functions, current: No"
    if (isolv==0) write(*,*) "-1 Toggle employing implicit solvation model, current: No"
    if (isolv==1) write(*,*) "-1 Toggle employing implicit solvation model, current: SMD, "//trim(solvname)
    if (isolv==2) write(*,*) "-1 Toggle employing implicit solvation model, current: CPCM, "//trim(solvname)
    if (itask==1) write(*,*) "0 Select task, current: Single point"
    if (itask==2) write(*,*) "0 Select task, current: Optimizing minimum"
    if (itask==3) write(*,*) "0 Select task, current: Frequency analysis"
    if (itask==4) write(*,*) "0 Select task, current: Optimizing minimum + Frequency"
    if (itask==5) write(*,*) "0 Select task, current: Optimizing TS + Frequency"
    if (itask==6) write(*,*) "0 Select task, current: Molecular dynamics"
    if (itask==7) write(*,*) "0 Select task, current: Interaction energy with counterpoise correction"
    if (itask==-7) write(*,*) "0 Select task, current: Interaction and complex energies with CP correction"
    if (itask==8) write(*,*) "0 Select task, current: NMR"
    if (itask==9) write(*,*) "0 Select task, current: Wavefunction stability test"
    !With b suffix, the internal index is original index plus 1000
    if (iORCAver==0) write(*,*) "1 B97-3c"
    if (iORCAver==1) write(*,*) "1 B97-3c      1b r2SCAN-3c"
    write(*,*) "2 RI-BLYP-D3(BJ)/def2-TZVP"
    write(*,*) "3 RI-B3LYP-D3(BJ)/def2-TZVP(-f)     4 RI-B3LYP-D3(BJ)/def2-TZVP"
    write(*,*) "5 RI-wB97M-V/def2-TZVP"
    write(*,*) "6 RI-PWPB95-D4/def2-TZVPP       7 RI-PWPB95-D4/def2-QZVPP"
    if (iORCAver==1) then
        write(*,*) "6b RI-wB97X-2-D3(BJ)/def2-TZVPP     7b RI-wB97X-2-D3(BJ)/def2-QZVPP"
        write(*,*) "6c RI-revDSD-PBEP86-D4/def2-TZVPP   7c RI-revDSD-PBEP86-D4/def2-QZVPP"
    end if
    write(*,*) "8 DLPNO-CCSD(T)/def2-TZVPP with normalPNO and RIJCOSX"
    write(*,*) "9 DLPNO-CCSD(T)/def2-TZVPP with tightPNO and RIJCOSX"
    write(*,*) "10 CCSD(T)/cc-pVTZ"
    write(*,*) "11 CCSD(T)-F12/cc-pVDZ-F12 with RI"
    write(*,*) "12 Approximated CCSD(T)/CBS with help of MP2 (cc-pVTZ->QZ extrapolation)"
    if (idiffuse==0) write(*,*) "13 DLPNO-CCSD(T)/CBS with RIJCOSX & tightPNO (def2-TZVPP->QZVPP extrapolation)"
    write(*,*) "14 CCSD(T)/CBS (cc-pVTZ->QZ extrapolation)"
    write(*,*) "20 sTD-DFT based on RI-wB97X-D3/def2-SV(P) orbitals"
    !write(*,*) "21 TDA-DFT RI-PBE0/def2-SV(P) with riints_disk (much faster than 22)" !riints_disk no longer supported by ORCA 6
    write(*,*) "22 TDDFT RI-PBE0/def2-SV(P)"
    write(*,*) "23 TDDFT RI-RSX-QIDH/def2-TZVP     231 TDDFT RI-DSD-PBEP86/def2-TZVP"
    write(*,*) "24 EOM-CCSD/cc-pVTZ                 25 STEOM-DLPNO-CCSD/def2-TZVP(-f)"
    read(*,"(a)") c80tmp
    if (c80tmp=="1b") then
        ilevel=1001
    else if (c80tmp=="6b") then
        ilevel=1006
    else if (c80tmp=="7b") then
        ilevel=1007
    else if (c80tmp=="6c") then
        ilevel=1008
    else if (c80tmp=="7c") then
        ilevel=1009
    else
        read(c80tmp,*) ilevel
    end if
    if (ilevel==-100) then
        do while(.true.)
            write(*,*) "Input the path of template input file, e.g. /sob/Bang/Dream.inp"
            write(*,"(a)") " Note: The template input file should contain all content of finally generated input file, &
            but coordinate part should be written as [geometry] or [GEOMETRY], which will be replaced with geometry of present system"
	        read(*,"(a)") c200tmp
	        inquire(file=c200tmp,exist=alive)
	        if (alive) exit
	        write(*,*) "Cannot find the file, input again!"
        end do
        open(ifileid,file=outname,status="replace")
        open(ifileid+1,file=c200tmp,status="old")
        do while(.true.)
            read(ifileid+1,"(a)",iostat=ierror) c200tmp
            if (ierror/=0) exit
            if (index(c200tmp,"[geometry]")/=0.or.index(c200tmp,"[GEOMETRY]")/=0) then
                do iatm=1,ncenter
                    if (atmbasname(iatm)==" ") then
                        write(ifileid,"(a,1x,3f14.8)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
                    else
                        write(ifileid,"(a,1x,3f14.8,a)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a," newGTO """//trim(atmbasname(iatm))//""" end"
                    end if
                end do
            else
                write(ifileid,"(a)") trim(c200tmp)
            end if
        end do
        close(ifileid)
        close(ifileid+1)
        write(*,"(a)") " Done! "//trim(outname)//" has been generated"
        return
    else if (ilevel==-11) then
        write(*,*) "Choose ORCA version compatibility"
        write(*,*) "0 I am using ORCA 4.x"
        write(*,*) "1 I am using ORCA 5.0 or newer"
        read(*,*) iORCAver
    else if (ilevel==0) then
        write(*,*) "1 Single point"
        write(*,*) "2 Optimization"
        write(*,*) "3 Frequency"
        write(*,*) "4 Optimization + Frequency"
        write(*,*) "5 Optimization for transition state + Frequency"
        write(*,*) "6 Molecular dynamics (MD)   6a MD with fixing all H related bond lengths"
        write(*,*) "7 Interaction energy with counterpoise correction"
        write(*,*) "-7 Interaction energy and complex energy with counterpoise correction"
        write(*,*) "8 NMR"
        write(*,*) "9 Wavefunction stability test"
        read(*,*) c80tmp
        ifixHbond=0
        if (c80tmp=="6a") then
            ifixHbond=1
            c80tmp(2:2)=" "
            write(*,"(a)") " Note: Lengths of all hydrogen related bonds will be constraint to initial values. &
            The existence of the bonds are in line with that observed in current status of main function 0"
        end if
        read(c80tmp,*) itask
        if (itask==7.or.itask==-7) then
            write(*,*) "Input indices of the atoms in fragment 1, e.g. 3,6-10,14"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,nCPfrag1,CPfrag1)
            write(*,*) "Input indices of the atoms in fragment 2, e.g. 1,2,4,5,11-13"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,nCPfrag2,CPfrag2)
        end if
    else if (ilevel==-1) then
        if (isolv/=0) then
            isolv=0
        else
            write(*,"(a)") " Input name of solvent for SMD solvation model, e.g. ethanol. The full list of available solvents can be found in ""The SMD Solvation Model"" of ORCA manual"
            write(*,*) "If press ENTER button directly, water will be employed as solvent"
            write(*,*) "If you want to use CPCM solvation model or use custom solvent, input ""c"""
            read(*,"(a)") solvname
            isolv=1
            if (solvname==" ") then
                solvname="water"
            else if (solvname=='c') then
                isolv=2
                write(*,"(a)") " Input name of one of supported solvents: Water, Acetone, Acetonitrile, Ammonia, &
                Benzene, CCl4, CH2Cl2, Chloroform, Cyclohexane, DMF, DMSO, Ethanol, Hexane, Methanol, Octanol, Pyridine, THF, Toluene"
                write(*,*) "If you want to defined custom solvent, input ""c"""
                read(*,*) solvname
                if (solvname=="c") then
                    write(*,*) "Input dielectic constant, e.g. 80"
                    read(*,*) diecons
                    write(*,*) "Input refractive index, e.g. 1.8"
                    write(*,*) "If actual value is known, usually 1.33 is an acceptable approximation"
                    read(*,*) refrac
                    solvname="custom"
                end if
            end if
        end if
    else if (ilevel==-2) then
        if (idiffuse==1) then
            idiffuse=0
        else
            idiffuse=1
        end if
    else if (ilevel==-3) then
        write(*,*) "Input the atoms to be frozen during optimization, e.g. 3,5-10,13,15"
        write(*,"(a)") " Note: The index starts from 1 (Multiwfn convention) rather than 0 (ORCA convention)"
        write(*,*) "If you only want to optimize hydrogens, input ""H"""
        read(*,"(a)") c2000tmp
        if (index(c2000tmp,'H')/=0) then
            ioptHonly=1
        else
            ioptHonly=0
            call str2arr(c2000tmp,nfrozenatm,frozenatm)
        end if
    else if (ilevel==-4) then
        do while(.true.)
            write(*,*)
            write(*,*) "0 Return"
            write(*,"(a,3f10.6,' a.u.')") " 1 Set external electric field, current:",efield(:)
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
                write(*,*) "Write X/Y/Z of external electric field vector in a.u., e.g. 0.002,0,0.0015"
                read(*,*) efield
            end if
        end do
    else if (ilevel==-5) then
        do while(.true.)
            write(*,*)
            if (all(atmbasname==" ")) then
                write(*,*) "Mixed basis set has not been set"
            else
                write(*,*) "Atoms with specific basis set (i.e. not identical to the calculation level):"
                do iatm=1,ncenter
                    if (atmbasname(iatm)/=" ") write(*,"(i6,'(',a,'):',1x,a)") iatm,ind2name(a(iatm)%index),trim(atmbasname(iatm))
                end do
            end if
            write(*,*)
            write(*,*) "Input indices of atoms for which specific basis set will be used"
            write(*,*) "For example, 2-8,13,16-31,89"
            write(*,*) "You can also input an element, e.g. Fe"
            write(*,*) "To exit, input ""q"""
            read(*,"(a)") c2000tmp
            if (index(c2000tmp,'q')/=0) then
                exit
            else
                if (iachar(c2000tmp(1:1))>=48.and.iachar(c2000tmp(1:1))<=57) then
                    call str2arr(c2000tmp,nsel,idxarr)
                else
                    call elename2idx(c2000tmp,iele)
                    if (iele/=0) then
                        nsel=0
                        do iatm=1,ncenter
                            if (a(iatm)%index==iele) then
                                nsel=nsel+1
                                idxarr(nsel)=iatm
                            end if
                        end do
                    else
                        write(*,*) "Error: Unable to recognize your input"
                        cycle
                    end if
                end if
            end if
            write(*,*) "Input basis set name, e.g. def2-TZVP"
            read(*,"(a)") c80tmp
            do idx=1,nsel
                atmbasname(idxarr(idx))=trim(c80tmp)
            end do
        end do
        
    else if (ilevel==-10) then
        write(*,*) "Input number of cores, e.g. 12"
        write(*,"(' If press ENTER button directly, current value',i4,' will be unchanged')") nprocs
        read(*,"(a)") c80tmp
        if (c80tmp/=" ") read(c80tmp,*) nprocs
        write(*,*) "Input usable memory per core, e.g. 1000 (in MB)"
        write(*,"(' If press ENTER button directly, current value',i6,' will be unchanged')") maxcore
        read(*,"(a)") c80tmp
        if (c80tmp/=" ") read(c80tmp,*) maxcore
    else
        exit
    end if
end do

grd4str=" grid4 gridx4"
grd5str=" grid5 gridx5"
if (iORCAver==1) then !Since ORCA 5.0, the default grid is defgrid2, and grid/gridx keyword is inacceptable
    grd4str=" "
    grd5str=" "
end if
if (idiffuse==0) then
    if (ilevel==1) c200tmp="! B97-3c"
    if (ilevel==1001) c200tmp="! r2SCAN-3c"
    if (ilevel==2) c200tmp="! BLYP D3 def2-TZVP def2/J" !For pure functional, RIJ is used by default even without def2/J
    if (ilevel==3) c200tmp="! B3LYP D3 def2-TZVP(-f) def2/J RIJCOSX"
    if (ilevel==4) c200tmp="! B3LYP D3 def2-TZVP def2/J RIJCOSX"
    if (ilevel==5) c200tmp="! wB97M-V def2-TZVP def2/J RIJCOSX strongSCF"//trim(grd5str)
    if (ilevel==6) c200tmp="! PWPB95 D4 def2-TZVPP def2/J def2-TZVPP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==7) c200tmp="! PWPB95 D4 def2-QZVPP def2/J def2-QZVPP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1006) c200tmp="! wB97X-2 D3 def2-TZVPP def2/J def2-TZVPP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1007) c200tmp="! wB97X-2 D3 def2-QZVPP def2/J def2-QZVPP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1008) c200tmp="! revDSD-PBEP86-D4/2021 def2-TZVPP def2/J def2-TZVPP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1009) c200tmp="! revDSD-PBEP86-D4/2021 def2-QZVPP def2/J def2-QZVPP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==8) c200tmp="! DLPNO-CCSD(T) normalPNO RIJCOSX def2-TZVPP def2/J def2-TZVPP/C tightSCF"
    if (ilevel==9) c200tmp="! DLPNO-CCSD(T) tightPNO RIJCOSX def2-TZVPP def2/J def2-TZVPP/C tightSCF"
    if (ilevel==10) c200tmp="! CCSD(T) cc-pVTZ tightSCF"
    if (ilevel==11) c200tmp="! CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C tightSCF"
    if (ilevel==12) c200tmp="! ExtrapolateEP2(3/4,cc,MP2) tightSCF"
    if (ilevel==13) c200tmp="! DLPNO-CCSD(T) tightPNO Extrapolate(3/4,def2) def2-QZVPP/C RIJCOSX tightSCF"
    if (ilevel==14) c200tmp="! CCSD(T) Extrapolate(3/4,cc) tightSCF"
    if (ilevel==20) c200tmp="! wB97X-D3 def2-SV(P) def2/J RIJCOSX"
    if (ilevel==21) c200tmp="! PBE0 def2-SV(P) def2/J def2-SVP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==22) c200tmp="! PBE0 def2-SV(P) def2/J RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==23) c200tmp="! RSX-QIDH def2-TZVP def2/J def2-TZVP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==231) c200tmp="! DSD-PBEP86 def2-TZVP def2/J def2-TZVP/C RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==24) c200tmp="! EOM-CCSD cc-pVTZ tightSCF"
    if (ilevel==25) c200tmp="! STEOM-DLPNO-CCSD RIJCOSX def2-TZVP(-f) def2/J def2-TZVP/C tightSCF"
else
    if (ilevel==1) c200tmp="! B97-3c"
    if (ilevel==1001) c200tmp="! r2SCAN-3c"
    if (ilevel==2) c200tmp="! BLYP D3 ma-def2-TZVP autoaux"
    if (ilevel==3) c200tmp="! B3LYP D3 ma-def2-TZVP(-f) autoaux RIJCOSX"
    if (ilevel==4) c200tmp="! B3LYP D3 ma-def2-TZVP autoaux RIJCOSX"
    if (ilevel==5) c200tmp="! wB97M-V ma-def2-TZVP autoaux RIJCOSX strongSCF"//trim(grd5str)
    if (ilevel==6) c200tmp="! PWPB95 D4 ma-def2-TZVPP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==7) c200tmp="! PWPB95 D4 ma-def2-QZVPP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1006) c200tmp="! wB97X-2 D3 ma-def2-TZVPP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1007) c200tmp="! wB97X-2 D3 ma-def2-QZVPP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1008) c200tmp="! revDSD-PBEP86-D4/2021 ma-def2-TZVPP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==1009) c200tmp="! revDSD-PBEP86-D4/2021 ma-def2-QZVPP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==8) c200tmp="! DLPNO-CCSD(T) normalPNO RIJCOSX ma-def2-TZVPP autoaux tightSCF"
    if (ilevel==9) c200tmp="! DLPNO-CCSD(T) tightPNO RIJCOSX ma-def2-TZVPP autoaux tightSCF"
    if (ilevel==10) c200tmp="! CCSD(T) aug-cc-pVTZ tightSCF"
    if (ilevel==11) c200tmp="! CCSD(T)-F12/RI cc-pVDZ-F12 cc-pVDZ-F12-CABS cc-pVTZ/C"
    if (ilevel==12) c200tmp="! ExtrapolateEP2(3/4,aug-cc,MP2) tightSCF"
    !if (ilevel==13) c200tmp="! DLPNO-CCSD(T) tightPNO Extrapolate(3/4,def2) tightSCF !This is unavailable when diffuse function is requested
    if (ilevel==14) c200tmp="! CCSD(T) Extrapolate(3/4,aug-cc) tightSCF"
    if (ilevel==20) c200tmp="! wB97X-D3 ma-def2-SV(P) autoaux RIJCOSX"
    if (ilevel==21) c200tmp="! PBE0 ma-def2-SV(P) autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==22) c200tmp="! PBE0 ma-def2-SV(P) autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==23) c200tmp="! RSX-QIDH ma-def2-TZVP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==231) c200tmp="! DSD-PBEP86 ma-def2-TZVP autoaux RIJCOSX tightSCF"//trim(grd4str)
    if (ilevel==24) c200tmp="! EOM-CCSD aug-cc-pVTZ tightSCF"
    if (ilevel==25) c200tmp="! STEOM-DLPNO-CCSD RIJCOSX ma-def2-TZVP(-f) autoaux tightSCF"
    
    if (ilevel==1.or.ilevel==11) then
        write(*,*) "Error: Diffuse functions cannot be added to this level!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
end if

if (abs(itask)==7.and.isolv/=0) then
    write(*,*) "Error: It is meaningless to perform counterpoise calculation with solvation model!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

!Generate keywords
if ((itask==2.or.itask==3.or.itask==4.or.itask==5.or.itask==8).and.index(c200tmp,"strongSCF")/=0) then !Remove strongSCF when tightSCF is needed
	itmp=index(c200tmp,"strongSCF")
    c200tmp(itmp:)=" "
end if
if (itask==2) then !opt task uses tightSCF by default
    keyword=trim(c200tmp)//" opt"
    if (ilevel==22.and.idiffuse==0) keyword=trim(keyword)//" def2-SVP/C" !RI-TDDFT opt must use /C
else if (itask==3) then !freq
    if (index(c200tmp,"tightSCF")==0) then
        keyword=trim(c200tmp)//" freq tightSCF"
    else
        keyword=trim(c200tmp)//" freq"
    end if
else if (itask==4) then !opt freq
    if (index(c200tmp,"tightSCF")==0) then
        keyword=trim(c200tmp)//" opt freq tightSCF"
    else
        keyword=trim(c200tmp)//" opt freq"
    end if
    if (ilevel==22.and.idiffuse==0) keyword=trim(keyword)//" def2-SVP/C" !RI-TDDFT opt must use /C
else if (itask==5) then !optTS
    if (index(c200tmp,"tightSCF")==0) then
        keyword=trim(c200tmp)//" optTS freq tightSCF"
    else
        keyword=trim(c200tmp)//" optTS freq"
    end if
    if (ilevel==22.and.idiffuse==0) keyword=trim(keyword)//" def2-SVP/C" !RI-TDDFT opt must use /C
else if (itask==8) then !NMR
    if (index(c200tmp,"tightSCF")==0) then
        keyword=trim(c200tmp)//" NMR tightSCF"
    else
        keyword=trim(c200tmp)//" NMR"
    end if
else
    keyword=trim(c200tmp)
end if
keyword=trim(keyword)//" noautostart miniprint"

if (loadcharge==-99) then !Not loaded from input file
    netcharge=nint(sum(a%charge)-nelec)
    if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
else
    netcharge=loadcharge
end if
if (loadmulti==-99) then !Not loaded from input file
    multi=nint(naelec-nbelec)+1
else
    multi=loadmulti
end if

open(ifileid,file=outname,status="replace")

!Counterpoise task is a special case. Directly write input file and return
if (itask==7.or.itask==-7) then
    write(ifileid,"('%pal nprocs',i4,' end')") nprocs
    write(ifileid,"(/,a)") trim(keyword)
    if (ilevel==1006.or.ilevel==1007) call ORCA_DFT_D3_parm(ifileid)
    !if (ilevel==1008.or.ilevel==1009) call ORCA_revDSD_PBEP86_D4_parm(ifileid)
    write(ifileid,"(a,i9)") "%maxcore",maxcore
    write(ifileid,"('* xyz',2i4)") netcharge,multi
    do i=1,ncenter
        if (atmbasname(i)==" ") then
	        write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
        else
            write(ifileid,"(a,1x,3f14.8,a)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a," newGTO """//trim(atmbasname(i))//""" end"
        end if
    end do
    write(ifileid,"(a)") "*"
    do itime=1,2
        write(ifileid,"(/,a)") "$new_job"
        write(ifileid,"(a)") trim(keyword)//" Pmodel" !Use Pmodel to regenerate new initial guess rather than use the previous wavefunction
        if (ilevel==1006.or.ilevel==1007) call ORCA_DFT_D3_parm(ifileid)
        !if (ilevel==1008.or.ilevel==1009) call ORCA_revDSD_PBEP86_D4_parm(ifileid)
        write(ifileid,"(a,i9)") "%maxcore",maxcore
        write(ifileid,"('* xyz',2i4)") netcharge,multi
        do i=1,ncenter
	        if ((itime==1.and.any(CPfrag2==i)).or.(itime==2.and.any(CPfrag1==i))) then
                if (atmbasname(i)==" ") then
                    write(ifileid,"(a,':',1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
                else    
                    write(ifileid,"(a,':',1x,3f14.8,a)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a," newGTO """//trim(atmbasname(i))//""" end"
                end if
            else
                if (atmbasname(i)==" ") then
                    write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
                else
                    write(ifileid,"(a,1x,3f14.8,a)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a," newGTO """//trim(atmbasname(i))//""" end"
                end if
            end if
        end do
        write(ifileid,"(a)") "*"
    end do
    if (itask==-7) then
        do itime=1,2
            write(ifileid,"(/,a)") "$new_job"
            write(ifileid,"(a)") trim(keyword)//" Pmodel" !Use Pmodel to regenerate new initial guess rather than use the previous wavefunction
            if (ilevel==1006.or.ilevel==1007) call ORCA_DFT_D3_parm(ifileid)
            !if (ilevel==1008.or.ilevel==1009) call ORCA_revDSD_PBEP86_D4_parm(ifileid)
            write(ifileid,"(a,i9)") "%maxcore",maxcore
            write(ifileid,"('* xyz',2i4)") netcharge,multi
            do i=1,ncenter
	            if ((itime==1.and.any(CPfrag1==i)) .or. (itime==2.and.any(CPfrag2==i))) then
                    write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
                end if
            end do
            write(ifileid,"(a)") "*"
        end do
    end if
    close(ifileid)
    write(*,"(a)") " Exporting ORCA input file finished!"
    if (itask==7) then
        write(*,"(a)") " If you run it, single point task will be carried out three times, and correspondingly, the ""FINAL SINGLE POINT ENERGY"" &
        data will be printed three times (EAB, EA and EB in turn), the interaction energy after counterpoise correction should be &
        manually calculated as EAB-EA-EB. You can also use the ""getEbind.sh"" script mentioned in http://sobereva.com/542 to automatically calculate this for all .out files in current folder."
    else
        write(*,"(a)") " If you run it, single point task will be carried out five times, and correspondingly, the ""FINAL SINGLE POINT ENERGY"" &
        data will be printed five times in this order: EAB, EA', EB', EA, EB. The interaction energy after counterpoise correction should be &
        manually calculated as EAB-EA'-EB'. BSSE correction energy is E_BSSE=(EA-EA')+(EB-EB'). BSSE corrected complex energy is EAB+E_BSSE. &
        You can also use the ""ORCA_CP.sh"" script mentioned in http://sobereva.com/542 to automatically calculate these terms for all .out files in current folder."
    end if
    write(*,"(a)") " In addition, please check charge and spin multiplicity to ensure their appropriateness"
    return
end if

!Start to write input file for normal cases
write(ifileid,"(a)") trim(keyword)
write(ifileid,"(a,i9)") "%maxcore",maxcore
write(ifileid,"('%pal nprocs',i4,' end')") nprocs
!Write solvation model information
if (isolv==2.and.solvname/="custom") then
    write(ifileid,"(a)")"! CPCM("//trim(solvname)//")"
else if (isolv>0) then
    write(ifileid,"(a)") "%cpcm"
    if (isolv==1) then
        write(ifileid,"(a)") "smd true"
        write(ifileid,"(a)") 'SMDsolvent "'//trim(solvname)//'"'
    else if (isolv==2) then
        write(ifileid,"(a,f8.2)") "epsilon",diecons
        write(ifileid,"(a,f8.2)") "refrac",refrac
    end if
    write(ifileid,"(a)") "end"
end if

!Write additional D3 parameter or functional definition
if (ilevel==1006.or.ilevel==1007) call ORCA_DFT_D3_parm(ifileid)
!if (ilevel==1008.or.ilevel==1009) call ORCA_revDSD_PBEP86_D4_parm(ifileid)

!Write task specific information
if (itask==5) then !Calculate Hessian for optTS
    write(ifileid,"(a)") "%geom Calc_Hess true end"
else if (itask==6) then !MD
    write(ifileid,"(a)") "%md"
    write(ifileid,"(a)") "#restart ifexists  # Continue MD by reading [basename].mdrestart if it exists. In this case ""initvel"" should be commented"
    write(ifileid,"(a)") "#minimize  # Do minimization prior to MD simulation"
    write(ifileid,"(a)") " timestep 0.5_fs  # This stepsize is safe at several hundreds of Kelvin"
    write(ifileid,"(a)") " initvel 100_K no_overwrite # Assign velocity according to temperature for atoms whose velocities are not available"
    if (iORCAver==0) then
        write(ifileid,"(a)") " thermostat berendsen 298.15_K timecon 30.0_fs  # Target temperature and coupling time constant"
    else if (iORCAver==1) then
        write(ifileid,"(a)") " thermostat CSVR 298.15_K timecon 30.0_fs  # Target temperature and coupling time constant"
    end if
    write(ifileid,"(a)") " dump position stride 1 format xyz filename ""pos.xyz""  # Dump position every ""stride"" steps"
    write(ifileid,"(a)") "#dump force stride 1 format xyz filename ""force.xyz""  # Dump force every ""stride"" steps"
    write(ifileid,"(a)") "#dump velocity stride 1 format xyz filename ""vel.xyz""  # Dump velocity every ""stride"" steps"
    write(ifileid,"(a)") "#dump gbw stride 20 filename ""wfn""  # Dump wavefunction to ""wfn[step].gbw"" files every ""stride"" steps"
    write(c80tmp,*) ncenter-1
    if (iORCAver==0) then
        write(ifileid,"(a)") " constraint add center 0.."//trim(adjustl(c80tmp))//"  # Fix center of mass at initial position"
        write(ifileid,"(a)") " run 2000  # Number of MD steps"
    else if (iORCAver==1) then
        write(ifileid,"(a)") " run 2000 CenterCOM # Number of MD steps. Remove motion of center of mass"
    end if
    if (ifixHbond==1) then
        do iatm=1,ncenter
            if (a(iatm)%index==1) then
                do jatm=1,ncenter
                    if (iatm==jatm) cycle
                    if (a(iatm)%index==1.and.a(jatm)%index==1.and.jatm<iatm) cycle !Avoiding double counting if H-H bond exists
                    iaddcons=0
                    dist=atomdist(iatm,jatm,0)
                    if (allocated(connmat)) then
                        if (connmat(iatm,jatm)/=0) iaddcons=1
                    else
                        if ( dist < ( covr(a(iatm)%index)+covr(a(jatm)%index) )*bondcrit) iaddcons=1
                    end if
                    if (iaddcons==1) write(ifileid,"(a,2i6,a,f9.5,'_A')") " constraint add distance",iatm-1,jatm-1,"  target",dist*b2a
                end do
            end if
        end do
    end if
    write(ifileid,"(a)") "end"
end if
if (ilevel==20) then
    write(ifileid,"(a)") "%tddft"
    write(ifileid,"(a)") "  Mode sTDDFT"
    write(ifileid,"(a)") "  Ethresh 7.0"
    write(ifileid,"(a)") "  PThresh 1e-4"
    write(ifileid,"(a)") "  PTLimit 30"
    write(ifileid,"(a)") "  triplets false" !When "triplets true", more excited states will be calculated
    write(ifileid,"(a)") "  maxcore 6000" !sTDDFT only support serial mode, therefore more memory could be assigned
    write(ifileid,"(a)") "end"
else if (ilevel==21) then
    write(ifileid,"(a)") "%tddft"
    write(ifileid,"(a)") "  nroots 10"
    write(ifileid,"(a)") "  mode riints_disk"
    write(ifileid,"(a)") "end"
else if (ilevel==22.or.ilevel==23.or.ilevel==231) then
    write(ifileid,"(a)") "%tddft"
    write(ifileid,"(a)") "  nroots 10"
    write(ifileid,"(a)") "  TDA false"
    write(ifileid,"(a)") "end"
else if (ilevel==24.or.ilevel==25) then
    write(ifileid,"(a)") "%mdci"
    write(ifileid,"(a)") "  nroots 3"
    write(ifileid,"(a)") "end"
end if

!Set frozen
if ((itask==2.or.itask==4.or.itask==5).and.(nfrozenatm>0.or.ioptHonly==1)) then
    write(ifileid,"(a)") "%geom"
    if (ioptHonly==0) then
        write(ifileid,"(a)") "  Constraints"
        do ifr=1,nfrozenatm
            write(ifileid,"('  { C',i7,' C }')") frozenatm(ifr)-1
        end do
        write(ifileid,"(a)") "  end"
    else
        write(ifileid,"(a)") "  optimizeHydrogens true"
    end if
    write(ifileid,"(a)") "end"
end if
!Wavefunction stability test, external field
if (itask==9.or.any(efield/=0)) then
    write(ifileid,"(a)") "%scf"
    if (itask==9) then
        write(ifileid,"(a)") "  STABPerform true"
        write(ifileid,"(a)") "  STABRestartUHFifUnstable true"
    end if
    if (any(efield/=0)) write(ifileid,"('  efield',f10.6,',',f10.6,',',f10.6)") efield(1),efield(2),efield(3)
    write(ifileid,"(a)") "end"
end if


!Write geometry
write(ifileid,"('* xyz',2i4)") netcharge,multi
do i=1,ncenter
    if (atmbasname(i)==" ") then
        	write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
    else    
        	write(ifileid,"(a,1x,3f14.8,a)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a," newGTO """//trim(atmbasname(i))//""" end"
    end if
end do
write(ifileid,"(a)") "*"
close(ifileid)
write(*,"(a)") " ORCA input file has been exported to "//trim(outname)//", you should properly check it and modify keywords"
end subroutine

!Write external DFT-D3(BJ) parameter for wB97X-2
subroutine ORCA_DFT_D3_parm(ifileid)
integer ifileid
write(ifileid,"(a)") "# DFT-D3(BJ) parameter for wB97X-2, see SI of PCCP, 20, 23175 (2018)"
write(ifileid,"(a)") "%method"
write(ifileid,"(a)") "  D3S6 0.547"
write(ifileid,"(a)") "  D3A1 3.520"
write(ifileid,"(a)") "  D3S8 0.0"
write(ifileid,"(a)") "  D3A2 7.795"
write(ifileid,"(a)") "end"
end subroutine
!Write functional definition for revDSD-PBEP86-D4. This is no longer needed, as ORCA 6 already supports revDSD-PBEP86-D4/2021 keyword
!See https://www.compchem.me/revdsd-pbep86-functional
!subroutine ORCA_revDSD_PBEP86_D4_parm(ifileid)
!integer ifileid
!write(ifileid,"(a)") "%method"
!write(ifileid,"(a)") "   Exchange X_PBE"
!write(ifileid,"(a)") "   Correlation C_P86"
!write(ifileid,"(a)") "   ScalHFX 0.69"
!write(ifileid,"(a)") "   ScalDFX 0.31"
!write(ifileid,"(a)") "   ScalGGAC 0.4210"
!write(ifileid,"(a)") "   ScalLDAC 0.4210"
!write(ifileid,"(a)") "   ScalMP2C 1.0"
!write(ifileid,"(a)") "   LDAOpt C_VWN5"
!write(ifileid,"(a)") "   D3S6 0.5132"
!write(ifileid,"(a)") "   D3S8 0.0"
!write(ifileid,"(a)") "   D3A1 0.44"
!write(ifileid,"(a)") "   D3A2 3.60"
!write(ifileid,"(a)") "end"
!write(ifileid,"(a)") "%mp2"
!write(ifileid,"(a)") "   DoSCS true"
!write(ifileid,"(a)") "   PS 0.5922"
!write(ifileid,"(a)") "   PT 0.0636"
!write(ifileid,"(a)") "end"
!end subroutine



!!---------- Output current coordinate to NWChem input file
subroutine outNWCheminp(outname,ifileid)
use defvar
character(len=*) outname
open(ifileid,file=outname,status="replace")
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
write(ifileid,"(a,/)") "start"
write(ifileid,"(a,i3)") "charge",netcharge
write(ifileid,"(a)") "GEOMETRY"
do i=1,ncenter
	write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,"(a)") "END"
write(ifileid,"(a)") "BASIS"
write(ifileid,"(a)") "  * library 6-31G*"
write(ifileid,"(a)") "END"
write(ifileid,"(a)") "DFT"
write(ifileid,"(a)") "  XC B3LYP"
write(ifileid,"(a,i3)") "  mult",nint(naelec-nbelec)+1
write(ifileid,"(a)") "END"
write(ifileid,"(a)") "TASK DFT ENERGY"
close(ifileid)
write(*,"(a)") " Exporting NWChem input file finished! It corresponds to single point task at B3LYP/6-31G* level"
end subroutine



!!---------- Interface of outputting MOPAC input file
subroutine outMOPACinp_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for generating MOPAC input file, e.g. C:\ltwd.mop"
write(*,"(a)") " If press ENTER button directly, will be exported to "//trim(c200tmp)//".mop"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".mop"
call outMOPACinp(outname,10)
end subroutine
!!---------- Output current coordinate to MOPAC input file
subroutine outMOPACinp(outname,ifileid)
use defvar
use util
character(len=*) outname
character multistr*80,chargestr*10,c80tmp*80,c2000tmp*2000,keyword*80
integer :: itask=2,isolv=0,ifopt(ncenter),frozenatm(ncenter),iMOZYME=0
real*8 :: eps
write(*,"(a)") " Note: You can select an option with positive index to generate input file with corresponding level, &
&other options can be used to change various calculation setting"

ifopt(:)=1
do while(.true.)
    write(*,*)
    if (itask==2) write(*,"(a,', current:',i6' atoms')") " -3 Set atom freeze in optimization",count(ifopt==0)
    if (iMOZYME==0) write(*,*) "-2 Toggle employing MOZYME to accelerate calculation, current: No"
    if (iMOZYME==1) write(*,*) "-2 Toggle employing MOZYME to accelerate calculation, current: Yes"
    if (isolv==1) write(*,"(a,f8.2)") " -1 Toggle employing COSMO solvation model, current: Yes, eps=",eps
    if (isolv==0) write(*,*) "-1 Toggle employing COSMO solvation model, current: No"
    if (itask==1) write(*,*) " 0 Select task, current: Single point"
    if (itask==2) write(*,*) " 0 Select task, current: Optimizing minimum"
    if (itask==3) write(*,*) " 0 Select task, current: Frequency analysis"
    write(*,*) " 1 PM6"
    write(*,*) " 2 PM6-DH+"
    write(*,*) " 3 PM6-D3H4X"
    write(*,*) " 4 PM7"
    read(*,*) isel
    
    if (isel>0) then
        exit
    else if (isel==0) then
        write(*,*) "1 Single point"
        write(*,*) "2 Optimizing minimum"
        write(*,*) "3 Frequency analysis"
        read(*,*) itask
    else if (isel==-1) then
        if (isolv==0) then
            write(*,*) "Input dielectric constant, e.g. 60.2"
            write(*,"(a)") " Note: If you press ENTER button directly, water (eps=78.3553) will be employed as solvent"
            read(*,"(a)") c80tmp
            if (c80tmp==" ") then
                eps=78.3553
            else
                read(c80tmp,*) eps
            end if
            isolv=1
        else
            isolv=0
        end if
    else if (isel==-2) then
        if (iMOZYME==0) then
            iMOZYME=1
        else
            iMOZYME=0
        end if
    else if (isel==-3) then
        write(*,*) "Input the atoms to be frozen during optimization, e.g. 3,5-10,13,15"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp,frozenatm)
        ifopt(frozenatm(1:ntmp))=0
    end if
end do

open(ifileid,file=outname,status="replace")

if (isel==1) keyword="PM6"
if (isel==2) keyword="PM6-DH+"
if (isel==3) keyword="PM6-D3H4X"
if (isel==4) keyword="PM7"
if (iMOZYME==1) keyword=trim(keyword)//" MOZYME"

if (itask==1) then
    ifopt=0
else if (itask==2) then
    keyword=trim(keyword)//" PRNT=2"
else if (itask==3) then
    keyword=trim(keyword)//" FORCE LET"
end if

if (isolv==1) then
    write(c80tmp,"(f12.4)") eps
    keyword=trim(keyword)//" eps="//adjustl(trim(c80tmp))
end if

netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
write(chargestr,"(i4)") netcharge
chargestr=adjustl(chargestr)
multival=nint(naelec-nbelec)+1

if (multival==1) then
	keyword=trim(keyword)//" precise charge="//trim(chargestr)
else
	if (multival==2) multistr="doublet"
	if (multival==3) multistr="triplet"
	if (multival==4) multistr="quartet"
	if (multival==5) multistr="quintet"
	if (multival==6) multistr="sextet"
	if (multival==7) multistr="septet"
	if (multival==8) multistr="octet"
	if (multival==9) multistr="nonet"
    keyword=trim(keyword)//" precise UHF charge="//trim(chargestr)//' '//multistr
end if

write(ifileid,"(a)") trim(keyword)
write(ifileid,"(a)") "molecule"
write(ifileid,"(a)") "All coordinates are Cartesian"

do i=1,ncenter
	write(ifileid,"(a,1x,f14.8,i2,f14.8,i2,f14.8,i2)") a(i)%name,a(i)%x*b2a,ifopt(i),a(i)%y*b2a,ifopt(i),a(i)%z*b2a,ifopt(i)
end do

close(ifileid)
write(*,*) "Exporting MOPAC input file finished!"
end subroutine



!!---------- Interface of outputting PSI4 input file
subroutine outPSI4inp_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,"(/,a)") " Note: If you benefit from this function when creating PSI4 input file, please cite Multiwfn original paper, thanks!"
write(*,*)
write(*,*) "Input path for generating PSI4 input file, e.g. C:\ltwd.inp"
write(*,"(a)") " If press ENTER button directly, will be exported to "//trim(c200tmp)//".inp"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".inp"
call outPSI4inp(outname,10)
end subroutine
!!---------- Output current coordinate to PSI4 input file
subroutine outPSI4inp(outname,ifileid)
use defvar
use util
character(len=*) outname
integer ifileid,atmlist1(ncenter),atmlist2(ncenter),ioutfch
character c200tmp*200,c2000tmp*2000,taskstr*80,methodstr*20
ioutfch=0
write(*,*)
write(*,*) "-1 Generate input file using a template input file" 
write(*,*) "0 Return"
write(*,*) "1 Intermolecular SAPT energy decomposition task"
write(*,*) "2 Single point energy or excitation energy task"
write(*,*) "3 Optimization task"
read(*,*) isel
if (isel==0) return

open(ifileid,file=outname,status="replace")

if (isel==-1) then !Special case
    write(*,*) "Input the path of the template file, e.g. C:\Sob\love\nico.inp"
    write(*,"(a)") " This file should contain everything of final input file, but the geometry part &
    should be denoted as [geometry] or [GEOMETRY], which will be automatically replaced by current coordinate"
    do while(.true.)
        read(*,"(a)") c200tmp
	    inquire(file=c200tmp,exist=alive)
	    if (alive) exit
	    write(*,*) "Cannot find the file, input again!"
    end do
    open(ifileid+1,file=c200tmp,status="old")
    do while(.true.)
        read(ifileid+1,"(a)",iostat=ierror) c200tmp
        if (ierror/=0) exit
        if (index(c200tmp,"[geometry]")/=0.or.index(c200tmp,"[GEOMETRY]")/=0) then
            do iatm=1,ncenter
	            write(ifileid,"(a,1x,3f14.8)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
            end do
        else
            write(ifileid,"(a)") trim(c200tmp)
        end if
    end do
    close(ifileid)
    close(ifileid+1)
    write(*,"(a)") " Done! The file has been exported to "//trim(outname)
    return
end if

write(ifileid,"('memory 4 GB')")
if (isel==1) then
    write(*,*) "1 Bronze (sSAPT0/jun-cc-pVDZ)"
    write(*,*) "2 Silver (SAPT2+/aug-cc-pVDZ)"
    write(*,*) "3 Gold (SAPT2+(3)dMP2/aug-cc-pVTZ)"
    write(*,*) "10 SAPT0/aug-cc-pVDZ with explicit charge-transfer energy"
    write(*,*) "11 SAPT2+(3)/aug-cc-pVTZ with explicit charge-transfer energy"
    write(*,*) "20 SAPT(DFT)/aug-cc-pVTZ"
    read(*,*) isel2
    
    write(*,*) "Input indices of the atoms in fragment 1, e.g. 3,8-15,20 "
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,nlist1,atmlist1)
    write(*,*) "Input indices of the atoms in fragment 2, e.g. 3,8-15,20 "
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,nlist2,atmlist2)
    
    write(ifileid,"(a)")
    write(ifileid,"(a)") "molecule dimer {"
    write(ifileid,"(a)") "0 1"
    do idx=1,nlist1
        iatm=atmlist1(idx)
	    write(ifileid,"(a,1x,3f14.8)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
    end do
    write(ifileid,"(a)") "--"
    write(ifileid,"(a)") "0 1"
    do idx=1,nlist2
        iatm=atmlist2(idx)
	    write(ifileid,"(a,1x,3f14.8)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
    end do
    write(ifileid,"(a)") "}"
    write(ifileid,"(a)")
    write(ifileid,"(a)") "set {"
    write(ifileid,"(a)") "    scf_type DF"
    write(ifileid,"(a)") "    freeze_core True"
    if (isel2==1) write(ifileid,"(a)") "    basis jun-cc-pVDZ"
    if (isel2==2) write(ifileid,"(a)") "    basis aug-cc-pVDZ"
    if (isel2==3) write(ifileid,"(a)") "    basis aug-cc-pVTZ"
    if (isel2==10) write(ifileid,"(a)") "    basis aug-cc-pVDZ"
    if (isel2==11) write(ifileid,"(a)") "    basis aug-cc-pVTZ"
    if (isel2==20) then
        write(ifileid,"(a)") "    basis aug-cc-pVTZ"
        write(ifileid,"(a)") "    SAPT_DFT_FUNCTIONAL PBE0"
        write(ifileid,"(a)") "    sapt_dft_grac_shift_a 0.03"
        write(ifileid,"(a)") "    sapt_dft_grac_shift_b 0.05"
        write(*,"(a)") " Note: You should set the functional name after ""SAPT_DFT_FUNCTIONAL"" to the actual one you want to use. &
        The ""sapt_dft_grac_shift_a"" and ""sapt_dft_grac_shift_b"" are shift parameters for monomer 1 and monomer 2, respectively, &
        they should be manually set (in Hartree) and evaluated as E(HOMO)+VIP, where E(HOMO) is the HOMO energy of original state, and VIP is vertical ionization potential, &
        which is calculated as E(N-1) - E(N) with N being the number of electrons of original state"
    end if
    write(ifileid,"(a)") "}"
    write(ifileid,"(a)")
    if (isel2==1) write(ifileid,"(a)") "energy('sapt0')"
    if (isel2==2) write(ifileid,"(a)") "energy('sapt2+')"
    if (isel2==3) write(ifileid,"(a)") "energy('sapt2+(3)dmp2')"
    if (isel2==10) write(ifileid,"(a)") "energy('sapt0-ct')"
    if (isel2==11) write(ifileid,"(a)") "energy('sapt2+(3)-ct')"
    if (isel2==20) write(ifileid,"(a)") "energy('sapt(dft)')"
    write(*,"(a)") " Exporting PSI4 input file finished! Note that charge of spin multiplicity &
    of the fragments are assumed to be 0 and 1, respectively"

else if (isel==2.or.isel==3) then !Single point, optimization
    !Note: At least for PSI 1.3.2, only DF-CCSD(T) gradient can use frozen core
    do while(.true.)
        write(*,*)
        if (isel==2) then
            if (ioutfch==0) write(*,*) "-1 Toggle outputting .fch file, current: No"
            if (ioutfch==1) write(*,*) "-1 Toggle outputting .fch file, current: Yes"
        end if
        write(*,*) "                         Select a calculation level:"
        if (isel==2) then !SP
            write(*,*) " 1 BLYP-D3(BJ)/def2-SVP with density fitting"
            write(*,*) " 2 B3LYP-D3(BJ)/def2-SVP with density fitting"
            write(*,*) " 3 PBEh-3c with density fitting"
            write(*,*) " 4 M06-2X/def2-TZVP with density fitting"
            write(*,*) " 5 PWPB95-D3(BJ)/cc-pVTZ with density fitting and frozen core"
            write(*,*) " 9 MP2/cc-pVTZ with density fitting and frozen core"
            write(*,*) "10 CCSD(T)/cc-pVTZ with frozen core"
            write(*,*) "11 CCSD(T)/cc-pVTZ with density fitting and frozen core"
            write(*,*) "12 Full-CI/cc-pVDZ"
            write(*,*) "19 ADC(2)/cc-pVTZ with frozen core"
            write(*,*) "20 EOM-CC2/cc-pVTZ with frozen core"
            write(*,*) "21 EOM-CCSD/cc-pVTZ with frozen core"
            write(*,*) "22 EOM-CC3/cc-pVTZ with frozen core"
        else if (isel==3) then !opt
            write(*,*) " 1 BLYP-D3(BJ)/def2-SVP with density fitting"
            write(*,*) " 2 B3LYP-D3(BJ)/def2-SVP with density fitting"
            write(*,*) " 3 PBEh-3c with density fitting"
            write(*,*) " 4 M06-2X/def2-TZVP with density fitting"
            write(*,*) " 5 PWPB95-D3(BJ)/cc-pVTZ with density fitting and frozen core (Numer. grad.)"
            write(*,*) " 9 MP2/cc-pVTZ with density fitting and frozen core"
            write(*,*) "10 CCSD(T)/cc-pVTZ without frozen core"
            write(*,*) "11 CCSD(T)/cc-pVTZ with density fitting and frozen core"
            write(*,*) "12 Full-CI/cc-pVDZ (Numerical gradient)"
            !write(*,*) "19 ADC(2)/cc-pVTZ with frozen core" !opt has not been supported, even numerical gradient
            write(*,*) "20 EOM-CC2/cc-pVTZ with frozen core (Numer. grad.)"
            write(*,*) "21 EOM-CCSD/cc-pVTZ" !EOM-CCSD gradient doesn't support frozen core
            write(*,*) "22 EOM-CC3/cc-pVTZ with frozen core (Numer. grad.)"
        end if
        read(*,*) isel2
        if (isel2==-1) then
            if (ioutfch==0) then
                ioutfch=1
            else
                ioutfch=0
            end if
        else
            exit
        end if
    end do
    
    netcharge=nint(sum(a%charge)-nelec)
    if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
    multi=nint(naelec-nbelec)+1
    
    !Write molecular geometry
    write(ifileid,"(a)") "molecule mol {"
    write(ifileid,"(2i4)") netcharge,multi
    do iatm=1,ncenter
	    write(ifileid,"(a,1x,3f14.8)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
    end do
    if (isel2==20.or.isel2==21.or.isel==22) write(ifileid,"(a)") "symmetry c1" !For EOM- calculation, assumed to be nosymm
    write(ifileid,"(a)") "}"
    
    !Set basis set and various settings
    write(ifileid,*)
    write(ifileid,"(a)") "set {"
    if (isel2==1.or.isel2==2) then
        write(ifileid,"(2x,a)") "basis def2-SVP"
    else if (isel2==4) then
        write(ifileid,"(2x,a)") "basis def2-TZVP"
    else if (isel2==5.or.isel2==9) then
        write(ifileid,"(2x,a)") "basis cc-pVTZ"
        write(ifileid,"(2x,a)") "freeze_core true"
        write(ifileid,"(2x,a)") "scf_type df"
    else if (isel2==10.or.isel2==11) then
        write(ifileid,"(2x,a)") "basis cc-pVTZ"
        if (.not.(isel==3.and.isel2==10)) write(ifileid,"(2x,a)") "freeze_core true" !Regular CCSD(T) cannot use FC
        if (isel2==11) then
            write(ifileid,"(2x,a)") "scf_type df"
            write(ifileid,"(2x,a)") "cc_type df"
        end if
    else if (isel2==12) then
        write(ifileid,"(2x,a)") "basis cc-pVDZ"
        write(ifileid,"(2x,a)") "#NUM_ROOTS 3"
    else if (isel2==20.or.isel2==21.or.isel2==22) then
        write(ifileid,"(2x,a)") "basis cc-pVTZ"
        if (.not.(isel==3.and.isel2==21)) write(ifileid,"(2x,a)") "freeze_core true"
        write(ifileid,"(2x,a)") "ROOTS_PER_IRREP [3]"
    end if
    if (multi>1) write(ifileid,"(2x,a)") "reference uhf"
    write(ifileid,"(a)") "}"
    
    !Method
    write(ifileid,*)
    if (isel==2) taskstr="energy"
    if (isel==3) taskstr="opt"
    if (isel2==1) then
        methodstr="blyp-d3bj"
    else if (isel2==2) then
        methodstr="b3lyp-d3bj"
    else if (isel2==3) then
        methodstr="pbeh3c"
    else if (isel2==4) then
        methodstr="m06-2x"
    else if (isel2==5) then
        methodstr="pwpb95-d3bj"
    else if (isel2==9) then
        methodstr="mp2"
    else if (isel2==10.or.isel2==11) then
        methodstr="ccsd(t)"
    else if (isel2==12) then
        methodstr="fci"
    else if (isel2==20) then
        methodstr="eom-cc2"
    else if (isel2==21) then
        methodstr="eom-ccsd"
    else if (isel2==22) then
        methodstr="eom-cc3"
    end if
    if (ioutfch==0) then
        write(ifileid,"(a)") trim(taskstr)//"('"//trim(methodstr)//"')"
    else if (ioutfch==1) then
        if (isel2<5) then
            write(ifileid,"(a)") "energy, wfn = "//trim(taskstr)//"('"//trim(methodstr)//"', return_wfn=True)"
        else
            write(ifileid,"(a)") "grad, wfn = gradient"//"('"//trim(methodstr)//"', return_wfn=True)"
        end if
        write(ifileid,"(a)") "fchk_writer = psi4.FCHKWriter(wfn)"
        itmp=index(outname,'.',back=.true.)
        if (itmp==0) itmp=len_trim(outname)+1
        write(ifileid,"(a)") "fchk_writer.write('"//outname(:itmp-1)//".fch')"
    end if
    write(*,"(a)") " Exporting PSI4 input file finished!"
end if

close(ifileid)
end subroutine



!!---------- Output current coordinate to MRCC input file
subroutine outMRCCinp(outname,ifileid)
use defvar
character(len=*) outname
open(ifileid,file=outname,status="replace")
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
write(ifileid,"(a)") "basis=cc-pVDZ"
write(ifileid,"(a)") "calc=CCSDT(Q)"
write(ifileid,"(a,i2)") "charge=",netcharge
write(ifileid,"(a,i2)") "mult=",nint(naelec-nbelec)+1
write(ifileid,"(a)") "mem=2500MB"
write(ifileid,*)
write(ifileid,"(a)") "geom=xyz"
write(ifileid,"(i5)") ncenter
write(ifileid,*)
do i=1,ncenter
	write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
close(ifileid)
write(*,"(a)") " Exporting MRCC input file finished! The file is named ""MINP"" in current folder. It corresponds to single point task at CCSDT(Q)/cc-pVDZ level"
end subroutine


!!---------- Output current coordinate to CFOUR input file
subroutine outCFOURinp(outname,ifileid)
use defvar
character(len=*) outname
open(ifileid,file=outname,status="replace")
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
write(ifileid,"(a)") "mol"
do i=1,ncenter
	write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,*)
write(ifileid,"(a,i3)") "CHARGE=",netcharge
write(ifileid,"(a,i3)") "MULTIPLICITY=",nint(naelec-nbelec)+1
write(ifileid,"(a)") "FROZEN_CORE=ON"
write(ifileid,"(a)") "*CFOUR(CALC=CCSD(T),BASIS=cc-pVTZ,COORD=CARTESIAN,CC_PROG=NCC,ABCDTYPE=AOBASIS)"
write(ifileid,*)
close(ifileid)
write(*,"(a)") " Exporting CFOUR input file finished! The file is named ""ZMAT"" in current folder. It corresponds to single point task at CCSD(T,FC)/cc-pVTZ level"
end subroutine


!!---------- Output current coordinate to Molpro input file
subroutine outMolproinp(outname,ifileid)
use defvar
character(len=*) outname
open(ifileid,file=outname,status="replace")
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
write(ifileid,"(a)") "angstrom"
write(ifileid,"(a)") "geometry={"
do i=1,ncenter
	write(ifileid,"(a,',,',f14.8,',',f14.8,',',f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,"(a)") "}"
write(ifileid,"(a)") "basis=vtz"
write(ifileid,"(a)") "{hf"
write(ifileid,"(a,i2,a,i2)") "wf, charge=",netcharge,", spin=",nint(naelec-nbelec)
write(ifileid,"(a)") "}"
write(ifileid,"(a)") "ccsd(T)"
close(ifileid)
write(*,"(a)") " Exporting Molpro input file finished! It corresponds to single point task at CCSD(T)/cc-pVTZ level"
end subroutine


!!---------- Output current coordinate to Molcas input file
subroutine outmolcasinp(outname,ifileid)
use defvar
character(len=*) outname
open(ifileid,file=outname,status="replace")
netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
write(ifileid,"(a)") "&GATEWAY"
write(ifileid,"(a)") "Coord"
write(ifileid,"(i5)") ncenter
write(ifileid,"(a)") "Generated by Multiwfn"
do i=1,ncenter
	write(ifileid,"(a,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,"(a)") "basis=cc-pVDZ"
write(ifileid,"(a)") "&SEWARD"
write(ifileid,"(a)") "&SCF"
if (naelec==nbelec) then
	write(ifileid,"(a)") "ksdft=b3lyp"
else
	write(ifileid,"(a)") "uhf;ksdft=b3lyp"
end if
write(ifileid,"(a,i3)") "charge=",netcharge
write(ifileid,"(a,i3)") "zspin=",nint(naelec-nbelec)
close(ifileid)
write(*,"(a)") " Exporting Molcas input file finished! It corresponds to single point task at B3LYP/cc-pVDZ level"
end subroutine


!!---------- Output current coordinate to Q-Chem input file
subroutine outQcheminp(outname,ifileid)
use defvar
character(len=*) outname
open(ifileid,file=outname,status="replace")
netcharge=nint(sum(a%charge)-nelec)
multi=nint(naelec-nbelec)+1
write(ifileid,"(a)") "$molecule"
write(ifileid,"(2i4)") netcharge,multi
do i=1,ncenter
	write(ifileid,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,"(a)") "$end"
write(ifileid,*)
write(ifileid,"(a)") "$rem"
write(ifileid,"(a)") "GUI 2    !Generate .fch file"
write(ifileid,"(a)") "EXCHANGE B3LYP"
write(ifileid,"(a)") "BASIS 6-31G*"
write(ifileid,"(a)") "DFT_D D3_BJ"
write(ifileid,"(a)") "$end"
close(ifileid)
write(*,"(a)") " Exporting Q-Chem input file finished! It corresponds to single point task at B3LYP-D3(BJ)/6-31G* level"
end subroutine


!!---------- Output current coordinate to Dalton input file (.dal and .mol)
!When symmetry is used, the number of atomtypes may not be consistent with the number of elements, such as O3.
!Therefore, for simplicity, we do not use symmetry (i.e. Nosymmetry), and thus assign atomtypes as the number of elements
!Spin multiplicity is not set in .dal, because the rule in Dalton is quite complicated
subroutine outDaltoninp(dalname,molname,ifileid)
use defvar
character(len=*) dalname,molname
character tmpstr*5,c20tmp*20,c20tmp2*20,selectyn
real*8 tmparr(nbasis),CObasa_reorder(nbasis,nbasis),CObasa_reordertmp(nbasis,nbasis)
integer ioutfmt
integer reorder_d_sph(5),reorder_f_sph(7),reorder_g_sph(9),reorder_h_sph(11)
integer reorder_d_car(6),reorder_f_car(10),reorder_g_car(15),reorder_h_car(21)

netcharge=nint(sum(a%charge)-nelec)
if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
if (dalname/=" ") then
	if (allocated(CObasa)) then
		write(*,*) "Basis function information is currently available, do you want to also write orbital expansion coefficients &
        into the .dal file, which can then be used as initial guess? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y'.or.selectyn=='Y') then
			!Although in the following code, the order of Cartesian shell is properly converted, finally I decide not to support
			!Because for >=D shell, the normalization rule of Cartesian basis function in Dalton is strange, which alters the coefficients that should be written
			if (any(bastype(1:nbasis)>=5)) then
				write(*,"(a)") " Error: Some basis functions are Cartesian type, however, only spherical-harmonic basis functions are supported by this function"
				write(*,*) "Press ENTER button to return"
				read(*,*)
				return
			end if
			write(*,*) "Choose the format of the coefficients to write:"
            write(*,"(a)") " 1 Using 4F18.14 format, which can be loaded without modifying Dalton source code; but &
            if a coefficient is very large, it cannot be properly written and loaded."
            write(*,"(a)") " 2 Using 4E20.12, which can properly record any coefficient, but you need to &
            modify ""DALTON/sirius/sirgp.F"" in DALTON source code folder, replacing ""PFMT = '(4F18.14)'"" with ""PFMT = '(4E20.12)'"", and then recompile."
            read(*,*) ioutfmt
			!Note: Order of pure shells of Dalton is D-2/-1/0/+1/+2, etc. Need conversion to CObasa_reorder
			reorder_d_sph(1:5)=(/ 5,3,1,2,4 /) !Mapping D0/+1/-1/+2/-2 to D-2/-1/0/+1/+2
			reorder_f_sph(1:7)=(/ 7,5,3,1,2,4,6 /) !Mapping F0/+1/-1/+2/-2/+3/-3 to F-3/-2/-1/0/+1/+2/+3
			reorder_g_sph(1:9)=(/ 9,7,5,3,1,2,4,6,8 /) !Mapping G0/+1/-1/+2/-2/+3/-3/+4/-4 to G-4/-3/-2/-1/0/+1/+2/+3/+4
			reorder_h_sph(1:11)=(/ 11,9,7,5,3,1,2,4,6,8,10 /) !Mapping H0/+1/-1/+2/-2/+3/-3/+4/-4/+5/-5 to G-5/-4/-3/-2/-1/0/+1/+2/+3/+4/+5
            !Mapping XX,YY,ZZ,XY,XZ,YZ to &
            !        XX,XY,XZ,YY,YZ,ZZ
            reorder_d_car(1:6)=(/ 1,4,5,2,6,3 /) 
            !Mapping XXX,YYY,ZZZ,XXY,XXZ, YYZ,XYY,XZZ,YZZ,XYZ to &
            !        XXX,XXY,XXZ,XYY,XYZ, XZZ,YYY,YYZ,YZZ,ZZZ
            reorder_f_car(1:10)=(/ 1,4,5,7,10,8,2,6,9,3 /) 
            !Mapping ZZZZ,YZZZ,YYZZ,YYYZ,YYYY, XZZZ,XYZZ,XYYZ,XYYY,XXZZ, XXYZ,XXYY,XXXZ,XXXY,XXXX to &
            !        XXXX,XXXY,XXXZ,XXYY,XXYZ, XXZZ,XYYY,XYYZ,XYZZ,XZZZ, YYYY,YYYZ,YYZZ,YZZZ,ZZZZ
            reorder_g_car(1:15)=(/ 15,14,13,12,11, 10,9,8,7,6, 5,4,3,2,1 /) 
            !Mapping ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ, YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ, XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY, XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY, XXXXX to &
            !        XXXXX,XXXXY,XXXXZ,XXXYY,XXXYZ, XXXZZ,XXYYY,XXYYZ,XXYZZ,XXZZZ, XYYYY,XYYYZ,XYYZZ,XYZZZ,XZZZZ, YYYYY,YYYYZ,YYYZZ,YYZZZ,YZZZZ, ZZZZZ
            reorder_h_car(1:21)=(/ 21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1 /) 
			ibas=1
			do while(.true.)
				if (bastype(ibas)>=0) then !S,P
                    CObasa_reorder(ibas,:)=CObasa(ibas,:)
                    ibas=ibas+1
                else if (bastype(ibas)>=-5.and.bastype(ibas)<=-1) then !D spherical
					do itmp=1,5
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_d_sph(itmp)-1,:)
                    end do
                    ibas=ibas+5
                else if (bastype(ibas)>=5.and.bastype(ibas)<=10) then !D Cartesian
					do itmp=1,6
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_d_car(itmp)-1,:)
                    end do
                    ibas=ibas+6
                else if (bastype(ibas)>=-12.and.bastype(ibas)<=-6) then !F spherical
					do itmp=1,7
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_f_sph(itmp)-1,:)
                    end do
                    ibas=ibas+7
                else if (bastype(ibas)>=11.and.bastype(ibas)<=20) then !F Cartesian
					do itmp=1,10
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_f_car(itmp)-1,:)
                    end do
                    ibas=ibas+10
                else if (bastype(ibas)>=-21.and.bastype(ibas)<=-13) then !G spherical
					do itmp=1,9
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_g_sph(itmp)-1,:)
                    end do
                    ibas=ibas+9
                else if (bastype(ibas)>=21.and.bastype(ibas)<=35) then !G Cartesian
					do itmp=1,15
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_g_car(itmp)-1,:)
                    end do
                    ibas=ibas+15
                else if (bastype(ibas)>=-32.and.bastype(ibas)<=-22) then !H spherical
					do itmp=1,11
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_h_sph(itmp)-1,:)
                    end do
                    ibas=ibas+11
                else if (bastype(ibas)>=36.and.bastype(ibas)<=56) then !H Cartesian
					do itmp=1,21
						CObasa_reorder(ibas+itmp-1,:)=CObasa(ibas+reorder_h_car(itmp)-1,:)
                    end do
                    ibas=ibas+21
                end if
                if (ibas>nbasis) exit
            end do
            !Reorder shell according to angular moment. In Dalton, shells in each atom are always in S, P, D, F... sequence
            CObasa_reordertmp=CObasa_reorder
            itmp=0
            do iatm=1,ncenter !Loop atoms
				do iang=0,5 !Loop angular moment from S to H
					do ibas=basstart(iatm),basend(iatm)
						itype=bastype(ibas)
						if ((iang==0.and.itype==1) .or. (iang==1.and.itype>=2.and.itype<=4) .or. (iang==2.and.itype>=-5.and.itype<=-1) &
                        .or. (iang==3.and.itype>=-12.and.itype<=-6) .or. (iang==4.and.itype>=-21.and.itype<=-13) .or. (iang==5.and.itype>=-32.and.itype<=-22)) then
                        itmp=itmp+1
                        CObasa_reorder(itmp,:)=CObasa_reordertmp(ibas,:)
                        end if
                    end do
				end do
            end do
            if (wfntype==1.or.wfntype==4) then
				write(*,"(/,a)") " Warning: Dalton is unable to perform unrestricted calculation, while present wavefunction is unrestricted. &
				&Only alpha orbitals are written to .dal file"
            end if
        end if
    end if
	open(ifileid,file=dalname,status="replace")
	write(ifileid,"(a)") "**DALTON INPUT"
	write(ifileid,"(a)") ".RUN WAVE FUNCTIONS"
	write(ifileid,"(a)") "**WAVE FUNCTIONS"
	write(ifileid,"(a)") ".DFT"
	write(ifileid,"(a)") " B3LYPg"
    if (selectyn=='y'.or.selectyn=='Y') then
		write(ifileid,"(a)") "*ORBITAL"
		write(ifileid,"(a)") ".MOSTART"
		write(ifileid,"(a)") " FORM18"
		write(ifileid,"(a)") "**MOLORB (Punched by Multiwfn)"
        !Note that the order of basis functions is dependent of occurrence order of atom types, which is considered here
        do imo=1,nbasis
			itmp=0
            ibas=1
			do iele=1,nelesupp
				do iatm=1,ncenter
					if (a(iatm)%index==iele) then
						itmp=itmp+1
                        nbastmp=basend(iatm)-basstart(iatm)+1
						tmparr(ibas:ibas+nbastmp-1)=CObasa_reorder(basstart(iatm):basend(iatm),imo)
                        ibas=ibas+nbastmp
					end if
				end do
			end do
			if (ioutfmt==1) then
				write(ifileid,"(4F18.14)") tmparr
			else
				write(ifileid,"(4E20.12)") tmparr
			end if
        end do
    end if
	write(ifileid,"(a)") "**END OF INPUT"
	close(ifileid)
	write(*,"(a)") " Exporting .dal file finished! It corresponds to single point task using B3LYPg functional"
end if

open(ifileid,file=molname,status="replace")
write(ifileid,"(a)") "ATOMBASIS"
write(ifileid,"(a)") "test molecule"
write(ifileid,"(a)") "Generated by Multiwfn"
natmtype=0
do iele=1,nelesupp
	if (any(a%index==iele)) natmtype=natmtype+1
end do
write(c20tmp,"(i5)") natmtype
c20tmp=adjustl(c20tmp)
write(c20tmp2,"(i5)") netcharge
c20tmp2=adjustl(c20tmp2)
write(ifileid,"(a,a,a,a)") "Atomtypes=",trim(c20tmp)," Angstrom Nosymmetry charge=",trim(c20tmp2)
do iele=1,nelesupp
	natmthis=count(a%index==iele)
	if (natmthis>0) then	
		write(c20tmp,"(f4.1)") dfloat(iele)
		c20tmp=adjustl(c20tmp)
		write(c20tmp2,"(i5)") natmthis
		c20tmp2=adjustl(c20tmp2)
		write(ifileid,"(a,a,a,a,a)") "Charge=",trim(c20tmp)," Atoms=",trim(c20tmp2)," Basis=6-31G*"
		itmp=0
		do iatm=1,ncenter
			if (a(iatm)%index/=iele) cycle
			itmp=itmp+1
			write(tmpstr,"(i5)") itmp
			tmpstr=adjustl(tmpstr)
			write(ifileid,"(a,3f14.8)") trim(a(iatm)%name)//trim(tmpstr),a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
		end do
	end if
end do
close(ifileid)
write(*,*) "Exporting .mol file finished! You need to properly set basis set in this file"
if (naelec/=nbelec) write(*,*) "NOTE: Electronic configuration in the .mol file should be manually set properly"
end subroutine




!!---------- Interface of outputting Quantum ESPRESSO input file
subroutine outQEinp_wrapper
use util
use defvar
character outname*200,c200tmp*200
call path2filename(filename,c200tmp)
write(*,*) "Input path for generating Quantum ESPRESSO input file, e.g. C:\ltwd.inp"
write(*,"(a)") " If press ENTER button directly, will export to "//trim(c200tmp)//".inp"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".inp"
call outQEinp(outname,10)
end subroutine
!!---------- Output current coordinate to Quantum ESPRESSO input file
subroutine outQEinp(outname,ifileid)
use defvar
use util
character(len=*) outname
integer ifileid,eleidx(ncenter)
character c80tmp*80,c200tmp*200,c2000tmp*2000

!Count number of elements
nele=1
eleidx(1)=a(1)%index
do iatm=2,ncenter
    if (all(eleidx(:nele)/=a(iatm)%index)) then
        nele=nele+1
        eleidx(nele)=a(iatm)%index
    end if
end do

open(ifileid,file=outname,status="replace")
call path2filename(filename,c200tmp)
write(ifileid,"(a)") "&control"
write(ifileid,"(a)") " calculation= 'scf'"
write(ifileid,"(a)") " prefix= '"//trim(c200tmp)//"'"
write(ifileid,"(a)") " pseudo_dir= './'"
write(ifileid,"(a)") "/"
write(ifileid,"(a)") "&system"
write(ifileid,"(a)") " ibrav= 0"
write(ifileid,"(a,i6)") " nat=",ncenter
write(ifileid,"(a,i5)") " ntyp=",nele
write(ifileid,"(a)") " ecutwfc= 44.0"
write(ifileid,"(a)") "/"

write(ifileid,"(a)") "&electrons"
write(ifileid,"(a)") "/"
write(ifileid,"(a)") "CELL_PARAMETERS angstrom"
write(ifileid,"(3f14.8)") cellv1(:)*b2a
write(ifileid,"(3f14.8)") cellv2(:)*b2a
write(ifileid,"(3f14.8)") cellv3(:)*b2a
write(ifileid,"(a)") "ATOMIC_SPECIES"
do itype=1,nele
    write(ifileid,"(a,f9.4,a)") ind2name(eleidx(itype)),atmwei(eleidx(itype)),"  Pseudopotential_file"
end do
write(ifileid,"(a)") "ATOMIC_POSITIONS angstrom"
do iatm=1,ncenter
    write(ifileid,"(a,3f12.6)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
end do
write(ifileid,"(a)") "K_POINTS automatic"
write(ifileid,"(a)") "3 3 3 0 0 0"

close(ifileid)

write(*,"(a)") " Quantum ESPRESSO input file has been exported to "//trim(outname)//", you need to manually set proper pseudopotential file"
if (ifPBC==0) write(*,"(a)") " Warning: Since the file loaded when Multiwfn booted up does not contain cell information, &
&all cell vectors in the created input file are set to be 0!"
end subroutine



!!---------- Interface of outputting VASP POSCAR file
subroutine outPOSCAR_wrapper
use util
use defvar
character outname*200,c200tmp*200
if (ifPBC==0) then
	write(*,"(a)") " Error: POSCAR cannot be generated, because the loaded file does not contain cell information!"
    write(*,"(a)") " See Section 2.9.3 of Multiwfn manual on which format can provide cell information to Multiwfn"
    write(*,*) "Press ENTER button to return"
    read(*,*)
	return
end if
call path2filename(filename,c200tmp)
write(*,*) "Input path for generating VASP POSCAR file, e.g. C:\ltwd.poscar"
write(*,"(a)") " If press ENTER button directly, will export to POSCAR in current folder"
read(*,"(a)") outname
if (outname==" ") outname="POSCAR"
call outPOSCAR(outname,10)
write(*,"(a)") " VASP POSCAR file has been exported to "//trim(outname)//" in current folder"
end subroutine
!!---------- Output current coordinate to VASP position file POSCAR
subroutine outPOSCAR(outname,ifileid)
use defvar
use util
character(len=*) outname
integer ifileid,natmele(nelesupp)
character c80tmp*80,c200tmp*200,c2000tmp*2000
real*8 fract(3),Cart(3)

call path2filename(filename,c200tmp)
open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") trim(c200tmp)//"; Created by Multiwfn"
write(ifileid,"(a)") "1.0"
write(ifileid,"(3f16.8)") cellv1*b2a
write(ifileid,"(3f16.8)") cellv2*b2a
write(ifileid,"(3f16.8)") cellv3*b2a
!Count how many atoms corresponding to each element
do iele=1,nelesupp
    natmele(iele)=count(a%index==iele)
end do
!Write element names
do iele=1,nelesupp
    if (natmele(iele)/=0) write(ifileid,"(a6)",advance="no") ind2name(iele)
end do
write(ifileid,*)
!Write number of atoms of each element
do iele=1,nelesupp
   if (natmele(iele)/=0) write(ifileid,"(i6)",advance="no") natmele(iele)
end do
write(ifileid,*)
write(ifileid,"(a)") "Direct"
!Write fractional atom coordinates
do iele=1,nelesupp
    do iatm=1,ncenter
        if (a(iatm)%index==iele) then
            Cart(1)=a(iatm)%x
            Cart(2)=a(iatm)%y
            Cart(3)=a(iatm)%z
            call Cart2fract(Cart,fract)
            write(ifileid,"(3f16.8)") fract(:)
        end if
    end do
end do
close(ifileid)
end subroutine



!!---------- Interface of outputting Gaussian type .cub file
subroutine outcube_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path of the new cube file, e.g. C:\Tree.cub"
write(*,"(a)") " If press ENTER button directly, will be exported to "//trim(c200tmp)//".cub"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".cub"
write(*,*) "Exporting, please wait..."
open(10,file=outname,status="replace")
call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
close(10)
write(*,*) "Done, cube file has been outputted"
end subroutine
!!-------- Output 3D matrix with property to a cube file. fileid must be opened before invoking this routine, and close it after that
!Example of calling: outcube(holegrid,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
subroutine outcube(matrix,numx,numy,numz,org_x,org_y,org_z,trans1,trans2,trans3,fileid)
use defvar
implicit real*8 (a-h,o-z)
integer numx,numy,numz,fileid
real*8 org_x,org_y,org_z,trans1(3),trans2(3),trans3(3)
real*8 matrix(numx,numy,numz)

write(fileid,"(' Generated by Multiwfn')")
write(fileid,"(' Totally ',i12,' grid points')") numx*numy*numz
if (ncenter>=1) then
	write(fileid,"(i5,3f12.6)") ncenter,org_x,org_y,org_z
else
	write(fileid,"(i5,3f12.6)") 1,org_x,org_y,org_z
end if
write(fileid,"(i5,3f12.6)") numx,trans1
write(fileid,"(i5,3f12.6)") numy,trans2
write(fileid,"(i5,3f12.6)") numz,trans3
if (ncenter>=1) then
	do i=1,ncenter
		write(fileid,"(i5,4f12.6)") a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z
	end do
else
	write(*,"(a)") " Note: Current system has no atom, in order to maximize compatibility of the generated .cub file, a hydrogen atom is added to 0,0,0"
	write(*,*)
	write(fileid,"(i5,4f12.6)") 1,1D0,0D0,0D0,0D0
end if
do i=1,numx
	do j=1,numy
		write(fileid,"(6(1PE14.5E3))",advance="no") matrix(i,j,1:numz)
		write(fileid,*)
	end do
end do
end subroutine



!!!------------- Output grid data to .vti file, which can be visualized by ParaView
subroutine outvti(outvtiname,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outvtiname
character outcmlname*200,selectyn
open(ifileid,file=outvtiname,status="replace")
write(ifileid,"(a)") '<?xml version="1.0"?>'
write(ifileid,"(a)") ' <VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
write(ifileid,"(a,6i5,a,3f12.6,a,3f12.6,a)") '  <ImageData WholeExtent="',0,nx-1,0,ny-1,0,nz-1,' " Origin="',&
orgx,orgy,orgz,' " Spacing="',gridv1(1),gridv2(2),gridv3(3),' ">'
write(ifileid,"(a,6i5,a)") '  <Piece Extent="',0,nx-1,0,ny-1,0,nz-1,' ">'
write(ifileid,"(a)") '  <PointData Scalars="scalars">'
write(ifileid,"(a)") '  <DataArray Name="scalars" type="Float64" NumberOfComponents="1" Format="ascii">'
write(ifileid,"(4E14.6)") (((cubmat(i,j,k),i=1,nx),j=1,ny),k=1,nz)
write(ifileid,"(a)") '    </DataArray>'
write(ifileid,"(a)") '    </PointData>'
write(ifileid,"(a)") '    </Piece>'
write(ifileid,"(a)") '    </ImageData>'
write(ifileid,"(a)") ' </VTKFile>'
close(ifileid)
write(*,*) "Exporting vti file finished!"
write(*,*)
write(*,"(a)") " Do you also want to export .cml file in Bohr, so that you can visualize grid data as well as molecular structure in ParaView? (y/n)"
read(*,*) selectyn
if (selectyn=='y') then
    write(*,*) "Input path of .cml file, e.g. C:\aqours.cml"
    read(*,"(a)") outcmlname
    call outcml(outcmlname,ifileid,1)
end if
end subroutine




!!------------- Output grid data in VASP format (like CHGCAR, LOCPOT)    
subroutine outVASPgrd(outname,ifileid)
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) outname

!If cell information is not available, Use grid data vectors as cell vectors
if (all(cellv1==0)) then
    cellv1(:)=gridv1(:)*nx
    cellv2(:)=gridv2(:)*ny
    cellv3(:)=gridv3(:)*nz
    itmp=1
end if

!The first part of CHGCAR/CHG is identical to POSCAR
!Because VASP grid data format does not record origin position, we need to temporarily translate atoms so that origin happens to be (0,0,0)
a(:)%x=a(:)%x-orgx
a(:)%y=a(:)%y-orgy
a(:)%z=a(:)%z-orgz
call outPOSCAR(outname,ifileid)
a(:)%x=a(:)%x+orgx
a(:)%y=a(:)%y+orgy
a(:)%z=a(:)%z+orgz

open(ifileid,file=outname,status="old",position="append")
write(ifileid,"(/,3i6)") nx,ny,nz
write(ifileid,"(5E18.11)") (((cubmat(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
close(ifileid)
write(*,"(a)") " Done, the grid data has been exported in VASP grid data format"

if (itmp==1) then
	cellv1(:)=0
	cellv2(:)=0
	cellv3(:)=0
end if
end subroutine




!!!--------------------- Output current wavefunction to a .wfn
!If isortatmind==1, then any atom without GTF posited on it will not be output, and the index is filled to assure contiguous
!Orbitals with zero occupiation will not be outputted
!If ioutinfo==1, output information
!ifileid: The fileid that can be temporarily in this subroutien
subroutine outwfn(outwfnname,isortatmind,ioutinfo,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outwfnname
integer isortatmind,ioutinfo,ifileid,indconv(ncenter)
!convGmul2wfn converts the g sequence used internally in Multiwfn (input) to commonly used g sequence in .wfn (the one outputted by Molden2AIM and g09 since B01 )
integer :: convGmul2wfn(35)=(/ (0,i=1,20), 23,29,32,27,22, 28,35,34,26,31, 33,30,25,24,21 /)
integer :: convGwfn2mul(35)=(/ (0,i=1,20), 35,25,21,34,33, 29,24,26,22,32, 30,23,31,28,27 /)

open(ifileid,file=outwfnname,status="replace")
write(ifileid,*) "Generated by Multiwfn"

if (isortatmind==1) then !Find real number of centers
	j=0
	do i=1,ncenter
		if (any(b%center==i)) j=j+1
	end do
	write(ifileid,"('GAUSSIAN',i15,' MOL ORBITALS',i7,' PRIMITIVES',i9,' NUCLEI')") count(MOocc(1:nmo)/=0D0),nprims,j
else
	write(ifileid,"('GAUSSIAN',i15,' MOL ORBITALS',i7,' PRIMITIVES',i9,' NUCLEI')") count(MOocc(1:nmo)/=0D0),nprims,ncenter
end if

j=1
do i=1,ncenter
	if (isortatmind==1) then
		if (all(b%center/=i)) cycle
	end if
	indconv(j)=i !The j actual atom corresponds to the i original atom
	write(ifileid,"(2x,a2,i4,4x,'(CENTRE',i3,')',1x,3f12.8,'  CHARGE =',f5.1)") a(i)%name,j,j,a(i)%x,a(i)%y,a(i)%z,a(i)%charge
	j=j+1
end do

if (isortatmind==1) then !Convert center of GTF to reordered index
	do i=1,j-1 !Cycle the centers with GTF
		where(b%center==indconv(i)) b%center=i
	end do
end if

write(ifileid,"('CENTRE ASSIGNMENTS  ',20i3)") b(1:nprims)%center

!Convert the g sequence to common sequence in .wfn, output them, and then convert back. This make output easier
do i=1,nprims
	if (b(i)%type>=21.and.b(i)%type<=35) b(i)%type=convGmul2wfn(b(i)%type)
end do
write(ifileid,"('TYPE ASSIGNMENTS    ',20i3)") b%type
do i=1,nprims
	if (b(i)%type>=21.and.b(i)%type<=35) b(i)%type=convGwfn2mul(b(i)%type)
end do

write(ifileid,"('EXPONENTS ',5D14.7)") b%exp
imo=0
nzeroocc=0
do i=1,nmo
	if (MOocc(i)==0D0) then
		nzeroocc=nzeroocc+1
		cycle
	end if
	imo=imo+1 !Use imo instead of i, this can make MO index contiguous
	write(ifileid,"('MO',I5,'     MO 0.0        OCC NO = ',f12.7,'  ORB. ENERGY =', f12.6)") imo,MOocc(i),MOene(i)
	write(ifileid,"(5D16.8)") (CO(i,j),j=1,nprims)
end do
write(ifileid,"('END DATA',/,' THE  HF ENERGY = ',f19.12,' THE VIRIAL(-V/T)= ',f12.8)") totenergy,virialratio
write(ifileid,"(/,a)") "$MOSPIN"
write(ifileid,"(30i2)") MOtype(:)
write(ifileid,"(a)") "$END"
if (ifPBC>0) then
    write(ifileid,"(/,a)") "[Cell]"
    write(ifileid,"(3f12.6)") cellv1(:)*b2a
    if (ifPBC>=2) write(ifileid,"(3f12.6)") cellv2(:)*b2a
    if (ifPBC==3) write(ifileid,"(3f12.6)") cellv3(:)*b2a
end if
close(ifileid)
if (nzeroocc>0.and.ioutinfo==1) write(*,"(a,i10,a)") " Note: Found",nzeroocc," zero occupied orbitals and have discarded them"
end subroutine



!!!----------------- Output current wavefunction to a .wfx file
!Orbitals with zero occupiation will not be outputted
!If ioutinfo==1, output information
subroutine outwfx(outwfxname,ioutinfo,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outwfxname
character c10tmp*10
integer ioutinfo,ifileid,indconv(ncenter)
! convGseq is used to convert g used in internal sequence of Multiwfn to .wfx sequence
! PS: spdfh sequence in .wfx is identical to Multiwfn
integer :: convGseq(35)
convGseq(21:35)=(/ 23,29,32,27,22,28,35,34,26,31,33,30,25,24,21 /) !g 21~35

nzeroocc=count(MOocc==0D0)
open(ifileid,file=outwfxname,status="replace")
write(ifileid,"(a)") "<Title>"
write(ifileid,*) "Generated by Multiwfn"
write(ifileid,"(a)") "</Title>"
write(ifileid,"(a)") "<Keywords>"
write(ifileid,"(a)") " GTO"
write(ifileid,"(a)") "</Keywords>"
write(ifileid,"(a)") "<Number of Nuclei>"
write(ifileid,"(i6)") ncenter
write(ifileid,"(a)") "</Number of Nuclei>"
write(ifileid,"(a)") "<Number of Occupied Molecular Orbitals>"
write(ifileid,"(i6)") nmo-nzeroocc
write(ifileid,"(a)") "</Number of Occupied Molecular Orbitals>"
write(ifileid,"(a)") "<Number of Perturbations>"
write(ifileid,"(i6)") 0
write(ifileid,"(a)") "</Number of Perturbations>"
write(ifileid,"(a)") "<Net Charge>"
write(ifileid,"(i6)") nint(sum(a%charge)-nelec)
write(ifileid,"(a)") "</Net Charge>"
write(ifileid,"(a)") "<Number of Electrons>"
write(ifileid,"(i6)") nint(nelec)
write(ifileid,"(a)") "</Number of Electrons>"
write(ifileid,"(a)") "<Number of Alpha Electrons>"
write(ifileid,"(i6)") nint(naelec)
write(ifileid,"(a)") "</Number of Alpha Electrons>"
write(ifileid,"(a)") "<Number of Beta Electrons>"
write(ifileid,"(i6)") nint(nbelec)
write(ifileid,"(a)") "</Number of Beta Electrons>"
write(ifileid,"(a)") "<Electronic Spin Multiplicity>"
write(ifileid,"(i6)") nint(naelec-nbelec)+1
write(ifileid,"(a)") "</Electronic Spin Multiplicity>"
write(ifileid,"(a)") "<Number of Core Electrons>"
write(ifileid,"(i6)") nEDFelec
write(ifileid,"(a)") "</Number of Core Electrons>"
write(ifileid,"(a)") "<Nuclear Names>"
do iatm=1,ncenter
	write(c10tmp,"(i6)") iatm
	write(ifileid,"(1x,a)") trim(a(iatm)%name)//trim(adjustl(c10tmp))
end do
write(ifileid,"(a)") "</Nuclear Names>"
write(ifileid,"(a)") "<Atomic Numbers>"
do iatm=1,ncenter
	write(ifileid,"(i3)") a(iatm)%index
end do
write(ifileid,"(a)") "</Atomic Numbers>"
write(ifileid,"(a)") "<Nuclear Charges>"
do iatm=1,ncenter
	write(ifileid,"(E20.12)") a(iatm)%charge
end do
write(ifileid,"(a)") "</Nuclear Charges>"
write(ifileid,"(a)") "<Nuclear Cartesian Coordinates>"
do iatm=1,ncenter
	write(ifileid,"(3E20.12)") a(iatm)%x,a(iatm)%y,a(iatm)%z
end do
write(ifileid,"(a)") "</Nuclear Cartesian Coordinates>"
write(ifileid,"(a)") "<Number of Primitives>"
write(ifileid,"(i6)") nprims
write(ifileid,"(a)") "</Number of Primitives>"
write(ifileid,"(a)") "<Primitive Centers>"
write(ifileid,"(5i20)") b%center
write(ifileid,"(a)") "</Primitive Centers>"
write(ifileid,"(a)") "<Primitive Types>"
icount=0
do iprim=1,nprims
	if (b(iprim)%type>=21.and.b(iprim)%type<=35) then
		write(ifileid,"(i20)",advance="no") convGseq(b(iprim)%type)
	else
		write(ifileid,"(i20)",advance="no") b(iprim)%type
	end if
	icount=icount+1
	if (icount==5.or.iprim==nprims) then
		write(ifileid,*)
		icount=0
	end if
end do
write(ifileid,"(a)") "</Primitive Types>"
write(ifileid,"(a)") "<Primitive Exponents>"
write(ifileid,"(5E20.12)") b%exp
write(ifileid,"(a)") "</Primitive Exponents>"
if (nEDFprims/=0) then
	write(ifileid,"(a)") "<Additional Electron Density Function (EDF)>"
	write(ifileid,"(a)") "<Number of EDF Primitives>"
	write(ifileid,"(i6)") nEDFprims
	write(ifileid,"(a)") "</Number of EDF Primitives>"
	write(ifileid,"(a)") "<EDF Primitive Centers>"
	write(ifileid,"(5i20)") b_EDF%center
	write(ifileid,"(a)") "</EDF Primitive Centers>"
	write(ifileid,"(a)") "<EDF Primitive Types>"
	write(ifileid,"(5i20)") b_EDF%type
	write(ifileid,"(a)") "</EDF Primitive Types>"
	write(ifileid,"(a)") "<EDF Primitive Exponents>"
	write(ifileid,"(5E20.12)") b_EDF%exp
	write(ifileid,"(a)") "</EDF Primitive Exponents>"
	write(ifileid,"(a)") "<EDF Primitive Coefficients>"
	write(ifileid,"(5E20.12)") CO_EDF
	write(ifileid,"(a)") "</EDF Primitive Coefficients>"
	write(ifileid,"(a)") "</Additional Electron Density Function (EDF)>"
end if

write(ifileid,"(a)") "<Molecular Orbital Occupation Numbers>"
do imo=1,nmo
	if (MOocc(imo)/=0D0) write(ifileid,"(E20.12)") MOocc(imo)
end do
write(ifileid,"(a)") "</Molecular Orbital Occupation Numbers>"

write(ifileid,"(a)") "<Molecular Orbital Energies>"
do imo=1,nmo
	if (MOocc(imo)/=0D0) write(ifileid,"(E20.12)") MOene(imo)
end do
write(ifileid,"(a)") "</Molecular Orbital Energies>"

write(ifileid,"(a)") "<Molecular Orbital Spin Types>"
do imo=1,nmo
	if (MOocc(imo)/=0D0) then
		if (MOtype(imo)==0) write(ifileid,"(a)") " Alpha and Beta"
		if (MOtype(imo)==1) write(ifileid,"(a)") " Alpha"
		if (MOtype(imo)==2) write(ifileid,"(a)") " Beta"
	end if
end do
write(ifileid,"(a)") "</Molecular Orbital Spin Types>"

write(ifileid,"(a)") "<Molecular Orbital Primitive Coefficients>"
i=0
nzeroocc=0
do imo=1,nmo
	if (MOocc(imo)==0D0) cycle
	i=i+1 !Use i instead of imo, this can make MO index contiguous
	write(ifileid,"(a)") "<MO Number>"
	write(ifileid,"(i6)") i
	write(ifileid,"(a)") "</MO Number>"
	write(ifileid,"(4E20.12)") CO(imo,:)
end do
write(ifileid,"(a)") "</Molecular Orbital Primitive Coefficients>"

write(ifileid,"(a)") "<Energy = T + Vne + Vee + Vnn>"
write(ifileid,"(E20.12)") totenergy
write(ifileid,"(a)") "</Energy = T + Vne + Vee + Vnn>"

write(ifileid,"(a)") "<Virial Ratio (-V/T)>"
write(ifileid,"(E20.12)") virialratio
write(ifileid,"(a)") "</Virial Ratio (-V/T)>"

close(ifileid)
if (nzeroocc>0.and.ioutinfo==1) write(*,"(a,i10,a)") " Note: Found",nzeroocc," zero occupied orbitals and have discarded them"
end subroutine



!!--------------- Interface of outputting .molden file
subroutine outmolden_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for generating .molden file, e.g. C:\sunshine\riko.molden"
write(*,"(a)") " If press ENTER button directly, will be exported to "//trim(c200tmp)//".molden"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".molden"
write(*,*) "Exporting, please wait..."
call outmolden(outname,10)
write(*,*) "Exporting Molden input file finished!"
end subroutine
!!!------------------------- Output current wavefunction to Molden input file
subroutine outmolden(outname,ifileid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
integer ifileid
character symbol,nowdate*20,nowtime*20

open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "[Molden Format]"
if (ifPBC>0) then
    write(ifileid,"(a)") "[Cell]"
    write(ifileid,"(3f12.6)") cellv1(:)*b2a
    if (ifPBC>=2) write(ifileid,"(3f12.6)") cellv2(:)*b2a
    if (ifPBC==3) write(ifileid,"(3f12.6)") cellv3(:)*b2a
end if
write(ifileid,"(a)") "[Title]"
call date_and_time(nowdate,nowtime)
write(ifileid,"(a,'  Date: ',a,'-',a,'-',a,' Time: ',a,':',a,':',a)") "Generated by Multiwfn",nowdate(1:4),nowdate(5:6),nowdate(7:8),nowtime(1:2),nowtime(3:4),nowtime(5:6)
write(ifileid,"(a)") "[Atoms] AU"
do i=1,ncenter
	write(ifileid,"(a,i7,i4,3f14.7)") a(i)%name,i,a(i)%index,a(i)%x,a(i)%y,a(i)%z
end do
write(ifileid,"(a)") "[GTO]"
do iatm=1,ncenter
	write(ifileid,"(2i6)") iatm,0
	do ish=1,nshell
		if (shcen(ish)==iatm) then
			symbol=shtype2name(shtype(ish))
			call struc2lc(symbol)
			write(ifileid,"(a,i4,' 1.0')") symbol,shcon(ish)
			if (ish==0) then
				istart=0
			else
				istart=sum(shcon(1:ish-1))
			end if
			do ipsh=istart+1,istart+shcon(ish)
				write(ifileid,"(2(1PE16.8))") primshexp(ipsh),primshcoeff(ipsh)
			end do
		end if
	end do
	write(ifileid,*)
end do
write(ifileid,*) 

if (any(shtype==-2).and.any(shtype==-3)) then !Default is 6d10f
	write(ifileid,"('[5D]')") !5d7f
else if (any(shtype==-2)) then
	write(ifileid,"('[5D10F]')") !5d10f
else if (any(shtype==-3)) then
	write(ifileid,"('[7F]')") !=6d7f
end if
if (any(shtype==-4)) write(ifileid,"('[9G]')") !Default is Cartesian G
if (any(shtype==-5)) write(ifileid,"('[11H]')") !Default is Cartesian H
write(ifileid,"(a)") "[MO]"
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Close shell
	do imo=1,nmo
		write(ifileid,"('Ene=',f16.8)") MOene(imo)
		write(ifileid,"('Spin= Alpha')")
		write(ifileid,"('Occup=',f10.6)") MOocc(imo)
		do ibas=1,nbasis
			write(ifileid,"(i6,f18.10)") ibas,CObasa(ibas,imo)
		end do
	end do
else !Open shell
	do isep=nmo,1,-1
		if (MOtype(isep)==1) exit
	end do
	do imo=1,isep !Alpha part
		write(ifileid,"('Ene=',f16.8)") MOene(imo)
		write(ifileid,"('Spin= Alpha')")
		write(ifileid,"('Occup=',f10.6)") MOocc(imo)
		do ibas=1,nbasis
			write(ifileid,"(i6,f18.10)") ibas,CObasa(ibas,imo)
		end do
	end do
	do imo=isep+1,nmo !Beta part
		write(ifileid,"('Ene=',f16.8)") MOene(imo)
		write(ifileid,"('Spin= Beta')")
		write(ifileid,"('Occup=',f10.6)") MOocc(imo)
		do ibas=1,nbasis
			write(ifileid,"(i6,f18.10)") ibas,CObasb(ibas,imo-isep)
		end do
	end do
end if

close(ifileid)
end subroutine




!!!--------------- Output current wavefunction to mkl file (old Molekel input file), then orca_2mkl can convert it to .gbw
!The format is exactly identical to the .mkl file produced by orca_2mkl (e.g. orca_2mkl test can generate test.mkl from test.gbw)
!Note that when if the gbw file converted from .mkl is used as initial guess, the number of orbitals recorded in the .mkl can be &
!smaller than number of basis functions. The only important point is that all occupied orbitals are recorded. The virtual orbitals &
!can be ignored in the .mkl, the coefficients can also be all zero. Therefore, when exporting .mkl file, it is not necessary to check &
!if the final several orbitals have all zero coeffienents due to the loaded wavefunction file has basis linear dependency
subroutine outmkl(outname,ifileid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
integer ifileid
character symbol,writeformat*20,selectyn
character(len=4) irrep(nmo)

write(*,*) "Will the generated .mkl file be used for ORCA? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.and.any(bastype(1:nbasis)>=5)) then
	write(*,*) "Error: Some basis functions are Cartesian type, however, ORCA only supports spherical-harmonic basis functions"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if
write(*,*) "Exporting, please wait..."

open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "$MKL"
write(ifileid,"(a)") "#"
write(ifileid,"(a)") "# MKL format file produced by Multiwfn"
write(ifileid,"(a)") "#"
write(ifileid,"(a)") "$CHAR_MULT"
write(ifileid,"(2i3)") nint(sum(a(:)%charge)-naelec-nbelec),nint(naelec-nbelec+1)
write(ifileid,"(a)") "$END"
write(ifileid,*)

write(ifileid,"(a)") "$COORD"
do i=1,ncenter
	write(ifileid,"(i4,1x,3f11.6)") a(i)%index,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
end do
write(ifileid,"(a)") "$END"
write(ifileid,*)

write(ifileid,"(a)") "$BASIS"
do iatm=1,ncenter
	do ish=1,nshell
		if (shcen(ish)==iatm) then
			symbol=shtype2name(shtype(ish))
			write(ifileid,"(i2,1x,a,f4.1)") shtype2nbas(shtype(ish)),symbol,1D0
			if (ish==0) then
				istart=0
			else
				istart=sum(shcon(1:ish-1))
			end if
			do ipsh=istart+1,istart+shcon(ish)
				write(ifileid,"(f19.9,f17.9)") primshexp(ipsh),primshcoeff(ipsh)
			end do
		end if
	end do
	if (iatm/=ncenter) write(ifileid,"(a)") "$$"
end do
write(ifileid,"(a)") "$END"
write(ifileid,*)

nframe=ceiling(nbasis/5D0)
!if (allocated(MOsym)) then
!    irrep=MOsym
!else
    irrep=" a1g" !orca_2mkl convention
!end if
write(ifileid,"(a)") "$COEFF_ALPHA"
do iframe=1,nframe
    if (iframe/=nframe) then
        ncol=5
    else
        ncol=nbasis-(nframe-1)*5
    end if
    ibeg=(iframe-1)*5+1
    iend=(iframe-1)*5+ncol
    write(writeformat,"('(',i1,'(a,1x))')") ncol
    write(ifileid,writeformat) irrep(ibeg:iend)
    write(writeformat,"('(',i1,'(f13.7,1x))')") ncol
    write(ifileid,writeformat) MOene(ibeg:iend)
    write(writeformat,"('(',i1,'(f12.7,1x))')") ncol
	do ibas=1,nbasis
        itmp=bastype(ibas)
        if ((selectyn=='y'.or.selectyn=='Y').and.(itmp==-6.or.itmp==-7.or.(itmp>=-16.and.itmp<=-13).or.(itmp>=-29.and.itmp<=-26))) then
            !write(*,"(2i5,1x,a)") ibas,itmp,GTFtype2name(itmp)
            write(ifileid,writeformat) -CObasa(ibas,ibeg:iend)
        else
            write(ifileid,writeformat) CObasa(ibas,ibeg:iend)
        end if
	end do
end do
write(ifileid,"(a)") " $END"
write(ifileid,*)

write(ifileid,"(a)") "$OCC_ALPHA"
write(ifileid,"(5f12.7)") MOocc(1:nbasis)
write(ifileid,"(a)") " $END"
    
if (wfntype==1.or.wfntype==4) then !Open shell
    write(ifileid,*)
    write(ifileid,"(a)") "$COEFF_BETA"
    do iframe=1,nframe
        if (iframe/=nframe) then
            ncol=5
        else
            ncol=nbasis-(nframe-1)*5
        end if
        ibeg=(iframe-1)*5+1
        iend=(iframe-1)*5+ncol
        write(writeformat,"('(',i1,'(a,1x))')") ncol
        write(ifileid,writeformat) irrep(ibeg+nbasis:iend+nbasis)
        write(writeformat,"('(',i1,'(f13.7,1x))')") ncol
        write(ifileid,writeformat) MOene(ibeg+nbasis:iend+nbasis)
        write(writeformat,"('(',i1,'(f12.7,1x))')") ncol
	    do ibas=1,nbasis
            itmp=bastype(ibas)
            if ((selectyn=='y'.or.selectyn=='Y').and.(itmp==-6.or.itmp==-7.or.(itmp>=-16.and.itmp<=-13).or.(itmp>=-29.and.itmp<=-26))) then
                write(ifileid,writeformat) -CObasb(ibas,ibeg:iend)
            else
                write(ifileid,writeformat) CObasb(ibas,ibeg:iend)
            end if
	    end do
    end do
    write(ifileid,"(a)") " $END"
    write(ifileid,*)

    write(ifileid,"(a)") "$OCC_BETA"
    write(ifileid,"(5f12.7)") MOocc(nbasis+1:nmo)
    write(ifileid,"(a)") " $END"
end if

close(ifileid)
write(*,*) "Exporting .mkl file finished!"
end subroutine



!!--------------- Interface of outputting .fch file
subroutine outfch_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for outputting .fch file, e.g. C:\sunshine\riko.fch"
write(*,"(a)") " If press ENTER button directly, the system will be exported to "//trim(c200tmp)//".fch in current folder"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".fch"
call outfch(outname,10,1)
end subroutine
!!!--------------- Output current wavefunction to .fch file
!informode=1: print prompt message, =0: do not print
subroutine outfch(outname,ifileid,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
integer ifileid,infomode
open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "Generated by Multiwfn"
if (wfntype==0.or.wfntype==3) write(ifileid,"(a10,a30,a30)") "SP        ","RB3LYP                        ","                      6-31G(d)"
if (wfntype==1.or.wfntype==4) write(ifileid,"(a10,a30,a30)") "SP        ","UB3LYP                        ","                      6-31G(d)"
if (wfntype==2)               write(ifileid,"(a10,a30,a30)") "SP        ","ROB3LYP                       ","                      6-31G(d)"
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of atoms                         ","I",ncenter
write(ifileid,"(A40,3X,A1,5X,I12)") "Charge                                  ","I",nint(sum(a%charge)-nelec)
write(ifileid,"(A40,3X,A1,5X,I12)") "Multiplicity                            ","I",nint(naelec-nbelec)+1
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of electrons                     ","I",nint(nelec)
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of alpha electrons               ","I",nint(naelec)
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of beta electrons                ","I",nint(nbelec)
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of basis functions               ","I",nbasis
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of independent functions         ","I",nbasis
if (ifPBC>0) write(ifileid,"(A40,3X,A1,5X,I12)") "Number of translation vectors           ","I",ifPBC
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Atomic numbers                          ","I",ncenter
write(ifileid,"(6I12)") a(:)%index
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Nuclear charges                         ","R",ncenter
write(ifileid,"(5(1PE16.8))") a(:)%charge
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Current cartesian coordinates           ","R",ncenter*3
write(ifileid,"(5(1PE16.8))") (a(i)%x,a(i)%y,a(i)%z,i=1,ncenter)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Translation vectors                     ","R",ifPBC*3
if (ifPBC==1) write(ifileid,"(3(1PE16.8))") cellv1(:)
if (ifPBC==2) write(ifileid,"(5(1PE16.8))") cellv1(:),cellv2(:)
if (ifPBC==3) write(ifileid,"(5(1PE16.8))") cellv1(:),cellv2(:),cellv3(:)
!Basis function definition
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of contracted shells             ","I",nshell
write(ifileid,"(A40,3X,A1,5X,I12)") "Number of primitive shells              ","I",nprimshell
write(ifileid,"(A40,3X,A1,5X,I12)") "Highest angular momentum                ","I",maxval(abs(shtype))
write(ifileid,"(A40,3X,A1,5X,I12)") "Largest degree of contraction           ","I",maxval(abs(shcon))
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Shell types                             ","I",nshell
write(ifileid,"(6I12)") shtype(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Number of primitives per shell          ","I",nshell
write(ifileid,"(6I12)") shcon(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Shell to atom map                       ","I",nshell
write(ifileid,"(6I12)") shcen(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Primitive exponents                     ","R",nprimshell
write(ifileid,"(5(1PE16.8))") primshexp(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Contraction coefficients                ","R",nprimshell
write(ifileid,"(5(1PE16.8))") primshcoeff(:)
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Coordinates of each shell               ","R",nshell*3
write(ifileid,"(5(1PE16.8))") (a(shcen(i))%x,a(shcen(i))%y,a(shcen(i))%z,i=1,nshell)
write(ifileid,"(A40,3X,A1,5X,1PE22.15)") "Virial Ratio                            ","R",virialratio
write(ifileid,"(A40,3X,A1,5X,1PE22.15)") "Total Energy                            ","R",totenergy
!After careful investigation, I found the following information is needed if you want to use unfchk &
!to convert it to chk and read guess from it, without it Gaussian will only read alpha orbitals as guess.
!According to Gaussian programmer's manual, the meaning of the first slot of ILSW is: 0=real RHF (default) 1=real UHF, etc.
!R and RO calculations do not need this
if (wfntype==1.or.wfntype==4) then
    write(ifileid,"(a)") "Num ILSW                                   I              1"
    write(ifileid,"(a)") "ILSW                                       I   N=         1"
    write(ifileid,"(a)") "           1"
end if
!Orbital informaiton
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Alpha Orbital Energies                  ","R",nbasis
write(ifileid,"(5(1PE16.8))") MOene(1:nbasis)
if (wfntype==1.or.wfntype==4) then
	write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Beta Orbital Energies                   ","R",nbasis
	write(ifileid,"(5(1PE16.8))") MOene(nbasis+1:2*nbasis)
end if
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Alpha MO coefficients                   ","R",nbasis*nbasis
write(ifileid,"(5(1PE16.8))") ((CObasa(ibasis,imo),ibasis=1,nbasis),imo=1,nbasis)
if (wfntype==1.or.wfntype==4) then
	write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Beta MO coefficients                    ","R",nbasis*nbasis
	write(ifileid,"(5(1PE16.8))") ((CObasb(ibasis,imo),ibasis=1,nbasis),imo=1,nbasis)
end if
write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Total SCF Density                       ","R",nbasis*(nbasis+1)/2
write(ifileid,"(5(1PE16.8))") ((Ptot(i,j),j=1,i),i=1,nbasis)
if (wfntype==1.or.wfntype==4) then
	write(ifileid,"(A40,3X,A1,3X,'N=',I12)") "Spin SCF Density                        ","R",nbasis*(nbasis+1)/2
	write(ifileid,"(5(1PE16.8))") ((Palpha(i,j)-Pbeta(i,j),j=1,i),i=1,nbasis)
end if
close(ifileid)
if (infomode==1) then
	if (wfntype==0.or.wfntype==2.or.wfntype==3) then
		write(*,"(a)") " Exporting .fch file finished! The ""Total SCF Density"" field in this file corresponds to density matrix at current level"
	else
		write(*,"(a)") " Exporting .fch file finished! The ""Total SCF Density"" and ""Spin SCF Density"" fields in this file correspond &
		&to total density matrix and difference between alpha and beta density matrices at current level, respectively"
	end if
end if
end subroutine


!!!---------- Output current wavefunction to .47 file, which is input file of NBO program
!Notice that although "UPPER" keyword is used here, in fact the matrix recorded in .47 is lower-triangular matrix (i.e. NBO's rule is confusing!)
subroutine out47(outname,ifileid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
character(80) c80tmp1,c80tmp2
real*8 primshcoefftmp(nprimshell)
real*8,allocatable :: halfmat(:)
real*8,external :: normgau
integer bastype2NBO(-12:50),nptr(nshell)
integer ifileid
bastype2NBO(-12)=351 !F 0
bastype2NBO(-11)=352 !F+1
bastype2NBO(-10)=353 !F-1
bastype2NBO(-9)=354 !F+2
bastype2NBO(-8)=355 !F-2
bastype2NBO(-7)=356 !F+3
bastype2NBO(-6)=357 !F-3
bastype2NBO(-5)=255 !z2
bastype2NBO(-4)=252 !xz
bastype2NBO(-3)=253 !yz
bastype2NBO(-2)=254 !x2-y2
bastype2NBO(-1)=251 !xy
bastype2NBO(1 )=1   !s
bastype2NBO(2 )=101 !x
bastype2NBO(3 )=102 !y
bastype2NBO(4 )=103 !z
bastype2NBO(5 )=201 !xx
bastype2NBO(8 )=202 !xy
bastype2NBO(9 )=203 !xz
bastype2NBO(6 )=204 !yy
bastype2NBO(10)=205 !yz
bastype2NBO(7 )=206 !zz
bastype2NBO(11)=301 !xxx
bastype2NBO(14)=302 !xxy
bastype2NBO(15)=303 !xxz
bastype2NBO(17)=304 !xyy
bastype2NBO(20)=305 !xyz
bastype2NBO(18)=306 !xzz
bastype2NBO(12)=307 !yyy
bastype2NBO(16)=308 !yyz
bastype2NBO(19)=309 !yzz
bastype2NBO(13)=310 !zzz
!Below g sequence comes from line 47384 in NBO_5 src
bastype2NBO(35)=401 !XXXX
bastype2NBO(34)=402 !XXXY
bastype2NBO(33)=403 !XXXZ
bastype2NBO(32)=404 !XXYY
bastype2NBO(31)=405 !XXYZ
bastype2NBO(30)=406 !XXZZ
bastype2NBO(29)=407 !XYYY
bastype2NBO(28)=408 !XYYZ
bastype2NBO(27)=409 !XYZZ
bastype2NBO(26)=410 !XZZZ
bastype2NBO(25)=411 !YYYY
bastype2NBO(24)=412 !YYYZ
bastype2NBO(23)=413 !YYZZ
bastype2NBO(22)=414 !YZZZ
bastype2NBO(21)=415 !ZZZZ

if (any(shtype<-3)) then
	write(*,"(a)") " Error: This function cannot be used if any spherical-harmonic basis function has >=g angular moment!"
	write(*,"(a)") " To make this function usable, there are three ways:"
    write(*,"(a)") " (1) Reduce quality of basis set so that highest angular moment does not exceed d. Note that most NBO analyses are quite insensitive to basis set"
    write(*,"(a)") " (2) Perform calculation based on Cartesian instead of spherical-harmonic basis functions"
    write(*,"(a)") " (3) Before booting up Multiwfn, set ""iloadasCart"" in settings.ini to 1. After that all spherical-harmonic &
    basis functions will be converted to Cartesian ones when loading input file. However, note that in this case Multiwfn is unable &
    to correctly generate Fock matrix based on orbital energies and coeffcients"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

write(*,*) "Exporting, please wait..."

open(ifileid,file=outname,status="replace")
write(c80tmp1,*) ncenter
write(c80tmp2,*) nbasis
if (wfntype==0.or.wfntype==3) then !Closed-shell
	write(ifileid,"(' $GENNBO NATOMS=',a,' NBAS=',a,' UPPER BODM $END')") trim(adjustl(c80tmp1)),trim(adjustl(c80tmp2))
else !Open-shell
	write(ifileid,"(' $GENNBO NATOMS=',a,' NBAS=',a,' UPPER BODM OPEN $END')") trim(adjustl(c80tmp1)),trim(adjustl(c80tmp2))
end if
write(ifileid,*) "$NBO $END"
write(ifileid,*) "$COORD"
write(ifileid,*) "Generated by Multiwfn"
do iatm=1,ncenter
	write(ifileid,"(2i6,3f12.6)") a(iatm)%index,int(a(iatm)%charge),a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
end do
write(ifileid,*) "$END"
!Basis function information
write(ifileid,*) "$BASIS"
write(ifileid,"(' CENTER =')")
write(ifileid,"(10i6)") bascen(1:nbasis)
write(ifileid,"(' LABEL =')")
write(ifileid,"(10i6)") bastype2NBO(bastype(1:nbasis))
write(ifileid,"(' $END')")

!Shell information
write(ifileid,*) "$CONTRACT"
write(ifileid,"(' NSHELL =',i6)") nshell
write(ifileid,"('   NEXP =',i6)") nprimshell
write(ifileid,"('  NCOMP =')")
write(ifileid,"(10i6)") shtype2nbas(shtype(:))
write(ifileid,"('  NPRIM =')")
write(ifileid,"(10i6)") shcon(:)
nptr(1)=1
do ish=2,nshell
	nptr(ish)=nptr(ish-1)+shcon(ish-1)
end do
write(ifileid,"('   NPTR =')")
write(ifileid,"(10i6)") nptr(:)
write(ifileid,"('    EXP =')")
write(ifileid,"(4E16.7)") primshexp(:)
!In standard .47 and .31 files, the shell contraction coefficients include normalization coefficients
!For d, f, g, the normalization coefficients are for XX/YY/ZZ, XXX/YYY/ZZZ, XXXX/YYYY/ZZZZ, respectively. &
!I think in fact this is problematic, because normalization coefficient of e.g. XX and XY is different, I don't know how NBO deal with this issue
do iang=0,maxval(abs(shtype))
	primshcoefftmp=0
	if (iang==0) then
		write(ifileid,"('     CS =')")
		do ish=1,nshell
			if (shtype(ish)==0) then
				do icon=1,shcon(ish)
					iprimshnow=nptr(ish)+icon-1
					primshcoefftmp(iprimshnow)=primshcoeff(iprimshnow)*normgau(1,primshexp(iprimshnow))
				end do
			end if
		end do
	else if (iang==1) then
		write(ifileid,"('     CP =')")
		do ish=1,nshell
			if (shtype(ish)==1) then
				do icon=1,shcon(ish)
					iprimshnow=nptr(ish)+icon-1
					primshcoefftmp(iprimshnow)=primshcoeff(iprimshnow)*normgau(2,primshexp(iprimshnow))
				end do
			end if
		end do
	else if (iang==2) then
		write(ifileid,"('     CD =')")
		do ish=1,nshell
			if (abs(shtype(ish))==2) then
				do icon=1,shcon(ish)
					iprimshnow=nptr(ish)+icon-1
					primshcoefftmp(iprimshnow)=primshcoeff(iprimshnow)*normgau(5,primshexp(iprimshnow))
				end do
			end if
		end do
	else if (iang==3) then
		write(ifileid,"('     CF =')")
		do ish=1,nshell
			if (abs(shtype(ish))==3) then
				do icon=1,shcon(ish)
					iprimshnow=nptr(ish)+icon-1
					primshcoefftmp(iprimshnow)=primshcoeff(iprimshnow)*normgau(11,primshexp(iprimshnow))
					!do itmp=11,20
					!	write(*,*) itmp,primshcoeff(iprimshnow)*normgau(itmp,primshexp(iprimshnow))
					!end do
				end do
			end if
		end do
	else if (iang==4) then
		write(ifileid,"('     CG =')")
		do ish=1,nshell
			if (abs(shtype(ish))==4) then
				do icon=1,shcon(ish)
					iprimshnow=nptr(ish)+icon-1
					primshcoefftmp(iprimshnow)=primshcoeff(iprimshnow)*normgau(21,primshexp(iprimshnow))
				end do
			end if
		end do
	end if
	write(ifileid,"(4E16.7)") primshcoefftmp(:)
end do
write(ifileid,"(' $END')")

allocate(halfmat(nbasis*(nbasis+1)/2))

!Overlap matrix
write(ifileid,"(' $OVERLAP')")
call mat2arr(Sbas,halfmat,2)
write(ifileid,"(5E15.7)") halfmat(:)
write(ifileid,"(' $END')")

!Density matrix
write(ifileid,"(' $DENSITY')")
if (wfntype==0.or.wfntype==3) then !Closed-shell
	call mat2arr(Ptot,halfmat,2)
	write(ifileid,"(5E15.7)") halfmat(:)
else !Open-shell
	call mat2arr(Palpha,halfmat,2)
	write(ifileid,"(5E15.7)") halfmat(:)
	call mat2arr(Pbeta,halfmat,2)
	write(ifileid,"(5E15.7)") halfmat(:)
end if
write(ifileid,"(' $END')")

!Fock matrix
if (allocated(FmatA).or.wfntype<=2) then
    write(*,*)
    write(*,*) "If writing Fock matrix to .47 file?"
    write(*,*) "0 Do not write"
    if (wfntype<=2) write(*,"(a)") " 1 Generate Fock matrix based on MO energies and coefficients via F=SCE(C^-1) and write it to .47 file"
    if (allocated(FmatA)) write(*,"(a)") " 2 Write the Fock matrix in memory to .47 file"
    read(*,*) iwriteF
    if (iwriteF==1) then
        call MOene2Fmat(istatus)
        if (istatus==1) write(*,"(a)") " Warning: Since Fock matrix was not successfully generated and thus will not be written to .47 file"
    end if
    if ((iwriteF==1.and.istatus==0).or.iwriteF==2) then
        write(ifileid,"(' $FOCK')")
	    call mat2arr(FmatA,halfmat,2)
	    write(ifileid,"(5E15.7)") halfmat(:)
        if (wfntype==1) then
	        call mat2arr(FmatB,halfmat,2)
	        write(ifileid,"(5E15.7)") halfmat(:)
        end if
        write(ifileid,"(' $END')")
    end if
end if

!LCAOMO matrix. Note that if "iloadasCart" is set to 1, when loading .molden and .fch where spherical harmonic basis functions are presented, &
!they will be converted to Cartesian type and retain this status. In this case some highest MOs have zero coefficients, but NBO can still work normally
! if (nmo==nbasis.or.nmo==2*nbasis) then !The input file must only contain Cartesian basis functions
write(ifileid,"(' $LCAOMO')")
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !R or RO
	do imo=1,nmo
		write(ifileid,"(5E15.7)") CObasa(:,imo)
	end do
	if (wfntype==2) then !RO, output orbitals twice
		do imo=1,nmo
			write(ifileid,"(5E15.7)") CObasa(:,imo)
		end do
	end if
else !U
	do imo=1,nbasis
		write(ifileid,"(5E15.7)") CObasa(:,imo)
	end do
	do imo=1,nbasis
		write(ifileid,"(5E15.7)") CObasb(:,imo)
	end do
end if
write(ifileid,"(' $END')")

!Dipole matrix
if (.not.allocated(Dbas)) then
    write(*,*) "Generating electric dipole moment integral matrix..."
    call genDbas_curr
end if

write(ifileid,"(' $DIPOLE')")
call mat2arr(Dbas(1,:,:),halfmat,2)
write(ifileid,"(5E15.7)") halfmat(:)*(-b2a) !Must be converted from Bohr to Angstrom
call mat2arr(Dbas(2,:,:),halfmat,2)
write(ifileid,"(5E15.7)") halfmat(:)*(-b2a)
call mat2arr(Dbas(3,:,:),halfmat,2)
write(ifileid,"(5E15.7)") halfmat(:)*(-b2a)
write(ifileid,"(' $END')")

close(ifileid)
write(*,*)
write(*,*) "Exporting .47 file finished!"
end subroutine



!!!------------ Interface of outputting .mwfn file
subroutine outmwfn_wrapper
use defvar
use util
character(len=200) outname,c200tmp
integer iexp
if (.not.allocated(CObasa)) then
	write(*,"(a)") " Error: This function works only when input file contains basis function information"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
else
    call path2filename(filename,c200tmp)
	write(*,*) "Input path for exporting file, e.g. C:\ltwd.mwfn"
    write(*,"(a)") " If press ENTER button directly, will be exported to "//trim(c200tmp)//".mwfn"
	read(*,"(a)") outname
    if (outname==" ") outname=trim(c200tmp)//".mwfn"
    write(*,*) "0 Only export wavefunction"
    write(*,*) "1 Export both wavefunction and density matrix"
    write(*,*) "2 Export wavefunction, density matrix and overlap matrix"
    read(*,*) iexp
	write(*,*) "Exporting, please wait..."
	call outmwfn(outname,10,iexp)
	write(*,*) "Exporting .mwfn file finished!"
end if
end subroutine
!!!------------ Output current wavefunction to .mwfn file
!ifileid is the file id that can be used during outputting
!iprint=0: Do not print any matrix, =1: Print density matrix, =2: Print density matrix and overlap matrix
!informode=1: Output prompt message, =0: Do not Output
subroutine outmwfn(outname,ifileid,iprint)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
integer ifileid,infomode
open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "# Generated by Multiwfn" !Comment
write(ifileid,"('Wfntype=',i4)") wfntype
write(ifileid,"('Charge=',f15.6)") sum(a%charge)-(naelec+nbelec)
write(ifileid,"('Naelec=',f15.6)") naelec
write(ifileid,"('Nbelec=',f15.6)") nbelec
write(ifileid,"('E_tot=',1PE16.8)") totenergy
write(ifileid,"('VT_ratio=',f12.8)") virialratio
if (ifPBC>0) then
    write(ifileid,"('Ndim=',i4)") ifPBC
    write(ifileid,"('cellv1=',3f12.8)") cellv1(:)*b2a
    if (ifPBC>1) write(ifileid,"('cellv2=',3f12.8)") cellv2(:)*b2a
    if (ifPBC>2) write(ifileid,"('cellv3=',3f12.8)") cellv3(:)*b2a
end if
write(ifileid,"(/,a)") "# Atom information"
write(ifileid,"('Ncenter=',i8)") ncenter
write(ifileid,"('$Centers')")
do icen=1,ncenter
    write(ifileid,"(i6,1x,a,i4,f6.1,3f16.8)") icen,a(icen)%name,a(icen)%index,a(icen)%charge,a(icen)%x*b2a,a(icen)%y*b2a,a(icen)%z*b2a
end do

write(ifileid,"(/,a)") "# Basis function information"
write(ifileid,"('Nbasis=    ',i8)") nbasis
if (nindbasis==0) then !nindbasis was not loaded, simply set the same as nbasis
    write(ifileid,"('Nindbasis= ',i8)") nbasis
else
    write(ifileid,"('Nindbasis= ',i8)") nindbasis
end if
write(ifileid,"('Nprims=    ',i8)") nprims
write(ifileid,"('Nshell=    ',i8)") nshell
write(ifileid,"('Nprimshell=',i8)") nprimshell
write(ifileid,"('$Shell types')")
write(ifileid,"(25i3)") shtype(:)
write(ifileid,"('$Shell centers')")
write(ifileid,"(10i8)") shcen(:)
write(ifileid,"('$Shell contraction degrees')")
write(ifileid,"(20i4)") shcon(:)
write(ifileid,"('$Primitive exponents')")
write(ifileid,"(5(1PE16.8))") primshexp(:)
write(ifileid,"('$Contraction coefficients')")
write(ifileid,"(5(1PE16.8))") primshcoeff(:)

!Output orbital information
if (wfntype==0.or.wfntype==3) then
    write(ifileid,"(/,a)") "# Orbital information (nindbasis orbitals)"
else
    write(ifileid,"(/,a)") "# Orbital information (2*nindbasis orbitals)"
end if
irealorb=0
nzero=0
do imo=1,nmo
    !Orbitals having all-zero coefficents, which is caused by linear-dependency, are not exported
    if (MOtype(imo)<=1) then !Spatial or alpha orbital
        if (all(CObasa(:,imo)==0).and.nindbasis<nbasis) then
            nzero=nzero+1
            cycle
        end if
    else
        if (all(CObasb(:,imo-nbasis)==0).and.nindbasis<nbasis) then
            nzero=nzero+1
            cycle
        end if
    end if
    
    irealorb=irealorb+1
    write(ifileid,*)
    write(ifileid,"('Index=',i10)") irealorb
	write(ifileid,"('Type=',i2)") MOtype(imo)
	write(ifileid,"('Energy=',1PE16.8)") MOene(imo)
	write(ifileid,"('Occ=',f12.8)") MOocc(imo)
	if (allocated(MOsym)) then
        if (MOsym(imo)/=" ") then
            write(ifileid,"('Sym=',1x,a)") MOsym(imo)
        else
            write(ifileid,"('Sym= ?')")
        end if
    else
        write(ifileid,"('Sym= ?')")
    end if
    write(ifileid,"('$Coeff')")
    if (MOtype(imo)<=1) then !Spatial or alpha orbital
		write(ifileid,"(5(1PE16.8))") CObasa(:,imo)
    else
		write(ifileid,"(5(1PE16.8))") CObasb(:,imo-nbasis)
    end if
end do
if (nzero>0) write(*,"(a,i5,a)") " Note: There are",nzero," orbitals with all-zero coefficients, which were not exported"

if (iprint>=1) then
    write(ifileid,"(/,a)") "# Various matrices"
    if (wfntype==0.or.wfntype==3) then
        write(ifileid,"(/,'$Total density matrix, dim=',2i5,'  lower= 1')") nbasis,nbasis !nbasis*(nbasis+1)/2 elements
        write(ifileid,"(5(1PE16.8))") ((Ptot(i,j),j=1,i),i=1,nbasis)
    else
        write(ifileid,"(/,'$Alpha density matrix, dim=',2i5,'  lower= 1')") nbasis,nbasis !nbasis*(nbasis+1)/2 elements
        write(ifileid,"(5(1PE16.8))") ((Palpha(i,j),j=1,i),i=1,nbasis)
        write(ifileid,"(/,'$Beta density matrix, dim=',2i5,'  lower= 1')") nbasis,nbasis !nbasis*(nbasis+1)/2 elements
        write(ifileid,"(5(1PE16.8))") ((Pbeta(i,j),j=1,i),i=1,nbasis)
    end if
    if (iprint==2) then
		call ask_sbas_PBC
        write(ifileid,"(/,'$Overlap matrix, dim=',2i5,'  lower= 1')") nbasis,nbasis !nbasis*(nbasis+1)/2 elements
        write(ifileid,"(5(1PE16.8))") ((Sbas(i,j),j=1,i),i=1,nbasis)
    end if
end if
close(ifileid)
end subroutine




!!!------------ Read wavefunction from .mwfn file
!infomode=0 means output info, =1 silent
subroutine readmwfn(name,infomode)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) name
!Temporary arrays used in process data
integer,allocatable :: shelltype(:),shell2atom(:),shellcon(:)
real*8,allocatable :: primexp(:),concoeff(:),amocoeff(:,:),bmocoeff(:,:)
integer :: s2f(-5:5,21)=0
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
character selectyn,c80tmp*80,c80tmp2*80
!For backing up spherical basis functions
integer,allocatable :: shelltype5D(:),MOtype5D(:)
real*8,allocatable :: CObasa5D(:,:),CObasb5D(:,:),MOocc5D(:),MOene5D(:),CO5D(:,:)
real*8,external :: normgau

ifiletype=14
imodwfn=0
s2f(-5,1:11)=(/ -32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22 /)
s2f(-4,1:9)=(/ -21,-20,-19,-18,-17,-16,-15,-14,-13 /)
s2f(-3,1:7)=(/ -12,-11,-10,-9,-8,-7,-6 /)
s2f(-2,1:5)=(/ -5,-4,-3,-2,-1 /)
s2f(-1,1:4)=(/ 1,2,3,4 /)
s2f(0,1)=1
s2f(1,1:3)=(/ 2,3,4 /)
s2f(2,1:6)=(/ 5,6,7,8,9,10 /)
s2f(3,1:10)=(/ 11,12,13,17,14,15,18,19,16,20 /) !Note: The sequence of f functions in Multiwfn is not identical to .fch, so convert here. While spdgh are identical
s2f(4,1:15)=(/ 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35 /)
s2f(5,1:21)=(/ 36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56 /)
call gensphcartab(1,conv5d6d,conv7f10f,conv9g15g,conv11h21h)

open(10,file=name,status="old")
if (infomode==0) write(*,*) "Loading various information of the wavefunction"

Nframe=1
call loclabel(10,'@Nframe',ifound,maxline=10)
if (ifound==1) then
    read(10,*) c80tmp,Nframe
    rewind(10)
    write(*,"(/,' There are',i8,' frames, load which one? e.g. 5')") Nframe
    read(*,*) iframe
end if

call loclabel(10,'Wfntype=')
read(10,*) c80tmp,wfntype
!call loclabel(10,'Charge=') !This information is not explicitly used by Multiwfn
!read(10,*) c80tmp,Charge
call loclabel(10,'Naelec=')
read(10,*) c80tmp,naelec
call loclabel(10,'Nbelec=')
read(10,*) c80tmp,nbelec
call loclabel(10,'E_tot=',ifound,maxline=100)
totenergy=0
if (ifound==1) read(10,*) c80tmp,totenergy
call loclabel(10,'VT_ratio=',ifound,maxline=100)
virialratio=0
if (ifound==1) read(10,*) c80tmp,virialratio
call loclabel(10,'Ndim',ifound,maxline=100)
ifPBC=0
if (ifound==1) read(10,*) c80tmp,ifPBC
if (ifPBC>0) then
    call loclabel(10,'cellv1',ifound,0)
    read(10,*) c80tmp,cellv1
end if
if (ifPBC>1) then
    call loclabel(10,'cellv2',ifound,0)
    read(10,*) c80tmp,cellv2
end if
if (ifPBC>2) then
    call loclabel(10,'cellv3',ifound,0)
    read(10,*) c80tmp,cellv3
end if

!Load atom information
call loclabel(10,'Ncenter=')
read(10,*) c80tmp,ncenter
allocate(a(ncenter))
call loclabel(10,'$Centers');read(10,*)
do icen=1,ncenter
    read(10,*) inouse,a(icen)%name,a(icen)%index,a(icen)%charge,a(icen)%x,a(icen)%y,a(icen)%z
end do

!Load basis function information
if (infomode==0) write(*,*) "Loading basis and primitive shell information..."
call loclabel(10,'Nbasis=')
read(10,*) c80tmp,nbasis
call loclabel(10,'Nindbasis=')
read(10,*) c80tmp,nindbasis
call loclabel(10,'Nprims=')
read(10,*) c80tmp,nprims
call loclabel(10,'Nshell=')
read(10,*) c80tmp,nshell
call loclabel(10,'Nprimshell=')
read(10,*) c80tmp,nprimshell
allocate(shelltype(nshell),shell2atom(nshell),shellcon(nshell))
call loclabel(10,'$Shell types');read(10,*)
read(10,*) shelltype
call loclabel(10,'$Shell centers');read(10,*)
read(10,*) shell2atom
call loclabel(10,'$Shell contraction degrees');read(10,*)
read(10,*) shellcon
allocate(primexp(nprimshell),concoeff(nprimshell))
call loclabel(10,'$Primitive exponents');read(10,*)
read(10,*) primexp
call loclabel(10,'$Contraction coefficients');read(10,*)
read(10,*) concoeff

!Copy information from local arrays to global arrays
allocate(shtype(nshell),shcen(nshell),shcon(nshell),primshexp(nprimshell),primshcoeff(nprimshell))
shtype=shelltype
shcen=shell2atom
shcon=shellcon
primshexp=primexp
primshcoeff=concoeff

!Multiwfn allows Cartesian and spherical basis functions mixed together. If any basis function is spherical, then isphergau=1
!Only the spherical ones will be treated specially
if (infomode==0) write(*,"(' The highest angular moment basis functions is ',a)") shtype2name(maxval(abs(shelltype))) 
isphergau=0
if (any(shelltype<=-2)) isphergau=1
if (any(abs(shelltype)>5).and.infomode==0) then
	write(*,"(a)") " Error: GTFs with angular moment higher than h are found, while Multiwfn currently only supports up to h. Press ENTER button to exit"
	read(*,*)
	stop
end if

if (nframe>1) then !Locate to specific frame
    write(c80tmp,*) iframe
    call loclabel(10,"@Frame= "//trim(adjustl(c80tmp)),ifound)
    if (ifound==0) then
        write(*,"(' Error: Unable to locate to frame',i8,' in the mwfn file!')") iframe
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        stop
    else
        do while(.true.)
            read(10,"(a)") c80tmp
            if (index(c80tmp,'Naelec=')/=0) then
                read(c80tmp,*) c80tmp2,naelec
            else if (index(c80tmp,'Nbelec=')/=0) then
                read(c80tmp,*) c80tmp2,nbelec
            else if (index(c80tmp,'E_tot=')/=0) then
                read(c80tmp,*) c80tmp2,totenergy
            else if (index(c80tmp,'VT_ratio=')/=0) then
                read(c80tmp,*) c80tmp2,virialratio
            else if (index(c80tmp,'cellv1=')/=0) then
                read(c80tmp,*) c80tmp2,cellv1
            else if (index(c80tmp,'cellv2=')/=0) then
                read(c80tmp,*) c80tmp2,cellv2
            else if (index(c80tmp,'cellv3=')/=0) then
                read(c80tmp,*) c80tmp2,cellv3
            else if (index(c80tmp,'Nindbasis=')/=0) then
                read(c80tmp,*) c80tmp2,nindbasis
            else if (index(c80tmp,'$Centers')/=0) then
                do icen=1,ncenter
                    read(10,*) inouse,a(icen)%name,a(icen)%index,a(icen)%charge,a(icen)%x,a(icen)%y,a(icen)%z
                end do
            else if (index(c80tmp,'Index=         1')/=0) then
                backspace(10)
                exit
            end if
        end do
    end if
end if

!Some treatments
nelec=naelec+nbelec
a%x=a%x/b2a
a%y=a%y/b2a
a%z=a%z/b2a
cellv1=cellv1/b2a
cellv2=cellv2/b2a
cellv3=cellv3/b2a
if (any(a%index==0).and.infomode==0) then
	write(*,*)
	write(*,*) "One or more dummy atoms (X) are found. Do you want to load them? (y/n)"
	write(*,"(a)") " Note: If all of them have corresponding basis functions, then they can be safely loaded. &
	&However, if some of them do not have basis functions, in general they should not be loaded, otherwise &
	&problems or crashes may occur when performing analyses based on wavefunction"
	read(*,*) selectyn
	if (selectyn=='n'.or.selectyn=='N') then
        do itime=1,ncenter
            backspace(10)
        end do
        ncenter=count(a%index/=0)
        deallocate(a);allocate(a(ncenter))
        icen=0
        do while(icen<=ncenter)
            read(10,*) inouse,c80tmp,itest
            if (itest/=0) then
                icen=icen+1
                backspace(10)
                read(10,*) inouse,a(icen)%name,a(icen)%index,a(icen)%charge,a(icen)%x,a(icen)%y,a(icen)%z
            end if
        end do
    end if
end if

!Load orbital information
!Note: Some basis maybe removed by linear dependence checking, hence the number of orbitals (nindbasis or 2*nindbasis) may be less than nbasis
!The information of those undefined orbitals (i.e. whose index is between nindbasis+1 and nbasis) are all set to zero
if (infomode==0) write(*,*) "Loading orbitals..."
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Restricted closed-shell or restricted open-shell
	nmo=nbasis
	allocate(amocoeff(nmo,nbasis))
    amocoeff=0
    norbload=nindbasis
    ires=1
else !Unrestricted open-shell wavefunction
    nmo=2*nbasis
	allocate(amocoeff(nbasis,nbasis),bmocoeff(nbasis,nbasis))
    amocoeff=0;bmocoeff=0
    norbload=2*nindbasis
    ires=0
end if
allocate(MOene(nmo),MOocc(nmo),MOtype(nmo),MOsym(nmo))
MOtype=0
MOocc=0
MOene=0
MOsym="?"
call loclabel(10,'Index=         1',irewind=0)
iaorb=0;iborb=0
do iorb=1,norbload
    if (iorb>1) read(10,*,iostat=ierror) !Skip space line between each orbital
    read(10,*,iostat=ierror) c80tmp !Skip index line
    if (index(c80tmp,"Index=")==0.or.ierror/=0) then !norbload is larger than actual number of recorded orbitals, often because linearly dependent basis functions are eliminated
        write(*,"(i6,' orbitals were actually loaded')") iorb-1
        exit
    else
        read(10,*) c80tmp,itype
        if (ires==1) then !Restricted
            MOtype(iorb)=itype
            read(10,*) c80tmp,MOene(iorb)
            read(10,*) c80tmp,MOocc(iorb)
            read(10,*) c80tmp,MOsym(iorb)
            read(10,*)
            read(10,*) amocoeff(:,iorb)
        else !Unrestricted
            if (itype==1) then
                iaorb=iaorb+1
                MOtype(iaorb)=itype
                read(10,*) c80tmp,MOene(iaorb)
                read(10,*) c80tmp,MOocc(iaorb)
                read(10,*) c80tmp,MOsym(iaorb)
                read(10,*)
                read(10,*) amocoeff(:,iaorb)
            else if (itype==2) then
                iborb=iborb+1
                MOtype(nbasis+iborb)=itype
                read(10,*) c80tmp,MOene(nbasis+iborb)
                read(10,*) c80tmp,MOocc(nbasis+iborb)
                read(10,*) c80tmp,MOsym(nbasis+iborb)
                read(10,*)
                read(10,*) bmocoeff(:,iborb)
            end if
        end if
    end if
end do
if (wfntype==1.or.wfntype==4) then !Unrestricted case, assign correct spin type for artificially filled orbitals
    MOtype(iaorb+1:nbasis)=1
    MOtype(nbasis+iborb+1:2*nbasis)=2
end if

close(10)

!!!!!! Loading have finished, now generate basis function information

!Backup spherical basis information (some of them may be Cartesian ones) with 5D suffix (of course, may be actually 7f, 9g, 11h...),
!convert them to Cartesian type temporarily to assign GTF information and calculate various integral matrices, at final stage recovery them back
if (isphergau==1) then
	allocate(shelltype5D(nshell))
	shelltype5D=shelltype
	where (shelltype<=-2) shelltype=-shelltype !Convert to Cartesian type
	nbasis5D=nbasis
	nbasis=0
	do i=1,nshell
		nbasis=nbasis+shtype2nbas(shelltype(i))
	end do
end if
allocate(shtypeCar(nbasis)) !Store shell information of Cartesian basis into global array, which may be used later
shtypeCar=shelltype
nbasisCar=nbasis

!Allocate space for arrays
allocate(b(nprims),CO(nmo,nprims),basshell(nbasis),bascen(nbasis),bastype(nbasis),primstart(nbasis),&
primend(nbasis),primconnorm(nprims),basstart(ncenter),basend(ncenter))
!Fill CObasa and CObasb
if (isphergau==0) then
	allocate(CObasa(nbasis,nbasis))
	CObasa=amocoeff
	if (ires==0) then
		allocate(CObasb(nbasis,nbasis))
		CObasb=bmocoeff
	end if
else if (isphergau==1) then
	allocate(CObasa(nbasis,nbasis),CObasa5D(nbasis5D,nbasis5D))
	CObasa5D=amocoeff
	CObasa=0D0
	if (ires==0) then
		allocate(CObasb(nbasis,nbasis),CObasb5D(nbasis5D,nbasis5D))
		CObasb5D=bmocoeff
		CObasb=0D0
	end if
	!Map 5D coefficient to 6D coefficient
	ipos5D=1
	ipos6D=1
	do ish=1,nshell
		ishtyp5D=shelltype5D(ish)
		ishtyp6D=shelltype(ish)
		numshorb5D=shtype2nbas(ishtyp5D)
		numshorb6D=shtype2nbas(ishtyp6D)
		if (ishtyp5D>=-1) then !S or P or other Cartesian shells
			CObasa(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasa5D(ipos5D:ipos5D+numshorb5D-1,:)
			if (ires==0) CObasb(ipos6D:ipos6D+numshorb6D-1,1:nbasis5D)=CObasb5D(ipos5D:ipos5D+numshorb5D-1,:)			
		else if (ishtyp5D==-2) then
			!5D->6D
			CObasa(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasa5D(ipos5D:ipos5D+4,:))
			if (ires==0) CObasb(ipos6D:ipos6D+5,1:nbasis5D)=matmul(conv5d6d,CObasb5D(ipos5D:ipos5D+4,:))
		else if (ishtyp5D==-3) then
			!7F->10F
			CObasa(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasa5D(ipos5D:ipos5D+6,:))
			if (ires==0) CObasb(ipos6D:ipos6D+9,1:nbasis5D)=matmul(conv7f10f,CObasb5D(ipos5D:ipos5D+6,:))
		else if (ishtyp5D==-4) then
			!9G->15G
			CObasa(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasa5D(ipos5D:ipos5D+8,:))
			if (ires==0) CObasb(ipos6D:ipos6D+14,1:nbasis5D)=matmul(conv9g15g,CObasb5D(ipos5D:ipos5D+8,:))
		else if (ishtyp5D==-5) then
			!11H->21H
			CObasa(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasa5D(ipos5D:ipos5D+10,:))
			if (ires==0) CObasb(ipos6D:ipos6D+20,1:nbasis5D)=matmul(conv11h21h,CObasb5D(ipos5D:ipos5D+10,:))
		end if
		ipos5D=ipos5D+numshorb5D
		ipos6D=ipos6D+numshorb6D
	end do
end if

if (infomode==0) write(*,*) "Converting basis function information to GTF information..."
!Distribute exponent, functype to every GTF and generate CO(:,:) from amocoeff/bmocoeff
!Fill: b,basshell,bascen,bastype,co,primstart,primend,primconnorm
k=1 !Current position of GTF
iexp=1
ibasis=1 !Current position of basis function
!Note: Below commented with !!! means the line associated to setting basis information
do i=1,nshell !Cycle each shell
	!write(*,*) i,shelltype(i),shellcon(i),shtype2nbas(shelltype(i)),k
	b(k:k+shellcon(i)*shtype2nbas(shelltype(i))-1)%center=shell2atom(i)
	basshell(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=i !!!Set basis attributed to which shell
	bascen(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=shell2atom(i) !!!Set basis attributed to which center
	do j=1,shtype2nbas(shelltype(i)) !Cycle each basis(orbital) in each shell
		b(k:k+shellcon(i)-1)%type=s2f(shelltype(i),j)
		bastype(ibasis)=s2f(shelltype(i),j) !!!Set basis type
		primstart(ibasis)=k !!!From where the GTFs attributed to ibasis'th basis
		primend(ibasis)=k+shellcon(i)-1 !!!To where the GTFs attributed to ibasis'th basis
		do l=1,shellcon(i) !Cycle each GTF in each basis in each shell
			b(k)%exp=primexp(iexp+l-1)
			tnormgau=normgau(b(k)%type,b(k)%exp)  !!!Normalization coefficient of GTFs
			temp=concoeff(iexp+l-1)  !!!Contraction coefficient of GTFs
			primconnorm(k)=temp*tnormgau !Combines contraction and normalization coefficient
			do imo=1,nmo
				if (ires==1) then !R or RO
					CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
				else !U
					if (isphergau==1) then
						if (imo<=nbasis5D) CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
						if (imo>nbasis5D) CO(imo,k)=CObasb(ibasis,imo-nbasis5D)*temp*tnormgau
					else
						if (imo<=nbasis) CO(imo,k)=CObasa(ibasis,imo)*temp*tnormgau
						if (imo>nbasis) CO(imo,k)=CObasb(ibasis,imo-nbasis)*temp*tnormgau
					end if
				end if
			end do
			k=k+1
		end do
		ibasis=ibasis+1
	end do
	iexp=iexp+shellcon(i)
end do

!Convert information from Cartesian basis to spherical basis
if (isphergau==1) then
	if (iloadasCart==1) then !For special purpose, keep Cartesian basis functions, e.g. converting spherical .fch/.molden to .47 file
		!Currently nbasis and dimension of all matrix (except for CO) correspond to full Cartesian case, but nmo &
		!and related arrays as well as CO still correspond to spherical harmonic case and thus need to be "expanded", the MO energies are artifically set to 100
		allocate(MOocc5D(nmo),MOene5D(nmo),MOtype5D(nmo),CO5D(nmo,nprims))
		MOocc5D=MOocc
		MOene5D=MOene
		MOtype5D=MOtype
		CO5D=CO
		deallocate(MOocc,MOene,MOtype,CO)
		if (ires==1) then
            nmo=nbasis !R, RO
		else
            nmo=2*nbasis !U
        end if
		allocate(MOocc(nmo),MOene(nmo),MOtype(nmo),CO(nmo,nprims))
		MOocc=0
		MOene=100
		CO=0
		if (ires==1) then !R, RO
			MOtype=0
			MOocc(1:nbasis5D)=MOocc5D
			MOene(1:nbasis5D)=MOene5D
			MOtype(1:nbasis5D)=MOtype5D
			CO(1:nbasis5D,:)=CO5D
		else !U
			MOtype(:nbasis)=1
			MOtype(nbasis+1:)=2
			MOocc(:nbasis5D)=MOocc5D(:nbasis5D)
			MOocc(nbasis+1:nbasis+nbasis5D)=MOocc5D(nbasis5D+1:)
			MOene(:nbasis5D)=MOene5D(:nbasis5D)
			MOene(nbasis+1:nbasis+nbasis5D)=MOene5D(nbasis5D+1:)
			MOtype(:nbasis5D)=MOtype5D(:nbasis5D)
			MOtype(nbasis+1:nbasis+nbasis5D)=MOtype5D(nbasis5D+1:)
			CO(:nbasis5D,:)=CO5D(:nbasis5D,:)
			CO(nbasis+1:nbasis+nbasis5D,:)=CO5D(nbasis5D+1:,:)
		end if
		isphergau=0
		
	else !Common case, transform to spherical basis
		if (infomode==0) write(*,*) "Back converting basis function information from Cartesian to spherical type..."
		!Recover spherical Gaussian basis function information
		nbasis=nbasis5D
		shelltype=shelltype5D
		ibasis=1
		do i=1,nshell
			basshell(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=i
			bascen(ibasis:ibasis+shtype2nbas(shelltype(i))-1)=shell2atom(i)
			do j=1,shtype2nbas(shelltype(i))
				bastype(ibasis)=s2f(shelltype(i),j)
				ibasis=ibasis+1
			end do
		end do
		deallocate(CObasa)
		allocate(CObasa(nbasis,nbasis))
		CObasa=CObasa5D
		if (ires==0) then
			deallocate(CObasb)
			allocate(CObasb(nbasis,nbasis))
			CObasb=CObasb5D
		end if
	end if
end if

!Generate basstart and basend
call bascen2basstart_end

!Generate one-particle density matrix for basis functions
if (igenP==1) then
	if (infomode==0) write(*,*) "Generating density matrix based on orbitals and occupation ..."
	call genP
end if

if (ifPBC==0) then !For PBC, Sbas will be generated after invoking "init_PBC" in "readinfile"
    if (infomode==0) write(*,*) "Generating overlap matrix..."
    call genSbas_curr
end if

!Output summary of present wavefunction
if (infomode==0) then
	write(*,*)
	write(*,"(' Total/Alpha/Beta electrons:',3f12.4)") nelec,naelec,nbelec
	write(*,"(' Net charge:',f12.5,'      Expected multiplicity:',i5)") sum(a(:)%charge)-nelec,nint(naelec-nbelec)+1
	write(*,"(' Atoms:',i7,',  Basis functions:',i7,',  GTFs:',i7)") ncenter,nbasis,nprims
	write(*,"(' Total energy:',f21.12,' Hartree,   Virial ratio:',f12.8)") totenergy,virialratio
	if (wfntype==0) then
		write(*,"(' This is a restricted single-determinant wavefunction')")
		write(*,"(' Orbitals from 1 to',i6,' are occupied')") nint(nelec/2)
	else if (wfntype==1) then
		write(*,"(' This is an unrestricted single-determinant wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are occupied')") 1,nbasis,1,nint(naelec)
		write(*,"(' Orbitals from ',i6,' to',i6,' are beta,  from',i6,' to',i6,' are occupied')") nbasis+1,nmo,nbasis+1,nbasis+nint(nbelec)
	else if (wfntype==2) then
		write(*,"(' This is a restricted open-shell wavefunction')")
		write(*,"(' Orbitals from',i6,' to',i6,' are doubly occupied')") 1,nint(nbelec)
		write(*,"(' Orbitals from',i6,' to',i6,' are singly occupied')") nint(nbelec)+1,nint(naelec)
	else if (wfntype==3) then
		write(*,"(' This is a restricted multiconfiguration wavefunction')")
	else if (wfntype==4) then
		write(*,"(' This is an unrestricted multiconfiguration wavefunction')")
		write(*,"(' Orbitals from ',i6,' to',i6,' are alpha, from',i6,' to',i6,' are beta')") 1,nbasis,nbasis+1,nmo
	end if
end if

call getHOMOidx !Find out index of HOMO, will be used in some cases
end subroutine




!!------------------- Read .cif file -------------------
! infomode=0: Output summary, =1: do not
!Sometimes symmetry operations are not explicitly given, but it can be inferred from _symmetry_space_group_name_H-M, &
!However Multiwfn ignores this currently. In principle symmetry operations can always be generated via http://cci.lbl.gov/sginfo/ or via spacegroupdata.py in cif2cell.py
!Ref. https://fossies.org/diffs/cp2k/3.0_vs_4.1/src/topology_cif.F-diff.html
subroutine readcif(name,infomode)
use deftype
use defvar
use util
use fparser
implicit real*8 (a-h,o-z)
integer infomode
character(len=*) name
character :: c200tmp*200,c80tmp*80,c80tmp2*80,c80tmp3*80,symopstr(500)*80,strarr(20)*30,var(3)=(/'x','y','z'/),testchar
real*8 vec1(3),vec2(3),vec1car(3),vec2car(3),val(3)
real*8,allocatable :: atmocc(:)
type(atomtype),allocatable :: a_dup(:)

ifiletype=17
ifPBC=3
open(10,file=name,status="old")

!!! Load cell information
if (infomode==0) write(*,*) "Loading cell information"
call loclabel(10,"_cell_length_a",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to find _cell_length_a field! Press ENTER button to exit"
    read(*,*);stop
end if
read(10,"(a)") c200tmp
call remove_parentheses(c200tmp)
read(c200tmp,*) c80tmp,asize
asize=asize/b2a
call loclabel(10,"_cell_length_b",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to find _cell_length_b field! Press ENTER button to exit"
    read(*,*);stop
end if
read(10,"(a)") c200tmp
call remove_parentheses(c200tmp)
read(c200tmp,*) c80tmp,bsize
bsize=bsize/b2a
call loclabel(10,"_cell_length_c",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to find _cell_length_c field! Press ENTER button to exit"
    read(*,*);stop
end if
read(10,"(a)") c200tmp
call remove_parentheses(c200tmp)
read(c200tmp,*) c80tmp,csize
csize=csize/b2a
call loclabel(10,"_cell_angle_alpha",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to find _cell_angle_alpha field! Press ENTER button to exit"
    read(*,*);stop
end if
read(10,"(a)") c200tmp
call remove_parentheses(c200tmp)
read(c200tmp,*) c80tmp,alpha
call loclabel(10,"_cell_angle_beta",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to find _cell_angle_beta field! Press ENTER button to exit"
    read(*,*);stop
end if
read(10,"(a)") c200tmp
call remove_parentheses(c200tmp)
read(c200tmp,*) c80tmp,beta
call loclabel(10,"_cell_angle_gamma",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to find _cell_angle_gamma field! Press ENTER button to exit"
    read(*,*);stop
end if
read(10,"(a)") c200tmp
call remove_parentheses(c200tmp)
read(c200tmp,*) c80tmp,gamma
call abc2cellv(asize,bsize,csize,alpha,beta,gamma)

!!! Load information of unique atoms
if (infomode==0) write(*,*) "Loading information of unique atoms"
call loclabel(10,"_atom_site_label",ifound)
if (ifound==0) then
    write(*,*) "Error: Unable to find _atom_site_label field! Press ENTER button to exit"
    read(*,*);stop
end if
do while(.true.)
    read(10,*) c80tmp
    if (c80tmp(1:1)/='_') exit
end do
ncenter_tmp=0
backspace(10)
do while(.true.) !Start loading atom part
    read(10,"(a)",iostat=ierror) c80tmp
    nunderline=strcharnum(c80tmp,'_') !I assume that if a line containing three underline, it must not atom name
    if (c80tmp==" ".or.index(c80tmp,'#')/=0.or.index(c80tmp,'loop_')/=0.or.nunderline>=3.or.ierror/=0) exit !.or.index(c80tmp,'_')/=0   Some atom names may contain _, so remove this criterion
    ncenter_tmp=ncenter_tmp+1
end do
if (allocated(a_tmp)) deallocate(a_tmp)
allocate(a_tmp(ncenter_tmp),atmocc(ncenter_tmp))
atmocc=1
if (infomode==0) write(*,"(a,i7)") " Number of symmetrically unique atoms:",ncenter_tmp

!Locate to "loop_" just before _atom_site_label
call loclabel(10,"_atom_site_label",ifound)
do while(.true.)
    backspace(10)
    read(10,"(a)") c80tmp
    if (index(c80tmp,"loop_")/=0) exit
    backspace(10)
end do
!Count how many atomic properties, and record index position of important labels
nlab=0
iocclab=0
do while(.true.)
    read(10,*) c80tmp
    if (index(c80tmp(1:1),'_')==0) exit
    nlab=nlab+1
    if (index(c80tmp,'_atom_site_label')/=0) iatmsitelab=nlab
    if (index(c80tmp,'_atom_site_type_symbol')/=0) iatmsitelab=nlab !If _atom_site_type_symbol is also given, use it prior to _atom_site_label
    if (index(c80tmp,'_atom_site_fract_x')/=0) ifrtxlab=nlab
    if (index(c80tmp,'_atom_site_fract_y')/=0) ifrtylab=nlab
    if (index(c80tmp,'_atom_site_fract_z')/=0) ifrtzlab=nlab
    if (index(c80tmp,'_atom_site_occupancy')/=0) iocclab=nlab
end do
backspace(10)
!Actual loading atomic information
do iatm=1,ncenter_tmp
    read(10,*) strarr(1:nlab)
    !Remove number and character after it in atom label if any
    call remove_parentheses(strarr(iatmsitelab))
    c80tmp=trim(strarr(iatmsitelab))
    do ichar=1,len_trim(c80tmp)
        if (iachar(c80tmp(ichar:ichar))>=48.and.iachar(c80tmp(ichar:ichar))<=57) then !This character is a digit
            c80tmp(ichar:)=" "
            exit
        end if
    end do
    a_tmp(iatm)%name=trim(c80tmp)
    !Detect element index
    call elename2idx(a_tmp(iatm)%name,a_tmp(iatm)%index)
    a_tmp(iatm)%charge=a_tmp(iatm)%index
    !Read fractional coordinates
    call remove_parentheses(strarr(ifrtxlab))
    call remove_parentheses(strarr(ifrtylab))
    call remove_parentheses(strarr(ifrtzlab))
    !write(*,"(i5,1x,a,1x,a,1x,a,1x,a)") iatm,trim(strarr(iatmsitelab)),trim(strarr(ifrtxlab)),trim(strarr(ifrtylab)),trim(strarr(ifrtzlab))
    read(strarr(ifrtxlab),*) a_tmp(iatm)%x
    read(strarr(ifrtylab),*) a_tmp(iatm)%y
    read(strarr(ifrtzlab),*) a_tmp(iatm)%z
    if (iocclab/=0) then
        call remove_parentheses(strarr(iocclab))
        read(strarr(iocclab),*,iostat=ierror) atmocc(iatm) !Some strange cif use ? for occupancy, in this case use 1.0 instead
        if (ierror/=0) atmocc(iatm)=1
    end if
end do
if (any(atmocc/=1)) then
    write(*,"(a)") " Warning: At least one atom has occupancy smaller than 1! All atoms will be loaded regardless of their occupancies"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
end if

!!! Read symmetry opteration and replicate atoms
!Only do this according to explicitly recorded symmetry operation, while gview and VESTA also loads point group like follows and generate symmetry opterations
!_symmetry_space_group_name_H-M   'P 21/c'
if (infomode==0) write(*,*) "Loading symmetry opteration and replicate atoms"
call loclabel(10,"_symmetry_equiv_pos_as_xyz",ifound)
if (ifound==0) call loclabel(10,"_space_group_symop_operation_xyz",ifound)
if (ifound==1) then
    do while(.true.) !Load to loop_ prior to symmetry operation field
        backspace(10)
        read(10,"(a)") c80tmp
        if (index(c80tmp,"loop_")/=0) exit
        backspace(10)
    end do
    nlab=0
    do while(.true.)
        read(10,*) c80tmp
        if (index(c80tmp,'_')==0) exit
        nlab=nlab+1
        if (index(c80tmp,'_symmetry_equiv_pos_as_xyz')/=0.or.index(c80tmp,'_space_group_symop_operation_xyz')/=0) isymop=nlab
    end do
    backspace(10)
    !Read symmetry operation string
    nsymopstr=0 !Number of symmetry operations
    do while(.true.)
        read(10,"(a)",iostat=ierror) c80tmp
        if (c80tmp==" ".or.index(c80tmp,'#')/=0.or.index(c80tmp,'_')/=0.or.ierror/=0) exit
        nsymopstr=nsymopstr+1
        do ichar=1,len_trim(c80tmp) !Replace all ' with space, make coordinate variables to lower case
            if (c80tmp(ichar:ichar)=="'") c80tmp(ichar:ichar)=" "
            if (c80tmp(ichar:ichar)=="X") c80tmp(ichar:ichar)="x"
            if (c80tmp(ichar:ichar)=="Y") c80tmp(ichar:ichar)="y"
            if (c80tmp(ichar:ichar)=="Z") c80tmp(ichar:ichar)="z"
        end do
        c80tmp=adjustl(c80tmp)
        if (nlab==1) then
            symopstr(nsymopstr)=trim(c80tmp)
        else !Load the second term, assume at most there are two terms in a line
            itmp=index(c80tmp,' ')
            symopstr(nsymopstr)=trim(adjustl(c80tmp(itmp+1:)))
        end if
    end do
    if (infomode==0) write(*,"(a,i4)") " Number of symmetry operations:",nsymopstr

    allocate(a_dup(ncenter_tmp*nsymopstr))
    call initf(3)
    ncenter_dup=ncenter_tmp
    a_dup(1:ncenter_dup)=a_tmp
    ifdoPBCx=1;ifdoPBCy=1;ifdoPBCz=1 !Because some PBC routines must be used in following codes now, PBC must be forced to be enabled at this stage
    !Cycle each unique atom and apply symmetry operation in turn, add it as actual atom if it doesn't have superposition to existing atom
    do isymop=1,nsymopstr
        itmp=index(symopstr(isymop),",")
        jtmp=index(symopstr(isymop),",",back=.true.)
        c80tmp=symopstr(isymop)(:itmp-1)
        c80tmp2=symopstr(isymop)(itmp+1:jtmp-1)
        c80tmp3=symopstr(isymop)(jtmp+1:)
        !Operation may be "1/2x+1/2y+z", however, to make parsing feasible, it should be replaced with "1/2*x+1/2*y+z", now we add * between number and x/y/z
        do itmp=2,80
			if (c80tmp(itmp:itmp)=='x'.or.c80tmp(itmp:itmp)=='y'.or.c80tmp(itmp:itmp)=='z') then
				testchar=c80tmp(itmp-1:itmp-1)
                if (iachar(testchar)>=48.and.iachar(testchar)<=57) then !Digit
					c80tmp(itmp+1:80)=c80tmp(itmp:79)
                    c80tmp(itmp:itmp)='*'
                end if
            end if
			if (c80tmp2(itmp:itmp)=='x'.or.c80tmp2(itmp:itmp)=='y'.or.c80tmp2(itmp:itmp)=='z') then
				testchar=c80tmp2(itmp-1:itmp-1)
                if (iachar(testchar)>=48.and.iachar(testchar)<=57) then !Digit
					c80tmp2(itmp+1:80)=c80tmp2(itmp:79)
                    c80tmp2(itmp:itmp)='*'
                end if
            end if
			if (c80tmp3(itmp:itmp)=='x'.or.c80tmp3(itmp:itmp)=='y'.or.c80tmp3(itmp:itmp)=='z') then
				testchar=c80tmp3(itmp-1:itmp-1)
                if (iachar(testchar)>=48.and.iachar(testchar)<=57) then !Digit
					c80tmp3(itmp+1:80)=c80tmp3(itmp:79)
                    c80tmp3(itmp:itmp)='*'
                end if
            end if
        end do
        !write(*,"(' Symmetry operation',i4,': ',a,1x,a,1x,a)") isymop,trim(c80tmp),trim(c80tmp2),trim(c80tmp3)
        call parsef(1,trim(c80tmp),var) !The string defines formula, var(1/2/3)=x/y/z specify variables in the formula. After that, we can use function evalf to evaluate result of the formula
        call parsef(2,trim(c80tmp2),var)
        call parsef(3,trim(c80tmp3),var)
        do iatm=1,ncenter_tmp !Duplicate unique atoms via symmetry operations
            vec1(1)=a_tmp(iatm)%x !Unique atom in fractional coordinate
            vec1(2)=a_tmp(iatm)%y
            vec1(3)=a_tmp(iatm)%z
            vec2(1)=evalf(1,vec1) !Obtain fractional coordinate of duplicated atom (vec2) by evaluating the formula defined above based on coordinate of the unique atom (vec1)
            vec2(2)=evalf(2,vec1)
            vec2(3)=evalf(3,vec1)
            !write(*,"(' Atom',i5,/,' Org fract.:',3f12.6,' Bohr',/,' New fract.:',3f12.6,' Bohr')") iatm,vec1,vec2
            call fract2Cart(vec2,vec2car)
            call move_to_cell(vec2car,vec2car) !Duplicated atom in Cartesian coordinate in the cell
        
            !Compare with existing atoms
            do jatm=1,ncenter_dup
                vec1(1)=a_dup(jatm)%x
                vec1(2)=a_dup(jatm)%y
                vec1(3)=a_dup(jatm)%z
                call fract2Cart(vec1,vec1car)
                call move_to_cell(vec1car,vec1car) !Existing atom in Cartesian coordinate in the cell
                call nearest_dist(vec1car,vec2car,dist)
                if (dist<0.05D0) exit
            end do
            if (jatm==ncenter_dup+1) then
                ncenter_dup=ncenter_dup+1
                !write(*,"(' Added new atom, current number:',i5)") ncenter_dup
                a_dup(ncenter_dup)=a_tmp(iatm)
                a_dup(ncenter_dup)%x=vec2(1)
                a_dup(ncenter_dup)%y=vec2(2)
                a_dup(ncenter_dup)%z=vec2(3)
            end if
        end do
    end do
    ncenter=ncenter_dup
    allocate(a(ncenter))
    a=a_dup(1:ncenter_dup)
    
else !Symmetry operations are not explicitly given
    if (infomode==0) write(*,"(a)") " Note: Symmetry operations are not explicitly given in the .cif file, therefore symmetry operations will not be taken into account"
    ncenter=ncenter_tmp
    allocate(a(ncenter))
    a=a_tmp
end if
deallocate(a_tmp)
close(10)

!!! Make all atoms to Cartesian coordinate and wrap into the cell
do iatm=1,ncenter
    vec1(1)=a(iatm)%x
    vec1(2)=a(iatm)%y
    vec1(3)=a(iatm)%z
    call fract2Cart(vec1,vec2)
    call move_to_cell(vec2,vec2)
    a(iatm)%x=vec2(1)
    a(iatm)%y=vec2(2)
    a(iatm)%z=vec2(3)
end do

if (infomode==0) write(*,"(' Totally',i8,' atoms')") ncenter
call guessnelec
end subroutine



!!---------- Interface of outputting cif file
subroutine outcif_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for outputting cif file, e.g. C:\ltwd.cif"
write(*,"(a)") " If press ENTER button directly,the system will be exported to "//trim(c200tmp)//".cif in current folder"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".cif"
call outcif(outname,10)
end subroutine
!!---------- Output current coordinate to cif file
subroutine outcif(outcifname,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outcifname
integer i,ifileid
real*8 rcoord(3),fcoord(3)
open(ifileid,file=outcifname,status="replace")
write(ifileid,"(a)") "#Generated by Multiwfn"
write(ifileid,"(a,f8.3)") "data_system"
call getcellabc(asize,bsize,csize,alpha,beta,gamma)
write(ifileid,"(a,f8.3)") "_cell_angle_alpha",alpha
write(ifileid,"(a,f8.3)") "_cell_angle_beta ",beta
write(ifileid,"(a,f8.3)") "_cell_angle_gamma",gamma
write(ifileid,"(a,f12.6)") "_cell_length_a",asize
write(ifileid,"(a,f12.6)") "_cell_length_b",bsize
write(ifileid,"(a,f12.6)") "_cell_length_c",csize
write(ifileid,"(a)") "loop_"
write(ifileid,"(a)") "_symmetry_equiv_pos_as_xyz"
write(ifileid,"(a)") "x,y,z"
write(ifileid,"(a)") "loop_"
write(ifileid,"(a)") "_atom_site_label"
write(ifileid,"(a)") "_atom_site_fract_x"
write(ifileid,"(a)") "_atom_site_fract_y"
write(ifileid,"(a)") "_atom_site_fract_z"
do i=1,ncenter
    rcoord(1)=a(i)%x
    rcoord(2)=a(i)%y
    rcoord(3)=a(i)%z
    call Cart2fract(rcoord,fcoord)
	write(ifileid,"(a,3f16.8)") ind2name(a(i)%index),fcoord(:)
end do

close(ifileid)
write(*,*) "Exporting cif file finished!"
end subroutine



!!---------- Interface of outputting gro file
subroutine outgro_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,*) "Input path for outputting gro file, e.g. C:\ltwd.gro"
write(*,"(a)") " If press ENTER button directly,the system will be exported to "//trim(c200tmp)//".gro in current folder"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".gro"
call outgro(outname,10)
end subroutine
!!---------- Output current coordinate to gro file
subroutine outgro(outgroname,ifileid)
use defvar
implicit real*8 (a-h,o-z)
character(len=*) outgroname
character tmpresname*3
integer ifileid

open(ifileid,file=outgroname,status="replace")
write(ifileid,"(a)") "Generated by Multiwfn"
write(ifileid,*) ncenter
do i=1,ncenter
	tmpresname="MOL"
	if (a(i)%resname/=" ") tmpresname=a(i)%resname
	write(ifileid,"(i5,a3,a7,i5,3f8.3)") a(i)%resid,tmpresname,a(i)%name,i,a(i)%x*b2a/10,a(i)%y*b2a/10,a(i)%z*b2a/10
end do
write(ifileid,"(9f12.6)") cellv1(1)*b2a/10,cellv2(2)*b2a/10,cellv3(3)*b2a/10,cellv1(2)*b2a/10,cellv1(3)*b2a/10,&
cellv2(1)*b2a/10,cellv2(3)*b2a/10,cellv3(1)*b2a/10,cellv3(2)*b2a/10
close(ifileid)
write(*,*) "Exporting gro file finished!"
end subroutine





!=======================================================================
!=======================================================================
!!!!!!!! Load parameters in settings.ini when boot up multiwfn !!!!!!!!!
!=======================================================================
!=======================================================================
subroutine loadsetting
use defvar
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,c200tmp2*200,settingpath*200
!Set default color of atomic spheres
atm3Dclr(:,1)=0.85D0
atm3Dclr(:,2)=0.6D0
atm3Dclr(:,3)=0.5D0
atm3Dclr(0,:)=(/0.3D0,  0.8D0,  0.85D0 /) !Bq
atm3Dclr(1,:)=(/0.95D0, 0.95D0, 0.95D0/) !H
atm3Dclr(5,:)=(/0.95D0, 0.7D0,  0.7D0 /) !B
atm3Dclr(6,:)=(/0.85D0, 0.85D0, 0.55D0/) !C
atm3Dclr(7,:)=(/0.5D0,  0.5D0,  1.0D0 /) !N
atm3Dclr(8,:)=(/1.0D0,  0.2D0,  0.2D0 /) !O
atm3Dclr(9,:)=(/0.6D0,  0.9D0,  0.9D0 /) !F
atm3Dclr(15,:)=(/0.9D0, 0.4D0,  0.0D0 /) !P
atm3Dclr(16,:)=(/0.9D0, 0.7D0,  0.1D0 /) !S
atm3Dclr(17,:)=(/0.1D0, 0.9D0,  0.1D0 /) !Cl

call getarg_str("-set",settingpath,ifound)

if (ifound==0) then
    inquire(file="settings.ini",exist=alive)
    if (alive) then
	    settingpath="settings.ini"
    else if (.not.alive) then
	    call getenv("Multiwfnpath",c80tmp)
	    if (isys==1) then
		    settingpath=trim(c80tmp)//"\settings.ini"
	    else if (isys==2) then
		    settingpath=trim(c80tmp)//"/settings.ini"
	    end if
	    inquire(file=settingpath,exist=alive)
	    if (.not.alive) then
		    write(*,"(a)") " WARNING: ""settings.ini"" was found neither in current folder nor in the path defined by ""Multiwfnpath"" &
		    environment variable. Now using default settings instead, and all your modifications (if any) on settings.ini will not take effect!"
            if (isys==1) then
				write(*,*) "Press ENTER button to continue"
                read(*,*)
            else
			    write(*,*)
            end if
		    return
	    end if
    end if
end if

open(20,file=settingpath,status="old")

! Below are the parameters can affect calculation results
call getarg_int("-uf",iuserfunc,ifound)
if (ifound==0) then
    call get_option_str(20,'iuserfunc=',c80tmp)
    if (c80tmp/=" ") read(c80tmp,*) iuserfunc
end if
call get_option_str(20,'refxyz=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) refx,refy,refz
call get_option_str(20,'iDFTxcsel=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iDFTxcsel
call get_option_str(20,'iKEDsel=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iKEDsel
call get_option_str(20,'uservar=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) uservar
call get_option_str(20,'uservar2=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) uservar2
call get_option_str(20,'ivdwprobe=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ivdwprobe
call get_option_str(20,'paircorrtype=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) paircorrtype
call get_option_str(20,'pairfunctype=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) pairfunctype
call get_option_str(20,'iautointgrid=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iautointgrid
call get_option_str(20,'radpot=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) radpot
call get_option_str(20,'sphpot=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) sphpot
call get_option_str(20,'radcut=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) radcut
call get_option_str(20,'expcutoff=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) expcutoff
call get_option_str(20,'expcutoff_PBC=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) expcutoff_PBC
call get_option_str(20,'RDG_maxrho=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) RDG_maxrho
call get_option_str(20,'RDGprodens_maxrho=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) RDGprodens_maxrho
call get_option_str(20,'IRI_rhocut=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) IRI_rhocut
call get_option_str(20,'ELF_addminimal=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ELF_addminimal
call get_option_str(20,'ELFLOL_type=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ELFLOL_type
call get_option_str(20,'iALIEdecomp=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iALIEdecomp
call get_option_str(20,'srcfuncmode=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) srcfuncmode
call get_option_str(20,'atomdenscut=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) atomdenscut
call get_option_str(20,'aug1D=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) aug1D
call get_option_str(20,'aug2D=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) aug2D
call get_option_str(20,'aug3D=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) aug3D
call get_option_str(20,'num1Dpoints=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) num1Dpoints
call get_option_str(20,'nprevorbgrid=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) nprevorbgrid
call get_option_str(20,'bndordthres=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) bndordthres
call get_option_str(20,'compthres=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) compthres
call get_option_str(20,'compthresCDA=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) compthresCDA
call get_option_str(20,'iCDAcomp=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iCDAcomp
call get_option_str(20,'ispheratm=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ispheratm
call get_option_str(20,'laplfac=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) laplfac
call get_option_str(20,'ipolarpara=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ipolarpara
call get_option_str(20,'ishowchgtrans=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ishowchgtrans
call get_option_str(20,'uESEinp=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) uESEinp
call get_option_str(20,'SpherIVgroup=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) SpherIVgroup
call get_option_str(20,'MCvolmethod=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) MCvolmethod
call get_option_str(20,'readEDF=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) readEDF
call get_option_str(20,'isupplyEDF=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) isupplyEDF
call get_option_str(20,'idelvirorb=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) idelvirorb
call get_option_str(20,'ifchprog=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ifchprog
call get_option_str(20,'ishowptESP=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ishowptESP
call get_option_str(20,'imolsurparmode=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) imolsurparmode
call get_option_str(20,'nPGmaxatm=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) nPGmaxatm
call get_option_str(20,'steric_addminimal=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) steric_addminimal
call get_option_str(20,'steric_potcutrho=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) steric_potcutrho
call get_option_str(20,'steric_potcons=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) steric_potcons
call get_option_str(20,'NICSnptlim=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) NICSnptlim
call get_option_str(20,'iplaneextdata=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iplaneextdata
call get_option_str(20,'igenP=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) igenP
call get_option_str(20,'iloadasCart=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iloadasCart
call get_option_str(20,'iloadGaugeom=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iloadGaugeom
call get_option_str(20,'maxloadexc=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) maxloadexc
call get_option_str(20,'iprintLMOorder=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iprintLMOorder
call get_option_str(20,'iMCBOtype=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iMCBOtype
call get_option_str(20,'ibasinlocmin=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ibasinlocmin
call get_option_str(20,'cfgcrossthres=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) cfgcrossthres
call get_option_str(20,'amIGMvdwscl=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) amIGMvdwscl
call get_option_str(20,'PBCnxnynz=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) PBCnx_in,PBCny_in,PBCnz_in
call get_option_str(20,'ifdoPBCxyz=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ifdoPBCx_in,ifdoPBCy_in,ifdoPBCz_in

call getarg_float("-ESPrhoiso",ESPrhoiso,ifound)
if (ifound==0) then
    call get_option_str(20,'ESPrhoiso=',c80tmp)
    if (c80tmp/=" ") read(c80tmp,*) ESPrhoiso
end if
call getarg_int("-ESPrhonlay",ESPrhonlay,ifound)
if (ifound==0) then
    call get_option_str(20,'ESPrhonlay=',c80tmp)
    if (c80tmp/=" ") read(c80tmp,*) ESPrhonlay
end if

!Below are the parameters involved in plotting
call get_option_str(20,'plotwinsize3D=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) plotwinsize3D
call get_option_str(20,'imodlayout=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) imodlayout
call get_option_str(20,'symbolsize=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) symbolsize
call get_option_str(20,'pleatmlabsize=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) pleatmlabsize
call get_option_str(20,'disshowlabel=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) disshowlabel
call get_option_str(20,'iatom_on_plane_far=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iatom_on_plane_far
call get_option_str(20,'iatmlabtype=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iatmlabtype
call get_option_str(20,'iatmlabtype3D=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iatmlabtype3D
call get_option_str(20,'graphformat=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) graphformat
call get_option_str(20,'graph1Dsize=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) graph1Dwidth,graph1Dheight
call get_option_str(20,'graph2Dsize=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) graph2Dwidth,graph2Dheight
call get_option_str(20,'graph3Dsize=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) graph3Dwidth,graph3Dheight
call get_option_str(20,'numdigxyz=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) numdigx,numdigy,numdigz
call get_option_str(20,'numdiglinexy=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) numdiglinex,numdigliney
call get_option_str(20,'numdigctr=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) numdigctr
call get_option_str(20,'fillcoloritpxy=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) fillcoloritpx,fillcoloritpy
call get_option_str(20,'itransparent=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) itransparent
call get_option_str(20,'isurfstyle=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) isurfstyle
call get_option_str(20,'bondRGB=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) bondclrR,bondclrG,bondclrB
call get_option_str(20,'atmlabRGB=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) atmlabclrR,atmlabclrG,atmlabclrB
call get_option_str(20,'CP_RGB=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) CP3n3RGB,CP3n1RGB,CP3p1RGB,CP3p3RGB
call get_option_str(20,'CP_RGB_2D=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) CP3n3RGB_2D,CP3n1RGB_2D,CP3p1RGB_2D,CP3p3RGB_2D
call get_option_str(20,'interbasin_RGB=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) IBSclrR,IBSclrG,IBSclrB
call get_option_str(20,'isoRGB_same=',c80tmp)
if (c80tmp/=" ") then
    read(c80tmp,*) clrRcub1same,clrGcub1same,clrBcub1same
    clrRcub1samemeshpt=clrRcub1same;clrGcub1samemeshpt=clrGcub1same;clrBcub1samemeshpt=clrBcub1same
end if
call get_option_str(20,'isoRGB_oppo=',c80tmp)
if (c80tmp/=" ") then
    read(c80tmp,*) clrRcub1oppo,clrGcub1oppo,clrBcub1oppo
    clrRcub1oppomeshpt=clrRcub1oppo;clrGcub1oppomeshpt=clrGcub1oppo;clrBcub1oppomeshpt=clrBcub1oppo
end if

!Set atom 3D color either according to external file or default setting
call get_option_str(20,'atmcolorfile=',c200tmp)
if (c200tmp/=" ".and.c200tmp/="none") then
    call remove_quotemark(c200tmp)
	inquire(file=c200tmp,exist=alive)
	if (alive) then
		write(*,"(' Note: Loading atom color settings from ',a)") trim(c200tmp)
		open(21,file=c200tmp,status="old")
		do iele=0,nelesupp
			read(21,*) inouse,atm3Dclr(iele,:)
		end do
		close(21)
        write(*,"(a,/)") " Done!"
    else
        write(*,"(a,/)") " Unable to find atomic color file: "//trim(c200tmp)
	end if
end if

!Load true type font filepath
call get_option_str(20,'ttfontfile=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) ttfontfile

!Load vdW radii from external file if the file is specified
call get_option_str(20,'vdwrfile=',c200tmp)
if (c200tmp/=" ".and.c200tmp/="none") then
    call remove_quotemark(c200tmp)
	inquire(file=c200tmp,exist=alive)
	if (alive) then
		write(*,"(' Note: Loading van der Waals radii from ',a)") trim(c200tmp)
		open(21,file=c200tmp,status="old")
        do while(.true.)
			read(21,"(a)",iostat=ierror) c200tmp
            if (ierror/=0.or.c200tmp==" ") exit
            if (index(c200tmp,'#')==0) then
				read(c200tmp,*) idx,tmpval
                vdwr(idx)=tmpval
            end if
        end do
		close(21)
        write(*,"(a,/)") " Done!"
    else
        write(*,"(a,/)") " Unable to find van der Waals radii file: "//trim(c200tmp)
	end if
end if

!Below are parameters about system
call getarg_int("-nt",nthreads,ifound)
if (ifound==0) then
    call get_option_str(20,'nthreads=',c80tmp)
    if (c80tmp/=" ") read(c80tmp,*) nthreads
end if
call get_option_str(20,'ompstacksize=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) ompstacksize
call get_option_str(20,'gaupath=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) gaupath
call get_option_str(20,'cubegenpath=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) cubegenpath
call get_option_str(20,'cubegendenstype=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) cubegendenstype
call get_option_str(20,'formchkpath=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) formchkpath
call get_option_str(20,'orcapath=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) orcapath
call get_option_str(20,'orca_2mklpath=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) orca_2mklpath
call get_option_str(20,'dftd3path=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) dftd3path
call testarg("-silent",isilent)
if (isilent==0) then
    call get_option_str(20,'isilent=',c80tmp)
    if (c80tmp/=" ") read(c80tmp,*) isilent
end if
call get_option_str(20,'iESPcode=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iESPcode
call get_option_str(20,'outmedinfo=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) outmedinfo
call get_option_str(20,'iwfntmptype=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iwfntmptype
call get_option_str(20,'iaddprefix=',c80tmp)
if (c80tmp/=" ") read(c80tmp,*) iaddprefix

call getarg_int("-ispecial",ispecial,ifound)
if (ifound==0) then
    call get_option_str(20,'ispecial=',c80tmp)
    if (c80tmp/=" ") read(c80tmp,*) ispecial
end if

!The last opened file name
call get_option_str(20,'lastfile=',c200tmp)
if (c200tmp/=" ") read(c200tmp,*) lastfile

close(20)

end subroutine



!!------- Delete files, cannot delete folder
subroutine delfile(delname)
use defvar
character(len=*) delname
character command*200
if (isys==1) then
    command="del /Q "//trim(delname)//" 2> NUL" !/Q means do not need to confirm when there is *. "2> NUL" redirects error (e.g. unable to find file) to black hole
else if (isys==2) then
    command="rm -f "//trim(delname)
end if
write(*,*) "Deleting "//trim(delname)
call system(trim(command))
end subroutine



!!------- Delete directory
subroutine deldir(delname)
use defvar
character(len=*) delname
character command*200
if (isys==1) then
    command="rmdir /Q /S "//trim(delname)//" 2> NUL"
else if (isys==2) then
    command="rm -rf "//trim(delname)
end if
write(*,*) "Deleting "//trim(delname)
call system(trim(command))
end subroutine



!!------- Copy a given file to another one
subroutine copyfile(file1,file2)
use defvar
character(len=*) file1,file2
character command*300

if (isys==1) then
	command="copy "//trim(file1)//' '//trim(file2)//" /Y > NUL"
else if (isys==2) then
	command="cp -f"//trim(file1)//' '//trim(file2)
end if
write(*,*) "Running: "//trim(command)
call system(trim(command))
end subroutine
