!=========================================================================================
!=========================================================================================
!============= Loading opened Gaussian+NBO or GENNBO output file (fileid=10) =============
!=========================================================================================
!=========================================================================================
!Usual using steps:
!use NAOmod
!open(10,file=filename,status="old")
!call checkNPA(ifound);if (ifound==0) return
!call loadNAOinfo
!call loadNAOatminfo  !Optional
!
!call checkDMNAO(ifound);if (ifound==0) return
!call loadDMNAO
!
!call checkNAOMO(ifound);if (ifound==0) return
!call loadNAOMO
!
!call checkAONAO(ifound);if (ifound==0) return
!call loadAONAO
!close(10)


!------ Check if natural population analysis information is presented
!ifound=0: NPA information is not found   ifound=1: NPA information is found
subroutine checkNPA(ifound)
use NAOmod
use util
integer ifound
call loclabel(10,"NATURAL POPULATIONS",ifound)
if (ifound==0) then
    close(10)
	write(*,"(a)") " Error: Cannot found natural population analysis information in the input file!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if
end subroutine



!------ Load detailed information about NAOs
subroutine loadNAOinfo
use NAOmod
use util
implicit real*8 (a-h,o-z)
character :: c80tmp*80,c80tmp2*80
character(len=3),allocatable :: shname_NAO_tmp(:),shset_NAO_tmp(:,:)
integer,allocatable :: shcen_NAO_tmp(:)

!Set iopshNAO
call loclabel(10,"*******         Alpha spin orbitals         *******",iopshNAO,1)
if (iopshNAO==0) write(*,*) "This is a closed shell calculation"
if (iopshNAO==1) write(*,*) "This is an open shell calculation"

!Set numNAO and ncenter_NAO
call loclabel(10,"NATURAL POPULATIONS")
read(10,*);read(10,*);read(10,*);read(10,*)
ilastspc=0
do while(.true.) !Find how many centers and how many NAOs. We need to carefully check where is ending
	read(10,"(a)") c80tmp
    !write(*,*) trim(c80tmp)
	if (c80tmp==' '.or.index(c80tmp,"low occupancy")/=0.or.index(c80tmp,"Population inversion found")/=0.or.index(c80tmp,"effective core potential")/=0) then
		if (ilastspc==1) then
			ncenter_NAO=iatm
			numNAO=iNAO
			exit
		end if
		ilastspc=1 !last line is space
	else
		read(c80tmp,*,iostat=ierror) iNAO,c80tmp2,iatm
        !Sometimes the content is " 940    H127  s      Ryd( 2s)...", there is no spacing between element symbol and atom index, use special treatment
        if (ierror/=0) then
            read(c80tmp,*) iNAO,c80tmp2
            do i=1,len_trim(c80tmp2)
                if ( iachar(c80tmp2(i:i))>=48 .and. iachar(c80tmp2(i:i))<=57 ) then !This is the first digit
                    read(c80tmp2(i:),*) iatm
                    exit
                end if
            end do
        end if
		ilastspc=0
	end if
end do
write(*,"(' The number of atoms:',i10)") ncenter_NAO
write(*,"(' The number of NAOs: ',i10)") numNAO

call dealloNAONBO(1)
allocate(NAOinit(ncenter_NAO),NAOend(ncenter_NAO),NAOcen(numNAO),NAOcenname(numNAO))
allocate(NAOset(numNAO,0:2),NAOtype(numNAO),NAOshname(numNAO),NAOocc(numNAO,0:2),NAOene(numNAO,0:2),atmname_NAO(ncenter_NAO))
NAOene=0 !For multiconfiguration method, there is no NAO energy, and thus NAO energies will be 0

!Load initial and ending index of NAOs corresponding to various atoms
call loclabel(10,"NATURAL POPULATIONS")
read(10,*);read(10,*);read(10,*);read(10,*)
ilastspc=1
do while(.true.)
	read(10,"(a)") c80tmp
	if (c80tmp/=' ') then
		read(c80tmp,*,iostat=ierror) iNAO,c80tmp2,iatm
        !Sometimes the content is " 940    H127  s      Ryd( 2s)...", there is no spacing between element symbol and atom index, use special treatment
        if (ierror/=0) then
            read(c80tmp,*) iNAO,c80tmp2
            do i=1,len_trim(c80tmp2)
                if ( iachar(c80tmp2(i:i))>=48 .and. iachar(c80tmp2(i:i))<=57 ) then !This is the first digit
                    read(c80tmp2(i:),*) iatm
                    exit
                end if
            end do
        end if
		NAOcen(iNAO)=iatm
		if (ilastspc==1) NAOinit(iatm)=iNAO
		ilastspc=0
	else
		NAOend(iatm)=iNAO
		if (iatm==ncenter_NAO) exit
		ilastspc=1
	end if
end do

!Load NAO information deriving from total density matrix
call loclabel(10,"NATURAL POPULATIONS")
read(10,*);read(10,*);read(10,*);read(10,*)
ispin=0
do iatm=1,ncenter_NAO
	do iNAO=NAOinit(iatm),NAOend(iatm)
		read(10,"(a)") c80tmp
		if (index(c80tmp,"Cor")/=0) then
			NAOset(iNAO,ispin)="Cor"
		else if (index(c80tmp,"Val")/=0) then
			NAOset(iNAO,ispin)="Val"
		else if (index(c80tmp,"Ryd")/=0) then
			NAOset(iNAO,ispin)="Ryd"
		end if
        !This code is even compatible with complicate case:   63   Th  1  f(0)   Val( 5f)     0.08546       0.08559
		read(c80tmp,*) c80tmp2,NAOcenname(iNAO),c80tmp2,NAOtype(iNAO)
        atmname_NAO(iatm)=NAOcenname(iNAO)
        do ic=1,len_trim(NAOtype(iNAO)) !Make sure that type name is always in lower case
            call uc2lc(NAOtype(iNAO)(ic:ic))
        end do
        !For total electron part, load orbital energies for closed shell case, while for open shell case the NAO energies will be loaded later
        i2=index(c80tmp,')',back=.true.)
		NAOshname(iNAO)=c80tmp(i2-3:i2-1)
        do ic=1,len_trim(NAOshname(iNAO)) !Make sure that shell name is always in lower case
            call uc2lc(NAOshname(iNAO)(ic:ic))
        end do
        if (iopshNAO==0) then
            read(c80tmp(i2+1:),*,iostat=ierror) NAOocc(iNAO,ispin),NAOene(iNAO,ispin) !For multiconfiguration method, there is no NAO energy
            if (ierror/=0) read(c80tmp(i2+1:),*,iostat=ierror) NAOocc(iNAO,ispin)
        else
		    read(c80tmp(i2+1:),*) NAOocc(iNAO,ispin)
            NAOene(iNAO,ispin)=0
        end if
        !mimic output of NPA output of NBO code
        !write(*,"(i6,a5,i5,2x,a7,a3,'(',a,')',f12.5,f14.5)") &
        !iNAO,NAOcenname(iNAO),NAOcen(iNAO),NAOtype(iNAO),NAOset(iNAO,ispin),NAOshname(iNAO),NAOocc(iNAO,ispin),NAOene(iNAO,ispin)
	end do
	read(10,*)
end do

if (iopshNAO==1) then
    !Load NAOset, NAOocc and NAOene of alpha spin
    do ispin=1,2
        call loclabel(10,"NATURAL POPULATIONS",ifound,0)
        read(10,*);read(10,*);read(10,*);read(10,*)
        do iatm=1,ncenter_NAO
	        do iNAO=NAOinit(iatm),NAOend(iatm)
		        read(10,"(a)") c80tmp
		        if (index(c80tmp,"Cor")/=0) then
			        NAOset(iNAO,ispin)="Cor"
		        else if (index(c80tmp,"Val")/=0) then
			        NAOset(iNAO,ispin)="Val"
		        else if (index(c80tmp,"Ryd")/=0) then
			        NAOset(iNAO,ispin)="Ryd"
		        end if
                read(c80tmp(i2+1:),*,iostat=ierror) NAOocc(iNAO,ispin),NAOene(iNAO,ispin) !For multiconfiguration method, there is no NAO energy
                if (ierror/=0) read(c80tmp(i2+1:),*,iostat=ierror) NAOocc(iNAO,ispin)
                !write(*,"(i6,a5,i5,2x,a7,a3,'(',a,')',f12.5,f14.5)") &
                !iNAO,NAOcenname(iNAO),NAOcen(iNAO),NAOtype(iNAO),NAOset(iNAO,ispin),NAOshname(iNAO),NAOocc(iNAO,ispin),NAOene(iNAO,ispin)
            end do
	        read(10,*)
        end do
    end do
end if

!Construct NAO shell information
numNAOsh=0
allocate(bassh_NAO(numNAO),shname_NAO_tmp(numNAO),shcen_NAO_tmp(numNAO),shset_NAO_tmp(numNAO,0:2))
do icen=1,ncenter_NAO
    ibeg=NAOinit(icen)
    iend=NAOend(icen)
    do iNAO=ibeg,iend
        if (iNAO>ibeg) then
            iold=0
            do ish=numNAOsh,1,-1 !Check if there is a old shell in this atom has identical name as current NAO
                if (NAOshname(iNAO)==shname_NAO_tmp(ish).and.shcen_NAO_tmp(ish)==icen) then
                    bassh_NAO(iNAO)=ish
                    iold=1 !There have been identical shell in this atom
                    exit
                end if
            end do
            if (iold==1) cycle
        end if
        numNAOsh=numNAOsh+1
        bassh_NAO(iNAO)=numNAOsh !The first NAO in this atom, must correspond to a new shell
		shname_NAO_tmp(numNAOsh)=NAOshname(iNAO)
        shset_NAO_tmp(numNAOsh,:)=NAOset(iNAO,:)
        shcen_NAO_tmp(numNAOsh)=icen
    end do
end do
allocate(shcen_NAO(numNAOsh),shname_NAO(numNAOsh),shset_NAO(numNAOsh,0:2))
shname_NAO=shname_NAO_tmp(1:numNAOsh)
shcen_NAO=shcen_NAO_tmp(1:numNAOsh)
shset_NAO(:,:)=shset_NAO_tmp(1:numNAOsh,:)

!do ish=1,numNAOsh
!    write(*,*) ish,shcen_NAO(ish),shname_NAO(ish)," ",shset_NAO(ish,0)
!end do
!write(*,*) bassh_NAO
end subroutine




!------ Load atom information from NAO information part
!ncenter will be set to ncenter_NAO, element and index property of "a" array will be filled
subroutine loadNAOatminfo
use NAOmod
use util
use defvar
implicit real*8 (a-h,o-z)
character :: c80tmp*80,c80tmp2*80

ncenter=ncenter_NAO
if (allocated(a)) deallocate(a)
allocate(a(ncenter))

call loclabel(10,"NATURAL POPULATIONS")
read(10,*);read(10,*);read(10,*);read(10,*)
do while(.true.)
	read(10,"(a)") c80tmp
	if (c80tmp/=' ') then
		read(c80tmp,*,iostat=ierror) inao,c80tmp2,iatm
        
        !Sometimes the content is " 940    H127  s      Ryd( 2s)...", there is no spacing between element symbol and atom index, use special treatment
        if (ierror/=0) then
            read(c80tmp,*) iNAO,c80tmp2
            do i=1,len_trim(c80tmp2)
                if ( iachar(c80tmp2(i:i))>=48 .and. iachar(c80tmp2(i:i))<=57 ) then !This is the first digit
                    read(c80tmp2(i:),*) iatm
                    c80tmp2(i:)=" "
                    exit
                end if
            end do
        end if
        
        call elename2idx(c80tmp2(1:2),iele)
        if (iele/=0) then
			a(iatm)%name=c80tmp2(1:2)
			a(iatm)%index=iele
        else
			write(*,"(a)") " Warning: Detected unrecognizable element name "//c80tmp2(1:2)//" !"
		end if
	else
		if (iatm==ncenter) exit
	end if
end do
end subroutine




!------ Check if DMNAO matrix is presented
!ifound=0: DMNAO is not found   ifound=1: DMNAO is found
subroutine checkDMNAO(ifound)
use NAOmod
use util
integer ifound
call loclabel(10,"NAO density matrix:",ifound)
if (ifound==0) then
    close(10)
	write(*,"(a)") " Error: Cannot found density matrix in NAO basis in the input file! You must use ""DMNAO"" keyword in the NBO"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if
end subroutine




!------ Load density matrix in NAO basis
subroutine loadDMNAO
use util
use NAOmod
implicit real*8 (a-h,o-z)
character c80tmp*80

call dealloNAONBO(2)
allocate(DMNAO(numNAO,numNAO))
call loclabel(10,"NAO density matrix:")
read(10,*);read(10,*);read(10,*)
read(10,"(a)") c80tmp
nskipcol=index(c80tmp,"- -")
backspace(10);backspace(10);backspace(10);backspace(10)
call readmatgau(10,DMNAO,0,"f8.4 ",nskipcol,8,3) !For open-shell case, the firstly printed density matrix is for alpha
call loclabel(10," Alpha spin orbitals ",iopenshell)
if (iopenshell==1) then
	allocate(DMNAOa(numNAO,numNAO),DMNAOb(numNAO,numNAO))
    DMNAOa=DMNAO
	write(*,*) "This is an open-shell calculation"
	call loclabel(10,"*******         Beta  spin orbitals         *******",ifound)
	call loclabel(10,"NAO density matrix:",ifound,0)
	call readmatgau(10,DMNAOb,0,"f8.4 ",nskipcol,8,3)
	DMNAO=DMNAOa+DMNAOb
end if
!call showmatgau(DMNAO,"Density matrix in NAO basis ",0,"f14.8",6)
end subroutine




!------ Check if NAOMO matrix is presented
!ifound=0: NAOMO is not found   ifound=1: NAOMO is found
subroutine checkNAOMO(ifound)
use NAOmod
use util
integer ifound
call loclabel(10,"MOs in the NAO basis:",ifound)
if (ifound==0) then
    close(10)
	write(*,"(a)") " Error: Cannot found coefficient matrix of NAOs in MOs from the input file, you should use ""NAOMO"" keyword in NBO module"
    write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
end subroutine




!------ Load NAOMO matrix
!numorb: The number of actual MOs (equals to NBsUse of Gaussian) 
subroutine loadNAOMO(numorb)
use util
use NAOmod
implicit real*8 (a-h,o-z)
character c80tmp*80
integer numorb

call dealloNAONBO(3)
allocate(NAOMO(numNAO,numorb))

call loclabel(10,"MOs in the NAO basis:")
!Check the columns that should be skipped during matrix reading, then return to title line
read(10,*);read(10,*);read(10,*)
read(10,"(a)") c80tmp
nskipcol=index(c80tmp,"- -")
backspace(10);backspace(10);backspace(10);backspace(10)
call readmatgau(10,NAOMO,0,"f8.4 ",nskipcol,8,3)

if (iopshNAO==1) then !Also load beta part
    allocate(NAOMOb(numNAO,numorb))
    call loclabel(10,"MOs in the NAO basis:",ifound,0)
    call readmatgau(10,NAOMOb,0,"f8.4 ",nskipcol,8,3)
end if

!call showmatgau(NAOMO,"MOs in the NAO basis",0,"f14.8",6)

end subroutine




!------ Check if AONAO matrix is presented, if presented, move to the line "NAOs in the AO basis:"
!ifound=0: AONAO is not found   ifound=1: AONAO is found
subroutine checkAONAO(ifound)
use NAOmod
use util
integer ifound
call loclabel(10,"NAOs in the AO basis:",ifound)
if (ifound==0) then
	close(10)
	write(*,"(a)") " Error: Cannot found coefficient matrix of NAOs in AO basis in the input file! &
    You must use ""AONAO"" keyword in NBO program!"
    write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
end subroutine




!----- Load AONAO matrix
!checkAONAO routine should be called to locate to the line containing "NAOs in the AO basis:" first
!numbas: Should be number of basis function before elimination of lineary dependency
!Since NAO orbitals is always generated by total density matrix, NBO only print it once even for open shell system
subroutine loadAONAO(numbas)
use defvar
use util
use NAOmod
implicit real*8 (a-h,o-z)
character c120tmp*120,c10tmp*10,form*10
integer numbas

write(*,"(a)") " Loading transformation matrix between original basis and NAO from input file..."
call dealloNAONBO(4)
allocate(AONAO(numbas,numNAO))
!Note: AONAO matrix in NBO output may be any number of columns (so that to ensure long data can be fully recorded), 
!so we must try to determine the actual number of rows and then use correct format to load it
!I assume that at least 5 columns and at most 8 columns
!Check number of columns that should be skipped and determining reading format, then return to title line
read(10,*);read(10,*);read(10,*)
read(10,"(a)") c120tmp
nskipcol=index(c120tmp,"- -")
nwid=len_trim(c120tmp)-index(trim(c120tmp)," ",back=.true.) !Test length of -------, may be different in different system
write(c10tmp,"(i2)") nwid+1
form='f'//trim(adjustl(c10tmp))//'.4' !The actual reading format
backspace(10);backspace(10);backspace(10);backspace(10)
!Start loading
read(10,*)
read(10,*)
read(10,"(a)") c120tmp
backspace(10);backspace(10);backspace(10)
if (index(c120tmp,'8')/=0) then !8 columns
    call readmatgau(10,AONAO,0,trim(form),nskipcol,8,3)
else if (index(c120tmp,'7')/=0) then !7 columns
	call readmatgau(10,AONAO,0,trim(form),nskipcol,7,3)
else if (index(c120tmp,'6')/=0) then !6 columns
	call readmatgau(10,AONAO,0,trim(form),nskipcol,6,3)
else if (index(c120tmp,'5')/=0) then !5 columns
	call readmatgau(10,AONAO,0,trim(form),nskipcol,5,3)
else
    write(*,*) "Error: Unable to load AONAO matrix!"
end if
!call showmatgau(AONAO)
end subroutine




!----- Deallocate all arrays involved in NAO/NBO analysis
!itype=0: Deallocate all NAO/NBO related arrays
!itype=1: Deallocate NAO information
!itype=2: Deallocate DMNAO matrix
!itype=3: Deallocate NAOMO matrix
!itype=4: Deallocate AONAO matrix
subroutine dealloNAONBO(itype)
use NAOmod
if (itype==0.or.itype==1) then
    if (allocated(atmname_NAO)) deallocate(atmname_NAO)
    if (allocated(NAOinit   )) deallocate(NAOinit   )
    if (allocated(NAOend    )) deallocate(NAOend    )
    if (allocated(NAOcen    )) deallocate(NAOcen    )
    if (allocated(NAOcenname)) deallocate(NAOcenname)
    if (allocated(NAOset    )) deallocate(NAOset    )
    if (allocated(NAOshname )) deallocate(NAOshname )
    if (allocated(NAOtype   )) deallocate(NAOtype   )
    if (allocated(NAOocc    )) deallocate(NAOocc    )
    if (allocated(NAOene    )) deallocate(NAOene    )
    if (allocated(shname_NAO)) deallocate(shname_NAO)
    if (allocated(shcen_NAO )) deallocate(shcen_NAO )
    if (allocated(bassh_NAO )) deallocate(bassh_NAO )
    if (allocated(shset_NAO )) deallocate(shset_NAO )
end if
if (itype==0.or.itype==2) then
    if (allocated(DMNAO     )) deallocate(DMNAO     )
    if (allocated(DMNAOa    )) deallocate(DMNAOa    )
    if (allocated(DMNAOb    )) deallocate(DMNAOb    )
end if
if (itype==0.or.itype==3) then
    if (allocated(NAOMO     )) deallocate(NAOMO     )
    if (allocated(NAOMOb    )) deallocate(NAOMOb    )
end if
if (itype==0.or.itype==4) then
    if (allocated(AONAO     )) deallocate(AONAO     )
end if
end subroutine