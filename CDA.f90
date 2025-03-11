!!!----------  Charge decomposition analysis (CDA) and extended CDA (ECDA)
! Closed-shell, unrestricted open-shell, natural orbitals are supported. Restricted open-shell is not supported
! For either complex or fragments, nbasis=nmo is always assumed! For open-shell, the nmo(or nmoCDA) used in this module is the number of either alpha or beta MOs, rather than their sum!
subroutine CDA
use util
use defvar
implicit real*8 (a-h,o-z)
!The matrices with b suffix are for beta part. Without it, denote total or alpha part (for open-shell cases)
!For the arrays contain multiple fragments, the 0th "fragment" is complex!!!
character,allocatable :: fragfilename(:)*200
real*8,allocatable :: ovlpbasmat(:,:) !Overlap matrix between all AO basis functions, loaded from complex output file
real*8,allocatable :: ovlpbasmatblk(:,:) !Off-diagonal block version of ovlpbasmat, diagonal blocks are zero
integer,allocatable :: elemidx(:,:) !Record elements index of all atoms, used for sanity check. The last element denotes frag#
real*8,allocatable :: atmpos(:,:,:) !Record coordinate of all atoms, used for sanity check. 1/2/3 of the second index corresponds to x/y/z. The last element denotes frag#
integer,allocatable :: natmCDA(:),nmoCDA(:) !The number of atoms and MOs in each fragment
integer,allocatable :: naelecCDA(:),nbelecCDA(:) !The number of alpha and beta electrons in each fragment
integer,allocatable :: iopsh(:),inatorb(:) !=1 means corresponding fragment is unrestricted open-shell / recording natural orbitals
integer,allocatable :: iRO(:) !=1 means corresponding fragment is RO-SCF fragment
!Below arrays relate to orbital information
real*8,allocatable :: occCDA(:,:),occCDAb(:,:) !Occupation numbers. The last index denotes frag#
real*8,allocatable :: eneCDA(:,:),eneCDAb(:,:) !MO energies, convert to eV after reading. The last index denotes frag#
real*8,allocatable :: cobasCDA(:,:,:),cobasCDAb(:,:,:) !Coefficient matrices in AO basis. The last index denotes frag#. The "CDA" in the name is to avoid conflict with global variable
!Below arrays will be computed rather than read
real*8,allocatable :: tmpmat(:,:),tmpmat2(:,:)
real*8,allocatable :: coFO(:,:),coFOb(:,:) !(i,j) denotes coefficient of FO i in complex orbitals j. FO index is sorted as 1~nmo1:~nmo2~nmo3...
real*8,allocatable :: FOcomp(:,:),FOcompb(:,:) !(i,j) denotes composition of FO i in complex orbitals j. Calculated by Mulliken or SCPA method. FO index is sorted as 1~nmo1~nmo2~nmo3...
real*8,allocatable :: FOovlpmat(:,:),FOovlpmatb(:,:) !Overlap matrix between all FOs. FO index is sorted as 1~nmo1~nmo2~nmo3...
real*8,allocatable :: dterm(:),bterm(:),rterm(:), dtermb(:),btermb(:),rtermb(:) !d,b,r defined in original paper
real*8 :: conncritleft=0.1D0,conncritright=0.1D0,degencrit=0.1D0,eneshiftA=0D0,eneshiftB=0D0,eneshiftcomp=0D0,eneintv=2D0
integer :: idrawMObar=1,iconnlogi=1,iout=6
character c80tmp*80,c80tmp2*80,selectyn
real*8 :: outthres=0 !Output threshold of CDA terms
!Some options relating to orbital interaction diagram
real*8 :: eneplotlow=-20,eneplothigh=5,complabshift=0.5D0
integer :: ilabelorbidx=1,ilabelcomp=1,labsize=40,ticknamesize=40,ispinplot=1

write(*,*) "Citation of generalized CDA method used in Multiwfn and original CDA method"
write(*,"(a)") " GCDA: Meng Xiao, Tian Lu, Generalized Charge Decomposition Analysis (GCDA) Method, J. Adv. Phys. Chem., 4, 111-124 (2015), http://dx.doi.org/10.12677/JAPC.2015.44013"
write(*,"(a)") " CDA: Stefan Dapprich, Gernot Frenking, J. Phys. Chem., 99, 9352-9362 (1995)"
write(*,*)

if (allocated(CObasa)) then !The input file contains basis function information
    igauout=0
else !Assume that the input file is Gaussian output file
    igauout=1
	open(10,file=filename,status="old")
	call loclabel(10,"NAtoms=",ifound) !Number of atoms
	read(10,*) c80tmp,ncenter
	do while(.true.)
		!Number of basis functions. However, sometimes it may be e.g. "Enter ECPInt, NBasis= 126.", so additional determination is needed
		call loclabel(10,"NBasis=",ifound,0)
		read(10,"(a)") c80tmp
		if (index(c80tmp,"ECP")==0) then
			backspace(10) 
			read(10,*) c80tmp,nbasis
			exit
		end if
	end do
	close(10)
	nmo=nbasis
end if

!====== Load complex (ifrag=0) and fragments
do while(.true.)
	write(*,*) "How many fragments do you want to define?  e.g. 2"
	read(*,*) nCDAfrag
	if (nCDAfrag>=2) exit
	write(*,*) "Error: The number of fragments must >= 2, input again"
end do
allocate(fragfilename(0:nCDAfrag),ovlpbasmat(nbasis,nbasis),elemidx(ncenter,0:nCDAfrag),atmpos(ncenter,3,0:nCDAfrag),&
natmCDA(0:nCDAfrag),nmoCDA(0:nCDAfrag),naelecCDA(0:nCDAfrag),nbelecCDA(0:nCDAfrag),iopsh(0:nCDAfrag),inatorb(0:nCDAfrag),iRO(0:nCDAfrag))
iRO=0
elemidx=0
atmpos=0D0
fragfilename=" "
! fragfilename(1)="examples\CDA\COBH3\CO.fch"
! fragfilename(2)="examples\CDA\COBH3\BH3.fch"

do ifrag=0,nCDAfrag !Here we first gather basic information of complex(ifrag=0) and fragments(ifrag>0)
	if (ifrag==0) then
		fragfilename(ifrag)=filename !The file loaded when Multiwfn boots up
		write(*,*) "Loading basic information of complex... Please wait"
	else
		if (igauout==0) then
            write(*,"(/,a,i4)") " Input .mwfn/.fch/.molden/.gms file of fragment",ifrag
            write(*,*) "e.g. C:\mol\CO.fch"
		else
            write(*,"(/,a,i4)") " Input Gaussian output file of fragment",ifrag
        end if
		do while(.true.)
			if (fragfilename(ifrag)/=" ") exit
			read(*,"(a)") fragfilename(ifrag)
			inquire(file=fragfilename(ifrag),exist=alive)
			if (alive) exit
			write(*,*) "Error: File not found, input again"
		end do
		write(*,*) "Loading basic information of this fragment... Please wait"
	end if
	if (igauout==1) then !Gaussian output file
		open(10,file=fragfilename(ifrag),status="old")
		call loclabel(10,"NAtoms=",ifound) !Number of atoms
		read(10,*) c80tmp,natmCDA(ifrag)
		call loclabel(10,"Input orientation:",ifound)
		if (ifound==0) call loclabel(10,"Z-Matrix orientation:",ifound) !Sometimes "Input orientation" doesn't occur, but "Z-Matrix orientation" occur, I don't know why
		read(10,*)
		read(10,*)
		read(10,*)
		read(10,*)
		read(10,*)
		do iatm=1,natmCDA(ifrag) !Load atom information
			ntmp=0
			if (ifrag>1) ntmp=sum(natmCDA(1:ifrag-1))
			jatm=ntmp+iatm
			read(10,*) nouse,elemidx(jatm,ifrag),nouse,atmpos(jatm,1:3,ifrag)
		end do
		do while(.true.)
			!Number of basis functions. However, sometimes it may be e.g. "Enter ECPInt, NBasis= 126.", so additional determination is needed
			call loclabel(10,"NBasis=",ifound,0)
			read(10,"(a)") c80tmp
			if (index(c80tmp,"ECP")==0) then
				backspace(10) 
				read(10,*) c80tmp,nbasistmp
				exit
			end if
		end do
		call loclabel(10,"NBsUse=",ifound,0) !NbsUse always equals to the actual number of MOs
		read(10,*) c80tmp,nmoCDA(ifrag)
		if (nmoCDA(ifrag)/=nbasistmp) then
			write(*,"(a)") " Error: The number of basis functions is unequal to the number of orbitals! Some basis functions may be &
			&eliminated during Gaussian calculation due to linear dependence problem. Please redo the Gaussian task with IOp(3/32=2) specified in route section"
			read(*,*)
			return
		end if
		call loclabel(10,"alpha electrons",ifound,1)
		read(10,*) naelecCDA(ifrag),c80tmp,c80tmp,nbelecCDA(ifrag)
		write(*,"(' Alpha electrons:',i7,'     Beta electrons:',i7,'     Multiplicity:',i3)") naelecCDA(ifrag),nbelecCDA(ifrag),naelecCDA(ifrag)-nbelecCDA(ifrag)+1
		if (ifrag==0) then
			write(*,*) "Loading overlap matrix of the complex, please wait..."
			call loclabel(10,"*** Overlap ***",ifound,0)
			call readmatgau(10,ovlpbasmat,1,"D14.6",7,5)
		end if
		call loclabel(10,"Orbital Coefficients",ifound,0)
		if (ifound==0) then
			write(*,*) "Error: Unable to find orbital coefficients in the file!"
			write(*,*) "Press ENTER button to exit"
			read(*,*)
			return
		end if
		read(10,*)
		call loclabel(10,"Orbital Coefficients",iopsh(ifrag),0) !If we can find "Orbital Coefficients" twice, that means this is unrestricted open-shell system
		call loclabel(10,"Natural Orbital Coefficients",inatorb(ifrag),1)
		close(10)
		if (iopsh(ifrag)==0.and.naelecCDA(ifrag)/=nbelecCDA(ifrag)) iRO(ifrag)=1 !This is a ROHF fragment
		
	else !The input file contains basis function information
		if (ifrag>0) then
			call dealloall(0)
			call readinfile(fragfilename(ifrag),1)
		end if
		natmCDA(ifrag)=ncenter
		do iatm=1,natmCDA(ifrag) !Load atom information
			ntmp=0
			if (ifrag>1) ntmp=sum(natmCDA(1:ifrag-1))
			jatm=ntmp+iatm
			elemidx(jatm,ifrag)=a(iatm)%index
			atmpos(jatm,1,ifrag)=a(iatm)%x
			atmpos(jatm,2,ifrag)=a(iatm)%y
			atmpos(jatm,3,ifrag)=a(iatm)%z
		end do
		nmoCDA(ifrag)=nbasis
		naelecCDA(ifrag)=nint(naelec)
		nbelecCDA(ifrag)=nint(nbelec)
		write(*,"(' Alpha electrons:',i7,'     Beta electrons:',i7,'     Multiplicity:',i3)") naelecCDA(ifrag),nbelecCDA(ifrag),naelecCDA(ifrag)-nbelecCDA(ifrag)+1
        if (ifrag==0) then
			call ask_Sbas_PBC !Compatible for periodic case
			ovlpbasmat=Sbas
        end if
		if (wfntype==1.or.wfntype==4) then
			iopsh(ifrag)=1
		else
			iopsh(ifrag)=0
		end if
		if (wfntype==3.or.wfntype==4) then
			inatorb(ifrag)=1
		else
			inatorb(ifrag)=0
		end if
		if (wfntype==2) iRO(ifrag)=1 !This is a ROHF fragment
	end if
	if (ifrag==0) then
		write(*,"(' The number of atoms in complex:',i10)") natmCDA(ifrag)
		write(*,"(' The number of basis functions in complex:',i10)") nmoCDA(ifrag)
		nbasisall=nbasis !Backup the nbasis and nmo of complex, they will be used later during reading fragments
		nmoall=nbasisall
	else
		write(*,"(' The number of basis functions in this fragment:',i10)") nmoCDA(ifrag)
		write(*,"(' The number of atoms in this fragment:',i10)") natmCDA(ifrag)
	end if
end do

!=== Check sanity of fragments
if (sum(natmCDA(1:))/=natmCDA(0)) then
	write(*,*)
	write(*,"(a)") " Error: The sum of the number of atoms in the fragments is not identical to complex!"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	return
end if
!Check atom consistency between complex and fragments
!write(*,*) "total"
!do iatm=1,ncenter_org
!    write(*,*) elemidx(iatm,0),ind2name(elemidx(iatm,0))
!end do
!do ifrag=1,nCDAfrag
!    ntmp=0
!    if (ifrag>1) ntmp=sum(natmCDA(1:ifrag-1))
!    write(*,*) "Frag",ifrag,"natom",natmCDA(ifrag)
!    do iatm=1,natmCDA(ifrag)
!	    jatm=ntmp+iatm
!        write(*,*) elemidx(jatm,ifrag),ind2name(elemidx(jatm,ifrag))
!    end do
!end do
do iatm=1,ncenter_org
	if (sum(elemidx(iatm,1:))/=elemidx(iatm,0)) then
		write(*,"(/,a)") " Error: The sequence of the atoms in the fragments is not consistent with that in complex, the result will be meaningless! Possible reasons:"
		write(*,"(a)") " 1 The fragment coordinates were not directly extracted from complex coordinate"
		write(*,"(a)") " 2 The loading sequence of the fragments is not identical to occurrence sequence of the fragments in complex"
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		return
	end if
end do
!Check consistency in coordinate
devmax=0
do ifrag=1,nCDAfrag
    ntmp=0
    if (ifrag>1) ntmp=sum(natmCDA(1:ifrag-1))
    do iatm=ntmp+1,natmCDA(ifrag)
        devtmp=dsqrt(sum((atmpos(iatm,:,0)-atmpos(iatm,:,ifrag))**2)) !Distance between fragment atom and complex atom
		!write(*,*) ifrag,iatm,devtmp
        if (devtmp>devmax) devmax=devtmp
    end do
end do
if (devmax>0.01D0) then
    write(*,"(/,' Maximum deviation of atomic coordinate between fragment and complex:',f10.3,' Angstrom')") devmax*b2a
    write(*,"(/,a)") " Warning: The coordinate of the fragments deviates from the complex distinctly! The result may be fully meaningless. &
    &Please check input files of your quantum chemistry code to make the coordinate of the fragments fully consistent with the complex. &
    &If you are a Gaussian user, please do not forget to add ""nosymm"" to avoid automatically translating and rotating the system to standard orientation"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
end if
!Check if the number of basis functions in complex is equal to the sum of that in all fragments
if (nmoCDA(0)/=sum(nmoCDA(1:))) then
	write(*,"(/,a)") " Error: The sum of the number of basis functions in all fragments is inconsistent with complex! &
	&Please carefully check the basis set you used in each calculation"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	return
end if

!Even for fragments, for simplicity, the size of the arrays are identical to complex, but of course only subblock will be filled. This will not waste too much memory
allocate(occCDA(nmoall,0:nCDAfrag),eneCDA(nmoall,0:nCDAfrag),cobasCDA(nbasisall,nmoall,0:nCDAfrag))
occCDA=0D0
eneCDA=0D0
cobasCDA=0D0
iopshCDA=0
!If any fragment is open-shell, the CDA will be seen as open-shell, and the closed-shell fragment will activate its beta part, namely copy information from alpha part
if (any(iopsh==1)) iopshCDA=1
if (iopshCDA==1) then
	allocate(occCDAb(nmoall,0:nCDAfrag),eneCDAb(nmoall,0:nCDAfrag),cobasCDAb(nbasisall,nmoall,0:nCDAfrag))
	occCDAb=0D0
	eneCDAb=0D0
	cobasCDAb=0D0
end if
!If any fragment use natural orbitals, the CDA will be performed in natural orbital mode
inatorbCDA=0
if (any(inatorb==1)) inatorbCDA=1

!=== Load information relating to orbitals
do ifrag=0,nCDAfrag
	if (ifrag==0) then
		write(*,"(/,a)") " Loading orbitals information for complex..."
	else if (ifrag>0) then
		write(*,"(/,a,i4,a)") " Loading orbitals information for fragment",ifrag,"..."
	end if
	if (igauout==1) then
		open(10,file=fragfilename(ifrag),status="old")
	else
		call dealloall(0)
		call readinfile(fragfilename(ifrag),1)
! 		if (nmo==2*nbasis) nmo=nbasis !The nmo used in present module is the number of either alpha or beta MOs rather than their sum as in other modules
	end if
	nmotmp=nmoCDA(ifrag)
	istart=1
	if (ifrag>1) istart=sum(nmoCDA(1:ifrag-1))+1
	iend=istart+nmotmp-1
		
	if (iopsh(ifrag)==0) then !Closed-shell complex or fragment
		if (igauout==1) then
			call loclabel(10,"Orbital Coefficients",ifound,1)
			call readmatgau(10,cobasCDA(istart:iend,istart:iend,ifrag),0,"f10.5",21,5,3-inatorb(ifrag)) !nbasis1+1:nbasis are empty
			!Set occupation number and energies
			if (inatorb(ifrag)==0) then !SCF orbital case
				occCDA(1:naelecCDA(ifrag),ifrag)=2
				call loclabel(10,"Orbital Coefficients",ifound,1)
				call readgauorbeig(eneCDA(1:nmotmp,ifrag),nmotmp,1) !Read MO energies
				eneCDA(:,ifrag)=eneCDA(:,ifrag)*au2eV
			else !Natural orbital case
				call loclabel(10,"Orbital Coefficients",ifound,1) 
				call readgauorbeig(occCDA(1:nmotmp,ifrag),nmotmp,0)
			end if
		else
			cobasCDA(istart:iend,istart:iend,ifrag)=CObasa
			eneCDA(1:nmotmp,ifrag)=MOene(:)*au2eV
			occCDA(1:nmotmp,ifrag)=MOocc(:)
		end if
! 		write(*,*) "Its orbital occupation numbers:"
! 		write(*,"(10f7.4)") occCDA(1:nmotmp,ifrag)
		if (iopshCDA==1) then !The whole CDA calculation is open-shell, so we duplicate this closed-shell fragment as alpha and beta part
			occCDA(:,ifrag)=occCDA(:,ifrag)/2D0
			occCDAb(:,ifrag)=occCDA(:,ifrag)
			eneCDAb(:,ifrag)=eneCDA(:,ifrag)
			cobasCDAb(:,:,ifrag)=cobasCDA(:,:,ifrag)
		end if
		
	else !Open-shell complex or fragment
		if (ifrag==0) write(*,*) "Note: The complex is an open-shell system"
		if (ifrag>0) write(*,*) "Note: This fragment is an open-shell system"
		if (igauout==1) then !Gaussian output file
			write(*,*) "Loading alpha MO cofficients..."
			call loclabel(10,"Orbital Coefficients",ifound,1)
			call readmatgau(10,cobasCDA(istart:iend,istart:iend,ifrag),0,"f10.5",21,5,3-inatorb(ifrag)) !For fragments, nmotmp+1:nmo are empty
			write(*,*) "Loading beta MO cofficients..."
			call loclabel(10,"Orbital Coefficients",ifound,0) !Beta
			call readmatgau(10,cobasCDAb(istart:iend,istart:iend,ifrag),0,"f10.5",21,5,3-inatorb(ifrag))
			!Set occupation number and energies
			if (inatorb(ifrag)==0) then !SCF orbital case
				occCDA(1:naelecCDA(ifrag),ifrag)=1
				occCDAb(1:nbelecCDA(ifrag),ifrag)=1
				call loclabel(10,"Orbital Coefficients",ifound,1) !Read MO energies
				call readgauorbeig(eneCDA(1:nmotmp,ifrag),nmotmp,1)
				call loclabel(10,"Orbital Coefficients",ifound,0) !Beta
				call readgauorbeig(eneCDAb(1:nmotmp,ifrag),nmotmp,1)
				eneCDA(:,ifrag)=eneCDA(:,ifrag)*au2eV
				eneCDAb(:,ifrag)=eneCDAb(:,ifrag)*au2eV
			else !Natural orbital case
				call loclabel(10,"Orbital Coefficients",ifound,1)
				call readgauorbeig(occCDA(1:nmotmp,ifrag),nmotmp,0)
				call loclabel(10,"Orbital Coefficients",ifound,0) !Beta
				call readgauorbeig(occCDAb(1:nmotmp,ifrag),nmotmp,0)
			end if
		else !Input file contains basis function information
			cobasCDA(istart:iend,istart:iend,ifrag)=CObasa
			cobasCDAb(istart:iend,istart:iend,ifrag)=CObasb
			eneCDA(1:nmotmp,ifrag)=MOene(1:nmotmp)*au2eV
			eneCDAb(1:nmotmp,ifrag)=MOene(nmotmp+1:2*nmotmp)*au2eV
			occCDA(1:nmotmp,ifrag)=MOocc(1:nmotmp)
			occCDAb(1:nmotmp,ifrag)=MOocc(nmotmp+1:2*nmotmp)
		end if
! 		write(*,*) "Its alpha orbital occupation numbers:"
! 		write(*,"(10f7.4)") occCDA(1:nmotmp,ifrag)
! 		write(*,*) "Its beta orbital occupation numbers:"
! 		write(*,"(10f7.4)") occCDAb(1:nmotmp,ifrag)
		if (ifrag==0) then
			write(*,*) "Flip electron spin of complex? (y/n)"
		else
			write(*,"(a,i4,a)") " Flip electron spin of fragment",ifrag,"? (y/n)"
		end if
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then !Swap the number of electrons and occupation numbers between alpha and beta orbitals
			tmpval=nbelecCDA(ifrag)
			nbelecCDA(ifrag)=naelecCDA(ifrag)
			naelecCDA(ifrag)=tmpval
			do imo=1,nmotmp
				tmpval=occCDAb(imo,ifrag)
				occCDAb(imo,ifrag)=occCDA(imo,ifrag)
				occCDA(imo,ifrag)=tmpval
				tmpval=eneCDAb(imo,ifrag)
				eneCDAb(imo,ifrag)=eneCDA(imo,ifrag)
				eneCDA(imo,ifrag)=tmpval
			end do
			allocate(tmpmat(nbasisall,nbasisall))
			tmpmat=cobasCDA(:,:,ifrag)
			cobasCDA(:,:,ifrag)=cobasCDAb(:,:,ifrag)
			cobasCDAb(:,:,ifrag)=tmpmat
			deallocate(tmpmat)
		end if	
	end if
	
	if (igauout==1) then
		close(10)
	else
		call dealloall(0)
		if (ifrag==nCDAfrag) then
			write(*,"(/,a,a)") " Reloading ",trim(firstfilename)
			call readinfile(firstfilename,1) !Recover to the first file
			if (nmo==2*nbasis) nmo=nbasis !The nmo used in present module is the number of either alpha or beta MOs rather than their sum as in other modules
		end if
	end if
end do

!=== Check sanity of electron number, if all fragments are not ROHF
if ( all(iRO==0) .and. (sum(naelecCDA(1:))/=naelecCDA(0) .or. sum(nbelecCDA(1:))/=nbelecCDA(0)) ) then
	write(*,"(/,a)") " Error: The sum of the number of alpha (beta) electrons of all fragments is unequal to the number of alpha (beta) electrons of the complex! &
	&Please check if electron spin flipping option was incorrectly set"
	read(*,*)
	return
end if

write(*,*)
write(*,*) "Calculating, please wait..."

!=== Calculate overlap matrix (FOovlpmat) between fragment orbitals
!The matrix multiplication formula of FOovlpmat is abstract, but absolutely correct
!The scheme of FOovlpmat
!      nmo1  nmo2  nmo3...
! nmo1
! nmo2
! nmo3
! ...

!Alpha or total part
allocate(tmpmat(nbasis,nbasis),tmpmat2(nbasis,nbasis),FOovlpmat(nbasis,nbasis),ovlpbasmatblk(nbasis,nbasis))
ovlpbasmatblk=ovlpbasmat !Off-diagonal block version of basis function overlap matrix, use this so that the FOovlpmat can be generated by only one line of code
tmpmat=0D0 !The nature of tmpmat is (nbasis1~nbasis2~nbasis3...nbasisX : nmo1~nmo2~nmo3...nmoX), it is a block diagonal matrix
itmp=0
do ifrag=1,nCDAfrag
	tmpmat=tmpmat+cobasCDA(:,:,ifrag)
	ovlpbasmatblk(itmp+1:itmp+nmoCDA(ifrag),itmp+1:itmp+nmoCDA(ifrag))=0
	itmp=itmp+nmoCDA(ifrag)
end do
!Transform overlap integral matrix from AO basis to FO basis, i.e. FOovlpmat=matmul(transpose(tmpmat),matmul(ovlpbasmatblk,tmpmat))
!The resultant FOovlpmat is off-diagonal block matrix to describe overlap integrals between fragment orbitals
write(*,*) "Calculating overlap integral matrix between fragment orbitals"
tmpmat2=matmul_blas(ovlpbasmatblk,tmpmat,nbasis,nbasis)
FOovlpmat=matmul_blas(tmpmat,tmpmat2,nbasis,nbasis,1,0)
forall (i=1:nbasis) FOovlpmat(i,i)=1 !All FO integrals within each fragment are zero, here set them to ideal value
write(*,*) "Calculating coefficient matrix of complex orbitals in fragment orbital basis"
allocate(coFO(nbasis,nmo))
call invmatsub(tmpmat,nbasis)
coFO=matmul_blas(tmpmat,cobasCDA(:,:,0),nbasis,nbasis)
deallocate(tmpmat)

!Beta part
if (iopshCDA==1) then
	!Calculate overlap matrix (FOovlpmatb) between beta fragment orbitals
	allocate(tmpmat(nbasis,nbasis),FOovlpmatb(nbasis,nbasis))
	tmpmat=0D0
	do ifrag=1,nCDAfrag
		tmpmat=tmpmat+cobasCDAb(:,:,ifrag)
	end do
	tmpmat2=matmul_blas(ovlpbasmatblk,tmpmat,nbasis,nbasis)
	FOovlpmatb=matmul_blas(tmpmat,tmpmat2,nbasis,nbasis,1,0)
	forall (i=1:nbasis) FOovlpmatb(i,i)=1 
	!Generate coefficient matrix of beta complex orbitals in beta fragment orbital basis
	allocate(coFOb(nbasis,nmo))
	call invmatsub(tmpmat,nbasis)
	coFOb=matmul_blas(tmpmat,cobasCDAb(:,:,0),nbasis,nbasis)
	deallocate(tmpmat)
end if
deallocate(ovlpbasmatblk,tmpmat2)

!===== Calculate complex orbital composition in FO basis method. Note: nmo=nbasis even for unrestricted cases in this module
write(*,*) "Calculating composition of complex orbitals..."
if (iCDAcomp==1) write(*,"(a)") " Note: Mulliken method is used. The method can be chosen by ""iCDAcomp"" parameter in settings.ini"
if (iCDAcomp==2) write(*,"(a)") " Note: SCPA method is used. The method can be chosen by ""iCDAcomp"" parameter in settings.ini"
allocate(FOcomp(nmo,nmo),tmpmat(nmo,nmo))

!$OMP PARALLEL DO SHARED(FOcomp) PRIVATE(imo,iFO,jFO,tmpval) schedule(dynamic) NUM_THREADS(nthreads)
do imo=1,nmo
	if (iCDAcomp==1) then !Mulliken
		do iFO=1,nmo
			tmpval=0
			do jFO=1,nmo
				tmpval=tmpval+coFO(jFO,imo)*FOovlpmat(iFO,jFO)
			end do
			FOcomp(iFO,imo)=coFO(iFO,imo)*tmpval
		end do
    else if (iCDAcomp==2) then !SCPA
		tmpval=sum(coFO(:,imo)**2)
		do iFO=1,nmo
			FOcomp(iFO,imo)=coFO(iFO,imo)**2/tmpval
        end do
    end if
end do
!$OMP END PARALLEL DO
if (iopshCDA==1) then !Beta part
	allocate(FOcompb(nmo,nmo))
	!$OMP PARALLEL DO SHARED(FOcompb) PRIVATE(imo,iFO,jFO,tmpval) schedule(dynamic) NUM_THREADS(nthreads)
	do imo=1,nmo
		if (iCDAcomp==1) then !Mulliken
			do iFO=1,nmo
				tmpval=0
				do jFO=1,nmo
					tmpval=tmpval+coFOb(jFO,imo)*FOovlpmatb(iFO,jFO)
				end do
				FOcompb(iFO,imo)=coFOb(iFO,imo)*tmpval
			end do
		else if (iCDAcomp==2) then !SCPA
			tmpval=sum(coFOb(:,imo)**2)
			do iFO=1,nmo
				FOcompb(iFO,imo)=coFOb(iFO,imo)**2/tmpval
			end do
		end if    
    end do
	!$OMP END PARALLEL DO
end if

if (nCDAfrag==2) then
	ifrag=1
	jfrag=2
	isel=0
	goto 10 !For two fragments cases, directly output CDA result
end if
1 continue


!============= Post-processing interface, but in fact CDA and ECDA are calculated by option 0
do while(.true.)
	write(*,*)
	write(*,"(a,f8.5)") " -3 Set threshold of printing CDA result in option 0, current:",outthres
	if (iout==6) write(*,*) "-2 Switch output destination (for options 0, 1, 6), current: Screen"
	if (iout==10) write(*,*) "-2 Switch output destination (for options 0, 1, 6), current: CDA.txt"
	write(*,*) "-1 Return to main menu"
	write(*,*) "0 Print CDA result and ECDA result"
	write(*,*) "1 Print full CDA result (All high-lying orbitals will be shown)" !For SCF, all high-lying result must be 0 and thus meaningless
	write(*,*) "2 Show fragment orbital contributions to specific complex orbital"
	write(*,*) "3 Export coefficient matrix of complex orbitals in fragment orbital basis"
	write(*,*) "4 Export overlap matrix between fragment orbitals"
	write(*,*) "5 Plot orbital interaction diagram" !Not available for natural orbitals
	write(*,*) "6 Decompose complex orbital contribution to CDA"
	read(*,*) isel
	
    if (isel==-3) then
		write(*,*) "Input threshold of absolute value of d or b or r term in printing, e.g. 0.01"
        read(*,*) outthres
	else if (isel==-2) then
		if (iout==6) then
			iout=10
		else if (iout==10) then
			iout=6
		end if
		write(*,*)
	else if (isel==-1) then
		return
	else if (isel==0.or.isel==1) then
		if (nCDAfrag==2) then
			ifrag=1
			jfrag=2
		else
			write(*,*) "Perform the analysis between which two fragments?  e.g. 1,3"
			read(*,*) ifrag,jfrag
			if (ifrag>nCDAfrag.or.jfrag>nCDAfrag) then
				write(*,*) "Error: The index exceeded valid range!"
				cycle
			end if
		end if
		if (iout==10.and.(isel==0.or.isel==1)) open(iout,file="CDA.txt",status="replace")
		!Calculate d,b,r terms
10		if (iopshCDA==0) refocc=2D0 !Reference orbital occupation number
		if (iopshCDA==1) refocc=1D0
		write(*,*)
		write(iout,*) "   ============= Charge decomposition analysis (CDA) result ============="
		write(iout,"(' d = The number of electrons donated from fragment',i3,' to fragment',i3)") ifrag,jfrag
		write(iout,"(' b = The number of electrons back donated from fragment',i3,' to fragment',i3)") jfrag,ifrag
		write(iout,*) "r = The number of electrons involved in repulsive polarization"
		!Alpha part or total
		if (iopshCDA==1) write(iout,"(/,a)") "                     **** Result for alpha electrons ****"
		if (allocated(dterm)) deallocate(dterm,bterm,rterm)
		allocate(dterm(nmo),bterm(nmo),rterm(nmo))
		dterm=0
		bterm=0
		rterm=0
		write(iout,*)
		write(iout,*) "   Orb.      Occ.          d           b        d - b          r"
		do iorb=1,nmo
			do iAidx=1,nmoCDA(ifrag) !Scan each FO pairs
				iAtmp=0
				if (ifrag>1) iAtmp=sum(nmoCDA(1:ifrag-1)) !Actual index of FO of the first fragment in all FOs
				iA=iAtmp+iAidx
				do iBidx=1,nmoCDA(jfrag)
					iBtmp=0
					if (jfrag>1) iBtmp=sum(nmoCDA(1:jfrag-1))
					iB=iBtmp+iBidx
					!occfac is defined by Tian Lu. =0 means both are virtual or occupied, or have identical occ.   >0 means donor->acceptor   <0 means donor<-acceptor
					!In particular, For SCF wavefunction, =1 occ_A->vir_B,  =-1 vir_A<-occ_B
					occfac=(occCDA(iAidx,ifrag)-occCDA(iBidx,jfrag))/refocc
					tmpval=occCDA(iorb,0)*occfac*coFO(iA,iorb)*coFO(iB,iorb)*FOovlpmat(iA,iB)
					!The overlap population between the two FOs in current complex orbital, is the twice of original r
					tmpval2=occCDA(iorb,0)*2*min(occCDA(iAidx,ifrag),occCDA(iBidx,jfrag))/refocc*coFO(iA,iorb)*coFO(iB,iorb)*FOovlpmat(iA,iB)
					if (occfac>0) dterm(iorb)=dterm(iorb)+tmpval
					if (occfac<0) bterm(iorb)=bterm(iorb)-tmpval  !The minus sign before tmpval cancelled the negative of occfac
					rterm(iorb)=rterm(iorb)+tmpval2
				end do
			end do
            if (abs(dterm(iorb))>outthres.or.abs(bterm(iorb))>outthres.or.abs(rterm(iorb))>outthres) then
				if (isel==1) then !Output all complex orbitals
					write(iout,"(i8,5f12.6)") iorb,occCDA(iorb,0),dterm(iorb),bterm(iorb),dterm(iorb)-bterm(iorb),rterm(iorb)
				else if (isel==0) then
					if (iorb<naelecCDA(0)+5.or.iorb==nmo) then
						write(iout,"(i8,5f12.6)") iorb,occCDA(iorb,0),dterm(iorb),bterm(iorb),dterm(iorb)-bterm(iorb),rterm(iorb)
					else
						write(iout,*) " ......"
						exit
					end if
				end if
            end if
		end do
		write(iout,*) "-------------------------------------------------------------------"
		write(iout,"(' Sum:   ',5f12.6)") sum(occCDA(:,0)),sum(dterm),sum(bterm),sum(dterm)-sum(bterm),sum(rterm)
        if (outthres/=0) write(*,*) "Note: The ""Sum"" includes all terms including those not printed above"
		!Beta part
		if (iopshCDA==1) then
			write(iout,"(/,a)") "                     **** Result for beta electrons ****"
			if (allocated(dtermb)) deallocate(dtermb,btermb,rtermb)
			allocate(dtermb(nmo),btermb(nmo),rtermb(nmo))
			dtermb=0
			btermb=0
			rtermb=0
			write(iout,*)
			write(iout,*) "   Orb.      Occ.          d           b        d - b          r"
			do iorb=1,nmo
				do iAidx=1,nmoCDA(ifrag) !Scan each FO pairs
					iAtmp=0
					if (ifrag>1) iAtmp=sum(nmoCDA(1:ifrag-1)) !Actual index of FO of the first fragment in all FOs
					iA=iAtmp+iAidx
					do iBidx=1,nmoCDA(jfrag)
						iBtmp=0
						if (jfrag>1) iBtmp=sum(nmoCDA(1:jfrag-1))
						iB=iBtmp+iBidx
						occfac=(occCDAb(iAidx,ifrag)-occCDAb(iBidx,jfrag))/refocc
						tmpval=occCDAb(iorb,0)*occfac*coFOb(iA,iorb)*coFOb(iB,iorb)*FOovlpmatb(iA,iB)
						if (occfac>0) dtermb(iorb)=dtermb(iorb)+tmpval
						if (occfac<0) btermb(iorb)=btermb(iorb)-tmpval
						tmpval=occCDAb(iorb,0)*2*min(occCDAb(iAidx,ifrag),occCDAb(iBidx,jfrag))/refocc*coFOb(iA,iorb)*coFOb(iB,iorb)*FOovlpmatb(iA,iB)
						rtermb(iorb)=rtermb(iorb)+tmpval
					end do
				end do
				if (abs(dtermb(iorb))>outthres.or.abs(btermb(iorb))>outthres.or.abs(rtermb(iorb))>outthres) then
					if (isel==1) then
						write(iout,"(i8,5f12.6)") iorb,occCDAb(iorb,0),dtermb(iorb),btermb(iorb),dtermb(iorb)-btermb(iorb),rtermb(iorb)
					else if (isel==0) then
						if (iorb<nbelecCDA(0)+5.or.iorb==nmo) then
							write(iout,"(i8,5f12.6)") iorb,occCDAb(iorb,0),dtermb(iorb),btermb(iorb),dtermb(iorb)-btermb(iorb),rtermb(iorb)
						else
							write(iout,*) " ......"
							exit
						end if
					end if
                end if
			end do
			write(iout,*) "-------------------------------------------------------------------"
			write(iout,"(' Sum:   ',5f12.6)") sum(occCDAb(:,0)),sum(dtermb),sum(btermb),sum(dtermb)-sum(btermb),sum(rtermb)
			if (outthres/=0) write(*,*) "Note: The ""Sum"" includes all terms including those not printed above"
			write(iout,*)
			write(iout,*) "Result for all electrons:"
			write(iout,"(' d=',f10.6,'  b=',f10.6,'  d - b =',f10.6,'  r=',f10.6)") sum(dterm)+sum(dtermb),sum(bterm)+sum(btermb),sum(dterm)-sum(bterm)+sum(dtermb)-sum(btermb),sum(rterm)+sum(rtermb)
		end if
		if (isel==0.and.inatorbCDA==1) write(iout,"(/,a)") " Note: The summation for d,b,r only covers the contribution of the orbitals shown above, to obtain full result, select option 1"
		write(iout,*)
		!ECDA
		if (isel==0) then
			if (inatorbCDA==1) then
				write(iout,"(a)") " Note: Complex or fragments are described by natural orbitals, ECDA analysis is thus skipped"
			else if (nCDAfrag>2) then
				write(iout,"(a)") " Note: ECDA is not applicable to more than two fragments cases, ECDA analysis is thus skipped"
			else if (any(iRO==1)) then
				write(iout,"(a)") " Note: ECDA is not applicable to restricted open-shell case and thus the analysis is skipped"
            else if (iCDAcomp==2) then !I found ECDA in combination with SCPA result in evidently wrong result, even for closed-shell
				write(iout,"(a)") " Note: ECDA is not applicable when SCPA is used to compute orbital composition, thus this analysis is skipped"
			else
				write(iout,*)
				write(iout,*) "     ========== Extended Charge decomposition analysis (ECDA) =========="
				if (iopshCDA==1) write(*,"(/,a)") "                     **** Result for alpha electrons ****"
				FOiocc_occ=0 !Contribution from all occupied orbitals of fragment i to all occupied complex orbitals
				FOivir_occ=0
				FOjocc_occ=0
				FOjvir_occ=0
				FOiocc_vir=0
				FOivir_vir=0
				FOjocc_vir=0
				FOjvir_vir=0
				itmp=0
				if (ifrag>1) itmp=sum(nmoCDA(1:ifrag-1)) !Actual index of FO of the first fragment in all FOs
				jtmp=0
				if (jfrag>1) jtmp=sum(nmoCDA(1:jfrag-1))
				naeleci=naelecCDA(ifrag)
				naelecj=naelecCDA(jfrag)
				do imo=1,nmo
					if (imo<=naelecCDA(0)) then !Complex orbital is occupied
						FOiocc_occ=FOiocc_occ+sum(FOcomp(itmp+1:itmp+naeleci,imo)) !Runs occupied orbitals of fragment i
						FOivir_occ=FOivir_occ+sum(FOcomp(itmp+naeleci+1:itmp+nmoCDA(ifrag),imo)) !Runs virtual orbitals of fragment i
						FOjocc_occ=FOjocc_occ+sum(FOcomp(jtmp+1:jtmp+naelecj,imo)) !Runs occupied orbitals of fragment j
						FOjvir_occ=FOjvir_occ+sum(FOcomp(jtmp+naelecj+1:jtmp+nmoCDA(jfrag),imo)) !Runs virtual orbitals of fragment j
					else !Complex orbital is virtual
						FOiocc_vir=FOiocc_vir+sum(FOcomp(itmp+1:itmp+naeleci,imo))
						FOivir_vir=FOivir_vir+sum(FOcomp(itmp+naeleci+1:itmp+nmoCDA(ifrag),imo))
						FOjocc_vir=FOjocc_vir+sum(FOcomp(jtmp+1:jtmp+naelecj,imo))
						FOjvir_vir=FOjvir_vir+sum(FOcomp(jtmp+naelecj+1:jtmp+nmoCDA(jfrag),imo))
					end if
				end do
				write(iout,*) "  Contribution to all occupied complex orbital:"
				write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") ifrag,FOiocc_occ*100,FOivir_occ*100
				write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") jfrag,FOjocc_occ*100,FOjvir_occ*100
				write(iout,*) "  Contribution to all virtual complex orbital:"
				write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") ifrag,FOiocc_vir*100,FOivir_vir*100
				write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") jfrag,FOjocc_vir*100,FOjvir_vir*100
				if (iopshCDA==0) ntmp=2
				if (iopshCDA==1) ntmp=1
					write(iout,"(' PL(',i2,') + CT(',i2,'->',i2,') =',f10.4,'      PL(',i2,') + CT(',i2,'->',i2,') =',f10.4)") ifrag,ifrag,jfrag,ntmp*FOiocc_vir,ifrag,jfrag,ifrag,ntmp*FOivir_occ
					write(iout,"(' PL(',i2,') + CT(',i2,'->',i2,') =',f10.4,'      PL(',i2,') + CT(',i2,'->',i2,') =',f10.4)") jfrag,ifrag,jfrag,ntmp*FOjvir_occ,jfrag,jfrag,ifrag,ntmp*FOjocc_vir
					write(iout,"(' The net electrons obtained by frag.',i2,' = CT(',i2,'->',i2,') - CT(',i2,'->',i2,') =',f10.4)") jfrag,ifrag,jfrag,jfrag,ifrag,ntmp*(FOiocc_vir-FOivir_occ)
				!Beta part
				if (iopshCDA==1) then
					write(*,"(/,a)") "                     **** Result for beta electrons ****"
					FOiocc_occb=0 !Contribution from all occupied orbitals of fragment i to all occupied complex orbitals
					FOivir_occb=0
					FOjocc_occb=0
					FOjvir_occb=0
					FOiocc_virb=0
					FOivir_virb=0
					FOjocc_virb=0
					FOjvir_virb=0
					itmp=0
					if (ifrag>1) itmp=sum(nmoCDA(1:ifrag-1)) !Actual index of FO of the first fragment in all FOs
					jtmp=0
					if (jfrag>1) jtmp=sum(nmoCDA(1:jfrag-1))
					nbeleci=nbelecCDA(ifrag)
					nbelecj=nbelecCDA(jfrag)
					do imo=1,nmo
						if (imo<=nbelecCDA(0)) then !Complex orbital is occupied
							FOiocc_occb=FOiocc_occb+sum(FOcompb(itmp+1:itmp+nbeleci,imo)) !Runs occupied orbitals of fragment i
							FOivir_occb=FOivir_occb+sum(FOcompb(itmp+nbeleci+1:itmp+nmoCDA(ifrag),imo)) !Runs virtual orbitals of fragment i
							FOjocc_occb=FOjocc_occb+sum(FOcompb(jtmp+1:jtmp+nbelecj,imo)) !Runs occupied orbitals of fragment j
							FOjvir_occb=FOjvir_occb+sum(FOcompb(jtmp+nbelecj+1:jtmp+nmoCDA(jfrag),imo)) !Runs virtual orbitals of fragment j
						else !Complex orbital is virtual
							FOiocc_virb=FOiocc_virb+sum(FOcompb(itmp+1:itmp+nbeleci,imo))
							FOivir_virb=FOivir_virb+sum(FOcompb(itmp+nbeleci+1:itmp+nmoCDA(ifrag),imo))
							FOjocc_virb=FOjocc_virb+sum(FOcompb(jtmp+1:jtmp+nbelecj,imo))
							FOjvir_virb=FOjvir_virb+sum(FOcompb(jtmp+nbelecj+1:jtmp+nmoCDA(jfrag),imo))
						end if
					end do
					write(iout,*) "  Contribution to all occupied complex orbital:"
					write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") ifrag,FOiocc_occb*100,FOivir_occb*100
					write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") jfrag,FOjocc_occb*100,FOjvir_occb*100
					write(iout,*) "  Contribution to all virtual complex orbital:"
					write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") ifrag,FOiocc_virb*100,FOivir_virb*100
					write(iout,"(' Occupied, virtual orbitals of fragment',i3,': ',f12.4,'%   ',f12.4,'%')") jfrag,FOjocc_virb*100,FOjvir_virb*100
					write(iout,"(' PL(',i2,') + CT(',i2,'->',i2,') =',f10.4,'      PL(',i2,') + CT(',i2,'->',i2,') =',f10.4)") ifrag,ifrag,jfrag,FOiocc_virb,ifrag,jfrag,ifrag,FOivir_occb
					write(iout,"(' PL(',i2,') + CT(',i2,'->',i2,') =',f10.4,'      PL(',i2,') + CT(',i2,'->',i2,') =',f10.4)") jfrag,ifrag,jfrag,FOjvir_occb,jfrag,jfrag,ifrag,FOjocc_virb
					write(iout,"(' The net electrons obtained by frag.',i2,' = CT(',i2,'->',i2,') - CT(',i2,'->',i2,') =',f10.4)") jfrag,ifrag,jfrag,jfrag,ifrag,FOiocc_virb-FOivir_occb
					write(iout,*)
					write(iout,"(' CT(',i2,'->',i2,') - CT(',i2,'->',i2,') for all electrons:',f10.4)") ifrag,jfrag,jfrag,ifrag,(FOiocc_vir-FOivir_occ)+(FOiocc_virb-FOivir_occb)
				end if
			end if
			write(iout,*)
		end if
		if (iout==10) then
			close(10)
			write(*,*) "Done! The result has been outputted to CDA.txt in current folder"
			write(*,*)
		end if 
		if (nCDAfrag==2) goto 1

	else if (isel==2) then !Compositions have been computed by Mulliken or SCPA method before
		do while(.true.)
			write(*,*) "Input the index of complex orbital you are interested in, e.g. 6"
			write(*,*) "To exit, input 0"
			write(*,"(a)") " Note: If you want to obtain the composition of an orbital (e.g. orbital 5) of a fragment (e.g. fragment 2) in all complex orbitals, you can input e.g. 2,5"
			read(*,"(a)") c80tmp
			if (index(c80tmp,',')==0) then !Output composition of a complex orbital
				read(c80tmp,*) iorb
				if (iorb==0) then
					exit
				else if (iorb<0.or.iorb>nmo) then
					write(*,"(a,i7,a,i7,/)") " Error: Orbital index should between",1," and",nmo
					cycle
				else
					write(*,"(a,f5.1,a,/)") " Note: Only the fragment orbitals with contribution >",compthresCDA," % will be shown below, &
                    &the threshold can be changed by ""compthresCDA"" in settings.ini"
					if (iopshCDA==0) then
						write(*,"(' Occupation number of orbital',i6,' of the complex:',f12.8)") iorb,occCDA(iorb,0)
                        sumcontri=0
						do ifrag=1,nCDAfrag
							do iFOidx=1,nmoCDA(ifrag)
								itmp=0
								if (ifrag>1) itmp=sum(nmoCDA(1:ifrag-1))
								iFO=itmp+iFOidx
								if (abs(FOcomp(iFO,iorb))*100>=compthresCDA) then
                                    write(*,"(' Orbital',i6,' of fragment',i3,', Occ:',f8.5,'    Contribution:',f8.2,' %')") &
								    iFOidx,ifrag,occCDA(iFOidx,ifrag),FOcomp(iFO,iorb)*100
                                    sumcontri=sumcontri+FOcomp(iFO,iorb)
                                end if
							end do
						end do
                        write(*,"(' Sum of values shown above:',f10.2,' %')") sumcontri*100
					else if (iopshCDA==1) then
                        sumcontri=0
						write(*,*) "                           **** Alpha orbitals ****"
						write(*,"(' Occupation number of alpha orbital',i6,' of the complex:',f12.8)") iorb,occCDA(iorb,0)
						do ifrag=1,nCDAfrag
							do iFOidx=1,nmoCDA(ifrag)
								itmp=0
								if (ifrag>1) itmp=sum(nmoCDA(1:ifrag-1))
								iFO=itmp+iFOidx
								if (abs(FOcomp(iFO,iorb))*100>=compthresCDA) then
                                    write(*,"(' Alpha orbital',i6,' of fragment',i3,', Occ:',f8.5,'    Contribution:',f8.2,' %')") &
								    iFOidx,ifrag,occCDA(iFOidx,ifrag),FOcomp(iFO,iorb)*100
                                    sumcontri=sumcontri+FOcomp(iFO,iorb)
                                end if
							end do
						end do
                        write(*,"(' Sum of values shown above:',f10.2,' %',/)") sumcontri*100
                        
                        sumcontri=0
						write(*,*) "                           **** Beta orbitals ****"
						write(*,"(' Occupation number of beta orbital ',i6,' of the complex:',f12.8)") iorb,occCDAb(iorb,0)
						do ifrag=1,nCDAfrag
							do iFOidx=1,nmoCDA(ifrag)
								itmp=0
								if (ifrag>1) itmp=sum(nmoCDA(1:ifrag-1))
								iFO=itmp+iFOidx
								if (abs(FOcompb(iFO,iorb))*100>=compthresCDA) then
                                    write(*,"(' Beta orbital ',i6,' of fragment',i3,', Occ:',f8.5,'    Contribution:',f8.2,' %')") &
								    iFOidx,ifrag,occCDAb(iFOidx,ifrag),FOcompb(iFO,iorb)*100
                                    sumcontri=sumcontri+FOcompb(iFO,iorb)
                                end if
							end do
						end do
                        write(*,"(' Sum of values shown above:',f10.2,' %')") sumcontri*100
					end if
					write(*,*)
				end if
			else !Output composition of a fragment orbital in all complex orbitals
				read(c80tmp,*) ifrg,ifrgorb
				if (ifrg>nCDAfrag) then
					write(*,"(a,i3,/)") " Error: The fragment index must <=",nCDAfrag
				else if (ifrgorb<=0.or.ifrgorb>nmoCDA(ifrg)) then
					write(*,"(a,/)") " Error: The index of the fragment orbital exceeded valid range!"
				else
					write(*,"(' Composition of orbital',i6,' of fragment',i3,' in all complex orbitals:')") ifrgorb,ifrg
					itmp=0
					if (ifrg>1) itmp=sum(nmoCDA(1:ifrg-1))
					ifrgorb=ifrgorb+itmp
					if (iopshCDA==0) then
						do iorb=1,nmo
							write(*,"(' Contribution to complex orbital',i6,' (Occ:',f8.5,'):',f10.4,'%')") iorb,occCDA(iorb,0),FOcomp(ifrgorb,iorb)*100
						end do
					else if (iopshCDA==1) then
						write(*,*) "                           **** Alpha orbitals ****"
						do iorb=1,nmo
							write(*,"(' Contribution to complex orbital',i6,' (Occ:',f8.5,'):',f10.4,'%')") iorb,occCDA(iorb,0),FOcomp(ifrgorb,iorb)*100
						end do
						write(*,*)
						write(*,*) "                           **** Beta orbitals ****"
						do iorb=1,nmo
							write(*,"(' Contribution to complex orbital',i6,' (Occ:',f8.5,'):',f10.4,'%')") iorb,occCDAb(iorb,0),FOcompb(ifrgorb,iorb)*100
						end do
					end if
					write(*,*)
				end if
			end if
		end do
		
	else if (isel==3) then !Output coefficient matrix of complex orbitals in FO basis
		open(10,file="coFO.txt",status="replace")
		if (iopshCDA==1) then
			call showmatgau(coFO,"Alpha coefficient matrix",0,"f14.8",10)
			call showmatgau(coFOb,"Beta coefficient matrix",0,"f14.8",10)
		else
			call showmatgau(coFO,"",0,"f14.8",10)
		end if
		close(10)
		write(*,"(a)") "Done! The coefficient matrix has been outputted to coFO.txt in current folder. &
		&The element (i,n) is the coefficient of fragment orbital i in complex orbital n."
		do ifrag=1,nCDAfrag
			istart=1
			if (ifrag>1) istart=sum(nmoCDA(1:ifrag-1))+1
			iend=istart+nmoCDA(ifrag)-1
			write(*,"(' FO from',i7,' to',i7,' correspond to the orbitals of fragment',i3)") istart,iend,ifrag
		end do
		write(*,*)
	else if (isel==4) then
		open(10,file="ovlpint.txt",status="replace")
		if (iopshCDA==1) then
			call showmatgau(FOovlpmat,"Overlap matrix between alpha FOs",0,"f14.8",10)
			call showmatgau(FOovlpmatb,"Overlap matrix between beta FOs",0,"f14.8",10)
		else
			call showmatgau(FOovlpmat,"Overlap matrix between FOs",0,"f14.8",10)
		end if
		close(10)
		write(*,"(a)") " Done! The matrix has been outputted to ovlpint.txt in current folder. &
		&The element (i,j) is the overlap integral between fragment orbital i and j"
		do ifrag=1,nCDAfrag
			istart=1
			if (ifrag>1) istart=sum(nmoCDA(1:ifrag-1))+1
			iend=istart+nmoCDA(ifrag)-1
			write(*,"(' FO from',i7,' to',i7,' correspond to the orbitals of fragment',i3)") istart,iend,ifrag
		end do
		write(*,*)
		
	else if (isel==5) then
		do while(.true.)
			if (inatorbCDA==1) then
				write(*,"(a,/)") " Error: This function is unavailable for natural orbitals!"
				exit
			end if
			write(*,*)
			write(*,*) "                  ----- Orbital interaction diagram -----"
			write(*,*) "0 Return"
			write(*,*) "1 Plot the diagram now!"
			write(*,*) "2 Save the diagram to current folder"
			write(*,"(a,2f12.4)") " 3 Set the energy range to be plotted, current (eV):",eneplotlow,eneplothigh
			if (iconnlogi==1) write(*,"(' 4 Set the rule for linking and drawing orbital bars, current: ',f5.1,',or,',f5.1)") conncritleft*100,conncritright*100
			if (iconnlogi==2) write(*,"(' 4 Set the rule for linking and drawing orbital bars, current: ',f5.1,',and,',f5.1)") conncritleft*100,conncritright*100
			if (iopshCDA==1.and.ispinplot==1) write(*,*) "5 Switch which type of orbitals to be plotted, current: Alpha"
			if (iopshCDA==1.and.ispinplot==2) write(*,*) "5 Switch which type of orbitals to be plotted, current: Beta"
			if (ilabelorbidx==0) write(*,*) "6 Enable labelling orbital indices"
			if (ilabelorbidx==1) write(*,"(a,i4)") " 6 Disable labelling orbital indices"
			if (ilabelcomp==0) write(*,*) "7 Enable labelling orbital composition"
			if (ilabelcomp==1) write(*,"(a,i4)") " 7 Disable labelling orbital composition"
			write(*,"(' 8 Set label size, current:',i4)") labsize
			write(*,"(' 9 Set shifting of composition labels, current:',f6.3)") complabshift
			write(*,*) "10 Print molecular orbital energies"
			write(*,"(a,f9.4,a)") " 11 Set the criterion for determining degeneration, current:",degencrit," eV"
			write(*,*) "12 Set orbital energy shifting value"
			write(*,"(a,f6.2,a)") " 13 Set energy interval in the axis, current:",eneintv," eV"
            write(*,"(a,i3)") " 14 Set size of ticks and axis names, current:",ticknamesize
			read(*,*) isel2
			
			if (isel2==0) then
				exit
			else if (isel2==1.or.isel2==2) then
				if (nCDAfrag==2) then
					ifrag=1
					jfrag=2
				else
					write(*,*) "Plot the diagram for which two fragments?  e.g. 1,3"
					read(*,*) ifrag,jfrag
					if (ifrag>nCDAfrag.or.jfrag>nCDAfrag) then
						write(*,*) "Error: The index exceeded valid range!"
						cycle
					end if
				end if
				if (isel2==1) then
                    c80tmp="show"
				else if (isel2==2) then
                    c80tmp="save"
                    c80tmp2=graphformat
                    write(*,*) "Hint: Using pdf format is recommended for orbital interaction diagram"
                    write(*,*)
                    call setgraphformat
                end if
				if (iopshCDA==0.or.(iopshCDA==1.and.ispinplot==1)) call plotintdiag(trim(c80tmp),ifrag,jfrag,nCDAfrag,nmoCDA,&
				FOcomp,nmo,nmoCDA(ifrag),nmoCDA(jfrag),occCDA(:,0),occCDA(:,ifrag),occCDA(:,jfrag),&
				eneCDA(:,0)+eneshiftcomp,eneCDA(:,ifrag)+eneshiftA,eneCDA(:,jfrag)+eneshiftB,eneplotlow,eneplothigh,eneintv,conncritleft,conncritright,&
				idrawMObar,iconnlogi,ilabelorbidx,ilabelcomp,labsize,ticknamesize,complabshift,degencrit,eneshiftA,eneshiftB,eneshiftcomp)
				if (iopshCDA==1.and.ispinplot==2) call plotintdiag(trim(c80tmp),ifrag,jfrag,nCDAfrag,nmoCDA,&
				FOcompb,nmo,nmoCDA(ifrag),nmoCDA(jfrag),occCDAb(:,0),occCDAb(:,ifrag),occCDAb(:,jfrag),&
				eneCDAb(:,0)+eneshiftcomp,eneCDAb(:,ifrag)+eneshiftA,eneCDAb(:,jfrag)+eneshiftB,eneplotlow,eneplothigh,eneintv,conncritleft,conncritright,&
				idrawMObar,iconnlogi,ilabelorbidx,ilabelcomp,labsize,ticknamesize,complabshift,degencrit,eneshiftA,eneshiftB,eneshiftcomp)
				if (isel2==2) then
                    write(*,*) "Done! The graph has been saved to current folder with ""dislin"" prefix"
                    graphformat=c80tmp2
                end if
			else if (isel2==3) then
				write(*,*) "Input the lower and upper limits of the MO energy to be plotted (in eV)"
				write(*,*) "e.g. -70.6,8.5 (0,0 corresponds to the full energy range)"
				read(*,*) eneplotlowtmp,eneplothightmp
				if (eneplotlowtmp==0.and.eneplothightmp==0) then
					tmpval=maxval(eneCDA(:,0))-minval(eneCDA(:,0))
					eneplotlow=minval(eneCDA(:,0))-0.1D0*tmpval
					eneplothigh=maxval(eneCDA(:,0))+0.1D0*tmpval
				else if (eneplothightmp>eneplotlowtmp) then
					eneplothigh=eneplothightmp
					eneplotlow=eneplotlowtmp
				else
					write(*,*) "Error: Invalid input"
                    cycle
				end if
                eneintv=(eneplothigh-eneplotlow)/10
                degentmp=(eneplothigh-eneplotlow)/20
                if (degentmp<0.1D0) degencrit=degentmp
			else if (isel2==4) then
				write(*,*) "Input the criterion and rule for connecting orbital bars"
				write(*,"(a)") " Example 1: ""15,or,20"" means the criterion for connecting FO of fragment A (B) and complex MO is >=15% (>=20%)."
				write(*,"(a)") " Example 2: ""15,and,20"" means if a FO of fragment A and a FO of fragment B have contribution to a complex MO >=15% &
				&and >=20%, simultaneously and respectively, then the complex MO will be connected to these two FOs"
				write(*,"(a)") " Note: Inputting ""k"" can keep the current criterion unchanged. Using a criterion larger than 100 can nullify the connection."
				do while(.true.)
					read(*,"(a)") c80tmp
					if (index(c80tmp,'k')/=0) goto 40
					read(c80tmp,*,iostat=ierror) conncritleft,c80tmp2,conncritright
					if (ierror==0) exit
					write(*,*) "Error: Unrecognized input! Input again"
				end do
				if (c80tmp2=="or") iconnlogi=1
				if (c80tmp2=="and") iconnlogi=2
				conncritleft=conncritleft/100
				conncritright=conncritright/100
40				write(*,*) "How to draw the orbital bars?"
				write(*,*) "1 Draw bars for all FOs and all complex MOs"
				write(*,*) "2 Draw bars only for the FOs and the complex MOs connected by line"
				write(*,*) "3 Draw bars for all complex MOs and for the FOs connected by line"
				write(*,*) "4 Draw bars for all FOs and for the complex MOs connected by line"
				read(*,*) idrawMObar
			else if (isel2==5) then
				if (ispinplot==1) then
					ispinplot=2
				else if (ispinplot==2) then
					ispinplot=1
				end if
			else if (isel2==6) then
				if (ilabelorbidx==1) then
					ilabelorbidx=0
				else if (ilabelorbidx==0) then
					ilabelorbidx=1
				end if
			else if (isel2==7) then
				if (ilabelcomp==1) then
					ilabelcomp=0
				else if (ilabelcomp==0) then
					ilabelcomp=1
				end if
			else if (isel2==8) then
				write(*,*) "Input a value, e.g. 40"
				read(*,*) labsize
			else if (isel2==9) then
				write(*,*) "Input a value, e.g. 0.5"
				write(*,"(a)") " Note: The value should within 0 and 1. The more the value close to 0, the more the &
				&composition labels in the lines will close to the bars of complex orbitals"
				read(*,*) complabshift
			else if (isel2==10) then
				do ifrag=0,nCDAfrag
					if (ifrag==0) write(*,*) "Energy of molecular orbitals of the complex, in eV:"
					if (ifrag>0) write(*,"(a,i3,a)") " Energy of molecular orbitals of fragment",ifrag,", in eV:"
					if (iopshCDA==0) then
						write(*,"(7f11.4)") eneCDA(1:nmoCDA(ifrag),ifrag)
					else if (iopshCDA==1) then
						write(*,*) "                           **** Alpha orbitals ****"
						write(*,"(7f11.4)") eneCDA(1:nmoCDA(ifrag),ifrag)
						write(*,*) "                           **** Beta orbitals ****"
						write(*,"(7f11.4)") eneCDAb(1:nmoCDA(ifrag),ifrag)
					end if
					write(*,*)
				end do
			else if (isel2==11) then
				write(*,*) "Input the criterion in eV, e.g. 0.5"
				read(*,*) degencrit
			else if (isel2==12) then
				write(*,"(a)") " Note: Fragment A and B means the fragment at the left and right side of the orbital interaction diagram, respectively"
				do while(.true.)
					write(*,*) "0 Return"
					write(*,"(a,f10.5)") " 1 Set energy shifting value for all orbitals of fragment A, current:",eneshiftA
					write(*,"(a,f10.5)") " 2 Set energy shifting value for all orbitals of fragment B, current:",eneshiftB
					write(*,"(a,f10.5)") " 3 Set energy shifting value for all orbitals of complex, current:",eneshiftcomp
					read(*,*) itmp
					if (itmp==0) then
						exit
					else
						write(*,*) "Input the value in eV, e.g. 2.4"
						if (itmp==1) read(*,*) eneshiftA
						if (itmp==2) read(*,*) eneshiftB
						if (itmp==3) read(*,*) eneshiftcomp
					end if
				end do
			else if (isel2==13) then
				write(*,*) "Input energy interval in eV, e.g. 1.5"
				read(*,*) eneintv
            else if (isel2==14) then
				write(*,*) "Input the text size, e.g. 50"
                read(*,*) ticknamesize
			end if
		end do
			
	else if (isel==6) then !Output detail of a complex orbital contribution to d,b,r
		if (nCDAfrag==2) then
			ifrag=1
			jfrag=2
		else
			write(*,*) "Perform the analysis between which two fragments?  e.g. 1,3"
			read(*,*) ifrag,jfrag
			if (ifrag>nCDAfrag.or.jfrag>nCDAfrag) then
				write(*,*) "Error: The index exceeded valid range!"
				cycle
			end if
		end if
		do while(.true.)
			write(*,*) "Input the index of complex orbital you are interested in, e.g. 6"
			write(*,*) "Input 0 can return"
			read(*,*) iorb
			if (iorb<0.or.iorb>nmo) then
				write(*,"(a,i7,a,i7,/)") " Error: Orbital index should between",1," and",nmo
				cycle
			else if (iorb==0) then
				exit
			end if
			write(*,"(a)") " Set output threshold, e.g. 0.001. If contribution of a pair of fragment orbitals to any of d,b,r is larger than this value then it will be printed."
			read(*,*) thres
			if (iout==10) open(iout,file="CDA.txt",status="replace")
			write(*,"(' Occupation number of orbital',i6,' of the complex:',f12.8)") iorb,occCDA(iorb,0)
			if (iopshCDA==0) refocc=2D0 !Reference orbital occupation number
			if (iopshCDA==1) refocc=1D0
			if (iopshCDA==1) write(iout,"(/,a)") "                     **** Result for alpha electrons ****"
			write(iout,*) "FragA Orb(Occ.)  FragB Orb(Occ.)      d           b        d - b          r"
            sumd=0
            sumb=0
            sumr=0
			do iAidx=1,nmoCDA(ifrag) !Scan each FO pairs
				iAtmp=0
				if (ifrag>1) iAtmp=sum(nmoCDA(1:ifrag-1)) !Actual index of FO of the first fragment in all FOs
				iA=iAtmp+iAidx
				do iBidx=1,nmoCDA(jfrag)
					iBtmp=0
					if (jfrag>1) iBtmp=sum(nmoCDA(1:jfrag-1))
					iB=iBtmp+iBidx
					dtmp=0
					btmp=0
					!occfac is defined by Tian Lu. =0 means both are virtual or occupied, or have identical occ.   >0 means donor->acceptor   <0 means donor<-acceptor
					!In particular, For SCF wavefunction, =1 occ_A->vir_B,  =-1 vir_A<-occ_B
					occfac=(occCDA(iAidx,ifrag)-occCDA(iBidx,jfrag))/refocc
					tmpval=occCDA(iorb,0)*occfac*coFO(iA,iorb)*coFO(iB,iorb)*FOovlpmat(iA,iB)
					if (occfac>0) dtmp=tmpval
					if (occfac<0) btmp=-tmpval  !The minus sign before tmpval cancelled the negative of occfac
					!The overlap population between the two FOs in current complex orbital, is the twice of original r
					rtmp=occCDA(iorb,0)*2*min(occCDA(iAidx,ifrag),occCDA(iBidx,jfrag))/refocc*coFO(iA,iorb)*coFO(iB,iorb)*FOovlpmat(iA,iB)
					if (abs(dtmp)>thres.or.abs(btmp)>thres.or.abs(rtmp)>thres) then
						write(iout,"(i5,'(',f7.4,')',3x,i5,'(',f7.4,')',4f12.6)") iAidx,occCDA(iAidx,ifrag),iBidx,occCDA(iBidx,jfrag),dtmp,btmp,dtmp-btmp,rtmp
                        sumd=sumd+dtmp
                        sumb=sumb+btmp
                        sumr=sumr+rtmp
                    end if
				end do
			end do
            write(iout,"(' Sum of above terms:',11x,4f12.6)") sumd,sumb,sumd-sumb,sumr
			if (iopshCDA==1) then
				write(iout,"(/,a)") "                     **** Result for beta electrons ****"
				write(iout,*) "FragA Orb(Occ.)  FragB Orb(Occ.)      d           b        d - b          r"
				sumd=0
				sumb=0
				sumr=0
				do iAidx=1,nmoCDA(ifrag) !Scan each FO pairs
					iAtmp=0
					if (ifrag>1) iAtmp=sum(nmoCDA(1:ifrag-1)) !Actual index of FO of the first fragment in all FOs
					iA=iAtmp+iAidx
					do iBidx=1,nmoCDA(jfrag)
						iBtmp=0
						if (jfrag>1) iBtmp=sum(nmoCDA(1:jfrag-1))
						iB=iBtmp+iBidx
						dtmp=0
						btmp=0
						occfac=(occCDAb(iAidx,ifrag)-occCDAb(iBidx,jfrag))/refocc
						tmpval=occCDAb(iorb,0)*occfac*coFOb(iA,iorb)*coFOb(iB,iorb)*FOovlpmatb(iA,iB)
						if (occfac>0) dtmp=tmpval
						if (occfac<0) btmp=-tmpval
						rtmp=occCDAb(iorb,0)*2*min(occCDAb(iAidx,ifrag),occCDAb(iBidx,jfrag))/refocc*coFOb(iA,iorb)*coFOb(iB,iorb)*FOovlpmatb(iA,iB)
						if (abs(dtmp)>thres.or.abs(btmp)>thres.or.abs(rtmp)>thres) then
							write(iout,"(i5,'(',f7.4,')',3x,i5,'(',f7.4,')',4f12.6)") iAidx,occCDAb(iAidx,ifrag),iBidx,occCDAb(iBidx,jfrag),dtmp,btmp,dtmp-btmp,rtmp
							sumd=sumd+dtmp
							sumb=sumb+btmp
							sumr=sumr+rtmp
                        end if
					end do
				end do
				write(iout,"(' Sum of above terms:',11x,4f12.6)") sumd,sumb,sumd-sumb,sumr
			end if
			write(*,*)
			if (iout==10) then
				close(10)
				write(*,*) "Done! Information has been output to CDA.txt in current folder"
				write(*,*)
			end if
		end do

	end if !End post-processing interface
end do
end subroutine


!--Read occupation number of natural orbitals from pop=full output. Assume that nbasis=nmo. Must use loclabel to locate to title line first
!ispace=0 used to read natural orbital occupation numbers; ispace=1 is used to read MO energies
subroutine readgauorbeig(eigvec,ndim,ispace)
character c80tmp*80
integer ipos
real*8 eigvec(ndim)
read(10,*)
read(10,*)
if (ispace==1) read(10,*)
read(10,"(a21)",advance='no') c80tmp
ipos=0
do while(.true.)
	ipos=ipos+1
	if (ipos>ndim) exit
	read(10,"(f10.5)",advance='no') eigvec(ipos)
	if (mod(ipos,5)==0) then
		do itmp=1,ndim+2+ispace
			read(10,*)
		end do
		read(10,"(a21)",advance='no') c80tmp
	end if
end do
backspace(10)
backspace(10)
if (ispace==1) backspace(10)
end subroutine


!----------- Plot orbital interaction diagram
!ilabelorbidx=0 means don't label orbital indices, =1 means label them
!ilabelcomp=0 means don't label orbital composition, =1 means label them
subroutine plotintdiag(status,ifrag,jfrag,nCDAfrag,nmoCDA,FOcp,nmo0,nmo1,nmo2,occ0,occ1,occ2,ene0,ene1,ene2,&
eneplotlow,eneplothigh,eneintv,conncritleft,conncritright,idrawMObar,iconnlogi,ilabelorbidx,ilabelcomp,labsize,ticknamesize,&
complabshift,degencrit,eneshiftA,eneshiftB,eneshiftcomp)
use dislin
use defvar
implicit real*8 (a-h,o-z)
integer nmo0,nmo1,nmo2,labsize,ticknamesize,idrawMObar,iconnlogi,ilabelorbidx,ilabelcomp,nmoCDA(0:nCDAfrag)
real*8 FOcp(nmo0,nmo0),occ0(nmo0),occ1(nmo1),occ2(nmo2),ene0(nmo0),ene1(nmo1),ene2(nmo2)
real*8 eneplotlow,eneplothigh,conncritleft,conncritright
real*8 movetextx,movetexty,complabshift,degencrit,eneshiftA,eneshiftB,eneshiftcomp
logical FO1involveconn(nmo1),FO2involveconn(nmo2)
character c80tmp*80,status*4
barsep=0.10D0
xlow1=1
xhigh1=2
xlow0=4
xhigh0=5.5D0
xlow2=7.5D0
xhigh2=8.5D0
xaxislow=xlow1-1
xaxishigh=xhigh2+1
movetextx=labsize/2D0*(xaxishigh-xaxislow)/1500
movetexty=labsize/2D0*(eneplothigh-eneplotlow)/2300
shiftdegen=0.70D0

call SCRMOD('REVERSE')
CALL PAGE(2700,2700)
CALL IMGFMT("RGB")
CALL setxid(0,'NONE') !If we don't set this, after we draw a graph embedded in GUI(e.g. relif map), curve map will not be shown 
if (status=="show") then
	call METAFL('xwin')
	call window(200,100,720,720)
else if (status=="save") then
	call METAFL(graphformat)
	call winsiz(1500,1500)
end if
CALL DISINI
if (status=="show") then
	call WINTIT("Click right mouse button to close...")
	call hwfont
else if (status=="save") then
	if (ttfontfile=="none") then
		CALL HELVES
    else
		CALL TTFONT(ttfontfile)
    end if
	CALL SHDCHA
end if
call center
call AXSLEN(1500,2300)
call height(ticknamesize)
call hname(ticknamesize)
write(c80tmp,"(a,i3)") "Frag.",ifrag
call messag(c80tmp,870,2600)
call messag("Complex",1330,2600)
write(c80tmp,"(a,i3)") "Frag.",jfrag
call messag(c80tmp,1900,2600)
CALL NAME('Orbital energy (eV)','Y')
CALL LABDIG(2,"Y")
CALL NAMDIS(60,'Y')
call ERRMOD("ALL","OFF")
call frame(0)
call setgrf('none','name','none','ticks')
CALL GRAF(xaxislow,xaxishigh,xaxislow,xaxishigh, eneplotlow,eneplothigh,eneplotlow,eneintv)

!X range: xaxislow~xaxishigh,  FO1 bar: xlow1~xhigh1  MO bar: xlow0~xhigh0  FO2 bar: xlow2~xhigh2
FO1involveconn=.false. !First assume all FOs are not involved in connecting
FO2involveconn=.false.

!Plot MOs of complex
degenmovex=0
do iorb0=1,nmo0
	involveconn=0 !First assume this complex MO is not involved in connecting
	eneval=ene0(iorb0)
	if (eneval<eneplotlow.or.eneval>eneplothigh) cycle
	call solid
	call height(nint(labsize*0.8D0)) !Use smaller size of text to label composition
	!Plot lines between complex and fragment A
	itmp=0
	if (ifrag>1) itmp=sum(nmoCDA(1:ifrag-1))
	jtmp=0
	if (jfrag>1) jtmp=sum(nmoCDA(1:jfrag-1))
	do iorb1=1,nmo1
		iorb1act=itmp+iorb1 !Actual orbital index in all FOs
		!Link orbital bars
		if (FOcp(iorb1act,iorb0)<conncritleft.or.ene1(iorb1)<eneplotlow.or.ene1(iorb1)>eneplothigh) cycle
		if (iconnlogi==2) then !"and" rule, so we must test if there is a FO of fragment B having enough large contribution
			do iorb2=1,nmo2
				eneval2=ene2(iorb2)
				if (eneval2<eneplotlow.or.eneval2>eneplothigh) cycle
				if (FOcp(jtmp+iorb2,iorb0)>=conncritright) exit
			end do
			if (iorb2==nmo2+1) cycle !No FO of fragment B has large enough contribution
		end if
		involveconn=1
		FO1involveconn(iorb1)=.true.
		CALL SETRGB(1D0,0D0,0D0)
		call rline(xlow0-barsep,eneval,xhigh1+barsep,ene1(iorb1))
		!Label orbital composition
		if (ilabelcomp==1) then
			CALL SETRGB(0.6D0,0D0,0D0)
			write(c80tmp,"(i3)") nint(FOcp(iorb1act,iorb0)*100)
			call rlmess(trim(adjustl(c80tmp))//'%',xlow0*(1-complabshift)+xhigh1*complabshift-movetextx,eneval*(1-complabshift)+ene1(iorb1)*complabshift+movetexty)
		end if
	end do
	!Plot lines between complex and fragment B
	do iorb2=1,nmo2
		iorb2act=jtmp+iorb2 !Actual orbital index in all FOs
		if (FOcp(iorb2act,iorb0)<conncritright.or.ene2(iorb2)<eneplotlow.or.ene2(iorb2)>eneplothigh) cycle
		if (iconnlogi==2) then !"and" rule, so we must test if there is a FO of fragment A having enough large contribution
			do iorb1=1,nmo1
				eneval1=ene1(iorb1)
				if (eneval1<eneplotlow.or.eneval1>eneplothigh) cycle
				if (FOcp(itmp+iorb1,iorb0)>=conncritleft) exit
			end do
			if (iorb1==nmo1+1) cycle !No FO of fragment B has large enough contribution
		end if
		involveconn=1
		FO2involveconn(iorb2)=.true.
		CALL SETRGB(1D0,0D0,0D0)
		call rline(xhigh0+barsep,eneval,xlow2-barsep,ene2(iorb2))
		if (ilabelcomp==1) then
			CALL SETRGB(0.6D0,0D0,0D0)
			write(c80tmp,"(i3)") nint(FOcp(iorb2act,iorb0)*100) 
			call rlmess(trim(adjustl(c80tmp))//'%',xlow2*complabshift+xhigh0*(1-complabshift)-movetextx,eneval*(1-complabshift)+ene2(iorb2)*complabshift+movetexty)
		end if
	end do
	
	!Plot orbital bar
	if ((idrawMObar==2.or.idrawMObar==4).and.involveconn==0) then
		ilastplotted=0
		cycle
	end if
	call height(labsize)
	CALL SETRGB(0D0,0D0,0D0)
	call solid
	if (occ0(iorb0)==0) call dash
	call rline(xlow0,eneval,xhigh0,eneval)
	if (ilabelorbidx==1) then
		CALL SETRGB(0D0,0D0,1D0)
		write(c80tmp,*) iorb0
		if (iorb0>=2) then !When near degenerate, several labels will superposition, so slightly shift text by degenmovex to avoid this problem
			if (abs(eneval-ene0(iorb0-1))<degencrit.and.ilastplotted==1) then
				degenmovex=degenmovex+shiftdegen
			else
				degenmovex=0D0
			end if
		end if
		call rlmess(trim(adjustl(c80tmp)),(xlow0+xhigh0)/2-0.3D0+degenmovex-movetextx,eneval+movetexty)
	end if
	ilastplotted=1 !This time the label is plotted
end do

!Plot FOs of fragment 1
call height(labsize)
degenmovex=0
do iorb=1,nmo1
	CALL SETRGB(0D0,0D0,0D0)
	eneval=ene1(iorb)
	if (eneval<eneplotlow.or.eneval>eneplothigh) then
		ilastplotted=0
		cycle
	end if
	if ((idrawMObar==2.or.idrawMObar==3).and.(.not.FO1involveconn(iorb))) cycle
	call solid !Use solid line to plot occupied orbital bars, use dashed line to plot virtual orbital bars
	if (occ1(iorb)==0) call dash
	call rline(xlow1,eneval,xhigh1,eneval)
	if (ilabelorbidx==1) then
		CALL SETRGB(0D0,0D0,1D0)
		write(c80tmp,*) iorb
		if (iorb>=2) then !When near degenerate, several labels will superposition, so slightly shift text by degenmovex to avoid this circumstance
			if (abs(eneval-ene1(iorb-1))<degencrit.and.ilastplotted==1) then
				degenmovex=degenmovex+shiftdegen
			else
				degenmovex=0D0
			end if
		end if
		call rlmess(trim(adjustl(c80tmp)),xlow1-0.2D0+degenmovex-movetextx,eneval+movetexty) !movetext*pix2usrx/y is used to shift text so that its center can occur at expected position
	end if
	ilastplotted=1 !This time the label is plotted
end do

!Plot FOs of fragment 2
degenmovex=0
do iorb=1,nmo2
	CALL SETRGB(0D0,0D0,0D0)
	eneval=ene2(iorb)
	if (eneval<eneplotlow.or.eneval>eneplothigh) then
		ilastplotted=0
		cycle
	end if
	if ((idrawMObar==2.or.idrawMObar==3).and.(.not.FO2involveconn(iorb))) cycle
	call solid
	if (occ2(iorb)==0) call dash
	call rline(xlow2,eneval,xhigh2,eneval)
	if (ilabelorbidx==1) then
		CALL SETRGB(0D0,0D0,1D0)
		write(c80tmp,*) iorb
		if (iorb>=2) then !When near degenerate, several labels will superposition, so slightly shift text by degenmovex to avoid this circumstance
			if (abs(eneval-ene2(iorb-1))<degencrit.and.ilastplotted==1) then
				degenmovex=degenmovex+shiftdegen
			else
				degenmovex=0D0
			end if
		end if
		call rlmess(trim(adjustl(c80tmp)),xlow2+0.3D0+degenmovex-movetextx,eneval+movetexty)
	end if
	ilastplotted=1 !This time the label is plotted
end do

CALL DISFIN
end subroutine