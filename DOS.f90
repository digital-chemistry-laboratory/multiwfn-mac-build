!!----------------- Plot various kinds of DOS map
!For .out or plain text file, only one type of spin MOs will be loaded and processed, and then we will not consider spin type
!For .mwfn/.fch/.molden/.gms, etc., user can switch spin type anytime

!  Note on arrays and open-shell case:
!  MOene_dos and MOocc_dos are working horse for present module, they are completely identical to MOene and MOocc but may be changed in this module (e.g. between a.u. and eV).
!Both of them have size of "nmo", the first and last half correspond to alpha and beta.
!  Other arrays such as str, FWHM, PDOSliney and compfrag also have size of nmo. For closed shell case, imoend=nmo, and all slots are meaningful.
!When only one spin is taken into account, namely alpha (ispin=1) or beta (ispin=2) case, imoend=nbasis, and their (1:nbasis) part correspond
!to alpha or beta MOs, respectively, and (nbasis+1:nmo) part is meaningless.
!  Therefore, to cycle all MOs in current set, looping "imo=1,imoend", the irealmo=imo for ispin=0/1 and irealmo=imo+nbasis for ispin==2 should be used for MOene_dos and MOocc_dos,
!while imo should be used for other arrays

subroutine DOS
use defvar
use util
use dislin
use functions
implicit real*8 (a-h,o-z)
integer,parameter :: nfragmax=10
integer,parameter :: num2Dpoints=200 !The number of points constituting the X-axis of 2D LDOS
integer :: iCOHP=0 !If=1, plot COHP instead of various DOS. In this case idoPDOS=0 and idoOPDOS=0
!Temporary variables
real*8,pointer :: tmpmat(:,:),tmpmat2(:,:)
character unitstr*5,c80tmp*80,c200tmp*200,c2000tmp*2000,selectyn
real*8 atmcomp(ncenter,nmo) !Use to store all atomic contributions to MOs by Hirshfeld/Becke methods
!Curve data
real*8 :: curvexpos(num1Dpoints),TDOScurve(num1Dpoints),OPDOScurve(num1Dpoints),PDOScurve(num1Dpoints,nfragmax),LDOScurve(num1Dpoints),COHPcurve(num1Dpoints)
real*8 :: LDOSxpos(num2Dpoints)
!All ?DOSliney share linexpos(:) as X axis data
real*8,allocatable :: linexpos(:),TDOSliney(:),TDOSliney_unocc(:),PDOSliney(:,:),OPDOSliney(:),LDOSliney(:),COHPliney(:),COHPliney_unocc(:) !TDOSliney_unocc only records TDOS of unoccupied MOs
real*8,allocatable :: compfrag(:,:) !i,k element is the composition of fragment k in MO i. compfrag(:,0) is used for recording degeneracy for TDOS
real*8,allocatable :: OPfrag12(:) !Overlap population between fragment 1 and 2
real*8,allocatable :: LDOScomp(:) !Composition at a point of each orbital
real*8,allocatable :: LDOSptscomp(:,:) !Composition of each MO, ipt in a given line
real*8,allocatable :: LDOS2Dmap(:,:) !LDOS curve, ipt in a given line
real*8 HOMOlevx(2),HOMOlevy(2)
real*8,allocatable :: MOene_dos(:),MOocc_dos(:) !Using the ene/occ in this to plot DOS, the values are scaled when changing between a.u. and eV. The original MOene is remain unchanged
integer,allocatable :: selorbarr(:)
character clegend*960 !Legend strings. (10+2) lines * 80 character per line
integer :: legendx=400,legendy=160
integer :: ticksize=45,height_axis=45,legtextsize=42
character :: TDOSstring*80="TDOS",OPDOSstring*80="OPDOS",COHPstring*80="COHP",graphformat_old*4
character :: PDOSstring(nfragmax)*80=(/ character(len=80) :: "PDOS frag.1","PDOS frag.2","PDOS frag.3","PDOS frag.4","PDOS frag.5","PDOS frag.6","PDOS frag.7","PDOS frag.8","PDOS frag.9","PDOS frag.10"/)
integer :: ishowPDOSline(nfragmax),ishowPDOScurve(nfragmax),icurvewidth=4,ilinewidth=2,intarr(2)
integer :: iclrPDOS(nfragmax)=(/ 1,3,10,14,12,9,13,11,6,7 /)
!Below are used for defining fragments. For Mulliken/SCPA, they correspond to basis function, while for Hirshfeld/Becke, they correspond to atom indices
integer :: nfragDOS(nfragmax) !The number of terms in the fragments (0=undefined)
integer,allocatable :: fragDOS(:,:) !The index of basis functions/atoms or MOs (#1) in fragments (#2)
!Below are used for photoelectron spectra (PES). In order to avoid confusing, they are not share with the ones for DOS
real*8 :: PES_shift=0,PES_Xlow=0,PES_Xhigh=5,PES_Xstep=0.5D0,scalePEScurve=0.2D0
real*8,allocatable :: PESlinex(:),PESliney(:),PES_str(:),PES_FWHM(:)
real*8 :: PEScurve(num1Dpoints)
real*8,allocatable :: bindene(:)
integer :: ishowPEScurve=1,ishowPESline=1,iusersetPES_Y=0,invPES_X=0,ilinebottom=0

if (.not.(ifiletype==0.or.allocated(CObasa))) then
	write(*,"(a,/)") " Error: This function is only available for input file containing basis function information &
	&(i.e. .mwfn/.fch/.molden/.gms) and plain text file with energy levels!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

!Initialize parameters and some arrays
if (allocated(FWHM)) deallocate(FWHM) !Global array
if (allocated(str)) deallocate(str) !Global array
iCOHP=0
defFWHM=0.05D0 !Default FWHM
icompmethod=1 !The method for calculate PDOS, =1 Mulliken, =2 SCPA, =3 Hirshfeld, =4 Becke
ibroadfunc=2 !Default is Gaussian function
scalecurve=0.1D0 !Multiply curves with this value
enelow=-0.8D0 !Energy range, a.u.
enehigh=0.2D0
eneshift=0
stepx=0.1D0
stepyleft=2
gauweigh=0.5D0 !The weight of Gaussian in Pseudo-Voigt function
iPDOStype=0 !=0 Undefined.  =1: The PDOS to be plotted is for basis function/atom.  =2: MO-PDOS
nfragDOS=0
ishowTDOScurve=1
ishowTDOSline=1
ishowPDOSline=0
ishowPDOScurve=0
ishowOPDOScurve=0
ishowOPDOSline=0
ishowCOHPline=0
ishowlegend=1
ishowHOMOlev=1
ishowvertline=0
ishowYlab=1
nlabdigX=2
nlabdigY=2
nlabdigY_OPDOS=2
nlabdigY_LDOS=3
nlabdigleftY_COHP=2
nlabdigrightY_COHP=2
iunitx=1
unitstr=" a.u."
idegen=0
!ispin=0: Restricted wavefunction, >0: Unrestricted wavefunction (=1: alpha, =2: beta, =3: mix alpha+beta)
ispin=0
if (wfntype==1.or.wfntype==4) ispin=1 !For U, alpha part by default
iusersetY=0 !If user has set lower and upper range of Y axis by himself
iusersetcolorscale=0 !If user has set color scale of 2D LDOS by himself
Yrightsclfac=0.5D0 !Scale factor relative to left Y-axis of OPDOS (right Y-axis)
yxratio=1D0
graphformat_old=graphformat !User may change graphformat, backup it
graphformat="pdf"
call setfil("dislin."//trim(graphformat))

ireadgautype=1
if (ifiletype==0) then
	!Read energy level information from text file, the first number in first row define how many energy levels
	!in there, the second number in first row if equals to 1, means below data are only energies, if equals to 2,
	!means both strength and FWHM also present.
	open(10,file=filename,status="old")
	call loclabel(10,"Gaussian, Inc",igauout,maxline=100)
	rewind(10)
	if (igauout==1) then
		write(*,*) "This is Gaussian output file"
		if (ireadgautype==1) then !Read energy level from Gaussian output
			call loclabel(10,"NBsUse=")
			read(10,*) c200tmp,nbasis
			nmo=nbasis
			allocate(MOene(nmo),MOocc(nmo),str(nmo),FWHM(nmo))
			call loclabel(10,"Orbital energies and kinetic energies",ifound) !First assume this is closed-shell
			if (ifound==1) then
				read(10,*)
				read(10,"(a)") c200tmp
				do i=1,nbasis
					read(10,"(a21)",advance="no") c200tmp
					read(10,*) MOene(i)
					MOocc(i)=0
					if (index(c200tmp,'O')/=0) MOocc(i)=2
				end do
				call loclabel(10,"Orbital energies and kinetic energies (beta)",ifound)
				if (ifound==1) then
					where(MOocc==2) MOocc=1
					write(*,*) "Read which type of orbitals? 1:alpha 2:beta"
					read(*,*) inp
					if (inp==2) then !Read beta energies, overlay read alpha counterpart
						read(10,*)
						read(10,*)
						do i=1,nbasis
							read(10,"(a21)",advance="no") c200tmp
							read(10,*) MOene(i)
							MOocc(i)=0
							if (index(c200tmp,'O')/=0) MOocc(i)=1
						end do
					end if
				end if			
				write(*,*) "Read orbital energy from the file"
			else
				write(*,"(a)") " Error: Cannot find orbital energies from this file, don't forget using pop=full keyword"
				write(*,*)
				return
			end if
		end if
		str=1D0
		FWHM=defFWHM
	else !Plain text file
		read(10,*) nmo,inp
		allocate(MOene(nmo),MOocc(nmo),str(nmo),FWHM(nmo))
		if (inp==1.or.inp==3) then
			do imo=1,nmo
				read(10,*) MOene(imo),MOocc(imo)
			end do
			str=1D0
			FWHM=defFWHM
		else if (inp==2.or.inp==4) then
			do imo=1,nmo
				read(10,*) MOene(imo),MOocc(imo),str(imo),FWHM(imo)
			end do
		end if
        if (inp==3.or.inp==4) then !Input file contains energies in eV, change unit
            MOene=MOene/au2eV !MOene only records a.u. Unit will be changed to eV when assigning MOene_dos
            if (inp==4) FWHM=FWHM/au2eV
            iunitx=2
            enelow=enelow*au2eV
		    enehigh=enehigh*au2eV
		    unitstr=" eV"
		    scalecurve=scalecurve*au2eV
		    stepx=stepx*au2eV
		    stepyleft=stepyleft/au2eV
        end if
	end if
	close(10)
else if (allocated(CObasa)) then !For ispin=1 or 2, only 1:nbasis is used, while for ispin==3 (both spin), all nmo slots will be used
	allocate(str(nmo),FWHM(nmo))
	str=1D0
	FWHM=defFWHM
else
	write(*,*) "Error: Your input file does not contain basis function information!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

if (all(MOene==0)) then
	write(*,*) "Error: All orbitals have zero energy! In this case DOS cannot be plotted!"
    if (ifiletype==9) write(*,"(a)") " If the molden file was produced by CP2K, note that OT should not be used, which does not produce orbital energies"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

!Allocate all arrays that may be used, do not consider if they will actually be used, because memory consuming is very little
allocate(linexpos(3*nmo),TDOSliney(3*nmo),TDOSliney_unocc(3*nmo),PDOSliney(3*nmo,nfragmax),OPDOSliney(3*nmo),LDOSliney(3*nmo),COHPliney(3*nmo),COHPliney_unocc(3*nmo))
allocate(compfrag(nmo,0:nfragmax),OPfrag12(nmo))
allocate(fragDOS(nbasis,nfragmax+1)) !The last slot is used to exchange fragment
allocate(LDOScomp(nmo))

!Set from where to where are active energy levels
if (ispin==0.or.ispin==3) imoend=nmo !Text file or restricted .fch, or unrestricted but consider both spins
if (ispin==1.or.ispin==2) imoend=nbasis !For unrestricted fch or Gaussian output file while only consider alpha or beta

!MOene_dos is the working horse, recording energies in actually employed unit
!MOocc_dos records occupation numbers for DOS plotting, which may be modified by users in present function
if (allocated(MOene_dos)) deallocate(MOene_dos,MOocc_dos)
allocate(MOene_dos(nmo),MOocc_dos(nmo))
MOene_dos=MOene
MOocc_dos=MOocc
if (ifiletype==0.and.(inp==3.or.inp==4)) MOene_dos=au2eV*MOene

!Suggested setting for periodic wavefunction, which may contain numerous densely distributed orbitals
if (ifPBC>0) then
    ilinebottom=1
    nlabdigX=1
	nlabdigY=1
	nlabdigY_OPDOS=1
    !Use eV unit
	iunitx=2
	unitstr=" eV"
    FWHM=0.5D0 !Usually suitable
	MOene_dos=MOene_dos*au2eV
	str=str/au2eV
	scalecurve=scalecurve*au2eV
	stepx=1
	stepyleft=1
	eneHOMO=-1D99 !Determine HOMO energy
	do imo=1,imoend
		irealmo=imo
		if (ispin==2) irealmo=imo+nbasis
		if (MOocc_dos(irealmo)>0) then
			if (MOene_dos(irealmo)>eneHOMO) eneHOMO=MOene_dos(irealmo)
        else
			eneLUMO=MOene_dos(irealmo)
            exit
        end if
	end do
    enelow=eneHOMO-8
    enehigh=eneLUMO+5
    icompmethod=2 !Use SCPA for PDOS by default, which doesn't rely on Sbas; evaluation of Sbas is quite time-consuming for large cell
end if



!!!!! ***** Main loop ***** !!!!!!
!!!!! ***** Main loop ***** !!!!!!
!!!!! ***** Main loop ***** !!!!!!
do while(.true.)

idoPDOS=0
idoOPDOS=0
if (iCOHP==0) then
	if (any(nfragDOS>0)) idoPDOS=1
	if (all(nfragDOS(1:2)>0).and.icompmethod<=2.and.iPDOStype==1) idoOPDOS=1
end if

!Unknow text file does not contain wavefunction info, couldn't define fragment
write(*,*)
write(*,"(a)") " Hint: You can input ""s"" to save current plotting status to a file, or input ""l"" to load status from a file"
write(*,*)
if (iCOHP==0) then
	call menutitle("Plot density-of-states (DOS)",10,2)
else
	call menutitle("Plot COHP",10,2)
end if
write(*,*) "-10 Return to main menu"
if (iCOHP==0) write(*,*) "-7 Change to COHP plotting mode"
write(*,"(a,f10.4,1x,a)") " -6 Set shift of energy levels, current:",eneshift,trim(unitstr)
write(*,*) "-5 Customize energy levels, occupations, strengths and FWHMs for specific MOs"
write(*,*) "-4 Show all orbital information"
write(*,*) "-3 Export energy levels, occupations, strengths and FWHMs to plain text file"
if (iCOHP==0) then
	write(*,*) "-2 Define MO fragments for MO-PDOS"
	if (allocated(CObasa)) write(*,*) "-1 Define fragments for PDOS/OPDOS"
	if (idoOPDOS==1) then
		write(*,*) "0 Draw TDOS+PDOS+OPDOS graph!    -0 Draw TDOS+PDOS!"
	else if (idoPDOS==1) then
		if (iPDOStype==1) then
			write(*,*) "0 Draw TDOS+PDOS graph!"
		else if (iPDOStype==2) then
			write(*,*) "0 Draw TDOS with MO-PDOS graph!"
		end if
	else
		write(*,*) "0 Draw TDOS graph!" !Reading text file can only draw spinless TDOS, because they impossible to define fragment
	end if
    if (nfragDOS(1)==0.and.nfragDOS(2)==0) write(*,*) "00 Draw TDOS and OPDOS between nearest atoms!"
else
	write(*,*) "-1 Define fragments for plotting COHP between fragments"
    if (nfragDOS(1)==0.and.nfragDOS(2)==0) then
		write(*,*) "0 Draw COHP between nearest atoms!"
    else
		write(*,*) "0 Draw COHP between the two defined fragments!"
    end if
end if
if (ibroadfunc==1) write(*,*) "1 Select broadening function, current: Lorentzian"
if (ibroadfunc==2) write(*,*) "1 Select broadening function, current: Gaussian"
if (ibroadfunc==3) write(*,*) "1 Select broadening function, current: Pseudo-Voigt"
write(*,"(a,f11.3,a,f11.3,a,a,f7.3)") " 2 Set energy range and step, current:",enelow," to",enehigh,trim(unitstr),", step:",stepx
if (maxval(FWHM)==minval(FWHM)) then
	write(*,"(a,f10.5,a)") " 3 Set full width at half maximum (FWHM), current:",FWHM(1),unitstr
else
	write(*,"(a,f10.5)") " 3 Set full width at half maximum (FWHM), current: Orbital dependent"
end if
if (iCOHP==0) then
	write(*,"(a,f10.5)") " 4 Set scale ratio for DOS curve, current:",scalecurve
else
	write(*,"(a,f10.5)") " 4 Set scale ratio for COHP curve, current:",scalecurve
end if
if (ibroadfunc==3) write(*,"(a,f10.5)") " 5 Set Gaussian-weighting coefficient, current:",gauweigh
if (ispin==1) write(*,*) "6 Choose orbital spin, current: Alpha"
if (ispin==2) write(*,*) "6 Choose orbital spin, current: Beta"
if (ispin==3) write(*,*) "6 Choose orbital spin, current: Both"
if (iCOHP==0) then
	if (allocated(CObasa)) then
		if (icompmethod==1) write(*,*) "7 Set the method for calculating PDOS, current: Mulliken"
		if (icompmethod==2) write(*,*) "7 Set the method for calculating PDOS, current: SCPA"
		if (icompmethod==3) write(*,*) "7 Set the method for calculating PDOS, current: Hirshfeld"
		if (icompmethod==4) write(*,*) "7 Set the method for calculating PDOS, current: Becke"
	end if
end if
if (iCOHP==0) write(*,*) "8 Switch unit between a.u. and eV, current:"//unitstr !Do now allow user to change to a.u. in COHP case to avoid cumbersome things
if (idegen==1) write(*,"(a,f6.3,' eV')") " 9 Toggle considering degenerate, current: Yes, with threshold of ",degencrit
if (idegen==0.and.iCOHP==0) write(*,"(a)") " 9 Toggle using line height to show orbital degeneracy, current: No"
if (iCOHP==0) then
	write(*,*) "10 Draw local DOS for a point"
	write(*,*) "11 Draw local DOS along a line"
	write(*,*) "12 Enter interface for plotting photoelectron spectrum (PES)"
end if

read(*,*) c80tmp

if (index(c80tmp,'s')/=0) then
    write(*,"(a)") " Input file path for saving current status, e.g. C:\Bang_Dream\RAS.dat"
    write(*,"(a)") " Note: If you press ENTER button directly, status will be saved to DOS.dat in current folder"
    read(*,"(a)") c200tmp
    if (c200tmp==" ") c200tmp="DOS.dat"
    open(10,file=c200tmp,status="replace")
    write(10,*) nmo !Used to determine consistency between the status file and the current system
    !Parameters that can be set in initial menu
    write(10,"(a)") unitstr
    write(10,*) iunitx
    write(10,*) ibroadfunc
    write(10,*) scalecurve
    write(10,*) gauweigh
    write(10,*) ispin
    write(10,*) icompmethod
    write(10,*) idegen
    write(10,*) degencrit
    write(10,*) iPDOStype
    write(10,*) iCOHP
    write(10,*) enelow
    write(10,*) enehigh
    write(10,*) Yrightsclfac
    !Parameters that can be set in post-processing menu
    write(10,"(a)") graphformat
    write(10,*) ylowerleft
    write(10,*) yupperleft
    write(10,*) stepx
    write(10,*) stepyleft
    write(10,*) ishowTDOScurve
    write(10,*) ishowTDOSline
    write(10,"(10i4)") ishowPDOSline(:)
    write(10,"(10i4)") ishowPDOScurve(:)
    write(10,*) ishowOPDOScurve
    write(10,*) ishowOPDOSline
    write(10,*) ishowlegend
    write(10,*) ishowHOMOlev
    write(10,*) ishowvertline,vertline_X
    write(10,*) ishowYlab
    write(10,*) legendx
    write(10,*) legendy
    write(10,*) icurvewidth
    write(10,*) ilinewidth
    write(10,"(10i4)") iclrPDOS(:)
    write(10,"(a)") TDOSstring
    write(10,"(a)") OPDOSstring
    do ifrag=1,nfragmax
        write(10,"(a)") PDOSstring(ifrag)
    end do
    !Only for 2D-LDOS
    write(10,*) yxratio
    write(10,*) clrscllow
    write(10,*) clrsclhigh
    write(10,*) clrsclstep
    write(10,*) Yleftstep
    !Status
    write(10,*) iusersetcolorscale
    write(10,*) iusersetY
    !System specific information
    write(10,*) nbasis
    write(10,*) imoend
    do imo=1,nmo
        write(10,"(2E16.8,2f12.6)") MOene_dos(imo),MOocc_dos(imo),FWHM(imo),str(imo)
    end do
    do ifrag=1,nfragmax
        write(10,*) nfragDOS(ifrag)
        if (nfragDOS(ifrag)>0) write(10,"(12i6)") fragDOS(:nfragDOS(ifrag),ifrag)
    end do
    !New parameters, recording with labels
    write(10,"('ilinebottom',i5)") ilinebottom
    write(10,"('eneshift',f12.6)") eneshift
    close(10)
    write(*,*) "Done!"
    cycle
else if (index(c80tmp,'l')/=0) then
    write(*,"(a)") " Input file path to load status from it, e.g. C:\Bang_Dream\RAS.dat"
    write(*,"(a)") " Note: If you press ENTER button directly, status will be load from DOS.dat in current folder"
    read(*,"(a)") c200tmp
    if (c200tmp==" ") c200tmp="DOS.dat"
	inquire(file=c200tmp,exist=alive)
	if (.not.alive) then
	    write(*,*) "Error: Cannot find the file! Press ENTER button to return"
        read(*,*)
        cycle
    end if
    open(10,file=c200tmp,status="old")
    read(10,*) nmotmp
    if (nmotmp/=nmo) then
        write(*,"(a)") " Warning: The status file seems to be inconsistent with current system, do you really want to continue to load it? (y/n)"
        read(*,*) selectyn
        if (selectyn=='n'.or.selectyn=='N') cycle
    end if
    nmo=nmotmp
    !Parameters that can be set in initial menu
    read(10,"(a)") unitstr
    read(10,*) iunitx
    read(10,*) ibroadfunc
    read(10,*) scalecurve
    read(10,*) gauweigh
    read(10,*) ispin
    read(10,*) icompmethod
    read(10,*) idegen
    read(10,*) degencrit
    read(10,*) iPDOStype
    read(10,*) iCOHP
    read(10,*) enelow
    read(10,*) enehigh
    read(10,*) Yrightsclfac
    !Parameters that can be set in post-processing menu
    read(10,"(a)") graphformat
    read(10,*) ylowerleft
    read(10,*) yupperleft
    read(10,*) stepx
    read(10,*) stepyleft
    read(10,*) ishowTDOScurve
    read(10,*) ishowTDOSline
    read(10,"(10i4)") ishowPDOSline(:)
    read(10,"(10i4)") ishowPDOScurve(:)
    read(10,*) ishowOPDOScurve
    read(10,*) ishowOPDOSline
    read(10,*) ishowlegend
    read(10,*) ishowHOMOlev
    read(10,*) ishowvertline,vertline_X
    read(10,*) ishowYlab
    read(10,*) legendx
    read(10,*) legendy
    read(10,*) icurvewidth
    read(10,*) ilinewidth
    read(10,"(10i4)") iclrPDOS(:)
    read(10,"(a)") TDOSstring
    read(10,"(a)") OPDOSstring
    do ifrag=1,nfragmax
        read(10,"(a)") PDOSstring(ifrag)
    end do
    !Only for 2D-LDOS
    read(10,*) yxratio
    read(10,*) clrscllow
    read(10,*) clrsclhigh
    read(10,*) clrsclstep
    read(10,*) Yleftstep
    !Status
    read(10,*) iusersetcolorscale
    read(10,*) iusersetY
    !System specific information
    read(10,*) nbasis
    read(10,*) imoend
    do imo=1,nmo
        read(10,"(2E16.8,2f12.6)") MOene_dos(imo),MOocc_dos(imo),FWHM(imo),str(imo)
    end do
    do ifrag=1,nfragmax
        read(10,*) nfragDOS(ifrag)
        if (nfragDOS(ifrag)>0) read(10,"(12i6)") fragDOS(:nfragDOS(ifrag),ifrag)
    end do
    call readoption_int(10,"ilinebottom",' ',ilinebottom)
    call readoption_float(10,"eneshift",' ',eneshift)
    close(10)
    write(*,*) "Loading finished!"
    cycle
else
	if (c80tmp=="-0") then
		idoOPDOS=0
    else if (c80tmp=="00") then
		idoOPDOS=1
        idoPDOS=0
    end if
    read(c80tmp,*) isel
end if


if (isel==-10) then
    graphformat=graphformat_old
	exit
    
else if (isel==-7) then
	if (iCOHP==0) then
		iCOHP=1
        icompmethod=1 !Basis function based fragment like Mulliken partition
        scalecurve=1
		ishowTDOScurve=0
		ishowTDOSline=0
		ishowOPDOScurve=0
		ishowOPDOSline=0
		ishowCOHPline=1
        if (.not.allocated(FmatA)) then
			do while(.true.)
				write(*,*)
				write(*,*) "How to provide the Fock/KS matrix used for plotting COHP"
				write(*,"(a)") " 1 Generating Fock/KS matrix by MO energies and coefficients as well as overlap matrix"
				write(*,*) "2 Loading Fock/KS matrix from a file"
				read(*,*) isel2
				if (isel2==1) then
					call MOene2Fmat(istatus)
				else if (isel2==2) then
					call loadFockfile(istatus)
				end if
				if (istatus==0) exit
			end do
        end if
		!Use eV unit
        if (iunitx==1) then
			iunitx=2
			unitstr=" eV"
			MOene_dos=MOene_dos*au2eV
			!str=str/au2eV     !I hope only for COHP, when X-axis use (eV), one unit of eV broadens to a normalized integral curve. So do not scale str
        end if
        if (ifPBC>0) str=str*au2eV !PBC case uses eV as X-axis unit by default and has scaled str by 1/au2eV, so scale it back
        ilinebottom=0 !Incompatible with showing lines at bottom
		FWHM=0.5D0 !Usually suitable
        call getHOMOidx
        !Set proper X-range for plotting COHP
        enelow=MOene_dos(idxHOMO)-10
        enehigh=MOene_dos(idxHOMO+1)+3
		stepx=(enehigh-enelow)/10
        if (ifPBC>0) then
			write(*,*)
			write(*,*) "Do you want to perform Lowdin orthonormalization? (y/n)"
			write(*,"(a)") " Note: For periodic systems calculated with a minimal basis set, choosing ""y"" will lead to more meaningful result"
			read(*,*) c80tmp
			if (c80tmp=='y'.or.c80tmp=='Y') then
				call ask_Sbas_PBC !Sbas is needed in transformation to orthonormal basis
				write(*,*) "Transforming coefficient and Fock/KS matrices to Lowdin orthonormal basis..."
				call symmortho(1)
				write(*,*) "Done!"
			end if
        end if
    else
		write(*,"(a)") " NOTE: Cannot change to DOS mode from COHP mode, please re-enter present module. Press ENTER button to continue"
        read(*,*)
    end if
    
else if (isel==-6) then !Shift energy levels
	eneshift_old=eneshift
	!Remove previous shift
	MOene_dos=MOene_dos-eneshift
    write(*,*) "Input shift value of orbital energies in "//trim(unitstr)//", e.g. 0.23"
    write(*,*) "If input ""H"", then orbital energies will be added by negative of HOMO energy"
    read(*,*) c80tmp
    if (c80tmp=='h'.or.c80tmp=='H') then
		eneHOMO=-1D99 !Determine HOMO energy
		do imo=1,imoend
			irealmo=imo
			if (ispin==2) irealmo=imo+nbasis
			if (MOocc_dos(irealmo)>0.and.MOene_dos(irealmo)>eneHOMO) eneHOMO=MOene_dos(irealmo)
		end do
        eneshift=-eneHOMO
    else
		read(c80tmp,*) eneshift
    end if
	MOene_dos=MOene_dos+eneshift
    !Correspondingly shift X-range
    enelow=enelow+(eneshift-eneshift_old)
    enehigh=enehigh+(eneshift-eneshift_old)
    write(*,*) "Done!"
    
else if (isel==-5) then
	do while(.true.)
		write(*,*)
		write(*,*) "0 Return"
		write(*,*) "1 Set orbital energies for specific orbitals"
		write(*,*) "2 Set occupation numbers for specific orbitals"
		write(*,*) "3 Set strengths for specific orbitals"
		write(*,*) "4 Set FWHMs for specific orbitals"
		read(*,*) isel
		
		if (isel==0) then
			exit
		else 
			write(*,"(a)") " Input orbital indices. e.g. 1,3-6,8,10-11 means orbitals 1,3,4,5,6,8,10,11 will be selected"
			if (ispin==1) write(*,*) "Note: The indices are for alpha ones"
			if (ispin==2) write(*,*) "Note: The indices are for beta ones, starting from 1"
			if (ispin==3) write(*,"(a,i6)") " Note: Alpha orbital index starts from 1, beta index starts from",nbasis+1
			read(*,"(a)") c200tmp
			call str2arr(c200tmp,nselorb)
			if (allocated(selorbarr)) deallocate(selorbarr)
			allocate(selorbarr(nselorb))
			call str2arr(c200tmp,nselorb,selorbarr)
			if (isel==1) then
				write(*,*) "Set their energies to which value? e.g. -0.13"
				write(*,*) "Note: The value should be given in"//unitstr
				read(*,*) enetmp
				if (ispin==2) selorbarr=selorbarr+nbasis
				do imoidx=1,nselorb
					imo=selorbarr(imoidx)
					MOene_dos(imo)=enetmp
				end do
			else if (isel==2) then
				write(*,*) "Set their occupation numbers to which value? e.g. 2.0"
				read(*,*) occtmp
				if (ispin==2) selorbarr=selorbarr+nbasis
				do imoidx=1,nselorb
					imo=selorbarr(imoidx)
					MOocc_dos(imo)=occtmp
				end do
			else if (isel==3) then
				write(*,*) "Set their strength to which value? e.g. 1.0"
				read(*,*) strtmp
				do imoidx=1,nselorb
					imo=selorbarr(imoidx)
					str(imo)=strtmp
				end do
			else if (isel==4) then
				write(*,*) "Set their FWHM to which value? e.g. 0.05"
				read(*,*) FWHMtmp
				do imoidx=1,nselorb
					imo=selorbarr(imoidx)
					FWHM(imo)=FWHMtmp
				end do
			end if
		end if
	end do
    
else if (isel==-4) then
	if (ispin==1.or.ispin==3) write(*,*) "Below orbitals are Alpha type"
	if (ispin==2) write(*,*) "Below orbitals are Beta type"
	do imo=1,imoend
		irealmo=imo
		if (ispin==2) irealmo=imo+nbasis
		if (ispin==3.and.imo==nbasis+1) write(*,"(/,a)") " Below orbitals are Beta type"
		write(*,"(' #',i6,'  Energy (',a,'):',f15.4,'  Occ:',f6.3,'  Str:',f6.3,'  FWHM:',f8.4)") &
        imo,trim(unitstr),MOene_dos(irealmo),MOocc_dos(irealmo),str(imo),FWHM(imo)
	end do
	write(*,*)
    
else if (isel==-3) then
	open(10,file="orginfo.txt",status="replace")
	if (iunitx==1) write(10,"(2i6)") imoend,2 !Current is a.u.
    if (iunitx==2) write(10,"(2i6)") imoend,4 !Current is eV
	do imo=1,imoend
		irealmo=imo
		if (ispin==2) irealmo=imo+nbasis
		write(10,"(f18.6,3f12.6)") MOene_dos(irealmo),MOocc_dos(irealmo),str(imo),FWHM(imo)
	end do
	close(10)
	write(*,"(a)") " The energy levels, occupation numbers, strengths, FWHMs have been exported to orginfo.txt in current directory, &
	&you can modify it and then load it into Multiwfn again"
	if (iunitx==1) write(*,*) "Note: The unit of energy levels and FWHMs in this file is a.u."
	if (iunitx==2) write(*,*) "Note: The unit of energy levels and FWHMs in this file is eV"
	write(*,*)

else if (isel==-2) then !Define MO sets for MO-PDOS
    if (iPDOStype==1) then
        write(*,"(a)") " Warning: You have defined atom or basis function fragments, which conflict with MO fragments. &
        &To proceed, these fragments will be cleaned, OK? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y') then
            nfragDOS=0
            iPDOStype=0
        else
            cycle
        end if
    end if
    write(*,"(a)") " Note 1: You can use options 1~10 to define up to 10 MO sets for plotting respective PDOS"
    write(*,"(a)") " Note 2: The first published paper employing MO-PDOS map is: Zeyu Liu, Tian Lu, Qinxue Chen, &
    &Carbon, 165, 461 (2020) DOI: 10.1016/j.carbon.2020.05.023. Please cite this paper if MO-PDOS map is employed in your work, thank you!"
    do while(.true.)
        write(*,*)
		write(*,*) "     ---------------- Define fragments for MO-PDOS map ----------------"
        write(*,*) " 0 Return"
        do iset=1,nfragmax
            if (nfragDOS(iset)==0) then
                write(*,"(i3,' This set is undefined')") iset
            else
                write(*,"(i3,' There are',i6,' orbitals')") iset,nfragDOS(iset)
            end if
        end do
        read(*,*) iset
        if (iset==0) then
            if (any(nfragDOS>0)) then
                iPDOStype=2
            else
                iPDOStype=0
            end if
            exit
        else
            write(*,*) "Input index of the MOs, e.g. 2,3,7-10,23"
            if (nfragDOS(iset)>0) write(*,*) "If input 0, then this set will be unset"
            read(*,"(a)") c2000tmp
            if (c2000tmp(1:1)=="0") then
                nfragDOS(iset)=0
            else
                call str2arr(c2000tmp,nfragDOS(iset),fragDOS(:,iset))
            end if
        end if
    end do
    
else if (isel==-1) then !Define fragments in common sense
    if (iPDOStype==2) then
        write(*,"(a)") " Warning: You have defined MO fragments, which conflicts with present fragment setting. To proceed, &
        &the MO fragments will be cleaned, OK? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y') then
            nfragDOS=0
            iPDOStype=0
        else
            cycle
        end if
    end if
    if (idegen==1) then
        write(*,"(a)") " Note: Degeneracy cannot be shown on the map when plotting PDOS, therefore the degeneracy setting has been disabled"
        idegen=0
    end if
    write(*,*)
	write(*,*) "           ----------------- Define fragments -----------------"
    if (iCOHP==0) then
		write(*,"(a)") " Note: Up to 10 fragments can be defined for plotting PDOS, but OPDOS will only be plotted for fragments 1 and 2 (they can have overlap, while local terms will not be calculated)"
    else
		write(*,"(a)") " Note: To plot COHP between fragments, you need to define fragments 1 and 2 (they can have overlap, while local terms will not be calculated)"
    end if
	do while(.true.)
		write(*,*)
        if (iCOHP==0) then
			do ifrag=1,nfragmax
				if (nfragDOS(ifrag)==0) then
					write(*,"(' Fragment',i5,', has not been defined')") ifrag
				else
					if (icompmethod<=2) write(*,"(' Fragment',i5,', number of basis functions:',i6)") ifrag,nfragDOS(ifrag)
					if (icompmethod>=3) write(*,"(' Fragment',i5,', number of atoms:',i6)") ifrag,nfragDOS(ifrag)
				end if
			end do
        else !COHP
			do ifrag=1,2
				if (nfragDOS(ifrag)==0) then
					write(*,"(' Fragment',i5,', has not been defined')") ifrag
				else
					write(*,"(' Fragment',i5,', number of basis functions:',i6)") ifrag,nfragDOS(ifrag)
				end if
            end do
        end if
		write(*,*) "Input fragment index to define it, e.g. 2"
		write(*,*) "Input a negative index will unset the fragment, e.g. -2"
		write(*,*) "Input two indices can exchange the two fragments, e.g. 1,4"
		write(*,*) "Input ""e"" can export current fragment setting to DOSfrag.txt in current folder"
		write(*,*) "Input ""i"" can import fragment setting from DOSfrag.txt in current folder"
		write(*,*) "To return to the last menu, input 0 or q"
		read(*,"(a)") c80tmp
		if (index(c80tmp(1:len_trim(c80tmp)),' ')/=0.or.c80tmp==" ") then
			write(*,*) "Inputting error!"
        else if (c80tmp=='0'.or.c80tmp=='q') then
            if (any(nfragDOS>0)) then
                iPDOStype=1
            else
                iPDOStype=0
            end if
            exit
		else if (c80tmp=='e') then !Export fragment definition
			open(10,file="DOSfrag.txt",status="replace")
			do ifrag=1,nfragmax
				write(10,*)
				if (icompmethod<=2) write(10,"(' #Fragment:',i4,'   nbasis:',i8)") ifrag,nfragDOS(ifrag)
				if (icompmethod>=3) write(10,"(' #Fragment:',i4,'   natoms:',i8)") ifrag,nfragDOS(ifrag)
				write(10,"(8i8)") fragDOS(1:nfragDOS(ifrag),ifrag)
			end do
			close(10)
			write(*,*) "Exporting finished!"
		else if (c80tmp=='i') then !Import fragment definition
			open(10,file="DOSfrag.txt",status="old")
			do ifrag=1,nfragmax
				read(10,*)
				read(10,*) c80tmp,inouse,c80tmp,nfragDOS(ifrag)
				read(10,"(8i8)") fragDOS(1:nfragDOS(ifrag),ifrag)
			end do
			close(10)
			write(*,*) "Importing finished!"
		else if (index(c80tmp,',')==0) then !Define a fragment
			read(c80tmp,*,iostat=ierror) ifragsel
            if (ierror/=0) then
			    write(*,*) "Inputting error!"
			    cycle
            end if
			if (ifragsel>nfragmax) then
				write(*,*) "Error: The index exceeded upper limit of fragments!"
			else if (ifragsel<0) then
				nfragDOS(abs(ifragsel))=0
			else !deffrag routine is only able to deal with global array frag1 and frag2, so we use frag1 as intermediate array
                if (icompmethod==1.or.icompmethod==2) then !Set basis functions in specific fragment
				    allocate(frag1(nfragDOS(ifragsel)))
				    frag1(:)=fragDOS(1:nfragDOS(ifragsel),ifragsel)
				    call deffrag(1)
				    if (allocated(frag1)) then
					    nfragDOS(ifragsel)=size(frag1)
					    fragDOS(1:nfragDOS(ifragsel),ifragsel)=frag1(:)
					    deallocate(frag1)
				    else
					    nfragDOS(ifragsel)=0
			    	end if
                else if (icompmethod==3.or.icompmethod==4) then !Set atoms in specific fragment
                    write(*,*) "Input index of the atoms comprising the fragment, e.g. 2,3,7-10"
                    write(*,*) "You can also input an element to choose all corresponding atoms, e.g. Fe"
                    read(*,"(a)") c2000tmp
                    if (iachar(c2000tmp(1:1))<48.or.iachar(c2000tmp(1:1))>57) then !Input an element
						call lc2uc(c2000tmp(1:1))
						call uc2lc(c2000tmp(2:2))
                        nfragDOS(ifragsel)=count(a%name==c2000tmp(1:2))
						if (nfragDOS(ifragsel)==0) then
							write(*,*) "Error: No atom is selected"
                        else
							itmp=0
							do iatm=1,ncenter
								if (a(iatm)%name==c2000tmp(1:2)) then
									itmp=itmp+1
                                    fragDOS(itmp,ifragsel)=iatm
                                end if
                            end do
                            write(*,"(i8,' atoms were selected')") itmp
                        end if
                    else
						call str2arr(c2000tmp,ntmp)
						nfragDOS(ifragsel)=ntmp
						call str2arr(c2000tmp,ntmp,fragDOS(1:ntmp,ifragsel))
						if (any(fragDOS(1:ntmp,ifragsel)<0).or.any(fragDOS(1:ntmp,ifragsel)>ncenter)) then
							write(*,*) "Error: The atom index exceeded valid range! You must redefine it"
							nfragDOS(ifragsel)=0
						end if
                    end if
                end if
			end if
		else !Exchange fragments
			read(c80tmp,*,iostat=ierror) ifragsel,jfragsel
            if (ierror/=0) then
			    write(*,*) "Inputting error!"
			    cycle
            end if
			ntmp=nfragDOS(jfragsel)
			nfragDOS(jfragsel)=nfragDOS(ifragsel)
			nfragDOS(ifragsel)=ntmp
			fragDOS(:,nfragmax+1)=fragDOS(:,jfragsel)
			fragDOS(:,jfragsel)=fragDOS(:,ifragsel)
			fragDOS(:,ifragsel)=fragDOS(:,nfragmax+1)
			write(*,*) "Exchanging finished!"
			write(*,*)
		end if
	end do
    
else if (isel==1) then
    write(*,*) "Choose one of broadening functions:"
	write(*,*) "1 Lorentzian"
	write(*,*) "2 Gaussian"
	write(*,*) "3 Pseudo-Voigt (i.e. Mixed Lorentzian and Gaussian)"
	read(*,*) ibroadfunc
    
else if (isel==2) then
	if (iunitx==1) then
		write(*,*) "Input lower, upper limits and stepsize between labels (in a.u.)"
		write(*,*) "e.g. -1.5,0.2,0.3"
	else if (iunitx==2) then
		write(*,*) "Input lower, upper limits and stepsize between labels (in eV)"
		write(*,*) "e.g. -20,5,2"
	end if
	read(*,*) enelow,enehigh,stepx
    
else if (isel==3) then
	write(*,*) "Input a value, e.g. 0.03"
	read(*,*) FWHMtmp
	if (FWHMtmp<0D0) write(*,*) "Error: The value should larger than zero, input again"
	FWHM=FWHMtmp
else if (isel==4) then
	write(*,*) "Input a value, e.g. 0.2"
	read(*,*) scalecurve
	if (scalecurve<0D0) write(*,*) "Error: The value should larger than zero, input again"
    
else if (isel==5) then
	write(*,*) "Input a value, e.g. 0.3"
	read(*,*) gauweigh
	if (gauweigh<0D0) write(*,*) "Error: The value should larger than zero, input again"
    
else if (isel==6) then
	write(*,*) "1 Alpha spin"
	write(*,*) "2 Beta spin"
	write(*,*) "3 Both spins"
	read(*,*) ispin
	if (ispin==1.or.ispin==2) then
		imoend=nbasis
	else if (ispin==3) then
		imoend=nmo
    end if

else if (isel==7) then
    ioldmethod=icompmethod
    write(*,"(a,/)") " Hint: Options 1 and 2 are very fast, however the result is not robust for unoccupied MOs, &
    &and even useless when diffuse functions are employed"
    write(*,*) "1 Mulliken method"
    write(*,*) "2 SCPA method"
    write(*,*) "3 Hirshfeld method"
    write(*,*) "4 Becke method"
    read(*,*) icompmethod
    if ( any(nfragDOS>0) .and. ((ioldmethod<=2.and.icompmethod>=3).or.(ioldmethod>=3.and.icompmethod<=2)) ) then !Need redefine fragment
        write(*,*) "Fragment definitions will be cleaned, OK? (y/n)"
        read(*,*) selectyn
        if (selectyn=='n'.or.selectyn=='N') cycle
        nfragDOS=0
    end if
    if (icompmethod==3) then !Use cheap grid to compute atomic contributions to all MOs
        call gen_orbatmcomp_space(1,atmcomp(:,:),1,nmo,1,0)
    else if (icompmethod==4) then
        call gen_orbatmcomp_space(2,atmcomp(:,:),1,nmo,1,0)
    end if
    
else if (isel==8) then
	if (iunitx==1) then !a.u.->eV
		iunitx=2
		MOene_dos=MOene_dos*au2eV
		FWHM=FWHM*au2eV
		enelow=enelow*au2eV
		enehigh=enehigh*au2eV
		unitstr=" eV"
        eneshift=eneshift*au2eV
		!After change the unit, in principle, the curve (and hence Y-range) will be automatically reduced by 27.2114.& 
		!str should also be reduced by 27.2114 so that the discrete line can be properly shown in the graph range &
		!To compensate the reduce of str, scalecurve thus be augmented by corresponding factor
		str=str/au2eV
		scalecurve=scalecurve*au2eV
		stepx=stepx*au2eV
		stepyleft=stepyleft/au2eV
	else if (iunitx==2) then !eV->a.u.
		iunitx=1
		MOene_dos=MOene_dos/au2eV
		FWHM=FWHM/au2eV
		enelow=enelow/au2eV
		enehigh=enehigh/au2eV
        eneshift=eneshift/au2eV
		unitstr=" a.u."
		str=str*au2eV
		scalecurve=scalecurve/au2eV
		stepx=stepx/au2eV
		stepyleft=stepyleft*au2eV
	end if
    
else if (isel==9) then
    if (idegen==1) then
        idegen=0
    else
        if (iPDOStype==1) then
            write(*,*) "Error: This feature cannot be enabled when plotting PDOS"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            cycle
        end if
        write(*,*) "Input threshold for determining degenerate in eV, e.g. 0.01"
        write(*,*) "If press ENTER button directly, 0.005 eV will be used"
        read(*,"(a)") c80tmp
        if (c80tmp==" ") then
            degencrit=0.005D0
        else
            read(c80tmp,*) degencrit
        end if
        idegen=1
    end if
	

!!!!!!!!-------------- Now we compute and plot DOS/COHP (isel==0) or LDOS curve (isel==10)
else if (isel==0.or.isel==10) then
	if (iCOHP==1) then
		if (nfragDOS(1)>0.and.nfragDOS(2)>0) then
			write(*,"(a)") " NOTE: COHP between the two fragments you defined will be plotted"
        else if (nfragDOS(1)==0.and.nfragDOS(2)==0) then
			write(*,"(a)") " NOTE: COHP between nearest atoms will be plotted. While if you want to plot COHP between two fragments, you need to define them by option -1 first"
        else
			write(*,"(a)") " Error: You should not only define one fragment! Please unset it (to plot COHP between all atoms) or &
            &also define another fragment (to plot COHP between two fragments) in option -1"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
        end if
    end if
	if (isel==10) then
		write(*,*) "Input X/Y/Z coordinate (in Bohr), e.g. 1.0,1.5,0.2"
		read(*,*) x,y,z
	end if
    
	!======Generate fragment composition, overlap population, COHP=======
	
	!Reset display setting
	ishowPDOScurve(:)=0
	ishowPDOSline(:)=0
    if (iCOHP==0) then
		do ifrag=1,nfragmax
			if (nfragDOS(ifrag)>0) then
				ishowPDOScurve(ifrag)=1
				ishowPDOSline(ifrag)=1
			end if
		end do
    end if
	ishowOPDOScurve=0
	ishowOPDOSline=0
	if (idoOPDOS==1) then
		ishowOPDOScurve=1
		ishowOPDOSline=1
	end if
	ishowLDOScurve=0
	ishowLDOSline=0
	if (isel==10) then
		ishowLDOScurve=1
		ishowLDOSline=1
	end if
    
    FWHMmax=maxval(FWHM)
    if (iunitx==1) then !a.u.
        enediff=degencrit/au2eV !The degencrit, which is used to detect orbital degeneracy, is always in eV, so convert to a.u.
    else if (iunitx==2) then !eV
        enediff=degencrit
    end if
	if (idoPDOS==1.or.idoOPDOS==1.or.iCOHP==1) then !Calculate composition used for plotting PDOS, OPDOS, or COHP data
        !if (iCOHP==1.and.nfragDOS(1)==0.and.nfragDOS(2)==0) call genconnmat(1,1) !Generate connectivity matrix for calculating COHP between neighbouring atoms
        compfrag=0
		if (iPDOStype/=2) then !Not MO-PDOS fragment, namely basis function or atom fragment for PDOS, OPDOS, COHP
			if (ispin/=3) then
				ntime=1 !One set of orbital will be calculated
			else
				ntime=2 !Unrestricted wavefunction and both spins are considered, twice calculation respectively using CObas of different spin
			end if
			if (icompmethod<=2) then !Mulliken/SCPA orbital composition analysis method, including COHP case
				if ((idoPDOS==1.and.icompmethod==1).or.idoOPDOS==1) call ask_Sbas_PBC !Generate overlap matrix if not available, which is needed by Mulliken or OPDOS
				OPfrag12=0
				call walltime(iwalltime1)
				do itime=1,ntime
					if (ispin/=3) then !ispin=0 (closed-shell), =1/2 (alpha/beta)
						if (iCOHP==0) write(*,*) "Calculating orbital composition, please wait..."
						if (iCOHP==1) write(*,*) "Calculating COHP, please wait..."
						tmpmat=>CObasa !Calculate total or alpha
						if (ispin==2) tmpmat=>CObasb !Calculate beta
                        if (iCOHP==1) then
							tmpmat2=>FmatA
                            if (ispin==2) tmpmat2=>FmatB
                        end if
					else !ispin==3, mix alpha and beta
						if (itime==1) then !Alpha spin
							if (iCOHP==0) write(*,*) "Calculating alpha orbital composition, please wait..."
							tmpmat=>CObasa
							if (iCOHP==1) then
								write(*,*) "Calculating alpha COHP, please wait..."
								tmpmat2=>FmatA
                            end if
						else if (itime==2) then !Beta spin
							if (iCOHP==0) write(*,*) "Calculating beta orbital composition, please wait..."
							tmpmat=>CObasb
							if (iCOHP==1) then
								write(*,*) "Calculating beta COHP, please wait..."
								tmpmat2=>FmatB
                            end if
						end if
					end if
					!Count how many orbitals will be evaluated this time
					ncalc=0
					do imo=1,nbasis
						imoall=imo !The orbital index defined from 1 to nmo
						if (ispin==2.or.itime==2) imoall=imo+nbasis
						if (MOene_dos(imoall)<enelow-3*FWHMmax.or.MOene_dos(imoall)>enehigh+3*FWHMmax.or.MOene(imoall)==0) cycle
						ncalc=ncalc+1 
					end do
					write(*,"(i8,' orbitals will be calculated')") ncalc
					!Calculate orbital compositions
					iprog=0
					!$OMP PARALLEL DO SHARED(compfrag,OPfrag12) PRIVATE(ifrag,imo,imoall,imoslot,allsqr,i,j,ibas,jbas,iatm,jatm,tmpval,iclose) schedule(dynamic) NUM_THREADS(nthreads)
					do imo=1,nbasis !Cycle MOs
						imoslot=imo !The index to be placed into arrays (used to generate PDOS). It starts from 1 even only beta is considered
						if (itime==2) imoslot=imo+nbasis !Both spin case
						imoall=imo !The orbital index defined from 1 to nmo
						if (ispin==2.or.itime==2) imoall=imo+nbasis
						if (MOene_dos(imoall)<enelow-3*FWHMmax.or.MOene_dos(imoall)>enehigh+3*FWHMmax) cycle
						if (MOene(imoall)==0) cycle !They are not actual orbitals recorded in inputted wavefunction file, which are abundant in CP2K molden file
						if (idoPDOS==1) then !Calculate fragment composition in MOs
							if (icompmethod==2) allsqr=sum(tmpmat(:,imo)**2) !Used by SCPA
							do ifrag=1,nfragmax
								if (nfragDOS(ifrag)==0) cycle
								do i=1,nfragDOS(ifrag) !Cycle each basis in the fragment
									ibas=fragDOS(i,ifrag)
									if (icompmethod==2) then !SCPA
										compfrag(imoslot,ifrag)=compfrag(imoslot,ifrag)+tmpmat(ibas,imo)**2/allsqr
									else !Mulliken
										do jbas=1,nbasis !Cycle all basis, included inner&external cross term and local term
											compfrag(imoslot,ifrag)=compfrag(imoslot,ifrag)+tmpmat(ibas,imo)*tmpmat(jbas,imo)*Sbas(jbas,ibas)
										end do
									end if
								end do
							end do
                        end if
 						!Calculate overlap population or COHP
						if (idoOPDOS==1.or.iCOHP==1) then
							if (nfragDOS(1)==0.and.nfragDOS(2)==0) then !Calculate OPDOS or COHP between nearest atoms
								tmpval=0
								do iatm=1,ncenter
									call nearest_atom(iatm,iclose)
									do ibas=basstart(iatm),basend(iatm)
										do jbas=basstart(iclose),basend(iclose)
											if (iCOHP==1) then !Overlap Hamilton
												tmpval=tmpval+tmpmat(ibas,imo)*tmpmat(jbas,imo)*tmpmat2(ibas,jbas)
                                            else !Overlap population
												tmpval=tmpval+tmpmat(ibas,imo)*tmpmat(jbas,imo)*Sbas(ibas,jbas)
                                            end if
										end do
									end do
                                    !Between all atoms
									!do jatm=1,ncenter
									!	if (jatm==iatm) cycle
									!	do ibas=basstart(iatm),basend(iatm)
									!		do jbas=basstart(jatm),basend(jatm)
									!			if (iCOHP==1) then !Overlap Hamilton
									!				tmpval=tmpval+tmpmat(ibas,imo)*tmpmat(jbas,imo)*tmpmat2(ibas,jbas)
									!			else !Overlap population
									!				tmpval=tmpval+tmpmat(ibas,imo)*tmpmat(jbas,imo)*Sbas(ibas,jbas)
									!			end if
									!		end do
									!	end do
         !                           end do
                                    !Consider all bonded atom pairs
									!do jatm=iatm+1,ncenter
									!	if (connmat(iatm,jatm)/=0) then
									!		do ibas=basstart(iatm),basend(iatm)
									!			do jbas=basstart(jatm),basend(jatm)
									!				tmpval=tmpval+tmpmat(ibas,imo)*tmpmat(jbas,imo)*tmpmat2(jbas,ibas)*2
									!			end do
									!		end do
         !                               end if
         !                           end do
                                end do
                                OPfrag12(imoslot)=OPfrag12(imoslot)+tmpval
                            else !Between two defined fragments
								do i=1,nfragDOS(1)
									ibas=fragDOS(i,1)
									do j=1,nfragDOS(2)
										jbas=fragDOS(j,2)
										if (ibas==jbas) cycle !Skip local term (on-site interaction)
										if (iCOHP==0) then !Overlap population
											OPfrag12(imoslot)=OPfrag12(imoslot)+tmpmat(ibas,imo)*tmpmat(jbas,imo)*Sbas(jbas,ibas)
										else !Overlap Hamilton
											OPfrag12(imoslot)=OPfrag12(imoslot)+tmpmat(ibas,imo)*tmpmat(jbas,imo)*tmpmat2(jbas,ibas)
										end if
									end do
								end do
                            end if
						end if
						!$OMP CRITICAL
						iprog=iprog+1
						call showprog(iprog,ncalc)
						!$OMP END CRITICAL
					end do
					!$OMP END PARALLEL DO
				end do
                if (iprog<ncalc) call showprog(ncalc,ncalc)
				if (idoOPDOS==1) then
					OPfrag12(:)=OPfrag12(:)*2 !Only (ibas,jbas) is considered above, now also include (jbas,ibas) contribution to overlap population
                else if (iCOHP==1) then
					OPfrag12(:)=OPfrag12(:)*2*au2eV !Also convert to eV unit to enlarge magnitude
                end if
				call walltime(iwalltime2)
				write(*,"(' Calculation of data took up wall clock time',i10,' s')") iwalltime2-iwalltime1
            
			else if (icompmethod>=3) then !Hirshfeld/Becke orbital composition analysis method, the compositions are not calculated here because they have already been calculated when switching to these methods in the interface
				do itime=1,ntime !1: alpha, 2: beta
					do imo=1,nbasis !Cycle MOs
						imoslot=imo !The index to be placed into arrays, count from 1 even for beta only
						if (itime==2) imoslot=imo+nbasis
						imoall=imo !The index defined from 1 to nmo
						if (ispin==2.or.itime==2) imoall=imo+nbasis
						do ifrag=1,nfragmax
							if (nfragDOS(ifrag)==0) cycle
							do i=1,nfragDOS(ifrag) !Cycle each atom in the fragment
								iatm=fragDOS(i,ifrag)
								compfrag(imoslot,ifrag)=compfrag(imoslot,ifrag)+atmcomp(iatm,imoall)
							end do
						end do
					end do
				end do
			end if
        
		else if (iPDOStype==2) then !MO-PDOS
			!For various fragments
			do iset=1,nfragmax
				if (nfragDOS(iset)==0) cycle
				if (idegen==1) then !Determine degeneracy and record as "compfrag"
					istart=1
					do while(istart<=nfragDOS(iset))
						ndegen=1
						imo=fragDOS(istart,iset)
						do jdx=istart+1,nfragDOS(iset)
							jmo=fragDOS(jdx,iset)
							if (abs(MOene_dos(jmo)-MOene_dos(imo))<enediff) ndegen=ndegen+1
						end do
						compfrag(imo,iset)=ndegen
						istart=istart+ndegen
					end do
				else
					do idx=1,nfragDOS(iset)
						imo=fragDOS(idx,iset)
						compfrag(imo,iset)=1
					end do
				end if
			end do
		end if
	end if
    
    !Set degeneracy for all MOs, can be used for MO-PDOS or normal TDOS map
    if (idegen==1) then
        imo=1
        do while(imo<=imoend)
            ndegen=1
            do jmo=imo+1,imoend
                if (abs(MOene_dos(imo)-MOene_dos(jmo))<enediff) ndegen=ndegen+1
            end do
            compfrag(imo,0)=ndegen
            if (ndegen>maxdegen) maxdegen=ndegen
            imo=imo+ndegen
        end do
    else
        compfrag(:,0)=1
    end if
    
    !Calculate LDOS for a point
	if (isel==10) then
		LDOScomp=0
		do imo=1,imoend
			imoall=imo
			if (ispin==2) imoall=imo+nbasis
			if (MOene_dos(imoall)<enelow-3*FWHMmax.or.MOene_dos(imoall)>enehigh+3*FWHMmax) cycle
			LDOScomp(imo)=fmo(x,y,z,imoall)**2
		end do
	end if
	
	!======Set X position of curves==========
	enestep=(enehigh-enelow)/(num1Dpoints-1) 
	do i=1,num1Dpoints
		curvexpos(i)=enelow+(i-1)*enestep
	end do
	
	!======Generate energy levels line=======
	TDOSliney=0
    TDOSliney_unocc=0
	LDOSliney=0
	OPDOSliney=0
    COHPliney=0
    COHPliney_unocc=0
    inow_unocc=0
    iwarnzero=0
	do imo=1,imoend
		inow=3*(imo-1)
		irealmo=imo
		if (ispin==2) irealmo=imo+nbasis
        !In the case of CP2K, virtual orbitals are not solved by default and energies are exactly zero, skip them. &
        !Use MOene rather than MOene_dos because the latter may be shifted
        if (MOene(irealmo)==0) then 
			if (iwarnzero==0) then
				write(*,*) "Note: There are orbitals with zero energy, they are automatically ignored"
				iwarnzero=1
            end if
			cycle
        end if
		linexpos(inow+1:inow+3)=MOene_dos(irealmo)
		if (isel==0) then
			if (iCOHP==0) then
				TDOSliney(inow+2)=str(imo)
				if (idegen==1) TDOSliney(inow+2)=compfrag(imo,0)
				if (MOocc_dos(irealmo)==0) then
					TDOSliney_unocc(inow+2)=str(imo)
					if (idegen==1) TDOSliney_unocc(inow+2)=compfrag(imo,0)
				end if
				do ifrag=1,nfragmax
					if (nfragDOS(ifrag)>0) then
						PDOSliney(inow+1,ifrag)=0D0
						PDOSliney(inow+2,ifrag)=str(imo)*compfrag(imo,ifrag)
						PDOSliney(inow+3,ifrag)=0D0
						!Using height to show degeneracy, do not scale because the axis will be shown at right
						if (idegen==1) PDOSliney(inow+2,ifrag)=compfrag(imo,ifrag)
					end if
				end do
				if (idoOPDOS==1) OPDOSliney(inow+2)=str(imo)*OPfrag12(imo)
            else
				COHPliney(inow+2)=str(imo)*OPfrag12(imo)
				if (MOocc_dos(irealmo)==0) COHPliney_unocc(inow+2)=str(imo)*OPfrag12(imo)
            end if
		else if (isel==10) then
			LDOSliney(inow+2)=str(imo)*LDOScomp(imo)
		end if
	end do
	
	!======Broaden to various kinds of DOS curves=======
	TDOScurve=0D0
	PDOScurve=0D0
	OPDOScurve=0D0
	LDOScurve=0D0
    COHPcurve=0D0
	if (ibroadfunc==1.or.ibroadfunc==3) then !Lorentzian function, see http://mathworld.wolfram.com/LorentzianFunction.html
		do imo=1,imoend !Cycle each orbital
			irealmo=imo
			if (ispin==2) irealmo=imo+nbasis
            if (MOene(irealmo)==0) cycle
			preterm=str(imo)*0.5D0/pi*FWHM(imo)
			do ipoint=1,num1Dpoints !Broaden imo as curve
				tmp=preterm/( (curvexpos(ipoint)-MOene_dos(irealmo))**2+0.25D0*FWHM(imo)**2 )
				if (isel==0) then
					if (iCOHP==0) then
						TDOScurve(ipoint)=TDOScurve(ipoint)+tmp
						do ifrag=1,nfragmax
							if (nfragDOS(ifrag)>0) PDOScurve(ipoint,ifrag)=PDOScurve(ipoint,ifrag)+tmp*compfrag(imo,ifrag)
						end do
						if (idoOPDOS==1) OPDOScurve(ipoint)=OPDOScurve(ipoint)+tmp*OPfrag12(imo)
                    else
						COHPcurve(ipoint)=COHPcurve(ipoint)+tmp*OPfrag12(imo)
                    end if
				else if (isel==10) then
					LDOScurve(ipoint)=LDOScurve(ipoint)+tmp*LDOScomp(imo)
				end if
			end do
		end do
	end if
	if (ibroadfunc==2.or.ibroadfunc==3) then !Gaussian function, see http://en.wikipedia.org/wiki/Gaussian_function
		if (ibroadfunc==3) TDOScurve=(1-gauweigh)*TDOScurve
		do imo=1,imoend !Cycle each orbital
			irealmo=imo
			if (ispin==2) irealmo=imo+nbasis
            if (MOene(irealmo)==0) cycle
			gauss_c=FWHM(imo)/2D0/sqrt(2*dlog(2D0))
			gauss_a=str(imo)/(gauss_c*sqrt(2D0*pi))
			do ipoint=1,num1Dpoints !Broaden imo as curve
				tmp=gauss_a*dexp( -(curvexpos(ipoint)-MOene_dos(irealmo))**2/(2*gauss_c**2) )
				if (ibroadfunc==3) tmp=gauweigh*tmp !Combine Lorentizan and Gaussian function
				if (isel==0) then
					if (iCOHP==0) then
						TDOScurve(ipoint)=TDOScurve(ipoint)+tmp
						do ifrag=1,nfragmax
							if (nfragDOS(ifrag)>0) PDOScurve(ipoint,ifrag)=PDOScurve(ipoint,ifrag)+tmp*compfrag(imo,ifrag)
						end do
						if (idoOPDOS==1) OPDOScurve(ipoint)=OPDOScurve(ipoint)+tmp*OPfrag12(imo)
                    else
						COHPcurve(ipoint)=COHPcurve(ipoint)+tmp*OPfrag12(imo)
                    end if
				else if (isel==10) then
					LDOScurve(ipoint)=LDOScurve(ipoint)+tmp*LDOScomp(imo)
				end if
			end do
		end do
	end if
	TDOScurve=TDOScurve*scalecurve
	PDOScurve=PDOScurve*scalecurve
	OPDOScurve=OPDOScurve*scalecurve
	LDOScurve=LDOScurve*scalecurve
	COHPcurve=COHPcurve*scalecurve
    
    !Calculate center of TDOS and PDOS curves. May be used to calcualte e.g. d-band center
    if (iCOHP==0) then
		TDOSnomin=0
		TDOSdenomin=0
		do ipoint=1,num1Dpoints
			TDOSnomin=TDOSnomin+curvexpos(ipoint)*TDOScurve(ipoint)*enestep
			TDOSdenomin=TDOSdenomin+TDOScurve(ipoint)*enestep
		end do
		write(*,"(/,' Center of TDOS:',f12.6,a)") TDOSnomin/TDOSdenomin,unitstr
		do ifrag=1,nfragmax
			if (nfragDOS(ifrag)>0) then
				PDOSnomin=0
				PDOSdenomin=0
				do ipoint=1,num1Dpoints
					PDOSnomin=PDOSnomin+curvexpos(ipoint)*PDOScurve(ipoint,ifrag)*enestep
					PDOSdenomin=PDOSdenomin+PDOScurve(ipoint,ifrag)*enestep
				end do
				write(*,"(' Center of PDOS',i3,':',f12.6,a)") ifrag,PDOSnomin/PDOSdenomin,unitstr
			end if
		end do
		write(*,*)
    end if
    
    !! All calculations have finished !!
    
	!Initialize some plotting statues, such as range of Y-axis
	idraw=1
    if (isilent==1) idraw=0 !Do not automatically draw at the first time entering the post-processing menu
	isavepic=0
    if (iCOHP==0) then
		if (iusersetY==0) then !Y axis range was not set by user, we automatically determine
			if (isel==0) then
				yupperleft=1.1D0*maxval(TDOScurve)
				ylowerleft=0
				if (idoPDOS==1) ylowerleft=minval(PDOScurve(:,:)) !PDOS may be negative
				if (idoOPDOS==1) ylowerleft=-yupperleft/2 !OPDOS may be large negative value
				if (ylowerleft>0) ylowerleft=0D0 !Don't allow lower plotting limit >0
				stepyleft=nint((yupperleft-ylowerleft)*10)/100D0
			else if (isel==10) then
				ylowerleft=0
				yupperleft=1.1D0*maxval(LDOScurve)
				stepyleft=(yupperleft-ylowerleft)/10
			end if
		end if
		ylowerright=ylowerleft*Yrightsclfac !Lower and upper limit of Y-axis for OPDOS
		yupperright=yupperleft*Yrightsclfac
    else !COHP
		if (iusersetY==0) then
			yupperleft=1.1D0*max(abs(maxval(COHPcurve)),abs(minval(COHPcurve)))
			ylowerleft=-yupperleft
			stepyleft=(yupperleft-ylowerleft)/10
			yupperright=1.3D0*max(abs(maxval(COHPliney)),abs(minval(COHPliney)))
			ylowerright=-yupperright
            stepyright=(yupperright-ylowerright)/10
        end if
    end if
	
	do while(.true.) !Cycle of plotting post-processing menu
    
        !======Draw various kinds of DOS now=======
		!If drawing lines at bottom, we first draw curves (without axis name and labels), and finally use additional section to draw lines
		if (idraw==1) then
			if (isavepic==0) then
				call METAFL('xwin')
				call window(200,100,1000,600)
			else if (isavepic==1) then
				call METAFL(graphformat)
				call winsiz(graph1Dwidth,graph1Dheight)
			end if
			call SCRMOD('REVERSE')
			CALL IMGFMT("RGB")
			CALL PAGE(3000,1800)
			CALL setxid(0,'NONE')
			call DISINI
            call HNAME(height_axis) !Height of axis name
            call height(ticksize) !Size of ticks
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
            leny_curve=1400
            iyaxspos_curve=1550
            leny_line=160 !When showing lines at bottom, leave space for plotting lines
            iyaxspos_line=1550
            if (ilinebottom==1) then
				leny_curve=leny_curve-leny_line
                iyaxspos_curve=iyaxspos_curve-leny_line
            end if
            lenx=2450
            if (iCOHP==1) then
				lenx=2100
			else if (ishowOPDOScurve==1.or.ishowOPDOSline==1) then !Because OPDOS axis will be shown at right side, use smaller X-axis length
				lenx=2300
            else if (idegen==1) then !Because degeneracy axis will be shown at right side, use smaller X-axis length
				lenx=2350
            end if
			call AXSLEN(lenx,leny_curve)
			call AXSPOS(350,iyaxspos_curve)
            if (iCOHP==1) call AXSPOS(450,iyaxspos_curve)
			if (isavepic==0) call WINTIT("Click right mouse button to close")
            
            !Set left axis. Right axis may be set later
            c80tmp="NAME"
            if (ilinebottom==1) c80tmp="LINE" !Do not show axis name if showing lines at bottom
			if (idoOPDOS==1.or.idegen==1.or.iCOHP==1) then
				call setgrf(trim(c80tmp),'NAME','TICKS','LINE')
			else
				call setgrf(trim(c80tmp),'NAME','TICKS','TICKS')
			end if
            if (ishowYlab==0) then
                call labels("NONE","Y")
                CALL TICKS(0,'Y')
                CALL TICKS(1,'X')
                CALL NAMDIS(70,'Y')
            else
                CALL TICKS(1,'XY')
                CALL NAMDIS(40,'Y')
            end if
			call ERRMOD("ALL","OFF")
            if (iCOHP==0) then
				CALL LABDIG(nlabdigX,"X")
				CALL LABDIG(nlabdigY,"Y")
            else
				CALL LABDIG(nlabdigleftY_COHP,"Y")
            end if
            if (ilinebottom==0) then
				if (iunitx==1) CALL NAME('Energy (a.u.)','X')
				if (iunitx==2) CALL NAME('Energy (eV)','X')
            end if
			if (iCOHP==0) then
				CALL NAME('Density-of-states','Y')
            else
				CALL NAME('-COHP curve (dimensionless)','Y')
            end if
			numleg=0
			ileg=1

			if (isel==0) then !Draw regular DOS or COHP
				!Sequence:
				!(1) Set legends
				!(2) Vertical dashed line
				!(3) TDOS (only plot curve if idegen==1)
				!(4) PDOS (only plot curve if idegen==1)
				!(5) TDOS&PDOS lines (idegen==1, using special axis)
				!(6) OPDOS
				!(7) COHP
                !(8) Lines at bottom (if ilinebottom=1, TDOS lines will be plotted in this stage, including TDOS and (MO-)PDOS with/without degen. Do not draw OPDOS lines)
            
				!Set legends. COHP is not involved
				if (ishowTDOScurve==1.or.ishowTDOSline==1) numleg=numleg+1
				do ifrag=1,nfragmax
					if (nfragDOS(ifrag)==0) cycle
					if (ishowPDOScurve(ifrag)==1.or.ishowPDOSline(ifrag)==1) numleg=numleg+1
				end do
				if (ishowOPDOScurve==1.or.ishowOPDOSline==1) numleg=numleg+1
				call legini(clegend,numleg,80)
				call legtit(' ')
				call legopt(2.5D0,0.5D0,1D0) !Decrease the length of legend color line
				call frame(0) !No box around legend
                call legpos(legendx,legendy) !Absolute position of legends
                
                if (iCOHP==1) then
                    if ( max(abs(yupperleft),abs(ylowerleft))>=1000 .or. min(abs(yupperleft),abs(ylowerleft))<0.001D0 ) call labels("FEXP","Y")
                end if
				CALL GRAF(enelow,enehigh,enelow,stepx, ylowerleft,yupperleft,ylowerleft,stepyleft) !Construct left axis!
				
				!Draw a vertical dashed line to highlight HOMO level
				if (ishowHOMOlev==1) then
					intarr(1)=10
					intarr(2)=20
					CALL MYLINE(intarr,2)
					enetmp=-1D99
					do imo=1,imoend
						irealmo=imo
						if (ispin==2) irealmo=imo+nbasis
						if (MOocc_dos(irealmo)>0.and.MOene_dos(irealmo)>enetmp) then
							enetmp=MOene_dos(irealmo)
							iHOMO=irealmo
						end if
					end do
					HOMOlevx=MOene_dos(iHOMO)
					write(*,"(a,f13.5,a)") " Note: The vertical dash line corresponds to HOMO level at",MOene_dos(iHOMO),unitstr
					HOMOlevy(1)=ylowerleft
					HOMOlevy(2)=yupperleft
					call linwid(ilinewidth) !Set to user-defined line width
					CALL CURVE(HOMOlevx,HOMOlevy,2)
					call solid !Restore default
				end if
                if (ishowvertline==1) then !Draw a vertical line at specific X
					HOMOlevx=vertline_X
					HOMOlevy(1)=ylowerleft
					HOMOlevy(2)=yupperleft
					call linwid(ilinewidth) !Set to user-defined line width
					CALL CURVE(HOMOlevx,HOMOlevy,2)
					call solid !Restore default
                end if
                
                !Draw various lines and curves
				if (iCOHP==0) then
					!Draw TDOS
					call linwid(icurvewidth) !Set to user-defined curve width
					if (ishowTDOScurve==1) CALL CURVE(curvexpos,TDOScurve,num1Dpoints) !Draw TDOS curve
					call linwid(ilinewidth) !Set to user-defined line width
					if (ishowTDOSline==1.and.idegen==0.and.ilinebottom==0) then
						call color('WHITE')
						CALL CURVE(linexpos,TDOSliney,3*imoend) !Draw TDOS lines
						call setcolor(6) !Gray
						CALL CURVE(linexpos,TDOSliney_unocc,3*imoend) !Draw TDOS lines for unoccupied MOs
						call color('WHITE')
					end if
					if (ishowTDOScurve==1.or.ishowTDOSline==1) then !Set legend
						call legpat(0,1,-1,-1,-1,ileg)
						CALL LEGLIN(clegend,trim(TDOSstring),ileg)
						ileg=ileg+1
					end if
				
					!Draw PDOS of each defined fragment
					do ifrag=1,nfragmax
						if (nfragDOS(ifrag)>0) then
							call setcolor(iclrPDOS(ifrag))
							call linwid(icurvewidth) !Set to user-defined curve width
							if (ishowPDOScurve(ifrag)==1) CALL CURVE(curvexpos,PDOScurve(:,ifrag),num1Dpoints)
							if (ilinebottom==0) then !If showing lines at bottom, PDOS lines will never be shown
								call linwid(ilinewidth) !Set to user-defined line width
								if (ishowPDOSline(ifrag)==1.and.idegen/=1) CALL CURVE(linexpos,PDOSliney(:,ifrag),3*imoend)
							end if
							if (ishowPDOScurve(ifrag)==1.or.ishowPDOSline(ifrag)==1) then !Set legend
								call legpat(0,1,-1,-1,-1,ileg)
								CALL LEGLIN(clegend,trim(PDOSstring(ifrag)),ileg)
								ileg=ileg+1
							end if
						end if
					end do
					call color('WHITE')
					call linwid(ilinewidth) !Set to user-defined curve width
					call XAXGIT !Draw a black line at Y=0 to cover the colored PDOS line
					call linwid(1) !Recover to default line width of dislin
					call height(legtextsize) !Define legend text size
					if (ishowlegend==1) call legend(clegend,3) !Draw the legends (for TDOS,PDOS), must before endgrf
					call endgrf !Finished drawing TDOS/PDOS
               
					!We do not draw lines before, now draw degeneracy of lines for TDOS and/or MO-PDOS using axis at right side (note: degeneracy is not supported for PDOS)
					if (idegen==1.and.ilinebottom==0) then
						CALL LABDIG(-1,"Y")
						CALL NAME('Degeneracy','Y')
						call height(ticksize)
						call setgrf('NONE','NONE','NONE','NAME')
						CALL GRAF(enelow,enehigh,enelow,stepx, 0D0,10D0,0D0,1D0)
						if (ishowTDOSline==1) then
							call color('WHITE')
							CALL CURVE(linexpos,TDOSliney,3*imoend) !Draw TDOS lines
							call setcolor(6) !Gray
							CALL CURVE(linexpos,TDOSliney_unocc,3*imoend) !Draw TDOS lines for unoccupied MOs
						end if
						do ifrag=1,nfragmax
							if (nfragDOS(ifrag)>0.and.ishowPDOSline(ifrag)==1) then
								call setcolor(iclrPDOS(ifrag))
								call linwid(ilinewidth) !Set to user-defined line width
								CALL CURVE(linexpos,PDOSliney(:,ifrag),3*imoend)
							end if
						end do
						call endgrf
					end if
					!Draw OPDOS
					if (ishowOPDOScurve==1.or.ishowOPDOSline==1) then
						CALL LABDIG(nlabdigY_OPDOS,"Y")
						CALL NAME('OPDOS','Y')
						call height(ticksize) !Size of ticks
						if (ishowYlab==1) CALL NAMDIS(45,'Y')
						call setgrf('NONE','NONE','NONE','NAME')
						CALL GRAF(enelow,enehigh,enelow,stepx, ylowerright,yupperright,ylowerright,(yupperright-ylowerright)/10D0)
						call setcolor(12) !Dark green
						call linwid(icurvewidth) !Set to user-defined curve width
						if (ishowOPDOScurve==1) CALL CURVE(curvexpos,OPDOScurve,num1Dpoints)
						if (ilinebottom==0) then !If showing lines at bottom, OPDOS lines will never be shown
							call linwid(ilinewidth) !Set to user-defined line width
							if (ishowOPDOSline==1) CALL CURVE(linexpos,OPDOSliney,3*imoend)
						end if
						call legpat(0,1,-1,-1,-1,ileg)
						CALL LEGLIN(clegend,trim(OPDOSstring),ileg)
						call color('WHITE')
						call XAXGIT !Draw a black line at Y=0 to cover the colored OPDOS line
						call linwid(1) !Recover to default line width of dislin
						call height(legtextsize) !Define legend text size
						if (ishowlegend==1) call legend(clegend,3) !Draw the legends (for TDOS,PDOS,OPDOS), must before endgrf
						call endgrf
					end if
                else !-COHP
					!Plot COHP curve using left axis
					call linwid(icurvewidth) !Set to user-defined curve width
					CALL CURVE(curvexpos,-COHPcurve,num1Dpoints) !Draw COHP curve
					call linwid(icurvewidth-1)
                    call XAXGIT !Draw a black line at Y=0
					call linwid(1) !Recover to default line width of dislin
                    call endgrf
					!Plot COHP line using right axis
					if (ishowCOHPline==1.and.idegen==0.and.ilinebottom==0) then
						CALL NAME('-COHP data (eV)','Y')
						call setgrf('NONE','NONE','NONE','NAME')
						CALL LABDIG(nlabdigrightY_COHP,"Y")
                        if ( max(abs(yupperright),abs(ylowerright))>=1000 .or. min(abs(yupperright),abs(ylowerright))<0.001D0 ) call labels("FEXP","Y")
						CALL GRAF(enelow,enehigh,enelow,stepx, ylowerright,yupperright,ylowerright,stepyright)
						call linwid(ilinewidth) !Set to user-defined line width
						call color('WHITE')
						CALL CURVE(linexpos,-COHPliney,3*imoend) !Draw COHP lines
						call setcolor(6) !Gray
						CALL CURVE(linexpos,-COHPliney_unocc,3*imoend) !Draw COHP lines for unoccupied MOs
						call color('WHITE')
						call XAXGIT !Draw a black line at Y=0, can also cover the colored COHP line
						call linwid(1) !Recover to default line width of dislin
						call endgrf
					end if
                end if
            
				!Draw lines at bottom
                !If current unit is eV (iunitx==2), the T/PDOSliney (when degeneracy is not considered) must be multiplied by au2eV, &
                !so that full height corresponds to 1.0, namely default str() array. While T/PDOSliney with consideration of degeneracy comes &
                !from compfrag() array, which is not affected by choice of unit of X-axis
				if (ilinebottom==1) then
					call AXSLEN(lenx,leny_line)
					call AXSPOS(350,iyaxspos_line)
					if (iunitx==1) CALL NAME('Energy (a.u.)','X')
					if (iunitx==2) CALL NAME('Energy (eV)','X')
					call height(ticksize)
					call color('WHITE')
					if (idegen==0) then !Do not consider degeneracy
						call setgrf('NAME','LINE','NONE','LINE')
						CALL GRAF(enelow,enehigh,enelow,stepx, 0D0,1.2D0,0D0,0.5D0)
						if (ishowTDOSline==1) then
                            if (iunitx==1) CALL CURVE(linexpos,TDOSliney,3*imoend) !Draw TDOS lines
                            if (iunitx==2) CALL CURVE(linexpos,TDOSliney*au2eV,3*imoend)
							call setcolor(6) !Gray
							if (iunitx==1) CALL CURVE(linexpos,TDOSliney_unocc,3*imoend) !Draw TDOS lines for unoccupied MOs
							if (iunitx==2) CALL CURVE(linexpos,TDOSliney_unocc*au2eV,3*imoend)
						end if
					else if (idegen==1) then !Draw degeneracy of lines for TDOS and/or MO-PDOS
						degenmax=1
						do imo=1,nmo !Get maximum degeneracy in the plotting energy range
							if (MOene_dos(imo)>enelow.and.MOene_dos(imo)<enehigh) then
								if (compfrag(imo,0)>degenmax) degenmax=compfrag(imo,0)
							end if
						end do
						CALL LABDIG(-1,"Y")
						CALL NAME('Degen.','Y')
						call setgrf('NAME','LINE','NONE','NAME')
		    				CALL TICKS(nint(degenmax),'Y')
						CALL GRAF(enelow,enehigh,enelow,stepx, 0D0,degenmax*1.2D0,0D0,degenmax*1.2D0)
						if (ishowTDOSline==1) then
							CALL CURVE(linexpos,TDOSliney,3*imoend) !Draw TDOS lines
							call setcolor(6) !Gray
							CALL CURVE(linexpos,TDOSliney_unocc,3*imoend) !Draw TDOS lines for unoccupied MOs
						end if
					end if
					do ifrag=1,nfragmax
						if (nfragDOS(ifrag)>0.and.ishowPDOSline(ifrag)==1) then
							call setcolor(iclrPDOS(ifrag))
							call linwid(ilinewidth) !Set to user-defined line width
							if (iunitx==1) CALL CURVE(linexpos,PDOSliney(:,ifrag),3*imoend)
							if (iunitx==2) CALL CURVE(linexpos,PDOSliney(:,ifrag)*au2eV,3*imoend)
						end if
					end do
					call endgrf
				end if
			
            !Draw LDOS
			else if (isel==10) then
				CALL LABDIG(nlabdigY_LDOS,"Y")
				CALL GRAF(enelow,enehigh,enelow,stepx, ylowerleft,yupperleft,ylowerleft,stepyleft)
                if (ishowLDOScurve==1) then
					call linwid(icurvewidth) !Set to user-defined line width
					CALL CURVE(curvexpos,LDOScurve,num1Dpoints)
                end if
                if (ishowLDOSline==1) then
					call linwid(ilinewidth) !Set to user-defined line width
					CALL CURVE(linexpos,LDOSliney,3*imoend)
                end if
				call linwid(1) !Recover to default line width of dislin
				call endgrf !Finished drawing LDOS
            end if

			call disfin
			if (isavepic==1) write(*,*) "Image file has been saved to current folder with ""dislin"" prefix"
		end if
		idraw=0
		
		!============= Post-processing menu ==============
		!============= Post-processing menu ==============
		!============= Post-processing menu ==============
		write(*,*)
		write(*,*) "                    -------- Post-processing menu --------"
		if (isel==0) then !T/P/OPDOS/COHP
			if (iCOHP==1) write(*,*) "-2 Show COHP raw data and ICOHP"
            write(*,*) "-1 Set format of saved image file, current: "//graphformat
			write(*,*) "0 Return"
			write(*,*) "1 Show graph again"
			write(*,*) "2 Save the graph to image file in current folder"
			write(*,*) "3 Export curve and line data to plain text file in current folder"
			if (idoPDOS==1) then
				write(*,"(' 4 Set Y-axis range and step for TDOS+PDOS, current:',f8.2,' to',f8.2,', step:',f6.2)") ylowerleft,yupperleft,stepyleft
			else
				if (iCOHP==0) then
					write(*,"(' 4 Set Y-axis range and step for TDOS, current:',f8.2,' to',f8.2,', step:',f6.2)") ylowerleft,yupperleft,stepyleft
                else
					if (abs(ylowerleft)<10000.and.abs(yupperleft)<10000) then
						write(*,"(' 4 Set left Y-axis range and step, current:',f9.3,' to',f9.3,', step:',f8.3)") ylowerleft,yupperleft,stepyleft
					else
						write(*,"(' 4 Set left Y-axis range and step, current:',1PE10.3,' to',1PE10.3,', step:',1PE10.3)") ylowerleft,yupperleft,stepyleft
					end if
					if (abs(ylowerright)<10000.and.abs(yupperright)<10000) then
						write(*,"('-4 Set right Y-axis range and step, current:',f9.3,' to',f9.3,', step:',f8.3)") ylowerright,yupperright,stepyright
					else
						write(*,"('-4 Set right Y-axis range and step, current:',1PE10.3,' to',1PE10.3,', step:',1PE10.3)") ylowerright,yupperright,stepyright
					end if
                end if
			end if
            if (iCOHP==0) then
				if (ishowTDOScurve==1) write(*,*) "5 Toggle showing TDOS curve, current: Yes"
				if (ishowTDOScurve==0) write(*,*) "5 Toggle showing TDOS curve, current: No"
				if (ishowTDOSline==1) write(*,*) "6 Toggle showing TDOS discrete lines, current: Yes"
				if (ishowTDOSline==0) write(*,*) "6 Toggle showing TDOS discrete lines, current: No"
            else
				if (ishowCOHPline==1) write(*,*) "6 Toggle showing COHP discrete lines, current: Yes"
				if (ishowCOHPline==0) write(*,*) "6 Toggle showing COHP discrete lines, current: No"
            end if
			if (idoPDOS==1) then
				write(*,*) "7 Toggle showing PDOS curves"
				write(*,*) "8 Toggle showing PDOS discrete lines"
			end if
			if (idoOPDOS==1) then
				if (ishowOPDOScurve==1) write(*,*) "9 Toggle showing OPDOS curve, current: Yes"
				if (ishowOPDOScurve==0) write(*,*) "9 Toggle showing OPDOS curve, current: No"
				if (ishowOPDOSline==1) write(*,*) "10 Toggle showing OPDOS lines, current: Yes"
				if (ishowOPDOSline==0) write(*,*) "10 Toggle showing OPDOS lines, current: No"
			end if
			if (idoPDOS==1) write(*,*) "11 Set color for PDOS curves and discrete lines"
            if (iCOHP==0) then
				if (ishowlegend==1) write(*,"(a,' X=',i5,', Y=',i5)") " 12 Set position of legends, current:",legendx,legendy
				if (ishowlegend==1) write(*,*) "13 Toggle showing legends, current: Yes"
				if (ishowlegend==0) write(*,*) "13 Toggle showing legends, current: No"
            end if
			if (idoOPDOS==1) write(*,"(a,f10.5)") " 14 Set scale factor of Y-axis range for OPDOS, current:",Yrightsclfac
			if (ishowHOMOlev==0) write(*,*) "15 Toggle showing a dashed line to highlight HOMO level, current: No"
			if (ishowHOMOlev==1) write(*,*) "15 Toggle showing a dashed line to highlight HOMO level, current: Yes"
			if (iCOHP==0) write(*,*) "16 Set the texts in the legends"
			write(*,"(a,i3)") " 17 Set width of curves, current:",icurvewidth
			write(*,"(a,i3)") " 18 Set width of discrete lines, current:",ilinewidth
		    if (ishowYlab==1) write(*,*) "19 Toggle showing labels and ticks on Y-axis, current: Yes"
		    if (ishowYlab==0) write(*,*) "19 Toggle showing labels and ticks on Y-axis, current: No"
            write(*,*) "20 Set number of decimal places for axes"
            write(*,*) "21 Set text sizes"
            if (iCOHP==0) then
				if (ilinebottom==0) write(*,*) "22 Toggle drawing discrete lines at bottom of curves, current: No"
				if (ilinebottom==1) write(*,*) "22 Toggle drawing discrete lines at bottom of curves, current: Yes"
            end if
            if (ishowvertline==0) write(*,*) "23 Toggle drawing a vertical line at a specific X position, current: No"
            if (ishowvertline==1) write(*,*) "23 Toggle drawing a vertical line at a specific X position, current: Yes"
			read(*,*) isel2

            if (isel2==-2) then
				sumocc=0
                sumunocc=0
				write(*,"(a,/)") " Note: Only the orbitals involved in present plot and with non-zero value are shown below"
				write(*,*) "   Orb.   Occ.      Ene.(eV)    Spin     -COHP(eV)"
				do imo=1,imoend
					if (OPfrag12(imo)==0) cycle
					irealmo=imo
					if (ispin==2) irealmo=imo+nbasis
                    write(*,"(i7,f8.3,f14.4,5x,a,f14.6)") imo,MOocc_dos(irealmo),MOene_dos(irealmo),orbtypename_short(MOtype(irealmo)),-OPfrag12(imo)
                    if (MOocc_dos(irealmo)==0) then
						sumunocc=sumunocc+OPfrag12(imo)
                    else
						sumocc=sumocc+OPfrag12(imo)
                    end if
                end do
                write(*,*)
                write(*,"(a,f16.4,' eV')") " Sum of COHP of occupied levels:  ",sumocc
                write(*,"(a,f16.4,' eV')") " Sum of COHP of unoccupied levels:",sumunocc
                write(*,"(a,f16.4,' eV')") " Sum of COHP of all levels:       ",sumunocc+sumocc
                write(*,*) "Note: Sum of COHP of occupied levels can be regarded as ICOHP"
            else if (isel2==-1) then
                call setgraphformat
			else if (isel2==0) then
				exit
			else if (isel2==1) then
				idraw=1
				isavepic=0
			else if (isel2==2) then
				idraw=1
				isavepic=1
			else if (isel2==3) then
				if (iCOHP==0) then
					open(10,file="DOS_line.txt",status="replace")
					if (idoOPDOS==1) then
						do i=1,3*imoend
							write(10,"(f12.5,12f12.6)") linexpos(i),TDOSliney(i),OPDOSliney(i),PDOSliney(i,:)
						end do
					else if (idoPDOS==1) then
						do i=1,3*imoend
							write(10,"(f12.5,11f12.6)") linexpos(i),TDOSliney(i),PDOSliney(i,:)
						end do
					else
						do i=1,3*imoend
							write(10,"(f12.5,f12.6)") linexpos(i),TDOSliney(i)
						end do
					end if
					close(10)
					open(10,file="DOS_curve.txt",status="replace")
					if (idoOPDOS==1) then
						do i=1,num1Dpoints
							write(10,"(f10.5,12f12.6)") curvexpos(i),TDOScurve(i),OPDOScurve(i),PDOScurve(i,:)
						end do
					else if (idoPDOS==1) then
						do i=1,num1Dpoints
							write(10,"(f10.5,11f12.6)") curvexpos(i),TDOScurve(i),PDOScurve(i,:)
						end do
					else
						do i=1,num1Dpoints
							write(10,"(f10.5,f12.6)") curvexpos(i),TDOScurve(i)
						end do
					end if
					close(10)
					write(*,*) "Curve data have been written to DOS_curve.txt in current folder"
					write(*,*) "Discrete line data have been written to DOS_line.txt in current folder"
					write(*,*) "Column 1: Energy ("//trim(unitstr)//")"
					write(*,*) "Column 2: TDOS"
					if (idoOPDOS==1) then
						write(*,*) "Column 3: OPDOS"
						write(*,*) "Column 4-13: PDOS of fragment 1-10"
					else if (idoPDOS==1) then
						write(*,*) "Column 3-12: PDOS of fragment 1-10"
					end if
                else !COHP
					open(10,file="COHP_line.txt",status="replace")
					do i=1,3*imoend
						write(10,"(f12.5,f20.6)") linexpos(i),-COHPliney(i)
					end do
					close(10)
					open(10,file="COHP_curve.txt",status="replace")
					do i=1,num1Dpoints
						write(10,"(f10.5,f20.6)") curvexpos(i),-COHPcurve(i)
					end do
					close(10)
					write(*,*) "Curve data have been written to COHP_curve.txt in current folder"
					write(*,*) "Discrete line data have been written to COHP_line.txt in current folder"
					write(*,*) "Column 1: Energy ("//trim(unitstr)//")"
					write(*,*) "Column 2: -COHP (curve in dimensionless, line in eV)"
                end if
			else if (isel2==4) then
				write(*,*) "Input lower and upper limit as well as stepsize, e.g. 0.0,2.4,0.3"
				read(*,*) ylowerleft,yupperleft,stepyleft
				iusersetY=1
                if (iCOHP==0) then
					ylowerright=ylowerleft*Yrightsclfac !Lower and upper limit for OPDOS. Set it here aims for immediately make effect
					yupperright=yupperleft*Yrightsclfac
                end if
			else if (isel2==-4) then !Only involved in the case of COHP
				write(*,*) "Input lower and upper limit as well as stepsize, e.g. 0.0,2.4,0.3"
				read(*,*) ylowerright,yupperright,stepyright
				iusersetY=1
			else if (isel2==5) then
				if (ishowTDOScurve==0) then
					ishowTDOScurve=1
				else
					ishowTDOScurve=0
				end if
			else if (isel2==6) then
				if (iCOHP==0) then
					if (ishowTDOSline==0) then
						ishowTDOSline=1
					else
						ishowTDOSline=0
					end if
                else
					if (ishowCOHPline==0) then
						ishowCOHPline=1
					else
						ishowCOHPline=0
					end if
                end if
			else if (isel2==7) then
				do while(.true.)
					write(*,*)
					write(*,*) "0 Return"
					do ifrag=1,nfragmax
						if (nfragDOS(ifrag)>0) then
							if (ishowPDOScurve(ifrag)==1) write(*,"(i2,' Toggle showing PDOS curve of fragment',i3,', current: Yes')") ifrag,ifrag
							if (ishowPDOScurve(ifrag)==0) write(*,"(i2,' Toggle showing PDOS curve of fragment',i3,', current: No')") ifrag,ifrag
						end if
					end do
					read(*,*) iselfrag
					if (iselfrag==0) exit
					if (ishowPDOScurve(iselfrag)==0) then
						ishowPDOScurve(iselfrag)=1
					else
						ishowPDOScurve(iselfrag)=0
					end if
				end do
			else if (isel2==8) then
				do while(.true.)
					write(*,*)
					write(*,*) "0 Return"
					do ifrag=1,nfragmax
						if (nfragDOS(ifrag)>0) then
							if (ishowPDOSline(ifrag)==1) write(*,"(i2,' Toggle showing PDOS line of fragment',i3,', current: Yes')") ifrag,ifrag
							if (ishowPDOSline(ifrag)==0) write(*,"(i2,' Toggle showing PDOS line of fragment',i3,', current: No')") ifrag,ifrag
						end if
					end do
					read(*,*) iselfrag
					if (iselfrag==0) exit
					if (ishowPDOSline(iselfrag)==0) then
						ishowPDOSline(iselfrag)=1
					else
						ishowPDOSline(iselfrag)=0
					end if
				end do
			else if (isel2==9) then
				if (ishowOPDOScurve==0) then
					ishowOPDOScurve=1
				else
					ishowOPDOScurve=0
				end if
			else if (isel2==10) then
				if (ishowOPDOSline==0) then
					ishowOPDOSline=1
				else
					ishowOPDOSline=0
				end if
			else if (isel2==11) then
				do while(.true.)
					write(*,*)
					write(*,*) "0 Return"
					do ifrag=1,nfragmax
						if (nfragDOS(ifrag)>0) write(*,"(i2,' Set color for fragment',i3,', current: ',a)") ifrag,ifrag,colorname(iclrPDOS(ifrag))
					end do
					read(*,*) iselfrag
					if (iselfrag<0.or.iselfrag>nfragmax) then
						write(*,*) "Error: Index exceeded valid range!"
					else if (iselfrag==0) then
						exit
					else
						call selcolor(iclrPDOS(iselfrag))
					end if
				end do
            else if (isel2==12) then
                write(*,*) "Input X position of the legends, e.g. 400"
                write(*,*) "The larger the value, the more the position of the legends is right"
                write(*,"(a,i6,a)") " If directly press ENTER button, current value",legendx," will be kept"
                read(*,"(a)") c80tmp
                if (c80tmp/=" ") read(c80tmp,*) legendx
                write(*,*) "Input Y position of the legends, e.g. 160"
                write(*,*) "The larger the value, the lower the position of the legends"
                write(*,"(a,i6,a)") " If directly ENTER button, current value",legendy," will be kept"
                read(*,"(a)") c80tmp
                if (c80tmp/=" ") read(c80tmp,*) legendy
			else if (isel2==13) then
				if (ishowlegend==0) then
					ishowlegend=1
				else
					ishowlegend=0
				end if
			else if (isel2==14) then
				write(*,*) "Input scale factor, e.g. 0.15"
				read(*,*) Yrightsclfac
				ylowerright=ylowerleft*Yrightsclfac !Lower and upper limit for OPDOS. Set it here aims for immediately make effect
				yupperright=yupperleft*Yrightsclfac
			else if (isel2==15) then
				if (ishowHOMOlev==0) then
					ishowHOMOlev=1
				else
					ishowHOMOlev=0
				end if
			else if (isel2==16) then
				do while(.true.)
					write(*,*) "Select a term to change its legend, e.g. 3"
					write(*,"(' -2 TDOS, current text: ',a)") trim(TDOSstring)
					if (idoOPDOS==1) write(*,"(' -1 OPDOS, current text: ',a)") trim(OPDOSstring)
					write(*,*) " 0 Return"
					do ifrag=1,nfragmax
						if (nfragDOS(ifrag)>0) write(*,"(i3,' PDOS of frag',i2,', current text: ',a)") ifrag,ifrag,trim(PDOSstring(ifrag))
					end do
					read(*,*) iseltmp
					if (iseltmp==0) exit
					write(*,*) "Input text for the legend"
					read(*,"(a)") c80tmp
					if (iseltmp==-2) then
						TDOSstring=c80tmp
					else if (iseltmp==-1) then
						OPDOSstring=c80tmp
					else
						PDOSstring(iseltmp)=c80tmp
					end if
				end do
			else if (isel2==17) then
				write(*,*) "Input width of curves, e.g. 5"
				read(*,*) icurvewidth
			else if (isel2==18) then
				write(*,*) "Input width of discrete lines, e.g. 3"
				read(*,*) ilinewidth
            else if (isel2==19) then
                if (ishowYlab==0) then
                    ishowYlab=1
                else
                    ishowYlab=0
                end if
            else if (isel2==20) then
                do while(.true.)
                    write(*,*)
                    write(*,*) "0 Return"
                    write(*,"(a,i3)") " 1 Set decimal places for X axis, current:",nlabdigX
                    if (iCOHP==0) then
						write(*,"(a,i3)") " 2 Set decimal places for Y axis of TDOS/PDOS, current:",nlabdigY
						write(*,"(a,i3)") " 3 Set decimal places for Y axis of OPDOS, current:",nlabdigY_OPDOS
						write(*,"(a,i3)") " 4 Set decimal places for Y axis of LDOS, current:",nlabdigY_LDOS
                    else
						write(*,"(a,i3)") " 5 Set decimal places for left Y axis of COHP, current:",nlabdigleftY_COHP
						write(*,"(a,i3)") " 6 Set decimal places for right Y axis of COHP, current:",nlabdigrightY_COHP
                    end if
                    read(*,*) itmp
                    if (itmp==0) then
                        exit
                    else
                        write(*,*) "Input the number of decimal places, e.g. 3"
                        if (itmp==1) read(*,*) nlabdigX
                        if (itmp==2) read(*,*) nlabdigY
                        if (itmp==3) read(*,*) nlabdigY_OPDOS
                        if (itmp==4) read(*,*) nlabdigY_LDOS
                        if (itmp==5) read(*,*) nlabdigleftY_COHP
                        if (itmp==6) read(*,*) nlabdigrightY_COHP
                    end if
                end do
            else if (isel2==21) then !Set text sizes
                do while(.true.)
                    write(*,*)
                    write(*,*) "0 Return"
                    write(*,"(a,i3)") " 1 Set text size of name of axes, current:",height_axis
                    write(*,"(a,i3)") " 2 Set text size of ticks, current:",ticksize
                    write(*,"(a,i3)") " 3 Set text size of legend, current:",legtextsize
                    read(*,*) itmp
                    if (itmp==0) then
                        exit
                    else if (itmp==1) then
                        write(*,*) "Input size, e.g. 45"
                        read(*,*) height_axis
                    else if (itmp==2) then
                        write(*,*) "Input size, e.g. 45"
                        read(*,*) ticksize
                    else if (itmp==3) then
                        write(*,*) "Input size, e.g. 40"
                        read(*,*) legtextsize
                    end if
                end do
            else if (isel2==22) then
				if (ilinebottom==1) then
					ilinebottom=0
                else
					ilinebottom=1
                end if
            else if (isel2==23) then
				write(*,*) "Input the X-position to plot a vertical line, e.g. -2.8"
                write(*,*) "If inputting ""q"", then no vertical line will be plotted"
                read(*,*) c80tmp
                if (index(c80tmp,'q')/=0) then
					ishowvertline=0
                else
					read(c80tmp,*) vertline_X
					ishowvertline=1
                end if
			end if
			
		else if (isel==10) then !LDOS in 1D
	    	write(*,*) "-1 Set format of saved image file, current: "//graphformat
			write(*,*) "0 Return"
			write(*,*) "1 Show graph again"
			write(*,*) "2 Save graphical file of the LDOS map in current folder"
			write(*,*) "3 Export curve data to plain text file in current folder"
			write(*,"(' 4 Set Y-axis range and step, current:',f9.3,' to',f9.3,', step:',f7.3)") ylowerleft,yupperleft,stepyleft
			if (ishowLDOSline==1) write(*,*) "6 Toggle showing LDOS line, current: Yes"
			if (ishowLDOSline==0) write(*,*) "6 Toggle showing LDOS line, current: No"
			write(*,"(a,i3)") " 7 Set width of curves, current:",icurvewidth
			write(*,"(a,i3)") " 8 Set width of lines, current:",ilinewidth
		    if (ishowYlab==1) write(*,*) "9 Toggle showing labels and ticks on Y-axis, current: Yes"
		    if (ishowYlab==0) write(*,*) "9 Toggle showing labels and ticks on Y-axis, current: No"
			read(*,*) isel2
            
            if (isel2==-1) then
                call setgraphformat
			else if (isel2==0) then
				iusersetY=0
				exit
			else if (isel2==1) then
				idraw=1
				isavepic=0
			else if (isel2==2) then
				idraw=1
				isavepic=1
			else if (isel2==3) then
				open(10,file="LDOS_line.txt",status="replace")
				do i=1,3*imoend
					write(10,"(f10.5,1PE15.6)") linexpos(i),LDOSliney(i)
				end do
				close(10)
				open(10,file="LDOS_curve.txt",status="replace")
				do i=1,num1Dpoints
					write(10,"(f10.5,1PE15.6)") curvexpos(i),LDOScurve(i)
				end do
				close(10)
				write(*,*) "Curve data have been written to LDOS_curve.txt in current folder"
				write(*,*) "Discrete line data have been written to LDOS_line.txt in current folder"
				write(*,*) "Column 1: Energy ("//trim(unitstr)//")"
				write(*,*) "Column 2: LDOS"
			else if (isel2==4) then
				write(*,*) "Input lower and upper limit as well as stepsize, e.g. 0.0,2.4,0.3"
				read(*,*) ylowerleft,yupperleft,stepyleft
				iusersetY=1
			else if (isel2==6) then
				if (ishowLDOSline==0) then
					ishowLDOSline=1
				else
					ishowLDOSline=0
				end if
			else if (isel2==7) then
				write(*,*) "Input width of curves, e.g. 5"
				read(*,*) icurvewidth
			else if (isel2==8) then
				write(*,*) "Input width of discrete lines, e.g. 3"
				read(*,*) ilinewidth
            else if (isel2==9) then
                if (ishowYlab==0) then
                    ishowYlab=1
                else
                    ishowYlab=0
                end if
			end if
		end if
	end do


!!================================================================
!!====== Plot local DOS along a line (2D color-filled map) =======
!!================================================================
else if (isel==11) then
	if (allocated(LDOS2Dmap)) deallocate(LDOS2Dmap,LDOSptscomp)
	write(*,*) "Input the starting point (in Bohr), e.g. 1.0,0,0.2"
	read(*,*) vorgx,vorgy,vorgz
	write(*,*) "Input the end point (in Bohr), e.g. 2.0,0,0.4"
	read(*,*) vendx,vendy,vendz
	write(*,*) "Input the number of points along the line, e.g. 200"
	read(*,*) numLDOSpt
	allocate(LDOS2Dmap(num2Dpoints,numLDOSpt),LDOSptscomp(nmo,numLDOSpt))
	write(*,*) "Calculating orbital composition, please wait..."
	xlen=vendx-vorgx
	vx=xlen/(numLDOSpt-1)
	ylen=vendy-vorgy
	vy=ylen/(numLDOSpt-1)
	zlen=vendz-vorgz
	vz=zlen/(numLDOSpt-1)
	totlen=dsqrt(xlen**2+ylen**2+zlen**2)
	dlen=dsqrt(vx**2+vy**2+vz**2)
	LDOSptscomp=0D0
	FWHMmax=maxval(FWHM)
	do ipt=1,numLDOSpt
		x=vorgx+vx*(ipt-1)
		y=vorgy+vy*(ipt-1)
		z=vorgz+vz*(ipt-1)
		do imo=1,imoend
			imoall=imo
			if (ispin==2) imoall=imo+nbasis
			if (MOene_dos(imoall)<enelow-3*FWHMmax.or.MOene_dos(imoall)>enehigh+3*FWHMmax) cycle
			LDOSptscomp(imo,ipt)=fmo(x,y,z,imoall)**2
		end do
	end do
	enestep=(enehigh-enelow)/(num2Dpoints-1) 
	do i=1,num2Dpoints
		LDOSxpos(i)=enelow+(i-1)*enestep
	end do
	
	LDOS2Dmap=0D0
	if (ibroadfunc==1.or.ibroadfunc==3) then !Lorentzian function, see http://mathworld.wolfram.com/LorentzianFunction.html
		do ilinept=1,numLDOSpt !Cycle each point in the line
			do imo=1,imoend !Cycle each orbital
				irealmo=imo
				if (ispin==2) irealmo=imo+nbasis
				preterm=str(imo)*0.5D0/pi*FWHM(imo)
				do imappt=1,num2Dpoints !Broaden imo as curve
					tmp=preterm/( (LDOSxpos(imappt)-MOene_dos(irealmo))**2+0.25D0*FWHM(imo)**2 )
					LDOS2Dmap(imappt,ilinept)=LDOS2Dmap(imappt,ilinept)+tmp*LDOSptscomp(imo,ilinept)
				end do
			end do
		end do
	end if
	if (ibroadfunc==2.or.ibroadfunc==3) then !Gaussian function, see http://en.wikipedia.org/wiki/Gaussian_function
		if (ibroadfunc==3) LDOS2Dmap=(1-gauweigh)*LDOS2Dmap
		do ilinept=1,numLDOSpt !Cycle each point in the line
			do imo=1,imoend !Cycle each orbital
				irealmo=imo
				if (ispin==2) irealmo=imo+nbasis
				gauss_c=FWHM(imo)/2D0/sqrt(2*dlog(2D0))
				gauss_a=str(imo)/(gauss_c*sqrt(2D0*pi))
				do imappt=1,num2Dpoints !Broaden imo as curve
					tmp=gauss_a*dexp( -(LDOSxpos(imappt)-MOene_dos(irealmo))**2/(2*gauss_c**2) )
					if (ibroadfunc==3) tmp=gauweigh*tmp !Combine Lorentizan and Gaussian function
					LDOS2Dmap(imappt,ilinept)=LDOS2Dmap(imappt,ilinept)+tmp*LDOSptscomp(imo,ilinept)
				end do
			end do
		end do
	end if
	LDOS2Dmap=LDOS2Dmap*scalecurve
	write(*,*)
	write(*,"(' Energy range: from',f12.5,' to',f12.5,1x,a)") enelow,enehigh,trim(unitstr)
	write(*,"(' Starting point:',3f12.6,' Bohr')") vorgx,vorgy,vorgz
	write(*,"(' End point:     ',3f12.6,' Bohr')") vendx,vendy,vendz
	write(*,"(' Stepsize:',f12.6,' Bohr, total length:',f12.6,' Bohr')") dlen,totlen
	
	idraw=1
	isavepic=0
    Yleftstep=totlen/10
	do while(.true.)
		if (isilent==0.and.idraw==1) then
	! 		call drawmatcolor(LDOS2Dmap,num2Dpoints,numLDOSpt,enelow,enehigh,0D0,totlen,0D0,maxval(LDOS2Dmap),(enehigh-enelow)/10,totlen/10,fillcoloritpx,2)
			lengthx=2300
			call SCRMOD('REVERSE')
			CALL setxid(0,'NONE')
			CALL PAGE(3200,2700)
			if (isavepic==0) then
				call METAFL('xwin')
				call window(200,100,900,770)
			else if (isavepic==1) then
				call METAFL(graphformat)
				call winsiz(graph2Dwidth,graph2Dheight)
				CALL IMGFMT('RGB')
			end if
			if (iusersetcolorscale==0) then
				clrscllow=0D0
				clrsclhigh=maxval(LDOS2Dmap)
				clrsclstep=maxval(LDOS2Dmap)/10
			end if
			call DISINI
			if (isavepic==0.and.isys==1) then
				call height(60)
				CALL HNAME(50)
			else !The text shown in graphic file is strangely larger than window, so slight decrease it
				call height(50)
				CALL HNAME(45)
			end if
			call ERRMOD("ALL","OFF") !If don't set this, when atom label in contour map is out of page range, DISLIN annoys users
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
			CALL LABDIG(2,"X")
			CALL LABDIG(3,"Y")
			CALL LABDIG(4,"Z")
			call height(50)
			call ticks(1,"XYZ")
			call WINTIT("Local DOS filled-color map")
			if (iunitx==1) CALL NAME('Energy (a.u.)','X')
			if (iunitx==2) CALL NAME('Energy (eV)','X')
			CALL NAME('Position in the line (Bohr)','Y')
			call center
			call AUTRES(num2Dpoints,numLDOSpt)
! 			call AX3LEN(lengthx,nint(lengthx*dfloat(numLDOSpt)/num2Dpoints),nint(lengthx*dfloat(numLDOSpt)/num2Dpoints))
			call AX3LEN(lengthx,nint(lengthx*yxratio),nint(lengthx*yxratio))
			call GRAF3(enelow,enehigh,enelow,(enehigh-enelow)/10,0D0,totlen,0D0,Yleftstep,clrscllow,clrsclhigh,clrscllow,clrsclstep)
			call CRVMAT(LDOS2Dmap,num2Dpoints,numLDOSpt,fillcoloritpx,fillcoloritpy)
			call DISFIN
		end if
		
		idraw=0
		isavepic=0
		write(*,*)
		write(*,*) "                  -------- Post-processing menu --------"
		write(*,*) "0 Return"
		write(*,*) "1 Show graph again"
		write(*,*) "2 Save the graph to image file in current folder"
		write(*,*) "3 Export curve data to plain text file in current folder"
		write(*,"(a,f8.3)") " 4 Set ratio of Y-axis relative to X-axis, current:",yxratio
		write(*,"(' 5 Set range and step of color scale, current:',f8.4,' to',f8.4,', step:',f7.4)") clrscllow,clrsclhigh,clrsclstep
        write(*,"(a,f6.3,' Bohr')") " 6 Set step of left Y-axis, current:",Yleftstep
		read(*,*) isel2

		if (isel2==0) then
			iusersetcolorscale=0
			exit
		else if (isel2==1) then
			idraw=1
			isavepic=0
		else if (isel2==2) then
			idraw=1
			isavepic=1
            write(*,*) "Done! The image has been saved to current folder with ""DISLIN"" prefix"
		else if (isel2==3) then
			open(10,file="LDOS.txt",status="replace")
			do imappt=1,num2Dpoints
				do ipt=1,numLDOSpt
					write(10,"(f8.3,f10.4,1PE16.8)") LDOSxpos(imappt),dlen*(ipt-1),LDOS2Dmap(imappt,ipt)
				end do
			end do
			close(10)
			write(*,*) "LDOS data have been written to LDOS.txt in current folder"
			write(*,*) "Column 1: Energy ("//trim(unitstr)//")"
			write(*,*) "Column 2: Coordinate in the line (Bohr)"
			write(*,*) "Column 3: LDOS value"
		else if (isel2==4) then
			write(*,*) "Input the ratio, e.g. 1.4"
			read(*,*) yxratio
		else if (isel2==5) then
			write(*,*) "Input lower and upper limits as well as stepsize of color scale"
			write(*,*) "e.g. 0.0,0.6,0.05"
			read(*,*) clrscllow,clrsclhigh,clrsclstep
			iusersetcolorscale=1
		else if (isel2==6) then
            write(*,*) "Input the step, e.g. 0.4"
            read(*,*) Yleftstep
		end if
	end do





!!=========================================
!!====== Plot photoelectron spectra =======
!!=========================================
!This is a highly independent module, all plotting parameters are irrelvant to the main interface of DOS
!This module only use eV as unit
else if (isel==12) then
	write(*,"(a)") " Note: Most plotting parameters of PES module are irrelevant to those in the DOS module. &
	&Only eV and Gaussian broadening function will be used for PES plotting. Spin types are not distinguished"
	nmoocc=count(MOocc/=0)
	if (allocated(bindene)) deallocate(bindene,PESlinex,PESliney,PES_str,PES_FWHM) !Binding energy
	allocate(PESlinex(3*nmoocc),PESliney(3*nmoocc),bindene(nmoocc),PES_str(nmoocc),PES_FWHM(nmoocc))
	!Only consider occupied MOs
	nmoocc=0
	do imo=1,nmo
		if (MOocc(imo)/=0) then
			nmoocc=nmoocc+1
			bindene(nmoocc)=MOene(imo)*au2eV
            PES_str(nmoocc)=str(imo)
            PES_FWHM(nmoocc)=FWHM(imo)*au2eV
		end if
	end do
    write(*,"(' HOMO level:',f18.4,' eV')") maxval(bindene)
	bindene=-bindene !Convert to binding energy
    if (maxval(PES_FWHM)==minval(PES_FWHM)) PES_FWHM=0.2D0
    
	do while(.true.)
		write(*,*)
		write(*,*) "                  -------- Photoelectron spectra (PES) --------"
		write(*,*) "-4 Set format of saved image file, current: "//graphformat
		write(*,*) "-3 Export occupied MO energies, strengths and FWHMs to plain text file"
		write(*,*) "-2 Show all binding energy level information"
		write(*,*) "-1 Export curve and line data to plain text file in current folder"
		write(*,*) "0 Return to last menu"
		write(*,*) "1 Show photoelectron spectra now!"
		write(*,*) "2 Save the graph to image file in current folder"
		write(*,"(a,f8.3,' eV')") " 3 Set shift value for binding energy, current:",PES_shift
		write(*,"(' 4 Set X-axis range and step, current:',f8.2,' to',f8.2,', step:',f6.2,' eV')") PES_Xlow,PES_Xhigh,PES_Xstep
		if (iusersetPES_Y==0) write(*,"(' 5 Set Y-axis range and step, current: Auto')")
		if (iusersetPES_Y==1) write(*,"(' 5 Set Y-axis range and step, current:',f8.2,' to',f8.2,', step:',f6.2)") PES_Ylow,PES_Yhigh,PES_Ystep
        if (maxval(PES_FWHM)==minval(PES_FWHM)) then
		    write(*,"(a,f8.3,' eV')") " 6 Set FWHM, current:",PES_FWHM(1)
        else
		    write(*,"(a)") " 6 Set FWHM, current: Orbital dependent"
        end if
		if (ishowPEScurve==1) write(*,*) "7 Toggle showing PES curve, current: Yes"
		if (ishowPEScurve==0) write(*,*) "7 Toggle showing PES curve, current: No"
		if (ishowPESline==1) write(*,*) "8 Toggle showing PES line, current: Yes"
		if (ishowPESline==0) write(*,*) "8 Toggle showing PES line, current: No"
		write(*,"(a,i3)") " 9 Set width of curves, current:",icurvewidth
		write(*,"(a,i3)") " 10 Set width of lines, current:",ilinewidth
		write(*,"(a,f10.5)") " 11 Set scale ratio for PES curve, current:",scalePEScurve
		if (invPES_X==1) write(*,*) "12 Toggle inverting X-axis, current: Yes"
		if (invPES_X==0) write(*,*) "12 Toggle inverting X-axis, current: No"
		if (ishowYlab==1) write(*,*) "13 Toggle showing labels and ticks on Y-axis, current: Yes"
		if (ishowYlab==0) write(*,*) "13 Toggle showing labels and ticks on Y-axis, current: No"
		read(*,*) isel2
        
        if (isel2==-4) then
            call setgraphformat
        else if (isel2==-3) then
	        open(10,file="PESinfo.txt",status="replace")
	        write(10,"(2i6)") nmoocc,4
	        do imo=1,nmoocc
		        write(10,"(f18.6,3f12.6)") -bindene(imo),1D0,PES_str(imo),PES_FWHM(imo)
	        end do
	        close(10)
	        write(*,"(a)") " The occupied MO energies (original), occupation numbers, strengths, FWHMs have been exported to &
            &PESinfo.txt in current directory, you can modify it and then load it into Multiwfn again"
	        write(*,*) "Note: The unit of energy levels and FWHMs in this file is eV"
	        write(*,*)
        else if (isel2==-2) then !The index is not necessary identical to original MOs
			do imo=1,nmoocc
				write(*,"(i6,'  Str:',f6.2,'  BE (org):',f12.3,' (shifted):',f12.3,'  FWHM:',f5.2,' eV')") &
                imo,PES_str(imo),bindene(imo),bindene(imo)+PES_shift,PES_FWHM(imo)
			end do
		else if (isel2==0) then
			exit
		else if (isel2==3) then
			write(*,*) "Input shift value for binding energies (in eV), e.g. 0.32"
			write(*,*) "Hint: To meet generalized Koopmans theorem, this value should be E(HOMO) + VIP"
			read(*,*) PES_shift
		else if (isel2==4) then
			write(*,*) "Input lower limit, upper limit and stepsize (in eV), e.g. 0.0,5.0,0.5"
			read(*,*) PES_Xlow,PES_Xhigh,PES_Xstep
		else if (isel2==5) then
			write(*,*) "Input lower limit, upper limit and stepsize, e.g. 0.0,2.4,0.3"
			read(*,*) PES_Ylow,PES_Yhigh,PES_Ystep
			iusersetPES_Y=1
		else if (isel2==6) then
			write(*,*) "Input FWHM (in eV), e.g. 0.5"
			read(*,*) tmp
            PES_FWHM(:)=tmp
		else if (isel2==7) then
			if (ishowPEScurve==0) then
				ishowPEScurve=1
			else
				ishowPEScurve=0
			end if
		else if (isel2==8) then
			if (ishowPESline==0) then
				ishowPESline=1
			else
				ishowPESline=0
			end if
		else if (isel2==9) then
			write(*,*) "Input width of curves, e.g. 5"
			read(*,*) icurvewidth
		else if (isel2==10) then
			write(*,*) "Input width of discrete lines, e.g. 3"
			read(*,*) ilinewidth
		else if (isel2==11) then
			write(*,*) "Input the scale ratio, e.g. 0.6"
			read(*,*) scalePEScurve
		else if (isel2==12) then
			if (invPES_X==0) then
				invPES_X=1
			else
				invPES_X=0
			end if
        else if (isel2==13) then
            if (ishowYlab==0) then
                ishowYlab=1
            else
                ishowYlab=0
            end if
		end if
		
		if (isel2==-1.or.isel2==1.or.isel2==2) then !Plot PES or export data
	        if (all(bindene>PES_Xhigh)) then
                write(*,"(a)") " Error: There is no molecular orbital in the current plotting range! &
                &You should increase upper limit of X-axis by option 4"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                cycle
            end if
        
			!Set X position of curves
			tmp=(PES_Xhigh-PES_Xlow)/(num1Dpoints-1) 
			do i=1,num1Dpoints
				curvexpos(i)=PES_Xlow+(i-1)*tmp
			end do
			
			bindene=bindene+PES_shift !Temporarily shift
	
			!Generate energy levels line
			do imo=1,nmoocc
				inow=3*(imo-1)
				PESlinex(inow+1:inow+3)=bindene(imo)
				PESliney(inow+1)=0D0
				PESliney(inow+2)=str(imo)
				PESliney(inow+3)=0D0
			end do
	
			!Generate PES curve
			PEScurve=0D0
			!Gaussian function, see http://en.wikipedia.org/wiki/Gaussian_function
			do imo=1,nmoocc !Cycle each orbital
				gauss_c=PES_FWHM(imo)/2D0/sqrt(2*dlog(2D0))
				gauss_a=PES_str(imo)/(gauss_c*sqrt(2D0*pi))
				do ipoint=1,num1Dpoints !Broaden imo as curve
					tmp=gauss_a*dexp( -(curvexpos(ipoint)-bindene(imo))**2/(2*gauss_c**2) )
					PEScurve(ipoint)=PEScurve(ipoint)+tmp
				end do
			end do
			PEScurve=PEScurve*scalePEScurve
			if (iusersetPES_Y==0) then !Y axis range was not set by user, automatically determine
				PES_Yhigh=1.1D0*maxval(PEScurve)
				PES_Ylow=0
				PES_Ystep=nint((PES_Yhigh-PES_Ylow)*10)/100D0
			end if
			
			bindene=bindene-PES_shift !Shift back
			
			if (isel2==1.or.isel2==2) then !Plot PES
				if (isel2==1) isavepic=0
				if (isel2==2) isavepic=1
				if (isavepic==0) then
					call METAFL('xwin')
					call window(200,100,1000,600)
				else if (isavepic==1) then
					call METAFL(graphformat)
					call winsiz(graph1Dwidth,graph1Dheight)
				end if
				call SCRMOD('REVERSE')
				CALL IMGFMT("RGB")
				CALL PAGE(3000,1800)
				CALL setxid(0,'NONE')
				call DISINI
                CALL HNAME(45)
				if (isavepic==0.and.isys==1) then
					call height(45)
				else
					call height(40) !The text shown in graphic file is strangely larger than window, so slight decrease it
				end if
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
				call AXSLEN(2450,1400)
				if (isavepic==0) call WINTIT("Click right mouse button to close")
                if (ishowYlab==0) then
				    call AXSPOS(310,1550)
                    call labels("NONE","Y")
				    CALL TICKS(0,'Y')
				    CALL TICKS (1,'X')
                    CALL NAMDIS(70,'Y')
                else
				    call AXSPOS(380,1550)
				    CALL TICKS (1,'XY')
                    CALL NAMDIS(50,'Y')
                end if
				call ERRMOD("ALL","OFF")
				CALL LABDIG(2,"X")
				CALL LABDIG(2,"Y")
				CALL NAME('Binding energy (eV)','X')
				CALL NAME('PES intensity (arb.)','Y')
				if (invPES_X==0) then
                    CALL GRAF(PES_Xlow,PES_Xhigh,PES_Xlow,PES_Xstep, PES_Ylow,PES_Yhigh,PES_Ylow,PES_Ystep)
				else if (invPES_X==1) then
                    CALL GRAF(PES_Xhigh,PES_Xlow,PES_Xhigh,-PES_Xstep, PES_Ylow,PES_Yhigh,PES_Ylow,PES_Ystep)
                end if
				call linwid(icurvewidth) !Set to user-defined line width
				if (ishowPEScurve==1) CALL CURVE(curvexpos,PEScurve,num1Dpoints)
				call linwid(ilinewidth) !Set to user-defined line width
				if (ishowPESline==1) CALL CURVE(PESlinex,PESliney,3*nmoocc)

				call disfin
				if (isavepic==1) write(*,*) "Graphical file has been saved to current folder with ""dislin"" prefix"
			
			else if (isel2==-1) then
				open(10,file="PES_line.txt",status="replace")
				do i=1,3*nmoocc
					write(10,"(f16.6,f12.6)") PESlinex(i),PESliney(i)
				end do
				close(10)
				open(10,file="PES_curve.txt",status="replace")
				do i=1,num1Dpoints
					write(10,"(f16.6,f12.6)") curvexpos(i),PEScurve(i)
				end do
				close(10)
				write(*,*) "Curve data have been written to PES_curve.txt in current folder"
				write(*,*) "Discrete line data have been written to PES_line.txt in current folder"
				write(*,*) "Column 1: Shifted binding energy (eV)"
				write(*,*) "Column 2: PES strength"
			end if
		end if
	end do
end if


end do !End of main loop
end subroutine