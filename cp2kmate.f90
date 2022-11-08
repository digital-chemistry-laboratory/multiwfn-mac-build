!!---------- cp2kmate is a collection of utilities of CP2K
subroutine cp2kmate
write(*,*) "1 Create CP2K input file"
write(*,*) "2 Convert band structure file of CP2K (.bs) to multiple column file"
!write(*,*) "3 Extract k-point wavefunction"
read(*,*) isel

if (isel==1) then
    call outCP2Kinp_wrapper
else if (isel==2) then
    call CP2K_BS
else if (isel==3) then
    !call CP2K_kptwfn
end if
end subroutine




!!---------- Interface of outputting CP2K input file
subroutine outCP2Kinp_wrapper
use util
use defvar
character(len=200) outname,c200tmp
call path2filename(filename,c200tmp)
write(*,"(/,a)") " Note: Please mention Multiwfn and cite original paper of Multiwfn if you benefits from this function in your study, thank you!"
write(*,*) 
write(*,*) "Input path for generating CP2K input file, e.g. C:\ltwd.inp"
write(*,"(a)") " If press ENTER button directly, will export to "//trim(c200tmp)//".inp"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".inp"
call outCP2Kinp(outname,10)
end subroutine
!!---------- Output current coordinate to CP2K input file
subroutine outCP2Kinp(outname,ifileid)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) outname
integer :: ifileid,ibas=2,tmparr(ncenter)
character selectyn,c80tmp*80,c80tmp2*80,c200tmp*200,c2000tmp*2000
character :: method*22="PBE",PBCdir*4="XYZ ",cellfix*4="NONE"
character(len=30) :: basname(-10:30)=" "
integer :: itask=1,idispcorr=0,imolden=0,ioutvibmol=1,ithermostat=0,ibarostat=0,inoSCFinfo=0,iSCCS=0,idipcorr=0,iwfc=0,iHFX=0,imoment=0,ioptmethod=1,iprintlevel=1
integer :: iTDDFT=0,nstates_TD=3,iTDtriplet=0,isTDA=0,iNTO=0,nADDED_MOS=0,icentering=0
integer :: iMDformat=1,nMDsavefreq=1,ioutcube=0,idiagOT=1,imixing=2,ismear=0,iatomcharge=0,ifineXCgrid=0,iouterSCF=1,iDFTplusU=0,NHOMO=0,NLUMO=0
integer :: natmcons=0,nthermoatm=0,ikpoint1=1,ikpoint2=1,ikpoint3=1,nrep1=1,nrep2=1,nrep3=1
integer,allocatable :: atmcons(:),thermoatm(:)
real*8 :: efieldvec(3)=0,vacsizex=5/b2a,vacsizey=5/b2a,vacsizez=5/b2a
real*8 :: frag1chg,frag2chg
integer :: frag1multi,frag2multi,totalmulti
integer :: iprestype=1,ioutSbas=0,ioutorbene=0,istate_force=1,idiaglib=1,iGAPW=0,iLSSCF=0,iLRIGPW=0,iPSOLVER=1
real*8 :: Piso=1.01325D0,Ptens(3,3)=reshape( [1.01325D0,0D0,0D0, 0D0,1.01325D0,0D0, 0D0,0D0,1.01325D0], shape=shape(Ptens))
integer :: CUTOFF=350,REL_CUTOFF=50

!Status information of current system. ",save" is used so that for the same system we can enter this interface multiple times to generate various input files
integer,save :: netchg,multispin
integer,allocatable,save :: atmkind(:) !The kind that atoms belonging to
integer,parameter :: nkindmax=200
integer,save :: nkind=0 !Current number of kinds
character(len=5),save :: kindname(nkindmax) !Name of each kind
integer,save :: kindeleidx(nkindmax) !Element idx of each kind
integer,save :: kindmag(nkindmax) !Magnetization of each kind
character,save :: lastinpname*200 !Input file of last time in this interface

!Defining array recording number of valence electrons for various elements of MOLOPT-GTH basis set, so that -q? can be automatically added to basis set name (except for 0)
!For some elements, there are small and large core versions. Here records the default one, which corresponds to large core for transition metals, and small core for others
!For CP2K 9.1, the information comes from GTH-PADE of GTH_POTENTIALS, except that Ln and Ac series come from LnPP1_POTENTIALS and AcPP1_POTENTIALS (medium-core version), respectively
integer :: Nval(0:nelesupp)=(/ 0, & !X
1, 2,                                             & !H ~He
3, 4,                                3,4,5,6,7,8, & !Li~Ne
9,10,                                3,4,5,6,7,8, & !Na~Ar
9,10,11,12,13,14,15,16,17,18,11,12, 13,4,5,6,7,8, & !K ~Kr
9,10,11,12,13,14,15,16,17,18,11,12, 13,4,5,6,7,8, & !Rb~Xe
9,10,11,12,13,14,15,16,17,18,29,30,31,32,33,34,35,& !Cs,Ba,La~Lu
        12,13,14,15,16,17,18,11,12, 13,4,5,6,7,8, & !Hf~Rn
0, 0,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,& !Fr,Ra,Ac~Lr
(0,iele=104,nelesupp) /)

iconvtest=0
ikpconvtest=0
natmcons=0 !Disable atom constraint setting if set last time
if (index(outname,"cutconv.inp")/=0) then
    iconvtest=1
    iprintlevel=2 !Medium printing level
    !ifineXCgrid=1 !I found this effectively often improves smoothness of convergence with respect to cutoff
    write(*,*) "Note: The generated file is dedicated to cutoff convergence testing purpose"
else if (index(outname,"kpconv.inp")/=0) then
    ikpconvtest=1
    write(*,*) "Note: The generated file is dedicated to k-point convergence testing purpose"
end if

!Conversion of basis set index to basis set name
basname(-1)="SZV-GTH"
basname(-2)="DZVP-GTH"
basname(-3)="TZVP-GTH"
basname(-4)="TZV2P-GTH"
basname(-5)="QZV2P-GTH"
basname(-6)="QZV3P-GTH"
basname(1)="SZV-MOLOPT-SR-GTH"
basname(2)="DZVP-MOLOPT-SR-GTH"
basname(3)="TZVP-MOLOPT-GTH"
basname(4)="TZV2P-MOLOPT-GTH"
basname(5)="TZV2PX-MOLOPT-GTH"
basname(10)="6-31G*"
basname(11)="6-311G**"
basname(12)="Ahlrichs-def2-TZVP"
basname(13)="pob-TZVP"
basname(14)="Ahlrichs-def2-QZVP"
basname(20)="cc-DZ with RI_DZ"
basname(21)="cc-TZ with RI_TZ"

!Construct kind information, charge and multiplicity
10 if (allocated(atmkind).and.filename/=lastinpname) deallocate(atmkind) !Has allocated last time, need to regenerate if the last system is different to the present one
if (.not.allocated(atmkind)) then !Haven't been constructed, namely first time use this function for present system, or this system it not identical to the last one
    write(*,*) "Generating KIND information..."
    allocate(atmkind(ncenter))
    !Build initial kind list
    nkind=1
    kindeleidx(1)=a(1)%index
    kindname(1)=a(1)%name
    do iatm=2,ncenter
        if (all(kindeleidx(:nkind)/=a(iatm)%index)) then
            nkind=nkind+1
            kindeleidx(nkind)=a(iatm)%index
            kindname(nkind)=a(iatm)%name
        end if
    end do
    !Assign each atom by a kind
    do iatm=1,ncenter
        do ikind=1,nkind
            if (a(iatm)%index==kindeleidx(ikind)) then
                atmkind(iatm)=ikind
                exit
            end if
        end do
    end do
    kindmag=0
    netchg=sum(a%charge)-nint(nelec)
    multispin=nint(naelec-nbelec)+1
    if (ncenter>500) ioptmethod=2 !LBFGS is more suitable for large system than BFGS
    lastinpname=filename
end if

do while(.true.)
    !do iatm=1,ncenter
    !    ikind=atmkind(iatm)
    !    write(*,"(' Atom',i6,1x,a,'  Kind index:',i3,'  Kind name: ',a,'  Kind elem idx:',i3)") iatm,a(iatm)%name,ikind,kindname(ikind),kindeleidx(ikind)
    !end do
    write(*,*)
    write(*,*) "-11 Enter the interface for geometry operations"
    write(*,*) "-10 Return"
    write(*,*) "-9 Other settings"
    write(*,*) "-7 Set direction(s) of applying periodic boundary condition, current: "//PBCdir
    if (itask==3.or.itask==4.or.itask==6.or.itask==7) then !Optimizing minimum, MD, TS search
        if (itask==6) write(*,"(a,i6)") " -6 Set frequency of writing molecular dynamics trajectory, current:",nMDsavefreq
        if (iMDformat==1) write(*,*) "-5 Choose format of outputted trajectory, current: xyz"
        if (iMDformat==2) write(*,*) "-5 Choose format of outputted trajectory, current: dcd"
        if (iMDformat==3) write(*,*) "-5 Choose format of outputted trajectory, current: pdb"
    end if
    if (iatomcharge==0) write(*,*) "-4 Calculate atomic charges, current: None"
    if (iatomcharge==1) write(*,*) "-4 Calculate atomic charges, current: Mulliken"
    if (iatomcharge==2) write(*,*) "-4 Calculate atomic charges, current: Lowdin"
    if (iatomcharge==3) write(*,*) "-4 Calculate atomic charges, current: Hirshfeld"
    if (iatomcharge==4) write(*,*) "-4 Calculate atomic charges, current: Hirshfeld-I"
    if (iatomcharge==5) write(*,*) "-4 Calculate atomic charges, current: Voronoi"
    if (iatomcharge==6) write(*,*) "-4 Calculate atomic charges, current: RESP"
    if (iatomcharge==7) write(*,*) "-4 Calculate atomic charges, current: REPEAT"
    if (ioutcube==0) write(*,*) "-3 Set exporting cube file, current: None"
    if (ioutcube==1) write(*,*) "-3 Set exporting cube file, current: Electron density"
    if (ioutcube==2) write(*,*) "-3 Set exporting cube file, current: ELF"
    if (ioutcube==3) write(*,*) "-3 Set exporting cube file, current: XC potential"
    if (ioutcube==4) write(*,*) "-3 Set exporting cube file, current: Hartree potential (negative of ESP)"
    if (ioutcube==5) write(*,*) "-3 Set exporting cube file, current: Electric field"
    if (ioutcube==6) write(*,"(a,i6,a,i6)") " -3 Set exporting cube file, current: MOs, with NHOMO=",NHOMO,", NLUMO=",NLUMO
    if (imolden==0) write(*,*) "-2 Toggle exporting .molden file for Multiwfn, current: No"
    if (imolden==1) write(*,*) "-2 Toggle exporting .molden file for Multiwfn, current: Yes"
    if (itask==1) write(*,*) "-1 Choose task, current: Energy"
    if (itask==2) write(*,*) "-1 Choose task, current: Energy + force"
    if (itask==3) write(*,*) "-1 Choose task, current: Optimizing structure"
    if (itask==4) write(*,*) "-1 Choose task, current: Optimizing structure and cell"
    if (itask==5) write(*,*) "-1 Choose task, current: Vibrational analysis"
    if (itask==6) write(*,*) "-1 Choose task, current: Molecular dynamics"
    if (itask==7) write(*,*) "-1 Choose task, current: Searching transition state"
    if (itask==8) write(*,*) "-1 Choose task, current: Nudge-elastic band"
    if (itask==9) write(*,*) "-1 Choose task, current: NMR"
    if (itask==10) write(*,*) "-1 Choose task, current: Polarizability"
    if (itask==11) write(*,*) "-1 Choose task, current: BSSE"
    if (itask==12) write(*,*) "-1 Choose task, current: BAND"
    if (itask==13) write(*,*) "-1 Choose task, current: Real-time propagation for electron dynamics"
    if (itask==14) write(*,*) "-1 Choose task, current: Path-integral molecular dynamics (PIMD)"
    write(*,*) " 0 Generate input file now!"
    write(*,*) " 1 Choose theoretical method, current: "//trim(method)
    if (method/="GFN1-xTB".and.method/="PM6".and.method/="SCC-DFTB") write(*,*) " 2 Choose basis set and pseudopotential, current: "//trim(basname(ibas))
    if (index(method,"MP2")==0.and.index(method,"RPA")==0.and.index(method,"GW")==0.and.method/="GFN1-xTB".and.method/="PM6".and.method/="SCC-DFTB".and.method/="BEEFVDW") then
        if (idispcorr==0) write(*,*) " 3 Set dispersion correction, current: None"
        if (idispcorr==1) write(*,*) " 3 Set dispersion correction, current: DFT-D3"
        if (idispcorr==2) write(*,*) " 3 Set dispersion correction, current: DFT-D3(BJ)"
        if (idispcorr==5) write(*,*) " 3 Set dispersion correction, current: rVV10"
    end if
    if (iLSSCF==0) then
        if (idiagOT==1) then
            write(*,*) " 4 Switching between diagonalization and OT, current: Diagonalization"
            if (method/="PM6") then !I found PM6 can only use Direct mixing
                if (imixing==1) write(*,*) " 5 Set density matrix mixing, current: Direct mixing + DIIS"
                if (imixing==2) write(*,*) " 5 Set density matrix mixing, current: Broyden mixing"
                if (imixing==3) write(*,*) " 5 Set density matrix mixing, current: Pulay mixing"
            end if
            if (ismear==0) write(*,*) " 6 Toggle smearing electron occupation, current: No"
            if (ismear==1) write(*,*) " 6 Toggle smearing electron occupation, current: Yes"
        else if (idiagOT==2) then
            write(*,*) " 4 Switching between diagonalization and OT, current: OT"
            if (iouterSCF==0) write(*,*) " 5 Toggle using outer SCF process, current: No"
            if (iouterSCF==1) write(*,*) " 5 Toggle using outer SCF process, current: Yes"
        end if
    end if
    if (iSCCS==0) write(*,*) " 7 Toggle using self-consistent continuum solvation (SCCS), current: No"
    if (iSCCS==1) write(*,*) " 7 Toggle using self-consistent continuum solvation (SCCS), current: Yes"
    if (ikpoint1==1.and.ikpoint2==1.and.ikpoint3==1) then
        write(*,*) " 8 Set k-points, current: GAMMA only"
    else
        write(*,"(a,3i3)") "  8 Set k-points, current: MONKHORST-PACK",ikpoint1,ikpoint2,ikpoint3
    end if
    !!!! Below are task specific options
    if (itask>=3.and.itask<=7) then
        if (natmcons==0) then
            write(*,*) " 9 Set atom position constraint, current: None"
        else
            write(*,"(a,i6)") "  9 Set atom position constraint, current:",natmcons
        end if
    end if
    if (itask==6) then
        if (ithermostat==0) write(*,*) "10 Set thermostat, current: None"
        if (ithermostat==1) write(*,*) "10 Set thermostat, current: Adaptive-Langevin"
        if (ithermostat==2) write(*,*) "10 Set thermostat, current: Canonical sampling through velocity rescaling"
        if (ithermostat==3) write(*,*) "10 Set thermostat, current: Generalized Langevin Equation (GLE)"
        if (ithermostat==4) write(*,*) "10 Set thermostat, current: Nose-Hoover"
        if (ithermostat>0) then
            !if (nthermoatm==ncenter) then !This option is fully misleading!
            !    write(*,*) "11 Set region for the thermostat, current: All atoms"
            !else
            !    write(*,"(a,i6,' atoms')") " 11 Set region for the thermostat, number of current atoms:",nthermoatm
            !end if
        end if
        if (ibarostat==0) write(*,*) "12 Set barostat, current: None"
        if (ibarostat==1) write(*,*) "12 Set barostat, current: Yes, flexible cell"
        if (ibarostat==2) write(*,*) "12 Set barostat, current: Yes, isotropic cell"
        if (inoSCFinfo==0) write(*,*) "13 Toggle suppressing printing SCF information during MD, current: No"
        if (inoSCFinfo==1) write(*,*) "13 Toggle suppressing printing SCF information during MD, current: Yes"
    else if (itask==3.or.itask==4) then
        if (ioptmethod==1) write(*,*) "10 Set optimization method, current: BFGS"
        if (ioptmethod==2) write(*,*) "10 Set optimization method, current: LBFGS"
        if (ioptmethod==3) write(*,*) "10 Set optimization method, current: CG"
    else if (itask==5) then
        if (ioutvibmol==0) write(*,*) "10 Toggle exporting Molden file recording vibrational modes, current: No"
        if (ioutvibmol==1) write(*,*) "10 Toggle exporting Molden file recording vibrational modes, current: Yes"
    end if
    if (itask==4) then
        write(*,*) "11 Set constraint of cell length(s) during optimization, current: "//trim(cellfix)
        if (iprestype==1) write(*,"(a,1PE13.5,' bar')") " 12 Set external pressure, current:",Piso
        if (iprestype==2) write(*,"(a)") " 12 Set external pressure, current: Anisotropic"
    end if
    if (iTDDFT==0) then
        write(*,*) "15 Toggle calculating excited states via TDDFT, current: No"
    else if (iTDDFT==1) then
        write(*,*) "15 Toggle calculating excited states via TDDFT, current: Yes"
        write(*,"(' 16 Set number of excited states to solve by TDDFT, current:',i5)") nstates_TD
        if (.not.(multispin>1.or.any(kindmag(1:nkind)/=0))) then !Current is closed-shell
            if (iTDtriplet==0) write(*,*) "17 Toggle spin of the excited states to be calculated, current: Singlet"
            if (iTDtriplet==1) write(*,*) "17 Toggle spin of the excited states to be calculated, current: Triplet"
        end if
        if (isTDA==0) write(*,*) "18 Toggle using sTDA approximation in TDDFT, current: No"
        if (isTDA==1) write(*,*) "18 Toggle using sTDA approximation in TDDFT, current: Yes"
        !At least for CP2K 8.1, 9.1, the occupied NTOs are wrong, so do not let user know this feature
        !if (iNTO==0) write(*,*) "19 Toggle if performing NTO analysis, current: No"
        !if (iNTO==1) write(*,*) "19 Toggle if performing NTO analysis, current: Yes"
        if (itask==2.or.itask==3.or.itask==4.or.itask==5.or.itask==6.or.itask==7.or.itask==8) then !Need force
            write(*,"(a,i5)") " 20 Choose the state to evaluate force, current:",istate_force
        end if
    end if
    read(*,*) isel
    
    if (isel==-11) then
        natmold=ncenter
        call geom_operation
        !User may use geometry operation interface to reorder atoms or construct supercell, so we need to rebuild atom kind information
        if (any(a(1:size(a_org))%index/=a_org%index).or.ncenter/=natmold) then
            deallocate(atmkind)
            if (allocated(atmcons)) deallocate(atmcons)
            natmcons=0
            !if (allocated(thermoatm)) deallocate(thermoatm)
            !nthermoatm=0
            ithermostat=0
            goto 10 
        end if
    else if (isel==-10) then
        return
    else if (isel==-9) then
        do while(.true.)
            write(*,*)
            write(*,*) "                     ---------- Other settings ----------"
            write(*,*) "0 Return"
            write(*,"(a,i5)") " 1 Set net charge, current:",netchg
            write(*,"(a,i5)") " 2 Set spin multiplicity, current:",multispin
            write(*,"(a,3i3)") " 3 Set number of repetitions of the cell in X, Y, Z, current:",nrep1,nrep2,nrep3
            if (ifineXCgrid==0) write(*,"(a)") " 4 Toggle using finer grid for exchange-correlation part, current: No"
            if (ifineXCgrid==1) write(*,"(a)") " 4 Toggle using finer grid for exchange-correlation part, current: Yes"
            write(*,"(a,i5,' and',i4,' Ry')") " 5 Set CUTOFF and REL_CUTOFF, current:",CUTOFF,REL_CUTOFF
            if (imoment==0) write(*,*) "6 Print moments, current: No"
            if (imoment==1) write(*,*) "6 Print moments, current: Yes"
            if (PBCdir=="NONE") write(*,"(a,3f8.3,' A')") " 7 Set vacuum size in both sides of X,Y,Z, current:",vacsizex*b2a,vacsizey*b2a,vacsizez*b2a
            if (PBCdir=="X") write(*,"(a,2f8.3,' A')") " 7 Set vacuum size in both sides of Y and Z, current:",vacsizey*b2a,vacsizez*b2a
            if (PBCdir=="Y") write(*,"(a,2f8.3,' A')") " 7 Set vacuum size in both sides of X and Z, current:",vacsizex*b2a,vacsizez*b2a
            if (PBCdir=="Z") write(*,"(a,2f8.3,' A')") " 7 Set vacuum size in both sides of X and Y, current:",vacsizex*b2a,vacsizey*b2a
            if (PBCdir=="XY") write(*,"(a,f8.3,' A')") " 7 Set vacuum size in both sides of Z, current:",vacsizez*b2a
            if (PBCdir=="XZ") write(*,"(a,f8.3,' A')") " 7 Set vacuum size in both sides of Y, current:",vacsizey*b2a
            if (PBCdir=="YZ") write(*,"(a,f8.3,' A')") " 7 Set vacuum size in both sides of X, current:",vacsizex*b2a
            if (iDFTplusU==0) write(*,*) "8 Toggle using DFT+U, current: No"
            if (iDFTplusU==1) write(*,*) "8 Toggle using DFT+U, current: Yes"
            if (all(kindmag(1:nkind)==0)) then
                write(*,*) "9 Define atomic magnetization"
            else
                nmagset=count(kindmag(1:nkind)/=0)
                write(*,"(a,i3,a)") " 9 Redefine atomic magnetization, current: Manually defined",nmagset," kinds"
            end if
            if (iprintlevel==0) write(*,*) "10 Choose printing level of output information, current: Silent"
            if (iprintlevel==1) write(*,*) "10 Choose printing level of output information, current: Low"
            if (iprintlevel==2) write(*,*) "10 Choose printing level of output information, current: Medium"
            if (iprintlevel==3) write(*,*) "10 Choose printing level of output information, current: High"
            if (all(efieldvec==0)) then
                write(*,"(a)") " 11 Set external electric field vector"
            else
                write(*,"(a,3f8.5,' a.u.')") " 11 Set external electric field vector, current:",efieldvec
            end if
            if (nADDED_MOS==-1) then
                write(*,*) "12 Set number of virtual orbitals to solve, current: All"
            else
                write(*,"(a,i6)") " 12 Set number of virtual orbitals to solve, current:",nADDED_MOS
            end if
            if (icentering==0) write(*,*) "13 Toggle centering the coordinates of the system in the box, current: No"
            if (icentering==1) write(*,*) "13 Toggle centering the coordinates of the system in the box, current: Yes"
            if (ioutorbene==0) write(*,*) "14 Toggle printing orbital energies and occupancies after SCF, current: No"
            if (ioutorbene==1) write(*,*) "14 Toggle printing orbital energies and occupancies after SCF, current: Yes"
            !if (ioutSbas==0) write(*,*) "15 Toggle outputting overlap matrix to a file, current: No"
            !if (ioutSbas==1) write(*,*) "15 Toggle outputting overlap matrix to a file, current: Yes"
            !if (idiaglib==1) write(*,*) "20 Choose diagonalization library, current: Default"
            !if (idiaglib==2) write(*,*) "20 Choose diagonalization library, current: ELPA"
            !if (idiaglib==3) write(*,*) "20 Choose diagonalization library, current: Scalapack"
            if (iLRIGPW==0) write(*,*) "20 Toggle using LRIGPW instead of GPW to accelerate calculation, current: No"
            if (iLRIGPW==1) write(*,*) "20 Toggle using LRIGPW instead of GPW to accelerate calculation, current: Yes"
            if (iLSSCF==0) write(*,*) "21 Toggle using Linear Scaling Self Consistent Field Method, current: No"
            if (iLSSCF==1) write(*,*) "21 Toggle using Linear Scaling Self Consistent Field Method, current: Yes"
            if (ifPBC=="XYZ") then
                if (iPSOLVER==1) write(*,*) "22 Set Poisson solver, current: PERIODIC"
                if (iPSOLVER==2) write(*,*) "22 Set Poisson solver, current: ANALYTIC"
                if (iPSOLVER==3) write(*,*) "22 Set Poisson solver, current: MT"
                if (iPSOLVER==4) write(*,*) "22 Set Poisson solver, current: WAVELET"
            else
                if (iPSOLVER==1) write(*,*) "22 Set Poisson solver and automatically set vacuum size, current: PERIODIC"
                if (iPSOLVER==2) write(*,*) "22 Set Poisson solver and automatically set vacuum size, current: ANALYTIC"
                if (iPSOLVER==3) write(*,*) "22 Set Poisson solver and automatically set vacuum size, current: MT"
                if (iPSOLVER==4) write(*,*) "22 Set Poisson solver and automatically set vacuum size, current: WAVELET"
            end if
            if (PBCdir=="XYZ".and.iPSOLVER==1) then
                if (idipcorr==0) write(*,*) "23 Set surface dipole correction, current: None"
                if (idipcorr==1) write(*,*) "23 Set surface dipole correction, current: X direction"
                if (idipcorr==2) write(*,*) "23 Set surface dipole correction, current: Y direction"
                if (idipcorr==3) write(*,*) "23 Set surface dipole correction, current: Z direction"
            end if
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
                write(*,*) "Input net charge, e.g. 1"
                read(*,*) netchg
            else if (isel2==2) then
                write(*,*) "Input spin multiplicity, e.g. 3"
                read(*,*) multispin
            else if (isel2==3) then
                write(*,*) "Input number of repetitions of the cell in X, Y, Z, e.g. 2,1,2"
                read(*,*) nrep1,nrep2,nrep3
            else if (isel2==4) then
                if (ifineXCgrid==0) then
                    ifineXCgrid=1
                else
                    ifineXCgrid=0
                end if
            else if (isel2==5) then
                write(*,*) "Input CUTOFF and REL_CUTOFF in Ry, e.g. 350,50"
                read(*,*) CUTOFF,REL_CUTOFF
            else if (isel2==6) then
                if (imoment==1) then
                    imoment=0
                else
                    imoment=1
                end if
            else if (isel2==7) then
                write(*,*) "Note: The inputted size will be applied to both sides of non-periodic direction"
                if (PBCdir=="NONE") then
                    write(*,*) "Input vacuum size in X,Y,Z in Angstrom, e.g. 5.0,6.5,5.0"
                    if (iPSOLVER==4) write(*,"(a)") " Note: WAVELET Poisson solver requires cubic cell, &
                    Multiwfn will determine the longest cell size and set it to all directions"
                    read(*,*) vacsizex,vacsizey,vacsizez
                    vacsizex=vacsizex/b2a;vacsizey=vacsizey/b2a;vacsizez=vacsizez/b2a
                else if (PBCdir=="X") then
                    write(*,*) "Input vacuum size in Y and Z in Angstrom, e.g. 5.0,6.5"
                    read(*,*) vacsizey,vacsizez
                    vacsizey=vacsizey/b2a;vacsizez=vacsizez/b2a
                else if (PBCdir=="Y") then
                    write(*,*) "Input vacuum size in X and Z in Angstrom, e.g. 5.0,6.5"
                    read(*,*) vacsizex,vacsizez
                    vacsizex=vacsizex/b2a;vacsizez=vacsizez/b2a
                else if (PBCdir=="Z") then
                    write(*,*) "Input vacuum size in X and Y in Angstrom, e.g. 5.0,6.5"
                    read(*,*) vacsizex,vacsizey
                    vacsizex=vacsizex/b2a;vacsizey=vacsizey/b2a
                else if (PBCdir=="XY") then
                    write(*,*) "Input vacuum size in Z in Angstrom, e.g. 5.0"
                    read(*,*) vacsizez
                    vacsizez=vacsizez/b2a
                else if (PBCdir=="XZ") then
                    write(*,*) "Input vacuum size in Y in Angstrom, e.g. 5.0"
                    read(*,*) vacsizey
                    vacsizey=vacsizey/b2a
                else if (PBCdir=="YZ") then
                    write(*,*) "Input vacuum size in X in Angstrom, e.g. 5.0"
                    read(*,*) vacsizex
                    vacsizex=vacsizex/b2a
                end if
            else if (isel2==8) then
                if (iDFTplusU==1) then
                    iDFTplusU=0
                else
                    iDFTplusU=1
                    write(*,"(a)") " IMPORTANT NOTE: DO NOT forget to manually replace the default DFT+U parameters in the generated input file with proper value!"
                end if
            else if (isel2==9) then
                !do ikind=1,nkind
                !    write(*,"(' #',i3,':  Kind name: ',a5,' Element: ',a,'  Magnetization:',i3,'   Natoms:',i5)") &
                !    ikind,kindname(ikind),ind2name(kindeleidx(ikind)),kindmag(ikind),count(atmkind(:)==ikind)
                !end do
                do while(.true.)
                    write(*,*)
                    write(*,*) "Current magnetization status:"
                    idx=0
                    do ikind=1,nkind
                        ncount=count(atmkind(:)==ikind)
                        if (ncount>0) then
                            idx=idx+1
                            write(*,"(' #',i3,':  Kind name: ',a5,' Element: ',a,'  Magnetization:',i3,'   Natoms:',i5)") &
                            idx,kindname(ikind),ind2name(kindeleidx(ikind)),kindmag(ikind),ncount
                        end if
                    end do
                    write(*,*)
                    write(*,*) "Input indices of the atoms to define magnetization, e.g. 1,5-10,13,19"
                    write(*,*) "To change kind name, input old and new name, e.g. Fe_1 Fe_B"
                    write(*,*) "To exit, inputting ""q"""
                    read(*,"(a)") c2000tmp
                    if (c2000tmp=="q") then
                        multispin=sum(kindmag(atmkind(:)))+1
                        exit
                    else if (iachar(c2000tmp(1:1))<48.or.iachar(c2000tmp(1:1))>57) then !The first character is not a number
                        read(c2000tmp,*) c80tmp,c80tmp2
                        do ikind=1,nkind
                            if (trim(kindname(ikind))==trim(c80tmp)) kindname(ikind)=trim(c80tmp2)
                        end do
                        cycle
                    end if
                    call str2arr(c2000tmp,ntmp,tmparr)
                    if (all(a(tmparr(1:ntmp))%index==a(tmparr(1))%index)) then
                        nkind=nkind+1
                        write(*,*) "Input magnetization (difference of alpha and beta electrons), e.g. 4"
                        read(*,*) kindmag(nkind)
                        atmkind(tmparr(1:ntmp))=nkind
                        iele=a(tmparr(1))%index
                        kindeleidx(nkind)=iele
                        isuffix=0
                        do while(.true.) !Test which suffix can be used
                            if (isuffix>0) then
                                write(c80tmp,"(i5)") isuffix
                                c80tmp=trim(ind2name(iele))//'_'//trim(adjustl(c80tmp))
                            else
                                c80tmp=ind2name(iele)
                            end if
                            if (all(kindname(1:nkind-1)/=trim(c80tmp))) then !New kind name
                                exit
                            else
                                ncount=count(kindname(atmkind(:))==trim(c80tmp))
                                if (ncount==0) exit
                            end if
                            isuffix=isuffix+1
                        end do
                        kindname(nkind)=trim(c80tmp)
                        write(*,*) "Done!"
                    else
                        write(*,*) "Error: Not all atoms you selected belong to the same element!"
                    end if
                end do
            else if (isel2==10) then
                write(*,*) "Choose printing level"
                write(*,*) "0 Silent"
                write(*,*) "1 Low"
                write(*,*) "2 Medium"
                write(*,*) "3 High"
                read(*,*) iprintlevel
            else if (isel2==11) then
                write(*,*) "Input external electric field vector in a.u., e.g. 0.0,0.0,0.025"
                read(*,*) efieldvec
                if (ifPBC>0.and.idiagOT==1) then
                    write(*,"(a)") " NOTE: Because this is a periodic system, &PERIODIC_EFIELD will be used in the generated input file to specify the field, &
                    in this case you must use OT rather than diagonalizaton!"
                end if
            else if (isel2==12) then
                write(*,*) "Input number of virtual orbitals to solve, e.g. 30"
                write(*,*) "If inputting -1, then all virtual orbitals will be solved" !Work since CP2K 9.1
                read(*,*) nADDED_MOS
            else if (isel2==13) then
                if (icentering==0) then
                    icentering=1
                else if (icentering==1) then
                    icentering=0
                end if
            else if (isel2==14) then
                if (ioutorbene==0) then
                    ioutorbene=1
                else
                    ioutorbene=0
                end if
            else if (isel2==15) then
                if (ioutSbas==0) then
                    ioutSbas=1
                else
                    ioutSbas=0
                end if
            else if (isel2==20) then
                !write(*,*) "Choose diagonalization library"
                !write(*,*) "1 Default"
                !write(*,*) "2 ELPA"
                !write(*,*) "3 Scalapack"
                !read(*,*) idiaglib
                if (iLRIGPW==0) then
                    iLRIGPW=1
                else
                    iLRIGPW=0
                end if
            else if (isel2==21) then
                if (iLSSCF==0) then
                    iLSSCF=1
                    if (idiagOT==2) then
                        idiagOT=1
                        write(*,"(a)") " Note: OT is disabled because it is incompatible with linear scaling self consistent field method"
                    end if
                else
                    iLSSCF=0
                end if
            else if (isel2==22) then
                write(*,*) "Choose Poisson solver:"
                write(*,*) "1 PERIODIC"
                write(*,*) "2 ANALYTIC"
                if (PBCdir=="NONE".or.PBCdir=="XY".or.PBCdir=="XZ".or.PBCdir=="YZ") write(*,*) "3 MT"
                if (PBCdir=="NONE".or.PBCdir=="XZ".or.PBCdir=="XYZ") write(*,*) "4 WAVELET"
                read(*,*) iPSOLVER
                call determine_vacuumsize(iPSOLVER,vacsizex,vacsizey,vacsizez,icentering) !Automatically set proper vacuum sizes
            else if (isel2==23) then
                write(*,*) "0 Do not use surface dipole correction"
                write(*,*) "1 Use surface dipole correction in X direction"
                write(*,*) "2 Use surface dipole correction in Y direction"
                write(*,*) "3 Use surface dipole correction in Z direction"
                read(*,*) idipcorr
            end if
        end do
    else if (isel==-7) then
        write(*,*) "Input one of following strings to specify periodic boundary condition (PBC)"
        write(*,*) "NONE, X, XY, XYZ, XZ, Y, YZ, Z"
        read(*,*) PBCdir
        if (PBCdir/="NONE".and.PBCdir/="X".and.PBCdir/="Y".and.PBCdir/="Z".and.PBCdir/="XY".and.PBCdir/="XY".and.PBCdir/="YZ") then
            write(*,*) "Error: Your input cannot be recognized!"
            PBCdir="XYZ"
            cycle
        end if
        !Automatically set proper Poisson solver and vacuum sizes
        if (PBCdir=="NONE") then !Use WAVELET for 0D, usually best choice
            iPSOLVER=4
            write(*,*) "Note: Poisson solver has been automatically changed to WAVELET"
        else if (PBCdir=="X".or.PBCdir=="Y".or.PBCdir=="Z".or.PBCdir=="XYZ") then !Use PERIODIC for 1D and 3D
            iPSOLVER=1
            write(*,*) "Note: Poisson solver has been automatically changed to PERIODIC"
        else if (PBCdir=="XY".or.PBCdir=="XZ".or.PBCdir=="YZ") then !Use MT for 2D. Needs vaccum size in each side >= half of system
            iPSOLVER=3
            write(*,*) "Note: Poisson solver has been automatically changed to MT"
        end if
        call determine_vacuumsize(iPSOLVER,vacsizex,vacsizey,vacsizez,icentering)
    else if (isel==-6) then
        write(*,*) "Input frequency of writing molecular dynamics trajectory, 1 means every step"
        read(*,*) nMDsavefreq
    else if (isel==-5) then
        if (itask==4) write(*,*) "Choose the format for recording trajectory of optimization"
        if (itask==6.or.itask==14) write(*,*) "Choose the format for recording trajectory of molecular dynamics"
        write(*,*) "1 xyz (Simplest. Does not contain cell information)"
        write(*,*) "2 dcd (Binary file, smallest size. Containing cell information)"
        write(*,*) "3 pdb (Containing cell information, but accuracy of coordinates is limited)"
        read(*,*) iMDformat
    else if (isel==-4) then
        write(*,*) "Printing which kind of atomic charge?"
        write(*,*) "0 None"
        write(*,*) "1 Mulliken"
        write(*,*) "2 Lowdin"
        write(*,*) "3 Hirshfeld"
        write(*,*) "4 Hirshfeld-I"
        write(*,*) "5 Voronoi"
        write(*,*) "6 RESP"
        write(*,*) "7 REPEAT"
        read(*,*) iatomcharge
        if (ifPBC==0.and.iatomcharge==7) then
            write(*,*) "Error: REPEAT can only be used for periodic system"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            iatomcharge=0
        end if
    else if (isel==-3) then
        write(*,*) "Output cube file for which real space function?"
        write(*,"(a)") " -1 Just for printing HOMO and LUMO energies as well as HOMO-LUMO gap (i.e. Outputting HOMO and LUMO cubes only)"
        write(*,*) "0 None"
        write(*,*) "1 Electron density (also with spin density for unrestricted calculation)"
        write(*,*) "2 Electron localization function (ELF)"
        write(*,*) "3 Exchange-correlation potential"
        write(*,*) "4 Hartree potential (negative of ESP)"
        write(*,*) "5 Each component of electric field"
        write(*,*) "6 Molecular orbital(s)"
        read(*,*) ioutcube
        if (ioutcube==6) then
            write(*,*) "Output how many highest occupied orbitals? e.g. 5"
            write(*,*) "If inputting -1, all occupied orbitals will be outputted"
            read(*,*) NHOMO
            write(*,*) "Output how many lowest unoccupied orbitals? e.g. 5"
            write(*,*) "If inputting -1, all unoccupied orbitals will be outputted"
            read(*,*) NLUMO
        else if (ioutcube==-1) then
            ioutcube=6
            NHOMO=1
            NLUMO=1
        end if
    else if (isel==-2) then
        if (imolden==0) then
            imolden=1
        else
            imolden=0
        end if
    else if (isel==-1) then
        write(*,*) "Please select task"
        write(*,*) "1 Energy"
        write(*,*) "2 Energy + force"
        write(*,*) "3 Optimizing structure (cell is fixed)"
        write(*,*) "4 Optimizing both structure and cell"
        write(*,*) "5 Vibrational analysis"
        write(*,*) "6 Molecular dynamics (MD)"
        write(*,*) "7 Searching transition state using dimer algorithm"
        !write(*,*) "8 Nudge-elastic band (NEB)"
        write(*,*) "9 NMR"
        write(*,*) "10 Polarizability"
        write(*,*) "11 Correct for basis set superposition error (BSSE)"
        !write(*,*) "12 BAND"
        write(*,*) "13 Real-time propagation for electron dynamics"
        write(*,*) "14 Path-integral molecular dynamics (PIMD)"
        read(*,*) itask
        if (itask==9.and.ibas<=5) then
            ibas=10 !Use 6-31G* if current basis set is pseudopotential basis set
        else if (itask==4) then !Use pdb to record variable cell during cell optimization
            iMDformat=3
        else if (itask==5) then !Vibrational analysis
            iprintlevel=2 !Use medium printing level
        else if (itask==11) then !BSSE
            write(*,*) "Input atoms in fragment 1, e.g. 1,3-6,10,14"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,nfrag1)
            if (allocated(frag1)) deallocate(frag1)
            allocate(frag1(nfrag1))
            call str2arr(c2000tmp,nfrag1,frag1)
            write(*,*) "Input charge for fragment 1, e.g. 0"
            read(*,*) frag1chg
            write(*,*) "Input spin multiplicity for fragment 1, e.g. 1"
            read(*,*) frag1multi
            write(*,*) "Input atoms in fragment 2, e.g. 2,7-9,11-13"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,nfrag2)
            if (allocated(frag2)) deallocate(frag2)
            allocate(frag2(nfrag2))
            call str2arr(c2000tmp,nfrag2,frag2)
            write(*,*) "Input charge for fragment 2, e.g. 0"
            read(*,*) frag2chg
            write(*,*) "Input spin multiplicity for fragment 2, e.g. 1"
            read(*,*) frag2multi
            write(*,*) "Input spin multiplicity for whole system, e.g. 1"
            read(*,*) totalmulti
        else if (itask==13) then !Real-time propagation of electron
            imoment=1
            !iatomcharge=1 !If enable this, Mulliken charge will print every step
        end if
        if (itask==1) then !Single point, medicore accuracy
            CUTOFF=350
            REL_CUTOFF=50
        else if (itask==6.or.itask==14) then !MD and PIMD, do not need high accuracy (assume without barostat)
            CUTOFF=300
            REL_CUTOFF=40
        else !If task involves energy derivative, or TDDFT, use higher cutoff
            CUTOFF=400
            REL_CUTOFF=55
        end if
    else if (isel==1) then !Functionals description: https://manual.cp2k.org/trunk/CP2K_INPUT/ATOM/METHOD/XC/XC_FUNCTIONAL.html
        !write(*,*) "-1 Molecular mechanism (MM)"
        write(*,*) "1 Pade (LDA)"
        write(*,*) "2 PBE        -2 revPBE     -3 PBEsol"
        write(*,*) "3 TPSS        4 BP86        5 BLYP"
        write(*,*) "6 PBE0       -6 PBE0 with ADMM"
        write(*,*) "7 B3LYP      -7 B3LYP with ADMM"
        write(*,*) "8 HSE06      -8 HSE06 with ADMM"
        write(*,*) "9 BEEF-vdW"
        write(*,*) "11 B97M-rV (via LibXC)       12 MN15L (via LibXC)"
        write(*,*) "13 SCAN (via LibXC)          14 r2SCAN (via LibXC)"
        write(*,*) "15 RPBE (via LibXC)          16 revTPSS (via LibXC)"
        write(*,*) "20 RI-MP2                    21 RI-SCS-MP2"
        write(*,*) "22 RI-(EXX+RPA)@PBE          23 GW@PBE"
        write(*,*) "25 RI-B2PLYP  26 RI-B2GP-PLYP  27 RI-DSD-BLYP  28 RI-revDSD-PBEP86 with ADMM"
        write(*,*) "30 GFN1-xTB      40 PM6      50 SCC-DFTB"
        read(*,*) isel2
        iwfc=0 !Do not involve wavefunction-based correlation
        iHFX=0 !Do not involve HF exchange
        if (isel2==1) method="Pade"
        if (isel2==2) method="PBE"
        if (isel2==-2) method="revPBE"
        if (isel2==-3) method="PBEsol"
        if (isel2==3) method="TPSS"
        if (isel2==4) method="BP"
        if (isel2==5) method="BLYP"
        if (isel2==6) method="PBE0"
        if (isel2==-6) method="PBE0_ADMM"
        if (isel2==7) method="B3LYP"
        if (isel2==-7) method="B3LYP_ADMM"
        if (isel2==8) method="HSE06"
        if (isel2==-8) method="HSE06_ADMM"
        if (isel2==9) then
            method="BEEFVDW"
            idispcorr=0
        end if
        if (isel2==11) then
            method="B97M-rV_LIBXC"
            idispcorr=5
        end if
        if (isel2==12) method="MN15L_LIBXC"
        if (isel2==13) method="SCAN_LIBXC"
        if (isel2==14) method="r2SCAN_LIBXC"
        if (isel2==15) method="RPBE_LIBXC"
        if (isel2==16) method="revTPSS_LIBXC"
        if (isel2>=20.and.isel2<30) then !Involve wavefunction-based correlation
            if (isel2==20) method="RI-MP2"
            if (isel2==21) method="RI-SCS-MP2"
            if (isel2==22) method="RI-(EXX+RPA)@PBE"
            if (isel2==23) method="GW@PBE"
            if (isel2==25) method="RI-B2PLYP"
            if (isel2==26) method="RI-B2GP-PLYP"
            if (isel2==27) method="RI-DSD-BLYP"
            if (isel2==28) method="RI-revDSD-PBEP86_ADMM"
            iwfc=1
            if (ibas/=20) ibas=21 !Default to cc-TZ with RI_TZ
        end if
        if (isel2==30) method="GFN1-xTB"
        if (isel2==40) method="PM6"
        if (isel2==50) method="SCC-DFTB"
        if (isel2==-6.or.isel2==-7.or.isel2==-8.or.isel2==30.or.isel2==40.or.isel2==50) idiagOT=2 !When ADMM is used, OT must also be used. OT is suggested for GFN-xTB and PM6 dealing with large system
        if (isel2==40) imixing=1
        if (index(method,"SCAN")/=0) then
            write(*,"(a)") " NOTE: If you are using CP2K >=9.1, in the generated CP2K input file, it is suggested to replace &
            ""POTENTIAL_FILE_NAME  POTENTIAL"" with ""POTENTIAL_FILE_NAME  POTENTIAL_UZH"", and replace ""BASIS_SET_FILE_NAME  BASIS_MOLOPT"" with ""BASIS_SET_FILE_NAME  BASIS_MOLOPT_UZH"", &
            and manually specify proper GTH potential and corresponding valence basis set optimized for SCAN calculation"
        else if (index(method,"PBE0")/=0) then
            write(*,"(a)") " NOTE: If you are using CP2K >=9.1, in the generated CP2K input file, it is suggested to replace &
            ""POTENTIAL_FILE_NAME  POTENTIAL"" with ""POTENTIAL_FILE_NAME  POTENTIAL_UZH"", and replace ""BASIS_SET_FILE_NAME  BASIS_MOLOPT"" with ""BASIS_SET_FILE_NAME  BASIS_MOLOPT_UZH"", &
            and manually specify proper GTH potential and corresponding valence basis set optimized for PBE0 calculation"
        end if
        if (abs(isel2)==6.or.abs(isel2)==7.or.abs(isel2)==8.or.iwfc==1) then !Hybrid functionals and wavefunction-based correlation methods need HF exchange
            iHFX=1
            if (method=="RI-(EXX+RPA)@PBE".or.method=="GW@PBE") iHFX=0 !Based on PBE orbitals, HFX is not involved in SCF process
        end if
        
    else if (isel==2) then !Select basis set
        write(*,"(a)") " Note: <=5, 20, 21 correspond to GPW calculation using GTH pseudopotential, the other ones correspond to full electron GAPW calculation"
        do i=-10,30
            if (basname(i)/=" ") write(*,"(1x,i2,1x,a)") i,trim(basname(i))
        end do
        read(*,*) ibassel
        if (iwfc==1.and.ibassel/=20.and.ibassel/=21) then
            write(*,"(a)") " Error: To perform RI calculation, you must choose 20 or 21, because they are accompanied by auxiliary basis set" 
            write(*,*) "Press ENTER button to continue"
            read(*,*)
        else
            ibas=ibassel
            if (ibas>=10.and.ibas<=14) then
                iGAPW=1
            else
                iGAPW=0
            end if
        end if
    else if (isel==3) then
        write(*,*) "Choose dispersion correction method"
        write(*,*) "0 None"
        write(*,*) "1 DFT-D3"
        write(*,*) "2 DFT-D3(BJ)"
        write(*,*) "5 rVV10"
        read(*,*) idispcorr
    else if (isel==4) then
        if (idiagOT==1) then
            if (ikpoint1/=1.or.ikpoint2/=1.or.ikpoint3/=1) then
                write(*,*) "Error: OT can only be used for Gamma point!"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                cycle
            end if
            if (ismear==1) then
                write(*,*) "Error: OT cannot be used in combination with smearing!"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                cycle
            end if
            idiagOT=2
        else if (idiagOT==2) then
            idiagOT=1
        end if
    else if (isel==5) then
        if (idiagOT==1) then
            write(*,*) "Choose how to mixing old and new density matrices"
            write(*,*) "1 Direct mixing with DIIS (default, usually poor)"
            write(*,*) "2 Broyden mixing"
            write(*,*) "3 Pulay mixing"
            read(*,*) imixing
        else if (idiagOT==2) then
            if (iouterSCF==0) then
                iouterSCF=1
            else
                iouterSCF=0
            end if
        end if
    else if (isel==6) then
        if (ismear==0) then
            ismear=1
            nADDED_MOS=30
            write(*,"(a)") " Note: The number of virtual orbitals to solve has been changed to 30, please properly adjust if needed"
        else
            ismear=0
            nADDED_MOS=0
            write(*,"(a)") " Note: The number of virtual orbitals to solve has been changed to 0"
        end if
    else if (isel==7) then
        if (iSCCS==0) then
            iSCCS=1
        else
            iSCCS=0
        end if
    else if (isel==8) then
        write(*,*) "Input number of k-points of MONKHORST-PACK in three directions, e.g. 8,6,2"
        read(*,*) ikpoint1,ikpoint2,ikpoint3
        if (idiagOT==2.and.(ikpoint1/=1.or.ikpoint2/=1.or.ikpoint3/=1)) then
            write(*,"(a)") " Warning: OT can be used for Gamma point only! Now diagonalization is used instead"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            idiagOT=1
        end if
    else if (isel==9) then
        do while(.true.)
            write(*,*) "Input indices of the atoms to be constraint (fixed), e.g. 1,5,9-12,14-18"
            write(*,"(a)") " If inputting ""optH"", then only hydrogens will be optimized while others will be fixed"
            read(*,"(a)") c2000tmp
            if (.not.allocated(atmcons)) allocate(atmcons(ncenter))
            if (index(c2000tmp,"optH")/=0) then
                natmcons=0
                do iatm=1,ncenter
                    if (a(iatm)%index==1) cycle
                    natmcons=natmcons+1
                    atmcons(natmcons)=iatm
                end do
            else
                call str2arr(c2000tmp,natmcons,atmcons)
                if (natmcons>ncenter) then
                    write(*,*) "Error: The indices you inputted is invalid!"
                    cycle
                end if
            end if
            write(*,"(i8,' atoms will be fixed')") natmcons
            write(*,"(a)") " Note that the direction(s) of fixing can be manually set by changing ""COMPONENTS_TO_FIX"" &
            in the generated input file, by default they are fixed in all directions"
            exit
        end do
    else if (isel==10) then
        if (itask==6) then
            write(*,*) "0 Do not use thermostat"
            write(*,*) "1 Adaptive-Langevin thermostat"
            write(*,"(a)") " 2 Canonical sampling through velocity rescaling (CSVR, also known as V-rescale, recommended!)"
            write(*,*) "3 Generalized Langevin Equation (GLE) thermostat"
            write(*,*) "4 Nose-Hoover thermostat"
            read(*,*) ithermostat
            !if (ithermostat>0) then
            !    if (.not.allocated(thermoatm)) allocate(thermoatm(ncenter))
            !    nthermoatm=ncenter
            !    forall(iatm=1:ncenter) thermoatm(iatm)=iatm
            !end if
        else if (itask==3.or.itask==4) then
            write(*,*) "Choose optimization method"
            write(*,*) "1 BFGS (Best choice for most situations)"
            write(*,*) "2 LBFGS (Suitable for very large systems)"
            write(*,"(a)") " 3 Conjugate gradient (More robust than BFGS and LBFGS especially when initial geometry &
            is far from minimum, unfortunately more expensive. Try it for difficult cases)"
            read(*,*) ioptmethod
        else if (itask==5) then
            if (ioutvibmol==0) then
                ioutvibmol=1
            else
                ioutvibmol=0
            end if
        end if
    else if (isel==11) then
        if (itask==6) then
            !write(*,*) "Input indices of the atoms to whom the thermostat will be applied"
            !write(*,*) "For example: 1,5,9-12,14-18"
            !read(*,"(a)") c2000tmp
            !if (.not.allocated(thermoatm)) allocate(thermoatm(ncenter))
            !call str2arr(c2000tmp,nthermoatm,thermoatm)
        else if (itask==4) then
            write(*,*) "Input constraint on cell optimization, please input one of following terms:"
            write(*,*) "NONE, X, Y, Z, XY, XZ, YZ"
            read(*,"(a)") cellfix
        end if
    else if (isel==12) then
        if (itask==6) then
            write(*,*) "0 Do not use barostat"
            write(*,*) "1 Use barostat, flexible cell"
            write(*,*) "2 Use barostat, isotropic cell"
            read(*,*) ibarostat
            if (ibarostat>0) iMDformat=2 !Use dcd format to record variable cell size during MD
        else if (itask==4) then
            write(*,*) "0 Return"
            write(*,*) "1 Set isotropic external pressure"
            write(*,*) "2 Set anisotropic external pressure by specifying 9 pressure tensor components"
            read(*,*) isel2
            if (isel2==0) then
                cycle
            else if (isel2==1) then
                write(*,*) "Input external pressure in bar, e.g. 150"
                read(*,*) Piso
                iprestype=1
            else if (isel2==2) then
                do while(.true.)
                    write(*,*)
                    write(*,*) "Current pressure tensor:"
                    call showmatgau(Ptens,fileid=6,form="1PE14.5")
                    write(*,*) "Input for example 1,3,250.4 to define component XZ as 250.4"
                    write(*,*) "Input ""q"" can return"
                    read(*,"(a)") c80tmp
                    if (index(c80tmp,'q')/=0) then
                        iprestype=2
                        exit
                    else
                        read(c80tmp,*,iostat=ierror) idx,jdx,presval
                        if (ierror/=0) then
                            write(*,*) "Error: Unable to recognize your input!"
                            cycle
                        end if
                        Ptens(idx,jdx)=presval
                    end if
                end do
            end if
        end if
    else if (isel==13) then
        if (itask==6) then
            if (inoSCFinfo==0) then
                inoSCFinfo=1
            else
                inoSCFinfo=0
            end if
        end if
    else if (isel==15) then
        if (iTDDFT==0) then
            iTDDFT=1
            write(*,*) "Note: The TDDFT realized by CP2K employs Tamm-Dancoff approximation"
            write(*,"(/,a)") " If outputting .molden file containing all occupied and a batch of virtual orbitals for post-processing analysis? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y') then
                write(*,"(a)") " How many virtual orbitals to solve and record in the .molden file? e.g. 40"
                write(*,*) "You can input -1 or a very large number to solve all virtual orbitals"
                read(*,*) nADDED_MOS
                imolden=1
                idiagOT=1
            end if
        else
            iTDDFT=0
        end if
    !Below are specific for TDDFT
    else if (isel==16) then
        write(*,*) "Input number of excited states to solve, e.g. 5"
        read(*,*) nstates_TD
    else if (isel==17) then
        if (iTDtriplet==0) then
            iTDtriplet=1
        else
            iTDtriplet=0
        end if
    else if (isel==18) then
        if (isTDA==0) then
            isTDA=1
        else
            isTDA=0
        end if
    else if (isel==19) then
        if (iNTO==0) then
            iNTO=1
        else
            iNTO=0
        end if
    else if (isel==20) then
        if (itask==2.or.itask==3.or.itask==4.or.itask==5.or.itask==6.or.itask==7.or.itask==8) then !Need force
            write(*,*) "Input the index of the excited state for which force will be evaluated, e.g. 2"
            read(*,*) istate_force
        end if
    
    else if (isel==0) then
        exit
    end if
end do

open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "#Generated by Multiwfn"
call path2filename(filename,c200tmp)
write(ifileid,"(a)") "&GLOBAL"
call path2filename(outname,c200tmp)
write(ifileid,"(a)") "  PROJECT "//trim(c200tmp)
if (idiaglib==2) write(ifileid,"(a)") "  PREFERRED_DIAG_LIBRARY ELPA #Library for diagonalization"
if (idiaglib==3) write(ifileid,"(a)") "  PREFERRED_DIAG_LIBRARY SL #Library for diagonalization"
if (iprintlevel==0) write(ifileid,"(a)") "  PRINT_LEVEL SILENT"
if (iprintlevel==1) write(ifileid,"(a)") "  PRINT_LEVEL LOW"
if (iprintlevel==2) write(ifileid,"(a)") "  PRINT_LEVEL MEDIUM"
if (iprintlevel==3) write(ifileid,"(a)") "  PRINT_LEVEL HIGH"
if (itask==1) write(ifileid,"(a)") "  RUN_TYPE ENERGY"
if (itask==2) write(ifileid,"(a)") "  RUN_TYPE ENERGY_FORCE"
if (itask==3) write(ifileid,"(a)") "  RUN_TYPE GEO_OPT"
if (itask==4) write(ifileid,"(a)") "  RUN_TYPE CELL_OPT"
if (itask==5) write(ifileid,"(a)") "  RUN_TYPE VIBRATIONAL_ANALYSIS"
if (itask==6) write(ifileid,"(a)") "  RUN_TYPE MD"
if (itask==7) write(ifileid,"(a)") "  RUN_TYPE GEO_OPT"
if (itask==9.or.itask==10) write(ifileid,"(a)") "  RUN_TYPE LR"
if (itask==11) write(ifileid,"(a)") "  RUN_TYPE BSSE"
if (itask==12) write(ifileid,"(a)") "  RUN_TYPE BAND"
if (itask==13) write(ifileid,"(a)") "  RUN_TYPE RT_PROPAGATION"
if (itask==14) write(ifileid,"(a)") "  RUN_TYPE PINT"
write(ifileid,"(a)") "&END GLOBAL"
write(ifileid,"(/,a)") "&FORCE_EVAL"
write(ifileid,"(a)") "  METHOD Quickstep"

write(ifileid,"(a)") "  &SUBSYS"
if (nrep1/=1.or.nrep2/=1.or.nrep3/=1) then !This is needed even not for classical force field calculation
    write(ifileid,"(a)") "    &TOPOLOGY"
    write(ifileid,"(a,3i3)") "      MULTIPLE_UNIT_CELL",nrep1,nrep2,nrep3
    write(ifileid,"(a)") "    &END TOPOLOGY"
end if

!---- &CELL
write(ifileid,"(a)") "    &CELL"
if (ifPBC==3.and.PBCdir=="XYZ") then !Real 3D system and request calculation as 3D periodicity, use original cell vectors
    write(ifileid,"(a,3f15.8)") "      A",cellv1(:)*b2a
    write(ifileid,"(a,3f15.8)") "      B",cellv2(:)*b2a
    write(ifileid,"(a,3f15.8)") "      C",cellv3(:)*b2a
    write(ifileid,"(a)") "      PERIODIC "//trim(PBCdir)//" #Direction of applied PBC (geometry aspect)"
else !Low-dimensional case, allow to set vacuum size. Periodic directions employ original cell vectors, while other direction(s) employ actual size with vacuum size, and they are parallel to the Cartesian axes
    !Get extended size (system+vacuum) in X,Y,Z w.r.t boundary atoms
    xdist=(maxval(a%x)-minval(a%x)+2*vacsizex)*b2a
    ydist=(maxval(a%y)-minval(a%y)+2*vacsizey)*b2a
    zdist=(maxval(a%z)-minval(a%z)+2*vacsizez)*b2a
    if (PBCdir=="NONE") then
        if (iPSOLVER==4) then !WAVELET, needs cubic box
            tmp=max(max(xdist,ydist),zdist)
            write(ifileid,"(a,3f10.3)") "      ABC",tmp,tmp,tmp
        else
            write(ifileid,"(a,3f10.3)") "      ABC",xdist,ydist,zdist
        end if
    else
        if (PBCdir=="X") then
            write(ifileid,"(a,3f15.8)") "      A",cellv1(:)*b2a
            write(ifileid,"(a,3f15.8)") "      B",0D0,ydist,0D0
            write(ifileid,"(a,3f15.8)") "      C",0D0,0D0,zdist
        else if (PBCdir=="Y") then
            write(ifileid,"(a,3f15.8)") "      A",xdist,0D0,0D0
            write(ifileid,"(a,3f15.8)") "      B",cellv2(:)*b2a
            write(ifileid,"(a,3f15.8)") "      C",0D0,0D0,zdist
        else if (PBCdir=="Z") then
            write(ifileid,"(a,3f15.8)") "      A",xdist,0D0,0D0
            write(ifileid,"(a,3f15.8)") "      B",0D0,ydist,0D0
            write(ifileid,"(a,3f15.8)") "      C",cellv3(:)*b2a
        else if (PBCdir=="XY") then
            write(ifileid,"(a,3f15.8)") "      A",cellv1(:)*b2a
            write(ifileid,"(a,3f15.8)") "      B",cellv2(:)*b2a
            write(ifileid,"(a,3f15.8)") "      C",0D0,0D0,zdist
        else if (PBCdir=="XZ") then
            write(ifileid,"(a,3f15.8)") "      A",cellv1(:)*b2a
            write(ifileid,"(a,3f15.8)") "      B",0D0,ydist,0D0
            write(ifileid,"(a,3f15.8)") "      C",cellv3(:)*b2a
        else if (PBCdir=="YZ") then
            write(ifileid,"(a,3f15.8)") "      A",xdist,0D0,0D0
            write(ifileid,"(a,3f15.8)") "      B",cellv2(:)*b2a
            write(ifileid,"(a,3f15.8)") "      C",cellv3(:)*b2a
        else if (PBCdir=="XYZ") then
            write(ifileid,"(a,3f15.8)") "      A",xdist,0D0,0D0
            write(ifileid,"(a,3f15.8)") "      B",0D0,ydist,0D0
            write(ifileid,"(a,3f15.8)") "      C",0D0,0D0,zdist
        end if
    end if
    if (iPSOLVER==1) then
        write(ifileid,"(a)") "      PERIODIC XYZ #Direction(s) of applied PBC (geometry aspect)"
        write(*,"(a)") " Note: PERIODIC in the generated input file is changed to XYZ since current PSOLVER is PERIODIC"
    else
        write(ifileid,"(a)") "      PERIODIC "//trim(PBCdir)//" #Direction(s) of applied PBC (geometry aspect)"
    end if
end if
if (nrep1/=1.or.nrep2/=1.or.nrep3/=1) write(ifileid,"(a,3i3)") "      MULTIPLE_UNIT_CELL",nrep1,nrep2,nrep3
write(ifileid,"(a)") "    &END CELL"
if (icentering==1) then
    write(ifileid,"(a)") "    &TOPOLOGY"
    write(ifileid,"(a)") "      &CENTER_COORDINATES #Centering the atoms in the box"
    write(ifileid,"(a)") "      &END CENTER_COORDINATES"
    write(ifileid,"(a)") "    &END TOPOLOGY"
end if

!---- &COORD
write(ifileid,"(a)") "    &COORD"
do iatm=1,ncenter
    write(ifileid,"(6x,a,3f14.8)") kindname(atmkind(iatm)),a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
end do
write(ifileid,"(a)") "    &END COORD"
if (itask==6) then
    write(ifileid,"(a)") "#   &VELOCITY #You can set initial atomic velocities in this section"
    write(ifileid,"(a)") "#   &END VELOCITY"
end if

!---- &KIND
if (method/="GFN1-xTB".and.method/="PM6".and.method/="SCC-DFTB") then !Semi-empirical methods do not need to define these
    ntime=1
    if (itask==11) ntime=2 !For BSSE task, the second time write same kinds but with _ghost suffix
    do itime=1,ntime
        do ikind=1,nkind
            if (count(atmkind(:)==ikind)==0) cycle
            if (itime==1) then
                write(ifileid,"(a)") "    &KIND "//kindname(ikind)
            else if (itime==2) then
                write(ifileid,"(a)") "    &KIND "//trim(kindname(ikind))//"_ghost"
            end if
            write(ifileid,"(a)") "      ELEMENT "//ind2name(kindeleidx(ikind))
            if (ibas==20) then
                write(ifileid,"(a)") "      BASIS_SET cc-DZ"
                write(ifileid,"(a)") "      BASIS_SET RI_AUX RI_DZ"
            else if (ibas==21) then
                write(ifileid,"(a)") "      BASIS_SET cc-TZ"
                write(ifileid,"(a)") "      BASIS_SET RI_AUX RI_TZ"
            else
                if ((ibas>=-6.and.ibas<=-1).or.(ibas>=1.and.ibas<=5)) then
                    write(c80tmp,"(i3)") Nval(kindeleidx(ikind))
                    write(ifileid,"(a)") "      BASIS_SET "//trim(basname(ibas))//'-q'//trim(adjustl(c80tmp))
                else
                    write(ifileid,"(a)") "      BASIS_SET "//trim(basname(ibas))
                end if
            end if
            if (index(method,"ADMM")/=0) then !Set auxiliary basis set file for ADMM
                write(c80tmp,"(i3)") Nval(kindeleidx(ikind))
                write(ifileid,"(a)") "      BASIS_SET AUX_FIT admm-dzp"//'-q'//trim(adjustl(c80tmp)) !Use ADMM
            end if
            if (iLRIGPW==1) then !Set auxiliary basis set file for LRIGPW
                write(ifileid,"(a)") "      BASIS_SET LRI_AUX LRI-DZVP-MOLOPT-GTH-MEDIUM"
            end if
            if (iGAPW==0) then !GPW
                if (method=="Pade") then
                    write(ifileid,"(a)") "      POTENTIAL GTH-PADE"
                else if (method=="BP") then          
                    write(ifileid,"(a)") "      POTENTIAL GTH-BP"
                else if (method=="BLYP".or.method=="B3LYP") then
                    write(ifileid,"(a)") "      POTENTIAL GTH-BLYP"
                !else if (index(method,"SCAN")/=0) then
                !    write(ifileid,"(a)") "      POTENTIAL GTH-SCAN"
                !else if (index(method,"PBE0")/=0) then
                !    write(ifileid,"(a)") "      POTENTIAL GTH-PBE0"
                else if (index(method,"MP2")/=0) then
                    write(ifileid,"(a)") "      POTENTIAL GTH-HF"
                else if (index(method,"revDSD-PBEP86")/=0) then
                    write(c80tmp,"(i3)") Nval(kindeleidx(ikind))
                    write(ifileid,"(a)") "      POTENTIAL GTH-PBE0"//'-q'//trim(adjustl(c80tmp))
                else                           
                    write(ifileid,"(a)") "      POTENTIAL GTH-PBE"
                end if
            else !GAPW
                write(ifileid,"(a)") "      POTENTIAL ALL"
            end if
            if (itime==2.and.itask==11) write(ifileid,"(a)") "      GHOST"
            if (iDFTplusU==1) then
                ie=kindeleidx(ikind)
                if ((ie>=21.and.ie<=28).or.(ie>=39.and.ie<=46).or.(ie>=57.and.ie<=78).or.ie>=89) then !Only +U for d or f element
                    if (ie/=24.and.ie/=25.and.ie/=42.and.ie/=43.and.ie/=46.and.ie/=75) then !Ignore half-occupied d shell elements
                        write(ifileid,"(a)") "      &DFT_PLUS_U T"
                        write(ifileid,"(a)") "        L 2 #Quantum number of angular momentum the atomic orbitals to +U. 2=d, 3=f"
                        write(ifileid,"(a)") "        U_MINUS_J [eV] 3.9 #Effective on-site Coulomb interaction parameter U(eff) = U - J"
                        write(ifileid,"(a)") "      &END DFT_PLUS_U"
                    end if
                end if
            end if
            if (kindmag(ikind)/=0) write(ifileid,"(a,i3)") "      MAGNETIZATION",kindmag(ikind)
            write(ifileid,"(a)") "    &END KIND"
        end do
    end do
end if

write(ifileid,"(a)") "  &END SUBSYS"

!---- &DFT
write(ifileid,"(/,a)") "  &DFT"
if (method/="GFN1-xTB".and.method/="PM6".and.method/="SCC-DFTB") then
    if (ibas<0) then
        write(ifileid,"(a)") "    BASIS_SET_FILE_NAME  GTH_BASIS_SETS"
    else if (ibas<=5) then
        write(ifileid,"(a)") "    BASIS_SET_FILE_NAME  BASIS_MOLOPT"
    else if (ibas==10.or.ibas==11.or.ibas==12.or.ibas==14) then
        write(ifileid,"(a)") "    BASIS_SET_FILE_NAME  EMSL_BASIS_SETS"
    else if (ibas==13) then
        write(ifileid,"(a)") "    BASIS_SET_FILE_NAME  BASIS_pob"
    else if (ibas==20.or.ibas==21) then
        write(ifileid,"(a)") "    BASIS_SET_FILE_NAME  BASIS_RI_cc-TZ"
    end if
    if (index(method,"ADMM")/=0) then !Set basis set file for ADMM
        write(ifileid,"(a)") "    BASIS_SET_FILE_NAME  BASIS_ADMM_UZH"
    end if
    if (iLRIGPW==1) then !Set basis set file for LRIGPW
        write(ifileid,"(a)") "    BASIS_SET_FILE_NAME  BASIS_LRIGPW_AUXMOLOPT"
    end if
    if (index(method,"MP2")/=0) then
        write(ifileid,"(a)") "    POTENTIAL_FILE_NAME  HF_POTENTIALS"
    else if (index(method,"revDSD-PBEP86")/=0) then
        write(ifileid,"(a)") "    POTENTIAL_FILE_NAME  POTENTIAL_UZH"
    !else if (index(method,"SCAN")/=0.or.index(method,"r2SCAN")/=0.or.index(method,"PBE0")/=0) then
    !    write(ifileid,"(a)") "    POTENTIAL_FILE_NAME  POTENTIAL_UZH"
    else
        write(ifileid,"(a)") "    POTENTIAL_FILE_NAME  POTENTIAL"
    end if
end if
write(ifileid,"(a)") "#   WFN_RESTART_FILE_NAME "//trim(c200tmp)//"-RESTART.wfn"
write(ifileid,"(a,i5,a)") "    CHARGE",netchg," #Net charge"
write(ifileid,"(a,i5,a)") "    MULTIPLICITY",multispin," #Spin multiplicity"
if (multispin>1.or.any(kindmag(1:nkind)/=0)) write(ifileid,"(a)") "    UKS"
if (any(efieldvec/=0)) then
    efieldmag=dsqrt(sum(efieldvec**2))
    if (ifPBC==0) then
        write(ifileid,"(a)") "    &EFIELD"
        write(ifileid,"(a,f8.5)") "      INTENSITY",efieldmag
        write(ifileid,"(a,3f8.5)") "      POLARISATION",efieldvec/efieldmag
        write(ifileid,"(a)") "    &END EFIELD"
    else
        write(ifileid,"(a)") "    &PERIODIC_EFIELD"
        write(ifileid,"(a,f8.5)") "      INTENSITY",efieldmag
        write(ifileid,"(a,3f8.5)") "      POLARISATION",efieldvec/efieldmag
        write(ifileid,"(a)") "    &END PERIODIC_EFIELD"
    end if
end if
if (ikpconvtest==1) then
    write(ifileid,"(a)") "    &KPOINTS"
    write(ifileid,"(a)") "      SCHEME MONKHORST-PACK kp_test"
    write(ifileid,"(a)") "    &END KPOINTS"
else if (ikpoint1/=1.or.ikpoint2/=1.or.ikpoint3/=1) then
    write(ifileid,"(a)") "    &KPOINTS"
    write(ifileid,"(a,3i3)") "      SCHEME MONKHORST-PACK",ikpoint1,ikpoint2,ikpoint3
    write(ifileid,"(a)") "    &END KPOINTS"
end if
if (iDFTplusU==1) write(ifileid,"(a)") "    PLUS_U_METHOD MULLIKEN #The method used in DFT+U. Can also be Lowdin"
if (iTDDFT==1.and.(itask==2.or.itask==3.or.itask==4.or.itask==5.or.itask==6.or.itask==7.or.itask==8)) then !Evaluate force for excited state
    write(ifileid,"(a)") "    &EXCITED_STATES"
    write(ifileid,"(a,i5,a)") "      STATE",istate_force," #For which excited state the force will be evaluated. Negative value indicates state following"
    write(ifileid,"(a)") "    &END EXCITED_STATES"
end if

!---- &QS
write(ifileid,"(a)") "    &QS"

if (iLSSCF==1) write(ifileid,"(a)") "      LS_SCF #Use linear scaling self consistent field method"
!Set proper EPS_SCF (default is 1E-5) and EPS_DEFAULT (default is 1E-10)
if (itask==1) then !Single point, do not need high accuracy
    eps_scf=5D-6
    eps_def=1D-11
    if (iwfc==1) then
        eps_scf=1D-6
        eps_def=1D-12
    end if
else if (itask==5) then !Vibration analysis is realized based on finite difference of force
    eps_scf=1D-7
    eps_def=1D-14
else if (itask==6.or.itask==14) then !For faster MD/PIMD, use even looser threshold
    eps_scf=1D-5
    eps_def=1D-10
else if (itask==9.or.itask==10) then !NMR and polar may need pretty tight threshold
    eps_scf=2D-8
    eps_def=1D-14
else !Other tasks involving energy derivative, use marginally tighter convergence
    eps_scf=1D-6
    eps_def=1D-12
end if
write(ifileid,"(a,1PE8.1,a)") "      EPS_DEFAULT",eps_def," #Set all EPS_xxx to values such that the energy will be correct up to this value"

if (iHFX==1.or.method=="RI-(EXX+RPA)@PBE") then
    write(ifileid,"(a)") "      EPS_PGF_ORB 1E-12 #If warning ""Kohn Sham matrix not 100% occupied"" occurs and meantime calculation is unstable, decrease it"
end if
if (itask==6) then
    write(ifileid,"(a)") "      EXTRAPOLATION ASPC #Extrapolation for wavefunction during e.g. MD. ASPC is default, PS also be used"
    write(ifileid,"(a)") "      EXTRAPOLATION_ORDER 3 #Order for PS or ASPC extrapolation. 3 is default"
end if
if (method=="GFN1-xTB") then
    write(ifileid,"(a)") "      METHOD xTB"
    write(ifileid,"(a)") "      &xTB"
    if (ifPBC>0) write(ifileid,"(a)") "        DO_EWALD T" !Default is Coulomb way to calculate electrostatic interaction
    write(ifileid,"(a)") "        CHECK_ATOMIC_CHARGES F #xTB calculation often crashes without setting this to false"
    write(ifileid,"(a)") "        &PARAMETER"
    write(ifileid,"(a)") "          DISPERSION_PARAMETER_FILE dftd3.dat"
    write(ifileid,"(a)") "          PARAM_FILE_NAME xTB_parameters"
    write(ifileid,"(a)") "        &END PARAMETER"
    write(ifileid,"(a)") "      &END xTB"
else if (method=="PM6") then
    write(ifileid,"(a)") "      METHOD PM6"
    if (ifPBC>0) then
        write(ifileid,"(a)") "      &SE"
        write(ifileid,"(a)") "        PERIODIC EWALD"
        write(ifileid,"(a)") "      &END SE"
    end if
else if (method=="SCC-DFTB") then
    write(ifileid,"(a)") "      METHOD DFTB"
    write(ifileid,"(a)") "      &DFTB"
    write(ifileid,"(a)") "        SELF_CONSISTENT  T"
    write(ifileid,"(a)") "        DISPERSION       T"
    if (ifPBC>0) write(ifileid,"(a)") "        DO_EWALD         T"
    write(ifileid,"(a)") "        &PARAMETER"
    write(ifileid,"(a)") "          PARAM_FILE_PATH  DFTB/scc"
    write(ifileid,"(a)") "          PARAM_FILE_NAME  scc_parameter"
    write(ifileid,"(a)") "          UFF_FORCE_FIELD  uff_table"
    write(ifileid,"(a)") "        &END PARAMETER"
    write(ifileid,"(a)") "      &END DFTB"
else if (iGAPW==1) then
    write(ifileid,"(a)") "      METHOD GAPW"
else !GPW with GTH pseudopotential
    if (iLRIGPW==1) then
        write(ifileid,"(a)") "      METHOD LRIGPW"
        write(ifileid,"(a)") "      &LRIGPW"
        write(ifileid,"(a)") "        LRI_OVERLAP_MATRIX AUTOSELECT #Choose automatically for each pair whether to calculate the inverse or pseudoinverse"
        write(ifileid,"(a)") "      &END LRIGPW"
    end if
end if
write(ifileid,"(a)") "    &END QS"

!---- &POISSON
write(ifileid,"(a)") "    &POISSON" !How to deal with electrostatic part
if (method=="PM6") then !Special for semi-empirical
    write(ifileid,"(a)") "      &EWALD"
    write(ifileid,"(a)") "        &MULTIPOLES"
    write(ifileid,"(a)") "          MAX_MULTIPOLE_EXPANSION QUADRUPOLE"
    write(ifileid,"(a)") "        &END MULTIPOLES"
    write(ifileid,"(a)") "        EWALD_TYPE EWALD"
    !write(ifileid,"(a)") "        ALPHA  0.5" !See e.g. https://github.com/misteliy/cp2k/blob/master/tests/SE/regtest-3/Al2O3.inp
    ngmax1=nint(cellv1(1)*nrep1*b2a) !Make gmax approximately 1 pt/Angstrom in each direction
    if (mod(ngmax1,2)==0) ngmax1=ngmax1+1 !Must be odd number
    ngmax2=nint(cellv2(2)*nrep2*b2a)
    if (mod(ngmax2,2)==0) ngmax2=ngmax2+1 !Must be odd number
    ngmax3=nint(cellv3(3)*nrep3*b2a)
    if (mod(ngmax3,2)==0) ngmax3=ngmax3+1
    write(ifileid,"(a,3i4)") "        GMAX",ngmax1,ngmax2,ngmax3
    write(ifileid,"(a)") "      &END EWALD"
else if (method=="SCC-DFTB") then !Special for DFTB
    write(ifileid,"(a)") "      &EWALD"
    write(ifileid,"(a)") "        EWALD_TYPE SPME"
    ngmax1=2*nint(cellv1(1)*nrep1*b2a) !Make gmax approximately 1 pt/Angstrom in each direction
    ngmax2=2*nint(cellv2(2)*nrep2*b2a)
    ngmax3=2*nint(cellv3(3)*nrep3*b2a)
    write(ifileid,"(a,3i4)") "        GMAX",ngmax1,ngmax2,ngmax3
    write(ifileid,"(a)") "      &END EWALD"
else !Common case
    if (iPSOLVER==1) then
        write(ifileid,"(a)") "      PERIODIC XYZ #Direction(s) of PBC for calculating electrostatics"
    else
        write(ifileid,"(a)") "      PERIODIC "//trim(PBCdir)//" #Direction(s) of PBC for calculating electrostatics"
    end if
    if (iPSOLVER==1) write(ifileid,"(a)") "      PSOLVER PERIODIC #The way to solve Poisson equation"
    if (iPSOLVER==2) write(ifileid,"(a)") "      PSOLVER ANALYTIC #The way to solve Poisson equation"
    if (iPSOLVER==3) write(ifileid,"(a)") "      PSOLVER MT #The way to solve Poisson equation"
    if (iPSOLVER==4) write(ifileid,"(a)") "      PSOLVER WAVELET #The way to solve Poisson equation"
end if
write(ifileid,"(a)") "    &END POISSON"

if (index(method,"ADMM")/=0) then !Use ADMM
    write(ifileid,"(a)") "    &AUXILIARY_DENSITY_MATRIX_METHOD"
    write(ifileid,"(a)") "      METHOD BASIS_PROJECTION #Method used for wavefunction fitting"
    if (iTDDFT==1) then
        write(ifileid,"(a)") "      ADMM_PURIFICATION_METHOD NONE #NONE is the only choice for TDDFT with ADMM"
    else
        write(ifileid,"(a)") "      ADMM_PURIFICATION_METHOD MO_DIAG"
    end if
    if (iTDDFT==1.and.index(method,"_ADMM")/=0) then !Needed otherwise cannot run. Suggested by https://www.cp2k.org/howto:tddft
        write(ifileid,"(a)") "      EXCH_SCALING_MODEL NONE"
        if (method=="B3LYP_ADMM") then
            write(ifileid,"(a)") "      EXCH_CORRECTION_FUNC BECKE88X"
        else
            write(ifileid,"(a)") "      EXCH_CORRECTION_FUNC PBEX"
        end if
    end if
    write(ifileid,"(a)") "    &END AUXILIARY_DENSITY_MATRIX_METHOD"
end if

!---- &XC
if (method=="PM6".or.method=="GFN1-xTB".or.method=="SCC-DFTB") goto 100
write(ifileid,"(a)") "    &XC"
if (index(method,"LIBXC")/=0) then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL"
    if (method=="B97M-rV_LIBXC") then !Non-separable XC
        write(ifileid,"(a)") "        &MGGA_XC_B97M_V"
        write(ifileid,"(a)") "        &END MGGA_XC_B97M_V"
    else !X-C separable
        if (method=="MN15L_LIBXC") then
            write(ifileid,"(a)") "        &MGGA_X_MN15_L"
            write(ifileid,"(a)") "        &END MGGA_X_MN15_L"
            write(ifileid,"(a)") "        &MGGA_C_MN15_L"
            write(ifileid,"(a)") "        &END MGGA_C_MN15_L"
        else if (method=="SCAN_LIBXC") then
            write(ifileid,"(a)") "        &MGGA_X_SCAN"
            write(ifileid,"(a)") "        &END MGGA_X_SCAN"
            write(ifileid,"(a)") "        &MGGA_C_SCAN"
            write(ifileid,"(a)") "        &END MGGA_C_SCAN"
        else if (method=="r2SCAN_LIBXC") then
            write(ifileid,"(a)") "        &MGGA_X_R2SCAN"
            write(ifileid,"(a)") "        &END MGGA_X_R2SCAN"
            write(ifileid,"(a)") "        &MGGA_C_R2SCAN"
            write(ifileid,"(a)") "        &END MGGA_C_R2SCAN"
        else if (method=="RPBE_LIBXC") then
            write(ifileid,"(a)") "        &GGA_X_RPBE"
            write(ifileid,"(a)") "        &END GGA_X_RPBE"
            write(ifileid,"(a)") "        &GGA_C_PBE"
            write(ifileid,"(a)") "        &END GGA_C_PBE"
        else if (method=="revTPSS_LIBXC") then
            write(ifileid,"(a)") "        &MGGA_X_REVTPSS"
            write(ifileid,"(a)") "        &END MGGA_X_REVTPSS"
            write(ifileid,"(a)") "        &MGGA_C_REVTPSS"
            write(ifileid,"(a)") "        &END MGGA_C_REVTPSS"
        end if
    end if
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (method=="GFN1-xTB".or.method=="PM6".or.method=="SCC-DFTB") then
    continue
else if (index(method,"PBE0")/=0) then
        write(ifileid,"(a)") "      &XC_FUNCTIONAL"
        write(ifileid,"(a)") "        &PBE"
        write(ifileid,"(a)") "          SCALE_X 0.75"
        write(ifileid,"(a)") "          SCALE_C 1.0"
        write(ifileid,"(a)") "        &END PBE"
        write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (index(method,"B3LYP")/=0) then
        write(ifileid,"(a)") "      &XC_FUNCTIONAL"
        write(ifileid,"(a)") "        &LYP"
        write(ifileid,"(a)") "          SCALE_C 0.81"
        write(ifileid,"(a)") "        &END"
        write(ifileid,"(a)") "        &BECKE88"
        write(ifileid,"(a)") "          SCALE_X 0.72"
        write(ifileid,"(a)") "        &END"
        write(ifileid,"(a)") "        &VWN"
        write(ifileid,"(a)") "          FUNCTIONAL_TYPE VWN3 #Gaussian's B3LYP definition"
        write(ifileid,"(a)") "          SCALE_C 0.19"
        write(ifileid,"(a)") "        &END"
        write(ifileid,"(a)") "        &XALPHA"
        write(ifileid,"(a)") "          SCALE_X 0.08"
        write(ifileid,"(a)") "        &END"
        write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (method=="revPBE".or.method=="PBEsol") then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL PBE"
    write(ifileid,"(a)") "        &PBE"
    if (method=="revPBE") write(ifileid,"(a)") "          PARAMETRIZATION REVPBE"
    if (method=="PBEsol") write(ifileid,"(a)") "          PARAMETRIZATION PBESOL"
    write(ifileid,"(a)") "        &END PBE"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (index(method,"HSE")/=0) then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL"
    write(ifileid,"(a)") "        &XWPBE"
    write(ifileid,"(a)") "          SCALE_X -0.25"
    write(ifileid,"(a)") "          SCALE_X0 1.0"
    write(ifileid,"(a)") "          OMEGA 0.11"
    write(ifileid,"(a)") "        &END XWPBE"
    write(ifileid,"(a)") "        &PBE"
    write(ifileid,"(a)") "          SCALE_X 0.0"
    write(ifileid,"(a)") "          SCALE_C 1.0"
    write(ifileid,"(a)") "        &END PBE"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (index(method,"MP2")/=0) then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL NONE"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (method=="RI-(EXX+RPA)@PBE".or.method=="GW@PBE") then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL PBE"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (index(method,"B2PLYP")/=0) then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL"
    write(ifileid,"(a)") "        &LYP"
    write(ifileid,"(a)") "          SCALE_C  0.73"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "        &BECKE88"
    write(ifileid,"(a)") "          SCALE_X  0.47"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (index(method,"B2GP-PLYP")/=0) then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL"
    write(ifileid,"(a)") "        &LYP"
    write(ifileid,"(a)") "          SCALE_C  0.64"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "        &BECKE88"
    write(ifileid,"(a)") "          SCALE_X  0.35"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (index(method,"DSD-BLYP")/=0) then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL"
    write(ifileid,"(a)") "        &LYP"
    write(ifileid,"(a)") "          SCALE_C  0.69"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "        &BECKE88"
    write(ifileid,"(a)") "          SCALE_X  0.31"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else if (index(method,"revDSD-PBEP86")/=0) then
    write(ifileid,"(a)") "      &XC_FUNCTIONAL"
    write(ifileid,"(a)") "        &GGA_X_PBE"
    write(ifileid,"(a)") "          SCALE 0.31"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "        &P86C"
    write(ifileid,"(a)") "          SCALE_C 0.4296"
    write(ifileid,"(a)") "        &END"
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
else !Common native GGA functionals
    write(ifileid,"(a)") "      &XC_FUNCTIONAL "//trim(method)
    write(ifileid,"(a)") "      &END XC_FUNCTIONAL"
end if
!HF part for hybrid functionals
if (iHFX==1) then !HFX potential
    write(ifileid,"(a)") "      &HF"
    if (index(method,"PBE0")/=0.or.index(method,"HSE")/=0) write(ifileid,"(a)") "        FRACTION 0.25 #HF composition"
    if (index(method,"B3LYP")/=0) write(ifileid,"(a)") "        FRACTION 0.2 #HF composition"
    if (index(method,"MP2")/=0) write(ifileid,"(a)") "        FRACTION 1.0 #HF composition"
    if (index(method,"B2PLYP")/=0) write(ifileid,"(a)") "        FRACTION 0.53 #HF composition"
    if (index(method,"B2GP-PLYP")/=0) write(ifileid,"(a)") "        FRACTION 0.65 #HF composition"
    if (index(method,"DSD-BLYP")/=0) write(ifileid,"(a)") "        FRACTION 0.69 #HF composition"
    if (index(method,"revDSD-PBEP86")/=0) write(ifileid,"(a)") "        FRACTION 0.69 #HF composition"
    write(ifileid,"(a)") "        &SCREENING"
    if (iwfc==1) then !For wavefunction-based correlation calculation, use tighter threshold for screening
        write(ifileid,"(a)") "          EPS_SCHWARZ 1E-7 #The larger the value, the lower the cost and lower the accuracy"
    else !Hybrid functionals
        write(ifileid,"(a)") "          EPS_SCHWARZ 1E-6 #The larger the value, the lower the cost and lower the accuracy"
    end if
    write(ifileid,"(a)") "          SCREEN_ON_INITIAL_P T #Screening ERI based on initial density matrix, need to provide wavefunction restart file"
    write(ifileid,"(a)") "        &END SCREENING"
    if (index(method,"HSE")/=0) then
        write(ifileid,"(a)") "        &INTERACTION_POTENTIAL"
        write(ifileid,"(a)") "          POTENTIAL_TYPE SHORTRANGE"
        write(ifileid,"(a)") "          OMEGA 0.11"
        write(ifileid,"(a)") "        &END INTERACTION_POTENTIAL"
    else
        if (ifPBC/=0) then !PBC system needs Coulomb truncation for common hybrid functionals
            write(ifileid,"(a)") "        &INTERACTION_POTENTIAL"
            write(ifileid,"(a)") "          POTENTIAL_TYPE TRUNCATED"
            !Set CUTOFF_RADIUS to 1/2.01 of shortest Cartesian length for small cell
            cellv1=cellv1*nrep1
            cellv2=cellv2*nrep2
            cellv3=cellv3*nrep3
            call cellxyzsize(xsize,ysize,zsize)
            trunc_rad=min(min(xsize,ysize),zsize)*b2a/2.01D0
            cellv1=cellv1/nrep1
            cellv2=cellv2/nrep2
            cellv3=cellv3/nrep3
            if (trunc_rad>6) then !If half of shortest box length is larger than 6 Angstrom, simply use 6, this is usually adequate
                write(ifileid,"(a,f8.4)") "          CUTOFF_RADIUS 6.0 #Cutoff radius for truncated 1/r Coulomb operator"
            else
                write(ifileid,"(a,f8.4,a)") "          CUTOFF_RADIUS",trunc_rad," #Cutoff radius for truncated 1/r Coulomb operator"
            end if
            write(ifileid,"(a)") "        &END INTERACTION_POTENTIAL"
        end if
    end if
    write(ifileid,"(a)") "        &MEMORY"
    write(ifileid,"(a)") "          MAX_MEMORY 3000 #Memory(MB) per MPI process for calculating HF exchange"
    !Scaling factor to scale EPS_SCHWARZ. Storage threshold for compression will be EPS_SCHWARZ*EPS_STORAGE_SCALING
    write(ifileid,"(a)") "          EPS_STORAGE_SCALING 0.1"
    write(ifileid,"(a)") "        &END MEMORY"
    write(ifileid,"(a)") "      &END HF"
end if

!Wavefunction-based correlation part
if (iwfc==1) then
    write(ifileid,"(a)") "      &WF_CORRELATION"
    if (method=="RI-(EXX+RPA)@PBE") then
        write(ifileid,"(a)") "        &RI_RPA"
        write(ifileid,"(a)") "          QUADRATURE_POINTS  10  #Number of quadrature points for the numerical integration in the RI-RPA method"
        write(ifileid,"(a)") "          MINIMAX  #Use Minimax quadrature scheme for performing the numerical integration"
        write(ifileid,"(a)") "          &HF"
        write(ifileid,"(a)") "            FRACTION 1.0"
        write(ifileid,"(a)") "            &SCREENING"
        write(ifileid,"(a)") "              EPS_SCHWARZ 1E-7"
        write(ifileid,"(a)") "            &END SCREENING"
        if (ifPBC/=0) then !PBC system needs Coulomb truncation for evaluating HF exchange of RPA energy
            write(ifileid,"(a)") "            &INTERACTION_POTENTIAL"
            write(ifileid,"(a)") "              POTENTIAL_TYPE TRUNCATED"
            !Set CUTOFF_RADIUS to 1/2.01 of shortest Cartesian length for small cell
            cellv1=cellv1*nrep1;cellv2=cellv2*nrep2;cellv3=cellv3*nrep3
            call cellxyzsize(xsize,ysize,zsize)
            trunc_rad=min(min(xsize,ysize),zsize)*b2a/2.01D0
            cellv1=cellv1/nrep1;cellv2=cellv2/nrep2;cellv3=cellv3/nrep3
            if (trunc_rad>6) then !If half of shortest box length is larger than 6 Angstrom, simply use 6, this is usually adequate
                write(ifileid,"(a,f8.4)") "              CUTOFF_RADIUS 6.0 #Cutoff radius for truncated 1/r Coulomb operator"
            else
                write(ifileid,"(a,f8.4,a)") "              CUTOFF_RADIUS",trunc_rad," #Cutoff radius for truncated 1/r Coulomb operator"
            end if
            write(ifileid,"(a)") "            &END INTERACTION_POTENTIAL"
        end if
        write(ifileid,"(a)") "          &END HF"
        write(ifileid,"(a)") "        &END RI_RPA"
    else if (method=="GW@PBE") then
        
    else !Other case, all involve MP2
        write(ifileid,"(a)") "        &RI_MP2"
        write(ifileid,"(a)") "        &END RI_MP2"
        if (index(method,"SCS-MP2")/=0) then
            write(ifileid,"(a)") "        SCALE_S 1.2"
            write(ifileid,"(a)") "        SCALE_T 0.3333333"
        else if (index(method,"B2PLYP")/=0) then
            write(ifileid,"(a)") "        SCALE_S 0.27"
            write(ifileid,"(a)") "        SCALE_T 0.27"
        else if (index(method,"B2GP-PLYP")/=0) then
            write(ifileid,"(a)") "        SCALE_S 0.36"
            write(ifileid,"(a)") "        SCALE_T 0.36"
        else if (index(method,"DSD-BLYP")/=0) then
            write(ifileid,"(a)") "        SCALE_S 0.46"
            write(ifileid,"(a)") "        SCALE_T 0.37"
        else if (index(method,"revDSD-PBEP86")/=0) then
            write(ifileid,"(a)") "        SCALE_S 0.0799"
            write(ifileid,"(a)") "        SCALE_T 0.5785"
        end if
    end if
    write(ifileid,"(a)") "        &INTEGRALS"
    write(ifileid,"(a)") "          &WFC_GPW"
    write(ifileid,"(a)") "            CUTOFF      300"
    write(ifileid,"(a)") "            REL_CUTOFF  50"
    write(ifileid,"(a)") "            EPS_FILTER  1E-12"
    write(ifileid,"(a)") "            EPS_GRID    1E-8"
    write(ifileid,"(a)") "          &END WFC_GPW"
    write(ifileid,"(a)") "        &END INTEGRALS"
    write(ifileid,"(a)") "        MEMORY    3000 #Maximum allowed total memory usage (MB) during wavefunction-based correlation"
    write(ifileid,"(a)") "        GROUP_SIZE  1 #Default. Also known as NUMBER_PROC"
    write(ifileid,"(a)") "      &END WF_CORRELATION"
end if

!--- Dispersion correction
if (idispcorr>0.or.method=="BEEFVDW") then
    write(ifileid,"(a)") "      &VDW_POTENTIAL"
    if (idispcorr==1.or.idispcorr==2) then
        write(ifileid,"(a)") "        POTENTIAL_TYPE PAIR_POTENTIAL"
        write(ifileid,"(a)") "        &PAIR_POTENTIAL"
        write(ifileid,"(a)") "          PARAMETER_FILE_NAME dftd3.dat"
        if (idispcorr==1) write(ifileid,"(a)") "          TYPE DFTD3"
        if (idispcorr==2) write(ifileid,"(a)") "          TYPE DFTD3(BJ)"
        !See qs_dispersion_pairpot.F on how to write functional name
        !Special cases are explicitly list here, used to change name, remove RI-, remove _LIBXC, etc.
        if (method=="BP") then
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL BP86"
        else if (method=="MN15L_LIBXC") then !i.e. Remove _LIBXC suffix
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL MN15L"
        else if (method=="SCAN_LIBXC") then !i.e. Remove _LIBXC suffix
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL SCAN"
        else if (method=="RPBE_LIBXC") then !i.e. Remove _LIBXC suffix
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL RPBE"
        else if (method=="revTPSS_LIBXC") then !i.e. Remove _LIBXC suffix
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL revTPSS"
        else if (index(method,"B2PLYP")/=0) then
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL B2PLYP"
        else if (index(method,"B2GP-PLYP")/=0) then
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL B2GPPLYP"
        else if (index(method,"DSD-BLYP")/=0) then
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL DSD-BLYP"
        else if (index(method,"revDSD-PBEP86")/=0) then
            write(ifileid,"(a)") "          D3BJ_SCALING 0.4377,0,0,5.5 #s6,a1,s8,a2"
        else
            c80tmp=trim(method)
            ipos=index(c80tmp,"_ADMM")
            if (ipos/=0) c80tmp(ipos:ipos+4)="" !Remove _ADMM suffix
            write(ifileid,"(a)") "          REFERENCE_FUNCTIONAL "//trim(c80tmp)
        end if
        !write(ifileid,"(a)") "          R_CUTOFF 10.5835442" !Default DFT-D potential range, cutoff will be 2 times this value
        write(ifileid,"(a)") "          #CALCULATE_C9_TERM T #Calculate C9-related three-body term, more accurate for large system"
        write(ifileid,"(a)") "        &END PAIR_POTENTIAL"
    else if (idispcorr==5.or.method=="BEEFVDW") then         
        write(ifileid,"(a)") "        POTENTIAL_TYPE NON_LOCAL"
        write(ifileid,"(a)") "        &NON_LOCAL"
        if (idispcorr==5) then
            write(ifileid,"(a)") "          TYPE RVV10"
            if (method=="B97M-rV_LIBXC") then !See: Ab initio molecular dynamics simulations of liquid water using high quality meta-GGA functionals
                write(ifileid,"(a)") "          PARAMETERS 6.0 0.01"
            else
                write(ifileid,"(a)") "    #The default rVV10 b and C parameters are given below. They should be replaced by proper values for current XC functional"
                write(ifileid,"(a)") "          PARAMETERS 6.3 9.3E-3"
            end if
            write(ifileid,"(a)") "          KERNEL_FILE_NAME rVV10_kernel_table.dat"
        else if (method=="BEEFVDW") then
            write(ifileid,"(a)") "          TYPE LMKLL"
            write(ifileid,"(a)") "          KERNEL_FILE_NAME vdW_kernel_table.dat"
        end if
        write(ifileid,"(a)") "        &END NON_LOCAL"
    end if
    write(ifileid,"(a)") "      &END VDW_POTENTIAL"
end if
if (ifineXCgrid==1) then
    write(ifileid,"(a)") "      &XC_GRID"
    write(ifileid,"(a)") "        USE_FINER_GRID T #Use finer grid for calculating XC"
    write(ifileid,"(a)") "      &END XC_GRID"
end if
write(ifileid,"(a)") "    &END XC"
100 continue

!--- &MGRID
if (method=="GFN1-xTB".or.method=="PM6".or.method=="SCC-DFTB") then
    continue !Semi-empirical methods do not need to set CUTOFF
else
    write(ifileid,"(a)") "    &MGRID"
    if (iconvtest==1) then
        write(ifileid,"(a)") "      CUTOFF LT_cutoff"
        write(ifileid,"(a)") "      REL_CUTOFF LT_rel_cutoff"
    else
        write(ifileid,"(a,i5)") "      CUTOFF",CUTOFF
        write(ifileid,"(a,i4)") "      REL_CUTOFF",REL_CUTOFF
    end if
    if (ibas==3.or.ibas==4.or.ibas==5) write(ifileid,"(a)") "      NGRIDS 5 #The number of multigrids to use. 5 is optimal for MOLOPT-GTH basis sets"
    write(ifileid,"(a)") "    &END MGRID"
end if

!--- &SCF or &LC_SCF
if (iLSSCF==1) then !&LS_SCF
    write(ifileid,"(a)") "    &LS_SCF"
    write(ifileid,"(a)") "      PURIFICATION_METHOD TRS4  #Scheme used to purify Kohn-Sham matrix into density matrix"
    write(ifileid,"(a)") "      #DYNAMIC_THRESHOLD T  #Should the threshold for the purification be chosen dynamically"
    write(ifileid,"(a)") "      EPS_FILTER 1E-7  #Threshold used to determine sparsity and thus speed and accuracy"
    write(ifileid,"(a)") "      EPS_SCF    5E-6  #Target accuracy for SCF convergence in terms of change of total energy per electron"
    write(ifileid,"(a)") "      MAX_SCF 40  #Maximum number of SCF iteration to be performed"
    write(ifileid,"(a)") "      S_PRECONDITIONER ATOMIC  #Preconditions S with some appropriate form. The default ATOMIC is suitable for most cases"
    write(ifileid,"(a)") "      #MIXING_FRACTION 0.45  #Fraction of mixing new density matrix. A value smaller than the default 0.45 may stablize SCF convergence"
    write(ifileid,"(a)") "      #MU -0.15  #Chemical potential in a.u., does not need to set if using TRS4"
    write(ifileid,"(a)") "    &END LS_SCF"
else !&SCF
    write(ifileid,"(a)") "    &SCF"
    if (iconvtest==1) then
        write(ifileid,"(a)") "      MAX_SCF 1"
    else
        if (idiagOT==1) then !Diagonalization
            write(ifileid,"(a)") "      MAX_SCF 128"
        else !OT, usually use more cycles
            if (iouterSCF==0) write(ifileid,"(a)") "      MAX_SCF 200 #Should be set to a small value (e.g. 20) if enabling outer SCF"
            if (iouterSCF==1) write(ifileid,"(a)") "      MAX_SCF 25 #Maximum number of steps of inner SCF"
        end if
    end if
    if (imixing==1) then
        write(ifileid,"(a)") "      MAX_DIIS 7 #Maximum number of DIIS vectors to be used" !The default 4 is too small
        write(ifileid,"(a)") "      EPS_DIIS 0.3 #Threshold on the convergence to start using DIAG/DIIS" !The default 0.1 is too small
    end if
    if (iouterSCF==0) then
        write(ifileid,"(a,1PE8.1,a)") "      EPS_SCF",eps_scf," #Convergence threshold of density matrix during SCF"
    else if (iouterSCF==1) then
        write(ifileid,"(a,1PE8.1,a)") "      EPS_SCF",eps_scf," #Convergence threshold of density matrix of inner SCF"
    end if
    !if (method=="GFN1-xTB".or.method=="PM6") write(ifileid,"(a)") "      SCF_GUESS MOPAC" !Seems they benefit from this
    write(ifileid,"(a)") "#     SCF_GUESS RESTART #Use wavefunction from WFN_RESTART_FILE_NAME file as initial guess"
    if (idiagOT==1) then
        write(ifileid,"(a)") "      &DIAGONALIZATION"
        write(ifileid,"(a)") "        ALGORITHM STANDARD #Algorithm for diagonalization. DAVIDSON is faster for large systems"
        write(ifileid,"(a)") "      &END DIAGONALIZATION"
    else if (idiagOT==2) then
        write(ifileid,"(a)") "      &OT"
        if (method=="GFN1-xTB".or.method=="PM6".or.method=="SCC-DFTB") then !Semi-empirical cannot use FULL_KINETIC. For large system FULL_SINGLE_INVERSE is the only good choice
            write(ifileid,"(a)") "        PRECONDITIONER FULL_SINGLE_INVERSE"
        else
            if (ncenter<300.or.iGAPW==1) then !GAPW should use FULL_ALL even if the system is large, because it converges much better than FULL_KINETIC according to my test and suggestion by Hutter 
                write(ifileid,"(a)") "        PRECONDITIONER FULL_ALL #Usually best but expensive for large system. Cheaper: FULL_SINGLE_INVERSE and FULL_KINETIC"
            else !GPW for large systems, using FULL_ALL will cause too high cost at the first step
                write(ifileid,"(a)") "        PRECONDITIONER FULL_KINETIC #FULL_SINGLE_INVERSE is also worth to try. FULL_ALL is better but quite expensive for large system"
            end if
        end if
        write(ifileid,"(a)") "        MINIMIZER DIIS #CG is worth to consider in difficult cases" !BROYDEN in fact can also be used, but quite poor!
        write(ifileid,"(a)") "        LINESEARCH 2PNT #1D line search algorithm for CG. 2PNT is default, 3PNT is better but more costly. GOLD is best but very expensive"
        write(ifileid,"(a)") "      &END OT"
        if (iouterSCF==0) then
            write(ifileid,"(a)") "      #Uncomment following lines can enable outer SCF, important for difficult convergence case"
            write(ifileid,"(a)") "      #&OUTER_SCF"
            write(ifileid,"(a)") "      #  MAX_SCF 20 #Maximum number of steps of outer SCF"
            write(ifileid,"(a,1PE8.1,a)") "      #  EPS_SCF",eps_scf," #Convergence threshold of outer SCF"
            write(ifileid,"(a)") "      #&END OUTER_SCF"
        else if (iouterSCF==1) then
            write(ifileid,"(a)") "      &OUTER_SCF"
            write(ifileid,"(a)") "        MAX_SCF 20 #Maximum number of steps of outer SCF"
            write(ifileid,"(a,1PE8.1,a)") "        EPS_SCF",eps_scf," #Convergence threshold of outer SCF"
            write(ifileid,"(a)") "      &END OUTER_SCF"
        end if
    end if
    !Case of using diagonalization
    if (idiagOT==1) then
        !--- &SCF \ &MIXING
        write(ifileid,"(a)") "      &MIXING #How to mix old and new density matrices"
        if (imixing==1) then !PM6 and only use this
            write(ifileid,"(a)") "        METHOD DIRECT_P_MIXING"
            write(ifileid,"(a)") "        ALPHA 0.4 #Default. Mixing 40% of new density matrix with the old one"
        else if (imixing==2) then
            write(ifileid,"(a)") "        METHOD BROYDEN_MIXING #PULAY_MIXING is also a good alternative"
            write(ifileid,"(a)") "        ALPHA 0.4 #Default. Mixing 40% of new density matrix with the old one"
            write(ifileid,"(a)") "        NBROYDEN 8 #Default is 4. Number of previous steps stored for the actual mixing scheme" !Equivalent to NBUFFER
        else if (imixing==3) then
            write(ifileid,"(a)") "        METHOD PULAY_MIXING #BROYDEN_MIXING is also a good alternative"
            write(ifileid,"(a)") "        NPULAY 8 #Default is 4. Number of previous steps stored for the actual mixing scheme" !Equivalent to NBUFFER
        end if
        write(ifileid,"(a)") "      &END MIXING"
        !--- &SCF \ &SMEAR
        if (ismear==1) then
            write(ifileid,"(a)") "      &SMEAR"
            write(ifileid,"(a)") "        METHOD FERMI_DIRAC" !Can also be ENERGY_WINDOW, LIST
            write(ifileid,"(a)") "        ELECTRONIC_TEMPERATURE 300 #Electronic temperature of Fermi-Dirac smearing in K"
            write(ifileid,"(a)") "      &END SMEAR"
        end if
        if (nADDED_MOS/=0) then
            if (multispin>1.or.any(kindmag(1:nkind)/=0)) then
                write(ifileid,"(a,i6,i6,a)") "      ADDED_MOS",nADDED_MOS,nADDED_MOS," #Number of virtual MOs to solve for alpha and beta spins"
            else
                write(ifileid,"(a,i6,a)") "      ADDED_MOS",nADDED_MOS," #Number of virtual MOs to solve"
            end if
        end if
    end if
    !--- &SCF \ &PRINT
    write(ifileid,"(a)") "      &PRINT"
    if (itask==5.or.itask==6.or.itask==14) then !Freq, MD, PIMD
        write(ifileid,"(a)") "        &RESTART OFF #Do not generate wfn file to suppress meaningless I/O cost"
        write(ifileid,"(a)") "        &END RESTART"
    else
        write(ifileid,"(a)") "        &RESTART #Note: Use ""&RESTART OFF"" can prevent generating .wfn file"
        write(ifileid,"(a)") "          BACKUP_COPIES 0 #Maximum number of backup copies of wfn file. 0 means never"
        write(ifileid,"(a)") "        &END RESTART"
    end if
    if (itask==6.and.inoSCFinfo==1) then !MD
        write(ifileid,"(a)") "        &PROGRAM_RUN_INFO"
        write(ifileid,"(a)") "          &EACH"
        write(ifileid,"(a)") "            MD 0 #Frequency of printing SCF process during MD. 0 means never"
        write(ifileid,"(a)") "          &END EACH"
        write(ifileid,"(a)") "        &END PROGRAM_RUN_INFO"
    end if
    write(ifileid,"(a)") "      &END PRINT"
    write(ifileid,"(a)") "    &END SCF"
end if

!--- &PRINT of DFT level, FORCE_EVAL/DFT/PRINT
if (imolden==1.or.ioutSbas==1.or.ioutcube>0.or.iatomcharge>0.or.itask==5.or.imoment==1.or.ioutorbene==1) then
    write(ifileid,"(a)") "    &PRINT"
    if (ioutSbas==1) then
        write(ifileid,"(a)") "      &S_CSR_WRITE #Exporting .csr file containing overlap matrix"
        write(ifileid,"(a)") "        REAL_SPACE T #Print the overlap matrix in real-space instead of k-space"
        write(ifileid,"(a)") "      &END S_CSR_WRITE"
    end if
    if (ioutorbene==1) then
        write(ifileid,"(a)") "      &MO"
        write(ifileid,"(a)") "        ENERGIES T"
        write(ifileid,"(a)") "        OCCUPATION_NUMBERS T"
        write(ifileid,"(a)") "        &EACH"
        write(ifileid,"(a)") "          QS_SCF 0"
        write(ifileid,"(a)") "        &END EACH"
        write(ifileid,"(a)") "      &END MO"
    end if
    if (imolden==1) then
        write(ifileid,"(a)") "      &MO_MOLDEN #Exporting .molden file containing wavefunction information"
        write(ifileid,"(a)") "        NDIGITS 9 #Output orbital coefficients if absolute value is larger than 1E-9"
        write(ifileid,"(a)") "      &END MO_MOLDEN"
    end if
    if (ioutcube>0) then
        if (ioutcube==1) then
            write(ifileid,"(a)") "      &E_DENSITY_CUBE"
            write(ifileid,"(a)") "        STRIDE 1 #Stride of exported cube file"
            write(ifileid,"(a)") "      &END E_DENSITY_CUBE"
        else if (ioutcube==2) then
            write(ifileid,"(a)") "      &ELF_CUBE"
            write(ifileid,"(a)") "        STRIDE 1 #Stride of exported cube file"
            write(ifileid,"(a)") "      &END ELF_CUBE"
        else if (ioutcube==3) then
            write(ifileid,"(a)") "      &V_XC_CUBE"
            write(ifileid,"(a)") "        STRIDE 1 #Stride of exported cube file"
            write(ifileid,"(a)") "      &END V_XC_CUBE"
        else if (ioutcube==4) then
            write(ifileid,"(a)") "      &V_HARTREE_CUBE"
            write(ifileid,"(a)") "        STRIDE 1 #Stride of exported cube file"
            write(ifileid,"(a)") "      &END V_HARTREE_CUBE"
        else if (ioutcube==5) then
            write(ifileid,"(a)") "      &EFIELD_CUBE"
            write(ifileid,"(a)") "        STRIDE 1 #Stride of exported cube file"
            write(ifileid,"(a)") "      &END EFIELD_CUBE"
        else if (ioutcube==6) then
            write(ifileid,"(a)") "      &MO_CUBES"
            write(ifileid,"(a)") "        STRIDE 1 #Stride of exported cube file"
            write(ifileid,"(a,i6)") "        NHOMO",NHOMO
            write(ifileid,"(a,i6)") "        NLUMO",NLUMO
            write(ifileid,"(a)") "      &END MO_CUBES"
        end if
    end if
    if (iatomcharge>0) then
        if (iatomcharge==1) then
            write(ifileid,"(a)") "      &MULLIKEN"
            write(ifileid,"(a)") "        PRINT_ALL F #If T, then printing full net AO and overlap population matrix"
            write(ifileid,"(a)") "      &END MULLIKEN"
        else if (iatomcharge==2) then
            write(ifileid,"(a)") "      &LOWDIN"
            write(ifileid,"(a)") "        PRINT_ALL F #If T, then printing full net AO and overlap population matrix"
            write(ifileid,"(a)") "      &END LOWDIN"
        else if (iatomcharge==3.or.iatomcharge==4) then
            write(ifileid,"(a)") "      &HIRSHFELD"
            write(ifileid,"(a)") "        SHAPE_FUNCTION DENSITY"
            if (iatomcharge==4) write(ifileid,"(a)") "        SELF_CONSISTENT T"
            write(ifileid,"(a)") "      &END HIRSHFELD"
        else if (iatomcharge==5) then
            write(ifileid,"(a)") "      &VORONOI"
            write(ifileid,"(a)") "        VORONOI_RADII Covalent" !Better than default of using vdW radii
            write(ifileid,"(a)") "      &END VORONOI"
        end if
    end if
    if (itask==5.or.imoment==1) then
        write(ifileid,"(a)") "      &MOMENTS"
        if (PBCdir=="NONE") then
            write(ifileid,"(a)") "        PERIODIC F #Use Berry phase formula (T) or simple operator (F), the latter normally applies to isolated systems"
            write(ifileid,"(a)") "        MAGNETIC F #If calculating magnetic moments. Can only be used for isolated system"
            write(ifileid,"(a)") "        REFERENCE COM #Reference point for calculating electric moment. COM=center of mass"
        else
            write(ifileid,"(a)") "        PERIODIC T #Use Berry phase formula (T) or simple operator (F), the latter normally applies to isolated systems"
        end if
        if (itask==13) then
            write(ifileid,"(a)") "        FILENAME electric"
            write(ifileid,"(a)") "        COMMON_ITERATION_LEVELS 5"
        end if
        write(ifileid,"(a)") "      &END MOMENTS"
    end if
    if (iSCCS==1) then
        write(ifileid,"(a)") "      @IF 1 #Printing SCCS information in each SCF iteration" !When print level is medium, will also print SCCS iteration information
        write(ifileid,"(a)") "      &SCCS"
        write(ifileid,"(a)") "        &EACH"
        write(ifileid,"(a)") "          QS_SCF 1"
        write(ifileid,"(a)") "        &END EACH"
        write(ifileid,"(a)") "      &END SCCS"
        write(ifileid,"(a)") "      @ENDIF"
    end if
    write(ifileid,"(a)") "    &END PRINT"
end if
if (iSCCS==1) then
    write(ifileid,"(a)") "    &SCCS"
    write(ifileid,"(a)") "      ALPHA [N*m^-1] 0.0"
    write(ifileid,"(a)") "      BETA [kbar] 0.0"
    write(ifileid,"(a)") "      GAMMA [mN/m] 0.0"
    write(ifileid,"(a)") "      DIELECTRIC_CONSTANT 78.36 #Water"
    write(ifileid,"(a)") "      EPS_SCCS 1E-6 #Default. Requested accuracy for the SCCS iteration cycle"
    write(ifileid,"(a)") "      EPS_SCF 0.5 #Default. SCCS iteration is activated only if SCF iteration is converged to this threshold"
    write(ifileid,"(a)") "      MAX_ITER 100 #Default. Maximum number of SCCS iteration steps"
    write(ifileid,"(a)") "      DERIVATIVE_METHOD FFT #Default. Method for calculation of numerical derivatives. Can also be CD3, CD5, CD7"
    write(ifileid,"(a)") "      &ANDREUSSI"
    write(ifileid,"(a)") "        RHO_MAX 0.0035 #Default"
    write(ifileid,"(a)") "        RHO_MIN 0.0001 #Default"
    write(ifileid,"(a)") "      &END ANDREUSSI"
    write(ifileid,"(a)") "    &END SCCS"
end if
if (idipcorr>0) then
    write(ifileid,"(a)") "    SURFACE_DIPOLE_CORRECTION T"
    if (idipcorr==1) write(ifileid,"(a)") "    SURF_DIP_DIR X"
    if (idipcorr==2) write(ifileid,"(a)") "    SURF_DIP_DIR Y"
    if (idipcorr==3) write(ifileid,"(a)") "    SURF_DIP_DIR Z"
end if
if (itask==13) then !Real-time propagation
    write(ifileid,"(a)") "    &REAL_TIME_PROPAGATION"
    write(ifileid,"(a)") "      INITIAL_WFN SCF_WFN #Initial wavefunction used for propagation is obtained by SCF. Can also be RESTART_WFN and RT_RESTART"
    write(ifileid,"(a)") "      EPS_ITER 1E-7 #Convergence criterion for the self consistent propagator loop. This is default value"
    write(ifileid,"(a)") "      MAX_ITER 50 #Maximal number of iterations for the self consistent propagator loop"
    write(ifileid,"(a)") "      APPLY_DELTA_PULSE #Applying a delta kick to the initial wavefunction"
    write(ifileid,"(a)") "      DELTA_PULSE_DIRECTION 0 0 1 #Direction of the applied electric field"
    write(ifileid,"(a)") "      &PRINT"
    write(ifileid,"(a)") "        &RESTART"
    write(ifileid,"(a)") "          BACKUP_COPIES 0 #Never generate backed up .rtpwfn files"
    write(ifileid,"(a)") "        &END RESTART"
    write(ifileid,"(a)") "      &END PRINT"
    write(ifileid,"(a)") "    &END REAL_TIME_PROPAGATION"
end if
write(ifileid,"(a)") "  &END DFT"

if (itask==2.or.inoSCFinfo==1) then !FORCE_EVAL/PRINT
    write(ifileid,"(a)") "  &PRINT"
    if (itask==2) then
        write(ifileid,"(a)") "    &FORCES ON"
        write(ifileid,"(a)") "    &END FORCES"
    end if
    if (inoSCFinfo==1) then
        write(ifileid,"(a)") "    &PROGRAM_RUN_INFO"
        write(ifileid,"(a)") "      &EACH"
        write(ifileid,"(a)") "        MD 0 #Frequency of printing evaluated energies during MD. 0 means never"
        write(ifileid,"(a)") "      &END EACH"
        write(ifileid,"(a)") "    &END PROGRAM_RUN_INFO"
    end if
    write(ifileid,"(a)") "  &END PRINT"
end if

if (itask==4.or.ibarostat>0) write(ifileid,"(a)") "  STRESS_TENSOR ANALYTICAL #Compute full stress tensor analytically" !By default not compute
if (itask==9.or.itask==10.or.iTDDFT==1) then !NMR, polar, TDDFT
    write(ifileid,"(a)") "  &PROPERTIES"
    if (itask==9.or.itask==10) then !NMR, polar
        write(ifileid,"(a)") "    &LINRES #Activate linear response calculation"
        if (ncenter>500) then
            write(ifileid,"(a)") "      PRECONDITIONER FULL_KINETIC #Preconditioner to be used with all minimization schemes"
        else
            write(ifileid,"(a)") "      PRECONDITIONER FULL_ALL #Preconditioner to be used with all minimization schemes"
        end if
        if (itask==9) then !NMR will do response calculation 3*Natoms times, quite expensive, so use relatively looser criterion
            write(ifileid,"(a)") "      EPS 1E-8 #Target accuracy for the convergence of the conjugate gradient" !Tigher than default 1E-6
        else
            write(ifileid,"(a)") "      EPS 1E-10 #Target accuracy for the convergence of the conjugate gradient" !Tigher than default 1E-6
        end if
        write(ifileid,"(a)") "      MAX_ITER 300 #Maximum number of conjugate gradient iteration to be performed for one optimization"
        if (itask==9) then
            write(ifileid,"(a)") "      &CURRENT"
            !write(ifileid,"(a)") "        GAUGE R_AND_STEP_FUNCTION #Default. The gauge used to compute the induced current within GAPW" !I found the default is the only useful choice
            write(ifileid,"(a)") "        ORBITAL_CENTER WANNIER #The orbital center. Can also be ATOM"
            write(ifileid,"(a)") "      &END CURRENT"
            write(ifileid,"(a)") "      &LOCALIZE"
            !I found CRAZY doesn't properly work, so do not explicitly mention, just use default JACOBI
            !write(ifileid,"(a)") "	      METHOD JACOBI #Localization optimization method. JACOBI=2x2 orbital rotations. CRAZY is less robust but usually much faster"
            write(ifileid,"(a)") "	      MAX_ITER 300 #Maximum number of iterations used for localization methods"
            !I found BOYS and PIPEK do not work, so do not explicitly mention them, just use default BERRY
            !write(ifileid,"(a)") "	      OPERATOR BERRY #The quantity to be minimized in localization. Can also be BOYS and PIPEK"
            write(ifileid,"(a)") "      &END LOCALIZE"
            write(ifileid,"(a)") "      &NMR"
            write(ifileid,"(a)") "      #  NICS T #Calculate NICS"
            write(ifileid,"(a)") "      #  NICS_FILE_NAME filepath #Path of the file containing NICS points coordinates"
            write(ifileid,"(a)") "      &END NMR"
        else if (itask==10) then
            write(ifileid,"(a)") "      &POLAR ON"
            write(ifileid,"(a)") "      &END POLAR"
        end if
        write(ifileid,"(a)") "    &END LINRES"
    end if
    if (iTDDFT==1) then !TDDFT
        write(ifileid,"(a)") "    &TDDFPT #TDDFT calculation with Tamm-Dancoff approximation"
        write(ifileid,"(a,i5,a)") "      NSTATES",nstates_TD," #Number of excited states to solve"
        if (isTDA==1) then
            write(ifileid,"(a)") "      KERNEL STDA #Using sTDA approximation"
            write(ifileid,"(a)") "      &STDA"
            if (ifPBC>0) write(ifileid,"(a)") "        DO_EWALD .T. #Use Ewald type method for periodic Coulomb interaction"
            write(ifileid,"(a)") "        FRACTION 0.2 #Fraction of TB Hartree-Fock exchange to use in the kernel"
            write(ifileid,"(a)") "      &END STDA"
        end if
        if (iTDtriplet==1) then
            write(ifileid,"(a)") "      RKS_TRIPLETS .T. #If calculating triplet rather than singlet excited states"
        else
            write(ifileid,"(a,i5)") "      RKS_TRIPLETS .F. #If calculating triplet rather than singlet excited states"
        end if
        write(ifileid,"(a)") "      CONVERGENCE [eV] 1E-4 #Convergence criterion of all excitation energies"
        write(ifileid,"(a)") "      MIN_AMPLITUDE 0.01 #The smallest excitation amplitude to print"
        write(ifileid,"(a)") "#     RESTART .T. #If restarting TDDFT calculation. If true, WFN_RESTART_FILE_NAME should be set to previous .tdwfn file"
        write(ifileid,"(a)") "#     WFN_RESTART_FILE_NAME "//trim(c200tmp)//"-RESTART.tdwfn"
        !if (ifPBC==3) then !The default VELOCITY is also correct for both periodic and isolated systems
        !    write(ifileid,"(a)") "      &DIPOLE_MOMENTS"
        !    write(ifileid,"(a)") "        DIPOLE_FORM BERRY"
        !    write(ifileid,"(a)") "      &END DIPOLE_MOMENTS"
        !end if
        if (index(method,"_ADMM")/=0.and.(itask==2.or.itask==3.or.itask==4.or.itask==5.or.itask==6.or.itask==7.or.itask==8)) then !Need force
            write(ifileid,"(a)") "      ADMM_KERNEL_CORRECTION_SYMMETRIC T"
        end if
        if (index(method,"_ADMM")==0) then !When ADMM is used for TDDFT, &XC should not appear, otherwise error shows: "ADMM is not implemented for a TDDFT kernel XC-functional which is different from the one used for the ground-state calculation"
            write(ifileid,"(a)") "      &XC"
            write(ifileid,"(a)") "        &XC_GRID"
            write(ifileid,"(a)") "          XC_DERIV SPLINE2_SMOOTH #The method used to compute the derivatives"
            write(ifileid,"(a)") "        &END XC_GRID"
            !XC functional for TDDFT
            if (index(method,"LIBXC")/=0) then
                write(ifileid,"(a)") "        &XC_FUNCTIONAL"
                if (method=="B97M-rV_LIBXC") then !Non-separable XC
                    write(ifileid,"(a)") "          &MGGA_XC_B97M_V"
                    write(ifileid,"(a)") "          &END MGGA_XC_B97M_V"
                else !X-C separable
                    if (method=="MN15L_LIBXC") then
                        write(ifileid,"(a)") "          &MGGA_X_MN15_L"
                        write(ifileid,"(a)") "          &END MGGA_X_MN15_L"
                        write(ifileid,"(a)") "          &MGGA_C_MN15_L"
                        write(ifileid,"(a)") "          &END MGGA_C_MN15_L"
                    else if (method=="SCAN_LIBXC") then
                        write(ifileid,"(a)") "          &MGGA_X_SCAN"
                        write(ifileid,"(a)") "          &END MGGA_X_SCAN"
                        write(ifileid,"(a)") "          &MGGA_C_SCAN"
                        write(ifileid,"(a)") "          &END MGGA_C_SCAN"
                    else if (method=="r2SCAN_LIBXC") then
                        write(ifileid,"(a)") "          &MGGA_X_R2SCAN"
                        write(ifileid,"(a)") "          &END MGGA_X_R2SCAN"
                        write(ifileid,"(a)") "          &MGGA_C_R2SCAN"
                        write(ifileid,"(a)") "          &END MGGA_C_R2SCAN"
                    else if (method=="RPBE_LIBXC") then
                        write(ifileid,"(a)") "          &GGA_X_RPBE"
                        write(ifileid,"(a)") "          &END GGA_X_RPBE"
                        write(ifileid,"(a)") "          &GGA_C_PBE"
                        write(ifileid,"(a)") "          &END GGA_C_PBE"
                    else if (method=="revTPSS_LIBXC") then
                        write(ifileid,"(a)") "          &MGGA_X_REVTPSS"
                        write(ifileid,"(a)") "          &END MGGA_X_REVTPSS"
                        write(ifileid,"(a)") "          &MGGA_C_REVTPSS"
                        write(ifileid,"(a)") "          &END MGGA_C_REVTPSS"
                    end if
                end if
                write(ifileid,"(a)") "        &END XC_FUNCTIONAL"
            else if (index(method,"PBE0")/=0) then
                write(ifileid,"(a)") "        &XC_FUNCTIONAL PBE0"
                write(ifileid,"(a)") "        &END XC_FUNCTIONAL"
            else if (index(method,"B3LYP")/=0) then
                write(ifileid,"(a)") "        &XC_FUNCTIONAL B3LYP"
                write(ifileid,"(a)") "        &END XC_FUNCTIONAL"
            else if (method=="revPBE".or.method=="PBEsol") then
                write(ifileid,"(a)") "        &XC_FUNCTIONAL PBE"
                write(ifileid,"(a)") "          &PBE"
                if (method=="revPBE") write(ifileid,"(a)") "          PARAMETRIZATION REVPBE"
                if (method=="PBEsol") write(ifileid,"(a)") "          PARAMETRIZATION PBESOL"
                write(ifileid,"(a)") "          &END PBE"
                write(ifileid,"(a)") "        &END XC_FUNCTIONAL"
            else if (index(method,"HSE")/=0) then
                write(ifileid,"(a)") "        &XC_FUNCTIONAL"
                write(ifileid,"(a)") "          &XWPBE"
                write(ifileid,"(a)") "            SCALE_X -0.25"
                write(ifileid,"(a)") "            SCALE_X0 1.0"
                write(ifileid,"(a)") "            OMEGA 0.11"
                write(ifileid,"(a)") "          &END XWPBE"
                write(ifileid,"(a)") "          &PBE"
                write(ifileid,"(a)") "            SCALE_X 0.0"
                write(ifileid,"(a)") "            SCALE_C 1.0"
                write(ifileid,"(a)") "          &END PBE"
                write(ifileid,"(a)") "        &END XC_FUNCTIONAL"
            else !Common native GGA functionals
                write(ifileid,"(a)") "        &XC_FUNCTIONAL "//trim(method)
                write(ifileid,"(a)") "        &END XC_FUNCTIONAL"
            end if
            write(ifileid,"(a)") "      &END XC"
        end if
        if (iNTO==1) then
            write(ifileid,"(a)") "      &PRINT"
            write(ifileid,"(a)") "        &NTO_ANALYSIS ON #Do NTO analysis for all excited states"
            !write(ifileid,"(a)") "          FILENAME NTO" !Seems not to affect name of actually exported .molden file
            write(ifileid,"(a)") "        &END NTO_ANALYSIS"
            write(ifileid,"(a)") "        &MOS_MOLDEN #Output .molden file containing NTO of the ""NSTATES""th state"
            write(ifileid,"(a)") "          NDIGITS 8"
            write(ifileid,"(a)") "          FILENAME NTO #Filename of NTO .molden file"
            write(ifileid,"(a)") "        &END MOS_MOLDEN"
            write(ifileid,"(a)") "      &END PRINT"
        end if
        write(ifileid,"(a)") "    &END TDDFPT"
    end if
    write(ifileid,"(a)") "  &END PROPERTIES"
end if
if (itask==11) then !BSSE setting
    write(ifileid,"(/,a)") "  &BSSE"
    write(ifileid,"(a)") "    &FRAGMENT"
    call outCP2K_LIST(ifileid,frag1,nfrag1,"      ")
    write(ifileid,"(a)") "    &END FRAGMENT"
    write(ifileid,"(a)") "    &FRAGMENT"
    call outCP2K_LIST(ifileid,frag2,nfrag2,"      ")
    write(ifileid,"(a)") "    &END FRAGMENT"
    write(ifileid,"(a)") "    &CONFIGURATION # real(A)+real(B)"
    write(ifileid,"(a)") "      GLB_CONF 1 1"
    write(ifileid,"(a)") "      SUB_CONF 1 1"
    write(ifileid,"(a,i3)") "      CHARGE",frag1chg+frag2chg
    write(ifileid,"(a,i3)") "      MULTIPLICITY ",totalmulti
    write(ifileid,"(a)") "    &END CONFIGURATION"
    write(ifileid,"(a)") "    &CONFIGURATION # real(A)"
    write(ifileid,"(a)") "      GLB_CONF 1 0"
    write(ifileid,"(a)") "      SUB_CONF 1 0"
    write(ifileid,"(a,i3)") "      CHARGE",frag1chg
    write(ifileid,"(a,i3)") "      MULTIPLICITY ",frag1multi
    write(ifileid,"(a)") "    &END CONFIGURATION"
    write(ifileid,"(a)") "    &CONFIGURATION # real(B)"
    write(ifileid,"(a)") "      GLB_CONF 0 1"
    write(ifileid,"(a)") "      SUB_CONF 0 1"
    write(ifileid,"(a,i3)") "      CHARGE",frag2chg
    write(ifileid,"(a,i3)") "      MULTIPLICITY ",frag2multi
    write(ifileid,"(a)") "    &END CONFIGURATION"
    write(ifileid,"(a)") "    &CONFIGURATION # real(A)+ghost(B)"
    write(ifileid,"(a)") "      GLB_CONF 1 1"
    write(ifileid,"(a)") "      SUB_CONF 1 0"
    write(ifileid,"(a,i3)") "      CHARGE",frag1chg
    write(ifileid,"(a,i3)") "      MULTIPLICITY ",frag1multi
    write(ifileid,"(a)") "    &END CONFIGURATION"
    write(ifileid,"(a)") "    &CONFIGURATION # ghost(A)+real(B)"
    write(ifileid,"(a)") "      GLB_CONF 1 1"
    write(ifileid,"(a)") "      SUB_CONF 0 1"
    write(ifileid,"(a,i3)") "      CHARGE",frag2chg
    write(ifileid,"(a,i3)") "      MULTIPLICITY ",frag2multi
    write(ifileid,"(a)") "    &END CONFIGURATION"
    write(ifileid,"(a)") "  &END BSSE"
end if
if (iatomcharge==6.or.iatomcharge==7) then !Atomic charges
    write(ifileid,"(/,a)") "  &PROPERTIES"
    write(ifileid,"(a)") "    &RESP"
    if (iatomcharge==7) write(ifileid,"(a)") "    USE_REPEAT_METHOD T"
    write(ifileid,"(a)") "      &SPHERE_SAMPLING"
    write(ifileid,"(a)") "        AUTO_VDW_RADII_TABLE CAMBRIDGE #vdW radii type. This is default. Can also be UFF"
    write(ifileid,"(a)") "        AUTO_RMIN_SCALE 1.0 #Scaled factor of vdW radii determining the inner boundary of sampling"
    write(ifileid,"(a)") "        AUTO_RMAX_SCALE 2.0 #Scaled factor of vdW radii determining the outer boundary of sampling"
    write(ifileid,"(a)") "      &END SPHERE_SAMPLING"
    if (ifPBC>0) then
        write(ifileid,"(a)") "      #Uncomment following lines can use slab sampling of fitting points"
        write(ifileid,"(a)") "      #&SLAB_SAMPLING #The fitting points will sampled above a slab"
        write(ifileid,"(a)") "      #  RANGE 2.0 4.0"
        write(ifileid,"(a)") "      #  LENGTH 3.0"
        write(ifileid,"(a)") "      #  ATOM_LIST 1..32 #List of considered atoms"
        write(ifileid,"(a)") "      #  SURF_DIRECTION Z #What above the surface means. Can also be e.g. -Z, X, Z..."
        write(ifileid,"(a)") "      #&END SLAB_SAMPLING"
    end if
    write(ifileid,"(a)") "      &PRINT"
    write(ifileid,"(a)") "        &COORD_FIT_POINTS"
    write(ifileid,"(a)") "        &END COORD_FIT_POINTS"
    write(ifileid,"(a)") "        &RESP_CHARGES_TO_FILE"
    write(ifileid,"(a)") "        &END RESP_CHARGES_TO_FILE"
    write(ifileid,"(a)") "      &END PRINT"
    write(ifileid,"(a)") "    &END RESP"
    write(ifileid,"(a)") "  &END PROPERTIES"
end if
write(ifileid,"(a)") "&END FORCE_EVAL"

!--- &MOTION
if (itask==3.or.itask==4.or.itask==5.or.itask==6.or.itask==7.or.itask==13.or.itask==14) then
    write(ifileid,"(/,a)") "&MOTION"
    if (itask==3.or.itask==7) then !Optimizing atoms for minimum or TS
        write(ifileid,"(a)") "  &GEO_OPT"
        if (itask==3) then
            write(ifileid,"(a)") "    TYPE MINIMIZATION #Search for minimum"
            write(ifileid,"(a)") "    KEEP_SPACE_GROUP F #If T, then space group will be detected and preserved"
            if (ioptmethod==1) then
                write(ifileid,"(a)") "    OPTIMIZER BFGS #Can also be CG (more robust for difficult cases) or LBFGS"
                write(ifileid,"(a)") "    &BFGS"
                write(ifileid,"(a)") "      TRUST_RADIUS 0.2 #Trust radius (maximum stepsize) in Angstrom"
                write(ifileid,"(a)") "#     RESTART_HESSIAN T #If read initial Hessian, uncomment this line and specify the file in the next line"
                write(ifileid,"(a)") "#     RESTART_FILE_NAME to_be_specified"
                write(ifileid,"(a)") "    &END BFGS"
            else if (ioptmethod==2) then
                write(ifileid,"(a)") "    OPTIMIZER LBFGS #Can also be CG (more robust for difficult cases) or BFGS"
                write(ifileid,"(a)") "    &LBFGS"
                write(ifileid,"(a)") "      TRUST_RADIUS 0.2 #Trust radius (maximum stepsize) in Angstrom"
                write(ifileid,"(a)") "      MAX_H_RANK 5 #Larger values (e.g. 30) will accelerate convergence behaviour at the cost of a larger memory consumption"
                write(ifileid,"(a)") "    &END LBFGS"
            else if (ioptmethod==3) then
                write(ifileid,"(a)") "    OPTIMIZER CG #Can also be BFGS or LBFGS"
                write(ifileid,"(a)") "    &CG"
                write(ifileid,"(a)") "      &LINE_SEARCH"
                write(ifileid,"(a)") "        TYPE 2PNT #Two-point extrapolation, cheap while acceptable. Can also be FIT, GOLD"
                write(ifileid,"(a)") "      &END LINE_SEARCH"
                write(ifileid,"(a)") "    &END CG"
            end if
        else if (itask==7) then
            write(ifileid,"(a)") "    TYPE TRANSITION_STATE #Optimizing TS using dimer algorithm"
            write(ifileid,"(a)") "    OPTIMIZER CG" !CG is the only choice for dimer
            write(ifileid,"(a)") "    &CG"
            write(ifileid,"(a)") "      &LINE_SEARCH"
            write(ifileid,"(a)") "        TYPE 2PNT"
            write(ifileid,"(a)") "      &END LINE_SEARCH"
            write(ifileid,"(a)") "    &END CG"
            write(ifileid,"(a)") "    &TRANSITION_STATE"
            write(ifileid,"(a)") "      &DIMER"
            write(ifileid,"(a)") "        DR 0.01 #Default. DR parameter"
            write(ifileid,"(a)") "        ANGLE_TOLERANCE [deg] 4.0 #Tolerance angle for line search performed to optimize dimer orientation"
            write(ifileid,"(a)") "        &ROT_OPT #How to optimizing dimer rotation"
            write(ifileid,"(a)") "          OPTIMIZER CG"
            write(ifileid,"(a)") "          MAX_ITER 50 #Maximum number of optimization steps, default is 200"
            write(ifileid,"(a)") "          #The following thresholds of dimer orientation are the default ones"
            write(ifileid,"(a)") "          MAX_DR 3E-3 #Maximum geometry change"
            write(ifileid,"(a)") "          RMS_DR 1.5E-3 #RMS geometry change"
            write(ifileid,"(a)") "          MAX_FORCE 4.5E-4 #Maximum force"
            write(ifileid,"(a)") "          RMS_FORCE 3E-4 #RMS force"
            write(ifileid,"(a)") "          &CG"
            write(ifileid,"(a)") "            &LINE_SEARCH"
            write(ifileid,"(a)") "              TYPE 2PNT"
            write(ifileid,"(a)") "            &END LINE_SEARCH"
            write(ifileid,"(a)") "          &END CG"
            write(ifileid,"(a)") "        &END ROT_OPT"
            write(ifileid,"(a)") "      &END DIMER"
            write(ifileid,"(a)") "    &END TRANSITION_STATE"
        end if
        write(ifileid,"(a)") "    MAX_ITER 400 #Maximum number of geometry optimization"
        write(ifileid,"(a)") "    #The following thresholds of geometry convergence are the default ones"
        write(ifileid,"(a)") "    MAX_DR 3E-3 #Maximum geometry change"
        write(ifileid,"(a)") "    RMS_DR 1.5E-3 #RMS geometry change"
        write(ifileid,"(a)") "    MAX_FORCE 4.5E-4 #Maximum force"
        write(ifileid,"(a)") "    RMS_FORCE 3E-4 #RMS force"
        write(ifileid,"(a)") "  &END GEO_OPT"
    else if (itask==4) then
        write(ifileid,"(a)") "  &CELL_OPT"
        write(ifileid,"(a)") "    MAX_ITER 400 #Maximum number of geometry optimization"
        if (iprestype==1) then
            write(ifileid,"(a,1PE13.5,a)") "    EXTERNAL_PRESSURE",Piso," #External pressure for cell optimization (bar)"
        else if (iprestype==2) then
            write(ifileid,"(a,9(1PE13.5),a)") "    EXTERNAL_PRESSURE",Ptens(1,1:3),Ptens(2,1:3),Ptens(3,1:3)," #External pressure for cell optimization (bar)"
        end if
        write(ifileid,"(a)") "    CONSTRAINT "//trim(cellfix)//" #Constraint of cell length, can be: NONE, X, Y, Z, XY, XZ, YZ"
        write(ifileid,"(a)") "    KEEP_ANGLES F #If T, then cell angles will be kepted"
        write(ifileid,"(a)") "    KEEP_SYMMETRY F #If T, then cell symmetry specified by &CELL / SYMMETRY will be kepted"
        write(ifileid,"(a)") "    KEEP_SPACE_GROUP F #If T, then space group will be detected and preserved"
        write(ifileid,"(a)") "    TYPE DIRECT_CELL_OPT #Geometry and cell are optimized at the same time. Can also be GEO_OPT, MD"
        write(ifileid,"(a)") "    #The following thresholds of optimization convergence are the default ones"
        write(ifileid,"(a)") "    MAX_DR 3E-3 #Maximum geometry change"
        write(ifileid,"(a)") "    RMS_DR 1.5E-3 #RMS geometry change"
        write(ifileid,"(a)") "    MAX_FORCE 4.5E-4 #Maximum force"
        write(ifileid,"(a)") "    RMS_FORCE 3E-4 #RMS force"
        write(ifileid,"(a)") "    PRESSURE_TOLERANCE 100 #Pressure tolerance (w.r.t EXTERNAL_PRESSURE)"
        if (ioptmethod==1) then
            write(ifileid,"(a)") "    OPTIMIZER BFGS #Can also be CG (more robust for difficult cases) or LBFGS"
            write(ifileid,"(a)") "    &BFGS"
            write(ifileid,"(a)") "      TRUST_RADIUS 0.2 #Trust radius (maximum stepsize) in Angstrom"
            write(ifileid,"(a)") "#     RESTART_HESSIAN T #If read initial Hessian, uncomment this line and specify the file in the next line"
            write(ifileid,"(a)") "#     RESTART_FILE_NAME to_be_specified"
            write(ifileid,"(a)") "    &END BFGS"
        else if (ioptmethod==2) then
            write(ifileid,"(a)") "    OPTIMIZER LBFGS #Can also be CG (more robust for difficult cases) or BFGS"
        else if (ioptmethod==3) then
            write(ifileid,"(a)") "    OPTIMIZER CG #Can also be BFGS or LBFGS"
            write(ifileid,"(a)") "    &CG"
            write(ifileid,"(a)") "      &LINE_SEARCH"
            write(ifileid,"(a)") "        TYPE 2PNT #Two-point extrapolation, cheap while acceptable. Can also be FIT, GOLD"
            write(ifileid,"(a)") "      &END LINE_SEARCH"
            write(ifileid,"(a)") "    &END CG"
        end if
        write(ifileid,"(a)") "  &END CELL_OPT"
    else if (itask==6.or.itask==13) then !MD or real-time TDDFT
        write(ifileid,"(a)") "  &MD"
        if (ithermostat==0.and.ibarostat==0) then
            write(ifileid,"(a)") "    ENSEMBLE NVE"
        else if (ithermostat>1.and.ibarostat==0) then
            write(ifileid,"(a)") "    ENSEMBLE NVT"
        else if (ithermostat==0) then
            if (ibarostat==1) write(ifileid,"(a)") "    ENSEMBLE NPE_F"
            if (ibarostat==2) write(ifileid,"(a)") "    ENSEMBLE NPE_I"
        else if (ithermostat>1) then
            if (ibarostat==1) write(ifileid,"(a)") "    ENSEMBLE NPT_F"
            if (ibarostat==2) write(ifileid,"(a)") "    ENSEMBLE NPT_I"
        end if
        if (itask==6) then !MD
            write(ifileid,"(a)") "    STEPS 200 #Number of steps to run"
            write(ifileid,"(a)") "    TIMESTEP 1.0 #Step size in fs. Decrease it properly for high temperature simulation"
        else if (itask==13) then !RT-TDDFT needs very small step
            write(ifileid,"(a)") "    STEPS 20000 #Number of steps to run"
            write(ifileid,"(a)") "    TIMESTEP 0.025 #Step size in fs"
        end if
        write(ifileid,"(a)") "    TEMPERATURE 298.15 #Initial and maintained temperature (K)"
        write(ifileid,"(a)") "#   COMVEL_TOL 0 #Uncomment this can remove translation motion of center-of-mass every step"
        if (ifPBC==0) then
            write(ifileid,"(a)") "#   ANGVEL_TOL 0 #Uncomment this can remove overall rotation every step"
            write(ifileid,"(a)") "    ANGVEL_ZERO T #Eliminate overall rotation component from initial velocity"
        end if
        if (ithermostat>0) then
            write(ifileid,"(a)") "    &THERMOSTAT"
            if (ithermostat==1) then
                write(ifileid,"(a)") "      TYPE AD_LANGEVIN"
            else if (ithermostat==2) then
                write(ifileid,"(a)") "      TYPE CSVR"
                write(ifileid,"(a)") "      &CSVR"
                write(ifileid,"(a)") "        TIMECON 200 #Time constant in fs. Smaller/larger results in stronger/weaker temperature coupling"
                write(ifileid,"(a)") "      &END CSVR"
            else if (ithermostat==3) then
                write(ifileid,"(a)") "      TYPE GLE"
            else if (ithermostat==4) then
                write(ifileid,"(a)") "      TYPE NOSE"
            end if
            !if (nthermoatm<ncenter) then !Misleading
            !    write(ifileid,"(a)") "      &DEFINE_REGION"
            !    call outCP2K_LIST(ifileid,thermoatm(1:nthermoatm),nthermoatm,"        ")
            !    write(ifileid,"(a)") "      &END DEFINE_REGION"
            !end if
            write(ifileid,"(a)") "    &END THERMOSTAT"
        end if
        if (ibarostat/=0) then
            write(ifileid,"(a)") "    &BAROSTAT"
            write(ifileid,"(a)") "      PRESSURE 1.01325 #Initial and maintained pressure (bar)"
            write(ifileid,"(a)") "      TIMECON 1000 #Barostat time constant (fs)"
            if (ibarostat==1) write(ifileid,"(a)") "      VIRIAL XYZ #Relax the cell along which cartesian axes"
            write(ifileid,"(a)") "    &END BAROSTAT"
        end if
        if (itask==6) then
            write(ifileid,"(a)") "    &PRINT"
            write(ifileid,"(a)") "      &PROGRAM_RUN_INFO"
            write(ifileid,"(a)") "        &EACH"
            write(ifileid,"(a,i6,a)") "          MD",1," #Output frequency of MD information, 0 means never"
            write(ifileid,"(a)") "        &END EACH"
            write(ifileid,"(a)") "      &END PROGRAM_RUN_INFO"
            write(ifileid,"(a)") "    &END PRINT"
        end if
        write(ifileid,"(a)") "  &END MD"
    else if (itask==14) then !PINT
        write(ifileid,"(a)") "  &PINT"
        write(ifileid,"(a)") "    PROPAGATOR PIMD #Type of propagator: PIMD, CMD, RPMD"
        write(ifileid,"(a)") "    DT 0.5 #Stepsize in fs"
        write(ifileid,"(a)") "    P 8 #Number of beads"
        write(ifileid,"(a)") "    NUM_STEPS 1000 #Number of dynamics steps"
        write(ifileid,"(a)") "    TEMP 298.15 #Simulation temperature (K)"
        write(ifileid,"(a)") "    T_TOL 50.0 #Threshold for the oscillations of the temperature excedeed which the temperature is rescaled. 0 means no rescaling"
        write(ifileid,"(a)") "    TRANSFORMATION NORMAL #Coordinate transformation method: NORMAL or STAGE"
        write(ifileid,"(a)") "    HARM_INT NUMERIC #Integrator scheme for integrating the harmonic bead springs: EXACT or NUMERIC"
        write(ifileid,"(a)") "    NRESPA 1 #Number of RESPA steps for the bead for each MD step"
        write(ifileid,"(a)") "    &NOSE #Use Nose-Hoover chain thermostat"
        write(ifileid,"(a)") "      NNOS 3 #Nose-Hoover chain length"
        write(ifileid,"(a)") "    &END NOSE"
        write(ifileid,"(a)") "  &END PINT"
    end if
    if (natmcons>0) then
        write(ifileid,"(a)") "  &CONSTRAINT"
        write(ifileid,"(a)") "    &FIXED_ATOMS #Set atoms to be fixed"
        write(ifileid,"(a)") "      COMPONENTS_TO_FIX XYZ #Which fractional components will be fixed, can be X, Y, Z, XY, XZ, YZ, XYZ"
        call outCP2K_LIST(ifileid,atmcons(1:natmcons),natmcons,"      ")
        write(ifileid,"(a)") "    &END FIXED_ATOMS"
        write(ifileid,"(a)") "  &END CONSTRAINT"
    end if
    
    !https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/PRINT.html
    
    !Control output frequency of various properties
    write(ifileid,"(a)") "  &PRINT"
    if (itask==3.or.itask==4.or.itask==7) then !Optimizing minimum, TS
        write(ifileid,"(a)") "    &TRAJECTORY"
        if (iMDformat==1) write(ifileid,"(a)") "      FORMAT xyz"
        if (iMDformat==2) write(ifileid,"(a)") "      FORMAT dcd"
        if (iMDformat==3) write(ifileid,"(a)") "      FORMAT pdb"
        write(ifileid,"(a)") "    &END TRAJECTORY"
    else if (itask==6.or.itask==14) then !MD or PIMD
        write(ifileid,"(a)") "    &TRAJECTORY"
        write(ifileid,"(a)") "      &EACH"
        if (itask==6) write(ifileid,"(a,i4,a)") "        MD",nMDsavefreq," #Output frequency of coordinates, 0 means never"
        if (itask==14) write(ifileid,"(a,i4,a)") "        PINT",nMDsavefreq," #Output frequency of coordinates, 0 means never"
        write(ifileid,"(a)") "      &END EACH"
        if (iMDformat==1) write(ifileid,"(a)") "      FORMAT xyz"
        if (iMDformat==2) write(ifileid,"(a)") "      FORMAT dcd"
        if (iMDformat==3) write(ifileid,"(a)") "      FORMAT pdb"
        write(ifileid,"(a)") "    &END TRAJECTORY"
        write(ifileid,"(a)") "    &VELOCITIES"
        write(ifileid,"(a)") "      &EACH"
        if (itask==6) write(ifileid,"(a,i6,a)") "        MD",0," #Output frequency of velocities, 0 means never"
        if (itask==14) write(ifileid,"(a,i6,a)") "        PINT",0," #Output frequency of velocities, 0 means never"
        write(ifileid,"(a)") "      &END EACH"
        write(ifileid,"(a)") "    &END VELOCITIES"
        write(ifileid,"(a)") "    &FORCES"
        write(ifileid,"(a)") "      &EACH"
        if (itask==6) write(ifileid,"(a,i6,a)") "        MD",0," #Output frequency of forces, 0 means never"
        if (itask==14) write(ifileid,"(a,i6,a)") "        PINT",0," #Output frequency of forces, 0 means never"
        write(ifileid,"(a)") "      &END EACH"
        write(ifileid,"(a)") "    &END FORCES"
    end if
    
    write(ifileid,"(a)") "    &RESTART"
    write(ifileid,"(a)") "      BACKUP_COPIES 0 #Maximum number of backing up restart file, 0 means never" !Do not generate annoying .restart.bak file
    !For other tasks, by default, restart file is updated every step. Only for MD it is default to 20, I explicitly provide option to change it
    if (itask==6.or.itask==14) then
        write(ifileid,"(a)") "      &EACH"
        if (itask==6) write(ifileid,"(a)") "        MD  1 #Frequency of updating last restart file, 0 means never"
        if (itask==14) write(ifileid,"(a)") "        PINT  1 #Frequency of updating last restart file, 0 means never"
        write(ifileid,"(a)") "      &END EACH"
    end if
    write(ifileid,"(a)") "    &END RESTART"
    
    !Control how to generate history .restart files
    !For GEO_OPT and MD, default is 500. For other tasks, default is every step
    !Because it is useless, so I completely suppress it
    if (itask==4.or.itask==6.or.itask==14) then !Cell opt, MD, PINT
        write(ifileid,"(a)") "    &RESTART_HISTORY OFF "
        write(ifileid,"(a)") "    &END RESTART_HISTORY"
    end if
    !write(ifileid,"(a)") "    &RESTART_HISTORY"
    !write(ifileid,"(a)") "      &EACH"
    !write(ifileid,"(a)") "        CELL_OPT 0 #How often a history .restart file is generated, 0 means never"
    !write(ifileid,"(a)") "      &END EACH"
    !write(ifileid,"(a)") "    &END RESTART_HISTORY"
    write(ifileid,"(a)") "  &END PRINT"
    write(ifileid,"(a)") "&END MOTION"
end if
if (itask==5) then
    write(ifileid,"(a)") "&VIBRATIONAL_ANALYSIS"
    write(ifileid,"(a)") "  DX 0.01 #Step size of finite difference. This is default (Bohr)"
    write(ifileid,"(a)") "  NPROC_REP 1 #Number of processors to be used per replica. This is default"
    write(ifileid,"(a)") "  TC_PRESSURE 101325 #1 atm. Pressure for calculate thermodynamic data (Pa)"
    write(ifileid,"(a)") "  TC_TEMPERATURE 298.15 #Temperature for calculate thermodynamic data (K)"
    write(ifileid,"(a)") "  THERMOCHEMISTRY #Print thermochemistry information (only valid for molecule in gas!)"
    write(ifileid,"(a)") "  INTENSITIES T #Calculate IR intensities"
    if (ifPBC==0) then
        write(ifileid,"(a)") "  FULLY_PERIODIC F #This is default. If T, avoiding to project out rotation component from Hessian matrix"
    else
        write(ifileid,"(a)") "  FULLY_PERIODIC T #Avoiding to project out rotation component from Hessian matrix"
    end if
    if (ioutvibmol==1) then
        write(ifileid,"(a)") "  &PRINT"
        write(ifileid,"(a)") "    &MOLDEN_VIB #Output .mol (Molden file) for visualization vibrational modes"
        write(ifileid,"(a)") "    &END MOLDEN_VIB"
        write(ifileid,"(a)") "  &END PRINT"
    end if
    write(ifileid,"(a)") "&END VIBRATIONAL_ANALYSIS"
end if

close(ifileid)

write(*,"(a)") " CP2K input file has been exported to "//trim(outname)
end subroutine

!----- Test if a given integer array contains contiguous numbers (index order is unimportant)
!If yes, ifconti=1, else =0
subroutine testidx_contiguous(array,narray,ifconti)
implicit real*8 (a-h,o-z)
integer narray,array(narray)
ifconti=0
do i=minval(array),maxval(array)
    if (all(array/=i)) return
end do
ifconti=1
end subroutine

!----- Output integer arrays in CP2K "LIST" convention
!list(1:nlist) are indices to be exported to ifileid. spacestr is string containing proper number of spaces to be outputted in front of data
subroutine outCP2K_LIST(ifileid,list,nlist,spacestr)
character(len=*) spacestr
character c80tmp*80
integer ifileid,nlist,ifconti
integer list(nlist)
call testidx_contiguous(list,nlist,ifconti)
if (ifconti==1) then
    write(c80tmp,*) maxval(list(1:nlist))
    write(ifileid,"(a,i8,'..',a)") spacestr//"LIST ",minval(list(1:nlist)),adjustl(c80tmp)
else
    write(ifileid,"(a)",advance='no') spacestr//"LIST "
    i=0
    do while(.true.) !Change line ended with \ when outputting every 12 terms
        write(ifileid,"(a)",advance='no') spacestr
        if (i*12<nlist-12) then !Not the last line, the remaining terms is more than 12
            write(ifileid,"(12i6,' \')") list(12*i+1:12*i+12)
            i=i+1
        else !The last line
            write(ifileid,"(12i6)") list(12*i+1:nlist)
            exit
        end if
    end do
end if
end subroutine



!!---------- Automatically set proper vacuum sizes
subroutine determine_vacuumsize(iPSOLVER,vacsizex,vacsizey,vacsizez,icentering)
use defvar
integer iPSOLVER,icentering
real*8 vacsizex,vacsizey,vacsizez
if (iPSOLVER==1) then !Usually adequate for PERIODIC
    vacsizex=5/b2a
    vacsizey=5/b2a
    vacsizez=5/b2a
else if (iPSOLVER==2) then !ANALYTIC converges quite slow
    vacsizex=10/b2a
    vacsizey=10/b2a
    vacsizez=10/b2a
else if (iPSOLVER==3) then !MT needs vaccum size in each side is >= half of system
    vacsizex=((maxval(a%x)-minval(a%x))*b2a+2*4D0)/2/b2a !4 A is extension distance for electron tail in each side
    vacsizey=((maxval(a%y)-minval(a%y))*b2a+2*4D0)/2/b2a
    vacsizez=((maxval(a%z)-minval(a%z))*b2a+2*4D0)/2/b2a
else if (iPSOLVER==4) then
    icentering=1
    vacsizex=3.5/b2a !WAVELET converges quite fast, 3.5 A is adequate for any case
    vacsizey=3.5/b2a
    vacsizez=3.5/b2a
end if
end subroutine



!!--------- Convert band structure file of CP2K (.bs) to multiple column file so that band map can be directly plotted by Origin or gnuplot
subroutine CP2K_BS
use defvar
use util
implicit real*8 (a-h,o-z)
character c200tmp*200,c10tmp*10
integer,parameter :: nSPmax=1000,nkpmax=10000
character SPlabel(nSPmax)*10 !Label of special points
real*8 SPvec(3,nSPmax) !XYZ of special points
integer*2 SPpath(nSPmax) !The path number that this special point belongs to. Each path consists of connected special points
integer*2 SPkp(nSPmax) !Index of k-point list of each special point
real*8 tmpvec(3),kpvec(3,nkpmax)
integer kpplot(nkpmax) !kpplot(i) is k-point index for plotting band map of real k-point i
E_VBT=-99999 !Valence band top
E_CBB=99999 !Conduction band bottom

write(*,*) "Input value for shifting energy levels, e.g. 3.42"
read(*,*) eshift

open(10,file=filename,status="old")

iopsh=0 !Closed-shell
call loclabel(10,"Spin 2",iopsh)
if (iopsh==0) then
    open(11,file="BS_occ.txt",status="replace")
    open(12,file="BS_vir.txt",status="replace")
else if (iopsh==1) then
    open(11,file="BS_A_occ.txt",status="replace")
    open(12,file="BS_A_vir.txt",status="replace")
    open(13,file="BS_B_occ.txt",status="replace")
    open(14,file="BS_B_vir.txt",status="replace")
end if

rewind(10)
nSP=0 !Number of currently loaded special points
nkp=0 !Current k-point index
nkpplot=0
ipath=1 !Current path index
xkp_last=-999
ykp_last=-999
zkp_last=-999
do while(.true.)
    !Loop all sets
    read(10,"(a)",iostat=ierror) c200tmp
    if (ierror/=0.or.c200tmp==" ") exit
    if (index(c200tmp,"# Set")/=0) then
        read(c200tmp,"(5x,i2,1x,i2,16x,i3,10x,i3)") iset,nSPthis,nkpthis,nlevelthis
        write(*,"(a,i3,'...')") " Loading set",iset
    end if
    
    !Read position and label of special points in current set
    do iSP=1,nSPthis
        read(10,"(a)") c200tmp
        read(c200tmp(24:),*) tmpvec(:),c10tmp
        if (iSP==1.and.nSP/=0) then !A new set but not the first set
            if (c10tmp/=SPlabel(nSP)) ipath=ipath+1 !The first special point in this set has different label to the last one, so begins a new path
        end if
        iadd=0
        if (nSP==0) then
            iadd=1
        else if (c10tmp/=SPlabel(nSP)) then
            iadd=1
        end if
        if (iadd==1) then !Add this special point to list
            nSP=nSP+1
            SPvec(:,nSP)=tmpvec(:)
            SPlabel(nSP)=c10tmp
            SPpath(nSP)=ipath
        end if
    end do
    
    !Loop k-points between neighbouring special points in current set
    do ikp=1,nkpthis
        !Alpha part
        read(10,"(a)") c200tmp
        read(10,*)
        read(c200tmp(25:),*) xkp,ykp,zkp !XYZ of k-point
        !write(*,"(2i5,6f10.6)") iset,ikp,xkp,ykp,zkp,xkp_last,ykp_last,zkp_last
        
        !The first point of this set is different to final point of the last set, make they share the same k-point index
        if (ikp==1.and.iset>1.and.(xkp/=xkp_last.or.ykp/=ykp_last.or.zkp/=zkp_last)) nkpplot=nkpplot-1
        
        if (nkp>1.and.xkp==xkp_last.and.ykp==ykp_last.and.zkp==zkp_last) then !This k-point is identical to the last one, so skip this one
            call skiplines(10,nlevelthis)
            if (iopsh==1) call skiplines(10,nlevelthis+2)
            cycle
        else !This k-point will be actually loaded
            nkp=nkp+1
            kpvec(1,nkp)=xkp
            kpvec(2,nkp)=ykp
            kpvec(3,nkp)=zkp
            nkpplot=nkpplot+1
            kpplot(nkp)=nkpplot
            xkp_last=xkp
            ykp_last=ykp
            zkp_last=zkp
        end if
        !write(11,"(i4,3f8.4)",advance="no") nkpplot,xkp,ykp,zkp !For checking xyz of k-point
        !write(12,"(i4,3f8.4)",advance="no") nkpplot,xkp,ykp,zkp
        write(11,"(i4)",advance="no") nkpplot
        write(12,"(i4)",advance="no") nkpplot
        !Loop energy levels in current k-point
        do ilevel=1,nlevelthis
            read(10,*) idx,ene,occ
            if (occ/=0) then
                write(11,"(f14.8)",advance="no") ene+eshift
                if (ene>E_VBT) E_VBT=ene
            else
                write(12,"(f14.8)",advance="no") ene+eshift
                if (ene<E_CBB) E_CBB=ene
            end if
        end do
        write(11,*)
        write(12,*)
        
        !Beta part
        if (iopsh==1) then
            read(10,*);read(10,*)
            write(13,"(i4)",advance="no") nkpplot
            write(14,"(i4)",advance="no") nkpplot
            do ilevel=1,nlevelthis
                read(10,*) idx,ene,occ
                if (occ/=0) then
                    write(13,"(f14.8)",advance="no") ene+eshift
                    if (ene>E_VBT) E_VBT=ene
                else
                    write(14,"(f14.8)",advance="no") ene+eshift
                    if (ene<E_CBB) E_CBB=ene
                end if
            end do
            write(13,*)
            write(14,*)
        end if
    end do
    nSP_last=nSP
end do

close(10)
close(11)
close(12)
if (iopsh==0) then
    write(*,"(/,a)") " Occupied and virtual levels for all k-points involved in band structure map have been written to BS_occ.txt and BS_vir.txt, respectively"
else if (iopsh==1) then
    close(13)
    close(14)
    write(*,"(/,a)") " Alpha occupied and virtual levels for all k-points involved in band structure map have been written to BS_A_occ.txt and BS_A_vir.txt, respectively"
    write(*,"(a)") " Beta occupied and virtual levels for all k-points involved in band structure map have been written to BS_B_occ.txt and BS_B_vir.txt, respectively"
end if

!write(*,*) "List of k-points:"
!do ikp=1,nkp
!    write(*,"(2i4,3f12.6)") ikp,kpplot(ikp),kpvec(:,ikp)
!end do

!Determine correspondence between special points and (actual) k-point index
SPkp(:)=0
do iSP=1,nSP
    if (iSP==1) then
        ibeg=1
    else
        ibeg=SPkp(iSP-1)+1
    end if
    do ikp=ibeg,nkp
        if (SPvec(1,iSP)==kpvec(1,ikp).and.SPvec(2,iSP)==kpvec(2,ikp).and.SPvec(3,iSP)==kpvec(3,ikp)) then
            SPkp(iSP)=ikp
            exit
        end if
    end do
end do
write(*,*)
write(*,"(' Number of total levels:',i7)") nlevelthis
write(*,*)
write(*,"(' Valence band top:      ',f14.8,' eV (unshifted)')") E_VBT
write(*,"(' Conduction band bottom:',f14.8,' eV (unshifted)')") E_CBB
write(*,"(' Band gap:              ',f14.8,' eV')") E_CBB-E_VBT
if (E_CBB<=E_VBT) then
    write(*,*) "Note: This is a metal, band gap in fact is zero"
end if
write(*,*)
write(*,*) "Special points for plotting band structure:"
write(*,*) "SP#   kp#   Label                Coordinates          Path"
do iSP=1,nSP
    if (iSP>1) then
        if (SPpath(iSP)/=SPpath(iSP-1)) write(*,*)
    end if
    write(*,"(i4,i5,4x,a,3f10.6,i5)") iSP,kpplot(SPkp(iSP)),SPlabel(iSP),SPvec(:,iSP),SPpath(iSP)
end do

end subroutine