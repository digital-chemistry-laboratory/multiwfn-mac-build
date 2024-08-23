!!---------- Extended Transition State - Natural Orbitals for Chemical Valence (ETS-NOCV)
!Terminology:
!Promolecular (pro) orbital/density/state: The total wavefunction directly combined from fragment wavefunctions
!Frozen state (frz) orbital/density: Wavefunction of antisymmetric product of occupied fragment orbitals &
!                                    via Lowdin orthogonalization (virtual orbitals are the same as original ones)
!If iopsh=1, means any of fragment or the whole system is unrestricted open-shell, thus alpha and beta are needed treated separately
!    About details of combination of fragment wavefunction:
!The CObasa_pro is a diagonal block matrix directly combined from CObasa of fragment, for example &
!(1:nbasis_A,1:nbasis_A) is CObasa of fragment A, (nbasis_A+1:nbasis_A+nbasis_B,nbasis_A+1:nbasis_A+nbasis_B) is CObasa of fragment B, etc.
!(1:nbasis_A,1:naelec_A) is coefficients of occupied orbitals of fragment A, (nbasis_A+1:nbasis_A+nbasis_B,nbasis_A+1:nbasis_A+nbelec) is that of fragment B, etc.
!The CObasa_frz is mix between (:,1:naelec_A) with (:,nbasis_A+1:nbasis_A+nbelec) using Lowdin orthogonalization
subroutine ETS_NOCV
use defvar
use util
use GUI
use functions
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,c2000tmp*200
integer :: nfrag=0 !Number of fragments
integer :: igridstat=0 !State of defining grid for visualization. 0/1 = Undefined/defined
character(len=200),allocatable :: fragfilename(:) !Filename of each fragment
real*8,allocatable :: tmpmat(:,:) !Temporary use
real*8 NOCVorb(nbasis,nbasis),NOCVeig(nbasis) !For backing up NOCV orbitals
real*8,allocatable :: NOCVorbB(:,:),NOCVeigB(:) !Same as above, for beta spin. Used only for open-shell case
real*8,allocatable :: NOCVene(:) !NOCV orbital energies
!Density matrices
real*8 Pfrz(nbasis,nbasis) !Frozen state density matrix of all electrons
real*8,allocatable :: PfrzA(:,:),PfrzB(:,:) !Same as above, for alpha and beta electrons. Used only for open-shell case
real*8 Pdiff(nbasis,nbasis) !Difference density matrix (real w.r.t. frozen state) of all electrons
real*8,allocatable :: PdiffA(:,:),PdiffB(:,:) !Same as above, for alpha and beta electrons. Used only for open-shell case
!Coefficient matrics
real*8 CObasapro(nbasis,nbasis) !Coefficient matrix of promolecular orbitals
real*8,allocatable :: CObasbpro(:,:) !Same as above, for beta spin. Used only for open-shell case
real*8 CObasafrz(nbasis,nbasis) !Coefficient matrix of frozen state orbitals
real*8,allocatable :: CObasbfrz(:,:) !Same as above, for beta spin. Used only for open-shell case
!Other arrays
real*8 occfrag(nbasis),occfragB(nbasis) !Occupation number of fragment orbitals, for all/alpha and beta spins
integer,allocatable :: ibasfrag(:) !ibasfrag(ifrag) is the starting basis function index of fragment i
integer,allocatable :: naelec_frag(:),nbelec_frag(:) !Number of alpha and beta electrons of each fragment 
integer,allocatable :: pairidx1(:),pairidx2(:) !Index of the two NOCV orbitals corresponding to a pair. If number of orbitals is odd, one orbital will not paired
integer,allocatable :: sellist(:)
real*8,allocatable :: rho_complex(:,:,:),rho_frz(:,:,:),rho_pro(:,:,:)
real*8 :: printthres=0.001D0
real*8,allocatable :: wfnval(:)

!Check input file
if (.not.allocated(CObasa)) then
    write(*,*) "Error: Your input file must contain basis function information!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if
if (wfntype==2) then
    write(*,*) "Error: Restricted open-shell wavefunction is not supported by this function!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
else if (wfntype==3.or.wfntype==4) then
    write(*,*) "Error: Only single-determinant wavefunction is supported by this function!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

write(*,*)
write(*,*) "How many fragments do you want to define for NOCV analysis?  e.g. 2"
do while(.true.)
    !nfrag=2 !For debug !!!
	read(*,*) nfrag
	if (nfrag>=2) exit
	write(*,*) "Error: The number of fragments must be >= 2, input again"
end do
allocate(ibasfrag(nfrag),fragfilename(nfrag),naelec_frag(nfrag),nbelec_frag(nfrag))

nbasis_org=nbasis
if (allocated(CObasa_org)) deallocate(CObasa_org)
allocate(CObasa_org(nbasis,nbasis))
CObasa_org=CObasa
nelec_org=nint(nelec)
naelec_org=nint(naelec)
nbelec_org=nint(nbelec)
nbasis_comb=0 !Information of combined fragments
nelec_comb=0
naelec_comb=0
nbelec_comb=0
ncenter_comb=0
iopsh=0
if (wfntype==1) iopsh=1

!Load each fragment
!There are two stages:
!(1) Load atom information and number of basis functions for checking sanity, and determine "iopsh"
!(2) Load orbital coefficients, construct promolecular density matrix
!fragfilename(1)="examples\ETS-NOCV\COBH3\CO.fch" !For debug !!!
!fragfilename(2)="examples\ETS-NOCV\COBH3\BH3.fch" !For debug !!!
do ifrag=1,nfrag
	write(*,"(/,a,i4)") " Input path of wavefunction file of fragment",ifrag
    write(*,*) "For example, D:\ll\A-Rise.mwfn"
    do while(.true.)
		read(*,"(a)") fragfilename(ifrag)
		inquire(file=fragfilename(ifrag),exist=alive)
		if (alive) exit
		write(*,*) "Cannot find the file, input again!"
    end do
	call dealloall(0)
	call readinfile(fragfilename(ifrag),1)
    !Check consistency of this fragment with respect to complex
    if (ncenter_comb+ncenter>ncenter_org) then
		write(*,"(a,i5,a,i5,a)") " Error: Total number of atoms of already loaded fragments (",ncenter_comb+ncenter,&
        ") is larger than that of the whole system (",ncenter_org,")!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
	    call dealloall(0)
	    call readinfile(firstfilename,1)
		return
    end if
    if (nbasis_comb+nbasis>nbasis_org) then
		write(*,"(a,i5,a,i5,a)") " Error: Total number of basis functions of already loaded fragments (",nbasis_comb+nbasis,&
        ") is larger than that of the whole system (",nbasis_org,")!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
	    call dealloall(0)
	    call readinfile(firstfilename,1)
		return
    end if
    do iatm=1,ncenter
		ncenter_comb=ncenter_comb+1
        xyzdev=dsqrt( (a(iatm)%x-a_org(ncenter_comb)%x)**2 + (a(iatm)%y-a_org(ncenter_comb)%y)**2 + (a(iatm)%z-a_org(ncenter_comb)%z)**2 )
        if (xyzdev>0.01D0) then
			write(*,"(a,i5,a,i5,a)") " Error: Atom",iatm," of this fragment has evidently different coordinate to the corresponding atom (",ncenter_comb," ) of the whole system!"
            write(*,"(' X,Y,Z of atom',i5,' in this fragment:',3f10.5,' Angstrom')") iatm,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
            write(*,"(' X,Y,Z of atom',i5,' in whole system: ',3f10.5,' Angstrom')") ncenter_comb,a_org(ncenter_comb)%x*b2a,a_org(ncenter_comb)%y*b2a,a_org(ncenter_comb)%z*b2a
			write(*,*) "Press ENTER button to return"
            read(*,*)
	        call dealloall(0)
	        call readinfile(firstfilename,1)
			return
        end if
    end do
    ibasfrag(ifrag)=nbasis_comb+1 !Starting basis function index of this fragment
    nbasis_comb=nbasis_comb+nbasis
    if (wfntype==1) then
        iopsh=1
    else if (wfntype==2) then
        write(*,*) "Error: Restricted open-shell wavefunction is not supported by this function!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
	    call dealloall(0)
	    call readinfile(firstfilename,1)
        return
    end if
end do

!Loading wavefunction information
CObasapro=0
if (iopsh==1) then
    allocate(CObasbpro(nbasis_org,nbasis_org))
    CObasbpro=0
end if
write(*,*)
do ifrag=1,nfrag
	call dealloall(0)
	call readinfile(fragfilename(ifrag),1)
    write(*,"(' Loading wavefunction information from ',a)") trim(fragfilename(ifrag))
    ibeg=ibasfrag(ifrag)
    !Load coefficient matrix of fragment
    CObasapro(ibeg:ibeg+nbasis-1,ibeg:ibeg+nbasis-1)=CObasa(:,:)
    occfrag(ibeg:ibeg+nbasis-1)=MOocc(1:nbasis)
    naelec_frag(ifrag)=nint(naelec)
    nbelec_frag(ifrag)=nint(nbelec)
    if (iopsh==1) then
        if (wfntype==0) then !This is a closed-shell fragment, let alpha share same information as beta
            occfrag(ibeg:ibeg+nbasis-1)=MOocc(1:nbasis)/2
            occfragB(ibeg:ibeg+nbasis-1)=MOocc(1:nbasis)/2
            CObasbpro(ibeg:ibeg+nbasis-1,ibeg:ibeg+nbasis-1)=CObasa(:,:)
        else if (wfntype==1) then
            occfragB(ibeg:ibeg+nbasis-1)=MOocc(nbasis+1:nmo)
            CObasbpro(ibeg:ibeg+nbasis-1,ibeg:ibeg+nbasis-1)=CObasb(:,:)
            !Ask how to flip spin
            write(*,"(' Fragment',i4,', number of alpha and beta electrons:',2i5)") ifrag,nint(naelec),nint(nbelec)
            write(*,*) "Do you want to flip spin (exchange information of alpha and beta spins)? (y/n)"
            read(*,*) c80tmp
            if (c80tmp=="y".or.c80tmp=="Y") then
                CObasapro(ibeg:ibeg+nbasis-1,ibeg:ibeg+nbasis-1)=CObasb(:,:)
                CObasbpro(ibeg:ibeg+nbasis-1,ibeg:ibeg+nbasis-1)=CObasa(:,:)
                occfrag(ibeg:ibeg+nbasis-1)=MOocc(nbasis+1:nmo)
                occfragB(ibeg:ibeg+nbasis-1)=MOocc(1:nbasis)
                naelec_frag(ifrag)=nint(nbelec)
                nbelec_frag(ifrag)=nint(naelec)
            end if
        end if
    end if
end do

write(*,"(/,' Alpha and beta electrons of all fragments:',2i7)") sum(naelec_frag),sum(nbelec_frag)

!Check consistency of number of electrons between complex and all fragments
if (sum(naelec_frag)/=naelec_org) then
    write(*,"(a,i5,a,i5,a)") " Total number of alpha electrons of all fragments (",sum(naelec_frag),") is different to that of the whole system (",naelec_org,")!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
	call dealloall(0)
	call readinfile(firstfilename,1)
    return
else if (sum(nbelec_frag)/=nbelec_org) then
    write(*,"(a,i5,a,i5,a)") " Total number of beta electrons of all fragments (",sum(nbelec_frag),") is different to that of the whole system (",nbelec_org,")!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
	call dealloall(0)
	call readinfile(firstfilename,1)
    return
end if

!Reloading the complex wavefunction
write(*,"(/,a,a)") " Reloading ",trim(firstfilename)
call dealloall(0)
call readinfile(firstfilename,1) !Recover to the first file

!Backup some arrays, and in the case of open-shell NOCV, split closed-shell complex wavefunction to alpha and beta
if (iopsh==1) then
    if (wfntype==0) then !Complex is closed-shell but open-shell analysis is needed, so transform it to unrestricted case
        nmo=2*nbasis
        deallocate(MOene_org,MOocc_org) !Backup MO energy/occupancy as two spins, used later when visualizing complex orbitals
        allocate(MOene_org(nmo),MOocc_org(nmo))
        MOene_org(1:nbasis)=MOene
        MOene_org(nbasis+1:nmo)=MOene
        MOocc_org(1:nbasis)=MOocc/2
        MOocc_org(nbasis+1:nmo)=MOocc/2
        if (allocated(CObasb_org)) deallocate(CObasb_org)
        allocate(CObasb(nbasis,nbasis),CObasb_org(nbasis,nbasis)) !Extent beta part as alpha part, and backup
        CObasb_org=CObasa
        deallocate(CO,MOocc,MOene,MOtype)
        allocate(CO(nmo,nprims),MOocc(nmo),MOene(nmo),MOtype(nmo)) !Make them have proper size
        MOtype(1:nbasis)=1
        MOtype(nbasis+1:nmo)=2
        allocate(Palpha(nbasis,nbasis),Pbeta(nbasis,nbasis))
        Palpha=Ptot/2
        Pbeta=Palpha
    else if (wfntype==1) then !Need to backup complex CObasb, may be used later
        allocate(CObasb_org(nbasis,nbasis))
        CObasb_org=CObasb
    end if
end if

!Perform Lowdin orthogonalization between occupied orbitals to generate frozen state wavefunction
if (iopsh==0) then
    ndimB=0
    allocate(CObasbpro(0,0),CObasbfrz(0,0),PfrzA(0,0),PfrzB(0,0))
else if (iopsh==1) then
    ndimB=nbasis
    allocate(CObasbfrz(nbasis,nbasis),PfrzA(nbasis,nbasis),PfrzB(nbasis,nbasis))
end if
call occorb_Lowdinorth(1,iopsh,ndimB,occfrag,occfragB,CObasapro,CObasbpro,CObasafrz,CObasbfrz,Pfrz,PfrzA,PfrzB)

!Generate difference density matrix
write(*,*)
write(*,*) "Calculating difference density matrix (delta-P_rho) ..."
if (iopsh==0) then
    Pdiff=Ptot-Pfrz
else if (iopsh==1) then
	PdiffA=Palpha-PfrzA
	PdiffB=Pbeta-PfrzB
end if

!Generate NOCV orbitals
!"gennatorb" is originally used to generate NOs, I borrow it to generate NOCVs by diagonalizing Pdiff
!  ETS-NOCV original paper represents delta-P in the basis of orthogonalized fragment orbitals, this is cumbersome. I just diagonalize &
!  delta-P in the basis of original basis functions; however, of course, the equation is (delta-P)C = SCv, so Lowdin orthogonalization &
!  is applied to basis functions first to convert it to (delta-P)C = Cv, this is what "gennatorb" does
write(*,*)
write(*,*) "Generating NOCV orbitals and eigenvalues ..."
if (iopsh==0) then !Diagonalizing total density matrix to generate NOCV
    Ptot=Pdiff
    call gennatorb(1,0) !This subroutine deals with Ptot (which is difference density matrix currently), generate CObasa (NOCVs currently) and MOocc (eigenvalues currently)
    call CObas2CO(1) !Convert CObasa to CO
    NOCVorb=CObasa !Backing up
    NOCVeig=MOocc
else if (iopsh==1) then !Diagonalizing alpha and beta matrices to generate NOCV of alpha and beta spins
    Palpha=PdiffA
    Pbeta=PdiffB
    call gennatorb(2,0) !Deal with Palpha and Pbeta, generate CObasa, CObasb, MOocc
    call CObas2CO(3) !Convert CObasa and CObasb to CO
    NOCVorb=CObasa !Backing up
    NOCVeig=MOocc
    allocate(NOCVorbB(nbasis,nbasis),NOCVeigB(nbasis))
    NOCVorbB=CObasb
    NOCVeigB=MOocc(nbasis+1:nmo)
end if
allocate(NOCVene(nmo))
NOCVene=0 !NOCV orbital energies, not evaluated by default
MOene=0 !Currently meaningless
write(*,"(/,a)") " Now orbitals in memory correspond to NOCV eigenvectors, and occupation numbers correspond to NOCV eigenvalues"
write(*,"(a)") " NOCV orbital energies are not calculated. To obtain them, you should select option -1 or -2 to load or generate Fock/KS matrix"

!Construct pair list of NOCV orbitals
!When number of basis functions is odd, then pairidx1 and pairidx2 of the last pair correspond to the same NOCV orbital
if (iopsh==0) then
    npair=ceiling(nbasis/2D0)
    allocate(pairidx1(npair),pairidx2(npair))
    do ipair=1,npair
        pairidx1(ipair)=ipair
        pairidx2(ipair)=nbasis+1-ipair
    end do
else if (iopsh==1) then
    npair=ceiling(nbasis/2D0)*2
    allocate(pairidx1(npair),pairidx2(npair))
    do ipair=1,npair/2
        !Alpha part
        pairidx1(ipair)=ipair
        pairidx2(ipair)=nbasis+1-ipair
        !Beta part
        pairidx1(ipair+npair/2)=nbasis+ipair
        pairidx2(ipair+npair/2)=2*nbasis+1-ipair
    end do
end if

!Show result
call showNOCV(6,pairidx1,pairidx2,npair,iopsh,printthres)

!Set default visualization setting to suitable status for showing NOCV related isosurface
drawisosurgui_SWGSTP_old=drawisosurgui_SWGSTP
drawisosurgui_lowlim_old=drawisosurgui_lowlim
drawisosurgui_highlim_old=drawisosurgui_highlim
sur_value_iso=0.005D0
sur_value_orb=0.05D0
aug3D=3D0 !NOCV related density mostly in internal side, so this extension distance is enough
aug3D_main0=4D0 !Used in showing NOCV orbitals
drawisosurgui_SWGSTP=0.001D0 !Visualize density difference need to use very small isovalue, so stepsize and range of drag bar should also be very small
drawisosurgui_lowlim=0D0
drawisosurgui_highlim=0.2D0

do while(.true.)
    write(*,*)
    write(*,*) "           ---------------- Post-processing menu ------------------"
    write(*,*) "-10 Return"
    write(*,*) "-6 Manually define correspondence between NOCV pairs and orbitals"
    if (igridstat==0) then
        write(*,*) "-5 Set grid for calculation of various densities, current: Undefined"
    else if (igridstat==1) then
        write(*,"(a,i10,a)") " -5 Set grid for calculation of various densities, current:",nx*ny*nz," points"
    end if
    write(*,*) "-4 Export NOCV information to a plain text file"
    write(*,"(a,f8.5)") " -3 Set printing threshold of NOCV eigenvalues, current:",printthres
    if (.not.allocated(FmatA)) then
        write(*,*) "-2 Generate Fock/KS matrix and evaluate NOCV orbital energies"
        write(*,*) "-1 Load Fock/KS matrix and evaluate NOCV orbital energies"
    else
        write(*,*) "-2 Generate Fock/KS matrix and re-evaluate NOCV orbital energies"
        write(*,*) "-1 Load Fock/KS matrix and re-evaluate NOCV orbital energies"
    end if
    write(*,*) " 0 Print NOCV information"
    write(*,*) " 1 Show isosurface of NOCV orbitals"
    write(*,*) " 2 Show isosurface of NOCV pair density"
    write(*,*) " 3 Show isosurface of Pauli deformation density"
    write(*,*) " 4 Show isosurface of orbital deformation density"
    write(*,*) " 5 Show isosurface of total deformation density"
    write(*,*) " 6 Export cube file of a NOCV orbital"
    write(*,*) " 7 Export cube file of NOCV pair density"
    write(*,*) " 8 Export cube file of Pauli deformation density"
    write(*,*) " 9 Export cube file of orbital deformation density"
    write(*,*) "10 Export cube file of total deformation density"
    write(*,*) "11 Visualize promolecular orbitals"
    write(*,*) "12 Visualize frozen state orbitals"
    write(*,*) "13 Visualize actual complex orbitals"
    write(*,*) "14 Calculate composition of NOCV orbitals and pairs"
    read(*,*) isel
    
    if (isel==-10) then !Return
        drawisosurgui_SWGSTP=drawisosurgui_SWGSTP_old
        drawisosurgui_lowlim=drawisosurgui_lowlim_old
        drawisosurgui_highlim=drawisosurgui_highlim_old
        write(*,"(/,a)") " Note that current orbitals in memory correspond to NOCV eigenvectors, &
        &and occupation numbers correspond to NOCV eigenvalues"
        exit
    
    else if (isel==-6) then !Manually define correspondence between NOCV pairs and orbitals
        do while(.true.)
            write(*,*) "Current NOCV pairs and orbitals are shown below:"
            call showNOCV(6,pairidx1,pairidx2,npair,iopsh,printthres)
            write(*,*)
            write(*,*) "Define which NOCV pair? e.g. 2"
            if (iopsh==1) write(*,"(a)") "If no suffix, alpha pair will be selected; if suffix is added, e.g. 3 b, beta pair will be selected"
            write(*,*) "Input ""q"" can return"
            read(*,"(a)") c80tmp
            if (index(c80tmp,"q")/=0) exit
            read(c80tmp,*) ipair
            write(*,*) "Input indices of two NOCV orbitals corresponding to this NOCV pair, e.g. 20,37"
            read(*,*) pairidx1(ipair),pairidx2(ipair)
            write(*,*) "Done!"
        end do
        
    else if (isel==-5.or.(isel>=2.and.isel<=10.and.igridstat==0)) then
        call setgrid(0,igridsel) !Haven't set grid, ask user to set now
        igridstat=1
        write(*,*)
        
    else if (isel==-4) then
        write(*,*) "Output NOCV information to which file? e.g. E:\K-ON\Azusa_Nakano.txt"
        write(*,"(a)") " If press ENTER button directly, NOCV information will be outputted to NOCV.txt in current folder"
        read(*,"(a)") c200tmp
        if (c200tmp=="") c200tmp="NOCV.txt"
        open(10,file=trim(c200tmp),status="replace")
        call showNOCV(10,pairidx1,pairidx2,npair,iopsh,printthres)
        close(10)
        write(*,*) "Outputting finished!"
        
    else if (isel==-3) then
        write(*,*) "Input printing threshold of NOCV eigenvalues, e.g. 0.05"
        read(*,*) printthres
        
    else if (isel==-1.or.isel==-2) then !Load or Generate Fock/KS matrix and evaluate NOCV orbital energies
        if (isel==-1) then !Load Fock
            call loadFockfile(istatus)
        else if (isel==-2) then !Generate Fock
            CObasa=CObasa_org
            if (iopsh==1) CObasb=CObasb_org
            MOene=MOene_org
            call MOene2Fmat(istatus)
            CObasa=NOCVorb !Restore to NOCV information
            if (iopsh==1) CObasb=NOCVorbB
        end if
        write(*,*)
        if (istatus==0) then
            write(*,*) "Calculating energies of NOCV orbitals ..."
            !Simply transforming Fock matrix from AO basis to NOCV basis, the diagonal terms are NOCV orbital energies
            allocate(tmpmat(nbasis,nbasis))
            !tmpmat=matmul(transpose(CObasa),matmul(FmatA,CObasa)) !Slower code
            tmpmat=matmul_blas(FmatA,CObasa,nbasis,nbasis)
            tmpmat=matmul_blas(CObasa,tmpmat,nbasis,nbasis,1,0)
            forall (i=1:nbasis) NOCVene(i)=tmpmat(i,i)
            if (iopsh==1) then
                tmpmat=matmul_blas(FmatB,CObasb,nbasis,nbasis)
                tmpmat=matmul_blas(CObasb,tmpmat,nbasis,nbasis,1,0)
                forall (i=1:nbasis) NOCVene(i+nbasis)=tmpmat(i,i)
            end if
            deallocate(tmpmat)
            MOene=NOCVene
            call showNOCV(6,pairidx1,pairidx2,npair,iopsh,printthres)
        else
            write(*,*) "Energies of NOCV orbitals are not calculated since Fock matrix is unavailable"
        end if
        
    else if (isel==0) then
        call showNOCV(6,pairidx1,pairidx2,npair,iopsh,printthres)
    end if
    
    !Options involving plot
    if (isel==1) then !Show isosurface of NOCV orbitals
        write(*,*) "Note: Current orbital occupations correspond to NOCV eigenvalues"
        ishoworbsel_prt=0 !When selecting an orbital, do not print its information to avoid affecting scroll bar,&
        !because user may be viewing NOCV information printed before
        sur_value=sur_value_orb
        call showNOCV(6,pairidx1,pairidx2,npair,iopsh,printthres)
        call drawmolgui
        ishoworbsel_prt=1
        sur_value_orb=sur_value
        iorbvis=0 !Recover its status. iorbvis=0 makes saved image file has DISLIN prefix
        
    else if (isel==2.or.isel==7) then !Show/export isosurface of NOCV pair density
        do while(.true.)
            call showNOCV(6,pairidx1,pairidx2,npair,iopsh,printthres)
            write(*,*)
            write(*,*) "Input the index of the NOCV pair of interest, e.g. 2"
            write(*,"(a)") " You can also input indices of multiple pairs, e.g. 3-6,8,10-12, their densities will be summed up"
            write(*,*) "To return to post-processing menu, input ""q"""
            read(*,"(a)") c2000tmp
            if (index(c2000tmp,'q')/=0) exit
            call str2arr(c2000tmp,nsellist)
            allocate(sellist(nsellist))
            call str2arr(c2000tmp,nsellist,sellist)
            if (allocated(cubmat)) deallocate(cubmat)
	 	    allocate(cubmat(nx,ny,nz))
            if (nsellist>1) then
                write(*,"(a,f10.5)") " Sum of positive eigenvalues of selected pairs:",sum(MOocc(pairidx1(sellist(1:nsellist))))
                if (any(MOene/=0)) then
                    totpairene=0
                    do ipairtmp=1,nsellist
                        ipair=sellist(ipairtmp)
                        call NOCVpairene(ipair,pairidx1,pairidx2,npair,pairene)
                        totpairene=totpairene+pairene
                    end do
                    write(*,"(a,f12.2,' kcal/mol')") " Sum of energies of selected pairs:",totpairene
                end if
                write(*,*)
            end if
            cubmat=0
            write(*,*) "Calculating grid data..."
            allocate(wfnval(nmo))
            ifinish=0;ishowprog=1
            ntmp=floor(ny*nz/100D0)
            !$OMP PARALLEL DO SHARED(cubmat,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz,ilow,ihigh,wfnval,idx,iorb,jorb) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
            do k=1,nz
	            do j=1,ny
		            do i=1,nx
			            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                        ilow=min(minval(pairidx1(sellist(:))),minval(pairidx2(sellist(:))))
                        ihigh=max(maxval(pairidx1(sellist(:))),maxval(pairidx2(sellist(:))))
                        call orbderv(1,ilow,ihigh,tmpx,tmpy,tmpz,wfnval)
                        do idx=1,nsellist
                            iorb=pairidx1(sellist(idx))
                            jorb=pairidx2(sellist(idx))
                            if (iorb/=jorb) then
			                    cubmat(i,j,k)=cubmat(i,j,k)+MOocc(iorb)*wfnval(iorb)**2+MOocc(jorb)*wfnval(jorb)**2
                            else
			                    cubmat(i,j,k)=cubmat(i,j,k)+MOocc(iorb)*wfnval(iorb)**2
                            end if
                        end do
		            end do
		            if (ntmp/=0) then
		                !$OMP CRITICAL
                        ifinish=ifinish+1
                        ishowprog=mod(ifinish,ntmp)
		                if (ishowprog==0) call showprog(floor(100D0*ifinish/(ny*nz)),100)
                        !$OMP END CRITICAL
		            end if
	            end do
            end do
            !$OMP END PARALLEL DO
            if (ishowprog/=0) call showprog(100,100)
            deallocate(wfnval)
            if (isel==2) then
                sur_value=sur_value_iso
                call drawisosurgui(1)
                sur_value_iso=sur_value
            else if (isel==7) then
                write(*,*)
                write(*,*) "Input the file path for outputting the cube file, e.g. D:\NOCV\pair.cub"
                write(*,"(a)") " If press ENTER button directly, NOCVpair.cub will be generated in current folder"
                read(*,"(a)") c200tmp
                if (c200tmp==" ") c200tmp="NOCVpair.cub"
                write(*,*) "Outputting NOCV pair density as cube file..."
		        open(10,file=trim(c200tmp),status="replace")
		        call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		        close(10)
		        write(*,*) "Done!"
            end if
            deallocate(cubmat,sellist)
        end do
        
    else if (isel==6) then !Export cube file of a NOCV orbital
        !call genmultiorbcube !Can export a batch, but require user to again define grid
        write(*,*) "Output which NOCV orbital? Input its index, e.g. 2"
        read(*,*) iNOCV
        if (iNOCV<1.or.iNOCV>nmo) then
            write(*,*) "Error: The index exceeded valid range!"
            cycle
        end if
        if (allocated(cubmat)) deallocate(cubmat)
	 	allocate(cubmat(nx,ny,nz))
        write(*,*) "Calculating grid data..."
        call savecubmat(4,1,iNOCV)
        write(*,*)
        write(c80tmp,*) iNOCV
        c2000tmp="NOCV_"//trim(adjustl(c80tmp))//".cub"
        write(*,*) "Input the file path for outputting the cube file, e.g. D:\Akiyama\MIO.cub"
        write(*,*) "If press ENTER button directly, "//trim(c2000tmp)//" will be generated in current folder"
        read(*,"(a)") c200tmp
        if (c200tmp==" ") c200tmp=trim(c2000tmp)
        write(*,*) "Outputting grid data..."
        open(10,file=c200tmp,status="replace")
        call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
        close(10)
        write(*,*) "Done!"
        deallocate(cubmat)
        
    else if (isel==3.or.isel==4.or.isel==5.or.isel==8.or.isel==9.or.isel==10) then !Show or export Pauli/orbital deform./total deform. density
        !isel=3/8: Pauli EDD,     namely rho(frz) - rho(pro)
        !isel=4/9: Orb. def. EDD, namely rho(complex) - rho(frz)
        !isel=5/10: Total EDD,     namely rho(complex) - rho(pro)
        !
        if (iopsh==0) then !Set occupancy to MO case
            MOocc=occfrag
        else if (iopsh==1) then
            MOocc(1:nbasis)=occfrag
            MOocc(nbasis+1:nmo)=occfragB
        end if
        if (allocated(cubmat)) deallocate(cubmat)
	 	allocate(cubmat(nx,ny,nz))
        !Calculate rho(pro)
        if (isel==3.or.isel==5.or.isel==8.or.isel==10) then
            if (allocated(rho_pro).and.size(rho_pro,1)==nx.and.size(rho_pro,2)==ny.and.size(rho_pro,3)==nz) then
                continue !Do not need to recompute
            else
                if (iopsh==0) then
                    CObasa=CObasapro
                    call CObas2CO(1) !Convert CObasa to CO
                else if (iopsh==1) then
                    CObasa=CObasapro
                    CObasb=CObasbpro
                    call CObas2CO(3) !Convert CObasa and CObasb to CO
                end if
                write(*,*)
                write(*,*) "Calculating grid data of promolecular density ..."
                call savecubmat(1,0,0) !Calculate rho
                if (allocated(rho_pro)) deallocate(rho_pro)
                allocate(rho_pro(nx,ny,nz))
                rho_pro=cubmat
            end if
        end if
        !Calculate rho(frz)
        if (isel==3.or.isel==4.or.isel==8.or.isel==9) then
            if (allocated(rho_frz).and.size(rho_frz,1)==nx.and.size(rho_frz,2)==ny.and.size(rho_frz,3)==nz) then
                continue !Do not need to recompute
            else
                if (iopsh==0) then
                    CObasa=CObasafrz
                    call CObas2CO(1) !Convert CObasa to CO
                else if (iopsh==1) then
                    CObasa=CObasafrz
                    CObasb=CObasbfrz
                    call CObas2CO(3) !Convert CObasa and CObasb to CO
                end if
                write(*,*)
                write(*,*) "Calculating grid data of frozen state density ..."
                call savecubmat(1,0,0) !Calculate rho
                if (allocated(rho_frz)) deallocate(rho_frz)
                allocate(rho_frz(nx,ny,nz))
                rho_frz=cubmat
            end if
        end if
        !Calculate rho(complex)
        if (isel==4.or.isel==5.or.isel==9.or.isel==10) then
            if (allocated(rho_complex).and.size(rho_complex,1)==nx.and.size(rho_complex,2)==ny.and.size(rho_complex,3)==nz) then
                continue !Do not need to recompute
            else
                if (iopsh==0) then
                    CObasa=CObasa_org
                    call CObas2CO(1) !Convert CObasa to CO
                else if (iopsh==1) then
                    CObasa=CObasa_org
                    CObasb=CObasb_org
                    call CObas2CO(3) !Convert CObasa and CObasb to CO
                end if
                MOocc=MOocc_org
                write(*,*)
                write(*,*) "Calculating grid data of complex density ..."
                call savecubmat(1,0,0) !Calculate rho
                if (allocated(rho_complex)) deallocate(rho_complex)
                allocate(rho_complex(nx,ny,nz))
                rho_complex=cubmat
            end if
        end if
        
        if (isel==3.or.isel==8) then !Calculate Pauli EDD
            cubmat = rho_frz - rho_pro
        else if (isel==4.or.isel==9) then !Calculate orb. def. EDD
            cubmat = rho_complex - rho_frz
        else if (isel==5.or.isel==10) then !Calculate total EDD
            cubmat = rho_complex - rho_pro
        end if
        
        if (isel==3.or.isel==4.or.isel==5) then
            sur_value=sur_value_iso
        		call drawisosurgui(1)
            sur_value_iso=sur_value
        else if (isel==8.or.isel==9.or.isel==10) then
            write(*,*)
            write(*,*) "Input the file path for outputting the cube file, e.g. D:\grid.cub"
            if (isel==8) write(*,"(a)") " If press ENTER button directly, Pauli.cub will be generated in current folder"
            if (isel==9) write(*,"(a)") " If press ENTER button directly, orb_def.cub will be generated in current folder"
            if (isel==10) write(*,"(a)") " If press ENTER button directly, EDD.cub will be generated in current folder"
            read(*,"(a)") c200tmp
            if (c200tmp==" ") then
                if (isel==8) c200tmp="Pauli.cub"
                if (isel==9) c200tmp="orb_def.cub"
                if (isel==10) c200tmp="EDD.cub"
            end if
            write(*,*) "Outputting cube file..."
		    open(10,file=trim(c200tmp),status="replace")
		    call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		    close(10)
		    write(*,*) "Done!"
        end if
        deallocate(cubmat)
        if (iopsh==0) then !Restore to NOCV orbitals
            CObasa=NOCVorb
            MOocc=NOCVeig
            call CObas2CO(1) !Convert CObasa to CO
            wfntype=3
        else if (iopsh==1) then
            CObasa=NOCVorb
            MOocc(1:nbasis)=NOCVeig
            CObasb=NOCVorbB
            MOocc(nbasis+1:nmo)=NOCVeigB
            call CObas2CO(3) !Convert CObasa and CObasb to CO
            wfntype=4
        end if
        
    else if (isel==11.or.isel==12.or.isel==13) then !Show isosurface of promolecular/frozen state/complex orbitals
        MOocc=occfrag
        if (isel==13) MOocc=MOocc_org
        if (iopsh==0) then
            if (isel==11) then !Promolecular orbitals
                CObasa=CObasapro
            else if (isel==12) then !Frozen state orbitals
                CObasa=CObasafrz
            else if (isel==13) then !Show complex orbital. For convenience checking, set energies to their actual energies
                CObasa=CObasa_org
                MOene=MOene_org
            end if
            call CObas2CO(1) !Convert CObasa to CO
            wfntype=0
        else if (iopsh==1) then
            if (isel==11) then
                CObasa=CObasapro
                CObasb=CObasbpro
            else if (isel==12) then
                CObasa=CObasafrz
                CObasb=CObasbfrz
            else if (isel==13) then
                CObasa=CObasa_org
                CObasb=CObasb_org
                MOene=MOene_org
            end if
            MOocc(1:nbasis)=occfrag
            MOocc(nbasis+1:nmo)=occfragB
            call CObas2CO(3) !Convert CObasa and CObasb to CO
            wfntype=1
        end if
        if (isel==11) write(*,"(a)") " Information of all occupied promolecular orbitals are shown below. Unoccupied ones are not meaningful in the present context"
        if (isel==12) write(*,"(a)") " Information of all occupied frozen state orbitals are shown below. Unoccupied ones are not meaningful in the present context"
        if (isel==13) write(*,"(a)") " Information of all occupied actual complex molecular orbitals are shown below"
        call showorbinfo3(1) !Show occupied orbital information so that user can easily check
        ishoworbsel_prt=0
        sur_value=sur_value_orb
        if (isel==13) write(*,"(a)") " Note: Current orbital energies have been temporarily set to actual energies of complex MOs"
        call drawmolgui
        sur_value_orb=sur_value
        ishoworbsel_prt=1
        if (iopsh==0) then !Restore to NOCV orbitals
            CObasa=NOCVorb
            MOocc=NOCVeig
            call CObas2CO(1) !Convert CObasa to CO
            wfntype=3
        else if (iopsh==1) then
            CObasa=NOCVorb
            MOocc(1:nbasis)=NOCVeig
            CObasb=NOCVorbB
            MOocc(nbasis+1:nmo)=NOCVeigB
            call CObas2CO(3) !Convert CObasa and CObasb to CO
            wfntype=4
        end if
        if (isel==13) MOene=NOCVene
        
    else if (isel==14) then !Calculate composition of NOCV orbitals and pairs
        do while(.true.)
            call showNOCV(6,pairidx1,pairidx2,npair,iopsh,printthres)
            write(*,*)
            write(*,*) "Input the NOCV pair index to obtain its composition, e.g. 5"
            write(*,*) "To return, input ""q"""
            read(*,"(a)") c80tmp
            if (c80tmp=="q") then
                exit
            else
                read(c80tmp,*) ipair
            end if
            call showNOCVcomposition(6,pairidx1,pairidx2,npair,iopsh,printthres,ipair)
        end do
        
    else if (isel==15) then !
        call ask_Sbas_PBC
        write(*,*) "Variation in electron population:"
        vartot=0
        do ifrag=1,nfrag
            ibeg=ibasfrag(ifrag)
            if (ifrag==nfrag) then
                iend=nbasis
            else
                iend=ibasfrag(ifrag+1)-1
            end if
            popval=sum(Sbas(ibeg:iend,ibeg:iend)*Pdiff(ibeg:iend,ibeg:iend))
            write(*,"(' Within fragment',i3,':',f12.6)") ifrag,popval
            vartot=vartot+popval
        end do
        write(*,*)
        do ifrag=1,nfrag
            ibeg=ibasfrag(ifrag)
            if (ifrag==nfrag) then
                iend=nbasis
            else
                iend=ibasfrag(ifrag+1)-1
            end if
            do jfrag=ifrag+1,nfrag
                jbeg=ibasfrag(jfrag)
                if (jfrag==nfrag) then
                    jend=nbasis
                else
                    jend=ibasfrag(jfrag+1)-1
                end if
                popval=2*sum(Sbas(ibeg:iend,jbeg:jend)*Pdiff(ibeg:iend,jbeg:jend))
                write(*,"(' Between fragment',i3,' and',i3,':',f12.6)") ifrag,jfrag,popval
                vartot=vartot+popval
            end do
        end do
        write(*,"(/,' Sum of variation:',f12.6)") vartot
    end if
end do
end subroutine


!------ Show NOCV information
subroutine showNOCV(ides,pairidx1,pairidx2,npair,iopsh,printthres)
use defvar
implicit real*8 (a-h,o-z)
integer ides,npair,iopsh,pairidx1(npair),pairidx2(npair)
real*8 printthres

if (ides==6) write(ides,*)
write(ides,"(a,a)") "        --------------- Pair and NOCV orbital information --------------"
write(ides,"(' There are totally',i5,' NOCV pairs and',i6,' NOCV orbitals')") npair,size(MOocc)
if (printthres>0) write(ides,"(a,1PE8.1,a)") " NOCV orbitals with absolute eigenvalues smaller than",printthres," are not shown"
if (allocated(FmatA)) then
    write(ides,*) "Note: All energies are given in kcal/mol"
else
    write(ides,*) "Note: Energies of NOCV orbitals have not been evaluated, so they are all zero"
end if
if (iopsh==0) then
    write(ides,*)
    write(ides,*) " Pair  Energy | Orbital  Eigenvalue    Energy  | Orbital  Eigenvalue    Energy"
    totpairene=0
    do ipair=1,npair
        imo=pairidx1(ipair)
        jmo=pairidx2(ipair)
        if (MOocc(imo)<printthres) cycle
        if (imo/=jmo) then
            pairene=(MOocc(imo)*MOene(imo)+MOocc(jmo)*MOene(jmo))*au2kcal
            write(ides,"(i6,f8.2,2x,i6,f13.5,f11.2,3x,i6,f13.5,f11.2)") ipair,pairene,imo,MOocc(imo),MOene(imo)*au2kcal,jmo,MOocc(jmo),MOene(jmo)*au2kcal
        else
            pairene=MOocc(imo)*MOene(imo)*au2kcal
            write(ides,"(i6,f8.2,2x,i6,f13.5,f11.2)") ipair,pairene,imo,MOocc(imo),MOene(imo)*au2kcal
        end if
        totpairene=totpairene+pairene
    end do
    write(ides,"(a,f10.5)") " Sum of NOCV eigenvalues:",sum(MOocc)
    if (allocated(FmatA)) write(ides,"(a,f11.2,' kcal/mol')") " Sum of pair energies:",totpairene
else if (iopsh==1) then
    write(ides,*) "                        ----- Alpha NOCV orbitals -----"
    write(ides,*) " Pair  Energy | Orbital  Eigenvalue    Energy  | Orbital  Eigenvalue    Energy"
    totpaireneA=0
    do ipair=1,npair/2
        imo=pairidx1(ipair)
        jmo=pairidx2(ipair)
        if (MOocc(imo)<printthres) cycle
        if (imo/=jmo) then
            pairene=(MOocc(imo)*MOene(imo)+MOocc(jmo)*MOene(jmo))*au2kcal
            write(ides,"(i6,f8.2,2x,i6,f13.5,f11.2,3x,i6,f13.5,f11.2)") ipair,pairene,imo,MOocc(imo),MOene(imo)*au2kcal,jmo,MOocc(jmo),MOene(jmo)*au2kcal
        else
            pairene=MOocc(imo)*MOene(imo)*au2kcal
            write(ides,"(i6,f8.2,2x,i6,f13.5,f11.2)") ipair,pairene,imo,MOocc(imo),MOene(imo)*au2kcal
        end if
        totpaireneA=totpaireneA+pairene
    end do
    write(ides,*) "                        ----- Beta NOCV orbitals -----"
    write(ides,*) " Pair  Energy | Orbital  Eigenvalue    Energy  | Orbital  Eigenvalue    Energy"
    totpaireneB=0
    do ipair=npair/2+1,npair
        imo=pairidx1(ipair)
        jmo=pairidx2(ipair)
        if (MOocc(imo)<printthres) cycle
        if (imo/=jmo) then
            pairene=(MOocc(imo)*MOene(imo)+MOocc(jmo)*MOene(jmo))*au2kcal
            write(ides,"(i6,f8.2,2x,i6,f13.5,f11.2,3x,i6,f13.5,f11.2)") ipair,pairene,imo,MOocc(imo),MOene(imo)*au2kcal,jmo,MOocc(jmo),MOene(jmo)*au2kcal
        else
            pairene=MOocc(imo)*MOene(imo)*au2kcal
            write(ides,"(i6,f8.2,2x,i6,f13.5,f11.2)") ipair,pairene,imo,MOocc(imo),MOene(imo)*au2kcal
        end if
        totpaireneB=totpaireneB+pairene
    end do
    write(ides,"(a,f10.5,a,f10.5)") " Sum of NOCV eigenvalues:  Alpha=",sum(MOocc(1:nbasis)),"  Beta=",sum(MOocc(nbasis+1:nmo))
    if (allocated(FmatA)) write(ides,"(a,f11.2,a,f11.2,a,f11.2)") " Sum of pair energies:  Alpha=",totpaireneA,"  Beta=",totpaireneB,"  Total=",totpaireneA+totpaireneB
end if
end subroutine



!!-------- Get NOCV pair energy in kcal/mol
subroutine NOCVpairene(ipair,pairidx1,pairidx2,npair,pairene)
use defvar
implicit real*8 (a-h,o-z)
integer npair,pairidx1(npair),pairidx2(npair)
real*8 pairene
integer ipair
imo=pairidx1(ipair)
jmo=pairidx2(ipair)
if (imo/=jmo) then
    pairene=(MOocc(imo)*MOene(imo)+MOocc(jmo)*MOene(jmo))*au2kcal
else
    pairene=MOocc(imo)*MOene(imo)*au2kcal
end if
end subroutine



!------ Show composition of a NOCV pair (ipair) using SCPA partition
subroutine showNOCVcomposition(ides,pairidx1,pairidx2,npair,iopsh,printthres,ipair)
use defvar
implicit real*8 (a-h,o-z)
integer ides,npair,iopsh,pairidx1(npair),pairidx2(npair),ipair
real*8 bascomp(nbasis,nmo),lcomp_i(0:10),lcomp_j(0:10),lcomp_pair(0:10)

iorb=pairidx1(ipair)
jorb=pairidx2(ipair)

call gen_allorbbascomp_SCPA(bascomp)

write(ides,"(a,f8.3,a)") " Note: Only the basis functions and shells having absolute composition in the NOCV pair or corresponding NOCV orbitals &
&larger than",compthres," % are printed below. The threshold can be changed via ""compthres"" parameter in settings.ini"

write(ides,"(/,' Contribution of each basis function to NOCV pair/orbitals:')")
write(ides,"('  Basis    Type    Atom    Shell    Orb.',i5,'   Orb.',i5,'   Pair',i5)") iorb,jorb,ipair
do ibas=1,nbasis
    orbcompi=bascomp(ibas,iorb)*100
    orbcompj=bascomp(ibas,jorb)*100
    if (iorb/=jorb) then
        paircomp=orbcompi*MOocc(iorb)+orbcompj*MOocc(jorb) !Note that MOocc(iorb)=-MOocc(jorb) as they are paired
    else
        paircomp=orbcompi*MOocc(iorb)
    end if
	if (abs(paircomp)>compthres.or.abs(orbcompi)>compthres.or.abs(orbcompj)>compthres) then
		write(ides,"(i6,5x,a,i5,a,i5,2x,f10.2,' %',f10.2,' %',f10.2,' %')") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),'('//a(bascen(ibas))%name//')',&
		basshell(ibas),orbcompi,orbcompj,paircomp
	end if
end do

write(ides,"(/,' Contribution of each basis function shell to NOCV pair/orbitals:')")
write(ides,"('  Shell    Type    Atom         Orb.',i5,'   Orb.',i5,'   Pair',i5)") iorb,jorb,ipair
lcomp_i=0
lcomp_j=0
lcomp_pair=0
do ish=1,nshell
	shellcom_i=0D0
	shellcom_j=0D0
	shellcom_pair=0D0
	do ibas=1,nbasis
		if (basshell(ibas)==ish) then
            compiorb=bascomp(ibas,iorb)*100
            compjorb=bascomp(ibas,jorb)*100
			shellcom_i=shellcom_i+compiorb
			shellcom_j=shellcom_j+compjorb
            if (iorb/=jorb) then
                shellcom_pair=shellcom_pair+compiorb*MOocc(iorb)+compjorb*MOocc(jorb)
            else
			    shellcom_pair=shellcom_pair+compiorb*MOocc(iorb)
            end if
			iatm=bascen(ibas)
		end if
	end do
	if ((abs(shellcom_i)>compthres).or.(abs(shellcom_j)>compthres).or.(abs(shellcom_pair)>compthres)) &
    write(ides,"(i6,5x,a,3x,i6,'(',a,')',3x,f10.2,' %',f10.2,' %',f10.2,' %')") ish,shtype2name(shtype(ish)),iatm,a(iatm)%name,shellcom_i,shellcom_j,shellcom_pair
    labs=abs(shtype(ish))
    lcomp_i(labs)=lcomp_i(labs)+shellcom_i
    lcomp_j(labs)=lcomp_j(labs)+shellcom_j
    lcomp_pair(labs)=lcomp_pair(labs)+shellcom_pair
end do

write(ides,"(/,' Contribution of various types of shells to NOCV pair/orbitals:')")
write(*,"(' Type     Orb.',i5,'   Orb.',i5,'   Pair',i5)") iorb,jorb,ipair
write(*,"(a,2x,f10.2,' %',f10.2,' %',f10.2,' %')") "  s:",lcomp_i(0),lcomp_j(0),lcomp_pair(0)
write(*,"(a,2x,f10.2,' %',f10.2,' %',f10.2,' %')") "  p:",lcomp_i(1),lcomp_j(1),lcomp_pair(1)
write(*,"(a,2x,f10.2,' %',f10.2,' %',f10.2,' %')") "  d:",lcomp_i(2),lcomp_j(2),lcomp_pair(2)
write(*,"(a,2x,f10.2,' %',f10.2,' %',f10.2,' %')") "  f:",lcomp_i(3),lcomp_j(3),lcomp_pair(3)
write(*,"(a,2x,f10.2,' %',f10.2,' %',f10.2,' %')") "  g:",lcomp_i(4),lcomp_j(4),lcomp_pair(4)
write(*,"(a,2x,f10.2,' %',f10.2,' %',f10.2,' %')") "  h:",lcomp_i(5),lcomp_j(5),lcomp_pair(5)

write(ides,"(/,' Contribution of each atom to NOCV pair/orbitals:')")
write(ides,"('   Atom         Orb.',i5,'   Orb.',i5,'   Pair',i5)") iorb,jorb,ipair
do iatm=1,ncenter
    if (basstart(iatm)==0) cycle
    atmcomp_i=sum( bascomp(basstart(iatm):basend(iatm),iorb) )*100
    atmcomp_j=sum( bascomp(basstart(iatm):basend(iatm),jorb) )*100
    if (iorb/=jorb) then
        atmcomp_pair=atmcomp_i*MOocc(iorb)+atmcomp_j*MOocc(jorb)
    else
        atmcomp_pair=atmcomp_i*MOocc(iorb)
    end if
	write(ides,"(i6,'(',a,'): ',f10.2,' %',f10.2,' %',f10.2,' %')") iatm,a(iatm)%name,atmcomp_i,atmcomp_j,atmcomp_pair
end do
end subroutine