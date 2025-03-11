!-----------------------------------------
!------ Fuzzy atomic space analysis ------
!-----------------------------------------
!Normally iwork=0. If iwork=1, directly choose isel==4 to calculate delocalization index in fuzzy atomic spase (namely fuzzy bond order, see statement in JPCA,110,5108 below Eq.9) and then return
!If iwork=2, directly choose isel==8 to calculate Laplacian bond order and then return
!
!For periodic system, a special subroutine is used for integrating via even grids. For isolated systems, &
!the integration grid is directly controlled by sphpot and radpot in settings.ini, since integrand may be not proportional to electron density,
!the grid will not be adjusted automatically as proposed by Becke for more efficient integration of XC functional
!
!Note on type of integration grid for isolated systems:
!  For integrating function in atomic spaces (subfunction 1) and calculating information-theoretic aromaticity index (subfunction 12), &
!molecular grid (iintgrid=2) is always used (looping all atoms in turn, use atomic integration grids to calculate integral contribution to every atom), &
!because it makes integrating e.g. Laplacian significantly more accurate when Hirshfeld(/-I) is used because of its over-extension behavior, while the increase of computational cost is basically negligible
!For subfunction 1, option "-5 Define the atoms to be considered" doesn't affect calculation process (all atoms are calculated), but only the integral of selected atoms will be printed finally
! For calculating AOM, using molecular grid makes calculation much more expensive for large systems than atomic grid, while improvement on result is only notable &
!when diffuse functions are heavily used, so atomic grid is employed by default and can be changed by users via option -6. This is determined by iAOMgrid
! Subfunction "2 Calculate atomic and molecular multipole moments and <r^2>" uses atomic integration grid, which is satisfactory
subroutine fuzzyana(iwork)
use functions
use util
use topo
implicit real*8 (a-h,o-z)
integer iwork
integer :: iintgrid=2 !Type of grids for integrating used in subfunctions 1 and 12. =1: Atomic grid =2: Molecular grid
integer atmcalclist(ncenter),natmcalclist !The atoms to be calculated will be recorded in the array
real*8 potx(sphpot),poty(sphpot),potz(sphpot),potw(sphpot)
type(content) gridatm(radpot*sphpot)
real*8 Pvec(ncenter),rintval(ncenter,10),funcval(radpot*sphpot) !rintval store integral, can record 10 integrand at the same time
real*8 atmspcweight(radpot*sphpot) !Selected weights of current atom on grids of current center
real*8 beckeweigrid(radpot*sphpot) !Becke weight of current atom on grids of current center
real*8,allocatable :: atmspcweiarr(:,:) !Weights of all atoms on grids of current center, (grid,center)
real*8,allocatable :: alldens(:,:) !Density of all free-atoms on grids of current center, (grid,center)
!--- AOM related
integer nmatsize,nmatsizeb !Number of lowest MOs considered for total/alpha and beta
integer :: iAOMgrid=1 !Type of grids for AOM integration. =1: Atomic grid   =2: Molecular grid
real*8,allocatable :: orbvalarr(:,:)
real*8 orbval(nmo)
real*8,allocatable :: AOMtmp(:,:)
!--- Fragment related
integer DIfrag1(ncenter),DIfrag2(ncenter),nDIfrag1,nDIfrag2 !Atom list used in evaluating interfragment DI and fragment LI
integer FOM1atm(ncenter),FOM2atm(ncenter),nFOM1atm,nFOM2atm !Atom list used in evaluating one or two FOMs
real*8,allocatable :: FOM1(:,:),FOM1b(:,:),FOM2(:,:),FOM2b(:,:) !Fragment overlap matrix of fragments 1 and 2
!--- Misc
real*8 rintvalp(ncenter,10) !Private for each OpenMP thread
real*8 promol(radpot*sphpot),atomdens(radpot*sphpot),selfdens(radpot*sphpot),selfdensgrad(3,radpot*sphpot),selfdenslapl(radpot*sphpot) !For Hirshfeld partition
real*8 specrho(radpot*sphpot),specrhograd2(radpot*sphpot) !Density and its gradient^2 of atom in specific state (user-provided atomic wavefunction). Used for taking Hirshfeld as reference to calculate relative Shannon and Fisher entropy
real*8 :: covr_becke(0:nelesupp)=0D0 !Covalent radii used for Becke partition
real*8 DI(ncenter,ncenter),DIa(ncenter,ncenter),DIb(ncenter,ncenter) !Delocalization index matrix
real*8 DI_tmp(ncenter,ncenter),DIa_tmp(ncenter,ncenter),DIb_tmp(ncenter,ncenter) !Temporarily used for parallel calculation of DI
real*8 LI(ncenter),LIa(ncenter),LIb(ncenter) !Localization index array
real*8 CLRK(ncenter,ncenter) !Condensed linear response kernel
real*8 ovlpinttot(ncenter,ncenter),ovlpintpos(ncenter,ncenter),ovlpintneg(ncenter,ncenter),ovlpintpostmp(ncenter,ncenter),ovlpintnegtmp(ncenter,ncenter) !Integration between fuzzy atoms, store positive part and negative part respectively
real*8 atmmono(ncenter) !Atomic monopole, filled during multipole integration task
integer :: ifunc=3,ipartition=1,PDIatom(6),FLUatom(ncenter),FLUorb(nmo),PLRatom(6)
real*8,allocatable :: AOM(:,:,:),AOMb(:,:,:) !Total/Alpha AOM, beta AOM
real*8,allocatable :: AOMsum(:,:),AOMsumb(:,:) !AOM(i,j,k) means overlap matrix of MO i,j in atom k space
real*8 :: FLUref(nelesupp,nelesupp)=-1D0
integer :: iraddefine=-1 !-1= Specific for Laplacian bond order. 0=Custom 1=CSD 2=Pyykko 3=Suresh
integer :: nbeckeiter=3
integer :: cenind(10) !Record atom index for multicenter DI
real*8 hess(3,3),rhogradw(3)
character :: radfilename*200,selectyn,c80inp*80,specatmfilename*80,c200tmp*200,c2000tmp*2000
character(len=200) :: atmvolwfn(nelesupp)=" " !Record the wavefunction file path of each element for evaluating atomic volume
integer,allocatable :: aromatatm(:)
real*8 atmpol(ncenter),atmC6(ncenter) !Atomic polarizability and C6 estimated by TS method
real*8 quadmom(3,3),tmpvec(3),tmpmat(3,3)

!Atomic polarizability table, 2020 version https://ctcp.massey.ac.nz/index.php?group=&page=dipole&menu=dipole
real*8 :: atmpol_free(1:nelesupp)=(/ &
4.50711D0,1.38375D0,164.1125D0,37.74D0,20.5D0,11.3D0,7.4D0,5.3D0,3.74D0,2.66110D0,& !H~Ne
162.7D0,71.2D0,57.8D0,37.3D0,25D0,19.4D0,14.6D0,11.083D0,& !Na~Ar
289.7D0,160.8D0,97D0,100D0,87D0,83D0,68D0,62D0,55D0,49D0,46.5D0,38.67D0,50D0,40D0,30D0,28.9D0,21D0,16.78D0,&  !K~Kr
319.8D0,197.2D0,162D0,112D0,98D0,87D0,79D0,72D0,66D0,26.14D0,55D0,46D0,65D0,53D0,43D0,38D0,32.9D0,27.32D0,& !Rb~Xe
400.9D0,272D0,215D0,205D0,216D0,208D0,200D0,192D0,184D0,158D0,170D0,163D0,156D0,150D0,144D0,139D0,137D0,& !Cs~Lu
103D0,74D0,68D0,62D0,57D0,54D0,48D0,36D0,33.91D0,50D0,47D0,48D0,44D0,42D0,35D0,& !& !Hf~Rn
317.8D0,246D0,203D0,217D0,154D0,129D0,151D0,132D0,131D0,144D0,125D0,122D0,118D0,113D0,109D0,110D0,320D0,& !Fr~Lr
112D0,42D0,40D0,38D0,36D0,34D0,32D0,32D0,28D0,29D0,31D0,70D0,67D0,62D0,58D0,169D0,159D0,(0D0,i=121,nelesupp)/) !Rf~Ubn (120)
!Free atomic C6 dispersion coefficient, from X. Chu and A. Dalgarno, JCP, 121, 4083 (2004), and H from TABLE I of PRL 102, 073005 (2009)
real*8 :: atmC6_free(1:nelesupp)=(/ &
6.5D0,1.42D0, 1392D0,227D0,99.5D0,46.6D0,24.2D0,15.6D0,9.52D0,6.20D0,& !H~Ne
1518D0,626D0, 528D0,305D0,185D0,134D0,94.6D0,64.2D0,& !Na~Ar
3923D0,2163D0, 1383D0,1044D0,832D0,602D0,552D0,482D0,408D0,373D0,253D0,284D0, 498D0,354D0,246D0,210D0,162D0,130D0,&  !K~Kr
4769D0,3175D0, (0D0,i=39,48),779D0,659D0,492D0,445D0,385D0,& !Ca,Sr...In-I
(0D0,i=54,nelesupp)/)

if (ispecial==2) then
	ipartition=2 !Use Hirshfeld for Shubin's 2nd project
	expcutoff=1 !Full accuracy
end if
if (ifPBC/=0) ipartition=3

do i=1,ncenter !Initialize the list of the atoms to be integrated
	atmcalclist(i)=i
end do
natmcalclist=ncenter

if (all(covr_becke==0D0)) then
	if (iraddefine==-1) covr_becke=covr_tianlu !The first time
	if (iraddefine==1) covr_becke=covr
	if (iraddefine==2) covr_becke=covr_pyy
	if (iraddefine==3) covr_becke=covr_Suresh
end if
if (all(FLUref==-1D0)) then !If =-1, means missing reference value
	FLUref(6,6)=1.468D0 !Calculated for benzene under HF/6-31G* opted
	FLUref(6,7)=1.566D0 !Pyridine
	FLUref(7,6)=1.566D0
	FLUref(5,7)=1.260D0 !Borazine
	FLUref(7,5)=1.260D0
end if
!Backup original grid setting, because calculating bond order may use lower grid quality
nradpotold=radpot
nsphpotold=sphpot


!==== Interface loop ====!
!==== Interface loop ====!
do while(.true.) 

!For some functions, e.g. calculate DI, it is safe to use relatively low grid quality for saving time,
!so sphpot and radpot may be adjusted automatically, but each time enter main interface we recover the ones set by users
radpot=nradpotold
sphpot=nsphpotold
if (iwork==0) then
	write(*,*)
	write(*,*) "                ======== Fuzzy atomic space analysis ========"
	if (numcp>0) write(*,*) "-11 Choose a critical point as reference point"
    if (ifPBC==0) then
		write(*,"(a,3f10.5,' Bohr')") " -10 Set X,Y,Z of reference point, current: ",refx,refy,refz
		if (iAOMgrid==1) then
			write(*,*) "-6 Choose type of integration grid for AOM, current: Atomic grid"
		else if (iAOMgrid==2) then
			write(*,*) "-6 Choose type of integration grid for AOM, current: Molecular grid"
		end if
		if (natmcalclist==ncenter) then
			write(*,*) "-5 Define the atoms to be considered in options 1, 2, 13, current: all atoms"
		else
			write(*,"(a,i5,a)") " -5 Define the atoms to be considered in options 1, 2 and 13, current:",natmcalclist," atoms"
		end if
		write(*,*) "-4 Adjust reference parameter for FLU"
		if (ipartition==1) then !For Becke
			write(*,"(' -3 Set the number of iterations for Becke partition, current:',i3)") nbeckeiter
			if (iraddefine==-1) write(*,*) "-2 Select radius definition for Becke partition, current: Modified CSD"
			if (iraddefine==0) write(*,*) "-2 Select radius definition for Becke partition, current: Custom"
			if (iraddefine==1) write(*,*) "-2 Select radius definition for Becke partition, current: CSD"
			if (iraddefine==2) write(*,*) "-2 Select radius definition for Becke partition, current: Pyykko"
			if (iraddefine==3) write(*,*) "-2 Select radius definition for Becke partition, current: Suresh"
			if (iraddefine==4) write(*,*) "-2 Select radius definition for Becke partition, current: Hugo"
		end if
    end if
	if (ipartition==1) write(*,*) "-1 Select method for partitioning atomic spaces, current: Becke"
	if (ipartition==2) write(*,*) "-1 Select method for partitioning atomic spaces, current: Hirshfeld"
	if (ipartition==3) write(*,*) "-1 Select method for partitioning atomic spaces, current: Hirshfeld*"
	if (ipartition==4) write(*,*) "-1 Select method for partitioning atomic spaces, current: Hirshfeld-I"
	if (ipartition==5) write(*,*) "-1 Select method for partitioning atomic spaces, current: MBIS"
	write(*,*) "0 Return"
	write(*,*) "1 Perform integration in fuzzy atomic spaces for a real space function"
    if (ifPBC==0) write(*,*) "2 Calculate atomic and molecular multipole moments and <r^2>"
	write(*,*) "3 Calculate and output atomic overlap matrix (AOM) in current folder"
	write(*,*) "33 Calculate and output fragment overlap matrix (FOM) in current folder"
	write(*,*) "4 Calculate localization (LI) and delocalization index (DI)"
	write(*,*) "44 Calculate fragment LI (FLI) and interfragment DI (IFDI)"
	write(*,*) "5 Calculate PDI (Para-delocalization index)"
	write(*,*) "6 Calculate FLU (Aromatic fluctuation index)"
	write(*,*) "7 Calculate FLU-pi"
    if (ifPBC==0) write(*,*) "8 Perform integration in fuzzy overlap region for a real space functions"
	if (allocated(CObasa)) write(*,*) "9 Calculate condensed linear response kernel (CLRK)" !Need virtual orbital information
	if (allocated(CObasa)) write(*,*) "10 Calculate PLR (Para linear response index)" !Need virtual orbital information
    if (ifPBC==0) then
		write(*,*) "11 Calculate multi-center delocalization index" !Only can be used for HF/DFT closed-shell wavefunction
		write(*,*) "12 Calculate information-theoretic aromaticity index (ACS Omega, 3, 18370)"
		write(*,"(a)") " 13 Calculate atomic effective volume, free volume, polarizability and C6 coefficient"
		!!If orbital occupancy has been modified before entering fuzzy analysis module, Hirshfeld-I cannot be chosen. &
		!This option makes orbital occupancy can be changed after generating H-I information, so that H-I partition can be used to integrate functions contributed by specific orbitals
		if (ipartition==4.or.ipartition==5) write(*,*) "26 Set occupation of some orbitals"
  		!write(*,*) "101 Integrate a function in Hirshfeld atomic space with molecular grid"
		if (ispecial==2) then
			write(*,*) "99 Calculate relative Shannon and Fisher entropy and 2nd-order term"
			write(*,"(a)") " 100 Calculate relative Shannon and Fisher entropy of specific state w.r.t. Hirshfeld density"
			write(*,*) "102 Obtain quadratic and cubic Renyi entropy"
			write(*,*) "103 Obtain quadratic and cubic Renyi relative entropy"
			write(*,*) "104 The same as 99, but also calculate relative g1,g2,g3"
		end if
    end if
	read(*,*) isel
	
else if (iwork==1) then !Calculate fuzzy bond order
	isel=4 !Directly calculate delocalization index
    if (ifPBC>0) write(*,*) "Note: The fuzzy bond order will be calculated under Hirshfeld partition"
else if (iwork==2) then
	isel=8 !Directly calculate Laplacian bond order
end if


!!===================================
!!--------- Adjust settings ---------
!!===================================
if (isel==0) then
	exit
	
else if (isel==101) then
	call intHirsh_molgrid
	
else if (isel==26) then
	call modorbocc
	
else if (isel==-11) then
	if (numcp>0) then
		write(*,*) "Summary of found CPs (in Bohr):"
        write(*,*) " Index         X               Y               Z         Type"
		do icp=1,numcp
			write(*,"(i6,3f16.9,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
		end do
		write(*,*) "Select a CP by inputting its index, e.g. 5"
		read(*,*) icp
		refx=CPpos(1,icp)
		refy=CPpos(2,icp)
		refz=CPpos(3,icp)
	else
		write(*,*) "Error: No CP has been found"
	end if
	
else if (isel==-10) then
	write(*,*) "Input X,Y,Z of reference point, e.g. 3.0,-4.12,0.0"
	read(*,*) refx,refy,refz
	write(*,*) "You inputted coordinate is in which unit?  1: Bohr  2: Angstrom"
	read(*,*) iunit
	if (iunit==2) then
		refx=refx/b2a
		refy=refy/b2a
		refz=refz/b2a
	end if
	
else if (isel==-6) then
	write(*,"(a)") " Choose type of grid for calculating AOM, all quantities related to it will also be affected"
    write(*,*) "1 Atomic grid"
    write(*,*) "2 Molecular grid"
	write(*,"(a)") " Note: Molecular grid is more accurate if virtual MOs are involved and meantime &
    &diffuse functions are presented, but it is significantly more expensive"
    read(*,*) iAOMgrid
	if (allocated(AOM)) deallocate(AOM,AOMsum)
    if (allocated(AOMb)) deallocate(AOMb,AOMsumb)
    
else if (isel==-5) then
	do while(.true.)
		write(*,*) "Input atom indices, e.g. 1,3-7,9,12"
        write(*,"(a)") " Note that this setting only determines which atoms will be printed and included into statistics, while the accuracy is not affected by this"
		read(*,"(a)") c200tmp
		call str2arr(c200tmp,natmcalclist,atmcalclist)
		if (any(atmcalclist(1:natmcalclist)>ncenter).or.any(atmcalclist(1:natmcalclist)<=0)) then
			write(*,*) "One or more atoms exceeded valid range!"
		else
			exit
		end if
	end do
	write(*,*) "Done! The atoms you chosen:"
	write(*,"(10i6)") atmcalclist(1:natmcalclist)
	write(*,*)
	
else if (isel==-4) then
	write(*,*) "Current FLU reference paramters:"
	do iref=1,nelesupp
		do jref=iref,nelesupp
			if (FLUref(iref,jref)/=-1) write(*,"(' ',a,a,a,a,f10.5)") ind2name(iref),'-',ind2name(jref),':',FLUref(iref,jref)
		end do
	end do
	do while(.true.)
		write(*,*) "Input two element indices and a new reference parameter"
		write(*,*) "e.g. 6,7,1.35 means set reference parameter for C-N to 1.35"
		write(*,*) "(Input q can return)"
		read(*,"(a)") c80inp
		if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
		read(c80inp,*) itmp,jtmp,refval
		FLUref(itmp,jtmp)=refval
		FLUref(jtmp,itmp)=refval
		write(*,*) "Done!"
	end do
	
else if (isel==-3) then
	nbeckeiterold=nbeckeiter
	write(*,*) "Do how many times of iteration? e.g. 3"
	write(*,*) "Note: Larger value gives rise to sharper atomic boundary"
	read(*,*) nbeckeiter
	if (nbeckeiter/=nbeckeiterold) then
		if (allocated(AOM)) deallocate(AOM,AOMsum)
		if (allocated(AOMb)) deallocate(AOMb,AOMsumb)
	end if
	
else if (isel==-2) then
	if (allocated(AOM)) deallocate(AOM,AOMsum)
	if (allocated(AOMb)) deallocate(AOMb,AOMsumb)
	do while(.true.)
		write(*,*)
		write(*,*) "-1 Use the modified version of CSD radii defined by Tian Lu"
		write(*,*) "0 Return"
		write(*,*) "1 Use CSD radii (Dalton Trans., 2008, 2832-2838)"
		write(*,*) "2 Use Pyykko radii (Chem. Eur.-J., 15, 186-197)"
		write(*,*) "3 Use Suresh radii (J. Phys. Chem. A, 105, 5940-5944)"
		write(*,*) "4 Use Hugo radii (Chem. Phys. Lett., 480, 127-131)"
		write(*,*) "10 Read radii from external file"
		write(*,*) "11 Modify current radii by manual input"
		write(*,*) "12 Print current radii list"
		read(*,*) iselrad
		if (iselrad==-1) then
			covr_becke=covr_TianLu
			iraddefine=-1
			write(*,*) "Done!"
		else if (iselrad==0) then
			exit
		else if (iselrad==1) then
			covr_becke=covr
			iraddefine=1
			write(*,*) "Done!"
		else if (iselrad==2) then
			covr_becke=covr_pyy
			iraddefine=2
			write(*,*) "Done!"
		else if (iselrad==3) then
			covr_becke=covr_Suresh
			iraddefine=3
			write(*,*) "Done!"
		else if (iselrad==4) then
			covr_becke=radii_hugo
			iraddefine=4
			write(*,*) "Done!"
		else if (iselrad==10) then
			iraddefine=0
			write(*,"(a)") " About the file format:"
			write(*,"(a)") " The first line should be the number of elements you want to modify, followed by element indices and radii (in Angstrom), for example:"
			write(*,"(a)") " 4"
			write(*,"(a)") " 1 0.35"
			write(*,"(a)") " 4 1.2"
			write(*,"(a)") " 5 1.12"
			write(*,"(a)") " 14 1.63"
			write(*,*)
			write(*,*) "Input file path, e.g. C:\radall.txt"
			read(*,"(a)") radfilename
			inquire(file=radfilename,exist=alive)
			if (alive.eqv..true.) then
				open(10,file=radfilename,status="old")
				read(10,*) nmodrad
				do irad=1,nmodrad
					read(10,*) indtmp,radtmp
					covr_becke(indtmp)=radtmp/b2a
				end do
				close(10)
				write(*,*) "Done!"
			else
				write(*,*) "Error: File cannot be found"
			end if
		else if (iselrad==11) then
			iraddefine=0
			write(*,*) "Input element index and radius (in Angstrom), e.g. 5,0.84"
			read(*,*) indtmp,radtmp
			covr_becke(indtmp)=radtmp/b2a
			write(*,*) "Done!"
		else if (iselrad==12) then
			do irad=0,nelesupp
				write(*,"(' Element:',i5,'(',a,')   Radius:',f8.3,' Angstrom')") irad,ind2name(irad),covr_becke(irad)*b2a
			end do
		end if
	end do
	
else if (isel==-1) then
	ipartitionold=ipartition
	write(*,*) "Select atomic space partition method"
    if (ifPBC==0) then
		write(*,*) "1 Becke"
		write(*,*) "2 Hirshfeld"
		write(*,*) "3 Hirshfeld (preferred over 2)"
		write(*,*) "4 Hirshfeld-I"
		write(*,*) "5 MBIS"
		write(*,"(a)") " Note: (2) uses atomic .wfn files to calculate Hirshfeld weights, they must be provided by yourself or let Multiwfn automatically &
		&invoke Gaussian to generate them. (3) evaluates the weights based on built-in radial atomic densities, thus is more convenient than (2)"
	else
		write(*,*) "3 Hirshfeld"
		write(*,*) "4 Hirshfeld-I"
		write(*,*) "5 MBIS"
    end if
    read(*,*) ipartition
	if (imodwfn==1.and.(ipartition==2.or.ipartition==4)) then !These two modes need reloading firstly loaded file, so they cannot be already modified
		write(*,"(a)") " Error: Since the wavefunction has been modified by you or by other functions, present function is unable to use. &
		&Please reboot Multiwfn and reload the file"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
		ipartition=ipartitionold
		cycle
	end if
	if (ipartition/=ipartitionold) then
		if (allocated(AOM)) deallocate(AOM,AOMsum)
		if (allocated(AOMb)) deallocate(AOMb,AOMsumb)
	end if
	if (ipartition==4) then !Generate radial density of all atoms by Hirshfeld-I
		if (ifPBC==0) then
			call Hirshfeld_I(2)
        else
			call Hirshfeld_I_evengrid(2,1)
        end if
    else if (ipartition==5) then !Generate radial density of all atoms by MBIS
		if (ifPBC==0) then
			call MBIS(2,0)
        else
			call MBIS(2,1)
        end if
	end if
end if
if (isel==26.or.isel==101.or.isel<0) cycle


!!=======================================
!!--------- Prepare calculation ---------
!!=======================================

if (isel==1.or.isel==8) then !Select which function to be integrated in single atomic space or overlap between two atomic spaces
	if (isel==8.and.ipartition/=1) then
		write(*,"(a)") " Error: Only the fuzzy atomic space defined by Becke can be used together with this function"
		cycle
	end if
	if (iwork==2) then  !When ==2, means calculate Laplacian bond order
		ifunc=3 !Laplacian of rho
		if (iautointgrid==1) then !Allow change integration grid adapatively. Do not use sphpot=230/266, the result will be frantic, I don't know why
			radpot=45
			sphpot=302
		end if
	else
		write(*,*) "-2 Deformation density"
		call selfunc_interface(1,ifunc)
	end if
else if (isel==2) then !Multipole moment integral needs electron density
	ifunc=1
    if (ispecial==1) then
		write(*,"(a)") " Note: Since ""ispecial"" in settings.ini was set to 1, the electron density involved in this function is replaced with user-defined function"
        ifunc=100
    end if
    write(*,*) "Select outputting destination of the result"
    write(*,*) "1 Output on screen"
    write(*,*) "2 Output to multipole.txt in current folder"
    write(*,"(a)") " If choose 2, atom_moment.txt will also be outputted in current folder, &
    &which contains electric dipole moments as well as eigenvalues and eigenvectors of quadrupole moment tensors of all atoms"
    read(*,*) itmp
    if (itmp==1) then
		iout=6
    else if (itmp==2) then
		iout=20
		open(20,file="multipole.txt",status="replace")
		open(21,file="atom_moment.txt",status="replace")
    end if
	write(iout,*) "Note: All units below are in a.u."
	write(iout,*)
	
!AOM,LI/DI,PDI,FLU/-pi/CLRK/PLR/Multicenter DI. Note: MO values will be generated when collecting data
else if (isel==3.or.isel==33.or.isel==4.or.isel==44.or.isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10.or.isel==11) then
	!Even (30,110) can be used for fuzzy bond order, so (45,170) is absolutely enough
	if (iautointgrid==1) then !Allow changing integration grid adapatively
		radpot=45
		sphpot=170
	end if
	!FLU,FLU-pi,PDI,PLR,CLRK are only applied to closed-shell system
	if (isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10) then
		if (wfntype/=0.and.wfntype/=3) then
			write(*,*) "Error: This function is only available for closed-shell system!"
			cycle
		end if
	else if (isel==11.and.wfntype/=0) then
		write(*,"(a)") " Error: This function is only available for single-determinant closed-shell system!"
		cycle
	end if
    if (isel==9.or.isel==10) then
		do imo=1,nmo
			if (MOocc(imo)==0.and.MOene(imo)==0) then
				write(*,*) "Error: Calculation of CLRK and PLR needs virtual orbitals, however currently one or &
                &more virtual orbitals have zero energy! Perhaps virtual orbitals were not solved"
                write(*,*) "Press ENTER button to continue, the result will be meaningless"
                read(*,*)
                exit
            end if
        end do
    end if
    if (isel==33.and.ifPBC==0) then !Calculate FOM. In the case of PBC, the fragment(s) is defined in subroutine AOMFOM_evengrid
		write(*,*)
		write(*,*) "Select the way of calculating fragment overlap matrix (FOM)"
		write(*,*) "1 Calculate and export FOM for one fragment"
		write(*,*) "2 Calculate and export FOM for two specific fragments"
		read(*,*) iFOMmode
		write(*,*) "Input indices of atoms for defining fragment 1, e.g. 3-9,14,19-20"
		read(*,"(a)") c2000tmp
		call str2arr(c2000tmp,nFOM1atm,FOM1atm)
		if (iFOMmode==2) then
			write(*,*) "Input indices of atoms for defining fragment 2, e.g. 1,2,10-13,15-18"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,nFOM2atm,FOM2atm)
		end if		
    end if
    if (isel==44) then !Calculate fragment LI and interfragment DI
		write(*,*) "Input atom indices to define fragment 1, e.g. 3,5-8,15-20"
		read(*,"(a)") c2000tmp
		call str2arr(c2000tmp,nDIfrag1,DIfrag1)
		write(*,*) "Input atom indices to define fragment 2, e.g. 1,2,4,9-14"
		read(*,"(a)") c2000tmp
		call str2arr(c2000tmp,nDIfrag2,DIfrag2)
    end if
	!Allocate space for AOM
	if ((isel==9.or.isel==10).and.allocated(AOM).and.size(AOM,1)==nmo) then !The AOM calculated for FLU/PDI is smaller than nmo, since virtual orbitals are not taken into account
		goto 10
	else if ((isel==5.or.isel==6.or.isel==7).and.(allocated(AOM))) then !Do not calculate AOM again. If AOM is calculated for PLR, the AOM is also applicable since it is more than necessary
		write(*,*) "Note: AOM has already been generated before, so skipping its calculation"
		goto 10
	else !Have not calculated proper AOM, hence needed to be calculated this time
		if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
			if (wfntype==3.or.isel==9.or.isel==10.or.ispecial==3) then !R-post-HF or CLRK or PLR, or requested by users, need to consider all orbitals
				nmatsize=nmo
                write(*,*) "Note: All orbitals are taken into account in AOM calculation"
			else !Single determinant, only consider occupied MOs
				!Notice that occupation number may be not contiguous, some low-lying orbital may have &
				!zero occupation due to modification by users, so we cannot simply use nelec to determine matrix size
				do nmatsize=nmo,1,-1
					if (MOocc(nmatsize)/=0) exit
				end do
				if (nmo-nmatsize>0) then
					write(*,"(' Note: The highest',i6,' virtual orbitals will not be taken into account')") nmo-nmatsize
                    if (isel==3.or.isel==33) write(*,*) "If you hope to consider all orbitals, set ""ispecial"" in settings.ini to 3"
                end if
			end if
            write(*,*) "Allocating memory..."
			if (allocated(AOM)) deallocate(AOM,AOMsum) !For PLR, the previous AOM and AOMsum allocated by PDI/FLU is too small, so here should be released
			allocate(AOM(nmatsize,nmatsize,ncenter),AOMsum(nmatsize,nmatsize))
            if (allocated(AOMtmp)) deallocate(AOMtmp) 
            allocate(AOMtmp(nmatsize,nmatsize))
			AOM=0
			AOMsum=0
		else if (wfntype==1.or.wfntype==4) then !UHF, U-post-HF
			do iendalpha=nmo,1,-1
				if (MOtype(iendalpha)==1) exit
			end do
			if (wfntype==4.or.ispecial==3) then !U-post-HF, or requested by users, need to consider all orbitals
				nmatsizea=iendalpha
				nmatsizeb=nmo-nmatsizea
                write(*,*) "Note: All orbitals are taken into account in AOM calculation"
			else !Single determinant, only consider occupied MOs
				do nmatsizea=iendalpha,1,-1
					if (MOocc(nmatsizea)/=0D0) exit
				end do
				if (nint(nbelec)==0) then
					nmatsizeb=0
				else
					do nmatsizeb=nmo,iendalpha+1,-1
						if (MOocc(nmatsizeb)/=0D0) exit
					end do
					nmatsizeb=nmatsizeb-iendalpha
				end if
                ntmpa=iendalpha-nmatsizea
                ntmpb=nmo-iendalpha-nmatsizeb
				if (ntmpa>0) write(*,"(' Note: The highest',i6,' alpha virtual orbitals will not be taken into account')") iendalpha-nmatsizea
				if (ntmpb>0) write(*,"('       The highest',i6,' beta virtual orbitals will not be taken into account')") nmo-iendalpha-nmatsizeb
                if ((isel==3.or.isel==33).and.(ntmpa>0.or.ntmpb>0)) write(*,*) "If you hope to consider all orbitals, set ""ispecial"" in settings.ini to 3"
			end if
			if (allocated(AOM)) deallocate(AOM,AOMsum) 
			if (allocated(AOMb)) deallocate(AOMb,AOMsumb)
			allocate(AOM(nmatsizea,nmatsizea,ncenter),AOMb(nmatsizeb,nmatsizeb,ncenter))
			allocate(AOMsum(nmatsizea,nmatsizea),AOMsumb(nmatsizeb,nmatsizeb))
            if (allocated(AOMtmp)) deallocate(AOMtmp) 
            allocate(AOMtmp(nmatsizea,nmatsizea)) !Because nmatsizea >= nmatsizeb
			AOM=0
			AOMb=0
			AOMsum=0
			AOMsumb=0
		end if
	end if
else if (isel==12) then
	write(*,*) "Choose one of information-theoretic quantity"
	write(*,*) "1 Shannon entropy"
	write(*,*) "2 Fisher information"
	write(*,*) "3 Ghosh-Berkowitz-Parr (GBP) entropy"
	read(*,*) infoaromat
	write(*,*) "Input the atoms constituing the ring, e.g. 2,3,7-10"
	read(*,"(a)") c2000tmp
	call str2arr(c2000tmp,naromatatm)
	if (allocated(aromatatm)) deallocate(aromatatm)
	allocate(aromatatm(naromatatm))
	call str2arr(c2000tmp,naromatatm,aromatatm)
else if (isel==13) then !Calculating atomic volume needs electron density
	ifunc=1
    do iele=1,nelesupp
		if (any(a(atmcalclist(1:natmcalclist))%index==iele)) then
			if (atmvolwfn(iele)==" ") then
				write(*,"(/,a)") " Input path of the file containing free-state wavefunction of element "//trim(ind2name(iele))
                write(*,*) "For example, D:\atmvol\"//trim(ind2name(iele))//".wfn (.fch, .mwfn, .molden... are also acceptable)"
                do while(.true.)
					read(*,"(a)") c200tmp
					inquire(file=c200tmp,exist=alive)
					if (alive) then
						atmvolwfn(iele)=c200tmp
						exit
                    end if
					write(*,*) "Error: Cannot find the file, input again!"
                end do
            else
				write(*,"(1x,a)") trim(atmvolwfn(iele))//" will be used for element "//trim(ind2name(iele))
            end if
        end if
    end do
end if



!!=======================================
!!---------------------------------------
!!--------- Start calculation -----------
!!---------------------------------------
!!=======================================

!For periodic case, use a specific subroutine to integrate based on even grids
if (ifPBC/=0) then
	if (isel==1) then !Integrate real space function in atomic spaces
		call intatmspace_evengrid(ipartition,ifunc,rintval(1:ncenter,1))
		goto 110 !Directly jump to result statistics
    else if (isel==3.or.isel==33.or.isel==4.or.isel==44.or.isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10) then
		if (isel==33.or.isel==44) then !isel==33: Calculate and export FOM using even grids. isel==44: Calculate AOM for atoms only involved in computing fragment LI/DI later for saving cost
			call AOMFOM_evengrid(ipartition,isel,AOM,AOMb,size(AOM,1),size(AOMb,1),iendalpha,DIfrag1,DIfrag2,nDIfrag1,nDIfrag2)
		else !Calculate full AOM
			call AOMFOM_evengrid(ipartition,iwork,AOM,AOMb,size(AOM,1),size(AOMb,1),iendalpha,DIfrag1,DIfrag2,nDIfrag1,nDIfrag2)
		end if
        goto 100 !Jump to the position just after cycling atoms for integration
    end if
end if

rintval=0D0 !Initialize accumulated variables
ovlpintpos=0D0
ovlpintneg=0D0
atmmono=0
xinttot=0D0;yinttot=0D0;zinttot=0D0 !For calculating system multiple moment
xxinttot=0D0;yyinttot=0D0;zzinttot=0D0;xyinttot=0D0;yzinttot=0D0;xzinttot=0D0
xxxinttot=0D0;yyyinttot=0D0;zzzinttot=0D0;yzzinttot=0D0;xzzinttot=0D0;xxzinttot=0D0;yyzinttot=0D0;xxyinttot=0D0;xyyinttot=0D0;xyzinttot=0D0

if (ipartition==2.or.ifunc==-2) call setpromol !In this routine reload first molecule at the end
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
write(*,*) "Please wait..."
write(*,*)
call walltime(nwalltime1)

call Lebedevgen(sphpot,potx,poty,potz,potw)

!Allocate arrays for molecular grid calculation
if ((iintgrid==2.or.iAOMgrid==2).and.(.not.allocated(atmspcweiarr))) allocate(atmspcweiarr(ncenter,radpot*sphpot),alldens(ncenter,radpot*sphpot))
if (iAOMgrid==2.and.(.not.allocated(orbvalarr))) allocate(orbvalarr(nmo,radpot*sphpot))

if ((isel==2.and.iout==20).or.(isel/=2.and.isel/=13)) then
	ifinish=0
	call showprog(0,natmcalclist)
end if

!! Cycle each atom !!!! Cycle each atom !!!! Cycle each atom !!!! Cycle each atom !!
!! Cycle each atom !!!! Cycle each atom !!!! Cycle each atom !!!! Cycle each atom !!
do iatm=1,ncenter

    if (isel==2.or.isel==13) then !These two functions employs atomic integration grid, nonrelated atoms are skipped and not printed during looping
		if (all(atmcalclist(1:natmcalclist)/=iatm)) cycle
    else if (isel==33.and.iAOMgrid==1) then !Calculate FOM using atomic grid
		if (iFOMmode==1) then
			if (all(FOM1atm(1:nFOM1atm)/=iatm)) cycle
		else
			if (all(FOM1atm(1:nFOM1atm)/=iatm).and.all(FOM2atm(1:nFOM2atm)/=iatm)) cycle
        end if
    else if (isel==44.and.iAOMgrid==1) then !Calculate interfragment DI and fragment LI
        if (all(DIfrag1(1:nDIfrag1)/=iatm).and.all(DIfrag2(1:nDIfrag2)/=iatm)) cycle
    end if
  !  if (iintgrid==1) then !When using atomic integration grid, non-involved atoms can be skipped for saving cost. This is not used, because currently always use molecular grid
		!if (isel==1.and.all(atmcalclist(1:natmcalclist)/=iatm)) cycle
  !  end if
    
	!!! Prepare grid points on current center
	iradcut=0 !Before where the radial points will be cut
	parm=1D0
	do i=1,radpot !Combine spherical point&weights with second kind Gauss-Chebyshev method for radial part
		radx=cos(i*pi/(radpot+1))
		radr=(1+radx)/(1-radx)*parm !Becke transform
		radw=2*pi/(radpot+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
		gridatm( (i-1)*sphpot+1:i*sphpot )%x=radr*potx
		gridatm( (i-1)*sphpot+1:i*sphpot )%y=radr*poty
		gridatm( (i-1)*sphpot+1:i*sphpot )%z=radr*potz
		gridatm( (i-1)*sphpot+1:i*sphpot )%value=radw*potw
		if (radcut/=0D0.and.iradcut==0.and.radr<radcut) iradcut=i-1
	end do
! 	if (radcut/=0) write(*,"(i8,' points per center were discarded')") iradcut*sphpot
	gridatm%x=gridatm%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatm%y+a(iatm)%y
	gridatm%z=gridatm%z+a(iatm)%z
    
    !! Calculate function value
	!For integrating real space function (1,8), calculate selected function value at each point here
	!For multipole moment integration (2) and evaluating atomic volume (13), calculate electron density here
    !For information-theoretic aromaticity index (12), calculate selected information density at each point here
	!For AOM, LI/DI, FLU/PDI calculation, wavefunction value will be computed in integration stage rather than here
	if (isel==1.or.isel==2.or.isel==8.or.isel==13.or.isel==102) then
		!$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
		do i=1+iradcut*sphpot,radpot*sphpot
			if (ifunc==-2.or.isel==102) then !Deformation density and Renyi entropy require molecular electron density
				funcval(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
			else
				funcval(i)=calcfuncall(ifunc,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
			end if
		end do
		!$OMP end parallel do
		
		!Calculate deformation density. We have calculated total density, thus now minusing it by each atomic density in free-state
		if ((isel==1.or.isel==8).and.ifunc==-2) then
			do jatm=1,ncenter_org !Calculate free atomic density
				call dealloall(0)
				call readwfn(custommapname(jatm),1)
				!$OMP parallel do shared(atomdens) private(ipt) num_threads(nthreads)
				do ipt=1+iradcut*sphpot,radpot*sphpot
					atomdens(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
				end do
				!$OMP end parallel do
				funcval=funcval-atomdens
			end do
			call dealloall(0)
			call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
		end if
    else if (isel==12) then !Information-theoretic aromaticity index
        !$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
		do i=1+iradcut*sphpot,radpot*sphpot
            if (infoaromat==1) then
                funcval(i)=infoentro(2,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
            else if (infoaromat==2) then
                funcval(i)=Fisherinfo(1,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
            else if (infoaromat==3) then
				funcval(i)=Ghoshentro(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,2)
			end if
		end do
		!$OMP end parallel do
	end if
	
	!!! Calculate atomic space weight function at all points around the atom (recorded in atmspcweight), which will be used later
	!Also integrate fuzzy overlap region here (only available for Becke partition)
	if (ipartition==1) then !Becke partition
		!$OMP parallel shared(atmspcweight,atmspcweiarr,ovlpintpos,ovlpintneg) private(i,rnowx,rnowy,rnowz,&
		!$OMP ii,jj,Pvec,tmpval,tmpval2,ovlpintpostmp,ovlpintnegtmp) num_threads(nthreads)
		ovlpintpostmp=0D0
		ovlpintnegtmp=0D0
		!$OMP do schedule(dynamic)
		do i=1+iradcut*sphpot,radpot*sphpot
			rnowx=gridatm(i)%x
			rnowy=gridatm(i)%y
			rnowz=gridatm(i)%z
            
			!Calculate Becke weight of all atoms (Pvec) at current point
            call BeckePvec(rnowx,rnowy,rnowz,Pvec,covr_becke,nbeckeiter)
            
			atmspcweight(i)=Pvec(iatm)
            if (allocated(atmspcweiarr)) atmspcweiarr(:,i)=Pvec(:)
			
			if (isel==8) then !Integration between two fuzzy atoms
				tmpval=Pvec(iatm)*funcval(i)*gridatm(i)%value
				do ii=1,ncenter !Note, ovlpint is lower triangular matrix, will be convert to full matrix during statistic stage
					do jj=ii,ncenter
						tmpval2=Pvec(jj)*Pvec(ii)*tmpval
						if (tmpval>0) then !ovlpinttot will be summed up to single value in statistic stage
							ovlpintpostmp(jj,ii)=ovlpintpostmp(jj,ii)+tmpval2
						else
							ovlpintnegtmp(jj,ii)=ovlpintnegtmp(jj,ii)+tmpval2
						end if
					end do
				end do
			end if
		end do
		!$OMP end do
		!$OMP CRITICAL
			ovlpintpos=ovlpintpos+ovlpintpostmp
			ovlpintneg=ovlpintneg+ovlpintnegtmp
		!$OMP end CRITICAL
		!$OMP end parallel
	else if (ipartition==2.or.ipartition==3.or.ipartition==4.or.ipartition==5) then !Hirshfeld(-/I), MBIS
		promol=0D0
		do jatm=1,ncenter_org !Calculate free atomic density of each atom and promolecular density
			if (ipartition==2) then !!Hirshfeld partition based on atomic .wfn files
				call dealloall(0)
				call readwfn(custommapname(jatm),1)
				!$OMP parallel do shared(atomdens,selfdensgrad,selfdenslapl) private(ipt,hess) num_threads(nthreads)
				do ipt=1+iradcut*sphpot,radpot*sphpot
					atomdens(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
					if (jatm==iatm) then !SPECIAL: For Shubin's project, evaluate derivative of rho
						if (isel==99.or.isel==100) then !Calculate rho and its gradient for free atom
							call calchessmat_dens(1,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,selfdens(ipt),selfdensgrad(1:3,ipt),hess)
						else if (isel==104) then !Calculate rho, its gradient and Hessian for free atom
							call calchessmat_dens(2,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,selfdens(ipt),selfdensgrad(1:3,ipt),hess)
							selfdenslapl(ipt)=hess(1,1)+hess(2,2)+hess(3,3)
						end if
					end if
				end do
				!$OMP end parallel do
            else if (ipartition==3) then !Hirshfeld based on interpolation of built-in atomic radius density
				!$OMP parallel do shared(atomdens) private(ipt) num_threads(nthreads)
				do ipt=1+iradcut*sphpot,radpot*sphpot
					atomdens(ipt)=calcatmdens(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
				end do
				!$OMP end parallel do
			else !Hirshfeld-I or MBIS based on refined atomic radial density
				!$OMP parallel do shared(atomdens) private(ipt) num_threads(nthreads)
				do ipt=1+iradcut*sphpot,radpot*sphpot
					atomdens(ipt)=fdens_rad(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
				end do
				!$OMP end parallel do
            end if
			promol(:)=promol(:)+atomdens(:)
			if (jatm==iatm) selfdens(:)=atomdens(:)
            if (iintgrid==2.or.iAOMgrid==2) alldens(jatm,:)=atomdens(:)
		end do
        !Note: Overlap between Rydberg orbitals can extend to extremely distant region, where promol is zero if interpolated density is used, &
        !in this case grid weights of present center are best set to 1 (though not rigorous), so that the orbital overlap can be integrated &
        !everywhere and thus satisfies orthogonal condition. Essentially, this corresponds to integrating zero-density region in Becke space of present center
		do i=1+iradcut*sphpot,radpot*sphpot !Get Hirshfeld（/-I) weight of present atom
			if (promol(i)/=0D0) then
				atmspcweight(i)=selfdens(i)/promol(i)
			else
				atmspcweight(i)=1D0
			end if
            if (iintgrid==2.or.iAOMgrid==2) then
				if (promol(i)/=0D0) then
					atmspcweiarr(:,i)=alldens(:,i)/promol(i)
				else
					atmspcweiarr(:,i)=0D0
                    atmspcweiarr(iatm,i)=1D0
				end if
            end if
		end do
        if (ipartition==2) then !Retrieve the firstly loaded file(whole molecule) in order to calculate real rho later
			call dealloall(0)
			call readinfile(firstfilename,1)
        end if
	end if
    
    !SPECIAL CASE: Calculate density of atom in specific-state, for Shubin's idea
    if (isel==100) then
		write(specatmfilename,"(a,i4.4,a)") "specwfn/",iatm,".wfn"
		write(*,*) "Processing "//trim(specatmfilename)
		call dealloall(0)
		call readwfn(specatmfilename,1)
		a=a_org(iatm) !Set atom position to actual atom position
		!$OMP parallel do shared(specrho,specrhograd2) private(ipt) num_threads(nthreads)
		do ipt=1+iradcut*sphpot,radpot*sphpot
		    specrho(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
            specrhograd2(ipt)=fgrad(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,'t')**2
		end do
		!$OMP end parallel do
		call dealloall(0)
		call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule) to calc real rho later
    end if
    
	!!! Perform integration
	if (isel==1.or.isel==12) then !Integration in atomic spaces, and information-theoretic aromaticity index
		if (iintgrid==1) then !Atomic grid
			do i=1+iradcut*sphpot,radpot*sphpot
				rintval(iatm,1)=rintval(iatm,1)+atmspcweight(i)*funcval(i)*gridatm(i)%value
			end do
		else if (iintgrid==2) then !Molecular grid. Calculate contribution of current grids to all atoms
			call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3) !Calculate Becke weight
			do i=1+iradcut*sphpot,radpot*sphpot
				rintval(:,1)=rintval(:,1)+atmspcweiarr(:,i)*funcval(i)*beckeweigrid(i)*gridatm(i)%value
			end do
        end if
	else if (isel==99.or.isel==100.or.isel==103.or.isel==104) then !SPECIAL SPECIAL SPECIAL
		!=99:  Calculate relative Shannon and Fisher entropy and 2nd-order term
		!=100: Calculate relative Shannon/Fisher entropy by taking Hirshfeld density as reference
		!=103: Calculate quadratic and cubic Renyi relative entropy
		!=104: Same as 99, but also calculate relative g1, g2, g3
		!$OMP parallel shared(rintval) private(i,rnowx,rnowy,rnowz,rhow,rhogradw,rhograd2w,hess,rholaplw,rintvalp,tmpx,tmpy,tmpz) num_threads(nthreads)
		rintvalp=0D0
		!$OMP do schedule(dynamic)
		do i=1+iradcut*sphpot,radpot*sphpot
			rnowx=gridatm(i)%x
			rnowy=gridatm(i)%y
			rnowz=gridatm(i)%z
			if (isel==99.or.isel==100) then
				call calchessmat_dens(1,rnowx,rnowy,rnowz,rhow,rhogradw(:),hess)
				rhow=atmspcweight(i)*rhow !rhoA
				rhogradw=atmspcweight(i)*rhogradw !grad_rhoA
				rhograd2w=sum(rhogradw(:)**2) !|grad_rhoA|^2
			else if (isel==104) then
				call calchessmat_dens(2,rnowx,rnowy,rnowz,rhow,rhogradw(:),hess)
				rhow=atmspcweight(i)*rhow !rhoA
				rhogradw=atmspcweight(i)*rhogradw !grad_rhoA
				rhograd2w=sum(rhogradw(:)**2) !|grad_rhoA|^2
                rholaplw=atmspcweight(i)*(hess(1,1)+hess(2,2)+hess(3,3))
			else
				rhow=atmspcweight(i)*fdens(rnowx,rnowy,rnowz) !rhoA at current point
			end if
			if (isel==99.or.isel==104) then
			    !Relative Shannon entropy w.r.t. free-state
			    rintvalp(iatm,1)=rintvalp(iatm,1)+rhow*log(rhow/selfdens(i))*gridatm(i)%value
			    !Relative Fisher information entropy w.r.t. free-state (old formula, incorrect)
			    rintvalp(iatm,2)=rintvalp(iatm,2)+(rhograd2w/rhow-sum(selfdensgrad(1:3,i)**2)/selfdens(i))*gridatm(i)%value
			    !Atomic Shannon
			    rintvalp(iatm,3)=rintvalp(iatm,3)-rhow*log(rhow)*gridatm(i)%value
			    !Atomic Fisher
			    rintvalp(iatm,4)=rintvalp(iatm,4)+(rhograd2w/rhow)*gridatm(i)%value
			    !1st-order term: rhoA-rho0
			    rintvalp(iatm,5)=rintvalp(iatm,5)+(rhow-selfdens(i))*gridatm(i)%value
			    !2nd-order term
			    rintvalp(iatm,6)=rintvalp(iatm,6)+(rhow-selfdens(i))**2/rhow/2D0 *gridatm(i)%value
			    !Relative Fisher information entropy w.r.t. free-state (new formula, correct)
			    tmpx=rhogradw(1)/rhow-selfdensgrad(1,i)/selfdens(i)
			    tmpy=rhogradw(2)/rhow-selfdensgrad(2,i)/selfdens(i)
			    tmpz=rhogradw(3)/rhow-selfdensgrad(3,i)/selfdens(i)
			    rintvalp(iatm,7)=rintvalp(iatm,7)+rhow*(tmpx**2+tmpy**2+tmpz**2)*gridatm(i)%value
                if (isel==104) then !Relative g1, g2, g3, respectively
                    rintvalp(iatm,8)=rintvalp(iatm,8) + rholaplw*log(rhow/selfdens(i))*gridatm(i)%value
                    rintvalp(iatm,9)=rintvalp(iatm,9) + rhow*(rholaplw/rhow-selfdenslapl(i)/selfdens(i))*gridatm(i)%value
                    tmpx=rhogradw(1)/rhow-selfdensgrad(1,i)/selfdens(i)
                    tmpy=rhogradw(2)/rhow-selfdensgrad(2,i)/selfdens(i)
                    tmpz=rhogradw(3)/rhow-selfdensgrad(3,i)/selfdens(i)
                    rintvalp(iatm,10)=rintvalp(iatm,10) + rhow*(tmpx**2+tmpy**2+tmpz**2)*gridatm(i)%value
                end if
			else if (isel==100) then
			    !Relative Shannon entropy of specific atomic state with Hirshfeld density as reference
			    rintvalp(iatm,1)=rintvalp(iatm,1)+specrho(i)*log(specrho(i)/rhow)*gridatm(i)%value
			    !Relative Fisher information entropy of specific atomic state with Hirshfeld density as reference
			    rintvalp(iatm,2)=rintvalp(iatm,2)+(specrhograd2(i)/specrho(i)-rhograd2w/rhow)*gridatm(i)%value
			else if (isel==103) then
				!Quadratic Renyi relative entropy 
				rintvalp(iatm,1)=rintvalp(iatm,1)+rhow**2/selfdens(i) *gridatm(i)%value
				!Cubic Renyi relative entropy 
				rintvalp(iatm,2)=rintvalp(iatm,2)+rhow**3/selfdens(i)**2 *gridatm(i)%value
			end if
		end do
		!$OMP end do
		!$OMP CRITICAL
			rintval=rintval+rintvalp
		!$OMP end CRITICAL
		!$OMP end parallel
            
	else if (isel==102) then !SPECIAL SPECIAL SPECIAL: Obtain quadratic and cubic Renyi entropy
		do i=1+iradcut*sphpot,radpot*sphpot
			rintval(iatm,1)=rintval(iatm,1)+atmspcweight(i)*funcval(i)**2*gridatm(i)%value
			rintval(iatm,2)=rintval(iatm,2)+atmspcweight(i)*funcval(i)**3*gridatm(i)%value
		end do
        
	else if (isel==2) then !Integrate multipole moments
		eleint=0D0
		xint=0D0;yint=0D0;zint=0D0
		xxint=0D0;yyint=0D0;zzint=0D0;xyint=0D0;yzint=0D0;xzint=0D0
		xxxint=0D0;yyyint=0D0;zzzint=0D0;yzzint=0D0;xzzint=0D0;xxzint=0D0;yyzint=0D0;xxyint=0D0;xyyint=0D0;xyzint=0D0
		do i=1+iradcut*sphpot,radpot*sphpot
			tmpmul=atmspcweight(i)*funcval(i)*gridatm(i)%value
			eleint=eleint-tmpmul !Monopole
            !Calculate atomic dipole/quadruple/octopole, the coordinates are w.r.t. nucleus
			rx=gridatm(i)%x-a(iatm)%x
			ry=gridatm(i)%y-a(iatm)%y
			rz=gridatm(i)%z-a(iatm)%z
			xint=xint+rx*tmpmul !Dipole moments
			yint=yint+ry*tmpmul
			zint=zint+rz*tmpmul
			xxint=xxint+rx*rx*tmpmul !Quadruple moments
			yyint=yyint+ry*ry*tmpmul
			zzint=zzint+rz*rz*tmpmul
			xyint=xyint+rx*ry*tmpmul
			yzint=yzint+ry*rz*tmpmul
			xzint=xzint+rx*rz*tmpmul
			xxxint=xxxint+rx*rx*rx*tmpmul !Octopole moments
			yyyint=yyyint+ry*ry*ry*tmpmul
			zzzint=zzzint+rz*rz*rz*tmpmul
			yzzint=yzzint+ry*rz*rz*tmpmul
			xzzint=xzzint+rx*rz*rz*tmpmul
			xxzint=xxzint+rx*rx*rz*tmpmul
			yyzint=yyzint+ry*ry*rz*tmpmul
			xxyint=xxyint+rx*rx*ry*tmpmul
			xyyint=xyyint+rx*ry*ry*tmpmul
			xyzint=xyzint+rx*ry*rz*tmpmul
            !Calculate molecular dipole/quadruple/octopole, the coordinates are w.r.t. origin
			rx=gridatm(i)%x
            ry=gridatm(i)%y
            rz=gridatm(i)%z
			xinttot=xinttot-rx*tmpmul
			yinttot=yinttot-ry*tmpmul
			zinttot=zinttot-rz*tmpmul
			xxinttot=xxinttot-rx*rx*tmpmul
			yyinttot=yyinttot-ry*ry*tmpmul
			zzinttot=zzinttot-rz*rz*tmpmul
			xyinttot=xyinttot-rx*ry*tmpmul
			yzinttot=yzinttot-ry*rz*tmpmul
			xzinttot=xzinttot-rx*rz*tmpmul
			xxxinttot=xxxinttot-rx*rx*rx*tmpmul
			yyyinttot=yyyinttot-ry*ry*ry*tmpmul
			zzzinttot=zzzinttot-rz*rz*rz*tmpmul
			yzzinttot=yzzinttot-ry*rz*rz*tmpmul
			xzzinttot=xzzinttot-rx*rz*rz*tmpmul
			xxzinttot=xxzinttot-rx*rx*rz*tmpmul
			yyzinttot=yyzinttot-ry*ry*rz*tmpmul
			xxyinttot=xxyinttot-rx*rx*ry*tmpmul
			xyyinttot=xyyinttot-rx*ry*ry*tmpmul
			xyzinttot=xyzinttot-rx*ry*rz*tmpmul
		end do
		rrint=xxint+yyint+zzint
		rrxint=xxxint+xyyint+xzzint
		rryint=xxyint+yyyint+yzzint
		rrzint=xxzint+yyzint+zzzint
		if (nEDFelec==0) then !Normal case, all electron basis or using pseudopotential but not accompanied by EDF
			atmchgtmp=a(iatm)%charge+eleint
		else !EDF is used for some atoms. Core electron density represented by EDF has been integrated, so nuclear charge should be augmented by nEDFelecatm
			atmchgtmp=a(iatm)%charge+nEDFelecatm(iatm)+eleint
		end if
		atmmono(iatm)=eleint
		write(iout,"('                           *****  Atom',i6,'(',a2,')  *****')") iatm,a(iatm)%name
		write(iout,"(' Atomic charge:',f12.6)") atmchgtmp
		write(iout,"(' Atomic monopole moment (from electrons):',f12.6)") eleint
		write(iout,"(' Atomic dipole moments:')")
		write(iout,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") -xint,-yint,-zint,dsqrt(xint**2+yint**2+zint**2)
		write(iout,"(' Contribution to molecular dipole moment:')")
		contridipx=atmchgtmp*a(iatm)%x-xint
		contridipy=atmchgtmp*a(iatm)%y-yint
		contridipz=atmchgtmp*a(iatm)%z-zint
		write(iout,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") contridipx,contridipy,contridipz,dsqrt(contridipx**2+contridipy**2+contridipz**2)
        if (ispecial==1) then
			write(iout,"(' Atomic quadrupole moments (Standard Cartesian form):')")
			write(iout,"(' XX=',f12.6,'  XY=',f12.6,'  XZ=',f12.6)") -xxint,-xyint,-xzint
			write(iout,"(' YX=',f12.6,'  YY=',f12.6,'  YZ=',f12.6)") -xyint,-yyint,-yzint
			write(iout,"(' ZX=',f12.6,'  ZY=',f12.6,'  ZZ=',f12.6)") -xzint,-yzint,-zzint
        end if
		write(iout,"(' Atomic quadrupole moments (Traceless Cartesian form):')")
		QXX=(-3*xxint+rrint)/2
		QYY=(-3*yyint+rrint)/2
		QZZ=(-3*zzint+rrint)/2
        quadmom(1,1)=QXX
        quadmom(1,2)=(-3*xyint)/2
        quadmom(1,3)=(-3*xzint)/2
        quadmom(2,1)=(-3*xyint)/2
        quadmom(2,2)=QYY
        quadmom(2,3)=(-3*yzint)/2
        quadmom(3,1)=(-3*xzint)/2
        quadmom(3,2)=(-3*yzint)/2
        quadmom(3,3)=QZZ
		write(iout,"(' XX=',f12.6,'  XY=',f12.6,'  XZ=',f12.6)") quadmom(1,:)
		write(iout,"(' YX=',f12.6,'  YY=',f12.6,'  YZ=',f12.6)") quadmom(2,:)
		write(iout,"(' ZX=',f12.6,'  ZY=',f12.6,'  ZZ=',f12.6)") quadmom(3,:)
		write(iout,"(' Magnitude of the traceless quadrupole moment tensor:',f12.6)") sqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
		R20=-(3*zzint-rrint)/2D0 !Notice that the negative sign, because electrons carry negative charge
		R2n1=-dsqrt(3D0)*yzint
		R2p1=-dsqrt(3D0)*xzint
		R2n2=-dsqrt(3D0)*xyint
		R2p2=-dsqrt(3D0)/2D0*(xxint-yyint)
		write(iout,"(' Atomic quadrupole moments (Spherical harmonic form):')")
		write(iout,"(' Q_2,0 =',f11.6,'   Q_2,-1=',f11.6,'   Q_2,1=',f11.6)") R20,R2n1,R2p1
		write(iout,"(' Q_2,-2=',f11.6,'   Q_2,2 =',f11.6)") R2n2,R2p2
		write(iout,"( ' Magnitude: |Q_2|=',f12.6)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)
        AESEx=xxint
        AESEy=yyint
        AESEz=zzint
        AESE=rrint
		write(iout,"(a,f14.6)") " Atomic electronic spatial extent <r^2>:",AESE
		write(iout,"(' Components of <r^2>:  X=',f14.6,'  Y=',f14.6,'  Z=',f14.6)") AESEx,AESEy,AESEz
		R30=-(5*zzzint-3*rrzint)/2D0
		R3n1=-dsqrt(3D0/8D0)*(5*yzzint-rryint)
		R3p1=-dsqrt(3D0/8D0)*(5*xzzint-rrxint)
		R3n2=-dsqrt(15D0)*xyzint
		R3p2=-dsqrt(15D0)*(xxzint-yyzint)/2D0
		R3n3=-dsqrt(5D0/8D0)*(3*xxyint-yyyint)
		R3p3=-dsqrt(5D0/8D0)*(xxxint-3*xyyint)
		write(iout,"(' Atomic octopole moments (Spherical harmonic form):')")
		write(iout,"(' Q_3,0 =',f11.6,'  Q_3,-1=',f11.6,'  Q_3,1 =',f11.6)") R30,R3n1,R3p1
		write(iout,"(' Q_3,-2=',f11.6,'  Q_3,2 =',f11.6,'  Q_3,-3=',f11.6,'  Q_3,3 =',f11.6)") R3n2,R3p2,R3n3,R3p3
		write(iout,"( ' Magnitude: |Q_3|=',f12.6)") dsqrt(R30**2+R3n1**2+R3p1**2+R3n2**2+R3p2**2+R3n3**2+R3p3**2)
		write(iout,*)
        !Output dipole and quadrupole moments
        if (iout==20) then
			write(21,"(' Atom',i6,' (',a,')')") iatm,a(iatm)%name
			write(21,"(' Atomic dipole moment:',3f12.6)") -xint,-yint,-zint
			call diagsymat(quadmom,tmpmat,tmpvec,istat)
			write(21,"(' Information of atomic quadrupole moment (Traceless Cartesian):')")
			write(21,"(' Eigenvalue 1:',f12.6,'  Eigenvector:',3f12.6)") tmpvec(1),tmpmat(:,1)
			write(21,"(' Eigenvalue 2:',f12.6,'  Eigenvector:',3f12.6)") tmpvec(2),tmpmat(:,2)
			write(21,"(' Eigenvalue 3:',f12.6,'  Eigenvector:',3f12.6)") tmpvec(3),tmpmat(:,3)
			write(21,*)
        end if
    
	else if (isel==13) then !Calculate atomic volume and related quantities
        !Calculate effective volume
		effV=0
		do i=1+iradcut*sphpot,radpot*sphpot
			rx=gridatm(i)%x-a(iatm)%x
			ry=gridatm(i)%y-a(iatm)%y
			rz=gridatm(i)%z-a(iatm)%z
            r3=dsqrt(rx**2+ry**2+rz**2)**3
			tmp=funcval(i)*atmspcweight(i) * r3 *gridatm(i)%value !Current funcval records electron density
            effV=effV+tmp
		end do
        !Calculate free volume
        freeV=0
        gridatm%x=gridatm%x-a(iatm)%x !Recover the integration points to (0,0,0) as center
		gridatm%y=gridatm%y-a(iatm)%y
		gridatm%z=gridatm%z-a(iatm)%z
		call dealloall(0)
		call readinfile(atmvolwfn(a_org(iatm)%index),1)
        gridatm%x=gridatm%x+a(1)%x !Move the points to current single atom as center
		gridatm%y=gridatm%y+a(1)%y
		gridatm%z=gridatm%z+a(1)%z
		do i=1+iradcut*sphpot,radpot*sphpot
			rx=gridatm(i)%x-a(1)%x
			ry=gridatm(i)%y-a(1)%y
			rz=gridatm(i)%z-a(1)%z
            r3=dsqrt(rx**2+ry**2+rz**2)**3
			tmp=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z) * r3 *gridatm(i)%value
            freeV=freeV+tmp
        end do
		call dealloall(0)
		call readinfile(firstfilename,1) !Retrieve to the whole molecule wavefunction
        write(*,"(' Atom',i5,'(',a,')  Effective V:',f10.3,'  Free V:',f10.3,' a.u.  Ratio:',f6.3)") &
        iatm,a(iatm)%name,effV,freeV,effV/freeV
        atmpol(iatm)=atmpol_free(a(iatm)%index)*effV/freeV
        atmC6(iatm)=atmC6_free(a(iatm)%index)*(effV/freeV)**2
	
	!Calculate atomic overlap matrix (AOM) for all tasks that require it
	else if (isel==3.or.isel==33.or.isel==4.or.isel==44.or.isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10.or.isel==11) then
		!Calculate total or alpha part
		if (wfntype==1.or.wfntype==4) nmatsize=nmatsizea !UHF,U-post-HF
		if (iAOMgrid==1) then !Atomic grid
			!$OMP parallel shared(AOM) private(i,imo,jmo,AOMtmp,orbval,tmpval) num_threads(nthreads)
			AOMtmp=0D0
			!$OMP do schedule(dynamic)
			do i=1+iradcut*sphpot,radpot*sphpot
				call orbderv(1,1,nmatsize,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,orbval)
				do jmo=1,nmatsize
					tmpval=atmspcweight(i)*orbval(jmo)*gridatm(i)%value
					do imo=jmo,nmatsize
						AOMtmp(imo,jmo)=AOMtmp(imo,jmo)+tmpval*orbval(imo)
					end do
				end do
			end do
			!$OMP end do
			!$OMP CRITICAL
				AOM(:,:,iatm)=AOM(:,:,iatm)+AOMtmp(1:nmatsize,1:nmatsize)
			!$OMP end CRITICAL
			!$OMP end parallel
        else if (iAOMgrid==2) then !Molecular grid
			call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3) !Calculate Becke weight
			!$OMP parallel do shared(orbvalarr) private(i) num_threads(nthreads)
			do i=1+iradcut*sphpot,radpot*sphpot
				call orbderv(1,1,nmatsize,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,orbvalarr(1:nmatsize,i))
			end do
			!$OMP end parallel do
			!$OMP parallel do shared(AOM) private(jatm,i,imo,jmo,tmpval,tmpval2) num_threads(nthreads)
			do jatm=1,ncenter
				do i=1+iradcut*sphpot,radpot*sphpot
					tmpval=beckeweigrid(i)*atmspcweiarr(jatm,i)*gridatm(i)%value
					do jmo=1,nmatsize
						tmpval2=tmpval*orbvalarr(jmo,i)
						do imo=jmo,nmatsize
							AOM(imo,jmo,jatm)=AOM(imo,jmo,jatm)+tmpval2*orbvalarr(imo,i)
						end do
					end do
				end do
			end do
			!$OMP end parallel do
        end if
		do jmo=1,nmatsize
			do imo=jmo+1,nmatsize
				AOM(jmo,imo,:)=AOM(imo,jmo,:)
			end do
		end do
        
        !Calculate Beta part for UHF,U-post-HF
		if ((wfntype==1.or.wfntype==4).and.nmatsizeb>0) then
			MOinit=iendalpha+1
			MOend=iendalpha+nmatsizeb
            if (iAOMgrid==1) then !Atomic grid
				!$OMP parallel shared(AOMb) private(i,imo,jmo,AOMtmp,orbval,tmpval) num_threads(nthreads)
				AOMtmp=0D0
				!$OMP do schedule(dynamic)
				do i=1+iradcut*sphpot,radpot*sphpot
					call orbderv(1,MOinit,MOend,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,orbval) !Calculate orbital wavefunction value of all MOs in current position and store to orbval
					do jmo=MOinit,MOend
						tmpval=atmspcweight(i)*orbval(jmo)*gridatm(i)%value
						do imo=jmo,MOend
							AOMtmp(imo-iendalpha,jmo-iendalpha)=AOMtmp(imo-iendalpha,jmo-iendalpha)+tmpval*orbval(imo)
						end do
					end do
				end do
				!$OMP end do
				!$OMP CRITICAL
					AOMb(1:nmatsizeb,1:nmatsizeb,iatm)=AOMb(1:nmatsizeb,1:nmatsizeb,iatm)+AOMtmp(MOinit-iendalpha:MOend-iendalpha,MOinit-iendalpha:MOend-iendalpha)
				!$OMP end CRITICAL
				!$OMP end parallel
			else if (iAOMgrid==2) then !Molecular grid
				call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3) !Calculate Becke weight
				!$OMP parallel do shared(orbvalarr) private(i) num_threads(nthreads)
				do i=1+iradcut*sphpot,radpot*sphpot
					call orbderv(1,MOinit,MOend,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,orbvalarr(1:nmatsizeb,i))
				end do
				!$OMP end parallel do
				!$OMP parallel do shared(AOMb) private(jatm,i,imo,jmo,imotmp,jmotmp,tmpval,tmpval2) num_threads(nthreads)
				do jatm=1,ncenter
					do i=1+iradcut*sphpot,radpot*sphpot
						tmpval=beckeweigrid(i)*atmspcweiarr(jatm,i)*gridatm(i)%value
						do jmo=MOinit,MOend
							tmpval2=tmpval*orbvalarr(jmo,i)
                            jmotmp=jmo-iendalpha
							do imo=jmo,MOend
								imotmp=imo-iendalpha
								AOMb(imotmp,jmotmp,jatm)=AOMb(imotmp,jmotmp,jatm)+tmpval2*orbvalarr(imo,i)
							end do
						end do
					end do
				end do
				!$OMP end parallel do
			end if
			do jmo=1,nmatsizeb
				do imo=jmo+1,nmatsizeb
					AOMb(jmo,imo,:)=AOMb(imo,jmo,:)
				end do
			end do
        end if
	end if
    
	!Show progress
    ifinish=ifinish+1
	if (isel==2.and.iout==20) then !Output multipole moments to multipole.txt and atom_moment.txt
        call showprog(ifinish,natmcalclist)
    else if (isel==33.and.iAOMgrid==1) then !Calculate FOM
		if (iFOMmode==1) then
			call showprog(ifinish,nFOM1atm)
		else
			call showprog(ifinish,nFOM1atm+nFOM2atm)
        end if
    else if (isel==44.and.iAOMgrid==1) then !Calculate interfragment DI and fragment LI
        call showprog(ifinish,nDIfrag1+nDIfrag2)
    else if (isel/=2.and.isel/=13) then !For evaluating atomic multipole moments (isel=2) and atomic volumes (isel=13), the process is not shown, because result is directly printed during looping
        call showprog(iatm,ncenter)
    end if
	
end do !End cycling atoms

call walltime(nwalltime2)
write(*,"(' Calculation took up',i8,' seconds wall clock time')") nwalltime2-nwalltime1

100 continue

!==== Check sanity of AOM ====!
if (isel==3.or.isel==4.or.isel==5.or.isel==6.or.isel==7.or.isel==9.or.isel==10.or.isel==11) then
	iwarn=0
	if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
		do iatm=1,ncenter
			AOMsum=AOMsum+AOM(:,:,iatm)
		end do
		AOMerror=identmaterr(AOMsum)/ncenter
		write(*,"(' Error of AOM is',f14.8)") AOMerror
		call identmatmaxerr(AOMsum,errdiag,idiag,errndiag,indiag,jndiag)
		write(*,"(' Maximum diagonal deviation to 1:   ',f10.6,' at orbital',i6)") errdiag,idiag
		write(*,"(' Maximum nondiagonal deviation to 0:',f10.6,' between orbitals',2i6)") errndiag,indiag,jndiag
        if (AOMerror>0.001D0) iwarn=1
	else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
		AOMerrorb=0D0
		do iatm=1,ncenter
			AOMsum=AOMsum+AOM(:,:,iatm)
			if (nmatsizeb>0) AOMsumb=AOMsumb+AOMb(:,:,iatm)
		end do
		AOMerrora=identmaterr(AOMsum)/ncenter
		write(*,"(' Error of alpha AOM is',f14.8)") AOMerrora
        if (AOMerrora>0.001D0) iwarn=1
		if (nmatsizeb>0) then
			AOMerrorb=identmaterr(AOMsumb)/ncenter
			write(*,"(' Error of Beta AOM is ',f14.8)") AOMerrorb
			if (AOMerrorb>0.001D0) iwarn=1
        end if
    end if
	if (iwarn==1) then
		write(*,"(/,a)") " Warning: The integration is not very accurate"
        if (ifPBC==0) then
			write(*,*) "To improve accuracy, please try one or some of following treatments:"
			write(*,*) "(1) Enlarge ""radpot"" and ""sphpot"" in settings.ini"
			write(*,*) "(2) Set ""radcut"" in settings.ini to 0"
			write(*,*) "(3) Choose option -6 to change to much more expensive molecular grid"
			write(*,*) "(4) If diffuse functions were heavily employed, remove them"
        else
			write(*,*) "To improve accuracy, please reduce grid spacing for integration"
        end if
		!open(10,file="AOMsum.txt",status="replace")
 	!	call showmatgau(AOMsum,"AOMsum",1,"f14.8",10)
		!close(10)
    end if
end if

!==== Generate DI, LI or condensed linear response kernel (CLRK) ====!
!DI-pi will be calculated for FLU-pi at later stage
!Multicenter DI will be calculated at later stage
10	if (isel==4.or.isel==44.or.isel==5.or.isel==6) then !For LI/DI, PDI and FLU
	if (any(MOocc<0)) then
		where(MOocc<0) MOocc=0
		write(*,"(a)") " Note: Some occupation numbers are negative. In order to make the calculation feasible, they have been set to zero"
		write(*,*) "Press ENTER button to continue"
		read(*,*)
	end if
    write(*,*)
    write(*,*) "Calculating LI and DI from AOM..."
	!RHF,R-post-HF, DI_A,B=2∑[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B     where i and j are non-spin orbitals
	if (wfntype==0.or.wfntype==3) then
		DI=0D0
        !$OMP parallel shared(DI) private(iatm,jatm,iorb,jorb,DI_tmp) num_threads(nthreads)
        DI_tmp=0
		!$OMP DO schedule(dynamic)
		do iatm=1,ncenter
			do jatm=iatm,ncenter
				do iorb=1,nmatsize
					do jorb=1,nmatsize
						DI_tmp(iatm,jatm)=DI_tmp(iatm,jatm)+dsqrt(MOocc(iorb)*MOocc(jorb))*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
					end do
				end do
			end do
		end do
		!$OMP END DO
		!$OMP CRITICAL
        DI=DI+DI_tmp
		!$OMP END CRITICAL
		!$OMP END PARALLEL
        forall(iatm=1:ncenter) LI(iatm)=DI(iatm,iatm)
		DI=2*(DI+transpose(DI))
		do iatm=1,ncenter !Diagonal terms are the sum of corresponding row or column
			DI(iatm,iatm)=0D0
			DI(iatm,iatm)=sum(DI(iatm,:))
		end do
	else if (wfntype==2) then !ROHF
		DIa=0D0
		DIb=0D0
		do nmoclose=nmatsize,1,-1
			if (MOtype(nmoclose)==0) exit
		end do
        !$OMP parallel shared(DIa,DIb) private(iatm,jatm,iorb,jorb,occi,occj,DIa_tmp,DIb_tmp) num_threads(nthreads)
        DIa_tmp=0
        DIb_tmp=0
		!$OMP DO schedule(dynamic)
		do iatm=1,ncenter
			do jatm=iatm,ncenter
				!Alpha
				do iorb=1,nmatsize !The number of close or alpha orbitals needed to be concerned
					occi=MOocc(iorb)
					if (MOtype(iorb)==0) occi=occi/2D0
					do jorb=1,nmatsize
						occj=MOocc(jorb)
						if (MOtype(jorb)==0) occj=occj/2D0
						DIa(iatm,jatm)=DIa(iatm,jatm)+dsqrt(occi*occj)*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
					end do
				end do
				!Beta
				do iorb=1,nmoclose !The number of close orbitals needed to be concerned
					do jorb=1,nmoclose
						DIb(iatm,jatm)=DIb(iatm,jatm)+dsqrt(MOocc(iorb)/2D0*MOocc(jorb)/2D0)*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
					end do
				end do
			end do
		end do
		!$OMP END DO
		!$OMP CRITICAL
        DIa=DIa+DIa_tmp
        DIb=DIb+DIb_tmp
		!$OMP END CRITICAL
		!$OMP END PARALLEL
        forall(iatm=1:ncenter) LIa(iatm)=DIa(iatm,iatm)
        forall(iatm=1:ncenter) LIb(iatm)=DIb(iatm,iatm)
		DIa=2*(DIa+transpose(DIa))
		DIb=2*(DIb+transpose(DIb))
		do iatm=1,ncenter !Diagonal terms are the sum of corresponding row or column
			DIa(iatm,iatm)=0D0
			DIb(iatm,iatm)=0D0
			DIa(iatm,iatm)=sum(DIa(iatm,:))
			DIb(iatm,iatm)=sum(DIb(iatm,:))
		end do
		!Combine alpha and Beta to total
		DI=DIa+DIb
		LI=LIa+LIb
	!UHF,U-post-HF   DI(A,B)=2∑[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B   where i and j are spin orbitals
	else if (wfntype==1.or.wfntype==4) then
		!Alpha
		DIa=0D0
        !$OMP parallel shared(DIa) private(iatm,jatm,iorb,jorb,DIa_tmp) num_threads(nthreads)
        DIa_tmp=0
		!$OMP DO schedule(dynamic)
		do iatm=1,ncenter
			do jatm=iatm,ncenter
				do iorb=1,nmatsizea
					do jorb=1,nmatsizea
						DIa(iatm,jatm)=DIa(iatm,jatm)+dsqrt(MOocc(iorb)*MOocc(jorb))*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
					end do
				end do
			end do
		end do
		!$OMP END DO
		!$OMP CRITICAL
        DIa=DIa+DIa_tmp
		!$OMP END CRITICAL
		!$OMP END PARALLEL
        forall(iatm=1:ncenter) LIa(iatm)=DIa(iatm,iatm)
		DIa=2*(DIa+transpose(DIa))
		!Beta
		if (nmatsizeb>0) then
			DIb=0D0
			MOinit=iendalpha+1 !Index range of beta orbitals
			MOend=iendalpha+nmatsizeb
			!$OMP parallel shared(DIb) private(iatm,jatm,iorb,iorbtmp,jorb,jorbtmp,DIb_tmp) num_threads(nthreads)
			DIb_tmp=0
			!$OMP DO schedule(dynamic)
			do iatm=1,ncenter
				do jatm=iatm,ncenter
					do iorb=MOinit,MOend
						iorbtmp=iorb-iendalpha
						do jorb=MOinit,MOend
							jorbtmp=jorb-iendalpha
							DIb(iatm,jatm)=DIb(iatm,jatm)+dsqrt(MOocc(iorb)*MOocc(jorb))*AOMb(iorbtmp,jorbtmp,iatm)*AOMb(iorbtmp,jorbtmp,jatm)
						end do
					end do
				end do
			end do
			!$OMP END DO
			!$OMP CRITICAL
			DIb=DIb+DIb_tmp
			!$OMP END CRITICAL
			!$OMP END PARALLEL
			forall(iatm=1:ncenter) LIb(iatm)=DIb(iatm,iatm)
			DIb=2*(DIb+transpose(DIb))
		end if
		do iatm=1,ncenter !Diagonal terms are the sum of corresponding row or column
			DIa(iatm,iatm)=0D0
			DIb(iatm,iatm)=0D0
			DIa(iatm,iatm)=sum(DIa(iatm,:))
			DIb(iatm,iatm)=sum(DIb(iatm,:))
		end do
		!Combine alpha and Beta to total
		DI=DIa+DIb
		LI=LIa+LIb
	end if
	
else if (isel==9.or.isel==10) then !Calculate condensed linear response kernel, PLR also uses it
	CLRK=0D0
	do iatm=1,ncenter
		do jatm=iatm,ncenter
			do iorb=1,nmo !Occupied MOs
				if (nint(MOocc(iorb))==2D0) then
					do jorb=idxHOMO+1,nmo !Virtual MOs
						if (nint(MOocc(jorb))==0D0) CLRK(iatm,jatm)=CLRK(iatm,jatm)+AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)/(MOene(iorb)-MOene(jorb))
					end do
				end if
			end do
		end do
	end do
	CLRK=CLRK*4D0
	CLRK=CLRK+transpose(CLRK)
	do iatm=1,ncenter
		CLRK(iatm,iatm)=CLRK(iatm,iatm)/2D0
	end do
end if




!!====================================================
!!------- Statistic results or post-processing -------
!!====================================================
110 write(*,*)
if (isel==1) then
	do iatm=1,ncenter
		if ( all(atmcalclist(1:natmcalclist)/=iatm) ) rintval(iatm,1)=0
    end do
	sumval=sum(rintval(atmcalclist(1:natmcalclist),1))
	sumabsval=sum(abs(rintval(atmcalclist(1:natmcalclist),1)))
	write(*,*) "  Atomic space        Value                % of sum            % of sum abs"
	if (any(abs(rintval(atmcalclist(1:natmcalclist),1))>1D9).or.all(abs(rintval(atmcalclist(1:natmcalclist),1))<1D-7)) then
		do idx=1,natmcalclist
			iatm=atmcalclist(idx)
			write(*,"(i6,'(',a2,')  ',E20.10,1x,f20.6,1x,f20.6)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,1)/sumval*100,rintval(iatm,1)/sumabsval*100
		end do
		write(*,"(' Summing up above values:',E20.10)") sumval
		write(*,"(' Summing up absolute value of above values:',E20.10)") sumabsval
    else
		do idx=1,natmcalclist
			iatm=atmcalclist(idx)
            write(*,"(i6,'(',a2,')  ',f20.8,1x,f20.6,1x,f20.6)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,1)/sumval*100,rintval(iatm,1)/sumabsval*100
		end do
		write(*,"(' Summing up above values:',f20.8)") sumval
		write(*,"(' Summing up absolute value of above values:',f20.8)") sumabsval
	end if
	
else if (isel==99.or.isel==104) then !SPECIAL: Relative Shannon and Fisher entropy and 2nd-order term
	write(*,*) "Relative Shannon entropy and relative Fisher information w.r.t. its free-state"
	write(*,*) "   Atom           Rel.Shannon       Rel.Fisher(old)   Rel.Fisher(new)"
	do iatm=1,ncenter
		write(*,"(i6,'(',a2,')  ',3f18.8)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2),rintval(iatm,7)
	end do
	write(*,"(' Summing up above values:',3f18.8)") sum(rintval(:,1)),sum(rintval(:,2)),sum(rintval(:,7))
	write(*,*)
	write(*,*) "Shannon and Fisher information entropy of each atom"
	write(*,*) "   Atom             Shannon            Fisher"
	do iatm=1,ncenter
		write(*,"(i6,'(',a2,')  ',2f18.8)") iatm,a(iatm)%name,rintval(iatm,3),rintval(iatm,4)
	end do
	write(*,"(' Summing up above values:',2f22.8)") sum(rintval(:,3)),sum(rintval(:,4))
	write(*,*)
	write(*,*) "1st and 2nd-order terms of each atom"
	write(*,*) "   Atom           1st           2nd"
	do iatm=1,ncenter
		write(*,"(i6,'(',a2,')  ',2f14.8)") iatm,a(iatm)%name,rintval(iatm,5),rintval(iatm,6)
	end do
	write(*,"(' Summing up above values:',2f16.8)") sum(rintval(:,5)),sum(rintval(:,6))
	write(*,*)
    if (isel==104) then
	    write(*,*) "Relative g1, g2 and g3"
        write(*,*) "    Atom            rel. g1           rel. g2           rel. g3"
	    do iatm=1,ncenter
		    write(*,"(i6,'(',a2,')  ',3f18.8)") iatm,a(iatm)%name,rintval(iatm,8),rintval(iatm,9),rintval(iatm,10)
	    end do
	    write(*,"(' Summing up above values:',3f18.8)") sum(rintval(:,8)),sum(rintval(:,9)),sum(rintval(:,10))
	end if
else if (isel==100) then !SPECIAL: Relative Shannon/Fisher by taking Hirshfeld density as reference
	write(*,*) "Relative Shannon and Fisher entropy of specific state w.r.t. Hirshfeld density"
	write(*,*) "   Atom         Relat_Shannon      Relat_Fisher"
	do iatm=1,ncenter
		write(*,"(i6,'(',a2,')  ',2f18.8)") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2)
	end do
else if (isel==102) then !SPECIAL: Quadratic and cubic Renyi entropy
	write(*,*) "Atomic contribution to int(rho^2) and int(rho^3) under Hirshfeld partition:"
	write(*,*) "   Atom            Quadratic             Cubic"
	do iatm=1,ncenter
		write(*,"(i6,'(',a2,')  ',2(1PE20.8))") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2)
	end do
	write(*,"('    Total   ',2(1PE20.8))") sum(rintval(:,1)),sum(rintval(:,2))
	write(*,*)
	write(*,"(' Molecular quadratic Renyi entropy:',f18.8)") -log10(sum(rintval(:,1)))
	write(*,"(' Molecular cubic Renyi entropy:    ',f18.8)") -log10(sum(rintval(:,2)))/2
else if (isel==103) then !SPECIAL: Quadratic and cubic Renyi relative entropy
	write(*,"(a)") " Note: rhoA=w_A(r)*rho(r) is density of A in molecule, rhoA0 is density of A in its free-state"
	write(*,*) "   Atom        int(rhoA^2/rhoA0)   int(rhoA^3/rhoA0^2)"
	do iatm=1,ncenter
		write(*,"(i6,'(',a2,')  ',2(1PE20.8))") iatm,a(iatm)%name,rintval(iatm,1),rintval(iatm,2)
	end do
	write(*,"('    Total   ',2(1PE20.8))") sum(rintval(:,1)),sum(rintval(:,2))
	write(*,*)
	write(*,"(' Molecular quadratic Renyi relative entropy:',f18.8)") -log10(sum(rintval(:,1)))
	write(*,"(' Molecular cubic Renyi relative entropy:    ',f18.8)") -log10(sum(rintval(:,2)))
		
else if (isel==2) then !Multipole moment
    if (natmcalclist/=ncenter) then
        write(iout,"(a)") " Note: The word ""Molecular"" in this context corresponds to the fragment you defined by option -5"
    else
        write(iout,"(a)") " Note: The word ""Molecular"" in this context corresponds to the entire current system"
    end if
    write(iout,*) "             *****  Molecular dipole and multipole moments  *****"
	write(iout,"(' Total number of electrons:',f14.6,'   Net charge:',f12.6)") -sum(atmmono),sum(a(atmcalclist(1:natmcalclist))%index)+sum(atmmono)
    
    ESEx=-xxinttot
    ESEy=-yyinttot
    ESEz=-zzinttot
    ESE=ESEx+ESEy+ESEz
	!Combine nuclear contribution and electron contribution to obtain molecular multiple moment
    do idx=1,natmcalclist
        iatm=atmcalclist(idx)
        xinttot=xinttot+a(iatm)%x*a(iatm)%index
        yinttot=yinttot+a(iatm)%y*a(iatm)%index
        zinttot=zinttot+a(iatm)%z*a(iatm)%index
        xxinttot=xxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%index
        yyinttot=yyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%index
        zzinttot=zzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%index
        xyinttot=xyinttot+a(iatm)%x*a(iatm)%y*a(iatm)%index
        yzinttot=yzinttot+a(iatm)%y*a(iatm)%z*a(iatm)%index
        xzinttot=xzinttot+a(iatm)%x*a(iatm)%z*a(iatm)%index
		xxxinttot=xxxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%index
		yyyinttot=yyyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%index
		zzzinttot=zzzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%index
		yzzinttot=yzzinttot+a(iatm)%y*a(iatm)%z*a(iatm)%z*a(iatm)%index
		xzzinttot=xzzinttot+a(iatm)%x*a(iatm)%z*a(iatm)%z*a(iatm)%index
		xxzinttot=xxzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%z*a(iatm)%index
		yyzinttot=yyzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%z*a(iatm)%index
		xxyinttot=xxyinttot+a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%index
		xyyinttot=xyyinttot+a(iatm)%x*a(iatm)%y*a(iatm)%y*a(iatm)%index
		xyzinttot=xyzinttot+a(iatm)%x*a(iatm)%y*a(iatm)%z*a(iatm)%index
    end do
	rrinttot=xxinttot+yyinttot+zzinttot
	rrxinttot=xxxinttot+xyyinttot+xzzinttot
	rryinttot=xxyinttot+yyyinttot+yzzinttot
	rrzinttot=xxzinttot+yyzinttot+zzzinttot
    
	write(iout,"(' Molecular dipole moment (a.u.): ',3f14.6)") xinttot,yinttot,zinttot
	write(iout,"(' Molecular dipole moment (Debye):',3f14.6)") xinttot*au2debye,yinttot*au2debye,zinttot*au2debye
	dipmag=sqrt(xinttot**2+yinttot**2+zinttot**2)
	write(iout,"(' Magnitude of molecular dipole moment (a.u.&Debye):',2f14.6)") dipmag,dipmag*au2debye
    write(iout,"(' Molecular quadrupole moments (Standard Cartesian form):')")
    fac=1
    !fac=au2debye*b2a !If using this factor, result will be identical to "Quadrupole moment (field-independent basis, Debye-Ang):" printed by Gaussian
	write(iout,"(' XX=',f12.4,'  XY=',f12.4,'  XZ=',f12.4)") xxinttot*fac,xyinttot*fac,xzinttot*fac
	write(iout,"(' YX=',f12.4,'  YY=',f12.4,'  YZ=',f12.4)") xyinttot*fac,yyinttot*fac,yzinttot*fac
	write(iout,"(' ZX=',f12.4,'  ZY=',f12.4,'  ZZ=',f12.4)") xzinttot*fac,yzinttot*fac,zzinttot*fac
    write(iout,"(' Molecular quadrupole moments (Traceless Cartesian form):')")
    !If removing the comment, the data will be identical to "Traceless Quadrupole moment (field-independent basis, Debye-Ang)" printed by Gaussian
	QXX=(3*xxinttot-rrinttot)/2 !*au2debye*b2a/1.5D0
	QYY=(3*yyinttot-rrinttot)/2 !*au2debye*b2a/1.5D0
	QZZ=(3*zzinttot-rrinttot)/2 !*au2debye*b2a/1.5D0
    QXY=3*xyinttot/2            !*au2debye*b2a/1.5D0
    QXZ=3*xzinttot/2            !*au2debye*b2a/1.5D0
    QYZ=3*yzinttot/2            !*au2debye*b2a/1.5D0
	write(iout,"(' XX=',f12.4,'  XY=',f12.4,'  XZ=',f12.4)") QXX,QXY,QXZ
	write(iout,"(' YX=',f12.4,'  YY=',f12.4,'  YZ=',f12.4)") QXY,QYY,QYZ
	write(iout,"(' ZX=',f12.4,'  ZY=',f12.4,'  ZZ=',f12.4)") QXZ,QYZ,QZZ
	write(iout,"(' Magnitude of the traceless quadrupole moment tensor:',f12.6)") sqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
	R20=(3*zzinttot-rrinttot)/2D0 !Notice that the negative sign, because electrons carry negative charge
	R2n1=dsqrt(3D0)*yzinttot
	R2p1=dsqrt(3D0)*xzinttot
	R2n2=dsqrt(3D0)*xyinttot
	R2p2=dsqrt(3D0)/2D0*(xxinttot-yyinttot)
	write(iout,"(' Molecular quadrupole moments (Spherical harmonic form):')")
	write(iout,"(' Q_2,0 =',f12.4,'   Q_2,-1=',f12.4,'   Q_2,1=',f12.4)") R20,R2n1,R2p1
	write(iout,"(' Q_2,-2=',f12.4,'   Q_2,2 =',f12.4)") R2n2,R2p2
	write(iout,"(' Magnitude: |Q_2|=',f12.4)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)
	write(iout,"(a,f16.6)") " Molecular electronic spatial extent <r^2>:",ESE
	write(iout,"(' Components of <r^2>:  X=',f15.6,'  Y=',f15.6,'  Z=',f15.6)") ESEx,ESEy,ESEz
	R30=(5*zzzinttot-3*rrzinttot)/2D0
	R3n1=dsqrt(3D0/8D0)*(5*yzzinttot-rryinttot)
	R3p1=dsqrt(3D0/8D0)*(5*xzzinttot-rrxinttot)
	R3n2=dsqrt(15D0)*xyzinttot
	R3p2=dsqrt(15D0)*(xxzinttot-yyzinttot)/2D0
	R3n3=dsqrt(5D0/8D0)*(3*xxyinttot-yyyinttot)
	R3p3=dsqrt(5D0/8D0)*(xxxinttot-3*xyyinttot)
	write(iout,"(' Molecular octopole moments (Cartesian form):')")
    fac=1
    !fac=au2debye*b2a*b2a !If using this factor, result will be identical to "Octapole moment (field-independent basis, Debye-Ang**2):" printed by Gaussian
	write(iout,"(' XXX=',f10.4,'  YYY=',f10.4,'  ZZZ=',f10.4,'  XYY=',f10.4,'  XXY=',f10.4)") &
    xxxinttot*fac,yyyinttot*fac,zzzinttot*fac,xyyinttot*fac,xxyinttot*fac
	write(iout,"(' XXZ=',f10.4,'  XZZ=',f10.4,'  YZZ=',f10.4,'  YYZ=',f10.4,'  XYZ=',f10.4)") &
    xxzinttot*fac,xzzinttot*fac,yzzinttot*fac,yyzinttot*fac,xyzinttot*fac
	write(iout,"(' Molecular octopole moments (Spherical harmonic form):')")
	write(iout,"(' Q_3,0 =',f11.4,'  Q_3,-1=',f11.4,'  Q_3,1 =',f11.4)") R30,R3n1,R3p1
	write(iout,"(' Q_3,-2=',f11.4,'  Q_3,2 =',f11.4,'  Q_3,-3=',f11.4,'  Q_3,3 =',f11.4)") R3n2,R3p2,R3n3,R3p3
	write(iout,"(' Magnitude: |Q_3|=',f12.4)") dsqrt(R30**2+R3n1**2+R3p1**2+R3n2**2+R3p2**2+R3n3**2+R3p3**2)
    if (iout==20) then
		close(20)
		close(21)
		write(*,"(a)") " Electric dipole and quadrupole moments of all atoms have been exported to multipole.txt and atom_moment.txt in current folder"
    end if
    
else if (isel==3) then !Output AOM
    write(*,*) "Exporting AOM.txt in current folder..."
	open(10,file="AOM.txt",status="replace")
	if (wfntype==0.or.wfntype==2.or.wfntype==3) then
		do iatm=1,ncenter
			write(10,"('Atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
			call showmatgau(AOM(:,:,iatm),"",1,"f14.8",10)
			write(10,*)
		end do
		write(10,"(a)") "Sum of atomic overlap matrices"
		call showmatgau(AOMsum,"",1,"f14.8",10)
	else if (wfntype==1.or.wfntype==4) then
		do iatm=1,ncenter
			write(10,"('Alpha part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
			call showmatgau(AOM(:,:,iatm),"",1,"f14.8",10)
			if (nmatsizeb>0) then
				write(10,"('Beta part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
				call showmatgau(AOMb(:,:,iatm),"",1,"f14.8",10)
			end if
			write(10,*)
		end do
		write(10,"(a)") "Sum of alpha part of atomic overlap matrices"
		call showmatgau(AOMsum,"",1,"f14.8",10)
		write(10,"(a)") "Sum of beta part of atomic overlap matrices"
		call showmatgau(AOMsumb,"",1,"f14.8",10)
	end if
	close(10)
	write(*,*) "Done, atomic overlap matrices have been exported to AOM.txt in current folder"
    
else if (isel==33.and.ifPBC==0) then !Construct and output FOM. For PBC case, FOM has been exported earlier in subroutine AOMFOM_evengrid
	!Generate FOMs
	allocate(FOM1(nmatsize,nmatsize),FOM2(nmatsize,nmatsize))
    FOM1=0
    FOM2=0
    do idx=1,nFOM1atm
		FOM1(:,:)=FOM1(:,:)+AOM(:,:,FOM1atm(idx))
    end do
    if (iFOMmode>1) then
		do idx=1,nFOM2atm
			FOM2(:,:)=FOM2(:,:)+AOM(:,:,FOM2atm(idx))
		end do
    end if
    if (wfntype==1.or.wfntype==2.or.wfntype==4) then
		allocate(FOM1b(nmatsizeb,nmatsizeb),FOM2b(nmatsizeb,nmatsizeb))
		FOM1b=0
		FOM2b=0
		do idx=1,nFOM1atm
			FOM1b(:,:)=FOM1b(:,:)+AOMb(:,:,FOM1atm(idx))
		end do
		if (iFOMmode>1) then
			do idx=1,nFOM2atm
				FOM2b(:,:)=FOM2b(:,:)+AOMb(:,:,FOM2atm(idx))
			end do
        end if
    end if
    !Export FOMs
	open(10,file="FOM.txt",status="replace")
	if (wfntype==0.or.wfntype==3) then !Closed-shell
		if (iFOMmode==1) then
			write(10,"('Fragment overlap matrix')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(*,*) "FOM has been exported to FOM.txt in current folder"
		else
			write(10,"('Fragment overlap matrix of fragment 1')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Fragment overlap matrix of fragment 2')")
			call showmatgau(FOM2(:,:),"",1,"f14.8",10)
			write(*,*) "FOM of fragments 1 and 2 has been exported to FOM.txt in current folder"
		end if
        deallocate(FOM1,FOM2)
	else !Open-shell
		if (iFOMmode==1) then
			write(10,"('Alpha part of fragment overlap matrix')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(*,*)
			write(10,"('Beta part of fragment overlap matrix')")
			call showmatgau(FOM1b(:,:),"",1,"f14.8",10)
			write(*,*) "FOM has been exported to FOM.txt in current folder"
		else
			write(10,"('Alpha part of fragment overlap matrix of fragment 1')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Beta part of fragment overlap matrix of fragment 1')")
			call showmatgau(FOM1b(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Alpha part of fragment overlap matrix of fragment 2')")
			call showmatgau(FOM2(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Beta part of fragment overlap matrix of fragment 2')")
			call showmatgau(FOM2b(:,:),"",1,"f14.8",10)
			write(*,*) "FOM of fragments 1 and 2 has been exported to FOM.txt in current folder"
		end if
        deallocate(FOM1,FOM1b,FOM2,FOM2b)
	end if
	close(10)
	
else if (isel==4) then !Show LI and DI or fuzzy bond order
	if (iwork==0) then !Output LI and DI
		write(*,"(a)") " Note: Delocalization index in fuzzy atomic space is also known as fuzzy bond order"
		!The strict definition of atomic valence in fuzzy space is Eq.18 in CPL,368,375, however in closed-shell case free valence is zero, so sum of bond order is just atomic valence
		write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements, for closed-shell cases, they are also known as atomic valence"
		write(*,*)
		selectyn='n'
		ioutid=6
		do while(.true.)
			if (wfntype==1.or.wfntype==2.or.wfntype==4) then !UHF,ROHF,U-post-HF, output each spin component first
				!Alpha
				call showmatgau(DIa,"Delocalization index matrix for alpha spin",0,"f14.8",ioutid)
				write(ioutid,*)
				if (iwork/=1) then
					write(ioutid,*) "Localization index for alpha spin:"
					do iatm=1,ncenter
						write(ioutid,"(i5,'(',a,'):',f7.3)",advance='no') iatm,ind2name(a(iatm)%index),LIa(iatm)
						if (mod(iatm,4)==0) write(ioutid,*)
					end do
					write(ioutid,*)
					write(ioutid,*)
				end if
				!Beta
				call showmatgau(DIb,"Delocalization index matrix for beta spin",0,"f14.8",ioutid)
				write(ioutid,*)
				if (iwork/=1) then
					write(ioutid,*) "Localization index for beta spin:"
					do iatm=1,ncenter
						write(ioutid,"(i5,'(',a,'):',f7.3)",advance='no') iatm,ind2name(a(iatm)%index),LIb(iatm)
						if (mod(iatm,4)==0) write(ioutid,*)
					end do
					write(ioutid,*)
					write(ioutid,*)
				end if
			end if
			!Alpha+Beta
			call showmatgau(DI,"Total delocalization index matrix",0,"f14.8",ioutid)
			write(ioutid,*)
			if (iwork/=1) then
				write(ioutid,*) "Localization index:"
				do iatm=1,ncenter
					write(ioutid,"(i5,'(',a,'):',f7.3)",advance='no') iatm,ind2name(a(iatm)%index),LI(iatm)
					if (mod(iatm,4)==0) write(ioutid,*)
				end do
				write(ioutid,*)
				write(ioutid,*)
			end if
			if (selectyn=='n') then !Just output result to screen, choose if output result to plain text file
				write(*,*) "If also outputting LI and DI to LIDI.txt in current folder? (y/n)"
				read(*,*) selectyn
				if (selectyn=='y') then
					open(10,file="LIDI.txt",status="replace")
					ioutid=10
				else if (selectyn=='n') then
					exit
				end if
			else if (selectyn=='y') then !Have already outputted result to LIDI.txt, exit cycle
				write(*,*) "Done, the LI and DI have been outputted to LIDI.txt in current folder"
				close(10)
				exit
			end if
		end do
		
	else if (iwork==1) then !Output fuzzy bond order
		write(*,"(' The total bond order >=',f10.6)") bndordthres
		itmp=0
		if (wfntype==1.or.wfntype==2.or.wfntype==4) then
			do i=1,ncenter
				do j=i+1,ncenter
					if (DIa(i,j)+DIb(i,j)>=bndordthres) then
						itmp=itmp+1
						write(*,"(' #',i5,':',i5,a,i5,a,' Alpha: ',f10.6,' Beta:',f10.6,' Total:',f10.6)") &
						itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',DIa(i,j),DIb(i,j),DIa(i,j)+DIb(i,j)
					end if
				end do
			end do
		else if (wfntype==0.or.wfntype==3) then
			itmp=0
			do i=1,ncenter
				do j=i+1,ncenter
					if (DI(i,j)>=bndordthres) then
						itmp=itmp+1
						write(*,"('#',i5,':',5x,i5,a,i5,a,f14.8)") itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',DI(i,j)
					end if
				end do
			end do
		end if
		if (allocated(frag1)) then
			bndordfraga=0
			bndordfragb=0
			bndordfragtot=0
			do i=1,size(frag1)
				do j=1,size(frag2)
					if (wfntype==1.or.wfntype==2.or.wfntype==4) then 
						bndordfraga=bndordfraga+DIa(frag1(i),frag2(j))
						bndordfragb=bndordfragb+DIb(frag1(i),frag2(j))
					else if (wfntype==0.or.wfntype==3) then
						bndordfragtot=bndordfragtot+DI(frag1(i),frag2(j))
					end if
				end do
			end do
			write(*,*)
			if (wfntype==1.or.wfntype==2.or.wfntype==4) then
				write(*,"(' The bond order between fragment 1 and 2:')")
				write(*,"(' Alpha:',f10.6,' Beta:',f10.6,' Total:',f10.6)") bndordfraga,bndordfragb,bndordfraga+bndordfragb
			else if (wfntype==0.or.wfntype==3) then
				write(*,"(' The bond order between fragment 1 and 2:',f12.6)") bndordfragtot
			end if
		end if
		write(*,*)
		write(*,*) "If outputting bond order matrix to bndmat.txt in current folder? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			open(10,file="bndmat.txt",status="replace")
			if (wfntype==1.or.wfntype==2.or.wfntype==4) then !UHF,ROHF,U-post-HF, output each spin component first
				call showmatgau(DIa,"Delocalization index matrix for alpha spin",0,"f14.8",10)
				write(10,*)
				call showmatgau(DIb,"Delocalization index matrix for beta spin",0,"f14.8",10)
				write(10,*)
			end if
			call showmatgau(DI,"Total delocalization index matrix",0,"f14.8",10)
			write(10,*)
			close(10)
			write(*,*) "Done, bond order matrix has been outputted to bndmat.txt in current folder"
			write(*,"(a)") " Note: Diagonal terms in the bond order matrix are the sum of corresponding row or column elements, for closed-shell cases, they are also known as atomic valence"
		end if
		radpot=nradpotold
		sphpot=nsphpotold
		return !Fuzzy bond order has been shown, now (normally) return to bond order analysis interface
	end if
    
else if (isel==44) then !FLI and IFDI
	write(*,*) "Calculating FLI and IFDI..."
    write(*,*)
	!Calculate IFDI
	sumDI=0
	sumDIa=0
    sumDIb=0
    do idx=1,nDIfrag1
		iatm=DIfrag1(idx)
		do jdx=1,nDIfrag2
			jatm=DIfrag2(jdx)
            if (wfntype==0.or.wfntype==3) then
				sumDI=sumDI+DI(iatm,jatm)
            else
				sumDIa=sumDIa+DIa(iatm,jatm)
				sumDIb=sumDIb+DIb(iatm,jatm)
            end if
        end do
    end do
    !Calculate FLI
    if (wfntype==0.or.wfntype==3) then
		sumLI1=0
		do idx=1,nDIfrag1
			sumLI1=sumLI1+LI(DIfrag1(idx))
			do jdx=idx+1,nDIfrag1
				sumLI1=sumLI1+DI(DIfrag1(idx),DIfrag1(jdx))
			end do
		end do
		sumLI2=0
		do idx=1,nDIfrag2
			sumLI2=sumLI2+LI(DIfrag2(idx))
			do jdx=idx+1,nDIfrag2
				sumLI2=sumLI2+DI(DIfrag2(idx),DIfrag2(jdx))
			end do
		end do
    else
		sumLI1a=0;sumLI1b=0
		do idx=1,nDIfrag1
			sumLI1a=sumLI1a+LIa(DIfrag1(idx))
			sumLI1b=sumLI1b+LIb(DIfrag1(idx))
			do jdx=idx+1,nDIfrag1
				sumLI1a=sumLI1a+DIa(DIfrag1(idx),DIfrag1(jdx))
				sumLI1b=sumLI1b+DIb(DIfrag1(idx),DIfrag1(jdx))
			end do
		end do
		sumLI2a=0;sumLI2b=0
		do idx=1,nDIfrag2
			sumLI2a=sumLI2a+LIa(DIfrag2(idx))
			sumLI2b=sumLI2b+LIb(DIfrag2(idx))
			do jdx=idx+1,nDIfrag2
				sumLI2a=sumLI2a+DIa(DIfrag2(idx),DIfrag2(jdx))
				sumLI2b=sumLI2b+DIb(DIfrag2(idx),DIfrag2(jdx))
			end do
		end do
    end if
    if (wfntype==0.or.wfntype==3) then
        write(*,"(' LI of fragment 1:',f14.6)") sumLI1
        write(*,"(' LI of fragment 2:',f14.6)") sumLI2
		write(*,"(' Interfragment DI:',f14.6)") sumDI
    else
		write(*,*) "Fragment 1:"
        write(*,"(' Fragment LI of alpha spin:',f14.6)") sumLI1a
        write(*,"(' Fragment LI of beta spin: ',f14.6)") sumLI1b
        write(*,"(' Total fragment LI:        ',f14.6)") sumLI1a+sumLI1b
        write(*,*)
		write(*,*) "Fragment 2:"
        write(*,"(' Fragment LI of alpha spin:',f14.6)") sumLI2a
        write(*,"(' Fragment LI of beta spin: ',f14.6)") sumLI2b
        write(*,"(' Total fragment LI:        ',f14.6)") sumLI2a+sumLI2b
        write(*,*)
		write(*,"(' Interfragment DI of alpha spin:',f14.6)") sumDIa
		write(*,"(' Interfragment DI of beta spin: ',f14.6)") sumDIb
		write(*,"(' Total interfragment DI:        ',f14.6)") sumDIa+sumDIb
    end if
	
else if (isel==5) then !PDI
	call showmatgau(DI,"Delocalization index matrix",0,"f14.8")
	write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements, for closed-shell cases, these also known as atomic valence"
	do while(.true.)
		write(*,"(/,a)") " Input indices of the six atoms constituting the ring, in clockwise or anti-clockwise order. e.g. 4,5,6,7,8,2"
		write(*,*) "(Input q can return)"
		read(*,"(a)") c80inp
		if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
		read(c80inp,*) PDIatom(:)
		write(*,"(' Delocalization index of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PDIatom(1),a(PDIatom(1))%name,PDIatom(4),a(PDIatom(4))%name,DI(PDIatom(1),PDIatom(4))
		write(*,"(' Delocalization index of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PDIatom(2),a(PDIatom(2))%name,PDIatom(5),a(PDIatom(5))%name,DI(PDIatom(2),PDIatom(5))
		write(*,"(' Delocalization index of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PDIatom(3),a(PDIatom(3))%name,PDIatom(6),a(PDIatom(6))%name,DI(PDIatom(3),PDIatom(6))
		write(*,"(' PDI value is',f12.6)") ( DI(PDIatom(1),PDIatom(4))+DI(PDIatom(2),PDIatom(5))+DI(PDIatom(3),PDIatom(6)) )/3D0
	end do
	
else if (isel==6) then !FLU
	call showmatgau(DI,"Delocalization index matrix",0,"f14.8")
	write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements, for closed-shell cases, these also known as atomic valence"
	write(*,*)
	write(*,*) "Current FLU reference parameters:"
	do iref=1,nelesupp
		do jref=iref,nelesupp
			if (FLUref(iref,jref)/=-1) write(*,"(' ',a,a,a,a,f10.5)") ind2name(iref),'-',ind2name(jref),':',FLUref(iref,jref)
		end do
	end do
	do while(.true.)
		write(*,"(/,a)") " Input indices of the atoms in the ring, in clockwise or anti-clockwise order"
		write(*,*) "e.g. 4,7,8,1,2,3      (Input q can exit)"
		read(*,"(a)") c80inp
		if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
		call str2arr(c80inp,nFLUatom,FLUatom)
		FLUval=0D0
        write(*,*)
		write(*,*) "        Atom pair         Contribution          DI"
		do iidx=1,nFLUatom
			jidx=iidx+1
			if (iidx==nFLUatom) jidx=1 !Return to the first element of the ring after a cycle
			iatm=FLUatom(iidx) !Actual atom index in present system
			jatm=FLUatom(jidx)
			iatmeleidx=a(iatm)%index !Index in periodic table
			jatmeleidx=a(jatm)%index
			refval=FLUref(iatmeleidx,jatmeleidx)
			if (refval==-1D0) then
				write(*,"(' Error: Missing reference parameter for',a,'-',a)") ind2name(iatmeleidx),ind2name(jatmeleidx)
				exit
			end if
			valenratio=DI(iatm,iatm)/DI(jatm,jatm) !DI(iatm,iatm) is the sum of corresponding row or column elements, namely atomic valence, rather than LI*2 of iatm
			if (valenratio<1) valenratio=1D0/valenratio
			FLUpair=(valenratio*( (DI(iatm,jatm)-refval)/refval ))**2/nFLUatom
			write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),FLUpair,DI(iatm,jatm)
			FLUval=FLUval+FLUpair
		end do
		write(*,"(' FLU value is',f12.6)") FLUval
	end do
	
else if (isel==7) then !FLU-pi
	write(*,*) "Which occupied orbitals are pi orbitals? Input their indices, e.g. 17,20,21"
	read(*,"(a)") c80inp
	call str2arr(c80inp,nFLUorb,FLUorb)
	!Generate DI for pi orbitals. DI_A,B=2∑[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B     where i and j are non-spin orbital
	DI=0D0
	do iatm=1,ncenter
		do jatm=iatm+1,ncenter
			tmpval=0D0
			do iidx=1,nFLUorb
				iorb=FLUorb(iidx)
				do jidx=1,nFLUorb
					jorb=FLUorb(jidx)
					tmpval=tmpval+dsqrt(MOocc(iorb)*MOocc(jorb))*AOM(iorb,jorb,iatm)*AOM(iorb,jorb,jatm)
				end do
			end do
			DI(iatm,jatm)=tmpval
		end do
	end do
	DI=2*(DI+transpose(DI))
	do iatm=1,ncenter !Calculate atomic valence
		DI(iatm,iatm)=sum(DI(iatm,:))
	end do
	call showmatgau(DI,"Delocalization index matrix for pi electrons",0,"f14.8")
	write(*,"(a)") " Note: Diagonal terms are the sum of corresponding row or column elements"
	do while(.true.)
		write(*,"(/,a)") " Input indices of the atoms in the ring, in clockwise or anti-clockwise"
		write(*,*) "e.g. 4,7,8,1,2,3      (Input q can exit)"
		read(*,"(a)") c80inp
		if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
		call str2arr(c80inp,nFLUatom,FLUatom)
		!Calculate average of DI-oi first
		avgDI=0D0
		do iidx=1,nFLUatom
			jidx=iidx+1
			if (iidx==nFLUatom) jidx=1
			avgDI=avgDI+DI(FLUatom(iidx),FLUatom(jidx))
		end do
		avgDI=avgDI/nFLUatom
		write(*,"(' Average of DI-pi is',f12.6)") avgDI
		FLUval=0D0
		write(*,*) "        Atom pair         Contribution          DI"
		do iidx=1,nFLUatom
			jidx=iidx+1
			if (iidx==nFLUatom) jidx=1
			iatm=FLUatom(iidx) !Actual atom index in present system
			jatm=FLUatom(jidx)
			iatmeleidx=a(iatm)%index !Index in periodic table
			jatmeleidx=a(jatm)%index
			valenratio=DI(iatm,iatm)/DI(jatm,jatm)
			if (valenratio<1) valenratio=1D0/valenratio
			FLUpair=(valenratio*(DI(iatm,jatm)-avgDI)/avgDI)**2/nFLUatom
			write(*,"(i5,'(',a,')  --',i5,'(',a,'):',f15.6,f16.6)") iatm,ind2name(iatmeleidx),jatm,ind2name(jatmeleidx),FLUpair,DI(iatm,jatm)
			FLUval=FLUval+FLUpair
		end do
		write(*,"(' FLU-pi value is',f12.6)") FLUval
	end do
	
else if (isel==8) then !Integral in overlap region
	ovlpintpos=ovlpintpos+transpose(ovlpintpos)
	ovlpintneg=ovlpintneg+transpose(ovlpintneg)
	sumdiagpos=0D0
	sumdiagneg=0D0
	do i=1,ncenter
		ovlpintpos(i,i)=ovlpintpos(i,i)/2D0
		sumdiagpos=sumdiagpos+ovlpintpos(i,i)
		ovlpintneg(i,i)=ovlpintneg(i,i)/2D0
		sumdiagneg=sumdiagneg+ovlpintneg(i,i)
	end do
	ovlpinttot=ovlpintpos+ovlpintneg
	if (iwork==2) then !Output Laplacian bond order
		write(*,"(' The bond orders >=',f10.6)") bndordthres
		itmp=0
		do i=1,ncenter
			do j=i+1,ncenter
				if (-10*ovlpintneg(i,j)>=bndordthres) then
					itmp=itmp+1
					write(*,"(' #',i5,':',i5,a,i5,a,':',f10.6)") &
					itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',-10*ovlpintneg(i,j)
				end if
			end do
		end do
		if (allocated(frag1)) then !Output interfragment bond order
			bndordfragtot=0
			do i=1,size(frag1)
				do j=1,size(frag2)
					bndordfragtot=bndordfragtot-10*ovlpintneg(frag1(i),frag2(j))
				end do
			end do
			write(*,*)
			write(*,"(' The bond order between fragment 1 and 2:',f12.6)") bndordfragtot
		end if
		write(*,*)
		write(*,*) "If outputting bond order matrix to bndmat.txt in current folder? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			open(10,file="bndmat.txt",status="replace")
			do i=1,ncenter
				ovlpintneg(i,i)=sum(ovlpintneg(i,:))-ovlpintneg(i,i) !Make diagonal terms are the sum of corresponding row elements, namely valence
			end do
			call showmatgau(-10*ovlpintneg,"Laplacian bond order matrix",0,"f14.8",10)
			close(10)
			write(*,*) "Done, bond order matrix has been outputted to bndmat.txt in current folder"
			write(*,"(a)") " Note: Diagonal terms in the bond order matrix are the sum of corresponding row or column elements"
			write(*,*)
		end if
		radpot=nradpotold
		sphpot=nsphpotold
		return !Laplacian bond order has been shown, now (normally) return to bond order analysis interface
	else
		call showmatgau(ovlpintpos,"Integration of positive values in overlap region",0,"f14.8")
		sumovlppos=sum(ovlpintpos)
		write(*,"(' Summing up diagonal matrix elements:     ',f20.8)") sumdiagpos
		write(*,"(' Summing up non-diagonal, matrix elements:',f20.8)") sumovlppos-sumdiagpos
		write(*,"(' Summing up all matrix elements:          ',f20.8)") sumovlppos
		write(*,*)
		sumovlpneg=sum(ovlpintneg)
		call showmatgau(ovlpintneg,"Integration of negative values in overlap region",0,"f14.8")
		write(*,"(' Summing up diagonal matrix elements:     ',f20.8)") sumdiagneg
		write(*,"(' Summing up non-diagonal, matrix elements:',f20.8)") sumovlpneg-sumdiagneg
		write(*,"(' Summing up all matrix elements:          ',f20.8)") sumovlpneg
		write(*,*)
		sumovlptot=sum(ovlpinttot)
		call showmatgau(ovlpinttot,"Integration of all values in overlap region",0,"f14.8")
		write(*,"(' Summing up diagonal matrix elements:     ',f20.8)") sumdiagpos+sumdiagneg
		write(*,"(' Summing up non-diagonal, matrix elements:',f20.8)") sumovlptot-sumdiagpos-sumdiagneg
		write(*,"(' Summing up all matrix elements:          ',f20.8)") sumovlptot
		write(*,*)
		write(*,*) "If also outputting above matrices to intovlp.txt in current folder? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			open(10,file="intovlp.txt",status="replace")
			call showmatgau(ovlpintpos,"Integration of positive values in overlap region",0,"f14.8",10)
			write(10,*)
			call showmatgau(ovlpintneg,"Integration of negative values in overlap region",0,"f14.8",10)
			write(10,*)
			call showmatgau(ovlpinttot,"Integration of all values in overlap region",0,"f14.8",10)
			write(10,*)
			close(10)
			write(*,*) "Done, the matrices have been outputted to intovlp.txt in current folder"
		end if
	end if
	
else if (isel==9) then !CLRK
	call showmatgau(CLRK,"Condensed linear response kernel (CLRK) matrix",0,"f14.8")
	write(*,*)
	write(*,*) "If also outputting CLRK to CLRK.txt in current folder? (y/n)"
	read(*,*) selectyn
	if (selectyn=='y'.or.selectyn=='Y') then
		open(10,file="CLRK.txt",status="replace")
		call showmatgau(CLRK,"Condensed linear response kernel (CLRK) matrix",0,"f14.8",10)
		close(10)
		write(*,*) "Done, the CLRK matrix has been outputted to CLRK.txt in current folder"
	end if
	
else if (isel==10) then !PLR
	call showmatgau(CLRK,"Condensed linear response kernel (CLRK) matrix",0,"f14.8")
	write(*,*)
	do while(.true.)
		write(*,"(/,34a)") " Input indices of the six atoms constituting the ring, in clockwise or anti-clockwise. e.g. 4,5,6,7,8,2"
		write(*,*) "(Input q can return)"
		read(*,"(a)") c80inp
		if (c80inp(1:1)=='q'.or.c80inp(1:1)=='Q') exit
		read(c80inp,*) PLRatom(:)
		write(*,"(' CLRK of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PLRatom(1),a(PLRatom(1))%name,PLRatom(4),a(PLRatom(4))%name,CLRK(PLRatom(1),PLRatom(4))
		write(*,"(' CLRK of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PLRatom(2),a(PLRatom(2))%name,PLRatom(5),a(PLRatom(5))%name,CLRK(PLRatom(2),PLRatom(5))
		write(*,"(' CLRK of ',i5,'(',a,')   --',i5,'(',a,'):',f12.6)") PLRatom(3),a(PLRatom(3))%name,PLRatom(6),a(PLRatom(6))%name,CLRK(PLRatom(3),PLRatom(6))
		write(*,"(' PLR index is',f12.6)") ( CLRK(PLRatom(1),PLRatom(4))+CLRK(PLRatom(2),PLRatom(5))+CLRK(PLRatom(3),PLRatom(6)) )/3D0
	end do
	
else if (isel==11) then !Multicenter DI
	do while(.true.)
		write(*,*) "Input atom indices, e.g. 3,4,7,8,10    (Up to 10 atoms)"
		write(*,*) "Input q can return to upper level menu"
		read(*,"(a)") c80inp
		if (c80inp(1:1)=='q') then
			exit
		else
			call str2arr(c80inp,nDIcen,cenind)
		end if
		DImulti=0D0
		write(*,*) "Please wait..."
		if (nDIcen==3) then
			do iorb=1,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
			DImulti=DImulti+&
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,iorb,cenind(3))
					end do
				end do
			end do
		else if (nDIcen==4) then
			do iorb=1,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
						do lorb=1,nmatsize
							if (MOocc(lorb)==0D0) cycle
			DImulti=DImulti+&
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,iorb,cenind(4))
						end do
					end do
				end do
			end do
		else if (nDIcen==5) then
			do iorb=1,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
						do lorb=1,nmatsize
							if (MOocc(lorb)==0D0) cycle
							do morb=1,nmatsize
								if (MOocc(morb)==0D0) cycle
			DImulti=DImulti+&
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,iorb,cenind(5))
							end do
						end do
					end do
				end do
			end do				
		else if (nDIcen==6) then
			do iorb=1,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
						do lorb=1,nmatsize
							if (MOocc(lorb)==0D0) cycle
							do morb=1,nmatsize
								if (MOocc(morb)==0D0) cycle
								do norb=1,nmatsize
									if (MOocc(norb)==0D0) cycle
			DImulti=DImulti+& !dsqrt(MOocc(iorb)*MOocc(jorb)*MOocc(korb)*MOocc(lorb)*MOocc(morb)*MOocc(norb))*
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*AOM(norb,iorb,cenind(6))
								end do
							end do
						end do
					end do
				end do
			end do
		else if (nDIcen==7) then
			do iorb=1,nmatsize
				write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
						do lorb=1,nmatsize
							if (MOocc(lorb)==0D0) cycle
							do morb=1,nmatsize
								if (MOocc(morb)==0D0) cycle
								do norb=1,nmatsize
									if (MOocc(norb)==0D0) cycle
									do iiorb=1,nmatsize
										if (MOocc(iiorb)==0D0) cycle
			DImulti=DImulti+&
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
			AOM(norb,iiorb,cenind(6))*AOM(iiorb,iorb,cenind(7))
									end do
								end do
							end do
						end do
					end do
				end do
			end do
		else if (nDIcen==8) then
			do iorb=1,nmatsize
				write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
						do lorb=1,nmatsize
							if (MOocc(lorb)==0D0) cycle
							do morb=1,nmatsize
								if (MOocc(morb)==0D0) cycle
								do norb=1,nmatsize
									if (MOocc(norb)==0D0) cycle
									do iiorb=1,nmatsize
										if (MOocc(iiorb)==0D0) cycle
										do jjorb=1,nmatsize
											if (MOocc(jjorb)==0D0) cycle
			DImulti=DImulti+&
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
			AOM(norb,iiorb,cenind(6))*AOM(iiorb,jjorb,cenind(7))*AOM(jjorb,iorb,cenind(8))
										end do
									end do
								end do
							end do
						end do
					end do
				end do
			end do
		else if (nDIcen==9) then
			do iorb=1,nmatsize
				write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
						do lorb=1,nmatsize
							if (MOocc(lorb)==0D0) cycle
							do morb=1,nmatsize
								if (MOocc(morb)==0D0) cycle
								do norb=1,nmatsize
									if (MOocc(norb)==0D0) cycle
									do iiorb=1,nmatsize
										if (MOocc(iiorb)==0D0) cycle
										do jjorb=1,nmatsize
											if (MOocc(jjorb)==0D0) cycle
											do kkorb=1,nmatsize
												if (MOocc(kkorb)==0D0) cycle
			DImulti=DImulti+&
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
			AOM(norb,iiorb,cenind(6))*AOM(iiorb,jjorb,cenind(7))*AOM(jjorb,kkorb,cenind(8))*AOM(kkorb,iorb,cenind(9))
											end do
										end do
									end do
								end do
							end do
						end do
					end do
				end do
			end do
		else if (nDIcen==10) then
			do iorb=1,nmatsize
				write(*,"(' Finished',i8,'/',i8)") iorb,nmatsize
				if (MOocc(iorb)==0D0) cycle
				do jorb=1,nmatsize
					if (MOocc(jorb)==0D0) cycle
					do korb=1,nmatsize
						if (MOocc(korb)==0D0) cycle
						do lorb=1,nmatsize
							if (MOocc(lorb)==0D0) cycle
							do morb=1,nmatsize
								if (MOocc(morb)==0D0) cycle
								do norb=1,nmatsize
									if (MOocc(norb)==0D0) cycle
									do iiorb=1,nmatsize
										if (MOocc(iiorb)==0D0) cycle
										do jjorb=1,nmatsize
											if (MOocc(jjorb)==0D0) cycle
											do kkorb=1,nmatsize
												if (MOocc(kkorb)==0D0) cycle
												do llorb=1,nmatsize
													if (MOocc(llorb)==0D0) cycle
			DImulti=DImulti+&
			AOM(iorb,jorb,cenind(1))*AOM(jorb,korb,cenind(2))*AOM(korb,lorb,cenind(3))*AOM(lorb,morb,cenind(4))*AOM(morb,norb,cenind(5))*&
			AOM(norb,iiorb,cenind(6))*AOM(iiorb,jjorb,cenind(7))*AOM(jjorb,kkorb,cenind(8))*AOM(kkorb,llorb,cenind(9))*AOM(llorb,iorb,cenind(10))
												end do
											end do
										end do
									end do
								end do
							end do
						end do
					end do
				end do
			end do
		end if
		DImulti=DImulti*2**(nDIcen-1)
		write(*,"(' Multicenter DI:',f13.7,/)") DImulti
		write(*,"(' Multicenter DI in normalized form: ',f13.7,/)") DImulti**(1D0/nDIcen)
	end do
	
else if (isel==12) then !Information-theoretic defined aromaticity
	valavg=0
	do idx=1,naromatatm
		iatm=aromatatm(idx)
		tmpval=rintval(iatm,1)
		write(*,"(' Atom',i5,'(',a,'):',f12.6)") iatm,a(iatm)%name,tmpval
		valavg=valavg+tmpval
	end do
	write(*,"(' The result (average of above data) is',f12.6)") valavg/naromatatm
    
else if (isel==13) then !Atomic polarizability
	write(*,*) "Atomic polarizabilities estimated using Tkatchenko-Scheffler method:"
	totpol=0
	do iatm=1,ncenter
		if (all(atmcalclist(1:natmcalclist)/=iatm)) cycle
        totpol=totpol+atmpol(iatm)
    end do
	do iatm=1,ncenter
		if (all(atmcalclist(1:natmcalclist)/=iatm)) cycle
        write(*,"(i5,'(',a,'):',f8.3,' a.u.  Contribution:',f6.2,' %  (Ref. data:',f8.3,' a.u.)')") &
        iatm,a(iatm)%name,atmpol(iatm),atmpol(iatm)/totpol*100,atmpol_free(a(iatm)%index)
    end do
    write(*,"(' Sum of atomic polarizabilities:',f10.3,' a.u.')") totpol
	write(*,*)
    write(*,*) "Atomic C6 coefficients estimated using Tkatchenko-Scheffler method:"
	do iatm=1,ncenter
		if (all(atmcalclist(1:natmcalclist)/=iatm)) cycle
        write(*,"(i5,'(',a,'):',f8.2,' a.u. (Ref. data:',f8.1,' a.u.)')") &
        iatm,a(iatm)%name,atmC6(iatm),atmC6_free(a(iatm)%index)
    end do
    write(*,*)
    write(*,*) "Note: Reference data denotes the built-in value of free-state atom"
    C6mol=0
    do iatm=1,ncenter
		C6mol=C6mol+atmC6(iatm)
		do jatm=iatm+1,ncenter
			C6AB=2*atmC6(iatm)*atmC6(jatm)/(atmpol(jatm)/atmpol(iatm)*atmC6(iatm) + atmpol(iatm)/atmpol(jatm)*atmC6(jatm))
			!write(*,"(i5,'(',a,') --',i5,'(',a,'):',f8.2,' a.u.')") iatm,a(iatm)%name,jatm,a(jatm)%name,C6AB
            C6mol=C6mol+2*C6AB
        end do
    end do
    write(*,"(a,f10.2,' a.u.')") " Homomolecular C6 coefficient:",C6mol
end if
	
end do !End interface loop

end subroutine





!!--------- Integrate atomic spaces for a function using evenly distributed grids
!mainly for periodic wavefunction only representing valence electrons
!ipartition: =3 Hirshfeld using built-in atomic density, =4 Hirshfeld-I, =5 MBIS
!ifunc: Index of the real space function to be integrated
subroutine intatmspace_evengrid(ipartition,ifunc,atmint)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8 atmrho(ncenter),tvec(3),atmint(ncenter),atmint_tmp(ncenter)

do iatm=1,ncenter
	if (a(iatm)%index>4.and.a(iatm)%index==nint(a(iatm)%charge)) then !MOLOPT is all-electron basis set for first <=Be
		write(*,"(a)") " Warning: This function employs evenly distributed grids for integration, it does not work well for all-electron wavefunction!"
		write(*,*) "Press ENTER button to continue"
		read(*,*)
		exit
    end if
end do

call setgrid_for_PBC(0.2D0,1)
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
call walltime(iwalltime1)
!Because uniform grid cannot integrate well core density, so temporarily disable EDFs
nEDFprims_org=nEDFprims
nEDFprims=0
call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
write(*,*) "Calculating grid data of selected real space function..."
call savecubmat(ifunc,0,1)
call delvirorb_back(1) !Restore to previous wavefunction
nEDFprims=nEDFprims_org
call calc_dvol(dvol)

write(*,*)
write(*,*) "Calculating atomic contributions..."
atmint(:)=0
ifinish=0;ishowprog=1
ntmp=floor(ny*nz/100D0)
!$OMP PARALLEL SHARED(atmint,ifinish,ishowprog) PRIVATE(i,j,k,tmpx,tmpy,tmpz,iatm,atmrho,prorho,atmint_tmp,&
!$OMP icell,jcell,kcell,tvec,dist2,tmprho,npt) NUM_THREADS(nthreads)
atmint_tmp(:)=0
!$OMP DO schedule(dynamic) collapse(2)
do k=1,nz
	do j=1,ny
		do i=1,nx
			if (cubmat(i,j,k)<1D-12) cycle
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
            atmrho(:)=0
            do icell=-PBCnx,PBCnx
                do jcell=-PBCny,PBCny
                    do kcell=-PBCnz,PBCnz
                        call tvec_PBC(icell,jcell,kcell,tvec)
                        do iatm=1,ncenter
                            dist2=(a(iatm)%x+tvec(1)-tmpx)**2+(a(iatm)%y+tvec(2)-tmpy)**2+(a(iatm)%z+tvec(3)-tmpz)**2
                            if (dist2>atmrhocutsqr(a(iatm)%index)) then
                                cycle
                            else
                                if (ipartition==3) then !Hirshfeld, using bulit-in atomic radial density to interpolate
                                    tmprho=eleraddens(a(iatm)%index,dsqrt(dist2),0)
                                else !Hirshfeld-I and MBIS. Refined atomic radial density of every atom has been available in atmraddens
									npt=atmradnpt(iatm)
									call lagintpol(atmradpos(1:npt),atmraddens(1:npt,iatm),npt,dsqrt(dist2),tmprho,rnouse,rnouse,1)
                                end if
                                atmrho(iatm)=atmrho(iatm)+tmprho
                            end if
                        end do
                    end do
                end do
            end do
            prorho=sum(atmrho(:))
            if (prorho>0) atmint_tmp(:)=atmint_tmp(:)+atmrho(:)/prorho*cubmat(i,j,k)
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
!$OMP END DO
!$OMP CRITICAL
atmint(:)=atmint(:)+atmint_tmp(:)*dvol
!$OMP END CRITICAL
!$OMP END PARALLEL
if (ishowprog/=0) call showprog(100,100)

call walltime(iwalltime2)
write(*,"(/,' Calculation totally took up wall clock time',i10,' s')") iwalltime2-iwalltime1
end subroutine





!!------------- Calculate atomic overlap matrix (AOM) or fragment overlap matrix (FOM) using evenly distributed grids. Fragment LI and interfragment DI can also be outputted
!  Calculating AOM is extremely expensive for large systems consisting of thousands of orbitals. I have try my best to optimize the code, it is however still very expensive
!ipartition: Partition method of atomic spaces
!iwork: 0 means calculating AOM normally; 1 means this subroutine is invoked for calculating fuzzy bond order
!       33 means calculating FOM and then exporting as .txt file
!       44 means calculating AOM only for the atoms involved in fragment LI and interfragment DI, which will be evaluated later by other code
!iendalpha: Index of last alpha orbital
!DIfrag1,DIfrag2: Atom indices of fragments 1 and 2 used to calculate interfragment DI
!nDIfrag1,nDIfrag2: Number of actual atoms in DIfrag1 and DIfrag2
subroutine AOMFOM_evengrid(ipartition,iwork,AOM,AOMb,nmatsize,nmatsizeb,iendalpha,DIfrag1,DIfrag2,nDIfrag1,nDIfrag2)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer ipartition,nmatsize,nmatsizeb,iendalpha,iwork
real*8 AOM(nmatsize,nmatsize,ncenter),AOMb(nmatsizeb,nmatsizeb,ncenter) !,tmpmat(nmatsize,nmatsize)
integer DIfrag1(ncenter),DIfrag2(ncenter),nDIfrag1,nDIfrag2
real*8 atmrho(ncenter),tvec(3),orbval(nmo)
real*8,allocatable :: AOM_tmp(:,:,:),AOMb_tmp(:,:,:)
real*8 :: thres=3D-6 !Threshold for ignoring an atom in loop atoms. This threshold is found to be good balance between accuracy and efficiency
integer FOM1atm(ncenter),FOM2atm(ncenter),nFOM1atm,nFOM2atm
real*8 bndmata(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmattot(ncenter,ncenter)
real*8,allocatable :: FOM1(:,:),FOM1b(:,:),FOM2(:,:),FOM2b(:,:)
logical atmdolist(ncenter)

do iatm=1,ncenter
	if (a(iatm)%index>4.and.a(iatm)%index==nint(a(iatm)%charge)) then !MOLOPT is all-electron basis set for first <=Be
		write(*,"(a)") " Warning: This function employs evenly distributed grids for integration, it does not work well for all-electron wavefunction!"
		write(*,*) "Press ENTER button to continue"
		read(*,*)
		exit
    end if
end do

if (iwork==33) then
	write(*,*)
	write(*,*) "Select the way of calculating fragment overlap matrix (FOM)"
    if (iwork==33) then
		write(*,*) "1 Calculate and export FOM for one fragment"
		write(*,*) "2 Calculate and export FOM for two specific fragments"
		write(*,*) "3 Calculate and export FOM for two fragments. The atoms in the first fragment is directly specified, &
		&the atoms in fragment 2 are those having Mayer bond order with any atom in fragment 1 larger than a specific threshold"
    else
		write(*,*) "2 Calculate FOM for two specific fragments and then calculate DI"
    end if
    read(*,*) iFOMmode
    write(*,*) "Input indices of atoms for defining fragment 1, e.g. 3-9,14,19-20"
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,nFOM1atm,FOM1atm)
    if (iFOMmode==2) then
		write(*,*) "Input indices of atoms for defining fragment 2, e.g. 1,2,10-13,15-18"
		read(*,"(a)") c2000tmp
		call str2arr(c2000tmp,nFOM2atm,FOM2atm)
    else if (iFOMmode==3) then
		write(*,*) "Input threshold of Mayer bond order, e.g. 0.001"
        write(*,*) "Note: If an atom not in fragment 1 has absolute value of total Mayer bond order with any atom &
        &in fragment 1 larger than this threshold, then this atom will be in fragment 2"
        read(*,*) FOM2thres
        call ask_Sbas_PBC
		call calcMayerbndord(bndmattot,bndmata,bndmatb)
        if (wfntype==1.or.wfntype==2.or.wfntype==4) bndmattot=bndmata+bndmatb !Make total bond order in open-shell case correspond to Mayer bond order definition
        nFOM2atm=0
        do iatm=1,ncenter
			if (any(FOM1atm(1:nFOM1atm)==iatm)) cycle
			if (any(abs(bndmattot(FOM1atm(1:nFOM1atm),iatm))>FOM2thres)) then
				nFOM2atm=nFOM2atm+1
                FOM2atm(nFOM2atm)=iatm
                cycle
            end if
        end do
        write(*,*) "Atoms in fragment 2:"
        call arr2str_2(FOM2atm(1:nFOM2atm),c2000tmp)
        write(*,*) trim(c2000tmp)
        sumbo=sum(bndmattot(FOM1atm(1:nFOM1atm),FOM2atm(1:nFOM2atm)))
        write(*,"(' Mayer bond order between fragments 1 and 2:',f8.3)") sumbo
    end if
	atmdolist=.false.
    do idx=1,nFOM1atm
		atmdolist(FOM1atm(idx))=.true.
    end do
    do idx=1,nFOM2atm
		atmdolist(FOM2atm(idx))=.true.
    end do
else if (iwork==44) then !Will calculate AOM only for the atoms involved in evaluating fragment LI and interfragment DI
	atmdolist=.false.
    do idx=1,nDIfrag1
		atmdolist(DIfrag1(idx))=.true.
    end do
    do idx=1,nDIfrag2
		atmdolist(DIfrag2(idx))=.true.
    end do
else !Will calculate full AOM
	atmdolist(:)=.true.
end if

write(*,*)
if (iwork==1) then
	call setgrid_for_PBC(0.35D0,1)
else
	call setgrid_for_PBC(0.2D0,1)
end if
call calc_dvol(dvol)

call walltime(iwalltime1)

call gen_neigh_GTF !Generate neighbouring GTFs list at reduced grids, for faster evaluation

!If each thread allocates AOM_tmp(nmatsize,nmatsize,ncenter), for large system it is very easily out-of-memory, so I decide &
!calculate AOM contributed by different range of atoms in different batches, so that each batch only use AOM_tmp(nmatsize,nmatsize,natm_per_batch). &
!The number of batches is determined dynamically according to present OpenMP stacksize
write(*,*)
write(*,*) "Calculating atomic overlap matrix (AOM)..."
write(*,"(' Number of atoms will be considered in atomic loop:',i7)") count(atmdolist.eqv..true.)
if (wfntype==0.or.wfntype==3) then
	write(*,"(' Number of orbitals will be considered:',i7)") nmatsize
else
	write(*,"(' Number of alpha orbitals will be considered:',i7)") nmatsize
	write(*,"(' Number of bete orbitals will be considered: ',i7)") nmatsizeb
end if
write(*,"(' Threshold of atomic weight in looping atoms:',f12.8)") thres
tmpmem_per_atm=8*(nmatsize*nmatsize+nmatsizeb*nmatsizeb) !The size (Bytes) of the matrix only recording contribution from one atom
write(*,"(' OpenMP stacksize for each thread: ',f10.2,' MB')") dfloat(ompstacksize)/1024/1024
write(*,"(a,f12.1,' MB memory')") " Each OpenMP thread needs at least",tmpmem_per_atm/1024/1024
natm_per_batch=min(floor(ompstacksize/tmpmem_per_atm),ncenter)
nbatch=ceiling(dfloat(ncenter)/natm_per_batch)
if (nbatch>1) write(*,"(a)") " Note: If you have enough physical memory, after properly increasing OpenMP stacksize, number of batches can be reduced, making total computational time lower"
write(*,"(' Maximum number of atoms per batch:',i6)") natm_per_batch
write(*,"(' Number of batches:',i6)") nbatch

AOM(:,:,:)=0
allocate(AOM_tmp(nmatsize,nmatsize,natm_per_batch))
if ((wfntype==1.or.wfntype==4).and.nmatsizeb>0) then
	AOMb(:,:,:)=0
	allocate(AOMb_tmp(nmatsizeb,nmatsizeb,natm_per_batch))
end if

ntmp=floor(ny*nz/100D0)
do ibatch=1,nbatch
    iatmbeg=(ibatch-1)*natm_per_batch+1
    iatmend=min(iatmbeg+natm_per_batch-1,ncenter)
	write(*,"(/,' Performing batch',i4,' /',i4,', from atom',i6,' to',i6)") ibatch,nbatch,iatmbeg,iatmend
    natmthis=iatmend-iatmbeg+1
	ifinish=0;ishowprog=1
	call showprog(0,100)
	!$OMP PARALLEL SHARED(AOM,AOMb,ifinish,ishowprog) PRIVATE(i,j,k,tmpx,tmpy,tmpz,iatm,atmrho,prorho,&
	!$OMP icell,jcell,kcell,tvec,dist2,tmprho,npt,AOM_tmp,AOMb_tmp,imo,jmo,MOinit,MOend,orbval,tmpval,wei,idx,icalcorb) NUM_THREADS(nthreads)
	AOM_tmp(:,:,:)=0
	if ((wfntype==1.or.wfntype==4).and.nmatsizeb>0) AOMb_tmp(:,:,:)=0
	!$OMP DO schedule(dynamic) collapse(2)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
				atmrho(:)=0
				do icell=-PBCnx,PBCnx
					do jcell=-PBCny,PBCny
						do kcell=-PBCnz,PBCnz
							call tvec_PBC(icell,jcell,kcell,tvec)
							do iatm=1,ncenter
								dist2=(a(iatm)%x+tvec(1)-tmpx)**2+(a(iatm)%y+tvec(2)-tmpy)**2+(a(iatm)%z+tvec(3)-tmpz)**2
								if (dist2>atmrhocutsqr(a(iatm)%index)) then
									cycle
								else
									if (ipartition==3) then !Hirshfeld, using built-in atomic radial density to interpolate
										tmprho=eleraddens(a(iatm)%index,dsqrt(dist2),0)
									else !Hirshfeld-I and MBIS. Refined atomic radial density of every atom has been available in atmraddens
										npt=atmradnpt(iatm)
										call lagintpol(atmradpos(1:npt),atmraddens(1:npt,iatm),npt,dsqrt(dist2),tmprho,rnouse,rnouse,1)
									end if
									atmrho(iatm)=atmrho(iatm)+tmprho
								end if
							end do
						end do
					end do
				end do
				prorho=sum(atmrho(:))
				if (prorho<1D-6) cycle !Accuracy loose is found to be negligible enough
				
                !If any atom whose AOM needs to calculate has weight on this grid larger than a threshold, then orbital values at this grid needs to be calculated
                icalcorb=0
                do idx=1,natmthis
					iatm=iatmbeg+idx-1
                    wei=atmrho(iatm)/prorho
                    if (wei>thres.and.atmdolist(iatm).eqv..true.) icalcorb=1
                end do
                if (icalcorb==0) cycle
                
				!Calculate total or alpha part. This is time-consuming for very large system
				call orbderv(1,1,nmatsize,tmpx,tmpy,tmpz,orbval)
                
                !Accumulate contribution of present integration grid to AOM of various atoms. This is most time-consuming for medium-sized system
                !The following code is found to be significant faster than the below one, at least for a tested small system
                do idx=1,natmthis
					iatm=iatmbeg+idx-1
                    if (atmdolist(iatm).eqv..false.) cycle !Skip atoms that do not need to evaluate AOM to reduce cost
                    wei=atmrho(iatm)/prorho
                    if (wei<thres) cycle
					do jmo=1,nmatsize
						tmpval=orbval(jmo)*wei
						do imo=jmo,nmatsize
							AOM_tmp(imo,jmo,idx)=AOM_tmp(imo,jmo,idx)+orbval(imo)*tmpval
						end do
					end do
                end do
                
                !The following code is slower
				!do jmo=1,nmatsize
				!	do imo=jmo,nmatsize
				!		AOM_tmp(imo,jmo,1:natmthis)=AOM_tmp(imo,jmo,1:natmthis)+orbval(imo)*orbval(jmo)*atmrho(iatmbeg:iatmend)/prorho
				!	end do
				!end do
                !The following code is extremely slower! Mostly because the matrix operation in atomic looping is too expensive
				!do jmo=1,nmatsize
				!	do imo=jmo,nmatsize
				!		tmpmat(imo,jmo)=orbval(imo)*orbval(jmo)
				!	end do
				!end do
				!do idx=1,natmthis
				!	iatm=iatmbeg+idx-1
                !	AOM_tmp(:,:,idx)=AOM_tmp(:,:,idx)+tmpmat(:,:)*atmrho(iatm)/prorho
				!end do
            
				!Calculate Beta part for UHF,U-post-HF
				if ((wfntype==1.or.wfntype==4).and.nmatsizeb>0) then
					MOinit=iendalpha+1
					MOend=iendalpha+nmatsizeb
                    call orbderv(1,MOinit,MOend,tmpx,tmpy,tmpz,orbval)
					do idx=1,natmthis
						iatm=iatmbeg+idx-1
						if (atmdolist(iatm).eqv..false.) cycle !Skip atoms that do not need to calculate to reduce cost
						wei=atmrho(iatm)/prorho
						if (wei<thres) cycle
						do jmo=MOinit,MOend
							tmpval=orbval(jmo)*wei
							do imo=jmo,MOend
								AOMb_tmp(imo-iendalpha,jmo-iendalpha,idx)=AOMb_tmp(imo-iendalpha,jmo-iendalpha,idx)+orbval(imo)*tmpval
							end do
						end do
                    end do
				end if
            
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
	!$OMP END DO
	!$OMP CRITICAL
	AOM(:,:,iatmbeg:iatmend)=AOM(:,:,iatmbeg:iatmend)+AOM_tmp(:,:,1:natmthis)
	if ((wfntype==1.or.wfntype==4).and.nmatsizeb>0) AOMb(:,:,iatmbeg:iatmend)=AOMb(:,:,iatmbeg:iatmend)+AOMb_tmp(:,:,1:natmthis)
	!$OMP END CRITICAL
	!$OMP END PARALLEL
	if (ishowprog/=0) call showprog(100,100)
end do

do jmo=1,nmatsize
	do imo=jmo+1,nmatsize
		AOM(jmo,imo,:)=AOM(imo,jmo,:)
	end do
end do
AOM=AOM*dvol
if ((wfntype==1.or.wfntype==4).and.nmatsizeb>0) then
	do jmo=1,nmatsizeb
		do imo=jmo+1,nmatsizeb
			AOMb(jmo,imo,:)=AOMb(imo,jmo,:)
		end do
	end do
	AOMb=AOMb*dvol
end if

call walltime(iwalltime2)
write(*,"(/,' Calculation totally took up wall clock time',i10,' s')") iwalltime2-iwalltime1

!Generate FOM and then export it
if (iwork==33) then
	!Generate FOMs
	allocate(FOM1(nmatsize,nmatsize),FOM2(nmatsize,nmatsize))
    FOM1=0
    FOM2=0
    do idx=1,nFOM1atm
		FOM1(:,:)=FOM1(:,:)+AOM(:,:,FOM1atm(idx))
    end do
    if (iFOMmode>1) then
		do idx=1,nFOM2atm
			FOM2(:,:)=FOM2(:,:)+AOM(:,:,FOM2atm(idx))
		end do
    end if
    if (wfntype==1.or.wfntype==2.or.wfntype==4) then
		allocate(FOM1b(nmatsizeb,nmatsizeb),FOM2b(nmatsizeb,nmatsizeb))
		FOM1b=0
		FOM2b=0
		do idx=1,nFOM1atm
			FOM1b(:,:)=FOM1b(:,:)+AOMb(:,:,FOM1atm(idx))
		end do
		if (iFOMmode>1) then
			do idx=1,nFOM2atm
				FOM2b(:,:)=FOM2b(:,:)+AOMb(:,:,FOM2atm(idx))
			end do
        end if
    end if
    !Export FOMs
	open(10,file="FOM.txt",status="replace")
	if (wfntype==0.or.wfntype==3) then !Closed-shell
		if (iFOMmode==1) then
			write(10,"('Fragment overlap matrix')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(*,*) "FOM has been exported to FOM.txt in current folder"
		else
			write(10,"('Fragment overlap matrix of fragment 1')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Fragment overlap matrix of fragment 2')")
			call showmatgau(FOM2(:,:),"",1,"f14.8",10)
			write(*,*) "FOM of fragments 1 and 2 has been exported to FOM.txt in current folder"
		end if
	else !Open-shell
		if (iFOMmode==1) then
			write(10,"('Alpha part of fragment overlap matrix')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(*,*)
			write(10,"('Beta part of fragment overlap matrix')")
			call showmatgau(FOM1b(:,:),"",1,"f14.8",10)
			write(*,*) "FOM has been exported to FOM.txt in current folder"
		else
			write(10,"('Alpha part of fragment overlap matrix of fragment 1')")
			call showmatgau(FOM1(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Beta part of fragment overlap matrix of fragment 1')")
			call showmatgau(FOM1b(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Alpha part of fragment overlap matrix of fragment 2')")
			call showmatgau(FOM2(:,:),"",1,"f14.8",10)
			write(10,*)
			write(10,"('Beta part of fragment overlap matrix of fragment 2')")
			call showmatgau(FOM2b(:,:),"",1,"f14.8",10)
			write(*,*) "FOM of fragments 1 and 2 has been exported to FOM.txt in current folder"
		end if
	end if
	close(10)
end if

end subroutine
