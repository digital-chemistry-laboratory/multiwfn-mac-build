!-------- Main interface of various other functions (part 1)
subroutine otherfunc_main
use defvar
implicit real*8 (a-h,o-z)
character c80tmp*80
do while(.true.)
	write(*,*)
	write(*,*) "              ============ Other functions (Part 1) ============ "
	write(*,*) "0 Return"
	write(*,*) "1 Draw scatter graph between two functions and generate their cube files"
	write(*,"(a)") " 2 Export various files (mwfn/pdb/xyz/wfn/wfx/molden/fch/47/mkl...) or generate input file of quantum chemistry programs"
	write(*,*) "3 Calculate molecular van der Waals volume"
	write(*,*) "4 Integrate a function in whole space"
	write(*,*) "5 Show overlap integral between alpha and beta orbitals"
	write(*,*) "6 Monitor SCF convergence process of Gaussian"
    write(*,*) "7 Auxiliary tools for CP2K (CP2Kmate)"
	write(*,*) "8 Generate Gaussian input file with initial guess from fragment wavefunctions"
	write(*,*) "9 Evaluate interatomic connectivity and atomic coordination number"
! 	write(*,*) "10 Generate spherically averaged atomic radial density" !Rarely used, so, hidden
	write(*,*) "11 Calculate overlap and centroid distance between two orbitals"
	write(*,*) "12 Biorthogonalization between alpha and beta orbitals"
	!write(*,*) "13 Calculate HOMA and Bird aromaticity index"
	write(*,*) "14 Calculate LOLIPOP (LOL Integrated Pi Over Plane)"
	write(*,*) "15 Calculate intermolecular orbital overlap"
    write(*,*) "17 Generate Fock/KS matrix based on orbital energies and coefficients"
	write(*,*) "18 Yoshizawa's electron transport route analysis"
	write(*,*) "19 Generate new wavefunction by combining fragment wavefunctions"
	write(*,*) "20 Calculate Hellmann-Feynman forces"
	write(*,*) "21 Calculate properties based on geometry information for specific atoms"
	write(*,*) "22 Detect pi orbitals, set occupation numbers and calculate pi composition"
	write(*,*) "23 Fit function distribution to atomic value"
	!write(*,*) "24 Obtain NICS_ZZ value for non-planar or tilted system"
	read(*,*) c80tmp

    if (c80tmp=="4a".or.c80tmp=="4b") then !For Chunying Rong's Fukui Shannon project
		if (c80tmp=="4a") call info_rhodiff
		if (c80tmp=="4b") call info_rhodiff_grid
        cycle
    else
		read(c80tmp,*) isel
    end if
	if (isel==0) then
		return
	else if (isel==1) then
		call funcvsfunc(0)
	else if (isel==2) then
		call outfile
	else if (isel==3) then
        call molvol_MC
	else if (isel==4) then
		if (ispecial/=1) then
			call selfunc_interface(1,ifunc)
			call intfunc(ifunc)
		else if (ispecial==1) then !For LSB
			call intfunc(1)
		end if
	else if (isel==-4) then
		call intdiff(1)
	else if (isel==-5) then
		call intdiff(2)
	else if (isel==5) then
		call aboverlap
	else if (isel==6) then
		call monitorscf
	else if (isel==7) then
		call cp2kmate
	else if (isel==8) then
		call fragguess
	else if (isel==9) then
		call conn_coordnum
	else if (isel==10) then
		call sphatmraddens
	else if (isel==11) then
		call ovlpdistorb
	else if (isel==12) then
		call biortho
	else if (isel==13) then !Hidden
		call HOMA_Bird
	else if (isel==14) then
		call LOLIPOP
	else if (isel==15) then
		call intmolovlp
	else if (isel==16) then
		write(*,*) " NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE"
		write(*,"(a)") " This function has been moved to main function 22. Please use main function 22 instead next time!!!"
		write(*,*) " NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE"
		call CDFT
	else if (isel==17) then
		call genFock
	else if (isel==18) then
		call Yoshieletrans
	else if (isel==19) then
		call gencombwfn
	else if (isel==20) then
		call hellmann_feynman
	else if (isel==21) then
		call calcgeomprop
	else if (isel==22) then
		call detectpiorb
	else if (isel==23) then
		call fitfunc
	else if (isel==24) then !Hidden
		call NICS_ZZ
	end if
end do
end subroutine




!! ----------------- Based on the equation for evaluating coordination number in DFT-D3 original paper
subroutine conn_coordnum
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 CNmat(ncenter,ncenter),k1,k2
character c80tmp*80,selectyn
write(*,*) "Input the threshold for printing connectivity index, e.g. 0.05"
write(*,*) "If you press ENTER button directly, 0.1 will be used"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
	printthres=0.1D0
else
	read(c80tmp,*) printthres
end if
k1=16D0
k2=4D0/3D0
CN=0D0
CNmat=0D0
do iatm=1,ncenter
	indi=a(iatm)%index
	do jatm=iatm+1,ncenter
		if (jatm==iatm) cycle
        if (ifPBC==0) then
    			r=atomdist(iatm,jatm,0)
        else
            call nearest_atmdistxyz(iatm,jatm,r,atmx,atmy,atmz)
        end if
		indj=a(jatm)%index
		sclcovsum=k2*(covr_pyy(indi)+covr_pyy(indj))
		CNmat(iatm,jatm)=1D0/( 1D0+dexp(-k1*(sclcovsum/r-1D0)) )
		if (CNmat(iatm,jatm)>=printthres) then
			write(*,"(i5,a,'  ---',i5,a,' :',f10.5,'   Nearest integer:',i3)") &
			iatm,ind2name(a(iatm)%index),jatm,ind2name(a(jatm)%index),CNmat(iatm,jatm),nint(CNmat(iatm,jatm))
		end if
	end do
end do
CNmat=CNmat+transpose(CNmat)
write(*,*)
do iatm=1,ncenter
	indi=a(iatm)%index
	write(*,"(i5,2x,a,2x,'Sum of connectivity:',f8.4,'   Sum of integer connectivity:',i3)") iatm,ind2name(indi),sum(CNmat(iatm,:)),sum(nint(CNmat(iatm,:)))
end do
write(*,*)
write(*,*) "If outputting connectivity matrix to connmat.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	open(10,file="connmat.txt",status="replace")
	call showmatgau(CNmat,"connectivity matrix",0,'f12.6',10)
	close(10)
	write(*,*) "Done! connmat.txt has been outputted to current folder"
end if
end subroutine



!! ----------- Output new file
subroutine outfile
use defvar
use util
implicit real*8 (a-h,o-z)
character c200tmp*200,c200tmp2*200,c200tmp3*200
write(*,*) "0 Return"
write(*,*) "                Export system to various formats of files:"
if (ifiletype==4) then !Used to convert .chg file to .pqr file
	write(*,*) "1 Output current structure and atomic charges to .pqr file"
	write(*,*) "-1 Output current structure and atomic charges to .pdb file"
else
	write(*,*) "1 Output current structure to .pdb file"
	!write(*,*) "1a Output current structure (with image atoms at boundary if any) to .pdb file" !Mainly used by AIM_PBC.bat
end if
write(*,*) "2 Output current structure to .xyz file"
write(*,*) "3 Output current structure and atomic charges to .chg file"
write(*,*) "4 Output current wavefunction as .wfx file"
write(*,*) "5 Output current wavefunction as .wfn file"
write(*,*) "6 Output current wavefunction as Molden input file (.molden)"
write(*,*) "7 Output current wavefunction as .fch file"
write(*,*) "8 Output current wavefunction as .47 file"
write(*,*) "9 Output current wavefunction as old Molekel input file (.mkl)"
write(*,*) "31 Output current structure to .cml file"
write(*,*) "32 Output current wavefunction as .mwfn file"
write(*,*) "33 Output current structure and cell information as .cif file"
write(*,*) "34 Output current structure and cell information as .gro file"
if (allocated(cubmat)) then
    write(*,*) "35 Output current grid data to Gaussian-type .cub file"
    write(*,*) "36 Output current grid data to .vti file"
    write(*,*) "37 Output current grid data to VASP grid data file"
end if
write(*,*) "             Generate input file of quantum chemistry codes:"
if (allocated(CObasa)) then
    write(*,*) "10 Gaussian"
else
    write(*,*) "10 Gaussian (with/without initial guess of wavefunction)"
end if
if (allocated(CObasa)) then
    write(*,*) "11 GAMESS-US"
else
    write(*,*) "11 GAMESS-US (with/without initial guess of wavefunction)"
end if
write(*,*) "12 ORCA         13 NWChem       14 MOPAC"
write(*,*) "15 PSI4         16 MRCC         17 CFOUR"
write(*,*) "18 Molpro       19 Dalton       20 Molcas"
write(*,*) "21 Q-Chem"
write(*,*) "25 CP2K         26 Quantum ESPRESSO      27 VASP (POSCAR)"
read(*,*) c200tmp
if (c200tmp=="1a") then
	call outpdb_PBC("mol.pdb",10)
    return
else
	read(c200tmp,*) isel
end if

if (isel==0) then
	return
else if (isel==1) then
	if (ifiletype==4) then
		call outpqr_wrapper
	else
		write(*,*) "Hint: You can also input ""pdb"" in main menu to quickly enter this function"
        write(*,*)
	    call outpdb_wrapper
	end if
else if (ifiletype==4.and.isel==-1) then
    write(*,*) "Hint: You can also input ""pdb"" in main menu to quickly enter this function"
    write(*,*)
	call outpdb_wrapper
else if (isel==2) then
    write(*,*) "Hint: You can also input ""xyz"" in main menu to quickly enter this function"
    write(*,*)
	call outxyz_wrapper
else if (isel==3) then
	write(*,*) "Input path for exporting file, e.g. C:\ltwd.chg"
	read(*,"(a)") c200tmp
	call outchg(c200tmp,10)
else if (isel==4) then
	if (.not.allocated(b)) then
		write(*,*) "Error: The input file you used does not contain GTF information!"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else
		write(*,*) "Input path for exporting file, e.g. C:\ltwd.wfx"
		read(*,"(a)") c200tmp
		call outwfx(c200tmp,1,10)
		write(*,*) "Done!"
	end if
else if (isel==5) then
	if (.not.allocated(b)) then
		write(*,*) "Error: The input file you used does not contain GTF information!"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else
		write(*,*) "Input path for exporting file, e.g. C:\ltwd.wfn"
		read(*,"(a)") c200tmp
		call outwfn(c200tmp,1,1,10)
		write(*,*) "Done!"
	end if
else if (isel==6) then
	if (.not.allocated(CObasa)) then
		write(*,"(a)") " Error: This function works only when input file contains basis function information"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else
        call outmolden_wrapper
	end if
else if (isel==7) then
	if (.not.allocated(CObasa)) then
		write(*,"(a)") " Error: This function works only when input file contains basis function information"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else
        call outfch_wrapper
	end if
else if (isel==8) then
	if (.not.allocated(CObasa)) then
		write(*,"(a)") " Error: This function works only when input file contains basis function information"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else
		write(*,*) "Input path for exporting file, e.g. C:\ltwd.47"
		read(*,"(a)") c200tmp
		call out47(c200tmp,10)
	end if
else if (isel==9) then
	if (.not.allocated(CObasa)) then
		write(*,"(a)") " Error: This function works only when input file contains basis function information"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else
		write(*,*) "Input path for exporting file, e.g. C:\ltwd.mkl"
		read(*,"(a)") c200tmp
		call outmkl(c200tmp,10)
	end if
else if (isel==10) then
	call outgjf_wrapper
else if (isel==11) then
	write(*,*) "Input path for generating file, e.g. C:\ltwd.inp"
	read(*,"(a)") c200tmp
	call outGAMESSinp(c200tmp,10)
else if (isel==12) then
    write(*,*) "Hint: You can also input ""oi"" in the main menu to enter this function"
    write(*,*)
	call outORCAinp_wrapper
else if (isel==13) then
	write(*,*) "Input path for generating file, e.g. C:\ltwd.nw"
	read(*,"(a)") c200tmp
	call outNWCheminp(c200tmp,10)
else if (isel==14) then
	call outMOPACinp_wrapper
else if (isel==15) then
    write(*,*) "Hint: You can also input ""pi"" in the main menu to enter this function"
    write(*,*)
	call outPSI4inp_wrapper
else if (isel==16) then
	c200tmp="MINP"
	call outMRCCinp(c200tmp,10)
else if (isel==17) then
	c200tmp="ZMAT"
	call outCFOURinp(c200tmp,10)
else if (isel==18) then
	write(*,*) "Input path for generating file, e.g. C:\ltwd.inp"
	read(*,"(a)") c200tmp
	call outmolproinp(c200tmp,10)
else if (isel==19) then
	c200tmp=" "
	write(*,"(a)") " Input path of .dal file, e.g. C:\DFT.dal (Directly pressing ENTER button if you do not need it)"
	read(*,"(a)") c200tmp
    call path2filename(filename,c200tmp3)
    c200tmp3=trim(c200tmp3)//".mol"
	write(*,*) "Input path of .mol file, e.g. C:\ltwd.mol"
    write(*,"(a)") " If pressing ENTER button directly, "//trim(c200tmp3)//" will be generated in current folder"
	read(*,"(a)") c200tmp2
    if (c200tmp2==" ") c200tmp2=c200tmp3
	call outDaltoninp(c200tmp,c200tmp2,10)
else if (isel==20) then
	write(*,*) "Input path for generating file, e.g. C:\ltwd.inp"
	read(*,"(a)") c200tmp
	call outmolcasinp(c200tmp,10)
else if (isel==21) then
	write(*,*) "Input path for generating file, e.g. C:\ltwd.inp"
	read(*,"(a)") c200tmp
	call outQcheminp(c200tmp,10)
else if (isel==25) then
    write(*,*) "Hint: You can also input ""cp2k"" in the main menu to enter this function"
    write(*,*)
	call outCP2Kinp_wrapper
else if (isel==26) then
    write(*,*) "Hint: You can also input ""QE"" in the main menu to enter this function"
    write(*,*)
	call outQEinp_wrapper
else if (isel==27) then
	call outPOSCAR_wrapper
else if (isel==31) then
	write(*,*) "Input path for exporting file, e.g. C:\ltwd.cml"
	read(*,"(a)") c200tmp
	call outcml(c200tmp,10,0)
else if (isel==32) then
    call outmwfn_wrapper
else if (isel==33) then
    call outcif_wrapper
else if (isel==34) then
    call outgro_wrapper
else if (isel==35) then
    call outcube_wrapper
else if (isel==36) then
	write(*,*) "Input path for exporting file, e.g. C:\ltwd.vti"
	read(*,"(a)") c200tmp
    write(*,*) "Exporting..."
	call outvti(c200tmp,10)
else if (isel==37) then
	write(*,*) "Input path for exporting file, e.g. D:\CHGCAR"
    write(*,"(a)") " If press ENTER button directly, the file will be exported to CHGCAR in current folder"
	read(*,"(a)") c200tmp
    if (c200tmp==" ") c200tmp="CHGCAR"
    write(*,*) "Exporting..."
    call outVASPgrd(c200tmp,10)
end if
end subroutine



!!----------- function vs. function
!iwork=0: General case
!iwork=1: NCI
!iwork=2: NCI based on promolecular approximation
!iwork=4: IRI
!iwork=5: DORI
subroutine funcvsfunc(iwork)
use plot
use defvar
use util
use GUI
implicit real*8 (a-h,o-z)
integer iwork
real*8,allocatable :: scatterx(:),scattery(:),exchangedata(:,:,:)
real*8,allocatable :: f2orgdata(:,:,:),scatteryorg(:) !Backup original data of function 2, because it may be modified by users
character c200tmp*200,f1name*20,f2name*20

f1name="function 1"
f2name="function 2"
if (iwork==1) then !NCI
	iselfunc1=15
	iselfunc2=13
    write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
    write(*,*) "  Original paper of NCI: J. Am. Chem. Soc., 132, 6498 (2010)"
	write(*,"(a)") "   Comprehensive reviews:"
    write(*,"(a)") " Tian Lu, Qinxue Chen, Visualization Analysis of &
	&Weak Interactions in Chemical Systems. In Comprehensive Computational Chemistry, vol. 2, pp. 240-264. Oxford: Elsevier (2024) DOI: 10.1016/B978-0-12-821978-2.00076-3"
    write(*,"(a)") " Tian Lu, Visualization Analysis of Covalent and Noncovalent Interactions in Real Space, Angew. Chem. Int. Ed., e202504895 (2025) DOI: 10.1002/anie.202504895"
else if (iwork==2) then !NCI based on promolecular approximation
	iselfunc1=16
	iselfunc2=14
    write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
    write(*,*) "  Original paper of NCI: J. Am. Chem. Soc., 132, 6498 (2010)"
	write(*,"(a)") "   Comprehensive reviews:"
    write(*,"(a)") " Tian Lu, Qinxue Chen, Visualization Analysis of &
	&Weak Interactions in Chemical Systems. In Comprehensive Computational Chemistry, vol. 2, pp. 240-264. Oxford: Elsevier (2024) DOI: 10.1016/B978-0-12-821978-2.00076-3"
    write(*,"(a)") " Tian Lu, Visualization Analysis of Covalent and Noncovalent Interactions in Real Space, Angew. Chem. Int. Ed., e202504895 (2025) DOI: 10.1002/anie.202504895"
else if (iwork==4) then !IRI
	iselfunc1=15
	iuserfunc_old=iuserfunc
	iselfunc2=100
	iuserfunc=99
    write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
    write(*,"(a)") "   Original paper of IRI: Tian Lu, Qinxue Chen, Interaction Region Indicator (IRI): A Simple Real Space Function &
    &Clearly Revealing Both Chemical Bonds and Weak Interactions, Chemistry-Methods, 1, 231-239 (2021) DOI: 10.1002/cmtd.202100007"   
	write(*,"(a)") "   Comprehensive reviews:"
    write(*,"(a)") " Tian Lu, Qinxue Chen, Visualization Analysis of &
	&Weak Interactions in Chemical Systems. In Comprehensive Computational Chemistry, vol. 2, pp. 240-264. Oxford: Elsevier (2024) DOI: 10.1016/B978-0-12-821978-2.00076-3"
    write(*,"(a)") " Tian Lu, Visualization Analysis of Covalent and Noncovalent Interactions in Real Space, Angew. Chem. Int. Ed., e202504895 (2025) DOI: 10.1002/anie.202504895"
else if (iwork==5) then !DORI
	iselfunc1=15
	iuserfunc_old=iuserfunc
	iselfunc2=100
	iuserfunc=20
    write(*,"(a,/)") " NOTE: Interaction Region Indicator (IRI) works much better than DORI, it is strongly suggested to use IRI instead!"
else
	call funclist
	if (allocated(b)) write(*,*) "Select function 1 (as X axis) and function 2 (as Y axis)  e.g. 15,13"
	if (.not.allocated(b)) write(*,*) "Select function 1 (as X axis) and function 2 (as Y axis)  e.g. 16,14"
	if (ifiletype==7) write(*,"(a)") " Note: If input 0,0, then the grid data in memory will be directly taken as function 1, &
	&and you need to input path of a cube file, whose data will be taken as function 2" 
	read(*,*) iselfunc1,iselfunc2
	if (iselfunc1==4) then
		write(*,*) "Select which orbital (for function 1)? Input the index"
		read(*,*) iorbsel1
	end if
	if (iselfunc2==4) then
		write(*,*) "Select which orbital (for function 2)? Input the index"
		read(*,*) iorbsel2
	end if
end if
if (iselfunc1==15.or.iselfunc1==16) f1name="sign(lambda2)rho"
if (iselfunc2==13.or.iselfunc2==14) f2name="RDG"
if (iselfunc2==24.or.(iselfunc2==100.and.iuserfunc==99)) f2name="IRI"
if (iselfunc2==100.and.iuserfunc==20) f2name="DORI"

if (iselfunc1==0.and.iselfunc2==0) then !Directly load grid data
	if (allocated(cubmattmp)) deallocate(cubmattmp)
	write(*,*) "Input file path for cube file of function 2, e.g. C:\test.cub"
	read(*,"(a)") c200tmp
	call readcubetmp(c200tmp,1,inconsis)
	if (inconsis==1) then
		write(*,"(a)") " Error: The grid setting of this cube file is inconsistent with that of present grid data, exit..."
		read(*,*)
		return
	end if
else !Calculate grid data
	if (iwork==1.or.iwork==2.or.iwork==4.or.iwork==5) aug3D=1.5D0 !Smaller than default value
	call setgrid(0,igridsel)
	if (allocated(cubmat)) deallocate(cubmat)
	if (allocated(cubmattmp)) deallocate(cubmattmp)
	allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
	call delvirorb(1)
	if (iselfunc1==15.and.iselfunc2==13) then !Since RDG and sign(lambda2)rho is often combined to study NCI, a special code is provided for reducing cost
		call savecubmat(1513,0,1)
	else if (iselfunc1==16.and.iselfunc2==14) then 
		call savecubmat(1614,0,1)
	else if (iselfunc1==15.and.(iselfunc2==24.or.(iselfunc2==100.and.iuserfunc==99))) then !Combinedly calculate IRI and sign(lambda2)rho to reduce cost by double
		call savecubmat(1599,0,1)
	else
		write(*,"(a)") " Calculating "//trim(f2name)//"..."
		call savecubmat(iselfunc2,0,iorbsel2)
		cubmattmp=cubmat
		write(*,"(a)") " Calculating "//trim(f1name)//"..."
		call savecubmat(iselfunc1,0,iorbsel1)
	end if
	call delvirorb_back(1)
end if

!Store grid data to array
allocate(scatterx(nx*ny*nz),scattery(nx*ny*nz))
ii=1
do k=1,nz
	do j=1,ny
		do i=1,nx
			scatterx(ii)=cubmat(i,j,k)
			scattery(ii)=cubmattmp(i,j,k)
			ii=ii+1
		end do
	end do
end do

if (iwork==0) then
	func1stddev=stddevarray(scatterx)
	func2stddev=stddevarray(scattery)
	write(*,"(a,2E18.10)") " Standard deviation of function 1 and 2:",func1stddev,func2stddev
	pearsoncoeff=covarray(scatterx,scattery)/func1stddev/func2stddev
	write(*,"(2(a,f12.6))") " Pearson correlation coefficient r:",pearsoncoeff,"  r^2:",pearsoncoeff**2
end if

xmin=minval(scatterx)
xmax=maxval(scatterx)
ymin=minval(scattery)
ymax=maxval(scattery)
if (iselfunc1==15.and.iselfunc2==13) then !sign(lambda2)*rho vs. RDG
	xmin=-RDG_maxrho
	xmax=RDG_maxrho
	if (RDG_maxrho==0D0) xmin=-2D0
	if (RDG_maxrho==0D0) xmax=2D0
	ymin=0D0
	ymax=2D0
else if (iselfunc1==16.and.iselfunc2==14) then !sign(lambda2)*rho vs. RDG based on promolecular density
	xmin=-RDGprodens_maxrho
	xmax=RDGprodens_maxrho
	if (RDGprodens_maxrho==0D0) xmin=-2D0
	if (RDGprodens_maxrho==0D0) xmax=2D0
	ymin=0D0
	ymax=2D0
else if (iselfunc1==15.and.(iselfunc2==24.or.(iselfunc2==100.and.iuserfunc==99))) then !sign(lambda2)*rho vs. IRI
	xmin=-0.4D0
	xmax=0.1D0
	ymin=0D0
	ymax=2.5D0
end if

do while (.true.)
	write(*,*)
	if (allocated(f2orgdata)) write(*,"(a)") " -4 Restore original "//trim(f2name)//" grid data"
	write(*,"(a)") " -3 Set "//trim(f2name)//" value where value of "//trim(f1name)//" is out of a certain range"
	write(*,"(a)") " -2 Set "//trim(f2name)//" value where value of "//trim(f1name)//" is within a certain range"
	write(*,*) "-1 Draw scatter graph"
	write(*,*) "0 Exit"
	write(*,*) "1 Save the scatter graph to file"
	write(*,*) "2 Output scatter points to output.txt in current folder"
	write(*,*) "3 Output cube files to func1.cub and func2.cub in current folder"
	write(*,"(' 4 Change range of X-axis of scatter graph, current:',1PE12.4,' to',1PE12.4)") xmin,xmax
	write(*,"(' 5 Change range of Y-axis of scatter graph, current:',1PE12.4,' to',1PE12.4)") ymin,ymax
	write(*,"(a)") " 6 Show isosurface of "//trim(f1name)
	write(*,"(a)") " 7 Show isosurface of "//trim(f2name)
	write(*,"(a)") " 8 Output "//trim(f1name)//" to output.txt where "//trim(f2name)//" is within in certain range"
    if (iwork==4) write(*,"(a)") " 9 Screen out covalent bond regions (set IRI to 100 for regions with sign(lambda2)rho < -0.1 a.u.)"
	read(*,*) isel
	if (isel==-4) then
		cubmattmp=f2orgdata
        scattery=scatteryorg
        deallocate(f2orgdata,scatteryorg)
        write(*,*) "Done! Original grid data and scatter data have been restored"
	else if (isel==-2.or.isel==-3) then
		if (.not.allocated(f2orgdata)) then
			allocate(f2orgdata(nx,ny,nz)) !Back up grid data of function 2
			f2orgdata=cubmattmp
            allocate(scatteryorg(nx*ny*nz)) !Back up scatter data of function 2
            scatteryorg=scattery
        end if
		write(*,"(a)") " Input lower and upper limit of the range of "//trim(f1name)//", e.g. 0.5,2.3"
		read(*,*) rlower,rupper
		write(*,"(a)") " Input the expected value of "//trim(f2name)//", e.g. 100"
		read(*,*) rsetfunc2
		if (isel==-2) then
			where (scatterx<=rupper.and.scatterx>=rlower) scattery=rsetfunc2
			where (cubmat<=rupper.and.cubmat>=rlower) cubmattmp=rsetfunc2
		else if (isel==-3) then
			where (scatterx>=rupper.or.scatterx<=rlower) scattery=rsetfunc2
			where (cubmat>=rupper.or.cubmat<=rlower) cubmattmp=rsetfunc2
		end if
		write(*,*) "Done! Then if you want to restore original data, you can use option -4"
	else if (isel==-1) then
		write(*,*) "Drawing graph, please wait..."
		if ((iselfunc1==15.and.iselfunc2==13).or.(iselfunc1==16.and.iselfunc2==14)) then
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","Reduced density gradient")
		else if (iselfunc1==15.and.iselfunc2==100.and.iuserfunc==20) then
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","DORI")
		else if (iselfunc1==15.and.(iselfunc2==24.or.(iselfunc2==100.and.iuserfunc==99))) then
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","IRI (a.u.)")
		else
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1)
		end if
	else if (isel==0) then
		if (iwork==4.or.iwork==5) iuserfunc=iuserfunc_old !IRI, DORI
		exit
	else if (isel==1) then
		isavepic=1
		write(*,*) "Exporting image file, please wait..."
		if ((iselfunc1==15.and.iselfunc2==13).or.(iselfunc1==16.and.iselfunc2==14)) then
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","Reduced density gradient")
		else if (iselfunc1==15.and.iselfunc2==100.and.iuserfunc==20) then
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","DORI")
		else if (iselfunc1==15.and.(iselfunc2==24.or.(iselfunc2==100.and.iuserfunc==99))) then
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","IRI (a.u.)")
		else
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1)
		end if
		isavepic=0
		write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""dislin"" prefix in current directory"
	else if (isel==2) then
		write(*,*) "Outputting output.txt in current folder..."
		open(10,file="output.txt",status="replace")
		do k=1,nz
			do j=1,ny
				do i=1,nx
					call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					write(10,"(3f12.6,2E16.8)") tmpx*b2a,tmpy*b2a,tmpz*b2a,cubmat(i,j,k),cubmattmp(i,j,k)
				end do
			end do
            call showprog(k,nz)
		end do
		close(10)
		write(*,"(a)") " Finished!"
        write(*,*) "Column 1/2/3: X/Y/Z in Angstrom"
        write(*,*) "Column 4/5: "//trim(f1name)//" and "//trim(f2name)//" in a.u."
		write(*,"(a)") " Obviously, if you will plot scatter map between "//trim(f1name)//" and "//trim(f2name)//" in external tools such as Origin, &
		&the last two columns should be taken as X and Y axes data"
	else if (isel==3) then
		write(*,*) "Exporting..."
		open(10,file="func1.cub",status="replace")
		call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,"(a)") " The cube file of "//trim(f1name)//" has been exported to func1.cub in current folder"
		open(10,file="func2.cub",status="replace")
		call outcube(cubmattmp,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,"(a)") " The cube file of "//trim(f2name)//" has been exported to func2.cub in current folder"
	else if (isel==4) then
		write(*,*) "Input lower limit and upper limit of X axis e.g. 0,1.5"
		read(*,*) xmin,xmax
	else if (isel==5) then
		write(*,*) "Input lower limit and upper limit of Y axis e.g. 0,1.5"
		read(*,*) ymin,ymax
	else if (isel==6) then
	 	write(*,*) "Input the value of isosurface, e.g. 0.02"
		read(*,*) sur_value
		call drawisosurgui(1)
	else if (isel==7) then
        allocate(exchangedata(nx,ny,nz))
		exchangedata=cubmat
		cubmat=cubmattmp
	 	write(*,*) "Input the value of isosurface, e.g. 0.5"
		read(*,*) sur_value
		call drawisosurgui(1)
		cubmat=exchangedata
        deallocate(exchangedata)
	else if (isel==8) then
		write(*,*) "Input range of "//trim(f2name)//", e.g. 0.0009, 0.0011"
		read(*,*) rlowlim,uplim
		rmin=1D200
		rmax=-1D200
		num=0
		open(10,file="output.txt",status="replace")
		do k=1,nz
			do j=1,ny
				do i=1,nx
					if (cubmattmp(i,j,k)<=uplim.and.cubmattmp(i,j,k)>=rlowlim) then
						num=num+1
						if (cubmattmp(i,j,k)>rmax) rmax=cubmattmp(i,j,k)
						if (cubmattmp(i,j,k)<rmin) rmin=cubmattmp(i,j,k)
                        call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
						write(10,"(3f12.6,2E16.8)") tmpx*b2a,tmpy*b2a,tmpz*b2a,cubmat(i,j,k),cubmattmp(i,j,k)
					end if
				end do
			end do
			call showprog(k,nz)
		end do
		close(10)
		write(*,"(a)") " Finished!"
        write(*,*) "Column 1/2/3: X/Y/Z in Angstrom"
        write(*,*) "Column 4/5: "//trim(f1name)//" and "//trim(f2name)//" in a.u."
		write(*,"(' Number of entries:',i10)") num
		if (num>=2) write(*,"(a,2E16.8)") " Min and Max of "//trim(f1name)//" in this range",rmin,rmax
	else if (isel==9) then
        where (scatterx<=-0.1D0) scattery=100
		where (cubmat<=-0.1D0) cubmattmp=100
        write(*,*) "Done!"
	end if
end do
end subroutine





!!---------- Intermolecular MO overlap
!Both alpha and beta are taken into account, the beta index is after alpha index
!Two cases of input files
!1 First load dimer Gaussian output file (must with iop(3/33=1), then load the two monomer Gaussian output files in turn (must with pop=full)
!2 Load dimer wavefunction, compute Sbas, then respectively load two monomer wavefunctions
subroutine intmolovlp
use util
use defvar
implicit real*8 (a-h,o-z)
real*8,allocatable :: cobas1(:,:),cobas2(:,:),ovlpbasmat(:,:),orbovlp(:,:) !ovlpmat is overlap matrix of basis functions
integer,allocatable :: elemidxcomb(:) !Element index array, used for sanity check of case 2
character monofile1*200,monofile2*200,c80tmp*80

if (ifiletype==0) then !Plain text file, assumed to be gaussian output file
    open(10,file=filename,status="old")
    call loclabel(10,"NBasis=",ifound) !Number of basis functions
    read(10,*) c80tmp,nbasis
    write(*,"(' The number of basis functions in the dimer',i10)") nbasis
    allocate(ovlpbasmat(nbasis,nbasis))
    write(*,*) "Loading overlap matrix of dimer, please wait..."
    call loclabel(10,"*** Overlap ***",ifound)
    call readmatgau(10,ovlpbasmat,1,"D14.6",7,5)
    close(10)

    write(*,*)
    write(*,*) "Input Gaussian output file of monomer 1, e.g. C:\monomer1.out"
    do while(.true.)
	    read(*,"(a)") monofile1
	    inquire(file=monofile1,exist=alive)
	    if (alive) exit
	    write(*,*) "Error: File cannot be found, input again"
    end do
    open(10,file=monofile1,status="old")
    call loclabel(10,"NBasis=",ifound)
    read(10,*) c80tmp,nbasis1 !Number of basis functions in monomer 1
    call loclabel(10,"NBsUse=",ifound) !NbsUse must equal to the number of MOs
    read(10,*) c80tmp,nmo1
    write(*,"(' The number of basis functions in monomer 1',i10)") nbasis1
    call loclabel(10,"Beta Molecular Orbital Coefficients:",iopsh1,0) !Determine if this monomer is open-shell
	allocate(cobas1(nbasis,nmo1))
	cobas1=0D0
    if (iopsh1==1) then
	    nmo1=nmo1*2
	    write(*,"(' MOs from',i8,' to',i8,' are Alpha orbitals')") 1,nmo1/2
	    write(*,"(' MOs from',i8,' to',i8,' are Beta orbitals')") nmo1/2+1,nmo1
	    write(*,*) "Loading molecular orbital coefficients of monomer 1, please wait..."
	    call loclabel(10,"Alpha Molecular Orbital Coefficients:",ifound)
	    call readmatgau(10,cobas1(1:nbasis1,1:nmo1/2),0,"f10.5",21,5,3) !nbasis1+1:nbasis are empty
	    call loclabel(10,"Beta Molecular Orbital Coefficients:",ifound)
	    call readmatgau(10,cobas1(1:nbasis1,nmo1/2+1:),0,"f10.5",21,5,3) !nbasis1+1:nbasis are empty
    else !Closed-shell
	    write(*,"(' The number of molecular orbitals in monomer 1',i10)") nmo1
	    write(*,*) "Loading molecular orbital coefficients of monomer 1, please wait..."
	    call loclabel(10,"Molecular Orbital Coefficients:",ifound)
	    call readmatgau(10,cobas1(1:nbasis1,:),0,"f10.5",21,5,3) !nbasis1+1:nbasis are empty
    end if
    close(10)

    write(*,*)
    write(*,*) "Input Gaussian output file of monomer 2, e.g. C:\monomer2.out"
    do while(.true.)
	    read(*,"(a)") monofile2
	    inquire(file=monofile2,exist=alive)
	    if (alive) exit
	    write(*,*) "Error: File cannot be found, input again"
    end do
    open(10,file=monofile2,status="old")
    call loclabel(10,"NBasis=",ifound)
    read(10,*) c80tmp,nbasis2 !Number of basis functions in monomer 1
    if (nbasis1+nbasis2/=nbasis) write(*,*) "Warning: The sum of the number of basis functions of the two monomers is unequal to dimer!"
    call loclabel(10,"NBsUse=",ifound) !NbsUse must equal to the number of MOs
    read(10,*) c80tmp,nmo2
    write(*,"(' The number of basis functions in monomer 2',i10)") nbasis2
    call loclabel(10,"Beta Molecular Orbital Coefficients:",iopsh2,0) !Determine if this monomer is open-shell
	allocate(cobas2(nbasis,nmo2))
	cobas2=0D0
    if (iopsh2==1) then
	    nmo2=nmo2*2
	    write(*,"(' MOs from',i8,' to',i8,' are Alpha orbitals')") 1,nmo2/2
	    write(*,"(' MOs from',i8,' to',i8,' are Beta orbitals')") nmo2/2+1,nmo2
	    write(*,*) "Loading molecular orbital coefficients of monomer 2, please wait..."
	    call loclabel(10,"Alpha Molecular Orbital Coefficients:",ifound)
	    call readmatgau(10,cobas2(nbasis1+1:,1:nmo2/2),0,"f10.5",21,5,3) !1:nbasis1 are empty
	    call loclabel(10,"Beta Molecular Orbital Coefficients:",ifound)
	    call readmatgau(10,cobas2(nbasis1+1:,nmo2/2+1:),0,"f10.5",21,5,3) !1:nbasis1 are empty
    else !Closed-shell
	    write(*,"(' The number of molecular orbitals in monomer 2',i10)") nmo2
	    write(*,*) "Loading molecular orbital coefficients of monomer 2, please wait..."
	    call loclabel(10,"Molecular Orbital Coefficients:",ifound)
	    call readmatgau(10,cobas2(nbasis1+1:,:),0,"f10.5",21,5,3) !1:nbasis1 are empty
    end if
    close(10)
    ! call showmatgau(cobas1,"111",0,"f14.8",6)

else !Using wavefunction file of dimer and monomer to do the analysis
    allocate(ovlpbasmat(nbasis,nbasis))
    ovlpbasmat=Sbas
    nmoall=nmo
    nbasisall=nbasis
    call dealloall(0)
    
    write(*,*) "Input wavefunction of monomer 1, e.g. C:\monomer1.fch"
    do while(.true.)
	    read(*,"(a)") monofile1
	    inquire(file=monofile1,exist=alive)
	    if (alive) exit
	    write(*,*) "File not found, input again"
    end do
    call readinfile(monofile1,1)
    write(*,"(' The number of basis functions in monomer 1',i10)") nbasis
    nmo1=nmo
    nbasis1=nbasis
	allocate(cobas1(nbasisall,nmo))
    cobas1=0
    iopsh1=0
    if (wfntype==1) then !Unrestricted
	    write(*,"(' Note: MOs from',i8,' to',i8,' are Alpha orbitals')") 1,nmo/2
	    write(*,"(' Note: MOs from',i8,' to',i8,' are Beta orbitals')") nmo/2+1,nmo
	    cobas1(1:nbasis,1:nmo/2)=CObasa
	    cobas1(1:nbasis,nmo/2+1:)=CObasb
        iopsh1=1
    else !R or RO
        cobas1(1:nbasis,:)=CObasa
        iopsh1=0
    end if
    
    allocate(elemidxcomb(ncenter_org))
    elemidxcomb(1:ncenter)=a%index
    ncenterfrag1=ncenter
    call dealloall(0)
    
    write(*,*)
    write(*,*) "Input wavefunction of monomer 2, e.g. C:\monomer2.fch"
    do while(.true.)
	    read(*,"(a)") monofile2
	    inquire(file=monofile2,exist=alive)
	    if (alive) exit
	    write(*,*) "File not found, input again"
    end do
    call readinfile(monofile2,1)
    write(*,"(' The number of basis functions in monomer 2',i10)") nbasis
    if (nbasisall/=nbasis1+nbasis) then
        write(*,"(a)") " Error: Sum of number of basis functions in the two monomers is inequivalent to that of the complex!"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        return
    end if
    
    nmo2=nmo
	allocate(cobas2(nbasisall,nmo))
    cobas2=0
    if (wfntype==1) then !Unrestricted
	    write(*,"(' Note: MOs from',i8,' to',i8,' are Alpha orbitals')") 1,nmo/2
	    write(*,"(' Note: MOs from',i8,' to',i8,' are Beta orbitals')") nmo/2+1,nmo
	    cobas2(nbasis1+1:,1:nmo/2)=CObasa
	    cobas2(nbasis1+1:,nmo/2+1:)=CObasb
        iopsh2=1
    else !R or RO
        cobas2(nbasis1+1:,:)=CObasa
        iopsh2=0
    end if
    elemidxcomb(ncenterfrag1+1:)=a%index
    call dealloall(0)
    
    !Check correspondence of fragments and dimer
    do iatm=1,ncenter_org
	    if (elemidxcomb(iatm)/=a_org(iatm)%index) then
		    write(*,"(/,a)") " Error: The sequence of the atoms in the fragments is not consistent with that in complex, the result will be meaningless! Possible reasons:"
		    write(*,"(a)") " 1 The fragment coordinates were not directly extracted from complex coordinate"
		    write(*,"(a)") " 2 The loading sequence of the fragments is not identical to occurrence sequence of the fragments in complex"
		    write(*,*) "Press ENTER button to exit"
		    read(*,*)
		    return
	    end if
    end do
    
    write(*,*)
    write(*,*) "Reloading the file initially loaded..."
    call readinfile(firstfilename,1)

end if

!Transform the overlap matrix in basis functions to the one between the MOs of the two monomers
write(*,*)
write(*,*) "Calculating the overlap matrix between MOs of the two monomers..."
allocate(orbovlp(nmo1,nmo2))
orbovlp=matmul(transpose(cobas1),matmul(ovlpbasmat,cobas2))

do while(.true.)
	write(*,*)
	write(*,"(a)") " Input e.g. 78,79 can print overlap integral between MO78 of monomer 1 and MO79 of monomer 2"
	write(*,"(a)") " Input o can output the whole overlap integral matrix to ovlpint.txt in current folder. Input q can exit"
	read(*,"(a)") c80tmp
	if (c80tmp(1:1)=='q') then
		exit
	else if (c80tmp(1:1)=='o') then
		open(10,file="ovlpint.txt",status="replace")
		call showmatgau(orbovlp,"",0,"f14.8",10)
		close(10)
		write(*,"(a)") " Done! The matrix has been outputted to ovlpint.txt in current folder"
		write(*,"(a)") " Note: The i,j element in the matrix corresponds to the overlap integral between MO i in monomer 1 and MO j in monomer 2"
	else
		read(c80tmp,*) idx1,idx2
		if (idx1>nmo1.or.idx1<1) then
			write(*,"(' Input error! The MO range of monomer 1 is between',i8,' and',i8)") 1,nmo1
			cycle
		end if
		if (idx2>nmo2.or.idx2<1) then
			write(*,"(' Input error! The MO range of monomer 2 is between',i8,' and',i8)") 1,nmo2
			cycle
		end if
		write(*,"(' Overlap integral is',f14.8)") orbovlp(idx1,idx2)
	end if
end do
end subroutine




!!----------- Generate fragment MO initial guess
subroutine fragguess
use util
use defvar
implicit real*8 (a-h,o-z)
integer,allocatable :: iflipfrag(:),numatom(:),numaelec(:),numbelec(:),numbas(:),ifunrestrict(:)
integer,allocatable :: wherefraga(:),wherefragb(:) !Each complex MO is come from which fragment
real*8,allocatable :: tmpcomat(:,:)
character selectyn*1,ctitle*80,c80tmp*80
character,allocatable :: namearray(:)*80

write(*,*) "How many fragments? (Including the fragment 1 that has been loaded)"
read(*,*) nfrag
allocate(iflipfrag(nfrag)) 
allocate(namearray(nfrag))
allocate(numatom(nfrag))
allocate(numbas(nfrag))
allocate(numaelec(nfrag))
allocate(numbelec(nfrag))
allocate(ifunrestrict(nfrag))
iflipfrag=0 !Don't flip spin by default
ifunrestrict=0 !All fragment are closed-shell by default
nchargetot=0 !charge of complex
nbasis=0
ncenter=0
naelec=0
nbelec=0

do i=1,nfrag
	if (i==1) then
		write(*,"(' Filename of fragment 1: ',a)") trim(filename)
	else if (i/=1) then
		do while(.true.)
			write(*,"(/,' Input Gaussian outputted filename of fragment',i4)") i
			read(*,"(a)") filename
			inquire(file=filename,exist=alive)
			if (alive) exit
			write(*,*) "Error: File not found, input again!"
		end do
	end if
	namearray(i)=filename

	open(10,file=filename,status="old") !Load other information of this fragment
	if (i==1) then
		call loclabel(10,"#")
		read(10,"(a)") ctitle
	end if
	call loclabel(10,"Charge =")
	read(10,"(9x,i3,15x,i2)") icharge,imulti
	call loclabel(10,"NBasis=")
	read(10,*) c80tmp,numbas(i)
	call loclabel(10,"NBsUse=")
	read(10,*) c80tmp,nbsuse
	if (nbsuse/=numbas(i)) then
		write(*,"(a)") " Error: Some linearly dependent basis functions were removed by Gaussian! You should regenerate the Gaussian output file with IOp(3/32=2)!"
		return
	end if
	call loclabel(10,"alpha electrons",ifound)
	read(10,*) numaelec(i),c80tmp,c80tmp,numbelec(i)
	call loclabel(10,"NAtoms",ifound)
	read(10,*) c80tmp,numatom(i)
	call loclabel(10,"Beta Molecular Orbital",ifunrestrict(i))
	close(10)

	write(*,"(' Charge and multiplicity of this fragment:',2i4)") icharge,numaelec(i)-numbelec(i)+1
	nchargetot=nchargetot+icharge
	if (imulti>1) then
		write(*,*) "Flip electron spin of this fragment? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			iflipfrag(i)=1
			itmp=numaelec(i)
			numaelec(i)=numbelec(i)
			numbelec(i)=itmp
		end if
	end if
end do

nbasis=sum(numbas)
naelec=sum(numaelec)
nbelec=sum(numbelec)
itotmulti=abs(naelec-nbelec)+1
write(*,"(' Total charge and multiplicity:',2i4)") nchargetot,itotmulti
write(*,"(' Total number of alpha and beta electrons:',2i6)") nint(naelec),nint(nbelec)
write(*,"(' Total number of basis functions:',i8)") nbasis
allocate(CObasa(nbasis,nbasis))
CObasa=0D0
allocate(wherefraga(nbasis)) !If wherefraga(i)=j means alpha orbital i in complex is contributed from fragment j
if (all(ifunrestrict==0)) then !Total is restricted
	allocate(MOocc(nbasis))
else !Since at least one fragment is unrestricted, so total is unrestricted
	allocate(MOocc(2*nbasis))
	allocate(CObasb(nbasis,nbasis))
	CObasb=0D0
	allocate(wherefragb(nbasis))
end if
ncenter=sum(numatom)
if (allocated(a)) deallocate(a)
allocate(a(ncenter))

iatm=1
locbas=1 !Current basis function position (index)
ioccmoa=1 !current alpha occupied orbital position (index)
ivirmoa=naelec+1 !current alpha virtual orbital position (index)
ioccmob=1
ivirmob=nbelec+1
do i=1,nfrag !Read atomic information and converged wavefunctions
	write(*,"(' Now loading fragment',i6,' ...')") i
	open(10,file=namearray(i),status="old")

	allocate(tmpcomat(numbas(i),numbas(i)))
	call loclabel(10,"Input orientation:",ifound) !I assume "nosymm" is used, so "standard orentation" does not appear
	if (ifound==0) call loclabel(10,"Coordinates in L301:",ifound)
 	if (ifound==0) call loclabel(10,"Z-Matrix orientation:",ifound) !sometimes the output of Gaussian is very weird
	if (ifound==0) then
	    write(*,"(a,/)") " Error: Atomic coordinates cannot be found!"
	    return
	end if
	do ii=1,5
		read(10,*)
	end do
	do ii=iatm,iatm+numatom(i)-1
		!Note: As usual, coordinate in Multiwfn is Bohr, but here is Angstrom for simplicity
		read(10,*) itmp,a(ii)%index,itmp2,a(ii)%x,a(ii)%y,a(ii)%z
	end do
	call loclabel(10,"Molecular Orbital Coefficients") !We first assume total is closed-shell, current is closed-shell too
	call readmatgau(10,tmpcomat,0,"f10.5",21,5,3) !i,j of tmpcomat is coefficient of i basis in j orbital
	if (iflipfrag(i)==0) then !fragment is closed-shell or needn't flip spin
		CObasa(locbas:locbas+numbas(i)-1,ioccmoa:ioccmoa+numaelec(i)-1)=tmpcomat(:,1:numaelec(i))
		CObasa(locbas:locbas+numbas(i)-1,ivirmoa:ivirmoa+(numbas(i)-numaelec(i))-1)=tmpcomat(:,numaelec(i)+1:)
	else if (iflipfrag(i)==1) then !flip spin, notice that now nbelec is already number of alpha electrons
		CObasb(locbas:locbas+numbas(i)-1,ioccmob:ioccmob+numbelec(i)-1)=tmpcomat(:,1:numbelec(i))
		CObasb(locbas:locbas+numbas(i)-1,ivirmob:ivirmob+(numbas(i)-numbelec(i))-1)=tmpcomat(:,numbelec(i)+1:)
	end if
	if (any(ifunrestrict==1)) then !Complex is open-shell, now fill beta part
		if (ifunrestrict(i)==0) then !fragment is closed-shell, beta coefficient is identical to alpha part
			CObasb(locbas:locbas+numbas(i)-1,ioccmob:ioccmob+numbelec(i)-1)=tmpcomat(:,1:numbelec(i))
			CObasb(locbas:locbas+numbas(i)-1,ivirmob:ivirmob+(numbas(i)-numbelec(i))-1)=tmpcomat(:,numbelec(i)+1:)
		else if (ifunrestrict(i)==1) then
			call loclabel(10,"Molecular Orbital Coefficients",ifound,0)
			call readmatgau(10,tmpcomat,0,"f10.5",21,5,3) !Load beta MOs information
			if (iflipfrag(i)==0) then
				CObasb(locbas:locbas+numbas(i)-1,ioccmob:ioccmob+numbelec(i)-1)=tmpcomat(:,1:numbelec(i))
				CObasb(locbas:locbas+numbas(i)-1,ivirmob:ivirmob+(numbas(i)-numbelec(i))-1)=tmpcomat(:,numbelec(i)+1:)
			else if (iflipfrag(i)==1) then !flipping spin
				CObasa(locbas:locbas+numbas(i)-1,ioccmoa:ioccmoa+numaelec(i)-1)=tmpcomat(:,1:numaelec(i))
				CObasa(locbas:locbas+numbas(i)-1,ivirmoa:ivirmoa+(numbas(i)-numaelec(i))-1)=tmpcomat(:,numaelec(i)+1:)
			end if
		end if
		wherefragb(ioccmob:ioccmob+numbelec(i)-1)=i
		wherefragb(ivirmob:ivirmob+(numbas(i)-numbelec(i))-1)=i
		ioccmob=ioccmob+numbelec(i)
		ivirmob=ivirmob+(numbas(i)-numbelec(i))
	end if
	
	wherefraga(ioccmoa:ioccmoa+numaelec(i)-1)=i
	wherefraga(ivirmoa:ivirmoa+(numbas(i)-numaelec(i))-1)=i
	ioccmoa=ioccmoa+numaelec(i) !Update ioccmoa
	ivirmoa=ivirmoa+(numbas(i)-numaelec(i)) !numbas(i)-numaelec(i) is number of virtual orbitals in current fragment

	deallocate(tmpcomat)
	iatm=iatm+numatom(i)
	locbas=locbas+numbas(i)
	close(10)
end do

MOocc=0D0
if (all(ifunrestrict==0)) then
	MOocc(1:nint(naelec))=2D0
else
	MOocc(1:nint(naelec))=1D0
	MOocc(nbasis+1:nbasis+nint(nbelec))=1D0
end if
!Generate new.gjf
open(10,file="new.gjf",status="replace")
write(10,"(a,/,/,a,/,/,2i3)") adjustl(trim(ctitle))//" guess=cards IOp(3/32=2)","Please check this file to ensure validity",nchargetot,itotmulti
do i=1,ncenter
	write(10,"(a,3f14.8)") ind2name(a(i)%index),a(i)%x,a(i)%y,a(i)%z
end do
write(10,"(/,'(5E16.9)',/,'-1')")
do i=1,nbasis !Cycle orbitals
	if (all(ifunrestrict==0)) then
		write(10,"('! Orbital:',i6,' Occ:',f10.6,' from fragment',i4)") i,MOocc(i),wherefraga(i)
	else
		write(10,"('! Alpha orbital:',i6,' Occ:',f10.6,' from fragment',i4)") i,MOocc(i),wherefraga(i)
	end if
	write(10,"(5E16.9)") (CObasa(j,i),j=1,nbasis)
end do
if (any(ifunrestrict==1)) then
	write(10,"('-1')")
	do i=1,nbasis
		write(10,"('! Beta orbital:',i6,' Occ:',f10.6,' from fragment',i4)") i,MOocc(nbasis+i),wherefragb(i)
		write(10,"(5E16.9)") (CObasb(j,i),j=1,nbasis)
	end do
end if
write(10,"('0',/)")
close(10)
write(*,*) "Input file with initial guess have been saved to new.gjf in current folder"
write(*,*) "Do not forget to manually check route section"
end subroutine



!!----------- Monitor SCF process
subroutine monitorscf
use defvar
use util
use dislin
implicit real*8 (a-h,o-z)
real*8,dimension(1000) :: DE,RMSDP,MaxDP,grad,stepnum(1000)=(/ (i,i=1,1000) /),constant
real*8 aimDE,aimRMSDP,aimMaxDP
character DEconv*3,RMSDPconv*3,MaxDPconv*3
itime=0
iscfqc=0
do while(.true.)
	if (itime/=0) then
		write(*,*)
		write(*,*) "1 Redraw convergence trend of all steps"
		write(*,*) "2 Redraw convergence trend of last 5 steps"
		write(*,*) "3 Redraw convergence trend of last 10 steps"
		write(*,*) "4 Redraw convergence trend of last specific number of steps"
		read(*,*) isel
		if (isel==4) then
			write(*,"(a,i5)") " Input a number, between 1 and",ifincyc
			read(*,*) ishow
		else if (isel==2.and.ifincyc>=5) then
			ishow=5
		else if (isel==3.and.ifincyc>=10) then
			ishow=10
		end if
	end if
	!Draw all points variation at first time get into this routine
	!I bracket filename with double quotation marks, otherwise if there have "+" in the filename then Multiwfn crash
	if (isys==1) call system("copy "//""""//trim(filename)//""""//" gauout.out /y")
	if (isys==2) call system("cp "//""""//trim(filename)//""""//" gauout.out -f")
	open(10,file="gauout.out",status="old")
	call loclabel(10," Quadratic Convergence",iscfqc)

	if (iscfqc==0) then !SCF=QC was not used
		call loclabel(10," Requested convergence on RMS density matrix")
		read(10,"(45x,D8.2)") aimRMSDP
		call loclabel(10," Requested convergence on MAX density matrix",ifound,0)
		read(10,"(45x,D8.2)") aimMaxDP
		call loclabel(10," Requested convergence on             energy",ifound,0)
		read(10,"(45x,D8.2)") aimDE

		i=1
		DE(1)=0D0 !In the first cycle DE is not present, we set it to an arbitrary value
		do while(.true.)
			call loclabel(10," RMSDP=",ifound,0)
			if (ifound==0) exit
			read(10,"(7x,1PD8.2)",advance="no") RMSDP(i)
			read(10,"(7x,1PD8.2)",advance="no") MaxDP(i)
			if (i/=1) read(10,"(4x,1PD9.2)") DE(i)
			i=i+1
		end do
		ifincyc=i-1
		if (ifincyc==0) then
			write(*,*) "Error: Unable to find SCF convergence information! Did you use #P as requested?"
			write(*,*) "Press ENTER button to exit program"
			stop
		end if
		if (isel==1.or.itime==0) ishow=ifincyc !At the first time, draw all point
		istart=ifincyc-ishow+1

		write(*,*) "Step#   RMSDP  Conv?   MaxDP  Conv?     DE    Conv?"
		do i=istart,ifincyc
			RMSDPconv="NO"
			MaxDPconv="NO"
			DEconv="NO"
			if (RMSDP(i)<aimRMSDP) RMSDPconv="YES"
			if (MaxDP(i)<aimMaxDP) MaxDPconv="YES"
			if (i==1) then
				DEconv="   "
			else if (abs(DE(i))<aimDE) then
				DEconv="YES"
			end if
			write(*,"(i5,2x,2(1PD8.2,a5,2x),1PD9.2,a5)") i,RMSDP(i),RMSDPconv,MaxDP(i),MaxDPconv,DE(i),DEconv
		end do
		write(*,"(' Goal  ',2(1PD8.2,7x),1PD9.2)") aimRMSDP,aimMaxDP,aimDE
	else if (iscfqc==1) then  !For SCF=QC
		i=1
		DE(1)=0D0
		do while(.true.)
			call loclabel(10,"Iteration ",ifound,0)
			if (ifound==0) exit
			read(10,"(50x)",advance="no")
			if (i/=1) read(10,*) DE(i)
			if (i==1) read(10,*)
			backspace(10)
			read(10,"(77x)",advance="no")
			read(10,*) grad(i)
			i=i+1
		end do
		ifincyc=i-1
		if (isel==1.or.itime==0) ishow=ifincyc
		istart=ifincyc-ishow+1

		write(*,*) "Step#       Delta-E         Grad"
		do i=istart,ifincyc
			write(*,"(i5,2x,2(1PE15.5))") i,DE(i),grad(i)
		end do
	end if
	call loclabel(10,"SCF Done",ifound)
	if (ifound==1) write(*,*) "SCF done!"
	if (ifound==0) write(*,*) "SCF failed or haven't converged"
	close(10)
    call delfile("gauout.out")

	call METAFL('xwin')
	call window(200,100,1000,700)
	call SCRMOD('REVERSE')
	CALL PAGE(3300,2310)
	call disini
	call hwfont
	call WINTIT("Monitor SCF process, variations at each step")
	CALL TICKS (1, 'XY')
	CALL TICPOS("REVERS","XYZ")
	call ERRMOD("ALL","OFF")
	CALL LABDIG(-1,"X")
	CALL LABDIG(2,"Y")
!  CALL AXSSCL('ELOG','Y')
	call labels('EXP','Y')
	stepx=ceiling((ifincyc-stepnum(istart))/20D0)
	if (iscfqc==0) then
		call AXSLEN(2300,600)
		!Draw RMSDP variation
		CALL NAME('RMS density matrix','Y')
		call AXSPOS(400,700)
		CALL GRAF(stepnum(istart),dfloat(ifincyc),stepnum(istart),stepx, 0D0,maxval(RMSDP(istart:ifincyc)),0D0,maxval(RMSDP(istart:ifincyc))/7)
		CALL CURVE(stepnum(istart:ifincyc),RMSDP(istart:ifincyc),ishow)
		call dash !Show where is goal
		call color('RED')
		constant=aimRMSDP
		CALL CURVE(stepnum(istart:ifincyc),constant(istart:ifincyc),ishow)
		call color('WHITE')
		call solid
		CALL ENDGRF
		!Draw MaxDP variation
		CALL NAME('MAX density matrix','Y')
		call AXSPOS(400,1400)
		CALL GRAF(stepnum(istart),dfloat(ifincyc),stepnum(istart),stepx, 0D0,maxval(MaxDP(istart:ifincyc)),0D0,maxval(MaxDP(istart:ifincyc))/7)
		CALL CURVE(stepnum(istart:ifincyc),MaxDP(istart:ifincyc),ishow)
		call dash !Show where is goal
		call color('RED')
		constant=aimMaxDP
		CALL CURVE(stepnum(istart:ifincyc),constant(istart:ifincyc),ishow)
		call color('WHITE')
		call solid
		CALL ENDGRF
		!Draw energy variation
	! 	CALL MARKER(21)
		CALL NAME('Energy','Y')
		call AXSPOS(400,2100)
		CALL NAME('Step number','X')
		fminDE=minval(DE(istart:ifincyc))
		fmaxDE=maxval(DE(istart:ifincyc))
		CALL GRAF(stepnum(istart),dfloat(ifincyc),stepnum(istart),stepx, fminDE,fmaxDE,fminDE,(fmaxDE-fminDE)/7)
		CALL CURVE(stepnum(istart:ifincyc),DE(istart:ifincyc),ishow)
		call dot !where is Y=0
		call XAXGIT
		call dash !Show where is goal
		call color('RED')
		constant=aimMaxDP
		CALL CURVE(stepnum(istart:ifincyc),constant(istart:ifincyc),ishow)
		constant=-aimMaxDP
		CALL CURVE(stepnum(istart:ifincyc),constant(istart:ifincyc),ishow)
		CALL ENDGRF
		! Show if have converged
		call color('RED')
		call height(50)
		if (RMSDP(ifincyc)<aimRMSDP) call messag("Done",2800,400)
		if (MaxDP(ifincyc)<aimMaxDP) call messag("Done",2800,1100)
		if (abs(DE(ifincyc))<aimDE) call messag("Done",2800,1800)
	else if (iscfqc==1) then
		call AXSLEN(2400,950)
		!Draw DE variation
		CALL NAME('Energy','Y')
		call AXSPOS(400,1050)
		fminDE=minval(DE(istart:ifincyc))
		fmaxDE=maxval(DE(istart:ifincyc))
		CALL GRAF(stepnum(istart),dfloat(ifincyc),stepnum(istart),stepx, fminDE,fmaxDE,fminDE,(fmaxDE-fminDE)/7)
		CALL CURVE(stepnum(istart:ifincyc),DE(istart:ifincyc),ishow)
		call endgrf
		!Draw gradient variation
		CALL NAME('Gradient','Y')
		call AXSPOS(400,2150)
		CALL NAME('Step number','X')
		CALL GRAF(stepnum(istart),dfloat(ifincyc),stepnum(istart),stepx, 0D0,maxval(grad(istart:ifincyc)),0D0,maxval(grad(istart:ifincyc))/7)
		CALL CURVE(stepnum(istart:ifincyc),grad(istart:ifincyc),ishow)
	end if
	call disfin
	itime=1
end do
end subroutine


!!------------- Show overlap matrix of alpha and beta orbitals
subroutine aboverlap
use util
use defvar
implicit real*8 (a-h,o-z)
character selectyn*1
real*8,allocatable :: ovlpmat(:,:)
real*8 :: GTFSintmat(nprims,nprims)

if (wfntype==0.or.wfntype==2.or.wfntype==3) then
	write(*,*) "Error: This function is only available for unrestricted wavefunction!"
	write(*,*)
	return
end if
write(*,*) "1 Calculate overlap between all alpha and all beta orbitals"
write(*,*) "2 Calculate overlap between alpha MOs and beta MOs with corresponding indices"
read(*,*) itype
do isplit=1,nmo
	if (MOtype(isplit)==2) exit
end do
numalphaMO=isplit-1
numbetaMO=nmo-isplit+1
allocate(ovlpmat(numalphaMO,numbetaMO))
ovlpmat=0D0
write(*,*) "Please wait..."
do ii=1,nprims
	do jj=1,nprims
		GTFSintmat(ii,jj)=doSint(ii,jj)
	end do
end do

if (itype==1) then
	!$OMP PARALLEL DO SHARED(ovlpmat) PRIVATE(iorba,iorbb,accum,ii,jj) schedule(dynamic) NUM_THREADS(nthreads)
	do iorba=1,numalphaMO !alpha orbitals
		do iorbb=isplit,nmo !beta orbitals
			accum=0D0
			if (allocated(CObasa)) then
				do ii=1,nbasis
					do jj=1,nbasis
						accum=accum+CObasa(ii,iorba)*CObasb(jj,iorbb-nbasis)*Sbas(ii,jj)
					end do
				end do
			else
				do ii=1,nprims
					do jj=1,nprims
						accum=accum+CO(iorba,ii)*CO(iorbb,jj)*GTFSintmat(ii,jj)
					end do
				end do
			end if
			ovlpmat(iorba,iorbb-isplit+1)=accum
		end do
	end do
	!$OMP end parallel do

	write(*,*)
	do iorb=1,numbetaMO
		write(*,"(' Overlap between the',i6,'th alpha and beta orbitals:',f12.6)") iorb,ovlpmat(iorb,iorb)
	end do

	if (wfntype==1) then !Calculate <S**2>
		tmp=0D0
		do iorba=1,naelec !Cycle occupied orbitals
			do iorbb=1,nbelec
				tmp=tmp+ovlpmat(iorba,iorbb)**2
			end do
		end do
	end if
	write(*,"(' <S**2> is ',f14.8)") (naelec-nbelec)/2*((naelec-nbelec)/2+1)+nbelec-tmp
	write(*,*)
	
	write(*,*) "Maximum pairing:"
	do i=1,numalphaMO
		imax=maxloc(abs(ovlpmat(i,:)),1)
		write(*,"(' Alpha',i6,'   Beta',i6,'   Overlap:',f12.6)") i,imax,abs(ovlpmat(i,imax))
	end do
	write(*,*)
	
	write(*,*) "If write overlap matrix to ovlpmat.txt in current folder? (y/n)"
	read(*,*) selectyn
	if (selectyn=='y'.or.selectyn=='Y') then
		open(10,file="ovlpmat.txt",status="replace")
		call showmatgau(ovlpmat,"Overlap matrix between alpha & beta orbitals",0,"1PE14.6",10)
		write(*,"(a)") " Done! The (i,j) element means the overlap integral between the ith alpha orbital and the jth beta orbital"
		close(10)
	end if
else if (itype==2) then
	do iorba=1,numbetaMO !Cycle all alpha orbitals, but the upper limit may be less than the beta ones, so use this limit
		accum=0D0
		iorbb=isplit-1+iorba
		do ii=1,nprims
			do jj=1,nprims
				accum=accum+CO(iorba,ii)*CO(iorbb,jj)*GTFSintmat(ii,jj)
			end do
		end do
		write(*,"(' Overlap between the ',i5,'th alpha and beta orbitals:',f12.6)") iorba,accum
	end do
end if
end subroutine



!!---------- Perform biorthogonalization between all alpha and beta orbitals
subroutine biortho
use defvar
character selectyn
real*8 Emat(nbasis,nbasis),tmparr(nbasis),orbsingval(nmo)

if (wfntype/=1) then
	write(*,*) "Error: This function is only available for unrestricted SCF wavefunction!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
if (.not.allocated(CObasa)) then
	write(*,*) "Error: This function requires basis function information!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

call ask_Sbas_PBC

!In fact, if choose option 1, the unoccupied alpha MO will keep unchanged, including its energy
!while beta MOs in the range of (naelec+1)~nbasis will be transformed at step 2 and energies are thus changed
write(*,*)
write(*,*) "Perform biorthogonalization for which orbitals?"
write(*,*) "1 Perform biorthogonalization only for occupied orbitals"
write(*,*) "2 Perform biorthogonalization for all orbitals"
write(*,*) "Note: Option 2 is time-consuming if there are large number of unoccupied MOs"
read(*,*) iorbrange

do while(.true.)
	write(*,*)
	write(*,*) "How to evaluate energies of biorthogonalized orbitals?"
	write(*,*) "0 Do not evaluate"
	write(*,*) "1 Evaluate, using the Fock matrix generated by MO energies and coefficients"
	write(*,*) "2 Evaluate, loading Fock matrix from a file"
	read(*,*) ievalene
	if (ievalene==0) then
		exit
	else if (ievalene==1) then
		call MOene2Fmat(istatus)
	else if (ievalene==2) then
		call loadFockfile(istatus)
	end if
	if (istatus==0) exit
end do

!Stage 1
call do_biortho(1,nint(naelec),1,nint(nbelec))
!Stage 2
if (nint(naelec)/=nint(nbelec)) then
	call do_biortho(nint(nbelec)+1,nint(naelec),nint(nbelec)+1,nbasis)
end if
!Stage 3
if (iorbrange==2) then
	call do_biortho(nint(naelec)+1,nbasis,nint(naelec)+1,nbasis)
else
	write(*,*) "Biorthogonalization between unoccupied orbitals is skipped as requested"
end if
orbsingval=MOene

!Evaluate orbital energies
if (ievalene>0) then
	write(*,*)
	write(*,*) "Evaluating energies of biorthogonalized orbitals..."
	Emat=matmul(matmul(transpose(CObasa),FmatA),CObasa)
	do iorb=1,nbasis
		MOene(iorb)=Emat(iorb,iorb)
	end do
	Emat=matmul(matmul(transpose(CObasb),FmatB),CObasb)
	do iorb=1,nbasis
		MOene(nbasis+iorb)=Emat(iorb,iorb)
	end do
	write(*,*) "Orbital energies have been successfully evaluated!"
	write(*,*)
	write(*,*) "Do you want to order each batch of orbitals according to their energies? (y/n)"
	read(*,*) selectyn
	if (selectyn=='y'.or.selectyn=='Y') then
		do itime=1,3
			if (itime==1) then
				nmobeg=1
				nmoend=nint(nbelec)
			else if (itime==2) then
				if (nint(naelec)==nint(nbelec)) cycle
				nmobeg=nint(nbelec)+1
				nmoend=nint(naelec)
			else if (itime==3) then
				if (nmoend==1) cycle
				nmobeg=nint(naelec)+1
				nmoend=nbasis
			end if
			do iorb=nmobeg,nmoend
				do jorb=iorb+1,nmoend
					avgene_i=(MOene(iorb)+MOene(iorb+nbasis))/2D0
					avgene_j=(MOene(jorb)+MOene(jorb+nbasis))/2D0
					if (avgene_i>avgene_j) then
						!alpha
						tmparr=CObasa(:,iorb)
						CObasa(:,iorb)=CObasa(:,jorb)
						CObasa(:,jorb)=tmparr
						tmpene=MOene(iorb)
						MOene(iorb)=MOene(jorb)
						MOene(jorb)=tmpene
						tmpsingval=orbsingval(iorb)
						orbsingval(iorb)=orbsingval(jorb)
						orbsingval(jorb)=tmpsingval
						!beta
						tmparr=CObasb(:,iorb)
						CObasb(:,iorb)=CObasb(:,jorb)
						CObasb(:,jorb)=tmparr
						tmpene=MOene(iorb+nbasis)
						MOene(iorb+nbasis)=MOene(jorb+nbasis)
						MOene(jorb+nbasis)=tmpene
						tmpsingval=orbsingval(iorb+nbasis)
						orbsingval(iorb+nbasis)=orbsingval(jorb+nbasis)
						orbsingval(jorb+nbasis)=tmpsingval
					end if
				end do
			end do
		end do
		write(*,*) "Ordering of biorthogonalized orbitals has finished!"
	end if
end if

!Output singular values and with/without energies in final sequence
write(*,*)
write(*,*) "Exporting biortho.txt..."
open(10,file="biortho.txt",status="replace")
nmoend=nbasis
if (iorbrange==1) nmoend=nint(naelec)
write(10,*) "S = Singular value, E = Energy (in eV), O= Occupancy, A=Alpha, B=Beta"
write(10,*)
if (ievalene>0) then
	do imo=1,nmoend
		write(10,"(' Orb:',i6,'   S=',f7.4,'  E(A)=',f11.3,'  O(A)=',f4.1,'  E(B)=',f11.3,'  O(B)=',f4.1)") &
		&imo,orbsingval(imo),MOene(imo)*au2eV,MOocc(imo),MOene(imo+nbasis)*au2eV,MOocc(imo+nbasis)
		if (imo==nint(nbelec).and.naelec/=nbelec) write(10,*) &
		"-------------------------------------------------------------------------------"
		if (iorbrange==2.and.imo==nint(naelec)) write(10,*) &
		"-------------------------------------------------------------------------------"
	end do
	write(*,"(a)") " Done! Singular values and energies of the biorthogonalized orbitals &
	&have been exported to biortho.txt in current folder"
else
	do imo=1,nmoend
		write(10,"(' Orb:',i6,'   S=',f7.4,'   O(A)=',f4.1,'   O(B)=',f4.1)") &
		imo,orbsingval(imo),MOocc(imo),MOocc(imo+nbasis)
		if (imo==nint(nbelec).and.naelec/=nbelec) &
		write(10,*) "-----------------------------------------------"
		if (iorbrange==2.and.imo==nint(naelec)) &
		write(10,*) "-----------------------------------------------"
	end do
	write(*,"(a)") " Done! Singular values of the biorthogonalized orbitals &
	&have been exported to biortho.txt in current folder"
end if
close(10)

write(*,*)
write(*,*) "Exporting biortho.fch..."
call outfch("biortho.fch",10,0)
write(*,"(a)",advance="no") " Done! biortho.fch has been exported in current folder, which contains biorthogonalized orbitals."
if (ievalene==0) then
	write(*,"(a)") " The orbital energies in a.u. correspond to singular values of orbital overlap matrix"
else
	write(*,"(a)") " The orbital energy information correspond to actual energy of the biorthogonalized orbitals"
end if
if (iorbrange==1) write(*,"(a,i6,a,i6,a)") " Note: Since biorthgonalization between unoccupied MOs was not performed, &
&the alpha and beta orbitals in the range",nint(naelec)+1," to",nbasis," are not meaningful"
write(*,*)
write(*,"(a)") " Do you want to load the biortho.fch now? If load, then you can directly visualize and analyze the biorthogonalized orbitals that just generated (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	write(*,*) "Loading biortho.fch..."
	call dealloall(0)
	call readinfile("biortho.fch",1)
	write(*,*) "Loading finished!"
else
	write(*,"(' Reloading ',a,'...')") trim(firstfilename)
	call dealloall(0)
	call readinfile(firstfilename,1)
end if
end subroutine
!!-------- Do biorthogonalization between specific range of alpha and beta orbitals
!Alpha range: na_beg~na_end; Beta range: nb_beg~nb_end
subroutine do_biortho(na_beg,na_end,nb_beg,nb_end)
use defvar
use util
implicit real*8 (a-h,o-z)
integer na_beg,na_end,nb_beg,nb_end
real*8,allocatable :: matU(:,:),matVT(:,:),matV(:,:),singval(:),ovlpmat(:,:)
! Since O = U * SIGMA * transpose(V) and U,V are unitary matrices, it clear that transpose(U) * O * V = SIGMA
write(*,"(/,' Doing biorthogonalization for alpha',i5,' to',i5,', Beta',i5,' to',i5,' ...')") na_beg,na_end,nb_beg,nb_end
na=na_end-na_beg+1
nb=nb_end-nb_beg+1
minnab=min(na,nb)
allocate(matU(na,na),matV(nb,nb),singval(minnab),ovlpmat(na,nb))
iprog=0
call showprog(iprog,na)
!$OMP PARALLEL DO SHARED(ovlpmat,iprog) PRIVATE(ia,imoa,ib,imob,ibas,jbas,tmpval) schedule(dynamic) NUM_THREADS(nthreads)
do ia=1,na
	imoa=na_beg+ia-1
	do ib=1,nb
		imob=nb_beg+ib-1
		tmpval=0
		do ibas=1,nbasis
			do jbas=1,nbasis
				tmpval=tmpval+CObasa(ibas,imoa)*CObasb(jbas,imob)*Sbas(ibas,jbas)
			end do
		end do
		ovlpmat(ia,ib)=tmpval
	end do
	!$OMP CRITICAL
	iprog=iprog+1
	call showprog(iprog,na)
	!$OMP end CRITICAL
end do
!$OMP END PARALLEL DO

call SVDmat(1,ovlpmat,matU,matV,singval,info)
if (info/=0) then
	write(*,*) "Error: SVD failed! The following result will be meaningless!"
	write(*,*) "Press ENTER button to continue"
	read(*,*)
	return
end if
write(*,*) "Singular values of orbital overlap matrix:"
write(*,"(8f9.4)") singval
CObasa(:,na_beg:na_end)=matmul(CObasa(:,na_beg:na_end),matU)
CObasb(:,nb_beg:nb_end)=matmul(CObasb(:,nb_beg:nb_end),matV)
MOene(na_beg:na_beg+minnab-1)=singval
MOene(nbasis+nb_beg:nbasis+nb_beg+minnab-1)=singval
end subroutine



!!---------- Integrating a function in whole space
!! ifunc: The real space function to be integrated
!The intval and funcval have 5 slots, the first one is used in normal case
!
!For isolated systems, atomic-center grids are used, &
!the integration grid is directly controlled by sphpot and radpot in settings.ini, since integrand may be not proportional to electron density,
!the grid will not be adjusted automatically (parm always=1) as proposed by Becke for more efficient integration of XC functional
!
!For periodic systems, evenly diestributed grids are used
subroutine intfunc(ifunc)
use functions
use util
implicit real*8 (a-h,o-z)
real*8 intval(5),intvalold(5),funcval(radpot*sphpot,5)
real*8 weigrid(radpot*sphpot) !Atom weighting function at grids
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 ELF2,ELF2r2 !Used to evaluate spherically symmetric average ELF (or LOL)

if (ifPBC==0) then !Isolated wavefunction
	call gen_GTFuniq(0) !Generate unique GTFs, for faster evaluation in orbderv
	if (outmedinfo==1) open(10,file="integrate.txt",status="replace")
	call walltime(iwalltime1)
	write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
	call gen1cintgrid(gridatmorg,iradcut)
	intval=0
	intvalold=0
	ELFsqr=0
	ELFsqrr2=0
	do iatm=1,ncenter
		write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		!$OMP parallel do shared(funcval) private(i,rnowx,rnowy,rnowz) num_threads(nthreads)
		do i=1+iradcut*sphpot,radpot*sphpot
			rnowx=gridatm(i)%x
			rnowy=gridatm(i)%y
			rnowz=gridatm(i)%z
			if (ispecial==0) then
				funcval(i,1)=calcfuncall(ifunc,rnowx,rnowy,rnowz) !This function automatically considers PBC
			else if (ispecial==1) then
				funcval(i,1)=infoentro(2,rnowx,rnowy,rnowz) !Shannon entropy density, see JCP,126,191107 for example
				funcval(i,2)=Fisherinfo(1,rnowx,rnowy,rnowz) !Fisher information density, see JCP,126,191107 for example
				funcval(i,3)=weizsacker(rnowx,rnowy,rnowz) !Steric energy
			end if
		end do
		!$OMP end parallel do
	
		!Generate atomic weighting function values at integration points
		!call gen1catmwei(iatm,iradcut,gridatm,weigrid,1) !Hirshfeld weightnig function
		call gen1cbeckewei(iatm,iradcut,gridatm,weigrid,covr_tianlu,3) !Becke weighting function
    
		do i=1+iradcut*sphpot,radpot*sphpot
			intval=intval+funcval(i,:)*weigrid(i)*gridatmorg(i)%value
			if (ifunc==9.or.ifunc==10) then !ELF and LOL
				tmp=funcval(i,1)**2*weigrid(i)*gridatmorg(i)%value
				r2=gridatm(i)%x**2+gridatm(i)%y**2+gridatm(i)%z**2
				ELFsqr=ELFsqr+tmp
				ELFsqrr2=ELFsqrr2+tmp*r2
			end if
 			if (outmedinfo==1) write(10,"(i7,3f12.5,3(1PE16.8))") i,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,funcval(i,1),gridatmorg(i)%value,weigrid(i)
		end do

		if (ispecial==0) write(*,"(' Accumulated value:',f20.10,'  Current center:',f20.10)") intval(1),intval(1)-intvalold(1)
		intvalold=intval
	end do

	call del_GTFuniq !Destory unique GTF informtaion
	call walltime(iwalltime2)
	write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

else !Periodic wavefunction
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
    write(*,*) "Calculating grid data..."
	!Because uniform grid cannot integrate well core density, so temporarily disable EDFs
    nEDFprims_org=nEDFprims
    nEDFprims=0
    call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
    call savecubmat(ifunc,0,1)
    call delvirorb_back(1) !Restore to previous wavefunction
    nEDFprims=nEDFprims_org
    call calc_dvol(dvol)
    intval(1)=sum(cubmat(:,:,:))*dvol
end if

write(*,*)
if (ispecial==0) then
	write(*,"(' Final result:',f20.10)") intval(1)
else if (ispecial==1) then
	write(*,"(' Shannon entropy:   ',f23.8)") intval(1)
	write(*,"(' Fisher information:',f23.8)") intval(2)
	write(*,"(' Steric energy:     ',f23.8)") intval(3)
end if
if (ifunc==9) then
	write(*,"(/,a,E14.6,' Bohr')") " int(ELF*r2*ELF):",ELFsqrr2
	write(*,"(a,E14.6)") " int(ELF*ELF):   ",ELFsqr
	write(*,"(a,f12.6,' Bohr')") " Spherically symmetric average ELF:",dsqrt(ELFsqrr2/ELFsqr)
else if (ifunc==10) then
	write(*,"(/,a,E14.6,' Bohr')") " int(LOL*r2*LOL):",ELFsqrr2
	write(*,"(a,E14.6)") " int(LOL*LOL):   ",ELFsqr
	write(*,"(a,f12.6,' Bohr')") " Spherically symmetric average LOL:",dsqrt(ELFsqrr2/ELFsqr)
end if

if (outmedinfo==1) then
	close(10)
    write(*,*) "integrate.txt has been exported to current folder"
    write(*,*) "Column 1: Index of integration points"
    write(*,*) "Columns 2~4: X, Y, Z of integration points in Bohr"
    write(*,*) "Column 5: Function value"
    write(*,*) "Column 6: Lebedev integration weighting"
    write(*,*) "Column 7: Atom weighting function"
end if
end subroutine


!!----------- Silent and simplied version of subroutine intfunc, directly return the integral. Currently used by EDA_SBL
subroutine intfunc_silent(ifunc,intval)
use functions
use util
implicit real*8 (a-h,o-z)
real*8 intval,funcval(radpot*sphpot),weigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)
intval=0
call showprog(0,ncenter)
do iatm=1,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	!$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
	do i=1+iradcut*sphpot,radpot*sphpot
		funcval(i)=calcfuncall(ifunc,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
	end do
	!$OMP end parallel do
	call gen1cbeckewei(iatm,iradcut,gridatm,weigrid,covr_tianlu,3)
	do i=1+iradcut*sphpot,radpot*sphpot
		intval=intval+funcval(i)*weigrid(i)*gridatmorg(i)%value
	end do
	call showprog(iatm,ncenter)
end do
end subroutine



!------ Integrate difference between two wavefunctions
!itype=1 : Integrate (f_wfn1 - f_wfn2)**2
!itype=2 : Integrate |f_wfn1 - f_wfn2|
!where f_wfn1 means real space function for wavefunction file 1, f_wfn2 is that for wavefunction file 2
!This function was specifically written for realizing Michael G. Medvedev's idea
subroutine intdiff(itype)
use functions
use util
implicit real*8 (a-h,o-z)
real*8 intval,intvalold,funcval1(radpot*sphpot),funcval2(radpot*sphpot),beckeweigrid(radpot*sphpot)
real*8 funcval1all(radpot*sphpot,ncenter) !For reuse data
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
character filename2*200,reusename*200

if (itype==1) write(*,*) "Note: The integrand is (f_wfn1 - f_wfn2)**2"
if (itype==2) write(*,*) "Note: The integrand is |f_wfn1 - f_wfn2|"
if (ifiletype/=2.and.ifiletype/=3) write(*,"(a)") " Hint: If you use .wfx or .wfn file as input instead, the calculation speed may be improved significantly!"
write(*,*)
write(*,*) "Select the function to be integrated over the whole space"
call selfunc_interface(1,ifunc)
write(*,*)
write(*,*) "Input the filename of another wavefunction file, e.g. C:\yuri.wfn"
read(*,"(a)") filename2
write(*,*)
write(*,"(a)") " Input density cutoff (e.g. 0.5), if electron density of the firstly loaded wavefunction at a point is &
&larger than this value, then corresponding integration grid will be ignored"
write(*,*) "If you do not want to enable this feature, input 0"
read(*,*) denscut

ireuse=0
write(reusename,"(a,'_',i3.3,'_',i4.4,'_',i4.4)") trim(firstfilename),ifunc,radpot,sphpot
inquire(file=reusename,exist=alive)
if (alive) then
	write(*,"(' Note: Data of reference system retrieved from ',a,' will be used')") trim(reusename)
	ireuse=1
	open(10,file=reusename,status="old")
	do iatm=1,ncenter
		read(10,*)
		read(10,*) funcval1all(:,iatm)
	end do
	close(10)
end if

write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)

call walltime(iwalltime1)

intval=0
do iatm=1,ncenter
	write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	
	if (ireuse==0) then !Calculate data for wfn1
		!$OMP parallel do shared(funcval1) private(i,rnowx,rnowy,rnowz) num_threads(nthreads) schedule(DYNAMIC)
		do i=1+iradcut*sphpot,radpot*sphpot
			rnowx=gridatm(i)%x
			rnowy=gridatm(i)%y
			rnowz=gridatm(i)%z
			funcval1(i)=calcfuncall(ifunc,rnowx,rnowy,rnowz)
		end do
		!$OMP end parallel do
		funcval1all(:,iatm)=funcval1
	else !Reuse data of wfn1 from external file
		funcval1=funcval1all(:,iatm)
	end if
	
	!Calculate data for wfn2
	call dealloall(0)
	call readinfile(filename2,1)
	!$OMP parallel do shared(funcval2) private(i,rnowx,rnowy,rnowz) num_threads(nthreads) schedule(DYNAMIC)
	do i=1+iradcut*sphpot,radpot*sphpot
		rnowx=gridatm(i)%x
		rnowy=gridatm(i)%y
		rnowz=gridatm(i)%z
		funcval2(i)=calcfuncall(ifunc,rnowx,rnowy,rnowz)
	end do
	!$OMP end parallel do
	!Recover to wfn1
	call dealloall(0)
	call readinfile(firstfilename,1)
	
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	do i=1+iradcut*sphpot,radpot*sphpot
		if (denscut/=0) then
			if (ifunc==1) then
				tmpdens=funcval1(i)
			else
				tmpdens=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
			end if
			if (tmpdens>denscut) cycle
		end if
		if (itype==1) then
			intval=intval+(funcval1(i)-funcval2(i))**2 *gridatmorg(i)%value*beckeweigrid(i)
		else
			intval=intval+abs(funcval1(i)-funcval2(i)) *gridatmorg(i)%value*beckeweigrid(i)
		end if
	end do
	write(*,"(' Accumulated value:',f20.10,'  Current center:',f20.10)") intval,intval-intvalold
	intvalold=intval
end do

call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,'s',/)") iwalltime2-iwalltime1
write(*,"(' Final result:',f24.12)") intval

!Write calculated wfn1 data to external file for reuse in the later calculation
if (ireuse==0) then
	open(10,file=reusename,status="replace")
	do iatm=1,ncenter
		write(10,*) iatm
		write(10,*) funcval1all(:,iatm)
	end do
	close(10)
	write(*,"(/,' Data of the firstly loaded file have been exported to ',a,' for possible later use')") trim(reusename)
end if
end subroutine




!------ Calculate overlap and distance between two orbitals
subroutine ovlpdistorb
use functions
use util
implicit real*8 (a-h,o-z)
real*8 intval(2),funcval(radpot*sphpot,2) !1/2=overlap of norm/sqr
real*8 cenpos(3,2),cenval(radpot*sphpot,3,2) !1st: point,  2st: x/y/z,  3st: iorb,jorb
real*8 beckeweigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
character selectyn

if (iautointgrid==1) then !This setting is good balance between cost and accuracy
	nradpotold=radpot
	nsphpotold=sphpot
	radcutold=radcut
	radpot=20
	sphpot=170
	radcut=18 !Enlarge radcut, because for Rydberg orbital, the default radcut 10 Bohr is not sufficient
end if

do while(.true.)
	write(*,*)
	write(*,*) "Input the index of the two orbitals, e.g. 32,35"
	write(*,*) "To exit, input 0,0"
	read(*,*) iorb,jorb
	if (iorb==0.and.jorb==0) then
		radpot=nradpotold
		sphpot=nsphpotold
		radcut=radcutold
		return
	end if
	
	if (iautointgrid==1) then
		write(*,"(a)") " Note: The default integration grids in general should be sufficient. If you want to change, &
		&set ""iautointgrid"" in settings.ini to 0, and set ""radpot"" and ""sphpot"" to expected values"
	end if
	write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
	call gen1cintgrid(gridatmorg,iradcut)

	intval=0
	cenpos=0
	do iatm=1,ncenter
        call showprog(iatm,ncenter)
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		!$OMP parallel do shared(funcval,cenval) private(i,vali,valj,vali2,valj2,rnowx,rnowy,rnowz) num_threads(nthreads)
		do i=1+iradcut*sphpot,radpot*sphpot
			rnowx=gridatm(i)%x
			rnowy=gridatm(i)%y
			rnowz=gridatm(i)%z
			vali=fmo(rnowx,rnowy,rnowz,iorb)
			vali2=vali**2 !rho of iorb
			valj=fmo(rnowx,rnowy,rnowz,jorb)
			valj2=valj**2 !rho of jorb
			funcval(i,1)=abs(vali)*abs(valj)
			funcval(i,2)=vali2*valj2
			cenval(i,1,1)=vali2*rnowx
			cenval(i,2,1)=vali2*rnowy
			cenval(i,3,1)=vali2*rnowz
			cenval(i,1,2)=valj2*rnowx
			cenval(i,2,2)=valj2*rnowy
			cenval(i,3,2)=valj2*rnowz
		end do
		!$OMP end parallel do
		
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
		do i=1+iradcut*sphpot,radpot*sphpot
			intval=intval+funcval(i,:)*gridatmorg(i)%value*beckeweigrid(i)
			cenpos=cenpos+cenval(i,:,:)*gridatmorg(i)%value*beckeweigrid(i)
		end do
		
	end do
	write(*,*)
	write(*,"(' X/Y/Z of centroid of electron density (Angstrom)')")
	write(*,"(' Orbital',i6,':',3f12.6)") iorb,cenpos(:,1)*b2a
	write(*,"(' Orbital',i6,':',3f12.6)") jorb,cenpos(:,2)*b2a
	write(*,"(' Centroid distance between the two orbitals:',f12.6,' Angstrom')") dsqrt(sum((cenpos(:,1)-cenpos(:,2))**2))*b2a
	write(*,"(' Overlap integral of norm of the two orbitals:',f16.10)") intval(1)
	write(*,"(' Overlap integral of square of the two orbitals:',f16.10)") intval(2)
	write(*,*)
	write(*,"(a)") " Do you want to add the two centroids as two dummy atoms? (y/n)"
	write(*,"(a)") " Note: Then you may enter main function 0 to visualize them along with corresponding orbital isosurfaces"
	read(*,*) selectyn
	if (selectyn=='y'.or.selectyn=='Y') then
		deallocate(a)
		allocate(a(ncenter+2))
		a(1:ncenter)=a_org
		a(ncenter+1)%name="Bq"
		a(ncenter+1)%index=0
		a(ncenter+1)%charge=0
		a(ncenter+1)%x=cenpos(1,1)
		a(ncenter+1)%y=cenpos(2,1)
		a(ncenter+1)%z=cenpos(3,1)
		a(ncenter+2)=a(ncenter+1)
		a(ncenter+2)%x=cenpos(1,2)
		a(ncenter+2)%y=cenpos(2,2)
		a(ncenter+2)%z=cenpos(3,2)
		ncenter=ncenter+2
		radpot=nradpotold
		sphpot=nsphpotold
		radcut=radcutold
		write(*,*) "Done!"
		return
	end if
end do
end subroutine



!!!------------ Interface of calculateing molecular volume by Monte Carlo method
subroutine molvol_MC
use defvar
implicit real*8 (a-h,o-z)
if (MCvolmethod==1) then !Using atomic vdW sphere superposition
	write(*,"(a)") " Note: The volume is defined as superposition of vdW spheres of atoms. 100*2^i points will be used to evaluate the volume by Monte Carlo method"
	do while(.true.)
		write(*,"(/,a)") " Please input i. Generally inputting 10 is recommended, however larger system requires larger i to ensure enough numerical accuracy"
		write(*,*) "Input 0 can return"
		read(*,*) pointexp
		if (pointexp==0) exit
		call calcvolume(1,pointexp,0D0,1D0)
	end do
else if (MCvolmethod==2) then !Using Bader's definition
    if (.not.allocated(b)) then
        write(*,"(a)") " Error: To calculate molecular volume based on isosurface of electron density, your input file must contain &
        &wavefunction information, namely you should use e.g. fch, molden, mwfn, wfn, wfx... file as input file"
        write(*,"(a)") " If you want to calculate molecular volume based on superposition of vdW spheres of atoms, please set ""MCvolmethod"" in settings.ini to 1"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        return
    end if
	write(*,*) "Note: 100*2^i points will be used to evaluate the volume by Monte Carlo method"
	write(*,"(a)") " The volume is defined as the region encompassed by the isosurface of electron density of x, &
	&The box used in the Monte Carlo calculation will be enlarged by vdW radii multiplied by k in each side"
	write(*,"(a)") " Hint: For evaluating the volume encompassed by rho=0.001 a.u. isosurface for small molecule, it is suggested to simply input 9,0.001,1.7"
	do while(.true.)
		write(*,*)
		write(*,*) "Please input i,x,k (input 0,0,0 can return)"
		read(*,*) pointexp,tmpisoval,enlarbox
		if (pointexp==0.and.tmpisoval==0.and.enlarbox==0) exit
		call calcvolume(2,pointexp,tmpisoval,enlarbox)
	end do
end if
end subroutine
!!!------------ Calculate molecular volume by Monte Carlo method
!If isoval == 0, then evaluate volume by approximate method, if >0, evaluate the volume based on rho
!imethod=1: Use superposition of vdW sphere to define vdW region, =2: use isosurface of electron density to define it
subroutine calcvolume(imethod,pointexp,mcisoval,enlarbox)
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 maxx,maxy,maxz,minx,miny,minz,lengthx,lengthy,lengthz,nowx,nowy,nowz,rx,ry,rz,r2,tmp
real*8 pointexp,mcisoval,enlarbox,boxvol
integer :: in,ntot,iatm,imethod,intmp
maxx=maxval( a(:)%x+enlarbox*vdwr(a(:)%index) )
maxy=maxval( a(:)%y+enlarbox*vdwr(a(:)%index) )
maxz=maxval( a(:)%z+enlarbox*vdwr(a(:)%index) )
minx=minval( a(:)%x-enlarbox*vdwr(a(:)%index) )
miny=minval( a(:)%y-enlarbox*vdwr(a(:)%index) )
minz=minval( a(:)%z-enlarbox*vdwr(a(:)%index) )
lengthx=maxx-minx
lengthy=maxy-miny
lengthz=maxz-minz
boxvol=lengthx*lengthy*lengthz
ntot=100*2**pointexp !More big more accurate

write(*,"(' Number of points used:',i12,', ',f10.3,' points per Bohr^3 in average')") ntot,ntot/(boxvol)
write(*,"(' Box size:',f12.3,' Bohr^3')") boxvol
write(*,*) "Please wait..."
in=0

if (imethod==1) then !I found if imethod=1 is parallelized too, the speed is much lowered!
	do i=1,ntot
		CALL RANDOM_NUMBER(nowx)
		CALL RANDOM_NUMBER(nowy)
		CALL RANDOM_NUMBER(nowz)
		nowx=nowx*lengthx+minx
		nowy=nowy*lengthy+miny
		nowz=nowz*lengthz+minz
		do iatm=1,ncenter
			rx=a(iatm)%x-nowx
			ry=a(iatm)%y-nowy
			rz=a(iatm)%z-nowz
			r2=rx*rx+ry*ry+rz*rz
			tmp=vdwr(a(iatm)%index)
			if (r2<=tmp*tmp) then
				in=in+1
				exit
			end if
		end do
	end do
else if (imethod==2) then
	!$OMP PARALLEL SHARED(in) PRIVATE(i,intmp,nowx,nowy,nowz) NUM_THREADS(nthreads)
	intmp=0
	!$OMP DO schedule(dynamic)
	do i=1,ntot
		CALL RANDOM_NUMBER(nowx)
		CALL RANDOM_NUMBER(nowy)
		CALL RANDOM_NUMBER(nowz)
		nowx=nowx*lengthx+minx
		nowy=nowy*lengthy+miny
		nowz=nowz*lengthz+minz
		if (fdens(nowx,nowy,nowz)>mcisoval) intmp=intmp+1
	end do
	!$OMP end do
	!$OMP CRITICAL
	in=in+intmp
	!$OMP end CRITICAL
	!$OMP END PARALLEL
end if

tmp=boxvol*(dfloat(in)/ntot)
write(*,"(' Molecular volume:',f10.3,' Bohr^3, (',f10.3,' Angstrom^3,',f9.3,' cm^3/mol)')") tmp,tmp*b2a**3,tmp*b2a**3*avogacst/1D24
end subroutine




!!!------------ Calculate LOLIPOP (LOL Integrated Pi Over Plane), see Chem. Commun., 48, 9239-9241 (2012)
!The selected plane is not necessarily parallel to XY plane, because grid data will be automatically projected
!------ Kawaii moe loli, pop, pop!
subroutine LOLIPOP
use defvar
use util
use GUI
implicit real*8 (a-h,o-z)
real*8 :: intradi=1.94D0,grdspc=0.08D0,LOLiso=0.55D0,disaway=0.5D0,vdwmulti=0.8D0
real*8 :: MOocc_old(nmo)
integer :: piorb(nmo),npiorb=0,ivisisosur=0
integer :: ringidx(100),nringidx=0
character c3000tmp*3000
integer :: ioutpt=0,iside=0
logical,allocatable :: gridconsider(:,:,:)

!Debug
!ivisisosur=1
!npiorb=7
!ioutpt=1
!piorb(1:7)=(/121,136,143,147,149,150,152/)
!grdspc=0.4D0

write(*,*) "## Kawaii moe loli, pop, pop!"
do while(.true.)
    write(*,*)
    write(*,*) "           =================  Calculate LOLIPOP  ================="
	write(*,*) "-1 Return"
	write(*,*) "0 Start calculation!"
	write(*,"(a,i7)") " 1 Choose pi orbitals, current number:",npiorb
	write(*,"(a,f8.4,a)") " 2 Set grid spacing, current:",grdspc," Bohr"
	write(*,"(a,f8.4,a)") " 3 Set integration radius, current:",intradi," Angstrom"
	write(*,"(a,f8.4,a)") " 4 Set the distance away the plane, current:",disaway," Angstrom"
    if (iside==0) write(*,*) "5 Choose side of the points to be taken into account, current: Both sides"
    if (iside==1) write(*,*) "5 Choose side of the points to be taken into account, current: Side 1"
    if (iside==2) write(*,*) "5 Choose side of the points to be taken into account, current: Side 2"
	if (ioutpt==0) write(*,*) "6 Toggle outputting actually considered points to pt.xyz, current: No"
	if (ioutpt==1) write(*,*) "6 Toggle outputting actually considered points to pt.xyz, current: Yes"
	if (ivisisosur==0) write(*,*) "7 Visualize LOL-pi isosurface after calculation, current: No"
	if (ivisisosur==1) write(*,*) "7 Visualize LOL-pi isosurface after calculation, current: Yes"
	read(*,*) isel

	if (isel==-1) then
		return
	else if (isel==1) then
		write(*,*) "Input indices of the pi orbitals, e.g. 17,20-25,36,37"
		read(*,"(a)") c3000tmp
		call str2arr(c3000tmp,npiorb,piorb)
	else if (isel==2) then
		write(*,*) "Input grid spacing in Bohr, e.g. 0.05"
		read(*,*) grdspc
	else if (isel==3) then
		write(*,*) "Input integration radius with respect to ring center in Angstrom, e.g. 1.94"
		read(*,*) intradi
	else if (isel==4) then
		write(*,*) "Input the distance away from the plane in Angstrom, e.g. 0.5"
		read(*,*) disaway
	else if (isel==5) then
        write(*,*) "Integrate LOL-pi in which side of the ring?"
        write(*,*) "0 Both sides"
        write(*,*) "1 Side 1"
        write(*,*) "2 Side 2"
        read(*,*) iside
    else if (isel==6) then
        if (ioutpt==1) then
            ioutpt=0
        else
            ioutpt=1
        end if
    else if (isel==7) then
        if (ivisisosur==0) then        
            ivisisosur=1
            sur_value=0.55D0
        else
            ivisisosur=0
        end if
		
	else if (isel==0) then
		if (npiorb==0) then
			write(*,*) "Error: You should use option 1 first to choose which orbitals are pi orbitals"
			write(*,*)
			cycle
		end if
		write(*,*) "Input indices of the atoms constituting the ring, e.g. 4,5,6,7,8,9"
		read(*,"(a)") c3000tmp
		call str2arr(c3000tmp,nringidx,ringidx)
		cenx=sum(a( ringidx(1:nringidx) )%x)/nringidx
		ceny=sum(a( ringidx(1:nringidx) )%y)/nringidx
		cenz=sum(a( ringidx(1:nringidx) )%z)/nringidx
		write(*,"(' Geometry center of the ring:',3f12.6,' Angstrom')") cenx*b2a,ceny*b2a,cenz*b2a
		endx=maxval( a(ringidx(1:nringidx))%x+vdwmulti*vdwr(a(ringidx(1:nringidx))%index) )
		endy=maxval( a(ringidx(1:nringidx))%y+vdwmulti*vdwr(a(ringidx(1:nringidx))%index) )
		endz=maxval( a(ringidx(1:nringidx))%z+vdwmulti*vdwr(a(ringidx(1:nringidx))%index) )
		orgx=minval( a(ringidx(1:nringidx))%x-vdwmulti*vdwr(a(ringidx(1:nringidx))%index) )
		orgy=minval( a(ringidx(1:nringidx))%y-vdwmulti*vdwr(a(ringidx(1:nringidx))%index) )
		orgz=minval( a(ringidx(1:nringidx))%z-vdwmulti*vdwr(a(ringidx(1:nringidx))%index) )
		dvol=grdspc**3
		write(*,"(' Spatial range of grid data to be calculated:')")
		write(*,"(' X is from',f10.4,'  to',f10.4,' Bohr')") orgx,endx
		write(*,"(' Y is from',f10.4,'  to',f10.4,' Bohr')") orgy,endy
		write(*,"(' Z is from',f10.4,'  to',f10.4,' Bohr')") orgz,endz
		write(*,"(' Differential element:',f12.6,' Bohr**3')") dvol
		xlength=endx-orgx
		ylength=endy-orgy
		zlength=endz-orgz
		dx=grdspc;gridv1=0;gridv1(1)=dx
		dy=grdspc;gridv2=0;gridv2(2)=dy
		dz=grdspc;gridv3=0;gridv3(3)=dz
		nx=nint(xlength/dx)+1
		ny=nint(ylength/dy)+1
		nz=nint(zlength/dz)+1
		if (allocated(cubmat)) deallocate(cubmat)
		allocate(cubmat(nx,ny,nz),gridconsider(nx,ny,nz))
        gridconsider=.false.
        
		write(*,"(' Number of points in x,y,z:',3i6,'  Total:',i10)") nx,ny,nz,nx*ny*nz
		write(*,*)
		write(*,"(' Pi orbitals:')")
		write(*,"(15i5)") piorb(1:npiorb)
        
		MOocc_old=MOocc !Backup
		do imo=1,nmo
			if (all(piorb(1:npiorb)/=imo)) MOocc(imo)=0D0 !Set occupation number of all orbitals to zero except for pi orbitals
		end do
        call delvirorb(0) !For saving time
		call savecubmat(10,0,1) !Calculate LOL
        call delvirorb_back(0)
        
		write(*,*)
		accum=0D0
		atm1x=a(ringidx(1))%x !Use 1,3,5 atoms in the ring to define the ring plane
		atm1y=a(ringidx(1))%y
		atm1z=a(ringidx(1))%z
		atm3x=a(ringidx(3))%x
		atm3y=a(ringidx(3))%y
		atm3z=a(ringidx(3))%z
		atm5x=a(ringidx(5))%x
		atm5y=a(ringidx(5))%y
		atm5z=a(ringidx(5))%z
		disple2crit=(disaway/b2a)**2
		discen2crit=(intradi/b2a)**2
        ncount=0
        
        !Construct two points for testing the current point is in which side by checking which one is closest to the current point
        call pointABCD(atm1x,atm1y,atm1z,atm3x,atm3y,atm3z,atm5x,atm5y,atm5z,pleA,pleB,pleC,pleD)
        valnorm=dsqrt(pleA**2+pleB**2+pleC**2)
        test1x=pleA/valnorm+cenx
        test1y=pleB/valnorm+ceny
        test1z=pleC/valnorm+cenz
        test2x=-pleA/valnorm+cenx
        test2y=-pleB/valnorm+ceny
        test2z=-pleC/valnorm+cenz
        
		do iz=1,nz
			do iy=1,ny
				do ix=1,nx
                    call getgridxyz(ix,iy,iz,xtmp,ytmp,ztmp)
					valtmp=cubmat(ix,iy,iz)
					call pointprjple(atm1x,atm1y,atm1z,atm3x,atm3y,atm3z,atm5x,atm5y,atm5z,xtmp,ytmp,ztmp,xprj,yprj,zprj) !Project grid point to the plane defined by atoms 1,3,5 in the ring
					if (valtmp>LOLiso) then
						disple2=(xtmp-xprj)**2+(ytmp-yprj)**2+(ztmp-zprj)**2 !The vertical distance**2 to the plane
						discen2=(xprj-cenx)**2+(yprj-ceny)**2+(zprj-cenz)**2 !The distance**2 of the projected point to ring center
						if (disple2>disple2crit.and.discen2<discen2crit) then
                            if (iside/=0) then
                                distest1=(xtmp-test1x)**2+(ytmp-test1y)**2+(ztmp-test1z)**2
                                distest2=(xtmp-test2x)**2+(ytmp-test2y)**2+(ztmp-test2z)**2
                                if ((iside==1.and.distest1>distest2).or.(iside==2.and.distest1<distest2)) cycle
                            end if
                            accum=accum+valtmp
                            ncount=ncount+1
                            gridconsider(ix,iy,iz)=.true.
                        end if
					end if
				end do
			end do
		end do
        
        write(*,"(' Number of points actually considered in the LOLIPOP integration:',i10)") count(gridconsider.eqv..true.)
        !Output points to pt.xyz for visual examination
        if (ioutpt==1) then
            open(10,file="pt.xyz",status="replace")
            write(10,*) ncount
            write(10,"(a,f12.6)") "Exported by Multiwfn, LOLIPOP value is",accum*dvol
		    do iz=1,nz
			    do iy=1,ny
				    do ix=1,nx
                        call getgridxyz(ix,iy,iz,xtmp,ytmp,ztmp)
					    if (gridconsider(ix,iy,iz)) write(10,"(a,3f12.6)") "C ",xtmp*b2a,ytmp*b2a,ztmp*b2a
				    end do
			    end do
		    end do
            close(10)
            write(*,"(a,/)") " Done! pt.xyz has been exported to current folder. You may use VMD program to load the file to visualize distribution of the points"
        end if
        
		write(*,"(' LOLIPOP value is',f12.6)") accum*dvol
		MOocc=MOocc_old !Recover original occupation number
        
        deallocate(gridconsider)
        if (ivisisosur==1) call drawisosurgui(1)
	end if
end do
end subroutine


!!---- Yoshizawa's electron transport route analysis
!Based on Eqs. 2 and 3 of Account of chemical research, 45, 1612
!Only closed-shell is supported
subroutine Yoshieletrans
use defvar
use util
use GUI
implicit real*8 (a-h,o-z)
character c80tmp*80,selectyn,c200tmp*200
real*8,allocatable :: NAOMO(:,:),mat(:,:),vallist(:)
integer,allocatable :: idx1(:),idx2(:) !Used to sort capacity
integer,allocatable :: atompi(:) !The ith element is the index of expected pi-AO of atom i
integer,allocatable :: piatm2atm(:) !Convert the index of the atom having pi-AO to absolute atom index
real*8 :: outcritval=0.01D0,outcritdistlow=0D0,outcritdisthigh=9999D0

!Load basic information about present system
if (ifiletype==0) then
    open(10,file=filename,status="old")
    call loclabel(10,"NBsUse=",ifound) !NbsUse must equal to the number of MOs
    if (ifound==0) then
	    write(*,"(a)") " Error: The input file you used does not meet requirement! Please carefully check Section 3.100.18 of the manual!"
	    write(*,*) "Press ENTER button to return"
	    read(*,*)
        close(10)
	    return
    end if
    read(10,*) c80tmp,nmo
    call loclabel(10,"NAtoms=",ifound)
    read(10,*) c80tmp,ncenter
    call loclabel(10,"alpha electrons",ifound)
    read(10,*) naelec,c80tmp,c80tmp,nbelec
    if (naelec/=nbelec) then
	    write(*,*) "Error: Only closed-shell wavefunction is supported!"
	    write(*,*) "Press ENTER button to return"
        read(*,*)
        close(10)
	    return
    end if
    !Load geometry
    if (allocated(a)) deallocate(a)
    allocate(a(ncenter))
    call loclabel(10,"Standard orientation:",ifound)
    if (ifound==0) then
	    write(*,*) "Error: Cannot found ""Standard orientation"" section!"
	    write(*,*) "Press ENTER button to return"
        read(*,*)
        close(10)
	    return
    end if
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    do iatm=1,ncenter !Assume that the coordinate is Angstrom
	    read(10,*) inouse,a(iatm)%index,inouse,a(iatm)%x,a(iatm)%y,a(iatm)%z
	    a(iatm)%name=ind2name(a(iatm)%index)
    end do
    a%x=a%x/b2a
    a%y=a%y/b2a
    a%z=a%z/b2a
    if (allocated(MOene)) deallocate(MOene)
    allocate(MOene(nmo))
    iHOMO=nint(naelec)
    iLUMO=iHOMO+1
    write(*,"(' Number of atoms:',i8)") ncenter
    write(*,"(' Total number of MOs:',i8)") nmo
    write(*,"(' Number of occupied MOs:',i8)") nint(naelec)
    write(*,"(' HOMO is orbital',i7,'           LUMO is orbital',i7)") iHOMO,iLUMO
    call loclabel(10,"occ. eigenvalues",ifound)
    do itime=1,ceiling(iHOMO/5D0) !Load occupied orbital energies
	    read(10,"(a)") c80tmp
	    ilow=(itime-1)*5+1
	    ihigh=itime*5
	    if (ihigh>=iHOMO) ihigh=iHOMO
	    read(c80tmp(29:),"(5f10.5)") MOene(ilow:ihigh)
    end do
    do itime=1,ceiling((nmo-iHOMO)/5D0) !Load unoccupied orbital energies
	    read(10,"(a)") c80tmp
	    ilow=(itime-1)*5+1+iHOMO
	    ihigh=itime*5+iHOMO
	    if (ihigh>=nmo) ihigh=nmo
	    read(c80tmp(29:),"(5f10.5)") MOene(ilow:ihigh)
    end do
     !write(*,*) "Energies of occupied MOs (a.u.):"
     !write(*,"(7f11.5)") MOene(:iHOMO)
     !write(*,*) "Energies of unoccupied MOs (a.u.):"
     !write(*,"(7f11.5)") MOene(iLUMO:)
else
    if (wfntype/=0) then
	    write(*,*) "Error: Only closed-shell wavefunction is supported!"
	    write(*,*) "Press ENTER button to return"
        read(*,*)
	    return
    end if
    write(*,"(a)") " Input file name containing information outputted by NBO program with NAOMO keyword, e.g. C:\test.out"
    do while(.true.)
	    read(*,"(a)") c200tmp
	    inquire(file=c200tmp,exist=alive)
	    if (alive) exit
	    write(*,*) "Cannot find the file, input again"
    end do
    open(10,file=c200tmp,status="old")
end if

write(*,*)
write(*,*) "The molecule is in which plane?  1=XY  2=YZ  3=XZ"
write(*,"(a)") " Note: This function cannot be used if all atoms are not in the same Cartesian plane"
read(*,*) iplesel

!Load information outputted by NBO program
allocate(atompi(ncenter))
atompi=0 !0 means the atom does not have corresponding pi orbital
call loclabel(10,"NATURAL POPULATIONS",ifound)
read(10,*)
read(10,*)
read(10,*)
read(10,*)
numNAO=0
do iatm=1,ncenter
	do while(.true.)
		read(10,"(a)") c80tmp
		if (c80tmp==" ") exit
		numNAO=numNAO+1
		if (iplesel==1.and.index(c80tmp,"Val")/=0.and.index(c80tmp,"pz")/=0) atompi(iatm)=numNAO
		if (iplesel==2.and.index(c80tmp,"Val")/=0.and.index(c80tmp,"px")/=0) atompi(iatm)=numNAO
		if (iplesel==3.and.index(c80tmp,"Val")/=0.and.index(c80tmp,"py")/=0) atompi(iatm)=numNAO
	end do
end do

write(*,"(' Number of natural atomic orbitals:',i8,/)") numNAO
numpiatom=count(atompi/=0)
write(*,"(' The number of atoms with expected pi atomic orbital:',i6)") numpiatom
numpair=numpiatom*(numpiatom-1)/2 !Number of pi-atom pairs
allocate(mat(numpiatom,numpiatom),piatm2atm(numpiatom))
itmp=0
do iatm=1,ncenter
	if (atompi(iatm)/=0) then
		write(*,"(' Atom:',i6,'       pi-NAO:',i6)") iatm,atompi(iatm)
		itmp=itmp+1
		piatm2atm(itmp)=iatm
	end if
end do

write(*,*) "Loading NAOMO matrix..."
allocate(NAOMO(numNAO,nmo))
call loclabel(10,"MOs in the NAO basis:",ifound,0) !Don't rewind
if (ifound==0) then
	write(*,"(a)") " Error: Cannot found ""MOs in NAO basis"" section! You should use ""NAOMO"" keyword in NBO module"
	write(*,*) "Press ENTER button to return"
    read(*,*)
    close(10)
	return
end if
!Check columns should be skipped during matrix reading, then return to title line
read(10,*);read(10,*);read(10,*)
read(10,"(a)") c80tmp
nskipcol=index(c80tmp,"- -")
backspace(10);backspace(10);backspace(10);backspace(10)
call readmatgau(10,NAOMO,0,"f8.4 ",nskipcol,8,3)

close(10)

iHOMO=nint(naelec)
iLUMO=iHOMO+1
Fermiene=(MOene(iHOMO)+MOene(iLUMO))/2D0
iorblow=1
iorbhigh=nmo
do while (.true.)
	write(*,*)
	write(*,*) "        ======= Yoshizawa's electron transport route analysis ======="
	write(*,*) "-10 Return"
	write(*,"(a,f10.4,a,f10.4,a)") " -4 Set distance criterion, current: From",outcritdistlow," to",outcritdisthigh," Angstrom"
	write(*,"(a,f12.8)") " -3 Set value criterion, current:",outcritval
	write(*,"(a,f12.6)") " -2 Set Fermi energy level, current (a.u.):",Fermiene
	write(*,"(a,i6,a,i6)") " -1 Select the range of MOs to be considered, current: from",iorblow,' to',iorbhigh
	write(*,*) "0 View molecular structure"
	write(*,*) "1 Output detail of electron transport probability between two atoms"
	write(*,*) "2 Output and rank all electron transport routes in the system"
	write(*,*) "3 Output and rank all electron transport routes for an atom"
	read(*,*) isel

	if (isel==-10) then
		return
	else if (isel==0) then
		call drawmolgui
	else if (isel==-1) then
		write(*,"(' Note: HOMO is orbital',i7,', LUMO is orbital',i7)") iHOMO,iLUMO
		write(*,*) "Input orbital range, e.g. 3,55 means all pi orbitals from orbital 3 to 55"
		write(*,*) "If input 0,0, then all MOs will be taken into account"
		read(*,*) iorblow,iorbhigh
		if (iorblow==0.and.iorbhigh==0) then
			iorblow=1
			iorbhigh=nmo
		else
			if (iorbhigh>nmo) iorbhigh=nmo
			if (iorblow<1) iorblow=1
		end if
	else if (isel==-2) then
		write(*,"(' Note: Energy of HOMO and LUMO is',2f12.6,' a.u., respectively')") MOene(iHOMO),MOene(iLUMO)
		write(*,*) "Input Fermi energy level in a.u., e.g. -0.005"
		read(*,*) Fermiene
	else if (isel==-3) then
		write(*,*) "Input the value criterion, e.g. 0.02"
		write(*,*) "Note: The routes whose value is smaller than this value will not be shown"
		read(*,*) outcritval
	else if (isel==-4) then
		write(*,*) "Input lower and upper limits of distance criterion (in Angstrom), e.g. 1.5,3.2"
		write(*,*) "Note: Only the routes within this range are possible to be shown"
		read(*,*) outcritdistlow,outcritdisthigh
	else if (isel==1) then
		write(*,*) "Input two atoms, e.g. 3,5"
		read(*,*) iatm,jatm
		transtot=0
		write(*,"(a,f10.6,a)") " Note: The MOs having contribution <",outcritval," will not be shown"
		write(*,*) "Note: HOMO is marked by asterisk"
		do iorb=iorblow,iorbhigh
			transtmp=NAOMO(atompi(iatm),iorb)*NAOMO(atompi(jatm),iorb)/(Fermiene-MOene(iorb))
			transtot=transtot+transtmp
			if (iorb==iHOMO) write(*,"('* MO:',i6,'   Energy(a.u.):',f12.6,'   Contribution:',f12.6)") iorb,MOene(iorb),transtmp
			if (iorb/=iHOMO.and.abs(transtmp)>=outcritval) write(*,"('  MO:',i6,'   Energy(a.u.):',f12.6,'   Contribution:',f12.6)") iorb,MOene(iorb),transtmp
! 			write(*,"(3f12.6,/)") NAOMO(atompi(iatm),iorb),NAOMO(atompi(jatm),iorb),Fermiene-MOene(iorb)
		end do
		write(*,"(' Total value is',f12.6)") transtot
		if (iorblow/=iHOMO.or.iorbhigh/=iLUMO) then
			contriHOMO=NAOMO(atompi(iatm),iHOMO)*NAOMO(atompi(jatm),iHOMO)/(Fermiene-MOene(iHOMO))
			contriLUMO=NAOMO(atompi(iatm),iLUMO)*NAOMO(atompi(jatm),iLUMO)/(Fermiene-MOene(iLUMO))
			write(*,"(/,' If only consider HOMO and LUMO, the value is',f12.6)") contriHOMO+contriLUMO
		end if
		write(*,"(' Distance of the route is',f12.6,' Angstrom')") dsqrt( (a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2 )*b2a
	else if (isel==2) then
		allocate(vallist(numpair),idx1(numpair),idx2(numpair))
		mat=0D0
		itmp=0
		do iatm=1,numpiatom !The index of pi-atom, must convert to actual atom index by piatm2atm
			do jatm=iatm+1,numpiatom
				itmp=itmp+1
				tmpval=0D0
				do iorb=iorblow,iorbhigh
					tmpval=tmpval+NAOMO(atompi(piatm2atm(iatm)),iorb)*NAOMO(atompi(piatm2atm(jatm)),iorb)/(Fermiene-MOene(iorb))
				end do
				mat(iatm,jatm)=tmpval
				vallist(itmp)=tmpval
				idx1(itmp)=iatm
				idx2(itmp)=jatm
			end do
		end do
		mat=mat+transpose(mat)
		!Sort
		do i=1,numpair
			do j=i+1,numpair
				if ( abs(vallist(i))<abs(vallist(j)) ) then
					temp=vallist(i)
					vallist(i)=vallist(j)
					vallist(j)=temp
					itemp=idx1(i)
					idx1(i)=idx1(j)
					idx1(j)=itemp
					itemp=idx2(i)
					idx2(i)=idx2(j)
					idx2(j)=itemp
				end if
			end do
		end do
		!Output ranked values
		write(*,*) "Electron transport route, ranked by transmission probability"
		write(*,"(' Note: The routes whose absolute value <',f10.6,' will not be shown')") outcritval
		write(*,"(' Note: The routes whose distance <',f10.4,' or >',f10.4,' Angstrom will not be shown')") outcritdistlow,outcritdisthigh
		do ipair=1,numpair
			if (abs(vallist(ipair))<outcritval) exit
			iatm=piatm2atm(idx1(ipair))
			jatm=piatm2atm(idx2(ipair))
			dist=dsqrt( (a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2 )*b2a
			if (dist>=outcritdistlow.and.dist<=outcritdisthigh) write(*,"(' Atom',i5,' -- Atom',i5,'  Value and distance:',2f12.6)") iatm,jatm,vallist(ipair),dist
		end do
		write(*,*) "Note: The units of the distances are Angstrom"
		write(*,*)
		write(*,*) "If outputting above data to result.txt in current folder? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			open(10,file="result.txt",status="replace")
			do ipair=1,numpair
				if (abs(vallist(ipair))<outcritval) exit
				iatm=piatm2atm(idx1(ipair))
				jatm=piatm2atm(idx2(ipair))
				dist=dsqrt( (a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2 )*b2a
				if (dist>=outcritdistlow.and.dist<=outcritdisthigh) write(10,"(2i6,2f12.6)") iatm,jatm,vallist(ipair),dist
			end do
			write(*,"(a)") " Done, the data has been outputted to result.txt in current folder. The first two columns correspond &
			&to atom indices, the third one corresponds to value, the fourth one is route distance (in Angstrom)."
			close(10)
		end if
		!Output matrix
		write(*,*)
		write(*,*) "If export electron transport probability matrix to current folder? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			open(10,file="transcapamat.txt",status="replace")
			call showmatgau(mat,"Electron transport probability matrix",0,"f14.8",10)
			close(10)
			write(*,*) "Done, the matrix has been outputted to transcapamat.txt in current folder"
			write(*,*) "Conversion between element in the matrix and actual atom index:"
			do iatm=1,numpiatom
				write(*,"(' Element',i6,'    -->   Atom',i6)") iatm,piatm2atm(iatm)
			end do
		end if
		deallocate(vallist,idx1,idx2)
	else if (isel==3) then
		write(*,*) "Input atom index, e.g. 5"
		read(*,*) iatm
		if (atompi(iatm)==0) then
			write(*,*) "Error: The atom does not have expected pi atomic orbital!"
			cycle
		end if
		allocate(vallist(numpiatom),idx1(numpiatom))
		do jatm=1,numpiatom !Calculate capacity to all other atoms
			tmpval=0D0
			do iorb=iorblow,iorbhigh
				tmpval=tmpval+NAOMO(atompi(piatm2atm(iatm)),iorb)*NAOMO(atompi(piatm2atm(jatm)),iorb)/(Fermiene-MOene(iorb))
			end do
			vallist(jatm)=tmpval
			idx1(jatm)=jatm
		end do
		vallist(iatm)=0D0 !The two sites are identical and meaningless
		!Sort
		do i=1,numpiatom
			do j=i+1,numpiatom
				if ( abs(vallist(i))<abs(vallist(j)) ) then
					temp=vallist(i)
					vallist(i)=vallist(j)
					vallist(j)=temp
					itemp=idx1(i)
					idx1(i)=idx1(j)
					idx1(j)=itemp
				end if
			end do
		end do
		write(*,*) "Electron transport route, ranked by transmission probability"
		write(*,"(' Note: The routes whose absolute value <',f10.6,' will not be shown')") outcritval
		write(*,"(' Note: The routes whose distance <',f10.4,' or >',f10.4,' Angstrom will not be shown')") outcritdistlow,outcritdisthigh
		do ipair=1,numpiatom
			jatm=piatm2atm(idx1(ipair))
			if (jatm==iatm) cycle
			if (abs(vallist(ipair))<outcritval) exit
			dist=dsqrt( (a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2 )*b2a
			if (dist>=outcritdistlow.and.dist<=outcritdisthigh) write(*,"(' To atom',i6,'    Value and distance (Angstrom):',2f12.6)") jatm,vallist(ipair),dist
		end do
		deallocate(vallist,idx1)
	end if
end do
end subroutine



!!----------- Detect pi orbital and set occupation number
subroutine detectpiorb
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
integer piorblist(nmo) !1 means this orbital is expected pi orbital
real*8,allocatable :: tmparr(:)
real*8 :: thresdens=0.01D0,thressingle=0.85D0
integer :: ionlyocc=1,idebug=0,icompmethod=1
integer,allocatable :: atmrange(:),tmpidx(:)
character c2000tmp*2000,c200tmp*200,selectyn
real*8 CObasa_LMO(nbasis,nbasis),CObasb_LMO(nbasis,nbasis),atmcomp(ncenter,nmo)

if (.not.allocated(b)) then
	write(*,*) "Error: Wavefunction information is not presented but needed!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
write(*,*)
write(*,"(a)") " Note: Please not only cite Multiwfn original paper but also cite below paper, in which the algorithm of this module &
&is very detailedly described and many interesting research examples are given:"
write(*,*) "Tian Lu, Qinxue Chen, Theor. Chem. Acc., 139, 25 (2020)"
write(*,*)
write(*,*) "Choose current situation:"
write(*,*) "-1: Orbitals are in localized form (e.g. LMO, NBO)"
write(*,"(a)") "  0: Orbitals are in delocalized form (e.g. MO, natural orbital, NTO). The system must be exactly planar"
read(*,*) iorbform

piorblist=0
pinelec=0D0
!tolerpara=0.1D0 !Too stringent, e.g. failed to detect all pi orbitals for cyclo[18]carbon under 0.029 a.u. field
tolerpara=0.15D0
tolerperp=80
if (iorbform==0) then !Delocalized case
    thres=0.05D0
    avgx=sum(a(:)%x)/ncenter
    avgy=sum(a(:)%y)/ncenter
    avgz=sum(a(:)%z)/ncenter
	if ( all(abs(a(:)%x-avgx)<thres) ) then
		iplane=2
		write(*,*) "This system is expected to be in YZ plane"
	else if ( all(abs(a(:)%y-avgy)<thres) ) then
		iplane=3
		write(*,*) "This system is expected to be in XZ plane"
	else if ( all(abs(a(:)%z-avgz)<thres) ) then
		iplane=1
		write(*,*) "This system is expected to be in XY plane"
        
	else !non-planar case
		write(*,"(a)") " Warning: Unable to detect the plane of the system! If you really want to carry out the pi-orbital detection, &
        &you need to manually choose an expected plane:"
        write(*,*) "0 Return"
        write(*,*) "1 XY plane"
        write(*,*) "2 YZ plane"
        write(*,*) "3 XZ plane"
        read(*,*) iplane
        if (iplane==0) return
        if (iplane==1) c200tmp="S, X and Y"
        if (iplane==2) c200tmp="S, Y and Z"
        if (iplane==3) c200tmp="S, X and Z"
        write(*,*) "Input a tolerance of expansion coefficient for "//trim(c200tmp)//" GTFs, e.g. 0.01"
        write(*,"(a)") " Note: If any above mentioned GTF has coefficient larger than this value, the orbital will not be identified as pi. &
        &Clearly, the larger the value, the higher the tendency that the orbitals will be determined as pi. &
        If you press ENTER button directly, 0.1 will be employed, which is usually suitable"
        read(*,"(a)") c200tmp
        if (c200tmp/=" ") read(c200tmp,*) tolerpara
        if (iplane==1) c200tmp="Z"
        if (iplane==2) c200tmp="X"
        if (iplane==3) c200tmp="Y"
        write(*,*) "Input a tolerance of percentage total contribution of "//trim(c200tmp)//" GTFs, e.g. 60"
        write(*,"(a)") " Note: If the contribution is lower than this value, the orbital will not be identified as pi. &
        &If you press ENTER button directly, 80% will be employed, which is usually suitable"
        read(*,"(a)") c200tmp
        if (c200tmp/=" ") read(c200tmp,*) tolerperp
        
	end if

    !SCPA method is used for calculating total contribution of P type of GTF that perpendicular to the plane (perpcontri)
	write(*,*) "Expected pi orbitals, occupation numbers and orbital energies (eV):"
	do imo=1,nmo
        perpcontri=0
		do iprim=1,nprims
			GTFtype=b(iprim)%type
			if (iplane==1) then !XY
				if ( (GTFtype==1.or.GTFtype==2.or.GTFtype==3).and.abs(CO(imo,iprim))>tolerpara ) exit !Orbital has S,X,Y component, so this is not pi-Z
                if (GTFtype==4.or.GTFtype==7.or.GTFtype==9.or.GTFtype==10) perpcontri=perpcontri+CO(imo,iprim)**2 !Z,ZZ,XZ,YZ
			else if (iplane==2) then !YZ
				if ( (GTFtype==1.or.GTFtype==3.or.GTFtype==4).and.abs(CO(imo,iprim))>tolerpara ) exit !Orbital has S,Y,Z component, so this is not pi-X
                if (GTFtype==2.or.GTFtype==5.or.GTFtype==8.or.GTFtype==9) perpcontri=perpcontri+CO(imo,iprim)**2 !X,XX,XY,XZ
                !if (imo==136) write(*,"(i5,2f12.6)") iprim,perpcontri,CO(imo,iprim)**2
			else if (iplane==3) then !XZ
				if ( (GTFtype==1.or.GTFtype==2.or.GTFtype==4).and.abs(CO(imo,iprim))>tolerpara ) exit !Orbital has S,X,Z component, so this is not pi-Y
                if (GTFtype==3.or.GTFtype==6.or.GTFtype==8.or.GTFtype==10) perpcontri=perpcontri+CO(imo,iprim)**2 !Y,YY,XY,YZ
			end if
			if (iprim==nprims) then
                testmag=sum(abs(CO(imo,:)))
                if (testmag<0.2D0) exit !The orbital may be core orbital of an atom, but GTF of this atom have been discarded using main function 6, therefore vanished
                !if (imo==136) write(*,*) perpcontri,sum(CO(imo,:)**2)
                perpcontri=perpcontri/sum(CO(imo,:)**2)*100 !Composition of 
                if (perpcontri<tolerperp) exit
				piorblist(imo)=1
				pinelec=pinelec+MOocc(imo)
				write(*,"(i6,2f14.6)") imo,MOocc(imo),MOene(imo)*au2ev
			end if
		end do
	end do
    isel=0
        
else if (iorbform==-1) then !LMO case
	if (.not.allocated(CObasa)) then
		write(*,"(a)") " Error: Basis function information is not presented but needed! See Section 2.5 of Multiwfn manual for detail"
		write(*,*) "Press ENTER button to return"
		read(*,*)
		return
	end if
	do while(.true.)
		write(*,*)
        write(*,"(a)") " -1 Detect pi orbitals and then evaluate pi composition for orbitals in another wavefunction file"
		write(*,*) "0 Detect pi orbitals and then set occupation numbers"
		write(*,"(a,f7.1,' %')") " 1 Set threshold for identifying single-center orbitals, current:",thressingle*100
		write(*,"(a,f8.4,' a.u.')") " 2 Set density threshold for identifying pi orbitals, current:",thresdens
		if (ionlyocc==1) write(*,*) "3 Switch the orbitals in consideration, current: Occupied localized orbitals"
		if (ionlyocc==0) write(*,*) "3 Switch the orbitals in consideration, current: All localized orbitals"
		if (idebug==1) write(*,*) "4 Switch outputting debug information, current: Yes"
		if (idebug==0) write(*,*) "4 Switch outputting debug information, current: No"
		if (.not.allocated(atmrange)) write(*,*) "5 Set constraint of atom range, current: undefined"
		if (allocated(atmrange)) write(*,"(a,i6,a)") " 5 Set constraint of atom range, current: ",natmrange," atoms"
        if (icompmethod==1) write(*,*) "6 Set the method for calculating orbital composition, current: Mulliken+SCPA"
        if (icompmethod==2) write(*,*) "6 Set the method for calculating orbital composition, current: Hirshfeld"
        if (icompmethod==3) write(*,*) "6 Set the method for calculating orbital composition, current: Becke"
		read(*,*) isel
		if (isel==0.or.isel==-1) then
			exit
		else if (isel==1) then
			write(*,*) "Input the threshold composition, e.g. 0.8"
			write(*,"(a)") " Note: If you input for example 0.8, then in an orbital, if an atom has contribution larger than 80%, &
			&then this orbital will be regarded as single-center orbital and will not be taken into account further"
			read(*,*) thressingle
		else if (isel==2) then
			write(*,*) "Input the threshold density in a.u., e.g. 0.02"
			write(*,"(a)") " Note: Assume that in an orbital, A and B are the two atoms having maximum contributions, &
			&the orbital will be regarded as pi orbital if its density at two representative points between A and B &
            &is both smaller than the threshold"
			read(*,*) thresdens
		else if (isel==3) then
			if (ionlyocc==1) then
				ionlyocc=0
			else if (ionlyocc==0) then
				ionlyocc=1
			end if
		else if (isel==4) then
			if (idebug==1) then
				idebug=0
			else if (idebug==0) then
				idebug=1
			end if
		else if (isel==5) then
			if (allocated(atmrange)) deallocate(atmrange)
			write(*,"(a)") " Input index range of the atoms, e.g. 2,3,7-10. &
			&Only LMOs with two largest contributing atoms belonging to this index range will be considered"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,natmrange)
			allocate(atmrange(natmrange))
			call str2arr(c2000tmp,natmrange,atmrange)
        else if (isel==6) then
            !Mulliken method is better than SCPA for LMOs, in particular for lone pair type of LMO, the SCPA shows it is too delocalized
            write(*,"(a,/)") " Hint: Option 1 is very fast, however the method is not robust, and it is even useless when diffuse functions are employed. &
            &Options 2 and 3 usually give similar result, they are more expensive, but fully compatible with diffuse functions"
            write(*,*) "1 Mulliken and SCPA methods for occupied and unoccupied orbitals, respectively"
            write(*,*) "2 Hirshfeld method"
            write(*,*) "3 Becke method"
            read(*,*) icompmethod
		end if
	end do
    
    !Evaluate compositions via Hirshfeld/Becke method, with cheap grid and silent mode
    if (icompmethod==2) call gen_orbatmcomp_space(1,atmcomp(:,:),1,nmo,0,0)
    if (icompmethod==3) call gen_orbatmcomp_space(2,atmcomp(:,:),1,nmo,0,0)
    write(*,*)
	write(*,*) "Expected pi orbitals, occupation numbers and orbital energies (eV):"
	allocate(tmparr(ncenter))
	do imo=1,nmo
		if (ionlyocc==1.and.MOocc(imo)==0) cycle
        !tmparr records composition of all atoms in current orbital
		if (icompmethod==1) then
            if (MOocc(imo)==0) then
                call gen_orbatmcomp_MMPA(1,imo,tmparr) !Mulliken method
            else
                call gen_orbatmcomp_MMPA(2,imo,tmparr) !SCPA method
            end if
        else
            tmparr=atmcomp(:,imo)
        end if
		!Find atom with maximum contribution (imax) and that with second maximum contribution (imax2)
		imax=maxloc(tmparr,1)
		tmpmax=-1D99
		do iatm=1,ncenter
			if (iatm==imax) cycle
			if (tmparr(iatm)>tmpmax) then
				imax2=iatm
				tmpmax=tmparr(iatm)
			end if
		end do
		!Use orbital density of at a few probe points between the atoms to determine if is pi orbital
        ratio1=0.7D0
        ratio2=0.3D0
        tmpx=ratio1*a(imax)%x+ratio2*a(imax2)%x
		tmpy=ratio1*a(imax)%y+ratio2*a(imax2)%y
		tmpz=ratio1*a(imax)%z+ratio2*a(imax2)%z
		dens1=fmo(tmpx,tmpy,tmpz,imo)**2
		tmpx=ratio2*a(imax)%x+ratio1*a(imax2)%x
		tmpy=ratio2*a(imax)%y+ratio1*a(imax2)%y
		tmpz=ratio2*a(imax)%z+ratio1*a(imax2)%z
		dens2=fmo(tmpx,tmpy,tmpz,imo)**2
		!tmpx=0.5D0*a(imax)%x+0.5D0*a(imax2)%x
		!tmpy=0.5D0*a(imax)%y+0.5D0*a(imax2)%y
		!tmpz=0.5D0*a(imax)%z+0.5D0*a(imax2)%z
		!dens3=fdens(tmpx,tmpy,tmpz)
		if (idebug==1) write(*,"(' Orb:',i5,'  max:',i5,f6.1,'%  max2:',i5,f6.1,'%  rho1:',f9.5,'  rho2:',f9.5)") &
        imo,imax,tmparr(imax)*100,imax2,tmparr(imax2)*100,dens1,dens2
		if (tmparr(imax)>thressingle) cycle !Pass single center LMO (Lone pair, inner-core)
		if (allocated(atmrange)) then
			if (all(atmrange/=imax).or.all(atmrange/=imax2)) cycle
		end if
		if (dens1<thresdens.and.dens2<thresdens) then
			piorblist(imo)=1
			pinelec=pinelec+MOocc(imo)
			write(*,"(i6,2f14.6)") imo,MOocc(imo),MOene(imo)*au2ev
		end if
	end do
	deallocate(tmparr)
end if

npiorb=count(piorblist==1)
if (npiorb==0) then
	write(*,*) "No pi orbital was found!"
	return
end if
write(*,"(' Total number of pi orbitals:',i6)") npiorb
write(*,"(' Total number of electrons in pi orbitals:',f12.6)") pinelec
if (iorbform==0.and.imodwfn==0) then !Only for MOs, one can safely separate inner and valence orbitals
	call getninnerele(ninnerele,0)
	ndelelec=ninnerele/2
	write(*,"(' Total number of inner-core electrons:',i6)") ninnerele
end if
if (npiorb/=0) then
    allocate(tmpidx(npiorb))
    nidx=0
    do imo=1,nmo
        if (piorblist(imo)==1.and.MOocc(imo)/=0) then
            nidx=nidx+1
            tmpidx(nidx)=imo
        end if
    end do
    call arr2str_2(tmpidx(1:nidx),c2000tmp)
    write(*,*) "Indices of occupied pi orbitals:"
    write(*,"(1x,a)") trim(c2000tmp)
    deallocate(tmpidx)
end if
write(*,*)

if (isel==0) then !Set occupation number
    write(*,*) "How to deal with these orbitals?"
    write(*,*) "0 Do nothing"
    write(*,*) "1 Set occupation number of these pi orbitals to zero"
    write(*,*) "2 Set occupation number of all other orbitals to zero"
    if (iorbform==0.and.imodwfn==0) then
	    write(*,*) "3 Set occupation number of valence pi orbitals to zero"
	    write(*,*) "4 Set occupation number of all except for valence pi orbitals to zero"
    end if
    read(*,*) isel

    if (isel/=0) then
	    if (isel==1) then
		    where (piorblist==1) MOocc=0
	    else if (isel==2) then
		    where (piorblist==0) MOocc=0
	    else if (isel==3.or.isel==4) then
		    if (wfntype==1.or.wfntype==4) then !UHF and U-post-HF wfn
			    do isplit=1,nmo !Where the first beta orbital appear now
				    if (motype(isplit)==2) exit
			    end do
			    do imo=1,nmo
				    if (isel==3) then
					    if (piorblist(imo)==1.and.imo<=isplit.and.imo>ndelelec) MOocc(imo)=0 !alpha part
					    if (piorblist(imo)==1.and.imo>isplit.and.imo>(isplit+ndelelec)) MOocc(imo)=0 !beta part
				    else if (isel==4) then
					    if (imo<=isplit) then !alpha part
						    if (piorblist(imo)==0.or.imo<=ndelelec) MOocc(imo)=0
					    else !beta part
						    if (piorblist(imo)==0.or.imo<=(isplit+ndelelec)) MOocc(imo)=0
					    end if
				    end if
			    end do
		    else if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Restricted(=0) or RO(=2) or post-R(=3) wavefunction	
			    do imo=1,nmo
				    if (isel==3.and.piorblist(imo)==1.and.imo>ndelelec) MOocc(imo)=0
				    if (isel==4.and.(piorblist(imo)==0.or.imo<=ndelelec)) MOocc(imo)=0
			    end do
		    end if
	    end if
	    imodwfn=1
	    call updatenelec
	    write(*,*) "Done!"
	    if (allocated(CObasa)) then
		    write(*,*) "Updating density matrix..."
		    call genP
		    write(*,*) "Density matrix has been updated"
	    end if
    end if
    write(*,*)

else if (isel==-1) then !Calculate pi composition for orbitals in another file
    CObasa_LMO=CObasa
    if (allocated(CObasb)) CObasb_LMO=CObasb
    call dealloall(0)
    write(*,*) "Input a file containing other set of orbitals, e.g. C:\riko.fch"
    do while(.true.)
	    read(*,"(a)") c200tmp
	    inquire(file=c200tmp,exist=alive)
	    if (alive) exit
	    write(*,*) "Cannot find the file, input again!"
    end do
    call readinfile(c200tmp,1)
    write(*,*) "Input threshold for printing (%), e.g. 85"
    write(*,*) "If pressing ENTER button directly, 50% will be employed as the threshold"
    read(*,"(a)") c200tmp
    if (c200tmp==" ") then
        thres=0.5D0
    else
        read(c200tmp,*) thres
        thres=thres/100
    end if
    if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Restricted orbitals
        do iorb=1,nmo
            if (ionlyocc==1.and.MOocc(iorb)==0) cycle !Only consider occupied orbitals
            picomp=0
            do ilmo=1,nmo
                if (piorblist(ilmo)==0) cycle !LMO is not pi
                coeff=sum(matmul(transpose(CObasa_lmo(:,ilmo:ilmo)),matmul(Sbas,CObasa(:,iorb:iorb))))
                picomp=picomp+coeff**2
            end do
            if (picomp>thres) write(*,"(' Orbital',i6,' (Occ=',f8.5')   pi composition:',f8.3,'%')") iorb,MOocc(iorb),picomp*100
        end do
    else if (wfntype==1.or.wfntype==4) then !Unrestricted orbitals
        write(*,*) "===== Alpha part ====="
        do iorb=1,nmo/2
            if (ionlyocc==1.and.MOocc(iorb)==0) cycle !Only consider occupied orbitals
            picomp=0
            do ilmo=1,nmo/2
                if (piorblist(ilmo)==0) cycle !LMO is not pi
                coeff=sum(matmul(transpose(CObasa_lmo(:,ilmo:ilmo)),matmul(Sbas,CObasa(:,iorb:iorb))))
                picomp=picomp+coeff**2
            end do
            if (picomp>thres) write(*,"(' Orbital',i6,' (Occ=',f8.5')   pi composition:',f8.3,'%')") iorb,MOocc(iorb),picomp*100
        end do
        write(*,*)
        write(*,*) "===== Beta part ====="
        do iorb=1,nmo/2
            if (ionlyocc==1.and.MOocc(iorb+nbasis)==0) cycle !Only consider occupied orbitals
            picomp=0
            do ilmo=1,nmo/2
                if (piorblist(ilmo+nbasis)==0) cycle !LMO is not pi
                coeff=sum(matmul(transpose(CObasb_lmo(:,ilmo:ilmo)),matmul(Sbas,CObasb(:,iorb:iorb))))
                picomp=picomp+coeff**2
            end do
            if (picomp>thres) write(*,"(' Orbital',i6,' (Occ=',f8.5')   pi composition:',f8.3,'%')") iorb,MOocc(iorb),picomp*100
        end do
    
    end if
    write(*,*)
    write(*,*) "Note: Current wavefunction corresponds to the file just loaded"
end if
end subroutine




!!------- A general routine used to fit atomic value from function value on vdW surface or on a set of given points
!Similar to routine "fitesp", but for universal purpose
subroutine fitfunc
use util
use defvar
use functions
implicit real*8 (a-h,o-z)
character addcenfile*200,extptfile*200
character selectyn,c80tmp
integer :: nlayer=4 !Number of fitting layers
real*8 :: funcfitvdwr(0:nelesupp)=-1D0,sclvdwlayer(100)=(/1.4D0,1.6D0,1.8D0,2D0,(2.2D0,i=5,100)/)
real*8,allocatable :: funcptval(:),funcptx(:),funcpty(:),funcptz(:),Bmat(:),Amat(:,:),Amatinv(:,:),atmval(:)
real*8,allocatable :: fitcenx(:),fitceny(:),fitcenz(:),fitcenvdwr(:),disptcen(:),origsphpt(:,:)
densperarea=6D0*b2a**2 !Point density per Angstrom**2 for MK, in order to convert to Bohr**2, multiply by b2a**2
iaddcen=0 !If give Additional center
iuseextpt=0 !If use external points
iskipfunccalc=0 !If read function value from external file directly rather than calculate here
iconstot=0

do while(.true.)
	write(*,*)
	write(*,*) "-2 Load additional fitting centers from external file"
	if (iuseextpt==0) write(*,*) "-1 Use fitting points recorded in external file instead of generating them"
	write(*,*) "0 Return"
	write(*,*) "1 Select a real space function and start calculation!"
	if (iuseextpt==0) then
		write(*,"(' 2 Set number of points per Angstrom^2, current:',f10.3)") densperarea/b2a**2 !Temporary convert to Angstrom**2 for convention
		write(*,"(' 3 Set number of layers per atom, current:',i4)") nlayer
		write(*,"(' 4 Set the scale factor of van der Waals radii in each layer')")
	end if
	if (iconstot==1) write(*,"(a,f14.7)") " 5 Set constraint for total value, current:",constotval
	if (iconstot==0) write(*,*) "5 Set constraint for total value, current: Not used"
	read(*,*) isel
	
	if (isel==-2) then
		iaddcen=1
		write(*,*) "Input the name of the file recording coordinates of additional fitting centers"
		read(*,"(a)") addcenfile
		write(*,*) "Done!"
	else if (isel==-1) then
		iuseextpt=1
		write(*,*) "Input the name of the file recording coordinates of fitting points"
		read(*,"(a)") extptfile
		write(*,*) "OK, the points recorded in this file will be used as fitting points"
	else if (isel==0) then
		Return
	else if (isel==1) then
		exit
	else if (isel==2) then
		write(*,*) "Input new value"
		read(*,*) densperarea
		densperarea=densperarea*b2a**2
	else if (isel==3) then
		write(*,*) "Input new value"
		read(*,*) nlayer
	else if (isel==4) then
		write(*,*) "Current values:"
		do ilayer=1,nlayer
			write(*,"(' Layer',i4,' :',f8.4)") ilayer,sclvdwlayer(ilayer)
		end do
		write(*,*)
		do ilayer=1,nlayer
			write(*,"(a,i4,',  e.g. 1.5')") " Input value for layer",ilayer
			read(*,*) sclvdwlayer(ilayer)
		end do
	else if (isel==5) then
		write(*,*) "Input a value, e.g. 3.0"
		write(*,*) "If input ""u"", constraint will not be applied to total value during fitting"
		read(*,"(a)") c80tmp
		if (c80tmp(1:1)=='u') then
			iconstot=0
		else
			iconstot=1
			read(c80tmp,*) constotval
		end if
	end if
end do

!Set vdW radius for MK, copied from GetvdW routine (utilam)
funcfitvdwr(1:17)=(/1.20d0,1.20d0,1.37d0,1.45d0,1.45d0,1.50d0,1.50d0,1.40d0,1.35d0,1.30d0,1.57d0,1.36d0,1.24d0,1.17d0,1.80d0,1.75d0,1.70d0/)
funcfitvdwr(1:17)=funcfitvdwr(1:17)/b2a
write(*,*) "Atomic radii used:"
do ielem=1,nelesupp
	if (any(a%index==ielem).and.funcfitvdwr(ielem)/=-1D0) write(*,"(' Element:',a,'     vdW radius (Angstrom):',f6.3)") ind2name(ielem),funcfitvdwr(ielem)*b2a
end do

!Check sanity and complete vdW radius table for all involved elements
do iatm=1,ncenter
	if (funcfitvdwr(a(iatm)%index)==-1D0) then
		write(*,"(' vdW radius used in fitting for element ',a,' is missing, input the radius (Bohr)')") ind2name(a(iatm)%index)
		write(*,"(a)") " Hint: If you do not know how to deal with the problem, simply input 3.4. (However, the radius of 3.4 Bohr may be not very appropriate for current element)" 
		read(*,*) funcfitvdwr(a(iatm)%index)
	end if
end do

!Check total number of fitting centers
naddcen=0
if (iaddcen==1) then
	open(10,file=addcenfile,status="old")
	read(10,*) naddcen
end if
nfitcen=ncenter+naddcen
allocate(fitcenx(nfitcen),fitceny(nfitcen),fitcenz(nfitcen),fitcenvdwr(nfitcen),disptcen(nfitcen))

!Generate information of fitting centers
do iatm=1,ncenter
	fitcenx(iatm)=a(iatm)%x
	fitceny(iatm)=a(iatm)%y
	fitcenz(iatm)=a(iatm)%z
	fitcenvdwr(iatm)=funcfitvdwr(a(iatm)%index) !vdW radius for each fitting center
end do
if (iaddcen==1) then
	do icen=ncenter+1,ncenter+naddcen
		read(10,*) fitcenx(icen),fitceny(icen),fitcenz(icen)
		fitcenvdwr(icen)=0D0
	end do
	close(10)
end if

write(*,*)
if (iuseextpt==0) then !Count number and generate coordinates of fitting points
	cutinnerscl=minval(sclvdwlayer(1:nlayer))
	write(*,"(' Note: If distance between a fitting point and any atom is smaller than',f6.3,' multiplied by corresponding vdW radius, then the point will be discarded')") cutinnerscl
	nfuncfitpt=0
	maxsphpt=nint(4D0*pi*(maxval(fitcenvdwr)*maxval(sclvdwlayer))**2 *densperarea) !Find maximal possible number of points in unit sphere to allocate temporary origsphpt
	allocate(origsphpt(3,maxsphpt))
	do icen=1,ncenter !Rather than nfitcen.   Count how many possible fitting points in total
		do ilayer=1,nlayer
			numsphpt=nint(4D0*pi*(fitcenvdwr(icen)*sclvdwlayer(ilayer))**2 *densperarea)
			nfuncfitpt=nfuncfitpt+numsphpt
		end do
	end do
	allocate(funcptval(nfuncfitpt),funcptx(nfuncfitpt),funcpty(nfuncfitpt),funcptz(nfuncfitpt)) !Currently nfuncfitpt is upper limit
	ifuncpt=0
	do icen=1,ncenter
		do ilayer=1,nlayer
			radius=fitcenvdwr(icen)*sclvdwlayer(ilayer)
			numsphpt=nint(4D0*pi*radius**2 *densperarea)
			call unitspherept(origsphpt,numsphpt) !Input expected number of point in unit sphere, return actual number of points
			origsphpt(:,1:numsphpt)=origsphpt(:,1:numsphpt)*radius
			origsphpt(1,1:numsphpt)=origsphpt(1,1:numsphpt)+fitcenx(icen) !Move unit sphere to atomic center
			origsphpt(2,1:numsphpt)=origsphpt(2,1:numsphpt)+fitceny(icen)
			origsphpt(3,1:numsphpt)=origsphpt(3,1:numsphpt)+fitcenz(icen)
			do ipt=1,numsphpt
				tmpx=origsphpt(1,ipt)
				tmpy=origsphpt(2,ipt)
				tmpz=origsphpt(3,ipt)
				iok=1
				do icen2=1,ncenter
					if (icen2==icen) cycle
					disptcensq=(fitcenx(icen2)-tmpx)**2+(fitceny(icen2)-tmpy)**2+(fitcenz(icen2)-tmpz)**2 !distance between point and center
					if (disptcensq<(fitcenvdwr(icen2)*cutinnerscl)**2) then !Less than vdW RADIUS*cutinner of atom icen2, it should be ommitted
						iok=0
						exit
					end if
				end do
				if (iok==1) then
					ifuncpt=ifuncpt+1
					funcptx(ifuncpt)=tmpx
					funcpty(ifuncpt)=tmpy
					funcptz(ifuncpt)=tmpz
				end if
			end do
		end do
	end do
	nfuncfitpt=ifuncpt
	deallocate(origsphpt)
	write(*,"(' Number of fitting points used:',i10)") nfuncfitpt
	
else if (iuseextpt==1) then !Directly use external fitting points
	open(10,file=extptfile,status="old")
	read(10,*) nfuncfitpt
	if (nfuncfitpt<0) then
		iskipfunccalc=1 !If the number of fitting points is negative, that means the fourth column records function value and needn't to be recalculated
		write(*,*) "Function value of all fitting points are read from external file directly"
	end if
	nfuncfitpt=abs(nfuncfitpt)
	write(*,"(' Number of fitting points used:',i10)") nfuncfitpt
	allocate(funcptval(nfuncfitpt),funcptx(nfuncfitpt),funcpty(nfuncfitpt),funcptz(nfuncfitpt))
	do i=1,nfuncfitpt
		if (iskipfunccalc==0) read(10,*) funcptx(i),funcpty(i),funcptz(i)
		if (iskipfunccalc==1) read(10,*) funcptx(i),funcpty(i),funcptz(i),funcptval(i)
	end do
	close(10)
end if

!Generate function value of fitting points
if (iskipfunccalc==0) then
	write(*,*) "Select the real space function to be fitted"
	call selfunc_interface(1,ifuncsel)
	write(*,*) "Calculating function value, please wait..."
    itmp=0
    call showprog(0,nfuncfitpt)
    !$OMP PARALLEL DO SHARED(funcptval,itmp) PRIVATE(ipt) schedule(dynamic) NUM_THREADS(nthreads)
	do ipt=1,nfuncfitpt
        if (ipt>=itmp*300) then
			call showprog(ipt,nfuncfitpt)
			itmp=itmp+1
		end if
		funcptval(ipt)=calcfuncall(ifuncsel,(funcptx(ipt)),(funcpty(ipt)),(funcptz(ipt)))
	end do
    !$OMP END PARALLEL DO
	call showprog(nfuncfitpt,nfuncfitpt)
	write(*,*) "Done!"
end if

matdim=nfitcen+1
allocate(Bmat(matdim),Amat(matdim,matdim),Amatinv(matdim,matdim),atmval(matdim))
Amat=0D0
do icen=1,nfitcen
	do jcen=icen,nfitcen
		do ipt=1,nfuncfitpt
			dis1=dsqrt( (funcptx(ipt)-fitcenx(icen))**2 + (funcpty(ipt)-fitceny(icen))**2 + (funcptz(ipt)-fitcenz(icen))**2 )
			dis2=dsqrt( (funcptx(ipt)-fitcenx(jcen))**2 + (funcpty(ipt)-fitceny(jcen))**2 + (funcptz(ipt)-fitcenz(jcen))**2 )
			Amat(icen,jcen)=Amat(icen,jcen)+1D0/dis1/dis2
		end do
	end do
end do
Amat=Amat+transpose(Amat)
do i=1,nfitcen
	Amat(i,i)=Amat(i,i)/2D0
end do
Amat(matdim,:)=1D0
Amat(:,matdim)=1D0
Amat(matdim,matdim)=0D0
Bmat=0D0
do icen=1,nfitcen
	do ipt=1,nfuncfitpt
		dis=dsqrt( (funcptx(ipt)-fitcenx(icen))**2 + (funcpty(ipt)-fitceny(icen))**2 + (funcptz(ipt)-fitcenz(icen))**2 )
		Bmat(icen)=Bmat(icen)+funcptval(ipt)/dis
	end do
end do
Bmat(matdim)=constotval !Constraint on the sum of all values
if (iconstot==1) then
	Amatinv=invmat(Amat,matdim)
	atmval=matmul(Amatinv,Bmat)
else
	Amatinv(1:nfitcen,1:nfitcen)=invmat(Amat(1:nfitcen,1:nfitcen),nfitcen)
	atmval(1:nfitcen)=matmul(Amatinv(1:nfitcen,1:nfitcen),Bmat(1:nfitcen))
end if

!Output summary
write(*,*) " Center       X           Y           Z             Value"
do i=1,ncenter
	write(*,"(i6,a,3f12.6,f16.6)") i,ind2name(a(i)%index),fitcenx(i),fitceny(i),fitcenz(i),atmval(i)
end do
do i=ncenter+1,ncenter+naddcen
	write(*,"(i6,2x,3f12.6,f16.6)") i,fitcenx(i),fitceny(i),fitcenz(i),atmval(i)
end do
write(*,"(' Sum of values:',f12.6)") sum(atmval(1:nfitcen))
!Calculate RMSE and RRMSE
RMSE=0D0
do ipt=1,nfuncfitpt
	atmvaleval=0D0 !Function value evaluated from atomic value by 1/r12
	do icen=1,nfitcen
		dis=dsqrt( (funcptx(ipt)-fitcenx(icen))**2 + (funcpty(ipt)-fitceny(icen))**2 + (funcptz(ipt)-fitcenz(icen))**2 )
		atmvaleval=atmvaleval+atmval(icen)/dis
	end do
	RMSE=RMSE+(funcptval(ipt)-atmvaleval)**2
end do
RRMSE=dsqrt(RMSE/sum(funcptval(1:nfuncfitpt)**2))
RMSE=dsqrt(RMSE/nfuncfitpt)
write(*,"(' RMSE:',f12.6,'   RRMSE:',f12.6)") RMSE,RRMSE

write(*,*)
write(*,"(a)") " If outputting coordinates and function value of all fitting points to funcfitpt.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=="Y") then
	open(10,file="funcfitpt.txt",status="replace")
	write(10,*) nfuncfitpt
	do ipt=1,nfuncfitpt
		write(10,"(3f13.7,E20.10)") funcptx(ipt),funcpty(ipt),funcptz(ipt),funcptval(ipt)
	end do
	write(*,*) "Data have been outputted to funcfitpt.txt in current folder"
	write(*,"(a)") " All units are in a.u. The first line shows the number of fitting points, &
	&the first three columns are X,Y,Z coordinates, the last column corresponds to function value"
	close(10)
end if
end subroutine



!!!------------------- Calculate properties based on atom geometry information
subroutine calcgeomprop
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer,allocatable :: calcatom(:)
do while(.true.)
	write(*,*) "Input indices of the atoms for which geometry information will be calculated"
	write(*,*) "e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will be considered"
	write(*,*) "Press ENTER button directly will analyze the whole system, input ""q"" will exit"
    write(*,*) "  Other commands:"
	write(*,*) "Input ""size"" will report size information of the whole system"
	write(*,*) "Input ""dist"" will report contact/distance between two specific fragments"
	write(*,*) "Input ""cav"" will report diameter of cavity enclosed by specific atoms"
	write(*,*) "Input ""ring"" will calculate area and perimeter of a specific ring"
	write(*,"(a)") " Input ""MPP"" will calculate molecular planarity parameter (MPP) and span of deviation from plane (SDP) for a fragment"
	read(*,"(a)") c2000tmp
	if (c2000tmp(1:1)=='q'.or.c2000tmp(1:1)=='Q') then
		exit
	else if (c2000tmp(1:4)=='size') then
		call calcmolsize
	else if (c2000tmp(1:3)=='cav') then
		call cavity_diameter
	else if (c2000tmp(1:4)=='dist') then
		call calcfragdist
	else if (c2000tmp(1:4)=='ring') then
		call calcringsize
	else if (index(c2000tmp,"MPP")/=0.or.index(c2000tmp,"mpp")/=0) then
		write(*,"(a,/)") " Hint: You can also directly enter this function by inputting ""MPP"" in main menu of Multiwfn"
		call calcMPP
	else
		if (allocated(calcatom)) deallocate(calcatom)
		if (c2000tmp==" ".or.index(c2000tmp,"all")/=0) then
			ncalcatom=ncenter
			allocate(calcatom(ncalcatom))
			do itmp=1,ncalcatom
				calcatom(itmp)=itmp
			end do
		else
			call str2arr(c2000tmp,ncalcatom)
			allocate(calcatom(ncalcatom))
			call str2arr(c2000tmp,ncalcatom,calcatom)
		end if
		call calcmolinfo(calcatom,ncalcatom)
	end if
	write(*,*)
end do
end subroutine
!!----- Show some molecular information based on geometry
!atmarray records which atoms will be taken into account, natmarr elements are there
subroutine calcmolinfo(atmarr,natmarr)
use util
use defvar
implicit real*8 (a-h,o-z)
integer atmarr(natmarr)
real*8 inertia(3,3),eigvalint(3),eigvecmatint(3,3)
totmass=sum(atmwei(a(atmarr(:))%index))
avgx=sum(a(atmarr(:))%x)/natmarr
avgy=sum(a(atmarr(:))%y)/natmarr
avgz=sum(a(atmarr(:))%z)/natmarr
cenmassx=sum(a(atmarr(:))%x*atmwei(a(atmarr(:))%index))/totmass
cenmassy=sum(a(atmarr(:))%y*atmwei(a(atmarr(:))%index))/totmass
cenmassz=sum(a(atmarr(:))%z*atmwei(a(atmarr(:))%index))/totmass
rgyr=dsqrt( sum( atmwei(a(atmarr(:))%index)* ((a(atmarr(:))%x-cenmassx)**2+(a(atmarr(:))%y-cenmassy)**2+(a(atmarr(:))%z-cenmassz)**2) ) / totmass )
totnucchg=sum(a(atmarr(:))%charge)
dipnucx=sum(a(atmarr(:))%x*a(atmarr(:))%charge)
dipnucy=sum(a(atmarr(:))%y*a(atmarr(:))%charge)
dipnucz=sum(a(atmarr(:))%z*a(atmarr(:))%charge)
dipnucnorm=dsqrt(dipnucx**2+dipnucy**2+dipnucz**2)
eleint=0D0
do iatmidx=1,natmarr
	iatm=atmarr(iatmidx)
    tmpval=0
	do jatmidx=iatmidx+1,natmarr
		jatm=atmarr(jatmidx)
        distval=atomdist(iatm,jatm,0)
		tmpval=tmpval+a(jatm)%charge/distval
	end do
    eleint=eleint+a(iatm)%charge*tmpval
    if (natmarr>20000) call showprog(iatmidx,natmarr)
end do
write(*,"(' Mass of these atoms:',f18.6,' amu')") totmass
write(*,"(' Geometry center (X/Y/Z):',3f14.8,' Angstrom')") avgx*b2a,avgy*b2a,avgz*b2a
write(*,"(' Center of mass (X/Y/Z): ',3f14.8,' Angstrom')") cenmassx*b2a,cenmassy*b2a,cenmassz*b2a
if (ifiletype==4) then !chg file
	write(*,"(' Sum of atomic charges:',f20.8)") totnucchg
	write(*,"(' Dipole from atomic charges (Norm): ',E12.5,' a.u.',E13.5,' Debye')") dipnucnorm,dipnucnorm*au2debye
	write(*,"(' Dipole from atomic charges (X/Y/Z):',3E12.5,' a.u.')") dipnucx,dipnucy,dipnucz
	write(*,"(' Electrostatic interaction energy between atomic charges:',/,f17.8,' a.u.',f20.5,' kcal/mol',f20.5,' KJ/mol')") eleint,eleint*au2kcal,eleint*au2KJ
else if (ifiletype/=4) then
	write(*,"(' Sum of nuclear charges:',f20.8)") totnucchg
	write(*,"(' Center of nuclear charges (X/Y/Z):',3f13.7,' Ang')") dipnucx*b2a/totnucchg,dipnucy*b2a/totnucchg,dipnucz*b2a/totnucchg
	write(*,"(' Dipole from nuclear charges (Norm):',E12.5,' a.u.',E13.5,' Debye')") dipnucnorm,dipnucnorm*au2debye
	write(*,"(' Dipole from nuclear charges (X/Y/Z):',3E12.5,' a.u.')") dipnucx,dipnucy,dipnucz
	write(*,"(' Electrostatic interaction energy between nuclear charges:',/,E18.10,' a.u.',E18.10,' kcal/mol',E18.10,' KJ/mol')") eleint,eleint*au2kcal,eleint*au2KJ
end if
write(*,"(' Radius of gyration:',f14.8,' Angstrom')") rgyr*b2a
xmin=a(atmarr(1))%x
ymin=a(atmarr(1))%y
zmin=a(atmarr(1))%z
xmax=xmin
ymax=ymin
zmax=zmin
ixmax=atmarr(1)
iymax=atmarr(1)
izmax=atmarr(1)
ixmin=atmarr(1)
iymin=atmarr(1)
izmin=atmarr(1)
do idx=1,natmarr
	iatm=atmarr(idx)
	if (a(iatm)%x>xmax) then
		xmax=a(iatm)%x
		ixmax=iatm
	end if
	if (a(iatm)%y>ymax) then
		ymax=a(iatm)%y
		iymax=iatm
	end if
	if (a(iatm)%z>zmax) then
		zmax=a(iatm)%z
		izmax=iatm
	end if
	if (a(iatm)%x<xmin) then
		xmin=a(iatm)%x
		ixmin=iatm
	end if
	if (a(iatm)%y<ymin) then
		ymin=a(iatm)%y
		iymin=iatm
	end if
	if (a(iatm)%z<zmin) then
		zmin=a(iatm)%z
		izmin=iatm
	end if
end do
write(*,"(' Minimum X is',f14.8,' Angstrom, at atom',i8,'(',a,')')") xmin*b2a,ixmin,a(ixmin)%name
write(*,"(' Minimum Y is',f14.8,' Angstrom, at atom',i8,'(',a,')')") ymin*b2a,iymin,a(iymin)%name
write(*,"(' Minimum Z is',f14.8,' Angstrom, at atom',i8,'(',a,')')") zmin*b2a,izmin,a(izmin)%name
write(*,"(' Maximum X is',f14.8,' Angstrom, at atom',i8,'(',a,')')") xmax*b2a,ixmax,a(ixmax)%name
write(*,"(' Maximum Y is',f14.8,' Angstrom, at atom',i8,'(',a,')')") ymax*b2a,iymax,a(iymax)%name
write(*,"(' Maximum Z is',f14.8,' Angstrom, at atom',i8,'(',a,')')") zmax*b2a,izmax,a(izmax)%name
if (natmarr>=2) then
	rmindist=1D50
	rmaxdist=0
	do iidx=1,natmarr
		i=atmarr(iidx)
		do jidx=iidx+1,natmarr
			j=atmarr(jidx)
            tmpval=atomdist(i,j,0)
			if (tmpval<rmindist) then
				rmindist=tmpval
				imindist=i
				jmindist=j
			end if
			if (tmpval>rmaxdist) then
				rmaxdist=tmpval
				imaxdist=i
				jmaxdist=j
			end if
		end do
        if (natmarr>20000) call showprog(iidx,natmarr)
	end do
	write(*,"(' Maximum distance is',f12.6,' Angstrom, between ',i8,'(',a,') and',i8,'(',a,')')") rmaxdist*b2a,imaxdist,a(imaxdist)%name,jmaxdist,a(jmaxdist)%name
	write(*,"(' Minimum distance is',f12.6,' Angstrom, between ',i8,'(',a,') and',i8,'(',a,')')") rmindist*b2a,imindist,a(imindist)%name,jmindist,a(jmindist)%name
end if
do iatmidx=1,natmarr
	iatm=atmarr(iatmidx)
	disttmp=dsqrt((avgx-a(iatm)%x)**2+(avgy-a(iatm)%y)**2+(avgz-a(iatm)%z)**2)
	if (iatmidx==1.or.disttmp>distmax) then
		distmax=disttmp
		idistmax=iatm
	end if
	if (iatmidx==1.or.disttmp<distmin) then
		distmin=disttmp
		idistmin=iatm
	end if
end do
write(*,"(' The atom closest to geometry center is',i9,'(',a,') Dist:',f12.6,' Angstrom')") idistmin,a(idistmin)%name,distmin*b2a
write(*,"(' The atom farthest to geometry center is',i8,'(',a,') Dist:',f12.6,' Angstrom')") idistmax,a(idistmax)%name,distmax*b2a
write(*,*)
inertia(1,1)=sum(atmwei(a(atmarr(:))%index)*( (a(atmarr(:))%y-cenmassy)**2+(a(atmarr(:))%z-cenmassz)**2) )*b2a*b2a
inertia(2,2)=sum(atmwei(a(atmarr(:))%index)*( (a(atmarr(:))%x-cenmassx)**2+(a(atmarr(:))%z-cenmassz)**2) )*b2a*b2a
inertia(3,3)=sum(atmwei(a(atmarr(:))%index)*( (a(atmarr(:))%x-cenmassx)**2+(a(atmarr(:))%y-cenmassy)**2) )*b2a*b2a
inertia(1,2)=-sum(atmwei(a(atmarr(:))%index)*(a(atmarr(:))%x-cenmassx)*(a(atmarr(:))%y-cenmassy))*b2a*b2a
inertia(2,1)=inertia(1,2)
inertia(1,3)=-sum(atmwei(a(atmarr(:))%index)*(a(atmarr(:))%x-cenmassx)*(a(atmarr(:))%z-cenmassz))*b2a*b2a
inertia(3,1)=inertia(1,3)
inertia(2,3)=-sum(atmwei(a(atmarr(:))%index)*(a(atmarr(:))%y-cenmassy)*(a(atmarr(:))%z-cenmassz))*b2a*b2a
inertia(3,2)=inertia(2,3)
call showmatgau(inertia,"Moments of inertia tensor (amu*Angstrom^2)")
write(*,"(' The moments of inertia relative to X,Y,Z axes (amu*Angstrom^2):',/,3E16.8)") inertia(1,1),inertia(2,2),inertia(3,3)
rotcstA=planckc/(8D0*pi**2*inertia(1,1)*amu2kg*1D-20)/1D9 !1D-20 used to convert Angstrom to meter, GHz=1D9Hz, 1Hz=1/s
rotcstB=planckc/(8D0*pi**2*inertia(2,2)*amu2kg*1D-20)/1D9
rotcstC=planckc/(8D0*pi**2*inertia(3,3)*amu2kg*1D-20)/1D9
write(*,"(' Rotational constant relative to X,Y,Z axes (GHz):',/,3f16.8)") rotcstA,rotcstB,rotcstC
call diagmat(inertia,eigvecmatint,eigvalint,300,1D-12)
call showmatgau(eigvecmatint,"Principal axes (each column vector)")
write(*,"(' The moments of inertia relative to principal axes (amu*Angstrom^2): ',/,3E16.8)") eigvalint
write(*,"(' Rotational constant relative to principal axes (GHz):',/,3f16.8)") planckc/(8D0*pi**2*eigvalint(1:3)*amu2kg*1D-20)/1D9
end subroutine

!!----- Calculate molecular diameter and length, width and height of the system
subroutine calcmolsize
use defvar
use util
use GUI
implicit real*8 (a-h,o-z)
real*8 inertia(3,3),eigvalarr(3),eigvecmat(3,3),vec3(1,3)

rmaxdist=0
do iatm=1,ncenter
	do jatm=iatm+1,ncenter
		tmpval=atomdist(iatm,jatm,0)
		if (tmpval>rmaxdist) then
			rmaxdist=tmpval
			imaxdist=iatm
			jmaxdist=jatm
		end if
	end do
end do
diameter=rmaxdist+vdwr(a(imaxdist)%index)+vdwr(a(jmaxdist)%index)
write(*,"(' Farthest distance:',i5,'(',a,')  ---',i5,'(',a,'):',f10.3, ' Angstrom')") imaxdist,a(imaxdist)%name,jmaxdist,a(jmaxdist)%name,rmaxdist*b2a
write(*,"(' vdW radius of',i5,'(',a,'):',f6.3,' Angstrom')") imaxdist,a(imaxdist)%name,vdwr(a(imaxdist)%index)*b2a
write(*,"(' vdW radius of',i5,'(',a,'):',f6.3,' Angstrom')") jmaxdist,a(jmaxdist)%name,vdwr(a(jmaxdist)%index)*b2a
write(*,"(' Diameter of the system:',f10.3,' Angstrom')") diameter*b2a
write(*,"(' Radius of the system:',f10.3,' Angstrom')") diameter*b2a/2
totmass=sum(atmwei(a%index))
cenmassx=sum(a%x*atmwei(a%index))/totmass
cenmassy=sum(a%y*atmwei(a%index))/totmass
cenmassz=sum(a%z*atmwei(a%index))/totmass
inertia(1,1)=sum(atmwei(a%index)*( (a%y-cenmassy)**2+(a%z-cenmassz)**2) )*b2a*b2a
inertia(2,2)=sum(atmwei(a%index)*( (a%x-cenmassx)**2+(a%z-cenmassz)**2) )*b2a*b2a
inertia(3,3)=sum(atmwei(a%index)*( (a%x-cenmassx)**2+(a%y-cenmassy)**2) )*b2a*b2a
inertia(1,2)=-sum(atmwei(a%index)*(a%x-cenmassx)*(a%y-cenmassy))*b2a*b2a
inertia(2,1)=inertia(1,2)
inertia(1,3)=-sum(atmwei(a%index)*(a%x-cenmassx)*(a%z-cenmassz))*b2a*b2a
inertia(3,1)=inertia(1,3)
inertia(2,3)=-sum(atmwei(a%index)*(a%y-cenmassy)*(a%z-cenmassz))*b2a*b2a
inertia(3,2)=inertia(2,3)
! call showmatgau(inertia,"Moments of inertia tensor (amu*Angstrom^2)")
call diagmat(inertia,eigvecmat,eigvalarr,300,1D-12)
! call showmatgau(eigvecmat,"Principal axes (each column vector)")
!Rotate the system to principal axis orientation
do iatm=1,ncenter
	vec3(1,1)=a(iatm)%x
	vec3(1,2)=a(iatm)%y
	vec3(1,3)=a(iatm)%z
	vec3(1,1:3)=matmul(vec3(1,1:3),eigvecmat)
	a(iatm)%x=vec3(1,1)
	a(iatm)%y=vec3(1,2)
	a(iatm)%z=vec3(1,3)
end do
!Make center of mass as origin
a%x=a%x-cenmassx
a%y=a%y-cenmassy
a%z=a%z-cenmassz
!Find minimum and maximum of the system, vdW radius has been taken into account
orgx=minval(a(:)%x-vdwr(a(:)%index))
endx=maxval(a(:)%x+vdwr(a(:)%index))
orgy=minval(a(:)%y-vdwr(a(:)%index))
endy=maxval(a(:)%y+vdwr(a(:)%index))
orgz=minval(a(:)%z-vdwr(a(:)%index))
endz=maxval(a(:)%z+vdwr(a(:)%index))
xlen=endx-orgx
ylen=endy-orgy
zlen=endz-orgz
write(*,"(' Length of the three sides:',3f10.3,' Angstrom')") xlen*b2a,ylen*b2a,zlen*b2a
do while(.true.)
	write(*,*)
	write(*,*) "0 Return"
	write(*,*) "1 Visualize the new orientation and molecular box"
	write(*,*) "2 Export the geometry in new orientation as new.pdb in current folder"
	read(*,*) isel
	if (isel==0) then
		a%x=a_org%x
		a%y=a_org%y
		a%z=a_org%z
		exit
	else if (isel==1) then
		ishowdatarange=1 !Draw grid data range
        !Define grid data information, use grid data framework to exhibit molecule size
        dtmp=0.25D0
        gridv1=0;gridv1(1)=dtmp
		gridv2=0;gridv2(2)=dtmp
		gridv3=0;gridv3(3)=dtmp
		nx=xlen/dtmp+1
		ny=ylen/dtmp+1
		nz=zlen/dtmp+1
		call miniGUI
		ishowdatarange=0
	else if (isel==2) then
		a%x=a%x-orgx !Temporarily make orgx,orgy,orgz as 0,0,0
		a%y=a%y-orgy
		a%z=a%z-orgz
		open(10,file="new.pdb",status="replace")
		write(10,"('REMARK   Generated by Multiwfn, Totally',i10,' atoms')") ncenter
		write(10,"('CRYST1',3f9.3,3f7.2,' P 1           1')") xlen*b2a,ylen*b2a,zlen*b2a,90D0,90D0,90D0
		do i=1,ncenter
			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
			"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,1.0,0.0,adjustr(ind2name_up(a(i)%index))
		end do
		close(10)
		write(*,*) "Exporting new.pdb file finished!"
		a%x=a%x+orgx
		a%y=a%y+orgy
		a%z=a%z+orgz
	end if
end do
end subroutine

!!--------- Calculate area for a given ring
subroutine calcringsize
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer,allocatable :: ringatmlist(:)
write(*,*) "Input the index of the atoms in the ring in clockwise manner, e.g. 2,3,4,6,7"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,nringatm)
allocate(ringatmlist(nringatm))
call str2arr(c2000tmp,nringatm,ringatmlist)
ringarea=0D0
ia=2
ib=nringatm
itri=0
do while(.true.)
	ntime=2
    if (ib-ia==1) ntime=1
    iatm1=ringatmlist(ia)
    iatm2=ringatmlist(ib)
    do itmp=1,ntime
        if (itmp==1) ic=ia-1
        if (itmp==2) ic=ib-1
        iatm3=ringatmlist(ic)
        itri=itri+1
        triarea=gettriangarea(a(iatm1)%x,a(iatm1)%y,a(iatm1)%z,a(iatm2)%x,a(iatm2)%y,a(iatm2)%z,a(iatm3)%x,a(iatm3)%y,a(iatm3)%z)*b2a**2
        ringarea=ringarea+triarea
        write(*,"(' Atoms in triangle',i4,':',3i7,'      Area:',f10.5,' Angstrom^2')") itri,iatm1,iatm2,iatm3,triarea
    end do
    ia=ia+1
    ib=ib-1
    if (ib<=ia) exit
end do
write(*,"(/,' The total ring area is',f12.6,' Angstrom^2')") ringarea
ringperi=0D0
do i=1,nringatm
    iatm1=ringatmlist(i)
    if (i==nringatm) then
		iatm2=ringatmlist(1)
	else
		iatm2=ringatmlist(i+1)
	end if
	ringperi=ringperi+dsqrt((a(iatm1)%x-a(iatm2)%x)**2+(a(iatm1)%y-a(iatm2)%y)**2+(a(iatm1)%z-a(iatm2)%z)**2)*b2a
end do
write(*,"(' The ring perimeter is',f12.6,' Angstrom',/)") ringperi
end subroutine

!!------ Calculate distance between two fragments
subroutine calcfragdist
use util
use defvar
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer fr1(ncenter),fr2(ncenter)
write(*,*) "Input atom indices for fragment 1, e.g. 3,5-8,15-20"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,nfr1,fr1)
write(*,*) "Input atom indices for fragment 2, e.g. 3,5-8,15-20"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,nfr2,fr2)

distmin=1D10
distmax=0
do idx=1,nfr1
    iatm=fr1(idx)
    do jdx=1,nfr2
        jatm=fr2(jdx)
        dist=atomdist(iatm,jatm,0)
        if (dist<distmin) then
            distmin=dist
            minatm1=iatm
            minatm2=jatm
        end if
        if (dist>distmax) then
            distmax=dist
            maxatm1=iatm
            maxatm2=jatm
        end if
    end do
end do

write(*,"(' Minimum distance:',f10.4,' Angstrom, between',i6,'(',a,') and',i6,'(',a,')')") distmin*b2a,minatm1,a(minatm1)%name,minatm2,a(minatm2)%name
write(*,"(' Maximum distance:',f10.4,' Angstrom, between',i6,'(',a,') and',i6,'(',a,')')") distmax*b2a,maxatm1,a(maxatm1)%name,maxatm2,a(maxatm2)%name

gc1x=sum(a(fr1(1:nfr1))%x)/nfr1
gc1y=sum(a(fr1(1:nfr1))%y)/nfr1
gc1z=sum(a(fr1(1:nfr1))%z)/nfr1
gc2x=sum(a(fr2(1:nfr2))%x)/nfr2
gc2y=sum(a(fr2(1:nfr2))%y)/nfr2
gc2z=sum(a(fr2(1:nfr2))%z)/nfr2
write(*,"(' Geometry center of fragment 1 (X/Y/Z):',3f10.4,' Angstrom')") gc1x*b2a,gc1y*b2a,gc1z*b2a
write(*,"(' Geometry center of fragment 2 (X/Y/Z):',3f10.4,' Angstrom')") gc2x*b2a,gc2y*b2a,gc2z*b2a
write(*,"(' Distance between the two geometry centers:',f10.4,' Angstrom')") dsqrt((gc1x-gc2x)**2+(gc1y-gc2y)**2+(gc1z-gc2z)**2)*b2a
fr1mass=sum(atmwei(a(fr1(1:nfr1))%index))
cm1x=sum(a(fr1(1:nfr1))%x*atmwei(a(fr1(1:nfr1))%index))/fr1mass
cm1y=sum(a(fr1(1:nfr1))%y*atmwei(a(fr1(1:nfr1))%index))/fr1mass
cm1z=sum(a(fr1(1:nfr1))%z*atmwei(a(fr1(1:nfr1))%index))/fr1mass
fr2mass=sum(atmwei(a(fr2(1:nfr2))%index))
cm2x=sum(a(fr2(1:nfr2))%x*atmwei(a(fr2(1:nfr2))%index))/fr2mass
cm2y=sum(a(fr2(1:nfr2))%y*atmwei(a(fr2(1:nfr2))%index))/fr2mass
cm2z=sum(a(fr2(1:nfr2))%z*atmwei(a(fr2(1:nfr2))%index))/fr2mass
write(*,"(' Mass center of fragment 1 (X/Y/Z):',3f10.4,' Angstrom')") cm1x*b2a,cm1y*b2a,cm1z*b2a
write(*,"(' Mass center of fragment 2 (X/Y/Z):',3f10.4,' Angstrom')") cm2x*b2a,cm2y*b2a,cm2z*b2a
write(*,"(' Distance between the two mass centers:',f10.4,' Angstrom')") dsqrt((cm1x-cm2x)**2+(cm1y-cm2y)**2+(cm1z-cm2z)**2)*b2a
end subroutine



!!------------ Calculate molecular planarity parameter (MPP) and span of deviation from plane (SDP)
subroutine calcMPP
use defvar
use util
implicit real*8 (a-h,o-z)
integer,allocatable :: atmsel(:),frameidx(:)
character c2000tmp*2000,c200tmp*200,outname*200,selectyn
real*8 MPP
write(*,"(a)") " If this function is employed in your work, please cite:"
write(*,"(a)") " Tian Lu, Simple, reliable, and universal metrics of molecular planarity, J. Mol. Model., 27, 263 (2021) DOI: 10.1007/s00894-021-04884-0"
write(*,*)

write(*,*) "Input index of the atoms for which data will be calculated, e.g. 2,3,7-10"
write(*,*) "Input ""a"" can choose all atoms, input ""h"" can choose all heavy atoms"
read(*,"(a)") c2000tmp
if (index(c2000tmp,'a')/=0) then
	nsel=ncenter
	allocate(atmsel(nsel))
	forall(i=1:nsel) atmsel(i)=i
else if (index(c2000tmp,'h')/=0) then
	nsel=count(a%index/=1)
	allocate(atmsel(nsel))
    nsel=0
    do iatm=1,ncenter
		if (a(iatm)%index==1) cycle
		nsel=nsel+1
        atmsel(nsel)=iatm
    end do
else
	call str2arr(c2000tmp,nsel)
	allocate(atmsel(nsel))
	call str2arr(c2000tmp,nsel,atmsel)
end if

!May be xyz trajectory file, use special code
itmp=len_trim(filename)
if (filename(itmp-3:itmp)==".xyz") then 
	!Test if is a .xyz trajectory file
	open(10,file=filename,status="old")
    do i=1,ncenter+2
		read(10,*)
    end do
    read(10,*,iostat=ierror) inouse
    close(10)
    if (ierror==0) then
		write(*,"(a)") " The loaded .xyz file seems to be a trajectory file. If you want to calculate MPP, SDP and atomic ds values &
        &for a range of frames, input the frame indices, e.g. 1,5-10,22-30. Inputting ""a"" can select all frames."
        write(*,*) "If you only want to calculate data for the first frame, input ""q"""
        read(*,"(a)") c2000tmp
        if (index(c2000tmp,'q')==0) then
			if (index(c2000tmp,'a')/=0) then
				open(10,file=filename,status="old")
				call xyzfile_nframe(10,nframe)
                close(10)
                write(*,"(a,i8 )") " Total number of frames:",nframe
				allocate(frameidx(nframe))
                forall(i=1:nframe) frameidx(i)=i
            else
				call str2arr(c2000tmp,nframe)
				allocate(frameidx(nframe))
				call str2arr(c2000tmp,nframe,frameidx)
            end if
            write(*,*) "Calculating, please wait..."
            open(11,file="MPP_SDP.txt",status="replace")
            open(12,file="ds.pqr",status="replace")
            !open(13,file="ds.pdb",status="replace") !If uncomment, ds.pdb will be exported, the B-factor corresponds to ds value
            open(10,file=filename,status="old")
			if (ifPBC>0) then !Write cell information
				call getcellabc(asize,bsize,csize,alpha,beta,gamma)
				write(12,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
			end if
			do iframe=1,maxval(frameidx)
				call readxyztrj(10)
                if (all(frameidx/=iframe)) cycle
				call ptsfitplane(atmsel,nsel,pleA,pleB,pleC,pleD,rmsfit) !Fit plane
                MPP=0
				iposmax=1
				inegmax=1
				distposmax=0
				distnegmax=0
				do idx=1,nsel
					iatm=atmsel(idx)
					call pointABCDdis(a(iatm)%x,a(iatm)%y,a(iatm)%z,pleA,pleB,pleC,pleD,dist,1)
					MPP=MPP+dist**2
					if (dist>0.and.dist>distposmax) then
						distposmax=dist
					else if (dist<0.and.dist<distnegmax) then
						distnegmax=dist
					end if
				end do
				MPP=dsqrt(MPP/nsel)
				SDP=distposmax-distnegmax
                write(11,"(i8,2f12.6)") iframe,MPP*b2a,SDP*b2a
                !Write pqr,pdb
				do i=1,ncenter
					devdis=0
					if (any(atmsel==i)) then
						call pointABCDdis(a(i)%x,a(i)%y,a(i)%z,pleA,pleB,pleC,pleD,dist,1)
						devdis=dist*b2a
					end if
					write(12,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
					"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,&
                    a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,devdis,vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
					!write(13,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
					!"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,1.0,devdis,adjustr(ind2name_up(a(i)%index))
				end do
                write(12,"(a)") "END"
                !write(13,"(a)") "END"
			end do
            write(*,"(a)") " Done!"
            write(*,"(a)") " Frame index, MPP and SDP have been exported to MPP_SDP.txt in current folder as the 1st, 2nd and 3rd columns, respectively"
            write(*,"(a)") " ds.pqr has been exported to current folder, the ""atomic charge"" column corresponds &
            &to the ds value of selected atoms in every frame (this value for unselected atoms is zero)"
            close(10)
            close(11)
            close(12)
            !close(13)
            return
		end if
    end if
end if

call ptsfitplane(atmsel,nsel,pleA,pleB,pleC,pleD,rmsfit) !Fit plane

write(*,*)
write(*,"(' Plane equation: A=',f10.5,'  B=',f10.5,'  C=',f10.5,'  D=',f10.5)") pleA,pleB,pleC,pleD
!write(*,"(' RMS of plane fitting:',f12.6)") rmsfit !This is different to MPP
write(*,*)

MPP=0
iposmax=1
inegmax=1
distposmax=0
distnegmax=0
do idx=1,nsel
	iatm=atmsel(idx)
	call pointABCDdis(a(iatm)%x,a(iatm)%y,a(iatm)%z,pleA,pleB,pleC,pleD,dist,1)
    write(*,"(' Deviation of atom',i5,'(',a,') to the plane:',f10.5,' Angstrom')") iatm,a(iatm)%name,dist*b2a
    MPP=MPP+dist**2
    if (dist>0.and.dist>distposmax) then
		distposmax=dist
        iposmax=iatm
    else if (dist<0.and.dist<distnegmax) then
		distnegmax=dist
        inegmax=iatm
    end if
end do
MPP=dsqrt(MPP/nsel)
SDP=distposmax-distnegmax
write(*,"(' Maximal positive deviation to the fitted plane is',i5,'(',a,'):',f10.5,' Angstrom')") iposmax,a(iposmax)%name,distposmax*b2a
write(*,"(' Maximal negative deviation to the fitted plane is',i5,'(',a,'):',f10.5,' Angstrom')") inegmax,a(inegmax)%name,distnegmax*b2a
write(*,"(/,' Molecular planarity parameter (MPP) is',f12.6,' Angstrom')") MPP*b2a
write(*,"(' Span of deviation from plane (SDP) is',f12.6,' Angstrom')") SDP*b2a
write(*,*)
write(*,"(a)") " Do you want to export .pqr of the system, in which ""charge"" atomic properties &
&of the selected atoms correspond to vertical distance to the fitted plane? (y/n)"
read(*,*) selectyn
if (selectyn=='Y'.or.selectyn=='y') then
	call path2filename(filename,outname)
    outname=trim(outname)//".pqr"
	write(*,*) "Input path of .pqr file, e.g. /sob/Akiyama/MIO.pqr"
    write(*,*) "If you press ENTER button directly, "//trim(outname)//" will be generated"
    read(*,"(a)") c200tmp
    if (c200tmp/=" ") outname=c200tmp
    open(10,file=trim(outname),status="replace")
	if (ifPBC>0) then !Write cell information
		call getcellabc(asize,bsize,csize,alpha,beta,gamma)
		write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
	end if
    do i=1,ncenter
		chgval=0
        if (any(atmsel==i)) then
			call pointABCDdis(a(i)%x,a(i)%y,a(i)%z,pleA,pleB,pleC,pleD,dist,1)
            chgval=dist*b2a
        end if
		write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
		"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,chgval,vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
	end do
    close(10)
    write(*,"(a)") " Done! In the newly generated .pqr file, ""charge"" property corresponds to the closest distance between &
    &atoms to fitted plane in Angstrom, positive and negative values correspond to whether the atom is above or below the plane. &
    &You can use e.g. VMD to color atoms corresponding this property"
end if
end subroutine




!!------------ Calculate diameter of cavity consisting of specific atoms in the system    
subroutine cavity_diameter
use defvar
use util
implicit real*8 (a-h,o-z)
integer,allocatable :: testatom(:)
character c2000tmp*2000
real*8 radmin,diamax,disp(3),coord(3),gvec(3),gvec_old(3),coordtmp(3),LSdisp(3),fract(3)

write(*,"(/,a)") " Input indices of a set of atoms to detect contact. e.g. 3-10,14,18-23"
write(*,*) "If pressing ENTER button directly, all atoms will be selected"
read(*,"(a)") c2000tmp
if (c2000tmp==" ") then
	ntestatom=ncenter
	allocate(testatom(ntestatom))
	forall(i=1:ncenter) testatom(i)=i
else
	call str2arr(c2000tmp,ntestatom)
	allocate(testatom(ntestatom))
	call str2arr(c2000tmp,ntestatom,testatom)
end if

write(*,*) "How to define initial sphere center?"
write(*,*) "1 Use geometric center of the atoms inputted in last step"
write(*,*) "2 Use a point by inputting its Cartesian coordinate"
if (ifPBC/=0) then
    write(*,*) "3 Use center of the cell"
	write(*,*) "4 Use a point by inputting its fractional coordinate"
end if
read(*,*) initcen

if (initcen==1) then
	coord(1)=sum(a(testatom(:))%x)/ntestatom
	coord(2)=sum(a(testatom(:))%y)/ntestatom
	coord(3)=sum(a(testatom(:))%z)/ntestatom
else if (initcen==2) then
	write(*,*) "Input X,Y,Z of the point in Angstrom, e.g. 2.5,0.79,0.2"
    read(*,*) coord(:)
    coord=coord/b2a
else if (initcen==3.or.initcen==4) then
	if (initcen==3) then
		fract=(/ 0.5D0,0.5D0,0.5D0 /)
    else
		write(*,*) "Input fractional coordinate of the point, e.g. 0,0.5,0.25"
		read(*,*) fract
    end if
	call fract2Cart(fract,coord)
end if

write(*,*) "If automatically adjusting the initial sphere center?"
write(*,*) "0 Do not adjust"
write(*,*) "1 Adjust in all directions"
write(*,*) "2 Adjust only in X"
write(*,*) "3 Adjust only in Y"
write(*,*) "4 Adjust only in Z"
write(*,*) "5 Adjust only in XY"
write(*,*) "6 Adjust only in XZ"
write(*,*) "7 Adjust only in YZ"
read(*,*) iadjust

if (iadjust==0) then
	write(*,"(' X/Y/Z of geometry center are',3f12.6,' Angstrom')") coord*b2a
    call getradmin(coord,testatom,ntestatom,radmin)
else
	!Iteratively update sphere center to find a position where radmin can be maximized
	!Barzilai¨CBorwein steepest ascent method cannot be used, because the object function is zigzag rather than smooth function
	write(*,"(' X/Y/Z of initial geometry center are',3f12.6,' Angstrom')") coord*b2a
	call getradmin(coord,testatom,ntestatom,radmin)
	write(*,"(' Initial sphere radius is',f12.6,' Angstrom')") radmin*b2a

	trustrad=0.3D0/b2a !Trust radius, 0.3 A
	dispconv=0.01D0/b2a !Displacement Threshold, 0.01 A
	stepscale=1D0
	nstepmax=100 !Maximum number of steps
	!Finite difference stepsize. Use a large value, not only because the object function is zigzag and thus derivative &
	!is invariant to stepsize in specific range, but also large stepsize can make position of sphere center can overcome slight barrier
	diff=0.3D0/b2a
    do istep=1,nstepmax
        write(*,*)
	    write(*,"(' Step',i5)") istep
    
        !Calculate gradient
        gvec=0
        do idir=1,3
            if ( (idir==1.and.(iadjust==1.or.iadjust==2.or.iadjust==5.or.iadjust==6)) &
            .or. (idir==2.and.(iadjust==1.or.iadjust==3.or.iadjust==5.or.iadjust==7)) &
            .or. (idir==3.and.(iadjust==1.or.iadjust==4.or.iadjust==6.or.iadjust==7)) ) then
                coordtmp=coord
                coordtmp(idir)=coord(idir)+diff
                call getradmin(coordtmp,testatom,ntestatom,tmpadd)
                coordtmp(idir)=coord(idir)-diff
                call getradmin(coordtmp,testatom,ntestatom,tmpmin)
                gvec(idir)=(tmpadd-tmpmin)/(2*diff)
            end if
        end do
        write(*,"(' Current coordinate:',3f12.6,' Angstrom')") coord*b2a
	    write(*,"(' Gradient:    ',3f12.6,'  Norm',f12.6)") gvec,dsqrt(sum(gvec**2))
        
        if (all(abs(gvec)<1D-6)) then
            disp=0
        else
            !Line search
            call getradmin(coord,testatom,ntestatom,radmin)
            LSdisp=0.5D0*gvec !Due to object function character, gvec must be a unit vector. So first step is 0.5 Bohr
            nmicro=0
            do nmicro=1,25
                coordtmp=coord+LSdisp
                call getradmin(coordtmp,testatom,ntestatom,radmintmp)
                !write(*,"(i4,2f12.6,3f16.10)") nmicro,radmin*b2a,radmintmp*b2a,LSdisp
                if (radmintmp>radmin) then
                    disp=coordtmp-coord
                    exit
                else
                    LSdisp=LSdisp*0.5D0
                end if
            end do
            if (nmicro==26) then !If finite difference stepsize is too large and thus inaccurate, then object functino cannot increase by following this direction
                diff=diff/2
				write(*,"(a,f10.6,' Bohr')") " Micro iteration in line search was not converged, reduce finite different stepsize by half to",diff
                cycle
            end if
            !write(*,"(' Number of micro steps in line search:',i6)") nmicro
        
            !Apply trust radius
            dispnorm=dsqrt(sum(disp**2))
            if (dispnorm>trustrad) then
                disp=disp*trustrad/dispnorm
            end if
    
	        coord=coord+stepscale*disp !Move coordinate
    
            if (ifPBC>0) call move_to_cell(coord,coord) !If moved to a position out of box, move it to central cell
        end if
    
        dispnorm=dsqrt(sum(disp**2))
	    write(*,"(' Displacement:',3f12.6,'  Norm',f12.6,' Angstrom')") disp*b2a,dispnorm*b2a
	    write(*,"(' Goal: displacement norm <',f12.8,' Angstrom')") dispconv*b2a
        if (dispnorm>dispconv) then
            write(*,"(' Not converged, new coordinate:',3f12.6,' Angstrom')") coord*b2a
        else
            write(*,"(/,' Converged after',i6,' iterations')") istep
            exit
        end if
        call getradmin(coord,testatom,ntestatom,radmin)
	    write(*,"(' Sphere radius at new coordinate:',f12.6,' Angstrom')") radmin*b2a
    end do
    write(*,"(/,' Final X/Y/Z of sphere center:',3f12.6,' Angstrom')") coord(:)*b2a
end if

diamax=2*radmin
write(*,"(' Radius is',f12.6,' Angstrom')") radmin*b2a
write(*,"(' Diameter is',f12.6,' Angstrom')") diamax*b2a
write(*,"(' Volume is',f14.6,' Angstrom^3')") 4D0/3D0*pi*(radmin*b2a)**3
write(*,*)
write(*,*) "Commands of drawing a sphere in VMD to show the cavity:"
write(*,*) "color Display Background white"
write(*,*) "draw material Transparent"
write(*,*) "draw color yellow"
write(*,"(a,3f9.3,a,f8.3,a)") " draw sphere {",coord(:)*b2a," } radius",radmin*b2a," resolution 100"
end subroutine

!------- Get minimum contact distance (w.r.t. atomic vdW sphere surface) between a point and all atoms in the list
subroutine getradmin(xyzA,testatom,ntestatom,radmin)
use defvar
implicit real*8 (a-h,o-z)
integer ntestatom,testatom(ntestatom)
real*8 radmin,xyzA(3),xyzB(3)
radmin=1D10
do idx=1,ntestatom
    iatm=testatom(idx)
    xyzB(1)=a(iatm)%x
    xyzB(2)=a(iatm)%y
    xyzB(3)=a(iatm)%z
    call nearest_dist(xyzA,xyzB,dist)
    contdist=dist-vdwr(a(iatm)%index)
    if (contdist<radmin) radmin=contdist
end do
end subroutine




!! ------------ Generate combined wavefunction by combining several fragment wavefunctions
!Any format of wavefunction file can be used as input, all orbitals including the virtual ones are stored to _all arrays first
subroutine gencombwfn
use defvar
implicit real*8 (a-h,o-z)
character selectyn*1,c200tmp*200
character,allocatable :: namearray(:)*200
type(atomtype),allocatable :: a_all(:)
type(primtype),allocatable :: b_all(:)
real*8,allocatable :: MOene_all(:),MOocc_all(:),CO_all(:,:),tmparr(:),tmparr2(:)
integer,allocatable :: shtype_all(:),shcon_all(:),shcen_all(:)
real*8,allocatable :: primshexp_all(:),primshcoeff_all(:),CObasa_all(:,:),CObasb_all(:,:)
integer,allocatable :: MOtype_all(:)
integer,allocatable :: iopshfrag(:),iflipspin(:)
real*8 naelec_all,nbelec_all

write(*,*) "0 Return"
write(*,*) "Generate which kind of file?"
write(*,*) "1: combine.wfn"
write(*,*) "2: combine.mwfn"
write(*,*) "3: combine.fch"
write(*,"(a)") " Note: For options 2 and 3, all inputted fragment wavefunction files must contain basis function information"
read(*,*) ifile
if (ifile==0) return

write(*,*) "How many fragments to combine? (Including the fragment 1 that has been loaded)"
write(*,"(a)") " Note: If the same basis functions and atoms were employed in every fragment calculation &
&(they only differ by the choice of ghost atoms and union of real atoms in each fragment corresponds to the whole system), input negative of number of fragments here. &
&In this case only single-determinant wavefunction is supported, and only occupied orbitals will be considered, while unoccupied orbitals will be ignored in the combined wavefunction file"
read(*,*) nfrag
isamebas=0
if (nfrag<0) isamebas=1
nfrag=abs(nfrag)
allocate(namearray(nfrag),iopshfrag(nfrag),iflipspin(nfrag))
do i=1,nfrag
	if (i==1) then
		write(*,"(' Filename of fragment 1: ',a)") trim(filename)
		namearray(1)=filename
	else if (i/=1) then
		do while(.true.)
			write(*,"(/,' Input wavefunction file of fragment',i4,', e.g. D:\combine\B.mwfn')") i
			if (ifile==1) then
				write(*,*) "(Any format of wavefunction file can be used, e.g. .wfn/wfx/mwfn/fch/molden)"
			else if (ifile==2.or.ifile==3) then
				write(*,*) "(The file must contain basis function information, e.g. .mwfn/fch/molden/gms)"
            end if
			read(*,"(a)") c200tmp
			inquire(file=c200tmp,exist=alive)
			if (alive) exit
			write(*,*) "File not found, input again"
		end do
		namearray(i)=c200tmp
	end if
end do

!Detect if need to treat this whole system as open-shell
iopsh=0
do i=1,nfrag
	call dealloall(0)
	call readinfile(namearray(i),1)
	if (wfntype==1.or.wfntype==4) then
		iopsh=1
		exit
	end if
end do

!Gain some basic information so that the arrays can be allocated
nprims_all=0
nshell_all=0
ncenter_all=0
if (ifile==2.or.ifile==3) then
	nbasis_all=0
    nprimshell_all=0
end if
nmo_all=0
nmoa_all=0
nmob_all=0
iopshfrag(:)=0 !Assume all fragments are closed-shell
naelec_all=0
nbelec_all=0
do i=1,nfrag
	call dealloall(0)
	write(*,"(/,' Loading ',a)") trim(namearray(i))
	call readinfile(namearray(i),1)
	ncenter_all=ncenter_all+ncenter
	nprims_all=nprims_all+nprims
    if (ifile==2.or.ifile==3) then
		if (nbasis==0) then
			write(*,*) "Error: This file does not contain basis function information!"
            write(*,*) "Press ENTER button to exit"
            read(*,*)
            stop
        end if
		nshell_all=nshell_all+nshell
        nprimshell_all=nprimshell_all+nprimshell
		nbasis_all=nbasis_all+nbasis
    end if
	if (iopsh==1) then !Open-shell treatment
		if (wfntype==1.or.wfntype==4) then !This is an open-shell fragment
			iopshfrag(i)=1
			nmoatmp=count(MOtype==1)
			nmobtmp=count(MOtype==2)
            write(*,"(' Number of alpha orbitals of this fragment:',i6)") nmoatmp
            write(*,"(' Number of beta orbitals of this fragment: ',i6)") nmobtmp
            write(*,"(' Number of alpha electrons of this fragment:',f10.6)") naelec
            write(*,"(' Number of beta electrons of this fragment: ',f10.6)") nbelec
            write(*,*)
			write(*,*) "If flipping electron spin for this fragment? (y/n)"
			read(*,*) selectyn
			if (selectyn=="y") then
				iflipspin(i)=1
				nmoa_all=nmoa_all+nmobtmp
				nmob_all=nmob_all+nmoatmp
				naelec_all=naelec_all+nbelec
				nbelec_all=nbelec_all+naelec
			else
				iflipspin(i)=0
				nmoa_all=nmoa_all+nmoatmp
				nmob_all=nmob_all+nmobtmp
				naelec_all=naelec_all+naelec
				nbelec_all=nbelec_all+nbelec
			end if
		else !This is a closed-shell fragment, separate it as equivalent alpha and beta parts
			nmoa_all=nmoa_all+nmo
			nmob_all=nmob_all+nmo
		end if
	else !Closed-shell treatment
		nmo_all=nmo_all+nmo
        naelec_all=naelec_all+naelec
        nbelec_all=nbelec_all+nbelec
	end if
end do
if (iopsh==1) nmo_all=nmoa_all+nmob_all
if (isamebas==1) then !This information is identical for all fragments, just take those of the last fragment
	ncenter_all=ncenter
	nprims_all=nprims
	nshell_all=nshell
    nprimshell_all=nprimshell
	nbasis_all=nbasis
    if (ifile==2.or.ifile==3) then
		if (iopsh==1) then !Open-shell treatment
			nmoa_all=nbasis
            nmob_all=nbasis
		else
			nmo_all=2*nbasis
		end if
    end if
end if

!Output brief information of combined wavefunction for check
write(*,*)
write(*,*) "Information of combined wavefunction:"
write(*,"(' Total number of atoms:',i6)") ncenter_all
write(*,"(' Total number of GTFs:',i6)") nprims_all
if (ifile==2.or.ifile==3) then
	write(*,"(' Total number of basis functions:',i6)") nbasis_all
	write(*,"(' Total number of basis function shells:',i6)") nshell_all
	write(*,"(' Total number of primitive shells:',i6)") nprimshell_all
end if
write(*,"(' Total number of orbitals:',i6)") nmo_all
if (iopsh==1) write(*,"(' Total number of alpha and beta orbitals:',2i6)") nmoa_all,nmob_all
write(*,*)

!Allocate array of combined wavefunction
allocate(a_all(ncenter_all),b_all(nprims_all),MOene_all(nmo_all),MOocc_all(nmo_all),MOtype_all(nmo_all),CO_all(nmo_all,nprims_all),tmparr(nprims_all))
CO_all=0
MOene_all=0
MOocc_all=0
MOtype_all=0
if (ifile==2.or.ifile==3) then
	allocate(shtype_all(nshell_all),shcen_all(nshell_all),shcon_all(nshell_all))
	allocate(primshexp_all(nprimshell_all),primshcoeff_all(nprimshell_all))
    allocate(CObasa_all(nbasis_all,nbasis_all))
    CObasa_all=0
    MOtype_all(1:nbasis)=1
    if (iopsh==1) then
		allocate(CObasb_all(nbasis_all,nbasis_all))
		CObasb_all=0
		MOtype_all(nbasis+1:nmo)=2
    end if
end if

!Read information from fragment wavefunction files and gradually construct combined wavefunction
icenter=1
iprim=1
ish=1
iprimsh=1
ibasis=1
imo=1
imoa=1
imob=nmoa_all+1
do i=1,nfrag
	write(*,"(a)") " Dealing with "//trim(namearray(i))//" ..."
	call dealloall(0)
	call readinfile(namearray(i),1)
    
    if (isamebas==0) then !Common case
		a_all(icenter:icenter+ncenter-1)=a(:)
		b_all(iprim:iprim+nprims-1)=b(:)
		b_all(iprim:iprim+nprims-1)%center=b_all(iprim:iprim+nprims-1)%center+(icenter-1)
		if (ifile==2.or.ifile==3) then
			shtype_all(ish:ish+nshell-1)=shtype(:)
			shcon_all(ish:ish+nshell-1)=shcon(:)
			shcen_all(ish:ish+nshell-1)=shcen(:)+(icenter-1)
			primshexp_all(iprimsh:iprimsh+nprimshell-1)=primshexp(:)
			primshcoeff_all(iprimsh:iprimsh+nprimshell-1)=primshcoeff(:)
		end if
		if (iopsh==0) then !Overall closed-shell situation
			MOene_all(imo:imo+nmo-1)=MOene
			MOocc_all(imo:imo+nmo-1)=MOocc
			MOtype_all(imo:imo+nmo-1)=0
			CO_all(imo:imo+nmo-1,iprim:iprim+nprims-1)=CO
			imo=imo+nmo
			if (ifile==2.or.ifile==3) CObasa_all(ibasis:ibasis+nbasis-1,ibasis:ibasis+nbasis-1)=CObasa
		else if (iopsh==1) then !Overall open-shell situation
			if (iopshfrag(i)==0) then !Closed-shell fragment
				!Alpha part
				MOene_all(imoa:imoa+nmo-1)=MOene
				MOocc_all(imoa:imoa+nmo-1)=MOocc/2D0
				MOtype_all(imoa:imoa+nmo-1)=0
				CO_all(imoa:imoa+nmo-1,iprim:iprim+nprims-1)=CO
				imoa=imoa+nmo
				if (ifile==2.or.ifile==3) CObasa_all(ibasis:ibasis+nbasis-1,ibasis:ibasis+nbasis-1)=CObasa
				!Beta part
				MOene_all(imob:imob+nmo-1)=MOene
				MOocc_all(imob:imob+nmo-1)=MOocc/2D0
				MOtype_all(imob:imob+nmo-1)=0
				CO_all(imob:imob+nmo-1,iprim:iprim+nprims-1)=CO
				imob=imob+nmo
				if (ifile==2.or.ifile==3) CObasb_all(ibasis:ibasis+nbasis-1,ibasis:ibasis+nbasis-1)=CObasa
			else !Open-shell fragment
				do isep=nmo,1,-1 !isep will be the last alpha MO
					if (MOtype(isep)==1) exit
				end do
				nmoatmp=count(MOtype==1)
				nmobtmp=count(MOtype==2)
				if (iflipspin(i)==0) then
					!Alpha part
					MOene_all(imoa:imoa+nmoatmp-1)=MOene(1:isep)
					MOocc_all(imoa:imoa+nmoatmp-1)=MOocc(1:isep)
					MOtype_all(imoa:imoa+nmoatmp-1)=1
					CO_all(imoa:imoa+nmoatmp-1,iprim:iprim+nprims-1)=CO(1:isep,:)
					imoa=imoa+nmoatmp
					if (ifile==2.or.ifile==3) CObasa_all(ibasis:ibasis+nbasis-1,ibasis:ibasis+nbasis-1)=CObasa
					!Beta part
					if (nmobtmp>0) then
						MOene_all(imob:imob+nmobtmp-1)=MOene(isep+1:nmo)
						MOocc_all(imob:imob+nmobtmp-1)=MOocc(isep+1:nmo)
						MOtype_all(imob:imob+nmobtmp-1)=2
						CO_all(imob:imob+nmobtmp-1,iprim:iprim+nprims-1)=CO(isep+1:nmo,:)
						imob=imob+nmobtmp
						if (ifile==2.or.ifile==3) CObasb_all(ibasis:ibasis+nbasis-1,ibasis:ibasis+nbasis-1)=CObasb
					end if
				else if (iflipspin(i)==1) then
					!Alpha part
					if (nmobtmp>0) then
						MOene_all(imoa:imoa+nmobtmp-1)=MOene(isep+1:nmo)
						MOocc_all(imoa:imoa+nmobtmp-1)=MOocc(isep+1:nmo)
						MOtype_all(imoa:imoa+nmobtmp-1)=1
						CO_all(imoa:imoa+nmobtmp-1,iprim:iprim+nprims-1)=CO(isep+1:nmo,:)
						imoa=imoa+nmobtmp
						if (ifile==2.or.ifile==3) CObasa_all(ibasis:ibasis+nbasis-1,ibasis:ibasis+nbasis-1)=CObasb
					end if
					!Beta part
					MOene_all(imob:imob+nmoatmp-1)=MOene(1:isep)
					MOocc_all(imob:imob+nmoatmp-1)=MOocc(1:isep)
					MOtype_all(imob:imob+nmoatmp-1)=2
					CO_all(imob:imob+nmoatmp-1,iprim:iprim+nprims-1)=CO(1:isep,:)
					imob=imob+nmoatmp
					if (ifile==2.or.ifile==3) CObasb_all(ibasis:ibasis+nbasis-1,ibasis:ibasis+nbasis-1)=CObasa
				end if
			end if
		end if
		icenter=icenter+ncenter
		iprim=iprim+nprims
		if (ifile==2.or.ifile==3) then
			ish=ish+nshell
			iprimsh=iprimsh+nprimshell
			ibasis=ibasis+nbasis
		end if
        
    else if (isamebas==1) then !Same basis functions were used for all fragments case. Only consider occupied orbitals
		!Add real atoms (nuclear charge is not zero) from fragment to a_all
		do iatm=1,ncenter
			if (a(iatm)%charge>0) a_all(iatm)=a(iatm)
        end do
		b_all(:)=b(:)
		if (ifile==2.or.ifile==3) then
			shtype_all(:)=shtype(:)
			shcon_all(:)=shcon(:)
			shcen_all(:)=shcen(:)
			primshexp_all(:)=primshexp(:)
			primshcoeff_all(:)=primshcoeff(:)
        end if
		if (iopsh==0) then !Overall closed-shell situation
			do itmp=1,nmo !Loop fragment orbitals to find occupied ones to include
				if (MOocc(itmp)>0) then
					MOene_all(imo)=MOene(itmp)
					MOocc_all(imo)=MOocc(itmp)
					MOtype_all(imo)=0
					CO_all(imo,:)=CO(itmp,:)
					if (ifile==2.or.ifile==3) CObasa_all(:,imo)=CObasa(:,itmp)
					imo=imo+1
                end if
            end do
		else if (iopsh==1) then !Overall open-shell situation
			if (iopshfrag(i)==0) then !Closed-shell fragment
				do itmp=1,nmo !Loop fragment orbitals to find occupied ones to include
					if (MOocc(itmp)>0) then
						MOene_all(imoa)=MOene(itmp)
						MOene_all(imob)=MOene(itmp)
						MOocc_all(imoa)=MOocc(itmp)/2
						MOocc_all(imob)=MOocc(itmp)/2
						MOtype_all(imoa)=1
						MOtype_all(imob)=2
						CO_all(imoa,:)=CO(itmp,:)
						CO_all(imob,:)=CO(itmp,:)
						if (ifile==2.or.ifile==3) then
							CObasa_all(:,imoa)=CObasa(:,itmp)
							CObasb_all(:,imob-nbasis)=CObasa(:,itmp)
                        end if
						imoa=imoa+1
						imob=imob+1
					end if
				end do
			else !Open-shell fragment
				do itmp=1,nmo !Loop fragment orbitals to find occupied ones to include
					if (MOocc(itmp)>0) then
                        if (iflipspin(i)==0) then
							if (MOtype(itmp)==1) then !Alpha
								MOene_all(imoa)=MOene(itmp)
								MOocc_all(imoa)=MOocc(itmp)
								MOtype_all(imoa)=1
								CO_all(imoa,:)=CO(itmp,:)
								if (ifile==2.or.ifile==3) CObasa_all(:,imoa)=CObasa(:,itmp)
								imoa=imoa+1
							else if (MOtype(itmp)==2) then !Beta
								MOene_all(imob)=MOene(itmp)
								MOocc_all(imob)=MOocc(itmp)
								MOtype_all(imob)=2
								CO_all(imob,:)=CO(itmp,:)
								if (ifile==2.or.ifile==3) CObasb_all(:,imob-nbasis)=CObasb(:,itmp-nbasis)
								imob=imob+1
							end if
                        else if (iflipspin(i)==1) then
							if (MOtype(itmp)==1) then !Alpha
								MOene_all(imob)=MOene(itmp)
								MOocc_all(imob)=MOocc(itmp)
								MOtype_all(imob)=2
								CO_all(imob,:)=CO(itmp,:)
								if (ifile==2.or.ifile==3) CObasb_all(:,imob-nbasis)=CObasa(:,itmp)
								imob=imob+1
							else if (MOtype(itmp)==2) then !Beta
								MOene_all(imoa)=MOene(itmp)
								MOocc_all(imoa)=MOocc(itmp)
								MOtype_all(imoa)=1
								CO_all(imoa,:)=CO(itmp,:)
								if (ifile==2.or.ifile==3) CObasa_all(:,imoa)=CObasb(:,itmp-nbasis)
								imoa=imoa+1
							end if
                        end if
					end if
				end do
			end if
        end if
    end if
end do

!Store the data to global arrays so that they can be exported
call dealloall(0)
allocate(a(ncenter_all),b(nprims_all),MOene(nmo_all),MOocc(nmo_all),MOtype(nmo_all),CO(nmo_all,nprims_all))
ncenter=ncenter_all
nprims=nprims_all
nmo=nmo_all
a=a_all
b=b_all
MOocc=MOocc_all
MOene=MOene_all
MOtype=MOtype_all
CO=CO_all
totenergy=0
virialratio=2
if (ifile==2.or.ifile==3) then
	nbasis=nbasis_all
	allocate(shtype(nshell_all),shcen(nshell_all),shcon(nshell_all),primshexp(nprimshell_all),primshcoeff(nprimshell_all))
    shtype=shtype_all
    shcen=shcen_all
    shcon=shcon_all
    nshell=nshell_all
    nprimshell=nprimshell_all
    primshexp=primshexp_all
    primshcoeff=primshcoeff_all
    naelec=naelec_all
    nbelec=nbelec_all
	allocate(CObasa(nbasis_all,nbasis_all))
	CObasa=CObasa_all
    if (iopsh==1) then
		allocate(CObasb(nbasis_all,nbasis_all))
		CObasb=CObasb_all
    end if
end if

!Determine wavefunction type
if (all(nint(MOocc)==MOocc)) then
	wfntype=0
	if (iopsh==0) then
		wfntype=0
		if (any(MOocc==1D0)) wfntype=2
	else if (iopsh==1) then
		wfntype=1
	end if
else
	wfntype=3
	if (iopsh==1) wfntype=4
end if

if (isamebas==0) then
	!Sort orbitals from occupancy from high to low. This guarantees that occupied orbitals occur prior to virtual ones
	!Sorting according to energies is meaningless, because energies of fragment orbitals evidently change when entering complex environment
	if (iopsh==0) ntime=1 !Closed-shell
	if (iopsh==1) ntime=2 !Open-shell
	allocate(tmparr2(nbasis))
	do itime=1,ntime
		if (iopsh==0) then
			ilow=1
			ihigh=nmo
		else if (iopsh==1) then !First time sort alpha orbitals, the second time sort beta orbitals
			if (itime==1) then !Alpha orbital range
				ilow=1
				ihigh=nmoa_all
			else if (itime==2) then !Beta orbital range
				ilow=nmoa_all+1
				ihigh=nmo_all
			end if
		end if
		do i=ilow,ihigh
			do j=i+1,ihigh
				if (MOocc(i)>=MOocc(j)) cycle
				temp=MOene(i)
				MOene(i)=MOene(j)
				MOene(j)=temp
				temp=MOocc(i)
				MOocc(i)=MOocc(j)
				MOocc(j)=temp
				itmp=MOtype(i)
				MOtype(i)=MOtype(j)
				MOtype(j)=itmp
				tmparr=CO(i,:)
				CO(i,:)=CO(j,:)
				CO(j,:)=tmparr
				if (ifile==2.or.ifile==3) then
					if (itime==1) then !Alpha or total
						tmparr2=CObasa(:,i)
						CObasa(:,i)=CObasa(:,j)
						CObasa(:,j)=tmparr2
					else if (itime==2) then !Beta
						tmparr2=CObasb(:,i-nbasis)
						CObasb(:,i-nbasis)=CObasb(:,j-nbasis)
						CObasb(:,j-nbasis)=tmparr2
					end if
				end if
			end do
		end do
	end do
end if

call updatenelec

if (ifile==1) then !Output .wfn file
	call outwfn("combine.wfn",1,1,10) !Note that the unoccupied MOs are automatically skipped
	write(*,*) "Combined wavefunction has been outputted to combine.wfn in current folder"
else if (ifile==2) then !Output .mwfn file
	call outmwfn("combine.mwfn",10,0)
	write(*,*) "Combined wavefunction has been outputted to combine.mwfn in current folder"
else if (ifile==3) then !Output .fch file
	call genP
	call outfch("combine.fch",10,0)
	write(*,*) "Combined wavefunction has been outputted to combine.fch in current folder"
end if
if (allocated(CObasa)) write(*,*) "The orbitals are sorted from occupation number from high to low"

!Recover to the first loaded system
call dealloall(0)
call readinfile(firstfilename,1)

end subroutine




!! ----------- Calculate Hellmann-Feynman forces at each nucleus
subroutine hellmann_feynman
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 HFforce_nuc(ncenter,3),HFforce_ele(ncenter,3),HFforce_tot(ncenter,3)
write(*,*) "Note: All units below are Hartree/Bohr"
write(*,*)
write(*,*) "Hellmann-Feynman forces contributed by electrons:"
write(*,*) "   Atom            X               Y               Z            Total"
diff=1D-5
do iatm=1,ncenter
	HFforce_ele(iatm,1)=(eleesp(a(iatm)%x+diff,a(iatm)%y,a(iatm)%z)-eleesp(a(iatm)%x-diff,a(iatm)%y,a(iatm)%z))/(2*diff)
	HFforce_ele(iatm,2)=(eleesp(a(iatm)%x,a(iatm)%y+diff,a(iatm)%z)-eleesp(a(iatm)%x,a(iatm)%y-diff,a(iatm)%z))/(2*diff)
	HFforce_ele(iatm,3)=(eleesp(a(iatm)%x,a(iatm)%y,a(iatm)%z+diff)-eleesp(a(iatm)%x,a(iatm)%y,a(iatm)%z-diff))/(2*diff)
	HFforce_ele(iatm,:)=-HFforce_ele(iatm,:)*a(iatm)%charge !Force is negative of gradient (1st derivative)
	write(*,"(i5,'(',a,')',4f16.8)") iatm,a(iatm)%name,HFforce_ele(iatm,:),dsqrt(sum(HFforce_ele(iatm,:)**2))
end do

write(*,*)
write(*,*) "Hellmann-Feynman forces contributed by nuclear charges:"
write(*,*) "   Atom            X               Y               Z            Total"
HFforce_nuc=0
do iatm=1,ncenter
	do jatm=1,ncenter
		if (jatm==iatm) cycle
		forcetmp=a(iatm)%charge*a(jatm)%charge/atomdist(iatm,jatm,0)**3
		HFforce_nuc(iatm,1)=HFforce_nuc(iatm,1)+forcetmp*(a(iatm)%x-a(jatm)%x)
		HFforce_nuc(iatm,2)=HFforce_nuc(iatm,2)+forcetmp*(a(iatm)%y-a(jatm)%y)
		HFforce_nuc(iatm,3)=HFforce_nuc(iatm,3)+forcetmp*(a(iatm)%z-a(jatm)%z)
	end do
	write(*,"(i5,'(',a,')',4f16.8)") iatm,a(iatm)%name,HFforce_nuc(iatm,:),dsqrt(sum(HFforce_nuc(iatm,:)**2))
end do

write(*,*)
HFforce_tot=HFforce_ele+HFforce_nuc
write(*,*) "Total Hellmann-Feynman forces:"
write(*,*) "   Atom            X               Y               Z            Total"
do iatm=1,ncenter
	write(*,"(i5,'(',a,')',4f16.8)") iatm,a(iatm)%name,HFforce_tot(iatm,:),dsqrt(sum(HFforce_tot(iatm,:)**2))
end do
end subroutine


!!------ Calculate attractive energy between an orbital and nuclei in a fragment
subroutine attene_orb_fragnuc
use functions
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000
real*8 intval,intvalold,funcval(radpot*sphpot),beckeweigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
write(*,*) "Input orbital index, e.g. 5"
read(*,*) iorb
write(*,"(a)") " Input atomic indices to define the fragment. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
read(*,"(a)") c2000tmp
if (allocated(frag1)) deallocate(frag1)
call str2arr(c2000tmp,nfragatm)
allocate(frag1(nfragatm))
call str2arr(c2000tmp,nfragatm,frag1)
write(*,"(i6,' atoms are selected',/)") nfragatm

write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)
intval=0
intvalold=0
do iatm=1,ncenter
	write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	!$OMP parallel do shared(funcval) private(i,jatmtmp,jatm,potnuc,rnowx,rnowy,rnowz) num_threads(nthreads)
	do i=1+iradcut*sphpot,radpot*sphpot
		rnowx=gridatm(i)%x
		rnowy=gridatm(i)%y
		rnowz=gridatm(i)%z
		potnuc=0
		do jatmtmp=1,nfragatm
			jatm=frag1(jatmtmp)
			potnuc=potnuc+a(jatm)%charge/dsqrt((rnowx-a(jatm)%x)**2+(rnowy-a(jatm)%y)**2+(rnowz-a(jatm)%z)**2)
		end do
 		funcval(i)=-potnuc*fmo(rnowx,rnowy,rnowz,iorb)**2
	end do
	!$OMP end parallel do
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	do i=1+iradcut*sphpot,radpot*sphpot
		intval=intval+funcval(i)*gridatmorg(i)%value*beckeweigrid(i)
	end do
	write(*,"(' Accumulated value:',f20.10,'  Current center:',f20.10)") intval,intval-intvalold
	intvalold=intval
end do
write(*,"(' Final result:',f14.6,' Hartree  ',f14.3,' kJ/mol')") intval,intval*au2kJ
end subroutine




!!----------- Output spherically averaged atomic radial density, can be used for generating promolecular density
subroutine sphatmraddens
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:),radpos(:),sphavgval(:)
truncrho=1D-8
rlow=0D0
rhigh=12
nsphpt=2030
nradpt=200 !Totally 200 radial points, but the number of point is truncated at truncrho
allocate(potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt),radpos(nradpt),sphavgval(nradpt))
call Lebedevgen(nsphpt,potx,poty,potz,potw)
ifinish=0
iprogstp=20
iprogcrit=iprogstp
write(*,*) "Calculating..."
!$OMP PARALLEL DO SHARED(sphavgval,radpos,ifinish,iprogcrit) PRIVATE(irad,radx,radr,isph,rnowx,rnowy,rnowz,tmpval) schedule(dynamic) NUM_THREADS(nthreads)
do irad=1,nradpt
	radx=cos(irad*pi/(nradpt+1))
	radr=(1+radx)/(1-radx) !Becke transform
	radpos(irad)=radr
	tmpval=0
	do isph=1,nsphpt
		rnowx=potx(isph)*radr
		rnowy=poty(isph)*radr
		rnowz=potz(isph)*radr
		tmpval=tmpval+fdens(rnowx,rnowy,rnowz)*potw(isph)
	end do
	sphavgval(irad)=tmpval !Spherically average density
    ifinish=ifinish+1
    if (ifinish==iprogcrit) then
		call showprog(ifinish,nradpt)
		iprogcrit=iprogcrit+iprogstp
	end if
end do
!$OMP END PARALLEL DO
open(10,file="sphavgval.txt",status="replace")
itmp=0
do irad=nradpt,1,-1
	if (sphavgval(irad)>truncrho) itmp=itmp+1
end do
write(10,"(a,i3,a)") "else if (iele==",a(1)%index,") then  !"
write(10,"('	npt=',i5)") itmp
itmp=0
do irad=nradpt,1,-1
	if (sphavgval(irad)>truncrho) then
		itmp=itmp+1
		write(10,"('	rhoarr(',i3,')=',f25.10,'D0')") itmp,sphavgval(irad)
	end if
end do
close(10)
write(*,*) "The result has been output to sphavgval.txt in current folder"
write(*,*) "The second column is radial distance (Bohr), the third column is value"
end subroutine




!!------- Load a .xyz trajectory, calculate distance between geometic center of two rings (fragments), and angle between the &
! normal vector of plane fitted for ring 1 and the line connecting the geometric centers of the two rings
subroutine ringring_geom
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer,allocatable :: tmparr(:)
real*8 geomcen1(3),geomcen2(3),vec1(3),vec2(3)

write(*,"(a)") " 1 Obtain distance between geometry centers of two fragments and angle between them"
write(*,*) "2 Angle between a fragment and Cartesian planes"
read(*,*) isel

write(*,*) "Input total number of frames, e.g. 4001"
read(*,*) nframetraj

allocate(frag1(ncenter),frag2(ncenter))

write(*,*) "Input index of the atoms in fragment 1, e.g. 2,3,7-10"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,nfrag1,frag1)

if (isel==1) then
    write(*,*) "Input index of the atoms in fragment 2, e.g. 2,3,7-10"
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,nfrag2,frag2)

    open(10,file=filename,status="old")
    open(11,file="distangle.txt",status="replace")
    !open(12,file="distangle2.txt",status="replace") !special

    call showprog(0,nframetraj)
    do iframe=1,nframetraj
	    call readxyztrj(10)
    
        geomcen1(1)=sum(a(frag1(1:nfrag1))%x)/nfrag1
        geomcen1(2)=sum(a(frag1(1:nfrag1))%y)/nfrag1
        geomcen1(3)=sum(a(frag1(1:nfrag1))%z)/nfrag1
        geomcen2(1)=sum(a(frag2(1:nfrag2))%x)/nfrag2
        geomcen2(2)=sum(a(frag2(1:nfrag2))%y)/nfrag2
        geomcen2(3)=sum(a(frag2(1:nfrag2))%z)/nfrag2
        cendist=dsqrt(sum((geomcen1-geomcen2)**2))*b2a
    
        call ptsfitplane(frag1,nfrag1,vec1(1),vec1(2),vec1(3),rnouse,rmsfit)
        vec2=geomcen2-geomcen1
        !write(*,"(' RMS error of the plane fitting:',f12.6,' Angstrom')") rmsfit*b2a
    
        !facnorm=sqrt(sum(vec1**2))
        !vec1=vec1/facnorm
        !write(*,"(' The vector 1 is',3f12.6)") vec1(:)
        !facnorm=sqrt(sum(vec2**2))
        !vec2=vec2/facnorm
        !write(*,"(' The vector 2 is',3f12.6,/)") vec2(:)
    
        angle=vecang(vec1(1),vec1(2),vec1(3),vec2(1),vec2(2),vec2(3))
        if (angle>90) angle=180-angle
    
        call ptsfitplane(frag2,nfrag2,vec2(1),vec2(2),vec2(3),rnouse,rmsfit)
        angle2=vecang(vec1(1),vec1(2),vec1(3),vec2(1),vec2(2),vec2(3))
        if (angle2>90) angle2=180-angle2
    
        write(11,"(i10,3f12.6)") iframe,cendist,angle,angle2
        !if (cendist<3) write(12,"(i10,2f12.6)") iframe,cendist,angle2 !special
    
        call showprog(iframe,nframetraj)
    end do
    write(*,*) "Results have been exported to distangle.txt in current folder"
    write(*,*) "Column 1: Frame index"
    write(*,*) "Column 2: Distance (Angstrom) between geometric centers of the two fragments"
    write(*,"(a)") " Column 3: Angle (degree) between normal vector of fragment 1 and linking line between geometric centers of the two fragments"
    write(*,"(a)") " Column 4: Angle (degree) between normal vector of fragment 1 and that of fragment 2"

else if (isel==2) then
    open(10,file=filename,status="old")
    open(11,file="angle.txt",status="replace")

    call showprog(0,nframetraj)
    do iframe=1,nframetraj
	    call readxyztrj(10)
        call ptsfitplane(frag1,nfrag1,vec1(1),vec1(2),vec1(3),rnouse,rmsfit)
        vec2=0;vec2(3)=1
        angleXY=vecang(vec1(1),vec1(2),vec1(3),vec2(1),vec2(2),vec2(3))
        if (angleXY>90) angleXY=180-angleXY
        vec2=0;vec2(1)=1
        angleYZ=vecang(vec1(1),vec1(2),vec1(3),vec2(1),vec2(2),vec2(3))
        if (angleYZ>90) angleYZ=180-angleYZ
        vec2=0;vec2(2)=1
        angleXZ=vecang(vec1(1),vec1(2),vec1(3),vec2(1),vec2(2),vec2(3))
        if (angleXZ>90) angleXZ=180-angleXZ
        write(11,"(i10,3f12.6)") iframe,angleXY,angleYZ,angleXZ
        call showprog(iframe,nframetraj)
    end do
    
    write(*,*) "Results have been exported to angle.txt in current folder"
    write(*,*) "Column 1: Frame index"
    write(*,*) "Column 2: Angle between fragment plane and XY plane"
    write(*,*) "Column 3: Angle between fragment plane and YZ plane"
    write(*,*) "Column 4: Angle between fragment plane and XZ plane"
end if

close(10)
close(11)
!close(12) !special
deallocate(frag1,frag2)
end subroutine



!!------- Calculate rotate angle velocity of C18 ring in OPP loop
subroutine ring_rotate
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer ringnatm
integer,allocatable :: ringatm(:)
type(atomtype) alast(ncenter)
real*8 cenloop(3),cen1(3),cenring(3),cenmove(3),vec1(3),vecc1(3,1),vecn(3),vecr(3),vecmove(3),vectan(3),mat(3,3)
integer loopnatm
integer,allocatable :: loopatm(:)

idetectout=1 !If skip frame if ring center deviates from loop center significantly (3 A)

write(*,*) "Input total number of frames, e.g. 4001"
!read(*,*) nframetraj
nframetraj=500001

write(*,*) "Input time interval between frames in ps, e.g. 2.0"
!read(*,*) timestep
timestep=0.2D0

write(*,*) "Input index of the atoms in ring, e.g. 2,3,7-10"
!read(*,"(a)") c2000tmp
c2000tmp="225-242"
!c2000tmp="243-260"
call str2arr(c2000tmp,ringnatm)
allocate(ringatm(ringnatm))
call str2arr(c2000tmp,ringnatm,ringatm)

if (idetectout==1) then
    write(*,*) "Input index of the atoms in the loop, e.g. 2,3,7-10"
    !read(*,"(a)") c2000tmp
    c2000tmp="16-18,20-22,43,49,62-63,82,97,112,127-128,145,149,154,179,186"
    !c2000tmp="23,28,33,38,52,57,64,69,100,105,159,164,169,174,197,202,207,212,217,222"
    call str2arr(c2000tmp,loopnatm)
    allocate(loopatm(loopnatm))
    call str2arr(c2000tmp,loopnatm,loopatm)
end if

write(*,*) "Input radius of the ring in Angstrom, e.g. 3.7"
!read(*,*) ringradius
ringradius=7.396D0/2 !wB97XD/6-311G* radius of C18
ringradius=ringradius/b2a

open(10,file=filename,status="old")
open(11,file="rotate.txt",status="replace")
!open(12,file="ring.xyz",status="replace")

nskip=0
iskiplast=0
call showprog(0,nframetraj)
do iframe=1,nframetraj
	call readxyztrj(10)
    
    if (iframe==1) then !Calculate geometry center of ring for first frame
        cen1(1)=sum(a(ringatm(:))%x)/ringnatm
        cen1(2)=sum(a(ringatm(:))%y)/ringnatm
        cen1(3)=sum(a(ringatm(:))%z)/ringnatm
    else !Calculate rotation with respect to last frame
        
        !Make normal vector coincide with first frame, so that rotation is removed
        !This is not used, because I found this treatment may cause artificial rotation and the magnitude may be quite large! (i.e. flip)
        !call ptsfitplane(ringatm,ringnatm,vec1(1),vec1(2),vec1(3),rnouse,rmsfit)
        !call rotmat_vec1_vec2(vec1,vecn1,mat)
        !do idx=1,ringnatm
        !    iatm=ringatm(idx)
        !    vecc1(1,1)=a(iatm)%x
        !    vecc1(2,1)=a(iatm)%y
        !    vecc1(3,1)=a(iatm)%z
        !    vecc1=matmul(mat,vecc1)
        !    a(iatm)%x=vecc1(1,1)
        !    a(iatm)%y=vecc1(2,1)
        !    a(iatm)%z=vecc1(3,1)
        !end do
        
        !Center of ring
        cenring(1)=sum(a(ringatm(:))%x)/ringnatm
        cenring(2)=sum(a(ringatm(:))%y)/ringnatm
        cenring(3)=sum(a(ringatm(:))%z)/ringnatm
    
        if (idetectout==1) then !Detect if ring is far from loop
            if (iskiplast==1) then !Last frame was skipped, so this frame cannot be calculated
                iskiplast=0
                alast=a
                nskip=nskip+1
                cycle
            end if
            cenloop(1)=sum(a(loopatm(:))%x)/loopnatm
            cenloop(2)=sum(a(loopatm(:))%y)/loopnatm
            cenloop(3)=sum(a(loopatm(:))%z)/loopnatm
            cendist=dsqrt(sum((cenring-cenloop)**2))
            if ( cendist*b2a > 3D0 ) then !Skip this frame
                iskiplast=1
                nskip=nskip+1
                cycle
            end if
        end if
        
        !Make ring center in current frame coincide with first frame, so that its translation is removed
        cenmove(:)=cenring(:)-cen1(:)
        do idx=1,ringnatm
            iatm=ringatm(idx)
            a(iatm)%x=a(iatm)%x-cenmove(1)
            a(iatm)%y=a(iatm)%y-cenmove(2)
            a(iatm)%z=a(iatm)%z-cenmove(3)
        end do
        
        call ptsfitplane(ringatm,ringnatm,vecn(1),vecn(2),vecn(3),rnouse,rmsfit) !Normal vector of the ring at first frame
        
        !Calculate tangential displacement for all atoms in the ring and then take average
        dtan_avg=0
        do idx=1,ringnatm
            iatm=ringatm(idx)
            !vecr is the vector pointing from ring center to reference atom
            vecr(1)=(a(iatm)%x+alast(iatm)%x)/2-cen1(1)
            vecr(2)=(a(iatm)%y+alast(iatm)%y)/2-cen1(2)
            vecr(3)=(a(iatm)%z+alast(iatm)%z)/2-cen1(3)
            !Construct tangential unit vector
            call crossprod(vecn,vecr,vectan)
            vectan=vectan/dsqrt(sum(vectan**2))
            !Movement vector of reference atom along tangential direction
            vecmove(1)=a(iatm)%x-alast(iatm)%x
            vecmove(2)=a(iatm)%y-alast(iatm)%y
            vecmove(3)=a(iatm)%z-alast(iatm)%z
            !Tangential displacement
            dtan=dot_product(vecmove,vectan)
            dtan_avg=dtan_avg+dtan
        end do
        dtan_avg=dtan_avg/ringnatm/timestep
        
        angvel=dtan_avg/ringradius
        write(11,"(i10,f12.2,4f10.4)") (iframe-1),(iframe-1)*timestep,dtan_avg*b2a,abs(dtan_avg)*b2a,angvel,abs(angvel)
        
    end if
    
    !Output processed ring geometry
    !write(12,*) ringnatm
    !write(12,"('Frame',i13)") iframe
    !do itmp=1,ringnatm
    !    iatm=ringatm(itmp)
	   ! write(12,"(a,3f16.8)") ind2name(a(iatm)%index),a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
    !end do
    
    alast=a
    if (mod(iframe,100)==0) call showprog(iframe,nframetraj)
end do
call showprog(nframetraj,nframetraj)

close(10)
close(11)
!close(12)

if (nskip>0) then
    write(*,"(/,' Number of skipped frames:',i10,' (',f8.3,' %)')") nskip,float(nskip)/nframetraj*100
else
    write(*,*) "No frame is skipped"
end if
write(*,*)
write(*,*) "Results have been exported to rotate.txt in current folder"
write(*,*) "Column 1: Frame"
write(*,*) "Column 2: Time (ps)"
write(*,*) "Column 3: Tangential velocity (signed, Angstrom*rad/ps)"
write(*,*) "Column 4: Tangential velocity (unsigned, Angstrom*rad/ps)"
write(*,*) "Column 5: Angle velocity (signed, rad/ps)"
write(*,*) "Column 6: Angle velocity (unsigned, rad/ps)"
!write(*,*)
!write(*,"(a)") " The ring coordinates with removal of overall rotation and translation with respect to the first frame has been exported to ring.xyz"
end subroutine




!!---------- Generate Fock/KS matrix based on orbital energies and coefficients
subroutine genFock
use defvar
implicit real*8 (a-h,o-z)
character outname*200
call MOene2Fmat(istatus)
if (istatus==0) then !Successfully generated
    write(*,*)
    write(*,*) "Input file path for exporting the generated matrix, e.g. C:\mol\Fock.txt"
    read(*,"(a)") outname
    open(10,file=outname,status="replace")
    write(10,"(4E18.8E3)") ((FmatA(i,j),j=1,i),i=1,nbasis)
    if (wfntype==1.or.wfntype==4) then
        write(10,"(4E18.8E3)") ((FmatB(i,j),j=1,i),i=1,nbasis)
    end if
    close(10)
    write(*,*) "Exporting finished!"
end if
end subroutine