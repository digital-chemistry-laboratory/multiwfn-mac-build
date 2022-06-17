!!!--------- Plot various kinds of spectra
!multiple.txt can records input file of multiple systems, but the system with the maximum number of transitions must be presented as the first entry
!If all(weight=1), then it is assumed that multiple.txt is not oriented for plotting weighted spectrum along with individual spectrum for each system, &
!but only for the latter, and in this case custom legend is allowed, which can be written in multiple.txt as the second column
!
!For single system, data is dataxall(1,1:numdata), nsystem=1, and numdataall(1)=numdata
!For n systems, nsystem=n, data of system i is dataxall(i,1:numdataall(i)), and in this case numdata is maximum length of numdataall
subroutine plotspectrum
use defvar
use dislin_d
use plot
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: weight(:) !Weight of various system for plotting mixed spectrum
real*8,allocatable :: dataxall(:,:),dataxall_org(:,:),strall(:,:),FWHMall(:,:) !Transition data loaded from multiple files. The first index corresponds to system index
integer,allocatable :: numdataall(:) !numdata is number of data, while for multiple system case, 
character(len=80),allocatable :: mollegend(:)
real*8,allocatable :: linexall(:,:),lineyall(:,:) !Array used to draw discrete lines for all systems. The first index corresponds to system index
real*8,allocatable :: curveyall(:,:) !The first index corresponds to system index. curvey in global array is used to record weighted curve
integer,allocatable :: tmparr(:),tmparr2(:)
real*8,allocatable :: tmparrr8(:),tmpmat(:,:)
real*8,allocatable :: indcurve(:,:) !Y value of curve of each individual band
real*8,allocatable :: indcontri(:) !Contribution of various transitions to a given X position
integer,allocatable :: indband2idx(:),idx2indband(:) !Used to map individual band index
character c80tmp*80,c200tmp*200,c200tmp2*200,strfmt*10,selectyn,graphformat_old*4,c2000tmp*2000
character clegend*2000 !Buffer for showing legends
integer :: icurveclr=1,ilineclr=5 !Default: Red for curve, black for discrete lines
integer :: thk_curve=3,thk_weighted=8,thk_Yeq0=2,thk_discrete=1,thk_axis=1,thk_grid=1,thk_PVS=3,thk_OPVS=3 !thickness
integer :: ishowlabelleft=1,ishowlabelright=1 !If showing labels on left and right Y-axes
integer :: ndecimalX=-1,ndecimalYleft=-1,ndecimalYright=-1 !Number of decimal places in axes, use auto by default
integer :: height_axis=36,ticksize=36,legtextsize=36,labtype_Yleft=1,ilegendpos=7
integer,parameter :: ncurrclr=15 !At most 15 individual colors. For more curves/lines, always use color of last one (8)
integer :: currclr(ncurrclr)=(/ 12,3,10,1,14,5,9,13,11,6,7,15,16,2,8 /) !Color index of current curve/line colors
integer :: iYeq0=1 !If drawing line corresponding to Y=0
real*8 :: degencrit=0.05D0,xlow=0,xhigh=1000,stepx=100
!Used for drawing spikes for indicating position of levels
integer,parameter :: maxspike=10
real*8,allocatable :: spikey(:) !Temporarily used for plotting spikes
integer,allocatable :: spikeidx(:,:) !Store level indices in each batch
integer :: spikenum(maxspike)=0 !The number of level indices in each batch
integer :: spikecolor(maxspike)=5 !Color of each spike batch, default to black
integer :: spikethick=3
!Used for showing extrema labels. The numbers are determined by numlocmax and numlocmin, which are initially zero and assigned when determining extrema
integer :: ishowextrema=0 !0=Do not show, 1=Show maxima, 2=Show minima, 3=Show both. Will be initialize to zero
integer minlabX(num1Dpoints),maxlabX(num1Dpoints) !Record the point index of the extrema in the curve
integer :: iextlabelrot=1,extlabeldecimal=1,extlabelsize=30,extlabelcontent=1,extlabelclr=3 !Blue
integer :: extmaxlabelshiftX=-16,extmaxlabelshiftY=20,extminlabelshiftX=20,extminlabelshiftY=-15 !Default label shifts, corresponding to rotated 90 degree case

!Partial vibrational spectrum (PVS) and overlap vibrational spectrum (OPVS) related arrays
!iVDOS, iPVSfragtype and iPVScomptype can combine freely, example:
!iVDOS=0,iPVScomptype=1,iPVSfragtype=1: PVS-NC(atom)
!iVDOS=1,iPVScomptype=1,iPVSfragtype=2: VDOS-NC(RIC)
!iVDOS=0,iPVScomptype=2,iPVSfragtype=1: PVS-I(atom)
!Nowever, combination between iPVScomptype=2 and iPVSfragtype=2 is not supported (calculating intensity for RICs is unsupported)
integer :: iPVSfragtype=1 !Type of PVS/PVDOS fragment, 1: Atoms  2: RIC
integer :: iPVScomptype=1 !Type of PVS composition, 1: Percentage contribution to normal coordinate (NC)  2: Percentage contribution to intensity (I. only available for IR now)
integer :: iVDOS=0 !1: Plot VDOS instead of PVS  0: PVS
integer :: OPVSidx1=0,OPVSidx2=0  !Indices of two PVS fragments for which OPVS will be plotted. 0,0 means do not plot OPVS
integer,parameter :: maxPVSfrag=10 !Maximum number of PVS fragments
integer :: PVScolor(maxPVSfrag)=(/ 1,3,14,9,13,15,16,11,7,8 /) !Color index of PVS/PVDOS curves
integer :: OPVScolor=12 !Color index of OPVS/OPVDOS curve
real*8,allocatable :: PVScomp(:,:),PVScompintra(:,:) !PVScomp(k,i) is composition of PVS fragment k in vibrational mode i, PVScompintra is intrafragment part of PVS-I
real*8,allocatable :: OPVScomp(:) !composition of OPVS in vibrational mode i
integer PVSnterm(maxPVSfrag) !Number of terms (atoms or RICs) in each PVS fragment
integer,allocatable :: PVSterm(:,:) !PVSterm(k,i) is index of term k in PVS fragment i
real*8,allocatable :: normmat_atm(:,:,:) !Temporarily record X/Y/Z component of all atoms in all normal vectors, normmat_atm(1:3,1:ncenter,1:numdata)
real*8 PVScurve(maxPVSfrag,num1Dpoints),OPVScurve(num1Dpoints) !Data point of PVS/PVDOS and OPVS/OPVDOS
real*8 PVScurveintra(maxPVSfrag,num1Dpoints) !Data point of intrafragment part of PVS-I
integer :: iPVSshow(maxPVSfrag)=1,iOPVSshow=1 !If showing PVS/PVDOS and OPVS/OPVDOS curves
character(len=80) PVSlegend(maxPVSfrag) !Legends of PVS/PVDOS curves
integer PVScart(maxPVSfrag) !Temporarily use, considered Cartesian components of atoms in fragment
integer :: nRIC=0 !Number of RICs
integer,allocatable :: RICatm(:,:) !RICatm(1:RICnatm(i),i) is the atom indices of the ith RIC loaded from Gaussian output
integer,allocatable :: RICnatm(:) !Number of atoms of various RIC
character,allocatable :: RICname(:)*6,RICdefname(:)*30 !Name and definition of RIC
real*8,allocatable :: strall_org(:) !Used for temporarily backup original strengths when plotting VDOS
integer :: ishowlocPVS=0
!Arrays used during calculating IR intensities
real*8,allocatable :: IRcomp(:,:),totIR(:) !Intensity of each mode of total IR
real*8,allocatable :: fragdipcomp(:,:,:,:) !x/y/z of fragment atoms contribution to x/y/z of dipole moment derivative w.r.t. q
integer,allocatable :: cartlist(:,:) !(3*ncenter,maxPVSfrag), if (j,i) is 1, then Cartesian coordinate j belongs to fragment i, otherwise belongs to others
real*8,allocatable :: dipder(:,:,:) !dipder(i,A,1/2/3) is derivative of x/y/z of dipole moment w.r.t i Cartesian coordinate of atom A
real*8,allocatable :: redmass(:) !Reduced masses of vibrational modes


if (ifiletype/=0) then
	if (ifiletype==1) then
		write(*,"(a)") " Error: As a Gaussian user, you must use output file (.out/.log) as input file for this function! &
		.fch/fchk does not contain information needed for plotting spectrum"
	else
		write(*,"(a)") " Error: The type of input file is wrong for plotting spectrum purpose! See Section 3.13.2 of manual for details"
	end if
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

!**** Spectrum types: 1=IR  2=Raman (or pre-resonance Raman)  3=UV-Vis  4=ECD  5=VCD  6=ROA
!
!Definition of units:  iunitx =0 cm^-1, =1 eV, =2 nm, =3 1000 cm^-1
!For vibrational spectrum, cm^-1 is always used, and used as internal unit; For electronic spectrum, eV (internal unit), nm, 1000 cm^-1 can be used 
!Note: For nm unit, we still store all data and FWHM in eV, and generate curve as usual in eV. Only at final stage, we scale the curve to get the one in nm
!If we choose 1000 cm^-1, we immediately convert all data and FWHM into 1000 cm^-1 before generating curve.
!When unit is changed, we reset lower and upper limit to auto rather than convert them to current unit to avoid problems.
!
!IR may use esu^2*cm^2 or km/mol as Y-axis unit, the data is always recorded in km/mol
!Unit conversion: 1 eV=8.0655447*1000 cm^-1    (1239.842/nm)eV    (1239.842/eV)nm     1esu^2*cm^2=2.5066km/mol

!! Initialize variables
gauweigh=0.5D0 !Gaussian weight used in Pseudo-Voigt broadening
iusersetY1=0 !User has not set the axes definition by himself
iusersetY2=0
iusersetX=0
orgy1=0;endy1=0;stepy1=0;orgy2=0;endy2=0;stepy2=0 !Temporarily used for scale left and right Y axes
isavepic=0
ishowline=1
ishowgrid=1
ishowlevel=0 !=1 means using spikes to indicate transition levels at bottom of the graph
iexportlevel=0
ishowweicurve=1 !0: Only show individual spectra, =1: Show both weighted spectrum and individual spectra, =2: Only show weighted spectrum
idegen=0
iweisyscurve=0
ishowextrema=0
iunitliney=1 !Only for IR
shiftx=0D0   !Shift value in X direction
iramantype=1 !1=Raman activities  2=Raman intensities
iROAtype=2   !1=Raman SCP(180)  2=ROA SCP(180)  3=Raman SCP(90)  4=ROA SCP(90)  5=Raman DCP(180)  6=ROA DCP(180)
graphformat_old=graphformat  !User may change format, backup the old
PVSnterm(:)=0

write(*,*) "Select type of the spectrum to plot"
write(*,*) "1:IR  2:Raman (or pre-resonance Raman)  3:UV-Vis  4:ECD  5:VCD  6:ROA  7:NMR"
read(*,*) ispectrum
if (ispectrum==7) then
    call NMRplot
    return
else if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then !IR, Raman, VCD, ROA
	ibroadfunc=1 !Use Lorentzian broadening
	iunitx=0 !cm^-1
else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
	ibroadfunc=2 !Use Gaussian broadening
	iunitx=2 !nm is default unit. But transition energies are loaded as eV
end if

if (ispectrum==2.or.ispectrum==4.or.ispectrum==5.or.ispectrum==6) then
!For ECD when eV is used, integrating the peak of a unit strength is assumed to be 1
!For Raman, VCD and ROA (cm^-1 is used), integrating the peak of a unit strength is assumed to be 1
	scalecurve=1D0
else if (ispectrum==1) then
!For IR when km/L is used, integrating the peak of a unit strength is 100
	scalecurve=100D0
else if (ispectrum==3) then
!1 unit oscillator strength can be broadened to 28700 area (X:eV Y:L/mol/cm)
!1 unit oscillator strength can be broadened to 1D0/4.32D-6 area (X:1000 cm^-1 Y:L/mol/cm)
!Their relationship: 1/(4.32D-9)/8065.5447=28700  1eV=8.0655447*1000*cm^-1 
!4.32D-9 can be found from Swizard manual (see also Review in C.C. vol.20 p168)
!The result is consistent with Gaussview
	if (iunitx==1.or.iunitx==2) scalecurve=28700 !nm,eV, both are recorded as eV internally
	if (iunitx==3) scalecurve=1D0/4.32D-6 !1000 cm^-1
end if


!! Load transition data from external files
nsystem=0
if (index(filename,"multiple.txt")/=0) then !Multiple file list with weights is recorded in multiple.txt
	open(11,file=filename,status="old") !Note that 10 is used by loadtransdata
	!Count total number of entries, find maximum number of data
    maxdata=0
	do while(.true.)
		read(11,"(a)",iostat=ierror) c200tmp
		if (ierror/=0.or.c200tmp==" ") exit
		read(c200tmp,*,iostat=ierror) c200tmp2
		nsystem=nsystem+1
        
		inquire(file=c200tmp2,exist=alive)
		if (.not.alive) then
			write(*,"(' Error: Cannot find ',a)") trim(c200tmp2)
			if (index(c200tmp,'/')/=0) then
				write(*,"(a)") " Reminder: Since the file path contains / symbol, you should add "" at the two ends of the path, so that the file can be properly loaded"
			end if
			write(*,*) "Press ENTER button to exit program"
			read(*,*)
			stop
        else !Obtain number of data
    			call loadtransdata(1,ispectrum,c200tmp2,numdata)
            if (numdata>maxdata) maxdata=numdata
		end if
	end do
	write(*,"(' There are',i4,' systems')") nsystem
	allocate(weight(nsystem),mollegend(nsystem))
    maxdata=maxdata*50 !Consider ORCA calculate both singlets and triplets, Gaussian may calculate anharmonic data, we need larger upperlimit
    allocate(dataxall(nsystem,maxdata),dataxall_org(nsystem,maxdata),strall(nsystem,maxdata),FWHMall(nsystem,maxdata),numdataall(nsystem))
	mollegend=" "
    !Actually load data
	rewind(11)
	do i=1,nsystem
		read(11,"(a)") c200tmp2
		read(c200tmp2,*,iostat=ierror) c200tmp,weight(i)
		if (ierror==0) then !The second field is weight
			write(mollegend(i),"(i3,' (',f5.1,'%)')") i,weight(i)*100
		else !The second field is legend rather than weight value
			ispc=index(c200tmp2," ")
            write(*,"(a)") trim(c200tmp2(:ispc-1))
            write(*,"(a)") trim(c200tmp2(ispc+1:))
			read(c200tmp2(:ispc-1),*) c200tmp
			read(c200tmp2(ispc+1:),"(a)") mollegend(i)
			weight(i)=1
			!If the first letter of the legend is $, it will be skipped
			if (mollegend(i)(1:1)=='$') mollegend(i)=mollegend(i)(2:)
		end if
		if (weight(i)==1) then
			write(*,"(' Loading ',a,'    Legend: ',a)") trim(c200tmp),trim(mollegend(i))
		else
			write(*,"(' Loading ',a,'    Weight:',f7.4)") trim(c200tmp),weight(i)
		end if
		call loadtransdata(0,ispectrum,c200tmp,numdata) !Data are loaded into datax,str,FWHM in global memory
		if (numdata>size(dataxall,2)) then !Very rare case, maxdata is still insufficiently large, user has to change order
			write(*,*)
			write(*,"(' The number of transitions in this file:',i6)") numdata
			write(*,"(' The size of data array:',i6)") size(dataxall,2)
			write(*,"(a)") " Error: You should put the system with maximum number of transitions to the first entry of multiple.txt"
			write(*,*) "Press ENTER button to exit program"
			read(*,*)
			stop
		end if
		dataxall(i,1:numdata)=datax
		strall(i,1:numdata)=str
		FWHMall(i,1:numdata)=FWHM
		numdataall(i)=numdata
	end do
	close(11)
	if (all(weight==1)) then !When all weights are unity, then no weighted spectrum will be plotted, but simply plotting all systems together
		ishowweicurve=0
		iweisyscurve=1
	end if
	!Below, numdata indicates actual maximum number of numdataall array
	numdata=maxval(numdataall)
else !Only one system
	nsystem=1
	allocate(weight(1))
	weight=1D0
	call loadtransdata(0,ispectrum,filename,numdata)
	allocate(dataxall(nsystem,numdata),dataxall_org(nsystem,numdata),strall(nsystem,numdata),FWHMall(nsystem,numdata),numdataall(nsystem))
	dataxall(1,:)=datax
	strall(1,:)=str
	FWHMall(1,:)=FWHM
	numdataall(1)=numdata
    allocate(strall_org(numdata))
    strall_org(:)=strall(1,:)
end if

dataxall_org=dataxall !Used to backup original data, because dataxall will be changed when setting scale factor

!! Allocate arrays properly
! curvey in global array is used to record weighted curve
! curveytmp is a temporary array used to record curve from each transition during generating curve of each system
if (allocated(curvex)) deallocate(curvex) !Global array
if (allocated(curvey)) deallocate(curvey) !Global array
if (allocated(curveytmp)) deallocate(curveytmp) !Global array
allocate(curvex(num1Dpoints),curvey(num1Dpoints),curveytmp(num1Dpoints),curveyall(nsystem,num1Dpoints))
allocate(linexall(nsystem,3*numdata),lineyall(nsystem,3*numdata))
 

!! Main interface
do while(.true.)
    write(*,"(/,a)") " Hint: You can input ""s"" to save current plotting settings to a file, or input ""l"" to load settings from a file"
	write(*,*)
    write(*,*) "               ================ Plot spectrum ==============="
	write(*,*) "-4 Set format of saving graphical file, current: "//graphformat
	write(*,*) "-3 Return to main menu"
	write(*,*) "-2 Export transition data to plain text file"
	write(*,*) "-1 Show transition data"
	write(*,*) "0 Plot spectrum"
	write(*,*) "1 Save graphical file of the spectrum in current folder"
	write(*,*) "2 Export X-Y data set of lines and curves to plain text file"
	if (iusersetX==0) write(*,*) "3 Set lower and upper limit of X-axis, current: Auto"
	if (iusersetX==1) write(*,"(a,f12.5,a,f12.5)") " 3 Set lower and upper limit of X-axis, current:",xlow," to",xhigh
	if (iusersetY1==0) write(*,*) "4 Set left Y-axis, current: Auto"
	if (iusersetY1==1) write(*,"(' 4 Set left Y-axis, current: low:',f12.3,' up:',f12.3,' step:',f11.3)") orgy1,endy1,stepy1
	if (iusersetY2==0) write(*,*) "5 Set right Y-axis, current: Auto"
	if (iusersetY2==1) write(*,"(' 5 Set right Y-axis, current: low:',f12.3,' up:',f12.3,' step:',f11.3)") orgy2,endy2,stepy2
	if (ibroadfunc==1) write(*,*) "6 Select broadening function, current: Lorentzian"
	if (ibroadfunc==2) write(*,*) "6 Select broadening function, current: Gaussian"
	if (ibroadfunc==3) write(*,*) "6 Select broadening function, current: Pseudo-Voigt"
	write(*,"(a,f20.5)") " 7 Set scale ratio for curve, current:",scalecurve
	do imol=1,nsystem
		if (any(FWHMall(imol,1:numdataall(imol))/=FWHMall(1,1))) exit
	end do
	if (imol==nsystem+1) then
		if (iunitx==0) then
			write(*,"(a,f10.5,' cm^-1')") " 8 Input full width at half maximum (FWHM), current:",FWHMall(1,1)
		else if (iunitx==1.or.iunitx==2) then
			!FWHM cannot be defined for nm, since it is not a linear unit, so what inputted is eV
			write(*,"(a,f10.5,' eV')") " 8 Input full width at half maximum (FWHM), current:",FWHMall(1,1)
		else if (iunitx==3) then
			write(*,"(a,f10.5,' 1000 cm^-1')") " 8 Input full width at half maximum (FWHM), current:",FWHMall(1,1)
		end if
	else
		write(*,*) "8 Set FWHM for all transitions, current: Loaded from input file"
	end if
	if (ishowline==1) write(*,*) "9 Toggle showing discrete lines, current: ON"
	if (ishowline==0) write(*,*) "9 Toggle showing discrete lines, current: OFF"
	if (ispectrum==1) then !IR allows using different strength units
		if (iunitliney==1) write(*,*) "10 Switch the unit of infrared intensity, current: km/mol"
		if (iunitliney==2) write(*,*) "10 Switch the unit of infrared intensity, current: esu^2*cm^2"
	else if (ispectrum==3.or.ispectrum==4) then !UV-Vis or ECD allows using different energy units
		if (iunitx==1) write(*,*) "10 Set the unit of excitation energy, current: eV"
		if (iunitx==2) write(*,*) "10 Set the unit of excitation energy, current: nm"
		if (iunitx==3) write(*,*) "10 Set the unit of excitation energy, current: 1000 cm^-1"
	end if
	if (ibroadfunc==3) write(*,"(a,f10.5)") " 11 Set Gaussian-weighting coefficient, current:",gauweigh
	write(*,"(a,f12.6)") " 12 Set shift value in X, current:",shiftx
	if (nsystem==1.or.any(weight/=1)) then
		write(*,*) "13 Set colors of curve and discrete lines"
	end if
	if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then
		write(*,*) "14 Set scale factor for vibrational frequencies"
	else if (ispectrum==3.or.ispectrum==4) then
		write(*,*) "14 Set scale factor for transition energies"
	end if
	if (nsystem==1) write(*,*) "15 Output contributions of individual transitions to the spectrum"
	if (.not.(nsystem>1.and.all(weight==1))) write(*,*) "16 Set status of showing labels of spectrum minima and maxima"
	write(*,*) "17 Other plotting settings"
	if (nsystem>1.and.any(weight/=1)) then
		if (iweisyscurve==1) write(*,*) "18 Toggle weighting spectrum of each system, current: ON"
		if (iweisyscurve==0) write(*,*) "18 Toggle weighting spectrum of each system, current: OFF"
	end if
	if (ispectrum==2) then !Raman
		if (iramantype==1) write(*,*) "19 Convert Raman activities to intensities"
		if (iramantype==2) write(*,*) "19 Convert Raman intensities to activities"
	else if (ispectrum==6) then !ROA
		if (iramantype==1) write(*,*) "19 Convert Raman or ROA activities to intensities"
		if (iramantype==2) write(*,*) "19 Convert Raman or ROA intensities to activities"
	end if
	if (ispectrum==1) then !IR
		write(*,*) "20 Modify IR strengths"
	else if (ispectrum==2) then !Raman
		if (iramantype==1) write(*,*) "20 Modify Raman activities"
		if (iramantype==2) write(*,*) "20 Modify Raman intensities"
	else if (ispectrum==3) then !UV-Vis
		write(*,*) "20 Modify oscillator strengths"
	else if (ispectrum==4.or.ispectrum==5) then !VCD,ECD
		write(*,*) "20 Modify rotatory strengths"
	else if (ispectrum==6) then !Frequency dependent Raman or ROA. The "intensity" is ambiguous, so use strength instead
		if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) write(*,*) "20 Modify Raman strengths"
		if (iROAtype==2.or.iROAtype==4.or.iROAtype==6) write(*,*) "20 Modify ROA strengths"
	end if
	if (any(weight/=1)) write(*,*) "21 Set status of showing weighted curve and curves of individual systems"
	write(*,*) "22 Set thickness of curves/lines/texts/axes/grid"
    if (nsystem==1) write(*,*) "23 Set status of showing spikes to indicate transition levels"
    if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then !IR, Raman, VCD, ROA
		write(*,*) "24 Set partial vibrational spectra (PVS) or vibrational DOS (VDOS)"
    end if
	read(*,"(a)") c80tmp
    
    if (index(c80tmp,'s')/=0) then
        write(*,"(a)") " Input file path for saving plotting settings, e.g. C:\Bang_Dream\RAS.dat"
        write(*,"(a)") " Note: If you press ENTER button directly, status will be saved to spectrum.dat in current folder"
        read(*,"(a)") c200tmp
        if (c200tmp==" ") c200tmp="spectrum.dat"
        open(10,file=c200tmp,status="replace")
        write(10,*) icurveclr
        write(10,*) ilineclr
        write(10,*) thk_curve
        write(10,*) thk_weighted
        write(10,*) thk_Yeq0
        write(10,*) thk_discrete
        write(10,*) thk_axis
        write(10,*) thk_grid
        write(10,*) spikethick
        write(10,*) spikecolor
        write(10,*) iextlabelrot
        write(10,*) extlabeldecimal
        write(10,*) extlabelsize
        write(10,*) extlabelcontent
        write(10,*) extlabelclr
        write(10,*) extmaxlabelshiftX
        write(10,*) extmaxlabelshiftY
        write(10,*) extminlabelshiftX
        write(10,*) extminlabelshiftY
        write(10,*) gauweigh
        write(10,*) iusersetY1
        write(10,*) iusersetY2
        write(10,*) iusersetX
        write(10,*) ishowline
        write(10,*) ishowgrid
        write(10,*) ishowlevel
        write(10,*) ishowweicurve
        write(10,*) idegen
        write(10,*) degencrit
        write(10,*) iweisyscurve
        write(10,*) ishowextrema
        write(10,*) iunitliney
        write(10,*) shiftx
        write(10,*) ibroadfunc
        write(10,*) iunitx
        write(10,*) scalecurve
        write(10,*) ispectrum
        write(10,*) FWHMall(1,1)
        write(10,*) xlow
        write(10,*) xhigh
        write(10,*) stepx
        write(10,*) orgy1
        write(10,*) endy1
        write(10,*) stepy1
        write(10,*) orgy2
        write(10,*) endy2
        write(10,*) stepy2
        write(10,"(a)") graphformat
        write(10,*) ishowlabelleft
        write(10,*) ishowlabelright
        write(10,"('height_axis',i5)") height_axis !Since 2020-Dec-23, settings are recorded with explicit labels
        write(10,"('ndecimalX',i5)") ndecimalX
        write(10,"('ndecimalYleft',i5)") ndecimalYleft
        write(10,"('ndecimalYright',i5)") ndecimalYright
        write(10,"('labtype_Yleft',i5)") labtype_Yleft
        write(10,"('ticksize',i5)") ticksize
        write(10,"('legtextsize',i5)") legtextsize
        write(10,"('ilegendpos',i5)") ilegendpos
        write(10,"('iYeq0',i5)") iYeq0
        write(10,"('thk_PVS',i5)") thk_PVS
        write(10,"('thk_OPVS',i5)") thk_OPVS
        write(10,"('iVDOS',i5)") iVDOS
        close(10)
        write(*,*) "Done!"
        cycle
    else if (index(c80tmp,'l')/=0)then
        write(*,"(a)") " Input file path to load plotting settings from it, e.g. C:\Bang_Dream\RAS.dat"
        write(*,"(a)") " Note: If you press ENTER button directly, status will be load from spectrum.dat in current folder"
        read(*,"(a)") c200tmp
        if (c200tmp==" ") c200tmp="spectrum.dat"
	    inquire(file=c200tmp,exist=alive)
	    if (.not.alive) then
	        write(*,*) "Error: Cannot find the file! Press ENTER button to return"
            read(*,*)
            cycle
        end if
        open(10,file=c200tmp,status="old")
        read(10,*) icurveclr
        read(10,*) ilineclr
        read(10,*) thk_curve
        read(10,*) thk_weighted
        read(10,*) thk_Yeq0
        read(10,*) thk_discrete
        read(10,*) thk_axis
        read(10,*) thk_grid
        read(10,*) spikethick
        read(10,*) spikecolor
        read(10,*) iextlabelrot
        read(10,*) extlabeldecimal
        read(10,*) extlabelsize
        read(10,*) extlabelcontent
        read(10,*) extlabelclr
        read(10,*) extmaxlabelshiftX
        read(10,*) extmaxlabelshiftY
        read(10,*) extminlabelshiftX
        read(10,*) extminlabelshiftY
        read(10,*) gauweigh
        read(10,*) iusersetY1
        read(10,*) iusersetY2
        read(10,*) iusersetX
        read(10,*) ishowline
        read(10,*) ishowgrid
        read(10,*) ishowlevel
        read(10,*) ishowweicurve
        read(10,*) idegen
        read(10,*) degencrit
        read(10,*) iweisyscurve
        read(10,*) ishowextrema
        read(10,*) iunitliney
        read(10,*) shiftx
        read(10,*) ibroadfunc
        read(10,*) iunitx
        read(10,*) scalecurve
        read(10,*) ispectrum
        read(10,*) FWHMall(1,1)
        read(10,*) xlow
        read(10,*) xhigh
        read(10,*) stepx
        read(10,*) orgy1
        read(10,*) endy1
        read(10,*) stepy1
        read(10,*) orgy2
        read(10,*) endy2
        read(10,*) stepy2
        read(10,"(a)") graphformat
        read(10,*) ishowlabelleft
        read(10,*) ishowlabelright
        call readoption_int(10,"height_axis",' ',height_axis)
        call readoption_int(10,"ndecimalX",' ',ndecimalX)
        call readoption_int(10,"ndecimalYleft",' ',ndecimalYleft)
        call readoption_int(10,"ndecimalYright",' ',ndecimalYright)
        call readoption_int(10,"ticksize",' ',ticksize)
        call readoption_int(10,"labtype_Yleft",' ',labtype_Yleft)
        call readoption_int(10,"legtextsize",' ',legtextsize)
        call readoption_int(10,"ilegendpos",' ',ilegendpos)
        call readoption_int(10,"iYeq0",' ',iYeq0)
        call readoption_int(10,"thk_PVS",' ',thk_PVS)
        call readoption_int(10,"thk_OPVS",' ',thk_OPVS)
        call readoption_int(10,"iVDOS",' ',iVDOS)
        close(10)
        write(*,*) "Loading finished!"
        cycle
    else if (c80tmp==" ") then
		cycle
    else
        read(c80tmp,*) isel
    end if
    
    if (isel==-4) then
        call setgraphformat
        
	else if (isel==-3) then
		istrtype=0
        graphformat=graphformat_old
		return
        
	else if (isel==-2) then !Export transition data
		if (nsystem==1) then
			open(10,file="transinfo.txt",status="replace")
			write(10,"(2i6)") numdata,2
			do i=1,numdata
				write(10,"(3f15.6)") dataxall(1,i),strall(1,i),FWHMall(1,i)
			end do
			close(10)
			write(*,"(a)") " The transition data have been exported to transinfo.txt in current directory, &
			this file can be directly used as input file of Multiwfn"
		else
			do imol=1,nsystem
				write(c200tmp,"(a,i3.3,a)") "transinfo",imol,".txt"
				open(10,file=c200tmp,status="replace")
				write(10,"(2i6)") numdataall(imol),2
				do i=1,numdataall(imol)
					write(10,"(3f15.6)") dataxall(imol,i),strall(imol,i),FWHMall(imol,i)
				end do
				close(10)
			end do
			write(*,"(a)") " The transition data have been exported to .txt with ""transinfo"" as prefix in current directory, &
			these files can be directly used as input file of Multiwfn"
		end if
        
	else if (isel==-1.or.isel==20) then !Show transition data and modify strengths
		do imol=1,nsystem
			if (nsystem>1) write(*,"(/,' Transition data of system',i5)") imol
			if (ispectrum==1) then !IR
				write(*,*) " Index  Freq.(cm^-1)  Intens.( km/mol   esu^2*cm^2)"
				do i=1,numdataall(imol)
					write(*,"(i6,1x,f12.5,7x,f12.5,f12.5)") i,dataxall(imol,i),strall(imol,i),strall(imol,i)/2.5066D0
				end do
			else if (ispectrum==2) then !Raman
				if (iramantype==1) write(*,*) " Index  Freq.(cm^-1)            Activities(A^4/amu)"
				if (iramantype==2) write(*,*) " Index  Freq.(cm^-1)                Intensity"
				do i=1,numdataall(imol)
					write(*,"(i6,3x,f12.5,7x,f18.5)") i,dataxall(imol,i),strall(imol,i)
				end do
			else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
				if (ispectrum==3) write(*,*) " Index  Excit.energy(eV       nm        1000 cm^-1)       Oscil.str."
				if (ispectrum==4) write(*,*) " Index  Excit.energy(eV       nm        1000 cm^-1)       Rotat.str."
				if (iunitx==3) dataxall=dataxall/8.0655447D0 !If unit is in 1000 cm^-1, temporarily convert to eV
				do i=1,numdataall(imol)
					write(*,"(i6,1x,4f15.5)") i,dataxall(imol,i),1239.842D0/dataxall(imol,i),8.0655447D0*dataxall(imol,i),strall(imol,i)
				end do
				if (iunitx==3) dataxall=dataxall*8.0655447D0 !Convert back from eV to 1000 cm^-1
			else if (ispectrum==5) then !VCD
				write(*,*) " Index   Freq.(cm^-1)        Rotat.str."
				do i=1,numdataall(imol)
					write(*,"(i6,3x,f12.5,7x,f12.5)") i,dataxall(imol,i),strall(imol,i)
				end do
			else if (ispectrum==6) then !ROA
				write(*,*) " Index   Freq.(cm^-1)         ROA str."
				do i=1,numdataall(imol)
					write(*,"(i6,3x,f12.5,7x,f12.5)") i,dataxall(imol,i),strall(imol,i)
				end do
			end if
		end do
		if (isel==20) then !Modify strengths
			imol=1
			if (nsystem>1) then
				write(*,*)
				write(*,*) "Input index of the system for which the strength(s) will be modified"
                write(*,*) "e.g. 2,4-6,10"
			    read(*,"(a)") c200tmp
			    call str2arr(c200tmp,ntmparr2)
			    allocate(tmparr2(ntmparr2))
			    call str2arr(c200tmp,ntmparr2,tmparr2)
            else
                ntmparr2=1
			    allocate(tmparr2(ntmparr2))
                tmparr2=1
			end if
            write(*,*)
			write(*,*) "Input index range of the transitions"
			write(*,*) "e.g. 1,3-6,22 means selecting modes 1,3,4,5,6,22"
            write(*,*) "Input ""q"" can return"
			read(*,"(a)") c200tmp
            if (c200tmp=="q") cycle
			call str2arr(c200tmp,ntmparr)
			allocate(tmparr(ntmparr))
			call str2arr(c200tmp,ntmparr,tmparr)
			write(*,*) "Input the strength to be set, e.g. 0.25"
			read(*,*) tmpval
            do i=1,ntmparr2
                imol=tmparr2(i)
			    do j=1,ntmparr
				    strall(imol,tmparr(j))=tmpval
			    end do
            end do
			deallocate(tmparr,tmparr2)
			write(*,*) "Done!"
		end if
	else if (isel==0) then !Draw curve
		isavepic=0
        
	else if (isel==1) then !Save curve picture
		isavepic=1
        
	else if (isel==3) then !Change X axis
		write(*,*) "Input lower limit, upper limit and step between ticks e.g. 200,1700,150"
		write(*,*) "Hint: If only input 0, the axis will be inverted"
		read(*,"(a)") c200tmp
		if (c200tmp=='0') then
			tmp=xlow
			xlow=xhigh
			xhigh=tmp
			stepx=-stepx
		else
			read(c200tmp,*,iostat=ierror) xlow,xhigh,stepx
            if (ierror/=0) then
                write(*,*) "Error: Unable to recognize the inputted content. Press ENTER button to continue"
                read(*,*)
                cycle
            end if
			if (xlow>xhigh.and.stepx>0) stepx=-stepx
		end if
		iusersetX=1 !User has modified it
        
	else if (isel==4) then !Change left Y axis
		if (orgy1==0.and.endy1==0.and.stepy1==0) then
			write(*,"(a)") " Note: To use this function, you should plot the graph at least once so that lower and upper limits of Y-axis could be initialized"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		orgy1old=orgy1
		endy1old=endy1
		stepy1old=stepy1
		write(*,*) "Input lower limit, upper limit and step between ticks e.g. 0,17000,2000"
		write(*,*) "Hint: If input 0, the axis will be inverted"
		read(*,"(a)") c200tmp
		if (c200tmp=='0') then
			tmp=orgy1
			orgy1=endy1
			endy1=tmp
			stepy1=-stepy1
		else
			read(c200tmp,*,iostat=ierror) orgy1,endy1,stepy1
            if (ierror/=0) then
                write(*,*) "Error: Unable to recognize the inputted content. Press ENTER button to continue"
                read(*,*)
                cycle
            end if
			if (orgy1>endy1.and.stepy1>0) stepy1=-stepy1
		end if
		iusersetY1=1
		write(*,"(a)") " Do you want to let program properly scale right Y axis so that its zero position exactly corresponds to left Y-axis? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			iusersetY2=1
			ratiotmp=(endy1old-orgy1old)/(endy2-orgy2)
			endy2=endy1/ratiotmp
			orgy2=orgy1/ratiotmp
			stepy2=stepy1/ratiotmp
		end if
        
	else if (isel==5) then !Change right Y axis
		if (orgy2==0.and.endy2==0.and.stepy2==0) then
			write(*,"(a)") " Note: To use this function, you should plot the graph at least once so that lower and upper limits of Y-axis could be initialized"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		orgy2old=orgy2
		endy2old=endy2
		stepy2old=stepy2
		write(*,*) "Input lower limit, upper limit and step between ticks e.g. -10,40,5"
		write(*,*) "Hint: If input 0, the axis will be inverted"
		read(*,"(a)") c200tmp
		if (c200tmp=='0') then
			tmp=orgy2
			orgy2=endy2
			endy2=tmp
			stepy2=-stepy2
		else
			read(c200tmp,*,iostat=ierror) orgy2,endy2,stepy2
            if (ierror/=0) then
                write(*,*) "Error: Unable to recognize the inputted content. Press ENTER button to continue"
                read(*,*)
                cycle
            end if
		end if
		iusersetY2=1
		write(*,"(a)") " Do you want to let program properly scale left Y axis so that its zero position exactly corresponds to right Y-axis? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			iusersetY1=1
			ratiotmp=(endy1-orgy1)/(endy2old-orgy2old)
			endy1=endy2*ratiotmp
			orgy1=orgy2*ratiotmp
			stepy1=stepy2*ratiotmp
		end if
        
	else if (isel==6) then !Set broadening function
        write(*,*) "Choose one of broadening functions:"
		write(*,*) "1 Lorentzian"
		write(*,*) "2 Gaussian"
		write(*,*) "3 Pseudo-Voigt"
		read(*,*) ibroadfunc
        
	else if (isel==7) then !Scale factor for curve
		write(*,*) "Input the scale factor, e.g. 0.8"
		read(*,*) scalecurve
        
	else if (isel==8) then !Set FWHM
		if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then !Always use cm^-1
			write(*,*) "Input the FWHM in cm^-1, e.g. 4"
		else if (ispectrum==3.or.ispectrum==4) then
			if (iunitx==2) write(*,"(a,/)") " NOTE: nm is not a linear unit of energy, so in principle one cannot define FWHM in nm. Nevertheless, in Multiwfn, when nm &
			is chosen as the unit, the curve will be generated in eV as X-axis first, and then convert to nm. Since current unit is nm, now you have to define the FWHM in eV."
			if (iunitx==1.or.iunitx==2) write(*,*) "Input the FWHM in eV, e.g. 0.5"
			if (iunitx==3) write(*,*) "Input the FWHM in 1000 cm^-1, e.g. 4"
		end if
		read(*,*) tmp
		FWHMall=tmp
        
	else if (isel==9) then !If show discrete lines
		if (ishowline==1) then
			ishowline=0
		else
			ishowline=1
		end if
        
	else if (isel==10) then !Change unit of X or Y axis. Not applied to VCD/ECD
		if (ispectrum==1) then !IR
			if (iunitliney==1) then
				iunitliney=2
			else
				iunitliney=1
			end if
		else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
			iusersetx=0 !Ensure auto X-axis is used, else will encounter problems 
			iold=iunitx
			write(*,*) "1: eV  2: nm  3: 1000 cm^-1"
			read(*,*) iunitx
			if (iunitx==1.or.iunitx==2) then
				if (iold==3) then !Convert data from 1000 cm^-1 to eV (For both eV and nm, transition energies are stored in eV)
					dataxall=dataxall/8.0655447D0
					FWHMall=FWHMall/8.0655447D0
					scalecurve=scalecurve/8.0655447D0 !1 unit oscillator strength can be broadened to 28700 area (X:eV Y:L/mol/cm)
				end if
			else if (iunitx==3) then
				if (iold==1.or.iold==2) then !Convert data from eV to 1000 cm^-1
					dataxall=dataxall*8.0655447D0
					FWHMall=FWHMall*8.0655447D0
					scalecurve=scalecurve*8.0655447D0
				end if
			end if
		end if
        
	else if (isel==11) then !Weight of Gaussian function
		write(*,*) "Input a value, e.g. 0.3"
		read(*,*) gauweigh
        
	else if (isel==12) then !Shift value in X
		write(*,*) "Input a value, e.g. 4.5"
		read(*,*) shiftx
        
	else if (isel==13) then !Set line and curve colors
		if (nsystem==1) write(*,*) "Use which color for curve?"
		if (nsystem>1) write(*,*) "Use which color for weighted curve?"
		call selcolor(icurveclr)
		write(*,*) "Use which color for discrete lines?"
		call selcolor(ilineclr)
        
	else if (isel==14) then !Set scale factor for transition energies or frequencies
		if (nsystem>1) then
			write(*,*) "Note: This operation will be applied to all systems loaded"
		end if
        if (any(dataxall/=dataxall_org)) then
			write(*,"(a)") " One or more data have been scaled previously, do you want to restore them to the original values? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y'.or.selectyn=='Y') then
				dataxall=dataxall_org !Recovered to initially loaded data
				write(*,*) "Original data have been restored"
            end if
        end if
        write(*,*)
		if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then !Vibrational spectra, mode selection is viable
			write(*,*) "Input the index range of the transitions you want to scaled"
			write(*,*) "e.g. 1,3-6,22 means selecting transitions 1,3,4,5,6,22"
            write(*,*) "If you want to select transitions according to frequency range, input ""f"" now"
			write(*,"(a)") " Note: Press ENTER button directly can select all modes. Input 0 can return"
            if (nsystem==1) write(*,"(a,i6,a)") " Note: There are",numdata," frequencies in total"
			read(*,"(a)") c200tmp
            if (index(c200tmp,'f')/=0) then !Select according to frequency, applied to all systems
				write(*,*) "Input lower limit of frequency in cm^-1, e.g. 1500"
                read(*,*) flow
				write(*,*) "Input upper limit of frequency in cm^-1, e.g. 4000"
                read(*,*) fup
				write(*,*) "Input scale factor, e.g. 0.97"
                read(*,*) tmpval
                nscl=0
                do isystem=1,nsystem
					do imode=1,numdataall(isystem)
						if (dataxall(isystem,imode)>=flow.and.dataxall(isystem,imode)<=fup) then
							dataxall(isystem,imode)=dataxall(isystem,imode)*tmpval
                            nscl=nscl+1
                        end if
					end do
                end do
                write(*,"(' Done!',i8,' frequencies have been scaled')") nsclall
			else if (c200tmp(1:1)=='0') then
                cycle
            else !Select according to index, or select 
				if (c200tmp==' '.or.index(c200tmp,"all")/=0) then !Selected all modes
					nmode=numdata
					allocate(tmparr(numdata))
					forall(itmp=1:numdata) tmparr(itmp)=itmp
				else
					call str2arr(c200tmp,nmode)
					allocate(tmparr(nmode))
					call str2arr(c200tmp,nmode,tmparr)
				end if
				write(*,"(i6,' frequencies are selected')") nmode
				write(*,"(/,a)") " Input scale factor, e.g. 0.97"
				write(*,"(a)") " Note: If pressing ENTER button directly, 0.9614 will be used, which is recommended for B3LYP/6-31G* level. &
				If inputting 1.0, frequencies will correspond to the ones originally loaded from input file"
				read(*,"(a)") c200tmp
				if (c200tmp==" ") then
					tmpval=0.9614D0
				else
					read(c200tmp,*) tmpval
				end if
				do idx=1,nmode
					dataxall(:,tmparr(idx))=dataxall(:,tmparr(idx))*tmpval
				end do
				deallocate(tmparr)
                write(*,*) "Done! Frequencies have been scaled"
            end if
		else !Electronic spectra, use universal scaling
			write(*,*) "Input the scale factor, e.g. 0.92"
			read(*,*) tmpval
			dataxall=dataxall*tmpval
			write(*,*) "Done! Transition energies have been scaled"
		end if
    
    else if (isel==16) then !Set showing method of extrema labels
        do while(.true.)
            write(*,*)
            write(*,*) "0 Return"
            if (ishowextrema==0) write(*,*) "1 Change displaying status of labels, current: Do not show"
            if (ishowextrema==1) write(*,*) "1 Change displaying status of labels, current: Show maxima"
            if (ishowextrema==2) write(*,*) "1 Change displaying status of labels, current: Show minima"
            if (ishowextrema==3) write(*,*) "1 Change displaying status of labels, current: Show minima and maxima"
            write(*,"(a,i3)") " 2 Set label size, current:",extlabelsize
            write(*,"(a,i3)") " 3 Set decimal digits, current:",extlabeldecimal
            if (iextlabelrot==0) write(*,*) "4 Switch label rotation, current: Do not rotate"
            if (iextlabelrot==1) write(*,*) "4 Switch label rotation, current: Rotated by 90 degree"
            write(*,*) "5 Set label color, current: "//trim(colorname(extlabelclr))
            if (extlabelcontent==1) write(*,*) "6 Switch label content, current: X-axis position"
            if (extlabelcontent==2) write(*,*) "6 Switch label content, current: Y-axis value"
            write(*,"(a,i4,',',i4)") " 7 Set shift of minimum labels in X and Y, current:",extminlabelshiftX,extminlabelshiftY
            write(*,"(a,i4,',',i4)") " 8 Set shift of maximum labels in X and Y, current:",extmaxlabelshiftX,extmaxlabelshiftY
            read(*,*) isellab
            if (isellab==0) then
                exit
            else if (isellab==1) then
                write(*,*) "0 Do not show extrema on the spectrum"
                write(*,*) "1 Show maxima on the spectrum"
                write(*,*) "2 Show minima on the spectrum"
                write(*,*) "3 Show both maxima and minima on the spectrum"
                read(*,*) ishowextrema
            else if (isellab==2) then
                isizeold=extlabelsize
                write(*,*) "Input label size, e.g. 25. The default is 30"
                read(*,*) extlabelsize
                if (iextlabelrot==0) extmaxlabelshiftY=extmaxlabelshiftY+(extlabelsize-isizeold)
            else if (isellab==3) then
                write(*,*) "Input decimal digits, e.g. 3. The default is 1"
                read(*,*) extlabeldecimal
            else if (isellab==4) then
                if (iextlabelrot==1) then
                    iextlabelrot=0
                    extmaxlabelshiftX=0
                    extmaxlabelshiftY=35
                    extmaxlabelshiftY=extmaxlabelshiftY+(extlabelsize-30) !30 is default Y shift
                    extminlabelshiftX=0
                    extminlabelshiftY=-15
                else if (iextlabelrot==0) then
                    iextlabelrot=1
                    extmaxlabelshiftX=-16
                    extmaxlabelshiftY=20
                    extminlabelshiftX=20
                    extminlabelshiftY=-15
                end if
            else if (isellab==5) then
                call selcolor(extlabelclr)
            else if (isellab==6) then
                if (extlabelcontent==1) then
                    extlabelcontent=2
                else if (extlabelcontent==2) then
                    extlabelcontent=1
                end if
            else if (isellab==7) then
                write(*,*) "Input shift of minimum labels in X and Y, respectively, e.g. 0,-20"
                read(*,*) extminlabelshiftX,extminlabelshiftY
            else if (isellab==8) then
                write(*,*) "Input shift of maximum labels in X and Y, respectively, e.g. 0,-20"
                read(*,*) extmaxlabelshiftX,extmaxlabelshiftY
            end if
        end do
        
	else if (isel==17) then !Other plotting settings
        do while(.true.)
            write(*,*)
            write(*,*) "0 Return"
	        if (ishowgrid==0) write(*,*) "1 Toggle showing dashed grid lines, current: No"
            if (ishowgrid==1) write(*,*) "1 Toggle showing dashed grid lines, current: Yes"
	        if (ishowlabelleft==0) write(*,*) "2 Toggle showing labels on left Y-axis, current: No"
            if (ishowlabelleft==1) write(*,*) "2 Toggle showing labels on left Y-axis, current: Yes"
	        if (ishowlabelright==0) write(*,*) "3 Toggle showing labels on right Y-axis, current: No"
            if (ishowlabelright==1) write(*,*) "3 Toggle showing labels on right Y-axis, current: Yes"
            if (ndecimalX==-1) then
                write(*,"(a)") " 4 Set number of decimal places in X-axis, current: Auto"
            else
                write(*,"(a,i2)") " 4 Set number of decimal places in X-axis, current:",ndecimalX
            end if
            if (ndecimalYleft==-1) then
                write(*,"(a)") " 5 Set number of decimal places in left Y-axis, current: Auto"
            else
                write(*,"(a,i2)") " 5 Set number of decimal places in left Y-axis, current:",ndecimalYleft
            end if
            if (ndecimalYright==-1) then
                write(*,"(a)") " 6 Set number of decimal places in right Y-axis, current: Auto"
            else
                write(*,"(a,i2)") " 6 Set number of decimal places in right Y-axis, current:",ndecimalYright
            end if
            write(*,"(a,i3)") " 7 Set text size of name of axes, current:",height_axis
            write(*,"(a,i3)") " 8 Set text size of ticks, current:",ticksize
            if (labtype_Yleft==1) write(*,"(a)") " 9 Set type of labels in left Y-axis, current: Float"
            if (labtype_Yleft==2) write(*,"(a)") " 9 Set type of labels in left Y-axis, current: Exponent"
            if (labtype_Yleft==3) write(*,"(a)") " 9 Set type of labels in left Y-axis, current: Scientific"
            if (nsystem>1.or.any(PVSnterm/=0)) then
                write(*,"(a,i3)") " 10 Set text size of legend, current:",legtextsize
                write(*,*) "11 Set position of legends"
            end if
            if (iYeq0==0) write(*,*) "12 Toggle drawing a line corresponding to Y-axis, current: No"
            if (iYeq0==1) write(*,*) "12 Toggle drawing a line corresponding to Y-axis, current: Yes"
            if (nsystem>1.and.ishowweicurve/=1) write(*,*) "13 Set color of line and curve of different systems"
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
		        if (ishowgrid==1) then
			        ishowgrid=0
		        else
			        ishowgrid=1
		        end if
            else if (isel2==2) then
		        if (ishowlabelleft==1) then
			        ishowlabelleft=0
		        else
			        ishowlabelleft=1
		        end if
            else if (isel2==3) then
		        if (ishowlabelright==1) then
			        ishowlabelright=0
		        else
			        ishowlabelright=1
		        end if
            else if (isel2==4) then
                write(*,*) "Input number of decimals to be kept, e.g. 2"
                read(*,*) ndecimalX
            else if (isel2==5) then
                write(*,*) "Input number of decimals to be kept, e.g. 2"
                read(*,*) ndecimalYleft
            else if (isel2==6) then
                write(*,*) "Input number of decimals to be kept, e.g. 2"
                read(*,*) ndecimalYright
            else if (isel2==7) then
                write(*,*) "Input size, e.g. 45"
                read(*,*) height_axis
            else if (isel2==8) then
                write(*,*) "Input size, e.g. 40"
                read(*,*) ticksize
            else if (isel2==9) then
                write(*,*) "1 Float"
                write(*,*) "2 Exponent"
                write(*,*) "3 Scientific"
                read(*,*) labtype_Yleft
            else if (isel2==10) then
                write(*,*) "Input size, e.g. 40"
                read(*,*) legtextsize
            else if (isel2==11) then
                write(*,*) "Choose position of legends"
                write(*,*) "5 Lower left corner"
                write(*,*) "6 Lower right corner"
                write(*,*) "7 Upper right corner"
                write(*,*) "8 Upper left corner"
                read(*,*) ilegendpos
            else if (isel2==12) then
                if (iYeq0==1) then
                    iYeq0=0
                else
                    iYeq0=1
                end if
            else if (isel2==13) then
				do while(.true.)
					write(*,"(a)") " The index and current color of various systems are listed below, &
                    now input index of a system to change its color, or input ""q"" to return"
					do isystem=1,nsystem
						if (isystem<=ncurrclr) then
							write(*,"(i4,2x,a)") isystem,colorname(currclr(isystem))
                        else
							write(*,"(i4,2x,a)") isystem,colorname(currclr(ncurrclr))
                        end if
                    end do
                    read(*,"(a)") c80tmp
                    if (index(c80tmp,'q')/=0) exit
                    read(c80tmp,*) isystem
                    write(*,*) "Set to which color?"
                    call selcolor(currclr(isystem))
                end do
            end if
        end do
        
	else if (isel==18) then !If weighting curve of each system
		if (iweisyscurve==1) then
			iweisyscurve=0
		else if (iweisyscurve==0) then
			iweisyscurve=1
		end if
        
	else if (isel==19) then !Convert between Raman activity and Raman intensity
		if (nsystem>1) write(*,*) "Note: This operation will be applied to all systems loaded"
		write(*,"(a)") " Reference of this Raman activity to intensity conversion: Chem. Asian J., 16, 56 (2021) DOI: 10.1002/asia.202001228"
		write(*,"(/,a)") " Input wavenumber of incident light, e.g. 532nm. If no unit is given, the unit will be default to cm^-1"
		read(*,"(a)") c200tmp
		ipos=index(c200tmp,'nm')
		if (ipos==0) then
			read(c200tmp,*) v0
		else
			read(c200tmp(:ipos-1),*) v0
			v0=1239.842D0/v0*8065.5447D0 !Convert nm to cm-1
		end if
		write(*,*) "Input temperature in K, e.g. 320"
		write(*,*) "Note: If pressing ENTER button directly, 298.15 K will be used"
        write(*,*) "      Inputting 0 will ignore the Boltzmann term"
		read(*,"(a)") c200tmp
		if (c200tmp==" ") then
			temper=298.15D0
		else
			read(c200tmp,*) temper
		end if
		Cfac=1D-12
		write(*,"(' Note: C factor of',1PE16.8,' is used')") Cfac
		do imol=1,nsystem
			do i=1,numdataall(imol)
				vi=dataxall(imol,i)
				if (temper==0) then
					Bfac=1D0
				else
					Bfac=1-exp( -vi*100*lightc*planckc/(boltzc*temper) )
				end if
				if (iramantype==1) strall(imol,i)= Cfac*(v0-vi)**4 /Bfac /vi * strall(imol,i) !Convert from activity to intensity
				if (iramantype==2) strall(imol,i)= strall(imol,i)*Bfac*vi /Cfac /(v0-vi)**4 !Convert from intensity to activity
			end do
		end do
		write(*,*) "Done!"
		if (iramantype==1) then
			iramantype=2
		else if (iramantype==2) then
			iramantype=1
		end if
        
	else if (isel==21) then !If weighting curve of each system
        write(*,*) "0: Only show spectra of individual systems"
        write(*,*) "1: Show both weighted spectrum and spectra of individual systems"
        write(*,*) "2: Only show weighted spectrum"
		read(*,*) ishowweicurve
        
	else if (isel==22) then
		do while(.true.)
			write(*,*)
			write(*,*) "0 Return"
			write(*,"(' 1 Set thickness of curves, current:',i3)") thk_curve
			write(*,"(' 2 Set thickness of weighted curves, current:',i3)") thk_weighted
			write(*,"(' 3 Set thickness of discrete lines, current:',i3)") thk_discrete
			if (iYeq0==1) write(*,"(' 4 Set thickness of Y=0 line, current:',i3)") thk_Yeq0
			write(*,"(' 5 Set thickness of axes, current:',i3)") thk_axis
			write(*,"(' 6 Set thickness of grid, current:',i3)") thk_grid
            if (any(PVSnterm/=0)) then
				write(*,"(' 7 Set thickness of PVS curves, current:',i3)") thk_PVS
				write(*,"(' 8 Set thickness of OPVS curves, current:',i3)") thk_OPVS
            end if
			read(*,*) isel2
			if (isel2==0) exit
			write(*,*) "Input the thickness, e.g. 3"
			if (isel2==1) read(*,*) thk_curve
			if (isel2==2) read(*,*) thk_weighted
			if (isel2==3) read(*,*) thk_discrete
			if (isel2==4) read(*,*) thk_Yeq0
			if (isel2==5) read(*,*) thk_axis
			if (isel2==6) read(*,*) thk_grid
			if (isel2==7) read(*,*) thk_PVS
			if (isel2==8) read(*,*) thk_OPVS
			write(*,*) "Done!"
		end do
        
    else if (isel==23) then !Set status of showing spikes
        if (.not.allocated(spikeidx)) allocate(spikeidx(maxspike,numdata))
        write(*,"(a)") " Note: You can use options 1~10 to define at most 10 sets of spikes"
        do while(.true.)
            write(*,*)
            if (iexportlevel==0) write(*,*) "-4 Toggle exporting spikes as spike.txt during plotting, current: No"
            if (iexportlevel==1) write(*,*) "-4 Toggle exporting spikes as spike.txt during plotting, current: Yes"
            if (idegen==1) write(*,"(a,f6.3)") " -3 Toggle considering degenerate, current: Yes, with threshold of ",degencrit
            if (idegen==0) write(*,"(a,f8.6)") " -3 Toggle considering degenerate, current: No"
            write(*,"(a,i3)") " -2 Set spike thick, current:",spikethick
            if (ispectrum==3.or.ispectrum==4) then
                if (ishowlevel==1) write(*,*) "-1 Toggle showing spikes to indicate transition energies, current: Yes"
                if (ishowlevel==0) write(*,*) "-1 Toggle showing spikes to indicate transition energies, current: No"
            else if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then
                if (ishowlevel==1) write(*,*) "-1 Toggle showing spikes to indicate vibrational levels, current: Yes"
                if (ishowlevel==0) write(*,*) "-1 Toggle showing spikes to indicate vibrational levels, current: No"
            end if
            write(*,*) " 0 Return"
            do ibatch=1,maxspike
                if (spikenum(ibatch)==0) then
                    write(*,"(i3,' Undefined spike set, choose to define')") ibatch
                else
                    write(*,"(i3,' Number of levels:',i6,',  Color: ',a)") ibatch,spikenum(ibatch),trim(colorname(spikecolor(ibatch)))
                end if
            end do
            read(*,*) isel2
            if (isel2==-4) then
                if (iexportlevel==0) then
                    iexportlevel=1
                    if (ishowlevel==0) write(*,"(a)") " Warning: This option takes effect only when status of option -1 has been charnged to ""Yes"""
                else
                    iexportlevel=0
                end if
            else if (isel2==-3) then
                if (idegen==1) then
                    idegen=0
                else
                    idegen=1
                    if (ispectrum==3.or.ispecturm==4) then
                        write(*,*) "Input threshold for determining degenerate in eV, e.g. 0.05"
                    else
                        write(*,*) "Input threshold for determining degenerate in cm^-1, e.g. 0.2"
                    end if
                    read(*,"(a)") c200tmp
                    if (c200tmp==" ") then
                        if (ispectrum==3.or.ispecturm==4) then
                            degencrit=0.05D0
                        else
                            degencrit=0.2D0
                        end if
                    else
                        read(c200tmp,*) degencrit
                    end if
                end if
            else if (isel2==-2) then
                write(*,*) "Input thick, e.g. 2"
                read(*,*) spikethick
            else if (isel2==-1) then
                if (ishowlevel==1) then
                    ishowlevel=0
                else
                    ishowlevel=1
                end if
            else if (isel2==0) then
                exit
            else
                nold=count(spikenum>0)
                write(*,*) "Input index of the levels to show, e.g. 2,3,7-10,44"
                if (spikenum(isel2)>0) then
                    write(*,*) "If input ""0"", this batch will be undefined"
                    write(*,*) "If press ENTER button directly, the definition will be kept unchanged"
                else
                    write(*,*) "If press ENTER button directly, all levels will be added"
                end if
                read(*,"(a)") c2000tmp
                if (spikenum(isel2)>0.and.c2000tmp(1:1)=="0") then
                    spikenum(isel2)=0
                    cycle
                end if
                if (spikenum(isel2)==0.and.c2000tmp==" ") then
                    spikenum(isel2)=numdata
                    forall (il=1:numdata) spikeidx(isel2,il)=il
                else if (spikenum(isel2)>0.and.c2000tmp==" ") then !Unchanged
                    continue
                else
                    call str2arr(c2000tmp,spikenum(isel2),spikeidx(isel2,:))
                    if (any(spikeidx(isel2,:)>numdata)) then
                        spikenum(isel2)=0
                        write(*,*) "Error: One or more indices exceeded valid range!"
                        cycle
                    end if
                end if
                write(*,*) "Use which color?"
		        call selcolor(spikecolor(isel2))
                if (nold==0.and.any(spikenum>0)) ishowlevel=1
            end if
        end do
        
    else if (isel==24) then !Set partial and overlap vibrational spectra (PVS and OPVS)
		if (nsystem>1) then
			write(*,*) "Error: This function is only available for single system!"
            write(*,*) "Press ENTER button to return"
			read(*,*)
            cycle
        end if
        if (ncenter==0) then !For xTB and CP2K cases, number of atoms has not been loaded from input file, guess it now so that arrays can be allocated
			if (mod((numdata+6),3)==0) then !Satisfy 3N-6
				ncenter=(numdata+6)/3
            else
				ncenter=(numdata+5)/3
            end if
			write(*,"(' Number of atoms is guessed to be',i8)") ncenter
        end if
        if (.not.allocated(PVSterm)) then
			allocate(PVScomp(maxPVSfrag,numdata),PVScompintra(maxPVSfrag,numdata),OPVScomp(numdata))
			allocate(PVSterm(10*ncenter,numdata)) !10*ncenter is large enough
            do ifrag=1,maxPVSfrag
                write(PVSlegend(ifrag),"('Fragment',i3)") ifrag
            end do
		end if
		icalccomp=0 !If 1, composition will be calculated when quit this option
		do while(.true.)
			write(*,*)
            call menutitle("Partial vibrational spectra (PVS) or vibrational DOS (VDOS)",10,1)
            if (any(PVSnterm/=0)) then
				if (icalccomp==1) then
					write(*,*) " q Generate composition data and return"
                else
					write(*,*) " q Return"
                end if
            else
				write(*,*) " q Return"
            end if
            if (iVDOS==0) then
				write(*,*) " d Set display status of PVS/OPVS curves"
				write(*,*) " c Set color of PVS/OPVS curves"
				write(*,*) " l Set legends of PVS curves"
				if (iPVScomptype==1.and.iPVSfragtype==1) write(*,*) " t Choose type of partial spectrum, current: PVS-NC(atom)"
				if (iPVScomptype==1.and.iPVSfragtype==2) write(*,*) " t Choose type of partial spectrum, current: PVS-NC(RIC)"
				if (iPVScomptype==2.and.iPVSfragtype==1) then
					if (ishowlocPVS==0) write(*,*) " t Choose type of partial spectrum, current: PVS-I(atom)"
					if (ishowlocPVS==1) write(*,"(a)") "  t Choose type of partial spectrum, current: PVS-I(atom) with intrafragment curves"
                end if
            else if (iVDOS==1) then
				write(*,*) " d Set display status of PVDOS/OPVDOS curves"
				write(*,*) " c Set color of PVDOS/OPVDOS curves"
				write(*,*) " l Set legends of PVDOS curves"
				if (iPVScomptype==1.and.iPVSfragtype==1) write(*,*) " t Choose type of partial spectrum, current: PVDOS-NC(atom)"
				if (iPVScomptype==1.and.iPVSfragtype==2) write(*,*) " t Choose type of partial spectrum, current: PVDOS-NC(RIC)"
            end if
            if (OPVSidx1==0.and.OPVSidx2==0) then
				if (iVDOS==0) write(*,*) " 0 Set OPVS, current: Undefined"
				if (iVDOS==1) write(*,*) " 0 Set OPVDOS, current: Undefined"
            else
				if (iVDOS==0) write(*,"(a,i3,a,i3)") "  0 Set OPVS, current: Between fragments",OPVSidx1," and",OPVSidx2
				if (iVDOS==1) write(*,"(a,i3,a,i3)") "  0 Set OPVDOS, current: Between fragments",OPVSidx1," and",OPVSidx2
            end if
            do iPVSfrag=1,maxPVSfrag
				if (PVSnterm(iPVSfrag)==0) then
					write(*,"(i3,' Define fragment  ',i3,', current: Undefined')") iPVSfrag,iPVSfrag
                else
					write(*,"(i3,' Redefine fragment',i3,', current:',i5,' terms, color: ',a)") iPVSfrag,iPVSfrag,PVSnterm(iPVSfrag),trim(colorname(PVScolor(iPVSfrag)))
                end if
            end do
            write(*,*) "Input e.g. -3 can unset fragment 3"
            read(*,"(a)") c80tmp
            
            if (index(c80tmp,'q')/=0) then !Load data, generate fragment compositions and then quit
				if (icalccomp==1.and.any(PVSnterm/=0)) then
                    open(10,file=filename,status="old")
                    call outputprog(10,iprog)
                    if (iprog==0) then !For xtb, load normal coordinates from g98.out
						call loclabel(10,"$vibrational spectrum",ixtb,maxline=100)
                        if (ixtb==1) then
							close(10)
							write(*,*) "Input path of the g98.out file outputted by xTB, e.g. D:\g98.out"
							do while(.true.)
								read(*,"(a)") c200tmp
								inquire(file=c200tmp,exist=alive)
								if (alive) exit
								write(*,*) "Cannot find the file, input again!"
							end do
                            open(10,file=c200tmp,status="old")
                        end if
                    end if
                    if (iPVSfragtype==1) then !Fragment based on atoms
                        if (iPVScomptype==1) then !Composition is computed according to contribution to normal coordinates
							write(*,*) "Loading normal vectors..."
							if (iprog==1.or.iprog==5.or.ixtb==1) then !Gaussian, CP2K, xTB (at least compatible with xTB 6.5)
								allocate(normmat_atm(3,ncenter,numdata))
								normmat_atm=0 !Note that vector of fixed atoms are not printed by Gaussian, so default to be zero
								ilackdata=numdata
								inow=1 !Normal mode index
								do while(.true.)
									if (ilackdata>3) then
										iread=3
									else
										iread=ilackdata
									end if
									if (iprog==1) then
										call loclabel(10,"Atom  AN",ifound,0)
									else if (iprog==5) then
										call loclabel(10,"ATOM  EL",ifound,0)
									else if (ixtb==1) then
										call loclabel(10,"Atom AN",ifound,0)
									end if
									read(10,*)
									do while(.true.)
										read(10,"(a)",iostat=ierror) c80tmp
										if (ierror/=0.or.c80tmp(1:6)==" ") exit
										read(c80tmp,"(i6)") iatm
										if (iread==1) read(c80tmp,*) itmp,c200tmp,normmat_atm(1:3,iatm,inow)
										if (iread==2) read(c80tmp,*) itmp,c200tmp,normmat_atm(1:3,iatm,inow),normmat_atm(1:3,iatm,inow+1)
										if (iread==3) read(c80tmp,*) itmp,c200tmp,normmat_atm(1:3,iatm,inow),normmat_atm(1:3,iatm,inow+1),normmat_atm(1:3,iatm,inow+2)
									end do
									if (ilackdata<=3) exit
									ilackdata=ilackdata-3
									inow=inow+3
								end do
							else if (iprog==2) then !ORCA. Note that ORCA outputs vectors of all normal modes including translations and rotations
								allocate(normmat_atm(3,ncenter,numdata))
								normmat_atm=0
								call loclabel(10,"The first frequency considered to be a vibration is")
								read(10,"(a)") c80tmp
								read(c80tmp(53:),*) itmp
								allocate(tmpmat(3*ncenter,numdata+itmp))
								call loclabel(10,"Thus, these vectors are normalized but",ierror)
								read(10,*)
								call readmatgau(10,tmpmat,inform="f11.6",inskipcol=11,inncol=6,iostat=ierror)
								if (ierror/=0) then
									write(*,*) "Error: Loading normal mode vectors was failed!"
									write(*,*) "Press ENTER button to exit program"
									read(*,*)
									stop
								end if
								do idata=1,numdata
									do iatm=1,ncenter
										do idir=1,3
											normmat_atm(idir,iatm,idata)=tmpmat(3*(iatm-1)+idir,itmp+idata)
										end do
									end do
								end do
								deallocate(tmpmat)
							end if
			  !              do imode=1,numdata !Check loaded normal mode
								!write(*,"(' Mode',i5)") imode
			  !                  do iatm=1,ncenter
								!	write(*,"(i5,3f12.6)") iatm,normmat_atm(1:3,iatm,imode)
			  !                  end do
			  !              end do
							!Generate compositions of fragments in various normal modes
							do idata=1,numdata !Cycle all vibrational modes
								facnorm=sum(normmat_atm(:,:,idata)**2)
								do ifrag=1,maxPVSfrag !Cycle all fragments
									if (PVSnterm(ifrag)/=0) then
										tmpval=0
										do idx=1,PVSnterm(ifrag)
											iatm=PVSterm(idx,ifrag)
											if (PVScart(ifrag)==0) then
												tmpval=tmpval+sum(normmat_atm(1:3,iatm,idata)**2)
											else if (PVScart(ifrag)==1) then
												tmpval=tmpval+normmat_atm(1,iatm,idata)**2
											else if (PVScart(ifrag)==2) then
												tmpval=tmpval+normmat_atm(2,iatm,idata)**2
											else if (PVScart(ifrag)==3) then
												tmpval=tmpval+normmat_atm(3,iatm,idata)**2
											else if (PVScart(ifrag)==4) then !XY
												tmpval=tmpval+normmat_atm(1,iatm,idata)**2+normmat_atm(2,iatm,idata)**2
											else if (PVScart(ifrag)==5) then !XZ
												tmpval=tmpval+normmat_atm(1,iatm,idata)**2+normmat_atm(3,iatm,idata)**2
											else if (PVScart(ifrag)==6) then !YZ
												tmpval=tmpval+normmat_atm(2,iatm,idata)**2+normmat_atm(3,iatm,idata)**2
											end if
										end do
										PVScomp(ifrag,idata)=tmpval/facnorm
									end if
								end do
                                OPVScomp(idata)=2*min(PVScomp(OPVSidx1,idata),PVScomp(OPVSidx2,idata))
							end do
							deallocate(normmat_atm)
                        else if (iPVScomptype==2) then !Composition is computed according to contribution to intensity
							!Load dipole derivatives, reduced masses and normal coordinates, which are needed for evaluating IR intensity
							if (.not.allocated(dipder)) then
								write(*,*) "Input path of fch/fchk file produced by frequency analysis task"
								write(*,*) "e.g. D:\test\ethanol.fch"
								do while(.true.)
									read(*,"(a)") c200tmp
									inquire(file=c200tmp,exist=alive)
									if (alive) exit
									write(*,*) "Cannot find the file, input again!"
								end do
								open(10,file=c200tmp,status="old")
                                write(*,*) "Loading electric dipole moment derivatives with respect to atoms..."
								!Load dipole derivatives
								allocate(dipder(3,ncenter,3))
								call loclabel(10,"Dipole Derivatives",ifound)
								if (ifound==0) then
									write(*,*) "Error: Unable to find ""Dipole Derivatives"" from the fch/fchk file!"
									write(*,*) "Press ENTER button to exit"
									read(*,*)
									stop
								end if
								allocate(tmparrr8(3*3*ncenter))
								read(10,*)
								read(10,*) tmparrr8(:)
								itmp=0
								do iatm=1,ncenter
									do j=1,3 !Atom Cartesian x/y/z component
										do idir=1,3 !Dipole moment direction
											itmp=itmp+1
											dipder(j,iatm,idir)=tmparrr8(itmp)
										end do
									end do
								end do
								deallocate(tmparrr8)
                                !Load reduced masses, they are recorded in part 2 in Vib-E2 field
                                write(*,*) "Loading reduced masses..."
								call loclabel(10,"Vib-E2",ifound,0)
								if (ifound==0) then
									write(*,*) "Error: Unable to find ""Vib-E2"" from the fch/fchk file!"
									write(*,*) "Press ENTER button to exit"
									read(*,*)
									stop
								end if
								allocate(tmparrr8(2*numdata),redmass(numdata))
								read(10,*)
								read(10,*) tmparrr8(:)
								redmass(:)=tmparrr8(numdata+1:2*numdata)
								deallocate(tmparrr8)
                                !Load normal coordinates
                                write(*,*) "Loading normal coordinates..."
								if (.not.allocated(normmat_atm)) allocate(normmat_atm(3,ncenter,numdata))
								call loclabel(10,"Vib-Modes",ifound,0)
								if (ifound==0) then
									write(*,*) "Error: Unable to find ""Vib-Modes"" from the fch/fchk file!"
									write(*,*) "Press ENTER button to exit"
									read(*,*)
									stop
								end if
								read(10,*)
								read(10,*) normmat_atm(:,:,:)
								close(10)
								!Convert normalized normal coordinates given by Gaussian to raw normal coordinates (in non-mass-weighted Cartesian form)
								do idata=1,numdata
                                    normmat_atm(:,:,idata)=normmat_atm(:,:,idata)/dsqrt(redmass(idata))
								end do
                            end if
                            write(*,*) "Calculating IR intensities..."
                            !See this about the coefficient 31.2231:
                            !https://mattermodeling.stackexchange.com/questions/5021/what-units-are-used-in-gaussian-16-for-dipole-derivatives-output
                            !Calculate total IR intensity
                            allocate(totIR(numdata),IRcomp(3,numdata),fragdipcomp(3,3,maxPVSfrag,numdata))
                            IRcomp=0
							do imode=1,numdata
								do idir=1,3
									do iatm=1,ncenter
										IRcomp(idir,imode)=IRcomp(idir,imode)+sum(dipder(:,iatm,idir)*normmat_atm(:,iatm,imode))
									end do
									IRcomp(idir,imode)=IRcomp(idir,imode)
								end do
								totIR(imode)=sum(IRcomp(:,imode)**2)
								!write(*,"(' Mode',i6,'   IR intensity:',f12.6,' km/mol')") imode,totIR(imode)*31.2231D0**2
							end do
                            !Construct a list recording which Cartesian coordinates belong to every fragment
                            !Order of the list: a1x,a1y,a1z,a2x,a2y,a2z...
                            allocate(cartlist(3*ncenter,maxPVSfrag))
                            cartlist=0
                            do ifrag=1,maxPVSfrag
                                if (PVSnterm(ifrag)==0) cycle
                                do iterm=1,PVSnterm(ifrag)
									itmp=3*(PVSterm(iterm,ifrag)-1)
									if (PVScart(ifrag)==0) then
										cartlist(itmp+1:itmp+3,ifrag)=1
									else if (PVScart(ifrag)==1) then !X
										cartlist(itmp+1:itmp+1,ifrag)=1
									else if (PVScart(ifrag)==2) then !Y
										cartlist(itmp+2:itmp+2,ifrag)=1
									else if (PVScart(ifrag)==3) then !Z
										cartlist(itmp+3:itmp+3,ifrag)=1
									else if (PVScart(ifrag)==4) then !XY
										cartlist(itmp+1:itmp+2,ifrag)=1
									else if (PVScart(ifrag)==5) then !XZ
										cartlist(itmp+1:itmp+1,ifrag)=1
										cartlist(itmp+3:itmp+3,ifrag)=1
									else if (PVScart(ifrag)==6) then !YZ
										cartlist(itmp+2:itmp+3,ifrag)=1
									end if
                                end do
                            end do
                            !Calculate fragment percentage contribution to IR intensity
                            fragdipcomp=0
                            do ifrag=1,maxPVSfrag !Cycle fragments
								do imode=1,numdata !Cycle vibration modes
                                    fragIR=0
                                    fragIRintra=0
									do idir=1,3 !Calculate fragment contribution to x,y,z component of total IR intensity
										do icart=1,3*ncenter !Loop Cartesian coordinates
											if (cartlist(icart,ifrag)==1) then
												iatm=ceiling(icart/3D0)
												icomp=icart-3*(iatm-1) !1/2/3=x,y,z of atom
                                                tmpval=dipder(icomp,iatm,idir)*normmat_atm(icomp,iatm,imode)
                                                fragdipcomp(icomp,idir,ifrag,imode)=fragdipcomp(icomp,idir,ifrag,imode)+tmpval
												do jcart=1,3*ncenter !Loop Cartesian coordinates
													jatm=ceiling(jcart/3D0)
													jcomp=jcart-3*(jatm-1)
													tmpval2=tmpval*dipder(jcomp,jatm,idir)*normmat_atm(jcomp,jatm,imode)
													fragIR=fragIR+tmpval2
                                                    if (cartlist(jcart,ifrag)==1) fragIRintra=fragIRintra+tmpval2
												end do
                                            end if
                                        end do
									end do
                                    if (totIR(imode)==0) then
										PVScomp(ifrag,imode)=0D0
                                    else
										PVScomp(ifrag,imode)=fragIR/totIR(imode)
										PVScompintra(ifrag,imode)=fragIRintra/totIR(imode)
                                    end if
                                end do
                            end do
                            !Calculate percentage coupling contribution between two fragments for plotting OPVS
                            if (OPVSidx1/=0.and.OPVSidx2/=0) then
								do imode=1,numdata !Cycle vibration modes
                                    coupIR=0
									do idir=1,3 !Calculate fragment contribution to x,y,z component of total IR intensity
										do icart=1,3*ncenter !Loop Cartesian coordinates
											if (cartlist(icart,OPVSidx1)==1) then
												iatm=ceiling(icart/3D0)
												icomp=icart-3*(iatm-1)
                                                tmpval=dipder(icomp,iatm,idir)*normmat_atm(icomp,iatm,imode)
												do jcart=1,3*ncenter !Loop Cartesian coordinates
													jatm=ceiling(jcart/3D0)
													jcomp=jcart-3*(jatm-1)
													tmpval2=tmpval*dipder(jcomp,jatm,idir)*normmat_atm(jcomp,jatm,imode)
													if (cartlist(jcart,OPVSidx2)==1) coupIR=coupIR+2*tmpval2
												end do
                                            end if
                                        end do
									end do
									if (totIR(imode)==0D0) then
										OPVScomp(imode)=0
                                    else
										OPVScomp(imode)=coupIR/totIR(imode)
                                    end if
                                end do
                            end if
                            write(*,"(a)") " If outputting information of dipole derivative w.r.t. atomic coordinate and atomic contribution to IRinten.txt in current folder? (y/n)"
                            read(*,*) selectyn
                            if (selectyn=='y') then
								open(10,file="IRinten.txt",status="replace")
								write(10,*) "Dipole moment derivative w.r.t. atomic Cartesian coordinates (a.u.)"
								write(10,"(a)") "             Dipole X       |      Dipole Y       |      Dipole Z"
								write(10,"(a)") "  Name  atom-X atom-Y atom-Z| atom-X atom-Y atom-Z| atom-X atom-Y atom-Z"
								do iatm=1,ncenter
									write(10,"(i5,a,3f7.3,'|',3f7.3,'|',3f7.3)") iatm,ind2name(a(iatm)%index),dipder(1:3,iatm,1),dipder(1:3,iatm,2),dipder(1:3,iatm,3)
								end do
                                do imode=1,numdata
									write(10,"(/,/,' *** Vibrational mode',i6,' (',f10.2,' cm^-1  Intens:',f10.3,' )')") imode,dataxall(1,imode),strall(1,imode)
									write(10,*) "X/Y/Z component of dipole derivative to normal coordinate q (e*amu^(-1/2)):"
									write(10,"(3f12.6)") IRcomp(1:3,imode)
									write(10,*) "X/Y/Z component of IR intensity:"
									write(10,"(3f12.6)") (IRcomp(1:3,imode)*31.2231D0)**2
									write(10,*)
									write(10,*) "Normal coordinate (amu^(-1/2)):"
									do iatm=1,ncenter
										write(10,"(i8,a,3f10.5)") iatm,ind2name(a(iatm)%index),normmat_atm(1:3,iatm,imode)
									end do
									write(10,*)
									write(10,*) "Atomic contribution to X/Y/Z component of dipole derivative to q:"
									write(10,"(a)") "             Dipole X       |      Dipole Y       |      Dipole Z"
									write(10,"(a)") "  Name  atom-X atom-Y atom-Z| atom-X atom-Y atom-Z| atom-X atom-Y atom-Z"
									do iatm=1,ncenter
										write(10,"(i8,a,3f7.3,'|',3f7.3,'|',3f7.3)") iatm,ind2name(a(iatm)%index),dipder(1:3,iatm,1)*normmat_atm(1:3,iatm,imode),dipder(1:3,iatm,2)*normmat_atm(1:3,iatm,imode),dipder(1:3,iatm,3)*normmat_atm(1:3,iatm,imode)
									end do
                                    tmp1=sum(dipder(1,:,1)*normmat_atm(1,:,imode))
                                    tmp2=sum(dipder(2,:,1)*normmat_atm(2,:,imode))
                                    tmp3=sum(dipder(3,:,1)*normmat_atm(3,:,imode))
                                    tmp4=sum(dipder(1,:,2)*normmat_atm(1,:,imode))
                                    tmp5=sum(dipder(2,:,2)*normmat_atm(2,:,imode))
                                    tmp6=sum(dipder(3,:,2)*normmat_atm(3,:,imode))
                                    tmp7=sum(dipder(1,:,3)*normmat_atm(1,:,imode))
                                    tmp8=sum(dipder(2,:,3)*normmat_atm(2,:,imode))
                                    tmp9=sum(dipder(3,:,3)*normmat_atm(3,:,imode))
                                    write(10,"(' Sum all: ',3f7.3,'|',3f7.3,'|',3f7.3)") tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
                                    do ifrag=1,maxPVSfrag
										if (PVSnterm(ifrag)>0) write(10,"(' Sum f',i2,': ',3f7.3,'|',3f7.3,'|',3f7.3)") ifrag,fragdipcomp(1:3,1,ifrag,imode),fragdipcomp(1:3,2,ifrag,imode),fragdipcomp(1:3,3,ifrag,imode)
                                    end do
                                    write(10,"('  Total:  ',5x,f10.5,6x,'|',5x,f10.5,6x,'|',5x,f10.5)") tmp1+tmp2+tmp3,tmp4+tmp5+tmp6,tmp7+tmp8+tmp9
                                end do
								close(10)
                                write(*,*) "IRinten.txt has been outputted to current folder"
                                write(*,*)
                            end if
                            deallocate(cartlist,IRcomp,totIR,fragdipcomp)
                        end if
                    else if (iPVSfragtype==2) then !Fragment based on redundant internal coordinates
						write(*,*) "Loading normal vectors..."
                        PVScomp(:,:)=0D0
						do idata=1,numdata !Cycle all vibrational modes
                            write(c80tmp,"('Normal Mode',i6)") idata
							call loclabel(10,trim(c80tmp),ifound,0)
                            if (ifound==0) then
								write(*,"(a)") " Error: Unable to find """//trim(c80tmp)//""" in from the Gaussian output file!"
                                write(*,"(a)") " Probably you forgot to add ""intmodes"" in ""freq"" keyword of your Gaussian input file"
                                write(*,*) "Press ENTER button to exit program"
                                read(*,*)
                                stop
                            end if
                            call skiplines(10,4)
                            do while(.true.)
								read(10,"(a)") c200tmp
                                if (index(c200tmp,'!')==0) exit
								do ifrag=1,maxPVSfrag !Cycle all fragments
									do iterm=1,PVSnterm(ifrag)
										if (index(c200tmp,RICname(PVSterm(iterm,ifrag)))/=0) then
											read(c200tmp(54:60),*) tmpval
                                            PVScomp(ifrag,idata)=PVScomp(ifrag,idata)+tmpval/100D0
                                        end if
                                    end do
                                end do
                            end do
                            OPVScomp(idata)=2*min(PVScomp(OPVSidx1,idata),PVScomp(OPVSidx2,idata))
						end do
                    end if
                    close(10)
		            do imode=1,numdata !Check loaded composition of various fragments in various normal modes
						write(*,"(' Vibrational mode',i6,' (',f10.2,' cm^-1  Intens:',f10.3,' )')") imode,dataxall(1,imode),strall(1,imode)
		                do ifrag=1,maxPVSfrag
							if (PVSnterm(ifrag)/=0) write(*,"('   Fragment',i3,'   Composition:',f10.4,' %')") ifrag,PVScomp(ifrag,imode)*100
		                end do
                        if (.not.(OPVSidx1==0.and.OPVSidx2==0)) write(*,"('   Overlap       Composition:',f10.4,' %')") OPVScomp(imode)*100
		            end do
					icurveclr=5 !Use black for total spectrum
					ishowline=0 !Do not show lines for clarity
                    thk_curve=5
                end if
				exit !Return to spectrum plotting interface
            else if (index(c80tmp,'0')/=0) then
				if (count(PVSnterm>0)<2) then
					write(*,*) "Error: You should define at least two fragments first!"
                    cycle
                end if
				write(*,"(a)") " Input indices of the two fragments for which overlap vibrational spectrum will be plotted, e.g. 1,3"
                write(*,*) "Note: If you input 0,0, then overlap vibrational spectrum will not be plotted"
                read(*,*) OPVSidx1,OPVSidx2
                if (OPVSidx1==0.and.OPVSidx2==0) cycle
                icalccomp=1
            else if (index(c80tmp,'c')/=0) then
				do while(.true.)
					write(*,*)
					write(*,*) " q Return"
					if (iVDOS==0) write(*,*) " 0 Set color for OPVS, current: "//trim(colorname(OPVScolor))
					if (iVDOS==1) write(*,*) " 0 Set color for OPVDOS, current: "//trim(colorname(OPVScolor))
					do ifrag=1,maxPVSfrag
						if (PVSnterm(ifrag)/=0) then
							write(*,"(i3,a,i3,a)") ifrag," Set color for fragment",ifrag,", current: "//trim(colorname(PVScolor(ifrag)))
						end if
					end do
					read(*,*) c80tmp
					if (index(c80tmp,'q')/=0) then
						exit
					else if (c80tmp(1:1)=='0') then
						write(*,*) "Use which color?"
						call selcolor(OPVScolor)
                    else
						read(c80tmp,*) ifrag
						write(*,*) "Use which color?"
						call selcolor(PVScolor(ifrag))
                    end if
                end do
            else if (index(c80tmp,'l')/=0) then
				do while(.true.)
					write(*,*)
					write(*,*) " q Return"
					do ifrag=1,maxPVSfrag
						if (PVSnterm(ifrag)/=0) then
							write(*,"(i3,a,i3,a)") ifrag," Set legend for fragment",ifrag,", current: "//trim(PVSlegend(ifrag))
						end if
					end do
					read(*,*) c80tmp
					if (index(c80tmp,'q')/=0) then
						exit
                    else
						read(c80tmp,*) ifrag
						write(*,*) "Input legend texts, e.g. love_nico"
						read(*,"(a)") PVSlegend(ifrag)
                    end if
                end do
            else if (index(c80tmp,'t')/=0) then
				iVDOS_old=iVDOS
                iPVSfragtype_old=iPVSfragtype
                write(*,*) " Meaning of partial spectrum:"
                write(*,*) "PVS:   Partial vibrational spectrum"
                write(*,*) "PVDOS: Partial vibrational density-of-states"
                write(*,*) "-NC: Composition corresponds to fragment contribution to normal coordinate"
                write(*,*) "-I:  Composition corresponds to fragment contribution to intensity"
                write(*,*) "(atom): Fragment is defined as a set of atoms"
                write(*,*) "(RIC):  Fragment is defined as a set of redundant internal coordinates"
                write(*,*)
                write(*,*) " Please choose the type of the partial spectrum you want to plot"
				write(*,*) "1 PVS-NC(atom)"
				write(*,*) "2 PVS-NC(RIC)"
				write(*,*) "3 PVS-I(atom)    -3 PVS-I(atom) with intrafragment curves"
				write(*,*) "4 PVDOS-NC(atom)"
				write(*,*) "5 PVDOS-NC(RIC)"
                read(*,*) itype
                if (abs(itype)==3) then
					if (ispectrum/=1) then
						write(*,*) "Error: PVS-I(atom) is currently only available for IR spectrum!"
						write(*,*) "Press ENTER button to return"
						read(*,*)
						cycle
					end if
                    open(10,file=filename,status="old")
                    call outputprog(10,iprog)
                    close(10)
                    if (iprog/=1) then
						write(*,*) "Error: PVS-I(atom) only supports vibration analysis of Gaussian!"
						write(*,*) "Press ENTER button to return"
						read(*,*)
						cycle
                    end if
                end if
                if (itype==1.or.itype==2.or.abs(itype)==3) then
					iVDOS=0
                else
					iVDOS=1
                end if
                if (abs(itype)==3) then
					iPVScomptype=2
					if (itype==3) then
						ishowlocPVS=0
					else if (itype==-3) then
						ishowlocPVS=1
					end if
                else
					iPVScomptype=1
                end if
                if (itype==1.or.abs(itype)==3.or.itype==4) then
					iPVSfragtype=1
                else
					iPVSfragtype=2
                end if
				if (iVDOS_old==0.and.iVDOS==1) then
                    strall(1,:)=1
                else if (iVDOS_old==1.and.iVDOS==0) then
                    strall(1,:)=strall_org(:)
                end if
                if (any(PVSnterm/=0)) then !Fragment type is changed, fragment must be cleaned
					if ((iPVSfragtype_old==1.and.iPVSfragtype==2).or.(iPVSfragtype_old==2.and.iPVSfragtype==1)) then
						PVSnterm(:)=0
                        write(*,*) "Note: Definition of previous fragments has been cleaned"
					end if
                end if
                if (iPVSfragtype_old==1.and.iPVSfragtype==2) then !Load RIC from Gaussian output file
					if (nRIC==0) then
						open(10,file=filename,status="old")
						call outputprog(10,iprog)
                        if (iprog/=1) then
							close(10)
							write(*,"(a)") " Error: This feature is only available for output file of Gaussian!"
							write(*,*) "Press ENTER button to return"
							read(*,*)
							cycle
                        end if
						write(*,*) "Loading redundant internal coordinates..."
						call loclabel(10,"! Name  Definition",ifound)
						if (ifound==0) then
							close(10)
							write(*,"(a)") " Error: Unable to find definition of redundant internal coordinates!"
							write(*,"(a)") " To guarantee that this information is always printed by Gaussian, the input file should contain ""freq=modr"" keyword"
							write(*,*) "Press ENTER button to return"
							read(*,*)
							cycle
						end if
						read(10,*);read(10,*)
						nRIC=0
						do while(.true.)
							read(10,"(a)") c200tmp
							if (index(c200tmp,"----")/=0) exit
							nRIC=nRIC+1
						end do
						allocate(RICatm(5,nRIC),RICnatm(nRIC),RICname(nRIC),RICdefname(nRIC)) !RIC may be linear bend, there involves five atoms
						do i=1,nRIC+1
							backspace(10)
						end do
                        do iRIC=1,nRIC
							read(10,"(a)") c200tmp
                            RICname(iRIC)=c200tmp(4:9)
                            RICdefname(iRIC)=c200tmp(10:30)
							ntmp=strcharnum(c200tmp,',')
							RICnatm(iRIC)=ntmp+1
							itmp=index(c200tmp,'(')
							jtmp=index(c200tmp,')')
							read(c200tmp(itmp+1:jtmp-1),*) RICatm(1:RICnatm(iRIC),iRIC)
                            !write(*,*) iRIC,RICname(iRIC),RICatm(1:RICnatm(iRIC),iRIC)
						end do
						close(10)
                        write(*,"(i6,' Redundant internal coordinates have been successfully loaded!')") nRIC
                    end if
                end if
            else if (index(c80tmp,'d')/=0) then
				do while(.true.)
					write(*,*) "q Return"
					if (iOPVSshow==0) then
						if (iVDOS==0) write(*,*) "0 Switch display status of OPVS, current: Do not show"
						if (iVDOS==1) write(*,*) "0 Switch display status of OPVDOS, current: Do not show"
					else if (iOPVSshow==1) then
						if (iVDOS==0) write(*,*) "0 Switch display status of OPVS, current: Show"
						if (iVDOS==1) write(*,*) "0 Switch display status of OPVDOS, current: Show"
                    end if
                    do ifrag=1,maxPVSfrag
						if (PVSnterm(ifrag)/=0) then
							if (iPVSshow(ifrag)==0) write(*,"(i2,a,i3,a)") ifrag," Switch display status of fragment",ifrag,", current: Do not show"
							if (iPVSshow(ifrag)==1) write(*,"(i2,a,i3,a)") ifrag," Switch display status of fragment",ifrag,", current: Show"
                        end if
                    end do
					read(*,*) c80tmp
					if (c80tmp=="q") then
						exit
					else if (c80tmp=="0") then
						if (iOPVSshow==0) then
							iOPVSshow=1
                        else
							iOPVSshow=0
                        end if
					else
						read(c80tmp,*) ifrag
						if (iPVSshow(ifrag)==0) then
							iPVSshow(ifrag)=1
                        else
							iPVSshow(ifrag)=0
                        end if
					end if
                end do
            else !Define PVS fragment
				icalccomp=1
				read(c80tmp,*) ifrag
				if (ifrag<0) then !Reset
					PVSnterm(abs(ifrag))=0
				else !Set fragment
					if (iPVSfragtype==1) then !PVS based on atoms
						write(*,"(a,i3,a)") " Input indices of the atoms in fragment",ifrag,", e.g. 2,3,7-10"
                        write(*,"(a)") " By default all Cartesian components are taken into account. You can also specify Cartesian component(s) as suffix, including X Y Z XY XZ YZ. &
                        For example, inputting ""2-4,9 XY"" will only define X and Y components of atoms 2,3,4,9 as the fragment"
						read(*,"(a)") c2000tmp
                        itmp=index(c2000tmp,' ')
						call str2arr(c2000tmp(1:itmp-1),PVSnterm(ifrag),PVSterm(:,ifrag))
                        if (c2000tmp(itmp+1:itmp+2)=="  ") then
							PVScart(ifrag)=0
                        else if (c2000tmp(itmp+1:itmp+2)=="X ") then
							PVScart(ifrag)=1
                        else if (c2000tmp(itmp+1:itmp+2)=="Y ") then
							PVScart(ifrag)=2
                        else if (c2000tmp(itmp+1:itmp+2)=="Z ") then
							PVScart(ifrag)=3
                        else if (c2000tmp(itmp+1:itmp+2)=="XY") then
							PVScart(ifrag)=4
                        else if (c2000tmp(itmp+1:itmp+2)=="XZ") then
							PVScart(ifrag)=5
                        else if (c2000tmp(itmp+1:itmp+2)=="YZ") then
							PVScart(ifrag)=6
                        end if
                    else if (iPVSfragtype==2) then !PVS based on RIC
						do while(.true.)
							write(*,*)
							if (PVSnterm(ifrag)/=0) then
								write(*,"(' There are',i6' redundant internal coordinates in this fragment')") PVSnterm(ifrag)
                            else
								write(*,*) "This fragment is empty currently"
                            end if
							write(*,*) "q Save and return"
                            write(*,*) "e Clean this fragment"
                            write(*,*) "p Print redundant internal coordinates in this fragment"
                            write(*,*) "c Add redundant internal coordinates to this fragment by conditions"
							write(*,"(a)") " Inputting 2/3/4 atoms can add corresponding bond/angle/dihedral coordinate"
                            write(*,*) "e.g. Inputting 5,8 will add bond (5,8) or (8,5)"
                            write(*,*) "     Inputting 5,8,1 will add angle (1,5,8) or (8,5,1)"
                            write(*,*) "     Inputting 5,8,1,2 will add angle (2,1,5,8) or (8,5,1,2)"
                            read(*,"(a)") c80tmp
                            if (index(c80tmp,'q')/=0) then
								exit
                            else if (index(c80tmp,'e')/=0) then
								PVSnterm(ifrag)=0
                                write(*,*) "This fragment has been cleaned"
                            else if (index(c80tmp,'p')/=0) then !Print RIC in current fragment
								do iterm=1,PVSnterm(ifrag)
									iRIC=PVSterm(iterm,ifrag)
									write(*,"(' Term',i6,':',2x,a,a)") iterm,RICname(iRIC),RICdefname(iRIC)
                                end do
                            else if (index(c80tmp,'c')/=0) then !Add RIC by conditions
								write(*,"(a)") " You will be asked to input conditions, the redundant internal coordinates (RICs) simultaneously satifying them will be added to this fragment"
                                write(*,*) "Which type of RICs may be added?"
                                write(*,*) "B: Bond    A: Angle   D: Dihedral   *: All"
                                read(*,*) c80tmp
                                if (c80tmp/='B'.and.c80tmp/='A'.and.c80tmp/='D'.and.c80tmp/='*') then
									write(*,*) "Error: Unable to recognize the type! It must be one of B A D *"
									cycle
                                end if
                                write(*,"(a)") " 1 RIC will be added if all its constituent atoms are in the inputted range"
                                write(*,"(a)") " 2 RIC will be added if any its constituent atom is in the inputted range"
                                read(*,*) icrit
                                write(*,"(a)") " Input the range of the atoms that make up the RICs, e.g. 3,9-12,18"
                                write(*,*) "If you press ENTER button directly, all atoms will be considered"
                                read(*,"(a)") c2000tmp
                                if (c2000tmp==" ") then
									allocate(tmparr(ncenter))
                                    forall(iatm=1:ncenter) tmparr(iatm)=iatm
                                else
									call str2arr(c2000tmp,ntmp)
									allocate(tmparr(ntmp))
									call str2arr(c2000tmp,ntmp,tmparr)
                                end if
                                numadd=0
                                do iRIC=1,nRIC
									iadd=0
									if ((c80tmp=='B'.or.c80tmp=='*').and.RICnatm(iRIC)==2) then
										if (icrit==1) then
											if (any(tmparr==RICatm(1,iRIC)).and.any(tmparr==RICatm(2,iRIC))) iadd=1
										else if (icrit==2) then
											if (any(tmparr==RICatm(1,iRIC)).or.any(tmparr==RICatm(2,iRIC))) iadd=1
                                        end if
                                    end if
									if ((c80tmp=='A'.or.c80tmp=='*').and.RICnatm(iRIC)==3) then
										if (icrit==1) then
											if (any(tmparr==RICatm(1,iRIC)).and.any(tmparr==RICatm(2,iRIC)).and.any(tmparr==RICatm(3,iRIC))) iadd=1
										else if (icrit==2) then
											if (any(tmparr==RICatm(1,iRIC)).or.any(tmparr==RICatm(2,iRIC)).or.any(tmparr==RICatm(3,iRIC))) iadd=1
                                        end if
                                    end if
									if ((c80tmp=='D'.or.c80tmp=='*').and.RICnatm(iRIC)==4) then
										if (icrit==1) then
											if (any(tmparr==RICatm(1,iRIC)).and.any(tmparr==RICatm(2,iRIC)).and.any(tmparr==RICatm(3,iRIC)).and.any(tmparr==RICatm(4,iRIC))) iadd=1
										else if (icrit==2) then
											if (any(tmparr==RICatm(1,iRIC)).or.any(tmparr==RICatm(2,iRIC)).or.any(tmparr==RICatm(3,iRIC)).or.any(tmparr==RICatm(4,iRIC))) iadd=1
                                        end if
									end if
                                    if (iadd==1) then
										if (any(PVSterm(1:PVSnterm(ifrag),ifrag)==iRIC)) cycle
										PVSnterm(ifrag)=PVSnterm(ifrag)+1
                                        PVSterm(PVSnterm(ifrag),ifrag)=iRIC
                                        numadd=numadd+1
                                    end if
                                end do
								deallocate(tmparr)
                                write(*,"(i6,' RICs have been added')") numadd
                            else !Add RIC by directly inputting atoms
                                ntmp=strcharnum(c80tmp,',')
                                if (ntmp==1) read(c80tmp,*) iatm,jatm
                                if (ntmp==2) read(c80tmp,*) iatm,jatm,katm
                                if (ntmp==3) read(c80tmp,*) iatm,jatm,katm,latm
                                do iRIC=1,nRIC
									i1=RICatm(1,iRIC);i2=RICatm(2,iRIC);i3=RICatm(3,iRIC);i4=RICatm(4,iRIC)
                                    iadd=0
									if (ntmp==1.and.RICnatm(iRIC)==2) then
										if ( (i1==iatm.and.i2==jatm).or.(i2==iatm.and.i1==jatm) ) iadd=1
									else if (ntmp==2.and.RICnatm(iRIC)==3) then
										if ( (i1==iatm.and.i2==jatm.and.i3==katm).or.(i3==iatm.and.i2==jatm.and.i1==katm) ) iadd=1
									else if (ntmp==3.and.RICnatm(iRIC)==4) then
										if ( (i1==iatm.and.i2==jatm.and.i3==katm.and.i4==latm).or.(i4==iatm.and.i3==jatm.and.i2==katm.and.i1==latm) ) iadd=1
                                    end if
                                    if (iadd==1) then
										PVSnterm(ifrag)=PVSnterm(ifrag)+1
                                        PVSterm(PVSnterm(ifrag),ifrag)=iRIC
                                        write(*,"(1x,a,' has been added to this fragment')") trim(RICname(iRIC))//' '//trim(RICdefname(iRIC))
										exit
                                    end if
                                end do
                                if (iRIC==nRIC+1) write(*,*) "No redundant internal coordinate was added"
                            end if
                        end do
                    end if
				end if
            end if
		end do
        
    end if
	
	if (isel==15.and.nsystem>1) then !Showing individual transition contribution is not possible when multiple files are involved
		write(*,*) "Error: This function is not available when multiple files are involved!"
		write(*,*) "Press ENTER button to continue"
		read(*,*)
		cycle
	end if
	
    
    
    
    

	!!=======================================================================!!
	!!=======================================================================!!
	!!=============== Below functions need calculation of curves ============!!
	!!=======================================================================!!
	!!=======================================================================!!
	if (isel==0.or.isel==1.or.isel==2.or.isel==15) then
		!====== Construct correspondence array if outputting individual bands. Only available when one file is loaded
		!This function is not available when multiple systems are considered
		if (isel==15) then
			write(*,"(a)") " Input criterion of strength, e.g. 0.2, then the contribution curves of the transitions &
            whose absolute strength larger than it will be exported to plain text files"
            write(*,"(a)") " If you input 0 now, then 10 maximum contributions to an inputted X position will be shown on screen"
			read(*,*) critindband
            if (critindband==0) then
                write(*,*) "Input the X position in current unit, e.g. 435.9"
                read(*,*) decompXpos
                if (iunitx==2) decompXpos=1239.842D0/decompXpos !Convert the inputted value in nm to internal unit eV
                allocate(indcontri(numdata))
                indcontri=0
            else
			    numindband=0 !The number of individual bands satisfying criterion
			    do idata=1,numdata
				    if (abs(strall(1,idata))>=critindband) numindband=numindband+1
			    end do
			    allocate(indband2idx(numindband),idx2indband(numdata),indcurve(num1Dpoints,numindband))
			    indcurve=0
			    itmp=0
			    idx2indband=0
			    do idata=1,numdata
				    if (abs(strall(1,idata))<critindband) cycle
				    itmp=itmp+1
				    indband2idx(itmp)=idata !Map index of outputted bands (indband) to actual index of all transitions (idx)
				    idx2indband(idata)=itmp !If not =0, the contribution from i band will be stored to idx2indband(i) slot of indcurve array
			    end do
            end if
		end if
		!====== Determine upper and lower limit of X axis =======
		if (iusersetx==0) then !Automatical scale, if user has not manually set the range
			tmpmin=minval(dataxall(1,1:numdataall(1)))
			tmpmax=maxval(dataxall(1,1:numdataall(1)))
			if (nsystem>1) then !Find upper and lower values for all systems
				do imol=2,nsystem
					tmpa=minval(dataxall(imol,1:numdataall(imol)))
					tmpb=maxval(dataxall(imol,1:numdataall(imol)))
					if (tmpa<tmpmin) tmpmin=tmpa
					if (tmpb>tmpmax) tmpmax=tmpb
				end do
			end if
			if (iunitx==0) then !cm^-1 for IR, Raman, VCD
				xlow=4000D0 !In common spectrum the energy is from high to low
				xhigh=0D0
			else if (iunitx==2) then !nm for UV-Vis, ECD, generate proper range in eV
				xhigh=1239.842D0/ (1239.842D0/tmpmin+40) !Note that when nm is used, in common spectrum the energy is from high to low, so we invert xlow and xhigh
				xlow=1239.842D0/ (1239.842D0/tmpmax-40)
			else
				rangetmp=tmpmax-tmpmin
				if (rangetmp==0) rangetmp=abs(dataxall(1,1)) !Only one data
				xlow=tmpmin-0.3D0*rangetmp
				xhigh=tmpmax+0.3D0*rangetmp
			end if
		else if (iusersetx==1) then !The range was defined by user
			!nm is selected unit, however we still using eV during broadening, so we convert user inputted range from nm to eV, after broadening then convert back
			if (iunitx==2) then
				xlow=1239.842D0/xlow
				xhigh=1239.842D0/xhigh
			end if
		end if
		!====== Set x position of curves ==========
		if (iunitx==2) then !For nm, which is not a linear energy unit, we generate points evenly distributed in X-axis with nm as unit
			xhighnm=1239.842D0/xhigh !nm->eV
			xlownm=1239.842D0/xlow
			xptstep=(xhighnm-xlownm)/(num1Dpoints-1) !Get proper spacing in nm
			do ipoint=1,num1Dpoints
				curvex(ipoint)=1239.842D0/(xlownm+(ipoint-1)*xptstep) !Now curvex is recorded in eV
			end do
		else !Common case, linear energy unit is used
			xptstep=(xhigh-xlow)/(num1Dpoints-1)
			do ipoint=1,num1Dpoints
				curvex(ipoint)=xlow+(ipoint-1)*xptstep
			end do
		end if
		!====== Generate energy levels line =======
		lineyall=0D0
		linexall=1D0 !To garantee that linexall will not be zero, otherwise may crash when converting unit
		do imol=1,nsystem
			do idata=1,numdataall(imol)
				inow=3*(idata-1)
				linexall(imol,inow+1:inow+3)=dataxall(imol,idata)
				lineyall(imol,inow+2)=weight(imol)*strall(imol,idata) !Line height is weighted! Otherwise plotting them is meaningless
			end do
		end do
		!===============================================
		!============== Broaden to curve ===============
		!===============================================
		!Under current X and Y axes units, below code guarantees that the integral of the peak broadened by one unit of strength is 1
		curveyall=0D0
        PVScurve=0D0
        OPVScurve=0D0
        PVScurveintra=0
		if (ibroadfunc==1.or.ibroadfunc==3) then !Lorentzian function or Pseudo-Voigt function, see http://mathworld.wolfram.com/LorentzianFunction.html
			do imol=1,nsystem
				do idata=1,numdataall(imol) !Cycle each transition
					preterm=strall(imol,idata)*0.5D0/pi*FWHMall(imol,idata) !Integral of the peak equals to str(idata)
					do ipoint=1,num1Dpoints
						curveytmp(ipoint)=preterm/( (curvex(ipoint)-dataxall(imol,idata))**2+0.25D0*FWHMall(imol,idata)**2 )
					end do
					curveyall(imol,:)=curveyall(imol,:)+curveytmp(:)
                    !Partial vibrational spectrum
                    do ifrag=1,maxPVSfrag
						if (PVSnterm(ifrag)/=0) PVScurve(ifrag,:)=PVScurve(ifrag,:)+curveytmp(:)*PVScomp(ifrag,idata)
                        if (iPVScomptype==2.and.ishowlocPVS==1) PVScurveintra(ifrag,:)=PVScurveintra(ifrag,:)+curveytmp(:)*PVScompintra(ifrag,idata)
                    end do
                    if ((OPVSidx1/=0.or.OPVSidx2/=0)) OPVScurve(:)=OPVScurve(:)+curveytmp(:)*OPVScomp(idata)
					if (isel==15) then !Individual contribution
                        if (critindband==0) then
                            indcontri(idata)=preterm/( (decompXpos-dataxall(imol,idata))**2+0.25D0*FWHMall(imol,idata)**2 )
                        else
						    if (idx2indband(idata)/=0) indcurve(:,idx2indband(idata))=curveytmp
                        end if
					end if
				end do
			end do
		end if
		if (ibroadfunc==2.or.ibroadfunc==3) then !Gaussian or Pseudo-Voigt function, see http://en.wikipedia.org/wiki/Gaussian_function
			if (ibroadfunc==3) then
				curveyall=(1-gauweigh)*curveyall !Scale Lorentzian function part of Pseudo-Voigt
				indcurve=(1-gauweigh)*indcurve
                indcontri=(1-gauweigh)*indcontri
			end if
			do imol=1,nsystem
				do idata=1,numdataall(imol)
					gauss_c=FWHMall(imol,idata)/2D0/sqrt(2*dlog(2D0))
					gauss_a=strall(imol,idata)/(gauss_c*sqrt(2D0*pi))
					do ipoint=1,num1Dpoints
						curveytmp(ipoint)=gauss_a*dexp( -(curvex(ipoint)-dataxall(imol,idata))**2/(2*gauss_c**2) )
					end do
					if (ibroadfunc==3) curveytmp(:)=gauweigh*curveytmp(:) !Scale Gaussian function part of Pseudo-Voigt
					curveyall(imol,:)=curveyall(imol,:)+curveytmp(:)
                    !Partial vibrational spectrum
                    do ifrag=1,maxPVSfrag
						if (PVSnterm(ifrag)/=0) PVScurve(ifrag,:)=PVScurve(ifrag,:)+curveytmp(:)*PVScomp(ifrag,idata)
                        if (iPVScomptype==2.and.ishowlocPVS==1) PVScurveintra(ifrag,:)=PVScurveintra(ifrag,:)+curveytmp(:)*PVScompintra(ifrag,idata)
                    end do
                    if (.not.(OPVSidx1==0.and.OPVSidx2==0)) OPVScurve(:)=OPVScurve(:)+curveytmp(:)*OPVScomp(idata)
					if (isel==15) then !Individual contribution
                        if (critindband==0) then
                            !write(*,*) gauss_a*dexp( -(decompXpos-dataxall(imol,idata))**2/(2*gauss_c**2) )
						    indcontri(idata)=indcontri(idata)+gauss_a*dexp( -(decompXpos-dataxall(imol,idata))**2/(2*gauss_c**2) )
                        else
						    if (idx2indband(idata)/=0) indcurve(:,idx2indband(idata))=indcurve(:,idx2indband(idata))+curveytmp
                        end if
					end if
				end do
			end do
		end if

		!Change units, scale curve, set axis
		if (ispectrum==1.and.iunitliney==2) lineyall=lineyall/2.5066D0 !For IR spectrum, convert strength unit from km/mol to esu^2*cm^2
		if (iunitx==2) then !eV->nm
			linexall=1239.842D0/linexall
			curvex=1239.842D0/curvex
			xlow=1239.842D0/xlow
			xhigh=1239.842D0/xhigh
		end if
		curveyall(:,:)=scalecurve*curveyall(:,:)
        PVScurve(:,:)=scalecurve*PVScurve(:,:)
        if (iPVScomptype==2.and.ishowlocPVS==1) PVScurveintra(:,:)=scalecurve*PVScurveintra(:,:)
        OPVScurve(:)=scalecurve*OPVScurve(:)
		if (isel==15) then !Scale individual contributions
            if (critindband==0) then
                indcontri(:)=scalecurve*indcontri(:)
            else
                indcurve(:,:)=scalecurve*indcurve(:,:)
            end if
        end if
		curvex(:)=curvex(:)+shiftx
		linexall(:,:)=linexall(:,:)+shiftx
		if (iusersetx==0) stepx=(xhigh-xlow)/10
		
		!Generate weighted curve from multiple curves
		curvey=0
		do imol=1,nsystem
			curvey(:)=curvey(:)+weight(imol)*curveyall(imol,:)
		end do
		!Weighting spectrum of each system
		if (iweisyscurve==1) then
			do imol=1,nsystem
				curveyall(imol,:)=weight(imol)*curveyall(imol,:)
			end do
		end if
	end if

	!================================================
	!======= Output data to external text file ======
	!================================================
    if (isel==15.and.critindband==0) then !Show contribution at the given position and then return
        !Sort from large to small
		allocate(tmparr(numdata))
		forall(itmp=1:numdata) tmparr(itmp)=itmp
        call sortr8(indcontri,"abs",tmparr)
        call invarrr8(indcontri,tmparr)
        totval=sum(abs(indcontri))
        write(*,"(a,f20.5)") " Sum of absolute values of all transitions:",totval
        write(*,*) "The individual terms are ranked by magnitude of contribution:"
        write(*,*) "  #Transition     Contribution      %"
        do itmp=1,min(10,numdata)
            write(*,"(i10,f20.5,f11.3)") tmparr(itmp),indcontri(itmp),indcontri(itmp)/totval*100
        end do
        if (ispectrum==4.or.ispectrum==5.or.ispectrum==6) then
            write(*,"(a)") " Note: The % shown above is calculated via dividing &
            absolute value of a transition by sum of absolute values of all transitions"
        end if
        write(*,*)
        deallocate(indcontri,tmparr)
        cycle
    end if
	if (isel==2.or.isel==15) then !Output curve for total and individual contributions, respectively
		if (isel==2) then !Output regular spectrum
			if (nsystem==1) then
				open(10,file="spectrum_curve.txt",status="replace")
				do ipt=1,num1Dpoints
					write(10,"(f13.5,1PE18.8E3)") curvex(ipt),curvey(ipt)
				end do
				close(10)
				write(*,*) "Curve data has been written to spectrum_curve.txt in current folder"
			else !Also output curve for all systems
				if (any(weight/=1)) then !Output weighted spectrum 
					open(10,file="spectrum_curve.txt",status="replace")
					do ipt=1,num1Dpoints
						write(10,"(f13.5,1PE18.8E3)") curvex(ipt),curvey(ipt)
					end do
					close(10)
					write(*,"(a)") " The curve data corresponding to weighted spectrum has been written to spectrum_curve.txt in current folder"
				end if
				open(10,file="curveall.txt",status="replace")
				do ipt=1,num1Dpoints
					write(10,"(f13.5)",advance="no") curvex(ipt)
					do imol=1,nsystem
						write(10,"(1PE18.8E3)",advance="no") curveyall(imol,ipt)
					end do
					write(10,*)
				end do
				close(10)
				write(*,"(a)") " Curve data of all systems have been exported to curveall.txt in current folder as different columns"
			end if
            if (any(PVSnterm/=0)) then
				!Output PVS curve
				open(10,file="PVScurve.txt",status="replace")
                write(10,"(a)",advance="no") "#      X     "
				do ifrag=1,maxPVSfrag !Write header
					if (PVSnterm(ifrag)/=0) then
						write(c80tmp,"(i3)") ifrag
                        c80tmp=adjustl(c80tmp)
                        write(10,"('      Frag',a,'      ')",advance="no") c80tmp(1:2)
                    end if
				end do
                write(10,*)
                do ipt=1,num1Dpoints
                    write(10,"(f13.5)",advance="no") curvex(ipt)
					do ifrag=1,maxPVSfrag
						if (PVSnterm(ifrag)/=0) write(10,"(1PE18.8E3)",advance="no") PVScurve(ifrag,ipt)
					end do
                    write(10,*)
                end do
                close(10)
                write(*,"(a)") " Curves of PVS spectra have been exported to PVScurve.txt in current folder, see first line of this file for meaning"
				!Output OPVS curve
				if (OPVSidx1/=0.and.PVSnterm(OPVSidx1)/=0.and.OPVSidx2/=0.and.PVSnterm(OPVSidx2)/=0) then
					open(10,file="OPVScurve.txt",status="replace")
					do ipt=1,num1Dpoints
						write(10,"(f13.5,1PE18.8E3)") curvex(ipt),OPVScurve(ipt)
					end do
					close(10)
					write(*,"(a)") " Curve of OPVS has been exported to OPVScurve.txt in current folder"
				end if
            end if
		else if (isel==15) then !Output individual band contributions
			open(10,file="spectrum_curve.txt",status="replace")
			do ipt=1,num1Dpoints
				write(10,"(f13.5)",advance="no") curvex(ipt)
				write(10,"(1PE18.8E3)",advance="no") curvey(ipt)
				do iindband=1,numindband
					write(10,"(1PE18.8E3)",advance="no") indcurve(ipt,iindband)
				end do
				write(10,*)
			end do
			close(10)
			write(*,"(a,i5,a)") " The total spectrum and the contributions from",numindband," transitions have been outputted to spectrum_curve.txt in current folder"
		end if
		!Explain meaning of each column
		if (ispectrum==1) then !IR
			write(*,*) "Column 1: Wavenumber (cm^-1)"
            if (nsystem==1) then
				write(*,*) "Column 2: Molar absorption coefficient (L/mol/cm)"
            else
				write(*,*) "Other columns: Molar absorption coefficient (L/mol/cm) of each curve"
            end if
		else if (ispectrum==2) then !Raman
			write(*,*) "Column 1: Wavenumber (cm^-1)"
            if (nsystem==1) then
				write(*,*) "Column 2: Relative Raman intensity"
            else
				write(*,*) "Other columns: Relative Raman intensity of each curve"
            end if
		else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
			if (iunitx==1) write(*,*) "Column 1: Excitation energy (eV)"
			if (iunitx==2) write(*,*) "Column 1: Wavelength (nm)"
			if (iunitx==3) write(*,*) "Column 1: Wavenumber (1000 cm^-1)"
            if (nsystem==1) then
				if (ispectrum==3) write(*,*) "Column 2: Molar absorption coefficient (L/mol/cm)"
				if (ispectrum==4) write(*,*) "Column 2: Delta molar absorption coefficient (arb.)"
            else
				if (ispectrum==3) write(*,*) "Other columns: Molar absorption coefficient (L/mol/cm) of each curve"
				if (ispectrum==4) write(*,*) "Other columns: Delta molar absorption coefficient (arb.) of each curve"
            end if
		else if (ispectrum==5) then !VCD
			write(*,*) "Column 1: Wavenumber (cm^-1)"
            if (nsystem==1) then
				write(*,*) "Column 2: Delta molar absorption coefficient (L/mol/cm)"
            else
				write(*,*) "Other columns: Delta molar absorption coefficient (L/mol/cm) of each curve"
            end if
		else if (ispectrum==6) then !ROA
			write(*,*) "Column 1: Wavenumber (cm^-1)"
            if (nsystem==1) then
				if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) write(*,*) "Column 2: IR+IL"
				if (iROAtype==2.or.iROAtype==4.or.iROAtype==6) write(*,*) "Column 2: IR-IL"
            else
				if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) write(*,*) "Other columns: IR+IL of each curve"
				if (iROAtype==2.or.iROAtype==4.or.iROAtype==6) write(*,*) "Other columns: IR-IL of each curve"
            end if
		end if
		if (isel==15) then
			write(*,*) "Correspondence between the columns and individual transitions in the file:"
			write(*,*) "     Column#   Transition#"
			do itmp=1,numindband
				write(*,"(2i12)") itmp+2,indband2idx(itmp)
			end do
			deallocate(indband2idx,idx2indband,indcurve)
		end if
		
		!Output discrete lines
		open(10,file="spectrum_line.txt",status="replace")
		do imol=1,nsystem
			do ipt=1,3*numdataall(imol)
				write(10,"(f16.5,1PE18.8E3)") linexall(imol,ipt),lineyall(imol,ipt)
			end do
			write(10,*)
		end do
		close(10)
		write(*,*)
		if (nsystem>1) then
			write(*,"(a)") " Discrete line data of all systems have been written together to spectrum_line.txt &
            in current folder, data of each system is separated by a blank line"
			if (any(weight/=1)) write(*,*) "Note: The height of discrete lines in this file have been weighted"
		else
			write(*,*) "Discrete line data have been written to spectrum_line.txt in current folder"
		end if
		if (ispectrum==1) then !IR
			write(*,*) "Column 1: Frequency (cm^-1)"
			if (iunitliney==1) write(*,*) "Column 2: IR intensities (km/mol)"
			if (iunitliney==2) write(*,*) "Column 2: IR intensities (esu^2*cm^2)"
		else if (ispectrum==2) then !Raman
			write(*,*) "Column 1: Frequency (cm^-1)"
			if (iramantype==1) write(*,*) "Column 2: Raman scattering activities (A^4/amu)"
			if (iramantype==2) write(*,*) "Column 2: Raman scattering intensities"
		else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
			if (iunitx==1) write(*,*) "Column 1: Excitation energy (eV)"
			if (iunitx==2) write(*,*) "Column 1: Wavelength (nm)"
			if (iunitx==3) write(*,*) "Column 1: Wavenumber (1000 cm^-1)"
			if (ispectrum==3) write(*,*) "Column 2: Oscillator strength"
			if (ispectrum==4) write(*,*) "Column 2: Rotatory strength in cgs (10^-40 erg-esu-cm/Gauss)"
		else if (ispectrum==5) then !VCD
			write(*,*) "Column 1: Frequency (cm^-1)"
			write(*,*) "Column 2: Rotatory strength (10^-44 esu^2 cm^2)"
		else if (ispectrum==6) then !ROA
			write(*,*) "Column 1: Frequency (cm^-1)"
			if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) write(*,*) "Column 2: Raman intensity (K)"
			if (iROAtype==2.or.iROAtype==4.or.iROAtype==6) write(*,*) "Column 2: ROA intensity (10^4 K)"
		end if
	end if	
	
    !Find and print minimum/maximum positions
    if (nsystem==1.or.any(weight/=1)) then !If simply plot multiple systems, extrema are meaningless
        if (isel==0) write(*,*) "Extrema on the spectrum curve:"
	    numlocmax=0
	    do ipoint=2,num1Dpoints-1
		    gradold=curvey(ipoint)-curvey(ipoint-1)
		    gradnew=curvey(ipoint+1)-curvey(ipoint)
		    if (gradold*gradnew<0D0.and.gradold>gradnew) then
			    numlocmax=numlocmax+1
                maxlabX(numlocmax)=ipoint
			    if (isel==0) write(*,"(' Maximum',i5,'   X:',f15.4,'   Value:',f15.4)") numlocmax,curvex(ipoint),curvey(ipoint)
		    end if
	    end do
	    numlocmin=0
	    do ipoint=2,num1Dpoints-1
		    gradold=curvey(ipoint)-curvey(ipoint-1)
		    gradnew=curvey(ipoint+1)-curvey(ipoint)
		    if (gradold*gradnew<0D0.and.gradold<gradnew) then
			    numlocmin=numlocmin+1
                minlabX(numlocmin)=ipoint
                if (ispectrum>2.and.isel==0) then !UV-Vis, ECD, VCD, ROA
			        write(*,"(' Minimum',i5,'   X:',f15.3,'   Value:',f15.4)") numlocmin,curvex(ipoint),curvey(ipoint)
                end if
		    end if
	    end do
        if (isel==0) then
			if (ispectrum<=2) write(*,"(a)") " Position of minima are not reported since they are commonly not of interest for this kind of spectrum"
			if (nsystem>1) write(*,*) "Note: The extrema reported above correspond to weighted spectrum"
        end if
        if (iPVScomptype==2.and.ishowlocPVS==1) write(*,"(a)") " Note: The dashed curves in the map correspond to intrafragment contribution to PVS-I spectrum"
    end if
    
    
	!========================================
	!============ Draw spectrum =============
	!========================================
	if (isel==0.or.isel==1) then

		if (iusersetY1==0) then !Set default lower and upper limit of left Y axis
			endy1=1.1D0*max(maxval(abs(curvey)),maxval(abs(curveyall)))
			if (nsystem>1.and.all(weight==1)) endy1=1.1D0*maxval(abs(curveyall))
			orgy1=-endy1/30D0 !Slightly lower it to avoid curve touching bottom
			!Positive and negative of VCD, ECD and ROA curves may have similar height
			if (ispectrum==4.or.ispectrum==5.or.(ispectrum==6.and.(iROAtype==2.or.iROAtype==4.or.iROAtype==6))) orgy1=-endy1
			stepy1=(endy1-orgy1)/10
		end if
		if (iusersetY2==0) then !Set default lower and upper limit of right Y axis
			endy2=1.1D0*maxval(abs(lineyall))
			orgy2=-endy2/30D0
			if (ispectrum==4.or.ispectrum==5.or.(ispectrum==6.and.(iROAtype==2.or.iROAtype==4.or.iROAtype==6))) orgy2=-endy2
			stepy2=(endy2-orgy2)/10
		end if
		if (isilent==1.and.isavepic==0) cycle
		
		if (isavepic==0) then
			call METAFL('xwin')
			call window(200,100,1200,770)
		else if (isavepic==1) then
			call METAFL(graphformat)
			call window(200,100,graph1Dwidth,graph1Dheight)
		end if
		call SCRMOD('REVERSE')
		CALL IMGFMT("RGB")
		CALL PAGE(3000,1800)
        if (isavepic==1) CALL PAGE(3050,1900) !Make page larger to avoid truncation when text size is set to relatively large
		call disini
		if (isavepic==1.and.graphformat=="png") then
			call TRIPLX 
		else
			CALL HWFONT
		end if
        nxpixel=2150
		if (ishowlevel==0) nypixel=1500
        if (ishowlevel==1) nypixel=1400 !When showing spikes, compress the spectrum region to leave space for plotting spikes
        call AXSLEN(nxpixel,nypixel)
        if (ishowline==1) then
		    if (ishowlevel==0) call axspos(400,1640) !Leave more space at right side to show axis
		    if (ishowlevel==1) call axspos(400,1540)
        else
		     if (ishowlevel==0) call axspos(510,1640)
		     if (ishowlevel==1) call axspos(510,1540)
        end if
! 		call center
		if (isavepic==0) call WINTIT("Click right mouse button to close")
		if (ishowlevel==0) CALL TICKS(1,'XY')
		if (ishowlevel==1) CALL TICKS(0,'X')
		call ERRMOD("ALL","OFF")
        call HNAME(height_axis) !Height of axis name
        call height(ticksize) !Size of ticks
        if (ndecimalX==-1) then
		    CALL LABDIG(1,"X")
		    if (iunitx==1) CALL LABDIG(2,"X") !eV, use more digits
        else
            if (ndecimalX==0) then
                CALL LABDIG(-1,"X") !Only show integer
            else
    		    CALL LABDIG(ndecimalX,"X")
            end if
        end if
        
		!Set name of X-axis
        if (ishowlevel==0) then !Do not use spikes at bottom to show specific levels, so show axis name here
		    if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then
			    call TEXMOD("ON")
			    CALL NAME('Wavenumber (cm$^{-1}$)','X')
		    end if
		    if (iunitx==1) CALL NAME('Excitation energy (eV)','X')
		    if (iunitx==2) CALL NAME('Wavelength (nm)','X')
		    if (iunitx==3) CALL NAME('Wavenumber (1000 cm$^{-1}$)','X')
        end if
        
		!Set name of Y-axis
		call TEXMOD("ON")
		if (ispectrum==1.or.ispectrum==3) then
			CALL NAME('Molar absorption coefficient  $\epsilon$ (L mol$^{-1}$cm$^{-1}$)','Y')
		else if (ispectrum==2) then
			CALL NAME('Relative Raman intensity','Y')
		else if (ispectrum==4.or.ispectrum==5) then
			CALL NAME('$\Delta\epsilon$ (arb.)','Y')
		else if (ispectrum==6) then
			if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) then
				CALL NAME('$I_R+I_L$','Y')
			else if (iROAtype==2.or.iROAtype==4.or.iROAtype==6) then
				CALL NAME('$I_R-I_L$','Y')
			end if
		end if
        if (iVDOS==1) CALL NAME('Vibrational density-of-states (arb. unit)','Y')
        
		!Set plotting parameters of axes and legends
		tmprange1=abs(endy1-orgy1)
        if (ndecimalYleft==-1) then
		    if (tmprange1>5) then
			    CALL LABDIG(1,"Y")
		    else if (tmprange1>0.5D0) then
			    CALL LABDIG(2,"Y")
		    else if (tmprange1>0.05D0) then
			    CALL LABDIG(3,"Y")
		    else
			    CALL LABDIG(4,"Y")
		    end if
        else
            if (ndecimalYleft==0) then
                CALL LABDIG(-1,"Y") !Only show integer
            else
    		    CALL LABDIG(ndecimalYleft,"Y")
            end if
        end if
		if (ishowline==1) then
			if (ishowlevel==0) call setgrf('NAME',"NAME",'TICKS','NONE') !If show discrete lines, leave right axis empty
            if (ishowlevel==1) call setgrf('TICKS',"NAME",'TICKS','NONE')
		else
			if (ishowlevel==0) call setgrf('NAME',"NAME",'TICKS','TICKS')
			if (ishowlevel==1) call setgrf('TICKS',"NAME",'TICKS','TICKS')
		end if
        if (ishowlabelleft==0) then
            call labels("NONE","Y")
            call namdis(60,'Y') !Use larger distance between name and axis
        else
            if (labtype_Yleft==1) call labels("FLOAT","Y")
            if (labtype_Yleft==2) call labels("XEXP","Y")
            if (labtype_Yleft==3) call labels("FEXP","Y")
        end if
		ileg=0 !Initialize legend information
        numleg=1
        if (nsystem>1) then
			numleg=1+nsystem
        else if (any(PVSnterm/=0)) then
			do ifrag=1,maxPVSfrag
				if (PVSnterm(ifrag)>0.and.iPVSshow(ifrag)==1) numleg=numleg+1
            end do
            if (OPVSidx1/=0.and.OPVSidx2/=0.and.iOPVSshow==1) numleg=numleg+1
        end if
		call legini(clegend,numleg,50)
		call legtit(' ')
		call frame(0) !No box around legend
		
        !Set axis and draw curves
		call LINWID(thk_axis)
		CALL GRAF(xlow+shiftx,xhigh+shiftx,xlow+shiftx,stepx, orgy1,endy1,orgy1,stepy1)
		if (ishowgrid==1) then !Draw shallow gray dashed lines
			call SETRGB(0.8D0,0.8D0,0.8D0) !Shallow gray
			call dash
			call LINWID(thk_grid)
			call grid(1,1)
		end if
		call solid
		if (ishowweicurve==1.or.ishowweicurve==2) then !Draw weighted curve (including only one system case) or total curve of PVS
			ileg=ileg+1
			call setcolor(icurveclr)
			CALL LINWID(thk_weighted) !Use very thick line for weighted curve when there are multiple systems
			if (nsystem==1) CALL LINWID(thk_curve)
			CALL CURVE(curvex,curvey,num1Dpoints)
			call legpat(0,1,-1,-1,-1,ileg)
            if (any(PVSnterm/=0)) then
				CALL LEGLIN(clegend,"Total",ileg)
            else
				CALL LEGLIN(clegend," Weighted",ileg)
            end if
		end if
		if (nsystem>1.and.ishowweicurve/=2) then !Draw curve for each system, and meantime set corresponding legend
			do imol=1,nsystem
                if (ishowweicurve==1) then !Conformation weighted and individual spectra
				    iclrtmp=imol+1 !The 1st color is red, which has been used by weighted curve
				    if (iclrtmp==4) iclrtmp=10 !4 corresponds to black, however, due to reverse, it corresponds to white, which is unable to use
				    if (iclrtmp==2) then !2 corresponds to green, which is too bright
					    iclrtmp=12 !Change to dark green
				    else if (iclrtmp==12) then
					    iclrtmp=2
				    end if
                else !Simply draw multiple system spectra together
                    iclrtmp=currclr(imol)
                    if (iclrtmp>ncurrclr) iclrtmp=ncurrclr
                end if
				call setcolor(iclrtmp)
				CALL LINWID(thk_curve)
				CALL CURVE(curvex,curveyall(imol,:),num1Dpoints)
				ileg=ileg+1
				call legpat(0,1,-1,-1,-1,ileg)
				CALL LEGLIN(clegend,trim(mollegend(imol)),ileg)
			end do
		end if
        
        !Draw (overlap) partial vibrational spectra
		do ifrag=1,maxPVSfrag
			if (PVSnterm(ifrag)>0.and.iPVSshow(ifrag)==1) then
				call setcolor(PVScolor(ifrag))
				CALL LINWID(thk_PVS)
				CALL CURVE(curvex,PVScurve(ifrag,:),num1Dpoints)
				ileg=ileg+1
				call legpat(0,1,-1,-1,-1,ileg)
				CALL LEGLIN(clegend,trim(PVSlegend(ifrag)),ileg)
				if (iPVScomptype==2.and.ishowlocPVS==1) then
					call dash
					CALL CURVE(curvex,PVScurveintra(ifrag,:),num1Dpoints)
                    call solid
                end if
            end if
        end do
        if (OPVSidx1/=0.and.OPVSidx2/=0.and.iOPVSshow==1) then
			call setcolor(OPVScolor)
            CALL LINWID(thk_OPVS)
			CALL CURVE(curvex,OPVScurve(:),num1Dpoints)
			ileg=ileg+1
			call legpat(0,1,-1,-1,-1,ileg)
			CALL LEGLIN(clegend,"OPVS",ileg)
        end if
        
		call color("WHITE")
		if (iYeq0==1) then
            CALL LINWID(thk_Yeq0)
			call xaxgit !Draw a line corresponding to Y=0
        end if
		call box2d !The dashed line of "call grid" overlaied frame of axis, so redraw frame
		call legopt(2.5D0,0.5D0,1D0) !Decrease the length of legend color line
        call height(legtextsize) !Define legend text size
		if (nsystem>1.and.ishowweicurve/=2) call legend(clegend,ilegendpos)
		if (any(PVSnterm/=0)) call legend(clegend,ilegendpos)
        call height(ticksize) !Size of ticks
		call endgrf
        call labels("FLOAT","Y") !Restore to default
		
        !Set axis and draw discrete lines
		if (ishowline==1) then
			if (ispectrum==1) then
				if (iunitliney==1) CALL NAME('IR intensities (km mol$^{-1}$)','Y')
				if (iunitliney==2) CALL NAME('IR intensities (esu$^2$ cm$^2$)','Y')
			else if (ispectrum==2) then
				if (iramantype==1) CALL NAME('Raman scattering activities (A$^4$ amu$^{-1}$)','Y')
				if (iramantype==2) CALL NAME('Raman scattering intensities','Y')
			else if (ispectrum==3) then
				call name("Oscillator strength",'Y')
			else if (ispectrum==4) then
				call name("Rotatory strength (cgs)",'Y')
			else if (ispectrum==5) then
				call name("Rotatory strength",'Y')
			else if (ispectrum==6) then
				if (iramantype==1) then !Directly loaded from Gaussian output file, formally they are known as intensity
					!Raman SCP(180) loaded from Gaussian output file is simply 4 times of pre-resonance Raman activity, clearly it is not "real" intensity
					if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) call name("Raman intensity (K)",'Y')
					if (iROAtype==2.or.iROAtype==4.or.iROAtype==6) call name("ROA intensity (10$^4$ K)",'Y')
				else if (iramantype==2) then !Convert to "real" intensity
					if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) call name("Raman intensity",'Y')
					if (iROAtype==2.or.iROAtype==4.or.iROAtype==6) call name("ROA intensity",'Y')
				end if
			end if
			tmprange2=abs(endy2-orgy2)
            if (ndecimalYright==-1) then
			    if (tmprange2>5) then
				    CALL LABDIG(1,"Y")
			    else if (tmprange2>0.5D0) then
				    CALL LABDIG(2,"Y")
			    else if (tmprange2>0.05D0) then
				    CALL LABDIG(3,"Y")
			    else
				    CALL LABDIG(4,"Y")
			    end if
            else
                if (ndecimalYright==0) then
                    CALL LABDIG(-1,"Y") !Only show integer
                else
    		        CALL LABDIG(ndecimalYright,"Y")
                end if
            end if
            call setgrf('NONE','NONE','NONE','NAME')
            if (ishowlabelright==0) then
                call labels("NONE","Y")
                call namdis(60,'Y')
            end if
			call LINWID(thk_axis)
			CALL GRAF(xlow+shiftx,xhigh+shiftx,xlow+shiftx,stepx, orgy2,endy2,orgy2,stepy2)
			CALL LINWID(thk_discrete)
            if (nsystem==1) then
			    call setcolor(ilineclr)
                imol=1
                CALL CURVE(linexall(imol,1:3*numdataall(imol)),lineyall(imol,1:3*numdataall(imol)),3*numdataall(imol))
            else
			    do imol=1,nsystem
                    if (any(weight/=1)) then !Conformation weighted spectrum
				        if (iweisyscurve==1) then !Individual spectral curves are weighted
					        iclrtmp=imol+1 !The 1st color is red, which has been used by weighted curve
					        if (iclrtmp==4) iclrtmp=10 !4 corresponds to black, however, due to reverse, it corresponds to white, which is unable to use
					        if (iclrtmp==2) then !2 corresponds to green, which is too bright
					        	iclrtmp=12 !Change to dark green
					        else if (iclrtmp==12) then
					        	iclrtmp=2
					        end if
					        call setcolor(iclrtmp)
                        else !Individual spectral curves are not weighted, all spikes share same color
			                call setcolor(ilineclr)
				        end if
                    else !Simply draw multiple conformation spectra together
                        iclrtmp=currclr(imol)
                        if (iclrtmp>ncurrclr) iclrtmp=ncurrclr
					    call setcolor(iclrtmp)
                    end if
				    CALL CURVE(linexall(imol,1:3*numdataall(imol)),lineyall(imol,1:3*numdataall(imol)),3*numdataall(imol))
			    end do
            end if
			call color("WHITE")
			call xaxgit !Draw a line corresponding to Y=0 using same thick as lines. This is important, otherwise the bottom line at the left side of the first spike will be invisible
            call endgrf
		end if
        
        !Show labels to indicate located extrema
        call setgrf('NONE','NONE','NONE','NONE')
		CALL GRAF(xlow+shiftx,xhigh+shiftx,xlow+shiftx,stepx, orgy1,endy1,orgy1,stepy1)
        pix2usrX=(xhigh-xlow)/nxpixel !used to shift label position in X
        pix2usrY=(endy1-orgy1)/nypixel !used to shift label position in Y
        call height(extlabelsize)
        call setcolor(extlabelclr)
        if (extlabeldecimal/=0) write(strfmt,"(a,i1,a)") "(f20.",extlabeldecimal,")"
        if (ishowextrema==1.or.ishowextrema==3) then !Show maxima labels
            if (iextlabelrot==1) call angle(90)
            do imax=1,numlocmax
                imaxpt=maxlabX(imax)
                if (extlabeldecimal==0) then
                    if (extlabelcontent==1) write(c200tmp,"(i20)") nint(curvex(imaxpt))
                    if (extlabelcontent==2) write(c200tmp,"(i20)") nint(curvey(imaxpt))
                else
                    if (extlabelcontent==1) write(c200tmp,strfmt) curvex(imaxpt)
                    if (extlabelcontent==2) write(c200tmp,strfmt) curvey(imaxpt)
                end if
                call rlmess(trim(adjustl(c200tmp)),curvex(imaxpt)+extmaxlabelshiftX*pix2usrX,curvey(imaxpt)+extmaxlabelshiftY*pix2usrY)
            end do
        end if
        if (ishowextrema==2.or.ishowextrema==3) then !Show minima labels
            if (iextlabelrot==1) call angle(-90)
            do imin=1,numlocmin
                iminpt=minlabX(imin)
                if (extlabeldecimal==0) then
                    if (extlabelcontent==1) write(c200tmp,"(i20)") nint(curvex(iminpt))
                    if (extlabelcontent==2) write(c200tmp,"(i20)") nint(curvey(iminpt))
                else
                    if (extlabelcontent==1) write(c200tmp,strfmt) curvex(iminpt)
                    if (extlabelcontent==2) write(c200tmp,strfmt) curvey(iminpt)
                end if
                call rlmess(trim(adjustl(c200tmp)),curvex(iminpt)+extminlabelshiftX*pix2usrX,curvey(iminpt)+extminlabelshiftY*pix2usrY)
            end do
        end if
        call setcolor(5) !Recover to black
        call height(36) !Recover to default
        call angle(0) !Recover to default
		call endgrf
        
        !Draw spikes at bottom to show all or specific transition levels, may or may not consider degeneracy
        if (ishowlevel==1) then
            call AXSLEN(nxpixel,90)
            if (ishowline==0) call axspos(510,1630)
            if (ishowline==1) call axspos(400,1630)
            CALL TICKS(1,'X')
            CALL TICKS(0,'Y')
            if (idegen==1) then !Determine maximum degeneracy of all transition levels
                maxdegen=1
                do idata=1,numdata
                    degentest=1
                    enei=dataxall(1,idata)
                    do jdata=idata+1,numdata
                        enej=dataxall(1,jdata)
                        if (abs(enej-enei)<degencrit) degentest=degentest+1
                    end do
                    if (degentest>maxdegen) maxdegen=degentest
                end do
                CALL TICKS(1,'Y')
            end if
            if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then
			    call TEXMOD("ON")
			    CALL NAME('Wavenumber (cm$^{-1}$)','X')
		    end if
		    if (iunitx==1) CALL NAME('Excitation energy (eV)','X')
		    if (iunitx==2) CALL NAME('Wavelength (nm)','X')
		    if (iunitx==3) CALL NAME('Wavenumber (1000 cm$^{-1}$)','X')
            call setgrf('NAME','TICKS','NONE','TICKS')
		    if (idegen==0) CALL GRAF(xlow+shiftx,xhigh+shiftx,xlow+shiftx,stepx, 0D0,1D0,0D0,1D0)
            if (idegen==1) CALL GRAF(xlow+shiftx,xhigh+shiftx,xlow+shiftx,stepx, 0D0,dfloat(maxdegen),0D0,1D0)
            allocate(spikey(3*numdata))
            call LINWID(spikethick)
            
            do iset=1,maxspike !Plotting all spike sets
                if (spikenum(iset)==0) cycle !This spike set has not been defined
                spikey=0
                if (idegen==0) then !Do not consider degeneracy
                    do idata=1,numdata
				        inow=3*(idata-1)
				        if (any(spikeidx(iset,1:spikenum(iset))==idata)) spikey(inow+2)=0.9D0
                    end do
                else !Consider degeneracy in current set
                    istart=1
                    do while(istart<spikenum(iset)) !Cycle all levels in this set, and compare with latter ones
                        do idata=istart,spikenum(iset)
                            ireal=spikeidx(iset,idata)
                            enei=dataxall(1,ireal)
                            degentest=1
                            do jdata=idata+1,spikenum(iset)
                                jreal=spikeidx(iset,jdata)
                                enej=dataxall(1,jreal)
                                if (abs(enej-enei)<degencrit) degentest=degentest+1
                            end do
				            inow=3*(ireal-1)
                            spikey(inow+2)=degentest
                        end do
                        istart=istart+degentest
                    end do
                end if
                call setcolor(spikecolor(iset))
			    CALL CURVE(linexall(1,1:3*numdata),spikey,3*numdata)
                if (iexportlevel==1) then !Export as .txt file
                    write(c200tmp,"('spike',i2.2,'.txt')") iset
                    open(10,file=c200tmp,status="replace")
                    do itmp=1,3*numdata
                        write(10,"(f12.4,f4.1)") linexall(1,itmp),spikey(itmp)
                    end do
                    close(10)
                end if
            end do
            deallocate(spikey)
		    CALL ENDGRF
			call color("WHITE")
            if (iexportlevel==1) write(*,"(a)") " Line data corresponding to various spike sets have been &
            exported to spike[index].txt in current folder"
        end if
        
		call disfin
		if (isavepic==1) write(*,*) "Graphical file has been saved to current folder with ""dislin"" prefix"
	end if

end do
end subroutine




!!----- Load transition data from a file to datax,str,FWHM in global memory
! ispectrum decides spectrum type, "filename" indicate the file to be loaded, numdata is the number of transitions in the file
! imode=0: Load data, imode=1: Only return numdata
! For electronic spectra, transition energies are loaded as eV
! istrtype is a global variable, it records type of strength. When loading multiple files, this variable avoids select type of strength multiple times
subroutine loadtransdata(imode,ispectrum,loadspecname,numdata)
use defvar
use util
implicit real*8 (a-h,o-z)
character(len=*) loadspecname
character ctest,ctest2,ctest3,c80tmp*80,c200tmp*200
integer ispectrum,imode
integer :: nrdfreq=0 ! >0 means pre-resonance raman, which loads external field frequency
real*8,allocatable :: rdfreq(:),tmparr(:)

if (allocated(datax)) deallocate(datax,str,FWHM)
open(10,file=loadspecname,status="old")

!Check if is sTDA output file
if (index(loadspecname,".dat")/=0) then
	if (imode==0) write(*,*) "Recognized as sTDA program output file"
	if (ispectrum==4.and.istrtype==0.and.imode==0) then
		write(*,*)
		write(*,*) "Read the rotatory strengths in which representation?"
		write(*,*) "1: Length representation     2: Velocity representation"
		write(*,*) "3: Mixed-form representation (recommended)"
		read(*,*) istrtype
	end if
	call loclabel(10,"DATXY")
	read(10,*)
	numdata=0
	do while(.true.)
		read(10,"(a)",iostat=ierror) c80tmp
		if (c80tmp==" ".or.ierror/=0) exit
		numdata=numdata+1
	end do
    if (imode==1) then !Have obtained number of data, return
        close(10)
        return
    end if
	allocate(datax(numdata),str(numdata),FWHM(numdata))
	if (ispectrum==3) FWHM=2D0/3D0
	if (ispectrum==4) FWHM=0.2D0
	call loclabel(10,"DATXY",ifound)
	read(10,*)
	do i=1,numdata
		!Note: The last four columns of tda.dat correspond to f_length, f_velocity, R_length, R_velocity
		read(10,*) inouse,datax(i),fl,fv,Rl,Rv
		if (ispectrum==3) then
			str(i)=fl
		else if (ispectrum==4) then
			if (istrtype==1) str(i)=Rl
			if (istrtype==2) str(i)=Rv
			if (istrtype==3) then
				str(i)=Rv
				if (fv/=0) str(i)=Rv*fl/fv
			end if
		end if
	end do
	close(10)
	return
end if

!Check if is Gaussian output file
call loclabel(10,"Gaussian, Inc",igauout,maxline=100)
rewind(10)
if (igauout==1) then
	if (imode==0) write(*,*) "Recognized as a Gaussian output file"

	!IR, Raman, VCD, ROA
	if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then
	
		!Raman and ROA, detect how many incident frequencies are there, and load which one
        if (imode==0) then
		    if (ispectrum==2.or.ispectrum==6) then
			    call loclabel(10,"NFrqRd=",ifound,0)
			    if (ifound==1) then
				    read(10,"(a)") c200tmp
				    itmp=index(c200tmp,"NFrqRd=")
				    read(c200tmp(itmp+7:),*) nrdfreq
				    if (nrdfreq>0) then
					    allocate(rdfreq(nrdfreq))
					    do itmp=0,nrdfreq,5
						    nleft=nrdfreq-itmp
						    read(10,"(a)") c200tmp
						    if (nleft>5) then
							    read(c200tmp(14:),*) rdfreq(itmp+1:itmp+5)
						    else
							    read(c200tmp(14:),*) rdfreq(itmp+1:nrdfreq)
						    end if
					    end do
					    if (ispectrum==2) write(*,*) "This is a pre-resonance Raman calculation, frequencies of incident light:"
					    if (ispectrum==6) write(*,*) "This is a ROA calculation, frequencies of incident light:"
					    do itmp=1,nrdfreq
						    if (rdfreq(itmp)==0D0) then
							    write(*,"(i5,':',f16.8,' a.u.')") itmp,rdfreq(itmp)
						    else
							    write(*,"(i5,':',f16.8,' a.u. (',f10.3,' nm )')") itmp,rdfreq(itmp),au2nm/rdfreq(itmp)
						    end if
					    end do
					    if (nrdfreq>1) then
						    write(*,*) "Load data for which frequency? Input the index, e.g. 3"
						    read(*,*) irdfreq
					    else
						    irdfreq=1
						    nrdfreq=1
					    end if
				    end if
			    end if
			    rewind(10)
		    end if
		    if (ispectrum==6) then
			    write(*,*)
			    write(*,*) "Load which type of data?" !For static incident light, no ROA is outputted
			    write(*,*) "1 Raman SCP(180)"
			    if (rdfreq(irdfreq)/=0D0) write(*,*) "2 ROA SCP(180), i.e. SCP backscattered ROA spectrum"
			    write(*,*) "3 Raman SCP(90)"
			    if (rdfreq(irdfreq)/=0D0) write(*,*) "4 ROA SCP(90)"
			    write(*,*) "5 Raman DCP(180)"
			    if (rdfreq(irdfreq)/=0D0) write(*,*) "6 ROA DCP(180)"
			    read(*,*) iROAtype
		    end if
        end if
		
		!Find how many frequencies in the file
		do while(.true.)
			call loclabel(10,"Frequencies -- ",ifound,0) !HPmodes is also compatible, because in this manner we locate to the traditional output section
			if (ifound==1) then
				i1=0
				i2=0
				i3=0
				backspace(10)
				backspace(10)
				read(10,*,iostat=ierror) i1,i2,i3
				if (ierror/=0) then
					read(10,*,iostat=ierror) i1,i2
					if (ierror/=0) then
						read(10,*,iostat=ierror) i1
					end if
				end if
				read(10,*)
				read(10,*)
				if (i1==0.or.i2==0.or.i3==0) exit
			else
				exit
			end if
		end do
		numdata=max(i1,i2,i3)
        if (imode==1) then !Have obtained number of data, return
            close(10)
            return
        end if
		rewind(10)

		!Load harmonic data
		allocate(datax(numdata),str(numdata),FWHM(numdata))
		FWHM=8D0
		ilackdata=numdata
		inow=1
		do while(.true.)
			if (ilackdata>3) then
				iread=3
			else
				iread=ilackdata
			end if
			call loclabel(10,"Frequencies -- ",ifound,0)
			read(10,"(16x)",advance="no")
			if (iread==1) read(10,*) datax(inow)
			if (iread==2) read(10,*) datax(inow),datax(inow+1)
			if (iread==3) read(10,*) datax(inow),datax(inow+1),datax(inow+2)
			if (ispectrum==1) then !IR
				call loclabel(10,"IR Inten    --",ifound,0)
                if (ifound==0) then
					write(*,*) "Error: Unable to find IR intensity information from this file!"
                    write(*,*) "Please double check your keyword in Gaussian input file"
                    write(*,*) "Press ENTER button to exit"
                    read(*,*)
                    stop
                end if
			else if (ispectrum==2) then !Raman
				if (nrdfreq==0) then !Normal Raman
					call loclabel(10,"Raman Activ --",ifound,0)
				else !Pre-resonance Raman
					write(c200tmp,"('RamAct Fr=',i2)") irdfreq
					call loclabel(10,trim(c200tmp),ifound,0)
				end if
                if (ifound==0) then
					write(*,*) "Error: Unable to find Raman activity information from this file!"
                    write(*,*) "Please double check your keyword in Gaussian input file"
                    write(*,*) "Press ENTER button to exit"
                    read(*,*)
                    stop
                end if
			else if (ispectrum==5) then !VCD
				call loclabel(10,"Rot. str.",ifound,0)
                if (ifound==0) then
					write(*,*) "Error: Unable to find rotatory strength information from this file!"
                    write(*,*) "Please double check your keyword in Gaussian input file"
                    write(*,*) "Press ENTER button to exit"
                    read(*,*)
                    stop
                end if
			else if (ispectrum==6) then !ROA
				if (iROAtype==1) write(c200tmp,"('Raman1 Fr=',i2)") irdfreq
				if (iROAtype==2) write(c200tmp,"('ROA1   Fr=',i2)") irdfreq
				if (iROAtype==3) write(c200tmp,"('Raman2 Fr=',i2)") irdfreq
				if (iROAtype==4) write(c200tmp,"('ROA2   Fr=',i2)") irdfreq
				if (iROAtype==5) write(c200tmp,"('Raman3 Fr=',i2)") irdfreq
				if (iROAtype==6) write(c200tmp,"('ROA3   Fr=',i2)") irdfreq
				call loclabel(10,trim(c200tmp),ifound,0)
                if (ifound==0) then
					if (iROAtype==1.or.iROAtype==3.or.iROAtype==5) then
						write(*,*) "Error: Unable to find Raman activity information information from this file!"
                    else
						write(*,*) "Error: Unable to find ROA information information from this file!"
                    end if
                    write(*,*) "Please double check your keyword in Gaussian input file"
                    write(*,*) "Press ENTER button to exit"
                    read(*,*)
                    stop
                end if
			end if
			read(10,"(16x)",advance="no")
			if (iread==1) read(10,*,iostat=ierror) str(inow)
			if (iread==2) read(10,*,iostat=ierror) str(inow),str(inow+1)
			if (iread==3) read(10,*,iostat=ierror) str(inow),str(inow+1),str(inow+2)
            if (ierror/=0) then
				backspace(10)
                read(10,"(a)") c80tmp
                write(*,*) "Error: Unable to load data successfully from this line:"
                write(*,*) trim(c80tmp)
                if (index(c80tmp,'*')/=0) write(*,"(a)") " Probably corresponding incident frequency is too close to an electron excitation energy, causing extremely large Raman activity"
                write(*,*) "Press ENTER button to exit"
                read(*,*)
                stop
            end if
			if (ilackdata<=3) exit
			ilackdata=ilackdata-3
			inow=inow+3
		end do
		
		!Load anharmonic data. Current Gaussian is unable to calculate anharmonic ROA data
		if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then  
			if (ispectrum==1) call loclabel(10,"Anharmonic Infrared Spectroscopy",ifound,0)
			if (ispectrum==2) call loclabel(10,"Anharmonic Raman Spectroscopy",ifound,0)
			if (ispectrum==5) call loclabel(10,"Anharmonic VCD Spectroscopy",ifound,0)
			if (ifound==1) then
				write(*,"(a)") " Found anharmonic data, if load them instead of the harmonic one? (y/n)"
				read(*,*) ctest
				if (ctest=='y'.or.ctest=='Y') then
					write(*,"(a)") " If also load overtone band frequencies? (y/n)"
					read(*,*) ctest2
					write(*,"(a)") " If also load combination band frequencies? (y/n)"
					read(*,*) ctest3
                    
                    !Detect total number of fund. over. comb. bands.
                    !This cannot be easily infer from number of harmonic frequency for non-linear system, i.e. over.=numdata, comb.=nummode*(nummode-1)/2
                    !However for linear system, the case is complicated, so we directly derect the number of outputted lines
                    numfund=numdata !Number of Anharm. is always identical to harm
                    write(*,"(' Number of fundamental frequencies:',i6)") numfund
                    numover=0
                    if (ctest2=='y'.or.ctest2=='Y') then
                        call loclabel(10,"Overtones",ifound,0)
                        read(10,*);read(10,*);read(10,*)
                        do while(.true.)
                            read(10,"(a)") c80tmp
                            if (c80tmp==" ") exit
                            numover=numover+1
                        end do
                        write(*,"(' Number of Overtone frequencies:   ',i6)") numover
                    end if
                    numcomb=0
                    if (ctest3=='y'.or.ctest3=='Y') then
                        call loclabel(10,"Combination Bands",ifound,0)
                        read(10,*);read(10,*);read(10,*)
                        do while(.true.)
                            read(10,"(a)") c80tmp
                            if (c80tmp==" ") exit
                            numcomb=numcomb+1
                        end do
                        write(*,"(' Number of combination frequencies:',i6)") numcomb
                    end if
                    numdata=numfund+numover+numcomb
					
					deallocate(datax,str,FWHM)
					allocate(datax(numdata),str(numdata),FWHM(numdata))
			        if (ispectrum==1) call loclabel(10,"Anharmonic Infrared Spectroscopy")
			        if (ispectrum==2) call loclabel(10,"Anharmonic Raman Spectroscopy")
			        if (ispectrum==5) call loclabel(10,"Anharmonic VCD Spectroscopy")
					FWHM=8D0
					idata=0
					call loclabel(10,"Fundamental Bands",ifound,0)
					read(10,*);read(10,*);read(10,*)
					do itmp=1,numfund
						idata=idata+1
                        read(10,"(22x)",advance='no')
						read(10,*) harmfreq,datax(idata),rnouse,str(idata)
						if (ispectrum==2) str(idata)=0.059320323D0*harmfreq*str(idata) !The conversion coefficient can be found in output file
					end do
					if (ctest2=='y'.or.ctest2=='Y') then
						call loclabel(10,"Overtones",ifound,0)
						read(10,*);read(10,*);read(10,*)
						do itmp=1,numover
							idata=idata+1
                            read(10,"(22x)",advance='no')
							read(10,*) harmfreq,datax(idata),str(idata)
							if (ispectrum==2) str(idata)=0.059320323D0*harmfreq*str(idata) !The conversion coefficient can be found in output file
						end do
					end if
					if (ctest3=='y'.or.ctest3=='Y') then
						call loclabel(10,"Combination Bands",ifound,0)
						read(10,*);read(10,*);read(10,*)
                        !In order to compatible with different version, find the position of the last character of the two labels
                        read(10,"(a)") c80tmp
                        iskip=index(c80tmp,')',back=.true.)
                        backspace(10)
						do itmp=1,numcomb
							idata=idata+1
                            read(10,"(a)") c80tmp
                            read(c80tmp(iskip+1:),*) harmfreq,datax(idata),str(idata)
							if (ispectrum==2) str(idata)=0.059320323D0*harmfreq*str(idata) !The conversion coefficient can be found in output file
						end do
					end if
				end if
			end if
		end if
	
	!UV-Vis, ECD
	else if (ispectrum==3.or.ispectrum==4) then
		!Because this may be an excited state optimization task, we need to determine how many steps are there, so that we can locate to the last output
        !Besides, for EOM-CCSD task, the EOM-CCSD energies and transition moments are outputted after CIS part, this also needs to find the last time of output
		numopt=0
		do while(.true.)
			call loclabel(10,"Excitation energies and oscillator strengths",ifound,0)
			numopt=numopt+ifound
			if (ifound==0) exit
			read(10,*)
		end do
		!Check how many states
		rewind(10)
		do i=1,numopt !Locate to the last time of output
			call loclabel(10,"Excitation energies and oscillator strengths",ifound,0)
			read(10,*)
		end do
		numdata=0
		do while(.true.)
			call loclabel(10," Excited State",ifound,0)
			if (ifound==0) exit
			read(10,*)
			numdata=numdata+1
		end do
        if (imode==1) then !Have obtained number of data, return
            close(10)
            return
        end if
		allocate(datax(numdata),str(numdata),FWHM(numdata))
		FWHM=2/3D0
		!Locate to the last time of output
		rewind(10)
		do i=1,numopt
			call loclabel(10,"Excitation energies and oscillator strengths",ifound,0)
			read(10,*)
		end do
		do i=1,numdata !Gaussian output is too flexible to use fixed format to read in
			call loclabel(10," Excited State",ifound,0)
			do while(.true.)
				read(10,"(a)",advance="no") ctest
				if (ctest=='-') exit
			end do
			read(10,"(5x)",advance="no")
			read(10,*) datax(i) !Read excitation energy (eV)
			backspace(10)
			do while(.true.)
				read(10,"(a)",advance="no") ctest
				if (ctest=='=') exit
			end do
			read(10,*) str(i) !Read oscillator strength
		end do
		if (ispectrum==4) then !Read ECD rotatory strength
			if (istrtype==0) then
				write(*,*) "Read the rotatory strengths in which representation?"
				write(*,*) "1: Length representation     2: Velocity representation (Recommended)"
				read(*,*) istrtype
			end if
            rewind(10)
			if (istrtype==1) then
		        do i=1,numopt !Locate to the last time of output
				    call loclabel(10,"R(length)",ifound,0)
				    read(10,*)
                end do
				do i=1,numdata
					read(10,*) inouse,rnouse,rnouse,rnouse,str(i)
				end do
			else
		        do i=1,numopt !Locate to the last time of output
				    call loclabel(10,"R(velocity)",ifound,0)
				    read(10,*)
                end do
				do i=1,numdata
					read(10,*) inouse,rnouse,rnouse,rnouse,str(i)
					do while(.true.)
						read(10,"(a)") c80tmp
						!Sometimes Gaussian output some additional info.
						if (index(c80tmp,"Total R(velocity) tensor")/=0.or.index(c80tmp,"R(velocity) tensor in inp. orien.")/=0) then
							do itmp=1,4
								read(10,*)
							end do
						else
							backspace(10)
							exit
						end if
					end do
				end do
			end if
		end if
	end if
	close(10)
	return
end if	

!Check if is ORCA output file
call loclabel(10,"O   R   C   A",iORCAout,maxline=100)
if (iORCAout==1) then
	if (imode==0) write(*,*) "Recognized as an ORCA output file"
	isTDA=0
	if (ispectrum==3.or.ispectrum==4) call loclabel(10,"ORCA sTD",isTDA) !When plotting UV-Vis or ECD, check if this is a sTDA or sTD-DFT calculation
	if (isTDA==0) then !Regular calculation
		if (ispectrum==1.or.ispectrum==2) then !IR, Raman
			if (ispectrum==1) call loclabel(10,"IR SPECTRUM")
			if (ispectrum==2) call loclabel(10,"RAMAN SPECTRUM")
			call loclabel(10," Mode    freq (cm**-1)",ifound,0)
			if (ifound==1) then !ORCA 4
				read(10,*);read(10,*)
			else !ORCA 5
				call loclabel(10," Mode   freq  ",ifound,1)
				read(10,*);read(10,*);read(10,*)
			end if
			numdata=0
			do while(.true.)
				read(10,"(a)") c80tmp
				if (c80tmp==" ") exit
				numdata=numdata+1
			end do
            if (imode==1) then !Have obtained number of data, return
                close(10)
                return
            end if
			allocate(datax(numdata),str(numdata),FWHM(numdata))
			if (ispectrum==1) call loclabel(10,"IR SPECTRUM",ifound,1)
			if (ispectrum==2) call loclabel(10,"RAMAN SPECTRUM",ifound,1)
			call loclabel(10,"Mode    freq (cm**-1)",ifound,0)
			if (ifound==1) then !ORCA 4
				read(10,*);read(10,*)
			else !ORCA 5
				call loclabel(10," Mode   freq  ",ifound,1)
				read(10,*);read(10,*);read(10,*)
			end if
			!The IR data under T**2 is in KM/mole
			do i=1,numdata
				read(10,*) c80tmp,datax(i),str(i)
			end do
			FWHM=8D0
		else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
			call loclabel(10,"Number of roots to be determined",ifound)
            if (ifound==1) then !TDDFT/TDA-DFT
				read(10,"(50x,i7)") numdata
				if (imode==1) then !Have obtained number of data, return
					close(10)
					return
				end if
				allocate(datax(numdata),str(numdata),FWHM(numdata))
				if (ispectrum==3) call loclabel(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",ifound,0)
				if (ispectrum==4) call loclabel(10,"CD SPECTRUM",ifound,0)
				call skiplines(10,5)
                datax=0;str=0
				do i=1,numdata
					read(10,*,iostat=ierror) itmp,t1,rnouse,t2
                    if (ierror/=0) exit !For SF-TDDFT, number of states recorded in this field is less than nstates by 1, because one of SF-TDDFT states is viewed as ground state
                    datax(itmp)=t1
                    str(itmp)=t2
				end do
				call loclabel(10,"SOC CORRECTED ABSORPTION",ifound,0)
				if (ifound==1) then
					write(*,"(a)") " Spin-orbit coupling corrected spectra information was found, &
					would you like to plot this kind of spectrum instead of the one without correction? (y/n)"
					read(*,*) ctest
					if (ctest=='y'.or.ctest=='Y') then
						numdata=4*numdata !if root=n, then there will be n singlet states and 3n triplet sublevels
						deallocate(datax,str,FWHM)
						allocate(datax(numdata),str(numdata),FWHM(numdata))
						if (ispectrum==3) continue
						if (ispectrum==4) call loclabel(10,"CD SPECTRUM",ifound,0)
						read(10,*)
						read(10,*)
						read(10,*)
						read(10,*)
						read(10,*)
						do i=1,numdata
							read(10,*) tmp1,tmp2,tmp3,tmp4,tmp5 !Since 5.0, the first column is always 0, strange!
							if (tmp1==0) then !ORCA >=5.0
								datax(i)=tmp3
								str(i)=tmp5
							else !ORCA 4.x
								datax(i)=tmp2
								str(i)=tmp4
							end if
						end do
					end if
				end if
            else !Should be (DLPNO-)(ST)EOM-CCSD or others
				!Count how many data are there
				call loclabel(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS")
                call skiplines(10,5)
                numdata=0
                do while(.true.)
					read(10,"(a)") c80tmp
                    if (c80tmp==" ") exit
					numdata=numdata+1
                end do
				allocate(datax(numdata),str(numdata),FWHM(numdata))
				if (ispectrum==3) call loclabel(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS")
				if (ispectrum==4) call loclabel(10,"CD SPECTRUM")
				call skiplines(10,5)
				do i=1,numdata
					read(10,*) rnouse,datax(i),tmpval,str(i)
                    if (tmpval>datax(i)) then !I noticed ORCA 5.0.1 has a bug, the values in cm-1 and nm are inversed
						ttt=tmpval
                        datax(i)=tmpval
                        tmpval=ttt
                    end if
				end do
            end if
			FWHM=2D0/3D0
			datax=datax/8065.5447D0 !Convert from cm-1 to eV
		end if
	else if (isTDA==1) then !sTDA or sTD-DFT calculation
        if (imode==0) then
		    write(*,*) "This is a sTDA or sTD-DFT calculation"
		    if (ispectrum==4.and.istrtype==0) then
			    write(*,*)
			    write(*,*) "Read the rotatory strengths in which representation?"
			    write(*,*) "1: Length representation     2: Velocity representation"
			    write(*,*) "3: Mixed-form representation (recommended)"
			    read(*,*) istrtype
		    end if
        end if
		call loclabel(10,"roots found,",ifound,0)
		read(10,*) numdata
        if (imode==1) then !Have obtained number of data, return
            close(10)
            return
        end if
		allocate(datax(numdata),str(numdata),FWHM(numdata))
		call loclabel(10,"state   eV        nm        fL",ifound,0)
		read(10,*)
		do i=1,numdata
			read(10,*) inouse,datax(i),rnouse,fl,fv,Rl,Rv
			if (ispectrum==3) then
				str(i)=fl
			else if (ispectrum==4) then
				if (istrtype==1) str(i)=Rl
				if (istrtype==2) str(i)=Rv
				if (istrtype==3) then
					str(i)=Rv
					if (fv/=0) str(i)=Rv*fl/fv
				end if
			end if
		end do
		if (ispectrum==3) FWHM=2D0/3D0
		if (ispectrum==4) FWHM=0.2D0
	end if
    close(10)
    return
end if

call loclabel(10,"$vibrational spectrum",ixtb,maxline=100)
if (ixtb==1) then
    if (imode==0) then
	    write(*,*) "Recognized as a vibspectrum file generated by Grimme's xtb program"
	    write(*,*)
    end if
    rewind(10)
    numdata=0
    do while(.true.)
		read(10,"(a)",iostat=ierror) c80tmp
        if (ierror/=0) exit
        if (index(c80tmp,'#')==0.and.index(c80tmp,'$')==0.and.index(c80tmp,' - ')==0) numdata=numdata+1
    end do
	allocate(FWHM(numdata),datax(numdata),str(numdata))
	FWHM=8D0
    rewind(10)
    itmp=0
    do while(.true.)
		read(10,"(a)",iostat=ierror) c80tmp
        if (ierror/=0.or.c80tmp==" ") exit
        if (index(c80tmp,'#')==0.and.index(c80tmp,'$')==0.and.index(c80tmp,' - ')==0) then
			itmp=itmp+1
			read(c80tmp,*) inouse,c200tmp,datax(itmp),str(itmp)
        end if
    end do
	close(10)
	return
end if

!Check if is CP2K output file
call loclabel(10,"CP2K|",iCP2K,maxline=300)
if (iCP2K==1) then
    if (imode==0) then
	    write(*,*) "Recognized as a CP2K output file"
	    write(*,*)
    end if
    rewind(10)
    
    if (ispectrum==3) then !UV-Vis
		call loclabelfinal(10,"number   energy (eV)",ifound)
        if (ifound==0) then
			write(*,*) "Error: Unable to find electronic excitation information!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            return
        end if
        read(10,*)
        read(10,*)
		numdata=0
		do while(.true.)
			read(10,"(a)") c80tmp
            if (c80tmp==" ") exit
			numdata=numdata+1
		end do
        if (imode==1) then !Have obtained number of data, return
            close(10)
            return
        end if
		allocate(datax(numdata),str(numdata),FWHM(numdata))
		FWHM=2/3D0
		do i=1,numdata+1
			backspace(10)
        end do
		do i=1,numdata !Read excitation energy (eV) and oscillator strength
			read(10,*) c80tmp,c80tmp,datax(i),c80tmp,c80tmp,c80tmp,str(i)
		end do
    
    else if (ispectrum==1) then !IR
		!Find how many frequencies in the file
		i1=0
		i2=0
		i3=0
		do while(.true.)
			call loclabel(10,"VIB|Frequency (cm^-1)",ifound,0)
			if (ifound==1) then
				i1=0;i2=0;i3=0
				backspace(10)
				read(10,"(a)") c200tmp
				read(c200tmp(10:),*,iostat=ierror) i1,i2,i3
				if (ierror/=0) then
					read(c200tmp(10:),*,iostat=ierror) i1,i2
					if (ierror/=0) then
						read(c200tmp(10:),*,iostat=ierror) i1
					end if
				end if
				read(10,*);read(10,*)
				if (i1==0.or.i2==0.or.i3==0) exit
			else
				exit
			end if
		end do
		rewind(10)
		numdata=max(i1,i2,i3)

		allocate(datax(numdata),str(numdata),FWHM(numdata))
		FWHM=8D0
		ilackdata=numdata
		inow=1
		do while(.true.)
			if (ilackdata>3) then
				iread=3
			else
				iread=ilackdata
			end if
			call loclabel(10,"VIB|Frequency (cm^-1)",ifound,0)
			read(10,"(22x)",advance="no")
			if (iread==1) read(10,*) datax(inow)
			if (iread==2) read(10,*) datax(inow),datax(inow+1)
			if (iread==3) read(10,*) datax(inow),datax(inow+1),datax(inow+2)
			read(10,"(22x)",advance="no")
			if (iread==1) read(10,*) str(inow)
			if (iread==2) read(10,*) str(inow),str(inow+1)
			if (iread==3) read(10,*) str(inow),str(inow+1),str(inow+2)
			if (ilackdata<=3) exit
			ilackdata=ilackdata-3
			inow=inow+3
		end do
    end if
	close(10)
    return
end if

!Plain text file
if (imode==0) write(*,*) "Recognized as a plain text file"
rewind(10)
read(10,*) numdata,inptype
if (imode==1) then !Have obtained number of data, return
    close(10)
    return
end if
allocate(datax(numdata),str(numdata),FWHM(numdata))
if (inptype==1) then !Only x-position and strengths
	if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) FWHM=8D0
	if (ispectrum==3) FWHM=2D0/3D0
	if (ispectrum==4) FWHM=0.2D0
	do i=1,numdata
		read(10,*) datax(i),str(i)
	end do
else if (inptype==2) then !also with FWHM
	do i=1,numdata
		read(10,*) datax(i),str(i),FWHM(i)
	end do
end if
close(10)

end subroutine






!!--------------- Define a module for exchanging data between subroutine NMRplot and subroutine loadNMR
module NMRmod
real*8,allocatable :: atmshd(:) !Isotropic shielding (ppm) loaded temporarily
real*8,allocatable :: atmshdall(:,:) !atmshd of all systems. (ncenter,nsystem)
end module

    
!!--------------- Plotting NMR spectrum
subroutine NMRplot
use defvar
use NMRmod
use dislin_d
use plot
use util
implicit real*8 (a-h,o-z)
!Plotting settings
character(len=3) :: elemplot="all" !The element to be plotted. "all" means all elements
integer :: igetshift=0 !The way of determining chemical shift. =0: None (only show shielding), =1: Use specific reference value, =2: Scaling method
real*8 :: NMRref=0 !Reference value of NMR in ppm, used for igetshift=1
real*8 :: NMRslope=1,NMRinter=0 !Slope and intercept of scaling method, used for igetshift=2
real*8 :: FWHM_NMR=0.2D0 !FWHM for broadening NMR data to curve. 0.2 for heavy atom, 0.02 for H.
  !Note that for reproducing experimental NMR, FWHM=0.01 is seemingly best choice. However, in order to make broadening effect clearly visible, we use larger FWHM
real*8 :: degentol=0.05D0 !Degeneracy tolerance in ppm
integer :: ishowline=1,ishowcurve=1,ishowgrid=1,ishowweighted=1
integer :: thk_curve=3,thk_weicurve=3,thk_line=3,thk_weiline=3,thk_legend=2,thk_axis=1,thk_grid=1 !thickness
integer :: nticksX=2,nticksY=2 !Number of ticks in X-axis, in signal strength axis
integer,allocatable :: line_clr(:),curve_clr(:) !Color list for discrete lines and curves. 0 corresponds to weighted data, 1/2/3... corresponds to 1/2/3th system
integer :: legposx=2200,legposy=155 !Legend positions in X and Y
integer :: ishowlegend=1 !Show legends when nsystem>1
integer :: ishowlabel=1,ilabelele=0,labelsize=45,legendsize=40,ilabelclr=3,labelshiftX=-14,labelshiftY=0,ionlyatmlabwei=1,ialllabtop=0 !Parameters of showing atomic labels
integer :: ishowdataright=1 !If show number on the Y-axis corresponding to NMR signal strength
!System information
real*8,allocatable :: weight(:) !Weight of various system for plotting mixed spectrum
character(len=80),allocatable :: mollegend(:) !Legends for multiple systems
real*8,allocatable :: atmshdall_org(:,:) !atmshdall will be converted to chemical shift before plotting, atmshdall_org backs up original values
!Data for final plotting. The "term" denotes the spike with degeneracy of 1 or more
integer,allocatable :: shdnum(:) !shdnum(r) is number of terms of system r
real*8,allocatable :: shdval(:,:) !shdval(j,r) is shielding value of j term for system r
integer,allocatable :: shdnatm(:,:) !shdnatm(j,r) is degeneracy of j term for system r
real*8,allocatable :: shdeffnatm(:,:) !shdeffnatm(j,r) is effective degeneracy of j term for system r (i.e. "strength" is taken into account)
integer,allocatable :: shdatm(:,:,:) !shdatm(1/2/3...,j,r) is index of 1/2/3...th atom in j term of system r
real*8,allocatable :: linexall(:,:),lineyall(:,:) !Array used to draw discrete lines for all systems. The 1st index corresponds to system index
real*8,allocatable :: curveyall(:,:) !Curves of all systems. The 1st index corresponds to system index. Global variable "curvex" is shared by all curves
real*8,allocatable :: atmstr(:) !(i) is strength of atom i, default to 1
!Weighted spectrum
real*8,allocatable :: atmshdwei(:),atmshdwei_org(:) !Weighting mixed isotropic shielding (ppm) and the array for backing up
integer shdnumwei !Number of terms of weighted spectrum
real*8,allocatable :: shdvalwei(:) !Shielding value of weighted spectrum
integer,allocatable :: shdnatmwei(:) !(j) is number of atoms of j term of weighted spectrum
real*8,allocatable :: shdeffnatmwei(:) !(j) is effective number of atoms of j term of weighted spectrum (i.e. "strength" is taken into account)
integer,allocatable :: shdatmwei(:,:) !(1/2/3...,j,r) is index of 1/2/3...th atom in j term of system r
real*8,allocatable :: linexwei(:),lineywei(:),curveywei(:)
!Other variables
character selectyn,c80tmp*80,c200tmp*200,c200tmp2*200,c2000tmp*2000,clegend*2000 !Buffer for showing legends
integer,allocatable :: tmparrint(:)


!Sequence and number of atoms in multiple systems must be the same!!!
if (index(filename,"multiple.txt")/=0) then !Load multiple systems
    nsystem=0
	open(11,file=filename,status="old") !Note that fileid=10 will be used by loadNMR
	!Count total number of entries
	do while(.true.)
		read(11,*,iostat=ierror) c200tmp
		if (ierror/=0.or.c200tmp==" ") exit
		nsystem=nsystem+1
	end do
	write(*,"(' There are',i4,' systems')") nsystem
	allocate(weight(nsystem),mollegend(nsystem))
	mollegend=" "
	rewind(11)
	do i=1,nsystem
		read(11,"(a)") c200tmp
		c200tmp2=c200tmp
		read(c200tmp2,*,iostat=ierror) c200tmp,weight(i)
		if (ierror==0) then !The second field is weight
			write(mollegend(i),"(i3,' (',f5.1,'%)')") i,weight(i)*100
		else !The second field is legend rather than weight value
			ispc=index(c200tmp," ")
			read(c200tmp2(:ispc-1),*) c200tmp
			read(c200tmp2(ispc+1:),"(a)") mollegend(i)
			weight(i)=1
			!If the first letter of the legend is $, it will be skipped
			if (mollegend(i)(1:1)=='$') mollegend(i)=mollegend(i)(2:)
		end if
		inquire(file=c200tmp,exist=alive)
		if (.not.alive) then
			write(*,"(' Error: Cannot find ',a)") trim(c200tmp)
			if (index(c200tmp,'/')/=0) then
				write(*,*) "Reminder: Since the file path contains / symbol, you should add "" at the two ends of the path, so that the file can be properly loaded"
			end if
			write(*,*) "Press ENTER button to exit program"
			read(*,*)
			stop
		end if
		if (weight(i)==1) then
			write(*,"(' Loading ',a,'    Legend: ',a)") trim(c200tmp),trim(mollegend(i))
		else
			write(*,"(' Loading ',a,'    Weight:',f7.4)") trim(c200tmp),weight(i)
		end if
        
        call loadNMR(trim(c200tmp)) !Load atomic shielding tensor to atmshd
        if (.not.allocated(atmshdall)) allocate(atmshdall(ncenter,nsystem))
        atmshdall(:,i)=atmshd(:)
        deallocate(atmshd)
	end do
	close(11)
	if (all(weight==1)) then !When all weights are unity, then no weighted spectrum will be plotted, but simply plotting all systems together
		ishowweighted=0
	end if
else !Only load one system
    nsystem=1
    call loadNMR(filename) !Load atomic shielding tensor to atmshd
    allocate(atmshdall(ncenter,1),weight(ncenter))
    atmshdall(:,1)=atmshd(:)
    weight=1
    ishowweighted=0
end if

!Allocate arrays used in plotting. ncenter is upper limit of (non)degenerate terms
allocate(shdnum(nsystem))
allocate(shdval(ncenter,nsystem))
allocate(shdnatm(ncenter,nsystem),shdeffnatm(ncenter,nsystem))
allocate(shdatm(20,ncenter,nsystem)) !Degeneracy is assumed to be at most 20
allocate(atmshdall_org(ncenter,nsystem))
allocate(atmstr(ncenter))
atmstr=1
!Generate weighted data. ncenter is upper limit
if (nsystem>1) then
    allocate(shdvalwei(ncenter),shdnatmwei(ncenter),shdeffnatmwei(ncenter),shdatmwei(20,ncenter))
    !If user defined weights, then calculated mixed shielding value
    if (any(weight/=1)) then
        allocate(atmshdwei(ncenter),atmshdwei_org(ncenter))
        sumwei=sum(weight)
        if (abs(sumwei-1)>0.001D0) then
            write(*,"(/,a,f7.2,a)") " Warning: The sum of all conformation weights is",sumwei*100,&
            " %, which deviates from 100% evidently, this may make the resulting map meaningless. Do you want to normalize the weights? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y'.or.selectyn=='Y') weight(:)=weight(:)/sumwei
        end if
        do iatm=1,ncenter
            atmshdwei(iatm)=sum(weight(:)*atmshdall(iatm,:))
        end do
    end if
end if

!Initialize colors
allocate(line_clr(0:nsystem),curve_clr(0:nsystem))
if (nsystem==1) then
    line_clr(1)=5 !Black
    curve_clr(1)=1 !Red
else
    !Black for weighted spectrum
    line_clr(0)=5
    curve_clr(0)=5
    idx=1
    do imol=1,nsystem !Assign colors for every system. black (corresponding to goodcolor(1)) is skipped
        idx=idx+1
        if (idx>ngoodcolor) idx=2
        iclr=goodcolor(idx)
        line_clr(imol)=iclr
        curve_clr(imol)=iclr
    end do
end if

!curveytmp is a temporary array used to record curve during generating curve of each system
if (allocated(curvex)) deallocate(curvex) !Global array
if (allocated(curveytmp)) deallocate(curveytmp) !Global array
allocate(curvex(num1Dpoints),curvey(num1Dpoints),curveytmp(num1Dpoints),curveyall(nsystem,num1Dpoints))
allocate(linexall(nsystem,3*ncenter),lineyall(nsystem,3*ncenter))
if (ishowweighted==1) allocate(linexwei(3*ncenter),lineywei(3*ncenter),curveywei(num1Dpoints))

!Set default element to be shown
if (any(a%name=="C")) then
    elemplot="C"
else if (any(a%name=="H")) then
    elemplot="H"
else
    elemplot=a(1)%name
end if

iusersetX=0
iusersetY1=0
iusersetY2=0

do while(.true.)
    write(*,"(/,a)") " Hint: You can input ""s"" to save current plotting settings to a file, or input ""l"" to load settings from a file"
	write(*,*)
    write(*,*) "        ==================== Plotting NMR spectrum ===================="
	write(*,*) "-10 Return to main menu"
	write(*,*) "-3 Set format of saving graphical file, current: "//graphformat
	write(*,*) "-1 Show NMR data on screen        -2 Export NMR data to a plain text file"
	write(*,*) "0 Plot NMR spectrum now!"
	write(*,*) "1 Save graphical file of NMR spectrum in current folder"
	write(*,*) "2 Export X-Y data set of spikes and curves to a plain text file"
    if (iusersetX==0) write(*,*) "3 Set lower and upper limits of X-axis, current: Auto"
	if (iusersetX==1) write(*,"(a,f8.1,a,f8.1,' ppm')") " 3 Set lower and upper limits of X-axis, current:",xlow," to",xhigh
	if (iusersetY1==0) write(*,*) "4 Set range of Y-axis corresponding to degeneracy, current: Auto"
	if (iusersetY1==1) write(*,"(' 4 Set Y-axis of degeneracy, current: low:',f8.2,' up:',f8.2,' step:',f6.2)") orgy1,endy1,stepy1
	if (iusersetY2==0) write(*,*) "5 Set range of Y-axis corresponding to signal strength, current: Auto"
	if (iusersetY2==1) write(*,"(' 5 Set Y-axis of signal strength, current: low:',f8.2,' up:',f8.2,' step:',f6.2)") orgy2,endy2,stepy2
    write(*,*) "6 Choose the element considered in plotting, current: "//elemplot
    if (igetshift==0) write(*,"(a)") " 7 Set how to determine chemical shifts, current: None"
    if (igetshift==1) write(*,"(a,f9.3,' ppm')") " 7 Set how to determine chemical shifts, current: Using reference value of",NMRref
    if (igetshift==2) write(*,"(a,f9.4,a,f9.4)") " 7 Set how to determine chemical shifts, current: Scaling method, slope:",NMRslope,", intercept:",NMRinter
    write(*,"(a,f6.2,' ppm')") " 8 Set full width at half maximum (FWHM) for broadening, current:",FWHM_NMR
    write(*,"(a,f6.3,' ppm')") " 9 Set tolerance for determining NMR degeneracy, current:",degentol
    write(*,*) "10 Average shielding values of specific atoms"
    write(*,*) "11 Set strength of specific atoms"
	if (ishowline==1) write(*,*) "12 Toggle showing spikes, current: ON"
	if (ishowline==0) write(*,*) "12 Toggle showing spikes, current: OFF"
	if (ishowcurve==1) write(*,*) "13 Toggle showing curves, current: ON"
	if (ishowcurve==0) write(*,*) "13 Toggle showing curves, current: OFF"
    write(*,*) "14 Set colors of curves and spikes"
	write(*,*) "15 Set thickness of curves/lines/texts/axes/grid"
	write(*,*) "16 Change setting of labelling atoms"
	if (any(weight/=1)) write(*,*) "17 Set status of showing weighted spectrum and that of individual systems"
    write(*,*) "18 Other plotting settings"
    
    read(*,*) c80tmp
    
    if (index(c80tmp,'s')/=0) then
        write(*,"(a)") " Input file path for saving plotting settings, e.g. C:\Bang_Dream\RAS.dat"
        write(*,"(a)") " Note: If you press ENTER button directly, status will be saved to NMR.dat in current folder"
        read(*,"(a)") c200tmp
        if (c200tmp==" ") c200tmp="NMR.dat"
        open(10,file=c200tmp,status="replace")
        if (iusersetX==1) then
            write(10,*) "iusersetX"
            write(10,*) xlow
            write(10,*) xhigh
            write(10,*) stepx
        end if
        if (iusersetY1==1) then
            write(10,*) "iusersetY1"
            write(10,*) orgy1
            write(10,*) endy1
            write(10,*) stepy1
        end if
        if (iusersetY2==1) then
            write(10,*) "iusersetY2"
            write(10,*) orgy2
            write(10,*) endy2
            write(10,*) stepy2
        end if
        write(10,*) "Other"
        write(10,*) elemplot
        write(10,*) igetshift
        write(10,*) NMRref
        write(10,*) NMRslope
        write(10,*) NMRinter
        write(10,*) FWHM_NMR
        write(10,*) degentol
        write(10,*) ishowline
        write(10,*) ishowcurve
        write(10,*) ishowgrid
        write(10,*) ishowweighted
        write(10,*) thk_curve
        write(10,*) thk_weicurve
        write(10,*) thk_line
        write(10,*) thk_weiline
        write(10,*) thk_legend
        write(10,*) thk_axis
        write(10,*) thk_grid
        write(10,*) nticksX
        write(10,*) nticksY
        write(10,*) ishowlabel
        write(10,*) ilabelele
        write(10,*) labelsize
        write(10,*) ilabelclr
        write(10,*) labelshiftX
        write(10,*) labelshiftY
        write(10,*) ionlyatmlabwei
        write(10,*) ialllabtop
        write(10,*) legposx
        write(10,*) legposy
        write(10,*) ishowlegend
        write(10,*) graphformat
        write(10,*) "color",nsystem
        do iclr=0,nsystem
            write(10,*) line_clr(iclr),curve_clr(iclr)
        end do
        close(10)
        write(*,*) "Done!"
        cycle
    else if (index(c80tmp,'l')/=0)then
        write(*,"(a)") " Input file path to load plotting settings from it, e.g. C:\Bang_Dream\RAS.dat"
        write(*,"(a)") " Note: If you press ENTER button directly, status will be load from NMR.dat in current folder"
        read(*,"(a)") c200tmp
        if (c200tmp==" ") c200tmp="NMR.dat"
	    inquire(file=c200tmp,exist=alive)
	    if (.not.alive) then
	        write(*,*) "Error: Cannot find the file! Press ENTER button to return"
            read(*,*)
            cycle
        end if
        open(10,file=c200tmp,status="old")
        call loclabel(10,"iusersetX",iusersetX)
        if (iusersetX==1) then
            read(10,*)
            read(10,*) xlow
            read(10,*) xhigh
            read(10,*) stepx
        end if
        call loclabel(10,"iusersetY1",iusersetY1)
        if (iusersetY1==1) then
            read(10,*)
            read(10,*) orgy1
            read(10,*) endy1
            read(10,*) stepy1
        end if
        call loclabel(10,"iusersetY2",iusersetY2)
        if (iusersetY2==1) then
            read(10,*)
            read(10,*) orgy2
            read(10,*) endy2
            read(10,*) stepy2
        end if
        call loclabel(10,"Other")
        read(10,*)
        read(10,*) elemplot
        read(10,*) igetshift
        read(10,*) NMRref
        read(10,*) NMRslope
        read(10,*) NMRinter
        read(10,*) FWHM_NMR
        read(10,*) degentol
        read(10,*) ishowline
        read(10,*) ishowcurve
        read(10,*) ishowgrid
        read(10,*) ishowweighted
        read(10,*) thk_curve
        read(10,*) thk_weicurve
        read(10,*) thk_line
        read(10,*) thk_weiline
        read(10,*) thk_legend
        read(10,*) thk_axis
        read(10,*) thk_grid
        read(10,*) nticksX
        read(10,*) nticksY
        read(10,*) ishowlabel
        read(10,*) ilabelele
        read(10,*) labelsize
        read(10,*) ilabelclr
        read(10,*) labelshiftX
        read(10,*) labelshiftY
        read(10,*) ionlyatmlabwei
        read(10,*) ialllabtop
        read(10,*) legposx
        read(10,*) legposy
        read(10,*) ishowlegend
        read(10,*) graphformat
        call loclabel(10,"color")
        read(10,*) c80tmp,nsystemtmp
        do iclr=0,min(nsystem,nsystemtmp)
            read(10,*) line_clr(iclr),curve_clr(iclr)
        end do
        close(10)
        if (all(weight==1).and.ishowweighted/=0) ishowweighted=0
        write(*,*) "Loading finished!"
        cycle
    else
        read(c80tmp,*) isel
    end if
    
    if (isel==-10) then
		if (allocated(atmshd)) deallocate(atmshd)
		if (allocated(atmshdall)) deallocate(atmshdall)
        return
        
    else if (isel==-3) then
        call setgraphformat
        
    else if (isel==-1.or.isel==-2) then
        if (isel==-1) then
            ides=6
        else if (isel==-2) then
            ides=10
            open(10,file="NMRdata.txt",status="replace")
        end if
        write(ides,*) "Note: The data shown below are in ppm"
        do isystem=1,nsystem
            if (nsystem>1) then
		        if (weight(isystem)==1) then
                    write(ides,"(/,' System',i5,'  (',a,'):')") isystem,trim(mollegend(isystem))
		        else
                    write(ides,"(/,' System',i5,'  (weight=',f7.4,'):')") isystem,weight(isystem)
		        end if
            else
                write(ides,*)
            end if
            if (igetshift==0) then
                write(ides,*) "   Atom       Shielding(iso)"
            else
                write(ides,*) "   Atom       Shielding(iso)      Chemical Shift"
            end if
            do iatm=1,ncenter
                if (a(iatm)%name==elemplot.or.elemplot=="all") then
                    sval=atmshdall(iatm,isystem)
                    if (igetshift==0) then
                        write(ides,"(i5,'(',a,')',4x,f12.3)") iatm,a(iatm)%name,sval
                    else if (igetshift==1) then
                        write(ides,"(i5,'(',a,')',4x,f12.3,9x,f12.3)") iatm,a(iatm)%name,sval,NMRref-sval
                    else if (igetshift==2) then
                        write(ides,"(i5,'(',a,')',4x,f12.3,9x,f12.3)") iatm,a(iatm)%name,sval,(sval-NMRinter)/NMRslope
                    end if
                end if
            end do
        end do
        if (any(weight/=1)) then
            write(ides,*)
            write(ides,*) "Weighted data:"
            if (igetshift==0) then
                write(ides,*) "   Atom       Shielding(iso)"
            else
                write(ides,*) "   Atom       Shielding(iso)      Chemical Shift"
            end if
            do iatm=1,ncenter
                if (a(iatm)%name==elemplot.or.elemplot=="all") then
                    sval=atmshdwei(iatm)
                    if (igetshift==0) then
                        write(ides,"(i5,'(',a,')',4x,f12.3)") iatm,a(iatm)%name,sval
                    else if (igetshift==1) then
                        write(ides,"(i5,'(',a,')',4x,f12.3,9x,f12.3)") iatm,a(iatm)%name,sval,NMRref-sval
                    else if (igetshift==2) then
                        write(ides,"(i5,'(',a,')',4x,f12.3,9x,f12.3)") iatm,a(iatm)%name,sval,(sval-NMRinter)/NMRslope
                    end if
                end if
            end do
        end if
        if (isel==-2) then
            write(*,*) "Done! Data have been exported to NMRdata.txt in current folder"
            close(10)
        end if
    else if (isel==3) then
		write(*,*) "Input lower limit, upper limit and step between ticks"
        write(*,*) "Example 1: -20,85,10     Example 2: 150,0,-20"
		write(*,*) "Hint: If only input 0, the axis will be inverted"
		read(*,"(a)") c200tmp
		if (c200tmp=='0') then
			tmp=xlow
			xlow=xhigh
			xhigh=tmp
			stepx=-stepx
		else
			read(c200tmp,*,iostat=ierror) xlow,xhigh,stepx
            if (ierror/=0) then
                write(*,*) "Error: Unable to recognize the inputted content. Press ENTER button to continue"
                read(*,*)
                cycle
            end if
            if (xlow>xhigh.and.stepx>0) stepx=-stepx
		end if
		iusersetX=1 !User has modified it
    
    else if (isel==4) then
		write(*,*) "Input lower limit, upper limit and step between ticks e.g. 0,3.5,1"
		read(*,"(a)") c200tmp
		read(c200tmp,*,iostat=ierror) orgy1,endy1,stepy1
        if (ierror/=0) then
            write(*,*) "Error: Unable to recognize the inputted content. Press ENTER button to continue"
            read(*,*)
            cycle
        end if
		iusersetY1=1
        
    else if (isel==5) then
		write(*,*) "Input lower limit, upper limit and step between ticks e.g. 0,8,0.5"
		read(*,"(a)") c200tmp
		read(c200tmp,*,iostat=ierror) orgy2,endy2,stepy2
        if (ierror/=0) then
            write(*,*) "Error: Unable to recognize the inputted content. Press ENTER button to continue"
            read(*,*)
            cycle
        end if
		iusersetY2=1
    
    else if (isel==6) then
        write(*,*) "Input the element, e.g. Cl"
        write(*,*) "You can input ""all"" to choose all elements"
        read(*,*) elemplot
        igetshift=0
        iusersetX=0
        if (elemplot=="all") then
            ilabelele=1
            FWHM_NMR=0.2D0
        else
            ilabelele=0
            if (elemplot=='H') then
                FWHM_NMR=0.02D0
            else
                FWHM_NMR=0.2D0
            end if
        end if
        
    else if (isel==7) then
        write(*,*) "0 Do not calculate chemical shift"
        write(*,*) "1 Set reference shielding value to determine chemical shift"
        write(*,*) "2 Set slope and intercept to determine chemical shift by scaling method"
        read(*,*) igetshift
        if (igetshift==1) then
            write(*,*) "Input reference value in ppm, e.g. 120.5"
            write(*,"(a)") " You can also input corresponding letter to use built-in values of TMS. &
            a~f were all calculated under chloroform represented by SMD model, geometries were optimized at B3LYP/def2-SVP in vaccum"
            write(*,*) "a: B97-2/def2-TZVP G09 (186.8707 for C and 31.5143 for H)"
            write(*,*) "b: B97-2/pcSseg-1 G09 (184.4144 for C and 31.2876 for H)"
            write(*,*) "c: MP2/pcSseg-1 G09 (195.6338 for C and 31.2732 for H)"
            write(*,*) "d: B97-2/def2-TZVP G16 (186.7595 for C and 31.5058 for H)"
            write(*,*) "e: B97-2/pcSseg-1 G16 (184.2007 for C and 31.2802 for H)"
            write(*,*) "f: MP2/pcSseg-1 G16 (195.5706 for C and 31.2652 for H)"
            write(*,"(a)") " The following one was calculated under chloroform represented by IEFPCM model, geometry was optimized at B3LYP/6-31G* in vaccum"
            write(*,*) "g: revTPSS/pcSseg-1 G16 (183.6902 for C and 31.7306 for H)"
            read(*,*) c80tmp
            if (iachar(c80tmp(1:1))>=48.and.iachar(c80tmp(1:1))<=57) then
                read(c80tmp,*) NMRref
            else
                if (elemplot=="C".or.elemplot=="H") then
                    if (c80tmp=="a") then
                        if (elemplot=="C") NMRref=186.8707D0
                        if (elemplot=="H") NMRref=31.5143D0
                    else if (c80tmp=="b") then
                        if (elemplot=="C") NMRref=184.4144D0
                        if (elemplot=="H") NMRref=31.2876D0
                    else if (c80tmp=="c") then
                        if (elemplot=="C") NMRref=195.6338D0
                        if (elemplot=="H") NMRref=31.2732D0
                    else if (c80tmp=="d") then
                        if (elemplot=="C") NMRref=186.7595D0
                        if (elemplot=="H") NMRref=31.5058D0
                    else if (c80tmp=="e") then
                        if (elemplot=="C") NMRref=184.2007D0
                        if (elemplot=="H") NMRref=31.2802D0
                    else if (c80tmp=="f") then
                        if (elemplot=="C") NMRref=195.5706D0
                        if (elemplot=="H") NMRref=31.2652D0
                    else if (c80tmp=="g") then
                        if (elemplot=="C") NMRref=183.6902D0
                        if (elemplot=="H") NMRref=31.7306D0
                    end if
                else
                    write(*,"(a)") " Error: The current element is neither C nor H! You should use option 6 to &
                    properly select the element to be considered"
                end if
            end if
        else if (igetshift==2) then
            write(*,*) "Input slope and intercept, e.g. -0.9,120.5"
            write(*,"(a)") " If inputting ""a"", the parameters for C and H fitted at B3LYP/6-31G* level with chloroform represented by SMD will be used"
            read(*,"(a)") c80tmp
            if (index(c80tmp,'a')==0) then
                read(c80tmp,*) NMRslope,NMRinter
            else
                if (elemplot=="C".or.elemplot=="H") then
                    if (elemplot=="C") then
                        NMRslope=-0.9449D0
                        NMRinter=188.4418D0
                    else if (elemplot=="H") then
                        NMRslope=-1.0157D0
                        NMRinter=32.210D0
                    end if
                    write(*,*) "Note: The parameters are taken from http://cheshirenmr.info" 
                else
                    write(*,"(a)") " Error: The current element is neither C nor H! You should use option 6 to &
                    properly select the element to be considered"
                end if
            end if
        end if
        
    else if (isel==8) then
        write(*,*) "Input full width at half maximum (FWHM) in ppm, e.g. 0.8"
        read(*,*) FWHM_NMR
        
    else if (isel==9) then
        write(*,*) "Input tolerance for determining NMR degeneracy, e.g. 0.1"
        read(*,*) degentol
        
    else if (isel==10) then
        write(*,*) "Input the indices of the atom for which the shieldings will be averaged"
        write(*,*) "e.g. 3,6,8-12,44"
        if (nsystem>1) write(*,*) "Average will be applied to all systems"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp)
        allocate(tmparrint(ntmp))
        call str2arr(c2000tmp,ntmp,tmparrint)
        do isystem=1,nsystem
            avgval=sum(atmshdall(tmparrint(:),isystem))/ntmp
            atmshdall(tmparrint(:),isystem)=avgval
            if (allocated(atmshdwei)) then
                avgval=sum(atmshdwei(tmparrint(:)))/ntmp
                atmshdwei(tmparrint(:))=avgval
            end if
        end do
        deallocate(tmparrint)
        write(*,*) "Done!"
        
    else if (isel==11) then
        if (all(atmstr==1)) then
            write(*,*) "Note: All atoms have strength of 1"
        else
            write(*,*) "Atom(s) having strength other than 1:"
            do iatm=1,ncenter
                if (atmstr(iatm)/=1) then
                    write(*,"(' Atom',i5,' (',a,')  Strength:',f10.5)") iatm,a(iatm)%name,atmstr(iatm)
                end if
            end do
        end if
        write(*,*)
        write(*,*) "Input indices of the atoms whose strengths will be set, e.g. 12,14-17,20"
        write(*,*) "Input ""q"" can return"
        read(*,"(a)") c2000tmp
        if (c2000tmp=='q') cycle
        write(*,*) "Set to which value? e.g. 0.5"
        write(*,*) "Note: If set to 0, then these atoms will be fully invisible in NMR spectrum"
        read(*,*) strtmp
        call str2arr(c2000tmp,ntmp)
        allocate(tmparrint(ntmp))
        call str2arr(c2000tmp,ntmp,tmparrint)
        do itmp=1,ntmp
            atmstr(tmparrint(itmp))=strtmp
        end do
        deallocate(tmparrint)
        
    else if (isel==12) then
        if (ishowline==1) then
            ishowline=0
        else
            ishowline=1
        end if
        
    else if (isel==13) then
        if (ishowcurve==1) then
            ishowcurve=0
        else
            ishowcurve=1
        end if
        
    else if (isel==14) then
        do while(.true.)
            write(*,*)
            if (nsystem==1) then
                write(*,*) "0 Return"
                write(*,"(a)") " 1 Set color for spikes, current: "//trim(colorname(line_clr(1)))
                write(*,"(a)") " 2 Set color for curve, current: "//trim(colorname(curve_clr(1)))
                read(*,*) isel2
                if (isel2==0) then
                    exit
                else if (isel2==1) then
                    call selcolor(line_clr(1))
                else if (isel2==2) then
                    call selcolor(curve_clr(1))
                end if
            else
                write(*,"(a)") " -2 Set color for weighted curve, current: "//trim(colorname(curve_clr(0)))
                write(*,"(a)") " -1 Set color for weighted spikes, current: "//trim(colorname(line_clr(0)))
                write(*,*) "0 Return"
                do imol=1,nsystem
                    if (weight(imol)==1) then
                        write(*,"(i2,' Set color for system',i2,' (',a,'), line: ',a,', curve: ',a)") &
                        imol,imol,trim(mollegend(imol)),trim(colorname(line_clr(imol))),trim(colorname(curve_clr(imol)))
                    else
                        write(*,"(i2,' Set color for system',i2,' (weight=',f7.4,'), line: ',a,', curve: ',a)") &
                        imol,imol,weight(imol),trim(colorname(line_clr(imol))),trim(colorname(curve_clr(imol)))
                    end if
                end do
                read(*,*) isel2
                if (isel2==0) then
                    exit
                else if (isel2==-1) then
                    call selcolor(line_clr(0))
                else if (isel2==-2) then
                    call selcolor(curve_clr(0))
                else
                    write(*,*) "Set color for spikes:"
                    call selcolor(line_clr(isel2))
                    write(*,*) "Set color for curves:"
                    call selcolor(curve_clr(isel2))
                end if
            end if
        end do
    
    else if (isel==15) then
		do while(.true.)
			write(*,*)
			write(*,*) "0 Return"
			write(*,"(' 1 Set thickness of curves, current:',i3)") thk_curve
			if (any(weight/=1)) write(*,"(' 2 Set thickness of weighted curve, current:',i3)") thk_weicurve
			write(*,"(' 3 Set thickness of spikes, current:',i3)") thk_line
			if (any(weight/=1)) write(*,"(' 4 Set thickness of weighted spike, current:',i3)") thk_weiline
			write(*,"(' 6 Set thickness of axes, current:',i3)") thk_axis
			write(*,"(' 7 Set thickness of grid, current:',i3)") thk_grid
			read(*,*) isel2
			if (isel2==0) exit
			write(*,*) "Input the thickness, e.g. 3"
			if (isel2==1) read(*,*) thk_curve
			if (isel2==2) read(*,*) thk_weicurve
			if (isel2==3) read(*,*) thk_line
			if (isel2==4) read(*,*) thk_weiline
			if (isel2==6) read(*,*) thk_axis
			if (isel2==7) read(*,*) thk_grid
			write(*,*) "Done!"
		end do
    
    else if (isel==16) then
        do while(.true.)
            write(*,*)
            write(*,*) "0 Return"
            if (ishowlabel==0) write(*,*) "1 Toggle showing atomic labels, current: No"
            if (ishowlabel==1) write(*,*) "1 Toggle showing atomic labels, current: Yes"
            write(*,"(a,i3)") " 2 Set label size, current:",labelsize
            write(*,*) "3 Set label color, current: "//trim(colorname(ilabelclr))
            write(*,"(a,i4)") " 4 Set shift of labels in X, current:",labelshiftX
            write(*,"(a,i4)") " 5 Set shift of labels in Y, current:",labelshiftY
            if (ialllabtop==0) write(*,*) "6 Toggle showing all atomic labels at top of spikes, current: No"
            if (ialllabtop==1) write(*,*) "6 Toggle showing all atomic labels at top of spikes, current: Yes"
            if (ilabelele==0) write(*,*) "7 Toggle showing element in labels, current: No"
            if (ilabelele==1) write(*,*) "7 Toggle showing element in labels, current: Yes"
            if (any(weight/=1)) then
                if (ionlyatmlabwei==0) write(*,*) "8 Only show labels for weighted data, current: No"
                if (ionlyatmlabwei==1) write(*,*) "8 Only show labels for weighted data, current: Yes"
            end if
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
                if (ishowlabel==0) then
                    ishowlabel=1
                else
                    ishowlabel=0
                end if
            else if (isel2==2) then
                write(*,*) "Input label size, e.g. 30"
                read(*,*) labelsize
            else if (isel2==3) then
                call selcolor(ilabelclr)
            else if (isel2==4) then
                write(*,*) "Input shift of labels in X, e.g. 5"
                read(*,*) labelshiftX
            else if (isel2==5) then
                write(*,*) "Input shift of labels in Y, e.g. 5"
                read(*,*) labelshiftY
            else if (isel2==6) then
                if (ialllabtop==0) then
                    ialllabtop=1
                else
                    ialllabtop=0
                end if
            else if (isel2==7) then
                if (ilabelele==1) then
                    ilabelele=0
                else
                    ilabelele=1
                end if
            else if (isel2==8) then
                if (ionlyatmlabwei==0) then
                    ionlyatmlabwei=1
                else
                    ionlyatmlabwei=0
                end if
            end if
        end do
    
    else if (isel==17) then
        write(*,*) "0: Only show spectra of individual systems"
        write(*,*) "1: Show both weighted spectrum and spectra of individual systems"
        write(*,*) "2: Only show weighted spectrum"
		read(*,*) ishowweighted
    
    else if (isel==18) then
        do while(.true.)
			write(*,*)
            write(*,*) "0 Return"
            write(*,"(a,i3)") " 1 Set number of ticks in X-axis, current:",nticksX
            write(*,"(a,i3)") " 2 Set number of ticks in the axis of signal strength, current:",nticksY
	        if (ishowgrid==1) write(*,*) "3 Toggle showing dashed grid lines, current: ON"
	        if (ishowgrid==0) write(*,*) "3 Toggle showing dashed grid lines, current: OFF"
            if (nsystem>1) then
                if (ishowlegend==0) write(*,*) "4 Toggle showing legends, current: No"
                if (ishowlegend==1) write(*,*) "4 Toggle showing legends, current: Yes"
                write(*,"(a,i5)") " 5 Set X position of legends, current:",legposx
                write(*,"(a,i5)") " 6 Set Y position of legends, current:",legposy
                write(*,"(a,i5)") " 7 Set text size of legends, current:",legendsize
            end if
            if (ishowdataright==0) write(*,*) "8 Toggle showing labels on the Y-axis of signal strength, current: No"
            if (ishowdataright==1) write(*,*) "8 Toggle showing labels on the Y-axis of signal strength, current: Yes"
            read(*,*) isel2
        
            if (isel2==0) then
                exit
            else if (isel2==1) then
                write(*,*) "Input number of ticks, e.g. 3"
                read(*,*) nticksX
            else if (isel2==2) then
                write(*,*) "Input number of ticks, e.g. 3"
                read(*,*) nticksY
            else if (isel2==3) then
                if (ishowgrid==1) then
                    ishowgrid=0
                else
                    ishowgrid=1
                end if
            else if (isel2==4) then
                if (ishowlegend==1) then
                    ishowlegend=0
                else
                    ishowlegend=1
                end if
            else if (isel2==5) then
                write(*,*) "Input the X position, e.g. 1900"
                write(*,*) "Note: The smaller the number, the more left of the position of the legends"
                read(*,*) legposx
            else if (isel2==6) then
                write(*,*) "Input the Y position, e.g. 200"
                write(*,*) "Note: The smaller the number, the higher the position of the legends"
                read(*,*) legposy
            else if (isel2==7) then
                write(*,*) "Input text size of legends, e.g. 50"
                read(*,*) legendsize
            else if (isel2==8) then
                if (ishowdataright==1) then
                    ishowdataright=0
                else
                    ishowdataright=1
                end if
            end if
        end do
        
	!!================================================================!!
	!!================================================================!!
	!!=============== Functions need generation of curves ============!!
	!!================================================================!!
	!!================================================================!!
	else if (isel==0.or.isel==1.or.isel==2) then
        if (isel==0) then
            isavepic=0
        else if (isel==1) then
            isavepic=1
        end if
        
        !====== Convert absolute shielding to chemical shift. After plotting, convert back
        atmshdall_org=atmshdall !Back up
        if (ishowweighted/=0) atmshdwei_org=atmshdwei !Back up
        if (igetshift==0) then
            continue
        else if (igetshift==1) then
            atmshdall=NMRref-atmshdall
            if (ishowweighted/=0) atmshdwei=NMRref-atmshdwei
        else if (igetshift==2) then
            atmshdall=(atmshdall-NMRinter)/NMRslope
            if (ishowweighted/=0) atmshdwei=(atmshdwei-NMRinter)/NMRslope
        end if
        
        !====== Generate (non)degenerate terms for plotting from atomic shielding values
        shdnatm=0
        shdeffnatm=0
        shdnum=0
        do imol=1,nsystem
            do iatm=1,ncenter
                if (a(iatm)%name/=elemplot.and.elemplot/="all") cycle
                shdthis=atmshdall(iatm,imol)
                effstr=atmstr(iatm)
                if (shdnum(imol)==0) then
                    shdnum(imol)=1
                    shdval(1,imol)=shdthis
                    shdnatm(1,imol)=1
                    shdeffnatm(1,imol)=effstr
                    shdatm(1,1,imol)=iatm
                else !Compare shielding of this atom with existing ones, if very close, then merge it
                    do iterm=1,shdnum(imol)
                        diff=abs(shdval(iterm,imol)-shdthis)
                        if (diff<degentol) then !Degenerate to old term
                            shdnatm(iterm,imol)=shdnatm(iterm,imol)+1
                            shdeffnatm(iterm,imol)=shdeffnatm(iterm,imol)+effstr
                            shdatm(shdnatm(iterm,imol),iterm,imol)=iatm
                            exit
                        end if
                    end do
                    if (iterm==shdnum(imol)+1) then !New term
                        nterm=shdnum(imol)+1
                        shdnum(imol)=nterm
                        shdval(nterm,imol)=shdthis
                        shdnatm(nterm,imol)=1
                        shdeffnatm(nterm,imol)=effstr
                        shdatm(1,nterm,imol)=iatm
                    end if
                end if
            end do
            !Output to screen
            if (isel==0.or.isel==1) then
                if (nsystem>1) write(*,"(/,' System',i5,':')") imol
                do iterm=1,shdnum(imol)
                    if (igetshift==0) then
                        write(*,"(' Term:',i5,'   Shielding(iso):',f10.3,' ppm   Atom:')",advance="no") iterm,shdval(iterm,imol)
                    else
                        write(*,"(' Term:',i5,'   Chemical shift:',f10.3,' ppm   Atom:')",advance="no") iterm,shdval(iterm,imol)
                    end if
                    do idx=1,shdnatm(iterm,imol)
                        iatm=shdatm(idx,iterm,imol)
                        write(*,"(i5,'(',a,')')",advance="no") iatm,a(iatm)%name
                    end do
                    if (any(atmstr/=1)) then
                        write(*,"(' Strength:',f6.3)") shdeffnatm(iterm,imol)
                    else
                        write(*,*)
                    end if
                end do
            end if
        end do
        
        !====== Generate (non)degenerate terms for plotting weighted spectrum
        if (ishowweighted/=0) then
            shdnatmwei=0
            shdeffnatmwei=0
            shdnumwei=0
            do iatm=1,ncenter
                if (a(iatm)%name/=elemplot.and.elemplot/="all") cycle
                shdthis=atmshdwei(iatm)
                effstr=atmstr(iatm)
                if (shdnumwei==0) then
                    shdnumwei=1
                    shdvalwei(1)=shdthis
                    shdnatmwei(1)=1
                    shdeffnatmwei(1)=effstr
                    shdatmwei(1,1)=iatm
                else
                    do iterm=1,shdnumwei
                        diff=abs(shdvalwei(iterm)-shdthis)
                        if (diff<degentol) then !Degenerate to old term
                            shdnatmwei(iterm)=shdnatmwei(iterm)+1
                            shdeffnatmwei(iterm)=shdeffnatmwei(iterm)+effstr
                            shdatmwei(shdnatmwei(iterm),iterm)=iatm
                            exit
                        end if
                    end do
                    if (iterm==shdnumwei+1) then !New term
                        shdnumwei=shdnumwei+1
                        shdvalwei(shdnumwei)=shdthis
                        shdnatmwei(shdnumwei)=1
                        shdeffnatmwei(shdnumwei)=effstr
                        shdatmwei(1,shdnumwei)=iatm
                    end if
                end if
            end do
            !Output to screen
            if (isel==0.or.isel==1) then
                write(*,"(/,' Weighted data:')")
                do iterm=1,shdnumwei
                    if (igetshift==0) then
                        write(*,"(' Term:',i5,'   Shielding(iso):',f10.3,' ppm   Atom:')",advance="no") iterm,shdvalwei(iterm)
                    else
                        write(*,"(' Term:',i5,'   Chemical shift:',f10.3,' ppm   Atom:')",advance="no") iterm,shdvalwei(iterm)
                    end if
                    do idx=1,shdnatmwei(iterm)
                        iatm=shdatmwei(idx,iterm)
                        write(*,"(i5,'(',a,')')",advance="no") iatm,a(iatm)%name
                    end do
                    if (any(atmstr/=1)) then
                        write(*,"(' Strength:',f6.3)") shdeffnatmwei(iterm)
                    else
                        write(*,*)
                    end if
                end do
            end if
        end if
        
		!====== Determine upper and lower limit of X axis
		if (iusersetX==0) then !Automatical scale, if user has not manually set the range
			tmpmin=minval(shdval(1:shdnum(1),1))
			tmpmax=maxval(shdval(1:shdnum(1),1))
			if (nsystem>1) then !Find upper and lower values for all systems
				do imol=2,nsystem
					tmpa=minval(shdval(1:shdnum(imol),imol))
					tmpb=maxval(shdval(1:shdnum(imol),imol))
					if (tmpa<tmpmin) tmpmin=tmpa
					if (tmpb>tmpmax) tmpmax=tmpb
				end do
			end if
            if (tmpmax==tmpmin) then
                xlow=floor(tmpmin-1)
                xhigh=ceiling(tmpmax+1)
            else
                extend=(tmpmax-tmpmin)/5
                xlow=floor(tmpmin-extend)
                xhigh=ceiling(tmpmax+extend)
            end if
            tmp=(xhigh-xlow)/10
            if (tmp>1) then !Make step close to 1,2,5,10,20,50
                if (tmp<1.5D0) then
                    stepx=1
                else if (tmp<3.5D0) then
                    stepx=2
                else if (tmp<7.5D0) then
                    stepx=5
                else if (tmp<15) then
                    stepx=10
                else if (tmp<35) then
                    stepx=20
                else if (tmp<75) then
                    stepx=50
                else
                    stepx=nint(tmp)
                end if
            else
                stepx=tmp
            end if
            if (igetshift>0) then !The data to be plotted is chemical shift rather than absolute shielding, lower limit of X-axis corresponds to larger shift
                tmp=xlow
                xlow=xhigh
                xhigh=tmp
                stepx=-stepx
            end if
		end if
        
		!====== Set x positions of curves
		xptstep=(xhigh-xlow)/(num1Dpoints-1)
		do ipoint=1,num1Dpoints
			curvex(ipoint)=xlow+(ipoint-1)*xptstep
		end do
        
		!====== Generate discrete line
		lineyall=0D0
		do imol=1,nsystem
			do iterm=1,shdnum(imol)
				inow=3*(iterm-1)
				linexall(imol,inow+1:inow+3)=shdval(iterm,imol)
				lineyall(imol,inow+2)=shdeffnatm(iterm,imol)
			end do
		end do
        if (ishowweighted/=0) then
			do iterm=1,shdnumwei
				inow=3*(iterm-1)
				linexwei(inow+1:inow+3)=shdvalwei(iterm)
                lineywei(inow+1)=0D0
				lineywei(inow+2)=shdeffnatmwei(iterm)
                lineywei(inow+3)=0D0
			end do
        end if
        
		!====== Broaden to curve
		!Under current X and Y axes units, below code guarantees that the integral of the peak broadened by one unit of strength is 1
		!Lorentzian function, see http://mathworld.wolfram.com/LorentzianFunction.html
		curveyall=0D0
		do imol=1,nsystem
			do iterm=1,shdnum(imol)
				preterm=shdeffnatm(iterm,imol)*0.5D0/pi*FWHM_NMR !Integral of the peak equals to degeneracy
				do ipoint=1,num1Dpoints
					curveytmp(ipoint)=preterm/( (curvex(ipoint)-shdval(iterm,imol))**2+0.25D0*FWHM_NMR**2 )
				end do
				curveyall(imol,:)=curveyall(imol,:)+curveytmp
			end do
		end do
        curveywei=0D0
        if (ishowweighted/=0) then
            do iterm=1,shdnumwei
				preterm=shdeffnatmwei(iterm)*0.5D0/pi*FWHM_NMR !Integral of the peak equals to degeneracy
				do ipoint=1,num1Dpoints
					curveytmp(ipoint)=preterm/( (curvex(ipoint)-shdvalwei(iterm))**2+0.25D0*FWHM_NMR**2 )
				end do
				curveywei(:)=curveywei(:)+curveytmp
			end do
        end if
        
        atmshdall=atmshdall_org !Restore to absolute shielding
        if (ishowweighted/=0) atmshdwei=atmshdwei_org !Restore to absolute shielding
        
        !====== Export NMR curves and discrete lines to text files
		if (isel==2) then
			if (nsystem==1) then
				open(10,file="NMR_line.txt",status="replace")
				do ipt=1,3*shdnum(1)
					write(10,"(2f13.5)") linexall(1,ipt),lineyall(1,ipt)
				end do
				close(10)
				write(*,*) "Discrete line data has been written to NMR_line.txt in current folder"
				open(10,file="NMR_curve.txt",status="replace")
				do ipt=1,num1Dpoints
					write(10,"(2f13.5)") curvex(ipt),curveyall(1,ipt)
				end do
				close(10)
				write(*,*) "Curve data has been written to NMR_curve.txt in current folder"
			else !Also output curve for all systems
				if (ishowweighted/=0) then !Output weighted spectrum
					open(10,file="NMR_linewei.txt",status="replace")
					do ipt=1,3*shdnumwei
						write(10,"(2f13.5)") linexwei(ipt),lineywei(ipt)
					end do
					close(10)
					write(*,"(a)") " The discrete line data corresponding to weighted spectrum has been written to NMR_linewei.txt in current folder"
					open(10,file="NMR_curvewei.txt",status="replace")
					do ipt=1,num1Dpoints
						write(10,"(2f13.5)") curvex(ipt),curveywei(ipt)
					end do
					close(10)
					write(*,"(a)") " The curve data corresponding to weighted spectrum has been written to NMR_curvewei.txt in current folder"
				end if
				open(10,file="NMR_curveall.txt",status="replace")
				do ipt=1,num1Dpoints
					write(10,"(f13.5)",advance="no") curvex(ipt)
					do imol=1,nsystem
						write(10,"(f13.5)",advance="no") curveyall(imol,ipt)
					end do
					write(10,*)
				end do
				close(10)
				write(*,"(a)") " Curve data of all systems have been exported to NMR_curveall.txt in current folder as different columns"
			end if
        
        !=================================
        !========= Draw NMR map ==========
        !=================================
        else if (isel==0.or.isel==1) then
            !Set default lower and upper limit
		    if (iusersetY1==0) then !Left Y axis (degeneracy)
                orgy1=0
			    endy1=maxval(shdeffnatm(:,:))+0.5D0
			    stepy1=0.5D0
		    end if
		    if (iusersetY2==0) then !Y axis of signal strength
			    orgy2=0
			    endy2=1.1D0*maxval(curveyall(:,:))
			    stepy2=(endy2-orgy2)/10
		    end if
            
            !Initialize DISLIN
            if (isavepic==0) then
			    call METAFL('xwin')
			    call window(200,100,1200,750) !1.6:1
		    else if (isavepic==1) then
			    call METAFL(graphformat)
			    call window(200,100,graph1Dwidth,graph1Dheight) !Default is 1.6:1
		    end if
		    call SCRMOD('REVERSE')
		    CALL IMGFMT("RGB")
		    CALL PAGE(3000,1875) !1.6:1
		    call disini
		    if (isavepic==0) call WINTIT("Click right mouse button to close")
		    call ERRMOD("ALL","OFF")
		    if (isavepic==1.and.graphformat=="png") then
			    call TRIPLX
		    else
			    CALL HWFONT
		    end if
            call axspos(350,1600)
            nxpixel=2300
            nypixel=1470
            call AXSLEN(nxpixel,nypixel)
            
		    !!! Set and draw X-axis and frame
            call namdis(45,'X')
		    CALL TICKS(nticksX,'X')
            CALL TICKS(0,'Y')
		    CALL LABDIG(-1,"X")
            if (abs(stepx)<1) CALL LABDIG(1,"X")
            if (abs(stepx)<0.1D0) CALL LABDIG(2,"X")
            CALL HNAME(40)
            if (igetshift==0) then
                CALL NAME('Shielding (ppm)','X')
            else
                CALL NAME('Chemical shift (ppm)','X')
            end if
            call height(40)
		    call setgrf('NAME','TICKS','NONE','TICKS')
		    call LINWID(thk_axis)
		    CALL GRAF(xlow,xhigh,xlow,stepx, orgy1,endy1,orgy1,stepy1)
            call endgrf !End of plotting X-axis and frame, as well as legends
            
            !!! Draw discrete lines !!!
            if (ishowline==1) then
		        !Set properties of left Y-axis (degeneracy)
		        call setgrf('NONE','NAME','NONE','NONE')
		        CALL TICKS(1,'Y')
		        CALL LABDIG(1,"Y")
                call namdis(60,'Y')
                CALL HNAME(40)
		        CALL NAME('Degeneracy','Y')
		        call LINWID(thk_axis)
		        CALL GRAF(xlow,xhigh,xlow,stepx, orgy1,endy1,orgy1,stepy1)
                
                !Draw shallow gray dashed lines
		        if (ishowgrid==1) then 
			        call SETRGB(0.85D0,0.85D0,0.85D0) !Shallow gray
			        call dash
			        call LINWID(thk_grid)
			        call grid(1,1)
		            call solid
		        end if
		        !Draw weighted discrete lines
		        if (ishowweighted/=0) then
			        call setcolor(line_clr(0))
			        CALL LINWID(thk_weiline)
			        CALL CURVE(linexwei(:),lineywei(:),3*shdnumwei)
		        end if
		        !Draw discrete line for every system
                if (ishowweighted/=2) then
			        CALL LINWID(thk_line)
			        do imol=1,nsystem
                        call setcolor(line_clr(imol))
                        ndata=3*shdnum(imol)
				        CALL CURVE(linexall(imol,1:ndata),lineyall(imol,1:ndata),ndata)
			        end do
                end if
			    call xaxgit !Draw a line corresponding to Y=0
		        call LINWID(thk_axis)
                call endgrf !End of plotting discrete lines
			    call color("WHITE")
		    end if
            
            !!! Draw curves !!!
            if (ishowcurve==1) then
                if (ishowline==0) then !Show at left side
                    call setgrf('NONE','NAME','NONE','NONE')
			    else if (ishowline==1) then !Show at right side
                    call setgrf('NONE','NONE','NONE','NAME')
                end if
		        CALL TICKS(nticksY,'Y')
                CALL LABDIG(1,"Y")
                if (ishowdataright==0) then
                    call labels("NONE","Y")
		            CALL TICKS(0,'Y')
                end if
                call namdis(60,'Y')
                CALL HNAME(40)
			    CALL NAME('Signal strength','Y')
			    CALL GRAF(xlow,xhigh,xlow,stepx, orgy2,endy2,orgy2,stepy2)
		        !Draw weighted curves
                if (ishowweighted/=0) then
			        CALL LINWID(thk_weicurve)
			        call setcolor(curve_clr(0))
                    CALL CURVE(curvex(:),curveywei(:),num1Dpoints)
                end if
		        !Draw curves for every system
                if (ishowweighted/=2) then
			        CALL LINWID(thk_curve)
			        do imol=1,nsystem
			            call setcolor(curve_clr(imol))
				        CALL CURVE(curvex(:),curveyall(imol,:),num1Dpoints)
			        end do
                end if
		        call LINWID(thk_axis)
                call endgrf
			    call color("WHITE")
            end if
            
            !!! Show atomic labels !!!
            if (ishowlabel==1) then
		        !Set properties of left Y-axis (degeneracy)
		        call setgrf('NONE','NONE','NONE','NONE')
		        CALL GRAF(xlow,xhigh,xlow,stepx, orgy1,endy1,orgy1,stepy1)
                pix2usrX=(xhigh-xlow)/nxpixel !Used to shift label position in X
                pix2usrY=(endy1-orgy1)/nypixel !Used to shift label position in Y
                call height(labelsize)
                call setcolor(ilabelclr)
                !Atoms labels on weighted lines
                if (ishowweighted/=0) then !Weighted lines are showing
			        do iterm=1,shdnumwei
                        if (shdeffnatmwei(iterm)==0) cycle
                        Xpos=shdvalwei(iterm)
                        if (Xpos<min(xlow,xhigh).or.Xpos>max(xlow,xhigh)) cycle
                        do idx=1,shdnatmwei(iterm)
                            if (ialllabtop==1) then
				                Ypos=shdeffnatmwei(iterm)+idx*labelsize*1*pix2usrY
                            else
    				            Ypos=sum(atmstr(shdatmwei(1:idx,iterm)))+labelsize*1.4*pix2usrY
                            end if
                            iatm=shdatmwei(idx,iterm)
                            write(c80tmp,*) iatm
                            shiftXtmp=labelshiftX
                            if (iatm>=10) shiftXtmp=shiftXtmp-20
                            if (ilabelele==1) then
                                c80tmp=trim(adjustl(c80tmp))//a(iatm)%name
                                shiftXtmp=shiftXtmp-20
                            end if
                            call rlmess(trim(adjustl(c80tmp)),Xpos+shiftXtmp*pix2usrX,Ypos+labelshiftY*pix2usrY)
                        end do
			        end do
                end if
                !Atoms labels on individual systems
                if (nsystem==1.or.(nsystem>1.and.ishowweighted/=2)) then
                    if (.not.(any(weight/=1).and.ionlyatmlabwei==1)) then
                        do imol=1,nsystem
			                do iterm=1,shdnum(imol)
                                if (shdeffnatm(iterm,imol)==0) cycle
                                Xpos=shdval(iterm,imol)
                                if (Xpos<min(xlow,xhigh).or.Xpos>max(xlow,xhigh)) cycle
                                do idx=1,shdnatm(iterm,imol)
                                    if (ialllabtop==1) then
				                        Ypos=shdeffnatm(iterm,imol)+idx*labelsize*1.4*pix2usrY
                                    else
                                        Ypos=sum(atmstr(shdatm(1:idx,iterm,imol)))+labelsize*1*pix2usrY
                                    end if
                                    iatm=shdatm(idx,iterm,imol)
                                    write(c80tmp,*) iatm
                                    shiftXtmp=labelshiftX
                                    if (iatm>=10) shiftXtmp=shiftXtmp-20
                                    if (ilabelele==1) then
                                        c80tmp=trim(adjustl(c80tmp))//a(iatm)%name
                                        shiftXtmp=shiftXtmp-20
                                    end if
                                    call rlmess(trim(adjustl(c80tmp)),Xpos+shiftXtmp*pix2usrX,Ypos+labelshiftY*pix2usrY)
                                end do
			                end do
                        end do
                    end if
                end if
                call endgrf
                call setcolor(5) !Recover to black
            end if
            
            !Show legends
            if (nsystem>1.and.ishowlegend==1) then
		        call setgrf('NONE','NONE','NONE','NONE')
		        CALL GRAF(xlow,xhigh,xlow,stepx, orgy1,endy1,orgy1,stepy1)
		        numleg=1+nsystem
		        call legini(clegend,numleg,50)
		        call legtit(' ')
		        call frame(0) !No box around legend
                call height(legendsize)
                !Set legend position
                legposxtmp=legposx
		        if (allocated(mollegend)) then
			        maxlen=maxval(len_trim(mollegend))
			        legposxtmp=legposxtmp-maxlen*22
		        end if
		        if (isavepic==1) legposxtmp=legposxtmp-100
		        call legpos(legposxtmp,legposy)
                ileg=0
                !Show legend for weighted spectrum
                if (ishowweighted==1) then
			        ileg=ileg+1
                    call setcolor(line_clr(0))
                    call legpat(0,1,-1,-1,-1,ileg)
                    CALL LEGLIN(clegend," Weighted",ileg)
                end if
                !Show legend for every spectrum
                if (ishowweighted/=2) then
			        do imol=1,nsystem
			            ileg=ileg+1
                        call setcolor(line_clr(imol))
			            call legpat(0,1,-1,-1,-1,ileg)
			            CALL LEGLIN(clegend,trim(mollegend(imol)),ileg)
                    end do
                end if
                call legopt(2.5D0,0.5D0,1D0) !Decrease the length of legend color line
			    call color("WHITE")
		        if (ileg>0) call legend(clegend,3)
                call endgrf !End of plotting X-axis and frame, as well as legends
            end if
            
            call setcolor(5) !Recover to black
            call height(36) !Recover to default
		    call disfin
		    if (isavepic==1) write(*,*) "Graphical file has been saved to current folder with ""dislin"" prefix"
            
        end if
        
    end if
end do
end subroutine



!--------- Load NMR data
subroutine loadNMR(infilename)
use defvar
use util
use NMRmod
implicit real*8 (a-h,o-z)
character c80tmp*80,c80tmp2*80
character(len=*) infilename

open(10,file=infilename,status="old")
call outputprog(10,iprog)

if (iprog==1) then !Gaussian
    call loclabel(10,"NAtoms=")
    read(10,*) c80tmp,ncenter
    if (.not.allocated(a)) allocate(a(ncenter))
    allocate(atmshd(ncenter))
    call loclabelfinal(10,"Magnetic shielding tensor (ppm)",ifound)
    if (ifound==0) then
        write(*,"(a)") " Error: Unable to find magnetic shielding tensors! Please check your Gaussian keywords"
        write(*,*) "Press ENTER button to exit program"
        read(*,*)
        stop
    end if
    read(10,*)
    do iatm=1,ncenter
        read(10,"(a)") c80tmp
        read(c80tmp,*) itmp,a(iatm)%name
        read(c80tmp(26:),*) atmshd(iatm)
        read(10,*);read(10,*);read(10,*);read(10,*)
    end do
else if (iprog==2) then !ORCA
	!Note that ORCA allows to only print shielding for user-specific atoms, so we have to load both atom index and value
	!Number of centers is determined when loading ORCA output file after Multiwfn just boots up
	if (ncenter==0) then
		call loclabel(10,"CARTESIAN COORDINATES (ANGSTROEM)",ifound)
		call skiplines(10,2)
		ncenter=0
		do while(.true.)
			read(10,"(a)") c80tmp
			if (c80tmp==" ") exit
			ncenter=ncenter+1
		end do
        allocate(a(ncenter))
        a%name=" "
    end if
    allocate(atmshd(ncenter))
    atmshd=0
    call loclabelfinal(10,"  Nucleus  Element",ifound)
    call skiplines(10,2)
    do while(.true.)
		read(10,"(a)") c80tmp
        if (c80tmp==" ") exit
        read(c80tmp,*) idx,c80tmp2,tmpshd
        iatm=idx+1
        a(iatm)%name=trim(c80tmp2)
        atmshd(iatm)=tmpshd
    end do
else if (iprog==0) then !Undetermined, viewed as plain text file
    ncenter=totlinenum(10,1)
    write(*,*) ncenter
    if (.not.allocated(a)) allocate(a(ncenter))
    allocate(atmshd(ncenter))
    do iatm=1,ncenter
        read(10,*) a(iatm)%name,atmshd(iatm)
    end do
else
    write(*,*) "Error: The program of generating this file is unsupported!"
    write(*,*) "Press ENTER button to exit program"
    read(*,*)
    stop
end if
close(10)

end subroutine