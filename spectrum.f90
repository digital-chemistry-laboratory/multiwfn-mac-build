!!!--------- Plot various kinds of spectra
!multiple.txt can records input file of multiple systems, but the system with the maximum number of transitions must be presented as the first entry
!If all(weight=1), then it is assumed that multiple.txt is not oriented for plotting weighted spectrum along with individual spectrum for each system, &
!but only for the latter, and in this case custom legend is allowed, which can be written in multiple.txt as the second column
!
!For single system, data is dataxall(1,1:numdata), nsystem=1, and numdataall(1)=numdata
!For n systems, nsystem=n, data of system i is dataxall(i,1:numdataall(i)), and in this case numdata is maximum length of numdataall
subroutine plotspectrum
use defvar
use dislin
use plot
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: weight(:) !Weight of various system for plotting mixed spectrum
real*8,allocatable :: dataxall(:,:),dataxall_org(:,:),strall(:,:),FWHMall(:,:) !Transition data loaded from multiple files. The first index corresponds to system index
integer,allocatable :: numdataall(:) !numdata is number of data, while for multiple system case
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
integer :: thk_curve=3,thk_weighted=8,thk_Yeq0=2,thk_discrete=1,thk_axis=1,thk_grid=1,thk_PVS=4,thk_OPVS=4 !thickness
integer :: ishowlabelleft=1,ishowlabelright=1 !If showing labels on left and right Y-axes
integer :: ndecimalX=-1,ndecimalYleft=-1,ndecimalYright=-1 !Number of decimal places in axes, use auto by default
integer :: height_axis=36,ticksize=36,legtextsize=36,labtype_Yleft=1,labtype_Yright=1,ilegendpos=7
integer :: legendx=400,legendy=160
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
real*8,allocatable :: PVScurve(:,:),OPVScurve(:) !Data point of PVS/PVDOS and OPVS/OPVDOS
real*8,allocatable :: PVScurveintra(:,:) !Data point of intrafragment part of PVS-I
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
		&.fch/fchk does not contain information needed for plotting spectrum"
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
write(*,*) " 1:IR  2:Raman (or pre-resonance Raman)  3:UV-Vis  4:ECD  5:VCD  6:ROA  7:NMR"
write(*,*) "-3: Directional UV-Vis"
write(*,*) " 0: Predicting color based on the UV-Vis spectrum recorded in text file"
read(*,*) ispectrum
if (ispectrum==-3) then
	ispectrum=3
    write(*,*) "Take which direction(s) into account?"
    write(*,*) "0 X+Y+Z (standard UV-Vis spectrum)"
    write(*,*) "1 X"
    write(*,*) "2 Y"
    write(*,*) "3 Z"
    write(*,*) "4 X+Y"
    write(*,*) "5 X+Z"
    write(*,*) "6 Y+Z"
    write(*,*) "7 Along specific direction by inputting a vector"
    read(*,*) iUVdir
    if (iUVdir==7) then
		write(*,*) "Input vector of oscillating electric field of incident light, e.g. 1,0,0.333"
        write(*,*) "Note: The vector will be automatically normalized"
        read(*,*) UVdirvec(:)
    end if
end if
if (ispectrum==7) then
    call NMRplot
    return
else if (ispectrum==1.or.ispectrum==2.or.ispectrum==5.or.ispectrum==6) then !IR, Raman, VCD, ROA
	ibroadfunc=1 !Use Lorentzian broadening
	iunitx=0 !cm^-1
else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
	ibroadfunc=2 !Use Gaussian broadening
	iunitx=2 !nm is default unit. But transition energies are loaded as eV
else if (ispectrum==0) then
	call gen_vis_curve
	call predict_color
	return
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
			if (isys==2.and.index(c200tmp,'\')/=0) then
				write(*,"(a)") " NOTE: ""\"" was found in the path, however in Linux/Mac, the separator in path should be ""/"", you need to manually fix it"
            end if
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
allocate(PVScurve(maxPVSfrag,num1Dpoints),OPVScurve(num1Dpoints),PVScurveintra(maxPVSfrag,num1Dpoints))
 

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
    if (ispectrum==3) write(*,*) "25 Evaluate color based on the spectrum in visible range"
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
        write(10,"('labtype_Yright',i5)") labtype_Yright
        write(10,"('num1Dpoints',i8)") num1Dpoints
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
        call readoption_int(10,"labtype_Yright",' ',labtype_Yright)
        call readoption_int(10,"num1Dpoints",' ',num1Dpoints)
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
			&this file can be directly used as input file of Multiwfn"
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
			&these files can be directly used as input file of Multiwfn"
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
			&is chosen as the unit, the curve will be generated in eV as X-axis first, and then convert to nm. Since current unit is nm, now you have to define the FWHM in eV."
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
				&If inputting 1.0, frequencies will correspond to the ones originally loaded from input file"
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
            if (labtype_Yright==1) write(*,"(a)") " -9 Set type of labels in right Y-axis, current: Float"
            if (labtype_Yright==2) write(*,"(a)") " -9 Set type of labels in right Y-axis, current: Exponent"
            if (labtype_Yright==3) write(*,"(a)") " -9 Set type of labels in right Y-axis, current: Scientific"
            if (nsystem>1.or.any(PVSnterm/=0)) then
                write(*,"(a,i3)") " 10 Set text size of legend, current:",legtextsize
                write(*,*) "11 Set position of legends"
            end if
            if (iYeq0==0) write(*,*) "12 Toggle drawing a line corresponding to Y=0, current: No"
            if (iYeq0==1) write(*,*) "12 Toggle drawing a line corresponding to Y=0, current: Yes"
            if (nsystem>1.and.ishowweicurve/=1) write(*,*) "13 Set color of line and curve of different systems"
            write(*,"(a,i8)") " 14 Set number of points in the curve, current:",num1Dpoints
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
            else if (isel2==-9) then
                write(*,*) "1 Float"
                write(*,*) "2 Exponent"
                write(*,*) "3 Scientific"
                read(*,*) labtype_Yright
            else if (isel2==10) then
                write(*,*) "Input size, e.g. 40"
                read(*,*) legtextsize
            else if (isel2==11) then
                write(*,*) "Choose position of legends"
                write(*,*) "0 Input coordinate of legend"
                write(*,*) "5 Lower left corner"
                write(*,*) "6 Lower right corner"
                write(*,*) "7 Upper right corner"
                write(*,*) "8 Upper left corner"
                read(*,*) ilegendpos
                if (ilegendpos==0) then
                    write(*,*) "Input X position of the legends with respect to upper left corner, e.g. 400"
					write(*,*) "The larger the value, the more the position of the legends is right"
					write(*,"(a,i6,a)") " If directly press ENTER button, current value",legendx," will be kept"
					read(*,"(a)") c80tmp
					if (c80tmp/=" ") read(c80tmp,*) legendx
					write(*,*) "Input Y position of the legends with respect to upper left corner, e.g. 160"
					write(*,*) "The larger the value, the lower the position of the legends"
					write(*,"(a,i6,a)") " If directly ENTER button, current value",legendy," will be kept"
					read(*,"(a)") c80tmp
					if (c80tmp/=" ") read(c80tmp,*) legendy
                end if
            else if (isel2==12) then
                if (iYeq0==1) then
                    iYeq0=0
                else
                    iYeq0=1
                end if
            else if (isel2==13) then
				do while(.true.)
					write(*,"(a)") " The index and current color of various systems are listed below, &
                    &now input index of a system to change its color, or input ""q"" to return"
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
            else if (isel2==14) then
				write(*,*) "Input number of points in the curve(s), e.g. 350"
                write(*,"(a)") " Note: Assume that the energy range for plotting is 380 to 760 nm, &
                &if you set number of points to (760-380)+1=381, then spacing between adjacent two points will be exactly 1 nm"
                read(*,*) num1Dpoints
                deallocate(curvex,curvey,curveytmp,curveyall,PVScurve,OPVScurve,PVScurveintra)
                allocate(curvex(num1Dpoints),curvey(num1Dpoints),curveytmp(num1Dpoints),curveyall(nsystem,num1Dpoints))
                allocate(PVScurve(maxPVSfrag,num1Dpoints),OPVScurve(num1Dpoints),PVScurveintra(maxPVSfrag,num1Dpoints))
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
                                if (OPVSidx1/=0.and.OPVSidx2/=0) OPVScomp(idata)=2*min(PVScomp(OPVSidx1,idata),PVScomp(OPVSidx2,idata))
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
                                        tmpall=0;tmpintra=0
										do icart=1,3*ncenter !Loop Cartesian coordinates
											iatm=ceiling(icart/3D0)
											icomp=icart-3*(iatm-1) !1/2/3=x,y,z of atom
                                            tmpval=dipder(icomp,iatm,idir)*normmat_atm(icomp,iatm,imode)
                                            tmpall=tmpall+tmpval
											if (cartlist(icart,ifrag)==1) then
												tmpintra=tmpintra+tmpval
												fragdipcomp(icomp,idir,ifrag,imode)=fragdipcomp(icomp,idir,ifrag,imode)+tmpval
                                            end if
                                        end do
                                        fragIRintra=fragIRintra+tmpintra**2 !Add intrafragment contribution in current direction
                                        fragIR=fragIR+tmpintra*tmpall !Add fragment contribution in current direction
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
										tmpintra1=0
                                        tmpintra2=0
										do icart=1,3*ncenter !Loop Cartesian coordinates
											iatm=ceiling(icart/3D0)
											icomp=icart-3*(iatm-1)
                                            tmpval=dipder(icomp,iatm,idir)*normmat_atm(icomp,iatm,imode)
											if (cartlist(icart,OPVSidx1)==1) then
												tmpintra1=tmpintra1+tmpval
											else if (cartlist(icart,OPVSidx2)==1) then
												tmpintra2=tmpintra2+tmpval
                                            end if
                                        end do
                                        coupIR=coupIR+2*tmpintra1*tmpintra2
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
                    thk_curve=4
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
                        &For example, inputting ""2-4,9 XY"" will only define X and Y components of atoms 2,3,4,9 as the fragment"
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
		end do !End menu loop of PVS
    
    end if
	
	if (isel==15.and.nsystem>1) then !Showing individual transition contribution is not possible when multiple files are involved
		write(*,*) "Error: This function is not available when multiple files are involved!"
		write(*,*) "Press ENTER button to continue"
		read(*,*)
		cycle
	end if
	if (isel==25) then !For predicting color, making range and point spacing in line with CIE1931 tristimulus functions
        xlow=360
		xhigh=830
		stepx=50
		iusersetX=1
		num1Dpoints=471
        deallocate(curvex,curvey,curveytmp,curveyall,PVScurve,OPVScurve,PVScurveintra)
        allocate(curvex(num1Dpoints),curvey(num1Dpoints),curveytmp(num1Dpoints),curveyall(nsystem,num1Dpoints))
        allocate(PVScurve(maxPVSfrag,num1Dpoints),OPVScurve(num1Dpoints),PVScurveintra(maxPVSfrag,num1Dpoints))
    end if
    
    
    

	!!=======================================================================!!
	!!=======================================================================!!
	!!=============== Below functions need calculation of curves ============!!
	!!=======================================================================!!
	!!=======================================================================!!
	if (isel==0.or.isel==1.or.isel==2.or.isel==15.or.isel==25) then
		!====== Construct correspondence array if outputting individual bands. Only available when one file is loaded
		!This function is not available when multiple systems are considered
		if (isel==15) then
			write(*,"(a)") " Input criterion of strength, e.g. 0.2, then the contribution curves of the transitions &
            &whose absolute strength larger than it will be exported to plain text files"
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
            &absolute value of a transition by sum of absolute values of all transitions"
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
                    if (ispecial==0) then
						write(10,"(f13.5,1PE18.8E3)") curvex(ipt),curvey(ipt)
                    else !Format is compatible with requirement of input file of CIE1931xy.V.1.6.0.2.exe
						write(c200tmp,"(i8)") nint(curvex(ipt))
						write(c200tmp2,"(f16.3)") curvey(ipt)
                        write(10,"(a)") trim(adjustl(c200tmp))//' '//trim(adjustl(c200tmp2))
                    end if
				end do
				close(10)
				write(*,*) "Curve data has been exported to spectrum_curve.txt in current folder"
			else !Also output curve for all systems
				if (any(weight/=1)) then !Output weighted spectrum 
					open(10,file="spectrum_curve.txt",status="replace")
					do ipt=1,num1Dpoints
                        write(10,"(f13.5,1PE18.8E3)") curvex(ipt),curvey(ipt)
					end do
					close(10)
					write(*,"(a)") " The curve data corresponding to weighted spectrum has been exported to spectrum_curve.txt in current folder"
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
            &in current folder, data of each system is separated by a blank line"
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
	if (isel==0.or.isel==1.or.isel==25) then

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
        if (ilegendpos==0) then
			call legpos(legendx,legendy) !Absolute position of legends
            call legend(clegend,3)
        else
			if (nsystem>1.and.ishowweicurve/=2) call legend(clegend,ilegendpos)
			if (any(PVSnterm/=0)) call legend(clegend,ilegendpos)
        end if
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
				call namdis(60,'Y') !Use larger distance between name and axis
			else
				if (labtype_Yright==1) call labels("FLOAT","Y")
				if (labtype_Yright==2) call labels("XEXP","Y")
				if (labtype_Yright==3) call labels("FEXP","Y")
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
			call height(ticksize) !Size of ticks
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
			call height(36) !Recover to default
		    CALL ENDGRF
			call color("WHITE")
            if (iexportlevel==1) write(*,"(a)") " Line data corresponding to various spike sets have been &
            &exported to spike[index].txt in current folder"
        end if
        
		call disfin
		if (isavepic==1) write(*,*) "Graphical file has been saved to current folder with ""dislin"" prefix"
        
        !Show color based on spectrum
        if (isel==25) call predict_color
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
character ctest,ctest2,ctest3,c3tmp*3,c80tmp*80,c200tmp*200
integer ispectrum,imode
integer :: nrdfreq=0 ! >0 means pre-resonance raman, which loads external field frequency
real*8 tmpvec(3)
real*8,allocatable :: rdfreq(:),tmparr(:)

if (allocated(datax)) deallocate(datax,str,FWHM)
open(10,file=loadspecname,status="old")

!Check if is sTDA output file
if (index(loadspecname,".dat")/=0) then
	if (imode==0) write(*,*) "Recognized as sTDA program output file"
    if (ispectrum==3.and.iUVdir/=0) then
		write(*,*) "Error: Plotting directional UV-Vis spectrum does not support sTDA code"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        stop
    end if
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
if (igauout==0) call loclabel(10,"Entering Gaussian System",igauout,maxline=200)
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
			if (ispectrum==1) then
				call loclabel(10,"Anharmonic Infrared Spectroscopy",ifound,0)
			else if (ispectrum==2) then
				call loclabel(10,"Anharmonic Raman Spectroscopy",ifound,0)
			else if (ispectrum==5) then
				call loclabel(10,"Anharmonic VCD Spectroscopy",ifound,0)
            end if
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
			        if (ispectrum==1) then
						call loclabel(10,"Anharmonic Infrared Spectroscopy")
			        else if (ispectrum==2) then
						call loclabel(10,"Anharmonic Raman Spectroscopy")
			        else if (ispectrum==5) then
						call loclabel(10,"Anharmonic VCD Spectroscopy")
                    end if
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
        if (iUVdir/=0) then !Load transition dipole moments and calculate Cartesian component of oscillator strength
			rewind(10)
			do i=1,numopt
				call loclabel(10,"Ground to excited state transition electric dipole moments",ifound,0)
				read(10,*)
			end do
			read(10,*)
            do i=1,numdata
				read(10,*) itmp,xdip,ydip,zdip
                call calc_trdipsqr(xdip,ydip,zdip,tmp)
                str(i)=2D0/3D0*tmp*datax(i)/au2eV
            end do
        end if
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
		if (ispectrum==1.or.ispectrum==2.or.ispectrum==5) then !IR, Raman
			if (ispectrum==1) call loclabel(10,"IR SPECTRUM")
			if (ispectrum==2) call loclabel(10,"RAMAN SPECTRUM")
			if (ispectrum==5) call loclabel(10,"VCD-Intensity")
			call loclabel(10," Mode    freq (cm**-1)",ioldformat,0)
			if (ioldformat==1) then !ORCA 4
				read(10,*);read(10,*)
			else !ORCA >=5
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
			if (ispectrum==1) call loclabelfinal(10,"IR SPECTRUM",nfound)
			if (ispectrum==2) call loclabelfinal(10,"RAMAN SPECTRUM",nfound)
			if (ispectrum==5) call loclabelfinal(10,"VCD-Intensity",nfound)
			if (ioldformat==1) then !ORCA 4
                call loclabel(10," Mode    freq (cm**-1)",ifound,0)
				read(10,*);read(10,*)
			else !ORCA >=5
				if (ispectrum==1.or.ispectrum==2) call loclabel(10," Mode   freq  ",ifound,0)
				read(10,*);read(10,*);read(10,*)
			end if
			do i=1,numdata
				read(10,*) c80tmp,datax(i),str(i)
                write(*,*) c80tmp,datax(i),str(i)
			end do
			FWHM=8D0
		else if (ispectrum==3.or.ispectrum==4) then !UV-Vis, ECD
			call loclabelfinal(10,"Number of roots to be determined",nfound)
            if (nfound>=1) then !TDDFT/TDA-DFT
				read(10,"(50x,i7)") numdata
				call loclabel(10,"Entering triplet calculation",itriplet)
                if (itriplet==1) numdata=numdata*2
				if (imode==1) then !Have obtained number of data, return
					close(10)
					return
				end if
				allocate(datax(numdata),str(numdata),FWHM(numdata))
				if (ispectrum==3) then
                    call loclabelfinal(10,"  ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",nfound) !Two spaces before the label make sure it will not locate to SOC corrected part
                else if (ispectrum==4) then
                    call loclabelfinal(10,"  CD SPECTRUM",nfound) !ORCA doesn't distinguish rotatory strength forms before 6.0, since ORCA 6, this way will load velocity form, which is good
                end if
                if (nfound>1) write(*,"(a)") " Note: Excited state information was outputted by ORCA multiple times, the finally outputted ones will be loaded"
				call skiplines(10,5)
                read(10,"(a)") c80tmp
                if (index(c80tmp,'->')/=0) then
					iORCAver=6 !>=6.0.x
                else
					iORCAver=5 !<=5.0.x
                end if
                backspace(10)
                datax=0;str=0
				do i=1,numdata
					if (ispectrum==3) then !t1: wavenumber, t2: f, t3: D2
						if (iORCAver==5) then
							read(10,*,iostat=ierror) itmp,t1,rnouse,t2,t3,xdip,ydip,zdip
                        else if (iORCAver==6) then
							itmp=i
							read(10,*,iostat=ierror) c80tmp,c80tmp,c80tmp,c80tmp,t1,rnouse,t2,t3,xdip,ydip,zdip
                        end if
					else if (ispectrum==4) then !t1: wavenumber, t2: R
						if (iORCAver==5) then
							read(10,*,iostat=ierror) itmp,t1,rnouse,t2
                        else if (iORCAver==6) then
							itmp=i
							read(10,*,iostat=ierror) c80tmp,c80tmp,c80tmp,c80tmp,t1,rnouse,t2
                        end if
                    end if
                    if (ierror/=0) exit !For SF-TDDFT, number of states recorded in this field is less than nstates by 1, because one of SF-TDDFT states is viewed as ground state
                    datax(itmp)=t1
                    str(itmp)=t2
					if (iUVdir/=0.and.ispectrum==3) then
						xdipsq=xdip**2
						ydipsq=ydip**2
						zdipsq=zdip**2
						if (iUVdir==1) then
							tmp=xdipsq
						else if (iUVdir==2) then
							tmp=ydipsq
						else if (iUVdir==3) then
							tmp=zdipsq
						else if (iUVdir==4) then
							tmp=xdipsq+ydipsq
						else if (iUVdir==5) then
							tmp=xdipsq+zdipsq
						else if (iUVdir==6) then
							tmp=ydipsq+zdipsq
						else if (iUVdir==7) then
							tmpvec=UVdirvec/dsqrt(sum(UVdirvec**2))
							tmp=tmpvec(1)**2*xdipsq+tmpvec(2)**2*ydipsq+tmpvec(3)**2*zdipsq
						end if
						str(i)=2D0/3D0*tmp*datax(i)/au2cm
					end if
				end do
				call loclabel(10,"SOC CORRECTED ABSORPTION",ifound,0)
				if (ifound==1) then
					numdata=numdata/2 !Because of itriplet==1, numdata has been doubled, now recover it to number of singlet states
					write(*,"(a)") " Spin-orbit coupling corrected spectra information was found, &
					&would you like to plot this kind of spectrum instead of the one without correction? (y/n)"
					read(*,*) ctest
					if (ctest=='y'.or.ctest=='Y') then
						numdata=4*numdata !If root=n, then there will be n singlet states and 3n triplet sublevels
						deallocate(datax,str,FWHM)
						allocate(datax(numdata),str(numdata),FWHM(numdata))
                        if (ispectrum==4) call loclabelfinal(10,"SOC CORRECTED CD SPECTRUM",nfound)
						call skiplines(10,5)
						do i=1,numdata
							if (iORCAver==5) then
								read(10,*) tmp1,tmp2,tmp3,tmp4,tmp5 !Since 5.0, the first column is always 0, strange!
								if (tmp1==0) then !ORCA >=5.0
									datax(i)=tmp3
									str(i)=tmp5
								else !ORCA 4.x
									datax(i)=tmp2
									str(i)=tmp4
								end if
                            else if (iORCAver==6) then
								if (ispectrum==3) then
									read(10,*) c80tmp,c80tmp,c80tmp,c80tmp,datax(i),c80tmp,str(i),D2sqr,xdipnorm,ydipnorm,zdipnorm
                                else if (ispectrum==4) then
									read(10,*) c80tmp,c80tmp,c80tmp,c80tmp,datax(i),c80tmp,str(i)
                                end if
                            end if
							if (iUVdir/=0.and.ispectrum==3) then
								if (iORCAver==6) then
									xdipsq=xdipnorm**2
									ydipsq=ydipnorm**2
									zdipsq=zdipnorm**2
                                else
									backspace(10)
									read(10,"(a)") c200tmp
									i1=strcharpos(c200tmp,'(',1)+1
									i2=strcharpos(c200tmp,')',1)-1
									read(c200tmp(i1:i2),*) xdipr,xdipi !Real and imaginary parts of transition dipole moment in X
									xdipsq=xdipr**2+xdipi**2
									i1=strcharpos(c200tmp,'(',2)+1
									i2=strcharpos(c200tmp,')',2)-1
									read(c200tmp(i1:i2),*) ydipr,ydipi !Real and imaginary parts of transition dipole moment in Y
									ydipsq=ydipr**2+ydipi**2
									i1=strcharpos(c200tmp,'(',3)+1
									i2=strcharpos(c200tmp,')',3)-1
									read(c200tmp(i1:i2),*) zdipr,zdipi !Real and imaginary parts of transition dipole moment in Z
									zdipsq=zdipr**2+zdipi**2
                                end if
								if (iUVdir==1) then
									tmp=xdipsq
								else if (iUVdir==2) then
									tmp=ydipsq
								else if (iUVdir==3) then
									tmp=zdipsq
								else if (iUVdir==4) then
									tmp=xdipsq+ydipsq
								else if (iUVdir==5) then
									tmp=xdipsq+zdipsq
								else if (iUVdir==6) then
									tmp=ydipsq+zdipsq
								else if (iUVdir==7) then
									tmpvec=UVdirvec/dsqrt(sum(UVdirvec**2))
									tmp=tmpvec(1)**2*xdipsq+tmpvec(2)**2*ydipsq+tmpvec(3)**2*zdipsq
								end if
								str(i)=2D0/3D0*tmp*datax(i)/au2cm
							end if
						end do
                    else
						if (iORCAver==6) then !>=6.0.x. Because in this case ORCA also outputs S to T transitions and mixed with S to S together, we need to load all of them, 2*numdata
							deallocate(datax,str,FWHM)
                            numdata=numdata*2
							allocate(datax(numdata),str(numdata),FWHM(numdata))
							if (ispectrum==3) then
								call loclabelfinal(10,"  ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",nfound)
							else if (ispectrum==4) then
								call loclabelfinal(10,"  CD SPECTRUM",nfound)
							end if
							call skiplines(10,5)
							datax=0;str=0
							do i=1,numdata
								if (ispectrum==3) then !t1: wavenumber, t2: f, t3: D2
									read(10,*,iostat=ierror) c80tmp,c80tmp,c80tmp,c80tmp,t1,rnouse,t2,t3,xdip,ydip,zdip
								else if (ispectrum==4) then !t1: wavenumber, t2: R
									read(10,*,iostat=ierror) c80tmp,c80tmp,c80tmp,c80tmp,t1,rnouse,t2
								end if
								datax(i)=t1
								str(i)=t2
							end do
                        end if
					end if
				end if
            else !Should be (DLPNO-)(ST)EOM-CCSD or others
				!Count how many data are there
				call loclabelfinal(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",nfound)
                if (nfound>1) write(*,"(a)") " Note: Excited state information was outputted by ORCA multiple times, the finally outputted ones will be loaded"
                call skiplines(10,5)
                read(10,"(a)") c80tmp
                if (index(c80tmp,'->')/=0) then
					iORCAver=6 !>=6.0.x
                else
					iORCAver=5 !<=5.0.x
                end if
                backspace(10)
                numdata=0
                do while(.true.)
					read(10,"(a)") c80tmp
                    if (c80tmp==" ".or.c80tmp(1:1)=='*') exit !Since 6.0, final line will be "*The positivity of oscillator strengths is ..."
					numdata=numdata+1
                end do
				allocate(datax(numdata),str(numdata),FWHM(numdata))
				if (ispectrum==3) then
					call loclabelfinal(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",nfound) !Since 6.0, the finally printed are the case of LEFT-RIGHT TRANSITION MOMENTS
				else if (ispectrum==4) then
					call loclabel(10,"CD SPECTRUM")
                end if
				call skiplines(10,5)
				do i=1,numdata
					if (iORCAver==5) then
						read(10,*) rnouse,datax(i),tmpval,str(i),tmpval,xdip,ydip,zdip
                    else
						read(10,*) c80tmp,c80tmp,c80tmp,c80tmp,datax(i),tmpval,str(i),tmpval,xdip,ydip,zdip
                    end if
                    if (tmpval>datax(i)) then !I noticed ORCA at least 5.0.1~5.0.3 has a bug, in output of CD SPECTRUM, the values in cm-1 and nm are inversed
						ttt=tmpval
                        datax(i)=tmpval
                        tmpval=ttt
                    end if
                    if (iUVdir/=0.and.ispectrum==3) then
						call calc_trdipsqr(xdip,ydip,zdip,tmp)
						str(i)=2D0/3D0*tmp*datax(i)/au2cm
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
		if (ispectrum==3) then
			FWHM=2D0/3D0
		else if (ispectrum==4) then
			FWHM=0.2D0
        end if
		if (ispectrum==3.and.iUVdir/=0) then
			call loclabel(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",ifound,0)
			call skiplines(10,5)
			do i=1,numdata
				read(10,*) rnouse,tmpene,tmpval,tmpval,tmpval,xdip,ydip,zdip
                call calc_trdipsqr(xdip,ydip,zdip,tmp)
				str(i)=2D0/3D0*tmp*tmpene/au2cm
			end do
        end if
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
		iSOC=0
        ctest='n'
        call loclabelfinal(10,"SOC-corrected exc.",iSOC)
        if (iSOC==1) then
			write(*,"(a)") " Spin-orbit coupling corrected spectra information was found, &
			&would you like to plot this kind of spectrum instead of the one without correction? (y/n)"
			read(*,*) ctest
			if (ctest=='n') then
				write(*,*) "Information of singlet excited states will be loaded"
				call loclabelfinal(10,"states of multiplicity 1",ifound)
				call loclabel(10,"number   energy (eV)",ifound,0)
			end if
        else
            call loclabelfinal(10,"number   energy (eV)",ifound)
			if (ifound==0) then
				write(*,*) "Error: Unable to find electronic excitation information!"
				write(*,*) "Press ENTER button to return"
				read(*,*)
				return
			end if
        end if
        read(10,*)
        read(10,*)
		numdata=0
		do while(.true.)
			read(10,"(a)") c80tmp
            if (c80tmp==" ".or.index(c80tmp,"SOC")/=0) exit
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
            if (iSOC==0.or.(iSOC==1.and.ctest=='n')) then
				read(10,*) c80tmp,c80tmp,datax(i),xdip,ydip,zdip,str(i)
				if (iUVdir/=0) then
					call calc_trdipsqr(xdip,ydip,zdip,tmp)
					str(i)=2D0/3D0*tmp*datax(i)/au2eV
				end if
            else
				read(10,*) c80tmp,datax(i),str(i)
            end if
		end do
    
    else if (ispectrum==1.or.ispectrum==2) then !IR and Raman
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
            !Read frequencies
            read(10,"(a)") c80tmp
			if (iread==1) then
				read(c80tmp(23:),*,iostat=ierror) datax(inow)
			else if (iread==2) then
				read(c80tmp(23:),*,iostat=ierror) datax(inow),datax(inow+1)
			else if (iread==3) then
				read(c80tmp(23:),*,iostat=ierror) datax(inow),datax(inow+1),datax(inow+2)
            end if
            if (ierror/=0) then
				write(*,*) "Error encountered while loading the following line:"
                write(*,"(a)") trim(c80tmp)
                write(*,*) "Press ENTER button to exit"
                read(*,*)
                stop
            end if
            !Read IR/Raman intensities
			if (ispectrum==2) read(10,*) !Skip IR line
			do itmp=1,iread
				read(10,"(a)") c80tmp
				if (itmp==1) then
					idxmode=inow
					read(c80tmp(22:),*,iostat=ierror) str(idxmode)
                else if (itmp==2) then
					idxmode=inow+1
					read(c80tmp(36:),*,iostat=ierror) str(idxmode)
                else if (itmp==3) then
					idxmode=inow+2
					read(c80tmp(57:),*,iostat=ierror) str(idxmode)
                end if
                if (ierror/=0) then
					if (ispectrum==1) write(*,"(a,i6,a)") " Error: Cannot properly load IR intensity of mode",idxmode,"!"
					if (ispectrum==2) write(*,"(a,i6,a)") " Error: Cannot properly load Raman activity of mode",idxmode,"!"
                    write(*,"(a)") " Please check corresponding data in the input file, perhaps the value is too large, &
                    &making the data be recorded as ************"
                    if (ispectrum==1) then
						write(*,*) "If it is the case, please manually input the IR intensity in kM/mol, e.g. 89.64"
                        write(*,"(a)") "Hint: If you request CP2K to export .mol (Molden) file containing vibrational information, &
                        &you can find the intensity in [INT] field. The value to be inputted here should be that in [INT] multiplied by 974.85541. &
                        &Alternatively, you can directly use the .mol file as input file" 
                    else if (ispectrum==2) then
						write(*,*) "If it is the case, please manually input the Raman activity, e.g. 89.64"
                    end if
                    read(*,*) str(idxmode)
					backspace(10)
                end if
            end do
			if (ilackdata<=3) exit
			ilackdata=ilackdata-3
			inow=inow+3
		end do
    end if
	close(10)
    return
end if

!Check if is BDF output file. This part of code was contributed by Cong Wang, 2022-Dec-1
call loclabel(10,"BDF",iBDFout,maxline=100)
rewind(10)
if (iBDFout==1) then
    if (imode==0) write(*,*) "Recognized as a BDF output file"
    !UV-Vis
    if (ispectrum==3) then
        call loclabel(10,"[algorithm]",ifound,0)
        call skiplines(10,2)
        numdata=0
        do while(.true.)
            read(10,"(a)") c200tmp
            if (index(c200tmp,"[dvdson_parameters]")/=0) exit
            read(c200tmp,*) itmp1,c3tmp,itmp2,itmp3
            numdata=numdata+itmp2
        end do
        if (imode==1) then !Have obtained number of data, return
            close(10)
            return
        end if
        allocate(datax(numdata),str(numdata),FWHM(numdata)) !Transition energy, strength and FWHM loaded from only one file
        call loclabel(10,"List of excitations",ifound,0)
        if (ifound==1) then
            do while(.true.)
                call loclabel(10,"List of excitations",ifound,0)
                if (ifound==1) then
                    call loclabel(10,"isf=",isfout,0,10)
                    if (isfout==1) then
                        read(10,"(a)") c200tmp
                        read(c200tmp(index(c200tmp,"=")+1:),*) isf
                        if (isf==0) exit
                    end if
                end if
            end do
        end if
        call loclabel(10,"No. Pair   ExSym   ExEnergies     Wavelengths      f",ifound,0,10)
        call skiplines(10,2)
        datax=0;str=0
        do iec=1,numdata
            read(10,"(a)") c200tmp
            read(c200tmp(index(c200tmp,"eV")-10:index(c200tmp,"eV")-1),*) datax(iec)
            read(c200tmp(index(c200tmp,"nm")+2:index(c200tmp,"nm")+10),*) str(iec)
        end do
        if (iUVdir/=0.and.ispectrum==3) then
            call loclabel(10,"Ground to excited state Transition electric dipole moments",ifound,0,10)
            call skiplines(10,2)
            !write(*,*) "TAG iUVdir: ", iUVdir
            do iec=1,numdata
                read(10,*) istate,xdip,ydip,zdip
                call calc_trdipsqr(xdip,ydip,zdip,tmp)
                str(iec)=2D0/3D0*tmp*datax(iec)/au2eV
            end do
        end if
        !write(*,*) "TAG -str: ", str
        call loclabel(10,"List of SOC-SI results",ifound,0)
        if (ifound==1) then
            write(*,"(a)") " Spin-orbit coupling corrected spectra information was found, &
            &would you like to plot this kind of spectrum instead of the one without correction? (y/n)"
            read(*,*) ctest
            if (ctest=='y'.or.ctest=='Y') then
                numdata=4*numdata !If root=n, then there will be n singlet states and 3n triplet sublevels
                deallocate(datax,str,FWHM)
                allocate(datax(numdata),str(numdata),FWHM(numdata))
                do isoc=1,numdata
                    call loclabel(10,"E_J-E_I         fosc",ifound,0)
                    call skiplines(10,2)
                    read(10,*) tmp1,tmp2,tmp3,tmp4,datax(isoc),str(isoc)
                    if (iUVdir/=0.and.ispectrum==3) then
                        call loclabel(10,"Norm=",ifound,0)
                        read(10,"(a)") c200tmp
                        read(c200tmp(index(c200tmp,"=")+1:),*) xdip,ydip,zdip
						call calc_trdipsqr(xdip,ydip,zdip,tmp)
                        str(isoc)=2D0/3D0*tmp*datax(isoc)/au2eV
                    end if
                end do
            end if
        end if
        FWHM=2D0/3D0
    end if
    close(10)
    return
end if

!Check if is Molden file containing vibrational frequencies and IR intensities
call loclabel(10,"[Molden Format]",imolden,maxline=100)
rewind(10)
if (imolden==1) then
    if (imode==0) write(*,*) "Recognized as a Molden input file"
    if (ispectrum/=1) then
		write(*,*) "Error: Only IR spectrum is supported when Molden file is used"
        close(10)
        return
    end if
    call loclabel(10,"[FREQ]")
    read(10,*)
    numdata=0
    do while(.true.)
		read(10,"(a)",iostat=ierror) c200tmp
        if (index(c200tmp,'[')/=0.or.c200tmp==" ".or.ierror/=0) exit
        numdata=numdata+1
    end do
    if (imode==1) then !Have obtained number of data, return
        close(10)
        return
    end if
	allocate(datax(numdata),str(numdata),FWHM(numdata))
	FWHM=8D0
    call loclabel(10,"[FREQ]")
    read(10,*)
	do i=1,numdata
		read(10,*) datax(i)
	end do
    call loclabel(10,"[INT]")
    read(10,*)
	do i=1,numdata
		read(10,*) str(i)
	end do
    str=str*974.85541D0
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




!!-------- Calculate finally used square of transition electric dipole moment
!Global variables: iUVdir=Type of direction(s) to consider, UVdirvec=The only direction to consider when iUVdir=7
!xdip,ydip,zdip: X/Y/Z component of transition electric dipole moment
!trdipsqr: Final result
subroutine calc_trdipsqr(xdip,ydip,zdip,trdipsqr)
use defvar
real*8 xdip,ydip,zdip,trdipsqr,UVdirvectmp(3)
if (iUVdir==1) then
	trdipsqr=xdip**2
else if (iUVdir==2) then
	trdipsqr=ydip**2
else if (iUVdir==3) then
	trdipsqr=zdip**2
else if (iUVdir==4) then
	trdipsqr=xdip**2+ydip**2
else if (iUVdir==5) then
	trdipsqr=xdip**2+zdip**2
else if (iUVdir==6) then
	trdipsqr=ydip**2+zdip**2
else if (iUVdir==7) then
    UVdirvectmp=UVdirvec/dsqrt(sum(UVdirvec**2))
    trdipsqr=(UVdirvectmp(1)*xdip)**2+(UVdirvectmp(2)*ydip)**2+(UVdirvectmp(3)*zdip)**2
end if
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
use dislin
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
			if (isys==2.and.index(c200tmp,'\')/=0) then
				write(*,"(a)") " NOTE: ""\"" was found in the path, however in Linux/Mac, the separator in path should be ""/"", you need to manually fix it"
            end if
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
allocate(shdatm(200,ncenter,nsystem)) !Degeneracy is assumed to be at most 200
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
if (allocated(curvey)) deallocate(curvey) !Global array
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
            &a~f were all calculated under chloroform represented by SMD model, geometries were optimized at B3LYP/def2-SVP in vaccum"
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
                    &properly select the element to be considered"
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
                    &properly select the element to be considered"
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
else if (iprog==7) then !BDF. This part of code was contributed by Cong Wang, 2022-Dec-1
    call loclabel(10,"Number of atoms",ifound)
    read(10,"(a)") c80tmp
    read(c80tmp(index(c80tmp,":")+1:),*) ncenter
    if (.not.allocated(a)) allocate(a(ncenter))
    allocate(atmshd(ncenter))
    rewind(10)
    do while(.true.)
        call loclabel(10,"Nuclear Magnetic shielding result in PPM",ifound,0)
        if (ifound==0) then
            write(*,"(a)") " Error: Unable to find magnetic shielding tensors!"
            write(*,*) "Press ENTER button to exit program"
            read(*,*)
            stop
        end if
        read(10,"(a)") c80tmp
        if(index(c80tmp,"giao")/=0) then
            backspace(10)
            exit
        end if
    end do
    do iatm=1,ncenter
        call loclabel(10,"NMR shielding tensor and constant of nucleus",irewind=0)
        read(10,"(a)") c80tmp
        read(c80tmp(index(c80tmp,"atom")+4:),*) a(iatm)%name
    end do
    call loclabel(10,"Isotropic/anisotropic constant by atom order",irewind=0)
    read(10,*)
    do iatm=1,ncenter
        read(10,*) atmshd(iatm), tmp
    end do
else if (iprog==0) then !Undetermined
	call loclabel(10,"Shielding atom at atomic positions",ifound)
    if (ifound==1) then !CP2K .data file of NMR
		call loclabelfinal(10,"ISOTROPY =",ncenter)
		if (.not.allocated(a)) allocate(a(ncenter))
		allocate(atmshd(ncenter))
        rewind(10)
        read(10,*)
        do iatm=1,ncenter
			read(10,*) c80tmp,a(iatm)%name
            call skiplines(10,12)
            read(10,*) c80tmp,c80tmp,atmshd(iatm)
        end do
    else !Viewed as a plain text file
		ncenter=totlinenum(10,1)
		if (.not.allocated(a)) allocate(a(ncenter))
		allocate(atmshd(ncenter))
		do iatm=1,ncenter
			read(10,*) a(iatm)%name,atmshd(iatm)
		end do
    end if
end if
close(10)

end subroutine




!!---- Load spectrum curve X-Y data from input file, and generate curve from 360 to 830 with 1 as spacing by interpolation
!The new curve is stored in global arrays curvex and curvey
subroutine gen_vis_curve
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: rawx(:),rawy(:)

if (allocated(curvex)) deallocate(curvex,curvey)
num1Dpoints=471
allocate(curvex(num1Dpoints),curvey(num1Dpoints))

write(*,*) "Loading data points..."
open(10,file=filename,status="old")
nraw=totlinenum(10,1)
allocate(rawx(nraw),rawy(nraw))
do idata=1,nraw
    read(10,*,iostat=ierror) rawx(idata),rawy(idata)
    if (ierror/=0) then
		write(*,"(a)") " Error: Unable to load data! The input file should contain two columns recording wavelengths and absorption strengths, respectively"
        write(*,*) "Press ENTER button to exit program"
        read(*,*)
        stop
    end if
    if (idata>1) then
	    if (rawx(idata)<=rawx(idata-1)) then
			write(*,"(a)") " Error: The X data in the input file must be ordered from small to large! Please manually fix your input file"
            write(*,*) "Press ENTER button to exit Multiwfn"
            read(*,*)
            stop
        end if
    end if
end do
write(*,"(' Number of data points:',i8)") nraw
write(*,"(' Range of X data:',f8.2,' to',f8.2)") minval(rawx),maxval(rawx)
close(10)

!Use Lagrangian interpolation to generate needed curve
do ipt=1,num1Dpoints
    inm=360+ipt-1
    curvex(ipt)=inm
    call lagintpol(rawx(:),rawy(:),nraw,dfloat(inm),curvey(ipt),der1,der2,1)
    !call linintpol(rawx(:),rawy(:),nraw,dfloat(inm),curvey(ipt))
    !write(11,*) curvex(ipt),curvey(ipt)
end do
end subroutine




!!---------- Predict color based on spectrum curve (curvex and curvey, from 360 to 830 nm with 1 nm spacing)
subroutine predict_color
use defvar
use plot
implicit real*8 (a-h,o-z)
character selectyn
real*8 CIE_Xfunc(830),CIE_Yfunc(830),CIE_Zfunc(830)

call gen_CIE_tristimulus_func(CIE_Xfunc,CIE_Yfunc,CIE_Zfunc)
CIE_X=0
CIE_Y=0
CIE_Z=0
do ipoint=1,num1Dpoints !In the present case, the X of spectrum curve starts from 360 and ends to 830 with spacing of 1 nm, which in line with tristimulus functions
	inm=ipoint-1+360
	CIE_X=CIE_X+curvey(ipoint)*CIE_Xfunc(inm)
	CIE_Y=CIE_Y+curvey(ipoint)*CIE_Yfunc(inm)
	CIE_Z=CIE_Z+curvey(ipoint)*CIE_Zfunc(inm)
end do
write(*,"(/,' CIE1931 XYZ:       ',3f18.6)") CIE_X,CIE_Y,CIE_Z
tmp=max(max(CIE_X,CIE_Y),CIE_Z)
if (tmp<1D-8) then
	write(*,"(a,1PE16.8,a)") " Error: Maximum component is merely",tmp,", it is colorless and thus predicting color is meaningless!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if
sCIE_X=CIE_X/tmp
sCIE_Y=CIE_Y/tmp
sCIE_Z=CIE_Z/tmp
write(*,"(' Fractional CIE1931 XYZ:',3f18.10)") sCIE_X,sCIE_Y,sCIE_Z
CIE_smallx=CIE_X/(CIE_X+CIE_Y+CIE_Z)
CIE_smally=CIE_Y/(CIE_X+CIE_Y+CIE_Z)
write(*,"(' CIE1931 xy:        ',2f18.10)") CIE_smallx,CIE_smally

!I also tried to take complementary color for X,Y,Z and then convert to RGB, in principle this is better. &
!But my test found this doesn't bring notable advantage compared to first convering to RGB and then take complementary color, the result is also not worsen
!sCIE_X=1-sCIE_X
!sCIE_Y=1-sCIE_Y
!sCIE_Z=1-sCIE_Z

Rcomp =  3.2404542D0*sCIE_X  -1.5371385D0*sCIE_Y  -0.4985314D0*sCIE_Z
Gcomp = -0.9692660D0*sCIE_X + 1.8760108D0*sCIE_Y + 0.0415560D0*sCIE_Z
Bcomp =  0.0556434D0*sCIE_X  -0.2040259D0*sCIE_Y + 1.0572252D0*sCIE_Z
write(*,"(a)") " Note the R,G,B values show below correspond to standard RGB (sRGB) color space"
write(*,"(' RGB (0-1):  ',3f10.6)") Rcomp,Gcomp,Bcomp
intRcomp=nint(Rcomp*255);intGcomp=nint(Gcomp*255);intBcomp=nint(Bcomp*255)
write(*,"(' RGB (0-255):',3i6)") intRcomp,intGcomp,intBcomp

!Gamma correction. I think this is redundant and often makes result worse
!if (Rcomp<=0.0031308D0) then
!	Rcomp=12.92D0*Rcomp
!else
!	Rcomp=1.055D0*Rcomp**(1D0/2.4D0)-0.055D0
!end if
!if (Gcomp<=0.0031308D0) then
!	Gcomp=12.92D0*Gcomp
!else
!	Gcomp=1.055D0*Gcomp**(1D0/2.4D0)-0.055D0
!end if
!if (Bcomp<=0.0031308D0) then
!	Bcomp=12.92D0*Bcomp
!else
!	Bcomp=1.055D0*Bcomp**(1D0/2.4D0)-0.055D0
!end if
        
if (Rcomp<0.or.Rcomp>1.or.Gcomp<0.or.Gcomp>1.or.Bcomp<0.or.Bcomp>1) then
	write(*,"(a)") " Note: The color exceeds sRGB color space! Now the R,G,B values are scaled into valid range:"
	if (Rcomp<0) Rcomp=0;if (Gcomp<0) Gcomp=0;if (Bcomp<0) Bcomp=0
	if (Rcomp>1) Rcomp=1;if (Gcomp>1) Gcomp=1;if (Bcomp>1) Bcomp=1
	write(*,"(' RGB (0-1):  ',3f10.6)") Rcomp,Gcomp,Bcomp
	intRcomp=nint(Rcomp*255);intGcomp=nint(Gcomp*255);intBcomp=nint(Bcomp*255)
	write(*,"(' RGB (0-255):',3i6)") intRcomp,intGcomp,intBcomp
end if
Rcomp_op=1-Rcomp !Opposite color
Gcomp_op=1-Gcomp
Bcomp_op=1-Bcomp
write(*,"(' RGB of complementary color:',3i6)") nint(Rcomp_op*255),nint(Gcomp_op*255),nint(Bcomp_op*255)
shiftval=1-max(max(Rcomp,Gcomp),Bcomp)
Rcomp_maxbr=Rcomp+shiftval
Gcomp_maxbr=Gcomp+shiftval
Bcomp_maxbr=Bcomp+shiftval
shiftval=1-max(max(Rcomp_op,Gcomp_op),Bcomp_op)
Rcomp_op_maxbr=Rcomp_op+shiftval
Gcomp_op_maxbr=Gcomp_op+shiftval
Bcomp_op_maxbr=Bcomp_op+shiftval
write(*,"(' RGB of original color (maximum brightness):     ',3i6)") nint(Rcomp_maxbr*255),nint(Gcomp_maxbr*255),nint(Bcomp_maxbr*255)
write(*,"(' RGB of complementary color (maximum brightness):',3i6)") nint(Rcomp_op_maxbr*255),nint(Gcomp_op_maxbr*255),nint(Bcomp_op_maxbr*255)
!Show color on screen
do itime=1,2
	CALL SCRMOD('REVERS')
	if (itime==1) CALL METAFL('CONS')
	if (itime==2) call METAFL(graphformat)
	call window(200,100,848,600)
	!CALL PAGE(1600,1000)
	CALL DISINI
	call ERRMOD("ALL","OFF")
	if (itime==1) then
		call WINTIT("Click right mouse button to close...")
		call hwfont
    else
		call COMPLX
    end if
	call shdpat(16)
	call SETRGB(Rcomp,Gcomp,Bcomp)
	call rectan(700,500,600,400)
	call SETRGB(Rcomp_op,Gcomp_op,Bcomp_op)
	call rectan(1600,500,600,400)
	call SETRGB(Rcomp_maxbr,Gcomp_maxbr,Bcomp_maxbr)
	call rectan(700,1200,600,400)
	call SETRGB(Rcomp_op_maxbr,Gcomp_op_maxbr,Bcomp_op_maxbr)
	call rectan(1600,1200,600,400)
	call SETRGB(0D0,0D0,0D0)
	call height(60)
	call messag("color",900,400)
    call messag("complementary color",1500,400)
    if (itime==1) then
		call messag("Maximum brightness of above colors:",750,1100)
    else
		call messag("Maximum brightness of above colors:",550,1100)
    end if
	CALL DISFIN
    if (itime==1) then
		write(*,*)
		write(*,*) "Do you want to save the colors to image file? (y/n)"
		read(*,*) selectyn
		if (selectyn=='n'.or.selectyn=='N') exit
    else
		write(*,*) "Done! The image file has been saved to current folder with ""DISLIN"" as prefix"
    end if
end do
end subroutine




!!---------- Generate CIE 1931 two-degree tristimulus functions, element index corresponds to nm value. Obtained from http://cvrl.ucl.ac.uk/cmfs.htm
subroutine gen_CIE_tristimulus_func(CIE_Xfunc,CIE_Yfunc,CIE_Zfunc)
real*8 CIE_Xfunc(830),CIE_Yfunc(830),CIE_Zfunc(830)
CIE_Xfunc(360)=0.000129900000D0;CIE_Yfunc(360)=0.000003917000D0;CIE_Zfunc(360)=0.000606100000D0 
CIE_Xfunc(361)=0.000145847000D0;CIE_Yfunc(361)=0.000004393581D0;CIE_Zfunc(361)=0.000680879200D0 
CIE_Xfunc(362)=0.000163802100D0;CIE_Yfunc(362)=0.000004929604D0;CIE_Zfunc(362)=0.000765145600D0 
CIE_Xfunc(363)=0.000184003700D0;CIE_Yfunc(363)=0.000005532136D0;CIE_Zfunc(363)=0.000860012400D0 
CIE_Xfunc(364)=0.000206690200D0;CIE_Yfunc(364)=0.000006208245D0;CIE_Zfunc(364)=0.000966592800D0 
CIE_Xfunc(365)=0.000232100000D0;CIE_Yfunc(365)=0.000006965000D0;CIE_Zfunc(365)=0.001086000000D0 
CIE_Xfunc(366)=0.000260728000D0;CIE_Yfunc(366)=0.000007813219D0;CIE_Zfunc(366)=0.001220586000D0 
CIE_Xfunc(367)=0.000293075000D0;CIE_Yfunc(367)=0.000008767336D0;CIE_Zfunc(367)=0.001372729000D0 
CIE_Xfunc(368)=0.000329388000D0;CIE_Yfunc(368)=0.000009839844D0;CIE_Zfunc(368)=0.001543579000D0 
CIE_Xfunc(369)=0.000369914000D0;CIE_Yfunc(369)=0.000011043230D0;CIE_Zfunc(369)=0.001734286000D0 
CIE_Xfunc(370)=0.000414900000D0;CIE_Yfunc(370)=0.000012390000D0;CIE_Zfunc(370)=0.001946000000D0 
CIE_Xfunc(371)=0.000464158700D0;CIE_Yfunc(371)=0.000013886410D0;CIE_Zfunc(371)=0.002177777000D0 
CIE_Xfunc(372)=0.000518986000D0;CIE_Yfunc(372)=0.000015557280D0;CIE_Zfunc(372)=0.002435809000D0 
CIE_Xfunc(373)=0.000581854000D0;CIE_Yfunc(373)=0.000017442960D0;CIE_Zfunc(373)=0.002731953000D0 
CIE_Xfunc(374)=0.000655234700D0;CIE_Yfunc(374)=0.000019583750D0;CIE_Zfunc(374)=0.003078064000D0 
CIE_Xfunc(375)=0.000741600000D0;CIE_Yfunc(375)=0.000022020000D0;CIE_Zfunc(375)=0.003486000000D0 
CIE_Xfunc(376)=0.000845029600D0;CIE_Yfunc(376)=0.000024839650D0;CIE_Zfunc(376)=0.003975227000D0 
CIE_Xfunc(377)=0.000964526800D0;CIE_Yfunc(377)=0.000028041260D0;CIE_Zfunc(377)=0.004540880000D0 
CIE_Xfunc(378)=0.001094949000D0;CIE_Yfunc(378)=0.000031531040D0;CIE_Zfunc(378)=0.005158320000D0 
CIE_Xfunc(379)=0.001231154000D0;CIE_Yfunc(379)=0.000035215210D0;CIE_Zfunc(379)=0.005802907000D0 
CIE_Xfunc(380)=0.001368000000D0;CIE_Yfunc(380)=0.000039000000D0;CIE_Zfunc(380)=0.006450001000D0 
CIE_Xfunc(381)=0.001502050000D0;CIE_Yfunc(381)=0.000042826400D0;CIE_Zfunc(381)=0.007083216000D0 
CIE_Xfunc(382)=0.001642328000D0;CIE_Yfunc(382)=0.000046914600D0;CIE_Zfunc(382)=0.007745488000D0 
CIE_Xfunc(383)=0.001802382000D0;CIE_Yfunc(383)=0.000051589600D0;CIE_Zfunc(383)=0.008501152000D0 
CIE_Xfunc(384)=0.001995757000D0;CIE_Yfunc(384)=0.000057176400D0;CIE_Zfunc(384)=0.009414544000D0 
CIE_Xfunc(385)=0.002236000000D0;CIE_Yfunc(385)=0.000064000000D0;CIE_Zfunc(385)=0.010549990000D0 
CIE_Xfunc(386)=0.002535385000D0;CIE_Yfunc(386)=0.000072344210D0;CIE_Zfunc(386)=0.011965800000D0 
CIE_Xfunc(387)=0.002892603000D0;CIE_Yfunc(387)=0.000082212240D0;CIE_Zfunc(387)=0.013655870000D0 
CIE_Xfunc(388)=0.003300829000D0;CIE_Yfunc(388)=0.000093508160D0;CIE_Zfunc(388)=0.015588050000D0 
CIE_Xfunc(389)=0.003753236000D0;CIE_Yfunc(389)=0.000106136100D0;CIE_Zfunc(389)=0.017730150000D0 
CIE_Xfunc(390)=0.004243000000D0;CIE_Yfunc(390)=0.000120000000D0;CIE_Zfunc(390)=0.020050010000D0 
CIE_Xfunc(391)=0.004762389000D0;CIE_Yfunc(391)=0.000134984000D0;CIE_Zfunc(391)=0.022511360000D0 
CIE_Xfunc(392)=0.005330048000D0;CIE_Yfunc(392)=0.000151492000D0;CIE_Zfunc(392)=0.025202880000D0 
CIE_Xfunc(393)=0.005978712000D0;CIE_Yfunc(393)=0.000170208000D0;CIE_Zfunc(393)=0.028279720000D0 
CIE_Xfunc(394)=0.006741117000D0;CIE_Yfunc(394)=0.000191816000D0;CIE_Zfunc(394)=0.031897040000D0 
CIE_Xfunc(395)=0.007650000000D0;CIE_Yfunc(395)=0.000217000000D0;CIE_Zfunc(395)=0.036210000000D0 
CIE_Xfunc(396)=0.008751373000D0;CIE_Yfunc(396)=0.000246906700D0;CIE_Zfunc(396)=0.041437710000D0 
CIE_Xfunc(397)=0.010028880000D0;CIE_Yfunc(397)=0.000281240000D0;CIE_Zfunc(397)=0.047503720000D0 
CIE_Xfunc(398)=0.011421700000D0;CIE_Yfunc(398)=0.000318520000D0;CIE_Zfunc(398)=0.054119880000D0 
CIE_Xfunc(399)=0.012869010000D0;CIE_Yfunc(399)=0.000357266700D0;CIE_Zfunc(399)=0.060998030000D0 
CIE_Xfunc(400)=0.014310000000D0;CIE_Yfunc(400)=0.000396000000D0;CIE_Zfunc(400)=0.067850010000D0 
CIE_Xfunc(401)=0.015704430000D0;CIE_Yfunc(401)=0.000433714700D0;CIE_Zfunc(401)=0.074486320000D0 
CIE_Xfunc(402)=0.017147440000D0;CIE_Yfunc(402)=0.000473024000D0;CIE_Zfunc(402)=0.081361560000D0 
CIE_Xfunc(403)=0.018781220000D0;CIE_Yfunc(403)=0.000517876000D0;CIE_Zfunc(403)=0.089153640000D0 
CIE_Xfunc(404)=0.020748010000D0;CIE_Yfunc(404)=0.000572218700D0;CIE_Zfunc(404)=0.098540480000D0 
CIE_Xfunc(405)=0.023190000000D0;CIE_Yfunc(405)=0.000640000000D0;CIE_Zfunc(405)=0.110200000000D0 
CIE_Xfunc(406)=0.026207360000D0;CIE_Yfunc(406)=0.000724560000D0;CIE_Zfunc(406)=0.124613300000D0 
CIE_Xfunc(407)=0.029782480000D0;CIE_Yfunc(407)=0.000825500000D0;CIE_Zfunc(407)=0.141701700000D0 
CIE_Xfunc(408)=0.033880920000D0;CIE_Yfunc(408)=0.000941160000D0;CIE_Zfunc(408)=0.161303500000D0 
CIE_Xfunc(409)=0.038468240000D0;CIE_Yfunc(409)=0.001069880000D0;CIE_Zfunc(409)=0.183256800000D0 
CIE_Xfunc(410)=0.043510000000D0;CIE_Yfunc(410)=0.001210000000D0;CIE_Zfunc(410)=0.207400000000D0 
CIE_Xfunc(411)=0.048995600000D0;CIE_Yfunc(411)=0.001362091000D0;CIE_Zfunc(411)=0.233692100000D0 
CIE_Xfunc(412)=0.055022600000D0;CIE_Yfunc(412)=0.001530752000D0;CIE_Zfunc(412)=0.262611400000D0 
CIE_Xfunc(413)=0.061718800000D0;CIE_Yfunc(413)=0.001720368000D0;CIE_Zfunc(413)=0.294774600000D0 
CIE_Xfunc(414)=0.069212000000D0;CIE_Yfunc(414)=0.001935323000D0;CIE_Zfunc(414)=0.330798500000D0 
CIE_Xfunc(415)=0.077630000000D0;CIE_Yfunc(415)=0.002180000000D0;CIE_Zfunc(415)=0.371300000000D0 
CIE_Xfunc(416)=0.086958110000D0;CIE_Yfunc(416)=0.002454800000D0;CIE_Zfunc(416)=0.416209100000D0 
CIE_Xfunc(417)=0.097176720000D0;CIE_Yfunc(417)=0.002764000000D0;CIE_Zfunc(417)=0.465464200000D0 
CIE_Xfunc(418)=0.108406300000D0;CIE_Yfunc(418)=0.003117800000D0;CIE_Zfunc(418)=0.519694800000D0 
CIE_Xfunc(419)=0.120767200000D0;CIE_Yfunc(419)=0.003526400000D0;CIE_Zfunc(419)=0.579530300000D0 
CIE_Xfunc(420)=0.134380000000D0;CIE_Yfunc(420)=0.004000000000D0;CIE_Zfunc(420)=0.645600000000D0 
CIE_Xfunc(421)=0.149358200000D0;CIE_Yfunc(421)=0.004546240000D0;CIE_Zfunc(421)=0.718483800000D0 
CIE_Xfunc(422)=0.165395700000D0;CIE_Yfunc(422)=0.005159320000D0;CIE_Zfunc(422)=0.796713300000D0 
CIE_Xfunc(423)=0.181983100000D0;CIE_Yfunc(423)=0.005829280000D0;CIE_Zfunc(423)=0.877845900000D0 
CIE_Xfunc(424)=0.198611000000D0;CIE_Yfunc(424)=0.006546160000D0;CIE_Zfunc(424)=0.959439000000D0 
CIE_Xfunc(425)=0.214770000000D0;CIE_Yfunc(425)=0.007300000000D0;CIE_Zfunc(425)=1.039050100000D0 
CIE_Xfunc(426)=0.230186800000D0;CIE_Yfunc(426)=0.008086507000D0;CIE_Zfunc(426)=1.115367300000D0 
CIE_Xfunc(427)=0.244879700000D0;CIE_Yfunc(427)=0.008908720000D0;CIE_Zfunc(427)=1.188497100000D0 
CIE_Xfunc(428)=0.258777300000D0;CIE_Yfunc(428)=0.009767680000D0;CIE_Zfunc(428)=1.258123300000D0 
CIE_Xfunc(429)=0.271807900000D0;CIE_Yfunc(429)=0.010664430000D0;CIE_Zfunc(429)=1.323929600000D0 
CIE_Xfunc(430)=0.283900000000D0;CIE_Yfunc(430)=0.011600000000D0;CIE_Zfunc(430)=1.385600000000D0 
CIE_Xfunc(431)=0.294943800000D0;CIE_Yfunc(431)=0.012573170000D0;CIE_Zfunc(431)=1.442635200000D0 
CIE_Xfunc(432)=0.304896500000D0;CIE_Yfunc(432)=0.013582720000D0;CIE_Zfunc(432)=1.494803500000D0 
CIE_Xfunc(433)=0.313787300000D0;CIE_Yfunc(433)=0.014629680000D0;CIE_Zfunc(433)=1.542190300000D0 
CIE_Xfunc(434)=0.321645400000D0;CIE_Yfunc(434)=0.015715090000D0;CIE_Zfunc(434)=1.584880700000D0 
CIE_Xfunc(435)=0.328500000000D0;CIE_Yfunc(435)=0.016840000000D0;CIE_Zfunc(435)=1.622960000000D0 
CIE_Xfunc(436)=0.334351300000D0;CIE_Yfunc(436)=0.018007360000D0;CIE_Zfunc(436)=1.656404800000D0 
CIE_Xfunc(437)=0.339210100000D0;CIE_Yfunc(437)=0.019214480000D0;CIE_Zfunc(437)=1.685295900000D0 
CIE_Xfunc(438)=0.343121300000D0;CIE_Yfunc(438)=0.020453920000D0;CIE_Zfunc(438)=1.709874500000D0 
CIE_Xfunc(439)=0.346129600000D0;CIE_Yfunc(439)=0.021718240000D0;CIE_Zfunc(439)=1.730382100000D0 
CIE_Xfunc(440)=0.348280000000D0;CIE_Yfunc(440)=0.023000000000D0;CIE_Zfunc(440)=1.747060000000D0 
CIE_Xfunc(441)=0.349599900000D0;CIE_Yfunc(441)=0.024294610000D0;CIE_Zfunc(441)=1.760044600000D0 
CIE_Xfunc(442)=0.350147400000D0;CIE_Yfunc(442)=0.025610240000D0;CIE_Zfunc(442)=1.769623300000D0 
CIE_Xfunc(443)=0.350013000000D0;CIE_Yfunc(443)=0.026958570000D0;CIE_Zfunc(443)=1.776263700000D0 
CIE_Xfunc(444)=0.349287000000D0;CIE_Yfunc(444)=0.028351250000D0;CIE_Zfunc(444)=1.780433400000D0 
CIE_Xfunc(445)=0.348060000000D0;CIE_Yfunc(445)=0.029800000000D0;CIE_Zfunc(445)=1.782600000000D0 
CIE_Xfunc(446)=0.346373300000D0;CIE_Yfunc(446)=0.031310830000D0;CIE_Zfunc(446)=1.782968200000D0 
CIE_Xfunc(447)=0.344262400000D0;CIE_Yfunc(447)=0.032883680000D0;CIE_Zfunc(447)=1.781699800000D0 
CIE_Xfunc(448)=0.341808800000D0;CIE_Yfunc(448)=0.034521120000D0;CIE_Zfunc(448)=1.779198200000D0 
CIE_Xfunc(449)=0.339094100000D0;CIE_Yfunc(449)=0.036225710000D0;CIE_Zfunc(449)=1.775867100000D0 
CIE_Xfunc(450)=0.336200000000D0;CIE_Yfunc(450)=0.038000000000D0;CIE_Zfunc(450)=1.772110000000D0 
CIE_Xfunc(451)=0.333197700000D0;CIE_Yfunc(451)=0.039846670000D0;CIE_Zfunc(451)=1.768258900000D0 
CIE_Xfunc(452)=0.330041100000D0;CIE_Yfunc(452)=0.041768000000D0;CIE_Zfunc(452)=1.764039000000D0 
CIE_Xfunc(453)=0.326635700000D0;CIE_Yfunc(453)=0.043766000000D0;CIE_Zfunc(453)=1.758943800000D0 
CIE_Xfunc(454)=0.322886800000D0;CIE_Yfunc(454)=0.045842670000D0;CIE_Zfunc(454)=1.752466300000D0 
CIE_Xfunc(455)=0.318700000000D0;CIE_Yfunc(455)=0.048000000000D0;CIE_Zfunc(455)=1.744100000000D0 
CIE_Xfunc(456)=0.314025100000D0;CIE_Yfunc(456)=0.050243680000D0;CIE_Zfunc(456)=1.733559500000D0 
CIE_Xfunc(457)=0.308884000000D0;CIE_Yfunc(457)=0.052573040000D0;CIE_Zfunc(457)=1.720858100000D0 
CIE_Xfunc(458)=0.303290400000D0;CIE_Yfunc(458)=0.054980560000D0;CIE_Zfunc(458)=1.705936900000D0 
CIE_Xfunc(459)=0.297257900000D0;CIE_Yfunc(459)=0.057458720000D0;CIE_Zfunc(459)=1.688737200000D0 
CIE_Xfunc(460)=0.290800000000D0;CIE_Yfunc(460)=0.060000000000D0;CIE_Zfunc(460)=1.669200000000D0 
CIE_Xfunc(461)=0.283970100000D0;CIE_Yfunc(461)=0.062601970000D0;CIE_Zfunc(461)=1.647528700000D0 
CIE_Xfunc(462)=0.276721400000D0;CIE_Yfunc(462)=0.065277520000D0;CIE_Zfunc(462)=1.623412700000D0 
CIE_Xfunc(463)=0.268917800000D0;CIE_Yfunc(463)=0.068042080000D0;CIE_Zfunc(463)=1.596022300000D0 
CIE_Xfunc(464)=0.260422700000D0;CIE_Yfunc(464)=0.070911090000D0;CIE_Zfunc(464)=1.564528000000D0 
CIE_Xfunc(465)=0.251100000000D0;CIE_Yfunc(465)=0.073900000000D0;CIE_Zfunc(465)=1.528100000000D0 
CIE_Xfunc(466)=0.240847500000D0;CIE_Yfunc(466)=0.077016000000D0;CIE_Zfunc(466)=1.486111400000D0 
CIE_Xfunc(467)=0.229851200000D0;CIE_Yfunc(467)=0.080266400000D0;CIE_Zfunc(467)=1.439521500000D0 
CIE_Xfunc(468)=0.218407200000D0;CIE_Yfunc(468)=0.083666800000D0;CIE_Zfunc(468)=1.389879900000D0 
CIE_Xfunc(469)=0.206811500000D0;CIE_Yfunc(469)=0.087232800000D0;CIE_Zfunc(469)=1.338736200000D0 
CIE_Xfunc(470)=0.195360000000D0;CIE_Yfunc(470)=0.090980000000D0;CIE_Zfunc(470)=1.287640000000D0 
CIE_Xfunc(471)=0.184213600000D0;CIE_Yfunc(471)=0.094917550000D0;CIE_Zfunc(471)=1.237422300000D0 
CIE_Xfunc(472)=0.173327300000D0;CIE_Yfunc(472)=0.099045840000D0;CIE_Zfunc(472)=1.187824300000D0 
CIE_Xfunc(473)=0.162688100000D0;CIE_Yfunc(473)=0.103367400000D0;CIE_Zfunc(473)=1.138761100000D0 
CIE_Xfunc(474)=0.152283300000D0;CIE_Yfunc(474)=0.107884600000D0;CIE_Zfunc(474)=1.090148000000D0 
CIE_Xfunc(475)=0.142100000000D0;CIE_Yfunc(475)=0.112600000000D0;CIE_Zfunc(475)=1.041900000000D0 
CIE_Xfunc(476)=0.132178600000D0;CIE_Yfunc(476)=0.117532000000D0;CIE_Zfunc(476)=0.994197600000D0 
CIE_Xfunc(477)=0.122569600000D0;CIE_Yfunc(477)=0.122674400000D0;CIE_Zfunc(477)=0.947347300000D0 
CIE_Xfunc(478)=0.113275200000D0;CIE_Yfunc(478)=0.127992800000D0;CIE_Zfunc(478)=0.901453100000D0 
CIE_Xfunc(479)=0.104297900000D0;CIE_Yfunc(479)=0.133452800000D0;CIE_Zfunc(479)=0.856619300000D0 
CIE_Xfunc(480)=0.095640000000D0;CIE_Yfunc(480)=0.139020000000D0;CIE_Zfunc(480)=0.812950100000D0 
CIE_Xfunc(481)=0.087299550000D0;CIE_Yfunc(481)=0.144676400000D0;CIE_Zfunc(481)=0.770517300000D0 
CIE_Xfunc(482)=0.079308040000D0;CIE_Yfunc(482)=0.150469300000D0;CIE_Zfunc(482)=0.729444800000D0 
CIE_Xfunc(483)=0.071717760000D0;CIE_Yfunc(483)=0.156461900000D0;CIE_Zfunc(483)=0.689913600000D0 
CIE_Xfunc(484)=0.064580990000D0;CIE_Yfunc(484)=0.162717700000D0;CIE_Zfunc(484)=0.652104900000D0 
CIE_Xfunc(485)=0.057950010000D0;CIE_Yfunc(485)=0.169300000000D0;CIE_Zfunc(485)=0.616200000000D0 
CIE_Xfunc(486)=0.051862110000D0;CIE_Yfunc(486)=0.176243100000D0;CIE_Zfunc(486)=0.582328600000D0 
CIE_Xfunc(487)=0.046281520000D0;CIE_Yfunc(487)=0.183558100000D0;CIE_Zfunc(487)=0.550416200000D0 
CIE_Xfunc(488)=0.041150880000D0;CIE_Yfunc(488)=0.191273500000D0;CIE_Zfunc(488)=0.520337600000D0 
CIE_Xfunc(489)=0.036412830000D0;CIE_Yfunc(489)=0.199418000000D0;CIE_Zfunc(489)=0.491967300000D0 
CIE_Xfunc(490)=0.032010000000D0;CIE_Yfunc(490)=0.208020000000D0;CIE_Zfunc(490)=0.465180000000D0 
CIE_Xfunc(491)=0.027917200000D0;CIE_Yfunc(491)=0.217119900000D0;CIE_Zfunc(491)=0.439924600000D0 
CIE_Xfunc(492)=0.024144400000D0;CIE_Yfunc(492)=0.226734500000D0;CIE_Zfunc(492)=0.416183600000D0 
CIE_Xfunc(493)=0.020687000000D0;CIE_Yfunc(493)=0.236857100000D0;CIE_Zfunc(493)=0.393882200000D0 
CIE_Xfunc(494)=0.017540400000D0;CIE_Yfunc(494)=0.247481200000D0;CIE_Zfunc(494)=0.372945900000D0 
CIE_Xfunc(495)=0.014700000000D0;CIE_Yfunc(495)=0.258600000000D0;CIE_Zfunc(495)=0.353300000000D0 
CIE_Xfunc(496)=0.012161790000D0;CIE_Yfunc(496)=0.270184900000D0;CIE_Zfunc(496)=0.334857800000D0 
CIE_Xfunc(497)=0.009919960000D0;CIE_Yfunc(497)=0.282293900000D0;CIE_Zfunc(497)=0.317552100000D0 
CIE_Xfunc(498)=0.007967240000D0;CIE_Yfunc(498)=0.295050500000D0;CIE_Zfunc(498)=0.301337500000D0 
CIE_Xfunc(499)=0.006296346000D0;CIE_Yfunc(499)=0.308578000000D0;CIE_Zfunc(499)=0.286168600000D0 
CIE_Xfunc(500)=0.004900000000D0;CIE_Yfunc(500)=0.323000000000D0;CIE_Zfunc(500)=0.272000000000D0 
CIE_Xfunc(501)=0.003777173000D0;CIE_Yfunc(501)=0.338402100000D0;CIE_Zfunc(501)=0.258817100000D0 
CIE_Xfunc(502)=0.002945320000D0;CIE_Yfunc(502)=0.354685800000D0;CIE_Zfunc(502)=0.246483800000D0 
CIE_Xfunc(503)=0.002424880000D0;CIE_Yfunc(503)=0.371698600000D0;CIE_Zfunc(503)=0.234771800000D0 
CIE_Xfunc(504)=0.002236293000D0;CIE_Yfunc(504)=0.389287500000D0;CIE_Zfunc(504)=0.223453300000D0 
CIE_Xfunc(505)=0.002400000000D0;CIE_Yfunc(505)=0.407300000000D0;CIE_Zfunc(505)=0.212300000000D0 
CIE_Xfunc(506)=0.002925520000D0;CIE_Yfunc(506)=0.425629900000D0;CIE_Zfunc(506)=0.201169200000D0 
CIE_Xfunc(507)=0.003836560000D0;CIE_Yfunc(507)=0.444309600000D0;CIE_Zfunc(507)=0.190119600000D0 
CIE_Xfunc(508)=0.005174840000D0;CIE_Yfunc(508)=0.463394400000D0;CIE_Zfunc(508)=0.179225400000D0 
CIE_Xfunc(509)=0.006982080000D0;CIE_Yfunc(509)=0.482939500000D0;CIE_Zfunc(509)=0.168560800000D0 
CIE_Xfunc(510)=0.009300000000D0;CIE_Yfunc(510)=0.503000000000D0;CIE_Zfunc(510)=0.158200000000D0 
CIE_Xfunc(511)=0.012149490000D0;CIE_Yfunc(511)=0.523569300000D0;CIE_Zfunc(511)=0.148138300000D0 
CIE_Xfunc(512)=0.015535880000D0;CIE_Yfunc(512)=0.544512000000D0;CIE_Zfunc(512)=0.138375800000D0 
CIE_Xfunc(513)=0.019477520000D0;CIE_Yfunc(513)=0.565690000000D0;CIE_Zfunc(513)=0.128994200000D0 
CIE_Xfunc(514)=0.023992770000D0;CIE_Yfunc(514)=0.586965300000D0;CIE_Zfunc(514)=0.120075100000D0 
CIE_Xfunc(515)=0.029100000000D0;CIE_Yfunc(515)=0.608200000000D0;CIE_Zfunc(515)=0.111700000000D0 
CIE_Xfunc(516)=0.034814850000D0;CIE_Yfunc(516)=0.629345600000D0;CIE_Zfunc(516)=0.103904800000D0 
CIE_Xfunc(517)=0.041120160000D0;CIE_Yfunc(517)=0.650306800000D0;CIE_Zfunc(517)=0.096667480000D0 
CIE_Xfunc(518)=0.047985040000D0;CIE_Yfunc(518)=0.670875200000D0;CIE_Zfunc(518)=0.089982720000D0 
CIE_Xfunc(519)=0.055378610000D0;CIE_Yfunc(519)=0.690842400000D0;CIE_Zfunc(519)=0.083845310000D0 
CIE_Xfunc(520)=0.063270000000D0;CIE_Yfunc(520)=0.710000000000D0;CIE_Zfunc(520)=0.078249990000D0 
CIE_Xfunc(521)=0.071635010000D0;CIE_Yfunc(521)=0.728185200000D0;CIE_Zfunc(521)=0.073208990000D0 
CIE_Xfunc(522)=0.080462240000D0;CIE_Yfunc(522)=0.745463600000D0;CIE_Zfunc(522)=0.068678160000D0 
CIE_Xfunc(523)=0.089739960000D0;CIE_Yfunc(523)=0.761969400000D0;CIE_Zfunc(523)=0.064567840000D0 
CIE_Xfunc(524)=0.099456450000D0;CIE_Yfunc(524)=0.777836800000D0;CIE_Zfunc(524)=0.060788350000D0 
CIE_Xfunc(525)=0.109600000000D0;CIE_Yfunc(525)=0.793200000000D0;CIE_Zfunc(525)=0.057250010000D0 
CIE_Xfunc(526)=0.120167400000D0;CIE_Yfunc(526)=0.808110400000D0;CIE_Zfunc(526)=0.053904350000D0 
CIE_Xfunc(527)=0.131114500000D0;CIE_Yfunc(527)=0.822496200000D0;CIE_Zfunc(527)=0.050746640000D0 
CIE_Xfunc(528)=0.142367900000D0;CIE_Yfunc(528)=0.836306800000D0;CIE_Zfunc(528)=0.047752760000D0 
CIE_Xfunc(529)=0.153854200000D0;CIE_Yfunc(529)=0.849491600000D0;CIE_Zfunc(529)=0.044898590000D0 
CIE_Xfunc(530)=0.165500000000D0;CIE_Yfunc(530)=0.862000000000D0;CIE_Zfunc(530)=0.042160000000D0 
CIE_Xfunc(531)=0.177257100000D0;CIE_Yfunc(531)=0.873810800000D0;CIE_Zfunc(531)=0.039507280000D0 
CIE_Xfunc(532)=0.189140000000D0;CIE_Yfunc(532)=0.884962400000D0;CIE_Zfunc(532)=0.036935640000D0 
CIE_Xfunc(533)=0.201169400000D0;CIE_Yfunc(533)=0.895493600000D0;CIE_Zfunc(533)=0.034458360000D0 
CIE_Xfunc(534)=0.213365800000D0;CIE_Yfunc(534)=0.905443200000D0;CIE_Zfunc(534)=0.032088720000D0 
CIE_Xfunc(535)=0.225749900000D0;CIE_Yfunc(535)=0.914850100000D0;CIE_Zfunc(535)=0.029840000000D0 
CIE_Xfunc(536)=0.238320900000D0;CIE_Yfunc(536)=0.923734800000D0;CIE_Zfunc(536)=0.027711810000D0 
CIE_Xfunc(537)=0.251066800000D0;CIE_Yfunc(537)=0.932092400000D0;CIE_Zfunc(537)=0.025694440000D0 
CIE_Xfunc(538)=0.263992200000D0;CIE_Yfunc(538)=0.939922600000D0;CIE_Zfunc(538)=0.023787160000D0 
CIE_Xfunc(539)=0.277101700000D0;CIE_Yfunc(539)=0.947225200000D0;CIE_Zfunc(539)=0.021989250000D0 
CIE_Xfunc(540)=0.290400000000D0;CIE_Yfunc(540)=0.954000000000D0;CIE_Zfunc(540)=0.020300000000D0 
CIE_Xfunc(541)=0.303891200000D0;CIE_Yfunc(541)=0.960256100000D0;CIE_Zfunc(541)=0.018718050000D0 
CIE_Xfunc(542)=0.317572600000D0;CIE_Yfunc(542)=0.966007400000D0;CIE_Zfunc(542)=0.017240360000D0 
CIE_Xfunc(543)=0.331438400000D0;CIE_Yfunc(543)=0.971260600000D0;CIE_Zfunc(543)=0.015863640000D0 
CIE_Xfunc(544)=0.345482800000D0;CIE_Yfunc(544)=0.976022500000D0;CIE_Zfunc(544)=0.014584610000D0 
CIE_Xfunc(545)=0.359700000000D0;CIE_Yfunc(545)=0.980300000000D0;CIE_Zfunc(545)=0.013400000000D0 
CIE_Xfunc(546)=0.374083900000D0;CIE_Yfunc(546)=0.984092400000D0;CIE_Zfunc(546)=0.012307230000D0 
CIE_Xfunc(547)=0.388639600000D0;CIE_Yfunc(547)=0.987418200000D0;CIE_Zfunc(547)=0.011301880000D0 
CIE_Xfunc(548)=0.403378400000D0;CIE_Yfunc(548)=0.990312800000D0;CIE_Zfunc(548)=0.010377920000D0 
CIE_Xfunc(549)=0.418311500000D0;CIE_Yfunc(549)=0.992811600000D0;CIE_Zfunc(549)=0.009529306000D0 
CIE_Xfunc(550)=0.433449900000D0;CIE_Yfunc(550)=0.994950100000D0;CIE_Zfunc(550)=0.008749999000D0 
CIE_Xfunc(551)=0.448795300000D0;CIE_Yfunc(551)=0.996710800000D0;CIE_Zfunc(551)=0.008035200000D0 
CIE_Xfunc(552)=0.464336000000D0;CIE_Yfunc(552)=0.998098300000D0;CIE_Zfunc(552)=0.007381600000D0 
CIE_Xfunc(553)=0.480064000000D0;CIE_Yfunc(553)=0.999112000000D0;CIE_Zfunc(553)=0.006785400000D0 
CIE_Xfunc(554)=0.495971300000D0;CIE_Yfunc(554)=0.999748200000D0;CIE_Zfunc(554)=0.006242800000D0 
CIE_Xfunc(555)=0.512050100000D0;CIE_Yfunc(555)=1.000000000000D0;CIE_Zfunc(555)=0.005749999000D0 
CIE_Xfunc(556)=0.528295900000D0;CIE_Yfunc(556)=0.999856700000D0;CIE_Zfunc(556)=0.005303600000D0 
CIE_Xfunc(557)=0.544691600000D0;CIE_Yfunc(557)=0.999304600000D0;CIE_Zfunc(557)=0.004899800000D0 
CIE_Xfunc(558)=0.561209400000D0;CIE_Yfunc(558)=0.998325500000D0;CIE_Zfunc(558)=0.004534200000D0 
CIE_Xfunc(559)=0.577821500000D0;CIE_Yfunc(559)=0.996898700000D0;CIE_Zfunc(559)=0.004202400000D0 
CIE_Xfunc(560)=0.594500000000D0;CIE_Yfunc(560)=0.995000000000D0;CIE_Zfunc(560)=0.003900000000D0 
CIE_Xfunc(561)=0.611220900000D0;CIE_Yfunc(561)=0.992600500000D0;CIE_Zfunc(561)=0.003623200000D0 
CIE_Xfunc(562)=0.627975800000D0;CIE_Yfunc(562)=0.989742600000D0;CIE_Zfunc(562)=0.003370600000D0 
CIE_Xfunc(563)=0.644760200000D0;CIE_Yfunc(563)=0.986444400000D0;CIE_Zfunc(563)=0.003141400000D0 
CIE_Xfunc(564)=0.661569700000D0;CIE_Yfunc(564)=0.982724100000D0;CIE_Zfunc(564)=0.002934800000D0 
CIE_Xfunc(565)=0.678400000000D0;CIE_Yfunc(565)=0.978600000000D0;CIE_Zfunc(565)=0.002749999000D0 
CIE_Xfunc(566)=0.695239200000D0;CIE_Yfunc(566)=0.974083700000D0;CIE_Zfunc(566)=0.002585200000D0 
CIE_Xfunc(567)=0.712058600000D0;CIE_Yfunc(567)=0.969171200000D0;CIE_Zfunc(567)=0.002438600000D0 
CIE_Xfunc(568)=0.728828400000D0;CIE_Yfunc(568)=0.963856800000D0;CIE_Zfunc(568)=0.002309400000D0 
CIE_Xfunc(569)=0.745518800000D0;CIE_Yfunc(569)=0.958134900000D0;CIE_Zfunc(569)=0.002196800000D0 
CIE_Xfunc(570)=0.762100000000D0;CIE_Yfunc(570)=0.952000000000D0;CIE_Zfunc(570)=0.002100000000D0 
CIE_Xfunc(571)=0.778543200000D0;CIE_Yfunc(571)=0.945450400000D0;CIE_Zfunc(571)=0.002017733000D0 
CIE_Xfunc(572)=0.794825600000D0;CIE_Yfunc(572)=0.938499200000D0;CIE_Zfunc(572)=0.001948200000D0 
CIE_Xfunc(573)=0.810926400000D0;CIE_Yfunc(573)=0.931162800000D0;CIE_Zfunc(573)=0.001889800000D0 
CIE_Xfunc(574)=0.826824800000D0;CIE_Yfunc(574)=0.923457600000D0;CIE_Zfunc(574)=0.001840933000D0 
CIE_Xfunc(575)=0.842500000000D0;CIE_Yfunc(575)=0.915400000000D0;CIE_Zfunc(575)=0.001800000000D0 
CIE_Xfunc(576)=0.857932500000D0;CIE_Yfunc(576)=0.907006400000D0;CIE_Zfunc(576)=0.001766267000D0 
CIE_Xfunc(577)=0.873081600000D0;CIE_Yfunc(577)=0.898277200000D0;CIE_Zfunc(577)=0.001737800000D0 
CIE_Xfunc(578)=0.887894400000D0;CIE_Yfunc(578)=0.889204800000D0;CIE_Zfunc(578)=0.001711200000D0 
CIE_Xfunc(579)=0.902318100000D0;CIE_Yfunc(579)=0.879781600000D0;CIE_Zfunc(579)=0.001683067000D0 
CIE_Xfunc(580)=0.916300000000D0;CIE_Yfunc(580)=0.870000000000D0;CIE_Zfunc(580)=0.001650001000D0 
CIE_Xfunc(581)=0.929799500000D0;CIE_Yfunc(581)=0.859861300000D0;CIE_Zfunc(581)=0.001610133000D0 
CIE_Xfunc(582)=0.942798400000D0;CIE_Yfunc(582)=0.849392000000D0;CIE_Zfunc(582)=0.001564400000D0 
CIE_Xfunc(583)=0.955277600000D0;CIE_Yfunc(583)=0.838622000000D0;CIE_Zfunc(583)=0.001513600000D0 
CIE_Xfunc(584)=0.967217900000D0;CIE_Yfunc(584)=0.827581300000D0;CIE_Zfunc(584)=0.001458533000D0 
CIE_Xfunc(585)=0.978600000000D0;CIE_Yfunc(585)=0.816300000000D0;CIE_Zfunc(585)=0.001400000000D0 
CIE_Xfunc(586)=0.989385600000D0;CIE_Yfunc(586)=0.804794700000D0;CIE_Zfunc(586)=0.001336667000D0 
CIE_Xfunc(587)=0.999548800000D0;CIE_Yfunc(587)=0.793082000000D0;CIE_Zfunc(587)=0.001270000000D0 
CIE_Xfunc(588)=1.009089200000D0;CIE_Yfunc(588)=0.781192000000D0;CIE_Zfunc(588)=0.001205000000D0 
CIE_Xfunc(589)=1.018006400000D0;CIE_Yfunc(589)=0.769154700000D0;CIE_Zfunc(589)=0.001146667000D0 
CIE_Xfunc(590)=1.026300000000D0;CIE_Yfunc(590)=0.757000000000D0;CIE_Zfunc(590)=0.001100000000D0 
CIE_Xfunc(591)=1.033982700000D0;CIE_Yfunc(591)=0.744754100000D0;CIE_Zfunc(591)=0.001068800000D0 
CIE_Xfunc(592)=1.040986000000D0;CIE_Yfunc(592)=0.732422400000D0;CIE_Zfunc(592)=0.001049400000D0 
CIE_Xfunc(593)=1.047188000000D0;CIE_Yfunc(593)=0.720003600000D0;CIE_Zfunc(593)=0.001035600000D0 
CIE_Xfunc(594)=1.052466700000D0;CIE_Yfunc(594)=0.707496500000D0;CIE_Zfunc(594)=0.001021200000D0 
CIE_Xfunc(595)=1.056700000000D0;CIE_Yfunc(595)=0.694900000000D0;CIE_Zfunc(595)=0.001000000000D0 
CIE_Xfunc(596)=1.059794400000D0;CIE_Yfunc(596)=0.682219200000D0;CIE_Zfunc(596)=0.000968640000D0 
CIE_Xfunc(597)=1.061799200000D0;CIE_Yfunc(597)=0.669471600000D0;CIE_Zfunc(597)=0.000929920000D0 
CIE_Xfunc(598)=1.062806800000D0;CIE_Yfunc(598)=0.656674400000D0;CIE_Zfunc(598)=0.000886880000D0 
CIE_Xfunc(599)=1.062909600000D0;CIE_Yfunc(599)=0.643844800000D0;CIE_Zfunc(599)=0.000842560000D0 
CIE_Xfunc(600)=1.062200000000D0;CIE_Yfunc(600)=0.631000000000D0;CIE_Zfunc(600)=0.000800000000D0 
CIE_Xfunc(601)=1.060735200000D0;CIE_Yfunc(601)=0.618155500000D0;CIE_Zfunc(601)=0.000760960000D0 
CIE_Xfunc(602)=1.058443600000D0;CIE_Yfunc(602)=0.605314400000D0;CIE_Zfunc(602)=0.000723680000D0 
CIE_Xfunc(603)=1.055224400000D0;CIE_Yfunc(603)=0.592475600000D0;CIE_Zfunc(603)=0.000685920000D0 
CIE_Xfunc(604)=1.050976800000D0;CIE_Yfunc(604)=0.579637900000D0;CIE_Zfunc(604)=0.000645440000D0 
CIE_Xfunc(605)=1.045600000000D0;CIE_Yfunc(605)=0.566800000000D0;CIE_Zfunc(605)=0.000600000000D0 
CIE_Xfunc(606)=1.039036900000D0;CIE_Yfunc(606)=0.553961100000D0;CIE_Zfunc(606)=0.000547866700D0 
CIE_Xfunc(607)=1.031360800000D0;CIE_Yfunc(607)=0.541137200000D0;CIE_Zfunc(607)=0.000491600000D0 
CIE_Xfunc(608)=1.022666200000D0;CIE_Yfunc(608)=0.528352800000D0;CIE_Zfunc(608)=0.000435400000D0 
CIE_Xfunc(609)=1.013047700000D0;CIE_Yfunc(609)=0.515632300000D0;CIE_Zfunc(609)=0.000383466700D0 
CIE_Xfunc(610)=1.002600000000D0;CIE_Yfunc(610)=0.503000000000D0;CIE_Zfunc(610)=0.000340000000D0 
CIE_Xfunc(611)=0.991367500000D0;CIE_Yfunc(611)=0.490468800000D0;CIE_Zfunc(611)=0.000307253300D0 
CIE_Xfunc(612)=0.979331400000D0;CIE_Yfunc(612)=0.478030400000D0;CIE_Zfunc(612)=0.000283160000D0 
CIE_Xfunc(613)=0.966491600000D0;CIE_Yfunc(613)=0.465677600000D0;CIE_Zfunc(613)=0.000265440000D0 
CIE_Xfunc(614)=0.952847900000D0;CIE_Yfunc(614)=0.453403200000D0;CIE_Zfunc(614)=0.000251813300D0 
CIE_Xfunc(615)=0.938400000000D0;CIE_Yfunc(615)=0.441200000000D0;CIE_Zfunc(615)=0.000240000000D0 
CIE_Xfunc(616)=0.923194000000D0;CIE_Yfunc(616)=0.429080000000D0;CIE_Zfunc(616)=0.000229546700D0 
CIE_Xfunc(617)=0.907244000000D0;CIE_Yfunc(617)=0.417036000000D0;CIE_Zfunc(617)=0.000220640000D0 
CIE_Xfunc(618)=0.890502000000D0;CIE_Yfunc(618)=0.405032000000D0;CIE_Zfunc(618)=0.000211960000D0 
CIE_Xfunc(619)=0.872920000000D0;CIE_Yfunc(619)=0.393032000000D0;CIE_Zfunc(619)=0.000202186700D0 
CIE_Xfunc(620)=0.854449900000D0;CIE_Yfunc(620)=0.381000000000D0;CIE_Zfunc(620)=0.000190000000D0 
CIE_Xfunc(621)=0.835084000000D0;CIE_Yfunc(621)=0.368918400000D0;CIE_Zfunc(621)=0.000174213300D0 
CIE_Xfunc(622)=0.814946000000D0;CIE_Yfunc(622)=0.356827200000D0;CIE_Zfunc(622)=0.000155640000D0 
CIE_Xfunc(623)=0.794186000000D0;CIE_Yfunc(623)=0.344776800000D0;CIE_Zfunc(623)=0.000135960000D0 
CIE_Xfunc(624)=0.772954000000D0;CIE_Yfunc(624)=0.332817600000D0;CIE_Zfunc(624)=0.000116853300D0 
CIE_Xfunc(625)=0.751400000000D0;CIE_Yfunc(625)=0.321000000000D0;CIE_Zfunc(625)=0.000100000000D0 
CIE_Xfunc(626)=0.729583600000D0;CIE_Yfunc(626)=0.309338100000D0;CIE_Zfunc(626)=0.000086133330D0 
CIE_Xfunc(627)=0.707588800000D0;CIE_Yfunc(627)=0.297850400000D0;CIE_Zfunc(627)=0.000074600000D0 
CIE_Xfunc(628)=0.685602200000D0;CIE_Yfunc(628)=0.286593600000D0;CIE_Zfunc(628)=0.000065000000D0 
CIE_Xfunc(629)=0.663810400000D0;CIE_Yfunc(629)=0.275624500000D0;CIE_Zfunc(629)=0.000056933330D0 
CIE_Xfunc(630)=0.642400000000D0;CIE_Yfunc(630)=0.265000000000D0;CIE_Zfunc(630)=0.000049999990D0 
CIE_Xfunc(631)=0.621514900000D0;CIE_Yfunc(631)=0.254763200000D0;CIE_Zfunc(631)=0.000044160000D0 
CIE_Xfunc(632)=0.601113800000D0;CIE_Yfunc(632)=0.244889600000D0;CIE_Zfunc(632)=0.000039480000D0 
CIE_Xfunc(633)=0.581105200000D0;CIE_Yfunc(633)=0.235334400000D0;CIE_Zfunc(633)=0.000035720000D0 
CIE_Xfunc(634)=0.561397700000D0;CIE_Yfunc(634)=0.226052800000D0;CIE_Zfunc(634)=0.000032640000D0 
CIE_Xfunc(635)=0.541900000000D0;CIE_Yfunc(635)=0.217000000000D0;CIE_Zfunc(635)=0.000030000000D0 
CIE_Xfunc(636)=0.522599500000D0;CIE_Yfunc(636)=0.208161600000D0;CIE_Zfunc(636)=0.000027653330D0 
CIE_Xfunc(637)=0.503546400000D0;CIE_Yfunc(637)=0.199548800000D0;CIE_Zfunc(637)=0.000025560000D0 
CIE_Xfunc(638)=0.484743600000D0;CIE_Yfunc(638)=0.191155200000D0;CIE_Zfunc(638)=0.000023640000D0 
CIE_Xfunc(639)=0.466193900000D0;CIE_Yfunc(639)=0.182974400000D0;CIE_Zfunc(639)=0.000021813330D0 
CIE_Xfunc(640)=0.447900000000D0;CIE_Yfunc(640)=0.175000000000D0;CIE_Zfunc(640)=0.000020000000D0 
CIE_Xfunc(641)=0.429861300000D0;CIE_Yfunc(641)=0.167223500000D0;CIE_Zfunc(641)=0.000018133330D0 
CIE_Xfunc(642)=0.412098000000D0;CIE_Yfunc(642)=0.159646400000D0;CIE_Zfunc(642)=0.000016200000D0 
CIE_Xfunc(643)=0.394644000000D0;CIE_Yfunc(643)=0.152277600000D0;CIE_Zfunc(643)=0.000014200000D0 
CIE_Xfunc(644)=0.377533300000D0;CIE_Yfunc(644)=0.145125900000D0;CIE_Zfunc(644)=0.000012133330D0 
CIE_Xfunc(645)=0.360800000000D0;CIE_Yfunc(645)=0.138200000000D0;CIE_Zfunc(645)=0.000010000000D0 
CIE_Xfunc(646)=0.344456300000D0;CIE_Yfunc(646)=0.131500300000D0;CIE_Zfunc(646)=0.000007733333D0 
CIE_Xfunc(647)=0.328516800000D0;CIE_Yfunc(647)=0.125024800000D0;CIE_Zfunc(647)=0.000005400000D0 
CIE_Xfunc(648)=0.313019200000D0;CIE_Yfunc(648)=0.118779200000D0;CIE_Zfunc(648)=0.000003200000D0 
CIE_Xfunc(649)=0.298001100000D0;CIE_Yfunc(649)=0.112769100000D0;CIE_Zfunc(649)=0.000001333333D0 
CIE_Xfunc(650)=0.283500000000D0;CIE_Yfunc(650)=0.107000000000D0;CIE_Zfunc(650)=0.000000000000D0 
CIE_Xfunc(651)=0.269544800000D0;CIE_Yfunc(651)=0.101476200000D0;CIE_Zfunc(651)=0.000000000000D0 
CIE_Xfunc(652)=0.256118400000D0;CIE_Yfunc(652)=0.096188640000D0;CIE_Zfunc(652)=0.000000000000D0 
CIE_Xfunc(653)=0.243189600000D0;CIE_Yfunc(653)=0.091122960000D0;CIE_Zfunc(653)=0.000000000000D0 
CIE_Xfunc(654)=0.230727200000D0;CIE_Yfunc(654)=0.086264850000D0;CIE_Zfunc(654)=0.000000000000D0 
CIE_Xfunc(655)=0.218700000000D0;CIE_Yfunc(655)=0.081600000000D0;CIE_Zfunc(655)=0.000000000000D0 
CIE_Xfunc(656)=0.207097100000D0;CIE_Yfunc(656)=0.077120640000D0;CIE_Zfunc(656)=0.000000000000D0 
CIE_Xfunc(657)=0.195923200000D0;CIE_Yfunc(657)=0.072825520000D0;CIE_Zfunc(657)=0.000000000000D0 
CIE_Xfunc(658)=0.185170800000D0;CIE_Yfunc(658)=0.068710080000D0;CIE_Zfunc(658)=0.000000000000D0 
CIE_Xfunc(659)=0.174832300000D0;CIE_Yfunc(659)=0.064769760000D0;CIE_Zfunc(659)=0.000000000000D0 
CIE_Xfunc(660)=0.164900000000D0;CIE_Yfunc(660)=0.061000000000D0;CIE_Zfunc(660)=0.000000000000D0 
CIE_Xfunc(661)=0.155366700000D0;CIE_Yfunc(661)=0.057396210000D0;CIE_Zfunc(661)=0.000000000000D0 
CIE_Xfunc(662)=0.146230000000D0;CIE_Yfunc(662)=0.053955040000D0;CIE_Zfunc(662)=0.000000000000D0 
CIE_Xfunc(663)=0.137490000000D0;CIE_Yfunc(663)=0.050673760000D0;CIE_Zfunc(663)=0.000000000000D0 
CIE_Xfunc(664)=0.129146700000D0;CIE_Yfunc(664)=0.047549650000D0;CIE_Zfunc(664)=0.000000000000D0 
CIE_Xfunc(665)=0.121200000000D0;CIE_Yfunc(665)=0.044580000000D0;CIE_Zfunc(665)=0.000000000000D0 
CIE_Xfunc(666)=0.113639700000D0;CIE_Yfunc(666)=0.041758720000D0;CIE_Zfunc(666)=0.000000000000D0 
CIE_Xfunc(667)=0.106465000000D0;CIE_Yfunc(667)=0.039084960000D0;CIE_Zfunc(667)=0.000000000000D0 
CIE_Xfunc(668)=0.099690440000D0;CIE_Yfunc(668)=0.036563840000D0;CIE_Zfunc(668)=0.000000000000D0 
CIE_Xfunc(669)=0.093330610000D0;CIE_Yfunc(669)=0.034200480000D0;CIE_Zfunc(669)=0.000000000000D0 
CIE_Xfunc(670)=0.087400000000D0;CIE_Yfunc(670)=0.032000000000D0;CIE_Zfunc(670)=0.000000000000D0 
CIE_Xfunc(671)=0.081900960000D0;CIE_Yfunc(671)=0.029962610000D0;CIE_Zfunc(671)=0.000000000000D0 
CIE_Xfunc(672)=0.076804280000D0;CIE_Yfunc(672)=0.028076640000D0;CIE_Zfunc(672)=0.000000000000D0 
CIE_Xfunc(673)=0.072077120000D0;CIE_Yfunc(673)=0.026329360000D0;CIE_Zfunc(673)=0.000000000000D0 
CIE_Xfunc(674)=0.067686640000D0;CIE_Yfunc(674)=0.024708050000D0;CIE_Zfunc(674)=0.000000000000D0 
CIE_Xfunc(675)=0.063600000000D0;CIE_Yfunc(675)=0.023200000000D0;CIE_Zfunc(675)=0.000000000000D0 
CIE_Xfunc(676)=0.059806850000D0;CIE_Yfunc(676)=0.021800770000D0;CIE_Zfunc(676)=0.000000000000D0 
CIE_Xfunc(677)=0.056282160000D0;CIE_Yfunc(677)=0.020501120000D0;CIE_Zfunc(677)=0.000000000000D0 
CIE_Xfunc(678)=0.052971040000D0;CIE_Yfunc(678)=0.019281080000D0;CIE_Zfunc(678)=0.000000000000D0 
CIE_Xfunc(679)=0.049818610000D0;CIE_Yfunc(679)=0.018120690000D0;CIE_Zfunc(679)=0.000000000000D0 
CIE_Xfunc(680)=0.046770000000D0;CIE_Yfunc(680)=0.017000000000D0;CIE_Zfunc(680)=0.000000000000D0 
CIE_Xfunc(681)=0.043784050000D0;CIE_Yfunc(681)=0.015903790000D0;CIE_Zfunc(681)=0.000000000000D0 
CIE_Xfunc(682)=0.040875360000D0;CIE_Yfunc(682)=0.014837180000D0;CIE_Zfunc(682)=0.000000000000D0 
CIE_Xfunc(683)=0.038072640000D0;CIE_Yfunc(683)=0.013810680000D0;CIE_Zfunc(683)=0.000000000000D0 
CIE_Xfunc(684)=0.035404610000D0;CIE_Yfunc(684)=0.012834780000D0;CIE_Zfunc(684)=0.000000000000D0 
CIE_Xfunc(685)=0.032900000000D0;CIE_Yfunc(685)=0.011920000000D0;CIE_Zfunc(685)=0.000000000000D0 
CIE_Xfunc(686)=0.030564190000D0;CIE_Yfunc(686)=0.011068310000D0;CIE_Zfunc(686)=0.000000000000D0 
CIE_Xfunc(687)=0.028380560000D0;CIE_Yfunc(687)=0.010273390000D0;CIE_Zfunc(687)=0.000000000000D0 
CIE_Xfunc(688)=0.026344840000D0;CIE_Yfunc(688)=0.009533311000D0;CIE_Zfunc(688)=0.000000000000D0 
CIE_Xfunc(689)=0.024452750000D0;CIE_Yfunc(689)=0.008846157000D0;CIE_Zfunc(689)=0.000000000000D0 
CIE_Xfunc(690)=0.022700000000D0;CIE_Yfunc(690)=0.008210000000D0;CIE_Zfunc(690)=0.000000000000D0 
CIE_Xfunc(691)=0.021084290000D0;CIE_Yfunc(691)=0.007623781000D0;CIE_Zfunc(691)=0.000000000000D0 
CIE_Xfunc(692)=0.019599880000D0;CIE_Yfunc(692)=0.007085424000D0;CIE_Zfunc(692)=0.000000000000D0 
CIE_Xfunc(693)=0.018237320000D0;CIE_Yfunc(693)=0.006591476000D0;CIE_Zfunc(693)=0.000000000000D0 
CIE_Xfunc(694)=0.016987170000D0;CIE_Yfunc(694)=0.006138485000D0;CIE_Zfunc(694)=0.000000000000D0 
CIE_Xfunc(695)=0.015840000000D0;CIE_Yfunc(695)=0.005723000000D0;CIE_Zfunc(695)=0.000000000000D0 
CIE_Xfunc(696)=0.014790640000D0;CIE_Yfunc(696)=0.005343059000D0;CIE_Zfunc(696)=0.000000000000D0 
CIE_Xfunc(697)=0.013831320000D0;CIE_Yfunc(697)=0.004995796000D0;CIE_Zfunc(697)=0.000000000000D0 
CIE_Xfunc(698)=0.012948680000D0;CIE_Yfunc(698)=0.004676404000D0;CIE_Zfunc(698)=0.000000000000D0 
CIE_Xfunc(699)=0.012129200000D0;CIE_Yfunc(699)=0.004380075000D0;CIE_Zfunc(699)=0.000000000000D0 
CIE_Xfunc(700)=0.011359160000D0;CIE_Yfunc(700)=0.004102000000D0;CIE_Zfunc(700)=0.000000000000D0 
CIE_Xfunc(701)=0.010629350000D0;CIE_Yfunc(701)=0.003838453000D0;CIE_Zfunc(701)=0.000000000000D0 
CIE_Xfunc(702)=0.009938846000D0;CIE_Yfunc(702)=0.003589099000D0;CIE_Zfunc(702)=0.000000000000D0 
CIE_Xfunc(703)=0.009288422000D0;CIE_Yfunc(703)=0.003354219000D0;CIE_Zfunc(703)=0.000000000000D0 
CIE_Xfunc(704)=0.008678854000D0;CIE_Yfunc(704)=0.003134093000D0;CIE_Zfunc(704)=0.000000000000D0 
CIE_Xfunc(705)=0.008110916000D0;CIE_Yfunc(705)=0.002929000000D0;CIE_Zfunc(705)=0.000000000000D0 
CIE_Xfunc(706)=0.007582388000D0;CIE_Yfunc(706)=0.002738139000D0;CIE_Zfunc(706)=0.000000000000D0 
CIE_Xfunc(707)=0.007088746000D0;CIE_Yfunc(707)=0.002559876000D0;CIE_Zfunc(707)=0.000000000000D0 
CIE_Xfunc(708)=0.006627313000D0;CIE_Yfunc(708)=0.002393244000D0;CIE_Zfunc(708)=0.000000000000D0 
CIE_Xfunc(709)=0.006195408000D0;CIE_Yfunc(709)=0.002237275000D0;CIE_Zfunc(709)=0.000000000000D0 
CIE_Xfunc(710)=0.005790346000D0;CIE_Yfunc(710)=0.002091000000D0;CIE_Zfunc(710)=0.000000000000D0 
CIE_Xfunc(711)=0.005409826000D0;CIE_Yfunc(711)=0.001953587000D0;CIE_Zfunc(711)=0.000000000000D0 
CIE_Xfunc(712)=0.005052583000D0;CIE_Yfunc(712)=0.001824580000D0;CIE_Zfunc(712)=0.000000000000D0 
CIE_Xfunc(713)=0.004717512000D0;CIE_Yfunc(713)=0.001703580000D0;CIE_Zfunc(713)=0.000000000000D0 
CIE_Xfunc(714)=0.004403507000D0;CIE_Yfunc(714)=0.001590187000D0;CIE_Zfunc(714)=0.000000000000D0 
CIE_Xfunc(715)=0.004109457000D0;CIE_Yfunc(715)=0.001484000000D0;CIE_Zfunc(715)=0.000000000000D0 
CIE_Xfunc(716)=0.003833913000D0;CIE_Yfunc(716)=0.001384496000D0;CIE_Zfunc(716)=0.000000000000D0 
CIE_Xfunc(717)=0.003575748000D0;CIE_Yfunc(717)=0.001291268000D0;CIE_Zfunc(717)=0.000000000000D0 
CIE_Xfunc(718)=0.003334342000D0;CIE_Yfunc(718)=0.001204092000D0;CIE_Zfunc(718)=0.000000000000D0 
CIE_Xfunc(719)=0.003109075000D0;CIE_Yfunc(719)=0.001122744000D0;CIE_Zfunc(719)=0.000000000000D0 
CIE_Xfunc(720)=0.002899327000D0;CIE_Yfunc(720)=0.001047000000D0;CIE_Zfunc(720)=0.000000000000D0 
CIE_Xfunc(721)=0.002704348000D0;CIE_Yfunc(721)=0.000976589600D0;CIE_Zfunc(721)=0.000000000000D0 
CIE_Xfunc(722)=0.002523020000D0;CIE_Yfunc(722)=0.000911108800D0;CIE_Zfunc(722)=0.000000000000D0 
CIE_Xfunc(723)=0.002354168000D0;CIE_Yfunc(723)=0.000850133200D0;CIE_Zfunc(723)=0.000000000000D0 
CIE_Xfunc(724)=0.002196616000D0;CIE_Yfunc(724)=0.000793238400D0;CIE_Zfunc(724)=0.000000000000D0 
CIE_Xfunc(725)=0.002049190000D0;CIE_Yfunc(725)=0.000740000000D0;CIE_Zfunc(725)=0.000000000000D0 
CIE_Xfunc(726)=0.001910960000D0;CIE_Yfunc(726)=0.000690082700D0;CIE_Zfunc(726)=0.000000000000D0 
CIE_Xfunc(727)=0.001781438000D0;CIE_Yfunc(727)=0.000643310000D0;CIE_Zfunc(727)=0.000000000000D0 
CIE_Xfunc(728)=0.001660110000D0;CIE_Yfunc(728)=0.000599496000D0;CIE_Zfunc(728)=0.000000000000D0 
CIE_Xfunc(729)=0.001546459000D0;CIE_Yfunc(729)=0.000558454700D0;CIE_Zfunc(729)=0.000000000000D0 
CIE_Xfunc(730)=0.001439971000D0;CIE_Yfunc(730)=0.000520000000D0;CIE_Zfunc(730)=0.000000000000D0 
CIE_Xfunc(731)=0.001340042000D0;CIE_Yfunc(731)=0.000483913600D0;CIE_Zfunc(731)=0.000000000000D0 
CIE_Xfunc(732)=0.001246275000D0;CIE_Yfunc(732)=0.000450052800D0;CIE_Zfunc(732)=0.000000000000D0 
CIE_Xfunc(733)=0.001158471000D0;CIE_Yfunc(733)=0.000418345200D0;CIE_Zfunc(733)=0.000000000000D0 
CIE_Xfunc(734)=0.001076430000D0;CIE_Yfunc(734)=0.000388718400D0;CIE_Zfunc(734)=0.000000000000D0 
CIE_Xfunc(735)=0.000999949300D0;CIE_Yfunc(735)=0.000361100000D0;CIE_Zfunc(735)=0.000000000000D0 
CIE_Xfunc(736)=0.000928735800D0;CIE_Yfunc(736)=0.000335383500D0;CIE_Zfunc(736)=0.000000000000D0 
CIE_Xfunc(737)=0.000862433200D0;CIE_Yfunc(737)=0.000311440400D0;CIE_Zfunc(737)=0.000000000000D0 
CIE_Xfunc(738)=0.000800750300D0;CIE_Yfunc(738)=0.000289165600D0;CIE_Zfunc(738)=0.000000000000D0 
CIE_Xfunc(739)=0.000743396000D0;CIE_Yfunc(739)=0.000268453900D0;CIE_Zfunc(739)=0.000000000000D0 
CIE_Xfunc(740)=0.000690078600D0;CIE_Yfunc(740)=0.000249200000D0;CIE_Zfunc(740)=0.000000000000D0 
CIE_Xfunc(741)=0.000640515600D0;CIE_Yfunc(741)=0.000231301900D0;CIE_Zfunc(741)=0.000000000000D0 
CIE_Xfunc(742)=0.000594502100D0;CIE_Yfunc(742)=0.000214685600D0;CIE_Zfunc(742)=0.000000000000D0 
CIE_Xfunc(743)=0.000551864600D0;CIE_Yfunc(743)=0.000199288400D0;CIE_Zfunc(743)=0.000000000000D0 
CIE_Xfunc(744)=0.000512429000D0;CIE_Yfunc(744)=0.000185047500D0;CIE_Zfunc(744)=0.000000000000D0 
CIE_Xfunc(745)=0.000476021300D0;CIE_Yfunc(745)=0.000171900000D0;CIE_Zfunc(745)=0.000000000000D0 
CIE_Xfunc(746)=0.000442453600D0;CIE_Yfunc(746)=0.000159778100D0;CIE_Zfunc(746)=0.000000000000D0 
CIE_Xfunc(747)=0.000411511700D0;CIE_Yfunc(747)=0.000148604400D0;CIE_Zfunc(747)=0.000000000000D0 
CIE_Xfunc(748)=0.000382981400D0;CIE_Yfunc(748)=0.000138301600D0;CIE_Zfunc(748)=0.000000000000D0 
CIE_Xfunc(749)=0.000356649100D0;CIE_Yfunc(749)=0.000128792500D0;CIE_Zfunc(749)=0.000000000000D0 
CIE_Xfunc(750)=0.000332301100D0;CIE_Yfunc(750)=0.000120000000D0;CIE_Zfunc(750)=0.000000000000D0 
CIE_Xfunc(751)=0.000309758600D0;CIE_Yfunc(751)=0.000111859500D0;CIE_Zfunc(751)=0.000000000000D0 
CIE_Xfunc(752)=0.000288887100D0;CIE_Yfunc(752)=0.000104322400D0;CIE_Zfunc(752)=0.000000000000D0 
CIE_Xfunc(753)=0.000269539400D0;CIE_Yfunc(753)=0.000097335600D0;CIE_Zfunc(753)=0.000000000000D0 
CIE_Xfunc(754)=0.000251568200D0;CIE_Yfunc(754)=0.000090845870D0;CIE_Zfunc(754)=0.000000000000D0 
CIE_Xfunc(755)=0.000234826100D0;CIE_Yfunc(755)=0.000084800000D0;CIE_Zfunc(755)=0.000000000000D0 
CIE_Xfunc(756)=0.000219171000D0;CIE_Yfunc(756)=0.000079146670D0;CIE_Zfunc(756)=0.000000000000D0 
CIE_Xfunc(757)=0.000204525800D0;CIE_Yfunc(757)=0.000073858000D0;CIE_Zfunc(757)=0.000000000000D0 
CIE_Xfunc(758)=0.000190840500D0;CIE_Yfunc(758)=0.000068916000D0;CIE_Zfunc(758)=0.000000000000D0 
CIE_Xfunc(759)=0.000178065400D0;CIE_Yfunc(759)=0.000064302670D0;CIE_Zfunc(759)=0.000000000000D0 
CIE_Xfunc(760)=0.000166150500D0;CIE_Yfunc(760)=0.000060000000D0;CIE_Zfunc(760)=0.000000000000D0 
CIE_Xfunc(761)=0.000155023600D0;CIE_Yfunc(761)=0.000055981870D0;CIE_Zfunc(761)=0.000000000000D0 
CIE_Xfunc(762)=0.000144621900D0;CIE_Yfunc(762)=0.000052225600D0;CIE_Zfunc(762)=0.000000000000D0 
CIE_Xfunc(763)=0.000134909800D0;CIE_Yfunc(763)=0.000048718400D0;CIE_Zfunc(763)=0.000000000000D0 
CIE_Xfunc(764)=0.000125852000D0;CIE_Yfunc(764)=0.000045447470D0;CIE_Zfunc(764)=0.000000000000D0 
CIE_Xfunc(765)=0.000117413000D0;CIE_Yfunc(765)=0.000042400000D0;CIE_Zfunc(765)=0.000000000000D0 
CIE_Xfunc(766)=0.000109551500D0;CIE_Yfunc(766)=0.000039561040D0;CIE_Zfunc(766)=0.000000000000D0 
CIE_Xfunc(767)=0.000102224500D0;CIE_Yfunc(767)=0.000036915120D0;CIE_Zfunc(767)=0.000000000000D0 
CIE_Xfunc(768)=0.000095394450D0;CIE_Yfunc(768)=0.000034448680D0;CIE_Zfunc(768)=0.000000000000D0 
CIE_Xfunc(769)=0.000089023900D0;CIE_Yfunc(769)=0.000032148160D0;CIE_Zfunc(769)=0.000000000000D0 
CIE_Xfunc(770)=0.000083075270D0;CIE_Yfunc(770)=0.000030000000D0;CIE_Zfunc(770)=0.000000000000D0 
CIE_Xfunc(771)=0.000077512690D0;CIE_Yfunc(771)=0.000027991250D0;CIE_Zfunc(771)=0.000000000000D0 
CIE_Xfunc(772)=0.000072313040D0;CIE_Yfunc(772)=0.000026113560D0;CIE_Zfunc(772)=0.000000000000D0 
CIE_Xfunc(773)=0.000067457780D0;CIE_Yfunc(773)=0.000024360240D0;CIE_Zfunc(773)=0.000000000000D0 
CIE_Xfunc(774)=0.000062928440D0;CIE_Yfunc(774)=0.000022724610D0;CIE_Zfunc(774)=0.000000000000D0 
CIE_Xfunc(775)=0.000058706520D0;CIE_Yfunc(775)=0.000021200000D0;CIE_Zfunc(775)=0.000000000000D0 
CIE_Xfunc(776)=0.000054770280D0;CIE_Yfunc(776)=0.000019778550D0;CIE_Zfunc(776)=0.000000000000D0 
CIE_Xfunc(777)=0.000051099180D0;CIE_Yfunc(777)=0.000018452850D0;CIE_Zfunc(777)=0.000000000000D0 
CIE_Xfunc(778)=0.000047676540D0;CIE_Yfunc(778)=0.000017216870D0;CIE_Zfunc(778)=0.000000000000D0 
CIE_Xfunc(779)=0.000044485670D0;CIE_Yfunc(779)=0.000016064590D0;CIE_Zfunc(779)=0.000000000000D0 
CIE_Xfunc(780)=0.000041509940D0;CIE_Yfunc(780)=0.000014990000D0;CIE_Zfunc(780)=0.000000000000D0 
CIE_Xfunc(781)=0.000038733240D0;CIE_Yfunc(781)=0.000013987280D0;CIE_Zfunc(781)=0.000000000000D0 
CIE_Xfunc(782)=0.000036142030D0;CIE_Yfunc(782)=0.000013051550D0;CIE_Zfunc(782)=0.000000000000D0 
CIE_Xfunc(783)=0.000033723520D0;CIE_Yfunc(783)=0.000012178180D0;CIE_Zfunc(783)=0.000000000000D0 
CIE_Xfunc(784)=0.000031464870D0;CIE_Yfunc(784)=0.000011362540D0;CIE_Zfunc(784)=0.000000000000D0 
CIE_Xfunc(785)=0.000029353260D0;CIE_Yfunc(785)=0.000010600000D0;CIE_Zfunc(785)=0.000000000000D0 
CIE_Xfunc(786)=0.000027375730D0;CIE_Yfunc(786)=0.000009885877D0;CIE_Zfunc(786)=0.000000000000D0 
CIE_Xfunc(787)=0.000025524330D0;CIE_Yfunc(787)=0.000009217304D0;CIE_Zfunc(787)=0.000000000000D0 
CIE_Xfunc(788)=0.000023793760D0;CIE_Yfunc(788)=0.000008592362D0;CIE_Zfunc(788)=0.000000000000D0 
CIE_Xfunc(789)=0.000022178700D0;CIE_Yfunc(789)=0.000008009133D0;CIE_Zfunc(789)=0.000000000000D0 
CIE_Xfunc(790)=0.000020673830D0;CIE_Yfunc(790)=0.000007465700D0;CIE_Zfunc(790)=0.000000000000D0 
CIE_Xfunc(791)=0.000019272260D0;CIE_Yfunc(791)=0.000006959567D0;CIE_Zfunc(791)=0.000000000000D0 
CIE_Xfunc(792)=0.000017966400D0;CIE_Yfunc(792)=0.000006487995D0;CIE_Zfunc(792)=0.000000000000D0 
CIE_Xfunc(793)=0.000016749910D0;CIE_Yfunc(793)=0.000006048699D0;CIE_Zfunc(793)=0.000000000000D0 
CIE_Xfunc(794)=0.000015616480D0;CIE_Yfunc(794)=0.000005639396D0;CIE_Zfunc(794)=0.000000000000D0 
CIE_Xfunc(795)=0.000014559770D0;CIE_Yfunc(795)=0.000005257800D0;CIE_Zfunc(795)=0.000000000000D0 
CIE_Xfunc(796)=0.000013573870D0;CIE_Yfunc(796)=0.000004901771D0;CIE_Zfunc(796)=0.000000000000D0 
CIE_Xfunc(797)=0.000012654360D0;CIE_Yfunc(797)=0.000004569720D0;CIE_Zfunc(797)=0.000000000000D0 
CIE_Xfunc(798)=0.000011797230D0;CIE_Yfunc(798)=0.000004260194D0;CIE_Zfunc(798)=0.000000000000D0 
CIE_Xfunc(799)=0.000010998440D0;CIE_Yfunc(799)=0.000003971739D0;CIE_Zfunc(799)=0.000000000000D0 
CIE_Xfunc(800)=0.000010253980D0;CIE_Yfunc(800)=0.000003702900D0;CIE_Zfunc(800)=0.000000000000D0 
CIE_Xfunc(801)=0.000009559646D0;CIE_Yfunc(801)=0.000003452163D0;CIE_Zfunc(801)=0.000000000000D0 
CIE_Xfunc(802)=0.000008912044D0;CIE_Yfunc(802)=0.000003218302D0;CIE_Zfunc(802)=0.000000000000D0 
CIE_Xfunc(803)=0.000008308358D0;CIE_Yfunc(803)=0.000003000300D0;CIE_Zfunc(803)=0.000000000000D0 
CIE_Xfunc(804)=0.000007745769D0;CIE_Yfunc(804)=0.000002797139D0;CIE_Zfunc(804)=0.000000000000D0 
CIE_Xfunc(805)=0.000007221456D0;CIE_Yfunc(805)=0.000002607800D0;CIE_Zfunc(805)=0.000000000000D0 
CIE_Xfunc(806)=0.000006732475D0;CIE_Yfunc(806)=0.000002431220D0;CIE_Zfunc(806)=0.000000000000D0 
CIE_Xfunc(807)=0.000006276423D0;CIE_Yfunc(807)=0.000002266531D0;CIE_Zfunc(807)=0.000000000000D0 
CIE_Xfunc(808)=0.000005851304D0;CIE_Yfunc(808)=0.000002113013D0;CIE_Zfunc(808)=0.000000000000D0 
CIE_Xfunc(809)=0.000005455118D0;CIE_Yfunc(809)=0.000001969943D0;CIE_Zfunc(809)=0.000000000000D0 
CIE_Xfunc(810)=0.000005085868D0;CIE_Yfunc(810)=0.000001836600D0;CIE_Zfunc(810)=0.000000000000D0 
CIE_Xfunc(811)=0.000004741466D0;CIE_Yfunc(811)=0.000001712230D0;CIE_Zfunc(811)=0.000000000000D0 
CIE_Xfunc(812)=0.000004420236D0;CIE_Yfunc(812)=0.000001596228D0;CIE_Zfunc(812)=0.000000000000D0 
CIE_Xfunc(813)=0.000004120783D0;CIE_Yfunc(813)=0.000001488090D0;CIE_Zfunc(813)=0.000000000000D0 
CIE_Xfunc(814)=0.000003841716D0;CIE_Yfunc(814)=0.000001387314D0;CIE_Zfunc(814)=0.000000000000D0 
CIE_Xfunc(815)=0.000003581652D0;CIE_Yfunc(815)=0.000001293400D0;CIE_Zfunc(815)=0.000000000000D0 
CIE_Xfunc(816)=0.000003339127D0;CIE_Yfunc(816)=0.000001205820D0;CIE_Zfunc(816)=0.000000000000D0 
CIE_Xfunc(817)=0.000003112949D0;CIE_Yfunc(817)=0.000001124143D0;CIE_Zfunc(817)=0.000000000000D0 
CIE_Xfunc(818)=0.000002902121D0;CIE_Yfunc(818)=0.000001048009D0;CIE_Zfunc(818)=0.000000000000D0 
CIE_Xfunc(819)=0.000002705645D0;CIE_Yfunc(819)=0.000000977058D0;CIE_Zfunc(819)=0.000000000000D0 
CIE_Xfunc(820)=0.000002522525D0;CIE_Yfunc(820)=0.000000910930D0;CIE_Zfunc(820)=0.000000000000D0 
CIE_Xfunc(821)=0.000002351726D0;CIE_Yfunc(821)=0.000000849251D0;CIE_Zfunc(821)=0.000000000000D0 
CIE_Xfunc(822)=0.000002192415D0;CIE_Yfunc(822)=0.000000791721D0;CIE_Zfunc(822)=0.000000000000D0 
CIE_Xfunc(823)=0.000002043902D0;CIE_Yfunc(823)=0.000000738090D0;CIE_Zfunc(823)=0.000000000000D0 
CIE_Xfunc(824)=0.000001905497D0;CIE_Yfunc(824)=0.000000688110D0;CIE_Zfunc(824)=0.000000000000D0 
CIE_Xfunc(825)=0.000001776509D0;CIE_Yfunc(825)=0.000000641530D0;CIE_Zfunc(825)=0.000000000000D0 
CIE_Xfunc(826)=0.000001656215D0;CIE_Yfunc(826)=0.000000598090D0;CIE_Zfunc(826)=0.000000000000D0 
CIE_Xfunc(827)=0.000001544022D0;CIE_Yfunc(827)=0.000000557575D0;CIE_Zfunc(827)=0.000000000000D0 
CIE_Xfunc(828)=0.000001439440D0;CIE_Yfunc(828)=0.000000519808D0;CIE_Zfunc(828)=0.000000000000D0 
CIE_Xfunc(829)=0.000001341977D0;CIE_Yfunc(829)=0.000000484612D0;CIE_Zfunc(829)=0.000000000000D0 
CIE_Xfunc(830)=0.000001251141D0;CIE_Yfunc(830)=0.000000451810D0;CIE_Zfunc(830)=0.000000000000D0 
end subroutine