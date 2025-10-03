program multiwfn
use defvar
use util
use GUI
use mod_2F2, only: set_alpha_level
!USE fparser
!use libreta
!use function
implicit real*8 (a-h,o-z)
character nowdate*20,nowtime*20,c200tmp*200,c2000tmp*2000,lovername*80,settingpath*200,strtmp*3,selectyn
real*8,allocatable :: tmparr(:),tmparr2(:),tmpmat(:,:),tmpmat2(:,:),tmpmat3D(:,:,:) !For debug purpose
integer,allocatable :: tmparri(:),tmparr2i(:),tmpmati(:,:),tmpmat2i(:,:)
real*8 tmpv1(3),tmpv2(3)
!Special treatment for Intel compiler
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
call kmp_set_warnings_off() !In rare case, "Cannot open message catalog "1041\libiomp5ui.dll"" may occurs, this calling avoid this problem, or user should set KMP_WARNINGS environment variable to 0
#endif

!Try to get input file name from argument, which should be the first argument
filename=" "
narg=command_argument_count()
if (narg>0) then
    call getarg(1,filename)
    inquire(file=filename,exist=alive)
    if (.not.alive) then
        write(*,*) "Error: Unable to find the input file you specified in argument"
        write(*,*)
        filename=" "
    end if
end if

10 call loadsetting
write(*,*) "Multiwfn -- A Multifunctional Wavefunction Analyzer"
write(*,*) "Version 3.8(dev), update date: 2025-Oct-2"
write(*,*) "Developer: Tian Lu (Beijing Kein Research Center for Natural Sciences)"
write(*,*) "Multiwfn official website: http://sobereva.com/multiwfn"
write(*,*) "Multiwfn English forum: http://sobereva.com/wfnbbs"
write(*,*) "Multiwfn Chinese forum: http://bbs.keinsci.com/wfn"

call date_and_time(nowdate,nowtime)
write(*,"(' ( Number of parallel threads:',i4,'  Current date: ',a,'-',a,'-',a,'  Time: ',a,':',a,':',a,' )')") &
nthreads,nowdate(1:4),nowdate(5:6),nowdate(7:8),nowtime(1:2),nowtime(3:4),nowtime(5:6)
write(*,*)
write(*,*) "Both following papers ***MUST BE CITED IN MAIN TEXT*** if Multiwfn is used:"
write(*,*) " Tian Lu, Feiwu Chen, J. Comput. Chem., 33, 580 (2012) DOI: 10.1002/jcc.22885"
write(*,*) " Tian Lu, J. Chem. Phys., 161, 082503 (2024) DOI: 10.1063/5.0216272"
write(*,*) "See ""How to cite Multiwfn.pdf"" in Multiwfn binary package for more information"


!!-------- Set up hardware resource information
!For Windows version of ifort, use KMP_SET_STACKSIZE_S() to directly set stacksize of OpenMP threads according to ompstacksize in settings.ini, &
!for other cases, the stacksize is determined by OMP_STACKSIZE environment variable, and we check if it has been defined here
!  ompstacksize variable is always actually used by Multiwfn to determine available memory
if (isys==1) then !Windows
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
    call KMP_SET_STACKSIZE_S(ompstacksize)
#else
    CALL getenv('OMP_STACKSIZE',c200tmp)
    if (c200tmp==" ") write(*,"(/,a)") " Warning: You should set OMP_STACKSIZE environment variable in Windows system to define stacksize of OpenMP threads!"
#endif
else if (isys==2) then !Linux/MacOS
    CALL getenv('OMP_STACKSIZE',c200tmp)
    if (c200tmp/=" ") call read_ompstacksize(c200tmp)
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
    if (c200tmp==" ") then
		CALL getenv('KMP_STACKSIZE',c200tmp) !OpenMP stacksize may also be defined by KMP_STACKSIZE by users if compiled with ifort
		if (c200tmp/=" ") call read_ompstacksize(c200tmp)
    end if
#endif
    if (c200tmp==" ") write(*,"(/,a)") " Warning: You should set OMP_STACKSIZE environment variable as mentioned in Section 2.1.2 of Multiwfn manual!"
end if
!write(*,"(' OpenMP stacksize for each thread: ',f10.2,' MB')") dfloat(ompstacksize)/1024/1024 !Use this if you want to monitor
!write(*,"(' OpenMP stacksize for each thread: ',f10.2,' MB')") dfloat(KMP_GET_STACKSIZE_S())/1024/1024 !Extended subroutine by Intel compiler

!Set number of cores used by calculation of MKL library (e.g. function matmul_blas)
#if defined(INTEL_MKL)
call mkl_set_num_threads(nthreads)
#endif


!!-------- Load input file
write(*,*)
if (trim(filename)==" ") then !Haven't defined filename variable
	call mylover(lovername)
	write(*,"(a,a,a)") " Now input file path, for example, E:\",trim(lovername),".mwfn"
	write(*,*) "(.wfn/wfn/wfx/fch/molden/pdb/xyz/mol2/cif/cub... see Section 2.5 of manual)"
	write(*,"(a)") " Hint: Pressing ENTER button directly can select a file in a GUI window. To reload the past file, inputting ""o"". &
	&Input such as ?miku.fch can open the miku.fch in the same folder as the past file"
	do while(.true.)
		read(*,"(a)") filename
		if (filename=='o') then
			write(*,"(' The file last time used: ',a)") trim(lastfile)
			filename=lastfile
		else if (filename==' ') then
			call selfileGUI
			if (filename==' ') then
				write(*,*) "You did not select a file, input path again!"
				cycle !User didn't select a file
			end if
			write(*,"(' Selected file: ',a)") trim(adjustl(filename))
		end if
        call remove_quotemark(filename) !Remove possible " or ' symbol
		if (filename(1:1)=='?') then
			do itmp=len_trim(lastfile),1,-1
				if (isys==1.and.lastfile(itmp:itmp)=='\') exit
				if (isys==2.and.lastfile(itmp:itmp)=='/') exit
			end do
			filename=lastfile(1:itmp)//trim(filename(2:))
		end if
		inquire(file=filename,exist=alive)
		if (alive.eqv..true.) exit
		write(*,"('""',a,'"" ',a)") trim(filename),"cannot be found, input again"
	end do
	!Write current opened file to "lastfile" in settings.ini
	inquire(file="settings.ini",exist=alive)
	if (alive) then
		settingpath="settings.ini"
	else if (.not.alive) then
		call getenv("Multiwfnpath",c200tmp)
		if (isys==1) then
			settingpath=trim(c200tmp)//"\settings.ini"
		else if (isys==2) then
			settingpath=trim(c200tmp)//"/settings.ini"
		end if
	end if
	inquire(file=settingpath,exist=alive)
	if (alive) then
		open(20,file=settingpath,status="old")
		call loclabel(20,"lastfile")
		write(20,"(a)") "lastfile= """//trim(filename)//""""
		close(20)
	end if
else
	inquire(file=filename,exist=alive)
	if (alive.eqv..false.) then
		write(*,*) "Error: File not found, exit program..."
		read(*,*)
		stop
	end if
end if

call readinfile(filename,0)
write(*,"(/,3a)") " Loaded ",trim(filename)," successfully!"


!!-------- Backup various information of first loaded (meanwhile unmodified) system
firstfilename=filename
if (allocated(a)) then
	allocate(a_org(ncenter))
	a_org=a
	ncenter_org=ncenter
end if
if (allocated(b)) then
	allocate(b_org(nprims))
	allocate(CO_org(nmo,nprims))
	allocate(MOocc_org(nmo))
	allocate(MOene_org(nmo))
	b_org=b
	CO_org=CO
	MOocc_org=MOocc
	MOene_org=MOene
	nprims_org=nprims
	nmo_org=nmo
end if
ifPBC_org=ifPBC
cellv1_org=cellv1
cellv2_org=cellv2
cellv3_org=cellv3


!!-------- Generate backed up fragment
nfragatm_org=ncenter
allocate(fragatm_org(nfragatm_org))
forall (i=1:nfragatm_org) fragatm_org(i)=i


!!-------- Call some routines only once
atmrhocutsqr(:)=atmrhocut(:)**2
atmrhocutsqr_1En5(:)=atmrhocut_1En5(:)**2
!Generate coordinate of atomic radial positions, which corresponds to atmraddens(:)
call genatmradpos(atmradpos(:))
!Convert prebuilt radii from Angstrom to Bohr. But some radii such as radii_hugo will remain unchanged since it is recorded as Bohr
if (ifirstMultiwfn==1) then
	vdwr=vdwr/b2a
	vdwr_tianlu=vdwr_tianlu/b2a
	vdwr_UFF=vdwr_UFF/b2a
	covr=covr/b2a
	covr_Suresh=covr_Suresh/b2a
	covr_pyy=covr_pyy/b2a
	covr_tianlu=covr_tianlu/b2a
end if
if (isys==2) call swgfnt("fixed",12) !In Linux, if font is not set for widget, the texts may not be displayed properly. Text size (12) is meaningless, since the size is bonded to font.
!Only get into dislin level 1 to get width and height of screen in pixels, don't do any other things
call METAFL("GKSL")
call disini
CALL GETSCR(iscrwidth,iscrheight)
call ERRMOD("ALL","OFF")
call disfin
CALL SYMFIL("NONE","DELETE")
if (itransparent==1) then
    CALL PNGMOD("ON",'TRANSPARENCY')
    CALL GIFMOD("ON",'TRANSPARENCY')
end if


!!-------- Show cell information
call showcellinfo
call getcellabc(asize,bsize,csize,alpha,beta,gamma)
crit=7D0/b2a
if (allocated(b)) then !Warning user if any direction is less than 7 Angstrom
	if (asize<crit.and.PBCnx==1) write(*,"(/,a,i3)") " Warning: Because size of a axis of the cell is relatively small, &
	&in order to guarantee analysis accuracy, it is suggested to increase the first index of ""PBCnxnynz"" in settings.ini to",ceiling(crit/asize)
	if (bsize<crit.and.PBCny==1) write(*,"(/,a,i3)") " Warning: Because size of b axis of the cell is relatively small, &
	&in order to guarantee analysis accuracy, it is suggested to increase the second index of ""PBCnxnynz"" in settings.ini to",ceiling(crit/bsize)
	if (csize<crit.and.PBCnz==1) write(*,"(/,a,i3)") " Warning: Because size of c axis of the cell is relatively small, &
	&in order to guarantee analysis accuracy, it is suggested to increase the third index of ""PBCnxnynz"" in settings.ini to",ceiling(crit/csize)
end if


!!-------- Show brief information of present system
if (allocated(a)) then
	write(*,*)
	call showformula
	totmass=sum(atmwei(a%index))
	write(*,"(' Molecule weight:',f16.5,' Da')") totmass
    !-- Show point group
    if (ifPBC==0.and.ncenter<nPGmaxatm.and.all(a%index>0)) then !Too large system will take evidently cost
        allocate(tmpmat(3,ncenter),tmpmat2i(ncenter,ncenter),tmparri(ncenter))
        tmpmat(1,:)=a%x*b2a;tmpmat(2,:)=a%y*b2a;tmpmat(3,:)=a%z*b2a
        !Tolerance of 0.0025 is suitable for most systems. Though it may be too tight to detect point group for difficult case, &
        !if the criterion is loosen, in rare case (e.g. C32) the SYVA routine will ceaselessly show &
        !"ERROR: Too many symmetry operations. Try a lower tolerance" and doesn't work or even make Multiwfn crash
        call PG_eqvatm(ncenter,a%index,tmpmat,0.0025D0,strtmp,ncls,tmparri,tmpmat2i)
        if (strtmp==" ".and.ncenter<50) then
            do i=1,20
                !write(*,"(' Loosen the tolerance of determining point group to',f12.6)") i*0.003D0
                call PG_eqvatm(ncenter,a%index,tmpmat,i*0.003D0,strtmp,ncls,tmparri,tmpmat2i)
                if (strtmp/=" ") exit
            end do
        end if
        if (strtmp==" ") then
            write(*,*) "Failed to detect point group"
        else
            write(*,"(' Point group: ',a)") strtmp
        end if
        !Test ability of detecting point group
        !do i=1,500
        !    val=i*0.00002D0
        !    call PG_eqvatm(ncenter,a%index,tmpmat,val,strtmp,ncls,tmparri,tmpmat2i)
        !    write(*,"(i6,f12.6,2x,a)") i,val,trim(strtmp)
        !end do
        deallocate(tmpmat,tmpmat2i,tmparri)
    end if
end if



!Special treatment and test new code
!call ORCAjson_gennatorb

!!!--------------------- Now everything start ---------------------!!!
do while(.true.) !Main loop

	write(*,*)
	if (allocated(cubmat)) write(*,*) "Note: A set of grid data is presented in memory"
    write(*,*) """q"": Exit program gracefully          ""r"": Load a new file"
	write(*,*) "                   ************ Main function menu ************"
	if (ifiletype==7) then
        write(*,*) "0 Show molecular structure and view isosurface"
    else if (ifiletype==8) then
        write(*,*) "0 View isosurface"
    else
        write(*,*) "0 Show molecular structure and view orbitals"
    end if
	write(*,*) "1 Output all properties at a point       2 Topology analysis"
	write(*,*) "3 Output and plot specific property in a line"
	write(*,*) "4 Output and plot specific property in a plane"
	write(*,*) "5 Output and plot specific property within a spatial region (calc. grid data)"
	write(*,*) "6 Check & modify wavefunction"
	write(*,*) "7 Population analysis and calculation of atomic charges"
	write(*,*) "8 Orbital composition analysis           9 Bond order analysis"
	write(*,*) "10 Plot total DOS, PDOS, OPDOS, local DOS, COHP and photoelectron spectrum"
	write(*,*) "11 Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum"
	write(*,*) "12 Quantitative analysis of molecular surface"
	if (allocated(cubmat)) write(*,*) "13 Process grid data"
	if (.not.allocated(cubmat)) write(*,*) "13 Process grid data (No grid data is presented currently)"
	write(*,*) "14 Adaptive natural density partitioning (AdNDP) analysis"
	write(*,*) "15 Fuzzy atomic space analysis"
	write(*,*) "16 Charge decomposition analysis (CDA) and plot orbital interaction diagram"
	write(*,*) "17 Basin analysis                       18 Electron excitation analysis"
	write(*,*) "19 Orbital localization analysis        20 Visual study of weak interaction"
	write(*,*) "21 Energy decomposition analysis        22 Conceptual DFT (CDFT) analysis"
    write(*,*) "23 ETS-NOCV analysis                    24 (Hyper)polarizability analysis"
    write(*,*) "25 Electron delocalization and aromaticity analyses"
    write(*,*) "26 Structure and geometry related analyses"
	write(*,*) "100 Other functions (Part 1)            200 Other functions (Part 2)"
	write(*,*) "300 Other functions (Part 3)"
	! write(*,*) "1000 Special functions"
	read(*,*) c200tmp
    
    if (c200tmp=="q".or.c200tmp=="-10") then !Exit program
        stop
	else if (c200tmp=="r".or.c200tmp=="-11") then !Load a new file
	    call dealloall(0)
	    call dealloall_org
	    filename=""
	    ifirstMultiwfn=0
	    goto 10
    else if (c200tmp=="oi") then
	    call outORCAinp_wrapper
    else if (c200tmp=="gi") then
        call outgjf_wrapper
    else if (c200tmp=="pi") then
        call outPSI4inp_wrapper
    else if (c200tmp=="cp") then
        call cp2kmate
    else if (c200tmp=="cp2k") then
        call outCP2Kinp_wrapper
    else if (c200tmp=="QE") then
        call outQEinp_wrapper
    else if (c200tmp=="iu") then
        write(*,*) "Input the index of the user-defined function you want to use, e.g. 5"
        read(*,*) iuserfunc
        call init_func
        write(*,*) "Done!"
    else if (c200tmp=="geomparm") then
		iallowPBC=0
		if (ifPBC>0) then
			write(*,*) "Consider periodic boundary condition? (y/n)"
            read(*,*) selectyn
            if (index(selectyn,'y')/=0) iallowPBC=1
        end if
		write(*,*) "Input path of the file to be exported, e.g. /sob/geomparm.txt"
        read(*,"(a)") c200tmp
		call showgeomparam(c200tmp,iallowPBC)
        write(*,"(a)") " Done! The internal coordinates have been exported to "//trim(c200tmp)
    else if (c200tmp=="MPP".or.c200tmp=="mpp") then
		call calcMPP
    else if (c200tmp=="cav") then
		call cavity_diameter
    else if (c200tmp=="pdb") then
	    call outpdb_wrapper
    else if (c200tmp=="xyz") then
	    call outxyz_wrapper
    else if (c200tmp=="cif") then
	    call outcif_wrapper
    else
        read(c200tmp,*,iostat=ierror) isel
        if (ierror/=0) then
			write(*,*) "Error: Unable to recognize the command you inputted"
            cycle
        end if

	    !!!---------------------------------------------------------------------------------------------
	    !1!------- Show system structure and view isosurface of MOs or the grid data read from cube file
	    if (isel==0) then
		    if (.not.(allocated(a).or.allocated(cubmat))) then
			    write(*,*) "Error: Data needed by this function is not presented! Check your input file!"
			    write(*,*) "Press ENTER button to continue"
			    read(*,*)
			    cycle
		    end if
		    if (ncenter>0) then
				if (ncenter<=300) then
					write(*,*) "Atom list:"
					do i=1,ncenter
						write(*,"(i5,'(',a2,')',' --> Charge:',f10.6,'  x,y,z(Bohr):',3f11.6)") i,a(i)%name,a(i)%charge,a(i)%x,a(i)%y,a(i)%z
					end do
                else
					write(*,"(a)") " Note: There are more than 300 atoms, so their information is not shown here now. &
                    &To print, in the manu bar please select ""Tools"" - ""Print XYZ coordinates"""
                    ishowatmlab=0
				end if
            end if
		    if (allocated(CObasa).and.imodwfn==0) then !fch and occupation number hasn't been modified
			    if (wfntype==0) then
				    write(*,"(' Note: Orbital',i6,' is HOMO, energy:',f12.6,' a.u.',f12.6,' eV')") nint(nelec/2),MOene(nint(nelec/2)),MOene(nint(nelec/2))*au2eV
				    if (nint(nelec/2)+1<=nmo) then
					    write(*,"('       Orbital',i6,' is LUMO, energy:',f12.6' a.u.',f12.6,' eV')") nint(nelec/2)+1,MOene(nint(nelec/2)+1),MOene(nint(nelec/2)+1)*au2eV
					    gapene=MOene(nint(nelec/2)+1)-MOene(nint(nelec/2))
					    write(*,"('       HOMO-LUMO gap:',f12.6,' a.u.',f12.6,' eV',f14.6,' kJ/mol')") gapene,gapene*au2eV,gapene*au2kJ
				    end if
			    else if (wfntype==1) then
				    write(*,"(' Range of alpha orbitals:',i5,' -',i5,'      Range of beta orbitals:',i5,' -',i5)") 1,nbasis,nbasis+1,nmo
				    write(*,"(' Note: Orbital',i6,' is alpha-HOMO, energy:',f12.6,' a.u.',f12.6,' eV')") nint(naelec),MOene(nint(naelec)),MOene(nint(naelec))*au2eV
				    write(*,"('       Orbital',i6,' is beta-HOMO, energy: ',f12.6,' a.u.',f12.6,' eV')") nbasis+nint(nbelec),MOene(nbasis+nint(nbelec)),MOene(nbasis+nint(nbelec))*au2eV
				    if (nbasis>=nint(naelec)+1) then
					    write(*,"('       Orbital',i6,' is alpha-LUMO, energy:',f12.6,' a.u.',f12.6,' eV')") nint(naelec)+1,MOene(nint(naelec)+1),MOene(nint(naelec)+1)*au2eV
					    write(*,"('       Orbital',i6,' is beta-LUMO, energy: ',f12.6,' a.u.',f12.6,' eV')") nbasis+nint(nbelec)+1,MOene(nbasis+nint(nbelec)+1),MOene(nbasis+nint(nbelec)+1)*au2eV
					    gapenea=MOene(nint(naelec)+1)-MOene(nint(naelec))
					    write(*,"('       HOMO-LUMO gap of alpha orbitals:',f12.6,' a.u.',f12.6,' eV')") gapenea,gapenea*au2eV
					    gapeneb=MOene(nbasis+nint(nbelec)+1)-MOene(nbasis+nint(nbelec))
					    write(*,"('       HOMO-LUMO gap of beta orbitals: ',f12.6,' a.u.',f12.6,' eV')") gapeneb,gapeneb*au2eV
				    end if
			    else if (wfntype==2) then
				    write(*,"(' Index of SOMO orbitals:',10i6)") (i,i=nint(nbelec+1),nint(naelec))
			    end if
		    end if
		    if (ifiletype==7.or.ifiletype==8) then !Visualize grid data
			    if (isilent==0) call drawisosurgui(1)
		    else
			    if (isilent==0) call drawmolgui
		    end if
            iorbvis=0 !Recover its status. iorbvis=0 makes saved image file has DISLIN prefix
            ishowdatarange=0 !Do not show data range if showing it
            call setfil("dislin."//trim(graphformat)) !The file name of saved image file may have been modified (e.g. equal to orbital index), so recover to default

            
	    !!!-------------------------------------------------
	    !1!-------------------- Output properties at a point
	    else if (isel==1) then
		    call study0dim


	    !!!--------------------------------------
	    !2!-------------------- Topology analysis
	    else if (isel==2) then
		    call topo_main


	    !!!--------------------------------------------
	    !3!-------------------- Draw property in a line
	    else if (isel==3) then
		    call study1dim


	    !!!-------------------------------------
	    !4!-------------------- Draw plane graph
	    else if (isel==4) then
		    call study2dim(0,0,0)

            
	    !!!--------------------------------------------------------
	    !5!------------------- Calculate, show and output grid file
	    else if (isel==5) then
		    call study3dim
            

	    !!!--------------------------------------------------------------------------------
	    !6!!------------------- Check & Modify wavefunction or show GTF/Orbital information
	    else if (isel==6) then
		    call modwfn


	    !!!---------------------------------------
	    !7!!------------------- Population analysis
	    else if (isel==7) then
		    call population_main


	    !!!-------------------------------------------------
	    !8!!------------------- Orbital composition analysis
	    else if (isel==8) then
		    call orbcomp_main


	    !!!----------------------------------------
	    !9!!------------------- Bond order analysis
	    else if (isel==9) then
		    call bondorder_main


	    !!!------------------------------
	    !10!!------------------- Plot DOS
	    else if (isel==10) then
		    call DOS
		
	    !!!----------------------------------
	    !11!!------------------- Plot spectra
	    else if (isel==11) then
		    call plotspectrum


	    !!!------------------------------------------------
	    !12!!------------------- Molecular surface analysis
	    else if (isel==12) then
		    call surfana

		
	    !!!---------------------------------------
	    !13!!------------------- Process grid data
	    else if (isel==13) then
		    call procgriddata


	    !!!-------------------------------------------------------------------
	    !14!!------------------- Adaptive natural density partitioning (AdNDP)
	    else if (isel==14) then
		    call AdNDP
		
		
	    !!!--------------------------------------------------
	    !15!!------------------- Integrate fuzzy atomic space
	    else if (isel==15) then
		    call fuzzyana(0)
	    else if (isel==-15) then
		    call fuzzySBL	


	    !!!---------------------------------------------------
	    !16!!------------------- Charge decomposition analysis
	    else if (isel==16) then
		    call CDA


	    !!!---------------------------------------
	    !17!!------------------- Basin integration
	    else if (isel==17) then
		    call basinana


	    !!!--------------------------------------------------
	    !18!!------------------- Electron excitation analysis
	    else if (isel==18) then
		    call excittrans_main


	    !!!---------------------------------------------------
	    !19!!------------------- Orbital localization analysis
	    else if (isel==19) then
		    call orbloc
		
		
	    !!!------------------------------------------------------
	    !20!!------------------- Visual study of weak interaction
	    else if (isel==20) then
		    call visweak_main
		
		
	    !!!---------------------------------------------------
	    !21!!------------------- Energy decomposition analysis
	    else if (isel==21) then
		    call EDA_main
		
		
	    !!!---------------------------------------------
	    !22!!------------------- Conceptual DFT analysis
	    else if (isel==22) then
		    call CDFT
		
		
	    !!!-----------------------------------
	    !23!!------------------- NOCV analysis
	    else if (isel==23) then
		    call ETS_NOCV
		
		
	    !!!----------------------------------------------------
	    !24!!------------------- (Hyper)polarizability analyses
	    else if (isel==24) then
		    call hyper_polarizability

            
	    !!!-------------------------------------------------------------
	    !25!!------------------- Delocalization and aromaticity analyses
	    else if (isel==25) then
		    call deloc_aromat

            
	    !!!-------------------------------------------------------------
	    !26!!------------------- Structure and geometry related analyses
	    else if (isel==26) then
			do while(.true.)
				write(*,*)
				call menutitle("Geometry related analysis",10,1)
				write(*,*) "0 Return"
				write(*,*) "1 Calculate properties based on geometry information for specific atoms"
				write(*,*) "2 Various geometry operation on present system"
				write(*,"(a)") " 3 Molecular planarity parameter (MPP) and span of deviation from plane (SDP)"
				write(*,*) "4 Evaluate interatomic connectivity and atomic coordination number"
				write(*,*) "5 Calculate average bond length and average coordinate number"
				write(*,*) "6 Calculate bond length/order alternation (BLA/BOA)"
				write(*,*) "7 Calculate kinetic diameter"
				write(*,*) "8 Calculate molecular diameter and length/width/height"
				write(*,*) "9 Visualize molecular cavity and calculate its volume"
				write(*,*) "10 Plot surface distance projection map"
				write(*,*) "11 View free regions and calculating free volume in a cell"
				write(*,*) "12 Calculate area of vdW surface of whole system or individual fragments"
				write(*,*) "13 Calculate vdW volume via Monte Carlo method"
				write(*,*) "14 Calculate vdW volume via Marching Tetrahedron algorithm based on electron density isosurface"
				write(*,*) "15 Calculate cavity diameter of molecule and crystal, as well as graphically representing it as sphere"
				write(*,*) "16 Calculate area and perimeter of a specific ring"
				write(*,"(a)") " 17 Calculate minimum/maximum and geometry/mass center distances between two fragments"
				read(*,*) isel
				if (isel==0) then
					exit
				else if (isel==1) then
					call calcgeomprop
				else if (isel==2) then
					call geom_operation
				else if (isel==3) then
					write(*,"(a,/)") " Hint: You can also directly enter this function by inputting ""MPP"" in main menu of Multiwfn"
					call calcMPP
				else if (isel==4) then
					call conn_coordnum
				else if (isel==5) then
					call atmavgdist
				else if (isel==6) then
					call BLABOA
				else if (isel==7) then
					write(*,"(a)") " Please check Section 4.12.12 of Multiwfn manual on how to use main function 12 to calculate kinetic diameter"
					write(*,*) "Press ENTER button to return"
					read(*,*)
				else if (isel==8) then
					call calcmolsize
				else if (isel==9) then
					write(*,"(a)") " Please check Section 4.200.14.2 of Multiwfn manual on how to use domain analysis to visualize molecular cavity and calculate its volume"
					write(*,*) "Press ENTER button to return"
					read(*,*)
				else if (isel==10) then
					call molsurf_distmap
				else if (isel==11) then
					call freeregion
				else if (isel==12) then
					write(*,"(a)") " Please check Section 4.12.9 of Multiwfn manual on how to use main function 12 to calculate area of van der Waals surface of whole system or individual fragments"
					write(*,*) "Press ENTER button to return"
					read(*,*)
				else if (isel==13) then
					call molvol_MC
				else if (isel==14) then
					write(*,"(a)") " Please check Section 4.12.1 of Multiwfn manual on how to use main function 12 to &
					&calculate van der Waals volume based on electron density isosurface"
					write(*,*) "Briefly speaking, you should simply enter main function 12 and choose option 6 to obtain the volume"
					write(*,*) "Press ENTER button to return"
					read(*,*)
				else if (isel==15) then
					call cavity_diameter
				else if (isel==16) then
					call calcringsize
				else if (isel==17) then
					call calcfragdist
				end if
			end do
		

	    !!!------------------------------------------------------------------
	    !100!!------------------- Misc and some unimportant functions, Part 1
	    else if (isel==100) then
		    call otherfunc_main
		
		
	    !!!------------------------------------------------------------------
	    !200!!------------------- Misc and some unimportant functions, Part 2
	    else if (isel==200) then
		    call otherfunc2_main
		
		
	    !!!------------------------------------------------------------------
	    !200!!------------------- Misc and some unimportant functions, Part 3
	    else if (isel==300) then
		    call otherfunc3_main


	    !!!-----------------------------------------
	    !1000!!------------------- Special functions
	    else if (isel==1000) then
		    write(*,*)
            write(*,*) " ---------------------------- Special functions ----------------------------"
		    write(*,*) "0 Return to main menu"
		    write(*,"(a,3f12.6,' Bohr')") " 1 Set reference point, current:",refx,refy,refz
		    write(*,"(a,i5,a,i5)") " 2 Set iuserfunc, current:",iuserfunc,"            3 Set iskipnuc, current:",iskipnuc
		    if (pleA==0D0.and.pleB==0D0.and.pleC==0D0.and.pleD==0D0) then
			    write(*,"(a)") " 4 Set the plane for user-defined function 38 (Not defined)"
		    else
			    write(*,"(a)") " 4 Set the plane for user-defined function 38 (Defined)"
		    end if
		    write(*,"(a,1PE18.8)") " 5 Set global temporary variable, current:",globaltmp
		    write(*,"(a,i3)") " 10 Set number of threads of parallel calculation, current:",nthreads
		    write(*,*) "11 Reload settings.ini file"
            write(*,*) "12 Add a Bq atom to specific position"
            write(*,*) "13 Convert bndmat.txt in current folder to Gaussian .gjf file with bond orders"
            write(*,*) "14 Convert current wavefunction to .rad file"
            write(*,*) "15 Make orbitals equivalent to basis functions"
            write(*,*) "-15 Make orbitals equivalent to Lowdin orthogonalized basis functions"
		    write(*,*) "16 Define one or two fragments for special purpose"
            write(*,*) "17 Generate promolecular wavefunction by calculating and combining atomic ones"
            if (allocated(cubmat)) write(*,*) "18 Set box information of grid data as cell information"
            write(*,*) "88/89 Calculate two-electron integral between for four PGTFs/orbitals"
		    write(*,*) "90 Calculate nuclear attractive energy between a fragment and an orbital"
		    write(*,*) "91 Exchange orbital energies and occupations"
		    write(*,*) "92 Calculate result of various kinetic energy functionals"
            write(*,*) "93 Output all Becke's integration points to intpt.txt in current folder"
		    write(*,*) "97 Generate natural orbitals based on density matrix outputted by MRCC program"
		    write(*,*) "98 Generate natural orbitals based on density matrix outputted by ORCA program"
		    write(*,*) "99 Show EDF information (if any)"
		    write(*,*) "100 Check the sanity of present wavefunction"
            write(*,*) "201 Ring-ring distance&angle statistical analysis for a trajectory" !Only used by Sobereva in C18 work
            write(*,*) "202 Ring rotation statistical analysis for a trajectory" !Only used by Sobereva in C18 work
            write(*,*) "1303 Setup fractional calculus"
		    read(*,*) i
		    if (i==1) then
			    write(*,*) "Input x,y,z in Bohr, e.g. 3.0,0.0,1.3"
                write(*,*) "Note: with ""A"" suffix, the unit will be in Angstrom, e.g. 2,0.12,3.5 A"
                read(*,"(a)") c200tmp
			    read(c200tmp,*) refx,refy,refz
                if (index(c200tmp,'A')/=0) then
					refx=refx/b2a
                    refy=refy/b2a
                    refz=refz/b2a
                end if
			    write(*,*) "Done!"
		    else if (i==2) then
			    write(*,*) "Input index of the user-defined function, e.g. 24"
			    read(*,*) iuserfunc
				call init_func
			    write(*,*) "Done!"
		    else if (i==3) then
			    write(*,*) "Input the index of the nucleus, e.g. 24"
			    read(*,*) iskipnuc
			    write(*,*) "Done!"
		    else if (i==4) then
			    write(*,*) "1 Input index of three atoms to define the plane"
			    write(*,*) "2 Input XYZ coordinate of three points to define the plane"
			    read(*,*) iseldef
			    if (iseldef==1) then
				    write(*,*) "Input three indices, e.g. 2,4,5"
				    read(*,*) i1,i2,i3
				    call pointABCD(a(i1)%x,a(i1)%y,a(i1)%z,a(i2)%x,a(i2)%y,a(i2)%z,a(i3)%x,a(i3)%y,a(i3)%z,pleA,pleB,pleC,pleD)
			    else if (iseldef==2) then
				    write(*,*) "Input coordinate for point 1 (in Bohr), e.g. 1.0,-0.2,0.3"
				    read(*,*) xtmp1,ytmp1,ztmp1
				    write(*,*) "Input coordinate for point 2 (in Bohr), e.g. 2.0,-0.3,0.1"
				    read(*,*) xtmp2,ytmp2,ztmp2
				    write(*,*) "Input coordinate for point 3 (in Bohr), e.g. 1.3,-1.2,0.33"
				    read(*,*) xtmp3,ytmp3,ztmp3
				    call pointABCD(xtmp1,ytmp1,ztmp1,xtmp2,ytmp2,ztmp2,xtmp3,ytmp3,ztmp3,pleA,pleB,pleC,pleD)
			    end if
			    tmpval=dsqrt(pleA**2+pleB**2+pleC**2)
			    write(*,"(' The unit vector normal to the plane is:',3f10.5)") pleA/tmpval,pleB/tmpval,pleC/tmpval
		    else if (i==5) then
			    write(*,*) "Input the value, e.g. 0.3"
			    read(*,*) globaltmp
			    write(*,*) "Done!"
		    else if (i==10) then
			    write(*,*) "Input an integer, e.g. 8"
			    read(*,*) nthreads
			    write(*,*) "Done!"
		    else if (i==11) then
			    call loadsetting
			    write(*,*) "Done!"
		    else if (i==12) then
                do while(.true.)
                    write(*,*) "Write the X,Y,Z of the Bq atom to be added in Bohr, e.g. 0.2,0,-3.5"
                    write(*,*) "Input ""q"" can exit"
                    read(*,"(a)") c200tmp
                    if (c200tmp=="q") then
                        exit
                    else
                        read(c200tmp,*) tmpx,tmpy,tmpz
                        call addbq(tmpx,tmpy,tmpz)
			            write(*,*) "Done!"
                    end if
                end do
		    else if (i==13) then
                allocate(tmpmat(ncenter,ncenter))
                inquire(file="bndmat.txt",exist=alive)
	            if (.not.alive) then
	                write(*,*) "Error: Cannot find the bndmat.txt in current folder!"
                    write(*,*) "Press ENTER button to return"
                    read(*,*)
                    cycle
                end if
                open(10,file="bndmat.txt",status="old")
                call loclabel(10,"Bond order matrix for all electron",ifound)
                if (ifound==1) then
					write(*,"(a)") " Note: Found ""Bond order matrix for all electron"" section, this part of bond order will be loaded"
                else
					rewind(10)
					call loclabel(10,"***")
                end if
                call readmatgau(10,tmpmat,0,"f14.8",6,5)
                close(10)
                open(10,file="gau.gjf",status="replace")
                write(10,"(a,/,/,a,/)") "#P B3LYP/6-31G* geom=connectivity","Generated by Multiwfn"
                netcharge=nint(sum(a%charge)-nelec)
                if (nelec==0) netcharge=0 !nelec==0 means no electron information, e.g. pdb file
                write(10,"(2i3)") netcharge,nint(naelec-nbelec)+1
                do i=1,ncenter
	                write(10,"(a,1x,3f14.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
                end do
                write(10,*)
                if (.not.allocated(connmat)) call genconnmat(1,0) !Generate connectivity matrix
                do iatm=1,ncenter
                    write(10,"(i8)",advance="no") iatm
                    do jatm=iatm+1,ncenter
                        if (connmat(iatm,jatm)==0) cycle
                        write(10,"(i8,f8.4)",advance="no") jatm,tmpmat(iatm,jatm)
                    end do
                    write(10,*)
                end do
                close(10)
                deallocate(tmpmat)
                write(*,*) "Done! gau.gjf has been generated in current folder"
            else if (i==14) then
                ipos=index(filename,'.',back=.true.)
                write(*,"(a)") " Converting to "//filename(1:ipos)//"rad"
                call atmwfn2atmrad(filename,filename(1:ipos)//"rad")
                write(*,*) "Done!"
            else if (abs(i)==15) then
                if (.not.allocated(CObasa)) then
                    write(*,*) "Error: The input file must contain basis function information!"
                    write(*,*) "Press ENTER button to return"
                    read(*,*)
                    cycle
                end if
                write(*,*) "Please wait..."
                CObasa=0
                do ibas=1,nbasis
                    CObasa(ibas,ibas)=1
                end do
				if (i==-15) then
					call ask_Sbas_PBC
					call symmortho(0)
                end if
                call CObas2CO(1)
                if (wfntype==0.or.wfntype==2.or.wfntype==3) then
                    write(*,"(a)") " Done! Now each orbital corresponds to a basis function, index of orbitals is identical to index of basis functions"
                else
                    write(*,"(a)") " Done! Now each alpha orbital corresponds to a basis function, index of alpha orbitals is identical to index of basis functions"
                end if
            else if (i==16) then
                if (allocated(frag1)) deallocate(frag1)
                if (allocated(frag2)) deallocate(frag2)
                write(*,*) "1 Define one fragment"
                write(*,*) "2 Define two fragments"
                read(*,*) ndeffrag
                do ifrag=1,ndeffrag
                    write(*,"(/,a,i2,a)") " Input indices of the atoms in fragment",ifrag,", e.g. 1-5,13,19"
                    read(*,"(a)") c2000tmp
                    if (ifrag==1) then
                        call str2arr(c2000tmp,nfrag1)
                        allocate(frag1(nfrag1))
                        call str2arr(c2000tmp,ntmp,frag1)
                    else if (ifrag==2) then
                        call str2arr(c2000tmp,nfrag2)
                        allocate(frag2(nfrag2))
                        call str2arr(c2000tmp,ntmp,frag2)
                    end if
                    write(*,"(' Number of inputted atoms in this fragment:',i6)") ntmp
                end do
            else if (i==17) then
				call generate_promolwfn(1)
            else if (i==18) then
				call grid2cellinfo
                write(*,*) "Done!"
                call showcellinfo
            else if (i==88) then
				call showGTF_ERI
            else if (i==89) then
				call showorb_ERI
		    else if (i==90) then
			    call attene_orb_fragnuc
		    else if (i==91) then
			    do iorb=1,nmo
				    tmp=MOocc(iorb)
				    MOocc(iorb)=MOene(iorb)
				    MOene(iorb)=tmp
			    end do
			    imodwfn=1
			    write(*,*) "Done!"
			    if (allocated(CObasa)) then
				    write(*,*) "Updating density matrix..."
				    call genP
				    write(*,*) "Density matrix has been updated"
			    end if
		    else if (i==92) then
			    call intKED
		    else if (i==93) then
                call outBeckeintpt
		    else if (i==97) then
			    call MRCC_gennatorb
		    else if (i==98) then
				call ORCAjson_gennatorb
		    else if (i==99) then
			    if (nEDFprims==0) then
				    write(*,*) "EDF field was not loaded"
			    else
				    write(*,"( ' The number of inner-core electrons represented by EDFs:',i6)") nEDFelec
				    write(*,*) "Information of EDF primitives:"
				    write(*,*) "Column 1: Index"
				    write(*,*) "Column 2: Atom"
				    write(*,*) "Column 3: Function type"
				    write(*,*) "Column 4: Exponent"
				    write(*,*) "Column 5: Coefficient"
				    do iEDFprim=1,nEDFprims
					    write(*,"(3i6,2f20.8)") iEDFprim,b_EDF(iEDFprim)%center,b_EDF(iEDFprim)%type,b_EDF(iEDFprim)%exp,CO_EDF(iEDFprim)
				    end do
			    end if
		    else if (i==100) then
			    call wfnsanity
            else if (i==201) then
                call ringring_geom
            else if (i==202) then
                call ring_rotate
            else if (i==1303) then
				call set_alpha_level()
		    end if
            
        !!!-----------------------------------------
	    !2000!!------------------- Very special functions
	    else if (isel==2000) then
			call energy_info_project
	    end if
        
    end if
    
end do !End main cycle

end program