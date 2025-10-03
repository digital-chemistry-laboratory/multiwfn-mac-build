!!! This module records electronic excitation information, shared by almost all electronic excitation analysis functions
!The excitfilename,ifiletypeexc,nstates,allexcene,allexcf,allexcmulti,allexcnorb are filled by "loadallexcinfo" routine
!The allexcdir,allorbleft,allorbright,allexccoeff are filled by "loadallexccoeff" routine
!All variables and arrays of user selected state is filled by "loadexccoeff"
!Notice that "loadallexcinfo" routine must be invoked at least one time prior to using loadallexccoeff and loadexccoeff
!
!When 50-50 is used in Gaussian, both singlet and triplet states are taken into account together. However, &
!when generation of triplet is specified in ORCA, only one set of excited state can be taken into account
module excitinfo
character :: excitfilename*200=" " !If nstates=0, then Multiwfn asks user to input
integer :: ifiletypeexc !1=Gaussian output; 2=ORCA output; 3=Plain text file; 4=Firefly output file
integer :: iORCAsTD=0 !1=ORCA with sTDA/sTDDFT, 0: Common ORCA TDA/TDDFT
integer :: iwarnORCA_TD=1 !If need to show warn when loading TDDFT coefficiens of ORCA, =1: Warn. Once warn once
integer :: nstates=0 !The total number of excited states (for ORCA, if triplet keyword is used, only one set of spin multiplicity states is loaded)
integer numexctime !The number of times of excited state printing of Gaussian/ORCA (1 corresponds to only once). >1 may due to geometry optimization or using state specific model
real*8,allocatable :: allexcene(:) !Excitation energies in eV
real*8,allocatable :: allexcf(:) !Oscillator strength
integer,allocatable :: allexcmulti(:) !Multiplicity of the states. 0 means the multiplicity is undefined (i.e. unrestricted reference state)
integer,allocatable :: allexcnorb(:) !The number of MO pairs in the states
!The last index in below arrays is the state index
integer,allocatable :: allexcdir(:,:) !1 means ->, 2 means <-
integer,allocatable :: allorbleft(:,:),allorbright(:,:) !Denote the MO at the left/right side in the excitation data (beta MO i is recorded as nbasis+i)
real*8,allocatable :: allexccoeff(:,:) !Coefficients of MO pairs
!Below is information of an interesting excited state selected by user
real*8 excene,excf
integer excmulti,excnorb
integer,allocatable :: excdir(:)
integer,allocatable :: orbleft(:),orbright(:)
real*8,allocatable :: exccoeff(:)
!Transition density matrix in basis function representation of alpha/total and beta electrons
real*8,allocatable :: tdmata(:,:),tdmatb(:,:)
end module


!---------- Deallocate excited state information
subroutine deallo_excitinfo
use excitinfo
nstates=0
excitfilename=" "
if (allocated(allexcene)) deallocate(allexcene)
if (allocated(allexcf)) deallocate(allexcf)
if (allocated(allexcmulti)) deallocate(allexcmulti)
if (allocated(allexcnorb)) deallocate(allexcnorb)
if (allocated(allexcdir)) deallocate(allexcdir)
if (allocated(allorbleft)) deallocate(allorbleft)
if (allocated(allorbright)) deallocate(allorbright)
if (allocated(allexccoeff)) deallocate(allexccoeff)
if (allocated(excdir)) deallocate(excdir)
if (allocated(orbleft)) deallocate(orbleft)
if (allocated(orbright)) deallocate(orbright)
if (allocated(exccoeff)) deallocate(exccoeff)
if (allocated(tdmata)) deallocate(tdmata)
if (allocated(tdmatb)) deallocate(tdmatb)
end subroutine


!-------- Main interface of various electron excitation analyses
subroutine excittrans_main
implicit real*8 (a-h,o-z)
do while(.true.)
	write(*,*)
	write(*,*) "           ============ Electronic excitation analyses ============ "
	write(*,*) "-1 Check, modify and export configuration coefficients of an excitation"
	write(*,*) "0 Return"
	write(*,"(a)") " 1 Analyze and visualize hole&electron distribution, transition density, and transition electric/magnetic dipole moment density"
	write(*,*) "2 Plot atom/fragment transition matrix of various kinds as heat map"
	write(*,*) "3 Analyze charge-transfer based on density difference grid data (JCTC,7,2498)"
	write(*,*) "4 Calculate delta_r index to measure charge-transfer length (JCTC,9,3118)"
	write(*,"(a)") " 5 Calculate transition electric/magnetic dipole moments between all states and for each state"
	write(*,*) "6 Generate natural transition orbitals (NTOs)"
	write(*,*) "7 Calculate ghost-hunter index (JCC,38,2151)"
	write(*,*) "8 Calculate interfragment charge transfer via IFCT method"
	write(*,*) "9 Generate and export transition density matrix"
	write(*,"(a)") " 10 Decompose transition electric/magnetic dipole moment as molecular orbital pair contributions"
	write(*,"(a)") " 11 Decompose transition electric/magnetic dipole moment as basis function and atom contributions"
	write(*,*) "12 Calculate Mulliken atomic transition charges"
	write(*,*) "13 Generate natural orbitals of specific excited states"
	write(*,*) "14 Calculate lambda index to characterize electron excitation (JCP,128,044118)"
	write(*,*) "15 Print major MO transitions in all excited states"
	write(*,*) "16 Charge-transfer spectrum (CTS) analysis (Carbon,187,78)"
    write(*,*) "17 Electron density polarization analysis based on electron excitations"

	read(*,*) isel
	if (isel==-1) then
		call modexccoeff
	else if (isel==0) then
		return
	else if (isel==1) then
		call hole_electron
	else if (isel==2) then
		call TDMplot
	else if (isel==3) then
		call CTanalyze
	else if (isel==4) then
		call delta_r
	else if (isel==5) then
		call exctransdip
	else if (isel==6) then
		call NTO
	else if (isel==7) then
		write(*,"(a)") " Note: To calculate the ghost-hunter index proposed in J. Comput. Chem., 38, 2151 (2017), you should use option 1 of subfunction 1 of main function 18 &
		&to calculate hole-electron distribution, then the index will be automatically printed. See Section 3.21.7 of the manual for more details."
! 		write(*,"(a)") " PS: The index calculated in this way is somewhat different to the original paper, &
! 		in which the 1/D term is calculated based on expensive relaxed density. If you really want to reproduce it, you can use subfunction 3 of main function 18&
! 		to calculate the 1/D term corresponding to relaxed density, and then manually calculate ghost-hunter index."
		write(*,*) "Press ENTER button to continue"
	else if (isel==8) then
		call IFCT
	else if (isel==9) then
		call exportTDM
	else if (isel==10) then
		call transdip_orbpair
	else if (isel==11) then
		call transdip_basatm
	else if (isel==12) then
		call transcharge
	else if (isel==13) then
		call genexcitNO
	else if (isel==14) then
		call lambda_excit
    else if (isel==15) then
        call majorMOtrans
    else if (isel==16) then
        call CTspectrum
    else if (isel==17) then
        call denspolar_excit
	end if
end do
end subroutine





!=============== Load and show basic information of all excited states
!ioutinfo=1: Output summary  =0: Do not output summary
!Loaded content: Total number of states, as well as spin multiplicity, energy, the number of MO pairs of each excited state
!If this information has already been loaded previously, then this routine only show summary
subroutine loadallexcinfo(ioutinfo)
use defvar
use excitinfo
use util
implicit real*8 (a-h,o-z)
integer ioutinfo
character c80tmp*80,c80tmp2*80,transmodestr*200,selectyn,c200tmp*200
!ifiletypeexc=1: Gaussian output file
!ifiletypeexc=2: ORCA output file
!ifiletypeexc=3: plain text file
!ifiletypeexc=4: Firefly TDDFT output file
!ifiletypeexc=5: GAMESS-US TDDFT output file
!ifiletypeexc=6: CP2K TDDFT output file
if (nstates>0) then
	write(*,"(' Note: Basic information of all excited states have been previously loaded from ',a)") trim(excitfilename)
else !The [excitfilename/=" ".and.nstates=0] case is involved in TDMplot
		
	if (excitfilename==" ") then
		if (ifiletype==10) then
			excitfilename=filename
		else
			write(*,"(a)") " Input path of Gaussian/ORCA/CP2K output file or plain text file, electron excitation information will be loaded from this file, see Section 3.21.A of manual for detail"
			write(*,*) "e.g. C:\lovelive\sunshine\yosoro.out"
            if (ifiletype==1.or.ifiletype==9.or.ifiletype==14) write(*,"(a)") " Hint: If pressing ENTER button directly, the file with same name as input file but &
            &with .out or .log suffix will be loaded"
			do while(.true.)
				read(*,"(a)") excitfilename
                if (excitfilename==" ") then
                    ipos=index(filename,'.',back=.true.)
                    excitfilename=trim(filename(:ipos))//"out"
				    inquire(file=excitfilename,exist=alive)
                    if (alive) then
                        write(*,"(a)") " Found "//trim(excitfilename)
                        exit
                    end if
                    excitfilename=trim(filename(:ipos))//"log"
				    inquire(file=excitfilename,exist=alive)
                    if (alive) then
                        write(*,"(a)") " Found "//trim(excitfilename)
                        exit
                    end if
                    write(*,"(' Error: Unable to find either ',a,' or ',a)") trim(filename(:ipos))//"out",trim(filename(:ipos))//"log"
                else
				    inquire(file=excitfilename,exist=alive)
				    if (alive) exit
				    write(*,*) "Cannot find this file, input again"
                end if
			end do
		end if
	end if

	open(10,file=excitfilename,status="old")
    
	call outputprog(10,iprog,1) !Determine type of output file and show information
    if (iprog==1) then !Gaussian output
		ifiletypeexc=1
    else if (iprog==2) then !ORCA output
		ifiletypeexc=2
    else if (iprog==3) then !GAMESS-US output
		ifiletypeexc=5
    else if (iprog==4) then !Firefly output
		ifiletypeexc=4
    else if (iprog==5) then !CP2K output
		ifiletypeexc=6
    else if(iprog==7) then !BDF output
        ifiletypeexc=7
    else if (iprog==0) then !Plain text file
		ifiletypeexc=3
    end if

	!Determine the number of excited states, so that proper size of arrays can be allocated
	nstates=0
	if (ifiletypeexc==1.or.ifiletypeexc==3) then !Gaussian output file or plain text file
		if (ifiletypeexc==1) then !Gaussian output file
			call loclabel(10,"Excitation energies and oscillator strengths:",ifound)
			if (ifound==0) then
				write(*,"(a)") " Error: This file is not output file of CIS/TDHF/TDDFT/TDA-DFT task, &
				&therefore cannot be used for present analysis. Please read Multiwfn manual Section 3.21 carefully"
				write(*,*) "Press ENTER button to return"
				read(*,*)
				return
			end if
            call loclabelfinal(10,"Excitation energies and oscillator strengths",numexctime)
            if (numexctime>1) then
			    write(*,"(a,i4,a)") " Note: Electron excitation information can be found",numexctime," times in the file, &
                &only the information printed last time will be loaded"
            end if
            read(10,*)
		else if (ifiletypeexc==3) then !Plain text file
			call loclabel(10,"Excited State",ifound)
			if (ifound==0) then
				write(*,"(a)") " Error: Unable to locate ""Excited State"" label"
				write(*,*) "Press ENTER button to return"
				read(*,*)
				return
			end if
		end if
		do while(.true.)
			call loclabel(10,"Excited State",ifound,0)
			if (ifound==1) then
				read(10,"(a)") c80tmp
                if ((index(c80tmp,'f=')/=0.and.ifiletypeexc==1).or.ifiletypeexc==3) nstates=nstates+1 !For Gaussian case, the line must contain f= to avoiding count irrelevant lines
			else
				exit
			end if
		end do
	else if (ifiletypeexc==2) then !ORCA output file
		call loclabel(10,"Number of roots to be determined",ifound)
		if (ifound==0) then
            call loclabel(10,"spectral range up to (eV)",ifound)
            if (ifound==1) then
                call loclabel(10,"roots found,",ifound,0)
                read(10,*) nstates
                iORCAsTD=1
            else
			    write(*,"(a)") "Error: This file is not output file of CIS/TDHF/TDDFT/TDA-DFT/SF-TDDFT/sTDA/sTDDFT task, &
			    &therefore cannot be used for present analysis. Please read Multiwfn manual Section 3.21 carefully"
			    write(*,*) "Press ENTER button to return"
			    read(*,*)
			    return
            end if
        else
            read(10,"(50x,i7)") nstates
		end if
        !This label is the most safest way of determining how many times of electronic excitation calculations have done
        !If singlet and triplet are both calculated, will be counted once
        call loclabelfinal(10,"TD-DFT XC SETUP",numexctime)
        if (numexctime>1) then
			write(*,"(a,i4,a)") " Note: Electron excitation information can be found",numexctime," times in the file, &
            &only the information printed last time will be loaded"
        end if
	else if (ifiletypeexc==4) then !Firefly output file
		call loclabel(10,"NUMBER OF STATES REQUESTED =",ifound)
		if (ifound==0) then
			write(*,*) "Error: It seems that this is not an electron excitation task!"
			write(*,*) "Press ENTER button to return"
			read(*,*)
			return
		end if
		read(10,"(29x,i6)") nstates
	else if (ifiletypeexc==5) then !GAMESS-US output file
		call loclabel(10,"NSTATE   =",ifound)
		if (ifound==0) then
			write(*,*) "Error: It seems that this is not an electron excitation task!"
			write(*,*) "Press ENTER button to return"
			read(*,*)
			return
		end if
		read(10,"(14x,i22)") nstates
	else if (ifiletypeexc==6) then !CP2K output file
		call loclabelfinal(10,"number   energy (eV)",ifound)
        if (ifound==0) then
			write(*,*) "Error: Unable to find electronic excitation information!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            return
        end if
        read(10,*)
        read(10,*)
		nstates=0
		do while(.true.)
			read(10,"(a)") c80tmp
            if (c80tmp==" ") exit
			nstates=nstates+1
		end do
    else if (ifiletypeexc==7) then !BDF output file. This part of code was contributed by Cong Wang, 2022-Dec-1
        call loclabel(10,"D<Pab>",ifound)
        if (ifound==0) then
            write(*,"(a)") " Error: Unable to locate ""D<Pab>"" label"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            return
        end if
        do while(.true.)
            call loclabel(10,"D<Pab>",ifound,0)
            if (ifound==1) then
                nstates=nstates+1
                read(10,*)
            else
                exit
            end if
        end do
	end if
    
	if (nstates>1) then
        if (maxloadexc==0) then
            write(*,"(' There are',i5,' excited states, loading basic information...')") nstates
        else if (maxloadexc<nstates) then
            write(*,"(' There are',i5,' excited states, however as ""maxloadexc"" in settings.ini has been set to',i5,', &
            &only the first these excited states are recognized')") nstates,maxloadexc
            nstates=maxloadexc
        end if
    end if
    write(*,*) 
	allocate(allexcene(nstates),allexcf(nstates),allexcmulti(nstates),allexcnorb(nstates))
	allexcnorb=0
	
	!Load excitation energy, multiplicity, oscillator strength, the number of MO pairs of each excited state
	if (ifiletypeexc==1) then !Gaussian output file
        rewind(10)
        do igeom=1,numexctime
		    call loclabel(10,"Excitation energies and oscillator strengths:",ifound,0)
            read(10,*)
        end do
		do iexc=1,nstates
			call loclabel(10,"Excited State",ifound,0)
			read(10,"(a)") transmodestr
			if (index(transmodestr,"Singlet")/=0) then
				allexcmulti(iexc)=1
			else if (index(transmodestr,"Triplet")/=0) then
				allexcmulti(iexc)=3
			else !May be unrestricted calculation, set as 0
				allexcmulti(iexc)=0
			end if
			do i=10,70
				if (transmodestr(i:i+1)=="eV") exit
			end do
			read(transmodestr(i-10:i-1),*) allexcene(iexc)
            itmp=index(transmodestr,'f=')
            read(transmodestr(itmp+2:),*) allexcf(iexc)
			!Count how many orbital pairs are involved in this excitation
			do while(.true.)
				read(10,"(a)") c80tmp
				if (index(c80tmp,'>')/=0) then
					allexcnorb(iexc)=allexcnorb(iexc)+1
				else if (index(c80tmp,'<')/=0) then
					allexcnorb(iexc)=allexcnorb(iexc)+1
				else
					exit
				end if
			end do
		end do
        
	else if (ifiletypeexc==2) then !ORCA output file
        if (iORCAsTD==0) then !Not sTD case
		    imultisel=1
		    if (wfntype==0.or.wfntype==3) then
			    call loclabel(10,"Generation of triplets") !When triplets=on, ORCA calculate both singlet and triplet excited state
			    read(10,"(a)") c80tmp
			    if (index(c80tmp," on ")/=0) then
				    write(*,*) "Load which kind of excited states?"
				    write(*,*) "1: Singlet   3: Triplet"
				    read(*,*) imultisel
				    allexcmulti=imultisel
                else
					allexcmulti=1 !Triplet excited states are not requested
			    end if
            else
				allexcmulti=0 !Spin multiplicity is undefined due to unrestricted reference
		    end if
            rewind(10)
            do igeom=1,numexctime
				call loclabel(10,"TD-DFT XC SETUP",ifound,0)
                read(10,*)
            end do
		    call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
		    if (imultisel==3) then
			    read(10,*)
			    call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
		    end if
		    do iexc=1,nstates
			    call loclabel(10,"STATE ",ifound,0)
			    read(10,"(a)") transmodestr
			    do i=10,70
				    if (transmodestr(i:i+1)=="eV") exit
			    end do
			    read(transmodestr(i-10:i-1),*) allexcene(iexc)
			    !Count how many orbital pairs are involved in this excitation
			    do while(.true.)
				    read(10,"(a)") c80tmp
				    if (index(c80tmp,'>')/=0) then
					    allexcnorb(iexc)=allexcnorb(iexc)+1
				    else if (index(c80tmp,'<')/=0) then
					    allexcnorb(iexc)=allexcnorb(iexc)+1
				    else
					    exit
				    end if
			    end do
		    end do
            if (imultisel==3) then !Triplet, spin-forbidden transition
                allexcf=0
            else !Singlet
                call loclabel(10,"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",ifound) !Load oscillator strengths
                call skiplines(10,5)
                read(10,"(a)") c80tmp
                if (index(c80tmp,'->')/=0) then
					iORCAver=6 !>=6.0.x
                    itmp=0
                else
					iORCAver=5 !<=5.0.x
                end if
                backspace(10)
                allexcf=0
				if (iORCAver==5) then !<=5.0.x
					do iexc=1,nstates
						read(10,*,iostat=ierror) itmp,rnouse,rnouse,tmpval
						if (ierror/=0) exit !For SF-TDDFT, number of states recorded in this field is less than nstates by 1, because one of SF-TDDFT states is viewed as ground state
						allexcf(itmp)=tmpval
					end do
                else if (iORCAver==6) then !>=6.0.x
					do iexc=1,2*nstates !Singlet and triplet states are outputted together, sorted by energy, so we need to examine 2*nstates
						read(10,*,iostat=ierror) c80tmp,c80tmp,c80tmp2,c80tmp,c80tmp,c80tmp,tmpval
						if (ierror/=0) exit !For SF-TDDFT, number of states recorded in this field is less than nstates by 1, because one of SF-TDDFT states is viewed as ground state
                        if (index(c80tmp2,"-1")/=0) then !e.g. 0-1A  ->  4-1A, where c80tmp2="4-1A"
							itmp=itmp+1
							allexcf(itmp)=tmpval
                        end if
					end do
                end if
            end if
        else if (iORCAsTD==1) then !sTDA or sTDDFT
            !If triplets=true, at least for ORCA 4.2, singlet and triplets are not outputted with clear labels, and transition strengths are always zero
            !In this case their information is loaded together
		    allexcmulti=1
            if (wfntype==0) then
                call loclabel(10,"calculate triplets",ifound)
                if (ifound==1) allexcmulti=0 !Spin multiplicity is undefined
            else
                allexcmulti=0 !Spin multiplicity is undefined due to unrestricted reference
            end if
            call loclabel(10,"state   eV")
            read(10,*)
		    do iexc=1,nstates
                read(10,"(a)") c200tmp
                read(c200tmp,*) inouse,allexcene(iexc),rnouse,allexcf(iexc)
                allexcnorb(iexc)=strcharnum(c200tmp,'(')
            end do
        end if
        
	else if (ifiletypeexc==4) then !Firefly output file
		call loclabel(10,"EXCITED STATE   1 ",ifound)
		do iexc=1,nstates
			read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)
			do while(.true.)
				read(10,"(a)") c80tmp
				if (c80tmp(11:11)=="-") then
					exit
				else
					allexcnorb(iexc)=allexcnorb(iexc)+1
				end if
			end do
			read(10,*);read(10,*)
		end do
		call loclabel(10,"STATE       HARTREE        EV",ifound,0)
		read(10,*);read(10,*)
		do iexc=1,nstates
			read(10,"(a)") c80tmp
			read(c80tmp(1:3),*) allexcmulti(iexc)
			read(c80tmp(26:33),*) allexcene(iexc)
			read(c80tmp(71:),*) allexcf(iexc)
		end do
        
	else if (ifiletypeexc==5) then !GAMESS-US output file
		if (wfntype==0) then
			allexcmulti=1
			call loclabel(10,"TRIPLET EXCITATIONS",ifound)
			if (ifound==1) allexcmulti=3
		else
			allexcmulti=0 !Undefined excited state multiplicity due to unrestricted reference
		end if
		call loclabel(10,"STATE #   1",ifound)
		do iexc=1,nstates
			read(10,"(22x,f12.6)") allexcene(iexc)
			read(10,"(22x,f12.6)") allexcf(iexc)
			read(10,*);read(10,*);read(10,*);read(10,*)
			if (wfntype==0) then
				read(10,*);read(10,*)
			end if
			do while(.true.)
				read(10,"(a)") c80tmp
				if (c80tmp==" ") then
					exit
				else
					allexcnorb(iexc)=allexcnorb(iexc)+2 !For TDDFT, excitation and de-excitation are outputted as a single line
				end if
			end do
		end do
        
	else if (ifiletypeexc==6) then !CP2K output file
        call loclabelfinal(10,"TDDFPT states of multiplicity",ifound)
        read(10,"(a)") c80tmp
        if (index(c80tmp,"U-TDDFPT")/=0) then
			allexcmulti=0 !Excited states are not spin pure states
        else
			allexcmulti=1
			if (index(c80tmp,'3')/=0) then
				imultisel=3
				call loclabel(10,"TDDFPT states of multiplicity 1",ifound) !SOC-TDDFT calculate both singlet and triplet excited states
				if (ifound==1) then
					write(*,*) "Load which kind of excited states?"
					write(*,*) "1: Singlet   3: Triplet"
					read(*,*) imultisel
				end if
                allexcmulti=imultisel
				if (imultisel==1) then
					call loclabelfinal(10,"TDDFPT states of multiplicity 1",ifound)
                else
					call loclabelfinal(10,"TDDFPT states of multiplicity 3",ifound)
                end if
			end if
        end if
        call loclabel(10,"TDDFPT|",ifound,0) !Load oscillator strength
		do iexc=1,nstates
			read(10,*) c80tmp,itmp,allexcene(iexc),r1,r2,r3,allexcf(iexc)
        end do
        call loclabel(10,"number             orbital",ifound,0)
        read(10,*);read(10,*)
		do iexc=1,nstates
			read(10,*)
			do while(.true.)
				read(10,"(a)") c80tmp
                if (c80tmp(2:6)=="-----") exit
				if (index(c80tmp,"eV")/=0) then
					backspace(10)
					exit
                end if
                allexcnorb(iexc)=allexcnorb(iexc)+1
			end do
		end do

    else if (ifiletypeexc==7) then !BDF output file. This part of code was contributed by Cong Wang, 2022-Dec-1
        rewind(10)
        isf=111
        iexc=0
        do while(.true.)
            read(10,"(a)",iostat=ierror) c200tmp
            if (ierror/=0) exit
            if (index(c200tmp,"Start TD-DFT Calculations for isf")/=0) then
                read(c200tmp(45:45),*) isf
            else if (index(c200tmp,"D<Pab>")/=0) then
                iexc=iexc+1
                if(iexc>nstates) then
                    write(*,*) "iexc>nstates"
                    stop
                end if
                if (isf==0) then
                    allexcmulti(iexc)=1
                else if (isf==1) then
                    allexcmulti(iexc)=3
                else !May be unrestricted calculation, set as 0
                    allexcmulti(iexc)=0
                end if
                do i=1,50
                    if (c200tmp(i:i+1)=="eV") exit
                end do
                read(c200tmp(i-10:i-1),*) allexcene(iexc)
                itmp=index(c200tmp,'f=')
                read(c200tmp(itmp+6:),*) allexcf(iexc)
                !Count how many orbital pairs are involved in this excitation
                do while(.true.)
                    read(10,"(a)") c200tmp
                    if (index(c200tmp,'>')/=0) then
                        allexcnorb(iexc)=allexcnorb(iexc)+1
                    else if (index(c200tmp,'<')/=0) then
                        allexcnorb(iexc)=allexcnorb(iexc)+1
                    else
                        exit
                    end if
                end do
            end if
        end do    
	else if (ifiletypeexc==3) then !Plain text file
		rewind(10)
		do iexc=1,nstates
			call loclabel(10,"Excited State",ifound,0)
			read(10,*) c80tmp,c80tmp,inouse,allexcmulti(iexc),allexcene(iexc)
			!Count how many orbital pairs are involved in this excitation
			do while(.true.)
				read(10,"(a)",iostat=ierror) c80tmp
				if (c80tmp==" ".or.ierror/=0) exit
				if (index(c80tmp,'>')/=0) then
					allexcnorb(iexc)=allexcnorb(iexc)+1
				else if (index(c80tmp,'<')/=0) then
					allexcnorb(iexc)=allexcnorb(iexc)+1
				end if
			end do
		end do
	end if
	
	close(10)
end if

if (ioutinfo==1) then
	if (nstates>1) then
		write(*,*) "Summary of excited states:"
		do iexc=1,nstates
			if (allexcmulti(iexc)>0) then
				write(*,"(' State:',i5,'    Exc. Energy:',f8.3,' eV   Multi.:',i2,'    MO pairs:',i8)") iexc,allexcene(iexc),allexcmulti(iexc),allexcnorb(iexc)
			else
				write(*,"(' State:',i5,'    Exc. Energy:',f8.3,' eV   Multi.: ?    MO pairs:',i8)") iexc,allexcene(iexc),allexcnorb(iexc)
			end if
		end do
	else
		write(*,*) "Information of the only excited state:"
		if (allexcmulti(1)>0) then
			write(*,"(' Excitation Energy:',f8.3,' eV   Multiplicity:',i2,'    MO pairs:',i8)") allexcene(1),allexcmulti(1),allexcnorb(1)
		else
			write(*,"(' Excitation Energy:',f8.3,' eV   Multiplicity: ?    MO pairs:',i8)") allexcene(1),allexcnorb(1)
		end if
	end if
	write(*,*)
end if

end subroutine



!================ Load details of all excited states
!If information has already been loaded previously, they will not be loaded again
!Information to be loaded: MO indices, excitation directions, determinant coefficients
!ioutinfo=1: Output summary  =0: Do not output summary  =2: Print nothing
subroutine loadallexccoeff(ioutinfo)
use defvar
use excitinfo
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,leftstr*80,rightstr*80

if (allocated(allexcdir)) then
	write(*,"(a)") " Detailed information of all excited states have already been loaded previously, now directly employ them"
else
	maxpair=maxval(allexcnorb)
	allocate(allexcdir(maxpair,nstates),allorbleft(maxpair,nstates),allorbright(maxpair,nstates),allexccoeff(maxpair,nstates))

	open(10,file=excitfilename,status="old")
	write(*,*) "Loading configuration coefficients..."

	!Notice that for unrestricted case, A and B are combined as single index, namely if orbital index is larger than nbasis, then it is B, else A
	if (ifiletypeexc==1.or.ifiletypeexc==3) then !Gaussian output file or plain text file
		if (ifiletypeexc==1) then !Gaussian output file
            do igeom=1,numexctime
		        call loclabel(10,"Excitation energies and oscillator strengths:",ifound,0)
                read(10,*)
            end do
        end if
		do iexc=1,nstates
			call loclabel(10,"Excited State",ifound,0)
			read(10,*)
			do itmp=1,allexcnorb(iexc)
				!Determine excitation direction
				read(10,"(a)") c80tmp
				ileft=index(c80tmp,'->')
				allexcdir(itmp,iexc)=1 !means ->
				isign=ileft
				if (ileft==0) then
					iright=index(c80tmp,'<-')
					allexcdir(itmp,iexc)=2 !means <-
					isign=iright
				end if
				!Process left side
				leftstr=c80tmp(:isign-1)
				if (wfntype==0.or.wfntype==3) then
					read(leftstr,*) allorbleft(itmp,iexc)
				else
					if (index(leftstr,'A')/=0) then !Alpha
						read(leftstr(:len_trim(leftstr)-1),*) allorbleft(itmp,iexc)
					else !Beta
						read(leftstr(:len_trim(leftstr)-1),*) allorbleft(itmp,iexc)
						allorbleft(itmp,iexc)=allorbleft(itmp,iexc)+nbasis
					end if
				end if
				!Process right side
				rightstr=c80tmp(isign+2:)
				if (wfntype==0.or.wfntype==3) then
					read(rightstr,*) allorbright(itmp,iexc),allexccoeff(itmp,iexc)
				else
					do isplit=1,80
						if (rightstr(isplit:isplit)=='A'.or.rightstr(isplit:isplit)=='B') exit
					end do
					read(rightstr(:isplit-1),*) allorbright(itmp,iexc)
					read(rightstr(isplit+1:),*) allexccoeff(itmp,iexc)
					if (rightstr(isplit:isplit)=='B') allorbright(itmp,iexc)=allorbright(itmp,iexc)+nbasis
				end if
			end do
		end do
		
	else if (ifiletypeexc==2) then !ORCA output file
        if (iORCAsTD==0) then !Regular case
		    !Worthnotingly, in at least ORCA 4.0, de-excitation is not separately outputted as <-, but combined into ->
		    !Here we still check <-, because hopefully Neese may change the convention of ORCA output in the future...
            do igeom=1,numexctime
                call loclabel(10,"TD-DFT XC SETUP",ifound,0)
                read(10,*)
            end do
		    call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
		    if (allexcmulti(1)==3) then !When triplets=on, ORCA calculates both singlet and triplet excited states, now move to the latter
			    read(10,*)
			    call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
		    end if
		    do iexc=1,nstates
			    call loclabel(10,"STATE ",ifound,0)
			    read(10,*)
			    do itmp=1,allexcnorb(iexc)
				    read(10,"(a)") c80tmp
				    if (index(c80tmp,'->')/=0) then
					    allexcdir(itmp,iexc)=1
				    else
					    allexcdir(itmp,iexc)=2
				    end if
				    do isign=1,80 !Find position of <- or ->
					    if (c80tmp(isign:isign)=='-'.or.c80tmp(isign:isign)=='<') exit
				    end do
				    !Process left side of <- or ->
				    read(c80tmp(:isign-1),"(a)") leftstr
				    read(leftstr(:len_trim(leftstr)-1),*) allorbleft(itmp,iexc)
				    allorbleft(itmp,iexc)=allorbleft(itmp,iexc)+1 !ORCA counts orbital from 0 rather than 1!!!
				    if (index(leftstr,'b')/=0) allorbleft(itmp,iexc)=allorbleft(itmp,iexc)+nbasis
				    !Process right side of <- or ->
				    read(c80tmp(isign+2:),*) rightstr
				    read(rightstr(:len_trim(rightstr)-1),*) allorbright(itmp,iexc)
				    allorbright(itmp,iexc)=allorbright(itmp,iexc)+1
				    if (index(rightstr,'b')/=0) allorbright(itmp,iexc)=allorbright(itmp,iexc)+nbasis
				    iTDA=index(c80tmp,'c=')
				    if (iTDA/=0) then !CIS, TDA task, both configuration contribution and coefficients are presented, e.g. 2a ->   5a  :     0.985689 (c=  0.99281847)
					    read(c80tmp(iTDA+2:iTDA+13),*) allexccoeff(itmp,iexc)
				    else !TD task. Positive contribution of i->a and negative contribution a<-i are summed up and printed, e.g. 2a ->   6a  :     0.968777
					    if (iexc==1.and.itmp==1) then
						    write(*,"(a)") " Warning: For TD task, ORCA does not print configuration coefficients but only print contributions of each MO pair to excitation, &
						    &in this case Multiwfn guesses configuration coefficients by calculating square root of the contributions. &
                            &However, this treatment may lead to fully misleading result, you should consider using TDA instead, which is perfectly supported."
						    write(*,*) "If you really want to proceed, press ENTER button"
						    read(*,*)
					    end if
					    read(c80tmp(23:32),*) tmpval
					    if (tmpval<0) allexcdir(itmp,iexc)=2 !Negative contribution is assumed to be significant de-excitation (of course this is not strict since -> and <- have been combined together)
					    allexccoeff(itmp,iexc)=dsqrt(abs(tmpval))
				    end if
				    !Although for closed-shell reference state, ORCA still outputs coefficients as normalization to 1.0, &
				    !However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
				    if (wfntype==0.or.wfntype==3) allexccoeff(itmp,iexc)=allexccoeff(itmp,iexc)/dsqrt(2D0)
			    end do
		    end do
        else if (iORCAsTD==1) then
            !sTDA/sTDDFT ignores very low lying occupied MOs and very high lying virtual MOs, &
            !we need to know how many core orbitals are ignored, so that actual MO index can be obtained
            call loclabel(10,"occ. MOs in")
            read(10,"(a)") c80tmp
            itmp=index(c80tmp,'...')
            if (wfntype==0) then
                read(c80tmp(itmp+3:),*) ntmp
                ncor=nint(naelec)-ntmp
            else
                read(c80tmp(itmp+3:),*) natmp,nbtmp
                nacor=nint(naelec)-natmp
                nbcor=nint(nbelec)-nbtmp
            end if
            call loclabel(10,"state   eV",ifound,0)
            read(10,*)
		    do iexc=1,nstates
                read(10,"(a)") c200tmp
                itmp=index(c200tmp,'(')
                ipos=itmp-9
                do ipair=1,allexcnorb(iexc)
                    if (wfntype==0) then
                        c80tmp=c200tmp(ipos:ipos+20)
                    else
                        c80tmp=c200tmp(ipos:ipos+22)
                    end if
                    read(c80tmp,*) allexccoeff(ipair,iexc)
                    if (wfntype==0) allexccoeff(ipair,iexc)=allexccoeff(ipair,iexc)/dsqrt(2D0)
                    if (index(c80tmp,'>')/=0) then
                        allexcdir(ipair,iexc)=1
                    else
                        allexcdir(ipair,iexc)=2
                    end if
                    itmp=index(c80tmp,'(') !Load left orbital index
                    if (wfntype==0) then
                        read(c80tmp(itmp+1:itmp+4),*) allorbleft(ipair,iexc)
                        allorbleft(ipair,iexc)=allorbleft(ipair,iexc)+ncor+1  !Note that the MO is 0 based index in sTDA/sTDDFT
                    else
                        read(c80tmp(itmp+1:itmp+4),*) allorbleft(ipair,iexc)
                        if (c80tmp(itmp+5:itmp+5)=='a') then
                            allorbleft(ipair,iexc)=allorbleft(ipair,iexc)+nacor+1
                        else
                            allorbleft(ipair,iexc)=allorbleft(ipair,iexc)+nbcor+1+nbasis
                        end if
                    end if
                    itmp=index(c80tmp,')') !Load right orbital index
                    if (wfntype==0) then
                        read(c80tmp(itmp-4:itmp-1),*) allorbright(ipair,iexc)
                        allorbright(ipair,iexc)=allorbright(ipair,iexc)+ncor+1
                    else
                        read(c80tmp(itmp-5:itmp-2),*) allorbright(ipair,iexc)
                        if (c80tmp(itmp-1:itmp-1)=='a') then
                            allorbright(ipair,iexc)=allorbright(ipair,iexc)+nacor+1
                        else
                            allorbright(ipair,iexc)=allorbright(ipair,iexc)+nbcor+1+nbasis
                        end if
                    end if
                    if (wfntype==0) then
                        ipos=ipos+21
                    else
                        ipos=ipos+23
                    end if
                end do
            end do
		end if
        
	else if (ifiletypeexc==4) then !Firefly output file
		call loclabel(10,"EXCITED STATE   1 ",ifound)
		do iexc=1,nstates
			read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)
			do itmp=1,allexcnorb(iexc)
				read(10,*) imo,jmo,allexccoeff(itmp,iexc)
				if (jmo>imo) then !Excitation
					allorbleft(itmp,iexc)=imo
					allorbright(itmp,iexc)=jmo
					allexcdir(itmp,iexc)=1
				else !De-excitation
					allorbleft(itmp,iexc)=jmo
					allorbright(itmp,iexc)=imo
					allexcdir(itmp,iexc)=2
				end if
			end do
			read(10,*);read(10,*);read(10,*)
		end do
		!Although for closed-shell reference state, Firefly still outputs coefficients as normalization to 1.0, &
		!However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
		if (wfntype==0.or.wfntype==3) allexccoeff=allexccoeff/dsqrt(2D0)
		
	else if (ifiletypeexc==5) then !GAMESS-US output file
		call loclabel(10," STATE #   1",ifound)
		do iexc=1,nstates
			if (wfntype==0) then !Restricted reference state
				read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)
				itmp=0
				do ipair=1,allexcnorb(iexc)/2
					read(10,*) imo,jmo,coeff1,coeff2
					itmp=itmp+1
					allorbleft(itmp,iexc)=imo
					allorbright(itmp,iexc)=jmo
					allexcdir(itmp,iexc)=1
					allexccoeff(itmp,iexc)=coeff1
					itmp=itmp+1
					allorbleft(itmp,iexc)=imo
					allorbright(itmp,iexc)=jmo
					allexcdir(itmp,iexc)=2
					allexccoeff(itmp,iexc)=coeff2
				end do
			else !Unrestricted reference state
				read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)
				itmp=0
				do ipair=1,allexcnorb(iexc)/2
					read(10,"(a)") c80tmp
					read(c80tmp,*) imo,jmo
					read(c80tmp(22:48),*) coeff1,coeff2
					ibeta=0
					if (index(c80tmp,"BET")/=0) ibeta=1
					itmp=itmp+1
					allorbleft(itmp,iexc)=imo
					allorbright(itmp,iexc)=jmo
					allexcdir(itmp,iexc)=1
					allexccoeff(itmp,iexc)=coeff1
					itmp=itmp+1
					allorbleft(itmp,iexc)=imo
					allorbright(itmp,iexc)=jmo
					allexcdir(itmp,iexc)=2
					allexccoeff(itmp,iexc)=coeff2
					if (ibeta==1) then
						allorbleft(itmp-1:itmp,iexc)=allorbleft(itmp-1:itmp,iexc)+nbasis
						allorbright(itmp-1:itmp,iexc)=allorbright(itmp-1:itmp,iexc)+nbasis
					end if
				end do
			end if
			read(10,*)
		end do
		!Even for closed-shell reference state, GAMESS-US still outputs coefficients as normalization to 1.0, &
		!However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
		if (wfntype==0.or.wfntype==3) allexccoeff=allexccoeff/dsqrt(2D0)
        
	else if (ifiletypeexc==6) then !CP2K output file
		allexcdir(:,:)=1 !CP2K is fully based on TDA
        if (allexcmulti(1)==1) then
			call loclabelfinal(10,"TDDFPT states of multiplicity 1",ifound)
			call loclabel(10,"number             orbital",ifound,0)
        else if (allexcmulti(1)==3) then
			call loclabelfinal(10,"TDDFPT states of multiplicity 3",ifound)
			call loclabel(10,"number             orbital",ifound,0)
		else
			call loclabelfinal(10,"number             orbital",ifound)
        end if
		read(10,*);read(10,*)
		do iexc=1,nstates
			read(10,*)
			do ipair=1,allexcnorb(iexc)
				read(10,"(a)") c80tmp
				read(c80tmp(27:34),*) allorbleft(ipair,iexc)
				read(c80tmp(49:55),*) allorbright(ipair,iexc)
				read(c80tmp(69:78),*) allexccoeff(ipair,iexc)
                if (c80tmp(37:39)=="bet") allorbleft(ipair,iexc)=allorbleft(ipair,iexc)+nbasis
                if (c80tmp(58:60)=="bet") allorbright(ipair,iexc)=allorbright(ipair,iexc)+nbasis
			end do
		end do
		!Even for closed-shell reference state, CP2K still outputs coefficients as normalization to 1.0, &
		!However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
		if (wfntype==0.or.wfntype==3) allexccoeff=allexccoeff/dsqrt(2D0)

    else if (ifiletypeexc==7) then !BDF output file. This part of code was contributed by Cong Wang, 2022-Dec-1
        do iexc=1,nstates
            call loclabel(10,"D<Pab>",ifound,0)
            read(10,*)
            do itmp=1,allexcnorb(iexc)
                !Determine excitation direction
                read(10,"(a)") c200tmp
                ileft=index(c200tmp,"->")
                allexcdir(itmp,iexc)=1 !means ->
                isign=ileft
                if (ileft==0) then
                    iright=index(c200tmp,"<-")
                    allexcdir(itmp,iexc)=2 !means <-
                    isign=iright
                end if
                !Process left side
                leftstr=c200tmp(:isign-1)
                do ileft=len_trim(leftstr),1,-1
                    if (leftstr(ileft:ileft)=="(") exit
                end do
                read(leftstr(ileft:len_trim(leftstr)-1),*) allorbleft(itmp,iexc)
                !Process right side
                rightstr=c200tmp(isign+2:)
                ir1=0
                ir2=0
                do iright=1,len_trim(rightstr)
                    if (rightstr(iright:iright)=="(") ir1=iright
                    if (rightstr(iright:iright)==")") then
                        ir2=iright
                        exit
                    end if
                end do
                read(rightstr(ir1+1:ir2-1),*) allorbright(itmp,iexc)
                ir1=0
                ir2=0
                do iright=1,len_trim(rightstr)
                    if (rightstr(iright:iright+3)=="Per:") ir1=iright
                    if (rightstr(iright:iright)=="%") then
                        ir2=iright
                        exit
                    end if
                end do
                read(rightstr(ir1+4:ir2-1),*) fper
                allexccoeff(itmp,iexc) = sqrt(fper/200D0)
            end do
        end do
	end if
	
	close(10)
end if

if (ioutinfo==1) then
	write(*,*) "Summary of excited states:"
	write(*,*) "Exc.state#     Exc.energy(eV)     Multi.   MO pairs    Normalization"
	do iexc=1,nstates
		sumsqrexc=0
		sumsqrdeexc=0
		do itmp=1,allexcnorb(iexc)
			if (allexcdir(itmp,iexc)==1) sumsqrexc=sumsqrexc+allexccoeff(itmp,iexc)**2
			if (allexcdir(itmp,iexc)==2) sumsqrdeexc=sumsqrdeexc-allexccoeff(itmp,iexc)**2
		end do
		sumsqrall=sumsqrexc+sumsqrdeexc
		if (allexcmulti(iexc)/=0) then
			write(*,"(i8,f18.5,4x,i8,i13,f16.6)") iexc,allexcene(iexc),allexcmulti(iexc),allexcnorb(iexc),sumsqrall
		else
			write(*,"(i8,f18.5,4x,'    ?   ',i13,f16.6)") iexc,allexcene(iexc),allexcnorb(iexc),sumsqrall
		end if
	end do
	write(*,*)
end if
end subroutine



!================ Load details of a selected excited states
!Loaded content: MO index, excitation direction, determinant coefficient
!ioutinfo=1: Print intermediate information and summary  =0: Do not print
subroutine loadexccoeff(istate,ioutinfo)
use defvar
use excitinfo
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,leftstr*80,rightstr*80
integer istate,ioutinfo
real*8,allocatable :: coefarr1(:),coefarr2(:)

excene=allexcene(istate)
excf=allexcf(istate)
excmulti=allexcmulti(istate)
excnorb=allexcnorb(istate)
if (allocated(excdir)) deallocate(excdir,orbleft,orbright,exccoeff)
allocate(excdir(excnorb),orbleft(excnorb),orbright(excnorb),exccoeff(excnorb))

!Directly take data from previously loaded detailed information of all states. Rarely used feature, so comment it
! if (allocated(allexcmulti)) then
! orbleft(:)=allorbleft(:,istate)
! orbright(:)=allorbright(:,istate)
! excdir(:)=allexcdir(:,istate)
! exccoeff(:)=allexccoeff(:,istate)
! return
! end if

open(10,file=excitfilename,status="old")
if (ioutinfo==1) write(*,"(' Loading configuration coefficients of excited state',i4,'...')") istate

!Notice that for unrestricted case, A and B are combined as single index, namely if orbital index is larger than nbasis, then it is B, else A
if (ifiletypeexc==1.or.ifiletypeexc==3) then !Gaussian output file or plain text file
	if (ifiletypeexc==1) then
        do igeom=1,numexctime
		    call loclabel(10,"Excitation energies and oscillator strengths:",ifound,0)
            read(10,*)
        end do
	else
		rewind(10)
	end if
	do iexc=1,istate
		call loclabel(10,"Excited State",ifound,0)
		read(10,*)
	end do
	do itmp=1,excnorb
		!Determine excitation direction
		read(10,"(a)") c80tmp
		ileft=index(c80tmp,'->')
		excdir(itmp)=1 !means ->
		isign=ileft
		if (ileft==0) then
			iright=index(c80tmp,'<-')
			excdir(itmp)=2 !means <-
			isign=iright
		end if
		!Process left side
		leftstr=c80tmp(:isign-1)
		if (wfntype==0.or.wfntype==3) then
			read(leftstr,*) orbleft(itmp)
		else
			if (index(leftstr,'A')/=0) then !Alpha
				read(leftstr(:len_trim(leftstr)-1),*) orbleft(itmp)
			else !Beta
				read(leftstr(:len_trim(leftstr)-1),*) orbleft(itmp)
				orbleft(itmp)=orbleft(itmp)+nbasis
			end if
		end if
		!Process right side
		rightstr=c80tmp(isign+2:)
		if (wfntype==0.or.wfntype==3) then
			read(rightstr,*) orbright(itmp),exccoeff(itmp)
		else
			do isplit=1,80
				if (rightstr(isplit:isplit)=='A'.or.rightstr(isplit:isplit)=='B') exit
			end do
			read(rightstr(:isplit-1),*) orbright(itmp)
			read(rightstr(isplit+1:),*) exccoeff(itmp)
			if (rightstr(isplit:isplit)=='B') orbright(itmp)=orbright(itmp)+nbasis
		end if
	end do
	
else if (ifiletypeexc==2) then !ORCA output file
    if (iORCAsTD==0) then !Regular case
		iORCAjson=0
		do itmp=len(excitfilename),1,-1
			if (excitfilename(itmp:itmp)=='.') exit
        end do
        c200tmp=trim(excitfilename(1:itmp))//"json"
		inquire(file=c200tmp,exist=alive)
		if (alive.eqv..true.) then
			write(*,"(1x,a,a)") trim(c200tmp)," was found in current folder, do you want to load configuration coefficients from this file instead of from ORCA output file? (y/n)"
            read(*,*) c80tmp
            if (c80tmp=='y'.or.c80tmp=='Y') iORCAjson=1
        end if
		if (iORCAjson==1) then !Load from json file
			close(10)
            write(*,*) "Loading "//trim(c200tmp)//"..."
            open(10,file=c200tmp,status="old")
			if (istate>1) then
				do iexc=1,istate-1
					call loclabel(10,"IRoot",ifound,0)
					if (ifound==0) then
						write(*,"(a)") " Error: Unable to find ""IRoot"" from the json file! The json file was not prepared in proper way, &
                        &please check Section 3.21.A of Multiwfn manual carefully"
						write(*,*) "Press ENTER button to exit program"
						read(*,*)
                        stop
					end if
					read(10,*)
				end do
            end if
            critkeep=1E-14
            !We first load all coefficients
            call loclabel(10,"IRoot",ifound,0)
			call loclabel(10,"OrbWin",ifound,0)
            read(10,*)
            read(10,*) iocclow
            read(10,*) ioccup
            read(10,*) ivirlow
            read(10,*) ivirup
            iocclow=iocclow+1 !ORCA counts from 0
            ioccup=ioccup+1
            ivirlow=ivirlow+1
            ivirup=ivirup+1
            ntotcoef=(ioccup-iocclow+1)*(ivirup-ivirlow+1)
            allocate(coefarr1(ntotcoef))
			call loclabel(10,"X",ifound,0)
            read(10,*)
            itmp=0
            read(10,*) coefarr1(:)
			call loclabel(10,"X-Y",ifTD,0)
            if (ifTD==0) then
				write(*,"(i8,' configuration coefficients were loaded')") ntotcoef
                excnorb=count(abs(coefarr1)>critkeep)
            else if (ifTD==1) then
				write(*,"(i8,' X+Y coefficients were loaded')") ntotcoef
				allocate(coefarr2(ntotcoef))
				read(10,*)
				read(10,*) coefarr2(:)
				write(*,"(i8,' X-Y coefficients were loaded')") ntotcoef
                !facnorm=dsqrt(sum(coefarr2(:)*coefarr1(:))) !The coefficients in json file were not normalized, we renormalize
                !write(*,"(' Normalization factor:',f16.10)") facnorm
                !coefarr1 and coefarr2 records X+Y and X-Y coefficients, respectively. Now convert to X and Y coefficients
                do itmp=1,ntotcoef
					XpY=coefarr1(itmp)
					XnY=coefarr2(itmp)
                    coefarr1(itmp)=(XpY+XnY)/2D0
                    coefarr2(itmp)=(XpY-XnY)/2D0
                end do
                excnorb=count(abs(coefarr1)>critkeep) + count(abs(coefarr2)>critkeep)
			end if
            write(*,"(' Absolute configuration coefficients larger than',1PE16.8,' will be kept')") critkeep
			deallocate(excdir,orbleft,orbright,exccoeff)
			allocate(excdir(excnorb),orbleft(excnorb),orbright(excnorb),exccoeff(excnorb))
            !Only retain nonnegligible configuration coefficients to final arrays
			excnorb=0
            itmp=0
			do iocc=iocclow,ioccup
				do ivir=ivirlow,ivirup
					itmp=itmp+1
                    !write(*,*) iocc,ivir,coefarr1(itmp)
					if (abs(coefarr1(itmp))>critkeep) then
						excnorb=excnorb+1
						exccoeff(excnorb)=coefarr1(itmp)
						orbleft(excnorb)=iocc
						orbright(excnorb)=ivir
						excdir(excnorb)=1
					end if
				end do
			end do
            write(*,"(' Number of kept excitation configuration coefficients:',i8)") excnorb
            ntmp=excnorb
            if (ifTD==1) then
				itmp=0
				do iocc=iocclow,ioccup
					do ivir=ivirlow,ivirup
						itmp=itmp+1
                        !write(*,*) iocc,ivir,coefarr2(itmp)
						if (abs(coefarr2(itmp))>critkeep) then
							excnorb=excnorb+1
							exccoeff(excnorb)=coefarr2(itmp)
							orbleft(excnorb)=iocc
							orbright(excnorb)=ivir
							excdir(excnorb)=2
						end if
					end do
				end do
				write(*,"(' Number of kept de-excitation configuration coefficients:',i8)") excnorb-ntmp
            end if
            
            facnorm=dsqrt(sum(coefarr1(:)**2)-sum(coefarr2(:)**2)) !The coefficients in json file were not normalized, we renormalize
            write(*,"(' Normalization factor:',f16.10)") facnorm
            exccoeff(:)=exccoeff(:)/facnorm
            
            if (wfntype==0.or.wfntype==3) exccoeff(:)=exccoeff(:)/dsqrt(2D0)
		else !Load from output file
			!Worthnotingly, in at least ORCA 4.0, de-excitation is not separately outputted as <-, but combined into ->
			!Here we still check <-, because hopefully Neese may change the convention of ORCA output in the future...
			do igeom=1,numexctime
				call loclabel(10,"TD-DFT XC SETUP",ifound,0)
				read(10,*)
			end do
			call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
			if (allexcmulti(1)==3) then !When triplets=on, ORCA calculates both singlet and triplet excited states, now move to the latter
				read(10,*)
				call loclabel(10,"the weight of the individual excitations are printed",ifound,0)
			end if
			do iexc=1,istate
				call loclabel(10,"STATE ",ifound,0)
				read(10,*)
				if (iexc==istate) then
					do itmp=1,excnorb
						read(10,"(a)") c80tmp
						if (index(c80tmp,'->')/=0) then
							excdir(itmp)=1
						else
							excdir(itmp)=2
						end if
						do isign=1,80 !Find position of <- or ->
							if (c80tmp(isign:isign)=='-'.or.c80tmp(isign:isign)=='<') exit
						end do
						!Process left side of <- or ->
						read(c80tmp(:isign-1),"(a)") leftstr
						read(leftstr(:len_trim(leftstr)-1),*) orbleft(itmp)
						orbleft(itmp)=orbleft(itmp)+1 !ORCA counts orbital from 0 rather than 1!!!
						if (index(leftstr,'b')/=0) orbleft(itmp)=orbleft(itmp)+nbasis
						!Process right side of <- or ->
						read(c80tmp(isign+2:),*) rightstr
						read(rightstr(:len_trim(rightstr)-1),*) orbright(itmp)
						orbright(itmp)=orbright(itmp)+1
						if (index(rightstr,'b')/=0) orbright(itmp)=orbright(itmp)+nbasis
						iTDA=index(c80tmp,'c=')
						if (iTDA/=0) then !CIS, TDA task, configuration coefficients are presented
							read(c80tmp(iTDA+2:iTDA+13),*) exccoeff(itmp)
						else !TD task, configuration coefficients are not presented. Contribution of i->a and i<-a are summed up and outputted as i->a
							if (iwarnORCA_TD==1.and.itmp==1) then
								write(*,"(a)") " Warning: For TD task, ORCA does not print configuration coefficients but only print corresponding contributions of each orbital pair, &
								&in this case Multiwfn determines configuration coefficients simply as square root of contribution values. However, this treatment is &
								&evidently inappropriate and the result is nonsense when de-excitation is significant (In this situation you have to use TDA-DFT instead)"
								write(*,*) "If you really want to proceed, press ENTER button to continue"
								read(*,*)
								iwarnORCA_TD=0
							end if
							read(c80tmp(23:32),*) tmpval
							if (tmpval<0) excdir(itmp)=2 !Negative contribution is assumed to be de-excitation (of course this is not strict since -> and <- have been combined together)
							exccoeff(itmp)=dsqrt(abs(tmpval))
						end if
						!Although for closed-shell ground state, ORCA still outputs coefficients as normalization to 1.0, &
						!However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
						if (wfntype==0.or.wfntype==3) exccoeff(itmp)=exccoeff(itmp)/dsqrt(2D0)
					end do
				end if
			end do
        end if
    else if (iORCAsTD==1) then
        !sTDA/sTDDFT ignores very low lying occupied MOs and very high lying virtual MOs, &
        !we need to know how many core orbitals are ignored, so that actual MO index can be obtained
        call loclabel(10,"occ. MOs in")
        read(10,"(a)") c80tmp
        itmp=index(c80tmp,'...')
        if (wfntype==0) then
            read(c80tmp(itmp+3:),*) ntmp
            ncor=nint(naelec)-ntmp
        else
            read(c80tmp(itmp+3:),*) natmp,nbtmp
            nacor=nint(naelec)-natmp
            nbcor=nint(nbelec)-nbtmp
        end if
        call loclabel(10,"state   eV",ifound,0)
        read(10,*)
		do iexc=1,nstates
            read(10,"(a)") c200tmp
            if (iexc==istate) then
                itmp=index(c200tmp,'(')
                ipos=itmp-9
                do ipair=1,excnorb
                    if (wfntype==0) then
                        c80tmp=c200tmp(ipos:ipos+20)
                    else
                        c80tmp=c200tmp(ipos:ipos+22)
                    end if
                    read(c80tmp,*) exccoeff(ipair)
                    if (wfntype==0) exccoeff(ipair)=exccoeff(ipair)/dsqrt(2D0)
                    if (index(c80tmp,'>')/=0) then
                        excdir(ipair)=1
                    else
                        excdir(ipair)=2
                    end if
                    itmp=index(c80tmp,'(') !Load left orbital index
                    if (wfntype==0) then
                        read(c80tmp(itmp+1:itmp+4),*) orbleft(ipair)
                        orbleft(ipair)=orbleft(ipair)+ncor+1  !Note that the MO is 0 based index in sTDA/sTDDFT
                    else
                        read(c80tmp(itmp+1:itmp+4),*) orbleft(ipair)
                        if (c80tmp(itmp+5:itmp+5)=='a') then
                            orbleft(ipair)=orbleft(ipair)+nacor+1
                        else
                            orbleft(ipair)=orbleft(ipair)+nbcor+1+nbasis
                        end if
                    end if
                    itmp=index(c80tmp,')') !Load right orbital index
                    if (wfntype==0) then
                        read(c80tmp(itmp-4:itmp-1),*) orbright(ipair)
                        orbright(ipair)=orbright(ipair)+ncor+1
                    else
                        read(c80tmp(itmp-5:itmp-2),*) orbright(ipair)
                        if (c80tmp(itmp-1:itmp-1)=='a') then
                            orbright(ipair)=orbright(ipair)+nacor+1
                        else
                            orbright(ipair)=orbright(ipair)+nbcor+1+nbasis
                        end if
                    end if
                    if (wfntype==0) then
                        ipos=ipos+21
                    else
                        ipos=ipos+23
                    end if
                end do
            end if
        end do
    end if
	
else if (ifiletypeexc==4) then !Firefly output file
	write(c80tmp,"(' EXCITED STATE',i4)") istate
	call loclabel(10,trim(c80tmp),ifound)
	read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)
	do itmp=1,excnorb
		read(10,*) imo,jmo,exccoeff(itmp)
		if (jmo>imo) then !Excitation
			orbleft(itmp)=imo
			orbright(itmp)=jmo
			excdir(itmp)=1
		else !De-excitation
			orbleft(itmp)=jmo
			orbright(itmp)=imo
			excdir(itmp)=2
		end if
	end do
	!Although for closed-shell reference state, Firefly still outputs coefficients as normalization to 1.0, &
	!However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
	if (wfntype==0.or.wfntype==3) exccoeff=exccoeff/dsqrt(2D0)
	
else if (ifiletypeexc==5) then !GAMESS-US output file
	write(c80tmp,"(' STATE #',i4)") istate
	call loclabel(10,trim(c80tmp),ifound)
	if (wfntype==0) then !Restricted reference state
		read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)
		itmp=0
		do ipair=1,excnorb/2
			read(10,*) imo,jmo,coeff1,coeff2
			itmp=itmp+1
			orbleft(itmp)=imo
			orbright(itmp)=jmo
			excdir(itmp)=1
			exccoeff(itmp)=coeff1
			itmp=itmp+1
			orbleft(itmp)=imo
			orbright(itmp)=jmo
			excdir(itmp)=2
			exccoeff(itmp)=coeff2
		end do
		!Even for closed-shell reference state, GAMESS-US still outputs coefficients as normalization to 1.0, &
		!However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
		if (wfntype==0.or.wfntype==3) exccoeff=exccoeff/dsqrt(2D0)
	else !Unrestricted reference state
		read(10,*);read(10,*);read(10,*);read(10,*);read(10,*);read(10,*)
		itmp=0
		do ipair=1,excnorb/2
			read(10,"(a)") c80tmp
			read(c80tmp,*) imo,jmo
			read(c80tmp(22:48),*) coeff1,coeff2
			ibeta=0
			if (index(c80tmp,"BET")/=0) ibeta=1
			itmp=itmp+1
			orbleft(itmp)=imo
			orbright(itmp)=jmo
			excdir(itmp)=1
			exccoeff(itmp)=coeff1
			itmp=itmp+1
			orbleft(itmp)=imo
			orbright(itmp)=jmo
			excdir(itmp)=2
			exccoeff(itmp)=coeff2
			if (ibeta==1) then
				orbleft(itmp-1:itmp)=orbleft(itmp-1:itmp)+nbasis
				orbright(itmp-1:itmp)=orbright(itmp-1:itmp)+nbasis
			end if
		end do
	end if
    
else if (ifiletypeexc==6) then !CP2K output file
	excdir(:)=1 !CP2K is fully based on TDA
    if (allexcmulti(1)==1) then
		call loclabelfinal(10,"TDDFPT states of multiplicity 1",ifound)
		call loclabel(10,"number             orbital",ifound,0)
    else if (allexcmulti(1)==3) then
		call loclabelfinal(10,"TDDFPT states of multiplicity 3",ifound)
		call loclabel(10,"number             orbital",ifound,0)
	else
		call loclabelfinal(10,"number             orbital",ifound)
    end if
	do iexc=1,istate
        call loclabel(10,"eV",ifound,0)
		read(10,*)
	end do
	do ipair=1,excnorb
		read(10,"(a)") c80tmp
		read(c80tmp(27:34),*) orbleft(ipair)
		read(c80tmp(49:55),*) orbright(ipair)
		read(c80tmp(69:78),*) exccoeff(ipair)
        if (c80tmp(37:39)=="bet") orbleft(ipair)=orbleft(ipair)+nbasis
        if (c80tmp(58:60)=="bet") orbright(ipair)=orbright(ipair)+nbasis
	end do
	!Even for closed-shell reference state, CP2K still outputs coefficients as normalization to 1.0, &
	!However, in order to follow the Gaussian convention, we change the coefficient as normalization to 0.5
	if (wfntype==0.or.wfntype==3) exccoeff=exccoeff/dsqrt(2D0)

else if (ifiletypeexc==7) then !BDF output file. This part of code was contributed by Cong Wang, 2022-Dec-1
    rewind(10)
    do iexc=1,istate
        call loclabel(10,"D<Pab>",ifound,0)
        read(10,*)
    end do
    do itmp=1,excnorb
        !Determine excitation direction
        read(10,"(a)") c200tmp
        ileft=index(c200tmp,"->")
        excdir(itmp)=1 !means ->
        isign=ileft
        if (ileft==0) then
            iright=index(c200tmp,"<-")
            excdir(itmp)=2 !means <-
            isign=iright
        end if
        !Process left side
        leftstr=c200tmp(:isign-1)
        do ileft=len_trim(leftstr),1,-1
            if (leftstr(ileft:ileft)=="(") exit
        end do
        read(leftstr(ileft:len_trim(leftstr)-1),*) orbleft(itmp)
        !Process right side
        rightstr=c200tmp(isign+2:)
        ir1=0
        ir2=0
        do iright=1,len_trim(rightstr)
            if (rightstr(iright:iright)=="(") ir1=iright
            if (rightstr(iright:iright)==")") then
                ir2=iright
                exit
            end if
        end do
        read(rightstr(ir1+1:ir2-1),*) orbright(itmp)
        ir1=0
        ir2=0
        do iright=1,len_trim(rightstr)
            if (rightstr(iright:iright+3)=="Per:") ir1=iright
            if (rightstr(iright:iright)=="%") then
                ir2=iright
                exit
            end if
        end do
        read(rightstr(ir1+4:ir2-1),*) fper
        exccoeff(itmp) = sqrt(fper/200D0)
    end do
end if
close(10)

if (ioutinfo==1) then
	sumsqrexc=0
	sumsqrdeexc=0
	do ipair=1,excnorb
		if (excdir(ipair)==1) sumsqrexc=sumsqrexc+exccoeff(ipair)**2
		if (excdir(ipair)==2) sumsqrdeexc=sumsqrdeexc-exccoeff(ipair)**2
	end do
	sumsqrall=sumsqrexc+sumsqrdeexc
	write(*,"(' Sum of square of excitation coefficients:',f10.6)") sumsqrexc
	write(*,"(' Negative of the sum of square of de-excitation coefficients:',f10.6)") sumsqrdeexc
	write(*,"(' Sum of above two values:',f10.6)") sumsqrall
    if (wfntype==0.or.wfntype==3) then
        dev=abs(sumsqrall-0.5D0)
        write(*,"(' Deviation to expected normalization value (0.5) is',f10.6)") dev
        if (dev>0.05D0) write(*,"(a)") " Warning: The deviation is too obvious, in this case the analysis result is not reliable or even fully misleading!"
    else
        dev=abs(sumsqrall-1D0)
        write(*,"(' Deviation to expected normalization value (1.0) is',f10.6)") dev
        if (dev>0.05D0) write(*,"(a)") " Warning: The deviation is too obvious, in this case the analysis result is not reliable or even fully misleading!"
    end if
    ihighvirMOA=0
    ihighvirMOB=0
    call getHOMOidx
	do ipair=1,excnorb
		iorb=orbright(ipair)
		if (iorb<=nbasis.and.iorb>ihighvirMOA) then
			ihighvirMOA=iorb
        else if (iorb>nbasis.and.iorb>ihighvirMOB) then
			ihighvirMOB=iorb
        end if
		iorb=orbleft(ipair)
		if (iorb<=nbasis.and.iorb>ihighvirMOA) then
			ihighvirMOA=iorb
        else if (iorb>nbasis.and.iorb>ihighvirMOB) then
			ihighvirMOB=iorb
        end if
	end do
	if (wfntype==0.or.wfntype==3) then
		write(*,"(' Involved highest unoccupied orbital:',i6,' (HOMO+',i6,')')") ihighvirMOA,ihighvirMOA-idxHOMO
	else if (wfntype==1.or.wfntype==4) then
		write(*,"(' Involved highest unoccupied alpha orbital:',i6,' (alpha HOMO+',i6,')')") ihighvirMOA,ihighvirMOA-idxHOMO
		write(*,"(' Involved highest unoccupied beta orbital: ',i6,' (beta HOMO+ ',i6,')')") ihighvirMOB-nbasis,ihighvirMOB-idxHOMOb
    end if
end if
end subroutine



!============= The interface used to select an excited state 
subroutine selexcit(istate)
use excitinfo
integer istate
if (nstates>1) then
	write(*,*) "Perform the analysis for which excited state? e.g. 2"
	do while(.true.)
		read(*,*) istate
		if (istate<1) then
			write(*,*) "Error: The index of the selected excited state must be larger than zero! Input again"
		else if (istate>nstates) then
			write(*,"(a,i4,a)") " Error: The index of the selected excited state should <=",nstates," ! Input again"
		else
			exit
		end if
	end do
else
	write(*,*) "The only excited state will be studied"
	istate=1
end if
end subroutine





!=====================================================================================================================
!=====================================================================================================================
!!------- Analyze or visualize hole-electron/transition density/transition dipole moment density distribution --------
!=====================================================================================================================
!=====================================================================================================================
!ROHF,RODFT are assumed to be impossible to be ground state
subroutine hole_electron
use excitinfo
use defvar
use GUI
use util
use functions
implicit real*8 (a-h,o-z)
integer :: idomag=0
real*8 orbval(nmo),wfnderv(3,nmo)
logical,allocatable :: skippair(:) !Record which orbital pairs will be ignored due to negligible coefficient
real*8,allocatable :: holegrid(:,:,:),elegrid(:,:,:),Sm(:,:,:),Sr(:,:,:),transdens(:,:,:),holecross(:,:,:),elecross(:,:,:),Cele(:,:,:),Chole(:,:,:),magtrdens(:,:,:,:)
real*8,allocatable :: cubx(:),cuby(:),cubz(:) !Used to calculate Coulomb attractive energy
character cubsuff*12
iaddstateidx=0
cubsuff=".cub"

if (.not.allocated(CObasa)) then
	write(*,*) "Error: The input file does not contain basis function information!"
	write(*,"(a)") " Please carefully read Section 3.21 of manual to make clear which kinds of input files could be used!"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	return
end if

call loadallexcinfo(1)
call selexcit(istate)
call loadexccoeff(istate,1)
write(*,*)
write(*,"(a)") " !!! Please cite hole-electron analysis as follows, the supplemental material presents a concise introduction of this analysis"
write(*,*) "Tian Lu, et al., Carbon, 165, 461-467 (2020) DOI: 10.1016/j.carbon.2020.05.023"

10 do while(.true.)
	write(*,*)
	write(*,*) "   ------------- Hole-electron and transition density analysis -----------"
	if (idomag==0) write(*,*) "-1 Toggle calculating transit. magnetic dip. density in option 1, current: No"
	if (idomag==1) write(*,*) "-1 Toggle calculating transit. magnetic dip. density in option 1, current: Yes"
	write(*,*) "0 Return"
	write(*,"(a)") " 1 Visualize and analyze hole, electron, transition density and so on"
	write(*,*) "2 Show molecular orbital contribution to hole and electron distribution"
	write(*,"(a)") " 3 Show atom or fragment contribution to hole and electron and plot the contributions as heat map"
	write(*,"(a)") " 4 Show basis function contribution to hole and electron"
	read(*,*) isel
	
	if (isel==-1) then
		if (idomag==0) then
			idomag=1
		else if (idomag==1) then
			idomag=0
		end if
	else if (isel==0) then
		return
	else if (isel==1) then
		exit
	else if (isel==2) then
		call hole_ele_MOcontri
	else if (isel==3) then
		call hole_ele_atmcontri_heatmap
	else if (isel==4) then
		call hole_ele_bascontri
	end if
end do

!Below we will calculate grid data
!Set up grid first
write(*,*)
call setgrid(0,igridsel)
if (allocated(holegrid)) deallocate(holegrid,elegrid,Sm,Sr,transdens,holecross,elecross)
allocate(holegrid(nx,ny,nz),elegrid(nx,ny,nz),Sm(nx,ny,nz),Sr(nx,ny,nz),transdens(nx,ny,nz),holecross(nx,ny,nz),elecross(nx,ny,nz))
holegrid=0D0
elegrid=0D0
holecross=0D0
elecross=0D0
transdens=0D0
if (idomag==1) then !Will also calculate transition magnetic dipole moment density
	if (allocated(magtrdens)) deallocate(magtrdens)
	allocate(magtrdens(nx,ny,nz,3))
	magtrdens=0D0
end if

if (cfgcrossthres/=0) then
	write(*,"(/,a,f8.5,a)") " Note: When calculating cross term of hole and electron, configurations with absolute value of coefficient <",cfgcrossthres,&
    " will be ignored to reduce cost. The threshold is determined by ""cfgcrossthres"" in settings.ini"
end if
allocate(skippair(excnorb))
skippair=.false.
do iexcitorb=1,excnorb
	if (abs(exccoeff(iexcitorb))<0.01D0) skippair(iexcitorb)=.true.
end do

if (ifPBC==0) then
	call gen_GTFuniq(1) !Generate unique GTFs, for faster evaluation
else
	call gen_neigh_GTF !Generate neighbouring GTFs list at reduced grids, for faster evaluation
end if
write(*,*) "Calculating grid data..."
call walltime(iwalltime1)

ifinish=0;ishowprog=1
ntmp=floor(ny*nz/100D0)
!$OMP PARALLEL DO SHARED(ifinish,ishowprog,holegrid,elegrid,transdens,holecross,elecross,magtrdens) &
!$OMP PRIVATE(i,j,k,tmpx,tmpy,tmpz,orbval,wfnderv,imo,jmo,excwei,iexcitorb,jexcitorb,ileft,jleft,iright,jright,tmpleft,tmpright,idir,jdir,tmpval) &
!$OMP schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			if (idomag==1) then !Study transition magnetic dipole moment density requests orbital 1st-derivative
				call orbderv(2,1,nmo,tmpx,tmpy,tmpz,orbval,wfnderv)
			else
				call orbderv(1,1,nmo,tmpx,tmpy,tmpz,orbval)
			end if
			!Calculate local term of hole and electron
            !If i->l,i<-l should be taken into account is unclear, its is ignored currently
			do iexcitorb=1,excnorb
				imo=orbleft(iexcitorb)
				jmo=orbright(iexcitorb)
				excwei=exccoeff(iexcitorb)**2
				if (excdir(iexcitorb)==1) then !Excitation
					holegrid(i,j,k)=holegrid(i,j,k)+excwei*orbval(imo)**2
					elegrid(i,j,k)=elegrid(i,j,k)+excwei*orbval(jmo)**2
				else !De-excitation
					holegrid(i,j,k)=holegrid(i,j,k)-excwei*orbval(imo)**2
					elegrid(i,j,k)=elegrid(i,j,k)-excwei*orbval(jmo)**2
				end if
				transdens(i,j,k)=transdens(i,j,k)+exccoeff(iexcitorb)*orbval(imo)*orbval(jmo)
				if (idomag==1) then
					if (excdir(iexcitorb)==1) then
						magtrdens(i,j,k,1)=magtrdens(i,j,k,1)+exccoeff(iexcitorb)*orbval(imo)*(tmpy*wfnderv(3,jmo)-tmpz*wfnderv(2,jmo))
						magtrdens(i,j,k,2)=magtrdens(i,j,k,2)+exccoeff(iexcitorb)*orbval(imo)*(tmpz*wfnderv(1,jmo)-tmpx*wfnderv(3,jmo))
						magtrdens(i,j,k,3)=magtrdens(i,j,k,3)+exccoeff(iexcitorb)*orbval(imo)*(tmpx*wfnderv(2,jmo)-tmpy*wfnderv(1,jmo))
					else !De-excitation has important influence on the transition magnetic dipole moment density, so must be considered explicitly
						magtrdens(i,j,k,1)=magtrdens(i,j,k,1)-exccoeff(iexcitorb)*orbval(imo)*(tmpy*wfnderv(3,jmo)-tmpz*wfnderv(2,jmo))
						magtrdens(i,j,k,2)=magtrdens(i,j,k,2)-exccoeff(iexcitorb)*orbval(imo)*(tmpz*wfnderv(1,jmo)-tmpx*wfnderv(3,jmo))
						magtrdens(i,j,k,3)=magtrdens(i,j,k,3)-exccoeff(iexcitorb)*orbval(imo)*(tmpx*wfnderv(2,jmo)-tmpy*wfnderv(1,jmo))
					end if
				end if
			end do
			!Calculate cross term of hole and electron. This part takes most majority of computational time
			do iexcitorb=1,excnorb
				!Currently only take below two cases into account:
				! Cross term of hole (do <i|j>):     (i->l,j->l), substract (i<-l,j<-l)
				! Cross term of electron (do <l|m>): (i->l,i->m), substract (i<-l,i<-m)
				!Below cases are skipped:
				! i->l,i->l and i->l,i<-l and i<-l,i<-l, since they do not correspond to cross term
				! i->l,j->m, since differ by two electrons, according to Slater-condon rule the result is zero
				! i->l,i<-m and i<-l,j->l, since excitation directions are different
				if (skippair(iexcitorb)) cycle
				idir=excdir(iexcitorb)
				ileft=orbleft(iexcitorb)
				iright=orbright(iexcitorb)
				tmpleft=exccoeff(iexcitorb)*orbval(ileft) !Use temporary variable to save the time for locating element
				tmpright=exccoeff(iexcitorb)*orbval(iright)
				do jexcitorb=iexcitorb+1,excnorb
					if (skippair(jexcitorb)) cycle
					jdir=excdir(jexcitorb)
					if (idir/=jdir) cycle
					jleft=orbleft(jexcitorb)
					jright=orbright(jexcitorb)
					if (ileft==jleft) then !do <l|m>
						if (iright==jright) cycle
						tmpval=tmpright*exccoeff(jexcitorb)*orbval(jright) !Originally virtual orbital
						if (idir==1) then !->
							elecross(i,j,k)=elecross(i,j,k)+tmpval
						else !<-
							elecross(i,j,k)=elecross(i,j,k)-tmpval
						end if
					else if (iright==jright) then !do <i|j>
						tmpval=tmpleft*exccoeff(jexcitorb)*orbval(jleft) !Originally occupied orbital
						if (idir==1) then !->
							holecross(i,j,k)=holecross(i,j,k)+tmpval
						else !<-
							holecross(i,j,k)=holecross(i,j,k)-tmpval
						end if
					end if
				end do
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

holecross=holecross*2 !Because "jexcitorb=iexcitorb+1,excnorb" does not considered duplicated cross terms
elecross=elecross*2
deallocate(skippair)
call del_GTFuniq !Destory unique GTF informtaion

if (wfntype==0.or.wfntype==3) then !For closed-shell wavefunction, the weights are normalized to 0.5 (or say the orbitals are doubly occupied), so correct it
	holegrid=holegrid*2
	elegrid=elegrid*2
	transdens=transdens*2
	holecross=holecross*2
	elecross=elecross*2
	if (idomag==1) magtrdens=magtrdens*2
end if
!Combine local term and cross term of hole to holegrid, that of electron to elegrid. Then local term will be holegrid-holecross and elegrid-elecross
holegrid=holegrid+holecross
elegrid=elegrid+elecross
! write(*,*) minval(holegrid),minval(elegrid) !In some places, holegrid or elegrid may be negative due to numerical inaccuracy, here examining its magnitude
where (holegrid<0) holegrid=0 !Force hole and electron to be non-negative, otherwise Sr cannot be calculated
where (elegrid<0) elegrid=0
Sm=min(holegrid,elegrid)
Sr=dsqrt(holegrid*elegrid)

call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,' s',/)") iwalltime2-iwalltime1

!Check normalization
call calc_dvol(dvol)
rnormhole=0
rnormele=0
Sm_index=0
Sr_index=0
rinttransdens=0
rtransdipx=0
rtransdipy=0
rtransdipz=0
centholex=0
centholey=0
centholez=0
centelex=0
centeley=0
centelez=0
rtransmagx=0
rtransmagy=0
rtransmagz=0
ele_deloc=0
hole_deloc=0
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			rnormhole=rnormhole+holegrid(i,j,k)
			rnormele=rnormele+elegrid(i,j,k)
			Sm_index=Sm_index+Sm(i,j,k)
			Sr_index=Sr_index+Sr(i,j,k)
			rinttransdens=rinttransdens+transdens(i,j,k)
			rtransdipx=rtransdipx-tmpx*transdens(i,j,k)
			rtransdipy=rtransdipy-tmpy*transdens(i,j,k)
			rtransdipz=rtransdipz-tmpz*transdens(i,j,k)
			centholex=centholex+holegrid(i,j,k)*tmpx
			centholey=centholey+holegrid(i,j,k)*tmpy
			centholez=centholez+holegrid(i,j,k)*tmpz
			centelex=centelex+elegrid(i,j,k)*tmpx
			centeley=centeley+elegrid(i,j,k)*tmpy
			centelez=centelez+elegrid(i,j,k)*tmpz
            ele_deloc=ele_deloc+elegrid(i,j,k)**2
            hole_deloc=hole_deloc+holegrid(i,j,k)**2
			if (idomag==1) then
				rtransmagx=rtransmagx+magtrdens(i,j,k,1)
				rtransmagy=rtransmagy+magtrdens(i,j,k,2)
				rtransmagz=rtransmagz+magtrdens(i,j,k,3)
			end if
		end do
	end do
end do
rnormhole=rnormhole*dvol
rnormele=rnormele*dvol
Sm_index=Sm_index*dvol
Sr_index=Sr_index*dvol
rinttransdens=rinttransdens*dvol
rtransdipx=rtransdipx*dvol
rtransdipy=rtransdipy*dvol
rtransdipz=rtransdipz*dvol
centholex=centholex*dvol/rnormhole
centholey=centholey*dvol/rnormhole
centholez=centholez*dvol/rnormhole
centelex=centelex*dvol/rnormele
centeley=centeley*dvol/rnormele
centelez=centelez*dvol/rnormele
ele_deloc=dsqrt(ele_deloc*dvol)*100
hole_deloc=dsqrt(hole_deloc*dvol)*100
rtransmagx=rtransmagx*dvol
rtransmagy=rtransmagy*dvol
rtransmagz=rtransmagz*dvol
if (excmulti==3) write(*,"(a,/)") " Note: The transition moments shown below correspond to spatial part, while orthogonality of spin part of wavefunction is ignored"

if (rnormele==0) then
	write(*,"(a,/)") " Warning: Integral of electron distribution is exactly zero, it is highly likely that &
    &wavefunctions of virtual orbitals were not recorded in your inputted wavefunction file! In this case hole-electron analysis is meaningless"
end if

write(*,"(' Integral of hole:    ',f12.6)") rnormhole
write(*,"(' Integral of electron:',f12.6)") rnormele
write(*,"(' Integral of transition density:',f12.6)") rinttransdens
write(*,"(' Transition dipole moment in X/Y/Z:',3f11.6,' a.u.')") rtransdipx,rtransdipy,rtransdipz
if (idomag==1) write(*,"(' Transition magnetic dipole moment in X/Y/Z:',3f10.6,' a.u.')") rtransmagx,rtransmagy,rtransmagz
write(*,"(' Sm index (integral of Sm function):',f10.5,' a.u.')") Sm_index
write(*,"(' Sr index (integral of Sr function):',f10.5,' a.u.')") Sr_index
write(*,"(' Centroid of hole in X/Y/Z:    ',3f12.6,' Angstrom')") centholex*b2a,centholey*b2a,centholez*b2a
write(*,"(' Centroid of electron in X/Y/Z:',3f12.6,' Angstrom')") centelex*b2a,centeley*b2a,centelez*b2a
disx=abs(centelex-centholex)
disy=abs(centeley-centholey)
disz=abs(centelez-centholez)
disnorm=dsqrt(disx**2+disy**2+disz**2)
write(*,"(' D_x:',f8.3,'  D_y:',f8.3,'  D_z:',f8.3,'    D index:',f8.3,' Angstrom')") disx*b2a,disy*b2a,disz*b2a,disnorm*b2a
write(*,"(' Variation of dipole moment with respect to ground state:')")
avgtransval=(rnormele+rnormhole)/2D0 !Ideal value is 1.0
dipvarx=-(centelex-centholex)*avgtransval
dipvary=-(centeley-centholey)*avgtransval
dipvarz=-(centelez-centholez)*avgtransval
dipvarnorm=dsqrt(dipvarx**2+dipvary**2+dipvarz**2)
write(*,"(' X:',f12.6,'  Y:',f12.6,'  Z:',f12.6,'    Norm:',f12.6,' a.u.')") dipvarx,dipvary,dipvarz,dipvarnorm

sigxele=0D0;sigyele=0D0;sigzele=0D0
sigxhole=0D0;sigyhole=0D0;sigzhole=0D0
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			sigxele=sigxele+elegrid(i,j,k)*(tmpx-centelex)**2
			sigyele=sigyele+elegrid(i,j,k)*(tmpy-centeley)**2
			sigzele=sigzele+elegrid(i,j,k)*(tmpz-centelez)**2
			sigxhole=sigxhole+holegrid(i,j,k)*(tmpx-centholex)**2
			sigyhole=sigyhole+holegrid(i,j,k)*(tmpy-centholey)**2
			sigzhole=sigzhole+holegrid(i,j,k)*(tmpz-centholez)**2
		end do
	end do
end do
sigxele=dsqrt(sigxele*dvol/rnormele)
sigyele=dsqrt(sigyele*dvol/rnormele)
sigzele=dsqrt(sigzele*dvol/rnormele)
signormele=dsqrt(sigxele**2+sigyele**2+sigzele**2)
sigxhole=dsqrt(sigxhole*dvol/rnormhole)
sigyhole=dsqrt(sigyhole*dvol/rnormhole)
sigzhole=dsqrt(sigzhole*dvol/rnormhole)
signormhole=dsqrt(sigxhole**2+sigyhole**2+sigzhole**2)
write(*,"(' RMSD of hole in X/Y/Z:    ',3f8.3,'   Norm:',f8.3,' Angstrom')") sigxhole*b2a,sigyhole*b2a,sigzhole*b2a,signormhole*b2a
write(*,"(' RMSD of electron in X/Y/Z:',3f8.3,'   Norm:',f8.3,' Angstrom')") sigxele*b2a,sigyele*b2a,sigzele*b2a,signormele*b2a
delta_sigma_x=sigxele-sigxhole
delta_sigma_y=sigyele-sigyhole
delta_sigma_z=sigzele-sigzhole
delta_sigma_index=signormele-signormhole
write(*,*) "Difference between RMSD of hole and electron (delta sigma):"
write(*,"(' X:',f7.3,'  Y:',f7.3,'  Z:',f7.3,'    Overall:',f7.3,' Angstrom')") delta_sigma_x*b2a,delta_sigma_y*b2a,delta_sigma_z*b2a,delta_sigma_index*b2a
Hx=(sigxele+sigxhole)/2D0
Hy=(sigyele+sigyhole)/2D0
Hz=(sigzele+sigzhole)/2D0
H_CT= ( dsqrt ( (Hx*(centelex-centholex))**2 + (Hy*(centeley-centholey))**2 + (Hz*(centelez-centholez))**2 ) ) / disnorm
H_index=(signormele+signormhole)/2D0
write(*,"(' H_x:',f7.3,'  H_y:',f7.3,'  H_z:',f7.3,'  H_CT:',f7.3,'  H index:',f7.3,' Angstrom')") Hx*b2a,Hy*b2a,Hz*b2a,H_CT*b2a,H_index*b2a
t_index=disnorm-H_CT
write(*,"(' t index:',f7.3,' Angstrom')") t_index*b2a
write(*,"(' Hole delocalization index (HDI):    ',f7.2)") hole_deloc
write(*,"(' Electron delocalization index (EDI):',f7.2)") ele_deloc

!---- Calculate ghost state diagnostic index proposed by Adamo (JCC,38,2151) and modified by Tian Lu
sumCsqr=0
do itmp=1,excnorb
	if (excdir(itmp)==2) cycle !Skip de-excitation
	sumCsqr=sumCsqr+exccoeff(itmp)**2
end do
ghostp1=0
do itmp=1,excnorb
	if (excdir(itmp)==2) cycle !Skip de-excitation
	ghostp1=ghostp1+exccoeff(itmp)**2/sumCsqr*(MOene(orbright(itmp))-MOene(orbleft(itmp)))
end do
ghostp2=1/disnorm !in a.u.
ghostidx=ghostp1-ghostp2
write(*,"(' Ghost-hunter index:',f11.3,' eV, 1st term:',f7.3,' eV, 2nd term:',f10.3,' eV')") ghostidx*au2eV,ghostp1*au2eV,ghostp2*au2eV
write(*,"(' Excitation energy of this state:',f10.3,' eV')") excene
if (excene<ghostidx*au2eV) then
	write(*,"(a)") " Note: Probably this is a ghost state, because excitation energy is lower than ghost-hunter index"
    write(*,"(a)") " Comment by Multiwfn author: It is frequently found that the condition of determining ghost state by the ghost-hunder &
    &index is too stringent, it is only suitable to be viewed as a necessary condition of ghost state rather than a sufficient condition. If you used a &
    &DFT functional having high Hartree-Fock composition in long range of electronic interaction, such as wB97XD, CAM-B3LYP, and M06-2X, &
    &you do not need to worry about presence of ghost state at all. In addition, ghost state should have excitation wavelength around 1000 nm or more and have &
    &negligible oscillator strength, if this is not the current case, you also do not to worry about the possibility of occurrence of ghost state"
end if

if (allocated(Cele)) deallocate(Cele,Chole)
allocate(Cele(nx,ny,nz),Chole(nx,ny,nz))
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,rnowx,rnowy,rnowz)
			Cele(i,j,k)=exp( -(rnowx-centelex)**2/(2*sigxele**2) -(rnowy-centeley)**2/(2*sigyele**2) -(rnowz-centelez)**2/(2*sigzele**2))
			Chole(i,j,k)=exp( -(rnowx-centholex)**2/(2*sigxhole**2) -(rnowy-centholey)**2/(2*sigyhole**2) -(rnowz-centholez)**2/(2*sigzhole**2))
		end do
	end do
end do
Cele=Cele*rnormele/(sum(Cele)*dvol)
Chole=Chole*rnormhole/(sum(Chole)*dvol)

if (ifPBC>0) write(*,"(/,a)") " Note: RMSD, t, H and ghost-hunter indices, as well as Cele and Chole, may be meaningless if your system is periodic!"

do while(.true.)
	write(*,*)
	write(*,*) "             -------------- Post-processing menu --------------"
    if (iaddstateidx==1) write(*,*) "-1 Toggle if adding state index at end of exported cube filename, current: Yes"
    if (iaddstateidx==0) write(*,*) "-1 Toggle if adding state index at end of exported cube filename, current: No"
	write(*,*) "0 Return"
	write(*,*) "1 Show isosurface of hole distribution"
	write(*,*) "2 Show isosurface of electron distribution"
	write(*,*) "3 Show isosurface of hole and electron distribution simultaneously"
	write(*,*) "4 Show isosurface of overlap of hole-electron (i.e. Sm or Sr function)"
	write(*,*) "5 Show isosurface of transition density"
	write(*,*) "6 Show isosurface of transition dipole moment density"
	write(*,*) "7 Show isosurface of charge density difference"
	write(*,*) "8 Show isosurface of Cele and Chole functions simultaneously"
	write(*,*) "9 Show isosurface of transition magnetic dipole moment density"
	write(*,*) "10 Output cube file of hole distribution to current folder"
	write(*,*) "11 Output cube file of electron distribution to current folder"
	write(*,*) "12 Output cube file of overlap of hole-electron to current folder"
	write(*,*) "13 Output cube file of transition density to current folder"
	write(*,*) "14 Output cube file of transition dipole moment density to current folder"
	write(*,*) "15 Output cube file of charge density difference to current folder"
	write(*,*) "16 Output cube file of Cele and Chole functions to current folder"
	write(*,*) "17 Output cube file of transition magnetic dipole moment density"
	write(*,*) "18 Calculate hole-electron Coulomb attractive energy (exciton binding energy)"
	read(*,*) isel
    if (isel==-1) then
        if (iaddstateidx==0) then
            iaddstateidx=1
            write(cubsuff,"('_',i5.5,'.cub')") istate
        else if (iaddstateidx==1) then
            iaddstateidx=0
            cubsuff=".cub"
        end if
	else if (isel==0) then
		goto 10
	else if (isel==1.or.isel==2.or.isel==4.or.isel==5.or.isel==6.or.isel==7.or.isel==9) then
	 	if (allocated(cubmat)) deallocate(cubmat)
	 	allocate(cubmat(nx,ny,nz))
	 	if (isel==1) then
! 	 		write(*,*) "Select the type of hole distribution to be shown. 1 is commonly used"
! 	 		write(*,*) "1 Total (local term + cross term)"
! 	 		write(*,*) "2 Local term only"
! 	 		write(*,*) "3 Cross term only"
! 	 		read(*,*) isel2
! 	 		if (isel2==1) then
! 	 			cubmat=holegrid
! 	 		else if (isel2==2) then
! 	 			cubmat=holegrid-holecross
! 	 		else if (isel2==3) then
! 	 			cubmat=holecross
! 	 		end if
	 		cubmat=holegrid
	 		sur_value=0.002D0
	 	else if (isel==2) then
! 	 		write(*,*) "Select the type of electron distribution to be shown. 1 is commonly used"
! 	 		write(*,*) "1 Total (local term + cross term)"
! 	 		write(*,*) "2 Local term only"
! 	 		write(*,*) "3 Cross term only"
! 	 		read(*,*) isel2
! 	 		if (isel2==1) then
! 	 			cubmat=elegrid
! 	 		else if (isel2==2) then
! 	 			cubmat=elegrid-elecross
! 	 		else if (isel2==3) then
! 	 			cubmat=elecross
! 	 		end if
	 		cubmat=elegrid
	 		sur_value=0.002D0
	 	else if (isel==4) then
	 		write(*,*) "Select the function that measures overlap of hole and electron distribution"
	 		write(*,*) "1 Sm function"
	 		write(*,*) "2 Sr function (recommended)"
	 		read(*,*) iovlptype
	 		if (iovlptype==1) then
	 			cubmat=Sm
	 			sur_value=0.002D0
	 		else if (iovlptype==2) then
	 			cubmat=Sr
	 			sur_value=0.002D0
	 		end if
	 	else if (isel==5) then
	 		cubmat=transdens
			isosur1style=4
	 		sur_value=0.001D0
		else if (isel==6) then
			write(*,*) "Select the component of transition dipole moment density"
			write(*,*) "1: X component  2: Y component  3: Z component  4: Norm, sqrt(x^2+y^2+z^2)"
	 		read(*,*) ifac
			isosur1style=4
	 		sur_value=0.001D0
	 		do k=1,nz
				do j=1,ny
					do i=1,nx
                        call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
						if (ifac==1) cubmat(i,j,k)=-tmpx*transdens(i,j,k)
						if (ifac==2) cubmat(i,j,k)=-tmpy*transdens(i,j,k)
						if (ifac==3) cubmat(i,j,k)=-tmpz*transdens(i,j,k)
                        if (ifac==4) cubmat(i,j,k)=dsqrt(tmpx**2+tmpy**2+tmpz**2)*abs(transdens(i,j,k))
					end do
				end do
			end do
	 	else if (isel==7) then
	 		cubmat=elegrid-holegrid
	 		sur_value=0.002D0
	 	else if (isel==9) then
            if (idomag==0) then
                write(*,"(a)") " Error: In order to visualize transition magnetic dipole moment density, you must select &
                &""-1 Toggle calculating transit. magnetic dip. density in option 1"" to switch status to ""Yes"" before calculating grid data!"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                cycle
            end if
			write(*,*) "Select the component of transition magnetic dipole moment density"
			write(*,*) "1: X component  2: Y component  3: Z component  4: Norm, sqrt(x^2+y^2+z^2)"
	 		read(*,*) ifac
			isosur1style=4
	 		sur_value=0.007D0
			if (ifac==1) cubmat=magtrdens(:,:,:,1)
			if (ifac==2) cubmat=magtrdens(:,:,:,2)
			if (ifac==3) cubmat=magtrdens(:,:,:,3)
			if (ifac==4) cubmat=dsqrt(magtrdens(:,:,:,1)**2+magtrdens(:,:,:,2)**2+magtrdens(:,:,:,3)**2)
	 	end if
		call drawisosurgui(1)
		deallocate(cubmat)
	else if (isel==3.or.isel==8) then
		if (isel==3) write(*,"(a)") " Note: Blue and green isosurfaces represent hole and electron distributions, respectively"
		if (isel==8) write(*,"(a)") " Note: Blue and green isosurfaces represent Chole and Cele functions, respectively"
		sur_value=0.002D0
		isosursec=1
		clrRcub2sameold=clrRcub2same !Backup previous color setting
		clrGcub2sameold=clrGcub2same
		clrBcub2sameold=clrBcub2same
		clrRcub2samemeshptold=clrRcub2samemeshpt
		clrGcub2samemeshptold=clrGcub2samemeshpt
		clrBcub2samemeshptold=clrBcub2samemeshpt
		clrRcub2same=0.3D0
		clrGcub2same=0.45D0
		clrBcub2same=0.9D0
		clrRcub2samemeshpt=0.3D0
		clrGcub2samemeshpt=0.45D0
		clrBcub2samemeshpt=0.9D0
		if (allocated(cubmat)) deallocate(cubmat)
		if (allocated(cubmattmp)) deallocate(cubmattmp)
		allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
		if (isel==3) then
			cubmat=elegrid
			cubmattmp=holegrid
		else if (isel==8) then
			cubmat=Cele
			cubmattmp=Chole
		end if
		isosur1style=4
		isosur2style=4
		call drawisosurgui(2)
		deallocate(cubmat,cubmattmp)
		clrRcub2same=clrRcub2sameold !Recover previous color setting
		clrGcub2same=clrGcub2sameold
		clrBcub2same=clrBcub2sameold
		clrRcub2samemeshpt=clrRcub2samemeshptold
		clrGcub2samemeshpt=clrGcub2samemeshptold
		clrBcub2samemeshpt=clrBcub2samemeshptold
	else if (isel==10) then
 		write(*,*) "Select the type of hole distribution to be exported. 1 is commonly used"
 		write(*,*) "1 Total (local term + cross term)"
 		write(*,*) "2 Local term only"
 		write(*,*) "3 Cross term only"
 		read(*,*) isel2
		write(*,*) "Outputting hole distribution to hole"//trim(cubsuff)//" in current folder"
		open(10,file="hole"//trim(cubsuff),status="replace")
		if (isel2==1) call outcube(holegrid,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		if (isel2==2) call outcube(holegrid-holecross,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		if (isel2==3) call outcube(holecross,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
	else if (isel==11) then
 		write(*,*) "Select the type of electron distribution to be exported. 1 is commonly used"
 		write(*,*) "1 Total (local term + cross term)"
 		write(*,*) "2 Local term only"
 		write(*,*) "3 Cross term only"
 		read(*,*) isel2
		write(*,*) "Outputting electron distribution to electron"//trim(cubsuff)//" in current folder"
		open(10,file="electron"//trim(cubsuff),status="replace")
		if (isel2==1) call outcube(elegrid,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		if (isel2==2) call outcube(elegrid-elecross,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		if (isel2==3) call outcube(elecross,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
	else if (isel==12) then
 		write(*,*) "Select the function that measures overlap of hole and electron distribution"
 		write(*,*) "1 Sm function"
 		write(*,*) "2 Sr function (recommended)"
 		read(*,*) iovlptype
 		if (iovlptype==1) then
			write(*,*) "Outputting Sm function to Sm"//trim(cubsuff)//" in current folder"
			open(10,file="Sm"//trim(cubsuff),status="replace")
			call outcube(Sm,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
 		else if (iovlptype==2) then
			write(*,*) "Outputting Sr function to Sr"//trim(cubsuff)//" in current folder"
			open(10,file="Sr"//trim(cubsuff),status="replace")
			call outcube(Sr,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
 		end if
		close(10)
		write(*,*) "Done!"
	else if (isel==13) then
		write(*,"(a)") " Outputting transition density to transdens"//trim(cubsuff)//" in current folder"
		open(10,file="transdens"//trim(cubsuff),status="replace")
		call outcube(transdens,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
 	else if (isel==14) then
		write(*,*) "Select the component of transition dipole moment density"
		write(*,*) "1: X component  2: Y component  3: Z component  4: Norm, sqrt(x^2+y^2+z^2)"
 		read(*,*) ifac
 		if (allocated(cubmat)) deallocate(cubmat)
	 	allocate(cubmat(nx,ny,nz))
 		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					if (ifac==1) cubmat(i,j,k)=-tmpx*transdens(i,j,k)
					if (ifac==2) cubmat(i,j,k)=-tmpy*transdens(i,j,k)
					if (ifac==3) cubmat(i,j,k)=-tmpz*transdens(i,j,k)
					if (ifac==4) cubmat(i,j,k)=dsqrt(tmpx**2+tmpy**2+tmpz**2)*abs(transdens(i,j,k))
				end do
			end do
		end do
		write(*,*) "Outputting to transdipdens"//trim(cubsuff)//" in current folder"
		open(10,file="transdipdens"//trim(cubsuff),status="replace")
		call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
		deallocate(cubmat)
	else if (isel==15) then
		write(*,*) "Outputting charge density difference to CDD"//trim(cubsuff)//" in current folder"
		open(10,file="CDD"//trim(cubsuff),status="replace")
		call outcube(elegrid-holegrid,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
	else if (isel==16) then
		open(10,file="Cele"//trim(cubsuff),status="replace")
		call outcube(Cele,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,"(' Cele function has been outputted to Cele"//trim(cubsuff)//" in current folder')")
		open(10,file="Chole"//trim(cubsuff),status="replace")
		call outcube(Chole,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,"(' Chole function has been outputted to Chole"//trim(cubsuff)//" in current folder')")
 	else if (isel==17) then
        if (idomag==0) then
            write(*,"(a)") " Error: In order to export transition magnetic dipole moment density, you must select &
            &""-1 Toggle calculating transit. magnetic dip. density in option 1"" to switch status to ""Yes"" before calculating grid data!"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            cycle
        end if
		write(*,*) "Select the component of transition magnetic dipole moment density"
		write(*,*) "1: X component  2: Y component  3: Z component  4: Norm, sqrt(x^2+y^2+z^2)"
 		read(*,*) ifac
		write(*,*) "Outputting to magtrdipdens"//trim(cubsuff)//" in current folder"
		open(10,file="magtrdipdens"//trim(cubsuff),status="replace")
		if (ifac==1) call outcube(magtrdens(:,:,:,1),nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		if (ifac==2) call outcube(magtrdens(:,:,:,2),nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		if (ifac==3) call outcube(magtrdens(:,:,:,3),nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		if (ifac==4) call outcube(dsqrt(magtrdens(:,:,:,1)**2+magtrdens(:,:,:,2)**2+magtrdens(:,:,:,3)**2),nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
 	else if (isel==18) then
		call walltime(iwalltime1)
		write(*,*) "Calculating, please wait..."
		allocate(cubx(nx),cuby(ny),cubz(nz))
		do i=1,nx
			cubx(i)=orgx+(i-1)*dx
		end do
		do i=1,ny
			cuby(i)=orgy+(i-1)*dy
		end do
		do i=1,nz
			cubz(i)=orgz+(i-1)*dz
		end do
 		coulene=0
		do k=1,nz
            cubzk=cubz(k)
			do j=1,ny
                cubyj=cuby(j)
				do i=1,nx
					tmpval=Cele(i,j,k)
                    cubxi=cubx(i)
					if (Cele(i,j,k)<1D-6) cycle !Typically leads to error at 0.001 magnitude
					!$OMP parallel shared(coulene) private(ii,jj,kk,distx2,disty2,distz2,dist,coulenetmp) num_threads(nthreads)
					coulenetmp=0
					!$OMP do schedule(DYNAMIC)
					do kk=1,nz
                        distz2=(cubzk-cubz(kk))**2
						do jj=1,ny
							disty2=(cubyj-cuby(jj))**2
							do ii=1,nx
								if (i==ii.and.j==jj.and.k==kk) cycle
								distx2=(cubxi-cubx(ii))**2
								dist=dsqrt(distx2+disty2+distz2)
								coulenetmp=coulenetmp+tmpval*Chole(ii,jj,kk)/dist
							end do
						end do
					end do
					!$OMP END DO
					!$OMP CRITICAL
					coulene=coulene+coulenetmp
					!$OMP END CRITICAL
					!$OMP END PARALLEL
				end do
			end do
			call showprog(k,nz)
		end do
		call calc_dvol(dvol)
		coulene=coulene*dvol*dvol
		call walltime(iwalltime2)
		write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
		write(*,*)
		write(*,"(' Coulomb attractive energy:',f12.6,' a.u.  (',f12.6,' eV )')") coulene,coulene*au2eV
		deallocate(cubx,cuby,cubz)
	end if
end do
end subroutine



!--------------------------------------------------------------------------------------
!----------- Show contribution of each MO to hole and electron distribution -----------
!Invoked by hole-electron analysis module
subroutine hole_ele_MOcontri
use defvar
use excitinfo
implicit real*8 (a-h,o-z)
real*8,allocatable :: tmparr1(:),tmparr2(:) !Arrays for temporary use

write(*,*) "Input the printing criterion"
write(*,*) "e.g. 0.5 means only the MOs having contribution >=0.5% will be printed"
read(*,*) printcrit
printcrit=printcrit/100
allocate(tmparr1(nmo),tmparr2(nmo))
tmparr1=0 !Record contribution of each MO to hole
tmparr2=0 !Record contribution of each MO to electron
do iexcitorb=1,excnorb !Cycle each excitation pair
	imo=orbleft(iexcitorb)
	jmo=orbright(iexcitorb)
	excwei=exccoeff(iexcitorb)**2
	if (excdir(iexcitorb)==1) then ! ->
		tmparr1(imo)=tmparr1(imo)+excwei
		tmparr2(jmo)=tmparr2(jmo)+excwei
	else ! <-
		tmparr1(imo)=tmparr1(imo)-excwei
		tmparr2(jmo)=tmparr2(jmo)-excwei
	end if
end do
if (wfntype==0.or.wfntype==3) then !Ground state is closed-shell
	tmparr1=tmparr1*2
	tmparr2=tmparr2*2
	do imo=1,nmo
		if (tmparr1(imo)>=printcrit.or.tmparr2(imo)>=printcrit) &
		write(*,"(' MO',i5,', Occ:',f10.5,'    Hole:',f9.3,' %    Electron:',f9.3,' %')") imo,MOocc(imo),tmparr1(imo)*100,tmparr2(imo)*100
	end do
else !Unrestricted reference state
	do imo=1,nmo
		if (tmparr1(imo)>=printcrit.or.tmparr2(imo)>=printcrit) then
			if (imo<=nbasis) then
				write(*,"(' MO',i5,'A, Occ:',f10.5,'    Hole:',f9.3,' %    Electron:',f9.3,' %')") imo,MOocc(imo),tmparr1(imo)*100,tmparr2(imo)*100
			else
				write(*,"(' MO',i5,'B, Occ:',f10.5,'    Hole:',f9.3,' %    Electron:',f9.3,' %')") imo-nbasis,MOocc(imo),tmparr1(imo)*100,tmparr2(imo)*100
			end if
		end if
	end do
end if
write(*,"(' Sum of hole:',f9.3,' %    Sum of electron:',f9.3,' %')") sum(tmparr1)*100,sum(tmparr2)*100
end subroutine





!------------------------------------------------------------------------------------------------------
! Calculate contribution of atom and fragment to hole and electron distribution and plot it as heat map
!------------------------------------------------------------------------------------------------------
subroutine hole_ele_atmcontri_heatmap
use defvar
use util
use dislin
implicit real*8 (a-h,o-z)
character c200tmp*200,c2000tmp*2000
real*8 :: atmcomp(ncenter,nmo),xyratio=0.15D0
real*8 he_atm(ncenter,3) !hole-electron composition of atoms. column 1,2,3 = hole,electron,overlap
real*8,allocatable :: he_atmnoh(:,:) !Similar to he_atm, but hydrogens are not taken into accounts
real*8,allocatable :: he_frag(:,:) !hole-electron composition of fragments. row 1,2,3 = hole,electron,overlap
integer,allocatable :: frag(:,:),fragnatm(:),nohlist(:)

write(*,*)
write(*,*) "Select the method for calculating fragment contributions to hole and electron"
write(*,*) "1 Mulliken (very fast, but incompatible with diffuse functions)"
write(*,*) "2 Hirshfeld partition (slower but very robust)"
read(*,*) imethod

if (imethod==1) then
	call atmcontri_holeele_Mulliken(he_atm(:,1),he_atm(:,2),1) !Calculate atomic contribution to hole and electron via Mulliken-like partition
else if (imethod==2) then
	call atmcontri_holeele_Hirshfeld(he_atm(:,1),he_atm(:,2),1) !Calculate atomic contribution to hole and electron via Hirshfeld partition
end if

do iatm=1,ncenter !Evaluate overlap in atom space. If contribution to hole or electron is unphysical negative, then the overlap will be regarded as zero
! 	he_atm(iatm,3)=minval(he_atm(iatm,1:2)) !Not as good as below
	if (he_atm(iatm,1)>0.and.he_atm(iatm,2)>0) he_atm(iatm,3)=dsqrt(he_atm(iatm,1)*he_atm(iatm,2))
end do

!Contract he_atm to he_atmnoh, which is he_atm without hydrogens
ncennoh=count(a(:)%index/=1)
write(*,"(' The number of non-hydrogen atoms:',i10)") ncennoh
allocate(he_atmnoh(ncennoh,3),nohlist(ncennoh))
itmp=0
do iatm=1,ncenter
	if (a(iatm)%index/=1) then
		itmp=itmp+1
		nohlist(itmp)=iatm
		he_atmnoh(itmp,:)=he_atm(iatm,:)
	end if
end do

write(*,*) "Contribution of each non-hydrogen atom to hole and electron:"
do idx=1,ncennoh
	iatm=nohlist(idx)
	write(*,"(i5,'(',a,')  Hole:',f6.2,' %  Electron:',f6.2,' %  Overlap:',f6.2,' %  Diff.:',f7.2,' %')") &
	iatm,a(iatm)%name,he_atmnoh(idx,:)*100,(he_atmnoh(idx,2)-he_atmnoh(idx,1))*100
end do
write(*,"(' Sum of hole shown above:',f7.2,'%   Sum of electron shown above:',f7.2,'%')") sum(he_atmnoh(:,1))*100,sum(he_atmnoh(:,2))*100

ifhydrogen=0 !By default, hydrogens are not taken into account
clrlimlow=0
clrlimhigh=maxval(he_atmnoh)
stepz=(clrlimhigh-clrlimlow)/5D0
nstepsize=ceiling(ncennoh/10D0)
do while(.true.)
	write(*,*)
	write(*,*) "     --------------- Plotting hole/electron composition ---------------"
	if (allocated(frag)) then
		write(*,"(a,i3,a)") " -1 Reload or clean fragment definition, current:",nfrag," fragments"
	else
		write(*,"(a)") " -1 Load fragment definition"
	end if
	write(*,*) "0 Return"
	write(*,*) "1 Plot hole/electron composition as heat map"
	write(*,*) "2 Export the data as he_atm.txt in current folder"
	write(*,*) "3 Save the heat map as graphic file in current folder"
	write(*,"(a,i3)") " 4 Set interval between labels in X axis, current:",nstepsize
	write(*,"(a,f7.3,a,f7.3)") " 5 Change lower and upper limit of color scale, current:",clrlimlow," to",clrlimhigh
	write(*,"(' 6 Set ratio between X and Y axes, current:',f5.2)") xyratio
	if (.not.allocated(frag)) then
		if (ifhydrogen==0) write(*,*) "7 Toggle if taking hydrogens into account, current: No"
		if (ifhydrogen==1) write(*,*) "7 Toggle if taking hydrogens into account, current: Yes"
	end if
	read(*,*) isel
	if (isel==0) then
		return
		
	else if (isel==-1) then
		write(*,*)
		write(*,*) "How many fragments to be defined? e.g. 3"
		write(*,"(a)") " Note: If you input 0, then fragment definition will be loaded from an external plain text file"
		if (allocated(frag)) write(*,*) "If input ""clean"", existing fragments will be discarded"
		read(*,"(a)") c200tmp
		if (c200tmp=="clean") then !Restore to non-hydrogen atom based case
			deallocate(frag,fragnatm,he_frag)
			nfrag=0
			ifhydrogen=0
			clrlimlow=0
			clrlimhigh=maxval(he_atmnoh)
			stepz=(clrlimhigh-clrlimlow)/5D0
			xyratio=0.15D0
			nstepsize=ceiling(ncennoh/10D0)
			cycle
		else
			if (allocated(frag)) deallocate(frag,fragnatm,he_frag)
			read(c200tmp,*) nfrag
			if (nfrag==0) then
				write(*,*) "Input the file containing fragment definition, e.g. C:\fragdef.txt"
				write(*,"(a)") " Note: If pressing ENTER button directly, fragdef.txt in current folder will be used"
				do while(.true.)
					read(*,"(a)") c200tmp
					if (c200tmp==" ") c200tmp="fragdef.txt"
					inquire(file=c200tmp,exist=alive)
					if (alive) then
						open(10,file=c200tmp,status="old")
						nfrag=0
						do while(.true.)
							read(10,"(a)",iostat=ierror) c200tmp
							if (c200tmp==" ".or.ierror/=0) exit
							nfrag=nfrag+1
						end do
						write(*,"(' There are',i5,' fragments')") nfrag
						allocate(frag(nfrag,ncenter),fragnatm(nfrag))
						rewind(10)
						do ifrag=1,nfrag
							read(10,"(a)") c200tmp
							call str2arr(c200tmp,fragnatm(ifrag),frag(ifrag,:))
							write(*,"(' Atoms in fragment',i3,':')") ifrag
							write(*,"(15i5)") frag(ifrag,1:fragnatm(ifrag))
						end do
						close(10)
						exit
					else
						write(*,*) "Unable to find the file, input again"
					end if
				end do
			else
				allocate(frag(nfrag,ncenter),fragnatm(nfrag))
				do ifrag=1,nfrag
					write(*,"(' Input atomic indices for fragment',i4,', e.g. 1,4,8-12,15')") ifrag
					read(*,"(a)") c2000tmp
					call str2arr(c2000tmp,fragnatm(ifrag),frag(ifrag,:))
				end do
			end if
		end if
	
		allocate(he_frag(nfrag,3))
		he_frag=0
		write(*,*) "Contribution of each fragment to hole and electron:"
		do ifrag=1,nfrag
			do idx=1,fragnatm(ifrag)
				iatm=frag(ifrag,idx)
				he_frag(ifrag,1:2)=he_frag(ifrag,1:2)+he_atm(iatm,1:2)
			end do
! 			he_frag(ifrag,3)=minval(he_frag(ifrag,1:2)) !Not as good as below
			if (he_frag(ifrag,1)>0.and.he_frag(ifrag,2)>0) he_frag(ifrag,3)=dsqrt(he_frag(ifrag,1)*he_frag(ifrag,2))
			write(*,"(' #',i3,'   Hole:',f6.2,' %  Electron:',f6.2,' %  Overlap:',f6.2,' %  Diff.:',f7.2,' %')") &
			ifrag,he_frag(ifrag,:)*100,(he_frag(ifrag,2)-he_frag(ifrag,1))*100
		end do
		clrlimlow=0
		clrlimhigh=maxval(he_frag)
		stepz=(clrlimhigh-clrlimlow)/5D0
		xyratio=0.4D0
		nstepsize=1
		write(*,*) "Constructing fragment contribution to hole/electron finished!"
		
	else if (isel==1.or.isel==3) then
		if (isel==3) isavepic=1
		lengthx=2000
		lengthy=nint(xyratio*lengthx)
		if (allocated(frag)) then
			numx=nfrag
		else
			if (ifhydrogen==0) numx=ncennoh
			if (ifhydrogen==1) numx=ncenter
		end if
		numy=3
		call SCRMOD('REVERSE')
		CALL setxid(0,'NONE')
		CALL PAGE(3200,2700)
		if (isavepic==0) then
			call METAFL('xwin')
			call window(200,100,900,770)
		else if (isavepic==1) then
			call METAFL("pdf")
			call winsiz(graph2Dwidth,graph2Dheight) !Because color bar need to draw, so width is bigger than height
			CALL IMGFMT('RGB')
		end if
		call DISINI
		call ERRMOD("ALL","OFF") !If do not set this, when atom label in contour map is out of page range, DISLIN annoys users
		itexth=45
        call HWFONT
		if (isavepic==0.and.allocated(frag)) itexth=55 !Since the number of fragments is small, use larger labels
		call height(itexth)
		CALL LABDIG(-1,"X")
		CALL LABDIG(-1,"Y")
		CALL LABDIG(3,"Z")
		call ticks(1,"XYZ")
		call WINTIT("Composition of hole and electron")
		ielelabx=70-(itexth-45)
		ielelaby=2700/2-itexth-20
		call messag("Hole",ielelabx,ielelaby-lengthy/3)
		call messag("Electron",ielelabx,ielelaby)
		call messag("Overlap",ielelabx,ielelaby+lengthy/3)
		call center
		call AUTRES(numx,numy)
		call AX3LEN(lengthx,lengthy,lengthy) !lengthx is 2300, lengthy&z is 300
		call sursze(1D0,dfloat(numx),1D0,dfloat(numy)) !Manually set center position of starting and ending grids to ensure boundary grids have the same size as internal grids
		call labels("NONE","Y")
		call GRAF3(0.5D0,numx+0.5D0,1D0,dfloat(nstepsize), numy+0.5D0,0.5D0,dfloat(numy),-1D0, clrlimlow,clrlimhigh,clrlimlow,stepz)
		if (allocated(frag)) then
			call CRVMAT(he_frag(1:numx,1:numy),numx,numy,1,1)
		else
			if (ifhydrogen==0) call CRVMAT(he_atmnoh(1:numx,1:numy),numx,numy,1,1)
			if (ifhydrogen==1) call CRVMAT(he_atm(1:numx,1:numy),numx,numy,1,1)
		end if
		call DISFIN
		if (isel==3) then
			isavepic=0
			write(*,"(a)") " Done, the image has been saved as pdf format in current folder with ""dislin"" prefix"
		end if
	else if (isel==2) then
		if (allocated(frag)) then
			open(10,file="he_frag.txt",status="replace")
			do ifrag=1,nfrag
				write(10,"(i5,3f10.3)") ifrag,he_frag(ifrag,1:3)*100
			end do
			write(*,"(a)") " he_frag.txt has been exported to current folder! The column 1,2,3,4 correspond to &
			&fragment index, percentage contribution to hole, electron and overlap, respectively"
			close(10)
		else
			if (ifhydrogen==1) then
				open(10,file="he_atm.txt",status="replace")
				do iatm=1,ncenter
					write(10,"(i5,3f10.3)") iatm,he_atm(iatm,1:3)*100
				end do
				write(*,"(a)") " he_atm.txt has been exported to current folder! The column 1,2,3,4 correspond to &
				&atom index, percentage contribution to hole, electron and overlap, respectively"
				close(10)
			else if (ifhydrogen==0) then
				open(10,file="he_atmnoh.txt",status="replace")
				do idx=1,ncennoh
					write(10,"(i5,3f10.3)") nohlist(idx),he_atmnoh(idx,1:3)*100
				end do
				write(*,"(a)") " he_atmnoh.txt has been exported to current folder! The column 1,2,3,4 correspond to &
				&atom index, percentage contribution to hole, electron and overlap, respectively"
				close(10)
			end if
		end if
	else if (isel==4) then
		write(*,*) "Please input a number, e.g. 5"
		read(*,*) nstepsize
	else if (isel==5) then
		write(*,*) "Input lower and upper limits, e.g. 0,0.5"
		read(*,*) clrlimlow,clrlimhigh
		stepz=(clrlimhigh-clrlimlow)/5D0
	else if (isel==6) then
		write(*,*) "Input the ratio"
		write(*,*) "e.g. 0.2 means length of Y axis is 0.2 times the length of X axis"
		read(*,*) xyratio
	else if (isel==7) then
		if (ifhydrogen==0) then
			ifhydrogen=1
			clrlimlow=0
			clrlimhigh=maxval(he_atm)
			write(*,*) "Contribution of each atom to hole and electron:"
			do iatm=1,ncenter
				write(*,"(i5,'(',a,')  Hole:',f6.2,' %  Electron:',f6.2,' %  Overlap:',f6.2,' %  Diff.:',f7.2,' %')") &
				iatm,a(iatm)%name,he_atm(iatm,:)*100,(he_atm(iatm,2)-he_atm(iatm,1))*100
			end do
		else if (ifhydrogen==1) then
			ifhydrogen=0
			clrlimlow=0
			clrlimhigh=maxval(he_atmnoh)
		end if
	end if
end do
end subroutine


!!------- Calculate atomic contribution to hole and electron via Mulliken-like partition
!ioutinfo=1 means output prompts, 0 means silent
subroutine atmcontri_holeele_Mulliken(atmhole,atmele,ioutinfo)
use defvar
use util
use excitinfo
implicit real*8 (a-h,o-z)
real*8 atmhole(ncenter),atmele(ncenter)
real*8 Tmat(nbasis,nbasis) !Used to store intermediate data during calculating atomic contribution

call ask_Sbas_PBC !For PBC case, calculate Sbas if it is not currently available

atmhole=0;atmele=0
call walltime(iwalltime1)
if (ioutinfo==1) write(*,*) "Evaluating atomic contributions..."
thres=0.001D0 !Threshold of product of configuration coefficient for ignoring inter-configuration contributions
do iexcitorb=1,excnorb
	ileft=orbleft(iexcitorb)
	iright=orbright(iexcitorb)
	coeffi=exccoeff(iexcitorb)
	idir=excdir(iexcitorb)
	call showprog(iexcitorb,excnorb)
	do jexcitorb=1,excnorb
		jleft=orbleft(jexcitorb)
		jright=orbright(jexcitorb)
		coeffj=exccoeff(jexcitorb)
		jdir=excdir(jexcitorb)
		coeffprod=coeffi*coeffj
		if (idir/=jdir) cycle
		if (abs(coeffprod)<thres) cycle !Ignore too small terms
		
		!Contribution to hole
		if (iright==jright) then
			if (ileft<=nbasis) then
				Tmat=matmul(CObasa(:,ileft:ileft),transpose(CObasa(:,jleft:jleft)))
			else if (ileft>nbasis) then
				Tmat=matmul(CObasb(:,ileft-nbasis:ileft-nbasis),transpose(CObasb(:,jleft-nbasis:jleft-nbasis)))
			end if
			Tmat=Tmat*Sbas
			!$OMP PARALLEL DO SHARED(atmhole) PRIVATE(iatm,ibeg,iend,tmphole) schedule(dynamic) NUM_THREADS(nthreads)
			do iatm=1,ncenter
				ibeg=basstart(iatm)
                if (ibeg==0) cycle
				iend=basend(iatm)
				tmphole=( sum(Tmat(ibeg:iend,:))+sum(Tmat(:,ibeg:iend)) )*coeffprod
				if (idir==1) then ! ->
					atmhole(iatm)=atmhole(iatm)+tmphole
				else ! <-
					atmhole(iatm)=atmhole(iatm)-tmphole
				end if
			end do
			!$OMP END PARALLEL DO
		end if
		
		!Contribution to electron
		if (ileft==jleft) then
			if (iright<=nbasis) then
				Tmat=matmul(CObasa(:,iright:iright),transpose(CObasa(:,jright:jright)))
			else if (iright>nbasis) then
				Tmat=matmul(CObasb(:,iright-nbasis:iright-nbasis),transpose(CObasb(:,jright-nbasis:jright-nbasis)))
			end if
			Tmat=Tmat*Sbas
			!$OMP PARALLEL DO SHARED(atmele) PRIVATE(iatm,ibeg,iend,tmpele) schedule(dynamic) NUM_THREADS(nthreads)
			do iatm=1,ncenter
				ibeg=basstart(iatm)
                if (ibeg==0) cycle
				iend=basend(iatm)
				tmpele=( sum(Tmat(ibeg:iend,:))+sum(Tmat(:,ibeg:iend)) )*coeffprod
				if (idir==1) then ! ->
					atmele(iatm)=atmele(iatm)+tmpele
				else ! <-
					atmele(iatm)=atmele(iatm)-tmpele
				end if
			end do
			!$OMP END PARALLEL DO
		end if
		
	end do
end do
atmhole=atmhole/2
atmele=atmele/2
if (wfntype==0.or.wfntype==3) then
	atmhole=atmhole*2
	atmele=atmele*2
end if
!Normalize the data, because minor terms have been ignored, the sum of all data is not exactly 100%
atmhole=atmhole/sum(atmhole)
atmele=atmele/sum(atmele)

call walltime(iwalltime2)
if (ioutinfo==1) write(*,"(' Calculation took up time',i10,' s')") iwalltime2-iwalltime1
end subroutine




!!------- Calculate atomic contribution to hole and electron via Hirshfeld partition
!ioutinfo=1 means output prompts, 0 means silent
subroutine atmcontri_holeele_Hirshfeld(atmhole,atmele,ioutinfo)
use defvar
use util
use functions
use excitinfo
implicit real*8 (a-h,o-z)
real*8 atmhole(ncenter),atmele(ncenter)
real*8 hole(radpot*sphpot),ele(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),selfdens(radpot*sphpot)
real*8 atmrho(ncenter),tvec(3),atmhole_tmp(ncenter),atmele_tmp(ncenter)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

atmhole=0
atmele=0
call walltime(iwalltime1)
write(*,*) "Evaluating atomic contributions..."

if (ifPBC==0) then !Atomic-center grids
	if (ioutinfo==1) write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot

	!Generate quadrature point and weights by combination of Gauss-Chebyshev and Lebedev grids
	call gen1cintgrid(gridatmorg,iradcut)

	do iatm=1,ncenter
		call showprog(iatm,ncenter)
		gridatm%value=gridatmorg%value !Weight in this grid point
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		!Calculate hole and electron at all atom grids
		!$OMP parallel do shared(hole,ele) private(i) num_threads(nthreads)
		do i=1+iradcut*sphpot,radpot*sphpot
			call calc_holeele(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,0.01D0,hole(i),ele(i))
		end do
		!$OMP end parallel do
		!Calculate free atomic density to obtain promolecule density
		promol=0D0
		do jatm=1,ncenter
			!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
			do ipt=1+iradcut*sphpot,radpot*sphpot
				tmpdens(ipt)=calcatmdens(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
			end do
			!$OMP end parallel do
			promol=promol+tmpdens
			if (jatm==iatm) selfdens=tmpdens
		end do
		do i=1+iradcut*sphpot,radpot*sphpot
			if (promol(i)/=0D0) then
				tmpv=selfdens(i)/promol(i)*gridatm(i)%value
				atmhole(iatm)=atmhole(iatm)+tmpv*hole(i)
				atmele(iatm)=atmele(iatm)+tmpv*ele(i)
			end if
		end do
	end do
    
else !PBC case, use evenly distributed integration grids
    write(*,"(' Grid spacing of',f6.3,' Bohr is used in integration')") 0.35D0
	call setgrid_for_PBC(0.35D0,2) !This is fully adequate for crude estimation
	call calc_dvol(dvol)
	ifinish=0;ishowprog=1
	ntmp=floor(ny*nz/100D0)
	!$OMP PARALLEL SHARED(atmhole,atmele,ifinish,ishowprog) PRIVATE(atmhole_tmp,atmele_tmp,i,j,k,tmpx,tmpy,tmpz,&
    !$OMP icell,jcell,kcell,tvec,iatm,dist2,atmrho,prorho,holeval,eleval) NUM_THREADS(nthreads)
	atmhole_tmp(:)=0
    atmele_tmp(:)=0
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
									atmrho(iatm)=atmrho(iatm)+eleraddens(a(iatm)%index,dsqrt(dist2),0)
								end if
							end do
						end do
					end do
				end do
				prorho=sum(atmrho(:))
                if (prorho>0) then
					call calc_holeele(tmpx,tmpy,tmpz,0.01D0,holeval,eleval)
                    atmhole_tmp(:)=atmhole_tmp(:)+atmrho(:)/prorho*holeval
                    atmele_tmp(:)=atmele_tmp(:)+atmrho(:)/prorho*eleval
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
    atmhole(:)=atmhole(:)+atmhole_tmp(:)*dvol
    atmele(:)=atmele(:)+atmele_tmp(:)*dvol
	!$OMP END CRITICAL
	!$OMP END PARALLEL
	if (ishowprog/=0) call showprog(100,100)
end if

!Normalize the data, because minor terms have been ignored, the sum of all data is not exactly 100%
if (ioutinfo==1) write(*,"(' Sum of hole:',f10.7, '   Sum of electron:',f10.7)") sum(atmhole),sum(atmele)
atmhole=atmhole/sum(atmhole)
atmele=atmele/sum(atmele)

if (ioutinfo==1) then
    call walltime(iwalltime2)    
    write(*,"(' Calculation took up time',i10,' s')") iwalltime2-iwalltime1
end if
end subroutine




!---------------------------------------------------------------------------------------------
!--------------- Calculate contribution of basis function to hole and electron ---------------
!---------------------------------------------------------------------------------------------
subroutine hole_ele_bascontri
use defvar
use util
use excitinfo
implicit real*8 (a-h,o-z)
real*8 bashole(nbasis),basele(nbasis),Tmat(nbasis,nbasis)

write(*,"(a)") " Input the printing criterion, e.g. 0.5 means only the basis functions having contribution >=0.5% will be printed"
read(*,*) printcrit

!if (ifPBC==3) then !Generate Sbas for PBC case
!	call ask_Sbas_PBC
!end if

call walltime(iwalltime1)
write(*,*) "Evaluating basis function contributions..."
thres=0.001D0 !Threshold of product of configuration coefficient for ignoring inter-configuration contributions
bashole=0
basele=0
do iexcitorb=1,excnorb
	ileft=orbleft(iexcitorb)
	iright=orbright(iexcitorb)
	coeffi=exccoeff(iexcitorb)
	idir=excdir(iexcitorb)
	call showprog(iexcitorb,excnorb)
	do jexcitorb=1,excnorb
		jleft=orbleft(jexcitorb)
		jright=orbright(jexcitorb)
		coeffj=exccoeff(jexcitorb)
		jdir=excdir(jexcitorb)
		coeffprod=coeffi*coeffj
		if (idir/=jdir) cycle
		if (abs(coeffprod)<thres) cycle !Ignore too small terms
		
		!Contribution to hole
		if (iright==jright) then
			if (ileft<=nbasis) then
				Tmat=matmul(CObasa(:,ileft:ileft),transpose(CObasa(:,jleft:jleft)))
			else if (ileft>nbasis) then
				Tmat=matmul(CObasb(:,ileft-nbasis:ileft-nbasis),transpose(CObasb(:,jleft-nbasis:jleft-nbasis)))
			end if
			Tmat=Tmat*Sbas
			do ibas=1,nbasis
				tmphole=( sum(Tmat(ibas,:))+sum(Tmat(:,ibas)) )*coeffprod
				if (idir==1) then ! ->
					bashole(ibas)=bashole(ibas)+tmphole
				else ! <-
					bashole(ibas)=bashole(ibas)-tmphole
				end if
			end do
		end if
		
		!Contribution to electron
		if (ileft==jleft) then
			if (iright<=nbasis) then
				Tmat=matmul(CObasa(:,iright:iright),transpose(CObasa(:,jright:jright)))
			else if (iright>nbasis) then
				Tmat=matmul(CObasb(:,iright-nbasis:iright-nbasis),transpose(CObasb(:,jright-nbasis:jright-nbasis)))
			end if
			Tmat=Tmat*Sbas
			do ibas=1,nbasis
				tmpele=( sum(Tmat(ibas,:))+sum(Tmat(:,ibas)) )*coeffprod
				if (idir==1) then ! ->
					basele(ibas)=basele(ibas)+tmpele
				else ! <-
					basele(ibas)=basele(ibas)-tmpele
				end if
			end do
		end if
		
	end do
end do
bashole=bashole/2
basele=basele/2
call walltime(iwalltime2)
write(*,"(' Calculation took up time',i10,' s')") iwalltime2-iwalltime1

if (wfntype==0.or.wfntype==3) then
	bashole=bashole*2
	basele=basele*2
end if
!Normalize the data, because minor terms have been ignored, the sum of all data is not exactly 100%
holenorm=sum(bashole)
elenorm=sum(basele)
write(*,"(/,' Sum of hole contributions:    ',f10.6)") holenorm
write(*,"(' Sum of electron contributions:',f10.6)") elenorm
bashole=bashole/holenorm
basele=basele/elenorm

write(*,"(/,a)") "  Basis  Type    Atom    Shell     Hole      Electron     Overlap      Diff."
accumhole=0
accumele=0
do ibas=1,nbasis
	if (abs(bashole(ibas))*100>printcrit.or.abs(basele(ibas))*100>printcrit) then
		ovlp=0
		if (bashole(ibas)>0.and.basele(ibas)>0) ovlp=dsqrt(bashole(ibas)*basele(ibas))*100
		write(*,"(i6,4x,a,i5,a,i5,f10.2,' %',f10.2,' %',f10.2,' %',f10.2,' %')") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),&
		'('//a(bascen(ibas))%name//')',basshell(ibas),bashole(ibas)*100,basele(ibas)*100,ovlp,basele(ibas)*100-bashole(ibas)*100
		accumhole=accumhole+bashole(ibas)*100
		accumele=accumele+basele(ibas)*100
	end if
end do
write(*,"(' Sum of above printed terms: ',f10.2,' %',f10.2,' %',12x,f10.2,' %')") accumhole,accumele,accumhole-accumele
end subroutine





!--------------------------------------------------------------------------------
!----------------- Check and modify configuration coefficients ------------------
!--------------------------------------------------------------------------------
subroutine modexccoeff
use excitinfo
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: exccoeffbackup(:) !Used to backup coefficient
character c80tmp*80,strdir*3,strspini,strspinj,strtmp1*10,strtmp2*10,c200tmp*200
integer,allocatable :: idxorder(:)

call loadallexcinfo(1)
if (nstates>1) then
	write(*,*) "Load configuration coefficients for which excited state? e.g. 2"
	read(*,*) istate
end if
call loadexccoeff(istate,1)
allocate(exccoeffbackup(excnorb))
exccoeffbackup=exccoeff !Backup information, so that they can be retrieved by option -1

!Print 10 MO transitions with highest contributions
allocate(idxorder(excnorb))
forall(i=1:excnorb) idxorder(i)=i
call sortr8(exccoeff,"abs",idxorder)
call invarri4(idxorder)
exccoeff=exccoeffbackup
write(*,*)
write(*,*) "Some MO transitions sorted by absolute contributions:"
do idx=1,min(excnorb,10)
	iexcitorb=idxorder(idx)
	imo=orbleft(iexcitorb)
	jmo=orbright(iexcitorb)
	strdir=" ->"
	contrisign=1
	if (excdir(iexcitorb)==2) then
		strdir=" <-"
		contrisign=-1
	end if
	if (wfntype==0.or.wfntype==3) then
		write(*,"(' #',i7,2x,i6,a,i6,'   Coeff.:',f10.5,'   Contri.:',f10.4,'%')") iexcitorb,imo,strdir,&
		jmo,exccoeff(iexcitorb),contrisign*exccoeff(iexcitorb)**2/0.5D0*100
	else
		strspini="A"
		if (imo>nbasis) then
			strspini="B"
			imo=imo-nbasis
		end if
		strspinj="A"
		if (jmo>nbasis) then
			strspinj="B"
			jmo=jmo-nbasis
		end if
		write(*,"(' #',i7,2x,i5,a,a,i5,a,'   Coeff.:',f10.5,'   Contri.:',f10.4,'%')") iexcitorb,imo,strspini,strdir,&
		jmo,strspinj,exccoeff(iexcitorb),contrisign*exccoeff(iexcitorb)**2*100
	end if
end do

do while(.true.)
	write(*,*)
	write(*,*) "    ------------ Check and modify configuration coefficients -------------"
	write(*,*) "-3 Export current excitation information to a plain text file"
	write(*,*) "-2 Print coefficient (and contribution to excitation) of some orbital pairs"
	write(*,*) "-1 Retrieve original coefficient of all orbital pairs"
	write(*,*) "0 Return"
	write(*,*) "1 Set coefficient of an orbital pair"
	write(*,*) "2 Set coefficient for specific range of orbital pairs"
	read(*,*) isel
	
	if (isel==-3) then
		write(*,*) "Input the file path to export information, e.g. C:\ltwd\state2.txt"
		read(*,"(a)") c200tmp
		open(10,file=c200tmp,status="replace")
		write(10,"(' Excited State 1',i3,f12.6)") excmulti,excene
		do iexcitorb=1,excnorb
			imo=orbleft(iexcitorb)
			jmo=orbright(iexcitorb)
			strdir=" ->"
			if (excdir(iexcitorb)==2) strdir=" <-"
			if (wfntype==0.or.wfntype==3) then
				write(10,"(i6,a,i6,f12.6)") imo,strdir,jmo,exccoeff(iexcitorb)
			else
				strspini="A"
				if (imo>nbasis) then
					strspini="B"
					imo=imo-nbasis
				end if
				strspinj="A"
				if (jmo>nbasis) then
					strspinj="B"
					jmo=jmo-nbasis
				end if
				write(10,"(i5,a,a,i5,a,f12.6)") imo,strspini,strdir,jmo,strspinj,exccoeff(iexcitorb)
			end if
		end do
		close(10)
		write(*,*) "Done! Exporting finished!"
		write(*,"(a)") " Note: This file can then be directly used to provide electronic excitation information for &
		&many functions, such as hole-electron analysis, NTO analysis and delta_r calculation"
		cycle
	else if (isel==-2) then
		write(*,*) "Input the threshold for printing, e.g. 0.01"
		write(*,"(a)") " Note: The orbital pairs whose absolute value of coefficient >= this value will be printed. &
		&If input 0, all orbital pairs will be printed"
		read(*,*) printthres
		ntmp=0
		do iexcitorb=1,excnorb
			if (abs(exccoeff(iexcitorb))<printthres) cycle
			imo=orbleft(iexcitorb)
			jmo=orbright(iexcitorb)
			strdir=" ->"
			contrisign=1
			if (excdir(iexcitorb)==2) then
				strdir=" <-"
				contrisign=-1
			end if
			if (wfntype==0.or.wfntype==3) then
				write(*,"(' #',i7,2x,i6,a,i6,'   Coeff.:',f10.5,'   Contri.:',f10.4,'%')") iexcitorb,imo,strdir,&
				jmo,exccoeff(iexcitorb),contrisign*exccoeff(iexcitorb)**2/0.5D0*100
			else
				strspini="A"
				if (imo>nbasis) then
					strspini="B"
					imo=imo-nbasis
				end if
				strspinj="A"
				if (jmo>nbasis) then
					strspinj="B"
					jmo=jmo-nbasis
				end if
				write(*,"(' #',i7,2x,i5,a,a,i5,a,'   Coeff.:',f10.5,'   Contri.:',f10.4,'%')") iexcitorb,imo,strspini,strdir,&
				jmo,strspinj,exccoeff(iexcitorb),contrisign*exccoeff(iexcitorb)**2*100
			end if
			ntmp=ntmp+1
		end do
		write(*,"(' There are',i7,' orbital pairs met the criterion you set')") ntmp
		cycle
	else if (isel==-1) then
		exccoeff=exccoeffbackup
		write(*,*) "Done, all coefficents have been reset to original values!"
	else if (isel==0) then
		return
	else if (isel==1) then !Modify a single coefficient
		if (wfntype==0.or.wfntype==3) then
			write(*,*) "Input index of the MOs orginally occupied and unoccupied, e.g. 12,26"
			read(*,*) iMO,jMO
		else
			write(*,*) "Input index of the MOs orginally occupied and unoccupied"
			write(*,*) "e.g. 12A,46A    or    23B,39B"
			read(*,*) strtmp1,strtmp2
			read(strtmp1(1:len_trim(strtmp1)-1),*) iMO
			read(strtmp2(1:len_trim(strtmp2)-1),*) jMO
			if (index(strtmp1,'B')/=0) then
				iMO=iMO+nbasis
				jMO=jMO+nbasis
			end if
		end if
		idir=1
		if (any(excdir==2)) then
			write(*,*) "Which is the type of the orbital transition?"
			if (isel==1) write(*,*) "1: Excitation (namely ->)  2: De-excitation (namely <-)"
			read(*,*) idir
		end if
		do iexcitorb=1,excnorb
			if (orbleft(iexcitorb)==iMO.and.orbright(iexcitorb)==jMO.and.excdir(iexcitorb)==idir) then
				ipair=iexcitorb
				exit
			end if
		end do
		if (iexcitorb==excnorb+1) then
			write(*,*) "Warning: Cannot find corresponding orbital pair!"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		write(*,"(' Original value is',f11.7,', now input the new value, e.g. 0.0143')") exccoeff(ipair)
		read(*,*) exccoeff(ipair)
		write(*,*) "Done, the coefficient has been modified"
	else if (isel==2) then !Set coefficient for specific range of orbital pairs
		!Select left part of orbital pair
		if (wfntype==0.or.wfntype==3) then
			write(*,*) "Input the index range of the MOs orginally occupied, e.g. 14,17"
			write(*,"(' Note: Range of occupied MOs is 1 to',i6)") nint(naelec)
		else
			write(*,*) "Input the index range of the MOs orginally occupied"
			write(*,*) "e.g. 12A,14A    or    22B,28B"
			write(*,"(' Note: Range of occupied alpha MOs is 1 to',i6)") nint(naelec)
			write(*,"('       Range of occupied Beta MOs is  1 to',i6)") nint(nbelec)
		end if
		write(*,*) "Note: Simply input 0 can select all orginally occupied MOs"
		read(*,"(a)") c80tmp
		if (c80tmp(1:1)=='0') then
			istart=1
			iend=nmo !This treatment is completely safe, do not worry
		else
			if (wfntype==0.or.wfntype==3) then
				read(c80tmp,*) istart,iend
			else
				read(c80tmp,*) strtmp1,strtmp2
				read(strtmp1(1:len_trim(strtmp1)-1),*) istart
				read(strtmp2(1:len_trim(strtmp2)-1),*) iend
				if (index(strtmp1,'B')/=0) then
					istart=istart+nbasis
					iend=iend+nbasis
				end if
			end if
		end if
		!Select right part of orbital pair
		if (wfntype==0.or.wfntype==3) then
			write(*,*) "Input the index range of the MOs orginally unoccupied, e.g. 72,93"
			write(*,"(' Note: Range of unoccupied MOs is',i6,' to',i6)") nint(naelec)+1,nbasis
		else
			write(*,*) "Input the index range of the MOs orginally unoccupied"
			write(*,*) "e.g. 72A,84A    or    62B,98B"
			write(*,"(' Note: Range of unoccupied alpha MOs is',i6,' to',i6)") nint(naelec)+1,nbasis
			write(*,"('       Range of unoccupied Beta MOs  is',i6,' to',i6)") nint(nbelec)+1,nbasis
		end if
		write(*,*) "Note: Simply input 0 can select all orginally unoccupied MOs"
		read(*,"(a)") c80tmp
		if (c80tmp(1:1)=='0') then
			jstart=1
			jend=nmo
		else
			if (wfntype==0.or.wfntype==3) then
				read(c80tmp,*) jstart,jend
			else
				read(c80tmp,*) strtmp1,strtmp2
				read(strtmp1(1:len_trim(strtmp1)-1),*) jstart
				read(strtmp2(1:len_trim(strtmp2)-1),*) jend
				if (index(strtmp1,'B')/=0) then
					jstart=jstart+nbasis
					jend=jend+nbasis
				end if
			end if
		end if
		idir=1
		if (any(excdir==2)) then
			write(*,*) "Which is the type of these orbital pairs?"
			write(*,*) "1: Excitation (namely ->)  2: De-excitation (namely <-)  3: Any"
			read(*,*) idir
		end if
		write(*,*) "Set the coefficients to which value?  e.g. 0.00148"
		read(*,*) tmpval
		iset=0
		do iexcitorb=1,excnorb
			if (orbleft(iexcitorb)>=istart.and.orbleft(iexcitorb)<=iend.and.orbright(iexcitorb)>=jstart.and.orbright(iexcitorb)<=jend) then
				if (idir==3.or.excdir(iexcitorb)==idir) then
					exccoeff(iexcitorb)=tmpval
					iset=iset+1
				end if
			end if
		end do
		write(*,"(' Coefficient of',i7,' orbital pairs have been set to',f12.7)") iset,tmpval
	end if
	
	!Coefficients have been modified, output current status
	sumsqrexc=0
	sumsqrdeexc=0
	do itmp=1,excnorb
		if (excdir(itmp)==1) sumsqrexc=sumsqrexc+exccoeff(itmp)**2
		if (excdir(itmp)==2) sumsqrdeexc=sumsqrdeexc-exccoeff(itmp)**2
	end do
	write(*,"(' Sum of current coeff.^2 of excitation and de-excitation:',2f11.6)") sumsqrexc,sumsqrdeexc
	nzerocoeff=count(exccoeff==0)
	if (nzerocoeff>0) write(*,"(' The number of pairs having zero coefficient:',i7)") nzerocoeff
end do
end subroutine







!===========================================================================
!------ Calculate delta_r index, see J. Chem. Theory Comput., 9, 3118 (2013)
!excitation and de-excitation cofficients are summed together, according to Eq.9 of the original paper
subroutine delta_r
use excitinfo
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: exccoefftot(:) !Store the coefficient combined from excitation and de-excitation of the same pair
real*8,allocatable :: orbcenx(:),orbceny(:),orbcenz(:) !Store centroid of MOs
real*8,allocatable :: GTFdipint(:,:)
integer,allocatable :: stateidx(:)
character strtmp1,strtmp2,selectyn,c2000tmp*2000

if (ifPBC/=0) then
	write(*,"(a)") " Error: This function does not support periodic case yet! Press ENTER button to return"
    read(*,*)
    return
end if

write(*,*) "Initializing data, please wait..."
write(*,*)
!Calculate dipole moment integral matrix
nsize=nprims*(nprims+1)/2
allocate(GTFdipint(3,nsize))
call genGTFDmat(GTFdipint,nsize)
!Calculate centroid of MOs
allocate(orbcenx(nmo),orbceny(nmo),orbcenz(nmo))
orbcenx=0
orbceny=0
orbcenz=0
!$OMP PARALLEL DO SHARED(orbcenx,orbceny,orbcenz) PRIVATE(iGTF,jGTF,ides,imo) schedule(dynamic) NUM_THREADS(nthreads)
do imo=1,nmo
	do iGTF=1,nprims
		do jGTF=1,nprims
			if (iGTF>=jGTF) then
				ides=iGTF*(iGTF-1)/2+jGTF
			else
				ides=jGTF*(jGTF-1)/2+iGTF
			end if
			orbcenx(imo)=orbcenx(imo)+CO(imo,iGTF)*CO(imo,jGTF)*GTFdipint(1,ides)
			orbceny(imo)=orbceny(imo)+CO(imo,iGTF)*CO(imo,jGTF)*GTFdipint(2,ides)
			orbcenz(imo)=orbcenz(imo)+CO(imo,iGTF)*CO(imo,jGTF)*GTFdipint(3,ides)
		end do
	end do
end do
!$OMP END PARALLEL DO
deallocate(GTFdipint)

call loadallexcinfo(1) !Print summary of available excited states

write(*,*) "Please input index of the excited states that you want to calculate delta_r"
write(*,*) "e.g. 1-6,10,15-20"
write(*,"(a)") " Note: If only one state is selected, the delta_r can be further decomposed as orbital pair contribution"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,nstatecalc)
allocate(stateidx(nstatecalc))
call str2arr(c2000tmp,nstatecalc,stateidx)
if (any(stateidx>nstates).or.any(stateidx<1)) then
	write(*,*) "Error: The index you inputted exceeded valid range!"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	return
end if

write(*,*) "Calculating, please wait..."
write(*,*)

do idx=1,nstatecalc
	istate=stateidx(idx)
	if (nstatecalc>1) then
		call loadexccoeff(istate,0) !Do not print summary of this state
	else
		call loadexccoeff(istate,1)
	end if

	!Combine coefficient square of excitations and de-excitations
	allocate(exccoefftot(excnorb))
	exccoefftot=exccoeff
	coeffsumsqr=0D0
	do itmp=1,excnorb
		if (excdir(itmp)==1) then !->, find corresponding <- pair and sum its coefficient to here
			do jtmp=1,excnorb
				if (excdir(jtmp)==2.and.orbleft(itmp)==orbleft(jtmp).and.orbright(itmp)==orbright(jtmp)) exccoefftot(itmp)=exccoefftot(itmp)+exccoeff(jtmp)
			end do
		end if
		coeffsumsqr=coeffsumsqr+exccoefftot(itmp)**2
	end do

	delta_r_value=0
	do itmp=1,excnorb
		if (excdir(itmp)==2) cycle
		imo=orbleft(itmp)
		jmo=orbright(itmp)
		delta_r_value=delta_r_value+exccoefftot(itmp)**2/coeffsumsqr *dsqrt((orbcenx(imo)-orbcenx(jmo))**2+(orbceny(imo)-orbceny(jmo))**2+(orbcenz(imo)-orbcenz(jmo))**2)
	end do
	write(*,"(' Excited state',i5,':   Delta_r =',f12.6,' Bohr,',f12.6,' Angstrom')") istate,delta_r_value,delta_r_value*b2a
	
	!Only selected on excited state, decompose result
	if (nstatecalc==1) then
		write(*,*)
		write(*,*) "If printing orbital pair contributions to delta_r? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			write(*,*) "Input threshold for printing contribution (Angstrom), e.g. 0.05"
			write(*,"(a)") " Note: If input -1, then all contributions will be exported to delta_r.txt in current folder"
			read(*,*) printthres
			iout=6
			if (printthres==-1) then
				iout=10
				open(10,file="delta_r.txt",status="replace")
			end if
			write(iout,"(a)") " Note: The configuration coefficients shown below have combined both excitation and de-excitation parts"
			write(iout,"(' Sum of square of configuration coefficients:',f12.6)") coeffsumsqr
			write(iout,*) "   #Pair     Orbitals      Coefficient     Contribution (Bohr and Angstrom)"
			do itmp=1,excnorb
				if (excdir(itmp)==2) cycle
				imo=orbleft(itmp)
				jmo=orbright(itmp)
				contrival=exccoefftot(itmp)**2/coeffsumsqr *dsqrt((orbcenx(imo)-orbcenx(jmo))**2+(orbceny(imo)-orbceny(jmo))**2+(orbcenz(imo)-orbcenz(jmo))**2)
				if (contrival*b2a<printthres) cycle
				if (wfntype==0.or.wfntype==3) then
					write(iout,"(i8,2i7,f16.7,3x,2f16.7)") itmp,imo,jmo,exccoefftot(itmp),contrival,contrival*b2a
				else
					strtmp1="A"
					strtmp2="A"
					if (imo>nbasis) then
						imo=imo-nbasis
						strtmp1="B"
					end if
					if (jmo>nbasis) then
						jmo=jmo-nbasis
						strtmp2="B"
					end if
					write(iout,"(i8,i6,a,i6,a,f16.7,3x,2f16.7)") itmp,imo,strtmp1,jmo,strtmp2,exccoefftot(itmp),contrival,contrival*b2a
				end if
			end do
			if (printthres==-1) then
				close(10)
				write(*,*) "Done, outputting finished"
			end if
		end if
	end if
	deallocate(exccoefftot)
end do
end subroutine





!===========================================================================================
!------ Calculate lambda index to study CT character, see J. Chem. Phys., 128, 044118 (2008)
!A lot of codes inherit from delta_r routine
subroutine lambda_excit
use excitinfo
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8,allocatable :: exccoefftot(:) !Store the coefficient combined from excitation and de-excitation of the same pair
integer,allocatable :: stateidx(:)
real*8 orbovlp(nmo,nmo) !Only (occ,vir) blocks are to be calculated
real*8 beckeweigrid(radpot*sphpot),orbval(nmo,radpot*sphpot),orbvalpt(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
character strtmp1,strtmp2,selectyn,c2000tmp*2000
! excitfilename="x\excittrans\4-Nitroaniline\4-Nitroaniline.out"

if (ifPBC/=0) then
	write(*,"(a)") " Error: This function does not support periodic case yet! Press ENTER button to return"
    read(*,*)
    return
end if

call walltime(iwalltime1)
write(*,*) "Calculating overlap matrix between MOs, please wait patiently..."
if (iautointgrid==1) then !Allow change integration grid adapatively
	nradpotold=radpot
	nsphpotold=sphpot
	radcutold=radcut
	!Although the automatic settings is seemingly quite low, my lots of tests showed that the lose of accuracy was found to be marginal, &
	!Because MOs are smooth and easy to be integrated
	radpot=15
	sphpot=170
	radcut=18 !Enlarge radcut, because for Rydberg orbital, the default radcut 10 Bohr is not sufficient
	write(*,"(a)") " Note: The default number of integration grids in general should be sufficient. If you want to change, &
	&set ""iautointgrid"" in settings.ini to 0, and set ""radpot"" and ""sphpot"" to proper values"
end if
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)
orbovlp=0

do iatm=1,ncenter
	call showprog(iatm,ncenter)
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	!$OMP parallel do shared(orbval) private(ipt) num_threads(nthreads)
	do ipt=1+iradcut*sphpot,radpot*sphpot
		call orbderv(1,1,nmo,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,orbval(:,ipt))
	end do
	!$OMP end parallel do
	orbval=abs(orbval)
	
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	do ipt=1+iradcut*sphpot,radpot*sphpot
		orbvalpt=orbval(:,ipt)
		wei=gridatmorg(ipt)%value*beckeweigrid(ipt)
		!Total or alpha
		!$OMP parallel do shared(orbovlp) private(imo,jmo) num_threads(nthreads)
		do imo=1,nint(naelec)
			do jmo=nint(naelec)+1,nbasis
				orbovlp(imo,jmo)=orbovlp(imo,jmo)+orbvalpt(imo)*orbvalpt(jmo)*wei
			end do
		end do
		!$OMP end parallel do
		!Beta
		if (wfntype==1) then
			!$OMP parallel do shared(orbovlp) private(imo,jmo) num_threads(nthreads)
			do imo=nbasis+1,nbasis+nint(nbelec)
				do jmo=nbasis+nint(nbelec)+1,nmo
					orbovlp(imo,jmo)=orbovlp(imo,jmo)+orbvalpt(imo)*orbvalpt(jmo)*wei
				end do
			end do
			!$OMP end parallel do
		end if
	end do
end do
call walltime(iwalltime2)
write(*,"(' Preparation finished! Calculation took up time',i10,' s')") iwalltime2-iwalltime1
write(*,*)

do while(.true.)
	call loadallexcinfo(1) !Print summary of available excited states

	write(*,*) "Please input index of the excited states that you want to calculate lambda"
	write(*,*) "e.g. 1-6,10,15-20"
	write(*,"(a)") " Note: If only one state is selected, the lambda can be further decomposed as orbital pair contribution"
	read(*,"(a)") c2000tmp
	call str2arr(c2000tmp,nstatecalc)
	if (allocated(stateidx)) deallocate(stateidx)
	allocate(stateidx(nstatecalc))
	call str2arr(c2000tmp,nstatecalc,stateidx)
	if (any(stateidx>nstates).or.any(stateidx<1)) then
		write(*,*) "Error: The index you inputted exceeded valid range!"
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		return
	end if

	write(*,*) "Calculating, please wait..."
	write(*,*)
	do idx=1,nstatecalc
		istate=stateidx(idx)
		if (nstatecalc>1) then
			call loadexccoeff(istate,0) !Do not print summary of this state
		else
			call loadexccoeff(istate,1)
		end if

		!Combine coefficient square of excitations and de-excitations
		allocate(exccoefftot(excnorb))
		exccoefftot=exccoeff
		coeffsumsqr=0D0
		do itmp=1,excnorb
			if (excdir(itmp)==1) then !->, find corresponding <- pair and sum its coefficient to here
				do jtmp=1,excnorb
					if (excdir(jtmp)==2.and.orbleft(itmp)==orbleft(jtmp).and.orbright(itmp)==orbright(jtmp)) exccoefftot(itmp)=exccoefftot(itmp)+exccoeff(jtmp)
				end do
			end if
			coeffsumsqr=coeffsumsqr+exccoefftot(itmp)**2
		end do

		rlambda=0
		do itmp=1,excnorb
			if (excdir(itmp)==2) cycle
			imo=orbleft(itmp)
			jmo=orbright(itmp)
			rlambda=rlambda+exccoefftot(itmp)**2/coeffsumsqr *orbovlp(imo,jmo)
		end do
		write(*,"(' Excited state',i5,':   lambda =',f12.6)") istate,rlambda
		
		!Only selected on excited state, decompose result
		if (nstatecalc==1) then
			write(*,*)
			write(*,*) "If printing orbital pair contributions to lambda? (y/n)"
			read(*,*) selectyn
			if (selectyn=='y'.or.selectyn=='Y') then
				write(*,*) "Input threshold for printing contribution"
				write(*,*) "e.g. 0.02 means MO pairs with contribution higher than 0.02 will be printed"
				write(*,"(a)") " Note: If input -1, then all contributions will be exported to lambda.txt in current folder"
				read(*,*) printthres
				iout=6
				if (printthres==-1) then
					iout=10
					open(10,file="lambda.txt",status="replace")
				end if
				write(iout,"(a)") " Note: The configuration coefficients shown below have combined both excitation and de-excitation parts"
				write(iout,"(' Sum of square of configuration coefficients:',f12.6)") coeffsumsqr
				write(iout,*) "   #Pair     Orbitals      Coefficient     Contribution"
				do itmp=1,excnorb
					if (excdir(itmp)==2) cycle
					imo=orbleft(itmp)
					jmo=orbright(itmp)
					contrival=exccoefftot(itmp)**2/coeffsumsqr *orbovlp(imo,jmo)
					if (contrival<printthres) cycle
					if (wfntype==0.or.wfntype==3) then
						write(iout,"(i8,2i7,f16.7,3x,f16.7)") itmp,imo,jmo,exccoefftot(itmp),contrival
					else
						strtmp1="A"
						strtmp2="A"
						if (imo>nbasis) then
							imo=imo-nbasis
							strtmp1="B"
						end if
						if (jmo>nbasis) then
							jmo=jmo-nbasis
							strtmp2="B"
						end if
						write(iout,"(i8,i6,a,i6,a,f16.7,3x,f16.7)") itmp,imo,strtmp1,jmo,strtmp2,exccoefftot(itmp),contrival
					end if
				end do
				if (printthres==-1) then
					close(10)
					write(*,*) "Done, outputting finished"
				end if
			end if
		end if
	
		deallocate(exccoefftot)	
	end do
	
	write(*,*)
	write(*,*) "Do the lambda analysis again? (y/n)"
	read(*,*) selectyn
	if (selectyn=='n'.or.selectyn=='N') then
		radpot=nradpotold
		sphpot=nsphpotold
		radcut=radcutold
		exit
	end if
end do
end subroutine






!===========================================================================
!---- Generate NTOs, original paper of NTO: J. Chem. Phys., 118, 4775 (2003)
!The NTO eigenvalues for unrestricted wavefunction is 1/2 of the ones outputted by Gaussian, &
!which is incorrect (i.e. the sum is 2.0 rather than 1.0), the result must be that Gaussian used incorrect factor
subroutine NTO
use defvar
use util
use excitinfo
implicit real*8 (a-h,o-z)
real*8,allocatable :: T_MO(:,:),TT(:,:),NTOvec(:,:),NTOval(:)
real*8 tmparr(nbasis)
character c200tmp*200

if (.not.allocated(CObasa)) then
	write(*,*) "Error: The input file does not contain basis function information!"
	write(*,*) "Please read manual to make clear which kinds of input files could be used!"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	return
end if

call loadallexcinfo(1)
call selexcit(istate)
call loadexccoeff(istate,1)

NTOvalcoeff=2
if (allocated(CObasb)) NTOvalcoeff=1

!*** Alpha part or restricted case
nocc=nint(naelec)
nvir=nbasis-nocc
allocate(T_MO(nocc,nvir))
T_MO=0
do iexcitorb=1,excnorb
	iocc=orbleft(iexcitorb)
	ivir=orbright(iexcitorb)-nocc
	if (iocc>nbasis) cycle !Here we only process alpha part, index of Beta orbitals are higher than nbasis
	!Ignoring de-excitation, this treatment makes the result identical to that ouputted by Gaussian
	if (excdir(iexcitorb)==1) T_MO(iocc,ivir)=exccoeff(iexcitorb)
end do
!Occupied part
allocate(TT(nocc,nocc),NTOvec(nocc,nocc),NTOval(nocc))
TT=matmul(T_MO,transpose(T_MO))
call diagsymat(TT,NTOvec,NTOval,ierror)
NTOval=NTOval*NTOvalcoeff
MOene(1:nocc)=NTOval !By default, the diagsymat gives result from low to high
CObasa(:,1:nocc)=matmul(CObasa(:,1:nocc),NTOvec)
deallocate(TT,NTOvec,NTOval)
!Virtual part
allocate(TT(nvir,nvir),NTOvec(nvir,nvir),NTOval(nvir))
TT=matmul(transpose(T_MO),T_MO)
call diagsymat(TT,NTOvec,NTOval,ierror)
NTOval=NTOval*NTOvalcoeff
MOene(nocc+1:nbasis)=NTOval
CObasa(:,nocc+1:nbasis)=matmul(CObasa(:,nocc+1:nbasis),NTOvec)
do itmp=1,int(nvir/2D0) !Exchange array, so that the sequence will be high->low rather than the default low->high
	i=nocc+itmp
	j=nbasis+1-itmp
	tmpval=MOene(i)
	MOene(i)=MOene(j)
	MOene(j)=tmpval
	tmparr=CObasa(:,i)
	CObasa(:,i)=CObasa(:,j)
	CObasa(:,j)=tmparr
end do

write(*,*)
if (nvir>10) then
	if (allocated(CObasb)) then !Open-shell case
		write(*,*) "The highest 10 eigenvalues of alpha NTO pairs:"
    else !Closed-shell case
		write(*,*) "The highest 10 eigenvalues of NTO pairs:"
    end if
	write(*,"(5f12.6)") MOene(nocc+1:nocc+10)
    sum1=sum(MOene(nocc+1:nocc+10))
else
	if (allocated(CObasb)) then !Open-shell case
		write(*,*) "Eigenvalues of alpha NTO pairs:"
    else !Closed-shell case
		write(*,*) "Eigenvalues of NTO pairs:"
    end if
	write(*,"(5f12.6)") MOene(nocc+1:nbasis)
    sum1=sum(MOene(nocc+1:nbasis))
end if
if (allocated(CObasb)) then
	write(*,"(' Sum of all alpha eigenvalues:',f10.6)") sum1
else
	write(*,"(' Sum of all eigenvalues:',f10.6)") sum1
end if
deallocate(TT,NTOvec,NTOval,T_MO)

!*** Beta part
if (allocated(CObasb)) then
	nocc=nint(nbelec)
	nvir=nbasis-nocc
	allocate(T_MO(nocc,nvir))
	T_MO=0
	do iexcitorb=1,excnorb
		iocc=orbleft(iexcitorb)
		ivir=orbright(iexcitorb)-nocc
		if (iocc>nbasis) then !Only process transition of Beta orbitals
			iocc=iocc-nbasis
			ivir=ivir-nbasis
			if (excdir(iexcitorb)==1) T_MO(iocc,ivir)=exccoeff(iexcitorb)
		end if
	end do
	!Occupied part
	allocate(TT(nocc,nocc),NTOvec(nocc,nocc),NTOval(nocc))
	TT=matmul(T_MO,transpose(T_MO))
	call diagsymat(TT,NTOvec,NTOval,ierror)
	NTOval=NTOval*NTOvalcoeff
	MOene(nbasis+1:nbasis+nocc)=NTOval !By default, the diagsymat gives result from low to high
	CObasb(:,1:nocc)=matmul(CObasb(:,1:nocc),NTOvec)
	deallocate(TT,NTOvec,NTOval)
	!Virtual part
	allocate(TT(nvir,nvir),NTOvec(nvir,nvir),NTOval(nvir))
	TT=matmul(transpose(T_MO),T_MO)
	call diagsymat(TT,NTOvec,NTOval,ierror)
	NTOval=NTOval*NTOvalcoeff
	MOene(nbasis+nocc+1:nbasis+nbasis)=NTOval
	CObasb(:,nocc+1:nbasis)=matmul(CObasb(:,nocc+1:nbasis),NTOvec)
	do itmp=1,int(nvir/2D0) !Exchange array, so that the sequence will be high->low rather than the default low->high
		i=nocc+itmp
		j=nbasis+1-itmp
		tmpval=MOene(nbasis+i)
		MOene(nbasis+i)=MOene(nbasis+j)
		MOene(nbasis+j)=tmpval
		tmparr=CObasb(:,i)
		CObasb(:,i)=CObasb(:,j)
		CObasb(:,j)=tmparr
	end do
	write(*,*)
	if (nvir>10) then
		write(*,*) "The highest 10 eigenvalues of beta NTO pairs:"
		write(*,"(5f12.6)") MOene(nbasis+nocc+1:nbasis+nocc+10)
        sum2=sum(MOene(nbasis+nocc+1:nbasis+nocc+10))
	else
		write(*,*) "Eigenvalues of beta NTO pairs:"
		write(*,"(5f12.6)") MOene(nbasis+nocc+1:nbasis+nbasis)
        sum2=sum(MOene(nbasis+nocc+1:nbasis+nbasis))
	end if
	write(*,"(' Sum of all beta eigenvalues:',f10.6)") sum2
    write(*,*)
	write(*,"(' Sum of all alpha and beta NTO eigenvalues:',f10.6)") sum1+sum2
	deallocate(TT,NTOvec,NTOval,T_MO)
end if
write(*,*)
write(*,*) "0 Return"
write(*,*) "1 Output NTO orbitals to .molden file"
write(*,*) "2 Output NTO orbitals to .fch file"
write(*,*) "3 Output NTO orbitals to .mwfn file (recommended)"
read(*,*) iselNTO
if (iselNTO==1) then
	write(*,*) "Input the file path to output, e.g. C:\S1.molden"
	read(*,"(a)") c200tmp
	call outmolden(c200tmp,10)
	write(*,"(a)") " Exporting .molden file finished!"
else if (iselNTO==2) then
	write(*,*) "Input the file path to output, e.g. C:\S1.fch"
	read(*,"(a)") c200tmp
	call outfch(c200tmp,10,1)
	write(*,"(a)") " Exporting .fch file finished!"
else if (iselNTO==3) then
	write(*,*) "Input the file path to output, e.g. C:\S1.mwfn"
	read(*,"(a)") c200tmp
	call outmwfn(c200tmp,10,0)
	write(*,"(a)") " Exporting .mwfn file finished!"
end if
write(*,"(a)") " Now you can load the newly generated file to visualize NTOs. Note that in this file the orbital energies correspond to NTO eigenvalues"
write(*,*)
write(*,"(a)") " Reloading "//trim(filename)//" to recover initial status..."
call dealloall(0)
call readinfile(filename,1)
write(*,*) "Loading finished!"
end subroutine




!!---------- Calculate all transition electric dipole moments between all states and for each state
!Also, this function is able to produce electric dipole moment of each state
!Note that this function cannot borrow the code in hole_electron for loading data, because this function &
!need transition information of all states as well as excitation energies
subroutine exctransdip
use defvar
use util
use functions
use excitinfo
implicit real*8 (a-h,o-z)
real*8,allocatable :: GTFdipint(:,:) !Dipole moment integral between GTFs, use compressed index. The first index is x,y,z
real*8,allocatable :: MOdipint(:,:,:) !Dipole moment integral between all MOs. The first index is x,y,z
!The tdvecmat(:,i,j) will record transition electric dipole moment between i and j states.&
!i and j start from 0 and end at nstates. Only electronic contribution is taken into account
real*8,allocatable :: tdvecmat(:,:,:),Xnorm(:)
real*8 grounddip(3),statedip(3),nucdip(3) !Total ground state dipole moment, electronic contribution, nuclear contribution
real*8 tdvec(3),tmpvec(3)
character,allocatable :: allexclab(:)*5 !Label of each states
integer :: idiptype=1 !1=Electric transition dipole moment, 2= Magnetic, 0=Both
integer :: iGSESonly=0

call loadallexcinfo(0)
call loadallexccoeff(1)
allocate(tdvecmat(3,0:nstates,0:nstates),Xnorm(nstates))
tdvecmat=0

do while(.true.)
	write(*,*)
    if (iGSESonly==0) write(*,*) "-1: Toggle if only calculating between ground and excited states, current: No"
    if (iGSESonly==1) write(*,*) "-1: Toggle if only calculating between ground and excited states, current: Yes"
	!if (idiptype==0) write(*,*) "0 Choose type of (transition) dipole moment, current: Electric & Magnetic"
	if (idiptype==1) write(*,*) "0 Choose type of (transition) dipole moment, current: Electric"
	if (idiptype==2) write(*,*) "0 Choose type of (transition) dipole moment, current: Magnetic"
	write(*,*) "1 Output (transition) dipole moments to screen"
	write(*,*) "2 Output (transition) dipole moments to transdipmom.txt in current folder"
    if (idiptype==1) then
		write(*,*) "3 Generate input file of SOS module of Multiwfn as SOS.txt in current folder"
		write(*,"(a)") " 4 Output electric dipole moment (including both electron and nuclear contributions) of each excited state to dipmom.txt in current folder"
    end if
	read(*,*) isel
    if (isel==0) then
		write(*,*) "Calculate which type of (transition) dipole moment?"
		!write(*,*) "0 Both electric and magnetic (transition) dipole moments"
		write(*,*) "1 Electric"
		write(*,*) "2 Magnetic"
        read(*,*) idiptype
    else if (isel==-1) then
		if (iGSESonly==1) then
			iGSESonly=0
        else
			iGSESonly=1
        end if
    else
		exit
    end if
end do

write(*,*)

if (any(allexcmulti/=allexcmulti(1))) then
	if (isel==3) then
		write(*,"(a)") " Error: SOS.txt cannot be generated when singlet and triplet excited states are simultaneously occurred"
		write(*,*) "Press ENTER button to return"
		read(*,*)
		return
	else if (isel==4) then
		write(*,"(a)") " Error: This function is unavailable when singlet and triplet excited states are simultaneously occurred"
		write(*,*) "Press ENTER button to return"
		read(*,*)
		return
	end if
end if

!Calculate dipole moment integral matrix between all MOs
!If only GTF information is available, we calculate <MO|-r|MO> based on <GTF|-r|GTF>. While if basis function is available, &
!we calculate MO dipole moment integral matrix by unitary transformation, which is much faster
if (allocated(CObasa)) then
    write(*,"(a)") " Stage 1: Generating dipole moment integral matrix between basis functions..."
    if (idiptype==1.and..not.allocated(Dbas)) then
        call genDbas_curr
    else if (idiptype==2.and..not.allocated(Magbas)) then
        call genMagbas_curr
    else
        write(*,*) "This stage is skipped since the matrix is already available"
    end if
    
	call walltime(iwalltime1)
	write(*,*) "Stage 2: Calculating dipole moment integrals between all MOs..."
	allocate(MOdipint(3,nmo,nmo))
    if (idiptype==1) then !Electric
		allocate(DorbA(3,nbasis,nbasis))
		if (allocated(CObasb)) then
			allocate(DorbB(3,nbasis,nbasis))
			call genDorb
			MOdipint=0
			MOdipint(:,1:nbasis,1:nbasis)=DorbA
			MOdipint(:,nbasis+1:nmo,nbasis+1:nmo)=DorbB
			deallocate(DorbA,DorbB)
		else
			call genDorb
			MOdipint=DorbA
			deallocate(DorbA)
		end if
    else if (idiptype==2) then !Magnetic
		allocate(MagorbA(3,nbasis,nbasis))
		if (allocated(CObasb)) then
			allocate(MagorbB(3,nbasis,nbasis))
			call genMagorb
			MOdipint=0
			MOdipint(:,1:nbasis,1:nbasis)=MagorbA
			MOdipint(:,nbasis+1:nmo,nbasis+1:nmo)=MagorbB
			deallocate(MagorbA,MagorbB)
		else
			call genMagorb
			MOdipint=MagorbA
			deallocate(MagorbA)
		end if
    end if
else
	write(*,*) "Stage 1: Calculating dipole moment integrals between all GTFs..."
	nsize=nprims*(nprims+1)/2
	allocate(GTFdipint(3,nsize))
	if (idiptype==1) then
		call genGTFDmat(GTFdipint,nsize)
	else if (idiptype==2) then
		call genGTFMmat(GTFdipint,nsize)
    end if
	call walltime(iwalltime1)
	write(*,*) "Stage 2: Calculating dipole moment integrals between all MOs..."
	allocate(MOdipint(3,nmo,nmo))
	!MOdipint will record dipole moment integrals between all MOs, including all occ+vir alpha and occ+vir beta ones
	iprog=0
	!$OMP PARALLEL DO SHARED(MOdipint,iprog) PRIVATE(imo,jmo,iGTF,jGTF,ides,tmpvec) schedule(dynamic) NUM_THREADS(nthreads)
	do imo=1,nmo
		do jmo=imo,nmo
			tmpvec=0
			do iGTF=1,nprims
				do jGTF=1,nprims
					if (iGTF>=jGTF) then
						ides=iGTF*(iGTF-1)/2+jGTF
					else
						ides=jGTF*(jGTF-1)/2+iGTF
					end if
					tmpvec=tmpvec+CO(imo,iGTF)*CO(jmo,jGTF)*GTFdipint(:,ides)
				end do
			end do
			MOdipint(:,imo,jmo)=tmpvec
		end do
		if (nprims>300) then
			!$OMP CRITICAL
			iprog=iprog+1
			call showprog(iprog,nmo)
			!$OMP END CRITICAL
		end if
	end do
	!$OMP END PARALLEL DO
end if

!Fill lower triangle part
do imo=1,nmo
	do jmo=imo+1,nmo
		if (idiptype==1) then
			MOdipint(:,jmo,imo)=MOdipint(:,imo,jmo)
        else if (idiptype==2) then
			MOdipint(:,jmo,imo)=-MOdipint(:,imo,jmo)
        end if
	end do
end do
call walltime(iwalltime2)
write(*,"(' (Stage 2 took up wall clock time',i10,' s)')") iwalltime2-iwalltime1

if (isel==1) then
	iout=6
else if (isel==2) then
	iout=10
	open(iout,file="transdipmom.txt",status="replace")
else if (isel==3) then
	iout=10
	open(iout,file="SOS.txt",status="replace")
else if (isel==4) then
	iout=10
	open(iout,file="dipmom.txt",status="replace")
end if

fac=1
if (wfntype==0.or.wfntype==3) fac=2

!Calculate dipole moment of ground state
statedip=0
do imo=1,nmo
	statedip(:)=statedip(:)+MOocc(imo)*MOdipint(:,imo,imo)
end do
tdvecmat(:,0,0)=statedip(:)

!Calculate electric dipole moment of ground state
if (idiptype==1) then
	nucdip(1)=sum(a%charge*a%x)
	nucdip(2)=sum(a%charge*a%y)
	nucdip(3)=sum(a%charge*a%z)
	grounddip(:)=statedip(:)+nucdip(:)
end if

!Transition dipole moment between ground state and excited states
!Note that the formulae for electric and magnetic transition dipole moments are different, see Eqs. 22 and 24 in &
!"Ab initio calculations of oscillator and rotatory strengths in the random-phase approximation"
do iexc=1,nstates
	tdvec=0
	do ipair=1,allexcnorb(iexc)
		imo=allorbleft(ipair,iexc)
		lmo=allorbright(ipair,iexc)
		wei=allexccoeff(ipair,iexc)
        if (idiptype==1) then !Electronic
			tdvec(:)=tdvec(:)+wei*MOdipint(:,imo,lmo)
        else if (idiptype==2) then !Magnetic
			if (allexcdir(ipair,iexc)==1) then !Excitation
				tdvec(:)=tdvec(:)+wei*MOdipint(:,imo,lmo)
            else !De-excitation
				tdvec(:)=tdvec(:)-wei*MOdipint(:,imo,lmo)
            end if
        end if
	end do
	tdvecmat(:,0,iexc)=tdvec(:)*fac
end do

!Write SOS.txt. Output index, excitation energies and transition dipole moment between &
!ground state and excited states, the format is in line with SOS module
if (isel==3) then
	write(iout,*) nstates
	do i=1,nstates
		write(iout,"(i6,f12.6)") i,allexcene(i)
	end do
	do iexc=0,nstates
		write(iout,"(2i6,3(1PE15.6))") 0,iexc,tdvecmat(:,0,iexc)
	end do
end if

!write(*,*) statedip
!do imo=1,nmo
!write(*,"(i6,3f12.8)") imo,MOdipint(:,imo,imo)
!end do

if (isel<=3) then
	write(*,*) "Stage 3: Calculating transition dipole moment between excited states..."
else if (isel==4) then
	write(*,*) "Stage 3: Calculating dipole moment of all excited states..."
end if
!Below code works for both R and U reference states; for the latter, the alpha MO---beta MO case is naturally avoided, because QC program never gives such a orbital transition
!See comment in subroutine "genTDM_2exc" on how transition dipole moment between two TD excited states is calculated
call walltime(iwalltime1)
iprog=0
!$OMP PARALLEL DO SHARED(tdvecmat,iprog) PRIVATE(iexc,jexc,tdvec,ipair,jpair,imo,lmo,jmo,kmo,wei) schedule(dynamic) NUM_THREADS(nthreads)
do iexc=1,nstates
	do jexc=iexc,nstates
		if (iGSESonly==1) cycle
		if (allexcmulti(iexc)/=allexcmulti(jexc)) cycle !For 50-50 calculation, only calculate same-spin case
		if (isel==4.and.jexc/=iexc) cycle !Only calculate the moment between excited states with same index
		tdvec=0
		do ipair=1,allexcnorb(iexc)
			imo=allorbleft(ipair,iexc)
			lmo=allorbright(ipair,iexc)
			do jpair=1,allexcnorb(jexc)
				jmo=allorbleft(jpair,jexc)
				kmo=allorbright(jpair,jexc)
				wei=allexccoeff(ipair,iexc)*allexccoeff(jpair,jexc)
				if (allexcdir(ipair,iexc)==allexcdir(jpair,jexc)) then
					if (allexcdir(ipair,iexc)==2) wei=-wei
					if (imo==jmo.and.lmo/=kmo) then
						tdvec(:)=tdvec(:)+wei*MOdipint(:,lmo,kmo)
					else if (imo/=jmo.and.lmo==kmo) then
						tdvec(:)=tdvec(:)-wei*MOdipint(:,jmo,imo) !The order of j,i is correct here and important in the case of transition magnetic dipole moment
					else if (imo==jmo.and.lmo==kmo) then
						tdvec(:)=tdvec(:)+wei*(statedip(:)-MOdipint(:,imo,imo)+MOdipint(:,lmo,lmo))
					end if
				end if
			end do
		end do
		tdvecmat(:,iexc,jexc)=tdvec(:)*fac
	end do
	
	if (nprims>300) then
		!$OMP CRITICAL
		iprog=iprog+1
		call showprog(iprog,nstates)
		!$OMP END CRITICAL
	end if
end do
!$OMP END PARALLEL DO
call walltime(iwalltime2)
write(*,"(' (Stage 3 took up wall clock time',i10,' s)',/)") iwalltime2-iwalltime1
if (isel==1.or.isel==2) write(iout,"(a,/)") " Note: The transition dipole moments reported below only correspond to spatial part, the spin part is not taken into account!"

if (all(allexcmulti==allexcmulti(1))) then !All states have same spin, in this case all options are available
	if (idiptype==1) then !Electric
		if (isel==1.or.isel==2) then !Output transition dipole moment between various states
			!The ground state dipole moment shown below include both nuclear charge and electronic contributions
			write(iout,"(' Ground state electric dipole moment in X,Y,Z:',3f12.6,' a.u.',/)") grounddip
			write(iout,"(' Transition electric dipole moment between ground state (0) and excited states (a.u.)')")
			write(iout,*) "    i     j         X             Y             Z        Diff.(eV)   Oscil.str"
			do iexc=1,nstates
				oscillstr=2D0/3D0*allexcene(iexc)/au2eV*sum(tdvecmat(:,0,iexc)**2)
				write(iout,"(2i6,3f14.7,2f12.5)") 0,iexc,tdvecmat(:,0,iexc),allexcene(iexc),oscillstr
			end do
            if (iGSESonly==0) then
				write(iout,*)
				write(iout,"(' Note: In below output the case of i=j corresponds to contribution of electron to dipole moment of excited state i (contribution of nuclear charges is not taken into account!)')")
				write(iout,"(' Transition electric dipole moment between excited states (a.u.):')")
				write(iout,*) "    i     j         X             Y             Z        Diff.(eV)   Oscil.str"
            end if
		else if (isel==4) then !Output dipole moment for each state
			write(iout,"(a)") " Note: The electric dipole moments shown below include both nuclear charge and electronic contributions"
			write(iout,"(' Ground state electric dipole moment in X,Y,Z:',3f12.6,' a.u.',/)") grounddip
			write(iout,"(' Excited state electric dipole moments (a.u.):')")
			write(iout,*) " State         X             Y             Z        exc.(eV)    exc.(nm)"
		end if
		do iexc=1,nstates
			if (isel<=3.and.iGSESonly==0) then
				do jexc=iexc,nstates
					if (isel==1.or.isel==2) then
						enediff=abs(allexcene(jexc)-allexcene(iexc))
						oscillstr=2D0/3D0*enediff/au2eV*sum(tdvecmat(:,iexc,jexc)**2)
						write(iout,"(2i6,3f14.7,2f12.5)") iexc,jexc,tdvecmat(:,iexc,jexc),enediff,oscillstr
					else if (isel==3) then
						write(iout,"(2i6,3(1PE15.6))") iexc,jexc,tdvecmat(:,iexc,jexc)
					end if
				end do
			else if (isel==4) then
				write(iout,"(i6,3f14.6,f12.4,f12.2)") iexc,tdvecmat(:,iexc,iexc)+nucdip(:),allexcene(iexc),1239.842D0/allexcene(iexc)
			end if
		end do
	else if (idiptype==2) then !Magnetic
		write(iout,"(' Transition magnetic dipole moment between ground state (0) and excited states (a.u.)')")
		write(iout,*) "    i     j         X             Y             Z        Diff.(eV)"
		do iexc=1,nstates
			write(iout,"(2i6,3f14.7,2f12.5)") 0,iexc,tdvecmat(:,0,iexc),allexcene(iexc)
		end do
		write(iout,*)
		write(iout,"(' Transition magnetic dipole moment between excited states (a.u.):')")
		write(iout,*) "    i     j         X             Y             Z        Diff.(eV)"
		do iexc=1,nstates
			do jexc=iexc,nstates
				if (isel==1.or.isel==2) then
					enediff=abs(allexcene(jexc)-allexcene(iexc))
					write(iout,"(2i6,3f14.7,2f12.5)") iexc,jexc,tdvecmat(:,iexc,jexc),enediff
				end if
			end do
		end do
    end if
else !Not all states have the same spin, 50-50 with singlet ground state is assumed. isel=3 and isel=4 are not available in this case
	allocate(allexclab(0:2*nstates))
	iS=0
	iT=0
	allexclab(0)=" S  0"
	do iexc=1,nstates
		if (allexcmulti(iexc)==1) then
			iS=iS+1
			write(allexclab(iexc),"(' S',i3)") iS
		else if (allexcmulti(iexc)==3) then
			iT=iT+1
			write(allexclab(iexc),"(' T',i3)") iT
		end if
	end do
	write(iout,"(' Note: In below output the case of i=j corresponds to contribution of electron to dipole moment of excited state i (contribution of nuclear charges is not taken into account!)')")
    if (idiptype==1) then
		write(iout,"(' Transition electric dipole moment between singlet states (a.u.):')")
    else if (idiptype==2) then
		write(iout,"(' Transition magnetic dipole moment between singlet states (a.u.):')")
    end if
	write(iout,*) "  i        j          X             Y             Z         Diff.(eV)  Osc.str."
    !Output S0 state
    write(iout,"(a,' -- ',a,3f14.7,f12.4,f12.5)") allexclab(0),allexclab(0),statedip(:),0D0
    !Output between S0 and various singlet excited states
	do iexc=1,nstates
		if (allexcmulti(iexc)/=1) cycle
        oscillstr=2D0/3D0*allexcene(iexc)/au2eV*sum(tdvecmat(:,0,iexc)**2)
		write(iout,"(a,' -- ',a,3f14.7,f12.4,f12.5)") allexclab(0),allexclab(iexc),tdvecmat(:,0,iexc),allexcene(iexc),oscillstr
	end do
    if (iGSESonly==0) then
		!Output between various singlet excited states
		do iexc=1,nstates
			do jexc=iexc,nstates
				oscillstr=2D0/3D0*abs(allexcene(jexc)-allexcene(iexc))/au2eV*sum(tdvecmat(:,iexc,jexc)**2)
				if (allexcmulti(iexc)==1.and.allexcmulti(jexc)==1) &
				write(iout,"(a,' -- ',a,3f14.7,f12.4,f12.5)") allexclab(iexc),allexclab(jexc),tdvecmat(:,iexc,jexc),abs(allexcene(jexc)-allexcene(iexc)),oscillstr
			end do
		end do
		!Output between various triplet excited states
		if (idiptype==1) then
			write(iout,"(/,' Transition electric dipole moment between triplet states (a.u.):')")
		else if (idiptype==2) then
			write(iout,"(/,' Transition magnetic dipole moment between triplet states (a.u.):')")
		end if
		write(iout,*) "  i        j          X             Y             Z         Diff.(eV)  Osc.str."
		do iexc=1,nstates
			do jexc=iexc,nstates
				oscillstr=2D0/3D0*abs(allexcene(jexc)-allexcene(iexc))/au2eV*sum(tdvecmat(:,iexc,jexc)**2)
				if (allexcmulti(iexc)==3.and.allexcmulti(jexc)==3) &
				write(iout,"(a,' -- ',a,3f14.7,f12.4,f12.5)") allexclab(iexc),allexclab(jexc),tdvecmat(:,iexc,jexc),abs(allexcene(jexc)-allexcene(iexc)),oscillstr
			end do
		end do
		write(iout,"(/,' Excitation energies (eV):')")
		do iexc=1,nstates
			if (allexcmulti(iexc)==1) write(iout,"(a,f12.4)") allexclab(iexc),allexcene(iexc)
		end do
		do iexc=1,nstates
			if (allexcmulti(iexc)==3) write(iout,"(a,f12.4)") allexclab(iexc),allexcene(iexc)
		end do
    end if
end if

if (isel==2) then
	close(iout)
	write(*,*) "Done! The result has been outputted to transdipmom.txt in current folder"
else if (isel==3) then
	close(iout)
	write(*,*) "Done! The result has been outputted to SOS.txt in current folder"
else if (isel==4) then
	close(iout)
	write(*,*) "Done! The result has been outputted to dipmom.txt in current folder"
end if
end subroutine




!!!------------- Analyze charge transfer based on grid data of density difference, See J. Chem. Theory Comput., 7, 2498
!Some indices in the original paper have been modified by me
subroutine CTanalyze
use GUI
use defvar
implicit real*8 (a-h,o-z)
real*8,allocatable :: Cpos(:,:,:),Cneg(:,:,:),tmpmat(:,:,:)

if (.not.allocated(cubmat)) then
	write(*,"(a)") " Error: No grid data is presented, grid data of electron density difference must be &
    &firstly calculated by main function 5 or loaded from external file when Multiwfn boots up"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

sumpos=0D0
sumneg=0D0
Rxpos=0D0;Rypos=0D0;Rzpos=0D0
Rxneg=0D0;Ryneg=0D0;Rzneg=0D0
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			if (cubmat(i,j,k)>0) then
				sumpos=sumpos+cubmat(i,j,k)
				Rxpos=Rxpos+cubmat(i,j,k)*tmpx
				Rypos=Rypos+cubmat(i,j,k)*tmpy
				Rzpos=Rzpos+cubmat(i,j,k)*tmpz
			else
				sumneg=sumneg+cubmat(i,j,k)
				Rxneg=Rxneg+cubmat(i,j,k)*tmpx
				Ryneg=Ryneg+cubmat(i,j,k)*tmpy
				Rzneg=Rzneg+cubmat(i,j,k)*tmpz
			end if
		end do
	end do
end do
call calc_dvol(dvol)
sumpos=sumpos*dvol
sumneg=sumneg*dvol
Rxpos=Rxpos*dvol/sumpos
Rypos=Rypos*dvol/sumpos
Rzpos=Rzpos*dvol/sumpos
Rxneg=Rxneg*dvol/sumneg
Ryneg=Ryneg*dvol/sumneg
Rzneg=Rzneg*dvol/sumneg
disx=abs(Rxpos-Rxneg)
disy=abs(Rypos-Ryneg)
disz=abs(Rzpos-Rzneg)
disnorm=dsqrt(disx**2+disy**2+disz**2)
write(*,"(' q_CT (positive and negative parts):',2f8.3,' a.u.')") sumpos,sumneg
write(*,"(' Barycenter of positive part in x,y,z (Angstrom):',3f8.3)") Rxpos*b2a,Rypos*b2a,Rzpos*b2a
write(*,"(' Barycenter of negative part in x,y,z (Angstrom):',3f8.3)") Rxneg*b2a,Ryneg*b2a,Rzneg*b2a
write(*,"(' Distance of CT in x,y,z (Angstrom):',3f8.3,'  D index:',f8.3)") disx*b2a,disy*b2a,disz*b2a,disnorm*b2a
dipx=-(Rxpos-Rxneg)*sumpos
dipy=-(Rypos-Ryneg)*sumpos
dipz=-(Rzpos-Rzneg)*sumpos
dipnorm=disnorm*sumpos
write(*,"(' Dipole moment variation (a.u.) :',3f8.3,' Norm:',f8.3)") dipx,dipy,dipz,dipnorm
write(*,"(' Dipole moment variation (Debye):',3f8.3,' Norm:',f8.3)") dipx*au2debye,dipy*au2debye,dipz*au2debye,dipnorm*au2debye

sigxpos=0D0;sigypos=0D0;sigzpos=0D0
sigxneg=0D0;sigyneg=0D0;sigzneg=0D0
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			if (cubmat(i,j,k)>0) then
				sigxpos=sigxpos+cubmat(i,j,k)*(tmpx-Rxpos)**2
				sigypos=sigypos+cubmat(i,j,k)*(tmpy-Rypos)**2
				sigzpos=sigzpos+cubmat(i,j,k)*(tmpz-Rzpos)**2
			else
				sigxneg=sigxneg+cubmat(i,j,k)*(tmpx-Rxneg)**2
				sigyneg=sigyneg+cubmat(i,j,k)*(tmpy-Ryneg)**2
				sigzneg=sigzneg+cubmat(i,j,k)*(tmpz-Rzneg)**2
			end if
		end do
	end do
end do
sigxpos=dsqrt(sigxpos/(sumpos/dvol))
sigypos=dsqrt(sigypos/(sumpos/dvol))
sigzpos=dsqrt(sigzpos/(sumpos/dvol))
signormpos=dsqrt(sigxpos**2+sigypos**2+sigzpos**2)
sigxneg=dsqrt(sigxneg/(sumneg/dvol))
sigyneg=dsqrt(sigyneg/(sumneg/dvol))
sigzneg=dsqrt(sigzneg/(sumneg/dvol))
signormneg=dsqrt(sigxneg**2+sigyneg**2+sigzneg**2)
write(*,"(' RMSD of positive part in x,y,z (Angstrom):',3f7.3,' Total:',f8.3)") sigxpos*b2a,sigypos*b2a,sigzpos*b2a,signormpos*b2a
write(*,"(' RMSD of negative part in x,y,z (Angstrom):',3f7.3,' Total:',f8.3)") sigxneg*b2a,sigyneg*b2a,sigzneg*b2a,signormneg*b2a
diffsigx=sigxpos-sigxneg
diffsigy=sigypos-sigyneg
diffsigz=sigzpos-sigzneg
delta_sigCT=signormpos-signormneg
write(*,*) "Difference between RMSD of positive and negative parts (Angstrom):"
write(*,"(' X:',f8.3,'  Y:',f8.3,'  Z:',f8.3,'  delta_sigma index:',f8.3)") diffsigx*b2a,diffsigy*b2a,diffsigz*b2a,delta_sigCT*b2a
Hx=(sigxpos+sigxneg)/2D0
Hy=(sigypos+sigyneg)/2D0
Hz=(sigzpos+sigzneg)/2D0
H_CT= ( dsqrt ( (Hx*(Rxpos-Rxneg))**2 + (Hy*(Rypos-Ryneg))**2 + (Hz*(Rzpos-Rzneg))**2 ) ) / disnorm
H_index=(signormpos+signormneg)/2D0
write(*,"(' H_x:',f7.3,'  H_y:',f7.3,'  H_z:',f7.3,'  H_CT:',f7.3,'  H index:',f7.3,' Angstrom')") Hx*b2a,Hy*b2a,Hz*b2a,H_CT*b2a,H_index*b2a
t_index=disnorm-H_CT
write(*,"(' t index:',f8.3,' Angstrom')") t_index*b2a

allocate(Cpos(nx,ny,nz),Cneg(nx,ny,nz),tmpmat(nx,ny,nz))
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,rnowx,rnowy,rnowz)
			Cpos(i,j,k)=exp( -(rnowx-Rxpos)**2/(2*sigxpos**2) -(rnowy-Rypos)**2/(2*sigypos**2) -(rnowz-Rzpos)**2/(2*sigzpos**2))
			Cneg(i,j,k)=exp( -(rnowx-Rxneg)**2/(2*sigxneg**2) -(rnowy-Ryneg)**2/(2*sigyneg**2) -(rnowz-Rzneg)**2/(2*sigzneg**2))
		end do
	end do
end do
Cpos=Cpos*sumpos/(sum(Cpos)*dvol)
Cneg=Cneg*sumneg/(sum(Cneg)*dvol)
write(*,"(' Overlap integral between C+ and C- (i.e. S+- index):',f10.6)") sum(dsqrt(Cpos/sumpos)*dsqrt(Cneg/sumneg))*dvol

sur_value=0.001D0
do while(.true.)
	write(*,*)
	write(*,*) "0 Return"
	write(*,*) "1 Show isosurface of C+ and C- functions simultaneously"
	write(*,*) "2 Export C+ and C- functions to cube file in current folder"
	read(*,*) isel

	if (isel==0) then
		return
	else if (isel==1) then
		isosursec=1
		clrRcub1sameold=clrRcub1same !Backup previous color setting
		clrGcub1sameold=clrGcub1same
		clrBcub1sameold=clrBcub1same
		clrRcub2oppoold=clrRcub2oppo
		clrGcub2oppoold=clrGcub2oppo
		clrBcub2oppoold=clrBcub2oppo
		clrRcub2samemeshptold=clrRcub2samemeshpt
		clrGcub2samemeshptold=clrGcub2samemeshpt
		clrBcub2samemeshptold=clrBcub2samemeshpt
		clrRcub2oppomeshptold=clrRcub2oppomeshpt
		clrGcub2oppomeshptold=clrGcub2oppomeshpt
		clrBcub2oppomeshptold=clrBcub2oppomeshpt
		clrRcub1same=0.3D0 !Set color to that Cpos is green, Cneg is blue. (Cpos/Cneg function is positive/negative everywhere due to sumpos and sumneg)
		clrGcub1same=0.75D0
		clrBcub1same=0.3D0
		clrRcub2oppo=0.3D0
		clrGcub2oppo=0.45D0
		clrBcub2oppo=0.9D0
		clrRcub2samemeshpt=0.3D0
		clrGcub2samemeshpt=0.75D0
		clrBcub2samemeshpt=0.3D0
		clrRcub2oppomeshpt=0.3D0
		clrGcub2oppomeshpt=0.45D0
		clrBcub2oppomeshpt=0.9D0
		if (allocated(cubmattmp)) deallocate(cubmattmp)
		allocate(cubmattmp(nx,ny,nz))
		tmpmat=cubmat
		cubmat=Cpos
		cubmattmp=Cneg
		!Since drawisosurgui is mainly designed for plotting one isosurface, while here we use it to draw both cubmat and cubmattmp, so disable change of style to avoid troubles
		call drawisosurgui(2)
		cubmat=tmpmat !Recover original data, namely density difference between two states
		deallocate(cubmattmp)
		clrRcub1same=clrRcub1sameold !Recover previous color setting
		clrGcub1same=clrGcub1sameold
		clrBcub1same=clrBcub1sameold
		clrRcub2oppo=clrRcub2oppoold
		clrGcub2oppo=clrGcub2oppoold
		clrBcub2oppo=clrBcub2oppoold
		clrRcub2samemeshpt=clrRcub2samemeshptold
		clrGcub2samemeshpt=clrGcub2samemeshptold
		clrBcub2samemeshpt=clrBcub2samemeshptold
		clrRcub2oppomeshpt=clrRcub2oppomeshptold
		clrGcub2oppomeshpt=clrGcub2oppomeshptold
		clrBcub2oppomeshpt=clrBcub2oppomeshptold
	else if (isel==2) then
		open(10,file="Cpos.cub",status="replace")
		call outcube(Cpos,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,"(' C+ function has been outputted to ""Cpos.cub"" in current folder')")
		open(10,file="Cneg.cub",status="replace")
		call outcube(Cneg,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,"(' C- function has been outputted to ""Cneg.cub"" in current folder')")
	end if
end do
end subroutine





!-------------------------------------------------------------------------------------------
!---------- Plot atom/fragment transition matrix of various kinds as heat map -----------
!-------------------------------------------------------------------------------------------
subroutine TDMplot
use defvar
use util
use plot
use excitinfo
implicit real*8 (a-h,o-z)
character tdmatfilename*200,c200tmp*200,c2000tmp*2000
real*8 tdmatbas(nbasis,nbasis),tmatatmtmp(ncenter,ncenter),tmatatm(ncenter,ncenter)
real*8,allocatable :: tmatatmnoh(:,:),tmatfrag(:,:)
integer,allocatable :: frag(:,:),fragnatm(:)

!Load density transition matrix in basis representation
if (excitfilename==" ") then
	write(*,*) "Below kinds of files are acceptable"
	write(*,"(a)") " (1) Output file of electron excitation task of Gaussian, ORCA, CP2K, GAMESS-US/Firefly, CP2K, xtb, BDF, TDM will be calculated based on configuration and MO coefficients"
	write(*,*) "(2) Gaussian output file of electron excitation task with IOp(6/8=3)"
	write(*,*) "(3) tdmat.txt, which contains transition density matrix"
	write(*,"(a)") " (4) One of AAtrdip.txt, AAtrdipX.txt, AAtrdipY.txt, AAtrdipZ.txt, which contains atom transition dipole moment matrix"
	write(*,*) "(5) atmCTmat.txt, which contains atom-atom charge transfer matrix"
	write(*,*) "Please input the file path, e.g. C:\ltwd\excit.out. Input 0 can return" 
	do while(.true.)
		read(*,"(a)") tdmatfilename
		if (tdmatfilename=='0') return
		inquire(file=tdmatfilename,exist=alive)
		if (alive) exit
		write(*,*) "Error: Cannot find the file, please input again"
	end do
else
	tdmatfilename=excitfilename
end if

tdmatbas=0D0
open(10,file=tdmatfilename,status="old")
!Directly load and use atom transition dipole moment matrix (commonly generated by option 11 of main function 18)
if (index(tdmatfilename,"AAtrdip")/=0.or.index(tdmatfilename,"aatrdip")/=0) then
	write(*,*) "Loading atom transition dipole moment matrix..."
	call readmatgau(10,tmatatm,0,"f14.8",6,5,1)
!Directly load and use atom-atom charge transfer matrix (commonly generated by IFCT module of main function 18)
else if (index(tdmatfilename,"atmCTmat.txt")/=0.or.index(tdmatfilename,"atmctmat.txt")/=0) then
	write(*,*) "Loading atom-atom charge transfer matrix..."
	call readmatgau(10,tmatatm,0,"f14.8",6,5,1)
else if (index(tdmatfilename,"tdmat.txt")/=0) then !Directly load and use transition density matrix
	write(*,*) "Loading transition density matrix..."
	call readmatgau(10,tdmatbas,0,"f14.8",6,5,1)
else !Load or evaluate basis function based transition density matrix and then contract it to atom transition matrix
	call loclabel(10,"6/8=3",igauTDM,maxline=100)
	rewind(10)
	if (igauTDM==1) then !Gaussian output file with TDM printed by IOp(6/8=3), directly load the TDM
		write(*,*) "This is a Gaussian output file with transition density matrix"
		call loclabel(10,"Alpha Density Matrix:",ifound)
		if (ifound==1) then !Open-shell
			write(*,*) "Use which type of transition density matrix?"
			write(*,*) "1=Alpha    2=Beta"
			read(*,*) iTDMtype
			write(*,*) "Loading transition density matrix..."
			if (iTDMtype==1) then
				call readmatgau(10,tdmatbas,1,"f10.5",21,5,1)
			else if  (iTDMtype==2) then
				call loclabel(10,"Beta Density Matrix:",ifound)
				call readmatgau(10,tdmatbas,1,"f10.5",21,5,1)
			end if
		else !Closed-shell
			call loclabel(10,"Density Matrix:",ifound)
			if (ifound==0) call loclabel(10,"DENSITY MATRIX.",ifound)
			if (ifound==0) then
				write(*,"(a,/)") " Error: Cannot found transition density matrix information from the Gaussian output file"
				return
			end if
			write(*,*) "Loading transition density matrix..."
			call readmatgau(10,tdmatbas,1,"f10.5",21,5,1)
		end if
	else !Load electronic excitation information and then generate TDM
		write(*,*) "The input file does not contain transition density matrix, so generate it"
		excitfilename=tdmatfilename
		call loadallexcinfo(1)
		write(*,*) "Generate transition density matrix for which electron excitation? e.g. 2"
		read(*,*) istate
		call loadexccoeff(istate,1)
		call genTDM(1,1) !Generate tdmat, which is a matrix in module excitinfo
		if (allocated(tdmatb)) then !Unrestricted reference
			write(*,*)
			write(*,*) "Use which type of transition density matrix?"
			write(*,*) "0=Total (i.e. Alpha+Beta)   1=Alpha   2=Beta"
			read(*,*) iTDMsel
			if (iTDMsel==0) tdmatbas=tdmata+tdmatb
			if (iTDMsel==1) tdmatbas=tdmata
			if (iTDMsel==2) tdmatbas=tdmatb
		else !Closed-shell reference
			tdmatbas=tdmata
		end if
	end if
end if
close(10)

!Contract TDM to atom transition matrix
if (index(tdmatfilename,"AAtrdip")==0.and.index(tdmatfilename,"aatrdip")==0.and.index(tdmatfilename,"atmCTmat.txt")==0.and.index(tdmatfilename,"atmctmat.txt")==0) then
	write(*,"(a)") " Select the way of constructing atom transition matrix based on transition density matrix (TDM), see manual for detail"
	write(*,*) "1: Sum of square of TDM elements"
	write(*,*) "2: Square root of sum of square of TDM elements"
	write(*,*) "3: Sum of absolute value of TDM elements"
	write(*,*) "4: Proposed in Chem. Rev., 102, 3171 (2002)"
	read(*,*) iatmTMtype
	tmatatm=0D0
	if (iatmTMtype==1.or.iatmTMtype==2.or.iatmTMtype==3) then !Contraction by summing up of square of elements
		do iatm=1,ncenter
            if (basstart(iatm)==0) cycle
			do jatm=1,ncenter
				if (basstart(jatm)==0) cycle
				do ibas=basstart(iatm),basend(iatm)
					do jbas=basstart(jatm),basend(jatm)
						if (iatmTMtype==1.or.iatmTMtype==2) then
							tmatatm(iatm,jatm)=tmatatm(iatm,jatm)+tdmatbas(ibas,jbas)**2
						else
							tmatatm(iatm,jatm)=tmatatm(iatm,jatm)+abs(tdmatbas(ibas,jbas))
						end if
					end do
				end do
			end do
		end do
		if (iatmTMtype==2) tmatatm=dsqrt(tmatatm)
	else if (iatmTMtype==4) then !Eqs 2.20, 2.21 in Chem. Rev., 102, 3171 (2002)
		!Diagonal terms
		do iatm=1,ncenter
            if (basstart(iatm)==0) cycle
			do ibas=basstart(iatm),basend(iatm)
				tmatatm(iatm,iatm)=tmatatm(iatm,iatm)+abs(tdmatbas(ibas,ibas))
			end do
		end do
		!Off-diagonal terms
		do iatm=1,ncenter
            if (basstart(iatm)==0) cycle
			do jatm=1,ncenter
				if (basstart(jatm)==0) cycle
				do ibas=basstart(iatm),basend(iatm)
					do jbas=basstart(jatm),basend(jatm)
						tmatatm(iatm,jatm)=tmatatm(iatm,jatm)+tdmatbas(ibas,jbas)**2
					end do
				end do
				tmatatm(iatm,jatm)=dsqrt(tmatatm(iatm,jatm))
			end do
		end do
	end if
end if

tmatatmtmp=tmatatm
!Convert tmatatm to tmatatmnoh (the tmatatm without hydrogens), tmatatmtmp is used as a intermediate
ncennoh=count(a(:)%index/=1)
write(*,"(' The number of non-hydrogen atoms:',i10)") ncennoh
allocate(tmatatmnoh(ncennoh,ncennoh))
itmp=0
do iatm=1,ncenter
	if (a(iatm)%index/=1) then
		itmp=itmp+1
		tmatatmtmp(:,itmp)=tmatatmtmp(:,iatm)
		tmatatmtmp(itmp,:)=tmatatmtmp(iatm,:)
	end if
end do
tmatatmnoh(:,:)=tmatatmtmp(1:ncennoh,1:ncennoh)

write(*,"(' Sum of all elements (including hydrogens):',f16.8)") sum(tmatatm)
write(*,"(' Maximum and minimum (including hydrogens):',2f16.8)") maxval(tmatatm),minval(tmatatm)
write(*,"(' Sum of all elements (without hydrogens):  ',f16.8)") sum(tmatatmnoh)
write(*,"(' Maximum and minimum (without hydrogens):  ',2f16.8)") maxval(tmatatmnoh),minval(tmatatmnoh)

clrlimlow=0 !minval(tmatatmnoh)
clrlimhigh=maxval(tmatatmnoh)
ifhydrogen=0
ninterpo=10
nstepsize=ceiling(dfloat(ncennoh)/10D0)
ifnormsum=0
facnorm=1D0 !Normalization factor, default is 1, namely do not do normalization
stepsizez=(clrlimhigh/facnorm-clrlimlow/facnorm)/10
nlabdigz=4
labsize=55
if (allocated(frag)) then
	labsize=65
    nlabdigz=3
end if

do while(.true.)
	write(*,*)
	write(*,*) "     ============ Plotting atom or fragment transition matrix ============"
	if (allocated(frag)) then
		write(*,"(a,i3,a)") " -1 Clean fragment definition, current:",nfrag," fragments"
	else
		write(*,*) "-1 Define fragments"
	end if
	write(*,*) "0 Return"
	write(*,*) "1 Show transition matrix map"
	write(*,*) "2 Save transition matrix map to a graphical file in current folder"
	write(*,*) "3 Export transition density matrix to tmat.txt in current folder"
	if (.not.allocated(frag)) then
		if (ifhydrogen==0) write(*,*) "4 Toggle if taking hydrogens into account, current: No"
		if (ifhydrogen==1) write(*,*) "4 Toggle if taking hydrogens into account, current: Yes"
	end if
	write(*,"(a,f10.4,a,f10.4)") " 5 Change lower and upper limit of color scale, current:",clrlimlow/facnorm," to",clrlimhigh/facnorm
	if (.not.allocated(frag)) then
		write(*,"(a,i3)") " 6 Set the number of interpolation steps between grids, current:",ninterpo
		write(*,"(a,i3,a,f9.5)") " 7 Set stepsize between labels, current:  X/Y is ",nstepsize,", Z is",stepsizez
		if (ifnormsum==0) write(*,*) "8 Switch if normalizing the sum of all elements to unity, current: No"
		if (ifnormsum==1) write(*,*) "8 Switch if normalizing the sum of all elements to unity, current: Yes"
    else
		write(*,"(a,f9.5)") " 7 Set stepsize between labels in Z axis, current:",stepsizez
	end if
    write(*,"(a,a)") " 9 Set color transition, current: ",trim(clrtransname(iclrtrans))
    write(*,"(a,i4)") " 10 Set label size, current:",labsize
    write(*,"(a,i2)") " 11 Set number of decimal places in labels of Z-axis, current:",nlabdigz
	read(*,*) isel

	if (isel==0) then
		return
	else if (isel==-1) then
		if (allocated(frag)) then !Unset fragments
			deallocate(frag,fragnatm,tmatfrag)
			nfrag=0
			ninterpo=10
			clrlimlow=0
			clrlimhigh=maxval(tmatatmnoh)
			ifhydrogen=0
			stepsizez=(clrlimhigh-clrlimlow)/10 !Note that when entering fragment TDM mode, facnorm has been set to 1 and ifnormsum=0
			cycle
		else !Define fragments
			write(*,*)
			write(*,*) "How many fragments to be defined? e.g. 3"
			write(*,"(a)") " Note: If you input 0, then fragment definition will be loaded from an external file"
			read(*,*) nfrag
			if (nfrag==0) then
				write(*,*) "Input the file containing fragment definition, e.g. C:\myfrag.txt"
				write(*,"(a)") " If pressing ENTER button directly, fragdef.txt in current folder will be loaded"
				read(*,"(a)") c200tmp
				if (c200tmp==" ") c200tmp="fragdef.txt"
				inquire(file=c200tmp,exist=alive)
				if (alive) then
					open(10,file=c200tmp,status="old")
					nfrag=0
					do while(.true.)
						read(10,"(a)",iostat=ierror) c200tmp
						if (c200tmp==" ".or.ierror/=0) exit
						nfrag=nfrag+1
					end do
					write(*,"(' There are',i5,' fragments')") nfrag
					allocate(frag(nfrag,ncenter),fragnatm(nfrag))
					rewind(10)
					do ifrag=1,nfrag
						read(10,"(a)") c2000tmp
						call str2arr(c2000tmp,fragnatm(ifrag),frag(ifrag,:))
					end do
					close(10)
				else
					write(*,"(' Error: Unable to find ',a,/)") trim(c200tmp)
					cycle
				end if
			else
				allocate(frag(nfrag,ncenter),fragnatm(nfrag))
				do ifrag=1,nfrag
					write(*,"(' Input atomic indices for fragment',i4,', e.g. 1,4,8-12,15')") ifrag
					read(*,"(a)") c2000tmp
					call str2arr(c2000tmp,fragnatm(ifrag),frag(ifrag,:))
				end do
			end if
			
			do ifrag=1,nfrag
				write(*,"(' Atoms in fragment',i3,':')") ifrag
				write(*,"(15i5)") frag(ifrag,1:fragnatm(ifrag))
			end do
					
			allocate(tmatfrag(nfrag,nfrag))
			tmatfrag=0
			do ifrag=1,nfrag
				do jfrag=1,nfrag
					do idx=1,fragnatm(ifrag)
						iatm=frag(ifrag,idx)
						do jdx=1,fragnatm(jfrag)
							jatm=frag(jfrag,jdx)
							tmatfrag(ifrag,jfrag)=tmatfrag(ifrag,jfrag)+tmatatm(iatm,jatm)
						end do
					end do
				end do
			end do
			clrlimlow=0
			clrlimhigh=maxval(tmatfrag)
			write(*,"(a)") " Constructing fragment transition matrix has finished! Now you can plot it with options 1 and 2 or export it via option 3"
			ninterpo=1
			facnorm=1D0
			ifnormsum=0
			stepsizez=(clrlimhigh-clrlimlow)/10
		end if
		
	else if (isel==1.or.isel==2) then
		if (isel==2) isavepic=1
		if (allocated(frag)) then !Fragment transition matrix
			call drawmatcolor(tmatfrag/facnorm,nfrag,nfrag,1D0,dfloat(nfrag),1D0,dfloat(nfrag),&
            clrlimlow/facnorm,clrlimhigh/facnorm,1D0,1D0,ninterpo,-1,labsize,nlabdigz,stepsizez)
		else !Atom transition matrix
			if (ifhydrogen==1) then
				call drawmatcolor(tmatatm/facnorm,ncenter,ncenter,1D0,dfloat(ncenter),1D0,dfloat(ncenter),&
                clrlimlow/facnorm,clrlimhigh/facnorm,dfloat(nstepsize),dfloat(nstepsize),ninterpo,-1,labsize,nlabdigz,stepsizez)
			else if (ifhydrogen==0) then
				call drawmatcolor(tmatatmnoh/facnorm,ncennoh,ncennoh,1D0,dfloat(ncennoh),1D0,dfloat(ncennoh),&
                clrlimlow/facnorm,clrlimhigh/facnorm,dfloat(nstepsize),dfloat(nstepsize),ninterpo,-1,labsize,nlabdigz,stepsizez)
			end if
		end if
		if (isel==2) then
			isavepic=0
			write(*,*) "Done, the image has been saved to current folder with ""dislin"" prefix"
		end if
	else if (isel==3) then
		open(10,file="tmat.txt",status="replace")
		if (allocated(frag)) then !Fragment transition matrix
			do ifrag=1,nfrag
				do jfrag=1,nfrag
					write(10,"(2i8,f12.6)") ifrag,jfrag,tmatfrag(ifrag,jfrag)/facnorm
				end do
			end do
		else !Atom transition matrix
			if (ifhydrogen==0) then
				do iatm=1,ncennoh
					do jatm=1,ncennoh
						write(10,"(2i8,f12.6)") iatm,jatm,tmatatmnoh(iatm,jatm)/facnorm
					end do
				end do
			else if (ifhydrogen==1) then
				do iatm=1,ncenter
					do jatm=1,ncenter
						write(10,"(2i8,f12.6)") iatm,jatm,tmatatm(iatm,jatm)/facnorm
					end do
				end do
			end if
		end if
		close(10)
		write(*,*) "Done, the data have been exported to tmat.txt in current folder"
	else if (isel==4) then
		if (ifhydrogen==0) then
			ifhydrogen=1
			clrlimlow=0
			clrlimhigh=maxval(tmatatm)
		else if (ifhydrogen==1) then
			ifhydrogen=0
			clrlimlow=0
			clrlimhigh=maxval(tmatatmnoh)
		end if
	else if (isel==5) then
		write(*,*) "Input lower and upper limits, respectively, e.g. 0,1.5"
		read(*,*) clrlimlow,clrlimhigh
		clrlimlow=clrlimlow*facnorm
		clrlimhigh=clrlimhigh*facnorm
        stepsizez=(clrlimhigh/facnorm-clrlimlow/facnorm)/10
	else if (isel==6) then
		write(*,*) "Please input interpolation steps between grids, e.g. 2"
		write(*,"(a)") " Note: Larger value gives rise to smoother graph, 1 &
		&means do not perform interpolation, in this case each grid in the map corresponds to a matrix element"
		read(*,*) ninterpo
	else if (isel==7) then
		if (.not.allocated(frag)) then
			write(*,*) "Please input stepsize for X and Y axes, e.g. 4"
			read(*,*) nstepsize
        end if
		write(*,*) "Please input stepsize for Z axis, e.g. 0.05"
		read(*,*) stepsizez
	else if (isel==8) then
		if (ifnormsum==0) then
			ifnormsum=1
			facnorm=sum(tmatatm(:,:)) !The sum of all elements
			write(*,"(a,f16.8)") " The normalization factor is",facnorm
		else if (ifnormsum==1) then
			ifnormsum=0
			facnorm=1D0
		end if
	else if (isel==9) then
        call selcolortable
	else if (isel==10) then
		write(*,*) "Input text size, e.g. 55"
        read(*,*) labsize
    else if (isel==11) then
		write(*,*) "Input number of decimal places in labels of Z-axis, e.g. 3"
        read(*,*) nlabdigz
	end if
end do
end subroutine






!--------------------------------------------------------------------
!---------- Interfragment charger transfer (IFCT) analysis ----------
!--------------------------------------------------------------------
subroutine IFCT
use defvar
use excitinfo
use util
implicit real*8 (a-h,o-z)
real*8 atmhole(ncenter),atmele(ncenter),CTmatatm(ncenter,ncenter)
integer,allocatable :: frag(:,:),fragnatm(:)
real*8,allocatable :: CTmatfrag(:,:),fraghole(:),fragele(:)
character c2000tmp*2000,c200tmp*200,selectyn
logical alivehole,aliveele

inquire(file="hole.cub",exist=alivehole)
inquire(file="electron.cub",exist=aliveele)
write(*,*) "Select the method for calculating hole and electron composition"
write(*,*) "1 Mulliken-like partition (very fast, but incompatible with diffuse functions)"
write(*,*) "2 Hirshfeld partition (slow, but very robust)"
if (alivehole.and.aliveele) write(*,"(a)") " 3 Becke partition via interpolation based on hole.cub and electron.cub in current folder"
read(*,*) icompmethod

if (icompmethod==1.or.icompmethod==2) then !Mulliken-like or Hirshfeld partition for selected excitation
	call loadallexcinfo(1)
	call selexcit(istate)
	call loadexccoeff(istate,1)
	if (icompmethod==1) then
		call atmcontri_holeele_Mulliken(atmhole,atmele,1)
	else if (icompmethod==2) then
		call atmcontri_holeele_Hirshfeld(atmhole,atmele,1)
    end if
else if (icompmethod==3) then !Employ fuzzy partition to obtain atomic contributions based on interpolated function of existing cube files of hole and electron
	iuserfuncold=iuserfunc
	iuserfunc=-1
	call readcube("hole.cub",1,1)
	write(*,*) "Evaluating composition of hole, please wait..."
	call atmcontrifunc(1,atmhole,100)
	call readcube("electron.cub",1,1)
	write(*,*) "Evaluating composition of electron, please wait..."
	call atmcontrifunc(1,atmele,100)
	iuserfunc=iuserfuncold
end if

do while (.true.)
	write(*,*)
	write(*,*) "How many fragments to be defined? e.g. 3"
	write(*,"(a)") " Note: If you input 0, then fragment definition will be loaded from an external file. &
	&If you input -1, then atom-atom charge transfer matrix will be exported to atmCTmat.txt in current folder"
	read(*,*) nfrag
	if (nfrag==0) then
		write(*,*) "Input the file containing fragment definition, e.g. C:\fragdef.txt"
		write(*,"(a)") " Note: If pressing ENTER button directly, fragdef.txt in current folder will be used"
		do while(.true.)
			read(*,"(a)") c200tmp
			if (c200tmp==" ") c200tmp="fragdef.txt"
			inquire(file=c200tmp,exist=alive)
			if (alive) then
				open(10,file=c200tmp,status="old")
				nfrag=0
				do while(.true.)
					read(10,"(a)",iostat=ierror) c200tmp
					if (c200tmp==" ".or.ierror/=0) exit
					nfrag=nfrag+1
				end do
				write(*,"(' There are',i5,' fragments')") nfrag
				allocate(frag(nfrag,ncenter),fragnatm(nfrag))
				rewind(10)
				do ifrag=1,nfrag
					read(10,"(a)") c200tmp
					call str2arr(c200tmp,fragnatm(ifrag),frag(ifrag,:))
					write(*,"(' Atoms in fragment',i3,':')") ifrag
					write(*,"(15i5)") frag(ifrag,1:fragnatm(ifrag))
				end do
				close(10)
				exit
			else
				write(*,*) "Unable to find the file, input again"
			end if
		end do
	else if (nfrag==-1) then
		write(*,*) "Generating and exporting atmCTmat.txt..."
		write(*,*) "Note: Negative matrix elements will be set to zero (if any)"
		do iatm=1,ncenter
			do jatm=1,ncenter
				CTmatatm(iatm,jatm)=atmhole(iatm)*atmele(jatm)
				if (CTmatatm(iatm,jatm)<0) CTmatatm(iatm,jatm)=0
			end do
		end do
		open(10,file="atmCTmat.txt",status="replace")
		call showmatgau(CTmatatm,"atom-atom charge transfer matrix",0,"f14.8",10)
		close(10)
		excitfilename=" "
		write(*,"(a)") " Done! atom-atom charge transfer matrix has been exported to atmCTmat.txt in current folder, &
		&rows and columns correspond to hole and electron, respectively."
		write(*,"(a)") " Hint: You can use subfunction 2 of main function 18 to plot this matrix as heat map by using this file as input file"
		return
	else
		allocate(frag(nfrag,ncenter),fragnatm(nfrag))
		do ifrag=1,nfrag
			write(*,"(' Input atom indices for fragment',i4,', e.g. 1,4,8-12,15')") ifrag
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,fragnatm(ifrag),frag(ifrag,:))
		end do
	end if
	
	allocate(CTmatfrag(nfrag,nfrag),fraghole(nfrag),fragele(nfrag))
	write(*,*) "Contribution of each fragment to hole and electron:"
	do ifrag=1,nfrag
		fraghole(ifrag)=sum(atmhole(frag(ifrag,1:fragnatm(ifrag))))
		fragele(ifrag)=sum(atmele(frag(ifrag,1:fragnatm(ifrag))))
		write(*,"(i3,'  Hole:',f7.2,' %     Electron:',f7.2,' %')") ifrag,fraghole(ifrag)*100,fragele(ifrag)*100
	end do
	if (any(fraghole<0).or.any(fragele<0)) then
		write(*,"(a)") " Warning: There are unphysical negative contributions, &
		&they are regarded as zero during the calculation of interfragment charger-transfer matrix"
        if (icompmethod==1) write(*,"(a)") " Using Hirshfeld partition instead of Mulliken partition may be able to obtain more reliable result"
		where (fraghole<0) fraghole=0
		where (fragele<0) fragele=0
	end if
    totalLE_intrinsic=0
    totalCT_intrinsic=0
	do ifrag=1,nfrag
		do jfrag=1,nfrag
			CTmatfrag(ifrag,jfrag)=fraghole(ifrag)*fragele(jfrag)
            if (ifrag==jfrag) then
				totalLE_intrinsic=totalLE_intrinsic+CTmatfrag(ifrag,jfrag)
            else
				totalCT_intrinsic=totalCT_intrinsic+CTmatfrag(ifrag,jfrag)
            end if
		end do
	end do
	write(*,*) "Construction of interfragment charger-transfer matrix has finished!"
	write(*,*)

	do ifrag=1,nfrag
		varpop=sum(CTmatfrag(:,ifrag))-sum(CTmatfrag(ifrag,:))
		write(*,"(' Variation of population number of fragment',i3,':',f10.5)") ifrag,varpop
        if (ifrag==1) totalCT_apparent=abs(varpop)
	end do
	
	write(*,*)
	do ifrag=1,nfrag
		write(*,"(' Intrafragment electron redistribution of fragment',i3,':',f10.5)") ifrag,CTmatfrag(ifrag,ifrag)
	end do

	write(*,*)
	write(*,*) "Transferred electrons between fragments:"
	do ifrag=1,nfrag
		do jfrag=ifrag+1,nfrag
			write(*,"(i3,' ->',i3,':',f10.5,5x,i3,' <-',i3,':',f10.5,'     Net',i3,' ->',i3,':',f10.5)") &
			ifrag,jfrag,CTmatfrag(ifrag,jfrag),ifrag,jfrag,CTmatfrag(jfrag,ifrag),ifrag,jfrag,CTmatfrag(ifrag,jfrag)-CTmatfrag(jfrag,ifrag)
		end do
	end do
	deallocate(CTmatfrag,fragele,fraghole)
    write(*,*)
    write(*,"(a,f10.3,' %')") " Intrinsic charge transfer percentage, CT(%): ",totalCT_intrinsic*100
    write(*,"(a,f10.3,' %')") " Intrinsic local excitation percentage, LE(%):",totalLE_intrinsic*100
    if (nfrag==2) then
		write(*,"(a,f10.3,' %')") " Apparent charge transfer percentage, CT(%):  ",totalCT_apparent*100
		write(*,"(a,f10.3,' %')") " Apparent local excitation percentage, LE(%): ",100-totalCT_apparent*100
		write(*,*) "Note: See Section 3.21.8 of manual on the difference of the two types of CT(%)"
    end if
	
	write(*,*)
	write(*,*) "If you want to redefine fragments and recalculate data, input 1"
	write(*,*) "If you want to exit this function, input 0"
	read(*,*) isel
	if (isel==0) then
		return
	else
		deallocate(frag,fragnatm)
	end if
end do
end subroutine





!------------------------------------------------------------------------------------------------
!---------- Decompose transition dipole moment as molecular orbital pair contributions ----------
!------------------------------------------------------------------------------------------------
subroutine transdip_orbpair
use defvar
use util
use excitinfo
implicit real*8 (a-h,o-z)
real*8,allocatable :: tmparr(:)
integer,allocatable :: idxlist(:)
real*8,allocatable :: dipcontri(:,:) !(1/2/3,iexc) contribution of orbital pairs "iexc" to transition dipole moment in X/Y/Z
character strdir*3,strspini,strspinj,c80tmp*80

write(*,*) "Choose type of transition dipole moment to decompose:"
write(*,*) "1 Electric"
write(*,*) "2 Magnetic"
read(*,*) idiptype

if (idiptype==1) then
	if (.not.allocated(Dbas)) then
		write(*,*) "Generating electric dipole moment integral matrix between basis functions..."
		call genDbas_curr
	end if
	if (.not.allocated(DorbA)) then
		write(*,*) "Generating electric dipole moment integral matrix between MOs..."
		allocate(DorbA(3,nbasis,nbasis))
		if (allocated(CObasb)) allocate(DorbB(3,nbasis,nbasis))
		call genDorb
	end if
else if (idiptype==2) then
	if (.not.allocated(Magbas)) then
		write(*,*) "Generating magnetic dipole moment integral matrix between basis functions..."
		call genMagbas_curr
	end if
	if (.not.allocated(MagorbA)) then
		write(*,*) "Generating magnetic dipole moment integral matrix between MOs..."
		allocate(MagorbA(3,nbasis,nbasis))
		if (allocated(CObasb)) allocate(MagorbB(3,nbasis,nbasis))
		call genMagorb
	end if
end if

write(*,*)

!Load electronic excitation information
call loadallexcinfo(1)
call selexcit(istate)
call loadexccoeff(istate,1)

if (idiptype==1) write(*,*) "Calculating orbital pair contributions to transition electric dipole moment..."
if (idiptype==2) write(*,*) "Calculating orbital pair contributions to transition magnetic dipole moment..."
allocate(dipcontri(3,excnorb)) ; dipcontri=0
if (wfntype==0.or.wfntype==3) then
	fac=2
else
	fac=1
end if
do iexcitorb=1,excnorb
	imo=orbleft(iexcitorb)
	jmo=orbright(iexcitorb)
    if (idiptype==1) then !Electric
		if (imo<=nbasis) then
			dipcontri(:,iexcitorb)=DorbA(:,imo,jmo)*exccoeff(iexcitorb)*fac
		else !beta->beta
			dipcontri(:,iexcitorb)=DorbB(:,imo-nbasis,jmo-nbasis)*exccoeff(iexcitorb)*fac
		end if
    else if (idiptype==2) then !Magnetic, see Eq 24 of "Abinitio calculations of oscillator and rotatory strengths in the random‐phase approximatio-Twisted mono‐olefins"
		tmp=1
        if (excdir(iexcitorb)==2) tmp=-1
		if (imo<=nbasis) then
            dipcontri(:,iexcitorb)=tmp*MagorbA(:,imo,jmo)*exccoeff(iexcitorb)*fac
		else !beta->beta
			dipcontri(:,iexcitorb)=tmp*MagorbB(:,imo-nbasis,jmo-nbasis)*exccoeff(iexcitorb)*fac
		end if
    end if
end do

xdipall=sum(dipcontri(1,:))
ydipall=sum(dipcontri(2,:))
zdipall=sum(dipcontri(3,:))
dipallnorm=dsqrt(xdipall**2+ydipall**2+zdipall**2)
write(*,*)
write(*,"(a,/)") " Note: Carbon, 165, 461 (2020) employed this function to study the nature of very strong absorption of cyclo[18]carbon, &
&you are suggested to look at this paper and cite it along with Multiwfn original paper"
if ((naelec==nbelec).and.excmulti==3) write(*,"(a,/)") " Notice: Since the spin multiplicity between ground state and the excited state &
&is different (spin-forbidden), the transition dipole moment analyzed in this function only considers spatial part"
write(*,"(' Transition dipole moment in X/Y/Z: ',3f11.6,' a.u.')") xdipall,ydipall,zdipall
write(*,"(' Norm of transition dipole moment:  ',f11.6,' a.u.')") dipallnorm
if (idiptype==1) then
	oscillstr=2D0/3D0*excene/au2eV*(xdipall**2+ydipall**2+zdipall**2)
	write(*,"(' Oscillator strength:',f12.7)") oscillstr
end if

do while(.true.)
	write(*,*)
    call menutitle("Show orbital pair contributions",10,1)
	write(*,*) "0 Return"
	write(*,*) "1 Print orbital pairs having contribution larger than a threshold"
	write(*,*) "2 Print orbital pairs in the order of absolute contribution to X component"
	write(*,*) "3 Print orbital pairs in the order of absolute contribution to Y component"
	write(*,*) "4 Print orbital pairs in the order of absolute contribution to Z component"
	write(*,*) "5 Print orbital pairs in the order of norm of contribution vector"
	write(*,*) "10 Export contribution of all orbital pairs to transdip.txt in current folder"
	read(*,*) isel
	
	if (isel==0) then
		return
	else if (isel==1) then
		write(*,"(/,a)") " Input the threshold for printing, e.g. 0.01 means the orbital pairs having contribution &
		&to any component of transition dipole moment >= 0.01 will be printed"
		read(*,*) printthres
		xdipsum=0; ydipsum=0; zdipsum=0
		nshownpair=0
		if ( any(abs(dipcontri)>printthres) ) then
			write(*,*) " #Pair   Orbital trans. Coefficient      Transition dipole X/Y/Z   Norm (a.u.)"
			do iexcitorb=1,excnorb
				if ( any(abs(dipcontri(:,iexcitorb))>printthres) ) then
					nshownpair=nshownpair+1
					xdipsum=xdipsum+dipcontri(1,iexcitorb)
					ydipsum=ydipsum+dipcontri(2,iexcitorb)
					zdipsum=zdipsum+dipcontri(3,iexcitorb)
					imo=orbleft(iexcitorb)
					jmo=orbright(iexcitorb)
					strdir=" ->"
					if (excdir(iexcitorb)==2) strdir=" <-"
					if (wfntype==0.or.wfntype==3) then
						write(*,"(i7,i7,a,i7,f11.6,4f11.6)") iexcitorb,imo,strdir,jmo,exccoeff(iexcitorb),dipcontri(:,iexcitorb),dsqrt(sum(dipcontri(:,iexcitorb)**2))
					else
						strspini="A"
						if (imo>nbasis) then
							strspini="B"
							imo=imo-nbasis
						end if
						strspinj="A"
						if (jmo>nbasis) then
							strspinj="B"
							jmo=jmo-nbasis
						end if
						write(*,"(i7,i6,a,a,i6,a,f11.6,4f11.6)") iexcitorb,imo,strspini,strdir,jmo,strspinj,exccoeff(iexcitorb),dipcontri(:,iexcitorb),dsqrt(sum(dipcontri(:,iexcitorb)**2))
					end if
				end if
			end do
			write(*,"(' Sum of the above',i8,' pairs:   ',3f11.6)") nshownpair,xdipsum,ydipsum,zdipsum
		else
			write(*,*) "No orbital pair statisfying the condition was found!"
		end if
		
	else if (isel==2.or.isel==3.or.isel==4.or.isel==5) then
		write(*,*) "Sorting, please wait..."
		allocate(tmparr(excnorb),idxlist(excnorb))
		if (isel==2) tmparr=dipcontri(1,:)
		if (isel==3) tmparr=dipcontri(2,:)
		if (isel==4) tmparr=dipcontri(3,:)
		if (isel==5) then
			do iexcitorb=1,excnorb
				tmparr(iexcitorb)=dsqrt(sum(dipcontri(:,iexcitorb)**2))
            end do
        end if
		forall (i=1:excnorb) idxlist(i)=i
		call sortr8(tmparr,"abs",idxlist)
		call invarri4(idxlist)
		if (isel==2) write(*,"(a,/)") " The orbital pairs have been sorted according to absolute value of X component of transition dipole moment"
		if (isel==3) write(*,"(a,/)") " The orbital pairs have been sorted according to absolute value of Y component of transition dipole moment"
		if (isel==4) write(*,"(a,/)") " The orbital pairs have been sorted according to absolute value of Z component of transition dipole moment"
		if (isel==5) write(*,"(a,/)") " The orbital pairs have been sorted according to norm of contribution vector to transition dipole moment"
		write(*,*) "How many orbital pairs to output? e.g. 20"
		read(*,*) noutpair
		if (noutpair>excnorb) noutpair=excnorb
		write(*,*) " #Pair   Orbital trans. Coefficient      Transition dipole X/Y/Z   Norm (a.u.)"
		do idx=1,noutpair
			iexcitorb=idxlist(idx)
			imo=orbleft(iexcitorb)
			jmo=orbright(iexcitorb)
			strdir=" ->"
			if (excdir(iexcitorb)==2) strdir=" <-"
			if (wfntype==0.or.wfntype==3) then
				write(*,"(i7,i7,a,i7,f11.6,4f11.6)") iexcitorb,imo,strdir,jmo,exccoeff(iexcitorb),dipcontri(:,iexcitorb),dsqrt(sum(dipcontri(:,iexcitorb)**2))
			else
				strspini="A"
				if (imo>nbasis) then
					strspini="B"
					imo=imo-nbasis
				end if
				strspinj="A"
				if (jmo>nbasis) then
					strspinj="B"
					jmo=jmo-nbasis
				end if
				write(*,"(i7,i6,a,a,i6,a,f11.6,4f11.6)") iexcitorb,imo,strspini,strdir,jmo,strspinj,exccoeff(iexcitorb),dipcontri(:,iexcitorb),dsqrt(sum(dipcontri(:,iexcitorb)**2))
			end if
		end do
		deallocate(tmparr,idxlist)
		
	else if (isel==10) then
		open(10,file="transdip.txt",status="replace")
		write(10,*) " #Pair   Orbital trans. Coefficient      Transition dipole X/Y/Z   Norm (a.u.)"
		do iexcitorb=1,excnorb
			imo=orbleft(iexcitorb)
			jmo=orbright(iexcitorb)
			strdir=" ->"
			if (excdir(iexcitorb)==2) strdir=" <-"
			if (wfntype==0.or.wfntype==3) then
				write(10,"(i7,i7,a,i7,f11.6,4f11.6)") iexcitorb,imo,strdir,jmo,exccoeff(iexcitorb),dipcontri(:,iexcitorb),dsqrt(sum(dipcontri(:,iexcitorb)**2))
			else
				strspini="A"
				if (imo>nbasis) then
					strspini="B"
					imo=imo-nbasis
				end if
				strspinj="A"
				if (jmo>nbasis) then
					strspinj="B"
					jmo=jmo-nbasis
				end if
				write(10,"(i7,i6,a,a,i6,a,f11.6,4f11.6)") iexcitorb,imo,strspini,strdir,jmo,strspinj,exccoeff(iexcitorb),dipcontri(:,iexcitorb),dsqrt(sum(dipcontri(:,iexcitorb)**2))
			end if
		end do
		close(10)
		write(*,*) "Done! The data have been exported to transdip.txt in current folder!"
	end if
end do

end subroutine






!--------------------------------------------------------------------------------------------------
!---------- Generate transition density matrix between ground state and an excited state ----------
!--------------------------------------------------------------------------------------------------
!There are two ways to construct TDM for TD method, see eqs. 22~24 in &
!JCP,66,3460 (Abinitio calculations of oscillator and rotatory strengths in the random‐phase approximatio-Twisted mono‐olefins) for detail
!iTDMtype=1: Correct for transition electric dipole moment, excitation and de-excitation are not distinguished
!iTDMtype=2: Correct for transition velocity/magnetic dipole moment, excitation and de-excitation are considered individually
!isym=1: Ask user if symmetrizing the TDM.  =3: Skip symmetrizing
!For SF-TDDFT, this routine is meaningless because TDM between reference state and SF-TDDFT states must be zero due to spin orthogonal
subroutine genTDM(iTDMtype,isym)
use defvar
use excitinfo
implicit real*8 (a-h,o-z)
character c80tmp*80
real*8,allocatable :: tmpmat(:,:),tmparr(:) !Arrays for temporary use

if (.not.allocated(CObasa)) then
	write(*,"(a)") " Error: TDM is unable to be generated because basis function information is not available! &
	&Please carefully read corresponding sections of the manual to make clear which kinds of input files should be used!"
	write(*,*) "Press ENTER button to exit program"
	read(*,*)
	stop
end if

if (iTDMtype==1) write(*,*) "Generating transition density matrix of type 1..."
if (iTDMtype==2) write(*,*) "Generating transition density matrix of type 2..."
if (.not.allocated(tdmata)) allocate(tdmata(nbasis,nbasis))
if ((wfntype==1.or.wfntype==4).and.(.not.allocated(tdmatb))) allocate(tdmatb(nbasis,nbasis))

!In order to significantly accelerate generation speed of TDM, for each occupied MO, we first accumulate coefficient vector of all related virtual MOs as tmparr
allocate(tmparr(nbasis),tmpmat(1,nbasis))
tdmata=0
do imo=1,nint(naelec) !Closed-shell or alpha part of open-shell
	tmparr=0
	if (iTDMtype==1) then
		do iexcitorb=1,excnorb
			if (orbleft(iexcitorb)==imo) tmparr=tmparr+exccoeff(iexcitorb)*CObasa(:,orbright(iexcitorb))
		end do
	else if (iTDMtype==2) then
		do iexcitorb=1,excnorb
			if (orbleft(iexcitorb)==imo) then
				if (excdir(iexcitorb)==1) then
					tmparr=tmparr+exccoeff(iexcitorb)*CObasa(:,orbright(iexcitorb))
				else
					tmparr=tmparr-exccoeff(iexcitorb)*CObasa(:,orbright(iexcitorb))
				end if
			end if
		end do
	end if
	tmpmat(1,:)=tmparr(:)
	tdmata=tdmata+matmul(CObasa(:,imo:imo),tmpmat)
end do

if (wfntype==0.or.wfntype==3) then
	tdmata=tdmata*2 !Closed-shell, double the current TDM
else if (wfntype==1.or.wfntype==4) then !Beta part of open-shell
	tdmatb=0
	do imo=nbasis+1,nbasis+nint(nbelec)
		tmparr=0
		if (iTDMtype==1) then
			do iexcitorb=1,excnorb
				if (orbleft(iexcitorb)==imo) tmparr=tmparr+exccoeff(iexcitorb)*CObasb(:,orbright(iexcitorb)-nbasis)
			end do
		else if (iTDMtype==2) then
			do iexcitorb=1,excnorb
				if (orbleft(iexcitorb)==imo) then
					if (excdir(iexcitorb)==1) then
						tmparr=tmparr+exccoeff(iexcitorb)*CObasb(:,orbright(iexcitorb)-nbasis)
					else
						tmparr=tmparr-exccoeff(iexcitorb)*CObasb(:,orbright(iexcitorb)-nbasis)
					end if
				end if
			end do
		end if
		tmpmat(1,:)=tmparr(:)
		imob=imo-nbasis
		tdmatb=tdmatb+matmul(CObasb(:,imob:imob),tmpmat)
	end do
end if
deallocate(tmparr,tmpmat)

!! Below codes are used to check transition properties based on transition density matrix and corresponding integral matrix
!Check transition electric dipole moment. The result is correct when iTDMtype==1, see above
! 		if (.not.allocated(Dbas)) call genDbas_curr
! 		Teledipx=sum(tdmata*Dbas(1,:,:))
! 		Teledipy=sum(tdmata*Dbas(2,:,:))
! 		Teledipz=sum(tdmata*Dbas(3,:,:))
! 		Teledipx=sum((tdmata+tdmatb)*Dbas(1,:,:))
! 		Teledipy=sum((tdmata+tdmatb)*Dbas(2,:,:))
! 		Teledipz=sum((tdmata+tdmatb)*Dbas(3,:,:))
! 		write(*,"(' Transition electric dipole moment:',3f12.6)") Teledipx,Teledipy,Teledipz

! 		!Check transition velocity dipole moment. The result is correct when iTDMtype==2, see above
! 		if (.not.allocated(Velbas)) call genvelbas_curr
! 		Tvdipx=sum(tdmata*Velbas(1,:,:))
! 		Tvdipy=sum(tdmata*Velbas(2,:,:))
! 		Tvdipz=sum(tdmata*Velbas(3,:,:))
! 		write(*,"(' Transition velocity dipole moment:',3f12.6)") Tvdipx,Tvdipy,Tvdipz

! 		!Check transition magnetic dipole moment. The result is correct when iTDMtype==2, see above
! 		if (.not.allocated(Magbas)) call genMagbas_curr
! 	    end if
! 	    call showmatgau(tdmata,"tdmata matrix",1)
! 		Tmagdipx=sum(tdmata*Magbas(1,:,:))
! 		Tmagdipy=sum(tdmata*Magbas(2,:,:))
! 		Tmagdipz=sum(tdmata*Magbas(3,:,:))
! 		write(*,"(' Transition magnetic dipole moment:',3f12.6)") Tmagdipx,Tmagdipy,Tmagdipz

! 		!Below formula is absolutely correct, however the printed result is not correct here because 
! 		!we cannot obtain correct transition electric and magnetic dipole moments based on the same TDM.
! 		!If we directly use the value outputted by Gaussian, you will find the result is in line with the rotatory strength printed by Gaussian
! 		!Four points:
! 		!1) The negative sign: Because we ignored the negative sign when evaluating magnetic integrals
! 		!2) Diveded by two: Necessary by definition
! 		!3) 2.54174619D-018: convert electric dipole moment from a.u. to cgs, see gabedit build-in converter
! 		!4) 1.85480184D-020: convert magnetic dipole moment from a.u. to cgs, see gabedit build-in converter
! 		Rlen=-(Teledipx*Tmagdipx+Teledipy*Tmagdipy+Teledipz*Tmagdipz)/2D0*2.54174619D-018*1.85480184D-020 *1D40
! 		write(*,"(' Rotatory strength in length representation:',f14.8,' 10^-40 cgs')") Rlen
! 		cycle

!The TDM generated above is absolutely correct, irrespective of R and U references
!The TDM symmetrized as follows is exactly identical to the Gaussian output using IOp(6/8=3), irrespective of R and U references
!The reason of Gaussian (and many papers) meantime introducing the 1/sqrt(2) (rather than 1/2 as expected) into &
!the symmetrized TDM expression is completely unclear, it should be unresonable, however here we follow this convention
if (isym==1) then
	write(*,*)
	write(*,*) "If symmetrizing the transition density matrix?"
	write(*,*) "0 or n: Do not symmetrize"
	write(*,*) "1 or y: Symmetrize as TDM_sym(i,j)=[TDM(i,j)+TDM(j,i)]/2"
	write(*,*) "2:      Symmetrize as TDM_sym(i,j)=[TDM(i,j)+TDM(j,i)]/sqrt(2)"
	read(*,*) c80tmp
	if (c80tmp(1:1)=='1'.or.c80tmp(1:1)=='y'.or.c80tmp(1:1)=='2') then
		if (c80tmp(1:1)=='1'.or.c80tmp(1:1)=='y') tmpfac=2
		if (c80tmp(1:1)=='2') tmpfac=dsqrt(2D0)
		do ibas=1,nbasis
			do jbas=ibas,nbasis
				tdmata(ibas,jbas)=(tdmata(ibas,jbas)+tdmata(jbas,ibas))/tmpfac
				tdmata(jbas,ibas)=tdmata(ibas,jbas)
				if (wfntype==1.or.wfntype==4) then
					tdmatb(ibas,jbas)=(tdmatb(ibas,jbas)+tdmatb(jbas,ibas))/tmpfac
					tdmatb(jbas,ibas)=tdmatb(ibas,jbas)
				end if
			end do
		end do
		write(*,*) "The transition density matrix has been symmetrized"
	end if
end if

end subroutine





!---------------------------------------------------------------------------------------------------------------
!---------- Generate transition density matrix between two excited states and store to tdmata/tdmatab ----------
!---------------------------------------------------------------------------------------------------------------
!isymmetry controls how to symmetrize the matrix. If =-1, then user will be asked to choose a symmetrization method
!SF-TDDFT case is not supported yet (in that case, closed-shell code in this routine should be used, while both CObasa and CObasb are involved)
subroutine genTDM_2exc(istate,jstate,isymmetry)
use defvar
use util
use excitinfo
implicit real*8 (a-h,o-z)
integer istate,jstate,isymmetry
real*8 CObasa_tr(nbasis,nbasis),CObasb_tr(nbasis,nbasis)
character c80tmp*80

!1E-5 is lowest acceptable threshold. Cost of 1E-6 will be higher than 1E-5 by one order of magnitude
write(*,"(a)") " Input the threshold of product of two configuration coefficients for skipping configurations, e.g. 1E-6"
write(*,"(a)") " Note: If you press ENTER button directly, then 0.00001 will be used, which is a good compromise between cost and accuracy"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
	thres=1D-5
else
	read(c80tmp,*) thres
end if

write(*,*) "Calculating, please wait patiently..."
call walltime(iwalltime1)

!I have tried to parallize this section in various ways, however the calculation cost is weirdly higher,
!the reason is unclear, probably due to e.g. intensive memory use
nterm=0
if (.not.allocated(tdmata)) allocate(tdmata(nbasis,nbasis))
tdmata=0
CObasa_tr=transpose(CObasa)
if (wfntype==0.or.wfntype==3) then !Closed-shell case
	do ipair=1,allexcnorb(istate)
		call showprog(ipair,allexcnorb(istate))
		imo=allorbleft(ipair,istate)
		lmo=allorbright(ipair,istate)
		do jpair=1,allexcnorb(jstate)
			jmo=allorbleft(jpair,jstate)
			kmo=allorbright(jpair,jstate)
			wei=allexccoeff(ipair,istate)*allexccoeff(jpair,jstate)
			if (abs(wei)<thres) cycle !For saving time
			!For two |X+Y> wavefunctions, X only acts with X, and Y only acts with Y
			!Assume what we calculated is <X-Y|-r|X+Y> = <X|-r|X> - <Y|-r|Y>. XX and YY are treated in the same way 
			!For TDA case, Y is nonexistent, the result for the same two states is exactly identical to density=rhoci outputted by Gaussian, &
			!while for TD case, the wavefunction is not rigorously defined, and the result is not completely in line with Gaussian, &
			!but the difference is usually marginal (I am not sure which formula is employed by Gaussian)
			!The way of evaluating transition dipole momemnt used here follows the rule in Section 3.21.9 of manual.
			if (allexcdir(ipair,istate)/=allexcdir(jpair,jstate)) cycle
			nterm=nterm+1
			if (allexcdir(ipair,istate)==2) wei=-wei
			if (imo==jmo.and.lmo/=kmo) then
				tdmata=tdmata+wei*matmul(CObasa(:,lmo:lmo),CObasa_tr(kmo:kmo,:))
			else if (imo/=jmo.and.lmo==kmo) then
				tdmata=tdmata-wei*matmul(CObasa(:,jmo:jmo),CObasa_tr(imo:imo,:))
			else if (imo==jmo.and.lmo==kmo) then
				tdmata=tdmata+wei*( Ptot - matmul(CObasa(:,imo:imo),CObasa_tr(imo:imo,:)) + matmul(CObasa(:,lmo:lmo),CObasa_tr(lmo:lmo,:)) )
			end if
		end do
	end do
	tdmata=tdmata*2
	
else !Open-shell case
	if (.not.allocated(tdmatb)) allocate(tdmatb(nbasis,nbasis))
	tdmatb=0
	CObasb_tr=transpose(CObasb)
	do ipair=1,allexcnorb(istate)
		if (ipair>iprog) then
			call showprog(ipair,allexcnorb(istate))
			iprog=iprog+nstep
		end if
		imo=allorbleft(ipair,istate)
		lmo=allorbright(ipair,istate)
		do jpair=1,allexcnorb(jstate)
			jmo=allorbleft(jpair,jstate)
			kmo=allorbright(jpair,jstate)
			wei=allexccoeff(ipair,istate)*allexccoeff(jpair,jstate)
			if (abs(wei)<thres) cycle
			if (allexcdir(ipair,istate)/=allexcdir(jpair,jstate)) cycle
			if (imo<=nbasis.and.jmo>nbasis) cycle !Skip interaction between A->A and B->B due to spin forbidden
			nterm=nterm+1
			if (allexcdir(ipair,istate)==2) wei=-wei
			if (imo<=nbasis) then !A->A interacts with A->A
				if (imo==jmo.and.lmo/=kmo) then
					tdmata=tdmata+wei*matmul(CObasa(:,lmo:lmo),CObasa_tr(kmo:kmo,:))
				else if (imo/=jmo.and.lmo==kmo) then
					tdmata=tdmata-wei*matmul(CObasa(:,imo:imo),CObasa_tr(jmo:jmo,:))
				else if (imo==jmo.and.lmo==kmo) then
					tdmata=tdmata+wei*( Palpha - matmul(CObasa(:,imo:imo),CObasa_tr(imo:imo,:)) + matmul(CObasa(:,lmo:lmo),CObasa_tr(lmo:lmo,:)) )
				end if
			else !B->B interacts with B->B
				if (imo==jmo.and.lmo/=kmo) then
					llmo=lmo-nbasis
					kkmo=kmo-nbasis
					tdmatb=tdmatb+wei*matmul(CObasb(:,llmo:llmo),CObasb_tr(kkmo:kkmo,:))
				else if (imo/=jmo.and.lmo==kmo) then
					iimo=imo-nbasis
					jjmo=jmo-nbasis
					tdmatb=tdmatb-wei*matmul(CObasb(:,iimo:iimo),CObasb_tr(jjmo:jjmo,:))
				else if (imo==jmo.and.lmo==kmo) then
					iimo=imo-nbasis
					llmo=lmo-nbasis
					tdmatb=tdmatb+wei*( Pbeta - matmul(CObasb(:,iimo:iimo),CObasb_tr(iimo:iimo,:)) + matmul(CObasb(:,llmo:llmo),CObasb_tr(llmo:llmo,:)) )
				end if
			end if
		end do
	end do
end if

call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
write(*,"(' Totally',i12,' terms were calculated')") nterm

if (isymmetry==-1) then
    write(*,*)
    write(*,*) "If symmetrizing the transition density matrix?"
    write(*,*) "0: Do not symmetrize"
    write(*,*) "1: Symmetrize as TDM_sym(i,j)=[TDM(i,j)+TDM(j,i)]/2"
    write(*,*) "2: Symmetrize as TDM_sym(i,j)=[TDM(i,j)+TDM(j,i)]/sqrt(2)" !Used by Gaussian, weird!
    read(*,*) isymmmethod
else
	isymmmethod=isymmetry
end if
if (isymmmethod/=0) then
	if (isymmmethod==1) tmpfac=2
	if (isymmmethod==2) tmpfac=dsqrt(2D0)
	do ibas=1,nbasis
		do jbas=ibas,nbasis
			tdmata(ibas,jbas)=(tdmata(ibas,jbas)+tdmata(jbas,ibas))/tmpfac
			tdmata(jbas,ibas)=tdmata(ibas,jbas)
			if (wfntype==1.or.wfntype==4) then
				tdmatb(ibas,jbas)=(tdmatb(ibas,jbas)+tdmatb(jbas,ibas))/tmpfac
				tdmatb(jbas,ibas)=tdmatb(ibas,jbas)
			end if
		end do
	end do
	write(*,*) "The transition density matrix has been symmetrized"
end if
end subroutine






!-------------------------------------------------------------------------
!---------- Export transition density matrix to plain text file ----------
!-------------------------------------------------------------------------
subroutine exportTDM
use defvar
use excitinfo
use util
implicit real*8 (a-h,o-z)
character selectyn

write(*,*) "1 Generate transition density matrix between ground state and excited state"
write(*,*) "2 Generate transition density matrix between two excited states"
read(*,*) isel

call loadallexcinfo(1)
if (isel==1) then
	write(*,*) "Please input index of the excited state, e.g. 2"
	read(*,*) istate
	call loadexccoeff(istate,1)
	call genTDM(1,1)
else if (isel==2) then
	call loadallexccoeff(0)
	write(*,*) "Input index of the two excited states, e.g. 2,5"
	read(*,*) istate,jstate
	call genTDM_2exc(istate,jstate,-1)
end if
write(*,*)
if (wfntype==0.or.wfntype==3) then
	write(*,*) "Exporting tdmat.txt..."
	open(10,file="tdmat.txt",status="replace")
	call showmatgau(tdmata,"Transition density matrix",0,"f14.8",10)
	close(10)
	write(*,"(a)") " Done! The transition density matrix has been outputted to tdmat.txt in current folder"
	write(*,"(a)") " Hint: You can use function 2 of main function 18 to plot the transition density matrix by using this file as input file"
else
	write(*,*) "Exporting tdmat.txt, tdmata.txt and tdmatb.txt..."
	open(10,file="tdmat.txt",status="replace")
	call showmatgau(tdmata+tdmatb,"Total transition density matrix",0,"f14.8",10)
	close(10)
	open(10,file="tdmata.txt",status="replace")
	call showmatgau(tdmata,"Alpha transition density matrix",0,"f14.8",10)
	close(10)
	open(10,file="tdmatb.txt",status="replace")
	call showmatgau(tdmatb,"Beta transition density matrix",0,"f14.8",10)
	close(10)
	write(*,"(a)") " Done! The total, alpha and beta transition density matrix has been outputted to tdmat.txt, tdmata.txt and tdmatb.txt in current folder, respectively"
	write(*,"(a)") " Hint: You can use function 2 of main function 18 to plot the transition density matrix by using these files as input file"
end if

write(*,"(/,a)") " If outputting TDM.fch file in current folder, whose density matrix field will correspond to the newly generated transition density matrix? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	if (wfntype==0.or.wfntype==2.or.wfntype==3) then !Closed-shell case
		Ptot=tdmata
	else
		Palpha=tdmata
		Pbeta=tdmatb
		Ptot=tdmata+tdmatb
	end if
	call outfch("TDM.fch",10,1)
	write(*,"(' Reloading: ',a)") trim(firstfilename)
	call dealloall(0)
	call readinfile(firstfilename,1)
end if
end subroutine





!--------------------------------------------------------------------------------------------------------------------
!  Decompose transition electric/magnetic dipole moment as basis function and atom contributions (Mulliken partition)
!--------------------------------------------------------------------------------------------------------------------
subroutine transdip_basatm
use defvar
use excitinfo
use util
implicit real*8 (a-h,o-z)
character selectyn
real*8 bastrdip(3,nbasis) !Contribution from each basis function to transition dipole moment, the first index 1/2/3=X/Y/Z
real*8 trdipmatbas(3,nbasis,nbasis),trdipmatatm(3,ncenter,ncenter) !Transition dipole moment matrix in basis function / atom representation, the first index 1/2/3=X/Y/Z
real*8,pointer :: tmpbas(:,:,:)
character c80tmp*80

call loadallexcinfo(1)
write(*,*) "Calculate the transition dipole moment between which two states?"
write(*,*) "If you input e.g. 3, then ground state and excited state 3 will be calculated"
write(*,*) "If you input e.g. 2,5, then excited states 2 and 5 will be calculated"
read(*,"(a)") c80tmp
write(*,*)
write(*,*) "Decompose which type of transition dipole moment?"
write(*,*) "1: Electric      2: Magnetic"
read(*,*) idecomptype

!Generate TDM without symmetrization, because symmetrization will be done later
if (index(c80tmp,',')/=0) then
    if (idecomptype==2) then
        write(*,"(a)") " Error: Calculating magnetic transition dipole moment between two excited states is not supported!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
	read(c80tmp,*) istate,jstate
	call loadallexccoeff(0)
	call genTDM_2exc(istate,jstate,0) !Generate TDM between two excited states, without symmetrization
else
    read(c80tmp,*) istate
	call loadexccoeff(istate,1)
    call genTDM(idecomptype,3) !Generate TDM using configuration coefficients in memory, without symmetrization
end if

!Haven't calculated dipole moment integral matrix, so reload the input file and calculate it
if (idecomptype==1.and.(.not.allocated(Dbas))) then
    write(*,*) "Generating electric dipole moment integral matrix..."
	call genDbas_curr
else if (idecomptype==2.and.(.not.allocated(Magbas))) then
    write(*,*) "Generating magnetic dipole moment integral matrix..."
	call genMagbas_curr
end if

if (idecomptype==1) then
	tmpbas=>Dbas
else if (idecomptype==2) then
	tmpbas=>Magbas
end if

if (wfntype==0.or.wfntype==3) then !Closed-shell
    trdipmatbas(1,:,:)=tdmata(:,:)*tmpbas(1,:,:)
    trdipmatbas(2,:,:)=tdmata(:,:)*tmpbas(2,:,:)
    trdipmatbas(3,:,:)=tdmata(:,:)*tmpbas(3,:,:)
else
    write(*,*) "Analyze which spin?   1=Alpha  2=Beta  3=Total"
    read(*,*) isel
    if (isel==1) then
	    trdipmatbas(1,:,:)=tdmata(:,:)*tmpbas(1,:,:)
	    trdipmatbas(2,:,:)=tdmata(:,:)*tmpbas(2,:,:)
	    trdipmatbas(3,:,:)=tdmata(:,:)*tmpbas(3,:,:)
    else if (isel==2) then
	    trdipmatbas(1,:,:)=tdmatb(:,:)*tmpbas(1,:,:)
	    trdipmatbas(2,:,:)=tdmatb(:,:)*tmpbas(2,:,:)
	    trdipmatbas(3,:,:)=tdmatb(:,:)*tmpbas(3,:,:)
    else if (isel==3) then
	    trdipmatbas(1,:,:)=(tdmata(:,:)+tdmatb(:,:))*tmpbas(1,:,:)
	    trdipmatbas(2,:,:)=(tdmata(:,:)+tdmatb(:,:))*tmpbas(2,:,:)
	    trdipmatbas(3,:,:)=(tdmata(:,:)+tdmatb(:,:))*tmpbas(3,:,:)
    end if
end if

open(10,file="trdipcontri.txt",status="replace")
write(10,"(a)") " Contribution of basis functions to transition dipole moment (X,Y,Z, in a.u.) derived by Mulliken partition:"
bastrdip=0
do ibas=1,nbasis !Note that trdipmatbas is not strictly a symmetry matrix, so we get average value for off-diagonal elements
	do jbas=1,nbasis
		if (ibas==jbas) then
			bastrdip(:,ibas)=bastrdip(:,ibas)+trdipmatbas(:,ibas,jbas)
		else
			bastrdip(:,ibas)=bastrdip(:,ibas)+(trdipmatbas(:,ibas,jbas)+trdipmatbas(:,jbas,ibas))/2
		end if
	end do
	write(10,"(i5,'#    Shell:',i5,'   Center:',i5,'(',a2,')   Type:',a,3f12.6)") &
	ibas,basshell(ibas),bascen(ibas),a(bascen(ibas))%name,GTFtype2name(bastype(ibas)),bastrdip(:,ibas)
end do
write(10,*)
write(10,"(a)") " Contribution from atoms to transition dipole moment (X,Y,Z, in a.u.) derived by Mulliken partition:"
do iatm=1,ncenter
	if (basstart(iatm)==0) then
		xval=0;yval=0;zval=0
    else
		xval=sum(bastrdip(1,basstart(iatm):basend(iatm)))
		yval=sum(bastrdip(2,basstart(iatm):basend(iatm)))
		zval=sum(bastrdip(3,basstart(iatm):basend(iatm)))
    end if
	write(10,"(i5,'(',a2,'):',3f12.6)") iatm,a(iatm)%name,xval,yval,zval
end do
write(10,*)
write(10,"(a,3f12.6,a)") " Transition dipole moment in X/Y/Z",sum(bastrdip(1,:)),sum(bastrdip(2,:)),sum(bastrdip(3,:))," a.u."
close(10)
write(*,*)
write(*,*) "Done! The result has been outputted to trdipcontri.txt in current folder"

write(*,*)
write(*,"(a)") " Would you also like to output atom transition dipole moment matrix in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	trdipmatatm=0
	do iatm=1,ncenter
		if (basstart(iatm)==0) cycle
		do jatm=1,ncenter
			if (basstart(jatm)==0) cycle
			trdipmatatm(1,iatm,jatm)=sum( trdipmatbas(1,basstart(iatm):basend(iatm),basstart(jatm):basend(jatm)) )
			trdipmatatm(2,iatm,jatm)=sum( trdipmatbas(2,basstart(iatm):basend(iatm),basstart(jatm):basend(jatm)) )
			trdipmatatm(3,iatm,jatm)=sum( trdipmatbas(3,basstart(iatm):basend(iatm),basstart(jatm):basend(jatm)) )
		end do
	end do
	open(10,file="AAtrdipX.txt",status="replace")
	call showmatgau(trdipmatatm(1,:,:),"Atom transition dipole moment matrix (X component)",0,"f14.8",10)
	close(10)
	open(10,file="AAtrdipY.txt",status="replace")
	call showmatgau(trdipmatatm(2,:,:),"Atom transition dipole moment matrix (Y component)",0,"f14.8",10)
	close(10)
	open(10,file="AAtrdipZ.txt",status="replace")
	call showmatgau(trdipmatatm(3,:,:),"Atom transition dipole moment matrix (Z component)",0,"f14.8",10)
	close(10)
	do iatm=1,ncenter
		do jatm=1,ncenter
			trdipmatatm(1,iatm,jatm)=trdipmatatm(1,iatm,jatm)**2+trdipmatatm(2,iatm,jatm)**2+trdipmatatm(3,iatm,jatm)**2
		end do
	end do
	open(10,file="AAtrdip.txt",status="replace")
	call showmatgau(trdipmatatm(1,:,:),"Atom transition dipole moment matrix (total)",0,"f14.8",10)
	close(10)
	write(*,"(a)") " Done! X/Y/Z components of the matrix have been outputted to AAtrdipX.txt, AAtrdipY.txt and AAtrdipZ.txt, &
	&respectively. The total matrix (sum of square of X,Y,Z components) has been outputted to AAtrdip.txt in current folder"
	write(*,"(a)") " Hint: These file can be plotted as heat map via subfunction 2 of main function 18"
end if

end subroutine





!------------------------------------------------------------------
!---------- Calculate Mulliken atomic transition charges ----------
!------------------------------------------------------------------
subroutine transcharge
use defvar
use excitinfo
implicit real*8 (a-h,o-z)
real*8 :: bastrpopa(nbasis),bastrpopb(nbasis) !Transition population of each basis function

call loadallexcinfo(1)
call selexcit(istate)
call loadexccoeff(istate,1)
call genTDM(1,3)
call ask_Sbas_PBC

bastrpopa=0
if (wfntype==1.or.wfntype==4) then
	bastrpopb=0
	open(10,file="atmtrchga.chg",status="replace")
	open(11,file="atmtrchgb.chg",status="replace")
else
	open(10,file="atmtrchg.chg",status="replace")
end if
!Note that the tdmat we currently used is not a symmetry matrix, so we get average value for off-diagonal elements
do ibas=1,nbasis
	do jbas=1,nbasis
		if (ibas==jbas) then
			bastrpopa(ibas)=bastrpopa(ibas)+tdmata(ibas,jbas)*Sbas(ibas,jbas)
			if (wfntype==1.or.wfntype==4) bastrpopb(ibas)=bastrpopb(ibas)+tdmatb(ibas,jbas)*Sbas(ibas,jbas)
		else
			bastrpopa(ibas)=bastrpopa(ibas)+(tdmata(ibas,jbas)+tdmata(jbas,ibas))*Sbas(ibas,jbas)/2
			if (wfntype==1.or.wfntype==4) bastrpopb(ibas)=bastrpopb(ibas)+(tdmatb(ibas,jbas)+tdmatb(jbas,ibas))*Sbas(ibas,jbas)/2
		end if
	end do
end do
do iatm=1,ncenter
	tmpvala=0D0
	tmpvalb=0D0
    if (basstart(iatm)/=0) then
		tmpvala=-sum(bastrpopa(basstart(iatm):basend(iatm)))
        tmpvalb=-sum(bastrpopb(basstart(iatm):basend(iatm)))
    end if
	write(10,"(1x,a,4f12.6)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a,tmpvala
	if (wfntype==1.or.wfntype==4) write(11,"(1x,a,4f12.6)") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a,tmpvalb
end do

if (wfntype==0.or.wfntype==3) then
	close(10)
	write(*,"(/,a)") " Done! Mulliken atomic transition charges have been outputted to atmtrchg.chg in current folder"
else
	close(10)
	close(11)
	write(*,"(/,a)") " Done! Mulliken atomic transition charges for alpha and beta electrons have been outputted &
	&to atmtrchga.chg and atmtrchgb.chg in current folder, respectively"
end if

end subroutine





!--------------------------------------------------------------------------------
!------------- Generate natural orbitals of specific excited states -------------
!--------------------------------------------------------------------------------
!Also compatible with SF-TDDFT case
subroutine genexcitNO
use defvar
use util
use excitinfo
implicit real*8 (a-h,o-z)
character c2000tmp*2000,c200tmp*200
integer,allocatable :: excarr(:)
real*8 orbwei(nmo),orbDM(nbasis,nbasis),tmpmat(nbasis,nbasis)
logical,allocatable :: skippair(:)

if (.not.allocated(CObasa)) then
	write(*,"(a)") " Error: TDM is unable to be generated because basis function information is not available! &
	&Please carefully read corresponding sections of the manual to make clear which kinds of input files should be used!"
	write(*,*) "Press ENTER button to exit program"
	read(*,*)
	stop
end if

!CObasa and CObasb will be modified during gennatorb, so back up
if (.not.allocated(CObasa_org)) allocate(CObasa_org(nbasis,nbasis))
CObasa_org=CObasa
if ((wfntype==1.or.wfntype==4).and.(.not.allocated(CObasb_org))) then
	allocate(CObasb_org(nbasis,nbasis))
	CObasb_org=CObasb
end if
iwfntype_org=wfntype !wfntype should be backed up since subroutine "gennatorb" will alter it

call loadallexcinfo(1)
write(*,"(a)") " Input indices of the excited states, for which the natural orbitals will be generated and exported to .mwfn file, e.g. 1,4,6-8"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,nexcsel)
allocate(excarr(nexcsel))
call str2arr(c2000tmp,nexcsel,excarr)

call ask_Sbas_PBC
if (cfgcrossthres/=0) then
	write(*,"(/,a,f9.6,a)") " Note: When calculating cross term of density matrix of excited state, configurations with absolute value of coefficient <",&
    cfgcrossthres," will be ignored to reduce cost. The threshold is determined by ""cfgcrossthres"" in settings.ini"
end if

do iexc=1,nexcsel
	istate=excarr(iexc)
	write(*,"(/,' Dealing with excited state',i5)") istate
	call loadexccoeff(istate,1)
	!In order to accelerate generation speed of DM, we first get contribution of all orbitals (sum of square of its coefficients in all CSFs), &
	!then cycle orbitals and calculate DM of the orbital, and correspondingly modify the ground state DM
	orbwei=0
	do iexcitorb=1,excnorb
		imo=orbleft(iexcitorb)
		jmo=orbright(iexcitorb)
		wei=exccoeff(iexcitorb)**2
		if (excdir(iexcitorb)==1) then !Excitation
			orbwei(imo)=orbwei(imo)-wei
			orbwei(jmo)=orbwei(jmo)+wei
		else !De-excitation
			orbwei(imo)=orbwei(imo)+wei
			orbwei(jmo)=orbwei(jmo)-wei
		end if
	end do
    
    !Construct an array used in determining which configurations will be ignored during calculating the expensive cross term
    !0.01 is found to be a good compromise between cost and accuracy
	allocate(skippair(excnorb))
	skippair=.false.
	do iexcitorb=1,excnorb
		if (abs(exccoeff(iexcitorb))<cfgcrossthres) skippair(iexcitorb)=.true.
	end do
    
	if (wfntype==0.or.wfntype==3) then !Closed-shell
        write(*,*) "Calculating local term of density matrix..."
		! Old slow code:
		!do iorb=1,nmo !This loop is most time-consuming part
		!	if (orbwei(iorb)/=0) then
		!		orbDM=matmul( CObasa(:,iorb:iorb),transpose(CObasa(:,iorb:iorb)) )
		!		Ptot=Ptot+2*orbwei(iorb)*orbDM
		!	end if
		!end do
		! New fast code:
		!A temporary matrix (tmpmat) is utilized, which effectively includes orbwei(:) into the coefficient matrix, &
		!this way is significantly cheaper than the above equivalent commented code. The same trick is used in subroutine genP
        do imo=1,nmo
			tmpmat(:,imo)=2*orbwei(imo)*CObasa(:,imo)
        end do
        Ptot=Ptot+matmul_blas(tmpmat,CObasa,nbasis,nbasis,0,1)
        
        write(*,*) "Calculating cross term of density matrix..."
		!Currently only take below two cases into account (the same as hole-electron analysis):
		! Calculate occupied MO coupling |i><j|, destination virtual MOs are the same (l): add (i->l,j->l) and substract (i<-l,j<-l) to Ptot
		! Calculate virtual MO coupling  |l><m|, source occupied MOs are the same (i):     add (i->l,i->m) and substract (i<-l,i<-m) to Ptot
		do iexcitorb=1,excnorb
			call showprog(iexcitorb,excnorb)
			if (skippair(iexcitorb)) cycle
			idir=excdir(iexcitorb)
			ileft=orbleft(iexcitorb)
			iright=orbright(iexcitorb)
			do jexcitorb=1,excnorb
				if (skippair(jexcitorb)) cycle !Skip local term
				jdir=excdir(jexcitorb)
				if (idir/=jdir) cycle !Two configurations must have the same excitation direction
				jleft=orbleft(jexcitorb)
				jright=orbright(jexcitorb)
				if (ileft==jleft) then !Source MO is the same, calculate virtual MO coupling |l><m| contribution to P
					if (iright==jright) cycle
					tmpmat=exccoeff(iexcitorb)*exccoeff(jexcitorb)*matmul_blas(CObasa(:,iright:iright),CObasa(:,jright:jright),nbasis,nbasis,0,1)
					if (idir==1) then !->
						Ptot=Ptot+tmpmat
					else !<-
						Ptot=Ptot-tmpmat
					end if
				else if (iright==jright) then !Destination MO is the same, calculate occupied MO coupling |i><j| contribution to P
					tmpmat=exccoeff(iexcitorb)*exccoeff(jexcitorb)*matmul_blas(CObasa(:,ileft:ileft),CObasa(:,jleft:jleft),nbasis,nbasis,0,1)
					if (idir==1) then !->
						Ptot=Ptot-tmpmat
					else !<-
						Ptot=Ptot+tmpmat
					end if
				end if
			end do
		end do
        
        write(*,*) "Generating natural orbitals..."
		call gennatorb(1,0) !Generate spatial NOs
        
	else if (wfntype==1.or.wfntype==4) then !Open-shell
        write(*,*) "Calculating local term of density matrix..."
		!Alpha part
        do imo=1,nbasis
			tmpmat(:,imo)=orbwei(imo)*CObasa(:,imo)
        end do
        Palpha=Palpha+matmul_blas(tmpmat,CObasa,nbasis,nbasis,0,1)
		!Beta part
        do imo=1,nbasis
			tmpmat(:,imo)=orbwei(imo+nbasis)*CObasb(:,imo)
        end do
        Pbeta=Pbeta+matmul_blas(tmpmat,CObasb,nbasis,nbasis,0,1)
        
        write(*,*) "Calculating cross term of density matrix..."
		do iexcitorb=1,excnorb
			call showprog(iexcitorb,excnorb)
			if (skippair(iexcitorb)) cycle
			idir=excdir(iexcitorb)
			ileft=orbleft(iexcitorb)
			iright=orbright(iexcitorb)
			do jexcitorb=1,excnorb
				if (skippair(jexcitorb)) cycle
				jdir=excdir(jexcitorb)
				if (idir/=jdir) cycle
				jleft=orbleft(jexcitorb)
				jright=orbright(jexcitorb)
				if (ileft==jleft) then !Source MO is the same, calculate virtual MO coupling |l><m| contribution to P
					if (iright==jright) cycle !Skip local term
                    if (iright<=nbasis) then !Transition to an alpha MOs
						tmpmat=exccoeff(iexcitorb)*exccoeff(jexcitorb)*matmul_blas(CObasa(:,iright:iright),CObasa(:,jright:jright),nbasis,nbasis,0,1)
						if (idir==1) then !->
							Palpha=Palpha+tmpmat
						else !<-
							Palpha=Palpha-tmpmat
						end if
                    else !Transition to a beta MOs
						tmpmat=exccoeff(iexcitorb)*exccoeff(jexcitorb)*matmul_blas(CObasb(:,iright-nbasis:iright-nbasis),CObasb(:,jright-nbasis:jright-nbasis),nbasis,nbasis,0,1)
						if (idir==1) then !->
							Pbeta=Pbeta+tmpmat
						else !<-
							Pbeta=Pbeta-tmpmat
						end if
                    end if
				else if (iright==jright) then !Destination MO is the same, calculate occupied MO coupling |i><j| contribution to P
                    if (ileft<=nbasis) then !Transition from an alpha MO
						tmpmat=exccoeff(iexcitorb)*exccoeff(jexcitorb)*matmul_blas(CObasa(:,ileft:ileft),CObasa(:,jleft:jleft),nbasis,nbasis,0,1)
                        if (idir==1) then !->
							Palpha=Palpha-tmpmat
						else !<-
							Palpha=Palpha+tmpmat
						end if
                    else !Transition from a beta MO
						tmpmat=exccoeff(iexcitorb)*exccoeff(jexcitorb)*matmul_blas(CObasb(:,ileft-nbasis:ileft-nbasis),CObasb(:,jleft-nbasis:jleft-nbasis),nbasis,nbasis,0,1)
                        if (idir==1) then !->
							Pbeta=Pbeta-tmpmat
						else !<-
							Pbeta=Pbeta+tmpmat
						end if
                    end if
				end if
			end do
		end do
        
        write(*,*) "Generating natural orbitals of alpha and beta spins..."
		call gennatorb(2,0) !Generate alpha and beta NOs
	end if
    
	deallocate(skippair)

	write(c200tmp,"('NO_',i4.4,'.mwfn')") istate
	write(*,"(' Exporting ',a,'...')") trim(c200tmp) 
	call outmwfn(c200tmp,10,0)
	write(*,"(a,i4,a)") " Natural orbitals of excited state",istate," have been exported to "//trim(c200tmp)//" in current folder"
	
    !Restore to original status and regenerate ground state density matrix
	MOocc=MOocc_org
    MOene=MOene_org
    CObasa=CObasa_org
    if (allocated(CObasb)) CObasb=CObasb_org
    wfntype=iwfntype_org
	call genP
end do

deallocate(CObasa_org)
if (allocated(CObasb_org)) deallocate(CObasb_org)
write(*,"(/,a)") " Done! .mwfn files containing natural orbitals of selected excited state(s) have been successfully generated!"
end subroutine




!------------------------------------------------------------------
!-------- Print major MO transitions in all excited states --------
!------------------------------------------------------------------
!This routine is independent of wavefunction information. Output file of excited state calculation should be used as input file
subroutine majorMOtrans
use defvar
use excitinfo
use util
implicit real*8 (a-h,o-z)
character c2tmp*2,c80tmp*80,c200tmp*200,outstr*80,outfilename*200,selectyn
integer,allocatable :: idxorder(:)
real*8,allocatable :: exccoefftmp(:)

compout=10*compthres !Composition criterion of outputting MO transitions

!Determine "iopsh", "wfntype" and "nbasis", the latter two are needed by the routines of loading excited state information
open(10,file=filename,status="old")
call outputprog(10,iprog)
if (iprog==1) then !Gaussian
    call loclabel(10,"  Beta  occ. eigenvalues",iopsh)
    if (iopsh==0) call loclabel(10,"S**2 before",iopsh) !When pop=no is used, even for open-shell case, there is no "  Beta  occ. eigenvalues"
    if (iopsh==0) then
		wfntype=0
    else if (iopsh==1) then
		wfntype=1
    end if
    call loclabel(10,"basis functions,")
    read(10,*) nbasis
else if (iprog==2) then !ORCA
    iopsh=0
    call loclabel(10,"SPIN DOWN ORBITALS",iopsh)
    if (iopsh==0) wfntype=0
    if (iopsh==1) wfntype=1
    call loclabel(10,"NO   OCC          E(Eh)") !Count number of MOs to determine nbasis (in "miniprint" case doesn't explicitly print this)
    read(10,*)
    nbasis=0
    do while(.true.)
        read(10,*,iostat=ierror) inouse
        if (ierror/=0) exit
        nbasis=nbasis+1
    end do
else if (iprog==3.or.iprog==4) then !GAMESS-US/Firefly output file
    !The inputted .gms file already contains "nbasis" and "wfntype"
    iopsh=0
    if (wfntype==1.or.wfntype==4) iopsh=1
else if (iprog==5) then !CP2K
    call loclabel(10,"Spherical basis functions",ifound)
	if (ifound==0) then
		write(*,*) "Error: Unable to find ""Spherical basis functions"""
        write(*,*) "Press ENTER button to return"
        close(10)
		return
    end if
    read(10,"(a)") c80tmp
    call readaftersign_int(c80tmp,':',nbasis)
    iopsh=0
    call loclabel(10,"U-TDDFPT",iopsh)
else if (iprog==7) then !BDF. This part of code was contributed by Cong Wang, 2022-Dec-1
    iopsh=0
    call loclabel(10,"Beta MOs",iopsh)
    if (iopsh==0) wfntype=0
    if (iopsh==1) wfntype=1
    call loclabel(10,"Basis Functions =")
    read(10,"(a)") c200tmp
    read(c200tmp(index(c200tmp,"=")+1:len_trim(c200tmp)),*) nbasis
else
    write(*,*) "Error: Unable to determine the content of this file"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if
if (iopsh==0) write(*,*) "The reference state is closed-shell"
if (iopsh==1) write(*,*) "The reference state is open-shell"
write(*,"(' The number of basis functions:',i6)") nbasis
close(10)

!Load all excited state information
excitfilename=filename
call loadallexcinfo(0)
call loadallexccoeff(1)

!Determine HOMO and LUMO index
if (iprog==1) then !Gaussian
    open(10,file=filename,status="old")
    call loclabel(10," alpha electrons")
    read(10,"(a)") c80tmp
    read(c80tmp,*) iHOMO_A
    read(c80tmp(24:),*) iHOMO_B
    close(10)
else if (iprog==2) then !ORCA
    open(10,file=filename,status="old")
    if (iORCAsTD==0) then !Normal case
        call loclabel(10,"N(Alpha)",ifound)
        if (ifound==1) then
            read(10,"(a)") c80tmp
            itmp=index(c80tmp,':')
            read(c80tmp(itmp+1:),*) naelec
            read(10,"(a)") c80tmp
            read(c80tmp(itmp+1:),*) nbelec
        else !Determine number of electrons from orbital information
            call loclabel(10,"ORBITAL ENERGIES",ifound)
            if (ifound==1) then
                read(10,*);read(10,*)
                read(10,"(a)") c80tmp
                iop=0
                if (index(c80tmp,"SPIN UP")/=0) iop=1 !Open-shell case
                read(10,*)
                naelec=0
                do while(.true.)
                    read(10,"(a)") c80tmp
                    if (c80tmp==" ".or.index(c80tmp,"time")/=0) exit
                    read(c80tmp,*) inouse,tmpval
                    naelec=naelec+tmpval
                end do
                if (iop==1) then
                    nbelec=0
                    read(10,*);read(10,*)
                    do while(.true.)
                        read(10,"(a)") c80tmp
                        if (c80tmp==" ".or.index(c80tmp,"time")/=0) exit
                        read(c80tmp,*) inouse,tmpval
                        nbelec=nbelec+tmpval
                    end do
                else
                    naelec=naelec/2
                    nbelec=naelec
                end if
            else
                write(*,"(a)") " Unable to determine the number of alpha and beta electrons from the input file, please input them respectively, e.g. 15,14"
                read(*,*) naelec,nbelec
            end if
        end if
        iHOMO_A=nint(naelec)
        iHOMO_B=nint(nbelec)
    else if (iORCAsTD==1) then !sTDA/sTDDFT
        !The situation is tricky: In this case the printed data happens to correct, because when use this function, during loading
        !configuration coefficients the naelec and nbelec are both 0
        iHOMO_A=0
        iHOMO_B=0
        if (iopsh==1) then !This situation is more complicated. Ignore
            write(*,"(a)") " Error: This function cannot be used for sTDDFT/sTDA task when reference state is open-shell"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            return
        end if
    end if
    close(10)
else if (iprog==3.or.iprog==4) then !GAMESS-US, Firefly
    iHOMO_A=nint(naelec)
    iHOMO_B=nint(nbelec)
else if (iprog==5) then !CP2K
    open(10,file=filename,status="old")
	call loclabel(10,"Number of occupied orbitals:",ifound)
    read(10,"(a)") c80tmp
    call readaftersign_int(c80tmp,':',iHOMO_A)
    if (iopsh==1) then
		call loclabel(10,"Number of occupied orbitals:",ifound,0)
		read(10,"(a)") c80tmp
		call readaftersign_int(c80tmp,':',iHOMO_B)
    end if
else if (iprog==7) then !BDF. This part of code was contributed by Cong Wang, 2022-Dec-1
    open(10,file=filename,status="old")
    call loclabel(10,"alpha electrons :")
    read(10,"(a)") c200tmp
    read(c200tmp(index(c200tmp,":")+1:len_trim(c200tmp)),*) iHOMO_A
    read(10,"(a)") c200tmp
    read(c200tmp(index(c200tmp,":")+1:len_trim(c200tmp)),*) iHOMO_B
    close(10)
end if

iLUMO_A=iHOMO_A+1
if (iopsh==0) then
    write(*,"(a,i6)") " HOMO index:",iHOMO_A
    write(*,"(a,i6)") " LUMO index:",iLUMO_A
else
	iLUMO_B=iHOMO_B+1
    write(*,"(a,2i6)") " HOMO and LUMO indices of alpha electron:",iHOMO_A,iLUMO_A
    write(*,"(a,2i6)") " HOMO and LUMO indices of beta electron: ",iHOMO_B,iLUMO_B
end if

write(*,*)
write(*,"(a,f5.1,a,/)") " Only MO transitions with absolute contribution >=",compout," % &
&are shown below. It corresponds to 10 times of ""compthres"" parameter in settings.ini"
call path2filename(filename,outfilename)
outfilename=trim(outfilename)//"_exc.txt"
do itime=1,2 !=1: Print on screen, =2: Print to file
    if (itime==1) then
        iout=6
    else if (itime==2) then
        iout=10
        open(10,file=outfilename,status="replace")
    end if
    do istat=1,nstates
        if (allexcmulti(istat)==0) then
            c2tmp=" ?"
        else
            write(c2tmp,"(i2)") allexcmulti(istat)
        end if
        write(iout,"(' #',i4,f9.4,' eV',f10.2,' nm   f=',f9.5,'   Spin multiplicity=',a,':')") istat,allexcene(istat),1239.842D0/allexcene(istat),allexcf(istat),c2tmp
        write(iout,"(2x)",advance="no")
        npair=allexcnorb(istat)
        allocate(excdir(npair),orbleft(npair),orbright(npair),exccoeff(npair))
        excdir=allexcdir(1:npair,istat)
        orbleft=allorbleft(1:npair,istat)
        orbright=allorbright(1:npair,istat)
        exccoeff=allexccoeff(1:npair,istat)
        !Sort
        allocate(idxorder(npair),exccoefftmp(npair))
        exccoefftmp=exccoeff
        forall(i=1:npair) idxorder(i)=i
        call sortr8(exccoefftmp,"abs",idxorder)
        deallocate(exccoefftmp)
        call invarri4(idxorder)
        !Print
        ifirst=1
        do idx=1,npair
	        ipair=idxorder(idx)
            comp=exccoeff(ipair)**2
            if (iopsh==0) comp=comp*2
            if (comp*100>compout) then
                !Output string at left side
                if (ifirst==1) then
                    outstr=" H"
                    ifirst=0
                else
                    outstr=", H"
                end if
                ileft=orbleft(ipair)
                c80tmp=" "
                if (iopsh==0) then
                    if (ileft/=iHOMO_A) write(c80tmp,*) ileft-iHOMO_A
                else
                    if (ileft<=nbasis) then
                        outstr=trim(outstr)//'a'
                        if (ileft/=iHOMO_A) write(c80tmp,*) ileft-iHOMO_A
                    else
                        outstr=trim(outstr)//'b'
                        if (ileft/=iHOMO_B+nbasis) write(c80tmp,*) ileft-(iHOMO_B+nbasis)
                    end if
                end if
                outstr=trim(outstr)//trim(adjustl(c80tmp))
                !Output excitation sign
                if (excdir(ipair)==1) outstr=trim(outstr)//' -> '
                if (excdir(ipair)==2) outstr=trim(outstr)//' <- '
                !Output right side
                iright=orbright(ipair)
                outstr=trim(outstr)//" L"
                c80tmp=" "
                if (iopsh==0) then
                    if (iright/=iLUMO_A) then
                        write(c80tmp,*) iright-iLUMO_A
                        outstr=trim(outstr)//'+'//trim(adjustl(c80tmp))
                    end if
                else
                    if (iright<=nbasis) then
                        outstr=trim(outstr)//'a'
                        if (iright/=iLUMO_A) then
                            write(c80tmp,*) iright-iLUMO_A
                            outstr=trim(outstr)//'+'//trim(adjustl(c80tmp))
                        end if
                    else
                        outstr=trim(outstr)//'b'
                        if (iright/=iLUMO_B+nbasis) then
                            write(c80tmp,*) iright-(iLUMO_B+nbasis)
                            outstr=trim(outstr)//'+'//trim(adjustl(c80tmp))
                        end if
                    end if
                end if
                !Output composition
                write(c80tmp,"(f5.1)") comp*100
                if (excdir(ipair)==1) outstr=trim(outstr)//' '//trim(adjustl(c80tmp))//'%'
                if (excdir(ipair)==2) outstr=trim(outstr)//' -'//trim(adjustl(c80tmp))//'%'
                !Output string
                write(iout,"(a)",advance='no') trim(outstr)
            end if
        end do
        write(iout,*)
        deallocate(excdir,orbleft,orbright,exccoeff,idxorder)
    end do
    if (itime==1) then
        write(*,"(/,a)") " Do you want to export above information to "//trim(outfilename)//" in current folder? (y/n)"
        read(*,*) selectyn
        if (selectyn=='n'.or.selectyn=='N') exit
    else  
        close(10)
        write(*,*) "Output finished!"
    end if
end do
end subroutine




!----------------------------------------------------------------------
!-------- Calculate data for plotting charge-transfer spectrum --------
!----------------------------------------------------------------------
subroutine CTspectrum
use excitinfo
use defvar
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c2000tmp*2000,tmpdir*12
integer,allocatable :: fragnatm(:),frag(:,:) !Indices of atoms in fragments, frag(1:fragnatm(i),i) are atoms in fragment i
real*8 atmhole(ncenter),atmele(ncenter)
real*8,allocatable :: fraghole(:,:),fragele(:,:)
integer sphpot_bk,radpot_bk

write(*,"(/,a,/)") " NOTE: The charge-transfer spectrum was first proposed in &
&Carbon, 187, 78 (2022) DOI: 10.1016/j.carbon.2021.11.005, please cite this paper in your work"
write(*,*) "How many fragments to define? At least 2 fragments should be defined"
write(*,*) "Note: Input 0 can exit"
read(*,*) nfrag
if (nfrag==0) return
allocate(fragnatm(nfrag),frag(ncenter,nfrag))

do ifrag=1,nfrag
    write(*,"(' Input atom indices in fragment',i3,', e.g. 3,5-8,15-20')") ifrag
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,fragnatm(ifrag),frag(:,ifrag))
end do

!Test examples\excit\NH2_C8_NO2\NH2_C8_NO2.out
!nfrag=3
!allocate(fragnatm(nfrag),frag(ncenter,nfrag))
!c2000tmp="1,13,14"
!call str2arr(c2000tmp,fragnatm(1),frag(:,1))
!c2000tmp="2-9,15-22"
!call str2arr(c2000tmp,fragnatm(2),frag(:,2))
!c2000tmp="10-12"
!call str2arr(c2000tmp,fragnatm(3),frag(:,3))

call loadallexcinfo(0)
call loadallexccoeff(0)

write(*,*) "Summary of excited states:"
do iexc=1,nstates
	sumsqrexc=0
	sumsqrdeexc=0
	do itmp=1,allexcnorb(iexc)
		if (allexcdir(itmp,iexc)==1) sumsqrexc=sumsqrexc+allexccoeff(itmp,iexc)**2
		if (allexcdir(itmp,iexc)==2) sumsqrdeexc=sumsqrdeexc-allexccoeff(itmp,iexc)**2
	end do
	sumsqrall=sumsqrexc+sumsqrdeexc
    write(*,"(' State:',i5,'   E_exc:',f8.3,' eV   f:',f7.4,'   Norm:',f8.5,'   MO pairs:',i7)") &
    iexc,allexcene(iexc),allexcf(iexc),sumsqrall,allexcnorb(iexc)
end do

write(*,*)
write(*,*) "Select the method for calculating fragment contributions to hole and electron"
write(*,*) "1 Mulliken (very fast, but incompatible with diffuse functions)"
write(*,*) "2 Hirshfeld partition (slower but very robust)"
read(*,*) imethod
if (imethod==2) then
    write(*,"(a)") " Note: During the calculation of cross term of hole and electron, &
	&the orbital pairs whose absolute value of coefficient <0.01 will be ignored to save time"
    if (iautointgrid==1) then !Though this grid setting is quite poor, I found it works equally well as (75,434)
        radpot_bk=radpot
        sphpot_bk=sphpot
        radpot=30
        sphpot=170
        write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
    end if
end if
    
allocate(fraghole(nfrag,nstates),fragele(nfrag,nstates))
open(11,file="IFCTdata.txt",status="replace")
write(11,"(a)",advance='no') "state "
do ifrag=1,nfrag
    write(11,"(1x,'hole(',i1,')')",advance='no') ifrag
    write(11,"(1x,'ele(',i1,') ')",advance='no') ifrag
end do
do ifrag=1,nfrag
    write(11,"(1x,'redis(',i1,')')",advance='no') ifrag
end do
do ifrag=1,nfrag
    do jfrag=ifrag+1,nfrag
        write(11,"(2x,i1,'->',i1,1x)",advance='no') ifrag,jfrag
        write(11,"(2x,i1,'<-',i1,1x)",advance='no') ifrag,jfrag
    end do
end do
write(11,*)

open(12,file="IFCTmajor.txt",status="replace")
write(12,*) "Contribution of IFCT terms larger than 5% are shown below"
write(12,"(/,a)") "state    f     nm"

call walltime(iwalltime1)
write(*,*)
do istate=1,nstates
    write(*,"(' Calculating excited state',i6,'  of',i6)") istate,nstates
    call loadexccoeff(istate,0)
    
	if (imethod==1) then
        call atmcontri_holeele_Mulliken(atmhole,atmele,0) !Calculate atomic contribution to hole and electron by Mulliken method
	else if (imethod==2) then
        call atmcontri_holeele_Hirshfeld(atmhole,atmele,0) !Calculate atomic contribution to hole and electron by Hirshfeld method
    end if

	do ifrag=1,nfrag
		fraghole(ifrag,istate)=sum(atmhole(frag(1:fragnatm(ifrag),ifrag)))
		fragele(ifrag,istate)=sum(atmele(frag(1:fragnatm(ifrag),ifrag)))
		write(*,"(' Fragment ',i3,'    Hole:',f7.2,' %     Electron:',f7.2,' %')") ifrag,fraghole(ifrag,istate)*100,fragele(ifrag,istate)*100
	end do
    !Unphysical negative contributions may occur when using Mulliken partition, simply set them to zero
    if (imethod==1) then
        if (any(fraghole(:,istate)<-0.01D0).or.any(fragele(:,istate)<-0.01D0)) then
            write(*,"(a)") " Warning: Evident negative fragment contribution is detected! &
            &The result may be unreliable, please use Hirshfeld partition instead! Now all negative contributions are set to 0"
        end if
	    where (fraghole(:,istate)<0) fraghole(:,istate)=0
	    where (fragele(:,istate)<0) fragele(:,istate)=0
    end if
    
    !Print IFCT data to IFCTdata.txt
    write(11,"(i4,1x)",advance='no') istate
    do ifrag=1,nfrag
        write(11,"(f8.4)",advance='no') fraghole(ifrag,istate)
        write(11,"(f8.4)",advance='no') fragele(ifrag,istate)
    end do
    do ifrag=1,nfrag
        write(11,"(f9.4)",advance='no') fraghole(ifrag,istate)*fragele(ifrag,istate)
    end do
    write(11,"(1x)",advance='no')
    do ifrag=1,nfrag
        do jfrag=ifrag+1,nfrag
            write(11,"(f7.4)",advance='no') fraghole(ifrag,istate)*fragele(jfrag,istate)
            write(11,"(f7.4)",advance='no') fragele(ifrag,istate)*fraghole(jfrag,istate)
        end do
    end do
    write(11,*)
    
    !Print major IFCT terms (those > 5%) to IFCTmajor.txt
    write(12,"(i4,1x,f7.4,f7.1':')",advance='no') istate,allexcf(istate),eV2nm/allexcene(istate)
    do ifrag=1,nfrag
        tmp=fraghole(ifrag,istate)*fragele(ifrag,istate)
        if (tmp>0.05D0) write(12,"('  Redis(',i1,')',f5.1,' %')",advance='no') ifrag,tmp*100
    end do
    do ifrag=1,nfrag
        do jfrag=ifrag+1,nfrag
            tmp1=fraghole(ifrag,istate)*fragele(jfrag,istate)
            if (tmp1>0.05D0) write(12,"(2x,i1,'->',i1,f5.1,' %')",advance='no') ifrag,jfrag,tmp1*100
            tmp2=fragele(ifrag,istate)*fraghole(jfrag,istate)
            if (tmp2>0.05D0) write(12,"(2x,i1,'<-',i1,f5.1,' %')",advance='no') ifrag,jfrag,tmp2*100
        end do
    end do
    write(12,*)
end do
close(11)
close(12)

call walltime(iwalltime2)
write(*,"(/,' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

if (imethod==2.and.iautointgrid==1) then
    radpot=radpot_bk
    sphpot=sphpot_bk
end if

!Generate CTspectrum.txt
call inquire_dir("CT_multiple",alive)
if (alive) then
	if (isys==1) then !delete old wfntmp folder
		write(*,*) "Running: rmdir /S /Q CT_multiple"
		call system("rmdir /S /Q CT_multiple")
	else if (isys==2) then
		write(*,*) "Running: rm -rf CT_multiple"
		call system("rm -rf CT_multiple")
	end if
end if
call system("mkdir CT_multiple")
if (isys==1) tmpdir="CT_multiple\"
if (isys==2) tmpdir="CT_multiple/"

!Generate files included in CTspectrum\CTspectrum.txt
open(12,file=tmpdir//"CT_multiple.txt",status="replace")
write(12,"(a)") '"'//tmpdir//"total_spectrum.txt"" UV-Vis"
do ifrag=1,nfrag
    write(12,"(a,i1,a,i1)") '"'//tmpdir//"Redis_",ifrag,".txt"" Redistribution ",ifrag
end do
do ifrag=1,nfrag
    do jfrag=ifrag+1,nfrag
        write(12,"(a,i1,a,i1,a,i1,a,i1)") '"'//tmpdir//"ET_",ifrag,"to",jfrag,".txt"" Electron transfer ",ifrag,"->",jfrag
        write(12,"(a,i1,a,i1,a,i1,a,i1)") '"'//tmpdir//"ET_",jfrag,"to",ifrag,".txt"" Electron transfer ",ifrag,"<-",jfrag
    end do
end do
close(12)

!Generate data in CT_multiple\total_spectrum.txt
open(12,file=tmpdir//"total_spectrum.txt",status="replace")
write(12,*) nstates,1
do istate=1,nstates
    write(12,"(f10.5,f10.5)") allexcene(istate),allexcf(istate)
end do
close(12)

!Generate data in CT_multiple\Redis_xxx.txt
do ifrag=1,nfrag
    write(c80tmp,"(a,i1,a)") tmpdir//"Redis_",ifrag,".txt"
    open(12,file=trim(c80tmp),status="replace")
    write(12,*) nstates,1
    do istate=1,nstates
        write(12,"(f10.5,f10.5)") allexcene(istate),allexcf(istate)*fraghole(ifrag,istate)*fragele(ifrag,istate)
    end do
    close(12)
end do

!Generate data in CT_multiple\ET_xxx.txt
do ifrag=1,nfrag
    do jfrag=ifrag+1,nfrag
        write(c80tmp,"(a,i1,a,i1,a)") tmpdir//"ET_",ifrag,"to",jfrag,".txt"
        open(12,file=trim(c80tmp),status="replace")
        write(12,*) nstates,1
        do istate=1,nstates
            write(12,"(f10.5,f10.5)") allexcene(istate),allexcf(istate)*fraghole(ifrag,istate)*fragele(jfrag,istate)
        end do
        close(12)
        
        write(c80tmp,"(a,i1,a,i1,a)") tmpdir//"ET_",jfrag,"to",ifrag,".txt"
        open(12,file=trim(c80tmp),status="replace")
        write(12,*) nstates,1
        do istate=1,nstates
            write(12,"(f10.5,f10.5)") allexcene(istate),allexcf(istate)*fragele(ifrag,istate)*fraghole(jfrag,istate)
        end do
        close(12)
    end do
end do

write(*,*)
write(*,*) "The following files have been generated in current folder:"
write(*,"(a)") "   IFCTdata.txt: IFCT data of all excited states, including fragment contributions to hole and electron, as well as intrafragment &
&redistribution and interfragment electron transfer terms"
write(*,"(a)") "   IFCTmajor.txt: IFCT terms with contribution larger than 5% for every excited state"
write(*,"(a)") "   Files in ""CT_multiple"" subfolder: The CT_multiple.txt can be used as input file of Multiwfn to plot UV-Vis charge-transfer spectrum via main function 11"
end subroutine




!!--------- Calculate hole and electron at a point
!If a coefficient has magnitude less than skipthres (e.g. 0.01), then this configuration will be skipped
!This code uses the same calculation method as subroutine hole_electron, but removed irrelevant codes
subroutine calc_holeele(x,y,z,skipthres,hole,ele)
use defvar
use excitinfo
use functions
implicit real*8 (a-h,o-z)
real*8 x,y,z,skipthres,hole,ele,orbval(nmo)

hole=0
ele=0
elecross=0
holecross=0

call orbderv(1,1,nmo,x,y,z,orbval)
!Calculate local term of hole and electron
do iexcitorb=1,excnorb
	imo=orbleft(iexcitorb)
	jmo=orbright(iexcitorb)
	excwei=exccoeff(iexcitorb)**2
	if (excdir(iexcitorb)==1) then ! ->
		hole=hole+excwei*orbval(imo)**2
		ele=ele+excwei*orbval(jmo)**2
	else ! <-
		hole=hole-excwei*orbval(imo)**2
		ele=ele-excwei*orbval(jmo)**2
	end if
end do
!Calculate cross term of hole and electron, this part takes most computational time
do iexcitorb=1,excnorb
	!Below cases are skipped:
	!i->l,i->l and i->l,i<-l and i<-l,i<-l, since ileft==jleft.and.iright==jright
	!i->l,j->m, since ileft/=jleft.and.iright/=jright
	!i->l,i<-m and i<-l,j->l, since excdir(iexcitorb)/=excdir(jexcitorb)
	!**If i->l,i<-l should be taken into account is unsolved
	! Currently only take below cases into account:
	! Cross term of hole (do <i|j>):     i->l,j->l substract i<-l,j<-l
	! Cross term of electron (do <l|m>): i->l,i->m substract i<-l,i<-m
    if (abs(exccoeff(iexcitorb))<skipthres) cycle
	ileft=orbleft(iexcitorb)
	iright=orbright(iexcitorb)
	tmpleft=exccoeff(iexcitorb)*orbval(ileft) !Use temporary variable to save the time for locating element
	tmpright=exccoeff(iexcitorb)*orbval(iright)
	idir=excdir(iexcitorb)
	do jexcitorb=1,excnorb
        if (abs(exccoeff(jexcitorb))<skipthres) cycle
		jleft=orbleft(jexcitorb)
		jright=orbright(jexcitorb)
		jdir=excdir(jexcitorb)
		if (idir/=jdir) cycle
		if (ileft==jleft) then !do <l|m>
			if (iright==jright) cycle
			tmpval=tmpright*exccoeff(jexcitorb)*orbval(jright) !Originally virtual orbital
			if (idir==1) then !->
				elecross=elecross+tmpval
			else !<-
				elecross=elecross-tmpval
			end if
		else if (iright==jright) then !do <i|j>
			tmpval=tmpleft*exccoeff(jexcitorb)*orbval(jleft) !Originally occupied orbital
			if (idir==1) then !->
				holecross=holecross+tmpval
			else !<-
				holecross=holecross-tmpval
			end if
		end if
	end do
end do

!For closed-shell wavefunction, the weights are normalized to 0.5 (or say the orbitals are doubly occupied), so correct it
if (wfntype==0.or.wfntype==3) then
	hole=hole*2
	ele=ele*2
	holecross=holecross*2
	elecross=elecross*2
end if

hole=hole+holecross
ele=ele+elecross
end subroutine




!!----------- Calculate grid data of transition density of specific excited state
!istate is the index of the excited state. transdens is the grid data
subroutine calc_transdens(istate,transdens)
use defvar
use excitinfo
use functions
implicit real*8 (a-h,o-z)
integer istate
real*8 transdens(nx,ny,nz),orbval(nmo)

call loadexccoeff(istate,0)    
transdens=0
ifinish=0
ishowprog=1
ntmp=floor(ny*nz/100D0)
!$OMP PARALLEL DO SHARED(ifinish,ishowprog,transdens) &
!$OMP PRIVATE(i,j,k,tmpx,tmpy,tmpz,orbval,imo,jmo,iexcitorb) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			call orbderv(1,1,nmo,tmpx,tmpy,tmpz,orbval)
			do iexcitorb=1,excnorb
				imo=orbleft(iexcitorb)
				jmo=orbright(iexcitorb)
				transdens(i,j,k)=transdens(i,j,k)+exccoeff(iexcitorb)*orbval(imo)*orbval(jmo)
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
if (wfntype==0.or.wfntype==3) then !For closed-shell wavefunction, the weights are normalized to 0.5 (or say the orbitals are doubly occupied), so correct it
	transdens=transdens*2
end if
end subroutine



    
!!-------- Electron density polarization analysis based on electron excitations
subroutine denspolar_excit
use defvar
use excitinfo
use util
use GUI
implicit real*8 (a-h,o-z)
integer nchg !Number of point charges
real*8,allocatable :: chgx(:),chgy(:),chgz(:),chgval(:) !XYZ coordinate (Bohr) and value (e) of point charges
real*8,allocatable :: extpot(:,:,:),transdens(:,:,:),denspolar(:,:,:)
real*8,allocatable :: intval(:),E2(:),ck(:),ck_abs(:) !Various terms and contribution of excited states
integer,allocatable :: maplist(:)
character c80tmp*80

if (ifPBC/=0) then
	write(*,*) "Error: Periodic system is not supported!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if
if (.not.allocated(CObasa)) then
	write(*,*) "Error: The input file does not contain basis function information!"
	write(*,"(a)") " Please carefully read Section 3.21 of manual to make clear which kinds of input files could be used!"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	return
end if

write(*,*)
call menutitle("Electron density polarization analysis based on electron excitations",10,2)

write(*,*) "Input number of point charges to place. e.g. 3"
write(*,*) "If pressing ENTER button directly, only one point charge will be placed"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
    nchg=1
else
    read(c80tmp,*) nchg
end if
allocate(chgx(nchg),chgy(nchg),chgz(nchg),chgval(nchg))

do ichg=1,nchg
    write(*,"(a,i5)") " Input X,Y,Z (in Angstrom) and value (in e) of point charge",ichg
    write(*,*) "For example, 0.2,-0.13,3,0.05"
    read(*,*) chgx(ichg),chgy(ichg),chgz(ichg),chgval(ichg)
end do
!chgx(1)=0.00001900D0
!chgy(1)=-4.13786800D0
!chgz(1)=0
!chgval(1)=-0.1D0
chgx=chgx/b2a
chgy=chgy/b2a
chgz=chgz/b2a

write(*,*) "Placed point charges:"
do ichg=1,nchg
    write(*,"(i5,'   X=',f10.5,'   Y=',f10.5,'   Z=',f10.5,' Angstrom, Q=',f10.5,' e')") ichg,chgx(ichg)*b2a,chgy(ichg)*b2a,chgz(ichg)*b2a,chgval(ichg)
end do

write(*,*)
call setgridfixspc

allocate(extpot(nx,ny,nz),transdens(nx,ny,nz),denspolar(nx,ny,nz))
if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))

write(*,*)
call loadallexcinfo(1)

allocate(intval(nstates),ck(nstates),ck_abs(nstates),E2(nstates),maplist(nstates))
call gen_GTFuniq(1) !Generate unique GTFs, for faster evaluation
call calc_dvol(dvol)

call walltime(iwalltime1)

write(*,*) "Calculating grid data of external potential..."
extpot=0
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
            do ichg=1,nchg
				tmpr=dsqrt( (tmpx-chgx(ichg))**2 + (tmpy-chgy(ichg))**2 + (tmpz-chgz(ichg))**2 )
				extpot(i,j,k)=extpot(i,j,k)-chgval(ichg)/tmpr !Beware that electron carries negative charge, so it is -chgval(ichg)/tmpr rather than chgval(ichg)/tmpr
            end do
		end do
	end do
end do

!Loop excited states
denspolar=0
do istate=1,nstates
    write(*,"(a,i5,'...')") " Calculating excited state",istate
    call calc_transdens(istate,transdens)
    intval(istate)=sum(transdens*extpot)*dvol
    ck(istate)=intval(istate)/(-allexcene(istate)/au2eV)
    E2(istate)=ck(istate)**2 *(-allexcene(istate)/au2eV)
    denspolar=denspolar+2*ck(istate)*transdens
end do

call del_GTFuniq !Destory unique GTF informtaion
call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,' s',/)") iwalltime2-iwalltime1

write(*,*) "Excited state contributions:"
!write(*,*) " State #    Exc. Ene (Ha)      Int(a.u.)             c_k"
!do istate=1,nstates
!	write(*,"(i8,f16.6,2(1PE20.8))") istate,allexcene(istate)/au2eV,intval(istate),ck(istate)
!end do
write(*,*) " State #    Exc. Ene (Ha)       c_k          E2(kJ/mol)"
do istate=1,nstates
	write(*,"(i8,f16.6,2f16.8)") istate,allexcene(istate)/au2eV,ck(istate),E2(istate)*au2kJ
end do

ck_abs=abs(ck)
forall(i=1:nstates) maplist(i)=i
call sortr8(ck_abs,"abs",maplist)
write(*,*)
write(*,*) "Excited state contributions sorted by |c_k|:"
!write(*,*) " State #    Exc. Ene (Ha)      Int(a.u.)             c_k"
!do istate=nstates,1,-1
!	write(*,"(i8,f16.6,2(1PE20.8))") maplist(istate),allexcene(maplist(istate))/au2eV,intval(maplist(istate)),ck(maplist(istate))
!end do
write(*,*) " State #    Exc. Ene (Ha)       c_k          E2(kJ/mol)"
do istate=nstates,1,-1
	write(*,"(i8,f16.6,2f16.8)") maplist(istate),allexcene(maplist(istate))/au2eV,ck(maplist(istate)),E2(maplist(istate))*au2kJ
end do

write(*,*)
totalE2=sum(E2)
write(*,"(' Total E2:',f12.6,' eV',' (',f10.4,' kJ/mol)')") totalE2*au2eV,totalE2*au2kJ
write(*,"(' Integral of density polarization:',1PE18.8)") sum(denspolar)*dvol
write(*,"(' Amount of electrons polarized:',f12.8)") sum(abs(denspolar))*dvol/2

!Place point charges as Bq for visualization
deallocate(a)
ncenter=ncenter+nchg
allocate(a(ncenter))
a(1:ncenter_org)=a_org(:)
do ichg=1,nchg
    a(ncenter_org+ichg)%name="Bq"
    a(ncenter_org+ichg)%index=0
    a(ncenter_org+ichg)%charge=0
    a(ncenter_org+ichg)%x=chgx(ichg)
    a(ncenter_org+ichg)%y=chgy(ichg)
    a(ncenter_org+ichg)%z=chgz(ichg)
end do

do while(.true.)
	write(*,*)
	call menutitle("Post-processing menu",10,1)
	write(*,*) "0 Return"
    write(*,*) "1 Visualize isosurface of density polarization"
    write(*,*) "2 Visualize isosurface of transition density of an excited state"
    write(*,*) "3 Visualize isosurface of external potential"
    write(*,*) "10 Export grid data of density polarization"
    write(*,*) "11 Export grid data of transition density of an excited state"
    write(*,*) "12 Export grid data of external potential"
    read(*,*) isel
    if (isel==0) then
        deallocate(a)
        ncenter=ncenter_org
        allocate(a(ncenter))
        a=a_org
		return
    else if (isel==1) then
	 	cubmat=denspolar
	 	sur_value=1E-4
        drawisosurgui_lowlim=0
        drawisosurgui_highlim=0.003D0
        drawisosurgui_SWGSTP=5E-5
        isosur1style=1
		call drawisosurgui(1)
    else if (isel==2) then
		write(*,*) "Input the index of the excited state of interest, e.g. 8"
        read(*,*) istate
        call calc_transdens(istate,transdens)
	 	cubmat=transdens
	 	sur_value=2E-3
		drawisosurgui_lowlim=0
		drawisosurgui_highlim=0.1D0
        drawisosurgui_SWGSTP=5E-4
		isosur1style=1
		call drawisosurgui(1)
    else if (isel==3) then
	 	cubmat=extpot
	 	sur_value=0.1D0
        drawisosurgui_lowlim=-10
        drawisosurgui_highlim=10
        drawisosurgui_SWGSTP=0.01D0
        isosur1style=1
		call drawisosurgui(1)
    else if (isel==10) then
		write(*,*) "Exporting denspolar.cub, please wait..."
		open(10,file="denspolar.cub",status="replace")
		call outcube(denspolar,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done! Density polarization has been exported to denspolar.cub in current folder"
    else if (isel==11) then
		write(*,*) "Input the index of the excited state of interest, e.g. 8"
        read(*,*) istate
        call calc_transdens(istate,transdens)
		write(*,*) "Exporting transdens.cub, please wait..."
		open(10,file="transdens.cub",status="replace")
		call outcube(transdens,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done! Transition density has been exported to transdens.cub in current folder"
    else if (isel==12) then
		write(*,*) "Exporting extpot.cub, please wait..."
		open(10,file="extpot.cub",status="replace")
		call outcube(extpot,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done! External potential has been exported to extpot.cub in current folder"
    end if
end do
end subroutine