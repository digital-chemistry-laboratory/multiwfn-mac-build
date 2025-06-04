!!!-------- Summary of routines for calculating charges
!Bickelhaupt !Print Bickelhaupt charges
!MMPA !Population analysis of SCPA and Stout & Politzer
!fitESP !Calculate MK and CHELPG charges, including interface
!RESP !Calculate RESP charges, including interface
!spacecharge and spacecharge_evengrid !Calculate Hirshfeld, VDD, ADCH, CM5 charges and so on, including interface
!Hirshfeld_I and Hirshfeld_I_evengrid !Calculate Hirshfeld-I charges, including interface
!doADC !Calculate ADC type of charges based on existing atomic charges, invoked by spacecharge
!doCM5 !Calculate CM5 charges based on existing Hirshfeld charges, invoked by spacecharge
!genHirshfeld !A routine directly return Hirshfeld charge based on built-in density
!MBIS !Calculate MBIS charges, including interface
!EEM !Calculate EEM charges, including interface
!gasteiger !Gasteiger (PEOE) charges
    
    
!----------- Interface of various population analyses methods
subroutine population_main
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000,c80tmp*80
character :: MPApath*200=" "

write(*,*) "NOTE: There is a review comprehensively introducing various atomic charges:"
write(*,*) "Tian Lu, Qinxue Chen, Partial Charges, In Exploring Chemical Concepts Through Theory and Computation. &
&WILEY-VCH GmbH: Weinheim (2024); pp. 161-187. DOI: 10.1002/9783527843435.ch6"
if (ifragcontri==1) then
	write(*,*) "Population analysis function could not be used combining with self-defined fragment"
else
	do while(.true.)
		imodwfnold=imodwfn
		write(*,*)
		write(*,*) "     ============== Population analysis and atomic charges =============="
		write(*,*) "-2 Calculate interaction energy between fragments based on atomic charges"
		if (.not.allocated(frag1)) then
			write(*,*) "-1 Define fragment"
		else
			write(*,"(a,i5)") " -1 Redefine fragment, current number of atoms:",size(frag1)
		end if
		write(*,*) "0 Return"
		write(*,*) "1 Hirshfeld atomic charge"
		write(*,*) "2 Voronoi deformation density (VDD) atom population"
		!Not available because integration errors of below two methods by means of Becke integration are too large
	! 		write(*,*) "3 Integrate electron density in voronoi cell"
	! 		write(*,*) "4 Adjusted method 3 by Rousseau et al."
		if (allocated(CObasa)) then
			write(*,*) "5 Mulliken atom & basis function population analysis"
			write(*,*) "6 Lowdin atom & basis function population analysis"
			write(*,*) "7 Modified Mulliken atom population defined by Ros & Schuit (SCPA)"
			write(*,*) "8 Modified Mulliken atom population defined by Stout & Politzer"
			write(*,*) "9 Modified Mulliken atom population defined by Bickelhaupt"
		end if
		write(*,*) "10 Becke atomic charge with atomic dipole moment correction"
		write(*,*) "11 Atomic dipole corrected Hirshfeld atomic charge (ADCH) (recommended)"
		write(*,*) "12 CHELPG ESP fitting atomic charge"
		write(*,*) "13 Merz-Kollmann (MK) ESP fitting atomic charge"
		write(*,*) "14 AIM atomic charge"
		write(*,*) "15 Hirshfeld-I atomic charge"
		write(*,*) "16 CM5 atomic charge    -16 Generate 1.2*CM5 atomic charge"
		write(*,*) "17 Electronegativity Equalization Method (EEM) atomic charge"
		write(*,*) "18 Restrained ElectroStatic Potential (RESP) atomic charge"
        write(*,*) "19 Gasteiger (PEOE) charge"
        write(*,*) "20 Minimal Basis Iterative Stockholder (MBIS) charge"
		!write(*,*) "50 Generate input file of uESE code"
		read(*,*) ipopsel
		
        if (ifPBC/=0.and.(ipopsel==2.or.ipopsel==10.or.ipopsel==11.or.ipopsel==12.or.ipopsel==13.or.ipopsel==18)) then
			write(*,*) "Error: This method currently does not support periodic wavefunction!"
            write(*,*) "Press ENTER button to return"
			read(*,*)
            cycle
        end if
		if (ipopsel==0) then
			if (allocated(frag1)) deallocate(frag1)
			return
		else if (ipopsel==-2) then
			if (ifiletype/=4) then
				write(*,*) "Error: You must use .chg file for this function!"
				write(*,*) "Press ENTER button to continue"
				read(*,*)
				cycle
			end if
			call coulint_atmchg
		else if (ipopsel==-1) then
			if (allocated(frag1)) then
				write(*,*) "Atoms in current fragment:"
				write(*,"(13i6)") frag1
				write(*,"(a)") " Input 0 to keep unchanged, or redefine fragment, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
			else
				write(*,"(a)") " Input atomic indices to define the fragment. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
			end if
			read(*,"(a)") c2000tmp
			if (c2000tmp(1:1)/='0') then
				if (allocated(frag1)) deallocate(frag1)
				call str2arr(c2000tmp,nfrag1)
				allocate(frag1(nfrag1))
				call str2arr(c2000tmp,nfrag1,frag1)
				if (any(frag1>ncenter)) then
					write(*,*) "Error: Some atomic indices exceeded valid range! Please define again"
					write(*,*)
					deallocate(frag1)
                else
				    write(*,*) "Done!"
                    if (ifiletype==4) then !chg or pqr
						write(*,"(/,' Sum of atomic charges loaded from input file for the atoms in the fragment is',f16.8)") sum(a(frag1(:))%charge)
                    end if
				end if
			end if
		else if (ipopsel==1) then
			write(*,*) "Citation of Hirshfeld method: Theor. Chim. Acta. (Berl), 44, 129-138 (1977)"
            icalcmode=2 !Atomic center grids
            if (ifiletype==7) then
				icalcmode=0  !Evenly distributed grids, and utilizing loaded electron density
            else if (ifPBC/=0) then !PBC wavefunction or .cub/VASP grid data, choose the way to calculate
				write(*,*)
				write(*,*) "Choose integration algorithm"
                write(*,*) "1: Evenly distributed grids"
                write(*,*) "2: Atomic center grids"
                write(*,"(a)") " NOTE: For periodic wavefunction only representating valence electrons, 1 is the best choice. 2 is also able to deal with core electrons, but significantly slower"
                write(*,*) "If directly pressing ENTER button, 1 will be used"
                read(*,"(a)") c80tmp
                if (c80tmp==" ") then
					icalcmode=1
                else
	                read(c80tmp,*) icalcmode
                end if
            end if
            if (icalcmode==0) then
				call spacecharge_evengrid(1,2)
            else if (icalcmode==1) then
				call spacecharge_evengrid(1,1)
            else if (icalcmode==2) then
				call spacecharge(1)
            end if
		else if (ipopsel==2) then
			write(*,*) "Citation of VDD method: J. Comput. Chem., 25, 189-210 (2004)"
			call spacecharge(2)
		else if (ipopsel==3) then
			call spacecharge(3)
		else if (ipopsel==4) then
			write(*,*) "Citation: J. Mol. Struct.(Theochem), 538, 235-238 (2001)"
			call spacecharge(4)
		else if (ipopsel==5) then
            call ask_Sbas_PBC !For PBC case, calculate Sbas if it is not currently available
			do while(.true.)
				write(*,*)
				write(*,*) "              ---------- Mulliken population analysis ----------"
				if (MPApath==" ") then
					write(*,*) "-1 Choose output destination, current: Screen"
				else
					write(*,"(a)") " -1 Choose output destination, current: "//trim(MPApath)
				end if
				write(*,*) "0 Return"
				write(*,*) "1 Output Mulliken population and atomic charges"
				write(*,*) "2 Output gross atomic population matrix and decompose it"
				write(*,*) "3 Output gross basis function population matrix and decompose it"
				write(*,*) "4 Output orbital contributions to atomic populations to atmpopdcp.txt"
				read(*,*) ipopsel2
				if (ipopsel2==0) then
					exit
				else if (ipopsel2==-1) then
					write(*,*) "Input the path for printing population analysis result, e.g. C:\ACG.txt"
					write(*,*) "Note: If press ENTER button directly, result will be printed on screen"
					read(*,"(a)") MPApath
				else
					call MPA(ipopsel2,MPApath)
				end if
			end do
		else if (ipopsel==6) then
            call ask_Sbas_PBC
			write(*,*) "Performing Lowdin orthogonalization, please wait..."
			call symmortho(0)
            write(*,*)
			write(*,*) "Input the path for printing population analysis result, e.g. C:\ACG.txt"
			write(*,*) "Note: If press ENTER button directly, result will be printed on screen"
			read(*,"(a)") MPApath
			call MPA(1,MPApath)
			call dealloall(1)
			call readinfile(firstfilename,1) !Current wavefunction has been altered, recover the initial state
		else if (ipopsel==7) then
			call MMPA(1)
		else if (ipopsel==8) then
            call ask_Sbas_PBC
			call MMPA(2)
		else if (ipopsel==9) then
            call ask_Sbas_PBC
			call Bickelhaupt
		else if (ipopsel==10) then
			call spacecharge(5)
		else if (ipopsel==11) then
			write(*,"(a)") " Citation of ADCH: Tian Lu, Feiwu Chen, Atomic dipole moment corrected Hirshfeld population method, J. Theor. Comput. Chem., 11, 163 (2012)"
			write(*,*)
			call spacecharge(6)
		else if (ipopsel==12) then
			call fitESP(2)
		else if (ipopsel==13) then
			call fitESP(1)
		else if (ipopsel==14) then !QTAIM
			write(*,"(a)") " NOTE: AIM charges cannot be calculated in present module but can be calculated in basin analysis module, &
			&please check the example given in Section 4.17.1 of the manual on how to do this"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
		else if (ipopsel==15) then !Hirshfeld-I
            if (ifiletype==7) then !Evenly distributed grids, and utilizing loaded electron density
				call Hirshfeld_I_evengrid(1,2)
            else if (ifPBC/=0) then !Evenly distributed grids, based on periodic wavefunction representing valence electrons
				call Hirshfeld_I_evengrid(1,1)
			else !Isolated system
				call Hirshfeld_I_wrapper(1)
            end if
		else if (abs(ipopsel)==16) then !CM5 or 1.2*CM5
            icalcmode=2 !Atomic center grids
            if (ifiletype==7) then
				icalcmode=0  !Evenly distributed grids, and utilizing loaded electron density
            else if (ifPBC/=0) then !PBC wavefunction or .cub/VASP grid data, choose the way to calculate
				write(*,*)
				write(*,*) "Choose integration algorithm"
                write(*,*) "1: Evenly distributed grids"
                write(*,*) "2: Atomic center grids"
                write(*,"(a)") " NOTE: For periodic wavefunction only representating valence electrons, 1 is the best choice. 2 is also able to deal with core electrons, but significantly slower"
                write(*,*) "If directly pressing ENTER button, 1 will be used"
                read(*,"(a)") c80tmp
                if (c80tmp==" ") then
					icalcmode=1
                else
	                read(c80tmp,*) icalcmode
                end if
            end if
            if (icalcmode==2) then
				if (ipopsel==16) call spacecharge(7)
				if (ipopsel==-16) call spacecharge(-7)
            else if (icalcmode==1) then
				if (ipopsel==16) call spacecharge_evengrid(7,1)
				if (ipopsel==-16) call spacecharge_evengrid(-7,1)
            else if (icalcmode==0) then
				if (ipopsel==16) call spacecharge_evengrid(7,2)
				if (ipopsel==-16) call spacecharge_evengrid(-7,2)
            end if
		else if (ipopsel==17) then
			call EEM
		else if (ipopsel==18) then
			call RESP
        else if (ipopsel==19) then
            call gasteiger
        else if (ipopsel==20) then !MBIS
            if (ifiletype==7) then !Evenly distributed grids, and utilizing loaded electron density
				call MBIS_wrapper(1,2)
            else if (ifPBC/=0) then !Evenly distributed grids, based on periodic wavefunction representing valence electrons
				call MBIS_wrapper(1,1)
			else !Isolated system
				call MBIS_wrapper(1,0)
            end if
		end if
		if (imodwfnold==1.and.(ipopsel==1.or.ipopsel==2.or.ipopsel==6.or.ipopsel==11)) then !1,2,6,11 are the methods need to reload the initial wavefunction
			write(*,"(a)") " Note: The wavefunction file has been reloaded, your previous modifications on occupation number will be ignored"
		end if
	end do
end if
end subroutine


!!---------- Calculate Coulomb interaction between two fragment based on atomic charges in .chg file
subroutine coulint_atmchg
use defvar
use util
character c2000tmp*2000
do while(.true.)
	write(*,*) "Input atom list for fragment 1, e.g. 1,4,6-9"
	write(*,*) "Input 0 can exit"
	read(*,"(a)") c2000tmp
	if (c2000tmp(1:1)=='0') exit
	call str2arr(c2000tmp,nfrag1)
	if (allocated(frag1)) deallocate(frag1)
	allocate(frag1(nfrag1))
	call str2arr(c2000tmp,nfrag1,frag1)
	write(*,*) "Input atom list for fragment 2, e.g. 1,4,6-9"
	read(*,"(a)") c2000tmp
	call str2arr(c2000tmp,nfrag2)
	if (allocated(frag2)) deallocate(frag2)
	allocate(frag2(nfrag2))
	call str2arr(c2000tmp,nfrag2,frag2)
	eleint=0
	do iatmidx=1,nfrag1
		iatm=frag1(iatmidx)
		do jatmidx=1,nfrag2
			jatm=frag2(jatmidx)
			eleint=eleint+a(iatm)%charge*a(jatm)%charge/atomdist(iatm,jatm,1)
		end do
	end do
	write(*,"(' Electrostatic interaction energy:',f12.2,' kJ/mol',f12.2,' kcal/mol',/)") eleint*au2kJ,eleint*au2kcal
end do
end subroutine



!!---------- Modified Mulliken population analysis defined by Bickelhaupt
subroutine Bickelhaupt
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 spinpop(ncenter)
real*8,target :: atmeletot(ncenter),atmelea(ncenter)
real*8,pointer :: tmpmat(:,:),tmpele(:)
do itime=1,2
	if (itime==1) then
		tmpele=>atmeletot
		tmpmat=>Ptot
	else if (itime==2) then
		tmpele=>atmelea
		tmpmat=>Palpha
	end if
	tmpele=0D0
	do i=1,ncenter
		if (basstart(i)==0) cycle
		do j=basstart(i),basend(i)
			cross=0D0
			do k=1,nbasis !This method equalvalent to use diagonal element of density matrix to partition nondiagonal element of P*S matrix
				if (k/=j) then
					if (tmpmat(j,j)+tmpmat(k,k)/=0D0) then
						cross=cross+tmpmat(j,j)/(tmpmat(j,j)+tmpmat(k,k))*tmpmat(j,k)*Sbas(j,k)
					else
						cross=cross+tmpmat(j,k)*Sbas(j,k)/2D0  !Use equivalent partition, when denominator is zero
					end if
				end if
			end do
			tmpele(i)=tmpele(i)+tmpmat(j,j)+2*cross !Plus electrons localized in basis and partitioned cross term
		end do
	end do
	if (wfntype==0.or.wfntype==3) exit
end do

if (wfntype==0.or.wfntype==3) then
	do iatm=1,ncenter
		write(*,"(' Atom',i6,'(',a2,')','  Population:',f12.8,'  Atomic charge:',f12.8)") iatm,a(iatm)%name,atmeletot(iatm),a(iatm)%charge-atmeletot(iatm)
	end do
	write(*,"(' Total net charge:',f14.8)") sum(a(:)%charge)-sum(atmeletot(:))
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	write(*,*) "    Atom      Alpha pop.   Beta pop.    Spin pop.     Atomic charge"
	totbetapop=0D0
	do iatm=1,ncenter
		betapop=atmeletot(iatm)-atmelea(iatm)
		totbetapop=totbetapop+betapop
		spinpop(iatm)=atmelea(iatm)-betapop
		write(*,"(i6,'(',a2,')',4f14.8)") iatm,a(iatm)%name,atmelea(iatm),betapop,spinpop(iatm),a(iatm)%charge-atmeletot(iatm)
	end do
	write(*,"(' Total net charge:',f14.8,'      Total spin electrons:',f12.8)") sum(a(:)%charge)-sum(atmeletot(:)),sum(atmelea(:))-totbetapop
end if

!Show fragment information
if (allocated(frag1)) then
	write(*,"(/,' Fragment charge:',f14.8)") sum(a(frag1)%charge-atmeletot(frag1))
	write(*,"(' Fragment population:',f14.8)") sum(atmeletot(frag1))
	if (wfntype==1.or.wfntype==2.or.wfntype==4) write(*,"(' Fragment spin population:',f12.8)") sum(spinpop(frag1))
end if

call outatmchg(10,a(:)%charge-atmeletot(:))
end subroutine



!!---------- Modified Mulliken population analysis
! isel=1 :Defined by Ros & Schuit (SCPA)"
! isel=2 :Defined by Stout & Politzer
subroutine MMPA(isel)
use defvar
use util
implicit real*8 (a-h,o-z)
integer isel
real*8,target :: atmelea(ncenter),atmeleb(ncenter)
real*8,pointer :: tmpmat(:,:),tmpele(:)
atmelea=0D0
atmeleb=0D0
do itime=1,2 !1=Alpha part or total electron, 2=Beta part
	do imo=1,nbasis
		if (itime==1) then
			irealmo=imo
			tmpele=>atmelea
			tmpmat=>CObasa
		else if (itime==2) then
			if (wfntype==1.or.wfntype==4) then
				irealmo=imo+nbasis
				tmpmat=>CObasb
			else !RO
				irealmo=imo
				tmpmat=>CObasa
			end if
			tmpele=>atmeleb
		end if
		if (MOocc(irealmo)==0D0) cycle
		if (isel==1) allbassqr=sum(tmpmat(:,imo)**2)
		do i=1,ncenter
			if (basstart(i)==0) cycle
			atmbassqr=sum(tmpmat(basstart(i):basend(i),imo)**2)
			if (isel==1) then !SCPA
				if (wfntype==2) then !RO
					if (MOocc(irealmo)==2D0.or.(MOtype(irealmo)==1.and.itime==1)) tmpele(i)=tmpele(i)+atmbassqr/allbassqr
				else
					tmpele(i)=tmpele(i)+MOocc(irealmo)*atmbassqr/allbassqr
				end if
			else if (isel==2) then !Stout & Politzer
				cross=0D0
				do ii=basstart(i),basend(i)
					do jj=1,nbasis
						denomin=tmpmat(ii,imo)**2+tmpmat(jj,imo)**2
						if (jj/=ii.and.denomin>=1D-120) cross=cross+tmpmat(ii,imo)**2/denomin*tmpmat(ii,imo)*tmpmat(jj,imo)*Sbas(ii,jj)
					end do
				end do
				if (wfntype==2) then !RO
					if (MOocc(irealmo)==2D0.or.(MOtype(irealmo)==1.and.itime==1)) tmpele(i)=tmpele(i)+(atmbassqr+2*cross)
				else
					tmpele(i)=tmpele(i)+MOocc(irealmo)*(atmbassqr+2*cross)
				end if
			end if
		end do
	end do
	if (wfntype==0.or.wfntype==3) exit
end do
if (wfntype==0.or.wfntype==3) then
	do iatm=1,ncenter
		write(*,"(' Atom',i6,'(',a2,')','  Population:',f12.8,'  Atomic charge:',f12.8)") iatm,a(iatm)%name,atmelea(iatm),a(iatm)%charge-atmelea(iatm)
	end do
	write(*,"(' Total net charge:',f14.8)") sum(a(:)%charge)-sum(atmelea(:))
	if (allocated(frag1)) then
        write(*,"(/,' Fragment charge:',f14.8)") sum(a(frag1)%charge-atmelea(frag1))
        write(*,"(' Fragment population:',f14.8)") sum(atmelea(frag1))
    end if
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	write(*,*) "    Atom      Alpha pop.    Beta pop.     Spin pop.      Charge"
	do iatm=1,ncenter
		write(*,"(i6,'(',a2,')',4f14.8)") iatm,a(iatm)%name,atmelea(iatm),atmeleb(iatm),atmelea(iatm)-atmeleb(iatm),a(iatm)%charge-atmelea(iatm)-atmeleb(iatm)
	end do
	write(*,"(' Total net charge:',f14.8,'   Total spin electrons:',f12.8)") sum(a(:)%charge)-sum(atmelea(:))-sum(atmeleb(:)),sum(atmelea(:))-sum(atmeleb(:))
	if (allocated(frag1)) then
		write(*,"(/,' Fragment charge:',f14.8)") sum(a(frag1)%charge-atmelea(frag1)-atmeleb(frag1))
        write(*,"(' Fragment population:',f14.8)") sum(atmelea(frag1)+atmeleb(frag1))
		write(*,"(' Fragment spin population:',f14.8)") sum(atmelea(frag1)-atmeleb(frag1))
	end if
end if

if (wfntype==0.or.wfntype==3) then
    call outatmchg(10,a(:)%charge-atmelea(:))
else
    call outatmchg(10,a(:)%charge-atmelea(:)-atmeleb(:))
end if

end subroutine




!!--------- Mulliken/Lowdin population analysis & decompose to orbital contributions
! isel=1 Output Mulliken/Lowdin population
! isel=2 Output gross atomic population matrix and decompose it to orbital contributions
! isel=3 Output gross basis function population matrix and decompose it to orbital contributions
! isel=4 Decompose atomic population to orbital contributions
! MPApath: Destination of outputting. " " means output on screen, otherwise output to specific file and the path is given here
!Note: If doing Lowdin population, density matrix and overlap matrix should be transformed first before invoking this routine
subroutine MPA(isel,MPApath)
use defvar
use util
implicit real*8 (a-h,o-z)
integer isel
real*8 MOcenmat(ncenter,nbasis),groatmmat(ncenter+1,ncenter),atmele(ncenter),charge(ncenter),spinpop(ncenter)
real*8,pointer :: ptmat(:,:)
real*8,allocatable :: tmpmat(:,:),basmata(:,:),angorbpop(:,:),angorbpopa(:,:),angorbpopb(:,:)
character selectyn,corbnum*6,cOcc*12
character(len=*) MPApath

if (MPApath==" ") then
	ides=6
else if (isel/=4) then
	open(10,file=trim(MPApath),status="replace")
	ides=10
end if

if (isel==1.or.isel==2.or.isel==4) then
	allocate(tmpmat(nbasis,nbasis),basmata(nbasis,nbasis)) !basmata stores basis gross population of alpha part
	do itime=1,3 !Total, alpha, beta
		if (itime==1) then
			tmpmat=Ptot*Sbas
		else if (itime==2) then
			tmpmat=Palpha*Sbas
			basmata=tmpmat !Backup
		else if (itime==3) then
			tmpmat=Pbeta*Sbas
		end if
		!Calculate gross atomic population matrix
        groatmmat=0
		do i=1,ncenter
			if (basstart(i)==0) cycle
			do j=1,ncenter
				if (basstart(j)==0) cycle
				accum=0D0
				do ii=basstart(i),basend(i)
					do jj=basstart(j),basend(j)
						accum=accum+tmpmat(ii,jj)
					end do
				end do
				groatmmat(i,j)=accum
			end do
		end do
		do i=1,ncenter !Stored atom populations to the last row
			groatmmat(ncenter+1,i)=sum(groatmmat(1:ncenter,i))
		end do
		totelec=0D0
		
		if (isel==1) then !Contract gross atomic population matrix and output population in each basis function/shell/atom
			if (wfntype==0.or.wfntype==3) then !Notice that only perform once (itime=1)
				allocate(angorbpop(ncenter,0:5)) !Record the population number in each angular moment orbitals, up to H
				angorbpop=0D0
				write(ides,*) "Population of basis functions:"
				write(ides,"('  Basis Type    Atom    Shell   Population')")
				do ibas=1,nbasis
					write(ides,"(i6,3x,a,i5,a,i5,f13.5)") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),'('//a(bascen(ibas))%name//')',basshell(ibas),sum(tmpmat(ibas,:))
				end do
				write(ides,*)
				write(ides,*) "Population of shells of basis functions:"
				do ish=1,nshell
					shellpop=0D0
					do ibas=1,nbasis
						if (basshell(ibas)==ish) then
							iatm=bascen(ibas) !Which atom this shell attribute to
							shellpop=shellpop+sum(tmpmat(ibas,:))
						end if
					end do
					write(ides,"(' Shell',i6,' Type: ',a,'    in atom',i5,'(',a,') :',f9.5)") ish,shtype2name(shtype(ish)),iatm,a(iatm)%name,shellpop
					iangtmp=abs(shtype(ish))
					angorbpop(iatm,iangtmp)=angorbpop(iatm,iangtmp)+shellpop
				end do
				write(ides,*)
				write(ides,*) "Population of each type of angular moment orbitals:"
				do iatm=1,ncenter
					write(ides,"(' Atom',i6,'(',a2,')',' s:',f7.4,' p:',f7.4,' d:',f7.4,' f:',f7.4,' g:',f7.4,' h:',f7.4)") iatm,a(iatm)%name,angorbpop(iatm,:)
				end do
				write(ides,"(' Sum  s:',f9.4,' p:',f9.4,' d:',f9.4,' f:',f9.4,' g:',f9.4,' h:',f9.4)") &
				sum(angorbpop(:,0)),sum(angorbpop(:,1)),sum(angorbpop(:,2)),sum(angorbpop(:,3)),sum(angorbpop(:,4)),sum(angorbpop(:,5))
				write(ides,*)
				write(ides,*) "Population of atoms:"
				do iatm=1,ncenter
					charge(iatm)=a(iatm)%charge-groatmmat(ncenter+1,iatm)
					write(ides,"(' Atom',i6,'(',a2,')','    Population:',f12.8,'    Net charge:',f12.8)") iatm,a(iatm)%name,groatmmat(ncenter+1,iatm),charge(iatm)
				end do
				write(ides,"(' Total net charge:',f14.8)") sum(a(:)%charge)-sum(groatmmat(ncenter+1,:))
			else if ((wfntype==1.or.wfntype==2.or.wfntype==4).and.itime==3) then !For unrestrict wfn, at last "itime" cycle print result
				allocate(angorbpopa(ncenter,0:5),angorbpopb(ncenter,0:5))
				angorbpopa=0D0
				angorbpopb=0D0
				write(ides,*) "Population of basis functions:"
				write(ides,"('  Basis Type    Atom    Shell   Alpha pop.   Beta pop.  Total pop.   Spin pop.')")
				do ibas=1,nbasis !Note: Currently, tmpmat stores basis gross population of beta part, basmata stores alpha part
					baspopa=sum(basmata(ibas,:))
					baspopb=sum(tmpmat(ibas,:))
					write(ides,"(i6,3x,a,i5,a,i5,1x,4f12.5)") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),&
					'('//a(bascen(ibas))%name//')',basshell(ibas),baspopa,baspopb,baspopa+baspopb,baspopa-baspopb
				end do
				write(ides,*)
				write(ides,*) "Population of shells:"
				write(ides,*) "Shell  Type     Atom     Alpha pop.  Beta pop.   Total pop.  Spin pop."
				do ish=1,nshell
					shellpopa=0D0
					shellpopb=0D0
					do ibas=1,nbasis
						if (basshell(ibas)==ish) then
							iatm=bascen(ibas) !Which atom this shell attribute to
							shellpopa=shellpopa+sum(basmata(ibas,:))
							shellpopb=shellpopb+sum(tmpmat(ibas,:))
						end if
					end do
					write(ides,"(i5,5x,a,i7,'(',a,')' ,4f12.5)") ish,shtype2name(shtype(ish)),iatm,a(iatm)%name,shellpopa,shellpopb,shellpopa+shellpopb,shellpopa-shellpopb
					iangtmp=abs(shtype(ish))
					angorbpopa(iatm,iangtmp)=angorbpopa(iatm,iangtmp)+shellpopa
					angorbpopb(iatm,iangtmp)=angorbpopb(iatm,iangtmp)+shellpopb
				end do
				write(ides,*)
				write(ides,*) "Population of each type of angular moment atomic orbitals:"
				write(ides,*) "    Atom    Type   Alpha pop.   Beta pop.    Total pop.   Spin pop."
				do iatm=1,ncenter
					if (angorbpopa(iatm,0)/=0D0.or.angorbpopb(iatm,0)/=0D0) write(ides,"(i6,'(',a2,')    s',4f13.5)") &
					iatm,a(iatm)%name,angorbpopa(iatm,0),angorbpopb(iatm,0),angorbpopa(iatm,0)+angorbpopb(iatm,0),angorbpopa(iatm,0)-angorbpopb(iatm,0)
					if (angorbpopa(iatm,1)/=0D0.or.angorbpopb(iatm,1)/=0D0) write(ides,"('              p',4f13.5)") &
					angorbpopa(iatm,1),angorbpopb(iatm,1),angorbpopa(iatm,1)+angorbpopb(iatm,1),angorbpopa(iatm,1)-angorbpopb(iatm,1)
					if (angorbpopa(iatm,2)/=0D0.or.angorbpopb(iatm,2)/=0D0) write(ides,"('              d',4f13.5)") &
					angorbpopa(iatm,2),angorbpopb(iatm,2),angorbpopa(iatm,2)+angorbpopb(iatm,2),angorbpopa(iatm,2)-angorbpopb(iatm,2)
					if (angorbpopa(iatm,3)/=0D0.or.angorbpopb(iatm,3)/=0D0) write(ides,"('              f',4f13.5)") &
					angorbpopa(iatm,3),angorbpopb(iatm,3),angorbpopa(iatm,3)+angorbpopb(iatm,3),angorbpopa(iatm,3)-angorbpopb(iatm,3)
					if (angorbpopa(iatm,4)/=0D0.or.angorbpopb(iatm,4)/=0D0) write(ides,"('              g',4f13.5)") &
					angorbpopa(iatm,4),angorbpopb(iatm,4),angorbpopa(iatm,4)+angorbpopb(iatm,4),angorbpopa(iatm,4)-angorbpopb(iatm,4)
					if (angorbpopa(iatm,5)/=0D0.or.angorbpopb(iatm,5)/=0D0) write(ides,"('              h',4f13.5)") &
					angorbpopa(iatm,5),angorbpopb(iatm,5),angorbpopa(iatm,5)+angorbpopb(iatm,5),angorbpopa(iatm,5)-angorbpopb(iatm,5)
				end do
				write(ides,*)
				if (sum(angorbpopa(:,0))/=0D0.or.sum(angorbpopb(:,0))/=0D0) write(ides,"('     Total    s',4f13.5)") &
				sum(angorbpopa(:,0)),sum(angorbpopb(:,0)),sum(angorbpopa(:,0))+sum(angorbpopb(:,0)),sum(angorbpopa(:,0))-sum(angorbpopb(:,0))
				if (sum(angorbpopa(:,1))/=0D0.or.sum(angorbpopb(:,1))/=0D0) write(ides,"('              p',4f13.5)") &
				sum(angorbpopa(:,1)),sum(angorbpopb(:,1)),sum(angorbpopa(:,1))+sum(angorbpopb(:,1)),sum(angorbpopa(:,1))-sum(angorbpopb(:,1))
				if (sum(angorbpopa(:,2))/=0D0.or.sum(angorbpopb(:,2))/=0D0) write(ides,"('              d',4f13.5)") &
				sum(angorbpopa(:,2)),sum(angorbpopb(:,2)),sum(angorbpopa(:,2))+sum(angorbpopb(:,2)),sum(angorbpopa(:,2))-sum(angorbpopb(:,2))
				if (sum(angorbpopa(:,3))/=0D0.or.sum(angorbpopb(:,3))/=0D0) write(ides,"('              f',4f13.5)") &
				sum(angorbpopa(:,3)),sum(angorbpopb(:,3)),sum(angorbpopa(:,3))+sum(angorbpopb(:,3)),sum(angorbpopa(:,3))-sum(angorbpopb(:,3))
				if (sum(angorbpopa(:,4))/=0D0.or.sum(angorbpopb(:,4))/=0D0) write(ides,"('              g',4f13.5)") &
				sum(angorbpopa(:,4)),sum(angorbpopb(:,4)),sum(angorbpopa(:,4))+sum(angorbpopb(:,4)),sum(angorbpopa(:,4))-sum(angorbpopb(:,4))
				if (sum(angorbpopa(:,5))/=0D0.or.sum(angorbpopb(:,5))/=0D0) write(ides,"('              h',4f13.5)") &
				sum(angorbpopa(:,5)),sum(angorbpopb(:,5)),sum(angorbpopa(:,5))+sum(angorbpopb(:,5)),sum(angorbpopa(:,5))-sum(angorbpopb(:,5))
				write(ides,*)
				write(ides,*) "Population of atoms:"
				write(ides,*) "    Atom      Alpha pop.   Beta pop.    Spin pop.     Atomic charge"
				totspinpop=0D0
				do iatm=1,ncenter
					alphaele=atmele(iatm)
					betaele=groatmmat(ncenter+1,iatm)
					charge(iatm)=a(iatm)%charge-(alphaele+betaele)
					spinpop(iatm)=alphaele-betaele
					write(ides,"(i6,'(',a2,')',3f13.5,f16.5)") iatm,a(iatm)%name,alphaele,betaele,spinpop(iatm),charge(iatm)
					totspinpop=totspinpop+alphaele-betaele
				end do
				write(ides,"(' Total net charge:',f12.8,'      Total spin electrons:',f12.8)") sum(charge),totspinpop
			end if
			if (itime==2) atmele(:)=groatmmat(ncenter+1,:) !Store alpha occupation of each atom to a temporary array
			
		else if (isel==2) then !Output gross atomic population matrix
			if (itime==1) call showmatgau(groatmmat,"Total gross atomic population matrix",0,"f14.8",ides)
			if (itime==2) call showmatgau(groatmmat,"Alpha gross atomic population matrix",0,"f14.8",ides)
			if (itime==3) call showmatgau(groatmmat,"Beta gross atomic population matrix",0,"f14.8",ides)
			write(ides,*)
		end if
		
		if (wfntype==0.or.wfntype==3) exit !RHF or ROHF, don't continue to process alpha & beta respectively
	end do
	
	!Show fragment information or some prompts
	if (isel==1) then
		write(*,*)
		if (allocated(frag1)) then
			write(ides,"(' Fragment charge:',f14.8)") sum(charge(frag1))
            write(ides,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(charge(frag1))
			if (wfntype==1.or.wfntype==2.or.wfntype==4) write(ides,"(' Fragment spin population:',f12.8)") sum(spinpop(frag1))
		end if
	else if (isel==2) then
		write(ides,*) "The last row is the sum of corresponding column elements (atomic population)"
	end if
	if (MPApath/=" ".and.isel/=4) then
		close(10)
		write(*,"(' Done! Data have been outputted to ',a)") trim(MPApath) 
	end if
	
    if (isel==1) call outatmchg(10,charge(:))
    
	!Decompose to orbital contributions
	selectyn='n'
	if (isel==2) then
		write(*,"(a)") " Decompose gross atomic population matrix to orbital contributions and write to groatmdcp.txt in current folder? (y/n)"
		read(*,*) selectyn
	else if (isel==4) then
		selectyn='y'
	end if
	if (selectyn=='y'.or.selectyn=='Y') then
		if (isel==2) then
			open(10,file="groatmdcp.txt",status="replace")
			write(10,*) "The last row is the sum of corresponding column elements"
			write(10,*)
		else if (isel==4) then
			open(10,file="atmpopdcp.txt",status="replace")
			write(10,*) "(i,j) elements is contribution to the ith atoms from the jth orbital"
			write(10,*)
		end if
		do itime=1,2
			MOcenmat=0D0
			if (itime==1) ptmat=>CObasa
			if (itime==2) ptmat=>CObasb
			do imo=1,nbasis
				if (itime==1) irealmo=imo
				if (itime==2) irealmo=imo+nbasis
				write(corbnum,"(i6)") imo
				write(cOcc,"(f12.8)") MOocc(irealmo)
				if (MOocc(irealmo)==0D0) cycle
				!Construct gross atomic population matrix
                groatmmat=0
				do i=1,ncenter
					if (basstart(i)==0) cycle
					do j=1,ncenter
						if (basstart(j)==0) cycle
						accum=0D0
						do ii=basstart(i),basend(i)
							do jj=basstart(j),basend(j)
								accum=accum+MOocc(irealmo)*ptmat(ii,imo)*ptmat(jj,imo)*Sbas(ii,jj)
							end do
						end do
						groatmmat(i,j)=accum
					end do
				end do
				do i=1,ncenter
					groatmmat(ncenter+1,i)=sum(groatmmat(1:ncenter,i))
				end do
				if (isel==2) then
					if (wfntype==0.or.wfntype==2.or.wfntype==3) then
						call showmatgau(groatmmat,"Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
					else if (itime==1.and.(wfntype==1.or.wfntype==4)) then
						call showmatgau(groatmmat,"Alpha Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
					else if (itime==2.and.(wfntype==1.or.wfntype==4)) then
						call showmatgau(groatmmat,"Beta Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
					end if
					write(10,*)
				else if (isel==4) then
					MOcenmat(:,imo)=groatmmat(ncenter+1,:)
				end if
			end do
			if (isel==4) then
				if (wfntype==0.or.wfntype==2) then
					call showmatgau(MOcenmat(:,1:nint(naelec)),"",0,"f14.8",10)
				else if (wfntype==1) then
					if (itime==1) call showmatgau(MOcenmat(:,1:nint(naelec)),"Alpha part",0,"f14.8",10)
					if (itime==2) call showmatgau(MOcenmat(:,1:nint(nbelec)),"Beta part",0,"f14.8",10)
				else if (wfntype==3) then
					call showmatgau(MOcenmat,"",0,"f14.8",10)
				else if (wfntype==4) then
					if (itime==1) call showmatgau(MOcenmat,"Alpha part",0,"f14.8",10)
					if (itime==2) call showmatgau(MOcenmat,"Beta part",0,"f14.8",10)
				end if
				write(10,*)
			end if
			if (wfntype==0.or.wfntype==2.or.wfntype==3) exit !ROHF needn't to separate to alpha and beta
		end do
		close(10)
		if (isel==2) then
			write(*,"(a)") " Done! The matrices have been outputted to groatmdcp.txt in current folder"
		else if (isel==4) then
			write(*,"(a)") " Done! The result have been outputted to atmpopdcp.txt in current folder"
		end if
	end if
		
else if (isel==3) then
	write(ides,*) "(i,j) element corresponds to P(i,j)*S(i,j)"
	write(ides,*)
	call showmatgau(Ptot*Sbas,"Total gross basis function population matrix",0,"f14.8",ides)
	write(ides,*)
	if (wfntype==1.or.wfntype==2.or.wfntype==4) then
		call showmatgau(Palpha*Sbas,"Alpha gross basis function population matrix",0,"f14.8",ides)
		write(ides,*)
		call showmatgau(Pbeta*Sbas,"Beta gross basis function population matrix",0,"f14.8",ides)
		write(ides,*)
	end if
	if (MPApath/=" ") then
		close(10)
		write(*,"(' Done! Data have been outputted to ',a)") trim(MPApath)
		write(*,*)
	end if
	write(*,"(a)") " Decompose gross basis function population matrix to orbital contributions and write to grobasdcp.txt in current folder? (y/n)"
	read(*,*) selectyn
	if (selectyn=='y'.or.selectyn=='Y') then
		open(10,file="grobasdcp.txt",status="replace")
		write(10,*) "Note:"
		write(10,*) "(i,j) element corresponds to occ(imo)*C(i,imo)*C(j,imo)*S(i,j)"
		write(10,*) "The last row is the sum of corresponding column elements"
		write(10,*)
		allocate(tmpmat(nbasis+1,nbasis))
		do itime=1,2
			if (itime==1) ptmat=>CObasa
			if (itime==2) ptmat=>CObasb
			do imo=1,nbasis
				if (itime==1) irealmo=imo
				if (itime==2) irealmo=imo+nbasis
				write(corbnum,"(i6)") imo
				write(cOcc,"(f12.8)") MOocc(irealmo)
				if (MOocc(irealmo)==0D0) cycle
				do ii=1,nbasis
					do jj=1,nbasis
						tmpmat(ii,jj)=MOocc(irealmo)*ptmat(ii,imo)*ptmat(jj,imo)*Sbas(ii,jj)
					end do
				end do
				do j=1,nbasis
					tmpmat(nbasis+1,j)=sum(tmpmat(1:nbasis,j))
				end do
				if (wfntype==0.or.wfntype==2.or.wfntype==3) then
					call showmatgau(tmpmat,"Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
				else if (itime==1.and.(wfntype==1.or.wfntype==4)) then
					call showmatgau(tmpmat,"Alpha Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
				else if (itime==2.and.(wfntype==1.or.wfntype==4)) then
					call showmatgau(tmpmat,"Beta Orbital"//corbnum//"  Occ:"//cOcc,0,"f14.8",10)
				end if
				write(10,*)
			end do
			if (wfntype==0.or.wfntype==2.or.wfntype==3) exit
		end do
		close(10)
		write(*,*) "Done!"
	end if
end if
end subroutine




!!========== Calculate atomic charges based on space partition methods
!chgtype: 1=Hirshfeld, 2=VDD, 3=Integrate electron density in Voronoi cell
!4=Adjusted method 3 by Rousseau et al., 5= Becke with/without ADC, 6= ADCH, 7= CM5, -7= 1.2*CM5
!Note: Use "subroutine spacecharge_evengrid" for periodic systems
subroutine spacecharge(chgtype)
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
integer chgtype
real*8 molrho(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),beckeweigrid(radpot*sphpot),selfdens(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 atmdipx(ncenter),atmdipy(ncenter),atmdipz(ncenter),charge(ncenter)
real*8 :: covr_becke(0:nelesupp) !Covalent radii used for Becke population
integer :: nbeckeiter=3
character radfilename*200

if (.not.allocated(b)) then
    write(*,"(a)") " Error: Your input file does not contain wavefunction information which is needed by this function! &
    &Please carefully check Section 2.5 of Multiwfn manual to understand which kind of input files can be used. You should use e.g. .wfn/.mwfn/.molden/.fch..."
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

if (chgtype==5) then !Select atomic radii for Becke population
	covr_becke=covr_TianLu
	iraddefine=2
	do while(.true.)
		write(*,*)
        call menutitle("Becke charge",10,1)
		write(*,*) "-1 Return"
		write(*,*) "0 Start calculation of Becke charge!"
		if (iraddefine==0) write(*,*) "1 Select the definition of atomic radii, current: Custom"
		if (iraddefine==1) write(*,*) "1 Select the definition of atomic radii, current: CSD"
		if (iraddefine==2) write(*,*) "1 Select the definition of atomic radii, current: Modified CSD"
		if (iraddefine==3) write(*,*) "1 Select the definition of atomic radii, current: Pyykko"
		if (iraddefine==4) write(*,*) "1 Select the definition of atomic radii, current: Suresh"
		if (iraddefine==5) write(*,*) "1 Select the definition of atomic radii, current: Hugo"
		write(*,"(a,i2)") " 2 Set the number of iterations for defining Becke atomic space, current:",nbeckeiter
		write(*,*) "10 Read radii from external file"
		write(*,*) "11 Modify current radii by manual input"
		write(*,*) "12 Print current radii list"
		read(*,*) isel
		if (isel==-1) then
			return
		else if (isel==0) then
			exit
		else if (isel==1) then
			write(*,*) "1 Use CSD radii (Dalton Trans., 2008, 2832-2838)"
			write(*,*) "2 Use the modified version of CSD radii defined by Tian Lu (Recommended)"
			write(*,*) "3 Use Pyykko radii (Chem. Eur.-J., 15, 186-197)"
			write(*,*) "4 Use Suresh radii (J. Phys. Chem. A, 105, 5940-5944)"
			write(*,*) "5 Use Hugo radii (Chem. Phys. Lett., 480, 127-131)"
			read(*,*) iselrad
			if (iselrad==1) then
				covr_becke=covr
				iraddefine=1
			else if (iselrad==2) then
				covr_becke=covr_TianLu
				iraddefine=2
			else if (iselrad==3) then
				covr_becke=covr_pyy
				iraddefine=3
			else if (iselrad==4) then
				covr_becke=covr_Suresh
				iraddefine=4
			else if (iselrad==5) then
				covr_becke=radii_hugo
				iraddefine=5
			end if
		else if (isel==2) then
			write(*,*) "Input a number, e.g. 3"
			read(*,*) nbeckeiter
		else if (isel==10) then
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
		else if (isel==11) then
			iraddefine=0
			write(*,*) "Input element index and radius (in Angstrom), e.g. 5,0.84"
			read(*,*) indtmp,radtmp
			covr_becke(indtmp)=radtmp/b2a
			write(*,*) "Done!"
		else if (isel==12) then
			do irad=0,nelesupp
				write(*,"(' Element:',i5,'(',a,')   Radius:',f8.3,' Angstrom')") irad,ind2name(irad),covr_becke(irad)*b2a
			end do
			write(*,*)
		end if
	end do
end if

!Generate quadrature point and weighs by combination of Gauss-Chebyshev and Lebedev grids
call gen1cintgrid(gridatmorg,iradcut)

!***** 1=Hirshfeld, 2=VDD, 6=ADCH, 7=CM5, -7=1.2*CM5
if (chgtype==1.or.chgtype==2.or.chgtype==6.or.chgtype==7.or.chgtype==-7) then
    if (ifPBC==0) then
	    write(*,*) "This task requests atomic densities, please select how to obtain them"
	    write(*,*) "1 Use build-in sphericalized atomic densities in free-states (more convenient)"
	    write(*,"(a)") " 2 Provide wavefunction file of involved elements by yourself or invoke Gaussian to automatically calculate them"
	    read(*,*) iatmdensmode
	    if (iatmdensmode==2) call setpromol !In this routine reload first molecule at the end
    else
        iatmdensmode=1
        write(*,"(a)") " Note: Build-in sphericalized atomic densities in free-states will be used in the calculation"
    end if
	write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
	write(*,*) "Calculating, please wait..."
	write(*,*)
	call walltime(nwalltime1)
	do iatm=1,ncenter !Cycle each atom to calculate their charges and dipole
		call delvirorb(0) !For faster calculation, remove high-lying virtual MOs in whole system, do not affect result
        call gen_GTFuniq(1) !Generate unique GTFs, for faster evaluation in orbderv
		atmx=a(iatm)%x
		atmy=a(iatm)%y
		atmz=a(iatm)%z
		gridatm%value=gridatmorg%value !Weight in this grid point
		gridatm%x=gridatmorg%x+atmx !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+atmy
		gridatm%z=gridatmorg%z+atmz
		!Calculate molecular density first
		!$OMP parallel do shared(molrho) private(i) num_threads(nthreads)
        do i=1+iradcut*sphpot,radpot*sphpot
			molrho(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
		end do
		!$OMP end parallel do
		call del_GTFuniq !Destory unique GTF informtaion
		call delvirorb_back(0) !Restore to wavefunction before using delvirorb (if used)
		!Calculate free atomic density to obtain promolecule density
		promol=0D0
		if (iatmdensmode==1) then !Use built-in atomic densities
			do jatm=1,ncenter
				!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
                do ipt=1+iradcut*sphpot,radpot*sphpot
					tmpdens(ipt)=calcatmdens(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
				end do
				!$OMP end parallel do
				promol=promol+tmpdens
				if (jatm==iatm) selfdens=tmpdens
			end do
            !In order to calculate Hirshfeld charge for periodic systems, the promolecular density should be calculated considering PBC,
            !as already done above (calcatmdens accounts for PBC). The density of current atom should be calculated without PBC, otherwise the function to be
            !integrated will not converge to zero as the integration points go away from current atomic center, and thus it is not integratabel
            if (ifPBC>0) then
                ifPBCorg=ifPBC
                ifPBC=0
                !$OMP parallel do shared(selfdens) private(ipt) num_threads(nthreads)
                do ipt=1+iradcut*sphpot,radpot*sphpot
				    selfdens(ipt)=calcatmdens(iatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
			    end do
                !$OMP end parallel do
                ifPBC=ifPBCorg
            end if
		else if (iatmdensmode==2) then !Calculate atomic densities based on atom .wfn files
			do jatm=1,ncenter
				call dealloall(1)
				call readwfn(custommapname(jatm),1)
				!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
                do ipt=1+iradcut*sphpot,radpot*sphpot
					tmpdens(ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
				end do
				!$OMP end parallel do
				promol=promol+tmpdens
				if (jatm==iatm) selfdens=tmpdens
			end do
			call dealloall(1)
			call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
		end if
		!Now we have needed data in hand, calculate atomic charges and atomic dipole moments
		tmpcharge=0D0
		dipx=0D0
		dipy=0D0
		dipz=0D0
		if (chgtype==1.or.chgtype==6.or.chgtype==7.or.chgtype==-7) then !Hirshfeld, ADCH charge, CM5 charge, 1.2*CM5 charge
            do i=1+iradcut*sphpot,radpot*sphpot
				if (promol(i)/=0D0) then
                    !write(10,"(3f14.6)") gridatm(i)%x,gridatm(i)%y,gridatm(i)%z
                    !write(10,"(5E14.6)") selfdens(i),promol(i),molrho(i),gridatm(i)%value,tmpcharge
					tmpv=selfdens(i)/promol(i)*molrho(i)*gridatm(i)%value
					tmpcharge=tmpcharge-tmpv
					dipx=dipx-(gridatm(i)%x-atmx)*tmpv
					dipy=dipy-(gridatm(i)%y-atmy)*tmpv
					dipz=dipz-(gridatm(i)%z-atmz)*tmpv
				end if
			end do
			if (nEDFelec==0) then
				charge(iatm)=a(iatm)%charge+tmpcharge
			else !EDF is used for some atoms. Core electron density represented by EDF has been integrated, so nuclear charge should be augmented by nEDFelecatm
				charge(iatm)=a(iatm)%charge+nEDFelecatm(iatm)+tmpcharge
            end if
		else if (chgtype==2) then !VDD charge
            do i=1+iradcut*sphpot,radpot*sphpot !Cycle each grid point of iatm, if the distance between the grid point and other atom is shorter than iatm, weight=0
				vddwei=1D0
				discen2=(gridatm(i)%x-atmx)**2+(gridatm(i)%y-atmy)**2+(gridatm(i)%z-atmz)**2 !Distance between this grid and current center atom
				do jatm=1,ncenter_org !Note: Current wfn is atomic wfn, so use _org suffix
					if (jatm==iatm) cycle
					disother2=(gridatm(i)%x-a_org(jatm)%x)**2+(gridatm(i)%y-a_org(jatm)%y)**2+(gridatm(i)%z-a_org(jatm)%z)**2
					if (disother2<discen2) then
						vddwei=0D0 !Using this weight is equivalent to using Voronoi cell
						exit
					end if
				end do
				tmpv=vddwei*(molrho(i)-promol(i))*gridatm(i)%value
				tmpcharge=tmpcharge-tmpv
				dipx=dipx-(gridatm(i)%x-atmx)*tmpv
				dipy=dipy-(gridatm(i)%y-atmy)*tmpv
				dipz=dipz-(gridatm(i)%z-atmz)*tmpv
			end do
			charge(iatm)=tmpcharge
		end if
		atmdipx(iatm)=dipx
		atmdipy(iatm)=dipy
		atmdipz(iatm)=dipz
		if (chgtype==1.or.chgtype==6.or.chgtype==7.or.chgtype==-7) then
            write(*,"(' Hirshfeld charge of atom ',i5,'(',a2,')',' is',f12.8)") iatm,a_org(iatm)%name,charge(iatm)
        else if (chgtype==2) then
            write(*,"(' VDD charge of atom ',i5,'(',a2,')',' is',f12.8)") iatm,a_org(iatm)%name,charge(iatm)
        end if
	end do
	
!***** 3=Integrate electron density in Voronoi cell, 4=Adjusted method 3 by Rousseau et al
else if (chgtype==3.or.chgtype==4) then
	write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
	write(*,*) "Calculating, please wait..."
	write(*,*)
	call walltime(nwalltime1)
	if (chgtype==4) then !vdW radius From J.Mol.Stru.(Theo.) 538,235-238 is not identical to original definition
		vdwr(1)=0.68D0/b2a
		!B,C,N,O,F
		vdwr(5)=1.46D0/b2a
		vdwr(6)=1.46D0/b2a
		vdwr(7)=1.39D0/b2a
		vdwr(8)=1.35D0/b2a
		vdwr(9)=1.29D0/b2a
		!P S Cl
		vdwr(15)=1.78D0/b2a
		vdwr(16)=1.74D0/b2a
		vdwr(17)=1.69D0/b2a
	end if
	do iatm=1,ncenter
		tmpcharge=0D0
		dipx=0D0
		dipy=0D0
		dipz=0D0
		atmx=a(iatm)%x
		atmy=a(iatm)%y
		atmz=a(iatm)%z
		gridatm%value=gridatmorg%value
		gridatm%x=gridatmorg%x+atmx !Move quadrature point with center of current atom
		gridatm%y=gridatmorg%y+atmy
		gridatm%z=gridatmorg%z+atmz
		do i=1+iradcut*sphpot,radpot*sphpot
			vorwei=1D0
			discen2=(gridatm(i)%x-atmx)**2+(gridatm(i)%y-atmy)**2+(gridatm(i)%z-atmz)**2 !Distance between this grid and current center atom
			do jatm=1,ncenter !Determine the boundary of cell
				if (jatm==iatm) cycle
				disother2=(gridatm(i)%x-a(jatm)%x)**2+(gridatm(i)%y-a(jatm)%y)**2+(gridatm(i)%z-a(jatm)%z)**2
				if (chgtype==3) then
					if (disother2<discen2) then
						vorwei=0D0 !Use this weights equivalent to use voronoi cell
						exit
					end if
				else if (chgtype==4) then !Adjusted voronoi
					vdwra=vdwr(a(iatm)%index)
					vdwrb=vdwr(a(jatm)%index)
					RAB=atomdist(iatm,jatm,1)
					rhoval=(RAB**2+discen2-disother2)/2D0/RAB
					rhoa=vdwra/(vdwra+vdwrb)*RAB
					if (rhoval>rhoa) then
						vorwei=0D0
						exit
					end if
				end if
			end do
			if (vorwei/=0D0) then
				tmpv=vorwei*fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)*gridatm(i)%value
				tmpcharge=tmpcharge-tmpv
				dipx=dipx-(gridatm(i)%x-atmx)*tmpv
				dipy=dipy-(gridatm(i)%y-atmy)*tmpv
				dipz=dipz-(gridatm(i)%z-atmz)*tmpv
			end if
		end do
		charge(iatm)=tmpcharge+a(iatm)%charge
		atmdipx(iatm)=dipx
		atmdipy(iatm)=dipy
		atmdipz(iatm)=dipz
		write(*,"(' Charge of atom ',i5,'(',a2,')',' is',f12.8)") iatm,a(iatm)%name,charge(iatm)
	end do
	
!***** Becke population
else if (chgtype==5) then
	write(*,"(' Radial grids:',i5,'    Angular grids:',i5,'   Total:',i10)") radpot,sphpot,radpot*sphpot
	write(*,*) "Calculating, please wait..."
	write(*,*)
	call walltime(nwalltime1)
	do iatm=1,ncenter !Cycle each atom to calculate their charges and dipole
		gridatm%value=gridatmorg%value !Weight in this grid point
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		!$OMP parallel do shared(tmpdens) private(i) num_threads(nthreads)
	    do i=1+iradcut*sphpot,radpot*sphpot !Calc molecular density first
		    tmpdens(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
	    end do
		!$OMP end parallel do
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_becke,nbeckeiter)
		tmpcharge=0D0
		dipx=0D0
		dipy=0D0
		dipz=0D0
		do i=1+iradcut*sphpot,radpot*sphpot
			tmpv=tmpdens(i)*beckeweigrid(i)*gridatm(i)%value
			tmpcharge=tmpcharge-tmpv
			dipx=dipx-(gridatm(i)%x-a(iatm)%x)*tmpv
			dipy=dipy-(gridatm(i)%y-a(iatm)%y)*tmpv
			dipz=dipz-(gridatm(i)%z-a(iatm)%z)*tmpv
		end do
		if (nEDFelec==0) then
			charge(iatm)=a(iatm)%charge+tmpcharge
		else !EDF is used for some atoms. Core electron density represented by EDF has been integrated, so nuclear charge should be augmented by nEDFelecatm
			charge(iatm)=a(iatm)%charge+nEDFelecatm(iatm)+tmpcharge
        end if
		atmdipx(iatm)=dipx
		atmdipy(iatm)=dipy
		atmdipz(iatm)=dipz
		write(*,"(' Becke charge of atom ',i5,'(',a2,')',' is',f12.8)") iatm,a(iatm)%name,charge(iatm)
	end do
end if

write(*,"(' Summing up all charges:',f15.8)") sum(charge)
write(*,*)
xmoldip=0D0
ymoldip=0D0
zmoldip=0D0
do i=1,ncenter
	xmoldip=xmoldip+a(i)%x*charge(i)
	ymoldip=ymoldip+a(i)%y*charge(i)
	zmoldip=zmoldip+a(i)%z*charge(i)
end do
totdip=dsqrt(xmoldip**2+ymoldip**2+zmoldip**2)
write(*,"(' Total dipole moment from atomic charges:',f12.6,' a.u.')") totdip
write(*,"(' X/Y/Z of dipole moment vector:',3f12.6,' a.u.')") xmoldip,ymoldip,zmoldip

if (chgtype==5.or.chgtype==6) then
	write(*,*)
	write(*,*) "Atomic dipole moments (a.u.):"
	do iatm=1,ncenter
		write(*,"(' Atom ',i5,'(',a2,')',' in X/Y/Z:',3f11.6,' Norm:',f11.6)") &
		iatm,a(iatm)%name,atmdipx(iatm),atmdipy(iatm),atmdipz(iatm),dsqrt(atmdipx(iatm)**2+atmdipy(iatm)**2+atmdipz(iatm)**2)
	end do
	totatmdip=dsqrt(sum(atmdipx)**2+sum(atmdipy)**2+sum(atmdipz)**2)
	write(*,"(' Total atomic dipole moment:',f12.6,' a.u.')") totatmdip
	write(*,"(' X/Y/Z of total atomic dipole:',3f12.6,' a.u.')") sum(atmdipx),sum(atmdipy),sum(atmdipz)
	corrdipx=xmoldip+sum(atmdipx) !Corresponding to actual molecular dipole moment derived from molecular density
	corrdipy=ymoldip+sum(atmdipy)
	corrdipz=zmoldip+sum(atmdipz)
	realdip=dsqrt(corrdipx**2+corrdipy**2+corrdipz**2)
	if (chgtype==5) call doADC(atmdipx,atmdipy,atmdipz,charge,realdip,5) !Becke with ADC
	if (chgtype==6) call doADC(atmdipx,atmdipy,atmdipz,charge,realdip,6) !ADCH
else if (chgtype==7) then
	call doCM5(charge,0,1)
else if (chgtype==-7) then
	call doCM5(charge,1,1)
end if

!Calculate and print normalized charge
call normalize_atmchg(charge(:))
call printatmchg(charge(:))

if (allocated(frag1)) then
    write(*,"(/,' Fragment charge:',f14.8)") sum(charge(frag1))
    write(*,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(charge(frag1))
end if

call walltime(nwalltime2)
write(*,"(/,' Calculation took up',i8,' seconds wall clock time')")  nwalltime2-nwalltime1
if (chgtype==7.and.uESEinp==1) call gen_uESE_input(10,charge)
call outatmchg(10,charge(:))
end subroutine





!!========== Calculate atomic charges based on space partition methods, using evenly distributed grid, mainly for GPW periodic wavefunctions
!chgtype: 1=Hirshfeld  7= CM5, -7= 1.2*CM5
!imode: 1=Calculate actual density from periodic wavefunction    2=Actual density is directly taken from grid data in memory
!Note: Use "subroutine spacecharge" for isolated systems
subroutine spacecharge_evengrid(chgtype,imode)
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
integer chgtype,imode
real*8 tvec(3),atmrho(ncenter),atmpop(ncenter),atmpop_tmp(ncenter),charge(ncenter)

if (imode==1) then !Calculate density from periodic wavefunction
    call setgrid_for_PBC(0.2D0,1)
    if (allocated(cubmat)) deallocate(cubmat)
    allocate(cubmat(nx,ny,nz))
	call walltime(iwalltime1)
	!Because uniform grid cannot integrate well core density, so temporarily disable EDFs
    nEDFprims_org=nEDFprims
    nEDFprims=0
    call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
    write(*,*) "Calculating electron density grid data..."
    call savecubmat(1,0,1)
    call delvirorb_back(1) !Restore to previous wavefunction
    nEDFprims=nEDFprims_org
else !Directly using loaded electron density from cub/VASP grid data, and transforming grid data information to cell information
	if (all(a%charge==0)) then
		write(*,*) "Error: All nuclear charges are zero! If this file was exported by CP2K, it is a bug. You need to manually &
        &edit the file so that effective nuclear charges (column 2 since line 8) are correctly recorded, otherwise atomic charges cannot be calculated"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
	end if
    call grid2cellinfo
    !call showcellinfo
    call walltime(iwalltime1)
end if

call calc_dvol(dvol)

write(*,*) "Calculating Hirshfeld charges..."
atmpop(:)=0
ifinish=0;ishowprog=1
ntmp=floor(ny*nz/100D0)
!$OMP PARALLEL SHARED(atmpop,ifinish,ishowprog) PRIVATE(i,j,k,tmpx,tmpy,tmpz,iatm,atmrho,prorho,atmpop_tmp,ic,jc,kc,icell,jcell,kcell,tvec,atmx,atmy,atmz,dist2,tmprho) NUM_THREADS(nthreads)
atmpop_tmp(:)=0
!$OMP DO schedule(dynamic) collapse(2)
do k=1,nz
	do j=1,ny
		do i=1,nx
			if (abs(cubmat(i,j,k))<1D-11) cycle !Note that the electron density around core produced by VASP PAW calculation can be negative, so use abs()
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
            atmrho(:)=0
            call getpointcell(tmpx,tmpy,tmpz,ic,jc,kc)
            do icell=ic-PBCnx,ic+PBCnx
                do jcell=jc-PBCny,jc+PBCny
                    do kcell=kc-PBCnz,kc+PBCnz
                        call tvec_PBC(icell,jcell,kcell,tvec)
                        do iatm=1,ncenter
                            atmx=a(iatm)%x+tvec(1)
                            atmy=a(iatm)%y+tvec(2)
                            atmz=a(iatm)%z+tvec(3)
                            dist2=(atmx-tmpx)**2+(atmy-tmpy)**2+(atmz-tmpz)**2
                            if (dist2>atmrhocutsqr(a(iatm)%index)) then
                                cycle
                            else
                                tmprho=eleraddens(a(iatm)%index,dsqrt(dist2),0)
                                atmrho(iatm)=atmrho(iatm)+tmprho
                            end if
                        end do
                    end do
                end do
            end do
            prorho=sum(atmrho(:))
            if (prorho>0) atmpop_tmp(:)=atmpop_tmp(:)+atmrho(:)/prorho*cubmat(i,j,k)*dvol
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
atmpop(:)=atmpop(:)+atmpop_tmp(:)
!$OMP END CRITICAL
!$OMP END PARALLEL
if (ishowprog/=0) call showprog(100,100)

write(*,*)
do iatm=1,ncenter
    write(*,"(' Hirshfeld population of atom ',i5,'(',a2,')',' is',f14.8)") iatm,a(iatm)%name,atmpop(iatm)
end do
write(*,"(' Sum of atomic populations:',f20.10)") sum(atmpop)

!Evaluate final Hirshfeld charge. Note that EDFs were not involved in evaluating system density
charge(:)=a(:)%charge-atmpop(:)

!Normalize atomic charges. This is not feasible if only grid data is available, &
!because in this case the nelec used in "normalize_atmchg" is simply guessed by assuming system is neutral
if (imode==1) call normalize_atmchg(charge(:))

if (chgtype==1) then
	call printatmchg(charge(:))
else if (chgtype==7) then
	call doCM5(charge,0,0)
else if (chgtype==-7) then
	call doCM5(charge,1,0)
end if
!write(*,"(' Sum of all atomic charges:',f12.6)") sum(charge)

if (allocated(frag1)) then
    write(*,"(/,' Fragment charge:',f14.8)") sum(charge(frag1))
    write(*,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(charge(frag1))
end if

call walltime(iwalltime2)
write(*,"(/,' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

call outatmchg(10,charge(:))
end subroutine





!!--------- Generate input file of uESE for evaluating solvation energy
!"charge" contains CM5 charges of all atoms. ifileid is the file id that can be used.
subroutine gen_uESE_input(ifileid,charge)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 charge(ncenter)
integer ifileid
character outname*200,c200tmp*200

call path2filename(filename,c200tmp)
write(*,*)
write(*,*) "Input path for generating uESE input file, e.g. C:\ltwd.inp"
write(*,"(a)") " If press ENTER button directly, will be exported to "//trim(c200tmp)//".inp"
read(*,"(a)") outname
if (outname==" ") outname=trim(c200tmp)//".inp"

open(ifileid,file=outname,status="replace")
write(ifileid,"(a)") "Generated by Multiwfn (http://sobereva.com/multiwfn)"
write(ifileid,"(a)") " Atomic coordinates"
do iatm=1,ncenter
	write(ifileid,"(i6,3f12.6)") a(iatm)%index,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
end do
write(ifileid,*)
write(ifileid,"(a)") " Ground state charges"
do iatm=1,ncenter
	write(ifileid,"(f8.4)",advance="no") charge(iatm)
    if (mod(iatm,10)==0) write(ifileid,*)
end do
write(ifileid,*)
write(ifileid,*)
close(ifileid)
end subroutine



!!------ Calculate atomic dipole moment corrected charge based on existing atomic charge (charge) and atomic dipole moments (dipx/y/z)
!This routine is previously specific for ADCH, but can be extended to any other types of atomic charges
!The "charge" is inputted Hirshfeld charge, finally it is replaced by ADC charge 
!chgtype 5= Becke with/without ADC, 6= ADCH
subroutine doADC(dipx,dipy,dipz,charge,realdip,chgtype)
use defvar
use util
implicit real*8 (a-h,o-z)
integer chgtype
real*8 gammamat(3,3),mat(3,3),avgr(3,1),avgrr(3,3),r(3,1),dip(3,1),tmp(1,1),eigval(3),eigvecmat(3,3)
real*8 w(ncenter),chargecorr(ncenter),charge(ncenter)
real*8 dipx(ncenter),dipy(ncenter),dipz(ncenter),realdip

write(*,*)
write(*,*) "Calculating atomic dipole moment corrected charge..."
write(*,*)
chargecorr=charge

do i=1,ncenter
	if (ishowchgtrans==1) write(*,"(' Atom: ',i4,a)") i,a(i)%name !ishowchgtrans==1 means output detail of charge transferation process during atomic dipole moment correction
	!Initialize variables
	totq=0D0
	tottmpdipx=0D0
	tottmpdipy=0D0
	tottmpdipz=0D0
	avgr=0D0
	avgrr=0D0
	dip(1,1)=dipx(i)
	dip(2,1)=dipy(i)
	dip(3,1)=dipz(i)

	!Calculate weight of every atom
	do j=1,ncenter
		r(1,1)=a(j)%x-a(i)%x
		r(2,1)=a(j)%y-a(i)%y
		r(3,1)=a(j)%z-a(i)%z
		r2=r(1,1)**2+r(2,1)**2+r(3,1)**2
		distij=dsqrt(r2)
		
		!Use modified Becke weight function with vdW radii criterion
        !Charge transfer is considered only when the distance between two atoms is smaller than sum of their vdW radii
		rmaxdist=vdwr(a(i)%index)+vdwr(a(j)%index)
		tr=distij/(rmaxdist/2D0)-1 !Transform variable so that it can in 0~rmaxdist range
		tr=1.5*tr-0.5*tr**3
		tr=1.5*tr-0.5*tr**3
		w(j)=0.5*(1-(1.5*tr-0.5*tr**3))
		if (distij>rmaxdist) w(j)=0D0

		avgr=avgr+w(j)*r
		avgrr=avgrr+w(j)*matmul(r,transpose(r))
	end do

	wtot=sum(w)
	avgr=avgr/wtot !Now avgr is <r> column vector
	avgrr=avgrr/wtot !Now avgrr is <r r^T> matrix
	gammamat=avgrr-matmul(avgr,transpose(avgr))
! 	call showmatgau(gammamat,form="f12.6")
	call diagmat(gammamat,eigvecmat,eigval,500,1D-10)
	if (outmedinfo==1) write(*,"(i5,a,3f14.10)") i,a(i)%name,eigval !Test eigenvalue of gamma matrix
	
	!Idea of the treatment:
	!The gammamat is a symmetric matrix, however some of its eigenvalues may be too close to zero to get reasonable inversed matrix,
	!therefore we diagonalize it, get its eigenvalues and eigenvectors, then all following steps are operated in new local coordinate.
	!In the new local coordinate the inverse of the gammamat is imply a diagonal matrix, whose elements are inverse of eigenvalues, hence in this case
	!the eigenvalues can be manually added by a minor value so that is inverse it not extremely large (now we simply omit very small eigenvalues)
	!position or dipole moment vectors in the new local coordinate and in old (Cartesian) coordinate can be transformed via the eigvecmat
	mat=0D0
! 	tmpmin=1D-5
! 	addtmp=tmpmin*(maxval(eigval)+tmpmin)
	do ii=1,3
		if (abs(eigval(ii))>1D-5) mat(ii,ii)=1/eigval(ii) !Ignore ADC for component less than 1E-4 to avoid numerical unstablity. This treatment is more meaningful
! 		mat(ii,ii)=1D0/(eigval(ii)+addtmp) !Original ADC implementation, still has numerical unstability problem in rare case and may worse result
	end do
	
	!Use transform matrix to transform r in old coordinate to r' in new coordinate, and transform P to P'
	!r=matmul(eigvecmat,r'), so r'=matmul(eigvecmat^(-1),r), because eigvecmat is unitary matrix, r'=matmul(transepose(eigvecmat),r)
	avgr=matmul(transpose(eigvecmat),avgr)
	dip=matmul(transpose(eigvecmat),dip)
	
! 	ishowchgtrans=1
! 	if (i==10.or.i==12) then
! 		ishowchgtrans=1
! 		write(*,"(f20.10)") addtmp
! 		write(*,"(3f20.10)") mat(1,1),mat(2,2),mat(3,3)
! 		call showmatgau(eigvecmat,form="f12.6")
! 		write(*,"(3f14.10)") avgr
! 		write(*,"(3f14.10)") dip
! 	end if

	!All values need have been calculated, now calculate final result
	do j=1,ncenter
		r(1,1)=a(j)%x-a(i)%x
		r(2,1)=a(j)%y-a(i)%y
		r(3,1)=a(j)%z-a(i)%z
		r=matmul(transpose(eigvecmat),r) ! Get r(i,j) vector in new coordinate
		tmp=w(j)/wtot*matmul(matmul(transpose(r-avgr),mat) ,dip) !delta q, namely the charge which atom i gives atom j
		chargecorr(j)=chargecorr(j)+tmp(1,1) !Charge after corrected
		if (ishowchgtrans==1) write(*,"(' Give atom ',i4,a4,f15.10,'  Weight',2f15.12)") j,a(j)%name,tmp(1,1),w(j)
		totq=totq+tmp(1,1)
		tottmpdipx=tottmpdipx+(a(j)%x-a(i)%x)*tmp(1,1)
		tottmpdipy=tottmpdipy+(a(j)%y-a(i)%y)*tmp(1,1)
		tottmpdipz=tottmpdipz+(a(j)%z-a(i)%z)*tmp(1,1)
	end do
	
	if (ishowchgtrans==1) write(*,*)
end do

call menutitle("Atomic dipole moment corrected (ADC) charges",10,1)
do i=1,ncenter
	write(*,"(' Atom: ',i4,a,'  Corrected charge:',f12.6,'  Before:',f12.6)") i,a(i)%name,chargecorr(i),charge(i)
end do
write(*,"(' Summing up all corrected charges:',f12.7)") sum(chargecorr)
if (chgtype==5) write(*,"(a)") " Note: The values shown after ""Corrected charge"" are atomic dipole moment corrected Becke charges, the ones after ""Before"" are normal Becke charges"
if (chgtype==6) write(*,"(a)") " Note: The values shown after ""Corrected charge"" are ADCH charges, the ones after ""Before"" are Hirshfeld charges"
ADCdipx=sum(a%x*chargecorr)
ADCdipy=sum(a%y*chargecorr)
ADCdipz=sum(a%z*chargecorr)
ADCdip=sqrt(ADCdipx**2+ADCdipy**2+ADCdipz**2)
write(*,*)
write(*,"(' Total dipole from ADC charges (a.u.)',f12.7,'  Error:',f12.7)") ADCdip,abs(ADCdip-realdip)
write(*,"(' X/Y/Z of dipole moment from the charge (a.u.)',3f12.6)") ADCdipx,ADCdipy,ADCdipz
charge=chargecorr !Overlay charge array, then return to Hirshfeld module and output result to .chg file
end subroutine


!!--------- Calculate CM5 or 1.2*CM5 charges by correcting inputted Hirshfeld charges
!itype=0: CM5, itype=1: 1.2*CM5
!ishowdip=1: Show dipole moment calculated by atomic charges, =0: Do not
!PBC is supported, because "atomdist" returns distance to nearest mirror
subroutine doCM5(charge,itype,ishowdip)
use defvar
use util
implicit real*8 (a-h,o-z)
integer itype,ishowdip
real*8 charge(ncenter),CMcharge(ncenter),radius(118),tvec(3)
alpha=2.474D0
!As shown in CM5 paper, the covalent radii used in CM5 equation are tabulated in CRC book 91th, where they are obtained as follows:
!For Z=1~96, the radii are the average of CSD radii (For Fe, Mn, Co the low-spin is used) and Pyykko radii
!For Z=97~118, the radii are Pyykko radii
radius(1:96)=(covr(1:96)+covr_pyy(1:96))/2D0
radius(97:118)=covr_pyy(97:118)
radius=radius*b2a !Because the radii have already been converted to Bohr, so we convert them back to Angstrom

if (ishowchgtrans==1) write(*,"(/,a)") " Details of CM5 charge correction (only terms > 1E-5 will be shown):"

if (ifPBC==0) then !Isolated system
	do iatm=1,ncenter
		CMcorr=0
		iZ=a(iatm)%index
		do jatm=1,ncenter
			if (iatm==jatm) cycle
			jZ=a(jatm)%index
			Bval=exp( -alpha*(atomdist(iatm,jatm,1)*b2a-radius(iZ)-radius(jZ)) )
			call getCM5Tval(iZ,jZ,Tval)
			CMcorr=CMcorr+Tval*Bval
			if (ishowchgtrans==1.and.abs(Tval*Bval)>1D-5) write(*,"(i4,a,i4,a,'  B_term:',f10.5,'  T_term:',f10.5,'  Corr. charge:',f10.5)") &
			iatm,a(iatm)%name,jatm,a(jatm)%name,Bval,Tval,Tval*Bval
		end do
		CMcharge(iatm)=charge(iatm)+CMcorr
	end do
else !PBC case
	do iatm=1,ncenter
		iZ=a(iatm)%index
		CMcorr=0
		do jatm=1,ncenter
			if (iatm==jatm) cycle
			jZ=a(jatm)%index
            call getCM5Tval(iZ,jZ,Tval)
            corrtmp=0
            do icell=-PBCnx,PBCnx
                do jcell=-PBCny,PBCny
                    do kcell=-PBCnz,PBCnz
                        call tvec_PBC(icell,jcell,kcell,tvec)
                        xtmp2=(a(jatm)%x+tvec(1)-a(iatm)%x)**2
                        ytmp2=(a(jatm)%y+tvec(2)-a(iatm)%y)**2
                        ztmp2=(a(jatm)%z+tvec(3)-a(iatm)%z)**2
                        dist=dsqrt(xtmp2+ytmp2+ztmp2)
						Bval=exp( -alpha*(dist*b2a-radius(iZ)-radius(jZ)) )
						corrtmp=corrtmp+Tval*Bval
                    end do
                end do
            end do
            CMcorr=CMcorr+corrtmp
			if (ishowchgtrans==1.and.abs(corrtmp)>1D-5) write(*,"(i4,a,'  Corr. charge due to',i4,a,':',f10.5)") iatm,a(iatm)%name,jatm,a(jatm)%name,corrtmp
		end do
		CMcharge(iatm)=charge(iatm)+CMcorr
	end do
end if

if (itype==1) CMcharge=CMcharge*1.2D0 !Convert to 1.2*CM5 charge

write(*,*)
if (itype==0) call menutitle("CM5 charges",10,1)
if (itype==1) call menutitle("1.2*CM5 charges",10,1)
do i=1,ncenter
	if (itype==0) write(*,"(' Atom: ',i4,a,'  CM5 charge:',f12.6,'  Hirshfeld charge:',f12.6)") i,a(i)%name,CMcharge(i),charge(i)
	if (itype==1) write(*,"(' Atom: ',i4,a,'  1.2*CM5 charge:',f12.6,'  Hirshfeld charge:',f12.6)") i,a(i)%name,CMcharge(i),charge(i)
end do
if (itype==0) write(*,"(' Summing up all CM5 charges:',f15.8)") sum(CMcharge)
if (itype==1) write(*,"(' Summing up all 1.2*CM5 charges:',f15.8)") sum(CMcharge)
CM5dipx=sum(a%x*CMcharge)
CM5dipy=sum(a%y*CMcharge)
CM5dipz=sum(a%z*CMcharge)
CM5dip=sqrt(CM5dipx**2+CM5dipy**2+CM5dipz**2)
if (ishowdip==1) then
	write(*,*)
	if (itype==0) then
		write(*,"(' Total dipole moment from CM5 charges',f12.7,' a.u.')") CM5dip
		write(*,"(' X/Y/Z of dipole moment from CM5 charges',3f10.5, ' a.u.')") CM5dipx,CM5dipy,CM5dipz
	else if (itype==1) then
		write(*,"(' Total dipole moment from 1.2*CM5 charges',f12.7,' a.u.')") CM5dip
		if (itype==1) write(*,"(' X/Y/Z of dipole moment from 1.2*CM5 charges',3f10.5, ' a.u.')") CM5dipx,CM5dipy,CM5dipz
	end if
end if
charge=CMcharge
end subroutine

!---- Return Tval used in CM5 charge
subroutine getCM5Tval(iZ,jZ,Tval)
integer iZ,jZ
real*8 Dparm(118),Tval
Dparm=0D0
Dparm(1)=0.0056D0
Dparm(2)=-0.1543D0
Dparm(4)=0.0333D0
Dparm(5)=-0.1030D0
Dparm(6)=-0.0446D0
Dparm(7)=-0.1072D0
Dparm(8)=-0.0802D0
Dparm(9)=-0.0629D0
Dparm(10)=-0.1088D0
Dparm(11)=0.0184D0
Dparm(13)=-0.0726D0
Dparm(14)=-0.0790D0
Dparm(15)=-0.0756D0
Dparm(16)=-0.0565D0
Dparm(17)=-0.0444D0
Dparm(18)=-0.0767D0
Dparm(19)=0.0130D0
Dparm(31)=-0.0512D0
Dparm(32)=-0.0557D0
Dparm(33)=-0.0533D0
Dparm(34)=-0.0399D0
Dparm(35)=-0.0313D0
Dparm(36)=-0.0541D0
Dparm(37)=0.0092D0
Dparm(49)=-0.0361D0
Dparm(50)=-0.0393D0
Dparm(51)=-0.0376D0
Dparm(52)=-0.0281D0
Dparm(53)=-0.0220D0
Dparm(54)=-0.0381D0
Dparm(55)=0.0065D0
Dparm(81)=-0.0255D0
Dparm(82)=-0.0277D0
Dparm(83)=-0.0265D0
Dparm(84)=-0.0198D0
Dparm(85)=-0.0155D0
Dparm(86)=-0.0269D0
Dparm(87)=0.0046D0
Dparm(113)=-0.0179D0
Dparm(114)=-0.0195D0
Dparm(115)=-0.0187D0
Dparm(116)=-0.0140D0
Dparm(117)=-0.0110D0
Dparm(118)=-0.0189D0
if (iZ==1.and.jZ==6) then
	Tval=0.0502D0
else if (iZ==6.and.jZ==1) then
	Tval=-0.0502D0
else if (iZ==1.and.jZ==7) then
	Tval=0.1747D0
else if (iZ==7.and.jZ==1) then
	Tval=-0.1747D0
else if (iZ==1.and.jZ==8) then
	Tval=0.1671D0
else if (iZ==8.and.jZ==1) then
	Tval=-0.1671D0
else if (iZ==6.and.jZ==7) then
	Tval=0.0556D0
else if (iZ==7.and.jZ==6) then
	Tval=-0.0556D0
else if (iZ==6.and.jZ==8) then
	Tval=0.0234D0
else if (iZ==8.and.jZ==6) then
	Tval=-0.0234D0
else if (iZ==7.and.jZ==8) then
	Tval=-0.0346D0
else if (iZ==8.and.jZ==7) then
	Tval=0.0346D0
else
	Tval=Dparm(iZ)-Dparm(jZ)
end if
end subroutine











!!============================ ESP charge ============================!!
!!============================ ESP charge ============================!!
!!============================ ESP charge ============================!!
!!============================ ESP charge ============================!!
!!============================ ESP charge ============================!!


!!------------ Calculate Restrained ElectroStatic Potential (RESP) charge
subroutine RESP
use util
use defvar
use functions
implicit real*8 (a-h,o-z)
character molfilepath*200,gauoutfilepath*200,eqvconsfilepath*200,addcenfilepath*200,c200tmp*200,c2000tmp*2000
real*8 :: hyper_a=0.0005D0,hyper_a_1=0.0005D0,hyper_a_2=0.001D0,hyper_b=0.1D0 !Hyperbolic restraint parameters
integer :: ideterbond=1,igridtype=1,iradiisel=1,iESPtype=1
real*8 tmpmat(1,1)
integer :: maxRESPiter=300
real*8 :: RESPconv=0.000001D0
!Charge constraint
character chgconsfilepath*200
integer nchgcons !The number of charge constraint terms
integer,allocatable :: chgconsnatm(:) !The i index is number of atoms in charge constraint i
integer,allocatable :: chgconsatm(:,:) !chgconsatm(1:chgconsnatm(i),i) is the atom indices of charge constraint i
real*8,allocatable :: chgconsval(:) !Value of charge constraint
!Conformation information
integer nconf
character(len=200),allocatable :: conffilepath(:) !File path of each conformation
real*8,allocatable :: confweight(:)
real*8,allocatable :: fitcen(:,:,:) !x/y/z, atom index, conformer index
!About distribution of fitting points
integer :: nMKlayer=4
real*8 :: espfitvdwr(nelesupp),sclvdwlayer(100)=(/1.4D0,1.6D0,1.8D0,2D0,(0D0,i=5,100)/)
real*8 :: fitspc=0.566917796573677D0 !0.3/b2a, spacing between grid for CHELPG
real*8 :: extdis=5.29123276802099D0 !2.8/b2a, extend 2.8 Angstrom to each side for CHELPG
real*8 :: MKptdens=1.68017136515525D0 !6D0*b2a**2, 6.0 Angstrom**2 point density per for MK. Multiply by b2a**2 to convert to Bohr**2
!Arrays used in ESP fitting
integer :: maxESPpt !The maximum of number of fitting points among all conformers (used to allocate ESPpt and ESPptval)
integer,allocatable :: nESPpt(:) !The actual number of fitting points of each conformer
real*8,allocatable :: ESPptval(:,:) !ESP values, conformer index
real*8,allocatable :: ESPpt(:,:,:) !x/y/z, fitting point index, conformer index
real*8,allocatable :: fitcenvdwr(:) !vdW radius of each fitting center
integer :: MKatmlist(ncenter) !Record index of atoms used to construct MK fitting points
integer :: CHELPGatmlist(ncenter) !Record index of atoms used to construct CHELPG fitting points
real*8,allocatable :: atmchg(:) !Final result
real*8,allocatable :: atmchg_stage1(:) !Record stage 1 result of standard RESP
!eqvlist records equivalence relationship. neqvlist is the number of equivalence constraints, eqvlistlen is the number of atoms in each equivalence constraint
!If e.g. eqvlist(1:eqvlistlen(3),3) contains 5,6, that means atoms 5 and 6 should be contrainted to be equivalent.
!If eqvlistlen(i) is 1, that means no equivalence constraint is imposed to atom i
integer neqvlist
integer,allocatable :: eqvlistlen(:),eqvlist(:,:)
!Arrays used for constructing standard RESP type of constraint 
real*8 bondedmat(ncenter,ncenter) !1/0 means the two atoms are bonded / not bonded
integer H_list(5) !Temporarily record index
integer nCHlist,CHlist(ncenter) !Record index of sp3 carbons, methyl carbons, and hydrogens attached to them, they are atoms to be fitted in RESP stage 2. nCHlist is actual length
integer neqvlist_H,eqvlistlen_H(ncenter),eqvlist_H(10,ncenter) !Constraint hydrogens in -CH3, =CH2, -CH2- to be equivalent

if (.not.allocated(b).and.ifiletype/=4) then
    write(*,"(a)") " Error: Your input file does not contain wavefunction information at all, &
	&evidently the RESP charges cannot be calculated! You should use &
	&e.g. .wfn/wfx/mwfn/fch/molden... as input file, see Section 2.5 of Multiwfn manual"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

!By default, only one conformer
nconf=1
allocate(confweight(nconf),conffilepath(nconf))
confweight(1)=1
conffilepath(1)=firstfilename

ifloadconflist=0
iloadgau=0
ieqvcons=2
ichgcons=0 !=0/1 apply/disable charge constraint
naddcen=0 !The number of additional fitting centers

maincyc: do while(.true.) !Main loop

do while(.true.) !Interface loop
	write(*,*)
	write(*,*) "            ------------ Calculation of RESP charges ------------"
	if (ifloadconflist==0) write(*,*) "-1 Load list of conformer and weights from external file"
	if (ifloadconflist==1) write(*,"(a,i4,a)") "-1 Reload list of conformers from external file, current:",nconf," conformers"
	write(*,*) "0 Return"
	write(*,*) "1 Start standard two-stage RESP fitting calculation"
	write(*,*) "2 Start one-stage ESP fitting calculation with constraints"
	if (iloadgau==0) then
		if (igridtype==1) write(*,*) "3 Set method and parameters for distributing fitting points, current: MK"
		if (igridtype==2) write(*,*) "3 Set method and parameters for distributing fitting points, current: CHELPG"
	end if
	write(*,*) "4 Set hyperbolic penalty and various other running parameters"
	if (ieqvcons==0) write(*,*) "5 Set equivalence constraint in fitting, current: No constraint"
	if (ieqvcons==1) write(*,*) "5 Set equivalence constraint in fitting, current: Customized"
	if (ieqvcons==2) write(*,*) "5 Set equivalence constraint in fitting, current: H in CH2 and CH3"
	if (ichgcons==0) write(*,*) "6 Set charge constraint in fitting, current: No constraint"
	if (ichgcons==1) write(*,*) "6 Set charge constraint in fitting, current: Customized"
	if (ideterbond==1) write(*,*) "7 Set the way of determining connectivity, current: Guess from bond length"
	if (ideterbond==2) write(*,*) "7 Set the way of determining connectivity, current: Load from .mol"
	if (iloadgau==0) write(*,"(a)") " 8 Toggle if loading fitting points and ESP values &
    &from Gaussian output file of pop=MK/CHELPG task with IOp(6/33=2) during the calculation, current: No"
	if (iloadgau==1) write(*,"(a)") " 8 Toggle if loading fitting points and ESP values &
    &from Gaussian output file of pop=MK/CHELPG task with IOp(6/33=2) during the calculation, current: Yes"
    if (naddcen==0) then
        write(*,*) "9 Load additional fitting centers, current: None"
    else
        write(*,"(a,i6,a)") " 9 Load additional fitting centers, current:",naddcen," additional centers"
    end if
	if (iradiisel==1) write(*,*) "10 Choose the atomic radii used in fitting, current: Automatic"
	if (iradiisel==2) write(*,*) "10 Choose the atomic radii used in fitting, current: Scaled UFF"
	if (iradiisel==3) write(*,*) "10 Choose the atomic radii used in fitting, current: Customized"
	if (iESPtype==1) write(*,*) "11 Choose ESP type, current: Nuclear + Electronic"
	if (iESPtype==2) write(*,*) "11 Choose ESP type, current: Electronic"
	if (iESPtype==3) write(*,*) "11 Choose ESP type, current: Transition electronic"
	read(*,*) isel
	
	if (isel==-1) then
		write(*,*) "Input path of the file containing conformer list, e.g. C:\conflist.txt"
		do while(.true.)
			read(*,"(a)") c2000tmp
			inquire(file=c2000tmp,exist=alive)
			if (alive) exit
			write(*,*) "Error: Cannot find the file, input again"
		end do
		open(10,file=c2000tmp,status="old")
		nconf=totlinenum(10,1)
		write(*,"(' There are',i5,' conformers')") nconf
		deallocate(conffilepath,confweight)
		allocate(conffilepath(nconf),confweight(nconf))
		rewind(10)
		do iconf=1,nconf
			read(10,"(a)") c200tmp
			isep=index(trim(c200tmp),' ',back=.true.)
			read(c200tmp(:isep-1),"(a)") conffilepath(iconf)
			read(c200tmp(isep+1:),*) confweight(iconf)
			inquire(file=conffilepath(iconf),exist=alive)
			if (.not.alive) then
				write(*,"(/,a)") " Error: Cannot find "//trim(conffilepath(iconf))
				write(*,*) "Press ENTER button to cancel the loading"
				read(*,*)
				nconf=1
				close(10)
				cycle maincyc
			end if
		end do
		close(10)
		totwei=sum(confweight)
		write(*,"(' Sum of weights:',f12.6)") totwei
		if (abs(totwei-1)>0.001D0) write(*,*) "Warning: The sum of weights deviates from 1.0 evidently!"
		ifloadconflist=1
	else if (isel==0) then
		Return
	else if (isel==1.or.isel==2) then !Start calculation
		exit
	else if (isel==3) then
		write(*,*) "Use which kind of fitting points?"
		write(*,*) "1 MK grid"
		write(*,*) "2 CHELPG grid"
		read(*,*) igridtype
		if (igridtype==1) then
			do while(.true.)
				write(*,*)
				write(*,*) "0 Finished!"
				write(*,"(' 1 Set number of points per Angstrom^2, current:',f10.3)") MKptdens/b2a**2 !Temporary convert to Angstrom**2 for convention
				write(*,"(' 2 Set number of layers per atom, current:',i4)") nMKlayer
				write(*,"(' 3 Set the value times van der Waals radius in each layer')")
				read(*,*) isel2
				if (isel2==0) then
					exit
				else if (isel2==1) then
					write(*,*) "Input a number (in per Angstrom^2), e.g. 6.0"
					read(*,*) MKptdens
					MKptdens=MKptdens*b2a**2
				else if (isel2==2) then
					write(*,*) "Input a value, e.g. 5"
					read(*,*) nMKlayer
				else if (isel2==3) then
					write(*,*) "Current values:"
					do ilayer=1,nMKlayer
						write(*,"(' Layer',i3,':',f8.4)") ilayer,sclvdwlayer(ilayer)
					end do
					write(*,*)
					do ilayer=1,nMKlayer
						write(*,"(a,i3,', e.g. 1.5')") " Input value for layer",ilayer
						read(*,*) sclvdwlayer(ilayer)
					end do
				end if
			end do
		else if (igridtype==2) then
			do while(.true.)
				write(*,*)
				write(*,*) "0 Finished!"
				write(*,"(' 1 Set grid spacing, current:',f7.3,' Bohr (',f7.3,' Angstrom)')") fitspc,fitspc*b2a
				write(*,"(' 2 Set box extension, current:',f7.3,' Bohr (',f7.3,' Angstrom)')") extdis,extdis*b2a
				read(*,*) isel2
				if (isel2==0) then
					exit
				else if (isel2==1) then
					write(*,*) "Input a value in Bohr, e.g. 0.5"
					read(*,*) fitspc
				else if (isel2==2) then
					write(*,*) "Input a value in Bohr, e.g. 6.5"
					read(*,*) extdis
				end if
			end do
		end if
		
	else if (isel==4) then
		do while(.true.)
			write(*,*)
			write(*,*) "0 Return"
			write(*,"(' 1 Set tightness parameter (b), current:',f8.5)") hyper_b
			write(*,"(' 2 Set restraint strength (a) for one-stage fitting, current:',f8.5)") hyper_a
			write(*,"(' 3 Set restraint strength in stage 1 of standard RESP, current:',f8.5)") hyper_a_1
			write(*,"(' 4 Set restraint strength in stage 2 of standard RESP, current:',f8.5)") hyper_a_2
			write(*,"(' 5 Set maximum number of RESP iterations, current:',i4)") maxRESPiter
			write(*,"(' 6 Set convergence threshold of RESP charge, current:',f14.10)") RESPconv
			read(*,*) isel2
			if (isel2==0) then
				exit
			else if (isel2==1) then
				write(*,*) "Input a value, e.g. 0.1"
				read(*,*) hyper_b
			else if (isel2==2) then
				write(*,*) "Input a value, e.g. 0.0005"
				read(*,*) hyper_a
			else if (isel2==3) then
				write(*,*) "Input a value, e.g. 0.0005"
				read(*,*) hyper_a_1
			else if (isel2==4) then
				write(*,*) "Input a value, e.g. 0.002"
				read(*,*) hyper_a_2
			else if (isel2==5) then
				write(*,*) "Input maximum number of RESP iterations, e.g. 35"
				read(*,*) maxRESPiter
			else if (isel2==6) then
				write(*,*) "Input charge convergence threshold, e.g. 0.0001"
				read(*,*) RESPconv
			end if
			write(*,*) "Done!"
		end do
		
	else if (isel==5) then
        write(*,*)
        write(*,"(a)") " Please select options 1~3. You can also use options 10 or 11 to generate file containing &
        &equivalence constraint, which can then be utilized by option 1"
        write(*,"(a)") " Note: For standard two-stage RESP fitting, options 0 and 1 only take effect for the first stage"
        do while(.true.)
            write(*,*)
		    write(*,*) "0 No equivalence constraint will be imposed"
		    write(*,*) "1 Load equivalence constraint setting from external plain text file"
		    write(*,"(a)") " 2 Constraint H in each =CH2, -CH2-, CH3 to be equivalent in one-stage fitting"
            write(*,"(a)") " 10 Export equivalence constraint corresponding to ""H in each =CH2, -CH2-, CH3"" to eqvcons_H.txt in current folder"
            write(*,"(a)") " 11 Generate equivalence constraint according to point group of global or local geometry and write to eqvcons_PG.txt in current folder" 
            read(*,*) ieqvcons
            if (ieqvcons==0) then
                write(*,*) "OK, no equivalence constraint will be imposed"
		    else if (ieqvcons==1) then
			    write(*,*) "Input path of the plain text file, e.g. C:\eqvcons.txt"
                write(*,*) "If pressing ENTER button directly, eqvcons.txt in current folder will be loaded"
			    do while(.true.)
				    read(*,"(a)") eqvconsfilepath
                    if (eqvconsfilepath==" ") eqvconsfilepath="eqvcons.txt"
				    inquire(file=eqvconsfilepath,exist=alive)
				    if (alive) exit
				    write(*,*) "Error: Cannot find the file, input again"
			    end do
			    write(*,*) "OK, equivalence constraint will be loaded from it during calculation"
            else if (ieqvcons==2) then
                write(*,*) "OK, this equivalence constraint will be employed during calculation"
            else if (ieqvcons==10) then
                ieqvold=ieqvcons
                ieqvcons=3
                goto 20
30              write(*,"(a)") " Done! The equivalence constraint has been exported to eqvcons_H.txt in current folder"
                ieqvcons=ieqvold
            else if (ieqvcons==11) then
                if (naddcen>0) write(*,"(a)") " Note: This function cannot be used to detect point group for additional fitting centers"
                call genPGeqvcons
		    end if
            if (ieqvcons==0.or.ieqvcons==1.or.ieqvcons==2) exit
        end do
		
	else if (isel==6) then
        write(*,"(a)") " Note: For standard two-stage RESP fitting, below options only take effect for the first stage"
        write(*,*)
		write(*,*) "0 No charge constraint will be imposed" 
		write(*,*) "1 Load charge constraint setting from external plain text file"
		read(*,*) ichgcons
		if (ichgcons==1) then
			write(*,*) "Input path of the plain text file, e.g. C:\chgcons.txt"
			do while(.true.)
				read(*,"(a)") chgconsfilepath
				inquire(file=chgconsfilepath,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			write(*,*) "OK, charge constraints will be loaded from it during calculation"
		end if
		
	else if (isel==7) then
		write(*,*) "1 Guess connectivity based on atomic covalent radii and interatomic distance"
		write(*,*) "2 Load connectivity from a .mol or .mol2 file"
		read(*,*) ideterbond
		if (ideterbond==2) then
			write(*,*) "Input file path of the .mol or .mol2 file, e.g. C:\tsushima_yoshiko.mol"
			do while(.true.)
				read(*,"(a)") molfilepath
				inquire(file=molfilepath,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			write(*,*) "OK, the connectivity will be loaded from this file during calculation"
		end if
		
	else if (isel==8) then
		if (iloadgau==1) then
			iloadgau=0
		else if (iloadgau==0) then
			iloadgau=1
			write(*,"(a)") " OK, ESP fitting points with ESP values will be directly loaded from specified Gaussian output file during the calculation"
		end if
        
    else if (isel==9) then
        write(*,*) "Input path of the file containing additional fitting centers, e.g. C:\tmp.txt"
        if (naddcen>0) write(*,*) "If you press ENTER button directly, the already loaded ones will be cleaned"
        read(*,"(a)") addcenfilepath
        if (addcenfilepath==" ") then
            naddcen=0
            write(*,*) "The current additional fitting centers have been cleaned"
        else
            inquire(file=addcenfilepath,exist=alive)
            if (.not.alive) then
                write(*,*) "Error: Cannot find the file! Press ENTER button to cancel"
                read(*,*)
                cycle
            end if
	        open(10,file=addcenfilepath,status="old")
            read(10,*) naddcen
            close(10)
            write(*,"(i6,' additional centers will be loaded from this file before fitting')") naddcen
        end if
		
	else if (isel==10) then
		write(*,*) "Please choose the way of determining atomic radii:"
		write(*,"(a)") " 1 Automatic (The radii actually used will depend on the way of distributing fitting points)"
		write(*,"(a)") " 2 UFF radii scaled by 1/1.2, defined for entire periodic table"
		write(*,"(a)") " 3 Customized (The radii will be loaded from external file during fitting)"
		read(*,*) iradiisel
		
	else if (isel==11) then
		write(*,*) "Choose ESP type:"
		write(*,*) "1 Nuclear + Electronic"
		write(*,*) "2 Electronic only"
		write(*,*) "3 Transition electronic (Mainly used for deriving TrEsp charges)"
		read(*,*) iESPtype
		if (iESPtype==3) then
			hyper_a=0
			write(*,*) "Note: Restraint strength (a) for one-stage fitting has been set to zero"
			if (ieqvcons/=0) then
				ieqvcons=0
				write(*,*) "Note: Atom equivalence constraint has been removed"
			end if
		end if
	end if
end do !end interface cycle


!*** Prepare calculation ***
nfitcen=ncenter+naddcen !Actual number of fitting centers
allocate(fitcen(3,nfitcen,nconf),fitcenvdwr(nfitcen),nESPpt(nconf))
allocate(atmchg(nfitcen),atmchg_stage1(nfitcen),eqvlistlen(nfitcen))
allocate(eqvlist(500,nfitcen)) !The size is sufficiently large

!Generate fitting centers, fitting points and calculate ESP value using internal code
if (iloadgau==0) then
	!Set and check sanity of vdW radii used in fitting
	espfitvdwr=-1
	if (iradiisel==1) then !Radii defined in MK or CHELPG 
		call setESPfitvdwr(igridtype,espfitvdwr)
	else if (iradiisel==2) then !Scaled UFF radii
		call setESPfitvdwr(0,espfitvdwr)
	else if (iradiisel==3) then !Customized radii
		call setESPfitvdwr(-1,espfitvdwr)
	end if
    !Set vdW radius for each fitting center
	do iatm=1,ncenter
		fitcenvdwr(iatm)=espfitvdwr(a(iatm)%index)
	end do
    if (naddcen>0) then
        fitcenvdwr(ncenter+1:nfitcen)=0
        write(*,*) "Note: Radii of all additional fitting centers are set to zero"
    end if
	
	!Cycle conformers. nconf may be 1
	do iconf=1,nconf
		if (ifloadconflist==1) then
			write(*,"(a,i5)") " Generating fitting points and calculate ESP for conformer",iconf
			call dealloall(1)
			call readinfile(conffilepath(iconf),1)
		end if
		!Generate information of fitting centers
		do iatm=1,ncenter !Real atoms
			fitcen(1,iatm,iconf)=a(iatm)%x
			fitcen(2,iatm,iconf)=a(iatm)%y
			fitcen(3,iatm,iconf)=a(iatm)%z
		end do
        if (naddcen>0) then !Additional fitting centers
            write(*,*) "Loading additional fitting centers"
            open(10,file=addcenfilepath,status="old")
            read(10,*)
            if (iconf>1) then
                nskip=(naddcen+1)*(iconf-1)
                do iskip=1,nskip
                    read(10,*)
                end do
            end if
            write(*,*) "Coordinate of loaded additional fitting centers (Angstrom):"
            do iaddcen=1,naddcen
                read(10,*) xtmp,ytmp,ztmp
                write(*,"(i6,3f12.6)") iaddcen,xtmp,ytmp,ztmp
                fitcen(1,ncenter+iaddcen,iconf)=xtmp/b2a
                fitcen(2,ncenter+iaddcen,iconf)=ytmp/b2a
                fitcen(3,ncenter+iaddcen,iconf)=ztmp/b2a
            end do
            close(10)
        end if
		
		if (igridtype==1) then !MK
			nMKatoms=ncenter !All atoms are used to construct MK fitting points
			forall(i=1:ncenter) MKatmlist(i)=i
			call setMKpt(1,nfitcen,fitcen(:,:,iconf),fitcenvdwr,nptthis,tmpmat,sclvdwlayer,MKptdens,nMKlayer,nMKatoms,MKatmlist) !The returned number of points is upper limit
			if (iconf==1) then !At the first conformer, allocate large enough size for ESPptval and ESPpt
				maxESPpt=nptthis*2
				allocate(ESPptval(maxESPpt,nconf),ESPpt(3,maxESPpt,nconf))
			end if
			call setMKpt(2,nfitcen,fitcen(:,:,iconf),fitcenvdwr,nptthis,ESPpt(:,:,iconf),sclvdwlayer,MKptdens,nMKlayer,nMKatoms,MKatmlist)
		else if (igridtype==2) then !CHELPG
			nCHELPGatoms=ncenter !All atoms are used to construct CHELPG fitting points
			forall(i=1:ncenter) CHELPGatmlist(i)=i
			call setCHELPGpt(1,nfitcen,fitcen(:,:,iconf),fitcenvdwr,nptthis,tmpmat,extdis,fitspc,nCHELPGatoms,CHELPGatmlist) !Return actual number of points
			if (iconf==1) then !At the first conformer, allocate large enough size for ESPptval and ESPpt
				maxESPpt=nptthis*2
				allocate(ESPptval(maxESPpt,nconf),ESPpt(3,maxESPpt,nconf))
			end if
			call setCHELPGpt(2,nfitcen,fitcen(:,:,iconf),fitcenvdwr,nptthis,ESPpt(:,:,iconf),extdis,fitspc,nCHELPGatoms,CHELPGatmlist)
		end if
		nESPpt(iconf)=nptthis

		!Calculate ESP value at fitting points
		if (nconf==1) then !Show prompts
			call fitESP_calcESP(1,iESPtype,nptthis,ESPpt(:,1:nptthis,iconf),ESPptval(1:nptthis,iconf),conffilepath(iconf))
		else !Don't show prompts
			call fitESP_calcESP(0,iESPtype,nptthis,ESPpt(:,1:nptthis,iconf),ESPptval(1:nptthis,iconf),conffilepath(iconf))
		end if
		
	end do
	
	if (ifloadconflist==1) then
		write(*,*) "Reloading the first file when Multiwfn boots up..."
		call dealloall(1)
		call readinfile(firstfilename,1)
	end if
	
!Reading ESP and coordinates of fitting points from Gaussian Iop(6/33=2) output. File containing geometry must be loaded to provide atom coordinates
else
	if (nconf==1) then
		write(*,"(a)") " Input path of the output file of Gaussian pop=MK/CHELPG task with IOp(6/33=2), e.g. C:\tsushima_yoshiko.out. &
		&Note that the geometry used in the Gaussian calculation must be exactly identical to current geometry"
		do while(.true.)
			read(*,"(a)") gauoutfilepath
			inquire(file=gauoutfilepath,exist=alive)
			if (alive) exit
			write(*,*) "Error: Cannot find the file, input again"
		end do
		write(*,"(a)") " Loading ESP data from "//trim(gauoutfilepath)
		iconf=1
		call loadgauESP_num(gauoutfilepath,nfitcen,nptthis) !Get number of points
		nESPpt(iconf)=nptthis
		maxESPpt=nptthis
		allocate(ESPptval(nptthis,iconf),ESPpt(3,nptthis,iconf))
		call loadgauESP(gauoutfilepath,nfitcen,fitcen(:,:,iconf),nptthis,ESPpt(:,:,iconf),ESPptval(:,iconf))
	else
		do iconf=1,nconf
			write(*,"(a)") " Loading ESP data from "//trim(conffilepath(iconf))
			call loadgauESP_num(conffilepath(iconf),nfitcen,nptthis) !Get number of points
			if (iconf==1) then !At the first conformer, allocate large enough size for ESPptval and ESPpt
				maxESPpt=nptthis*2
				allocate(ESPptval(maxESPpt,nconf),ESPpt(3,maxESPpt,nconf))
			end if
			nESPpt(iconf)=nptthis
			call loadgauESP(conffilepath(iconf),nfitcen,fitcen(:,:,iconf),nptthis,ESPpt(:,1:nptthis,iconf),ESPptval(1:nptthis,iconf))
		end do
	end if
end if
write(*,*)


!For standard two-stage RESP fitting, or one-stage fitting with eqv. cons. on H in -CH3, -CH2-, =CH2, &
!we need to identify bonding and generate relevant arrays
if (isel==1.or.ieqvcons==2) then
	!Generate bonding matrix, used to identify -CH3, -CH2-, =CH2 groups
20	bondedmat=0
	if (ideterbond==1) then !Guess bonding relationship
		do i=1,ncenter
			do j=i+1,ncenter
				if ( atomdist(i,j,1)<( covr(a(i)%index)+covr(a(j)%index) )*bondcrit ) bondedmat(i,j)=1
				bondedmat(j,i)=bondedmat(i,j)
			end do
		end do
	else if (ideterbond==2) then !Generate bonding matrix based on loaded connectivity matrix in specified .mol or .mol2 file
		call readmolconn(molfilepath)
		do i=1,ncenter
			do j=i+1,ncenter
				if (connmat(i,j)>0) bondedmat(i,j)=1
				bondedmat(j,i)=bondedmat(i,j)
			end do
		end do
	end if
	
	!Identify active C and H in stage 2 of standard RESP fitting, meantime generate equivalence &
    !list considering hydrogens in each -CH3, -CH2-, =CH2 is equivalent
	neqvlist_H=0
	eqvlist_H=0
	nCHlist=0
	do iatm=1,ncenter
		if (a(iatm)%index==6) then !Carbon
			nbond=sum(bondedmat(iatm,:))
			natt_H=0
			natt_C=0
			do jatm=1,ncenter
				if (bondedmat(iatm,jatm)==1) then
					if (a(jatm)%index==1) then
						natt_H=natt_H+1
						H_list(natt_H)=jatm
					else if (a(jatm)%index==6) then
						natt_C=natt_C+1
					end if
				end if
			end do
			!sp3 or methylene carbon, its charge should be fitted at stage 2 of standard RESP
			if (nbond==4.or.(natt_H==2.and.natt_C==1)) then
				nCHlist=nCHlist+1
				CHlist(nCHlist)=iatm
				neqvlist_H=neqvlist_H+1
				eqvlistlen_H(neqvlist_H)=1
				eqvlist_H(1,neqvlist_H)=iatm
				if (natt_H>0) then !Hydrogens attached to such carbon should have equivalent constraint
					neqvlist_H=neqvlist_H+1
					eqvlistlen_H(neqvlist_H)=natt_H
					do iatt_H=1,natt_H
						nCHlist=nCHlist+1
						CHlist(nCHlist)=H_list(iatt_H)
						eqvlist_H(iatt_H,neqvlist_H)=H_list(iatt_H)
					end do
				end if
			end if
		end if 
	end do
	!If there are atoms have not appeared in the equivalence constraint list, add them to a slot
	do icen=1,ncenter
		if (all(eqvlist_H/=icen)) then
			neqvlist_H=neqvlist_H+1
			eqvlistlen_H(neqvlist_H)=1
			eqvlist_H(1,neqvlist_H)=icen
		end if
	end do
    if (ieqvcons==3) then !Export the generated eqv. const.
        open(20,file="eqvcons_H.txt",status="replace")
        do itmp=1,neqvlist_H
            if (eqvlistlen_H(itmp)>1) then
                do jtmp=1,eqvlistlen_H(itmp)-1
                     write(20,"(i6,',')",advance='no') eqvlist_H(jtmp,itmp)
                end do
                write(20,"(i6)") eqvlist_H(jtmp,itmp)
            end if
        end do
        close(20)
        goto 30
    end if
end if


!Setting up charge constraint
if (ichgcons==1) then !Load setting from plain text file
	write(*,"(' Loading charge constraint setting from ',a)") trim(chgconsfilepath)
	open(10,file=chgconsfilepath,status="old")
	nchgcons=totlinenum(10,1) !The number of charge constraint terms
	allocate(chgconsnatm(nchgcons),chgconsval(nchgcons),chgconsatm(nfitcen,nchgcons))
	rewind(10)
	do icons=1,nchgcons
		read(10,"(a)") c2000tmp
		isep=index(trim(c2000tmp),' ',back=.true.)
		read(c2000tmp(isep+1:),*) chgconsval(icons)
		call str2arr(c2000tmp(:isep-1),ntmp)
		chgconsnatm(icons)=ntmp
		call str2arr(c2000tmp(:isep-1),ntmp,chgconsatm(1:ntmp,icons))
		write(*,"(' Charge constraint',i4,':',i4,' atoms, charge:',f12.6)") icons,ntmp,chgconsval(icons)
	end do
	close(10)
else !No charge constraint
	nchgcons=0
	allocate(chgconsnatm(nchgcons),chgconsval(nchgcons),chgconsatm(nfitcen,nchgcons))
	write(*,*) "No charge constraint is imposed in this stage"
end if


!Setting up equivalence constraint
!First assume that no equivalence constraint is imposed
neqvlist=nfitcen
eqvlistlen=1
forall(i=1:neqvlist) eqvlist(1,i)=i
if (ieqvcons==1) then !Employ equivalence setting from plain text file
	eqvlist=0
	write(*,"(' Loading equivalence constraint setting from ',a)") trim(eqvconsfilepath)
	open(10,file=eqvconsfilepath,status="old")
	neqvlist=totlinenum(10,1)
	rewind(10)
	do ieqv=1,neqvlist
		read(10,"(a)") c2000tmp
		call str2arr(c2000tmp,ntmp)
		eqvlistlen(ieqv)=ntmp
		call str2arr(c2000tmp,ntmp,eqvlist(1:ntmp,ieqv))
	end do
	close(10)
	!If there are atoms have not appeared in the equivalence constraint list, add them to a slot
	do icen=1,nfitcen
		if (all(eqvlist/=icen)) then
			neqvlist=neqvlist+1
			eqvlistlen(neqvlist)=1
			eqvlist(1,neqvlist)=icen
		end if
	end do
else if (ieqvcons==2) then !Impose equivalence constraint on hydrogens in each CH2 and CH3 group for one-stage fitting
    if (isel==1) then !Standard two-stage fitting
        continue !Such a constraint will be imposed later in the second fitting stage
    else if (isel==2) then !One-stage fitting
        neqvlist=neqvlist_H
	    eqvlistlen(1:neqvlist_H)=eqvlistlen_H(1:neqvlist_H)
	    eqvlist(1:size(eqvlist_H,1),1:neqvlist_H)=eqvlist_H(:,1:neqvlist_H)
        if (naddcen>0) then !No constraint for additional centers
            do ieqv=1,naddcen
                neqvlist=neqvlist+1
                eqvlistlen(neqvlist)=1
                eqvlist(1,neqvlist)=ncenter+ieqv
            end do
        end if
    end if
end if

!write(*,*) neqvlist,nfitcen
!do ieqv=1,neqvlist
!    write(*,*) eqvlist(1:eqvlistlen(ieqv),ieqv)
!end do
!pause
call showeqvcons(neqvlist,nfitcen,eqvlistlen,eqvlist)
write(*,*)


!------ Start charge fitting now!
!------ Start charge fitting now!
!Do first step of standard RESP fitting
if (isel==1) then
	write(*,*) "**** Stage 1: RESP fitting under weak hyperbolic penalty"
	call RESPiter(hyper_a_1,hyper_b,nconf,confweight,nfitcen,fitcen,maxESPpt,nESPpt,ESPpt,ESPptval,&
	atmchg_stage1,nchgcons,chgconsnatm,chgconsatm,chgconsval,neqvlist,eqvlistlen,eqvlist,iESPtype,maxRESPiter,RESPconv)
else if (isel==2) then !One-stage ESP fitting	
	write(*,*) "One-stage restrainted ESP fitting iteration has started"
	call RESPiter(hyper_a,hyper_b,nconf,confweight,nfitcen,fitcen,maxESPpt,nESPpt,ESPpt,ESPptval,&
	atmchg,nchgcons,chgconsnatm,chgconsatm,chgconsval,neqvlist,eqvlistlen,eqvlist,iESPtype,maxRESPiter,RESPconv)
end if

!Do second step of standard RESP fitting
if (isel==1) then
	write(*,*)
	write(*,*) "**** Stage 2: RESP fitting under strong hyperbolic penalty"
	if (nCHlist>0) then
		!Equivalence constraint: Hydrogens in -CH3, =CH2, -CH2-
        neqvlist=neqvlist_H
	    eqvlistlen(1:neqvlist_H)=eqvlistlen_H(1:neqvlist_H)
	    eqvlist(1:size(eqvlist_H,1),1:neqvlist_H)=eqvlist_H(:,1:neqvlist_H)
        if (naddcen>0) then !No constraint for additional centers
            do ieqv=1,naddcen
                neqvlist=neqvlist+1
                eqvlistlen(neqvlist)=1
                eqvlist(1,neqvlist)=ncenter+ieqv
            end do
        end if
		call showeqvcons(neqvlist,nfitcen,eqvlistlen,eqvlist)
		write(*,*) "Fitting objects: sp3 carbons, methyl carbons and hydrogens attached to them"
		write(*,*) "Indices of these atoms:"
		do i=1,nCHlist
			write(*,"(i5,a)",advance='no') CHlist(i),a(CHlist(i))%name
			if (mod(i,10)==0) write(*,*)
		end do
		
		deallocate(chgconsatm,chgconsnatm,chgconsval)
		nchgcons=nfitcen-nCHlist
		allocate(chgconsnatm(nchgcons),chgconsval(nchgcons),chgconsatm(nfitcen,nchgcons))
		chgconsnatm=1
		icons=0
		do icen=1,nfitcen
			if (all(CHlist(1:nCHlist)/=icen)) then !Charge of this atom should be fixed
				icons=icons+1
				chgconsval(icons)=atmchg_stage1(icen)
				chgconsatm(1,icons)=icen
			end if
		end do
		!If current folder has chgcons_stage2.txt, then export constrainted charges in the RESP stage 2 to this file
		inquire(file="chgcons_stage2.txt",exist=alive)
		if (alive) then
			open(10,file="chgcons_stage2.txt",status="replace")
			do idx=1,icons
				write(10,"(i6,f12.6)") chgconsatm(1,idx),chgconsval(idx)
			end do
			close(10)
			write(*,"(/,a)") " The charges that kept fixed at stage 2 of RESP fitting &
			&has been written to chgcons_stage2.txt in current folder"
		end if
		
		write(*,*)
		call RESPiter(hyper_a_2,hyper_b,nconf,confweight,nfitcen,fitcen,maxESPpt,nESPpt,ESPpt,ESPptval,&
		atmchg,nchgcons,chgconsnatm,chgconsatm,chgconsval,neqvlist,eqvlistlen,eqvlist,iESPtype,maxRESPiter,RESPconv)
	else
		write(*,*) "Stage 2 of standard RESP fitting is skipped since no atom needs to be fitted"
		atmchg=atmchg_stage1
	end if
	
end if


!Output summary
write(*,*)
write(*,*) "  Center       Charge"
do icen=1,nfitcen
    if (icen<=ncenter) then
    	write(*,"(i6,'(',a,')',f15.10)") icen,ind2name(a(icen)%index),atmchg(icen)
    else
        write(*,"(i6,'(',a,')',f15.10)") icen,"X ",atmchg(icen)
    end if
end do
write(*,"(' Sum of charges:',f15.10)") sum(atmchg)

!Calculate RMSE and RRMSE
weiRMSE=0;weiRRMSE=0
do iconf=1,nconf
	RMSE=0D0
	do ipt=1,nESPpt(iconf)
		atmchgesp=0D0
		do icen=1,nfitcen
			dis=dsqrt( sum((ESPpt(:,ipt,iconf)-fitcen(:,icen,iconf))**2) )
			atmchgesp=atmchgesp+atmchg(icen)/dis
		end do
		RMSE=RMSE+(ESPptval(ipt,iconf)-atmchgesp)**2
	end do
	RRMSE=dsqrt(RMSE/sum(ESPptval(1:nESPpt(iconf),iconf)**2))
	RMSE=dsqrt(RMSE/nESPpt(iconf))
	if (nconf>1) then
		write(*,"(' Conformer:',i5,'   RMSE:',f12.6,'   RRMSE:',f12.6)") iconf,RMSE,RRMSE
		weiRMSE=weiRMSE+RMSE*confweight(iconf)
		weiRRMSE=weiRRMSE+RRMSE*confweight(iconf)
	else
		write(*,"(' RMSE:',f12.6,'   RRMSE:',f12.6)") RMSE,RRMSE
	end if
end do
if (nconf>1) write(*,"(' Weighted RMSE:',f12.6,'   Weighted RRMSE',f12.6)") weiRMSE,weiRRMSE

!Show fragment charge
if (allocated(frag1)) then
    write(*,"(/,' Fragment charge:',f14.8)") sum(atmchg(frag1))
    write(*,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(atmchg(frag1))
end if

!Export .chg file
if (nconf==1) then
    if (nfitcen==ncenter) then
        call outatmchg(10,atmchg(:))
    else
        call outallchg(10,atmchg(:),fitcen,nfitcen)    
    end if
else
    write(*,"(/,a)") " Note: Because present calculation involves multiple conformers, the result cannot be exported to .chg file"
end if

deallocate(fitcen,fitcenvdwr,nESPpt,ESPptval,ESPpt,chgconsatm,chgconsnatm,chgconsval)
deallocate(atmchg,atmchg_stage1,eqvlistlen,eqvlist)

end do maincyc

end subroutine



!!------ Perform RESP iteration
subroutine RESPiter(hyper_a,hyper_b,nconf,confweight,nfitcen,fitcen,maxESPpt,nESPpt,ESPpt,ESPptval,&
atmchg,nchgcons,chgconsnatm,chgconsatm,chgconsval,neqvlist,eqvlistlen,eqvlist,iESPtype,maxiter,convcrit)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 hyper_a,hyper_b
integer nconf,nfitcen,maxESPpt,nchgcons,neqvlist,iESPtype
integer nESPpt(nconf)
real*8 confweight(nconf),fitcen(3,nfitcen,nconf),ESPpt(3,maxESPpt,nconf),ESPptval(maxESPpt,nconf),atmchg(nfitcen),chgconsval(nchgcons)
integer chgconsnatm(nchgcons),chgconsatm(nfitcen,nchgcons),eqvlistlen(nfitcen),eqvlist(500,nfitcen)
real*8,allocatable :: Bvec(:),Amat(:,:),Amat_bk(:,:),Amat_tmp(:,:),Amatinv(:,:),qvec(:),qvec_old(:)
real*8,allocatable :: Bveceqv(:),Amateqv(:,:),Amateqvinv(:,:),qveceqv(:) !The counterpart of Amat,Bvec,qvec when considering equivalence contraint
integer maxiter
real*8 convcrit

if (hyper_a==0) then
	write(*,*) "Since restraint strength was set to zero, no iteration will be carried out"
else
	write(*,"(' Convergence criterion:',f14.10)") convcrit
	write(*,"(' Hyperbolic restraint strength (a):',f9.6,'    Tightness (b):',f9.6)") hyper_a,hyper_b
end if

! write(*,*) nESPpt,nconf
! write(*,*) neqvlist,nchgcons,nfitcen
! do ieqv=1,neqvlist
! 	write(*,*) ieqv,eqvlist(1:eqvlistlen(ieqv),ieqv)
! end do
! do icons=1,nchgcons
! 	write(*,*) icons,chgconsatm(1:chgconsnatm(icons),icons),chgconsval(icons)
! end do
! do iconf=1,nconf
! 	write(*,*) iconf,confweight(iconf)
! end do

matdim=nfitcen+nchgcons+1
allocate(Bvec(matdim),Amat(matdim,matdim),Amat_bk(matdim,matdim),Amatinv(matdim,matdim),qvec(matdim),qvec_old(matdim))

!Forming Amat
Amat=0D0
do icen=1,nfitcen
	do jcen=icen,nfitcen
		do iconf=1,nconf
			tmp=0
			do ipt=1,nESPpt(iconf)
				dis1=dsqrt( sum((ESPpt(:,ipt,iconf)-fitcen(:,icen,iconf))**2) )
				dis2=dsqrt( sum((ESPpt(:,ipt,iconf)-fitcen(:,jcen,iconf))**2) )
				tmp=tmp+1D0/dis1/dis2
			end do
			Amat(icen,jcen)=Amat(icen,jcen)+tmp*confweight(iconf)
		end do
		Amat(jcen,icen)=Amat(icen,jcen)
	end do
end do
Amat(nfitcen+1,:nfitcen)=1D0
Amat(:nfitcen,nfitcen+1)=1D0
do icons=1,nchgcons
	do idx=1,chgconsnatm(icons)
		icen=chgconsatm(idx,icons)
		Amat(icen,nfitcen+1+icons)=1
		Amat(nfitcen+1+icons,icen)=1
	end do
end do
Amat_bk=Amat !Backup, its diagonal terms will be used during RESP calculation

!Forming Bvec
Bvec=0D0
do icen=1,nfitcen
	do iconf=1,nconf
		tmp=0
		do ipt=1,nESPpt(iconf)
			dis=dsqrt( sum((ESPpt(:,ipt,iconf)-fitcen(:,icen,iconf))**2) )
			tmp=tmp+ESPptval(ipt,iconf)/dis
		end do
		Bvec(icen)=Bvec(icen)+tmp*confweight(iconf)
	end do
end do
if (iESPtype==1) then !Take nuclei into account
	Bvec(nfitcen+1)=sum(a(:)%charge)-nelec
else if (iESPtype==2) then !Do not take nuclei into account
	Bvec(nfitcen+1)=-nelec
else if (iESPtype==3) then !Electronic transition density
	Bvec(nfitcen+1)=0
end if
do icons=1,nchgcons
	Bvec(nfitcen+1+icons)=chgconsval(icons)
end do

matdimeqv=neqvlist+nchgcons+1
allocate(Amateqv(matdimeqv,matdimeqv),Amateqvinv(matdimeqv,matdimeqv),Bveceqv(matdimeqv),qveceqv(matdimeqv))
allocate(Amat_tmp(matdimeqv,matdim)) !Used to temporarily store Amat with contracted rows due to eqv. cons.
qvec=0

do iter=1,maxiter
	qvec_old=qvec
	
	!Only diagonal terms have an additional term due to hyperbolic restraint
	do icen=1,nfitcen !Penalty function is only applied to non-hydrogen atoms
        if (icen<=ncenter) then
		    if (a(icen)%index==1) cycle
		end if
        Amat(icen,icen)=Amat_bk(icen,icen)+ hyper_a/dsqrt(qvec(icen)**2+hyper_b**2)
	end do
	
	!Construct matrices with consideration of equivalence constraint
	Amat_tmp=0
	Amateqv=0
	Bveceqv=0
	!Contract row of Amat to form Amat_tmp and meantime contract Bvec to form Bveceqv
	do ieqv=1,neqvlist
		do idx=1,eqvlistlen(ieqv)
			irow=eqvlist(idx,ieqv)
			Amat_tmp(ieqv,:)=Amat_tmp(ieqv,:)+Amat(irow,:)
			Bveceqv(ieqv)=Bveceqv(ieqv)+Bvec(irow)
		end do
	end do
	do idx=1,nchgcons+1
		Amat_tmp(neqvlist+idx,:)=Amat(nfitcen+idx,:)
		Bveceqv(neqvlist+idx)=Bvec(nfitcen+idx)
	end do
	!Contract column of Amat_tmp to form Amateqv
	do ieqv=1,neqvlist
		do idx=1,eqvlistlen(ieqv)
			icol=eqvlist(idx,ieqv)
			Amateqv(:,ieqv)=Amateqv(:,ieqv)+Amat_tmp(:,icol)
		end do
	end do
	do idx=1,nchgcons+1
		Amateqv(:,neqvlist+idx)=Amat_tmp(:,nfitcen+idx)
	end do
	
	!Evaluate and update charges
	Amateqvinv=invmat(Amateqv,matdimeqv)
	qveceqv=matmul(Amateqvinv,Bveceqv)
	do ieqv=1,neqvlist
		do idx=1,eqvlistlen(ieqv)
			qvec(eqvlist(idx,ieqv))=qveceqv(ieqv)
		end do
	end do
	
	if (hyper_a==0) exit
	deltamax=maxval(abs(qvec(1:nfitcen)-qvec_old(1:nfitcen)))
	write(*,"(' Iter:',i4,'   Maximum charge variation:',f16.10)") iter,deltamax
	if (deltamax<convcrit) exit
	
end do
if (iter==maxiter+1) then
	write(*,*) "Error: Convergence failed!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
else if (hyper_a/=0) then
	write(*,*) "Successfully converged!"
end if

atmchg=qvec(1:nfitcen)
end subroutine


!!----- Show atom equivalence constraint in ESP fitting
subroutine showeqvcons(neqvlist,nfitcen,eqvlistlen,eqvlist)
use defvar
implicit real*8 (a-h,o-z)
integer neqvlist,nfitcen
integer eqvlistlen(nfitcen),eqvlist(500,nfitcen)
if (any(eqvlistlen(1:neqvlist)>1)) then
	write(*,*) "Atom equivalence constraint imposed in this fitting stage:"
	icons=0
	do ieqv=1,neqvlist
		if (eqvlistlen(ieqv)>1) then
			icons=icons+1
			write(*,"(' Constraint',i4,':')",advance='no') icons
			do idx=1,eqvlistlen(ieqv)
				iatm=eqvlist(idx,ieqv)
                if (iatm<=ncenter) then
				    write(*,"(i5,'(',a,')')",advance='no') iatm,a(iatm)%name
                else
				    write(*,"(i5,'(',a,')')",advance='no') iatm,"X "
                end if
			end do
			write(*,*)
		end if
	end do
else
	write(*,*) "No atom equivalence constraint is imposed in this fitting stage"
end if
end subroutine



!!----------- Generate equivalence list according to point group of provided fragment geometry and write to eqvcons_PG.txt
subroutine genPGeqvcons
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000,pglabel*3,c80tmp*80,selectyn
real*8 coord(3,ncenter)
integer atmlist(ncenter),atmindex(ncenter),classnatm(ncenter),classidx(ncenter,ncenter)
toler=0.03D0 !loose
open(10,file="eqvcons_PG.txt",status="replace")
write(*,"(a)") " Note: You can change the tolerance for detecting point group to e.g. 0.02 by inputting ""t 0.02"". The default tolerance is 0.03. &
&If Multiwfn shows ""ERROR: Too many symmetry operations"", you should set the tolerance to a smaller value and retry"
do while(.true.)
    write(*,*)
    write(*,*) "Input indices of the atoms in the fragment, e.g. 3,6-10,12,14"
    write(*,*) "To add all atoms in the system, input ""a"""
    write(*,*) "To exit, input ""q"""
    read(*,"(a)") c2000tmp
    if (index(c2000tmp,'q')/=0) then
        exit
    else if (index(c2000tmp,'t')/=0) then
        read(c2000tmp,*) c80tmp,toler
        write(*,"(' The tolerance has been set to',f8.4)") toler
    else
        if (index(c2000tmp,'a')/=0) then
            nselatm=ncenter
            forall(iatm=1:ncenter) atmlist(iatm)=iatm
        else
            call str2arr(c2000tmp,nselatm,atmlist)
        end if
        atmindex(1:nselatm)=a(atmlist(1:nselatm))%index
        coord(1,1:nselatm)=a(atmlist(1:nselatm))%x*b2a
        coord(2,1:nselatm)=a(atmlist(1:nselatm))%y*b2a
        coord(3,1:nselatm)=a(atmlist(1:nselatm))%z*b2a
        call PG_eqvatm(nselatm,atmindex(1:nselatm),coord(:,1:nselatm),toler,pglabel,nclass,classnatm(1:nselatm),classidx(1:nselatm,1:nselatm))
        nwrite=count(classnatm(1:nclass)>1)
        !Using loose tolerance, the equivalence atoms can also be detected, but sometimes point group is failed to recognized
        if (pglabel/=" ") then
            write(*,"(' Detected point group: ',a)") pglabel
        else
            !Even in many cases point group cannot be recognized, the symmetry-equivalence atoms are still correctly recognized
            !write(*,*) "Failed to detect point group"
        end if
        if (nwrite==0) then
            write(*,*) "No symmetry-equivalence atoms were found"
            cycle
        end if
        write(*,"(' Number of symmetry-equivalence classes:',i5)") nwrite
        ic=0
        do iclass=1,nclass
            ntmp=classnatm(iclass)
            if (ntmp==1) cycle
            ic=ic+1
            write(*,"(' Class',i5,a,':',i5,' atoms')") ic,' ('//a(classidx(1,iclass))%name//')',ntmp
            do iatm=1,ntmp
                write(*,"(i5)",advance="no") atmlist(classidx(iatm,iclass))
                if (iatm/=ntmp) write(*,"(a)",advance="no") ','
            end do
            write(*,*)
        end do
        write(*,*) "Accept and append to eqvcons_PG.txt in current folder? (y/n)"
        read(*,*) selectyn
        if (selectyn=='y') then
            do iclass=1,nclass
                ntmp=classnatm(iclass)
                if (ntmp==1) cycle
                do iatm=1,ntmp
                    write(10,"(i5)",advance="no") atmlist(classidx(iatm,iclass))
                    if (iatm/=ntmp) write(10,"(a)",advance="no") ','
                end do
                write(10,*)
            end do
            write(*,*) "The constraints have been appended to eqvcons_PG.txt in current folder"
        end if
    end if
end do
close(10)
write(*,*) "Generating eqvcons_PG.txt is finalized!"
end subroutine




!!------------ Calculate MK and CHELPG charges
! igridtype=1:MK   igridtype=2:CHELPG
! This module does not have wealth of features as RESP routine, but it is very clear and have these special points:
! (1) Support additional fitting center (2) Support external fitting points (3) Fitting points could be exported
! (4) Can calculate TrEsp (transition charge from electrostatic potential), see manual
subroutine fitESP(igridtype)
use util
use defvar
use functions
implicit real*8 (a-h,o-z)
integer igridtype
character(len=200) addcenfile,extptfile,chgfile,gauoutfilepath,outchgfilepath
character selectyn,c80tmp*80,c2000tmp*2000
integer :: iESPtype=1,ioutfitptval=0,iradiisel=1
!About distribution of fitting points
integer :: nMKlayer=4
integer :: MKatmlist(ncenter) !Record index of atoms used to construct MK fitting points
integer :: CHELPGatmlist(ncenter) !Record index of atoms used to construct CHELPG fitting points
real*8 :: espfitvdwr(nelesupp),sclvdwlayer(100)=(/1.4D0,1.6D0,1.8D0,2D0,(0D0,i=5,100)/)
real*8 :: fitspc=0.566917796573677D0 !0.3/b2a, spacing between grid for CHELPG
real*8 :: extdis=5.29123276802099D0 !2.8/b2a, extend 2.8 Angstrom to each side for CHELPG
real*8 :: MKptdens=1.68017136515525D0 !6D0*b2a**2, 6.0 Angstrom**2 point density per for MK. Multiply by b2a**2 to convert to Bohr**2
!Arrays used in ESP fitting
real*8,allocatable :: ESPptval(:),ESPpt(:,:),fitcen(:,:) !x/y/z,index
real*8,allocatable :: Bvec(:),Amat(:,:),Amatinv(:,:),qvec(:)
real*8,allocatable :: fitcenvdwr(:)
real*8,allocatable :: cenchg(:) !Final result
real*8,allocatable :: ESPerr(:)

if (.not.allocated(b).and.ifiletype/=4) then
    write(*,"(a)") " Error: Your input file does not contain wavefunction information at all, &
    &evidently the electrostatic potential fitted charges cannot be calculated! You should use &
    &e.g. .wfn/wfx/mwfn/fch/molden... as input file, see Section 2.5 of Multiwfn manual"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

iaddcen=0 !If give additional centers
naddcen=0
iuseextpt=0 !If use external points
iskipespcalc=0 !If read ESP from external file directly rather than calculate here
iloadgau=0
iloadchg=0
nMKatoms=ncenter !The number of atoms used to construct MK fitting points
nCHELPGatoms=ncenter !The number of atoms used to construct CHELPG fitting points
forall(i=1:ncenter) MKatmlist(i)=i
forall(i=1:ncenter) CHELPGatmlist(i)=i

10 do while(.true.)
	write(*,*)
	if (igridtype==1) write(*,*) "            -------------- Calculation of MK charges --------------"
	if (igridtype==2) write(*,*) "            ------------ Calculation of CHELPG charges ------------"
	if (iloadchg==0) write(*,*) "-3 Toggle using atomic charges in external file, current: No"
	if (iloadchg==1) write(*,*) "-3 Toggle using atomic charges in external file, current: Yes"
	if (iaddcen==0) write(*,*) "-2 Toggle loading additional fitting centers from external file, current: No"
	if (iaddcen==1) write(*,*) "-2 Toggle loading additional fitting centers from external file, current: Yes"
	if (iuseextpt==0) write(*,*) "-1 Toggle using fitting points in external file, current: No"
	if (iuseextpt==1) write(*,*) "-1 Toggle using fitting points in external file, current: Yes"
	write(*,*) "0 Return"
	write(*,*) "1 Start calculation!"
	if (iuseextpt==0) then
		if (igridtype==1) then
			write(*,"(' 2 Set number of points per Angstrom^2, current:',f10.3)") MKptdens/b2a**2 !Temporary convert to Angstrom**2 for convention
			write(*,"(' 3 Set number and scale factors of layers of MK fitting points')")
			if (nMKatoms==ncenter) then
				write(*,*) "4 Choose the atoms used for generating fitting points, current: All atoms"
			else
				write(*,"(a,i5,a)") " 4 Choose the atoms used for generating fitting points, current:",nMKatoms," atoms"
			end if
		else if (igridtype==2) then
			write(*,"(' 2 Set grid spacing, current:',f7.3,' Bohr (',f7.3,' Angstrom)')") fitspc,fitspc*b2a
			write(*,"(' 3 Set box extension, current:',f7.3,' Bohr (',f7.3,' Angstrom)')") extdis,extdis*b2a
			if (nCHELPGatoms==ncenter) then
				write(*,*) "4 Choose the atoms used for generating fitting points, current: All atoms"
			else
				write(*,"(a,i5,a)") " 4 Choose the atoms used for generating fitting points, current:",nCHELPGatoms," atoms"
			end if
		end if
	end if
	if (iESPtype==1) write(*,*) "5 Choose ESP type, current: Nuclear + Electronic"
	if (iESPtype==2) write(*,*) "5 Choose ESP type, current: Electronic"
	if (iESPtype==3) write(*,*) "5 Choose ESP type, current: Transition electronic"
	if (ioutfitptval==0) write(*,*) "6 Toggle if exporting fitting points with ESP after the task, current: No"
	if (ioutfitptval==1) write(*,*) "6 Toggle if exporting fitting points with ESP after the task, current: Yes"
 	if (iloadgau==0) write(*,"(a)") " 7 Toggle if reading fitting points and ESP values from Gaussian output file of pop=MK/CHELPG task with IOp(6/33=2), current: No"
 	if (iloadgau==1) write(*,"(a)") " 7 Toggle if reading fitting points and ESP values from Gaussian output file of pop=MK/CHELPG task with IOp(6/33=2), current: Yes"
	if (iradiisel==1) write(*,*) "10 Choose the atomic radii used in fitting, current: Automatic"
	if (iradiisel==2) write(*,*) "10 Choose the atomic radii used in fitting, current: Scaled UFF"
	if (iradiisel==3) write(*,*) "10 Choose the atomic radii used in fitting, current: Customized"
	read(*,*) isel
	
	if (isel==-3) then
		if (iloadchg==1) then
			iloadchg=0
		else
			write(*,*) "Input the path of the .chg file containing atomic charges, e.g. C:\H2O.chg"
			do while(.true.)
				read(*,"(a)") chgfile
				inquire(file=chgfile,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			iloadchg=1
			write(*,"(a)") " OK, the atomic charges contained in this file will be directly used and no ESP fitting charges will be yielded"
		end if
	else if (isel==-2) then
		if (iaddcen==1) then
			iaddcen=0
		else
			write(*,"(a)") " Input the path of the file recording coordinates of additional fitting centers, e.g. C:\ll_sunshine\Riko.txt"
			do while(.true.)
				read(*,"(a)") addcenfile
				inquire(file=addcenfile,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			iaddcen=1
			write(*,*) "Done!"
		end if
	else if (isel==-1) then
		if (iuseextpt==1) then
			iuseextpt=0
		else
			write(*,"(a)") " Input the path of the file recording coordinates of ESP fitting points, e.g. C:\ll_sunshine\You.txt"
			do while(.true.)
				read(*,"(a)") extptfile
				inquire(file=extptfile,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			iuseextpt=1
			write(*,*) "OK, the points recorded in this file will be used as fitting points"
		end if
	else if (isel==0) then
		Return
	else if (isel==1) then
		exit
	else if (isel==2) then
		if (igridtype==1) then
			write(*,*) "Input a number (in per Angstrom^2), e.g. 6.0"
			read(*,*) MKptdens
			MKptdens=MKptdens*b2a**2 !Convert to per Bohr^2
		else if (igridtype==2) then
			write(*,*) "Input a value in Bohr, e.g. 0.5"
			read(*,*) fitspc
		end if
	else if (isel==3) then
		if (igridtype==1) then
			write(*,*) "Current layers of MK fitting points:"
			do ilayer=1,nMKlayer
				write(*,"(' Layer',i3,':',f8.4)") ilayer,sclvdwlayer(ilayer)
			end do
			nMKlayer=0
			write(*,*)
			do while(.true.)
				nMKlayer=nMKlayer+1
				write(*,"(' Input scale factor (w.r.t atomic vdW radius) for layer',i4)") nMKlayer
				write(*,*) "If input ""q"", the setting will be finished"
				if (sclvdwlayer(nMKlayer)/=0) write(*,"(' If pressing ENTER button directly, current value',f8.4,' will be retained')") sclvdwlayer(nMKlayer)
				read(*,"(a)") c80tmp
				if (index(c80tmp,'q')/=0) then
					sclvdwlayer(nMKlayer:)=0
					nMKlayer=nMKlayer-1
					exit
				end if
				if (c80tmp/=" ") read(c80tmp,*) sclvdwlayer(nMKlayer)
			end do
			write(*,*) "Current layer of MK fitting points:"
			do ilayer=1,nMKlayer
				write(*,"(' Layer',i3,'  Scale factor:',f8.4)") ilayer,sclvdwlayer(ilayer)
			end do
		else if (igridtype==2) then
			write(*,*) "Input a value in Bohr, e.g. 6.5"
			read(*,*) extdis
		end if
	else if (isel==4) then
        if (igridtype==1) then
		    write(*,*) "Input the indices of the atoms used to construct MK fitting points"
		    write(*,*) "e.g. 1-5,8,10-12"
		    read(*,"(a)") c2000tmp
		    call str2arr(c2000tmp,nMKatoms,MKatmlist)
        else if (igridtype==2) then
		    write(*,*) "Input the indices of the atoms used to construct CHELPG fitting points"
		    write(*,*) "e.g. 1-5,8,10-12"
		    read(*,"(a)") c2000tmp
		    call str2arr(c2000tmp,nCHELPGatoms,CHELPGatmlist)
        end if
		write(*,*) "Done!"
	else if (isel==5) then
		write(*,*) "Choose ESP type:"
		write(*,*) "1 Nuclear + Electronic"
		write(*,*) "2 Electronic only"
		write(*,*) "3 Transition electronic (Mainly used for deriving TrEsp charges)"
		read(*,*) iESPtype
	else if (isel==6) then
		if (ioutfitptval==1) then
			ioutfitptval=0
		else if (ioutfitptval==0) then
			ioutfitptval=1
		end if
	else if (isel==7) then
		if (iloadgau==1) then
			iloadgau=0
		else if (iloadgau==0) then
			iloadgau=1
			write(*,"(a)") " Input file path of the Gaussian pop=MK/CHELPG task with IOp(6/33=2) keyword, e.g. C:\tsushima_yoshiko.out. &
			&Note that the geometry used in the Gaussian calculation must be exactly identical to current geometry"
			do while(.true.)
				read(*,"(a)") gauoutfilepath
				inquire(file=gauoutfilepath,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			write(*,"(a)") " OK, ESP fitting points with ESP values will be directly loaded from this file during calculation"
		end if
	else if (isel==10) then
		write(*,*) "Please choose the way of determining atomic radii:"
		write(*,"(a)") " 1 Automatic (The radii actually used will depend on the way of distributing fitting points)"
		write(*,"(a)") " 2 UFF radii scaled by 1/1.2, defined for entire periodic table"
		write(*,"(a)") " 3 Customized (The radii will be loaded from external file during fitting)"
		read(*,*) iradiisel
	end if
end do


!Generate fitting centers, fitting points and calculate ESP value using internal code
if (iloadgau==0) then
	!Set and check sanity of vdW radii used in fitting
	espfitvdwr=-1
	if (iradiisel==1) then
		call setESPfitvdwr(igridtype,espfitvdwr)
	else if (iradiisel==2) then !Scaled UFF
		call setESPfitvdwr(0,espfitvdwr)
	else if (iradiisel==3) then !Customize radii
		call setESPfitvdwr(-1,espfitvdwr)
	end if

	!Set total number of fitting centers
	naddcen=0
	if (iaddcen==1) then !Load additional fitting centers
		open(10,file=addcenfile,status="old")
		read(10,*) naddcen
	end if
	nfitcen=ncenter+naddcen
	allocate(fitcen(3,nfitcen),fitcenvdwr(nfitcen))
	!Generate positions and vdW radii of fitting centers
	do iatm=1,ncenter
		fitcen(1,iatm)=a(iatm)%x
		fitcen(2,iatm)=a(iatm)%y
		fitcen(3,iatm)=a(iatm)%z
		fitcenvdwr(iatm)=espfitvdwr(a(iatm)%index) !vdW radius for each fitting center
	end do
	if (iaddcen==1) then
		do icen=ncenter+1,ncenter+naddcen
			read(10,*) fitcen(:,icen)
			fitcen(:,icen)=fitcen(:,icen)/b2a
			fitcenvdwr(icen)=0 !The important thing is reproducing ESP on molecular surface, therefore extra point should not affect construction of fitting points
! 			write(*,"(' Please input radius for extra fitting center',i4,' (Bohr), e.g. 0.4')") icen
! 			read(*,*) fitcenvdwr(icen)
		end do
		close(10)
        write(*,*) "Coordinate of loaded additional fitting centers (Angstrom):"
        do icen=1,naddcen
            write(*,"(i6,3f12.6)") icen,fitcen(:,ncenter+icen)*b2a
        end do
	end if

	!Generate position of fitting points
	write(*,*)
	if (iuseextpt==0) then !Count number and generate coordinates of fitting points
		if (igridtype==1) then !MK
			allocate(ESPpt(3,0)) !Temporarily assign a minimal length
			call setMKpt(1,nfitcen,fitcen,fitcenvdwr,nESPpt,ESPpt,sclvdwlayer,MKptdens,nMKlayer,nMKatoms,MKatmlist) !The returned nESPpt is upper limit
			deallocate(ESPpt);allocate(ESPptval(nESPpt),ESPpt(3,nESPpt))
			call setMKpt(2,nfitcen,fitcen,fitcenvdwr,nESPpt,ESPpt,sclvdwlayer,MKptdens,nMKlayer,nMKatoms,MKatmlist)
		else if (igridtype==2) then !CHELPG
			allocate(ESPpt(3,0)) !Temporarily assign a minimal length
			call setCHELPGpt(1,nfitcen,fitcen,fitcenvdwr,nESPpt,ESPpt,extdis,fitspc,nCHELPGatoms,CHELPGatmlist) !Return actual nESPpt
			deallocate(ESPpt);allocate(ESPptval(nESPpt),ESPpt(3,nESPpt))
			call setCHELPGpt(2,nfitcen,fitcen,fitcenvdwr,nESPpt,ESPpt,extdis,fitspc,nCHELPGatoms,CHELPGatmlist)
		end if
	else if (iuseextpt==1) then !Directly use fitting points in external file. The ESP values may or may not be provided in the file
		open(10,file=extptfile,status="old")
		read(10,*) nESPpt
		if (nESPpt<0) then
			iskipespcalc=1 !If the number of fitting points is negative, that means the fourth column records ESP value and needn't to be recalculated
			write(*,*) "ESP value of all fitting points are read from external file directly"
		end if
		nESPpt=abs(nESPpt)
		write(*,"(' Total number of fitting points used:',i10)") nESPpt
		allocate(ESPptval(nESPpt),ESPpt(3,nESPpt))
		do i=1,nESPpt
			if (iskipespcalc==0) read(10,*) ESPpt(:,i)
			if (iskipespcalc==1) read(10,*) ESPpt(:,i),ESPptval(i)
		end do
		close(10)
	end if

	!Generate ESP value of fitting points
	if (iskipespcalc==0) call fitESP_calcESP(1,iESPtype,nESPpt,ESPpt,ESPptval,filename)

else !Reading ESP and coordinates of fitting points from Gaussian Iop(6/33=2) output. File containing geometry must be loaded to provide atom coordinates
	write(*,"(a)") " Loading ESP data from "//trim(gauoutfilepath)
	call loadgauESP_num(gauoutfilepath,nfitcen,nESPpt)
	allocate(fitcen(3,nfitcen),ESPptval(nESPpt),ESPpt(3,nESPpt))
	call loadgauESP(gauoutfilepath,nfitcen,fitcen,nESPpt,ESPpt,ESPptval)
end if

allocate(cenchg(nfitcen))
if (iloadchg==0) then
	!Calculate ESP fitting charges. See original paper of MK for detail of algorithem
	matdim=nfitcen+1
	allocate(Bvec(matdim),Amat(matdim,matdim),Amatinv(matdim,matdim),qvec(matdim))
	!Forming Amat
	Amat=0D0
	do icen=1,nfitcen
		do jcen=icen,nfitcen
			do ipt=1,nESPpt
				dis1=dsqrt( sum((ESPpt(:,ipt)-fitcen(:,icen))**2) )
				dis2=dsqrt( sum((ESPpt(:,ipt)-fitcen(:,jcen))**2) )
				Amat(icen,jcen)=Amat(icen,jcen)+1D0/dis1/dis2
			end do
			Amat(jcen,icen)=Amat(icen,jcen)
		end do
	end do
	Amat(matdim,:nfitcen)=1D0
	Amat(:nfitcen,matdim)=1D0
	!Forming Bvec
	Bvec=0D0
	do icen=1,nfitcen
		do ipt=1,nESPpt
			dis=dsqrt( sum((ESPpt(:,ipt)-fitcen(:,icen))**2) )
			Bvec(icen)=Bvec(icen)+ESPptval(ipt)/dis
		end do
	end do
	if (iESPtype==1) then !Take nuclei into account
		Bvec(matdim)=sum(a(:)%charge)-nelec !Net charge of the system
	else if (iESPtype==2) then !Do not take nuclei into account
		Bvec(matdim)=-nelec
	else if (iESPtype==3) then !Electronic transition density
		Bvec(matdim)=0
	end if
	Amatinv=invmat(Amat,matdim)
	qvec=matmul(Amatinv,Bvec)
	cenchg=qvec(1:nfitcen)
	deallocate(Bvec,Amat,Amatinv,qvec)
else if (iloadchg==1) then
	!Directly load charges of fitting centers from external files
	write(*,"(a)") " Loading charges of fitting centers from "//trim(chgfile)
	open(10,file=chgfile,status="old")
	do icen=1,nfitcen
		read(10,*) c80tmp,tmp1,tmp2,tmp3,cenchg(icen)
	end do
	close(10)
end if

!Output summary
write(*,*)
write(*,*) "  Center       Charge"
do i=1,ncenter+naddcen
	if (i<=ncenter) then
        write(*,"(i6,'(',a,')',f15.10)") i,ind2name(a(i)%index),cenchg(i)
    else
        write(*,"(i6,'(',a,')',f15.10)") i,"X ",cenchg(i)
    end if
end do
write(*,"(' Sum of charges:',f15.10)") sum(cenchg(1:nfitcen))

!Calculate RMSE and RRMSE
if (allocated(ESPerr)) deallocate(ESPerr)
allocate(ESPerr(nESPpt))
RMSE=0D0
do ipt=1,nESPpt
	atmchgesp=0D0
	do icen=1,nfitcen
		dis=dsqrt( sum((ESPpt(:,ipt)-fitcen(:,icen))**2) )
		atmchgesp=atmchgesp+cenchg(icen)/dis
	end do
	ESPerr(ipt)=abs(ESPptval(ipt)-atmchgesp)
	RMSE=RMSE+(ESPptval(ipt)-atmchgesp)**2
end do
RRMSE=dsqrt(RMSE/sum(ESPptval(1:nESPpt)**2))
RMSE=dsqrt(RMSE/nESPpt)
write(*,"(' RMSE:',f12.6,'   RRMSE:',f12.6)") RMSE,RRMSE

!Show fragment charge
if (allocated(frag1)) then
    write(*,"(/,' Fragment charge:',f14.8)") sum(cenchg(frag1))
    write(*,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(cenchg(frag1))
end if

!Exporting fitting points with ESP values
if (ioutfitptval==1) then
	do while(.true.)
		write(*,*)
		write(*,*) "0 Continue"
		write(*,*) "1 Export fitting points with ESP value to ESPfitpt.txt in current folder"
		write(*,*) "2 Export fitting points with ESP value to ESPfitpt.pqr in current folder"
		write(*,"(a)") " 3 Export fitting points with ESP reproduction error to ESPerr.pqr in current folder"
		read(*,*) ides
		if (ides==0) then
			exit
		else if (ides==1) then
			open(10,file="ESPfitpt.txt",status="replace")
			write(10,*) nESPpt
			do ipt=1,nESPpt
				write(10,"(3f12.6,f14.8)") ESPpt(:,ipt),ESPptval(ipt)
			end do
			write(*,*) "Done! Fitting points have been exported to ESPfitpt.txt in current folder"
			write(*,"(a)") " All units are in a.u. The first line shows the number of fitting points, &
			&the first three columns are X,Y,Z coordinates, the last column corresponds to ESP value"
			close(10)
		else if (ides==2) then
			open(10,file="ESPfitpt.pqr",status="replace")
			do ipt=1,nESPpt
				write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.3,a2)") &
				"HETATM",ipt,' '//"O "//' ',"MOL",'A',1,ESPpt(:,ipt)*b2a,ESPptval(ipt)*au2kcal,0.1D0,"O "
			end do
			write(*,"(a)") " Done! Fitting points have been exported to ESPfitpt.pqr in current folder. &
			&The ""charge"" column in this file corresponds to ESP value in kcal/mol. The radius column is meaningless"
			close(10)
		else if (ides==3) then
			open(10,file="ESPerr.pqr",status="replace")
			do ipt=1,nESPpt
				write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.3,a2)") &
                "HETATM",ipt,' '//"O "//' ',"MOL",'A',1,ESPpt(:,ipt)*b2a,ESPerr(ipt)*au2kcal,0.1D0,"O "
			end do
			write(*,"(a)") " Done! Fitting points have been exported to ESPerr.pqr in current folder. &
			&The ""charge"" column in this file corresponds to absolute different (in kcal/mol) between the exactly evaluated ESP &
			&and that evaluated based on atomic charges. The radius column is meaningless"
			close(10)
		end if
	end do
end if


!Export .chg file
if (nfitcen==ncenter) then
    call outatmchg(10,cenchg(:))
else
    call outallchg(10,cenchg(:),fitcen,nfitcen)    
end if

deallocate(ESPpt,ESPptval,fitcen)
if (allocated(fitcenvdwr)) deallocate(fitcenvdwr)
deallocate(cenchg)
goto 10 !Return to ESP fitting interface
end subroutine


!!------- Set MK fitting points for present system
!imode=1: Used to get upper limit of number of fitting points (nESPpt), the nESPpt and ESPpt that passed in could be any length
!imode=2: Fill coordinate of ESP fitting points to ESPpt array, in this case the passed-in ESPpt must have enough size
!MKatmlist records the real atoms used to constructing fitting points, nMKatoms is its length
subroutine setMKpt(imode,nfitcen,fitcen,fitcenvdwr,nESPpt,ESPpt,sclvdwlayer,MKptdens,nMKlayer,nMKatoms,MKatmlist)
use defvar
implicit real*8 (a-h,o-z)
integer nfitcen,nESPpt,nMKlayer,MKatmlist(nMKatoms)
real*8 fitcen(3,nfitcen),fitcenvdwr(nfitcen),ESPpt(3,nESPpt),sclvdwlayer(100),MKptdens
real*8,allocatable :: origsphpt(:,:)
if (imode==1) then !Count how many possible ESP points in total, the number is upper limit because some points will be pruned
	nESPpt=0
	do idx=1,nMKatoms !Note that only real atoms will be used to construct fitting points
		icen=MKatmlist(idx)
		do ilayer=1,nMKlayer
			numsphpt=nint(4D0*pi*(fitcenvdwr(icen)*sclvdwlayer(ilayer))**2 *MKptdens)
			nESPpt=nESPpt+numsphpt
		end do
	end do
else
	cutinnerscl=minval(sclvdwlayer(1:nMKlayer)) !If distance between a ESP point and any atom is smaller than this, the point will be discarded
	maxsphpt=nint(4D0*pi*(maxval(fitcenvdwr)*maxval(sclvdwlayer))**2 *MKptdens) !Find maximal possible number of points in unit sphere to allocate temporary origsphpt
	allocate(origsphpt(3,maxsphpt))
	iESPpt=0
	do idx=1,nMKatoms
		icen=MKatmlist(idx)
		do ilayer=1,nMKlayer
			radius=fitcenvdwr(icen)*sclvdwlayer(ilayer)
			numsphpt=nint(4D0*pi*radius**2 *MKptdens)
			call unitspherept(origsphpt,numsphpt) !Input expected number of point in unit sphere, return actual number of points
			origsphpt(:,1:numsphpt)=origsphpt(:,1:numsphpt)*radius
			do idir=1,3
				origsphpt(idir,1:numsphpt)=origsphpt(idir,1:numsphpt)+fitcen(idir,icen) !Move unit sphere to atomic center
			end do
			!Prune out the fitting points lying inside the intermost shell
			do ipt=1,numsphpt
				tmpx=origsphpt(1,ipt)
				tmpy=origsphpt(2,ipt)
				tmpz=origsphpt(3,ipt)
				iok=1
				do icen2=1,ncenter
					if (icen2==icen) cycle
					disptcensq=(fitcen(1,icen2)-tmpx)**2+(fitcen(2,icen2)-tmpy)**2+(fitcen(3,icen2)-tmpz)**2 !Distance between point and center
					if (disptcensq<(fitcenvdwr(icen2)*cutinnerscl)**2) then !Less than vdW RADIUS*cutinner of atom icen2, it should be ommitted
						iok=0
						exit
					end if
				end do
				if (iok==1) then
					iESPpt=iESPpt+1
					ESPpt(1,iESPpt)=tmpx
					ESPpt(2,iESPpt)=tmpy
					ESPpt(3,iESPpt)=tmpz
				end if
			end do
		end do
	end do
	nESPpt=iESPpt
	write(*,"(' Number of MK fitting points used:',i10)") nESPpt
end if
end subroutine

!!------- Set CHELPG fitting points for present system
!imode=1: Used to get number of fitting piints (nESPpt), the nESPpt and ESPpt that passed in could be any length
!imode=2: Fill coordinate of ESP fitting points to ESPpt array, in this case the passed-in ESPpt must have enough size
subroutine setCHELPGpt(imode,nfitcen,fitcen,fitcenvdwr,nESPpt,ESPpt,extdis,fitspc,nCHELPGatoms,CHELPGatmlist)
use defvar
implicit real*8 (a-h,o-z)
integer nfitcen,nESPpt,nCHELPGatoms,CHELPGatmlist(nCHELPGatoms)
real*8 fitcen(3,nfitcen),fitcenvdwr(nfitcen),disptcen(nfitcen)
real*8 ESPpt(3,nESPpt)
xlow=minval(fitcen(1,:))-extdis
xhigh=maxval(fitcen(1,:))+extdis
ylow=minval(fitcen(2,:))-extdis
yhigh=maxval(fitcen(2,:))+extdis
zlow=minval(fitcen(3,:))-extdis
zhigh=maxval(fitcen(3,:))+extdis
xlen=xhigh-xlow
ylen=yhigh-ylow
zlen=zhigh-zlow
nxfit=int(xlen/fitspc)+1
nyfit=int(ylen/fitspc)+1
nzfit=int(zlen/fitspc)+1
if (imode==1) then
	nESPpt=0
	do ix=0,nxfit
		do iy=0,nyfit
            do iz=0,nzfit
				tmpx=xlow+ix*fitspc
				tmpy=ylow+iy*fitspc
				tmpz=zlow+iz*fitspc
                iadd=0
                !If the point is not in the layer, ignore it
				do icen=1,nfitcen
					disptcen(icen)=dsqrt( (fitcen(1,icen)-tmpx)**2+(fitcen(2,icen)-tmpy)**2+(fitcen(3,icen)-tmpz)**2 )
					if (disptcen(icen)<=fitcenvdwr(icen)) exit
				end do
                if (icen==nfitcen+1) then
                    if (any(disptcen<=extdis)) iadd=1
                end if
                !If the point is belonging to an atom that in the given atom list, finally take it into account
                if (iadd==1) then
                    call pointcloseatom(tmpx,tmpy,tmpz,icloseatm)
                    if (any(CHELPGatmlist==icloseatm)) nESPpt=nESPpt+1
                end if
			end do
		end do
	end do
	write(*,"(' Number of CHELPG fitting points used:',i10)") nESPpt
else
	iESPpt=0
	do ix=0,nxfit
		do iy=0,nyfit
			do iz=0,nzfit
				tmpx=xlow+ix*fitspc
				tmpy=ylow+iy*fitspc
				tmpz=zlow+iz*fitspc
                iadd=0
                !If the point is not in the layer, ignore it
				do icen=1,nfitcen
					disptcen(icen)=dsqrt( (fitcen(1,icen)-tmpx)**2+(fitcen(2,icen)-tmpy)**2+(fitcen(3,icen)-tmpz)**2 )
					if (disptcen(icen)<=fitcenvdwr(icen)) exit
				end do
                if (icen==nfitcen+1) then
                    if (any(disptcen<=extdis)) iadd=1
                end if
                !If the point is belonging to an atom that in the given atom list, finally take it into account
                if (iadd==1) then
                    call pointcloseatom(tmpx,tmpy,tmpz,icloseatm)
                    if (any(CHELPGatmlist==icloseatm)) then
					    iESPpt=iESPpt+1
					    ESPpt(1,iESPpt)=tmpx
					    ESPpt(2,iESPpt)=tmpy
					    ESPpt(3,iESPpt)=tmpz
                    end if
                end if
			end do
		end do
	end do
end if
end subroutine


!!---- Return index of (real) atom that closest to a point
subroutine pointcloseatom(x,y,z,iclose)
use defvar
implicit real*8 (a-h,o-z)
real*8 x,y,z
integer iclose
dist2min=1E10
do iatm=1,ncenter
    if (a(iatm)%index<1) cycle
    dist2=(x-a(iatm)%x)**2+(y-a(iatm)%y)**2+(z-a(iatm)%z)**2
    if (dist2<dist2min) then
        iclose=iatm
        dist2min=dist2
    end if
end do
end subroutine



!!------ Set vdW radius and check sanity and complete vdW radius table for all involved elements
!iradtype=-1: Load radii from external file
!iradtype=0: Use UFF radii scaled by 1/1.2
!iradtype=1: Use radius specific for MK
!iradtype=2: Use radius specific for CHELPG
!espfitvdwr is in Bohr
subroutine setESPfitvdwr(iradtype,espfitvdwr)
use defvar
use util
implicit real*8 (a-h,o-z)
integer iradtype
real*8 espfitvdwr(nelesupp)
character c80*80,c200tmp*200

if (iradtype>0) then
	if (iradtype==1) then !For MK, copied from GetvdW routine (utilam)
		espfitvdwr(1:17)=(/1.20d0,1.20d0,&
		1.37d0,1.45d0,1.45d0,1.50d0,1.50d0,1.40d0,1.35d0,1.30d0,&
		! 1.57d0,1.36d0,1.24d0,1.17d0,1.80d0,1.75d0,1.70d0/) !From Gaussian code
		1.57d0,1.65d0,1.65d0,1.8d0,1.80d0,1.75d0,1.70d0/) !Modified according to my chemical intuition with consulting Bondi and UFF radii
		espfitvdwr(1:17)=espfitvdwr(1:17)/b2a !to Bohr
	else if (iradtype==2) then !For CHELPG
		espfitvdwr(1:2)=1.45D0 !vdW radius copied from GetvdW routine (utilam), some of them are given in CHELPG original paper
		espfitvdwr(3:6)=1.5D0
		espfitvdwr(7:10)=1.7D0
		espfitvdwr(11:18)=2D0
		espfitvdwr(1:18)=espfitvdwr(1:18)/b2a !to Bohr
	end if
	write(*,*) "Atomic radii used:"
	do ielem=1,nelesupp
		if (any(a%index==ielem).and.espfitvdwr(ielem)/=-1D0) write(*,"(' Element:',a,'     vdW radius (Angstrom):',f6.3)") ind2name(ielem),espfitvdwr(ielem)*b2a
	end do
	do iatm=1,ncenter
		idxtmp=a(iatm)%index
		if (espfitvdwr(idxtmp)==-1D0) then
			if (ispecial==0) then
				write(*,"(/,' vdW radius used in fitting for ',a,' is missing, input it in Angstrom, e.g. 1.7')") ind2name(idxtmp)
				write(*,"(a)") " Hint: If pressing ENTER button directly, corresponding UFF radii scaled by 1/1.2 will be used, &
                &which usually is a proper workaround. (If you want to automatically employ this treatment, setting ""ispecial"" in settings.ini to 1)"
				read(*,"(a)") c80
				if (c80==" ") then
					espfitvdwr(idxtmp)=vdwr_UFF(idxtmp)/1.2D0
					write(*,"(' Radius of ',a,' has been set to',f8.4,' Angstrom')") ind2name(idxtmp),espfitvdwr(idxtmp)*b2a
				else
					read(c80,*) tmpval
					espfitvdwr(idxtmp)=tmpval/b2a
				end if
            else if (ispecial==1) then
                espfitvdwr(idxtmp)=vdwr_UFF(idxtmp)/1.2D0
				write(*,"(/,' NOTE: vdW radius used in fitting for ',a,' is missing!')") ind2name(idxtmp)
				write(*,"(a,f6.3,a)") " Because ""ispecial"" parameter has been set to 1, therefore UFF radii of this element scaled by &
                &1/1.2 (",espfitvdwr(idxtmp)*b2a," Angstrom) is directly used, which usually is a proper workaround"
            end if
		end if
	end do
	
else if (iradtype==0) then !UFF radii
	do iatm=1,ncenter
		idxtmp=a(iatm)%index
		if (espfitvdwr(idxtmp)==-1D0) then
			espfitvdwr(idxtmp)=vdwr_UFF(idxtmp)/1.2D0
			write(*,"(' Radius of ',a,' has been set to',f8.4,' Angstrom')") ind2name(idxtmp),espfitvdwr(idxtmp)*b2a
		end if
	end do

else if (iradtype==-1) then !From external file
2737	write(*,*) "Input the path of the file containing element radii, e.g. C:\elerad.txt"
	do while(.true.)
		read(*,"(a)") c200tmp
		inquire(file=c200tmp,exist=alive)
		if (alive) exit
		write(*,*) "Cannot find the file, input again"
	end do
	open(10,file=c200tmp,status="old")
	do iatm=1,ncenter
		idxtmp=a(iatm)%index
		if (espfitvdwr(idxtmp)==-1D0) then
			call loclabel(10,ind2name(idxtmp),ifound)
			if (ifound==1) then
				read(10,*) c80,espfitvdwr(idxtmp)
				espfitvdwr(idxtmp)=espfitvdwr(idxtmp)/b2a
				write(*,"(' Radius of ',a,' has been set to',f8.4,' Angstrom')") ind2name(idxtmp),espfitvdwr(idxtmp)*b2a
			else
				write(*,"(' Error: Unable to find radius of element ',a,'!')") ind2name(idxtmp)
				close(10)
				goto 2737
			end if
		end if
	end do
	close(10)
end if
end subroutine


!!-------- Load the number of ESP fitting centers and fitting points from Gaussian pop=MK/CHELPG
subroutine loadgauESP_num(gauoutfilepath,nfitcen,nESPpt)
use util
character c80*80,gauoutfilepath*200
open(10,file=gauoutfilepath,status="old")
call loclabel(10,"NAtoms")
read(10,*) c80,nfitcen
call loclabel(10,"points will be used for")
read(10,*) nESPpt
write(*,"(' Total number of fitting points used:',i10)") nESPpt
close(10)
end subroutine
!!-------- Load ESP fitting centers, fitting points with ESP values from Gaussian pop=MK/CHELPG with IOp(6/33=2) output file
! Using these data, the result will be exactly identical to those outputted by Gaussian
! Note: For Gaussian pop=MK task, even if IOp(6/42=6) has been specified to make the density of fitting point identical to the default value of Multiwfn,
! the result of Multiwfn and Gaussian is different, because after my careful visual comparison, I found that Gaussian automatically eliminate fitting points
! in highly aggregated region. While for single atom, the number of fitting points used by Gaussian and Multiwfn are completely identical.
subroutine loadgauESP(gauoutfilepath,nfitcen,fitcen,nESPpt,ESPpt,ESPptval)
use util
use defvar
implicit real*8 (a-h,o-z)
character gauoutfilepath*200
real*8 fitcen(3,nfitcen),ESPpt(3,nESPpt),ESPptval(nESPpt)
open(10,file=gauoutfilepath,status="old")
call loclabel(10,"Atomic Center    1 is at")
do icen=1,nfitcen
	read(10,"(32x,3f10.6)") fitcen(:,icen)
end do
fitcen=fitcen/b2a
call loclabel(10,"ESP Fit Center",ifound,0)
do ipt=1,nESPpt
	read(10,"(32x,3f10.6)") ESPpt(:,ipt)
end do
ESPpt=ESPpt/b2a
call loclabel(10," Fit ",ifound,0)
do ipt=1,nESPpt
	read(10,"(14x,f10.6)") ESPptval(ipt)
end do
close(10)
end subroutine

!!---------- Calculating ESP at fitting points during calculation of ESP fitting charges
!iESPtype=1: Take nuclear charge into account  /=1: Ignore nuclear contribution
!ishowprompt=1: Show prompts  =0: Do not show
subroutine fitESP_calcESP(ishowprompt,iESPtype,nESPpt,ESPpt,ESPptval,calcfilepath)
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
integer ishowprompt,nESPpt,iESPtype
real*8 ESPpt(3,nESPpt),ESPptval(nESPpt)
character c400tmp*400,filename_tmp*200
character(len=*) calcfilepath

call walltime(iwalltime1)
if (ishowprompt==1) write(*,*) "Calculating ESP at fitting points, please wait..."
!If possible, use cubegen to calculate ESP to reduce computational time
alive=.false.
if (cubegenpath/=" ".and.ifiletype==1) then
	inquire(file=cubegenpath,exist=alive)
	if ((.not.alive).and.ishowprompt==1) then
		write(*,"(a)") " Note: Albeit current file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been defined, &
		&the cubegen cannot be found, therefore electrostatic potential will still be calculated using internal code of Multiwfn"
	end if
end if

if (alive.and.ifiletype==1) then !Use cubegen to calculate ESP
	if (ishowprompt==1) write(*,"(a)") " Since the input file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been properly defined, &
	&now Multiwfn directly invokes cubegen to calculate electrostatic potential"
	
	!Generate cubegen input file
	open(10,file="cubegenpt.txt",status="replace")
	do ipt=1,nESPpt
		write(10,"(3f16.8)") ESPpt(:,ipt)*b2a
	end do
	close(10)
	ncubegenthreads=1 !Parallel implementation prior to G16 is buggy, so test here
	if (index(cubegenpath,"G16")/=0.or.index(cubegenpath,"g16")/=0) ncubegenthreads=nthreads
    
	filename_tmp=calcfilepath
	if (index(filename,".chk")/=0) call chk2fch(filename_tmp)
	write(c400tmp,"(a,i5,a)") """"//trim(cubegenpath)//"""",ncubegenthreads," potential="//trim(cubegendenstype)//" "//&
	""""//trim(filename_tmp)//""""//" ESPresult.cub -5 h < cubegenpt.txt > nouseout"
    call runcommand(c400tmp)
	if (index(filename,".chk")/=0) call delfile(filename_tmp)
	
	!Load ESP data from cubegen resulting file
	open(10,file="ESPresult.cub",status="old")
	do iskip=1,6+ncenter
		read(10,*)
	end do
	do ipt=1,nESPpt
		read(10,*) xtmp,ytmp,ztmp,ESPptval(ipt)
		if (iESPtype==2.or.iESPtype==3) ESPptval(ipt)=ESPptval(ipt)-nucesp(xtmp/b2a,ytmp/b2a,ztmp/b2a) !Remove nuclear contribution
	end do
	close(10)
    
	!Delete intermediate files
    call delfile("cubegenpt.txt ESPresult.cub nouseout")
    
else !Use internal code to evaluate ESP
    nESPthreads=nthreads
    if (iESPcode==2.or.iESPcode==3) then
        call doinitlibreta(1)
        if (isys==1.and.nESPthreads>12) nESPthreads=12
    end if
    write(*,*)
    ifinish=0;ishowprog=1
    call showprog(0,nESPpt)
    ntmp=floor(nESPpt/100D0)
	!$OMP PARALLEL DO SHARED(ifinish,ESPptval,ishowprog) PRIVATE(ipt) schedule(dynamic) NUM_THREADS(nESPthreads)
	do ipt=1,nESPpt
		if (iESPtype==1) then !Take nuclear charge into account
			ESPptval(ipt)=totesp(ESPpt(1,ipt),ESPpt(2,ipt),ESPpt(3,ipt))
		else if (iESPtype==2.or.iESPtype==3) then !Do not take nuclear charge into account
			ESPptval(ipt)=eleesp(ESPpt(1,ipt),ESPpt(2,ipt),ESPpt(3,ipt))
		end if
		if (ntmp/=0) then
			!$OMP CRITICAL
			ifinish=ifinish+1
			ishowprog=mod(ifinish,ntmp)
			if (ishowprog==0) call showprog(floor(100D0*ifinish/nESPpt),100)
			!$OMP END CRITICAL
        end if
	end do    !$OMP END PARALLEL DO
    if (ishowprog/=0) call showprog(100,100)
end if
call walltime(iwalltime2)
if (ishowprompt==1) write(*,"(' Calculation of ESP took up wall clock time',i10,' s')") iwalltime2-iwalltime1
end subroutine

!!--------- Generate numpt points scattered evenly on an unit sphere, used by e.g. module of obtaining MK charge
! Input argument numpt is the expected number of points, while the return value is actual number
! ptcrd store coordinates of the points 
subroutine unitspherept(ptcrd,numpt)
implicit real*8 (a-h,o-z)
real*8 ptcrd(3,numpt)
integer numpt
pi=3.141592653589793D0
!The average number of equator points in all XY layes is numequ*2/pi, and there are numvert=numequ/2 layers
!Solve (numequ*2/pi)*numequ/2=numpt one can get numequ=sqrt(numpt*pi)
numequ=int(sqrt(numpt*pi)) !Maximal number of point in each XY layer
numvert=numequ/2
ipt=0
do ivert=0,numvert
	angz=dfloat(ivert)/numvert*pi
	scalexy=sin(angz)
	z=cos(angz)
	numxy=int(numequ*scalexy)
	if (numxy==0) numxy=1
	do ihori=1,numxy
		ipt=ipt+1
		if (ipt>numpt) then
			numpt=ipt-1
			return
		end if
		angxy=2D0*pi*ihori/numxy
		ptcrd(1,ipt)=cos(angxy)*scalexy
		ptcrd(2,ipt)=sin(angxy)*scalexy
		ptcrd(3,ipt)=z
	end do
end do
numpt=ipt
end subroutine











!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!!============================ Hirshfeld-I ============================!!
!Wrapper of Hirshfeld-I module to automatically set radpot and sphpot to proper values
subroutine Hirshfeld_I_wrapper(itype)
use defvar
implicit real*8 (a-h,o-z)
nradpotold=radpot
nsphpotold=sphpot
if (iautointgrid==1) then
	radpot=30
	sphpot=170
	if (any(a%index>18)) radpot=40
	if (any(a%index>36)) radpot=50
	if (any(a%index>54)) radpot=60
end if
call Hirshfeld_I(itype)
if (iautointgrid==1) then
	radpot=nradpotold
	sphpot=nsphpotold
end if
end subroutine

!!========== Calculate Hirshfeld-I charge or atomic radial density
!I've compared this module with hipart, this module is faster than hipart, and the accuracy under default setting is at least never lower than hipart
!Note: Use "subroutine Hirshfeld_I_evengrid" for periodic systems
!itype: 1=Calculate and print charges =2: Only generate atomic spaces, namely filling "atmraddens" global array by final radial density of each atom
subroutine Hirshfeld_I(itype)
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
integer itype
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 molrho(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),selfdens(radpot*sphpot),molrhoall(ncenter,radpot*sphpot)
real*8 charge(ncenter),lastcharge(ncenter) !Atomic charge of current iter. and last iter.
real*8 radrholow(200),radrhohigh(200)
character sep,c80tmp*80
character(len=2) :: statname(-4:4)=(/ "-4","-3","-2","-1","_0","+1","+2","+3","+4" /)
integer :: maxcyc=50,ioutmedchg=0
real*8 :: crit=0.0002D0
!Used for mode 2. e.g. atmstatgrid(iatm,igrid,jatm,-1) means density of jatm with -1 charge state at igrid around iatm
real*8,allocatable :: atmstatgrid(:,:,:,:)

!Ignore jatm contribution to iatm centered grids if distance between iatm and jatm is larger than 1.5 times of sum of their vdwr
!This can reduce lots of time for large system, the lose of accuracy can be ignored (error is ~0.0001 per atom)
integer :: ignorefar=1
real*8 :: vdwsumcut=2D0

!Mode 1 use very low memory but expensive, because most data is computed every iteration
!Mode 2 use large memory but fast, because most data is only computed once at initial stage
!The result of the two modes differ with each other marginally, probably because in mode 1 radial density is related to max(npthigh,nptlow), which is not involved in mode 2
!In principle, result of mode 2 is slightly better
integer :: imode=2

ntotpot=radpot*sphpot
write(*,"(/,a,/)") " IMPORTANT HINT: If your system does not contain lanthanides and actinides and meantime you want to directly &
&perform Hirshfeld-I calculation without letting Multiwfn to automatically invoke Gaussian &
&to generate atomic .wfn files for various charged states, it is strongly suggested to copy ""atmrad"" folder from ""examples"" directory &
&to current directory, then the atomic radial density files in the ""atmrad"" folder will be directly utilized. See Section &
&3.9.13 of Multiwfn manual for more detail about the underlying mechanism, and see Section 4.7.4 for example of Hirshfeld-I calculation"

do while(.true.)
    if (itype==1) write(*,*) "     =============== Iterative Hirshfeld (Hirshfeld-I) ==============="
    if (itype==2) write(*,*) "     ============== Generate Hirshfeld-I atomic weights =============="
	if (ignorefar==1) write(*,"(a,f6.3)") " -3 Switch if speeding up calculation using distance cutoff, current: Yes, ratio factor is",vdwsumcut
	if (ignorefar==0) write(*,*) "-3 Switch if speeding up calculation using distance cutoff, current: No"
	if (imode==1) write(*,*) "-2 Switch algorithm, current: Slow & low memory requirement"
	if (imode==2) write(*,*) "-2 Switch algorithm, current: Fast & large memory requirement"
	if (itype==1) then
		if (ioutmedchg==0) write(*,*) "-1 Switch if outputting intermediate results, current: No"
		if (ioutmedchg==1) write(*,*) "-1 Switch if outputting intermediate results, current: Yes"
		write(*,*) "0 Return"
	end if
	write(*,*) "1 Start calculation!"
	write(*,"(a,i4)") " 2 Set the maximum number of iterations, current:",maxcyc
	write(*,"(a,f10.6)") " 3 Set convergence criterion of atomic charges, current:",crit
	read(*,*) isel
	if (isel==-3) then
        if (ignorefar==1) then
            ignorefar=0
        else
            ignorefar=1
            write(*,*) "Input ratio factor of cutoff, e.g. 2.5"
            write(*,*) "Note: The higher the value, the more accurate the result and the more robust &
            &the calculation will be, however the computational cost will be correspondingly higher. The default value is 2.0"
            read(*,*) vdwsumcut
        end if
	else if (isel==-2) then
		if (imode==1) then
			imode=2
		else
			imode=1
			crit=0.001 !mode 1 is more time-consuming, use loose criterion
		end if
	else if (isel==-1) then
		if (ioutmedchg==1) then
			ioutmedchg=0
		else
			ioutmedchg=1
		end if
	else if (isel==0) then
		return
	else if (isel==1) then
		exit
	else if (isel==2) then
		write(*,*) "Input maximum number of iterations, e.g. 30"
		read(*,*) maxcyc
	else if (isel==3) then
		write(*,*) "Input convergence criterion of atomic charges, e.g. 0.001"
		read(*,*) crit
	end if
end do

!Generate all needed .rad files if not provided
call genatmradfile

call walltime(iwalltime1)

!Generate single center integration grids
call gen1cintgrid(gridatmorg,iradcut)
write(*,"(' Radial grids:',i5,'  Angular grids:',i5,'  Total:',i7,'  After pruning:',i7)") radpot,sphpot,radpot*sphpot,radpot*sphpot-iradcut*sphpot

!Calculate molecular density
write(*,*) "Calculating actual density of system at all grids..."
do iatm=1,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
    !$OMP parallel do shared(molrhoall) private(ipt) num_threads(nthreads)
	do ipt=1+iradcut*sphpot,ntotpot
		molrhoall(iatm,ipt)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
	end do
    !$OMP end parallel do
    call showprog(iatm,ncenter)
end do

if (allocated(atmradnpt)) deallocate(atmradnpt)
if (allocated(atmraddens)) deallocate(atmraddens)
allocate(atmradnpt(ncenter),atmraddens(200,ncenter))
sep='/' !Separation symbol of directory
if (isys==1) sep='\'

ichgmin=-3 !Actually allowed charged range
ichgmax=3
!Calculate contribution of all atoms in every state to each atomic centered grids
if (imode==2) then
	allocate(atmstatgrid(ncenter,ntotpot,ncenter,ichgmin:ichgmax))
	atmstatgrid=0
    rmem=dfloat(ncenter)**2*dfloat(ntotpot)*(ichgmax-ichgmin+1)*8/1024D0/1024D0 !Must change integer to float, otherwise intermediate integer may exceed upper limit of recording
    write(*,"(a,f10.1,a)") " Memory requirement for storing atomic densities on grids:",rmem," MB"
	write(*,*) "Calculating atomic density contribution to grids..."
	do iatm=1,ncenter !The center of grids
		gridatm%value=gridatmorg%value !Weight in this grid point
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		do istat=ichgmin,ichgmax !Charge state
			do jatm=1,ncenter
				if (ignorefar==1) then
					if (atomdist(iatm,jatm,1)>(vdwr(a(iatm)%index)+vdwr(a(jatm)%index))*vdwsumcut) cycle
				end if
				if (a(jatm)%index==1.and.istat==1) cycle !H+ doesn't contains electron and cannot compute density
				c80tmp="atmrad"//sep//trim(a(jatm)%name)//statname(istat)//".rad"
				inquire(file=c80tmp,exist=alive)
				if (.not.alive) cycle !If the .rad file of jatm in charge state "istat" is not available, skip calculation
				open(10,file=c80tmp,status="old")
				read(10,*) atmradnpt(jatm)
				do ipt=1,atmradnpt(jatm)
					read(10,*) rnouse,atmraddens(ipt,jatm)
				end do
				close(10)
                !I have made great effort to try to parallelize this part, however after doing this, the speed is even significantly lowered
                !So I decided not to parallelize it
				do ipt=1+iradcut*sphpot,ntotpot
					atmstatgrid(iatm,ipt,jatm,istat)=fdens_rad(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
				end do
			end do
		end do
        call showprog(iatm,ncenter)
	end do
end if

!Set atomic initial radial density as neutral state, which is loaded from corresponding .rad file
atmraddens=0
do iatm=1,ncenter
    c80tmp="atmrad"//sep//trim(a(iatm)%name)//"_0.rad"
    inquire(file=c80tmp,exist=alive)
    if (.not.alive) then
        write(*,*)
        write(*,*) "Error: The file "//trim(c80tmp)//" cannot be found!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
	open(10,file=c80tmp,status="old")
	read(10,*) atmradnpt(iatm)
	do ipt=1,atmradnpt(iatm)
		read(10,*) rnouse,atmraddens(ipt,iatm)
	end do
	close(10)
end do

write(*,*)
write(*,*) "Performing Hirshfeld-I iteration to refine atomic spaces..."
lastcharge=0
!Cycle each atom to calculate their charges
do icyc=1,maxcyc
	if (ioutmedchg==1) write(*,*)
	if (icyc==1) then
		write(*,"(' Cycle',i5)") icyc
	else
		write(*,"(' Cycle',i5,'   Maximum change:',f10.6)") icyc,varmax
	end if
	
	do iatm=1,ncenter
		gridatm%value=gridatmorg%value !Weight in this grid point
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		
		!Molecular density
		molrho=molrhoall(iatm,:)
		
		!Calculate promolecular and proatomic density 
		promol=0D0
		do jatm=1,ncenter
			if (ignorefar==1) then
				if (atomdist(iatm,jatm,1)>(vdwr(a(iatm)%index)+vdwr(a(jatm)%index))*vdwsumcut) cycle
			end if
			if (imode==1) then
				!$OMP parallel do shared(tmpdens) private(ipt) schedule(dynamic) num_threads(nthreads)
				do ipt=1+iradcut*sphpot,ntotpot
					tmpdens(ipt)=fdens_rad(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
				end do
				!$OMP end parallel do
			else if (imode==2) then
				if (icyc==1) then
					tmpdens=atmstatgrid(iatm,:,jatm,0)
				else
					ichglow=floor(lastcharge(jatm))	
					ichghigh=ceiling(lastcharge(jatm))
					tmpdens=(lastcharge(jatm)-ichglow)*atmstatgrid(iatm,:,jatm,ichghigh) + (ichghigh-lastcharge(jatm))*atmstatgrid(iatm,:,jatm,ichglow)
				end if
			end if
			promol=promol+tmpdens
			if (jatm==iatm) selfdens=tmpdens
		end do
		
		!Calculate atomic charge
		electmp=0D0
		do ipt=1+iradcut*sphpot,ntotpot
			if (promol(ipt)/=0D0) electmp=electmp+selfdens(ipt)/promol(ipt)*molrho(ipt)*gridatm(ipt)%value
		end do
		if (nEDFelec==0) then
			charge(iatm)=a(iatm)%charge-electmp
		else !EDF is used for some atoms. Core electron density represented by EDF has been integrated, so nuclear charge should be augmented by nEDFelecatm
			charge(iatm)=a(iatm)%charge+nEDFelecatm(iatm)-electmp
        end if
		if (ioutmedchg==1) write(*,"(' Charge of atom',i5,'(',a2,')',': ',f12.6,'  Delta:',f12.6)") &
		iatm,a(iatm)%name,charge(iatm),charge(iatm)-lastcharge(iatm)
	end do
	
	!Check convergence
	varmax=maxval(abs(charge-lastcharge))
	if (varmax<crit) then
		if (itype==1) then
			write(*,"(a,f10.6)") " All atomic charges have converged to criterion of",crit
			write(*,"(' Sum of all charges:',f14.8)") sum(charge)
            !Calculate and print normalized charge
			call normalize_atmchg(charge(:))
            call printatmchg(charge(:))
			exit
		else
			write(*,*) "Hirshfeld-I atomic spaces converged successfully!"
			return
		end if
	else
		if (icyc==maxcyc) then
			write(*,"(/,' Convergence failed within',i4,' cycles!')") maxcyc
			exit
		end if
	end if
	
	!Update atomic radial density by means of interpolation of adjacent charge state
	do iatm=1,ncenter
		!Read radial density of lower limit state
		ichglow=floor(charge(iatm))
		radrholow=0
		c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(ichglow)//".rad"
		inquire(file=c80tmp,exist=alive)
		if (.not.alive) then
			write(*,"(' Error: ',a,' is needed but was not prepared!')") trim(c80tmp)
            write(*,"(' Current charge of atom',i5,'(',a,'):',f12.8)") iatm,a(iatm)%name,charge(iatm)
            write(*,"(a)") " Note: This error implies that this atom has unusual charge. You should manually provide the corresponding .rad file &
            &in ""atmrad"" folder prior to the calculation. See Section 3.9.13 of Multiwfn manual for detail."
			return
		end if
		open(10,file=c80tmp,status="old")
		read(10,*) nptlow
		do ipt=1,nptlow
			read(10,*) rnouse,radrholow(ipt)
		end do
		close(10)
		!Read radial density of upper limit state
		ichghigh=ceiling(charge(iatm))
		radrhohigh=0
		c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(ichghigh)//".rad"
		inquire(file=c80tmp,exist=alive)
		if (.not.alive) then
			write(*,"(' Error: ',a,' is needed but was not prepared!')") trim(c80tmp)
            write(*,"(' Current charge of atom',i5,'(',a,'):',f12.8)") iatm,a(iatm)%name,charge(iatm)
            write(*,"(a)") " Note: This error implies that this atom has unusual charge. You should manually provide the corresponding .rad file &
            &in ""atmrad"" folder prior to the calculation. See Section 3.9.13 of Multiwfn manual for detail."
			return
		end if
		open(10,file=c80tmp,status="old")
		read(10,*) npthigh
		do ipt=1,npthigh
			read(10,*) rnouse,radrhohigh(ipt)
		end do
		close(10)
		!Update current radial density
		atmraddens(:,iatm)=(charge(iatm)-ichglow)*radrhohigh(:) + (ichghigh-charge(iatm))*radrholow(:)
		atmradnpt(iatm)=max(npthigh,nptlow)
	end do
	
	lastcharge=charge
end do

if (allocated(frag1)) then
    write(*,"(/,' Fragment charge:',f14.8)") sum(charge(frag1))
    write(*,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(charge(frag1))
end if
call walltime(iwalltime2)
write(*,*)
write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

call outatmchg(10,charge(:))
end subroutine





!!========== Calculate Hirshfeld-I atomic charges or atomic radial density, using evenly distributed grid, mainly for GPW periodic wavefunctions
!itype: 1=Calculate and print charges =2: Only generate atomic spaces, namely filling "atmraddens" global array by final radial density of each atom
!imode: 1=Calculate actual density from periodic wavefunction    2=Actual density is directly taken from grid data in memory
!Note: Use "subroutine Hirshfeld_I" for isolated systems
subroutine Hirshfeld_I_evengrid(itype,imode)
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
integer itype,imode
real*8 charge(ncenter),lastcharge(ncenter) !Atomic charge of current iter. and last iter.
real*8 radrholow(200),radrhohigh(200)
real*8 tvec(3),atmrho(ncenter),atmpop(ncenter),atmpop_tmp(ncenter)
character sep,c80tmp*80
character(len=2) :: statname(-4:4)=(/ "-4","-3","-2","-1","_0","+1","+2","+3","+4" /)
integer :: maxcyc=50,ioutmedchg=0
real*8 :: crit=0.0002D0

write(*,"(/,a,/)") " IMPORTANT HINT: If your system does not contain lanthanides and actinides and meantime you want to directly &
&perform Hirshfeld-I calculation without letting Multiwfn to automatically invoke Gaussian &
&to generate atomic .wfn files for various charged states, it is strongly suggested to copy ""atmrad"" folder from ""examples"" directory &
&to current directory, then the atomic radial density files in the ""atmrad"" folder will be directly utilized. See Section &
&3.9.13 of Multiwfn manual for more detail about the underlying mechanism, and see Section 4.7.4 for example of Hirshfeld-I calculation"

do while(.true.)
    if (itype==1) write(*,*) "     =============== Iterative Hirshfeld (Hirshfeld-I) ==============="
    if (itype==2) write(*,*) "     ============== Generate Hirshfeld-I atomic weights =============="
	if (itype==1) then
		if (ioutmedchg==0) write(*,*) "-1 Switch if outputting intermediate results, current: No"
		if (ioutmedchg==1) write(*,*) "-1 Switch if outputting intermediate results, current: Yes"
		write(*,*) "0 Return"
	end if
	write(*,*) "1 Start calculation!"
	write(*,"(a,i4)") " 2 Set the maximum number of iterations, current:",maxcyc
	write(*,"(a,f10.6)") " 3 Set convergence criterion of atomic charges, current:",crit
	read(*,*) isel
	if (isel==-1) then
		if (ioutmedchg==1) then
			ioutmedchg=0
		else
			ioutmedchg=1
		end if
	else if (isel==0) then
		return
	else if (isel==1) then
		exit
	else if (isel==2) then
		write(*,*) "Input maximum number of iterations, e.g. 30"
		read(*,*) maxcyc
	else if (isel==3) then
		write(*,*) "Input convergence criterion of atomic charges, e.g. 0.001"
		read(*,*) crit
	end if
end do

!Generate all needed .rad files if not provided
call genatmradfile

!Prepare density of the actual system
if (imode==1) then !Calculate density from periodic wavefunction
    call setgrid_for_PBC(0.2D0,1)
    if (allocated(cubmat)) deallocate(cubmat)
    allocate(cubmat(nx,ny,nz))
    call walltime(iwalltime1)
	!Because uniform grid cannot integrate well core density, so temporarily disable EDFs
    nEDFprims_org=nEDFprims
    nEDFprims=0
    call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
    write(*,*) "Calculating electron density grid data..."
    call savecubmat(1,0,1)
    call delvirorb_back(1) !Restore to previous wavefunction
    nEDFprims=nEDFprims_org
else !Directly using loaded electron density from cub/VASP grid data, and transforming grid data information to cell information
	if (all(a%charge==0)) then
		write(*,*) "Error: All nuclear charges are zero! If this file was exported by CP2K, it is a bug. You need to manually &
        &edit the file so that effective nuclear charges (column 2 since line 8) are correctly recorded, otherwise atomic charges cannot be calculated"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
	end if
    call grid2cellinfo
    !call showcellinfo
    call walltime(iwalltime1)
end if

call calc_dvol(dvol)
if (allocated(atmradnpt)) deallocate(atmradnpt)
if (allocated(atmraddens)) deallocate(atmraddens)
allocate(atmradnpt(ncenter),atmraddens(200,ncenter))
sep='/' !Separation symbol of directory
if (isys==1) sep='\'
ichgmin=-3 !Actually allowed charged range
ichgmax=3

!Set atomic initial radial density as neutral state, which is loaded from corresponding .rad file
atmraddens=0
do iatm=1,ncenter
    c80tmp="atmrad"//sep//trim(a(iatm)%name)//"_0.rad"
    inquire(file=c80tmp,exist=alive)
    if (.not.alive) then
        write(*,*)
        write(*,*) "Error: The file "//trim(c80tmp)//" cannot be found!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
	open(10,file=c80tmp,status="old")
	read(10,*) atmradnpt(iatm)
	do ipt=1,atmradnpt(iatm)
		read(10,*) rnouse,atmraddens(ipt,iatm)
	end do
	close(10)
end do

write(*,*)
write(*,*) "Performing Hirshfeld-I iteration to refine atomic spaces..."
lastcharge=0
do icyc=1,maxcyc
	if (ioutmedchg==1) write(*,*)
	if (icyc==1) then
		write(*,"(' Cycle',i5)") icyc
	else
		write(*,"(' Cycle',i5,'   Maximum change:',f10.6)") icyc,varmax
	end if

	atmpop(:)=0
	ifinish=0;ishowprog=1
	ntmp=floor(ny*nz/100D0)
	!$OMP PARALLEL SHARED(atmpop,ifinish,ishowprog) PRIVATE(i,j,k,tmpx,tmpy,tmpz,iatm,npt,atmrho,prorho,atmpop_tmp,&
    !$OMP ic,jc,kc,icell,jcell,kcell,tvec,atmx,atmy,atmz,dist2,tmprho) NUM_THREADS(nthreads)
	atmpop_tmp(:)=0
	!$OMP DO schedule(dynamic) collapse(2)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				if (abs(cubmat(i,j,k))<1D-11) cycle !Note that the electron density around core produced by VASP PAW calculation can be negative, so use abs()
				call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
				atmrho(:)=0
				call getpointcell(tmpx,tmpy,tmpz,ic,jc,kc)
				do icell=ic-PBCnx,ic+PBCnx
					do jcell=jc-PBCny,jc+PBCny
						do kcell=kc-PBCnz,kc+PBCnz
							call tvec_PBC(icell,jcell,kcell,tvec)
							do iatm=1,ncenter
								atmx=a(iatm)%x+tvec(1)
								atmy=a(iatm)%y+tvec(2)
								atmz=a(iatm)%z+tvec(3)
								dist2=(atmx-tmpx)**2+(atmy-tmpy)**2+(atmz-tmpz)**2
								if (dist2>atmrhocutsqr(a(iatm)%index)) then
									cycle
								else
									npt=atmradnpt(iatm)
									call lagintpol(atmradpos(1:npt),atmraddens(1:npt,iatm),npt,dsqrt(dist2),tmprho,rnouse,rnouse,1)
									atmrho(iatm)=atmrho(iatm)+tmprho
								end if
							end do
						end do
					end do
				end do
				prorho=sum(atmrho(:))
                if (prorho>0) atmpop_tmp(:)=atmpop_tmp(:)+atmrho(:)/prorho*cubmat(i,j,k)*dvol
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
	atmpop(:)=atmpop(:)+atmpop_tmp(:)
	!$OMP END CRITICAL
	!$OMP END PARALLEL
	if (ishowprog/=0) call showprog(100,100)

	!Evaluate current charge. Note that EDFs were not involved in evaluating system density
    charge(:)=a(:)%charge-atmpop(:)
	if (ioutmedchg==1) then
		do iatm=1,ncenter
			write(*,"(' Charge of atom',i5,'(',a2,')',': ',f12.6,'  Delta:',f12.6)") iatm,a(iatm)%name,charge(iatm),charge(iatm)-lastcharge(iatm)
        end do
    end if

	!Check convergence
	varmax=maxval(abs(charge(:)-lastcharge(:)))
	if (varmax<crit.or.icyc==maxcyc) then
        if (varmax<crit) write(*,"(/,a,f10.6)") " All atomic charges have converged to criterion of",crit
        if (icyc==maxcyc) write(*,"(/,' Convergence failed within',i4,' cycles!')") maxcyc
		if (itype==2) then
			write(*,*) "Construction of Hirshfeld-I atomic spaces has finished"
			return
        end if
		exit
	end if
	
	!Update atomic radial density by means of interpolation of adjacent charge state
	do iatm=1,ncenter
		!Read radial density of lower limit state
		ichglow=floor(charge(iatm))
		radrholow=0
		c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(ichglow)//".rad"
		inquire(file=c80tmp,exist=alive)
		if (.not.alive) then
			write(*,"(' Error: ',a,' is needed but was not prepared!')") trim(c80tmp)
            write(*,"(' Current charge of atom',i5,'(',a,'):',f12.8)") iatm,a(iatm)%name,charge(iatm)
            write(*,"(a)") " Note: This error implies that this atom has unusual charge. You should manually provide the corresponding .rad file &
            &in ""atmrad"" folder prior to the calculation. See Section 3.9.13 of Multiwfn manual for detail."
			return
		end if
		open(10,file=c80tmp,status="old")
		read(10,*) nptlow
		do ipt=1,nptlow
			read(10,*) rnouse,radrholow(ipt)
		end do
		close(10)
		!Read radial density of upper limit state
		ichghigh=ceiling(charge(iatm))
		radrhohigh=0
		c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(ichghigh)//".rad"
		inquire(file=c80tmp,exist=alive)
		if (.not.alive) then
			write(*,"(' Error: ',a,' is needed but was not prepared!')") trim(c80tmp)
            write(*,"(' Current charge of atom',i5,'(',a,'):',f12.8)") iatm,a(iatm)%name,charge(iatm)
            write(*,"(a)") " Note: This error implies that this atom has unusual charge. You should manually provide the corresponding .rad file &
            &in ""atmrad"" folder prior to the calculation. See Section 3.9.13 of Multiwfn manual for detail."
			return
		end if
		open(10,file=c80tmp,status="old")
		read(10,*) npthigh
		do ipt=1,npthigh
			read(10,*) rnouse,radrhohigh(ipt)
		end do
		close(10)
		!Update current radial density
		atmraddens(:,iatm)=(charge(iatm)-ichglow)*radrhohigh(:) + (ichghigh-charge(iatm))*radrholow(:)
		atmradnpt(iatm)=max(npthigh,nptlow)
	end do
	
    !Update atomic charges
	lastcharge(:)=charge(:)
end do

write(*,"(' Sum of all raw charges:',f14.8)") sum(charge(:))
!Normalize atomic charges. This is not feasible if only grid data is available, &
!because in this case the nelec used in "normalize_atmchg" is simply guessed by assuming system is neutral
if (imode==1) call normalize_atmchg(charge(:))
!Print final atomic charges
call printatmchg(charge(:))

if (allocated(frag1)) then
	write(*,"(/,' Fragment charge:',f14.8)") sum(charge(frag1))
	write(*,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(charge(frag1))
end if

call walltime(iwalltime2)
write(*,"(/,' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

call outatmchg(10,charge(:))
end subroutine





!!------- Generate atomic radial density files at different states, used for e.g. Hirshfeld-I
!"atmrad" in current folder is used as working directory
!Various charge states of elements in current system will be calculated to produce atomic .wfn file by Gaussian, then radial density file (.rad) will be generated
!If needed atomic wfn file is already existed, calculation will be skipped
!Radial distance values are the same as built-in atomic density, i.e. those in atmraddens.f90
subroutine genatmradfile
use defvar
use util
implicit real*8 (a-h,o-z)
character c80tmp*80,c200tmp*200,calclevel*80,radname*200,sep
character(len=2) :: statname(-3:3)=(/ "-3","-2","-1","_0","+1","+2","+3" /)
integer :: chgmulti(nelesupp,-3:3)=0 !Ground state multiplicity of each charge state of each element. If value=0, means undefined

!Define chgmulti for elements for possible states
!Charge states unlikely involved in Hirshfeld-I calculation will not be defined here and calculated later
!H,Li,Na,K,Rb,Cs
chgmulti(1,0)=2
chgmulti(1,1)=1
chgmulti(1,-1)=1
chgmulti(3,:)=chgmulti(1,:)
chgmulti(11,:)=chgmulti(1,:)
chgmulti(19,:)=chgmulti(1,:)
chgmulti(37,:)=chgmulti(1,:)
chgmulti(55,:)=chgmulti(1,:)
!He,Ne,Ar,Kr,Xe,Rn
chgmulti(2,0)=1
chgmulti(2,1)=2
chgmulti(2,-1)=2
chgmulti(10,:)=chgmulti(2,:)
chgmulti(18,:)=chgmulti(2,:)
chgmulti(36,:)=chgmulti(2,:)
chgmulti(54,:)=chgmulti(2,:)
chgmulti(86,:)=chgmulti(2,:)
!Be,Mg,Ca,Sr,Ba
chgmulti(4,0)=1
chgmulti(4,1)=2
chgmulti(4,2)=1
chgmulti(4,-1)=2
chgmulti(12,:)=chgmulti(4,:)
chgmulti(20,:)=chgmulti(4,:)
chgmulti(38,:)=chgmulti(4,:)
chgmulti(56,:)=chgmulti(4,:)
!B,Al,Ga,In,Tl
chgmulti(5,0)=2
chgmulti(5,1)=1
chgmulti(5,2)=2
chgmulti(5,-1)=3
chgmulti(5,-2)=4
chgmulti(13,:)=chgmulti(5,:)
chgmulti(31,:)=chgmulti(5,:)
chgmulti(49,:)=chgmulti(5,:)
chgmulti(81,:)=chgmulti(5,:)
!C,Si,Ge,Sn,Pb
chgmulti(6,0)=3
chgmulti(6,1)=2
chgmulti(6,2)=1
chgmulti(6,-1)=4
chgmulti(6,-2)=3
chgmulti(14,:)=chgmulti(6,:)
chgmulti(32,:)=chgmulti(6,:)
chgmulti(50,:)=chgmulti(6,:)
chgmulti(82,:)=chgmulti(6,:)
!N,P,As,Sb,Bi
chgmulti(7,0)=4
chgmulti(7,1)=3
chgmulti(7,2)=2
chgmulti(7,-1)=3
chgmulti(7,-2)=2
chgmulti(15,:)=chgmulti(7,:)
chgmulti(33,:)=chgmulti(7,:)
chgmulti(51,:)=chgmulti(7,:)
chgmulti(83,:)=chgmulti(7,:)
!O,S,Se,Te,Po
chgmulti(8,0)=3
chgmulti(8,1)=4
chgmulti(8,2)=3
chgmulti(8,-1)=2
chgmulti(8,-2)=1
chgmulti(16,:)=chgmulti(8,:)
chgmulti(34,:)=chgmulti(8,:)
chgmulti(52,:)=chgmulti(8,:)
chgmulti(84,:)=chgmulti(8,:)
!F,Cl,Br,I,At
chgmulti(9,0)=2
chgmulti(9,1)=3
chgmulti(9,2)=4
chgmulti(9,-1)=1
chgmulti(17,:)=chgmulti(9,:)
chgmulti(35,:)=chgmulti(9,:)
chgmulti(53,:)=chgmulti(9,:)
chgmulti(85,:)=chgmulti(9,:)
!Spin multiplicity of transition metal for each state is determined by chemical intuition as well as a few single point energy data
!For simplicity, I assume that later elements in each row has identical configuration, of course this is not always correct but not too bad
!Sc (3d1,4s2)
chgmulti(21,0)=2
chgmulti(21,1)=3
chgmulti(21,2)=2
chgmulti(21,-1)=3
chgmulti(39,:)=chgmulti(21,:) !Y
chgmulti(57,:)=chgmulti(21,:) !La
!Ti (3d2,4s2)
chgmulti(22,0)=3
chgmulti(22,1)=4
chgmulti(22,2)=3
chgmulti(22,-1)=4
chgmulti(40,:)=chgmulti(22,:) !Zr
chgmulti(72,:)=chgmulti(22,:) !Hf
!V  (3d3,4s2)
chgmulti(23,0)=4
chgmulti(23,1)=5
chgmulti(23,2)=4
chgmulti(23,-1)=5
chgmulti(41,:)=chgmulti(23,:) !Nb
chgmulti(73,:)=chgmulti(23,:) !Ta
!Cr (3d5,4s1)
chgmulti(24,0)=7
chgmulti(24,1)=6
chgmulti(24,2)=5
chgmulti(24,-1)=6
chgmulti(42,:)=chgmulti(24,:) !Mo
chgmulti(74,:)=chgmulti(24,:) !W
!Mn (3d5,4s2)
chgmulti(25,0)=6
chgmulti(25,1)=7
chgmulti(25,2)=6
chgmulti(25,-1)=5
chgmulti(43,:)=chgmulti(25,:) !Tc
chgmulti(75,:)=chgmulti(25,:) !Re
!Fe (3d6,4s2)
chgmulti(26,0)=5
chgmulti(26,1)=6
chgmulti(26,2)=5
chgmulti(26,-1)=4
chgmulti(44,:)=chgmulti(26,:) !Ru
chgmulti(76,:)=chgmulti(26,:) !Os
!Co (3d7,4s2)
chgmulti(27,0)=4
chgmulti(27,1)=5
chgmulti(27,2)=4
chgmulti(27,-1)=3
chgmulti(45,:)=chgmulti(27,:) !Rh
chgmulti(77,:)=chgmulti(27,:) !Ir
!Ni (3d8,4s2)
chgmulti(28,0)=3
chgmulti(28,1)=4
chgmulti(28,2)=3
chgmulti(28,-1)=2
chgmulti(46,:)=chgmulti(28,:) !Pd
chgmulti(78,:)=chgmulti(28,:) !Pt
!Cu (3d10,4s1)
chgmulti(29,0)=2
chgmulti(29,1)=1
chgmulti(29,2)=2
chgmulti(29,-1)=1
chgmulti(47,:)=chgmulti(29,:) !Ag
chgmulti(79,:)=chgmulti(29,:) !Au
!Zn (3d10,4s2)
chgmulti(30,0)=1
chgmulti(30,1)=2
chgmulti(30,2)=1
chgmulti(30,-1)=2
chgmulti(48,:)=chgmulti(30,:) !Cd
chgmulti(80,:)=chgmulti(30,:) !Hg

sep='/' !Separation symbol of directory
if (isys==1) sep='\'
calclevel=" "
icalcatmwfn=0

!Cycle each charge state of each atom. Each element is only calculated once. If the .rad file is already existent, don't calculate again
do iatm=1,ncenter
	iele=a(iatm)%index
	do istat=-3,3
		if (chgmulti(iele,istat)==0) cycle !Undefined state
		radname="atmrad"//sep//trim(a(iatm)%name)//statname(istat)//".rad"
		inquire(file=radname,exist=alive)
		if (alive) cycle
		radname="atmrad"//sep//trim(a(iatm)%name)//statname(istat)//".wfn"
		inquire(file=radname,exist=alive)
		if (alive) cycle
		
        icalcatmwfn=1
		!Check Gaussian path
		inquire(file=gaupath,exist=alive)
		if (.not.alive) then
			write(*,*) "Could not find Gaussian path defined in ""gaupath"" variable in settings.ini"
			if (isys==1) write(*,*) "Input the path of Gaussian executable file, e.g. ""D:\study\g09w\g09.exe"""
			if (isys==2) write(*,*) "Input the path of Gaussian executable file, e.g. ""/sob/g09/g09"""
			do while(.true.)
				read(*,"(a)") gaupath
				inquire(file=gaupath,exist=alive)
				if (alive) exit
				write(*,*) "Could not find Gaussian executable file, input again"
			end do
		end if
		
		!Input calculation level
		if (calclevel==" ") then
			write(*,*) "Some atomic .wfn files are not found in ""atmrad"" folder in current directory"
			write(*,"(a)") " Now input the level for calculating these .wfn files, e.g. B3LYP/def2SVP"
			write(*,"(a)") " You can also add other keywords at the same time, e.g. M062X/6-311G(2df,2p) scf=xqc int=ultrafine"
			read(*,"(a)") calclevel
		end if
		
		!Generate .gjf file 
		call inquire_dir("atmrad",alive)
		if (.not.alive) call system("mkdir atmrad")
		c200tmp="atmrad"//sep//trim(a(iatm)%name)//statname(istat)//".gjf"
		open(10,file=c200tmp,status="replace")
		write(10,"(a)") "# "//trim(calclevel)//" out=wfn"
		write(10,*)
		write(10,"(a)") trim(a(iatm)%name)//statname(istat)
		write(10,*)
		write(10,"(2i3)") istat,chgmulti(iele,istat)
		write(10,"(a)") a(iatm)%name
		write(10,*)
		c200tmp="atmrad"//sep//trim(a(iatm)%name)//statname(istat)//".wfn"
		write(10,"(a)") trim(c200tmp)
		write(10,*)
		write(10,*)
		close(10)
		
		!Start calculation
		c80tmp="atmrad"//sep//trim(a(iatm)%name)//statname(istat)
		call runcommand('"'//trim(gaupath)//'" "'//trim(c80tmp)//'.gjf" "'//trim(c80tmp)//'"')
		
		!Check if Gaussian task was successfully finished
		if (isys==1) then
			inquire(file=trim(c80tmp)//".out",exist=alive)
		else
			inquire(file=trim(c80tmp)//".log",exist=alive)
		end if
		if (alive) then
			if (isys==1) then
				open(10,file=trim(c80tmp)//".out",status="old")
			else
				open(10,file=trim(c80tmp)//".log",status="old")
			end if
			call loclabel(10,"Normal termination",igaunormal)
			close(10)
			if (igaunormal==0) then
				write(*,"(a)") " Gaussian running may be failed! Please manually check Gaussian input and output files in atmrad folder"
				write(*,*) "Press ENTER button to continue"
				read(*,*)
			end if
		else
			write(*,"(a)") " Gaussian running may be failed! Please manually check Gaussian input and output files in atmrad folder"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
		end if
	end do
end do

!All element wfn files have been generated, now calculate corresponding radial density file (.rad)
!Existing .rad file will not be recalculated
write(*,*)
!write(*,*) "Generating atomic radial density from atomic wfn file..."
do iatm=1,ncenter
	iele=a_org(iatm)%index
	do istat=-3,3
		if (chgmulti(iele,istat)==0) cycle !Undefined state, does not need to calculate
		c80tmp="atmrad"//sep//trim(a_org(iatm)%name)//statname(istat)
		inquire(file=trim(c80tmp)//".rad",exist=alive)
		if (alive) cycle
		inquire(file=trim(c80tmp)//".wfn",exist=alive)
		if (.not.alive) then
			write(*,"(' Error: ',a,' was not found!')") trim(c80tmp)//".wfn"
			write(*,*) "If you want to skip, press ENTER button directly"
			read(*,*)
			cycle
		end if
		write(*,"(' Converting ',a,' to ',a)") trim(c80tmp)//".wfn",trim(c80tmp)//".rad"
		call atmwfn2atmrad(trim(c80tmp)//".wfn",trim(c80tmp)//".rad")
	end do
end do

!Recover to the firstly loaded file if atomic wavefunction is calculated above
if (icalcatmwfn==1) then
	call dealloall(1)
	call readinfile(firstfilename,1)
end if
end subroutine




!!------- Generate atomic radial density from atomic .wfn file
!The code is adapted from sphatmraddens
subroutine atmwfn2atmrad(infile,outfile)
use defvar
use functions
implicit real*8 (a-h,o-z)
character(len=*) infile,outfile
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:),radpos(:),sphavgval(:)
call dealloall(0)
call readinfile(infile,1)
truncrho=1D-8
rlow=0D0
rhigh=12
nsphpt=974
nradpt=200 !Totally 200 radial points, but the number of point is truncated at truncrho (because the interpolation routine doesn't work well for very low value)
allocate(potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt),radpos(nradpt),sphavgval(nradpt))
sphavgval=0
call Lebedevgen(nsphpt,potx,poty,potz,potw)
!$OMP PARALLEL DO SHARED(sphavgval,radpos) PRIVATE(irad,radx,radr,isph,rnowx,rnowy,rnowz) schedule(dynamic) NUM_THREADS(nthreads)
do irad=1,nradpt
	radx=cos(irad*pi/(nradpt+1))
	radr=(1+radx)/(1-radx) !Becke transform
	radpos(irad)=radr
	do isph=1,nsphpt
		rnowx=potx(isph)*radr
		rnowy=poty(isph)*radr
		rnowz=potz(isph)*radr
		sphavgval(irad)=sphavgval(irad)+fdens(rnowx,rnowy,rnowz)*potw(isph)
	end do
end do
!$OMP END PARALLEL DO
open(10,file=outfile,status="replace")
write(10,*) count(sphavgval>truncrho)
do irad=nradpt,1,-1
	if (sphavgval(irad)>truncrho) write(10,"(f20.12,E18.10)") radpos(irad),sphavgval(irad)
end do
close(10)
end subroutine





!!============================ MBIS ============================!!
!!============================ MBIS ============================!!
!!============================ MBIS ============================!!
!!============================ MBIS ============================!!
!!============================ MBIS ============================!!
!A wrapper of subroutine MBIS to automatically set efficient radpot and sphpot for imode=0
subroutine MBIS_wrapper(itype,imode)
use defvar
implicit real*8 (a-h,o-z)
integer itype,imode
nradpotold=radpot
nsphpotold=sphpot
if (iautointgrid==1) then
  	radpot=30
  	sphpot=302 !I carefully tested, difference between 170 and 434 in atomic charge is at most 0.003, 302 is safer
 	if (any(a%index>18)) radpot=40
 	if (any(a%index>36)) radpot=50
 	if (any(a%index>54)) radpot=60
end if
call MBIS(itype,imode)
if (iautointgrid==1) then
	radpot=nradpotold
	sphpot=nsphpotold
end if
end subroutine

!!--------- Calculate MBIS charge or atomic radial density. Suitable for both isolated and periodic systems
!itype: 1=Calculate and print charges =2: Only generate atomic spaces, namely filling "atmraddens" global array by final radial density of each atom
!imode:
!0=Atomic center grid, only for isolated systems
!1=Evenly distributed grid, calculate actual density from periodic wavefunction
!2=Evenly distributed grid, actual density is directly taken from grid data in memory
subroutine MBIS(itype,imode)
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
integer itype,imode
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 charge(ncenter),lastcharge(ncenter) !Atomic charges of current iter. and last iter.
real*8 beckeweigrid(radpot*sphpot),tvec(3)
integer,parameter :: maxshell=6
real*8 shpop(maxshell,ncenter),shsig(maxshell,ncenter) !Shell populations and shell sigma (width)
real*8 shpopnew(maxshell,ncenter),shsignew(maxshell,ncenter) !New shell populations and shell sigma during iteration
real*8 shpopnew_tmp(maxshell,ncenter),shsignew_tmp(maxshell,ncenter)
real*8 rho0sh(maxshell,ncenter) !Shell density at current grid
real*8 tmpdens(radpot*sphpot,ncenter) !tmpdens(ipt,iatm) corresponds to contribution of iatm to molecular density at grid ipt, and meantime multiplied by single-center integration weight at that point
real*8 atmdis2min(ncenter)
integer mshell(ncenter) !Actual number of shells of atoms
integer :: maxcyc=500,ioutmedchg=0,ioutshell=0,ignorefar=1
real*8 :: crit=0.0001D0,eps=1D-14,dencut=1D-10

if (any(a%index>86)) then
    write(*,*) "Error: MBIS for elements beyond Rn is not supported"
    write(*,*) "Press ENTER button to exit"
    read(*,*)
    return
end if
if (imode==2.and.any(cubmat<0)) then
	write(*,"(a)") " Warning: Negative electron density has been found on some grids! In this case MBIS result will be inaccurate. &
    Note that negative electron density at core region can be caused by PAW calculation of VASP, if it is the present case, please consider not to use PAW or use other method to calculate atomic charges"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
end if

do while(.true.)
    write(*,*)
    call menutitle("MBIS",15,2)
    if (ignorefar==1) write(*,*) "-4 Toggle if reducing cost by ignoring atoms far from grid, current: Yes"
    if (ignorefar==0) write(*,*) "-4 Toggle if reducing cost by ignoring atoms far from grid, current: No"
    if (imode==0) write(*,*) "-3 Enter frj implementation of MBIS code"
    if (ioutshell==1) write(*,*) "-2 Toggle if outputting population and width of shells, current: Yes"
    if (ioutshell==0) write(*,*) "-2 Toggle if outputting population and width of shells, current: No"
    if (ioutmedchg==1) write(*,*) "-1 Toggle if outputting atomic charges during iterations, current: Yes"
    if (ioutmedchg==0) write(*,*) "-1 Toggle if outputting atomic charges during iterations, current: No"
    write(*,*) "0 Return"
	write(*,*) "1 Start calculation!"
	write(*,"(a,i4)") " 2 Set the maximum number of iterations, current:",maxcyc
	write(*,"(a,f10.6)") " 3 Set convergence criterion of atomic charges, current:",crit
	read(*,*) isel
    if (isel==0) then
        return
	else if (isel==-4) then
        if (ignorefar==1) then
            ignorefar=0
        else
            ignorefar=1
        end if
    else if (isel==-3) then
        call mbis_frj
	else if (isel==-2) then
		if (ioutshell==1) then
			ioutshell=0
		else
			ioutshell=1
		end if
	else if (isel==-1) then
		if (ioutmedchg==1) then
			ioutmedchg=0
		else
			ioutmedchg=1
		end if
	else if (isel==0) then
		return
	else if (isel==1) then
		exit
	else if (isel==2) then
		write(*,*) "Input maximum number of iterations, e.g. 30"
		read(*,*) maxcyc
	else if (isel==3) then
		write(*,*) "Input convergence criterion of atomic charges, e.g. 0.001"
		read(*,*) crit
	end if
end do

!Prepare actual density of present system at integration points
if (imode==0) then !Atomic center grids, only for isolated systems
	call walltime(iwalltime1)
	ntotpot=radpot*sphpot
	call gen1cintgrid(gridatmorg,iradcut)
	write(*,"(' Radial grids:',i4,'  Angular grids:',i5,'  Total:',i7,'  After pruning:',i7)") radpot,sphpot,radpot*sphpot,radpot*sphpot-iradcut*sphpot
	write(*,"(a)") " Calculating atomic contribution to electron density of present system on grid points..."
    ifinish=0
    call showprog(ifinish,ncenter)
	!$OMP PARALLEL DO SHARED(tmpdens,ifinish) PRIVATE(iatm,gridatm,beckeweigrid,dtmp) schedule(dynamic) NUM_THREADS(nthreads)
	do iatm=1,ncenter
		gridatm%value=gridatmorg%value
		gridatm%x=gridatmorg%x+a(iatm)%x
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
		do ipt=1+iradcut*sphpot,ntotpot
			dtmp = fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
			tmpdens(ipt,iatm) = dtmp*gridatm(ipt)%value*beckeweigrid(ipt)
		end do
		!$OMP CRITICAL
		ifinish=ifinish+1
		call showprog(ifinish,ncenter)
		!$OMP END CRITICAL
	end do
	!$OMP END PARALLEL DO
else if (imode==1) then !Calculate density from periodic wavefunction
    call setgrid_for_PBC(0.2D0,1)
	call calc_dvol(dvol)
    if (allocated(cubmat)) deallocate(cubmat)
    allocate(cubmat(nx,ny,nz))
	call walltime(iwalltime1)
	!Because uniform grid cannot integrate well core density, so temporarily disable EDFs
    nEDFprims_org=nEDFprims
    nEDFprims=0
    call delvirorb(1) !Delete high-lying virtual orbitals for faster calculation
    write(*,*) "Calculating electron density grid data..."
    call savecubmat(1,0,1)
    call delvirorb_back(1) !Restore to previous wavefunction
    nEDFprims=nEDFprims_org
else !Directly using loaded electron density from cub/VASP grid data, and transforming grid data information to cell information
	if (all(a%charge==0)) then
		write(*,*) "Error: All nuclear charges are zero! If this file was exported by CP2K, it is a bug. You need to manually &
        &edit the file so that effective nuclear charges (column 2 since line 8) are correctly recorded, otherwise atomic charges cannot be calculated"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
	end if
    call grid2cellinfo
	call calc_dvol(dvol)
	call walltime(iwalltime1)
end if

!Set initial sigma and population of various shells
shpop(:,:)=0
icore=1 !If consider core shells. If =0, initial population of core shells will be 0, and core shells will not be utilized during iteration (population and sigma will be zero throughout iterations)
!For density only representing valence electrons, I found ignoring core shells do not improve convergence. After first several iterations, core population automatically decreases to nearly zero
do iatm=1,ncenter
    iele = a(iatm)%index
    if (iele==0) then !Ghost atom, initialize as shsig=1 and with a tiny population. This scheme is defined by frj
        mshell(iatm)=1
        shsig(1,iatm) = 1
        shpop(1,iatm)=1D-3
    else if (iele<=2) then
        mshell(iatm) = 1
        shsig(1,iatm) = 1D0/(2*iele)
        shpop(1,iatm)=iele
    else if (iele<=10) then
        mshell(iatm) = 2
        shsig(1,iatm) = 1D0/(2*iele)
        shsig(2,iatm) = 1D0/2
        if (icore==1) shpop(1,iatm)=2
        shpop(2,iatm)=iele-2
    else if (iele<=18) then
        mshell(iatm) = 3
        shsig(1,iatm) = 1D0/(2*iele)
        shsig(2,iatm) = 1D0/(2*sqrt(dfloat(iele)))
        shsig(3,iatm) = 1D0/2
        if (icore==1) shpop(1,iatm)=2
        if (icore==1) shpop(2,iatm)=8
        shpop(3,iatm)=iele-10
    else if (iele<=36) then
        mshell(iatm) = 4
        shsig(1,iatm) = 1D0/(2*iele)
        do ishell=2,3
            shsig(ishell,iatm) = 1D0/(2*iele**(1-dfloat(ishell-1)/(mshell(iatm)-1)))
        end do
        shsig(4,iatm) = 1D0/2
        if (icore==1) shpop(1,iatm)=2
        if (icore==1) shpop(2,iatm)=8
        if (icore==1) shpop(3,iatm)=8
        shpop(4,iatm)=iele-18
    else if (iele<=54) then
        mshell(iatm) = 5
        shsig(1,iatm) = 1D0/(2*iele)
        do ishell=2,4
            shsig(ishell,iatm) = 1D0/(2*iele**(1-dfloat(ishell-1)/(mshell(iatm)-1)))
        end do
        shsig(5,iatm) = 1D0/2
        if (icore==1) shpop(1,iatm)=2
        if (icore==1) shpop(2,iatm)=8
        if (icore==1) shpop(3,iatm)=8
        if (icore==1) shpop(4,iatm)=18
        shpop(5,iatm)=iele-36
    else if (iele<=86) then
        mshell(iatm) = 6
        shsig(1,iatm) = 1D0/(2*iele)
        do ishell=2,5
            shsig(ishell,iatm) = 1D0/(2*iele**(1-dfloat(ishell-1)/(mshell(iatm)-1)))
        end do
        shsig(6,iatm) = 1D0/2
        if (icore==1) shpop(1,iatm)=2
        if (icore==1) shpop(2,iatm)=8
        if (icore==1) shpop(3,iatm)=8
        if (icore==1) shpop(4,iatm)=18
        if (icore==1) shpop(5,iatm)=18
        shpop(6,iatm)=iele-54
    end if
end do

write(*,*)
write(*,*) "Performing MBIS iterations to refine atomic spaces..."
lastcharge=0

do icyc=1,maxcyc
	if (ioutmedchg==1) write(*,*)
	if (icyc==1) then
		write(*,"(' Cycle',i5)") icyc
	else
		write(*,"(' Cycle',i5,'   Maximum change:',f12.8)") icyc,varmax
	end if
    
    !Monitor population and width of shells
	!write(*,*) "Population of each shell"
	!do iatm=1,ncenter
	!	write(*,"(i5,'(',a,'):',6f10.6,' q(atm):',f11.6)") iatm,a(iatm)%name,(shpop(ish,iatm),ish=1,6),a(iatm)%charge-sum(shpop(1:mshell(iatm),iatm))
	!end do
	!write(*,*) "Width (sigma) of each shell in Bohr"
	!do iatm=1,ncenter
	!	write(*,"(i5,'(',a,'):',6f12.6)") iatm,a(iatm)%name,(shsig(ish,iatm),ish=1,mshell(iatm))
	!end do
	
	shpopnew(:,:)=0 !New population of shells of various atoms
	shsignew(:,:)=0 !New sigma of shells of various atoms
    if (imode==0) then !Using multicenter integration to evaluate population of various shells of various atoms based on present sigma (Eq. 18 of MBIS paper)
		do iatm=1,ncenter
			gridatm%x=gridatmorg%x+a(iatm)%x
			gridatm%y=gridatmorg%y+a(iatm)%y
			gridatm%z=gridatmorg%z+a(iatm)%z
			do ipt=1+iradcut*sphpot,ntotpot
				rho0sh(:,:)=0 !Record {rho_0_Ai} at present grid, namely density of shell i of atom A at this integration point, constructed by Eq. 7 of MBIS paper
				rho0=0 !rho_0, namely reference density at this integration point, constructed by Eq. 6 of MBIS paper
				do jatm=1,ncenter
					dx = gridatm(ipt)%x - a(jatm)%x
					dy = gridatm(ipt)%y - a(jatm)%y
					dz = gridatm(ipt)%z - a(jatm)%z
					dis2 = dx*dx + dy*dy + dz*dz
                    if (ignorefar==1.and.dis2>atmrhocutsqr(a(jatm)%index)) cycle !My tested showed that this reduce cost nearly half, while accuracy lost is negligible
                    dis=dsqrt(dis2)
					do ishell=1,mshell(jatm)
						sigval = shsig(ishell,jatm)
						tmp = shpop(ishell,jatm)/sigval**3/8/pi*exp(-dis/sigval) !Eq. 7 of MBIS paper
						if (tmp<dencut) tmp = 0 !I don't know why frj introduced this criterion. Seems that this can make insignificant grid ignored and reduce cost (because of wtot>0)?
						rho0sh(ishell,jatm) = tmp
						rho0 = rho0 + tmp
					end do
				end do
				!Accumulate contribution of this integration grid to new population and sigma of shells
				tmpden = tmpdens(ipt,iatm)
				if (rho0>0.and.tmpden>eps) then
					do jatm=1,ncenter
						dx = gridatm(ipt)%x - a(jatm)%x
						dy = gridatm(ipt)%y - a(jatm)%y
						dz = gridatm(ipt)%z - a(jatm)%z
						dis2 = dx*dx + dy*dy + dz*dz
						if (ignorefar==1.and.dis2>atmrhocutsqr(a(jatm)%index)) cycle !My tested showed that this reduce cost nearly half, while accuracy lost is negligible (<0.0004)
                        dis=dsqrt(dis2)
						do ishell=1,mshell(jatm)
							shpopnew(ishell,jatm) = shpopnew(ishell,jatm) + tmpden*rho0sh(ishell,jatm)/rho0 !Eq. 18 of MBIS paper
							shsignew(ishell,jatm) = shsignew(ishell,jatm) + tmpden*rho0sh(ishell,jatm)/rho0*dis !Integral part of Eq. 19 of MBIS paper
						end do
					end do
				end if
			end do
		end do
    
    else !Using evenly distributed grids
		ifinish=0;ishowprog=1
		ntmp=floor(ny*nz/100D0)
		!$OMP PARALLEL SHARED(shpopnew,shsignew,ifinish,ishowprog) PRIVATE(shpopnew_tmp,shsignew_tmp,i,j,k,rho0sh,rho0,tmpx,tmpy,tmpz,tvec, &
		!$OMP ic,jc,kc,icell,jcell,kcell,iatm,dx,dy,dz,dis,dis2,dis2min,ishell,sigval,tmp,tmp2,tmp3,tmpden,atmdis2min) NUM_THREADS(nthreads)
		shpopnew_tmp(:,:)=0
        shsignew_tmp(:,:)=0
		!$OMP DO schedule(dynamic) collapse(2)
		do k=1,nz
			do j=1,ny
				do i=1,nx
					if (cubmat(i,j,k)<1D-11) cycle !VASP PAW density is intrinsically incompatible with MBIS because of negative density around core. However, ignoring negative density can make calculation still feasible though inaccurate
					rho0sh(:,:)=0 !Record {rho_0_Ai} at present grid, namely density of shell i of atom A at this integration point, constructed by Eq. 7 of MBIS paper
					rho0=0 !rho_0, namely reference density at this integration point, constructed by Eq. 6 of MBIS paper
					call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					!call getpointcell(tmpx,tmpy,tmpz,ic,jc,kc)
                    atmdis2min(:)=1D10
					do icell=-PBCnx,+PBCnx
						do jcell=-PBCny,+PBCny
							do kcell=-PBCnz,+PBCnz
								call tvec_PBC(icell,jcell,kcell,tvec)
								do iatm=1,ncenter
									dx=a(iatm)%x+tvec(1)-tmpx
									dy=a(iatm)%y+tvec(2)-tmpy
									dz=a(iatm)%z+tvec(3)-tmpz
									dis2=dx*dx+dy*dy+dz*dz
                                    if (dis2<atmdis2min(iatm)) atmdis2min(iatm)=dis2
									if (dis2>atmrhocutsqr(a(iatm)%index)) cycle !Ignore atoms that do not contribute notably to present grid
                                    dis=dsqrt(dis2)
									do ishell=1,mshell(iatm)
										sigval = shsig(ishell,iatm)
                                        if (sigval==0) cycle
										tmp = shpop(ishell,iatm)/sigval**3/8/pi*exp(-dis/sigval) !Eq. 7 of MBIS paper
										rho0sh(ishell,iatm) = rho0sh(ishell,iatm) + tmp
										rho0 = rho0 + tmp
									end do
								end do
							end do
						end do
					end do
                    
					!Accumulate contribution of this integration grid to new population and sigma of shells
					tmpden = cubmat(i,j,k)*dvol
					if (rho0>0.and.tmpden>eps) then
						do iatm=1,ncenter
							tmp2=tmpden/rho0
                            tmp3=tmp2*dsqrt(atmdis2min(iatm))
							do ishell=1,mshell(iatm)
								shpopnew_tmp(ishell,iatm) = shpopnew_tmp(ishell,iatm) + tmp2*rho0sh(ishell,iatm) !Eq. 18 of MBIS paper
								shsignew_tmp(ishell,iatm) = shsignew_tmp(ishell,iatm) + tmp3*rho0sh(ishell,iatm) !Integral part of Eq. 19 of MBIS paper
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
		shpopnew(:,:)=shpopnew(:,:)+shpopnew_tmp(:,:)
		shsignew(:,:)=shsignew(:,:)+shsignew_tmp(:,:)
		!$OMP END CRITICAL
		!$OMP END PARALLEL
		if (ishowprog/=0) call showprog(100,100)
    end if
    
	!write(*,*) "Population of each shell"
	!do iatm=1,ncenter
	!	write(*,"(i5,'(',a,'):',6f10.6,' q(atm):',f11.6)") iatm,a(iatm)%name,(shpopnew(ish,iatm),ish=1,6),a(iatm)%charge-sum(shpopnew(1:mshell(iatm),iatm))
	!end do
	!write(*,*) "Width (sigma) of each shell in Bohr"
	!do iatm=1,ncenter
	!	write(*,"(i5,'(',a,'):',6f12.6)") iatm,a(iatm)%name,(shsignew(ish,iatm),ish=1,mshell(iatm))
	!end do
 !   write(*,*) "--------------------------"
    
    !Include prefix part of Eq. 19 of MBIS paper
	do iatm=1,ncenter
		do ish=1,mshell(iatm)
			if (shpopnew(ish,iatm)>0) shsignew(ish,iatm)=shsignew(ish,iatm)/(3*shpopnew(ish,iatm))
        end do
    end do
    
    !Summing up shell populations to atomic population and get atomic charge
    do iatm=1,ncenter
		tmppop=sum(shpopnew(1:mshell(iatm),iatm)) !Atomic population
        if (nEDFelec==0.or.imode>0) then !Note that EDFs were not involved in evaluating system density when using even grids (imode>0)
            charge(iatm) = a(iatm)%charge - tmppop
        else !EDF is used for some atoms. Core electron density represented by EDF has been integrated, so nuclear charge should be augmented by nEDFelecatm
            charge(iatm) = a(iatm)%charge+nEDFelecatm(iatm) - tmppop
        end if
        if (ioutmedchg==1) write(*,"(i5,'(',a,')   charge:',f12.6)") iatm,a(iatm)%name,charge(iatm)
    end do

    !Check convergence
    varmax=maxval(abs(charge(:)-lastcharge(:)))
	if (varmax<crit.or.icyc==maxcyc) then
        if (varmax<crit) write(*,"(/,a,f10.6)") " All atomic charges have converged to criterion of",crit
        if (icyc==maxcyc) write(*,"(/,' Convergence failed within',i4,' cycles!')") maxcyc
		exit
	end if
    
    !Update atomic charges, shell population and sigma
	lastcharge(:)=charge(:)
    shpop(:,:)=shpopnew(:,:)
    shsig(:,:)=shsignew(:,:)
end do

if (itype==1) then
	write(*,"(' Sum of all raw charges:',f14.8)") sum(charge(:))
	!Normalize atomic charges. This is not feasible if only grid data is available, &
	!because in this case the nelec used in "normalize_atmchg" is simply guessed by assuming system is neutral
	if (imode==1) call normalize_atmchg(charge(:))
	!Print final atomic charges
	call printatmchg(charge(:))

	if (allocated(frag1)) then
		write(*,"(/,' Fragment charge:',f14.8)") sum(charge(frag1))
		write(*,"(' Fragment population:',f14.8)") sum(a(frag1)%charge) - sum(charge(frag1))
	end if

	if (ioutshell==1) then
		write(*,*)
		write(*,*) "Population of each shell"
		do iatm=1,ncenter
		   write(*,"(i5,'(',a,'):',6f12.6)") iatm,a(iatm)%name,(shpop(ish,iatm),ish=1,mshell(iatm))
		end do
		write(*,*)
		write(*,*) "Width (sigma) of each shell in Bohr"
		do iatm=1,ncenter
		   write(*,"(i5,'(',a,'):',6f12.6)") iatm,a(iatm)%name,(shsig(ish,iatm),ish=1,mshell(iatm))
		end do
	end if
end if

call walltime(iwalltime2)
write(*,"(/,' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

if (itype==1) then !Output charges
    call outatmchg(10,charge(:))
else if (itype==2) then !Generate radial density of every atom
    if (allocated(atmradnpt)) deallocate(atmradnpt)
    if (allocated(atmraddens)) deallocate(atmraddens)
    allocate(atmradnpt(ncenter),atmraddens(200,ncenter))
    do iatm=1,ncenter
        do ipt=1,200
            tmprho=0
            do ishell=1,mshell(iatm)
				sigval=shsig(ishell,iatm)
                tmprho = tmprho + shpop(ishell,iatm)/sigval**3/8/pi*exp(-atmradpos(ipt)/sigval)
            end do
            atmraddens(ipt,iatm)=tmprho
            if (tmprho<1D-8) then !Electron density truncation
                atmradnpt(iatm)=ipt
                exit
            end if
        end do
    end do
    write(*,*) "Construction of MBIS atomic spaces has been finished!"
end if

end subroutine








!!----------------------------------------
!!--------- Calculate EEM charge ---------
!!----------------------------------------
subroutine EEM
use defvar
use util
implicit real*8 (a-h,o-z)
integer,parameter :: maxBO=3 !Maximum of possible bond order
character c200tmp*200
real*8 EEMmat(ncenter+1,ncenter+1),EEMarr(ncenter+1),qarr(ncenter+1)
real*8 kappa,Aparm(nelesupp,maxBO),Bparm(nelesupp,maxBO) !If parameter is -1, means undefined parameter
real*8 :: chgnet=0

if (ifiletype/=11.and.ifiletype/=13) then
	write(*,"(/,a)") " Warning: Commonly MDL Molfile (.mol) or .mol2 file should be used as input file, &
	&because it contains bond information, which is needed by present function. If you want to let Multiwfn guess interatomic connectivity &
	&and then calculate EEM charges, you can input ""g"", however bond multiplicity cannot be determined in this way, and thus &
	&the EEM charges may be problematic if the employed EEM parameters explicitly involve bond multiplicity"
	write(*,*) "If you simply want to return, press ENTER button"
	read(*,"(a)") c200tmp
	if (index(c200tmp,'g')/=0) then
		call genconnmat(1,0)
	else
		return
	end if
end if

iparmset=2
call genEEMparm(iparmset,kappa,Aparm,Bparm)
isel2=-10
	
EEMcyc: do while(.true.)
	write(*,*)
	write(*,*) "           ------ Electronegativity Equalization Method (EEM) ------"
	write(*,*) "-1 Return"
	write(*,*) "0 Start calculation"
	write(*,*) "1 Choose EEM parameters"
	write(*,"(a,f4.1)") " 2 Set net charge, current:",chgnet
	read(*,*) isel
	if (isel==-1) then
		return
	else if (isel==1) then
		do while(.true.)
			if (isel2/=-1) then
				write(*,*)
				write(*,*) "Present EEM parameters:"
				write(*,"(' kappa',f12.6)") kappa
				do iele=1,nelesupp
					do imulti=1,maxBO
						if (Aparm(iele,imulti)/=-1) write(*,"(1x,a,'  Multiplicity:',i2,'    A:',f9.5, '    B:',f9.5)") ind2name(iele),imulti,Aparm(iele,imulti),Bparm(iele,imulti)
					end do
				end do
			end if
			write(*,*)
			write(*,*) "-2 Return"
			write(*,*) "-1 Export present parameters to external file"
			write(*,*) "0 Load parameters from external file"
			write(*,*) "1 Use parameters fitted to HF/STO-3G Mulliken charge, IJMS, 8, 572 (2007)"
			write(*,*) "2 Use parameters fitted to B3LYP/6-31G* CHELPG charge, JCC, 30, 1174 (2009)"
			write(*,*) "3 Use parameters fitted to HF/6-31G* CHELPG charge, JCC, 30, 1174 (2009)"
			write(*,*) "4 Use parameters fitted to B3LYP/6-311G* NPA charge, J Cheminform, 8, 57(2016)"
			read(*,*) isel2
			if (isel2==-2) then
				exit
			else if (isel2==-1) then
				open(10,file="EEMparm.txt",status="replace")
				write(10,"(f12.6)") kappa
				do iele=1,nelesupp
					do imulti=1,maxBO
						if (Aparm(iele,imulti)/=-1) write(10,"(1x,a,i3,2f9.5)") ind2name(iele),imulti,Aparm(iele,imulti),Bparm(iele,imulti)
					end do
				end do
				close(10)
				write(*,*) "Parameters have been exported to EEMparm.txt in current folder"
			else if (isel2==0) then
				write(*,*) "Input path of parameter file, e.g. C:\aqours.txt"
				do while(.true.)
					read(*,"(a)") c200tmp
					inquire(file=c200tmp,exist=alive)
					if (alive) exit
					write(*,*) "Cannot find the file, input again"
				end do
				Aparm=-1
				Bparm=-1
				open(10,file=c200tmp,status="old")
				read(10,*) kappa
				nload=0
				do while(.true.)
					read(10,*,iostat=ierror) c200tmp,imulti,tmpA,tmpB
					if (ierror/=0) exit
                    call elename2idx(c200tmp,iele)
					Aparm(iele,imulti)=tmpA
					Bparm(iele,imulti)=tmpB
					nload=nload+1
				end do
				write(*,"(' Loaded',i5,' entries')") nload
				close(10)
			else
				call genEEMparm(isel2,kappa,Aparm,Bparm)
			end if
		end do
	else if (isel==2) then
		write(*,*) "Input net charge of the system, e.g. -1"
		read(*,*) chgnet
	else if (isel==0) then
		write(*,*) "Calculating..."
		write(*,*)
		!Construct EEM array
		EEMarr(ncenter+1)=chgnet
		do iatm=1,ncenter
			imulti=maxval(connmat(iatm,:))
			if (imulti>maxBO) then
				write(*,"(' Error: Multiplicity of atom',i5,' (',i2,') exceeded upper limit (',i2,')!')") iatm,imulti,maxBO
                write(*,"(a)") " The present EEM parameters do not support such bonding status, or connectivity in your input file is wrong"
				cycle EEMcyc
            else if (imulti==0) then
				write(*,"(' Error: Atom',i5,'(',a,') is not bonded to any atom!')") iatm,a(iatm)%name
                write(*,"(a)") " The present EEM parameters do not support such bonding status, or connectivity in your input file is wrong"
				cycle EEMcyc
			end if
			tmpval=Aparm(a(iatm)%index,imulti)
			if (tmpval==-1) then
				write(*,"(' Error: Parameter for atom',i5,'(',a,') is unavailable!')") iatm,a(iatm)%name
				cycle EEMcyc
			else
				EEMarr(iatm)=-tmpval
			end if
		end do

		!Construct EEM matrix
		EEMmat=0
		EEMmat(ncenter+1,1:ncenter)=1
		EEMmat(1:ncenter,ncenter+1)=-1
		do i=1,ncenter
			imulti=maxval(connmat(i,:))
			do j=1,ncenter
				if (i==j) then
					EEMmat(i,j)=Bparm(a(i)%index,imulti)
				else
					EEMmat(i,j)=kappa/(atomdist(i,j,1)*b2a)
				end if
			end do
		end do
		
		!Solve EEM equation
		qarr=matmul(invmat(EEMmat,ncenter+1),EEMarr)
		do iatm=1,ncenter
			write(*,"(' EEM charge of atom',i8,'(',a,'):',f15.10)") iatm,a(iatm)%name,qarr(iatm)
		end do
		write(*,"(' Electronegativity:',f12.6)") qarr(ncenter+1)
        call outatmchg(10,qarr(:))
	end if

end do EEMcyc
end subroutine

!---- Generate EEM parameters
subroutine genEEMparm(iset,kappa,Aparm,Bparm)
use defvar
integer,parameter :: maxBO=3 !Maximum of bond order
real*8 kappa,Aparm(nelesupp,maxBO),Bparm(nelesupp,maxBO)
Aparm=-1
Bparm=-1
if (iset==1) then !Parameters fitted to Mulliken charge at HF/STO-3G, Int. J. Mol. Sci., 8, 572-582 (2007)
	write(*,"(a)") " Parameters have been set to those fitted to Mulliken charges at HF/STO-3G, see Int. J. Mol. Sci., 8, 572-582 (2007)"
	kappa=0.44D0
	Aparm(1,1)= 2.396D0  !H
	Bparm(1,1)= 0.959D0
	Aparm(6,1)= 2.459D0  !C,multi=1
	Bparm(6,1)= 0.611D0
	Aparm(7,1)= 2.597D0  !N,multi=1
	Bparm(7,1)= 0.790D0
	Aparm(8,1)= 2.625D0  !O,multi=1
	Bparm(8,1)= 0.858D0
	Aparm(16,1)= 2.407D0  !S,multi=1
	Bparm(16,1)= 0.491D0
	Aparm(6,2)= 2.464D0  !C,multi=2
	Bparm(6,2)= 0.565D0
	Aparm(7,2)= 2.554D0  !N,multi=2
	Bparm(7,2)= 0.611D0
	Aparm(8,2)= 2.580D0  !O,multi=2
	Bparm(8,2)= 0.691D0
else if (iset==2) then !Parameters fitted to CHELPG charges at B3LYP/6-31G*, J. Comput. Chem., 30, 1174 (2009)
	write(*,"(a)") " Parameters have been set to those fitted to CHELPG charges at B3LYP/6-31G*, see J. Comput. Chem., 30, 1174 (2009)"
	kappa=0.302D0
	Aparm(35,1)= 2.659D0  !Br,multi=1
	Bparm(35,1)= 1.802D0
	Aparm(6,1)= 2.482D0  !C,multi=1
	Bparm(6,1)= 0.464D0
	Aparm(17,1)= 2.519D0  !Cl,multi=1
	Bparm(17,1)= 1.450D0
	Aparm(9,1)= 3.577D0  !F,multi=1
	Bparm(9,1)= 3.419D0
	Aparm(1,1)= 2.385D0  !H,multi=1
	Bparm(1,1)= 0.737D0
	Aparm(7,1)= 2.595D0  !N,multi=1
	Bparm(7,1)= 0.468D0
	Aparm(8,1)= 2.825D0  !O,multi=1
	Bparm(8,1)= 0.844D0
	Aparm(16,1)= 2.452D0  !S,multi=1
	Bparm(16,1)= 0.362D0
	Aparm(30,1)= 2.298D0  !Zn,multi=1
	Bparm(30,1)= 0.420D0
	Aparm(6,2)= 2.464D0  !C,multi=2
	Bparm(6,2)= 0.392D0
	Aparm(7,2)= 2.556D0  !N,multi=2
	Bparm(7,2)= 0.377D0
	Aparm(8,2)= 2.789D0  !O,multi=2
	Bparm(8,2)= 0.834D0
else if (iset==3) then !Parameters fitted to CHELPG charges at HF/6-31G*, J. Comput. Chem., 30, 1174 (2009)
	write(*,"(a)") " Parameters have been set to those fitted to CHELPG charges at HF/6-31G*, see J. Comput. Chem., 30, 1174 (2009)"
	kappa=0.227D0
	Aparm(35,1)= 2.615D0  !Br,multi=1
	Bparm(35,1)= 1.436D0
	Aparm(6,1)= 2.481D0  !C,multi=1
	Bparm(6,1)= 0.373D0
	Aparm(17,1)= 2.517D0  !Cl,multi=1
	Bparm(17,1)= 1.043D0
	Aparm(9,1)= 3.991D0  !F,multi=1
	Bparm(9,1)= 3.594D0
	Aparm(1,1)= 2.357D0  !H,multi=1
	Bparm(1,1)= 0.688D0
	Aparm(7,1)= 2.585D0  !N,multi=1
	Bparm(7,1)= 0.329D0
	Aparm(8,1)= 2.870D0  !O,multi=1
	Bparm(8,1)= 0.717D0
	Aparm(16,1)= 2.450D0  !S,multi=1
	Bparm(16,1)= 0.269D0
	Aparm(30,1)= 2.185D0  !Zn,multi=1
	Bparm(30,1)= 0.375D0
	Aparm(6,2)= 2.475D0  !C,multi=2
	Bparm(6,2)= 0.292D0
	Aparm(7,2)= 2.556D0  !N,multi=2
	Bparm(7,2)= 0.288D0
	Aparm(8,2)= 2.757D0  !O,multi=2
	Bparm(8,2)= 0.621D0
else if (iset==4) then !Parameters fitted to NPA charges at B3LYP/6-311G*, see J. Cheminform., 8, 57 (2016)
!The data were taken from "13321_2016_171_MOESM5_ESM.rar Additional file 5" of supplmental material
!13321_2016_171_MOESM5_ESM\neemp\ideal_q1\set3_DE_RMSD_B3LYP_6311Gd_NPA_cross_ideal\output_set3_DE_RMSD_B3LYP_6311Gd_NPA_cross_ideal_all
	write(*,"(a)") " Parameters have been set to those fitted to NPA charges at B3LYP/6-311G*, they were extracted from SI of J. Cheminform., 8, 57 (2016)"
	kappa=0.4024D0
	Aparm(1,1)= 2.4598D0  !H,multi=1
	Bparm(1,1)= 0.9120D0
	Aparm(6,1)= 2.5957D0  !C,multi=1
	Bparm(6,1)= 0.5083D0
	Aparm(6,2)= 2.6029D0  !C,multi=2
	Bparm(6,2)= 0.5021D0
	Aparm(6,3)= 2.5326D0  !C,multi=3
	Bparm(6,3)= 0.5932D0
	Aparm(7,1)= 2.7802D0  !N,multi=1
	Bparm(7,1)= 0.7060D0
	Aparm(7,2)= 2.7141D0  !N,multi=2
	Bparm(7,2)= 0.5366D0
	Aparm(7,3)= 2.6391D0  !N,multi=3
	Bparm(7,3)= 0.5171D0
	Aparm(8,1)= 2.9496D0  !O,multi=1
	Bparm(8,1)= 0.8264D0
	Aparm(8,2)= 2.8595D0  !O,multi=2
	Bparm(8,2)= 0.6589D0
	Aparm(9,1)= 2.9165D0  !F,multi=1
	Bparm(9,1)= 0.8427D0
	Aparm(15,2)= 2.1712D0  !P,multi=2
	Bparm(15,2)= 0.4802D0
	Aparm(16,1)= 2.5234D0  !S,multi=1
	Bparm(16,1)= 0.3726D0
	Aparm(16,2)= 2.5334D0  !S,multi=2
	Bparm(16,2)= 0.3519D0
	Aparm(17,1)= 2.5625D0  !Cl,multi=1
	Bparm(17,1)= 0.9863D0
	Aparm(35,1)= 2.4772D0  !Br,multi=1
	Bparm(35,1)= 1.2131D0
end if
end subroutine





!!----- A routine directly return Hirshfeld charges based on built-in density, adapted from "spacecharge" routine
!Mainly used by obtaining the charges used in conceptual density function theory module
subroutine genHirshfeld(charge)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8 charge(ncenter)
real*8 molrho(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),selfdens(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

!Generate quadrature point and weighs by combination of Gauss-Chebyshev and Lebedev grids
call gen1cintgrid(gridatmorg,iradcut)

do iatm=1,ncenter
	gridatm%value=gridatmorg%value !Weight in this grid point
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	!Calculate molecular density first
	!$OMP parallel do shared(molrho) private(i) num_threads(nthreads)
	do i=1,radpot*sphpot
		molrho(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
	end do
	!$OMP end parallel do
	!Calc free atomic density to obtain promolecule density
	promol=0D0
	do jatm=1,ncenter
		!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
		do ipt=1,radpot*sphpot
			tmpdens(ipt)=calcatmdens(jatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,0)
		end do
		!$OMP end parallel do
		promol=promol+tmpdens
		if (jatm==iatm) selfdens=tmpdens
	end do
	!Now we have needed data in hand, calculate atomic charges
	tmpcharge=0D0
	do i=1,radpot*sphpot
		if (promol(i)/=0D0) then
			tmpv=selfdens(i)/promol(i)*molrho(i)*gridatm(i)%value
			tmpcharge=tmpcharge-tmpv
		end if
	end do
	if (nEDFelec==0) then
		charge(iatm)=a(iatm)%charge+tmpcharge
	else !EDF is used for some atoms. Core electron density represented by EDF has been integrated, so nuclear charge should be augmented by nEDFelecatm
		charge(iatm)=a(iatm)%charge+nEDFelecatm(iatm)+tmpcharge
    end if
    call showprog(iatm,ncenter)
end do
end subroutine




!!--------------- Calculate Gasteiger charges
!Parameter comes from its original paper: Tetrahedron, 36, 3219 (1980)
!
!Can compare: obabel 1.mol -osmi --partialcharge gasteiger --print (result is always identical to Avogadro)
!But for conjugated systems, the result of present code is somewhat different to those.
!
!Can also compare: antechamber -i 1.mol2 -fi mol2 -o test.prepin -fo prepi -c gas
!The charges produced by present code are identical and directly comparable to those in the intermediate file ANTECHAMBER_PREP.AC
subroutine gasteiger
use defvar
use util
implicit real*8 (a-h,o-z)
character selectyn
integer,parameter :: maxbond=4
integer eleidx
!The a,b,c and initial charge. Row is element index, column corresponds to number of bonds (hybridzation state)
real*8 parma(nelesupp,maxbond),parmb(nelesupp,maxbond),parmc(nelesupp,maxbond),initchg(nelesupp,maxbond)
real*8 atom_a(ncenter),atom_b(ncenter),atom_c(ncenter),delta_q(ncenter),eleneg(ncenter),elenegQ1(ncenter),charge(ncenter)
integer nbond(ncenter) !Number of bonds of each atom
integer bondlist(ncenter,maxbond) !(i,j) element is the index of the jth atom that connected to atom i

!Parameters are extracted from original paper, while a few are supplemented by consulting GASPARM.dat in AmberTools
!Despite nearly ionic atom type is available in GASPARM.dat, they are not taken into account, because they only occur in very rare cases such H2SO4, H3PO4...
parma=0;parmb=0;parmc=0;initchg=0
!H
parma(1,1)=7.17D0; parmb(1,1)=6.24D0; parmc(1,1)=-0.56D0 !One bond
!C
parma(6,1:2)=10.39D0; parmb(6,1:2)=9.45D0; parmc(6,1:2)=0.73D0 !1 or 2 bonds
parma(6,3)=8.79D0;  parmb(6,3)=9.32D0; parmc(6,3)=1.51D0 !3 bonds
parma(6,4)=7.98D0;  parmb(6,4)=9.18D0; parmc(6,4)=1.88D0 !4 bonds
!N
parma(7,1)=15.68D0; parmb(7,1)=11.70D0; parmc(7,1)=-0.27D0 !One bond
parma(7,2)=12.87D0; parmb(7,2)=11.15D0; parmc(7,2)=0.85D0 !Two bonds
parma(7,3)=11.54D0; parmb(7,3)=10.82D0; parmc(7,3)=1.36D0 !Three bonds. The parameter of N in sp2 state will be assigned later
parma(7,4)= 0.00D0; parmb(7,4)=11.86D0; parmc(7,4)=11.86D0 !Four bonds
!P
parma(15,:)=8.90D0; parmb(15,:)=8.24D0; parmc(15,:)=0.96D0
!O
parma(8,1)=17.07D0; parmb(8,1)=13.79D0; parmc(8,1)=0.47D0 !One bond
parma(8,2)=14.18D0; parmb(8,2)=12.92D0; parmc(8,2)=1.39D0 !Two bonds
!S
parma(16,1)=10.88D0; parmb(16,1)=9.485D0; parmc(16,1)=1.325D0 !One bond
parma(16,2:3)=10.14D0; parmb(16,2:3)=9.13D0; parmc(16,2:3)=1.38D0 !Two or three bonds
parma(16,4)=12.00D0; parmb(16,4)=10.805D0; parmc(16,4)=1.195D0 !Four bonds
!F
parma(9,1)=14.66D0; parmb(9,1)=13.85D0; parmc(9,1)=2.31D0 !One bond
!Cl
parma(17,1)=11.00D0; parmb(17,1)=9.69D0; parmc(17,1)=1.35D0 !One bond
!Br
parma(35,1)=10.08D0; parmb(35,1)=8.47D0; parmc(35,1)=1.16D0 !One bond
!I
parma(53,1)=9.90D0; parmb(53,1)=7.96D0; parmc(53,1)=0.96D0 !One bond

!Below are GASPARM.dat in AmberTools. The ones not explicitly taken into account in Multiwfn is labelled with *
!Antechamber first determines Amber atom types, then convert to Gasgeiter types according to ATOMTYPE_GAS.DEF,&
!then checks GASPARM.DAT to determine actual parameters, finally invokes "charge.c" to carry out Gasteiger calculation.&
!The resulting Gasteiger type can be found from intermediate file ANTECHAMBER_GAS_AT.AC
!The O1 in below information corresponds to oxygen in =O case
!                  a       b       c      d     formal_charge
! GASPARM	h	  7.17	  6.24	 -0.56	 20.02	  0.00  !numbond=1
! GASPARM	c1	 10.39	  9.45	  0.73	 20.57	  0.00  !numbond=1 or 2
! GASPARM	c2	  8.79	  9.32	  1.51	 19.62	  0.00  !numbond=3
! GASPARM	c3	  7.98	  9.18	  1.88	 19.04	  0.00  !numbond=4
!*GASPARM	cg	  8.79	  9.32	  1.51	 19.62	  0.04  !numbond=3 (N3,N3,N3)
! GASPARM	n1	 15.68	 11.70	 -0.27	 27.11	  0.00  !numbond=1
! GASPARM	n2	 12.87	 11.15	  0.85	 24.87	  0.00  !numbond=2
! GASPARM	n3	 11.54	 10.82	  1.36	 23.72	  0.00  !numbond=3
! GASPARM	na 	 12.32	 11.20	  1.34	 24.86	  0.00  !numbond=3 Sp2 N with three connected atoms, e.g. piptide bond
!*GASPARM	na+	 12.32	 11.20	  1.34	 24.86	  1.00  !numbond=3 charged na
!*GASPARM	ng	 12.32	 11.20	  1.34	 24.86	  0.32  !numbond=3 (C3(N3,N3))
! GASPARM	n4	  0.00	 11.86	 11.86	 23.72	  1.00  !numbond=4
! GASPARM	o2	 17.07	 13.79	  0.47	 31.33	  0.00  !numbond=1
! GASPARM	o3	 14.18	 12.92	  1.39	 28.49	  0.00  !numbond=2
!*GASPARM	o-1	 17.07	 13.79	  0.47	 31.33	 -1.00  !numbond=1 (C4)
!*GASPARM	o-2	 17.07	 13.79	  0.47	 31.33	 -0.50  !numbond=1 (C3(O1))
! GASPARM	os	 17.07	 13.79	  0.47	 31.33	 -1.00  !numbond=1 (S4) or (S3)
!*GASPARM	op#	 17.07	 13.79	  0.47	 31.33	 -1.00  !numbond=1 (P4(O1,O1,O1))
!*GASPARM	op	 17.07	 13.79	  0.47	 31.33	 -0.50  !numbond=1 (P4(O1))
!*GASPARM	op=	 17.07	 13.79	  0.47	 31.33	 -0.67  !numbond=1 (P4(O1,O1))
! GASPARM	s	 10.14	  9.13	  1.38	 20.65	  0.00  !General, numbond=2
! GASPARM	s3	 10.14	  9.13	  1.38	 20.65	  0.00  !General, numbond=2
! GASPARM	s2	 10.88	  9.485	  1.325	 21.69	  0.00  !numbond=1
!*GASPARM	s-1	 10.88	  9.485	  1.325	 21.69	 -1.00  !numbond=1 (C4)
! GASPARM	so	 10.14	  9.13	  1.38	 20.65	  1.00  !numbond=3 (O1)
! GASPARM	so1	 12.00	 10.805	  1.195	 24.00	  1.00  !numbond=4 (O1)
! GASPARM	so2	 12.00	 10.805	  1.195	 24.00	  2.00  !numbond=4 (O1,O1)
! GASPARM	so3	 12.00	 10.805	  1.195	 24.00	  3.00  !numbond=4 (O1,O1,O1)
! GASPARM	so4	 12.00	 10.805	  1.195	 24.00	  4.00  !numbond=4 (O1,O1,O1,O1)
! GASPARM	f	 14.66	 13.85	  2.31	 30.82	  0.00  !numbond=1
! GASPARM	cl	 11.00	  9.69	  1.35	 22.04	  0.00  !numbond=1
! GASPARM	br	 10.08	  8.47	  1.16	 19.71	  0.00  !numbond=1
! GASPARM	i	  9.90	  7.96	  0.96	 18.82	  0.00  !numbond=1
! GASPARM	p	  8.90	  8.24	  0.96	 18.10	  0.00  !General P
!*GASPARM	p#	  8.90	  8.24	  0.96	 18.10	  1.00  !numbond=4 (O1,O1,O1,O1)
!*GASPARM	p=	  8.90	  8.24	  0.96	 18.10	  0.01  !numbond=4 (O1,O1,O1)
!*GASPARM	pn	  8.90	  8.24	  0.96	 18.10	  0.00  !numbond=4 (O1,O1)
!Note: It is not necessary to consider p=, pn, p#, op#, op op=. &
!I tested H3PO4, H2CH3PO4, Antechamber identify all phosporous as p, &
!hydroxyl oxygen is identified as o3, =O as o2, the result is identical to Antechamber. It seems that only when P connected to multiple =O &
!the case will be come more complicated, but this situation is quite rare and thus not needed to consider

call genconnmat(1,1)
write(*,*)

nbond=0
charge=0
!Assign a,b,c parameters to each atom
!nbond array is fully constructed after this looping
write(*,*) "Determined parameters:"
do iatm=1,ncenter
    do jatm=1,ncenter
        if (connmat(iatm,jatm)/=0) then
            nbond(iatm)=nbond(iatm)+1
            if (nbond(iatm)>maxbond) then
                write(*,"(a,i5,a)") " Error: Number of bonds of atom",iatm," exceeded 4, in this case Gasteiger charge cannot be calculated!"
                write(*,*) "Press ENTER button to return"
                read(*,*)
                return
            end if
            bondlist(iatm,nbond(iatm))=jatm
        end if
    end do
    eleidx=a(iatm)%index
    atom_a(iatm)=parma(eleidx,nbond(iatm))
    atom_b(iatm)=parmb(eleidx,nbond(iatm))
    atom_c(iatm)=parmc(eleidx,nbond(iatm))
    
    !Special case: Sp2 N with three bonds, i.e. "GASPARM na"
    !Determine this according to distance of N to the plane defined by the connected three atoms
    !This rule is different to Antechamber, which determines based on complicated connectivity and atomic property, meantime .mol2 must be used
    !Evidently, the current rule is more strict and any format could be used, however the geometry should be optimized first, e.g. PM7
    if (eleidx==7.and.nbond(iatm)==3) then
        i1=bondlist(iatm,1)
        i2=bondlist(iatm,2)
        i3=bondlist(iatm,3)
        dist=potpledis(a(i1)%x,a(i1)%y,a(i1)%z,a(i2)%x,a(i2)%y,a(i2)%z,a(i3)%x,a(i3)%y,a(i3)%z,a(iatm)%x,a(iatm)%y,a(iatm)%z)
        if (dist<0.5D0) then !The nitrogen is close to the local plane, regarding it as sp2
            atom_a(iatm)=12.32D0
            atom_b(iatm)=11.20D0
            atom_c(iatm)=1.34D0
        end if
    end if
    !Special case: N with four bonds, i.e. "GASPARM n4", make it carry +1 charge
    if (eleidx==7.and.nbond(iatm)==4) charge(iatm)=1D0
    !Special case: O only connected to sulfur, i.e. "GASPARM os"
    if (eleidx==8.and.nbond(iatm)==1) then
        if (a(bondlist(iatm,1))%index==16) charge(iatm)=-1D0
    end if
    !Special case, -COO, make each O has initial charge of -1 and C has +1
    if (eleidx==6.and.nbond(iatm)==3) then
		!Check if there are two terminal oxygens
		ncheck=0
		do itmp=1,3
			idx=bondlist(iatm,itmp)
			if (a(idx)%index==8.and.count(connmat(:,idx)/=0)==1) ncheck=ncheck+1
        end do
        if (ncheck==2) then
			charge(iatm)=1D0
			do itmp=1,3
				idx=bondlist(iatm,itmp)
				if (a(idx)%index==8) charge(idx)=-1D0
            end do
        end if
    end if
    
    write(*,"(i5,'(',a,')  numbond=',i2,'   a=',f8.3,'   b=',f8.3,'   c=',f8.3)") &
    iatm,ind2name(eleidx),nbond(iatm),atom_a(iatm),atom_b(iatm),atom_c(iatm)
    if (atom_a(iatm)+atom_b(iatm)+atom_c(iatm)==0) then
        write(*,*) "Error: Parameter is missing for this atom!"
        write(*,*) "Press ENTER button to exit"
        read(*,*)
        return
    end if
end do

!Special case: Assign initial charge for S according to the number of coordinates and connected =O atoms
do iatm=1,ncenter
    if (a(iatm)%index/=16) cycle
    noxy=0 !Count the number of connected =O atoms
    do jdx=1,nbond(iatm)
        jatm=bondlist(iatm,jdx)
        if (a(jatm)%index==8.and.nbond(jatm)==1) noxy=noxy+1
    end do
    if (nbond(iatm)==3.and.noxy==1) then !"GASPARM so"
        charge(iatm)=1D0
    else if (nbond(iatm)==4) then !"GASPARM so1,so2,so3,so4"
        charge(iatm)=noxy
    end if
end do

write(*,*)
inquire(file="PEOEinit.txt",exist=alive)
if (alive) then
	write(*,"(a)") " PEOEinit.txt has been found in current folder, read initial charges from it to override default ones? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=="Y") then
		open(10,file="PEOEinit.txt",status="old")
		nread=0
		do while(.true.)
			read(10,*,iostat=ierror) iatm,tmp
			if (ierror/=0) exit
			charge(iatm)=tmp
			nread=nread+1
		end do
		close(10)
		write(*,"(i8,a)") nread," initial charges have been loaded"
    end if
else
	write(*,"(a)") " Note: If you want to manually set initial charges, you can prepare PEOEinit.txt &
    &in current folder, see Section 3.9.17 of manual"
end if

write(*,*)
if (any(charge/=0)) then
	write(*,*) "Nonzero initial charges:"
	do iatm=1,ncenter
		if (charge(iatm)/=0) write(*,"(i5,'(',a,')  q=',f12.6)") iatm,ind2name(a(iatm)%index),charge(iatm)
	end do
else
	write(*,"(a)") " All initial charges are zero"
end if

!Calculate atom electronegativity at q=1 state
do iatm=1,ncenter
    if (a(iatm)%index==1) then !X(q=1) of hydrogen is special, as mentioned in right side of page 3220 of Tetrahedron, 36, 3219 (1980)
        elenegQ1(iatm) = 20.02D0
    else
        elenegQ1(iatm) = atom_a(iatm) + atom_b(iatm) + atom_c(iatm)
    end if
end do

fdamp=0.5D0
convcrit=0.0001D0 !This criterion is more strict than OpenBabel and Avogadro. It seems that their criteria are about 0.002
maxcyc=50
if (outmedinfo==0) write(*,"(/,a)") " Note: If you want to print atomic charges, amount of charge transfers and &
&atomic electronegativity in each cycle, set ""outmedinfo"" in settings.ini to 1"
write(*,*)
write(*,"(' Max cycles:',i3,'  Charge convergence criterion:',f8.5,'  Damping factor:',f6.3)") maxcyc,convcrit,fdamp
write(*,*)
do icyc=1,maxcyc
    write(*,"(' Cycle',i5)",advance="no") icyc
    
    !Calculate atom electronegativity
    do iatm=1,ncenter
        eleneg(iatm) = atom_a(iatm) + atom_b(iatm)*charge(iatm) + atom_c(iatm)*charge(iatm)**2
    end do
    
    !Calculate transfered charge
    delta_q=0
    do iatm=1,ncenter
        do jatm=1,ncenter
            if (connmat(iatm,jatm)/=0) then
                if (eleneg(jatm)>eleneg(iatm)) then
                    delta_q(iatm)=delta_q(iatm) + (eleneg(jatm)-eleneg(iatm))/elenegQ1(iatm)
                else
                    delta_q(iatm)=delta_q(iatm) + (eleneg(jatm)-eleneg(iatm))/elenegQ1(jatm)
                end if
            end if
        end do
    end do
    delta_q=delta_q*fdamp**icyc
    
    deltamax=maxval(abs(delta_q))
    write(*,"('    Maximum change of charges:',f12.6)") deltamax
    charge=charge+delta_q
    if (outmedinfo==1) then
        do iatm=1,ncenter
	        write(*,"(i6,'(',a,')  Δq=',f10.6,'  X=',f10.3,' eV  q=',f12.6)") iatm,ind2name(a(iatm)%index),delta_q(iatm),eleneg(iatm),charge(iatm)
        end do
    end if
    if (deltamax<convcrit) exit
end do
if (icyc==maxcyc+1) then
    write(*,*) "Error: The convergence is failed!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
else
    write(*,"(' Convergence succeeded after',i4,' cycles')") icyc
    write(*,*)
    write(*,*) "   Atom        Charge"
    do iatm=1,ncenter
	    write(*,"(i6,'(',a,')',f14.8)") iatm,ind2name(a(iatm)%index),charge(iatm)
    end do
    write(*,"(/,a,f12.6)") " Total charge:",sum(charge(:))
	if (allocated(frag1)) then
		write(ides,"(/,' Fragment charge:',f14.8)") sum(charge(frag1))
	end if
    call outatmchg(10,charge(:))
end if

end subroutine




!------ Normalize inputted atomic charges
!Normalization makes sum of atomic populations equal to actual number of electrons. Correct a%charge and nelec must be available
!This is mainly used for getting rid of integration inaccuracy of the atomic charges derived by real space integration
subroutine normalize_atmchg(charge)
use defvar
implicit real*8 (a-h,o-z)
real*8 charge(ncenter)
write(*,*)
write(*,*) "Normalize sum of atomic populations to actual number of electrons..."
totnumelec=sum(a%charge-charge)
facnorm=nelec/totnumelec
do iatm=1,ncenter
	charge(iatm)=a(iatm)%charge-facnorm*(a(iatm)%charge-charge(iatm))
end do
end subroutine


!---- Print atomic charges on screen
subroutine printatmchg(charge)
use defvar
real*8 charge(ncenter)
write(*,*)
write(*,*) "Final atomic charges:"
do iatm=1,ncenter
	write(*,"(' Atom',i5,'(',a2,')',': ',f14.8)") iatm,a(iatm)%name,charge(iatm)
end do
end subroutine



!---- Output atomic charges to a .chg file
!ifileid is the file id that can be used in this subroutine
subroutine outatmchg(ifileid,charge)
use defvar
use util
integer ifileid,i
real*8 charge(ncenter)
character selectyn,chgfilename*200

call path2filename(firstfilename,chgfilename)
write(*,"(/,a)") " If outputting atom coordinates with charges to "//trim(chgfilename)//".chg in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=="y".or.selectyn=="Y") then
	open(ifileid,file=trim(chgfilename)//".chg",status="replace")
	do i=1,ncenter
		write(ifileid,"(a,3f12.6,f15.10)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,charge(i)
	end do
	close(ifileid)
	write(*,"(a)") " Result have been saved to "//trim(chgfilename)//".chg in current folder"
	write(*,"(a)") " Columns 1 to 5 are name,X,Y,Z,charge respectively, coordinates are in Angstrom"
end if
end subroutine



!---- Output all calculated charges (including additional fitting centers) to a .chg file
!ifileid is the file id that can be used in this routine
subroutine outallchg(ifileid,charge,fitcen,nfitcen)
use defvar
use util
integer ifileid,i,nfitcen
real*8 charge(nfitcen),fitcen(3,nfitcen)
character selectyn,chgfilename*200
call path2filename(firstfilename,chgfilename)
write(*,"(a)") " If outputting coordinates of all fitting centers with their charges to "//trim(chgfilename)//".chg in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=="y".or.selectyn=="Y") then
	open(ifileid,file=trim(chgfilename)//".chg",status="replace")
    do icen=1,nfitcen
		if (icen<=ncenter) then
			write(10,"(a,3f12.6,f15.10)") a(icen)%name,fitcen(:,icen)*b2a,charge(icen)
		else
			write(10,"(a,3f12.6,f15.10)") "X ",fitcen(:,icen)*b2a,charge(icen)
		end if
	end do
	close(ifileid)
	write(*,"(a)") " Result have been saved to "//trim(chgfilename)//".chg in current folder"
	write(*,"(a)") " Columns 1 to 5 are name,X,Y,Z,charge respectively, unit is Angstrom"
    write(*,*) "X correspond to additional fitting centers"
end if
end subroutine
