!-------- Main interface of various energy decomposition analyses
subroutine EDA_main
use defvar
implicit real*8 (a-h,o-z)
do while(.true.)
	write(*,*)
	write(*,*) "           ============ Energy decomposition analysis ============ "
	write(*,*) "0 Return"
	write(*,*) "1 Energy decomposition analysis based on molecular forcefield (EDA-FF)"
	write(*,*) "2 Shubin Liu's energy decomposition analysis (Gaussian is needed)"
	write(*,*) "3 sobEDA and sobEDAw energy decomposition analysis"
! 		write(*,*) "2 Mayer energy decomposition analysis"
! 		write(*,*) "3 Fuzzy space based energy decomposition analysis"
	read(*,*) infuncsel2
	if (infuncsel2==0) then
		return
	else if (infuncsel2==1) then
		call EDA_forcefield
	else if (infuncsel2==2) then	
		call EDA_SBL
	else if (infuncsel2==3) then
		write(*,"(a)") " This kind of analysis needs using shell script. Please check detailed sobEDA/sobEDAw tutorial: http://sobereva.com/soft/sobEDA_tutorial.zip"
        write(*,*) "Also see original paper: https://doi.org/10.26434/chemrxiv-2023-n79rz"
        write(*,*) "Press ENTER to return"
		read(*,*)
	end if
end do
end subroutine








!==============================================================================
!==== Energy decomposition analysis based on molecular forcefield (EDA-FF) ====
!==============================================================================

!------ Module for EDA-FF
module EDA_FF_mod
integer :: ielemode=1,ivdwmode=2
real*8,allocatable :: parmA(:),parmB(:)
end module

subroutine EDA_forcefield
use defvar
use util
use EDA_FF_mod
implicit real*8 (a-h,o-z)
integer :: nfrag=0,ishowatmpair=0,ioutatmpqr=0
character c200tmp*200,c200tmp2*200,c2000tmp*2000
integer,allocatable :: frag(:,:),fragnatm(:)
real*8,allocatable :: elemat(:,:),repmat(:,:),dispmat(:,:),totmat(:,:) !Interfragment interaction matrix of electrostatic, repulsive, dispersion and total interaction
real*8 eleatm(ncenter),repatm(ncenter),dispatm(ncenter),totatm(ncenter)
!vdW parameter of each atom, for UFF they are well depths and vdW distance; For AMBER/GAFF they are well depths and vdW radii
character(len=2) c2tmp,FFtype(ncenter) !Force field atom types, only needed by AMBER and GAFF

FFtype="?"
if (ifiletype/=4) a%charge=0 !If the input file is not chg or pqr format, assume atomic charges to be zero
allocate(parmA(ncenter),parmB(ncenter))

10 do while(.true.)
	write(*,*)
    write(*,"(a)") " If this analysis is employed in your work, please cite the following paper, which briefly described and employed EDA-FF analysis"
    write(*,"(a)") "   Mat. Sci. Eng. B-Adv., 273, 115425 (2021) DOI: 10.1016/j.mseb.2021.115425"
    write(*,*)
	write(*,*) "---- Energy decomposition analysis based on molecular forcefield (EDA-FF) ----"
	if (ioutatmpqr==1) write(*,*) "-4 Toggle if outputting atom contributions to .pqr files, current: Yes"
	if (ioutatmpqr==0) write(*,*) "-4 Toggle if outputting atom contributions to .pqr files, current: No"
	if (ishowatmpair==1) write(*,*) "-3 Toggle if outputting atom pairwise interactions to interatm.txt,current:Yes"
	if (ishowatmpair==0) write(*,*) "-3 Toggle if outputting atom pairwise interactions to interatm.txt, current:No"
	if (ielemode==1) write(*,*) "-2 Select model for electrostatic interaction, current: 1/r"
	if (ielemode==2) write(*,*) "-2 Select model for electrostatic interaction, current: 1/r^2"
	if (ivdwmode==1) write(*,*) "-1 Select model for van der Waals interaction, current: UFF"
	if (ivdwmode==2) write(*,*) "-1 Select model for van der Waals interaction, current: AMBER99 & GAFF"
	write(*,*) "0 Return"
	write(*,*) "1 Start analysis"
	if (allocated(frag)) then
		write(*,"(a,i4,a)") " 2 Redefine fragments, current:",nfrag," fragments"
	else
		write(*,"(a)") " 2 Define fragments (undefined currently)"
	end if
	if (ivdwmode==1) then !UFF
		write(*,*) "3 Load atomic charges for current system"
	else !AMBER&GAFF
		write(*,*) "3 Load atomic types and charges for current system"
	end if
	if (ivdwmode==1) then !UFF
		write(*,*) "4 Show current atomic charges"
	else !AMBER&GAFF
		write(*,*) "4 Show current atomic types and charges"
	end if
	read(*,*) isel
	
	if (isel==-4) then
		if (ioutatmpqr==1) then
			ioutatmpqr=0
		else
			ioutatmpqr=1
		end if
	else if (isel==-3) then
		if (ishowatmpair==1) then
			ishowatmpair=0
		else
			ishowatmpair=1
		end if
	else if (isel==-2) then
		write(*,*) "1: 1/r interaction potential"
		write(*,*) "2: 1/r^2 interaction potential"
		read(*,*) ielemode
	else if (isel==-1) then
		write(*,*) "1: UFF van der Waals model"
		write(*,*) "2: AMBER99 & GAFF van der Waals model"
		read(*,*) ivdwmode
		
	else if (isel==0) then
		return
		
	else if (isel==1) then
		if (nfrag==0) then
			write(*,*) "Error: You should use option 2 to define fragments first!"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		exit
		
	else if (isel==2) then !Define fragments
		if (allocated(frag)) deallocate(frag,fragnatm)
		write(*,*) "How many fragments to be defined? e.g. 3"
		write(*,"(a)") " Note: If you input 0, then fragment definition will be loaded from an external plain text file"
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
		else
			allocate(frag(nfrag,ncenter),fragnatm(nfrag))
			do ifrag=1,nfrag
				write(*,"(' Input atomic indices for fragment',i4,', e.g. 1,4,8-12,15')") ifrag
				read(*,"(a)") c2000tmp
				call str2arr(c2000tmp,fragnatm(ifrag),frag(ifrag,:))
			end do
		end if
		
	else if (isel==3) then !Load atomic charges and types
		write(*,*) "Input the path of the file containing list of molecule data"
		write(*,*) "e.g. C:\hanoko\water2_ethanol3.txt"
		write(*,"(a)") " Note: If pressing ENTER button directly, mollist.txt in current folder will be loaded"
		do while(.true.)
			read(*,"(a)") c200tmp
			if (c200tmp==" ") c200tmp="mollist.txt"
			inquire(file=c200tmp,exist=alive)
			if (alive) then
				iatm=0 !Position of the atom in whole system
				open(10,file=c200tmp,status="old")
				do while(.true.) !Read each line of fragment list
					read(10,"(a)",iostat=ierror) c200tmp
					if (c200tmp==" ".or.ierror/=0) exit
					isep=index(trim(c200tmp),' ',back=.true.)
                    if (isep==0) then
                        write(*,*) "Error: You didn't specify number of molecules in the molecular list file!"
                        write(*,*) "Press ENTER button to exit"
                        read(*,*)
                        stop
                    end if
					read(c200tmp(isep+1:),*) nthisfrag
					c200tmp=c200tmp(:isep-1)
					inquire(file=c200tmp,exist=alive)
					if (.not.alive) then
						write(*,"(' Error: Unable to find ',a,', which was specified in the loaded molecular list file')") trim(c200tmp)
						write(*,*) "Press ENTER button to return"
						read(*,*)
						close(10)
						goto 10
					end if
					if (ivdwmode==1) then !UFF
						write(*,"(' Loading atomic charges from ',a,' ...')") trim(c200tmp)
					else !AMBER99 & GAFF
						write(*,"(' Loading atomic charges and types from ',a,' ...')") trim(c200tmp)
					end if
					open(11,file=c200tmp,status="old")
					!Test how many atoms in this fragment .chg file
					natmthis=0
					do while(.true.)
						read(11,*,iostat=ierror) c200tmp2
						if (c200tmp2==" ".or.ierror/=0) exit
						natmthis=natmthis+1
					end do
					!Fill atom charge with/without atom type of the fragment into current system
					do itime=1,nthisfrag
						rewind(11)
						do jatm=1,natmthis
							iatm=iatm+1
							if (iatm>ncenter) then
								write(*,"(a)") " Error: The index of the atom to be loaded exceeded actual number of atoms in the whole system! &
								Please carefully check your molecular list and molecule type file"
								write(*,*) "Press ENTER button to cancel loading"
								read(*,*)
								goto 10
							end if
							if (ivdwmode==1) then !UFF
								read(11,*,iostat=ierror) a(iatm)%charge
							else !AMBER99 & GAFF
								read(11,*,iostat=ierror) FFtype(iatm),a(iatm)%charge
							end if
							if (ierror/=0) then
								write(*,"(a,i5,a)") " Error encountered while loading atom",jatm," in this file!"
								write(*,*) "Press ENTER button to cancel loading"
								read(*,*)
								goto 10
							end if
						end do
					end do
					close(11)
				end do
				close(10)
				exit
			else
				write(*,*) "Unable to find the file, input again"
			end if
		end do
		write(*,*) "Loading finished!"
		if (iatm/=ncenter) write(*,"(' Warning: Only',i5,' atomic charges have been loaded, however there are'i5' atoms in current system!')") iatm,ncenter
		
	else if (isel==4) then !Show current atomic charges and types
		if (allocated(frag)) then
			do ifrag=1,nfrag
				write(*,"(' *** Fragment',i4,':')") ifrag
				do idx=1,fragnatm(ifrag)
					iatm=frag(ifrag,idx)
					if (ivdwmode==1) then
						c2tmp="UF"
					else
						c2tmp=FFtype(iatm)
					end if
					write(*,"(' Atom:',i5,'(',a,')    Charge:',f12.6,'    Type: ',a)") iatm,ind2name(a(iatm)%index),a(iatm)%charge,c2tmp
				end do
			end do
		else
			write(*,*) "No fragment has been defined!"
			do iatm=1,ncenter
				if (ivdwmode==1) then
					c2tmp="UF"
				else
					c2tmp=FFtype(iatm)
				end if
				write(*,"(' Atom:',i5,'(',a,')    Charge:',f12.6,'    Type: ',a)") iatm,ind2name(a(iatm)%index),a(iatm)%charge,c2tmp
			end do
		end if
	end if
	
end do

!================== Now start analysis!
allocate(elemat(nfrag,nfrag),repmat(nfrag,nfrag),dispmat(nfrag,nfrag),totmat(nfrag,nfrag))

!! Setup vdW parameters
call setvdWparm(ivdwmode,FFtype,parmA,parmB,istatus)
if (istatus==1) then !Some parameters are missing
	deallocate(elemat,repmat,dispmat,totmat)
	goto 10
end if

!! Evaluate interactions between fragments
elemat=0;repmat=0;dispmat=0
eleatm=0;repatm=0;dispatm=0
do ifrag=1,nfrag
	do jfrag=ifrag+1,nfrag
		do idx=1,fragnatm(ifrag)
			iatm=frag(ifrag,idx)
			do jdx=1,fragnatm(jfrag)
				jatm=frag(jfrag,jdx)
                call calcAAele(iatm,jatm,eleval)
				elemat(ifrag,jfrag)=elemat(ifrag,jfrag)+eleval
                eleatm(iatm)=eleatm(iatm)+eleval
                eleatm(jatm)=eleatm(jatm)+eleval
                call calcAAvdW(iatm,jatm,repval,dispval)
                repmat(ifrag,jfrag)=repmat(ifrag,jfrag)+repval
                dispmat(ifrag,jfrag)=dispmat(ifrag,jfrag)+dispval
                repatm(iatm)=repatm(iatm)+repval
                repatm(jatm)=repatm(jatm)+repval
                dispatm(iatm)=dispatm(iatm)+dispval
                dispatm(jatm)=dispatm(jatm)+dispval
			end do
		end do
		elemat(jfrag,ifrag)=elemat(ifrag,jfrag)
        repmat(jfrag,ifrag)=repmat(ifrag,jfrag)
        dispmat(jfrag,ifrag)=dispmat(ifrag,jfrag)
	end do
end do
eleatm=eleatm/2
repatm=repatm/2
dispatm=dispatm/2
totmat=elemat+repmat+dispmat
totatm=eleatm+repatm+dispatm

!! Final stage: Output result
write(*,*) "Note: All energies shown below are in kJ/mol!"

!Output interatomic total interaction energies
if (ishowatmpair==1) then
	open(10,file="interatm.txt",status="replace")
	write(10,*) "Note: All energies shown below are in kJ/mol!"
	do ifrag=1,nfrag
		do jfrag=ifrag+1,nfrag
			write(10,"(/,' ******* Between fragment',i4,' and fragment',i4,':')") ifrag,jfrag
			write(10,*) " Atom_i  Atom_j  Dist(Ang) Electrostatic   Repulsive    Dispersion     Total"
			do idx=1,fragnatm(ifrag)
				iatm=frag(ifrag,idx)
				do jdx=1,fragnatm(jfrag)
					jatm=frag(jfrag,jdx)
                    call calcAAele(iatm,jatm,eleval)
                    call calcAAvdW(iatm,jatm,repval,dispval)
					write(10,"(2i7,':  ',f8.3,4f13.2)") iatm,jatm,atomdistA(iatm,jatm,0),eleval,repval,dispval,eleval+repval+dispval
				end do
			end do
		end do
	end do
	close(10)
	write(*,"(a)") " Interatomic total interaction energies have been outputted to interatm.txt in current folder!"
end if

!Output contribution of each atom to total interfragment interaction energy
write(*,"(/,a)") " Contribution of each atom in defined fragments to overall interfragment interaction energies:"
do ifrag=1,nfrag
	do idx=1,fragnatm(ifrag)
		iatm=frag(ifrag,idx)
		write(*,"(' Atom',i5,'(',a,')   Elec:',f8.2,'  Rep:',f8.2,'  Disp:',f8.2,'  Total:',f8.2)") &
        iatm,ind2name(a(iatm)%index),eleatm(iatm),repatm(iatm),dispatm(iatm),totatm(iatm)
	end do
end do

!Output .pqr files
if (ioutatmpqr==1) then
	open(10,file="atmint_tot.pqr",status="replace")
	write(10,"('REMARK   Generated by Multiwfn, totally',i10,' atoms')") ncenter
	do i=1,ncenter
		write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.4,f9.4,a2)") &
		"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,totatm(i),vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
	end do
	write(10,"('END')")
	open(10,file="atmint_ele.pqr",status="replace")
	write(10,"('REMARK   Generated by Multiwfn, totally',i10,' atoms')") ncenter
	do i=1,ncenter
		write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.4,f9.4,a2)") &
		"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,eleatm(i),vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
	end do
	write(10,"('END')")
	open(10,file="atmint_rep.pqr",status="replace")
	write(10,"('REMARK   Generated by Multiwfn, totally',i10,' atoms')") ncenter
	do i=1,ncenter
		write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.4,f9.4,a2)") &
		"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,repatm(i),vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
	end do
	write(10,"('END')")
	open(10,file="atmint_disp.pqr",status="replace")
	write(10,"('REMARK   Generated by Multiwfn, totally',i10,' atoms')") ncenter
	do i=1,ncenter
		write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.4,f9.4,a2)") &
		"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,dispatm(i),vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
	end do
	write(10,"('END')")
	open(10,file="atmint_vdW.pqr",status="replace")
	write(10,"('REMARK   Generated by Multiwfn, totally',i10,' atoms')") ncenter
	do i=1,ncenter
		write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.4,f9.4,a2)") &
		"HETATM",i,' '//ind2name_up(a(i)%index)//' ',"MOL",'A',1,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a,repatm(i)+dispatm(i),vdwr(a(i)%index)*b2a,adjustr(ind2name_up(a(i)%index))
	end do
	write(10,"('END')")
	write(*,*)
	write(*,"(a)") " atmint_tot.pqr, atmint_ele.pqr, atmint_rep.pqr, atmint_disp.pqr, atmint_vdW.pqr have been exported to current folder. &
	Their atomic charge fields record atom contribution to total/electrostatic/repulsive/dispersion/vdW interaction energy between all fragments, &
	respectively. the radius column corresponds to Bondi vdW radii"
	close(10)
end if

!Output interfragment interaction energies
write(*,*)
write(*,*) "Interaction energy components between all fragments:" 
write(*,*) "                        Electrostatic   Repulsion   Dispersion     Total"
do ifrag=1,nfrag
	do jfrag=ifrag+1,nfrag
		write(*,"(' Frag',i4,' -- Frag',i4,':',4f13.2)") ifrag,jfrag,elemat(ifrag,jfrag),repmat(ifrag,jfrag),dispmat(ifrag,jfrag),totmat(ifrag,jfrag)
	end do
end do

deallocate(elemat,repmat,dispmat,totmat)
goto 10

end subroutine


!!----- Calculate atom-atom electrostatic interaction using forcefield
subroutine calcAAele(iatm,jatm,eleval)
use defvar
use util
use EDA_FF_mod
implicit real*8 (a-h,o-z)
integer iatm,jatm
real*8 eleval
if (ielemode==1) then
	eleval=a(iatm)%charge*a(jatm)%charge/atomdist(iatm,jatm,1)
else
	eleval=a(iatm)%charge*a(jatm)%charge/atomdist(iatm,jatm,1)**2
end if
eleval=eleval*au2kJ
end subroutine

!!----- Calculate atom-atom exchange-repulsion (repval) and dispersion interaction (dispval) using forcefield
subroutine calcAAvdW(iatm,jatm,repval,dispval)
use defvar
use util
use EDA_FF_mod
implicit real*8 (a-h,o-z)
integer iatm,jatm
real*8 repval,dispval
if (ivdwmode==1) then !UFF. The unit obtained here is temporarily in kcal/mol
	Dij=dsqrt(parmA(iatm)*parmA(jatm))*cal2J !Well depth in kJ/mol
	Xij=dsqrt(parmB(iatm)*parmB(jatm)) !vdW distance
else if (ivdwmode==2) then !AMBER99 & GAFF
	Dij=dsqrt(parmA(iatm)*parmA(jatm))*cal2J !Well depth in kJ/mol
	Xij=parmB(iatm)+parmB(jatm) !vdW distance obtained by sum of vdW radius
end if
tmpval=(Xij/atomdistA(iatm,jatm,0))**6
repval=Dij*tmpval**2 !Exchange-repulsion
dispval=-2*Dij*tmpval !Dispersion
end subroutine


!!------ Define UFF parameter, because UFF parameters has more than one use in Multiwfn
subroutine defineUFFparm(UFF_A_in,UFF_B_in)
real*8 :: UFF_A_in(103),UFF_B_in(103)
real*8 :: UFF_A(103)=(/ & !UFF well depth (D) in kcal/mol
0.044D0,0.056D0,& !1,2 (H,He)
0.025D0,0.085D0,0.18D0,0.105D0,0.069D0,0.06D0,0.05D0,0.042D0,& !3~10 (Li~Ne)
0.03D0,0.111D0,0.505D0,0.402D0,0.305D0,0.274D0,0.227D0,0.185D0,& !11~18 (Na~Ar)
0.035D0,0.238D0,0.019D0,0.017D0,0.016D0,0.015D0,0.013D0,0.013D0,0.014D0,0.015D0,0.005D0,0.124D0,& !19~30 (K~Zn)
0.415D0,0.379D0,0.309D0,0.291D0,0.251D0,0.22D0,& !31~36 (Ga~Kr)
0.04D0,0.235D0,0.072D0,0.069D0,0.059D0,0.056D0,0.048D0,0.056D0,0.053D0,0.048D0,0.036D0,0.228D0,& !37~48 (Rb~Cd)
0.599D0,0.567D0,0.449D0,0.398D0,0.339D0,0.332D0,& !49~54 (In~Xe)
0.045D0,0.364D0,0.017D0,0.013D0,0.010D0,0.010D0,0.009D0,0.008D0,0.008D0,0.009D0,0.007D0,0.007D0,0.007D0,0.007D0,0.006D0,0.228D0,0.041D0,& !55~71 (Cs~Lu)
0.072D0,0.081D0,0.067D0,0.066D0,0.037D0,0.073D0,0.080D0,0.039D0,0.385D0,0.680D0,0.663D0,0.518D0,0.325D0,0.284D0,0.248D0,& !72~86 (Hf~Rn)
0.050D0,0.404D0,0.033D0,0.026D0,0.022D0,0.022D0,0.019D0,0.016D0,0.014D0,0.013D0,0.013D0,0.013D0,0.012D0,0.012D0,0.011D0,0.011D0,0.011D0 /) !87~103 (Fr~Lr)
real*8 :: UFF_B(103)=(/ & !UFF vdW distance (x) in Angstrom
2.886D0,2.362D0,& !1,2 (H,He)
2.451D0,2.745D0,4.083D0,3.851D0,3.660D0,3.500D0,3.364D0,3.243D0,& !3~10 (Li~Ne)
2.983D0,3.021D0,4.499D0,4.295D0,4.147D0,4.035D0,3.947D0,3.868D0,& !11~18 (Na~Ar)
3.812D0,3.399D0,3.295D0,3.175D0,3.144D0,3.023D0,2.961D0,2.912D0,2.872D0,2.834D0,3.495D0,2.763D0,& !19~30 (K~Zn)
4.383D0,4.280D0,4.230D0,4.205D0,4.189D0,4.141D0,& !31~36 (Ga~Kr)
4.114D0,3.641D0,3.345D0,3.124D0,3.165D0,3.052D0,2.998D0,2.963D0,2.929D0,2.899D0,3.148D0,2.848D0,& !37~48 (Rb~Cd)
4.463D0,4.392D0,4.420D0,4.470D0,4.50D0,4.404D0,& !49~54 (In~Xe)
4.517D0,3.703D0,3.522D0,3.556D0,3.606D0,3.575D0,3.547D0,3.520D0,3.493D0,3.368D0,3.451D0,3.428D0,3.409D0,3.391D0,3.374D0,3.355D0,3.640D0,& !55~71 (Cs~Lu)
3.141D0,3.170D0,3.069D0,2.954D0,3.120D0,2.840D0,2.754D0,3.293D0,2.705D0,4.347D0,4.297D0,4.370D0,4.709D0,4.750D0,4.765D0,& !72~86 (Hf~Rn)
4.90D0,3.677D0,3.478D0,3.396D0,3.424D0,3.395D0,3.424D0,3.424D0,3.381D0,3.326D0,3.339D0,3.313D0,3.299D0,3.286D0,3.274D0,3.248D0,3.236D0 /) !87~103 (Fr~Lr)
UFF_A_in=UFF_A
UFF_B_in=UFF_B
end subroutine


!!------ Assign van der Waals parameters for all atoms according to their types, used for FF based EDA
!FFtype, parmA and parmB should have dimension of ncenter.
!ivdwmode=1: UFF, =2: AMBER99 & GAFF
!Return value: istatus=0 means all parameters are successfully assigned.  =1: Some parameters are missing
subroutine setvdWparm(ivdwmode,FFtype,parmA,parmB,istatus)
use defvar
use util
implicit real*8 (a-h,o-z)
integer ivdwmode,istatus
real*8 parmA(ncenter),parmB(ncenter)
character(len=2) FFtype(ncenter)
real*8 :: UFF_A(103),UFF_B(103)

!AMBER99 atomic information is taken from parm99.dat of AMBERtools package, the lone pair (LP) type is not taken into account
!OM is not standard AMBER atomic type, however gview assigns carboxyl oxygen as OM, while in standard AMBER naming it should be O
integer,parameter :: nAMBERtype=62
character(len=2) :: AMBERname(nAMBERtype)=(/ "H ","HO","HS","HC","H1","H2","H3","HP","HA",&
"H4","H5","HW","HZ","O ","O2","OW","OH","OS","C*","CT","C ","N ","N3","S ","SH",&
"P ","IM","Li","IP","Na","K ","Rb","Cs","MG","C0","Zn","F ","Cl","Br","I ","IB",&
"NA","N2","N*","NC","NB","NT","NY",&  !7 types equivalent to "N "
"CA","CB","CC","CD","CK","CM","CN","CQ","CR","CV","CW","CY","CZ",&  !13 types equivalent to "C*"
"OM"/)  !Special atomic type assigned by gview. Use the same vdW parameter as "O"

real*8 :: AMBER_A(nAMBERtype)=(/ & !AMBER well depth in kcal/mol
0.0157D0,0.0000D0,0.0157D0,0.0157D0,0.0157D0,0.0157D0,0.0157D0,0.0157D0,0.0150D0,0.0150D0,0.0150D0,&
0.0000D0,0.0150D0,0.2100D0,0.2100D0,0.1520D0,0.2104D0,0.1700D0,0.0860D0,0.1094D0,0.0860D0,0.1700D0,&
0.1700D0,0.2500D0,0.2500D0,0.2000D0,0.1D0,0.0183D0,0.00277D0,0.00277D0,0.000328D0,0.00017D0,&
0.0000806D0,0.8947D0,0.459789D0,0.0125D0,0.061D0,0.265D0,0.320D0,0.40D0,0.1D0,&
0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,&
0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,&
0.2100D0 /)
real*8 :: AMBER_B(nAMBERtype)=(/ & !AMBER vdW radii in Angstrom
0.6000D0,0.0000D0,0.6000D0,1.4870D0,1.3870D0,1.2870D0,1.1870D0,1.1000D0,1.4590D0,1.4090D0,1.3590D0,0.0000D0,1.4590D0,&
1.6612D0,1.6612D0,1.7683D0,1.7210D0,1.6837D0,1.9080D0,1.9080D0,1.9080D0,1.8240D0,1.8240D0,2.0000D0,2.0000D0,2.1000D0,&
2.47D0,1.1370D0,1.8680D0,1.8680D0,2.6580D0,2.9560D0,3.3950D0,0.7926D0,1.7131D0,1.10D0,1.75D0,1.948D0,2.22D0,2.35D0,5.0D0,&
1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,&
1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,&
1.6612D0 /)

!GAFF atomic information is taken from gaff.dat of AMBERtools package
integer,parameter :: nGAFFtype=83
character(len=2) :: GAFFname(nGAFFtype)=(/ &
"h1","h2","h3","h4","h5","ha","hc","hn","ho","hp","hs","hw","hx","o ",&
"oh","os","op","oq","ow","c ","c1","c2","c3","ca","cc","cd","ce","cf",&
"cg","ch","cp","cq","cu","cv","cx","cy","cz","n ","ni","nj","n1","n2",&
"n3","np","nq","n4","nk","nl","na","nb","nc","nd","ne","nf","nh","nm",&
"nn","no","s ","s2","s4","s6","sx","sy","sh","ss","sp","sq","p2","p3",&
"p4","p5","pb","pc","pd","pe","pf","px","py","f ","cl","br","i " /)
character(len=2) c2tmp
real*8 :: GAFF_A(nGAFFtype)=(/ & !GAFF well depth in kcal/mol
0.0157D0,0.0157D0,0.0157D0,0.0150D0,0.0150D0,0.0150D0,0.0157D0,0.0157D0,&
0.0000D0,0.0157D0,0.0157D0,0.0000D0,0.0157D0,0.2100D0,0.2104D0,0.1700D0,&
0.1700D0,0.1700D0,0.1520D0,0.0860D0,0.2100D0,0.0860D0,0.1094D0,0.0860D0,&
0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.2100D0,0.2100D0,0.0860D0,0.0860D0,&
0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.0860D0,0.1700D0,0.1700D0,0.1700D0,&
0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,&
0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,0.1700D0,&
0.1700D0,0.1700D0,0.2500D0,0.2500D0,0.2500D0,0.2500D0,0.2500D0,0.2500D0,&
0.2500D0,0.2500D0,0.2500D0,0.2500D0,0.2000D0,0.2000D0,0.2000D0,0.2000D0,&
0.2000D0,0.2000D0,0.2000D0,0.2000D0,0.2000D0,0.2000D0,0.2000D0,0.061D0,0.265D0,0.420D0,0.50D0 /)
real*8 :: GAFF_B(nGAFFtype)=(/ & !GAFF vdW radii in Angstrom
1.3870D0,1.2870D0,1.1870D0,1.4090D0,1.3590D0,1.4590D0,1.4870D0,0.6000D0,0.0000D0,0.6000D0,0.6000D0,&
0.0000D0,1.1000D0,1.6612D0,1.7210D0,1.6837D0,1.6837D0,1.6837D0,1.7683D0,1.9080D0,1.9080D0,1.9080D0,&
1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.9080D0,&
1.9080D0,1.9080D0,1.9080D0,1.9080D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,&
1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,1.8240D0,&
1.8240D0,1.8240D0,1.8240D0,2.0000D0,2.0000D0,2.0000D0,2.0000D0,2.0000D0,2.0000D0,2.0000D0,2.0000D0,&
2.0000D0,2.0000D0,2.1000D0,2.1000D0,2.1000D0,2.1000D0,2.1000D0,2.1000D0,2.1000D0,2.1000D0,2.1000D0,&
2.1000D0,2.1000D0,1.75D0,1.948D0,2.02D0,2.15D0 /)

!Fill UFF parameter
call defineUFFparm(UFF_A,UFF_B)

if (ivdwmode==1) then
	do iatm=1,ncenter
		parmA(iatm)=UFF_A(a(iatm)%index)
		parmB(iatm)=UFF_B(a(iatm)%index)
	end do
else if (ivdwmode==2) then !AMBER99 & GAFF
atmcyc:	do iatm=1,ncenter
		if (FFtype(iatm)=="UF") then !Borrow UFF parameter
			parmA(iatm)=UFF_A(a(iatm)%index)
			parmB(iatm)=UFF_B(a(iatm)%index)/2
			cycle atmcyc
		end if
		!Check if is AMBER atom type
		ntottype=nAMBERtype
		do itype=1,ntottype
			if (FFtype(iatm)==AMBERname(itype)) then
				parmA(iatm)=AMBER_A(itype)
				parmB(iatm)=AMBER_B(itype)
				cycle atmcyc
			end if
		end do
		!Check if is GAFF atom type
		ntottype=nGAFFtype
		do itype=1,ntottype
			if (FFtype(iatm)==GAFFname(itype)) then
				parmA(iatm)=GAFF_A(itype)
				parmB(iatm)=GAFF_B(itype)
				cycle atmcyc
			end if
		end do
        !Final check. GaussView always uses uppercase, if the second character is convert to lowercase, may be assignment can be successful
		ntottype=nAMBERtype
		do itype=1,ntottype
            c2tmp=FFtype(iatm)
            call uc2lc(c2tmp(2:2))
			if (c2tmp==AMBERname(itype)) then
				parmA(iatm)=AMBER_A(itype)
				parmB(iatm)=AMBER_B(itype)
				cycle atmcyc
			end if
		end do
		write(*,"(' Error: Unable to find parameter for atom',i5,' with type ',a)") iatm,FFtype(iatm)
		write(*,"(a)") " If the element is indeed defined by the forcefield, please check if spelling and capitalization is correct"
		write(*,*) "Press ENTER button to cancel calculation"
		read(*,*)
		istatus=1
		return
	end do atmcyc
end if

istatus=0
end subroutine





!====================================================
!==== Shubin Liu's energy decomposition analysis ====
!====================================================
subroutine EDA_SBL
use defvar
use util
implicit real*8 (a-h,o-z)
character c200tmp*200

write(*,*) "Citation: J. Chem. Phys., 126, 244103 (2007)"
write(*,*)

write(*,*) "Input path of Gaussian output file with ""ExtraLinks=L608"", e.g. C:\Yohane.out"
do while(.true.)
	read(*,"(a)") c200tmp
	inquire(file=c200tmp,exist=alive)
	if (alive) exit
	write(*,*) "Cannot find the file, input again"
end do

open(10,file=c200tmp,status="old")
call loclabel(10,"ET=",ifound)
if (ifound==0) then
	write(*,"(a)") " Error: Unable to find ""ET="" term, please double check your Gaussian input file and Multiwfn manual"
	write(*,*) "Press ENTER button to exit"
	read(*,*)
	return
end if

read(10,"(a)") c200tmp
read(c200tmp(6:),*) ET
read(c200tmp(37:),*) EJ
read(c200tmp(71:),*) ENuc
!Frequently, the EV is quite large making the corresponding output is *****. Therefore, EV will be obtained as EV=ENTVJ-ET-EJ-ENuc

call loclabel(10,"Ex= ")
read(10,"(18x,f12.6,4x,f12.6,4x,f12.6)") ENTVJ,Ex,Ec
EV=ENTVJ-ET-EJ-ENuc
close(10)
write(*,*) "Data has been successfully loaded from this file"

write(*,*)
write(*,*) "Calculating E_steric (Weizsacker kinetic energy)"
iuserfunc_old=iuserfunc
iuserfunc=5
call intfunc_silent(100,TW)
iuserfunc=iuserfunc_old

E_elst=EJ+EV+ENuc
E_quan=ET-TW+Ex+Ec
E_tot=TW+E_elst+E_quan !Equivalent to ET+Ex+Ec+EJ+EV+ENuc
write(*,"(/,a,f16.6,' Hartree')") " Electronic kinetic energy (ET):",ET
write(*,"(a,f16.6,' Hartree')") " Weizsacker kinetic energy (TW):",TW
write(*,"(a,f16.6,' Hartree')") " Interelectronic Coulomb repulsion energy (EJ):",EJ
write(*,"(a,f16.6,' Hartree')") " Internuclear Coulomb repulsion energy (ENuc):",ENuc
write(*,"(a,f16.6,' Hartree')") " Nuclear-electronic Coulomb attraction energy (EV):",EV
write(*,"(a,f16.6,' Hartree')") " Energy without electronic correlation (ET+EV+EJ+ENuc):",ET+EV+EJ+ENuc !Corresponding to the "ENTVJ" term in output file
write(*,"(a,f16.6,' Hartree')") " Exchange correlation energy (Ex):",Ex
write(*,"(a,f16.6,' Hartree')") " Coulomb correlation energy (Ec):",Ec
write(*,"(a,f16.6,' Hartree')") " Pauli kinetic energy (ET-TW):",ET-TW
write(*,*)
write(*,*) "----- EDA-SBL energy decomposition terms:"
write(*,"(a,f16.6,' Hartree')") " E_steric:       ",TW
write(*,"(a,f16.6,' Hartree')") " E_electrostatic:",E_elst
write(*,"(a,f16.6,' Hartree')") " E_quantum:      ",E_quan
write(*,"(/,a,f16.6,' Hartree')") " E_total:        ",E_tot

!I found for G16 A.03 with M06-2X, the E_tot shown above, which is exactly equals to "ETot=" printed by L608, is remarkably
!different to the single point energy printed in either output file or fch file. By comparing with G09, I found the reason is
!that the "Ec=" shown by G16 with M06-2X is incorrect
!Same problem was found for G03 E.01+B3LYP and G16 A.03+B3LYP, while G09 D.01 works normally with B3LYP and M06-2X
diff=totenergy-E_tot
if (abs(diff)>1D-4) then
	write(*,"(/,a,f14.6,a)") " Warning: The total energy shown above is detectably different to the total energy (",totenergy," a.u.) in fch/fchk file! &
	This issue is known in some versions of Gaussian for certain DFT functionals, the reason is that the &
	correlation energy printed by Link 608 is inaccurate. You are suggested to either try to use other Gaussian version (G09 D.01 seems always gives correct result), or try to use other functionals"
	write(*,*) "Press ENTER button to continue"
	read(*,*)
end if
end subroutine
