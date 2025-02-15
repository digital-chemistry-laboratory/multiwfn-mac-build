!Index of routines:
!orbcomp_main: Main interface of various orbital composition analyses
!orballcomp_MMPA: Interface for calculating basis function, shell and atom contribution via Mulliken, Stout-Politzer, SCPA
!orbfragcomp_MMPA: Interface for calculating fragment contribution via Mulliken, Stout-Politzer, SCPA
!orbatmcomp_space: Interface for calculating atom contribution via fuzzy partition methods
!orballcomp_NAO: Interface for calculating atom contribution via NAOMO method
!gen_orbatmcomp_space: Returning atomic contributions to specific range of orbitals by Hirshfeld (built-in density) or Becke method
!gen_orbatmcomp_MMPA: Returning all atom contributions to one orbital via Mulliken and SCPA
!gen_allorbbascomp_SCPA: Returning all basis function contributions to all orbitals via SCPA
!gen_allorbatmcomp_SCPA: Returning all atom contributions to all orbitals via SCPA
    
    
    

!!--------- Interface of various orbital composition analyses
subroutine orbcomp_main
use defvar
implicit real*8 (a-h,o-z)

do while(.true.)
    write(*,"(/,a)") " If this module is used in your work, please cite the following paper together with &
    &Multiwfn original paper (J. Comput. Chem., 33, 580 (2012) and J. Chem. Phys., 161, 082503 (2024)), as the methods employed in this module are comprehensively described in this paper:"
    write(*,"(a)") " Tian Lu, Feiwu Chen, Acta Chimica Sinica, 69, 2393-2406 (2011) http://sioc-journal.cn/Jwk_hxxb/CN/abstract/abstract340458.shtml"
    write(*,*)
	write(*,*) "       ================ Orbital composition analysis ==============="
	write(*,*) "-10 Return"
	if (allocated(CObasa)) then
		write(*,*) "-2 Define fragment 2 (for option 4,5)"
		write(*,*) "-1 Define fragment 1 (for option 1~6)"
		write(*,*) "1 Orbital composition analysis with Mulliken partition"
		write(*,*) "2 Orbital composition analysis with Stout-Politzer partition"
		write(*,*) "3 Orbital composition analysis with Ros-Schuit (SCPA) partition"
		write(*,*) "4 Print frag. 1 & inter-fragment compositions in all orbitals (Mulliken)"
		write(*,*) "5 Print frag. 1 & inter-fragment compositions in all orbitals (Stout-Politzer)"
		write(*,*) "6 Print frag. 1 compositions in all orbitals (SCPA)"
	end if
	write(*,*) "7 Orbital composition analysis by natural atomic orbital (NAO) method"
	if (allocated(b)) write(*,*) "8 Calculate atom and fragment contributions by Hirshfeld method"
	if (allocated(b)) write(*,*) "9 Calculate atom and fragment contributions by Becke method"
	if (allocated(b)) write(*,*) "10 Calculate atom and fragment contributions by Hirshfeld-I method"
	if (allocated(b)) write(*,*) "11 Calculate atom and fragment contributions by AIM method"
	if (allocated(CObasa)) write(*,*) "100 Evaluate oxidation state by LOBA/mLOBA method"
	read(*,*) icompana

    if (ifPBC/=0.and.(icompana==9.or.icompana==10)) then
		write(*,*) "Error: This analysis does not support periodic wavefunctions!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        cycle
    end if
    
	if (icompana==-10) then
		if (allocated(frag1)) deallocate(frag1)
		if (allocated(frag2)) deallocate(frag2)
		exit
	else if (icompana==-1) then
		call deffrag(1)
	else if (icompana==-2) then
		call deffrag(2)
	else if (icompana==1.or.icompana==2.or.icompana==3) then
        if (icompana==1.or.icompana==2) call ask_Sbas_PBC
		if (icompana==1) call orballcomp_MMPA(1)
		if (icompana==2) call orballcomp_MMPA(2)
		if (icompana==3) call orballcomp_MMPA(3)
	else if (icompana==4.or.icompana==5.or.icompana==6) then
		if (.not.allocated(frag1)) then
			write(*,*) "Error: Please use option -1 to define fragment 1"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
			cycle
		end if
		if (.not.allocated(frag2).and.(icompana==4.or.icompana==5)) write(*,*) "Note: Fragment 2 was not be defined"
        if (icompana==4.or.icompana==5) call ask_Sbas_PBC
		write(*,*)
		if (icompana==4) call orbfragcomp_MMPA(1)
		if (icompana==5) call orbfragcomp_MMPA(2)
		if (icompana==6) call orbfragcomp_MMPA(3)
	else if (icompana==7) then
		call orballcomp_NAO
	else if (icompana==8) then
        call orbatmcomp_space(1)
	else if (icompana==9) then
		call orbatmcomp_space(2)
	else if (icompana==10) then
		call Hirshfeld_I(2)
		call orbatmcomp_space(3)
	else if (icompana==11) then
        write(*,"(a)") " Note: To compute orbital composition under AIM partition, you should use basin analysis module (main function 17). &
        &see Section 4.8.6 of Multiwfn manual for example"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
	else if (icompana==100) then
		call LOBA
	end if
end do
end subroutine



!!-------------- Define fragment 1 & 2 and store to global array frag1 and frag2. "isel" denote define which fragment
!If the fragment is empty, then will not allcoate frag1/2 when exit
subroutine deffrag(isel)
use defvar
use util
implicit real*8 (a-h,o-z)
!fragtmp stores temporary basis index before exit setting interface, by default they are all zero,
!if basis function 4,89,32 are in this fragment, then value 4,89,32 will be in uncertain three slots of fragtmp
integer fragtmp(nbasis)
integer termtmp(nbasis) !Store each time read serial number
integer vectmp(nbasis) !used by "inv" command
integer atmsellist(ncenter),bassellist(nbasis)
character c10000tmp*10000,c80tmp*80,lchar
fragtmp=0
if (isel==1.and.allocated(frag1)) then
	fragtmp(1:size(frag1))=frag1(:)
	deallocate(frag1)
else if (isel==2.and.allocated(frag2)) then
	fragtmp(1:size(frag2))=frag2(:)
	deallocate(frag2)
end if
10 write(*,*) "Commands and examples:"
write(*,*) "q: Save fragment and exit"
write(*,*) "clean: Clean current fragment"
write(*,*) "list: List basis functions in current fragment"
write(*,*) "all: List all basis functions of current system"
write(*,*) "addall: Add all basis functions to fragment"
write(*,*) "inv: Invert selection, viz. delete existing basis functions add all other ones"
write(*,*) "help: Print these help information again"
write(*,*) "cond: Add basis functions satisfying given conditions to fragment"
write(*,*) "a 2,5-8,12: Add all basis functions in atoms 2,5,6,7,8,12 to fragment"
write(*,*) "s 2,5-8,12: Add all basis functions in shells 2,5,6,7,8,12 to fragment"
write(*,*) "b 2,5-8,12: Add basis functions 2,5,6,7,8,12 to fragment"
write(*,*) "l s,d: Add basis functions of s and d angular moments to fragment"
write(*,*) "e Fe: Add all atoms of Fe element to fragment"
write(*,*) "da 5,7,11-13: Delete all basis functions in atoms 5,7,11,12,13 from fragment"
write(*,*) "db 5,7,11-13: Delete basis functions 5,7,11,12,13 from fragment"

do while(.true.)
	read(*,"(a)") c10000tmp
	if (c10000tmp(1:4)=="help") then
		goto 10
	else if (c10000tmp(1:6)=="addall") then
		forall(i=1:nbasis) fragtmp(i)=i
		write(*,*) "Done!"
	else if (c10000tmp(1:1)=="q") then
		numbas=count(fragtmp/=0)
		if (isel==1.and.numbas>=1) allocate(frag1(numbas))
		if (isel==2.and.numbas>=1) allocate(frag2(numbas))
		write(*,*) "Basis function index in current fragment:"
		if (numbas==0) then
			write(*,*) "None"
		else
			do i=1,nbasis
				if (fragtmp(i)/=0) write(*,"(i6)",advance="no") fragtmp(i)
			end do
		end if
		write(*,*)
		write(*,*) "Fragment is saved"
		write(*,*)
		if (numbas>=1) then
			if (isel==1) frag1=pack(fragtmp,fragtmp/=0)
			if (isel==2) frag2=pack(fragtmp,fragtmp/=0)
		end if
		exit
	else if (c10000tmp(1:5)=="clean") then
		fragtmp=0
		write(*,*) "Done!"
	else if (c10000tmp(1:3)=="all") then
		write(*,*) "The basis functions with asterisk are those presented in current fragment"
		do i=1,nbasis
			if (any(fragtmp==i)) then
				write(*,"('* Basis:',i6,'    Shell:',i5,'    Center:',i5,'(',a2,')    Type: ',a)")&
				 i,basshell(i),bascen(i),a(bascen(i))%name,GTFtype2name(bastype(i))
			else
				write(*,"('  Basis:',i6,'    Shell:',i5,'    Center:',i5,'(',a2,')    Type: ',a)")&
				 i,basshell(i),bascen(i),a(bascen(i))%name,GTFtype2name(bastype(i))
			end if
		end do
	else if (c10000tmp(1:4)=="list") then
		write(*,*) "Basis functions in current fragment:"
		if (all(fragtmp==0)) then
			write(*,*) "None"
		else
            ntmp=0
			do i=1,nbasis
				if (fragtmp(i)/=0) then
                    write(*,"(i5)",advance="no") fragtmp(i)
                    ntmp=ntmp+1
                end if
			end do
            write(*,"(/,' Totally',i8,' basis functions')") ntmp
			write(*,*)
		end if
	else if (c10000tmp(1:3)=="inv") then
		vectmp=fragtmp
		fragtmp=0
		itmp=0
		do ibas=1,nbasis
			if (all(vectmp/=ibas)) then
				itmp=itmp+1
				fragtmp(itmp)=ibas
			end if 
		end do
		write(*,*) "Done!"
	else if (c10000tmp(1:4)=="cond") then
		write(*,"(a)") " Note: You will be prompted to input three conditions in turn, &
		&the basis functions satisfying all conditions will be added to current fragment"
		write(*,*)
		write(*,*) "Condition 1: Input range of atoms, e.g. 2,5-8,12"
        write(*,*) "You can also input one element name (case sensitive), e.g. Fe"
		write(*,*) "If you press ENTER button directly, the atom will be arbitrary"
		read(*,"(a)") c10000tmp
		if (c10000tmp==" ".or.index(c10000tmp,'a')/=0) then
			nselatm=ncenter
			forall(i=1:nselatm) atmsellist(i)=i
		else
			if (iachar(c10000tmp(1:1))>=48.and.iachar(c10000tmp(1:1))<=57) then !Inputted number
				call str2arr(c10000tmp,nselatm,atmsellist)
				if (any(atmsellist(1:nselatm)>ncenter)) then
					write(*,*) "Error: One or more atom indices exceeded valid range!"
					cycle
				end if
            else !Inputted element
				nselatm=0
				do iatm=1,ncenter
					if (a(iatm)%name==c10000tmp) then
						nselatm=nselatm+1
                        atmsellist(nselatm)=iatm
                    end if
                end do
            end if
		end if
		write(*,*) "Condition 2: Input range of basis functions, e.g. 2,5-8,12"
		write(*,*) "If you press ENTER button directly, the basis function index will be arbitrary"
		read(*,"(a)") c10000tmp
		if (c10000tmp==" ".or.index(c10000tmp,'a')/=0) then
			nselbas=nbasis
			forall(i=1:nselbas) bassellist(i)=i
		else
			call str2arr(c10000tmp,nselbas,bassellist)
			if (any(bassellist(1:nselbas)>nbasis)) then
				write(*,*) "Error: One or more basis function indices exceeded valid range!"
				cycle
			end if
		end if
		write(*,*) "Condition 3: Choose basis function type"
        write(*,*) "Could be one of such as S, Y, Z, XY, YY, ZZZ, D+1, D 0 ..."
		write(*,*) "You can also choose shell type, one of S, P, D, F, G, H"
		write(*,*) "If you press ENTER button directly, the type will be arbitrary"
		read(*,"(a)") c10000tmp
		naddbas=0
		do ibasidx=1,nselbas !Examine all basis functions
			ibas=bassellist(ibasidx)
			iadd=0
			if ( any(fragtmp==ibas) ) cycle !Skip if already presented in current fragment
			if ( all(atmsellist(1:nselatm)/=bascen(ibas)) ) cycle !Atom index condition
			if ( index(c10000tmp,'a')/=0.or.c10000tmp==" ") then
				iadd=1
			else !Basis function type condition is defined
				itype=bastype(ibas)
				if ( trim(c10000tmp)=='P' .and. ( itype>=2.and.itype<=4 ) ) iadd=1
				if ( trim(c10000tmp)=='D' .and. ( (itype>=-5 .and.itype<=-1 ).or.(itype>=5 .and.itype<=10) ) ) iadd=1
				if ( trim(c10000tmp)=='F' .and. ( (itype>=-12.and.itype<=-6 ).or.(itype>=11.and.itype<=20) ) ) iadd=1
				if ( trim(c10000tmp)=='G' .and. ( (itype>=-21.and.itype<=-13).or.(itype>=21.and.itype<=35) ) ) iadd=1
				if ( trim(c10000tmp)=='H' .and. ( (itype>=-32.and.itype<=-22).or.(itype>=36.and.itype<=56) ) ) iadd=1
				if ( trim(c10000tmp)==GTFtype2name(bastype(ibas)) ) iadd=1 !Inputted is detailed type
			end if
			if (iadd==1) then
				do i=1,nbasis !Find space slot to save this basis function
					if (fragtmp(i)==0) then
						fragtmp(i)=ibas
						naddbas=naddbas+1
						exit
					end if
				end do
			end if
		end do
		write(*,"(' Done!',i6,' new basis functions have been added to current fragment')") naddbas
		
	else if (c10000tmp(1:2)=="a ".or.c10000tmp(1:2)=="s ".or.c10000tmp(1:2)=="b ".or.c10000tmp(1:2)=="db".or.c10000tmp(1:2)=="da") then
		call str2arr(c10000tmp(3:),nterm,termtmp)
		if (c10000tmp(1:2)=="a ") then
			if (any(termtmp(1:nterm)<=0).or.any(termtmp(1:nterm)>ncenter)) then
				write(*,*) "Error: Atom index exceeded valid range! Ignored"
				cycle
			end if
			do iterm=1,nterm
				iatm=termtmp(iterm)
                if (basstart(iatm)==0) cycle
				do ibas=basstart(iatm),basend(iatm)
					if (all(fragtmp/=ibas)) then
						do i=1,nbasis !Find an empty slot to record this basis function
							if (fragtmp(i)==0) then
								fragtmp(i)=ibas
								exit
							end if
						end do
					end if
				end do
			end do
		else if (c10000tmp(1:2)=="s ") then
			do ibas=1,nbasis
				if (any(termtmp(1:nterm)==basshell(ibas)).and.all(fragtmp/=ibas)) then
					do j=1,nbasis !Find an empty slot to record this basis function
						if (fragtmp(j)==0) then
							fragtmp(j)=ibas
							exit
						end if
					end do
				end if
			end do
		else if (c10000tmp(1:2)=="b ") then
			do i=1,nterm
				if (all(fragtmp/=termtmp(i)).and.termtmp(i)<=nbasis.and.termtmp(i)>0) then
					do j=1,nbasis !Find empty slot to save this basis function
						if (fragtmp(j)==0) then
							fragtmp(j)=termtmp(i)
							exit
						end if
					end do
				end if
			end do
		else if (c10000tmp(1:2)=="db") then
			do i=1,nterm
				where(fragtmp==termtmp(i)) fragtmp=0
			end do
		else if (c10000tmp(1:2)=="da") then
			do idx=1,nterm
				iatm=termtmp(idx)
				if (basstart(iatm)==0) cycle
				do ibas=basstart(iatm),basend(iatm)
					where(fragtmp==ibas) fragtmp=0
				end do
			end do
		end if
		write(*,*) "Done!"
    else if (c10000tmp(1:2)=="e ") then
		read(c10000tmp(3:),*) c80tmp
        call elename2idx(c80tmp,idx)
        if (idx==0) then
			write(*,*) "Error: Unable to identify the element!"
        else
			naddatm=0
			do iatm=1,ncenter
				if (a(iatm)%index==idx) then
					naddatm=naddatm+1
					if (basstart(iatm)==0) cycle
					do ibas=basstart(iatm),basend(iatm)
						if (all(fragtmp/=ibas)) then
							do i=1,nbasis !Find an empty slot to record this basis function
								if (fragtmp(i)==0) then
									fragtmp(i)=ibas
									exit
								end if
							end do
						end if
					end do
                end if
			end do
            write(*,"(' Done!',i6,' atoms have been added')") naddatm
        end if
	else if (c10000tmp(1:2)=="l ") then
		naddbas=0
		do ibas=1,nbasis
			ido=0
			itype=bastype(ibas)
			if ( (index(c10000tmp,'s')/=0.or.index(c10000tmp,'S')/=0).and.itype==1 ) ido=1
			if ( (index(c10000tmp,'p')/=0.or.index(c10000tmp,'P')/=0).and.( itype>=2.and.itype<=4) ) ido=1
			if ( (index(c10000tmp,'d')/=0.or.index(c10000tmp,'D')/=0).and.( (itype>=-5 .and.itype<=-1 ).or.(itype>=5 .and.itype<=10) ) ) ido=1
			if ( (index(c10000tmp,'f')/=0.or.index(c10000tmp,'F')/=0).and.( (itype>=-12.and.itype<=-6 ).or.(itype>=11.and.itype<=20) ) ) ido=1
			if ( (index(c10000tmp,'g')/=0.or.index(c10000tmp,'G')/=0).and.( (itype>=-21.and.itype<=-13).or.(itype>=21.and.itype<=35) ) ) ido=1
			if ( (index(c10000tmp,'h')/=0.or.index(c10000tmp,'H')/=0).and.( (itype>=-32.and.itype<=-22).or.(itype>=36.and.itype<=56) ) ) ido=1
			if (ido==1.and.all(fragtmp/=ibas)) then
				do j=1,nbasis !Find an empty slot to record this basis function
					if (fragtmp(j)==0) then
						fragtmp(j)=ibas
						naddbas=naddbas+1
						exit
					end if
				end do
			end if
		end do
		write(*,"(' Done!',i6,' new basis functions have been added to current fragment')") naddbas
	else
		write(*,*) "Error: Unrecognized command! Ignored"
	end if
end do
if (allocated(frag1).and.allocated(frag2)) then
	j=0
	do i=1,size(frag2)
		if (any(frag1==frag2(i))) then
			write(*,"(' Warning: basis function',i6,' in fragment 2 also present in fragment 1!')") frag2(i)
			j=1
		end if
	end do
	if (j==1) then
		write(*,*) "You must redefine fragment 1 or 2 to eliminate intersection set!"
		write(*,*)
	end if
end if
end subroutine



!!--------- fragment and inter-fragment composition analysis
!isel==1 fragment 1 & inter-fragment composition with Mulliken partition"
!isel==2 fragment 1 & inter-fragment composition with Stout & Politzer partition"
!isel==3 fragment 1 composition with Ros & Schuit (SCPA) partition"
subroutine orbfragcomp_MMPA(isel)
use defvar
use util
implicit real*8 (a-h,o-z)
integer isel
real*8,pointer :: tmpmat(:,:)
real*8 :: ovpfrg12(nmo),ovpfrg12_1(nmo)
character orbtype*2
ovpfrg12=0D0
ovpfrg12_1=0D0
if (isel==1.or.isel==2) then
	write(*,"(a)") " ""c^2"": The square of coefficients of all basis functions within fragment 1"
	write(*,"(a)") " ""Int.cross"": Cross term within fragment 1"
	write(*,"(a,/)") " ""Ext.cross"": The fragment 1 part of the total cross term between fragment 1 and all other basis functions"
	write(*,"(' Orb#  Type  Ene(a.u.)    Occ        c^2     Int.cross    Ext.cross     Total')")
else if (isel==3) then
	write(*,"(' Orb#   Type   Energy(a.u.)     Occ       Composition')")
end if

do imo=1,nmo
	if (imo<=nbasis) then
		irealmo=imo
		tmpmat=>CObasa
	else if (imo>nbasis) then
		irealmo=imo-nbasis
		tmpmat=>CObasb
	end if
	floc=0D0 !Local term in fragment 1
	crossint=0D0
	crossext=0D0 !fragment 1 part of all overlap between fragment 1 and external atoms
	if (isel==3) allsqr=sum(tmpmat(:,irealmo)**2)
    if (all(CObasa(:,imo)==0)) cycle !This orbital should be an artificially filled one
	do i=1,size(frag1)
		ibas=frag1(i)
		if (isel==1.or.isel==2) floc=floc+tmpmat(ibas,irealmo)**2
		if (isel==3) floc=floc+tmpmat(ibas,irealmo)**2/allsqr
		!Calculate cross term
		if (isel==1.or.isel==2) then
			do jbas=1,nbasis
				if (jbas==ibas) cycle
				crossij=tmpmat(ibas,irealmo)*tmpmat(jbas,irealmo)*Sbas(ibas,jbas)
				if (any(frag1==jbas)) then
					crossint=crossint+crossij !Internal cross term
				else if (isel==1) then
					crossext=crossext+crossij
				else if (isel==2) then
					tmpdenom=tmpmat(ibas,irealmo)**2+tmpmat(jbas,irealmo)**2
					if (tmpdenom>1D-30) crossext=crossext+tmpmat(ibas,irealmo)**2/tmpdenom*crossij*2
				end if
				!Cross term between fragment 1 and 2. We assume there is no intersection set between frag1 and frag2
				!Note: any(frag2==jbas) is a subset of .not.any(frag1==jbas), so ovpfrg12 is part of crossext
				if (allocated(frag2).and.any(frag2==jbas)) then
					ovpfrg12(imo)=ovpfrg12(imo)+2*crossij
					if (isel==1) then
						ovpfrg12_1(imo)=ovpfrg12_1(imo)+crossij
					else if (isel==2) then
						tmpdenom=tmpmat(ibas,irealmo)**2+tmpmat(jbas,irealmo)**2
						if (tmpdenom>1D-30) ovpfrg12_1(imo)=ovpfrg12_1(imo)+tmpmat(ibas,irealmo)**2/tmpdenom*crossij*2
					end if
				end if
			end do
		end if
	end do
	if (MOtype(imo)==0) orbtype="AB"
	if (MOtype(imo)==1) orbtype="A "
	if (MOtype(imo)==2) orbtype="B "
	if (isel==1.or.isel==2) write(*,"(i6,2x,a,f12.4,f9.5,f10.3,' %',f10.3,' %',f10.3,' %',f10.3,' %')") &
	imo,orbtype,MOene(imo),MOocc(imo),floc*100,crossint*100,crossext*100,(floc+crossint+crossext)*100
	if (isel==3) write(*,"(i6,5x,a,f16.5,f11.5,f14.5,' %')") imo,orbtype,MOene(imo),MOocc(imo),floc*100
end do

if (isel/=3.and.allocated(frag2)) then !Print cross term between fragment 1 and 2
	write(*,*)
	write(*,*) "Cross term between fragment 1 and 2 and their individual parts:"
	write(*,"(' Orb#  Type   Ene(a.u.)   Occ       Frag.1 part     Frag.2 part        Total')")
	do imo=1,nmo
		if (MOtype(imo)==0) orbtype="AB"
		if (MOtype(imo)==1) orbtype="A "
		if (MOtype(imo)==2) orbtype="B "
		write(*,"(i6,2x,a,f12.4,f9.5,f14.4,' %',f14.4,' %',f14.4,' %')") &
		imo,orbtype,MOene(imo),MOocc(imo),ovpfrg12_1(imo)*100,(ovpfrg12(imo)-ovpfrg12_1(imo))*100,ovpfrg12(imo)*100
	end do
end if
end subroutine




!!---------- Calculate basis function, shell and atom contribution to orbitals
!imethod=1 : Mulliken partition    =2:Stout & Politzer partition   =3: Ros & Schuit (SCPA) partition
subroutine orballcomp_MMPA(imethod)
use defvar
implicit real*8 (a-h,o-z)
real*8 basloc(nbasis),bascross(nbasis),bastot(nbasis)
real*8,pointer :: tmpmat(:,:)
real*8 atmcomp(ncenter)
integer imethod
character orbtype*10,c80tmp*80

do while(.true.)
    write(*,*)
    write(*,*) "Input the orbital index to print composition, e.g. 4"
    if (allocated(CObasb)) write(*,"(a)") " Note: Positive index and negative index correspond to alpha orbital and beta orbital, respectively"
    write(*,*) "Input ""a"" can print basic information of all orbitals, input 0 to return"
    if (wfntype==0.or.wfntype==1.or.wfntype==3) call orblabsel_prompt
    read(*,"(a)") c80tmp
    
	if (c80tmp=="0") then
		exit
	else if (c80tmp=="a") then !Show all orbital information
        call showorbinfo(1,nmo)
	else
		if (index(c80tmp,'h')==0.and.index(c80tmp,'l')==0) then 
			read(c80tmp,*) ishowmo
			if (ishowmo<0.and.allocated(CObasb)) ishowmo=abs(ishowmo)+nbasis
			if (ishowmo<=0.or.ishowmo>nmo) then
				write(*,"(' Error: Orbital index should be in the range of 1 to',i6)") nmo
				cycle
			end if
        else
			call orblabsel(c80tmp,ishowmo)
            if (ishowmo==0) then
				write(*,*) "Error: The orbital label you inputted is wrong! Please double check"
				cycle
            end if
        end if
		if (MOtype(ishowmo)==0) orbtype="Alpha&Beta"
		if (MOtype(ishowmo)==1) orbtype="Alpha     "
		if (MOtype(ishowmo)==2) orbtype="Beta      "
		if (ishowmo<=nbasis) tmpmat=>CObasa
		if (ishowmo>nbasis) tmpmat=>CObasb
		write(*,"(' Orbital:',i6,'  Energy(a.u.):',f14.6,'  Occ:',f10.6,'  Type: ',a)") ishowmo,MOene(ishowmo),MOocc(ishowmo),orbtype
		write(*,"(' Printing threshold of absolute value:  >',f10.5,' %')") compthres
        write(*,*)
		if (imethod==1.or.imethod==2) write(*,"('  Basis Type    Atom    Shell      Local       Cross term        Total')")
		if (imethod==3) write(*,"('  Basis Type    Atom    Shell   Composition')")
		if (ishowmo>nbasis) ishowmo=ishowmo-nbasis !For wfntype==1.or.wfntype==4, change to #beta for CObasb
		if (imethod==3) allsqr=sum(tmpmat(:,ishowmo)**2)
		bascross=0D0
		do ibas=1,nbasis
			basloc(ibas)=tmpmat(ibas,ishowmo)**2
			if (imethod==3) then
				basloc(ibas)=basloc(ibas)/allsqr
			else if (imethod==1.or.imethod==2) then
				do jbas=1,nbasis
					if (jbas==ibas) cycle
					if (imethod==1) then
						bascross(ibas)=bascross(ibas)+tmpmat(ibas,ishowmo)*tmpmat(jbas,ishowmo)*Sbas(ibas,jbas)
					else if (imethod==2) then
						tmp=tmpmat(ibas,ishowmo)**2+tmpmat(jbas,ishowmo)**2
						if (tmp>1D-30) bascross(ibas)=bascross(ibas)+tmpmat(ibas,ishowmo)**2/tmp*2*tmpmat(ibas,ishowmo)*tmpmat(jbas,ishowmo)*Sbas(ibas,jbas)
					end if
				end do
			end if
			bastot(ibas)=basloc(ibas)+bascross(ibas)
			if (abs(bastot(ibas))*100>compthres) then
				if (imethod==1.or.imethod==2) write(*,"(i6,3x,a,i5,a,i5,f13.5,' %',f13.5,' %',f13.5,' %')") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),'('//a(bascen(ibas))%name//')',&
				basshell(ibas),basloc(ibas)*100,bascross(ibas)*100,bastot(ibas)*100
				if (imethod==3) write(*,"(i6,3x,a,i5,a,i5,f13.5,' %')") ibas,GTFtype2name(bastype(ibas)),bascen(ibas),'('//a(bascen(ibas))%name//')',&
				basshell(ibas),bastot(ibas)*100
			end if
		end do
		aboveloc=sum(basloc(:),abs(bastot)*100>compthres)*100
		abovecross=sum(bascross(:),abs(bastot)*100>compthres)*100
		if (imethod==1.or.imethod==2) write(*,"(' Sum up those listed above: ',f13.5,' %',f13.5,' %',f13.5,' %')") aboveloc,abovecross,aboveloc+abovecross
		if (imethod==1.or.imethod==2) write(*,"(' Sum up all basis functions:',f13.5,' %',f13.5,' %',f13.5,' %')") sum(basloc(:))*100,sum(bascross(:))*100,sum(bastot(:))*100
		if (imethod==3) write(*,"(' Sum up those listed above: ',f13.5,' %')") sum(bastot(:),abs(bastot)*100>compthres)*100
		if (imethod==3) write(*,"(' Sum up all basis functions:',f13.5,' %')") sum(bastot(:))*100
		write(*,*)
		write(*,"(' Composition of each shell')")
        s_comp=0;p_comp=0;d_comp=0;f_comp=0;g_comp=0;h_comp=0
		do i=1,nshell
			shellcom=0D0
			do j=1,nbasis
				if (basshell(j)==i) then
					shellcom=shellcom+bastot(j)
					iatm=bascen(j)
				end if
			end do
			if (abs(shellcom)*100>compthres) write(*,"(' Shell',i6,' Type: ',a,'    in atom',i5,'(',a,') :',f12.5,' %')") i,shtype2name(shtype(i)),iatm,a(iatm)%name,shellcom*100
            if (shtype(i)==0) s_comp=s_comp+shellcom
            if (shtype(i)==1) p_comp=p_comp+shellcom
            if (abs(shtype(i))==2) d_comp=d_comp+shellcom
            if (abs(shtype(i))==3) f_comp=f_comp+shellcom
            if (abs(shtype(i))==4) g_comp=g_comp+shellcom
            if (abs(shtype(i))==5) h_comp=h_comp+shellcom
		end do
		write(*,*)
		write(*,*) "Composition of different types of shells (%):"
        write(*,"(' s:',f8.3,'  p:',f8.3,'  d:',f8.3,'  f:',f8.3,'  g:',f8.3,'  h:',f8.3)") s_comp*100,p_comp*100,d_comp*100,f_comp*100,g_comp*100,h_comp*100
		write(*,*)
		write(*,*) "Composition of each atom:"
		do i=1,ncenter
			if (basstart(i)==0) then
				atmcomp(i)=0
			else
	            atmcomp(i)=sum(bastot(basstart(i):basend(i)))*100
            end if
			write(*,"(' Atom',i6,'(',a,') :',f12.5,' %')") i,a(i)%name,atmcomp(i)
		end do
        write(*,"(/,' Orbital delocalization index:',f8.2)") sum(atmcomp**2)/100
		if (allocated(frag1)) then
			write(*,*)
			write(*,"(' Composition of the fragment:',f12.5,' %')") sum(bastot(frag1(:)))*100
		end if
	end if
end do
end subroutine


!!-------- A routine to calculate all atomic contributions to specific orbital, utilized by a few routines (e.g. detecting pi orbitals)
!iorb is the index of the interesting orbital, atmcomp is array with dimension of ncenter to return atomic contribution
!imethod=1: Mulliken    =2: SCPA
subroutine gen_orbatmcomp_MMPA(imethod,iorbin,atmcomp)
use defvar
implicit real*8 (a-h,o-z)
real*8 atmcomp(ncenter),bastot(nbasis)
integer iorbin,imethod
real*8,pointer :: tmpmat(:,:)
if (iorbin<=nbasis) then
	tmpmat=>CObasa
	iorb=iorbin
else if (iorbin>nbasis) then
	tmpmat=>CObasb
	iorb=iorbin-nbasis
end if

atmcomp=0
if (imethod==1) then !Mulliken
	do ibas=1,nbasis
		bascross=0D0
		do jbas=1,nbasis
			if (jbas==ibas) cycle
			bascross=bascross+tmpmat(ibas,iorb)*tmpmat(jbas,iorb)*Sbas(ibas,jbas)
		end do
		bastot(ibas)=tmpmat(ibas,iorb)**2+bascross
	end do
	do iatm=1,ncenter
		if (basstart(iatm)==0) cycle
		atmcomp(iatm)=sum(bastot(basstart(iatm):basend(iatm)))
	end do
else if (imethod==2) then !SCPA
	sumsqr=sum(tmpmat(:,iorb)**2)
	do iatm=1,ncenter
		if (basstart(iatm)==0) cycle
		atmcomp(iatm)=sum(tmpmat(basstart(iatm):basend(iatm),iorb)**2)/sumsqr
	end do
end if
end subroutine



!!!-------- Partition orbital to atomic contributions by space partition methods
!Periodic wavefunction is supported based on evenly distributed grids. Only Hirshfeld is supported currently
!itype: 1=Hirshfeld, 2=Becke, 3=Hirshfeld-I (must have used "call Hirshfeld_I(2)", which generates atomic radial density)
subroutine orbatmcomp_space(itype)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
type(content) gridorg(radpot*sphpot),gridatm(radpot*sphpot)
real*8 resultvec(ncenter)
real*8 allpotx(ncenter,radpot*sphpot),allpoty(ncenter,radpot*sphpot),allpotz(ncenter,radpot*sphpot),allpotw(ncenter,radpot*sphpot)
real*8 tmpdens(radpot*sphpot),selfdens(radpot*sphpot),promol(radpot*sphpot),orbval(nmo)
real*8 atmcomp(ncenter),orbcomp(ncenter,nmo)
integer,allocatable :: idxarr(:)
character c80tmp*80,c2000tmp*2000

if (ifPBC==0) then !Using atomic-center grids to calculate intermediate array
	if (iautointgrid==1) then
		nradpotold=radpot
		nsphpotold=sphpot
		radcutold=radcut
		radpot=30
		sphpot=302
		radcut=15
	end if

	!Generate allpotx/y/z/w
	!allpotx(iatm,j) is x coordinate of the jth point for integrating center iatm
	call gen1cintgrid(gridorg,iradcut)
	do iatm=1,ncenter
		allpotx(iatm,:)=gridorg(:)%x+a(iatm)%x
		allpoty(iatm,:)=gridorg(:)%y+a(iatm)%y
		allpotz(iatm,:)=gridorg(:)%z+a(iatm)%z
	end do

	!allpotw combines Becke multi-center integration weight with Becke/Hirshfeld/Hirshfeld-I weight
	allpotw=0D0
	write(*,"(i6,' quadrature points are used for each atom to compute orbital compositions')") radpot*sphpot
	write(*,"(a)") " Note: You can manually define the number of radial and angular points by &
	&setting ""iautointgrid"" in settings.ini to 0 and setting ""radpot"" and ""sphpot"""
	call walltime(iwalltime1)

	if (itype==1.or.itype==3) then !Hirshfeld or Hirshfeld-I partition
		write(*,*)
		if (itype==1) then
			write(*,*) "Hirshfeld analysis requests atomic densities, please select how to obtain them"
			write(*,*) "1 Use build-in sphericalized atomic densities in free-states (recommended)"
			write(*,"(a)") " 2 Provide wavefunction file of involved elements by yourself or invoke Gaussian to automatically calculate them"
			read(*,*) ihirshmode
		end if
    
		if ((itype==1.and.ihirshmode==1).or.itype==3) then !Hirshfeld or Hirshfeld-I based on interpolation density
			if (itype==1.and.ihirshmode==1) write(*,*) "Generating Hirshfeld atomic weighting functions at all grids..."
			if (itype==3) write(*,*) "Generating Hirshfeld-I atomic weighting functions at all grids..."
			do iatm=1,ncenter
				promol=0D0
				do jatm=1,ncenter
					if (itype==1.and.ihirshmode==1) then !Hirshfeld based on interpolation of built-in atomic radius density
						!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
						do ipt=1+iradcut*sphpot,radpot*sphpot
							tmpdens(ipt)=calcatmdens(jatm,allpotx(iatm,ipt),allpoty(iatm,ipt),allpotz(iatm,ipt),0)
						end do
						!$OMP end parallel do
					else if (itype==3) then !Hirshfeld-I based on refined atomic radial density
						!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
						do ipt=1+iradcut*sphpot,radpot*sphpot
							tmpdens(ipt)=fdens_rad(jatm,allpotx(iatm,ipt),allpoty(iatm,ipt),allpotz(iatm,ipt))
						end do
						!$OMP end parallel do
					end if
					promol=promol+tmpdens
					if (jatm==iatm) selfdens=tmpdens
				end do
				do i=1+iradcut*sphpot,radpot*sphpot !Get Hirshfeld weight of present atom
					if (promol(i)/=0D0) then
						allpotw(iatm,i)=selfdens(i)/promol(i)
					else
						allpotw(iatm,i)=0D0
					end if
				end do
				allpotw(iatm,:)=allpotw(iatm,:)*gridorg(:)%value !Combine Hirshfeld weight with single-center integration weight
            
				call showprog(iatm,ncenter)
			end do
		
		else if (itype==1.and.ihirshmode==2) then !Hirshfeld based on atomic .wfn file
			call setpromol
			write(*,*) "Generating Hirshfeld atomic weighting functions at all grids..."
			do iatm=1,ncenter_org
				promol=0D0
				do jatm=1,ncenter_org
					call dealloall(0)
					call readwfn(custommapname(jatm),1)
					!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
					do ipt=1+iradcut*sphpot,radpot*sphpot
						tmpdens(ipt)=fdens(allpotx(iatm,ipt),allpoty(iatm,ipt),allpotz(iatm,ipt))
					end do
					!$OMP end parallel do
					promol=promol+tmpdens
					if (jatm==iatm) selfdens=tmpdens
				end do
				do i=1+iradcut*sphpot,radpot*sphpot !Get Hirshfeld weight of present atom
					if (promol(i)/=0D0) then
						allpotw(iatm,i)=selfdens(i)/promol(i)
					else
						allpotw(iatm,i)=0D0
					end if
				end do
				allpotw(iatm,:)=allpotw(iatm,:)*gridorg(:)%value !Combine Hirshfeld weight with single-center integration weight
            
				call showprog(iatm,ncenter_org)
			end do
			call dealloall(0)
			call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule) in order to calculate real rho later
		end if

	else if (itype==2) then !Becke partition
		write(*,*) "Generating Becke atomic weighting functions at all grids..."
		do iatm=1,ncenter !Cycle points of each atom
			gridatm%x=gridorg%x+a(iatm)%x
			gridatm%y=gridorg%y+a(iatm)%y
			gridatm%z=gridorg%z+a(iatm)%z
			call gen1cbeckewei(iatm,iradcut,gridatm,allpotw(iatm,:),covr_tianlu,3)
			allpotw(iatm,:)=allpotw(iatm,:)*gridorg%value !Combine Becke weight with single-center integration weight
			call showprog(iatm,ncenter)
		end do
	end if

else !Using evenly distributed grids to calculate composition of all orbitals
	if (itype==1) then
		call gen_orbatmcomp_space(1,orbcomp,1,nmo,1,0)
    else
		write(*,*) "Error: This method is not supported for periodic wavefunctions"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
end if

call walltime(iwalltime2)
write(*,"(' Done! Initialization took up wall clock time',i10,' s')") iwalltime2-iwalltime1

if (allocated(frag1)) deallocate(frag1)
do while(.true.)
    write(*,*)
	write(*,*) "Now input the orbital index to print orbital composition, e.g. 5"
    if (wfntype==1.or.wfntype==4) write(*,"(a)") " Note: If you want to input index of beta orbital, add ""b"" suffix, e.g. 39b"
    if (wfntype==0.or.wfntype==1.or.wfntype==3) call orblabsel_prompt
	write(*,"(a)") " You can also input:"
	if (.not.allocated(frag1)) then
	    write(*,"(a,i6)") " -9: Define fragment, current: undefined"
	else
    	write(*,"(a,i6)") " -9: Redefine fragment, current number of atoms:",size(frag1)
    end if
	write(*,"(a)") "  0: Return"
	write(*,"(a)") " -1: Print basic information of all orbitals"
	write(*,"(a)") " -2: Print atom contribution to a batch of orbitals"
	write(*,"(a)") " -3: Print fragment contribution to a batch of orbitals"
	write(*,"(a)") " -4: Export composition of every atom in every orbital to orbcomp.txt"
	write(*,"(a)") " -5: Print orbital delocalization index (ODI) for a batch of orbitals"
	read(*,"(a)") c80tmp
    	if (index(c80tmp,'h')/=0.or.index(c80tmp,'l')/=0) then 
		call orblabsel(c80tmp,ishowmo)
        if (ishowmo==0) then
			write(*,*) "Error: The orbital label you inputted is wrong! Please double check"
			cycle
        end if
    else if (index(c80tmp,'b')/=0) then !Convert beta orbital index to global orbital index
        itmp=index(c80tmp,'b')
        read(c80tmp(1:itmp-1),*) ishowmo
        do istart=1,nmo !Find starting orbital of beta
	        if (MOtype(istart)==2) exit
        end do
        ishowmo=istart-1+ishowmo
    else
        read(c80tmp,*) ishowmo
    end if
	if (ishowmo>nmo) then
		write(*,"(' Error: Orbital index should between 1 and',i6,/)") nmo
		cycle
	else if (ishowmo==0) then
        if (ifPBC==0.and.iautointgrid==1) then
	        radpot=nradpotold
	        sphpot=nsphpotold
            radcut=radcutold
        end if
        if (allocated(frag1)) deallocate(frag1)
		exit
	else if (ishowmo==-1) then !Show all orbital information
        call showorbinfo(1,nmo)
	else if (ishowmo==-2.or.ishowmo==-3) then !Print atom/fragment contribution to specific range of orbitals
		if ((.not.allocated(frag1)).and.ishowmo==-3) then
			write(*,*) "Error: You must define the fragment by option -9 first!"
			write(*,*)
			cycle
		end if
		if (ishowmo==-2) then
			write(*,*) "Input atom index, e.g. 4"
			read(*,*) iatm
            if (iatm>ncenter) then
                write(*,*) "Error: The atom index exceeded upper limit!"
                cycle
            end if
		end if
        write(*,*) "Input orbital indices, e.g. 3-10,13,18"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp)
        if (allocated(idxarr)) deallocate(idxarr)
        allocate(idxarr(ntmp))
        call str2arr(c2000tmp,ntmp,idxarr)
        if (any(idxarr>nmo)) then
            write(*,*) "Error: One or more orbital indices exceeded upper limit!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            cycle
        end if
        
		write(*,"(' Orb#    Type       Ene(a.u.)     Occ    Composition    Population')")
		pop=0D0
		do idx=1,ntmp
            imo=idxarr(idx)
			if (ishowmo==-2) then !Calculate only one atom
				if (ifPBC==0) then
					tmp=0D0
					!$OMP parallel shared(tmp) private(ipot,value,tmpprivate) num_threads(nthreads)
					tmpprivate=0D0
					!$OMP do schedule(dynamic)
					do ipot=1+iradcut*sphpot,radpot*sphpot
						if (allpotw(iatm,ipot)<1D-8) cycle !May lose 0.001% accuracy
						value=fmo(allpotx(iatm,ipot),allpoty(iatm,ipot),allpotz(iatm,ipot),imo)**2
						tmpprivate=tmpprivate+value*allpotw(iatm,ipot)
					end do
					!$OMP end do
					!$OMP CRITICAL
					tmp=tmp+tmpprivate
					!$OMP end CRITICAL
					!$OMP end parallel
                else
					tmp=orbcomp(iatm,imo)
                end if
			else if (ishowmo==-3) then !Calculate fragment
				if (ifPBC==0) then
					tmp=0D0
					do itmp=1,nfrag1
						iatm=frag1(itmp)
						!$OMP parallel shared(tmp) private(ipot,value,tmpprivate) num_threads(nthreads)
						tmpprivate=0D0
						!$OMP do schedule(dynamic)
						do ipot=1+iradcut*sphpot,radpot*sphpot
							if (allpotw(iatm,ipot)<1D-8) cycle !May lose 0.001% accuracy
							value=fmo(allpotx(iatm,ipot),allpoty(iatm,ipot),allpotz(iatm,ipot),imo)**2
							tmpprivate=tmpprivate+value*allpotw(iatm,ipot)
						end do
						!$OMP end do
						!$OMP CRITICAL
						tmp=tmp+tmpprivate
						!$OMP end CRITICAL
						!$OMP end parallel
					end do
                else
					tmp=sum(orbcomp(frag1(1:nfrag1),imo))
                end if
			end if
			write(*,"(i5,1x,a,f13.4,f9.3,f11.3,' %',f15.6)") imo,orbtypename(MOtype(imo)),MOene(imo),MOocc(imo),tmp*100,MOocc(imo)*tmp
			pop=pop+MOocc(imo)*tmp
		end do
		if (ishowmo==-2) write(*,"(a,f12.6)") " Population of this atom in these orbitals:",pop
		if (ishowmo==-3) write(*,"(a,f12.6)") " Population of this fragment in these orbitals:",pop
		
	else if (ishowmo==-4) then !Export composition of every atom in every orbital to orbcomp.txt in current folder
		if (ifPBC==0) then
			write(*,*) "Calculating, please wait..."
			orbcomp=0
			!$OMP parallel do shared(orbcomp) private(iatm,ipot,orbval) num_threads(nthreads) schedule(dynamic)
			do iatm=1,ncenter
				do ipot=1+iradcut*sphpot,radpot*sphpot
					if (allpotw(iatm,ipot)<1D-9) cycle !May lose 0.001% accuracy
					call orbderv(1,1,nmo,allpotx(iatm,ipot),allpoty(iatm,ipot),allpotz(iatm,ipot),orbval)
					orbcomp(iatm,:)=orbcomp(iatm,:)+orbval(:)**2*allpotw(iatm,ipot)
				end do
			end do
			!$OMP end parallel do
        end if
		open(10,file="orbcomp.txt",status="replace")
		do imo=1,nmo
			write(10,"(' Orbital',i6)") imo
            valnorm=sum(orbcomp(:,imo))
            if (valnorm/=0) then !Some orbitals may have all-zero coefficients, composition is thus all-zero
				do iatm=1,ncenter
					write(10,"(i6,f11.3,' %')") iatm,orbcomp(iatm,imo)/valnorm*100
				end do
            end if
		end do
		close(10)
		write(*,*) "Done! orbcomp.txt has been exported to current folder"
        
	else if (ishowmo==-5) then !Print orbital delocalization index (ODI) for a range of orbitals
        write(*,*) "Input orbital indices, e.g. 3-10,13,18"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp)
        if (allocated(idxarr)) deallocate(idxarr)
        allocate(idxarr(ntmp))
        call str2arr(c2000tmp,ntmp,idxarr)
        if (any(idxarr>nmo)) then
            write(*,*) "Error: One or more orbital indices exceeded upper limit!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            cycle
        end if
        if (ifPBC==0) then
			write(*,*) "Please wait..."
			orbcomp=0
			imax=maxval(idxarr)
			imin=minval(idxarr)
			!$OMP parallel do shared(orbcomp) private(iatm,ipot,orbval) num_threads(nthreads) schedule(dynamic)
			do iatm=1,ncenter
				do ipot=1+iradcut*sphpot,radpot*sphpot
					if (allpotw(iatm,ipot)<1D-8) cycle !May lose 0.001% accuracy
					call orbderv(1,imin,imax,allpotx(iatm,ipot),allpoty(iatm,ipot),allpotz(iatm,ipot),orbval)
					orbcomp(iatm,imin:imax)=orbcomp(iatm,imin:imax)+orbval(imin:imax)**2*allpotw(iatm,ipot)
				end do
			end do
			!$OMP end parallel do
        end if
        if (allocated(frag1)) write(*,*) "ODI of the whole system:"
        do idx=1,ntmp
            imo=idxarr(idx)
            orbcomp(:,imo)=orbcomp(:,imo)/sum(orbcomp(:,imo)) !Normalize
            odi=sum((orbcomp(:,imo)*100)**2)/100
            write(*,"(' Orb:',i5,'  E(a.u.):',f14.6,'  Occ:',f8.4,'  Type: ',a,'  ODI:',f7.2)") &
            imo,MOene(imo),MOocc(imo),orbtypename(MOtype(imo)),odi
        end do
        if (allocated(frag1)) then
            write(*,*)
            write(*,*) "ODI of the fragment you defined:"
            do idx=1,ntmp
                imo=idxarr(idx)
                fragcomp=sum(orbcomp(frag1(:),imo))
                odi=sum((orbcomp(frag1(:),imo)*100/fragcomp)**2)/100
                write(*,"(' Orb:',i5,'  E(a.u.):',f14.6,'  Occ:',f8.4,'  Type: ',a,'  ODI:',f7.2)") &
                imo,MOene(imo),MOocc(imo),orbtypename(MOtype(imo)),odi
            end do
        end if
    
	else if (ishowmo==-9) then !Define fragment
        call definefragment
	    
	else !Print composition for an orbital
		write(*,"(' Orbital:',i5,'  Energy(a.u.):',f14.6,'  Occ:',f10.5,'  Type: ',a)") ishowmo,MOene(ishowmo),MOocc(ishowmo),orbtypename(MOtype(ishowmo))
        if (ifPBC==0) then
			write(*,*) "Please wait..."
			accum=0D0
			do iatm=1,ncenter
				tmp=0D0
				!$OMP parallel shared(tmp) private(ipot,value,tmpprivate) num_threads(nthreads)
				tmpprivate=0
				!$OMP do schedule(dynamic)
				do ipot=1+iradcut*sphpot,radpot*sphpot
					if (allpotw(iatm,ipot)<1D-8) cycle !May lose 0.001% accuracy
					value=fmo(allpotx(iatm,ipot),allpoty(iatm,ipot),allpotz(iatm,ipot),ishowmo)**2
					tmpprivate=tmpprivate+value*allpotw(iatm,ipot)
				end do
				!$OMP end do
				!$OMP CRITICAL
				tmp=tmp+tmpprivate
				!$OMP end CRITICAL
				!$OMP end parallel
				accum=accum+tmp
				resultvec(iatm)=tmp
			end do
			write(*,"(' The sum of contributions before normalization',11x,f12.6,' %',/)") accum*100
        end if
		write(*,*) "Contributions after normalization:"
		do iatm=1,ncenter
			if (ifPBC==0) then
	            atmcomp(iatm)=resultvec(iatm)/accum*100
            else
				atmcomp(iatm)=orbcomp(iatm,ishowmo)*100
            end if
			write(*,"(' Atom',i6,'(',a,') :',f11.3,' %')") iatm,a(iatm)%name,atmcomp(iatm)
		end do
        orbdeloc=sum(atmcomp**2)/100
        write(*,"(/,' Orbital delocalization index:',f8.2)") orbdeloc
		if (allocated(frag1)) then
		    fragcomp=sum(atmcomp(frag1))
            write(*,"(/,' Fragment contribution:',f11.3,' %')") fragcomp
            orbdeloc=sum((atmcomp(frag1(:))/(fragcomp/100))**2)/100
            write(*,"(' Orbital delocalization index of the fragment:',f8.2)") orbdeloc
        end if
	end if
end do
end subroutine





!!!-------- Calculate atomic contributions to specific range of orbitals by Hirshfeld (built-in density) or Becke method
!In fact this is a simplified version of subroutine orbatmcomp_space
!itype=1: Hirshfeld (supports PBC), itype=2: Becke (only isolated systems)
!atmcomp: The array returned, atmcomp(A,i) is contribution of atom A to orbital ibeg+i-1. The second index has length of iend-ibeg+1
!ibeg,iend: The beginning and ending index of the orbitals to be computed, ranging from 1 to nmo
!info=0: silent mode, info=1: print intermediate prompts
!igrid=0: Use lowest acceptable grid (accurate to one decimal place), =1: Use medium quality grid
subroutine gen_orbatmcomp_space(itype,atmcomp,ibeg,iend,info,igrid)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8 atmcomp(ncenter,iend-ibeg+1),atmcomp_tmp(iend-ibeg+1,ncenter),atmcomp_tmpsum(iend-ibeg+1,ncenter),comptmp(iend-ibeg+1)
type(content) gridorg(radpot*sphpot),gridatm(radpot*sphpot)
real*8 allpotx(radpot*sphpot,ncenter),allpoty(radpot*sphpot,ncenter),allpotz(radpot*sphpot,ncenter),allpotw(radpot*sphpot,ncenter)
real*8 tmpdens(radpot*sphpot),selfdens(radpot*sphpot),promol(radpot*sphpot),orbval(nmo)
real*8 atmrho(ncenter),tvec(3)

if (ifPBC/=0.and.itype==2) then
	write(*,*) "Error: Becke method was not implemented for periodic case!"
    write(*,*) "Press ENTER button to exit"
    read(*,*)
    stop
end if

if (itype==1) write(*,*) "Calculating atomic contributions to orbitals by Hirshfeld method..."
if (itype==2) write(*,*) "Calculating atomic contributions to orbitals by Becke method..."
atmcomp=0

if (ifPBC==0) then !Using atomic-center grids
	if (iautointgrid==1) then
		nradpotold=radpot
		nsphpotold=sphpot
		radcutold=radcut
		if (igrid==0) then
			radpot=25
			sphpot=170
		else if (igrid==1) then
			radpot=30
			sphpot=302
		end if
		radcut=15
	end if
	if (info==1) then
		write(*,"(i7,' integration points are used for each atom')") radpot*sphpot
		write(*,"(a)") " Note: You can manually define the number of radial and angular points by &
		&setting ""iautointgrid"" in settings.ini to 0 and setting ""radpot"" and ""sphpot"""
	end if

	call gen1cintgrid(gridorg,iradcut)
	do iatm=1,ncenter
		allpotx(:,iatm)=gridorg(:)%x+a(iatm)%x
		allpoty(:,iatm)=gridorg(:)%y+a(iatm)%y
		allpotz(:,iatm)=gridorg(:)%z+a(iatm)%z
	end do

	!allpotw combines Becke multi-center integration weight with Becke/Hirshfeld weight
	allpotw=0D0
	if (itype==1) then !Hirshfeld based on built-in atomic density
		if (info==1) write(*,*) "Generating Hirshfeld atomic weighting functions at all grids..."
		do iatm=1,ncenter
			promol=0D0
			do jatm=1,ncenter
				!$OMP parallel do shared(tmpdens) private(ipt) num_threads(nthreads)
				do ipt=1+iradcut*sphpot,radpot*sphpot
					tmpdens(ipt)=calcatmdens(jatm,allpotx(ipt,iatm),allpoty(ipt,iatm),allpotz(ipt,iatm),0)
				end do
				!$OMP end parallel do
				promol=promol+tmpdens
				if (jatm==iatm) selfdens=tmpdens
			end do
			do i=1+iradcut*sphpot,radpot*sphpot !Get Hirshfeld weight of present atom
				if (promol(i)/=0D0) then
					allpotw(i,iatm)=selfdens(i)/promol(i)
				else
					allpotw(i,iatm)=0D0
				end if
			end do
			allpotw(:,iatm)=allpotw(:,iatm)*gridorg(:)%value !Combine Hirshfeld weight with single-center integration weight
		end do
	else if (itype==2) then !Becke partition
		if (info==1) write(*,*) "Generating Becke weights..."
		do iatm=1,ncenter !Cycle points of each atom
			gridatm%x=gridorg%x+a(iatm)%x
			gridatm%y=gridorg%y+a(iatm)%y
			gridatm%z=gridorg%z+a(iatm)%z
			call gen1cbeckewei(iatm,iradcut,gridatm,allpotw(:,iatm),covr_tianlu,3)
			allpotw(:,iatm)=allpotw(:,iatm)*gridorg%value !Combine Becke weight with single-center integration weight
		end do
	end if

	if (info==1) write(*,*) "Calculating orbital compositions, please wait..."
	do iatm=1,ncenter
		!$OMP parallel shared(atmcomp) private(ipt,orbval,comptmp) num_threads(nthreads)
		comptmp=0
		!$OMP do schedule(dynamic)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			if (allpotw(ipt,iatm)<1D-9) cycle !May lose 0.001% accuracy
			call orbderv(1,ibeg,iend,allpotx(ipt,iatm),allpoty(ipt,iatm),allpotz(ipt,iatm),orbval)
			comptmp(:)=comptmp(:)+orbval(ibeg:iend)**2*allpotw(ipt,iatm)
		end do
		!$OMP END DO
		!$OMP CRITICAL
		atmcomp(iatm,:)=atmcomp(iatm,:)+comptmp(:)
		!$OMP END CRITICAL
		!$OMP end parallel
		call showprog(iatm,ncenter)
	end do

	if (iautointgrid==1) then
		radpot=nradpotold
		sphpot=nsphpotold
		radcut=radcutold
	end if

else !Using evenly distributed grids for peridic systems
	!call walltime(iwalltime1)
    if (igrid==0) spcgrd=0.35D0
    if (igrid==1) spcgrd=0.25D0
    if (info==1) write(*,"(' Grid spacing of',f6.3,' Bohr is used in integration')") spcgrd
	call setgrid_for_PBC(spcgrd,2) !This is fully adequate for crude estimation
	call calc_dvol(dvol)
	call gen_neigh_GTF !Generate neighbouring GTFs list at reduced grids, for faster evaluation in orbderv
	ifinish=0;ishowprog=1
	ntmp=floor(ny*nz/100D0)
	!$OMP PARALLEL SHARED(atmcomp_tmpsum,ifinish,ishowprog) PRIVATE(atmcomp_tmp,i,j,k,tmpx,tmpy,tmpz,icell,jcell,kcell,tvec,iatm,dist2,atmrho,prorho,orbval) NUM_THREADS(nthreads)
	atmcomp_tmp(:,:)=0
    atmcomp_tmpsum(:,:)=0
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
                    call orbderv_PBC(1,ibeg,iend,tmpx,tmpy,tmpz,orbval)
					do iatm=1,ncenter
						atmcomp_tmp(:,iatm)=atmcomp_tmp(:,iatm)+atmrho(iatm)/prorho*orbval(ibeg:iend)**2
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
    atmcomp_tmpsum(:,:)=atmcomp_tmpsum(:,:)+atmcomp_tmp(:,:)
	!$OMP END CRITICAL
	!$OMP END PARALLEL
	if (ishowprog/=0) call showprog(100,100)
    atmcomp(:,:)=dvol*transpose(atmcomp_tmpsum(:,:))
	!call walltime(iwalltime2)
	!write(*,"(' Calculation of orbital compositions took',i10,' s')") iwalltime2-iwalltime1
end if

!Do normalization
do imo=ibeg,iend
	!write(*,"(' Orbital',i6)") imo
	!do iatm=1,ncenter
	!	write(*,"(i6,f12.6,' %')") iatm,atmcomp(iatm,imo)*100
	!end do
 !   write(*,"(i6,f12.6,/)") imo,sum(atmcomp(:,imo))*100
    itmp=imo-ibeg+1
    allsum=sum(atmcomp(:,itmp))
    if (allsum==0) cycle !When MO is larger than basis function, highest several MO will be filled by blank information to make them equal
    atmcomp(:,itmp)=atmcomp(:,itmp)/allsum
end do

end subroutine





!!------------ Orbital composition analysis by NAO method
subroutine orballcomp_NAO
use defvar
use NAOmod
use util
implicit real*8 (a-h,o-z)
integer :: ioutmode=3
integer :: ispinmode=0 !The spin under study. 0/1/2 = closed shell/alpha/beta
integer,allocatable :: fragidx(:),termtmp(:),tmparr(:)
character :: c80tmp*80,c200*200,c3tmp*3
real*8 :: outcrit=0.5D0
character c2000tmp*2000
real*8,allocatable :: atmcomp(:),shcomp(:)

open(10,file=filename,status="old")

call checkNPA(ifound);if (ifound==0) return
call loadNAOinfo

!Get actual number of MOs
!Gaussian may eliminate some linear dependency basis functions, so MO may be smaller than numNAO. NBsUse always equals to actual number of MO
call loclabel(10,"NBsUse=",ifound)
if (ifound==1) then
	read(10,"(a)") c80tmp
	!For DKH case, G09 may output such as RelInt: Using uncontracted basis, NUniq=   180 NBsUse=   180 0.100E-13, this nbsuse is meaningless, use next nbsuse
	if (index(c80tmp,"RelInt")/=0) call loclabel(10,"NBsUse=",ifound)
	itmp=index(c80tmp,'NBsUse=')
	read(c80tmp(itmp+8:),*) numorb
else !This make this analysis could be independent of Gaussian
    numorb=numNAO
end if
write(*,"(' The number of MOs:',i10)") numorb

!Note that in the current context, numorb=nbasis even for open shell case, because we explicitly distinguish spin
call checkNAOMO(ifound);if (ifound==0) return
call loadNAOMO(numorb)

allocate(fragidx(numNAO),termtmp(numNAO))
ifragend=0 !ifragend is actual ending position (number of elements) of fragidx array
if (iopshNAO==1) ispinmode=1

!Interface
NAOcompmaincyc: do while(.true.)
do while(.true.)
	write(*,*)
    write(*,*) "------- Orbital composition analysis based on natural atomic orbitals -------"
	write(*,*) "-10 Return"
	write(*,*) "-1 Define fragment for option 1"
	write(*,*) " 0 Show composition of an orbital"
	write(*,*) " 1 Show fragment contribution to a batch of orbitals"
	if (ioutmode==0) write(*,"(a)") "  2 Select output mode for option 0, current: All terms"
	if (ioutmode==1) write(*,"(a)") "  2 Select output mode for option 0, current: Non-Rydberg terms"
	if (ioutmode==2) write(*,"(a,f6.2,'%')") "  2 Select output mode for option 0, current: All terms with contribution >",outcrit
	if (ioutmode==3) write(*,"(a,f6.2,'%')") "  2 Select output mode for option 0, current: Non-Rydberg terms with contribution >",outcrit
	if (ispinmode==1) write(*,*) " 3 Switch spin type, current: Alpha"
	if (ispinmode==2) write(*,*) " 3 Switch spin type, current: Beta"
	read(*,*) isel

	if (isel==-10) then
		close(10)
		return
	else if (isel==-1) then
20		write(*,*) "Commands and examples:"
		write(*,*) "q: Save setting and exit"
		write(*,*) "help: Print these help content again"
		write(*,*) "clean: Clean current fragment"
		write(*,*) "list: List NAOs in current fragment"
		write(*,*) "all: Print all NAOs of current system"
		write(*,*) "addall: Add all NAOs to fragment"
		write(*,*) "a 2,5-8,12: Add all NAOs in atoms 2,5,6,7,8,12 to fragment"
		write(*,*) "b 2,5-8,12: Add NAOs 2,5,6,7,8,12 to fragment"
		write(*,*) "da 5,7,11-13: Delete all NAOs in atoms 5,7,11,12,13 from fragment"
		write(*,*) "db 5,7,11-13: Delete NAOs 5,7,11,12,13 from fragment"
		do while(.true.)
            write(*,*)
            write(*,*) "Please input command. Input ""help"" can show help. Input ""q"" can save and exit"
			read(*,"(a)") c200
			if (index(c200,"help")/=0) then
				goto 20
			else if (index(c200,"addall")/=0) then
				forall(i=1:numNAO) fragidx(i)=i
				ifragend=numNAO
				write(*,*) "Done!"
			else if (c200(1:1)=="q") then
				write(*,*)
				write(*,*) "Fragment is saved, indices of NAOs in current fragment:"
				if (ifragend==0) then
					write(*,*) "None"
				else
					write(*,"(12i6)") fragidx(1:ifragend)
				end if
				exit
			else if (index(c200,"clean")/=0) then
				fragidx=0
				ifragend=0
				write(*,*) "Done!"
			else if (index(c200,"all")/=0) then
				write(*,*) "Note: The ones with asterisks are those presented in current fragment"
                write(*,*) "      NAO#   Atom&Index    Type   Set&Shell    Occupancy   Energy (a.u.)"
                do iNAO=1,numNAO
                    c3tmp=" "
                    if (any(fragidx(1:ifragend)==iNAO)) c3tmp=" * "
                    write(*,"(a,i7,3x,a5,i5,5x,a7,1x,a3,'(',a,')',f12.5,f14.5)") &
                    c3tmp,iNAO,NAOcenname(iNAO),NAOcen(iNAO),NAOtype(iNAO),NAOset(iNAO,ispinmode),NAOshname(iNAO),NAOocc(iNAO,ispinmode),NAOene(iNAO,ispinmode)
                end do
			else if (index(c200,"list")/=0) then
				write(*,*) "NAO indices in current fragment:"
				if (ifragend==0) then
					write(*,*) "None"
				else
					write(*,"(12i6)") fragidx(1:ifragend)
                    write(*,*)
                    write(*,*) "Detailed information:"
                    write(*,*) "      NAO#   Atom&Index    Type   Set&Shell    Occupancy   Energy (a.u.)"
                    do idx=1,ifragend
                        iNAO=fragidx(idx)
                        write(*,"(a,i7,3x,a5,i5,5x,a7,1x,a3,'(',a,')',f12.5,f14.5)") &
                        "   ",iNAO,NAOcenname(iNAO),NAOcen(iNAO),NAOtype(iNAO),NAOset(iNAO,ispinmode),NAOshname(iNAO),NAOocc(iNAO,ispinmode),NAOene(iNAO,ispinmode)
                    end do
				end if
				
			else if (c200(1:2)=="a ".or.c200(1:2)=="s ".or.c200(1:2)=="b ".or.c200(1:2)=="db".or.c200(1:2)=="da") then
				call str2arr(c200(3:),nterm,termtmp)
				if (c200(1:2)=="a ".or.c200(1:2)=="da") then !Check sanity of the input
					if (any(termtmp(1:nterm)<=0).or.any(termtmp(1:nterm)>ncenter_NAO)) then
						write(*,*) "Error: Atom index exceeded valid range! Ignoring..."
						cycle
					end if
				else if (c200(1:2)=="b ".or.c200(1:2)=="db") then
					if (any(termtmp(1:nterm)<=0).or.any(termtmp(1:nterm)>numNAO)) then
						write(*,*) "Error: NAO index exceeded valid range! Ignoring..."
						cycle
					end if
				end if
				!Modify fragidx according to the selection
				if (c200(1:2)=="a ") then
					do iterm=1,nterm
						iatm=termtmp(iterm)
						do ibas=NAOinit(iatm),NAOend(iatm)
							if ( any(fragidx(1:ifragend)==ibas) ) cycle
							ifragend=ifragend+1
							fragidx(ifragend)=ibas
						end do
					end do
				else if (c200(1:2)=="b ") then
					do ibas=1,nterm
						if ( any(fragidx(1:ifragend)==termtmp(ibas)) ) cycle
						ifragend=ifragend+1
						fragidx(ifragend)=termtmp(ibas)
					end do
				else if (c200(1:2)=="db") then
					do itmp=1,nterm
						do iscan=1,ifragend
							if (fragidx(iscan)==termtmp(itmp)) then
								fragidx(iscan:ifragend-1)=fragidx(iscan+1:ifragend)
								ifragend=ifragend-1
								exit
							end if
							if (iscan==ifragend) write(*,"('Note: NAO with index of',i7,' does not present in current fragment')") termtmp(itmp)
						end do
					end do
				else if (c200(1:2)=="da") then
					do itmp=1,nterm
						iatm=termtmp(itmp)
						do iNAOincen=NAOinit(iatm),NAOend(iatm) !Cycle NAO index in itmp atom
							do iscan=1,ifragend !Scan fragidx array to find out iNAOincen
								if (fragidx(iscan)==iNAOincen) then
									fragidx(iscan:ifragend-1)=fragidx(iscan+1:ifragend)
									ifragend=ifragend-1
									exit
								end if
								if (iscan==ifragend) write(*,"('Note: NAO with index of',i7,' does not present in current fragment')") iNAOincen
							end do
						end do
					end do
				end if
				write(*,*) "Done!"
			else
				write(*,*) "Error: Cannot recognize the command you inputted"
			end if
			
		end do !End of fragment definition module
		
	else if (isel==0) then
		exit
		
	else if (isel==1) then
		if (ifragend==0) then
			write(*,*) "Error: You have not defined fragment or the fragment is empty!"
		else
			write(*,*) "Input orbital indices, for which the composition will be printed"
            write(*,*) "For example 3,5-8,15-20"
            read(*,"(a)") c2000tmp
            if (allocated(tmparr)) deallocate(tmparr)
            call str2arr(c2000tmp,ntmp)
            allocate(tmparr(ntmp))
            call str2arr(c2000tmp,ntmp,tmparr)
			if (any(tmparr<=0).or.any(tmparr>numorb)) then
				write(*,*) "Error: The range of orbital indices is invalid!"
				cycle
			else
				exit
			end if
		end if
		
	else if (isel==2) then
		write(*,*) "0 Show all NAOs/shells"
		write(*,*) "1 Show non-Rydberg NAOs/shells"
		write(*,*) "2 Show NAOs/shells whose contributions are larger than specific criterion"
		write(*,"(a)") " 3 Show non-Rydberg NAOs/shells whose contributions are larger than specific criterion"
		read(*,*) ioutmode
		if (ioutmode==2.or.ioutmode==3) then
			write(*,*) "Input printing criterion in %, e.g. 0.5"
			read(*,*) outcrit
		end if
		
	else if (isel==3) then
		if (ispinmode==1) then
			ispinmode=2
		else if (ispinmode==2) then
			ispinmode=1
		end if
	end if
end do

!Start analysis
if (isel==0) then !Analyze one orbital
    allocate(atmcomp(ncenter_NAO),shcomp(numNAOsh))
	do while(.true.)
		write(*,*) "Input the index of the orbital to be analyzed, e.g. 5"
        write(*,*) "Input 0 can return"
		read(*,*) iorboutcomp
		if (iorboutcomp==0) then
            if (allocated(atmcomp)) deallocate(atmcomp,shcomp)
            exit
		else if (iorboutcomp<=0.or.iorboutcomp>numorb) then
			write(*,"(a,i7)") " Error: The orbital index should between  1 and",numorb
			cycle
		end if
		
		if (ioutmode==2) write(*,"(a,f6.2,a)") " Note: All NAOs/shells whose contributions <=",outcrit," % will not be printed"
		if (ioutmode==3) write(*,"(a,f6.2,a)") " Note: All Rydberg NAOs/shells or contributions <=",outcrit," % will not be printed"
		if (ispinmode==1) write(*,"(a,i6)") " Below are composition of alpha orbitals",iorboutcomp
		if (ispinmode==2) write(*,"(a,i6)") " Below are composition of beta orbitals",iorboutcomp
        write(*,*)
		write(*,*) "   NAO#   Center   Label      Type    Composition"
        shcomp=0D0
        atmcomp=0D0
        Corcomp=0D0
        Valcomp=0D0
        Rydcomp=0D0
		do iNAO=1,numNAO
			tmpcomp=NAOMO(iNAO,iorboutcomp)**2*100
            if (ispinmode==2) tmpcomp=NAOMOb(iNAO,iorboutcomp)**2*100 !Use beta MO instead
            if (NAOset(iNAO,ispinmode)=="Cor") Corcomp=Corcomp+tmpcomp
            if (NAOset(iNAO,ispinmode)=="Val") Valcomp=Valcomp+tmpcomp
            if (NAOset(iNAO,ispinmode)=="Ryd") Rydcomp=Rydcomp+tmpcomp
            shcomp(bassh_NAO(iNAO))=shcomp(bassh_NAO(iNAO))+tmpcomp
            atmcomp(NAOcen(iNAO))=atmcomp(NAOcen(iNAO))+tmpcomp
			if ((ioutmode==1.or.ioutmode==3).and.NAOset(iNAO,ispinmode)=="Ryd") cycle !Skip showing Rydberg
			if ((ioutmode==2.or.ioutmode==3).and.tmpcomp<outcrit) cycle !Skip showing too small terms
			write(*,"(i8,i5,'(',a,')',4x,a,2x,a,'(',a,')',f10.3,' %' )") &
            iNAO,NAOcen(iNAO),NAOcenname(iNAO),NAOtype(iNAO),NAOset(iNAO,ispinmode),NAOshname(iNAO),tmpcomp
		end do
		write(*,*)
		write(*,*) "Condensed NAO terms to shells:"
        s_comp=0;p_comp=0;d_comp=0;f_comp=0;g_comp=0;h_comp=0
        do ish=1,numNAOsh
            icen=shcen_NAO(ish)
            if (index(shname_NAO(ish),'s')/=0) s_comp=s_comp+shcomp(ish)
            if (index(shname_NAO(ish),'p')/=0) p_comp=p_comp+shcomp(ish)
            if (index(shname_NAO(ish),'d')/=0) d_comp=d_comp+shcomp(ish)
            if (index(shname_NAO(ish),'f')/=0) f_comp=f_comp+shcomp(ish)
            if (index(shname_NAO(ish),'g')/=0) g_comp=g_comp+shcomp(ish)
            if (index(shname_NAO(ish),'h')/=0) h_comp=h_comp+shcomp(ish)
            if ((ioutmode==1.or.ioutmode==3).and.shset_NAO(ish,ispinmode)=="Ryd") cycle !Skip showing Rydberg
			if ((ioutmode==2.or.ioutmode==3).and.shcomp(ish)<outcrit) cycle !Skip showing too small terms
            write(*,"('   Atom:',i6,'(',a,')  Shell:',i6,'(',a,1x,a,')',f10.3,' %' )") &
            icen,atmname_NAO(shcen_NAO(ish)),ish,shname_NAO(ish),shset_NAO(ish,ispinmode),shcomp(ish)
        end do
        write(*,*)
		write(*,*) "Composition of different types of shells (%):"
        write(*,"(' s:',f8.3,'  p:',f8.3,'  d:',f8.3,'  f:',f8.3,'  g:',f8.3,'  h:',f8.3)") s_comp,p_comp,d_comp,f_comp,g_comp,h_comp
		write(*,*)
		write(*,*) "Condensed NAO terms to atoms:"
		write(*,*) "  Center   Composition"
		do icen=1,ncenter_NAO
            if ((ioutmode==2.or.ioutmode==3).and.atmcomp(icen)<outcrit) cycle
			write(*,"(i6,'(',a,')',f10.3,' %' )") icen,atmname_NAO(icen),atmcomp(icen)
		end do
		write(*,*)
		write(*,"(' Core composition:   ',f10.3,' %')") Corcomp
		write(*,"(' Valence composition:',f10.3,' %')") Valcomp
		write(*,"(' Rydberg composition:',f10.3,' %')") Rydcomp
		write(*,*)
        write(*,"(' Orbital delocalization index:',f8.2)") sum(atmcomp(:)**2)/100
		write(*,*)
	end do
    
else if (isel==1) then !Show fragment contribution in a range of orbitals
	if (ispinmode==1) write(*,*) "Below are composition of alpha orbitals"
	if (ispinmode==2) write(*,*) "Below are composition of beta orbitals"
	write(*,*) " Orb.#       Core      Valence     Rydberg      Total"
	do imoidx=1,ntmp
        imo=tmparr(imoidx)
		sumcompcor=0D0
		sumcompval=0D0
		sumcompryd=0D0
		do itmp=1,ifragend
			iNAO=fragidx(itmp)
			tmpcomp=NAOMO(iNAO,imo)**2*100
            if (ispinmode==2) tmpcomp=NAOMOb(iNAO,imo)**2*100
			if (NAOset(iNAO,ispinmode)=="Cor") sumcompcor=sumcompcor+tmpcomp
			if (NAOset(iNAO,ispinmode)=="Val") sumcompval=sumcompval+tmpcomp
			if (NAOset(iNAO,ispinmode)=="Ryd") sumcompryd=sumcompryd+tmpcomp
		end do
		sumcomptot=sumcompcor+sumcompval+sumcompryd
		write(*,"(i6,3x,4(f10.3,' %'))") imo,sumcompcor,sumcompval,sumcompryd,sumcomptot
	end do
end if

end do NAOcompmaincyc

end subroutine




!!---------- Localized orbital bonding analysis (LOBA) and modified LOBA (mLOBA)
!Input file should contains localized MO
subroutine LOBA
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 oxdstat(ncenter),atmcomp(ncenter,nmo)
integer,allocatable :: fragLOBA(:)
character c2000tmp*2000,c80tmp*80

if (wfntype==0) then
	call gen_orbatmcomp_space(1,atmcomp,1,nint(nelec/2),0,0)
else
	call gen_orbatmcomp_space(1,atmcomp,1,nmo,0,0)
end if
write(*,*)
write(*,"(a)") " Hint: If you want to make Multiwfn output attribution of LMOs, &
&you can set ""outmedinfo"" in settings.ini to 1 before booting up Multiwfn"
write(*,*)

nfragLOBA=0
do while(.true.)
	write(*,*) "Input percentage threshold to determine oxidation states (50 is commonly used)"
    write(*,"(a)") " Input ""m"" will determine oxidation states by assigning electron(s) in each orbital to the atom of maximal contribution, which corresponds to modified LOBA (mLOBA)"
	if (allocated(fragLOBA)) then
        write(*,"(a,i6,a)") " Input -1 can redefine the fragment, current:",nfragLOBA," atoms"
    else
        write(*,*) "Input -1 can define fragment, current: Undefined"
    end if
	write(*,*) "Input 0 can exit"
	read(*,"(a)") c80tmp
	if (c80tmp=="0") then
		return
	else if (c80tmp=="-1") then
		write(*,*) "Input indices of the atoms constituting the fragment"
		write(*,*) "e.g. 1,3-5,9"
		read(*,"(a)") c2000tmp
		if (allocated(fragLOBA)) deallocate(fragLOBA)
		call str2arr(c2000tmp,nfragLOBA)
		allocate(fragLOBA(nfragLOBA))
		call str2arr(c2000tmp,nfragLOBA,fragLOBA)
		write(*,*) "Done!"
		write(*,*)
		cycle
	else
        oxdstat(:)=a(:)%charge
		if (index(c80tmp,'m')==0) then
			read(c80tmp,*) thres
			thres=thres/100
			do iatm=1,ncenter
				do imo=1,nmo
					if (MOocc(imo)==0) cycle
					if (atmcomp(iatm,imo)>thres) then
						oxdstat(iatm)=oxdstat(iatm)-MOocc(imo)
 						if (outmedinfo==1) write(*,"(f8.4' electrons in LMO',i6,' are attributed to atom',i6)") MOocc(imo),imo,iatm
					end if
				end do
				write(*,"(' Oxidation state of atom',i4,'(',a,') :',i3)") iatm,a(iatm)%name,nint(oxdstat(iatm))
			end do
        else
			do imo=1,nmo
				if (MOocc(imo)==0) cycle
				iatm=maxloc(atmcomp(:,imo),1)
				oxdstat(iatm)=oxdstat(iatm)-MOocc(imo)
				if (outmedinfo==1) write(*,"(f8.4' electrons in LMO',i6,' are attributed to atom',i6)") MOocc(imo),imo,iatm
			end do
			do iatm=1,ncenter
				write(*,"(' Oxidation state of atom',i4,'(',a,') :',i3)") iatm,a(iatm)%name,nint(oxdstat(iatm))
			end do
        end if
        
		write(*,"(' The sum of oxidation states:',i6)") sum(nint(oxdstat))
        
		if (allocated(fragLOBA)) then
			if (index(c80tmp,'m')==0) then
				oxdfrag=0
				do iatmidx=1,nfragLOBA
					iatm=fragLOBA(iatmidx)
					oxdfrag=oxdfrag+a(iatm)%charge
				end do
				do imo=1,nmo
					compval=sum(atmcomp(fragLOBA(:),imo))
					if (compval>thres) then
						oxdfrag=oxdfrag-MOocc(imo)
						if (outmedinfo==1) write(*,"(f8.4' electrons in LMO',i6,' are attributed to the defined fragment')") MOocc(imo),imo
					end if
				end do
				write(*,"(' Oxidation state of the fragment:',i4)") nint(oxdfrag)
            else
				write(*,"(' Oxidation state of the fragment:',i4)") sum(nint(oxdstat(fragLOBA(:))))
            end if
		end if
		write(*,*)
	end if
end do

end subroutine






!!------- Generate all atomic composition of all orbitals by SCPA method, or load them from orbcomp.txt in current folder
!Currently not employed by any function
subroutine gen_allorbatmcomp_SCPA(atmcomp)
use defvar
implicit real*8 (a-h,o-z)
real*8 atmcomp(ncenter,nmo)

atmcomp=0
inquire(file="orbcomp.txt",exist=alive)
if (alive) then !Load atomic contribution from orbcomp.txt, which may be outputted by option -4 of Hirshfeld/Becke composition analysis
	write(*,"(a)") " orbcomp.txt was found in current folder, now load atomic contribution to all orbitals from this file..."
	open(10,file="orbcomp.txt",status="old")
	do imo=1,nmo
		read(10,*)
		do iatm=1,ncenter
			read(10,*,iostat=ierror) inouse,atmcomp(iatm,imo)
			if (ierror/=0) then
				write(*,"(a)") " Error: Problem was encountered while loading orbcomp.txt in current folder! Please carefully check this file!"
				write(*,*) "Press ENTER button to exit program"
				read(*,*)
				stop
			end if
		end do
	end do
	close(10)
	atmcomp=atmcomp/100
else
	write(*,"(a)") " orbcomp.txt was not found in current folder, therefore calculate atomic contributions to all orbitals by SCPA method..."
	do imo=1,nmo
		if (MOtype(imo)==0.or.MOtype(imo)==1) then !Closed-shell or alpha part of open-shell
			sumsqr=sum(CObasa(:,imo)**2)
			do iatm=1,ncenter
				if (basstart(iatm)==0) cycle
				atmcomp(iatm,imo)=sum(CObasa(basstart(iatm):basend(iatm),imo)**2)/sumsqr
			end do
		else !Beta part of open-shell
			iimo=imo-nbasis
			sumsqr=sum(CObasb(:,iimo)**2)
			do iatm=1,ncenter
				if (basstart(iatm)==0) cycle
				atmcomp(iatm,imo)=sum(CObasb(basstart(iatm):basend(iatm),iimo)**2)/sumsqr
			end do
		end if
	end do
end if
end subroutine



!!------- Generate all basis function composition in all orbitals by SCPA method
subroutine gen_allorbbascomp_SCPA(bascomp)
use defvar
implicit real*8 (a-h,o-z)
real*8 bascomp(nbasis,nmo)
do imo=1,nmo
	if (MOtype(imo)==0.or.MOtype(imo)==1) then !Closed-shell or alpha part of open-shell
		sumsqr=sum(CObasa(:,imo)**2)
		do ibas=1,nbasis
			bascomp(ibas,imo)=CObasa(ibas,imo)**2/sumsqr
		end do
	else !Beta part of open-shell
		iimo=imo-nbasis
		sumsqr=sum(CObasb(:,iimo)**2)
		do ibas=1,nbasis
			bascomp(ibas,imo)=CObasb(ibas,iimo)**2/sumsqr
		end do
	end if
end do
end subroutine