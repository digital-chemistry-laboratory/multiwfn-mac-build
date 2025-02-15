!-------- Main interface of various bond order analyses
subroutine bondorder_main
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000
do while(.true.)
	write(*,*)
	write(*,*) "           ================ Bond order analysis ==============="
	if (allocated(b)) then
		if (allocated(frag1)) then
			write(*,*) "-1 Redefine fragment 1 and 2 for options 1,3,4,7,8,10"
		else
			write(*,*) "-1 Define fragment 1 and 2 for options 1,3,4,7,8,10 (to be defined)"
		end if
	end if
	write(*,*) "0 Return"
	write(*,*) "1 Mayer bond order analysis"
	write(*,*) "2 Multicenter bond order analysis"
	write(*,*) "-2 Multicenter bond order analysis in NAO basis"
! 	write(*,*) "-3 Multicenter bond order analysis in Lowdin orthogonalized basis" !Can be used, but not very meaningful, so not shown
	write(*,*) "3 Wiberg bond order analysis in Lowdin orthogonalized basis"
	write(*,*) "4 Mulliken bond order (Mulliken overlap population) analysis"
	write(*,*) "5 Decompose Mulliken bond order between two atoms to orbital contributions"
	write(*,*) "6 Orbital occupancy-perturbed Mayer bond order"
	write(*,*) "7 Fuzzy bond order analysis (FBO)"
	write(*,*) "8 Laplacian bond order (LBO)"
	write(*,*) "9 Decompose Wiberg bond order in NAO basis as atomic orbital pair contribution"
    write(*,*) "10 Intrinsic bond strength index (IBSI)"
    write(*,*) "11 AV1245 index (approximate multicenter bond order for large rings) and AVmin"
    write(*,*) "20 Bond order density (BOD) and natural adaptive orbital (NAdO) analyses"
	read(*,*) ibondana
    
    if (ifPBC>0.and.(ibondana==-2.or.ibondana==8.or.ibondana==9.or.ibondana==10)) then
		write(*,*) "Error: This function does not support periodic systems currently!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        cycle
    end if
    
	if (.not.allocated(CObasa).and.((ibondana>=1.and.ibondana<=6).or.ibondana==20)) then
		write(*,"(a)") " Error: The input file you used does not contain basis function information! Please check Section 2.5 of the manual for explanation"
		write(*,*) "Press ENTER button to return"
		read(*,*)
		return
	else if (.not.allocated(b).and.(ibondana==7.or.ibondana==8)) then
		write(*,"(a)") " Error: The input file you used does not contain GTF information! Please check Section 2.5 of the manual for explanation"
		write(*,*) "Press ENTER button to return"
		read(*,*)
		return
	end if
	
	if (ibondana==-1) then
		!Define frag1, the size just accomodates content
		if (allocated(frag1)) then
			write(*,*) "Atoms in current fragment 1:"
			write(*,"(13i6)") frag1
			write(*,"(a)") " Input 0 to keep unchanged, or redefine fragment 1, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment 1"
		else
			write(*,"(a)") " Input atomic indices to define fragment 1. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment 1"
		end if
		read(*,"(a)") c2000tmp
		if (c2000tmp(1:1)/='0') then
			if (allocated(frag1)) deallocate(frag1)
			call str2arr(c2000tmp,nfragatm)
			allocate(frag1(nfragatm))
			call str2arr(c2000tmp,nfragatm,frag1)
		end if
		!Define frag2, the size just accomodates content
		if (allocated(frag2)) then
			write(*,*) "Atoms in current fragment 2:"
			write(*,"(13i6)") frag2
			write(*,"(a)") " Input 0 to keep unchanged, or redefine fragment 2, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute fragment 2"
		else
			write(*,"(a)") " Input atomic indices to define fragment 2. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment 2"
		end if
		read(*,"(a)") c2000tmp
		if (c2000tmp(1:1)/='0') then
			if (allocated(frag2)) deallocate(frag2)
			call str2arr(c2000tmp,nfragatm)
			allocate(frag2(nfragatm))
			call str2arr(c2000tmp,nfragatm,frag2)
		end if
		if (any(frag1>ncenter).or.any(frag2>ncenter)) then
			write(*,*) "Error: Some atomic indices exceeded valid range! Please define again"
			write(*,*)
			deallocate(frag1,frag2)
			cycle
		end if
		write(*,*) "Setting is saved"
		write(*,*) "Now the atoms in fragment 1 are"
		write(*,"(13i6)") frag1
		write(*,*) "Now the atoms in fragment 2 are"
		write(*,"(13i6)") frag2
		if (any(frag1<=0).or.any(frag1>ncenter).or.any(frag2<=0).or.any(frag2>ncenter)) write(*,*) "Warning: Indices of some atoms exceed valid range! Please redefine fragment"
		do i=1,size(frag1)
			if (any(frag2==frag1(i))) then
				write(*,"(a)") "Warning: Indices of some atoms are duplicated in the two fragments! Please redefine them"
				exit
			end if
		end do
		write(*,*)
		
	else if (ibondana==0) then
		if (allocated(frag1)) deallocate(frag1)
		if (allocated(frag2)) deallocate(frag2)
		exit
	else if (ibondana==1) then
        call ask_Sbas_PBC
		write(*,*) "Calculating, please wait..."
		call mayerbndord
	else if (ibondana==2) then
        call ask_Sbas_PBC
		call multicenter(2)
	else if (ibondana==-2) then
		call multicenterNAO
	else if (ibondana==3.or.ibondana==-3) then
        call ask_Sbas_PBC
	    !In symmortho the density matrix, CObasa/b and Sbas will change, so backup them
	    if (allocated(CObasb)) then !Open-shell
            allocate(CObasa_org(nbasis,nmo/2),CObasb_org(nbasis,nmo/2),Sbas_org(nbasis,nbasis))
	        CObasb_org=CObasb
	    else
			allocate(CObasa_org(nbasis,nmo),Sbas_org(nbasis,nbasis)) 
	    end if
	    CObasa_org=CObasa
	    Sbas_org=Sbas
	    write(*,*) "Performing Lowdin orthogonalization..."
 		call symmortho(0)
 		if (ibondana==3) then
			write(*,*) "Calculating Wiberg bond order..."
			call mayerbndord
		else if (ibondana==-3) then
			call multicenter(-3)
		end if
        write(*,*) "Regenerating original density matrix..."
        write(*,*)
        CObasa=CObasa_org
        Sbas=Sbas_org
        deallocate(CObasa_org,Sbas_org)
        if (allocated(CObasb_org)) then
            CObasb=CObasb_org
            deallocate(CObasb_org)
        end if
        call genP
	else if (ibondana==4) then
        call ask_Sbas_PBC
		call mullikenbndord
	else if (ibondana==5) then
        call ask_Sbas_PBC
		call decompMullikenBO
	else if (ibondana==6) then
        call ask_Sbas_PBC
		call OrbPertMayer
	else if (ibondana==7) then
		call fuzzyana(1)
	else if (ibondana==8) then
		write(*,"(/,a)") " Citation of Laplacian bond order (LBO):" 
		write(*,"(a,/)") " Tian Lu and Feiwu Chen, &
		&Bond Order Analysis Based on the Laplacian of Electron Density in Fuzzy Overlap Space, J. Phys. Chem. A, 117, 3100-3108 (2013)"
		call fuzzyana(2)
	else if (ibondana==9) then
		call decompWibergNAO
	else if (ibondana==10) then
        call IBSI
    else if (ibondana==11) then
        call AV1245
    else if (ibondana==20) then
        call BOD
	end if
end do
end subroutine


!! ----------------- Mayer/Generalized Wiberg 2-c bond order analysis
! If Lowdin orthogonalization has been performed, that is carry out Wiberg bond order analysis in Lowdin orthogonalized basis
! Note: For closed-shell, two methods give the same result. For open-shell, Mayer bond order for all electrons is the sum of
! alpha and beta bond order, while GWBO directly use total density matrix to generate
! total bond order, the "Mayer bond order" in Gaussian is actually GWBO!
subroutine Mayerbndord
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 bndmata(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmattot(ncenter,ncenter),PSmattot(nbasis,nbasis)
character selectyn

call calcMayerbndord(bndmattot,bndmata,bndmatb)

write(*,"(' Bond orders with absolute value >=',f10.6)") bndordthres
itmp=0
if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	do i=1,ncenter
		do j=i+1,ncenter
			if (abs(bndmata(i,j)+bndmatb(i,j))>=bndordthres) then
				itmp=itmp+1
				write(*,"(' #',i5,':',i5,a,i5,a,' Alpha: ',f10.6,' Beta:',f10.6,' Total:',f10.6)") &
				itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmata(i,j),bndmatb(i,j),bndmata(i,j)+bndmatb(i,j)
			end if
		end do
	end do
	write(*,*)    
    write(*,"(a)") " Note: The ""Total"" bond orders shown above are more meaningful than the below ones. If you are not familiar &
    &related theory, you can simply ignore below output"
    write(*,*)
	write(*,"(' Bond order from mixed alpha&beta density matrix >=',f10.6)") bndordthres
end if
itmp=0
do i=1,ncenter
	do j=i+1,ncenter
		if (abs(bndmattot(i,j))>=bndordthres) then
			itmp=itmp+1
			write(*,"(' #',i5,':',5x,i5,a,i5,a,f14.8)") itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmattot(i,j)
		end if
	end do
end do

write(*,*)
write(*,*) "Total valences and free valences defined by Mayer:"
PSmattot=matmul_blas(Ptot,Sbas,nbasis,nbasis)
do i=1,ncenter
    if (basstart(i)==0) cycle
	accum=0D0
	accum2=0D0
	do ii=basstart(i),basend(i)
		accum=accum+2*PSmattot(ii,ii)
		do jj=basstart(i),basend(i)
			accum2=accum2+PSmattot(ii,jj)*PSmattot(jj,ii)
		end do
	end do
	freeval=accum-accum2-(bndmata(i,i)+bndmatb(i,i))
	if (wfntype==0.or.wfntype==3) freeval=0D0
	write(*,"(' Atom',i6,'(',a,') :',2f14.8)") i,a(i)%name,accum-accum2,freeval
end do

!Between fragment
if (allocated(frag1)) then
	bndordfraga=0
	bndordfragb=0
	bndordfragtot=0
	do i=1,size(frag1)
		do j=1,size(frag2)
			bndordfraga=bndordfraga+bndmata(frag1(i),frag2(j))
			bndordfragb=bndordfragb+bndmatb(frag1(i),frag2(j))
			bndordfragtot=bndordfragtot+bndmattot(frag1(i),frag2(j))
		end do
	end do
	write(*,*)
	if (wfntype==1.or.wfntype==2.or.wfntype==4) then
		write(*,"(' The bond order between fragment 1 and 2:')")
		write(*,"(' Alpha:',f10.6,' Beta:',f10.6,' Total:',f10.6,' Mixed Alpha&Beta:',f10.6)") bndordfraga,bndordfragb,bndordfraga+bndordfragb,bndordfragtot
	else
		write(*,"(' The bond order between fragment 1 and 2:',f12.6)") bndordfragtot
	end if
end if

write(*,*)
write(*,*) "If outputting bond order matrix to bndmat.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	open(10,file="bndmat.txt",status="replace")
	write(10,*) "Note: The diagonal elements are the sum of corresponding row elements"
	if (wfntype==0.or.wfntype==3) then
		call showmatgau(bndmattot,"Bond order matrix",0,"f14.8",10)
	else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
		call showmatgau(bndmata,"Bond order matrix for alpha electrons",0,"f14.8",10)
		call showmatgau(bndmatb,"Bond order matrix for beta electrons",0,"f14.8",10)
		call showmatgau(bndmata+bndmatb,"Bond order matrix for all electrons",0,"f14.8",10)
		call showmatgau(bndmattot,"Bond order matrix from mixed density",0,"f14.8",10)
	end if
	close(10)
	write(*,*) "Result have been outputted to bndmat.txt in current folder"
	write(*,*)
end if
end subroutine



!!------ Return bond order matrix of Mayer bond order
!For restricted closed-shell wavefunction, only bndmattot is returned
!For unrestricted wavefunction, bndmattot corresponds to generalized Wiberg bond order, bndmata and bndmatb correspond to alpha and beta Mayer bond order, respectively
subroutine calcMayerbndord(bndmattot,bndmata,bndmatb)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 bndmata(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmattot(ncenter,ncenter),&
PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis),PSmattot(nbasis,nbasis)

bndmata=0D0
bndmatb=0D0
bndmattot=0D0

!Calculate total bond order for restricted closed-shell wavefunction (for open-shell calculate GWBO, namely using density matrix of Palpha+Pbeta)
PSmattot=matmul_blas(Ptot,Sbas,nbasis,nbasis)
do i=1,ncenter
	if (basstart(i)==0) cycle
	do j=i+1,ncenter
		if (basstart(j)==0) cycle
		accum=0D0
		do ii=basstart(i),basend(i)
			do jj=basstart(j),basend(j)
				accum=accum+PSmattot(ii,jj)*PSmattot(jj,ii)
			end do
		end do
		bndmattot(i,j)=accum
	end do
end do
bndmattot=bndmattot+transpose(bndmattot) !Because we only filled one triangular region, copy it to another
do i=1,ncenter
	bndmattot(i,i)=sum(bndmattot(:,i))
end do

!Unrestricted wavefunction
if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	PSmata=matmul_blas(Palpha,Sbas,nbasis,nbasis)
	PSmatb=matmul_blas(Pbeta,Sbas,nbasis,nbasis)
	bndmata=0
	bndmatb=0
	do i=1,ncenter
		if (basstart(i)==0) cycle
		do j=i+1,ncenter
			if (basstart(j)==0) cycle
			accuma=0D0
			accumb=0D0
			do ii=basstart(i),basend(i)
				do jj=basstart(j),basend(j)
					accuma=accuma+PSmata(ii,jj)*PSmata(jj,ii)
					accumb=accumb+PSmatb(ii,jj)*PSmatb(jj,ii)
				end do
			end do
			bndmata(i,j)=accuma
			bndmatb(i,j)=accumb
		end do
	end do
	bndmata=2*(bndmata+transpose(bndmata))
	bndmatb=2*(bndmatb+transpose(bndmatb))
	do i=1,ncenter
		bndmata(i,i)=sum(bndmata(:,i))
		bndmatb(i,i)=sum(bndmatb(:,i))
	end do
end if

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!--------- Calculate Mulliken bond order
subroutine Mullikenbndord
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis),bndmattot(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmata(ncenter,ncenter)
character selectyn
bndmattot=0D0
if (wfntype==0.or.wfntype==3) then
	PSmata=Sbas*Ptot !Condensed to basis function matrix
	do i=1,ncenter !Contract PSmata to Condensed to "Condensed to atoms" 
		if (basstart(i)==0) cycle
		do j=i+1,ncenter
			if (basstart(j)==0) cycle
			bndmattot(i,j)=sum(PSmata(basstart(i):basend(i),basstart(j):basend(j)))
		end do
	end do
	bndmattot=2*(bndmattot+transpose(bndmattot))
	forall (i=1:ncenter) bndmattot(i,i)=sum(bndmattot(i,:))
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	bndmata=0D0
	bndmatb=0D0
	PSmata=Palpha*Sbas
	PSmatb=Pbeta*Sbas
	do i=1,ncenter
		if (basstart(i)==0) cycle
		do j=i+1,ncenter
			if (basstart(j)==0) cycle
			bndmata(i,j)=sum(PSmata(basstart(i):basend(i),basstart(j):basend(j)))
			bndmatb(i,j)=sum(PSmatb(basstart(i):basend(i),basstart(j):basend(j)))
		end do
	end do
	bndmata=2*(bndmata+transpose(bndmata))
	bndmatb=2*(bndmatb+transpose(bndmatb))
	forall (i=1:ncenter) bndmata(i,i)=sum(bndmata(i,:))
	forall (i=1:ncenter) bndmatb(i,i)=sum(bndmatb(i,:))
	bndmattot=bndmata+bndmatb
end if

write(*,"(' Bond orders with absolute value >=',f10.6)") bndordthres
itmp=0
do i=1,ncenter
	do j=i+1,ncenter
		if (wfntype==0.or.wfntype==3) then
			if (abs(bndmattot(i,j))>=bndordthres) then
				itmp=itmp+1
				write(*,"(' #',i5,':',5x,i5,a,i5,a,f14.8)") itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmattot(i,j)
			end if
		else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
			if (abs(bndmata(i,j)+bndmatb(i,j))>=bndordthres) then
				itmp=itmp+1
				write(*,"(' #',i5,':',i5,a,i5,a,' Alpha: ',f10.6,' Beta:',f10.6,' Total:',f10.6)") &
				itmp,i,'('//a(i)%name//')',j,'('//a(j)%name//')',bndmata(i,j),bndmatb(i,j),bndmattot(i,j)
			end if
		end if
	end do
end do
write(*,*)

!Between fragment
if (allocated(frag1)) then
	bndordfraga=0
	bndordfragb=0
	bndordfragtot=0
	do i=1,size(frag1)
		do j=1,size(frag2)
			bndordfraga=bndordfraga+bndmata(frag1(i),frag2(j))
			bndordfragb=bndordfragb+bndmatb(frag1(i),frag2(j))
			bndordfragtot=bndordfragtot+bndmattot(frag1(i),frag2(j))
		end do
	end do
	if (wfntype==1.or.wfntype==2.or.wfntype==4) then
		write(*,"(' The Mulliken bond order between fragment 1 and 2:')")
		write(*,"(' Alpha:',f12.6,' Beta:',f12.6,' Total:',f12.6)") bndordfraga,bndordfragb,bndordfragtot
	else if (wfntype==0.or.wfntype==3) then
		write(*,"(' The Mulliken bond order between fragment 1 and 2:',f12.6)") bndordfragtot
	end if
	write(*,*)
end if

write(*,*) "If outputting bond order matrix to bndmat.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	open(10,file="bndmat.txt",status="replace")
	write(10,*) "Note:The diagonal elements are the sum of corresponding row elements"
	if (wfntype==0.or.wfntype==3) then
		call showmatgau(bndmattot,"Mulliken bond order matrix",0,"f14.8",10)
	else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
		call showmatgau(bndmata,"Mulliken bond order matrix for alpha electrons",0,"f14.8",10)
		call showmatgau(bndmatb,"Mulliken bond order matrix for beta electrons",0,"f14.8",10)
		call showmatgau(bndmattot,"Mulliken bond order matrix all electrons",0,"f14.8",10)
	end if
	close(10)
	write(*,*) "Result have been outputted to bndmat.txt in current folder"
	write(*,*)
end if
end subroutine



!!--------- Decompose Mulliken bond order to MO contribution
subroutine decompMullikenBO
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,pointer :: ptmat(:,:)

do while(.true.)
	write(*,*) "Input index of two atom (e.g. 3,5)"
	write(*,*) "Note: Input 0,0 can return to upper level menu"
	read(*,*) ind1,ind2

	if (ind1==0.and.ind2==0) exit
	bndorda=0D0
	bndordb=0D0
	do itime=1,2
		if (itime==1) ptmat=>CObasa
		if (itime==2) ptmat=>CObasb
		if (itime==1.and.(wfntype==1.or.wfntype==4)) write(*,*) "Alpha orbitals:"
		if (itime==2.and.(wfntype==1.or.wfntype==4)) write(*,*) "Beta orbitals:"
		do imo=1,nbasis
			if (itime==1) irealmo=imo
			if (itime==2) irealmo=imo+nbasis
			if (MOocc(irealmo)==0D0) cycle
			accum=0D0
			do i=basstart(ind1),basend(ind1)
				do j=basstart(ind2),basend(ind2)
					accum=accum+MOocc(irealmo)*ptmat(i,imo)*ptmat(j,imo)*Sbas(i,j)
				end do
			end do
			if (itime==1) bndorda=bndorda+accum*2
			if (itime==2) bndordb=bndordb+accum*2
			write(*,"(' Orbital',i6,' Occ:',f10.6,' Energy:',f12.6,' contributes',f14.8)") imo,MOocc(irealmo),MOene(irealmo),accum*2
		end do
		if (wfntype==0.or.wfntype==2.or.wfntype==3) then
			write(*,"(' Total Mulliken bond order:',f14.8,/)") bndorda
			exit
		else if (wfntype==1.or.wfntype==4) then
			if (itime==1) write(*,"(' Mulliken bond order of all alpha electrons:',f14.8,/)") bndorda
			if (itime==2) write(*,"(' Mulliken bond order of all beta electrons:',f14.8,/)") bndordb
			if (itime==2) write(*,"(' Total Mulliken bond order:',f14.8,/)") bndorda+bndordb
		end if
	end do
end do
end subroutine



!!----- Orbital occupancy-perturbed Mayer bond order (Decompose Mayer bond-order between two atoms to orbital contributions)
!--- J. Chem. Theory Comput. 2012, 8, 908, 914
!For simplicity, this routine only calculate Mayer bond for alpha and beta and then sum them up, don't concern mixed alpha+beta cases
subroutine OrbPertMayer
use defvar
use util
implicit real*8 (a-h,o-z)
character orbtypechar*2
real*8 bndmata(ncenter,ncenter),bndmatb(ncenter,ncenter),bndmattot(ncenter,ncenter)
real*8,allocatable :: PSmattot(:,:),Ptottmp(:,:)
real*8,allocatable :: PSmata(:,:),PSmatb(:,:),Palphatmp(:,:),Pbetatmp(:,:)
bndmata=0D0
bndmatb=0D0
do while(.true.)
write(*,*) "Input indices of two atoms, e.g. 3,5"
	read(*,*) iatm,jatm
	if (iatm>=1.and.iatm<=ncenter.and.jatm>=1.and.jatm<=ncenter.and.iatm/=jatm) exit
	write(*,*) "Error: Invalid input, please input again"
end do
if (wfntype==0.or.wfntype==3) then !Close shell
	allocate(PSmattot(nbasis,nbasis),Ptottmp(nbasis,nbasis))
	sumupvar=0D0
	do imo=0,nmo !Cycle all MOs
		Ptottmp=Ptot !Do not use Ptot to make troubles, because Ptot is a global array
		if (imo/=0) then !Calculate perturbed density. At the first time (imo=1), we don't pertube density matrix to yield original Mayer bond order
			if (MOocc(imo)<=1D-10) cycle
			do ibas=1,nbasis
				do jbas=1,nbasis
					Ptottmp(ibas,jbas)=Ptottmp(ibas,jbas)-MOocc(imo)*CObasa(ibas,imo)*CObasa(jbas,imo)
				end do
			end do
		end if
		PSmattot=matmul_blas(Ptottmp,Sbas,nbasis,nbasis) !Calculate Mayer bond order based on Ptottmp
		bndordtot=0D0
		do ii=basstart(iatm),basend(iatm)
			do jj=basstart(jatm),basend(jatm)
				bndordtot=bndordtot+PSmattot(ii,jj)*PSmattot(jj,ii)
			end do
		end do
		if (imo==0) then
			beforepert=bndordtot
			write(*,"(' Mayer bond order before orbital occupancy-perturbation:',f12.6)") beforepert
			write(*,*)
			write(*,"(' Mayer bond order after orbital occupancy-perturbation:')")
			write(*,*) "Orbital     Occ      Energy    Bond order   Variance"
		else
			bndordvar=bndordtot-beforepert
			write(*,"(i6,f12.5,f11.5,2f12.6)") imo,MOocc(imo),MOene(imo),bndordtot,bndordvar
			sumupvar=sumupvar+bndordvar
		end if
	end do
	write(*,"(' Summing up occupancy perturbation from all orbitals:',f10.5)") sumupvar
	
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open shell
	sumupvar=0D0
	allocate(PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis),Palphatmp(nbasis,nbasis),Pbetatmp(nbasis,nbasis))
	do imo=0,nmo
		Palphatmp=Palpha
		Pbetatmp=Pbeta
		if (imo/=0) then !The first time, we calculate actual Mayer bond order
			if (MOocc(imo)<=1D-10) cycle
			if (wfntype==1.or.wfntype==4) then
				if (imo<=nbasis) then !Alpha orbitals
					do ibas=1,nbasis
						do jbas=1,nbasis
							Palphatmp(ibas,jbas)=Palphatmp(ibas,jbas)-MOocc(imo)*CObasa(ibas,imo)*CObasa(jbas,imo)
						end do
					end do
				else !Beta orbitals, between nbasis+1 and nmo
					do ibas=1,nbasis
						do jbas=1,nbasis
							Pbetatmp(ibas,jbas)=Pbetatmp(ibas,jbas)-MOocc(imo)*CObasb(ibas,imo-nbasis)*CObasb(jbas,imo-nbasis)
						end do
					end do
				end if
			else if (wfntype==2) then !ROHF
				if (MOtype(imo)==0) then !Doubly occupied orbitals
					do ibas=1,nbasis
						do jbas=1,nbasis
							Palphatmp(ibas,jbas)=Palphatmp(ibas,jbas)-1D0*CObasa(ibas,imo)*CObasa(jbas,imo)
							Pbetatmp(ibas,jbas)=Pbetatmp(ibas,jbas)-1D0*CObasa(ibas,imo)*CObasa(jbas,imo) !For ROHF, CObasb==CObasa, and hence CObasb is not allocated
						end do
					end do
				else if (MOtype(imo)==1) then !Alpha orbitals
					do ibas=1,nbasis
						do jbas=1,nbasis
							Palphatmp(ibas,jbas)=Palphatmp(ibas,jbas)-1D0*CObasa(ibas,imo)*CObasa(jbas,imo)
						end do
					end do				
				end if
			end if
		end if
		PSmata=matmul_blas(Palphatmp,Sbas,nbasis,nbasis)
		PSmatb=matmul_blas(Pbetatmp,Sbas,nbasis,nbasis)
		bndorda=0D0
		bndordb=0D0
		do ii=basstart(iatm),basend(iatm)
			do jj=basstart(jatm),basend(jatm)
				bndorda=bndorda+PSmata(ii,jj)*PSmata(jj,ii)
				bndordb=bndordb+PSmatb(ii,jj)*PSmatb(jj,ii)
			end do
		end do
		bndorda=bndorda*2
		bndordb=bndordb*2
		if (imo==0) then
			beforepert=bndorda+bndordb
			write(*,"(' Mayer bond order before orbital occupancy-perturbation:')") 
			write(*,"(' Alpha:',f12.6,'  Beta:',f12.6,'  Total:',f12.6)") bndorda,bndordb,bndorda+bndordb
			write(*,*)
			write(*,"(' Mayer bond order after orbital occupancy-perturbation:')")
			write(*,*) "Orbital     Occ      Energy  Type     Alpha      Beta     Total      Variance"
		else
			bndordvar=bndorda+bndordb-beforepert
			if (MOtype(imo)==0) orbtypechar="AB"
			if (MOtype(imo)==1) orbtypechar="A "
			if (MOtype(imo)==2) orbtypechar="B "		
			write(*,"(i6,f12.6,f11.5,2x,a,2x,3f10.5,3x,f10.5)") imo,MOocc(imo),MOene(imo),orbtypechar,bndorda,bndordb,bndorda+bndordb,bndordvar
			sumupvar=sumupvar+bndordvar
		end if
	end do
	write(*,"(' Summing up occupancy perturbation from all orbitals:',f10.5)") sumupvar
end if
write(*,*)
end subroutine





!------ Decompose Wiberg bond order as NAO pair and NAO shell pair contributions
!NBO output file with DMNAO keyword should be used as input file
subroutine decompWibergNAO
use defvar
use NAOmod
use util
implicit real*8 (a-h,o-z)
character(len=3) icenshname(100),jcenshname(100) !Record all shell type names in centers i and j
character c80tmp*80,c2000tmp*2000
real*8,allocatable :: shcontri(:,:),shcontrib(:,:) !alpha/total part, beta part
integer NAOfrag1(ncenter),NAOfrag2(ncenter)
integer nNAOfrag1,nNAOfrag2
!Arrays for unique shells. Each unique shell corresponds to a unique kind of NAO shell in present system
character(len=3) uniqsh_name(100) !Name of unique shells
integer numuniqsh !Number of unique shells
real*8 uniqsh_val(100,100),uniqsh_valb(100,100) !uniqsh_val(i,j) is contribution to bond order due to interaction of unique shell i in fragment 1 and unique shell j in fragment 2 

!Load NAO and DMNAO information
open(10,file=filename,status="old")
call checkNPA(ifound);if (ifound==0) return
call loadNAOinfo
call checkDMNAO(ifound);if (ifound==0) return
call loadDMNAO
close(10)

!Construct shell list with unique shell name used to record interfragment shell interactions
numuniqsh=0
do ish=1,numNAOsh
    if ( all(uniqsh_name(1:numuniqsh)/=shname_NAO(ish)) ) then
        numuniqsh=numuniqsh+1
        uniqsh_name(numuniqsh)=shname_NAO(ish)
    end if
end do
nNAOfrag1=0
nNAOfrag2=0

write(*,"(a)") " Note: The threshold for printing contribution is controlled by ""bndordthres"" in settings.ini"
do while(.true.)
	write(*,*)
	write(*,*) "Input two atom indices to decompose their bond order, e.g. 3,4"
    write(*,*) "To decompose bond order between two fragments, input -1"
	write(*,*) "Input 0 can exit"
    
	read(*,"(a)") c80tmp
	if (c80tmp(1:1)=='0') then
        return
    else if (c80tmp=="-1") then
        write(*,*) "Input indices of the atoms in fragment 1, e.g. 4,8,9-12,18"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nNAOfrag1,NAOfrag1)
        write(*,*) "Input indices of the atoms in fragment 2, e.g. 1,3,13-16"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nNAOfrag2,NAOfrag2)
    else
	    read(c80tmp,*) iatm,jatm
    end if
    
    if (c80tmp/="-1") then !Interatom analysis
	    !Construct name list of all shells for iatm and jatm
	    numicensh=1
	    icenshname(1)=NAOshname(NAOinit(iatm))
	    do ibas=NAOinit(iatm)+1,NAOend(iatm)
		    if (all(icenshname(1:numicensh)/=NAOshname(ibas))) then
			    numicensh=numicensh+1
			    icenshname(numicensh)=NAOshname(ibas)
		    end if
	    end do
	    numjcensh=1
	    jcenshname(1)=NAOshname(NAOinit(jatm))
	    do jbas=NAOinit(jatm)+1,NAOend(jatm)
		    if (all(jcenshname(1:numjcensh)/=NAOshname(jbas))) then
			    numjcensh=numjcensh+1
			    jcenshname(numjcensh)=NAOshname(jbas)
		    end if
	    end do
	    allocate(shcontri(numicensh,numjcensh),shcontrib(numicensh,numjcensh))
	    shcontri=0
        shcontrib=0
	    bndord=0
	    bndordb=0
    
	    !Calculate Wiberg bond order and output worthnoting components
	    write(*,*) "Contribution from NAO pairs that larger than printing threshold:"
	    if (iopshNAO==0) write(*,*) " Contri.  NAO   Center   NAO type             NAO   Center   NAO type"
	    if (iopshNAO==1) write(*,*) "Spin   Contri.  NAO   Center   NAO type          NAO   Center   NAO type"
	    do iNAO=NAOinit(iatm),NAOend(iatm)
		    do ish=1,numicensh !Find the belonging shell index within this atom for iNAO
			    if (NAOshname(iNAO)==icenshname(ish)) exit
		    end do
		    do jNAO=NAOinit(jatm),NAOend(jatm)
			    do jsh=1,numjcensh
				    if (NAOshname(jNAO)==jcenshname(jsh)) exit
			    end do
                if (iopshNAO==0) then !Closed shell
			        contri=DMNAO(iNAO,jNAO)**2
			        if (contri>bndordthres) write(*,"(f8.4,1x,i5,i5,'(',a,')  ',a,'(',a,') ',a,'--- ',i5,i5,'(',a,')  ',a,'(',a,') ',a)") contri,&
			        iNAO,NAOcen(iNAO),NAOcenname(iNAO),NAOset(iNAO,0),NAOshname(iNAO),NAOtype(iNAO),&
			        jNAO,NAOcen(jNAO),NAOcenname(jNAO),NAOset(jNAO,0),NAOshname(jNAO),NAOtype(jNAO)
			        shcontri(ish,jsh)=shcontri(ish,jsh)+contri
			        bndord=bndord+contri
                else !Open shell
                    contri1=2*DMNAOa(iNAO,jNAO)**2
			        if (contri1>bndordthres) write(*,"(' Alpha',f8.4,1x,i5,i5,'(',a,') ',a,'(',a,') ',a,'--',i5,i5,'(',a,') ',a,'(',a,') ',a)") contri1,&
			        iNAO,NAOcen(iNAO),NAOcenname(iNAO),NAOset(iNAO,1),NAOshname(iNAO),NAOtype(iNAO),&
			        jNAO,NAOcen(jNAO),NAOcenname(jNAO),NAOset(jNAO,1),NAOshname(jNAO),NAOtype(jNAO)
                    contri2=2*DMNAOb(iNAO,jNAO)**2
			        if (contri2>bndordthres) write(*,"(' Beta ',f8.4,1x,i5,i5,'(',a,') ',a,'(',a,') ',a,'--',i5,i5,'(',a,') ',a,'(',a,') ',a)") contri2,&
			        iNAO,NAOcen(iNAO),NAOcenname(iNAO),NAOset(iNAO,2),NAOshname(iNAO),NAOtype(iNAO),&
			        jNAO,NAOcen(jNAO),NAOcenname(jNAO),NAOset(jNAO,2),NAOshname(jNAO),NAOtype(jNAO)
			        shcontri(ish,jsh)=shcontri(ish,jsh)+contri1
			        shcontrib(ish,jsh)=shcontrib(ish,jsh)+contri2
			        bndord=bndord+contri1
			        bndordb=bndordb+contri2
                end if
		    end do
	    end do
	    write(*,*)
        write(*,*) "Contribution from NAO shell pairs that larger than printing threshold:"
	    if (iopshNAO==0) write(*,*) " Contri. Shell  Center   Type        Shell  Center   Type"
	    if (iopshNAO==1) write(*,*) "Spin    Contri. Shell  Center   Type        Shell  Center   Type"
	    do ish=1,numicensh
		    do jsh=1,numjcensh
	            if (iopshNAO==0) then
			        if (shcontri(ish,jsh)>bndordthres) write(*,"(f8.4,1x,i5,i5,'(',a,')    ',a,'   --- ',i5,i5,'(',a,')    ',a)") shcontri(ish,jsh),&
			        ish,NAOcen(NAOinit(iatm)),NAOcenname(NAOinit(iatm)),icenshname(ish),&
			        jsh,NAOcen(NAOinit(jatm)),NAOcenname(NAOinit(jatm)),jcenshname(jsh)
                else
			        if (shcontri(ish,jsh)>bndordthres) write(*,"(' Alpha',f9.4,1x,i5,i5,'(',a,')    ',a,'   --- ',i5,i5,'(',a,')    ',a)") shcontri(ish,jsh),&
			        ish,NAOcen(NAOinit(iatm)),NAOcenname(NAOinit(iatm)),icenshname(ish),&
			        jsh,NAOcen(NAOinit(jatm)),NAOcenname(NAOinit(jatm)),jcenshname(jsh)
			        if (shcontrib(ish,jsh)>bndordthres) write(*,"(' Beta ',f9.4,1x,i5,i5,'(',a,')    ',a,'   --- ',i5,i5,'(',a,')    ',a)") shcontrib(ish,jsh),&
			        ish,NAOcen(NAOinit(iatm)),NAOcenname(NAOinit(iatm)),icenshname(ish),&
			        jsh,NAOcen(NAOinit(jatm)),NAOcenname(NAOinit(jatm)),jcenshname(jsh)
                end if
		    end do
	    end do
    
        !Total result    
        if (iopshNAO==0) then
	        write(*,"(/,a,f8.4)") " Total Wiberg bond order:",bndord
        else
	        write(*,"(/,a,f8.4)") " Total alpha Wiberg bond order:",bndord
	        write(*,"(a,f8.4)") " Total beta Wiberg bond order: ",bndordb
	        write(*,"(a,f8.4)") " Sum of alpha and beta Wiberg bond orders:",bndord+bndordb
        end if
	    deallocate(shcontri,shcontrib)
    
    else !Interfragment analysis
        if (nNAOfrag1>0.and.nNAOfrag2>0) then
            write(*,*)
            write(*,*) "Interfragment bond order analysis:"
            uniqsh_val=0
            uniqsh_valb=0
            do idx=1,nNAOfrag1 !Cycle fragment 1
                iatm=NAOfrag1(idx)
                do jdx=1,nNAOfrag2 !Cycle fragment 1
                    jatm=NAOfrag2(jdx)
                    do iNAO=NAOinit(iatm),NAOend(iatm)
			            do ish=1,numuniqsh
				            if (NAOshname(iNAO)==uniqsh_name(ish)) exit
			            end do
		                do jNAO=NAOinit(jatm),NAOend(jatm)
			                do jsh=1,numuniqsh
				                if (NAOshname(jNAO)==uniqsh_name(jsh)) exit
			                end do
                            if (iopshNAO==0) then !Closed shell
			                    contri=DMNAO(iNAO,jNAO)**2
			                    uniqsh_val(ish,jsh)=uniqsh_val(ish,jsh)+contri
                            else !Open shell
                                contri1=2*DMNAOa(iNAO,jNAO)**2
                                contri2=2*DMNAOb(iNAO,jNAO)**2
			                    uniqsh_val(ish,jsh)=uniqsh_val(ish,jsh)+contri1
			                    uniqsh_valb(ish,jsh)=uniqsh_valb(ish,jsh)+contri2
                            end if
		                end do
                    end do
                end do
            end do
	        if (iopshNAO==0) write(*,*) " Contribution   Fragment 1   Fragment 2"
	        if (iopshNAO==1) write(*,*) "  Spin   Contribution   Fragment 1   Fragment 2"
            do ish=1,numuniqsh
                do jsh=1,numuniqsh
                    if (iopshNAO==0) then !Closed shell
			            if (uniqsh_val(ish,jsh)>bndordthres) write(*,"(f12.5,8x,a,10x,a)") &
                        uniqsh_val(ish,jsh),uniqsh_name(ish),uniqsh_name(jsh)
                    else !Open shell
			            if (uniqsh_val(ish,jsh)>bndordthres) write(*,"('   Alpha',f12.5,8x,a,10x,a)") &
                        uniqsh_val(ish,jsh),uniqsh_name(ish),uniqsh_name(jsh)
			            if (uniqsh_valb(ish,jsh)>bndordthres) write(*,"('   Beta ',f12.5,8x,a,10x,a)") &
                        uniqsh_valb(ish,jsh),uniqsh_name(ish),uniqsh_name(jsh)
                    end if
                end do
            end do
            if (iopshNAO==0) then
	            write(*,"(/,a,f8.4)") " Interfragment Wiberg bond order:",sum(uniqsh_val)
            else
	            write(*,"(/,a,f8.4)") " Interfragment alpha Wiberg bond order: ",sum(uniqsh_val)
	            write(*,"(a,f8.4)") " Interfragment beta Wiberg bond order:  ",sum(uniqsh_valb)
	            write(*,"(a,f8.4)") " Interfragment total Wiberg bond orders:",sum(uniqsh_val)+sum(uniqsh_valb)
            end if
        end if
    end if
end do
end subroutine



!!-----------------------------------------------------------------
!!------------ Intrinsic bond strength index (IBSI) ---------------
!!-----------------------------------------------------------------
subroutine IBSI
use defvar
use util
implicit real*8 (a-h,o-z)
integer atmlist(ncenter),iIGMtype
real*8 atmpairdg(ncenter,ncenter),IBSImat(ncenter,ncenter),IBSIfrag
real*8 :: refval_IGM=0.558245D0,refval_IGMH=0.566653D0,refval_mIGM=0.439922D0 !Calculated under perfect grid for H2 at 0.74144 Angstrom, namely IBSI value obtained with reference value of 1.0
real*8 :: distprintthres=3.5D0
character c2000tmp*2000

natmlist=ncenter
forall(i=1:ncenter) atmlist(i)=i

write(*,*)
if (allocated(b)) then
    iIGMtype=2
    refval=refval_IGMH
    write(*,"(a)") " Note: The current reference value corresponds to H2 in &
    &experimental structure (0.74144 Ang) with density generated at B3LYP/6-311G** level"
else
    iIGMtype=1
    refval=refval_IGM
    write(*,"(a)") " Note: The current reference value corresponds to H2 in experimental structure (0.74144 Ang)"
end if

do while(.true.)
    write(*,*)
    write(*,*) "           ---------- Intrinsic bond strength index (IBSI) ----------"
    write(*,*) "0 Return"
    write(*,*) "1 Start calculation"
    if (iIGMtype==1) write(*,*) "2 Set type of IGM, current: IGM based on promolecular approximation"
    if (iIGMtype==2) write(*,"(a)") " 2 Set type of IGM, current: IGM based on Hirshfeld partition of electron density (IGMH)"
    if (iIGMtype==-1) write(*,*) "2 Set type of IGM, current: Modified IGM (mIGM)"
    if (natmlist==ncenter) then
        write(*,*) "3 Input the range of the atoms to be taken into account, current: all"
    else
        write(*,"(a,i5,' atoms')") " 3 Input the range of the atoms to be taken into account, current:",natmlist
    end if
    write(*,"(a,f10.6)") " 4 Set H2 reference value, current:",refval
    write(*,"(a,f6.3,' Angstrom')") " 5 Set distance threshold for printing result, current: <",distprintthres
    read(*,*) isel
    if (isel==0) then
        return
    else if (isel==2) then
		write(*,*) " 1 IGM based on promolecular approximation (original IGM)"
		write(*,*) " 2 IGM based on Hirshfeld partition of electron density (IGMH)"
		write(*,*) "-1 Modified IGM (mIGM)"
        read(*,*) iIGMtype
        if (iIGMtype==2) then
            if (.not.allocated(b)) then
                write(*,"(a)") " Error: Your input file does not contain wavefunction information, &
                &therefore only the IGM based on promolecular approximation can be used"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                cycle
            end if
            refval=refval_IGMH
            write(*,"(a)") " Note: The current reference value corresponds to H2 in &
            &experimental structure (0.74144 Ang) with density generated at B3LYP/6-311G** level"
        else if (iIGMtype==1) then
            refval=refval_IGM
            write(*,"(a)") " Note: The current reference value corresponds to H2 in experimental structure (0.74144 Ang)"
        else if (iIGMtype==-1) then
            refval=refval_mIGM
            write(*,"(a)") " Note: The current reference value corresponds to H2 in experimental structure (0.74144 Ang)"
        end if
    else if (isel==3) then
        write(*,*) "Input index of the atoms to be considered in the calculation, e.g. 2,3,7-10"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,natmlist,atmlist)
    else if (isel==4) then
        write(*,*) "Input reference value, e.g. 0.324"
        read(*,*) refval
    else if (isel==5) then
        write(*,*) "Input distance threshold for printing in Angstrom, e.g. 3.0"
        write(*,"(a)") " Note: If distance between two atoms is larger than this value, then the corresponding data will not be printed"
        read(*,*) distprintthres
        
    else if (isel==1) then
        call calcatmpairdg(iIGMtype,natmlist,atmlist,natmlist,atmlist,atmpairdg(1:natmlist,1:natmlist))
        write(*,*)
        write(*,"(a)") " Note: ""Dist"" is distance between the two atoms in Angstrom, Int(dg_pair) is the integral &
        &in the numerator of the IBSI formule (atomic pair delta-g index)"
        write(*,*)
        do idx=1,natmlist
            iatm=atmlist(idx)
            do jdx=idx+1,natmlist
                jatm=atmlist(jdx)
                dist=atomdist(iatm,jatm,1)*b2a
                if (dist>distprintthres) cycle
                IBSImat(idx,jdx)=atmpairdg(idx,jdx)/dist**2/refval
                write(*,"(i5,'(',a,')',i5,'(',a,')  Dist:',f8.4,'   Int(dg_pair):',f8.5,'   IBSI:',f8.5)") &
                iatm,a(iatm)%name,jatm,a(jatm)%name,dist,atmpairdg(idx,jdx),IBSImat(idx,jdx)
            end do
        end do
        
        !Between fragment
        if (allocated(frag1)) then
            if (natmlist==ncenter) then
	            IBSIfrag=0
	            do i=1,size(frag1)
		            do j=1,size(frag2)
			            IBSIfrag=IBSIfrag+IBSImat(frag1(i),frag2(j))
		            end do
	            end do
	            write(*,"(/,' The total IBSI between fragment 1 and 2:',f10.5)") IBSIfrag
            else
                write(*,"(/,a)") " Note: IBSI between the two defined fragments is not shown because the range &
                &of the atoms to be taken into account is not all atoms"
            end if
        end if
    end if
end do
end subroutine