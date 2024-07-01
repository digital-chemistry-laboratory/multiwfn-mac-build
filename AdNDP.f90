!! --------------- Perform Adaptive natural density partitioning
subroutine AdNDP
use defvar
use util
use GUI
use NAOmod
implicit real*8 (a-h,o-z)
integer,allocatable :: nNAOatm(:) !The number of NAOs of each atom
integer,allocatable :: atmcomb(:) !Store combination of specific number of atoms
integer,allocatable :: idxarray(:) !A contiguous array by default, used as underground numbering array for generating atom combinations
integer,allocatable :: searchlist(:) !Store atom indices for those will be exhaustively searched
real*8,allocatable :: orbeigval(:),orbeigvec(:,:),removemat(:,:),colvec(:,:),rowvec(:,:),DMNAOblk(:,:)
! real*8 :: bndcrit(30)=(/ 1.9D0,1.7D0,1.7D0,(1.8D0,i=4,30) /) 
integer :: numprint=100
real*8,allocatable :: candiocc(:),candivec(:,:) !Store candidate orbital information, occupation, eigenvector(row vector in NAOs)
integer,allocatable :: candinatm(:) !Number of atoms of candidate orbitals
integer,allocatable :: candiatmlist(:,:) !Store atom list of candidate orbitals
real*8,allocatable :: savedocc(:),savedvec(:,:) !Store saved orbital information, occupation, eigenvector(row vector in NAOs)
integer,allocatable :: savednatm(:) !Number of atoms of saved orbitals
integer,allocatable :: savedatmlist(:,:) !Store atom list of saved orbitals
real*8,allocatable :: oldsavedocc(:),oldsavedvec(:,:),oldDMNAO(:,:) !For temporarily store data
integer,allocatable :: oldsavednatm(:) !For temporarily store data
integer,allocatable :: oldsavedatmlist(:,:),eiguselist(:),tmparr(:)
real*8,allocatable :: adndpCObas(:,:),Fmat(:,:),Emat(:,:)
real*8,allocatable :: atmcomp(:),shcomp(:)
character :: c80tmp*80,c80tmp2*80,c200tmp*200,c2000tmp*2000,selectyn,fchfilename*200=' '

bndcrit=1.7D0  !Bond occupation threshold for different center-bonds

open(10,file=filename,status="old")

call checkNPA(ifound);if (ifound==0) return
write(*,*) "Loading NAO information..."
call loadNAOinfo
call loadNAOatminfo
allocate(nNAOatm(ncenter))
nNAOatm(:)=NAOend(:)-NAOinit(:)+1

call loclabel(10,'basis functions,',igauout)
if (igauout==1) then !Gaussian+NBO output file
	read(10,*) nbasis
	write(*,"(' The number of basis functions:',i6)") nbasis
    if (numNAO/=nbasis) then
        write(*,"(a)") " Warning: The number of basis functions is unequal to the number of NAOs, commonly this is because &
        &diffuse functions are used. If then Multiwfn fails to load data from the input file and crashes, or the result is weird, &
        &you should remove diffuse functions and regenerate the files, then try to redo the AdNDP analysis"
        write(*,*) "This situation is usually safe if NBO >=6.0 is used, while dangerous if NBO 3.1 is used"
        !I noticed that when linear dependency is occurred, although the number of NAOs is smaller than nbasis, &
        !the DMNAO printed by NBO 3.1 still has dimension of nbasis, making DMNAO incorrectly loaded by Multiwfn.
        !While for NBO 6, the dimension of DMNAO is always numNAO, as expected.
    end if
else !Cannot load from Gaussian output, assume the number of basis functions is identical to numNAO. However, when this is not true, cumbersome thing will occur...
	nbasis=numNAO
end if

write(*,*) "Loading DMNAO matrix..."
write(*,*)
call checkDMNAO(ifound);if (ifound==0) return
call loadDMNAO

if (iopshNAO==1) then !Open shell case
	write(*,*) "Use which density matrix for AdNDP analysis? 0=Total 1=Alpha 2=Beta"
	read(*,*) ispin
    if (ispin==1.or.ispin==2) bndcrit=bndcrit/2D0
    if (ispin==0) then
        continue !The current DMNAO is already total density matrix
    else if (ispin==1) then
        DMNAO=DMNAOa
    else if (ispin==2) then
        DMNAO=DMNAOb
    end if
else
    ispin=0
end if

close(10) !Loading finished

DMnele=0
do iNAO=1,numNAO
	!write(*,"(i5,3f12.6)") iNAO,DMNAO(iNAO,iNAO),DMNAOa(iNAO,iNAO),DMNAOb(iNAO,iNAO)
	DMnele=DMnele+DMNAO(iNAO,iNAO)
end do
write(*,"(' Number of electrons of this density matrix:',f12.5)") DMnele
!==== Remove core contribution from density matrix
removeocc=0
do iNAO=1,numNAO
	if (NAOset(iNAO,ispin)=="Cor") then
        DMNAO(iNAO,iNAO)=0D0
        removeocc=removeocc+NAOocc(iNAO,ispin)
    end if
end do
write(*,"(a,f10.3,a)") " Note: Contributions from core NAOs to density matrix (",removeocc," e) &
&have been eliminated by setting corresponding diagonal terms of density matrix to zero"
write(*,*) "Note: Default exhaustive search list is the entire system"
write(*,*)
!Initialization
nlencandi=100*ncenter !This is absolutely enough for any number of center searching
allocate(candiocc(nlencandi+1),candivec(nlencandi+1,numNAO),candiatmlist(nlencandi+1,ncenter),candinatm(nlencandi+1)) !The last element is used as temporary space to exchange information
allocate(eiguselist(nlencandi))
nlensaved=30*ncenter !Each atom can form at most four bonds, but we leave more space
allocate(savedocc(nlensaved),savedvec(nlensaved,numNAO),savedatmlist(nlensaved,ncenter),savednatm(nlensaved))
ncenana=1
ioutdetail=0
numsaved=0 !Number of picked AdNDP orbitals
numcandi=0 !Number of candidate orbitals
lensearchlist=ncenter !The list length of atom search range, default is entire system
allocate(colvec(numNAO,1),rowvec(1,numNAO))
allocate(atmcomb(ncenter),idxarray(ncenter),searchlist(ncenter)) !Note: Effective number of elements in searchlist is lensearchlist
forall (i=1:ncenter) searchlist(i)=i

isel=0
write(*,*) "      ======== Adaptive natural density partitioning (AdNDP) ========"
do while(.true.)
1	if (isel/=5.and.isel/=13.and.isel/=16) then
		!Sort candidate orbitals according to occupation varies from large to small and then print out them
		do i=1,numcandi
			do j=1,numcandi
				if (candiocc(i)>candiocc(j)) then !Exchange candidate orbital i and j, nlencandi+1 is used as temporary slot for exchanging, candinatm is needn't to be exchanged since they are the same
					candiocc(nlencandi+1)=candiocc(i)
					candivec(nlencandi+1,:)=candivec(i,:)
					candiatmlist(nlencandi+1,:)=candiatmlist(i,:)
					candiocc(i)=candiocc(j)
					candivec(i,:)=candivec(j,:)
					candiatmlist(i,:)=candiatmlist(j,:)
					candiocc(j)=candiocc(nlencandi+1)
					candivec(j,:)=candivec(nlencandi+1,:)
					candiatmlist(j,:)=candiatmlist(nlencandi+1,:)
				end if
			end do
		end do
		if (numcandi>0) then
			if (numcandi>numprint) then
                write(*,"(a,i6,a)") " Note: Only",numprint," candidate orbitals with highest occupancy are printed (This &
                &threshold can be changed via suboption 2 of option -2"
            end if
            write(*,*) "  ---- Current candidate orbital list, sorted according to occupation ----"
			do icandi=min(numcandi,numprint),1,-1 !Print from occupation of small to large, so index is decreased
				write(*,"(' #',i4,' Occ:',f7.4,' Atom:',9(i4,a))") icandi,candiocc(icandi),(candiatmlist(icandi,ii),a(candiatmlist(icandi,ii))%name,ii=1,candinatm(icandi))
			end do
			write(*,*)
 		else
 			write(*,*) "Note: Candidate orbital list is empty currently"
		end if
	end if

	!Check total number of electrons
	remainelec=0D0
	do iatm=1,lensearchlist
		do iNAO=NAOinit(searchlist(iatm)),NAOend(searchlist(iatm))
			remainelec=remainelec+DMNAO(iNAO,iNAO)
		end do
	end do
	write(*,"(' Residual valence electrons of all atoms in the search list:',f12.6)") remainelec

	write(*,*) "-10 Return to main menu"
    write(*,*) "-2 Various other settings and functions"
	write(*,*) "-1 Define exhaustive search list"
	if (numcandi>0) write(*,*) "0 Pick out some candidate orbitals and update occupations of others"
	write(*,"(' 1 Perform orbitals search for a specific atom combination')")
	write(*,"(' 2 Perform exhaustive search of ',i2,'-centers orbitals within the search list')") ncenana
	write(*,*) "3 Set the number of centers in the next exhaustive search"
	write(*,"(a,f8.3)") " 4 Set occupation threshold in the next exhaustive search, current:",bndcrit
	if (numsaved>0) write(*,"(' 5 Show information of AdNDP orbitals, current number:',i5)") numsaved
	if (numsaved>0) write(*,*) "6 Delete some AdNDP orbitals"
	if (numsaved>0) write(*,*) "7 Visualize AdNDP orbitals and molecular geometry"
	if (numsaved==0) write(*,*) "7 Visualize AdNDP orbitals (none) and molecular geometry"
	if (numcandi>0) write(*,*) "8 Visualize candidate orbitals and molecular geometry"
	if (numcandi==0) write(*,*) "8 Visualize candidate orbitals (none) and molecular geometry"
	if (numsaved>0) write(*,*) "9 Export some AdNDP orbitals to cube files"
	if (numcandi>0) write(*,"(a)") " 10 Export some candidate orbitals to Gaussian-type cube files"
	if (allocated(oldDMNAO)) write(*,"(a)") " 11 Save current density matrix and AdNDP orbital list again"
	if (.not.allocated(oldDMNAO)) write(*,"(a)") " 11 Save current density matrix and AdNDP orbital list (Unsaved)"
	if (allocated(oldDMNAO)) write(*,"(a)") " 12 Load saved density matrix and AdNDP orbital list"
	write(*,"(a)") " 13 Show residual density distributions on the atoms in the search list"
	if (numsaved>0) write(*,"(a)") " 14 Export AdNDP orbitals to AdNDP.mwfn file"
	write(*,"(a)") " 15 Evaluate and output composition of AdNDP orbitals"
	write(*,"(a)") " 16 Evaluate and output energy of AdNDP orbitals"
	read(*,*) isel
	
	
	if (isel==-10) then
		return
	else if (isel==-2) then
        do while(.true.)
            write(*,*)
            write(*,*) "0 Return"
	        if (ioutdetail==1) write(*,*) "1 Switch if output detail of exhaustive searching process, current: Yes"
	        if (ioutdetail==0) write(*,*) "1 Switch if output detail of exhaustive searching process, current: No"
            write(*,"(a,i5)") " 2 Set maximum number of candidate orbitals to be printed, current:",numprint
            write(*,*) "3 Output current density matrix to DMNAO.txt in current folder"
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
		        if (ioutdetail==1) then
			        ioutdetail=0
		        else
			        ioutdetail=1
		        end if
            else if (isel2==2) then
                write(*,*) "Input a number, e.g. 80"
                read(*,*) numprint
            else if (isel2==3) then
		        open(10,file="DMNAO.txt",status="replace")
		        call showmatgau(DMNAO,"Density matrix in NAO basis",0,"f12.6",10)
		        close(10)
		        write(*,*) "Done, density matrix in NAO basis has been saved to DMNAO.txt in current folder"
            end if
        end do
        
	else if (isel==-1) then
		call adndpdeflist(searchlist,lensearchlist)
        
	else if (isel==0) then
		write(*,*) "Input indices of the candidate orbitals to be picked out, e.g. 1,4-10,15,16"
		write(*,"(a)") " Note: If only inputting one number (N), then N orbitals with largest occupation will be picked out"
        write(*,"(a)") " To pick out e.g. the 7th orbital, input 7-7"
		read(*,"(a)") c2000tmp
		if (index(c2000tmp,',')/=0.or.index(c2000tmp,'-')/=0) then !Inputted indices
            call str2arr(c2000tmp,npickout)
			allocate(tmparr(npickout))
			call str2arr(c2000tmp,npickout,tmparr)
            write(*,"(i5,' orbitals are picked out')") npickout
		else
			read(c2000tmp,*) npickout
			allocate(tmparr(npickout))
            forall(i=1:npickout) tmparr(i)=i
		end if
		if (any(tmparr>numcandi)) then
			write(*,"(a)") " Error: Index of any orbital to be picked out should not be larger than the total number of candidate orbitals!"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            deallocate(tmparr)
			goto 1
		end if
        
		!Pick some orbitals from candidate list to permanent list
        do ipick=1,npickout
			iorb=tmparr(ipick)
			savedocc(numsaved+ipick)=candiocc(iorb)
			savedvec(numsaved+ipick,:)=candivec(iorb,:)
			savedatmlist(numsaved+ipick,:)=candiatmlist(iorb,:)
			savednatm(numsaved+ipick)=candinatm(iorb)
			!Deplete the density of the picked orbitals from DMNAO
			colvec(:,1)=candivec(iorb,:)
			rowvec(1,:)=candivec(iorb,:)
			DMNAO=DMNAO-candiocc(iorb)*matmul(colvec,rowvec)
        end do
		numsaved=numsaved+npickout
        
		!Make candidate list contiguous
        do ipick=1,npickout
			iorb=tmparr(ipick)
            if (iorb<numcandi) then
				candiocc(iorb:numcandi-1)=candiocc(iorb+1:numcandi)
				candivec(iorb:numcandi-1,:)=candivec(iorb+1:numcandi,:)
				candiatmlist(iorb:numcandi-1,:)=candiatmlist(iorb+1:numcandi,:)
				candinatm(iorb:numcandi-1)=candinatm(iorb+1:numcandi)
				where(tmparr>iorb) tmparr=tmparr-1
            end if
			numcandi=numcandi-1
        end do
        
        deallocate(tmparr)
        
		!Recalculate eigenval and eigenvec of remained candidate orbitals
		!Since an atom combination may have many orbitals exceed threshold, we first specify which eigenvalue of each combination will be used. If the combination is unique, then the largest one will be used
		nlenlist=candinatm(1) !The same to other candidate orbitals
		eiguselist(:)=0 !iwhichuse=n means the nth largest eigenvalue/eigenvector in this atom combination will be used. 0 means hasn't specified
		do icandi=1,numcandi
			iseleig=1
			if (eiguselist(icandi)/=0) cycle !Already specified
			do jcandi=1,numcandi
				if (all(candiatmlist(jcandi,1:nlenlist)==candiatmlist(icandi,1:nlenlist))) then
					eiguselist(jcandi)=iseleig
					iseleig=iseleig+1
				end if
			end do
		end do
		
		do icandi=1,numcandi
			nNAOblk=sum(nNAOatm(candiatmlist(icandi,1:candinatm(icandi)))) !Number of NAOs in current DMNAO block
			allocate(DMNAOblk(nNAOblk,nNAOblk),orbeigval(nNAOblk),orbeigvec(nNAOblk,nNAOblk))
			!Construct density matrix block
			irowed=0
			do idx1=1,candinatm(icandi) !Scan rows
				iatm=candiatmlist(icandi,idx1)
				irowbg=irowed+1
				irowed=irowbg+nNAOatm(iatm)-1
				icoled=0
				do idx2=1,candinatm(icandi) !Scan columns
					jatm=candiatmlist(icandi,idx2)
					icolbg=icoled+1
					icoled=icolbg+nNAOatm(jatm)-1
					DMNAOblk(irowbg:irowed,icolbg:icoled)=DMNAO(NAOinit(iatm):NAOend(iatm),NAOinit(jatm):NAOend(jatm))
				end do
			end do
			!Diagonalize the block
			call diagsymat(DMNAOblk,orbeigvec,orbeigval,istat)
			if (istat/=0) write(*,*) "Error: Diagonalization failed!"
			!Update candidate orbital
			ieiguse=nNAOblk-eiguselist(icandi)+1
			candiocc(icandi)=orbeigval(ieiguse) !The last element, namely nNAOblk, correspond to the largest occupation element
			candivec(icandi,:)=0D0 !Clean
			ied=0
			do idx=1,candinatm(icandi)
				iatm=candiatmlist(icandi,idx)
				ibg=ied+1
				ied=ibg+nNAOatm(iatm)-1
				candivec(icandi,NAOinit(iatm):NAOend(iatm))=orbeigvec(ibg:ied,ieiguse)
			end do
			deallocate(DMNAOblk,orbeigval,orbeigvec)
		end do
		write(*,"(i4,' candidate orbitals remain')") numcandi
		
	else if (isel==1.or.isel==2) then
		if (isel==1) then
			write(*,*) "Input atom indices, e.g. 3,4-6,7,12"
			write(*,*) "Note: Input ""a"" means all atoms in present system will be chosen"
			read(*,"(a)") c2000tmp
			if (index(c2000tmp,"a")/=0) then
				ncenana=ncenter
				forall(i=1:ncenter) atmcomb(i)=i
			else
                call str2arr(c2000tmp,ncenana,atmcomb)
				if (any(atmcomb(1:ncenana)>ncenter).or.any(atmcomb(1:ncenana)<=0)) then
					write(*,*) "Error: Some inputted atom indices exceeded valid range!"
					goto 1
				end if
			end if
		else if (isel==2) then
			forall (i=1:ncenana) idxarray(i)=i !Used as underground numbering index during generating combination
			atmcomb(1:ncenana)=searchlist(idxarray(1:ncenana))
			write(*,"(' Exhaustively searching ',i2,'-center orbitals, please wait...')") ncenana
		end if
		ipos=ncenana !Current position in the array
		numcandi=0 !Clean current candidate list
		ntotcomb=0 !Number of tried
		ioutcomb=1 !If do analysis for present combination
		cyccomb: do while(ipos>0)
			if (ioutcomb==1) then
				ntotcomb=ntotcomb+1
				!============Analyze atom combination in this time
				if (ioutdetail==1) write(*,"(' Searching atom',12i5)") atmcomb(1:ncenana)
				nNAOblk=sum(nNAOatm(atmcomb(1:ncenana))) !Number of NAOs in current DMNAO block
				allocate(DMNAOblk(nNAOblk,nNAOblk),orbeigval(nNAOblk),orbeigvec(nNAOblk,nNAOblk))
				!Construct density matrix block
				irowed=0
				do idx1=1,ncenana !Scan rows
					iatm=atmcomb(idx1)
					irowbg=irowed+1
					irowed=irowbg+nNAOatm(iatm)-1
					icoled=0
					do idx2=1,ncenana !Scan columns
						jatm=atmcomb(idx2)
						icolbg=icoled+1
						icoled=icolbg+nNAOatm(jatm)-1
						DMNAOblk(irowbg:irowed,icolbg:icoled)=DMNAO(NAOinit(iatm):NAOend(iatm),NAOinit(jatm):NAOend(jatm))
					end do
				end do
				!Diagonalize the block
				call diagsymat(DMNAOblk,orbeigvec,orbeigval,istat)
				if (istat/=0) write(*,*) "Error: Diagonalization failed!"
				if (isel==2.and.ioutdetail==1) then
					write(*,*) "Eigenvalues:"
					write(*,"(10f7.4)") orbeigval
				end if
				!Analyze result at this time
				do iNAO=nNAOblk,1,-1 !orbeigval varies from small to large, so cycle from large to small
					if (orbeigval(iNAO)>bndcrit.or.isel==1) then !When user specified combination, all orbitals will be outputted
						numcandi=numcandi+1
						if (isel==2.and.ioutdetail==1) write(*,"('Found the ',i4,'th candidate orbital with occupation:',f8.4)") numcandi,orbeigval(iNAO)
						if (numcandi>nlencandi) then
							write(*,"(a)") " Error: Candidate orbital list is overflowed! The direct way of solving this problem is &
                            &increasing ""nlencandi"" in AdNDP.f90 in Multiwfn source code package and then recompile the code. However, it is more likely that &
                            &you were searching AdNDP orbitals in an incorrect way, such as occupation threshold is too low, high-occupancy orbitals with lower &
                            &number of centers were not picked out, etc."
							write(*,*) "Press ENTER button to continue"
							read(*,*)
							deallocate(DMNAOblk,orbeigval,orbeigvec)
							numcandi=numcandi-1
							exit cyccomb
						end if
						!Move this orbital to candidate list
						candiocc(numcandi)=orbeigval(iNAO)
						candiatmlist(numcandi,1:ncenana)=atmcomb(1:ncenana)
						candinatm(numcandi)=ncenana
						candivec(numcandi,:)=0D0 !Clean
						ied=0
						do idx=1,ncenana
							iatm=atmcomb(idx)
							ibg=ied+1
							ied=ibg+nNAOatm(iatm)-1
							candivec(numcandi,NAOinit(iatm):NAOend(iatm))=orbeigvec(ibg:ied,iNAO)
						end do
					else
						exit
					end if
				end do
				deallocate(DMNAOblk,orbeigval,orbeigvec)
				if (ioutdetail==1) write(*,*)
				!============End analyze this combination
			end if
			if (isel==1) exit !isel==1 only do once for user inputted combination
			
			ioutcomb=0
			idxarray(ipos)=idxarray(ipos)+1
			if (idxarray(ipos)>lensearchlist) then
				ipos=ipos-1 !Go back to last position
				cycle
			end if
			if (ipos<ncenana) then
				ipos=ipos+1
				idxarray(ipos)=idxarray(ipos-1)
				cycle
			end if
			if (ipos==ncenana) then
				ioutcomb=1
				atmcomb(1:ncenana)=searchlist(idxarray(1:ncenana))
			end if
		end do cyccomb
		
		if (isel==2) write(*,"(' Tried',i9,' combinations, totally found',i9,' candidate orbitals',/)") ntotcomb,numcandi
		if (ncenana<lensearchlist.and.isel==2) ncenana=ncenana+1
        
	else if (isel==3) then
		write(*,"(a,i6)") " Input a number, should between 1 and",lensearchlist
		read(*,*) ncenanatmp
		if (ncenanatmp>lensearchlist) then
			write(*,*) "Error: The number of centers to be searched exceeds valid range!"
			goto 1
		end if
		ncenana=ncenanatmp
        
	else if (isel==4) then
		write(*,*) "Input a number, e.g. 1.9"
		read(*,*) bndcrit
        
	else if (isel==5) then
		write(*,*) "                         ---- AdNDP orbital list ----"
		do i=1,numsaved
			write(*,"(' #',i5,' Occ:',f8.4,' Atom:',9(i4,a))") i,savedocc(i),(savedatmlist(i,ii),a(savedatmlist(i,ii))%name,ii=1,savednatm(i))
		end do
		write(*,"(' Total occupation number in above orbitals:',f10.4,/)") sum(savedocc(1:numsaved))
        
	else if (isel==6) then
		write(*,*) "Input orbital index range that will be removed, e.g. 7,10"
		write(*,*) "Note: The density of these orbitals will not be returned to density matrix"
		read(*,*) ilow,ihigh
		if (ihigh<=numsaved) then
			numback=numsaved-ihigh
			savedocc(ilow:ilow+numback-1)=savedocc(ihigh+1:numsaved)
			savedvec(ilow:ilow+numback-1,:)=savedvec(ihigh+1:numsaved,:)
			savedatmlist(ilow:ilow+numback-1,:)=savedatmlist(ihigh+1:numsaved,:)
			savednatm(ilow:ilow+numback-1)=savednatm(ihigh+1:numsaved)
			numsaved=numsaved-(ihigh-ilow+1)
		else
			write(*,*) "Error: Index exceeded valid range"
		end if
		
	else if (isel==7.or.isel==8.or.isel==9.or.isel==10.or.isel==14) then !Visualize or export cube file for candidate or saved orbitals, or save them as .mwfn
		!Now we need basis functions information, load them from .fch file
        !At the same time, the CObas matrix generated by AONAO and AdNDP in NAO basis will be passed to "readfchadndp" to yield CO matrix used for real space visualization
		!If .fch or .fchk file with identical name in identical folder as initial input file can be found, then directly load it
		lenname=len_trim(filename)
		inquire(file=filename(1:lenname-3)//'fch',exist=alive)
		if (alive) then
			fchfilename=filename(1:lenname-3)//'fch'
		else
			inquire(file=filename(1:lenname-3)//'fchk',exist=alive)
			if (alive) fchfilename=filename(1:lenname-3)//'fchk'
		end if
		if (fchfilename==' ') then
			write(*,*) "Input path of corresponding .fch file, e.g. C:\test.fch"
			read(*,"(a)") fchfilename
			inquire(file=fchfilename,exist=alive)
			if (alive.eqv..false.) then
				write(*,*) "Error: File cannot be found! Hence orbitals cannot be visualized"
				write(*,*)
				fchfilename=' '
				goto 1
			end if
		end if
		if (.not.allocated(AONAO)) then
            open(10,file=filename,status="old")
            call checkAONAO(ifound);if (ifound==0) cycle
			call loadAONAO(nbasis)
            close(10)
		end if
		!Load mainbody of .fch file, and convert adndp orbitals (NAO basis) to CO matrix (GTF basis) so that fmo function can directly calculate orbital wavefunction value
		write(*,"(' Loading ',a)") trim(fchfilename)
		ifixorbsign=1 !Automatically fix sign of the isosurfaces generated by drawmolgui
        
		if (isel==7) then !Visualize saved orbitals
			allocate(adndpcobas(nbasis,numsaved))
			adndpcobas(:,:)=matmul(AONAO,transpose(savedvec(1:numsaved,:))) !cobasaadndp(i,j) means coefficient of basis function i in orbital j
			call readfchadndp(fchfilename,ispin,savedocc(1:numsaved),adndpcobas,numsaved)
			call drawmolgui
            write(*,*)
            
		else if (isel==8) then !Visualize candidate orbitals
			allocate(adndpcobas(nbasis,numcandi))
            !Original size: AONAO(nbasis,numNAO),candivec(nlencandi+1,numNAO)
			adndpcobas(:,:)=matmul( AONAO,transpose(candivec(1:numcandi,:)) )
			call readfchadndp(fchfilename,ispin,candiocc(1:numcandi),adndpcobas,numcandi)
			call drawmolgui
            write(*,*)
            
		else if (isel==14) then !Output all saved AdNDP orbitals as .mwfn file
			call dealloall(0)
			call readfch(fchfilename,1)
			wfntype=3
			CObasa=0
			MOocc=0
			MOene=0
			CObasa(:,1:numsaved)=matmul(AONAO,transpose(savedvec))
			MOocc(1:numsaved)=savedocc(1:numsaved)
			nmo=nbasis !AdNDP only performed for total density or single set of spin spin, therefore when nmo should be forced to equal to nbasis
			call outmwfn("AdNDP.mwfn",10,0)
			write(*,*) "Done! All AdNDP orbitals have been exported to AdNDP.mwfn in current folder"
			write(*,*)
            
		else if (isel==9.or.isel==10) then !Export saved or candidate AdNDP orbitals as cube file
			if (isel==9) then !Saved orbitals
				allocate(adndpcobas(nbasis,numsaved))
				adndpcobas(:,:)=matmul(AONAO,transpose(savedvec(1:numsaved,:))) !cobasaadndp(i,j) means coefficient of basis function i in orbital j
				call readfchadndp(fchfilename,ispin,savedocc(1:numsaved),adndpcobas,numsaved)
			else if (isel==10) then !Candidate orbitals
				allocate(adndpcobas(nbasis,numcandi))
				adndpcobas(:,:)=matmul(AONAO,transpose(candivec(1:numcandi,:)))
				call readfchadndp(fchfilename,ispin,candiocc(1:numcandi),adndpcobas,numcandi)
			end if
			!Set up grid
            call setgrid(0,igridsel)
			if (allocated(cubmat)) deallocate(cubmat)
			allocate(cubmat(nx,ny,nz))
				
			if (isel==9) then !Export saved AdNDP orbitals
                write(*,*)
				write(*,*) "Input index range of AdNDP orbitals to be exported"
                write(*,*) "e.g. 1-3,8,10-12 corresponds to 1,2,3,8,10,11,12"
                do while(.true.)
                    read(*,"(a)") c2000tmp
                    call str2arr(c2000tmp,ntmp)
                    allocate(tmparr(ntmp))
                    call str2arr(c2000tmp,ntmp,tmparr)
				    if (any(tmparr>numsaved)) then
					    write(*,*) "Error: Inputted index exceeded valid range! Input again"
					    deallocate(tmparr)
					else
                        exit
				    end if
                end do
				do itmp=1,ntmp
                    iorb=tmparr(itmp)
                    write(*,"(' Calculating grid data for orbital',i6,', please wait...')") iorb
					call savecubmat(4,1,iorb)
					if (sum(cubmat)<0) cubmat=-cubmat
					write(c80tmp,"('AdNDPorb',i4.4,'.cub')") iorb
					open(10,file=c80tmp,status="replace")
					call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
					close(10)
					write(*,"(1x,a,' has been exported to current folder')") trim(c80tmp)
				end do
			else if (isel==10) then !Export candidate orbitals
                write(*,*)
				write(*,*) "Input index range of candidate orbitals to be exported"
                write(*,*) "e.g. 1-3,8,10-12 corresponds to 1,2,3,8,10,11,12"
                do while(.true.)
                    read(*,"(a)") c2000tmp
                    call str2arr(c2000tmp,ntmp)
                    allocate(tmparr(ntmp))
                    call str2arr(c2000tmp,ntmp,tmparr)
				    if (any(tmparr>numcandi)) then
					    write(*,*) "Error: Inputted index exceeded valid range! Input again"
					    deallocate(tmparr)
					else
                        exit
				    end if
                end do
				do itmp=1,ntmp
                    iorb=tmparr(itmp)
                    write(*,"(' Calculating grid data for orbital',i6,', please wait...')") iorb
					call savecubmat(4,1,iorb)
					if (sum(cubmat)<0) cubmat=-cubmat
					write(c80tmp,"('candiorb',i4.4,'.cub')") iorb
					open(10,file=c80tmp,status="replace")
					call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
					close(10)
					write(*,"(1x,a,' has been exported to current folder')") trim(c80tmp)
				end do
			end if
			write(*,*)
			deallocate(cubmat,tmparr)
		end if
		if (allocated(adndpcobas)) deallocate(adndpcobas)
		
	else if (isel==11) then
		if (.not.allocated(oldDMNAO)) allocate(oldDMNAO(numNAO,numNAO))
		oldDMNAO=DMNAO
		if (allocated(oldsavedvec)) deallocate(oldsavedvec,oldsavedocc,oldsavedatmlist,oldsavednatm)
		allocate(oldsavedvec(numsaved,numNAO),oldsavedocc(numsaved),oldsavedatmlist(numsaved,ncenter),oldsavednatm(numsaved))
		oldsavedvec=savedvec(1:numsaved,:)
		oldsavedocc=savedocc(1:numsaved)
		oldsavedatmlist=savedatmlist(1:numsaved,:)
		oldsavednatm=savednatm(1:numsaved)
		noldorb=numsaved
		noldcenana=ncenana
		write(*,*) "Done, current density matrix in NAO basis and AdNDP orbital list has been saved"
        
	else if (isel==12) then
		if (.not.allocated(oldDMNAO)) then
			write(*,*) "Error: Density matrix in NAO basis has not been saved before!"
		else
			DMNAO=oldDMNAO
			numsaved=noldorb
			ncenana=noldcenana
			savedvec(1:numsaved,:)=oldsavedvec
			savedocc(1:numsaved)=oldsavedocc
			savedatmlist(1:numsaved,:)=oldsavedatmlist
			savednatm(1:numsaved)=oldsavednatm
			write(*,"(a)") " Done, the saved density matrix in NAO basis and AdNDP orbital list has been recovered"
            write(*,*)
		end if
        
	else if (isel==13) then
		write(*,*) "Residual valence electrons on each atom in the search list:"
		do idx=1,lensearchlist
			iatm=searchlist(idx)
			residatmdens=0
			do iNAO=NAOinit(iatm),NAOend(iatm)
				residatmdens=residatmdens+DMNAO(iNAO,iNAO)
			end do
			write(*,"(i8,a,':',f8.4)",advance='no') iatm,a(iatm)%name,residatmdens
			if (mod(iatm,4)==0) write(*,*)
		end do
		write(*,*)
        
	else if (isel==15) then !Output orbital composition of AdNDP orbitals. Adapted from NAOMO analysis module
        allocate(atmcomp(ncenter_NAO),shcomp(numNAOsh))
    	do while(.true.)
		    write(*,*)
		    write(*,*) "Analyze which orbital? e.g. 5"
            write(*,*) "Input 0 can return"
		    read(*,*) iorboutcomp
		    if (iorboutcomp==0) then
                if (allocated(atmcomp)) deallocate(atmcomp,shcomp)
                exit
		    else if (iorboutcomp<=0.or.iorboutcomp>numsaved) then
			    write(*,"(a,i7)") " Error: The orbital index should between  1 and",numsaved
			    cycle
		    end if
            
            write(*,"(a,f6.3,a)") " Terms whose absolute contribution > ",compthres," % are printed. &
            &This threshold can be changed by ""compthres"" in settings.ini"
            write(*,*)
		    write(*,*) "   NAO#   Center   Label      Type    Composition"
            shcomp=0D0
            atmcomp=0D0
            Rydcomp=0D0
		    do iNAO=1,numNAO
			    tmpcomp=savedvec(iorboutcomp,iNAO)**2*100
                if (NAOset(iNAO,ispin)=="Ryd") Rydcomp=Rydcomp+tmpcomp
                shcomp(bassh_NAO(iNAO))=shcomp(bassh_NAO(iNAO))+tmpcomp
                atmcomp(NAOcen(iNAO))=atmcomp(NAOcen(iNAO))+tmpcomp
			    if (abs(tmpcomp)<compthres) cycle !Skip showing too small terms
			    write(*,"(i8,i5,'(',a,')',4x,a,2x,a,'(',a,')',f10.3,'%' )") &
                iNAO,NAOcen(iNAO),NAOcenname(iNAO),NAOtype(iNAO),NAOset(iNAO,ispin),NAOshname(iNAO),tmpcomp
		    end do
		    write(*,*)
		    write(*,*) "Condensed NAO terms to shells:"
            do ish=1,numNAOsh
                icen=shcen_NAO(ish)
                if (abs(shcomp(ish))<compthres) cycle !Skip showing too small terms
                write(*,"('   Atom:',i6,'(',a,')  Shell:',i6,'(',a,1x,a,')',f10.3,'%' )") &
                icen,atmname_NAO(shcen_NAO(ish)),ish,shname_NAO(ish),shset_NAO(ish,ispin),shcomp(ish)
            end do
            write(*,*)
		    write(*,*) "Condensed NAO terms to atoms:"
		    write(*,*) "  Center   Composition"
		    do icen=1,ncenter_NAO
                if (abs(atmcomp(icen))<compthres) cycle !Skip showing too small terms
			    write(*,"(i6,'(',a,')',f10.3,'%' )") icen,atmname_NAO(icen),atmcomp(icen)
		    end do
		    write(*,*)
		    write(*,"(' Rydberg composition:',f10.3,'%')") Rydcomp
	    end do
    
	else if (isel==16) then !Evaluate orbital energy
		if (numsaved==0) then
			write(*,*) "Error: You need to pick out at least one orbital!"
			write(*,*)
			cycle
		end if
		!Load Fock matrix in AOs
		if (.not.allocated(Fmat)) then
			write(*,"(a)") " Input the file recording Fock matrix in original basis functions in lower triangular form, e.g. C:\fock.txt"
			write(*,*) "Note: If the suffix is .47, the Fock matrix will be directly loaded from it"
			read(*,"(a)") c200tmp
			inquire(file=c200tmp,exist=alive)
			if (.not.alive) then
				write(*,*) "Error: Unable to find this file!"
				cycle
			end if
			allocate(Fmat(nbasis,nbasis))
			open(10,file=c200tmp,status="old")
			if (index(c200tmp,".47")/=0) then
				call loclabel(10,"$FOCK",ifound)
				if (ifound==0) then
					write(*,*) "Error: Unable to find $FOCK field in this file!"
					close(10)
					cycle
				end if
				read(10,*)
				write(*,*) "Loading Fock matrix from .47 file..."
			end if
			read(10,*) ((Fmat(i,j),j=1,i),i=1,nbasis) !Load total or alpha Fock matrix
			if (iopshNAO==1) then !Open-shell
				if (ispin==0) then !User selected total density
					write(*,"(a,/)") " Error: This is an open-shell system but you selected analyzing total density, in this case orbital energy cannot be printed"
					cycle
				else if (ispin==2) then !User selected beta spin
					read(10,*) ((Fmat(i,j),j=1,i),i=1,nbasis) !Load beta Fock matrix
				end if
			end if
			do i=1,nbasis !Fill upper triangular part
 				do j=i+1,nbasis
 					Fmat(i,j)=Fmat(j,i)
 				end do
			end do
            close(10)
		end if
		if (.not.allocated(AONAO)) then
            open(10,file=filename,status="old")
            call checkAONAO(ifound);if (ifound==0) cycle
			call loadAONAO(nbasis)
            close(10)
		end if
		allocate(adndpcobas(nbasis,numsaved),Emat(numsaved,numsaved))
		!Note that savedvec is savedvec(numsaved,numNAO)
		adndpcobas=matmul(AONAO,transpose(savedvec(1:numsaved,:)))
		Emat=matmul(matmul(transpose(adndpcobas),Fmat),adndpcobas)
		if (allocated(MOene)) deallocate(MOene)
		allocate(MOene(numsaved))
		do iorb=1,numsaved
			MOene(iorb)=Emat(iorb,iorb)
		end do
		write(*,"(/,a)") " Energy of picked AdNDP orbitals:"
		do iorb=1,numsaved
			write(*,"(' Orbital:',i6,'  Energy (a.u./eV):',f12.6,f12.4)") iorb,MOene(iorb),MOene(iorb)*au2eV
		end do
		write(*,*)
		deallocate(adndpcobas,Emat,MOene)
	end if

end do

end subroutine




!!!------------- Define search list for exhaustive search
subroutine adndpdeflist(searchlist,lensearchlist)
use defvar
use util
implicit real*8 (a-h,o-z)
integer searchlisttmp(ncenter),searchlist(ncenter) !tmp verison is used to temporarily store index
integer lensearchlisttmp,lensearchlist !Effective length of the search list
integer tmparr(ncenter)
character cmd*200,elename*2

lensearchlisttmp=lensearchlist
searchlisttmp=searchlist
if (lensearchlisttmp>0) then
	write(*,"(' Currently',i5,' atoms are present in the search list:')") lensearchlisttmp
	do i=1,lensearchlisttmp
		write(*,"(i5,'(',a,')')",advance='no') searchlisttmp(i),a(searchlisttmp(i))%name
		if (mod(i,8)==0) write(*,*)
	end do
	write(*,*)
else
	write(*,*) "Current search list is empty"
end if
write(*,*)
write(*,*) "Exemplificative commands:"
write(*,*) "a 1,4,5,6 : Add atom 1,4,5,6 to the list"
write(*,*) "a 2-6     : Add atom 2,3,4,5,6 to the list"
write(*,*) "d 6,2,3   : Remove atom 6,2,3 from the list"
write(*,*) "d 2-6     : Remove atom 2,3,4,5,6 from the list"
write(*,*) "ae Al     : Add all aluminium atoms to the list"
write(*,*) "de H      : Remove all hydrogen atoms from the list"
write(*,*) "addall    : Add all atoms to the list"
write(*,*) "clean     : Clean the list"
write(*,*) "list      : Show current search list"
write(*,*) "help      : Show help information again"
write(*,*) "x         : Save the list and quit"
write(*,*) "q         : Quit without saving"

do while(.true.)
    write(*,"(/,a)") " Please input command. Press ENTER button directly can show help, inputting ""x"" can save and return"
	read(*,"(a)") cmd
	
    if (cmd==" ") then
		write(*,*) "Exemplificative commands:"
		write(*,*) "a 1,4,5,6 : Add atom 1,4,5,6 to the list"
		write(*,*) "a 2-6     : Add atom 2,3,4,5,6 to the list"
		write(*,*) "d 6,2,3   : Remove atom 6,2,3 from the list"
		write(*,*) "d 2-6     : Remove atom 2,3,4,5,6 from the list"
		write(*,*) "ae Al     : Add all aluminium atoms to the list"
		write(*,*) "de H      : Remove all hydrogen atoms from the list"
		write(*,*) "addall    : Add all atoms to the list"
		write(*,*) "clean     : Clean the list"
		write(*,*) "show      : Show current search list"
		write(*,*) "help      : Show help information again"
		write(*,*) "x         : Save the list and quit"
		write(*,*) "q         : Quit without saving"
	else if (cmd=="list") then
		if (lensearchlisttmp>0) then
			write(*,"(' Currently',i5,' atoms are present in the search list:')") lensearchlisttmp
			do i=1,lensearchlisttmp
				write(*,"(i5,'(',a,')')",advance='no') searchlisttmp(i),a(searchlisttmp(i))%name
				if (mod(i,8)==0) write(*,*)
			end do
			write(*,*)
		else
			write(*,*) "Current search list is empty"
		end if
		write(*,*)
	else if (cmd=='q') then
		write(*,*)
		exit
	else if (cmd=='x') then
		searchlist=searchlisttmp
		lensearchlist=lensearchlisttmp
		write(*,*) "Search list has been saved"
		if (lensearchlisttmp>0) then
			write(*,"('Currently',i5,' atoms are present in the search list:')") lensearchlisttmp
			do i=1,lensearchlisttmp
				write(*,"(i5,'(',a,')')",advance='no') searchlisttmp(i),a(searchlisttmp(i))%name
				if (mod(i,8)==0) write(*,*)
			end do
			write(*,*)
		else
			write(*,*) "Current search list is empty"
		end if
		write(*,*)
		exit
	else if (cmd=='addall') then
		lensearchlisttmp=ncenter
		forall (i=1:ncenter) searchlisttmp(i)=i
		write(*,*) "Done!"
	else if (cmd=='clean') then
		lensearchlisttmp=0
		write(*,*) "Done!"
	else if (cmd(1:2)=='ae'.or.cmd(1:2)=='de') then
			elename=cmd(4:5)
            call elename2idx(elename,iele)
			if (iele/=0) then
				if (cmd(1:2)=='ae') then !Add atom
					do icyclist=1,ncenter !Scan and find out corresponding atom from entire system
						if (a(icyclist)%index==iele) then
							if (any(searchlisttmp(1:lensearchlisttmp)==icyclist)) cycle !Check if it has presented in search list
							lensearchlisttmp=lensearchlisttmp+1
							searchlisttmp(lensearchlisttmp)=icyclist
						end if
					end do
				else !remove atom
					ipos=1
					do while(.true.)
						if (a(searchlisttmp(ipos))%index==iele) then
							if (lensearchlisttmp>=ipos+1) searchlisttmp(ipos:lensearchlisttmp-1)=searchlisttmp(ipos+1:lensearchlisttmp)
							lensearchlisttmp=lensearchlisttmp-1
						else
							ipos=ipos+1
						end if
						if (ipos>lensearchlisttmp) exit
					end do
				end if
				write(*,*) "Done!"
			else
				write(*,*) "Error: Unrecognizable element name "//elename
			end if
	else if (cmd(1:2)=='a '.or.cmd(1:2)=='d ') then
		if (index(cmd,'-')==0) then !Doesn't use range select for atoms
			iterm=1
			do i=1,len_trim(cmd)
				if (cmd(i:i)==',') iterm=iterm+1
			end do
			read (cmd(3:len_trim(cmd)),*) tmparr(1:iterm)
		else
			do i=1,len_trim(cmd) !Find position of -
				if (cmd(i:i)=='-') exit
			end do
			read(cmd(3:i-1),*) ilow
			read(cmd(i+1:),*) ihigh
			iterm=ihigh-ilow+1
			forall (i=1:iterm) tmparr(i)=i+ilow-1
		end if
		
		if (cmd(1:2)=='a') then
			do i=1,iterm
				if (any(searchlisttmp(1:lensearchlisttmp)==tmparr(i))) cycle
				lensearchlisttmp=lensearchlisttmp+1
				searchlisttmp(lensearchlisttmp)=tmparr(i)
			end do
		else
			ipos=1
			do while(.true.)
				if (any( tmparr(1:iterm)==searchlisttmp(ipos) )) then
					if (lensearchlisttmp>=ipos+1) searchlisttmp(ipos:lensearchlisttmp-1)=searchlisttmp(ipos+1:lensearchlisttmp)
					lensearchlisttmp=lensearchlisttmp-1
				else
					ipos=ipos+1
				end if
				if (ipos>lensearchlisttmp) exit
			end do
		end if
		write(*,*) "Done!"
	else
		write(*,*) "Error: Unrecognizable input"
	end if
end do
end subroutine