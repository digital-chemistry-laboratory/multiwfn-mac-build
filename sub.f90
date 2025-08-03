!------------- Modify & Check wavefunction
subroutine modwfn
use defvar
use util
implicit real*8 (a-h,o-z)
character seltmpc*10,selectyn,c1000tmp*1000,c2000tmp*2000
real*8 eigval(nbasis),eigvec(nbasis,nbasis),tmpmat(nbasis,nbasis)
real*8,allocatable :: tmparr(:)
integer orbarr(nmo)
integer,allocatable :: exclfragatm(:),tmparrint(:),idxsel(:)
character(len=3) :: orbtype(0:2)=(/ "A+B"," A "," B " /)
character symstr*6

write(*,*) "Note: ""GTF"" in this function refers to primitive Gaussian type function"
call updatenelec !Sometimes the naelec/nbelec/nelec recorded in the input file is not in line with sum of occupation numbers, so update
do while(.true.)
	write(*,*)
	write(*,*) "            ============ Modify & Check wavefunction ============ "
	write(*,"(' GTF:',i6,', Orb:',i6,', Atom:',i5,', A/B/tot ele:',3f10.3)") nprims,nmo,ncenter,naelec,nbelec,nelec
	if (ifragcontri/=1) write(*,*) "-4 Exclude contribution of some atoms to real space functions"
	if (ifragcontri/=1) write(*,*) "-3 Only retain contribution of some atoms to real space functions"
    if (ifPBC==3) write(*,"(a,3f10.6)") " -2 Set k-point, current:",kp1crd,kp2crd,kp3crd
	write(*,*) "-1 Return"
	write(*,*) "0 Save the present wavefunction to new.wfn file in current folder"
	if (allocated(CObasa)) then
		write(*,*) "1 List all GTFs                    2 List all basis functions"
	else
		write(*,*) "1 List all GTFs"
	end if
	write(*,*) "3 List all orbitals                4 Print detail information of an orbital"
	if (allocated(CObasa)) write(*,*) "5 Print coefficient matrix in basis functions"
	if (allocated(CObasa)) write(*,*) "6 Print density matrix in basis function"
	if (allocated(CObasa)) write(*,*) "7 Print various kinds of integral matrix between basis functions"
	write(*,*) "11 Exchange some information of two GTFs"
	write(*,*) "21 Set center of a GTF             22 Set type of a GTF"
	write(*,*) "23 Set exponent of a GTF           24 Set coefficient of a GTF in an orbital"
	if (allocated(CObasa)) then
		write(*,*) "25 Set coefficients of GTFs/basis functions that satisfied certain conditions"
	else
		write(*,*) "25 Set coefficients of GTFs that satisfied certain conditions"
	end if
	write(*,*) "26 Set occupation of some orbitals  27 Set type of some orbitals"
	write(*,*) "28 Set energy of some orbitals      29 Exchange information of two orbitals"
    write(*,*) "30 Exchange energies and occupation numbers for all orbitals"
	write(*,*) "31 Translate the system             32 Translate and duplicate the system"
	write(*,*) "33 Rotate wavefunction, namely X->Y, Y->Z, Z->X"
	if (imodwfn==0) write(*,*) "34 Set occupation number of inner orbitals to zero" !If occupation has been modified, don't do this to complicate things
	if (allocated(MOsym)) write(*,*) "35 Keep or discard orbital contributions according to irreducible rep."
	write(*,*) "36 Invert phase of some orbitals"
	write(*,*) "37 Split spatial orbitals as alpha and beta spin orbitals"
    write(*,*) "38 Make orbital occupations integer and satisfy Aufbau principle"
	read(*,*) isel
	
	if (isel==-1) then
		if (allocated(CObasa).and.imodwfn==1) then
			write(*,*) "Updating density matrix..."
			call genP
			write(*,*) "Density matrix has been updated"
		end if
        if (imodwfn==1) if_initlibreta=0 !LIBRETA should then be re-initialized
		exit
	else if (isel==-2) then
		write(*,*) "Input k-point coordinate in reciprocal space, e.g. 0.5,0.5,0"
        read(*,*) kp1crd,kp2crd,kp3crd
	else if (isel==-3.or.isel==-4) then
		if (allocated(fragatm)) deallocate(fragatm) !fragatm has been defined previously by default, fragatm contains all atoms
		if (isel==-3) then
			! "fragatm" is convertion relationship from fragment to the whole,
			! e.g. fragatm(4) is the actual atom index corresponding the 4th atom in fragment list
			write(*,"(a)") " Input atomic indices to define the fragment, e.g. 1,3-6,8,10-11 means atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,nfragatm)
			allocate(fragatm(nfragatm))
			call str2arr(c2000tmp,nfragatm,fragatm)
			call sort(fragatm,"val")
		else if (isel==-4) then
			write(*,*) "Input indices of the atoms to be excluded"
			write(*,*) "e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will be excluded"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,nexclatm)
			nfragatm=ncenter-nexclatm
			allocate(fragatm(nfragatm),exclfragatm(nexclatm))
			call str2arr(c2000tmp,nexclatm,exclfragatm)
			j=0
			do i=1,ncenter
				if (all(exclfragatm/=i)) then
					j=j+1
					fragatm(j)=i
				end if
			end do
		end if
		j=0
		do i=1,nprims
			if (any(fragatm==b(i)%center)) then
				j=j+1      !Move GTFs in the fragment to head of list
				CO(:,j)=CO(:,i)
				b(j)=b(i)
			end if
		end do
		ifragcontri=1 !Fragment has been defined by users
		write(*,"(' Done,',i8,' GTFs have been discarded,',i8,' GTFs reserved')") nprims-j,j
		nprims=j !Cut list at j, all functions after j seems nonexistent
		if (isel==-4) deallocate(exclfragatm)

		!Modification of wavefunction has finished, now reduce size of b, CO... to current nprims and nmo to avoid potential problems
		if (allocated(b)) then !Only for input file contains wavefunctions
			call resizebynmo(nmo,nprims) !Reduce size of CO, MOene, MOocc, MOtype
			allocate(b_tmp(nprims))
			b_tmp(:)=b(1:nprims)
			deallocate(b)
			allocate(b(nprims))
			b=b_tmp
			deallocate(b_tmp)
		end if
        
	else if (isel==0) then
		call outwfn("new.wfn",1,1,10)
		write(*,*) "Wavefunction has been outputted to new.wfn in current folder"
	
	else if (isel==1) then
		GTFlarge=0
		do i=1,nprims
			write(*,"(i6,' Center:',i5,'(',a2,')','   Type: ',a,'   Exponent:',E16.7)") i,b(i)%center,a(b(i)%center)%name,GTFtype2name(b(i)%type),b(i)%exp
            if (b(i)%exp>GTFlarge) then
				ilarge=i
                GTFlarge=b(i)%exp
            end if
		end do
        write(*,"(/,' Largest exponent is',E16.7,' of GTF',i6)") GTFlarge,ilarge
	
	else if (isel==2) then
		do i=1,nbasis
            if (isphergau==1) then
			    write(*,"(' Basis:',i5,'   Shell:',i5,'   Center:',i5,'(',a2,')   Type:',a)")&
                i,basshell(i),bascen(i),a(bascen(i))%name,GTFtype2name(bastype(i))
            else !The primstart/end is constructed w.r.t. Cartesian basis functions
			    write(*,"(' Basis:',i5,'   Shell:',i5,'   Center:',i5,'(',a2,')   Type:',a,'  GTF:',i6,' to',i6)")&
                i,basshell(i),bascen(i),a(bascen(i))%name,GTFtype2name(bastype(i)),primstart(i),primend(i)
            end if
		end do
	
	else if (isel==3) then
		write(*,*) "Basic information of all orbitals:"
		symstr=" "
		naorb=count(MOtype==1)
		do i=1,nmo
			if (allocated(MOsym)) symstr='('//MOsym(i)//')'
			if (wfntype==0.or.wfntype==2.or.wfntype==3) then
				write(*,"(' Orb:',i6,' Ene(au/eV):',f13.6,f13.4,' Occ:',f9.6,' Type:',a,1x,a)") &
				i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
			else
				if (MOtype(i)==1) then
					write(*,"(i6,9x,' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
					i,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
				else
					write(*,"(i6,' (',i6,')',' E(au/eV):',f12.5,f13.4,' Occ:',f9.6,' Typ:',a,1x,a)") &
					i,i-naorb,MOene(i),MOene(i)*au2eV,MOocc(i),orbtype(MOtype(i)),symstr
				end if
			end if
		end do
		if (any(MOtype==2)) write(*,"(a)") " Note: For beta orbitals, &
		&the index in the parenthese shown above is the index counted from the first beta orbital"
		
	else if (isel==4) then
		write(*,*) "Input the orbital index, e.g. 12"
		read(*,*) i
		if (i<1.or.i>nmo) then
			write(*,"(' Invalid orbital index, should within range of',i5,' and ',i5)") 1,nmo
		else
			write(*,"(' Occupation number is ',f12.7,'     Energy is',f12.6,' Hartree')") MOocc(i),MOene(i)
			if (MOtype(i)==0) write(*,*) "This is a closed-shell orbital"
			if (MOtype(i)==1) write(*,*) "This is an alpha orbital"
			if (MOtype(i)==2) write(*,*) "This is a beta orbital"
			write(*,*)
			do j=1,nprims
				write(*,"(' GTF:',i6,' Cen:',i5,'(',a2,')',' Type: ',a,' Coeff:',1PE16.8,' Exp: ',1PE13.7)") &
				j,b(j)%center,a(b(j)%center)%name,GTFtype2name(b(j)%type),CO(i,j),b(j)%exp
			end do
			write(*,"(a,/)") " Note: The ""coeff."" are expansion coefficients of orbitals with respect to GTFs, including normalization constant"
			if (allocated(b)) then
				do j=1,nbasis
					if (MOtype(i)==0.or.MOtype(i)==1) covalue=CObasa(j,i)
					if (MOtype(i)==2) covalue=CObasb(j,i-nbasis)
					write(*,"(' Basis func:',i6,'  Cen:',i5,'(',a2,')',' Shell:',i5,' Type: ',a,' Coeff:',f12.8)") &
					j,bascen(j),a(bascen(j))%name,basshell(j),GTFtype2name(bastype(j)),covalue
				end do
			end if
			write(*,"(a,/)") " Note: The ""coeff."" are expansion coefficients of orbitals with respect to basis functions, which are normalized functions"
		end if
	
	else if (isel==5) then
		write(*,*) "0 Return"
		write(*,*) "1 Print on screen"
		write(*,*) "2 Print to Cmat.txt in current folder"
		read(*,*) iseltmp
		if (iseltmp==1.or.iseltmp==2) then
			if (iseltmp==1) ides=6
			if (iseltmp==2) then
				ides=10
				open(ides,file="Cmat.txt",status="replace")
			end if
			write(ides,*) "Note: (i,j) element means coefficient of ith basis function in jth orbital"
			if (wfntype==0.or.wfntype==2.or.wfntype==3) then
				call showmatgau(CObasa,"Coefficient matrix",0,fileid=ides)
			else if (wfntype==1.or.wfntype==4) then
				call showmatgau(CObasa,"Alpha coefficient matrix",0,fileid=ides)
				call showmatgau(CObasb,"Beta coefficient matrix",0,fileid=ides)
			end if
			if (iseltmp==2) then
				write(*,*) "Done! The matrix has been outputted to Cmat.txt in current folder"
				close(ides)
			end if
		end if
	
	else if (isel==6) then
		write(*,*) "0 Return"
		write(*,*) "1 Print on screen"
		write(*,*) "2 Print to Pmat.txt in current folder"
		read(*,*) iseltmp
		if (iseltmp==1.or.iseltmp==2) then
			if (iseltmp==1) ides=6
			if (iseltmp==2) then
				ides=10
				open(ides,file="Pmat.txt",status="replace")
			end if
			call showmatgau(Ptot,"Total density matrix",1,fileid=ides)
			sumt=0
			do i=1,nbasis
				sumt=sumt+Ptot(i,i)
			end do
			write(ides,"(/,' Trace of density matrix:',f14.6)") sumt
			write(ides,"(' Trace of density matrix multiplied by overlap matrix:',f14.6)") sum(Ptot*Sbas)
			if (wfntype==1.or.wfntype==2.or.wfntype==4) then
				suma=0
				sumb=0
				do i=1,nbasis
					suma=suma+Palpha(i,i)
					sumb=sumb+Pbeta(i,i)
				end do
				write(ides,*)
				call showmatgau(Palpha-Pbeta,"Spin density matrix",1,fileid=ides)
				write(ides,*)
				call showmatgau(Palpha,"Alpha density matrix",1,fileid=ides)
				write(ides,*)
				call showmatgau(Pbeta,"Beta density matrix",1,fileid=ides)
				write(ides,*)
				write(ides,"(' Trace of alpha and beta density matrix:',2f14.6)") suma,sumb
				write(ides,"(' Trace of density matrix multiplied by overlap matrix:',/,' Alpha:',f14.6,'     Beta:',f14.6)") sum(Palpha*Sbas),sum(Pbeta*Sbas)
			end if
			if (iseltmp==2) then
				write(*,*) "Done! The matrix has been outputted to Pmat.txt in current folder"
				close(ides)
			end if
		end if
	
	!Print various kinds of integral matrix between basis functions
	else if (isel==7) then
		write(*,*) "Print which kind of integral matrix?"
		write(*,*) "0 Fock/KS matrix"
		write(*,*) "1 Overlap integral"
		write(*,*) "2 Electric dipole moment integral"
		write(*,*) "3 Magnetic dipole moment integral"
		write(*,*) "4 Velocity integral"
		write(*,*) "5 Kinetic energy integral"
		write(*,*) "6 Electric quadrupole moment integral"
		write(*,*) "7 Electric octopole moment integral"
		write(*,*) "8 Electric hexadecapole moment integral"
		read(*,*) imattype
        iout=2
        if (imattype<=5) then !Amount of quadrupole and octopole is too large, and only be printed to file
		    write(*,*) "Select destination of outputting"
		    write(*,*) "1 Print on screen"
		    write(*,*) "2 Print to intmat.txt in current folder"
		    read(*,*) iout
        end if
		if (iout==1) then
            ides=6
		else
			ides=10
			open(ides,file="intmat.txt",status="replace")
		end if
		if (imattype==0) then
			if (wfntype==3.or.wfntype==4) then
				write(*,*) "Error: Fock/KS matrix is only available for single-determinant wavefunction"
                cycle
            end if
			if (.not.allocated(FmatA)) then
				do while(.true.)
					write(*,*)
					write(*,*) "How to provide the Fock/KS matrix?"
					write(*,"(a)") " 1 Generating Fock/KS matrix by MO energies and coefficients as well as overlap matrix"
					write(*,*) "2 Loading Fock/KS matrix from a file"
					read(*,*) isel2
					if (isel2==1) then
						call MOene2Fmat(istatus)
					else if (isel2==2) then
						call loadFockfile(istatus)
					end if
					if (istatus==0) exit
				end do
			end if
			write(ides,*)
            if (wfntype==0.or.wfntype==2) then
				call showmatgau(FmatA(:,:),"Fock/KS matrix",1,fileid=ides)
            else if (wfntype==1) then
				call showmatgau(FmatA(:,:),"Fock/KS matrix of alpha spin",1,fileid=ides)
				write(ides,*)
				call showmatgau(FmatB(:,:),"Fock/KS matrix of beta spin",1,fileid=ides)
            end if
		else if (imattype==1) then
			call ask_Sbas_PBC
			call showmatgau(Sbas,"Overlap matrix",1,fileid=ides)
            tmpmat=Sbas
			call diagsymat(Sbas,eigvec,eigval,ierror)
			write(ides,*)
			write(ides,*) "Eigenvalues:"
			write(ides,"(6f12.8)") eigval
            Sbas=tmpmat
		else if (imattype==2) then
			if (.not.allocated(Dbas)) then
                write(*,*) "Calculating the matrix..."
                call genDbas_curr
            end if
			write(ides,*)
			call showmatgau(Dbas(1,:,:),"Electric dipole moment matrix (X component)",1,fileid=ides)
			write(ides,*)
			call showmatgau(Dbas(2,:,:),"Electric dipole moment matrix (Y component)",1,fileid=ides)
			write(ides,*)
			call showmatgau(Dbas(3,:,:),"Electric dipole moment matrix (Z component)",1,fileid=ides)
		else if (imattype==3) then
			if (.not.allocated(Magbas)) then
                write(*,*) "Calculating the matrix..."
                call genMagbas_curr
            end if
			write(ides,*)
			call showmatgau(Magbas(1,:,:),"Magnetic dipole moment matrix (X component)",0,fileid=ides)
			write(ides,*)
			call showmatgau(Magbas(2,:,:),"Magnetic dipole moment matrix (Y component)",0,fileid=ides)
			write(ides,*)
			call showmatgau(Magbas(3,:,:),"Magnetic dipole moment matrix (Z component)",0,fileid=ides)
		else if (imattype==4) then
			if (.not.allocated(Velbas)) then
                write(*,*) "Calculating the matrix..."
                call genVelbas_curr
            end if
			write(ides,*)
			call showmatgau(Velbas(1,:,:),"Velocity matrix (X component)",0,fileid=ides)
			write(ides,*)
			call showmatgau(Velbas(2,:,:),"Velocity matrix (Y component)",0,fileid=ides)
			write(ides,*)
			call showmatgau(Velbas(3,:,:),"Velocity matrix (Z component)",0,fileid=ides)
        else if (imattype==5) then
			if (.not.allocated(Tbas)) then
                write(*,*) "Calculating the matrix..."
                call genTbas_curr
            end if
			write(ides,*)
			call showmatgau(Tbas(:,:),"Kinetic energy matrix",1,fileid=ides)
        else if (imattype==6) then
			if (.not.allocated(Quadbas)) then
                write(*,*) "Calculating the matrix..."
                call genMultipolebas_curr
            end if
			              call showmatgau(Quadbas(1,:,:),"Quadrupole moment matrix (XX component)",1,fileid=ides)
			write(ides,*);call showmatgau(Quadbas(2,:,:),"Quadrupole moment matrix (YY component)",1,fileid=ides)
			write(ides,*);call showmatgau(Quadbas(3,:,:),"Quadrupole moment matrix (ZZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(Quadbas(4,:,:),"Quadrupole moment matrix (XY component)",1,fileid=ides)
			write(ides,*);call showmatgau(Quadbas(5,:,:),"Quadrupole moment matrix (YZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(Quadbas(6,:,:),"Quadrupole moment matrix (XZ component)",1,fileid=ides)
        else if (imattype==7) then
			if (.not.allocated(Octobas)) then
                write(*,*) "Calculating the matrix..."
                call genMultipolebas_curr
            end if
			              call showmatgau(Octobas(1,:,:),"Octopole moment matrix (XXX component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(2,:,:),"Octopole moment matrix (YYY component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(3,:,:),"Octopole moment matrix (ZZZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(4,:,:),"Octopole moment matrix (YZZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(5,:,:),"Octopole moment matrix (XZZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(6,:,:),"Octopole moment matrix (XXZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(7,:,:),"Octopole moment matrix (YYZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(8,:,:),"Octopole moment matrix (XXY component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(9,:,:),"Octopole moment matrix (XYY component)",1,fileid=ides)
			write(ides,*);call showmatgau(Octobas(10,:,:),"Octopole moment matrix (XYZ component)",1,fileid=ides)
        else if (imattype==8) then
			if (.not.allocated(hexdebas)) then
                write(*,*) "Calculating the matrix..."
                call genMultipolebas_curr
            end if
			              call showmatgau(hexdebas(1,:,:),"Hexadecapole moment matrix (XXXX component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(2,:,:),"Hexadecapole moment matrix (YYYY component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(3,:,:),"Hexadecapole moment matrix (ZZZZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(4,:,:),"Hexadecapole moment matrix (XXXY component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(5,:,:),"Hexadecapole moment matrix (XXXZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(6,:,:),"Hexadecapole moment matrix (YYYX component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(7,:,:),"Hexadecapole moment matrix (YYYZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(8,:,:),"Hexadecapole moment matrix (ZZZX component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(9,:,:),"Hexadecapole moment matrix (ZZZY component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(10,:,:),"Hexadecapole moment matrix (XXYY component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(11,:,:),"Hexadecapole moment matrix (XXZZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(12,:,:),"Hexadecapole moment matrix (YYZZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(13,:,:),"Hexadecapole moment matrix (XXYZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(14,:,:),"Hexadecapole moment matrix (YYXZ component)",1,fileid=ides)
			write(ides,*);call showmatgau(hexdebas(15,:,:),"Hexadecapole moment matrix (ZZXY component)",1,fileid=ides)
		end if
		if (iout==2) then
			write(*,*) "Done! The matrix has been outputted to intmat.txt in current folder"
			close(ides)
		end if
		
	else if (isel==11) then
		write(*,*) "Exchange information of which two GTFs? Input their indices, e.g. 18,21"
		read(*,*) i,j
		write(*,*) "Exchange which information for the two GTFs?"
		write(*,*) "1 Exchange all properties"
		write(*,*) "2 Exchange center"
		write(*,*) "3 Exchange function type"
		write(*,*) "4 Exchange exponent"
		write(*,*) "5 Exchange orbital expansion coefficient"
		read(*,*) iswapcontent
		if (iswapcontent==1) call swapGTF(i,j,"all")
		if (iswapcontent==2) call swapGTF(i,j,"cen")
		if (iswapcontent==3) call swapGTF(i,j,"typ")
		if (iswapcontent==4) call swapGTF(i,j,"exp")
		if (iswapcontent==5) call swapGTF(i,j,"MO ")
		write(*,*) "Exchanging finished!"
	
	else if (isel==21) then
		write(*,*) "Input index of the GTF, e.g. 13"
		read(*,*) i
		write(*,*) "Input the index of the center that you want the function centered at, e.g. 5"
		read(*,*) j
		if (j<=ncenter.and.j>0) then
			b(i)%center=j
            write(*,*) "Done!"
		else
			write(*,"(' Error: The center index should >0 and <=',i7)") ncenter
		end if
	
	else if (isel==22) then
		write(*,*) "Input index of the GTF, e.g. 13"
		read(*,*) i
		write(*,*) "Input the type you want to set, should be one of the following ones"
		write(*,*) "S,X,Y,Z,XX,YY,ZZ,XY,XZ,YZ,XXX,YYY,ZZZ,XXY,XXZ,YYZ,XYY,XZZ,YZZ,XYZ"
		write(*,*) "ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX"
		write(*,"(a)") " ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX"
		read(*,*) seltmpc
		do j=1,size(GTFtype2name)
			if (seltmpc==GTFtype2name(j)) then
				b(i)%type=j
				write(*,*) "Done!"
				exit
			end if
			if (j==20) write(*,*) "Error: Cannot recognize this type"
		end do
	
	else if (isel==23) then
		write(*,*) "Input index of the GTF, e.g. 13"
		read(*,*) i
		write(*,*) "Input exponent, e.g. 0.035"
		read(*,*) rexp
		b(i)%exp=rexp
        write(*,*) "Done!"
	
	else if (isel==24) then
		write(*,*) "Input index of the GTF, e.g. 13"
		read(*,*) iprm
		write(*,*) "Input orbital index, e.g. 12"
		read(*,*) imonum
		if (iprm<=nprims.and.iprm>0.and.imonum<=nmo.and.imonum>0) then
			write(*,*) "Input the coefficient, e.g. 2.83E-3"
			read(*,*) rcoeff
			CO(imonum,iprm)=rcoeff
            write(*,*) "Done!"
		else
			write(*,"(' Error: The index of the GTF or orbital exceeded valid range!')")
		end if
	
	else if (isel==25) then
		isetmode=1
		if (allocated(CObasa)) then
			write(*,*) "1 Set coefficients of some GTFs in some orbitals"
			write(*,*) "2 Set coefficients of some basis functions in some orbitals"
			read(*,*) isetmode
		end if
		if (isetmode==1) then
			allocate(idxsel(nprims)) 
            idxsel=1 !If idxsel(i)=1, then GTF i is selected
			write(*,"(a)") " You will be asked to input several conditions, the GTFs simultaneously satisfying them will be selected"
			write(*,*)
			write(*,*) "Input indices that the GTFs must be, e.g. 2,3,7-10"
            write(*,*) "Press ENTER button directly means ignoring this condition"
			read(*,"(a)") c2000tmp
            if (c2000tmp/=" ") then
				call str2arr(c2000tmp,ntmp)
				allocate(tmparrint(ntmp))
				call str2arr(c2000tmp,ntmp,tmparrint)
				do iGTF=1,nprims
					if (all(tmparrint/=iGTF)) idxsel(iGTF)=0
				end do
				deallocate(tmparrint)
            end if
			write(*,*) "Input indices of the atoms that the GTFs must be centered at, e.g. 3,9-13,18"
            write(*,*) "Press ENTER button directly means ignoring this condition"
			read(*,"(a)") c2000tmp
            if (c2000tmp/=" ") then
				call str2arr(c2000tmp,ntmp)
				allocate(tmparrint(ntmp))
				call str2arr(c2000tmp,ntmp,tmparrint)
				do iGTF=1,nprims
					if (all(tmparrint/=b(iGTF)%center)) idxsel(iGTF)=0
				end do
				deallocate(tmparrint)
            end if
			write(*,*) "Input the type of the GTFs must be, should be one of S,X,Y,Z,XX,XY"
			write(*,*) "You can also input one of S,P,D,F,G,H to select according to angular moment"
            write(*,"(a)") " Adding a negative sign means selecting all other types, for example, ""-P"" means selecting all GTFs other than P type"
            write(*,*) "Press ENTER button directly means ignoring this condition"
			read(*,"(a)") seltmpc
            if (seltmpc/=" ") then
				if (seltmpc(1:1)=='-') then
					seltmpc=seltmpc(2:)
					do iGTF=1,nprims
						if ( GTFtype2name(b(iGTF)%type)==trim(seltmpc).or.type2ang(b(iGTF)%type)==trim(seltmpc) ) idxsel(iGTF)=0
					end do
                else
					do iGTF=1,nprims
						if ( GTFtype2name(b(iGTF)%type)/=trim(seltmpc).and.type2ang(b(iGTF)%type)/=trim(seltmpc) ) idxsel(iGTF)=0
					end do
                end if
            end if
			write(*,"(' Coefficient of',i8,' GTFs have been selected')") count(idxsel==1)
			write(*,"(/,a)") " Input indices of the orbitals for which the coefficients of selected GTFs will be set, e.g. 2,9-13,18-25"
            write(*,*) "Press ENTER button directly selecting all orbitals"
			read(*,"(a)") c2000tmp
            if (c2000tmp/=" ") then
				call str2arr(c2000tmp,ntmp)
				allocate(tmparrint(ntmp))
				call str2arr(c2000tmp,ntmp,tmparrint)
            else
				ntmp=nmo
				allocate(tmparrint(ntmp))
                forall (iorb=1:nmo) tmparrint(iorb)=iorb
            end if
			write(*,*) "Input the expansion coefficient you want to set, e.g. 0.5"
			read(*,*) coval
            do iGTF=1,nprims
				if (idxsel(iGTF)==0) cycle
                CO(tmparrint(1:ntmp),iGTF)=coval
            end do
            deallocate(tmparrint,idxsel)
            write(*,*) "Done!"
		else if (isetmode==2) then
			allocate(idxsel(nbasis)) 
            idxsel=1 !If idxsel(i)=1, then basis function i is selected
			write(*,"(a)") " You will be asked to input several conditions, the basis function simultaneously satisfying them will be selected"
			write(*,*)
			write(*,*) "Input indices that the basis functions must be, e.g. 2,3,7-10"
            write(*,*) "Press ENTER button directly means ignoring this condition"
			read(*,"(a)") c2000tmp
            if (c2000tmp/=" ") then
				call str2arr(c2000tmp,ntmp)
				allocate(tmparrint(ntmp))
				call str2arr(c2000tmp,ntmp,tmparrint)
				do ibas=1,nbasis
					if (all(tmparrint/=ibas)) idxsel(ibas)=0
				end do
				deallocate(tmparrint)
            end if
			write(*,*) "Input indices of the atoms that the basis functions must be centered at"
            write(*,*) "For example, 3,9-13,18"
            write(*,*) "Press ENTER button directly means ignoring this condition"
			read(*,"(a)") c2000tmp
            if (c2000tmp/=" ") then
				call str2arr(c2000tmp,ntmp)
				allocate(tmparrint(ntmp))
				call str2arr(c2000tmp,ntmp,tmparrint)
				do ibas=1,nbasis
					if (all(tmparrint/=bascen(ibas))) idxsel(ibas)=0
				end do
				deallocate(tmparrint)
            end if
			write(*,*) "Input the type of the basis functions must be, should be one of S,X,Y,Z,XX,XY"
			write(*,*) "You can also input one of S,P,D,F,G,H to select according to angular moment"
            write(*,"(a)") " Adding a negative sign means selecting all other types, for example, ""-P"" means selecting all basis functions other than P type"
            write(*,*) "Press ENTER button directly means ignoring this condition"
			read(*,"(a)") seltmpc
            if (seltmpc/=" ") then
				if (seltmpc(1:1)=='-') then
					seltmpc=seltmpc(2:)
					do ibas=1,nbasis
						if ( GTFtype2name(bastype(ibas))==trim(seltmpc).or.type2ang(bastype(ibas))==trim(seltmpc) ) idxsel(ibas)=0
					end do
                else
					do ibas=1,nbasis
						if ( GTFtype2name(bastype(ibas))/=trim(seltmpc).and.type2ang(bastype(ibas))/=trim(seltmpc) ) idxsel(ibas)=0
					end do
                end if
            end if
			write(*,"(' Coefficient of',i8,' basis functions have been selected')") count(idxsel==1)
			write(*,"(/,a)") " Input indices of the orbitals for which the coefficients of selected basis functions will be set, e.g. 2,9-13,18-25"
            write(*,*) "Press ENTER button directly selecting all orbitals"
			read(*,"(a)") c2000tmp
            if (c2000tmp/=" ") then
				call str2arr(c2000tmp,ntmp)
				allocate(tmparrint(ntmp))
				call str2arr(c2000tmp,ntmp,tmparrint)
            else
				ntmp=nmo
				allocate(tmparrint(ntmp))
                forall (iorb=1:nmo) tmparrint(iorb)=iorb
            end if
			write(*,*) "Input the expansion coefficient you want to set, e.g. 0.5"
			read(*,*) coval
            do ibas=1,nbasis
				if (idxsel(ibas)==0) cycle
                do idxorb=1,ntmp
					iorb=tmparrint(idxorb)
                    if (iorb<=nbasis) then
						CObasa(ibas,iorb)=coval
                    else
						CObasb(ibas,iorb-nbasis)=coval
                    end if
                end do
            end do
            deallocate(tmparrint,idxsel)
            write(*,*) "Done!"
			imodwfn=1
		end if

	else if (isel==26) then
		call modorbocc
	
	else if (isel==27) then
		do while(.true.)
			write(*,*) "Set type for which range of orbitals?"
			write(*,*) "e.g. 2,4,13-16,20 means selecting orbitals 2,4,13,14,15,16,20"
			write(*,*) "Input 0 can select all orbitals, input q or 00 can return"
			read(*,"(a)") c1000tmp
			if (c1000tmp(1:1)=='q'.or.c1000tmp(1:2)=='00') exit
			if (c1000tmp(1:1)=='0') then
				numorbsel=nmo
				do i=1,nmo
					orbarr(i)=i
				end do
			else
				call str2arr(c1000tmp,numorbsel,orbarr)
				if ( any(orbarr(1:numorbsel)<1).or.any(orbarr(1:numorbsel)>nmo) ) then
					write(*,*) "Error: One or more orbital indices exceeded valid range!"
					cycle
				end if
			end if
			write(*,*) "Set to which type?  0=Alpha+Beta  1=Alpha  2=Beta"
			read(*,*) isettype
			MOtype(orbarr(1:numorbsel))=isettype
			!Recount alpha and beta electrons
			call updatenelec
			write(*,*) "Done!"
			!Update wavefunction type
			if (all(MOtype==0)) then
				if (all(MOocc==nint(MOocc))) then !All A+B orbital & integer occupation
					wfntype=0
					write(*,"(' Note: Now the wavefunction is recognized as a restricted closed-shell single-determinant wavefunction')")
				else !All A+B orbital & partial occupation
					wfntype=3
					write(*,"(' Note: Now the wavefunction is recognized as a restricted multiconfiguration wavefunction')")
				end if
			else
				if (any(MOocc/=nint(MOocc)).and.all(MOtype/=0)) then !Either A or B, and partial occupation
					wfntype=4
					write(*,"(' Note: Now the wavefunction is recognized as an unrestricted multiconfiguration wavefunction')")
				else if (all(MOocc==nint(MOocc)).and.all(MOtype/=0)) then !Integer occupation and either A or B
					wfntype=1
					write(*,"(' Note: Now the wavefunction is recognized as an unrestricted single-determinant wavefunction')")
				else if (all(MOocc==nint(MOocc)).and.any(MOtype==0).and.all(MOtype/=2)) then !Integer occupation and at least one orbital is A, and B is unexisted
					wfntype=2
					write(*,"(' Note: Now the wavefunction is recognized as a restricted open-shell wavefunction')")
				else
					write(*,"(' Warning: The type of present wavefunction cannot be identified! You need to reset orbital types')")
					write(*,*) "Press ENTER button to continue"
					read(*,*)
				end if
			end if
			imodwfn=1
		end do
		
	else if (isel==28) then
		do while(.true.)
			write(*,*) "Select the orbitals for which the energy are needed to be changed"
			write(*,*) "e.g. 2,4,13-16,20 means selecting orbitals 2,4,13,14,15,16,20"
			write(*,*) "Input 0 can select all orbitals, input q or 00 can return"
			read(*,"(a)") c1000tmp
			if (c1000tmp(1:1)=='q'.or.c1000tmp(1:2)=='00') exit
			if (c1000tmp(1:1)=='0') then
				numorbsel=nmo
				do i=1,nmo
					orbarr(i)=i
				end do
			else
				call str2arr(c1000tmp,numorbsel,orbarr)
				if ( any(orbarr(1:numorbsel)<1).or.any(orbarr(1:numorbsel)>nmo) ) then
					write(*,*) "Error: One or more orbital indices exceeded valid range!"
					cycle
				end if
			end if
			write(*,*) "0 Recover their initial orbital energies"
			write(*,*) "1 Set the orbital energies to a specific value"
			write(*,*) "2 Add a value to the orbital energies"
			write(*,*) "3 Minus a value from the orbital energies"
			write(*,*) "4 Multiply the orbital energies by a value"
			write(*,*) "5 Divide the orbital energies by a value"
			read(*,*) iselop
			if (iselop==0) then
				MOene(orbarr(1:numorbsel))=MOene_org(orbarr(1:numorbsel))
			else if (iselop==1) then
				write(*,*) "Input the energy in eV, e.g. -3.6"
				read(*,*) tmpval
				MOene(orbarr(1:numorbsel))=tmpval/au2eV
			else if (iselop==2) then
				write(*,*) "Input the value in eV, e.g. 0.6"
				read(*,*) tmpval
				MOene(orbarr(1:numorbsel))=MOene(orbarr(1:numorbsel))+tmpval/au2eV
			else if (iselop==3) then
				write(*,*) "Input the value in eV, e.g. 0.6"
				read(*,*) tmpval
				MOene(orbarr(1:numorbsel))=MOene(orbarr(1:numorbsel))-tmpval/au2eV
			else if (iselop==4) then
				write(*,*) "Input the factor, e.g. 0.9"
				read(*,*) tmpval
				MOene(orbarr(1:numorbsel))=MOene(orbarr(1:numorbsel))*tmpval
			else if (iselop==5) then
				write(*,*) "Input the factor, e.g. 1.2"
				read(*,*) tmpval
				MOene(orbarr(1:numorbsel))=MOene(orbarr(1:numorbsel))/tmpval
			end if
			write(*,*) "Done!"
		end do
	
	else if (isel==29) then
		write(*,"(a)") " Input indices of two orbitals, e.g. 3,8, then all of their information (energy, occupation, coefficients) will be exchanged"
        read(*,*) iorb,jorb
        if (MOtype(iorb)==MOtype(jorb)) then
			enetmp=MOene(iorb)
			MOene(iorb)=MOene(jorb)
			MOene(jorb)=enetmp
			occtmp=MOocc(iorb)
			MOocc(iorb)=MOocc(jorb)
			MOocc(jorb)=enetmp
			allocate(tmparr(nprims))
			tmparr(:)=CO(iorb,:)
			CO(iorb,:)=CO(jorb,:)
			CO(jorb,:)=tmparr(:)
			deallocate(tmparr)
			if (allocated(CObasa)) then
				allocate(tmparr(nbasis))
				if (iorb<=nbasis) then
					tmparr(:)=CObasa(:,iorb)
					CObasa(:,iorb)=CObasa(:,jorb)
					CObasa(:,jorb)=tmparr(:)
                else
					tmparr(:)=CObasb(:,iorb-nbasis)
					CObasb(:,iorb-nbasis)=CObasb(:,jorb-nbasis)
					CObasb(:,jorb-nbasis)=tmparr(:)
                end if
				deallocate(tmparr)
			end if
			imodwfn=1
			write(*,*) "Done!"
        else
			write(*,*) "Error: In order to alter information, the two orbitals must have the same type!"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
        end if
	
	else if (isel==30) then
		write(*,*) "1 Exchange orbital energies (in eV) with occupation numbers"
		write(*,*) "2 Exchange orbital energies (in Hartree) with occupation numbers"
        write(*,*) "Note: For EDDB user, you should choose 1"
        read(*,*) isel2
		allocate(tmparr(nmo))
        tmparr=MOocc
        if (isel2==1) then
			MOocc=MOene*au2eV
        else if (isel2==2) then
			MOocc=MOene
        end if
        MOene=tmparr
        deallocate(tmparr)
        write(*,*) "Orbital energies and occupation numbers have been exchanged with each other"
		if (any(MOocc/=int(MOocc))) then
			if (wfntype==0) then
				wfntype=3 !RHF-> Restricted multiconfiguration wavefunction
				write(*,"(a)") " Note: Now the wavefunction is recognized as a restricted multiconfiguration wavefunction"
			else if (wfntype==1.or.wfntype==2) then !UHF/ROHF-> Unrestricted multiconfiguration wavefunction
				wfntype=4
				write(*,"(a)") " Note: Now the wavefunction is recognized as an unrestricted multiconfiguration wavefunction"
			end if
		end if
        imodwfn=1
        call updatenelec
    
	else if (isel==31) then
		write(*,*) "Input X,Y,Z of translation vector (e.g. 3.2,1.0,0)"
		read(*,*) pbctransx,pbctransy,pbctransz
		write(*,*) "You inputted coordinates are in which unit?  1: Bohr  2: Angstrom"
		read(*,*) iunit
		if (iunit==2) then
			pbctransx=pbctransx/b2a
			pbctransy=pbctransy/b2a
			pbctransz=pbctransz/b2a
		end if
		do i=1,ncenter
			a(i)%x=a(i)%x+pbctransx
			a(i)%y=a(i)%y+pbctransy
			a(i)%z=a(i)%z+pbctransz
		end do
		imodwfn=1
	
	else if (isel==32) then
		write(*,*) "Input X,Y,Z of translation vector (e.g. 3.2,1.0,0)"
		read(*,*) pbctransx,pbctransy,pbctransz
		write(*,*) "You inputted coordinates are in which unit?  1: Bohr  2: Angstrom"
		read(*,*) iunit
		if (iunit==2) then
			pbctransx=pbctransx/b2a
			pbctransy=pbctransy/b2a
			pbctransz=pbctransz/b2a
		end if
		write(*,*) "Duplicate system how many times? e.g. 3"
		read(*,*) numdup
		!_tmp is for backing up current information
        if (allocated(a_tmp)) deallocate(a_tmp)
		allocate(a_tmp(ncenter))
		allocate(b_tmp(nprims))
		allocate(CO_tmp(nmo,nprims))
		a_tmp=a
		b_tmp=b
		CO_tmp=CO
		deallocate(a,b,CO)
		nprims_tmp=nprims
		ncenter_tmp=ncenter
		nprims=nprims*(numdup+1)
		ncenter=ncenter*(numdup+1)
		nelec=nelec*(numdup+1)
		naelec=naelec*(numdup+1)
		nbelec=nbelec*(numdup+1)
		allocate(a(ncenter))
		allocate(b(nprims))
		allocate(CO(nmo,nprims))
		do idup=0,numdup
			a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))=a_tmp(1:ncenter_tmp)
			a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))%x=a_tmp(1:ncenter_tmp)%x+pbctransx*idup
			a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))%y=a_tmp(1:ncenter_tmp)%y+pbctransy*idup
			a(ncenter_tmp*idup+1:ncenter_tmp*(idup+1))%z=a_tmp(1:ncenter_tmp)%z+pbctransz*idup
			b(nprims_tmp*idup+1:nprims_tmp*(idup+1))=b_tmp(1:nprims_tmp)
			b(nprims_tmp*idup+1:nprims_tmp*(idup+1))%center=b_tmp(1:nprims_tmp)%center+ncenter_tmp*idup
			CO(:,nprims_tmp*idup+1:nprims_tmp*(idup+1))=CO_tmp(:,1:nprims_tmp) !Notice that the orbitals do not satisify normalization condition any more, and the orbital occupation number will be artifical
		end do
		deallocate(a_tmp,b_tmp,CO_tmp)
		imodwfn=1
	
	else if (isel==33) then
		write(*,*) "Rotate which orbital? (Input 0 to rotate all orbitals)"
		read(*,*) iorb
		if (iorb/=0) then
			call orbcoeffrotate(iorb)
		else if (iorb==0) then
            write(*,*) "Please wait..."
			do imo=1,nmo
				call orbcoeffrotate(imo)
			end do
			write(*,*) "Also rotate atomic coordinates? (y/n)"
			read(*,*) selectyn
			if (selectyn=='y'.or.selectyn=='Y') then
				do iatm=1,ncenter
					tmpval=a(iatm)%x
					a(iatm)%x=a(iatm)%z
					a(iatm)%z=a(iatm)%y
					a(iatm)%y=tmpval
				end do
			end if
		end if
		write(*,*) "Done!"
	
	else if (isel==34) then
		call getninnerele(ninnerele,1)
		nelec=nelec-ninnerele
		naelec=naelec-ninnerele/2
		nbelec=nbelec-ninnerele/2
		if (wfntype==1.or.wfntype==4) then !UHF and U-post-HF wfn
			MOocc(1:ninnerele/2)=0D0
			do j=1,nmo !Where the first beta orbital appear now
				if (motype(j)==2) exit
			end do
			MOocc(j:j+ninnerele/2-1)=0D0
			write(*,"(' The occupation of',i7,' lowest energy orbitals have been set to zero')") ninnerele
		else if (wfntype==0.or.wfntype==2.or.wfntype==3) then !restricted(=0) or RO(=2) or post-R(=3) wavefunction
			MOocc(1:ninnerele/2)=0D0
			write(*,"(' The occupation of',i7,' lowest energy orbitals have been set to zero')") ninnerele/2
		end if
		if (wfntype==3.or.wfntype==4) write(*,"(' Warning: Discarding inner orbitals for multiconfiguration wavefunction will lead to unexpected result!')") 
		imodwfn=1
	
	else if (isel==35) then
		call selMO_IRREP
	
	else if (isel==36) then
		write(*,*) "Input index of the orbitals, e.g. 2,3,7-10"
		read(*,"(a)") c1000tmp
		call str2arr(c1000tmp,ntmp,orbarr)
		do idxtmp=1,ntmp
			idx=orbarr(idxtmp)
			CO(idx,:)=-CO(idx,:)
			if (allocated(CObasa)) then
				if (idx<=nbasis) then
					CObasa(:,idx)=-CObasa(:,idx)
				else
					CObasb(:,idx-nbasis)=-CObasb(:,idx-nbasis)
				end if
			end if
		end do
		write(*,*) "Done!"
		imodwfn=1
	
	else if (isel==37) then
		if (wfntype==1.or.wfntype==4) then
			write(*,*) "Error: The current wavefunction is already unrestricted!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            cycle
        end if
        !Back up current wavefunction
		if (allocated(CO_back)) deallocate(CO_back,MOene_back,MOocc_back)
		allocate(CO_back(nmo,nprims),MOene_back(nmo),MOocc_back(nmo))
		nmo_back=nmo
		CO_back=CO
		MOene_back=MOene
		MOocc_back=MOocc
        !Reallocate arrays with new size
        deallocate(CO,MOene,MOocc,MOtype)
        nmo=nmo_back*2
        allocate(CO(nmo,nprims),MOene(nmo),MOocc(nmo),MOtype(nmo))
        if (wfntype==0.or.wfntype==3) then !RHF, R-post-HF
			do imo=1,nmo_back
				MOene(imo)=MOene_back(imo)
				MOene(imo+nmo_back)=MOene_back(imo)
				MOocc(imo)=MOocc_back(imo)/2
				MOocc(imo+nmo_back)=MOocc_back(imo)/2
                MOtype(imo)=1
                MOtype(imo+nmo_back)=2
                CO(imo,:)=CO_back(imo,:)
                CO(imo+nmo_back,:)=CO_back(imo,:)
            end do
            if (wfntype==0) wfntype=1
            if (wfntype==3) wfntype=4
        else if (wfntype==2) then !ROHF
			MOocc=0
			do imo=1,nmo_back
				MOene(imo)=MOene_back(imo)
				MOene(imo+nmo_back)=MOene_back(imo)
                if (MOocc_back(imo)==2) then
					MOocc(imo)=1
					MOocc(imo+nmo_back)=1
                else if (MOocc_back(imo)==1) then
					MOocc(imo)=1
                end if
                MOtype(imo)=1
                MOtype(imo+nmo_back)=2
                CO(imo,:)=CO_back(imo,:)
                CO(imo+nmo_back,:)=CO_back(imo,:)
            end do
            wfntype=1
        end if
        deallocate(CO_back,MOene_back,MOocc_back)
        if (allocated(CObasa)) then
			allocate(CObasb(nbasis,nbasis))
            CObasb=CObasa
        end if
        imodwfn=1
        write(*,*) "Finished!"
	
	else if (isel==38) then
		call make_occ_integer_Aufbau
        imodwfn=1
	end if
    
end do
end subroutine



!!----------- Modify orbital occupancy
subroutine modorbocc
use defvar
use util
implicit real*8 (a-h,o-z)
character c1000tmp*1000
integer orbarr(nmo)

do while(.true.)
	write(*,*)
	write(*,*) "Select the orbitals for which the occupation numbers are needed to be changed"
	write(*,*) "e.g. 2,4,13-16,20 means selecting orbitals 2,4,13,14,15,16,20"
	write(*,*) "Input 0 can select all orbitals, input q or 00 can return"
	read(*,"(a)") c1000tmp
	if (c1000tmp(1:1)=='q'.or.c1000tmp(1:2)=='00') exit
	if (c1000tmp(1:1)=='0') then
		numorbsel=nmo
		do i=1,nmo
			orbarr(i)=i
		end do
	else
		call str2arr(c1000tmp,numorbsel,orbarr)
		if ( any(orbarr(1:numorbsel)<1).or.any(orbarr(1:numorbsel)>nmo) ) then
			write(*,*) "Error: One or more orbital indices exceeded valid range!"
			cycle
		end if
	end if
    write(*,*)
	write(*,*) "Set occupation number to which value? e.g. 1.2"
	write(*,*) "Note:"
	write(*,"(a)") " You can also input for example ""+1.1"" ""-1.1"" ""*1.1"" ""/1.1"" to add, minus, multiply and divide the occupation numbers by 1.1"
	write(*,"(a)") " To recover the initial occupation numbers, input ""i"""
	write(*,"(a)") " To generate occupancy state for calculating odd electron density, input ""odd"""
	write(*,"(a)") " To generate occupancy state for calculating FOD, input ""fod"""
	read(*,"(a)") c1000tmp
	if (index(c1000tmp,"odd")/=0) then
		if (wfntype==3) then
			oddnum=0
			do iorb=1,numorbsel
				MOocc(orbarr(iorb))=min(2-MOocc(orbarr(iorb)),MOocc(orbarr(iorb)))
				oddnum=oddnum+MOocc(orbarr(iorb))
			end do
			write(*,*) "Done!"
			write(*,"(a,f12.6)") " Sum of occupation numbers of selected orbitals:",oddnum
			if (nEDFprims/=0) then
				deallocate(b_EDF,CO_EDF,nEDFelecatm)
				nEDFprims=0
				nEDFelec=0
				write(*,"(/,a)") " NOTE: EDF information has been removed to avoid their unexpected influence on subsequent calculation of electron density (corresponding to odd electron density in the present context)"
			end if
        else
			write(*,"(a)") " Error: This function is usable only when the wavefunction is represented in terms of spatial (spinless) natural orbitals"
            write(*,*) "Press ENTER button to return"
            read(*,*)
        end if
	else if (index(c1000tmp,"fod")/=0) then
        fodnum=0
		if (wfntype==3) then
			idxHOMO=nint(nelec/2)
			do iorbidx=1,numorbsel
				iorb=orbarr(iorbidx)
                if (iorb<=idxHOMO) MOocc(iorb)=2-MOocc(iorb)
				fodnum=fodnum+MOocc(iorb)
			end do
        else
			do itmp=nmo,1,-1 !Find the last alpha MO
				if (MOtype(itmp)==1) exit
			end do
			do iorbidx=1,numorbsel
				iorb=orbarr(iorbidx)
                if (iorb<=itmp) then !Alpha
					if (iorb<=nint(naelec)) MOocc(iorb)=1-MOocc(iorb)
                else !Beta
					if (iorb-itmp<=nint(nbelec)) MOocc(iorb)=1-MOocc(iorb)
                end if
				fodnum=fodnum+MOocc(iorb)
			end do
        end if
        write(*,*) "Done!"
		if (nEDFprims/=0) then
			deallocate(b_EDF,CO_EDF,nEDFelecatm)
			nEDFprims=0
			nEDFelec=0
			write(*,"(/,a)") " NOTE: EDF information has been removed to avoid their unexpected influence on subsequent calculation of electron density (corresponding to FOD in the present context)"
		end if
        if (numorbsel==nmo) then
			write(*,"(a,f12.6)") " N_FOD is",fodnum
        else
			write(*,"(a,f12.6)") " Sum of occupation numbers of selected orbitals:",fodnum
        end if
	else if (index(c1000tmp,"i")/=0) then
		MOocc(orbarr(1:numorbsel))=MOocc_org(orbarr(1:numorbsel))
		write(*,*) "The occupation numbers have been recovered"
	else if (c1000tmp(1:1)=='+'.or.c1000tmp(1:1)=='-'.or.c1000tmp(1:1)=='*'.or.c1000tmp(1:1)=='/') then
		read(c1000tmp(2:),*) tmpval
		if (c1000tmp(1:1)=='+') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))+tmpval
		if (c1000tmp(1:1)=='-') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))-tmpval
		if (c1000tmp(1:1)=='*') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))*tmpval
		if (c1000tmp(1:1)=='/') MOocc(orbarr(1:numorbsel))=MOocc(orbarr(1:numorbsel))/tmpval
		write(*,*) "Done!"
	else
		read(c1000tmp,*) tmpval
		MOocc(orbarr(1:numorbsel))=tmpval
		write(*,*) "Done!"
	end if
	call updatenelec !Update the number of electrons
	imodwfn=1
	if (any(MOocc/=int(MOocc))) then
		if (wfntype==0.or.wfntype==2) then
            if (wfntype==2) MOtype=0
			wfntype=3 !RHF/ROHF-> Restricted multiconfiguration wavefunction
			write(*,"(a)") " Note: Now the wavefunction is recognized as a restricted multiconfiguration wavefunction"
		else if (wfntype==1) then !UHF-> Unrestricted multiconfiguration wavefunction
			wfntype=4
			write(*,"(a)") " Note: Now the wavefunction is recognized as an unrestricted multiconfiguration wavefunction"
		end if
	end if
end do
end subroutine




!!-------- Make orbital occupations integer and satisfy Aufbau principle
subroutine make_occ_integer_Aufbau
use defvar
implicit real*8 (a-h,o-z)
MOocc=0
if (wfntype==0.or.wfntype==3) then !Closed-shell or RO
	diffint=abs(nelec-nint(nelec)) !Due to recording accuracy, the difference is not exactly zero even total number of electron actually is integer
    if (diffint<1D-5) then
		write(*,"(a)") " Note: The difference between number of electrons and integer is less than 1E-5, so number of electrons is set to integer"
        nelec=nint(nelec)
    end if
    nintocc=floor(nelec/2D0)
	MOocc(1:nintocc)=2
    write(*,*) "Done!"
    if (diffint<1D-5.and.wfntype==3) then
		wfntype=0
		write(*,*) "Note: The wavefunction has been set to single-determinant closed-shell type"
    else !Distribute reminder electrons to LUMO
		MOocc(nintocc+1)=mod(nelec,2D0)
		if (nmo>=nintocc+2) MOocc(nintocc+2:)=0
    end if
else if (wfntype==1.or.wfntype==4) then !Unrestricted
	diffinta=abs(naelec-nint(naelec))
	diffintb=abs(nbelec-nint(nbelec))
    if (diffinta<1D-5.and.diffintb<1D-5) then
		write(*,"(a)") " Note: The difference between number of electrons and integer is less than 1E-5, so number of electrons is set to integer"
        naelec=nint(naelec)
        nbelec=nint(nbelec)
        nelec=nint(nelec)
    end if
    nintocca=floor(naelec)
    nintoccb=floor(nbelec)
    do ibeta=1,nmo
		if (MOtype(ibeta)==2) exit
    end do
	MOocc(1:nintocca)=1
	MOocc(ibeta:ibeta+nintoccb-1)=1
    write(*,*) "Done!"
    if (diffinta<1D-5.and.diffintb<1D-5.and.wfntype==4) then
		wfntype=1
		write(*,"(a)") " Note: The wavefunction has been set to single-determinant unrestricted open-shell type"
    else !Distribute reminder electrons to LUMO
		MOocc(nintocca+1)=mod(naelec,1D0)
		if (ibeta-1>=nintocca+2) MOocc(nintocca+2:ibeta-1)=0
		MOocc(ibeta+nintoccb)=mod(nbelec,1D0)
		if (nmo>=ibeta+nintoccb+1) MOocc(ibeta+nintoccb+1:)=0
    end if
end if
end subroutine




!!---------------- Select MOs according to irreducible representation
subroutine selMO_IRREP
use defvar
use util
character symlab(nmo)*4,c2000tmp*2000,symstat(nmo)*9 !Allocate the array lengths as upper limit
integer tmparr(nmo),symNorb(nmo) !Allocate the array lengths as upper limit

if (wfntype/=0.and.wfntype/=1) then
	write(*,"(a)") " Error: This function only works for RHF or UHF wavefunctions (or the DFT counterparts)"
	return
end if
nsym=0
do imo=1,nmo
	if (MOocc_org(imo)==0D0) cycle
	if (all(symlab(1:nsym)/=MOsym(imo))) then
		nsym=nsym+1
		symlab(nsym)=MOsym(imo)
	end if
end do
symNorb=0
do imo=1,nmo
	if (MOocc_org(imo)==0D0) cycle
	do isym=1,nsym
		if (MOsym(imo)==symlab(isym)) symNorb(isym)=symNorb(isym)+1
	end do
end do
symstat="Normal"
MOocc=MOocc_org
if (imodwfn==1) write(*,*) "Note: Original occupation status has been recovered"
write(*,*) "Note: Only the orbitals that originally occupied are taken into account here"
do while(.true.)
	write(*,*)
	write(*,*) "Information of various irreducible representations:"
	do isym=1,nsym
		write(*,"(i5,'  Sym: ',a,'  N_orb:',i5,'    Status: ',a)") isym,symlab(isym),symNorb(isym),symstat(isym)
	end do
	write(*,*)
	write(*,*) "0 Save and return"
	write(*,*) "1 Discard specific irreducible representations"
	write(*,*) "2 Recover original status"
	write(*,*) "3 Reverse status"
	read(*,*) isel
	
	if (isel==0) then
		call updatenelec
		imodwfn=1
		write(*,*) "The current orbital occupation status has been saved"
		write(*,*) "Updating density matrix..."
		call genP
		write(*,*) "Density matrix has been updated"
		exit
	else if (isel==2) then
		MOocc=MOocc_org
		symstat="Normal"
	else if (isel==1.or.isel==3) then
		if (isel==1) then
			write(*,*) "Input the index of the irreducible representations to be discarded, e.g. 1,3-5"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,nsymsel,tmparr)
			do isym=1,nsymsel
				symstat(tmparr(isym))="Discarded"
			end do
		else if (isel==3) then
			do isym=1,nsym
				if (symstat(isym)=="Normal") then
					symstat(isym)="Discarded"
				else
					symstat(isym)="Normal"
				end if
			end do
		end if
		do imo=1,nmo
			if (MOocc_org(imo)==0D0) cycle
			do isym=1,nsym
				if (MOsym(imo)==symlab(isym)) then
					if (symstat(isym)=="Normal") MOocc(imo)=MOocc_org(imo)
					if (symstat(isym)=="Discarded") MOocc(imo)=0D0
					exit
				end if
			end do
		end do
	end if
	write(*,*) "Done!"
end do
end subroutine



!!---------- Update the number of electrons
subroutine updatenelec
use defvar
integer imo
nelec=0
naelec=0
nbelec=0
do imo=1,nmo
	if (MOtype(imo)==0) then
		naelec=naelec+MOocc(imo)/2D0
		nbelec=nbelec+MOocc(imo)/2D0
	else if (MOtype(imo)==1) then
		naelec=naelec+MOocc(imo)
	else if (MOtype(imo)==2) then
		nbelec=nbelec+MOocc(imo)
	end if
end do
nelec=naelec+nbelec
end subroutine



!!-------- Generate nelec, naelec, nbelec by guessing when wavefunction information is not available
!Commonly invoked when reading file only contains geometry information
!The system is always regarded as neutral!!! If number of electrons is even, then naelec=nbelec=nelec/2. Else odd, naelec=nbelec+1
subroutine guessnelec
use defvar
nelec=sum(a%charge)
if (mod(nint(nelec),2)==0) then
    naelec=nint(nelec)/2
    nbelec=naelec
else
    nbelec=(nint(nelec)-1)/2
    naelec=nbelec+1
end if
end subroutine
			
			

!!-------- Check if present wavefunction is sanity, i.e. all orbital satisfies normalization condition
subroutine wfnsanity
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 GTFSmat(nprims*(nprims+1)/2)
real*8,allocatable :: tmpmat(:,:)

itype=1
if (allocated(CObasa)) then
	write(*,*) "1 Check normalization based on GTF information"
	write(*,*) "2 Check normalization based on basis function information"
    write(*,*) "3 Check Tr(P*S) (should be equal to number of electrons)"
	read(*,*) itype
end if

rmaxdev=0
rmaxdevint=0
if (itype==1) then
	if (ifPBC/=0) then
		write(*,*) "Error: This function is only available for non-periodic case"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
    call genGTFSmat(GTFSmat,nprims*(nprims+1)/2)
    do imo=1,nmo
	    tmp=0
	    !$OMP parallel shared(tmp) private(iGTF,jGTF,tmp1,tmp2,tmpprivate) num_threads(nthreads)
	    tmpprivate=0
	    !$OMP do schedule(dynamic)
	    do iGTF=1,nprims
			tmp1=CO(imo,iGTF)**2*GTFSmat(iGTF*(iGTF-1)/2+iGTF)
            tmp2=0
		    do jGTF=iGTF+1,nprims
			    tmp2=tmp2+CO(imo,jGTF)*GTFSmat(jGTF*(jGTF-1)/2+iGTF)
		    end do
		    tmpprivate=tmpprivate+tmp1+2*CO(imo,iGTF)*tmp2
	    end do
	    !$OMP END DO
	    !$OMP CRITICAL
	    tmp=tmp+tmpprivate
	    !$OMP END CRITICAL
	    !$OMP END PARALLEL
	    write(*,"(' Orbital',i7,'   Occ:',f10.6,'   Value:',f16.10)") imo,MOocc(imo),tmp
	    tmpt=abs(tmp-1)
	    if (tmpt>rmaxdev) rmaxdev=tmpt
	    tmpt=abs(tmp-nint(tmp))
	    if (tmpt>rmaxdevint) rmaxdevint=tmpt
    end do
	write(*,"(' Maximum deviation to 1:',f16.10)") rmaxdev
	write(*,"(' Maximum deviation to integer:',f16.10)") rmaxdevint

else if (itype==2) then
    call ask_Sbas_PBC
    allocate(tmpmat(nbasis,nbasis))
    if (allocated(CObasb)) write(*,*) "Note: The check is only conducted to alpha orbitals"
    do imo=1,nbasis
        if (all(CObasa(:,imo)==0)) cycle
        tmp=0
	    !$OMP parallel shared(tmp) private(ibas,jbas,tmp1,tmp2,tmpprivate) num_threads(nthreads)
	    tmpprivate=0
	    !$OMP do schedule(dynamic)
        do ibas=1,nbasis
			tmp1=CObasa(ibas,imo)**2*Sbas(ibas,ibas)
            tmp2=0
            do jbas=ibas+1,nbasis
                tmp2=tmp2+CObasa(jbas,imo)*Sbas(ibas,jbas)
            end do
            tmpprivate=tmpprivate+tmp1+2*CObasa(ibas,imo)*tmp2
        end do
	    !$OMP END DO
	    !$OMP CRITICAL
	    tmp=tmp+tmpprivate
	    !$OMP END CRITICAL
	    !$OMP END PARALLEL
	    write(*,"(' Orbital',i7,'   Occ:',f10.6,'   Value:',f16.10)") imo,MOocc(imo),tmp
	    tmpt=abs(tmp-1)
	    if (tmpt>rmaxdev) rmaxdev=tmpt
	    tmpt=abs(tmp-nint(tmp))
	    if (tmpt>rmaxdevint) rmaxdevint=tmpt
    end do
	write(*,"(' Maximum deviation to 1:',f16.10)") rmaxdev
	write(*,"(' Maximum deviation to integer:',f16.10)") rmaxdevint
else if (itype==3) then
    call ask_Sbas_PBC
	write(*,"(' Trace of density matrix multiplied by overlap matrix:',f14.6)") sum(Ptot*Sbas)
    write(*,"(' Expected number of electrons:',i10)") nint(nelec)
    if (wfntype==1) then
		write(*,"(' Trace of alpha density matrix multiplied by overlap matrix:',f14.6)") sum(Palpha*Sbas)
		write(*,"(' Expected number of alpha electrons:',i10)") nint(naelec)
		write(*,"(' Trace of beta density matrix multiplied by overlap matrix: ',f14.6)") sum(Pbeta*Sbas)
		write(*,"(' Expected number of beta electrons: ',i10)") nint(nbelec)
    end if
end if
write(*,*)
write(*,*) "Press ENTER button to continue"
read(*,*)
end subroutine




!!---------- Return normalization coefficient for specific type of Cartesian type GTF, see Levine 5ed p487
!The meaning of itype is defined in GTFtype2name
real*8 function normgau(itype,exp)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 exp
ix=type2ix(itype)
iy=type2iy(itype)
iz=type2iz(itype)
normgau=(2*exp/pi)**0.75D0*dsqrt( (8*exp)**(ix+iy+iz)*ft(ix)*ft(iy)*ft(iz)/(ft(2*ix)*ft(2*iy)*ft(2*iz)) )
end function


!---- Renormalizing shells (modifying contraction coefficients in basis shells) for Molden input file
!Various basis functions in a shell may have different normalization factor, however, since we are dealing with shell, 
!we can use any specific basis function to derive normalization factor and fix contraction coefficients
!Here we only consider (Lval,0,0) or (0,Lval,0) or (0,0,Lval) type, which have identical normalization factor
!Lval=0/1/2/3/4/5 corresponds to s/p/d/f/g/h, exp and con are exponents and contraction coefficients in the basis shell, respectively
subroutine renormmoldengau(nlen,Lval,exp,con)
implicit real*8 (a-h,o-z)
integer nlen,Lval
real*8 exp(nlen),con(nlen),ctmp(nlen)

pi=acos(-1D0)
call genn1n2nf(Lval,n1,n2,nf)
fc=2D0**n1/(pi**3*nf)
do i=1,nlen
	prmnormfac=sqrt(sqrt(fc*exp(i)**n2)) !Normalization factor of (Lval,0,0) type of GTF
	ctmp(i)=con(i)*prmnormfac
end do
facnorm=0D0 !Calculate <Lval,0,0|Lval,0,0> overlap integral for (Lval,0,0) type of contracted basis function
do i=1,nlen
	do j=1,i
	  expavg=(exp(i)+exp(j))/2D0
	  facadd=ctmp(i)*ctmp(j)/sqrt(fc*expavg**n2)
	  if (i/=j) facadd=facadd*2
	  facnorm=facnorm+facadd
	end do
end do
if (facnorm>1D-10) facnorm=1/sqrt(facnorm)
con=con*facnorm !Fix contraction coefficients
end subroutine
!!---- Produce normalization factor of (Lval,0,0) type of GTF with angular moment of Lval, &
!used to fix the contraction coefficient problem of the Molden input file generated by ORCA
!This routine can also be replaced by the slower routine "normgau", e.g. renormgau_ORCA(0.8D0,3) = normgau(21,0.8D0)
!cf. fnorm_lmn of m2a
real*8 function renormgau_ORCA(exp,Lval)
real*8 exp,pi
integer Lval,n1,n2,nf
pi=acos(-1D0)
call genn1n2nf(Lval,n1,n2,nf)
renormgau_ORCA=dsqrt(dsqrt(2**n1*exp**n2/(pi**3*nf))) !Get norm, the norm^4 = 2^n1 * a^n2 / (pi^3 * nf)
end function
!!--- Generate n1,n2,f for a given angular moment (up to h), the data only correspond to (Lval,0,0)=(0,Lval,0)=(0,0,Lval) case
!n1=3+4*(l+m+n)
!n2=3+2*(l+m+n)
!nf=[(2l-1)!!(2m-1)!!(2n-1)!!]^2  <--- e.g. f shell, putting XXXX, i.e. (4,0,0) into it returns 11025. Here 7!! means 7*5*3*1
subroutine genn1n2nf(Lval,n1,n2,nf)
integer Lval,n1,n2,nf
if (Lval==0) then
	n1=3
	n2=3
	nf=1
else if (Lval==1) then
	n1=7
	n2=5
	nf=1
else if (Lval==2) then
	n1=11
	n2=7
	nf=9
else if (Lval==3) then
	n1=15
	n2=9
	nf=225
else if (Lval==4) then
	n1=19
	n2=11
	nf=11025
else if (Lval==5) then
	n1=23
	n2=13
	nf=893025
end if
end subroutine



!!----- Use Lowdin orthogonalization method to transform density matrix and coefficient matrix &
!to orthonormal basis, meanwhile update Sbas to identity matrix. If itask=1, also update Fock matrix to orthonormal basis if available
!See Szabo p143 for details
!NOTICE that the matrix Xmat in this subroutine (Xmat=S^(1/2)) is inverse of that of p143 of Szabo book
subroutine symmortho(itask)
use defvar
use util
implicit real*8 (a-h,o-z)
integer itask
real*8 Umat(nbasis,nbasis),svalvec(nbasis),Xmat(nbasis,nbasis),tmpmat(nbasis,nbasis)
real*8,allocatable :: Xmatinv(:,:)

!call walltime(iwalltime1)
call diagsymat(Sbas,Umat,svalvec,ierror)
if (ierror/=0) write(*,*) "Error: Diagonalization of overlap matrix failed!"
!call walltime(iwalltime2)
!write(*,"(' Time cost for diagonalizing overlap matrix',i10,' s')") iwalltime2-iwalltime1

!Now Sbas is already diagonalized
forall (i=1:nbasis) Sbas(i,i)=dsqrt(svalvec(i)) !Use Sbas as temporary matrix here

!call walltime(iwalltime1)
tmpmat=matmul_blas(Umat,Sbas,nbasis,nbasis)
Xmat=matmul_blas(tmpmat,Umat,nbasis,nbasis,0,1)
!call walltime(iwalltime2)
!write(*,"(' Time cost for constructing X matrix',i10,' s')") iwalltime2-iwalltime1

!call walltime(iwalltime1)
tmpmat=matmul_blas(Xmat,Ptot,nbasis,nbasis)
Ptot=matmul_blas(tmpmat,Xmat,nbasis,nbasis)
if (allocated(Palpha)) then
	tmpmat=matmul_blas(Xmat,Palpha,nbasis,nbasis)
	Palpha=matmul_blas(tmpmat,Xmat,nbasis,nbasis)
	Pbeta=Ptot-Palpha
end if
!call walltime(iwalltime2)
!write(*,"(' Time cost for transforming density matrix',i10,' s')") iwalltime2-iwalltime1

!call walltime(iwalltime1)
CObasa=matmul_blas(Xmat,CObasa,nbasis,nbasis)
if (allocated(CObasb)) CObasb=matmul_blas(Xmat,CObasb,nbasis,nbasis)
!call walltime(iwalltime2)
!write(*,"(' Time cost for transforming coefficient matrix',i10,' s')") iwalltime2-iwalltime1

if (itask==1.and.allocated(FmatA)) then
	allocate(Xmatinv(nbasis,nbasis))
	forall (i=1:nbasis) Sbas(i,i)=1D0/Sbas(i,i)
	tmpmat=matmul_blas(Umat,Sbas,nbasis,nbasis) !Szabo, Eq. 3.167
	Xmatinv=matmul_blas(tmpmat,Umat,nbasis,nbasis,0,1)
	tmpmat=matmul_blas(transpose(Xmatinv),FmatA,nbasis,nbasis) !Szabo, Eq. 3.177
	FmatA=matmul_blas(tmpmat,Xmatinv,nbasis,nbasis)
	if (allocated(FmatB)) then
		tmpmat=matmul_blas(transpose(Xmatinv),FmatB,nbasis,nbasis)
		FmatB=matmul_blas(tmpmat,Xmatinv,nbasis,nbasis)
    end if
end if

forall(i=1:nbasis) Sbas(i,i)=1D0 !Set overlap matrix to identity matrix as expected
end subroutine


!!----- Input overlap matrix and return Lowdin orthogonalization transformation matrix Xmat=S^0.5 and Xmatinv=S^-0.5
!Smatin is inputted overlap matrix, which will not be modified
!NOTICE that the matrix Xmat in this subroutine (Xmat=S^(1/2)) is inverse of that of p143 of Szabo book
subroutine symmorthomat(Smatin,Xmat,Xmatinv)
use defvar
use util
real*8 Umat(nbasis,nbasis),svalvec(nbasis),Smatin(nbasis,nbasis),Smat(nbasis,nbasis),Xmat(nbasis,nbasis),Xmatinv(nbasis,nbasis),tmpmat(nbasis,nbasis)
Smat=Smatin
call diagsymat(Smat,Umat,svalvec,ierror)
if (ierror/=0) write(*,*) "Error: Diagonalization of overlap matrix is fail!"
forall (i=1:nbasis) Smat(i,i)=dsqrt(svalvec(i))
!Xmat=matmul(matmul(Umat,Smat),transpose(Umat))
tmpmat=matmul_blas(Umat,Smat,nbasis,nbasis)
Xmat=matmul_blas(tmpmat,Umat,nbasis,nbasis,0,1)
forall (i=1:nbasis) Smat(i,i)=1D0/Smat(i,i)
!Xmatinv=matmul(matmul(Umat,Smat),transpose(Umat))
tmpmat=matmul_blas(Umat,Smat,nbasis,nbasis)
Xmatinv=matmul_blas(tmpmat,Umat,nbasis,nbasis,0,1)
end subroutine




!!!------------ Return the number of inner-core orbitals
subroutine getninnerele(ninnerele,info)
use defvar
integer ninnerele,info
ninnerele=0
do i=1,ncenter
	if (int(a(i)%charge)/=a(i)%index) then
		if (info==1) write(*,"(' Note: Atom',i5,' is not taken into account since it utilizes pseudopotential')") i
		cycle
	end if
	if (a(i)%index>2.and.a(i)%index<=10) ninnerele=ninnerele+2
	if (a(i)%index>10.and.a(i)%index<=18) ninnerele=ninnerele+10
	if (a(i)%index>18.and.a(i)%index<=36) ninnerele=ninnerele+18
	if (a(i)%index>36.and.a(i)%index<=54) ninnerele=ninnerele+36
	if (a(i)%index>54.and.a(i)%index<=86) ninnerele=ninnerele+54
	if (a(i)%index>86) ninnerele=ninnerele+86
end do
end subroutine



!!!------------------------- Swap two GTF
subroutine swapGTF(i,j,swaptype)
use defvar
integer n,i,j
character swaptype*3
type(primtype) tempb !For exchanging basis functions' order
if (swaptype=="all") then
	tempb=b(i)
	b(i)=b(j)
	b(j)=tempb
else if (swaptype=="cen") then
	tempb%center=b(i)%center
	b(i)%center=b(j)%center
	b(j)%center=tempb%center
else if (swaptype=="typ") then
	tempb%type=b(i)%type
	b(i)%type=b(j)%type
	b(j)%type=tempb%type
else if (swaptype=="exp") then
	tempb%exp=b(i)%exp
	b(i)%exp=b(j)%exp
	b(j)%exp=tempb%exp
end if
if (swaptype=="all".or.swaptype=="MO ") then
	do n=1,nmo
		temp=CO(n,i)
		CO(n,i)=CO(n,j)
		CO(n,j)=temp
	end do
end if
end subroutine




!!!---- Rotate(exchange) GTF and basis function coefficients within all shell in different direction of specific orbital
! If using this three times, namely XYZ->ZXY->YZX->XYZ, the original coefficient will be recovered
! In detail, for example, for d-type will lead to such coefficient exchange: XX to YY, YY to ZZ, ZZ to XX, XY to YZ, XZ to XY, YZ to XZ
! exchange only involve the GTFs/basis func. in the same shell
subroutine orbcoeffrotate(orb) !orb=Rotate which orbital (from 1 to nmo)
use defvar
implicit real*8 (a-h,o-z)
integer orb
real*8 COorborg(nprims) !For backing up origin CO
real*8 CObasa_tmp(nbasis) !For backing up origin CObasa
real*8 CObasb_tmp(nbasis) !For backing up origin CObasb
COorborg(:)=CO(orb,:)
do i=1,nprims
	ixtmp=type2iz(b(i)%type)
	iytmp=type2ix(b(i)%type)
	iztmp=type2iy(b(i)%type)
	do j=1,nprims
		if (type2ix(b(j)%type)==ixtmp.and.type2iy(b(j)%type)==iytmp.and.&
		type2iz(b(j)%type)==iztmp.and.b(j)%exp==b(i)%exp.and.b(j)%center==b(i)%center) CO(orb,j)=COorborg(i)
	end do
end do
if (allocated(CObasa)) then
    if (orb<=nbasis) then
	    CObasa_tmp=CObasa(:,orb)
	else
        CObasb_tmp=CObasb(:,orb-nbasis)
    end if
	do iatm=1,ncenter
		if (basstart(iatm)==0) cycle
		do ibas=basstart(iatm),basend(iatm)
			ityp=bastype(ibas)
			ixtmp=type2iz(ityp)
			iytmp=type2ix(ityp)
			iztmp=type2iy(ityp)
			do jbas=basstart(iatm),basend(iatm)
				jtyp=bastype(jbas)
				if (type2ix(jtyp)==ixtmp.and.type2iy(jtyp)==iytmp.and.type2iz(jtyp)==iztmp.and.basshell(ibas)==basshell(jbas)) then
                    if (orb<=nbasis) then
					    CObasa(jbas,orb)=CObasa_tmp(ibas)
					else
                        CObasb(jbas,orb-nbasis)=CObasb_tmp(ibas)
                    end if
				end if
			end do
		end do
	end do
end if
end subroutine



!!!------ Output molecular formula
subroutine showformula
use defvar
implicit real*8 (a-h,o-z)
character tmp*6
write(*,"(' Formula: ')",advance="no")
do i=0,nelesupp
	n=0
	do iatm=1,ncenter
		if (a(iatm)%index==i) n=n+1
	end do
	write(tmp,"(i6)") n
	if (n/=0) write(*,"(a,a,' ')",advance="no") trim(ind2name(i)),trim(adjustl(tmp))
end do
write(*,"('     Total atoms:',i8)") ncenter
end subroutine




!!!----------- Resize number of orbitals of CO, MOene, MOtype, MOocc to "newnmo", also resize number of GTFs of CO to "newnprims"
subroutine resizebynmo(newnmo,newnprims)
use defvar
implicit real*8 (a-h,o-z)
real*8,allocatable :: CO_bk(:,:),MOene_bk(:),MOocc_bk(:)
integer,allocatable :: MOtype_bk(:)
integer newnmo,oldnmo,newnprims,oldnprims
oldnmo=size(CO,1)
oldnprims=size(CO,2)
allocate(CO_bk(oldnmo,oldnprims),MOene_bk(oldnmo),MOocc_bk(oldnmo),MOtype_bk(oldnmo))
CO_bk=CO
MOene_bk=MOene
MOocc_bk=MOocc
MOtype_bk=MOtype
deallocate(CO,MOene,MOocc,MOtype)
allocate(CO(newnmo,newnprims),MOene(newnmo),MOocc(newnmo),MOtype(newnmo))
if (newnmo>=oldnmo) then !Enlarge array size, don't forget to fill the gap afterwards
	if (newnprims>=oldnprims) CO(1:oldnmo,1:oldnprims)=CO_bk(:,:)
	if (newnprims<oldnprims) CO(1:oldnmo,:)=CO_bk(:,1:newnprims)
	MOene(1:oldnmo)=MOene_bk(:)
	MOocc(1:oldnmo)=MOocc_bk(:)
	MOtype(1:oldnmo)=MOtype_bk(:)
else if (newnmo<oldnmo) then !Reduce array size
	if (newnprims>=oldnprims) CO(:,1:oldnprims)=CO_bk(1:newnmo,:)
	if (newnprims<oldnprims) CO(:,:)=CO_bk(1:newnmo,1:newnprims)
	MOene(:)=MOene_bk(1:newnmo)
	MOocc(:)=MOocc_bk(1:newnmo)
	MOtype(:)=MOtype_bk(1:newnmo)
end if
deallocate(CO_bk,MOene_bk,MOocc_bk,MOtype_bk)
end subroutine




!!!------------ Generate gjf of atoms in molecule, and invoke Gaussian to get .wfn, then input them into custom list
subroutine setpromol
use defvar
use util
implicit real*8 (a-h,o-z)
integer :: itype=0
character(len=2) typename(100),nametmp
character basisset*80,tmpdir*80,c80tmp*80,outwfnname*80
logical alivegauout,alivewfntmp,aliveatomwfn
if (isys==1) call delfile("gxx.* fort.6 Gau*.inp") !Clean Gaussian scratch files in current folder

!The only difference between c80tmp and tmpdir is that the latter has \ or / separator at the end
if (iwfntmptype==1) then
	if (isys==1) tmpdir="wfntmp\"
	if (isys==2) tmpdir="wfntmp/"
	c80tmp="wfntmp"
    call inquire_dir("wfntmp",alivewfntmp)
	if (isys==1.and.alivewfntmp) then !Delete old wfntmp folder
		write(*,*) "Running: rmdir /S /Q wfntmp"
		call system("rmdir /S /Q wfntmp")
	else if (isys==2.and.alivewfntmp) then
		write(*,*) "Running: rm -rf wfntmp"
		call system("rm -rf wfntmp")
	end if
else if (iwfntmptype==2) then
	do i=1,9999 !Find a proper name of temporary folder
		write(c80tmp,"('wfntmp',i4.4)") i
		call inquire_dir(c80tmp,alivewfntmp)
		if (.not.alivewfntmp) exit
	end do
	if (isys==1) write(tmpdir,"('wfntmp',i4.4,'\')") i
	if (isys==2) write(tmpdir,"('wfntmp',i4.4,'/')") i
end if
write(*,*) "Running: mkdir "//trim(c80tmp) !Build new temporary folder
call system("mkdir "//trim(c80tmp))
call inquire_dir("atomwfn",aliveatomwfn)
if (isys==1.and.aliveatomwfn) then
	write(*,*) "Running: copy atomwfn\*.wfn "//trim(tmpdir)
	call system("copy atomwfn\*.wfn "//trim(tmpdir))
else if (isys==2.and.aliveatomwfn) then
	write(*,*) "Running: cp atomwfn/*.wfn "//trim(tmpdir)
	call system("cp atomwfn/*.wfn "//trim(tmpdir))
end if

noatmwfn=0 !Check if the atomic wfn file have pre-stored in atomwfn folder, if not, invoke gaussian to calc it
do i=1,nfragatm
	if (isys==1) inquire(file="atomwfn\"//a(fragatm(i))%name//".wfn",exist=alive)
	if (isys==2) inquire(file="atomwfn/"//a(fragatm(i))%name//".wfn",exist=alive)
	if (.not.alive) then
		noatmwfn=1
		exit
	end if
end do

if (noatmwfn==0) then
	write(*,"(a)") " All atom .wfn files needed have already presented in ""atomwfn"" folder, we will not calculate them"
else if (noatmwfn==1) then !Some or all atomic wfn don't exist, calc them
	!Select calculation level
	write(*,"(a)") " Note: Some or all atom .wfn files needed are not present in ""atomwfn"" folder, they must be calculated now. See Section 3.7.3 of the manual for detail."
	write(*,"(a)") " Now please input the level for calculating atom wfn files, theoretical method is optional."
	write(*,"(a)") " For example: B3LYP/6-31G*    You can also add other keywords at the same time, e.g. wB97XD/def2TZVP scf=xqc int=ultrafine"
	read(*,"(a)") basisset !Note: 6d 10f is not required for generating wfn files, since the work has been done in L607 internally
	!Check Gaussian path
	inquire(file=gaupath,exist=alive)
	if (.not.alive) then
		write(*,*) "Could not find Gaussian path defined in ""gaupath"" variable in settings.ini"
		if (isys==1) write(*,*) "Input the path of Gaussian executable file, e.g. ""D:\study\g16w\g16.exe"""
		if (isys==2) write(*,*) "Input the path of Gaussian executable file, e.g. ""/sob/g16/g16"""
		do while(.true.)
			read(*,"(a)") gaupath
			inquire(file=gaupath,exist=alive)
			if (alive) exit
			write(*,*) "Could not find Gaussian executable file, input again"
		end do
	end if
end if

!Generate .gjf file for all elements, regardless if their wfn file have already presented, meanwhile count the total number of elements
itype=0
do i=1,nfragatm
	inquire(file=trim(tmpdir)//a(fragatm(i))%name//".gjf",exist=alive)
	if (.not.alive) then
		itype=itype+1 !How many different types
		typename(itype)=a(fragatm(i))%name
				
		if (a_org(fragatm_org(i))%index>36) then
			inquire(file=trim(tmpdir)//a(fragatm(i))%name//".wfn",exist=alive)
			if (.not.alive) then !The wfn file of the heavy element hasn't been provided in "atomwfn" and hence cannot be found in "wfntmp" here
				write(*,"(a,a,a)") " Error: Multiwfn cannot invoke Gaussian to generate wavefunction file and sphericalize density for ",a(fragatm(i))%name,", since its &
				&index is larger than 36! You should provide corresponding atom .wfn files in ""atomwfn"" folder manually"
				write(*,*) "Press ENTER button to continue"
				read(*,*)
				return
			end if
		end if
		
		open(14,file=trim(tmpdir)//a(fragatm(i))%name//".gjf",status="replace")
		!If user inputted including "/" e.g. B3LYP/6-31g*, will replace default theoretical method
		if (index(basisset,'/')==0) then
			if (a(fragatm(i))%index<=20.or.a(fragatm(i))%index>=31) then
				write(14,"(a,/)") "#T out=wfn ROHF/"//trim(basisset) !Main group elements. If not use scf=sp, in g09, RO calculations for IIIA elements are to converge
				write(14,"(a,/)") "Temporary file for promolecule, ROHF"//trim(basisset)
			else
				write(14,"(a,/)") "#T out=wfn UB3LYP/"//trim(basisset) !Transition metals
				write(14,"(a,/)") "Temporary file for promolecule, UB3LYP"//trim(basisset)
			end if
		else
			if (a(fragatm(i))%index<=20.or.a(fragatm(i))%index>=31) then
				write(14,"(a,/)") "#T out=wfn RO"//trim(basisset) !Main group elements
				write(14,"(a,/)") "Temporary file for promolecule, RO"//trim(basisset)
			else
				write(14,"(a,/)") "#T out=wfn U"//trim(basisset) !Transition metals (RO may leads to convergence problem)
				write(14,"(a,/)") "Temporary file for promolecule, U"//trim(basisset)
			end if
		end if

		!Currently support up to the fourth row
		if (a(fragatm(i))%name=="H ".or.a(fragatm(i))%name=="Li".or.a(fragatm(i))%name=="Na".or.a(fragatm(i))%name=="K") write(14,*) "0 2"
		if (a(fragatm(i))%name=="Be".or.a(fragatm(i))%name=="Mg".or.a(fragatm(i))%name=="Ca") write(14,*) "0 1"
		if (a(fragatm(i))%name=="B ".or.a(fragatm(i))%name=="Al".or.a(fragatm(i))%name=="Ga") write(14,*) "0 2"
		if (a(fragatm(i))%name=="C ".or.a(fragatm(i))%name=="Si".or.a(fragatm(i))%name=="Ge") then
			if (SpherIVgroup==0) write(14,*) "0 5"
			if (SpherIVgroup==1) write(14,*) "0 3"
		end if
		if (a(fragatm(i))%name=="N ".or.a(fragatm(i))%name=="P ".or.a(fragatm(i))%name=="As") write(14,*) "0 4"
		if (a(fragatm(i))%name=="O ".or.a(fragatm(i))%name=="S ".or.a(fragatm(i))%name=="Se") write(14,*) "0 3"
		if (a(fragatm(i))%name=="F ".or.a(fragatm(i))%name=="Cl".or.a(fragatm(i))%name=="Br") write(14,*) "0 2"
		if (a(fragatm(i))%name=="He".or.a(fragatm(i))%name=="Ne".or.a(fragatm(i))%name=="Ar".or.a(fragatm(i))%name=="Kr") write(14,*) "0 1"
		if (a(fragatm(i))%name=="Sc") write(14,*) "0 2" !3d1 4s2
		if (a(fragatm(i))%name=="Ti") write(14,*) "0 3" !3d2 4s2
		if (a(fragatm(i))%name=="V ") write(14,*) "0 4" !3d3 4s2
		if (a(fragatm(i))%name=="Cr") write(14,*) "0 7" !3d5 4s1, needn't sphericalization
		if (a(fragatm(i))%name=="Mn") write(14,*) "0 6" !3d5 4s2, needn't sphericalization
		if (a(fragatm(i))%name=="Fe") write(14,*) "0 5" !3d6 4s2
		if (a(fragatm(i))%name=="Co") write(14,*) "0 4" !3d7 4s2
		if (a(fragatm(i))%name=="Ni") write(14,*) "0 3" !3d8 4s2
		if (a(fragatm(i))%name=="Cu") write(14,*) "0 2" !3d10 4s1, needn't sphericalization
		if (a(fragatm(i))%name=="Zn") write(14,*) "0 1" !3d10 4s2, needn't sphericalization
		write(14,*) a(fragatm(i))%name,0.0,0.0,0.0
		write(14,*)
		write(14,*) trim(tmpdir)//a(fragatm(i))%name//".wfn" !The output path of wfn file
		write(14,*)
		write(14,*)
		close(14)
	end if
end do

if (noatmwfn==0) then
	call delfile(trim(tmpdir)//"*.gjf") !The .gjf generated have valueless now, delete them for avoiding user's misunderstanding
else if (noatmwfn==1) then !Some wfn needs to be genereated by Gaussian and sphericalized here
	do i=1,nfragatm
		nametmp=a_org(fragatm_org(i))%name
		inquire(file=trim(tmpdir)//nametmp//".wfn",exist=alive)
		if (alive) cycle !If the .wfn file had copied from atomwfn folder, needn't recalculate

		call runcommand('"'//trim(gaupath)//'" "'//trim(tmpdir)//nametmp//'.gjf" "'//trim(tmpdir)//nametmp//'"')
		!Check if Gaussian task was successfully finished
		if (isys==1) inquire(file=trim(tmpdir)//trim(nametmp)//".out",exist=alivegauout)
		if (isys==2) inquire(file=trim(tmpdir)//trim(nametmp)//".log",exist=alivegauout)
		if (alivegauout) then
			if (isys==1) open(10,file=trim(tmpdir)//trim(nametmp)//".out",status="old")
			if (isys==2) open(10,file=trim(tmpdir)//trim(nametmp)//".log",status="old")
			call loclabel(10,"Normal termination",igaunormal)
			close(10)
			if (igaunormal==0) then
				write(*,"(a)") " Gaussian running may be failed! Please manually check Gaussian input and output files in wfntmp folder. Press ENTER button to continue"
				read(*,*)
			else if (igaunormal==1) then
				write(*,*) "Finished successfully!"
			end if
		else
			write(*,"(a)") " Gaussian running may be failed! Please manually check Gaussian input and output files in wfntmp folder"
			read(*,*)
		end if
	
		!Load and sphericalize electron density for the just generated wfn, and then save
		if (ispheratm==1.and.igaunormal==1) then
			call dealloall(1)
			call readwfn(trim(tmpdir)//nametmp//".wfn",1)
			!Main group, restrict open-shell
			if (nametmp=="H ".or.nametmp=="Li".or.nametmp=="Na".or.nametmp=="K") MOocc(nmo)=1D0
			if (nametmp=="B ".or.nametmp=="Al".or.nametmp=="Ga") then
				nmo=nmo+2
				call resizebynmo(nmo,nprims) !Enlarge nmo by 2, but don't interfere nprims
				MOene(nmo-1:nmo)=MOene(nmo-2)
				MOtype(nmo-2:nmo)=1 !actually no use, because we only use atomic wfn. files to get total density
				MOocc(nmo-2:nmo)=1D0/3D0
				call orbcoeffrotate(nmo-2) !XYZ->ZXY, note: nmo-2 is original single occupied orbital
				CO(nmo-1,:)=CO(nmo-2,:)
				call orbcoeffrotate(nmo-2) !ZXY->YZX
				CO(nmo,:)=CO(nmo-2,:)
				call orbcoeffrotate(nmo-2) !YZX->XYZ, namely recovered
				!Now nmo-2,nmo-1,nmo correspond XYZ,ZXY,YZX
			end if
			if (nametmp=="C ".or.nametmp=="Si".or.nametmp=="Ge") then
				if (SpherIVgroup==0) then
					MOocc(nmo-3:nmo)=1D0
				else if (SpherIVgroup==1) then
					nmo=nmo+1
					call resizebynmo(nmo,nprims)
					MOene(nmo)=MOene(nmo-2) !MOene(nmo-1) is degenerate to MOene(nmo-2)
					MOtype(nmo-2:nmo)=1
					MOocc(nmo-2:nmo)=2D0/3D0
					!Rotate and copy the first occupied p orbital (nmo-5)
					call orbcoeffrotate(nmo-2) !XYZ->ZXY
					CO(nmo-1,:)=CO(nmo-2,:) !Overlap the already occupied orbital
					call orbcoeffrotate(nmo-2) !ZXY->YZX
					CO(nmo,:)=CO(nmo-2,:)
					call orbcoeffrotate(nmo-2) !YZX->XYZ, namely recovered
				end if
			end if
			if (nametmp=="N ".or.nametmp=="P ".or.nametmp=="As") MOocc(nmo-2:nmo)=1D0
			if (nametmp=="O ".or.nametmp=="S ".or.nametmp=="Se") MOocc(nmo-2:nmo)=4D0/3D0
			if (nametmp=="F ".or.nametmp=="Cl".or.nametmp=="Br") MOocc(nmo-2:nmo)=5D0/3D0
			!Transition metals, unrestrict open-shell, find boundary of alpha and beta first
			do ibound=2,nmo
				if (MOene(ibound)<MOene(ibound-1)) exit !from ii is beta orbitals
			end do
			!For Sc, Ti and V, rotate and duplicate d orbitals in each diection to get *near* spherical density, as for III main group
			!Note: Don't use Hartree-Fock, because correct energy sequence couldn't be reproduced, so can't be sphericalized correctly!
			if (nametmp=="Sc".or.nametmp=="Ti".or.nametmp=="V ") then !3d1 4s2, 3d2 4s2, 3d3 4s2
				if (nametmp=="Sc") then
					ibeg=1 !alpha 4s orbital, because this s orbital shows very strong unequlitity
					iend=2
				else if (nametmp=="Ti") then
					ibeg=2
					iend=3
				else if (nametmp=="V") then
					ibeg=2
					iend=4
				end if
				ienlarge=(iend-ibeg+1)*2
				call resizebynmo(nmo+ienlarge,nprims)
				ipass=0
				do iavgorb=ibeg,iend
					call orbcoeffrotate(ibound-iavgorb) !rotate this orbital
					CO(nmo+1+ipass,:)=CO(ibound-iavgorb,:) !Duplicate this orbital
					call orbcoeffrotate(ibound-iavgorb)
					CO(nmo+2+ipass,:)=CO(ibound-iavgorb,:)
					call orbcoeffrotate(ibound-iavgorb) !recover
					MOocc(ibound-iavgorb)=1D0/3D0
					MOene(nmo+1+ipass:nmo+2+ipass)=MOene(ibound-iavgorb)
					ipass=ipass+2 !next time skip nmo+1 and nmo+2
				end do
				MOocc(nmo+1:nmo+ienlarge)=1D0/3D0
				MOtype(nmo+1:nmo+ienlarge)=1
				nmo=nmo+ienlarge
			else if (nametmp=="Fe") then !3d6 4s2
				MOocc(nmo-1)=0D0 !delete the only d-beta orbital, the "nmo"th orbital is 4s-beta
				MOocc(ibound-6:ibound-2)=1.2D0 !Scatter one electron in beta-d orbital to alpha orbitals evenly. MOocc(ibound) is 4s orbital
			else if (nametmp=="Co") then !3d7 4s2
				MOocc(nmo-2:nmo-1)=0D0
				MOocc(ibound-6:ibound-2)=1.4D0
			else if (nametmp=="Ni") then !3d8 4s2
				MOocc(nmo-3:nmo-1)=0D0
				MOocc(ibound-6:ibound-2)=1.6D0
			end if
			call outwfn(trim(tmpdir)//nametmp//".wfn",0,0,10)
		end if
	end do
end if
write(*,*)

!Setup custom operation array with current size
ncustommap=nfragatm_org
if (allocated(custommapname)) deallocate(custommapname)
if (allocated(customop)) deallocate(customop)
allocate(custommapname(ncustommap))
allocate(customop(ncustommap))

!Generate atomic wfn file from element wfn file, meanwhile take them into custom operation list
do i=1,itype !Scan each atomtype in current system
	call dealloall(1)
	call readwfn(trim(tmpdir)//typename(i)//".wfn",1)
	do j=1,nfragatm_org
		if (a_org(fragatm_org(j))%name==typename(i)) then !Find atoms attributed to current element
			a(1)%x=a_org(fragatm_org(j))%x !Modify the atomic .wfn, then output to new .wfn
			a(1)%y=a_org(fragatm_org(j))%y
			a(1)%z=a_org(fragatm_org(j))%z
			write(outwfnname,"(a2,i4,a4)") typename(i),fragatm_org(j),".wfn"
			call outwfn(trim(tmpdir)//outwfnname,0,0,10)
			custommapname(j)=trim(tmpdir)//outwfnname !Sequence is identical to atom in fragment
		end if
	end do
end do

call dealloall(1)
call readinfile(firstfilename,1)
end subroutine




!!------- Generate a promolecular wavefunction by calculating and then combining atomic .wfn files, &
!store to global arrays with _pmol, namely MOocc_pmol, MOene_pmol, MOtype_pmol, CO_pmol, and then recover the original wavefunction
!itask=0: Only generate promolecular wavefunction information in memory (*_pmol); =1 Also export the promolecular wavefunction to a .wfn file
!  If wavefunction information is available, the basis set for evaluating atoms must be identical to present system, otherwise error will be reported!
!  Mainly used to calculate information gain at a batch of points, which needs evaluation of promolecular density frequently
!  The density corresponding to this promolecular wavefunction is exactly identical to superposition of densities of all atomic .wfn files
!  Since spin density is not interest in this context, spin flip is not taken into account
subroutine generate_promolwfn(itask)
use defvar
implicit real*8 (a-h,o-z)
integer itask
character selectyn
type(primtype),allocatable :: btmp(:)

!Generate atomic .wfn files
call setPromol

nmo_pmol=0
if (allocated(MOocc_pmol)) deallocate(MOocc_pmol,MOtype_pmol,MOene_pmol,CO_pmol)

!Allocate proper array sizes
if (allocated(b)) then
	allocate(MOocc_pmol(2*nmo),MOtype_pmol(2*nmo),MOene_pmol(2*nmo),CO_pmol(2*nmo,nprims)) !2*nmo is large enough for storing combined MOs
else !The original system does not have wavefunction information. Estimate sufficient size of nmo and prims
	nmotmp=0
	nprimstmp=0
	do iatm=1,ncustommap
		call dealloall(0)
		write(*,"(a)") " Dealing with "//trim(custommapname(iatm))
		call readinfile(custommapname(iatm),1)
		nmotmp=nmotmp+nmo
		nprimstmp=nprimstmp+nprims
	end do
    allocate(btmp(nprimstmp),MOocc_pmol(2*nmotmp),MOtype_pmol(2*nmotmp),MOene_pmol(2*nmotmp),CO_pmol(2*nmotmp,nprimstmp))
end if

MOocc_pmol=0
MOene_pmol=0
MOtype_pmol=0
CO_pmol=0

!Cycling each atomic .wfn file and merge into combined wavefunction
iGTF=1
imo_pmol=1
do iatm=1,ncustommap
    call dealloall(0)
    write(*,"(a)") " Dealing with "//trim(custommapname(iatm))
    call readinfile(custommapname(iatm),1)
    if (nprims_org/=0) then !The original system has wavefunction information
		if (iGTF+nprims-1>nprims_org) then
			write(*,"(a)") " Error: The basis set used for calculating atoms must be different to that originally used for calculating molecule!"
			write(*,*) "This situation is not supported. Press ENTER button to exit"
			read(*,*)
			stop
		end if
    end if
    MOocc_pmol(imo_pmol:imo_pmol+nmo-1)=MOocc(:)
    MOene_pmol(imo_pmol:imo_pmol+nmo-1)=MOene(:)
    MOtype_pmol(imo_pmol:imo_pmol+nmo-1)=MOtype(:)
    CO_pmol(imo_pmol:imo_pmol+nmo-1,iGTF:iGTF+nprims-1)=CO(:,:)
    if (allocated(btmp)) then !Original system does not have GTF information, now accumulate
		btmp(iGTF:iGTF+nprims-1)=b(:)
        btmp(iGTF:iGTF+nprims-1)%center=iatm
    end if
    iGTF=iGTF+nprims
    imo_pmol=imo_pmol+nmo
end do
nmo_pmol=imo_pmol-1
if (nprims_org/=0) then !The original system has wavefunction information
	if (iGTF-1/=nprims_org) then
		write(*,"(/,a)") " Error: The basis set used for calculating atoms must be different to that originally used for calculating molecule!"
		write(*,"(' Number of GTFs of promolecular wavefunction:',i10)") iGTF-1
		write(*,"(' Number of GTFs of original wavefunction:    ',i10)") nprims_org
		write(*,*) "This situation is not supported. Press ENTER button to exit"
		read(*,*)
		stop
	end if
end if

write(*,"(/,a)") " Done! Promolecular wavefunction information has been successfully generated in memory!"

write(*,"(/,a)") " Reloading "//trim(firstfilename)
call dealloall(1)
call readinfile(firstfilename,1)

if (itask==1) then
	write(*,"(/,a)") " Do you want to make wavefunction information in memory correspond to the just generated promolecular wavefunction? (y/n)"
    read(*,*) selectyn
    if (selectyn=='y'.or.selectyn=='Y') then
		if (allocated(b)) then
			deallocate(MOocc,MOtype,MOene,CO)
        else
			nprims=nprimstmp
            allocate(b(nprims))
            b(:)=btmp(:)
        end if
        allocate(MOocc(nmo_pmol),MOene(nmo_pmol),MOtype(nmo_pmol),CO(nmo_pmol,nprims))
		nmo=nmo_pmol
		MOocc=MOocc_pmol
		MOene=MOene_pmol
		MOtype=MOtype_pmol
		CO=CO_pmol
		!call outwfn("promol.wfn",1,1,10) !Note that the unoccupied MOs are automatically skipped
		!write(*,*) "Promolecular wavefunction has been outputted to promol.wfn in current folder"
        write(*,*) "Done! The analysis you performed later will correspond to promolecular case. You can also use option 0 in main function 6 to export the promolecular wavefunction as a .wfn file"
    end if
end if
end subroutine




!!!------------------ Generate density matrix, can be used when basis function information is available
!PS: density matrics can also be directly loaded from .fch via subroutine readfchdensmat
subroutine genP
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 tmpmat(nbasis,nbasis)

if (allocated(Ptot)) deallocate(Ptot)
if (allocated(Palpha)) deallocate(Palpha)
if (allocated(Pbeta)) deallocate(Pbeta)
allocate(Ptot(nbasis,nbasis))
Ptot=0
if (wfntype==1.or.wfntype==2.or.wfntype==4) then !open-shell
	allocate(Palpha(nbasis,nbasis))
	allocate(Pbeta(nbasis,nbasis))
	Palpha=0D0
	Pbeta=0D0
end if

!For SCF wavefunction, if the wavefunction has not been modified (imodwfn==0), use fast way to construct it
!However, if the wavefunction has been modified, the case may be complicated, for example, there is a hole orbital. In these cases
!We use general way (as used for post-HF) to construct density matrix
if (wfntype==0.and.imodwfn==0) then !RHF
	!Ptot=2*matmul(CObasa(:,1:nint(naelec)),transpose(CObasa(:,1:nint(naelec))))
    Ptot=matmul_blas(CObasa(:,1:nint(naelec)),CObasa(:,1:nint(naelec)),nbasis,nbasis,0,1) !Parallel MKL, fastest
    !call matprod(2,Ptot,CObasa(:,1:nint(naelec)),CObasa(:,1:nint(naelec)))
    Ptot=Ptot*2
else if (wfntype==1.and.imodwfn==0) then !UHF
    Palpha=matmul_blas(CObasa(:,1:nint(naelec)),CObasa(:,1:nint(naelec)),nbasis,nbasis,0,1)
    Pbeta=matmul_blas(CObasb(:,1:nint(nbelec)),CObasb(:,1:nint(nbelec)),nbasis,nbasis,0,1)
	Ptot=Palpha+Pbeta
else if (wfntype==2.and.imodwfn==0) then !ROHF
	Palpha=matmul(CObasa(:,1:nint(naelec)),transpose(CObasa(:,1:nint(naelec))))
	Pbeta=matmul(CObasa(:,1:nint(nbelec)),transpose(CObasa(:,1:nint(nbelec))))
	Ptot=Palpha+Pbeta
else if (wfntype==3.or.((wfntype==0.or.wfntype==2).and.imodwfn==1)) then !Restricted post-HF
	!A temporary matrix (tmpmat) is utilized, which effectively includes orbital occupancies into the coefficient matrix, &
	!this way is significantly cheaper than the equivalent commented code
	!do imo=1,nmo
	!	if (MOocc(imo)==0D0) cycle
	!	Ptot=Ptot+MOocc(imo)*matmul(CObasa(:,imo:imo),transpose(CObasa(:,imo:imo)))
	!end do
    do imo=1,nbasis
		tmpmat(:,imo)=CObasa(:,imo)*MOocc(imo)
    end do
    Ptot=matmul_blas(tmpmat(:,:),CObasa(:,:),nbasis,nbasis,0,1)
else if (wfntype==4.or.(wfntype==1.and.imodwfn==1)) then !Unrestricted post-HF
    do imo=1,nbasis
		tmpmat(:,imo)=CObasa(:,imo)*MOocc(imo)
    end do
    Palpha=matmul_blas(tmpmat(:,:),CObasa(:,:),nbasis,nbasis,0,1)
    do imo=1,nbasis
		tmpmat(:,imo)=CObasb(:,imo)*MOocc(imo+nbasis)
    end do
    Pbeta=matmul_blas(tmpmat(:,:),CObasb(:,:),nbasis,nbasis,0,1)
	Ptot=Palpha+Pbeta
end if
end subroutine




!!!------ Generate density matrix based on GTF
subroutine genPprim
use defvar
use util
implicit real*8 (a-h,o-z)

if (allocated(Ptot_prim)) deallocate(Ptot_prim)
allocate(Ptot_prim(nprims,nprims))
Ptot_prim=0

if (wfntype==1.or.wfntype==2.or.wfntype==4) then
    if (allocated(Palpha_prim)) deallocate(Palpha_prim,Pbeta_prim)
    allocate(Palpha_prim(nprims,nprims),Pbeta_prim(nprims,nprims))
    Palpha_prim=0
    Pbeta_prim=0
end if

naint=nint(naelec)
nbint=nint(nbelec)
if (wfntype==0.and.imodwfn==0) then !R wavefunction
    Ptot_prim=2*matmul_blas(transpose(CO(1:naint,:)),CO(1:naint,:),nprims,nprims,0,0)
    
else if (wfntype==1.and.imodwfn==0) then !U wavefunction
    do ibbeg=1,nmo
	    if (MOtype(ibbeg)==2) exit
    end do
    ibend=ibbeg-1+nbint
    Palpha_prim=matmul_blas(transpose(CO(1:naint,:)),CO(1:naint,:),nprims,nprims,0,0)
    if (nbint>0) then
        Pbeta_prim=matmul_blas(transpose(CO(ibbeg:ibend,:)),CO(ibbeg:ibend,:),nprims,nprims,0,0)
    end if
    Ptot_prim=Palpha_prim+Pbeta_prim
    
else if (wfntype==2.and.imodwfn==0) then !RO wavefunction
    Palpha_prim=matmul_blas(transpose(CO(1:naint,:)),CO(1:naint,:),nprims,nprims,0,0)
    if (nbint>0) then
        Pbeta_prim=matmul_blas(transpose(CO(1:nbint,:)),CO(1:nbint,:),nprims,nprims,0,0)
    end if
    Ptot_prim=Palpha_prim+Pbeta_prim
    
else if (wfntype==3.or.((wfntype==0.or.wfntype==2).and.imodwfn==1)) then !Restricted post-HF
    !$OMP PARALLEL DO SHARED(Ptot_prim) PRIVATE(iGTF,jGTF,imo) schedule(auto) NUM_THREADS(nthreads)
    do iGTF=1,nprims
        do jGTF=1,nprims
            do imo=1,nmo
                Ptot_prim(iGTF,jGTF) = Ptot_prim(iGTF,jGTF)+MOocc(imo)*CO(imo,iGTF)*CO(imo,jGTF)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    
else if (wfntype==4.or.(wfntype==1.and.imodwfn==1)) then !Unrestricted post-HF
    do ibbeg=1,nmo
	    if (MOtype(ibbeg)==2) exit
    end do
    !$OMP PARALLEL DO SHARED(Palpha_prim) PRIVATE(iGTF,jGTF,imo) schedule(auto) NUM_THREADS(nthreads)
    do iGTF=1,nprims
        do jGTF=1,nprims
            do imo=1,ibbeg-1
                Palpha_prim(iGTF,jGTF) = Palpha_prim(iGTF,jGTF)+MOocc(imo)*CO(imo,iGTF)*CO(imo,jGTF)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SHARED(Pbeta_prim) PRIVATE(iGTF,jGTF,imo) schedule(auto) NUM_THREADS(nthreads)
    do iGTF=1,nprims
        do jGTF=1,nprims
            do imo=ibbeg,nmo
                Pbeta_prim(iGTF,jGTF) = Pbeta_prim(iGTF,jGTF)+MOocc(imo)*CO(imo,iGTF)*CO(imo,jGTF)
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    Ptot_prim=Palpha_prim+Pbeta_prim
end if
end subroutine




!!!------ Show system one-electron properties based on density matrix and integral matrix between basis functions
subroutine sys1eprop
use defvar
if (.not.allocated(Sbas)) allocate(Sbas(nbasis,nbasis))
call genSbas_curr
write(*,"(' Total number of electrons:',f16.8)") sum(Ptot*Sbas)
if (.not.allocated(Tbas)) allocate(Tbas(nbasis,nbasis))
call genTbas_curr
write(*,"(' Kinetic energy:',f18.9,' a.u.')") sum(Ptot*Tbas)
if (.not.allocated(Dbas)) allocate(Dbas(3,nbasis,nbasis))
call genDbas_curr
write(*,"(' Electric dipole moment in X/Y/Z:',3f13.7,' a.u.')") sum(Ptot*Dbas(1,:,:)),sum(Ptot*Dbas(2,:,:)),sum(Ptot*Dbas(3,:,:))
if (.not.allocated(Magbas)) allocate(Magbas(3,nbasis,nbasis))
call genMagbas_curr
write(*,"(' Magnetic dipole moment in X/Y/Z:',3f13.7,' a.u.')") sum(Ptot*Magbas(1,:,:)),sum(Ptot*Magbas(2,:,:)),sum(Ptot*Magbas(3,:,:))
if (.not.allocated(Velbas)) allocate(Velbas(3,nbasis,nbasis))
call genVelbas_curr
write(*,"(' Linear momentum in X/Y/Z:       ',3f13.7,' a.u.')") sum(Ptot*Velbas(1,:,:)),sum(Ptot*Velbas(2,:,:)),sum(Ptot*Velbas(3,:,:))
end subroutine




!!!------------- Show all properties at a point
!ifuncsel: Controls the gradient and Hessian for which function
!ifileid: Output to which file destination, of course 6=screen
subroutine showptprop(inx,iny,inz,ifuncsel,ifileid)
use util
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 inx,iny,inz
real*8 eigvecmat(3,3),eigval(3),elegrad(3),elehess(3,3),funcgrad(3),funchess(3,3),stresstensor(3,3)
real*8 tmparr(3,1),tmpmat(3,3),tmpgrad1(3),tmpgrad2(3)
integer ifuncsel,ifileid

if (allocated(b)) then !If loaded file contains wavefuntion information
	call gencalchessmat(2,1,inx,iny,inz,elerho,elegrad,elehess) !Generate electron density, gradient and hessian
	write(ifileid,"(' Density of all electrons:',E18.10)") elerho
	if (ipolarpara==0) then
		tmpval=fspindens(inx,iny,inz,'s')
		write(ifileid,"(' Density of Alpha electrons:',E18.10)") (elerho+tmpval)/2D0
		write(ifileid,"(' Density of Beta electrons:',E18.10)") (elerho-tmpval)/2D0
		write(ifileid,"(' Spin density of electrons:',E18.10)") tmpval
	else if (ipolarpara==1) then
		write(ifileid,"(' Spin polarization parameter function:',E18.10)") fspindens(inx,iny,inz,'s')
	end if
	valG=lagkin(inx,iny,inz,0)
	valGx=lagkin(inx,iny,inz,1)
	valGy=lagkin(inx,iny,inz,2)
	valGz=lagkin(inx,iny,inz,3)
	write(ifileid,"(' Lagrangian kinetic energy G(r):',E18.10)") valG
	write(ifileid,"(' G(r) in X,Y,Z:',3E18.10)") valGx,valGy,valGz
	valK=Hamkin(inx,iny,inz,0)
	write(ifileid,"(' Hamiltonian kinetic energy K(r):',E18.10)") valK
! 	valKx=Hamkin(inx,iny,inz,1);valKy=Hamkin(inx,iny,inz,2);valKz=Hamkin(inx,iny,inz,3)
! 	write(ifileid,"(' K(r) in X,Y,Z:',3E18.10)") valKx,valKy,valKz
	write(ifileid,"(' Potential energy density V(r):',E18.10)") -valK-valG !When without EDF, also equals to flapl(inx,iny,inz,'t')/4D0-2*valG
	write(ifileid,"(' Energy density E(r) or H(r):',E18.10)") -valK
	write(ifileid,"(' Laplacian of electron density:',E18.10)") laplfac*(elehess(1,1)+elehess(2,2)+elehess(3,3))
	write(ifileid,"(' Electron localization function (ELF):',E18.10)") ELF_LOL(inx,iny,inz,"ELF")
	write(ifileid,"(' Localized orbital locator (LOL):',E18.10)") ELF_LOL(inx,iny,inz,"LOL")
	write(ifileid,"(' Local information entropy:',E18.10)") infoentro(1,inx,iny,inz)
	write(ifileid,"(' Interaction region indicator (IRI):',E18.10)") IRIfunc(inx,iny,inz)
	write(ifileid,"(' Reduced density gradient (RDG):',E18.10)") fgrad(inx,iny,inz,'r')
	write(ifileid,"(' Reduced density gradient with promolecular approximation:',E18.10)") RDGprodens(inx,iny,inz)
	write(ifileid,"(' Sign(lambda2)*rho:',E18.10)") signlambda2rho(inx,iny,inz)
	write(ifileid,"(' Sign(lambda2)*rho with promolecular approximation:',E18.10)") signlambda2rho_prodens(inx,iny,inz)
	if (pairfunctype==1) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. hole for alpha, ref.:",refx,refy,refz,pairfunc(refx,refy,refz,inx,iny,inz)
	if (pairfunctype==2) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. hole for beta, ref.:",refx,refy,refz,pairfunc(refx,refy,refz,inx,iny,inz)
	if (pairfunctype==4) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. fac. for alpha, ref.:",refx,refy,refz,pairfunc(refx,refy,refz,inx,iny,inz)
	if (pairfunctype==5) write(ifileid,"(a,3f10.5,' :',E18.10)") " Corr. fac. for beta, ref.:",refx,refy,refz,pairfunc(refx,refy,refz,inx,iny,inz)
	if (pairfunctype==7) write(ifileid,"(a,3f10.5,' :',E18.10)") " Exc.-corr. dens. for alpha, ref:",refx,refy,refz,pairfunc(refx,refy,refz,inx,iny,inz)
	if (pairfunctype==8) write(ifileid,"(a,3f10.5,' :',E18.10)") " Exc.-corr. dens. for beta, ref:",refx,refy,refz,pairfunc(refx,refy,refz,inx,iny,inz)
	write(ifileid,"(' Source function, ref.:',3f10.5,' :',E18.10)") refx,refy,refz,srcfunc(inx,iny,inz,srcfuncmode)
	if (nmo/=0) write(ifileid,"(' Wavefunction value for orbital',i8,' :',E18.10)") iorbsel,fmo(inx,iny,inz,iorbsel)
	if (iALIEdecomp==0) then
		write(ifileid,"(' Average local ionization energy (ALIE):',E18.10)") avglocion(inx,iny,inz)
	else if (iALIEdecomp==1) then
		call avglociondecomp(ifileid,inx,iny,inz)
	end if
	write(ifileid,"(' van der Waals potential (probe atom: ',a,'):',E18.10,' kcal/mol')") ind2name(ivdwprobe),vdwpotfunc(inx,iny,inz,1)
	write(ifileid,"(' Delta-g (under promolecular approximation):',E18.10)") delta_g_promol(inx,iny,inz)
	write(ifileid,"(' Delta-g (under Hirshfeld partition):',E18.10)") delta_g_Hirsh(inx,iny,inz)
	write(ifileid,"(' User-defined real space function:',E18.10)") userfunc(inx,iny,inz)
    if (ifPBC==0) then
	    fesptmp=nucesp(inx,iny,inz)
	    if (ifiletype==4) then
		    write(ifileid,"(' ESP from atomic charges:',E18.10)") fesptmp
	    else
		    write(ifileid,"(' ESP from nuclear charges:',E18.10)") fesptmp
	    end if
	    if (ishowptESP==1) then
		    fesptmpelec=eleesp(inx,iny,inz)
		    write(ifileid,"(' ESP from electrons:',E18.10)") fesptmpelec
		    write(ifileid,"(' Total ESP:',E18.10,' a.u. (',E14.7,' eV,',E14.7,' kcal/mol)')") &
		    fesptmpelec+fesptmp,(fesptmpelec+fesptmp)*au2eV,(fesptmpelec+fesptmp)*au2kcal
            if (iskipnuc/=0) then
		        tmpval=totespskip(inx,iny,inz,iskipnuc)
                write(ifileid,"(' Total ESP without contribution from nuclear charge of &
		        &atom',i6,':',/,E18.10,' a.u. (',E15.7,' eV,',E15.7,' kcal/mol)')") iskipnuc,tmpval,tmpval*au2eV,tmpval*au2kcal
            end if
        end if
    else if (ifPBC>0) then
        write(ifileid,"(a)") " Electrostatic potential (ESP) information is not shown because it cannot be calculated for periodic systems currently"
	end if

else !Only loaded structure, use YWT promolecule density
	call calchessmat_prodens(inx,iny,inz,elerho,elegrad,elehess)
	write(ifileid,"(/,a,/)") " Note: The input file does not contain wavefunction information, so the following quantities &
    &that related to electron density are evaluated based on promolecular density, which is constructed by the built-in &
    &free-state atomic densities described in Appendix 3 of Multiwfn manual"
	write(ifileid,"(' Density of electrons:',E18.10)") elerho
	write(ifileid,"(' Reduced density gradient:',E18.10)") RDGprodens(inx,iny,inz)
	write(ifileid,"(' Sign(lambda2)*rho:',E18.10)") signlambda2rho_prodens(inx,iny,inz)
	if (ifiletype==4) then
		write(ifileid,"(' ESP from atomic charges:',E18.10)") nucesp(inx,iny,inz)
	else
		write(ifileid,"(' ESP from nuclear charges:',E18.10)") nucesp(inx,iny,inz)
	end if
	write(ifileid,"(' van der Waals potential (probe atom: ',a,'):',E18.10,' kcal/mol')") ind2name(ivdwprobe),vdwpotfunc(inx,iny,inz,1)
	write(ifileid,"(' User-defined real space function:',E18.10)") userfunc(inx,iny,inz)
end if

write(ifileid,*)
if (ifuncsel==1) then
	write(ifileid,*) "Note: Below information is for electron density"
	funchess=elehess
	funcgrad=elegrad
else
	if (ifuncsel==3) then
        write(ifileid,*) "Note: Below information is for Laplacian of electron density"
	else if (ifuncsel==4) then
        write(ifileid,*) "Note: Below information is for value of orbital wavefunction"
	else if (ifuncsel==9) then
        write(ifileid,*) "Note: Below information is for electron localization function"
	else if (ifuncsel==10) then
        write(ifileid,*) "Note: Below information is for localized orbital locator"
	else if (ifuncsel==12) then
        write(ifileid,*) "Note: Below information is for total ESP"
	else if (ifuncsel==100) then
        write(ifileid,*) "Note: Below information is for user-defined real space function"
	else
        write(ifileid,"(a,i4)") " Note: Below information is for real space function",ifuncsel
    end if
	call gencalchessmat(2,ifuncsel,inx,iny,inz,funcvalue,funcgrad,funchess)
end if
write(ifileid,*)
write(ifileid,*) "Components of gradient in x/y/z are:"
write(ifileid,"(3E18.10)") funcgrad(1),funcgrad(2),funcgrad(3)
write(ifileid,"(' Norm of gradient is:',E18.10)") dsqrt(sum(funcgrad**2))
write(ifileid,*)
write(ifileid,*) "Components of Laplacian in x/y/z are:"
write(ifileid,"(3E18.10)") funchess(1,1),funchess(2,2),funchess(3,3)
write(ifileid,"(' Total:',E18.10)") funchess(1,1)+funchess(2,2)+funchess(3,3)
write(ifileid,*)
write(ifileid,*) "Hessian matrix:"
write(ifileid,"(3E18.10)") funchess
call diagsymat(funchess,eigvecmat,eigval,idiagok)
if (idiagok/=0) write(*,*) "Note: Diagonization of Hessian matrix failed!"
write(ifileid,"(' Eigenvalues of Hessian:',3E18.10)") eigval(1:3)
write(ifileid,*) "Eigenvectors (columns) of Hessian:"
write(ifileid,"(3E18.10)") ((eigvecmat(i,j),j=1,3),i=1,3)
write(ifileid,"(' Determinant of Hessian:',E18.10)") detmat(funchess)

if (ifuncsel==1) then
	!Output ellipticity for rho
	call sort(eigval) !Sort eigenvalues from low to high
	write(ifileid,"(a,f12.6)") " Ellipticity of electron density:",eigval(1)/eigval(2)-1
	write(ifileid,"(a,f12.6)") " eta index:",abs(eigval(1))/eigval(3)
	write(ifileid,"(a,f12.6)") " Stiffness:",abs(eigval(2))/eigval(3)
    
	!Output stress tensor and related quantities
    call stress_tensor(inx,iny,inz,stresstensor)
    write(ifileid,*)
	write(ifileid,*) "Stress tensor:"
	write(ifileid,"(3E18.10)") stresstensor
	call diagsymat(stresstensor,eigvecmat,eigval,idiagok)
	if (idiagok/=0) write(*,*) "Note: Diagonization of Hessian matrix failed!"
	write(ifileid,"(' Eigenvalues of stress tensor:',3E18.10)") eigval(1:3)
	write(ifileid,*) "Eigenvectors (columns) of stress tensor:"
	write(ifileid,"(3E18.10)") ((eigvecmat(i,j),j=1,3),i=1,3)
	call sort(eigval) !Sort eigenvalues from low to high
	write(ifileid,"(a,f12.6)") " Stress tensor stiffness:",abs(eigval(1))/abs(eigval(3))
	write(ifileid,"(a,f12.6)") " Stress tensor polarizability:",abs(eigval(3))/abs(eigval(1))
end if
end subroutine




!!!------- Decompose property at a point as contribution from various orbitals
subroutine decompptprop(x,y,z)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
character c2000tmp*2000
real*8 MOocctmp(nmo)
real*8,allocatable :: valarr(:)
integer,allocatable :: orbidx(:)
write(*,*) "Select the function to be studied"
call funclist
read(*,*) ifunc
write(*,*) "Input range of orbitals, e.g. 3,6-8,12-15"
write(*,"(a)") " Note: If press ENTER button directly, all occupied orbitals will be taken into account, &
&and 10 orbitals having largest contributions will be shown"
read(*,"(a)") c2000tmp
if (c2000tmp==" ") then
	norb=count(MOocc(1:nmo)/=0D0)
	allocate(orbidx(norb),valarr(norb))
	idx=0
	do imo=1,nmo
		if (MOocc(imo)/=0) then
			idx=idx+1
			orbidx(idx)=imo
		end if
	end do
else
	call str2arr(c2000tmp,norb)
	allocate(orbidx(norb),valarr(norb))
	call str2arr(c2000tmp,norb,orbidx)
end if

totval=calcfuncall(ifunc,x,y,z)
MOocctmp=MOocc
sumval=0
do itmp=1,norb
	iorb=orbidx(itmp)
	MOocc=0
	MOocc(iorb)=MOocctmp(iorb)
	valarr(itmp)=calcfuncall(ifunc,x,y,z)
	sumval=sumval+valarr(itmp)
end do
MOocc=MOocctmp

call sortr8(valarr,"abs",orbidx)
call invarrr8(valarr,orbidx)

nout=norb
if (c2000tmp==" ".and.norb>10) nout=10
do itmp=1,nout
	iorb=orbidx(itmp)
	if (ifunc==1) then !For electron density, also print percentage contribution
		write(*,"(' Contribution from orbital',i6,' (occ=',f9.6,'):',f14.6,' a.u. (',f6.2,'% )')") iorb,MOocc(iorb),valarr(itmp),valarr(itmp)/totval*100
	else
		write(*,"(' Contribution from orbital',i6,' (occ=',f9.6,'):',1E16.8,' a.u.')") iorb,MOocc(iorb),valarr(itmp)
	end if
end do

if (ifunc==1) then
	write(*,"(' Sum of above values:',f16.8,' a.u. ( ',f6.2,'% )')") sumval,sumval/totval*100
	write(*,"(' Exact value:',f16.8,' a.u.')") totval
else
	write(*,"(' Sum of above values:',1E16.8,' a.u.')") sumval
	write(*,"(' Exact value:',1E16.8,' a.u.')") totval
end if
end subroutine




!!!------------------- Delete virtual orbitals higher than LUMO+10 for HF/DFT wavefunctions
!Each time "delvirob" has been called, "delvirorb_back" should be then called to recover previous wavefunction status
!infomode=1 means show prompt
subroutine delvirorb(infomode)
use defvar
implicit real*8 (a-h,o-z)
integer :: infomode,nvirsave=10 !Lowest "nvirsave" virtual orbitals will be reserved
if (.not.allocated(CObasa)) return !Only works when all orbitals are available, which implies CObasa is allocated
if (idelvirorb==0) return !Do not delete orbitals as explicitly requested in settings.ini
!Linear response kernel, local electron affinity, orbital-weighted Fukui/dual descriptor require all orbital information
if (iuserfunc==24.or.iuserfunc==27.or.iuserfunc==28.or.iuserfunc==29.or.iuserfunc==95.or.iuserfunc==96.or.iuserfunc==97.or.iuserfunc==98) return
if (imodwfn==1) return !Wavefunction has been modified, so don't let this subroutine make thing more complicated!
if (wfntype==3.or.wfntype==4) return !This routine doesn't work for post-HF cases
if (ifdelvirorb==1) return !This routine has already been called while delvirorb_back was not used after that

if (allocated(CO_back)) deallocate(CO_back,MOene_back,MOocc_back,MOtype_back)
allocate(CO_back(nmo,nprims),MOene_back(nmo),MOocc_back(nmo),MOtype_back(nmo))
nmo_back=nmo
CO_back=CO
MOene_back=MOene
MOocc_back=MOocc
MOtype_back=MOtype

call getHOMOidx !Return idxHOMO (R, RO), or idxHOMO and idxHOMOb (U)
if (wfntype==0.or.wfntype==2) then !RHF, ROHF
	if (nmo<=idxHOMO+nvirsave) return
	nmo=idxHOMO+nvirsave !Simply shield those virtual orbitals
else if (wfntype==1) then !Reserve up to LUMO+10 for alpha, and same number of orbitals for beta
	if (nmo/2<=idxHOMO+nvirsave) return !naelec is always >= nbelec
	nperserve=idxHOMO+nvirsave !Number of lowest preserved alpha orbitals
	!CObasa and CObasb needn't to be modified, because they are not directly involved in real space function calculation
	CO(nperserve+1:2*nperserve,:)=CO(nmo/2+1:nmo/2+nperserve,:)
	MOene(nperserve+1:2*nperserve)=MOene(nmo/2+1:nmo/2+nperserve)
	MOocc(nperserve+1:2*nperserve)=MOocc(nmo/2+1:nmo/2+nperserve)
	MOtype(nperserve+1:2*nperserve)=MOtype(nmo/2+1:nmo/2+nperserve)
	nmo=2*nperserve
end if

ifdelvirorb=1 !delvirorb has taken effect
if (infomode==1) write(*,"(a,/)") " Note: Virtual orbitals higher than LUMO+9 have been temporarily discarded for saving computational time"
end subroutine



!!!------------------- Recover the status of wavefunction before calling delvirorb
!infomode=1 means show prompt, =0 means do not show
subroutine delvirorb_back(infomode)
use defvar
implicit real*8 (a-h,o-z)
integer infomode
if (ifdelvirorb==1) then
    nmo=nmo_back
    CO=CO_back
    MOene=MOene_back
    MOocc=MOocc_back
    MOtype=MOtype_back
    deallocate(CO_back,MOene_back,MOocc_back,MOtype_back)
    ifdelvirorb=0 !Has been restored
    if (infomode==1) write(*,"(a)") " Note: Previous orbital information has been restored"
end if
end subroutine




!!------- Generate neighbouring list of GTFs at reduced grids
!The grid positions are determined based on present grid (orgx/y/z to endx/y/z)
subroutine gen_neigh_GTF
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 tvec(3)

if (.not.allocated(b)) return

spcred=1.5D0/b2a !Spacing of reduced grid, this is found to be optimal value
!Define a orthogonal box for reduced grid
call cellmaxxyz(xmax,ymax,zmax)
call cellminxyz(orgx_neigh,orgy_neigh,orgz_neigh)
nx_red=floor((xmax-orgx_neigh)/spcred) !after +1 is number of reduced grids
ny_red=floor((ymax-orgy_neigh)/spcred)
nz_red=floor((zmax-orgz_neigh)/spcred)
if (allocated(neighGTF)) then
    if (size(neighnGTF,1)==nx_red+1 .and. size(neighnGTF,2)==ny_red+1 .and. size(neighnGTF,3)==nz_red+1) return !Proper arrays have already been allocated
    deallocate(neighnGTF,neighGTF,neighGTFcell)
end if

!write(*,"(' X/Y/Z max:',3f12.6)") xmax,ymax,zmax
!write(*,"(' X/Y/Z min:',3f12.6)") orgx_neigh,orgy_neigh,orgz_neigh
!write(*,"(' Number of reduced grids in X, Y, Z:',3i5)") nx_red,ny_red,nz_red
allocate(neighnGTF(0:nx_red,0:ny_red,0:nz_red))
cencordist=spcred/2*dsqrt(3D0) !Distance between center position and corner of (cubic) reduced grid
cencordist2=cencordist**2

!itime=1: Count maximum number of neighbouring GTFs at all grids so that arrays can be located
!itime=2: Filling neighGTF,neighnGTF,neighGTFcell
write(*,*) "Constructing neighbouring list of GTFs at reduced grids..."
call walltime(iwalltime1)
maxneighGTF=0
do itime=1,2
    !Loop reduced grids
    !$OMP PARALLEL DO SHARED(neighGTF,neighGTFcell) PRIVATE(x,y,z,ix,iy,iz,neighnGTF_tmp,icell,jcell,kcell,ic,jc,kc,lastcen,tvec,iGTF,icen,sftx,sfty,sftz,rr,iadd,rrtmp,tmpval) &
    !$OMP schedule(dynamic) NUM_THREADS(nthreads)
    do iz=0,nz_red
        z=orgz_neigh+iz*spcred+spcred/2
        do iy=0,ny_red
            y=orgy_neigh+iy*spcred+spcred/2
            do ix=0,nx_red
                x=orgx_neigh+ix*spcred+spcred/2
                call getpointcell(x,y,z,ic,jc,kc)
                neighnGTF_tmp=0
                !Loop cells
                do icell=ic-PBCnx,ic+PBCnx
                    do jcell=jc-PBCny,jc+PBCny
                        do kcell=kc-PBCnz,kc+PBCnz
			                lastcen=-1 !Arbitrary value
                            call tvec_PBC(icell,jcell,kcell,tvec)
                            !Loop GTFs
                            do iGTF=1,nprims
					            icen=b(iGTF)%center
					            if (icen/=lastcen) then
						            sftx=x-(a(icen)%x+tvec(1))
						            sfty=y-(a(icen)%y+tvec(2))
						            sftz=z-(a(icen)%z+tvec(3))
                                    rr=sftx*sftx+sfty*sfty+sftz*sftz
					            end if
					            lastcen=icen
                                iadd=0
                                if (rr<cencordist2) then
                                    iadd=1
					            else
                                    rrtmp=(dsqrt(rr)-cencordist)**2
					                tmpval=-b(iGTF)%exp*rrtmp
                                    if (tmpval>expcutoff_PBC.or.expcutoff_PBC>0) iadd=1
					            end if
                                if (iadd==1) then
                                    neighnGTF_tmp=neighnGTF_tmp+1
                                    if (itime==2) then
                                        neighGTF(neighnGTF_tmp,ix,iy,iz)=iGTF
                                        neighGTFcell(1,neighnGTF_tmp,ix,iy,iz)=icell
                                        neighGTFcell(2,neighnGTF_tmp,ix,iy,iz)=jcell
                                        neighGTFcell(3,neighnGTF_tmp,ix,iy,iz)=kcell
                                    end if
                                end if
                            end do
                            !End loop GTFs
                        end do
                    end do
                end do
                !End loop cells
                if (itime==1) neighnGTF(ix,iy,iz)=neighnGTF_tmp
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    !End loop reduced grids
    if (itime==1) then
        maxneighGTF=maxval(neighnGTF(:,:,:))
        !write(*,"(' Maximum number of neighbouring GTFs at a reduced grid:',i10)") maxneighGTF
        allocate(neighGTF(maxneighGTF,0:nx_red,0:ny_red,0:nz_red))
        allocate(neighGTFcell(3,maxneighGTF,0:nx_red,0:ny_red,0:nz_red))
    end if
end do
call walltime(iwalltime2)
!write(*,"(' Constructing neighbouring list of GTFs took up',i10,' s')") iwalltime2-iwalltime1
end subroutine




!!----------- Generate unique GTF information (b_uniq, CO_uniq)
!This reduces number of GTFs if generally contracted basis set is used, then cost in subroutine "orbserv" can be lowered
!infomode=0: Print number of unique GTFs. =1: Do not print
!Note that after using it, del_GTFuniq should be invoked as early as possible, otherwise when CO is changed due to some reasons, CO_uniq will be out-of-dated
subroutine gen_GTFuniq(infomode)
use defvar
integer ifcombined(nprims)

if (.not.allocated(b)) return
if (allocated(b_uniq)) deallocate(b_uniq)
if (allocated(CO_uniq)) deallocate(CO_uniq)

do itime=1,2
	nprims_uniq=0
	ifcombined=0
	do iGTF=1,nprims
		if (ifcombined(iGTF)==0) then
			nprims_uniq=nprims_uniq+1
            if (itime==2) then
				b_uniq(nprims_uniq)=b(iGTF)
				CO_uniq(:,nprims_uniq)=CO(1:nmo,iGTF)
            end if
			do jGTF=iGTF+1,nprims
				if (b(jGTF)%center==b(iGTF)%center.and.b(jGTF)%type==b(iGTF)%type.and.b(jGTF)%exp==b(iGTF)%exp) then
					if (itime==2) CO_uniq(:,nprims_uniq)=CO_uniq(:,nprims_uniq)+CO(1:nmo,jGTF)
					ifcombined(jGTF)=1
				end if
			end do
        end if
	end do
    if (itime==1) allocate(b_uniq(nprims_uniq),CO_uniq(nmo,nprims_uniq))
end do

if (infomode==0) write(*,"(' Unique GTFs have been constructed. Number of unique GTFs:',i8)") nprims_uniq
end subroutine



!!----------- Destory unique GTF information (b_uniq, CO_uniq)
subroutine del_GTFuniq
use defvar
nprims_uniq=0
if (allocated(b_uniq)) deallocate(b_uniq)
if (allocated(CO_uniq)) deallocate(CO_uniq)
end subroutine



!!!-------- imode=1: Convert unit of grid/plane parameters from Bohr to Angstrom. =2: Convert them back
subroutine convgridlenunit(imode)
use defvar
implicit none
integer imode
real*8 scll
if (imode==1) scll=b2a
if (imode==2) scll=1/b2a
orgx=orgx*scll
orgy=orgy*scll
orgz=orgz*scll
orgx2D=orgx2D*scll
orgy2D=orgy2D*scll
orgz2D=orgz2D*scll
endx=endx*scll
endy=endy*scll
endz=endz*scll
dx=dx*scll;gridv1=gridv1*scll
dy=dy*scll;gridv2=gridv2*scll
dz=dz*scll;gridv3=gridv3*scll
v1x=v1x*scll
v1y=v1y*scll
v1z=v1z*scll
v2x=v2x*scll
v2y=v2y*scll
v2z=v2z*scll
a1x=a1x*scll
a1y=a1y*scll
a1z=a1z*scll
a2x=a2x*scll
a2y=a2y*scll
a2z=a2z*scll
a3x=a3x*scll
a3y=a3y*scll
a3z=a3z*scll
d1=d1*scll
d2=d2*scll
end subroutine



!!-------- Deallocate all arrays about wavefunction except for the _org ones, and re-initialize some variables
!imode=0: Deallocate all global arrays
!imode=1: Same as 1 but do not deallocate frag1, frag2 and fragatm, which may be used later
subroutine dealloall(imode)
use defvar
integer imode

call delvirorb_back(0) !If delvirorb has taken effect, use this routine to deallocate relevant arrays
call del_GTFuniq
call dealloEDF !Deallocate EDF information
if (allocated(a)) deallocate(a)
if (allocated(b)) deallocate(b)
if (allocated(CO)) deallocate(CO)
if (allocated(MOocc)) deallocate(MOocc)
if (allocated(MOsym)) deallocate(MOsym)
if (allocated(MOene)) deallocate(MOene)
if (allocated(MOtype)) deallocate(MOtype)
if (allocated(connmat)) deallocate(connmat)
!Related to basis functions
if (allocated(shtype)) deallocate(shtype,shcen,shcon,primshexp,primshcoeff)
if (allocated(basshell)) deallocate(basshell,bascen,bastype,basstart,basend,primstart,primend,primconnorm)
if (allocated(shtypeCar)) deallocate(shtypeCar)
if (allocated(CObasa)) deallocate(CObasa)
if (allocated(CObasb)) deallocate(CObasb)
if (allocated(Ptot)) deallocate(Ptot)
if (allocated(Palpha)) deallocate(Palpha)
if (allocated(Pbeta)) deallocate(Pbeta)
if (allocated(Sbas)) deallocate(Sbas)
if (allocated(Dbas)) deallocate(Dbas)
if (allocated(DorbA)) deallocate(DorbA)
if (allocated(DorbB)) deallocate(DorbB)
if (allocated(Magbas)) deallocate(Magbas)
if (allocated(MagorbA)) deallocate(MagorbA)
if (allocated(MagorbB)) deallocate(MagorbB)
if (allocated(Quadbas)) deallocate(Quadbas)
if (allocated(Octobas)) deallocate(Octobas)

if (allocated(Ptot_prim)) deallocate(Ptot_prim)
if (allocated(Palpha_prim)) deallocate(Palpha_prim)
if (allocated(Pbeta_prim)) deallocate(Pbeta_prim)
if (allocated(Dprim)) deallocate(Dprim)
if (allocated(Quadprim)) deallocate(Quadprim)
if (allocated(Octoprim)) deallocate(Octoprim)

if (allocated(neighGTF)) deallocate(neighGTF,neighnGTF,neighGTFcell)

if (imode==0) then
	if (allocated(frag1)) deallocate(frag1)
	if (allocated(frag2)) deallocate(frag2)
	if (allocated(fragatm)) deallocate(fragatm)
end if

!Arrays used by delvirorb for backing up original MO information
if (allocated(CO_back)) deallocate(CO_back,MOene_back,MOocc_back,MOtype_back)

loadmulti=-99
loadcharge=-99
totenergy=0
virialratio=2
nelec=0
naelec=0
nbelec=0
nmo=0
nprims=0
ncenter=0
iresinfo=0
nbasis=0
nindbasis=0
kp1crd=0;kp2crd=0;kp3crd=0

ifPBC=0
cellv1=0;cellv2=0;cellv3=0
if (allocated(Sbas_PBC)) deallocate(Sbas_PBC)

call deallo_excitinfo !Deallocate excited state information
call deallo_basinana(0) !Deallocate basin analysis information
call deallo_topo !Clean topology analysis information

!Make LIBRETA in uninitialized status
if_initlibreta=0
nstates=0
end subroutine


!!-------- Deallocate all arrays about the system for the _org ones
subroutine dealloall_org
use defvar
firstfilename=" "
ncenter_org=0
nmo_org=0
nprims_org=0
if (allocated(a_org)) deallocate(a_org)
if (allocated(b_org)) deallocate(b_org,CO_org,MOocc_org,MOene_org)
if (allocated(Sbas_org)) deallocate(Sbas_org)
if (allocated(CObasa_org)) deallocate(CObasa_org)
if (allocated(CObasb_org)) deallocate(CObasb_org)
if (allocated(fragatm_org)) deallocate(fragatm_org)
ifPBC_org=0
cellv1_org=0
cellv2_org=0
cellv3_org=0
end subroutine


!!-------- Deallocate EDF related information
subroutine dealloEDF
use defvar
nEDFprims=0
nEDFelec=0
if (allocated(CO_EDF)) deallocate(CO_EDF)
if (allocated(b_EDF)) deallocate(b_EDF)
if (allocated(nEDFelecatm)) deallocate(nEDFelecatm)
end subroutine



!!------- Generate fragment Hirshfeld weighting function (based on atomic .wfn) and store it to planemat, calculate free-atom/fragmental density and store it to planemattmp
!The atoms in the fragment is inputted as "selatm" array, nselatm is the number of its elements
!if itype=1, use atomic wavefunction to calculate Hirshfeld weight, and setpromol must have been invoked; if =2, use built-in atomic density to generate it
subroutine genHirshplanewei(selatm,nselatm,itype)
use defvar
use functions
implicit real*8 (a-h,o-z)
integer selatm(nselatm),nselatm,itype

if (allocated(planemat)) deallocate(planemat)
if (allocated(planemattmp)) deallocate(planemattmp)
allocate(planemat(ngridnum1,ngridnum2),planemattmp(ngridnum1,ngridnum2))
planemat=0D0
planemattmp=0D0
do iatm=1,ncenter_org !Calc free atomic density of each atom, get promolecular density and Hirshfeld weight of present atom
	iyes=0
	if (any(selatm==iatm)) iyes=1
	if (itype==1) then
		call dealloall(0)
		call readwfn(custommapname(iatm),1)
	end if
	!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz,tmpval) shared(planemat) schedule(dynamic) NUM_THREADS(nthreads)
	do i=1,ngridnum1 !First calculate promolecular density and store it to planemat
		do j=1,ngridnum2
            call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
			if (itype==1) then
				tmpval=fdens(rnowx,rnowy,rnowz)
			else
				tmpval=calcatmdens(iatm,rnowx,rnowy,rnowz,0)
			end if
			planemat(i,j)=planemat(i,j)+tmpval
			if (iyes==1) planemattmp(i,j)=planemattmp(i,j)+tmpval
		end do
	end do
	!$OMP END PARALLEL DO
end do
if (itype==1) then
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule)
end if

do i=1,ngridnum1 !Calculate Hirshfeld weighting function
	do j=1,ngridnum2
		if (planemat(i,j)/=0D0) then
			planemat(i,j)=planemattmp(i,j)/planemat(i,j)
		else
			planemat(i,j)=0D0
		end if
	end do
end do
end subroutine



!!----- Generate fragment Hirshfeld weighting function (based on atomic .wfn) and store it to cubmat
!The atoms in the fragment is inputted as "selatm" array, nselatm is the number of its elements
!if itype=1, use atomic wavefunction to calculate Hirshfeld weight, and setpromol must have been invoked; if =2, use built-in atomic density to generate it
subroutine genHirshcubewei(selatm,nselatm,itype)
use defvar
use functions
implicit real*8 (a-h,o-z)
integer selatm(nselatm),nselatm,itype

if (allocated(cubmat)) deallocate(cubmat)
if (allocated(cubmattmp)) deallocate(cubmattmp)
allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
cubmat=0D0
cubmattmp=0D0
do iatm=1,ncenter_org
	write(*,"(' Finished',i6,'  /',i6)") iatm,ncenter_org
	if (itype==1) then
		call dealloall(0)
		call readwfn(custommapname(iatm),1)
	end if
	!$OMP PARALLEL DO SHARED(cubmat,cubmattmp,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz,tmpval) schedule(dynamic) NUM_THREADS(nthreads)
	do k=1,nz !First calculate promolecular density and store it to cubmat
		do j=1,ny
			do i=1,nx
                call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
				if (itype==1) then
					tmpval=fdens(tmpx,tmpy,tmpz)
				else
					tmpval=calcatmdens(iatm,tmpx,tmpy,tmpz,0)
				end if
				cubmat(i,j,k)=cubmat(i,j,k)+tmpval
				if (any(selatm==iatm)) cubmattmp(i,j,k)=cubmattmp(i,j,k)+tmpval
			end do
		end do
	end do
	!$OMP END PARALLEL DO
end do
if (itype==1) then
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule)
end if

do k=1,nz !Calculate Hirshfeld weighting function
	do j=1,ny
		do i=1,nx
			if (cubmat(i,j,k)/=0D0) then
				cubmat(i,j,k)=cubmattmp(i,j,k)/cubmat(i,j,k)
			else
				cubmat(i,j,k)=0D0
			end if
		end do
	end do
end do
end subroutine



!!--- Generate single-center integration grid. Return iradcut and gridatm
subroutine gen1cintgrid(gridatm,iradcut)
use defvar
implicit real*8 (a-h,o-z)
integer iradcut
real*8 potx(sphpot),poty(sphpot),potz(sphpot),potw(sphpot)
type(content) gridatm(radpot*sphpot)

call Lebedevgen(sphpot,potx,poty,potz,potw)
iradcut=0 !Before where the radial points will be cut
parm=1D0 !Not adapted according to element, for balance description for all cases
do i=1,radpot !Combine spherical point&weights with second kind Gauss-Chebyshev method for radial part
	radx=cos(i*pi/(radpot+1))
	radr=(1+radx)/(1-radx)*parm !Becke transform
	radw=2*pi/(radpot+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
	gridatm( (i-1)*sphpot+1:i*sphpot )%x=radr*potx
	gridatm( (i-1)*sphpot+1:i*sphpot )%y=radr*poty
	gridatm( (i-1)*sphpot+1:i*sphpot )%z=radr*potz
	gridatm( (i-1)*sphpot+1:i*sphpot )%value=radw*potw
	if (radcut/=0D0.and.iradcut==0.and.radr<radcut) iradcut=i-1
end do
end subroutine



!!----- A general interface for generating weighting function value of an atom at its single-center integration grids
!iweitype=1: Hirshfeld atom weighting function
!iweitype=2: Becke atom weighting function with covr_tianlu radii
!iweitype=3: Becke atom weighting function with covr radii
!iweitype=4: Tian lu error function type atom weighting function
!iweitype=5: Tian lu Gaussian function type atom weighting function
subroutine gen1catmwei(iatm,iradcut,gridatm,weigrid,iweitype)
use defvar
implicit real*8 (a-h,o-z)
integer iatm,iradcut
real*8 weigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot)
if (ifPBC>0.and.(iweitype==2.or.iweitype==3)) then
	write(*,*) "Error: Becke atom weighting function does not support periodic system!"
    write(*,*) "Press ENTER button to exit program"
    read(*,*)
    stop
end if
!$OMP parallel do shared(weigrid) private(i,rnowx,rnowy,rnowz) num_threads(nthreads) schedule(dynamic)
do i=1+iradcut*sphpot,radpot*sphpot
	rnowx=gridatm(i)%x
	rnowy=gridatm(i)%y
	rnowz=gridatm(i)%z
    if (iweitype==1) then
		call Hirshatmwei(iatm,rnowx,rnowy,rnowz,weigrid(i))
    else if (iweitype==2) then
		call Beckeatmwei(iatm,rnowx,rnowy,rnowz,weigrid(i),covr_tianlu,3)
    else if (iweitype==3) then
		call Beckeatmwei(iatm,rnowx,rnowy,rnowz,weigrid(i),covr,3)
    else if (iweitype==4) then
		call TLatmwei(iatm,rnowx,rnowy,rnowz,weigrid(i),1)
    else if (iweitype==5) then
		call TLatmwei(iatm,rnowx,rnowy,rnowz,weigrid(i),2)
    end if
end do
!$OMP end parallel do
end subroutine



!!------ Generate Becke weight for a batch of points around iatm
!Input: iatm, iradcut, gridatm, radinp, nbeckeiter   Return: beckeweigrid
!  radinp(1:nsuppele) are radii used for determining Becke atomic weights. In principle, using covr_tianlu seems reasonable, &
!compared to using covr, this allows larger space for more electronegative atoms and thus covering their rich electron distribution. &
!However using covr is frequently found to result in better accuracy.
!  nbeckeiter is sharpness parameter, without special reason is should be 3
!In fact, this routine is equivalent to "call gen1catmwei(iatm,iradcut,gridatm,weigrid,2)", but present function allows to directly set the radii used
subroutine gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,radinp,nbeckeiter)
use defvar
implicit real*8 (a-h,o-z)
integer iatm,iradcut,nbeckeiter
real*8 beckeweigrid(radpot*sphpot),radinp(0:nelesupp)
type(content) gridatm(radpot*sphpot)

!$OMP parallel do shared(beckeweigrid) private(i) num_threads(nthreads) schedule(dynamic)
do i=1+iradcut*sphpot,radpot*sphpot
	call Beckeatmwei(iatm,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,beckeweigrid(i),radinp,nbeckeiter)
end do
!$OMP end parallel do
end subroutine



!!------- Calculate Becke weighting function of iatm at (x,y,z). Does not support PBC
subroutine Beckeatmwei(iatm,x,y,z,weight,radinp,nbeckeiter)
use defvar
integer iatm,nbeckeiter
real*8 x,y,z,weight,radinp(0:nelesupp),Pvec(ncenter)
call BeckePvec(x,y,z,Pvec,radinp,nbeckeiter)
weight=Pvec(iatm)
end subroutine
!!------- Calculate Becke weighting function of all atoms at (x,y,z), storing to Pvec(:). Does not support PBC
subroutine BeckePvec(x,y,z,Pvec,radinp,nbeckeiter)
use defvar
use util
implicit real*8 (a-h,o-z)
integer nbeckeiter
real*8 x,y,z,Pvec(ncenter)
real*8 smat(ncenter,ncenter),rdist(ncenter),atmrad(ncenter),radinp(0:nelesupp)

do i=1,ncenter
	atmrad(i)=radinp(a(i)%index) !Set actually used radii
	rdist(i)=dsqrt( (x-a(i)%x)**2+(y-a(i)%y)**2+(z-a(i)%z)**2 ) !Distance between current point to every atom
end do

smat(:,:)=1D0
    
!Calculate weight for atom pair (ii,jj). smat(ii,jj) is weight of ii, smat(jj,ii) is weight of jj, they sum to 1
do ii=1,ncenter
	do jj=ii+1,ncenter
        dist=atomdist(jj,ii,1)
        diff=rdist(ii)-rdist(jj)
            
        !Quick determination. If difference of distance between current point to ii and to jj is longer than their distance, atom of distant side should have weight of zero
        !However, practical tested showed that this doesn't accelerate computation
        if (diff>dist) then
			smat(ii,jj)=0
            cycle
        else if (-diff>dist) then
			smat(jj,ii)=0
            cycle
        end if
            
		rmiu=diff/dist
            
 		!Adjust for heteronuclear, if their distance is close. In principle this adjustment should apply to all atomic pairs,
        !however this will cause computational burden if the system is large. In fact, adjusting this for two distant atoms is basically meaningless
        if (dist<8.and.a(ii)%index/=a(jj)%index) then
			chi=atmrad(ii)/atmrad(jj)
			uij=(chi-1)/(chi+1)
			aij=uij/(uij*uij-1)
			if (aij>0.5D0) then
				aij=0.5D0
			else if (aij<-0.5D0) then
				aij=-0.5D0
			end if
			rmiu=rmiu+aij*(1-rmiu**2)
        end if
            
        !Sharpness parameter. This is major overhead
        if (nbeckeiter==3) then !Common case, faster than using small loop
			tmps=1.5D0*rmiu-0.5D0*rmiu**3
			tmps=1.5D0*tmps-0.5D0*tmps**3
			tmps=1.5D0*tmps-0.5D0*tmps**3
        else
			tmps=rmiu
			do iter=1,nbeckeiter
				tmps=1.5D0*tmps-0.5D0*tmps**3
			end do
        end if
            
		smat(ii,jj)=0.5D0*(1-tmps)
        smat(jj,ii)=1-smat(ii,jj)
	end do
end do
    
!Pvec(i) is unnormalized weight of atom i at this point
Pvec=1D0
do i=1,ncenter
	Pvec(:)=Pvec(:)*smat(:,i)
end do
!Normalize weights
tmp=sum(Pvec)
Pvec(:)=Pvec(:)/tmp
end subroutine



!!------- Calculate Hirshfeld weighting function of iatm at (x,y,z)
!PBC is taken into account because function "calcatmdens" considers PBC
subroutine Hirshatmwei(iatm,x,y,z,weight)
use defvar
use functions
implicit real*8 (a-h,o-z)
integer iatm
real*8 x,y,z,weight

promol=0
do jatm=1,ncenter
	tmpdens=calcatmdens(jatm,x,y,z,0)
	promol=promol+tmpdens
	if (jatm==iatm) selfdens=tmpdens
end do

!ifPBC=0
!selfdens=calcatmdens(iatm,x,y,z,0)
!ifPBC=3

if (promol==0) then
	weight=0
	!Unable to determine atomic weights at this point because promol is zero, mostly becaues this point is too far from any atom
    !In this case, find the atom closest to this point, make its weight to 1. However this makes sharp change of weighting function at boundary, so not employed
	!call closest_atm_pt(x,y,z,iatmclose,r)
 !   if (iatmclose==iatm) weight=1
else
	weight=selfdens/promol
end if
end subroutine



!!------- Calculate Tian Lu weighting function of iatm at (x,y,z). PBC is taken into account
!Calculate value of simple atomic decaying function for all atoms, and finally calculate weight of iatm using Hirshfeld-like manner
!itype=1: Error function type. xscale = 0.85, leading to modest sharpness. Weight of 0.5 equals to CSD covalent radii
!itype=2: Gaussian function type. FWHM is CSD covalent radii
!itype=3: Becke function type. Weight of 0.5 equals to CSD covalent radii. This is poor, because this function decays quickly to zero, making weighting function at distant region cannot be calculated
!See http://sobereva.com/539 for illustration of different weighting functions
!
!All the functions have deficiency, namely they become exactly zero at distant region, &
!in this case weight cannot be determined and will be default to zero, if they used in multi-center integration purpose, very distant regions will be omitted
!I found none 
subroutine TLatmwei(iatm,x,y,z,weight,itype)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
integer iatm,itype
real*8 x,y,z,weight,tvec(3)

promol=0
selfdens=0
xscale=0.85D0

if (ifPBC==0) then
    do jatm=1,ncenter
		r=dsqrt( (a(jatm)%x-x)**2 + (a(jatm)%y-y)**2 + (a(jatm)%z-z)**2 )
		if (itype==1) then
			tmpval=switch_erf(r,covr(a(jatm)%index),xscale)
		else if (itype==2) then
			tmpval=switch_Gauss(r,2*covr(a(jatm)%index))
		else if (itype==3) then
			tmpval=switch_Becke(r,covr(a(jatm)%index),3)
		end if
		promol=promol+tmpval
		if (jatm==iatm) selfdens=tmpval
	end do
else !Periodic case
	call getpointcell(x,y,z,ic,jc,kc)
	do icell=ic-PBCnx,ic+PBCnx
		do jcell=jc-PBCny,jc+PBCny
			do kcell=kc-PBCnz,kc+PBCnz
				call tvec_PBC(icell,jcell,kcell,tvec)
				do jatm=1,ncenter
					atmx=a(jatm)%x+tvec(1)
					atmy=a(jatm)%y+tvec(2)
					atmz=a(jatm)%z+tvec(3)
					r=dsqrt( (atmx-x)**2 + (atmy-y)**2 + (atmz-z)**2 )
					if (itype==1) then
						tmpval=switch_erf(r,covr(a(jatm)%index),xscale)
					else if (itype==2) then
						tmpval=switch_Gauss(r,2*covr(a(jatm)%index))
					else if (itype==3) then
						tmpval=switch_Becke(r,covr(a(jatm)%index),3)
					end if
					promol=promol+tmpval
					if (jatm==iatm) selfdens=selfdens+tmpval
				end do
			end do
		end do
	end do
end if

if (promol==0) then
	weight=0
	!Unable to determine atomic weights at this point because promol is zero, mostly becaues this point is too far from any atom
    !In this case, find the atom closest to this point, make its weight to 1. However this makes sharp change of weighting function at boundary, so not employed
	call closest_atm_pt(x,y,z,iatmclose,r)
    if (iatmclose==iatm) weight=1
else
	weight=selfdens/promol
end if
end subroutine



!!--------- A standalone routine to calculate atomic contribution to specific real space function. Employed by IFCT analysis
!iparttype=1: Becke partition
!atmcontri: Returned array containing atomic contribution
!ifunc: The function to be evaluated
subroutine atmcontrifunc(iparttype,atmcontri,ifunc)
use defvar
use functions
implicit real*8 (a-h,o-z)
integer ifunc
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
real*8 atmcontri(ncenter),beckeweigrid(radpot*sphpot)

atmcontri=0
! write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
if (iparttype==1) then !Becke partition
	call gen1cintgrid(gridatmorg,iradcut)
	!$OMP PARALLEL DO SHARED(atmcontri) PRIVATE(iatm,ipt,gridatm,beckeweigrid,funcval) schedule(dynamic) NUM_THREADS(nthreads)
	do iatm=1,ncenter
! 		write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			funcval=calcfuncall(ifunc,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
			atmcontri(iatm)=atmcontri(iatm)+funcval*gridatmorg(ipt)%value*beckeweigrid(ipt)
		end do
	end do
	!$OMP END PARALLEL DO
end if
end subroutine



!!--------- Output all Becke's integration points to intpt.txt in current folder
subroutine outBeckeintpt
use defvar
implicit real*8 (a-h,o-z)
real*8 beckeweigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
character outwei

write(*,*) "See Appendix 6 of Multiwfn manual for detail of this function"
write(*,*)
write(*,*) "Also export Becke's integration weights? (y/n)"
read(*,*) outwei

open(10,file="intpt.txt",status="replace")
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
if (outwei=='y') then
    write(*,*) "Calculating and exporting, please wait..."
else
    write(*,*) "Exporting, please wait..."
end if
call gen1cintgrid(gridatmorg,iradcut)
do iatm=1,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
    if (outwei=='n') then
        do ipt=1+iradcut*sphpot,radpot*sphpot
            write(10,"(i6,4E16.8)") iatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,gridatmorg(ipt)%value
        end do
    else
	    call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
        do ipt=1+iradcut*sphpot,radpot*sphpot
            write(10,"(i6,5E16.8)") iatm,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,gridatmorg(ipt)%value,beckeweigrid(ipt)
        end do
    end if
end do  
close(10)
write(*,"(a)") " Done! All Becke's integration points have been exported to intpt.txt in current folder"
write(*,*) "Column 1: Index of the atom that the point belongs to"
write(*,*) "Column 2/3/4: X/Y/Z coordinate of the point in Bohr"
write(*,"(a)") " Column 5: Single center integration weight (second kind Gauss-Chebyshev for radial part and Lebedev for angular part"
if (outwei=='y') write(*,*) "Column 6: Becke's integration weight"
end subroutine



!!--------- Convert 1RDM in MO basis outputted by MRCC program to natural orbitals
!In CCDENSITIES, the density matrix is represented in MO basis
!When frozen core is enabled, the indices are counted from the first correlated orbital
subroutine MRCC_gennatorb
use defvar
use util
implicit real*8 (a-h,o-z)
character c200tmp*200,selectyn,outname*200
real*8,allocatable :: eigvecmat(:,:),eigvalarr(:),tmparr(:)
do while(.true.)
	write(*,*) "Input the path of CCDENSITIES file, e.g. C:\lovelive\CCDENSITIES"
    write(*,*) "If press ENTER button directly, CCDENSITIES in current folder will be loaded"
	read(*,"(a)") c200tmp
    if (c200tmp==" ") c200tmp="CCDENSITIES"
	inquire(file=c200tmp,exist=alive)
	if (alive) exit
	write(*,*) "Error: Cannot find the file, input again"
end do
write(*,*)
write(*,*) "Input the number of frozen orbitals, e.g. 3"
write(*,*) "If no orbitals are frozen, simply input 0"
write(*,"(a)") " PS: For unrestricted reference, if you input n, the n lowest alpha and n lowest beta MOs will be regarded as frozen"
read(*,*) nfrz
write(*,*) "Please wait..."
if (wfntype==0.or.wfntype==2) then !RHF and ROHF reference
	open(10,file=c200tmp,status="old")
	Ptot=0
	do while(.true.)
		read(10,*,iostat=ierror) tmp,i,j,k,l
		if (ierror/=0) exit
		if (k==0.and.l==0) then !Only load 1RDM
			Ptot(i+nfrz,j+nfrz)=tmp
			Ptot(j+nfrz,i+nfrz)=tmp
		end if
	end do
	close(10)
	do ifrz=1,nfrz
		Ptot(ifrz,ifrz)=2D0
	end do
	allocate(eigvecmat(nbasis,nbasis),eigvalarr(nbasis),tmparr(nbasis))
	call diagsymat(Ptot,eigvecmat,eigvalarr,istat)
	MOocc=eigvalarr
	!Currently the occupation is from low to high, now invert the sequence
	do i=1,int(nmo/2D0)
		idx=i
		jdx=nmo+1-i
		tmp=MOocc(idx)
		MOocc(idx)=MOocc(jdx)
		MOocc(jdx)=tmp
		tmparr=eigvecmat(:,idx)
		eigvecmat(:,idx)=eigvecmat(:,jdx)
		eigvecmat(:,jdx)=tmparr
	end do
	CObasa=matmul(CObasa,eigvecmat) !Transform to AO basis
	wfntype=3
	write(*,*) "Occupation numbers:"
	write(*,"(6f12.8)") MOocc
else if (wfntype==1) then !UHF reference
	!In CCDENSITIES, the sequence is:
	!2RDM-alpha
	!  0.00000000000000000000E+00   0   0   0   0
	!2RDM-beta
	!  0.00000000000000000000E+00   0   0   0   0
	!Unknown
	!  0.00000000000000000000E+00   0   0   0   0
	!1RDM-alpha
	!  0.00000000000000000000E+00   0   0   0   0
	!1RDM-beta
	!  0.00000000000000000000E+00   0   0   0   0
	open(10,file=c200tmp,status="old")
	Palpha=0
	Pbeta=0
	itime=0
	do while(.true.)
		read(10,*) tmp,i,j,k,l
		if (i==0.and.j==0.and.k==0.and.l==0) then
			itime=itime+1
			if (itime==5) exit
			cycle
		end if
		if (itime==3) then
			Palpha(i+nfrz,j+nfrz)=tmp
			Palpha(j+nfrz,i+nfrz)=tmp
		else if (itime==4) then
			Pbeta(i+nfrz,j+nfrz)=tmp
			Pbeta(j+nfrz,i+nfrz)=tmp
		end if
	end do
	close(10)
	do ifrz=1,nfrz
		Palpha(ifrz,ifrz)=1D0
		Pbeta(ifrz,ifrz)=1D0
	end do
	allocate(eigvecmat(nbasis,nbasis),eigvalarr(nbasis),tmparr(nbasis))
	!Alpha part
	call diagsymat(Palpha,eigvecmat,eigvalarr,istat)
	MOocc(1:nbasis)=eigvalarr
	do i=1,int(nbasis/2D0)
		idx=i
		jdx=nbasis+1-i
		tmp=MOocc(idx)
		MOocc(idx)=MOocc(jdx)
		MOocc(jdx)=tmp
		tmparr=eigvecmat(:,idx)
		eigvecmat(:,idx)=eigvecmat(:,jdx)
		eigvecmat(:,jdx)=tmparr
	end do
	CObasa=matmul(CObasa,eigvecmat)
	write(*,*) "Occupation number of Alpha part:"
	write(*,"(6f12.6)") MOocc(1:nbasis)
	!Beta part
	call diagsymat(Pbeta,eigvecmat,eigvalarr,istat)
	MOocc(nbasis+1:nmo)=eigvalarr
	do i=1,int(nbasis/2D0)
		idx=nbasis+i
		jdx=nmo+1-i
		tmp=MOocc(idx)
		MOocc(idx)=MOocc(jdx)
		MOocc(jdx)=tmp
		tmparr=eigvecmat(:,i)
		eigvecmat(:,i)=eigvecmat(:,nbasis+1-i)
		eigvecmat(:,nbasis+1-i)=tmparr
	end do
	CObasb=matmul(CObasb,eigvecmat)
	write(*,*) "Occupation number of Beta part:"
	write(*,"(6f12.6)") MOocc(nbasis+1:nmo)
	wfntype=4
end if

call genP
MOene=0

call path2filename(filename,c200tmp)
outname=trim(c200tmp)//".mwfn"
call outmwfn(outname,10,0)
write(*,"(/,a)") " All natural orbitals has been exported to "//trim(outname)//" in current folder"
write(*,"(a)") " Do you want to load it now so that you can perform wavefunction analysis based on the natural orbitals at the corresponding level? (y/n)"
read(*,*) selectyn
call dealloall(0)
if (selectyn=='y'.or.selectyn=='Y') then
    write(*,*) "Loading "//trim(outname)
    call readinfile(outname,1)
    write(*,"(a)") " Loading finished!"
else
    write(*,"(a)") " Reloading "//trim(filename)//" to recover initial status..."
    call readinfile(filename,1)
    write(*,*) "Loading finished!"
end if

end subroutine




!!--------- Convert 1RDM in json file outputted by orca_2json in ORCA program to natural orbitals
!gennatorb is invoked in this subroutine
!This subroutine is very similar with gennatorb
subroutine ORCAjson_gennatorb
use util
use defvar
implicit real*8 (a-h,o-z)
real*8,allocatable :: Pspin(:,:)
character c200tmp*200,selectyn,denstype*10,locstr*40

if (ifiletype/=9) then
	write(*,*) "Error: Molden file of ORCA must be used as input file for this function!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

do while(.true.)
	write(*,*) "Input the path of json file, e.g. C:\Palaio\Faliro.json"
	read(*,"(a)") c200tmp
	inquire(file=c200tmp,exist=alive)
	if (alive) exit
	write(*,*) "Error: Cannot find the file, input again"
end do
open(10,file=c200tmp,status="old")

write(*,*) "Input the type of density matrix, e.g. pmp2re, autocipre..."
write(*,"(a)") " Hint: If you do not know how to input, please run e.g. ""orca_plot test.gbw -i"" and then choose ""1 - Enter type of plot""&
& to check which densities are available. For example, pmp2re corresponds to MP2 relaxed density, autocipre corresponds to AutoCI relaxed density"
do while(.true.)
	read(*,"(a)") denstype
	call loclabel(10,trim(denstype),ifound)
	if (ifound==1) exit
    write(*,"(a)") " Error: Unable to find "//trim(denstype)//" from the json file! Input again"
end do

write(*,"(' Found ',a,', loading it...')") trim(denstype)
read(10,*)
do ibas=1,nbasis
	read(10,*)
	read(10,*) Ptot(:,ibas)
	read(10,*)
end do
call ORCA_mat_reorder(nbasis,Ptot)

iNOtype=1
if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	write(*,*) "Select the type of natural orbitals you want to obtain"
	write(*,*) "1 Spatial natural orbitals (diagonalizing total density matrix)"
	write(*,*) "2 Alpha and beta natural orbitals (diagonalizing respective density matrix)"
	write(*,*) "3 Spin natural orbitals (diagonalizing spin density matrix)"
	read(*,*) iNOtype
    if (iNOtype>1) then
		write(*,*) "Input the label corresponding to the spin density matrix, e.g.:"
        write(*,*) "rmp2re corresponds to MP2 relaxed spin density"
        write(*,*) "autocirre corresponds to AutoCI relaxed spin density"
        do while(.true.)
			read(*,"(a)") denstype
			call loclabel(10,trim(denstype),ifound)
			if (ifound==1) exit
			write(*,"(a)") " Error: Unable to find "//trim(denstype)//" from the json file! Input again"
		end do
		write(*,"(' Found ',a,', loading it...')") trim(denstype)
		!Load spin density matrix to construct alpha and beta DM
		allocate(Pspin(nbasis,nbasis))
		read(10,*)
		do ibas=1,nbasis
			read(10,*)
			read(10,*) Pspin(:,ibas)
			read(10,*)
		end do
		call ORCA_mat_reorder(nbasis,Pspin)
		Palpha=(Ptot+Pspin)/2D0
		Pbeta=(Ptot-Pspin)/2D0
    end if
end if

close(10)
write(*,*) "Density matrix was successfully loaded from the json file"

call gennatorb(iNOtype,1)
write(*,*) "Done! Basis function information now correspond to natural orbitals"

write(*,"(/,a)") " If next you intend to analyze real space functions based on the NOs, you should export new.mwfn &
&in current folder and then reload it, so that GTF information will also correspond to NOs"
write(*,*) "Would you like to do this immediately? (y/n)"
read(*,*) selectyn
if (selectyn=='y') then
    call outmwfn("new.mwfn",10,0)
	write(*,*) "The NOs have been exported to new.mwfn in current folder"
	call dealloall(0)
	write(*,*) "Loading new.mwfn..."
	call readinfile("new.mwfn",1)
	write(*,"(a)") " Loading finished, now you can use main function 0 to visualize NOs as isosurfaces, or perform various wavefunction analyses based on the NOs"
end if
end subroutine




!!----------- Generate spherical harmonic -> Cartesian basis function conversion table for d,f,g,h.
!iprog=1: for readfch;  iprog=2: for readmolden
!The table comes from IJQC,54,83, which is used by Gaussian
!The sequence of d and f shell is also identical to .molden convention, but for g, another conversion table is used, &
!since in Multiwfn g Cartesian shell starts from ZZZZ, but that of .molden starts from xxxx
subroutine gensphcartab(iprog,matd,matf,matg,math)
real*8 matd(6,5),matf(10,7),matg(15,9),math(21,11)
integer iprog

matd=0D0
matf=0D0
matg=0D0
math=0D0
! From 5D: D 0,D+1,D-1,D+2,D-2
! To 6D:  1  2  3  4  5  6
!        XX,YY,ZZ,XY,XZ,YZ
!
! D0=-0.5*XX-0.5*YY+ZZ
matd(1:3,1)=(/ -0.5D0,-0.5D0,1D0 /)
! D+1=XZ
matd(5,2)=1D0
! D-1=YZ
matd(6,3)=1D0
! D+2=SQRT(3)/2*(XX-YY)
matd(1:2,4)=(/ sqrt(3D0)/2D0,-sqrt(3D0)/2D0 /)
! D-2=XY
matd(4,5)=1D0

! From 7F: F 0,F+1,F-1,F+2,F-2,F+3,F-3
! To 10F:  1   2   3   4   5   6   7   8   9  10      
!         XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ (Gaussian sequence, not identical to Multiwfn)
!
! F 0=-3/(2*√5)*(XXZ+YYZ)+ZZZ
matf(3,1)=1D0
matf(6,1)=-1.5D0/sqrt(5D0)
matf(9,1)=-1.5D0/sqrt(5D0)
! F+1=-√(3/8)*XXX-√(3/40)*XYY+√(6/5)*XZZ
matf(1,2)=-sqrt(3D0/8D0)
matf(4,2)=-sqrt(3D0/40D0)
matf(7,2)=sqrt(6D0/5D0)
! F-1=-√(3/40)*XXY-√(3/8)*YYY+√(6/5)*YZZ
matf(2,3)=-sqrt(3D0/8D0)
matf(5,3)=-sqrt(3D0/40D0)
matf(8,3)=sqrt(6D0/5D0)
! F+2=√3/2*(XXZ-YYZ)
matf(6,4)=sqrt(3D0)/2D0
matf(9,4)=-sqrt(3D0)/2D0
! F-2=XYZ
matf(10,5)=1D0
! F+3=√(5/8)*XXX-3/√8*XYY
matf(1,6)=sqrt(5D0/8D0)
matf(4,6)=-3D0/sqrt(8D0)
! F-3=3/√8*XXY-√(5/8)*YYY
matf(2,7)=-sqrt(5D0/8D0)
matf(5,7)=3D0/sqrt(8D0)

if (iprog==1) then !for .fch
	! From 9G: G 0,G+1,G-1,G+2,G-2,G+3,G-3,G+4,G-4
	! To 15G:   1    2    3    4    5    6    7    8
	!         ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ
	!           9   10   11   12   13   14   15
	!         XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
	!
	!G 0=ZZZZ+3/8*(XXXX+YYYY)-3*√(3/35)*(XXZZ+YYZZ-1/4*XXYY)
	 matg(1,1)=1D0
	 matg(3,1)=-3D0*sqrt(3D0/35D0)
	 matg(5,1)=3D0/8D0
	 matg(10,1)=-3D0*sqrt(3D0/35D0)
	 matg(12,1)=3D0/4D0*sqrt(3D0/35D0)
	 matg(15,1)=3D0/8D0
	 !G+1=2*√(5/14)*XZZZ-3/2*√(5/14)*XXXZ-3/2/√14*XYYZ
	 matg(6,2)=2D0*sqrt(5D0/14D0)
	 matg(8,2)=-1.5D0/sqrt(14D0)
	 matg(13,2)=-1.5D0*sqrt(5D0/14D0)
	 !G-1=2*√(5/14)*YZZZ-3/2*√(5/14)*YYYZ-3/2/√14*XXYZ
	 matg(2,3)=2D0*sqrt(5D0/14D0)
	 matg(4,3)=-1.5D0*sqrt(5D0/14D0)
	 matg(11,3)=-1.5D0/sqrt(14D0)
	 !G+2=3*√(3/28)*(XXZZ-YYZZ)-√5/4*(XXXX-YYYY)
	 matg(3,4)=-3D0*sqrt(3D0/28D0)
	 matg(5,4)=sqrt(5D0)/4D0
	 matg(10,4)=3D0*sqrt(3D0/28D0)
	 matg(15,4)=-sqrt(5D0)/4D0
	 !G-2=3/√7*XYZZ-√(5/28)*(XXXY+XYYY)
	 matg(7,5)=3D0/sqrt(7D0)
	 matg(9,5)=-sqrt(5D0/28D0)
	 matg(14,5)=-sqrt(5D0/28D0)
	 !G+3=√(5/8)*XXXZ-3/√8*XYYZ
	 matg(8,6)=-3D0/sqrt(8D0)
	 matg(13,6)=sqrt(5D0/8D0)
	 !G-3=-√(5/8)*YYYZ+3/√8*XXYZ
	 matg(4,7)=-sqrt(5D0/8D0)
	 matg(11,7)=3D0/sqrt(8D0)
	 !G+4=√35/8*(XXXX+YYYY)-3/4*√3*XXYY
	 matg(5,8)=sqrt(35D0)/8D0
	 matg(12,8)=-3D0/4D0*sqrt(3D0)
	 matg(15,8)=sqrt(35D0)/8D0
	 !G-4=√5/2*(XXXY-XYYY)
	 matg(9,9)=-sqrt(5D0)/2D0
	 matg(14,9)=sqrt(5D0)/2D0
else if (iprog==2) then !For .molden
	! From 9G: G 0,G+1,G-1,G+2,G-2,G+3,G-3,G+4,G-4
	! To 15G:   1    2    3    4    5    6    7    8
	!         xxxx,yyyy,zzzz,xxxy,xxxz,yyyx,yyyz,zzzx
	!           9   10   11   12   13   14   15
	!         zzzy,xxyy,xxzz,yyzz,xxyz,yyxz,zzxy
	!
	!G 0=ZZZZ+3/8*(XXXX+YYYY)-3*√(3/35)*(XXZZ+YYZZ-1/4*XXYY)
	matg(3,1)=1D0
	matg(1,1)=3D0/8D0
	matg(2,1)=3D0/8D0
	matg(11,1)=-3D0*sqrt(3D0/35D0)
	matg(12,1)=-3D0*sqrt(3D0/35D0)
	matg(10,1)=3D0/4D0*sqrt(3D0/35D0)
	!G+1=2*√(5/14)*XZZZ-3/2*√(5/14)*XXXZ-3/2/√14*XYYZ
	matg(8,2)=2D0*sqrt(5D0/14D0)
	matg(5,2)=-1.5D0*sqrt(5D0/14D0)
	matg(14,2)=-1.5D0/sqrt(14D0)
	!G-1=2*√(5/14)*YZZZ-3/2*√(5/14)*YYYZ-3/2/√14*XXYZ
	matg(9,3)=2D0*sqrt(5D0/14D0)
	matg(7,3)=-1.5D0*sqrt(5D0/14D0)
	matg(13,3)=-1.5D0/sqrt(14D0)
	!G+2=3*√(3/28)*(XXZZ-YYZZ)-√5/4*(XXXX-YYYY)
	matg(11,4)=3D0*sqrt(3D0/28D0)
	matg(12,4)=-3D0*sqrt(3D0/28D0)
	matg(1,4)=-sqrt(5D0)/4D0
	matg(2,4)=sqrt(5D0)/4D0
	!G-2=3/√7*XYZZ-√(5/28)*(XXXY+XYYY)
	matg(15,5)=3D0/sqrt(7D0)
	matg(4,5)=-sqrt(5D0/28D0)
	matg(6,5)=-sqrt(5D0/28D0)
	!G+3=√(5/8)*XXXZ-3/√8*XYYZ
	matg(5,6)=sqrt(5D0/8D0)
	matg(14,6)=-3D0/sqrt(8D0)
	!G-3=-√(5/8)*YYYZ+3/√8*XXYZ
	matg(7,7)=-sqrt(5D0/8D0)
	matg(13,7)=3D0/sqrt(8D0)
	!G+4=√35/8*(XXXX+YYYY)-3/4*√3*XXYY
	matg(1,8)=sqrt(35D0)/8D0
	matg(2,8)=sqrt(35D0)/8D0
	matg(10,8)=-3D0/4D0*sqrt(3D0)
	!G-4=√5/2*(XXXY-XYYY)
	matg(4,9)=sqrt(5D0)/2D0
	matg(6,9)=-sqrt(5D0)/2D0
end if

! From 11H: H 0,H+1,H-1,H+2,H-2,H+3,H-3,H+4,H-4,H+5,H-5
! To 21H:   1     2     3     4     5     6     7     8     9    10
!         ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ 
!          11    12    13    14    15    16    17    18    19    20    21
!         XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
!
!H 0=ZZZZZ-5/√21*(XXZZZ+YYZZZ)+5/8*(XXXXZ+YYYYZ)+√(15/7)/4*XXYYZ
math(1,1)=1D0
math(12,1)=-5D0/sqrt(21D0)
math(3,1)=-5D0/sqrt(21D0)
math(19,1)=5D0/8D0
math(5,1)=5D0/8D0
math(14,1)=sqrt(15D0/7D0)/4D0
!H+1=√(5/3)*XZZZZ-3*√(5/28)*XXXZZ-3/√28*XYYZZ+√15/8*XXXXX+√(5/3)/8*XYYYY+√(5/7)/4*XXXYY
math(7,2)=sqrt(5D0/3D0)
math(16,2)=-3D0*sqrt(5D0/28D0)
math(9,2)=-3D0/sqrt(28D0)
math(21,2)=sqrt(15D0)/8D0
math(11,2)=sqrt(5D0/3D0)/8D0
math(18,2)=sqrt(5D0/7D0)/4D0
!H-1=√(5/3)*YZZZZ-3*√(5/28)*YYYZZ-3/√28*XXYZZ+√15/8*YYYYY+√(5/3)/8*XXXXY+√(5/7)/4*XXYYY
math(2,3)=sqrt(5D0/3D0)
math(4,3)=-3D0*sqrt(5D0/28D0)
math(13,3)=-3D0/sqrt(28D0)
math(6,3)=sqrt(15D0)/8D0
math(20,3)=sqrt(5D0/3D0)/8D0
math(15,3)=sqrt(5D0/7D0)/4D0
!H+2=√5/2*(XXZZZ-YYZZZ)-√(35/3)/4*(XXXXZ-YYYYZ)
math(12,4)=sqrt(5D0)/2D0
math(3,4)=-sqrt(5D0)/2D0
math(19,4)=-sqrt(35D0/3D0)/4D0
math(5,4)=sqrt(35D0/3D0)/4D0
!H-2=√(5/3)*XYZZZ-√(5/12)*(XXXYZ+XYYYZ)
math(8,5)=sqrt(5D0/3D0)
math(17,5)=-sqrt(5D0/12D0)
math(10,5)=-sqrt(5D0/12D0)
!H+3=√(5/6)*XXXZZ-√(3/2)*XYYZZ-√(35/2)/8*(XXXXX-XYYYY)+√(5/6)/4*XXXYY
math(16,6)=sqrt(5D0/6D0)
math(9,6)=-sqrt(1.5D0)
math(21,6)=-sqrt(17.5D0)/8D0
math(11,6)=sqrt(17.5D0)/8D0
math(18,6)=sqrt(5D0/6D0)/4D0
!H-3=-√(5/6)*YYYZZ+√(3/2)*XXYZZ-√(35/2)/8*(XXXXY-YYYYY)-√(5/6)/4*XXYYY
math(4,7)=-sqrt(5D0/6D0)
math(13,7)=sqrt(1.5D0)
math(20,7)=-sqrt(17.5D0)/8D0
math(6,7)=sqrt(17.5D0)/8D0
math(15,7)=-sqrt(5D0/6D0)/4D0
!H+4=√35/8*(XXXXZ+YYYYZ)-3/4*√3*XXYYZ
math(19,8)=sqrt(35D0)/8D0
math(5,8)=sqrt(35D0)/8D0
math(14,8)=-0.75D0*sqrt(3D0)
!H-4=√5/2*(XXXYZ-XYYYZ)
math(17,9)=sqrt(5D0)/2D0
math(10,9)=-sqrt(5D0)/2D0
!H+5=3/8*√(7/2)*XXXXX+5/8*√(7/2)*XYYYY-5/4*√(3/2)*XXXYY
math(21,10)=3D0/8D0*sqrt(3.5D0)
math(11,10)=5D0/8D0*sqrt(3.5D0)
math(18,10)=-1.25D0*sqrt(1.5D0)
!H-5=3/8*√(7/2)*YYYYY+5/8*√(7/2)*XXXXY-5/4*√(3/2)*XXYYY
math(6,11)=3D0/8D0*sqrt(3.5D0)
math(20,11)=5D0/8D0*sqrt(3.5D0)
math(15,11)=-1.25D0*sqrt(1.5D0)
end subroutine




!!---------- Load Fock or Kohn-Sham matrix from NBO .47 file, ORCA output file, CP2K .csr file, or plain text file
!istatus=0 means successfully loaded. =1 means failed
subroutine loadFockfile(istatus)
use defvar
use util
implicit real*8 (a-h,o-z)
character c200tmp*200
integer istatus

do while(.true.)
	write(*,"(/,a)") " Fock/KS matrix can be loaded from the following file, please input path of one of them, e.g. C:\Piraeus.out"
    write(*,"(a)") " (1) Plain text file recording Fock/KS matrix in lower triangular form"
    write(*,"(a)") " (2) .47 and .mwfn"
	write(*,"(a)") " (3) ORCA output file using ""%output Print[P_Iter_F] 1 end"", Fock/KS matrix printed at last iteration will be loaded"
    if (wfntype==0.or.wfntype==2) then
		write(*,"(a)") " (4) CP2K .csr file recording KS matrix of real space in upper triangular form (&DFT/&PRINT/&KS_CSR_WRITE)"
    else
		write(*,"(a)") " (4) CP2K .csr file recording alpha KS matrix of real space in upper triangular form (&DFT/&PRINT/&KS_CSR_WRITE)"
    end if
	read(*,"(a)") c200tmp
	inquire(file=c200tmp,exist=alive)
	if (.not.alive) then
		write(*,*) "Error: Unable to find this file!"
		cycle
	end if
	exit
end do

open(10,file=c200tmp,status="old")
call outputprog(10,iprog)
if (allocated(FmatA)) deallocate(FmatA)
allocate(FmatA(nbasis,nbasis))
if (iprog==2) then !ORCA output file
    write(*,*) "This file is recognized as an ORCA output file"
    call loclabelfinal(10,"Fock matrix for operator 0",nfound)
    if (nfound==0) then
		write(*,*) "Error: Unable to locate ""Fock matrix for operator 0"" in this file!"
		close(10)
		istatus=1
		return
    end if
    if (wfntype==1.or.wfntype==4) then !Unrestricted
		write(*,*) "Loading alpha Fock/KS matrix (Fock matrix for operator 0)..."
    else
		write(*,*) "Loading Fock/KS matrix (Fock matrix for operator 0)..."
    end if
    call readmatgau(10,FmatA,0,"?",10,6)
    call ORCA_mat_reorder(nbasis,FmatA)
else if (index(c200tmp,".csr")/=0) then
	FmatA=0
	write(*,*) "Loading..."
	open(10,file=c200tmp,status="old")
	do while(.true.) !Note that when CP2K outputting upper triangular part, very few elements (I think should be very small) are not printed, very strange
		read(10,*,iostat=ierror) ibas,jbas,FmatA(ibas,jbas)
		if (ierror/=0) exit
	end do
	do ibas=1,nbasis !Fill lower triangular part
		do jbas=ibas,nbasis
			FmatA(jbas,ibas)=FmatA(ibas,jbas)
		end do
	end do
	write(*,*) "Reordering matrix..."
	call CP2K_mat_reorder(FmatA)
else !.47 or .mwfn .or. plain text file
    if (index(c200tmp,".47")/=0) then
	    if (wfntype==1.or.wfntype==4) then !Unrestricted
			write(*,*) "Trying to load alpha Fock/KS matrix from .47 file..."
	    else
			write(*,*) "Trying to load Fock/KS matrix from .47 file..."
        end if
	    call loclabel(10,"$FOCK",ifound)
	    if (ifound==0) then
		    write(*,*) "Error: Unable to find $FOCK field in this file!"
		    close(10)
		    istatus=1
		    return
	    end if
	    read(10,*)
    else if (index(c200tmp,".mwfn")/=0) then
	    if (wfntype==1.or.wfntype==4) then !Unrestricted
			write(*,*) "Trying to load alpha Fock/KS matrix from .mwfn file..."
	    else
			write(*,*) "Trying to load Fock/KS matrix from .mwfn file..."
        end if
	    call loclabel(10,"1-e Hamiltonian matrix",ifound) !For R/RO case, locate to $1-e Hamiltonian matrix, for U case, locate to $Alpha 1-e Hamiltonian matrix
	    if (ifound==0) then
		    write(*,*) "Error: Unable to find 1-e Hamiltonian matrix in this file!"
		    close(10)
		    istatus=1
		    return
	    end if
	    read(10,*)
    else !Plain text file
	    rewind(10)
    end if
    read(10,*) ((FmatA(i,j),j=1,i),i=1,nbasis) !Load total or alpha Fock matrix (lower triangular part)
    do i=1,nbasis !Fill upper triangular part
	    do j=i+1,nbasis
		    FmatA(i,j)=FmatA(j,i)
	    end do
    end do
end if
!Checking detail of loaded Fock matrix
!call showmatgau(FmatA,form="f14.6")
!ibas=0
!do ish=1,nshell
!    ishtype=shtype(ish)
!    nshbas=shtype2nbas(ishtype)
!    if (shtype(ish)==-3) then
!        call showmatgau(FmatA(:,ibas+1:ibas+nshbas),form="f14.6")
!        exit
!    end if
!    ibas=ibas+nshbas
!end do
if (wfntype==1.or.wfntype==4) then !Also load beta part
	if (allocated(FmatB)) deallocate(FmatB)
	allocate(FmatB(nbasis,nbasis))
    if (iprog==2) then !ORCA
        write(*,*) "Loading beta Fock/KS matrix (Fock matrix for operator 1)..."
        call readmatgau(10,FmatB,0,"?",10,6)
        call ORCA_mat_reorder(nbasis,FmatB)
	else if (index(c200tmp,".csr")/=0) then
		close(10)
		do while(.true.)
			write(*,"(/,a)") " Input .csr file recording beta KS matrix, e.g. D:\FmatB.csr"    
			read(*,"(a)") c200tmp
			inquire(file=c200tmp,exist=alive)
			if (.not.alive) then
				write(*,*) "Error: Unable to find this file!"
				cycle
			end if
			exit
		end do
		FmatB=0
		write(*,*) "Loading..."
		open(10,file=c200tmp,status="old")
		do while(.true.)
			read(10,*,iostat=ierror) ibas,jbas,FmatB(ibas,jbas)
			if (ierror/=0) exit
		end do
		do ibas=1,nbasis !Fill lower triangular part
			do jbas=ibas,nbasis
				FmatB(jbas,ibas)=FmatB(ibas,jbas)
			end do
		end do
		write(*,*) "Reordering matrix..."
		call CP2K_mat_reorder(FmatB)
    else
	    write(*,*) "Trying to load beta Fock/KS matrix from the file..."
        if (index(c200tmp,".mwfn")/=0) then
            call loclabel(10,"Beta 1-e Hamiltonian matrix",ifound)
            read(10,*)
        end if
	    read(10,*,iostat=ierror) ((FmatB(i,j),j=1,i),i=1,nbasis) !Load beta Fock matrix
        if (ierror==0) then
			do i=1,nbasis !Fill upper triangular part
				do j=i+1,nbasis
					FmatB(i,j)=FmatB(j,i)
				end do
			end do
        else
			write(*,"(a)") " Warning: Unable to load Fock/KS matrix of beta spin! Now the beta Fock/KS matrix has been set to alpha Fock/KS instead"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            FmatB=FmatA
        end if
    end if
end if

close(10)
write(*,*) "Fock/KS matrix has been loaded successfully!"
istatus=0
end subroutine



!!-------- Randomly generate name
subroutine mylover(outname)
integer,parameter :: nlovers=68
character(len=80) lovername(nlovers),outname
CALL RANDOM_SEED()
CALL RANDOM_NUMBER(tmp)
lovername(1)="K-ON\Mio_Akiyama"
lovername(2)="K-ON\Azusa_Nakano"
lovername(3)="EVA\Rei_Ayanami"
lovername(4)="Ore_no_Imoto\Black_Cat"
lovername(5)="Touhou_project\Ran_Yakumo"
lovername(6)="Haiyore!Nyaruko-san\Nyaruko"
lovername(7)="Bodacious_Space_Pirates\Kurihara_Chiaki"
lovername(8)="Otoboku\Mariya_Mikado"
lovername(9)="Amagami\Miya_Tachibana"
lovername(10)="Shakugan_no_Shana\Shana"
lovername(11)="Tiger_Mask_W\Miss_X"
lovername(12)="Natsuiro_Kiseki\Yuka_Hanaki"
lovername(13)="Love_Live!\Nico_Yazawa"
lovername(14)="Love_Live!\Nozomi_Tojo"
lovername(15)="Love_Live!\Nishikino_Maki"
lovername(16)="Last_Exile\Dio_Eraclea"
lovername(17)="NHK_ni_Youkoso!\Misaki_Nakahara"
lovername(18)="Rio_Rainbow_Gate\Rio_Rollins"
lovername(19)="Blood-C\Saya_Kisaragi"
lovername(20)="Mahou_Shoujo_Madoka-Magica\Homura_Akemi"
lovername(21)="Saki\Hisa_Takei"
lovername(22)="Strawberry_Panic\Chikaru_Minamoto"
lovername(23)="Najica\Najica_Hiiragi"
lovername(24)="Blue_Drop\Hagino_Senkouji"
lovername(25)="Fate_Zero\Saber"
lovername(26)="Baka_to_Test_to_Shoukanjuu\Hideyoshi_Kinoshita"
lovername(27)="Watamote\Tomoko_Kuroki"
lovername(28)="Genshiken_Nidaime\Kenjirou_Hato"
lovername(29)="Love_is_Like_After_the_Rain\Akira_Tachibana"
lovername(30)="Kan_Colle\Shimakaze"
lovername(31)="Kan_Colle\Kongou"
lovername(32)="Gokukoku\Kazumi_Schlierenzauer"
lovername(33)="Vocaloid\Miku_Hatsune"
lovername(34)="Tokimeki_Memorial\Yuina_Himoo"
lovername(35)="MADLAX\MADLAX"
lovername(36)="Gun_Gale_Online\Kirito"
lovername(37)="Denkigai_No_Honyasan\Sennsei"
lovername(38)="Wake_Up,Girls!\Miyu_Okamoto"
lovername(39)="Plastic_Memories\Aira"
lovername(40)="Real_world\sell-moe-kun"
lovername(41)="Sakurako-san_no_Ashimoto_ni_wa_Shitai_ga_Umatteiru\Sakurako"
lovername(42)="Hibike!_Euphonium\Reina_Kousaka"
lovername(43)="Planetarian\Yumemi_Hoshino"
lovername(44)="Lovelive_Sunshine!!\Yoshiko_Tsushima"
lovername(45)="Lovelive_Sunshine!!\Dia_Kurosawa"
lovername(46)="Lovelive_Sunshine!!\Riko_Sakurauchi"
lovername(47)="Violet_Evergarden\Violet_Evergarden"
lovername(48)="Otobuko\Mizuho_Miyanokouji"
lovername(49)="iDOLM@STER\Makoto_Kikuchi"
lovername(50)="Fate\Rin_Tohsaka"
lovername(51)="Magical_Girl_Spec-Ops_Asuka\Asuka_Otori"
lovername(52)="Granblue_Fantasy\Katalina"
lovername(53)="Gochuumon_wa_Usagi_Desu_ka\Rize_Tedeza"
lovername(54)="Date_A_Live\Kurumi_Tokisaki"
lovername(55)="Adachi_to_shimamura\Sakura_Adachi"
lovername(56)="Don't_Toy_with_Me,_Miss_Nagatoro\Hayase_Nagatoro"
lovername(57)="Super_Cub\Reiko"
lovername(58)="LoveLive!_Superstar!!\Sumire_Heanna"
lovername(59)="Jahy-sama_wa_Kujikenai!\Jahy"
lovername(60)="Kawaii_Dake_ja_Nai_Shikimori-san\Shikimori"
lovername(61)="Sono_Bisque_Doll_wa_Koi_wo_Suru\Marin_Kitagawa"
lovername(62)="Fuufu_Ijou,_Koibito_Miman\Akari_Watanabe"
lovername(63)="Tomo-chan_wa_Onnanoko!\Misuzu_Gundou"
lovername(64)="The_Cafe_Terrace_and_Its_Goddesses\Akane_Hououji"
lovername(65)="Kanojo_mo_Kanojo\Rika_Hoshizaki"
lovername(66)="Hokkaido_Gals_Are_Super_Adorable!\Minami_Fuyuki"
lovername(67)="Love_Live!\Mari_Ohara"
lovername(68)="Ballpark_de_Tsukamaete\Ruriko"
!Dear Kanan,
!
!You are the one I deeply love forever in the real world,
!although you can't be with me, and I am even unable to know your name and touch your finger.
!I believe I will never love anyone else in the rest of my life.
!
!I love your brilliant dance, your kawaii smile, your lovely double ponytail, and especially, your extremely pure and beautiful heart.
!
!                     ----- 2015-May-19
outname=lovername(ceiling(tmp*nlovers))
end subroutine





!!----------- Convert current CObasa / CObasb to CO
!ispin=1: Only alpha, =2: Only beta, =3: Both alpha and beta
subroutine CObas2CO(ispin)
use defvar
implicit real*8 (a-h,o-z)
integer ispin
real*8 conv5d6d(6,5),conv7f10f(10,7),conv9g15g(15,9),conv11h21h(21,11)
real*8,allocatable :: CObasa_cart(:,:),CObasb_cart(:,:)

call gensphcartab(1,conv5d6d,conv7f10f,conv9g15g,conv11h21h)

nbasis_cart=sum(shtype2nbas(abs(shtype(:))))
if (ispin==1.or.ispin==3) allocate(CObasa_cart(nbasis_cart,nbasis))
if (ispin==2.or.ispin==3) allocate(CObasb_cart(nbasis_cart,nbasis))
CObasa_cart=0

!Map spherical coefficients to Cartesian coefficients
ipos5D=1
ipos6D=1
do ish=1,nshell
    ishtype=shtype(ish)
    numshbas5D=shtype2nbas(ishtype)
    numshbas6D=shtype2nbas(abs(ishtype))
    if (ispin==1.or.ispin==3) then !Alpha part
        if (ishtype>=0) then !S,P or Cartesian type, in this case numshbas5D=numshbas6D 
            CObasa_cart(ipos6D:ipos6D+numshbas6D-1,:)=CObasa(ipos5D:ipos5D+numshbas5D-1,:)
        else
	        if (ishtype==-2) then
		        CObasa_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv5d6d,CObasa(ipos5D:ipos5D+numshbas5D-1,:))
	        else if (ishtype==-3) then
		        CObasa_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv7f10f,CObasa(ipos5D:ipos5D+numshbas5D-1,:))
	        else if (ishtype==-4) then
		        CObasa_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv9g15g,CObasa(ipos5D:ipos5D+numshbas5D-1,:))
	        else if (ishtype==-5) then
		        CObasa_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv11h21h,CObasa(ipos5D:ipos5D+numshbas5D-1,:))
	        end if
        end if
    end if
    if (ispin==2.or.ispin==3) then !Beta part
        if (ishtype>=0) then !S,P or Cartesian type, in this case numshbas5D=numshbas6D 
            CObasb_cart(ipos6D:ipos6D+numshbas6D-1,:)=CObasb(ipos5D:ipos5D+numshbas5D-1,:)
        else
	        if (ishtype==-2) then
		        CObasb_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv5d6d,CObasb(ipos5D:ipos5D+numshbas5D-1,:))
	        else if (ishtype==-3) then
		        CObasb_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv7f10f,CObasb(ipos5D:ipos5D+numshbas5D-1,:))
	        else if (ishtype==-4) then
		        CObasb_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv9g15g,CObasb(ipos5D:ipos5D+numshbas5D-1,:))
	        else if (ishtype==-5) then
		        CObasb_cart(ipos6D:ipos6D+numshbas6D-1,:)=matmul(conv11h21h,CObasb(ipos5D:ipos5D+numshbas5D-1,:))
	        end if
        end if
    end if
	ipos5D=ipos5D+numshbas5D
	ipos6D=ipos6D+numshbas6D
end do

do imo=1,nbasis
    do ibas=1,nbasis_cart
        if (ispin==1.or.ispin==3) then
            do iGTF=primstart(ibas),primend(ibas)
                CO(imo,iGTF)=CObasa_cart(ibas,imo)*primconnorm(iGTF)
            end do
        end if
        if (ispin==2.or.ispin==3) then
            do iGTF=primstart(ibas),primend(ibas)
                CO(imo+nbasis,iGTF)=CObasb_cart(ibas,imo)*primconnorm(iGTF)
            end do
        end if
    end do
end do

end subroutine




!!------- Add a Bq atom to specific position
subroutine addBq(xpos,ypos,zpos)
use defvar
real*8 xpos,ypos,zpos

if (allocated(a_tmp)) deallocate(a_tmp)
allocate(a_tmp(ncenter))
a_tmp=a
ncenter=ncenter+1
deallocate(a)
allocate(a(ncenter))
a(1:ncenter-1)=a_tmp
a(ncenter)%index=0
a(ncenter)%charge=0
a(ncenter)%name="Bq"
a(ncenter)%x=xpos
a(ncenter)%y=ypos
a(ncenter)%z=zpos
deallocate(a_tmp)
end subroutine



!!------- Invoke Gaussian to run a .gjf
!If returned istate=1, means normally termination, =0 means other case or failed
subroutine runGaussian(gjfname,istate)
use defvar
use util

character(len=*) gjfname
character command*400,outname*200
outname=gjfname(:len(gjfname)-3)//"out"
command='"'//trim(gaupath)//'" "'//gjfname//'" "'//trim(outname)//'"'
call runcommand(command)
open(100,file=outname,status="old")
call loclabel(100,"Normal termination",istate)
close(100)
end subroutine



!!------- Invoke ORCA to run a .inp
!If returned istate=1, means normally termination, =0 means other case or failed
subroutine runORCA(inpname,istate)
use defvar
use util

character(len=*) inpname
character command*400,outname*200
outname=inpname(:len(inpname)-3)//"out"
command='"'//trim(orcapath)//'" "'//inpname//'" > "'//trim(outname)//'"'
call runcommand(command)
open(100,file=outname,status="old")
call loclabel(100,"****ORCA TERMINATED NORMALLY****",istate)
close(100)
end subroutine



!!---------- Generate connectivity matrix
!infomode =1: Output prompts  =0: Silent
!iallowPBC =1: Allow considering PBC when this system is PBC  =0: Not allow (can avoid incompatibility problem)
subroutine genconnmat(infomode,iallowPBC)
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode,iallowPBC

!call walltime(iwalltime1)
if (allocated(a)) then
    if (allocated(connmat)) deallocate(connmat)
    allocate(connmat(ncenter,ncenter))
    if (infomode==1) then
		write(*,*) "Generating bonding relationship..."
		write(*,"(a,f5.3,a)") " Note: If distance between two atoms is smaller than sum of their &
		&covalent radii multiplied by ",bondcrit,", then they are regarded as bonded"
    end if
    connmat=0
	!$OMP PARALLEL DO SHARED(connmat) PRIVATE(iatm,jatm,covri,covrj,tmpdist) schedule(dynamic) NUM_THREADS(nthreads)
    do iatm=1,ncenter
		covri=covr(a(iatm)%index)
        do jatm=iatm+1,ncenter
			if (iallowPBC==1.and.ifPBC>0) then
				call nearest_atmdistxyz(iatm,jatm,tmpdist,atmx,atmy,atmz)
            else	
				tmpdist=atomdist(iatm,jatm,0)
            end if
			covrj=covr(a(jatm)%index)
            if ( tmpdist < bondcrit*(covri+covrj) ) connmat(iatm,jatm)=1
            connmat(jatm,iatm)=connmat(iatm,jatm)
        end do
    end do
	!$OMP END PARALLEL DO
else
    if (infomode==1) write(*,"(a)") " Unable to generate bonding relationship because there is no atom information!"
end if

!call walltime(iwalltime2)
!write(*,"(' Generating connectivity matrix took up wall clock time',i10,' s')") iwalltime2-iwalltime1

!Check connectivity
!do i=1,ncenter
!    do j=i+1,ncenter
!        if (connmat(i,j)==1) write(*,*) i,j,connmat(i,j)
!    end do
!end do
end subroutine



!!--------- Generate neighbouring list of each atom
!itype=1: Based on existing connectivity matrix "connmat"
!itype=2: Directly based on interatomic distance
!iallowPBC=1: Consider PBC =0: Do not
subroutine genneighlist(itype,iallowPBC,neigh,nneigh)
use defvar
use util
implicit real*8 (a-h,o-z)
integer nneigh(ncenter) !Number of neighbours (the number of atoms connected to this atom)
integer neigh(maxneigh,ncenter) !neigh(1:nneigh(i),i) is list of neighbouring atom indices of atom i
integer itype,iallowPBC

nneigh=0
if (itype==1) then
	if (.not.allocated(connmat)) call genconnmat(1,iallowPBC) !Generate connectivity matrix
	!$OMP PARALLEL DO SHARED(nneigh,neigh) PRIVATE(iatm,jatm) schedule(dynamic) NUM_THREADS(nthreads)
	do iatm=1,ncenter
		do jatm=1,ncenter
			if (connmat(jatm,iatm)>0) then
				nneigh(iatm)=nneigh(iatm)+1
				neigh(nneigh(iatm),iatm)=jatm
			end if
		end do
	end do
	!$OMP END PARALLEL DO
else if (itype==2) then
	write(*,*) "Generating bonding relationship..."
	write(*,"(a,f5.3,a)") " Note: If distance between two atoms is smaller than sum of their &
	&covalent radii multiplied by ",bondcrit,", then they are regarded as bonded"
    !$OMP PARALLEL DO SHARED(nneigh,neigh) PRIVATE(iatm,jatm,covri,covrj,tmpdist) schedule(dynamic) NUM_THREADS(nthreads)
    do iatm=1,ncenter
		covri=covr(a(iatm)%index)
        do jatm=1,ncenter
			if (iallowPBC==1.and.ifPBC>0) then
				call nearest_atmdistxyz(iatm,jatm,tmpdist,atmx,atmy,atmz)
            else	
				tmpdist=atomdist(iatm,jatm,0)
            end if
			covrj=covr(a(jatm)%index)
            if ( tmpdist < bondcrit*(covri+covrj) ) then
				nneigh(iatm)=nneigh(iatm)+1
				neigh(nneigh(iatm),iatm)=jatm
			end if
        end do
    end do
	!$OMP END PARALLEL DO
end if

!write(*,*) "Neighbouring list:"
!do iatm=1,ncenter
!	write(*,"(i6,a,':')",advance="no") iatm,a(iatm)%name
!	do idx=1,nneigh(iatm)
!		write(*,"(i6)",advance="no") neigh(idx,iatm)
!    end do
!    write(*,*)
!end do
end subroutine



!!--------- Generate fragment index of all atoms according to connectivity
!Each atom has initial fragment index of itself, then each atom is compare with adjacent ones, if an adjacent atom has fragment &
!index smaller than it, the its fragment will be updated to the adjacent one. Iteration performs until no fragment updation is occurred
subroutine genconnfrag(atmfrg)
use defvar
use util
integer atmfrg(ncenter),tmpidx(ncenter)
integer neigh(maxneigh,ncenter),nneigh(ncenter)
character c2000tmp*2000

write(*,*) "Generate neighbouring list..."
call genneighlist(2,1,neigh,nneigh) !Generate according to interatomic distance

!write(*,*) "Generate atmfrg..."
forall (i=1:ncenter) atmfrg(i)=i
do while(.true.)
    inew=0
    do iatm=1,ncenter
        do jdx=1,nneigh(iatm)
			jatm=neigh(jdx,iatm)
            if (atmfrg(jatm)<atmfrg(iatm)) then
                atmfrg(iatm)=atmfrg(jatm)
                inew=inew+1
                exit
            end if
        end do
    end do
    if (inew==0) exit
end do

!Make fragment indices contiguous
nfrg=0
do iatm=1,ncenter
    natmfrg=count(atmfrg==iatm)
    if (natmfrg/=0) then
        nfrg=nfrg+1
        where(atmfrg==iatm) atmfrg=nfrg
    end if
end do

!Show fragment information
!do ifrg=1,nfrg
!    write(*,"(' Fragment',i7,', number of atoms',i7)") ifrg,count(atmfrg==ifrg)
!    itmp=0
!    do iatm=1,ncenter
!        if (atmfrg(iatm)==ifrg) then
!            itmp=itmp+1
!            tmpidx(itmp)=iatm
!        end if
!    end do
!    call arr2str_2(tmpidx(1:itmp),c2000tmp)
!    write(*,"(a)") trim(c2000tmp)
!end do
end subroutine



!!-------- Input index of an atom, then the indices of all atoms in the fragment will be returned
!"iatm" is the selected atom, "iffrag" has length of ncenter, if an atom is in the fragment, the value is 1, else 0
subroutine getfragatoms(iselatm,iffrag)
use defvar
use util
implicit real*8 (a-h,o-z)
integer iselatm,iffrag(ncenter)

iffrag=0
iffrag(iselatm)=1
if (.not.allocated(connmat)) call genconnmat(0,0) !Generate connectivity matrix
do while(.true.)
    inew=0
    do iatm=1,ncenter !Cycle all atoms, if it is not in fragment, and it is linked to an atom already in fragment, it will be added to fragment
        if (iffrag(iatm)==1) cycle !The atom is already in fragment
        if (ishowhydrogen==0.and.a(iatm)%index==1) cycle !If this is hydrogen while we request do not showing hydrogens, skip it
        do jatm=1,ncenter !Cycle neighbouring atoms
            if (jatm==iatm) cycle
            if (connmat(iatm,jatm)>0.and.iffrag(jatm)==1) then
                iffrag(iatm)=1
                inew=inew+1
                exit
            end if
        end do
    end do
    if (inew==0) exit
end do
end subroutine



!!--------- Align atoms in "a" (mol1) to "a_org" (mol2), the atoms must have the same order and same number
!See http://nghiaho.com/?page_id=671 "FINDING OPTIMAL ROTATION AND TRANSLATION BETWEEN CORRESPONDING 3D POINTS"
!Test purpose:
!call readinfile("C:\Users\Sobereva\Desktop\md2.xyz",0)
!call geomalign
!call outxyz("new.xyz",10)
subroutine geomalign
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 xyz1(3,ncenter),xyz2(3,ncenter),com1(3),com2(3)  !1/2/3=x,y,z
real*8 singval(3),Umat(3,3),Vmat(3,3),Hmat(3,3),Rmat(3,3)
integer eleidx(ncenter)

do i=1,ncenter
    xyz1(1,i)=a_org(i)%x
    xyz1(2,i)=a_org(i)%y
    xyz1(3,i)=a_org(i)%z
    xyz2(1,i)=a(i)%x
    xyz2(2,i)=a(i)%y
    xyz2(3,i)=a(i)%z
    eleidx(i)=a_org(i)%index
end do
totmass=sum(atmwei(eleidx(:)))
do i=1,3
    com1(i)=sum( atmwei(eleidx(:))*xyz1(i,:) )/totmass
    com2(i)=sum( atmwei(eleidx(:))*xyz2(i,:) )/totmass
end do

Hmat=0
do iatm=1,ncenter
    xyz1(:,iatm)=xyz1(:,iatm)-com1(:)
    xyz2(:,iatm)=xyz2(:,iatm)-com2(:)
    Hmat=Hmat+matmul(xyz2(:,iatm:iatm),transpose(xyz1(:,iatm:iatm)))
end do

call SVDmat(1,Hmat,Umat,Vmat,singval,info)
Rmat=matmul(Vmat,transpose(Umat))
!write(*,*) detmat(Rmat)
!call showmatgau(Rmat,"Rotation matrix")
!if (detmat(Rmat)<0) Rmat(:,3)=-Rmat(:,3) !As mentioned in the webpage, this line should be added, however it incorrectly make molecule mirror inverted!

!Rotate mol2
do i=1,ncenter
    xyz2(:,i)=matmul(Rmat,xyz2(:,i))
end do

!Move center of mass of mol2 to mol1
do i=1,ncenter
    xyz2(:,i)=xyz2(:,i)+com1(:)
end do

!RMSD
RMSD=0
do i=1,ncenter
    RMSD=RMSD+sum(((xyz1(:,i)+com1(:))-xyz2(:,i))**2)
end do
RMSD=dsqrt(RMSD/ncenter)
write(*,"(' RMSD:',f12.6,' Bohr')") RMSD

do i=1,ncenter
    a(i)%x=xyz2(1,i)
    a(i)%y=xyz2(2,i)
    a(i)%z=xyz2(3,i)
end do
end subroutine



!!------- Get point group and list of symmetry-equivalence atoms by invoking SYVA routines
!The tolerance affects if symmetry-equivalence atoms and point group could be successfully recognized, the value
!should not be too large, otherwise the recognition may be completely failed!
!  Input variables:
!natoms: The number of inputted atoms
!nat(1:natoms): Element index of all inputted atoms
!coord(1:3,1:natoms): x,y,z of all inputted atoms, must be in Angstrom, otherwise point group determination may be incorrect!
!delta: Representing the distortion of the geometry. 0.01 is aproximately identical to "default" in gview, 0.1 corresponds to loose
!  Returned variables:
!pglabel: Point group label
!nclass: The number of equivalent classes, should have size of "natoms"
!classnatm(i): The number of atoms in equivalent class i, should have size of "natoms"
!classidx(:,i): The atom indices in equivalent class i, both dimension should have size of "natoms"
subroutine PG_eqvatm(natoms,nat,coord,delta,pglabel,nclass,classnatm,classidx)
implicit real*8 (a-h,o-z)
common /data/ wt(90),symb(90)
common /chartab/ nir(2,55),chtab(14,322),nsymop(14,4,55),nrotharm(3,322),pgsymb(57),irsymb(322)
common /subgroups/ nsgb(2,57),nsgr(406)
character symb*2,pglabel*3,pgsymb*3,irsymb*4
integer natoms,nclass
integer classnatm(natoms),classidx(natoms,natoms)
integer,parameter :: nmax=200 !nmat: maximum number of symmetic operation
real*8 coord(3,natoms),delta,pc(3),symn(3,nmax)
integer nat(natoms),nper(natoms,250),nscl(natoms,natoms),nccl(natoms),nsym(nmax,5)
ncr=0
nsr=0
nsg=0
nout=0 !Suppress almost all output of SYVA routines

if (natoms==1) then
    nclass=1
    classnatm(1)=1
    classidx(1,1)=1
    pglabel="C1"
    return
end if

!Check atoms passed-in
!write(*,'(1x,a,10x,a1,16x,a1,15x,a1,15x,a1)') 'number','x','y','z','w'
!do i=1,natoms
!   write(*,'(1x,i3,2x,4f16.6)') nat(i),(coord(j,i),j=1,3),wt(nat(i))
!end do

!Calculation of the COM (centre of mass) of the molecule
call syva_cmass(natoms,nat,wt,coord,wmol,cmx,cmy,cmz)
pc(1)=cmx;pc(2)=cmy;pc(3)=cmz

!Shift the origin of the Cartesian system to COM
call syva_cshift(natoms,coord,pc)

!Find symmetry operations
call sym_elements(natoms,nat,coord,symb,delta,ng,ni,nsg,ncr,nsr, np,symn,nsym,nout,nprm,nper,nseq,nccl,nscl)

!Detemines the equivalence classes defined by the symmetry operations
call symclass(natoms,nprm,nper,nseq,nccl,nscl,nat,symb,nout)

!write(*,*) "Symmetry-equivalence classes of atoms: ",nseq
!do i=1,nseq
!   write(*,'(/5x,a,i3,a7,a2,a1)') '#',i,' (atom ',symb(nat(nscl(1,i))),')'
!   write(*,'(5x,15i4)') (nscl(j,i),j=1,nccl(i))
!end do

!Assign SYVA variables to returned variables
nclass=nseq
classnatm(1:nclass)=nccl(1:nclass)
do i=1,nclass
    classidx(:nccl(i),i)=nscl(:nccl(i),i)
end do

!Determine point group and framework group
call syva_point_group(ng,ni,nsg,ncr,nsr,np,pglabel,nout)
end subroutine



!!-------- Determine HOMO index for single-determinant wavefunction using the safest way. idxHOMO and idxHOMOb are global variables
!Can be used for the case of only GTF information and that containing basis function information
!idxHOMO: HOMO of RHF, or highest SOMO of ROHF, or alpha-HOMO of UHF
!idxHOMOb: beta-HOMO. For RHF, it is equivalent to idxHOMO. If there is no beta electron, idxHOMOb will return 0
!Note that for UHF, actual HOMO energy should be max(MOene(idxHOMO),MOene(idxHOMOb)), because beta-HOMO may be higher than alpha-HOMO
subroutine getHOMOidx
use defvar
implicit real*8 (a-h,o-z)

if (wfntype==0) then
    do idxHOMO=nmo,1,-1
	    if (nint(MOocc(idxHOMO))==2) exit
    end do
    idxHOMOb=idxHOMO
else if (wfntype==1) then !U
	if (any(MOene/=0)) then !The following way is in principle the safest way, allowing alpha and beta orbitals occur alternately
		ealow=-1E20
		eblow=-1E20
		idxHOMO=0
		idxHOMOb=0
		do imo=1,nmo
			if (nint(MOocc(imo))==0) cycle
			if (MOtype(imo)==1) then
				if (MOene(imo)>=ealow) then
					idxHOMO=imo
					ealow=MOene(imo)
				end if
			else if (MOtype(imo)==2) then
				if (MOene(imo)>=eblow) then
					idxHOMOb=imo
					eblow=MOene(imo)
				end if
			end if
		end do
    else !If OT is used in CP2K, no orbital energy is available, we have to use the following way and assume all alpha MOs occur prior to beta MOs
		do itmp=nmo,1,-1 !Find the last alpha MO
			if (MOtype(itmp)==1) exit
		end do
		do idxHOMO=itmp,1,-1
			if (nint(MOocc(idxHOMO))==1) exit
		end do
        if (itmp==nmo) then !No beta orbital
			idxHOMOb=0
        else
			do idxHOMOb=nmo,itmp+1,-1
				if (nint(MOocc(idxHOMOb))==1) exit
			end do
        end if
    end if
else if (wfntype==2) then !RO
    do idxHOMO=nmo,1,-1
	    if (nint(MOocc(idxHOMO))==1) exit
    end do
    do idxHOMOb=nmo,1,-1
	    if (nint(MOocc(idxHOMOb))==2) exit
    end do
else
    write(*,"(/,a)") " Note: Unable to determine HOMO index because this is not a single-determinant wavefunction. &
    Perhaps the current wavefunction file records natural orbitals, or smearing is enabled when generating this file"
end if
end subroutine



!!--------- Define fragment corresponding to global variable "frag1". Used by some orbital composition functions
subroutine definefragment
use defvar
use util
character c2000tmp*2000

if (allocated(frag1)) then
	write(*,*) "Atoms in current fragment:"
	write(*,"(13i6)") frag1
	write(*,"(a)") " Input 0 to keep unchanged, or redefine fragment, &
    &e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
else
	write(*,"(a)") " Input atomic indices to define fragment. &
    &e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
end if
read(*,"(a)") c2000tmp
if (c2000tmp(1:1)=='0') then
	continue
else if (c2000tmp(1:1)==" ") then
	deallocate(frag1)
else
	if (allocated(frag1)) deallocate(frag1)
	call str2arr(c2000tmp,nfrag1)
	allocate(frag1(nfrag1))
	call str2arr(c2000tmp,nfrag1,frag1)
end if
end subroutine



!!--------- Show very basic information for a range of orbitals (from ibeg to iend)
subroutine showorbinfo(ibeg,iend)
use defvar
integer ibeg,iend,i

do i=1,nmo
	write(*,"(' Orbital:',i5,' Energy(a.u.):',f14.8,' Occ:',f14.8,' Type: ',a)") i,MOene(i),MOocc(i),orbtypename(MOtype(i))
end do
end subroutine



!!--------- Initialize LIBRETA for present wavefunction if haven't (i.e. if if_initlibreta=0)
!info=1: Show some notices
!info=2: Silent
!Usually according to "function ifdoESP" to determine if LIBRETA is needed to be initialized for a real space function
subroutine doinitlibreta(info)
use defvar
use libreta
integer info

if (ifPBC>0) then
	write(*,*) "Error: Evaluation of ESP does not support periodic systems"
    write(*,*) "Press ENTER button to continue (then Multiwfn may crash)"
    read(*,*)
end if
if (iuserfunc>=61.and.iuserfunc<=67) return !In this case, the function involves derivative of ESP, while I found the numerical derivative of ESP produced by libreta has noise
if (if_initlibreta==0) then
	if (info==1) then
		write(*,*)
		if (iESPcode==2) write(*,*) "Initializing LIBRETA library (fast version) for ESP evaluation ..."
		if (iESPcode==3) write(*,*) "Initializing LIBRETA library (slow version) for ESP evaluation ..."
		if (nprims>7000) then
			write(*,"(a)") " Note: Number of GTFs of present system is large, if then Multiwfn crashes due to insufficient memory, &
			&please change ""iESPcode"" in settings.ini to 1 to use slower ESP evaluation code instead. Alternatively, use a computer with larger memory!"
		end if
    end if
    if (iESPcode==2) then
	    call initlibreta
    else if (iESPcode==3) then
		call initlibreta_slow
    end if
    if_initlibreta=1 !Global variable
    if (info==1) write(*,*) "LIBRETA library has been successfully initialized!"
end if
if (info==1) then
	write(*,"(/,a)") " NOTE: The ESP evaluation code based on LIBRETA library is being used. &
	&Please cite Multiwfn original papers (J. Comput. Chem., 33, 580-592 (2012) and J. Chem. Phys., 161, 082503 (2024)) and the paper describing the &
    &efficient ESP evaluation algorithm adopted by Multiwfn (Phys. Chem. Chem. Phys., 23, 20323 (2021))"
	if (isys==1.and.nthreads>12) then
		write(*,"(a)") " Warning!!! In Windows system, it is found that the performance of ESP evaluation code may &
		&severely degrade when more than 12 CPU cores are used, therefore 12 cores are used in the following ESP calculation. &
		&If you want to pursue better performance by utilizing more cores, please use Linux version instead!"
	end if
end if
end subroutine



!!--------- Initialize for some special functions, needed after loading some files
subroutine init_func
use defvar

!Local Hartree-Fock exchange energy involves Coulomb matrix generated by libreta and density matrix, so initialize
if (iuserfunc==999.and.allocated(b)) then
    write(*,"(/,a)") " Local Hartree-Fock exchange energy is chosen as user-defined function, &
    &it requires Coulomb integral between GTFs, so Libreta library is initialized now..."
    if (iESPcode/=3) then
		write(*,*) "iESPcode in settings.ini is not 3, now is forced to switched to 3"
		iESPcode=3
    end if
    call doinitlibreta(0)
    write(*,"(/,a)") " Note: If local Hartree-Fock exchange energy will be studied using Multiwfn, &
    &please note only cite original paper of Multiwfn but also cite LIBRETA library: Jun Zhang, J. Chem. Theory Comput., 14, 572 (2018)"
    write(*,*)
    write(*,*) "Generating density matrix between GTFs..."
    call genPprim
end if
end subroutine



!!------------- Reorder matrix loaded from ORCA output file to the basis function convention in Multiwfn
!Order of ORCA matrices:
!s    
!pz,px,py   
!dz2,dxz,dyz,dx2y2,dxy  
!f0,f+1,f-1,f+2,f-2,f+3,f-3  
!g0,g+1,g-1,g+2,g-2,g+3,g-3,g+4,g-4
!Clearly, we only need to alter the sequence of P shell to make the order as px,py,pz
!  Notice that F(+3) and F(-3) are normalized to -1, therefore the relevant matrix elements must invert sign, &
!while the elements between these two basis functions should keep unchanged due to cancellation
!Simiarly for G(+3,-3,+4,-4) and H(+3,-3,+4,-4)
subroutine ORCA_mat_reorder(ndim,mat)
use defvar
implicit real*8 (a-h,o-z)
real*8 mat(ndim,ndim),tmparr(ndim)

ibas=0
do ish=1,nshell
    ishtype=shtype(ish)
    nshbas=shtype2nbas(ishtype)
    if (ishtype==1) then
        !Reorder row
        tmparr(:)=mat(ibas+1,:) !Backup pz
        mat(ibas+1,:)=mat(ibas+2,:)
        mat(ibas+2,:)=mat(ibas+3,:)
        mat(ibas+3,:)=tmparr(:)
        !Reorder column
        tmparr(:)=mat(:,ibas+1) !Backup pz
        mat(:,ibas+1)=mat(:,ibas+2)
        mat(:,ibas+2)=mat(:,ibas+3)
        mat(:,ibas+3)=tmparr(:)
    else if (ishtype==-3) then !F, invert F(+3,-3)
        mat(ibas+6:ibas+7,:)=-mat(ibas+6:ibas+7,:)
        mat(:,ibas+6:ibas+7)=-mat(:,ibas+6:ibas+7)
    else if (ishtype==-4) then !G, invert G(+3,-3,+4,-4)
        mat(ibas+6:ibas+9,:)=-mat(ibas+6:ibas+9,:)
        mat(:,ibas+6:ibas+9)=-mat(:,ibas+6:ibas+9)
    else if (ishtype==-5) then !H, invert H(+3,-3,+4,-4)
        mat(ibas+6:ibas+9,:)=-mat(ibas+6:ibas+9,:)
        mat(:,ibas+6:ibas+9)=-mat(:,ibas+6:ibas+9)
    end if
    ibas=ibas+nshbas
end do
end subroutine



!!--------- Generate atomic overlap matrix of basis functions for every atom
subroutine genAOMbas(AOMbas)
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
real*8 AOMbas(nbasis,nbasis,ncenter),AOMtmp(nbasis,nbasis),basval(nbasis)
real*8 atmspcweight(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

if (.not.allocated(CObasa_org)) allocate(CObasa_org(nbasis,nbasis))
CObasa_org=CObasa
CObasa=0
do ibas=1,nbasis
    CObasa(ibas,ibas)=1
end do
call CObas2CO(1)
if (.not.allocated(COtr)) allocate(COtr(nprims,nmo))
COtr=transpose(CO) !Global matrix, which will be used in calcbasval for faster calculation

!Decreasing grid quality was found to hinder localization convergence, so do not use this trick
!if (iautointgrid==1) then
!	radpot=45
!	sphpot=170
!end if

!Using grid distance cutoff may hinder localization convergence for e.g. Li6 cluster (deviation of normalization to 1 of this case is quite large)
!So, force to disable this trick. This is probably the contribution of basis function overlap far from atomic center is nonnegligible for sparse system
radcut_old=radcut
radcut=0
call gen1cintgrid(gridatmorg,iradcut) !Generate integration grid

write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call walltime(nwalltime1)

AOMbas=0
ifinish=0
call showprog(0,ncenter)
!Cycle each atom
do iatm=1,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
    gridatm%value=gridatmorg%value
	call gen1cbeckewei(iatm,iradcut,gridatm,atmspcweight,covr_tianlu,3)
    
    !$OMP parallel shared(AOMbas) private(ipt,ibas,jbas,AOMtmp,basval,weitmp,weitmp2) num_threads(nthreads)
    AOMtmp=0D0
    !$OMP do schedule(dynamic)
    do ipt=1+iradcut*sphpot,radpot*sphpot
        weitmp=atmspcweight(ipt)*gridatm(ipt)%value
	    !call orbderv(1,1,nbasis,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,basval)
        call calcbasval(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,basval) !Faster than using above line
        !Original version, slower
     !   do ibas=1,nbasis
     !       weitmp2=basval(ibas)*weitmp
		   ! do jbas=ibas,nbasis
			  !  AOMtmp(ibas,jbas)=AOMtmp(ibas,jbas)+basval(jbas)*weitmp2
		   ! end do
	    !end do
		do jbas=1,nbasis
            AOMtmp(jbas:nbasis,jbas)=AOMtmp(jbas:nbasis,jbas)+basval(jbas:nbasis)*basval(jbas)*weitmp
	    end do
    end do
    !$OMP end do
    !$OMP CRITICAL
	    AOMbas(:,:,iatm)=AOMbas(:,:,iatm)+AOMtmp(:,:)
    !$OMP end CRITICAL
    !$OMP end parallel
    
    ifinish=ifinish+1
    call showprog(ifinish,ncenter)
end do !End cycling atoms
    
!Original version, slower
!do ibas=1,nbasis
!    do jbas=ibas+1,nbasis
!        AOMbas(jbas,ibas,:)=AOMbas(ibas,jbas,:)
!    end do
!end do
do jbas=1,nbasis
    do ibas=jbas,nbasis
		AOMbas(jbas,ibas,:)=AOMbas(ibas,jbas,:)
	end do
end do

call walltime(nwalltime2)
write(*,"(' Generation of atomic overlap matrix took up',i8,' seconds wall clock time')") nwalltime2-nwalltime1

!Check quality of AOMbas
devmax=0
do ibas=1,nbasis
    tmp=sum(AOMbas(ibas,ibas,:))
    !write(*,"(' Basis function:',i6,'    Normalization:',f16.10)") ibas,tmp
    if (abs(1-tmp)>devmax) devmax=abs(1-tmp)
end do
write(*,"(a,f12.8)") " Maximal deviation of normalization of basis function to unity:",devmax

deallocate(COtr)
CObasa=CObasa_org
CO=CO_org
radcut=radcut_old
end subroutine



!!---- Return electric dipole moment vector, wavefunction information must be available in either b or CObasa
subroutine get_dipole_moment(vec)
use defvar
implicit real*8 (a-h,o-z)
real*8 vec(3)

if (allocated(CObasa)) then
    write(*,*) "Calculating electric dipole moment integral matrix..."
    call genDbas_curr
    xeledip=sum(Dbas(1,:,:)*Ptot(:,:))
    yeledip=sum(Dbas(2,:,:)*Ptot(:,:))
    zeledip=sum(Dbas(3,:,:)*Ptot(:,:))
else if (allocated(b)) then
    write(*,*) "Calculating density matrix based on GTFs..."
    call genPprim
    write(*,*) "Calculating electric dipole moment integral matrix..."
    call genDprim
    xeledip=sum(Dprim(1,:,:)*Ptot_prim(:,:))
    yeledip=sum(Dprim(2,:,:)*Ptot_prim(:,:))
    zeledip=sum(Dprim(3,:,:)*Ptot_prim(:,:))
end if
xnucdip=0
ynucdip=0
znucdip=0
do iatm=1,ncenter
    xnucdip=xnucdip+a(iatm)%x*a(iatm)%charge
    ynucdip=ynucdip+a(iatm)%y*a(iatm)%charge
    znucdip=znucdip+a(iatm)%z*a(iatm)%charge
end do
vec(1)=xeledip+xnucdip
vec(2)=yeledip+ynucdip
vec(3)=zeledip+znucdip
end subroutine



!!------ Generate "fragatm" array containing all atoms in the present system
subroutine genfragatm
use defvar

if (allocated(fragatm)) deallocate(fragatm)
allocate(fragatm(ncenter))
nfragatm=ncenter
forall (i=1:nfragatm) fragatm(i)=i
ifragcontri=0 !Fragment has not been manually defined by user
end subroutine


!!-------- Export all internal geometry parameters to int_coord.txt in current folder
subroutine showgeomparam(outname,iallowPBC)
use defvar
use util
character(len=*) outname
integer iallowPBC
integer nneigh(ncenter) !Number of neighbours (the number of atoms connect to this atom)
integer neigh(maxneigh,ncenter) !neigh(1:nneigh(i),i) is list of neighbouring atom indices of atom i
integer nbond !Number of bonds
integer,allocatable :: bond(:,:) !bond(1/2,i) are indices of the two atoms involved in bond i
integer nangle !Number of angles
integer,allocatable :: angle(:,:) !angle(1/2/3,i) are indices of the three atoms involved in angle i
integer ndih !Number of dihedrals
integer,allocatable :: dih(:,:) !dih(1/2/3/4,i) are indices of the four atoms involved in dihedral i

open(10,file=outname,status="replace")

!Determining neighbours
nneigh=0
do iatm=1,ncenter
    do jatm=1,ncenter
        if (jatm==iatm) cycle
        if (atomdist(iatm,jatm,iallowPBC)<( covr(a(iatm)%index)+covr(a(jatm)%index) )*bondcrit) then
            nneigh(iatm)=nneigh(iatm)+1
            neigh(nneigh(iatm),iatm)=jatm
        end if
    end do
end do

!Determining bond list
do itime=1,2
    nbond=0
    do iatm=1,ncenter
        do jatm=iatm+1,ncenter
            if (any(neigh(1:nneigh(iatm),iatm)==jatm)) then
                nbond=nbond+1
                if (itime==2) then
                    bond(1,nbond)=iatm
                    bond(2,nbond)=jatm
                end if
            end if
        end do
    end do
    if (itime==1) allocate(bond(2,nbond))
end do
write(10,"(' Number of bonds:',i8)") nbond
do ibond=1,nbond
    iatm=bond(1,ibond)
    jatm=bond(2,ibond)
    write(10,"(' #',i5,'      Atoms:',2i5,'    Distance:',f12.6,' Angstrom')") ibond,iatm,jatm,atomdist(iatm,jatm,iallowPBC)*b2a
end do

!Determining angle list
do itime=1,2
    nangle=0
    do icen=1,ncenter
        if (nneigh(icen)<2) cycle
        do idx=1,nneigh(icen)
            iatm=neigh(idx,icen)
            do jdx=idx+1,nneigh(icen)
                jatm=neigh(jdx,icen)
                nangle=nangle+1
                if (itime==2) then
                    angle(1,nangle)=iatm
                    angle(2,nangle)=icen
                    angle(3,nangle)=jatm
                end if
            end do
        end do
    end do
    if (itime==1) allocate(angle(3,nangle))
end do
!Sort angle array, making index from small to large, from left to right
call sortidxlist(angle,3,nangle)
write(10,*)
write(10,"(' Number of angles:',i8)") nangle
do iangle=1,nangle
    iatm=angle(1,iangle)
    jatm=angle(2,iangle)
    katm=angle(3,iangle)
    angval=atomang(iatm,jatm,katm,iallowPBC)
    write(10,"(' #',i5,'      Atoms:',3i5,'    Angle:',f12.6,' degree')") iangle,iatm,jatm,katm,angval
end do

!Determining dihedral list. iatm-icen-jcen-jatm
nlinear=0
do itime=1,2
    ndih=0
    do ibond=1,nbond
        icen=bond(1,ibond)
        jcen=bond(2,ibond)
        neiA=nneigh(icen)
        neiB=nneigh(jcen)
        if (neiA==0.or.neiB==0) cycle
        do idx=1,neiA
            iatm=neigh(idx,icen)
            if (iatm==jcen) cycle
            do jdx=1,neiB
                jatm=neigh(jdx,jcen)
                if (jatm==icen) cycle
                if (iatm==jatm) cycle
                angdev1=180-atomang(iatm,icen,jcen,iallowPBC)
                angdev2=180-atomang(icen,jcen,jatm,iallowPBC)
                if (angdev1<0.5D0.or.angdev2<0.5D0) then
                    nlinear=nlinear+1
                    cycle
                end if
                ndih=ndih+1
                if (itime==2) then
                    if (jatm>iatm) then
                        dih(1,ndih)=iatm
                        dih(2,ndih)=icen
                        dih(3,ndih)=jcen
                        dih(4,ndih)=jatm
                    else !Require the index of the fourth atom is larger than that of the first atom
                        dih(1,ndih)=jatm
                        dih(2,ndih)=jcen
                        dih(3,ndih)=icen
                        dih(4,ndih)=iatm
                    end if
                end if
            end do
        end do
    end do
    if (itime==1) allocate(dih(4,ndih))
end do
if (nlinear>0) then
    write(*,"(a,i5,a)") " Note:",nlinear," dihedrals deviate from linear less than 0.5 degree, they are ignored"
end if

!Sort dih array, making index from small to large, from left to right
call sortidxlist(dih,4,ndih)
write(10,*)
write(10,"(' Number of dihedrals:',i8)") ndih
do idih=1,ndih
    iatm=dih(1,idih)
    jatm=dih(2,idih)
    katm=dih(3,idih)
    latm=dih(4,idih)
    dihval=atomdih(iatm,jatm,katm,latm,iallowPBC)
    write(10,"(' #',i5,'      Atoms:',4i5,'    Dihedral:',f12.6,' degree')") idih,iatm,jatm,katm,latm,dihval
end do

close(10)

end subroutine




!!------ Generate Z-matrix. The each line of Zmat corresponds to an entry of the matrix, atomic indices are recorded
!If ierror=1, then the generation is failed (e.g. C2H2). =0 means normal
subroutine genZmat(Zmat,ierror)
use defvar
use util
implicit real*8 (a-h,o-z)
integer Zmat(ncenter,3),ierror
integer nneigh(ncenter) !Number of neighbours (the number of atoms connected to this atom)
integer neigh(maxneigh,ncenter) !neigh(1:nneigh(i),i) is list of neighbouring atom indices of atom i

call genneighlist(1,1,neigh,nneigh) !Determining neighbours according to connmat

ierror=0
Zmat=0 !Element of 0 means unavailable

!Content of a line: iatm i1 i2 i3
!Will determine suitable i1, then i2, then i3
!When determining the suitable one w.r.t. current atom p, first search an atom with index less than iatm and connects to p, &
!if nothing is found, use an atom with index less than iatm and nearest to p
!Three atoms should not lie nearly a straight line, the angle must be >0.5 and <179.5 degree

do iatm=2,ncenter
    do jatm=1,iatm-1 !If an atom has index less than iatm, use it
        if (any(neigh(1:nneigh(iatm),iatm)==jatm)) then
            i1=jatm
            exit
        end if
    end do
    if (jatm==iatm) then !No atom connected to iatm, find the atom nearest to the current one
        distmin=1E10
        do jatm=1,iatm-1
            dist=atomdist(iatm,jatm,0)
            if (dist<distmin) then
                i1=jatm
                distmin=dist
            end if
        end do
    end if
    Zmat(iatm,1)=i1
    
    if (iatm>=3) then
        do jatm=1,iatm-1 !If an atom has index less than i1, use it
            if (jatm==i1) cycle
            if (any(neigh(1:nneigh(i1),i1)==jatm)) then
				if (ncenter>3) then
					tmpang=atomang(iatm,i1,jatm,0)
					if (tmpang<0.5D0.or.tmpang>179.5D0) cycle
                end if
                i2=jatm
                exit
            end if
        end do
        if (jatm==iatm) then !No atom connected to iatm, find the atom nearest to the current one
            distmin=1E10
            do jatm=1,iatm-1
                if (jatm==i1) cycle
                dist=atomdist(i1,jatm,0)
                if (dist<distmin) then
					if (ncenter>3) then
						tmpang=atomang(iatm,i1,jatm,0)
						if (tmpang<0.5D0.or.tmpang>179.5D0) cycle
                    end if
                    i2=jatm
                    distmin=dist
                end if
            end do
            if (distmin==1E10) then !distmin was unchanged, no atom satisfies condition
				ierror=1
                return
            end if
        end if
        Zmat(iatm,2)=i2
    
        if (iatm>=4) then
            do jatm=1,iatm-1 !If an atom has index less than iatm, use it
                if (jatm==i1.or.jatm==i2) cycle
                if (any(neigh(1:nneigh(i2),i2)==jatm)) then
					tmpang=atomang(i1,i2,jatm,0)
					if (tmpang<0.5D0.or.tmpang>179.5D0) cycle
					i3=jatm
					exit
                end if
            end do
            if (jatm==iatm) then !No atom connected to iatm, find the atom nearest to the current one
                distmin=1E10
                do jatm=1,iatm-1
                    if (jatm==i1.or.jatm==i2) cycle
                    dist=atomdist(i2,jatm,0)
                    if (dist<distmin) then
						tmpang=atomang(i1,i2,jatm,0)
						if (tmpang<0.5D0.or.tmpang>179.5D0) cycle
                        i3=jatm
                        distmin=dist
                    end if
                end do
				if (distmin==1E10) then
					ierror=1
					return
				end if
            end if
            Zmat(iatm,3)=i3
        end if
    end if
end do
end subroutine




!!-------- Load OpenMP stacksize from a string, which may correspond to that set by OMP_STACKSIZE or KMP_STACKSIZE
!The string should not be empty. See https://www.openmp.org/spec-html/5.0/openmpse54.html for more information about format
subroutine read_ompstacksize(c200tmp)
use defvar
use util
character(len=*) c200tmp

iunit=2 !KB, this is default unit according to OpenMP standard when unit is not explicitly specified
if (index(c200tmp,'b')/=0.or.index(c200tmp,'B')/=0) iunit=1 !Bytes
if (index(c200tmp,'M')/=0.or.index(c200tmp,'n')/=0) iunit=3 !MBytes
if (index(c200tmp,'G')/=0.or.index(c200tmp,'g')/=0) iunit=4 !GBytes
call remove_char(c200tmp,'b')
call remove_char(c200tmp,'B')
call remove_char(c200tmp,'k')
call remove_char(c200tmp,'K')
call remove_char(c200tmp,'m')
call remove_char(c200tmp,'M')
call remove_char(c200tmp,'g')
call remove_char(c200tmp,'G')
read(c200tmp,*) ompstacksize
!ompstacksize is recorded in Bytes
if (iunit==2) ompstacksize=ompstacksize*1024
if (iunit==3) ompstacksize=ompstacksize*1024**2
if (iunit==4) ompstacksize=ompstacksize*1024**3
end subroutine




!!---- Print information to let user to input orbital label
!Only supports single-determinant wavefunction, except for RO, and basis function should be available
subroutine orblabsel_prompt
use defvar
if ((wfntype==0.or.wfntype==1).and.allocated(CObasa)) then
	if (wfntype==0) then
		write(*,*) "You can also input for example: (quotation marks should not be inputted)"
		write(*,*) """h"" for HOMO, ""h-5"" for HOMO-5, ""l"" for LUMO, ""l+3"" for LUMO+3"
	else
		write(*,*) "You can also input for example: (quotation marks should not be inputted)"
		write(*,*) """ha"" for HOMO(alpha), ""hb"" for HOMO(beta), ""hb-5"" for HOMO(beta)-5"
		write(*,*) """la"" for LUMO(alpha), ""lb"" for LUMO(beta), ""la+3"" for LUMO(alpha)+3"
	end if
end if
end subroutine

!!------- Return orbital index according to inputted orbital label
!Only supports single-determinant wavefunction, except for RO, and basis function should be available
!iorb is orbital index. If failed to determine, return 0
!Suggested usage:
!if (index(c200tmp,'h')==0.and.index(c200tmp,'l')==0) then 
!   ... directly load orbital index
!else
!   call orblabsel(c200tmp,iorb)
!	if (ishowmo==0) then
!		write(*,*) "Error: The orbital label you inputted is wrong! Please double check"
!		...cycle
!	end if
!end if
subroutine orblabsel(orblab,iorb)
use defvar
character(len=*) orblab
integer iorb,tmpval
if (index(orblab,',')/=0) then
    write(*,*) "Error: In this selection mode, you can select only one orbital each time"
    iorb=0
    return
end if
call getHOMOidx
tmpval=0
orblab=adjustl(orblab)
if (orblab(1:2)=="ha") then
    if (orblab(3:3)=='-') read(orblab(3:),*) tmpval
    iorb=idxHOMO+tmpval
else if (orblab(1:2)=="hb") then
    if (orblab(3:3)=='-') read(orblab(3:),*) tmpval
    iorb=idxHOMOb+tmpval
else if (orblab(1:2)=="h ".or.orblab(1:2)=="h-") then
    if (orblab(2:2)=='-') read(orblab(2:),*) tmpval
    iorb=idxHOMO+tmpval
else if (orblab(1:2)=="la") then
    if (orblab(3:3)=='+') read(orblab(3:),*) tmpval
    iorb=idxHOMO+1+tmpval
else if (orblab(1:2)=="lb") then
    if (orblab(3:3)=='+') read(orblab(3:),*) tmpval
    iorb=idxHOMOb+1+tmpval
else if (orblab(1:2)=="l ".or.orblab(1:2)=="l+") then
    if (orblab(2:2)=='+') read(orblab(2:),*) tmpval
    iorb=idxHOMO+1+tmpval
else
    iorb=0
end if
if (iorb<1.or.iorb>nmo) iorb=0
end subroutine



!!----------- Transform MO energies to Fock matrix. PBC is supported
!The numerical accuracy may be marginally lower than that directly loaded from e.g. 47
!istatus is returned variable. It is 0 by default, when one or more orbital energies are zero and thus cannot generate Fock, istatus=1
subroutine MOene2Fmat(istatus)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 tmpmat(nbasis,nbasis),Emat(nbasis,nbasis),Cinv(nbasis,nbasis)
integer istatus

istatus=0
if (.not.allocated(CObasa)) then
	write(*,*) "Error: Basis function information is not available!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    istatus=1
    return
end if
nzero=count(MOene==0)
if (nzero>0) then
    if (nzero==nbasis) then
        write(*,"(a)") " Error: All orbital energies are zero, hence Fock/KS matrix &
        &cannot be generated based on orbital energies and expansion coefficients. Commonly this is because the orbital energies were not &
        &solved at all during your quantum chemistry calculation. For example, the current orbitals are not molecular orbitals, or you enabled OT in CP2K"
    else
        write(*,"(a,i6,a)") " Error:",nzero," orbital energies are zero, hence Fock/KS matrix &
        &cannot be generated based on orbital energies and expansion coefficients. Commonly this is because there is linear dependency problem &
        &in your quantum chemistry calculation, and thus some basis functions are automatically removed by your quantum chemistry code. &
        &To tackle this problem, removing diffuse functions if they have been used; if you are a Gaussian user, also add IOp(3/32=2) to fully avoid &
        &eliminating linear-dependent basis functions"
    end if
    write(*,*) "Press ENTER button to return"
    read(*,*)
    istatus=1
    return
end if

call ask_Sbas_PBC

write(*,"(a)") " Generating Fock/KS matrix based on orbital energies, expansion coefficients and overlap matrix..."

!Because FC=SCE, thus F=SCE(C)^-1
if (allocated(FmatA)) deallocate(FmatA)
allocate(FmatA(nbasis,nbasis))
Cinv=invmat(CObasa,nbasis)
Emat=0
forall(i=1:nbasis) Emat(i,i)=MOene(i)
!FmatA=matmul(matmul(matmul(Sbas,CObasa),Emat),Cinv) !Slower code
tmpmat=matmul_blas(Sbas,CObasa,nbasis,nbasis,0,0) !Parallel MKL
tmpmat=matmul_blas(tmpmat,Emat,nbasis,nbasis,0,0)
FmatA=matmul_blas(tmpmat,Cinv,nbasis,nbasis,0,0)

if (wfntype==1.or.wfntype==4) then !Beta part
    if (allocated(FmatB)) deallocate(FmatB)
    allocate(FmatB(nbasis,nbasis))
    Cinv=invmat(CObasb,nbasis)
    Emat=0
    forall(i=1:nbasis) Emat(i,i)=MOene(i+nbasis)
    tmpmat=matmul_blas(Sbas,CObasb,nbasis,nbasis,0,0) !Parallel MKL
    tmpmat=matmul_blas(tmpmat,Emat,nbasis,nbasis,0,0)
    FmatB=matmul_blas(tmpmat,Cinv,nbasis,nbasis,0,0)
end if

write(*,*) "Fock/KS matrix has been successfully generated!"

!Compare with that loaded from .47
!tmpmat=FmatA
!call loadFockfile(istatus)
!!call showmatgau(FmatA,form="f14.8")
!write(*,"(' Maximum deviation to loaded Fock matrix',f20.10)") maxval(abs(tmpmat-FmatA))
end subroutine



!!----------- Find element index from element name. If idx returns 0, means the index was not found
subroutine elename2idx(name,idx)
use defvar
use util
character(len=*) name
integer idx

idx=0
call lc2uc(name(1:1)) !Convert to upper case
call uc2lc(name(2:2)) !Convert to lower case
do iele=1,nelesupp
	if ( name(1:2)==ind2name(iele) ) then
		idx=iele
		return
	end if
end do
end subroutine



!!--------- Show menu title with center alignment, e.g.          ------ ltwd ------
!str is the string of title, nsymbol_in is number - to be shown at each side. If it is 0, then print - as much as possible to fill all blank spaces
!itype=1: Show -   itype=2: Show =
!Example: call menutitle("Post-processing menu",10,1)
subroutine menutitle(str,nsymbol_in,itype)
implicit real*8 (a-h,o-z)
integer nsymbol,strlen,itype
character(len=*) str
character c80tmp*80

strlen=len_trim(str)
maxnsym=(80-(strlen+4))/2
if (nsymbol_in>maxnsym) then
	nsymbol=maxnsym
else
	nsymbol=nsymbol_in
end if
nchar=nsymbol*2+2+strlen !Number of characters to print without spaces at two sides
nspace=int((80-nchar)/2)
c80tmp=" "
do i=nspace+1,nspace+nsymbol
	if (itype==1) then
		c80tmp(i:i)='-'
    else
		c80tmp(i:i)='='
    end if
end do
c80tmp(nspace+nsymbol+2:nspace+nsymbol+1+strlen)=str
itmp=nspace+nsymbol+1+strlen
do i=itmp+2,itmp+1+nsymbol
	if (itype==1) then
		c80tmp(i:i)='-'
    else
		c80tmp(i:i)='='
    end if
end do
write(*,"(a)") trim(c80tmp)
end subroutine



!!--------- Find the atom (iatmclose) closest to a given point and return the distance (r)
!PBC is considered
subroutine closest_atm_pt(x,y,z,iatmclose,r)
use defvar
implicit real*8 (a-h,o-z)
real*8 x,y,z,r,tvec(3)
integer iatmclose

r2close=1D100
if (ifPBC>0) then
	call getpointcell(x,y,z,ic,jc,kc)
	do icell=ic-PBCnx,ic+PBCnx
		do jcell=jc-PBCny,jc+PBCny
			do kcell=kc-PBCnz,kc+PBCnz
				call tvec_PBC(icell,jcell,kcell,tvec)
				do jatm=1,ncenter
					atmx=a(jatm)%x+tvec(1)
					atmy=a(jatm)%y+tvec(2)
					atmz=a(jatm)%z+tvec(3)
					r2=(atmx-x)**2 + (atmy-y)**2 + (atmz-z)**2
                    if (r2<r2close) then
						r2close=r2
                        iatmclose=jatm
                    end if
				end do
			end do
		end do
	end do
else
	do jatm=1,ncenter
		r2=(a(jatm)%x-x)**2+(a(jatm)%y-y)**2+(a(jatm)%z-z)**2
        if (r2<r2close) then
			r2close=r2
            iatmclose=jatm
        end if
	end do
end if
r=dsqrt(r2close)
end subroutine



!!---- Setup move vectors and move lengths towards 26 neighbouring grid. Used in basin generation and domain analysis
!vec26x/y/z(i) is change of grid index in directions 1/2/3 of movement mode i
!len26(i) is move distance (Bohr) of movement mode i, will be used to determine gradient of corresponding movement
subroutine setupmovevec
use defvar
use basinintmod
real*8 vec(3)

vec26x=0
vec26y=0
vec26z=0
!The nearest neighbours:
vec26x(1)=1
vec26x(2)=-1
vec26y(3)=1
vec26y(4)=-1
vec26z(5)=1
vec26z(6)=-1
!On the edges:
vec26x(7)=1
vec26y(7)=1
vec26x(8)=-1
vec26y(8)=1
vec26x(9)=-1
vec26y(9)=-1
vec26x(10)=1
vec26y(10)=-1
vec26x(11)=1
vec26z(11)=1
vec26x(12)=-1
vec26z(12)=1
vec26x(13)=-1
vec26z(13)=-1
vec26x(14)=1
vec26z(14)=-1
vec26y(15)=1
vec26z(15)=1
vec26y(16)=1
vec26z(16)=-1
vec26y(17)=-1
vec26z(17)=-1
vec26y(18)=-1
vec26z(18)=1
!At the vertices:
vec26z(19:22)=1
vec26x(19)=1
vec26y(19)=1
vec26x(20)=-1
vec26y(20)=1
vec26x(21)=-1
vec26y(21)=-1
vec26x(22)=1
vec26y(22)=-1
vec26z(23:26)=-1
vec26x(23)=1
vec26y(23)=-1
vec26x(24)=1
vec26y(24)=1
vec26x(25)=-1
vec26y(25)=1
vec26x(26)=-1
vec26y(26)=-1
do i=1,26
	vec=vec26x(i)*gridv1+vec26y(i)*gridv2+vec26z(i)*gridv3
    len26(i)=dsqrt(sum(vec**2))
	!write(*,"(i3,3i5,f12.6)") i,vec26x(i),vec26y(i),vec26z(i),len26(i)
end do
end subroutine



!!------ Convert bascen to basstart and basend
!basstart and basend should have been allocated
subroutine bascen2basstart_end
use defvar
implicit real*8 (a-h,o-z)

nowcen=0
indcen=0
basstart=0 !Bq atom without basis function will have 0 for basstart and basend
basend=0
do ibasis=1,nbasis
	if (bascen(ibasis)/=nowcen) then
		nowcen=bascen(ibasis)
		indcen=indcen+1
		basstart(indcen)=ibasis
		if (indcen/=1) basend(indcen-1)=ibasis-1
	end if
end do
basend(indcen)=nbasis
!do iatm=1,ncenter
!	write(*,*) iatm,basstart(iatm),basend(iatm)
!end do
end subroutine




!!------- Show minimum and maximum of a curve
!npointcurve: Number of points of xdata and ydata
!ixunit=1: xdata is in Bohr, =2: xdata is in Angstrom
subroutine showcurveminmax(npointcurve,xdata,ydata,ixunit)
use defvar
implicit real*8 (a-h,o-z)
integer npointcurve
real*8 xdata(npointcurve),ydata(npointcurve)

numlocmin=0
numlocmax=0
do ipt=2,npointcurve-1
	gradold=ydata(ipt)-ydata(ipt-1)
	gradnew=ydata(ipt+1)-ydata(ipt)
	if (gradold*gradnew<0D0) then
		if (gradold>gradnew) then
			numlocmax=numlocmax+1
			if (ixunit==1) then
				write(*,"(' Maximum X (Bohr):',f12.6,'  Value:',E18.8)") xdata(ipt),ydata(ipt)
			else
				write(*,"(' Maximum X (Angstrom):',f12.6,'  Value:',E18.8)") xdata(ipt)*b2a,ydata(ipt)
			end if
		else if (gradold<gradnew) then
			numlocmin=numlocmin+1
			if (ixunit==1) then
				write(*,"(' Minimum X (Bohr):',f12.6,'  Value:',E18.8)") xdata(ipt),ydata(ipt)
			else
				write(*,"(' Minimum X (Angstrom):',f12.6,'  Value:',E18.8)") xdata(ipt)*b2a,ydata(ipt)
			end if
		end if
	end if
end do
write(*,"(' Totally found',i5,' minima,',i5,' maxima')") numlocmin,numlocmax
end subroutine




!!------- Interface of "occorb_Lowdinorth" for general purpose, CObasa/CObasb and Ptot/Palpha/Pbeta will be updated
subroutine occorb_Lowdinorth_wrapper
use defvar
implicit real*8 (a-h,o-z)
real*8 occfragB(nbasis),CObasapro(nbasis,nbasis)
real*8,allocatable :: CObasbpro(:,:)
real*8 PfrzA(0,0),PfrzB(0,0),CObasbfrz(0,0)

CObasapro=CObasa
if (wfntype==0) then
    allocate(CObasbpro(0,0))
    call occorb_Lowdinorth(0,0,0,MOocc,occfragB,CObasapro,CObasbpro,CObasa,CObasbfrz,Ptot,PfrzA,PfrzB)
    call CObas2CO(1)
else if (wfntype==1) then
    allocate(CObasbpro(nbasis,nbasis))
    CObasbpro=CObasb
    call occorb_Lowdinorth(0,1,nbasis,MOocc(1:nbasis),MOocc(nbasis+1:nmo),CObasapro,CObasbpro,CObasa,CObasb,Ptot,Palpha,Pbeta)
    call CObas2CO(3)
else
    write(*,"(a)") " Error: Lowdin orthogonalization between occupied orbitals can only be performed for single-determinant wavefunction"
end if
write(*,*) "Done! Orbital coefficients and density matrix have been updated"
end subroutine


!!------- Perform Lowdin orthogonalization between occupied orbitals. Density matrix is also updated. Mainly used in ETS-NOCV
!Sbas must has been available. Single-determinant wavefunction is assumed
!infomode=0: Silent   =1: Show information
!iopsh=0: closed-shell   =1: unrestricted open-shell
!occfrag,occfragB: List of occupation number of fragment orbitals, for spatial/alpha and beta spins
!CObasapro,CObasbpro*: Coefficient matrix of original orbitals
!CObasafrz,CObasbfrz*: Coefficient matrix of new orbitals
!Pfrz,PfrzA*,PfrzB*: Total, alpha, beta density matrix constructed based on new orbitals
!NOTE: For closed-shell case, ndimB should be 0, and in this case all matrics with asterisk above will have size of 0 and thus do not waste memory
subroutine occorb_Lowdinorth(infomode,iopsh,ndimB,occfrag,occfragB,CObasapro,CObasbpro,CObasafrz,CObasbfrz,Pfrz,PfrzA,PfrzB)
use defvar
use util
implicit real*8 (a-h,o-z)
integer infomode,iopsh
real*8 occfrag(nbasis),occfragB(nbasis)
integer,allocatable :: occidx(:),occidxB(:) !occidx(i)=j : ith occupied orbital corresponds to jth global orbital index. For all/alpha and beta spins
real*8 CObasapro(nbasis,nbasis),CObasbpro(ndimB,ndimB)
real*8 CObasafrz(nbasis,nbasis),CObasbfrz(ndimB,ndimB)
real*8 Pfrz(nbasis,nbasis),PfrzA(ndimB,ndimB),PfrzB(ndimB,ndimB)
real*8,allocatable :: Umat(:,:),svalvec(:),Xmat(:,:),tmpmat(:,:)
real*8 FOovlp(nbasis,nbasis) !Overlap matrix between fragment orbitals (spatial or alpha)
real*8,allocatable :: FOovlpB(:,:) !Same as above, for beta spin. Used only for open-shell case
real*8,allocatable :: FOovlpocc(:,:) !Overlap matrix between occupied fragment orbitals (spatial or alpha)
real*8,allocatable :: FOovlpoccB(:,:) !Same as above, for beta spin. Used only for open-shell case

!Generate correspondence of index between occupied orbitals and global orbitals
nocc=nint(naelec)
allocate(FOovlpocc(nocc,nocc),occidx(nocc))
iocc=0
do imo=1,nbasis
    if (occfrag(imo)==0) cycle
    iocc=iocc+1
    occidx(iocc)=imo
end do
if (iopsh==1) then !Beta part
    noccB=nint(nbelec)
    allocate(FOovlpoccB(noccB,noccB),occidxB(noccB))
    iocc=0
    do imo=1,nbasis
        if (occfragB(imo)==0) cycle
        iocc=iocc+1
        occidxB(iocc)=imo
    end do
end if

!Transform overlap integral matrix from AO basis to FO basis
if (infomode==1) then
    write(*,*)
    write(*,*) "Calculating overlap integral matrix between fragment orbitals ..."
end if
allocate(tmpmat(nbasis,nbasis))
tmpmat=matmul_blas(Sbas,CObasapro,nbasis,nbasis)
FOovlp=matmul_blas(CObasapro,tmpmat,nbasis,nbasis,1,0)
!Generate overlap integral matrix between occupied combined fragment orbitals
do imo=1,nocc
    do jmo=1,nocc
        FOovlpocc(imo,jmo)=FOovlp(occidx(imo),occidx(jmo))
    end do
end do
if (iopsh==1) then !Beta part
    tmpmat=matmul_blas(Sbas,CObasbpro,nbasis,nbasis)
    FOovlpB=matmul_blas(CObasbpro,tmpmat,nbasis,nbasis,1,0)
    do imo=1,noccB
        do jmo=1,noccB
            FOovlpoccB(imo,jmo)=FOovlpB(occidxB(imo),occidxB(jmo))
        end do
    end do
end if
deallocate(tmpmat)

!Lowdin orthogonalization between all occupied orbitals. See Eq. 7 of ETS-NOCV original paper
!Unoccupied fragment orbitals keep unchanged since they are not involved in construction of frozen state density matrix, &
! Do not follow ETS-NOCV paper's statement, it is redundant and unnecessary (they also use Lowdin orthogonalization between fragment unoccupied orbitals, &
! and then do Schmidt orthogonalization for unoccupied orbitals w.r.t. occupied ones)
if (infomode==1) then
    write(*,*)
    write(*,*) "Lowdin orthogonalization between occupied fragment orbitals ..."
end if
allocate(Umat(nocc,nocc),svalvec(nocc),Xmat(nocc,nocc),tmpmat(nocc,nocc))
call diagsymat(FOovlpocc,Umat,svalvec,ierror)
if (ierror/=0) write(*,*) "Error: Diagonalization of overlap matrix failed!"
tmpmat=0
forall (i=1:nocc) tmpmat(i,i)=dsqrt(1/svalvec(i)) !Use Sbas as temporary matrix here
Xmat=matmul(matmul(Umat,tmpmat),transpose(Umat)) !Then Xmat is S^(-1/2)
CObasafrz=CObasapro
do imo=1,nocc
    CObasafrz(:,occidx(imo))=0
    do jmo=1,nocc
        CObasafrz(:,occidx(imo))=CObasafrz(:,occidx(imo))+CObasapro(:,occidx(jmo))*Xmat(imo,jmo)
    end do
end do
deallocate(Umat,svalvec,Xmat,tmpmat)
if (iopsh==1) then !Beta part
    allocate(Umat(noccB,noccB),svalvec(noccB),Xmat(noccB,noccB),tmpmat(noccB,noccB))
    call diagsymat(FOovlpoccB,Umat,svalvec,ierror)
    if (ierror/=0) write(*,*) "Error: Diagonalization of overlap matrix failed!"
    tmpmat=0
    forall (i=1:noccB) tmpmat(i,i)=dsqrt(1/svalvec(i)) !Use Sbas as temporary matrix here
    Xmat=matmul(matmul(Umat,tmpmat),transpose(Umat)) !Then Xmat is S^-0.5
    CObasbfrz=CObasbpro
    do imo=1,noccB
        CObasbfrz(:,occidxB(imo))=0
        do jmo=1,noccB
            CObasbfrz(:,occidxB(imo))=CObasbfrz(:,occidxB(imo))+CObasbpro(:,occidxB(jmo))*Xmat(imo,jmo)
        end do
    end do
    deallocate(Umat,svalvec,Xmat,tmpmat)
end if

!Generate frozen state density matrix
if (infomode==1) then
    write(*,*)
    write(*,*) "Generating frozen state density matrix (P_frz) ..."
end if
if (iopsh==0) then !Closed-shell case
    Pfrz=0
    do iocc=1,nocc
        imo=occidx(iocc)
	    Pfrz=Pfrz+occfrag(imo)*matmul(CObasafrz(:,imo:imo),transpose(CObasafrz(:,imo:imo)))
    end do
    if (infomode==1) write(*,"(' Tr(P_frz*S):',f12.6)") sum(Pfrz*Sbas)
    !call showmatgau(Pfrz,label="Pfrz",form="f12.8")
else if (iopsh==1) then !Open shell case, generate alpha and beta respectively, and also combine to total
    PfrzA=0
    do iocc=1,nocc
        imo=occidx(iocc)
	    PfrzA=PfrzA+occfrag(imo)*matmul(CObasafrz(:,imo:imo),transpose(CObasafrz(:,imo:imo)))
    end do
    PfrzB=0
    do iocc=1,noccB
        imo=occidxB(iocc)
	    PfrzB=PfrzB+occfragB(imo)*matmul(CObasbfrz(:,imo:imo),transpose(CObasbfrz(:,imo:imo)))
    end do
    if (infomode==1) then
        write(*,"(' Tr(P_frz_A*S):',f12.6)") sum(PfrzA*Sbas)
        write(*,"(' Tr(P_frz_B*S):',f12.6)") sum(PfrzB*Sbas)
    end if
end if
end subroutine