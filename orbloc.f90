!!------------------------------------------------------
!! ----------- Orbital localization analysis ----------- 
!!------------------------------------------------------
!The final wavefunction can be exported as .fch, I don't select .molden because it doesn't record atomic charge, this will be problematic when ECP is used
subroutine orbloc
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
character c2000tmp*2000
integer :: maxcyc=80,ireload=1,idoene=0,idocore=1,imethod=1,domark(4),iPMexp=2,ilmocen=0
real*8 :: arrayi(nbasis),arrayj(nbasis),crit=1D-4,tmpbasarr(nbasis),tmpprimarr(nprims),bastot(nbasis)
real*8,pointer :: Cmat(:,:)
real*8,allocatable :: FLMOA(:,:),FLMOB(:,:),Xmat(:,:),Xmatinv(:,:),SC(:,:),AOMbas(:,:,:),tmpmat(:,:)
real*8 :: orbcomp(ncenter,nbasis) !Used in printing major orbital character
integer :: orbtype(nmo) !LMO type determined according to composition, 0=other, 1=one-center, 2=two-center
integer :: orbatm(2,nmo) !The index of the atom mainly involved in the LMO. one-center LMO has first element, two-center LMO has two elements
integer :: istatarr(nbasis),iatmarr(ncenter,nbasis) !Used in printing major orbital character
real*8 :: irj(3),iri(3),jrj(3) !Used for Boys localization, store dipole moment integral between two orbitals
real*8 :: crit1c=0.85D0,crit2c=0.80D0
integer :: thres_1c=90,thres_2c=85,icompmethod=2
character c200tmp*200,typestr*4,selectyn
integer orblist(nmo),norblist !Record the index of MOs included in the localization
integer orblist_did(nmo) !If the orbital is involved in orbital localization, it will be set to 1, otherwise 0
!Used by PM-Becke
real*8 tmprow(1,nbasis),mat11(1,1),veccoeffi(1,nbasis),veccoeffj(1,nbasis)
!Below for calculating LMO centers and dipole moments
real*8,allocatable :: orbvalpt(:,:)
real*8 LMOpos(3,nmo),tmpvec(3)
real*8 beckeweigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

if (wfntype==2.or.wfntype==3.or.wfntype==4) then
	write(*,*) "Error: This function only works for single-determinant SCF wavefunction!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
if (.not.allocated(CObasa)) then
	write(*,*) "Error: Basis function information was not provided by your input file!"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if

do while(.true.)
	write(*,*)
	write(*,*) "                ======== Orbital localization analysis ========"
	if (icompmethod>0) write(*,"(a,2i4,' %')") " -10 Set threshold for determining one- and two-center LMO, current:",nint(crit1c*100),nint(crit2c*100)
    if (icompmethod==0) write(*,*) "-9 Set the method for calculating orbital composition, current: Do not calculate"
    if (icompmethod==1) write(*,*) "-9 Set the method for calculating orbital composition, current: Mulliken+SCPA"
    if (icompmethod==2) write(*,*) "-9 Set the method for calculating orbital composition, current: Hirshfeld"
    if (icompmethod==3) write(*,*) "-9 Set the method for calculating orbital composition, current: Becke"
	if (ilmocen==1) write(*,*) "-8 If calculating center position and dipole moment of LMOs, current: Yes"
	if (ilmocen==0) write(*,*) "-8 If calculating center position and dipole moment of LMOs, current: No"
	!By default exponent of 2 is used for PM. Exponent of 4 can also be chosen by option -7.
    !According to my test, p=4 converges slower than p=2 for PM-Mulliken and PM-Lowdin,
    !and the degree of localization is not as substantial as p=2, so I finally decide not to expose this option to users
! 	if (imethod==1.or.imethod==2.or.imethod==3) write(*,"(a,i3)") " -7 Set exponent of Pipek-Mezey method, current:",iPMexp
	if (imethod==1) write(*,*) "-6 Set localization method, current: Pipek-Mezey with Mulliken population"
	if (imethod==2) write(*,*) "-6 Set localization method, current: Pipek-Mezey with Lowdin population"
	if (imethod==3) write(*,*) "-6 Set localization method, current: Pipek-Mezey with Becke population"
	if (imethod==10) write(*,*) "-6 Set localization method, current: Foster-Boys"
	if (idocore==1) write(*,*) "-5 If also localizing core orbitals, current: Yes"
	if (idocore==0) write(*,*) "-5 If also localizing core orbitals, current: No"
	if (idoene==1) write(*,*) "-4 If calculating and printing orbital energies, current: Yes"
	if (idoene==0) write(*,*) "-4 If calculating and printing orbital energies, current: No"
	if (ireload==1) write(*,*) "-3 If reloading newly generated .fch file, current: Yes"
	if (ireload==0) write(*,*) "-3 If reloading newly generated .fch file, current: No"
	write(*,"(a,f12.8)") " -2 Set criterion of convergence, current:",crit
	write(*,"(a,i4)") " -1 Set maximum cycles, current:",maxcyc
	write(*,*) "0 Return"
	write(*,*) "1 Localizing occupied orbitals only"
	write(*,*) "2 Localizing both occupied and unoccupied orbitals separately"
	write(*,*) "3 Localizing specific set of orbitals"
	read(*,*) isel
    
	if (isel==0) then
		return
	else if (isel==-1) then
		write(*,*) "Input maximum cycles, e.g. 30"
		read(*,*) maxcyc
	else if (isel==-2) then
		write(*,*) "Input criterion of convergence, e.g. 0.0001"
		read(*,*) crit
	else if (isel==-3) then
		if (ireload==1) then
			ireload=0
		else if (ireload==0) then
			ireload=1
		end if
	else if (isel==-4) then
		if (idoene==1) then
			idoene=0
		else if (idoene==0) then
			do while(.true.)
				write(*,*)
				write(*,*) "How to evaluate energy of localized orbitals?"
				write(*,*) "0 Do not evaluate"
				write(*,*) "1 Evaluate, using the Fock matrix generated by MO energies and coefficients"
				write(*,*) "2 Evaluate, loading Fock matrix from a file"
				read(*,*) ievalene
				if (ievalene==0) then
					exit
				else if (ievalene==1) then
					call MOene2Fmat(istatus)
				else if (ievalene==2) then
					call loadFockfile(istatus)
				end if
				if (istatus==0) then
					idoene=1
					exit
                end if
			end do
		end if
	else if (isel==-5) then
		if (idocore==1) then
			idocore=0
		else if (idocore==0) then
			idocore=1
		end if
	else if (isel==-6) then
		write(*,*) "Please select orbital localization method"
		if (ifPBC==0) write(*,"(a,/)") " Hint: 1 and 2 are very fast, but may not well work when diffuse functions are presented. 10 &
        &is evidently slower, and 3 is quite time-consuming for large systems, but they are but fully compatible &
        &with diffuse functions. Only 10 is unable to give sigma-pi separated LMOs"
		write(*,*) "1 Pipek-Mezey based on Mulliken type of population"
		write(*,*) "2 Pipek-Mezey based on Lowdin type of population"
        if (ifPBC==0) then
			write(*,*) "3 Pipek-Mezey based on Becke population"
			write(*,*) "10 Foster-Boys"
        end if
		read(*,*) imethod
	else if (isel==-7) then
		write(*,*) "Input exponent of Pipek-Mezey method. 2 or 4 is allowed"
		write(*,"(a)") " Note: Original paper of PM method use exponent of 2, while 4 is shown to give LMO with more localized character"
		read(*,*) iPMexp
		if (iPMexp/=2.and.iPMexp/=4) then
			write(*,*) "Input error! The value should be either 2 or 4!"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			iPMexp=2
		end if
	else if (isel==-8) then
		if (ilmocen==0) then
			if (ireload==0) then
				write(*,"(a)") " Error: To use this function, you must first switch the option ""If reloading newly generated .fch file"" to ""Yes"""
				write(*,*) "Press ENTER button to continue"
				read(*,*)
				cycle
			end if
			ilmocen=1
		else if (ilmocen==1) then
			ilmocen=0
		end if
    else if (isel==-9) then
        !Sometimes SCPA tends to overestimate delocalization character of lone pair type of LMO, &
        !while Mulliken method tends to result in extremely large atomic composition for virtual LMOs
        write(*,"(a,/)") " Hint: Option 1 is very fast, however the method is not robust, and it is even useless when diffuse functions are employed. &
        &Options 2 and 3 usually give similar result, they are more expensive, but fully compatible with diffuse functions"
        write(*,*) "0 Do not calculate orbital compositions"
        write(*,*) "1 Mulliken method for occupied LMOs and SCPA method for unoccupied LMOs"
        write(*,*) "2 Hirshfeld method"
        write(*,*) "3 Becke method"
        read(*,*) icompmethod
	else if (isel==-10) then
		write(*,"(a)") " Input thresholds for identifying one- and two-center LMOs. For example, &
		&inputting 0.9,0.85 means using 90% and 85%, respectively"
		read(*,*) crit1c,crit2c
        
	else if (isel==1.or.isel==2) then
		exit
        
	else if (isel==3) then
        write(*,*) "Select type of orbitals"
        if (wfntype==0) then
            write(*,*) "1 Occupied orbitals"
            write(*,*) "2 Unoccupied orbitals"
        else if (wfntype==1) then
            write(*,*) "1 Occupied alpha orbitals"
            write(*,*) "2 Unoccupied alpha orbitals"
            write(*,*) "3 Occupied beta orbitals"
            write(*,*) "4 Unoccupied beta orbitals"
        end if
        read(*,*) iorbtype
        write(*,"(a)") " Input indices of the orbitals to localize, e.g. 2,3,7-10. The range must belong to the type you selected"
        if (iorbtype==3.or.iorbtype==4) write(*,*) "Note: The beta indices are 1 based"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,norblist,orblist)
        exit
	end if
end do

!Will perform Alpha/total-occ, Alpha/total-vir, Beta-occ, Beta-vir, set a label to indicate which part will actually done
domark=0
if (wfntype==0) then
    if (isel==1) then
	    domark(1)=1
	else if (isel==2) then
        domark(1:2)=1
    else if (isel==3) then
        domark(iorbtype)=1
    end if
else if (wfntype==1) then
    if (isel==1) then
	    domark(1)=1
	    domark(3)=1
	else if (isel==2) then
        domark=1
    else if (isel==3) then
        domark(iorbtype)=1
    end if
end if

if (isel/=3.and.idocore==0) then
	call getninnerele(ninnerele,0) !Count the number of inner electrons
	write(*,"(' Note: Lowest',i5,' orbitals are regarded as core orbitals and will not be localized')") ninnerele/2
end if

!Preparation work for specific method
if (imethod==1) then !PM with Mulliken
    call ask_Sbas_PBC
	allocate(SC(nbasis,nbasis))
else if (imethod==2) then !PM with Lowdin
	write(*,*) "Performing Lowdin orthonormalization..."
    call ask_Sbas_PBC
	allocate(Xmat(nbasis,nbasis),Xmatinv(nbasis,nbasis))
	call symmorthomat(Sbas,Xmat,Xmatinv)
    CObasa=matmul_blas(Xmat,CObasa,nbasis,nbasis) !Transform coefficient matrix to orthogonalized basis
    if (allocated(CObasb)) CObasb=matmul_blas(Xmat,CObasb,nbasis,nbasis)
	Sbas_org=Sbas
	Sbas=0
	forall (i=1:nbasis) Sbas(i,i)=1
else if (imethod==3) then !PM with Becke
    allocate(AOMbas(nbasis,nbasis,ncenter))
    write(*,*) "Generating atomic overlap matrix of basis functions..."
    call genAOMbas(AOMbas)
else if (imethod==10) then !Foster-Boys
	if (.not.allocated(Dbas)) then
		write(*,*) "Generating electric dipole moment integral matrix..."
		call genDbas_curr
	end if
end if

call walltime(iwalltime1)

!Carry out localization
!Alpha-occ,Alpha-vir,Beta-occ,Beta-vir
orblist_did=0
do itime=1,4
	if (domark(itime)==0) cycle
	if (itime<=2) then
		Cmat=>CObasa
	else
		Cmat=>CObasb
	end if
    !Set orbital range to deal with. 1 based index for both alpha and beta
    if (isel==1.or.isel==2) then
	    if (itime==1) then
		    nmobeg=1
		    if (idocore==0) nmobeg=ninnerele/2+1
		    nmoend=naelec
	    else if (itime==2) then
		    nmobeg=naelec+1
		    nmoend=nbasis
	    else if (itime==3) then
		    nmobeg=1
		    if (idocore==0) nmobeg=ninnerele/2+1
		    nmoend=nbelec
	    else if (itime==4) then
		    nmobeg=nbelec+1
		    nmoend=nbasis
	    end if
        norblist=nmoend-nmobeg+1
        forall(i=1:norblist) orblist(i)=nmobeg-1+i
    else if (isel==3) then
        continue !orblist has already been defined by user
    end if
    orblist_did(orblist(1:norblist))=1
	if (imethod==1) then
		write(*,*)
        write(*,*) "Calculating intermediate matrix, please wait..."
        SC=matmul_blas(Sbas,Cmat,nbasis,nbasis,0,0) !Generate intermediate matrix SC for PM-Mulliken for lowering cost from N^4 to N^3
    end if
	
	if (wfntype==0) then
		if (itime==1) write(*,"(/,a)") " Localizing occupied orbitals..."
		if (itime==2) write(*,"(/,a)") " Localizing unoccupied orbitals..."
	else if (wfntype==1) then
		if (itime==1) write(*,"(/,a)") " Localizing alpha occupied orbitals..."
		if (itime==2) write(*,"(/,a)") " Localizing alpha unoccupied orbitals..."
		if (itime==3) write(*,"(/,a)") " Localizing beta occupied orbitals..."
		if (itime==4) write(*,"(/,a)") " Localizing beta unoccupied orbitals..."
	end if

	Pvalold=0
	do icyc=1,maxcyc
		do idx=1,norblist-1 !Cycle each orbital pair
            imo=orblist(idx)
			do jdx=idx+1,norblist
                jmo=orblist(jdx)
                !PM-Mulliken and PM-Lowdin. Working equation is Journal of Computational Chemistry, 14, 6, 736 (1993)
				if (imethod==1.or.imethod==2) then
				    Aval=0;Bval=0
					do iatm=1,ncenter
						is=basstart(iatm)
                        if (is==0) cycle
						ie=basend(iatm)
                        if (imethod==1) then !Mulliken population
						    Qij=0.5D0*sum(Cmat(is:ie,jmo)*SC(is:ie,imo)+Cmat(is:ie,imo)*SC(is:ie,jmo))
						    Qii=sum(Cmat(is:ie,imo)*SC(is:ie,imo))
						    Qjj=sum(Cmat(is:ie,jmo)*SC(is:ie,jmo))
                        else if (imethod==2) then !Lowdin population
						    Qij=sum(Cmat(is:ie,imo)*Cmat(is:ie,jmo))
						    Qii=sum(Cmat(is:ie,imo)*Cmat(is:ie,imo))
						    Qjj=sum(Cmat(is:ie,jmo)*Cmat(is:ie,jmo))
                        end if
						if (iPMexp==2) then
							Aval=Aval+( Qij**2-(Qii-Qjj)**2/4D0 )
							Bval=Bval+( Qij*(Qii-Qjj) )
						else if (iPMexp==4) then
							Aval=Aval-Qii**4-Qjj**4+6*(Qii**2+Qjj**2)*Qij**2+Qii**3*Qjj+Qjj**3*Qii
							Bval=Bval+4*Qij*(Qii**3-Qjj**3)
						end if
					end do
                    !The following parallization makes calculation evidently slower perhaps due to additional overhead due to create threads
     !               !$OMP parallel private(iatm,Avaltmp,Bvaltmp,is,ie,Qij,Qii,Qjj) num_threads(nthreads)
					!Avaltmp=0
					!Bvaltmp=0
     !               !$OMP do schedule(dynamic)
					!do iatm=1,ncenter
					!	is=basstart(iatm)
     !                   if (is==0) cycle
					!	ie=basend(iatm)
					!	Qij=0.5D0*sum(Cmat(is:ie,jmo)*SC(is:ie,imo)+Cmat(is:ie,imo)*SC(is:ie,jmo))
					!	Qii=sum(Cmat(is:ie,imo)*SC(is:ie,imo))
					!	Qjj=sum(Cmat(is:ie,jmo)*SC(is:ie,jmo))
					!	Avaltmp=Avaltmp+( Qij**2-(Qii-Qjj)**2/4D0 )
					!	Bvaltmp=Bvaltmp+( Qij*(Qii-Qjj) )
					!end do
     !               !$OMP end do
     !               !$OMP CRITICAL
	    !            Aval=Aval+Avaltmp
	    !            Bval=Bval+Bvaltmp
     !               !$OMP end CRITICAL
     !               !$OMP end parallel
                !PM-Becke. See Eqs. 14 and 17 in J. Chem. Theory Comput., 10, 642 (2014)
                else if (imethod==3) then
				    Aval=0;Bval=0
                    veccoeffi(1,:)=Cmat(:,imo)
                    veccoeffj(1,:)=Cmat(:,jmo)
                    !$OMP parallel private(iatm,Avaltmp,Bvaltmp,tmprow,mat11,Qij,Qii,Qjj) num_threads(nthreads)
					Avaltmp=0
					Bvaltmp=0
                    !$OMP do schedule(dynamic)
					do iatm=1,ncenter !Note: Constructing tmprow is the most time-consuming step of PM-Becke method, especially for large system
                        tmprow=matmul_blas(veccoeffi(1:1,:),AOMbas(:,:,iatm),1,nbasis,0,0)
                        mat11=matmul_blas(tmprow(1:1,:),Cmat(:,jmo:jmo),1,1,0,0)
                        Qij=mat11(1,1)
                        mat11=matmul_blas(tmprow(1:1,:),Cmat(:,imo:imo),1,1,0,0)
                        Qii=mat11(1,1)
                        tmprow=matmul_blas(veccoeffj(1:1,:),AOMbas(:,:,iatm),1,nbasis,0,0)
                        mat11=matmul_blas(tmprow(1:1,:),Cmat(:,jmo:jmo),1,1,0,0)
                        Qjj=mat11(1,1)
						if (iPMexp==2) then
							Avaltmp=Avaltmp+( Qij**2-(Qii-Qjj)**2/4D0 )
							Bvaltmp=Bvaltmp+( Qij*(Qii-Qjj) )
						else if (iPMexp==4) then
							Avaltmp=Avaltmp-Qii**4-Qjj**4+6*(Qii**2+Qjj**2)*Qij**2+Qii**3*Qjj+Qjj**3*Qii
							Bvaltmp=Bvaltmp+4*Qij*(Qii**3-Qjj**3)
						end if
					end do
                    !$OMP end do
                    !$OMP CRITICAL
	                Aval=Aval+Avaltmp
	                Bval=Bval+Bvaltmp
                    !$OMP end CRITICAL
                    !$OMP end parallel
                !Foster-Boys localization
				else if (imethod==10) then
					call boysdipint(iri,jrj,irj,imo,jmo,Cmat)
					Aval=sum(irj**2)-sum((iri-jrj)**2)/4
					Bval=sum(irj*(iri-jrj))
				end if
                
				if (Aval**2+Bval**2<1D-12) cycle
				gamma=sign(1D0,Bval)*acos(-Aval/dsqrt(Aval**2+Bval**2))/4D0
				arrayi=cos(gamma)*Cmat(:,imo)+sin(gamma)*Cmat(:,jmo) !Rotate imo and jmo by angle gamma to yield a new orbital
				arrayj=-sin(gamma)*Cmat(:,imo)+cos(gamma)*Cmat(:,jmo)
				Cmat(:,imo)=arrayi !Update current coefficient matrix by new mixed orbital
				Cmat(:,jmo)=arrayj
				if (imethod==1) then !For PM-Mulliken, also correspondingly update auxiliary matrix SC
					arrayi=cos(gamma)*SC(:,imo)+sin(gamma)*SC(:,jmo)
					arrayj=-sin(gamma)*SC(:,imo)+cos(gamma)*SC(:,jmo)
					SC(:,imo)=arrayi
					SC(:,jmo)=arrayj
				end if
                
			end do
		end do
        
		Pval=0
		do idx=1,norblist
            imo=orblist(idx)
			do iatm=1,ncenter
				if (basstart(iatm)==0) cycle
				Pval=Pval+Qval(Cmat,imo,imo,iatm)**iPMexp
			end do
		end do
		deltaPval=Pval-Pvalold
        
		write(*,"(' Cycle:',i5,'  P:',f16.8,'  Delta P:',f16.8)") icyc,Pval,deltaPval
		if (abs(deltaPval)<crit) exit
		Pvalold=Pval
	end do
    
	if (icyc==maxcyc+1) then
		write(*,*) "Warning: Convergence failed!"
	else
		write(*,"(a)") " Successfully converged!"
	end if
end do

if (imethod==2) then !PM with Lowdin. Back convert CObas from orthonormal basis to original basis
    CObasa=matmul_blas(Xmatinv,CObasa,nbasis,nbasis,0,0)
	if (allocated(CObasb)) CObasb=matmul_blas(Xmatinv,CObasb,nbasis,nbasis,0,0)
	Sbas=Sbas_org
end if

call walltime(iwalltime2)
write(*,"(/,' Orbital localization took up wall clock time',i10,' s')") iwalltime2-iwalltime1


!Print orbital energies, sort orbitals according to energies
if (idoene==1) then
	write(*,*) "Evaluating orbital energies..."
	nmobeg=1
	if (isel/=3.and.idocore==0) nmobeg=ninnerele/2+1
	!Do Alpha part or closed-shell orbitals
	allocate(FLMOA(nbasis,nbasis),tmpmat(nbasis,nbasis))
    tmpmat=matmul_blas(CObasa,FmatA,nbasis,nbasis,1,0)
    FLMOA=matmul_blas(tmpmat,CObasa,nbasis,nbasis,0,0)
	!FLMOA=matmul(matmul(transpose(CObasa),FmatA),CObasa)
	if (isel==1) then
        nmoend=naelec
	else if (isel==2.or.isel==3) then
        nmoend=nbasis
    end if
	do iorb=nmobeg,nmoend
		MOene(iorb)=FLMOA(iorb,iorb)
	end do
	do iorb=nmobeg,nmoend
		do jorb=iorb+1,nmoend
			if (MOene(iorb)>MOene(jorb)) then	
				tmpbasarr=CObasa(:,iorb)
				CObasa(:,iorb)=CObasa(:,jorb)
				CObasa(:,jorb)=tmpbasarr
				tmpprimarr=CO(iorb,:)
				CO(iorb,:)=CO(jorb,:)
				CO(jorb,:)=tmpprimarr
				tmpene=MOene(iorb)
				MOene(iorb)=MOene(jorb)
				MOene(jorb)=tmpene
			end if
		end do
	end do
    write(*,*)
	write(*,*) "Energies of localized orbitals:"
	do iorb=nmobeg,nmoend
        if (orblist_did(iorb)==0) cycle
		typestr="A+B"
		if (wfntype==1)	typestr="A"
		write(*,"(i6,'   Energy:',f13.7,' a.u.',f13.4,' eV   Type: ',a,'  Occ:',f4.1)") &
		iorb,MOene(iorb),MOene(iorb)*au2eV,typestr,MOocc(iorb)
	end do
	!Do beta part
	if (wfntype==1) then
		allocate(FLMOB(nbasis,nbasis))
		tmpmat=matmul_blas(CObasb,FmatB,nbasis,nbasis,1,0)
		FLMOB=matmul_blas(tmpmat,CObasb,nbasis,nbasis,0,0)
		!FLMOB=matmul(matmul(transpose(CObasb),FmatB),CObasb)
		if (isel==1) then
            nmoend=nbelec
		else if (isel==2.or.isel==3) then
            nmoend=nbasis
        end if
		do iorb=nmobeg,nmoend
			MOene(nbasis+iorb)=FLMOB(iorb,iorb)
		end do
		do iorb=nmobeg,nmoend
			do jorb=iorb+1,nmoend
				if (MOene(nbasis+iorb)>MOene(nbasis+jorb)) then	
					tmpbasarr=CObasb(:,iorb)
					CObasb(:,iorb)=CObasb(:,jorb)
					CObasb(:,jorb)=tmpbasarr
					tmpprimarr=CO(nbasis+iorb,:)
					CO(nbasis+iorb,:)=CO(nbasis+jorb,:)
					CO(nbasis+jorb,:)=tmpprimarr
					tmpene=MOene(nbasis+iorb)
					MOene(nbasis+iorb)=MOene(nbasis+jorb)
					MOene(nbasis+jorb)=tmpene
				end if
			end do
		end do
		typestr="B"
		do iorb=nmobeg,nmoend
            if (all(orblist(1:norblist)/=iorb)) cycle
			write(*,"(i6,'   Energy:',f13.7,' a.u.',f13.4,' eV   Type: ',a,'  Occ:',f4.1)") &
			iorb+nbasis,MOene(iorb+nbasis),MOene(iorb+nbasis)*au2eV,typestr,MOocc(iorb+nbasis)
		end do
	end if
	if (isel/=3.and.idocore==0)  write(*,*) "Energies of core orbitals are not updated since they were not localized"
	if (isel==1) write(*,*) "Energies of unoccupied orbitals are not updated since they were not localized"
	
	!Second-order perturbation analysis between occupied and virtual orbitals (like NBO E2), this only works when both of them have been localized
	!This part is commented since it don't print any useful result, because it is easy to proved that Fock element between occupied and virtual LMOs are exactly zero
! 	if (isel==1) then
! 		write(*,*) " Note: E(2) analysis is skipped since virtual orbitals were not localized"
! 	else if (isel==2) then
! 		write(*,*) "Second-order perturbation theory analysis of interaction energy:"
! 		!Regenerated Fock matrix in LMO, since they have been sorted
! 		FLMOA=matmul(matmul(transpose(CObasa),FmatA),CObasa)
! 		call showmatgau(FLMOA)
! 		read(*,*)
! 		coeff=2
! 		do iocc=1,naelec
! 			do ivir=naelec+1,nbasis
! 				E2val=-coeff* FLMOA(iocc,ivir)**2/(MOene(ivir)-MOene(iocc))
! 				qCT=coeff* ( FLMOA(iocc,ivir)/(MOene(ivir)-MOene(iocc)) )**2
! 				if (abs(E2val*au2kcal)>0.2D0) then
! 					write(*,"(' Donor:',i5,'  -  Acceptor:',i5,'   E(2):',f7.2,' kcal/mol   q_CT:',f10.5)") iocc,ivir,E2val*au2kcal,qCT
! 					write(*,"(3f16.10)") FLMOA(iocc,ivir),MOene(ivir)-MOene(iocc)
! 				end if
! 			end do
! 		end do
! 	end if
end if


!Calculate orbital composition and print major character of LMOs
if (icompmethod>0) then
	orbtype=0
	!At most four batches, 1: Alpha-occ 2: Alpha-vir 3: Beta-occ 4: Beta-vir
	do itime=1,4
        if (domark(itime)==0) cycle
        if (isel==1.or.isel==2) then
		    if (itime<=2) then !For batches 1 and 2, use alpha density matrix
			    Cmat=>CObasa
		    else !For batches 3 and 4, use beta density matrix
			    Cmat=>CObasb
		    end if
		    if (itime==1) then
			    ibeg=1
			    if (idocore==0) ibeg=ninnerele/2+1
			    iend=naelec
		    else if (itime==2) then
			    ibeg=naelec+1
			    iend=nbasis
		    else if (itime==3) then
			    ibeg=1
			    if (idocore==0) ibeg=ninnerele/2+1
			    iend=nbelec
		    else if (itime==4) then
			    ibeg=nbelec+1
			    iend=nbasis
		    end if
            norblist=iend-ibeg+1
            forall(i=1:norblist) orblist(i)=ibeg-1+i
        else if (isel==3) then
            !orblist has already been set, however ibeg and iend should be set since will be used later
            ibeg=minval(orblist(1:norblist))
            iend=maxval(orblist(1:norblist))
        end if
		
		!Calculate orbital composition and then sort from large to small. Size: orbcomp(1:ncenter,1:nbasis)
        write(*,*) "Calculating orbital compositions..."
        if (icompmethod==1) then !Mulliken+SCPA
			orbcomp=0
		    do idx=1,norblist
                iorb=orblist(idx)
			    if (itime==1.or.itime==3) then !For occupied LMOs, use Mulliken
				    do ibas=1,nbasis
					    bascross=0D0
					    do jbas=1,nbasis
						    if (jbas==ibas) cycle
						    bascross=bascross+Cmat(ibas,iorb)*Cmat(jbas,iorb)*Sbas(ibas,jbas)
					    end do
					    bastot(ibas)=Cmat(ibas,iorb)**2+bascross
				    end do
				    do iatm=1,ncenter
						if (basstart(iatm)==0) cycle
					    orbcomp(iatm,iorb)=sum(bastot(basstart(iatm):basend(iatm)))
				    end do
			    else !For virtual LMOs, use SCPA. Since it guarantees that the result is within 0~100%
 				    do iatm=1,ncenter
						if (basstart(iatm)==0) cycle
 					    orbcomp(iatm,iorb)=sum(Cmat(basstart(iatm):basend(iatm),iorb)**2)
 				    end do
 				    orbcomp(:,iorb)=orbcomp(:,iorb)/sum(orbcomp(:,iorb))
			    end if
            end do
        else if (icompmethod==2.or.icompmethod==3) then
            write(*,*)
            if (icompmethod==2) itmp=1 !Hirshfeld
            if (icompmethod==3) itmp=2 !Becke
            allocate(CO_tmp(nmo,nprims))
            CO_tmp=CO
            if (itime==1.or.itime==2) then !Alpha
                call CObas2CO(1) !Convert current CObasa to CO, which is needed by Hirshfeld/Becke method
                call gen_orbatmcomp_space(itmp,orbcomp(:,ibeg:iend),ibeg,iend,0,0)
            else if (itime==3.or.itime==4) then !Beta
                call CObas2CO(2) !Convert current CObasb to CO, which is needed by Hirshfeld/Becke method
                call gen_orbatmcomp_space(itmp,orbcomp(:,ibeg:iend),ibeg+nbasis,iend+nbasis,0,0) !Use cheap grid
            end if
            CO=CO_tmp
            deallocate(CO_tmp)
        end if
        
        !Sort atomic contributions
		do idx=1,norblist
            iorb=orblist(idx)
			forall(iatm=1:ncenter) iatmarr(iatm,iorb)=iatm
			call sort(orbcomp(:,iorb),"val",iatmarr(:,iorb)) !Sort atomic contributions from small to large
			call invarr(orbcomp(:,iorb),iatmarr(:,iorb)) !Then become from large to small
		end do
		
        if (itime==1) write(*,"(/,a)") " Hint: If you hope to print the LMOs in the order of atoms and atomic pairs, &
        &set ""iprintLMOorder"" in settings.ini to 1 prior to the analysis"
		if (wfntype==0) then
			if (itime==1) write(*,"(/,a)") " **** Major character of occupied LMOs:"
			if (itime==2) write(*,"(/,a)") " **** Major character of unoccupied LMOs:"
		else if (wfntype==1) then
			if (itime==1) write(*,"(/,a)") " **** Major character of alpha occupied LMOs:"
			if (itime==2) write(*,"(/,a)") " **** Major character of alpha unoccupied LMOs:"
			if (itime==3) write(*,"(/,a)") " **** Major character of beta occupied LMOs:"
			if (itime==4) write(*,"(/,a)") " **** Major character of beta unoccupied LMOs:"
		end if
		istatarr=0 !If =1, then the character of the LMO has been identified and printed
		write(*,"(' Almost single center LMOs: (An atom has contribution >',f5.1,'%)')") crit1c*100
		itmp=0
        if (iprintLMOorder==0) then !Print LMO in original order
		    do idx=1,norblist
                iorb=orblist(idx)
			    if (orbcomp(1,iorb)>crit1c) then
				    write(*,"(i5,':',i4,'(',a,')',f5.1,'%    ')",advance='no') iorb,iatmarr(1,iorb),a(iatmarr(1,iorb))%name,orbcomp(1,iorb)*100
				    itmp=itmp+1
				    if (mod(itmp,3)==0) write(*,*)
				    istatarr(iorb)=1
				    if (itime<=2) then
					    orbtype(iorb)=1
					    orbatm(1,iorb)=iatmarr(1,iorb)
				    else
					    orbtype(iorb+nbasis)=1
					    orbatm(1,iorb+nbasis)=iatmarr(1,iorb)
				    end if
			    end if
		    end do
        else if (iprintLMOorder==1) then !Print LMO in the sequence of major contribution atom
            do iatm=1,ncenter
		        do idx=1,norblist
                    iorb=orblist(idx)
			        if (orbcomp(1,iorb)>crit1c.and.iatmarr(1,iorb)==iatm) then
				        write(*,"(i5,':',i4,'(',a,')',f5.1,'%    ')",advance='no') iorb,iatmarr(1,iorb),a(iatmarr(1,iorb))%name,orbcomp(1,iorb)*100
				        itmp=itmp+1
				        if (mod(itmp,3)==0) write(*,*)
				        istatarr(iorb)=1
				        if (itime<=2) then
					        orbtype(iorb)=1
					        orbatm(1,iorb)=iatmarr(1,iorb)
				        else
					        orbtype(iorb+nbasis)=1
					        orbatm(1,iorb+nbasis)=iatmarr(1,iorb)
				        end if
			        end if
		        end do
            end do
        end if
		if (mod(itmp,3)/=0) write(*,*)
		if (itmp==0) write(*,*) "None!"
		
		if (all(istatarr(orblist(1:norblist))==1)) cycle
		write(*,*)
		write(*,"(' Almost two-center LMOs: (Sum of two largest contributions >',f5.1,'%)')") crit2c*100
		itmp=0
        if (iprintLMOorder==0) then !Print LMO in original order
		    do idx=1,norblist
                iorb=orblist(idx)
			    if (istatarr(iorb)==1) cycle
			    if (orbcomp(1,iorb)+orbcomp(2,iorb)>crit2c) then
				    write(*,"(i5,':',i4,'(',a,')',f5.1,'%',i4,'(',a,')',f5.1,'%     ')",advance='no') &
				    iorb,iatmarr(1,iorb),a(iatmarr(1,iorb))%name,orbcomp(1,iorb)*100,&
				    iatmarr(2,iorb),a(iatmarr(2,iorb))%name,orbcomp(2,iorb)*100
				    itmp=itmp+1
				    if (mod(itmp,2)==0) write(*,*)
				    istatarr(iorb)=1
				    if (itime<=2) then
					    orbtype(iorb)=2
					    orbatm(:,iorb)=iatmarr(1:2,iorb)
				    else
					    orbtype(iorb+nbasis)=2
					    orbatm(:,iorb+nbasis)=iatmarr(1:2,iorb)
				    end if
			    end if
		    end do
        else if (iprintLMOorder==1) then !Print LMO in the sequence of major contribution atomic pair
            do iatm=1,ncenter
                do jatm=iatm+1,ncenter
		            do idx=1,norblist
                        iorb=orblist(idx)
			            if (istatarr(iorb)==1) cycle
			            if (orbcomp(1,iorb)+orbcomp(2,iorb)>crit2c) then
                            if (iatmarr(1,iorb)==iatm.and.iatmarr(2,iorb)==jatm) then
				                write(*,"(i5,':',i4,'(',a,')',f5.1,'%',i4,'(',a,')',f5.1,'%     ')",advance='no') &
				                iorb,iatmarr(1,iorb),a(iatmarr(1,iorb))%name,orbcomp(1,iorb)*100,&
				                iatmarr(2,iorb),a(iatmarr(2,iorb))%name,orbcomp(2,iorb)*100
                            else if (iatmarr(1,iorb)==jatm.and.iatmarr(2,iorb)==iatm) then
				                write(*,"(i5,':',i4,'(',a,')',f5.1,'%',i4,'(',a,')',f5.1,'%     ')",advance='no') &
				                iorb,iatmarr(2,iorb),a(iatmarr(2,iorb))%name,orbcomp(2,iorb)*100,&
				                iatmarr(1,iorb),a(iatmarr(1,iorb))%name,orbcomp(1,iorb)*100
                            else
                                cycle
                            end if
				            itmp=itmp+1
				            if (mod(itmp,2)==0) write(*,*)
				            istatarr(iorb)=1
				            if (itime<=2) then
					            orbtype(iorb)=2
					            orbatm(:,iorb)=iatmarr(1:2,iorb)
				            else
					            orbtype(iorb+nbasis)=2
					            orbatm(:,iorb+nbasis)=iatmarr(1:2,iorb)
				            end if
			            end if
		            end do
                end do
            end do
        end if
		if (mod(itmp,2)/=0) write(*,*)
		if (itmp==0) write(*,*) "None!"
		
		if (all(istatarr(orblist(1:norblist))==1)) cycle
		write(*,*)
		write(*,*) "More delocalized LMOs: (Three largest contributions are printed)"
		do idx=1,norblist
            iorb=orblist(idx)
			if (istatarr(iorb)==1) cycle
			write(*,"(i5,':',i5,'(',a,')',f5.1,'%',i5,'(',a,')',f5.1,'%',i5,'(',a,')',f5.1,'%')") &
			iorb,iatmarr(1,iorb),a(iatmarr(1,iorb))%name,orbcomp(1,iorb)*100,&
			iatmarr(2,iorb),a(iatmarr(2,iorb))%name,orbcomp(2,iorb)*100,&
			iatmarr(3,iorb),a(iatmarr(3,iorb))%name,orbcomp(3,iorb)*100
		end do
	end do
end if

write(*,*)
write(*,*) "Exporting localized orbitals to new.fch in current folder..."
call outfch("new.fch",10,1)

if (ireload==1) then !Automatically reload the newly generated new.fch as requested
	call dealloall(0)
	write(*,*) "Loading new.fch..."
    call readinfile("new.fch",1)
	write(*,*) "Loading finished!"

	!Adding center of LMOs as Bq atoms
	if (ilmocen==1) then
        if (isel==3) then
            write(*,*) "Note: Calculating center of LMOs is not supported in current case"
            return
        end if
        write(*,*) "Generating electric dipole moment integral matrix..."
        call genDbas_curr
		write(*,*)
		write(*,*) "Calculating center of LMOs and meanwhile adding them as Bq atoms..."
		
		!Backup a to a_tmp, then add Bq to a_tmp during calculation, and finally copy a_tmp to a
		if (isel==1) then
			ncenter_new=ncenter+nint(naelec)
			if (wfntype==1) ncenter_new=ncenter_new+nint(nbelec)
		else if (isel==2) then
			ncenter_new=ncenter+nbasis
			if (wfntype==1) ncenter_new=ncenter_new+nbasis
		end if
		if (idocore==0) then
			if (wfntype==0) ncenter_new=ncenter_new-ninnerele/2
			if (wfntype==1) ncenter_new=ncenter_new-ninnerele
		end if
		if (allocated(a_tmp)) deallocate(a_tmp)
		allocate(a_tmp(ncenter_new))
		a_tmp(1:ncenter)=a
		itmp=ncenter
		
		LMOpos=0
		open(10,file="LMOcen.txt",status="replace")
		if (wfntype==0) ntime=1
		if (wfntype==1) ntime=2
		do itime=1,ntime !=1: Total or alpha, =2: Beta
			if (itime==1) then
				ibeg=1
				if (idocore==0) ibeg=ninnerele/2+1
				iend=nint(naelec)
				if (isel==2) iend=nbasis
			else if (itime==2) then
				ibeg=nbasis+1
				if (idocore==0) ibeg=nbasis+ninnerele/2+1
				iend=nbasis+nint(nbelec)
				if (isel==2) iend=2*nbasis
			end if
			do iorb=ibeg,iend
				!The Dbas contains negative sign of operator r, therefore we need to use "minus"
				if (itime==1) then
					do ibas=1,nbasis
						do jbas=1,nbasis
							LMOpos(:,iorb)=LMOpos(:,iorb)-Dbas(:,ibas,jbas)*CObasa(ibas,iorb)*CObasa(jbas,iorb)
						end do
					end do
				else
					do ibas=1,nbasis
						do jbas=1,nbasis
							LMOpos(:,iorb)=LMOpos(:,iorb)-Dbas(:,ibas,jbas)*CObasb(ibas,iorb-nbasis)*CObasb(jbas,iorb-nbasis)
						end do
					end do
				end if
				itmp=itmp+1
				a_tmp(itmp)%x=LMOpos(1,iorb)
				a_tmp(itmp)%y=LMOpos(2,iorb)
				a_tmp(itmp)%z=LMOpos(3,iorb)
				a_tmp(itmp)%index=0
				a_tmp(itmp)%charge=0
				a_tmp(itmp)%name="Bq"
				if (itime==1) then
					if (wfntype==0) write(10,"(' LMO',i6,' corresponds to Bq',i6,', X,Y,Z:',3f10.4,' Bohr')") iorb,itmp,LMOpos(:,iorb)
					if (wfntype==1) write(10,"(' Alpha LMO',i6,': Bq',i6,', X,Y,Z:',3f10.4,' Bohr')") iorb,itmp,LMOpos(:,iorb)
				else
					write(10,"(' Beta LMO ',i6,': Bq',i6,', X,Y,Z:',3f10.4,' Bohr')") iorb-nbasis,itmp,LMOpos(:,iorb)
				end if
			end do
		end do
		close(10)
		
		ncenter=ncenter_new
		deallocate(a)
		allocate(a(ncenter))
		a=a_tmp
		deallocate(a_tmp)
		!Set plotting parameters in main function 0 to the best status for showing LMO centers
		iatmlabtype3D=4
		ishowaxis=0
		bondradius=0.06D0
		ratioatmsphere=0.6D0
		iatmlabtype3D=5
		isosur1style=2
		isosur2style=2
		atmlabclrR=1D0;atmlabclrG=0D0;atmlabclrB=0D0
		write(*,"(/,a)") " Done! The Bq atoms in current system now correspond to center of LMOs. The LMO center coordinates, &
		&correspondence between LMO indices and Bq indices have been exported to LMOcen.txt in current folder. In addition, &
		&the plotting parameters in main function 0 have been set to the best status for showing LMO centers"
		write(*,"(a)") " Note: Since these Bq atoms do not have corresponding basis functions, &
		&the present wavefunction should not be subjected to wavefunction analyses, otherwise Multiwfn may crash"
		write(*,*)
		
		
		!---- Calculate dipole moment LMOs
		!LMOpos is <phi|r|phi>, therefore its negative is LMO electronic contribution to dipole moment
		write(*,*) "Would you like to perform dipole moment analysis for occupied LMOs? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			!Evaluate system total dipole moment
			x_nuc=sum(a%charge*a%x);y_nuc=sum(a%charge*a%y);z_nuc=sum(a%charge*a%z)
			x_ele=0;y_ele=0;z_ele=0
			if (wfntype==0) then
				occval=2D0
				do iLMO=1,nint(naelec)
					x_ele=x_ele-LMOpos(1,iLMO);y_ele=y_ele-LMOpos(2,iLMO);z_ele=z_ele-LMOpos(3,iLMO)
				end do
				x_ele=x_ele*occval;y_ele=y_ele*occval;z_ele=z_ele*occval
			else if (wfntype==1) then
				occval=1D0
				do iLMO=1,nint(naelec)
					x_ele=x_ele-LMOpos(1,iLMO);y_ele=y_ele-LMOpos(2,iLMO);z_ele=z_ele-LMOpos(3,iLMO)
				end do
				do iLMO=1,nint(nbelec)
					x_ele=x_ele-LMOpos(1,iLMO+nbasis);y_ele=y_ele-LMOpos(2,iLMO+nbasis);z_ele=z_ele-LMOpos(3,iLMO+nbasis)
				end do
			end if
			x_tot=x_nuc+x_ele;y_tot=y_nuc+y_ele;z_tot=z_nuc+z_ele
			
			open(10,file="LMOdip.txt",status="replace")
			write(10,"(' Dipole moment of the whole system (including nuclear and electronic contributions)',/,' X/Y/Z:',3f10.5,'  Norm:',f10.5,' a.u.')") &
			x_tot,y_tot,z_tot,dsqrt(x_tot**2+y_tot**2+z_tot**2)
			write(10,"(' Nuclear contribution:',/,' X/Y/Z:',3f10.5,'  Norm:',f10.5,' a.u.')") &
			x_nuc,y_nuc,z_nuc,dsqrt(x_nuc**2+y_nuc**2+z_nuc**2)
			write(10,"(' Electronic contribution:',/,' X/Y/Z:',3f10.5,'  Norm:',f10.5,' a.u.',/)") &
			x_ele,y_ele,z_ele,dsqrt(x_ele**2+y_ele**2+z_ele**2)
			!Total or alpha part
			if (wfntype==1) then
				write(10,*) "===================================="
				write(10,*) "============ Alpha LMOs ============"
				write(10,*) "===================================="
				write(10,*)
			end if
			write(10,*) "Single-center orbital dipole moments (a.u.):"
			tmpbx=0;tmpby=0;tmpbz=0
			do iLMO=1,nint(naelec)
				if (orbtype(iLMO)==1) then
					iatm=orbatm(1,iLMO)
					xdip=(a(iatm)%x-LMOpos(1,iLMO)) *occval
					ydip=(a(iatm)%y-LMOpos(2,iLMO)) *occval
					zdip=(a(iatm)%z-LMOpos(3,iLMO)) *occval
					write(10,"(i5,' (',i4,a,')  X/Y/Z:',3f10.5,'  Norm:',f10.5)") &
					iLMO,iatm,a(iatm)%name,xdip,ydip,zdip,dsqrt(xdip**2+ydip**2+zdip**2)
					tmpbx=tmpbx+xdip;tmpby=tmpby+ydip;tmpbz=tmpbz+zdip
				end if
			end do
			write(10,"(' Sum            X/Y/Z:',3f10.5,'  Norm:',f10.5)") tmpbx,tmpby,tmpbz,dsqrt(tmpbx**2+tmpby**2+tmpbz**2)
			write(10,*)
			!If you want to visualize the positive charge center, uncomment below lines
			!allocate(a_tmp(ncenter+count(orbtype==2)))
			!a_tmp(1:ncenter)=a
			write(10,*) "Two-center bond dipole moments (a.u.):"
			tmpbx=0;tmpby=0;tmpbz=0
			do iLMO=1,nint(naelec)
				if (orbtype(iLMO)==2) then
					iatm=orbatm(1,iLMO)
					jatm=orbatm(2,iLMO)
					ratioi=covr(a(iatm)%index) !Adjust positive charge center according to covalent radii
					ratioj=covr(a(jatm)%index)
					!ratioi=0.5D0 !Use bond midpoint as positive charge center, as the "dipole" keyword of NBO used
					!ratioj=0.5D0
					tmpsum=ratioi+ratioj
					ratioi=ratioi/tmpsum
					ratioj=ratioj/tmpsum
					posx=a(iatm)%x*ratioj+a(jatm)%x*ratioi
					posy=a(iatm)%y*ratioj+a(jatm)%y*ratioi
					posz=a(iatm)%z*ratioj+a(jatm)%z*ratioi
					xdip=occval*posx-occval*LMOpos(1,iLMO)
					ydip=occval*posy-occval*LMOpos(2,iLMO)
					zdip=occval*posz-occval*LMOpos(3,iLMO)
					!ncenter=ncenter+1
					!a_tmp(ncenter)%x=posx;a_tmp(ncenter)%y=posy;a_tmp(ncenter)%z=posz
					!a_tmp(ncenter)%index=0;a_tmp(ncenter)%name="Bq";a_tmp(ncenter)%charge=0
					!write(*,"(3i5,2f6.3,3f10.4)") ncenter,iatm,jatm,ratioi,ratioj,posx,posy,posz
					!write(*,"(3f10.4)") -LMOpos(:,iLMO)
					write(10,"(i5,' (',i4,a,' -',i4,a,')  X/Y/Z:',3f10.5,'  Norm:',f10.5)") &
					iLMO,iatm,a(iatm)%name,jatm,a(jatm)%name,xdip,ydip,zdip,dsqrt(xdip**2+ydip**2+zdip**2)
					tmpbx=tmpbx+xdip;tmpby=tmpby+ydip;tmpbz=tmpbz+zdip
				end if
			end do
			write(10,"(' Sum                    X/Y/Z:',3f10.5,'  Norm:',f10.5)") tmpbx,tmpby,tmpbz,dsqrt(tmpbx**2+tmpby**2+tmpbz**2)
			!deallocate(a)
			!allocate(a(ncenter))
			!a=a_tmp
			!deallocate(a_tmp)
			if (any(orbtype(1:nint(naelec))==0)) then
				write(10,"(/,a)") " Some more delocalized LMOs exist, but ignored here" 
			end if
			write(10,*)
			write(10,*) "Contributions of all occupied LMOs to system dipole moment (a.u.):"
			do iLMO=1,nint(naelec)
				write(10,"(i5,'  X/Y/Z:',3f10.5,'  Norm:',f10.5)") iLMO,-LMOpos(:,iLMO)*occval,dsqrt(sum(LMOpos(:,iLMO)**2))*occval
			end do
			x_LMOs=sum(-LMOpos(1,1:nint(naelec)))*occval
			y_LMOs=sum(-LMOpos(2,1:nint(naelec)))*occval
			z_LMOs=sum(-LMOpos(3,1:nint(naelec)))*occval
			write(10,"(' Sum   X/Y/Z:',3f10.5,'  Norm:',f10.5)") x_LMOs,y_LMOs,z_LMOs,dsqrt(x_LMOs**2+y_LMOs**2+z_LMOs**2)
			
			!Beta part, simply repeat above codes
			if (wfntype==1) then
				write(10,*)
				write(10,*) "===================================="
				write(10,*) "============ Beta LMOs ============="
				write(10,*) "===================================="
				write(10,*)
				write(10,*) "Single-center orbital dipole moments (a.u.):"
				tmpbx=0;tmpby=0;tmpbz=0
				do iLMO=1,nint(nbelec)
					if (orbtype(iLMO+nbasis)==1) then
						iatm=orbatm(1,iLMO+nbasis)
						xdip=(a(iatm)%x-LMOpos(1,iLMO+nbasis)) *occval
						ydip=(a(iatm)%y-LMOpos(2,iLMO+nbasis)) *occval
						zdip=(a(iatm)%z-LMOpos(3,iLMO+nbasis)) *occval
						write(10,"(i5,' (',i4,a,')  X/Y/Z:',3f10.5,'  Norm:',f10.5)") &
						iLMO,iatm,a(iatm)%name,xdip,ydip,zdip,dsqrt(xdip**2+ydip**2+zdip**2)
						tmpbx=tmpbx+xdip;tmpby=tmpby+ydip;tmpbz=tmpbz+zdip
					end if
				end do
				write(10,"(' Sum            X/Y/Z:',3f10.5,'  Norm:',f10.5)") tmpbx,tmpby,tmpbz,dsqrt(tmpbx**2+tmpby**2+tmpbz**2)
				write(10,*)
				write(10,*) "Two-center bond dipole	moments (a.u.):"
				tmpbx=0;tmpby=0;tmpbz=0
				do iLMO=1,nint(nbelec)
					if (orbtype(iLMO+nbasis)==2) then
						iatm=orbatm(1,iLMO+nbasis)
						jatm=orbatm(2,iLMO+nbasis)
						ratioi=covr(a(iatm)%index) !Adjust positive charge center according to covalent radii
						ratioj=covr(a(jatm)%index)
						!ratioi=0.5D0 !Use bond midpoint as positive charge center, this is the "dipole" keyword of NBO used
						!ratioj=0.5D0
						tmpsum=ratioi+ratioj;ratioi=ratioi/tmpsum;ratioj=ratioj/tmpsum
						posx=a(iatm)%x*ratioj+a(jatm)%x*ratioi
						posy=a(iatm)%y*ratioj+a(jatm)%y*ratioi
						posz=a(iatm)%z*ratioj+a(jatm)%z*ratioi
						xdip=occval*posx-occval*LMOpos(1,iLMO+nbasis)
						ydip=occval*posy-occval*LMOpos(2,iLMO+nbasis)
						zdip=occval*posz-occval*LMOpos(3,iLMO+nbasis)
						write(10,"(i5,' (',i4,a,' -',i4,a,')  X/Y/Z:',3f10.5,'  Norm:',f10.5)") &
						iLMO,iatm,a(iatm)%name,jatm,a(jatm)%name,xdip,ydip,zdip,dsqrt(xdip**2+ydip**2+zdip**2)
						tmpbx=tmpbx+xdip;tmpby=tmpby+ydip;tmpbz=tmpbz+zdip
					end if
				end do
				write(10,"(' Sum                    X/Y/Z:',3f10.5,'  Norm:',f10.5)") tmpbx,tmpby,tmpbz,dsqrt(tmpbx**2+tmpby**2+tmpbz**2)
				if (any(orbtype(1:nint(naelec))==0)) then
					write(10,"(/,a)") " Some more delocalized LMOs exist, but ignored here" 
				end if
				write(10,*)
				write(10,*) "Contributions of all occupied LMOs to system dipole moment (a.u.):"
				do iLMO=1,nint(nbelec)
					write(10,"(i5,'  X/Y/Z:',3f10.5,'  Norm:',f10.5)") iLMO,&
					-LMOpos(:,iLMO+nbasis)*occval,dsqrt(sum(LMOpos(:,iLMO+nbasis)**2))*occval
				end do
				x_LMOs=sum(-LMOpos(1,1:nint(nbelec)))*occval
				y_LMOs=sum(-LMOpos(2,1:nint(nbelec)))*occval
				z_LMOs=sum(-LMOpos(3,1:nint(nbelec)))*occval
				write(10,"(' Sum   X/Y/Z:',3f10.5,'  Norm:',f10.5)") x_LMOs,y_LMOs,z_LMOs,dsqrt(x_LMOs**2+y_LMOs**2+z_LMOs**2)
			end if
			
			close(10)
			write(*,"(a)") " Done! Dipole moment of occupied LMOs as well as their contribution to &
			&system dipole moment have been exported as LMOdip.txt in current folder"
		end if
	end if
end if
end subroutine



!------ Calculate Q value used for determining localization convergence
real*8 function Qval(Cmat,imo,jmo,iatm)
use defvar
real*8 Cmat(nbasis,nbasis)
integer imo,jmo,iatm,ibas
Qval=0
do ibas=basstart(iatm),basend(iatm)
	Qval=Qval+sum( (Cmat(:,imo)*Cmat(ibas,jmo)+Cmat(ibas,imo)*Cmat(:,jmo)) *Sbas(ibas,:) )
end do
Qval=Qval/2D0
end function


!------ Calculate dipole moment integral between orbitals that involved in Boys localization
subroutine boysdipint(iri,jrj,irj,imo,jmo,Cmat)
use defvar
implicit real*8 (a-h,o-z)
integer imo,jmo
real*8 :: iri(3),jrj(3),irj(3),iripriv(3),jrjpriv(3),irjpriv(3),Cmat(nbasis,nbasis)
iri=0
jrj=0
irj=0
!$OMP parallel shared(iri,jrj,irj) private(iripriv,jrjpriv,irjpriv) num_threads(nthreads)
iripriv=0
jrjpriv=0
irjpriv=0
!$OMP do schedule(DYNAMIC)
do ibas=1,nbasis
	do jbas=1,nbasis
		iripriv=iripriv+Dbas(:,ibas,jbas)*Cmat(ibas,imo)*Cmat(jbas,imo)
		jrjpriv=jrjpriv+Dbas(:,ibas,jbas)*Cmat(ibas,jmo)*Cmat(jbas,jmo)
		irjpriv=irjpriv+Dbas(:,ibas,jbas)*Cmat(ibas,imo)*Cmat(jbas,jmo)
	end do
end do
!$OMP END DO
!$OMP CRITICAL
iri=iri+iripriv
jrj=jrj+jrjpriv
irj=irj+irjpriv
!$OMP END CRITICAL
!$OMP END PARALLEL
end subroutine