!-------- Main interface of various other functions part 2
subroutine otherfunc2_main
implicit real*8 (a-h,o-z)
do while(.true.)
	write(*,*)
	write(*,*) "              ============ Other functions (Part 2) ============"
	write(*,*) "0 Return"
	write(*,*) "1 Calculate core-valence bifurcation (CVB) index and related quantities"
	write(*,*) "2 Calculate atomic and bond dipole moments in Hilbert space"
	write(*,*) "3 Generate cube file for multiple orbital wavefunctions"
	!write(*,*) "4 Generate iso-chemical shielding surfaces (ICSS) and related quantities"
	write(*,*) "5 Plot radial distribution function for a real space function"
	write(*,*) "6 Analyze correspondence between orbitals in two wavefunctions"
	!write(*,*) "7 Parse output of (hyper)polarizability task of Gaussian"
	!write(*,*) "8 Study (hyper)polarizability by sum-over-states (SOS) method"
	write(*,*) "9 Calculate average bond length and average coordinate number"
	write(*,*) "10 Output various kinds of integral between orbitals"
	write(*,*) "11 Calculate center, first/second moments and radius of gyration of a function"
	write(*,*) "12 Calculate energy index (EI) or bond polarity index (BPI)"
    write(*,*) "13 Evaluate orbital contributions to density difference or other grid data"
	write(*,*) "14 Domain analysis (Obtaining properties within isosurfaces of a function)"
	write(*,*) "15 Calculate electron correlation index (PCCP, 18, 24015)"
	write(*,*) "16 Generate natural orbitals based on the density matrix in .fch/.fchk file"
    write(*,*) "17 Calculate Coulomb and exchange integrals between two orbitals"
    write(*,*) "18 Calculate bond length/order alternation (BLA/BOA)"
    write(*,*) "19 Calculate spatial delocalization index (SDI) for orbitals or a function"
    write(*,*) "20 Bond order density (BOD) and natural adaptive orbital (NAdO) analyses"
    write(*,*) "21 Perform Lowdin orthogonalization between occupied orbitals"
	read(*,*) isel
	if (isel==0) then
		return
	else if (isel==1) then
		call CVB_index
	else if (isel==2) then
		call atmbonddip
	else if (isel==3) then
		call genmultiorbcube
	else if (isel==4) then !hidden
		call ICSS
	else if (isel==5) then
		call plotraddis
	else if (isel==6) then
		call orbcorres
	else if (isel==7) then !Hidden
		call parseGauPolar
	else if (isel==8) then !Hidden
		call SOS
	else if (isel==9) then
		call atmavgdist
	else if (isel==10) then
		call outorbint
	else if (isel==11) then
		call funcmoment
	else if (isel==12) then
		call calcEIBPI
    else if (isel==13) then
        call orbfitEDD
	else if (isel==14) then
		call domainana
	else if (isel==15) then
		call elecorridx
	else if (isel==16) then
		call fch_gennatorb
    else if (isel==17) then
        call orb_coulexcint
    else if (isel==18) then
        call BLABOA
    else if (isel==19) then
		call SDI
    else if (isel==20) then
        call BOD
    else if (isel==21) then
        call occorb_Lowdinorth_wrapper
	end if
end do
end subroutine




!!----- Calculate atomic and bond dipole moments in Hilbert space
!For derivation, see Ideas of Quantum Chemistry, p634
subroutine atmbonddip
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 xdipmat(nbasis,nbasis),ydipmat(nbasis,nbasis),zdipmat(nbasis,nbasis),Ptottmp(nbasis,nbasis)
character c80tmp*80

if (.not.allocated(CObasa)) then
	write(*,"(a)") " Error: No basis function information is provided in your input file! See Section 2.5 of Multiwfn manual for detail"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
if (.not.allocated(Dbas)) then
    call genDbas_curr
	write(*,*) "Generating electric dipole moment integral matrix..."
end if
xdipmat=Dbas(1,:,:)
ydipmat=Dbas(2,:,:)
zdipmat=Dbas(3,:,:)

!Calculate total dipole moment
xnucdip=sum(a(:)%charge*a(:)%x)
ynucdip=sum(a(:)%charge*a(:)%y)
znucdip=sum(a(:)%charge*a(:)%z)
write(*,"(' Molecular nuclear dipole moment (a.u.):')")
write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xnucdip,ynucdip,znucdip,dsqrt(xnucdip**2+ynucdip**2+znucdip**2)
xeledip=sum(Ptot*xdipmat)
yeledip=sum(Ptot*ydipmat)
zeledip=sum(Ptot*zdipmat)
write(*,"(' Molecular electron dipole moment (a.u.):')")
write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xeledip,yeledip,zeledip,dsqrt(xeledip**2+yeledip**2+zeledip**2)
xmoldip=xnucdip+xeledip
ymoldip=ynucdip+yeledip
zmoldip=znucdip+zeledip
write(*,"(' Molecular dipole moment (a.u.):')")
write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xmoldip,ymoldip,zmoldip,dsqrt(xmoldip**2+ymoldip**2+zmoldip**2)

do while(.true.)
	write(*,*)
	write(*,*) "         ----- Atomic and bond dipole moments in Hilbert space -----"
	write(*,*) "0 Return"
	write(*,*) "1 Output atomic dipole moment of specific atom"
	write(*,*) "2 Output bond dipole moment of specific atomic pair"
	write(*,*) "3 Output atomic overall dipole moment of specific atom (Mulliken partition)"
	write(*,*) "10 Export entire dipole moment matrix"
	read(*,*) isel
	if (isel==0) then
		exit
	else if (isel==1) then
		do while(.true.)
			write(*,*) "Input the atom index, e.g. 5"
			write(*,*) "Note: Input 0 can return, input -1 can output result for all atoms"
			read(*,*) isel2
			if (isel2==0) then
				exit
			else if (isel2==-1) then
				iatmstart=1
				iatmend=ncenter
			else
				if (isel2>ncenter) then
					write(*,*) "Atom index exceeded valid range! Input again"
					cycle
				end if
				iatmstart=isel2
				iatmend=isel2
			end if
			do iatm=iatmstart,iatmend
				if (basstart(iatm)==0) cycle
				write(*,"(' Result of atom',i8,' (',a2,')')") iatm,a(iatm)%name
				istart=basstart(iatm)
				iend=basend(iatm)
				atmelepop=sum(Ptot(istart:iend,istart:iend)*Sbas(istart:iend,istart:iend))
				xatmelediptot=sum(Ptot(istart:iend,istart:iend)*xdipmat(istart:iend,istart:iend))
				yatmelediptot=sum(Ptot(istart:iend,istart:iend)*ydipmat(istart:iend,istart:iend))
				zatmelediptot=sum(Ptot(istart:iend,istart:iend)*zdipmat(istart:iend,istart:iend))
				xatmeledip=xatmelediptot+atmelepop*a(iatm)%x
				yatmeledip=yatmelediptot+atmelepop*a(iatm)%y
				zatmeledip=zatmelediptot+atmelepop*a(iatm)%z
				xatmnucdip=a(iatm)%charge*a(iatm)%x
				yatmnucdip=a(iatm)%charge*a(iatm)%y
				zatmnucdip=a(iatm)%charge*a(iatm)%z
				xatmdiptot=xatmnucdip+xatmelediptot
				yatmdiptot=yatmnucdip+yatmelediptot
				zatmdiptot=zatmnucdip+zatmelediptot
				write(*,"(' Atomic local population number:',f12.6)") atmelepop
				write(*,"(' Atomic dipole moment (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmeledip,yatmeledip,zatmeledip,dsqrt(xatmeledip**2+yatmeledip**2+zatmeledip**2)
				write(*,"(' Contribution to system dipole moment due to nuclear charge (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmnucdip,yatmnucdip,zatmnucdip,dsqrt(xatmnucdip**2+yatmnucdip**2+zatmnucdip**2)
				write(*,"(' Contribution to system dipole moment due to electron (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmelediptot,yatmelediptot,zatmelediptot,dsqrt(xatmelediptot**2+yatmelediptot**2+zatmelediptot**2)
				write(*,"(' Contribution to system dipole moment (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmdiptot,yatmdiptot,zatmdiptot,dsqrt(xatmdiptot**2+yatmdiptot**2+zatmdiptot**2)
				write(*,*)
			end do
		end do
	else if (isel==2) then
		do while(.true.)
			write(*,*) "Input index of two atoms, e.g. 5,8"
			write(*,*) "Note: Input q can return. Input b can output result for all bonds"
			read(*,"(a)") c80tmp
			if (index(c80tmp,'q')/=0) then
				exit
			else if (index(c80tmp,'b')/=0) then
				write(*,*) "Notice that the bonds are determined according to distance criterion"
				write(*,*)
				bondcritval=1.15D0
			else
				read(c80tmp,*) iatomsel1,iatomsel2
				if (iatomsel1>ncenter.or.iatomsel2>ncenter) then
					write(*,*) "Atom index exceeded valid range! Input again"
					cycle
				end if
				if (iatomsel1>iatomsel2) then
					itmp=iatomsel2
					iatomsel2=iatomsel1
					iatomsel1=itmp
				end if
			end if
			do iatm=1,ncenter
				if (basstart(iatm)==0) cycle
				do jatm=iatm+1,ncenter
					if (basstart(jatm)==0) cycle
					bonddist=dsqrt((a(iatm)%x-a(jatm)%x)**2+(a(iatm)%y-a(jatm)%y)**2+(a(iatm)%z-a(jatm)%z)**2)
					if (index(c80tmp,'b')/=0) then
						if (bonddist>( covr(a(iatm)%index)+covr(a(jatm)%index) )*bondcritval) cycle
					else
						if (iatm/=iatomsel1.or.jatm/=iatomsel2) cycle
					end if
					xcen=(a(iatm)%x+a(jatm)%x)/2D0
					ycen=(a(iatm)%y+a(jatm)%y)/2D0
					zcen=(a(iatm)%z+a(jatm)%z)/2D0
					write(*,"(' Result between atom',i7,' (',a2,')  and atom',i7,' (',a2,'), distance:',f10.5,' Ang')") iatm,a(iatm)%name,jatm,a(jatm)%name,bonddist*b2a
					istart=basstart(iatm)
					iend=basend(iatm)
					jstart=basstart(jatm)
					jend=basend(jatm)
					bondpop=2*sum(Ptot(istart:iend,jstart:jend)*Sbas(istart:iend,jstart:jend)) !The matrix is symmetical, so multiplied by 2
					xbonddiptot=2*sum(Ptot(istart:iend,jstart:jend)*xdipmat(istart:iend,jstart:jend))
					ybonddiptot=2*sum(Ptot(istart:iend,jstart:jend)*ydipmat(istart:iend,jstart:jend))
					zbonddiptot=2*sum(Ptot(istart:iend,jstart:jend)*zdipmat(istart:iend,jstart:jend))
					xbonddip=xbonddiptot+bondpop*xcen
					ybonddip=ybonddiptot+bondpop*ycen
					zbonddip=zbonddiptot+bondpop*zcen
					write(*,"(' Bond population number (Overlap population):',f12.6)") bondpop
					write(*,"(' Bond dipole moment (a.u.):')")
					write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xbonddip,ybonddip,zbonddip,dsqrt(xbonddip**2+ybonddip**2+zbonddip**2)
					write(*,"(' Contribution to system dipole moment (a.u.):')")
					write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xbonddiptot,ybonddiptot,zbonddiptot,dsqrt(xbonddiptot**2+ybonddiptot**2+zbonddiptot**2)
					write(*,*)
				end do
			end do
		end do
	else if (isel==3) then
		do while(.true.)
			write(*,*) "Input the atom index, e.g. 5"
			write(*,*) "Note: Input 0 can return, input -1 can output result for all atoms"
			read(*,*) isel2
			if (isel2==0) then
				exit
			else if (isel2==-1) then
				iatmstart=1
				iatmend=ncenter
			else
				if (isel2>ncenter) then
					write(*,*) "Atom index exceeded valid range! Input again"
					cycle
				end if
				iatmstart=isel2
				iatmend=isel2
			end if
			do iatm=iatmstart,iatmend
				if (basstart(iatm)==0) cycle
				write(*,"(' Result of atom',i8,' (',a2,')')") iatm,a(iatm)%name
				atmelepop=0D0
				xatmelediptot=0D0
				yatmelediptot=0D0
				zatmelediptot=0D0
				istart=basstart(iatm)
				iend=basend(iatm)
				do jatm=1,ncenter
					if (basstart(jatm)==0) cycle
					jstart=basstart(jatm)
					jend=basend(jatm)
					atmelepop=atmelepop+sum(Ptot(istart:iend,jstart:jend)*Sbas(istart:iend,jstart:jend))
					xatmelediptot=xatmelediptot+sum(Ptot(istart:iend,jstart:jend)*xdipmat(istart:iend,jstart:jend))
					yatmelediptot=yatmelediptot+sum(Ptot(istart:iend,jstart:jend)*ydipmat(istart:iend,jstart:jend))
					zatmelediptot=zatmelediptot+sum(Ptot(istart:iend,jstart:jend)*zdipmat(istart:iend,jstart:jend))
				end do
				xatmeledip=xatmelediptot+atmelepop*a(iatm)%x
				yatmeledip=yatmelediptot+atmelepop*a(iatm)%y
				zatmeledip=zatmelediptot+atmelepop*a(iatm)%z
				xatmnucdip=a(iatm)%charge*a(iatm)%x
				yatmnucdip=a(iatm)%charge*a(iatm)%y
				zatmnucdip=a(iatm)%charge*a(iatm)%z
				xatmdiptot=xatmnucdip+xatmelediptot
				yatmdiptot=yatmnucdip+yatmelediptot
				zatmdiptot=zatmnucdip+zatmelediptot
				write(*,"(' Atomic Mulliken population number:',f12.6)") atmelepop
				write(*,"(' Atomic overall dipole moment (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmeledip,yatmeledip,zatmeledip,dsqrt(xatmeledip**2+yatmeledip**2+zatmeledip**2)
				write(*,"(' Contribution to system dipole moment due to nuclear charge (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmnucdip,yatmnucdip,zatmnucdip,dsqrt(xatmnucdip**2+yatmnucdip**2+zatmnucdip**2)
				write(*,"(' Contribution to system dipole moment due to electron (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmelediptot,yatmelediptot,zatmelediptot,dsqrt(xatmelediptot**2+yatmelediptot**2+zatmelediptot**2)
				write(*,"(' Contribution to system dipole moment (a.u.):')")
				write(*,"('  X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") xatmdiptot,yatmdiptot,zatmdiptot,dsqrt(xatmdiptot**2+yatmdiptot**2+zatmdiptot**2)
				write(*,*)
			end do
		end do
	else if (isel==10) then
		open(10,file="dipmatx.txt",status="replace")
		call showmatgau(Ptot*xdipmat,"",1,"f14.8",10)
		close(10)
		open(10,file="dipmaty.txt",status="replace")
		call showmatgau(Ptot*ydipmat,"",1,"f14.8",10)
		close(10)
		open(10,file="dipmatz.txt",status="replace")
		call showmatgau(Ptot*zdipmat,"",1,"f14.8",10)
		close(10)
		write(*,"(a)") " X, Y and Z components of electron dipole moment matrix have been outputted to dipmatx, dipmaty and dipmatz.txt in current folder, respectively"
	end if
end do
end subroutine


!!------------ Generate cube file for multiple orbitals
subroutine genmultiorbcube
use defvar
use util
implicit real*8 (a-h,o-z)
integer orbsellist(nmo)
integer tmparr(nmo+1)
character c1000tmp*1000,cubname*20
real*8,allocatable :: orbcubmat(:,:,:,:)
write(*,"(a)") " Input orbital index. e.g. 1,3-6,8,10-11 denotes 1,3,4,5,6,8,10,11"
call orblabsel_prompt
write(*,*) "Input q can return"
read(*,"(a)") c1000tmp

if (index(c1000tmp,'q')/=0) return
if (index(c1000tmp,'h')==0.and.index(c1000tmp,'l')==0) then !Use inputted a series of orbital indices
    call str2arr(c1000tmp,norbsel,orbsellist)
    if ( any(orbsellist(1:norbsel)<1) .or. any(orbsellist(1:norbsel)>nmo) ) then
	    write(*,*) "Error: The orbitals you selected exceeded valid range!"
	    return
    end if
    call setgrid(0,itmp)
    if (allocated(cubmat)) deallocate(cubmat)
    allocate(cubmat(nx,ny,nz))
    write(*,*)
    write(*,*) "1 Output the grid data of these orbitals as separate cube files"
    write(*,*) "2 Output the grid data of these orbitals as a single cube file"
    read(*,*) ioutmode

    if (ioutmode==1) then
	    do iorbidx=1,norbsel
		    iorb=orbsellist(iorbidx)
		    write(cubname,"('orb',i6.6,'.cub')") iorb
		    write(*,"(' Calculating and exporting orbital',i6)") iorb
		    call savecubmat(4,1,iorb)
		    open(10,file=cubname,status="replace")
		    call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		    close(10)
		    write(*,"(' Orbital',i7,' has been exported to ',a,' in current folder',/)") iorb,trim(cubname)
	    end do
	
    else if (ioutmode==2) then
	    allocate(orbcubmat(nx,ny,nz,norbsel))
	    do iorbidx=1,norbsel
		    iorb=orbsellist(iorbidx)
		    write(*,"(a,i6,a)") " Calculating grid data for orbital",iorb,"..."
		    call savecubmat(4,1,iorb)
		    orbcubmat(:,:,:,iorbidx)=cubmat(:,:,:)
	    end do
	    where (abs(orbcubmat)<=1D-99) orbcubmat=0D0 !Diminish too small value, otherwise the symbol "E" cannot be shown by 1PE13.5 format e.g. 9.39376-116, 
	    write(*,*)
	    write(*,*) "Exporting cube file, please wait..."
	    open(10,file="orbital.cub",status="replace")
	    write(10,"(' Generated by Multiwfn')")
	    write(10,"(' Totally ',i12,' grid points')") nx*ny*nz
	    write(10,"(i5,3f12.6)") -ncenter,orgx,orgy,orgz
	    write(10,"(i5,3f12.6)") nx,dx,0.0,0.0
	    write(10,"(i5,3f12.6)") ny,0.0,dy,0.0
	    write(10,"(i5,3f12.6)") nz,0.0,0.0,dz
	    do i=1,ncenter
		    write(10,"(i5,4f12.6)") a(i)%index,a(i)%charge,a(i)%x,a(i)%y,a(i)%z
	    end do
	    tmparr(1)=norbsel
	    tmparr(2:norbsel+1)=orbsellist(1:norbsel)
	    write(10,"(10i5)") tmparr(1:norbsel+1)
	    do ix=1,nx
		    do iy=1,ny
			    write(10,"(6(1PE13.5))",advance="no") ((orbcubmat(ix,iy,iz,iorbidx),iorbidx=1,norbsel),iz=1,nz)
			    write(10,*)
		    end do
	    end do
	    close(10)
	    write(*,*) "The grid data of the orbitals have been stored to orbital.cub in current folder"
	    deallocate(orbcubmat)
    end if
    deallocate(cubmat)
    
else
    if (index(c1000tmp,',')/=0) then
        write(*,*) "Error: In this selection mode, you can select only one orbital each time"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
    call orblabsel(c1000tmp,iorb)
    if (iorb==0) then
	    write(*,*) "Error: The orbital label you inputted is wrong! Please double check"
	    return
    end if
    call setgrid(0,itmp)
    if (allocated(cubmat)) deallocate(cubmat)
    allocate(cubmat(nx,ny,nz))
    cubname=trim(c1000tmp)//".cub"
	write(*,*) "Calculating and exporting orbital..."
	call savecubmat(4,1,iorb)
	open(10,file=cubname,status="replace")
	call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
	close(10)
	write(*,"(' Orbital',i7,' has been exported to ',a,' in current folder')") iorb,trim(cubname)
    deallocate(cubmat)
end if
end subroutine




!!----- Plot radial distribution function for a real space function
subroutine plotraddis
use defvar
use functions
use GUI
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:),radval(:),radpos(:),intradval(:),sphavgval(:)
ifunc=1
cenx=0D0
ceny=0D0
cenz=0D0
rlow=0D0
rhigh=5D0/b2a !5 Angstrom
nsphpt=2030
nradpt=500
do while(.true.)
	write(*,*)
	write(*,*) "  ====== Plot radial distribution function for a real space function ======"
	write(*,*) "-1 Exit"
	write(*,*) "0 Calculate radial distribution function and its integration curve"
	write(*,"(a,i4)") " 1 Select real space function, current:",ifunc
	write(*,"(a,3f10.4,' Ang')") " 2 Set sphere center, current",cenx*b2a,ceny*b2a,cenz*b2a
	write(*,"(a,2f9.4,' Ang')") " 3 Set lower and upper limit of radial distance, current:",rlow*b2a,rhigh*b2a
	write(*,"(a,i6)") " 4 Set the number of integration point in each shell, current:",nsphpt
	write(*,"(a,i6)") " 5 Set the number of radial points, current:",nradpt
	read(*,*) isel
	if (isel==-1) then
		return
	else if (isel==1) then
		call selfunc_interface(1,ifunc)
	else if (isel==2) then
		write(*,*) "Input sphere center (Angstrom), e.g. 0.0,1.2,-0.4"
		read(*,*) cenx,ceny,cenz
		cenx=cenx/b2a !Convert to Bohr
		ceny=ceny/b2a
		cenz=cenz/b2a
	else if (isel==3) then
		write(*,*) "Input lower and upper limit (Angstrom), e.g. 0.0,8.0"
		read(*,*) rlow,rhigh
		rlow=rlow/b2a !Convert to Bohr
		rhigh=rhigh/b2a
	else if (isel==4) then
		write(*,"(a)") " Input the number of integration point in each shell, the value must be one of &
		&110/170/230/266/302/434/590/770/974/1454/1730/2030/2354/2702/3074/3470/3890/4334/4802/5294/5810"
		read(*,*) nsphpt
	else if (isel==5) then
		write(*,*) "Input the number of radial points, e.g. 800"
		read(*,*) nradpt
	else if (isel==0) then
		allocate(potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt))
		allocate(radval(nradpt),radpos(nradpt),intradval(nradpt),sphavgval(nradpt)) !radval records RDF, radpot records r position, intradval records integration of RDF
		call Lebedevgen(nsphpt,potx,poty,potz,potw)
		radval=0D0
		radstp=(rhigh-rlow)/(nradpt-1)
		ifinish=0
		iprogstp=20
		iprogcrit=iprogstp
		write(*,*) "Calculating..."
		!$OMP PARALLEL DO SHARED(radval,radpos,ifinish,iprogcrit) PRIVATE(irad,radnow,isph,rnowx,rnowy,rnowz,tmpval) schedule(dynamic) NUM_THREADS(nthreads)
		do irad=1,nradpt
			radnow=rlow+(irad-1)*radstp
			radpos(irad)=radnow
			tmpval=0
			do isph=1,nsphpt
				rnowx=potx(isph)*radnow+cenx
				rnowy=poty(isph)*radnow+ceny
				rnowz=potz(isph)*radnow+cenz
				tmpval=tmpval+calcfuncall(ifunc,rnowx,rnowy,rnowz)*potw(isph)
			end do
			radval(irad)=4*pi*tmpval*radnow**2 !Multiplied by 4*pi is because the Lebedev integration routine produces unity rather than 4*pi
			sphavgval(irad)=tmpval !Spherically average function
	        ifinish=ifinish+1
	        if (ifinish==iprogcrit) then
				call showprog(ifinish,nradpt)
				iprogcrit=iprogcrit+iprogstp
			end if
		end do
		!$OMP END PARALLEL DO
		
		!Calculate integration of RDF
		intradval(1)=0D0
		do irad=2,nradpt
			intradval(irad)=intradval(irad-1)+radval(irad-1)*radstp
		end do
		write(*,"(a,f22.10)") " Integrating the RDF in the specified range is",intradval(nradpt)
		valrange=maxval(radval)-minval(radval)
		valrangeint=maxval(intradval)-minval(intradval)
		ilenunit1D=2
		do while(.true.)
			write(*,*)
			if (ilenunit1D==1) write(*,*) "-1 Switch the length unit for plotting, current: Bohr"
			if (ilenunit1D==2) write(*,*) "-1 Switch the length unit for plotting, current: Angstrom"
			write(*,*) "0 Return"
			write(*,*) "1 Plot the radial distribution function (RDF)"
			write(*,*) "2 Plot integration curve of the RDF"
			write(*,*) "3 Save the radial distribution function map in current folder"
			write(*,*) "4 Save integration curve of the RDF map in current folder"
			write(*,*) "5 Export the radial distribution function to RDF.txt in current folder"
			write(*,*) "6 Export integration curve of the RDF to intRDF.txt in current folder"
			write(*,*) "7 Export the spherically averaged function"
			read(*,*) isel
			if (isel==-1) then
				if (ilenunit1D==1) then
					ilenunit1D=2
				else if (ilenunit1D==2) then
					ilenunit1D=1
				end if
			else if (isel==0) then
				deallocate(potx,poty,potz,potw,radval,radpos,intradval,sphavgval)
				exit
			else if (isel==1.or.isel==3) then
				ylow=minval(radval)-0.1D0*valrange
				yhigh=maxval(radval)+0.1D0*valrange
				if (isel==1) then
					call drawcurve(radpos,radval,nradpt,rlow,rhigh,(rhigh-rlow)/10,ylow,yhigh,(yhigh-ylow)/10,"show")
				else
					call drawcurve(radpos,radval,nradpt,rlow,rhigh,(rhigh-rlow)/10,ylow,yhigh,(yhigh-ylow)/10,"save")
					write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""dislin"" prefix in current directory"
				end if
			else if (isel==2.or.isel==4) then
				ylow=minval(intradval)-0.1D0*valrangeint
				yhigh=maxval(intradval)+0.1D0*valrangeint
				if (isel==2) then
					call drawcurve(radpos,intradval,nradpt,rlow,rhigh,(rhigh-rlow)/10,ylow,yhigh,(yhigh-ylow)/10,"show")
				else if (isel==4) then
					call drawcurve(radpos,intradval,nradpt,rlow,rhigh,(rhigh-rlow)/10,ylow,yhigh,(yhigh-ylow)/10,"save")
					write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""dislin"" prefix in current directory"
				end if
			else if (isel==5) then
				open(10,file="RDF.txt",status="replace")
				do irad=1,nradpt
					write(10,"(i7,f12.4,f22.10)") irad,radpos(irad)*b2a,radval(irad)
				end do
				close(10)
				write(*,*) "The result has been output to RDF.txt in current folder"
				write(*,*) "The second column is radial distance (Angstrom), the third column is value"
			else if (isel==6) then
				open(10,file="intRDF.txt",status="replace")
				do irad=1,nradpt
					write(10,"(i7,f12.4,f22.10)") irad,radpos(irad)*b2a,intradval(irad)
				end do
				close(10)
				write(*,*) "The result has been output to intRDF.txt in current folder"
				write(*,*) "The second column is radial distance (Angstrom), the third column is value"
			else if (isel==7) then
				open(10,file="sphavgval.txt",status="replace")
				do irad=1,nradpt
					write(10,"(i7,f12.6,f18.7)") irad,radpos(irad),sphavgval(irad)
				end do
				close(10)
				write(*,*) "The result has been output to sphavgval.txt in current folder"
				write(*,*) "The second column is radial distance (Bohr), the third column is value"
			end if
		end do
	
	end if
end do
end subroutine



!!--------- Analyze correspondence between orbitals in two wavefunctions
subroutine orbcorres
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: convmat(:,:) !(i,j) is the coefficient of j MO of the second wavefunction in i MO of current wavefunction
real*8,allocatable :: Snormmat(:,:) !(i,j) is the overlap integral between norm of i MO of current wavefunction and norm of j MO of the second wavefunction
real*8,allocatable :: Snorm2mat(:,:) !(i,j) is the overlap integral between square of i MO of current wavefunction and norm of j MO of the second wavefunction
real*8,allocatable :: MOvalgrd(:,:),MOvalgrd2(:,:) !MOvalgrd(j,n),MOvalgrd2(j,n) means the the value of the nth MO of the first/second wavefunction at the jth grid
real*8,allocatable :: comparr(:),beckeweigrid(:)
integer,allocatable :: comparridx(:)
character filename2*200,c80tmp*80
type(content),allocatable :: gridatm(:),gridatmorg(:)
if (iautointgrid==1) then
	nsphpotold=sphpot
	nradpotold=radpot
    radcutold=radcut
	sphpot=302
	radpot=30
    radcut=15
end if
allocate(gridatm(radpot*sphpot),gridatmorg(radpot*sphpot),beckeweigrid(radpot*sphpot))

do isep=nmo,1,-1
	if (MOtype(isep)==1) exit
end do
if (wfntype==1.or.wfntype==4) write(*,"(' Note: The orbitals from',i6,' to',i6,' are alpha; from',i6,' to',i6,' are beta')") 1,isep,isep+1,nmo
write(*,"(a)") " Input lower and upper limits of the orbitals to be considered in the present wavefunction, e.g. 2,9"
write(*,*) "If press ENTER button directly, all orbitals will be taken into account"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
	istart1=1
	iend1=nmo
else
	read(c80tmp,*) istart1,iend1
end if

write(*,*)
write(*,*) "Input path of the second wavefunction, e.g. C:\ltwd.fch"
do while(.true.)
	read(*,"(a)") filename2
	inquire(file=filename2,exist=alive)
	if (alive) exit
	write(*,*) "Cannot find the file, input again"
end do
call dealloall(0)
call readinfile(filename2,1) !Get some knowledge about the second wavefunction
nmo2=nmo !The number of MOs in the wfn2
iwfntype2=wfntype
do isep=nmo2,1,-1
	if (MOtype(isep)==1) exit
end do
write(*,*)
if (iwfntype2==1.or.iwfntype2==4) write(*,"(' Note: The orbitals from',i6,' to',i6,' are alpha; from',i6,' to',i6,' are beta')") 1,isep,isep+1,nmo
write(*,"(a)") " Input lower and upper limits of the orbitals to be considered in the second wavefunction, e.g. 2,9"
write(*,*) "If press ENTER button directly, all orbitals will be taken into account"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
	istart2=1
	iend2=nmo2
else
	read(c80tmp,*) istart2,iend2
end if

call dealloall(0)
call readinfile(firstfilename,1)
allocate(MOvalgrd(radpot*sphpot,nmo),MOvalgrd2(radpot*sphpot,nmo2),convmat(nmo,nmo2),Snormmat(nmo,nmo2),Snorm2mat(nmo,nmo2))
convmat=0D0
Snormmat=0D0
Snorm2mat=0D0

write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
write(*,*) "Calculating, please wait..."
call gen1cintgrid(gridatmorg,iradcut)

call walltime(iwalltime1)

do iatm=1,ncenter
    call showprog(iatm,ncenter)
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	
	!Calculate value of all MOs of the first and second wavefunction at all grids
	call dealloall(0)
	call readinfile(filename2,1) !Load wfn2
	!$OMP parallel do shared(MOvalgrd2) private(ipt) num_threads(nthreads) schedule(dynamic)
	do ipt=1+iradcut*sphpot,radpot*sphpot
		call orbderv(1,istart2,iend2,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,MOvalgrd2(ipt,:))
	end do
	!$OMP end parallel do
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve to wfn1
	!$OMP parallel do shared(MOvalgrd) private(ipt) num_threads(nthreads) schedule(dynamic)
	do ipt=1+iradcut*sphpot,radpot*sphpot
		call orbderv(1,istart1,iend1,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,MOvalgrd(ipt,:))
	end do
	!$OMP end parallel do

	!Calculate Becke weight at all grids
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	
	!$OMP parallel do shared(convmat,Snormmat,Snorm2mat) private(imo,jmo,tmpval,tmpval2,tmpval2sqr,ipt) num_threads(nthreads) schedule(dynamic)
	do imo=istart1,iend1
		do jmo=istart2,iend2
			tmpval=0D0
			tmpval2=0D0
            tmpval2sqr=0
			do ipt=1+iradcut*sphpot,radpot*sphpot
				tmpval=tmpval+beckeweigrid(ipt)*gridatmorg(ipt)%value*MOvalgrd(ipt,imo)*MOvalgrd2(ipt,jmo)
				tmpval2=tmpval2+beckeweigrid(ipt)*gridatmorg(ipt)%value*abs(MOvalgrd(ipt,imo)*MOvalgrd2(ipt,jmo))
				tmpval2sqr=tmpval2sqr+beckeweigrid(ipt)*gridatmorg(ipt)%value*MOvalgrd(ipt,imo)**2*MOvalgrd2(ipt,jmo)**2
			end do
			convmat(imo,jmo)=convmat(imo,jmo)+tmpval
            Snormmat(imo,jmo)=Snormmat(imo,jmo)+tmpval2
            Snorm2mat(imo,jmo)=Snorm2mat(imo,jmo)+tmpval2sqr
		end do
	end do
	!$OMP end parallel do
end do
call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,'s',/)") iwalltime2-iwalltime1

! call showmatgau(convmat,"convmat",0,"f12.3")

devmax=0D0
idevmax=1
! write(*,*) "The sum of composition of each orbital of current wavefunction"
if ((wfntype==0.or.wfntype==2).and.(iwfntype2==0.or.iwfntype2==2)) then !Both of the two wavefunctions are R or RO types
	do imo=istart1,iend1 !Check normalization
		totcomp=sum(convmat(imo,:)**2)*100D0
		if (abs(totcomp-100D0)>devmax) then
			devmax=abs(totcomp-100D0)
			idevmax=imo
		end if
	end do
end if
write(*,"(' The maximum deviation to normalization condition is',f8.3,' % (Orbital',i6,')')") devmax,idevmax
write(*,"(a)") " Note: The first column below is the index of the orbitals in present wavefunction, the largest five contributions from &
&the orbitals in the second wavefunction are shown at right side. If the dominative index is inconsistent to the first column, the row will be marked by asterisk"
write(*,*)
allocate(comparr(nmo2),comparridx(nmo2))
do imo=istart1,iend1
	!Sort the composition array from small to large based on magnitude of coefficients
	comparr(:)=convmat(imo,:)
	do itmp=1,nmo2
		comparridx(itmp)=itmp
	end do
	do i=istart2,iend2
		do j=i+1,iend2
			if (abs(comparr(i))>abs(comparr(j))) then
				temp=comparr(i)
				comparr(i)=comparr(j)
				comparr(j)=temp
				itemp=comparridx(i)
				comparridx(i)=comparridx(j)
				comparridx(j)=itemp
			end if
		end do
	end do
	if (comparridx(iend2)/=imo) then
		write(*,"('*',i5,':  ')",advance="no") imo
	else
		write(*,"(' ',i5,':  ')",advance="no") imo
	end if
	do i=iend2,iend2-4,-1
		write(*,"(i5,'(',f6.2,'%)')",advance="no") comparridx(i),comparr(i)**2*100D0
	end do
	write(*,*)
end do

write(*,*)
do while(.true.)
	write(*,*) "Input the orbital index to print detail compositions and coefficients, e.g. 5"
	write(*,*) "  Input -1 can output all overlap integrals between the chosen orbitals"
	write(*,"(a)") "   Input -2 can output all overlap integrals of wavefunction norm between the chosen orbitals"
	write(*,"(a)") "   Input -3 can output all overlap integrals of square of wavefunction between the chosen orbitals"
	write(*,*) "  Input 0 can exit"
	read(*,*) imo
	if (imo==0) then
		exit
	else if (imo==-1) then
		open(10,file="convmat.txt",status="replace")
		do i=istart1,iend1
			do j=istart2,iend2
				write(10,"(2i7,f12.6)") i,j,convmat(i,j)
			end do
		end do
		close(10)
		write(*,"(a,/)") " The overlap integrals have been outputted to convmat.txt in current folder, &
        &the first and second columns correspond to the orbital indices in present and in the second wavefunctions, respectively"
	else if (imo==-2) then
		open(10,file="Snormmat.txt",status="replace")
		do i=istart1,iend1
			do j=istart2,iend2
				write(10,"(2i7,f12.6)") i,j,Snormmat(i,j)
			end do
		end do
		close(10)
		write(*,"(a,/)") " The overlap integrals of wavefunction norms have been outputted to Snormmat.txt in current folder, &
        &the first and second columns correspond to the orbital indices in present and in the second wavefunctions, respectively"
	else if (imo==-3) then
		open(10,file="Snorm2mat.txt",status="replace")
		do i=istart1,iend1
			do j=istart2,iend2
				write(10,"(2i7,f18.12)") i,j,Snorm2mat(i,j)
			end do
		end do
		close(10)
		write(*,"(a,/)") " The overlap integrals of square of wavefunctions have been outputted to Snorm2mat.txt in current folder, &
        &the first and second columns correspond to the orbital indices in present and in the second wavefunctions, respectively"
	else if (imo<istart1.or.imo>iend1) then
		write(*,"(a,i6,a,i6)") "Error: Exceeded valid range! The value should within",istart1," and",iend1
	else
		tmpval=0D0
		do jmo=istart2,iend2
			tmpval=tmpval+convmat(imo,jmo)**2*100D0
			write(*,"(i6,'   Contribution:',f10.3,' %    Coefficient:',f12.6)") jmo,convmat(imo,jmo)**2*100D0,convmat(imo,jmo)
		end do
		write(*,"(' Total:',f10.3,' %')") tmpval
		write(*,*)
	end if
end do

if (iautointgrid==1) then
	radpot=nradpotold
	sphpot=nsphpotold
    radcut=radcutold
end if
end subroutine





!!---------- Calculate average bond length between two elements and average coordinate number
subroutine atmavgdist
use defvar
use util
implicit real*8 (a-h,o-z)
character elesel1*2,elesel2*2,selectyn
if (all(a%name==a(1)%name)) then !Only contain one kind of atom
	elesel1=a(1)%name
	elesel2=a(1)%name
else
	write(*,*) "Input two elements for which their average bond length will be calculated"
	write(*,*) "For example, B,Al"
	read(*,*) elesel1,elesel2
	call lc2uc(elesel1(1:1))
	call uc2lc(elesel1(2:2))
	call lc2uc(elesel2(1:1))
	call uc2lc(elesel2(2:2))
end if
write(*,*) "Input distance cutoff in Angstrom, e.g. 3.2"
read(*,*) discrit
discrit=discrit/b2a

avgdist=0
iwithin=0
distmax=0
distmin=1000000
do iatm=1,ncenter
	do jatm=iatm+1,ncenter
		if ((a(iatm)%name==elesel1.and.a(jatm)%name==elesel2).or.(a(jatm)%name==elesel1.and.a(iatm)%name==elesel2)) then
			dist=atomdist(iatm,jatm,1)
			if (dist<=discrit) then
				iwithin=iwithin+1
				write(*,"(i6,'#   ',i6,'(',a')   --',i6,'(',a,')    Length:',f12.6,' Angstrom')") iwithin,iatm,a(iatm)%name,jatm,a(jatm)%name,dist*b2a
				avgdist=avgdist+dist
				if ((iatm==1.and.jatm==iatm+1).or.dist>distmax) then
					distmax=dist
					idistmax1=iatm
					idistmax2=jatm
				end if
				if ((iatm==1.and.jatm==iatm+1).or.dist<distmin) then
					distmin=dist
					idistmin1=iatm
					idistmin2=jatm
				end if
			end if
		end if
	end do
end do
if (iwithin>0) then
	avgdist=avgdist/iwithin
	write(*,"(' Average bond length between ',a,1x,a,' is',f12.6,' Angstrom')") elesel1,elesel2,avgdist*b2a
	write(*,"(' Minimum length is',f12.6,' Angstrom, between',i6,'(',a,') and',i6,'(',a,')')") distmin*b2a,idistmin1,elesel1,idistmin2,elesel2
	write(*,"(' Maximum length is',f12.6,' Angstrom, between',i6,'(',a,') and',i6,'(',a,')')") distmax*b2a,idistmax1,elesel1,idistmax2,elesel2
else
	write(*,*) "No bond that satisfied your criterion is found"
	return
end if
write(*,*)
write(*,*) "If also calculate average coordinate number? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	ncoordtot=0
	ntmp=0
	do iatm=1,ncenter
		if (a(iatm)%name/=elesel1) cycle
		ncoord=0
		do jatm=1,ncenter
			if (iatm==jatm.or.a(jatm)%name/=elesel2) cycle
			if (atomdist(iatm,jatm,1)<=discrit) ncoord=ncoord+1
		end do
		write(*,"(' The coordinate number of',i6,'(',a,') due to ',a,' - ',a,' bond:',i5)") iatm,elesel1,elesel1,elesel2,ncoord
		ncoordtot=ncoordtot+ncoord
		ntmp=ntmp+1
	end do
	write(*,"(/,' The average coordinate number of ',a,' due to ',a,' - ',a,' bond:',f10.5)") elesel1,elesel1,elesel2,dfloat(ncoordtot)/ntmp
end if
end subroutine





!!!------- Calculate electric/magnetic/velocity... integral between orbitals
subroutine outorbint
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: GTFint(:),GTFvecint(:,:)
real*8 vecint(3),vecinttmp(3),intval
write(*,*) "Output which kind of integral between orbitals?"
write(*,*) "1: Electric dipole moment  2: Magnetic dipole moment  3: Velocity"
write(*,*) "4: Kinetic energy   5: Overlap"
read(*,*) itype
if (wfntype==0.or.wfntype==1.or.wfntype==2) then
	write(*,*) "Output the integrals between which orbitals?"
	if (wfntype==0.or.wfntype==1.or.wfntype==2) then
		write(*,*) "1 Between all occupied orbitals"
		write(*,*) "2 Between all occupied and all unoccupied orbitals"
	end if
	write(*,*) "3 Between all orbitals"
	write(*,*) "4 Between the same orbitals"
	write(*,*) "5 Between specific range of orbitals"
	write(*,*) "6 Between two specific orbitals"
	read(*,*) irange
end if
if (irange==5) then
	write(*,*) "Input orbital range of the first index, e.g. 25,30"
	read(*,*) ibeg,iend
	write(*,*) "Input orbital range of the second index, e.g. 25,30"
	read(*,*) jbeg,jend
else if (irange==6) then
100	write(*,*) "Input index for the two orbitals, e.g. 144,340"
	write(*,*) "Input 0,0 can exit"
	read(*,*) iMOsel,jMOsel
	if (iMOsel==0.and.jMOsel==0) return
end if

call walltime(iwalltime1)

nsize=nprims*(nprims+1)/2
if (itype==1.or.itype==2.or.itype==3) then
	allocate(GTFvecint(3,nsize))
else
	allocate(GTFint(nsize))
end if
if (itype==1) then
	call genGTFDmat(GTFvecint,nsize)
else if (itype==2) then
	call genGTFMmat(GTFvecint,nsize) !Notice that this only generate lower triangular matrix, (i,j)=-(j,i)
else if (itype==3) then	
	call genGTFVelmat(GTFvecint,nsize)
else if (itype==4) then	
	call genGTFTmat(GTFint,nsize)
else if (itype==5) then	
	call genGTFSmat(GTFint,nsize)
end if

if (irange/=6) open(10,file="orbint.txt",status="replace")
do imo=1,nmo
	do jmo=1,nmo
		if (irange==1) then
			if (MOocc(imo)==0D0.or.MOocc(jmo)==0D0) cycle
		else if (irange==2) then
			if (MOocc(imo)==0D0.or.MOocc(jmo)/=0D0) cycle
		else if (irange==4) then
			if (imo/=jmo) cycle
		else if (irange==5) then
			if (imo<ibeg.or.imo>iend.or.jmo<jbeg.or.jmo>jend) cycle
		else if (irange==6) then
			if (imo/=iMOsel.or.jmo/=jMOsel) cycle
		end if
		
		if (itype==1.or.itype==2.or.itype==3) then !Vector integral
			vecint=0D0
			!$OMP PARALLEL SHARED(vecint) PRIVATE(iGTF,jGTF,ides,vecinttmp) NUM_THREADS(nthreads)
			vecinttmp=0D0
			!$OMP DO schedule(dynamic)
			do iGTF=1,nprims
				do jGTF=1,nprims
					if (iGTF>=jGTF) then
						ides=iGTF*(iGTF-1)/2+jGTF
					else
						ides=jGTF*(jGTF-1)/2+iGTF
					end if
					if ((itype==2.or.itype==3).and.iGTF>jGTF) then !Magnetic and velocity operators are Hermitean, so inverse sign
						vecinttmp=vecinttmp-CO(imo,iGTF)*CO(jmo,jGTF)*GTFvecint(:,ides)
					else
						vecinttmp=vecinttmp+CO(imo,iGTF)*CO(jmo,jGTF)*GTFvecint(:,ides)
					end if
				end do
			end do
			!$OMP END DO
			!$OMP CRITICAL
			vecint=vecint+vecinttmp
			!$OMP END CRITICAL
			!$OMP END PARALLEL
			valnorm=sqrt(sum(vecint(:)**2))
			if (irange==6) then
				write(*,"(' X, Y, Z:',3f18.10,' a.u.')") vecint(:)
				write(*,"(' Norm:',f18.10,' a.u.')") valnorm
			else
				write(10,"(2i8,4x,4f18.10)") imo,jmo,vecint(:),valnorm
			end if
		else !Scalar integral
			intval=0D0
			do iGTF=1,nprims
				do jGTF=1,nprims
					if (iGTF>=jGTF) then
						ides=iGTF*(iGTF-1)/2+jGTF
					else
						ides=jGTF*(jGTF-1)/2+iGTF
					end if
					intval=intval+CO(imo,iGTF)*CO(jmo,jGTF)*GTFint(ides)
				end do
			end do
			if (irange==6) then
				write(*,"(' Result:',f18.10)") intval
			else
				write(10,"(2i8,4x,f18.10)") imo,jmo,intval
			end if
		end if
	end do
	if (irange/=6) write(*,"(' Finished',i7,' /',i7)") imo,nmo
end do
if (irange==6) then
	if (allocated(GTFvecint)) deallocate(GTFvecint)
	if (allocated(GTFint)) deallocate(GTFint)
	write(*,*)
	goto 100
else
	close(10)
	call walltime(iwalltime2)
	write(*,"(' Done! Calculation took up wall clock time',i10,' s',/)") iwalltime2-iwalltime1
	write(*,*) "The integrals have been outputted to orbint.txt in current folder"
	if (itype==1.or.itype==2.or.itype==3) then
		write(*,"(a)") " The first and the second columns correspond to orbital indices, &
		&the next three columns correspond to the integral in X/Y/Z (a.u.), the final column is the norm"
	else
		write(*,"(a)") " The first and the second columns correspond to orbital indices. The last column corresponds to the integral value (a.u.)"
	end if
end if
end subroutine





!!----------- Calculate energy index (EI) or bond polarity index (BPI)
!!J. Phys. Chem., 94, 5602-5607 and J. Phys. Chem.,96, 157-164
subroutine calcEIBPI
use defvar
implicit real*8 (a-h,o-z)
if (wfntype==3.or.wfntype==4) then
	write(*,*) "Error: Multiconfiguration wavefunction has not been supported yet!"
	return
end if
call getninnerele(ninnerele,0)
write(*,"(' The number of inner electrons is assumed to be',i5,/)") ninnerele
do while(.true.)
	write(*,*) "Calculate EI index for which atom? e.g. 5"
	write(*,*) "Input 0 can exit"
	read(*,*) iatm
	if (iatm==0) exit
	val1=0D0
	val2=0D0
	do imo=ninnerele/2+1,nbasis
		if (MOocc(imo)==0D0) exit
		compos=0
		do ibas=basstart(iatm),basend(iatm)
			compos=compos+CObasa(ibas,imo)**2
			do jbas=1,nbasis
				if (jbas==ibas) cycle
				compos=compos+CObasa(ibas,imo)*CObasa(jbas,imo)*Sbas(ibas,jbas)
			end do
		end do
! 		write(*,"(i5,3f12.6)") imo,MOocc(imo),MOene(imo),compos
		val1=val1+MOocc(imo)*MOene(imo)*compos
		val2=val2+MOocc(imo)*compos
	end do
	if (wfntype==1) then !beta part
		do imo=nbasis+ninnerele/2+1,nmo
			if (MOocc(imo)==0D0) exit
			compos=0
			do ibas=basstart(iatm),basend(iatm)
				compos=compos+CObasb(ibas,imo-nbasis)**2
				do jbas=1,nbasis
					if (jbas==ibas) cycle
					compos=compos+CObasb(ibas,imo-nbasis)*CObasb(jbas,imo-nbasis)*Sbas(ibas,jbas)
				end do
			end do
! 			write(*,"(i5,3f12.6)") imo,MOocc(imo),MOene(imo),compos
			val1=val1+MOocc(imo)*MOene(imo)*compos
			val2=val2+MOocc(imo)*compos
		end do
	end if
	write(*,"(' The numerator:  ',f12.6,' a.u.')") val1
	write(*,"(' The denominator:',f12.6,' a.u.')") val2 !Corresponding to Mulliken occupation number of this atom in valence MOs
	write(*,"(' The EI index:   ',f12.6,' a.u.',/)") val1/val2
end do
end subroutine


!!------------ Domain analysis (Integrate real space function within isosurface of a real space function)
!I use the same data structure as basin analysis to illustrate definition of isosurfaces
subroutine domainana
use defvar
use GUI
use util
use functions
use basinintmod !Use its vec26 array
implicit real*8 (a-h,o-z)
integer :: ifunciso=13,ifPBCgrid=0
integer,allocatable :: mergelist(:),grididx(:,:,:),dogrid(:,:)
logical,allocatable :: boundgrid(:)
character :: defdomain*20="<0.5",c80tmp*80,c1000tmp*1000,c2000tmp*2000,cubname*200
integer,allocatable :: tmparr(:)

if (allocated(gridxyz)) deallocate(gridxyz)
if (allocated(domainsize)) deallocate(domainsize)
if (allocated(domaingrid)) deallocate(domaingrid)
if (ifPBC>0) ifPBCgrid=1

do while(.true.)
	write(*,*)
    write(*,*) "  ---------------------------- Domain analysis ----------------------------"
	if (allocated(cubmat)) write(*,*) "-1 Yield domains based on the grid data in memory"
	write(*,*) "0 Return"
	write(*,*) "1 Calculate grid data and yield domains"
	write(*,"(a,i5)") " 2 Select real space function to be calculated for option 1, current:",ifunciso
	write(*,"(a,a)") " 3 Set criterion for defining domain, current: ",trim(defdomain)
    if (ifPBCgrid==0) write(*,*) "4 Toggle considering periodicity during domain analysis, current: No"
    if (ifPBCgrid==1) write(*,*) "4 Toggle considering periodicity during domain analysis, current: Yes"
	read(*,*) isel
	if (isel==0) then
		return
	else if (isel==1.or.isel==-1) then
		exit
	else if (isel==2) then
		call selfunc_interface(1,ifunciso)
	else if (isel==3) then
		write(*,*) "Input the definition, e.g. <0.05"
		write(*,*) "Note: The first character must be < or >"
		read(*,"(a)") defdomain
	else if (isel==4) then
		if (ifPBCgrid==0) then
			ifPBCgrid=1
        else
			ifPBCgrid=0
        end if
	end if
end do

!Set grid and generate grid data
if (isel==1) then
	call setgridfixspc
	if (allocated(cubmat)) deallocate(cubmat)
	allocate(cubmat(nx,ny,nz))
	call savecubmat(ifunciso,0,iorbsel)
end if
call calc_dvol(dvol)

!Count the number of grids satisfying the criterion
read(defdomain(2:),*) valiso
if (defdomain(1:1)=='<') then
	ngrid=count(cubmat<valiso)
else if (defdomain(1:1)=='>') then
	ngrid=count(cubmat>valiso)
else
	write(*,*) "Error: Parsing of the domain definition failed!"
	return
end if
write(*,"(/,' The number of grids satisfying the criterion:',i10)") ngrid

!Clustering grids that satisfied criterion to domain
!The idea is very clever:
!For grids that meet isovalue criterion, I assign each grid with a different index, and perform iterations, in each iteration all of these grids &
!are looped over, index of each grid is set to that of one of the nearest 26 grids if its index is larger than current grid. Finally, &
!grids in each domain will have identical index
write(*,*) "Clustering domains..."
call walltime(iwalltime1)
!dogrid: ix/iy/iz index of grids satisfying condition
!gridxyz: XYZ coordinate of grids satisfying condition
!grididx: Firstly, records initial index of grids satisfying condition; then update so that grids in each domain share identical index; at final stage, it records domain index
!boundgrid: If the grid is boundary grid (grids around isosurfaces, or at cell boundary)
allocate(dogrid(3,ngrid),grididx(nx,ny,nz),gridxyz(3,ngrid),boundgrid(ngrid))
boundgrid=.false.

!Initialize grid indices. Make each grid statisfying the condition have unique index, and build gridxyz
write(*,*) "Initializing grid indices..."
grididx=-1 !Irrelevant grids have very small value
idx=0
do iz=1,nz
	do iy=1,ny
		do ix=1,nx
			if ((defdomain(1:1)=='<'.and.cubmat(ix,iy,iz)<valiso).or.(defdomain(1:1)=='>'.and.cubmat(ix,iy,iz)>valiso)) then
				idx=idx+1
				dogrid(1,idx)=ix
				dogrid(2,idx)=iy
				dogrid(3,idx)=iz
				grididx(ix,iy,iz)=idx
                call getgridxyz(ix,iy,iz,gridxyz(1,idx),gridxyz(2,idx),gridxyz(3,idx))
			end if
		end do
	end do
end do

call setupmovevec

!Determine if grid is boundary grid
write(*,*) "Determining which grids are boundary grids..."
do igrid=1,ngrid !Loop each grid statisfying condition
	ix=dogrid(1,igrid)
	iy=dogrid(2,igrid)
	iz=dogrid(3,igrid)
	if (ix==1.or.ix==nx.or.iy==1.or.iy==ny.or.iz==1.or.iz==nz) then !If current grid is at box boundary, it will be regarded as boundary grid
		boundgrid(igrid)=.true.
		cycle
	end if
	do imove=1,26 !Check each neighbouring grid. If function value of any neighbouring grid does not meet condition, present grid should be boundary grid
		ixtmp=ix+vec26x(imove);iytmp=iy+vec26y(imove);iztmp=iz+vec26z(imove) !ix,iy,iz index of neighbouring grid
		if (ixtmp<1.or.ixtmp>nx.or.iytmp<1.or.iytmp>ny.or.iztmp<1.or.iztmp>nz) cycle !Skip grids out of box boundary
		valtmp=cubmat(ixtmp,iytmp,iztmp)
		if ((defdomain(1:1)=='<'.and.valtmp>valiso).or.(defdomain(1:1)=='>'.and.valtmp<valiso)) then
			boundgrid(igrid)=.true.
			exit
		end if
	end do
end do

!Iteration to make indices of grids in each domain are identical
write(*,*) "Iteratively updating indices of grids in each domain..."
do while(.true.)
	iupdate=0 !If there are grids updated index in this cycle
	do itmp=1,ngrid !Cycle all grids that satisify condition
		ix=dogrid(1,itmp)
		iy=dogrid(2,itmp)
		iz=dogrid(3,itmp)
		do imove=1,26 !Cycle neighbouring grids
			ixtmp=ix+vec26x(imove)
			iytmp=iy+vec26y(imove)
			iztmp=iz+vec26z(imove)
            if (ifPBCgrid==0) then
				if (ixtmp<1.or.ixtmp>nx.or.iytmp<1.or.iytmp>ny.or.iztmp<1.or.iztmp>nz) cycle !Skip grids at box boundary
            else if (ifPBCgrid==1) then
				call PBCgrididx(ixtmp,iytmp,iztmp) !If neighbouring grid is outside cell, wrap its index into the cell
            end if
			if (grididx(ix,iy,iz)<grididx(ixtmp,iytmp,iztmp)) then
				grididx(ix,iy,iz)=grididx(ixtmp,iytmp,iztmp)
				iupdate=1
				exit
			end if
		end do
	end do
	if (iupdate==0) exit !This iteration does not update any index of the grids to be considered, so finished
end do

!After below step, grididx will record domain index of each grid that satisfies condition
write(*,*) "Building the array recording domain index of each grid..."
ndone=0
ndomain=0
do while(.true.)
	ndomain=ndomain+1
	mintmp=2147483647 !Largest value can be recorded by integer*4
	do itmp=1,ngrid
		idxtmp=grididx(dogrid(1,itmp),dogrid(2,itmp),dogrid(3,itmp))
		if (idxtmp<ndomain) cycle
		if (idxtmp<mintmp) mintmp=idxtmp
	end do
	do itmp=1,ngrid
		idxtmp=grididx(dogrid(1,itmp),dogrid(2,itmp),dogrid(3,itmp))
		if (idxtmp==mintmp) then
			grididx(dogrid(1,itmp),dogrid(2,itmp),dogrid(3,itmp))=ndomain
			ndone=ndone+1
		end if
	end do
	if (ndone==ngrid) exit
end do

!Generate domainsize and domaingrid (grid index that contained in each domain)
write(*,"(a)") " Counting number of grids in each domain and building list of grids in each domain..."
allocate(domainsize(ndomain),domaingrid(ngrid,ndomain))
do idom=1,ndomain
	j=0
	do itmp=1,ngrid
		if (grididx(dogrid(1,itmp),dogrid(2,itmp),dogrid(3,itmp))==idom) then
			j=j+1
			domaingrid(j,idom)=itmp
		end if
	end do
	domainsize(idom)=j
end do

write(*,*) "Clustering domains finished!"
call walltime(iwalltime2)
write(*,"(' Clustering took up wall clock time',i10,' s')") iwalltime2-iwalltime1

do idom=1,ndomain
	write(*,"(' Domain:',i7,'    Grids:',i9,'    Volume:',f11.4,' Angstrom^3')") idom,domainsize(idom),domainsize(idom)*dvol*b2a**3
end do

do while(.true.)
	write(*,*)
    write(*,*) " ----------------- Post-processing menu of domain analysis -----------------"
	write(*,*) "-1 Merge specific domains"
	write(*,*) "0 Exit"
	write(*,*) "1 Perform integration for a domain"
	write(*,*) "2 Perform integration for all domains"
	write(*,"(a)") " 2b Perform integration for subregion of some domains according to range of sign(lambda2)*rho"
	write(*,*) "3 Visualize domains"
    write(*,*) "4 Sort indices of domains according to their domain volumes"
	write(*,*) "5 Calculate q_bind index for a domain"
	write(*,*) "10 Export a domain as domain.cub file in current folder"
	write(*,*) "11 Export boundary grids of a domain to domain.pdb file in current folder"
    write(*,"(a)") " 12 Export X,Y,Z coordinates and value of all grids in each domain to domain.txt in current folder"
    read(*,"(a)") c80tmp
    if (index(c80tmp,'b')/=0) then
		isel2=-2
    else
		read(c80tmp,*) isel2
    end if
    
	if (isel2==0) then
		return
	else if (isel2==-1) then
		if (ndomain<2) then
			write(*,*) "Error: At least two domains must be presented!"
			cycle
		end if
		write(*,*) "Input indices of the domains you want to merge, e.g. 4,5,8-10"
		read(*,"(a)") c1000tmp
		call str2arr(c1000tmp,nmerge) !Find how many terms
		allocate(mergelist(nmerge))
		call str2arr(c1000tmp,nmerge,mergelist)
		call sort(mergelist)
		idom=mergelist(1)
		do jdx=nmerge,2,-1 !Gradually merge the last domain (jdom) in the list into the first domain (idom)
			jdom=mergelist(jdx)
			nsizei=domainsize(idom)
			nsizej=domainsize(jdom)
			domainsize(idom)=nsizei+nsizej
			domaingrid(nsizei+1:nsizei+nsizej,idom)=domaingrid(1:nsizej,jdom)
			ndomain=ndomain-1
			domainsize(jdom:ndomain)=domainsize(jdom+1:ndomain+1) !Move all after jdom (disappeared) forward
			domaingrid(:,jdom:ndomain)=domaingrid(:,jdom+1:ndomain+1)
			where (grididx(:,:,:)==jdom) grididx=idom
		end do
		deallocate(mergelist)
		idrawdomainidx=0 !Do not draw domain in the visualizer
		write(*,"(a,i6)") " Done! The domains you selected have been merged as domain",idom
		write(*,*) "Size of current domains:"
		do idom=1,ndomain
			write(*,"(' Domain:',i6,'     Grids:',i8)") idom,domainsize(idom)
		end do
        
	else if (isel2==1) then !Perform integration for a domain
		write(*,*) "Input the index of the domain to be integrated, e.g. 3"
		read(*,*) intdom
		if (intdom<1.or.intdom>ndomain) then
			write(*,"(a)") " Error: The index of the domain to be integrated is incorrect"
			cycle
		end if
        write(*,*) "Which data will be integrated in the selected domain?"
        write(*,*) "1 The grid data in memory"
        write(*,*) "2 Choose a real space function, which will be directly calculated"
        write(*,*) "3 Grid data loaded from a .cub file"
        read(*,*) inttype
        if (inttype==2) then
			write(*,*) "Select the real space function to be integrated, e.g. 1"
			call selfunc_interface(1,ifuncint)
        else if (inttype==3) then
			if (allocated(cubmattmp)) deallocate(cubmattmp)
            write(*,*) "Input path of the cube file, e.g. D:\ltwd\rho.cub"
            write(*,"(a)") " Note: Distribution of the grids in this file must be exactly identical to that of grid data in memory"
            do while(.true.)
				read(*,"(a)") cubname
				inquire(file=cubname,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again!"
			end do
            call readcubetmp(cubname,1,inconsis)
            if (inconsis==1) then
				write(*,*) "Error: Distribution of the grids in this file is not exactly identical to that of grid data in memory"
                write(*,*) "Press ENTER button to return"
                read(*,*)
                cycle
            end if
        end if
		valint=0
		volint=0
		valmin=1D100
		valmax=-1D100
		xmin=1D100;xmax=-1D100;ymin=1D100;ymax=-1D100;zmin=1D100;zmax=-1D100
		do igrd=1,domainsize(intdom)
			idx=domaingrid(igrd,intdom)
			xnow=gridxyz(1,idx)
			ynow=gridxyz(2,idx)
			znow=gridxyz(3,idx)
            if (inttype==1) then
				tmpval=cubmat(dogrid(1,idx),dogrid(2,idx),dogrid(3,idx))
            else if (inttype==2) then
				tmpval=calcfuncall(ifuncint,xnow,ynow,znow)
            else if (inttype==3) then
				tmpval=cubmattmp(dogrid(1,idx),dogrid(2,idx),dogrid(3,idx))
            end if
			valint=valint+tmpval
			volint=volint+1
			if (tmpval<valmin) valmin=tmpval
			if (tmpval>valmax) valmax=tmpval
			if (xnow<xmin) xmin=xnow
			if (xnow>xmax) xmax=xnow
			if (ynow<ymin) ymin=ynow
			if (ynow>ymax) ymax=ynow
			if (znow<zmin) zmin=znow
			if (znow>zmax) zmax=znow
		end do
		avgval=valint/domainsize(intdom)
		valint=valint*dvol
		volint=volint*dvol
		write(*,"(' Integration result:',E20.10,' a.u.')") valint
		write(*,"(' Volume:',f12.6,' Bohr^3  (',f12.6,' Angstrom^3 )')") volint,volint*b2a**3
		write(*,"(' Average:',E20.10)") avgval
		write(*,"(' Maximum:',E20.10,'   Minimum:',E20.10)") valmax,valmin
		write(*,"(/,' Position statistics for coordinates of domain points (Angstrom):')")
		write(*,"(' X minimum:',f10.4,'  X maximum:',f10.4,'  Span:',f10.4)") xmin*b2a,xmax*b2a,(xmax-xmin)*b2a
		write(*,"(' Y minimum:',f10.4,'  Y maximum:',f10.4,'  Span:',f10.4)") ymin*b2a,ymax*b2a,(ymax-ymin)*b2a
		write(*,"(' Z minimum:',f10.4,'  Z maximum:',f10.4,'  Span:',f10.4)") zmin*b2a,zmax*b2a,(zmax-zmin)*b2a
        
	else if (isel2==2) then !Perform integration for all domains
        write(*,*) "Which data will be integrated in the selected domain?"
        write(*,*) "1 The grid data in memory"
        write(*,*) "2 Choose a real space function, which will be directly calculated"
        write(*,*) "3 Grid data loaded from a .cub file"
        read(*,*) inttype
        if (inttype==2) then
			write(*,*) "Select the real space function to be integrated, e.g. 1"
			call selfunc_interface(1,ifuncint)
        else if (inttype==3) then
			if (allocated(cubmattmp)) deallocate(cubmattmp)
            write(*,*) "Input path of the cube file, e.g. D:\ltwd\rho.cub"
            write(*,"(a)") " Note: Distribution of the grids in this file must be exactly identical to that of grid data in memory"
            do while(.true.)
				read(*,"(a)") cubname
				inquire(file=cubname,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again!"
			end do
            call readcubetmp(cubname,1,inconsis)
            if (inconsis==1) then
				write(*,*) "Error: Distribution of the grids in this file is not exactly identical to that of grid data in memory"
                write(*,*) "Press ENTER button to return"
                read(*,*)
                cycle
            end if
        end if
		write(*,*) "Domain    Integral (a.u.)     Volume (Bohr^3)      Average"
		valinttot=0
		volinttot=0
		do intdom=1,ndomain
			valint=0
			volint=0
			do igrd=1,domainsize(intdom)
				idx=domaingrid(igrd,intdom)
				if (inttype==1) then
					tmpval=cubmat(dogrid(1,idx),dogrid(2,idx),dogrid(3,idx))
                else if (inttype==2) then
					tmpval=calcfuncall(ifuncint,gridxyz(1,idx),gridxyz(2,idx),gridxyz(3,idx))
				else if (inttype==3) then
					tmpval=cubmattmp(dogrid(1,idx),dogrid(2,idx),dogrid(3,idx))
                end if
				valint=valint+tmpval
				volint=volint+1
			end do
			avgval=valint/domainsize(intdom)
			valint=valint*dvol
			volint=volint*dvol
			write(*,"(i6,E20.10,f17.6,E20.10)") intdom,valint,volint,avgval
			valinttot=valinttot+valint
			volinttot=volinttot+volint
		end do
		write(*,"(' Integration result of all domains:',E20.10,' a.u.')") valinttot
		write(*,"(' Volume of all domains:',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volinttot,volinttot*b2a**3
        
	else if (isel2==-2) then !Perform integration for subregion of some domains according to sign(lambda)*rho
		write(*,*) "Input indices of the domains you want to integrate, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, all domains will be selected"
		read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
			ntmp=ndomain
			allocate(tmparr(ntmp))
            forall(i=1:ntmp) tmparr(i)=i
        else
			call str2arr(c2000tmp,ntmp)
			allocate(tmparr(ntmp))
			call str2arr(c2000tmp,ntmp,tmparr)
        end if
        write(*,"(a)") " Integrate which part of the selected domains?"
        write(*,*) "0 All region"
        write(*,*) "1 Subregion with positive sign(lambda_2)*rho"
        write(*,*) "2 Subregion with negative sign(lambda_2)*rho"
        write(*,*) "3 Subregion within specific range of sign(lambda_2)*rho"
        read(*,*) iregion
        if (iregion==0) then
			rlow=-1D200
            rhigh=1D200
        else if (iregion==1) then
			rlow=0
            rhigh=1D200
        else if (iregion==2) then
			rlow=-1D200
            rhigh=0
        else if (iregion==3) then
			write(*,*) "Please input lower and upper limits of sign(lambda_2)*rho in a.u."
            write(*,*) "For example, -0.03,0.02"
            read(*,*) rlow,rhigh
        end if
		write(*,*) "Select the real space function to be integrated, e.g. 3"
		call selfunc_interface(1,ifuncint)
        write(*,*)
        write(*,*) "Meaning of outputted terms:"
        write(*,*) "int(tot): Total integral in specific subregion of selected domains, in a.u."
        write(*,*) "vol(tot): Total volume in specific subregion of selected domains, in Bohr"
        write(*,"(a)") " int(pos) and int(neg): The integrals where lambda_2 is positive and negative in the subregion, respectively. They sum to int(tot)"
        write(*,"(a)") " vol(pos) and vol(neg): The volumes where lambda_2 is positive and negative in the subregion, respectively. They sum to vol(tot)"
        write(*,"(/,a)") " Domain     int(tot)    vol(tot)    int(pos)    vol(pos)    int(neg)    vol(neg)"
        sum_rintneg=0
        sum_rintpos=0
        sum_volneg=0
        sum_volpos=0
		do itmp=1,ntmp !Loops all selected domains
			intdom=tmparr(itmp)
			rintneg=0
			rintpos=0
			volneg=0
			volpos=0
			do igrd=1,domainsize(intdom) !Loops all grid of this domain
				idx=domaingrid(igrd,intdom)
                xnow=gridxyz(1,idx)
                ynow=gridxyz(2,idx)
                znow=gridxyz(3,idx)
				call signlambda2rho_RDG(xnow,ynow,znow,sl2r,RDG,rho)
                if (sl2r>=rlow.and.sl2r<=rhigh) then
					valint=calcfuncall(ifuncint,xnow,ynow,znow)
					if (sl2r<0) then
						rintneg=rintneg+valint
						volneg=volneg+1
					else
						rintpos=rintpos+valint
						volpos=volpos+1
					end if
                end if
			end do
			rintneg=rintneg*dvol
			rintpos=rintpos*dvol
			volneg=volneg*dvol
			volpos=volpos*dvol
            write(*,"(i6,1x,E15.6,f9.3,E15.6,f9.3,E15.6,f9.3)") intdom,rintneg+rintpos,volneg+volpos,rintpos,volpos,rintneg,volneg
            sum_rintneg=sum_rintneg+rintneg
            sum_rintpos=sum_rintpos+rintpos
            sum_volneg=sum_volneg+volneg
            sum_volpos=sum_volpos+volpos
		end do
        write(*,"('  sum  ',E15.6,f9.3,E15.6,f9.3,E15.6,f9.3)") sum_rintneg+sum_rintpos,sum_volneg+sum_volpos,sum_rintpos,sum_volpos,sum_rintneg,sum_volneg
        deallocate(tmparr)
        
	else if (isel2==3) then !Visualize domains
		idrawdomain=1
		aug3Dold=aug3D
		if (aug3D<3) aug3D=3 !Often we set extension distance to zero, e.g. RDG, in this case the molecule will be truncated, therefore here temporarily augment it
		call drawdomaingui
		aug3D=aug3Dold
		idrawdomain=0
        
	else if (isel2==4) then !Sort domains according to volumes from large to small
		write(*,*) "Sorting, please wait..."
		allocate(tmparr(ngrid))
		do idom=1,ndomain
			do jdom=1,ndomain-idom
				if (domainsize(jdom)<domainsize(jdom+1)) then
					tmp=domainsize(jdom)
                    domainsize(jdom)=domainsize(jdom+1)
                    domainsize(jdom+1)=tmp
                    tmparr(:)=domaingrid(:,jdom+1)
                    domaingrid(:,jdom+1)=domaingrid(:,jdom)
                    domaingrid(:,jdom)=tmparr(:)
                    where(grididx==jdom+1) grididx=-1
                    where(grididx==jdom) grididx=jdom+1
                    where(grididx==-1) grididx=jdom
                end if
            end do
        end do
		deallocate(tmparr)
		write(*,*) "Domain information after sorting:"
		do idom=1,ndomain
			write(*,"(' Domain:',i7,'    Grids:',i9,'    Volume:',f11.4,' Angstrom^3')") idom,domainsize(idom),domainsize(idom)*dvol*b2a**3
		end do
        
	else if (isel2==5) then !Calculate q_bind index for a domain
		write(*,*) "Input the index of the domain to be integrated, e.g. 3"
		read(*,*) intdom
		if (intdom<1.or.intdom>ndomain) then
			write(*,"(a)") " Error: The index of the domain to be integrated is incorrect"
			cycle
		end if
        write(*,"(a)") "Input exponent of electron density. For example, if you input 1.1, then rho^1.1 will be taken as the integrand"
        write(*,*) "If pressing ENTER button directly, then 4/3 will be used"
        read(*,"(a)") c80tmp
        if (c80tmp==" ") then
			expfac=4D0/3D0
		else
			read(c80tmp,*) expfac
        end if
		qatt=0
		qrep=0
		volneg=0
		volpos=0
		do igrd=1,domainsize(intdom)
			idx=domaingrid(igrd,intdom)
			call signlambda2rho_RDG(gridxyz(1,idx),gridxyz(2,idx),gridxyz(3,idx),sl2r,RDG,rho)
			if (sl2r<0) then
				qatt=qatt+rho**expfac
				volneg=volneg+1
			else
				qrep=qrep+rho**expfac
				volpos=volpos+1
			end if
		end do
		qatt=qatt*dvol
		qrep=qrep*dvol
		qbind=-(qatt-qrep)
		volneg=volneg*dvol
		volpos=volpos*dvol
		write(*,"(' q_att: ',f16.8,' a.u.')") qatt
		write(*,"(' q_rep: ',f16.8,' a.u.')") qrep
		write(*,"(' q_bind:',f16.8,' a.u.')") qbind
		write(*,"(' Volume (lambda2<0):',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volneg
		write(*,"(' Volume (lambda2>0):',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volpos
		write(*,"(' Volume (Total):    ',f13.6,' Bohr^3  ',f13.6,' Angstrom^3')") volneg+volpos
        
	else if (isel2==10) then !Export a domain as domain.cub file in current folder
		write(*,*) "Input the index of the domain to be exported, e.g. 4"
		read(*,*) idomain
		write(*,*) "Outputting domain.cub..."
		open(10,file="domain.cub",status="replace")
		write(10,"(' Generated by Multiwfn')")
		write(10,"(' Totally ',i12,' grid points')") nx*ny*nz
		write(10,"(i5,3f12.6)") ncenter,orgx,orgy,orgz
		write(10,"(i5,3f12.6)") nx,dx,0.0,0.0
		write(10,"(i5,3f12.6)") ny,0.0,dy,0.0
		write(10,"(i5,3f12.6)") nz,0.0,0.0,dz
		do icenter=1,ncenter
			write(10,"(i5,4f12.6)") a(icenter)%index,a(icenter)%charge,a(icenter)%x,a(icenter)%y,a(icenter)%z
		end do
		icount=0
		do ix=1,nx
			do iy=1,ny
				do iz=1,nz
					if (ix==1.or.ix==nx.or.iy==1.or.iy==ny.or.iz==1.or.iz==nz) then
						write(10,"(1PE13.5)",advance="no") 0D0 !Boundary grid is set to 0, so that the isosurface will always be closed
					else if (grididx(ix,iy,iz)==idomain) then
						write(10,"(1PE13.5)",advance="no") 1D0
					else
						write(10,"(1PE13.5)",advance="no") 0D0
					end if
					icount=icount+1
					if (icount==6) then
						write(10,*)
						icount=0
					end if
				end do
			end do
		end do
		close(10)
		write(*,"(a)") " Done! domain.cub has been outputted to current folder. &
		&The grids belonging and not belonging the domain have value of 1 and 0, respectively"
        
	else if (isel2==11) then !Export boundary grids of a domain to domain.pdb file in current folder
		write(*,*) "Input index of the domain, e.g. 4"
		read(*,*) idomain
		write(*,*) "Outputting domain.pdb..."
		open(10,file="domain.pdb",status="replace")
		if (ifPBC>0) then
			call getcellabc(asize,bsize,csize,alpha,beta,gamma)
			write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
		end if
		do igrd=1,domainsize(idomain)
			idx=domaingrid(igrd,idomain)
			if (boundgrid(idx)) then
				write(10,"(a,i5,1x,a4,14x,3f8.3)") "HETATM",igrd," C  ",gridxyz(1,idx)*b2a,gridxyz(2,idx)*b2a,gridxyz(3,idx)*b2a
			end if
		end do
		close(10)
		write(*,"(a)") " Done! domain.pdb has been outputted to current folder"
        
    else if (isel2==12) then !Export X,Y,Z coordinates and value of all grids in each domain
		write(*,*) "Input index of the domain, e.g. 4"
		read(*,*) idomain
		write(*,*) "Outputting domain.txt..."
		open(10,file="domain.txt",status="replace")
		do igrd=1,domainsize(idomain)
			idx=domaingrid(igrd,idomain)
            ix=dogrid(1,idx)
			iy=dogrid(2,idx)
			iz=dogrid(3,idx)
            write(10,"(3f11.5,1PE16.8)") gridxyz(1,idx)*b2a,gridxyz(2,idx)*b2a,gridxyz(3,idx)*b2a,cubmat(ix,iy,iz)
		end do
		close(10)
		write(*,"(a)") " Done! domain.txt has been outputted to current folder"
        write(*,*) "Column 1,2,3: X,Y,Z (in Angstrom) of grids in this domain"
        write(*,*) "Column 4: Value of the grid data recorded in memory at corresponding grids"
        write(*,"(a,f10.5,' Bohr^3 (',f10.5,' Angstrom^3)')") " Note: Volume of each grid is",dvol,dvol*b2a**3
	end if
end do
end subroutine


!------ Calculate electron correlation index proposed by Matito et al.
subroutine elecorridx
use defvar
real*8 occ(nmo),I_ND,I_D,I_T
write(*,*) "See Section 4.A.7 of Multiwfn manual for details"
I_ND=0
I_D=0
if (wfntype==3) then
	occ=MOocc/2
	where(occ>1) occ=1 !Remove unphysical occupation number larger than unity
	where(occ<0) occ=0 !Remove unphysical negative occupation number
	do i=1,nmo
		I_D=I_D+ dsqrt(occ(i)*(1-occ(i))) - 2*occ(i)*(1-occ(i))
	end do
	I_D=I_D/4
	do i=1,nmo
		I_ND=I_ND+ occ(i)*(1-occ(i))
	end do
	I_ND=I_ND/2
	!Above we only consider half part, another part is identical to that, so double the result
	I_D=I_D*2
	I_ND=I_ND*2
else if (wfntype==4) then
	occ=MOocc
	where(occ>1) occ=1
	where(occ<0) occ=0
	do i=1,nmo
		I_D=I_D+ dsqrt(occ(i)*(1-occ(i))) - 2*occ(i)*(1-occ(i))
	end do
	I_D=I_D/4
	do i=1,nmo
		I_ND=I_ND+ occ(i)*(1-occ(i))
	end do
	I_ND=I_ND/2
end if
I_T=I_ND+I_D
write(*,"(' Nondynamic correlation index:',f12.8)") I_ND
write(*,"(' Dynamic correlation index:   ',f12.8)") I_D
write(*,"(' Total correlation index:     ',f12.8)") I_T
end subroutine



!!------ Generate natural orbitals based on the density matrix loaded from .fch/.fchk file
!gennatorb is invoked in this routine
subroutine fch_gennatorb
use util
use defvar
implicit real*8 (a-h,o-z)
real*8,allocatable :: Pspin(:,:)
character selectyn,denstype*10,locstr*40

if (ifiletype/=1) then
	write(*,*) "Error: .fch/.fchk should be used as input file for this function"
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
write(*,*) "Input the type of density matrix, e.g. SCF, CI, CC, MP2, MP3, MP4..."
write(*,"(a)") " e.g. If the .fch was produced by Gaussian at MP2 level and ""density"" keyword has been used, you may input ""SCF"" or ""MP2"""
write(*,"(a)") " Hint: You can manually open .fch file by text editor and search ""density"" to see which density matrices are available in this file"
read(*,"(a)") denstype
write(locstr,"('Total ',a,' Density')") trim(denstype)
open(10,file=filename,status="old")
call loclabel(10,trim(locstr),ifoundDM)
if (ifoundDM==0) then
	write(*,"(' Error: Unable to find ""',a,'"" from the input file')") trim(locstr)
	write(*,*) "Press ENTER button to return"
	read(*,*)
	return
end if
iNOtype=1
if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	write(*,*) "Select the type of natural orbitals you want to obtain"
	write(*,*) "1 Spatial natural orbitals (diagonalizing total density matrix)"
	write(*,*) "2 Alpha and beta natural orbitals (diagonalizing respective density matrix)"
	write(*,*) "3 Spin natural orbitals (diagonalizing spin density matrix)"
	read(*,*) iNOtype
end if

write(*,*) "Loading density matrix..."
!Load total density matrix
Ptot=0D0
call loclabel(10,trim(locstr))
read(10,*)
read(10,"(5(1PE16.8))") ((Ptot(i,j),j=1,i),i=1,nbasis)
Ptot=Ptot+transpose(Ptot)
do i=1,nbasis
	Ptot(i,i)=Ptot(i,i)/2D0
end do
!Load spin density matrix to construct alpha and beta DM
if (iNOtype>1) then
	allocate(Pspin(nbasis,nbasis))
	Pspin=0
	read(10,*)
	read(10,"(5(1PE16.8))") ((Pspin(i,j),j=1,i),i=1,nbasis)
	Pspin=Pspin+transpose(Pspin)
	do i=1,nbasis
		Pspin(i,i)=Pspin(i,i)/2D0
	end do
	Palpha=(Ptot+Pspin)/2D0
	Pbeta=(Ptot-Pspin)/2D0
end if
close(10)
write(*,*) "Density matrix was loaded from .fch/.fchk file"

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



!!------ Generate natural orbitals based on the density matrix in memory, wavefunction information including wfntype will be updated to NO case
!iNOtype=1: Spatial NO, =2: Alpha and beta NO, =3: Spin NO
!ioutmode=1: Print intermediate information =0: Do not print
subroutine gennatorb(iNOtype,ioutmode)
use util
use defvar
implicit real*8 (a-h,o-z)
integer iNOtype
real*8,allocatable :: tmparr(:),Pspin(:,:)
real*8 Xmat(nbasis,nbasis),Xmatinv(nbasis,nbasis),tmpmat(nbasis,nbasis)

!To produce natural orbitals, we need to convert P to orthogonalized basis and then diagonalize it
allocate(tmparr(nbasis))
if (ioutmode==1) write(*,*)
if (iNOtype==1.or.iNOtype==3) then
	if (iNOtype==1) then
		if (ioutmode==1) write(*,*) "Generating NOs, please wait..."
		call symmorthomat(Sbas,Xmat,Xmatinv)
		!call diagsymat(matmul(matmul(transpose(Xmat),Ptot),Xmat),CObasa,MOocc,ierror) !CObasa now is NOs in orthogonalized basis
        tmpmat=matmul_blas(Xmat,Ptot,nbasis,nbasis,1,0)
        tmpmat=matmul_blas(tmpmat,Xmat,nbasis,nbasis)
		call diagsymat(tmpmat,CObasa,MOocc,ierror) !CObasa now is NOs in orthogonalized basis
	else
		allocate(Pspin(nbasis,nbasis))
		Pspin=Palpha-Pbeta !Construct spin density matrix
		if (ioutmode==1) write(*,*) "Generating SNOs, please wait..."
		call symmorthomat(Sbas,Xmat,Xmatinv)
		!call diagsymat(matmul(matmul(transpose(Xmat),Pspin),Xmat),CObasa,MOocc,ierror) !CObasa now is SNOs in orthogonalized basis
        tmpmat=matmul_blas(Xmat,Pspin,nbasis,nbasis,1,0)
        tmpmat=matmul_blas(tmpmat,Xmat,nbasis,nbasis)
		call diagsymat(tmpmat,CObasa,MOocc,ierror) !CObasa now is NOs in orthogonalized basis
	end if
	MOene=0
	CObasa=matmul_blas(Xmatinv,CObasa,nbasis,nbasis) !Back convert CObasa to original basis
	!Sort NOs according to occupation number
	do i=1,nbasis
		do j=i+1,nbasis
			if (MOocc(i)<MOocc(j)) then
				tmpocc=MOocc(i)
				MOocc(i)=MOocc(j)
				MOocc(j)=tmpocc
				tmparr=CObasa(:,i)
				CObasa(:,i)=CObasa(:,j)
				CObasa(:,j)=tmparr
			end if
		end do
	end do
	if (wfntype==1.or.wfntype==4) then !Then wfntype will be 3, deallocate useless arrays and resize arrays
		deallocate(CObasb,Palpha,Pbeta,MOene,tmparr)
		allocate(MOene(nbasis))
		MOene=0
		allocate(tmparr(nmo))
		tmparr=MOocc
		deallocate(MOocc)
		allocate(MOocc(nbasis))
		MOocc=tmparr(1:nbasis)
		nmo=nbasis
	end if
	if (ioutmode==1) then
		write(*,*) "Occupation numbers:"
		write(*,"(6f12.6)") MOocc
    end if
	wfntype=3
else if (iNOtype==2) then
	if (ioutmode==1) write(*,*) "Generating alpha and beta NOs, please wait..."
	call symmorthomat(Sbas,Xmat,Xmatinv)
	!call diagsymat(matmul(matmul(transpose(Xmat),Palpha),Xmat),CObasa,MOocc(1:nbasis),ierror)
    tmpmat=matmul_blas(Xmat,Palpha,nbasis,nbasis,1,0)
    tmpmat=matmul_blas(tmpmat,Xmat,nbasis,nbasis)
	call diagsymat(tmpmat,CObasa,MOocc(1:nbasis),ierror)
    CObasa=matmul_blas(Xmatinv,CObasa,nbasis,nbasis) !Back convert CObasa to original basis
	MOene(1:nbasis)=0
	do i=1,nbasis
		do j=i+1,nbasis
			if (MOocc(i)<MOocc(j)) then
				tmpocc=MOocc(i)
				MOocc(i)=MOocc(j)
				MOocc(j)=tmpocc
				tmparr=CObasa(:,i)
				CObasa(:,i)=CObasa(:,j)
				CObasa(:,j)=tmparr
			end if
		end do
	end do
	if (ioutmode==1) then
		write(*,*) "Occupation numbers of alpha NOs:"
		write(*,"(6f12.6)") MOocc(1:nbasis)
		write(*,*)
    end if
	call symmorthomat(Sbas,Xmat,Xmatinv)
	!call diagsymat(matmul(matmul(transpose(Xmat),Pbeta),Xmat),CObasb,MOocc(nbasis+1:nmo),ierror)
    tmpmat=matmul_blas(Xmat,Pbeta,nbasis,nbasis,1,0)
    tmpmat=matmul_blas(tmpmat,Xmat,nbasis,nbasis)
	call diagsymat(tmpmat,CObasb,MOocc(nbasis+1:nmo),ierror)
    CObasb=matmul_blas(Xmatinv,CObasb,nbasis,nbasis) !Back convert CObasb to original basis
	MOene(nbasis+1:nmo)=0
	do i=1,nbasis
		ii=i+nbasis
		do j=i+1,nbasis
			jj=j+nbasis
			if (MOocc(ii)<MOocc(jj)) then
				tmpocc=MOocc(ii)
				MOocc(ii)=MOocc(jj)
				MOocc(jj)=tmpocc
				tmparr=CObasb(:,i)
				CObasb(:,i)=CObasb(:,j)
				CObasb(:,j)=tmparr
			end if
		end do
	end do
	if (ioutmode==1) then
		write(*,*) "Occupation numbers of beta NOs:"
		write(*,"(6f12.6)") MOocc(nbasis+1:nmo)
		write(*,*)
    end if
	wfntype=4
end if

nindbasis=nbasis !After generating NOs, all orbitals are actually filled, so nindbasis must be set to nbasis
end subroutine




!!--------- Calculate core-valence bifurcation (CVB) index and related quantities
subroutine CVB_index
use defvar
use functions
implicit real*8 (a-h,o-z)
integer,parameter :: nptELFcurve=6000 !The number of points comprising the ELF curve,it is adequate to find exact ELF_CV and ELF_DHA
real*8 ELF_x(nptELFcurve),ELF_y(nptELFcurve)

write(*,*) "Original paper of CVB index: Theor. Chem. Acc., 104, 13 (2000)"
write(*,*)
write(*,*) "------ Calculating core-valence bifurcation (CVB) and related quantities -----"
write(*,*) "Input index of donor atom, hydrogen and acceptor atom in the H-bond (D-H...A)"
write(*,*) "For example: 1,3,4"
read(*,*) iD,iH,iA

!First time: calculate and analyze D<-H ELF curve
!Second time: calculate and analyze H->A ELF curve
ELF_DHA_x=0
do itime=1,2
	orgx1D=a(iH)%x
	orgy1D=a(iH)%y
	orgz1D=a(iH)%z
	if (itime==1) then
		endx1D=a(iD)%x
		endy1D=a(iD)%y
		endz1D=a(iD)%z
	else
		endx1D=a(iA)%x
		endy1D=a(iA)%y
		endz1D=a(iA)%z
	end if
	transx=(endx1D-orgx1D)/nptELFcurve
	transy=(endy1D-orgy1D)/nptELFcurve
	transz=(endz1D-orgz1D)/nptELFcurve
	transr=dsqrt(transx**2+transy**2+transz**2)
	!$OMP parallel do shared(ELF_x,ELF_y) private(ipt,rnowx,rnowy,rnowz) num_threads(nthreads)
	do ipt=1,nptELFcurve
		rnowx=orgx1D+(ipt-1)*transx
		rnowy=orgy1D+(ipt-1)*transy
		rnowz=orgz1D+(ipt-1)*transz
		ELF_x(ipt)=ipt*transr
		ELF_y(ipt)=ELF_LOL(rnowx,rnowy,rnowz,"ELF")
	end do
	!$OMP end parallel do
			
	!Find minimum
	do ipt=2,nptELFcurve-1
		gradold=ELF_y(ipt)-ELF_y(ipt-1)
		gradnew=ELF_y(ipt+1)-ELF_y(ipt)
		if (gradold*gradnew<0D0.and.gradnew>gradold) then !Find minimum
			if (itime==1) then !First minimum, C-V of donor atom
				ELF_CV_D=ELF_y(ipt)
				ELF_CV_x_D=ELF_x(ipt)
				exit
			else
				if (ELF_DHA_x==0) then !First minimum, bifurcation at DH-A
					ELF_DHA=ELF_y(ipt)
					ELF_DHA_x=ELF_x(ipt)
				else !Second minimum, C-V of acceptor atom
					ELF_CV_A=ELF_y(ipt)
					ELF_CV_x_A=ELF_x(ipt)
					exit
				end if
			end if
		end if
	end do
	
	if (itime==1) then
		write(*,"(' Core-valence bifurcation value at donor, ELF(C-V,D):',f8.4)") ELF_CV_D
		write(*,"(' Distance between corresponding minimum and the hydrogen:',f8.3,' Angstrom')") ELF_CV_x_D*b2a
		write(*,*)
	else
		write(*,"(' Core-valence bifurcation value at acceptor, ELF(C-V,A):',f8.4)") ELF_CV_A
		write(*,"(' Distance between corresponding minimum and the hydrogen:',f8.3,' Angstrom')") ELF_CV_x_A*b2a
		write(*,*)
		write(*,"(' Bifurcation value at H-bond, ELF(DH-A):',f8.4)") ELF_DHA
		write(*,"(' Distance between corresponding minimum and the hydrogen:',f8.3,' Angstrom')") ELF_DHA_x*b2a
		write(*,*)
	end if
end do

!ELF_CV=max(ELF_CV_D,ELF_CV_A)
!write(*,"(' ELF(C-V): ',f12.6)") ELF_CV
!write(*,*)
!write(*,"(' CVB index:',f12.6)") ELF_CV-ELF_DHA
write(*,"(' The CVB index, namely ELF(C-V,D) - ELF(DH-A):',f12.6)") ELF_CV_D - ELF_DHA
end subroutine







!!------------------------------------------------------------------------------------------
!!---- Fit Fukui function or other kind of density difference to orbital representation ----
!!------------------------------------------------------------------------------------------
!Note: Other kind of grid data is also acceptable
subroutine orbfitEDD
use defvar
use util
use GUI
use functions
implicit real*8 (a-h,o-z)
character c200tmp*200,c2000tmp*2000,selectyn
integer,allocatable :: orbidx(:)
real*8,allocatable :: Amat(:,:),Amatinv(:,:),Bvec(:),fvec(:),fval(:),EDD(:,:,:),orbgrid(:,:,:,:)
real*8 orbval(nmo)
integer,allocatable :: idxlist(:)
integer :: isetcons=1,imode=1,ioutfitcub=1
real*8 :: consval=1

write(*,*)
write(*,"(a)") " Input path of a cube file containing density difference (or other kind of grid data), e.g. C:\rize\f+.cub"
do while(.true.)
    c200tmp="C:\Users\Sobereva\Desktop\NBO_Fukui\H2CO\f-.cub" !!!!!!!!!!!!!!!!!!!!!!!!!!!!
	read(*,"(a)") c200tmp
	inquire(file=c200tmp,exist=alive)
	if (alive) exit
	write(*,*) "Cannot find the file, input again!"
end do
call readcube(c200tmp,1,1)
write(*,*) "Loading finished"

!Store EDD
allocate(EDD(nx,ny,nz))
EDD=cubmat

do while(.true.)
    write(*,*)
    write(*,*) "  ----------- Calculation of orbital contributions to grid data -----------"
    write(*,*) "-1 Return"
    write(*,*) " 0 Choose orbital range and start analysis!"
    if (isetcons==1) write(*,"(a,f12.3)") "  1 Set constraint on the sum of contributions, current:",consval
    if (isetcons==0) write(*,*) " 1 Set constraint on the sum of contributions, current: No constraint"
    if (imode==1) write(*,*) " 2 Switch calculation mode, current: fast (memory based)"
    if (imode==2) write(*,*) " 2 Switch calculation mode, current: slow (cube file based)"
    read(*,*) isel
    
    if (isel==-1) then
        exit
        
    else if (isel==1) then
        write(*,*) "1 Do not set constraint on the sum of contributions"
        write(*,*) "2 Set constraint on the sum of contributions to a specific value"
        read(*,*) isel2
        if (isel2==1) then
            isetcons=0
        else if (isel2==2) then
            isetcons=1
            write(*,*) "Input the value of constraint, e.g. 1.5"
            read(*,*) consval
        end if
        
    else if (isel==2) then
        if (imode==1) then
            imode=2
        else
            imode=1
        end if
        
    else if (isel==0) then !Start analysis!
        write(*,*)
        write(*,*) "Input index of the orbitals to be taken into account, e.g. 2,3,7-10"
        write(*,*) "If pressing ENTER button directly, all orbitals will be taken into account"
        write(*,"(a)") " If inputting ""o"", all orbitals with non-zero occupation will be taken into account"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            norb=nmo
            allocate(orbidx(nmo))
            forall(i=1:nmo) orbidx(i)=i
        else if (c2000tmp=="o") then
            norb=count(MOocc/=0)
            allocate(orbidx(norb))
            itmp=0
            do imo=1,nmo
                if (MOocc(imo)/=0) then
                    itmp=itmp+1
                    orbidx(itmp)=imo
                end if
            end do
        else
            call str2arr(c2000tmp,norb)
            allocate(orbidx(norb))
            call str2arr(c2000tmp,norb,orbidx)
        end if

        if (imode==1) then !Calculate grid data of |psi|^2 for all selected orbitals and store to memory
            if (allocated(orbgrid)) deallocate(orbgrid)
            allocate(orbgrid(nx,ny,nz,norb))
            ilow=minval(orbidx)
            ihigh=maxval(orbidx)
            write(*,*) "Calculating orbital density..."
            !$OMP PARALLEL DO SHARED(orbgrid) PRIVATE(i,j,k,tmpx,tmpy,tmpz,idx,iorb,orbval) schedule(dynamic) NUM_THREADS(nthreads)
            do k=1,nz
	            do j=1,ny
		            do i=1,nx
                        call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                        call orbderv(1,ilow,ihigh,tmpx,tmpy,tmpz,orbval)
                        do idx=1,norb
                            iorb=orbidx(idx)
			                orbgrid(i,j,k,idx)=orbval(iorb)**2
                        end do
		            end do
	            end do
            end do
            !$OMP END PARALLEL DO
        
        else if (imode==2) then !Calculate and export cube file of |psi|^2 for all selected orbitals
            do idx=1,norb
                iorb=orbidx(idx)
                write(c200tmp,"('rho_',i5.5,'.cub')") iorb
                inquire(file=c200tmp,exist=alive)
	            if (alive) then
                    write(*,"(1x,a)") trim(c200tmp)//" has already existed and thus will not be calculated"
                else
                    write(*,"(' Calculating orbital density for orbital',i6)") iorb
                    call savecubmat(4,1,iorb)
                    cubmat=cubmat**2
                    open(10,file=c200tmp,status="replace")
                    call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
                    close(10)
                end if
            end do
        end if

        !Construct A matrix and B vector
        if (isetcons==0) then
            ndim=norb
        else if (isetcons==1) then
            ndim=norb+1
        end if
        allocate(Amat(ndim,ndim),Amatinv(ndim,ndim),Bvec(ndim),fvec(ndim),fval(norb),idxlist(norb))
        write(*,*) "Calculating fitted coefficients..."
        call showprog(0,norb)
        do idx=1,norb
            iorb=orbidx(idx)
            
            if (imode==1) then !Memory based
                !Construct A matrix. Since it is symmetric, only construct upper-right part
                do jdx=idx,norb
                    tmpA=sum(orbgrid(:,:,:,idx)*orbgrid(:,:,:,jdx))
                    Amat(idx,jdx)=tmpA
                    Amat(jdx,idx)=tmpA
                end do
                !Construct B vector
                Bvec(idx)=sum(orbgrid(:,:,:,idx)*EDD)
                
            else if (imode==2) then !cube file based
                write(c200tmp,"('rho_',i5.5,'.cub')") iorb
                call readcube(c200tmp,2,1)
                !Construct A matrix. Since it is symmetric, only construct upper-right part
                do jdx=idx,norb
                    jorb=orbidx(jdx)
                    write(c200tmp,"('rho_',i5.5,'.cub')") jorb
                    call readcubetmp(c200tmp,2,itmp)
                    tmpA=sum(cubmat*cubmattmp)
                    Amat(idx,jdx)=tmpA
                    Amat(jdx,idx)=tmpA
                end do
                !Construct B vector
                Bvec(idx)=sum(cubmat*EDD)
            end if
    
            call showprog(idx,norb)
        end do

        !Finalize matrix and solve linear equation
        if (isetcons==1) then
            Amat(ndim,1:norb)=1
            Amat(1:norb,ndim)=1
            Amat(ndim,ndim)=0
            Bvec(ndim)=consval
        end if
        Amatinv=invmat(Amat,ndim)
        fvec=matmul(Amatinv,Bvec)

        !Show result
        forall(i=1:norb) idxlist(i)=i
        fval=fvec(1:norb)
        call sortr8(fval,list=idxlist)
        write(*,*)
        do idx=1,norb
            idxold=idxlist(idx)
            write(*,"(' Orbital',i6,'   Value:',f10.3)") orbidx(idxold),fval(idx)
        end do
        write(*,"(' Sum of all values:',f12.3)") sum(fval)

        !Show fitting error
        cubmat=0
        do idx=1,norb
            if (imode==1) then
                cubmat=cubmat+fvec(idx)*orbgrid(:,:,:,idx)
            else if (imode==2) then
                iorb=orbidx(idx)
                write(c200tmp,"('rho_',i5.5,'.cub')") iorb
                call readcubetmp(c200tmp,2,itmp)
                cubmat=cubmat+fvec(idx)*cubmattmp
            end if
        end do
        call calc_dvol(dvol)
        write(*,"(' Fitting error (definition 1):',f12.4)") sum(abs(EDD-cubmat))*dvol
        write(*,"(' Fitting error (definition 2):',f12.6)") sum(abs(EDD-cubmat)**2)*dvol

        do while(.true.)
            write(*,*)
            write(*,*) "0 Exit"
            write(*,*) "1 Output fitted grid data to fitted.cub in current folder"
            write(*,*) "2 Visualize isosurface of provided grid data"
            write(*,*) "3 Visualize isosurface of fitted grid data"
            write(*,"(a)") " 4 Visualize isosurface of difference between the provided grid data and the fitted grid data"
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==1) then
                open(10,file="fitted.cub",status="replace")
                call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
                close(10)
                write(*,*) "Outputting finished!"
            else if (isel2==2) then
                if (.not.allocated(cubmattmp)) allocate(cubmattmp(nx,ny,nz))
                cubmattmp=cubmat
                cubmat=EDD
                sur_value=0.01D0
                call drawisosurgui(1)
                cubmat=cubmattmp
            else if (isel2==3) then
                sur_value=0.01D0
                call drawisosurgui(1)
            else if (isel2==4) then
                if (.not.allocated(cubmattmp)) allocate(cubmattmp(nx,ny,nz))
                cubmattmp=cubmat
                cubmat=EDD-cubmat
                sur_value=0.005D0
                call drawisosurgui(1)
                cubmat=cubmattmp
            end if
        end do

        if (imode==2) then
            write(*,*) "Do you want to clean all .cub files involved in this run? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y') then
                do idx=1,norb
                    iorb=orbidx(idx)
                    write(c200tmp,"('rho_',i5.5,'.cub')") iorb
                    inquire(file=c200tmp,exist=alive)
	                if (alive) then
                        open(10,file=c200tmp,status="old")
                        close(10,status="delete")
                    end if
                end do
                write(*,*) "Done! All rho_xxxxx.cub files in current folder have been deleted"
            end if
        end if
        
        deallocate(orbidx,Amat,Amatinv,Bvec,fvec,fval,idxlist)
 
   end if
    
end do

end subroutine




!!--------------------------------------------------------------------------------
!!--------- Calculate Coulomb and exchange integral between two orbitals ---------
!!--------------------------------------------------------------------------------
subroutine orb_coulexcint
use defvar
use util
implicit real*8 (a-h,o-z)
character c200tmp*200
real*8,allocatable :: cubx(:),cuby(:),cubz(:),rhoii(:,:,:),rhojj(:,:,:),rhoij(:,:,:)
real*8 :: Coulcrit=1D-6,exccrit=1D-5 !Only leads to marginal error, speed may increase several times

if (allocated(b)) then !cubmat and cubmattmp will record orbital wavefunction grid data of j and i orbitals
    write(*,*) "Input index of two orbitals, e.g. 4,10"
    read(*,*) iorb,jorb
    call setgridfixspc
    write(*,*) "Calculating grid data of orbital wavefunction..."
	if (allocated(cubmat)) deallocate(cubmat)
	if (allocated(cubmattmp)) deallocate(cubmattmp)
	allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
	call savecubmat(4,1,iorb)
    cubmattmp=cubmat
	call savecubmat(4,1,jorb)
else
    write(*,*) "Input path of another cube file recording wavefunction of an orbital"
    write(*,*) "e.g. C:\otoboku\MO10.cub"
    do while(.true.)
	    read(*,"(a)") c200tmp
	    inquire(file=c200tmp,exist=alive)
	    if (alive) exit
	    write(*,*) "Cannot find the file, input again!"
    end do
    call readcubetmp(c200tmp,1,itmp)
end if
call calc_dvol(dvol)

allocate(rhoii(nx,ny,nz),rhojj(nx,ny,nz),rhoij(nx,ny,nz))
rhoii=cubmattmp**2
rhojj=cubmat**2
rhoij=cubmattmp*cubmat

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

do while(.true.)
    write(*,*)
    write(*,*) "------- Calculate Coulomb and exchange integrals based on uniform grid -------"
    write(*,*) "0 Return"
    write(*,*) "1 Calculate Coulomb integral"
    write(*,"(a,f12.8)") " 2 Set truncation value of Coulomb integral, current:",Coulcrit
    write(*,*) "3 Calculate exchange integral"
    write(*,"(a,f12.8)") " 4 Set truncation value of exchange integral, current:",exccrit
    read(*,*) isel
    if (isel==0) then
        exit
    else if (isel==2) then
        write(*,*) "Input a value, e.g. 1E-6"
        read(*,*) Coulcrit
    else if (isel==4) then
        write(*,*) "Input a value, e.g. 1E-5"
        read(*,*) exccrit
    end if

    if (isel==1.or.isel==3) then
        write(*,*) "Calculating integrals, please wait..."
        call walltime(iwalltime1)
        coulene=0
        excene=0
        do k=1,nz
			cubzk=cubz(k)
	        do j=1,ny
				cubyj=cuby(j)
		        do i=1,nx
					cubxi=cubx(i)
                    if (isel==1) then !Coulomb integral
			            if (rhoii(i,j,k)>Coulcrit) then
			                !$OMP parallel shared(coulene) private(ii,jj,kk,distx2,disty2,distz2,dist,coulenetmp) num_threads(nthreads)
			                coulenetmp=0
			                !$OMP do schedule(DYNAMIC)
					        do kk=1,nz
						        distz2=(cubzk-cubz(kk))**2
				                do jj=1,ny
					                disty2=(cubyj-cuby(jj))**2
			                        do ii=1,nx
				                        distx2=(cubxi-cubx(ii))**2
                                        dist=dsqrt(distx2+disty2+distz2)
                                        if (dist==0) cycle
						                coulenetmp=coulenetmp + rhojj(ii,jj,kk) /dist
					                end do
				                end do
			                end do
			                !$OMP END DO
			                !$OMP CRITICAL
			                coulene=coulene+rhoii(i,j,k)*coulenetmp
			                !$OMP END CRITICAL
			                !$OMP END PARALLEL
                        end if
                    else if (isel==3) then !Exchange integral
			            if (abs(rhoij(i,j,k))>exccrit) then
			                !$OMP parallel shared(excene) private(ii,jj,kk,distx2,disty2,distz2,dist,excenetmp) num_threads(nthreads)
                            excenetmp=0
			                !$OMP do schedule(DYNAMIC)
					        do kk=1,nz
						        distz2=(cubz(k)-cubz(kk))**2
				                do jj=1,ny
					                disty2=(cuby(j)-cuby(jj))**2
			                        do ii=1,nx
				                        distx2=(cubx(i)-cubx(ii))**2
						                dist=dsqrt(distx2+disty2+distz2)
                                        if (dist==0) cycle
						                excenetmp=excenetmp + rhoij(ii,jj,kk) /dist
					                end do
				                end do
			                end do
			                !$OMP END DO
			                !$OMP CRITICAL
                            excene=excene+rhoij(i,j,k)*excenetmp
			                !$OMP END CRITICAL
			                !$OMP END PARALLEL
                        end if
                    end if
                    
		        end do
	        end do
	        call showprog(k,nz)
        end do
        
        coulene=coulene*dvol*dvol
        excene=excene*dvol*dvol
        call walltime(iwalltime2)
        write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
        write(*,*)
        if (isel==1) write(*,"(' Coulomb integral (ii|jj): ',f12.6,' a.u.')") coulene
        if (isel==3) write(*,"(' Exchange integral (ij|ji):',f12.6,' a.u.')") excene
    end if
    
end do
end subroutine





!!---------------------------------------------------------------
!!------------- Calculate bond length/order alternation (BLA/BOA)
!!---------------------------------------------------------------
subroutine BLABOA
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000,selectyn
integer,allocatable :: chainatm(:),atmseq(:),atmtmp(:)
real*8,allocatable :: PSmat(:,:),PSmata(:,:),PSmatb(:,:)
integer :: cenind(2000)

write(*,*) "Input atom indices in the chain (the sequence is arbitrary)"
write(*,*) "e.g. 2,14,16-17,19,21,23-24"
write(*,*) "If pressing enter button directly, all atoms will be selected"
read(*,"(a)") c2000tmp
if (c2000tmp==" ") then
	nchainatm=ncenter
	allocate(chainatm(nchainatm),atmseq(nchainatm),atmtmp(ncenter))
    forall(i=1:ncenter) chainatm(i)=i
else
	call str2arr(c2000tmp,nchainatm)
	allocate(chainatm(nchainatm),atmseq(nchainatm),atmtmp(ncenter))
	call str2arr(c2000tmp,nchainatm,chainatm)
end if

write(*,*) "Input index of the two atoms at the two ends of the path, e.g. 13,24"
write(*,"(a)") " If the path is a closed path (e.g. a ring), input twice of starting atom index, e.g. 5,5"
read(*,*) ibeg,iend

if (.not.allocated(connmat)) call genconnmat(1,0) !Generate connectivity matrix

!Identify the atom sequence in the chain
!From ibeg, gradually add adjacent atom to the sequence, until the iend is encountered
atmtmp(:)=0 !All atoms have not been added to the sequence. If atmtmp(i)=1, that means the atom i has already been added to the sequence
atmtmp(ibeg)=1
atmseq(1)=ibeg !The ibeg is the first atom in the sequence
inow=ibeg !Atom index of current step
idx=1 !Current position in the sequence
do while(.true.)
    do itmp=1,nchainatm !Compare current atom to other atoms
        iatm=chainatm(itmp)
        if (atmtmp(iatm)==1) cycle !This atom has already added, skip it
        if (connmat(inow,iatm)/=0) then
            inow=iatm !Move to this atom
            atmtmp(iatm)=1
            idx=idx+1
            atmseq(idx)=inow
            exit
        end if
    end do
    if (ibeg/=iend) then !Open path
        if (inow==iend) exit
    else !Closed path
        if (idx>2.and.connmat(inow,ibeg)/=0) exit
    end if
end do

write(*,*)
write(*,*) "Sequence of the atoms in the chain from the beginning side to the ending side"
write(*,"(9i8)") atmseq

open(10,file="bondalter.txt",status="replace")

write(*,*)
if (allocated(CObasa)) then
	call ask_Sbas_PBC
    iBO=1
    if (wfntype==0.or.wfntype==3) then !Closed-shell
        allocate(PSmat(nbasis,nbasis))
        PSmat=matmul(Ptot,Sbas)
    else !Open-shell
        allocate(PSmata(nbasis,nbasis),PSmatb(nbasis,nbasis))
	    PSmata=matmul(Palpha,Sbas)
	    PSmatb=matmul(Pbeta,Sbas)
    end if
    write(*,*) " Bond     Atom1     Atom2   Length (Angstrom)   Mayer bond order"
else
    iBO=0
    write(*,*) " Bond     Atom1     Atom2   Length (Angstrom)"
end if

avglen_even=0
avglen_odd=0
avgBO_even=0
avgBO_odd=0
n_even=0
n_odd=0
idxend=nchainatm-1
if (ibeg==iend) idxend=nchainatm
do idx=1,idxend
    iatm=atmseq(idx)
    if (ibeg/=iend) then !Open path
        jatm=atmseq(idx+1)
    else !Closed path
        if (idx<idxend) then
            jatm=atmseq(idx+1)
        else if (idx==idxend) then
            jatm=atmseq(1)
        end if
    end if
    dist=atomdist(iatm,jatm,1)*b2a
    if (mod(idx,2)==1) then !odd
        avglen_odd=avglen_odd+dist
        n_odd=n_odd+1
    else !even
        avglen_even=avglen_even+dist
        n_even=n_even+1
    end if
    
    if (iBO==0) then
        write(*,"(i5,2i10,f14.3)") idx,iatm,jatm,dist
        write(10,"(i5,2i10,f14.3)") idx,iatm,jatm,dist
    else !Also calculate bond order
        cenind(1)=iatm
        cenind(2)=jatm
        if (wfntype==0.or.wfntype==3) then !Closed-shell
            call calcmultibndord(2,cenind,PSmat,nbasis,bondorder)
        else !Open-shell
            call calcmultibndord(2,cenind,PSmata,nbasis,resulta)
            call calcmultibndord(2,cenind,PSmatb,nbasis,resultb)
            bondorder=2*(resulta+resultb)
        end if
        write(*,"(i5,2i10,f14.4,f20.4)") idx,iatm,jatm,dist,bondorder
        write(10,"(i5,2i10,f14.4,f20.4)") idx,iatm,jatm,dist,bondorder
        if (mod(idx,2)==1) then !odd
            avgBO_odd=avgBO_odd+bondorder
        else
            avgBO_even=avgBO_even+bondorder
        end if
    end if
end do

write(*,*)
write(*,*) "The data shown above have also been exported to bondalter.txt in current folder"
write(*,"(a,i6)") " The number of even bonds:",n_even
write(*,"(a,i6)") " The number of odd bonds: ",n_odd
avglen_even=avglen_even/n_even
avglen_odd=avglen_odd/n_odd
write(*,"(a,f12.4,' Angstrom')") " Average length of even bonds: ",avglen_even
write(*,"(a,f12.4,' Angstrom')") " Average length of odd bonds:  ",avglen_odd
BLA=avglen_even-avglen_odd
write(*,"(a,f12.4,' Angstrom')") " Bond length alternation (BLA):",BLA
if (iBO==1) then
    avgBO_even=avgBO_even/n_even
    avgBO_odd=avgBO_odd/n_odd
    write(*,"(a,f12.4)") " Average bond order of even bonds:",avgBO_even
    write(*,"(a,f12.4)") " Average bond order of odd bonds: ",avgBO_odd
    BOA=avgBO_even-avgBO_odd
    write(*,"(a,f12.4)") " Bond order alternation (BOA):    ",BOA
end if

write(*,"(/,a)") " Do you also want to calculate variation of bond angle and dihedral along the path? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
    write(*,*) "Note The unit of printed values is degree"
    write(*,*)
    !Angles
	avgang_even=0
	avgang_odd=0
	n_even=0
	n_odd=0
    iangle=0
    idxend=nchainatm
    if (ibeg==iend) idxend=nchainatm+2
    do itmp=3,idxend
		iangle=iangle+1
        idx=itmp-2
        iatm=atmseq(idx)
        jdx=itmp-1
        if (jdx>nchainatm) jdx=jdx-nchainatm
        jatm=atmseq(jdx)
        kdx=itmp
        if (kdx>nchainatm) kdx=kdx-nchainatm
        katm=atmseq(kdx)
        angle=xyz2angle(a(iatm)%x,a(iatm)%y,a(iatm)%z,a(jatm)%x,a(jatm)%y,a(jatm)%z,a(katm)%x,a(katm)%y,a(katm)%z)
        write(*,"(' #',i5,'  Atoms:',3i6,'  Angle:',f10.3)") iangle,iatm,jatm,katm,angle
		if (mod(iangle,2)==1) then !odd
			avgang_odd=avgang_odd+angle
			n_odd=n_odd+1
		else !even
			avgang_even=avgang_even+angle
			n_even=n_even+1
		end if
    end do
	write(*,"(a,i6)") " The number of even angles:",n_even
	write(*,"(a,i6)") " The number of odd angles: ",n_odd
	avgang_even=avgang_even/n_even
	avgang_odd=avgang_odd/n_odd
	write(*,"(a,f12.4)") " Average degree of even angles: ",avgang_even
	write(*,"(a,f12.4)") " Average degree of odd angles:  ",avgang_odd
    
    !Dihedral
    write(*,*)
    idih=0
    if (ibeg==iend) idxend=nchainatm+3
    do itmp=4,idxend
		idih=idih+1
        idx=itmp-3
        iatm=atmseq(idx)
        jdx=itmp-2
        if (jdx>nchainatm) jdx=jdx-nchainatm
        jatm=atmseq(jdx)
        kdx=itmp-1
        if (kdx>nchainatm) kdx=kdx-nchainatm
        katm=atmseq(kdx)
        ldx=itmp
        if (ldx>nchainatm) ldx=ldx-nchainatm
        latm=atmseq(ldx)
        dih=abs(atomdih(iatm,jatm,katm,latm,1)) !Note that the returned value of is [0,180]
        dev=dih
        if (dih>90) dev=180-dih
        write(*,"(' #',i5,'  Atoms:',4i6,'  Dih.:',f7.2,', dev. to planar:',f7.2)") idih,iatm,jatm,katm,latm,dih,dev
    end do
end if

end subroutine





!!--------------------------------------------------------------------------------------
!!------  Bond order density (BOD) and natural adaptive orbitals (NAdOs) analyses ------
!!--------------------------------------------------------------------------------------
subroutine BOD
use defvar
use util
implicit real*8 (a-h,o-z)
real*8,allocatable :: BODmat(:,:),BODmatb(:,:)
real*8,allocatable :: BOM1(:,:),BOM2(:,:),AOM1(:,:),AOM2(:,:),FOM1(:,:),FOM2(:,:)
real*8,allocatable :: BOM1b(:,:),BOM2b(:,:),AOM1b(:,:),AOM2b(:,:),FOM1b(:,:),FOM2b(:,:)
real*8,allocatable :: eigvecmat(:,:),eigvalarr(:)
real*8 tmparr(nbasis),mat11(1,1)
integer frag1atm(ncenter),frag2atm(ncenter)
character c80tmp*80,c200tmp*200,c2000tmp*200,selectyn

if (wfntype>=2) then !Only for R and U SCF currently
    write(*,"(a)") " Error: Only restricted and unresticted single-determinant wavefunction are supported!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

if (.not.allocated(CObasa)) then
    write(*,"(a)") " Error: The input file does not contain basis function information!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

if (ifdelvirorb==1) call delvirorb_back(1)

if (all(MOene==0)) then
	write(*,*) "Warning: Energies of all orbitals are zero, so you need to manually input HOMO index"
	if (wfntype==0) then
		write(*,*) "Input index of HOMO, e.g. 300"
		read(*,*) idxHOMO
    else if (wfntype==1) then
		write(*,*) "Input index of HOMO of alpha spin, e.g. 300"
		read(*,*) idxHOMO
		write(*,*) "Input index of HOMO of beta spin, e.g. 299"
		read(*,*) idxHOMOb
        idxHOMOb=idxHOMOb+nbasis
    end if
else
	call getHOMOidx
end if
idoene=0

do while(.true.)
	write(*,*)
	write(*,*) "   ======= Bond order density and natural adaptive orbital analysis ======="
	if (idoene==0) write(*,*) "-1 Toggle if calculating energies for NAdOs, current: No"
	if (idoene==1) write(*,*) "-1 Toggle if calculating energies for NAdOs, current: Yes"
	write(*,*) " 0 Return"
	write(*,"(a)") "  1 Interatomic interaction analysis based on atomic overlap matrix (AOM)"
	write(*,"(a)") "  2 Interbasin interaction analysis based on basin overlap matrix (BOM)"
	write(*,"(a)") "  3 Interfragment interaction analysis based on the fragment overlap matrix (FOM) constructed from AOM"
	write(*,"(a)") "  4 Interfragment interaction analysis based on FOM directly provided in FOM.txt"
	read(*,*) isel

	if (isel==0) then
		return
    
	else if (isel==-1) then
		if (idoene==1) then
			idoene=0
		else if (idoene==0) then
			do while(.true.)
				write(*,*)
				write(*,*) "How to evaluate energy of NAdOs?"
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
    else
		exit
	end if
end do
    
if (isel==1) then !Interatomic interaction analysis based on atomic overlap matrix (AOM)
    write(*,*) "Input the path of the file containing AOM, e.g. C:\AOM.txt"
    write(*,*) "If press ENTER button directly, AOM.txt in current folder will be loaded"
    read(*,"(a)") c200tmp
    if (c200tmp==" ") c200tmp="AOM.txt"
    write(*,*) "Input index of the two atoms, e.g. 2,4"
    read(*,*) iatm1,iatm2
    allocate(AOM1(idxHOMO,idxHOMO),AOM2(idxHOMO,idxHOMO))
    open(10,file=c200tmp,status="old")
    if (wfntype==0) then !Closed-shell
        do itime=1,2 !Two atoms
            if (itime==1) iatm=iatm1
            if (itime==2) iatm=iatm2
            write(c80tmp,"('Atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call loclabel(10,c80tmp,ifound)
            if (ifound==0) then
                write(*,"(' Error: Unable to find AOM of atom',i6)") iatm
                write(*,*) "Press ENTER button to exit"
                read(*,*)
                close(10)
                return
            else
                if (itime==1) call readmatgau(10,AOM1(:,:),1,"f14.8",6,5)
                if (itime==2) call readmatgau(10,AOM2(:,:),1,"f14.8",6,5)
            end if
        end do
    else if (wfntype==1) then !Unrestricted wavefunction
        nborb=idxHOMOb-nbasis
        allocate(AOM1b(nborb,nborb),AOM2b(nborb,nborb))
        do itime=1,2 !Two atoms
            if (itime==1) iatm=iatm1
            if (itime==2) iatm=iatm2
            do ispin=1,2
                if (ispin==1) then !alpha
                    write(c80tmp,"('Alpha part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
                    call loclabel(10,c80tmp,ifound)
                    if (itime==1) call readmatgau(10,AOM1(:,:),1,"f14.8",6,5)
                    if (itime==2) call readmatgau(10,AOM2(:,:),1,"f14.8",6,5)
                else if (ispin==2) then !beta
                    write(c80tmp,"('Beta part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
                    call loclabel(10,c80tmp,ifound)
                    if (itime==1) call readmatgau(10,AOM1b(:,:),1,"f14.8",6,5)
                    if (itime==2) call readmatgau(10,AOM2b(:,:),1,"f14.8",6,5)
                end if
            end do
        end do
    end if
    close(10)
    write(*,*) "AOMs of the two atoms have been successfully loaded"

else if (isel==2) then !Interbasin interaction analysis based on basin overlap matrix (BOM)
    write(*,*) "Input the path of the file containing BOM, e.g. C:\BOM.txt"
    write(*,*) "If press ENTER button directly, BOM.txt in current folder will be loaded"
    read(*,"(a)") c200tmp
    if (c200tmp==" ") then
        c200tmp="BOM.txt"
    end if
    write(*,*) "Input index of the two basins, e.g. 2,4"
    read(*,*) ibasin1,ibasin2
    allocate(BOM1(idxHOMO,idxHOMO),BOM2(idxHOMO,idxHOMO))
    open(10,file=c200tmp,status="old")
    if (wfntype==0) then
        do itime=1,2
            if (itime==1) ibasin=ibasin1
            if (itime==2) ibasin=ibasin2
            write(c80tmp,"('Orbital overlap matrix of basin',i6)") ibasin
            call loclabel(10,c80tmp,ifound)
            if (ifound==0) then
                write(*,"(' Error: Unable to find BOM of basin',i6)") ibasin
                write(*,*) "Press ENTER button to exit"
                read(*,*)
                close(10)
                return
            else
                if (itime==1) call readmatgau(10,BOM1(:,:),1,"f14.8",6,5)
                if (itime==2) call readmatgau(10,BOM2(:,:),1,"f14.8",6,5)
            end if
        end do
    else if (wfntype==1) then
        nborb=idxHOMOb-nbasis
        allocate(BOM1b(nborb,nborb),BOM2b(nborb,nborb))
        do itime=1,2 !Two basins
            if (itime==1) ibasin=ibasin1
            if (itime==2) ibasin=ibasin2
            do ispin=1,2
                if (ispin==1) then !alpha
                    write(c80tmp,"('Alpha part of orbital overlap matrix of basin',i6)") ibasin
                    call loclabel(10,c80tmp,ifound)
                    if (itime==1) call readmatgau(10,BOM1(:,:),1,"f14.8",6,5)
                    if (itime==2) call readmatgau(10,BOM2(:,:),1,"f14.8",6,5)
                else if (ispin==2) then !beta
                    write(c80tmp,"('Beta part of orbital overlap matrix of basin',i6)") ibasin
                    call loclabel(10,c80tmp,ifound)
                    if (itime==1) call readmatgau(10,BOM1b(:,:),1,"f14.8",6,5)
                    if (itime==2) call readmatgau(10,BOM2b(:,:),1,"f14.8",6,5)
                end if
            end do
        end do
    end if
    close(10)
    write(*,*) "BOMs of the two basins have been successfully loaded"
    
else if (isel==3) then !Interfragment interaction analysis based on the fragment overlap matrix (FOM) constructed from AOM
    write(*,*) "Input the path of the file containing AOM, e.g. C:\AOM.txt"
    write(*,*) "If press ENTER button directly, AOM.txt in current folder will be loaded"
    read(*,"(a)") c200tmp
    if (c200tmp==" ") c200tmp="AOM.txt"
    write(*,*) "Input index of the atoms in fragment 1, e.g. 2,4-10,16"
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,nfrag1atm,frag1atm)
    write(*,*) "Input index of the atoms in fragment 2, e.g. 1,3,11-15,17"
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,nfrag2atm,frag2atm)
    
    allocate(AOM1(idxHOMO,idxHOMO),FOM1(idxHOMO,idxHOMO),FOM2(idxHOMO,idxHOMO))
    FOM1=0
    FOM2=0
    open(10,file=c200tmp,status="old")
    if (wfntype==0) then !Closed-shell
        !Load AOM (into AOM1 as temporary slot) and construct FOM1
        do idx=1,nfrag1atm
            iatm=frag1atm(idx)
            write(c80tmp,"('Atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call loclabel(10,c80tmp,ifound)
            if (ifound==0) then
                write(*,"(' Error: Unable to find AOM of atom',i6)") iatm
                write(*,*) "Press ENTER button to exit"
                read(*,*)
                close(10)
                return
            else
                write(*,"(' Loading atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
                call readmatgau(10,AOM1(:,:),1,"f14.8",6,5)
                FOM1=FOM1+AOM1
            end if
        end do
        !Load AOM and construct FOM2
        do idx=1,nfrag2atm
            iatm=frag2atm(idx)
            write(c80tmp,"('Atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call loclabel(10,c80tmp,ifound)
            if (ifound==0) then
                write(*,"(' Error: Unable to find AOM of atom',i6)") iatm
                write(*,*) "Press ENTER button to exit"
                read(*,*)
                close(10)
                return
            else
                write(*,"(' Loading atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
                call readmatgau(10,AOM1(:,:),1,"f14.8",6,5)
                FOM2=FOM2+AOM1
            end if
        end do
        
    else if (wfntype==1) then !Unrestricted wavefunction
        nborb=idxHOMOb-nbasis
        allocate(AOM1b(nborb,nborb),FOM1b(nborb,nborb),FOM2b(nborb,nborb))
        FOM1b=0
        FOM2b=0
        !Load AOM (into AOM1 as temporary slot) and construct FOM1
        do idx=1,nfrag1atm
            iatm=frag1atm(idx)
            write(*,"(' Loading atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            write(c80tmp,"('Alpha part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call loclabel(10,c80tmp,ifound)
            call readmatgau(10,AOM1(:,:),1,"f14.8",6,5)
            FOM1=FOM1+AOM1
            write(c80tmp,"('Beta part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call loclabel(10,c80tmp,ifound)
            call readmatgau(10,AOM1b(:,:),1,"f14.8",6,5)
            FOM1b=FOM1b+AOM1b
        end do
        !Load AOM and construct FOM2
        do idx=1,nfrag2atm
            iatm=frag2atm(idx)
            write(*,"(' Loading atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            write(c80tmp,"('Alpha part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call loclabel(10,c80tmp,ifound)
            call readmatgau(10,AOM1(:,:),1,"f14.8",6,5)
            FOM2=FOM2+AOM1
            write(c80tmp,"('Beta part of atomic overlap matrix of',i6,'(',a2,')')") iatm,a(iatm)%name
            call loclabel(10,c80tmp,ifound)
            call readmatgau(10,AOM1b(:,:),1,"f14.8",6,5)
            FOM2b=FOM2b+AOM1b
        end do
    end if
    close(10)
    write(*,*) "Fragment overlap matrices have been successfully constructed!"

else if (isel==4) then !Interfragment interaction analysis based on FOM directly provided in FOM.txt
	inquire(file="FOM.txt",exist=alive)
    if (alive.eqv..false.) then
		write(*,"(a)") " Cannot find FOM.txt in current folder. &
		&Please directly input path of the file containing FOM of the two fragments, e.g. D:\Nea_makri\FOM.txt"
		do while(.true.)
			read(*,"(a)") c200tmp
			inquire(file=c200tmp,exist=alive)
			if (alive) exit
			write(*,*) "Cannot find the file, input again!"
		end do
    else
		c200tmp="FOM.txt"
    end if
    open(10,file=c200tmp,status="old")
    allocate(FOM1(idxHOMO,idxHOMO),FOM2(idxHOMO,idxHOMO))
    if (wfntype==0) then !Closed-shell
		write(*,*) "Loading FOM of fragment 1..."
		call readmatgau(10,FOM1(:,:),1,"f14.8",6,5)
		write(*,*) "Loading FOM of fragment 2..."
        read(10,*)
		call readmatgau(10,FOM2(:,:),1,"f14.8",6,5)
    else !Open-shell
        nborb=idxHOMOb-nbasis
        allocate(FOM1b(nborb,nborb),FOM2b(nborb,nborb))
		write(*,*) "Loading alpha FOM of fragment 1..."
		call readmatgau(10,FOM1(:,:),1,"f14.8",6,5)
		write(*,*) "Loading beta FOM of fragment 1..."
        read(10,*)
		call readmatgau(10,FOM1b(:,:),1,"f14.8",6,5)
		write(*,*) "Loading alpha FOM of fragment 2..."
        read(10,*)
		call readmatgau(10,FOM2(:,:),1,"f14.8",6,5)
		write(*,*) "Loading beta FOM of fragment 2..."
        read(10,*)
		call readmatgau(10,FOM2b(:,:),1,"f14.8",6,5)
    end if
    close(10)
    write(*,*) "Loading finished"
end if

write(*,*)
write(*,*) "Generating natural adaptive orbitals (NAdOs)..."

!Deal with closed-shell case or alpha part of unrestricted wavefunction
allocate(BODmat(idxHOMO,idxHOMO))
if (isel==1) then
	BODmat=matmul(AOM1,AOM2)+matmul(AOM2,AOM1)
else if (isel==2) then
	BODmat=matmul(BOM1,BOM2)+matmul(BOM2,BOM1)
else if (isel==3.or.isel==4) then
	BODmat=matmul(FOM1,FOM2)+matmul(FOM2,FOM1)
end if
allocate(eigvecmat(idxHOMO,idxHOMO),eigvalarr(idxHOMO))
call diagsymat(BODmat,eigvecmat,eigvalarr,istat)
CObasa(:,1:idxHOMO)=matmul(CObasa(:,1:idxHOMO),eigvecmat)
MOene(1:idxHOMO)=0 !Energy is meaningless
!Now the eigenvalues are from low to high. Reorder them so that rank from high to low
do iorb=1,int(idxHOMO/2)
    tmparr=CObasa(:,iorb)
    CObasa(:,iorb)=CObasa(:,idxHOMO-iorb+1)
    CObasa(:,idxHOMO-iorb+1)=tmparr
end do
call invarrr8(eigvalarr)
if (wfntype==0) then
    MOocc(1:idxHOMO)=eigvalarr(:)*2
    write(*,"(a,f10.5,a)") " Eigenvalues of NAdOs: (sum=",2*sum(eigvalarr)," )"
    write(*,"(7f10.5)") eigvalarr*2
else
    MOocc(1:idxHOMO)=eigvalarr(:)
    write(*,"(a,f10.5,a)") " Eigenvalues of alpha NAdOs: (sum=",sum(eigvalarr)," )"
    write(*,"(7f10.5)") eigvalarr
end if

if (wfntype==1) then !Deal with beta part of unrestricted wavefunction
    !nborb is the number of occupied beta orbitals
    allocate(BODmatb(nborb,nborb))
    if (isel==1) then
		BODmatb=matmul(AOM1b,AOM2b)+matmul(AOM2b,AOM1b)
    else if (isel==2) then
		BODmatb=matmul(BOM1b,BOM2b)+matmul(BOM2b,BOM1b)
    else if (isel==3.or.isel==4) then
		BODmatb=matmul(FOM1b,FOM2b)+matmul(FOM2b,FOM1b)
    end if
    deallocate(eigvecmat,eigvalarr)
    allocate(eigvecmat(nborb,nborb),eigvalarr(nborb))
    call diagsymat(BODmatb,eigvecmat,eigvalarr,istat)
    CObasb(:,1:nborb)=matmul(CObasb(:,1:nborb),eigvecmat)
    MOene(nbasis+1:idxHOMOb)=0 !Energy is meaningless
    do iorb=1,int((nborb)/2)
        tmparr=CObasb(:,iorb)
        CObasb(:,iorb)=CObasb(:,nborb-iorb+1)
        CObasb(:,nborb-iorb+1)=tmparr
    end do
    call invarrr8(eigvalarr)
    MOocc(nbasis+1:idxHOMOb)=eigvalarr(:)
    write(*,"(a,f10.5,a)") " Eigenvalues of beta NAdOs: (sum=",sum(eigvalarr)," )"
    write(*,"(7f10.5)") eigvalarr
end if

!Generate NAdO energies
write(*,*)
if (idoene==0) then
	write(*,*) "Note: Energies of NAdOs are not generated, they are simply set to zero"
else
	write(*,*) "Generating energies for NAdOs..."
    do iorb=1,idxHOMO
		mat11=matmul(matmul(transpose(CObasa(:,iorb:iorb)),FmatA),CObasa(:,iorb:iorb))
		MOene(iorb)=mat11(1,1)
    end do
    if (wfntype==1) then
		do iorbtmp=nbasis+1,idxHOMOb
			iorb=iorbtmp-nbasis
            mat11=matmul(matmul(transpose(CObasb(:,iorb:iorb)),FmatB),CObasb(:,iorb:iorb))
			MOene(iorbtmp)=mat11(1,1)
		end do
    end if
    write(*,*) "Done!"
end if

call outmwfn("NAdOs.mwfn",10,0)
write(*,"(/,a)") " All NAdO orbitals has been exported to NAdOs.mwfn in current folder"
write(*,"(a)") " In this file, original occupied orbitals have been replaced with NAdOs, occupation numbers correspond to eigenvalues, &
&while the unoccupied orbitals are still the original ones"
write(*,*)
write(*,*) "Do you want to load it now so that you can then study NAdOs? (y/n)"
read(*,*) selectyn
call dealloall(0)
if (selectyn=='y'.or.selectyn=='Y') then
    write(*,*) "Loading NAdOs.mwfn..."
    call readinfile("NAdOs.mwfn",1)
    write(*,"(a)") " Done! Now you can study NAdOs, such as using main function 0 to visualize orbitals"
else
    write(*,"(a)") " Reloading "//trim(filename)//" to recover initial status..."
    call readinfile(filename,1)
    write(*,*) "Loading finished!"
end if

end subroutine



!!----------- Calculate center, first/second moments, radius of gyration, and <r^2> of a function
subroutine funcmoment
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8 intval,moment1(3),moment2(3,3),moment2nuc(3,3),funcval(radpot*sphpot),beckeweigrid(radpot*sphpot),eigvecmat(3,3),eigval(3)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
integer :: iabs=0
character selectyn

ifunc=1
cenx=0
ceny=0
cenz=0
do while(.true.)
	write(*,*)
    if (iabs==0) write(*,"(a)") " -1 Toggle using absolute function value instead of original function value when calculating various quantities via option 1, current: No"
    if (iabs==1) write(*,"(a)") " -1 Toggle using absolute function value instead of original function value when calculating various quantities via option 1, current: Yes"
	write(*,*) "0 Return"
	write(*,*) "1 Calculate various quantities of the selected function"
	write(*,*) "2 Calculate center and integral of the selected function"
	write(*,"(a,i5)") " 3 Select the function to be studied, current:",ifunc
	write(*,"(a,3f11.5,' Ang')") " 4 Set the center for option 1, current:",cenx*b2a,ceny*b2a,cenz*b2a
	write(*,*) "5 Calculate center and integral of the absolute of the selected function"
	read(*,*) isel
	
	if (isel==0) then
		return
    else if (isel==-1) then
		if (iabs==0) then
			iabs=1
        else
			iabs=0
        end if
		cycle
	else if (isel==3) then
		call selfunc_interface(1,ifunc)
		cycle
	else if (isel==4) then
		write(*,*) "Input X,Y,Z of the center in Angstrom, e.g. 2.0,0,1.5"
		read(*,*) cenx,ceny,cenz
		cenx=cenx/b2a !To Bohr
		ceny=ceny/b2a
		cenz=cenz/b2a
		cycle
	end if

	write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
	call gen1cintgrid(gridatmorg,iradcut)

	call walltime(iwalltime1)

	intval=0
	moment1=0
	moment2=0
	realcenx=0
	realceny=0
	realcenz=0
    !Calculate center, integral, first and second moments
	do iatm=1,ncenter
		write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
		gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		!$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
		do i=1+iradcut*sphpot,radpot*sphpot
			funcval(i)=calcfuncall(ifunc,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
		end do
		!$OMP end parallel do
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
        
        !The following code doesn't work normally if using ifort 19 with -O2 under Linux, see http://sobereva.com/wfnbbs/viewtopic.php?id=669
		!do i=1+iradcut*sphpot,radpot*sphpot
		!	tmpval=funcval(i)*gridatmorg(i)%value*beckeweigrid(i)
		!	xtmp=gridatm(i)%x-cenx
		!	ytmp=gridatm(i)%y-ceny
		!	ztmp=gridatm(i)%z-cenz
  !          if (isel==5) tmpval=abs(tmpval)
		!	intval=intval+tmpval
		!	moment1(1)=moment1(1)+xtmp*tmpval
		!	moment1(2)=moment1(2)+ytmp*tmpval
		!	moment1(3)=moment1(3)+ztmp*tmpval
		!	moment2(1,1)=moment2(1,1)+xtmp*xtmp*tmpval
		!	moment2(2,2)=moment2(2,2)+ytmp*ytmp*tmpval
		!	moment2(3,3)=moment2(3,3)+ztmp*ztmp*tmpval
		!	moment2(1,2)=moment2(1,2)+xtmp*ytmp*tmpval
		!	moment2(2,3)=moment2(2,3)+ytmp*ztmp*tmpval
		!	moment2(1,3)=moment2(1,3)+xtmp*ztmp*tmpval
		!	realcenx=realcenx+gridatm(i)%x*tmpval
		!	realceny=realceny+gridatm(i)%y*tmpval
		!	realcenz=realcenz+gridatm(i)%z*tmpval
		!end do
        
        !Contributed by i.s.ger:
        !$OMP parallel do num_threads(nthreads) &
        !$OMP default(none) &
        !$OMP private(i, xtmp, ytmp, ztmp, tmpval) &
        !$OMP shared(iabs, iradcut, radpot, sphpot, isel, funcval, gridatm, gridatmorg, beckeweigrid, cenx, ceny, cenz) &
        !$OMP reduction(+:intval, moment1, moment2, realcenx, realceny, realcenz)
        do i = 1 + iradcut*sphpot, radpot*sphpot
            tmpval = funcval(i)*gridatmorg(i)%value*beckeweigrid(i)
            xtmp = gridatm(i)%x - cenx
            ytmp = gridatm(i)%y - ceny
            ztmp = gridatm(i)%z - cenz
            if (isel == 5 .or. iabs==1) tmpval = abs(tmpval)
            intval = intval + tmpval
            moment1(1) = moment1(1) + xtmp*tmpval
            moment1(2) = moment1(2) + ytmp*tmpval
            moment1(3) = moment1(3) + ztmp*tmpval
            moment2(1, 1) = moment2(1, 1) + xtmp*xtmp*tmpval
            moment2(2, 2) = moment2(2, 2) + ytmp*ytmp*tmpval
            moment2(3, 3) = moment2(3, 3) + ztmp*ztmp*tmpval
            moment2(1, 2) = moment2(1, 2) + xtmp*ytmp*tmpval
            moment2(2, 3) = moment2(2, 3) + ytmp*ztmp*tmpval
            moment2(1, 3) = moment2(1, 3) + xtmp*ztmp*tmpval
            realcenx = realcenx + gridatm(i)%x*tmpval
            realceny = realceny + gridatm(i)%y*tmpval
            realcenz = realcenz + gridatm(i)%z*tmpval
        end do
        !$OMP end parallel do
	end do
    
	call walltime(iwalltime2)
	write(*,"(' Calculation took up wall clock time',i10,' s',/)") iwalltime2-iwalltime1
	
    !Print result
	if (isel==1) then !Various quantities
		moment2(3,1)=moment2(1,3)
		moment2(2,1)=moment2(1,2)
		moment2(3,2)=moment2(2,3)
		write(*,*) "Note: Unless otherwise specified, all data shown below are in a.u."
		write(*,"(/,' Integral over whole space:',1PE16.8,/)") intval
		write(*,"(' The first moment:')")
		write(*,"(' X= ',1PE16.8,'   Y= ',1PE16.8,'   Z= ',1PE16.8)") moment1
		write(*,"(' Norm= ',1PE16.8,/)") sum(moment1**2)
		write(*,"(' The second moment:')")
		write(*,"(' XX=',1PE16.8,'   XY=',1PE16.8,'   XZ=',1PE16.8)") moment2(1,:)
		write(*,"(' YX=',1PE16.8,'   YY=',1PE16.8,'   YZ=',1PE16.8)") moment2(2,:)
		write(*,"(' ZX=',1PE16.8,'   ZY=',1PE16.8,'   ZZ=',1PE16.8)") moment2(3,:)

		call diagmat(moment2,eigvecmat,eigval,300,1D-10)
		call sort(eigval)
		write(*,"(a,3(1PE16.8))") ' Eigenvalues:',eigval
        write(*,"(a,1PE16.8)") " Sum of eigenvalues (trace of the second moment tensor):",sum(eigval)
		write(*,"(' Anisotropy:',1PE16.8,/)") eigval(3)-(eigval(1)+eigval(2))/2D0
        rgyr=dsqrt((moment2(1,1)+moment2(2,2)+moment2(3,3))/intval)
		write(*,"(' Radius of gyration:',1PE16.8, ' Bohr ',1PE16.8, ' Angstrom')") rgyr,rgyr*b2a
        write(*,"(/,a,f16.6)") " Spatial extent of the function <r^2>:",sum(eigval)

        !If the selected function is electron density, also print electric moments
		if (ifunc==1) then
			moment2nuc=0
			do iatm=1,ncenter
				xtmp=a(iatm)%x-cenx
				ytmp=a(iatm)%y-ceny
				ztmp=a(iatm)%z-cenz
				tmpval=a(iatm)%charge
				moment2nuc(1,1)=moment2nuc(1,1)+xtmp*xtmp*tmpval
				moment2nuc(2,2)=moment2nuc(2,2)+ytmp*ytmp*tmpval
				moment2nuc(3,3)=moment2nuc(3,3)+ztmp*ztmp*tmpval
				moment2nuc(1,2)=moment2nuc(1,2)+xtmp*ytmp*tmpval
				moment2nuc(2,3)=moment2nuc(2,3)+ytmp*ztmp*tmpval
				moment2nuc(1,3)=moment2nuc(1,3)+xtmp*ztmp*tmpval
			end do
			moment2nuc(3,1)=moment2nuc(1,3)
			moment2nuc(2,1)=moment2nuc(1,2)
			moment2nuc(3,2)=moment2nuc(2,3)
			write(*,*)
			write(*,"(' The quadrupole moment of nuclear charges:')")
			write(*,"(' XX=',f16.8,'   XY=',f16.8,'   XZ=',f16.8)") moment2nuc(1,:)
			write(*,"(' YX=',f16.8,'   YY=',f16.8,'   YZ=',f16.8)") moment2nuc(2,:)
			write(*,"(' ZX=',f16.8,'   ZY=',f16.8,'   ZZ=',f16.8)") moment2nuc(3,:)
			write(*,*)
			write(*,"(' The quadrupole moment of the system:')")
			write(*,"(' XX=',f16.8,'   XY=',f16.8,'   XZ=',f16.8)") moment2nuc(1,:)-moment2(1,:)
			write(*,"(' YX=',f16.8,'   YY=',f16.8,'   YZ=',f16.8)") moment2nuc(2,:)-moment2(2,:)
			write(*,"(' ZX=',f16.8,'   ZY=',f16.8,'   ZZ=',f16.8)") moment2nuc(3,:)-moment2(3,:)
		end if
		
	else if (isel==2.or.isel==5) then !Function center
		realcenx=realcenx/intval
		realceny=realceny/intval
		realcenz=realcenz/intval
        if (isel==2) then
			write(*,"(' Integral of the function:',1PE16.8,' a.u.')") intval
			write(*,"(/,' Center of the function:')")
        else
			write(*,"(' Integral of the absolute of the function:',1PE16.8,' a.u.')") intval
			write(*,"(/,' Center of the absolute of the function:')")
        end if
		write(*,"(' X=',f16.8,' Y=',f16.8,' Z=',f16.8,' Angstrom',/)") realcenx*b2a,realceny*b2a,realcenz*b2a
		write(*,*) "Use this center for subsequent calculations? (y/n)"
		read(*,*) selectyn
		if (selectyn=='y') then
			cenx=realcenx
			ceny=realceny
			cenz=realcenz
		end if
	end if
end do
end subroutine




!!----------- Calculate spatial delocalization index (SDI)
subroutine SDI
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
character c2000tmp*2000
real*8 :: expfac=2D0
integer,allocatable :: tmparr(:)
real*8,allocatable :: SDIarr(:),densarr(:,:)
real*8 funcval(radpot*sphpot),beckeweigrid(radpot*sphpot),wfnval(nmo)
type(content) gridatmorg(radpot*sphpot),gridatm(radpot*sphpot)

if (ifPBC/=0) then
    write(*,*) "Error: PBC has not been supported by this function yet!"
    return
end if

do while(.true.)
	write(*,*)
	call menutitle("Calculate spatial delocalization index (SDI)",10,1)
	write(*,"(a,f8.4)") " -1 Set exponent factor, current:",expfac
	write(*,*) "0 Return"
	write(*,*) "1 Calcluate SDI for a real space function"
	write(*,*) "2 Calculate SDI for density of orbital wavefunctions"
	if (allocated(cubmat)) write(*,*) "3 Calculate SDI based on the grid data in memory"
	read(*,*) icalctype
	if (icalctype==0) then
		return
	else if (icalctype==-1) then
		write(*,*) "Input exponent factor, e.g. 1.5"
        read(*,*) expfac
	else
		exit
	end if
end do

if (icalctype==1) then !Calculate for real space function
    call selfunc_interface(1,ifunc)
	write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
	call gen1cintgrid(gridatmorg,iradcut)
	call walltime(iwalltime1)
	valintabs=0
	valint2=0
	do iatm=1,ncenter
		write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
		gridatm%x=gridatmorg%x+a(iatm)%x
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
		!$OMP parallel do shared(funcval) private(ipt) num_threads(nthreads)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			funcval(ipt)=calcfuncall(ifunc,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
		end do
		!$OMP end parallel do
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			valintabs=valintabs+abs(funcval(ipt))*gridatmorg(ipt)%value*beckeweigrid(ipt)
			valint2=valint2+abs(funcval(ipt))**expfac *gridatmorg(ipt)%value*beckeweigrid(ipt)
		end do
	end do
    valint2=valint2/valintabs**expfac
	call walltime(iwalltime2)
	write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
    SDIval=1D0/dsqrt(valint2)
    write(*,"(/,' Spatial delocalization index is',f12.6)") SDIval

else if (icalctype==2) then !Calculate based on orbital density
    write(*,*) "Input indices of the orbitals to calculate SDI, e.g. 2,3,7-10"
    read(*,"(a)") c2000tmp
    call str2arr(c2000tmp,ntmp)
    allocate(tmparr(ntmp),SDIarr(ntmp))
    call str2arr(c2000tmp,ntmp,tmparr)
    if (iautointgrid==1) then
        nradpotold=radpot
        nsphpotold=sphpot
        radcutold=radcut
        radpot=30
        sphpot=302
        radcut=15
    end if
    allocate(densarr(nmo,radpot*sphpot))
    write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
    call gen1cintgrid(gridatmorg,iradcut)
    call walltime(iwalltime1)
    densarr=0
    do iatm=1,ncenter
	    write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
	    gridatm%x=gridatmorg%x+a(iatm)%x
	    gridatm%y=gridatmorg%y+a(iatm)%y
	    gridatm%z=gridatmorg%z+a(iatm)%z
	    !$OMP parallel do shared(densarr) private(ipt,wfnval) num_threads(nthreads)
	    do ipt=1+iradcut*sphpot,radpot*sphpot
            call orbderv(1,1,nmo,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,wfnval(:))
		    densarr(:,ipt)=wfnval(:)**2
	    end do
	    !$OMP end parallel do
	    call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	    do ipt=1+iradcut*sphpot,radpot*sphpot
		    SDIarr(:)=SDIarr(:)+densarr(:,ipt)**expfac *gridatmorg(ipt)%value*beckeweigrid(ipt)
	    end do
    end do
    SDIarr(:)=1D0/dsqrt(SDIarr(:))
    call walltime(iwalltime2)
    write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
    write(*,*)
    do itmp=1,ntmp
        iorb=tmparr(itmp)
        write(*,"(' SDI of orbital',i6,':',f10.4)") iorb,SDIarr(iorb)
    end do
    write(*,*)
    write(*,*) "Note: The unit of SDI of orbital density is a.u."
    if (iautointgrid==1) then
	    radpot=nradpotold
	    sphpot=nsphpotold
        radcut=radcutold
    end if
    
else if (icalctype==3) then !Calculate based on grid data
    call calc_dvol(dvol)
    sumall=sum(abs(cubmat(:,:,:)))
    valint=sumall*dvol
    write(*,"(' Differential element:',f15.10,' Bohr^3')") dvol
	write(*,"(' Sum of all values of grid data:',1PE16.8)") sumall
	write(*,"(' Integral of absolute value of grid data:',1PE16.8)") valint
    if (allocated(cubmattmp)) deallocate(cubmattmp)
    allocate(cubmattmp(nx,ny,nz))
    cubmattmp=cubmat/valint
    SDIval=1/dsqrt(sum(abs(cubmattmp(:,:,:))**expfac)*dvol)
    write(*,"(/,' Spatial delocalization index is',f12.6)") SDIval
    deallocate(cubmattmp)
end if

end subroutine