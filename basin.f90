!!------ Forming basin for specific real space function, and integrate real space functions in the basins
!!------ The method is adapted from J. Phys.: Condens. Matter 21 (2009) 084204
!Grid must be ortho
subroutine basinana
use defvar
use basinintmod
use util
use GUI
implicit real*8 (a-h,o-z)
integer walltime1,walltime2
integer :: igridmethod=3
integer,allocatable :: attconvold(:),realattconv(:),usercluslist(:),tmparr(:)
character basinfilename*200,selectyn,c80tmp*80,ctmp1*20,ctmp2*20,ctmp3*20,ctmp4*20,c2000tmp*2000,c200tmp*200
real*8 :: threslowvalatt=1D-5
real*8,allocatable :: tmparrf(:)
logical alive1,alive2
ishowattlab=1
ishowatmlab=0
ishowatt=1
!Don't show topology analysis results to avoid confusion
ishowCPlab=0
ishowpathlab=1
ishow3n3=0
ishow3n1=0
ishow3p1=0
ishow3p3=0
!When enter this module, we reset the whole state, which is signified by numatt. Because the user may calculate grid data at external functions,
!and hence messed up grid setting (orgx/y/z...), so all of the functions in this module that rely on grid setting will be totally wrong
numatt=0

do while(.true.)
	write(*,*)
	write(*,*) "                 ============= Basin analysis ============="
	write(*,*) "-10 Return to main menu"
	if (numatt>0) write(*,*) "-6 Set parameter for attractor clustering or manually perform clustering"
	if (numatt==0) write(*,*) "-6 Set parameter for attractor clustering"
    if (numatt>0) write(*,*) "-45 Export attractor information and cube file of present grid data"
	if (numatt>0) write(*,*) "-5 Export basins as cube file"
	if (numatt>0) write(*,*) "-4 Export attractors as pdb/pqr/txt/gjf file"
	if (numatt>0) write(*,*) "-3 Show information of attractors"
	if (numatt>0) write(*,*) "-2 Measure distances, angles and dihedral angles between attractors or atoms"
    write(*,*) "-1 Select the method for generating basins"
	if (numatt>0) write(*,*) " 0 Visualize attractors and basins"
	if (numatt>0) then
		write(*,*) " 1 Regenerate basins and relocate attractors"
	else
		write(*,*) " 1 Generate basins and locate attractors"
	end if
	if (numatt>0) then
		write(*,*) " 2 Integrate real space functions in the basins"
		write(*,*) " 3 Calculate electric multipole moments and <r^2> for basins"
		write(*,*) " 4 Calculate localization index (LI) and delocalization index (DI) for basins"
		write(*,*) " 5 Output orbital overlap matrix in basins to BOM.txt in current folder"
! 		write(*,*) " 100 Integrate real space functions in the basins with multi-level refinement"
		if (ifuncbasin==1) then
		    write(*,*) " 6 Output orbital overlap matrix in atoms to AOM.txt in current folder"
			write(*,*) " 7 Integrate real space functions in AIM basins with mixed type of grids"
			write(*,*) " 8 Calculate electric multipole moments in AIM basins with mixed type of grids"
			write(*,*) " 9 Obtain atomic contribution to population of external basins"
        else if (ifuncbasin==9) then
            write(*,"(a)") " 10 Calculate high ELF localization domain population and volume (HELP, HELV)"
		end if
		write(*,*) "11 Calculate orbital compositions contributed by various basins"
        if (ifuncbasin==9) write(*,*) "12 Assign ELF basin labels"
	end if
	read(*,*) isel
	if (numatt==0.and.(isel==-5.or.isel==-4.or.isel==-3.or.isel==-2.or.isel==0.or.isel==2.or.isel==3.or.isel==4.or.isel==5)) then
		write(*,*) "Error: You should use option 1 to generate basins first!"
		cycle
	end if
    
    if ((isel>=1.and.isel<=3).or.(isel>=7.and.isel<=10)) call delvirorb(0) !For functions do not rely on virtual orbitals, remove high-lying ones to save cost

	if (isel==-10) then
		idrawbasinidx=-10 !Don't display interbasin in any other GUI
		return
		
	else if (isel==-6) then
		do while(.true.)
			write(*,*) "0 Return"
			write(*,"(a,f14.5,a)") " 1 Set relative value difference criterion, current:",valcritclus*100,"%"
			write(*,"(a,i3)") " 2 Set the multiplier for distance criterion, current:",mergeattdist
			if (numatt>=2) write(*,*) "3 Cluster specified attractors"
			read(*,*) isel2
			if (isel2==0) then
				exit
			else if (isel2==1) then
				write(*,"(a)") " Input a value in percentage, e.g. 0.5"
				write(*,"(a)") " Note: Inputting X means if the value difference between two attractors is less than X%, &
				&then they will be regarded as degenerate and may be clustered according to distance criterion"
				read(*,*) valcritclus
				valcritclus=valcritclus/100D0
			else if (isel2==2) then
				write(*,"(a)") " Input a value, e.g. 6"
				write(*,"(a)") " Note: Inputting P means for any two attractors that have relative value difference less than the one set by option 1, &
				&if their interval is less than P*sqrt(dx^2+dy^2+dz^2), &
				&where dx,dy,dz are grid spacing in X,Y,Z,then they will be clustered together. If you want to nullify the clustering step after generating &
				&basins, simply set this value to 0"
				read(*,*) mergeattdist
			else if (isel2==3) then
				if (numatt<2) then
					write(*,*) "Error: At least two attractors must be presented!"
					cycle
				end if
				write(*,*) "Input the attractors you want to cluster, e.g. 4,5,8,9,22,21"
				read(*,"(a)") c2000tmp
				call str2arr(c2000tmp,ntobeclus) !Find how many terms
				if (allocated(attconvold)) deallocate(attconvold,realattconv,usercluslist)
				allocate(attconvold(-2:numatt),realattconv(-2:numrealatt),usercluslist(ntobeclus))
				realattconv(-2)=-2
				realattconv(-1)=-1
				realattconv(0)=0
				call str2arr(c2000tmp,ntobeclus,usercluslist)
				call sort(usercluslist)
				attconvold=attconv
				idesrealatt=usercluslist(1)
				do itmp=2,ntobeclus
					irealatt=usercluslist(itmp)
					do jtmp=1,nrealatthas(irealatt)
						iatt=realatttable(irealatt,jtmp)
						attconv(iatt)=idesrealatt
					end do
				end do
				call clusdegenatt(1)
				realattconv(1)=1
				do itmp=1,numatt !Build conversion relationship between previous real attractors and new real attractors
					if (itmp/=1) then
						if (attconvold(itmp)/=attconvold(itmp-1)) realattconv(attconvold(itmp))=attconv(itmp)
					end if
				end do
				do iz=izlow,izup
					do iy=iylow,iyup
						do ix=ixlow,ixup
							gridbas(ix,iy,iz)=realattconv(gridbas(ix,iy,iz))
						end do
					end do
				end do
				call detectinterbasgrd(6) !Detect interbasin grids
				write(*,*) "Done!"
				write(*,*)
			end if
		end do
	
	else if (isel==-45) then !Export attractors+basins and cube file of present grid data
		call export_basinana_info
		
	else if (isel==-5) then !Output cube file for basins for visualizing interbasin surface
		write(*,"(a)") " Inputting e.g. 4,6-9,15,17 will output these basins to individual basin[index].cub files to current folder"
		write(*,"(a)") " Inputting ""a"" will output basin.cub in current folder, whose value corresponds to basin index"
		write(*,"(a)") "   If you want to export grid data of function value in the region of specific basins to basinsel.cub in current folder, input ""c"""
        write(*,"(a)") "   If you want to use VMD to plot ELF isosurface map colored by basin type &
        &(monosynaptic or disynaptic), input ""b"". See Section 4.17.10 of manual for detail"
		read(*,"(a)") c2000tmp
        if (index(c2000tmp,'a')/=0) then
			write(*,*) "Exporting..."
		    if (allocated(cubmattmp)) deallocate(cubmattmp)
            nxtmp=ixup-ixlow+1
            nytmp=iyup-iylow+1
            nztmp=izup-izlow+1
		    allocate(cubmattmp(nxtmp,nytmp,nztmp))
            cubmattmp=gridbas(ixlow:ixup,iylow:iyup,izlow:izup)
            open(10,file="basin.cub",status="replace")
            call outcube(cubmattmp,nxtmp,nytmp,nztmp,orgx+gridv1(1),orgy+gridv2(2),orgz+gridv3(3),gridv1,gridv2,gridv3,10)
            close(10)
		    deallocate(cubmattmp)
			write(*,"(a)") " Done! basin.cub has been outputted to current folder. The grid values correspond to basin index"
        else if (index(c2000tmp,'b')/=0) then
            if (allocated(tmparr)) deallocate(tmparr)
			write(*,*) "Input indices of monosynaptic basins, e.g. 1,4-9,12,15"
            write(*,"(a)") " Hint: You can visually identify indices of monosynaptic basins in option 0 or with help of option 12"
            write(*,*) "To return, input ""q"""
            read(*,"(a)") c2000tmp
            if (c2000tmp=="q") cycle
            call str2arr(c2000tmp,ntmp)
            allocate(tmparr(ntmp))
            call str2arr(c2000tmp,ntmp,tmparr)
		    if (allocated(cubmattmp)) deallocate(cubmattmp)
		    allocate(cubmattmp(nx,ny,nz))
            cubmattmp=0
            do itmp=1,ntmp
				ibas=tmparr(itmp)
                where(gridbas==ibas) cubmattmp=-1
            end do
            deallocate(tmparr)
			write(*,*) "Input indices of disynaptic basins, e.g. 2,3,10-11"
            write(*,"(a)") " Hint: You can visually identify indices of monosynaptic basins in option 0 or with help of option 12"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,ntmp)
            allocate(tmparr(ntmp))
            call str2arr(c2000tmp,ntmp,tmparr)
            do itmp=1,ntmp
				ibas=tmparr(itmp)
                where(gridbas==ibas) cubmattmp=1
            end do
			write(*,"(a)") " Exporting cube file containing basin type index to basinsyn.cub in current folder..."
			open(10,file="basinsyn.cub",status="replace")
            call outcube(cubmattmp,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
			close(10)
            write(*,"(a)") " Finished! In this file the values within the monosynaptic and disynaptic basin regions are -1 and 1, respectively, all other regions have value of 0. All other regions have value of 0"
            deallocate(cubmattmp)
			write(*,*) "Exporting the real space function used to generate basins to basinfunc.cub..."
			open(10,file="basinfunc.cub",status="replace")
            call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
			close(10)
            write(*,*) "Finished! This file records value of the function used to generate basins"
            write(*,"(a)") " Now you can use ""examples\scripts\basinsyn.tcl"" in Multiwfn package to draw basin type colored isosurface map in VMD"
        else if (index(c2000tmp,'c')/=0) then
            if (allocated(tmparr)) deallocate(tmparr)
			write(*,*) "Input indices of the basins of interest, e.g. 2,3,10-11"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,ntmp)
            allocate(tmparr(ntmp))
            call str2arr(c2000tmp,ntmp,tmparr)
            write(*,*) "Set the value of the grids that do not belong to the selected basins, e.g. 0"
            read(*,*) tmpval
            allocate(cubmattmp(nx,ny,nz))
            cubmattmp=tmpval
			do iz=izlow,izup
				do iy=iylow,iyup
					do ix=ixlow,ixup
						if (any(tmparr==gridbas(ix,iy,iz))) cubmattmp(ix,iy,iz)=cubmat(ix,iy,iz)
					end do
				end do
			end do
			write(*,"(a)") " Exporting cube file of selected basins to basinsel.cub in current folder..."
			open(10,file="basinsel.cub",status="replace")
            call outcube(cubmattmp,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
			close(10)
            write(*,*) "Finished!"
            deallocate(cubmattmp)
		else
            call str2arr(c2000tmp,ntmp)
            if (allocated(tmparr)) deallocate(tmparr)
            allocate(tmparr(ntmp))
            call str2arr(c2000tmp,ntmp,tmparr)
            if (any(tmparr<=0).or.any(tmparr>numatt)) then
                write(*,*) "Error: One or more basin indices exceeded valid range!"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                cycle
            end if
            write(*,"(a)") " Select how to output cube files. The value of the grids inside and outside the region of selected basins will be 1.0 and 0.0, respectively"
            write(*,*) "1 Output all basin grids"
            write(*,*) "2 Output boundary basin grids"
            write(*,*) "3 Output all basin grids where electron density > 0.001 a.u."
            write(*,"(a)") " If you want to export the current real space function in respective basin region &
            &(function value outside the corresponding basin is set to zero), input 0"
		    read(*,*) igrid
		    if (allocated(cubmattmp)) deallocate(cubmattmp)
		    allocate(cubmattmp(nx,ny,nz))
		    do idx=1,ntmp
                ibas=tmparr(idx)
			    write(basinfilename,"('basin',i4.4,'.cub')") ibas
			    write(*,"(' Outputting basin',i8,' as ',a)") ibas,trim(basinfilename)
                if (igrid==0) then
					cubmattmp=cubmat
					where(gridbas/=ibas) cubmattmp=0
                else
					cubmattmp=gridbas(:,:,:)
					where(cubmattmp/=ibas)
						cubmattmp=0
					elsewhere
						cubmattmp=1
					end where
					if (igrid==1) then
						continue
					else if (igrid==2) then
						where (interbasgrid.eqv..false.) cubmattmp=0
					else if (igrid==3) then
						if (ifuncbasin==1) then !The cubmat already records electron density
							where(cubmat<0.001D0) cubmattmp=0
						else
							call saverhocub
							where(rhocub<0.001D0) cubmattmp=0
						end if
					end if
                end if
			    open(10,file=basinfilename,status="replace")
                call outcube(cubmattmp,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
			    close(10)
		    end do
		    deallocate(cubmattmp)
			write(*,"(a)") " Done! Cube files for the selected basins have been outputted to current folder"
            if (igrid==2) then
			    write(*,"(a)") " You can plot isosurface with isovalue=0.5 to visualize the basin boundary surface"
            else if (igrid==1.or.igrid==3) then
			    write(*,"(a)") " Values 1 and 0 in these files indicate that the corresponding point belongs and not belongs to the basin, & 
			    &respectively. You can plot isosurface with isovalue=0.5 to visualize basin region"
            end if
		end if
		
	else if (isel==-4) then
        write(*,*) "0 Return"
        write(*,*) "1 Export coordinates of all attractors as attractors.pdb"
        write(*,"(a)") " 2 Export coordinates and function values of all attractors as attractors.pqr"
        write(*,"(a)") " 3 Export coordinates and function values of all attractors as attractors.txt"
        write(*,"(a)") " 4 Export atoms and all attractors as attractors.gjf"
        read(*,*) isel2
        if (isel2==0) then
            continue
        else if (isel2==1) then
		    open(10,file="attractors.pdb",status="replace")
			if (ifPBC>0) then
				call getcellabc(asize,bsize,csize,alpha,beta,gamma)
				write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
			end if
		    do iatt=1,numatt
			    write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",iatt,' '//"C "//' ',"ATT",'A',attconv(iatt),attxyz(:,iatt)*b2a,1.0,0.0,"C "
		    end do
		    close(10)
		    write(*,*) "Done, all attractors have been exported to attractors.pdb in current folder"
		    write(*,"(a)") " Note: The residue indices in the pdb file correspond to the attractor indices after clustering, &
            &while the atomic indices correspond to the raw attractor indices"
        else if (isel2==2) then
		    open(10,file="attractors.pqr",status="replace")
			if (ifPBC>0) then
				call getcellabc(asize,bsize,csize,alpha,beta,gamma)
				write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
			end if
		    do iatt=1,numatt
			    write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,E18.8E3,f5.2,1x,a2)") "HETATM",iatt,' '//"C "//' ',"ATT",'A',attconv(iatt),attxyz(:,iatt)*b2a,attval(iatt),1.0,"C "
		    end do
		    close(10)
		    write(*,*) "Done, all attractors have been exported to attractors.pqr in current folder"
		    write(*,"(a)") " Note: The residue indices in the pqr file correspond to the attractor indices after clustering, &
            &while the atomic indices correspond to the raw attractor indices. The atomic charge column corresponds to function value"
        else if (isel2==3) then
		    open(10,file="attractors.txt",status="replace")
		    do iatt=1,numatt
                write(10,"(3f12.6,1PE18.8E3)") attxyz(:,iatt),attval(iatt)
            end do
		    close(10)
		    write(*,*) "Done, all attractors have been exported to attractors.txt in current folder"
            write(*,"(a)") " In this file, the first three columns correspond to X,Y,Z coordinate in Bohr, the final column is function value"
        else if (isel2==4) then
		    open(10,file="attractors.gjf",status="replace")
            write(10,"(a,/,/,a,/,/,'0 1')") "# B3LYP/6-31G* geom=connectivity","Generated by basin analysis module of Multiwfn"
            do iatm=1,ncenter
                write(10,"(a,3f12.6,' H')") a(iatm)%name,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a
            end do
		    do iatt=1,numatt
                write(10,"('Bq',3f12.6,' L')") attxyz(:,iatt)*b2a
            end do
            call genconnmat(1,0)
            write(10,*)
            do iatm=1,ncenter
                write(10,"(i6)",advance="no") iatm
                do jatm=iatm+1,ncenter
                    if (connmat(iatm,jatm)==0) cycle
                    write(10,"(i6,f4.1)",advance="no") jatm,dfloat(connmat(iatm,jatm))
                end do
                write(10,*)
            end do
            do itmp=1,numatt
                write(10,"(i6)") itmp+ncenter
            end do
            write(10,*);write(10,*)
		    close(10)
		    write(*,"(a)") " Done, all atoms and attractors have been exported to attractors.gjf in current folder. &
            &You can directly use GaussView to visualize. It is suggested to choose &
            &""File""-""Preference""-""View""-""Display Format""-""Molecule"", and set low layer as ""Tube"""
        end if
		
	else if (isel==-3) then
		write(*,*) "  Attractor       X,Y,Z coordinate (Angstrom)            Value"
		do irealatt=1,numrealatt
			write(*,"(i8,3f14.8,1E18.9)") irealatt,realattxyz(:,irealatt)*b2a,realattval(irealatt)
		end do
		do irealatt=1,numrealatt
			if (nrealatthas(irealatt)/=1) then
				write(*,*)
				write(*,"(' The members of degenerate attractor',i6,':')") irealatt
				do itmp=1,nrealatthas(irealatt)
					iatt=realatttable(irealatt,itmp)
					write(*,"(i8,'  XYZ:',3f13.7,'  Value:',1E18.9)") iatt,attxyz(:,iatt)*b2a,attval(iatt)
				end do
			end if
		end do
        if (numrealatt>1) then
            write(*,"(/,a)") " Do you also want to print the attractor information after sorting according to values? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y'.or.selectyn=='Y') then
                if (allocated(tmparr)) deallocate(tmparr)
                allocate(tmparr(numrealatt))
                forall(i=1:numrealatt) tmparr(i)=i
                allocate(tmparrf(numrealatt))
                tmparrf=realattval
                call sortr8(tmparrf,"val",tmparr)
		        write(*,*) "  Attractor       X,Y,Z coordinate (Angstrom)            Value"
		        do idx=1,numrealatt
                    irealatt=tmparr(idx)
			        write(*,"(i8,3f14.8,1E18.9)") irealatt,realattxyz(:,irealatt)*b2a,realattval(irealatt)
                end do
                deallocate(tmparrf,tmparr)
            end if
        end if
		
	else if (isel==-2) then
		write(*,*) "q = Quit. Selection method:"
		write(*,*) "a? = Atom ?"
		write(*,*) "c? = Attractor ? (If the attractor is degenerate, average coordinate is used)"
		write(*,*) "d? = Attractor ? (Only for degenerate attractors, cycle each of its members)"
		write(*,*) "For example:"
		write(*,*) """a1 c3"" returns the distance between atom1 and att.3"
		write(*,*) """a4 a2"" returns the distance between atom4 and atom2"
		write(*,*) """c6 a2 a5"" returns the angle of att.6-atom2-atom5"
		write(*,*) """c2 c4 a3 c7"" returns the dihedral angle of att.2-att.4-atom3-att.7"
		write(*,*) """d3 a1"" returns the distance between members of att.3 and atom1"
		write(*,"(a)") " ""d3 a1 c2"" returns the the angle of (members of att.3)--atom1--att.2, &
		&meanwhile outputs the vertical distance from (members of att.3) to the line linking atom1 and att.2"
		do while(.true.)
			read(*,"(a)") c80tmp
			c80tmp=adjustl(c80tmp)
			imeasure=0
			do ichar=1,len_trim(c80tmp) !imeasure=1/2/3: measure distance,angle,dihedral
				if (c80tmp(ichar:ichar)==','.or.c80tmp(ichar:ichar)==' ') imeasure=imeasure+1
			end do
			nelement=0 !Validate "d" selection
			do ichar=1,len_trim(c80tmp)
				if (c80tmp(ichar:ichar)=='d') nelement=nelement+1
			end do
			if (nelement/=0) then
				if (c80tmp(1:1)/='d'.or.nelement>1) then
					write(*,*) "Error: ""d"" type of selection can only occur once and must be the first term"
					cycle
				else if (imeasure==3) then
					write(*,*) "Error: Dihedral angle calculation does not support ""d"" type of selection"
					cycle					
				end if
			end if
							
			if (c80tmp(1:1)=='q') then
				exit
			else if (imeasure==1.or.imeasure==2.or.imeasure==3) then
				if (imeasure==1) read(c80tmp,*) ctmp1,ctmp2 !Read two terms
				if (imeasure==2) read(c80tmp,*) ctmp1,ctmp2,ctmp3 !Read three terms
				if (imeasure==3) read(c80tmp,*) ctmp1,ctmp2,ctmp3,ctmp4 !Read four terms
				
				if (ctmp1(1:1)=='a') then
					read(ctmp1(2:),*) iatm
					tmpx1=a(iatm)%x
					tmpy1=a(iatm)%y
					tmpz1=a(iatm)%z
				else if (ctmp1(1:1)=='c') then
					read(ctmp1(2:),*) irealatt
					tmpx1=realattxyz(1,irealatt)
					tmpy1=realattxyz(2,irealatt)
					tmpz1=realattxyz(3,irealatt)
				else if (ctmp1(1:1)=='d') then
					read(ctmp1(2:),*) irealattde
				end if
				if (ctmp2(1:1)=='a') then
					read(ctmp2(2:),*) iatm
					tmpx2=a(iatm)%x
					tmpy2=a(iatm)%y
					tmpz2=a(iatm)%z
				else if (ctmp2(1:1)=='c') then
					read(ctmp2(2:),*) irealatt
					tmpx2=realattxyz(1,irealatt)
					tmpy2=realattxyz(2,irealatt)
					tmpz2=realattxyz(3,irealatt)
				end if
				
				if (imeasure==1) then
					if (ctmp1(1:1)/='d') then
						write(*,"(' The distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") &
						dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2),dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2)*b2a
					else if (ctmp1(1:1)=='d') then
						distmin=9999999
						distmax=-1
						distavg=0
						do itmp=1,nrealatthas(irealattde)
							iatt=realatttable(irealattde,itmp)
							tmpx1=attxyz(1,iatt)
							tmpy1=attxyz(2,iatt)
							tmpz1=attxyz(3,iatt)
							distnow=dsqrt((tmpx1-tmpx2)**2+(tmpy1-tmpy2)**2+(tmpz1-tmpz2)**2)
							distavg=distavg+distnow
							if (distnow<distmin) distmin=distnow
							if (distnow>distmax) distmax=distnow
						end do
						distavg=distavg/nrealatthas(irealattde)
						write(*,"(' Note: Attractor',i6,' has',i6,' members')") irealattde,nrealatthas(irealattde)
						write(*,"(' The minimum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmin,distmin*b2a
						write(*,"(' The maximum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmax,distmax*b2a
						write(*,"(' The average distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distavg,distavg*b2a
					end if
				end if

				if (imeasure==2.or.imeasure==3) then !Analyze one more term, then print angle
					if (ctmp3(1:1)=='a') then
						read(ctmp3(2:),*) iatm
						tmpx3=a(iatm)%x
						tmpy3=a(iatm)%y
						tmpz3=a(iatm)%z
					else if (ctmp3(1:1)=='c') then
						read(ctmp3(2:),*) irealatt
						tmpx3=realattxyz(1,irealatt)
						tmpy3=realattxyz(2,irealatt)
						tmpz3=realattxyz(3,irealatt)
					end if
				end if
				if (imeasure==2) then
					if (ctmp1(1:1)/='d') then
						write(*,"(' The angle is',f12.6,' degree')") xyz2angle(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
					else if (ctmp1(1:1)=='d') then
						angmin=180
						angmax=-1
						angavg=0
						distmin=999999
						distmax=-1
						distavg=0
						prjx=0
						prjy=0
						prjz=0
						do itmp=1,nrealatthas(irealattde)
							iatt=realatttable(irealattde,itmp)
							tmpx1=attxyz(1,iatt)
							tmpy1=attxyz(2,iatt)
							tmpz1=attxyz(3,iatt)
							angnow=xyz2angle(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
							angavg=angavg+angnow
							if (angnow<angmin) angmin=angnow
							if (angnow>angmax) angmax=angnow
							distnow=potlinedis(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
							distavg=distavg+distnow
							if (distnow<distmin) distmin=distnow
							if (distnow>distmax) distmax=distnow
							call pointprjline(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3,prjxtmp,prjytmp,prjztmp)
							prjx=prjx+prjxtmp
							prjy=prjy+prjytmp
							prjz=prjz+prjztmp
						end do
						angavg=angavg/nrealatthas(irealattde)
						distavg=distavg/nrealatthas(irealattde)
						prjx=prjx/nrealatthas(irealattde)
						prjy=prjy/nrealatthas(irealattde)
						prjz=prjz/nrealatthas(irealattde)
						write(*,"(' Note: Attractor',i6,' has',i6,' members')") irealattde,nrealatthas(irealattde)
						write(*,"(' The minimum angle is',f12.6,' degree')") angmin
						write(*,"(' The maximum angle is',f12.6,' degree')") angmax
						write(*,"(' The average angle is',f12.6,' degree')") angavg
						write(*,"(' The minimum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmin,distmin*b2a
						write(*,"(' The maximum distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distmax,distmax*b2a
						write(*,"(' The average distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") distavg,distavg*b2a
						write(*,"(' The average X,Y,Z coordinate by projecting the members of ',a,' to the line linking ',a,' and ',a)") trim(ctmp1),trim(ctmp2),trim(ctmp3)
						write(*,"(3f14.8,'   Angstrom')") prjx*b2a,prjy*b2a,prjz*b2a
					end if
				end if
				
				if (imeasure==3) then !Analyze one more term, then print dihedral angle
					if (ctmp4(1:1)=='a') then
						read(ctmp4(2:),*) iatm
						tmpx4=a(iatm)%x
						tmpy4=a(iatm)%y
						tmpz4=a(iatm)%z
					else if (ctmp4(1:1)=='c') then
						read(ctmp4(2:),*) irealatt
						tmpx4=realattxyz(1,irealatt)
						tmpy4=realattxyz(2,irealatt)
						tmpz4=realattxyz(3,irealatt)
					end if
					write(*,"(' The dihedral angle is',f12.6,' degree')") xyz2dih_sign(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3,tmpx4,tmpy4,tmpz4)
				end if
			else
				write(*,*) "Error: Invalid input!"
			end if
		end do
		
	else if (isel==-1) then
		do while(.true.)
			write(*,*)
			write(*,*) "0 Return"
			if (igridmethod==1) write(*,*) "1 Select algorithm for generating basins, current: On-grid"
			if (igridmethod==2) write(*,*) "1 Select algorithm for generating basins, current: Near-grid"
			if (igridmethod==3) write(*,*) "1 Select algorithm for generating basins, current: Near-grid with refinement"
			if (ibasinlocmin==0) write(*,*) "2 Switch the object to locate if all grid data are positive, current: Maxima"
			if (ibasinlocmin==1) write(*,*) "2 Switch the object to locate if all grid data are positive, current: Minima"
			read(*,*) isel2
			if (isel2==0) then
				exit
			else if (isel2==1) then
				write(*,*) "1: On-grid method, Comput .Mat. Sci., 36, 354 (2006)"
				write(*,*) "2: Near-grid method, J. Phys.: Condens. Matter, 21, 08420 (2009)"
				write(*,*) "3: Near-grid method with boundary refinement step"
				write(*,"(a)") " Note: Near-grid method (adapted by Tian Lu) is more accurate than On-grid method and thus &
				&s more recommended; with the boundary refinement step, the result will be better"
				read(*,*) igridmethod
			else if (isel2==2) then
				if (ibasinlocmin==0) then
					ibasinlocmin=1
				else
					ibasinlocmin=0
				end if
			end if
		end do
	else if (isel==0) then
		ioldtextheigh=textheigh
		textheigh=40 !Default textheigh is too small
		call drawbasinintgui
		textheigh=ioldtextheigh
		
	else if (isel==1) then !Regenerate basins and relocate attractors
		!When previously exported basinana.txt and basinana.cub are available, directly load rather than regenerate
		inquire(file="basinana.cub",exist=alive1)
		inquire(file="basinana.txt",exist=alive2)
		if ((alive1.eqv..true.).and.(alive2.eqv..true.)) then
			write(*,*) "basinana.txt and basinana.cub are found in current folder, &
            directly load attractors&basins information and grid data from them? (y/n)"
            read(*,*) selectyn
            if (selectyn=='y') then
				call load_basinana_info
				!Set range for looping over grids. For isolated case, grids at boundary should be ignored to avoid move outside
				if (ifPBC==0) then
					ixlow=2;iylow=2;izlow=2
					ixup=nx-1;iyup=ny-1;izup=nz-1
				else
					ixlow=1;iylow=1;izlow=1
					ixup=nx;iyup=ny;izup=nz
				end if
				cycle
            end if
        end if
		!Preparing grid data for partitioning basins
		isourcedata=1
		if (allocated(cubmat)) then
			write(*,"(a)") " Note: There has been a grid data in the memory, please select generating the basins by which manner"
			write(*,*) "0 Return"
			write(*,*) "1 Generate the basins by selecting a real space function"
			write(*,*) "2 Generate the basins by using the grid data stored in memory"
            if (ifiletype==7.or.ifiletype==8) write(*,*) "3: Same as 2, and perform analysis as three-dimension periodic system" !Can be used if grid data is loaded from input file
			read(*,*) isourcedata
            if (isourcedata==3) then
				isourcedata=2
				call grid2cellinfo !Add cell information from grid information
            end if
		end if
		if (isourcedata==0) then
			cycle
		else if (isourcedata==1) then
			write(*,*) "Select the real space function used for partitioning basins"
			call selfunc_interface(1,ifuncbasin)
			call setgridforbasin(ifuncbasin)
			if (allocated(cubmat)) deallocate(cubmat)
			allocate(cubmat(nx,ny,nz))
            !For basin analysis of IRI-pi, IRI_rhocut should be set to zero to avoid automatically setting it to 5 in low rho region, which will lead to huge number of artificial extrema
            ichange=0
            if (ifuncbasin==100.and.iuserfunc==99.and.IRI_rhocut/=0) then
				write(*,"(a)") " Note: IRI_rhocut parameter has been temporarily set to 0 during calculating grid data"
				tmpval=IRI_rhocut
				IRI_rhocut=0
                ichange=1
            end if
            !call gen_GTFuniq(0) !Generate unique GTFs, for faster evaluation in orbderv. This is automatically done in savecubmat
			call savecubmat(ifuncbasin,0,iorbsel)
            if (ichange==1) then
				write(*,*) "Note: Original IRI_rhocut parameter has been restored"
				IRI_rhocut=tmpval
            end if
		else if (isourcedata==2) then
			write(*,"(' Grid vector 1 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv1,dsqrt(sum(gridv1**2))
			write(*,"(' Grid vector 2 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv2,dsqrt(sum(gridv2**2))
			write(*,"(' Grid vector 3 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv3,dsqrt(sum(gridv3**2))
			write(*,"(' Number of points in three directions is',3i5,'  Total:',i12)") nx,ny,nz,nx*ny*nz
		end if
        
        !Set range for looping over grids. For isolated case, grids at boundary should be ignored to avoid move outside
        if (ifPBC==0) then
			ixlow=2;iylow=2;izlow=2
			ixup=nx-1;iyup=ny-1;izup=nz-1
        else
			ixlow=1;iylow=1;izlow=1
			ixup=nx;iyup=ny;izup=nz
        end if
		
        !For the regions with negative value, invert its sign to positive and locate attractor as usual, then minima of the negative part can be located
		allocate(grdposneg(nx,ny,nz)) !.true./.false means this grid has positive/negative value
		grdposneg=.true.
		do iz=1,nz
			do iy=1,ny
				do ix=1,nx
					if (cubmat(ix,iy,iz)<0D0) then
						grdposneg(ix,iy,iz)=.false.
						cubmat(ix,iy,iz)=-cubmat(ix,iy,iz)
					end if
				end do
			end do
		end do
! 		where (cubmat<0D0)  !DO NOT USE THIS, because I found that "where" will consuming vary large amount of memory!
! 			grdposneg=.false.
! 			cubmat=-cubmat !Invert negative values to positive, after basins are generated the values will be recovered
! 		end where
        if (ibasinlocmin==1) then
            if (all(cubmat>=0)) then
                write(*,"(/,a)") " Note: Since ""ibasinlocmin"" in settings.ini has been set to 1, and all grid data have non-negative value, &
                &therefore minima (repulsors) rather than maxima (attractors) will be located. The attractors reported subsequently in fact correspond to minima"
		        grdposneg(:,:,:)=.false.
		        cubmat(:,:,:)=-cubmat(:,:,:)
            end if
        end if
        
        !Initialize grid array of recording basin indices
		if (allocated(gridbas)) deallocate(gridbas)
		allocate(gridbas(nx,ny,nz))
		gridbas=0 !Unassigned state
        if (ifPBC==0) then
			gridbas(1,:,:)=-2 !Use index of -2 to identify box boundary grid
			gridbas(nx,:,:)=-2
			gridbas(:,1,:)=-2
			gridbas(:,ny,:)=-2
			gridbas(:,:,1)=-2
			gridbas(:,:,nz)=-2
        end if
        
		call setupmovevec !Generate movement vectors, will be used in generatebasin
        
        !If .cub file is used, ask if generate core density as corerhogrid(:,:,:), which can be added to current grid &
        !data of valence electron density, so that attractors occur at nuclear positions
        ifcorerho=0
        if (ifiletype==7.and.(index(filename,".cub")/=0.or.index(filename,".cube")/=0)) then
			if (any(a%index/=a%charge).and.all(a%charge/=0)) then
				write(*,*)
				write(*,*) "Consider core electron density during generating AIM basins? (y/n)"
                write(*,"(a)") " Note: To use this feature, effective charge information of atoms in your cube file must correspond to actual number of valence electrons"
                write(*,*) "If your grid data does not contain electron density, input ""n"""
                read(*,*) selectyn
                if (selectyn=='y') then
					call getcorerhogrid !Generate core density grid data, for AIM basin partition based on cub file recording valence electron density
					ifcorerho=1
                    ifuncbasin=1 !The function used for partitioning is electron density
                    write(*,*) "Summing up valence density and core density grid data..."
                    cubmat(:,:,:)=cubmat(:,:,:)+corerhogrid(:,:,:) !Construct all-electron density
                end if
            end if
        end if
        
        !Generate basins now!!!!!!
		write(*,*)
		write(*,*) "Generating basins, please wait..."
		call walltime(walltime1)
		call generatebasin(igridmethod)
		call walltime(walltime2)
		write(*,"(' Generating basins took up wall clock time',i10,' s')") walltime2-walltime1
        
        if (ifcorerho==1) then !Restore to valence electron density
            write(*,*) "Restoring valence density..."
			cubmat(:,:,:)=cubmat(:,:,:)-corerhogrid(:,:,:)
        end if
        
		numunassign=count(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==0)
		write(*,"(' The number of unassigned grids:',i12)") numunassign
		numgotobound=count(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==-1)
		write(*,"(' The number of grids travelled to box boundary:',i12)") numgotobound
		where (grdposneg.eqv..false.) cubmat=-cubmat !Recover original grid data
		deallocate(grdposneg)
		
        !Eliminate the attractors with very low value
		do iatt=1,numatt
			if ( abs(cubmat(attgrid(1,iatt),attgrid(2,iatt),attgrid(3,iatt)))<threslowvalatt ) then
				write(*,"(' Note: There are attractors having very low absolute value (<',1PE8.2,') and thus insignificant, how to deal with them?')") threslowvalatt
				write(*,*) "1 Do nothing"
				write(*,*) "2 Set corresponding grids as unassigned status"
				write(*,*) "3 Assign corresponding grids to the nearest significant attractors"
				write(*,*) "Hint: For most cases, option 3 is recommended"
				read(*,*) isel2
				call elimlowvalatt(threslowvalatt,isel2)
				exit
			end if
		end do
		
		!Generate actual coordinate of attractors, which will be used to plot in the local GUI, and in any other external GUIs
		if (allocated(attxyz)) deallocate(attxyz,attval)
		allocate(attxyz(3,numatt),attval(numatt))
		do iatt=1,numatt
			ix=attgrid(1,iatt)
			iy=attgrid(2,iatt)
			iz=attgrid(3,iatt)
            call getgridxyz(ix,iy,iz,attxyz(1,iatt),attxyz(2,iatt),attxyz(3,iatt))
			attval(iatt)=cubmat(ix,iy,iz)
		end do
		!Currently attractors have been finally determined, one should not perturb them anymore
		
        !Cluster degenerate attractors as "real attractors" and calculate average coordinate and value for the real attractors
		call clusdegenatt(0)

        !Detect interbasin grids
		call detectinterbasgrd(6)
		numinterbas=count(interbasgrid.eqv..true.)
		write(*,"(' The number of interbasin grids:',i12)") numinterbas
		
	else if (isel==2) then !Integrate real space function using uniform grid
		call integratebasin
        
    else if (isel==3) then !Output multipole moments
        call multipolebasin
		
	else if (isel==4.or.isel==5.or.isel==6.or.isel==7.or.isel==-7.or.isel==8) then
		if (.not.allocated(b)) then
			write(*,"(a)") " Note: No GTF (Gaussian type function) information is available in your input file, please input the file &
			&containing GTF information of your system, such as .wfn/.wfx/.mwfn/.fch/.molden file. e.g. C:\abc.wfn"
			read(*,"(a)") c200tmp
			call readinfile(c200tmp,1)
            if (ifPBC==0) then
				call gen_GTFuniq(1) !Generate unique GTFs, for faster evaluation in orbderv
            else
				call gen_neigh_GTF !Generate neighbouring GTFs list at reduced grids, for faster evaluation
            end if
		end if
		if (isel==4) then !Output LI and DI
			call LIDIbasin
		else if (isel==5) then !Output BOM
			call outBOMAOM(1)
		else if (isel==6) then !Output AOM
			call outBOMAOM(2)
		else if (isel==7) then !Integrate real space function using mixed grid
			write(*,*) "0 Return"
			write(*,*) "1 Integrate a specific function with atomic-center + uniform grids"
			write(*,*) "2 The same as 1, but with exact refinement of basin boundary"
			write(*,*) "3 The same as 2, but with approximate refinement of basin boundary"
			write(*,*) "Hint:"
			write(*,*) "Accuracy: 2>=3>>1     Time spent: 2>3>>1     Memory requirement: 3>2=1"
	! 		write(*,*) "Robust: 1=2>3"
			read(*,*) iseltmp
            if (iseltmp==0) cycle
			call integratebasinmix(iseltmp)
		else if (isel==-7) then
			call integratebasinmix_LSB
		else if (isel==8) then !Calculate multipole moment using mixed grid
			call integratebasinmix(10)
		end if
		
! 	else if (isel==100) then
! 		call integratebasinrefine
	else if (isel==9) then
		call atmpopinbasin
    else if (isel==10) then
        call HELP_HELV
    else if (isel==11) then
        call basinorbcomp
    else if (isel==12) then
		call assignELFbasinlab
	end if
    
    !Because delvirorb has been called before using some options, restore previous orbital information here
    if ((isel>=1.and.isel<=3).or.(isel>=7.and.isel<=10)) call delvirorb_back(0)
	!call del_GTFuniq !Destory unique GTF informtaion
end do

end subroutine



!!------- Generate basins from regular grid
! igridmethod=1: On grid, Comput.Mat.Sci.,36,354
! =2: Near-grid method (slower, but more accurate), see J.Phys.:Condens.Matter,21,084204
! =3: Near-grid with refinement
! The near-grid method was improved by Tian Lu, namely at later stage automatically switches to on-grid method to guarantee convergence
subroutine generatebasin(igridmethod)
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
integer,parameter :: nmaxtrjgrid=3000
integer igridmethod
integer ntrjgrid !Recording how many members are contained in trjgrid now
integer trjgrid(3,nmaxtrjgrid) !The trajectory contains which grids, record their indices sequentially, trjgrid(1/2/3,i) = ix/iy/iz of the ith grid
! real*8 trjval(nmaxtrjgrid),gradmaxval(nmaxtrjgrid) !******For debugging******
if (allocated(attgrid)) deallocate(attgrid)
allocate(attgrid(3,nint(nx*ny*nz/20D0))) !I think the number of attractors in general is impossible to exceeds nx*ny*nz/20
numatt=0
write(*,*) "  Attractor       X,Y,Z coordinate (Angstrom)                Value"
!Cycle all grids, but the box boundary ones are ignored (due to gridbas=-2), since can't evalute its gradients in all directions by finite difference

nsteplimit=min( nmaxtrjgrid,nint(dsqrt(dfloat(nx*nx+ny*ny+nz*nz))*2) )
nstepdiscorr=nint(nsteplimit/2D0)
ntime=1
if (igridmethod==3) ntime=2 !The second time performs refinement, emit steepest ascent trajectories from every boundary grid

do itime=1,ntime
!$OMP PARALLEL DO private(ix,iy,iz,corrx,corry,corrz,ntrjgrid,inowx,inowy,inowz,trjgrid,valnow,imove,gradtmp,igradmax,gradmax,iatt,tmpx,tmpy,tmpz,&
!$OMP itrjgrid,idtmp,icorrx,icorry,icorrz,gradx,grady,gradz,sclgrad,ineiidx,ixtmp,iytmp,iztmp,ixtmp2,iytmp2,iztmp2) shared(gridbas,numatt,attgrid) schedule(DYNAMIC) NUM_THREADS(nthreads)
do iz=izlow,izup
	do iy=iylow,iyup
		do ix=ixlow,ixup
			if (itime==1) then
				if (gridbas(ix,iy,iz)/=0) cycle
            else !Refinement
				if (.not.interbasgrid(ix,iy,iz)) cycle !Only refine boundary grids
				if (gridbas(ix,iy,iz)<=0) cycle !Ignored the ones unassigned or gone to box boundary
            end if
			ntrjgrid=0 !Number of grids passed by this searching trajectory
			inowx=ix !Grid index of the first point of trajectory
			inowy=iy
			inowz=iz
			corrx=0D0 !Correction vector
			corry=0D0
			corrz=0D0
			do while(.true.) !Emit steepest ascent trajectory from every grid
				ntrjgrid=ntrjgrid+1
				if (ntrjgrid>nsteplimit) exit !These gridbas for these unconverged grids will still be 0
				trjgrid(1,ntrjgrid)=inowx
				trjgrid(2,ntrjgrid)=inowy
				trjgrid(3,ntrjgrid)=inowz
				!Test all 26 directions to find out the maximal gradient direction
				valnow=cubmat(inowx,inowy,inowz)
				do imove=1,26
					ixtmp=inowx+vec26x(imove)
					iytmp=inowy+vec26y(imove)
					iztmp=inowz+vec26z(imove)
					if (ifPBC>0) call PBCgrididx(ixtmp,iytmp,iztmp) !If this grid is outside cell, wrap its index into the cell
					gradtmp=(cubmat(ixtmp,iytmp,iztmp)-valnow)/len26(imove)
					if (imove==1.or.gradtmp>gradmax) then
						igradmax=imove
						gradmax=gradtmp
					end if
				end do
				!Test if this is an attractor (maximum gradient is negative)
				if (gradmax<=0D0) then !Equal sign is important, because when system has symmetry, adjacent grid may be degenerate about mirrow plane, now the ascent should be terminated
					if (valnow==0D0) exit !The region is far beyond system, the value may be exactly zero due to cutoff of exponent, these grids should not be regarded as attractors
					!$OMP CRITICAL
cyciatt:            do iatt=1,numatt !Test if current grid corresponds to an existing attractor 
						if (inowx==attgrid(1,iatt).and.inowy==attgrid(2,iatt).and.inowz==attgrid(3,iatt)) then
							if (itime==1) then !Assign all grids in this trajectory
								do itrjgrid=1,ntrjgrid
									gridbas(trjgrid(1,itrjgrid),trjgrid(2,itrjgrid),trjgrid(3,itrjgrid))=iatt
								end do
                            else !During refinement, only update the boundary grid
								gridbas(ix,iy,iz)=iatt
                            end if
							exit cyciatt
						end if
						!Test if neighbour grid is attractor iatt, is yes, assign all grids in this trajectory
						!The reason I do this is because when the points are symmetric to Cartesian plane, many adjacent and degenerate attractors will occur,&
						!while this treatment can combine them as a single one
						do imove=1,26
							ixtmp=inowx+vec26x(imove)
							iytmp=inowy+vec26y(imove)
							iztmp=inowz+vec26z(imove)
							if (ifPBC>0) call PBCgrididx(ixtmp,iytmp,iztmp) !If this grid is outside cell, wrap its index into the cell
							if (ixtmp==attgrid(1,iatt).and.iytmp==attgrid(2,iatt).and.iztmp==attgrid(3,iatt)) then
								if (itime==1) then
									do itrjgrid=1,ntrjgrid
										gridbas(trjgrid(1,itrjgrid),trjgrid(2,itrjgrid),trjgrid(3,itrjgrid))=iatt
									end do
                                else
									gridbas(ix,iy,iz)=iatt
                                end if
								exit cyciatt
							end if
						end do
					end do cyciatt
					!Found a new attractor. Add to list, and assign the grids in the current trajectory
					if (iatt==numatt+1) then
						if (itime==1) then
							numatt=numatt+1
							do itrjgrid=1,ntrjgrid
								gridbas(trjgrid(1,itrjgrid),trjgrid(2,itrjgrid),trjgrid(3,itrjgrid))=numatt
							end do
							attgrid(1,numatt)=inowx
							attgrid(2,numatt)=inowy
							attgrid(3,numatt)=inowz
							call getgridxyz(inowx,inowy,inowz,tmpx,tmpy,tmpz)
							if (grdposneg(inowx,inowy,inowz)) then
								write(*,"(i8,3f14.8,f20.8)") numatt,tmpx*b2a,tmpy*b2a,tmpz*b2a,cubmat(inowx,inowy,inowz)
							else !This grid should have negative value
								write(*,"(i8,3f14.8,f20.8)") numatt,tmpx*b2a,tmpy*b2a,tmpz*b2a,-cubmat(inowx,inowy,inowz)
							end if
                        else if (itime==2) then
							write(*,*) "Warning: Found new attractor at refining process!"
                        end if
					end if
					!$OMP end CRITICAL
					exit
				end if

                !On-grid method
                if (igridmethod==1) then
					!Move to next grid
					inowx=inowx+vec26x(igradmax)
					inowy=inowy+vec26y(igradmax)
					inowz=inowz+vec26z(igradmax)
					if (ifPBC>0) call PBCgrididx(inowx,inowy,inowz) !If current grid is outside cell, wrap its index into the cell
					!Test if this grid has already been assigned, if yes, assign all grids in this trajectory with same basin index and terminate this trajectory
					idtmp=gridbas(inowx,inowy,inowz)
					if (idtmp>0) then
						do itrjgrid=1,ntrjgrid
							gridbas(trjgrid(1,itrjgrid),trjgrid(2,itrjgrid),trjgrid(3,itrjgrid))=idtmp
						end do
						exit !Terminate this trajectory
					end if
                    
                !Near-grid method without/with refinement
                else if (igridmethod==2.or.igridmethod==3) then
					!Note: Correction step may lead to oscillator when encountering circular ELF/LOL attractor, so if the current number of step is already large
					!(larger than half of upper limit of step number), then correction will be disabled (namely switch to on-grid method) to guarantee convergence
					!  Correction vector is already large, this time we do correction step
					if ( ntrjgrid<nstepdiscorr .and. (abs(corrx)>(dx/2D0).or.abs(corry)>(dy/2D0).or.abs(corrz)>(dz/2D0)) ) then
						if (abs(corrx)>(dx/2D0)) then
							icorrx=nint(corrx/abs(corrx)) !Get sign of corrx
							inowx=inowx+icorrx
                            call PBCgrididx1(inowx)
							corrx=corrx-icorrx*dx
						end if
						if (abs(corry)>(dy/2D0)) then
							icorry=nint(corry/abs(corry))
							inowy=inowy+icorry
                            call PBCgrididx2(inowy)
							corry=corry-icorry*dy
						end if
						if (abs(corrz)>(dz/2D0)) then
							icorrz=nint(corrz/abs(corrz))
							inowz=inowz+icorrz
                            call PBCgrididx3(inowz)
							corrz=corrz-icorrz*dz
						end if
					else !Move to next grid according to maximal gradient and then meantime update correction vector
						!Calculate true gradient
						ixtmp=inowx+1
						ixtmp2=inowx-1
						iytmp=inowy+1
						iytmp2=inowy-1
						iztmp=inowz+1
						iztmp2=inowz-1
                        if (ifPBC>0) then
							call PBCgrididx1(ixtmp)
							call PBCgrididx1(ixtmp2)
							call PBCgrididx2(iytmp)
							call PBCgrididx2(iytmp2)
							call PBCgrididx3(iztmp)
							call PBCgrididx3(iztmp2)
                        end if
						gradx=(cubmat(ixtmp,inowy,inowz)-cubmat(ixtmp2,inowy,inowz))/(2*dx)
						grady=(cubmat(inowx,iytmp,inowz)-cubmat(inowx,iytmp2,inowz))/(2*dy)
						gradz=(cubmat(inowx,inowy,iztmp)-cubmat(inowx,inowy,iztmp2))/(2*dz)
						sclgrad=min(dx/abs(gradx),dy/abs(grady),dz/abs(gradz))
						corrx=corrx+gradx*sclgrad-vec26x(igradmax)*dx
						corry=corry+grady*sclgrad-vec26y(igradmax)*dy
						corrz=corrz+gradz*sclgrad-vec26z(igradmax)*dz
						inowx=inowx+vec26x(igradmax)
						inowy=inowy+vec26y(igradmax)
						inowz=inowz+vec26z(igradmax)
						if (ifPBC>0) call PBCgrididx(inowx,inowy,inowz)
					end if
                    
					!Test if this grid has already been assigned and all of the neighbours were also assigned to the same attractor
					idtmp=gridbas(inowx,inowy,inowz)
					if (idtmp>0) then !This grid has already been assigned
						if (itime==1.or.(itime==2.and.ntrjgrid>60)) then
							do imove=1,26 !Check neighbours
								ixtmp=inowx+vec26x(imove)
								iytmp=inowy+vec26y(imove)
								iztmp=inowz+vec26z(imove)
								if (ifPBC>0) call PBCgrididx(ixtmp,iytmp,iztmp)
								ineiidx=gridbas(ixtmp,iytmp,iztmp)
								if (ineiidx/=idtmp) exit
							end do
							!All neighbours of current grid have same basin index as current grid, indicating current trajectory &
							!has fully enter a basin, so all grids in this trajectory can be assigned to that basin index, and the trjectory can be terminated
							if (imove==27) then
								if (itime==1) then
									do itrjgrid=1,ntrjgrid
										gridbas(trjgrid(1,itrjgrid),trjgrid(2,itrjgrid),trjgrid(3,itrjgrid))=idtmp
									end do
                                else !Refinement
									gridbas(ix,iy,iz)=idtmp
                                end if
								exit
							end if
                        end if
					end if
                end if
                    
				!Test if encountered box boundary for isolated system
				if (ifPBC==0.and.(inowx==1.or.inowx==nx.or.inowy==1.or.inowy==ny.or.inowz==1.or.inowz==nz)) then
					do itrjgrid=1,ntrjgrid
						gridbas(trjgrid(1,itrjgrid),trjgrid(2,itrjgrid),trjgrid(3,itrjgrid))=-1
					end do
					exit !Terminate trajectory
				end if
			end do !End steepest ascent trajectory
		end do !End cycle x
	end do !End cycle y
end do !End cycle z
!$OMP END PARALLEL DO
if (igridmethod==3.and.itime==1) then
	write(*,*) "Detecting boundary grids..."
	call detectinterbasgrd(6) !It seems that using 26 directions to determine boundary grids doesn't bring evident benefit
	write(*,"(' There are',i12,' grids at basin boundary')") count(interbasgrid.eqv..true.)
	write(*,*) "Refining basin boundary..."
end if
end do

! The following code for refining is not needed, since it is done via above code at itime=2
!!!! Refining basin boundaries. Emit steepest ascent trajectory from every boundary grid
!if (igridmethod==3) then
!! 	do itime=1,3 !Refine once in general is sufficient
!	write(*,*) "Detecting boundary grids..."
!	call detectinterbasgrd(6) !It seems that using 26 directions to determine boundary grids doesn't bring evident benefit
!	write(*,"(' There are',i12,' grids at basin boundary')") count(interbasgrid.eqv..true.)
!	write(*,*) "Refining basin boundary..."
!	!Below code is the adapted copy of above near-grid code
!	!$OMP PARALLEL DO private(ix,iy,iz,corrx,corry,corrz,ntrjgrid,inowx,inowy,inowz,valnow,imove,gradtmp,igradmax,gradmax,iatt,&
!	!$OMP ineiidx,idtmp,icorrx,icorry,icorrz,gradx,grady,gradz,sclgrad) shared(gridbas,attgrid) schedule(DYNAMIC) NUM_THREADS(nthreads)
!	do iz=izlow,izup
!		do iy=iylow,iyup
!			do ix=ixlow,ixup
!				if (.not.interbasgrid(ix,iy,iz)) cycle
!				if (gridbas(ix,iy,iz)<=0) cycle !Ignored the ones unassigned or gone to box boundary
!				ntrjgrid=0
!				inowx=ix
!				inowy=iy
!				inowz=iz
!				corrx=0D0 !Correction vector
!				corry=0D0
!				corrz=0D0
!				do while(.true.) !Steepest ascent
!					ntrjgrid=ntrjgrid+1
!					if (ntrjgrid>nsteplimit) exit
!					valnow=cubmat(inowx,inowy,inowz)
!					do imove=1,26
!						gradtmp=(cubmat(inowx+vec26x(imove),inowy+vec26y(imove),inowz+vec26z(imove))-valnow)/len26(imove)
!						if (imove==1.or.gradtmp>gradmax) then
!							igradmax=imove
!							gradmax=gradtmp
!						end if
!					end do
!					if (gradmax<=0) then !Equal sign is important, because when system has symmetry, adjacent grid may be degenerate about mirrow plane, now the ascent should be terminated
!cyciatt2:               do iatt=1,numatt
!							if (inowx==attgrid(1,iatt).and.inowy==attgrid(2,iatt).and.inowz==attgrid(3,iatt)) then
!								gridbas(ix,iy,iz)=iatt
!								exit cyciatt2
!							end if
!							do imove=1,26 !Test if neighbour grid (+/-x,+/-y,+/-z) is attractor iatt
!								if (inowx+vec26x(imove)==attgrid(1,iatt).and.inowy+vec26y(imove)==attgrid(2,iatt).and.inowz+vec26z(imove)==attgrid(3,iatt)) then
!									gridbas(ix,iy,iz)=iatt
!									exit cyciatt2
!								end if
!							end do
!						end do cyciatt2
!						if (iatt>numatt) then
!							write(*,*) "Warning: Found new attractor at refining process!"
!						end if
!						exit
!					end if
!					if ( ntrjgrid<nstepdiscorr .and. (abs(corrx)>(dx/2D0).or.abs(corry)>(dy/2D0).or.abs(corrz)>(dz/2D0)) ) then !This time we do correction step
!						if (abs(corrx)>(dx/2D0)) then
!							icorrx=nint(corrx/abs(corrx)) !Get sign of corrx
!							inowx=inowx+icorrx
!							corrx=corrx-icorrx*dx
!						end if
!						if (abs(corry)>(dy/2D0)) then
!							icorry=nint(corry/abs(corry))
!							inowy=inowy+icorry
!							corry=corry-icorry*dy
!						end if
!						if (abs(corrz)>(dz/2D0)) then
!							icorrz=nint(corrz/abs(corrz))
!							inowz=inowz+icorrz
!							corrz=corrz-icorrz*dz
!						end if
!					else !Move to next grid according to maximal gradient and then update correction vector
!						gradx=(cubmat(inowx+1,inowy,inowz)-cubmat(inowx-1,inowy,inowz))/(2*dx)
!						grady=(cubmat(inowx,inowy+1,inowz)-cubmat(inowx,inowy-1,inowz))/(2*dy)
!						gradz=(cubmat(inowx,inowy,inowz+1)-cubmat(inowx,inowy,inowz-1))/(2*dz)
!						sclgrad=min(dx/abs(gradx),dy/abs(grady),dz/abs(gradz))
!						inowx=inowx+vec26x(igradmax)
!						inowy=inowy+vec26y(igradmax)
!						inowz=inowz+vec26z(igradmax)
!						corrx=corrx+gradx*sclgrad-vec26x(igradmax)*dx
!						corry=corry+grady*sclgrad-vec26y(igradmax)*dy
!						corrz=corrz+gradz*sclgrad-vec26z(igradmax)*dz
!					end if
!					idtmp=gridbas(inowx,inowy,inowz)
!					if (ntrjgrid>60.and.idtmp>0) then !If enable this doesn't affect result detectably, but enabling it will evidently reduce computational cost
!						do imove=1,26
!							ineiidx=gridbas(inowx+vec26x(imove),inowy+vec26y(imove),inowz+vec26z(imove))
!							if (ineiidx/=idtmp) exit
!						end do
!						if (imove==27) then
!							gridbas(ix,iy,iz)=idtmp
!							exit
!						end if
!					end if
!					if (inowx==1.or.inowx==nx.or.inowy==1.or.inowy==ny.or.inowz==1.or.inowz==nz) exit
!				end do !End steepest ascent trajectory
!			end do !End cycle x
!		end do !End cycle y
!	end do !End cycle z
!	!$OMP END PARALLEL DO
!! 	end do
!end if

end subroutine



!!------- Eliminate low value attractors and the corresponding basins
!imode=1: Do nothing
!imode=2: Set corresponding grids as unassigned status
!imode=3: Assign corresponding grids to the nearest significant attractors
subroutine elimlowvalatt(threslowvalatt,imode)
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
integer imode
real*8 threslowvalatt
integer highvalatt(numatt)
integer att2att(-2:numatt) !Convert indices of old attractors to new ones
integer atttypelist(numatt)

if (imode==1) then
	return
else
	icounthigh=0
	atttypelist=0
	do iatt=1,numatt !Generate list
		if (abs(cubmat(attgrid(1,iatt),attgrid(2,iatt),attgrid(3,iatt)))>=threslowvalatt) then
			icounthigh=icounthigh+1
			highvalatt(icounthigh)=iatt
			atttypelist(iatt)=1 !high
		end if
	end do
	nlowvalatt=numatt-icounthigh
	if (imode==2) then
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					if (atttypelist(gridbas(ix,iy,iz))==0) gridbas(ix,iy,iz)=0
				end do
			end do
		end do
	else if (imode==3) then
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
					if (atttypelist(gridbas(ix,iy,iz))==0) then !Find out which grid belongs to insignificant attractors
						shortmaxsqr=1D20
						do iatt=1,numatt !Will be attribute to the nearest significant attractors
							if (atttypelist(iatt)==1) then
								call getgridxyz(attgrid(1,iatt),attgrid(2,iatt),attgrid(3,iatt),xatt,yatt,zatt)
								disttmpsqr=(xatt-rnowx)**2+(yatt-rnowy)**2+(zatt-rnowz)**2
								if (disttmpsqr<shortmaxsqr) then
									ishortmax=iatt
									shortmaxsqr=disttmpsqr
								end if
							end if
						end do
						gridbas(ix,iy,iz)=ishortmax
					end if
				end do
			end do
		end do
	end if
	!Build conversion table between old and new attractors, so that the indices of attractors could be contiguous
	att2att=-3
	att2att(-2)=-2
	att2att(-1)=-1
	att2att(0)=0
	icount=0
	do iatt=1,numatt
		if (atttypelist(iatt)==0) cycle
		icount=icount+1
		att2att(iatt)=icount
	end do
	numatt=numatt-nlowvalatt
	!Update attgrid
	do iatt=1,numatt
		attgrid(:,iatt)=attgrid(:,highvalatt(iatt))
	end do
	!Final update of grid attribution
	do iz=izlow,izup
		do iy=iylow,iyup
			do ix=ixlow,ixup
				gridbas(ix,iy,iz)=att2att(gridbas(ix,iy,iz))
			end do
		end do
	end do
	write(*,"(i6,' insignificant attractors have been eliminated')") nlowvalatt
end if
end subroutine



!!------- !Detect interbasin grids for real attractors
!ndir can be 6 and 26, meaning if only primary 6 directions or all 26 directions are used to determine boundary grids
subroutine detectinterbasgrd(ndir)
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
integer ndir
if (allocated(interbasgrid)) deallocate(interbasgrid)
allocate(interbasgrid(nx,ny,nz))
interbasgrid=.false.
do iz=izlow,izup
	do iy=iylow,iyup
		do ix=ixlow,ixup
			idxnow=gridbas(ix,iy,iz)
			do imove=1,ndir !To reduce the number of interbasin grids and thus speed up displaying, I don't use full 26 direction as criterion
				ixtmp=ix+vec26x(imove)
                iytmp=iy+vec26y(imove)
                iztmp=iz+vec26z(imove)
                call PBCgrididx(ixtmp,iytmp,iztmp)
				idxmove=gridbas(ixtmp,iytmp,iztmp)
				if ((idxmove/=idxnow).and.idxmove/=-2) then !This is a grid between two or more basins, and is not adjacent to box boundary
					interbasgrid(ix,iy,iz)=.true.
					exit
				end if
			end do
		end do
	end do
end do
end subroutine



!!------- Cluster degenerate and adjacent attractors together
!!Nearest neighbour method is used to cluster the attractors, see Leach book p493. Each cluster corresponds to a final attractor
!If imode=0, run the whole subroutine; if imode=1, then external attconv is provided, and only some code in this routine is used to make the index contiguous
subroutine clusdegenatt(imode)
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
integer imode
integer neleatt,eleatt(numatt) !eleatt(i)=j means the ith element in clustering stage corresponds to actual jth attractor. nclustlist is its total current elements
integer passlist(numatt) !if the ith element is 1, means i attractor has passed clustering stage
integer ncluster,ncluele(numatt),cluele(numatt,numatt) !cluele(i,4/5/9...) means cluster i has element 4,5,9..., ncluele(i) is current size of cluster i, ncluster is total number of cluster
real*8 attdismat(numatt,numatt),attvalmat(numatt,numatt) !Distance matrix and relative value difference matrix of the original attractors
real*8 posvec1(3),posvec2(3)
distcrit=mergeattdist*dsqrt(dx*dx+dy*dy+dz*dz) !Distance criterion
irefined=0 !If 1, that means this combining process takes effects
idebugclusdegenatt=0 !If output debug information

if (imode==1) goto 2

if (allocated(attconv)) deallocate(attconv)
allocate(attconv(-2:numatt))
do iatt=-2,numatt !First assume that each attractor is a real attractor
	attconv(iatt)=iatt
end do

!Generate difference matrix for distance and relative value 
attdismat=0D0
attvalmat=0D0
do iatt=1,numatt
	valiatt=attval(iatt)
    posvec1(:)=attxyz(:,iatt)
	do jatt=iatt+1,numatt
		valjatt=attval(jatt)
		posvec2(:)=attxyz(:,jatt)
        if (ifPBC==0) then
			distdiff=dsqrt( sum((posvec1(:)-posvec2(:))**2) )
        else
			call nearest_dist(posvec1,posvec2,distdiff)
        end if
		attdismat(iatt,jatt)=distdiff
		valdiff=abs((valiatt-valjatt)/valiatt)
		attvalmat(iatt,jatt)=valdiff
		if (valdiff<valcritclus.and.distdiff<distcrit) irefined=1
	end do
end do
attdismat=attdismat+transpose(attdismat)
attvalmat=attvalmat+transpose(attvalmat)
do iatt=1,numatt
	attdismat(iatt,iatt)=attdismat(iatt,iatt)/2D0
	attvalmat(iatt,iatt)=attvalmat(iatt,iatt)/2D0
end do
if (irefined==0) then !Nothing to do. Construct real attractor information and then directly return
	numrealatt=numatt
	if (allocated(nrealatthas)) deallocate(nrealatthas,realattval,realattxyz,realatttable)
	allocate(nrealatthas(numrealatt),realattval(numrealatt),realattxyz(3,numrealatt),realatttable(numrealatt,1))
	nrealatthas=1
	do iatt=1,numatt
		realattxyz(:,iatt)=attxyz(:,iatt)
		realattval(iatt)=attval(iatt)
		realatttable(iatt,1)=iatt
	end do
	return 
else
	write(*,*) "Degenerate attractors detected, clustering them..."
	write(*,"(' Criterion for clustering: Interval <',f8.5,' Bohr, value difference <',f8.5,'%')") distcrit,valcritclus*100
end if
! do i=1,numatt !Print distance matrix
! 	do j=1,numatt
! 		write(10,*) i,j,attdismat(i,j)
! 	end do
! end do

!Clustering via nearest neighbour method, and accordingly update conversion list
passlist=0
do iatt=1,numatt
	if (passlist(iatt)==1) cycle
	if (count(attdismat(iatt,:)<distcrit)==1) then !This attractor has no neighbour
		if (passlist(iatt)==1) cycle
		cycle
	end if
	!Construct mapping table of attractor, capacity is neleatt
	neleatt=0 !The number of degenerate attractors in this batch
	do jatt=iatt,numatt
		if (attvalmat(iatt,jatt)<valcritclus) then
			neleatt=neleatt+1
			eleatt(neleatt)=jatt
			passlist(jatt)=1
		end if
	end do
	if (neleatt==1) cycle !Although this attractor has one or more neighbours, but no other attractor is degenerate with itself, so pass
	if (idebugclusdegenatt==1) then
		write(*,"(a,i5)") " Attractors in this time of clustering triggered by attractor:",iatt
		write(*,"(' The number of initial clusters:',i5,', corresponding to these attractor:')") neleatt
		write(*,"(15i5)") eleatt(1:neleatt)
	end if

	!Start clustering. Now you should forget actual index of attractors. Only consider element 1,2,3...
	!Initialize, assume that each element in this batch of degenerate attractors is a cluster
	ncluster=neleatt
	ncluele(:)=1
	do iclu=1,ncluster
		cluele(iclu,1)=iclu
	end do
	ntimeclu=0
	do while(.true.) !Gradually clustering until closet distance between any cluster pair is larger than criteria
		ntimeclu=ntimeclu+1
		imergeclu=0
		do iclu=1,ncluster !Note that ncluster remain unchanged, eliminated cluster is simply emptyed but still reserve its index
			if (ncluele(iclu)==0) cycle !Has been incorporated to others
			shortmax=9999999D0
			ishortmax=0
			do jclu=iclu+1,ncluster
				if (ncluele(jclu)==0) cycle
				!Calculate shortest distance between cluster i and j
				shortmaxtmp=9999999D0
				do ieletmp=1,ncluele(iclu)
					do jeletmp=1,ncluele(jclu)
						iele=cluele(iclu,ieletmp)
						jele=cluele(jclu,jeletmp)
						distele=attdismat(eleatt(iele),eleatt(jele))
						if (distele<shortmaxtmp) shortmaxtmp=distele
					end do
				end do
				if (shortmaxtmp<shortmax) then !Find out the closest cluster to cluster iclu
					shortmax=shortmaxtmp
					ishortmax=jclu
				end if
			end do
			if (shortmax<distcrit) then !Combine cluster ishortmax to cluster iclu
				imergeclu=1
				newlen=ncluele(iclu)+ncluele(ishortmax)
				cluele(iclu,ncluele(iclu)+1:newlen)=cluele(ishortmax,1:ncluele(ishortmax))
				ncluele(iclu)=newlen
				ncluele(ishortmax)=0 !Empty
			end if
		end do
		if (imergeclu==0) exit !Between all cluster pair the shortmax is longer than distcrit, therefore exit
	end do
	do iclu=1,ncluster !Update attractor conversion list
		if (ncluele(iclu)==0) cycle
		if (idebugclusdegenatt==1) then
			write(*,"(' Clustering finished after',i5,' times of cycle')") ntimeclu
			write(*,"(' Cluster',i5,', deriving from  attractor',i5,', contains:')") iclu,eleatt(iclu)
	 		write(*,"(' Element:  ',15i5)") cluele(iclu,1:ncluele(iclu))
			write(*,"(' Attractor:',15i5)") eleatt(cluele(iclu,1:ncluele(iclu)))
		end if
		do itmp=1,ncluele(iclu)
			ielimele=cluele(iclu,itmp)
			ielimatt=eleatt(ielimele)
			itargetatt=eleatt(iclu)
			attconv(ielimatt)=itargetatt
		end do
	end do
end do

! do iatt=1,numatt !Print conversion list just after clustering
! 	write(*,*) iatt,attconv(iatt)
! end do
!Update attractor conversion list to make the attractor index contiguous, and hence get "real" attractor
!e.g. old conversion list, after slash is the new one
!   1           1/1
!   2           2/2
!   ...         2/2
!  13           2/2
!  14          14/3 <---After above clustering, 14,16,17 will convert to 14. The first new entry at right column always matched corresponding left column, this is guaranteed by above clustering code
!  15           2/2
!  16          14/3
!  17          14/3
!  18          18/4
2 numrealatt=0
passlist=0 !Record which attractor has already been converted
do iatt=1,numatt
	if (passlist(iatt)==1) cycle
	numrealatt=numrealatt+1 !Of course, numrealatt always <=iatt, so will not unexpectly overlap data during update conversion list
	idestmp=attconv(iatt)
	do jatt=iatt,numatt
		if (attconv(jatt)==idestmp) then
			passlist(jatt)=1
			attconv(jatt)=numrealatt
		end if
	end do
end do
! do iatt=1,numatt !Print conversion list after reordered
! 	write(*,*) iatt,attconv(iatt)
! end do

call updaterealattprop

if (imode==1) return
!Update basin attribution
do iz=izlow,izup
	do iy=iylow,iyup
		do ix=ixlow,ixup
			gridbas(ix,iy,iz)=attconv(gridbas(ix,iy,iz))
		end do
	end do
end do
!Output final attractors
write(*,*) "The attractors after clustering:"
write(*,*) "   Index      Average X,Y,Z coordinate (Angstrom)               Value"
do irealatt=1,numrealatt
	write(*,"(i8,3f15.8,f20.8)") irealatt,realattxyz(:,irealatt)*b2a,realattval(irealatt)
end do
end subroutine



!!------- Update real properties of attractors (average value, xyz, the number of members)
subroutine updaterealattprop
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
if (allocated(nrealatthas)) deallocate(nrealatthas,realattval,realattxyz,realatttable)
allocate(nrealatthas(numrealatt),realattval(numrealatt),realattxyz(3,numrealatt),realatttable(numrealatt,numatt))
nrealatthas=0
realattval=0D0
realattxyz=0D0
do iatt=1,numatt
	irealatt=attconv(iatt)
	nrealatthas(irealatt)=nrealatthas(irealatt)+1
	realatttable(irealatt,nrealatthas(irealatt))=iatt
	realattxyz(:,irealatt)=realattxyz(:,irealatt)+attxyz(:,iatt)
	realattval(irealatt)=realattval(irealatt)+attval(iatt)
end do
do irealatt=1,numrealatt
	realattval(irealatt)=realattval(irealatt)/nrealatthas(irealatt)
	realattxyz(:,irealatt)=realattxyz(:,irealatt)/nrealatthas(irealatt)
end do
end subroutine



!!------------------ Set grid for generating basin, adapted from the subroutine "setgrid"
subroutine setgridforbasin(ifuncsel)
use defvar
use GUI
implicit real*8 (a-h,o-z)
integer :: iselexttype=3,ifuncsel
real*8 :: molxlen,molylen,molzlen,tmpx,tmpy,tmpz,rhocrit=1D-6
real*8 :: gridextdist=5D0,enlarbox=2.1D0,spclowqual=0.2D0,spcmedqual=0.1D0,spchighqual=0.06D0,spclunaqual=0.04D0,tmparr6(6)
character c80tmp*80,cubefilename*200

if (.not.allocated(b)) iselexttype=2 !Unable to set box using rho
do while(.true.)
	if (iselexttype==1) then
		orgx=minval(a%x)-gridextdist
		orgy=minval(a%y)-gridextdist
		orgz=minval(a%z)-gridextdist
		endx=maxval(a%x)+gridextdist
		endy=maxval(a%y)+gridextdist
		endz=maxval(a%z)+gridextdist
	else if (iselexttype==2) then
		orgx=minval( a(:)%x-enlarbox*vdwr_tianlu(a(:)%index) )
		orgy=minval( a(:)%y-enlarbox*vdwr_tianlu(a(:)%index) )
		orgz=minval( a(:)%z-enlarbox*vdwr_tianlu(a(:)%index) )
		endx=maxval( a(:)%x+enlarbox*vdwr_tianlu(a(:)%index) )
		endy=maxval( a(:)%y+enlarbox*vdwr_tianlu(a(:)%index) )
		endz=maxval( a(:)%z+enlarbox*vdwr_tianlu(a(:)%index) )
	end if
	molxlen=endx-orgx
	molylen=endy-orgy
	molzlen=endz-orgz
	ntotlow=(nint(molxlen/spclowqual)+1)*(nint(molylen/spclowqual)+1)*(nint(molzlen/spclowqual)+1)
	ntotmed=(nint(molxlen/spcmedqual)+1)*(nint(molylen/spcmedqual)+1)*(nint(molzlen/spcmedqual)+1)
	ntothigh=(nint(molxlen/spchighqual)+1)*(nint(molylen/spchighqual)+1)*(nint(molzlen/spchighqual)+1)
	ntotluna=(nint(molxlen/spclunaqual)+1)*(nint(molylen/spclunaqual)+1)*(nint(molzlen/spclunaqual)+1)
	
	write(*,*) "Please select a method for setting up grid"
	if (iselexttype==1) write(*,"(a,f10.5,a)") " -10 Set grid extension distance for mode 1~6, current: Fixed,",gridextdist," Bohr"
	if (iselexttype==2) write(*,"(a)") " -10 Set grid extension distance for mode 1~6, current: Adaptive"
	if (iselexttype==3) write(*,"(a)") " -10 Set grid extension distance for mode 1~6, current: Detect rho isosurface"
	if (iselexttype==1.or.iselexttype==2) then
		write(*,"(a,f4.2,a,i14)") " 1 Low quality grid, spacing=",spclowqual," Bohr, number of grids:    ",ntotlow
		write(*,"(a,f4.2,a,i14)") " 2 Medium quality grid, spacing=",spcmedqual," Bohr, number of grids: ",ntotmed
		write(*,"(a,f4.2,a,i14)") " 3 High quality grid, spacing=",spchighqual," Bohr, number of grids:   ",ntothigh
		write(*,"(a,f4.2,a,i14)") " 4 Lunatic quality grid, spacing=",spclunaqual," Bohr, number of grids:",ntotluna
	else
		write(*,"(a,f4.2,a,i14)") " 1 Low quality grid, spacing=",spclowqual," Bohr, cost: 1x"
		write(*,"(a,f4.2,a,i14)") " 2 Medium quality grid, spacing=",spcmedqual," Bohr, cost: 8x"
		write(*,"(a,f4.2,a,i14)") " 3 High quality grid, spacing=",spchighqual," Bohr, cost: 36x"
		write(*,"(a,f4.2,a,i14)") " 4 Lunatic quality grid, spacing=",spclunaqual," Bohr, cost: 120x"
	end if
	write(*,*) "5 Only input grid spacing, automatically set other parameters"
	write(*,*) "6 Only input the number of points in X,Y,Z, automatically set other parameters"
	write(*,*) "7 Input original point, translation vector and the number of points"
	write(*,*) "8 Set center position, grid spacing and box length"
	write(*,*) "9 Use grid setting of another cube file"
	write(*,*) "10 Set box of grid data visually using a GUI window"
    if (ifPBC>0) write(*,"(a)") " 11 Use translation vectors of current cell, manually specify origin, box lengths and grid spacing"
	read(*,*) igridsel
    
    if (igridsel==-10) then
		write(*,*) "Please select the type of extension distance:"
		write(*,*) "1 Use fixed value in each direction"
		write(*,*) "2 Adaptively determined according to van der Waals radius"
		write(*,"(a,1PE7.1)") " 3 Detect rho to make the box just accommodate its isosurface, isoval: ",rhocrit
		write(*,*) "Hint: In general mode 3 is recommended, mode 1 should be avoided to be used"
		read(*,*) iselexttype
		if (iselexttype==1) then
			do while(.true.)
				write(*,*) "Input extension distance (Bohr) e.g. 6.5"
				read(*,*) gridextdist
				if (gridextdist>0D0) exit
				write(*,*) "Error: The value must be larger than 0!" 
			end do
		else if (iselexttype==2) then
			do while(.true.)
				write(*,*) "Input the ratio for scaling vdW radii, e.g. 1.7"
				write(*,"(' Current value is',f12.6)") enlarbox
				read(*,*) enlarbox
				if (enlarbox>0D0) exit
				write(*,*) "Error: The value must be larger than 0!" 
			end do
		else if (iselexttype==3) then
			write(*,*) "Select the isovalue of rho"
			write(*,*) "1: 0.0001       (10^-4)"
			write(*,*) "2: 0.00001      (10^-5)"
			write(*,*) "3: 0.000005   (5*10^-6)"
			write(*,*) "4: 0.000003   (3*10^-6)"
			write(*,*) "5: 0.000001     (10^-6) In general recommended"
			write(*,*) "6: 0.0000005  (5*10^-7)"
			write(*,*) "7: 0.0000001    (10^-7)"
			write(*,*) "8: 0.00000005 (5*10^-8)"
			write(*,*) "9: 0.00000001   (10^-8)"
			write(*,*) "10: Input by yourself"
			read(*,*) iselrho
			if (iselrho==1) rhocrit=1D-4
			if (iselrho==2) rhocrit=1D-5
			if (iselrho==3) rhocrit=5D-6
			if (iselrho==4) rhocrit=3D-6
			if (iselrho==5) rhocrit=1D-6
			if (iselrho==6) rhocrit=5D-7
			if (iselrho==7) rhocrit=1D-7
			if (iselrho==8) rhocrit=5D-8
			if (iselrho==9) rhocrit=1D-8
			if (iselrho==10) then
				write(*,*) "Input the value, e.g. 0.000005 or 5D-6"
				read(*,*) rhocrit
			end if
		end if
	else
		exit
	end if
end do

!Properly set aug3D for some modes
if (igridsel==1.or.igridsel==2.or.igridsel==3.or.igridsel==4.or.igridsel==5.or.igridsel==6) then
	if (iselexttype==1) then
		aug3D=gridextdist !Set aug3D, which only affects GUI display, so that the scope of the axis is wide enough to contain the whole grid region
	else if (iselexttype==2) then
		aug3D=maxval(enlarbox*vdwr_tianlu(a(:)%index))*1.35D0
	else if (iselexttype==3) then
		call detectboxrho(rhocrit) !Determine orgx/y/z and endx/y/z
		molxlen=endx-orgx
		molylen=endy-orgy
		molzlen=endz-orgz
		tmparr6(1)=abs(orgx-minval(a%x))
		tmparr6(2)=abs(orgy-minval(a%y))
		tmparr6(3)=abs(orgz-minval(a%z))
		tmparr6(4)=abs(endx-maxval(a%x))
		tmparr6(5)=abs(endy-maxval(a%y))
		tmparr6(6)=abs(endz-maxval(a%z))
		aug3D=maxval(tmparr6)+1
	end if
end if

!Input information for various grid setting modes
!Note: orgx,orgy,orgz,endx,endy,endz as well as molx/y/zlen for igridsel==1~6 have already been set above
if (igridsel==1.or.igridsel==2.or.igridsel==3.or.igridsel==4.or.igridsel==5) then
	if (igridsel==1) dx=spclowqual
	if (igridsel==2) dx=spcmedqual
	if (igridsel==3) dx=spchighqual
	if (igridsel==4) dx=spclunaqual
	if (igridsel==5) then
		write(*,*) "Input the grid spacing (Bohr), e.g. 0.08"
		read(*,*) dx
	end if
	dy=dx
	dz=dx
	nx=nint(molxlen/dx)+1
	ny=nint(molylen/dy)+1
	nz=nint(molzlen/dz)+1
else if (igridsel==6) then
	write(*,*) "Input the number of grid points in X,Y,Z direction, e.g. 139,59,80"
	read(*,*) nx,ny,nz
	dx=molxlen/(nx-1)
	dy=molylen/(ny-1)
	dz=molzlen/(nz-1)
else if (igridsel==7) then
	write(*,*) "Input X,Y,Z coordinate of original point (Bohr), e.g. 0.1,4,-1"
	read(*,*) orgx,orgy,orgz
	write(*,*) "Input X,Y,Z component of translation vector (Bohr), e.g. 0.1,0.1,0.15"
	read(*,*) dx,dy,dz
	write(*,*) "Input the number of points in X,Y,Z direction, e.g. 139,59,80"
	read(*,*) nx,ny,nz
else if (igridsel==8) then
	write(*,*) "Input X,Y,Z coordinate of box center (in Angstrom)"
	write(*,*) "or input such as a8 to take the coordinate of atom 8 as box center"
	write(*,*) "or input such as a3,a7 to take the midpoint of atom 3 and atom 7 as box center"
	read(*,"(a)") c80tmp
	if (c80tmp(1:1)=='a') then
		do ich=1,len_trim(c80tmp)
			if (c80tmp(ich:ich)==',') exit
		end do
		if (ich==len_trim(c80tmp)+1) then
			read(c80tmp(2:),*) itmp
			cenx=a(itmp)%x
			ceny=a(itmp)%y
			cenz=a(itmp)%z
		else
			read(c80tmp(2:ich-1),*) itmp
			read(c80tmp(ich+2:),*) jtmp			
			cenx=(a(itmp)%x+a(jtmp)%x)/2D0
			ceny=(a(itmp)%y+a(jtmp)%y)/2D0
			cenz=(a(itmp)%z+a(jtmp)%z)/2D0
		end if
	else
		read(c80tmp,*) cenx,ceny,cenz
		cenx=cenx/b2a
		ceny=ceny/b2a
		cenz=cenz/b2a
	end if
	write(*,*) "Input the grid spacing (Bohr), e.g. 0.08"
	read(*,*) dx
	dy=dx
	dz=dx
	write(*,*) "Input the box lengths in X,Y,Z direction (Bohr), e.g. 8.0,8.0,13.5"
	read(*,*) molxlen,molylen,molzlen
	orgx=cenx-molxlen/2D0
	orgy=ceny-molylen/2D0
	orgz=cenz-molzlen/2D0
	endx=orgx+molxlen
	endy=orgy+molylen
	endz=orgz+molzlen
	nx=nint(molxlen/dx)+1
	ny=nint(molylen/dy)+1
	nz=nint(molzlen/dz)+1
else if (igridsel==9) then
	write(*,*) "Input path of a cube file, e.g. C:\Fate\Saber.cub"
	do while(.true.)
		read(*,"(a)") cubefilename
		inquire(file=cubefilename,exist=alive)
		if (alive) then
			open(10,file=cubefilename,status="old")
			read(10,*)
			read(10,*)
			read(10,*) nouse,orgx,orgy,orgz
			read(10,*) nx,dx
			read(10,*) ny,rnouse,dy
			read(10,*) nz,rnouse,rnouse,dz
			close(10)
			exit
		else
			write(*,*) "Error: File cannot be found, input again"
		end if
	end do
else if (igridsel==10) then
	call setboxGUI
else if (igridsel==11) then
	call setgrid_for_PBC(0.1D0,0)
end if


!For some modes, if system is symmetric to Cartesian plane, then slightly adjust grid setting, so that the distribution of grids are symmetric to Cartesian plane
!This treatment will make the integrals have much better symmetry
if (igridsel>=1.and.igridsel<=6) then
	diffx=abs(orgx+endx)
	diffy=abs(orgy+endy)
	diffz=abs(orgz+endz)
	if (diffx<0.05D0) then !The system is symmetry to YZ plane
		distxmin=1D10
		do ix=1,nx
			rnowx=orgx+(ix-1)*dx
			if (rnowx>=0D0.and.rnowx<distxmin) distxmin=rnowx
		end do
		!1D-12 is a perturbation to avoid the grids are degenerate with respect to Cartesian plane, which results in cumbersome degenerate attractors
		!However if it is introduced, the integrals will not so close to symmetry
		orgx=orgx+(dx/2D0-distxmin) !+1D-15
	end if
	if (diffy<0.05D0) then !The system is symmetry to XZ plane
		distymin=1D10
		do iy=1,ny
			rnowy=orgy+(iy-1)*dy
			if (rnowy>=0D0.and.rnowy<distymin) distymin=rnowy
		end do
		orgy=orgy+(dy/2D0-distymin) !+1D-15
	end if
	if (diffz<0.05D0) then !The system is symmetry to XY plane
		distzmin=1D10
		do iz=1,nz
			rnowz=orgz+(iz-1)*dz
			if (rnowz>=0D0.and.rnowz<distzmin) distzmin=rnowz
		end do
		orgz=orgz+(dz/2D0-distzmin) !+1D-15
	end if
end if

!Make dx,dy,dz and gridv1/2/3 available for all cases
if (igridsel==9.or.igridsel==11) then
    dx=gridv1(1)
    dy=gridv2(2)
    dz=gridv3(3)
else
    gridv1=0;gridv1(1)=dx
    gridv2=0;gridv2(2)=dy
    gridv3=0;gridv3(3)=dz
end if

call getgridend !Generate endx,endy,endz

write(*,"(' Coordinate of origin in X,Y,Z is   ',3f12.6,' Bohr')") orgx,orgy,orgz
write(*,"(' Coordinate of end point in X,Y,Z is',3f12.6,' Bohr')") endx,endy,endz
if (igridsel==9.or.igridsel==11) then !Grid may be nonorthogonal
	write(*,"(' Grid vector 1 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv1,dsqrt(sum(gridv1**2))
	write(*,"(' Grid vector 2 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv2,dsqrt(sum(gridv2**2))
	write(*,"(' Grid vector 3 in X,Y,Z is',3f10.6,' Bohr, norm:',f10.6)") gridv3,dsqrt(sum(gridv3**2))
	write(*,"(' Number of points in three directions is',3i5,'  Total:',i12)") nx,ny,nz,nx*ny*nz
else
	write(*,"(' Grid spacing in X,Y,Z is',3f12.6,' Bohr')") dx,dy,dz
	write(*,"(' Number of points in X,Y,Z is',3i5,'   Total:',i12)") nx,ny,nz,nx*ny*nz
end if

! do i=1,nx
! 	write(c80tmp,"(D20.13)") orgx+(i-1)*dx
! 	read(c80tmp,*) tmpval
! 	write(*,*) i," x ",tmpval
! end do
! do i=1,ny
! 	write(c80tmp,"(D20.13)") orgy+(i-1)*dy
! 	read(c80tmp,*) tmpval
! 	write(*,*) i," y ",tmpval
! end do
! do i=1,nz
! 	write(c80tmp,"(D20.13)") orgz+dfloat(i-1)*dz
! 	read(c80tmp,*) tmpval
! 	write(*,*) i," z ",tmpval
! end do
! read(*,*)
end subroutine


!!-- Detect proper box (namely set orgx/y/z,endx/y/z) so that the box can just enclose the isosurface of rho=rhocrit
subroutine detectboxrho(rhocrit)
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 rhocrit
tmpfac=1D0
orgxtmp=minval( a(:)%x-tmpfac*vdwr_tianlu(a(:)%index) ) !Define a too small box, extend each side (to obtain orgx/y/z,endx/y/z) by detecting rho
orgytmp=minval( a(:)%y-tmpfac*vdwr_tianlu(a(:)%index) )
orgztmp=minval( a(:)%z-tmpfac*vdwr_tianlu(a(:)%index) )
endxtmp=maxval( a(:)%x+tmpfac*vdwr_tianlu(a(:)%index) )
endytmp=maxval( a(:)%y+tmpfac*vdwr_tianlu(a(:)%index) )
endztmp=maxval( a(:)%z+tmpfac*vdwr_tianlu(a(:)%index) )
scanstp=0.3D0 !Spacing of scan
nstpx=nint((endxtmp-orgxtmp)/scanstp)+1 !The number of detection grid in X direction
nstpy=nint((endytmp-orgytmp)/scanstp)+1 !The number of detection grid in Y direction
nstpz=nint((endztmp-orgztmp)/scanstp)+1 !The number of detection grid in Z direction
distmove=0.2D0
write(*,*) "Detecting proper box size..."
!Detect Z low. Scan XY plane, if rho at a point is larger than "rhocrit", lower down the plane by "distmove" and check again, until all points have rho<rhocrit
orgz=orgztmp
do while(.true.)
	iok=1
zl:do ix=1,nstpx
		rnowx=orgxtmp+(ix-1)*scanstp
		do iy=1,nstpy
			rnowy=orgytmp+(iy-1)*scanstp
			if (fdens(rnowx,rnowy,orgz)>rhocrit) then
				iok=0
				exit zl
			end if
		end do
	end do zl
	if (iok==1) exit
	orgz=orgz-0.2D0 !Lower down the layer
end do
!Detect Z high
endz=endztmp
do while(.true.)
	iok=1
zh:do ix=1,nstpx
		rnowx=orgxtmp+(ix-1)*scanstp
		do iy=1,nstpy
			rnowy=orgytmp+(iy-1)*scanstp
			if (fdens(rnowx,rnowy,endz)>rhocrit) then
				iok=0
				exit zh
			end if
		end do
	end do zh
	if (iok==1) exit
	endz=endz+0.2D0
end do
!Detect X low. Scan YZ plane
orgx=orgxtmp
do while(.true.)
	iok=1
xl:do iy=1,nstpy
		rnowy=orgytmp+(iy-1)*scanstp
		do iz=1,nstpz
			rnowz=orgztmp+(iz-1)*scanstp
			if (fdens(orgx,rnowy,rnowz)>rhocrit) then
				iok=0
				exit xl
			end if
		end do
	end do xl
	if (iok==1) exit
	orgx=orgx-0.2D0
end do
!Detect X high
endx=endxtmp
do while(.true.)
	iok=1
xh:do iy=1,nstpy
		rnowy=orgytmp+(iy-1)*scanstp
		do iz=1,nstpz
			rnowz=orgztmp+(iz-1)*scanstp
			if (fdens(endx,rnowy,rnowz)>rhocrit) then
				iok=0
				exit xh
			end if
		end do
	end do xh
	if (iok==1) exit
	endx=endx+0.2D0
end do
!Detect Y low. Scan XZ plane
orgy=orgytmp
do while(.true.)
	iok=1
yl:do ix=1,nstpx
		rnowx=orgxtmp+(ix-1)*scanstp
		do iz=1,nstpz
			rnowz=orgztmp+(iz-1)*scanstp
			if (fdens(rnowx,orgy,rnowz)>rhocrit) then
				iok=0
				exit yl
			end if
		end do
	end do yl
	if (iok==1) exit
	orgy=orgy-0.2D0
end do
!Detect Y high
endy=endytmp
do while(.true.)
	iok=1
yh:do ix=1,nstpx
		rnowx=orgxtmp+(ix-1)*scanstp
		do iz=1,nstpz
			rnowz=orgztmp+(iz-1)*scanstp
			if (fdens(rnowx,endy,rnowz)>rhocrit) then
				iok=0
				exit yh
			end if
		end do
	end do yh
	if (iok==1) exit
	endy=endy+0.2D0
end do
end subroutine



!!------- Integrate a real space function in the basins already partitioned
!It is not assumed that wavefunction information is available. This routine also works purely based on grid data
subroutine integratebasin
use defvar
use util
use functions
use basinintmod
implicit real*8 (a-h,o-z)
real*8 intval(-1:numatt),intvalpriv(-1:numatt),basinvol(-1:numatt),basinvolpriv(-1:numatt),basinvdwvol(-1:numatt),basinvdwvolpriv(-1:numatt),posvec(3)
integer att2atm(numrealatt) !The attractor corresponds to which atom. If =0, means this is a NNA
integer walltime1,walltime2
character grdfilename*200

if (ifuncbasin==1.and.ifPBC==0) then
	write(*,"(a)") " IMPORTANT NOTE: To integrate AIM basins, it is much better to use option 7 &
    &to integrate the basins based on uniform + atomic center grid, the accuracy is significantly better than using the present function!"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
end if

write(*,*)
write(*,*) "Please select the integrand:"
write(*,*) "-2 Return"
write(*,*) "-1 The values of the grid data stored in an external file (.cub/.grd)"
write(*,*) "0 The values of the grid data stored in memory"
call selfunc_interface(1,ifuncint)
if (ifuncint==-1) then !Load external grid data to cubmattmp
	do while(.true.)
		write(*,*) "Input another .cub or .grd file name, e.g. C:\CODE_GEASS\CC.cub"
		read(*,"(a)") grdfilename
		inquire(file=grdfilename,exist=alive)
		if (alive) exit
		write(*,*) "Error: File not found, input again"
		write(*,*)
	end do
	inamelen=len_trim(grdfilename)
	if (grdfilename(inamelen-2:inamelen)=="cub".or.grdfilename(inamelen-3:inamelen)=="cube") then
		call readcubetmp(grdfilename,1,inconsis)
	else if (grdfilename(inamelen-2:inamelen)=="grd") then
		call readgrdtmp(grdfilename,inconsis)
	end if
	if (inconsis==1) return
	write(*,*)
else if (ifuncint==-2) then
	return
end if

call walltime(walltime1)
write(*,*) "Integrating, please wait..."
intval=0D0
basinvol=0D0
basinvdwvol=0D0
ifinish=0
!$OMP PARALLEL private(ix,iy,iz,irealatt,rnowx,rnowy,rnowz,tmpval,intvalpriv,basinvolpriv,basinvdwvolpriv) shared(intval,basinvol,basinvdwvol,ifinish) NUM_THREADS(nthreads)
intvalpriv=0D0
basinvolpriv=0D0
basinvdwvolpriv=0D0
!$OMP do schedule(DYNAMIC)
do iz=izlow,izup
	do iy=iylow,iyup
		do ix=ixlow,ixup
			call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
			if (ifuncint==-1) then
				tmpval=cubmattmp(ix,iy,iz)
			else if (ifuncint==0) then
				tmpval=cubmat(ix,iy,iz)
			else
				tmpval=calcfuncall(ifuncint,rnowx,rnowy,rnowz)
			end if
			irealatt=gridbas(ix,iy,iz)
			intvalpriv(irealatt)=intvalpriv(irealatt)+tmpval
			basinvolpriv(irealatt)=basinvolpriv(irealatt)+1
			if (tmpval>=0.001D0) basinvdwvolpriv(irealatt)=basinvdwvolpriv(irealatt)+1
		end do
	end do
    !$OMP CRITICAL
    ifinish=ifinish+1
    if (ifuncint>=1) call showprog(ifinish,nz-2)
    !$OMP END CRITICAL
end do
!$OMP end do
!$OMP CRITICAL
    intval=intval+intvalpriv
    basinvol=basinvol+basinvolpriv
    basinvdwvol=basinvdwvol+basinvdwvolpriv
!$OMP end CRITICAL
!$OMP END PARALLEL
if (ifinish<nz-2) call showprog(nz-2,nz-2)
call calc_dvol(dvol)
intval=intval*dvol
basinvol=basinvol*dvol !Basin volume
basinvdwvol=basinvdwvol*dvol !Basin volume within isosurface of 0.001 a.u.
write(*,*) "  #Basin        Integral(a.u.)      Volume(a.u.^3)"
do irealatt=1,numrealatt
	write(*,"(i8,f22.10,f20.8)") irealatt,intval(irealatt),basinvol(irealatt)
end do
write(*,"(' Sum of above values:',f20.8)") sum(intval(1:numrealatt))
if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0)
if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1)

call walltime(walltime2)
write(*,"(' Integrating basins took up wall clock time',i10,' s')") walltime2-walltime1

if (ifuncbasin==1.and.(ifuncint==-1.or.ifuncint==0.or.ifuncint==1)) then
	!Wavefunction information may be not available, so the correspondence between attractors and atoms are only crudely determined
    att2atm=0
	do iatt=1,numrealatt !Cycle each attractors
		do iatm=1,ncenter !Try to find corresponding atom
			if (ifPBC==0) then
				disttest=dsqrt( (realattxyz(1,iatt)-a(iatm)%x)**2+(realattxyz(2,iatt)-a(iatm)%y)**2+(realattxyz(3,iatt)-a(iatm)%z)**2 )
            else
				posvec(1)=a(iatm)%x
				posvec(2)=a(iatm)%y
				posvec(3)=a(iatm)%z
				call nearest_dist(realattxyz(:,iatt),posvec,disttest)
            end if
			if (disttest<0.3D0) att2atm(iatt)=iatm !If distance between attractor and a nucleus is smaller than 0.3 Bohr, then the attractor will belong to the atom
        end do
    end do
    !Print atomic charges
	write(*,"(/,a)") " Because the basins are partitioned based on electron density, and the integrand is electron density, &
    &atomic charges together with basin volumes and the atomic vdW volumes (within rho=0.001 a.u. isosurface) are printed below"
    rnormfac=1 !It is not possible to determine normalization factor of number of electron, since if purely grid data is used, we do not know net charge at all
    write(*,*)
	write(*,*) "    Atom      Basin       Charge (e)     Volume (Bohr^3)   vdW Volume (Bohr^3)"
    do iatm=0,ncenter
		do iatt=1,numrealatt
			if (att2atm(iatt)==iatm.and.iatm==0) then
				write(*,"(7x,'(NNA)',i7,4x,f14.8,f16.3,3x,f16.3)") iatt,-intval(iatt)/rnormfac,basinvol(iatt),basinvdwvol(iatt)
			else if (att2atm(iatt)==iatm.and.iatm/=0) then
				if (nEDFelec==0) then !Normal case, all electron basis or using pseudopotential but not accompanied by EDF
					write(*,"(i7,' (',a,')',i7,4x,f14.8,f16.3,3x,f16.3)") iatm,a(iatm)%name,iatt,a(iatm)%charge-intval(iatt)/rnormfac,basinvol(iatt),basinvdwvol(iatt)
				else !EDF is used, so using a(iatm)%index instead of a(iatm)%charge
					write(*,"(i7,' (',a,')',i7,4x,f14.8,f16.3,3x,f16.3)") iatm,a(iatm)%name,iatt,a(iatm)%index-intval(iatt)/rnormfac,basinvol(iatt),basinvdwvol(iatt)
				end if
			end if
		end do
	end do
end if
end subroutine



!!------- Calculate orbital compositions contributed by various basins by uniform grid integration
subroutine basinorbcomp
use defvar
use util
use basinintmod
use functions
implicit real*8 (a-h,o-z)
real*8 intval(-1:numatt),intvalpriv(-1:numatt),orbcomp(ncenter)
real*8 comparr(nmo,numrealatt),comparrpriv(nmo,numrealatt),orbval(nmo)
integer att2atm(numrealatt),atm2att(ncenter)
character c80tmp*80,c2000tmp*2000

if (allocated(frag1)) deallocate(frag1)
if (ifuncbasin==1) then
    call calc_att2atm(att2atm)
    atm2att=0 !If an atom doesn't have corresponding attractor, the attractor index will be 0
    do irealatt=1,numrealatt
        atm2att(att2atm(irealatt))=irealatt
    end do
end if
call calc_dvol(dvol)

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
	!write(*,"(a)") " -2: Print atom contribution to a batch of orbitals"
	!write(*,"(a)") " -3: Print fragment contribution to a batch of orbitals"
	write(*,"(a)") " -4: Export composition of every basin in every orbital to orbcomp.txt"
	!write(*,"(a)") " -5: Print orbital delocalization index (ODI) for a batch of orbitals"
    read(*,"(a)") c80tmp
    
    if (c80tmp=="0") then
        if (allocated(frag1)) deallocate(frag1)
        return
    else if (c80tmp=="-1") then
        call showorbinfo(1,nmo)
        
    else if (c80tmp=="-4") then
        if (ifuncbasin/=1) then
            write(*,"(a)") " Error: This option is meaningful only when the function used to partition the basins is electron density!"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            cycle
        end if
        comparr=0D0
        ifinish=0
        !$OMP PARALLEL private(ix,iy,iz,irealatt,rnowx,rnowy,rnowz,comparrpriv,orbval) shared(ifinish,comparr) NUM_THREADS(nthreads)
        comparrpriv=0D0
        !$OMP do schedule(DYNAMIC)
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
				    call orbderv(1,1,nmo,rnowx,rnowy,rnowz,orbval(:))
			        irealatt=gridbas(ix,iy,iz)
			        comparrpriv(:,irealatt)=comparrpriv(:,irealatt)+orbval(:)**2
		        end do
	        end do
			!$OMP CRITICAL
            ifinish=ifinish+1
            call showprog(ifinish,nz-2)
			!$OMP END CRITICAL
        end do
        !$OMP end do
        !$OMP CRITICAL
        comparr=comparr+comparrpriv
        !$OMP end CRITICAL
        !$OMP END PARALLEL
        comparr=comparr*dvol
    
		open(10,file="orbcomp.txt",status="replace")
		do imo=1,nmo
			write(10,"(' Orbital',i6)") imo
            do iatm=1,ncenter
                orbcomp(iatm)=comparr(imo,atm2att(iatm))
            end do
            valnorm=sum(comparr(imo,:)) !Since there may be NNA, normalization factor should use sum of all basins
			do iatm=1,ncenter
				write(10,"(i6,f11.3,' %')") iatm,orbcomp(iatm)/valnorm*100
			end do
		end do
		close(10)
		write(*,*) "Done! orbcomp.txt has been exported to current folder"
        if (numrealatt/=ncenter) then
            write(*,"(a)") " Note: The number of attractors is unequal to number of atoms, therefore the sum of atom contributions is not equal to 100 %"
        end if
        
    else if (c80tmp=="-9") then
        if (ifuncbasin/=1) then
            write(*,"(a)") " Error: This option is meaningful only when the function used to partition the basins is electron density!"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            cycle
        end if
        call definefragment
        
    else
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
    
		write(*,"(' Orbital:',i5,'  Energy(a.u.):',f14.6,'  Occ:',f10.5,'  Type: ',a)") ishowmo,MOene(ishowmo),MOocc(ishowmo),orbtypename(MOtype(ishowmo))
        write(*,*) "Please wait..."
		intval=0D0
        ifinish=0
        !$OMP PARALLEL private(ix,iy,iz,irealatt,rnowx,rnowy,rnowz,tmpval,intvalpriv) shared(intval,ifinish) NUM_THREADS(nthreads)
        intvalpriv=0D0
        !$OMP do schedule(DYNAMIC)
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
                    call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
				    tmpval=fmo(rnowx,rnowy,rnowz,ishowmo)**2
			        irealatt=gridbas(ix,iy,iz)
			        intvalpriv(irealatt)=intvalpriv(irealatt)+tmpval
		        end do
	        end do
			!$OMP CRITICAL
            ifinish=ifinish+1
            call showprog(ifinish,nz-2)
			!$OMP END CRITICAL
        end do
        !$OMP end do
        !$OMP CRITICAL
        intval=intval+intvalpriv
        !$OMP end CRITICAL
        !$OMP END PARALLEL
        intval=intval*dvol
        sumcontri=sum(intval(1:numrealatt))
        write(*,"(' Sum of raw contributions:',f18.3,' %')") sumcontri*100
        intval=intval/sumcontri
        write(*,*)
        write(*,*) "Final data after normalization:"
        do irealatt=1,numrealatt
	        write(*,"(' Basin:',i8,'   Contribution:',f10.3,' %')") irealatt,intval(irealatt)*100
        end do
        !For AIM basins, also show atomic values
        if (ifuncbasin==1) then
            write(*,*)
            write(*,*) "Contributions from atoms:"
            do iatm=0,ncenter
			    do iatt=1,numrealatt
				    if (att2atm(iatt)==iatm) then
                        if (iatm==0) then
					        write(*,"(i7,' (NNA)   Contribution:',f8.3,' %')") iatt,intval(iatt)*100
				        else
				            write(*,"(i7,' (',a,')    Contribution:',f8.3,' %')") iatm,a(iatm)%name,intval(iatt)*100
				        end if
                    end if
			    end do
		    end do
            write(*,*)
            write(*,"(' Orbital delocalization index:',f8.2)") sum((intval(1:numrealatt)*100)**2)/100
		    if (allocated(frag1)) then
			    write(*,*)
                compfrag=sum(intval(atm2att(frag1)))
			    write(*,"(' Composition of the fragment:',f8.3,' %')") compfrag*100
                orbdeloc=sum((intval(atm2att(frag1))*100/compfrag)**2)
                write(*,"(' Orbital delocalization index of the fragment:',f8.2)") orbdeloc/100
		    end if
        end if
    end if
end do
end subroutine




!!------- Obtain atomic population in a basin region (e.g. ELF bond basin) via AIM partition
subroutine atmpopinbasin
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
inquire(file="basin.cub",exist=alive)
if (.not.alive) then
	write(*,*) "Error: basin.cub is not existed in current folder!"
	return
else
	if (allocated(cubmattmp)) deallocate(cubmattmp)
    write(*,*) "Loading basin.cub in current folder..."
	call readcubetmp("basin.cub",2,inconsis)
	if (inconsis==1) then
		write(*,*) "Error: The grid setting of basin.cub is inconsistent with present grid data"
		return
	end if
end if
do while(.true.)
	write(*,*)
	write(*,*) "Study population of which atom? Input the index of corresponding attractor"
	write(*,*) "For example, 3"
	write(*,*) "Input 0 can exit"
	read(*,*) iatt
	if (iatt==0) then
		exit
	else if (iatt<0.or.iatm>numrealatt) then
		write(*,*) "Error: Attractor index exceeded valid range!"
		return
	end if
	write(*,*) "Study its population in which basin of the basin.cub file? e.g. 6"
	read(*,*) ibasin
	atmpop=0
	do iz=izlow,izup
		do iy=iylow,iyup
			do ix=ixlow,ixup
				if (nint(cubmattmp(ix,iy,iz))==ibasin.and.gridbas(ix,iy,iz)==iatt) atmpop=atmpop+cubmat(ix,iy,iz)
			end do
		end do
	end do
    call calc_dvol(dvol)
	atmpop=atmpop*dvol
	write(*,"(' Population of attractor',i4,' in external basin',i4,' is',f12.5)") iatt,ibasin,atmpop
end do
deallocate(cubmattmp)
end subroutine



!!----------------- Calculate multipole moment in the basins via uniform grid
subroutine multipolebasin
use defvar
use basinintmod
use functions
implicit real*8 (a-h,o-z)

if (ifuncbasin==1.and.ifPBC==0.and.allocated(b)) then
	write(*,"(a)") " IMPORTANT NOTE: To calculate multipole moments for AIM basins, it is much better to use option 8 &
    &to employ uniform + atomic center grid, the accuracy is significantly better than the present function!"
    write(*,*) "Press ENTER button to continue"
    read(*,*)
end if

do while(.true.)
	write(*,*)
	write(*,*) "Input the index of the basin, e.g. 5"
	write(*,"(a)") " Note: -1 means printing all basin results on screen, -2 means printing to multipole.txt in current folder. Input 0 can return"
	read(*,*) itmp
    ioutid=6
	if (itmp==-1.or.itmp==-2) then
		ibegin=1
		iend=numrealatt
		if (itmp==-2) then
			ioutid=10
			open(10,file="multipole.txt",status="replace")
		end if
	else if (itmp==0) then
		return
	else
		ibegin=itmp
		iend=itmp
	end if
	write(*,*) "Calculating, please wait..."
	write(ioutid,*) "Note: All data shown below are in a.u.!"
	write(ioutid,*)
	
	dipelextot=0D0
	dipeleytot=0D0
	dipeleztot=0D0
	eleinttot=0D0
	do ibas=ibegin,iend
	! 	xcen=0 !Use centroid of basin as center, but this is a bad idea
	! 	ycen=0
	! 	zcen=0
	! 	nbasgrid=0
	!	do iz=izlow,izup
	!		do iy=iylow,iyup
	!			do ix=ixlow,ixup
	! 				if (gridbas(ix,iy,iz)==ibas) then
	! 					nbasgrid=nbasgrid+1
	! 					call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
	! 					tmpdens=ELF_LOL(rnowx,rnowy,rnowz,"ELF")
	! 					xcen=xcen+rnowx*tmpdens
	! 					ycen=ycen+rnowy*tmpdens
	! 					zcen=zcen+rnowz*tmpdens
	! 				end if
	! 			end do
	! 		end do
	! 	end do
	! 	xcen=xcen/nbasgrid
	! 	ycen=ycen/nbasgrid
	! 	zcen=zcen/nbasgrid
		xcen=realattxyz(1,ibas)
		ycen=realattxyz(2,ibas)
		zcen=realattxyz(3,ibas)
		call calc_dvol(dvol)
	! 	write(*,"(' The X,Y,Z of the center of the basin:')")
	! 	write(*,"(3f14.8,' Bohr',/)") xcen,ycen,zcen
	
		eleint=0D0
		xint=0D0
		yint=0D0
		zint=0D0
		xxint=0D0
		yyint=0D0
		zzint=0D0
		xyint=0D0
		yzint=0D0
		xzint=0D0
		!$OMP PARALLEL private(ix,iy,iz,rnowx,rnowy,rnowz,tmpmul,rx,ry,rz,eleintp,xintp,yintp,zintp,xxintp,yyintp,zzintp,xyintp,yzintp,xzintp) NUM_THREADS(nthreads)
		eleintp=0D0
		xintp=0D0
		yintp=0D0
		zintp=0D0
		xxintp=0D0
		yyintp=0D0
		zzintp=0D0
		xyintp=0D0
		yzintp=0D0
		xzintp=0D0
		!$OMP do schedule(DYNAMIC)
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					if (gridbas(ix,iy,iz)==ibas) then
                        call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
						if (ifuncbasin==1) then
							tmpmul=dvol*cubmat(ix,iy,iz) !The cubmat currently is just electron density
						else
							tmpmul=dvol*fdens(rnowx,rnowy,rnowz)
						end if
						rx=rnowx-xcen
						ry=rnowy-ycen
						rz=rnowz-zcen
						eleintp=eleintp+tmpmul !monopole
						xintp=xintp+rx*tmpmul
						yintp=yintp+ry*tmpmul
						zintp=zintp+rz*tmpmul
						xxintp=xxintp+rx*rx*tmpmul
						yyintp=yyintp+ry*ry*tmpmul
						zzintp=zzintp+rz*rz*tmpmul
						xyintp=xyintp+rx*ry*tmpmul
						yzintp=yzintp+ry*rz*tmpmul
						xzintp=xzintp+rx*rz*tmpmul
					end if
				end do
			end do
		end do
		!$OMP end do
		!$OMP CRITICAL
		eleint=eleint+eleintp
		xint=xint+xintp
		yint=yint+yintp
		zint=zint+zintp
		xxint=xxint+xxintp
		yyint=yyint+yyintp
		zzint=zzint+zzintp
		xyint=xyint+xyintp
		yzint=yzint+yzintp
		xzint=xzint+xzintP
		!$OMP end CRITICAL
		!$OMP END PARALLEL
		if (itmp==-1.or.itmp==-2) write(ioutid,"(/,' ***** Basin',i8)") ibas
		if (itmp==-2) write(*,"(' Outputting basin',i8)") ibas
		write(ioutid,"(' Basin monopole moment:',f12.6)") -eleint
		write(ioutid,"(' Basin dipole moment:')") 
		write(ioutid,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") -xint,-yint,-zint,sqrt(xint**2+yint**2+zint**2)
		eleinttot=eleinttot+eleint
		dipelex=-eleint*xcen+(-xint) !Contribution to molecular total dipole moment
		dipeley=-eleint*ycen+(-yint)
		dipelez=-eleint*zcen+(-zint)
		dipelextot=dipelextot+dipelex
		dipeleytot=dipeleytot+dipeley
		dipeleztot=dipeleztot+dipelez
		write(ioutid,"(' Basin electron contribution to molecular dipole moment:')")
		write(ioutid,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") dipelex,dipeley,dipelez,sqrt(dipelex**2+dipeley**2+dipelez**2)
		write(ioutid,"(' Basin quadrupole moment (Traceless Cartesian form):')")
		rrint=xxint+yyint+zzint
		QXX=-(3*xxint-rrint)/2
		QYY=-(3*yyint-rrint)/2
		QZZ=-(3*zzint-rrint)/2
		write(ioutid,"(' XX=',f12.6,'  XY=',f12.6,'  XZ=',f12.6)") QXX,-(3*xyint)/2,-(3*xzint)/2
		write(ioutid,"(' YX=',f12.6,'  YY=',f12.6,'  YZ=',f12.6)") -(3*xyint)/2,QYY,-(3*yzint)/2
		write(ioutid,"(' ZX=',f12.6,'  ZY=',f12.6,'  ZZ=',f12.6)") -(3*xzint)/2,-(3*yzint)/2,QZZ
		write(ioutid,"(' Magnitude of the traceless quadrupole moment tensor:',f12.6)") dsqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
		R20=-(3*zzint-rrint)/2D0 !Notice that the negative sign, because electrons carry negative charge
		R2n1=-dsqrt(3D0)*yzint
		R2p1=-dsqrt(3D0)*xzint
		R2n2=-dsqrt(3D0)*xyint
		R2p2=-dsqrt(3D0)/2D0*(xxint-yyint)
		write(ioutid,"(' Basin quadrupole moments (Spherical harmonic form):')")
		write(ioutid,"(' Q_2,0 =',f11.6,'   Q_2,-1=',f11.6,'   Q_2,1=',f11.6)") R20,R2n1,R2p1
		write(ioutid,"(' Q_2,-2=',f11.6,'   Q_2,2 =',f11.6)") R2n2,R2p2
		write(ioutid,"(' Magnitude: |Q_2|=',f12.6)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)
        ESEx=xxint;ESEy=yyint;ESEz=zzint
        ESE=ESEx+ESEy+ESEz
        	write(ioutid,"(a,f16.6)") " Basin electronic spatial extent <r^2>:",ESE
		write(ioutid,"(' Components of <r^2>:  X=',f15.6,'  Y=',f15.6,'  Z=',f15.6)") ESEx,ESEy,ESEz
	end do
	if (itmp==-1.or.itmp==-2) then !Output overall properties, most users are not interested in them
		dipnucx=sum(a(:)%x*a(:)%charge)
		dipnucy=sum(a(:)%y*a(:)%charge)
		dipnucz=sum(a(:)%z*a(:)%charge)
		write(ioutid,*)
		write(ioutid,"( ' Molecular net charge:',f12.6)") sum(a%charge)-eleinttot
		write(ioutid,"( ' Nuclear contribution to molecular dipole moment:')") 
		write(ioutid,"(' X=',f14.6,'  Y=',f14.6,'  Z=',f14.6,'  Norm=',f14.6)") dipnucx,dipnucy,dipnucz,sqrt(dipnucx**2+dipnucy**2+dipnucz**2)
		write(ioutid,"( ' Electron contribution to molecular dipole moment:')") 
		write(ioutid,"(' X=',f14.6,'  Y=',f14.6,'  Z=',f14.6,'  Norm=',f14.6)") dipelextot,dipeleytot,dipeleztot,sqrt(dipelextot**2+dipeleytot**2+dipeleztot**2)
		dipmolx=dipnucx+dipelextot
		dipmoly=dipnucy+dipeleytot
		dipmolz=dipnucz+dipeleztot
		write(ioutid,"( ' Molecular dipole moment:')")
		write(ioutid,"(' X=',f14.6,'  Y=',f14.6,'  Z=',f14.6,'  Norm=',f14.6)") dipmolx,dipmoly,dipmolz,sqrt(dipmolx**2+dipmoly**2+dipmolz**2)
	end if
	if (itmp==-2) then
		close(10)
		write(*,*) "Outputting finished!"
	end if
end do
end subroutine


!!--------------- Calculate localized index within and delocalization index between basins
subroutine LIDIbasin
use defvar
use basinintmod
use util
implicit real*8 (a-h,o-z)
real*8 LI(numrealatt),LIa(numrealatt),LIb(numrealatt) !Localization index array
real*8 DI(numrealatt,numrealatt),DIa(numrealatt,numrealatt),DIb(numrealatt,numrealatt) !Delocalization index matrix
real*8 atmDI(ncenter,ncenter),atmDIa(ncenter,ncenter),atmDIb(ncenter,ncenter)
integer maplist(ncenter)
character selectyn,label*20

write(*,*) "Generating basin overlap matrix (BOM)..."
call genBOM !Generate BOM
nmatsize=size(BOM,1)
nmatsizeb=size(BOMb,1)

if (any(MOocc<0)) then
	where(MOocc<0) MOocc=0
	write(*,"(a)") " Note: Some orbital occupation numbers are negative. In order to make the calculation feasible, they have been set to zero"
	write(*,*) "Press ENTER button to continue"
	read(*,*)
end if

write(*,*) "Calculating LI and DI..."
!RHF,R-post-HF, DI_A,B=2∑[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B     where i and j are spatial orbitals
if (wfntype==0.or.wfntype==3) then
	DI=0D0
	do ibas=1,numrealatt
		do jbas=ibas,numrealatt
			do iorb=1,nmatsize
				do jorb=1,nmatsize
					DI(ibas,jbas)=DI(ibas,jbas)+dsqrt(MOocc(iorb)*MOocc(jorb))*BOM(iorb,jorb,ibas)*BOM(iorb,jorb,jbas)
				end do
			end do
		end do
		LI(ibas)=DI(ibas,ibas)
	end do
	DI=2*(DI+transpose(DI))
	do ibas=1,numrealatt !Diagonal terms are the sum of corresponding row or column
		DI(ibas,ibas)=sum(DI(ibas,:))-DI(ibas,ibas)
	end do
!ROHF
else if (wfntype==2) then
	DIa=0D0
	DIb=0D0
	do nmoclose=nmatsize,1,-1
		if (MOtype(nmoclose)==0) exit
	end do
	do ibas=1,numrealatt
		do jbas=ibas,numrealatt
			!Alpha
			do iorb=1,nmatsize !The number of close or alpha orbitals needed to be concerned
				occi=MOocc(iorb)
				if (MOtype(iorb)==0) occi=occi/2D0
				do jorb=1,nmatsize
					occj=MOocc(jorb)
					if (MOtype(jorb)==0) occj=occj/2D0
					DIa(ibas,jbas)=DIa(ibas,jbas)+dsqrt(occi*occj)*BOM(iorb,jorb,ibas)*BOM(iorb,jorb,jbas)
				end do
			end do
			!Beta
			do iorb=1,nmoclose !The number of close orbitals needed to be concerned
				do jorb=1,nmoclose
					DIb(ibas,jbas)=DIb(ibas,jbas)+dsqrt(MOocc(iorb)/2D0*MOocc(jorb)/2D0)*BOM(iorb,jorb,ibas)*BOM(iorb,jorb,jbas)
				end do
			end do
		end do
		LIa(ibas)=DIa(ibas,ibas)
		LIb(ibas)=DIb(ibas,ibas)
	end do
	DIa=2*(DIa+transpose(DIa))
	DIb=2*(DIb+transpose(DIb))
	do ibas=1,numrealatt !Diagonal terms are the sum of corresponding row or column
		DIa(ibas,ibas)=sum(DIa(ibas,:))-DIa(ibas,ibas)
		DIb(ibas,ibas)=sum(DIb(ibas,:))-DIb(ibas,ibas)
	end do
	!Combine alpha and Beta to total
	DI=DIa+DIb
	LI=LIa+LIb
!UHF,U-post-HF   DI(A,B)=2∑[i,j]dsqrt(n_i*n_j)*S_i,j_A * S_i,j_B   where i and j are spin orbitals
else if (wfntype==1.or.wfntype==4) then
	!Alpha
	DIa=0D0
	do ibas=1,numrealatt
		do jbas=ibas,numrealatt
			do iorb=1,nmatsize
				do jorb=1,nmatsize
					DIa(ibas,jbas)=DIa(ibas,jbas)+dsqrt(MOocc(iorb)*MOocc(jorb))*BOM(iorb,jorb,ibas)*BOM(iorb,jorb,jbas)
				end do
			end do
		end do
		LIa(ibas)=DIa(ibas,ibas)
	end do
	DIa=2*(DIa+transpose(DIa))
	!Beta
	if (nmatsizeb>0) then
		DIb=0D0
		do iendalpha=nmo,1,-1
			if (MOtype(iendalpha)==1) exit
		end do
		MOinit=iendalpha+1 !Index range of beta orbitals
		MOend=iendalpha+nmatsizeb
		do ibas=1,numrealatt
			do jbas=ibas,numrealatt
				do iorb=MOinit,MOend
					iorbtmp=iorb-iendalpha
					do jorb=MOinit,MOend
						jorbtmp=jorb-iendalpha
						DIb(ibas,jbas)=DIb(ibas,jbas)+dsqrt(MOocc(iorb)*MOocc(jorb))*BOMb(iorbtmp,jorbtmp,ibas)*BOMb(iorbtmp,jorbtmp,jbas)
					end do
				end do
			end do
			LIb(ibas)=DIb(ibas,ibas)
		end do
		DIb=2*(DIb+transpose(DIb))
	end if
	do ibas=1,numrealatt !Diagonal terms are the sum of corresponding row or column
		DIa(ibas,ibas)=sum(DIa(ibas,:))-DIa(ibas,ibas)
		DIb(ibas,ibas)=sum(DIb(ibas,:))-DIb(ibas,ibas)
	end do
	!Combine alpha and beta to total
	DI=DIa+DIb
	LI=LIa+LIb
end if

!Output LI and DI in basin indices
label=""
if (ifuncbasin==1) label=" (basin index) "
ioutid=6
write(*,*)
100 write(ioutid,"(a)") " Note: Diagonal terms of the following matrices are the sum of corresponding row or column elements"
if (wfntype==1.or.wfntype==2.or.wfntype==4) then !UHF,ROHF,U-post-HF, output each spin component first
	!Alpha
    write(ioutid,*)
	call showmatgau(DIa,"Delocalization index matrix"//trim(label)//" for alpha spin",0,"f14.8",ioutid,formindex="5x,i5,4x")
	write(ioutid,*)
	write(ioutid,*) "Localization index"//trim(label)//" for alpha spin:"
	do ibas=1,numrealatt
		write(ioutid,"(i5,':',f9.5)",advance='no') ibas,LIa(ibas)
		if (mod(ibas,5)==0) write(ioutid,*)
	end do
    if (mod(numrealatt,5)/=0) write(ioutid,*)
	!Beta
	write(ioutid,*)
	call showmatgau(DIb,"Delocalization index matrix"//trim(label)//" for beta spin",0,"f14.8",ioutid,formindex="5x,i5,4x")
	write(ioutid,*)
	write(ioutid,*) "Localization index"//trim(label)//" for beta spin:"
	do ibas=1,numrealatt
		write(ioutid,"(i5,':',f9.5)",advance='no') ibas,LIb(ibas)
		if (mod(ibas,5)==0) write(ioutid,*)
	end do
	if (mod(numrealatt,5)/=0) write(ioutid,*)
end if
!Alpha+Beta
write(ioutid,*)
call showmatgau(DI,"Total delocalization index matrix"//trim(label),0,"f14.8",ioutid,formindex="5x,i5,4x")
write(ioutid,*)
write(ioutid,*) "Total localization index"//trim(label)//":"
do ibas=1,numrealatt
	write(ioutid,"(i5,':',f9.5)",advance='no') ibas,LI(ibas)
	if (mod(ibas,5)==0) write(ioutid,*)
end do
if (mod(numrealatt,5)/=0) write(ioutid,*)

!Output LI and DI in atom indices
if (ifuncbasin==1) then
    write(ioutid,*)
    call atmidx2attidx(maplist,1,ioutid)
    !UHF,ROHF,U-post-HF, output each spin component first
    if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	    !Alpha
        atmDIa=0
        atmDIb=0
        do iatm=1,ncenter
            do jatm=1,ncenter
                if (maplist(iatm)==0.or.maplist(jatm)==0) cycle !Corresponding basin index was not identified
                atmDIa(iatm,jatm)=DIa(maplist(iatm),maplist(jatm))
                atmDIb(iatm,jatm)=DIb(maplist(iatm),maplist(jatm))
            end do
        end do
		write(ioutid,*)
	    call showmatgau(atmDIa,"Delocalization index matrix (atom index) for alpha spin",0,"f14.8",ioutid,formindex="5x,i5,4x")
	    write(ioutid,*)
	    write(ioutid,*) "Localization index (atom index) for alpha spin:"
	    do iatm=1,ncenter
            if (maplist(iatm)==0) cycle
		    write(ioutid,"(i5,':',f9.5)",advance='no') iatm,LIa(maplist(iatm))
		    if (mod(iatm,5)==0) write(ioutid,*)
	    end do
		if (mod(ncenter,5)/=0) write(ioutid,*)
		write(ioutid,*)
		call showmatgau(atmDIb,"Delocalization index matrix (atom index) for beta spin",0,"f14.8",ioutid,formindex="5x,i5,4x")
		write(ioutid,*)
		write(ioutid,*) "Localization index (atom index) for beta spin:"
		do iatm=1,ncenter
			if (maplist(iatm)==0) cycle
			write(ioutid,"(i5,':',f9.5)",advance='no') iatm,LIb(maplist(iatm))
			if (mod(iatm,5)==0) write(ioutid,*)
		end do
		if (mod(ncenter,5)/=0) write(ioutid,*)
    end if
    
    !Alpha+Beta
    atmDI=0
    do iatm=1,ncenter
        do jatm=1,ncenter
            if (maplist(iatm)==0.or.maplist(jatm)==0) cycle !Corresponding basin index was not identified
            atmDI(iatm,jatm)=DI(maplist(iatm),maplist(jatm))
        end do
    end do
    write(ioutid,*)
    call showmatgau(atmDI,"Total delocalization index matrix (atom index)",0,"f14.8",ioutid,formindex="5x,i5,4x")
    write(ioutid,*)
    write(ioutid,*) "Total localization index (atom index):"
    do iatm=1,ncenter
        if (maplist(iatm)==0) cycle
	    write(ioutid,"(i5,':',f9.5)",advance='no') iatm,LI(maplist(iatm))
	    if (mod(iatm,5)==0) write(ioutid,*)
    end do
    if (mod(ncenter,5)/=0) write(ioutid,*)
end if

if (ioutid==10) then
	write(*,*) "Done!"
	close(10)
	return
end if

write(*,*)
write(*,*) "If also outputting above information to LIDI.txt in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	open(10,file="LIDI.txt",status="replace")
	ioutid=10
	goto 100
end if
end subroutine



!!------- Find basin indices corresponding to various atom indices in the case of AIM analysis
!maplist(i) is the attractor index that corresponds to atom i. If no corresponding can be found, the value is 0
!info: =1 show correspondence, =0 do not show
!ioutid: Destination of output information
!Distance <0.3 Bohr is the criterion used to identify correspondence
subroutine atmidx2attidx(maplist,info,ioutid)
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
integer maplist(ncenter),ioutid
distcrit=0.3D0
maplist=0
if (info==1) write(ioutid,*) "Detecting correspondence between basin and atom indices (criterion: <0.3 Bohr)"
do iatm=1,ncenter
	do ibas=1,numrealatt
        disttest=dsqrt( (realattxyz(1,ibas)-a(iatm)%x)**2+(realattxyz(2,ibas)-a(iatm)%y)**2+(realattxyz(3,ibas)-a(iatm)%z)**2 )
        if (disttest<distcrit) then
            if (info==1) write(ioutid,"(' Basin',i6,' corresponds to atom',i6,' (',a,')')") ibas,iatm,a(iatm)%name 
            maplist(iatm)=ibas
            exit
        end if
    end do
    if (ibas==numrealatt+1.and.info==1) write(ioutid,"(a,i6)") " Warning: Unable to find basin index that corresponds to atom",iatm
end do
end subroutine



!!------- Integrate a real space function in the basins with multi-level refinement, indenspensible for Laplacian
!DEPRECATED, since mixed type of grids perform better and faster
!Also, this routine is not compatible with nonorthogonal grid
!subroutine integratebasinrefine
!use defvar
!use util
!use function
!use basinintmod
!implicit real*8 (a-h,o-z)
!character c200tmp*200
!real*8 intval(-1:numatt),basinvol(-1:numatt),intvalpriv(-1:numatt),basinvolpriv(-1:numatt)
!integer walltime1,walltime2
!real*8 :: critlevel1=0.1D0,critlevel2=0.5D0,critlevel3=1D0
!integer :: nrefine1=1,nrefine2=2,nrefine3=5,nrefine4=7
!if (ifuncbasin/=1) then
!	write(*,*) "Error: This function is only applicable to AIM basins!"
!	return
!end if
!
!do while(.true.)
!	write(*,*) "-2 Return"
!	write(*,*) "-1 Print and set parameters for multi-level refinement"
!	call selfunc_interface(1,ifuncint)
!	if (ifuncint==-2) then
!		return
!	else if (ifuncint==-1) then
!		write(*,"(a)") " Note: A number n means each grid in corresponding value range will be transformed &
!		as n^3 grids around it during the integration to gain a higher integration accuracy. n=1 means the grids will remain unchanged. &
!		The ""value range"" referred here is the value range of the function used to generate basins"
!		write(*,*) "The value range and the times of refinement:"
!		write(*,"(' Smaller than ',f10.5,' :',i4)") critlevel1,nrefine1
!		write(*,"(' Between',f10.5,' and ',f10.5,' :',i4)") critlevel1,critlevel2,nrefine2
!		write(*,"(' Between',f10.5,' and ',f10.5,' :',i4)") critlevel2,critlevel3,nrefine3
!		write(*,"(' Larger than  ',f10.5,' :',i4)") critlevel3,nrefine4
!		write(*,*)
!		write(*,*) "Please input three thresholds to define the four ranges, e.g. 0.1,0.5,1"
!		write(*,*) "Note: Press ENTER button can retain current values unchanged"
!		read(*,"(a)") c200tmp
!		if (c200tmp/=' ') read(c200tmp,*) critlevel1,critlevel2,critlevel3
!		write(*,*) "Input the times of refinement for the grids in the four ranges, e.g. 1,2,5,7"
!		write(*,*) "Note: Press ENTER button can retain current values unchanged"
!		read(*,"(a)") c200tmp
!		if (c200tmp/=' ') read(c200tmp,*) nrefine1,nrefine2,nrefine3,nrefine4
!		write(*,*) "Done!"
!		write(*,*)
!	end if
!end do
!
!call walltime(walltime1)
!write(*,*) "Integrating, please wait..."
!intval=0D0
!basinvol=0D0
!ifinish=0
!!$OMP PARALLEL private(ix,iy,iz,ixref,iyref,izref,ndiv,irealatt,rnowx,rnowy,rnowz,rnowxtmp,rnowytmp,rnowztmp,orgxref,orgyref,orgzref,dxref,dyref,dzref,&
!!$OMP tmpval,tmpvalrefine,intvalpriv,basinvolpriv,nrefine) shared(intval,basinvol,ifinish) NUM_THREADS(nthreads)
!intvalpriv=0D0
!basinvolpriv=0D0
!!$OMP do schedule(DYNAMIC)
!do iz=izlow,izup
!	do iy=iylow,iyup
!		do ix=ixlow,ixup
!			if (cubmat(ix,iy,iz)<critlevel1) then
!				nrefine=nrefine1 !The number of point to represent each edge
!			else if (cubmat(ix,iy,iz)<critlevel2) then
!				nrefine=nrefine2
!			else if (cubmat(ix,iy,iz)<critlevel3) then
!				nrefine=nrefine3
!			else
!				nrefine=nrefine4
!			end if
!			ndiv=nrefine**3
!            call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
!			orgxref=rnowx-dx/2 !Take corner position as original point of microcycle
!			orgyref=rnowy-dy/2
!			orgzref=rnowz-dz/2
!			dxref=dx/nrefine
!			dyref=dy/nrefine
!			dzref=dz/nrefine
!			tmpval=0D0
!			do ixref=1,nrefine
!				do iyref=1,nrefine
!					do izref=1,nrefine
!						rnowxtmp=orgxref+(ixref-0.5D0)*dxref
!						rnowytmp=orgyref+(iyref-0.5D0)*dyref
!						rnowztmp=orgzref+(izref-0.5D0)*dzref
!						if (ifuncint==-1) then
!							tmpvalrefine=cubmattmp(ix,iy,iz)
!						else if (ifuncint==0) then
!							tmpvalrefine=cubmat(ix,iy,iz)
!						else
!							tmpvalrefine=calcfuncall(ifuncint,rnowxtmp,rnowytmp,rnowztmp)
!						end if
!						tmpval=tmpval+tmpvalrefine/ndiv
!					end do
!				end do
!			end do
!			irealatt=gridbas(ix,iy,iz)
!			intvalpriv(irealatt)=intvalpriv(irealatt)+tmpval
!			basinvolpriv(irealatt)=basinvolpriv(irealatt)+1
!		end do
!	end do
!    ifinish=ifinish+1
!    call showprog(ifinish,nz-2)
!end do
!!$OMP end do
!!$OMP CRITICAL
!    intval=intval+intvalpriv
!    basinvol=basinvol+basinvolpriv
!!$OMP end CRITICAL
!!$OMP END PARALLEL
!call calc_dvol(dvol)
!intval=intval*dvol
!basinvol=basinvol*dvol !Basin volume
!write(*,*) "  #Basin          Integral        Volume(a.u.^3)"
!do irealatt=1,numrealatt
!	write(*,"(i8,f22.10,f20.8)") irealatt,intval(irealatt),basinvol(irealatt)
!end do
!write(*,"(' Sum of above values:',f20.8)") sum(intval(1:numrealatt))
!if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0)
!if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1)
!
!call walltime(walltime2)
!write(*,"(' Integrating basins took up wall clock time',i10,' s')") walltime2-walltime1
!end subroutine



!!------- Integrate AIM basins using mixed atomic-center and uniform grids
! itype=1: Integrate specific real space function
! itype=2: Integrate specific real space function with exact refinement of basin boundary
! itype=3: Integrate specific real space function with approximate refinement of basin boundary by exact+linear interpolation
! itype=10: Produce electric multipole moments
!NNA and ECP are supported
!Notice that the grid generated must cover the whole space!
!It is assumed that wavefunction information is available
subroutine integratebasinmix(itype)
use defvar
use util
use functions
use basinintmod
use topo
implicit real*8 (a-h,o-z)
!p suffix means private variable for parallel mode
real*8 intval(-1:numrealatt,20),intvalp(-1:numrealatt,20) !Up to 20 functions can be evaluated and stored simultaneously
real*8 basinvol(-1:numrealatt),basinvolp(-1:numrealatt),basinvdwvol(-1:numrealatt),basinvdwvolp(-1:numrealatt) !vdW is used to obtain the basin volume enclosed by 0.001 isosurface of rho
real*8 trustrad(numrealatt),intbasinthread(numrealatt),intbasin(numrealatt)
real*8 dens,grad(3),hess(3,3),k1(3),k2(3),k3(3),k4(3),xarr(nx),yarr(ny),zarr(nz)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:)
real*8,allocatable :: rhogrid(:,:,:),rhogradgrid(:,:,:,:) !Used in shubin's 2nd project
real*8,allocatable :: prorhogrid(:,:,:) !Used in integrating deformation density
type(content),allocatable :: gridatt(:) !Record correspondence between attractor and grid
integer att2atm(numrealatt) !The attractor corresponds to which atom. If =0, means this is a NNA
real*8 eleint(-1:numrealatt),xint(-1:numrealatt),yint(-1:numrealatt),zint(-1:numrealatt),&
xxint(-1:numrealatt),yyint(-1:numrealatt),zzint(-1:numrealatt),xyint(-1:numrealatt),yzint(-1:numrealatt),xzint(-1:numrealatt)
real*8 eleintp(-1:numrealatt),xintp(-1:numrealatt),yintp(-1:numrealatt),zintp(-1:numrealatt),&
xxintp(-1:numrealatt),yyintp(-1:numrealatt),zzintp(-1:numrealatt),xyintp(-1:numrealatt),yzintp(-1:numrealatt),xzintp(-1:numrealatt)
integer radpotAIM,sphpotAIM
real*8 prodensgrad(0:4),gridval(100000,6) !Used for storing various information during integrating inside sphere
character c80tmp*80,selectyn
real*8 quadmom(3,3),tmpmat(3,3),tmpvec(3)

nbeckeiter=8
if (ifuncbasin/=1) then
	write(*,"(a)") " Error: This function is only applicable to AIM basins! That means in option 1, you should select electron density to construct the basins."
	return
end if

if (ispecial==0) then
	if (itype==1.or.itype==2.or.itype==3) then
		write(*,*) "Please select integrand:"
		write(*,*) "-2 Return"
		write(*,*) "-1 Deformation density"
		call selfunc_interface(1,ifuncint)
		if (ifuncint==-2) then
			return
		else if (ifuncint==-1) then
			call setpromol
		end if
    else if (itype==10) then
		write(*,"(a)") " Note: Atomic dipole/multipole moments will be calculated with respect to attractor of &
        &electron density rather than nuclear position, the discrepancy is usually very small"
	end if
else if (ispecial==1) then
	continue !Don't let user to select integrand
else if (ispecial==2) then
	call setpromol !Don't let user to select integrand
	expcutoff=1 !Use full accuracy for shubin's 2nd project
end if

numcp=0
att2atm=0
intval=0D0
basinvol=0D0
basinvdwvol=0D0
eleint=0D0
xint=0D0
yint=0D0
zint=0D0
xxint=0D0
yyint=0D0
zzint=0D0
xyint=0D0
yzint=0D0
xzint=0D0

call walltime(iwalltime1)

!Determine trust radius and then integrate in the trust sphere
!For each electron attractor determined by basin analysis, we identify its corresponding atom, then use Newton method to obtain exact NCP position
!The grid employed in this stage is centered from NCP
write(*,*)
write(*,*) "Integrating in trust sphere..."
do iatt=1,numrealatt !Cycle each attractors
	do iatm=1,ncenter !Try to find corresponding atom
		disttest=dsqrt( (realattxyz(1,iatt)-a(iatm)%x)**2+(realattxyz(2,iatt)-a(iatm)%y)**2+(realattxyz(3,iatt)-a(iatm)%z)**2 )
		if (disttest<0.3D0) then !If distance between attractor and a nucleus is smaller than 0.3 Bohr, then the attractor will belong to the atom
			att2atm(iatt)=iatm
			write(*,"(/,' Attractor',i6,' corresponds to atom',i6,' (',a,')')") iatt,iatm,a(iatm)%name
            !Find exact NCP position from the nucleus and record to CPpos
			numcpold=numcp
			call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,1) !If successfully converge to a NCP, the position is write to CPpos(:,numcp)
			if (numcp==numcpold) then !Failed to converge to a NCP
				write(*,*) "Note: Unable to locate exact CP position! Use nuclear position"
				numcp=numcp+1
				CPpos(1,numcp)=a(iatm)%x
				CPpos(2,numcp)=a(iatm)%y
				CPpos(3,numcp)=a(iatm)%z
			end if
! 			write(*,"(' Coordinate after refinement:',3f16.8)") CPpos(:,numcp)
			exit
		end if
	end do
	if (att2atm(iatt)==0) then !No real atom corresponds to this attractor
		write(*,"(/,a,i6,a)") " Warning: Unable to determine the attractor",iatt," belongs to which atom!"
		write(*,"(a)") " If this is a non-nuclear attractor, simply press ENTER button to continue. If you used pseudopotential &
		&and this attractor corresponds to the cluster of all maxima of its valence electron, then input the index of this atom (e.g. 9). &
		&Else you should input q to return and regenerate basins with smaller grid spacing"
		read(*,"(a)") c80tmp
		if (c80tmp=='q') then
			return
		else if (c80tmp==" ") then
			numcpold=numcp
			call findcp(realattxyz(1,iatt),realattxyz(2,iatt),realattxyz(3,iatt),1)
			if (numcp==numcpold) then
				write(*,*) "Unable to locate exact CP position! Use attractor position instead"
				numcp=numcp+1
				CPpos(:,numcp)=realattxyz(:,iatt)
			end if
		else !ECP, input corresponding atom by user and directly use its nuclear position as grid center
			read(c80tmp,*) iatmtmp
			att2atm(iatt)=iatmtmp
			numcp=numcp+1
			CPpos(1,numcp)=a(iatmtmp)%x
			CPpos(2,numcp)=a(iatmtmp)%y
			CPpos(3,numcp)=a(iatmtmp)%z
		end if
	end if
	
	!Set integration points and weights, and meantime determine trust radius 
	radpotAIM=200 !My test shows that increasing it to 400 may improve total integral by N*0.0001. This is unimportant, so only use 200, it is already good enough
	parm=1
	isettrustrad=0
	nintgrid=0 !Then number of integration grids within trust radius
	if (allocated(gridatt)) deallocate(gridatt) !Used to record grids in trust sphere of this attractor
	allocate(gridatt(radpotAIM*500))
	do ish=1,radpotAIM !Cycle each radial shell. Radius distance is from near to far
		!Becke, namely the second-kind Gauss-Chebyshev
		itmp=radpotAIM+1-ish !Invert ish to make radr from near to far
		radx=cos(itmp*pi/(radpotAIM+1D0))
		radr=(1+radx)/(1-radx)*parm
		radw=2*pi/(radpotAIM+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
! 		!Handy, also known as Euler-Maclaurin. See: Murray, C. W.; Handy, N. C.; Laming, G. J. Mol Phys 1993, 78, 997
! 		radx=dfloat(ish)/(radpotAIM+1D0)
! 		radr=radx**2/(1-radx)**2*parm
! 		radw=2*radx**5/dfloat(radpotAIM+1)/(1-radx)**7*parm**3 *4*pi
		
		!Set Lebedev grids according to shell radius
		!For more inner shell, the distribution is more akin to spherically symmetric, therefore lower number of grids could be used
		if (att2atm(iatt)==0) then !NNA
			sphpotAIM=302
		else
			radtmp=covr(a(att2atm(iatt))%index)
			if (radr<0.2D0*radtmp) then
				sphpotAIM=26
			else if (radr<0.5D0*radtmp) then
				sphpotAIM=74
			else if (radr<0.8D0*radtmp) then
				sphpotAIM=146
			else
				sphpotAIM=194
			end if
		end if
		if (allocated(potx)) deallocate(potx,poty,potz,potw)
		allocate(potx(sphpotAIM),poty(sphpotAIM),potz(sphpotAIM),potw(sphpotAIM))
		call Lebedevgen(sphpotAIM,potx,poty,potz,potw)
		!Combine radial point and weights with angular part, and make them centered at current attractor
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%x=radr*potx+CPpos(1,numcp)
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%y=radr*poty+CPpos(2,numcp)
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%z=radr*potz+CPpos(3,numcp)
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%value=radw*potw
		!Find out trust radius for present attractor
		!If in a shell, the angle between "linking line between nucleus and a shell point" and "gradient vector of this point" &
		!is larger than 45 degree, then this shell is trust radius
		angmax=0
		if (att2atm(iatt)==0) then
			radinit=0
		else
			radrinit=0.15D0
			if (a(att2atm(iatt))%index>2) radrinit=0.5D0
		end if
		if (isettrustrad==0.and.radr>radrinit) then
			do isphpt=1,sphpotAIM
				xtmp=gridatt(nintgrid+isphpt)%x
				ytmp=gridatt(nintgrid+isphpt)%y
				ztmp=gridatt(nintgrid+isphpt)%z
				call calchessmat_dens(1,xtmp,ytmp,ztmp,dens,grad,hess) !Only value and gradient
				dirx=CPpos(1,numcp)-xtmp
				diry=CPpos(2,numcp)-ytmp
				dirz=CPpos(3,numcp)-ztmp
				angtmp=vecang(dirx,diry,dirz,grad(1),grad(2),grad(3))
				if (angtmp>angmax) angmax=angtmp
				if (angtmp>45) then
					write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr')") iatt,trustrad(iatt)
					isettrustrad=1 !The radius of last shell should be the final trust radius. Now exit
					exit
				end if
			end do
			if (isettrustrad==0) trustrad(iatt)=radr !Passed this shell and temporarily set the radius as trust radius. Continue to enlarge the trust radius, until reached angmax>45 degree
		end if
		nintgrid=nintgrid+sphpotAIM
	end do
	if (isettrustrad==0) then !Trust radius was not set after run over all shells
		trustrad(iatt)=1000 !Infinite, for isolated atom
		if (ispecial==2) trustrad(iatt)=20 !For Shubin's 2nd project, should not be as large as 1000, because for a point very far from nucleus the relative entropy cannot be evaluated
		write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr')") iatt,trustrad(iatt)
	end if
	
	!Use DFT integration algorithm to integrate the region inside trust radius
	!$OMP PARALLEL private(ipt,ptx,pty,ptz,rx,ry,rz,dist,tmps,iter,switchwei,intvalp,&
	!$OMP eleintp,xintp,yintp,zintp,xxintp,yyintp,zzintp,xyintp,yzintp,xzintp,tmpval,tmpval2,tmpval3,tmpval4) &
	!$OMP shared(intval,gridval) NUM_THREADS(nthreads)
	intvalp=0D0
	eleintp=0D0
	xintp=0D0
	yintp=0D0
	zintp=0D0
	xxintp=0D0
	yyintp=0D0
	zzintp=0D0
	xyintp=0D0
	yzintp=0D0
	xzintp=0D0
	!$OMP do schedule(DYNAMIC)
	do ipt=1,nintgrid
		ptx=gridatt(ipt)%x
		pty=gridatt(ipt)%y
		ptz=gridatt(ipt)%z
		rx=ptx-CPpos(1,numcp) !The relative distance between current point to corresponding attractor
		ry=pty-CPpos(2,numcp)
		rz=ptz-CPpos(3,numcp)
		!Calculate switching function
		dist=dsqrt(rx*rx+ry*ry+rz*rz)
		tmps=dist-trustrad(iatt)
		if (tmps>1) then
			switchwei=0
		else if (tmps<-1) then
			switchwei=1
		else
			do iter=1,nbeckeiter
				tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
			end do
			switchwei=0.5D0*(1-tmps)
		end if
		gridval(ipt,3)=switchwei
! 		if (dist>trustrad(iatt)) cycle !Discrete separation between atomic-center and uniform integration
		if (switchwei<1D-7) cycle !For saving computational time
		
		if (itype==1.or.itype==2.or.itype==3) then !Integrate a function
			if (ifuncint==-1) then !Deformation density, store molecular density of present center temporarily
				gridval(ipt,1)=fdens(ptx,pty,ptz)
			else if (ispecial==0) then !Normal case
				tmpval=calcfuncall(ifuncint,ptx,pty,ptz)
				intvalp(iatt,1)=intvalp(iatt,1)+gridatt(ipt)%value*tmpval*switchwei
			else if (ispecial==1) then !For Shubin's project, simultaneously output many properties
				tmpval=infoentro(2,ptx,pty,ptz) !Shannon entropy density
				tmpval2=Fisherinfo(1,ptx,pty,ptz) !Fisher information density
				tmpval3=weizsacker(ptx,pty,ptz) !Steric energy
				tmpval4=fdens(ptx,pty,ptz) !Electron density
				intvalp(iatt,1)=intvalp(iatt,1)+gridatt(ipt)%value*tmpval*switchwei
				intvalp(iatt,2)=intvalp(iatt,2)+gridatt(ipt)%value*tmpval2*switchwei
				intvalp(iatt,3)=intvalp(iatt,3)+gridatt(ipt)%value*tmpval3*switchwei
				intvalp(iatt,4)=intvalp(iatt,4)+gridatt(ipt)%value*tmpval4*switchwei
			else if (ispecial==2) then !For Shubin's 2nd project
				call calchessmat_dens(1,ptx,pty,ptz,gridval(ipt,1),gridval(ipt,4:6),hess)
				gridval(ipt,2)=sum(gridval(ipt,4:6)**2)
			end if
		else if (itype==10) then !Calculate multipole moment
			tmpval=gridatt(ipt)%value*fdens(ptx,pty,ptz)*switchwei
			eleintp(iatt)=eleintp(iatt)+tmpval
			xintp(iatt)=xintp(iatt)+rx*tmpval
			yintp(iatt)=yintp(iatt)+ry*tmpval
			zintp(iatt)=zintp(iatt)+rz*tmpval
			xxintp(iatt)=xxintp(iatt)+rx*rx*tmpval
			yyintp(iatt)=yyintp(iatt)+ry*ry*tmpval
			zzintp(iatt)=zzintp(iatt)+rz*rz*tmpval
			xyintp(iatt)=xyintp(iatt)+rx*ry*tmpval
			yzintp(iatt)=yzintp(iatt)+ry*rz*tmpval
			xzintp(iatt)=xzintp(iatt)+rx*rz*tmpval
		end if
	end do
	!$OMP end do
	!$OMP CRITICAL
	if (itype==1.or.itype==2.or.itype==3) then
		intval=intval+intvalp
	else if (itype==10) then
		eleint=eleint+eleintp
		xint=xint+xintp
		yint=yint+yintp
		zint=zint+zintp
		xxint=xxint+xxintp
		yyint=yyint+yyintp
		zzint=zzint+zzintp
		xyint=xyint+xyintp
		yzint=yzint+yzintp
		xzint=xzint+xzintp
	end if
	!$OMP end CRITICAL
	!$OMP END PARALLEL
	
	!Some special cases:
	if (ispecial==2) then !Shubin's 2nd project, integrate relative Shannon and Fisher information
		!Current gridval content: 1=rho, 2=gradrho^2, 3=switchwei, 4=gradrho_x, 5=gradrho_y, 6=gradrho_z
		call dealloall(0)
		write(*,"(' Loading ',a,/)") trim(custommapname(att2atm(iatt)))
		call readwfn(custommapname(att2atm(iatt)),1)
		do ipt=1,nintgrid
			switchwei=gridval(ipt,3)
			if (switchwei<1D-7) cycle
			ptx=gridatt(ipt)%x
			pty=gridatt(ipt)%y
			ptz=gridatt(ipt)%z
			call calchessmat_dens(1,ptx,pty,ptz,prodensgrad(0),prodensgrad(1:3),hess)
			prodens=prodensgrad(0) !rho0_A
			prodensgrad2=sum(prodensgrad(1:3)**2)
			!Relative Shannon entropy, integrate rho*log(rho/rho0_A)]
			tmpval=gridval(ipt,1)*log(gridval(ipt,1)/prodens)
			intval(iatt,1)=intval(iatt,1)+gridatt(ipt)%value*tmpval*switchwei
			!Relative Fisher information (old and incorrect implementation)
			tmpval=gridval(ipt,2)/gridval(ipt,1)-prodensgrad2/prodens
			intval(iatt,2)=intval(iatt,2)+gridatt(ipt)%value*tmpval*switchwei
			!Relative Fisher information (new and correct implementation)
			tmpvalx=gridval(ipt,4)/gridval(ipt,1) - prodensgrad(1)/prodens
			tmpvaly=gridval(ipt,5)/gridval(ipt,1) - prodensgrad(2)/prodens
			tmpvalz=gridval(ipt,6)/gridval(ipt,1) - prodensgrad(3)/prodens
			tmpval=gridval(ipt,1)*(tmpvalx**2+tmpvaly**2+tmpvalz**2)
			intval(iatt,3)=intval(iatt,3)+gridatt(ipt)%value*tmpval*switchwei
		end do
		call dealloall(0)
		call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
	else if (ifuncint==-1) then !Integrate deformation density
		gridval(:,2)=0D0
		do iatm=1,ncenter_org !Cycle each atom to calculate deformation density at all integration grid
			call dealloall(0)
			call readwfn(custommapname(iatm),1)
			do ipt=1,nintgrid
				switchwei=gridval(ipt,3)
				if (switchwei<1D-7) cycle
				gridval(ipt,2)=gridval(ipt,2)+fdens(gridatt(ipt)%x,gridatt(ipt)%y,gridatt(ipt)%z)
			end do
		end do
		do ipt=1,nintgrid !Now gridval(ipt,1/2/3) records actual, promolecular density and weight at ipt
			defdens=gridval(ipt,1)-gridval(ipt,2)
			intval(iatt,1)=intval(iatt,1)+gridatt(ipt)%value*defdens*gridval(ipt,3)
		end do
		call dealloall(0)
		call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
	end if
end do !End cycle attractors

if (itype==1.or.itype==2.or.itype==3) then
	write(*,*) "Integration result inside trust spheres"
	if (ispecial/=2) then
		write(*,*) "  #Sphere       Integral(a.u.)"
		do iatt=1,numrealatt
			write(*,"(i8,f22.10)") iatt,intval(iatt,1)
		end do
		write(*,"(' Sum of above values:',f20.8)") sum(intval(1:numrealatt,1))
	else if (ispecial==2) then !Shubin's 2nd project
		write(*,*) "  #Sphere      Rel.Shannon         Rel.Fisher(old)     Rel.Fisher(new)"
		do iatt=1,numrealatt
			write(*,"(i8,3f20.10)") iatt,intval(iatt,1:3)
		end do
		write(*,"(' Sum of rel.Shannon:',f20.8)") sum(intval(1:numrealatt,1))
		write(*,"(' Sum of rel.Fisher(old): ',f20.8)") sum(intval(1:numrealatt,2))
		write(*,"(' Sum of rel.Fisher(new): ',f20.8)") sum(intval(1:numrealatt,3))
	end if
end if

!Set coordinate of uniform grids
call calc_dvol(dvol)
do ix=1,nx
	xarr(ix)=orgx+(ix-1)*dx
end do
do iy=1,ny
	yarr(iy)=orgy+(iy-1)*dy
end do
do iz=1,nz
	zarr(iz)=orgz+(iz-1)*dz
end do

!--------- Integrating uniform grids, basin boundary grids will be calculated in the later stage
!NOTE: For shubin's 2nd project or deformation density, do not integrate here. At next stage boundary grids will be updated, and then at the next stage,&
!they will be integrated by a special module. Because at each grid if we reload atomic .wfn file will be too time consuming
if (ispecial==2.or.ifuncint==-1) goto 10 
write(*,*)
write(*,*) "Integrating uniform grids..."
ifinish=0
!$OMP PARALLEL private(ix,iy,iz,iatt,icp,rnowx,rnowy,rnowz,rx,ry,rz,tmpval,tmpval2,tmpval3,intvalp,basinvolp,basinvdwvolp,dist,tmps,iter,switchwei,&
!$OMP eleintp,xintp,yintp,zintp,xxintp,yyintp,zzintp,xyintp,xzintp,yzintp) shared(intval,basinvol,basinvdwvol,ifinish) NUM_THREADS(nthreads)
intvalp=0D0
basinvolp=0D0
basinvdwvolp=0D0
eleintp=0D0
xintp=0D0
yintp=0D0
zintp=0D0
xxintp=0D0
yyintp=0D0
zzintp=0D0
xyintp=0D0
yzintp=0D0
xzintp=0D0
!$OMP do schedule(DYNAMIC)
do iz=izlow,izup
	do iy=iylow,iyup
		do ix=ixlow,ixup
			if ((itype==2.or.itype==3).and.interbasgrid(ix,iy,iz)) cycle !If refine boundary grid at next stage, we don't calculate them at present stage
            call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
			iatt=gridbas(ix,iy,iz)
! 			do icp=1,numcp
! 				dist=dsqrt( (rnowx-CPpos(1,icp))**2+(rnowy-CPpos(2,icp))**2+(rnowz-CPpos(3,icp))**2 )
! 				if (disttest<trustrad(icp)) cycle cycix !The function inside trust radius is integrated by DFT integration
! 			end do
			rx=rnowx-CPpos(1,iatt) !The relative distance between current point to exact position of corresponding attractor
			ry=rnowy-CPpos(2,iatt)
			rz=rnowz-CPpos(3,iatt)
			!Calculate switching function at current grid
			dist=dsqrt(rx*rx+ry*ry+rz*rz)
			tmps=dist-trustrad(iatt)
			if (tmps>1) then
				switchwei=0
			else if (tmps<-1) then
				switchwei=1
			else
				do iter=1,nbeckeiter
					tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
				end do
				switchwei=0.5D0*(1-tmps)
			end if
			switchwei=1-switchwei
			basinvolp(iatt)=basinvolp(iatt)+1 !Calculate basin volume
			if (cubmat(ix,iy,iz)>0.001D0) basinvdwvolp(iatt)=basinvdwvolp(iatt)+1
			if (switchwei<1D-7) cycle !For saving time
			if (itype==1.or.itype==2.or.itype==3) then
				if (ispecial==0) then
					if (ifuncint==1) then !Electron density on each grid has already been calculated
						tmpval=cubmat(ix,iy,iz)
					else
						tmpval=calcfuncall(ifuncint,rnowx,rnowy,rnowz)
					end if
					intvalp(iatt,1)=intvalp(iatt,1)+tmpval*switchwei
				else if (ispecial==1) then !For Chunying Rong
					tmpval=infoentro(2,rnowx,rnowy,rnowz) !Shannon entropy density, see JCP,126,191107 for example
					tmpval2=Fisherinfo(1,rnowx,rnowy,rnowz) !Fisher information density, see JCP,126,191107 for example
					tmpval3=weizsacker(rnowx,rnowy,rnowz) !Steric energy
					tmpval4=fdens(rnowx,rnowy,rnowz) !Electron density
					intvalp(iatt,1)=intvalp(iatt,1)+tmpval*switchwei
					intvalp(iatt,2)=intvalp(iatt,2)+tmpval2*switchwei
					intvalp(iatt,3)=intvalp(iatt,3)+tmpval3*switchwei
					intvalp(iatt,4)=intvalp(iatt,4)+tmpval4*switchwei
				end if
			else if (itype==10) then
				tmpval=cubmat(ix,iy,iz)*switchwei
				eleintp(iatt)=eleintp(iatt)+tmpval
				xintp(iatt)=xintp(iatt)+rx*tmpval
				yintp(iatt)=yintp(iatt)+ry*tmpval
				zintp(iatt)=zintp(iatt)+rz*tmpval
				xxintp(iatt)=xxintp(iatt)+rx*rx*tmpval
				yyintp(iatt)=yyintp(iatt)+ry*ry*tmpval
				zzintp(iatt)=zzintp(iatt)+rz*rz*tmpval
				xyintp(iatt)=xyintp(iatt)+rx*ry*tmpval
				yzintp(iatt)=yzintp(iatt)+ry*rz*tmpval
				xzintp(iatt)=xzintp(iatt)+rx*rz*tmpval
			end if
		end do
	end do
	!$OMP CRITICAL
    ifinish=ifinish+1
    if (ifuncint/=1) call showprog(ifinish,nz-2)
	!$OMP end CRITICAL
end do
!$OMP end do
!$OMP CRITICAL
if (itype==1.or.itype==2.or.itype==3) then
	intval=intval+intvalp*dvol
	basinvol=basinvol+basinvolp*dvol
	basinvdwvol=basinvdwvol+basinvdwvolp*dvol
else if (itype==10) then
	eleint=eleint+eleintp*dvol
	xint=xint+xintp*dvol
	yint=yint+yintp*dvol
	zint=zint+zintp*dvol
	xxint=xxint+xxintp*dvol
	yyint=yyint+yyintp*dvol
	zzint=zzint+zzintp*dvol
	xyint=xyint+xyintp*dvol
	yzint=yzint+yzintp*dvol
	xzint=xzint+xzintp*dvol
end if
!$OMP end CRITICAL
!$OMP END PARALLEL

10 continue
!---- Exact refinement with/without multi-level splitting of boundary grids
if (itype==2.or.itype==3) then
	if (itype==3) then !Calculate grid data of gradient of electron density used to linear interpolation to obtain the value at any point
		call gengradmat
	end if
	write(*,*) "Integrating grids at basin boundary..." 
	nrk4lim=100
	nrk4gradswitch=40
	hsizeinit=0.25D0
	ifinish=0
    gridd1=dsqrt(sum(gridv1**2))
    gridd2=dsqrt(sum(gridv2**2))
    gridd3=dsqrt(sum(gridv3**2))
	!$OMP PARALLEL private(ix,iy,iz,iatt,rnowx,rnowy,rnowz,rx,ry,rz,tmpval,tmpval2,tmpval3,intvalp,basinvolp,basinvdwvolp,dist,tmps,iter,switchwei,&
	!$OMP rnowxtmp,rnowytmp,rnowztmp,orgxref,orgyref,orgzref,dxref,dyref,dzref,ixref,iyref,izref,nrefine,ndiv,&
	!$OMP k1,k2,k3,k4,dens,denshold,grad,hess,iattref,xtmp,ytmp,ztmp,irk4,hsize,ixtest,iytest,iztest,tmpdist) shared(intval,basinvol,basinvdwvol,ifinish) NUM_THREADS(nthreads)
	intvalp=0D0
	basinvolp=0D0
	basinvdwvolp=0D0
	!$OMP do schedule(DYNAMIC)
	do iz=izlow,izup
		do iy=iylow,iyup
			do ix=ixlow,ixup
				if (.not.interbasgrid(ix,iy,iz)) cycle
                call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
! 				if (cubmat(ix,iy,iz)>0.001D0) then
! 					nrefine=2 !3 is the best
! 				else if (cubmat(ix,iy,iz)>0.001D0) then !0.0001 is the best
! 					nrefine=1
! 				else
! 					nrefine=1
! 				end if
 				nrefine=1
				ndiv=nrefine**3
				orgxref=rnowx-dx/2 !Take corner position as original point of microcycle
				orgyref=rnowy-dy/2
				orgzref=rnowz-dz/2
				dxref=dx/nrefine
				dyref=dy/nrefine
				dzref=dz/nrefine
				do ixref=1,nrefine
					do iyref=1,nrefine
						do izref=1,nrefine
							rnowxtmp=orgxref+(ixref-0.5D0)*dxref !Coordinate of current refined grid
							rnowytmp=orgyref+(iyref-0.5D0)*dyref
							rnowztmp=orgzref+(izref-0.5D0)*dzref
							if (cubmat(ix,iy,iz)<=0.001D0) then !Only refine the boundary inside vdW surface
								iattref=gridbas(ix,iy,iz)
							else
								xtmp=rnowxtmp !This point will continuously move in the iteration
								ytmp=rnowytmp
								ztmp=rnowztmp
								hsize=hsizeinit
								densold=0D0
								!** Tracing steepest ascent trajectory using 4-order Runge-Kutta (RK4)
		cycrk4:					do irk4=1,nrk4lim
									!For full accuracy refinement, or the first step, or when interpolation gradient works worse,&
									!namely has not converge until nrk4gradswitch, use exactly evaluated gradient
									if (itype==2.or.irk4==1.or.irk4==2.or.irk4>nrk4gradswitch) then 
										if (itype==3.and.irk4==nrk4gradswitch+1) then !Interpolated gradient doesn't work well, switch to full accuracy, reset the coordinate
											xtmp=rnowxtmp
											ytmp=rnowytmp
											ztmp=rnowztmp
											hsize=hsizeinit
										end if
										call calchessmat_dens(1,xtmp,ytmp,ztmp,dens,grad,hess) !Only value and gradient
										if (dens<densold-1D-10) then
											hsize=hsize*0.75D0 !Reduce step size if density decrease
										else if (dens>densold+1D-10) then
											hsize=hsizeinit !Recover to initial step size
										end if
										denshold=dens
										k1=grad/dsqrt(sum(grad**2))
										call calchessmat_dens(1,xtmp+hsize/2*k1(1),ytmp+hsize/2*k1(2),ztmp+hsize/2*k1(3),dens,grad,hess) !Only value and gradient
										k2=grad/dsqrt(sum(grad**2))
										call calchessmat_dens(1,xtmp+hsize/2*k2(1),ytmp+hsize/2*k2(2),ztmp+hsize/2*k2(3),dens,grad,hess) !Only value and gradient
										k3=grad/dsqrt(sum(grad**2))
										call calchessmat_dens(1,xtmp+hsize*k3(1),ytmp+hsize*k3(2),ztmp+hsize*k3(3),dens,grad,hess) !Only value and gradient
										k4=grad/dsqrt(sum(grad**2))
									else !Using the gradients evaluated by trilinear interpolation from pre-calculated grid data to save computational time
										call linintp3dvec(xtmp,ytmp,ztmp,grad) !Only value and gradient
! 										if (dens<densold-1D-10) then
! 											hsize=hsize*0.75D0
! 										else if (dens>densold+1D-10) then
! 											hsize=hsizeinit !Recover to initial step size
! 										end if
! 										denshold=dens
										k1=grad/dsqrt(sum(grad**2))
										call linintp3dvec(xtmp+hsize/2*k1(1),ytmp+hsize/2*k1(2),ztmp+hsize/2*k1(3),grad)
										k2=grad/dsqrt(sum(grad**2))
										call linintp3dvec(xtmp+hsize/2*k2(1),ytmp+hsize/2*k2(2),ztmp+hsize/2*k2(3),grad)
										k3=grad/dsqrt(sum(grad**2))
										call linintp3dvec(xtmp+hsize*k3(1),ytmp+hsize*k3(2),ztmp+hsize*k3(3),grad)
										k4=grad/dsqrt(sum(grad**2))
									end if
									xtmp=xtmp+hsize/6*(k1(1)+2*k2(1)+2*k3(1)+k4(1)) !Update current coordinate
									ytmp=ytmp+hsize/6*(k1(2)+2*k2(2)+2*k3(2)+k4(2))
									ztmp=ztmp+hsize/6*(k1(3)+2*k2(3)+2*k3(3)+k4(3))
									!Check if current position has entered trust radius of an attractor
									do iatttmp=1,numrealatt
										dist=dsqrt( (xtmp-CPpos(1,iatttmp))**2+(ytmp-CPpos(2,iatttmp))**2+(ztmp-CPpos(3,iatttmp))**2 )
										if (dist<trustrad(iatttmp)) then
											iattref=iatttmp
											exit cycrk4
										end if
									end do
									!Check if the closest grid and its 26 neighbours have the same attribution, if yes, employ its attribution then exit
                                    !Find the closest grid to current position as (ixtest,iytest,iztest)
									do ixtest=ixlow,ixup !Find closest grid to current in direction 1
										tmpdist=abs(xtmp-xarr(ixtest))
										if (tmpdist<1.01D0*gridd1/2D0) exit
									end do
									do iytest=iylow,iyup !Find closest grid to current in direction 2
										tmpdist=abs(ytmp-yarr(iytest))
										if (tmpdist<1.01D0*gridd2/2D0) exit
									end do
									do iztest=izlow,izup !Find closest grid to current in direction 3
										tmpdist=abs(ztmp-zarr(iztest))
										if (tmpdist<1.01D0*gridd3/2D0) exit
									end do
                                    if (ixtest==nx.or.iytest==ny.or.iztest==nz) cycle !Current position is far away from any grid, usually the initial case
									iattref=gridbas(ixtest,iytest,iztest)
									do imove=1,26
										if ( gridbas(ixtest+vec26x(imove),iytest+vec26y(imove),iztest+vec26z(imove))/=iattref ) exit
									end do
									if (imove==27) exit !Successfully passed neighbour test
								end do cycrk4
								if (irk4==nrk4lim+1) then !Didn't enter trust radius or didn't approach a grid who and whose neighbour have the same attribution
! 									write(*,*) "Warning: Exceeded the step limit of steepest ascent process!" !may frighten users, so comment
									iattref=gridbas(ix,iy,iz) !Use its original attribution
								end if
							end if
							gridbas(ix,iy,iz)=iattref !Update attribution of boundary grids
							!Calculate switching function at current grid
							rx=rnowxtmp-CPpos(1,iattref) !The relative distance between current point to corresponding attractor
							ry=rnowytmp-CPpos(2,iattref)
							rz=rnowztmp-CPpos(3,iattref)
							dist=dsqrt(rx*rx+ry*ry+rz*rz)
							tmps=dist-trustrad(iattref)
							if (tmps>1) then
								switchwei=0
							else if (tmps<-1) then
								switchwei=1
							else
								do iter=1,nbeckeiter
									tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
								end do
								switchwei=0.5D0*(1-tmps)
							end if
							switchwei=1-switchwei
							basinvolp(iattref)=basinvolp(iattref)+1D0/ndiv !Calculate boundary basin volume
							if (cubmat(ix,iy,iz)>0.001D0) basinvdwvolp(iattref)=basinvdwvolp(iattref)+1D0/ndiv
							if (ispecial==2.or.ifuncint==-1) then !Only for Shubin
								continue !Don't calculate function value, but only update attribution of boundary grids at this stage
							else if (ispecial==0) then
								tmpval=calcfuncall(ifuncint,rnowxtmp,rnowytmp,rnowztmp)
								intvalp(iattref,1)=intvalp(iattref,1)+tmpval*switchwei/ndiv
							else if (ispecial==1) then !Only for Shubin
								tmpval=infoentro(2,rnowxtmp,rnowytmp,rnowztmp) !Shannon entropy density, see JCP,126,191107 for example
								tmpval2=Fisherinfo(1,rnowxtmp,rnowytmp,rnowztmp) !Fisher information density, see JCP,126,191107 for example
								tmpval3=weizsacker(rnowxtmp,rnowytmp,rnowztmp) !Steric energy
								tmpval4=fdens(rnowxtmp,rnowytmp,rnowztmp) !Electron density
								intvalp(iattref,1)=intvalp(iattref,1)+tmpval*switchwei/ndiv
								intvalp(iattref,2)=intvalp(iattref,2)+tmpval2*switchwei/ndiv
								intvalp(iattref,3)=intvalp(iattref,3)+tmpval3*switchwei/ndiv
								intvalp(iattref,4)=intvalp(iattref,4)+tmpval4*switchwei/ndiv
							end if
						end do !End refine grid
					end do
				end do
				
			end do !End cycle ix grid
		end do
		!$OMP CRITICAL
		ifinish=ifinish+1
		call showprog(ifinish,nz-2)
		!$OMP end CRITICAL
	end do
	!$OMP end do
	!$OMP CRITICAL
	intval=intval+intvalp*dvol
	basinvol=basinvol+basinvolp*dvol
	basinvdwvol=basinvdwvol+basinvdwvolp*dvol
	!$OMP end CRITICAL
	!$OMP END PARALLEL
	call detectinterbasgrd(6)
	write(*,*) "Basin boundary has been updated"
	numinterbas=count(interbasgrid.eqv..true.)
	write(*,"(' The number of interbasin grids:',i12)") numinterbas	
end if

!Below are special modules for the cases when density of atom in free-state are involved
if (ispecial==2) then !Shubin's 2nd project, integrate relative Shannon and relative Fisher information
	allocate(rhogrid(nx,ny,nz),rhogradgrid(3,nx,ny,nz))
	ifinish=0
	write(*,*)
	write(*,*) "Calculating electron density and its gradient for actual system at each grid"
	!$OMP PARALLEL DO SHARED(rhogrid,rhogradgrid,ifinish) PRIVATE(ix,iy,iz,ptx,pty,ptz) schedule(dynamic) NUM_THREADS(nthreads)
	do iz=izlow,izup
		do iy=iylow,iyup
			do ix=ixlow,ixup
				call getgridxyz(ix,iy,iz,ptx,pty,ptz)
				call calchessmat_dens(1,ptx,pty,ptz,rhogrid(ix,iy,iz),rhogradgrid(:,ix,iy,iz),hess)
			end do
		end do
		!$OMP CRITICAL
		ifinish=ifinish+1
		write(*,"(' Finished',i6,' /',i6)") ifinish,izup-izlow+1
		!$OMP end CRITICAL
	end do
	!$OMP end PARALLEL DO
	write(*,*)
	write(*,*) "Calculating electron density and its gradient for free-state atom at each grid"
	do iatt=1,numrealatt !Cycle each attractors
		write(*,"(' Processing ',a)") trim(custommapname(att2atm(iatt)))
		call dealloall(0)
		call readwfn(custommapname(att2atm(iatt)),1)
		!$OMP PARALLEL private(intvalp,ix,iy,iz,ptx,pty,ptz,rx,ry,rz,dist,tmps,switchwei,&
		!$OMP rho,rhograd2,prodens,prodensgrad,prodensgrad2,tmpval1,tmpval2,tmpval3,tmpx,tmpy,tmpz) shared(intval) NUM_THREADS(nthreads)
		intvalp=0D0
		!$OMP do schedule(DYNAMIC)
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					if (gridbas(ix,iy,iz)==iatt) then
						!Calculate switching function at current grid
						call getgridxyz(ix,iy,iz,ptx,pty,ptz)
						rx=ptx-CPpos(1,iatt) !The relative distance between current point to corresponding attractor
						ry=pty-CPpos(2,iatt)
						rz=ptz-CPpos(3,iatt)
						dist=dsqrt(rx*rx+ry*ry+rz*rz)
						tmps=dist-trustrad(iatt)
						if (tmps>1) then
							switchwei=0
						else if (tmps<-1) then
							switchwei=1
						else
							do iter=1,nbeckeiter
								tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
							end do
							switchwei=0.5D0*(1-tmps)
						end if
						switchwei=1-switchwei
						rho=rhogrid(ix,iy,iz)
						rhograd2=sum(rhogradgrid(1:3,ix,iy,iz)**2)
						call calchessmat_dens(1,ptx,pty,ptz,prodens,prodensgrad(1:3),hess)
						prodensgrad2=sum(prodensgrad(1:3)**2)
						!Relative Shannon entropy
						tmpval1=rho*log(rho/prodens)
						intvalp(iatt,1)=intvalp(iatt,1)+tmpval1*switchwei
						!Relative Fisher information (old and incorrect form)
						tmpval2=rhograd2/rho-prodensgrad2/prodens
						intvalp(iatt,2)=intvalp(iatt,2)+tmpval2*switchwei
						!Relative Fisher information (new and correct form)
						tmpx=rhogradgrid(1,ix,iy,iz)/rho-prodensgrad(1)/prodens
						tmpy=rhogradgrid(2,ix,iy,iz)/rho-prodensgrad(2)/prodens
						tmpz=rhogradgrid(3,ix,iy,iz)/rho-prodensgrad(3)/prodens
						tmpval3=(tmpx**2+tmpy**2+tmpz**2)*rho
						intvalp(iatt,3)=intvalp(iatt,3)+tmpval3*switchwei
					end if
				end do
			end do
		end do
		!$OMP end do
		!$OMP CRITICAL
		intval=intval+intvalp*dvol
		!$OMP end CRITICAL
		!$OMP END PARALLEL
	end do
	deallocate(rhogrid,rhogradgrid)
	call dealloall(0)
	write(*,"(' Reloading ',a)") trim(firstfilename)
	call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule)
else if (ifuncint==-1) then !Deformation density
	allocate(prorhogrid(nx,ny,nz))
	write(*,*)
	write(*,*) "Calculating promolecular density at each grid"
	prorhogrid=0D0
	do iatm=1,ncenter_org !Cycle each atom
		write(*,"(' Processing atom',i6,a,'...')") iatm,a_org(iatm)%name
		call dealloall(0)
		call readwfn(custommapname(iatm),1)
		!$OMP PARALLEL DO SHARED(prorhogrid) PRIVATE(tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads)
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					call getgridxyz(ix,iy,iz,tmpx,tmpy,tmpz)
					prorhogrid(ix,iy,iz)=prorhogrid(ix,iy,iz)+fdens(tmpx,tmpy,tmpz)
				end do
			end do
		end do
		!$OMP end PARALLEL DO
	end do
	do iatt=1,numrealatt !Cycle each attractors
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					if (gridbas(ix,iy,iz)==iatt) then
						!Calculate switching function at current grid
						call getgridxyz(ix,iy,iz,tmpx,tmpy,tmpz)
						rx=tmpx-CPpos(1,iatt) !Relative distance between current point to corresponding attractor
						ry=tmpy-CPpos(2,iatt)
						rz=tmpz-CPpos(3,iatt)
						dist=dsqrt(rx*rx+ry*ry+rz*rz)
						tmps=dist-trustrad(iatt)
						if (tmps>1) then
							switchwei=0
						else if (tmps<-1) then
							switchwei=1
						else
							do iter=1,nbeckeiter
								tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
							end do
							switchwei=0.5D0*(1-tmps)
						end if
						switchwei=1-switchwei
						defdens=cubmat(ix,iy,iz)-prorhogrid(ix,iy,iz)
						intval(iatt,1)=intval(iatt,1)+defdens*switchwei*dvol
					end if
				end do
			end do
		end do
	end do
	deallocate(prorhogrid)
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule) to calc real rho again
end if

!!----------- Output SUMMARY
if (itype==1.or.itype==2.or.itype==3) then !Integrate specific real space function(s)
	write(*,*)
	if (ifuncint==1) then
		write(*,*) "Total result:"
		write(*,*) "  #Basin        Integral(a.u.)      Vol(Bohr^3)    Vol(rho>0.001)"
		do iatt=1,numrealatt
			write(*,"(i8,f22.10,2f16.3)") iatt,intval(iatt,1),basinvol(iatt),basinvdwvol(iatt)
		end do
		write(*,"(' Sum of above integrals:',f20.8)") sum(intval(1:numrealatt,1))
		write(*,"(' Sum of basin volumes (rho>0.001):',f12.3,' Bohr^3')") sum(basinvdwvol(1:numrealatt))
		if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0,1)
		if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1,1)
		write(*,*)
		rnormfac=sum(intval(1:numrealatt,1))/(nelec+nEDFelec) !The electrons represented by EDF must be taken into account!
		write(*,"(' Normalization factor of the integral of electron density is',f12.6)") rnormfac
		write(*,*) "The atomic charges after normalization and atomic volumes:"
		do iatm=0,ncenter
			do iatt=1,numrealatt
				if (att2atm(iatt)==iatm.and.iatm==0) then
					write(*,"(i7,' (NNA)   Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatt,-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
				else if (att2atm(iatt)==iatm.and.iatm/=0) then
				    if (nEDFelec==0) then !Normal case, all electron basis or using pseudopotential but not accompanied by EDF
					    write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%charge-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
					else !EDF is used, so using a(iatm)%index instead of a(iatm)%charge
					    write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%index-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
					end if
				end if
			end do
		end do
	else
		if (ispecial==0) then
			write(*,*) "Total result:"
			write(*,*) "    Atom       Basin       Integral(a.u.)   Vol(Bohr^3)   Vol(rho>0.001)"
			do iatm=0,ncenter
				do iatt=1,numrealatt
					if (att2atm(iatt)==iatm.and.iatm==0) then
						write(*,"('      NNA   ',i8,f20.8,2f14.3)") iatt,intval(iatt,1),basinvol(iatt),basinvdwvol(iatt)
					else if (att2atm(iatt)==iatm.and.iatm/=0) then
						write(*,"(i7,' (',a,')',i8,f20.8,2f14.3)") iatm,a(iatm)%name,iatt,intval(iatt,1),basinvol(iatt),basinvdwvol(iatt)
					end if
				end do
			end do
			write(*,"(' Sum of above integrals:',f23.8)") sum(intval(1:numrealatt,1))
			write(*,"(' Sum of basin volumes (rho>0.001):',f12.3,' Bohr^3')") sum(basinvdwvol(1:numrealatt))
			if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==0)) write(*,"(' Integral of unassigned grids:',f20.8)") intval(0,1)
			if (any(gridbas(ixlow:ixup,iylow:iyup,izlow:izup)==-1)) write(*,"(' Integral of the grids travelled to box boundary:',f20.8)") intval(-1,1)
		else if (ispecial==1) then
			write(*,*) "Total result:"
			write(*,*) "    Atom       Basin     Shannon       Fisher      Steric ene  Vol(rho>0.001)"
			do iatm=0,ncenter
				do iatt=1,numrealatt
					if (att2atm(iatt)==iatm.and.iatm==0) then
						write(*,"('      NNA   ',i8,3f14.7,f13.5)") iatt,intval(iatt,1),intval(iatt,1:3),basinvdwvol(iatt)
					else if (att2atm(iatt)==iatm.and.iatm/=0) then
						write(*,"(i7,' (',a,')',i8,3f14.7,f13.5)") iatm,a(iatm)%name,iatt,intval(iatt,1:3),basinvdwvol(iatt)
					end if
				end do
			end do
			write(*,"(' Total Shannon entropy:   ',f23.8)") sum(intval(1:numrealatt,1))
			write(*,"(' Total Fisher information:',f23.8)") sum(intval(1:numrealatt,2))
			write(*,"(' Total steric energy:     ',f23.8)") sum(intval(1:numrealatt,3))
			write(*,"(' Sum of basin volumes (rho>0.001):',f12.3,' Bohr^3')") sum(basinvdwvol(1:numrealatt))
			rnormfac=sum(intval(1:numrealatt,4))/(nelec+nEDFelec) !The electrons represented by EDF must be taken into account!
			write(*,"(/,' Normalization factor of the integral of electron density is',f12.6)") rnormfac
			write(*,*) "The atomic charges after normalization and atomic volumes:"
			do iatm=0,ncenter
				do iatt=1,numrealatt
					if (att2atm(iatt)==iatm.and.iatm==0) then
						write(*,"(i7,' (NNA)   Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatt,-intval(iatt,1)/rnormfac,basinvdwvol(iatt)
					else if (att2atm(iatt)==iatm.and.iatm/=0) then
						if (nEDFelec==0) then !Normal case, all electron basis or using pseudopotential but not accompanied by EDF
							write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%charge-intval(iatt,4)/rnormfac,basinvdwvol(iatt)
						else !EDF is used, so using a(iatm)%index instead of a(iatm)%charge
							write(*,"(i7,' (',a,')    Charge:',f12.6,'     Volume:',f10.3,' Bohr^3')") iatm,a(iatm)%name,a(iatm)%index-intval(iatt,4)/rnormfac,basinvdwvol(iatt)
						end if
					end if
				end do
			end do
		else if (ispecial==2) then
			write(*,*) "Total result:"
			write(*,*) "    Atom      Basin        Rel.Shannon     Rel.Fisher(old)     Rel.Fisher(new)"
			do iatm=1,ncenter
				do iatt=1,numrealatt
					if (att2atm(iatt)==iatm) write(*,"(i7,' (',a,')',i7,3f20.8)") iatm,a(iatm)%name,iatt,intval(iatt,1:3)
				end do
			end do
			write(*,"(' Sum of relat_Shannon:     ',f23.8)") sum(intval(1:numrealatt,1))
			write(*,"(' Sum of relat_Fisher(old): ',f23.8)") sum(intval(1:numrealatt,2))
			write(*,"(' Sum of relat_Fisher(new): ',f23.8)") sum(intval(1:numrealatt,3))
		end if
	end if
	write(*,*)
else if (itype==10) then !Electric multipole moment
	ioutid=6
101	eleinttot=0D0
	dipelextot=0D0
	dipeleytot=0D0
	dipeleztot=0D0
	if (ioutid==6) write(*,*)
	write(ioutid,*) "Note: All data shown below are in a.u.!"
	write(ioutid,*)
	do iatm=0,ncenter
		do iatt=1,numrealatt
			if (att2atm(iatt)/=iatm) cycle
			if (iatm==0) then
				write(ioutid,"(' *****  Result of NNA',i6)") iatt
			else
				write(ioutid,"(' *****  Result of atom',i6,' (',a,'), corresponding to basin',i6)") iatm,a(iatm)%name,iatt
			end if
			write(ioutid,"(' Basin monopole moments (from electrons):',f12.6)") -eleint(iatt)
			if (iatm/=0) then
				if (nEDFelec==0) then !Normal case, all electron basis or using pseudopotential but not accompanied by EDF
					atmchgtmp=a(iatm)%charge-eleint(iatt)
				else !EDF is used, so using a(iatm)%index instead of a(iatm)%charge
					atmchgtmp=a(iatm)%index-eleint(iatt)
				end if
				write(ioutid,"(' Atomic charge:',f12.6)") atmchgtmp
            end if
			write(ioutid,"(' Basin dipole moments:')") 
			write(ioutid,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") -xint(iatt),-yint(iatt),-zint(iatt),sqrt(xint(iatt)**2+yint(iatt)**2+zint(iatt)**2)
			eleinttot=eleinttot+eleint(iatt)
			dipelex=-eleint(iatt)*CPpos(1,iatt)+(-xint(iatt)) !Contribution to molecular total dipole moment
			dipeley=-eleint(iatt)*CPpos(2,iatt)+(-yint(iatt))
			dipelez=-eleint(iatt)*CPpos(3,iatt)+(-zint(iatt))
			dipelextot=dipelextot+dipelex
			dipeleytot=dipeleytot+dipeley
			dipeleztot=dipeleztot+dipelez
			write(ioutid,"(' Basin electron contribution to molecular dipole moment:')")
			write(ioutid,"(' X=',f12.6,'  Y=',f12.6,'  Z=',f12.6,'  Norm=',f12.6)") dipelex,dipeley,dipelez,sqrt(dipelex**2+dipeley**2+dipelez**2)
			write(ioutid,"(' Basin quadrupole moments (Traceless Cartesian form):')")
			rrint=xxint(iatt)+yyint(iatt)+zzint(iatt)
			QXX=-(3*xxint(iatt)-rrint)/2
			QYY=-(3*yyint(iatt)-rrint)/2
			QZZ=-(3*zzint(iatt)-rrint)/2
			quadmom(1,1)=QXX
			quadmom(1,2)=(-3*xyint(iatt))/2
			quadmom(1,3)=(-3*xzint(iatt))/2
			quadmom(2,1)=(-3*xyint(iatt))/2
			quadmom(2,2)=QYY
			quadmom(2,3)=(-3*yzint(iatt))/2
			quadmom(3,1)=(-3*xzint(iatt))/2
			quadmom(3,2)=(-3*yzint(iatt))/2
			quadmom(3,3)=QZZ
			write(ioutid,"(' XX=',f12.6,'  XY=',f12.6,'  XZ=',f12.6)") quadmom(1,:)
			write(ioutid,"(' YX=',f12.6,'  YY=',f12.6,'  YZ=',f12.6)") quadmom(2,:)
			write(ioutid,"(' ZX=',f12.6,'  ZY=',f12.6,'  ZZ=',f12.6)") quadmom(3,:)
			write(ioutid,"(' Magnitude of the traceless quadrupole moment tensor:',f12.6)") sqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
			R20=-(3*zzint(iatt)-rrint)/2D0 !Notice the negative sign, because electrons carry negative charge
			R2n1=-dsqrt(3D0)*yzint(iatt)
			R2p1=-dsqrt(3D0)*xzint(iatt)
			R2n2=-dsqrt(3D0)*xyint(iatt)
			R2p2=-dsqrt(3D0)/2D0*(xxint(iatt)-yyint(iatt))
			write(ioutid,"(' Basin quadrupole moments (Spherical harmonic form):')")
			write(ioutid,"(' Q_2,0 =',f11.6,'   Q_2,-1=',f11.6,'   Q_2,1=',f11.6)") R20,R2n1,R2p1
			write(ioutid,"(' Q_2,-2=',f11.6,'   Q_2,2 =',f11.6)") R2n2,R2p2
			write(ioutid,"(' Magnitude: |Q_2|=',f12.6)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)
			ESEx=xxint(iatt);ESEy=yyint(iatt);ESEz=zzint(iatt)
			ESE=ESEx+ESEy+ESEz
        		write(ioutid,"(a,f16.6)") " Basin electronic spatial extent <r^2>:",ESE
			write(ioutid,"(' Components of <r^2>:  X=',f15.6,'  Y=',f15.6,'  Z=',f15.6)") ESEx,ESEy,ESEz
			write(ioutid,*)
			!Output dipole and quadrupole moments
            if (ioutid==10.and.iatm/=0) then
				write(20,"(' Atom',i6,' (',a,')')") iatm,a(iatm)%name
				write(20,"(' Atomic dipole moment:',3f12.6)") -xint(iatt),-yint(iatt),-zint(iatt)
				call diagsymat(quadmom,tmpmat,tmpvec,istat)
				write(20,"(' Information of atomic quadrupole moment (Traceless Cartesian):')")
				write(20,"(' Eigenvalue 1:',f12.6,'  Eigenvector:',3f12.6)") tmpvec(1),tmpmat(:,1)
				write(20,"(' Eigenvalue 2:',f12.6,'  Eigenvector:',3f12.6)") tmpvec(2),tmpmat(:,2)
				write(20,"(' Eigenvalue 3:',f12.6,'  Eigenvector:',3f12.6)") tmpvec(3),tmpmat(:,3)
				write(20,*)
            end if
		end do
	end do
	!Output overall electric properties
	dipnucx=sum(a(:)%x*a(:)%charge)
	dipnucy=sum(a(:)%y*a(:)%charge)
	dipnucz=sum(a(:)%z*a(:)%charge)
	write(ioutid,"( ' Molecular net charge:',f12.6)") sum(a%charge)-eleinttot
	write(ioutid,"( ' Nuclear contribution to molecular dipole moment:')") 
	write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipnucx,dipnucy,dipnucz,sqrt(dipnucx**2+dipnucy**2+dipnucz**2)
	write(ioutid,"( ' Electron contribution to molecular dipole moment:')") 
	write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipelextot,dipeleytot,dipeleztot,sqrt(dipelextot**2+dipeleytot**2+dipeleztot**2)
	dipmolx=dipnucx+dipelextot
	dipmoly=dipnucy+dipeleytot
	dipmolz=dipnucz+dipeleztot
	write(ioutid,"( ' Molecular dipole moment:')")
	write(ioutid,"(' X=',f12.6,' Y=',f12.6,' Z=',f12.6,' Norm=',f12.6)") dipmolx,dipmoly,dipmolz,sqrt(dipmolx**2+dipmoly**2+dipmolz**2)
	if (ioutid==10) then
		close(10)
        close(20)
		write(*,*) "Outputting finished!"
		return
	end if
end if

call walltime(iwalltime2)
write(*,"(/,' Integrating basins took up wall clock time',i10,' s')") iwalltime2-iwalltime1

if (itype==10) then
	write(*,*)
	write(*,*) "If also outputting the above result to multipole.txt in current folder? (y/n)"
    write(*,"(a)") " If choose ""y"", atom_moment.txt will also be outputted in current folder, &
    &which contains electric dipole moments as well as eigenvalues and eigenvectors of quadrupole moment tensors of all atoms"
	read(*,*) selectyn
	if (selectyn=='y'.or.selectyn=='Y') then
		ioutid=10
		open(10,file="multipole.txt",status="replace")
        open(20,file="atom_moment.txt",status="replace")
		goto 101
	end if
end if
end subroutine





!!------ Generate grid data of gradient of electron density, used to refine basin boundary
subroutine gengradmat
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8 wfnval(nmo),wfnderv(3,nmo),gradrho(3),EDFgrad(3),sumgrad2
if (allocated(cubmatvec)) then
	if (size(cubmatvec,2)==nx.and.size(cubmatvec,3)==ny.and.size(cubmatvec,4)==nz) return !Already generated
	deallocate(cubmatvec)
end if
allocate(cubmatvec(3,nx,ny,nz))
ifinish=0
write(*,*) "Generating grid data of gradient of electron density..."
!$OMP PARALLEL DO SHARED(cubmatvec,ifinish) PRIVATE(ix,iy,iz,tmpx,tmpy,tmpz,wfnval,wfnderv,gradrho,imo) schedule(dynamic) NUM_THREADS(nthreads)
do iz=1,nz
	do iy=1,ny
		do ix=1,nx
            call getgridxyz(ix,iy,iz,tmpx,tmpy,tmpz)
			call orbderv(2,1,nmo,tmpx,tmpy,tmpz,wfnval,wfnderv)
			gradrho=0D0
			do imo=1,nmo
				gradrho(:)=gradrho(:)+MOocc(imo)*wfnval(imo)*wfnderv(:,imo)
			end do
			cubmatvec(:,ix,iy,iz)=2*gradrho(:)
		end do
	end do
    !$OMP CRITICAL
    ifinish=ifinish+1
    call showprog(ifinish,nz)
    !$OMP END CRITICAL
end do
!$OMP END PARALLEL DO
end subroutine



!!------- Calculate high ELF localization domain populations and volumes (HELP and HELV), ChemPhysChem, 14, 3714 (2013)
subroutine HELP_HELV
use defvar
use functions
use basinintmod
use util
implicit real*8 (a-h,o-z)
real*8 :: ELFthres=0.5D0,rhothres=0.001D0
real*8 HELP,HELV,HELP_priv,HELV_priv
character c2000tmp*2000
integer tmparr(numrealatt)

do while(.true.)
    write(*,*)
    write(*,*) "    -------------------- Calculation of HELP and HELV --------------------"
    write(*,*) "-1 Return"
    write(*,*) " 0 Select basins and calculate their HELP and HELV"
    write(*,"(a,f6.3)") "  1 Define threshold of ELF, current: >",ELFthres
    write(*,"(a,f8.5,' a.u.')") "  2 Define threshold of electron density, current: >",rhothres
    read(*,*) isel
    if (isel==-1) then
        return
    else if (isel==1) then
        write(*,*) "Input threshold of ELF, e.g. 0.5"
        read(*,*) ELFthres
    else if (isel==2) then
        write(*,*) "Input threshold of electron density in a.u., e.g. 0.001"
        read(*,*) rhothres
    else if (isel==0) then
        write(*,*) "Input indices of the ELF basins for which HELP and HELV will be calculated"
        write(*,*) "e.g. 2,3,7-10,14"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nselbasin,tmparr)

        write(*,*) "Integrating, please wait patiently..."
        do idx=1,nselbasin
            ibasin=tmparr(idx)
            HELP=0
            HELV=0
            basinpop=0
            basinvol=0
            !$OMP PARALLEL private(ix,iy,iz,irealatt,rnowx,rnowy,rnowz,HELP_priv,HELV_priv,basinpop_priv,basinvol_priv,rhotmp,ELFtmp) &
            !$OMP shared(HELP,HELV,basinpop,basinvol) NUM_THREADS(nthreads)
            HELP_priv=0
            HELV_priv=0
            basinpop_priv=0
            basinvol_priv=0
            !$OMP do schedule(DYNAMIC)
			do iz=izlow,izup
				do iy=iylow,iyup
					do ix=ixlow,ixup
			            call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
			            irealatt=gridbas(ix,iy,iz)
                        rhotmp=fdens(rnowx,rnowy,rnowz)
                        ELFtmp=ELF_LOL(rnowx,rnowy,rnowz,"ELF")
                        if (irealatt==ibasin) then
                            if (rhotmp>rhothres.and.ELFtmp>ELFthres) then
                                HELP_priv=HELP_priv+rhotmp
                                HELV_priv=HELV_priv+1
                            end if
                            basinpop_priv=basinpop_priv+rhotmp
                            basinvol_priv=basinvol_priv+1
                        end if
		            end do
	            end do
            end do
            !$OMP end do
            !$OMP CRITICAL
                HELP=HELP+HELP_priv
                HELV=HELV+HELV_priv
                basinpop=basinpop+basinpop_priv
                basinvol=basinvol+basinvol_priv
            !$OMP end CRITICAL
            !$OMP END PARALLEL
            call calc_dvol(dvol)
            HELP=HELP*dvol
            HELV=HELV*dvol
            basinpop=basinpop*dvol
            basinvol=basinvol*dvol
            write(*,*)
            if (nselbasin>1) write(*,"(a,i6)") " ------ Information of basin",ibasin
            write(*,*) "Basin information: (constraints are not taken into account)"
            write(*,"(' Population:',f8.4,'   Volume:',f10.4,' Bohr^3')") basinpop,basinvol
            write(*,*)
            write(*,"(' High ELF localization domain population (HELP):',f10.4)") HELP
            write(*,"(' High ELF localization domain volume (HELV):    ',f10.4,' Bohr^3')") HELV
        end do
    end if
end do
end subroutine



!!------------ Determine correspondence between attractor index and atom index
!The passed in array is att2atm(numrealatt)
!Returned arry att2atm(i)=k means real attractor i corresponds to atom k. If k=0, that means attractor i should be NNA
subroutine calc_att2atm(att2atm)
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
integer att2atm(numrealatt)
att2atm=0
crit=0.3D0 !If distance between attractor and a nucleus is smaller than 0.3 Bohr, then the attractor will belong to the atom
do iatt=1,numrealatt !Cycle each attractors
	do iatm=1,ncenter
        disttest=dsqrt( (realattxyz(1,iatt)-a(iatm)%x)**2+(realattxyz(2,iatt)-a(iatm)%y)**2+(realattxyz(3,iatt)-a(iatm)%z)**2 )
        if (disttest<0.3D0) then 
	        att2atm(iatt)=iatm
            exit
        end if
    end do
end do
end subroutine



!!------------ Deallocate information of basin analysis
subroutine deallo_basinana(info)
use basinintmod
integer info
numatt=0
numrealatt=0
if (info==1.and.numatt>0) then
    write(*,*) "Note: Basin analysis information has been discarded"
end if
if (allocated(gridbas)) deallocate(gridbas)
if (allocated(attgrid)) deallocate(attgrid)
if (allocated(attval)) deallocate(attval)
if (allocated(attxyz)) deallocate(attxyz)
if (allocated(attconv)) deallocate(attconv)
if (allocated(nrealatthas)) deallocate(nrealatthas)
if (allocated(realatttable)) deallocate(realatttable)
if (allocated(realattval)) deallocate(realattval)
if (allocated(realattxyz)) deallocate(realattxyz)
if (allocated(interbasgrid)) deallocate(interbasgrid)
end subroutine






!!--------- Calculate basin overlap matrix using uniform or uniform + atomic center grid
!For closed-shell, the matrix will be stored to BOM in module basinintmod
!For open-shell, the alpha and beta matrix will be stored to BOM and BOMb in module basinintmod
!NNA and ECP are supported for both kinds of integration grid
subroutine genBOM
use defvar
use util
use basinintmod
use topo
use functions
implicit real*8 (a-h,o-z)
character c80tmp*80
real*8 orbval(nmo)
real*8,allocatable :: BOMsum(:,:),BOMsumb(:,:),BOMtmp(:,:),BOMtmpb(:,:)
!Below are used by uniform + atomic center grid integration
real*8 trustrad(numrealatt)
real*8 dens,grad(3),hess(3,3)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:)
type(content),allocatable :: gridatt(:) !Record correspondence between attractor and grid
integer att2atm(numrealatt) !The attractor corresponds to which atom. If =0, means this is a NNA
integer radpotAIM,sphpotAIM

if (ifuncbasin==1) then !Basin analysis for electron density
	write(*,*) "Use which kind of integration grid to evaluate basin overlap matrix?"
	write(*,*) "1 Uniform grid"
	write(*,*) "2 Mixed grid (uniform + atomic center grid)"
    write(*,*) "Note: 2 is much more accurate than 1 for core orbitals, but more expensive"
	read(*,*) igridkind
else
	igridkind=1
end if

call calc_dvol(dvol)
if (allocated(BOM)) deallocate(BOM)
if (allocated(BOMb)) deallocate(BOMb)

!Allocate space for BOM
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
	if (wfntype==3) then !R-post-HF, need to consider all orbitals
		nmatsize=nmo
	else !RHF,ROHF
		!High-lying virtual orbitals will be deleted, especially for .fch case (In fact, before entering this function, those >= LUMO+10 has already been discarded)
		!Notice that occupation number may be not contiguous, some low-lying orbitals may have
		!zero occupation due to modification by users, so we can't simply use nelec to determine matrix size
		do nmatsize=nmo,1,-1
			if (MOocc(nmatsize)/=0) exit
		end do
		if (nmo-nmatsize>0) write(*,"(' Note: The highest',i6,' virtual orbitals will not be taken into account')") nmo-nmatsize
	end if
	allocate(BOM(nmatsize,nmatsize,numrealatt),BOMtmp(nmatsize,nmatsize))
    nmatsizeb=0
	BOM=0
else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
	do iendalpha=nmo,1,-1
		if (MOtype(iendalpha)==1) exit
	end do
	if (wfntype==4) then !U-post-HF
		nmatsize=iendalpha !Total number of alpha orbitals
		nmatsizeb=nmo-nmatsize !Total number of beta orbitals
	else if (wfntype==1) then !UHF
		do nmatsize=iendalpha,1,-1
			if (MOocc(nmatsize)/=0D0) exit
		end do
		if (nint(nbelec)==0) then
			nmatsizeb=0
		else
			do nmatsizeb=nmo,iendalpha+1,-1
				if (MOocc(nmatsizeb)/=0D0) exit
			end do
			nmatsizeb=nmatsizeb-iendalpha
		end if
		if (iendalpha-nmatsize>0) write(*,"(' Note: The highest',i6,' alpha virtual orbitals will not be taken into account')") iendalpha-nmatsize
		if (nmo-iendalpha-nmatsizeb>0) write(*,"(' Note: The highest',i6,' beta virtual orbitals will not be taken into account')") nmo-iendalpha-nmatsizeb
	end if
	allocate(BOM(nmatsize,nmatsize,numrealatt),BOMb(nmatsizeb,nmatsizeb,numrealatt))
    allocate(BOMtmp(nmatsize,nmatsize),BOMtmpb(nmatsizeb,nmatsizeb))
	BOM=0
	BOMb=0
	MOinit=iendalpha+1
	MOend=iendalpha+nmatsizeb
end if

call walltime(iwalltime1)

!Start calculation of BOM
!!!!!!! Using uniform grid
if (igridkind==1) then
	do ibas=1,numrealatt
		write(*,"(' Generating orbital overlap matrix for basin',i6,'  of',i6,' ...')") ibas,numrealatt
		!$OMP parallel shared(BOM,BOMb) private(ix,iy,iz,rnowx,rnowy,rnowz,imo,jmo,imotmp,jmotmp,BOMtmp,BOMtmpb,orbval) num_threads(nthreads)
		BOMtmp=0D0
		if (nmatsizeb>0) BOMtmpb=0D0
		!$OMP do schedule(DYNAMIC)
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
					if (gridbas(ix,iy,iz)==ibas) then
						call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
                        !Total or alpha part
						call orbderv(1,1,nmatsize,rnowx,rnowy,rnowz,orbval)
						do imo=1,nmatsize
							do jmo=imo,nmatsize
								BOMtmp(imo,jmo)=BOMtmp(imo,jmo)+orbval(imo)*orbval(jmo)
							end do
						end do
                        !Beta part of UHF,U-post-HF
                        if (nmatsizeb>0) then
							call orbderv(1,MOinit,MOend,rnowx,rnowy,rnowz,orbval)
							do imo=MOinit,MOend
								imotmp=imo-iendalpha !So that the index start from 1 to nbelec
								do jmo=imo,MOend
									jmotmp=jmo-iendalpha
									BOMtmpb(imotmp,jmotmp)=BOMtmpb(imotmp,jmotmp)+orbval(imo)*orbval(jmo)
								end do
							end do
                        end if
					end if
				end do
			end do
		end do
		!$OMP end do
		!$OMP CRITICAL
		BOM(:,:,ibas)=BOM(:,:,ibas)+BOMtmp(:,:)*dvol
        if (nmatsizeb>0) BOMb(:,:,ibas)=BOMb(:,:,ibas)+BOMtmpb(:,:)*dvol
		!$OMP end CRITICAL
		!$OMP end parallel
	end do
    
!!!!!!! Using uniform + atomic center grid
else if (igridkind==2) then
	nbeckeiter=8
	numcp=0
	att2atm=0

	!Determine trust radius and then integrate in the trust sphere
	!For each electron attractor determined by basin analysis, we identify its corresponding atom, then use Newton method to obtain exact NCP position
	!The grid employed in this stage is centered from NCP
	write(*,*) "Integrating in trust sphere..."
	do iatt=1,numrealatt !Cycle each electron density attractor
    
		do iatm=1,ncenter !Try to find corresponding atom
			disttest=dsqrt( (realattxyz(1,iatt)-a(iatm)%x)**2+(realattxyz(2,iatt)-a(iatm)%y)**2+(realattxyz(3,iatt)-a(iatm)%z)**2 )
			if (disttest<0.3D0) then !If distance between attractor and a nucleus is smaller than 0.3 Bohr, then the attractor corresponds to the atom
				att2atm(iatt)=iatm
				write(*,"(' Attractor',i6,' corresponds to atom',i6,' (',a,')')") iatt,iatm,a(iatm)%name
				!Find exact NCP position from the nucleus and record to CPpos
				numcpold=numcp
				call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,1) !If successfully converge to a NCP, the position is write to CPpos(:,numcp)
				if (numcp==numcpold) then !Failed to converge to a NCP
					write(*,*) "Note: Unable to locate exact CP position! Use nuclear position"
					numcp=numcp+1
					CPpos(1,numcp)=a(iatm)%x
					CPpos(2,numcp)=a(iatm)%y
					CPpos(3,numcp)=a(iatm)%z
				end if
				exit
			end if
		end do
		if (att2atm(iatt)==0) then !No real atom corresponds to this attractor
			write(*,"(a,i6,a)") " Warning: Unable to determine the attractor",iatt," belongs to which atom!"
			write(*,"(a)") " If this is a non-nuclear attractor, simply press ENTER button to continue. If you used pseudopotential &
			&and this attractor corresponds to the cluster of all maxima of its valence electron, then input the index of this atom (e.g. 9). &
			&Else you should input q to return and regenerate basins with smaller grid spacing"
			read(*,"(a)") c80tmp
			if (c80tmp=='q') then
				return
			else if (c80tmp==" ") then
				numcpold=numcp
				call findcp(realattxyz(1,iatt),realattxyz(2,iatt),realattxyz(3,iatt),1)
				if (numcp==numcpold) then
					write(*,*) "Unable to locate exact CP position! Use attractor position instead"
					numcp=numcp+1
					CPpos(:,numcp)=realattxyz(:,iatt)
				end if
			else !ECP, input corresponding atom by user and directly use its nuclear position as grid center
				read(c80tmp,*) iatmtmp
				att2atm(iatt)=iatmtmp
				numcp=numcp+1
				CPpos(1,numcp)=a(iatmtmp)%x
				CPpos(2,numcp)=a(iatmtmp)%y
				CPpos(3,numcp)=a(iatmtmp)%z
			end if
		end if
	
		!Set integration points and weights, and meantime determine trust radius
		radpotAIM=200
		parm=1
		isettrustrad=0
		nintgrid=0 !Then number of integration grids within trust radius
		if (allocated(gridatt)) deallocate(gridatt) !Used to record grids in trust sphere of this attractor
		allocate(gridatt(radpotAIM*500))
		do ish=1,radpotAIM !Cycle each radial shell. Radius distance is from near to far
			!Becke, namely the second-kind Gauss-Chebyshev
			itmp=radpotAIM+1-ish !Invert ish to make radr from near to far
			radx=cos(itmp*pi/(radpotAIM+1D0))
			radr=(1+radx)/(1-radx)*parm
			radw=2*pi/(radpotAIM+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
		
			!Set Lebedev grids according to shell radius
			!For more inner shell, the distribution is more akin to spherically symmetric, therefore lower number of grids could be used
			if (att2atm(iatt)==0) then !NNA
				sphpotAIM=302
			else
				radtmp=covr(a(att2atm(iatt))%index)
				if (radr<0.2D0*radtmp) then
					sphpotAIM=26
				else if (radr<0.5D0*radtmp) then
					sphpotAIM=74
				else if (radr<0.8D0*radtmp) then
					sphpotAIM=146
				else
					sphpotAIM=194
				end if
			end if
			if (allocated(potx)) deallocate(potx,poty,potz,potw)
			allocate(potx(sphpotAIM),poty(sphpotAIM),potz(sphpotAIM),potw(sphpotAIM))
			call Lebedevgen(sphpotAIM,potx,poty,potz,potw)
			!Combine radial point and weights with angular part, and make them centered at current attractor
			gridatt( nintgrid+1:nintgrid+sphpotAIM )%x=radr*potx+CPpos(1,numcp)
			gridatt( nintgrid+1:nintgrid+sphpotAIM )%y=radr*poty+CPpos(2,numcp)
			gridatt( nintgrid+1:nintgrid+sphpotAIM )%z=radr*potz+CPpos(3,numcp)
			gridatt( nintgrid+1:nintgrid+sphpotAIM )%value=radw*potw
			!Find out trust radius for present attractor
			!If in a shell, the angle between "linking line between nucleus and a shell point" and "gradient vector of this point" &
			!is larger than 45 degree, then this shell is trust radius
			angmax=0
			if (att2atm(iatt)==0) then
				radinit=0
			else
				radrinit=0.15D0
				if (a(att2atm(iatt))%index>2) radrinit=0.5D0
			end if
			if (isettrustrad==0.and.radr>radrinit) then
				do isphpt=1,sphpotAIM
					xtmp=gridatt(nintgrid+isphpt)%x
					ytmp=gridatt(nintgrid+isphpt)%y
					ztmp=gridatt(nintgrid+isphpt)%z
					call calchessmat_dens(1,xtmp,ytmp,ztmp,dens,grad,hess) !Only value and gradient
					dirx=CPpos(1,numcp)-xtmp
					diry=CPpos(2,numcp)-ytmp
					dirz=CPpos(3,numcp)-ztmp
					angtmp=vecang(dirx,diry,dirz,grad(1),grad(2),grad(3))
					if (angtmp>angmax) angmax=angtmp
					if (angtmp>45) then
						write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr',/)") iatt,trustrad(iatt)
						isettrustrad=1 !The radius of last shell should be the final trust radius. Now exit
						exit
					end if
				end do
				if (isettrustrad==0) trustrad(iatt)=radr !Passed this shell and temporarily set the radius as trust radius. Continue to enlarge the trust radius, until reached angmax>45 degree
			end if
			nintgrid=nintgrid+sphpotAIM
		end do
		if (isettrustrad==0) then !Trust radius was not set after run over all shells
			trustrad(iatt)=1000 !Infinite, for isolated atom
			write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr',/)") iatt,trustrad(iatt)
		end if
	
		!Use DFT integration algorithm to integrate the region inside trust radius
		!$OMP PARALLEL private(ipt,ptx,pty,ptz,rx,ry,rz,dist,tmps,iter,switchwei,orbval,BOMtmp,BOMtmpb,imo,jmo,imotmp,jmotmp) shared(BOM,BOMb) NUM_THREADS(nthreads)
        BOMtmp=0
        if (nmatsizeb>0) BOMtmpb=0
		!$OMP do schedule(DYNAMIC)
		do ipt=1,nintgrid
			ptx=gridatt(ipt)%x
			pty=gridatt(ipt)%y
			ptz=gridatt(ipt)%z
			rx=ptx-CPpos(1,numcp) !The relative distance between current point to corresponding attractor
			ry=pty-CPpos(2,numcp)
			rz=ptz-CPpos(3,numcp)
			!Calculate switching function
			dist=dsqrt(rx*rx+ry*ry+rz*rz)
			tmps=dist-trustrad(iatt)
			if (tmps>1) then
				switchwei=0
			else if (tmps<-1) then
				switchwei=1
			else
				do iter=1,nbeckeiter
					tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
				end do
				switchwei=0.5D0*(1-tmps)
			end if
			if (switchwei<1D-7) cycle !For saving computational time
            
            !Total or alpha part
			call orbderv(1,1,nmatsize,ptx,pty,ptz,orbval)
			do imo=1,nmatsize
				do jmo=imo,nmatsize
					BOMtmp(imo,jmo)=BOMtmp(imo,jmo)+orbval(imo)*orbval(jmo)*gridatt(ipt)%value*switchwei
				end do
			end do
            !Beta part of UHF,U-post-HF
            if (nmatsizeb>0) then
				call orbderv(1,MOinit,MOend,ptx,pty,ptz,orbval)
				do imo=MOinit,MOend
					imotmp=imo-iendalpha !So that the index start from 1 to nbelec
					do jmo=imo,MOend
						jmotmp=jmo-iendalpha
						BOMtmpb(imotmp,jmotmp)=BOMtmpb(imotmp,jmotmp)+orbval(imo)*orbval(jmo)*gridatt(ipt)%value*switchwei
					end do
				end do
            end if
		end do
		!$OMP end do
		!$OMP CRITICAL
		BOM(:,:,iatt)=BOM(:,:,iatt)+BOMtmp(:,:)
        if (nmatsizeb>0) BOMb(:,:,iatt)=BOMb(:,:,iatt)+BOMtmpb(:,:)
		!$OMP end CRITICAL
		!$OMP END PARALLEL
	end do !End cycle attractors

    	!Integrating uniform grids
	write(*,*) "Integrating uniform grids..."
	ifinish=0
    
    call showprog(0,numrealatt)
	do ibas=1,numrealatt
		!$OMP PARALLEL private(ix,iy,iz,iatt,rnowx,rnowy,rnowz,rx,ry,rz,dist,tmps,iter,switchwei,orbval,tmpval,imo,jmo,imotmp,jmotmp,BOMtmp,BOMtmpb) shared(BOM,BOMb) NUM_THREADS(nthreads)
		BOMtmp=0D0
		if (nmatsizeb>0) BOMtmpb=0D0
		!$OMP do schedule(DYNAMIC)
		do iz=izlow,izup
			do iy=iylow,iyup
				do ix=ixlow,ixup
                    call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
					iatt=gridbas(ix,iy,iz)
                    if (iatt/=ibas) cycle
					rx=rnowx-CPpos(1,iatt) !The relative distance between current point to corresponding attractor
					ry=rnowy-CPpos(2,iatt)
					rz=rnowz-CPpos(3,iatt)
					!Calculate switching function at current grid
					dist=dsqrt(rx*rx+ry*ry+rz*rz)
					tmps=dist-trustrad(iatt)
					if (tmps>1) then
						switchwei=0
					else if (tmps<-1) then
						switchwei=1
					else
						do iter=1,nbeckeiter
							tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
						end do
						switchwei=0.5D0*(1-tmps)
					end if
					switchwei=1-switchwei
					if (switchwei<1D-7) cycle !For saving time
                
					!Total or alpha part
					call orbderv(1,1,nmatsize,rnowx,rnowy,rnowz,orbval)
					do imo=1,nmatsize
						tmpval=orbval(imo)*switchwei
						do jmo=imo,nmatsize
							BOMtmp(imo,jmo)=BOMtmp(imo,jmo)+orbval(jmo)*tmpval
						end do
					end do
					!Beta part of UHF,U-post-HF
					if (nmatsizeb>0) then
						call orbderv(1,MOinit,MOend,rnowx,rnowy,rnowz,orbval)
						do imo=MOinit,MOend
							imotmp=imo-iendalpha !So that the index start from 1 to nbelec
							tmpval=orbval(imo)*switchwei
							do jmo=imo,MOend
								jmotmp=jmo-iendalpha
								BOMtmpb(imotmp,jmotmp)=BOMtmpb(imotmp,jmotmp)+orbval(jmo)*tmpval
							end do
						end do
					end if
				end do
			end do
		end do
		!$OMP end do
		!$OMP CRITICAL
		BOM(:,:,ibas)=BOM(:,:,ibas)+BOMtmp(:,:)*dvol
		if (nmatsizeb>0) BOMb(:,:,ibas)=BOMb(:,:,ibas)+BOMtmpb(:,:)*dvol
		!$OMP end CRITICAL
		!$OMP end parallel
		call showprog(ibas,numrealatt)
    end do
end if

!Generate another half of matrix
do imo=1,nmatsize
	do jmo=imo,nmatsize
		BOM(jmo,imo,:)=BOM(imo,jmo,:)
	end do
end do
if (nmatsizeb>0) then
	do imo=1,nmatsizeb
		do jmo=imo,nmatsizeb
			BOMb(jmo,imo,:)=BOMb(imo,jmo,:)
		end do
	end do
end if

call walltime(iwalltime2)
write(*,"(' Generating BOM took up wall clock time',i10,' s')") iwalltime2-iwalltime1

!Check sanity of BOM
write(*,*)
allocate(BOMsum(nmatsize,nmatsize))
BOMsum=0
if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
	do ibas=1,numrealatt
		BOMsum=BOMsum+BOM(:,:,ibas)
	end do
	BOMerror=identmaterr(BOMsum)/numrealatt
	write(*,"(' Error of BOM is',f14.8)") BOMerror
	if (BOMerror>0.02D0) write(*,"(a)") " Warning: The integration is not very accurate"
! 	call showmatgau(BOMsum,"BOMsum",0,"f14.8",6)
else if (wfntype==1.or.wfntype==4) then !UHF,U-post-HF
	do ibas=1,numrealatt
		BOMsum=BOMsum+BOM(:,:,ibas)
	end do
	BOMerrora=identmaterr(BOMsum)/numrealatt
	write(*,"(' Error of alpha BOM is',f14.8)") BOMerrora
	if (nmatsizeb>0) then
		allocate(BOMsumb(nmatsizeb,nmatsizeb))
		BOMsumb=0
		do ibas=1,numrealatt
			BOMsumb=BOMsumb+BOMb(:,:,ibas)
		end do
		BOMerrorb=identmaterr(BOMsumb)/numrealatt
		write(*,"(' Error of Beta BOM is ',f14.8)") BOMerrorb
    else
		BOMerrorb=0D0
    end if
	if (BOMerrora>0.02D0.or.BOMerrorb>0.02D0) write(*,"(a)") " Warning: The integration is not very accurate"
end if
end subroutine





!!-------- Output basin overlap matrix (BOM) or atomic overlap matrix (AOM)
!itask=1: Output BOM to BOM.txt in current folder
!itask=2: Convert BOM to AOM and then output to AOM.txt in current folder
subroutine outBOMAOM(itask)
use defvar
use basinintmod
use util
integer itask,maplist(ncenter)

write(*,*) "Generating basin overlap matrix (BOM)..."
call genBOM !Generate BOM
nmatsizeb=size(BOMb,1)

if (itask==1) then !Output BOM
	open(10,file="BOM.txt",status="replace")
	if (wfntype==0.or.wfntype==2.or.wfntype==3) then
		do ibas=1,numrealatt
			write(10,"('Orbital overlap matrix of basin',i6)") ibas
			call showmatgau(BOM(:,:,ibas),"",1,"f14.8",10)
			write(10,*)
		end do
	else if (wfntype==1.or.wfntype==4) then
		do ibas=1,numrealatt
			write(10,"('Alpha part of orbital overlap matrix of basin',i6)") ibas
			call showmatgau(BOM(:,:,ibas),"",1,"f14.8",10)
			if (nmatsizeb>0) then
				write(10,"('Beta part of orbital overlap matrix of basin',i6)") ibas
				call showmatgau(BOMb(:,:,ibas),"",1,"f14.8",10)
			end if
			write(10,*)
		end do
	end if
	close(10)
	write(*,*)
	write(*,*) "Done, the matrices have been exported to BOM.txt in current folder"
else if (itask==2) then !Output AOM
    write(*,*)
    call atmidx2attidx(maplist,1,6)
	open(10,file="AOM.txt",status="replace")
    do iatm=1,ncenter
        ibas=maplist(iatm)
	    if (wfntype==0.or.wfntype==2.or.wfntype==3) then
		    write(10,"('Atomic overlap matrix of',i6,'(',a,')')") iatm,a(iatm)%name
		    call showmatgau(BOM(:,:,ibas),"",1,"f14.8",10)
		    write(10,*)
	    else if (wfntype==1.or.wfntype==4) then
		    write(10,"('Alpha part of atomic overlap matrix of',i6,'(',a,')')") iatm,a(iatm)%name
		    call showmatgau(BOM(:,:,ibas),"",1,"f14.8",10)
		    if (nmatsizeb>0) then
			    write(10,"('Beta part of atomic overlap matrix of',i6,'(',a,')')") iatm,a(iatm)%name
			    call showmatgau(BOMb(:,:,ibas),"",1,"f14.8",10)
		    end if
		    write(10,*)
        end if
    end do
	close(10)
	write(*,*)
	write(*,*) "Done, the matrices have been exported to AOM.txt in current folder"
end if
end subroutine




!!-------- Assign ELF basin labels
!This function can be used only when wavefunction information is available
subroutine assignELFbasinlab
use defvar
use basinintmod
use functions
implicit real*8 (a-h,o-z)
character(len=200) :: basinlab(numrealatt) !Label of each basin
character c80tmp*80,c200tmp*200,c10000tmp*10000
integer :: basinnatm(numrealatt) !=-1: Core basin, =n: n-synapatic basin
integer :: basinatm(20,numrealatt) !basinatm(1:abs(basinnatm(i)),i): Atom indices of basin i
integer :: tmparr(20),arrtmp(numrealatt)
real*8 :: basinpop(-1:numrealatt),basinvol(-1:numrealatt),basinpop_priv(-1:numrealatt),basinvol_priv(-1:numrealatt) !Population and volume of basins
integer :: basinidx(numrealatt) !Record indices of sorted basins according to basin label
!ELFcore_r(i) is the radius (Bohr) of outermost ELF minimum of isolate atom of element i, manually determined for atomic wavefunctions in http://sobereva.com/235
!A few transition metals have this radius even smaller than 1.0, may be due to ELF is inadequate to &
!show shell structure for certain cases, for these elements the average of neighbouring elements is taken
!Pd: 0.614, replaced with (2.876+2.762)/2=2.819
!Pt: 0.734, replaced with (2.442+2.584)/2=2.513
!Au: 0.714, replaced with (2.442+2.584)/2=2.513
real*8 :: ELFcore_r(0:nelesupp)=(/ 0D0,0D0,1D0,& !Ghost,H,He(1~2). Arbitrarily determined, since there is no core
1.512D0,1.012D0,0.726D0,0.57D0,0.468D0,0.392D0,0.338D0,0.296D0,& !Li~Ne(3~10)
2.102D0,1.674D0,1.386D0,1.174D0,1.018D0,0.904D0,0.81D0,0.732D0,& !Na~Ar(11~18)
3.028D0,2.498D0,2.382D0,2.290D0,2.200D0,2.356D0,2.052D0,1.988D0,1.940D0,& !K~Co(19~27)
1.884D0,2.134D0,1.786D0,1.524D0,1.376D0,1.252D0,1.158D0,1.076D0,1.002D0,& !Ni~Kr(28~36)
3.396D0,2.878D0,2.668D0,2.574D0,2.746D0,2.940D0,2.462D0,2.852D0,2.876D0,& !Rb~Rh(37~45)
2.819D0,2.762D0,2.298D0,1.974D0,1.766D0,1.618D0,1.512D0,1.416D0,1.330D0,& !Pd~Xe(46~54)
3.906D0,3.330D0,2.998D0,2.938D0,3.084D0,3.036D0,2.990D0,2.948D0,2.912D0,2.682D0,2.830D0,& !Cs~Tb(55~65)
2.794D0,2.760D0,2.730D0,2.700D0,2.672D0,2.412D0,2.368D0,2.338D0,2.344D0,2.364D0,2.370D0,2.442D0,& !Dy~Ir(66~77)
2.513D0,2.513D0,2.584D0,2.310D0,1.942D0,1.736D0,1.634D0,1.540D0,1.454D0,4.112D0,3.402D0,& !Pt~Ra(78~88)
3.182D0,2.990D0,3.036D0,3.088D0,3.036D0,3.136D0,3.086D0,2.750D0,2.994D0,2.956D0,2.916D0,2.874D0,2.836D0,2.802D0,2.448D0,(3D0,i=104,nelesupp) /) !Ac~Lr(89~103),~all

if (any(a%index>103)) then
	write(*,"(/,a)") " Error: There are atom(s) with element index larger than 103, this function is unable to use in this case"
	return
end if

basinlab=" "
basinnatm=0

!Assigning labels for core basins. If distance between an attractor and a nucleus is smaller than ELFcore_r, then do the assignment
write(*,*) "Assigning labels for core basins..."
do irealatt=1,numrealatt !Cycle each real attractor
    cycatt: do idx=1,nrealatthas(irealatt) !Cycle each sub attractor
		iatt=realatttable(irealatt,idx)
		attx=attxyz(1,iatt)
		atty=attxyz(2,iatt)
		attz=attxyz(3,iatt)
		do iatm=1,ncenter
			dist=dsqrt((attx-a(iatm)%x)**2+(atty-a(iatm)%y)**2+(attz-a(iatm)%z)**2)
			thres=ELFcore_r(a(iatm)%index)
			if (dist<thres) then
				basinnatm(irealatt)=-1
				basinatm(1,irealatt)=iatm
				exit cycatt
			end if
		end do
    end do cycatt
end do

!Assigning labels for valence basins. If any grid (boundary, rho>0.001) belonging to an attractor is neighbouring to a grid of a core basin, &
!then the owner atom of this core basin grid will be added to member list of this attractor
write(*,*) "Assigning labels for valence basins..."
ntest=2 !Test two-layers neighbouring grid. If =1, then V(S,F) of SF6 will be incorrectly assigned as V(S), since there is a layer of V(F) grids between core basin of V(S) and C(F)
do irealatt=1,numrealatt !Cycle all attractors
	icycle=0
	if (basinnatm(irealatt)==-1) then
		icycle=1
		cycle !Already assigned as core basin
    end if
    if (nrealatthas(irealatt)==1) then !Check if this attractor correspond to any H or He, if yes, add to member list
		do iatm=1,ncenter
			if (a(iatm)%index<=2) then
				attx=realattxyz(1,irealatt)
				atty=realattxyz(2,irealatt)
				attz=realattxyz(3,irealatt)
				dist=dsqrt((attx-a(iatm)%x)**2+(atty-a(iatm)%y)**2+(attz-a(iatm)%z)**2)
                if (dist<0.2D0) then
					basinnatm(irealatt)=1
                    basinatm(1,irealatt)=iatm
					exit
                end if
			end if
        end do
    end if
	do iz=izlow,izup
		do iy=iylow,iyup
			do ix=ixlow,ixup
				!This is a boundary grid and belongs to this attractor
				if (interbasgrid(ix,iy,iz).and.gridbas(ix,iy,iz)==irealatt) then
					call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
                    if (fdens(rnowx,rnowy,rnowz)<0.001D0) cycle !More safe, avoiding considering very distant grids, which are noisy
					!Cycle surrounding grids
					cycsurr: do k=iz-ntest,iz+ntest
						do j=iy-ntest,iy+ntest
							do i=ix-ntest,ix+ntest
								if (i<1.or.i>nx.or.j<1.or.j>ny.or.k<1.or.k>nz) cycle
								neiatt=gridbas(i,j,k)
								if (neiatt<=0) cycle cycsurr !This surrounding grid was not assigned to an atom
								if (basinnatm(neiatt)==-1) then !This grid belongs to a core basin
									jatm=basinatm(1,neiatt) !The current basin is neightbouring to jatm
                                    nmember=basinnatm(irealatt)
                                    if (all(basinatm(1:nmember,irealatt)/=jatm)) then !Add jatm to member list of irealatt if it hasn't occur
										basinnatm(irealatt)=nmember+1
                                        basinatm(nmember+1,irealatt)=jatm
                                    end if
									exit cycsurr
                                end if
                            end do
                        end do
                    end do cycsurr
                end if
            end do
        end do
	end do
    call showprog(irealatt,numrealatt)
end do
if (icycle==1) call showprog(numrealatt,numrealatt)

!Sort member indices from small to large for multiple-synapatic basins
do iatt=1,numrealatt
	nmember=basinnatm(iatt)
	if (nmember>=2) then
		do i=1,nmember
			do j=i+1,nmember
				if (basinatm(i,iatt)>basinatm(j,iatt)) then
					itmp=basinatm(i,iatt)
                    basinatm(i,iatt)=basinatm(j,iatt)
                    basinatm(j,iatt)=itmp
                end if
            end do
        end do
    end if
end do

!Generate label strings
do iatt=1,numrealatt
	if (basinnatm(iatt)==-1) then
		iatm=basinatm(1,iatt)
		write(c80tmp,"(i5)") iatm
		basinlab(iatt)="C("//trim(a(iatm)%name)//trim(adjustl(c80tmp))//")"
    else
		c200tmp=" "
		do idx=1,basinnatm(iatt)
			iatm=basinatm(idx,iatt)
			write(c80tmp,"(i5)") iatm
            if (idx==1) then
				c200tmp=trim(a(iatm)%name)//trim(adjustl(c80tmp))
            else
				c200tmp=trim(c200tmp)//','//trim(a(iatm)%name)//trim(adjustl(c80tmp))
            end if
        end do
		basinlab(iatt)="V("//trim(adjustl(c200tmp))//")"
    end if
end do

!Evaluating basin populations and volumes
write(*,*) "Evaluating basin populations and volumes, please wait..."
basinpop=0D0
basinvol=0D0
ifinish=0
!$OMP PARALLEL private(ix,iy,iz,irealatt,rnowx,rnowy,rnowz,basinpop_priv,basinvol_priv) shared(basinpop,basinvol,ifinish) NUM_THREADS(nthreads)
basinpop_priv=0D0
basinvol_priv=0D0
!$OMP do schedule(DYNAMIC)
do iz=izlow,izup
	do iy=iylow,iyup
		do ix=ixlow,ixup
			call getgridxyz(ix,iy,iz,rnowx,rnowy,rnowz)
			irealatt=gridbas(ix,iy,iz)
			basinpop_priv(irealatt)=basinpop_priv(irealatt)+fdens(rnowx,rnowy,rnowz)
			basinvol_priv(irealatt)=basinvol_priv(irealatt)+1
		end do
	end do
	!$OMP CRITICAL
    ifinish=ifinish+1
    call showprog(ifinish,izup-izlow+1)
	!$OMP end CRITICAL
end do
!$OMP end do
!$OMP CRITICAL
basinpop=basinpop+basinpop_priv
basinvol=basinvol+basinvol_priv
!$OMP end CRITICAL
!$OMP END PARALLEL
call calc_dvol(dvol)
basinpop=basinpop*dvol
basinvol=basinvol*dvol

!Output result according to basin indices
write(*,*)
write(*,*) "The following information is printed according to basin indices"
write(*,*) "Basin indices, populations (e), volumes (Angstrom^3) and assigned labels:"
sumcore=0
sumval=0
do iatt=1,numrealatt
	write(*,"(' Basin',i5,'  Pop.:',f8.4,'  Vol.:',f9.3,'  Label: ',a)") iatt,basinpop(iatt),basinvol(iatt)*b2a**3,trim(basinlab(iatt))
    if (basinnatm(iatt)==-1) then
		sumcore=sumcore+basinpop(iatt)
    else
		sumval=sumval+basinpop(iatt)
    end if
end do
write(*,*)
write(*,"(' Sum of core basin populations:   ',f12.4)") sumcore
write(*,"(' Sum of valence basin populations:',f12.4)") sumval
write(*,"(' Sum of all basin populations:    ',f12.4)") sumcore+sumval

!Sort basins, core basin first, then valence ones according to number of members. For each kind, sort according to basin index of 1,2,3... member
write(*,*)
write(*,*) "Sorting basins according to labels..."
forall(i=1:numrealatt) basinidx(i)=i
!First, sort according to number of members
do i=1,numrealatt
	do j=i+1,numrealatt
		if (basinnatm(basinidx(i))>basinnatm(basinidx(j))) then
            itmp=basinidx(i)
            basinidx(i)=basinidx(j)
            basinidx(j)=itmp
        end if
    end do
end do
!Second, sort each kind
do nmember=-1,20 !At most consider 20 member case
    if (count(basinnatm(:)==nmember)<=1) cycle !Do not need to sort
    do ibeg=1,numrealatt
		if (basinnatm(basinidx(ibeg))==nmember) exit
    end do
    do iend=numrealatt,1,-1
		if (basinnatm(basinidx(iend))==nmember) exit
    end do
    do imember=1,abs(nmember) !Cycle each position of basin label and sort
		do i=ibeg,iend !Cycle the range of basins which have nmember
			iatt=basinidx(i)
			do j=i+1,iend
				jatt=basinidx(j)
				if (basinatm(imember,iatt)>basinatm(imember,jatt)) then
					if (imember>=2) then !If changing order of imember will reverse order of imember-1, then do not change order
						if (basinatm(imember-1,iatt)<basinatm(imember-1,jatt)) cycle
                    end if
					itmp=basinidx(i)
					basinidx(i)=basinidx(j)
					basinidx(j)=itmp
                end if
            end do
        end do
    end do
end do

!Output result according to basin labels
write(*,*) "The following information is printed according to order of basin labels"
write(*,*) "Basin indices, populations (e), volumes (Angstrom^3) and assigned labels"
ncore=0
n1syn=0
n2syn=0
do idx=1,numrealatt
	iatt=basinidx(idx)
	write(*,"(' #',i5,'  Basin',i5,'  Pop.:',f8.4,'  Vol.:',f9.3,'  Label: ',a)") idx,iatt,basinpop(iatt),basinvol(iatt)*b2a**3,trim(basinlab(iatt))
	!write(*,"(' Number of members: ',i6)") basinnatm(iatt)
end do
write(*,*)
do nmember=-1,20
	if (nmember==0) cycle
	nbasin=count(basinnatm(:)==nmember)
	if (nbasin>0) then
		if (nmember==-1) then
			write(*,"(a,i6,a)") " Number of core basins is",nbasin,", their indices:"
        else
			write(*,"(a,i2,a,i6,a)") " Number of ",nmember,"-synaptic basins is",nbasin,", their indices:"
        end if
        itmp=0
        do idx=1,numrealatt
			iatt=basinidx(idx)
			if (basinnatm(iatt)==nmember) then
				itmp=itmp+1
				arrtmp(itmp)=iatt
            end if
        end do
		call arr2str_2(arrtmp(1:nbasin),c10000tmp)
        write(*,"(1x,a)") trim(c10000tmp)
	end if
end do
end subroutine




!!-------- Find outermost minimum of spherically averaged ELF for current system
!Mainly used to construct ELFcore_r in subroutine assignELFbasinlab
subroutine get_ELF_outer_minimum
use functions
implicit real*8 (a-h,o-z)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:),radpos(:)
nsphpt=170
allocate(potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt))
call Lebedevgen(nsphpt,potx,poty,potz,potw)
radr=5
sphavgval=0
idecrease=0
write(*,*) "Calculating..."
do while (radr>0.1D0)
    valold=sphavgval
	sphavgval=0
	do isph=1,nsphpt
		rnowx=potx(isph)*radr
		rnowy=poty(isph)*radr
		rnowz=potz(isph)*radr
		sphavgval=sphavgval+ELF_LOL(rnowx,rnowy,rnowz,"ELF")*potw(isph)
	end do
    !write(*,"(' r=',f10.4,' Bohr, avg. ELF=',f10.4)") radr,sphavgval
    if (sphavgval>valold.and.idecrease==1) then
		write(*,"(' Outermost core minimum is',f10.6,' Bohr')") radr
		exit
    end if
    if (sphavgval<valold) idecrease=1
    radr=radr-0.002D0
end do
end subroutine




!!--------- Generate core electron density grid data, store to corerhogrid(:,:,:)
subroutine getcorerhogrid
use defvar
use functions
use basinintmod
use util
implicit real*8 (a-h,o-z)
call readEDFlib(0)
write(*,*) "Generating grid data of core density..."
if (allocated(corerhogrid)) deallocate(corerhogrid)
allocate(corerhogrid(nx,ny,nz))
ifinish=0
!$OMP PARALLEL DO SHARED(corerhogrid,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			call EDFrho(1,tmpx,tmpy,tmpz,corerhogrid(i,j,k))
		end do
	end do
	!$OMP CRITICAL
    ifinish=ifinish+1
    call showprog(ifinish,nz)
    !$OMP END CRITICAL
end do
!$OMP END PARALLEL DO
call dealloEDF
!open(12,file="corerhogrid.cub",status="replace")
!call outcube(corerhogrid,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,12)
!close(12)
end subroutine




!!---------- Export attractor/basin information as basinana.txt and grid data as basinana.cub in current folder
subroutine export_basinana_info
use defvar
use basinintmod
implicit real*8 (a-h,o-z)

write(*,*) "Exporting present grid data to basinana.cub in current folder..."
open(10,file="basinana.cub",status="replace")
call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
close(10)
write(*,*) "Done!"

write(*,*) "Exporting attractors to basinana.txt in current folder..."
open(10,file="basinana.txt",status="replace")
write(10,"('ifuncbasin',i5)") ifuncbasin
write(10,"('numatt',i10)") numatt
write(10,"('numrealatt',i10)") numrealatt
write(10,*) "#Pristine attractors (ix,iy,iz,x,y,z,value)"
do iatt=1,numatt
	write(10,"(3i5,4E16.8)") attgrid(:,iatt),attxyz(:,iatt),attval(iatt)
end do
write(10,*) "#Real attractors (x,y,z,value)"
do irealatt=1,numrealatt
	write(10,"(4E16.8)") realattxyz(:,irealatt),realattval(irealatt)
end do
write(10,*) "#Number of pristine attractors that each real attractor has"
do irealatt=1,numrealatt
	write(10,"(12i6)") nrealatthas(irealatt)
end do
write(10,*) "#Conversion table from real attractors to pristine attractors"
do irealatt=1,numrealatt
	write(10,"(12i6)") realatttable(irealatt,1:nrealatthas(irealatt))
end do
write(10,*) "#Conversion table from pristine attractors to real attractors"
do iatt=1,numatt
	write(10,"(12i6)") attconv(iatt)
end do
write(10,*) "#Correspondence between every grid and basin index (gridbas)"
write(10,"(15i5)") (((gridbas(i,j,k),i=1,nx),j=1,ny),k=1,nz)
write(10,*) "#If grid is interbasin (interbasgrid)"
write(10,"(30L)") (((interbasgrid(i,j,k),i=1,nx),j=1,ny),k=1,nz)
close(10)
write(*,*) "Done!"
end subroutine




!!---------- Load attractor/basin information from basinana.txt and grid data from basinana.cub in current folder
subroutine load_basinana_info
use defvar
use basinintmod
implicit real*8 (a-h,o-z)
character c80tmp*80

if (allocated(cubmat)) deallocate(cubmat)
write(*,*) "Loading grid data from basinana.cub..."
call readcube("basinana.cub",1,0)
write(*,*) "Loading finished!"

call deallo_basinana(0)
write(*,*) "Loading attractors and basins information from basinana.txt..."
open(10,file="basinana.txt",status="old")
read(10,*) c80tmp,ifuncbasin
read(10,*) c80tmp,numatt
read(10,*) c80tmp,numrealatt
allocate(attgrid(3,numatt),attxyz(3,numatt),attval(numatt),attconv(numatt)) !Pristine attractors
allocate(nrealatthas(numrealatt),realattxyz(3,numrealatt),realattval(numrealatt),realatttable(numrealatt,numatt)) !Real attractors
allocate(gridbas(nx,ny,nz),interbasgrid(nx,ny,nz))

read(10,*) !Pristine attractors (ix,iy,iz,x,y,z,value)
do iatt=1,numatt
	read(10,"(3i5,4E16.8)") attgrid(:,iatt),attxyz(:,iatt),attval(iatt)
end do
read(10,*) !Real attractors (x,y,z,value)
do irealatt=1,numrealatt
	read(10,"(4E16.8)") realattxyz(:,irealatt),realattval(irealatt)
end do
read(10,*) !Number of pristine attractors that each real attractor has
read(10,*) nrealatthas(1:numrealatt)
read(10,*) !Conversion table from real attractors to pristine attractors
do irealatt=1,numrealatt
	read(10,*) realatttable(irealatt,1:nrealatthas(irealatt))
end do
read(10,*) !Conversion table from pristine attractors to real attractors"
read(10,*) attconv(1:numatt)
read(10,*) !Correspondence between every grid and basin index (gridbas)"
read(10,*) (((gridbas(i,j,k),i=1,nx),j=1,ny),k=1,nz)
read(10,*) !If grid is interbasin (interbasgrid)"
read(10,*) (((interbasgrid(i,j,k),i=1,nx),j=1,ny),k=1,nz)
close(10)
write(*,*) "Loading finished!"
end subroutine