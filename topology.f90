!============== Topology analysis for various kind of real space functions
subroutine topo_main
use topo
use defvar
use util
use GUI
use functions
use plot
use boysfunc
implicit real*8 (a-h,o-z)
character selectyn,c80tmp*80,c200*200,c1000*1000,c2000tmp*2000,ctmp1*20,ctmp2*20,ctmp3*20,ctmp4*20,icp1text*12,icp2text*12
real*8,allocatable :: randptx(:),randpty(:),randptz(:) !x,y,z of the points in the sphere
real*8,allocatable :: bassurpathtmp(:,:,:,:) !Used for temporary store bassurpath
real*8,allocatable :: shanCPrho(:) !For calculate Shannon aromaticity
integer shanCPind(100),searchcenlist(500)
integer :: isearch=0 !Searching mode of CP (suboption 10 in option -1)
integer :: nsearchlist=0 !The number of atoms considered in the searching modes 2,3,4,5 
integer,allocatable :: searchlist(:) !The list of the atoms considered in the searching modes 2,3,4,5
logical inlist1(ncenter),inlist2(ncenter) !If inlist1(i)=.true./.false., then atom i is in atom list 1 for interfragment searching
integer,allocatable :: tmparr(:),tmparr2(:)
real*8 hesstmp(3,3),gradtmp(3) !Temporarily used for calculating curvature
real*8 tvec(3)
integer :: ishowsearchpath=0

!Do not need virtual orbitals, remove them for faster calculation
if (iuserfunc/=27) call delvirorb(1) !Local electron affinity is related to virtual orbitals
call gen_GTFuniq(0) !Generate unique GTFs, for faster evaluation in orbderv

!Initialize searching and plotting parameters
toposphrad=3D0 !Radius of searching sphere is 3 Bohr
numsearchpt=1000 !1000 points randomly scattered in the sphere
sphcenx=0D0 !Position of the sphere center
sphceny=0D0
sphcenz=0D0
!Initialize plot parameter
ZVU=5D0 !Closer than other case
idrawisosur=0
ishow3n3=1
ishow3n1=1
ishow3p1=1
ishow3p3=1
idrawpath=1
bondradius=0.07D0 !For showing CPs, we use thiner default bond to avoid overlay the (3,-1)
ratioatmsphere=0.6D0 !Use smaller radius of atom
textheigh=36
ishowCPlab=0
ishowatmlab=0
ishowpathlab=0
ishowsearchlevel=0

!Backup parameters related to IRI and RDG, because they may be modified. When leaving this module, restore
tmp_IRI_rhocut=IRI_rhocut
tmp_RDG_maxrho=RDG_maxrho

call deallo_basinana(1) !If basin analysis has been performed, information should be dellocated otherwise attractors will be shown when visualizing CPs

write(*,*)
write(*,*) "         !!! Note: All length units in this module are Bohr !!!"
do while(.true.)
	write(*,*)
	write(*,"(a)") "              ================ Topology analysis ==============="
	if (numcp>0.or.numpath>0.or.nple3n1path>0.or.numbassurf>0) then
		write(*,"(a,i5)") " -11 Delete results and reselect real space function, current:",ifunctopo
    else
		write(*,"(a,i5)") " -11 Select real space function, current:",ifunctopo
    end if
	write(*,*) "-10 Return to main menu"
	write(*,*) "-9 Measure distances, angles and dihedral angles between CPs or atoms"
	write(*,*) "-5 Modify or print detail or export paths, or plot property along a path"
	write(*,*) "-4 Modify or export CPs (critical points)"
	write(*,*) "-3 Set interbasin surface generating parameters"
	write(*,*) "-2 Set path generating parameters"
	write(*,*) "-1 Set CP searching parameters"
    if (numcp<500.and.numpath<500) then
		write(*,*) "0 Print and visualize all generated CPs, paths and interbasin surfaces"
    else
		write(*,*) "0 Visualize all generated CPs, paths and interbasin surfaces"
		write(*,*) "00 Print information of all generated CPs and paths"
    end if
	write(*,*) "1 Search CPs from given starting points"
	write(*,*) "2 Search CPs from nuclear positions"
	write(*,*) "3 Search CPs from midpoint of atomic pairs"
	write(*,*) "4 Search CPs from triangle center of three atoms"
	write(*,*) "5 Search CPs from pyramid center of four atoms"
	write(*,*) "6 Search CPs from a batch of points within sphere(s)" !cube, random in sphere
	write(*,*) "7 Show real space function values at specific CP or all CPs"
	write(*,*) "8 Generating the paths connecting (3,-3) and (3,-1) CPs"
	write(*,*) "9 Generating the paths connecting (3,+1) and (3,+3) CPs"
	write(*,*) "10 Add or delete interbasin surfaces"
	if (ifunctopo==1) write(*,*) "20 Calculate Shannon aromaticity index"
	write(*,*) "21 Calculate gradient and curvature of electron density along a direction"
    
    read(*,"(a)") c80tmp
    if (c80tmp=="00") then !Special case
		call showtoposummary
        cycle
    else if (c80tmp=="vmin".or.c80tmp=="Vmin") then !Set most suitable setting for searching ESP minima via steepest descent
		ifunctopo=12
		itopomethod=4
		numsearchpt=10
		gradconv=1D-4
		dispconv=1D-5
		if (iESPcode==2.or.iESPcode==3) then
			call doinitlibreta(1)
			iaccurateESP=1
			if (isys==1.and.nthreads>12) nthreads=12
		end if
        isel=6
    else
		read(c80tmp,*) isel
    end if

	if (isel==-11) then
        write(*,"(/,a)") "  The following functions support both analytical gradient and Hessian:"
        write(*,"(a)") " Electron density, gradient norm of electron density, orbital wavefunction, IRI, RDG, vdW potential, Shannon entropy density & local information entropy, relative shannon entropy density, Fisher information density"
        write(*,"(a)") "  The following functions support analytic gradient and semi-numerical Hessian:"
        write(*,"(a)") " Laplacian of electron density, ELF, LOL, second Fisher information density"
        write(*,"(a)") "  All other functions do not have analytic derivative and thus their topology analysis may be relatively slow and numerical accuracy cannot be well guaranteed"
        write(*,*)
		write(*,*) "0 Return"
        call selfunc_interface(1,ifunctopo)
        if (ifunctopo==0) cycle
        if (ifdoESP(ifunctopo).and.(iESPcode==2.or.iESPcode==3)) then
            call doinitlibreta(1)
            iaccurateESP=1 !Ask libreta to calculate boys function in high precision, otherwise numerical gradient/Hessian will be too poor to converge
            if (isys==1.and.nthreads>12) nthreads=12
        end if
		if (numcp>0.or.numpath>0.or.nple3n1path>0.or.numbassurf>0) write(*,*) "Note: All found CPs, paths, surfaces have been cleaned"
		call deallo_topo
		!Set special parameters for specific real space functions
		if (ifunctopo==1.or.ifunctopo==2.or.ifunctopo==4.or.(ifunctopo==100.and.(iuserfunc==49.or.iuserfunc==50.or.iuserfunc==51))) then !Tight criteria for functions with analytical Hessian
			gradconv=1D-7
			dispconv=1D-8
		else if (ifunctopo==3.or.ifunctopo==9.or.ifunctopo==10.or.(ifunctopo==100.and.iuserfunc==52)) then !Looser criteria for functions with semi-numerical Hessian
			gradconv=1D-5
			dispconv=1D-6
		else !Use low criteria for functions without any order analytic derivative. IRI and RDG have analytical Hessian but also use this, because difficulty to converge to their CPs
			gradconv=1D-4
			dispconv=1D-5
		end if
        !write(*,"(' Note: Convergence thresholds of gradient and displacement have been set to',1PE12.4,' and',1PE12.4' a.u., respectively')") gradconv,dispconv
        write(*,"(a)") " Note: Some CP searching and path generating parameters have been automatically reset"
        itopomethod=1
        toposphrad=3D0
        nsurfpt=100
		if (ifunctopo==1) then !For electron density
			maxpathpttry=451
			pathstepsize=0.03D0
			nsurfpathpercp=60
			surfpathstpsiz=0.03D0
		else
			if (ifunctopo==4) toposphrad=1.5D0 !CPs for orbital wavefunction is very close to nuclei, so use smaller radii
			maxpathpttry=901
			pathstepsize=0.015D0 !Curvature is much larger than paths for rho, so smaller differentiate step must be used
			if (ifunctopo==4) numsearchpt=5000 !Use large value since evaluation of orbital wavefunction is very fast, but hard to locate CPs
			nsurfpathpercp=200
			surfpathstpsiz=0.008D0
		end if
        
        !For topology analysis of RDG, RDG_maxrho should be set to zero, otherwise steepest descent algorithm cannot be used
        if (ifunctopo==13.and.RDG_maxrho/=0) then
			write(*,"(a)") " Note: RDG_maxrho parameter has been temporarily set to 0"
			RDG_maxrho=0
        end if
        !For topology analysis of IRI/IRI-pi, IRI_rhocut should be set to zero to avoid automatically setting it to 5 in low rho region, which will lead to huge number of artificial extrema
        if (ifunctopo==24.and.IRI_rhocut/=0) then
			write(*,"(a)") " Note: IRI_rhocut parameter has been temporarily set to 0"
			IRI_rhocut=0
        end if
		if (ifunctopo==13.or.ifunctopo==24) then !RDG and IRI
			write(*,"(a)") " Note: CP searching method has been changed to steepest descent method, because it is most suitable for this case. &
            &In addition, gradient convergence criterion has been changed to a very large value (1000) to deactivate its effect, &
            &because for IRI or RDG, it is almost impossible to use steepest descent method to converge very accurately to a position with small enough gradient"
			itopomethod=4
            dispconv=0.0001D0
            gradconv=1000
            numsearchpt=1000
        else if (ifunctopo==25) then !vdW potential
			write(*,"(a)") " Note: CP searching method has been changed to steepest descent method, because it is most suitable for this case"
			itopomethod=4
            numsearchpt=100
		end if
        
	else if (isel==-10) then
		call del_GTFuniq !Destory unique GTF informtaion
        call delvirorb_back(1)
        if (IRI_rhocut/=tmp_IRI_rhocut) then
			write(*,*) "Note: Original IRI_rhocut parameter has been restored"
			IRI_rhocut=tmp_IRI_rhocut
        end if
        if (RDG_maxrho/=tmp_RDG_maxrho) then
			write(*,*) "Note: Original RDG_maxrho parameter has been restored"
			RDG_maxrho=tmp_RDG_maxrho
        end if
		exit
        
!-9 -9 -9 -9 -9 -9 -9
	else if (isel==-9) then
		write(*,*) "q = quit"
		write(*,*) "Selection method: a? = Atom?, c? = Critical point ?"
		write(*,*) "e.g. ""a1 c3"" returns the distance between atom1 and CP3"
		write(*,*) "     ""a4 a2"" returns the distance between atom4 and atom2"
		write(*,*) "     ""c6 a2 a5"" returns the angle of CP6-atom2-atom5"
		write(*,*) "     ""c2 c4 a3 c7"" returns the dihedral angle of CP2-CP4-atom3-CP7"
		do while(.true.)
			read(*,"(a)") c200
			c200=adjustl(c200)
			imeasure=0
			do ichar=1,len_trim(c200) !imeasure=1/2/3: measure distance,angle,dihedral
				if (c200(ichar:ichar)==','.or.c200(ichar:ichar)==' ') imeasure=imeasure+1
			end do
			if (c200(1:1)=='q') then
				exit
			else if (imeasure==1.or.imeasure==2.or.imeasure==3) then
				if (imeasure==1) read(c200,*) ctmp1,ctmp2 !Read two terms
				if (imeasure==2) read(c200,*) ctmp1,ctmp2,ctmp3 !Read three terms
				if (imeasure==3) read(c200,*) ctmp1,ctmp2,ctmp3,ctmp4 !Read four terms
				
				if (ctmp1(1:1)=='a') then
					read(ctmp1(2:),*) iatm
					tmpx1=a(iatm)%x
					tmpy1=a(iatm)%y
					tmpz1=a(iatm)%z
				else if (ctmp1(1:1)=='c') then
					read(ctmp1(2:),*) icp
					tmpx1=CPpos(1,icp)
					tmpy1=CPpos(2,icp)
					tmpz1=CPpos(3,icp)
				end if
				if (ctmp2(1:1)=='a') then
					read(ctmp2(2:),*) iatm
					tmpx2=a(iatm)%x
					tmpy2=a(iatm)%y
					tmpz2=a(iatm)%z
				else if (ctmp2(1:1)=='c') then
					read(ctmp2(2:),*) icp
					tmpx2=CPpos(1,icp)
					tmpy2=CPpos(2,icp)
					tmpz2=CPpos(3,icp)
				end if
				if (imeasure==1) then
                    tmpdist=xyz2dist(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2)
                    write(*,"(' The distance is',f12.6,' Bohr (',f12.6 ' Angstrom)')") tmpdist,tmpdist*b2a
                end if

				if (imeasure==2.or.imeasure==3) then !Analyze one more term, then print angle
					if (ctmp3(1:1)=='a') then
						read(ctmp3(2:),*) iatm
						tmpx3=a(iatm)%x
						tmpy3=a(iatm)%y
						tmpz3=a(iatm)%z
					else if (ctmp3(1:1)=='c') then
						read(ctmp3(2:),*) icp
						tmpx3=CPpos(1,icp)
						tmpy3=CPpos(2,icp)
						tmpz3=CPpos(3,icp)
					end if
				end if
				if (imeasure==2) write(*,"(' The angle is',f12.6,' degree')") xyz2angle(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3)
				
				if (imeasure==3) then !Analyze one more term, then print dihedral angle
					if (ctmp4(1:1)=='a') then
						read(ctmp4(2:),*) iatm
						tmpx4=a(iatm)%x
						tmpy4=a(iatm)%y
						tmpz4=a(iatm)%z
					else if (ctmp4(1:1)=='c') then
						read(ctmp4(2:),*) icp
						tmpx4=CPpos(1,icp)
						tmpy4=CPpos(2,icp)
						tmpz4=CPpos(3,icp)
					end if
					write(*,"(' The dihedral angle is',f12.6,' degree')") xyz2dih_sign(tmpx1,tmpy1,tmpz1,tmpx2,tmpy2,tmpz2,tmpx3,tmpy3,tmpz3,tmpx4,tmpy4,tmpz4)
				end if
			else
				write(*,*) "Input error"
			end if
		end do
		
!-5 -5 -5 -5 -5 -5 -5
!-5 -5 -5 -5 -5 -5 -5
!-5 -5 -5 -5 -5 -5 -5
	else if (isel==-5) then
		do while(.true.)
			write(*,*)
			write(*,*) "                      -------- Process paths --------"
			write(*,*) "0 Return"
			write(*,*) "1 Print summary of paths"
			write(*,*) "2 Print detail of a path"
			write(*,*) "3 Delete some paths"
			write(*,*) "4 Save points of all paths to paths.txt in current folder"
			write(*,*) "5 Load paths from an external file"
			write(*,*) "6 Export paths as paths.pdb file in current folder"
			if (ifPBC>0) write(*,"(a)") " -6 Export paths (including images at cell boundary if any) as paths.pdb file in current folder"
			write(*,*) "7 Calculate and plot specific real space function along a path"
			write(*,"(a)") " 8 Only retain bond paths (and corresponding CPs) connecting two specific molecular fragments while remove all other bond paths and BCPs"
			read(*,*) isel2
			
			if (isel2==0) then
				exit
				
			else if (isel2==1) then
				if (numpath>0) then
					do i=1,numpath
						call path_cp(i,icp1,icp2,ipathtype)
						if (icp1==0) then
							icp1text="   Unknown  "
						else
							write(icp1text,"(i5,1x,a)") icp1,CPtyp2lab(CPtype(icp1))
						end if
						if (icp2==0) then
							icp2text="  Unknown   "
						else
							write(icp2text,"(i5,1x,a)") icp2,CPtyp2lab(CPtype(icp2))
						end if
                        call getpathlength(i,pathlength)
						write(*,"(' #',i5,5x,'CP:',a,' --->',' CP:',a,'   Length:',f9.5)") i,icp1text,icp2text,pathlength
					end do
				else
					write(*,*) "No path has been generated"
				end if
				
			else if (isel2==2) then
				write(*,*) "Input the index of the path, e.g. 2"
				read(*,*) ipath
				if (ipath>numpath.or.ipath<=0) then
					write(*,*) "Error: Invalid path index!"
				else
                    call getpathlength(ipath,pathlength)
					write(*,"(a,i6,a,f10.5,a,i5)") " Path:",ipath,"   Length:",pathlength," Bohr   Total points:",pathnumpt(ipath)
					write(*,"(' From',3f18.12,/,' to  ',3f18.12)") topopath(:,1,ipath),topopath(:,pathnumpt(ipath),ipath)
					write(*,*) "The X/Y/Z coordinate (Bohr) and length of points in the path:"
                    curlen=0
					do ipt=1,pathnumpt(ipath)
                        if (ipt>1) then
                            if (ifPBC==0) then
                                disp=dsqrt(sum( (topopath(:,ipt,ipath)-topopath(:,ipt-1,ipath))**2) )
                            else
                                call nearest_dist(topopath(:,ipt,ipath),topopath(:,ipt-1,ipath),disp)
                            end if
                            curlen=curlen+disp 
                        end if
						write(*,"(i6,3f16.10,f10.5)") ipt,topopath(:,ipt,ipath),curlen
					end do
				end if
				
			else if (isel2==3) then
				write(*,*) "Input the index of the paths that will be deleted, e.g. 3,5,7-10,20"
				write(*,*) "Note: Input 0 can delete all paths"
				read(*,"(a)") c1000
				numpathold=numpath
				if (c1000(1:1)=='0') then
					numpath=0
				else
					call str2arr(c1000,ntmp)
					allocate(tmparr(ntmp))
					call str2arr(c1000,ntmp,tmparr)
					do itmp=1,ntmp
						idx=tmparr(itmp)
						if (idx<1.or.idx>numpath) cycle
						topopath(:,:,idx:numpath-1)=topopath(:,:,idx+1:numpath)
						pathnumpt(idx:numpath-1)=pathnumpt(idx+1:numpath)
						where(tmparr(itmp+1:ntmp)>idx) tmparr(itmp+1:ntmp)=tmparr(itmp+1:ntmp)-1
						numpath=numpath-1
					end do
					deallocate(tmparr)
				end if
				write(*,"(i6,' paths were deleted, now there are',i6,' paths left')") numpathold-numpath,numpath
				
			else if (isel2==4) then
				open(10,file="paths.txt",status="replace")
				write(10,"(2i10)") numpath
				do ipath=1,numpath
					write(10,"(/,'Path index:',i10)") ipath
					write(10,"(i10)") pathnumpt(ipath)
					do ipt=1,pathnumpt(ipath)
						write(10,"(3E20.12)") topopath(:,ipt,ipath)
					end do
				end do
				close(10)
				write(*,*) "Done, path information has been saved to paths.txt in current folder"
				write(*,"(a)") " Units are in Bohr. The first line corresponds to the number of paths"
				
			else if (isel2==5) then
				write(*,"(a)") " Note: The format of the input file must be identical to the one outputted by option 4"
				if (numpath>0) write(*,*) "Note: After loading the file, all current paths will be cleaned"
				write(*,*) "Input file path, e.g. C:\paths.txt"
                write(*,*) "If press ENTER button directly, paths.txt in current folder will be loaded"
				read(*,"(a)") c200
                if (c200==" ") c200="paths.txt"
				inquire(file=c200,exist=alive)
				if (alive.eqv..false.) then
					write(*,*) "Error: File was not found!"
				else
					open(10,file=c200,status="old")
					read(10,*) numpath
					do ipath=1,numpath
						read(10,*)
						read(10,*)
						read(10,*) pathnumpt(ipath)
						do ipt=1,pathnumpt(ipath) 
							read(10,*) topopath(:,ipt,ipath)
						end do
					end do
					close(10)
					write(*,*) "Done, path information has been recovered from the file"
				end if
				
			else if (isel2==6.or.isel2==-6) then
				open(10,file="paths.pdb",status="replace")
				if (ifPBC>0) then
					call getcellabc(asize,bsize,csize,alpha,beta,gamma)
					write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
				end if
				itmp=0
                if (isel2==6) then
					do ipath=1,numpath
						do ipt=1,pathnumpt(ipath)
							itmp=itmp+1
							write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",itmp,&
							' '//"C "//' ',"PTH",'A',ipath,topopath(:,ipt,ipath)*b2a,1.0,0.0,"C "
						end do
						write(10,"('TER')")
					end do
                else
					numpath_tmp=numpath
					pathnumpt_tmp=pathnumpt
					topopath_tmp=topopath
					forall(i=1:numpath) path_tmp_idx(i)=i
					call construct_pathtmp_withbound(numpath_tmp) !Generate paths at cell boundary
					do ipath=1,numpath_tmp
						do ipt=1,pathnumpt_tmp(ipath)
							itmp=itmp+1
							write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",itmp,&
							' '//"C "//' ',"PTH",'A',path_tmp_idx(ipath),topopath_tmp(:,ipt,ipath)*b2a,1.0,0.0,"C "
						end do
						write(10,"('TER')")
					end do
                end if
				close(10)
				write(*,"(a)") " Done, paths have been saved to paths.pdb in current folder. Each atom in this file &
                &corresponds to a point in the paths. Residue index corresponds path index"
				
			else if (isel2==7) then
				call plotpathprop
				
			else if (isel2==8) then
				call retainpath
			end if
		end do
		
!-4 -4 -4 -4 -4 -4 -4
!-4 -4 -4 -4 -4 -4 -4
!-4 -4 -4 -4 -4 -4 -4
	else if (isel==-4) then
		do while(.true.)
			write(*,*)
			write(*,*) "             ------------ Modify or export found CPs ------------"
			write(*,*) "0 Return"
			write(*,*) "1 Print summary of CPs (in Bohr)    -1 Print summary of CPs (in Angstrom)"
			write(*,*) "2 Delete some CPs"
			write(*,*) "3 Add a CP artificially"
			write(*,*) "4 Save CPs to CPs.txt in current folder"
			write(*,*) "5 Load CPs from a file"
			write(*,*) "6 Export CPs as CPs.pdb file in current folder"
			if (ifPBC>0) write(*,"(a)") " -6 Export CPs (including images at cell boundary if any) as CPs.pdb file in current folder"
			read(*,*) isel2
			
			if (isel2==0) then
				exit
			else if (isel2==1.or.isel2==-1) then
				if (numcp>0) then
					write(*,*) "Summary of found CPs:"
                    write(*,*) " Index         X               Y               Z         Type"
					do icp=1,numcp
						if (isel2==1) write(*,"(i6,3f16.9,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
						if (isel2==-1) write(*,"(i6,3f16.9,3x,a)") icp,CPpos(:,icp)*b2a,CPtyp2lab(CPtype(icp))
					end do
				else
					write(*,*) "No CP has been found"
				end if
			else if (isel2==2) then
				if (numcp==0) then
					write(*,*) "There is no CP currently"
					cycle
				end if
				if (numbassurf>0) then
					write(*,"(a)") " Warning: If one or more CPs are deleted, all interbasin surfaces will be lost, continue? 0/1=No/Yes"
					read(*,*) ifok
					if (ifok==0) then
						cycle
					else
						numbassurf=0
						cp2surf=0
						deallocate(bassurpath)
					end if
				end if
				do while(.true.)
					!The type of the CPs to be deleted will be temporarily set to -1, finally they will be ripped out
					write(*,*)
					write(*,*) "-1 Delete all CPs"
					write(*,*) "0 Return"
					write(*,*) "1 Delete CPs by inputting index range"
					write(*,*) "2 Delete CPs by specifying electron density range"
					write(*,*) "3 Delete all (3,-3) CPs"
					write(*,*) "4 Delete all (3,-1) CPs"
					write(*,*) "5 Delete all (3,+1) CPs"
					write(*,*) "6 Delete all (3,+3) CPs"
					write(*,*) "7 Delete CPs that far from a specific fragment"
					read(*,*) idel
					if (idel==-1) then
						numCP=0
						write(*,*) "Done, all CPs have been deleted"
					else if (idel==0) then
						exit
					else if (idel==1) then
						write(*,*) "Input index of the CPs you want to delete, e.g. 2,3,7-10"
						read(*,"(a)") c2000tmp
						call str2arr(c2000tmp,ntmp)
						allocate(tmparr(ntmp))
						call str2arr(c2000tmp,ntmp,tmparr)
						do itmp=1,ntmp
							CPtype(tmparr(itmp))=-1
						end do
						deallocate(tmparr)
					else if (idel==2) then
						write(*,"(a)") " Input lower and upper limits of electron density. All CPs within this range will be deleted, e.g. 0,0.1"
						read(*,*) rholow,rhohigh
						do icp=1,numcp
							tmpval=fdens(CPpos(1,icp),CPpos(2,icp),CPpos(3,icp))
							if (tmpval>=rholow.and.tmpval<=rhohigh) CPtype(icp)=-1
						end do
					else if (idel==3) then
						where(CPtype==1) CPtype=-1
					else if (idel==4) then
						where(CPtype==2) CPtype=-1
					else if (idel==5) then
						where(CPtype==3) CPtype=-1
					else if (idel==6) then
						where(CPtype==4) CPtype=-1
					else if (idel==7) then
						write(*,"(a)") " Input a factor, e.g. 1.3. If distances between a CP and all atoms in a specific fragment &
						&are larger than corresponding atomic covalent radius multiplied by this factor, then the CP will be deleted"
						read(*,*) facdist
						write(*,*) "Input index of the atoms in the fragment, e.g. 2,3,7-10"
						read(*,"(a)") c2000tmp
						call str2arr(c2000tmp,ntmp)
						allocate(tmparr(ntmp))
						call str2arr(c2000tmp,ntmp,tmparr)
						do icp=1,numcp
							do idx=1,ntmp
								iatm=tmparr(idx)
								dist=dsqrt((CPpos(1,icp)-a(iatm)%x)**2+(CPpos(2,icp)-a(iatm)%y)**2+(CPpos(3,icp)-a(iatm)%z)**2)
								if (dist < covr(a(iatm)%index)*facdist) exit
							end do
							if (idx>ntmp) CPtype(icp)=-1
						end do
					end if
					!Actually delete the CPs that to be deleted from the CP list
					if (idel>0) then
						numcpold=numcp
						icp=1
						do while(icp<=numcp)
							if (CPtype(icp)==-1) then
								CPpos(:,icp:numcp-1)=CPpos(:,icp+1:numcp)
								CPtype(icp:numcp-1)=CPtype(icp+1:numcp)
								numcp=numcp-1
							end if
							if (numcp<icp) exit
							if (CPtype(icp)==-1) cycle
							icp=icp+1
						end do
						write(*,"(i6,' CPs were deleted, now there are',i6,' CPs left')") numcpold-numcp,numcp
					end if
				end do
			else if (isel2==3) then
				numcp=numcp+1
				write(*,*) "Input X/Y/Z coordinate in Bohr, e.g. 1.2,0.0,3.44"
				read(*,*) CPpos(1,numcp),CPpos(2,numcp),CPpos(3,numcp)
				write(*,*) "Input CP type, 1=(3,-3) 2=(3,-1) 3=(3,+1) 4=(3,+3)"
				read(*,*) CPtype(numcp)
			else if (isel2==4) then
				open(10,file="CPs.txt",status="replace")
				write(10,"(i6)") numcp
				do icp=1,numcp
					write(10,"(i6,3f12.6,3x,i4)") icp,CPpos(:,icp),CPtype(icp)
				end do
				close(10)
				write(*,*) "Done, CP information has been saved to CPs.txt in current folder"
				write(*,*) "Note: The last column is CP type, 1=(3,-3) 2=(3,-1) 3=(3,+1) 4=(3,+3)"
			else if (isel2==5) then
				write(*,*) "Input file path, e.g. C:\ltwd\CPs.txt"
				write(*,"(a)") " (The format of the file must be identical to the one outputted by option 4)"
				if (numcp>0) write(*,*) "Note: After loading the file, all found CPs will be cleaned"
                write(*,"(a)") " Note: If pressing ENTER button directly, CPs.txt in current folder will be loaded"
				read(*,"(a)") c200
                if (c200==" ") c200="CPs.txt"
				inquire(file=c200,exist=alive)
				if (alive.eqv..false.) then
					write(*,*) "Error: File was not found!"
				else
					open(10,file=c200,status="old")
					read(10,*) numcp
					do icp=1,numcp
						read(10,*) nouse,CPpos(:,icp),CPtype(icp)
					end do
					close(10)
					write(*,*) "Done, CP information has been retrieved from the file"
				end if
			else if (isel2==6.or.isel2==-6) then
				open(10,file="CPs.pdb",status="replace")
				if (ifPBC>0) then
					call getcellabc(asize,bsize,csize,alpha,beta,gamma)
					write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
				end if
				write(10,"('REMARK   C=(3,-3) N=(3,-1) O=(3,+1) F=(3,+3)')")
                if (isel2==6) then
					do icp=1,numcp
						if (CPtype(icp)==1) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"C "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"C "
						if (CPtype(icp)==2) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"N "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"N "
						if (CPtype(icp)==3) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"O "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"O "
						if (CPtype(icp)==4) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"F "//' ',"CPS",'A',1,CPpos(:,icp)*b2a,1.0,0.0,"F "
					end do
                else
					numcp_tmp=numcp
					CPpos_tmp=CPpos
					CPtype_tmp=CPtype
					forall(i=1:numCP) CP_tmp_idx(i)=i
					call construct_CPtmp_withbound(numcp_tmp)
					do icp=1,numcp_tmp
						if (CPtype_tmp(icp)==1) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"C "//' ',"CPS",'A',CP_tmp_idx(icp),CPpos_tmp(:,icp)*b2a,1.0,0.0,"C "
						if (CPtype_tmp(icp)==2) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"N "//' ',"CPS",'A',CP_tmp_idx(icp),CPpos_tmp(:,icp)*b2a,1.0,0.0,"N "
						if (CPtype_tmp(icp)==3) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"O "//' ',"CPS",'A',CP_tmp_idx(icp),CPpos_tmp(:,icp)*b2a,1.0,0.0,"O "
						if (CPtype_tmp(icp)==4) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") "HETATM",icp,' '//"F "//' ',"CPS",'A',CP_tmp_idx(icp),CPpos_tmp(:,icp)*b2a,1.0,0.0,"F "
					end do
                end if
				write(*,*) "Done, CP information has been saved to CPs.pdb in current folder"
				write(*,*) "Note: Element C/N/O/F correspond to (3,-3)/(3,-1)/(3,+1)/(3,+3), respectively"
                if (isel2==-6) write(*,"(a)") " For image CPs at cell boundary, their residue index corresponds to index of unique CP"
				close(10)
			end if
		end do
!-3 -3 -3 -3 -3 -3 -3
	else if (isel==-3) then
		do while(.true.)
			write(*,*)
			write(*,*) "     ------------ Set interbasin surface generating parameters ------------"
			write(*,"(a)")      " 0 Return"
			write(*,"(a,i6)")   " 1 Number of paths in each interbasin surface",nsurfpathpercp
			write(*,"(a,i6)")   " 2 Number of points in each interbasin surface path, current:",nsurfpt
			write(*,"(a,f8.4)") " 3 Stepsize, current:",surfpathstpsiz
			read(*,*) isel2

			if (isel2==0) then
				exit
			else
				if (allocated(bassurpath)) then
					write(*,"(a)") " If the parameter is changed, all already generated interbasin surfaces will be cleaned, OK? (y/n)"
					read(*,*) selectyn
					if (selectyn=='y'.or.selectyn=='Y') then
						deallocate(bassurpath)
						numbassurf=0
						CP2surf=0
					end if
				end if
				if (.not.allocated(bassurpath)) then
					if (isel2==1) then
				    	write(*,*) "Input the value, e.g. 80"
                        read(*,*) nsurfpathpercp
					else if (isel2==2) then
				    	write(*,*) "Input the value, e.g. 100"
                        read(*,*) nsurfpt
					else if (isel2==3) then
				    	write(*,*) "Input the value, e.g. 0.02"
                        read(*,*) surfpathstpsiz
                    end if
				end if
			end if
		end do
!-2 -2 -2 -2 -2 -2 -2
	else if (isel==-2) then
		do while(.true.)
			write(*,*)
			write(*,*) "           ------------ Set path generation parameters ------------"
			write(*,"(a)")      " 0 Return"
			write(*,"(a,i4)")   " 1 Maximum number of points of a path, current:",maxpathpttry
			write(*,"(a,f8.4)") " 2 Stepsize, current:",pathstepsize
			write(*,"(a,f8.4)") " 3 Stop generation if distance to any CP is smaller than:",discritpathfin
! 			write(*,"(a,i4)")   " 4 Maximal number of bisection, current:",npathtry
			if (ishowsearchpath==0) write(*,*) "5 Toggle if printing details of path generating procedure, current: No"
			if (ishowsearchpath==1) write(*,*) "5 Toggle if printing details of path generating procedure, current: Yes"
			read(*,*) isel2

			if (isel2==0) then
				exit
			else if (isel2==1) then
				write(*,*) "Input a value, e.g. 400"
				write(*,"(' Note: The value must be smaller than',i5)") maxpathpt
				read(*,*) maxpathpttry
			else if (isel2==2) then
				selectyn='y'
				if (numpath>0) then
					write(*,*) "Warning: All generated paths will be cleaned, OK? (y/n)"
					read(*,*) selectyn
				end if
				if (selectyn=='y'.or.selectyn=='Y') then
					write(*,*) "Input a value in Bohr, e.g. 0.03"
					read(*,*) pathstepsize
					numpath=0
					write(*,*) "Done! Generated paths have been discarded"
				end if
			else if (isel2==3) then
				write(*,*) "Input a value in Bohr, e.g. 0.05"
				read(*,*) discritpathfin
			else if (isel2==4) then
				write(*,*) "Input an integer, must >=2"
				read(*,*) npathtry
			else if (isel2==5) then
				if (ishowsearchpath==0) then
					ishowsearchpath=1
					if (nthreads>1) write(*,"(a)") " NOTE: The printed path generation details will be messed up since more than one thread is utilized"
				else if (ishowsearchpath==1) then
					ishowsearchpath=0
				end if
			end if
		end do
!-1 -1 -1 -1 -1 -1 -1
	else if (isel==-1) then
		do while(.true.)
			write(*,*)
			write(*,*)"              ------------ Set CP searching parameters ------------"
			write(*,"(a)") " -1 Restore to default settings"
			write(*,"(a)") " 0 Return"
			write(*,"(a,i5)") " 1 Set maximal iterations:",topomaxcyc
			write(*,"(a,f12.6)") " 2 Set scale factor of stepsize:",CPstepscale
			write(*,"(a,1PE12.5)") " 3 Criteria for gradient norm convergence:",gradconv
			write(*,"(a,1PE12.5)") " 4 Criteria for displacement convergence:",dispconv
			write(*,"(a,f12.6)") " 5 Minimal distance between CPs:",minicpdis
			write(*,"(a,f8.2)") " 6 Skip search if distance between atoms is longer than the sum of their vdW radius multiplied by:",vdwsumcrit
			if (ishowsearchlevel==0) write(*,"(a)") " 7 Set printing level of details of CP searching, current: None"
			if (ishowsearchlevel==1) write(*,"(a)") " 7 Set printing level of details of CP searching: Minor detail"
			if (ishowsearchlevel==2) write(*,"(a)") " 7 Set printing level of details of CP searching: Some detail"
			if (ishowsearchlevel==3) write(*,"(a)") " 7 Set printing level of details of CP searching: All detail"
			write(*,"(a,1PE15.8)") " 8 Set criterion for determining if Hessian matrix is singular:",singularcrit
			if (CPsearchlow==CPsearchhigh) then
				write(*,*) "9 Set value range for reserving CPs, current: Reserve All CPs"
			else
				write(*,"(a,1PE12.4,a,1PE12.4)") " 9 Set value range for reserving CPs, current:",CPsearchlow,' to',CPsearchhigh
			end if
			if (isearch==0) then
				write(*,*) "10 Set the atoms to be considered in searching modes 2, 3, 4, 5: Undefined"
			else if (isearch==1) then
				write(*,"(a,i5,a)") " 10 Set the atoms to be considered in searching modes 2, 3, 4, 5, current: Within",nsearchlist," atoms"
			else if (isearch==2) then
				write(*,"(a,i5,a,i5,a)") " 10 Set the atoms to be considered in searching modes 2, 3, 4, 5, current: Betweem",count(inlist1.eqv..true.)," atoms and",count(inlist2.eqv..true.)," atoms"
			end if
            if (topotrustrad==0) then
                write(*,*) "11 Set trust radius of searching, current: Undefined"
            else
                write(*,"(a,f12.6)") " 11 Set trust radius of searching, current:",topotrustrad
            end if
            if (itopomethod==1) write(*,*) "12 Choose searching algorithm, current: Newton"
            if (itopomethod==2) write(*,*) "12 Choose searching algorithm, current: Barzilai-Borwein"
            if (itopomethod==3) write(*,*) "12 Choose searching algorithm, current: Steepest ascent"
            if (itopomethod==4) write(*,*) "12 Choose searching algorithm, current: Steepest descent"
			read(*,*) isel2

			if (isel2==-1) then
				topomaxcyc=120
				CPstepscale=1D0
				gradonv=1D-7
				dispconv=1D-9
				minicpdis=0.03D0
				vdwsumcrit=1.2D0
				ishowsearchlevel=0
				CPsearchlow=0D0
				CPsearchhigh=0D0
				singularcrit=5D-22
                topotrustrad=0.5D0
                itopomethod=1
				nsearchlist=0
				if (allocated(searchlist)) deallocate(searchlist)
			else if (isel2==0) then
				exit
			else if (isel2==1) then
				write(*,*) "Input an integer, e.g. 120"
				read(*,*) topomaxcyc
			else if (isel2==2) then
				write(*,*) "Input a value, default is 1.0"
				read(*,*) CPstepscale
			else if (isel2==3) then
				write(*,*) "Input a value, e.g. 1E-7"
				read(*,*) gradconv
			else if (isel2==4) then
				write(*,*) "Input a value, e.g. 1E-9"
				read(*,*) dispconv
			else if (isel2==5) then
				write(*,*) "Input a value, e.g. 0.03"
				read(*,*) minicpdis
			else if (isel2==6) then
				write(*,*) "Input a value, e.g. 1.2"
				read(*,*) vdwsumcrit
			else if (isel2==7) then
				write(*,*) "0 Do not print details"
				write(*,*) "1 Print minor details"
				write(*,*) "2 Print some details"
				write(*,*) "3 Print all details"
				read(*,*) ishowsearchlevel
				if (nthreads>1) write(*,*) "Warning: The printed details may be messed up since parallel running is enabled"
			else if (isel2==8) then
				write(*,"(a)") " Input a value, if absolute value of determinant of Hessiant matrix is lower than this value, &
                &then the Hessian will be regarded as singular, and the CP search will stop. e.g. 1E-15"
				read(*,*) singularcrit
			else if (isel2==9) then
				write(*,"(a)") " Input lower and upper limits. For example, if you input 0.05,0.22, &
				&then during CP searching, only when the real space function of a new CP is between 0.05 and 0.22 then it will be reserved"
				write(*,*) "Note: If the two values are identical, all CPs will be reserved (default case)"
				read(*,*) CPsearchlow,CPsearchhigh
			else if (isel2==10) then
				write(*,*) "0 Do not set"
                write(*,*) "1 Input an atom list, the searching will be confined in this region"
                write(*,"(a)") " 2 Input two atom lists without any overlap, the searching will be carried out between the two regions"
                read(*,*) isearch
                if (isearch==0) then
					nsearchlist=0
                    if (allocated(searchlist)) deallocate(searchlist)
                else if (isearch==1) then
					if (allocated(searchlist)) deallocate(searchlist)
					write(*,*) "Input index of the atoms, e.g. 2,3,7-10"
					read(*,"(a)") c2000tmp
					call str2arr(c2000tmp,nsearchlist)
					allocate(searchlist(nsearchlist))
					call str2arr(c2000tmp,nsearchlist,searchlist)
					if (any(searchlist>ncenter).or.any(searchlist<0)) then
						write(*,*) "Error: The range exceeded valid range!"
                        isearch=0
						nsearchlist=0
						deallocate(searchlist)
					else
						write(*,"(' Done! There are',i6,' atoms in the range')") nsearchlist
					end if
                else if (isearch==2) then
					write(*,*) "Input index of the atoms in atom list 1, e.g. 2,3,7-10"
					read(*,"(a)") c2000tmp
					call str2arr(c2000tmp,natmlist1)
					allocate(tmparr(natmlist1))
					call str2arr(c2000tmp,natmlist1,tmparr)
					if (any(tmparr>ncenter).or.any(tmparr<0)) then
						write(*,*) "Error: The range exceeded valid range!"
						deallocate(tmparr)
                        isearch=0
					else
						write(*,"(' Done! There are',i6,' atoms in the range')") natmlist1
                        inlist1=.false.
                        do itmp=1,natmlist1
							inlist1(tmparr(itmp))=.true.
                        end do
					end if
                    write(*,*)
					write(*,*) "Input index of the atoms in atom list 2, e.g. 1,4-6,11,12"
					read(*,"(a)") c2000tmp
					call str2arr(c2000tmp,natmlist2)
					allocate(tmparr2(natmlist2))
					call str2arr(c2000tmp,natmlist2,tmparr2)
					if (any(tmparr2>ncenter).or.any(tmparr2<0)) then
						write(*,*) "Error: The range exceeded valid range!"
						deallocate(tmparr2)
                        isearch=0
					else
						write(*,"(' Done! There are',i6,' atoms in the range')") natmlist2
                        inlist2=.false.
                        do itmp=1,natmlist2
							inlist2(tmparr2(itmp))=.true.
                        end do
					end if
                    nsearchlist=natmlist1+natmlist2
					if (allocated(searchlist)) deallocate(searchlist)
                    allocate(searchlist(nsearchlist))
                    searchlist(1:natmlist1)=tmparr(:)
                    searchlist(natmlist1+1:)=tmparr2(:)
                    deallocate(tmparr,tmparr2)
                end if
            else if (isel2==11) then
                write(*,*) "Input the trust radius in Bohr, e.g. 0.2"
                write(*,"(a)") " Note: If norm of displacement is larger than this value, the displacement vector will be scaled to this value. &
                &To disable this consideration, input 0"
                read(*,*) topotrustrad
            else if (isel2==12) then
                write(*,"(a)") " Choose the CP searching algorithm"
                write(*,*) "1 Newton"
                write(*,*) "2 Barzilai-Borwein"
                write(*,*) "3 Steepest ascent"
                write(*,*) "4 Steepest descent"
                write(*,"(a)") " Note: 1 and 2 can be used to search all kinds of CPs, while 3 and 4 only locate maxima and minima, respectively"
                read(*,*) itopomethod
			end if
		end do
		
	!0000000000000000000 Visualize and inspect resulting topology information
	else if (isel==0) then
		if (numcp<500.and.numpath<500) then
			call showtoposummary
        else
			call showPHrelat
			write(*,"(/,a)") " Note: Since there are too many critical points and/or paths, their summary is not automatically shown. &
            &To print them, input 00 in the topology analysis menu"
        end if
		if (isilent==0) call drawmoltopogui
        
	!111111111111111111111 Locate CPs from one or more given starting points
	else if (isel==1) then
        do while(.true.)
            write(*,*)
            write(*,*) "0 Return"
            write(*,*) "1 Manually input a starting point"
            write(*,*) "2 Using an atom or midpoint of two atoms as starting point"
            write(*,*) "3 Using all atom coordinates from a given .pdb or .pqr file as starting points"
            write(*,*) "4 Using starting points from a .txt file"
            read(*,*) isel2
            if (isel2==0) exit
		    numcpold=numcp
            if (isel2==1) then
		        write(*,*) "Input X,Y,Z of starting point"
                write(*,*) "e.g. ""2.0,3.1,-0.5"" means (2.0,3.1,-0.5) Bohr"
                write(*,*) "     ""2.5,3.2,-0.1 A"" means (2.5,3.2,-0.1) Angstrom"
                read(*,"(a)") c200
                read(c200,*) x,y,z
                if (index(c200,'A')/=0.or.index(c200,'a')/=0) then
                    x=x/b2a;y=y/b2a;z=z/b2a
                end if
		        call findcp(x,y,z,ifunctopo)
            else if (isel2==2) then
		        write(*,"(a)") " Input index of atom to use its nuclear position as starting point, e.g. 5"
                write(*,*) "  or indices of two atoms to use their midpoint as starting point, e.g. 3,6"
                read(*,"(a)") c200
                if (index(c200,',')/=0) then
    			    read(c200,*) iatm,jatm
			        x=(a(iatm)%x+a(jatm)%x)/2
			        y=(a(iatm)%y+a(jatm)%y)/2
			        z=(a(iatm)%z+a(jatm)%z)/2
                else
                    read(c200,*) iatm
			        x=a(iatm)%x
			        y=a(iatm)%y
			        z=a(iatm)%z
                end if
		        call findcp(x,y,z,ifunctopo)
            else if (isel2==3.or.isel2==4) then
                if (isel2==3) then
                    write(*,*) "Input path of the .pdb or .pqr file, e.g. C:\Cutie_panther.pdb"
                else if (isel2==4) then
                    write(*,*) "Input path of the .txt file, e.g. C:\maki.txt"
                    write(*,*) "Each line of this file should contain X,Y,Z of a starting point in Bohr"
                    write(*,"(a)") " If pressing ENTER button directly, attractors.txt in current folder will be loaded if it exists"
                end if
                do while(.true.)
	                read(*,"(a)") c200
                    if (c200==" ") then
						c200="attractors.txt"
                        exit
                    else
						inquire(file=c200,exist=alive)
						if (alive) exit
						write(*,*) "Cannot find the file, input again!"
                    end if
                end do
                open(10,file=c200,status="old")
                ncalc=0
                do while(.true.)
	                read(10,"(a)",iostat=ierror) c200
	                if (ierror/=0) exit
                    if (isel2==3.and.(c200(1:6)=="HETATM".or.c200(1:6)=="ATOM  ")) then
                        ncalc=ncalc+1
		                read(c200,"(30x,3f8.3)") x,y,z
                        write(*,"(' Doing starting point:',i6,'   X,Y,Z:',3f8.3,' Angstrom')") ncalc,x,y,z
		                call findcp(x/b2a,y/b2a,z/b2a,ifunctopo)
                    else if (isel2==4.and.c200/=" ") then
                        ncalc=ncalc+1
		                read(c200,*) x,y,z
                        write(*,"(' Doing starting point:',i6,'   X,Y,Z:',3f12.6,' Bohr')") ncalc,x,y,z
		                call findcp(x,y,z,ifunctopo)
                    end if
                end do
                close(10)
            end if
		    if (numcp==numcpold) then
			    write(*,*) "No new critical point was found"
		    else
			    write(*,"(' Found',i5,' new critical points')") numcp-numcpold
				if ((numcp-numcpold)/=0) then
					call sortCP(numcpold+1)
					write(*,*) "                            ---- Summary ----"
					write(*,*) " Index                       Coordinate               Type"
					do icp=numcpold+1,numcp
						write(*,"(i6,3f15.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
					end do
				end if
		    end if
        end do
		
	!Typical searching modes (mainly for AIM)
	else if (isel==2.or.isel==3.or.isel==4.or.isel==5) then
		if (isearch==0) then !Construct a temporary list, including all atoms
			nsearchlist=ncenter
			allocate(searchlist(nsearchlist))
			forall(i=1:ncenter) searchlist(i)=i
        else
			continue !searchlist has already been defined by user
		end if
		
		!222222222222222222222, search CPs from nuclear positions
		if (isel==2) then
			numcpold=numcp
			inow=0
            gradconv_bk=gradconv
            gradconv=1  !Remove gradient convergence criterion, because for very heavy atoms, the nucleus has too sharp density to converge
			!$OMP PARALLEL DO SHARED(itime) PRIVATE(iatm) schedule(dynamic) NUM_THREADS(nthreads)
			do idx=1,nsearchlist
				iatm=searchlist(idx)
				if (ishowsearchlevel>=1) then
					write(*,"(' #',i5,' /',i5,a,i5,'(',a,')')") iatm,ncenter,": Trying from nuclear position of ",iatm,a(iatm)%name
				else
					!$OMP CRITICAL
					inow=inow+1
					call showprog(inow,nsearchlist)
					!$OMP END CRITICAL
				end if
				call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,ifunctopo)
			end do
			!$OMP END PARALLEL DO
			if ((numcp-numcpold)/=0) then
				call sortCP(numcpold+1)
				write(*,*) "                            ---- Summary ----"
				write(*,*) " Index                       Coordinate               Type"
				do icp=numcpold+1,numcp
					write(*,"(i6,3f15.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
				end do
			end if
			write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
            gradconv=gradconv_bk
			if (ifunctopo==1.and.count(CPtype(1:numcp)==1)<nsearchlist) write(*,*) "Warning: Some (3,-3) may missing, try to search again with different parameters"
		
		!333333333333333333333, midpoint between two atoms as initial guess
		else if (isel==3) then
			if (ifPBC>0) write(*,*) "Note: Looping atoms in this mode considers image ones"
			numcpold=numcp
			itime=0
			ntime=0
            if (ifPBC==0) then
				do idx=1,nsearchlist !Test how many iterations will be done
					iatm=searchlist(idx)
					do jdx=idx+1,nsearchlist
						jatm=searchlist(jdx)
						if (isearch==2) then !Search will be performed between two fragments
							if ( (inlist1(iatm).and.inlist1(jatm)).or.(inlist2(iatm).and.inlist2(jatm)) ) cycle !Two atoms are in the same fragment
						end if
						disttmp=atomdist(iatm,jatm,1)
						if ( disttmp <= vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) ntime=ntime+1
					end do
				end do
				write(*,"(' Number of starting points:',i8)") ntime
				!$OMP PARALLEL DO SHARED(itime) PRIVATE(idx,iatm,jdx,jatm,disttmp) schedule(dynamic) NUM_THREADS(nthreads)
				do idx=1,nsearchlist
					iatm=searchlist(idx)
					do jdx=idx+1,nsearchlist
						jatm=searchlist(jdx)
						if (isearch==2) then !Search will be performed between two fragments
							if ( (inlist1(iatm).and.inlist1(jatm)).or.(inlist2(iatm).and.inlist2(jatm)) ) cycle !Two atoms are in the same fragment
						end if
						disttmp=atomdist(iatm,jatm,1)
						if ( disttmp > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
						!$OMP CRITICAL
						itime=itime+1
						if (ishowsearchlevel>=1) then
							write(*,"(' #',i5,' /',i5,a,i5,'(',a,')',a,i5,'(',a,')')") &
							itime,ntime,": Trying from midpoint between ",iatm,a(iatm)%name," and",jatm,a(jatm)%name
						else
							call showprog(itime,ntime)
						end if
						!$OMP end CRITICAL
						call findcp( (a(iatm)%x+a(jatm)%x)/2D0,(a(iatm)%y+a(jatm)%y)/2D0,(a(iatm)%z+a(jatm)%z)/2D0, ifunctopo)
					end do
				end do
				!$OMP END PARALLEL DO
            else !PBC case
				do idx=1,nsearchlist !Test how many iterations will be done
					iatm=searchlist(idx)
					do jdx=idx,nsearchlist
						jatm=searchlist(jdx)
						if (isearch==2) then !Search will be performed between two fragments
							if ( (inlist1(iatm).and.inlist1(jatm)).or.(inlist2(iatm).and.inlist2(jatm)) ) cycle !Two atoms are in the same fragment
						end if
						ntime=ntime+1
					end do
				end do
				write(*,"(' Number of atom pairs:',i8)") ntime
				!$OMP PARALLEL DO SHARED(itime) PRIVATE(idx,iatm,jdx,jatm,icell,jcell,kcell,disttmp,atmx,atmy,atmz) schedule(dynamic) NUM_THREADS(nthreads)
				do idx=1,nsearchlist
					iatm=searchlist(idx)
					do jdx=idx,nsearchlist
						jatm=searchlist(jdx)
						if (isearch==2) then !Search will be performed between two fragments
							if ( (inlist1(iatm).and.inlist1(jatm)).or.(inlist2(iatm).and.inlist2(jatm)) ) cycle !Two atoms are in the same fragment
						end if
						!$OMP CRITICAL
						itime=itime+1
						if (ishowsearchlevel>=1) then
							write(*,"(' #',i5,' /',i5,a,i5,'(',a,')',a,i5,'(',a,')')") &
							itime,ntime,": Trying between ",iatm,a(iatm)%name," and",jatm,a(jatm)%name
						else
							call showprog(itime,ntime)
						end if
						!$OMP end CRITICAL
                        do icell=-1,1
							do jcell=-1,1
								do kcell=-1,1
									if (jdx==idx.and.icell==0.and.jcell==0.and.kcell==0) cycle
									call tvec_PBC(icell,jcell,kcell,tvec)
                                    atmx=a(jatm)%x+tvec(1)
                                    atmy=a(jatm)%y+tvec(2)
                                    atmz=a(jatm)%z+tvec(3)
									disttmp=dsqrt( (atmx-a(iatm)%x)**2 + (atmy-a(iatm)%y)**2 + (atmz-a(iatm)%z)**2 )
									if ( disttmp > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
									call findcp( (a(iatm)%x+atmx)/2D0,(a(iatm)%y+atmy)/2D0,(a(iatm)%z+atmz)/2D0, ifunctopo)
                                end do
                            end do
                        end do
					end do
				end do
				!$OMP END PARALLEL DO
            end if
			if ((numcp-numcpold)/=0) then
				call sortCP(numcpold+1)
				write(*,*) "                            ---- Summary ----"
				write(*,*) " Index                       Coordinate               Type"
				do icp=numcpold+1,numcp
					write(*,"(i6,3f15.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
				end do
			end if
			write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
		
		!4444444444444444444, search CPs from triangle center of three atoms
        !If distance between any atomic pair is longer than threshold, the search will be skipped. For PBC case, only consider the other two atoms closest to it
		else if (isel==4) then
			if (ifPBC>0) write(*,*) "Note: Looping atoms in this mode does not consider image ones"
			numcpold=numcp
			itime=0
			ntime=0
			do idx=1,nsearchlist !Test how many iterations will be done
				iatm=searchlist(idx)
				do jdx=idx+1,nsearchlist
					jatm=searchlist(jdx)
                    call nearest_atmdistxyz(iatm,jatm,distij,atmjx,atmjy,atmjz)
					if ( distij > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
					do kdx=jdx+1,nsearchlist
						katm=searchlist(kdx)
						if (isearch==2) then !Search will be performed between two fragments
							!Three atoms are in the same fragment, skip
							if ( (inlist1(iatm).and.inlist1(jatm).and.inlist1(katm)).or.(inlist2(iatm).and.inlist2(jatm).and.inlist2(katm)) ) cycle
						end if
                        call nearest_atmdistxyz(iatm,katm,distki,atmkx,atmky,atmkz)
                        distkj=xyz2dist(atmkx,atmky,atmkz,atmjx,atmjy,atmjz)
						if ( distki > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(iatm)%index)).or.&
						distkj > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
						ntime=ntime+1
					end do
				end do
			end do
            write(*,"(' Number of starting points:',i8)") ntime
			!$OMP PARALLEL DO SHARED(itime) PRIVATE(iatm,jatm,katm,distij,distki,distkj,atmjx,atmjy,atmjz,atmkx,atmky,atmkz) schedule(dynamic) NUM_THREADS(nthreads)
			do idx=1,nsearchlist
				iatm=searchlist(idx)
				do jdx=idx+1,nsearchlist
					jatm=searchlist(jdx)
                    call nearest_atmdistxyz(iatm,jatm,distij,atmjx,atmjy,atmjz)
					if ( distij > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
					do kdx=jdx+1,nsearchlist
						katm=searchlist(kdx)
						if (isearch==2) then !Search will be performed between two fragments
							!Three atoms are in the same fragment, skip
							if ( (inlist1(iatm).and.inlist1(jatm).and.inlist1(katm)).or.(inlist2(iatm).and.inlist2(jatm).and.inlist2(katm)) ) cycle
						end if
                        call nearest_atmdistxyz(iatm,katm,distki,atmkx,atmky,atmkz)
                        distkj=xyz2dist(atmkx,atmky,atmkz,atmjx,atmjy,atmjz)
						if ( distki > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(iatm)%index)).or.&
						distkj > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
						!$OMP CRITICAL
						itime=itime+1
						if (ishowsearchlevel>=1) then
							write(*,"(' #',i5,' /',i5,a ,i5,'(',a,')' ,i5,'(',a,')' ,i5,'(',a,')' )") &
							itime,ntime,": Trying from triangle center of ",iatm,a(iatm)%name,jatm,a(jatm)%name,katm,a(katm)%name
						else
							call showprog(itime,ntime)
						end if
						!$OMP end CRITICAL
						call findcp( (a(iatm)%x+atmjx+atmkx)/3D0,(a(iatm)%y+atmjy+atmky)/3D0,(a(iatm)%z+atmjz+atmkz)/3D0, ifunctopo)
					end do
				end do
			end do
			!$OMP END PARALLEL DO
			if ((numcp-numcpold)/=0) then
				call sortCP(numcpold+1)
				write(*,*) "                            ---- Summary ----"
				write(*,*) " Index                       Coordinate               Type"
				do icp=numcpold+1,numcp
					write(*,"(i6,3f15.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
				end do
			end if
			write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
	
		!5555555555555555555, search CPs from pyramid center of four atoms
		else if (isel==5) then
			if (ifPBC>0) write(*,*) "Note: Looping atoms in this mode does not consider image ones"
			numcpold=numcp
			itime=0 
			ntime=0
			do idx=1,nsearchlist !Test how many iterations will be done; ij,jk,kl,li,lj,ik
				iatm=searchlist(idx)
				do jdx=idx+1,nsearchlist
					jatm=searchlist(jdx)
                    call nearest_atmdistxyz(iatm,jatm,distij,atmjx,atmjy,atmjz)
					if ( distij > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
					do kdx=jdx+1,nsearchlist
						katm=searchlist(kdx)
                        call nearest_atmdistxyz(iatm,katm,distik,atmkx,atmky,atmkz)
                        distkj=xyz2dist(atmjx,atmjy,atmjz,atmkx,atmky,atmkz)
						if ( distkj > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
						do ldx=kdx+1,nsearchlist
							latm=searchlist(ldx)
							if (isearch==2) then !Search will be performed between two fragments
								!Four atoms are in the same fragment, skip
								if ( (inlist1(iatm).and.inlist1(jatm).and.inlist1(katm).and.inlist1(latm)).or.(inlist2(iatm).and.inlist2(jatm).and.inlist2(katm).and.inlist2(latm)) ) cycle
							end if
                            call nearest_atmdistxyz(iatm,latm,distli,atmlx,atmly,atmlz)
                            distlk=xyz2dist(atmlx,atmly,atmlz,atmkx,atmky,atmkz)
                            distlj=xyz2dist(atmlx,atmly,atmlz,atmjx,atmjy,atmjz)
							if ( distlk > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(katm)%index)).or.&
							distli > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(iatm)%index)).or.&
							distlj > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(jatm)%index)).or.&
							distik > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(katm)%index))) cycle
							ntime=ntime+1
						end do
					end do
				end do
			end do
            write(*,"(' Number of starting points:',i8)") ntime
			!$OMP PARALLEL DO SHARED(itime) PRIVATE(iatm,jatm,katm,latm,distij,distik,distkj,distli,distlk,distlj, &
            !$OMP atmjx,atmjy,atmjz,atmkx,atmky,atmkz,atmlx,atmly,atmlz) schedule(dynamic) NUM_THREADS(nthreads)
			do idx=1,nsearchlist
				iatm=searchlist(idx)
				do jdx=idx+1,nsearchlist
					jatm=searchlist(jdx)
                    call nearest_atmdistxyz(iatm,jatm,distij,atmjx,atmjy,atmjz)
					if ( distij > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(jatm)%index)) ) cycle
					do kdx=jdx+1,nsearchlist
						katm=searchlist(kdx)
                        call nearest_atmdistxyz(iatm,katm,distik,atmkx,atmky,atmkz)
                        distkj=xyz2dist(atmjx,atmjy,atmjz,atmkx,atmky,atmkz)
						if ( distkj > vdwsumcrit*(vdwr(a(katm)%index)+vdwr(a(jatm)%index)) ) cycle
						do ldx=kdx+1,nsearchlist
							latm=searchlist(ldx)
							if (isearch==2) then !Search will be performed between two fragments
								!Four atoms are in the same fragment, skip
								if ( (inlist1(iatm).and.inlist1(jatm).and.inlist1(katm).and.inlist1(latm)).or.(inlist2(iatm).and.inlist2(jatm).and.inlist2(katm).and.inlist2(latm)) ) cycle
							end if
                            call nearest_atmdistxyz(iatm,latm,distli,atmlx,atmly,atmlz)
                            distlk=xyz2dist(atmlx,atmly,atmlz,atmkx,atmky,atmkz)
                            distlj=xyz2dist(atmlx,atmly,atmlz,atmjx,atmjy,atmjz)
							if ( distlk > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(katm)%index)).or.&
							distli > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(iatm)%index)).or.&
							distlj > vdwsumcrit*(vdwr(a(latm)%index)+vdwr(a(jatm)%index)).or.&
							distik > vdwsumcrit*(vdwr(a(iatm)%index)+vdwr(a(katm)%index))) cycle
							!$OMP CRITICAL
							itime=itime+1
							if (ishowsearchlevel>=1) then
								write(*,"(' #',i5,' /',i5,a ,i5,'(',a,')' ,i5,'(',a,')' ,i5,'(',a,')' ,i5,'(',a,')')") &
								itime,ntime,": Trying from center of ",&
								iatm,a(iatm)%name,jatm,a(jatm)%name,katm,a(katm)%name,latm,a(latm)%name
							else
								call showprog(itime,ntime)
							end if
							!$OMP end CRITICAL
							call findcp( (a(iatm)%x+atmjx+atmkx+atmlx)/4D0,&
							(a(iatm)%y+atmjy+atmky+atmly)/4D0,(a(iatm)%z+atmjz+atmkz+atmlz)/4D0, ifunctopo)
						end do
					end do
				end do
			end do
			!$OMP END PARALLEL DO
			if ((numcp-numcpold)/=0) then
				call sortCP(numcpold+1)
				write(*,*) "                            ---- Summary ----"
				write(*,*) " Index                       Coordinate               Type"
				do icp=numcpold+1,numcp
					write(*,"(i6,3f15.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
				end do
			end if
			write(*,"(' Totally found',i6,' new critical points')") numcp-numcpold
		end if
		
		if (isearch==0) then !Discard the temporarily constructed list for CP searching mode 0
			nsearchlist=0
			deallocate(searchlist)
		end if
		
	!6666666666666666666. Force search using a large number of starting points
	else if (isel==6) then
		do while(.true.)
			write(*,*) "   --------------- Distribute starting points in sphere(s) ---------------"
			write(*,"(' Sphere center:',3f10.5,'   Radius:',f6.2,'   Points:',i8)") sphcenx,sphceny,sphcenz,toposphrad,numsearchpt
			write(*,"(a)") " -9 Return"
			write(*,"(a)") " -2 Start the search using some nuclei as sphere center in turn"
			write(*,"(a)") " -1 Start the search using each nucleus as sphere center in turn"
			write(*,"(a)") " 0 Start the search using the defined sphere center"
			write(*,"(a)") " 1 Input coordinate of the sphere center"
			write(*,"(a)") " 2 Set the sphere center at a nuclear position"
			write(*,"(a)") " 3 Set the sphere center at midpoint between two atoms"
			write(*,"(a)") " 4 Set the sphere center at triangle center of three atoms"
			write(*,"(a)") " 5 Set the sphere center at a CP"
			write(*,"(a)") " 6 Set the sphere center at midpoint between two CPs"			
			write(*,"(a)") " 10 Set the sphere radius"
			write(*,"(a)") " 11 Set the number of points in the sphere"
			read(*,*) isel2

			if (isel2==-9) then
				exit
			else if (isel2==0.or.isel2==-1.or.isel2==-2) then
				if (isel2==-2) then
					write(*,*) "Input the indices of the atoms, e.g. 3,4,5,12, at most 1000 characters"
					read(*,"(a)") c1000
					call str2arr(c1000,nsearchcen,searchcenlist)
				end if
				if (isel2==0) nsearchcen=1
				if (isel2==-1) nsearchcen=ncenter
				!4.189=4/3*pi, this approximately assesses upper limit of points need to be generated in the cube region, &
                !so that numsearchpt points will occur in the sphere
				numcpold=numcp
				numsearchpt_tmp=nint(8D0/4.189D0*numsearchpt)
				allocate(randptx(numsearchpt_tmp),randpty(numsearchpt_tmp),randptz(numsearchpt_tmp))
                sphcenx_old=sphcenx
                sphceny_old=sphceny
                sphcenz_old=sphcenz
				itime=0
				ioutcount=0
				call walltime(iwalltime1)
                call showprog(0,numsearchpt*nsearchcen)
				do icenidx=1,nsearchcen
					icen=icenidx
					if (isel2==-2) icen=searchcenlist(icenidx)
					if (isel2==-1.or.isel2==-2) then !Cycle each atom center
						sphcenx=a(icen)%x
						sphceny=a(icen)%y
						sphcenz=a(icen)%z
					end if
					CALL RANDOM_NUMBER(randptx)
					CALL RANDOM_NUMBER(randpty)
					CALL RANDOM_NUMBER(randptz)
					randptx=randptx*2*toposphrad+(sphcenx-toposphrad) !Move distribution center of random point to sphere center
					randpty=randpty*2*toposphrad+(sphceny-toposphrad)
					randptz=randptz*2*toposphrad+(sphcenz-toposphrad)
					randptx(1)=sphcenx !The first try point is set to sphere center, this faciliates to locate CP at nuclei
					randpty(1)=sphceny
					randptz(1)=sphcenz
					!$OMP PARALLEL DO SHARED(itime) PRIVATE(i) schedule(dynamic) NUM_THREADS(nthreads)
					do i=1,numsearchpt_tmp
						dispt_cen=dsqrt( (randptx(i)-sphcenx)**2+(randpty(i)-sphceny)**2+(randptz(i)-sphcenz)**2 )
						if (dispt_cen>toposphrad) cycle !Outside the expected sphere, ignore this starting point
						call findcp(randptx(i),randpty(i),randptz(i),ifunctopo)
						!$OMP CRITICAL
						itime=itime+1
                        if (itime<numsearchpt*nsearchcen) call showprog(itime,numsearchpt*nsearchcen)
						!$OMP end CRITICAL
					end do
					!$OMP END PARALLEL DO
				end do
				deallocate(randptx,randpty,randptz)
                call showprog(numsearchpt*nsearchcen,numsearchpt*nsearchcen)
				call walltime(iwalltime2)
				write(*,"(' Searching CPs took up',i8,' seconds wall clock time')") iwalltime2-iwalltime1
                sphcenx=sphcenx_old
                sphceny=sphceny_old
                sphcenz=sphcenz_old
				
				if ((numcp-numcpold)/=0) then
! 					call sortCP(numcpold+1) !Nonsense here, because the starting points occur randomly
					write(*,*)
					write(*,*) "                            ---- Summary ----"
					write(*,*) " Index                       Coordinate               Type"
					do icp=numcpold+1,numcp
						write(*,"(i6,3f15.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(CPtype(icp))
					end do
				end if
				write(*,"(' Totally find',i6,' new critical points')") numcp-numcpold
				write(*,*)
			else if (isel2==1) then
				write(*,*) "Input x,y,z of sphere center in Bohr, e.g. 1.2,0.2,-0.44"
				read(*,*) sphcenx,sphceny,sphcenz
			else if (isel2==2) then
				write(*,*) "Input atom index, e.g. 3"
				read(*,*) iatm
				if (iatm>ncenter.or.iatm<=0) then
					write(*,*) "Error: Invalid input"
				else
					sphcenx=a(iatm)%x
					sphceny=a(iatm)%y
					sphcenz=a(iatm)%z
				end if
			else if (isel2==3) then
				write(*,*) "Input index of the two atoms, e.g. 3,7"
				read(*,*) iatm,jatm
				if ( iatm>ncenter.or.iatm<=0.or.jatm>ncenter.or.jatm<=0 ) then
					write(*,*) "Invalid input"
				else
					sphcenx=(a(iatm)%x+a(jatm)%x)/2D0
					sphceny=(a(iatm)%y+a(jatm)%y)/2D0
					sphcenz=(a(iatm)%z+a(jatm)%z)/2D0
				end if
			else if (isel2==4) then
				write(*,*) "Input index of the three atoms, e.g. 2,3,7"
				read(*,*) iatm,jatm,katm
				if ( iatm>ncenter.or.iatm<=0.or.jatm>ncenter.or.jatm<=0.or.katm>ncenter.or.katm<=0 ) then
					write(*,*) "Invalid input"
				else
					sphcenx=(a(iatm)%x+a(jatm)%x+a(katm)%x)/3D0
					sphceny=(a(iatm)%y+a(jatm)%y+a(katm)%y)/3D0
					sphcenz=(a(iatm)%z+a(jatm)%z+a(katm)%z)/3D0
				end if
			else if (isel2==5) then
				write(*,*) "Input the CP index, e.g. 3"
				read(*,*) icp
				if (icp>numcp.or.icp<=0) then
					write(*,*) "Invalid input"
				else
					sphcenx=CPpos(1,icp)
					sphceny=CPpos(2,icp)
					sphcenz=CPpos(3,icp)
				end if
			else if (isel2==6) then
				write(*,*) "Input index of the two CPs, e.g. 3,7"
				read(*,*) icp,jcp
				if ( icp>numcp.or.icp<=0.or.jcp>numcp.or.jcp<=0 ) then
					write(*,*) "Invalid input"
				else
					sphcenx=(CPpos(1,icp)+CPpos(1,jcp))/2D0
					sphceny=(CPpos(2,icp)+CPpos(2,jcp))/2D0
					sphcenz=(CPpos(3,icp)+CPpos(3,jcp))/2D0
				end if			
			else if (isel2==10) then
				write(*,*) "Input a radius in Bohr, e.g. 3.5"
				read(*,*) toposphrad
			else if (isel2==11) then
				write(*,*) "Input a number, e.g. 200"
				read(*,*) numsearchpt
			end if
		end do
		
!7777777777777777777
	else if (isel==7) then
		write(*,*) "Input the index of the CP that you are interested in, e.g. 3"
		write(*,"(a)") " Note 1: If input 0, then properties of all CPs will be outputted to CPprop.txt in current folder &
		&(and if you feel the output speed is slow, you can input -1 to avoid outputting ESP, which is the most expensive one)"
		write(*,"(a)") " Note 2: If input CP index with ""d"" suffix, e.g. 17d, then property of this CP can be decomposed into orbital contribution"
		read(*,*) c200
		if (index(c200,"d")/=0) then
			read(c200(1:index(c200,"d")-1),*) indcp
			call decompptprop(CPpos(1,indcp),CPpos(2,indcp),CPpos(3,indcp))
		else
			read(c200,*) indcp
			if (indcp==0.or.indcp==-1) then
				write(*,*) "Please wait..."
				iback=ishowptESP
				if (indcp==-1) ishowptESP=0 !Avoid outputting ESP
                c200="CPprop.txt"
                if (iaddprefix==1) call addprefix(c200)
				open(10,file=c200,status="replace")
				do icp=1,numcp
					write(*,"(' Outputting properties of CP',i6,'  /',i6)") icp,numcp
					write(10,"(' ----------------   CP',i6,',     Type ',a,'   ----------------')") icp,CPtyp2lab(CPtype(icp))
                    if (ifunctopo==1) then
						if (CPtype(icp)==1) then !NCP of AIM, output its corresponding nucleus if possible
							call NCP2atmidx(icp,iatm)
							if (iatm==0) then
								write(10,"(' Corresponding nucleus:  Unknown')")
							else
								write(10,"(' Corresponding nucleus:',i6,'(',a,')')") iatm,a(iatm)%name
							end if
						else if (CPtype(icp)==2) then !BCP of AIM, output its corresponding nucleus if possible
							call BCP2atmidx(icp,iatm,jatm)
							if (iatm==0.or.jatm==0) then
								!write(10,"(' Connected atoms:  Unknown')") 
							else
								write(10,"(' Connected atoms:',i6,'(',a')','   --',i6,'(',a')')") iatm,a(iatm)%name,jatm,a(jatm)%name
							end if
                        end if
                    end if
					write(10,"(' Position (Bohr):    ',3f18.12)") CPpos(:,icp)
					write(10,"(' Position (Angstrom):',3f18.12)") CPpos(:,icp)*b2a
					call showptprop(CPpos(1,icp),CPpos(2,icp),CPpos(3,icp),ifunctopo,10)
					write(10,*)
				end do
				close(10)
				if (indcp==-1) ishowptESP=iback
				write(*,"(a)") " Done! The results have been outputted to "//trim(c200)//" in current folder"
				write(*,*) "Note: Unless otherwise specified, all units are in a.u."
			else if (indcp>0.and.indcp<=numcp) then
				write(*,*) "Note: Unless otherwise specified, all units are in a.u."
				write(*,"(' Position (Bohr):    ',3f18.12)") CPpos(:,indcp)
				write(*,"(' Position (Angstrom):',3f18.12)") CPpos(:,indcp)*b2a
				write(*,"(' CP type: ',a)") CPtyp2lab(CPtype(indcp))
				call showptprop(CPpos(1,indcp),CPpos(2,indcp),CPpos(3,indcp),ifunctopo,6)
			else
				write(*,*) "Error: Invalid input"
			end if
		end if
		
!8888888888888888888, generate the paths connecting (3,-3) and (3,-1)
	else if (isel==8) then
		numpathold=numpath
		num3n1=count(CPtype(1:numcp)==2)
		if (num3n1==0) then
			write(*,"(a)") " Error: You must at least locate one (3,-1) CP before trying to generate topology paths!"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		if (count(CPtype(1:numcp)==1)==0) then
			write(*,"(a)") " Error: You must locate (3,-3) CPs first before trying to generate topology paths!"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		write(*,"(' Number of (3,-1) CPs:',i6,'    Generating topology paths...')") num3n1
		inow=0
        call showprog(0,num3n1)
		!$OMP PARALLEL DO SHARED(inow) PRIVATE(i) schedule(dynamic) NUM_THREADS(nthreads)
		do i=1,numcp
			if (CPtype(i)==2) then
				if (ishowsearchpath==0) then
					call findpath(i,1,ifunctopo,0)
					!$OMP CRITICAL
					inow=inow+1
					call showprog(inow,num3n1)
					!$OMP END CRITICAL
				else
					call findpath(i,1,ifunctopo,1)
				end if
			end if
		end do
		!$OMP END PARALLEL DO
		call sortpath
		write(*,"(' Totally found',i6,' new paths')") numpath-numpathold
		if ((numpath-numpathold)>0) idrawmol=0 !Avoid atoms and bonds blocking newly generated paths
		
!9999999999999999999
	else if (isel==9) then
		num3p1=count(CPtype(1:numcp)==3)
		if (num3p1==0) then
			write(*,"(a)") " Error: You must at least locate one (3,+1) CP before trying to generate topology paths!"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		if (count(CPtype(1:numcp)==4)==0) then
			write(*,"(a)") " Error: You must at least locate one (3,+3) CP before trying to generate topology paths!"
			write(*,*) "Press ENTER button to continue"
			read(*,*)
			cycle
		end if
		numpathold=numpath
		write(*,"(' Number of (3,+1) CPs:',i6,'    Generating topology paths...')") num3p1
		inow=0
		!$OMP PARALLEL DO SHARED(inow) PRIVATE(i) schedule(dynamic) NUM_THREADS(nthreads)
		do i=1,numcp
			if (CPtype(i)==3) then !Cycle each (3,+1)
				if (ishowsearchpath==0) then
					!$OMP CRITICAL
					inow=inow+1
					call showprog(inow,num3p1)
					!$OMP END CRITICAL
					call findpath(i,2,ifunctopo,0)
				else
					call findpath(i,2,ifunctopo,1)
				end if
			end if
		end do
		!$OMP END PARALLEL DO
		call sortpath
		write(*,"(' Totally found',i6,' new paths')") numpath-numpathold
		
!10 10 10 10 10 10 10
	else if (isel==10.and.count(CPtype(1:numcp)==2)==0) then
		write(*,*) "Error: You have to find at least one (3,-1) critical point"
	else if (isel==10) then
		do while(.true.)
			write(*,*)
			write(*,*) "Generate or delete interbasin surface for which (3,-1)?"
			write(*,*) "e.g. 5 means generate the surface from the (3,-1) with index of 5"
			if (numbassurf> 0) write(*,*) "     -3 means delete the surface from the (3,-1) with index of 3" 
			if (numbassurf> 0) write(*,*) "If input 0, surfaces from all (3,-1) will be deleted"
			if (numbassurf==0) write(*,*) "If input 0, surfaces from all (3,-1) will be generated"
			if (numbassurf> 0) write(*,*) "To list all generated surfaces, input the letter ""l"""
			if (numbassurf> 0) write(*,*) "To export the surfaces from the (3,-1) with index of 4, input ""o 4"""
			write(*,*) "To return, input ""q"""
			read(*,"(a)") c200

			if (c200(1:1)=='q') then
				exit
			else if (c200(1:1)=='l') then
				do isurf=1,numbassurf
					do icp=1,numcp
						if (cp2surf(icp)==isurf) exit
					end do
					write(*,"(' Index of interbasin surface:',i8,'  Index of corresponding (3,-1):',i8)") isurf,icp
				end do
			else if (c200(1:1)=='o') then
				read(c200(3:),*) icp
				if (icp>numcp.or.icp<=0) then
					write(*,*) "Error: The index of the surface is nonexistent"
				else if (cp2surf(icp)==0) then
					write(*,*) "Error: This CP is not (3,-1), input again"
				else
					open(10,file="surpath.txt",status="replace")
					do ipath=1,nsurfpathpercp
						write(10,"(' Path',i8)") ipath
						do ipt=1,nsurfpt
							write(10,"(i6,3f14.8)") ipt,bassurpath(:,ipt,ipath,cp2surf(icp))
						end do
					end do
					close(10)
					write(*,"(a)") "The coordinates of the paths of the surface have been exported to surpath.txt in current folder"
				end if
			else
				read(c200,*) isel2
				if (isel2==0) then
					if (numbassurf>0) then
						numbassurf=0
						cp2surf=0 !If cps2surf(i)==0, means the (3,-1) with total index of i hasn't been given surface
						deallocate(bassurpath)
					else if (numbassurf==0) then !Generate all interbasin surfaces
						numbassurf=count(CPtype(1:numcp)==2)
						allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
						isurf=1
						write(*,"(i8,' surfaces will be generated, please wait...')") numbassurf
						do icp=1,numcp
							if (CPtype(icp)==2) then
								cp2surf(icp)=isurf
								call genbassurf(icp,isurf,ifunctopo)
								isurf=isurf+1
								write(*,"(a,i6)") " Finished the surface generation from (3,-1) with index of",icp
							end if
						end do
					end if
				else if (abs(isel2)>numcp) then
					write(*,*) "Error: This CP is nonexisted"
				else if (CPtype(abs(isel2))/=2) then
					write(*,*) "This CP is not (3,-1), input again"
				else if (isel2>0) then !Add a surface
					if (cp2surf(isel2)/=0) then
						write(*,*) "The interbasin surface has already been generated" 
					else
						write(*,*) "Please wait..."
						if (numbassurf>0) then
							allocate(bassurpathtmp(3,nsurfpt,nsurfpathpercp,numbassurf))
							bassurpathtmp=bassurpath
							deallocate(bassurpath)
							numbassurf=numbassurf+1
							allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
							bassurpath(:,:,:,1:numbassurf-1)=bassurpathtmp(:,:,:,:)
							deallocate(bassurpathtmp)
						else if (numbassurf==0) then
							numbassurf=numbassurf+1
							allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
						end if
						call genbassurf(isel2,numbassurf,ifunctopo)
						cp2surf(isel2)=numbassurf
						write(*,*) "Done!"
					end if
				else if (isel2<0) then !Delete a surface
					idelsurf=cp2surf(abs(isel2))
					if (idelsurf==0) then
						write(*,*) "Surface has not been generated from this CP"
					else
						allocate(bassurpathtmp(3,nsurfpt,nsurfpathpercp,numbassurf))
						bassurpathtmp=bassurpath
						deallocate(bassurpath)
						numbassurf=numbassurf-1
						allocate(bassurpath(3,nsurfpt,nsurfpathpercp,numbassurf))
						bassurpath(:,:,:,1:idelsurf-1)=bassurpathtmp(:,:,:,1:idelsurf-1)
						bassurpath(:,:,:,idelsurf:)=bassurpathtmp(:,:,:,idelsurf+1:)
						deallocate(bassurpathtmp)
						cp2surf(abs(isel2))=0
						where(cp2surf>idelsurf) cp2surf=cp2surf-1
					end if
				end if
			end if
		end do
	else if (isel==20.and.ifunctopo==1) then
		do while(.true.)
            write(*,*)
			write(*,*) "Input the indices of the CPs in the ring, e.g. 22,23,25,28,32,11"
			write(*,*) "(Input q can exit)"
			read(*,"(a)") c200
			if (c200(1:1)=='q'.or.c200(1:1)=='Q') exit
			call str2arr(c200,nshanaromat,shanCPind)
			allocate(shanCPrho(nshanaromat))
			totdens=0D0
			do ishan=1,nshanaromat
				shanCPrho(ishan)=fdens(CPpos(1,shanCPind(ishan)),CPpos(2,shanCPind(ishan)),CPpos(3,shanCPind(ishan)))
			end do
			totrho=sum(shanCPrho(1:nshanaromat))
			shant=0D0
            write(*,*)
			do ishan=1,nshanaromat
				shanentropy=-shanCPrho(ishan)/totrho*log(shanCPrho(ishan)/totrho)
				write(*,"(' Electron density at CP',i3,':',f15.10,'  Local entropy:',f15.10)") shanCPind(ishan),shanCPrho(ishan),shanentropy
				shant=shant+shanentropy
			end do
			shanmax=log(dfloat(nshanaromat))
			write(*,"(' Total electron density:',f15.10)") totrho
			write(*,"(' Total Shannon entropy:',f15.10)") shant
			write(*,"(' Expected maximum Shannon entropy:',f15.10)") shanmax		
			write(*,*)
			write(*,"(' Shannon aromaticity index:',f16.10)") shanmax-shant 
			deallocate(shanCPrho)
		end do
	else if (isel==21.and.ifunctopo==1) then
		write(*,*) "Input XYZ coordinate of the point to be calculated, e.g. 2.0,2.4,1.1"
        write(*,*) "or input index of a CP, e.g. 4"
		read(*,"(a)") c200
		if ( index(c200,',')==0 .and. index(trim(c200),' ')==0 ) then
			read(c200,*) ithisCP
			tmpx=CPpos(1,ithisCP)
			tmpy=CPpos(2,ithisCP)
			tmpz=CPpos(3,ithisCP)
		else
			read(c200,*) tmpx,tmpy,tmpz
			write(*,*) "The coordinate you inputted is in which unit?  1=Bohr  2=Angstrom"
			read(*,*) iunit
			if (iunit==2) then
				tmpx=tmpx/b2a
				tmpy=tmpy/b2a
				tmpz=tmpz/b2a
			end if
		end if
        write(*,*) "Calculate gradient and curvature of electron density in which direction?"
        write(*,*) "1 Along the direction by specifying a vector"
        write(*,*) "2 Along the direction perpendicular to a given plane"
        read(*,*) idirtype
        if (idirtype==1) then
			write(*,*) "Input X,Y,Z of the direction vector, e.g. 3.2,1.0,5"
            read(*,*) xnor,ynor,znor
        else if (idirtype==2) then
			write(*,*) "Input indices of atoms to define a fitting plane, e.g. 3,4,9-12"
			write(*,*) "Evidently, at least three atoms should be inputted"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,ntmp)
			allocate(tmparr(ntmp))
			call str2arr(c2000tmp,ntmp,tmparr)
			call ptsfitplane(tmparr,ntmp,xnor,ynor,znor,planeD,rmsfit)
            deallocate(tmparr)
        end if
		facnorm=sqrt(xnor**2+ynor**2+znor**2) !Normalize the vector, then (xnor,ynor,znor) is the unit vector
		xnor=xnor/facnorm
		ynor=ynor/facnorm
		znor=znor/facnorm
		if (allocated(b)) then
			call gencalchessmat(2,1,tmpx,tmpy,tmpz,densvalue,gradtmp,hesstmp)
			densgrad=xnor*gradtmp(1)+ynor*gradtmp(2)+znor*gradtmp(3)
			denscurvature=xnor*xnor*hesstmp(1,1)+xnor*ynor*hesstmp(1,2)+xnor*znor*hesstmp(1,3)+&
						  ynor*xnor*hesstmp(2,1)+ynor*ynor*hesstmp(2,2)+ynor*znor*hesstmp(2,3)+&
						  znor*xnor*hesstmp(3,1)+znor*ynor*hesstmp(3,2)+znor*znor*hesstmp(3,3)
			write(*,"(' X,Y,Z of unit direction vector is',3f14.8)") xnor,ynor,znor
			write(*,"(' Electron density is          ',f30.10,' a.u.')") densvalue
			write(*,"(' Electron density gradient is ',f30.10,' a.u.')") densgrad
			write(*,"(' Electron density curvature is',f30.10,' a.u.')") denscurvature
		end if
		write(*,*)
		write(*,"(a)") " BTW: The X,Y,Z coordinates (row) of current point, the points below and above 1 Angstrom of the plane from current point, respectively (in Angstrom)"
		write(*,"(3f16.10)") tmpx*b2a,tmpy*b2a,tmpz*b2a
		write(*,"(3f16.10)") (tmpx-xnor/b2a)*b2a,(tmpy-ynor/b2a)*b2a,(tmpz-znor/b2a)*b2a
		write(*,"(3f16.10)") (tmpx+xnor/b2a)*b2a,(tmpy+ynor/b2a)*b2a,(tmpz+znor/b2a)*b2a
		write(*,*)
	end if
end do
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!--------- Generate interbasin surface from (3,-1)
subroutine genbassurf(ithisCP,ithissurf,ifunc)
use defvar
use functions
use util
use topo
implicit real*8 (a-h,o-z)
integer ifunc,ithisCP,ithissurf
real*8 initstpsize,hess(3,3),grad(3),eigvecmat(3,3),eigval(3),basvec1(3),basvec2(3),k1(3),k2(3)
initstpsize=surfpathstpsiz/4D0 !Smaller than pathstepsize(0.02)

call gencalchessmat(2,ifunc,CPpos(1,ithisCP),CPpos(2,ithisCP),CPpos(3,ithisCP),value,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)
!Generate normalized basis vectors from two negative eignvectors of hessian matrix
itmp=0
do i=1,3
	if (eigval(i)<0) then
		if (itmp==0) then
			rnorm=dsqrt(eigvecmat(1,i)**2+eigvecmat(2,i)**2+eigvecmat(3,i)**2)
			basvec1=eigvecmat(:,i)/rnorm
			itmp=1
		else if (itmp==1) then
			rnorm=dsqrt(eigvecmat(1,i)**2+eigvecmat(2,i)**2+eigvecmat(3,i)**2)
			basvec2=eigvecmat(:,i)/rnorm
		end if
	end if
end do
!Use the two basis vectors to generate a circle of initial points around (3,-1)
angstp=360D0/nsurfpathpercp
do i=1,nsurfpathpercp
	ang=(i-1)*angstp/180D0*pi !Convert to arc unit
	bassurpath(:,1,i,ithissurf)=initstpsize*( sin(ang)*basvec1(:)+cos(ang)*basvec2(:) ) + CPpos(:,ithisCP)
end do

!Generate gradient path from each initial point to comprise interbasin surface
do ipath=1,nsurfpathpercp
	do ipt=2,nsurfpt
		!Move point, RK2 method. Only calculate function value and gradient
		xtmp=bassurpath(1,ipt-1,ipath,ithissurf)
		ytmp=bassurpath(2,ipt-1,ipath,ithissurf)
		ztmp=bassurpath(3,ipt-1,ipath,ithissurf)
		call gencalchessmat(1,ifunc,xtmp,ytmp,ztmp,value,grad,hess)
		k1=grad/dsqrt(sum(grad**2))
		call gencalchessmat(1,ifunc,xtmp+surfpathstpsiz/2*k1(1),ytmp+surfpathstpsiz/2*k1(2),ztmp+surfpathstpsiz/2*k1(3),value,grad,hess)
		k2=grad/dsqrt(sum(grad**2))
		bassurpath(:,ipt,ipath,ithissurf)=bassurpath(:,ipt-1,ipath,ithissurf)-surfpathstpsiz*k2
	end do
end do
end subroutine



!!--------- Generate interbasin path from (3,-1) on a given plane
!ifunc=which real space function, ithisCP: The total index of (3,-1), ithispath: Store to which slot of ple3n1path array
subroutine gen3n1plepath(ifunc,ithisCP,ithispath)
use defvar
use topo
use functions
use util
implicit real*8 (a-h,o-z)
integer ifunc,ithisCP,ithispath
real*8 initstpsize,hess(3,3),grad(3),eigvecmat(3,3),eigval(3),initvec(3),plenormvec(3),k1(3),k2(3)
initstpsize=ple3n1pathstpsiz/4D0 !smaller than pathstepsize(0.02)

call gencalchessmat(2,ifunc,CPpos(1,ithisCP),CPpos(2,ithisCP),CPpos(3,ithisCP),value,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)
do iposvec=1,3
	if (eigval(iposvec)>0) exit
end do
if (plesel==1) then !XY
	plenormvec(1:2)=0D0
	plenormvec(3)=1D0
else if (plesel==2) then !XZ
	plenormvec(1:3)=0D0
	plenormvec(2)=1D0
else if (plesel==3) then !YZ
	plenormvec(2:3)=0D0
	plenormvec(1)=1D0
else
	call pointABCD(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,plenormvec(1),plenormvec(2),plenormvec(3),tmp)
end if

call vecprod(plenormvec(1),plenormvec(2),plenormvec(3), eigvecmat(1,iposvec),eigvecmat(2,iposvec),eigvecmat(3,iposvec), initvec(1),initvec(2),initvec(3))
rnorm=dsqrt(sum(initvec**2))
initvec=initvec/rnorm
ple3n1path(:,1,1,ithispath)=CPpos(:,ithisCP)+initvec(:)*initstpsize
ple3n1path(:,1,2,ithispath)=CPpos(:,ithisCP)-initvec(:)*initstpsize

!Generate gradient path from the given (3,-1) to comprise interbasin path on the plane
do idir=1,2
	do ipt=2,n3n1plept
		!Move point, RK2 method. Only calculate function value and gradient, don't calculate Hessian for saving time
		xtmp=ple3n1path(1,ipt-1,idir,ithispath)
		ytmp=ple3n1path(2,ipt-1,idir,ithispath)
		ztmp=ple3n1path(3,ipt-1,idir,ithispath)
		call gencalchessmat(1,ifunc,xtmp,ytmp,ztmp,value,grad,hess)
		k1=grad/dsqrt(sum(grad**2))
		call gencalchessmat(1,ifunc,xtmp+surfpathstpsiz/2*k1(1),ytmp+surfpathstpsiz/2*k1(2),ztmp+surfpathstpsiz/2*k1(3),value,grad,hess)
		k2=grad/dsqrt(sum(grad**2))
		ple3n1path(:,ipt,idir,ithispath)=ple3n1path(:,ipt-1,idir,ithispath)-surfpathstpsiz*k2
	end do
end do
end subroutine



!!!--------- Find path from a critical point at x,y,z
!itype=1: from (3,-1) to (3,-3)  =2: from (3,+1) to (3,+3)  =3: between (3,+1) and (3,-1) (not implemented)
!If info==1, print intermediate searching information, else do not print
subroutine findpath(ithisCP,itype,ifunc,info)
use topo
use functions
use util
implicit real*8 (a-h,o-z)
!Global variable "npathtry" is the maximum number of paths that can emerge from this CP
!In each time invoking this routine, for itype=1 and 2, search two times (npathtry=2 in fact is sufficient); for itype=3, search npathtry times
integer ifunc,info,ithisCP
real*8 grad(3),hess(3,3),k1(3),k2(3)
real*8 eigvecmat(3,3),eigval(3),tmpvec(3)
real*8,allocatable :: tmparr(:,:,:)
real*8 pathtmp(3,maxpathpt,npathtry) !XYZ coordinate of the points in the newly generated paths

!Determine eigenvalue and eigenvector of Hessian at initial point
call gencalchessmat(2,ifunc,CPpos(1,ithisCP),CPpos(2,ithisCP),CPpos(3,ithisCP),valueCP,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)

if (itype==1.or.itype==2) then !From (3,-1) to (3,-3), or from (3,+1) to (3,+3)
	pathtmp(:,1,1)=CPpos(:,ithisCP) !Set first point as input coordinate
	pathtmp(:,1,2)=CPpos(:,ithisCP)
	if (itype==1) then
		do iposi=1,3
			if (eigval(iposi)>0) exit !Find positive eigenvalue of Hessian of (3,-1)
		end do
	else
		do iposi=1,3
			if (eigval(iposi)<0) exit !Find negative eigenvalue of Hessian of (3,+1)
		end do
	end if
iterdir:	do idir=1,2
		if (info==1) then
			if (idir==1) write(*,"(' Go forward from CP: ',i6,1x,a,' Position:',3f11.5)") ithisCP,CPtyp2lab(CPtype(ithisCP)),CPpos(:,ithisCP)
			if (idir==2) write(*,"(' Go backward from CP:',i6,1x,a,' Position:',3f11.5)") ithisCP,CPtyp2lab(CPtype(ithisCP)),CPpos(:,ithisCP)
		end if
		posvecnorm=dsqrt(sum(eigvecmat(:,iposi)**2))
		if (idir==1) pathtmp(:,2,idir)=pathtmp(:,1,idir)+pathstepsize*eigvecmat(:,iposi)/posvecnorm !Move forwards along eigenvector with positive eigenvalue
		if (idir==2) pathtmp(:,2,idir)=pathtmp(:,1,idir)-pathstepsize*eigvecmat(:,iposi)/posvecnorm !Move backwards along eigenvector with positive eigenvalue
        if (ifPBC>0) call move_to_cell(pathtmp(:,2,idir),pathtmp(:,2,idir)) !If move to a position out of box, move it to central cell
        !Check if the path has already presented by comparing corresponding first two points
		do ickpath=1,numpath
			if ( dsqrt(sum( (topopath(:,1,ickpath)-pathtmp(:,1,idir))**2 ))<0.01D0 &
			.and.dsqrt(sum( (topopath(:,2,ickpath)-pathtmp(:,2,idir))**2 ))<0.01D0) then
				if (info==1) then
					write(*,*) "The path has already presented, skip search"
					write(*,*)
				end if
				cycle iterdir
			end if
		end do

		value=valueCP
iterpt:	do ipt=2,maxpathpttry
			!Check if the distance between current point and any other CP is smaller than threshold (discritpathfin)
			do icp=1,numcp
				if (icp==ithisCP) cycle
                if (ifPBC==0) then
    				distcp=dsqrt(sum( (pathtmp(:,ipt,idir)-CPpos(:,icp))**2 ))
                else !PBC case
                    call nearest_mirror(pathtmp(:,ipt,idir),CPpos(:,icp),tmpvec)
    				distcp=dsqrt(sum( (pathtmp(:,ipt,idir)-tmpvec(:))**2 ))
                end if
				if (distcp<discritpathfin) then !Terminal of this path is close enough to a CP, end the path generation
					!$OMP CRITICAL
					numpath=numpath+1
					pathnumpt(numpath)=ipt
					topopath(:,:,numpath)=pathtmp(:,:,idir)
                    !Add the nearly reached CP as the final point of the path. In the case of PBC, use the nearest image of the CP as this point
                    pathnumpt(numpath)=pathnumpt(numpath)+1
                    if (ifPBC==0) then
                        topopath(:,pathnumpt(numpath),numpath)=CPpos(:,icp)
                    else !PBC case
                        call nearest_mirror(pathtmp(:,ipt,idir),CPpos(:,icp),tmpvec)
                        topopath(:,pathnumpt(numpath),numpath)=tmpvec
                    end if
					!$OMP END CRITICAL
					if (info==1) then
						write(*,"(a,i6,a,f8.4)") " Found new path after",ipt," iterations, path length:",(ipt-1)*pathstepsize
						write(*,"(' Reached CP',i6,1x,a,'    Position:',3f12.6)") icp,CPtyp2lab(CPtype(icp)),CPpos(:,icp)
					end if
					exit iterpt
				end if
			end do
			if (ipt==maxpathpttry) then
				if (info==1) write(*,"(a)") " Path generation failed, upper limit of number of iterations is exceeded before reaching a CP"
				exit
			end if
			
			!Move point, RK2 method. Only calculate function value and gradient
			xtmp=pathtmp(1,ipt,idir)
			ytmp=pathtmp(2,ipt,idir)
			ztmp=pathtmp(3,ipt,idir)
			valueold=value
			call gencalchessmat(1,ifunc,xtmp,ytmp,ztmp,value,grad,hess)
            !For electron density, also check if value decrease during generating bond path
            !While for other functions like LOL, when path curvature at CP is large, this check often make path generation failed, &
            !because at initial a few steps, the value may marginally decrease, but in fact this doesn't matter
            if (ifunc==1.and.itype==1.and.value<valueold) then 
				if (info==1) write(*,"(a)") " Path generation failed, function value is lower than the last step but no located (3,-3) has reached"
				exit
            end if
			k1=grad/dsqrt(sum(grad**2))
			call gencalchessmat(1,ifunc,xtmp+pathstepsize/2*k1(1),ytmp+pathstepsize/2*k1(2),ztmp+pathstepsize/2*k1(3),value,grad,hess)
			k2=grad/dsqrt(sum(grad**2))
			if (itype==1) pathtmp(:,ipt+1,idir)=pathtmp(:,ipt,idir)+pathstepsize*k2
			if (itype==2) pathtmp(:,ipt+1,idir)=pathtmp(:,ipt,idir)-pathstepsize*k2
            if (ifPBC>0) then !If the point has moved to a position out of box and it is not on a plane (avoid jumping of boundary points), moving it to central cell
				call check_on_plane(pathtmp(:,ipt+1,idir),0.001D0,ionplane)
				if (ionplane==0) call move_to_cell(pathtmp(:,ipt+1,idir),pathtmp(:,ipt+1,idir))
            end if
		end do iterpt

		if (info==1) write(*,*)
	end do iterdir
	
else if (itype==3) then !between (3,+1) and (3,-1)
	!Not implemented
end if

end subroutine




!!--------- Determine path type and connected to which two CPs, 0 means unknown (not connected a found CP)
!ipathtype =0: other   =1: (3,-1)->(3,-3) =2: (3,+1)->(3,+3) =3: (3,-1)<-->(3,+1)
subroutine path_cp(ipath,icp1,icp2,ipathtype)
use defvar
use topo
implicit real*8 (a-h,o-z)
integer ipath,icp1,icp2
real*8 tmpvec(3)

iunknown=0
do icp1=1,numcp
    if (ifPBC==0) then
        tmpvec=CPpos(:,icp1)
    else if (ifPBC>0) then
        call nearest_mirror(topopath(:,1,ipath),CPpos(:,icp1),tmpvec)
    end if
	if (sum( (topopath(:,1,ipath)-tmpvec(:))**2 )<discritpathfin**2) exit !Test the first point in the path
	if (icp1==numcp) iunknown=1
end do
if (iunknown==1) icp1=0

iunknown=0
do icp2=1,numcp
    if (ifPBC==0) then
        tmpvec=CPpos(:,icp2)
    else if (ifPBC>0) then
        call nearest_mirror(topopath(:,pathnumpt(ipath),ipath),CPpos(:,icp2),tmpvec)
    end if
	if (sum( (topopath(:,pathnumpt(ipath),ipath)-tmpvec(:))**2 )<discritpathfin**2) exit !Test the last point
	if (icp2==numcp) iunknown=1
end do
if (iunknown==1) icp2=0
if (icp1==0.or.icp2==0) then
	ipathtype=0
else
	icp1type=CPtype(icp1)
	icp2type=CPtype(icp2)
	if (icp1type==2.and.icp2type==1) then
		ipathtype=1
	else if (icp1type==3.and.icp2type==4) then
		ipathtype=2
	else if ( (icp1type==2.and.icp2type==3).or.(icp1type==3.and.icp2type==2) ) then
		ipathtype=3
	else
		ipathtype=0
	end if
end if
end subroutine



!!!----------- Find critical points from initial guess at X,Y,Z using Newton method
!ifunc is index of real space functions
!ishowsearchlevel=0/1/2/3:  Print none/minor/some/all detail. Notice that in parallel mode, the outputted details are messed up
subroutine findcp(x,y,z,ifunc)
use topo
use functions
use util
implicit real*8 (a-h,o-z)
integer ifunc
real*8 x,y,z
real*8 coord(3,1),grad(3,1),hess(3,3),disp(3,1),disptmp(3),gvec(3),gvec_old(3),dvec(3),tmpvec(3)
real*8 eigvecmat(3,3),eigval(3) !,tmpmat(3,3)

coord(1,1)=x
coord(2,1)=y
coord(3,1)=z
if (ifPBC>0) call move_to_cell(coord,coord) !Initial position may already be out of cell, move it to central cell
if (ishowsearchlevel>1) write(*,"(' Starting point:',3f12.6)") coord(1:3,1)

do i=1,topomaxcyc
	if (ishowsearchlevel==3) then
		write(*,"(/,' Iteration',i5)") i
        write(*,"(' Coordinate:',3f16.10)") coord
    end if
    
	!Determine function value and gradient at present geometry, and calculate raw displacement
    if (itopomethod==1) then !Newton method
	    call gencalchessmat(2,ifunc,coord(1,1),coord(2,1),coord(3,1),value,grad(1:3,1),hess) !Obtain gradient and Hessian
	    singulartest=abs(detmat(hess))
	    if (singulartest<singularcrit) then
		    if (ishowsearchlevel>1) then
			    write(*,*) "Hessian matrix is singular at current position, stop iteration"
			    write(*,"(' Absolute of determinant of Hessian matrix:',1PE18.10)") singulartest
			    write(*,"(' Criterion for detecting singular:',1PE18.10)") singularcrit
		    end if
		    exit
	    end if
	    disp=-matmul(invmat(hess,3),grad)
    else if (itopomethod==2) then !Barzilai-Borwein method
		!Use BB2 (see "STABILIZED BARZILAI-BORWEIN METHOD") to determine stepsize, corresponding to Barzilai-Borwein in https://en.wikipedia.org/wiki/Gradient_descent
		!Although this method was originally designed for searching minimum, in fact it mimics Newton method, it converges to stationary point (not necessarily minimum or maximum)
		!In otherwords, BB2 is somewhat like quasi-Newton method, but Hessian is not explicitly needed to constructed and updated
        gvec_old=gvec
        call gencalchessmat(1,ifunc,coord(1,1),coord(2,1),coord(3,1),value,gvec(:),hess) !Obtain gradient
        if (i>1) then 
            dvec(:)=disp(:,1)
            val1=sum(dvec(:)*(gvec-gvec_old))
            val2=sum((gvec-gvec_old)**2)
            disp(:,1)=-val1/val2*gvec(:)
        else if (i==1) then !Steepest descent with small fixed step at first step
            stepinit=0.01D0
            gvecnorm=dsqrt(sum(gvec(:)**2))
            if (gvecnorm<1D-20) then
				disp=0
            else
				disp(:,1)=-gvec(:)*stepinit/gvecnorm
            end if
        end if
        grad(:,1)=gvec(:)
    else if (itopomethod==3.or.itopomethod==4) then !Steepest ascent/descent
        call gencalchessmat(1,ifunc,coord(1,1),coord(2,1),coord(3,1),value,gvec(:),hess(:,:)) !Obtain gradient
		gradnorm=dsqrt(sum(gvec**2))
		tmpvec=gvec(:)/gradnorm !Unit vector along gradient
		sclfac=0.1D0 !First stepsize
        micromax=30
        if (ishowsearchlevel==3) then
			write(*,*) "** Start micro iterations of line search"
			write(*,"(' Initial function value:',f24.13)") value
			write(*,"(' Gradient:     ',3f16.10)") gvec(:)
            write(*,*) "Iter.             X,Y,Z of attempt displacement                 New value"
        end if
		do imicro=1,micromax
			if (itopomethod==3) then !Steepest ascent
				disptmp(:)=sclfac*tmpvec(:)
            else !Steepest descent
				disptmp(:)=-sclfac*tmpvec(:)
            end if
			xtest=coord(1,1)+disptmp(1)
			ytest=coord(2,1)+disptmp(2)
			ztest=coord(3,1)+disptmp(3)
			tmpval=calcfuncall(ifunc,xtest,ytest,ztest)
            if (ishowsearchlevel==3) write(*,"(i3,3f19.14,1PE18.10)") imicro,sclfac*tmpvec(:),tmpval
			if ((itopomethod==3.and.tmpval>value).or.(itopomethod==4.and.tmpval<value)) then
				if (ishowsearchlevel==3) write(*,*) "** Displacement accepted!"
				disp(:,1)=disptmp(:)
				exit
			else
				sclfac=sclfac/2.5D0
			end if
		end do
		if (imicro==micromax+1) then !Usually when line search was failed, displacement is already quite small, while gradient is not quite small
			disp(:,1)=0D0
            gvec(:)=0D0
			if (ishowsearchlevel==3) write(*,*) "Warning: Line search was failed, assumed to be converged"
		end if
        grad(:,1)=gvec(:)
    end if
    
	gradnorm=dsqrt(sum(grad**2))
	if (ishowsearchlevel==3) then
		write(*,"(' Function value:',f24.13)") value
		write(*,"(' Grad:',3f16.10,'  Norm:',f18.10)") grad,gradnorm
	end if
    
    !Apply scale factor
    disp=CPstepscale*disp
    !Apply trust radius
    if (topotrustrad>0) then
        dispnorm=dsqrt(sum(disp**2))
        if (dispnorm>topotrustrad) disp=disp*topotrustrad/dispnorm
    end if
	dispnorm=dsqrt(sum(disp**2))
    
    !Update coordinate
	coord=coord+disp
    if (ifPBC>0) call move_to_cell(coord,coord) !If moved to a position out of box, move it to central cell
    
	if (ishowsearchlevel==3) then
		write(*,"(' Disp:',3f16.10,'  Norm:',f18.10)") disp,dispnorm
		write(*,"(' Goal: |disp|<',E18.8,'    |grad|<',E18.8)") dispconv,gradconv
        if (dispnorm>dispconv.or.gradnorm>gradconv) write(*,*) "Not converged"
	end if
    
	if (dispnorm<dispconv.and.gradnorm<gradconv) then
		if (ishowsearchlevel>1) write(*,"(' Converged after',i6,' iterations')") i
		if (ishowsearchlevel==3) write(*,"(/,a)") "         ---------------------- Iteration ended ----------------------"
		if (itopomethod==2.or.itopomethod==3.or.itopomethod==4) then !When using methods other than Newton, only gradient is calculated, here we calculate Hessian for determining CP type
			call gencalchessmat(2,ifunc,coord(1,1),coord(2,1),coord(3,1),value,gvec(:),hess)
		end if
        !$OMP CRITICAL
		inewcp=1
		do icp=1,numcp
            if (ifPBC==0) then
    				r=dsqrt(sum( (coord(:,1)-CPpos(:,icp))**2 ))
            else
                call nearest_dist(coord(:,1),CPpos(:,icp),r)
            end if
			if (r<=minicpdis) then
				if (ishowsearchlevel>1) write(*,"(a,i6,a)") " This CP is too close to CP",icp,", ignored"
				inewcp=0
				exit
			end if
		end do
		if (CPsearchlow/=CPsearchhigh.and.(value<CPsearchlow.or.value>CPsearchhigh)) then !Check value of CP
			if (ishowsearchlevel>1) write(*,"(a,1PE12.5,a)") " The value of this CP is ",value,", which exceeded user-defined range and thus ignored"
			inewcp=0
		end if
        if (inewcp==1) then !Check Hessian of CP
            call diagsymat(hess,eigvecmat,eigval,idiagok)
			igt0=count(eigval>0)
            if (ishowsearchlevel>1) write(*,"(' Eigenvalues:',3(1PE20.10))") eigval
 			if (idiagok/=0) then
				if (ishowsearchlevel>1) write(*,*) "Note: Diagonization of Hessian matrix failed! This CP is ignored"
				inewcp=0
			end if
			if (all(eigval==0)) then !Occur in very rare case, such as (0,0,0) of thiophene-pi.mwfn of IRI official tutorial
				if (ishowsearchlevel>1) write(*,*) "All eigenvalues of Hessian are zero, this CP is ignored"
				inewcp=0
			end if
        end if
        if (itopomethod==3) then !Steepest ascent only accepts (3,-3)
			if (igt0/=0) then
				if (ishowsearchlevel>1) write(*,*) "Steepest ascent method only accepts (3,-3) CP but this is not, so ignored"
				inewcp=0
            end if
        else if (itopomethod==4) then !Steepest descent only accepts (3,+3)
			if (igt0/=3) then
				inewcp=0
				if (ishowsearchlevel>1) write(*,*) "Steepest descent method only accepts (3,+3) CP but this is not, so ignored"
            end if
        end if
		if (inewcp==1) then !Finally add CP to list
			numcp=numcp+1
			CPpos(:,numcp)=coord(:,1)
			if (igt0==3) then
				if (ishowsearchlevel>1) write(*,"(' Found new (3,+3) at',3f15.10)") coord
				CPtype(numcp)=4
			else if (igt0==2) then
				if (ishowsearchlevel>1) write(*,"(' Found new (3,+1) at',3f15.10)") coord
				CPtype(numcp)=3
			else if (igt0==1) then
				if (ishowsearchlevel>1) write(*,"(' Found new (3,-1) at',3f15.10)") coord
				CPtype(numcp)=2
				call sort(eigval)
				if (ishowsearchlevel>1) write(*,"(' Bond ellipticity is',f15.10)") eigval(1)/eigval(2)-1D0
			else if (igt0==0) then
				if (ishowsearchlevel>1) write(*,"(' Found new (3,-3) at',3f15.10)") coord
				CPtype(numcp)=1
			end if
		end if
        !$OMP end CRITICAL
		exit
	end if
	if (i==topomaxcyc.and.(ishowsearchlevel>1)) write(*,"(/,a)") "          !! Exceeded maximal cycles until find a stationary point !!"
end do
if (ishowsearchlevel>1) write(*,*)
end subroutine



!!----- Sort newly found CPs according to coordinates. This is mainly used to garantee that the CP indices are identical in each time of execution under parallel mode
!numcpoldp1: The number of CPs before this search + 1
subroutine sortCP(numcpoldp1)
use topo
implicit real*8 (a-h,o-z)
integer numcpoldp1,typetmp
real*8 tmparr(3)
do itmp=numcpoldp1,numcp
	do jtmp=itmp+1,numcp
		tmpvali=CPpos(1,itmp)*0.234134D0+CPpos(2,itmp)*1.9837322D0-CPpos(3,itmp)*0.5413578924D0 !Use three random number to generate unique code
		tmpvalj=CPpos(1,jtmp)*0.234134D0+CPpos(2,jtmp)*1.9837322D0-CPpos(3,jtmp)*0.5413578924D0 
		if (tmpvali>tmpvalj) then
			typetmp=CPtype(itmp)
			CPtype(itmp)=CPtype(jtmp)
			CPtype(jtmp)=typetmp
			tmparr=CPpos(:,itmp)
			CPpos(:,itmp)=CPpos(:,jtmp)
			CPpos(:,jtmp)=tmparr
		end if
	end do
end do
end subroutine



!!----- Sort newly generated paths. This is mainly used to garantee that the path indices are identical in each time of execution under parallel mode
subroutine sortpath
use topo
implicit real*8 (a-h,o-z)
integer cp1(numpath),cp2(numpath)
real*8 tmparr(3,maxpathpt)
do ipath=1,numpath
	call path_cp(ipath,cp1(ipath),cp2(ipath),ipathtype)
end do
!Sorting according to index of starting CP of the paths, from small to large
do i=1,numpath
	do j=i+1,numpath
		if (cp1(i)>cp1(j)) then
			ntmp=pathnumpt(i)
			pathnumpt(i)=pathnumpt(j)
			pathnumpt(j)=ntmp
			tmparr=topopath(:,:,i)
			topopath(:,:,i)=topopath(:,:,j)
			topopath(:,:,j)=tmparr
			itmp=cp1(i)
			cp1(i)=cp1(j)
			cp1(j)=itmp
			itmp=cp2(i)
			cp2(i)=cp2(j)
			cp2(j)=itmp
		end if
	end do
end do
!Sorting according to index of ending CP of the paths, from small to large
do i=1,numpath
	do j=i+1,numpath
		if (cp1(i)==cp1(j).and.cp2(i)>cp2(j)) then
			ntmp=pathnumpt(i)
			pathnumpt(i)=pathnumpt(j)
			pathnumpt(j)=ntmp
			tmparr=topopath(:,:,i)
			topopath(:,:,i)=topopath(:,:,j)
			topopath(:,:,j)=tmparr
		end if
	end do
end do
end subroutine



!!------- Find corresponding atomic index of NCPs (icp)
!iatm is corresponding atom index, if return 0, that means this NCP doesn't corresponds to a nuclear
!If NCP and an atom is less than 0.15 Angstrom, then they will be regarded as corresponded
!The threshold should not be too small. For example, NCP and nuclear position of the hydrogen connected to nitrogen may be as large 0.12 Angstrom
subroutine NCP2atmidx(icp,iatm)
use defvar
use topo
implicit real*8 (a-h,o-z)
real*8 xyzA(3),xyzB(3)
xyzA(1)=CPpos(1,icp)
xyzA(2)=CPpos(2,icp)
xyzA(3)=CPpos(3,icp)
do iatm=1,ncenter
    xyzB(1)=a(iatm)%x
    xyzB(2)=a(iatm)%y
    xyzB(3)=a(iatm)%z
    call nearest_dist(xyzA,xyzB,dist)
	if (dist<0.15D0/b2a) exit
end do
if (iatm==ncenter+1) iatm=0
end subroutine



!!------- Return the two atoms connecting to the given BCP (iBCP) according to generated bond path
!iatm/jatm is corresponding atom index; if return 0, that means the atom connecting to the BCP cannot be determined
!Rule: Staring point of the bond path emerges from the BCP, and ending point is close enough to a nucleus 
subroutine BCP2atmidx(iBCP,iatm,jatm)
use defvar
use topo
implicit real*8 (a-h,o-z)
integer iBCP,iatm,jatm
iatm=0
jatm=0
do ipath=1,numpath
	call path_cp(ipath,icp1,icp2,ipathtype) !Find indices of the two CPs that connected by this path
    if (icp1==iBCP.and.ipathtype==1) then !Staring from the iBCP, and path is (3,-1)->(3,-3)
		if (iatm==0) then
			call NCP2atmidx(icp2,iatm)
        else
			call NCP2atmidx(icp2,jatm)
        end if
    end if
    if (iatm/=0.and.jatm/=0) exit
end do
end subroutine



!!-------- Only retain bond paths (and corresponding BCPs) connecting &
! two specific molecular fragments while remove all others bond paths and BCPs
subroutine retainpath
use util
use defvar
use topo
implicit real*8 (a-h,o-z)
character c2000tmp*2000,selectyn
integer delCP(numcp),delpath(numpath),ndelcp,ndelpath
write(*,*) "Input atom list for fragment 1, e.g. 3,5-8,15-20"
write(*,*) "Note: Input 0 can exit"
read(*,"(a)") c2000tmp
if (c2000tmp(1:1)=="0") return
call str2arr(c2000tmp,nterm1)
if (allocated(frag1)) deallocate(frag1)
allocate(frag1(nterm1))
call str2arr(c2000tmp,nterm1,frag1)
write(*,*) "Input atom list for fragment 2, e.g. 1,2,4,9-14,18"
read(*,"(a)") c2000tmp
call str2arr(c2000tmp,nterm2)
if (allocated(frag2)) deallocate(frag2)
allocate(frag2(nterm2))
call str2arr(c2000tmp,nterm2,frag2)

ndelcp=0
ndelpath=0
write(*,*) "Please wait..."
do icp=1,numcp !Cycle all (3,-1)
    call showprog(icp,numcp)
	if (CPtype(icp)/=2) cycle
	ifound=0
	do ipath=1,numpath !Find path that emitted from icp
		call path_cp(ipath,icp1,icp2,ipathtype) !Find indices of the two CPs that connected by this path
		if (icp1==icp) then
			ifound=ifound+1
			if (ifound==1) then
				ipath1=ipath 
				call NCP2atmidx(icp2,iatm1)
			else if (ifound==2) then
				ipath2=ipath
				call NCP2atmidx(icp2,iatm2)
				exit
			end if
		end if
	end do
	if ( (any(frag1==iatm1).and.any(frag2==iatm2)) .or. (any(frag1==iatm2).and.any(frag2==iatm1)) ) then
		continue !Retain this path
	else
		ndelcp=ndelcp+1
		delcp(ndelcp)=icp
        !In rare case, the BCP between the two fragments may only have one path because (3,-3) in atom of one fragment is missing or generation of path is failed
        !It is also possible that a BCP is not linked to two fragments simultaneously. Therefore conditions must be applied
        if (ifound==1.or.ifound==2) then
		    ndelpath=ndelpath+1
		    delpath(ndelpath)=ipath1
        end if
        if (ifound==2) then
		    ndelpath=ndelpath+1
		    delpath(ndelpath)=ipath2
        end if
	end if
end do
deallocate(frag1,frag2)

do idx=1,ndelpath
	ipath=delpath(idx)
    maxpt=maxval(pathnumpt)
	topopath(:,1:maxpt,ipath:numpath-1)=topopath(:,1:maxpt,ipath+1:numpath)
	pathnumpt(ipath:numpath-1)=pathnumpt(ipath+1:numpath)
    do jdx=idx+1,ndelpath
        if (delpath(jdx)>ipath) delpath(jdx)=delpath(jdx)-1
    end do
	numpath=numpath-1
end do
write(*,"(i6,' bond paths have been deleted, also delete corresponding BCPs? (y/n)')") ndelpath
read(*,*) selectyn
if (selectyn=='y'.or.selectyn=='Y') then
	do idx=1,ndelCP
		icp=delCP(idx)
		CPpos(:,icp:numcp-1)=CPpos(:,icp+1:numcp)
		CPtype(icp:numcp-1)=CPtype(icp+1:numcp)
		where(delCP(idx+1:ndelCP)>icp) delCP(idx+1:ndelCP)=delCP(idx+1:ndelCP)-1
		numcp=numcp-1
	end do
	write(*,"(i6,' BCPs have been deleted')") ndelCP
end if
end subroutine



!!--------- Plot real space function along topology path
subroutine plotpathprop
use defvar
use topo
use plot
use functions
implicit real*8 (a-h,o-z)
character c200*200

do while(.true.)
	write(*,*)
	write(*,*) "Input index of a path, e.g. 3     Input ""q"" can return"
	write(*,"(a)") " Hint: If input index of two paths (e.g. 6,7) emitted from the same (3,-1) CP, &
	&then the real space function along the combined paths will be outputted"
	if (ifunctopo==1) write(*,"(a)") " Hint: If you input e.g. c6 and meantime the CP6 is BCP, then the two bond paths emitted from it will be chosen"
	if (allocated(MOsym)) write(*,"(a)") " Hint: You can input ""s"" to choose which irreducible &
	&representations will be taken into account in the real space function evaluation"
	read(*,"(a)") c200
					
	if (index(c200,'q')/=0) then
		return
	else if (index(c200,'s')/=0) then
		call selMO_IRREP
		cycle
	end if
					
	itwopath=0
	if (index(c200,'c')==0.and.index(c200,'C')==0) then
		if (index(c200,',')==0) then !Only one path
			read(c200,*) ipath
		else !Two paths
			itwopath=1
			read(c200,*) ipath,jpath
		end if
		if (ipath>numpath.or.ipath<=0 .or. (itwopath==1.and.(jpath>numpath.or.jpath<=0)) ) then
			write(*,*) "Error: Invalid index"
			cycle
		end if
		if (itwopath==1.and.( topopath(1,1,ipath)/=topopath(1,1,jpath) .or. &
			topopath(2,1,ipath)/=topopath(2,1,jpath) .or. topopath(3,1,ipath)/=topopath(3,1,jpath) )) then
			write(*,*) "Error: The two paths are not emitted from the same (3,-1) critical point!"
			cycle
		end if
	else
		read(c200(2:),*) iCP
		if (CPtype(iCP)/=2) then
			write(*,*) "Error: This CP is not a BCP!"
			cycle
		end if
		ipath=0;jpath=0
		do kpath=1,numpath
			call path_cp(kpath,icp1,icp2,ipathtype)
			if (ipathtype==1) then
				if (icp1==iCP.or.icp2==iCP) then
					if (ipath==0) then
						ipath=kpath
					else
						jpath=kpath
						exit
					end if
				end if
			end if
		end do
		if (kpath==numpath+1) then
			write(*,*) "Error: Unable to find two bond paths emitted from the BCP you chosen!"
			cycle
		else
			write(*,"(' Index of the two bond paths is',2i6)") ipath,jpath
			itwopath=1
			write(*,*)
		end if
	end if
	
	write(*,*) "Select the real space function to calculate along the path"
	call selfunc_interface(1,iselfunc)
	!Store the values along the path into curvex and curvey, show them as text and curve map
	npointcurve=pathnumpt(ipath)
	if (itwopath==1) npointcurve=pathnumpt(ipath)+pathnumpt(jpath)-1
	if (allocated(curvex)) deallocate(curvex)
	if (allocated(curvey)) deallocate(curvey)
	allocate(curvex(npointcurve),curvey(npointcurve))
	if (itwopath==0) then
		do ipt=1,pathnumpt(ipath)
            if (ipt==1) then
                curvex(ipt)=0
            else
			    curvex(ipt)=curvex(ipt-1) + dsqrt(sum( (topopath(:,ipt,ipath)-topopath(:,ipt-1,ipath))**2) )
            end if
			curvey(ipt)=calcfuncall(iselfunc,topopath(1,ipt,ipath),topopath(2,ipt,ipath),topopath(3,ipt,ipath))
		end do
	else if (itwopath==1) then
		ipttmp=0
		do ipt=pathnumpt(ipath),1,-1
			ipttmp=ipttmp+1
            if (ipt==pathnumpt(ipath)) then
                curvex(ipttmp)=0
            else
			    curvex(ipttmp)=curvex(ipttmp-1) + dsqrt(sum( (topopath(:,ipt,ipath)-topopath(:,ipt+1,ipath))**2) )
            end if
			curvey(ipttmp)=calcfuncall(iselfunc,topopath(1,ipt,ipath),topopath(2,ipt,ipath),topopath(3,ipt,ipath))
		end do
		do jpt=2,pathnumpt(jpath) !The first point of jpath is (3,-1), which has been included in ipath above
			ipttmp=ipttmp+1
			curvex(ipttmp)=curvex(ipttmp-1) + dsqrt(sum( (topopath(:,jpt,jpath)-topopath(:,jpt-1,jpath))**2) )
			curvey(ipttmp)=calcfuncall(iselfunc,topopath(1,jpt,jpath),topopath(2,jpt,jpath),topopath(3,jpt,jpath))
		end do
	end if	
	write(*,*) " Index           X/Y/Z Coordinate (Bohr)          Dist.       Value"
	if (itwopath==0) then
		do ipt=1,pathnumpt(ipath)
			write(*,"(i6,3f14.8,f9.4,E16.8)") ipt,topopath(:,ipt,ipath),curvex(ipt),curvey(ipt)
		end do
		icurve_vertlinex=0
	else if (itwopath==1) then
		ipttmp=0
		do ipt=pathnumpt(ipath),1,-1
			ipttmp=ipttmp+1
			write(*,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,ipt,ipath),curvex(ipttmp),curvey(ipttmp)
		end do
		ibcptmp=ipttmp
		do jpt=2,pathnumpt(jpath) !The first point of jpath is (3,-1), which has been included in ipath above
			ipttmp=ipttmp+1
			write(*,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,jpt,jpath),curvex(ipttmp),curvey(ipttmp)
		end do
		write(*,"(' Note: The point',i5,' corresponds to (3,-1) critical point')") ibcptmp
		icurve_vertlinex=1 !Draw a vertical line to highlight BCP point
		curve_vertlinex=curvex(ibcptmp)
		write(*,*) "The dashed line on the graph corresponds to the position of BCP"
	end if
	!Show the result as curve map
	if (minval(curvey)>=0) then
		curveymin=0
		curveymax=1.1D0*maxval(curvey)
		steplaby=curveymax/11
	else
		exty=(maxval(curvey)-minval(curvey))/10
		curveymin=minval(curvey)-exty
		curveymax=maxval(curvey)+exty
		steplaby=exty
	end if
	steplabx=0.5D0
	ilog10y=0
	if (isilent==0) call drawcurve(curvex,curvey,npointcurve,0D0,maxval(curvex),steplabx,curveymin,curveymax,steplaby,"show")
	do while(.true.)
		write(*,*)
		write(*,*) "0 Return"
		write(*,*) "1 Plot the graph again"
		write(*,*) "2 Save the graph in current folder"
		if (ilog10y==0) then
			write(*,"(' 3 Set range of Y-axis, current:',f13.5,' to',f13.5,' Step:',f12.5)") curveymin,curveymax,steplaby
			write(*,*) "4 Switch Y-axis type, current: linear scaling"
		else if (ilog10y==1) then
			write(*,"(' 3 Set range of Y-axis, current:',1PE14.6,' to',1PE14.6)") 10**curveymin,10**curveymax
			write(*,*) "4 Switch Y-axis type, current: logarithmic scaling"
		end if
		write(*,*) "5 Export the data of the path shown above as pathvalue.txt in current folder"
		read(*,*) iselcurve
						
		if (iselcurve==0) then
			exit
		else if (iselcurve==1) then
			call drawcurve(curvex,curvey,npointcurve,0D0,maxval(curvex),steplabx,curveymin,curveymax,steplaby,"show")
		else if (iselcurve==2) then
			call drawcurve(curvex,curvey,npointcurve,0D0,maxval(curvex),steplabx,curveymin,curveymax,steplaby,"save")
		else if (iselcurve==3) then
			if (ilog10y==0) then
				write(*,*) "Input lower and upper limits and step between labels, e.g. 0,1.5,0.2"
				read(*,*) curveymin,curveymax,steplaby
			else if (ilog10y==1) then
				write(*,*) "Input minimum and maximum value of Y axis  e.g. -1,5 means from 10^-1 to 10^5"
				read(*,*) curveymin,curveymax
			end if
		else if (iselcurve==4) then
			if (ilog10y==0) then
				ilog10y=1
				curveymin=log10(curveymin)
				curveymax=log10(curveymax)
			else
				ilog10y=0
				curveymin=10**curveymin
				curveymax=10**curveymax
			end if
		else if (iselcurve==5) then				
			open(10,file="pathvalue.txt",status="replace")
			if (itwopath==0) then
				do ipt=1,pathnumpt(ipath)
					write(10,"(i6,3f14.8,f9.4,E16.8)") ipt,topopath(:,ipt,ipath),curvex(ipt),curvey(ipt)
				end do
			else if (itwopath==1) then
				ipttmp=0
				do ipt=pathnumpt(ipath),1,-1
					ipttmp=ipttmp+1
					write(10,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,ipt,ipath),curvex(ipttmp),curvey(ipttmp)
				end do
				do jpt=2,pathnumpt(jpath) !The first point of jpath is (3,-1), which has been included in ipath above
					ipttmp=ipttmp+1
					write(10,"(i6,3f14.8,f9.4,E16.8)") ipttmp,topopath(:,jpt,jpath),curvex(ipttmp),curvey(ipttmp)
				end do
			end if
			close(10)
			write(*,*) "Done! The data has been outputted to pathvalue.txt in current folder"
		end if
	end do
	icurve_vertlinex=0
	deallocate(curvex,curvey)
end do
end subroutine



!!----------- Get length of a path
subroutine getpathlength(ipath,pathlength)
use topo
use defvar
implicit real*8 (a-h,o-z)
integer ipath
pathlength=0
do ipt=2,pathnumpt(ipath)
    if (ifPBC==0) then
        disp=dsqrt(sum( (topopath(:,ipt,ipath)-topopath(:,ipt-1,ipath))**2 ))
    else !PBC
        call nearest_dist(topopath(:,ipt,ipath),topopath(:,ipt-1,ipath),disp)
    end if
    pathlength=pathlength+disp
end do
end subroutine



!!---------- Show topology information summary
subroutine showtoposummary
use defvar
use topo
implicit real*8 (a-h,o-z)
character icp1text*12,icp2text*12

if (numpath>0) then
	write(*,*) "Summary of generated topology paths:"
	do i=1,numpath
		call path_cp(i,icp1,icp2,ipathtype)
		if (icp1==0) then
			icp1text="   Unknown  "
		else
			write(icp1text,"(i5,1x,a)") icp1,CPtyp2lab(CPtype(icp1))
		end if
		if (icp2==0) then
			icp2text="  Unknown   "
		else
			write(icp2text,"(i5,1x,a)") icp2,CPtyp2lab(CPtype(icp2))
		end if
        call getpathlength(i,pathlength)
		write(*,"(' Path',i5,',',4x,'CP:',a,' -->',' CP:',a,'   Length:',f9.5,' Bohr')") i,icp1text,icp2text,pathlength
	end do
else
	write(*,*) "No path has been generated"
end if
write(*,*)
if (numcp>0) then
	write(*,*) "Summary of found CPs:"
	write(*,*) " Index            XYZ Coordinate (Bohr)            Type"
	do icp=1,numcp
		icptype=CPtype(icp)
		if (ifunctopo==1.and.icptype==1) then !NCP of AIM
			call NCP2atmidx(icp,iatm)
			if (iatm==0) then
				write(*,"(i6,3f14.8,3x,a,'   Nucleus: Unknown')") icp,CPpos(:,icp),CPtyp2lab(icptype)
			else
				write(*,"(i6,3f14.8,3x,a,'   Nucleus:',i5,'(',a')')") icp,CPpos(:,icp),CPtyp2lab(icptype),iatm,a(iatm)%name
			end if
		else if (ifunctopo==1.and.icptype==2) then !BCP of AIM
			call BCP2atmidx(icp,iatm,jatm)
			if (iatm==0.or.jatm==0) then
				write(*,"(i6,3f14.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(icptype)
			else
				write(*,"(i6,3f14.8,3x,a,i5,'(',a')',' --',i5,'(',a')')") icp,CPpos(:,icp),CPtyp2lab(icptype),iatm,a(iatm)%name,jatm,a(jatm)%name
			end if
		else
			write(*,"(i6,3f14.8,3x,a)") icp,CPpos(:,icp),CPtyp2lab(icptype)
		end if
	end do
    call showPHrelat
	if (numbassurf>0) write(*,"(' The number of generated interbasin surfaces:',i8)") numbassurf
else
	write(*,*) "No CP has been found"
end if
end subroutine



!!---------- Show satisfaction of Poincare-Hopf relationship
subroutine showPHrelat
use defvar
use topo
implicit real*8 (a-h,o-z)
NumCPtype1=count(CPtype(1:numcp)==1)
NumCPtype2=count(CPtype(1:numcp)==2)
NumCPtype3=count(CPtype(1:numcp)==3)
NumCPtype4=count(CPtype(1:numcp)==4)
write(*,*) "The number of critical points of each type:"
write(*,"(' (3,-3):',i6,',   (3,-1):',i6,',   (3,+1):',i6,',   (3,+3):',i6)") NumCPtype1,NumCPtype2,NumCPtype3,NumCPtype4
itestPH=NumCPtype1-NumCPtype2+NumCPtype3-NumCPtype4 !Poincare-Hopf relationship
write(*,"(' Poincare-Hopf verification:',i6,'  -',i6,'  +',i6,'  -',i6,'  =',i6)") NumCPtype1,NumCPtype2,NumCPtype3,NumCPtype4,itestPH
if (ifPBC==0) then
	if (itestPH==1) then
        write(*,*) "Fine, Poincare-Hopf relationship is satisfied, all CPs may have been found"
    else
        write(*,*) "Warning: Poincare-Hopf relationship is not satisfied, some CPs may be missing"
    end if
else if (ifPBC>0) then
	if (itestPH==0) then
        write(*,*) "Fine, Poincare-Hopf relationship is satisfied, all CPs may have been found"
    else
        write(*,*) "Warning: Poincare-Hopf relationship is not satisfied, some CPs may be missing"
    end if
end if
end subroutine



!!--------- Clean topology analysis information
subroutine deallo_topo
use topo
numcp=0
numcp_tmp=0
numpath=0
numpath_tmp=0
nple3n1path=0
numbassurf=0
CPtype=0 !Clean relationship
cp2surf=0
cp2ple3n1path=0
CPsearchlow=0D0
CPsearchhigh=0D0
lab_oneCP=0
if (allocated(bassurpath)) deallocate(bassurpath)
if (allocated(ple3n1path)) deallocate(ple3n1path)
end subroutine





!!!--------- Construct global array CPpos_tmp and CPtype_tmp, which contains real CPs and replicated boundary CPs, as well as CP_tmp_idx
!numcp_tmp is returned variable containing number of all CPs including the replicated boundary ones
subroutine construct_CPtmp_withbound(numcp_tmp)
use defvar
use topo
implicit real*8 (a-h,o-z)
real*8 Cart(3),fract(3),ftest(3)
integer numcp_tmp

call getcellabc(asize,bsize,csize,alpha,beta,gamma)
devthres=1D-2 !If distance to cell boundary is less than 0.01 Bohr, then it will be regarded as boundary CP
devthres2=devthres**2
fdeva=devthres/asize
fdevb=devthres/bsize
fdevc=devthres/csize

do iCP=1,numcp
    Cart(1:3)=CPpos(1:3,iCP)
    call Cart2fract(Cart,fract)
    xdiff=asize*abs(fract(1)-nint(fract(1)))
    ydiff=bsize*abs(fract(2)-nint(fract(2)))
    zdiff=csize*abs(fract(3)-nint(fract(3)))
    if (xdiff<devthres.or.ydiff<devthres.or.zdiff<devthres) then !Existing boundary CP is at least very close to one of cell walls
        !Try to replicate the existing boundary CP to all possible neighbouring mirror sites
        do ix=-1,1
            do iy=-1,1
                do iz=-1,1
                    if (ix==0.and.iy==0.and.iz==0) cycle
                    ftest(1)=fract(1)+ix
                    ftest(2)=fract(2)+iy
                    ftest(3)=fract(3)+iz
                    !Mirror boundary must be within or quasi within current cell
                    if (ftest(1)>-fdeva.and.ftest(1)<1+fdeva.and.ftest(2)>-fdevb.and.ftest(2)<1+fdevb.and.ftest(3)>-fdevc.and.ftest(3)<1+fdevc) then
                        call fract2Cart(ftest,Cart)
                        !Check if the mirror boundary CP is too close to existing CPs
                        iadd=1
                        do jCP=1,numcp
                            if (jCP==iCP) cycle
                            dist2=sum((Cart(:)-CPpos(:,jCP))**2)
                            if (dist2<devthres2) then !Too close, skip it
                                iadd=0
                                exit
                            end if
                        end do
                        if (iadd==1) then
                            numcp_tmp=numcp_tmp+1
                            CP_tmp_idx(numcp_tmp)=iCP
                            CPtype_tmp(numcp_tmp)=CPtype(iCP)
                            CPpos_tmp(:,numcp_tmp)=Cart(:)
                        end if
                    end if
                end do
            end do
        end do
    end if
end do
end subroutine





!!!--------- Construct global array pathnumpt_tmp and topopath_tmp, which contains real paths and replicated boundary paths, as well as path_tmp_idx
!numpath_tmp is returned variable containing number of all CPs including the replicated boundary ones
subroutine construct_pathtmp_withbound(numpath_tmp)
use defvar
use topo
implicit real*8 (a-h,o-z)
real*8 Cart(3),fract(3),fract2(3),fracttmp(3),ftest(3),ftest2(3)
integer numpath_tmp

call getcellabc(asize,bsize,csize,alpha,beta,gamma)
devthres=1D-2 !If distance to cell boundary is less than 0.01 Bohr, then it will be regarded as boundary CP
fdeva=devthres/asize
fdevb=devthres/bsize
fdevc=devthres/csize

do ipath=1,numpath
    call Cart2fract(topopath(1:3,1,ipath),fract)
    xdiff1=asize*abs(fract(1)-nint(fract(1)))
    ydiff1=bsize*abs(fract(2)-nint(fract(2)))
    zdiff1=csize*abs(fract(3)-nint(fract(3)))
    call Cart2fract(topopath(1:3,pathnumpt(ipath),ipath),fract2)
    xdiff2=asize*abs(fract2(1)-nint(fract2(1)))
    ydiff2=bsize*abs(fract2(2)-nint(fract2(2)))
    zdiff2=csize*abs(fract2(3)-nint(fract2(3)))
    !Both first and last points of the path are very close to one of cell walls, so this is a boundary path
    if ((xdiff1<devthres.or.ydiff1<devthres.or.zdiff1<devthres).and.(xdiff2<devthres.or.ydiff2<devthres.or.zdiff2<devthres)) then
        !Try to replicate the existing boundary CP to all possible neighbouring mirror sites
        do ix=-1,1
            do iy=-1,1
                do iz=-1,1
                    if (ix==0.and.iy==0.and.iz==0) cycle
                    ftest(1)=fract(1)+ix
                    ftest(2)=fract(2)+iy
                    ftest(3)=fract(3)+iz
                    ftest2(1)=fract2(1)+ix
                    ftest2(2)=fract2(2)+iy
                    ftest2(3)=fract2(3)+iz
                    !Mirror boundary must be within or quasi within current cell
                    if (ftest(1)>-fdeva.and.ftest(1)<1+fdeva.and.ftest(2)>-fdevb.and.ftest(2)<1+fdevb.and.ftest(3)>-fdevc.and.ftest(3)<1+fdevc) then
                        !Check if the mirror boundary path is too close to existing path
                        iadd=1
                        !do jCP=1,numcp !To be modified
                        !    if (jCP==iCP) cycle
                        !    dist2=sum((Cart(:)-CPpos(:,jCP))**2)
                        !    if (dist2<devthres2) then !Too close, skip it
                        !        iadd=0
                        !        exit
                        !    end if
                        !end do
                        if (iadd==1) then
                            numpath_tmp=numpath_tmp+1
                            path_tmp_idx(numpath_tmp)=ipath
                            pathnumpt_tmp(numpath_tmp)=pathnumpt(ipath)
                            do ipt=1,pathnumpt(ipath)
                                call Cart2fract(topopath(1:3,ipt,ipath),fracttmp)
                                fracttmp(1)=fracttmp(1)+ix
                                fracttmp(2)=fracttmp(2)+iy
                                fracttmp(3)=fracttmp(3)+iz
                                call fract2Cart(fracttmp,topopath_tmp(:,ipt,numpath_tmp))
                            end do
                        end if
                    end if
                end do
            end do
        end do
    end if
end do
end subroutine
