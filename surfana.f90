!!!============== Molecular surface analysis, output a variety of properties
subroutine surfana
use defvar
use surfvertex
use util
use functions
use GUI
implicit real*8 (a-h,o-z)
integer*2,allocatable :: corpos(:,:,:) !corner position
logical,allocatable :: ifbndcub(:,:,:) !if true, means this is a boundary cub
integer,allocatable :: mergerelat(:),HirBecatm(:)
integer tmpintarr3(3),wrtnumedr
character pdbfilename*200,filename_tmp*200,c80tmp*80,c200tmp*200,c400tmp*400,c2000tmp*2000,c10000tmp*10000,selectyn,grdfilename*200,char1tmp
character FPAfile1*200,FPAfile2*200
real*8 fragsurarea(ncenter,3),fragsuravg(ncenter,3),fragsurvar(ncenter,3) !Area, average value and variance of each atom surface. 1,2,3 corresponds to all,positive,negative part
real*8 fragsurmax(ncenter),fragsurmin(ncenter),fragsurchgsep(ncenter)
integer surfrag(ncenter),ifatmfrag(ncenter) !User-defined fragment contain which atoms; ifatmfrag(iatm)=1/0 means iatm belong / doesn't belong to user-defined fragment
integer,allocatable :: surtrifrag(:) !Each surface triangle belongs to which fragment
real*8 nucchgbackup(ncenter) !Backup nuclear charge, because which may be flushed by .chg file
integer isurftype !How to define the surface. 1=Isosurface of electron density, 2=A certain real space function, 5/6=Hirshfeld/Becke surface, 10=Isosurface of existing grid data
real*8 Pvec(ncenter)
integer,allocatable :: tmpintarr(:)

isurftype=1
surfisoval=0.001D0
imapfunc=1
ifelim=1 !If elimnate redundant surface vertices
ireadextmapval=0
grdspc=0.25D0
spcmergeratio=0.5D0
critmerge=grdspc*spcmergeratio !If the distance between two surface vertices smaller than this value, merge them
vdwmulti=1.7D0
nbisec=3
iFPA=0 !By default do not use focal-point approximation (FPA)

surfanaloop: do while(.true.)
do while(.true.)
	tetravol0=0D0
	tetravol1=0D0
	tetravol2=0D0
	tetravol3=0D0
	write(*,*)
	write(*,"(a)") " NOTE: If this module is used in your research, please cite both the Multiwfn original papers (those shown when Multiwfn boots up) and the following paper, &
	&the latter detailedly described the underlying algorithm employed in this module"
	write(*,"(a)") " Tian Lu, Feiwu Chen, Quantitative analysis of molecular surface based on improved Marching Tetrahedra algorithm, &
	&J. Mol. Graph. Model., 38, 314-323 (2012)"
	write(*,*)
	write(*,*) "     ============= Quantitative Molecular Surface Analysis ============="
	write(*,*) "-1 Return to main menu"
	write(*,*) "0 Start analysis now!"
	if (isurftype==1) write(*,"(a,f9.5)") " 1 Select the way to define surface, current: Electron density, iso:",surfisoval
	if (isurftype==2) write(*,"(a,i4,a,f9.5)") " 1 Select the way to define surface, current: Function",ifuncintp,", iso:",surfisoval
	if (isurftype==5) write(*,"(a,i5,a)") " 1 Select the way to define surface, current: Hirshfeld surface, for",nHirBecatm," atoms"
	if (isurftype==6) write(*,"(a,i5,a)") " 1 Select the way to define surface, current: Becke surface, for",nHirBecatm," atoms"
	if (isurftype==10) write(*,"(a,f9.5)") " 1 Select the way to define surface, current: Existing grid data, iso:",surfisoval
	
	if (imapfunc==-1) write(*,*) "2 Select mapped function, current: User-defined function"
	if (imapfunc==0) write(*,*) "2 Select mapped function, current: Function from external file"
	if (imapfunc==1) write(*,*) "2 Select mapped function, current: Electrostatic potential (ESP)"
	if (imapfunc==2) write(*,*) "2 Select mapped function, current: Average local ionization energy (ALIE)"
	if (imapfunc==3) write(*,"(a)") " 2 Select mapped function, current: Electrostatic potential from atomic charges"
	if (imapfunc==4) write(*,*) "2 Select mapped function, current: Local electron affinity"
	if (imapfunc==-4) write(*,*) "2 Select mapped function, current: Local electron attachment energy (LEAE)"
	if (imapfunc==5) write(*,*) "2 Select mapped function, current: Electron delocalization range function EDR(r;d)" 
	if (imapfunc==6) write(*,*) "2 Select mapped function, current: Orbital overlap distance function D(r)"	
	if (imapfunc==10) write(*,*) "2 Select mapped function, current: Pair density"
	if (imapfunc==11) then
		if (allocated(b)) then
			write(*,*) "2 Select mapped function, current: Electron density"
		else
			write(*,*) "2 Select mapped function, current: Promolecular electron density"
        end if
    end if
	if (imapfunc==12) write(*,*) "2 Select mapped function, current: Sign(lambda2)*rho"
	if (imapfunc==20) write(*,*) "2 Select mapped function, current: d_i"
	if (imapfunc==21) write(*,*) "2 Select mapped function, current: d_e"
	if (imapfunc==22) write(*,*) "2 Select mapped function, current: d_norm"
	
	if (isurftype/=10) write(*,"(a,f10.6)") " 3 Spacing of grid points for generating molecular surface:",grdspc
	write(*,*) "4 Advanced options"
	if (ireadextmapval==0) write(*,*) "5 Loading mapped function values from external file, current: No"
	if (ireadextmapval==1) write(*,*) "5 Loading mapped function values from external file, current: Vertices data"
	if (ireadextmapval==2) write(*,*) "5 Loading mapped function values from external file, current: Cubegen output"
	if (ireadextmapval==3) write(*,*) "5 Loading mapped function values from external file, current: Cube file"
    write(*,*) "6 Start analysis without considering mapped function"
    read(*,*) isel
    
	if (isel==-1) then
		exit surfanaloop
		
	else if (isel==0) then
        iskipmapfunc=0
		exit
        
	else if (isel==6) then
        iskipmapfunc=1
		exit
		
	else if (isel==1) then !Select the way to define surface
		write(*,*) "How to define the surface to be analyzed?"
		write(*,*) "1 Isosurface of electron density"
		write(*,*) "2 Isosurface of a specific real space function "
		write(*,*) "5 Hirshfeld surface (isosurface of Hirshfeld weight) of a fragment"
		write(*,*) "6 Becke surface (isosurface of Becke weight) of a fragment"
		write(*,*) "10 Isosurface of a grid data loaded from external file"
		if (allocated(cubmat)) write(*,*) "11 The isosurface of the grid data in memory"
		read(*,*) isurftypetmp
		if (isurftypetmp==1) then
			isurftype=1
		else if (isurftypetmp==2) then
			isurftype=2
			call selfunc_interface(1,ifuncintp)
		else if (isurftypetmp==5.or.isurftypetmp==6) then
			if (isurftypetmp==5) then
				isurftype=5
                grdspc=0.2D0 !Make points in finger plot denser
			else if (isurftypetmp==6) then
				isurftype=6
            end if
            do while(.true.)
				write(*,"(a)") " Input atomic indices. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will be selected"
				read(*,"(a)") c10000tmp
				call str2arr(c10000tmp,nHirBecatm)
				if (allocated(HirBecatm)) deallocate(HirBecatm)
				allocate(HirBecatm(nHirBecatm))
				call str2arr(c10000tmp,nHirBecatm,HirBecatm)
                if (nHirBecatm==ncenter) then
					write(*,"(a)") " Error: You should not define all atoms as the fragment for Hirshfeld or Becke surface analysis! &
                    &Please carefully read Multiwfn manual to correctly understand the idea of this kind of analysis. The defined fragment &
                    &should correspond to a subset of the whole system"
                    write(*,*) "Press ENTER button to re-input the atom indices"
                    read(*,*)
                else
					write(*,*) "Done! The fragment has been defined"
					exit
                end if
            end do
			imapfunc=11
		else if (isurftypetmp==10) then
			isurftype=10
			write(*,*) "Input the path of .cub or .grd file, e.g. C:\ltwd.cub"
			do while(.true.)
				read(*,"(a)") grdfilename
				inquire(file=grdfilename,exist=alive)
				if (alive) exit
				write(*,*) "Cannot find the file, input again"
			end do
			inamelen=len_trim(grdfilename)
			if (grdfilename(inamelen-2:inamelen)=="cub".or.grdfilename(inamelen-3:inamelen)=="cube") then
				call readcube(grdfilename,1,1)
			else if (grdfilename(inamelen-2:inamelen)=="grd") then
				call readgrd(grdfilename,1,1)
			end if
			if (dx/=dy.or.dy/=dz) write(*,*) "Warning: The grid in this file is not cubic! The result may be problematic"
			critmerge=min(dx,dy,dz)*spcmergeratio
		else if (isurftypetmp==11) then
			isurftype=10
			critmerge=min(dx,dy,dz)*spcmergeratio
		end if
		!Set isovalue
		if (isurftype==5.or.isurftype==6) then
			surfisoval=0.5D0 !Hirshfeld and Becke surface the isoval must be 0.5
		else
			write(*,*) "Input the isovalue for defining the isosurface, e.g. 0.5"
			if (isurftype==1) write(*,*) "Hint: Isovalue of 0.001 is commonly used to define molecular vdW surface"
			read(*,*) surfisoval
		end if
		!Set additional parameters
		if (isurftype==1) then !Only remove redundant vertices for electron density isosurface, for other types, elimination method may work poorly
			ifelim=1
			critmerge=grdspc*0.5D0
			nbisec=3
		else if (isurftype==2) then
			ifelim=0
			nbisec=3
		else if (isurftype==5.or.isurftype==6.or.isurftype==10) then
			ifelim=1
			nbisec=0 !Linear interpolation
		end if
		if (.not.allocated(b).and.isurftype/=5.and.isurftype/=6) imapfunc=-1
		
	else if (isel==2) then !Select mapped function
		if (imapfunc==3) a%charge=nucchgbackup !If .chg file was loaded previously by option 3 shown below, so this time firstly recovery actual nuclear charges
		write(*,*) "Select to real space function to be mapped on the molecular surface"
		write(*,"(a,i4)") "-1 User-defined real space function, iuserfunc=",iuserfunc
		write(*,*) "0 Function from external file"
		write(*,*) "1 Electrostatic potential (ESP)"
		write(*,*) "2 Average local ionization energy (ALIE)"
		write(*,*) "3 Electrostatic potential from atomic charges in a .chg file"
		write(*,*) "4 Local electron affinity (LEA)"
		write(*,*) "-4 Local electron attachment energy (LEAE)"
		write(*,*) "5 Electron delocalization range function EDR(r;d)"	
		write(*,*) "6 Orbital overlap length function D(r) which maximizes EDR(r;d)"	
! 		if (allocated(b)) write(*,*) "10 Pair density" !Rarely used by normal users, so comment it
		if (allocated(b)) then
			write(*,*) "11 Electron density"
			write(*,*) "12 Sign(lambda2)*rho"
		else
			write(*,*) "11 Electron density with promolecular approximation"
			write(*,*) "12 Sign(lambda2)*rho with promolecular approximation"
		end if
		if (nHirBecatm>0) then
			write(*,*) "20 d_i: Distance in Angstrom from the nearest nucleus inside the surface"
			write(*,*) "21 d_e: Distance in Angstrom from the nearest nucleus outside the surface"
			write(*,*) "22 d_norm: Normalized contact distance in Angstrom"
		end if
		read(*,*) imapfunc
		
		if (imapfunc==0) then !Determine appropriate grid spacing
			ireadextmapval=1
			grdspc=0.2D0
		else
			ireadextmapval=0
			if (imapfunc==1.or.imapfunc==11.or.imapfunc==12) grdspc=0.25D0
			if (imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4.or.imapfunc==5.or.imapfunc==6.or.imapfunc==-1.or.imapfunc==10) grdspc=0.2D0
		end if
		critmerge=grdspc*spcmergeratio
        if (imapfunc==-4) then
			surfisoval=0.004D0 !Commonly LEAE is studied on rho=0.004 a.u. isosurface
            write(*,*) "Note: Isovalue has been changed to 0.004 a.u."
		else if (imapfunc==3) then
			write(*,"(a)") " Please input the path of the .chg file which contains the atomic charges of present system. e.g. C:\sob.chg"
			do while(.true.)
				read(*,"(a)") c200tmp
				inquire(file=c200tmp,exist=alive)
				if (alive) exit
				write(*,*) "Error: Cannot find the file, input again"
			end do
			nucchgbackup=a%charge
			open(10,file=c200tmp,status="old")
			do icen=1,ncenter
				read(10,*) c80tmp,rnouse,rnouse,rnouse,a(icen)%charge
				if (trim(c80tmp)/=trim(a(icen)%name)) write(*,"(' Warning: The name of atom',i7,' in .chg file is not consistent with present system',i7)") icen
			end do
			close(10)
			write(*,*) "The atomic charges have been loaded"
		else if (imapfunc==5) then !Input length scale to evaluate EDR(r;d)
			write(*,*) "The EDR(r;d) computing code was contributed by Arshad Mehmood"
			write(*,"(a,/)") " References: J. Chem. Phys., 141, 144104 (2014); J. Chem. Theory Comput., 12, 79 (2016); Angew. Chem. Int. Ed., 56, 6878 (2017)"
			write(*,*) " Input length scale d (Bohr), e.g. 0.85"
			read(*,*) dedr
		else if (imapfunc==6) then !Input parameters to evaluate D(r)
			write(*,*) "The D(r) computing code was contributed by Arshad Mehmood"
			write(*,"(a,/)") " References: J. Chem. Theory Comput., 12, 3185 (2016); Phys. Chem. Chem. Phys., 17, 18305 (2015)"
			write(*,*) "1 Manually input total number, start and increment in EDR exponents"
			write(*,*) "2 Use default values   i.e. 20,2.50,1.50"
			read(*,*) edrmaxpara
			if (edrmaxpara==1) then  
				write(*,*) "Please input in order: exponents start increment, e.g. 20,2.5,1.5"
				write(*,*) "Note: Max. allowed exponents are 50 and min. allowed increment is 1.01"
				read (*,*) nedr,edrastart,edrainc
				if (nedr<1) then
					write(*,*) "Error: Bad number of EDR exponents. Should be between 1 to 50"
					write(*,*) "Press ENTER button to exit"
					read(*,*)
					stop
				else if (nedr>50) then
					write(*,*) "Error: Bad number of EDR exponents. Should be between 1 to 50"
					write(*,*) "Press ENTER button to exit"
					read(*,*)
					stop
				end if
				if (edrainc<1.01d0) then
					write(*,*) "Error: Bad increment in EDR exponents. Should not be less than 1.01"
					write(*,*) "Press ENTER button to exit"
					read(*,*)
					stop
				end if
			else if (edrmaxpara==2) then
				nedr=20
				edrastart=2.5d0
				edrainc=1.5d0
			end if
			write(*,*) "The following EDR exponents will be used in calculation:"
			wrtstart=edrastart
			do wrtnumedr=1,nedr
				wrtexpo(wrtnumedr)=wrtstart
				wrtstart=wrtstart/edrainc
				write(*,"(E13.5)") wrtexpo(wrtnumedr) 
			end do
		else if (imapfunc==20.or.imapfunc==21) then
			write(*,*) "NOTE: ALL VALUES OF THIS FUNCTION SHOWN IN LATER STAGE WILL BE BOHR!"
		end if
		
	else if (isel==3) then
		write(*,*) "Input a value (in Bohr), e.g. 0.2"
		if (imapfunc==0.or.imapfunc==1.or.imapfunc==20.or.imapfunc==21.or.imapfunc==22) write(*,*) "Note: In general 0.25 is enough. For higher accuracy, 0.15~0.20 is recommended"
		if (imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4.or.imapfunc==-1) write(*,*) "Note: In general 0.20 is enough. For higher accuracy, 0.13~0.17 is recommended"
		read(*,*) grdspc
		critmerge=grdspc*spcmergeratio
		
	else if (isel==4) then
		do while(.true.)
			write(*,*)
			write(*,*) "0 Return to upper level menu"
			write(*,"(a,f7.4)") " 1 Set ratio of vdW radii used to extend spatial region for grids:",vdwmulti
			if (ifelim==0) write(*,*) "2 Toggle if eliminating redundant vertices: No"
			if (ifelim==1) write(*,"(a,f6.3,a)") " 2 Toggle if eliminating redundant vertices: Yes, criterion is",critmerge," Bohr"
			write(*,"(' 3 Number of bisections before linear interpolation, current:',i5)") nbisec
			if (iFPA==0) write(*,*) "4 Toggle using focal-point approximation to evaluate ESP, current: No"
			if (iFPA==1) write(*,*) "4 Toggle using focal-point approximation to evaluate ESP, current: Yes"
			read(*,*) isel2
			
			if (isel2==0) then
				exit
			else if (isel2==1) then
				write(*,*) "Input ratio of vdW radii, e.g. 1.5"
				if (isurftype==1) write(*,"(a)") " Note: 1.7 is enough for the case of isovalue=0.001, for lower isovalue, a larger value is needed"
				read(*,*) vdwmulti
			else if (isel2==2) then
				if (ifelim==1) then
					ifelim=0
				else if (ifelim==0) then
					write(*,"(a)") " Input a value (in Bohr). If the distance between any two vertices is smaller than this value, one of them will be merged to the other."
					write(*,"(a)") " Hint: 0.5 times grid spacing is recommended in general"
					read(*,*) critmerge
					ifelim=1
				end if
			else if (isel2==3) then
				write(*,*) "Perform how many times biections before linear interpolation?"
				write(*,*) "Note: Input 0 means do interpolation directly"
				read(*,*) nbisec
            else if (isel2==4) then
				if (iFPA==0) then
					write(*,*) "Input the wavefunction file for ""low method / small basis-set"""
                    write(*,*) "e.g. /sob/MP2_def2_TZVP.mwfn"
                    write(*,"(a)") " Note: Any format containing recognizable wavefunction information can be used. The coordinate must be identical to present system"
                    do while(.true.)
						read(*,"(a)") FPAfile1
						inquire(file=FPAfile1,exist=alive)
						if (alive) exit
						write(*,*) "Cannot find the file, input again!"
					end do
					write(*,*) "Input the wavefunction file for ""low method / large basis-set"""
                    write(*,*) "e.g. /sob/MP2_def2_QZVP.mwfn"
                    do while(.true.)
						read(*,"(a)") FPAfile2
						inquire(file=FPAfile2,exist=alive)
						if (alive) exit
						write(*,*) "Cannot find the file, input again!"
					end do
                    write(*,*) "OK, focal-point approximation will be used in the calculation"
                    iFPA=1
                else
					iFPA=0
                end if
			end if
		end do
		
	else if (isel==5) then
		write(*,*) "0 Do not load mapped function but directly calculate by Multiwfn"
		write(*,*) "1 Load mapped function at all surface vertices from plain text file"
		write(*,*) "2 Similar to 1, but specific for the case of using cubegen utility of Gaussian"
		write(*,*) "3 Mapped function will be interpolated from an external cube file"
		read(*,*) ireadextmapval
	end if
	
end do


!======== Start calculation ========!
!======== Start calculation ========!
!======== Start calculation ========!

if (imapfunc/=0.and.imapfunc/=4.and.imapfunc/=20.and.imapfunc/=21.and.imapfunc/=22) call delvirorb(1) ! Delete high-lying virtual orbitals to speed up calculation
if (allocated(b)) then
	if (ifPBC==0) then
		call gen_GTFuniq(0) !Generate unique GTFs, for faster evaluation in orbderv
	else
		call gen_neigh_GTF !Generate neighbouring GTFs list at reduced grids, for faster evaluation
	end if
end if

call walltime(iclktime1)
if (isurftype==1.or.isurftype==2.or.isurftype==5.or.isurftype==6) then !Calculate grid data for determining isosurface
	if (isurftype==5.or.isurftype==6) then !Hirshfeld/Becke surface, only used selected atoms to define the box
		orgx=minval( a(HirBecatm(:))%x-vdwmulti*vdwr(a(HirBecatm(:))%index) )
		orgy=minval( a(HirBecatm(:))%y-vdwmulti*vdwr(a(HirBecatm(:))%index) )
		orgz=minval( a(HirBecatm(:))%z-vdwmulti*vdwr(a(HirBecatm(:))%index) )
		endx=maxval( a(HirBecatm(:))%x+vdwmulti*vdwr(a(HirBecatm(:))%index) )
		endy=maxval( a(HirBecatm(:))%y+vdwmulti*vdwr(a(HirBecatm(:))%index) )
		endz=maxval( a(HirBecatm(:))%z+vdwmulti*vdwr(a(HirBecatm(:))%index) )
	else
		orgx=minval( a(:)%x-vdwmulti*vdwr(a(:)%index) )
		orgy=minval( a(:)%y-vdwmulti*vdwr(a(:)%index) )
		orgz=minval( a(:)%z-vdwmulti*vdwr(a(:)%index) )
		endx=maxval( a(:)%x+vdwmulti*vdwr(a(:)%index) )
		endy=maxval( a(:)%y+vdwmulti*vdwr(a(:)%index) )
		endz=maxval( a(:)%z+vdwmulti*vdwr(a(:)%index) )
	end if
	
	write(*,"(' Spatial range of grid data:')")
	write(*,"(' X is from',f10.4,'  to',f10.4,' Bohr')") orgx,endx
	write(*,"(' Y is from',f10.4,'  to',f10.4,' Bohr')") orgy,endy
	write(*,"(' Z is from',f10.4,'  to',f10.4,' Bohr')") orgz,endz
	xlength=endx-orgx
	ylength=endy-orgy
	zlength=endz-orgz
	dx=grdspc;gridv1=0;gridv1(1)=dx
	dy=grdspc;gridv2=0;gridv2(2)=dy
	dz=grdspc;gridv3=0;gridv3(3)=dz
	nx=nint(xlength/dx)+1
	ny=nint(ylength/dy)+1
	nz=nint(zlength/dz)+1
	write(*,"(' The number of point in x,y,z:',3i6,'  Total:',i10)") nx,ny,nz,nx*ny*nz
	write(*,*)
	
	if (isurftype==1.or.isurftype==2) then
		if (isurftype==1) then
            write(*,*) "Calculating grid data of electron density..."
            if (.not.allocated(b)) write(*,*) "Note: Promolecular approximation is used"
		else if (isurftype==2) then
            write(*,*) "Calculating grid data of the real space function..."
        end if
		if (allocated(cubmat)) deallocate(cubmat)
		allocate(cubmat(nx,ny,nz))
		if (isurftype==1) then
			call savecubmat(1,0,1)
		else if (isurftype==2) then
			call savecubmat(ifuncintp,0,iorbsel)
		end if
		!Calculate MICC(molecular intrinsic characteristic contour), however it seems that the result obtained in this manner is not in line with the one in original paper
		! iuserfunc=10
		! call savecubmat(100,1,1)
		
	else if (isurftype==5) then !Hirshfeld analysis, calculate weighting distribution of a specific set of atom
		if (allocated(cubmat)) deallocate(cubmat)
		if (allocated(cubmattmp)) deallocate(cubmattmp)
		allocate(cubmat(nx,ny,nz),cubmattmp(nx,ny,nz))
		cubmat=0D0
		cubmattmp=0D0
		!I found ihirshmode=1 is always a satisfactory choice, so I decide not bother users to select
		!write(*,*) "Hirshfeld analysis requests atomic densities, please select how to obtain them"
		!write(*,*) "1 Use build-in sphericalized atomic densities in free-states (recommended)"
		!write(*,"(a)") " 2 Provide wavefunction file of involved elements by yourself or invoke Gaussian to automatically calculate them"
		!read(*,*) ihirshmode
		!if (ihirshmode==1) then
			do iatm=1,ncenter
				!I attempted to use truncation, which reduce cost significantly for >1000 atoms system, unfortunately this will cost Hirshfeld surface occur at molecular open boundary...
				!xtmp=a(iatm)%x;ytmp=a(iatm)%y;ztmp=a(iatm)%z
    !            crit2=(vdwr(a(iatm)%index)*3)**2
                iadd=0
                if (any(HirBecatm==iatm)) iadd=1
				!$OMP PARALLEL DO SHARED(cubmat,cubmattmp) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp,dist2) schedule(dynamic) NUM_THREADS(nthreads)
				do k=1,nz
					do j=1,ny
						do i=1,nx
							call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                            !dist2=(xtmp-tmpx)**2+(ytmp-tmpy)**2+(ztmp-tmpz)**2
                            !if (dist2>crit2) cycle
							denstmp=calcatmdens(iatm,tmpx,tmpy,tmpz,-1) !Use STO fitted atomic density, can cover rather long distance to avoid zero denominator of Hirshfeld weight
							cubmattmp(i,j,k)=cubmattmp(i,j,k)+denstmp
							if (iadd==1) cubmat(i,j,k)=cubmat(i,j,k)+denstmp !Density of specified fragment
						end do
					end do
				end do
				!$OMP END PARALLEL DO
				call showprog(iatm,ncenter)
			end do
		!else if (ihirshmode==2) then
		!	call setpromol
		!	do iatm=1,ncenter_org
		!		call dealloall(0)
		!		call readwfn(custommapname(iatm),1)
		!		!$OMP PARALLEL DO SHARED(cubmat,cubmattmp) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp) schedule(dynamic) NUM_THREADS(nthreads)
		!		do k=1,nz
		!			do j=1,ny
		!				do i=1,nx
		!					call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
		!					denstmp=fdens(tmpx,tmpy,tmpz)
		!					cubmattmp(i,j,k)=cubmattmp(i,j,k)+denstmp
		!					if (any(HirBecatm==iatm)) cubmat(i,j,k)=cubmat(i,j,k)+denstmp !Density of specified fragment
		!				end do
		!			end do
		!		end do
		!		!$OMP END PARALLEL DO
		!		call showprog(iatm,ncenter_org)
		!	end do
		!	call dealloall(0)
		!	write(*,"(' Reloading ',a)") trim(firstfilename)
		!	call readinfile(firstfilename,1) !Retrieve to first loaded file(whole molecule)
		!end if
		do k=1,nz
			do j=1,ny
				do i=1,nx
					if (cubmattmp(i,j,k)==0) then !At very distant region, promolecular density can be zero
						cubmat(i,j,k)=0
                    else
						cubmat(i,j,k)=cubmat(i,j,k)/cubmattmp(i,j,k)
                    end if
                end do
            end do
        end do

	else if (isurftype==6) then !Becke analysis, calculate weighting distribution of a specific set of atom
		if (allocated(cubmat)) deallocate(cubmat)
		allocate(cubmat(nx,ny,nz))
		cubmat=0D0
		ifinish=0
        call showprog(0,nz)
		!We calculate Becke weight for all atoms, but only summing up the value of we selected atoms to cubmat
		!$OMP PARALLEL DO SHARED(cubmat) PRIVATE(i,j,k,rnowx,rnowy,rnowz,Pvec) schedule(dynamic) NUM_THREADS(nthreads)
		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,rnowx,rnowy,rnowz)
                    !Calculate Becke weight of all atoms (Pvec) at current point
					call BeckePvec(rnowx,rnowy,rnowz,Pvec,covr,3)
					cubmat(i,j,k)=cubmat(i,j,k)+sum(Pvec(HirBecatm(:)))
				end do
			end do
			!$OMP CRITICAL
			ifinish=ifinish+1
			call showprog(ifinish,nz)
			!$OMP end CRITICAL
		end do
		!$OMP END PARALLEL DO
	end if
end if

allocate(corpos(nx,ny,nz),ifbndcub(nx-1,ny-1,nz-1))
write(*,*)

!Initialize
numcubx=nx-1
numcuby=ny-1
numcubz=nz-1
ifbndcub=.false.
nsurvtx=0
nsurtri=0
tetravol1=0D0 !Volume in generated surface tetrahedron of type 1,2,3, will be accumulated in polygonization process
tetravol2=0D0
tetravol3=0D0
allocate(abs2suridx(nx,ny,nz))
abs2suridx=0
surlocminidx=0
surlocmaxidx=0

!Which corners(nx*ny*nz in total) is in internal (=0) or external (=1) of specified isosurface
do ix=1,nx
	do iy=1,ny
		do iz=1,nz
			if (cubmat(ix,iy,iz)>=surfisoval) then !internal corner
				corpos(ix,iy,iz)=0
			else !external corner
				corpos(ix,iy,iz)=1
			end if
		end do
	end do
end do

nintcub=0
nextcub=0
nbndcub=0
icorind=0
do ix=1,numcubx
	do iy=1,numcuby
		do iz=1,numcubz
			icubtest=corpos(ix,iy,iz)+corpos(ix+1,iy,iz)+corpos(ix,iy+1,iz)+corpos(ix,iy,iz+1)+corpos(ix+1,iy+1,iz)+corpos(ix,iy+1,iz+1)+corpos(ix+1,iy,iz+1)+corpos(ix+1,iy+1,iz+1)
			if (icubtest==0) then !internal cube
				nintcub=nintcub+1
			else if (icubtest==8) then !external cube
				nextcub=nextcub+1
			else !boundary cube
				ifbndcub(ix,iy,iz)=.true.
				nbndcub=nbndcub+1
				!Give each corner of boundary cube a unique index
				if (abs2suridx(ix,iy,iz)==0) then
					icorind=icorind+1
					abs2suridx(ix,iy,iz)=icorind
				end if
				if (abs2suridx(ix+1,iy,iz)==0) then
					icorind=icorind+1
					abs2suridx(ix+1,iy,iz)=icorind
				end if
				if (abs2suridx(ix,iy+1,iz)==0) then
					icorind=icorind+1
					abs2suridx(ix,iy+1,iz)=icorind
				end if
				if (abs2suridx(ix,iy,iz+1)==0) then
					icorind=icorind+1
					abs2suridx(ix,iy,iz+1)=icorind
				end if
				if (abs2suridx(ix+1,iy+1,iz)==0) then
					icorind=icorind+1
					abs2suridx(ix+1,iy+1,iz)=icorind
				end if
				if (abs2suridx(ix,iy+1,iz+1)==0) then
					icorind=icorind+1
					abs2suridx(ix,iy+1,iz+1)=icorind
				end if
				if (abs2suridx(ix+1,iy,iz+1)==0) then
					icorind=icorind+1
					abs2suridx(ix+1,iy,iz+1)=icorind
				end if
				if (abs2suridx(ix+1,iy+1,iz+1)==0) then
					icorind=icorind+1
					abs2suridx(ix+1,iy+1,iz+1)=icorind
				end if
			end if
		end do
	end do
end do
numsurfcubcor=count(abs2suridx/=0)
deallocate(corpos)
write(*,"(' The number of boundary cubes:        ',i10)") nbndcub  !debug info
write(*,"(' The number of corners of boundary cubes:',i10)") numsurfcubcor  !debug info

nassumsurfvtx=7*nbndcub !I assume can generate up to 7*nbndcub surface vertices. In principle, each cube can generate up to 19 vertices, since there are 19 edges in main-axis decomposition
nassumfacet=12*nbndcub !each boundary cube can generate up to 2*6=12 facets in principle
allocate(survtx(nassumsurfvtx)) !, so in principle should be 19*nbndcub... 
allocate(surtriang(nassumfacet))
allocate(surfcor2vtx(numsurfcubcor,14),surcor2vtxpos(numsurfcubcor)) !14 is max number that each corner can link other corner for main-axis decomposition of cube
surfcor2vtx(:,:)%athcor=0
surcor2vtxpos=0
allocate(vtxconn(nassumsurfvtx,100),vtxconnpos(nassumsurfvtx)) !I assume that each surface vertex can connect up to 100 other surface vertices
vtxconn=0
vtxconnpos=0

write(*,*) "Generating isosurface by Marching Tetrahedra algorithm, please wait..."
!DO NOT use parallel method in this step, otherwise the numbering will be disordered, and we can't read external mapped function data, since indices cannot be consistent
do ix=1,numcubx
	do iy=1,numcuby
		do iz=1,numcubz
			if (ifbndcub(ix,iy,iz)) then !Numbering of cube corner is identical to figure 3 of WFA original paper
				call marchtetra(ix,iy,iz)
! 				call marchcube(ix,iy,iz)
			end if
		end do
	end do
end do
deallocate(ifbndcub,abs2suridx,surfcor2vtx,surcor2vtxpos) !Discard arrays used in polygonization

nsuredge=sum(vtxconnpos(1:nsurvtx))/2
write(*,"(' The number of surface vertices (V): ',i10)") nsurvtx
write(*,"(' The number of surface edges (E):    ',i10)") nsuredge
write(*,"(' The number of triangular facets (F):',i10)") nsurtri
iVFEtest=nsurvtx+nsurtri-nsuredge-2
if (imapfunc/=20.and.imapfunc/=21.and.imapfunc/=22) then !The surface of Hirshfeld/Becke is usually open and thus must break this rule, so don't check this
	write(*,"(' V+F-E-2=',i10)") iVFEtest  !debug info
	if (iVFEtest/=0) write(*,"(' Warning: V+F-E-2=0 is violated! Probably grid spacing is too large or the isosurface is not closed')")
end if
call calc_dvol(volcub)
totintvol=nintcub*volcub !Volume of internal cubes
totvol=tetravol0+tetravol1+tetravol2+tetravol3+totintvol !Volume of boundary tetrahedra
! write(*,"('Vol. of type-0,1,2,3:',4f11.7,' Ang.^3')") tetravol0*b2a**3,tetravol1*b2a**3,tetravol2*b2a**3,tetravol3*b2a**3
write(*,"(' Volume enclosed by the isosurface:',f12.5,' Bohr^3  (',f10.5,' Angstrom^3)')") totvol,totvol*b2a**3
 write(*,*) "Among all surface vertices:"
 write(*,"(' Min-X:',f12.4,'  Max-X:',f10.4,' Angstrom')") minval(survtx(1:nsurvtx)%x)*b2a,maxval(survtx(1:nsurvtx)%x)*b2a
 write(*,"(' Min-Y:',f12.4,'  Max-Y:',f10.4,' Angstrom')") minval(survtx(1:nsurvtx)%y)*b2a,maxval(survtx(1:nsurvtx)%y)*b2a
 write(*,"(' Min-Z:',f12.4,'  Max-Z:',f10.4,' Angstrom')") minval(survtx(1:nsurvtx)%z)*b2a,maxval(survtx(1:nsurvtx)%z)*b2a
write(*,*)

!Merge redundant vertices
allocate(elimvtx(nsurvtx),elimtri(nint(nsurtri*1.1D0)),mergerelat(nsurvtx)) !The reason for *1.1D0, is because during elimination of vertices with three connections, we need add new triangle
elimvtx=0 !If it is 1, means this surface vertex has been eliminated
elimtri=0
forall(i=1:nsurvtx) mergerelat(i)=i !If mergerelat(i)=j, means surface vertex i has been absorbed into vertex j
if (ifelim==1) then
	write(*,*) "Eliminating redundant surface vertices..."
	do ivt1=1,nsurvtx
		if (elimvtx(ivt1)==1) cycle
		vt1x=survtx(ivt1)%x
		vt1y=survtx(ivt1)%y
		vt1z=survtx(ivt1)%z
		do j=1,vtxconnpos(ivt1)
			ivt2=vtxconn(ivt1,j)
			if (elimvtx(ivt2)==1) cycle
			vt2x=survtx(ivt2)%x
			vt2y=survtx(ivt2)%y
			vt2z=survtx(ivt2)%z
			dist2=(vt1x-vt2x)**2+(vt1y-vt2y)**2+(vt1z-vt2z)**2
			if (dist2<critmerge**2) then !merge vt2 to vt1
				!Scan and update connectivity of neighbour of vt2 (except vt1), also add new connections to vt1 from vt2
				do k=1,vtxconnpos(ivt2)
					inei=vtxconn(ivt2,k)
					if (elimvtx(inei)==1) cycle
					if (inei==ivt1) cycle
					ihasvt1=0
					do itmp=1,vtxconnpos(inei)
						ineinei=vtxconn(inei,itmp)
						if (elimvtx(ineinei)==1) cycle
						if (ineinei==ivt1) ihasvt1=1
					end do
					if (ihasvt1==0) then
						vtxconnpos(inei)=vtxconnpos(inei)+1 !if this neighbour vertex of vt2 hasn't connected to vt1, add the connection
						vtxconn(inei,vtxconnpos(inei))=ivt1
						vtxconnpos(ivt1)=vtxconnpos(ivt1)+1
						vtxconn(ivt1,vtxconnpos(ivt1))=inei
					end if
				end do
				mergerelat(ivt2)=ivt1
				elimvtx(ivt2)=1 !eliminate ivt2
				!Change coordinate of ivtx1 to the average of ivtx1 and ivtx2
				survtx(ivt1)%x=(vt1x+vt2x)/2D0
				survtx(ivt1)%y=(vt1y+vt2y)/2D0
				survtx(ivt1)%z=(vt1z+vt2z)/2D0
				vt1x=survtx(ivt1)%x
				vt1y=survtx(ivt1)%y
				vt1z=survtx(ivt1)%z
			end if
		end do
	end do
! 	do i=1,nsurvtx
! 		write(1,*) i,mergerelat(i),elimvtx(i)
! 	end do
	
	!Update vertices of facets
	do itri=1,nsurtri
		do itmp=1,3
			do while(.true.) !Update vertex indices of each facet according to mergence relation map
				if (elimvtx(surtriang(itri)%idx(itmp))==1) then
					surtriang(itri)%idx(itmp)=mergerelat(surtriang(itri)%idx(itmp))
				else
					exit
				end if
			end do
		end do
outdo:	do itmp=1,3 !If a facet has two identical vertices, then delete this facet
			do jtmp=itmp+1,3
				if (surtriang(itri)%idx(itmp)==surtriang(itri)%idx(jtmp)) then
					elimtri(itri)=1
					exit outdo
				end if
			end do
		end do outdo
	end do
! 	write(*,"('Eliminated surface vertices in merge step:',i10)") count(elimvtx==1)
! 	write(*,"('Eliminated surface facets in merge step:',i10)") count(elimtri==1)
	
	!Eliminate surface vertices with two neighbours
	nelimtwoconvtx=0
	do ivt1=1,nsurvtx
		if (elimvtx(ivt1)==1) cycle
		nlink=0
		do j=1,vtxconnpos(ivt1)
			ivt2=vtxconn(ivt1,j)
			if (elimvtx(ivt2)==0) nlink=nlink+1
			if (nlink>2) exit
		end do
		if (nlink==2) then
			nelimtwoconvtx=nelimtwoconvtx+1
			elimvtx(ivt1)=1
		end if
	end do
	!Eliminate surface vertices with three neighbours
	nelimthreeconvtx=0
	do ivt1=1,nsurvtx
		if (elimvtx(ivt1)==1) cycle
		nlink=0
		do j=1,vtxconnpos(ivt1)
			ivt2=vtxconn(ivt1,j)
			if (elimvtx(ivt2)==0) then
				nlink=nlink+1
				if (nlink>3) exit
				tmpintarr3(nlink)=ivt2
			end if
		end do
		if (nlink==3) then
			nelimthreeconvtx=nelimthreeconvtx+1
			elimvtx(ivt1)=1
			nsurtri=nsurtri+1
			surtriang(nsurtri)%idx(:)=tmpintarr3(:)
			elimtri(nsurtri)=0
		end if
	end do
	!Remove the triangles which has at least one vertex has been discarded
! 	iz=0
	do itri=1,nsurtri
		if (elimtri(itri)==1) cycle
		do itmp=1,3
			if (elimvtx(surtriang(itri)%idx(itmp))==1) then
				elimtri(itri)=1
! 				iz=iz+1
! 				write(*,"(5i10)") iz,surtriang(itri)%idx(itmp),surtriang(itri)%idx(:)
				exit
			end if
		end do
	end do
! 	write(*,"('Eliminated surface vertices with two connections:',i10)") nelimtwoconvtx
! 	write(*,"('Eliminated surface vertices with three connections:',i10)") nelimthreeconvtx
	
	!Elimination has finished, count how many surface vertices, edges and facets currently
	ncurrtri=count(elimtri(1:nsurtri)==0)
	ncurrvtx=count(elimvtx(1:nsurvtx)==0)
	ncurredge=0
	do ivtx=1,nsurvtx
		if (elimvtx(ivtx)==1) cycle
		do itmp=1,vtxconnpos(ivtx)
			if (elimvtx(vtxconn(ivtx,itmp))==1) cycle
			ncurredge=ncurredge+1
		end do
	end do
	ncurredge=ncurredge/2
	write(*,"(' After elimination, V=',i9,',  E=',i9,',  F=',i9,',  V+F-E-2=',i6)") ncurrvtx,ncurredge,ncurrtri,ncurrvtx+ncurrtri-ncurredge-2
! 	write(*,*) maxval(vtxconnpos(1:nsurvtx))
else
	ncurrtri=nsurtri !If not do vertices elimnation
	ncurrvtx=nsurvtx
	ncurredge=nsuredge
end if

!Check deviation of surface vertices from specified isovalue
! valtotdev=0D0
! valRMSD=0D0
! valmaxdev=0D0
! num=0
! do i=1,nsurvtx
! 	if (elimvtx(i)==1) cycle
! 	num=num+1
! 	devtmp=abs(fdens(survtx(i)%x,survtx(i)%y,survtx(i)%z)-surfisoval)
! 	valtotdev=valtotdev+devtmp
! 	valRMSD=valRMSD+devtmp**2
! 	if (devtmp>valmaxdev) valmaxdev=devtmp
! end do
! valtotdev=valtotdev/num
! valRMSD=dsqrt(valRMSD/num)
! write(*,"('Average deviation from specified isovalue:',f20.16)") valtotdev
! write(*,"('RMSD from specified isovalue:',f20.16)") valRMSD
! write(*,"('Maximum deviation from specified isovalue:',f20.16)") valmaxdev
! read(*,*)

!Calculate isosurface area and geometry center of each facet
surfareaall=0D0
triareamin=100D0
triareamax=0D0
do icyc=1,nsurtri
	if (elimtri(icyc)==1) cycle
	idx1=surtriang(icyc)%idx(1)
	idx2=surtriang(icyc)%idx(2)
	idx3=surtriang(icyc)%idx(3)
	xidx1=survtx(idx1)%x
	yidx1=survtx(idx1)%y
	zidx1=survtx(idx1)%z
	xidx2=survtx(idx2)%x
	yidx2=survtx(idx2)%y
	zidx2=survtx(idx2)%z
	xidx3=survtx(idx3)%x
	yidx3=survtx(idx3)%y
	zidx3=survtx(idx3)%z
	surtriang(icyc)%area=gettriangarea( xidx1,yidx1,zidx1,xidx2,yidx2,zidx2,xidx3,yidx3,zidx3 )
	surfareaall=surfareaall+surtriang(icyc)%area
	if (surtriang(icyc)%area>triareamax) triareamax=surtriang(icyc)%area
	if (surtriang(icyc)%area<triareamin) triareamin=surtriang(icyc)%area
end do

! write(*,"('Minimum and maximum facet area:',E18.10,f15.10,' Angstrom^2')") triareamin*b2a*b2a,triareamax*b2a*b2a  !debug info
write(*,"(' Isosurface area:',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareaall,surfareaall*b2a*b2a
!https://en.wikipedia.org/wiki/Sphericity
write(*,"(' Sphericity:',f8.4)")  pi**(1D0/3D0)*(6*totvol)**(2D0/3D0)/surfareaall

if (iskipmapfunc==1) goto 10 !Skip dealing with mapped function

!Calculate mapped function at each surface vertex, save to survtx%value
if (ireadextmapval==0) then !Directly calculate
	write(*,*)
	if (imapfunc==-1) then
		write(*,*) "Calculating user-defined real space function at surface vertices..."
	else if (imapfunc==1) then
		write(*,*) "Calculating electrostatic potential at surface vertices, please wait patiently"
	else if (imapfunc==2) then
		write(*,*) "Calculating average local ionization energy at surface vertices..."
	else if (imapfunc==3) then
		write(*,*) "Calculating electrostatic potential based on atomic charges at surface vertices"
	else if (imapfunc==4) then
		write(*,*) "Calculating local electron affinity at surface vertices..."
	else if (imapfunc==-4) then
		write(*,*) "Calculating local electron attachment energy at surface vertices..."
	else if (imapfunc==5) then		
		write(*,*) "Calculating EDR(r;d) at surface vertices..."
	else if (imapfunc==6) then		
		write(*,*) "Calculating orbital overlap length D(r) at surface vertices..."
	else
		write(*,*) "Calculating mapped function value at surface vertices..."
	end if
	
	call walltime(iwalltime1)
	
	!If the mapped function is ESP, explore the possibility of using cubegen to accelerate calculation
	alive=.false.
	if (cubegenpath/=" ".and.ifiletype==1.and.imapfunc==1) then
		inquire(file=cubegenpath,exist=alive)
		if (.not.alive.and.cubegenpath/="none") then
			write(*,"(a)") " Note: Albeit current file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been defined, &
			&the cubegen cannot be found, therefore electrostatic potential will still be calculated using internal code of Multiwfn"
		end if
	end if
	if (alive.and.ifiletype==1.and.imapfunc==1) then !Use cubegen to calculate ESP
		write(*,"(a)") " Since the input file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been properly defined, &
		&now Multiwfn directly invokes cubegen to calculate electrostatic potential"
		
		!Generate cubegen input file
		open(10,file="cubegenpt.txt",status="replace")
		do icyc=1,nsurvtx
			if (elimvtx(icyc)==1) cycle
			write(10,"(3f16.8)") survtx(icyc)%x*b2a,survtx(icyc)%y*b2a,survtx(icyc)%z*b2a
		end do
		close(10)
		ncubegenthreads=1 !Parallel implementation prior to G16 is buggy, so test here
		if (index(cubegenpath,"G16")/=0.or.index(cubegenpath,"g16")/=0) ncubegenthreads=nthreads
		
		filename_tmp=filename
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
		do icyc=1,nsurvtx
			if (elimvtx(icyc)==1) cycle
			read(10,*) rnouse,rnouse,rnouse,survtx(icyc)%value
		end do
		close(10)
		
		!Delete intermediate files
		call delfile("cubegenpt.txt ESPresult.cub nouseout")
		
	else !Normal case, use Multiwfn internal code to evaluate mapped functions
		ii=0
		iprog=0
		iprogstep=1
		if (.not.allocated(HirBecatm)) then !allocate a dummy array
			nHirBecatm=0
			allocate(HirBecatm(nHirBecatm))
		end if
        
        noldthreads=nthreads
        if (imapfunc==1) then !ESP
            if (iESPcode==2.or.iESPcode==3) then
                call doinitlibreta(1)
                if (isys==1.and.nthreads>12) nthreads=12
            end if
        end if
        write(*,*)
        
        ndosurvtx=count(elimvtx(:)==0)
        ifinish=0;ishowprog=1
        call showprog(0,100)
        ntmp=floor(ndosurvtx/100D0)
		!$OMP PARALLEL DO SHARED(survtx,ifinish,ishowprog) PRIVATE(icyc) schedule(dynamic) NUM_THREADS(nthreads)
		do icyc=1,nsurvtx
			if (elimvtx(icyc)==1) cycle
			survtx(icyc)%value=calcmapfunc(imapfunc,survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,nHirBecatm,HirBecatm)
			if (ntmp/=0) then
				!$OMP CRITICAL
				ifinish=ifinish+1
				ishowprog=mod(ifinish,ntmp)
				if (ishowprog==0) call showprog(floor(100D0*ifinish/ndosurvtx),100)
				!$OMP END CRITICAL
            end if
		end do
		!$OMP END PARALLEL DO
		if (ishowprog/=0) call showprog(100,100)
        
        !Focal-point approximation for evaluating ESP
        !"High method / large BS" = "High method / small BS" (current survtx%value) + "Low method / large BS" (FPAfile2) - "Low method / small BS" (FPAfile1)
		if (iFPA==1) then
			if (imapfunc==1) then
				write(*,*)
				write(*,*) "Doing Focal-point approximation for evaluating ESP..."
				write(*,"(a)") " Loading "//trim(FPAfile1)
				write(*,*) "Calculating ESP..."
				call dealloall(0)
                call readinfile(FPAfile1,1)
                !$OMP PARALLEL DO SHARED(survtx) PRIVATE(icyc) schedule(dynamic) NUM_THREADS(nthreads)
				do icyc=1,nsurvtx
					if (elimvtx(icyc)==1) cycle
					survtx(icyc)%value=survtx(icyc)%value - calcmapfunc(imapfunc,survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,nHirBecatm,HirBecatm)
				end do
				!$OMP END PARALLEL DO
				write(*,"(a)") " Loading "//trim(FPAfile2)
				write(*,*) "Calculating ESP..."
				call dealloall(0)
                call readinfile(FPAfile2,1)
                !$OMP PARALLEL DO SHARED(survtx) PRIVATE(icyc) schedule(dynamic) NUM_THREADS(nthreads)
				do icyc=1,nsurvtx
					if (elimvtx(icyc)==1) cycle
					survtx(icyc)%value=survtx(icyc)%value + calcmapfunc(imapfunc,survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,nHirBecatm,HirBecatm)
				end do
				!$OMP END PARALLEL DO
				call dealloall(0)
                write(*,*) "Using focal-point approximation to evaluate ESP has been finished"
                write(*,"(a)") " Reloading "//trim(firstfilename)
                call readinfile(firstfilename,1)
            else
				write(*,"(a)") " Error: Focal-point approximation can be used only when mapped function is ESP. This calculation is skipped"
            end if
        end if
        nthreads=noldthreads
	end if
	call walltime(iwalltime2)
	write(*,"(' Calculation of mapped function took up wall clock time',i10,' s')") iwalltime2-iwalltime1

else if (ireadextmapval==1) then !Load mapped function from external file
	open(10,file="surfptpos.txt",status="replace")
	write(10,"(i10)") ncurrvtx
	do icyc=1,nsurvtx
		if (elimvtx(icyc)==1) cycle
		write(10,"(3f13.7)") survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z
	end do
	close(10)
	write(*,"(/,a)") " The X/Y/Z coordinates of the surface vertices have been exported to surfptpos.txt in current folder, unit is in Bohr. &
	&Now you can use your favourite program to calculate mapped function values at this points"
	write(*,"(/,a)") " Now input the path of the file containing calculated mapped function values. The data will be read in free format. The first row is the number of points, &
	&in following lines the 1/2/3/4 column should correspond to X,Y,Z and mapped function value, respectively. All units must be in a.u."
	do while(.true.)
		read(*,"(a)") c200tmp
		inquire(file=c200tmp,exist=alive)
		if (alive) exit
		write(*,*) "File cannot be found, input again"
	end do
	open(10,file=c200tmp,status="old")
	read(10,*) npttmp !How many points in the file
	if (npttmp/=ncurrvtx) then
		write(*,"(a,i10,a,i10)") " Error: The number of points in this file",npttmp," is inconsistent with the number of vertices on the molecular surface",ncurrvtx
		exit surfanaloop
	end if
	do icyc=1,nsurvtx
		if (elimvtx(icyc)==1) cycle
		read(10,*) rnouse,rnouse,rnouse,survtx(icyc)%value
	end do
	close(10)
	write(*,*) "Loading data finished!"
	
else if (ireadextmapval==2) then !Output positions to intermediate file, let user to calculate data via cubegen and then load data
	open(10,file="cubegenpt.txt",status="replace")
	do icyc=1,nsurvtx
		if (elimvtx(icyc)==1) cycle
		write(10,"(3f13.7)") survtx(icyc)%x*b2a,survtx(icyc)%y*b2a,survtx(icyc)%z*b2a
	end do
	close(10)
	write(*,"(/,a)") " The coordinate of the surface vertices has been outputted to cubegenpt.txt"
	write(*,"(a,/)") " Use for example ""cubegen 0 potential CNT.fch result.cub -5 h < cubegenpt.txt"" to get cubegen output file"
	write(*,*) "Now input the file name of cubegen output file, e.g. D:\g09\result.cub"
	do while(.true.)
		read(*,"(a)") c200tmp
		inquire(file=c200tmp,exist=alive)
		if (alive) exit
		write(*,*) "File cannot be found, input again"
	end do
	open(10,file=c200tmp,status="old")
	do iskip=1,6+ncenter
		read(10,*)
	end do
	do icyc=1,nsurvtx
		if (elimvtx(icyc)==1) cycle
		read(10,*) rnouse,rnouse,rnouse,survtx(icyc)%value
	end do
	close(10)
	write(*,*) "Loading data finished!"
	
else if (ireadextmapval==3) then !Will calculate mapped function by interpolating from external cube file (cubmattmp)
	write(*,*)
	write(*,*) "Outputting template cube file..."
	open(10,file="template.cub",status="replace")
	call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
	close(10)
	write(*,"(/,a)") " The template cube file has been outputted to template.cub in current folder"
	write(*,"(a)") " Now input the name of the cube file representing mapped function, e.g. C:\t.cub. &
	&The grid setting in this cube file must be exactly identical to the template cube file"
	do while(.true.)
		read(*,"(a)") c200tmp
		inquire(file=c200tmp,exist=alive)
		if (alive) exit
		write(*,*) "File cannot be found, input again"
	end do
	call readcubetmp(c200tmp,1,inconsis)
	write(*,*) "Loading data finished!"
	if (inconsis==1) then
		write(*,"(a)") " Warning: The grid setting in the cube file you inputted is not identical to the template cube file! The analysis result may be meaningless!"
		read(*,*)
	end if
	do icyc=1,nsurvtx
		if (elimvtx(icyc)==1) cycle
		survtx(icyc)%value=linintp3d(survtx(icyc)%x,survtx(icyc)%y,survtx(icyc)%z,2)
	end do
end if

!Find maximum and minimum value of surface vertex
itime=0
do icyc=1,nsurvtx
	if (elimvtx(icyc)==1) cycle
	itime=itime+1
	if (itime==1) then
		indsurfmin=icyc
		indsurfmax=icyc
	end if
	if (survtx(icyc)%value>survtx(indsurfmax)%value) indsurfmax=icyc
	if (survtx(icyc)%value<survtx(indsurfmin)%value) indsurfmin=icyc
end do

!Approximately evaluate mapped function value at each facet, use average of values at three verices as the value of facet 
do icyc=1,nsurtri
	if (elimtri(icyc)==1) cycle
	idx1=surtriang(icyc)%idx(1)
	idx2=surtriang(icyc)%idx(2)
	idx3=surtriang(icyc)%idx(3)
	surtriang(icyc)%value=(survtx(idx1)%value+survtx(idx2)%value+survtx(idx3)%value)/3D0
end do
if (size(survtx)>0) then !The number of vertex is zero for empty grid data
	write(*,*)
	write(*,"(' Global surface minimum:',f10.6,' a.u. at',3f11.6,' Ang')") survtx(indsurfmin)%value,survtx(indsurfmin)%x*b2a,survtx(indsurfmin)%y*b2a,survtx(indsurfmin)%z*b2a
	write(*,"(' Global surface maximum:',f10.6,' a.u. at',3f11.6,' Ang')") survtx(indsurfmax)%value,survtx(indsurfmax)%x*b2a,survtx(indsurfmax)%y*b2a,survtx(indsurfmax)%z*b2a
end if
write(*,*)

!Find local minimum of mapped function on surface
!To assign a minimum, not only vertices that directly contact to it must have smaller value, &
!but also the secondary neighbouring vertices must be lower than it
nsurlocmin=0
cmin: do ivtx=1,nsurvtx
	if (elimvtx(ivtx)==1) cycle
	vali=survtx(ivtx)%value
	do jtmp=1,vtxconnpos(ivtx) !Compare value with directly connected neighbours
		jvtx=vtxconn(ivtx,jtmp)
		if (elimvtx(jvtx)==1) cycle
		if (vali>=survtx(jvtx)%value) cycle cmin
		do ktmp=1,vtxconnpos(jvtx) !Compare value with secondary neighbours
			kvtx=vtxconn(jvtx,ktmp)
			if (elimvtx(kvtx)==1) cycle
			if (kvtx==ivtx) cycle
			if (vali>=survtx(kvtx)%value) cycle cmin !If want to retain pseudo-minima, comment this line
		end do
	end do
	nsurlocmin=nsurlocmin+1
	surlocminidx(nsurlocmin)=ivtx
end do cmin
write(*,"(a,i6)") " The number of surface minima:",nsurlocmin
if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
	write(*,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
else
	write(*,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
end if
do i=1,nsurlocmin
	idx=surlocminidx(i)
	char1tmp=' '
	if (idx==indsurfmin) char1tmp='*' !Mark minimum term by asterisk
	if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
		write(*,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
	else
		write(*,"(a,i5,f19.7,3x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
	end if
end do
write(*,*)

!Find local maximum of mapped function on surface
nsurlocmax=0
cmax: do ivtx=1,nsurvtx
	if (elimvtx(ivtx)==1) cycle
	vali=survtx(ivtx)%value
	do jtmp=1,vtxconnpos(ivtx)
		jvtx=vtxconn(ivtx,jtmp)
		if (elimvtx(jvtx)==1) cycle
		if (vali<=survtx(jvtx)%value) cycle cmax
		do ktmp=1,vtxconnpos(jvtx)
			kvtx=vtxconn(jvtx,ktmp)
			if (elimvtx(kvtx)==1) cycle
			if (kvtx==ivtx) cycle
			if (vali<=survtx(kvtx)%value) cycle cmax !If want to retain pseudo-maxima, comment this line
		end do
	end do
	nsurlocmax=nsurlocmax+1
	surlocmaxidx(nsurlocmax)=ivtx
end do cmax
write(*,"(a,i6)") " The number of surface maxima:",nsurlocmax
if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
	write(*,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
else
	write(*,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
end if
do i=1,nsurlocmax
	idx=surlocmaxidx(i) !Convert to absolute surface vertice index
	char1tmp=' '
	if (idx==indsurfmax) char1tmp='*' !Mark maximum term by asterisk
	if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
		write(*,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
	else
		write(*,"(a,i5,f19.7,3x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
	end if
end do


!============ Perform statistic
10 write(*,*)
write(*,*) "      ================= Summary of surface analysis ================="
write(*,*)
totvalang3=totvol*b2a**3
totmass=sum(atmwei(a%index))
densmass=totmass/avogacst*(1D24/totvalang3)
write(*,"(' Volume:',f12.5,' Bohr^3  (',f10.5,' Angstrom^3)')") totvol,totvalang3
write(*,"(' Estimated density according to mass and volume (M/V):',f10.4,' g/cm^3')") densmass

if (iskipmapfunc==1) then !Only show surface area, do not print any quantity related to mapped function
    write(*,"(' Overall surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareaall,surfareaall*b2a*b2a
    goto 20
end if

if (imapfunc==1.or.imapfunc==3) then !ESP or ESP from atomic charges
	write(*,"(' Minimal value:',f13.5,' kcal/mol   Maximal value:',f13.5,' kcal/mol')") survtx(indsurfmin)%value*au2kcal,survtx(indsurfmax)%value*au2kcal
else if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) then !ALIE or LEA or LEAE
	write(*,"(' Minimal value:',f13.5,' eV,   Maximal value:',f13.5,' eV')") survtx(indsurfmin)%value*au2eV,survtx(indsurfmax)%value*au2eV
else
	write(*,"(' Minimal value:',f16.8,'    Maximal value:',f16.8)") survtx(indsurfmin)%value,survtx(indsurfmax)%value
end if
surfareapos=0D0
surfareaneg=0D0
surfareapol=0D0
surfareanonpol=0D0
do i=1,nsurtri
	if (elimtri(i)==1) cycle
	if (surtriang(i)%value>=0) then
		surfareapos=surfareapos+surtriang(i)%area
	else
		surfareaneg=surfareaneg+surtriang(i)%area
	end if
    if (abs(surtriang(i)%value*au2kcal)<=10) then
        surfareanonpol=surfareanonpol+surtriang(i)%area
    else
        surfareapol=surfareapol+surtriang(i)%area
    end if
end do
write(*,"(' Overall surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareaall,surfareaall*b2a*b2a
write(*,"(' Positive surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareapos,surfareapos*b2a*b2a
write(*,"(' Negative surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") surfareaneg,surfareaneg*b2a*b2a

sumposval=0D0
sumnegval=0D0
do i=1,nsurtri
	if (elimtri(i)==1) cycle
	valtmp=surtriang(i)%value
	if (valtmp>=0D0) then
		sumposval=sumposval+valtmp*surtriang(i)%area
	else
		sumnegval=sumnegval+valtmp*surtriang(i)%area
	end if
end do
sumallval=sumposval+sumnegval
avgpos=sumposval/surfareapos
avgneg=sumnegval/surfareaneg
avgall=sumallval/surfareaall
if (imapfunc==1.or.imapfunc==3) then !ESP or ESP from atomic charges
	write(*,"(' Overall average value: ',f13.8,' a.u. (',f13.5,' kcal/mol)')") avgall,avgall*au2kcal
	write(*,"(' Positive average value:',f13.8,' a.u. (',f13.5,' kcal/mol)')") avgpos,avgpos*au2kcal
	write(*,"(' Negative average value:',f13.8,' a.u. (',f13.5,' kcal/mol)')") avgneg,avgneg*au2kcal
else if (imapfunc==2) then !ALIE is always positive
	write(*,"(' Average value: ',f13.8,' a.u. (',f13.5,' eV,',f13.5,' kcal/mol)')") avgpos,avgpos*au2eV,avgpos*au2kcal
else if (imapfunc==4.or.imapfunc==-4) then !LEA or LEAE
	write(*,"(' Overall average value: ',f11.7,' a.u. (',f10.5,' eV,',f13.5,' kcal/mol)')") avgall,avgall*au2eV,avgall*au2kcal
	write(*,"(' Positive average value:',f11.7,' a.u. (',f10.5,' eV,',f13.5,' kcal/mol)')") avgpos,avgpos*au2eV,avgpos*au2kcal
	write(*,"(' Negative average value:',f11.7,' a.u. (',f10.5,' eV,',f13.5,' kcal/mol)')") avgneg,avgneg*au2eV,avgneg*au2kcal
else
	write(*,"(' Overall average value: ',f20.8)") avgall
	write(*,"(' Positive average value:',f20.8)") avgpos
	write(*,"(' Negative average value:',f20.8)") avgneg
end if

chgsep=0D0
varipos=0D0
varineg=0D0
polaridx=0D0
do i=1,nsurtri
	if (elimtri(i)==1) cycle
	valtmp=surtriang(i)%value
    areatmp=surtriang(i)%area
	if (valtmp>=0D0) then
		varipos=varipos+areatmp*(valtmp-avgpos)**2
	else
		varineg=varineg+areatmp*(valtmp-avgneg)**2
	end if
	chgsep=chgsep+abs(valtmp-avgall)*areatmp
    polaridx=polaridx+abs(valtmp)*areatmp
end do
if (surfareapos/=0D0) varipos=varipos/surfareapos
if (surfareaneg/=0D0) varineg=varineg/surfareaneg
variall=varipos+varineg
chgsep=chgsep/surfareaall
polaridx=polaridx/surfareaall
if (imapfunc==1.or.imapfunc==3) then
	balencechg=varipos*varineg/(varipos+varineg)**2
	write(*,"(' Overall variance (sigma^2_tot):',f12.8,' a.u.^2 (',f12.5,' (kcal/mol)^2)')") variall,variall*au2kcal**2
	write(*,"(' Positive variance:     ',f13.8,' a.u.^2 (',f13.5,' (kcal/mol)^2)')") varipos,varipos*au2kcal**2
	write(*,"(' Negative variance:     ',f13.8,' a.u.^2 (',f13.5,' (kcal/mol)^2)')") varineg,varineg*au2kcal**2
	write(*,"(' Balance of charges (nu):',f13.8)") balencechg
	write(*,"(' Product of sigma^2_tot and nu: ',f12.8,' a.u.^2 (',f11.5,' (kcal/mol)^2)')") balencechg*variall,balencechg*variall*au2kcal**2
	write(*,"(' Internal charge separation (Pi):',f13.8,' a.u. (',f13.5,' kcal/mol)')") chgsep,chgsep*au2kcal
    write(*,"(' Molecular polarity index (MPI):',f13.8,' eV (',f13.5,' kcal/mol)')") polaridx*au2eV,polaridx*au2kcal
    write(*,"(' Nonpolar surface area (|ESP| <= 10 kcal/mol):',f10.2,' Angstrom^2  (',f6.2,' %)')") surfareanonpol*b2a*b2a,surfareanonpol/surfareaall*100
    write(*,"(' Polar surface area (|ESP| > 10 kcal/mol):    ',f10.2,' Angstrom^2  (',f6.2,' %)')") surfareapol*b2a*b2a,surfareapol/surfareaall*100
else if (imapfunc==2) then
	write(*,"(' Variance:  ',f13.8,' a.u.^2  (',f13.5,' eV^2,',E13.5,' kcal/mol^2)')") variall,variall*au2eV**2,variall*au2kcal**2
else if (imapfunc==4.or.imapfunc==-4) then
	write(*,"(' Overall variance: ',f10.6,' a.u.^2  (',f10.5,' eV^2,',E12.5,' kcal/mol^2)')") variall,variall*au2eV**2,variall*au2kcal**2
	write(*,"(' Positive variance:',f10.6,' a.u.^2  (',f10.5,' eV^2,',E12.5,' kcal/mol^2)')") varipos,varipos*au2eV**2,varipos*au2kcal**2
	write(*,"(' Negative variance:',f10.6,' a.u.^2  (',f10.5,' eV^2,',E12.5,' kcal/mol^2)')") varineg,varineg*au2eV**2,varineg*au2kcal**2
else
	write(*,"(' Overall variance:  ',f20.10)") variall
	write(*,"(' Positive variance: ',f20.10)") varipos
	write(*,"(' Negative variance: ',f20.10)") varineg
end if

!Calculate Skewness, as requested in http://sobereva.com/wfnbbs/viewtopic.php?id=925
skewall=0
skewpos=0
skewneg=0
do i=1,nsurtri
	if (elimtri(i)==1) cycle
	valtmp=surtriang(i)%value
    areatmp=surtriang(i)%area
    skewall=skewall+areatmp*(valtmp-avgall)**3
	if (valtmp>=0D0) then
		skewpos=skewpos+areatmp*(valtmp-avgpos)**3
	else
		skewneg=skewneg+areatmp*(valtmp-avgneg)**3
	end if
end do
skewall=skewall/surfareaall/variall**(3D0/2D0)
skewpos=skewpos/surfareapos/varipos**(3D0/2D0)
skewneg=skewneg/surfareaneg/varineg**(3D0/2D0)
write(*,"(' Overall skewness: ',f20.10)") skewall
if (surfareapos>0) write(*,"(' Positive skewness:',f20.10)") skewpos
if (surfareaneg>0) write(*,"(' Negative skewness:',f20.10)") skewneg

write(*,*)
write(*,*) "Surface analysis finished!"

call walltime(iclktime2)
write(*,"(' Total wall clock time passed during this task:',i6,' s')") iclktime2-iclktime1
if (imapfunc/=0.and.imapfunc/=4.and.imapfunc/=20.and.imapfunc/=21.and.imapfunc/=22) call delvirorb_back(1)
if (allocated(b)) call del_GTFuniq !Destory unique GTF informtaion
if (imapfunc==1) write(*,"(a)") " Citation of molecular polarity index (MPI): Carbon, 171, 514 (2021) DOI: 10.1016/j.carbon.2020.09.048"


!============= Post-processing stage
!For charged system, magnitude of ESP at vertices may be quite large, therefore use eV instead of kcal/mol when outputting .pdb files
if (imapfunc==1.or.imapfunc==3) then
	iESPev=0
	if (survtx(indsurfmin)%value*au2kcal<=-99.99D0.or.survtx(indsurfmax)%value*au2kcal>999.99D0) iESPev=1
end if
20 textheigh=30D0
ratioatmsphere=1D0
bondradius=0.2D0
ishowatmlab=1
ishowaxis=1
ishowlocminlab=0
ishowlocmaxlab=0
ishowlocminpos=1
ishowlocmaxpos=1
do while(.true.)
	write(*,*)
	write(*,*) "                   ---------- Post-processing menu ----------"
	write(*,*) "-3 Visualize the surface"
	write(*,*) "-2 Export the grid data to surf.cub in current folder"
	write(*,*) "-1 Return to upper level menu"
	write(*,*) "0 View molecular structure, surface minima and maxima"
	write(*,*) "1 Export surface extrema as surfanalysis.txt in current folder"
	write(*,*) "2 Export surface extrema as surfanalysis.pdb in current folder"
	write(*,*) "3 Discard surface minima in certain value range"
	write(*,*) "4 Discard surface maxima in certain value range"
	write(*,*) "5 Export molecule as pdb format file"
	write(*,*) "6 Export all surface vertices to vtx.pdb in current folder" !if 66, also output connectivity
	write(*,*) "7 Export all surface vertices to vtx.txt in current folder"
	write(*,*) "8 Export all surface vertices and surface extrema as vtx.pqr and extrema.pqr"
	write(*,*) "9 Output surface area in specific value range of mapped function"
	write(*,*) "10 Output the closest and farthest distance between the surface and a point"
	write(*,*) "11 Output surface properties of each atom"
	write(*,*) "12 Output surface properties of specific fragment"
	write(*,*) "13 Calculate grid data of mapped function and export it to mapfunc.cub"
    write(*,*) "14 Calculate area and function average in a region around a surface extreme"
    write(*,*) "15 Basin-like partition of surface and calculate areas"
! 	write(*,*) "16 Export center of surface facets as pdb file" !Can also output to xyz file
    write(*,*) "18 Discard some surface extrema by inputting their indices"
    write(*,*) "19 Merge some surface extrema and take their average position"
	if (isurftype==5.or.isurftype==6) write(*,*) "20 Fingerprint plot and local contact analyses"
	read(*,*) isel
	
	if (isel==-3) then
		sur_value=surfisoval
		call drawisosurgui(1)
		
	else if (isel==-2) then
		write(*,*) "Exporting surf.cub..."
		open(10,file="surf.cub",status="replace")
		call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done, the grid data has been exported to surf.cub in current folder"
		
	else if (isel==-1) then
		exit
		
	else if (isel==0.or.isel==1) then
		if (isel==0) then
			ides=6
		else if (isel==1) then
			ides=10
			open(ides,file="surfanalysis.txt",status="replace")
		end if
		write(ides,"(a,i6)") " Number of surface minima:",count(surlocminidx(1:nsurlocmin)/=0)
		if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
			write(ides,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
		else
			write(ides,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
		end if
		do i=1,nsurlocmin
			idx=surlocminidx(i)
			if (idx==0) cycle !The extreme has already been discarded
			char1tmp=" "
			if (idx==indsurfmin) char1tmp="*"
			if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
				write(ides,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,&
				survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
			else
				write(ides,"(a,i5,f19.7,3x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
			end if
		end do
		write(ides,*)
		write(ides,"(a,i6)") " Number of surface maxima:",count(surlocmaxidx(1:nsurlocmax)/=0)
		if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
			write(ides,*) "  #       a.u.         eV      kcal/mol           X/Y/Z coordinate(Angstrom)"
		else
			write(ides,*) "  #             Value           X/Y/Z coordinate(Angstrom)"
		end if
		do i=1,nsurlocmax
			idx=surlocmaxidx(i)
			if (idx==0) cycle !The extreme has already been discarded
			char1tmp=" "
			if (idx==indsurfmax) char1tmp="*"
			if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
				write(ides,"(a,i5,f12.8,f12.6,f12.6,4x,3f11.6)") char1tmp,i,survtx(idx)%value,&
				survtx(idx)%value*au2eV,survtx(idx)%value*au2kcal,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
			else
				write(ides,"(a,i5,f19.7,3x,3f11.6)") char1tmp,i,survtx(idx)%value,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
			end if
		end do
		if (isel==0) then
			call drawsurfanalysis
		else if (isel==1) then
			close(10)
			write(*,"(a)") " Results have been outputted to surfanalysis.txt in current folder"
		end if
		
	else if (isel==2) then
		!Output positions of local maximum (as carbon) and local minimum (as oxygen) to pdb file
		open(10,file="surfanalysis.pdb",status="replace")
		if (ifPBC>0) then
			call getcellabc(asize,bsize,csize,alpha,beta,gamma)
			write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
		end if
		if (imapfunc==1.or.imapfunc==3) then !ESP
			if (iESPev==0) then
				write(10,"(a)") "REMARK   Unit of B-factor field (i.e. ESP) is kcal/mol"
			else if (iESPev==1) then
				write(10,"(a)") "REMARK   Unit of B-factor field (i.e. ESP) is eV"
			end if
        else if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) then
            if (imapfunc==2) write(10,"(a)") "REMARK   Unit of B-factor field (ALIE) is eV"
            if (imapfunc==4) write(10,"(a)") "REMARK   Unit of B-factor field (LEA) is eV"
            if (imapfunc==-4) write(10,"(a)") "REMARK   Unit of B-factor field (LEAE) is eV"
		end if
        write(10,"(a)") "REMARK   Carbon: Surface maximum    Oxygen: surface minimum"
		do i=1,nsurlocmax
			idx=surlocmaxidx(i)
			if (idx==0) cycle
			if (imapfunc==1.or.imapfunc==3) then
				if (iESPev==0) then
					tmpval=survtx(idx)%value*au2kcal
				else if (iESPev==1) then
					tmpval=survtx(idx)%value*au2eV
				end if
			else if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) then
				tmpval=survtx(idx)%value*au2eV
			else
				tmpval=survtx(idx)%value
			end if
			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
			"HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a,1.0,tmpval,"C "
		end do
		do i=1,nsurlocmin
			idx=surlocminidx(i)
			if (idx==0) cycle !The extreme has already been discarded
			if (imapfunc==1.or.imapfunc==3) then
				if (iESPev==0) then
					tmpval=survtx(idx)%value*au2kcal
				else if (iESPev==1) then
					tmpval=survtx(idx)%value*au2eV
				end if
			else if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) then
				tmpval=survtx(idx)%value*au2eV
			else
				tmpval=survtx(idx)%value
			end if
			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
			"HETATM",i,' '//"O "//' ',"MOL",'A',1,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a,1.0,tmpval,"O "
		end do
		write(10,"('END')")
		close(10)
		write(*,"(a)") " Results have been outputted to surfanalysis.pdb in current folder"
		write(*,"(a)",advance='no') " Carbons and oxygens correspond to local maximum and minimum points respectively."
		if (imapfunc==1.or.imapfunc==3) then
			if (iESPev==0) then
				write(*,"(a)") " B-factor field records mapped function value in kcal/mol"
			else if (iESPev==1) then
				write(*,"(a)") " B-factor field records mapped function value in eV"
			end if
		else if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) then
			write(*,"(a)") " B-factor field records mapped function value in eV"
		else
			write(*,"(a)") " B-factor field records mapped function value in original unit"
		end if
		
	else if (isel==3) then
		write(*,*) "Input value range (in a.u.),  e.g. 0.8,90"
        write(*,"(a)") " Hint: If you simply input letter ""d"", all surface minima with positive value will be removed"
        read(*,"(a)") c200tmp
        if (index(c200tmp,'d')/=0) then
            ncount=0
		    do i=1,nsurlocmin
			    idx=surlocminidx(i)
			    if (idx/=0.and.survtx(idx)%value>0) then
                    surlocminidx(i)=0
                    ncount=ncount+1
                end if
		    end do
            write(*,"(i6,a)") ncount," surface minima have been removed"
        else
		    read(c200tmp,*) vallowlim,valhighlim
		    do i=1,nsurlocmin
			    idx=surlocminidx(i)
			    if (idx==0) cycle !The extreme has already been discarded
			    if (survtx(idx)%value>vallowlim.and.survtx(idx)%value<valhighlim) surlocminidx(i)=0
		    end do
		    write(*,"(' Surface minima with value from',f10.5,' to',f10.5,' have been discarded')") vallowlim,valhighlim
        end if
		
	else if (isel==4) then
		write(*,*) "Input value range (in a.u.),  e.g. -10,20.5"
        write(*,"(a)") " Hint: If you simply input letter ""d"", all surface maxima with negative value will be removed"
        read(*,"(a)") c200tmp
        if (index(c200tmp,'d')/=0) then
            ncount=0
		    do i=1,nsurlocmax
			    idx=surlocmaxidx(i)
			    if (idx/=0.and.survtx(idx)%value<0) then
                    surlocmaxidx(i)=0
                    ncount=ncount+1
                end if
		    end do
            write(*,"(i6,a)") ncount," surface maxima have been removed"
        else
		    read(c200tmp,*) vallowlim,valhighlim
		    do i=1,nsurlocmax
			    idx=surlocmaxidx(i)
			    if (idx==0) cycle !The extreme has already been discarded
			    if (survtx(idx)%value>vallowlim.and.survtx(idx)%value<valhighlim) surlocmaxidx(i)=0
		    end do
		    write(*,"(' Surface maxima with value from',f10.5,' to',f10.5,' have been discarded')") vallowlim,valhighlim
        end if
		
	else if (isel==5) then
		write(*,*) "Input the filename you want to save to, e.g. C:\K-ON\Mio.pdb"
		read(*,"(a)") pdbfilename
		call outpdb(pdbfilename,10)
		
	else if (isel==6.or.isel==66) then !If 66, also outputting connectivity
		!The density value in intermolecular region is always very small, and beta column in .pdb file is very narrow, we need scale it by a factor
		if (imapfunc==11.or.imapfunc==12) then
			write(*,*) "Input the factor for scaling the density, e.g. 100"
			read(*,*) tmpfac
		end if
		open(10,file="vtx.pdb",status="replace")
		write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface vertices')") nsurvtx
		if (ifPBC>0) then
			call getcellabc(asize,bsize,csize,alpha,beta,gamma)
			write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
		end if
		if (imapfunc==1.or.imapfunc==3) then
			if (iESPev==0) then
				write(10,"(a)") "REMARK   Unit of B-factor field (i.e. ESP) is kcal/mol"
			else if (iESPev==1) then
				write(10,"(a)") "REMARK   Unit of B-factor field (i.e. ESP) is eV"
			end if
		end if
		do i=1,nsurvtx
			if (imapfunc==1.or.imapfunc==3) then
				if (iESPev==0) then
					tmpfuncval=survtx(i)%value*au2kcal
				else if (iESPev==1) then
					tmpfuncval=survtx(i)%value*au2eV
				end if
			else if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) then
				tmpfuncval=survtx(i)%value*au2eV
			else if (imapfunc==11.or.imapfunc==12) then
				tmpfuncval=survtx(i)%value*tmpfac
			else
				tmpfuncval=survtx(i)%value
			end if
			!Scale data to avoid exceeding limit, because B-factor field only have three integer positions
			if (tmpfuncval>999.99D0) tmpfuncval=999.99D0
			if (tmpfuncval<=-99.99D0) tmpfuncval=-99.99D0
			if (elimvtx(i)==0) then !Output normal vertices as carbon atoms
				write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
				"HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,1.0,tmpfuncval,"C "
			else !Output eliminated vertices as oxygens when hidden option 66 is selected
			    if (isel==66) then
				    write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
				    "HETATM",i,' '//"O "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,1.0,tmpfuncval,"O "
				end if
			end if
		end do
		if (isel==66) then !Hidden option, also output connectivity of vertices, which can be visualized in VMD using the way introduced in http://sobereva.com/121
		    do i=1,nsurvtx
			    if (elimvtx(i)==1) cycle
			    write(10,"('CONECT',i6)",advance='no') i 
			    do j=1,vtxconnpos(i)
				    if (elimvtx(vtxconn(i,j))==0) write(10,"(i6)",advance='no') vtxconn(i,j)
			    end do
			    write(10,*)
		    end do
		end if
		write(10,"('END')")
		close(10)
		write(*,"(a)") " Surface vertices have been outputted to vtx.pdb in current folder"
		if (imapfunc==1.or.imapfunc==3) then
			if (iESPev==0) write(*,"(a)") " B-factor field records mapped function value in kcal/mol"
			if (iESPev==1) write(*,"(a)") " B-factor field records mapped function value in eV"
		else if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) then
			write(*,"(a)") " B-factor field records mapped function value in eV"
		else
			write(*,"(a)") " B-factor field records mapped function value in original unit"
		end if
		if (isel==6) then
! 		    write(*,"(a)") " Note: The eliminated vertices and connectivity are not outputted. If you would like to output them you should select option 66 (a hidden option)"
		else if (isel==66) then
		    write(*,"(a)") " Note: Carbons and oxygens correspond to actually used and eliminated vertices respecitvely. CONECT field records connectivity between the vertices."
		end if
		
		!Output vertices with three connnections
! 		open(10,file="vtx3.pdb",status="replace")
! 		do i=1,nsurvtx
! 			if (elimvtx(i)==1) cycle
! 			numconn=0
! 			do j=1,vtxconnpos(i)
! 				if (elimvtx(vtxconn(i,j))==0) numconn=numconn+1
! 			end do
! 			if (numconn==3) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
! 				"HETATM",i,' '//"N "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,1.0,0.0,"N "
! 		end do
! 		close(10)

	else if (isel==7) then
		open(10,file="vtx.txt",status="replace")
		write(10,"(i10)") ncurrvtx
		do i=1,nsurvtx
			if (elimvtx(i)==0) then
				if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
					write(10,"(3f13.7,4x,3f16.10)") survtx(i)%x,survtx(i)%y,survtx(i)%z,survtx(i)%value,survtx(i)%value*au2eV,survtx(i)%value*au2kcal
				else
					write(10,"(3f13.7,4x,f18.10)") survtx(i)%x,survtx(i)%y,survtx(i)%z,survtx(i)%value
				end if
			end if
		end do
		write(*,"(a)") " Done, all surface vertices (not including the ones have been eliminated) have been outputted to vtx.txt in current folder"
		if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
			write(*,"(a)") " The first line is the number of points. Column 1,2,3,4,5,6 respectively correspond to coordinate of X/Y/Z in Bohr, mapped function in a.u., eV and kcal/mol"
		else
			write(*,"(a)") " The first line is the number of points. Column 1,2,3,4 respectively correspond to coordinate of X/Y/Z in Bohr and mapped function in original unit"
		end if
		close(10)
        
    else if (isel==8) then !Export all surface vertices and surface extrema as vtx.pqr and extrema.pqr
    
    	open(10,file="extrema.pqr",status="replace")
		write(10,"(a)") "REMARK   The third last column is function values in a.u."
		do i=1,nsurlocmax
			idx=surlocmaxidx(i)
            if (idx==0) cycle !The extreme has already been discarded
			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
			"HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a,survtx(idx)%value,1.0," C"
		end do
		do i=1,nsurlocmin
			idx=surlocminidx(i)
            if (idx==0) cycle !The extreme has already been discarded
			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
			"HETATM",i,' '//"O "//' ',"MOL",'A',1,survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a,survtx(idx)%value,1.0," O"
		end do
		write(10,"('END')")
		close(10)
		write(*,"(a)") " Surface extrema have been outputted to extrema.pqr in current folder. &
        &Carbons and oxygens correspond to local maximum and minimum points respectively. &
        &The atomic charge column (i.e. third last column) corresponds to function value in a.u."
        
		open(10,file="vtx.pqr",status="replace")
		write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface vertices')") nsurvtx
		write(10,"(a)") "REMARK   The third last column is function values in a.u."
		do i=1,nsurvtx
            if (elimvtx(i)==0) write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
			"HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,survtx(i)%value,1.0," C"
		end do
		write(10,"('END')")
		close(10)
		write(*,"(/,a)") " Surface vertices have been outputted to vtx.pqr in current folder. &
        &The atomic charge column (i.e. third last column) corresponds to function value in a.u."
		
	else if (isel==9) then
		write(*,"(a)") " Input atomic indices to define the fragment. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
		write(*,"(a)") " If input ""all"" or press ENTER button directly, then the whole system will be analyzed"
		read(*,"(a)") c2000tmp
	    if (allocated(surtrifrag)) deallocate(surtrifrag)
	    allocate(surtrifrag(nsurtri))
		if (index(c2000tmp,'all')/=0.or.c2000tmp==" ") then
		    surtrifrag=1 !All vertices are taken into account
		else
		    call str2arr(c2000tmp,nsurfragatm,surfrag) !surfrag contains nsurfragatm elements, they are the member of the user-defined fragment
		    do iatm=1,ncenter !Generate table ifatmfrag, determine if an atom is belong to the user-defined fragment
			    if ( any(surfrag(1:nsurfragatm)==iatm) ) then
				    ifatmfrag(iatm)=1
			    else
				    ifatmfrag(iatm)=0
			    end if
		    end do
		    !Determine the triangles belong to which atom, use Voronoi-like partition
		    surtrifrag=0 !The eliminated vertices will be attributed to null (0)
		    do itri=1,nsurtri
			    if (elimtri(itri)==1) cycle
			    effdistmax=-1D50
			    surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
			    surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
			    surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
			    do iatm=1,ncenter
				    disttmp=dsqrt((surtrix-a(iatm)%x)**2+(surtriy-a(iatm)%y)**2+(surtriz-a(iatm)%z)**2)
		 		    if (imolsurparmode==1) then
		 			    vdwrtmp=1D0 !Use original Voronoi partition
		 		    else
					    vdwrtmp=vdwr(a(iatm)%index)
				    end if
				    effdisttmp=1-disttmp/vdwrtmp
				    if (effdisttmp>effdistmax) then
					    surtrifrag(itri)=iatm
					    effdistmax=effdisttmp
				    end if
			    end do
			    surtrifrag(itri)=ifatmfrag(surtrifrag(itri)) ! Then surtrifrag(itri)=1/0 means this itri belongs / doesn't belong to user-defined fragment
		    end do
	    end if
		write(*,"(a)") " Input lower and upper limits of the mapped function for statistical analysis (unit will be selected later), e.g. -45,50.5"
		read(*,*) rangelow,rangehigh
		write(*,*) "Input the number of statistical intervals, e.g. 5"
		read(*,*) nintval
		convunit=1D0
		if (imapfunc==1.or.imapfunc==2.or.imapfunc==3.or.imapfunc==4.or.imapfunc==-4) then
			write(*,*) "The inputted range is in which unit?"
			write(*,*) "1: a.u.     2: eV     3: kcal/mol"
			read(*,*) iunit
			if (iunit==2) convunit=au2eV
			if (iunit==3) convunit=au2kcal
		end if
		step=(rangehigh-rangelow)/nintval
		write(*,*) "Note: Area unit is in Angstrom^2"
		write(*,*) "     Begin        End       Center       Area         %"
		areatmptot=0D0
		do istep=1,nintval
			tmplow=rangelow+(istep-1)*step
			tmphigh=rangelow+istep*step
			areatmp=0D0
			do itri=1,nsurtri
				if (elimtri(itri)==1) cycle
				if (surtrifrag(itri)==0) cycle
				if (surtriang(itri)%value*convunit>=tmplow.and.surtriang(itri)%value*convunit<tmphigh) areatmp=areatmp+surtriang(itri)%area
			end do
			areatmp=areatmp*b2a**2
			areatmptot=areatmptot+areatmp
			write(*,"(5f12.4)") tmplow,tmphigh,(tmphigh+tmplow)/2D0,areatmp,areatmp/(surfareaall*b2a**2)*100
		end do
		write(*,"(' Sum:',31x,2f12.4)") areatmptot,areatmptot/(surfareaall*b2a**2)*100
        write(*,*)
        write(*,"(a)") " Citing below paper is highly suggested, in which the surface area analysis for mapped function was proposed and employed for the first time"
        write(*,*) "Sergio Manzetti, Tian Lu, J. Phys. Org. Chem., 26, 473 (2013)"
        write(*,*) "  More examples of this kind of analysis you may also cite:"
        write(*,*) "Struct. Chem., 25, 1521¨C1533 (2014)"
        write(*,*) "Carbon, 171, 514-523 (2021)"
		
	else if (isel==10) then
		write(*,*) "Input XYZ coordinate of the point (in Angstrom), e.g. 0.2,-4.3,1.66"
		write(*,*) "or input atomic index to use its nuclear position as the point, e.g. 5" 
		write(*,*) "or input ""g"" to set geometry center of present system as the point"
		write(*,"(a)") " If input ""f"", then the farthest distance between all surface points will be outputted"
		read(*,"(a)") c200tmp
		if (index(c200tmp,'f')/=0) then
			write(*,*) "Calculating, please wait..."
			distmax2=0
            idx=0;jdx=0
			do ivtx=1,nsurvtx
				do jvtx=ivtx+1,nsurvtx
					dist2=(survtx(ivtx)%x-survtx(jvtx)%x)**2+(survtx(ivtx)%y-survtx(jvtx)%y)**2+(survtx(ivtx)%z-survtx(jvtx)%z)**2
					if (dist2>distmax2) then
                        distmax2=dist2
                        idx=ivtx
                        jdx=jvtx
                    end if
				end do
			end do
			write(*,"(' The farthest distance is:',f12.6,' Bohr (',f12.6,' Angstrom)')") dsqrt(distmax2),dsqrt(distmax2)*b2a
            write(*,"(' X/Y/Z of point 1: ',3f12.6,' Angstrom')") survtx(idx)%x*b2a,survtx(idx)%y*b2a,survtx(idx)%z*b2a
            write(*,"(' X/Y/Z of point 2: ',3f12.6,' Angstrom')") survtx(jdx)%x*b2a,survtx(jdx)%y*b2a,survtx(jdx)%z*b2a
            write(*,"(' X/Y/Z of midpoint:',3f12.6,' Angstrom')") &
            (survtx(idx)%x+survtx(jdx)%x)*b2a/2,(survtx(idx)%y+survtx(jdx)%y)*b2a/2,(survtx(idx)%z+survtx(jdx)%z)*b2a/2
		else
			if (index(c200tmp,',')/=0) then
				read(c200tmp,*) tmpx,tmpy,tmpz
				tmpx=tmpx/b2a
				tmpy=tmpy/b2a
				tmpz=tmpz/b2a
			else if (index(c200tmp,'g')/=0) then
				tmpx=sum(a(:)%x)/ncenter
				tmpy=sum(a(:)%y)/ncenter
				tmpz=sum(a(:)%z)/ncenter
			else
				read(c200tmp,*) iatm
				tmpx=a(iatm)%x
				tmpy=a(iatm)%y
				tmpz=a(iatm)%z
			end if
			do i=1,nsurvtx
				dist=dsqrt( (survtx(i)%x-tmpx)**2+(survtx(i)%y-tmpy)**2+(survtx(i)%z-tmpz)**2 )
				if (i==1.or.dist<distmin) distmin=dist
				if (i==1.or.dist>distmax) distmax=dist
			end do
			if (index(c200tmp,',')==0) write(*,"(' The XYZ coordinate of the point you chosen:',/,3f12.6,' Angstrom')") tmpx*b2a,tmpy*b2a,tmpz*b2a
			write(*,"(a,f12.6,a,f12.6,a)") " The closest distance to the point:",distmin," Bohr (",distmin*b2a," Angstrom)"
			write(*,"(a,f12.6,a,f12.6,a)") " The farthest distance to the point:",distmax," Bohr (",distmax*b2a," Angstrom)"
		end if
		
	else if (isel==11.or.isel==12) then !Separate the properties of the whole molecular surface into atomic or fragment local surface
		if (isel==11) then
			nsurfrag=ncenter !Each fragment corresponds to an atom
		else if (abs(isel)==12) then
			nsurfrag=1 !User define only one fragment
			write(*,"(a)") " Input atomic indices to define the fragment. e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 will constitute the fragment"
			read(*,"(a)") c2000tmp
			call str2arr(c2000tmp,nsurfragatm,surfrag) !surfrag contains nsurfragatm elements, they are the member of the user-defined fragment
			do iatm=1,ncenter !Generate table ifatmfrag, determine if an atom is belong to the user-defined fragment
				if ( any(surfrag(1:nsurfragatm)==iatm) ) then
					ifatmfrag(iatm)=1
				else
					ifatmfrag(iatm)=0
				end if
			end do
		end if
		!Determine the triangles belong to which atom, use Voronoi-like partition
		if (allocated(surtrifrag)) deallocate(surtrifrag)
		allocate(surtrifrag(nsurtri))
		surtrifrag=0 !The eliminated vertices will be attributed to null (0)
		do itri=1,nsurtri
			if (elimtri(itri)==1) cycle
			effdistmax=-1D50
			surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
			surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
			surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0

			do iatm=1,ncenter
				disttmp=dsqrt((surtrix-a(iatm)%x)**2+(surtriy-a(iatm)%y)**2+(surtriz-a(iatm)%z)**2)
		 		if (imolsurparmode==1) then
		 			vdwrtmp=1D0 !Use original Voronoi partition
		 		else
					vdwrtmp=dsqrt(vdwr(a(iatm)%index))
				end if
				effdisttmp=1-disttmp/vdwrtmp
				if (effdisttmp>effdistmax) then
					surtrifrag(itri)=iatm
					effdistmax=effdisttmp
				end if
			end do
			if (abs(isel)==12) surtrifrag(itri)=ifatmfrag(surtrifrag(itri)) ! Then surtrifrag(itri)=1/0 means this itri belongs / doesn't belong to user-defined fragment
		end do
! 		do ifrag=1,nsurfrag	!Show the number of constituent facets for each fragment
! 			write(*,*) ifrag,count(surtrifrag==ifrag)
! 		end do
		!Obtain atomic or fragmental values
		fragsurarea=0D0
		fragsuravg=0D0
		fragsurvar=0D0
		fragsurmin=1D50
		fragsurmax=-1D50
		do itri=1,nsurtri !Get area and average
			if (elimtri(itri)==1) cycle
			iattfrag=surtrifrag(itri)
			if (iattfrag==0) cycle
			tmpval=surtriang(itri)%value
			tmpvalmin=minval(survtx(surtriang(itri)%idx(1:3))%value) !Find minimal and maximal value from its three vertices
			tmpvalmax=maxval(survtx(surtriang(itri)%idx(1:3))%value)
			fragsurarea(iattfrag,1)=fragsurarea(iattfrag,1)+surtriang(itri)%area
			fragsuravg(iattfrag,1)=fragsuravg(iattfrag,1)+surtriang(itri)%area*tmpval
			if (tmpval>0) then
				fragsurarea(iattfrag,2)=fragsurarea(iattfrag,2)+surtriang(itri)%area
				fragsuravg(iattfrag,2)=fragsuravg(iattfrag,2)+surtriang(itri)%area*tmpval
			else if (tmpval<0) then
				fragsurarea(iattfrag,3)=fragsurarea(iattfrag,3)+surtriang(itri)%area
				fragsuravg(iattfrag,3)=fragsuravg(iattfrag,3)+surtriang(itri)%area*tmpval
			end if
			if (tmpvalmax>fragsurmax(iattfrag)) fragsurmax(iattfrag)=tmpvalmax
			if (tmpvalmin<fragsurmin(iattfrag)) fragsurmin(iattfrag)=tmpvalmin
		end do
		fragsuravg=fragsuravg/fragsurarea
		fragsurchgsep=0D0
        fragpolaridx=0 !Only output for fragment, not every atom case
        fragareapol=0
        fragareanonpol=0
		do itri=1,nsurtri !Get variance
			if (elimtri(itri)==1) cycle
			iattfrag=surtrifrag(itri)
			if (iattfrag==0) cycle
			tmpval=surtriang(itri)%value
            tmparea=surtriang(itri)%area
			if (tmpval>0) then
				fragsurvar(iattfrag,2)=fragsurvar(iattfrag,2)+tmparea*(tmpval-fragsuravg(iattfrag,2))**2
			else if (tmpval<0) then
				fragsurvar(iattfrag,3)=fragsurvar(iattfrag,3)+tmparea*(tmpval-fragsuravg(iattfrag,3))**2
			end if
			fragsurchgsep(iattfrag)=fragsurchgsep(iattfrag)+tmparea*abs(tmpval-fragsuravg(iattfrag,1))
            fragpolaridx=fragpolaridx+tmparea*abs(tmpval)
            if (abs(tmpval)*au2kcal>10D0) then
                fragareapol=fragareapol+tmparea
            else
                fragareanonpol=fragareanonpol+tmparea
            end if
		end do
		fragsurvar(:,2:3)=fragsurvar(:,2:3)/fragsurarea(:,2:3)
 		do ifrag=1,nsurfrag
			if (fragsurarea(ifrag,2)>0) fragsurvar(ifrag,1)=fragsurvar(ifrag,1)+fragsurvar(ifrag,2)
			if (fragsurarea(ifrag,3)>0) fragsurvar(ifrag,1)=fragsurvar(ifrag,1)+fragsurvar(ifrag,3)
 		end do
		fragsurchgsep=fragsurchgsep/fragsurarea(:,1)
        fragpolaridx=fragpolaridx/fragsurarea(1,1)
		!Print atomic values
		write(*,*)
		if (isel==11) write(*,"(a,/)") " Note: The atoms having zero surface area (i.e. buried) are not shown below"
		if (abs(isel)==12) write(*,*) "Properties on the surface of this fragment:"
		if (imapfunc==1.or.imapfunc==3) then !ESP
			if (isel==11) then
				write(*,*) "Note: Minimal and maximal value below are in kcal/mol"
				write(*,*) " Atom#    All/Positive/Negative area (Ang^2)  Minimal value   Maximal value"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,1)==0D0) cycle
					write(*,"(i7,1x,3f12.5,2f16.8)") ifrag,fragsurarea(ifrag,:)*b2a*b2a,fragsurmin(ifrag)*au2kcal,fragsurmax(ifrag)*au2kcal
				end do
				write(*,*)
				write(*,*) "Note: Average and variance below are in kcal/mol and (kcal/mol)^2 respectively"
				write(*,*) " Atom#    All/Positive/Negative average       All/Positive/Negative variance"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,1)==0D0) cycle
					write(*,"(i7,1x,3f11.5,3x,3f11.5)") ifrag,fragsuravg(ifrag,:)*au2kcal,fragsurvar(ifrag,:)*au2kcal**2
				end do
				write(*,*)
				write(*,*) "Note: Internal charge separation (Pi) is in kcal/mol, nu = Balance of charges"
				write(*,*) " Atom#           Pi              nu         nu*sigma^2"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,1)==0D0) cycle
					balencechg=fragsurvar(ifrag,2)*fragsurvar(ifrag,3)/fragsurvar(ifrag,1)**2
					write(*,"(i7,1x,3f16.6)") ifrag,fragsurchgsep(ifrag)*au2kcal,balencechg,balencechg*fragsurvar(ifrag,1)*au2kcal**2
				end do
			else if (abs(isel)==12) then
				balencechg=fragsurvar(1,2)*fragsurvar(1,3)/fragsurvar(1,1)**2
				write(*,"(' Minimal value:',f13.6,' kcal/mol   Maximal value:',f13.6,' kcal/mol')") fragsurmin(1)*au2kcal,fragsurmax(1)*au2kcal
				write(*,"(' Overall surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,1),fragsurarea(1,1)*b2a*b2a
				write(*,"(' Positive surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
				write(*,"(' Negative surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,3),fragsurarea(1,3)*b2a*b2a
				write(*,"(' Overall average value: ',f13.8,' a.u. (',f14.8,' kcal/mol)')") fragsuravg(1,1),fragsuravg(1,1)*au2kcal
				write(*,"(' Positive average value:',f13.8,' a.u. (',f14.8,' kcal/mol)')") fragsuravg(1,2),fragsuravg(1,2)*au2kcal
				write(*,"(' Negative average value:',f13.8,' a.u. (',f14.8,' kcal/mol)')") fragsuravg(1,3),fragsuravg(1,3)*au2kcal
				write(*,"(' Overall variance (sigma^2_tot):',f12.8,' a.u.^2 (',f12.5,' (kcal/mol)^2)')") fragsurvar(1,1),fragsurvar(1,1)*au2kcal**2
				write(*,"(' Positive variance:     ',f13.8,' a.u.^2 (',f13.5,' (kcal/mol)^2)')") fragsurvar(1,2),fragsurvar(1,2)*au2kcal**2
				write(*,"(' Negative variance:     ',f13.8,' a.u.^2 (',f13.5,' (kcal/mol)^2)')") fragsurvar(1,3),fragsurvar(1,3)*au2kcal**2
				write(*,"(' Balance of charges (nu):',f13.8)") balencechg
				write(*,"(' Product of sigma^2_tot and nu: ',f12.8,' a.u.^2 (',f11.5,' (kcal/mol)^2)')") balencechg*fragsurvar(1,1),balencechg*fragsurvar(1,1)*au2kcal**2
				write(*,"(' Internal charge separation (Pi):',f13.8,' a.u. (',f13.5,' kcal/mol)')") fragsurchgsep(1),fragsurchgsep(1)*au2kcal
                write(*,"(' Molecular polarity index (MPI):',f13.8,' eV (',f13.5,' kcal/mol)')") fragpolaridx*au2eV,fragpolaridx*au2kcal
                write(*,"(' Nonpolar surface area (|ESP| <= 10 kcal/mol):',f10.2,' Angstrom^2  (',f6.2,' %)')") fragareanonpol*b2a*b2a,fragareanonpol/fragsurarea(1,1)*100
                write(*,"(' Polar surface area (|ESP| > 10 kcal/mol):    ',f10.2,' Angstrom^2  (',f6.2,' %)')") fragareapol*b2a*b2a,fragareapol/fragsurarea(1,1)*100
			end if
		else if (imapfunc==2) then !ALIE
			if (isel==11) then
				write(*,*) "Minimal, maximal and average value are in eV, variance is in eV^2"
				write(*,*) " Atom#      Area(Ang^2)  Min value   Max value       Average        Variance"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,2)==0D0) cycle
					write(*,"(i7,f16.5,2f12.6,2f15.6)") ifrag,fragsurarea(ifrag,2)*b2a*b2a,fragsurmin(ifrag)*au2eV,fragsurmax(ifrag)*au2eV,fragsuravg(ifrag,2)*au2eV,fragsurvar(ifrag,2)*au2eV**2
				end do
			else if (abs(isel)==12) then
				write(*,"(' Minimal value:',f13.6,' eV   Maximal value:',f13.6,' eV')") fragsurmin(1)*au2eV,fragsurmax(1)*au2eV
				write(*,"(' Surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
				write(*,"(' Average value: ',f13.8,' a.u. (',f13.8,' eV,',f14.8,' kcal/mol)')") fragsuravg(1,2),fragsuravg(1,2)*au2eV,fragsuravg(1,2)*au2kcal
				write(*,"(' Variance:  ',f13.8,' a.u.^2  (',f13.8,' eV^2,',E14.6,' kcal/mol^2)')") fragsurvar(1,2),fragsurvar(1,2)*au2eV**2,fragsurvar(1,2)*au2kcal**2
			end if
		else if (imapfunc==4.or.imapfunc==-4) then
			if (isel==11) then
				write(*,*) "Note: Below minimal and maximal values are in eV"
				write(*,*) " Atom#    All/Positive/Negative area (Ang^2)  Minimal value   Maximal value"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,1)==0D0) cycle
					write(*,"(i7,1x,3f12.5,2f16.8)") ifrag,fragsurarea(ifrag,:)*b2a*b2a,fragsurmin(ifrag)*au2eV,fragsurmax(ifrag)*au2eV
				end do
				write(*,*)
				write(*,*) "Note: Average and variance below are in eV and eV^2 respectively"
				write(*,*) " Atom#    All/Positive/Negative average       All/Positive/Negative variance"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,1)==0D0) cycle
					write(*,"(i7,1x,3f11.5,3x,3f11.5)") ifrag,fragsuravg(ifrag,:)*au2eV,fragsurvar(ifrag,:)*au2eV**2
				end do
				write(*,*)
			else if (abs(isel)==12) then
				write(*,"(' Minimal value:',f14.8,' eV,   Maximal value:',f14.8,' eV')") fragsurmin(1)*au2eV,fragsurmax(1)*au2eV
				write(*,"(' Overall surface area: ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,1),fragsurarea(1,1)*b2a*b2a
				write(*,"(' Positive surface area:',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
				write(*,"(' Negative surface area:',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,3),fragsurarea(1,3)*b2a*b2a
				write(*,"(' Overall average value: ',f13.8,' a.u. (',f14.8,' eV)')") fragsuravg(1,1),fragsuravg(1,1)*au2eV
				write(*,"(' Positive average value:',f13.8,' a.u. (',f14.8,' eV)')") fragsuravg(1,2),fragsuravg(1,2)*au2eV
				write(*,"(' Negative average value:',f13.8,' a.u. (',f14.8,' eV)')") fragsuravg(1,3),fragsuravg(1,3)*au2eV
				write(*,"(' Overall variance: ',f13.8,' a.u.^2 (',f14.8,' eV^2)')") fragsurvar(1,1),fragsurvar(1,1)*au2eV**2
				write(*,"(' Positive variance:',f13.8,' a.u.^2 (',f14.8,' eV^2)')") fragsurvar(1,2),fragsurvar(1,2)*au2eV**2
				write(*,"(' Negative variance:',f13.8,' a.u.^2 (',f14.8,' eV^2)')") fragsurvar(1,3),fragsurvar(1,3)*au2eV**2
			end if
		else !Other or unknown real space function provided by user
			if (isel==11) then
				write(*,*) " Atom#    All/Positive/Negative area (Ang^2)  Minimal value   Maximal value"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,1)==0D0) cycle
					write(*,"(i7,1x,3f12.5,2f16.9)") ifrag,fragsurarea(ifrag,:)*b2a*b2a,fragsurmin(ifrag),fragsurmax(ifrag)
				end do
				write(*,*)
				write(*,*) " Atom#   All/Positive/Negative average       All/Positive/Negative variance"
				do ifrag=1,nsurfrag
					if (fragsurarea(ifrag,1)==0D0) cycle
					write(*,"(i7,3(1PE13.5),3(1PE11.4))") ifrag,fragsuravg(ifrag,:),fragsurvar(ifrag,:)
				end do
			else if (abs(isel)==12) then
				write(*,"(' Minimal value:',f16.8,'    Maximal value:',f16.8)") fragsurmin(1),fragsurmax(1)
				write(*,"(' Overall surface area:      ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,1),fragsurarea(1,1)*b2a*b2a
				write(*,"(' Positive surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,2),fragsurarea(1,2)*b2a*b2a
				write(*,"(' Negative surface area:     ',f12.5,' Bohr^2  (',f10.5,' Angstrom^2)')") fragsurarea(1,3),fragsurarea(1,3)*b2a*b2a
				write(*,"(' Overall average value: ',f16.8)") fragsuravg(1,1)
				write(*,"(' Positive average value:',f16.8)") fragsuravg(1,2)
				write(*,"(' Negative average value:',f16.8)") fragsuravg(1,3)
				write(*,"(' Overall variance:  ',f16.8)") fragsurvar(1,1)
				write(*,"(' Positive variance: ',f16.8)") fragsurvar(1,2)
				write(*,"(' Negative variance: ',f16.8)") fragsurvar(1,3)
			end if
		end if
		! write(*,"(' Sum of area of total/pos./neg.:',3f12.6,' Bohr^2')") sum(fragsurarea(:,1)),sum(fragsurarea(:,2)),sum(fragsurarea(:,3))
		write(*,"(/,a)") " If outputting the surface facets to locsurf.pqr in current folder? By which you can visualize local surface via VMD program (y/n)"
		read(*,*) selectyn
		if (selectyn=='y'.or.selectyn=='Y') then
			open(10,file="locsurf.pqr",status="replace")
			write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface triangles')") nsurtri
			write(10,"(a)") "REMARK   The third last column is function values in a.u."
			if (ifPBC>0) then
				call getcellabc(asize,bsize,csize,alpha,beta,gamma)
				write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
			end if
			do itri=1,nsurtri
				if (elimtri(itri)==1) cycle
				surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
				surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
				surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
				write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
				"HETATM",itri,' '//"C "//' ',"MOL",'A',surtrifrag(itri),surtrix*b2a,surtriy*b2a,surtriz*b2a,surtriang(itri)%value,1.0," C"
			end do
			write(*,"(a)",advance="no") " Coordinate of the geometry center of the surface facets have been outputted to locsurf.pqr in current folder. &
            &The atomic charge column (i.e. the third last column) corresponds to function value in a.u."
			if (isel==11) write(*,"(a)") " Residue indices correspond to the attribution of the surface facets."
			if (abs(isel)==12) write(*,"(a)") " Residue index of 1/0 means that the corresponding surface facet belongs / does not belong to the fragment you defined."
			close(10)
		end if
		
	else if (isel==13) then !Calculate and export cube file of a mapped function
		write(*,*) "Calculating grid data of mapped function..."
		if (allocated(cubmattmp)) deallocate(cubmattmp)
        allocate(cubmattmp(nx,ny,nz))
        ifinish=0;ishowprog=1
        ntmp=floor(ny*nz/100D0)
		!$OMP PARALLEL DO SHARED(cubmattmp,ifinish,ishowprog) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
		do k=1,nz
			do j=1,ny
				do i=1,nx
                    call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					cubmattmp(i,j,k)=calcmapfunc(imapfunc,tmpx,tmpy,tmpz,nHirBecatm,HirBecatm)
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
		!$OMP END PARALLEL DO
        if (ishowprog/=0) call showprog(100,100)
		write(*,*) "Exporting mapfunc.cub..."
		open(10,file="mapfunc.cub",status="replace")
		call outcube(cubmattmp,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done! grid data has been exported to mapfunc.cub in current folder"
        deallocate(cubmattmp)
        
	else if (isel==14) then !Calculate area of the region around a specific surface extreme
        call extreme_area
        
    else if (isel==15) then !Basin-like partition of surface and calculate area
        call surfbasin
		
	else if (isel==16) then
		!Output facets via xyz file for debugging
! 		open(10,file="fac.xyz",status="replace")
! 		write(10,*) count(elimtri==0)
! 		write(10,*)
! 		do i=1,nsurtri
! 			if (elimtri(i)==1) cycle
			!surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
			!surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
			!surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
! 			write(10,"(a,3f14.8)") "O   ",surtrix*b2a,surtriy*b2a,surtriz*b2a
! 		end do
! 		close(10)
! 		write(*,*) "Center of surface facets have been outputted to fac.xyz in current folder"

		open(10,file="tri.pdb",status="replace")
		write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface triangles')") nsurtri
		if (ifPBC>0) then
			call getcellabc(asize,bsize,csize,alpha,beta,gamma)
			write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
		end if
		do itri=1,nsurtri
			if (elimtri(itri)/=0) cycle
			surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
			surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
			surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
			if (imapfunc==-1.or.imapfunc==0) tmpfuncval=surtriang(itri)%value
			if (imapfunc==1.or.imapfunc==3) tmpfuncval=surtriang(itri)%value*au2kcal
			if (imapfunc==2.or.imapfunc==4.or.imapfunc==-4) tmpfuncval=surtriang(itri)%value*au2eV
			if (tmpfuncval>999.99D0) tmpfuncval=999.99D0 !Avoid excess limit then become, because B-factor field only have three integer position
			if (tmpfuncval<-99.99D0) tmpfuncval=-99.99D0
			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
			"HETATM",itri,' '//"C "//' ',"MOL",'A',1,surtrix*b2a,surtriy*b2a,surtriz*b2a,1.0,tmpfuncval,"C "
		end do
		write(*,*) "Center of surface triangles have been outputted to tri.pdb in current folder"
		close(10)
        
	else if (isel==18) then  !Discard some surface extrema by inputting indices
        write(*,*) "Select the type of surface extrema"
        write(*,*) "1: Surface maxima   2: Surface minima"
        write(*,*) "Input 0 can return"
        read(*,*) itype
        if (itype==0) cycle
        write(*,*) "Input index range, e.g. 2,6-8,10"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nremove)
        allocate(tmpintarr(nremove))
        call str2arr(c2000tmp,nremove,tmpintarr)
        if (itype==1) surlocmaxidx(tmpintarr)=0
        if (itype==2) surlocminidx(tmpintarr)=0
        write(*,"(' Done!',i5,' extrema have been removed')") nremove
        deallocate(tmpintarr)
	
	else if (isel==19) then  !Merge some surface extrema and take their average position
        write(*,"(a)") " In this function, you need to select a set of surface extrema. The position of the first selected extreme &
        &will be replaced with average coordinate of all selected extrema, other extrema will be removed"
        write(*,*)
        write(*,*) "Select the type of surface extrema"
        write(*,*) "1: Surface maxima   2: Surface minima"
        write(*,*) "Input 0 can return"
        read(*,*) itype
        if (itype==0) cycle
        write(*,*) "Input index range, e.g. 2,6-8,10"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,nselsur)
        allocate(tmpintarr(nselsur))
        call str2arr(c2000tmp,nselsur,tmpintarr)
        avgx=0;avgy=0;avgz=0
        ncount=0 !The number of not discarded selected extrema
        ifirst=0 !The first not discarded selected extrema
        do itmp=1,nselsur
            isur=tmpintarr(itmp)
            if (itype==1) idx=surlocmaxidx(isur)
            if (itype==2) idx=surlocminidx(isur)
            if (idx==0) cycle !Already discarded
            avgx=avgx+survtx(idx)%x
            avgy=avgy+survtx(idx)%y
            avgz=avgz+survtx(idx)%z
            ncount=ncount+1
            if (ifirst==0) then
                ifirst=isur
                ifirstvtx=idx
            else
                if (itype==1) surlocmaxidx(isur)=0
                if (itype==2) surlocminidx(isur)=0
            end if
        end do
        survtx(ifirstvtx)%x=avgx/ncount
        survtx(ifirstvtx)%y=avgy/ncount
        survtx(ifirstvtx)%z=avgz/ncount
        deallocate(tmpintarr)
        if (itype==1) write(*,"(i5,' positive extrema have been merged as positive extreme',i5)") ncount,ifirst
        if (itype==2) write(*,"(i5,' negative extrema have been merged as negative extreme',i5)") ncount,ifirst
        write(*,*) "Note that the value of mapped function at this extreme is not updated"
	else if (isel==20) then !Fingerprint plot
        call fingerprt(HirBecatm,nHirBecatm,imapfunc)
	end if
end do

deallocate(elimvtx,elimtri,mergerelat)
deallocate(survtx,surtriang,vtxconn,vtxconnpos)
end do surfanaloop
end subroutine






!!----------------- Fingerprint plot analysis
subroutine fingerprt(HirBecatm,nHirBecatm,imapfunc)
use plot
use surfvertex
use util
use functions
implicit real*8 (a-h,o-z)
integer nHirBecatm,HirBecatm(nHirBecatm),tmparr(ncenter),ifcontactvtx(nsurvtx),imapfunc
real*8 dens(nsurvtx),densall(nsurvtx),elecontact(nelesupp,nelesupp)
real*8 d_i(nsurvtx),d_e(nsurvtx) !In Angstrom
real*8 :: rlow=0.6D0,rhigh=2.6D0,rstep=0.2D0
integer :: ptsize=15
integer,allocatable :: inarr(:),outarr(:),notHirBecatm(:)
real*8,allocatable :: xarr(:),yarr(:)
character c2000tmp*2000,c2tmp*2,c80tmp*80
!Set default inside and outside fragment
ninarr=nHirBecatm
nnotHirBecatm=ncenter-nHirBecatm !The number of atoms do not belong to Hirshfeld/Becke fragment
noutarr=nnotHirBecatm
allocate(inarr(ninarr),outarr(noutarr),notHirBecatm(nnotHirBecatm))
inarr=HirBecatm
itmp=0
do i=1,ncenter
	if (all(HirBecatm/=i)) then
		itmp=itmp+1
		notHirBecatm(itmp)=i
	end if
end do
outarr=notHirBecatm

do while(.true.)
    write(*,*)
    call menutitle("Fingerprint plot and local contact analyses",10,1)
	write(*,*) "-1 Return"
	write(*,*) "0 Start analysis"
	write(*,"(a,i8)") " 1 Set inside atoms for option 0, current number is",ninarr
	write(*,"(a,i8)") " 2 Set outside atoms for option 0, current number is",noutarr
    write(*,*) "3 Calculate contact area between different elements"
	read(*,*) isel

	if (isel==-1) then
		return
	else if (isel==1) then
		tmparr=0 !If tmparr(i)=1, i will be in this fragment
		write(*,*) "The indices of current inside atoms"
		write(*,"(10i7)") inarr
		write(*,"(a)") " Now input two conditions (atom index range and then element), the inside atoms will be their intersection"
		write(*,*)
		write(*,"(a)") " Input index range of atoms, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 are selected"
		write(*,"(a)") " Note: If press ENTER button directly, all atoms in the Hirshfeld/Becke fragment will be taken into consideration"
		read(*,"(a)") c2000tmp
		if (c2000tmp==" ") then
			nelement=nHirBecatm
			tmparr(1:nelement)=HirBecatm
		else
			call str2arr(c2000tmp,nelement)
			call str2arr(c2000tmp,nelement,tmparr(1:nelement))
		end if
		write(*,*) "Input element as filter condition, e.g. Cl"
		write(*,*) "Note: If press ENTER button directly, this condition will be ignored"
		read(*,"(a)") c2tmp
		deallocate(inarr)
		if (c2tmp==" ") then
			allocate(inarr(nelement))
			inarr=tmparr(1:nelement)
			ninarr=nelement
		else
			ninarr=count(a(tmparr(1:nelement))%name==c2tmp)
			allocate(inarr(ninarr))
			itmp=0
			do i=1,nelement
				if (a(tmparr(i))%name==c2tmp) then
					itmp=itmp+1
					inarr(itmp)=tmparr(i)
				end if
			end do
		end if
		write(*,*) "The inside atoms you chosen are"
		write(*,"(10i7)") inarr
		write(*,*)
		do i=1,ninarr
			if (all(HirBecatm/=inarr(i))) write(*,"(' Warning: atom',i7,' does not belong to the Hirshfeld/Becke fragment you previously set')") i
		end do
	else if (isel==2) then
		tmparr=0 !If tmparr(i)=1, i will be in this fragment
		write(*,*) "Index of current outside atoms"
		write(*,"(10i7)") outarr
		write(*,"(a)") " Now input two conditions (atom index range and then element), the inside atoms will be their intersection"
		write(*,*)
		write(*,"(a)") " Input index range of atoms, e.g. 1,3-6,8,10-11 means the atoms 1,3,4,5,6,8,10,11 are selected"
		write(*,"(a)") " Note: If press ENTER button directly, all atoms do not belong to the Hirshfeld/Becke fragment will be taken into consideration"
		read(*,"(a)") c2000tmp
		if (c2000tmp==" ") then
			nelement=nnotHirBecatm
			tmparr(1:nelement)=notHirBecatm
		else
			call str2arr(c2000tmp,nelement)
			call str2arr(c2000tmp,nelement,tmparr(1:nelement))
		end if
		write(*,*) "Input element as filter condition, e.g. Cl"
		write(*,*) "Note: If press ENTER button directly, this condition will be ignored"
		read(*,"(a)") c2tmp
		deallocate(outarr)
		if (c2tmp==" ") then
			allocate(outarr(nelement))
			outarr=tmparr(1:nelement)
			noutarr=nelement
		else
			noutarr=count(a(tmparr(1:nelement))%name==c2tmp)
			allocate(outarr(noutarr))
			itmp=0
			do i=1,nelement
				if (a(tmparr(i))%name==c2tmp) then
					itmp=itmp+1
					outarr(itmp)=tmparr(i)
				end if
			end do
		end if
		write(*,*) "The outside atoms you chosen are"
		write(*,"(10i7)") outarr
		write(*,*)
	else if (isel==0) then
		write(*,*) "Calculating points on contact surface..."
        !Finding vertices corresponding to contact surface between the two atom sets. If ivtx belongs to, then ifcontactvtx(ivtx)=1
        !Also determine d_i, d_e, d_norm of all vertices
		ifcontactvtx=0
		ncontactvtx=0 !Number of vertices on local contact surface
		ncurrvtx=0 !Number of non-eliminated vertices
		do ivtx=1,nsurvtx
			if (elimvtx(ivtx)==1) cycle
			ncurrvtx=ncurrvtx+1
			dist2minin=1D100
			dist2minout=1D100
			iminin=0
			iminout=0
			do iatm=1,ncenter
				dist2=(a(iatm)%x-survtx(ivtx)%x)**2+(a(iatm)%y-survtx(ivtx)%y)**2+(a(iatm)%z-survtx(ivtx)%z)**2
				if (any(HirBecatm==iatm)) then !Find the closest atom in Hirshfeld/Becke fragment to the surface point
					if (dist2<dist2minin) then
						dist2minin=dist2
						iminin=iatm
					end if
				else !Find the closest atom that does not belong to Hirshfeld/Becke fragment to the surface point
					if (dist2<dist2minout) then
						dist2minout=dist2
						iminout=iatm
					end if
				end if
			end do
			d_i(ivtx)=dsqrt(dist2minin)*b2a
			d_e(ivtx)=dsqrt(dist2minout)*b2a
			if (any(inarr==iminin).and.any(outarr==iminout)) then !This point is on local contact surface
				ncontactvtx=ncontactvtx+1
                ifcontactvtx(ivtx)=1 !This is a point on local contact surface
			end if
		end do
        
        !Calculate distribution density in fingerprint plot of every point on local contact surface
        !For each point, calculate number of other points in specific radius and divide by circle area
		radtest=0.03D0 !Compare radius in Angstrom
        dens=0
        densall=0
        !$OMP PARALLEL DO SHARED(dens,densall) PRIVATE(ivtx,di1,de1,jvtx,di2,de2,dist) schedule(dynamic) NUM_THREADS(nthreads)
		do ivtx=1,nsurvtx
			if (elimvtx(ivtx)==1) cycle
            di1=d_i(ivtx)
            de1=d_e(ivtx)
            !Statisticize entire contact surface
			do jvtx=1,nsurvtx
				if (elimvtx(jvtx)==1) cycle
				di2=d_i(jvtx)
				de2=d_e(jvtx)
                dist=dsqrt((di1-di2)**2+(de1-de2)**2)
                if (dist<radtest) densall(ivtx)=densall(ivtx)+1
            end do
            !Statisticize local contact surface (may also correspond to entire surface if user do not modify definition of inside/outside atoms)
			if (ifcontactvtx(ivtx)==0) cycle
			do jvtx=1,nsurvtx
				if (ifcontactvtx(jvtx)==0) cycle
				di2=d_i(jvtx)
				de2=d_e(jvtx)
                dist=dsqrt((di1-di2)**2+(de1-de2)**2)
                if (dist<radtest) dens(ivtx)=dens(ivtx)+1
			end do
        end do
        !$OMP END PARALLEL DO
        dens=dens/(pi*radtest**2)
        densmax=maxval(densall)/(pi*radtest**2) !Used to define upper limit of color scale
        
        !Loop over all center of surface facets, if a facet belongs to a contact surface, sum it area to arealoc
        !Also, calculate Hirshfeld surface area here
        arealoc=0
        areaall=0D0
        do itri=1,nsurtri
	        if (elimtri(itri)==1) cycle
			surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
			surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
			surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
			dist2minin=1D100
			dist2minout=1D100
			iminin=0
			iminout=0
			do iatm=1,ncenter
				dist2=(a(iatm)%x-surtrix)**2+(a(iatm)%y-surtriy)**2+(a(iatm)%z-surtriz)**2
				if (any(HirBecatm==iatm)) then !Find the closest atom in Hirshfeld/Becke fragment to the surface point
					if (dist2<dist2minin) then
						dist2minin=dist2
						iminin=iatm
					end if
				else !Find the closest atom that does not belong to Hirshfeld/Becke fragment to the surface point
					if (dist2<dist2minout) then
						dist2minout=dist2
						iminout=iatm
					end if
				end if
			end do
			if (any(inarr==iminin).and.any(outarr==iminout)) arealoc=arealoc+surtriang(itri)%area
	        areaall=areaall+surtriang(itri)%area
        end do
        arealoc=arealoc*b2a*b2a
        areaall=areaall*b2a*b2a
        
        write(*,*)
		write(*,"(' The area of the local contact surface is',f10.3,' Angstrom^2')") arealoc
		write(*,"(' The area of the total contact surface is',f10.3,' Angstrom^2')") areaall
		write(*,"(' The local surface occupies',f8.2,'% of the total surface')") arealoc/areaall*100D0
        write(*,*)
		write(*,"(' The number of points on the contact surface:        ',i10)") ncontactvtx
		write(*,"(' The number of points on total Hirshfeld/Becke surface:',i10)") ncurrvtx
		!write(*,"(' They occupy',f10.3,'% of total Hirshfeld/Becke surface points')") dfloat(ncontactvtx)/ncurrvtx*100D0
        
		do while(.true.)
		    write(*,*)
            call menutitle("Post process menu of fingerprint map and local contact analyses",6,1)
			write(*,*) "-1 Return"
			write(*,*) "0 Show fingerprint plot (quality is much poorer than that saved by option 1)"
			write(*,*) "1 Save fingerprint plot to an image file in current folder"
			write(*,"(a,i3)") " 2 Set size of points, current:",ptsize
			write(*,"(a,f8.3,a,f8.3,a,f5.2,a)") " 3 Set range of axes, current: from",rlow," to",rhigh,', with step',rstep,' Angstrom'
			write(*,*) "4 Export surface points to .pqr file in current folder"
			write(*,*) "5 Export d_i and d_e of surface points to .txt file in current folder"
			read(*,*) isel3
            
			if (isel3==-1) then
				exit
			else if (isel3==0.or.isel3==1) then
				call SCRMOD('REVERSE')
				CALL PAGE(3000,3000)
				if (isel3==0) then
					call METAFL('xwin')
					call window(100,100,800,800)
				else if (isel3==1) then
					call METAFL("pdf")
					call winsiz(graph1Dwidth,graph1Dwidth) !Ensure 1:1
				end if
				CALL setxid(0,'NONE')
				CALL DISINI
				if (isel3==0) call WINTIT("Click right mouse button to continue")
				call ERRMOD("ALL","OFF")
				call hwfont
				call center
                call AXSLEN(2400,2400)
				CALL HNAME(60) !Axis name size
				CALL height(65) !Axis tick size
				call TEXMOD("ON")
				CALL NAME('$d{_i}$ (Angstrom)','X')
				CALL NAME('$d{_e}$ (Angstrom)','Y')
				CALL LABDIG(1,"X")
				CALL LABDIG(1,"Y")
				CALL NAMDIS(70,"X")
				CALL NAMDIS(85,"Y")
                CALL TICKS(1,"XYZ")
				CALL GRAF(rlow,rhigh,rlow,rstep,rlow,rhigh,rlow,rstep)
                call setRGB(0.8D0,0.8D0,0.8D0) !Light gray
                CALL GRID(1,1)
				CALL INCMRK(-1) !Plot every point and do not connect them by line
				CALL MARKER(21) !Symbol of point
				CALL HSYMBL(ptsize) !Size of point
                !Plotting the surface points not on current contact surface
                nnow=ncurrvtx-ncontactvtx
                allocate(xarr(nnow),yarr(nnow))
                itmp=0
				do ivtx=1,nsurvtx
					if (elimvtx(ivtx)==1) cycle
                    if (ifcontactvtx(ivtx)==0) then
						itmp=itmp+1
						xarr(itmp)=d_i(ivtx)
						yarr(itmp)=d_e(ivtx)
                    end if
                end do
				call curve(xarr,yarr,nnow)
                deallocate(xarr,yarr)
                !Plot points on local contact surface
                do ivtx=1,nsurvtx
					if (ifcontactvtx(ivtx)==1) then
                        call percent2RGB(1,dens(ivtx)/densmax,Rval,Gval,Bval)
                        call setRGB(Rval,Gval,Bval)
						CALL MARKER(21)
						call curve(d_i(ivtx),d_e(ivtx),1)
                    end if
                end do
                !Finish
				call color("WHITE") !Restore to default (black)
				CALL DISFIN
				if (isel3==1) write(*,"(a)") " The image has been saved to a pdf file in current folder with DISLIN prefix"
            else if (isel3==2) then
				write(*,*) "Input size of points, e.g. 10"
				read(*,*) ptsize
			else if (isel3==3) then
				write(*,*) "Input lower, upper limits of axis in Angstrom, e.g. 0.4,2.8"
                read(*,*) rlow,rhigh
                write(*,*) "Input step size, e.g. 0.2"
                read(*,*) rstep
			else if (isel3==4) then
				open(10,file="finger.pqr",status="replace")
				write(10,"('REMARK   Generated by Multiwfn, totally',i10,' points on the surface')") ncontactvtx
				if (ifPBC>0) then
					call getcellabc(asize,bsize,csize,alpha,beta,gamma)
					write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
				end if
				do ivtx=1,nsurvtx
					if (ifcontactvtx(ivtx)==0) cycle
                    idx=ivtx
                    if (ivtx>99999) idx=99999
                    write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
					"HETATM",idx,' '//"C "//' ',"MOL",'A',1,survtx(ivtx)%x*b2a,survtx(ivtx)%y*b2a,survtx(ivtx)%z*b2a,survtx(ivtx)%value,1D0," C"
				end do
				write(10,"('END')")
				close(10)
				write(*,"(a)") " The points on the local contact surface have been outputted to finger.pqr in current folder"
				open(10,file="finger_all.pqr",status="replace")
				write(10,"('REMARK   Generated by Multiwfn, totally',i10,' points on the surface')") ncurrvtx
				if (ifPBC>0) then
					call getcellabc(asize,bsize,csize,alpha,beta,gamma)
					write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
				end if
				do ivtx=1,nsurvtx
					if (elimvtx(ivtx)==1) cycle
                    idx=ivtx
                    if (ivtx>99999) idx=99999
					write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f9.4,a2)") &
					"HETATM",idx,' '//"C "//' ',"MOL",'A',1,survtx(ivtx)%x*b2a,survtx(ivtx)%y*b2a,survtx(ivtx)%z*b2a,survtx(ivtx)%value,1D0," C"
				end do
				write(10,"('END')")
				close(10)
				write(*,"(a)") " The points on the entire contact surface have been outputted to finger_all.pqr in current folder"
                if (imapfunc==11) then
					write(*,*) """Charge"" field in these files correspond to electron density (a.u.)"
                else if (imapfunc==22) then
					write(*,*) """Charge"" field in these files correspond to d_norm in Angstrom"
                else
					write(*,*) """Charge"" field in these files correspond to the mapped function value"
                end if
			else if (isel3==5) then
				open(10,file="di_de.txt",status="replace")
				do ivtx=1,nsurvtx
					if (ifcontactvtx(ivtx)==0) cycle
					write(10,"(2f12.6)") d_i(ivtx),d_e(ivtx)
				end do
				close(10)
				write(*,"(a)") " d_i and d_e of the points on the local contact surface been exported to di_de.txt in current folder"
				open(10,file="di_de_all.txt",status="replace")
				do ivtx=1,nsurvtx
					if (elimvtx(ivtx)==1) cycle
					write(10,"(2f12.6)") d_i(ivtx),d_e(ivtx)
				end do
				close(10)
				write(*,"(a)") " d_i and d_e of the points on the entire Hirshfeld/Becke surface been exported to di_de_all.txt in current folder"
                write(*,*) "Columns 1 and 2 in these files are d_i and d_e, respectively"
			end if
		end do

    else if (isel==3) then !Calculate contact area between different elements 
        !Loop over all center of surface facets, attribute its area to element pair
		elecontact=0 !elecontact(8,1) means contact area between oxygen of inside atoms and hydrogen of outside atoms
        do itri=1,nsurtri
	        if (elimtri(itri)==1) cycle
			surtrix=sum(survtx(surtriang(itri)%idx(1:3))%x)/3D0 !Center of triangles
			surtriy=sum(survtx(surtriang(itri)%idx(1:3))%y)/3D0
			surtriz=sum(survtx(surtriang(itri)%idx(1:3))%z)/3D0
			dist2minin=1D100
			dist2minout=1D100
			iminin=0
			iminout=0
			do iatm=1,ncenter
				dist2=(a(iatm)%x-surtrix)**2+(a(iatm)%y-surtriy)**2+(a(iatm)%z-surtriz)**2
				if (any(HirBecatm==iatm)) then !Find the closest atom in Hirshfeld/Becke fragment to the surface point
					if (dist2<dist2minin) then
						dist2minin=dist2
						iminin=iatm
					end if
				else !Find the closest atom that does not belong to Hirshfeld/Becke fragment to the surface point
					if (dist2<dist2minout) then
						dist2minout=dist2
						iminout=iatm
					end if
				end if
			end do
            iele=a(iminin)%index
            jele=a(iminout)%index
			elecontact(iele,jele)=elecontact(iele,jele)+surtriang(itri)%area
        end do
        elecontact=elecontact*b2a*b2a
        areaall=sum(elecontact)
        write(*,"(a)") " Inside element, outside element, their contact area (Angstrom^2) and percentage (%)"
        do iele=1,nelesupp !Loop inside elements
			do jele=1,nelesupp !Loop outside elements
				if (elecontact(iele,jele)/=0) then
					write(c80tmp,"(a,'-',a)") trim(ind2name(iele)),trim(ind2name(jele))
					write(*,"(2x,a,f12.3,f12.3)") c80tmp(1:5),elecontact(iele,jele),elecontact(iele,jele)/areaall*100
                end if
            end do
        end do
        write(*,*)
        write(*,*) "The same as above, but do not distinguish inside and outside elements"
        do iele=1,nelesupp
			do jele=iele,nelesupp
				if (iele==jele) then
					if (elecontact(iele,iele)/=0) then
						write(c80tmp,"(a,'-',a)") trim(ind2name(iele)),trim(ind2name(jele))
						write(*,"(2x,a,f12.3,f12.3)") c80tmp(1:11),elecontact(iele,iele),elecontact(iele,iele)/areaall*100
					end if
                else
					tmparea=elecontact(iele,jele)+elecontact(jele,iele)
					if (tmparea/=0) then
						write(c80tmp,"(a,'-',a,'/',a,'-',a)") trim(ind2name(iele)),trim(ind2name(jele)),trim(ind2name(jele)),trim(ind2name(iele))
						write(*,"(2x,a,f12.3,f12.3)") c80tmp(1:11),tmparea,tmparea/areaall*100
					end if
                end if
            end do
        end do
        write(*,*)
		write(*,"(' Area of total contact surface is',f10.3,' Angstrom^2')") areaall
	end if
end do
end subroutine



!!-------- Use Marching Tetrahedra algorithem, decompose cube to several tetrahedra
subroutine marchtetra(ix,iy,iz)
implicit real*8 (a-h,o-z)
integer ix,iy,iz
!    7-------8
!   /|      /|
!  6-+-----4 |
!  | |     | |
!  | 3-----+-1
!  |/      |/
!  5-------2
!
!   Z
!   |
!   0---Y    
!  / 
! X

! Five tetrahedra, may lead too big spacing somewhere
! call genvertex(ix,iy,iz,1,2,5,4) 
! call genvertex(ix,iy,iz,4,5,6,7)
! call genvertex(ix,iy,iz,1,3,5,7)
! call genvertex(ix,iy,iz,1,4,7,8)
! call genvertex(ix,iy,iz,1,4,7,5)

! WFA original paper
! call genvertex(ix,iy,iz,3,2,1,8) 
! call genvertex(ix,iy,iz,3,2,4,8)
! call genvertex(ix,iy,iz,3,7,4,8)
! call genvertex(ix,iy,iz,3,4,5,2)
! call genvertex(ix,iy,iz,3,4,5,7)
! call genvertex(ix,iy,iz,6,4,5,7)

! Main-axis decomposition, all tetrahedra share 4-3 axis, as in http://paulbourke.net/geometry/polygonise/
call genvertex(ix,iy,iz,4,3,5,2) 
call genvertex(ix,iy,iz,4,3,5,6)
call genvertex(ix,iy,iz,4,3,7,6)
call genvertex(ix,iy,iz,4,3,7,8)
call genvertex(ix,iy,iz,4,3,1,8)
call genvertex(ix,iy,iz,4,3,1,2)
end subroutine

!!-------- Use marching cube algorithem, interpolate each edge. Unfinished routine, can only generated surface vertices but not connecitivity
! subroutine marchcube(ix,iy,iz)
! use surfvertex
! use defvar
! implicit real*8 (a-h,o-z)
! integer ix,iy,iz
! itestc1=0
! itestc2=0
! itestc3=0
! itestc4=0
! itestc5=0
! itestc6=0
! itestc7=0
! itestc8=0
! !Test function value at each corner of cube, =0/1 means lower/higher than isovalue
! if (cubmat(ix,iy+1,iz)>=surfisoval) itestc1=1
! if (cubmat(ix+1,iy+1,iz)>=surfisoval) itestc2=1
! if (cubmat(ix,iy,iz)>=surfisoval) itestc3=1
! if (cubmat(ix+1,iy+1,iz+1)>=surfisoval) itestc4=1
! if (cubmat(ix+1,iy,iz)>=surfisoval) itestc5=1
! if (cubmat(ix+1,iy,iz+1)>=surfisoval) itestc6=1
! if (cubmat(ix,iy,iz+1)>=surfisoval) itestc7=1
! if (cubmat(ix,iy+1,iz+1)>=surfisoval) itestc8=1
! !interpolate each cube edge
! ! 1-2,2-5,5-3,3-1
! if (itestc1+itestc2==1) call vertexinterpolate(ix,iy+1,iz,ix+1,iy+1,iz)
! if (itestc2+itestc5==1) call vertexinterpolate(ix+1,iy+1,iz,ix+1,iy,iz)
! if (itestc5+itestc3==1) call vertexinterpolate(ix+1,iy,iz,ix,iy,iz)
! if (itestc3+itestc1==1) call vertexinterpolate(ix,iy,iz,ix,iy+1,iz)
! ! 1-8,2-4,5-6,3-7
! if (itestc1+itestc8==1) call vertexinterpolate(ix,iy+1,iz,ix,iy+1,iz+1)
! if (itestc2+itestc4==1) call vertexinterpolate(ix+1,iy+1,iz,ix+1,iy+1,iz+1)
! if (itestc5+itestc6==1) call vertexinterpolate(ix+1,iy,iz,ix+1,iy,iz+1)
! if (itestc3+itestc7==1) call vertexinterpolate(ix,iy,iz,ix,iy,iz+1)
! ! 8-4,4-6,6-7,7-8
! if (itestc8+itestc4==1) call vertexinterpolate(ix,iy+1,iz+1,ix+1,iy+1,iz+1)
! if (itestc4+itestc6==1) call vertexinterpolate(ix+1,iy+1,iz+1,ix+1,iy,iz+1)
! if (itestc6+itestc7==1) call vertexinterpolate(ix+1,iy,iz+1,ix,iy,iz+1)
! if (itestc7+itestc8==1) call vertexinterpolate(ix,iy,iz+1,ix,iy+1,iz+1)
! end subroutine


!!---------- Generate surface vertices from tetrahedra. inum is the vertex index within each tetrahedron 
subroutine genvertex(baseix,baseiy,baseiz,inum1,inum2,inum3,inum4)
use util
use defvar
use surfvertex
implicit real*8 (a-h,o-z)
integer baseix,baseiy,baseiz,inum1,inum2,inum3,inum4
integer itestv(4),vix(4),viy(4),viz(4),newvtxind(4) !1~4 is the four vertices index within current tetrahedron. newvtxind is absolute surface vertex index
real*8 vtxval(4)
!Convert index within cube to absolute cubic index
call getvertind(vix(1),viy(1),viz(1),baseix,baseiy,baseiz,inum1)
call getvertind(vix(2),viy(2),viz(2),baseix,baseiy,baseiz,inum2)
call getvertind(vix(3),viy(3),viz(3),baseix,baseiy,baseiz,inum3)
call getvertind(vix(4),viy(4),viz(4),baseix,baseiy,baseiz,inum4)
!Values of vertices in current tetrahedron
vtxval(1)=cubmat(vix(1),viy(1),viz(1))
vtxval(2)=cubmat(vix(2),viy(2),viz(2))
vtxval(3)=cubmat(vix(3),viy(3),viz(3))
vtxval(4)=cubmat(vix(4),viy(4),viz(4))
!Test function value at each vertex of tetrahedron, =0/1 means lower/higher than isovalue
itestv=0
where (vtxval>surfisoval) itestv=1 !If the tetrahedron vertex is inside isosurface, then its itestv=1, else 0
!    2
!    |
!    |
!    1
!   / \
!  /   \
! 3     4
!Interpolate surface vertices, if it is new, then accumulate it to "survtx" array
itesttot=sum(itestv)
newvtxind=0
inewvtxind=0
if (itesttot==0) then ! external tetrahedron, return directly
	return
else if (itesttot==4) then !Internal tetrahedron, calculate its volume and then return
	pax=orgx+(vix(1)-1)*dx
	pay=orgy+(viy(1)-1)*dy
	paz=orgz+(viz(1)-1)*dz
	pbx=orgx+(vix(2)-1)*dx
	pby=orgy+(viy(2)-1)*dy
	pbz=orgz+(viz(2)-1)*dz
	pcx=orgx+(vix(3)-1)*dx
	pcy=orgy+(viy(3)-1)*dy
	pcz=orgz+(viz(3)-1)*dz
	pdx=orgx+(vix(4)-1)*dx
	pdy=orgy+(viy(4)-1)*dy
	pdz=orgz+(viz(4)-1)*dz
	tetravol0=tetravol0+gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz) !Volume type 0
else if (itesttot==1) then
	do iint=1,4 !Find the only internal tetrahedron vertex
		if (itestv(iint)==1) exit
	end do
	do iplt=1,4
		if (iplt==iint) cycle
		inewvtxind=inewvtxind+1
		call vertexinterpolate(vix(iint),viy(iint),viz(iint),vix(iplt),viy(iplt),viz(iplt),newvtxind(inewvtxind))
	end do
	call addvtxconn(newvtxind(1),newvtxind(2))
	call addvtxconn(newvtxind(1),newvtxind(3))
	call addvtxconn(newvtxind(2),newvtxind(3))
	nsurtri=nsurtri+1
	surtriang(nsurtri)%idx(1:3)=newvtxind(1:3)
	!Calculate volume of newly generated tetrahedron
	
	pax=orgx+(vix(iint)-1)*dx
	pay=orgy+(viy(iint)-1)*dy
	paz=orgz+(viz(iint)-1)*dz
	pbx=survtx(newvtxind(1))%x
	pby=survtx(newvtxind(1))%y
	pbz=survtx(newvtxind(1))%z
	pcx=survtx(newvtxind(2))%x
	pcy=survtx(newvtxind(2))%y
	pcz=survtx(newvtxind(2))%z
	pdx=survtx(newvtxind(3))%x
	pdy=survtx(newvtxind(3))%y
	pdz=survtx(newvtxind(3))%z
	tetravol1=tetravol1+gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz) !Volume type 1
else if (itesttot==3) then
	do iext=1,4 !Find the only external tetrahedron vertex
		if (itestv(iext)==0) exit
	end do
	do iplt=1,4
		if (iplt==iext) cycle
		inewvtxind=inewvtxind+1
		call vertexinterpolate(vix(iext),viy(iext),viz(iext),vix(iplt),viy(iplt),viz(iplt),newvtxind(inewvtxind))
	end do
	call addvtxconn(newvtxind(1),newvtxind(2))
	call addvtxconn(newvtxind(1),newvtxind(3))
	call addvtxconn(newvtxind(2),newvtxind(3))
	nsurtri=nsurtri+1
	surtriang(nsurtri)%idx(1:3)=newvtxind(1:3)
	!Calculate volume of newly generated tetrahedron
	pax=orgx+(vix(iext)-1)*dx
	pay=orgy+(viy(iext)-1)*dy
	paz=orgz+(viz(iext)-1)*dz
	pbx=survtx(newvtxind(1))%x
	pby=survtx(newvtxind(1))%y
	pbz=survtx(newvtxind(1))%z
	pcx=survtx(newvtxind(2))%x
	pcy=survtx(newvtxind(2))%y
	pcz=survtx(newvtxind(2))%z
	pdx=survtx(newvtxind(3))%x
	pdy=survtx(newvtxind(3))%y
	pdz=survtx(newvtxind(3))%z
	voidvol=gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz)
	!Calculate the entire volume of original tetrahedron
	
	pax=orgx+(vix(1)-1)*dx
	pay=orgy+(viy(1)-1)*dy
	paz=orgz+(viz(1)-1)*dz
	pbx=orgx+(vix(2)-1)*dx
	pby=orgy+(viy(2)-1)*dy
	pbz=orgz+(viz(2)-1)*dz
	pcx=orgx+(vix(3)-1)*dx
	pcy=orgy+(viy(3)-1)*dy
	pcz=orgz+(viz(3)-1)*dz
	pdx=orgx+(vix(4)-1)*dx
	pdy=orgy+(viy(4)-1)*dy
	pdz=orgz+(viz(4)-1)*dz
	wholetetravol=gettetravol(pax,pay,paz,pbx,pby,pbz,pcx,pcy,pcz,pdx,pdy,pdz)
	tetravol2=tetravol2+(wholetetravol-voidvol) !Volume type 2
else if (itesttot==2) then
	itv1=0
	itv2=0
	do iint=1,4
		if (itestv(iint)==1) then  !Find the two internal tetrahedron vertices
			if (itv1==0) then
				itv1=iint
			else
				itv2=iint
			end if
			do iother=1,4 !Find external vertices (namely itestv=0)
				if (itestv(iother)==1) cycle !If this is also interal one, skip it
				inewvtxind=inewvtxind+1
				call vertexinterpolate(vix(iint),viy(iint),viz(iint),vix(iother),viy(iother),viz(iother),newvtxind(inewvtxind))
			end do
		end if
	end do
	! 1--2
	! |  |
	! 3--4
	! 1-2-4 is the first triangle, 1-3-4 is the second one
	call addvtxconn(newvtxind(1),newvtxind(2))
	call addvtxconn(newvtxind(1),newvtxind(3))
	call addvtxconn(newvtxind(2),newvtxind(4))
	call addvtxconn(newvtxind(3),newvtxind(4))
	call addvtxconn(newvtxind(1),newvtxind(4))
	nsurtri=nsurtri+1
	surtriang(nsurtri)%idx(1:2)=newvtxind(1:2)
	surtriang(nsurtri)%idx(3)=newvtxind(4)
	nsurtri=nsurtri+1
	surtriang(nsurtri)%idx(1)=newvtxind(1)
	surtriang(nsurtri)%idx(2:3)=newvtxind(3:4)
	!Now split the widget with six vertices to three tetrahedra, so that volume can be calculated
	!decompose to tv1-sv1-sv3-sv4,tv1-tv2-sv1-sv4,tv2-sv1-sv2-sv4  (tv: tetrahedron vertex with value larger than isovalue;  sv: square vertex)
!
! 	      tv1-tv2
! 	
!      sv1-------sv2
!      /         /
!    sv3-------sv4

	tv1x=orgx+(vix(itv1)-1)*dx
	tv1y=orgy+(viy(itv1)-1)*dy
	tv1z=orgz+(viz(itv1)-1)*dz
	tv2x=orgx+(vix(itv2)-1)*dx
	tv2y=orgy+(viy(itv2)-1)*dy
	tv2z=orgz+(viz(itv2)-1)*dz
	sv1x=survtx(newvtxind(1))%x
	sv1y=survtx(newvtxind(1))%y
	sv1z=survtx(newvtxind(1))%z
	sv2x=survtx(newvtxind(2))%x
	sv2y=survtx(newvtxind(2))%y
	sv2z=survtx(newvtxind(2))%z
	sv3x=survtx(newvtxind(3))%x
	sv3y=survtx(newvtxind(3))%y
	sv3z=survtx(newvtxind(3))%z
	sv4x=survtx(newvtxind(4))%x
	sv4y=survtx(newvtxind(4))%y
	sv4z=survtx(newvtxind(4))%z
	tetravol3=tetravol3+gettetravol(tv1x,tv1y,tv1z,sv1x,sv1y,sv1z,sv3x,sv3y,sv3z,sv4x,sv4y,sv4z) !Volume type 3
	tetravol3=tetravol3+gettetravol(tv1x,tv1y,tv1z,tv2x,tv2y,tv2z,sv1x,sv1y,sv1z,sv4x,sv4y,sv4z)
	tetravol3=tetravol3+gettetravol(tv2x,tv2y,tv2z,sv1x,sv1y,sv1z,sv2x,sv2y,sv2z,sv4x,sv4y,sv4z)
end if
	
! if (itestv(1)+itestv(2)==1) call vertexinterpolate(vix(1),viy(1),viz(1),vix(2),viy(2),viz(2),newvtxindt)
! if (itestv(1)+itestv(3)==1) call vertexinterpolate(vix(1),viy(1),viz(1),vix(3),viy(3),viz(3),newvtxindt)
! if (itestv(1)+itestv(4)==1) call vertexinterpolate(vix(1),viy(1),viz(1),vix(4),viy(4),viz(4),newvtxindt)
! if (itestv(2)+itestv(3)==1) call vertexinterpolate(vix(2),viy(2),viz(2),vix(3),viy(3),viz(3),newvtxindt)
! if (itestv(2)+itestv(4)==1) call vertexinterpolate(vix(2),viy(2),viz(2),vix(4),viy(4),viz(4),newvtxindt)
! if (itestv(3)+itestv(4)==1) call vertexinterpolate(vix(3),viy(3),viz(3),vix(4),viy(4),viz(4),newvtxindt)
end subroutine

!!--------- Convert numbering of each vertex of tetrahedron to absolute corner-index
subroutine getvertind(indx,indy,indz,baseix,baseiy,baseiz,inum)
implicit real*8 (a-h,o-z)
integer indx,indy,indz,baseix,baseiy,baseiz,inum
indx=baseix
indy=baseiy
indz=baseiz
if (inum==1) then
	indy=baseiy+1
else if (inum==2) then
	indx=baseix+1
	indy=baseiy+1
else if (inum==3) then
	continue
else if (inum==4) then
	indx=baseix+1
	indy=baseiy+1
	indz=baseiz+1
else if (inum==5) then
	indx=baseix+1
else if (inum==6) then
	indx=baseix+1
	indz=baseiz+1
else if (inum==7) then
	indz=baseiz+1
else if (inum==8) then
	indy=baseiy+1
	indz=baseiz+1
end if
end subroutine

!!--------- Interpolate surface vertex from input index of vertex, then save to "survtx" array
!newind is the index of this surface vertex
!ifuncint define use which real space function to do the interpolation
subroutine vertexinterpolate(iax,iay,iaz,ibx,iby,ibz,newind)
use defvar
use surfvertex
use functions
implicit real*8 (a-h,o-z)
integer iax,iay,iaz,ibx,iby,ibz,newind
icora=abs2suridx(iax,iay,iaz)
icorb=abs2suridx(ibx,iby,ibz)
do i=1,surcor2vtxpos(icora)
	if (surfcor2vtx(icora,i)%athcor==icorb) then !Have already interpolated for these two corner, return interpolated vertex index directly
		newind=surfcor2vtx(icora,i)%itpvtx
		return
	end if
end do

nsurvtx=nsurvtx+1
newind=nsurvtx
surcor2vtxpos(icora)=surcor2vtxpos(icora)+1
surfcor2vtx(icora,surcor2vtxpos(icora))%athcor=icorb
surfcor2vtx(icora,surcor2vtxpos(icora))%itpvtx=nsurvtx
surcor2vtxpos(icorb)=surcor2vtxpos(icorb)+1
surfcor2vtx(icorb,surcor2vtxpos(icorb))%athcor=icora
surfcor2vtx(icorb,surcor2vtxpos(icorb))%itpvtx=nsurvtx

vala=cubmat(iax,iay,iaz)
valb=cubmat(ibx,iby,ibz)

if (nbisec==0) then !Linear interpolate directly
	scl=(vala-surfisoval)/(vala-valb)
	survtx(nsurvtx)%x=(1-scl)*(orgx+(iax-1)*dx)+scl*(orgx+(ibx-1)*dx)
	survtx(nsurvtx)%y=(1-scl)*(orgy+(iay-1)*dy)+scl*(orgy+(iby-1)*dy)
	survtx(nsurvtx)%z=(1-scl)*(orgz+(iaz-1)*dz)+scl*(orgz+(ibz-1)*dz)
else !First perform bisect several times
	ax=orgx+(iax-1)*dx
	ay=orgy+(iay-1)*dy
	az=orgz+(iaz-1)*dz
	bx=orgx+(ibx-1)*dx
	by=orgy+(iby-1)*dy
	bz=orgz+(ibz-1)*dz
	do i=1,nbisec
		bisecx=(ax+bx)/2D0
		bisecy=(ay+by)/2D0
		bisecz=(az+bz)/2D0
		bisecval=calcfuncall(ifuncintp,bisecx,bisecy,bisecz)
		if ((bisecval-surfisoval)*(vala-surfisoval)<0) then
			bx=bisecx
			by=bisecy
			bz=bisecz
			valb=bisecval
		else
			ax=bisecx
			ay=bisecy
			az=bisecz
			vala=bisecval
		end if
	end do
	!Last time, use linear interpolation
	if ((bisecval-surfisoval)*(vala-surfisoval)<0) then !interpolate for a and bisec
		scl=(vala-surfisoval)/(vala-bisecval)
		survtx(nsurvtx)%x=(1-scl)*ax+scl*bisecx
		survtx(nsurvtx)%y=(1-scl)*ay+scl*bisecy
		survtx(nsurvtx)%z=(1-scl)*az+scl*bisecz
	else !interpolate for b and bisec
		scl=(valb-surfisoval)/(valb-bisecval)
		survtx(nsurvtx)%x=(1-scl)*bx+scl*bisecx
		survtx(nsurvtx)%y=(1-scl)*by+scl*bisecy
		survtx(nsurvtx)%z=(1-scl)*bz+scl*bisecz
	end if
end if
end subroutine


!!----------Add a connection between two surface vertices
subroutine addvtxconn(ind1,ind2)
use surfvertex
implicit real*8 (a-h,o-z)
integer ind1,ind2
do i=1,vtxconnpos(ind1)
	if (vtxconn(ind1,i)==ind2) return !Have already existed this connection entry
end do
vtxconnpos(ind1)=vtxconnpos(ind1)+1
vtxconn(ind1,vtxconnpos(ind1))=ind2
vtxconnpos(ind2)=vtxconnpos(ind2)+1
vtxconn(ind2,vtxconnpos(ind2))=ind1
end subroutine


!!---------- Calculate mapped function value
real*8 function calcmapfunc(imapfunc,x,y,z,nHirBecatm,HirBecatm)
use defvar
use functions
integer imapfunc,nHirBecatm,HirBecatm(nHirBecatm)
real*8 x,y,z
if (imapfunc==1) then
	calcmapfunc=totesp(x,y,z)
else if (imapfunc==2) then
	calcmapfunc=avglocion(x,y,z)
else if (imapfunc==3) then
	calcmapfunc=nucesp(x,y,z)
else if (imapfunc==4) then
	calcmapfunc=loceleaff(x,y,z)
else if (imapfunc==-4) then
	calcmapfunc=loceleatt(x,y,z)
else if (imapfunc==5) then
	calcmapfunc=edr(x,y,z)
else if (imapfunc==6) then
	calcmapfunc=edrdmax(x,y,z)
else if (imapfunc==-1) then
	calcmapfunc=userfunc(x,y,z)
else if (imapfunc==10) then
	calcmapfunc=pairfunc(refx,refy,refz,x,y,z)  !Calculate pair density, not for normal users
else if (imapfunc==11) then  !Calculate electron density
	if (allocated(b)) then
		calcmapfunc=fdens(x,y,z)
	else !with promolecular approximation
		calcmapfunc=calcprodens(x,y,z,0)
	end if
else if (imapfunc==12) then  !Calculate Sign(lambda2)*rho
	if (allocated(b)) then
		calcmapfunc=signlambda2rho(x,y,z)
	else !with promolecular approximation
		calcmapfunc=signlambda2rho_prodens(x,y,z)
	end if
else if (imapfunc==20) then
	calcmapfunc=surfana_di(x,y,z,nHirBecatm,HirBecatm)
else if (imapfunc==21) then
	calcmapfunc=surfana_de(x,y,z,nHirBecatm,HirBecatm)
else if (imapfunc==22) then
	calcmapfunc=surfana_norm(x,y,z,nHirBecatm,HirBecatm)
end if
end function





!!------ Calculate area of the region around a specific surface extreme
subroutine extreme_area
use defvar
use surfvertex
implicit real*8 (a-h,o-z)
integer iselvtx(nsurvtx) !=1 means this vertex belongs to selected region, =0 means not

write(*,*) "Select type of surface extreme, 1=minimum  2=maximum"
read(*,*) iexttype
write(*,*) "Input index of the surface extreme, e.g. 3"
read(*,*) iext

if (iexttype==1) then
    idx=surlocminidx(iext)
    valext=survtx(idx)%value
    write(*,"(' Value of this surface minimum:',f12.6,' a.u.',/)") valext
    write(*,*) "Input criterion of determining the region, e.g. 0.05"
    write(*,"(a)") " Note: If a vertex is directly or indirectly connected to the selected surface minimum and &
    &its value is lower than this criterion, then this vertex will belong to the region"
else
    idx=surlocmaxidx(iext)
    valext=survtx(idx)%value
    write(*,"(' Value of this surface maximum:',f12.6,' a.u.',/)") valext
    write(*,*) "Input criterion of determining the region, e.g. 0.05"
    write(*,"(a)") " Note: If a vertex is directly or indirectly connected to the selected surface maximum and &
    &its value is higher than this criterion, then this vertex will belong to the region"
end if
read(*,*) critval

nselvtx=1 !The number of selected vertices
iselvtx(:)=0
iselvtx(idx)=1 !Initially, only the vertex corresponding to extreme is selected
do while(.true.)
    nselvtx_old=nselvtx
    !Cycle all surface vertex, if it is connected to a vertex that direct/indirect contact with &
    !the selected extreme and satisfies criterion, then this vertex will also be labelled as "contact"
    !The loop is infinite, until no further update can be realized
    do ivtx=1,nsurvtx
		if (elimvtx(ivtx)==1) cycle
        vtxval=survtx(ivtx)%value
        if ( (iexttype==1.and.vtxval>critval).or.((iexttype==2.and.vtxval<critval)) ) cycle
        do jtmp=1,vtxconnpos(ivtx) !Cycle neighbouring vertices
            jvtx=vtxconn(ivtx,jtmp)
            if (elimvtx(jvtx)==1) cycle
            if (iselvtx(jvtx)==1) then
                iselvtx(ivtx)=1
                exit
            end if
        end do
    end do
    nselvtx=count(iselvtx==1)
    if (nselvtx==nselvtx_old) exit
end do

write(*,"(' Number of surface vertices in selected surface region:',i10)") nselvtx

areasel=0
avgvalsel=0
do icyc=1,nsurtri !If a triangle is composed of three selected vertices, the triangle should be taken into account
	if (elimtri(icyc)==1) cycle
	idx1=surtriang(icyc)%idx(1)
	idx2=surtriang(icyc)%idx(2)
	idx3=surtriang(icyc)%idx(3)
    if (iselvtx(idx1)==1.and.iselvtx(idx2)==1.and.iselvtx(idx3)==1) then
	    areasel=areasel+surtriang(icyc)%area
        avgvalsel=avgvalsel+surtriang(icyc)%value*surtriang(icyc)%area
    end if
end do
avgvalsel=avgvalsel/areasel
write(*,"(' Area of selected surface region:',f10.3,' Angstrom^2')") areasel*b2a**2
write(*,"(' Average value of selected surface region:',f12.5,' a.u.')") avgvalsel
write(*,"(' Product of above two values:',f16.5,' a.u.*Angstrom^2')") avgvalsel* areasel*b2a**2

open(10,file="selsurf.pqr",status="replace")
write(10,"('REMARK   Generated by Multiwfn, totally',i10,' surface vertices')") nselvtx
write(10,"(a)") "REMARK   The third last column is function values in a.u."
do i=1,nsurvtx
	if (elimvtx(i)==0.and.iselvtx(i)==1) then !Output as carbon atoms
        write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f14.8,f9.4,a2)") &
		"HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,survtx(i)%value,1.0," C"
	end if
end do
write(10,"('END')")
close(10)
write(*,"(a)") " Vertices of selected surface region have been exported to selsurf.pqr in current folder so that you can visualize them"

end subroutine





!!------ Perform basin-like partition on molecular surface
subroutine surfbasin
use defvar
use surfvertex
use util
implicit real*8 (a-h,o-z)
integer vtxatt(nsurvtx) !Attribution of a vertex to extreme index. =0 means not attributed
integer,allocatable :: extidx(:)
integer :: pseudoext(-10000:10000) !List of pseudo-extrema
integer :: mappseudo(-10000:10000) !Mapping pseudo-extreme indices to actual extreme indices
real*8,allocatable :: basinarea(:),basinavgval(:)

!Construct vertex index list of all surface extrema. Positive/Negative position index records maxima/minima vertex index
allocate(extidx(-nsurlocmin:nsurlocmax))
extidx(0)=0
do imax=1,nsurlocmax
    extidx(imax)=surlocmaxidx(imax)
end do
do imin=1,nsurlocmin
    extidx(-imin)=surlocminidx(imin)
end do

!do i=-nsurlocmin,nsurlocmax
!    write(*,*) i,extidx(i)
!end do

!The pseudo-minimum means the value at this point is smaller than directly connected vertices, but it is not
!smaller than the secondary neighbouring vertices, therefore it is not regarded as actual minimum
!However, vertex may finally move to these pseudo-minima, therefore we construct a list to map the pseudo-minima
!to actual minima, the rule is their distance is smaller than a threshold
npseudomin=0
npseudomax=0
pmin: do ivtx=1,nsurvtx
	if (elimvtx(ivtx)==1) cycle
    if (any(extidx==ivtx)) cycle
    vali=survtx(ivtx)%value
    do jtmp=1,vtxconnpos(ivtx) !Compare value with secondary neighbours
		jvtx=vtxconn(ivtx,jtmp)
		if (elimvtx(jvtx)==1) cycle
		if (vali>survtx(jvtx)%value) cycle pmin
	end do
    npseudomin=npseudomin+1
    pseudoext(-npseudomin)=ivtx
    !Calculate distance between this pseudo-minimum with real minimum, it will be mapped to the closest one
    dist2min=1D10
    do imin=1,nsurlocmin
        iminvtx=surlocminidx(imin)
        xext=survtx(iminvtx)%x;yext=survtx(iminvtx)%y;zext=survtx(iminvtx)%z
        xvtx=survtx(ivtx)%x;yvtx=survtx(ivtx)%y;zvtx=survtx(ivtx)%z
        dist2=(xext-xvtx)**2+(yext-yvtx)**2+(zext-zvtx)**2
        if (dist2<dist2min) then
            iclose=imin
            dist2min=dist2
        end if
    end do
    mappseudo(-npseudomin)=-iclose
end do pmin
!Similarly, process pseudo-maxima
pmax: do ivtx=1,nsurvtx
	if (elimvtx(ivtx)==1) cycle
    if (any(extidx==ivtx)) cycle
    vali=survtx(ivtx)%value
    do jtmp=1,vtxconnpos(ivtx) !Compare value with secondary neighbours
		jvtx=vtxconn(ivtx,jtmp)
		if (elimvtx(jvtx)==1) cycle
		if (vali<survtx(jvtx)%value) cycle pmax
	end do
    npseudomax=npseudomax+1
    pseudoext(npseudomax)=ivtx
    !Calculate distance between this pseudo-maximum with real maximum, it will be mapped to the closest one
    dist2min=1D10
    do imax=1,nsurlocmax
        imaxvtx=surlocmaxidx(imax)
        xext=survtx(imaxvtx)%x;yext=survtx(imaxvtx)%y;zext=survtx(imaxvtx)%z
        xvtx=survtx(ivtx)%x;yvtx=survtx(ivtx)%y;zvtx=survtx(ivtx)%z
        dist2=(xext-xvtx)**2+(yext-yvtx)**2+(zext-zvtx)**2
        if (dist2<dist2min) then
            iclose=imax
            dist2min=dist2
        end if
    end do
    mappseudo(npseudomax)=iclose
end do pmax
write(*,"(' Number of pseudo-maxima:',i6,'    Number of pseudo-minima:',i6)") npseudomin,npseudomax
!do i=-npseudomin,npseudomax
!    write(*,*) i,pseudoext(i),mappseudo(i)
!end do

!Determining attribution of each vertex
vtxatt=0
maxcyc=500
notelim=count(elimvtx==0)
iprog=0
nfailed=0
do ivtxcyc=1,nsurvtx
	if (elimvtx(ivtxcyc)==1) cycle
    ivtx=ivtxcyc !ivtx is a temporary variable recording current index of ivtxcyc
    icyc=0
    
movei:   do while(.true.) !Update position of ivtx until find an extrema
        idxmax=ivtx
        absmax=abs(survtx(ivtx)%value)
        !Cycle neighbouring vertices and update index
        do jtmp=1,vtxconnpos(ivtx)
		    jvtx=vtxconn(ivtx,jtmp)
            if (elimvtx(jvtx)==1) cycle 
            absjval=abs(survtx(jvtx)%value)
            if (absjval>absmax) then !Finding vertex carrying largest abs value around current vertex
                absmax=absjval
                idxmax=jvtx
            end if
        end do
        ivtx=idxmax
        icyc=icyc+1
        !Check if attribution of current vertex has already been determined, if yes, directly use its attribution
        if (vtxatt(ivtx)/=0) then
            vtxatt(ivtxcyc)=vtxatt(ivtx)
            exit movei
        end if
        !Check if current index is an actual extrema
        do icheck=-nsurlocmin,nsurlocmax
            if (icheck==0) cycle
            if (idxmax==extidx(icheck)) then
                vtxatt(ivtxcyc)=icheck !Attribution of ivtxcyc now is determined
                exit movei
            end if
        end do
        !Check if current index is a pseudo-extrema
        do icheck=-npseudomin,npseudomax
            if (icheck==0) cycle
            if (idxmax==pseudoext(icheck)) then
                vtxatt(ivtxcyc)=mappseudo(icheck) !Attribution of ivtxcyc now is determined
                exit movei
            end if
        end do
        if (icyc==maxcyc) then
            nfailed=nfailed+1
            exit
        end if
    end do movei
    
    iprog=iprog+1
    call showprog(iprog,notelim)
end do

!open(10,file="failed.txt",status="replace")
!do ivtx=1,nsurvtx
!	if (elimvtx(ivtx)==1) cycle
!    if (vtxatt(ivtx)==0) write(10,*) ivtx
!end do
!close(10)
if (nfailed/=0) write(*,"(' Determination of',i6,' vertices is failed')") nfailed
write(*,*)

allocate(basinarea(-nsurlocmin:nsurlocmax),basinavgval(-nsurlocmin:nsurlocmax))
basinarea=0
basinavgval=0
!!If the three vertices of a triangle share the same extrema index, then value of this triangle will be &
!attributed to corresponding basin. Clearly boundary triangle is not attributed to any basin
do itri=1,nsurtri
	if (elimtri(itri)==1) cycle
	idx1=surtriang(itri)%idx(1)
	idx2=surtriang(itri)%idx(2)
	idx3=surtriang(itri)%idx(3)
    iext=vtxatt(idx1)
    if (vtxatt(idx2)==iext.and.vtxatt(idx3)==iext) then
	    basinarea(iext)=basinarea(iext)+surtriang(itri)%area
        basinavgval(iext)=basinavgval(iext)+surtriang(itri)%value*surtriang(itri)%area
    end if
end do
basinavgval=basinavgval/basinarea
basinarea=basinarea*b2a*b2a

do imin=1,nsurlocmin
    Nvert=count(vtxatt==-imin)
    if (Nvert==0) cycle
    write(*,"(' Minimum',i4,'  N_vert:',i6,',',f8.3,' Angstrom^2  Avg. value:',f12.6,' a.u.')") &
    imin,Nvert,basinarea(-imin),basinavgval(-imin)
end do
write(*,*)
do imax=1,nsurlocmax
    Nvert=count(vtxatt==imax)
    if (Nvert==0) cycle
    write(*,"(' Maximum',i4,'  N_vert:',i6,',',f8.3,' Angstrom^2  Avg. value:',f12.6,' a.u.')") &
    imax,Nvert,basinarea(imax),basinavgval(imax)
end do

open(10,file="surfbasin.pdb",status="replace")
write(10,"('REMARK   Generated by Multiwfn')")
if (ifPBC>0) then
	call getcellabc(asize,bsize,csize,alpha,beta,gamma)
	write(10,"('CRYST1',3f9.3,3f7.2)") asize,bsize,csize,alpha,beta,gamma
end if
do i=1,nsurvtx
	if (elimvtx(i)==0) then !Output as carbon atoms
        betaval=vtxatt(i)
		write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
		"HETATM",i,' '//"C "//' ',"MOL",'A',1,survtx(i)%x*b2a,survtx(i)%y*b2a,survtx(i)%z*b2a,1.0,betaval,"C "
	end if
end do
write(10,"('END')")
close(10)
write(*,"(/,a)") " Surface basins have been exported to surfbasin.pdb in current folder, &
&beta values correspond to indices of extrema, negative/positive index corresponds to minima/maxima"

end subroutine