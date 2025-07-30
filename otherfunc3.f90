!-------- Main interface of various other functions part 2
subroutine otherfunc3_main
implicit real*8 (a-h,o-z)

do while(.true.)
	write(*,*)
	write(*,*) "              ============ Other functions (Part 3) ============ "
	write(*,*) "0 Return"
	write(*,*) "1 Viewing free regions and calculating free volume in a cell"
	write(*,*) "2 Fitting atomic radial density as linear combination of multiple STOs or GTFs"
    !write(*,*) "3 Visualize (hyper)polarizability via unit sphere and vector representations"
    write(*,*) "4 Simulating scanning tunneling microscope (STM) image"
    write(*,*) "5 Calculate electric dipole/multipole moments and electronic spatial extent"
    write(*,*) "6 Calculate energies of present orbitals by inputting Fock matrix"
    write(*,*) "7 Geometry operation on the present system"
    write(*,*) "8 Plot surface distance projection map"
    write(*,*) "9 Determine Fermi level"
	read(*,*) isel
	if (isel==0) then
		return
	else if (isel==1) then
        call freeregion
	else if (isel==2) then
        call fitatmdens
	else if (isel==3) then !Hidden
        call vis_hypol
	else if (isel==4) then
        call STM
    else if (isel==5) then
        call calc_multipole
    else if (isel==6) then
        call calc_orb_energy
    else if (isel==7) then
        call geom_operation
    else if (isel==8) then
        call molsurf_distmap
    else if (isel==9) then
        call calc_Fermi_level
	end if
end do
end subroutine



!!--------- Viewing free regions and calculating free volume in a box
!  Error function with scale factor of 1.0 is used. Because Becke function even with n_iter=1 is to sharp, isosurface will even occur &
!in the center region between two stacking aromatic rings. While Gaussian function decays to slow to be employed as a general choice, &
!it is quite artificial when lots of atoms are crowded in small area, because their effect are accumulated too much in the region far away from nuclei
subroutine freeregion
use defvar
use GUI
use util
implicit real*8 (a-h,o-z)
character c80tmp*80
integer :: ismooth=1,iclosebound=1,ismoothmethod=3,nBeckeiter=1
real*8 :: sclrad=1,xyzA(3),xyzB(3),Gauss_radmulti=2D0,erfscl=1D0
real*8,allocatable :: cubmat_smooth(:,:,:)

if (ifPBC>0) then
    iPBC=1
else
    iPBC=0
end if

do while(.true.)
    do while(.true.)
        write(*,*)
        write(*,*) " -------- Viewing free regions and calculating free volume in a cell --------"
        write(*,*) "0 Return"
        write(*,*) "1 Set grid and start calculation"
        if (iPBC==1) write(*,*) "2 Toggle considering periodic boundary condition, current: Yes"
        if (iPBC==0) write(*,*) "2 Toggle considering periodic boundary condition, current: No"
        if (ismooth==1) then
            write(*,*) "3 Toggle calculating smoothed grid data of free regions, current: Yes"
            if (ismoothmethod==1) write(*,"(a,f4.1)") " 4 Set method of smoothing, current: Gaussian function, FWHM is vdW radius multiplied by",Gauss_radmulti
            if (ismoothmethod==2) write(*,"(a,i2)") " 4 Set method of smoothing, current: Becke function, number of iterations is",nBeckeiter
            if (ismoothmethod==3) write(*,"(a,f5.2)") " 4 Set method of smoothing, current: Error function, scale factor:",erfscl
        else
            write(*,*) "3 Toggle calculating smoothed grid data of free regions, current: No"
        end if
        if (iclosebound==0) write(*,"(' 5 Toggle making isosurface closed at boundary, current: No')")
        if (iclosebound==1) write(*,"(' 5 Toggle making isosurface closed at boundary, current: Yes')")
        write(*,"(' 6 Set scale factor of vdW radii for calculating free volume and raw free region, current:',f6.3)") sclrad
    
        read(*,*) isel
        if (isel==0) then
            return
        else if (isel==1) then
            exit
        else if (isel==2) then
            if (iPBC==1) then
                iPBC=0
            else
                iPBC=1
            end if
        else if (isel==3) then
            if (ismooth==1) then
                ismooth=0
            else
                ismooth=1
            end if
        else if (isel==4) then
            write(*,*) "Choose the switching function used to calculate smoothed grid data"
            write(*,*) "1 Gaussian function"
            write(*,*) "2 Becke function (transformed)"
            write(*,*) "3 Error function (transformed)"
            write(*,"(a)") " Note: 2 and 3 have very similar decaying character and computational cost, see Section 3.300.1 of Multiwfn manual"
            read(*,*) ismoothmethod
            if (ismoothmethod==1) then
                write(*,"(a)") " Set full width at half maximum (FWHM) of Gaussian function. If inputting for example 1.8, &
                &then the FWHM will be 1.8 times of van der Waals radius of corresponding atom"
                write(*,*) "Hint: Usually the value should be <=2.0"
                read(*,*) Gauss_radmulti
            else if (ismoothmethod==2) then
                write(*,*) "Input number of iterations of Becke function, e.g. 2"
                write(*,"(a)") " Hint: Should be >=0. The smaller the number, the smoother the surface, &
                &but the farther region an atom will affect and thus more artificial it is. Usually 1 is used"
                read(*,*) nBeckeiter
            else if (ismoothmethod==3) then
                write(*,*) "Input scale factor of the error function, e.g. 0.3"
                write(*,"(a)") " Hint: The larger the value, the smoother the surface, &
                &but the farther region an atom will affect and thus more artificial it is. Usually 1.0 is a good compromise"
                read(*,*) erfscl
            end if
        else if (isel==5) then
            if (iclosebound==1) then
                iclosebound=0
            else
                iclosebound=1
            end if
        else if (isel==6) then
            write(*,*) "Input scale factor for vdW radii, e.g. 1.8"
            write(*,*) "Note: This setting does not affect smoothed grid data"
            read(*,*) sclrad
        end if
    end do

    !Set grid
    write(*,*) "Input X, Y, Z of origin of grid data in Angstrom, e.g. 0.3,4.5,0"
    write(*,*) "If press ENTER button directly, (0,0,0) will be used"
    read(*,"(a)") c80tmp
    if (c80tmp==" ") then
        orgx=0;orgy=0;orgz=0
    else
        read(c80tmp,*) orgx,orgy,orgz
        orgx=orgx/b2a
        orgy=orgy/b2a
        orgz=orgz/b2a
    end if
    write(*,*) "Input lengths of the three sides of the grid data range in Angstrom"
    write(*,*) "For example, 31.06,31.15,31.092"
    if (ifPBC==3) write(*,*) "If press ENTER button directly, length of three sides of the cell will be used"
    read(*,"(a)") c80tmp
    if (c80tmp==" ") then
        v1len=dsqrt(sum(cellv1**2))
        v2len=dsqrt(sum(cellv2**2))
        v3len=dsqrt(sum(cellv3**2))
    else
        read(c80tmp,*) v1len,v2len,v3len
        v1len=v1len/b2a
        v2len=v2len/b2a
        v3len=v3len/b2a
    end if
    write(*,*) "Input grid spacing in Angstrom, e.g. 0.2"
    write(*,*) "If press ENTER button directly, 0.25 Angstrom will be used"
    write(*,"(a)") " Note: The smaller the value, the more accurate the result and the smoother of the isosurface, but the higher the computational cost"
    read(*,"(a)") c80tmp
    if (c80tmp==" ") then
        grdspc=0.25D0
    else
        read(c80tmp,*) grdspc
    end if
    grdspc=grdspc/b2a
    nx=nint(v1len/grdspc)
    ny=nint(v2len/grdspc)
    nz=nint(v3len/grdspc)
    !Marginally re-adjust grid spacing so that the number of grids could be integer
    grdspcv1=v1len/(nx-1)
    grdspcv2=v2len/(ny-1)
    grdspcv3=v3len/(nz-1)
    if (ifPBC==3) then
        gridv1(:)=cellv1(:)/v1len*grdspcv1
        gridv2(:)=cellv2(:)/v2len*grdspcv2
        gridv3(:)=cellv3(:)/v3len*grdspcv3
    else
        gridv1(:)=(/ grdspcv1,0D0,0D0 /)
        gridv2(:)=(/ 0D0,grdspcv2,0D0 /)
        gridv3(:)=(/ 0D0,0D0,grdspcv3 /)
    end if
    call getgridend
	write(*,"(' Coordinate of origin in X,Y,Z is   ',3f12.6,' Angs')") orgx*b2a,orgy*b2a,orgz*b2a
	write(*,"(' Coordinate of end point in X,Y,Z is',3f12.6,' Angs')") endx*b2a,endy*b2a,endz*b2a
	write(*,"(' Grid vector 1 in X,Y,Z is',3f10.6,' Angs, norm:',f10.6)") gridv1*b2a,dsqrt(sum(gridv1**2))*b2a
	write(*,"(' Grid vector 2 in X,Y,Z is',3f10.6,' Angs, norm:',f10.6)") gridv2*b2a,dsqrt(sum(gridv2**2))*b2a
	write(*,"(' Grid vector 3 in X,Y,Z is',3f10.6,' Angs, norm:',f10.6)") gridv3*b2a,dsqrt(sum(gridv3**2))*b2a
	write(*,"(' Number of points in three directions is',3i5,'  Total:',i12)") nx,ny,nz,nx*ny*nz

    call walltime(iwalltime1)
    allocate(cubmat(nx,ny,nz))
    if (ismooth==1) allocate(cubmat_smooth(nx,ny,nz))
    cubmat=1 !Free region has value of 1, occupied region has value of 0
    cubmat_smooth=1
    write(*,*)
    write(*,*) "Calculating, please wait..."
    ifinish=0
    !$OMP PARALLEL DO SHARED(cubmat,cubmat_smooth,ifinish) PRIVATE(i,j,k,tmpx,tmpy,tmpz,iatm,atmvdwr,sclvdwr,xyzA,xyzB,dist,tmpval,parmc,iter,tmps) schedule(dynamic) NUM_THREADS(nthreads)
    do k=1,nz
	    do j=1,ny
		    do i=1,nx
                call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                do iatm=1,ncenter
                    atmvdwr=vdwr(a(iatm)%index)
                    sclvdwr=sclrad*atmvdwr
                    if (iPBC==1) then
                        xyzA(1)=tmpx
                        xyzA(2)=tmpy
                        xyzA(3)=tmpz
                        xyzB(1)=a(iatm)%x
                        xyzB(2)=a(iatm)%y
                        xyzB(3)=a(iatm)%z
                        call nearest_dist(xyzA,xyzB,dist)
                    else
                        dist=dsqrt((tmpx-a(iatm)%x)**2+(tmpy-a(iatm)%y)**2+(tmpz-a(iatm)%z)**2)
                    end if
                    if (dist<sclvdwr) then
                        cubmat(i,j,k)=0
                        if (ismooth==0) exit
                    end if
                    !Generate smoothed grid data
                    if (ismooth==1) then
                        if (ismoothmethod==1) then !Use (unnormalized) Gaussian function to smooth, namely broadening atoms as Gaussian distributions
                            !When atoms are very crowd, e.g. C60 .cif, the smoothed function will be negative everywhere, because Gaussian is too extent
                            !If FWHM = 2*atmvdwr, the value of Gaussian at atmvdwr is 0.5
                            if (dist>3*atmvdwr*(Gauss_radmulti/2)) cycle !Ignore very far atom from this point. This was found to be very safe
                            tmpval=switch_Gauss(dist,Gauss_radmulti*atmvdwr)
                        else if (ismoothmethod==2) then !Use Becke switching function for smoothing, vdW radius position corresponds to 0.5
                            if (dist>2*atmvdwr) cycle !Ignore very far atom from this point. This was found to be very safe
                            tmpval=switch_Becke(dist,atmvdwr,nBeckeiter)
                        else if (ismoothmethod==3) then !Use switching function based on error function for smoothing
                            if (dist>2*atmvdwr) cycle
                            tmpval=switch_erf(dist,atmvdwr,erfscl)
                        end if
                        cubmat_smooth(i,j,k)=cubmat_smooth(i,j,k)-tmpval
                    end if
                end do
		    end do
	    end do
        !$OMP CRITICAL
	    ifinish=ifinish+1
        call showprog(ifinish,nz)
        !$OMP END CRITICAL
    end do
    !$OMP END PARALLEL DO

    where(cubmat_smooth<0) cubmat_smooth=0

    !Calculate free volume. Strickly speaking, center of grid should be used to determine occupancy status, however for simplicity,
    !we do below way, corresponding to use one vertex of a grid to determine occupancy
    call calc_dvol(dvol)
    freevol=sum(cubmat(1:nx-1,1:ny-1,1:nz-1))*dvol*b2a**3
    voltot=(nx-1)*(ny-1)*(nz-1)*dvol*b2a**3
    write(*,"(' Volume of entire box:',f12.3,' Angstrom^3')") voltot
    write(*,"(' Free volume:',f12.3,' Angstrom^3, corresponding to',f8.2,' % of whole space')") freevol,freevol/voltot*100

    call walltime(iwalltime2)
    write(*,"(/,' Calculation took up time',i10,' s')") iwalltime2-iwalltime1

    !Set value of boundary grids to make isosurface map close
    if (iclosebound==1) then
        !XY layers
        do i=1,nx
            do j=1,ny
                cubmat(i,j,1)=0
                cubmat(i,j,nz)=0
                if (ismooth==1) then
                    cubmat_smooth(i,j,1)=0
                    cubmat_smooth(i,j,nz)=0
                end if
            end do
        end do
        !XZ layers
        do i=1,nx
            do k=1,nz
                cubmat(i,1,k)=0
                cubmat(i,ny,k)=0
                if (ismooth==1) then
                    cubmat_smooth(i,1,k)=0
                    cubmat_smooth(i,ny,k)=0
                end if
            end do
        end do
        !YZ layers
        do j=1,ny
            do k=1,nz
                cubmat(1,j,k)=0
                cubmat(nx,j,k)=0
                if (ismooth==1) then
                    cubmat_smooth(1,j,k)=0
                    cubmat_smooth(nx,j,k)=0
                end if
            end do
        end do
    end if

    idrawmol=1
    ishowatmlab=0
    ishowaxis=0
    sur_value=0.5D0
    do while(.true.)
        write(*,*)
        write(*,*) "                     ------- Post-processing menu -------"
        write(*,*) "0 Return"
        write(*,*) "1 Visualize isosurface of raw free region"
        write(*,*) "2 Export raw grid data as free_prim.cub in current folder"
        if (ismooth==1) then
            write(*,*) "3 Visualize isosurface of smoothed free region"
            write(*,*) "4 Export smoothed grid data as free_smooth.cub in current folder"
        end if
        read(*,*) isel
        if (isel==0) then
            deallocate(cubmat)
            if (ismooth==1) deallocate(cubmat_smooth)
            exit
        else if (isel==1) then
            call drawisosurgui(1)
        else if (isel==2) then
            open(10,file="free_prim.cub",status="replace")
            call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
            close(10)
            write(*,"(' Done! Grid data has been exported to free_prim.cub in current folder')")
        else if (isel==3) then
            allocate(cubmattmp(nx,ny,nz))
            cubmattmp=cubmat
            cubmat=cubmat_smooth
            call drawisosurgui(1)
            cubmat=cubmattmp
            deallocate(cubmattmp)
        else if (isel==4) then
            open(10,file="free_smooth.cub",status="replace")
            call outcube(cubmat_smooth,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
            close(10)
            write(*,"(' Done! Grid data has been exported to free_smooth.cub in current folder')")
        end if
    end do
end do
end subroutine



!!------- Module used by fitatmdens and related routines
module fitatmdens_mod
real*8,allocatable :: radr(:),radrho(:) !Radial distance and sphericallized density
integer :: ifittype=2 !=1: Minimize absolute error, =2: Minimize relative error, =3: Minimize RDF error
integer :: ifunctype=1 !=1: STO, =2: GTF
integer,parameter :: maxfitfunc=1000 !Maximum number of fitting functions
real*8 exp_fit(maxfitfunc) !If exponents are requested to be fixed, use this array to pass exponents to the calculation routine
integer :: ifixexp=0 !=0: Fit both coefficients and exponents, =1: Fix exponents unchanged
integer :: nfitfunc=0 !Number of STO/GTF used for fitting
end module

!!------------ Use Levenberg-Marquardt algorithm to fit sphericallized atomic radial density as multiple STOs or GTFs
subroutine fitatmdens
use defvar
use fitatmdens_mod
use functions
use util
use plot
implicit real*8 (a-h,o-z)
character(len=3) :: funclab(2)=(/ "STO","GTF" /)
character clegend*80,c200tmp*200,c80tmp*80
real*8 parm(maxfitfunc*2) !Up to maxfitfunc fitting functions. The first half is coefficient, the latter half is exponents
real*8,allocatable :: fiterr(:),fitrho(:) !Fitting error and fitted density at each fitting point
integer :: npoint_CB=0
real*8,allocatable :: radr_CB(:),radw_CB(:),rho_CB(:) !Position, weight and sphericalized density at second kind Gauss-Chebyshev points
real*8 :: tol=1D-5 !Fitting tolerance. Should not be too small, otherwise it is too difficult to converge until reach maximum of function calls
integer :: iscale=1 !=1: Scale coefficients so that integral equals to actual number of electrons, =0: Do not scale
integer :: isort=1,idelredun=1
external :: atmdens_fiterr
integer,parameter :: nsphpt=170 !Number of points used to calculate sphericalized radial density
real*8 potx(nsphpt),poty(nsphpt),potz(nsphpt),potw(nsphpt)
integer seqidx(maxfitfunc)
integer,allocatable :: tmpidxarr(:)
real*8,allocatable :: tmpvalarr(:)

if (ncenter/=1.or.a(1)%x/=0.or.a(1)%y/=0.or.a(1)%z/=0) then
    write(*,*) "Error: In order to use this function, there must be only one atom and at (0,0,0) point!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
end if

maxcall=10000
expcutoff=1 !Disable exponent truncation
radstep=0.001D0/b2a
npoint=4000
!radstep=0.002D0/b2a
!npoint=2000
call delvirorb(1)
call gen_GTFuniq(0) !Generate unique GTFs, for faster evaluation in orbderv

do while(.true.)
    write(*,*)
    write(*,*) "--------------- Fitting atomic radial density as STOs or GTFs ---------------"
    write(*,*) "0 Return"
    write(*,*) "1 Start fitting"
    write(*,"(a,a)") " 2 Switch type of fitting functions, current: ",funclab(ifunctype)
    write(*,*) "3 Check or set initial guess of coefficients and exponents"
    write(*,"(a,1PE10.2)") " 4 Set fitting tolerance, current:",tol
    write(*,"(a,i5)") " 5 Set number of evenly placed fitting points, current:",npoint
    write(*,"(a,f6.3,' Angstrom')") " 6 Set spacing between fitting points, current:",radstep*b2a
    if (iscale==0) write(*,*) "7 Toggle scaling coefficients to actual number of electrons, current: No"
    if (iscale==1) write(*,*) "7 Toggle scaling coefficients to actual number of electrons, current: Yes"
    if (ifittype==1) write(*,*) "8 Select fitting type, current: Minimizing absolute error"
    if (ifittype==2) write(*,*) "8 Select fitting type, current: Minimizing relative error"
    if (ifittype==3) write(*,*) "8 Select fitting type, current: Minimizing RDF error"
    if (ifixexp==0) write(*,*) "9 Toggle fixing exponents, current: No"
    if (ifixexp==1) write(*,*) "9 Toggle fixing exponents, current: Yes"
    if (isort==0) write(*,*) "10 Toggle sorting functions according to exponents, current: No"
    if (isort==1) write(*,*) "10 Toggle sorting functions according to exponents, current: Yes"
    if (idelredun==0) write(*,*) "11 Toggle removing redundant fitting functions, current: No"
    if (idelredun==1) write(*,*) "11 Toggle removing redundant fitting functions, current: Yes"
    write(*,"(a,i4)") " 12 Set number of second kind Gauss-Chebyshev fitting points, current:",npoint_CB
    write(*,"(a,i10)") " 13 Set maximum number of function calls, current:",maxcall
    read(*,*) isel
    
    if (isel==0) then
        call del_GTFuniq !Destory unique GTF informtaion
	    call delvirorb_back(1)
        return
    else if (isel==2) then
        if (ifunctype==1) then
            ifunctype=2
        else
            ifunctype=1
        end if
    else if (isel==3) then
        do while(.true.)
            if (nfitfunc>0) then
                write(*,*)
                write(*,*) "Current initial guess:"
                write(*,*) "           Coefficient       Exponent"
                do ifunc=1,nfitfunc
                    write(*,"(1x,a,i3,':',2(1PE16.6))") funclab(ifunctype),ifunc,parm(ifunc),parm(nfitfunc+ifunc)
                end do
            else
                write(*,*) "Note: Fitting functions have not been set"
            end if
            write(*,*)
            write(*,*) "0 Return"
            write(*,*) "1 Load initial guess from text file"
            write(*,*) "  Prebuilt settings:"
            write(*,*) "2 Crude fitting by a few STOs with variable exponents"
            write(*,*) "3 Ideal fitting by 30 GTFs with fixed exponents"
            write(*,*) "4 Fine fitting by 15 GTFs with variable exponents"
            write(*,*) "5 Fine fitting by 10 GTFs with variable exponents"
            !write(*,*) "6 Set initial guess as ""fine fitting by 20 GTFs with variable exponents"""
            !GTF with variable exponents higher than 15 is not a good idea, it is extremely expensive &
            !and often the lmdif1 routine returns unoptimized result or convergence tolerance cannot be reached, &
            !and the result is never detectably better than 15 GTFs
            write(*,*) "7 Relatively fine fitting by no more than 10 GTFs with variable exponents"
            write(*,*) "Hint: Accuracy & number of functions: 3>4>5>7>2"
            if (nfitfunc>0) write(*,*) "10 Combine two fitting functions together"
            read(*,*) isel2
            
            if (isel2==0) then
                exit
            else if (isel2==1) then
                write(*,"(/,a)") " Input path of the text file containing initial guess of coefficients and exponents, e.g. C:\Popipa.txt"
                do while(.true.)
                    read(*,"(a)") c200tmp
	                inquire(file=c200tmp,exist=alive)
	                if (alive) exit
	                write(*,*) "Cannot find the file, input again!"
                end do
                open(10,file=c200tmp,status="old")
                nfitfunc=totlinenum(10,1)
                if (nfitfunc>maxfitfunc) then
                    write(*,"(' Error: Number of fitting functions should not exceed',i6)") maxfitfunc
                    write(*,*) "Press ENTER button to cancel loading"
                    read(*,*)
                else
                    do ifunc=1,nfitfunc
                        read(10,*) parm(ifunc),parm(nfitfunc+ifunc)
                    end do
                end if
                close(10)
            else if (isel2==2) then !Initial guess of few STOs
                ifixexp=0
                ifunctype=1
                idelredun=1
                if (a(1)%index<=2) then
                    nfitfunc=1
                    parm(1)=1D0
                    parm(2)=2D0
                else if (a(1)%index<=10) then !For second row, using more STOs does not improve result
                    nfitfunc=2
                    parm(1:nfitfunc)=(/ 100,1 /)
                    parm(nfitfunc+1:2*nfitfunc)=(/ 10,2 /)
                else if (a(1)%index<=36) then
                    !nfitfunc=3 !Qualtiy is not as good as nfitunc=4
                    !parm(1:nfitfunc)=(/ 1000,100,1 /)
                    !parm(nfitfunc+1:2*nfitfunc)=(/ 20,5,1 /)
                    nfitfunc=4 !For Mn, this setting may result in negative value at some points
                    parm(1:nfitfunc)=(/ 1000,300,20,1 /)
                    parm(nfitfunc+1:2*nfitfunc)=(/ 27,9,3,1 /)
                else
                    nfitfunc=6
                    do ifunc=nfitfunc,1,-1
                        parm(ifunc)=1*15**(ifunc-1)
                    end do
                    do ifunc=1,nfitfunc
                        parm(nfitfunc+ifunc)=0.5*3**(ifunc)
                    end do
                end if
            else if (isel2==3) then
                ifixexp=1
                ifunctype=2
                nfitfunc=30
                idelredun=1
                do ifunc=1,nfitfunc
                    !parm(nfitfunc+ifunc)=0.1D0*2D0**(ifunc-1) !Using this is safer, but tail cannot be represented as well as below setting
                    parm(nfitfunc+ifunc)=0.05D0*2D0**(ifunc-1)
                end do
                parm(1:nfitfunc)=1
            else if (isel2==4.or.isel2==5.or.isel2==6) then
                ifixexp=0
                ifunctype=2
                idelredun=1
                if (isel2==4) then
                    nfitfunc=15
                    do ifunc=1,nfitfunc
                        parm(nfitfunc+ifunc)=0.1D0*2D0**(ifunc-1)
                    end do
                else if (isel2==5) then
                    nfitfunc=10
                    do ifunc=1,nfitfunc
                        parm(nfitfunc+ifunc)=0.1D0*2.5D0**(ifunc-1)
                    end do
                !else if (isel2==6) then
                !    nfitfunc=20
                !    do ifunc=1,nfitfunc
                !        parm(nfitfunc+ifunc)=0.05D0*1.6D0**(ifunc-1)
                !    end do
                end if
                parm(1:nfitfunc)=1
                tol=1E-4 !Use more loose tolerance than default make convergence easier while the quality is not detectably lowered
            else if (isel2==7) then !Initial guess of few GTFs
                ifixexp=0
                ifunctype=2
                idelredun=1
                if (a(1)%index<=18) then
                    nfitfunc=6
                    do ifunc=1,nfitfunc
                        parm(nfitfunc+ifunc)=0.15D0*2.5D0**(ifunc-1)
                    end do
                else
                    nfitfunc=10
                    do ifunc=1,nfitfunc
                        parm(nfitfunc+ifunc)=0.1D0*2.5D0**(ifunc-1)
                    end do
                end if
                parm(1:nfitfunc)=1
            else if (isel2==10) then
                write(*,"(a)") " Input indices of the functions to combine (e.g. 1,3), whose average exponent will be taken as the new exponent, &
                &and their coefficients will be summed up"
                read(*,*) ifunc1,ifunc2
                parm(ifunc1)=parm(ifunc1)+parm(ifunc2)
                parm(nfitfunc+ifunc1)=(parm(nfitfunc+ifunc1)+parm(nfitfunc+ifunc2))/2
                parm(ifunc2:ifunc2+2*nfitfunc)=parm(ifunc2+1:ifunc2+2*nfitfunc+1)
                nfitfunc=nfitfunc-1
                parm(nfitfunc+ifunc2:nfitfunc+ifunc2+nfitfunc)=parm(nfitfunc+ifunc2+1:nfitfunc+ifunc2+1+nfitfunc)
                write(*,*) "Done!"
            end if
        end do
        maxcall=10000*2*nfitfunc
    else if (isel==4) then
        write(*,*) "Input fitting tolerance, e.g. 1E-7"
        write(*,*) "The smaller the value, the better the fitting accuracy while higher the cost"
        read(*,*) tol
    else if (isel==5) then
        write(*,*) "Input number of fitting points, e.g. 300"
        read(*,*) npoint
    else if (isel==6) then
        write(*,*) "Input spacing between fitting points (in Angstrom), e.g. 0.02"
        read(*,*) radstep
        radstep=radstep/b2a
    else if (isel==7) then
        if (iscale==1) then
            iscale=0
        else
            iscale=1
        end if
    else if (isel==8) then
        write(*,*) "1 Minimizing absolute error"
        write(*,*) "2 Minimizing relative error"
        write(*,*) "3 Minimizing radial distribution function (RDF) error"
        read(*,*) ifittype
    else if (isel==9) then
        if (ifixexp==1) then
            ifixexp=0
        else
            ifixexp=1
        end if
    else if (isel==10) then
        if (isort==1) then
            isort=0
        else
            isort=1
        end if
    else if (isel==11) then
        if (idelredun==1) then
            idelredun=0
        else
            idelredun=1
        end if
    else if (isel==12) then
        write(*,*) "Set the number of second kind Gauss-Chebyshev points used in fitting, e.g. 80"
        write(*,*) "If input 0, then this kind of points will not be included in fitting"
        read(*,*) npoint_CB
    else if (isel==13) then
        write(*,*) "Input maximum number of function calls, e.g. 300000"
        write(*,"(a)") " Note: If the error minimization does not converge when reaches this condition, the minimization will stop and unconverged result will be reported"
        read(*,*) maxcall
    else if (isel==1) then !Do fitting
    
        if (nfitfunc==0) then
            write(*,*) "Error: You should use option 3 to set initial guess of fitting functions first!"
            write(*,*) "Press ENTER button to return"
            read(*,*)
            cycle
        end if
    
        !Output initial guess
        write(*,*) "Initial guess of fitting functions:"
        write(*,*) "           Coefficient       Exponent"
        do ifunc=1,nfitfunc
            write(*,"(1x,a,i3,':',2(1PE16.6))") funclab(ifunctype),ifunc,parm(ifunc),parm(nfitfunc+ifunc)
        end do
        write(*,*)
        
        !Calculate actual sphericalized radial density at evenly distributed fitting points
        allocate(radr(npoint),radrho(npoint),fiterr(npoint),fitrho(npoint))
        write(*,*) "Calculating sphericalized radial density..."
        call Lebedevgen(nsphpt,potx,poty,potz,potw)
        !$OMP PARALLEL DO SHARED(radr,radrho) PRIVATE(ipt,tmpdens,isph,xtmp,ytmp,ztmp) schedule(dynamic) NUM_THREADS(nthreads)
        do ipt=1,npoint
            radr(ipt)=radstep*ipt
	        tmpdens=0
	        do isph=1,nsphpt
		        xtmp=potx(isph)*radr(ipt)
		        ytmp=poty(isph)*radr(ipt)
		        ztmp=potz(isph)*radr(ipt)
		        tmpdens=tmpdens+fdens(xtmp,ytmp,ztmp)*potw(isph)
	        end do
            radrho(ipt)=tmpdens
        end do
        !$OMP END PARALLEL DO
        
        !Take second kind Gauss-Chebyshev fitting points into account
        if (npoint_CB>0) then
            allocate(radr_CB(npoint_CB),radw_CB(npoint_CB),rho_CB(npoint_CB))
            !Calculate density at second kind Gauss-Chebyshev points
            parmbk=1D0
            do i=1,npoint_CB !Combine spherical point&weights with second kind Gauss-Chebyshev method for radial part
	            radx=cos(i*pi/(npoint_CB+1))
	            radr_CB(i)=(1+radx)/(1-radx)*parmbk !Becke transform
	            radw_CB(i)=2*pi/(npoint_CB+1)*parmbk**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
            end do
            !$OMP PARALLEL DO SHARED(rho_CB) PRIVATE(ipt,tmpdens,isph,xtmp,ytmp,ztmp) schedule(dynamic) NUM_THREADS(nthreads)
            do ipt=1,npoint_CB
	            tmpdens=0
	            do isph=1,nsphpt
		            xtmp=potx(isph)*radr_CB(ipt)
		            ytmp=poty(isph)*radr_CB(ipt)
		            ztmp=potz(isph)*radr_CB(ipt)
		            tmpdens=tmpdens+fdens(xtmp,ytmp,ztmp)*potw(isph)
	            end do
                rho_CB(ipt)=tmpdens
            end do
            !$OMP END PARALLEL DO
            !Merge second kind Gauss-Chebyshev points into fitting points
            nadd=(count(radr_CB<(10/b2a))) !Ignore points farther than 10 Angstrom, their densities are negligible
            allocate(tmpvalarr(npoint))
            tmpvalarr=radr
            deallocate(radr);allocate(radr(npoint+nadd))
            radr(1:npoint)=tmpvalarr
            tmpvalarr=radrho
            deallocate(radrho);allocate(radrho(npoint+nadd))
            radrho(1:npoint)=tmpvalarr
            do ipt=1,npoint_CB !Combine spherical point&weights with second kind Gauss-Chebyshev method for radial part
                if (radr_CB(ipt)>10/b2a) cycle !Ignore > 10 Angstrom points
                npoint=npoint+1
	            radr(npoint)=radr_CB(ipt)
                radrho(npoint)=rho_CB(ipt)
            end do
            deallocate(fitrho,fiterr)
            allocate(fitrho(npoint),fiterr(npoint))
            !Sort according to r
            allocate(tmpidxarr(npoint))
            forall(i=1:npoint) tmpidxarr(i)=i
            call sortr8(radr,list=tmpidxarr)
            fitrho(:)=fitrho(tmpidxarr(:))
            radrho(:)=radrho(tmpidxarr(:))
            deallocate(tmpvalarr,tmpidxarr)
        end if
        
        !!!!!! Start fitting !!!!!!
        nparm=nfitfunc*2
        write(*,"(/,' Maximum number of function calls:',i8)") maxcall
        write(*,"(' Convergence tolerance:',f16.8)") tol
        if (ifixexp==1) write(*,*) "Exponents are kept fixed as requested"
        if (ifittype==1) write(*,*) "Fitting type: Minimizing absolute error"
        if (ifittype==2) write(*,*) "Fitting type: Minimizing of relative error"
        if (ifittype==3) write(*,*) "Fitting type: Minimizing radial distribution function error"
        write(*,*) "Fitting via Levenberg-Marquardt algorithm..."
        write(*,*)
        if (idelredun==0) then !Simple fitting
            if (ifixexp==0) then !Fit both coefficients and exponents
                call lmdif1(atmdens_fiterr,npoint,nparm,parm(1:nparm),fiterr,tol,maxcall,info)
            else !Do not fit exponents but keep initial values
                exp_fit(1:nfitfunc)=parm(nfitfunc+1:nparm)
                call lmdif1(atmdens_fiterr,npoint,nfitfunc,parm(1:nfitfunc),fiterr,tol,maxcall,info)
            end if
        else if (idelredun==1) then !Fitting with automatical removal of redundant fitting functions
            nremove=0
            do while(.true.) !Repeat until no function is to be removed
                if (ifixexp==0) then !Fit both coefficients and exponents
                    call lmdif1(atmdens_fiterr,npoint,nparm,parm(1:nparm),fiterr,tol,maxcall,info)
                else !Do not fit exponents but keep initial values
                    exp_fit(1:nfitfunc)=parm(nfitfunc+1:nparm)
                    call lmdif1(atmdens_fiterr,npoint,nfitfunc,parm(1:nfitfunc),fiterr,tol,maxcall,info)
                end if
                iremove=0
                !Remove redundant fitting functions according to exponent and coefficient, adapted according to the ChkRed routine in denfit.f90 in molden2aim
                do ifunc=1,nfitfunc
                    coefftmp=parm(ifunc)
                    exptmp=parm(nfitfunc+ifunc)
                    if ((exptmp>1D5.and.abs(coefftmp)<5).or.(exptmp<3.and.abs(coefftmp)<1D-4)) then
                        iremove=ifunc
                        exit
                    end if    
                end do
                !Remove redundant fitting function having maximum negative contribution if this point has negative fitted density
                if (iremove==0) then
                    do ipt=1,npoint_CB+2*npoint
                        if (ipt<=npoint_CB) then !Using second kind Gauss-Chebyshev part of fitting points to check
                            rtmp=radr_CB(ipt)
                        else !Use double dense grids including r=0 to check. The evenly distributed fitting points are subset of this set
                            rtmp=radstep/2*(ipt-npoint_CB-1)
                        end if
                        contrimin=0
                        tmpval=0
                        imin=0
                        do ifunc=1,nfitfunc
                            coeff=parm(ifunc)
                            expon=parm(nfitfunc+ifunc)
                            if (ifunctype==1) then !STO
                                contrival=coeff*exp(-expon*rtmp)
                            else if (ifunctype==2) then !GTF
                                contrival=coeff*exp(-expon*rtmp**2)
                            end if
                            tmpval=tmpval+contrival
                            if (ifunc==1) then
                                imin=1
                                contrimin=contrival
                            else
                                if (contrival<contrimin) then
                                    imin=ifunc
                                    contrimin=contrival
                                end if
                            end if
                        end do
                        if (tmpval<0) then
                            iremove=imin
                            write(*,"(a,f8.4,a)") " Negative fitted density is found at r=",rtmp*b2a," Angstrom"
                            exit
                        end if
                    end do
                end if
                !Check if variation of fitted density is monotonically decrease. If not, the function with largest exponent may be harmful and should be removed
                !Use double dense grid to perform check 
                if (iremove==0) then
                    do ipt=1,npoint
                        rtmp=radstep/2*(ipt-1)
                        rhotmp=calcfitdens(rtmp,nparm,parm(1:nparm))
                        if (ipt==1) then
                            rhoold=rhotmp
                            cycle
                        else
                            if (rhotmp>rhoold) then !Function is increased with increasing r
                                write(*,"(a,f8.4,a)") " Fitted density is not monotonically decreased at r=",rtmp*b2a," Angstrom"
                                iremove=maxloc(parm(nfitfunc+1:nparm),dim=1)
                                exit
                            else
                                rhoold=rhotmp
                            end if
                        end if
                    end do
                end if
                if (iremove==0) then
                    exit
                else
                    write(*,"(' Delete redundant function (coeff=',1PE12.5,' exp=',1PE12.5,'), refitting...')") parm(iremove),parm(nfitfunc+iremove)
                    parm(nfitfunc+iremove:nfitfunc+nfitfunc-1)=parm(nfitfunc+iremove+1:nfitfunc+nfitfunc)
                    parm(iremove:nfitfunc+nfitfunc-2)=parm(iremove+1:nfitfunc+nfitfunc-1)
                    nfitfunc=nfitfunc-1
                    nparm=nfitfunc*2
                    nremove=nremove+1
                end if
            end do
            if (nremove==0) then
                write(*,*) "No redundant fitting functions were found"
            else
                write(*,"(' Totally',i3,' redundant fitting functions have been eliminated')") nremove
            end if
            write(*,*)
        end if
        
        !Output fitting status
        if (info==1.or.info==2.or.info==3) then
            write(*,*) "Fitting has successfully finished!"
        else if (info==5) then
            write(*,"(a,i7)") " Warning: Convergence tolerance has not met while the maximum number of function calls has reached",maxcall
        else if (info==6.or.info==7) then
            write(*,*) "Error: Tolerance is too small, unable to reach the tolerance!"
        end if
        
        !Sort according to exponents from small to large
        if (isort==1) then
            write(*,*) "Sorting fitting functions according to their exponents..."
            forall(i=1:maxfitfunc) seqidx(i)=i
            call sortr8(parm(nfitfunc+1:nparm),list=seqidx(1:nfitfunc))
            parm(1:nfitfunc)=parm(seqidx(1:nfitfunc))
        end if
        
        !Check integral and scale
        rhoint=fitdensint(100,nparm,parm(1:nparm))
        write(*,"(/,' Integral of fitted density calculated using 100 points:',f14.8)") rhoint
        if (iscale==1) then !Scaling fitted coefficients to actual number of electrons
            write(*,"(' Fitted coefficients are scaled by',f16.8)") nelec/rhoint
            parm(1:nfitfunc)=parm(1:nfitfunc)/rhoint*nelec
        end if
        
        !Show final fitted parameters
        write(*,*)
        if (iscale==0) write(*,*) "Fitted parameters (a.u.):"
        if (iscale==1) write(*,*) "Fitted parameters (a.u.) after scaling:"
        write(*,*) "           Coefficient       Exponent"
        do ifunc=1,nfitfunc
            write(*,"(1x,a,i3,':',2(1PE16.6))") funclab(ifunctype),ifunc,parm(ifunc),parm(nfitfunc+ifunc)
        end do
        write(*,*)
        
        !Error statistics
        call calcfitdens_arr(npoint,nparm,parm(1:nparm),fitrho) !Calculate fitted density at fitting points using final (may be scaled) parameters
        write(*,"(' RMSE of fitting error at all points:',f20.6,' a.u.^2')") dsqrt(sum((fitrho-radrho)**2)/npoint)
	    pearsoncoeff=covarray(radrho,fitrho)/stddevarray(radrho)/stddevarray(fitrho)
	    write(*,"(2(a,f12.6))") " Pearson correlation coefficient r:",pearsoncoeff,"  r^2:",pearsoncoeff**2
        if (any(fitrho<0)) write(*,"(/,a)") " Warning: Fitted density at one or more fitting points is negative!"
        
        do while(.true.)
            write(*,*)
            write(*,*) "       -------------------- Quality check & Others --------------------"
            write(*,*) "-2 Export fitted function parameters to fitparm.txt in current folder"
            write(*,*) "-1 Print fitted function parameters again"
            write(*,*) "0 Return"
            write(*,*) "1 Print actual&fitted density and error at fitting points on screen"
            write(*,*) "2 Export actual&fitted density and error at fitting points to radfit.txt"
            write(*,*) "3 Visualize actual density and fitted density curves using logarithmic scaling"
            write(*,*) "4 Visualize actual density and fitted density curves using linear scaling"
            write(*,"(a)") " 5 Export fitted density from 0 to 10 Angstrom with double dense grid to fitdens.txt in current folder"
            write(*,*) "6 Check integral of fitted density"
            write(*,*) "7 Check fitted density at a given radial distance"
            write(*,*) "8 Output coefficients and exponents as Fortran code to a .txt file" !For my personal use
            read(*,*) isel2
            if (isel2==0) then
                exit
            else if (isel2==2) then
                open(10,file="fitparm.txt",status="replace")
                do ifunc=1,nfitfunc
                    write(10,"(2(1PE16.6))") parm(ifunc),parm(nfitfunc+ifunc)
                end do
                close(10)
                write(*,*) "Done! fitparm.txt has been exported in current folder"
                write(*,*) "Column 1: Coefficients (a.u.)"
                write(*,*) "Column 2: Exponents (a.u.)"
            else if (isel2==-1) then
                write(*,*) "Note: Units are in a.u."
                write(*,*) "           Coefficient       Exponent"
                do ifunc=1,nfitfunc
                    write(*,"(1x,a,i3,':',2(1PE16.6))") funclab(ifunctype),ifunc,parm(ifunc),parm(nfitfunc+ifunc)
                end do
            else if (isel2==1) then
                write(*,"(a)") " Radial distance (Angstrom), actual density (a.u.), &
                &difference between fitted and actual density (a.u.) as well as relative difference"
                do ipt=1,npoint
                    write(*,"(' #',i5,'  r:',f8.5,'  rho:',f18.8,'  Diff:',f16.8,' (',f8.2,' %)')") &
                    ipt,radr(ipt)*b2a,radrho(ipt),fitrho(ipt)-radrho(ipt),(fitrho(ipt)-radrho(ipt))/radrho(ipt)*100
                end do
            else if (isel2==2) then
                open(10,file="radfit.txt",status="replace")
                do ipt=1,npoint
                    write(10,"(f8.4,3f20.10)") radr(ipt)*b2a,radrho(ipt),fitrho(ipt),fitrho(ipt)-radrho(ipt)
                end do
                close(10)
                write(*,*)
                write(*,*) "Data has been exported to radfit.txt in current folder. Content:"
                write(*,*) "Column 1: Radial distance (Angstrom)"
                write(*,*) "Column 2: Actual density (a.u.)"
                write(*,*) "Column 3: Fitted density (a.u.)"
                write(*,*) "Column 4: Fitting error (a.u.)"
            else if (isel2==3.or.isel2==4) then
				call METAFL('xwin')
				call window(100,0,1200,720)
				call SCRMOD('REVERSE')
				CALL PAGE(3000,1800)
				call disini
				call height(40)
				CALL HNAME(45)
				call hwfont
				call AXSLEN(2450,1400)
				call WINTIT("Click right mouse button to close")
				call ERRMOD("ALL","OFF")
				call AXSPOS(380,1550)
                !Set axis style
                CALL NAMDIS(40,'Y')
				CALL NAME('Radial distance (Angstrom)','X')
				CALL NAME('Electron density (a.u.)','Y')
				CALL TICKS(1,'XY') !1 tick between two labels
                CALL LABDIG(1,"X")
                if (isel2==3) then
                    call AXSSCL('log','Y')
                    Ymin=-5 !Start from 1E-5
                    Ymax=ceiling(log10(maxval(radrho)))
                    Ystep=1
                    call labels('log','Y')
				    CALL LABDIG(-1,"Y") !Do not show digit in Y
				    CALL GRAF(0D0,4D0,0D0,0.2D0, Ymin,Ymax,Ymin,Ystep) !Plot r=0~4 Angstrom
                else if (isel2==4) then
                    call AXSSCL('lin','Y')
                    Ymin=0
                    Ymax=min(max(maxval(radrho),maxval(fitrho)),2D0)
                    Ystep=(Ymax-Ymin)/10
                    call labels('float','Y')
				    CALL LABDIG(2,"Y")
				    CALL GRAF(0D0,2D0,0D0,0.2D0, Ymin,Ymax,Ymin,Ystep) !Plot r=0~2 Angstrom
                end if
                call SETRGB(0.8D0,0.8D0,0.8D0) !Shallow gray grid
			    call dash
			    call LINWID(1)
                CALL GRID(1,1)
                call solid
                !Draw legends and curves
                call legini(clegend,2,40)
				call legtit(' ') !Do now show legend title
				call frame(0) !No box around legend
				call linwid(5) !Use thick line
                call setcolor(5) !Black to plot actual density
                CALL LEGLIN(clegend,"Actual density",1)
				CALL CURVE(radr*b2a,radrho,npoint)
                call setcolor(3) !Blue to plot fitted density
                CALL LEGLIN(clegend,"Fitted density",2)
				CALL CURVE(radr*b2a,fitrho,npoint)
                call setcolor(5) !Black legend text
                call legend(clegend,7)
				call disfin
            else if (isel2==5) then
                open(10,file="fitdens.txt",status="replace")
                rtmp=0
                do while(.true.)
                    write(10,"(f8.4,f35.12)") rtmp*b2a,calcfitdens(rtmp,nparm,parm)
                    rtmp=rtmp+radstep/2
                    if (rtmp*b2a>10) exit
                end do
                close(10)
                write(*,*) "Done! fitdens.txt has been exported in current folder"
                write(*,*) "Column 1: Radial distance (Angstrom)"
                write(*,*) "Column 2: Fitted density (a.u.)"
            else if (isel2==6) then
                do ntmp=40,300,20
                    rhoint=fitdensint(ntmp,nparm,parm(1:nparm))
                    write(*,"(' Number of integration points:',i5,'    Integral:',f16.8)") ntmp,rhoint
                end do
            else if (isel2==7) then
                write(*,*) "Input radial distance in Angstrom, e.g. 3.8"
                read(*,*) rtmp
                write(*,"(' Fitted density is',1PE16.8,' a.u.')") calcfitdens(rtmp/b2a,nparm,parm)
            else if (isel2==8) then
                write(c80tmp,"(i3.3,'_',a,'.txt')") a(1)%index,trim(ind2name(a(1)%index))
                open(10,file=trim(c80tmp),status="replace")
                write(10,"(a,i3,a)") "case(",a(1)%index,")"
                if (ifunctype==1) write(10,"(a,i2)") "    nSTO=",nfitfunc
                if (ifunctype==2) write(10,"(a,i2)") "    nGTF=",nfitfunc
                do ifunc=1,nfitfunc
                    write(10,"('    atomcoeff(',i2,')=',1PD16.8)") ifunc,parm(ifunc)
                end do
                do ifunc=1,nfitfunc
                    write(10,"('    atomexp(',i2,')=',1PD16.8)") ifunc,parm(nfitfunc+ifunc)
                end do
                close(10)
                write(*,*) "Data has been outputted to "//trim(c80tmp)
            end if
        end do
        
        deallocate(radr,radrho,fiterr,fitrho)
        if (npoint_CB>0) deallocate(radr_CB,radw_CB,rho_CB)
    end if
end do

end subroutine


!!---- Input coefficients and exponents, return array containing difference between actual density and fitted density at positions in radr(:) array
!radr, radrho, ifunctype and ifittype in fitatmdens_mod are involved
!If iflag=-999, the returned "diff" will correspond to fitted density
!The arguments are required by lmdif1 routine in MINPACK
subroutine atmdens_fiterr(npoint,nparm,parm,diff,iflag)
use fitatmdens_mod
implicit real*8 (a-h,o-z)
integer npoint,nparm,iflag
real*8 :: parm(nparm),diff(npoint)
do ipt=1,npoint
    tmpval=0
    do ifunc=1,nfitfunc
        coeff=parm(ifunc)
        if (ifixexp==0) then
            expon=parm(nfitfunc+ifunc)
        else
            expon=exp_fit(ifunc)
        end if
        if (ifunctype==1) then !STO
            tmpval=tmpval+coeff*exp(-expon*radr(ipt))
        else if (ifunctype==2) then !GTF
            tmpval=tmpval+coeff*exp(-expon*radr(ipt)**2)
        end if
    end do
    if (iflag==-999) then
        diff(ipt)=tmpval
    else
        if (ifittype==1) then
            diff(ipt)=tmpval-radrho(ipt)
        else if (ifittype==2) then
            diff(ipt)=abs(tmpval-radrho(ipt))/radrho(ipt)
        else if (ifittype==3) then
            diff(ipt)=(tmpval-radrho(ipt))*radr(ipt)**2
        end if
    end if
end do
end subroutine

!!---- Input coefficients and exponents, return the array containing fitted density at positions in radr(:)
subroutine calcfitdens_arr(npoint,nparm,parm,fitrho)
integer npoint,nparm
real*8 parm(nparm),fitrho(npoint)
call atmdens_fiterr(npoint,nparm,parm,fitrho,-999)
end subroutine

!!---- Input coefficients and exponents, return value of fitted density at a given radial distance
real*8 function calcfitdens(r,nparm,parm)
use fitatmdens_mod
implicit real*8 (a-h,o-z)
real*8 parm(nparm),r
calcfitdens=0
do ifunc=1,nfitfunc
    if (ifunctype==1) then !STO
        calcfitdens=calcfitdens+parm(ifunc)*exp(-parm(nfitfunc+ifunc)*r)
    else if (ifunctype==2) then !GTF
        calcfitdens=calcfitdens+parm(ifunc)*exp(-parm(nfitfunc+ifunc)*r**2)
    end if
end do
end function

!!---- Return integral of fitted density based on input number of integration points and parameters
!ifunctype in fitatmdens_mod is involved
real*8 function fitdensint(nradpt,nparm,parm)
use defvar
use fitatmdens_mod
implicit real*8 (a-h,o-z)
integer nradpt,nparm
real*8 parm(nparm),radr_int(nradpt),radw_int(nradpt)
parmbk=1D0
do i=1,nradpt !Combine spherical point&weights with second kind Gauss-Chebyshev method for radial part
	radx=cos(i*pi/(nradpt+1))
	radr_int(i)=(1+radx)/(1-radx)*parmbk
	radw_int(i)=2*pi/(nradpt+1)*parmbk**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
end do
fitdensint=0
do i=nradpt,1,-1 !From close to far
    rhotmp=0
    do ifunc=1,nfitfunc
        if (ifunctype==1) then !STO
            rhotmp=rhotmp+parm(ifunc)*exp(-parm(nfitfunc+ifunc)*radr_int(i))
        else if (ifunctype==2) then !GTF
            rhotmp=rhotmp+parm(ifunc)*exp(-parm(nfitfunc+ifunc)*radr_int(i)**2)
        end if
    end do
    fitdensint=fitdensint+rhotmp*radw_int(i)
end do
end function





!!----------- Simulating scanning tunneling microscope (STM) image
!Main ref: https://en.wikipedia.org/wiki/Scanning_tunneling_microscope
!Partial ref: Tersoff and Hamann, Theory of the scanning tunneling microscope, PRB, 31, 805 (1985)
subroutine STM
use defvar
use GUI
use functions
use util
implicit real*8 (a-h,o-z)
integer :: imode=2
real*8 :: bias=0
character c80tmp*80
real*8,external :: LDOS_STM

if (allocated(cubmat)) deallocate(cubmat)
nx=200;ny=200
!Set initial range, in Bohr
if (ifPBC==0) then
    orgx=minval(a%x)-3
    endx=maxval(a%x)+3
    orgy=minval(a%y)-3
    endy=maxval(a%y)+3
else
    call cellmaxxyz(endx,endy,zmax)
    call cellminxyz(orgx,orgy,zmin)
end if
orgz=maxval(a%z)+0.7D0/b2a !Scan Z=0.7~2.5 Angstrom with respect to top atom
endz=orgz+1.8D0/b2a

if (.not.allocated(b)) then
    write(*,"(a)") " Error: In order to use this function, the input file must at least contain GTF information! See Section 2.5 of manual for detail."
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

if (wfntype==3.or.wfntype==4) then
    write(*,"(a)") " Error: This function does not formally support wavefunction with non-integer orbital occupancy!"
    write(*,"(a)") " Please enter subfunction 9 of main function 300 to determine Fermi level, which can be used in the present function, and &
    &at the same time, all orbital occupancies will be set to integer"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
else if (wfntype==2) then
    write(*,"(a)") " Error: Restricted open-shell wavefunction is not directly supported by this function. You should first use subfunction 37 in &
    &main function 6 to transform the wavefunction to equivalent unrestricted open-shell form!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
else !Single-determinant wavefunction
    if (allocated(CObasa)) then
        call getHOMOidx
        Ef=(MOene(idxHOMO)+MOene(idxHOMO+1))/2
        bias=MOene(idxHOMO)-Ef
    
        if (wfntype==0) then
            write(*,*) "Note: The default Fermi level has been set to average of E(HOMO) and E(LUMO)"
        else
            write(*,*) "Note: The default Fermi level has been set to average of E(HOMO) and E(LUMO) of alpha spin. In this case, &
            &the result will be problematic if bias voltage is set to positive value (electron flows from tip to sample)"
        end if
        write(*,"(a)") " The default bias voltage has been set to the difference between E(HOMO) and Fermi level, &
        &therefore under default setting only HOMO will be imaged"
    else
        Ef=maxval(MOene)
        write(*,*) "Note: The default Fermi level has been set to HOMO"
        write(*,"(a)") " Note: Since there is no unoccupied MO, the bias voltage must be set to negative value (electron flows from sample to tip)"
    end if
end if

do while(.true.)
    write(*,*)
    write(*,*) " ----------- Simulating scanning tunneling microscope (STM) image -----------"
    write(*,*) "-1 Return"
    if (imode==1) write(*,*) "0 Calculating grid data of tunneling current!"
    if (imode==2) write(*,*) "0 Calculating tunneling current on the plane!"
    if (imode==1) write(*,*) "1 Toggle mode of STM image, current: Constant current"
    if (imode==2) write(*,*) "1 Toggle mode of STM image, current: Constant distance"
    write(*,"(a,f10.3,' V')") " 2 Set bias voltage, current:",bias*au2eV
    write(*,"(a,f10.3,' eV')") " 3 Set Fermi level, current:",Ef*au2eV
    !write(*,"(a,f6.3,' eV')") " 3 Set FWHM for Gaussian broadening, current:",STM_FWHM
    if (imode==1) write(*,"(a,3i5)") " 4 Set number of grid points in X,Y,Z, current:",nx,ny,nz
    if (imode==2) write(*,"(a,3i5)") " 4 Set number of grid points in X and Y, current:",nx,ny
    write(*,"(a,f8.3,a,f8.3,a)") " 5 Set range in X direction, current: From",orgx*b2a," to",endx*b2a," Angstrom"
    write(*,"(a,f8.3,a,f8.3,a)") " 6 Set range in Y direction, current: From",orgy*b2a," to",endy*b2a," Angstrom"
    if (imode==1) write(*,"(a,f8.3,a,f8.3,a)") " 7 Set range in Z direction, current: From",orgz*b2a," to",endz*b2a," Angstrom"
    if (imode==2) write(*,"(a,f8.3,a)") " 7 Set Z coordinate of the XY plane, current:",orgz*b2a," Angstrom"
    read(*,*) isel
    
    if (isel==-1) then
        return
    else if (isel==1) then
        if (imode==1) then
            imode=2
            nx=200;ny=200
        else if (imode==2) then
            imode=1
            nx=150;ny=150;nz=80
        end if
    else if (isel==2) then
        write(*,*) "Input bias voltage in V, e.g. -3.5"
        write(*,"(a)") " Note: Negative value lets electron flow from sample to tip, thus density of occupied MOs are imaged. &
        &Positive value lets electron flow from tip to sample, thus unoccupied MOs are imaged."
        read(*,*) bias
        bias=bias/au2eV
    else if (isel==3) then
        write(*,*) "Input Fermi energy in eV, e.g. -5.82"
        read(*,*) Ef
        Ef=Ef/au2eV
    else if (isel==4) then
        if (imode==1) then
            write(*,*) "1 Coarse grid (100*100*40)"
            write(*,*) "2 Medium grid (150*150*70)"
            write(*,*) "3 Fine grid (200*200*100)"
            write(*,*) "or, directly input number of grid points in X,Y,Z, e.g. 80,80,30"
            read(*,"(a)") c80tmp
            read(c80tmp,*,iostat=ierror) nx,ny,nz
            if (ierror/=0) then
                read(c80tmp,*) isel2
                if (isel2==1) then
                    nx=100;ny=100;nz=40
                else if (isel2==2) then
                    nx=150;ny=150;nz=70
                else if (isel2==3) then
                    nx=200;ny=200;nz=100
                end if
            end if
        else if (imode==2) then
            write(*,*) "Input number of grid points in X and Y, e.g. 80,80"
            read(*,*) nx,ny
        end if
    else if (isel==5) then
        write(*,*) "Input lower and upper limit of X in Angstrom, e.g. -5.8,6.4"
        write(*,"(a)") " If you only input a number, it will be employed as extension distance (in Angstrom) in X direction to properly determine the X range"
        read(*,"(a)") c80tmp
        read(c80tmp,*,iostat=ierror) orgx,endx
        if (ierror==0) then
            orgx=orgx/b2a
            endx=endx/b2a
        else
            read(c80tmp,*) ext
            orgx=minval(a%x)-ext/b2a
            endx=maxval(a%x)+ext/b2a
        end if
    else if (isel==6) then
        write(*,*) "Input lower and upper limit of Y in Angstrom, e.g. -5.8,6.4"
        write(*,"(a)") " If you only input a number, it will be employed as extension distance (in Angstrom) in Y direction to properly determine the X range"
        read(*,"(a)") c80tmp
        read(c80tmp,*,iostat=ierror) orgy,endy
        if (ierror==0) then
            orgy=orgy/b2a
            endy=endy/b2a
        else
            read(c80tmp,*) ext
            orgy=minval(a%y)-ext/b2a
            endy=maxval(a%y)+ext/b2a
        end if
    else if (isel==7) then
        if (imode==1) then
            write(*,*) "Input lower and upper limit of Z in Angstrom, e.g. 0,2.5"
            read(*,*) orgz,endz
            orgz=orgz/b2a
            endz=endz/b2a
        else if (imode==2) then
            write(*,*) "Input Z coordinate of the XY plane in Angstrom, e.g. 2.2"
            read(*,*) orgz
            orgz=orgz/b2a
        end if
    
    else if (isel==0) then !Start calculation !!!
        
        !Show which MOs will be taken into account
        if (bias<=0) then
            Elow=Ef+bias
            Ehigh=Ef
        else
            Elow=Ef
            Ehigh=Ef+bias
        end if
        write(*,"(/,' Lower limit of MO energy considered in the calculation:',f12.3,' eV')") Elow*au2eV
        write(*,"(' Upper limit of MO energy considered in the calculation:',f12.3,' eV')") Ehigh*au2eV
        write(*,*) "The MOs taken into account in the current STM simulation:"
        ialphabeg=0 !Record range of alpha (or spatial) MOs to be considered
        ialphaend=0
        ibetabeg=0 !Record range of beta MOs to be considered
        ibetaend=0
        nconsider=0
        do imo=1,nmo
            iadd=0
            if (bias<=0) then !Electron flows from sample to tip
                if (MOocc(imo)/=0.and.MOene(imo)>=Elow.and.MOene(imo)<=Ehigh) then
                    write(*,"(' MO',i6,'   Occ=',f6.3,'   Energy=',f12.4,' eV   Type: ',a)") imo,MOocc(imo),MOene(imo)*au2eV,trim(orbtypename(MOtype(imo)))
                    iadd=1
                end if
            else if (bias>0) then !Electron flows from tip to sample
                if (MOocc(imo)==0.and.MOene(imo)>=Elow.and.MOene(imo)<=Ehigh) then
                    write(*,"(' MO',i6,'   Occ=',f6.3,'   Energy=',f12.4,' eV   Type: ',a)") imo,MOocc(imo),MOene(imo)*au2eV,trim(orbtypename(MOtype(imo)))
                    iadd=1
                end if
            end if
            if (iadd==1) then
                if (MOtype(imo)<=1) then
                    if (ialphabeg==0) ialphabeg=imo
                    ialphaend=imo
                else
                    if (ibetabeg==0) ibetabeg=imo
                    ibetaend=imo
                end if
                nconsider=nconsider+1
            end if
        end do
        if (nconsider==0) then
            write(*,*) "None. Therefore the calculation is canceled"
            cycle
        else
            if (wfntype==0) then
                write(*,"(' Range of MOs to be taken into account:',2i8)") ialphabeg,ialphaend
            else
                write(*,"(' Range of alpha MOs to be taken into account:',2i8)") ialphabeg,ialphaend
                write(*,"(' Range of beta MOs to be taken into account: ',2i8)") ibetabeg,ibetaend
            end if
            write(*,"(' Totally',i5,' MOs are taken into account',/)") nconsider
        end if
        
        dx=(endx-orgx)/(nx-1)
        dy=(endy-orgy)/(ny-1)
        !Prepare settings for plane plot, they are utilized by "drawplane" routine via "planemap_interface" routine
        call gencontour(0,0D0,0D0,0) !Generate contour lines
        ngridnum1=nx
        ngridnum2=ny
        plesel=1 !XY plane
        disshowlabel=100 !Very broad threshold to make sure showing all atom labels
        ilenunit2D=2 !Use Angstrom
        planestpx=(endx-orgx)*b2a/7
        planestpy=(endy-orgy)*b2a/7
        iclrtrans=6 !Grey transition in color-filled map
        iatom_on_plane=1 !Show atom label on map
        numdigz=4
        
        if (imode==1) then !Constant current STM
            dz=(endz-orgz)/(nz-1)
            gridv1=(/ dx,0D0,0D0 /)
            gridv2=(/ 0D0,dy,0D0 /)
            gridv3=(/ 0D0,0D0,dz /)
	        write(*,"(' Grid spacings in X,Y,Z are',3f12.6,' Bohr')") dx,dy,dz
            if (ifPBC==0) then
	            call gen_GTFuniq(0) !Generate unique GTFs, for faster evaluation in orbderv
            else
	            call gen_neigh_GTF !Generate neighbouring GTFs list at reduced grids, for faster evaluation
            end if
            write(*,*) "Calculating, please wait..."
            allocate(cubmat(nx,ny,nz))
            call walltime(iwalltime1)
            ifinish=0;ishowprog=1
            ntmp=floor(ny*nz/100D0)
            !$OMP PARALLEL DO SHARED(cubmat,ifinish) PRIVATE(ix,xpos,iy,ypos,iz,zpos) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
            do iz=1,nz
                do iy=1,ny
                    do ix=1,nx
                        call getgridxyz(ix,iy,iz,xpos,ypos,zpos)
                        cubmat(ix,iy,iz)=LDOS_STM(xpos,ypos,zpos,ialphabeg,ialphaend)
                        if (wfntype==0) then
                            cubmat(ix,iy,iz)=cubmat(ix,iy,iz)*2
                        else if (wfntype==1.and.ibetaend/=0) then !Beta part
                            cubmat(ix,iy,iz)=cubmat(ix,iy,iz)+LDOS_STM(xpos,ypos,zpos,ibetabeg,ibetaend)
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
            !$OMP END PARALLEL DO
            if (ishowprog/=0) call showprog(100,100)
            call walltime(iwalltime2)
            write(*,"(' Calculation took up time',i10,' s')") iwalltime2-iwalltime1
            valmax=maxval(cubmat)
            write(*,"(' Maximal value (LDOS) is',f12.6,' a.u.')") valmax
            sur_value=valmax/2 !For isosurface plot
            
            do while(.true.)
                write(*,*)
                write(*,*) "          ------------------- Post-processing menu -------------------"
                write(*,*) "0 Return"
                write(*,*) "1 Visualize isosurface of current"
                write(*,*) "2 Export grid data of current as STM.cub in current folder"
                write(*,*) "3 Calculate and visualize constant current STM image"
                read(*,*) isel2
                if (isel2==0) then
                    deallocate(cubmat)
                    exit
                else if (isel2==1) then
                    call drawisosurgui(1)
                else if (isel2==2) then
                    open(10,file="STM.cub",status="replace")
			        call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
                    close(10)
                    write(*,*) "Exporting finished!"
                else if (isel2==3) then
                    write(*,*) "Input constant current value, e.g. 0.004"
                    write(*,"(' Note: The value should be larger than 0 and smaller than',f10.6)") valmax
                    read(*,*) constcurr
                    write(*,*) "Calculating constant current map, please wait..."
                    allocate(planemat(nx,ny))
                    ifinish=0
                    do iy=1,ny
                        do ix=1,nx
                            !Gradually decrease tip from high z to low z, until isosurface is encountered, and using linear interpolation to determine its position
                            do iz=nz,2,-1
                                call getgridxyz(ix,iy,iz,xpos,ypos,zpos)
                                zposnext=orgz+(iz-2)*dz
                                val=cubmat(ix,iy,iz)
                                valnext=cubmat(ix,iy,iz-1)
                                if (val<constcurr.and.valnext>constcurr) then
                                    zdiffmin=zposnext+(constcurr-valnext)/(val-valnext)*dz
                                    exit
                                end if
                            end do
                            if (iz==1) zdiffmin=orgz !Failed to determine isosurface position
                            planemat(ix,iy)=zdiffmin*b2a
                        end do
                    end do
                    
                    write(*,"(' Minimal Z is',f12.6,' Angstrom')") minval(planemat)
                    write(*,"(' Maximal Z is',f12.6,' Angstrom')") maxval(planemat)
                    clrlow=minval(planemat)*0.99999D0;clrhigh=maxval(planemat)*1.00001D0 !Avoid a few points marginally exceed upper limit
                    planestpz=(clrhigh-clrlow)/10
                    orgz2D=endz !In fact this is meaningless for present case, but can be used to determine distance for plotting atomic labels, namely distance between atom and maximum z of calculated region
                    idrawtype=1
                    idrawcontour=0
                    call gencontour(2,clrlow,clrhigh,10) !Generate contour lines evenly covering lower and upper limits
                    call planemap_interface("constant current STM","STM",orgx,endx,orgy,endy,clrlow,clrhigh)
                    deallocate(planemat)
                end if
            end do
        
        else if (imode==2) then !Constant height STM
	        write(*,"(' Grid spacings in X and Y are',2f12.6,' Bohr')") dx,dy
            write(*,*) "Calculating, please wait..."
            allocate(planemat(nx,ny))
            !$OMP PARALLEL DO SHARED(planemat) PRIVATE(ix,xpos,iy,ypos) schedule(dynamic) NUM_THREADS(nthreads)
            do iy=1,ny
                do ix=1,nx
                    xpos=orgx+(ix-1)*dx
                    ypos=orgy+(iy-1)*dy
                    planemat(ix,iy)=LDOS_STM(xpos,ypos,orgz,ialphabeg,ialphaend)
                    if (wfntype==0) then
                        planemat(ix,iy)=planemat(ix,iy)*2
                    else if (wfntype==1.and.ibetaend/=0) then !Beta part
                        planemat(ix,iy)=planemat(ix,iy)+LDOS_STM(xpos,ypos,orgz,ibetabeg,ibetaend)
                    end if
                end do
            end do
            !$OMP END PARALLEL DO
            write(*,"(' Maximal value (LDOS) is',f12.6,' a.u.')") maxval(planemat)
            clrlow=0D0
            clrhigh=maxval(planemat)
            planestpz=(clrhigh-clrlow)/10
            orgz2D=orgz
            idrawtype=1
            idrawcontour=0
            call planemap_interface("constant height STM","STM",orgx,endx,orgy,endy,clrlow,clrhigh)
            deallocate(planemat)
        end if
        
    end if
end do

end subroutine


!!-------- Return LDOS at x,y,z (Bohr) contributed by MO from ibeg to iend. Invoked by subroutine STM
real*8 function LDOS_STM(x,y,z,ibeg,iend)
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 x,y,z,wfnval(nmo)
integer ibeg,iend

call orbderv(1,ibeg,iend,x,y,z,wfnval)
LDOS_STM=0
do imo=ibeg,iend
    LDOS_STM=LDOS_STM+wfnval(imo)**2
end do
end function




!!--------- Calculate electric dipole/multipole moments and electronic spatial extent based on analytic integrals
subroutine calc_multipole
use defvar
use util
implicit real*8 (a-h,o-z)

!Calculate based on atomic charges, and print dipole moment information
if (ifiletype==4) then
    write(*,"(a)") " Note: Wavefunction information is not available, the data will be evaluated based on atomic charges"
    xdip=0
    ydip=0
    zdip=0
    do iatm=1,ncenter
        xdip=xdip+a(iatm)%x*a(iatm)%charge
        ydip=ydip+a(iatm)%y*a(iatm)%charge
        zdip=zdip+a(iatm)%z*a(iatm)%charge
    end do
    write(*,"(/,' Dipole moment (a.u.): ',3f14.6)") xdip,ydip,zdip
    write(*,"(' Dipole moment (Debye):',3f14.6)") xdip*au2debye,ydip*au2debye,zdip*au2debye
    dipmag=dsqrt(xdip**2+ydip**2+zdip**2)
    write(*,"(' Magnitude of dipole moment:',f14.6,' a.u.',f14.6,' Debye')") dipmag,dipmag*au2debye
    
    xxinttot=0;yyinttot=0;zzinttot=0;xyinttot=0;yzinttot=0;xzinttot=0
    xxxinttot=0;yyyinttot=0;zzzinttot=0;yzzinttot=0
    xzzinttot=0;xxzinttot=0;yyzinttot=0;xxyinttot=0;xyyinttot=0;xyzinttot=0
    xxxxinttot=0;yyyyinttot=0;zzzzinttot=0;xxxyinttot=0;xxxzinttot=0
    yyyxinttot=0;yyyzinttot=0;zzzxinttot=0;zzzyinttot=0;xxyyinttot=0
    xxzzinttot=0;yyzzinttot=0;xxyzinttot=0;yyxzinttot=0;zzxyinttot=0
    do iatm=1,ncenter
        xxinttot=xxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%charge
        yyinttot=yyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%charge
        zzinttot=zzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%charge
        xyinttot=xyinttot+a(iatm)%x*a(iatm)%y*a(iatm)%charge
        yzinttot=yzinttot+a(iatm)%y*a(iatm)%z*a(iatm)%charge
        xzinttot=xzinttot+a(iatm)%x*a(iatm)%z*a(iatm)%charge
	    xxxinttot=xxxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%charge
	    yyyinttot=yyyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%charge
	    zzzinttot=zzzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%charge
	    yzzinttot=yzzinttot+a(iatm)%y*a(iatm)%z*a(iatm)%z*a(iatm)%charge
	    xzzinttot=xzzinttot+a(iatm)%x*a(iatm)%z*a(iatm)%z*a(iatm)%charge
	    xxzinttot=xxzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%z*a(iatm)%charge
	    yyzinttot=yyzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%z*a(iatm)%charge
	    xxyinttot=xxyinttot+a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%charge
	    xyyinttot=xyyinttot+a(iatm)%x*a(iatm)%y*a(iatm)%y*a(iatm)%charge
	    xyzinttot=xyzinttot+a(iatm)%x*a(iatm)%y*a(iatm)%z*a(iatm)%charge
        xxxxinttot=xxxxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%charge
        yyyyinttot=yyyyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%charge
        zzzzinttot=zzzzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%charge
        xxxyinttot=xxxyinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%charge
        xxxzinttot=xxxzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%z*a(iatm)%charge
        yyyxinttot=yyyxinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%x*a(iatm)%charge
        yyyzinttot=yyyzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%z*a(iatm)%charge
        zzzxinttot=zzzxinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%x*a(iatm)%charge
        zzzyinttot=zzzyinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%y*a(iatm)%charge
        xxyyinttot=xxyyinttot+a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%y*a(iatm)%charge
        xxzzinttot=xxzzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%z*a(iatm)%z*a(iatm)%charge
        yyzzinttot=yyzzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%z*a(iatm)%z*a(iatm)%charge
        xxyzinttot=xxyzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%z*a(iatm)%charge
        yyxzinttot=yyxzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%x*a(iatm)%z*a(iatm)%charge
        zzxyinttot=zzxyinttot+a(iatm)%z*a(iatm)%z*a(iatm)%x*a(iatm)%y*a(iatm)%charge
    end do

else !Calculation based on wavefunction information, and print dipole moment information
    if (ispecial==1) then
        xnucdip=0
        ynucdip=0
        znucdip=0
        do iatm=1,ncenter
            xnucdip=xnucdip+a(iatm)%x*a(iatm)%charge
            ynucdip=ynucdip+a(iatm)%y*a(iatm)%charge
            znucdip=znucdip+a(iatm)%z*a(iatm)%charge
        end do
        write(*,"(/,' Dipole moment from nuclear charges (a.u.): ',3f11.6)") xnucdip,ynucdip,znucdip
        write(*,"(a)") " Because ispecial=1, now displacing nuclear coordinates to make their contributions to dipole moment vanishing"
        sumnuc=sum(a%charge)
        do iatm=1,ncenter
            a(iatm)%x=a(iatm)%x-xnucdip/sumnuc
            a(iatm)%y=a(iatm)%y-ynucdip/sumnuc
            a(iatm)%z=a(iatm)%z-znucdip/sumnuc
        end do
        write(*,*) "Done!"
    end if

    if (allocated(CObasa)) then
        write(*,"(a)") " Calculating electric dipole, quadruple, octopole and Hexadecapole moment integral matrix..."
        call genMultipolebas_curr

        xinttot=sum(Dbas(1,:,:)*Ptot(:,:))
        yinttot=sum(Dbas(2,:,:)*Ptot(:,:))
        zinttot=sum(Dbas(3,:,:)*Ptot(:,:))

        xxinttot=sum(Quadbas(1,:,:)*Ptot(:,:))
        yyinttot=sum(Quadbas(2,:,:)*Ptot(:,:))
        zzinttot=sum(Quadbas(3,:,:)*Ptot(:,:))
        xyinttot=sum(Quadbas(4,:,:)*Ptot(:,:))
        yzinttot=sum(Quadbas(5,:,:)*Ptot(:,:))
        xzinttot=sum(Quadbas(6,:,:)*Ptot(:,:))

        xxxinttot=sum(Octobas(1,:,:)*Ptot(:,:))
        yyyinttot=sum(Octobas(2,:,:)*Ptot(:,:))
        zzzinttot=sum(Octobas(3,:,:)*Ptot(:,:))
        yzzinttot=sum(Octobas(4,:,:)*Ptot(:,:))
        xzzinttot=sum(Octobas(5,:,:)*Ptot(:,:))
        xxzinttot=sum(Octobas(6,:,:)*Ptot(:,:))
        yyzinttot=sum(Octobas(7,:,:)*Ptot(:,:))
        xxyinttot=sum(Octobas(8,:,:)*Ptot(:,:))
        xyyinttot=sum(Octobas(9,:,:)*Ptot(:,:))
        xyzinttot=sum(Octobas(10,:,:)*Ptot(:,:))
    
        xxxxinttot=sum(Hexdebas(1,:,:)*Ptot(:,:))
        yyyyinttot=sum(Hexdebas(2,:,:)*Ptot(:,:))
        zzzzinttot=sum(Hexdebas(3,:,:)*Ptot(:,:))
        xxxyinttot=sum(Hexdebas(4,:,:)*Ptot(:,:))
        xxxzinttot=sum(Hexdebas(5,:,:)*Ptot(:,:))
        yyyxinttot=sum(Hexdebas(6,:,:)*Ptot(:,:))
        yyyzinttot=sum(Hexdebas(7,:,:)*Ptot(:,:))
        zzzxinttot=sum(Hexdebas(8,:,:)*Ptot(:,:))
        zzzyinttot=sum(Hexdebas(9,:,:)*Ptot(:,:))
        xxyyinttot=sum(Hexdebas(10,:,:)*Ptot(:,:))
        xxzzinttot=sum(Hexdebas(11,:,:)*Ptot(:,:))
        yyzzinttot=sum(Hexdebas(12,:,:)*Ptot(:,:))
        xxyzinttot=sum(Hexdebas(13,:,:)*Ptot(:,:))
        yyxzinttot=sum(Hexdebas(14,:,:)*Ptot(:,:))
        zzxyinttot=sum(Hexdebas(15,:,:)*Ptot(:,:))
    
    else if (allocated(b)) then
        write(*,*) "Calculating density matrix based on GTFs..."
        call genPprim
        write(*,"(a)") " Calculating electric dipole, quadruple, octopole and Hexadecapole moment integral matrix..."
        call genMultipoleprim

        xinttot=sum(Dprim(1,:,:)*Ptot_prim(:,:))
        yinttot=sum(Dprim(2,:,:)*Ptot_prim(:,:))
        zinttot=sum(Dprim(3,:,:)*Ptot_prim(:,:))

        xxinttot=sum(Quadprim(1,:,:)*Ptot_prim(:,:))
        yyinttot=sum(Quadprim(2,:,:)*Ptot_prim(:,:))
        zzinttot=sum(Quadprim(3,:,:)*Ptot_prim(:,:))
        xyinttot=sum(Quadprim(4,:,:)*Ptot_prim(:,:))
        yzinttot=sum(Quadprim(5,:,:)*Ptot_prim(:,:))
        xzinttot=sum(Quadprim(6,:,:)*Ptot_prim(:,:))

        xxxinttot=sum(Octoprim(1,:,:)*Ptot_prim(:,:))
        yyyinttot=sum(Octoprim(2,:,:)*Ptot_prim(:,:))
        zzzinttot=sum(Octoprim(3,:,:)*Ptot_prim(:,:))
        yzzinttot=sum(Octoprim(4,:,:)*Ptot_prim(:,:))
        xzzinttot=sum(Octoprim(5,:,:)*Ptot_prim(:,:))
        xxzinttot=sum(Octoprim(6,:,:)*Ptot_prim(:,:))
        yyzinttot=sum(Octoprim(7,:,:)*Ptot_prim(:,:))
        xxyinttot=sum(Octoprim(8,:,:)*Ptot_prim(:,:))
        xyyinttot=sum(Octoprim(9,:,:)*Ptot_prim(:,:))
        xyzinttot=sum(Octoprim(10,:,:)*Ptot_prim(:,:))

        xxxxinttot=sum(Hexdeprim(1,:,:)*Ptot_prim(:,:))
        yyyyinttot=sum(Hexdeprim(2,:,:)*Ptot_prim(:,:))
        zzzzinttot=sum(Hexdeprim(3,:,:)*Ptot_prim(:,:))
        xxxyinttot=sum(Hexdeprim(4,:,:)*Ptot_prim(:,:))
        xxxzinttot=sum(Hexdeprim(5,:,:)*Ptot_prim(:,:))
        yyyxinttot=sum(Hexdeprim(6,:,:)*Ptot_prim(:,:))
        yyyzinttot=sum(Hexdeprim(7,:,:)*Ptot_prim(:,:))
        zzzxinttot=sum(Hexdeprim(8,:,:)*Ptot_prim(:,:))
        zzzyinttot=sum(Hexdeprim(9,:,:)*Ptot_prim(:,:))
        xxyyinttot=sum(Hexdeprim(10,:,:)*Ptot_prim(:,:))
        xxzzinttot=sum(Hexdeprim(11,:,:)*Ptot_prim(:,:))
        yyzzinttot=sum(Hexdeprim(12,:,:)*Ptot_prim(:,:))
        xxyzinttot=sum(Hexdeprim(13,:,:)*Ptot_prim(:,:))
        yyxzinttot=sum(Hexdeprim(14,:,:)*Ptot_prim(:,:))
        zzxyinttot=sum(Hexdeprim(15,:,:)*Ptot_prim(:,:))
    else
        write(*,*) "Error: The current input file contain neither wavefunction information nor atomic charges, this function cannot be used!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if

    ESEx=-xxinttot
    ESEy=-yyinttot
    ESEz=-zzinttot

    !Combine nuclear contribution and electron contribution to obtain multiple moments
    xnucdip=0
    ynucdip=0
    znucdip=0
    do iatm=1,ncenter
        xnucdip=xnucdip+a(iatm)%x*a(iatm)%charge
        ynucdip=ynucdip+a(iatm)%y*a(iatm)%charge
        znucdip=znucdip+a(iatm)%z*a(iatm)%charge
        xxinttot=xxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%charge
        yyinttot=yyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%charge
        zzinttot=zzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%charge
        xyinttot=xyinttot+a(iatm)%x*a(iatm)%y*a(iatm)%charge
        yzinttot=yzinttot+a(iatm)%y*a(iatm)%z*a(iatm)%charge
        xzinttot=xzinttot+a(iatm)%x*a(iatm)%z*a(iatm)%charge
	    xxxinttot=xxxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%charge
	    yyyinttot=yyyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%charge
	    zzzinttot=zzzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%charge
	    yzzinttot=yzzinttot+a(iatm)%y*a(iatm)%z*a(iatm)%z*a(iatm)%charge
	    xzzinttot=xzzinttot+a(iatm)%x*a(iatm)%z*a(iatm)%z*a(iatm)%charge
	    xxzinttot=xxzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%z*a(iatm)%charge
	    yyzinttot=yyzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%z*a(iatm)%charge
	    xxyinttot=xxyinttot+a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%charge
	    xyyinttot=xyyinttot+a(iatm)%x*a(iatm)%y*a(iatm)%y*a(iatm)%charge
	    xyzinttot=xyzinttot+a(iatm)%x*a(iatm)%y*a(iatm)%z*a(iatm)%charge
        xxxxinttot=xxxxinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%charge
        yyyyinttot=yyyyinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%charge
        zzzzinttot=zzzzinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%charge
        xxxyinttot=xxxyinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%charge
        xxxzinttot=xxxzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%x*a(iatm)%z*a(iatm)%charge
        yyyxinttot=yyyxinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%x*a(iatm)%charge
        yyyzinttot=yyyzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%y*a(iatm)%z*a(iatm)%charge
        zzzxinttot=zzzxinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%x*a(iatm)%charge
        zzzyinttot=zzzyinttot+a(iatm)%z*a(iatm)%z*a(iatm)%z*a(iatm)%y*a(iatm)%charge
        xxyyinttot=xxyyinttot+a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%y*a(iatm)%charge
        xxzzinttot=xxzzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%z*a(iatm)%z*a(iatm)%charge
        yyzzinttot=yyzzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%z*a(iatm)%z*a(iatm)%charge
        xxyzinttot=xxyzinttot+a(iatm)%x*a(iatm)%x*a(iatm)%y*a(iatm)%z*a(iatm)%charge
        yyxzinttot=yyxzinttot+a(iatm)%y*a(iatm)%y*a(iatm)%x*a(iatm)%z*a(iatm)%charge
        zzxyinttot=zzxyinttot+a(iatm)%z*a(iatm)%z*a(iatm)%x*a(iatm)%y*a(iatm)%charge
    end do

    write(*,"(/,' X, Y, Z of center of positive charges (nuclear charges) in Angstrom',/,3f12.6)") &
    xnucdip/sum(a%charge)*b2a,ynucdip/sum(a%charge)*b2a,znucdip/sum(a%charge)*b2a
    write(*,"(' X, Y, Z of center of negative charges (electronic charges) in Angstrom',/,3f12.6)") &
    -xinttot/nelec*b2a,-yinttot/nelec*b2a,-zinttot/nelec*b2a

    write(*,"(/,' Dipole moment from nuclear charges (a.u.): ',3f11.6)") xnucdip,ynucdip,znucdip
    write(*,"(' Dipole moment from electrons (a.u.):       ',3f11.6)") xinttot,yinttot,zinttot
    xinttot=xinttot+xnucdip
    yinttot=yinttot+ynucdip
    zinttot=zinttot+znucdip
    write(*,*)
    write(*,"(' Dipole moment (a.u.): ',3f14.6)") xinttot,yinttot,zinttot
    write(*,"(' Dipole moment (Debye):',3f14.6)") xinttot*au2debye,yinttot*au2debye,zinttot*au2debye
    dipmag=sqrt(xinttot**2+yinttot**2+zinttot**2)
    write(*,"(' Magnitude of dipole moment:',f14.6,' a.u.',f14.6,' Debye')") dipmag,dipmag*au2debye
end if

!Print multipole moments
rrinttot=xxinttot+yyinttot+zzinttot
rrxinttot=xxxinttot+xyyinttot+xzzinttot
rryinttot=xxyinttot+yyyinttot+yzzinttot
rrzinttot=xxzinttot+yyzinttot+zzzinttot
write(*,*)
write(*,*) "Note: All units given below are in a.u."
write(*,"(/,' Quadrupole moments (Standard Cartesian form):')")
fac=1
!fac=au2debye*b2a !If using this factor, result will be identical to "Quadrupole moment (field-independent basis, Debye-Ang):" printed by Gaussian
write(*,"(' XX=',f12.6,'  XY=',f12.6,'  XZ=',f12.6)") xxinttot*fac,xyinttot*fac,xzinttot*fac
write(*,"(' YX=',f12.6,'  YY=',f12.6,'  YZ=',f12.6)") xyinttot*fac,yyinttot*fac,yzinttot*fac
write(*,"(' ZX=',f12.6,'  ZY=',f12.6,'  ZZ=',f12.6)") xzinttot*fac,yzinttot*fac,zzinttot*fac
write(*,"(' Quadrupole moments (Traceless Cartesian form):')")
!If removing the comment, the data will be identical to "Traceless Quadrupole moment (field-independent basis, Debye-Ang)" printed by Gaussian
QXX=(3*xxinttot-rrinttot)/2 !*au2debye*b2a/1.5D0
QYY=(3*yyinttot-rrinttot)/2 !*au2debye*b2a/1.5D0
QZZ=(3*zzinttot-rrinttot)/2 !*au2debye*b2a/1.5D0
QXY=3*xyinttot/2            !*au2debye*b2a/1.5D0
QXZ=3*xzinttot/2            !*au2debye*b2a/1.5D0
QYZ=3*yzinttot/2            !*au2debye*b2a/1.5D0
write(*,"(' XX=',f12.6,'  XY=',f12.6,'  XZ=',f12.6)") QXX,QXY,QXZ
write(*,"(' YX=',f12.6,'  YY=',f12.6,'  YZ=',f12.6)") QXY,QYY,QYZ
write(*,"(' ZX=',f12.6,'  ZY=',f12.6,'  ZZ=',f12.6)") QXZ,QYZ,QZZ
write(*,"(' Magnitude of the traceless quadrupole moment tensor:',f12.6)") sqrt(2D0/3D0*(QXX**2+QYY**2+QZZ**2))
R20=(3*zzinttot-rrinttot)/2D0 !Notice that the negative sign, because electrons carry negative charge
R2n1=dsqrt(3D0)*yzinttot
R2p1=dsqrt(3D0)*xzinttot
R2n2=dsqrt(3D0)*xyinttot
R2p2=dsqrt(3D0)/2D0*(xxinttot-yyinttot)
write(*,"(' Quadrupole moments (Spherical harmonic form):')")
write(*,"(' Q_2,0 =',f11.6,'   Q_2,-1=',f11.6,'   Q_2,1=',f11.6)") R20,R2n1,R2p1
write(*,"(' Q_2,-2=',f11.6,'   Q_2,2 =',f11.6)") R2n2,R2p2
write(*,"( ' Magnitude: |Q_2|=',f12.6)") dsqrt(R20**2+R2n1**2+R2p1**2+R2n2**2+R2p2**2)

R30=(5*zzzinttot-3*rrzinttot)/2D0
R3n1=dsqrt(3D0/8D0)*(5*yzzinttot-rryinttot)
R3p1=dsqrt(3D0/8D0)*(5*xzzinttot-rrxinttot)
R3n2=dsqrt(15D0)*xyzinttot
R3p2=dsqrt(15D0)*(xxzinttot-yyzinttot)/2D0
R3n3=dsqrt(5D0/8D0)*(3*xxyinttot-yyyinttot)
R3p3=dsqrt(5D0/8D0)*(xxxinttot-3*xyyinttot)
write(*,"(/,' Octopole moments (Cartesian form):')")
fac=1
!fac=au2debye*b2a*b2a !If using this factor, result will be identical to "Octapole moment (field-independent basis, Debye-Ang**2):" printed by Gaussian
write(*,"(' XXX=',f10.4,'  YYY=',f10.4,'  ZZZ=',f10.4,'  XYY=',f10.4,'  XXY=',f10.4)") &
xxxinttot*fac,yyyinttot*fac,zzzinttot*fac,xyyinttot*fac,xxyinttot*fac
write(*,"(' XXZ=',f10.4,'  XZZ=',f10.4,'  YZZ=',f10.4,'  YYZ=',f10.4,'  XYZ=',f10.4)") &
xxzinttot*fac,xzzinttot*fac,yzzinttot*fac,yyzinttot*fac,xyzinttot*fac
write(*,"(' Octopole moments (Spherical harmonic form):')")
write(*,"(' Q_3,0 =',f11.4,'  Q_3,-1=',f11.4,'  Q_3,1 =',f11.4)") R30,R3n1,R3p1
write(*,"(' Q_3,-2=',f11.4,'  Q_3,2 =',f11.4,'  Q_3,-3=',f11.4,'  Q_3,3 =',f11.4)") R3n2,R3p2,R3n3,R3p3
write(*,"( ' Magnitude: |Q_3|=',f12.4)") dsqrt(R30**2+R3n1**2+R3p1**2+R3n2**2+R3p2**2+R3n3**2+R3p3**2)

!The outputting order is identical to Gaussian
fac=1
!fac=au2debye*b2a*b2a*b2a !If using this, result will be identical to "Hexadecapole moment (field-independent basis, Debye-Ang**3):" printed by Gaussian
write(*,"(/,' Hexadecapole moments:')")
write(*,"(' XXXX=',f16.4,'  YYYY=',f16.4,'  ZZZZ=',f16.4)") xxxxinttot*fac,yyyyinttot*fac,zzzzinttot*fac
write(*,"(' XXXY=',f16.4,'  XXXZ=',f16.4,'  YYYX=',f16.4)") xxxyinttot*fac,xxxzinttot*fac,yyyxinttot*fac
write(*,"(' YYYZ=',f16.4,'  ZZZX=',f16.4,'  ZZZY=',f16.4)") yyyzinttot*fac,zzzxinttot*fac,zzzyinttot*fac
write(*,"(' XXYY=',f16.4,'  XXZZ=',f16.4,'  YYZZ=',f16.4)") xxyyinttot*fac,xxzzinttot*fac,yyzzinttot*fac
write(*,"(' XXYZ=',f16.4,'  YYXZ=',f16.4,'  ZZXY=',f16.4)") xxyzinttot*fac,yyxzinttot*fac,zzxyinttot*fac

if (allocated(b)) then
    ESE=ESEx+ESEy+ESEz
    write(*,"(/,a,f16.6)") " Electronic spatial extent <r^2>:",ESE
    write(*,"(' Components of <r^2>:  X=',f15.6,'  Y=',f15.6,'  Z=',f15.6)") ESEx,ESEy,ESEz
end if
end subroutine




!!------------ A general routine for obtaining energies of present orbitals in memory based on loaded Fock matrix    
subroutine calc_orb_energy
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 orbene(nmo),Emat(nbasis,nbasis),tmpmat(nbasis,nbasis)

if (.not.allocated(CObasa)) then
    write(*,"(a)") " Error: To use this function, the input file must contain basis function information! See Section 2.5 of Multiwfn manual for detail"
    write(*,*) "Press ENTER button to return"
    return
end if

call loadFockfile(istat)
if (istat==1) then
	write(*,*) "Error: Unable to evaluate orbital energies!"
else
    write(*,*) "Evaluating orbital energies..."
	!Emat=matmul(matmul(transpose(CObasa),FmatA),CObasa)
    tmpmat=matmul_blas(CObasa,FmatA,nbasis,nbasis,1,0)
    Emat=matmul_blas(tmpmat,CObasa,nbasis,nbasis)
	do iorb=1,nbasis
		orbene(iorb)=Emat(iorb,iorb)
	end do
    if (allocated(CObasb)) then
	    !Emat=matmul(matmul(transpose(CObasb),FmatB),CObasb)
        tmpmat=matmul_blas(CObasb,FmatB,nbasis,nbasis,1,0)
        Emat=matmul_blas(tmpmat,CObasb,nbasis,nbasis)
	    do iorb=1,nbasis
		    orbene(nbasis+iorb)=Emat(iorb,iorb)
	    end do
    end if
	write(*,*) "Orbital energies have been successfully evaluated!"
    write(*,*)
    write(*,*) "0 Do nothing"
    write(*,*) "1 Export orbital energies to orbene.txt in current folder"
    write(*,*) "2 Replace the original orbital energies in memory by the newly evaluated ones"
    write(*,*) "3 Do both 1 and 2"
    read(*,*) isel
    if (isel==1.or.isel==3) then
        open(10,file="orbene.txt")
        if (wfntype==0.or.wfntype==2.or.wfntype==3) then
            do iorb=1,nbasis
                write(10,"(i7,'  Occ=',f8.4,'  E=',f16.8,' Hartree',f12.4,' eV')") iorb,MOocc(iorb),orbene(iorb),orbene(iorb)*au2eV
            end do
        else
            write(10,*) "===== Alpha orbitals ====="
            do iorb=1,nbasis
                write(10,"(i7,'  Occ=',f8.4,'  E=',f16.8,' Hartree',f12.4,' eV')") iorb,MOocc(iorb),orbene(iorb),orbene(iorb)*au2eV
            end do
            write(10,*)
            write(10,*) "===== Beta orbitals ====="
            do iorb=1,nbasis
                write(10,"(i7,'  Occ=',f8.4,'  E=',f16.8,' Hartree',f12.4,' eV')") iorb,MOocc(nbasis+iorb),orbene(nbasis+iorb),orbene(nbasis+iorb)*au2eV
            end do
        end if
        close(10)
        write(*,*) "The new orbital energies have been exported to orbene.txt in current folder!"
    end if
    if (isel==2.or.isel==3) then
        MOene=orbene
        write(*,*) "The original orbital energies in memory have been replaced by the new ones!"
    end if
end if
end subroutine





!!------- Generate randomly displaced geometries
subroutine displace_geom
use defvar
use util
implicit real*8 (a-h,o-z)
character c2000tmp*2000,c80tmp*80
integer atmlist(ncenter)

write(*,*)
write(*,*) " ------------------ Generate randomly displaced geometries ------------------"
write(*,*) "Input index of the atoms that you want to randomly displace, e.g. 2,3,7-10"
write(*,*) "To choose the whole system, press ENTER button directly"
write(*,*) "To exit, input ""q"""
read(*,"(a)") c2000tmp
if (c2000tmp==" ".or.index(c2000tmp,'a')/=0) then
    nsel=ncenter
    forall(i=1:nsel) atmlist(i)=i
else if (c2000tmp=='q') then
    return
else
    call str2arr(c2000tmp,nsel,atmlist)
end if

write(*,*) "Displace which Cartesian coordinates?"
write(*,*) "1 X coordinate"
write(*,*) "2 Y coordinate"
write(*,*) "3 Z coordinate"
write(*,*) "4 X and Y coordinates"
write(*,*) "5 Y and Z coordinates"
write(*,*) "6 X and Z coordinates"
write(*,*) "7 X, Y and Z coordinates"
read(*,*) itype

write(*,*) "Input standard variation of displacement in Angstrom, e.g. 0.01"
write(*,*) "If you press ENTER button directly, 0.03 Angstrom will be used"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
    stdvar=0.03D0
else
    read(c80tmp,*) stdvar
end if
stdvar=stdvar/b2a

write(*,*) "Generate how many geometries? e.g. 4"
write(*,*) "If you press ENTER button directly, only one geometry will be generated"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
    numgen=1
else
    read(c80tmp,*) numgen
end if

open(10,file="new.xyz",status="replace")
do igen=1,numgen
    do i=1,nsel
        iatm=atmlist(i)
        !Use basic form of Box–Muller transform to generate random numbers &
        !satisfying normal distribution (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
        CALL RANDOM_NUMBER(ran1)
        CALL RANDOM_NUMBER(ran2)
        CALL RANDOM_NUMBER(ran3)
        tmp=stdvar*dsqrt(-2*log(ran1))
        pertx=0
        perty=0
        pertz=0
        if (itype==1.or.itype==4.or.itype==6.or.itype==7) pertx=tmp*cos(2*pi*ran2)
        if (itype==2.or.itype==4.or.itype==5.or.itype==7) perty=tmp*sin(2*pi*ran2)
        if (itype==3.or.itype==5.or.itype==6.or.itype==7) pertz=tmp*cos(2*pi*ran3)
        a(iatm)%x=a(iatm)%x+pertx
        a(iatm)%y=a(iatm)%y+perty
        a(iatm)%z=a(iatm)%z+pertz
        if (numgen==1) write(*,"(' Moved',i6,'(',a,') by X:',f9.5,'  Y:',f9.5,'  Z:',f9.5,'  Total:',f9.5,' Ang')") &
        iatm,a(iatm)%name,pertx*b2a,perty*b2a,pertz*b2a,dsqrt(pertx**2+perty**2+pertz**2)*b2a
    end do
    write(10,"(i6)") ncenter
    write(10,"(' Geometry',i10,' Generated by Multiwfn')") igen
    do i=1,ncenter
	    write(10,"(a,3f16.8)") a(i)%name,a(i)%x*b2a,a(i)%y*b2a,a(i)%z*b2a
    end do
    a%x=a_org%x
    a%y=a_org%y
    a%z=a_org%z
end do
close(10)

write(*,*)
if (numgen==1) then
    write(*,*) "Done! The displaced geometry has been updated to new.xyz in current folder"
else
    write(*,*) "Done! The displaced geometries have been updated to new.xyz in current folder"
end if
end subroutine




!!---------- Geometry relevant operations on the present system
subroutine geom_operation
use defvar
use deftype
use GUI
use util
implicit real*8 (a-h,o-z)
integer fragsel(ncenter),nfragsel,iffrag(ncenter)
character c200tmp*200,c2000tmp*2000,selectyn
real*8 mat(3,3),matt1(3,3),matt2(3,3)
real*8 vec(3),vec2(3),vecr(1,3),vecc(3,1) !vector, row vector, column vector
type(atomtype),allocatable :: a_old(:)
type(atomtype) atmp
real*8 inertia(3,3),eigvecmat(3,3),eigvalarr(3)
real*8 rcoord(3),fcoord(3)
integer,allocatable :: tmparr(:),intarr(:)
real*8,allocatable :: fcoordall(:,:)

allocate(intarr(ncenter)) !It should be reallocated when constructing supercell

if (allocated(b)) then
    write(*,"(/,a)") " NOTE: This function only changes geometry, the wavefunction will not be correspondingly affected!"
end if
do while(.true.)
    write(*,*)
    write(*,*) "      ----------------- Geometry relevant operations -----------------"
    write(*,*) "-10 Return"
    write(*,*) "-9 Restore the original geometry"
    write(*,*) "-3 Output system to .gjf file   -4 Output system to .cif file"
    write(*,*) "-1 Output system to .xyz file   -2 Output system to .pdb file"
    write(*,*) " 0 Visualize current geometry"
    write(*,*) " 1 Translate selected atoms according to a translation vector"
    write(*,*) " 2 Translate the system such that the center of selected atoms is at origin"
    write(*,*) " 3 Rotate selected atoms around a Cartesian axis or a bond"
    write(*,*) " 4 Rotate selected atoms by applying a given rotation matrix"
    write(*,*) " 5 Make a bond parallel to a vector or Cartesian axis"
    write(*,*) " 6 Make a vector parallel to a vector or Cartesian axis"
    write(*,*) " 7 Make electric dipole moment parallel to a vector or Cartesian axis"
    write(*,*) " 8 Make longest axis of selected atoms parallel to a vector or Cartesian axis"
    write(*,*) " 9 Mirror inversion for selected atoms"
    write(*,*) "10 Center inversion for selected atoms"
    write(*,*) "11 Make the plane defined by selected atoms parallel to a Cartesian plane"
    if (ifPBC>0) then
        write(*,*) "12 Scale Cartesian coordinates of selected atoms and cell vectors"
    else
        write(*,*) "12 Scale Cartesian coordinates of selected atoms"
    end if
    write(*,*) "13 Reorder atom sequence"
    write(*,*) "15 Add an atom     16 Remove some atoms      17 Crop some atoms"
    write(*,*) "18 Generate randomly displaced geometries"
    write(*,*) "19 Translate and duplicate cell (constructing supercell)"
    write(*,*) "20 Make truncated molecules by cell boundary whole"
    write(*,*) "21 Scale cell length and atom coordinates correspondingly"
    write(*,*) "22 Wrap all atoms outside the cell into the cell"
    write(*,*) "23 Translate system along cell axes by given distances"
    write(*,*) "24 Translate system to center selected part in the cell"
    write(*,*) "25 Extract a molecular cluster (central molecule + neighbouring ones)"
    write(*,*) "26 Set cell information      27 Add boundary atoms     28 Axes interconversion"
    read(*,*) isel
    
    if (isel==-10) then
        return
    else if (isel==-9) then
        deallocate(a)
        if (allocated(connmat)) deallocate(connmat)
        ncenter=ncenter_org
        allocate(a(ncenter))
        a=a_org
        ifPBC=ifPBC_org
        cellv1=cellv1_org
        cellv2=cellv2_org
        cellv3=cellv3_org
        write(*,*) "Current geometry has been restored to the original one"
    else if (isel==-4) then
        if (ifPBC<3) then
            write(*,*) "Error: This function can only be used for three-dimension periodic systems"
            write(*,*) "Press ENTER button to continue"
            read(*,*)
            cycle
        else
            call outcif_wrapper
        end if
    else if (isel==-3) then
    	    call outgjf_wrapper
    else if (isel==-2) then
    	    call outpdb_wrapper
    else if (isel==-1) then
	    call outxyz_wrapper
    else if (isel==0) then
        if (ncenter<=300) then
            write(*,*) "Atom list:"
		    call showcoordA(0)
        else
			write(*,"(a)") " Note: There are more than 300 atoms, so their information is not shown here now. &
            &To print, in the manu bar please select ""Tools"" - ""Print XYZ coordinates"""
        end if
        call drawmolgui
        
    else if (isel==1) then !Translate selected part
        write(*,*) "Input indices of the atoms you want to translate, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        write(*,*) "Input X, Y, Z of translation vector in Angstrom, e.g. 2.3,1.05,-3.24"
        read(*,*) vec
        vec=vec/b2a
        do idx=1,nfragsel
            iatm=fragsel(idx)
            a(iatm)%x=a(iatm)%x+vec(1)
            a(iatm)%y=a(iatm)%y+vec(2)
            a(iatm)%z=a(iatm)%z+vec(3)
        end do
        write(*,*) "Done! The geometry as been updated"
        
    else if (isel==2) then !Move center of selected part to origin
        write(*,*) "Choose center type:"
        write(*,*) "1 Center of geometry"
        write(*,*) "2 Center of mass"
        write(*,*) "3 Center of nuclear charges"
        read(*,*) icentype
        
        write(*,*) "Input indices of the atoms for which the center will be defined, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        write(*,*) "Note: The whole system will be translated"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        
        if (icentype==1) then
            vec(1)=sum(a(fragsel(1:nfragsel))%x)/nfragsel
            vec(2)=sum(a(fragsel(1:nfragsel))%y)/nfragsel
            vec(3)=sum(a(fragsel(1:nfragsel))%z)/nfragsel
        else if (icentype==2) then
            totmass=sum(atmwei(a(fragsel(1:nfragsel))%index))
            vec(1)=sum(a(fragsel(1:nfragsel))%x*atmwei(a(fragsel(1:nfragsel))%index))/totmass
            vec(2)=sum(a(fragsel(1:nfragsel))%y*atmwei(a(fragsel(1:nfragsel))%index))/totmass
            vec(3)=sum(a(fragsel(1:nfragsel))%z*atmwei(a(fragsel(1:nfragsel))%index))/totmass
        else if (icentype==3) then
            totchg=sum(a(fragsel(1:nfragsel))%charge)
            vec(1)=sum(a(fragsel(1:nfragsel))%x*a(fragsel(1:nfragsel))%charge)/totchg
            vec(2)=sum(a(fragsel(1:nfragsel))%y*a(fragsel(1:nfragsel))%charge)/totchg
            vec(3)=sum(a(fragsel(1:nfragsel))%z*a(fragsel(1:nfragsel))%charge)/totchg
        end if
        write(*,"(' X, Y, Z coordinate of the center:',3f12.6,' Bohr')") vec
        
        do iatm=1,ncenter
            a(iatm)%x=a(iatm)%x-vec(1)
            a(iatm)%y=a(iatm)%y-vec(2)
            a(iatm)%z=a(iatm)%z-vec(3)
        end do
        write(*,*) "Done! The geometry as been updated"
            
    else if (isel==3) then !Rotate selected part
        write(*,*) "Input indices of the atoms you want to rotate, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        write(*,*) "How to rotate?"
        write(*,*) "1 Rotate around X axis"
        write(*,*) "2 Rotate around Y axis"
        write(*,*) "3 Rotate around Z axis"
        write(*,*) "4 Rotate around a bond"
        write(*,*) "5 Rotate around a specific vector"
        read(*,*) irotsel
        if (irotsel>=1.and.irotsel<=5) then
            write(*,*) "Rotate how many degrees? e.g. 132  (may be negative)"
            read(*,*) rotdeg
        end if
        
        if (irotsel==1) then
            call get_rotmat_aroundXYZ(0D0,0D0,rotdeg,mat)
        else if (irotsel==2) then
            call get_rotmat_aroundXYZ(0D0,rotdeg,0D0,mat)
        else if (irotsel==3) then
            call get_rotmat_aroundXYZ(rotdeg,0D0,0D0,mat)
        else if (irotsel==4) then
            write(*,*) "Input indices of two atoms to define the bond, e.g. 4,6"
            read(*,*) iatm,jatm
            vec(1)=a(jatm)%x-a(iatm)%x
            vec(2)=a(jatm)%y-a(iatm)%y
            vec(3)=a(jatm)%z-a(iatm)%z
            call get_rotmat_aroundvec(vec,rotdeg,mat)
        else if (irotsel==5) then
            write(*,*) "Input X,Y,Z component of the vector, e.g. -1.05,0.26,8.3"
            read(*,*) vec(:)
            call get_rotmat_aroundvec(vec,rotdeg,mat)
        end if
        write(*,*) "Rotation matrix applied to coordinates of selected atoms:"
        write(*,"(3f16.8)") mat(1,:)
        write(*,"(3f16.8)") mat(2,:)
        write(*,"(3f16.8)") mat(3,:)
        write(*,*)
        !Apply rotation matrix
        do idx=1,nfragsel
            iatm=fragsel(idx)
            vecc(1,1)=a(iatm)%x;vecc(2,1)=a(iatm)%y;vecc(3,1)=a(iatm)%z
            vecc=matmul(mat,vecc)
            a(iatm)%x=vecc(1,1);a(iatm)%y=vecc(2,1);a(iatm)%z=vecc(3,1)
        end do
        write(*,*) "Done! The geometry as been updated"
        
    else if (isel==4) then !Rotate the selected part by applying a given rotation matrix
        write(*,*) "Input indices of the atoms you want to apply rotation matrix, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        write(*,*) "Input XX,XY,XZ of rotation matrix"
        read(*,*) mat(1,:)
        write(*,*) "Input YX,YY,YZ of rotation matrix"
        read(*,*) mat(2,:)
        write(*,*) "Input ZX,ZY,ZZ of rotation matrix"
        read(*,*) mat(3,:)
        write(*,*)
        write(*,*) "How to apply the rotation matrix (mat)?"
        write(*,*) "1 v_new = mat * v_old  (v is column vector of coordinates)"
        write(*,*) "2 v_new = v_old * mat  (v is row vector of coordinates)"
        write(*,"(a)") " Note: Option 1 corresponds to common situation. However, if your purpose is to apply the &
        &eigenvector matrix of Hessian printed by Multiwfn as rotation matrix to eliminate its non-diagonal terms, you should choose 2"
        read(*,*) iapply
        
        do idx=1,nfragsel
            iatm=fragsel(idx)
            if (iapply==1) then
                vecc(1,1)=a(iatm)%x;vecc(2,1)=a(iatm)%y;vecc(3,1)=a(iatm)%z
                vecc=matmul(mat,vecc)
                a(iatm)%x=vecc(1,1);a(iatm)%y=vecc(2,1);a(iatm)%z=vecc(3,1)
            else
                vecr(1,1)=a(iatm)%x;vecr(1,2)=a(iatm)%y;vecr(1,3)=a(iatm)%z
                vecr=matmul(vecr,mat)
                a(iatm)%x=vecr(1,1);a(iatm)%y=vecr(1,2);a(iatm)%z=vecr(1,3)
            end if
        end do
        write(*,*) "Done! The geometry as been updated"
        
    else if (isel==5.or.isel==6.or.isel==7.or.isel==8) then !Reorientate the system to make a bond / vector / dipole moment / longest axis parallel to a vector
        if (isel==5) then
            write(*,*) "Input indices of two atoms to define the bond, e.g. 4,6"
            read(*,*) iatm,jatm
            vec(1)=a(jatm)%x-a(iatm)%x
            vec(2)=a(jatm)%y-a(iatm)%y
            vec(3)=a(jatm)%z-a(iatm)%z
            write(*,"(' X,Y,Z of the bond vector:',3f12.6)") vec
        else if (isel==6) then
            write(*,*) "Input X,Y,Z component of the vector, e.g. -1.05,0.26,8.3"
            read(*,*) vec(:)
        else if (isel==7) then
            if (.not.allocated(b)) then
                write(*,"(a)") " Error: Wavefunction information is not provided by your input file, this function does not work!"
                write(*,*) "Press ENTER button to continue"
                read(*,*)
                cycle
            end if
            call get_dipole_moment(vec)
            write(*,"(' X,Y,Z of the dipole moment:',3f12.6,' a.u.')") vec
        else if (isel==8) then
            write(*,*) "Input indices of the atoms for which the longest axis will be determined"
            write(*,*) "For example, 2,3,7-10,15"
            write(*,*) "If press ENTER button directly, the whole system will be selected"
            read(*,"(a)") c2000tmp
            if (c2000tmp==" ") then
                forall(iatm=1:ncenter) fragsel(iatm)=iatm
                nfragsel=ncenter
            else
                call str2arr(c2000tmp,nfragsel,fragsel)
            end if
            totmass=sum(atmwei(a(fragsel(1:nfragsel))%index))
            cenmassx=sum( a(fragsel(1:nfragsel))%x*atmwei(a(fragsel(1:nfragsel))%index) )/totmass
            cenmassy=sum( a(fragsel(1:nfragsel))%y*atmwei(a(fragsel(1:nfragsel))%index) )/totmass
            cenmassz=sum( a(fragsel(1:nfragsel))%z*atmwei(a(fragsel(1:nfragsel))%index) )/totmass
            inertia(1,1)=sum(atmwei(a(fragsel(1:nfragsel))%index)*( (a(fragsel(1:nfragsel))%y-cenmassy)**2+(a(fragsel(1:nfragsel))%z-cenmassz)**2) )*b2a*b2a
            inertia(2,2)=sum(atmwei(a(fragsel(1:nfragsel))%index)*( (a(fragsel(1:nfragsel))%x-cenmassx)**2+(a(fragsel(1:nfragsel))%z-cenmassz)**2) )*b2a*b2a
            inertia(3,3)=sum(atmwei(a(fragsel(1:nfragsel))%index)*( (a(fragsel(1:nfragsel))%x-cenmassx)**2+(a(fragsel(1:nfragsel))%y-cenmassy)**2) )*b2a*b2a
            inertia(1,2)=-sum(atmwei(a(fragsel(1:nfragsel))%index)*(a(fragsel(1:nfragsel))%x-cenmassx)*(a(fragsel(1:nfragsel))%y-cenmassy))*b2a*b2a
            inertia(2,1)=inertia(1,2)
            inertia(1,3)=-sum(atmwei(a(fragsel(1:nfragsel))%index)*(a(fragsel(1:nfragsel))%x-cenmassx)*(a(fragsel(1:nfragsel))%z-cenmassz))*b2a*b2a
            inertia(3,1)=inertia(1,3)
            inertia(2,3)=-sum(atmwei(a(fragsel(1:nfragsel))%index)*(a(fragsel(1:nfragsel))%y-cenmassy)*(a(fragsel(1:nfragsel))%z-cenmassz))*b2a*b2a
            inertia(3,2)=inertia(2,3)
            call diagmat(inertia,eigvecmat,eigvalarr,300,1D-12)
            tmpval=eigvalarr(1)
            vec=eigvecmat(:,1)
            do i=2,3
                if (eigvalarr(i)<tmpval) then
                    tmpval=eigvalarr(i)
                    vec=eigvecmat(:,i)
                end if
            end do
            write(*,"(' X,Y,Z of the longest axis of selected atoms:',3f12.6)") vec
        end if
        write(*,*)
        write(*,*) "1: Parallel to X axis"
        write(*,*) "2: Parallel to Y axis"
        write(*,*) "3: Parallel to Z axis"
        write(*,*) "4: Parallel to a given vector"
        read(*,*) iaxis
        vec2=0
        if (iaxis==1) then
            vec2(1)=1
        else if (iaxis==2) then
            vec2(2)=1
        else if (iaxis==3) then
            vec2(3)=1
        else if (iaxis==4) then
            write(*,*) "Input the vector, e.g. 2.5,0,-0.2"
            read(*,*) vec2
        end if
        
        !Generate rotation matrix
        vec=vec/dsqrt(sum(vec**2)) !Normalize
        vec2=vec2/dsqrt(sum(vec2**2)) !Normalize
        call rotmat_vec1_vec2(vec,vec2,mat)
        write(*,*) "Rotation matrix applied to coordinates of selected atoms:"
        write(*,"(3f16.8)") mat(1,:)
        write(*,"(3f16.8)") mat(2,:)
        write(*,"(3f16.8)") mat(3,:)
        write(*,*)
        !Apply rotation matrix
        do iatm=1,ncenter
            vecc(1,1)=a(iatm)%x;vecc(2,1)=a(iatm)%y;vecc(3,1)=a(iatm)%z
            vecc=matmul(mat,vecc)
            a(iatm)%x=vecc(1,1);a(iatm)%y=vecc(2,1);a(iatm)%z=vecc(3,1)
        end do
        write(*,*) "Done! The geometry as been updated"
        
    else if (isel==9) then !Mirror invertion for selected atoms
        write(*,*) "Input indices of the atoms you want to perform mirror inversion, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        write(*,*) "How to invert?"
        write(*,*) "1 Invert with respect to XY plane"
        write(*,*) "2 Invert with respect to YZ plane"
        write(*,*) "3 Invert with respect to XZ plane"
        read(*,*) iinvsel
        !Apply inversion
        do idx=1,nfragsel
            iatm=fragsel(idx)
            if (iinvsel==1) then
                a(iatm)%z=-a(iatm)%z
            else if (iinvsel==2) then
                a(iatm)%x=-a(iatm)%x
            else if (iinvsel==3) then
                a(iatm)%y=-a(iatm)%y
            end if
        end do
        write(*,*) "Done! The geometry as been updated"
        
    else if (isel==10) then !Center invertion for selected atoms
        write(*,*) "Input indices of the atoms you want to perform center inversion, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        do idx=1,nfragsel
            iatm=fragsel(idx)
            a(iatm)%x=-a(iatm)%x
            a(iatm)%y=-a(iatm)%y
            a(iatm)%z=-a(iatm)%z
        end do
        write(*,*) "Done! The geometry as been updated"
        
    else if (isel==11) then !Make the plane defined by selected atoms parallel to a Cartesian plane
        write(*,*) "Input indices of the atoms you want to fit a plane, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        call ptsfitplane(fragsel(1:nfragsel),nfragsel,xnor,ynor,znor,rnouse,rmsfit)
        write(*,"(' RMS error of the plane fitting:',f12.6,' Angstrom')") rmsfit*b2a
        facnorm=sqrt(xnor**2+ynor**2+znor**2)
        vec(1)=xnor/facnorm !Normalize normal vector
        vec(2)=ynor/facnorm
        vec(3)=znor/facnorm
        write(*,"(' X, Y and Z of the unit normal vector is',3f13.8)") xnor,ynor,znor
        write(*,*)
        write(*,*) "Choose the plane to which the fitted plane will be parallel to"
        write(*,*) "1: XY"
        write(*,*) "2: YZ"
        write(*,*) "3: XZ"
        read(*,*) iple
        vec2=0
        if (iple==1) then
            vec2(3)=1
        else if (iple==2) then
            vec2(1)=1
        else if (iple==3) then
            vec2(2)=1
        end if
        
        !Generate rotation matrix
        vec=vec/dsqrt(sum(vec**2)) !Normalize
        call rotmat_vec1_vec2(vec,vec2,mat)
        write(*,*) "Rotation matrix applied to coordinates of selected atoms:"
        write(*,"(3f16.8)") mat(1,:)
        write(*,"(3f16.8)") mat(2,:)
        write(*,"(3f16.8)") mat(3,:)
        write(*,*)
        !Apply rotation matrix
        do iatm=1,ncenter
            vecc(1,1)=a(iatm)%x;vecc(2,1)=a(iatm)%y;vecc(3,1)=a(iatm)%z
            vecc=matmul(mat,vecc)
            a(iatm)%x=vecc(1,1);a(iatm)%y=vecc(2,1);a(iatm)%z=vecc(3,1)
        end do
        write(*,*) "Done! The geometry as been updated"
    
    else if (isel==12) then !Scale Cartesian coordinate of selected atoms
        write(*,*) "Input indices of the atoms you want to scale coordinates, e.g. 2,3,7-10"
        write(*,*) "If press ENTER button directly, the whole system will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            forall(iatm=1:ncenter) fragsel(iatm)=iatm
            nfragsel=ncenter
        else
            call str2arr(c2000tmp,nfragsel,fragsel)
        end if
        write(*,*) "Scale which coordinates? Input the corresponding index, e.g. 2"
        write(*,*) "0: All"
        write(*,*) "1: X"
        write(*,*) "2: Y"
        write(*,*) "3: Z"
        read(*,*) icoord
        write(*,*) "Input scale factor, e.g. 0.95"
        read(*,*) scltmp
        do idx=1,nfragsel
            iatm=fragsel(idx)
            if (icoord==0.or.icoord==1) a(iatm)%x=scltmp*a(iatm)%x
            if (icoord==0.or.icoord==2) a(iatm)%y=scltmp*a(iatm)%y
            if (icoord==0.or.icoord==3) a(iatm)%z=scltmp*a(iatm)%z
        end do
        write(*,*) "Done! The geometry as been updated"
        if (ifPBC>0) then
            write(*,"(/,a)") " Also scale corresponding Cartesian component of all cell vectors? (y/n)"
            write(*,"(a)") " (For example, if you have chosen to scale X coordinates before, then X component of all cell vectors will also be scaled)"
            read(*,*) selectyn
            if (selectyn=='y'.or.selectyn=='Y') then
                if (icoord==0.or.icoord==1) then
                    cellv1(1)=cellv1(1)*scltmp
                    cellv2(1)=cellv2(1)*scltmp
                    cellv3(1)=cellv3(1)*scltmp
                end if
                if (icoord==0.or.icoord==2) then
                    cellv1(2)=cellv1(2)*scltmp
                    cellv2(2)=cellv2(2)*scltmp
                    cellv3(2)=cellv3(2)*scltmp
                end if
                if (icoord==0.or.icoord==3) then
                    cellv1(3)=cellv1(3)*scltmp
                    cellv2(3)=cellv2(3)*scltmp
                    cellv3(3)=cellv3(3)*scltmp
                end if
                write(*,*) "Cell vectors have been scaled"
            end if
        end if
        
    else if (isel==13) then !Reorder atom sequence
        if (allocated(b)) write(*,"(a,/)") " Notice: Wavefunction information will not be correspondingly updated if you use this function to reorder atoms"
        write(*,*) "How to reorder atom sequence?"
        write(*,*) "-1 Reverse order of all atoms"
        write(*,*) "0 Return"
        write(*,*) "1 Reorder according to X coordinate (small to large)"
        write(*,*) "2 Reorder according to Y coordinate (small to large)"
        write(*,*) "3 Reorder according to Z coordinate (small to large)"
        write(*,*) "4 Reorder according to bonding (making indices contiguous in every fragment)"
        write(*,*) "5 Reorder according to element index (large to small)"
        write(*,*) "6 Put hydrogens at the end"
        write(*,*) "7 Exchange two atoms"
        write(*,*) "8 Reorder atoms according to a given list"
        write(*,*) "9 Put specific atoms prior to others"
        read(*,*) imode
        
        if (imode==-1) then
            do iatm=1,int(ncenter/2)
                jatm=ncenter-iatm+1
                atmp=a(iatm)
                a(iatm)=a(jatm)
                a(jatm)=atmp
            end do
        else if (imode==1.or.imode==2.or.imode==3) then
            allocate(a_old(ncenter))
            a_old=a
            forall(i=1:ncenter) intarr(i)=i
            if (imode==1) call sortr8(a(:)%x,list=intarr)
            if (imode==2) call sortr8(a(:)%y,list=intarr)
            if (imode==3) call sortr8(a(:)%z,list=intarr)
            do iatm=1,ncenter
                a(iatm)=a_old(intarr(iatm))
            end do
            deallocate(a_old)
        else if (imode==4) then
            iatm=1
            nassign=0
            intarr=1 !intarr(i)=1 means atom i has not been assigned into fragment
            allocate(a_old(ncenter))
            a_old=a
            do while(.true.)
                call getfragatoms(iatm,iffrag)
                do jatm=1,ncenter
                    if (iffrag(jatm)==1) then !jatm is in the newly chosen fragment
                        nassign=nassign+1
                        a(nassign)=a_old(jatm)
                        intarr(jatm)=0
                    end if
                end do
                if (nassign==ncenter) exit
                do iatm=1,ncenter
                    if (intarr(iatm)==1) exit
                end do
            end do
            deallocate(a_old)
            if (allocated(connmat)) deallocate(connmat) !getfragatoms generated connmat, while after changing order, the relationship is no longer valid
        else if (imode==5) then
            do iatm=1,ncenter
                do jatm=iatm+1,ncenter
                    if (a(iatm)%index<a(jatm)%index) then
                        atmp=a(iatm)
                        a(iatm)=a(jatm)
                        a(jatm)=atmp
                    end if
                end do
            end do
        else if (imode==6) then
            allocate(a_old(ncenter))
            a_old=a
            itmp=0
            do iatm=1,ncenter
                if (a_old(iatm)%index/=1) then
                    itmp=itmp+1
                    a(itmp)=a_old(iatm)
                end if
            end do
            do iatm=1,ncenter
                if (a_old(iatm)%index==1) then
                    itmp=itmp+1
                    a(itmp)=a_old(iatm)
                end if
            end do
            deallocate(a_old)
        else if (imode==7) then
            write(*,*) "Input indices of two atoms to exchange them, e.g. 5,8"
            write(*,*) "You can successively enter many times, input ""q"" can exit"
            do while(.true.)
                read(*,"(a)") c200tmp
                if (index(c200tmp,'q')/=0) then
                    exit
                else
                    read(c200tmp,*) iatm,jatm
                    atmp=a(iatm)
                    a(iatm)=a(jatm)
                    a(jatm)=atmp
                    write(*,"(i5,' and',i5,' have been exchanged')") iatm,jatm
                end if
            end do
        else if (imode==8) then
            write(*,*) "Input the path of the file containing new order, e.g. C:\new.txt"
            write(*,"(a)") " Note: The content of this file should contain indices of all atoms in new order, for example (a system containing 6 atoms):"
            write(*,*) "5"
            write(*,*) "2-4"
            write(*,*) "6"
            write(*,*) "1"
            do while(.true.)
                read(*,"(a)") c200tmp
	            inquire(file=c200tmp,exist=alive)
	            if (alive) exit
	            write(*,*) "Error: Cannot find the file, input again!"
            end do
            open(10,file=c200tmp,status="old")
            ncenter=0
            allocate(a_old(ncenter))
            a_old=a
            do while(.true.)
                read(10,"(a)",iostat=ierror) c2000tmp
                if (ierror/=0.or.c2000tmp==" ") exit
                call str2arr(c2000tmp,ntmp,intarr)
                do itmp=1,ntmp
                    iatm=intarr(itmp)
                    ncenter=ncenter+1
                    a(ncenter)=a_old(iatm)
                end do
            end do
            close(10)
            if (ncenter/=ncenter_org) then
                write(*,"(a)") " Error: The atom list does not define order of all atoms! The original atom information is restored"
                ncenter=ncenter_org
                a=a_old
                cycle
            end if
            deallocate(a_old)
        else if (imode==9) then
            write(*,*) "Input indices of the atoms, they will appear in front of others. e.g. 2,3,7-10"
            read(*,"(a)") c2000tmp
            call str2arr(c2000tmp,ntmp)
            allocate(tmparr(ntmp))
            call str2arr(c2000tmp,ntmp,tmparr)
            allocate(a_old(ncenter))
            a_old=a
            ncenter_old=ncenter
            ncenter=0
            do itmp=1,ntmp
                ncenter=ncenter+1
                a(ncenter)=a_old(tmparr(itmp))
            end do
            do iatm=1,ncenter_old
                if (all(tmparr/=iatm)) then
                    ncenter=ncenter+1
                    a(ncenter)=a_old(iatm)
                end if
            end do
            deallocate(a_old)
        end if
        write(*,*) "Done!"
        
    else if (isel==15) then !Add an atom
        write(*,*) "Input element of the newly added atom, e.g. Fe"
        read(*,*) c200tmp
        call elename2idx(c200tmp,iele)
        if (iele==0) then
            write(*,*) "Error: The element cannot be recognized. Press ENTER button to return"
            read(*,*)
            cycle
        end if
        write(*,*) "Input X,Y,Z coordinate in Angstrom, e.g. 3.2,0.9,-1.5"
        read(*,*) xtmp,ytmp,ztmp
        xtmp=xtmp/b2a
        ytmp=ytmp/b2a
        ztmp=ztmp/b2a
        if (allocated(a_tmp)) deallocate(a_tmp)
        allocate(a_tmp(ncenter))
        a_tmp=a
        deallocate(a)
        ncenter=ncenter+1
        allocate(a(ncenter))
        a(1:ncenter-1)=a_tmp
        a(ncenter)%index=iele
        a(ncenter)%charge=iele
        a(ncenter)%name=ind2name(iele)
        a(ncenter)%x=xtmp
        a(ncenter)%y=ytmp
        a(ncenter)%z=ztmp
        a(ncenter)%resid=1
        a(ncenter)%resname=" "
        deallocate(a_tmp)
        write(*,*) "Done! The new atom has been added"
        
    else if (isel==16) then !Delete some atoms
        write(*,*) "Input indices of the atoms you want to remove, e.g. 2,3-5,7-10"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp,intarr)
        if (allocated(a_tmp)) deallocate(a_tmp)
        allocate(a_tmp(ncenter))
        a_tmp=a
        deallocate(a)
        ncenter=ncenter-ntmp
        allocate(a(ncenter))
        ncenter=0
        do iatm=1,size(a_tmp)
            if (any(intarr(1:ntmp)==iatm)) cycle
            ncenter=ncenter+1
            a(ncenter)=a_tmp(iatm)
        end do
        deallocate(a_tmp)
        write(*,"(a,i5,a)") " Done!",ntmp," atoms have been removed"
        
    else if (isel==17) then !Crop some atoms
        write(*,*) "Input indices of the atoms you want to crop, e.g. 2,3,7-10"
        write(*,*) "Note: The atoms not belonging to the selected range will be removed"
        read(*,"(a)") c2000tmp
        call str2arr(c2000tmp,ntmp,intarr)
        if (allocated(a_tmp)) deallocate(a_tmp)
        allocate(a_tmp(ncenter))
        a_tmp=a
        deallocate(a)
        ncenter=ntmp
        allocate(a(ncenter))
        do idx=1,ntmp
            iatm=intarr(idx)
            a(idx)=a_tmp(iatm)
        end do
        deallocate(a_tmp)
        write(*,"(a,i5,a)") " Done!",ntmp," atoms have been cropped"
        
    else if (isel==18) then !Generate randomly displaced geometries
        call displace_geom
    end if
    
    if (isel>=19.and.ifPBC==0) then !The following functions need PBC information
        write(*,"(a)") " Error: To use this function, cell information must be provided by your input file! See Section 2.9.3 of Multiwfn manual for detail"
        write(*,*) "Press ENTER button to continue"
        read(*,*)
        cycle
    end if
    if (isel==19) then !Construct supercell
        nrepli1=1;nrepli2=1;nrepli3=1
        if (ifPBC>=1) then
            write(*,*) "Input number of replicas in direction 1 (i.e. times of current cell), e.g. 3"
            write(*,*) "If inputting 1, the cell will remain unchanged in this direction"
            write(*,*) "If inputting -n, the cell will be duplicated in negative direction by n times"
            read(*,*) nrepli1
        end if
        if (ifPBC>=2) then
            write(*,*) "Input number of replicas in direction 2 (i.e. times of current cell), e.g. 2"
            write(*,*) "If inputting 1, the cell will remain unchanged in this direction"
            write(*,*) "If inputting -n, the cell will be duplicated in negative direction by n times"
            read(*,*) nrepli2
        end if
        if (ifPBC==3) then
            write(*,*) "Input number of replicas in direction 3 (i.e. times of current cell), e.g. 2"
            write(*,*) "If inputting 1, the cell will remain unchanged in this direction"
            write(*,*) "If inputting -n, the cell will be duplicated in negative direction by n times"
            read(*,*) nrepli3
        end if
        allocate(a_old(ncenter))
        a_old=a
        ncenter_old=ncenter
        deallocate(a)
        n1=nrepli1 !Actual number of replicas after processing in direction 1
        if (nrepli1<0) n1=1+abs(nrepli1)
        n2=nrepli2
        if (nrepli2<0) n2=1+abs(nrepli2)
        n3=nrepli3
        if (nrepli3<0) n3=1+abs(nrepli3)
        ncenter=ncenter_old*n1*n2*n3
        allocate(a(ncenter))
        icen=0
        do irepli1=1,n1
            do irepli2=1,n2
                do irepli3=1,n3
                    !i1,i2,i3 are translation times in direction 1,2,3
                    if (nrepli1>=1) then
                        i1=irepli1-1
                    else if (nrepli1<=-1) then
                        if (irepli1==1) then
                            i1=0 !Original cell
                        else
                            i1=-(irepli1-1)
                        end if
                    end if
                    if (nrepli2>=1) then
                        i2=irepli2-1
                    else if (nrepli2<=-1) then
                        if (irepli2==1) then
                            i2=0 !Original cell
                        else
                            i2=-(irepli2-1)
                        end if
                    end if
                    if (nrepli3>=1) then
                        i3=irepli3-1
                    else if (nrepli3<=-1) then
                        if (irepli3==1) then
                            i3=0 !Original cell
                        else
                            i3=-(irepli3-1)
                        end if
                    end if
                    call tvec_PBC(i1,i2,i3,vec)
                    a(icen+1:icen+ncenter_old)=a_old(:)
                    a(icen+1:icen+ncenter_old)%x=a_old(:)%x+vec(1)
                    a(icen+1:icen+ncenter_old)%y=a_old(:)%y+vec(2)
                    a(icen+1:icen+ncenter_old)%z=a_old(:)%z+vec(3)
                    icen=icen+ncenter_old
                end do
            end do
        end do
        cellv1=n1*cellv1
        cellv2=n2*cellv2
        cellv3=n3*cellv3
        nelec=nelec*n1*n2*n3
        deallocate(a_old)
        deallocate(fragatm,intarr)
        allocate(fragatm(ncenter),intarr(ncenter))
        nfragatm=ncenter
        forall (i=1:nfragatm) fragatm(i)=i
        write(*,*) "Done!"
        
    else if (isel==20) then !Make truncated molecules by cell boundary whole
        call makemolwhole
        
    else if (isel==21) then !Scale cell size and atom coordinates
        write(*,*) "0: Scale all cell lengths and atom coordinates correspondingly"
        if (ifPBC>=1) write(*,*) "1: Scale cell length ""a"" and atom coordinates correspondingly"
        if (ifPBC>=2) write(*,*) "2: Scale cell length ""b"" and atom coordinates correspondingly"
        if (ifPBC>=3) write(*,*) "3: Scale cell length ""c"" and atom coordinates correspondingly"
        read(*,*) iscl
        write(*,*) "Input scale factor, e.g. 0.95"
        read(*,*) sfac
        do iatm=1,ncenter
            rcoord(1)=a(iatm)%x
            rcoord(2)=a(iatm)%y
            rcoord(3)=a(iatm)%z
            call Cart2fract(rcoord,fcoord)
            if (iscl==0.or.iscl==1) fcoord(1)=fcoord(1)*sfac
            if (iscl==0.or.iscl==2) fcoord(2)=fcoord(2)*sfac
            if (iscl==0.or.iscl==3) fcoord(3)=fcoord(3)*sfac
            call fract2Cart(fcoord,rcoord)
            a(iatm)%x=rcoord(1)
            a(iatm)%y=rcoord(2)
            a(iatm)%z=rcoord(3)
        end do
        if (iscl==0.or.iscl==1) cellv1=cellv1*sfac
        if (iscl==0.or.iscl==2) cellv2=cellv2*sfac
        if (iscl==0.or.iscl==3) cellv3=cellv3*sfac
        write(*,*) "Done!"
        
    else if (isel==22) then !Wrap all atoms into the center cell
        do iatm=1,ncenter
            rcoord(1)=a(iatm)%x
            rcoord(2)=a(iatm)%y
            rcoord(3)=a(iatm)%z
            call Cart2fract(rcoord,fcoord)
            do idir=1,3
                if (fcoord(idir)>1) then
                    fcoord(idir)=fcoord(idir)-int(fcoord(idir))
                else if (fcoord(idir)<0) then
                    fcoord(idir)=fcoord(idir)+ceiling(abs(fcoord(idir)))
                end if
            end do
            call fract2Cart(fcoord,rcoord)
            a(iatm)%x=rcoord(1)
            a(iatm)%y=rcoord(2)
            a(iatm)%z=rcoord(3)
        end do
        write(*,*) "Done!"
        
    else if (isel==23) then !Translate system along cell axes
        write(*,*) "Input translation distance in direction 1 (in Angstrom), e.g. 0.6"
        write(*,*) "The value could be either positive or negative"
        read(*,*) dist1
        dist2=0
        if (ifPBC>=2) then
            write(*,*) "Input translation distance in direction 2 (in Angstrom), e.g. 0.35"
            write(*,*) "The value could be either positive or negative"
            read(*,*) dist2
        end if
        dist3=0
        if (ifPBC==3) then
            write(*,*) "Input translation distance in direction 3 (in Angstrom), e.g. -0.79"
            write(*,*) "The value could be either positive or negative"
            read(*,*) dist3
        end if
        dist1=dist1/b2a
        dist2=dist2/b2a
        dist3=dist3/b2a
        rnorm1=dsqrt(sum(cellv1**2))
        rnorm2=dsqrt(sum(cellv2**2))
        rnorm3=dsqrt(sum(cellv3**2))
        iout=0
        do iatm=1,ncenter
            a(iatm)%x = a(iatm)%x + dist1*cellv1(1)/rnorm1 + dist2*cellv2(1)/rnorm2 + dist3*cellv3(1)/rnorm3
            a(iatm)%y = a(iatm)%y + dist1*cellv1(2)/rnorm1 + dist2*cellv2(2)/rnorm2 + dist3*cellv3(2)/rnorm3
            a(iatm)%z = a(iatm)%z + dist1*cellv1(3)/rnorm1 + dist2*cellv2(3)/rnorm2 + dist3*cellv3(3)/rnorm3
            rcoord(1)=a(iatm)%x
            rcoord(2)=a(iatm)%y
            rcoord(3)=a(iatm)%z
            call Cart2fract(rcoord,fcoord)
            if (any(fcoord(:)>1)) iout=1
        end do
        write(*,*) "Done!"
        if (iout==1) write(*,"(a)") " Warning: One or more atoms are lying outside the cell, if you want to wrap them into cell, you can use option 22"
    
    else if (isel==24) then !Translate system to center selected part in the cell
        write(*,"(a)") " Input indices of the atoms, whose geometry center will be centered in the cell by translating the whole system. e.g. 1,3,7-10"
        write(*,*) "If you press ENTER button directly, all atoms will be selected"
        read(*,"(a)") c2000tmp
        if (c2000tmp==" ") then
            allocate(tmparr(ncenter))
            ntmp=ncenter
            forall(i=1:ncenter) tmparr(i)=i
        else
            call str2arr(c2000tmp,ntmp)
            allocate(tmparr(ntmp))
            call str2arr(c2000tmp,ntmp,tmparr)
        end if
        fcoord(:)=0.5D0
        call fract2Cart(fcoord,rcoord)
        xcencell=rcoord(1)
        ycencell=rcoord(2)
        zcencell=rcoord(3)
        write(*,"(' X/Y/Z of center of cell:',3f12.6,' Bohr')") xcencell,ycencell,zcencell
        xcenfrag=sum(a(tmparr(:))%x)/ntmp
        ycenfrag=sum(a(tmparr(:))%y)/ntmp
        zcenfrag=sum(a(tmparr(:))%z)/ntmp
        write(*,"(' X/Y/Z of center of selected atoms:',3f12.6,' Bohr')") xcenfrag,ycenfrag,zcenfrag
        rmove1=xcencell-xcenfrag
        rmove2=ycencell-ycenfrag
        rmove3=zcencell-zcenfrag
        iout=0
        do iatm=1,ncenter
            a(iatm)%x = a(iatm)%x + rmove1
            a(iatm)%y = a(iatm)%y + rmove2
            a(iatm)%z = a(iatm)%z + rmove3
            rcoord(1)=a(iatm)%x
            rcoord(2)=a(iatm)%y
            rcoord(3)=a(iatm)%z
            call Cart2fract(rcoord,fcoord)
            if (any(fcoord(:)>1)) iout=1
        end do
        write(*,*) "Done!"
        if (iout==1) write(*,"(a)") " Warning: One or more atoms are lying outside the cell, if you want to wrap them into cell, you can use option 22"
    
    else if (isel==25) then !Extract a molecular cluster (central molecule + neighbouring ones)
        call extract_molclus
        
    else if (isel==26) then !Set cell information
        do while(.true.)
            write(*,*)
            call getcellabc(asize,bsize,csize,alpha,beta,gamma) !The returned a/b/csize are in Angstrom
            write(*,*) "Current cell information:"
            write(*,"(' Cell vector 1,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Angstrom')") cellv1*b2a
            if (ifPBC>1) then
                write(*,"(' Cell vector 2,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Angstrom')") cellv2*b2a
                if (ifPBC>2) write(*,"(' Cell vector 3,  X=',f11.5,'  Y=',f11.5,'  Z=',f11.5,'  Angstrom')") cellv3*b2a
            end if
            if (ifPBC==3) then
                write(*,"(' Cell angles:  Alpha=',f9.4,'  Beta=',f9.4,'  Gamma=',f9.4,' degree')") alpha,beta,gamma
                write(*,"(' Cell size:        a=',f9.4,'     b=',f9.4,'      c=',f9.4,' Angstrom')") dsqrt(sum(cellv1**2))*b2a,dsqrt(sum(cellv2**2))*b2a,dsqrt(sum(cellv3**2))*b2a
                call calc_cellvol(cellvol)
                write(*,"(' Cell volume:',f12.3,' Angstrom^3')") cellvol*b2a**3
            end if
            write(*,*)
            write(*,*) "-1 Restore to original cell information"
            write(*,*) "0 Return"
            if (ifPBC==3) then
                write(*,*) "1 Set length of a"
                write(*,*) "2 Set length of b"
                write(*,*) "3 Set length of c"
                write(*,*) "4 Set angle of alpha (angle between b and c)"
                write(*,*) "5 Set angle of beta (angle between a and c)"
                write(*,*) "6 Set angle of gamma (angle between a and b)"
            end if
            write(*,*) "7 Set cell vector 1"
            if (ifPBC>1) then
                write(*,*) "8 Set cell vector 2"
                if (ifPBC>2) then
                    write(*,*) "9 Set cell vector 3"
                end if
            end if
            if (ifPBC==3) write(*,*) "10 Align cell (make ""a"" along X-axis, and make ""ab"" parallel to XY plane)"
            !if (ifPBC==3) write(*,*) "11 Transform to orthogonal cell"
            read(*,*) isel2
            if (isel2==-1) then
                cellv1=cellv1_org
                cellv2=cellv2_org
                cellv3=cellv3_org
            else if (isel2==0) then
                exit
            end if
            if (isel2>=1.and.isel2<=6) then
                if (isel2==1) then
                    write(*,*) "Input the length of a in Angstrom, e.g. 8.32"
                    read(*,*) asize
                else if (isel2==2) then
                    write(*,*) "Input the length of b in Angstrom, e.g. 8.32"
                    read(*,*) bsize
                else if (isel2==3) then
                    write(*,*) "Input the length of c in Angstrom, e.g. 8.32"
                    read(*,*) csize
                else if (isel2==4) then
                    write(*,*) "Input angle of alpha, e.g. 84.5"
                    read(*,*) alpha
                else if (isel2==5) then
                    write(*,*) "Input angle of beta, e.g. 84.5"
                    read(*,*) beta
                else if (isel2==6) then
                    write(*,*) "Input angle of gamma, e.g. 84.5"
                    read(*,*) gamma
                end if
                call abc2cellv(asize/b2a,bsize/b2a,csize/b2a,alpha,beta,gamma)
            end if
            if (isel2==7) then
                write(*,*) "Input X/Y/Z of cell vector 1 in Angstrom, e.g. 0.84,0.2,0"
                read(*,*) cellv1
                cellv1=cellv1/b2a
            else if (isel2==8) then
                write(*,*) "Input X/Y/Z of cell vector 2 in Angstrom, e.g. 0.2,0.84,0"
                read(*,*) cellv2
                cellv2=cellv2/b2a
            else if (isel2==9) then
                write(*,*) "Input X/Y/Z of cell vector 3 in Angstrom, e.g. 0,0.2,0.84"
                read(*,*) cellv3
                cellv3=cellv3/b2a
            else if (isel2==10) then
                call getcellabc(asize,bsize,csize,alpha,beta,gamma)
                asize=asize/b2a
                bsize=bsize/b2a
                csize=csize/b2a
                call abc2cellv(asize,bsize,csize,alpha,beta,gamma)
                write(*,*) "Done!"
            else if (isel2==11) then !Transform to orthogonal cell. This is useless, because although the volume will be unchanged, the translation symmetry will be broken
                !Make b point to Y axis
                cellv2(1)=0
                cellv2(3)=0
                !Make c point to Z axis
                cellv3(1:2)=0
                write(*,*) "Done!"
            end if
        end do
        
    else if (isel==27) then !Add boundary atoms
        ncenter_old=ncenter
        call construct_atmp_withbound(ncenter_tmp)
        deallocate(a)
        ncenter=ncenter_tmp
        allocate(a(ncenter))
        a=a_tmp
        deallocate(a_tmp)
        ntmp=ncenter-ncenter_old
        if (ntmp==0) then
            write(*,*) "No new boundary atoms were added"
        else
            write(*,"(a,i6,a)") " Done!",ntmp," boundary atoms have been added to present system"
        end if
        
    else if (isel==28) then !Axes interconversion
        write(*,"(a)") " Note: Selected two components of atomic fractional coordinates and cell vector lengths will be exchanged"
        write(*,*)
        write(*,*) "0 Return"
        write(*,*) "1 a<-->b interconversion"
        write(*,*) "2 a<-->c interconversion"
        write(*,*) "3 b<-->c interconversion"
        read(*,*) idir
        if (idir==0) cycle
        !Exchange fractional coordinates
        allocate(fcoordall(3,ncenter))
        do iatm=1,ncenter
            rcoord(1)=a(iatm)%x
            rcoord(2)=a(iatm)%y
            rcoord(3)=a(iatm)%z
            call Cart2fract(rcoord,fcoord)
            if (idir==1) then
                tmpval=fcoord(1)
                fcoord(1)=fcoord(2)
                fcoord(2)=tmpval
            else if (idir==2) then
                tmpval=fcoord(1)
                fcoord(1)=fcoord(3)
                fcoord(3)=tmpval
            else if (idir==3) then
                tmpval=fcoord(2)
                fcoord(2)=fcoord(3)
                fcoord(3)=tmpval
            end if
            fcoordall(:,iatm)=fcoord(:)
        end do
        !Exchange cell lengths
        anorm=dsqrt(sum(cellv1**2))
        bnorm=dsqrt(sum(cellv2**2))
        cnorm=dsqrt(sum(cellv3**2))
        if (idir==1) then
            cellv1=cellv1/anorm*bnorm
            cellv2=cellv2/bnorm*anorm
        else if (idir==2) then
            cellv1=cellv1/anorm*cnorm
            cellv3=cellv3/cnorm*anorm
        else if (idir==3) then
            cellv2=cellv2/bnorm*cnorm
            cellv3=cellv3/cnorm*bnorm
        end if
        !Obtain new atomic coordinates
        do iatm=1,ncenter
            call fract2Cart(fcoordall(:,iatm),rcoord)
            a(iatm)%x=rcoord(1)
            a(iatm)%y=rcoord(2)
            a(iatm)%z=rcoord(3)
        end do
        deallocate(fcoordall)
        write(*,*) "Done!"
    end if
end do
    
end subroutine



!!---------- Extract a molecular cluster (central molecule + neighbouring ones)
!The botteneck of cost is call genconnmat(1,0). In fact the cost of this routine can be reduced &
!if directly constructing neighbouring list and then atmfrg
subroutine extract_molclus
use defvar
use util
implicit real*8 (a-h,o-z)
type(atomtype),allocatable :: a_old(:)
character c2000tmp*2000,c80tmp*80
real*8 vec(3)
integer,allocatable :: cenidx(:),atmfrg(:),iffrag(:),idxarr(:)
integer,allocatable :: iext(:) !=1 need to extract, =0 do not need

write(*,"(/,a)") " Input index of an atom, e.g. 5. The whole molecule that the atom belongs to as well as the neighbouring molecules will be extracted"
read(*,*) iselatm
write(*,*) "Input threshold of detecting contact"
write(*,"(a)") " If you input e.g. 1.5, then if the closest distance between a neighbouring molecule &
&and current molecule is shorter than the sum of the corresponding atomic vdW radii multiplied by 1.5, the neighbouring molecule will be extracted"
write(*,*) "If pressing ENTER button directly, 1.2 will be used"
write(*,"(a)") " If you simply want to get the single molecule containing the atom you selected, input 0"
read(*,"(a)") c80tmp
if (c80tmp==" ") then
    crit=1.2D0
else
    read(c80tmp,*) crit
end if

allocate(a_old(ncenter))
a_old=a
ncenter_old=ncenter
deallocate(a)
nelec=0

call walltime(iwalltime1)

write(*,*) "Constructing 5*5*5 supercell..."
!Temporarily replicate original cell in each direction as (-2:2), where 0 is current cell
!This guarantee that we can extract a cluster containing selected molecule and all whole molecules surrounding it
!The original atoms keep their original indices
ncenter=ncenter_old*5*5*5
allocate(a(ncenter))
a(1:ncenter_old)=a_old
icen=ncenter_old
do i1=-2,2
    do i2=-2,2
        do i3=-2,2
            if (i1==0.and.i2==0.and.i3==0) cycle !Original cell
            call tvec_PBC(i1,i2,i3,vec)
            a(icen+1:icen+ncenter_old)=a_old(:)
            a(icen+1:icen+ncenter_old)%x=a_old(:)%x+vec(1)
            a(icen+1:icen+ncenter_old)%y=a_old(:)%y+vec(2)
            a(icen+1:icen+ncenter_old)%z=a_old(:)%z+vec(3)
            icen=icen+ncenter_old
        end do
    end do
end do
deallocate(a_old)

!PBC information is useless anymore
ifPBC=0
cellv1=0
cellv2=0
cellv3=0

write(*,*) "Generating fragments according to connectivity..."
allocate(atmfrg(ncenter))
call genconnfrag(atmfrg) !atmfrg records fragment index of every atom
ncenidx=count(atmfrg==atmfrg(iselatm)) !Number of atoms in central molecule
allocate(iext(ncenter),cenidx(ncenidx))
iext=0
itmp=0
do iatm=1,ncenter
    if (atmfrg(iatm)==atmfrg(iselatm)) then 
        itmp=itmp+1
        cenidx(itmp)=iatm !cenidx array records indices of central atoms
        iext(iatm)=1 !This atom will be extracted
    end if
end do

!Determine which atoms should be kepted
write(*,*) "Determining neighbouring molecules..."
do iatm=1,ncenter
    if (iext(iatm)==1) cycle !This atom is already to be extracted
    vdwr_i=vdwr(a(iatm)%index)
    do idx=1,ncenidx !Cycle atoms in central molecule
        jatm=cenidx(idx)
        sumvdwrad=vdwr_i+vdwr(a(jatm)%index)
        if (atomdist(iatm,jatm,0)<sumvdwrad*crit) then !If an atom is close enough to an atom in central molecule, all atoms having the same fragment index as it will be extracted
            iatmfrg=atmfrg(iatm)
            do katm=1,ncenter
                if (atmfrg(katm)==iatmfrg) iext(katm)=1
            end do
            exit
        end if
    end do
end do

allocate(a_old(ncenter))
a_old=a
ncenter_old=ncenter
ncenter=count(iext==1)
write(*,"(a,i8)") " Done! Number of extracted atoms:",ncenter
deallocate(a)
allocate(a(ncenter))
icen=0
do iatm=1,ncenter_old
    if (iext(iatm)==1) then
        icen=icen+1
        a(icen)=a_old(iatm)
        if (iatm==iselatm) iselnew=icen
    end if
end do

call walltime(iwalltime2)
write(*,"(' Extraction totally took up wall clock time',i10,' s')") iwalltime2-iwalltime1

write(*,"(/,' Index of you previously selected atom in current system:',i8)") iselnew
allocate(iffrag(ncenter),idxarr(ncenter))
call getfragatoms(iselnew,iffrag)
icen=0
do iatm=1,ncenter
    if (iffrag(iatm)==1) then
        icen=icen+1
        idxarr(icen)=iatm
    end if
end do
call arr2str_2(idxarr(1:icen),c2000tmp)
write(*,*) "Atoms in the central molecule:"
write(*,"(a)") trim(c2000tmp)

!Finalize
deallocate(fragatm)
allocate(fragatm(ncenter))
nfragatm=ncenter
forall (i=1:nfragatm) fragatm(i)=i
end subroutine



!!----------- Plot surface distance projection map
subroutine molsurf_distmap
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
real*8 :: rhoiso=0.05D0,zstep=0.05D0/b2a
integer :: isurftype=1
character c80tmp*80,c80tmp2*80

nx=300
ny=300
xlow=minval(a%x)-3
xhigh=maxval(a%x)+3
ylow=minval(a%y)-3
yhigh=maxval(a%y)+3
zstart=maxval(a%z)
zend=minval(a%z)
if (abs(zstart-zend)<0.1D0) zstart=zstart+2/b2a !For molecule in exactly XY plane, zstart is higher than the plane by 2 Angstrom

!Plotting settings
plesel=1 !XY plane
planestpx=1
planestpy=1
planestpz=1
ilenunit2D=2
numdigx=0
numdigy=0
numdigz=1
disshowlabel=2
iatom_on_plane=1 !Show atom label on map
iclrindatmlab=3 !Blue label
iatom_on_plane_far=1
iclrtrans=3
ctrnegstyle=1;ctrnegstyle=0 !Use solid line for negative contour
ngenctr=25
idrawtype=1
idrawcontour=1

if (.not.allocated(a)) then
    write(*,"(a)") " Error: Your input file must contain atom information! Please read Section 2.5 of Multiwfn manual carefully!"
    write(*,*) "Press ENTER button to return"
    read(*,*)
    return
end if

do while(.true.)
    write(*,*)
    write(*,*) " ====================== Surface distance projection map ====================="
    write(*,*) "-1 Return"
    write(*,*) " 0 Start calculation"
    if (isurftype==1) write(*,"(a,f9.5,' a.u.')") "  1 Set definition of system surface, current: Promolecular density with isovalue of",rhoiso
    if (isurftype==2) write(*,"(a,f9.5,' a.u.')") "  1 Set definition of system surface, current: Electron density with isovalue of",rhoiso
    if (isurftype==3) write(*,"(a,f7.3)") "  1 Set definition of system surface, current: Superposition of vdW spheres multiplied by",vdwsphscl
    write(*,"(a,f8.3,a,f8.3,a)") "  2 Choose range of X axis, current:",xlow*b2a," to",xhigh*b2a," Angstrom"
    write(*,"(a,f8.3,a,f8.3,a)") "  3 Choose range of Y axis, current:",ylow*b2a," to",yhigh*b2a," Angstrom"
    write(*,"(a,i4)") "  4 Set number of grids in X, current:",nx
    write(*,"(a,i4)") "  5 Set number of grids in Y, current:",ny
    write(*,"(a,f8.4,' Angstrom')") "  6 Set stepsize for detection isosurface in Z, current:",zstep*b2a
    write(*,"(a,f9.4,' Angstrom')") "  7 Set Z position of starting detection, current:",zstart*b2a
    write(*,"(a,f9.4,' Angstrom')") "  8 Set Z position of ending detection, current:",zend*b2a
    write(*,"(a,i5)") "  9 Set number of contour lines to generate, current:",ngenctr
    read(*,*) isel
    
    if (isel==-1) then
        return
    else if (isel==1) then
        write(*,*) "Use which definition of system surface?"
        write(*,*) "1 Promolecular electron density"
        if (allocated(b)) write(*,*) "2 Electron density calculated based on wavefunction"
        write(*,*) "3 Superposition of atomic van der Waals spheres scaled by a factor"
        read(*,*) isurftype
        if (isurftype==1.or.isurftype==2) then
            write(*,*) "Input isovalue of electron density in a.u., e.g. 0.01"
            read(*,*) rhoiso
        else if (isurftype==3) then
            write(*,*) "Input scale factor for vdW radii, e.g. 0.8"
            read(*,*) vdwsphscl
        end if
    else if (isel==2) then
        write(*,*) "Input lower and upper limits of X in Angstrom, e.g. -2.6,8.3"
        read(*,*) xlow,xhigh
        xlow=xlow/b2a
        xhigh=xhigh/b2a
    else if (isel==3) then
        write(*,*) "Input lower and upper limits of Y in Angstrom, e.g. -2.6,8.3"
        read(*,*) ylow,yhigh
        ylow=ylow/b2a
        yhigh=yhigh/b2a
    else if (isel==4) then
        write(*,*) "Input number of grids in X, e.g. 150"
        read(*,*) nx
    else if (isel==5) then
        write(*,*) "Input number of grids in Y, e.g. 150"
        read(*,*) ny
    else if (isel==6) then
        write(*,*) "Input stepsize for detecting isosurface in Z (in Angstrom), e.g. 0.01"
        read(*,*) zstep
        zstep=zstep/b2a
    else if (isel==7) then
        write(*,*) "Input Z position of starting detection (in Angstrom), e.g. 6.25"
        write(*,*) "You can also input e.g. ""a 5"" to use Z coordinate of atom 5 as this position"
        read(*,"(a)") c80tmp
        if (index(c80tmp,'a')/=0) then
            read(c80tmp,*) c80tmp2,iatm
            zstart=a(iatm)%z
        else
            read(c80tmp,*) zstart
            zstart=zstart/b2a
        end if
    else if (isel==8) then
        write(*,*) "Input Z position of ending detection (in Angstrom), e.g. -3.2"
        read(*,*) zend
        zend=zend/b2a
    else if (isel==9) then
        write(*,*) "Number of contour lines to generate, e.g. 30"
        read(*,*) ngenctr
        
    else if (isel==0) then !Start calculation!
        allocate(planemat(nx,ny))
        dx=(xhigh-xlow)/(nx-1)
        dy=(yhigh-ylow)/(ny-1)
        call walltime(iwalltime1)
        call showprog(0,nx)
        planemat=zend-zstart
        ifinish=0
        !$OMP PARALLEL DO SHARED(planemat,ifinish) PRIVATE(ix,iy,xtmp,ytmp,ztmp,tmpval,ifin,tmpvallast,ztmplast) schedule(dynamic) NUM_THREADS(nthreads)
        do ix=1,nx
            xtmp=xlow+dx*(ix-1)
            do iy=1,ny
                ytmp=ylow+dy*(iy-1)
                ztmp=zstart
                do while(ztmp>=zend)
                    if (isurftype==1.or.isurftype==2) then
                        if (isurftype==1) then
                            tmpval=calcprodens(xtmp,ytmp,ztmp,0)
                        else if (isurftype==2) then
                            tmpval=fdens(xtmp,ytmp,ztmp)
                        end if
                        if (tmpval>rhoiso) then
                            if (ztmp==zstart) then
                                planemat(ix,iy)=0
                            else !Use linear interpolation to improve smoothness
                                planemat(ix,iy)=ztmp+(rhoiso-tmpval)/(tmpvallast-tmpval)*zstep - zstart
                            end if
                            exit
                        end if
                        tmpvallast=tmpval
                        ztmplast=ztmp
                    else if (isurftype==3) then
                        call inside_vdWsph(xtmp,ytmp,ztmp,vdwsphscl,ifin)
                        if (ifin==1) then
                            planemat(ix,iy)=ztmp-zstart
                            exit
                        end if
                    end if
                    ztmp=ztmp-zstep
                end do
            end do
            !$OMP CRITICAL
            ifinish=ifinish+1
            call showprog(ifinish,nx)
            !$OMP END CRITICAL
        end do
        !$OMP END PARALLEL DO
        planemat=planemat*b2a
        call walltime(iwalltime2)
        write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
        
        ngridnum1=nx
        ngridnum2=ny
        orgz2D=zstart
        clrhigh=0
        clrlow=ceiling((zend-zstart)*b2a)-1E-5
        call gencontour(2,clrlow,clrhigh,25) !Generate contour lines evenly covering lower and upper limits
        call planemap_interface("surface projection","distmap",xlow,xhigh,ylow,yhigh,clrlow,clrhigh)
        deallocate(planemat)
    end if
end do
end subroutine

!------- Test if a point is inside molecular vdW surface defined by superposition of Bondi vdW spheres scaled by vdwsphscl
!ifin=1: inside  ifin=0: outside
subroutine inside_vdWsph(xtmp,ytmp,ztmp,vdwsphscl,ifin)
use defvar
implicit real*8 (a-h,o-z)
real*8 xtmp,ytmp,ztmp,vdwsphscl,xyzA(3),xyzB(3)
integer ifin
ifin=0
do iatm=1,ncenter
	rx=a(iatm)%x-xtmp
	ry=a(iatm)%y-ytmp
	rz=a(iatm)%z-ztmp
    if (ifPBC==0) then
    	dist=dsqrt(rx*rx+ry*ry+rz*rz)
    else
        xyzA(1)=xtmp;xyzA(2)=ytmp;xyzA(3)=ztmp
        xyzB(1)=a(iatm)%x;xyzB(2)=a(iatm)%y;xyzB(3)=a(iatm)%z
        call nearest_dist(xyzA,xyzB,dist)
    end if
    if (dist<vdwr(a(iatm)%index)*vdwsphscl) then
        ifin=1
        return
    end if
end do
end subroutine




!!-------- Calculate Fermi-level based on Fermi-Dirac distribution using energies of molecular orbitals
subroutine calc_Fermi_level
use defvar
implicit real*8 (a-h,o-z)
character c80tmp*80
real*8 numelec,numelec_tmp
real*8 :: thres=1D-6

if (allocated(b)) then
    if (wfntype==3.or.wfntype==4) then
        write(*,"(a)") " Note: Because orbital occupation numbers are not all integer, now make orbital occupations integer and satisfy Aufbau principle"
        call make_occ_integer_Aufbau
    end if
    call getHOMOidx
    if (idxHOMO==nmo) then
        write(*,*)
        write(*,*) "Error: This function cannot be used because there is no unoccupied orbitals!"
        return
    end if
else !Load orbital energies from plain text file
    if (allocated(MOene)) deallocate(MOene)
    open(10,file=filename,status="old")
    read(10,"(a)") c80tmp
    read(c80tmp,*,iostat=ierror) nocc,nvir,noccB,nvirB
    if (ierror==0) then !Unrestricted
        wfntype=1
        nelec=nocc+noccB
        nmo=nocc+nvir+noccB+nvirB
        idxHOMO=nocc
        idxHOMOb=nocc+nvir+noccB
    else !Restricted
        read(c80tmp,*) nocc,nvir
        wfntype=0
        nelec=2*nocc
        nmo=nocc+nvir
        idxHOMO=nocc
    end if
    allocate(MOene(nmo))
    read(10,*) MOene(1:nocc)
    read(10,*) MOene(nocc+1:nocc+nvir)
    if (wfntype==1) then
        read(10,*) MOene(nocc+nvir+1:nocc+nvir+noccB)
        read(10,*) MOene(nocc+nvir+noccB+1:nocc+nvir+noccB+nvirB)
    end if
    close(10)
    write(*,*) "Loading orbital energies from input file finished!"
end if

!call nelec_Ef_T(0.14448977718533D0,300D0,numelec)
!write(*,*) numelec

do while(.true.)
    write(*,*)
    write(*,*) "Input temperature (K) for determining Fermi level, e.g. 400"
    write(*,*) "If press ENTER button directly, 298.15 K will be used"
    write(*,*) "Input ""q"" can exit. Input ""thres"" can modify convergence threshold"
    read(*,"(a)") c80tmp
    if (index(c80tmp,'q')/=0) exit
    if (c80tmp==" ") then
        T=298.15D0
    else if (c80tmp=="thres") then
        write(*,*) "Input threshold, e.g. 1E-8"
        read(*,*) thres
        cycle
    else
        read(c80tmp,*) T
    end if
    if (T==0) then
        write(*,"(a)") " Error: Fermi level is ill-defined at 0 K! Any value between E_HOMO and E_LUMO may be acceptable. Press ENTER button to return"
        read(*,*)
        cycle
    end if
    !Use bisection method to determine Fermi level
    if (wfntype==0) then
        vallow=MOene(idxHOMO)
        valhigh=MOene(idxHOMO+1)
    else
        vallow=min(MOene(idxHOMO),MOene(idxHOMOb))
        valhigh=max(MOene(idxHOMO+1),MOene(idxHOMOb+1))
    end if
    if (valhigh==0) then
        write(*,"(a)") " Warning: This file does not contain actually calculated unoccupied orbitals! The result is meaningless!"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        exit
    end if
    iter=0
    write(*,*) "Starting bisection method to determine Fermi level"
    write(*,"(a,f10.1)") " Expected number of electrons:",nelec
    write(*,"(a,1PE10.1)") " Convergence threshold of deviation to expected number of electrons:",thres
    write(*,"(a,2f14.8,a)") " Initial lower and upper limits:",vallow,valhigh," Hartree"
    call nelec_Ef_T(vallow,T,tmp_low)
    call nelec_Ef_T(valhigh,T,tmp_high)
    write(*,"(a,2f18.8)") " Corresponding number of electrons:",tmp_low,tmp_high
    write(*,*)
    do while(.true.)
        iter=iter+1
        Ef=(vallow+valhigh)/2D0
        call nelec_Ef_T(Ef,T,numelec)
        write(*,"(' Iter:',i5,'  Nelec:',f16.8,'  Dev.:',D13.5,'  Ef:',f14.8,' a.u.')") iter,numelec,numelec-nelec,Ef
        if (iter>=1000) then
            write(*,*) "Error: Number of electrons failed to converge to 1E-6 within 1000 iterations!"
            exit
        else if (abs(numelec-nelec)<thres) then
            write(*,"(' Converged! Fermi level is',f14.8,' Hartree',f14.6,' eV')") Ef,Ef*au2eV
            exit
        end if
        if (numelec>nelec) then
            valhigh=Ef
            call nelec_Ef_T(vallow,T,numelec_tmp)
            if (numelec_tmp>nelec) then
                vallow=vallow-0.01D0
                write(*,"(' Decrease lower limit by 0.01 to ',f16.8,' a.u.')") vallow
            end if
        else
            vallow=Ef
            call nelec_Ef_T(valhigh,T,numelec_tmp)
            if (numelec_tmp<nelec) then
                valhigh=valhigh+0.01D0
                write(*,"(' Increase upper limit by 0.01 to ',f16.8,' a.u.')") valhigh
            end if
        end if
    end do
end do
end subroutine
!--------- Get number of electrons at given Fermi level (a.u.) and temperature (K) based on Fermi-Dirac distribution using current orbitals
subroutine nelec_Ef_T(Ef,T,numelec)
use defvar
real*8 occtmp(nmo),Ef,T,numelec
occtmp=0
do imo=1,nmo
    if (MOene(imo)==0) cycle !Skip falsely filled MOs
    !write(*,"(i5,4E16.8)") imo,(MOene(imo)-Ef),(boltzcau*T),(MOene(imo)-Ef)/(boltzcau*T),1D0/( 1+exp( (MOene(imo)-Ef)/(boltzcau*T) ) )
    occtmp(imo)= 1D0/( 1+exp( (MOene(imo)-Ef)/(boltzcau*T) ) )
    !write(*,"(i5,f16.8)") imo,occtmp(imo)
end do
numelec=sum(occtmp)
if (wfntype==0.or.wfntype==2) numelec=numelec*2 !For spatial orbitals, the occupation should be doubled
end subroutine