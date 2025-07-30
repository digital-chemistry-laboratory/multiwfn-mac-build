!-------- Main interface of visual study of weak interaction
subroutine visweak_main
use defvar
implicit real*8 (a-h,o-z)

write(*,"(a)") " Reviews of functions in this module:"
write(*,"(a)") " Tian Lu and Qinxue Chen, Visualization Analysis of Weak Interactions in Chemical Systems. &
&In Comprehensive Computational Chemistry, vol. 2, pp. 240-264. Oxford: Elsevier (2024) DOI: 10.1016/B978-0-12-821978-2.00076-3"
write(*,"(a)") " Tian Lu, Visualization Analysis of Covalent and Noncovalent Interactions in Real Space, Angew. Chem. Int. Ed., 137, e202504895 (2025) DOI: 10.1002/anie.202504895"
do while(.true.)
	write(*,*)
	write(*,*) "           ============ Visual study of weak interaction ============ "
	write(*,*) " 0 Return"
	write(*,*) " 1 NCI analysis (also known as RDG analysis)"
	write(*,*) " 2 NCI analysis based on promolecular density"
	write(*,*) " 3 aNCI: Averaged NCI analysis"
	write(*,*) " 4 IRI: Interaction region indicator analysis (Chemistry-Methods, 1, 231)"
	write(*,*) " 5 DORI analysis"
	write(*,*) " 6 Visualization of van der Waals potential (JMM, 26, 315)"
	write(*,*) " 9 Becke/Hirshfeld surface analysis"
	write(*,*) " 10 IGM analysis"
	write(*,*) "-10 mIGM: Modified IGM analysis"
	write(*,"(a)") "  11 IGMH: IGM analysis based on Hirshfeld partition of molecular density (JCC, 43, 539)"
	write(*,*) " 12 aIGM: Averaged IGM analysis"
	write(*,*) "-12 amIGM: Averaged mIGM analysis"
	read(*,*) isel
    
	if ((isel==1.or.isel==4.or.isel==5.or.isel==11).and.(.not.allocated(b))) then
		write(*,*) "Error: Wavefunction information is not available!"
        write(*,"(a)") " You must use a file containing wavefunction information (e.g. .wfn/fch/mwfn/wfx...) as input file to perform the analysis. &
        &See Section 2.5 of Multiwfn manual to understand which files contain wavefunction information"
        write(*,*) "Press ENTER button to return"
        read(*,*)
        return
    end if
	if (isel==0) then
		return
	else if (isel==1) then
		call funcvsfunc(1)
	else if (isel==2) then
		call funcvsfunc(2)
	else if (isel==3) then
		call aNCI
	else if (isel==4) then
		call funcvsfunc(4)
	else if (isel==5) then
		call funcvsfunc(5)
    else if (isel==6) then
        call vdWpotential
	else if (isel==9) then
		write(*,"(a)") " Note: To perform Becke or Hirshfeld surface analysis, you should use main function 12. &
		&Please check Section 3.15.5 of Multiwfn manual on how to do that. Corresponding examples are given as Sections 4.12.5 and 4.12.6."
		write(*,*) "Press ENTER button to continue"
		read(*,*)
	else if (isel==10) then
		call IGM(1)
	else if (isel==-10) then
		call IGM(-1)
	else if (isel==11) then
		call IGM(2)
	else if (isel==12) then
		call aIGM(1)
	else if (isel==-12) then
		call aIGM(-1)
	end if
end do
end subroutine



!!-------- Averaged RDG/NCI analysis for MD, use promolecular approximation
subroutine aNCI
use defvar
use util
use GUI
use functions
implicit real*8 (a-h,o-z)
!The first index of avggrad and the first two indices of avghess correspond to components of gradient and Hessian, respectively
real*8,allocatable :: avgdens(:,:,:),avggrad(:,:,:,:),avghess(:,:,:,:,:)
real*8,allocatable :: avgRDG(:,:,:),avgsl2r(:,:,:),thermflu(:,:,:)
real*8,allocatable :: scatterx(:),scattery(:)

write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
write(*,*) "  Original paper of aNCI: J. Chem. Theory Comput., 9, 2226 (2013)"
write(*,"(a)") "   Review: Tian Lu, Qinxue Chen, Visualization Analysis of &
&Weak Interactions in Chemical Systems. In Comprehensive Computational Chemistry, vol. 2, pp. 240-264. Oxford: Elsevier (2024) DOI: 10.1016/B978-0-12-821978-2.00076-3"

write(*,*)
write(*,"(a)") " NOTE: amIGM is the significantly better choice than aNCI! Please consider to use amIGM instead"

write(*,*)
write(*,*) "Input range of the frames to be analyzed, e.g. 150,400 means 150 to 400 frames"
write(*,*) "Note: The frame index starts from 1"
read(*,*) ifpsstart,ifpsend
nfps=ifpsend-ifpsstart+1
write(*,"(' Selected',i8,' frames, frames from',i8,' to',i8,' will be processed',/)") nfps,ifpsstart,ifpsend

call setgrid(0,igridsel)

call walltime(iwalltime1)
write(*,"(/,a)") " Calculating averaged density, gradient and Hessian of density"
allocate(avgdens(nx,ny,nz),avggrad(3,nx,ny,nz),avghess(3,3,nx,ny,nz))
call avg_rhogradhess(avgdens,avggrad,avghess,ifpsstart,ifpsend)
write(*,"(a)") " Calculating averaged RDG and averaged sign(lambda2)*rho..."
allocate(avgRDG(nx,ny,nz),avgsl2r(nx,ny,nz))
call avg_RDG_sl2r(avgdens,avggrad,avghess,avgRDG,avgsl2r)
call walltime(iwalltime2)
write(*,"(' Calculation totally took up wall clock time',i10,' s')") iwalltime2-iwalltime1
deallocate(avggrad,avghess) !Will not be used further, so release its memory

allocate(scatterx(nx*ny*nz),scattery(nx*ny*nz))
ii=1
do k=1,nz
	do j=1,ny
		do i=1,nx
			scatterx(ii)=avgsl2r(i,j,k)
			scattery(ii)=avgRDG(i,j,k)
			ii=ii+1
		end do
	end do
end do
!Default axis range of scatter plot
xmin=-RDGprodens_maxrho
xmax=RDGprodens_maxrho
if (RDGprodens_maxrho==0D0) xmin=-2D0
if (RDGprodens_maxrho==0D0) xmax=2D0
ymin=0D0
ymax=1D0

do while (.true.)
    write(*,*)
	write(*,*) "0 Return"
	write(*,*) "1 Draw scatter graph between averaged RDG and sign(lambda2)*rho"
	write(*,*) "2 Save the scatter graph to file"
	write(*,*) "3 Change range of X-axis of the scatter graph"
	write(*,*) "4 Change range of Y-axis of the scatter graph"
	write(*,*) "5 Export scatter points to output.txt in current folder"
	write(*,"(a)") " 6 Export averaged RDG and sign(lambda2)*rho to avgRDG.cub and avgsl2r.cub in current folder, respectively"
	write(*,"(a)") " 7 Compute thermal fluctuation index (TFI) and export to thermflu.cub in current folder"
	write(*,*) "8 Show isosurface of averaged RDG"
    write(*,*) "9 Export averaged density to avgdens.cub in current folder"
	read(*,*) isel
	
	if (isel==0) then
		exit
	else if (isel==1) then
		write(*,*) "Drawing graph, please wait..."
		call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"Averaged $sign({\lambda}_2)\rho$ (a.u.)","Averaged reduced density gradient")
	else if (isel==2) then
		isavepic=1
		call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"Averaged $sign({\lambda}_2)\rho$ (a.u.)","Averaged reduced density gradient")
		isavepic=0
		write(*,"(a,a,a)") " Graph have been saved to ",trim(graphformat)," file with ""dislin"" prefix in current directory"
	else if (isel==3) then
		write(*,*) "Input lower limit and upper limit of X axis e.g. 0,1.5"
		read(*,*) xmin,xmax
	else if (isel==4) then
		write(*,*) "Input lower limit and upper limit of Y axis e.g. 0,1.5"
		read(*,*) ymin,ymax
	else if (isel==5) then
		write(*,*) "Outputting output.txt in current folder..."
		open(10,file="output.txt",status="replace")
        do k=1,nz
			do j=1,ny
				do i=1,nx
					call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
					write(10,"(3f12.6,2E16.8)") tmpx*b2a,tmpy*b2a,tmpz*b2a,avgsl2r(i,j,k),avgRDG(i,j,k)
				end do
			end do
			call showprog(k,nz)
		end do
		close(10)
		write(*,*) "Finished!"
        write(*,*) "Column 1/2/3: X/Y/Z in Angstrom"
        write(*,*) "Column 4/5: Averaged sign(lambda2)*rho in a.u. and averaged RDG"
	else if (isel==6) then
		write(*,*) "Outputting averaged reduced density gradient to avgRDG.cub in current folder"
		open(10,file="avgRDG.cub",status="replace")
		call outcube(avgRDG,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
		write(*,*)
		write(*,*) "Outputting averaged sign(lambda2)*rho to avgsl2r.cub in current folder"
		open(10,file="avgsl2r.cub",status="replace")
		call outcube(avgsl2r,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
	else if (isel==7) then
        call calcexport_TFI(avgdens,ifpsstart,ifpsend)
	else if (isel==8) then
		write(*,*) "Please wait..."
	 	sur_value=0.3D0
	 	if (allocated(cubmat)) deallocate(cubmat)
	 	allocate(cubmat(nx,ny,nz))
	 	cubmat=avgRDG
		call drawisosurgui(1)
		deallocate(cubmat)
    else if (isel==9) then
		write(*,*) "Outputting averaged density to avgdens.cub in current folder"
		open(10,file="avgdens.cub",status="replace")
		call outcube(avgdens,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "Done!"
	end if
end do

end subroutine




!!----- Generate averaged density, averaged gradient and averaged Hessian
subroutine avg_rhogradhess(avgdens,avggrad,avghess,ifpsstart,ifpsend)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8 avgdens(nx,ny,nz),avggrad(3,nx,ny,nz),avghess(3,3,nx,ny,nz)
real*8 gradtmp(3),hesstmp(3,3)
avgdens=0D0
avggrad=0D0
avghess=0D0
nfps=ifpsend-ifpsstart+1
open(10,file=filename,status="old")
do ifps=1,ifpsend
	!call readxyz(filename,1,0)
	call readxyztrj(10)
	if (ifps<ifpsstart) cycle
    call showprog(ifps,nfps)
	!Global variable "fragatm" must be defined for each frame, because each frame may have different content, and calchessmat_prodens will use it
	nfragatm=ncenter
	if (allocated(fragatm)) deallocate(fragatm)
	allocate(fragatm(nfragatm))
    forall(itmp=1:nfragatm) fragatm(itmp)=itmp
	!$OMP PARALLEL DO SHARED(avgdens,avggrad,avghess) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp,gradtmp,hesstmp) schedule(dynamic) NUM_THREADS(nthreads)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
				call calchessmat_prodens(tmpx,tmpy,tmpz,denstmp,gradtmp,hesstmp)
				avgdens(i,j,k)=avgdens(i,j,k)+denstmp
				avggrad(:,i,j,k)=avggrad(:,i,j,k)+gradtmp
				avghess(:,:,i,j,k)=avghess(:,:,i,j,k)+hesstmp
			end do
		end do
	end do
	!$OMP END PARALLEL DO
end do
close(10)
avgdens=avgdens/nfps
avggrad=avggrad/nfps
avghess=avghess/nfps
end subroutine



!!----- Use averaged density, averaged gradient and averaged Hessian to compute averaged RDG and averaged Sign(lambda2)*rho
subroutine avg_RDG_sl2r(avgdens,avggrad,avghess,avgRDG,avgsl2r)
use defvar
use util
implicit real*8 (a-h,o-z)
real*8 avgdens(nx,ny,nz),avggrad(3,nx,ny,nz),avghess(3,3,nx,ny,nz),avgRDG(nx,ny,nz),avgsl2r(nx,ny,nz)
real*8 eigvecmat(3,3),eigval(3)
!$OMP PARALLEL DO SHARED(avgRDG,avgsl2r) PRIVATE(i,j,k,tmpx,tmpy,tmpz,avggradnormtmp,eigval) schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)

			avggradnormtmp=dsqrt(sum(avggrad(:,i,j,k)**2))
			if (RDGprodens_maxrho/=0D0.and.avgdens(i,j,k)>=RDGprodens_maxrho) then
				avgRDG(i,j,k)=100D0
			else if (avggradnormtmp==0D0.or.avgdens(i,j,k)==0D0) then
				avgRDG(i,j,k)=999D0
			else
				avgRDG(i,j,k)=0.161620459673995D0*avggradnormtmp/avgdens(i,j,k)**(4D0/3D0) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
			end if
			
			call diagmat(avghess(:,:,i,j,k),eigvecmat,eigval,100,1D-6)
			call sort(eigval)
			if (eigval(2)/=0D0) then
				avgsl2r(i,j,k)=avgdens(i,j,k)*eigval(2)/abs(eigval(2)) !At nuclei of single atom system, Hessian returned may be zero matrix
			else
				avgsl2r(i,j,k)=-avgdens(i,j,k) !Around nuclei, eigval(2)/abs(eigval(2)) always be negative
			end if
			
		end do
	end do
end do
!$OMP END PARALLEL DO
end subroutine



!!-------- Calculate and export thermal fluctuation index (TFI)
!avgdens is inputted averaged electron density
subroutine calcexport_TFI(avgdens,ifpsstart,ifpsend)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
real*8 thermflu(nx,ny,nz),avgdens(nx,ny,nz)

nfps=ifpsend-ifpsstart+1
call walltime(iwalltime1)
write(*,*) "Calculating thermal fluctuation index..."
thermflu=0D0
open(10,file=filename,status="old")
!Calculate fluctuation square term of density and temporarily store it to thermflu array
do ifps=1,ifpsend
	!call readxyz(filename,1,0)
	call readxyztrj(10)
	if (ifps<ifpsstart) cycle
    call showprog(ifps,nfps)
	!$OMP PARALLEL DO SHARED(thermflu) PRIVATE(i,j,k,tmpx,tmpy,tmpz,denstmp) schedule(dynamic) NUM_THREADS(nthreads)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
				call calchessmat_prodens(tmpx,tmpy,tmpz,denstmp)
				thermflu(i,j,k)=thermflu(i,j,k)+(denstmp-avgdens(i,j,k))**2
			end do
		end do
	end do
	!$OMP END PARALLEL DO
end do
close(10)
do k=1,nz
	do j=1,ny
		do i=1,nx
            if (avgdens(i,j,k)/=0) then
                thermflu(i,j,k)=dsqrt(thermflu(i,j,k)/nfps)/avgdens(i,j,k)
            else
                thermflu(i,j,k)=0D0
            end if
		end do
	end do
end do
call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,' s',/)") iwalltime2-iwalltime1
write(*,*) "Outputting thermal fluctuation index to thermflu.cub in current folder"
open(10,file="thermflu.cub",status="replace")
call outcube(thermflu,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
close(10)
write(*,*) "Done!"
end subroutine




!!------------------------------------------------------------------------------------------------
!! ----------- Independent Gradient Model (IGM) analysis based on promolecular density -----------
!!------------------------------------------------------------------------------------------------
!iIGMtype=1: Based on promolecular approximation (original IGM)
!iIGMtype=-1: Modified IGM (mIGM), namely based on Hirshfeld partition of promolecular density
!iIGMtype=2: Based on Hirshfeld partition of actual density (IGMH)
subroutine IGM(iIGMtype)
use functions
use util
use defvar
use GUI
implicit real*8 (a-h,o-z)
character c2000tmp*2000,selectyn
real*8 grad(3),IGM_gradnorm,IGM_gradnorm_inter,gradtmp(3)
integer iIGMtype
integer,allocatable :: IGMfrag(:,:),IGMfragsize(:) !Definition of each fragment used in IGM, and the number of atoms in each fragment
real*8,allocatable :: dg_intra(:,:,:) !delta-g_intra of fragments
real*8,allocatable :: dg_inter(:,:,:) !delta-g_inter between fragment 1 and 2
real*8,allocatable :: dg(:,:,:) !delta-g
real*8,allocatable :: sl2r(:,:,:) !sign(lambda2)rho
real*8,allocatable :: rhogrid(:,:,:) !real density
real*8,allocatable :: gradgrid(:,:,:,:) !real density gradient (1/2/3,i,j,k)
real*8,allocatable :: scatterx(:),scattery(:),tmparr(:),scattery2(:)
integer,allocatable :: tmpidx1(:),tmpidx2(:),allatm(:)
real*8,allocatable :: atmpairdg(:,:) !(i,j) is integral of dg between ith atom in a fragment and jth atom in another fragment
real*8,allocatable :: IBSIWmat(:,:)
!write(*,*) "Citation: Phys. Chem. Chem. Phys., 19, 17928 (2017)"
write(*,*) "Note: Atomic unit is used for all outputs of this function"

write(*,*)
if (iIGMtype==1) then
	write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
    write(*,*) "  Original paper of IGM:"
    write(*,*) "Phys. Chem. Chem. Phys., 19, 17928 (2017)"
    write(*,*) "  Implementation of IGM in Multiwfn:"
    write(*,*) "J. Comput. Chem., 43, 539 (2022) DOI: 10.1002/jcc.26812"
else if (iIGMtype==-1) then !mIGM
	write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
    write(*,*) "J. Comput. Chem., 43, 539 (2022) DOI: 10.1002/jcc.26812"
else if (iIGMtype==2) then
	write(*,*)
	write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
    write(*,*) "  Original paper of IGMH:"
    write(*,*) "Tian Lu, Qinxue Chen, J. Comput. Chem., 43, 539 (2022) DOI: 10.1002/jcc.26812"
    write(*,*) "  Erratum to the IGMH paper:"
    write(*,*) "Tian Lu, Qinxue Chen, ChemRxiv (2022) DOI: 10.26434/chemrxiv-2022-g1m34"
end if
write(*,"(a)") "   Comprehensive reviews:"
write(*,"(a)") " Tian Lu, Qinxue Chen, Visualization Analysis of &
&Weak Interactions in Chemical Systems. In Comprehensive Computational Chemistry, vol. 2, pp. 240-264. Oxford: Elsevier (2024) DOI: 10.1016/B978-0-12-821978-2.00076-3"
write(*,"(a)") " Tian Lu, Visualization Analysis of Covalent and Noncovalent Interactions in Real Space, Angew. Chem. Int. Ed., 137, e202504895 (2025) DOI: 10.1002/anie.202504895"

!----- Define fragments
write(*,*)
write(*,*) "How many fragments will be defined? e.g. 3"
write(*,"(a)") " Note: At least one fragment should be defined. Union set of all fragment atoms will be taken into calculation"
read(*,*) nIGMfrag
if (nIGMfrag<1) then
	write(*,*) "Are you crazy? At least one fragment should be defined!!!"
	return
end if
allocate(IGMfrag(nIGMfrag,ncenter),IGMfragsize(nIGMfrag))
IGMfrag=0
do ifrag=1,nIGMfrag
	write(*,"(a,i3,a)") " Input atom indices for fragment",ifrag,", e.g. 3,5-8,15-20"
	if (nIGMfrag==1) write(*,*) "Note: If input ""a"", the fragment will correspond to the whole system"
	if (nIGMfrag>1.and.ifrag==nIGMfrag) write(*,"(a)") " Note: If input ""c"", the atoms complementary to already defined fragments will be selected"
	read(*,"(a)") c2000tmp
	if (index(c2000tmp,'a')/=0) then
		IGMfragsize(ifrag)=ncenter
		forall (i=1:ncenter) IGMfrag(ifrag,i)=i
	else if (index(c2000tmp,'c')/=0) then
		IGMfragsize(ifrag)=0
		do iatm=1,ncenter
			if (all(IGMfrag(1:ifrag-1,:)/=iatm)) then
				IGMfragsize(ifrag)=IGMfragsize(ifrag)+1
                IGMfrag(ifrag,IGMfragsize(ifrag))=iatm
            end if
        end do
        write(*,"(i5,a,i3)") IGMfragsize(ifrag)," atoms have been selected as fragment",ifrag
	else
		call str2arr(c2000tmp,IGMfragsize(ifrag),IGMfrag(ifrag,:))
	end if
end do

!Set "allatm" array, which contains index of all atoms involved in defined fragments (the coverage of all user-defined fragments may be not equal to the whole system)
allocate(tmpidx1(ncenter))
tmpidx1=0
do ifrag=1,nIGMfrag
	tmpidx1(IGMfrag(ifrag,1:IGMfragsize(ifrag)))=1
end do
ntmp=count(tmpidx1==1)
allocate(allatm(ntmp))
itmp=0
do iatm=1,ncenter
	if (tmpidx1(iatm)==1) then
		itmp=itmp+1
		allatm(itmp)=iatm
	end if
end do
deallocate(tmpidx1)

if (abs(iIGMtype)==1) then
    isl2r=2
    if (allocated(b)) then
	    write(*,"(a)") " Your input file contains wavefunction information, which kind of sign(lambda2)rho would you like to use?"
	    write(*,*) "1 sign(lambda2)rho based on actual electron density"
	    write(*,*) "2 sign(lambda2)rho based on promolecular density"
        write(*,*) "Note: 1 is more accurate and ideal but more expensive"
	    read(*,*) isl2r
    end if
else if (iIGMtype==2) then
    isl2r=1
end if

!----- Set grid
aug3D=2D0 !Smaller than default value
call setgrid(0,igridsel)
allocate(dg_intra(nx,ny,nz),dg_inter(nx,ny,nz),dg(nx,ny,nz),sl2r(nx,ny,nz))
if (iIGMtype==2) allocate(rhogrid(nx,ny,nz),gradgrid(3,nx,ny,nz))

!----- Calculate grid data
call delvirorb(1)
if (ifPBC==0) then
	call gen_GTFuniq(0) !Generate unique GTFs, for faster evaluation in orbderv
else
	call gen_neigh_GTF !Generate neighbouring GTFs list at reduced grids, for faster evaluation
end if
call walltime(iwalltime1)
write(*,*) "Calculating sign(lambda2)rho..."
ifinish=0;ishowprog=1
ntmp=floor(ny*nz/100D0)
!$OMP PARALLEL DO SHARED(ifinish,ishowprog,sl2r,rhogrid,gradgrid) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			if (isl2r==1) then !sign(lambda2)rho based on actual electron density
                if (abs(iIGMtype)==1) then
				    sl2r(i,j,k)=signlambda2rho(tmpx,tmpy,tmpz)
                else
                    !Store density and gradient to "rhogrid" and "gradgrid", which will be passed into IGMgrad_Hirsh, thus avoiding recalculate them later
                    !This treatment is not applied to IGM, because it is quite cheap and often employed for quite large system, this strategy will double memory consuming
                    call signlambda2rho_RDG(tmpx,tmpy,tmpz,sl2r(i,j,k),RDG,rhogrid(i,j,k),gradgrid(:,i,j,k))
                end if
			else !sign(lambda2)rho based on promolecular density
				sl2r(i,j,k)=signlambda2rho_prodens(tmpx,tmpy,tmpz)
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

write(*,*) "Calculating delta-g, delta-g_inter and delta-g_intra..."
ifinish=0
ishowprog=1
dg_inter=0
ntmp=floor(ny*nz/100D0)
!$OMP PARALLEL DO SHARED(ifinish,ishowprog,dg,dg_inter) PRIVATE(i,j,k,tmpx,tmpy,tmpz,grad,gradnorm,IGM_gradnorm,IGM_gradnorm_inter,gradtmp) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			!Calculate gradient vector and IGM gradient norm of whole system and get dg
			if (iIGMtype==1) then !IGM
                call IGMgrad_promol(tmpx,tmpy,tmpz,allatm,grad,IGM_gradnorm)
            else if (iIGMtype==-1) then !mIGM
                call IGMgrad_Hirshpromol(tmpx,tmpy,tmpz,allatm,grad,IGM_gradnorm)
            else if (iIGMtype==2) then !IGMH
                call IGMgrad_Hirsh(tmpx,tmpy,tmpz,allatm,grad,IGM_gradnorm,rhogrid(i,j,k),gradgrid(:,i,j,k))
            end if
            gradnorm=dsqrt(sum(grad**2))
			dg(i,j,k)=IGM_gradnorm-gradnorm
            !Calculate delta-g_inter for defined fragments
            grad=0
			IGM_gradnorm_inter=0
			do ifrag=1,nIGMfrag
                if (iIGMtype==1) then
				    call IGMgrad_promol(tmpx,tmpy,tmpz,IGMfrag(ifrag,1:IGMfragsize(ifrag)),gradtmp(:),IGM_gradnorm) !The IGM_gradnorm given here is useless
                else if (iIGMtype==-1) then
				    call IGMgrad_Hirshpromol(tmpx,tmpy,tmpz,IGMfrag(ifrag,1:IGMfragsize(ifrag)),gradtmp(:),IGM_gradnorm) !The IGM_gradnorm given here is useless
                else if (iIGMtype==2) then
				    call IGMgrad_Hirsh(tmpx,tmpy,tmpz,IGMfrag(ifrag,1:IGMfragsize(ifrag)),gradtmp(:),IGM_gradnorm,rhogrid(i,j,k),gradgrid(:,i,j,k)) !The IGM_gradnorm given here is useless
                end if
                grad=grad+gradtmp
                IGM_gradnorm=dsqrt(sum(gradtmp**2))
				IGM_gradnorm_inter=IGM_gradnorm_inter+IGM_gradnorm
			end do
			dg_inter(i,j,k)=IGM_gradnorm_inter-dsqrt(sum(grad**2))
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

dg_intra=dg-dg_inter

call del_GTFuniq !Destory unique GTF informtaion
call delvirorb_back(1)
call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1

call calc_dvol(dvol)
dginter_int=dvol*sum(dg_inter)
dg_int=dvol*sum(dg)
dgintra_int=dg_int-dginter_int
write(*,*)
write(*,"(a,f12.6,' a.u.')") " Integral of delta-g over whole space:      ",dg_int
write(*,"(a,f12.6,' a.u.')") " Integral of delta-g_inter over whole space:",dginter_int
write(*,"(a,f12.6,' a.u.')") " Integral of delta-g_intra over whole space:",dgintra_int

ymin=0D0
ymax=maxval(dg)
xmin=-0.6D0
xmax=0.2D0

do while (.true.)
	write(*,*)
	write(*,*) "             -------------- Post-processing menu --------------"
	write(*,"(' -3 Change range of Y-axis of scatter graph, current:',f11.6,' to',f11.6)") ymin,ymax
	write(*,"(' -2 Change range of X-axis of scatter graph, current:',f11.6,' to',f11.6)") xmin,xmax
	write(*,*) "-1 Draw scatter graph           1 Save the scatter graph to file"
	write(*,*) "0 Exit"
	write(*,*) "2 Output scatter points to output.txt in current folder"
	write(*,*) "3 Output cube files to current folder"
	write(*,*) "4 Show isosurface of grid data"
	write(*,*) "5 Screen delta-g_intra in high density region"
	write(*,"(a)") " 6 Evaluate contribution of atomic pairs and atoms to interfragment interaction (atomic and atomic pair delta-g indices as well as IBSIW index)"
	write(*,*) "7 Set delta-g where value of sign(lambda2)rho is out of a certain range"
	write(*,*) "8 Set delta-g_inter where value of sign(lambda2)rho is out of a certain range"
	read(*,*) isel
	if (isel==1.or.isel==-1) then
		write(*,*) "1 delta-g_inter vs. sign(lambda2)rho"
		write(*,*) "2 delta-g_intra vs. sign(lambda2)rho"
		write(*,*) "3 delta-g vs. sign(lambda2)rho"
		write(*,*) "4 delta-g_inter + delta-g_intra vs. sign(lambda2)rho"
		read(*,*) itype
		if (.not.allocated(scatterx)) allocate(scatterx(nx*ny*nz),scattery(nx*ny*nz))
		if (itype==4.and.(.not.allocated(scattery2))) allocate(scattery2(nx*ny*nz))
		write(*,*) "Drawing graph, please wait..."
		ii=1
		do k=1,nz
			do j=1,ny
				do i=1,nx
					scatterx(ii)=sl2r(i,j,k)
					if (itype==1) then
						scattery(ii)=dg_inter(i,j,k)
					else if (itype==2) then
						scattery(ii)=dg_intra(i,j,k)
					else if (itype==3) then
						scattery(ii)=dg(i,j,k)
					else
						scattery(ii)=dg_inter(i,j,k)
						scattery2(ii)=dg_intra(i,j,k)
					end if
					ii=ii+1
				end do
			end do
		end do
		if (isel==-1) then
			if (itype==1) call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g^{inter}$ (a.u.)")
			if (itype==2) call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g^{intra}$ (a.u.)")
			if (itype==3) call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g$ (a.u.)")
			if (itype==4) call drawscatter(scatterx,scattery2,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g^{inter/intra}$ (a.u.)",scatterx,scattery,nx*ny*nz)
		else
			isavepic=1
			if (itype==1) call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g^{inter}$ (a.u.)")
			if (itype==2) call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g^{intra}$ (a.u.)")
			if (itype==3) call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g$ (a.u.)")
			if (itype==4) call drawscatter(scatterx,scattery2,nx*ny*nz,xmin,xmax,ymin,ymax,1,"$sign({\lambda}_2)\rho$ (a.u.)","${\delta}g^{inter/intra}$ (a.u.)",scatterx,scattery,nx*ny*nz)
			isavepic=0
			write(*,"(a)") " Figure has been saved to "//trim(graphformat)//" file with ""dislin"" prefix in current directory"
		end if
	else if (isel==-2) then
		write(*,*) "Input lower limit and upper limit of X axis  e.g. 0,1.5"
		read(*,*) xmin,xmax
	else if (isel==-3) then
		write(*,*) "Input lower limit and upper limit of Y axis  e.g. 0,1.5"
		read(*,*) ymin,ymax
	else if (isel==0) then
		exit
	else if (isel==2) then
		open(10,file="output.txt",status="replace")
		write(*,*) "Outputting output.txt..."
		write(10,"(4E16.8)") ((( dg_inter(i,j,k),dg_intra(i,j,k),dg(i,j,k),sl2r(i,j,k),k=1,nz),j=1,ny),i=1,nx)
		close(10)
		write(*,*) "Finished! Column 1/2/3/4 = delta-g_inter/delta-g_intra/delta-g/sign(lambda2)rho"
	else if (isel==3) then
		write(*,*) "Exporting..."
		open(10,file="dg_inter.cub",status="replace")
		call outcube(dg_inter,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "delta-g_inter has been exported to dg_inter.cub in current folder"
		open(10,file="dg_intra.cub",status="replace")
		call outcube(dg_intra,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "delta-g_intra has been exported to dg_intra.cub in current folder"
		open(10,file="dg.cub",status="replace")
		call outcube(dg,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "delta-g has been exported to dg.cub in current folder"
		open(10,file="sl2r.cub",status="replace")
		call outcube(sl2r,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
		write(*,*) "sign(lambda2)rho has been exported to sl2r.cub in current folder"
	else if (isel==4) then
		write(*,*) "1 delta-g_inter"
		write(*,*) "2 delta-g_intra"
		write(*,*) "3 delta-g"
		write(*,*) "4 sign(lambda2)rho"
		read(*,*) itype
		if (allocated(cubmat)) deallocate(cubmat)
		allocate(cubmat(nx,ny,nz))
		if (itype==1) cubmat=dg_inter
		if (itype==2) cubmat=dg_intra
		if (itype==3) cubmat=dg
		if (itype==4) cubmat=sl2r
	 	if (itype==1) then
	 		write(*,*) "Input the value of isosurface, e.g. 0.02"
	 	else
	 		write(*,*) "Input the value of isosurface, e.g. 0.5"
	 	end if
		read(*,*) sur_value
		call drawisosurgui(1)
	else if (isel==5) then
		write(*,*) "Input range of sign(lambda2)rho, e.g. -0.12,0.08"
		write(*,"(a)") " delta-g_intra will be set to 0 if sign(lambda2)rho is out of this range"
		read(*,*) denslow,denshigh
		do k=1,nz
			do j=1,ny
				do i=1,nx
					if (sl2r(i,j,k)<denslow.or.sl2r(i,j,k)>denshigh) dg_intra(i,j,k)=0
				end do
			end do
		end do
		write(*,*) "Done!"
        
	else if (isel==6) then !Calculate and print atomic and atomic pair dg indices as well as IBSIW
		if (nIGMfrag==2) then
			ifrag=1
			jfrag=2
		else
			write(*,*) "Input fragment index to select two fragments, e.g. 1,3"
			read(*,*) ifrag,jfrag
		end if
		ni=IGMfragsize(ifrag)
		nj=IGMfragsize(jfrag)
        
        !Calculate atomic pair delta-g indices
		allocate(atmpairdg(ni,nj))
        call calcatmpairdg(iIGMtype,ni,IGMfrag(ifrag,:),nj,IGMfrag(jfrag,:),atmpairdg)
        
		open(10,file="atmdg.txt",status="replace")
		!Output the first fragment
		allocate(tmparr(ni),tmpidx1(ni))
		do iatmtmp=1,ni
			tmparr(iatmtmp)=sum(atmpairdg(iatmtmp,:))
			tmpidx1(iatmtmp)=IGMfrag(ifrag,iatmtmp)
		end do
		call sort(tmparr,list=tmpidx1) !Sort from small to large
		write(10,"(' Atomic delta-g indices of fragment',i3,' and percentage contributions')") ifrag
        totval=sum(tmparr(:))
		do idx=ni,1,-1
			write(10,"(' Atom',i5,' :',f12.6,'  (',f7.2,' % )')") tmpidx1(idx),tmparr(idx),tmparr(idx)/totval*100
		end do
		deallocate(tmparr,tmpidx1)
		!Output the second fragment
		allocate(tmparr(nj),tmpidx1(nj))
		do jatmtmp=1,nj
			tmparr(jatmtmp)=sum(atmpairdg(:,jatmtmp))
			tmpidx1(jatmtmp)=IGMfrag(jfrag,jatmtmp)
		end do
		call sort(tmparr,list=tmpidx1) !Sort from small to large
		write(10,"(/,' Atomic delta-g indices of fragment',i3,' and percentage contributions')") jfrag
		do jdx=nj,1,-1
			write(10,"(' Atom',i5,' :',f12.6,'  (',f7.2,' % )')") tmpidx1(jdx),tmparr(jdx),tmparr(jdx)/totval*100
		end do
		deallocate(tmparr,tmpidx1)
		!Output atomic pair dg indices
		allocate(tmparr(ni*nj),tmpidx1(ni*nj),tmpidx2(ni*nj))
		itmp=0
		do iatmtmp=1,ni
			do jatmtmp=1,nj
				itmp=itmp+1
				tmpidx1(itmp)=IGMfrag(ifrag,iatmtmp)
				tmpidx2(itmp)=IGMfrag(jfrag,jatmtmp)
				tmparr(itmp)=atmpairdg(iatmtmp,jatmtmp)
			end do
		end do
		call sort(tmparr,list=tmpidx1,list2=tmpidx2) !Sort from small to large
		write(10,"(/,' Atomic pair delta-g indices and percentage contributions (zero terms are not shown)')")
		do idx=ni*nj,1,-1
			if (tmparr(idx)==0) cycle
			write(10,"(2i5,' :',f12.6,'  (',f7.2,' % )')") tmpidx1(idx),tmpidx2(idx),tmparr(idx),tmparr(idx)/totval*100
		end do
        write(10,"(/,a,f12.6)") " Sum of all atomic pair delta-g indices:",sum(atmpairdg(:,:))
		close(10)
        deallocate(tmparr,tmpidx1,tmpidx2)
		write(*,"(a)") " Atomic and atomic pair delta-g indices have been outputted to atmdg.txt in current folder"
        
        !Calculate IBSIW index
        allocate(IBSIWmat(ni,nj))
		open(10,file="IBSIW.txt",status="replace")
        !Convert atomic pair dg indices to IBSIW
		do idx=1,ni
			do jdx=1,nj
                dist=atomdist(IGMfrag(ifrag,idx),IGMfrag(jfrag,jdx),1)*b2a
                IBSIWmat(idx,jdx)=atmpairdg(idx,jdx)/dist**2*100
			end do
		end do
		!Output the first fragment
		allocate(tmparr(ni),tmpidx1(ni))
		do iatmtmp=1,ni
			tmparr(iatmtmp)=sum(IBSIWmat(iatmtmp,:))
			tmpidx1(iatmtmp)=IGMfrag(ifrag,iatmtmp)
		end do
		call sort(tmparr,list=tmpidx1) !Sort from small to large
		write(10,"(' Sum of related IBSIW for each atom in fragment',i3)") ifrag
		do idx=ni,1,-1
			write(10,"(' Atom',i5,' :',f12.6)") tmpidx1(idx),tmparr(idx)
		end do
		deallocate(tmparr,tmpidx1)
		!Output the second fragment
		allocate(tmparr(nj),tmpidx1(nj))
		do jatmtmp=1,nj
			tmparr(jatmtmp)=sum(IBSIWmat(:,jatmtmp))
			tmpidx1(jatmtmp)=IGMfrag(jfrag,jatmtmp)
		end do
		call sort(tmparr,list=tmpidx1) !sort from small to large
		write(10,"(/,' Sum of related IBSIW for each atom in fragment',i3)") jfrag
		do jdx=nj,1,-1
			write(10,"(' Atom',i5,' :',f12.6)") tmpidx1(jdx),tmparr(jdx)
		end do
		deallocate(tmparr,tmpidx1)
		!Output atomic pair dg indices
		allocate(tmparr(ni*nj),tmpidx1(ni*nj),tmpidx2(ni*nj))
		itmp=0
		do iatmtmp=1,ni
			do jatmtmp=1,nj
				itmp=itmp+1
				tmpidx1(itmp)=IGMfrag(ifrag,iatmtmp)
				tmpidx2(itmp)=IGMfrag(jfrag,jatmtmp)
				tmparr(itmp)=IBSIWmat(iatmtmp,jatmtmp)
			end do
		end do
		call sort(tmparr,list=tmpidx1,list2=tmpidx2) !Sort from small to large
		write(10,"(/,' IBSIW index (zero terms are not shown)')")
		do idx=ni*nj,1,-1
			if (tmparr(idx)==0) cycle
			write(10,"(2i5,' :',f12.6)") tmpidx1(idx),tmpidx2(idx),tmparr(idx)
		end do
		close(10)
        deallocate(IBSIWmat,tmparr,tmpidx1,tmpidx2)
		write(*,*) "IBSIW have been outputted to IBSIW.txt in current folder"
        
        iouttype=1 !1: atmdg.pdb. 2: atmdg.pqr & atmdg%.pqr. pqr is not good choice because VMD is unable to load elements from it, making determination of bonding often incorrect
        if (iouttype==1) then
            write(*,"(/,a)") " If outputting the two fragments as atmdg.pdb in current folder, whose ""Beta"" and ""Occupancy"" &
            &atomic properties correspond to atomic delta-g index multiplied by 10 and percentage atomic delta-g index, respectively? (y/n)"
			read(*,*) selectyn
			if (selectyn=='y') then
				open(10,file="atmdg.pdb",status="replace")
                write(10,"('REMARK   Generated by Multiwfn, Totally',i10,' atoms')") ni+nj
				i=1
				do iatmtmp=1,ni
					val=sum(atmpairdg(iatmtmp,:))*10 !Atomic delta-g multiplied by 10
					val2=val/totval*10 !Percentage contribution
					iatm=IGMfrag(ifrag,iatmtmp)
					write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
					"HETATM",i,' '//ind2name_up(a(iatm)%index)//' ',"MOL",'A',1,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a,val2,val,adjustr(ind2name_up(a(iatm)%index))
					i=i+1
				end do
				do jatmtmp=1,nj
					val=sum(atmpairdg(:,jatmtmp))*10
					val2=val/totval*10
					jatm=IGMfrag(jfrag,jatmtmp)
					write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,2f6.2,10x,a2)") &
					"HETATM",i,' '//ind2name_up(a(jatm)%index)//' ',"MOL",'B',1,a(jatm)%x*b2a,a(jatm)%y*b2a,a(jatm)%z*b2a,val2,val,adjustr(ind2name_up(a(jatm)%index))
					i=i+1
				end do
				write(10,"('END')")
				close(10)
                write(*,*) "Done! atmdg.pdb has been outputted to current folder"
			end if
        else if (iouttype==2) then !atmdg.pqr & atmdg%.pqr
   !     		write(*,"(/,a)") " If outputting the two fragments as atmdg.pqr and atmdg%.pqr in current folder, whose ""Charge"" atomic property carries &
			!atom delta-g index and its percentage contribution, respectively? (y/n)"
			!read(*,*) selectyn
			!if (selectyn=='y') then
			!	do itime=1,2
			!		if (itime==1) open(10,file="atmdg.pqr",status="replace")
			!		if (itime==2) open(10,file="atmdg%.pqr",status="replace")
			!		write(10,"('REMARK   Generated by Multiwfn, Totally',i10,' atoms')") ni+nj
			!		i=1
			!		do iatmtmp=1,ni
			!			val=sum(atmpairdg(iatmtmp,:))
			!			if (itime==2) val=val/totval
			!			iatm=IGMfrag(ifrag,iatmtmp)
			!			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f8.3,1x,a2)") &
			!			"HETATM",i,' '//ind2name_up(a(iatm)%index)//' ',"MOL",'A',1,a(iatm)%x*b2a,a(iatm)%y*b2a,a(iatm)%z*b2a,val,vdwr(a(iatm)%index)*b2a,adjustr(ind2name_up(a(iatm)%index))
			!			i=i+1
			!		end do
			!		do jatmtmp=1,nj
			!			val=sum(atmpairdg(:,jatmtmp))
			!			if (itime==2) val=val/totval
			!			jatm=IGMfrag(jfrag,jatmtmp)
			!			write(10,"(a6,i5,1x,a4,1x,a3, 1x,a1,i4,4x,3f8.3,f13.8,f8.3,1x,a2)") &
			!			"HETATM",i,' '//ind2name_up(a(jatm)%index)//' ',"MOL",'B',1,a(jatm)%x*b2a,a(jatm)%y*b2a,a(jatm)%z*b2a,val,vdwr(a(iatm)%index)*b2a,adjustr(ind2name_up(a(jatm)%index))
			!			i=i+1
			!		end do
			!		write(10,"('END')")
			!		close(10)
			!	end do
			!	write(*,*) "Done! atmdg.pqr and atmdg%.pqr have been outputted to current folder"
			!end if
        end if
        
		deallocate(atmpairdg)
		
	else if (isel==7.or.isel==8) then
		write(*,"(a)") " Input lower and upper limit of range of sign(lambda2)rho, e.g. -0.04,-0.025"
		read(*,*) rlower,rupper
		if (rupper<rlower) then
			write(*,*) "Error: Upper limit is smaller than lower limit!"
			write(*,*) "Press ENTER button to continue. The data will not be changed"
			read(*,*)
			cycle
		end if
		if (isel==7) write(*,*) "Input expected value of delta-g, e.g. 0"
		if (isel==8) write(*,*) "Input expected value of delta-g_inter, e.g. 0"
		read(*,*) tmpval
		if (isel==7) where (sl2r>rupper.or.sl2r<rlower) dg=tmpval
		if (isel==8) where (sl2r>rupper.or.sl2r<rlower) dg_inter=tmpval
		write(*,*) "Done!"
	end if !end of menu
end do
end subroutine




!!------ Calculate atomic pair delta-g index of atomic pairs between two lists atmlist1 and atmlist2
! iIGMtype=1: IGM; =2: IGM based on Hirshfeld partition (IGMH); iIGMtype=-1: mIGM
subroutine calcatmpairdg(iIGMtype,natm1,atmlist1,natm2,atmlist2,atmpairdg)
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
integer iIGMtype,natm1,natm2,atmlist1(natm1),atmlist2(natm2),atmpair(2)
real*8 atmpairdg(natm1,natm2),atmpairdg_add(natm1,natm2),atmpairdg_addtmp(natm1,natm2),grad(3),IGM_grad(3)
type(content),allocatable :: gridatmorg(:),gridatm(:)
real*8,allocatable :: beckeweigrid(:)
real*8 atmrho(ncenter),atmgrad(3,ncenter),gradtmp(3),gradnormIGMtmp
real*8 atmprorho(ncenter),atmprograd(3,ncenter)
real*8 prorho,prograd(3),realrho,realgrad(3),hess(3,3),Hirshwei(ncenter)

if (ifPBC>0) then
    write(*,*)
	write(*,"(a)") " Warning: This function currently does not fully support periodic systems. &
    &To use it for a periodic system, the atoms of interest should be far from box boundary!"
    write(*,*)
end if

iapprox=1 !Enable use approximation to accelerate calculation

nradpot_bk=radpot
nsphpot_bk=sphpot
write(*,*) "Please select grid for integrating the delta-g_pair functions"
!IGMH and mIGM have higher requirement on integration grid since distribution region is narrow
if (iIGMtype==2) write(*,*) "Note: Option 1 is deprecated since numerical accuracy is too low for IGMH"
if (iIGMtype==-1) write(*,*) "Note: Option 1 is deprecated since numerical accuracy is too low for mIGM"
write(*,*) "1 Medium quality (radial=30, angular=110. Cost=1.0 x)"
write(*,*) "2 High quality (radial=40, angular=170. Cost=2.1 x)"
write(*,*) "3 Ultrafine quality (radial=60, angular=302. Cost=5.5 x)"
write(*,*) "4 Perfect quality (radial=75, angular=434. Cost=9.9 x)"
write(*,*) "5 Use ""radpot"" and ""sphpot"" defined in settings.ini"
read(*,*) isel
if (isel==1) then !In fact this is already enough to guarantee quantitative accuracy
    radpot=30
    sphpot=110
else if (isel==2) then
    radpot=40
    sphpot=170
else if (isel==3) then
    radpot=60
    sphpot=302
else if (isel==4) then
    radpot=75
    sphpot=434
end if
allocate(beckeweigrid(radpot*sphpot),gridatm(radpot*sphpot),gridatmorg(radpot*sphpot))
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call walltime(iwalltime1)
call gen1cintgrid(gridatmorg,iradcut)
call showprog(0,ncenter)
    
!Becke's multi-center integration, cycle all atoms
!For each atom grids, calculate atmpairdg_add(:,:), which records all dg_pair contributed by these grids, and added to atmpairdg
atmpairdg=0
do icen=1,ncenter
    !If this center is far from any pair of closely contacted atoms among the two fragments, &
	!then grids centered on this atom should not contribute to any dg_pair, so this center will be skipped
    if (iapprox==1) then
        distmin=1E10
        do itmp=1,natm1
            iatm=atmlist1(itmp)
            do jtmp=1,natm2
                jatm=atmlist2(jtmp)
                if (atomdist(iatm,jatm,1)>9) cycle !The two atoms are closely contacted
                xmid=(a(iatm)%x+a(jatm)%x)/2
                ymid=(a(iatm)%y+a(jatm)%y)/2
                zmid=(a(iatm)%z+a(jatm)%z)/2
                dist=dsqrt((xmid-a(icen)%x)**2+(ymid-a(icen)%y)**2+(zmid-a(icen)%z)**2)
                if (dist<distmin) distmin=dist
            end do
        end do
        if (distmin>6) cycle !This is found to be lowest acceptable threshold. The accuracy loss in this case is fully negligible
    end if
        
	gridatm%x=gridatmorg%x+a(icen)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(icen)%y
	gridatm%z=gridatmorg%z+a(icen)%z
	call gen1cbeckewei(icen,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
    atmpairdg_add=0
	!$OMP parallel shared(atmpairdg_add) private(atmpairdg_addtmp,ipt,rnowx,rnowy,rnowz,atmrho,atmgrad,itmp,jtmp,iatm,jatm,&
    !$OMP gradtmp,gradnormIGMtmp,atmprorho,atmprograd,realrho,realgrad,prorho,prograd,Hirshwei,idir,t1,t2,t3,tmp,dgtmp) num_threads(nthreads)
	atmpairdg_addtmp=0
    !$OMP DO schedule(dynamic)
	do ipt=1+iradcut*sphpot,radpot*sphpot
		rnowx=gridatm(ipt)%x
        rnowy=gridatm(ipt)%y
        rnowz=gridatm(ipt)%z
        
        !Calculate gradient of every isolated atom (atmgrad)
        if (iIGMtype==1) then !IGM
            !Calculate atomic gradient in this grid if it is involved in either atmlist1 or atmlist2
            do iatm=1,ncenter
                if (any(atmlist1==iatm).or.any(atmlist2==iatm)) then
                    call proatmgrad(1,iatm,rnowx,rnowy,rnowz,atmrho(iatm),atmgrad(:,iatm))
                end if
            end do
        else if (iIGMtype==2) then !IGMH
            !Calculate molecular density and its gradient. This is the major overhead (>80%)
            call calchessmat_dens(1,rnowx,rnowy,rnowz,realrho,realgrad,hess)
            !Calculate atom density and gradient in free state
            do iatm=1,ncenter
                if (iapprox==1) then !Ignore atom farther than 8 Bohr from current point, this is quite safe
                    tmp=(rnowx-a(iatm)%x)**2+(rnowy-a(iatm)%y)**2+(rnowz-a(iatm)%z)**2
                    if (tmp>64) then
                        atmprorho(iatm)=0
                        atmprograd(:,iatm)=0
                        cycle
                    end if
                end if
                call proatmgrad(1,iatm,rnowx,rnowy,rnowz,atmprorho(iatm),atmprograd(:,iatm))
            end do
            !Calculate promolecular density and gradient of molecule
            prorho=sum(atmprorho(:))
            do idir=1,3
                prograd(idir)=sum(atmprograd(idir,:))
            end do
            !Calculate Hirshfeld weight
            do iatm=1,ncenter
                if (prorho==0) then
                    Hirshwei(iatm)=0
                else
                    Hirshwei(iatm)=atmprorho(iatm)/prorho
                end if
            end do
            !Calculate gradient of atomic density partitioned by Hirshfeld
            do iatm=1,ncenter
                if (any(atmlist1==iatm).or.any(atmlist2==iatm)) then
                    do idir=1,3
                        t1=Hirshwei(iatm)*realgrad(idir)
                        if (prorho==0) then
                            t2=0;t3=0
                        else
                            t2=realrho/prorho*atmprograd(idir,iatm)
                            t3=-realrho*atmprorho(iatm)/prorho**2 * prograd(idir)
                        end if
                        !atmgrad(idir,iatm)=t1+t2+t3 !Mathematically correct free-state atomic gradient but poor effect for IGMH. Older than 2022-Sep-18
                        atmgrad(idir,iatm)=t1-t2-t3 !IGMH-type "special" free-state atomic gradient
                    end do
                end if
            end do
        else if (iIGMtype==-1) then !mIGM
            !Calculate atom density and gradient in free state
            do iatm=1,ncenter
                if (iapprox==1) then !Ignore atom farther than 8 Bohr from current point, this is quite safe
                    tmp=(rnowx-a(iatm)%x)**2+(rnowy-a(iatm)%y)**2+(rnowz-a(iatm)%z)**2
                    if (tmp>64) then
                        atmprorho(iatm)=0
                        atmprograd(:,iatm)=0
                        cycle
                    end if
                end if
                call proatmgrad(2,iatm,rnowx,rnowy,rnowz,atmprorho(iatm),atmprograd(:,iatm))
            end do
            !Calculate promolecular density and gradient of molecule
            prorho=sum(atmprorho(:))
            do idir=1,3
                prograd(idir)=sum(atmprograd(idir,:))
            end do
			!Calculate gradient of atomic density partitioned by Hirshfeld
            do iatm=1,ncenter
                if (any(atmlist1==iatm).or.any(atmlist2==iatm)) then
					if (prorho==0) then
						atmgrad(:,iatm)= -atmprograd(:,iatm)
                    else
						atmgrad(:,iatm)=2*atmprorho(iatm)/prorho*prograd(:) - atmprograd(:,iatm)
                    end if
                end if
			end do
        end if
        
        !Calculate atomic pair delta-g matrix contributed by this point
        do itmp=1,natm1
            iatm=atmlist1(itmp)
            do jtmp=1,natm2
                jatm=atmlist2(jtmp)
                gradtmp(:)=atmgrad(:,iatm)+atmgrad(:,jatm)
                gradnormIGMtmp=dsqrt(sum(atmgrad(:,iatm)**2))+dsqrt(sum(atmgrad(:,jatm)**2))
                dgtmp=gradnormIGMtmp-dsqrt(sum(gradtmp**2))
                atmpairdg_addtmp(itmp,jtmp)=atmpairdg_addtmp(itmp,jtmp)+dgtmp*gridatmorg(ipt)%value*beckeweigrid(ipt)
            end do
        end do
	end do
	!do ipt=1+iradcut*sphpot,radpot*sphpot
	!	atmpairdg(:,:)=atmpairdg(:,:)+dgmat(:,:,ipt)*gridatmorg(ipt)%value*beckeweigrid(ipt)
	!end do
	!$OMP END DO
	!$OMP CRITICAL
	atmpairdg_add=atmpairdg_add+atmpairdg_addtmp
	!$OMP END CRITICAL
	!$OMP END PARALLEL
    atmpairdg=atmpairdg+atmpairdg_add
    call showprog(icen,ncenter)
end do
call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
radpot=nradpot_bk
sphpot=nsphpot_bk
end subroutine




!!-------- Averaged IGM/mIGM. Supports PBC
!iIGMtype=1: aIGM   =-1: amIGM
!iIGMtype=3 (obsolete): Calculate averaged (regular/IGM) density gradient first, then calculate aIGM
subroutine aIGM(iIGMtype)
use defvar
use util
use GUI
use functions
implicit real*8 (a-h,o-z)
real*8 gradtmp(3),grad_inter(3),IGM_gradnorm_inter,vec1(3),vec2(3)
integer iIGMtype
integer,allocatable :: IGMfrag(:,:),IGMfragsize(:) !Definition of each fragment used in IGM, and the number of atoms in each fragment
real*8,allocatable :: frag_grad(:,:,:,:,:) !frag_grad(1:3,nx,ny,nz,nfrag), gradient vector of each fragment at every point
real*8,allocatable :: dg_inter(:,:,:),TFI_IGM(:,:,:)
logical,allocatable :: dogrid(:,:,:)
!The first index of avggrad and the first two indices of avghess correspond to components of gradient and Hessian, respectively
real*8,allocatable :: avgdens(:,:,:),avggrad(:,:,:,:),avghess(:,:,:,:,:)
real*8,allocatable :: avgRDG(:,:,:),thermflu(:,:,:),avgsl2r(:,:,:)
real*8,allocatable :: scatterx(:),scattery(:)
character c2000tmp*2000,selectyn
real*8 prorho,prograd(3),atmprorho(ncenter),atmprograd(3,ncenter),atmgrad(3),atomcoeff(10),atomexp(10)

write(*,*) "*** Please cite the following papers along with Multiwfn original papers ***"
write(*,"(a)") "   Original paper of aIGM: Tian Lu, Qinxue Chen, Visualization Analysis of &
&Weak Interactions in Chemical Systems. In Comprehensive Computational Chemistry, vol. 2, pp. 240-264. Oxford: Elsevier (2024) DOI: 10.1016/B978-0-12-821978-2.00076-3"

if (iIGMtype==1) then
	write(*,*)
	write(*,"(a)") " Warning: amIGM is the significantly better choice than aIGM! Please consider to use amIGM instead"
end if

write(*,*)
write(*,*) "How many fragments will be defined? e.g. 2"
write(*,"(a)") " Note: At least two fragments should be defined"
read(*,*) nIGMfrag
allocate(IGMfrag(nIGMfrag,ncenter),IGMfragsize(nIGMfrag))
do ifrag=1,nIGMfrag
	write(*,"(a,i3,a)") " Input atom indices for fragment",ifrag,", e.g. 3,5-8,15-20"
	if (ifrag==nIGMfrag) write(*,"(a)") " Note: If input ""c"", the atoms complementary to already defined fragments will be selected"
	read(*,"(a)") c2000tmp
	if (index(c2000tmp,'c')/=0) then
		IGMfragsize(ifrag)=0
		do iatm=1,ncenter
			if (all(IGMfrag(1:ifrag-1,:)/=iatm)) then
				IGMfragsize(ifrag)=IGMfragsize(ifrag)+1
                IGMfrag(ifrag,IGMfragsize(ifrag))=iatm
            end if
        end do
        write(*,"(i5,a,i3)") IGMfragsize(ifrag)," atoms have been selected as fragment",ifrag
	else
		call str2arr(c2000tmp,IGMfragsize(ifrag),IGMfrag(ifrag,:))
	end if
end do

write(*,*)
write(*,*) "Input range of the frames to be analyzed, e.g. 150,400 means from 150 to 400 frames"
write(*,*) "Note: The frame index starts from 1"
read(*,*) ifpsstart,ifpsend
nfps=ifpsend-ifpsstart+1
write(*,"(' Selected',i8,' frames, frames from',i8,' to',i8,' will be processed',/)") nfps,ifpsstart,ifpsend

call setgrid(0,igridsel)

call walltime(iwalltime1)

write(*,"(/,a)") " Calculating averaged density, gradient and Hessian of density..."
allocate(avgdens(nx,ny,nz),avggrad(3,nx,ny,nz),avghess(3,3,nx,ny,nz))
call avg_rhogradhess(avgdens,avggrad,avghess,ifpsstart,ifpsend)

write(*,"(a)") " Calculating averaged RDG and averaged sign(lambda2)*rho..."
allocate(avgRDG(nx,ny,nz),avgsl2r(nx,ny,nz))
call avg_RDG_sl2r(avgdens,avggrad,avghess,avgRDG,avgsl2r) !RDG is a byproduct
deallocate(avggrad,avghess) !Will not be used further, so release its memory
call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time until now',i10,' s')") iwalltime2-iwalltime1

if (abs(iIGMtype)==1) then !Calculate IGM/mIGM for each frame, then take average
    allocate(dg_inter(nx,ny,nz))
    open(10,file=filename,status="old")
    if (iIGMtype==1) write(*,*) "Calculating grid data of averaged IGM..."
    if (iIGMtype==-1) write(*,*) "Calculating grid data of averaged mIGM..."
    dg_inter=0
    
    !Prescreening grids. If a grid is within scaled vdW radius of any atom of fragment 1, then this grid will be calculated. 
    !This treatment can safely reduce certain computational cost if scale factor is set to 2
    if (amIGMvdwscl/=0) then
		write(*,"(a,f6.2)") " Prescreening grids with amIGMvdwscl parameter:",amIGMvdwscl
		allocate(dogrid(nx,ny,nz))
		dogrid=.false.
		do k=1,nz
			do j=1,ny
				do i=1,nx
					call getgridxyz(i,j,k,vec1(1),vec1(2),vec1(3))
					do idx=1,IGMfragsize(1)
						iatm=IGMfrag(1,idx)
						vec2(1)=a(iatm)%x
						vec2(2)=a(iatm)%y
						vec2(3)=a(iatm)%z
						call nearest_dist(vec1,vec2,dist)
						if (dist < amIGMvdwscl*vdwr(a(iatm)%index)) then
							dogrid(i,j,k)=.true.
							exit
						end if
					end do
				end do
			end do
		end do
        write(*,"(' Percent of screened grids:',f10.2,'%')") dfloat(count(dogrid.eqv..false.))/(nx*ny*nz)*100
    end if
    
    !Loop frames
    do ifps=1,ifpsend
	    call readxyztrj(10)
	    if (ifps<ifpsstart) cycle
        call showprog(ifps,nfps)
        !$OMP PARALLEL DO SHARED(dg_inter) PRIVATE(i,j,k,ifrag,gradtmp,grad_inter,IGM_gradnorm_inter,tmpx,tmpy,tmpz,tmpval,&
        !$OMP idx,iatm,atmprorho,atmprograd,prorho,prograd,idir,&
        !$OMP rx,ry,rz,rx2,ry2,rz2,r,r2,iele,nSTO,atomcoeff,atomexp,iSTO,tmp,term) schedule(dynamic) NUM_THREADS(nthreads)
        do k=1,nz
	        do j=1,ny
		        do i=1,nx
					if (amIGMvdwscl/=0) then
						if (dogrid(i,j,k).eqv..false.) cycle
                    end if
					call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                    grad_inter=0
                    IGM_gradnorm_inter=0
                    if (iIGMtype==1) then !aIGM
						do ifrag=1,nIGMfrag
							call IGMgrad_promol(tmpx,tmpy,tmpz,IGMfrag(ifrag,1:IGMfragsize(ifrag)),gradtmp(:),tmpval)
							grad_inter(:)=grad_inter(:)+gradtmp(:)
							IGM_gradnorm_inter=IGM_gradnorm_inter+dsqrt(sum(gradtmp**2))
						end do
                    else !amIGM
						!  In fact, one can simply replacing "IGMgrad_promol" in above code by "IGMgrad_Hirshpromol" to realize amIGM
						!However, the cost is high, because proatomic density for all atoms will be evaluated in this subroutine when looping each fragment
						!So, in the following code, promolecular is calculated first, and record proatomic gradient at the same time
                    
						!Obtain density and gradient of all atoms in their isolated states
                        atmprorho=0
                        atmprograd=0
						do iatm=1,ncenter
							!call proatmgrad(2,iatm,tmpx,tmpy,tmpz,atmprorho(iatm),atmprograd(:,iatm))
							!The following code is equivalent to the above line, while more efficient
							rx=tmpx-a(iatm)%x
							ry=tmpy-a(iatm)%y
							rz=tmpz-a(iatm)%z
							rx2=rx*rx
							ry2=ry*ry
							rz2=rz*rz
							r2=rx2+ry2+rz2
							iele=a(iatm)%index
							if (r2<atmrhocutsqr_1En5(iele)) then
								r=dsqrt(r2)
								call genatmraddens_STOfitparm(iele,nSTO,atomcoeff,atomexp)
								do iSTO=1,nSTO
									term=atomcoeff(iSTO)*dexp(-r*atomexp(iSTO))
									atmprorho(iatm)=atmprorho(iatm)+term
									if (r/=0) then
										tmp=term*atomexp(iSTO)/r
										atmprograd(1,iatm)=atmprograd(1,iatm)-tmp*rx
										atmprograd(2,iatm)=atmprograd(2,iatm)-tmp*ry
										atmprograd(3,iatm)=atmprograd(3,iatm)-tmp*rz
									end if
								end do
							end if
						end do
                        
						prorho=sum(atmprorho(:))
						!Cycle fragments
                        if (prorho/=0) then
							!Calculate gradient of promolecular density
							do idir=1,3
								prograd(idir)=sum(atmprograd(idir,:))
							end do
							do ifrag=1,nIGMfrag
								!Calculate gradient and IGM gradient of this fragment in the mIGM way
                                gradtmp=0
								do idx=1,IGMfragsize(ifrag)
									iatm=IGMfrag(ifrag,idx)
									atmgrad(:)=2*atmprorho(iatm)/prorho*prograd(:) - atmprograd(:,iatm)
									gradtmp=gradtmp+atmgrad
								end do
								grad_inter(:)=grad_inter(:)+gradtmp(:)
								IGM_gradnorm_inter=IGM_gradnorm_inter+dsqrt(sum(gradtmp**2))
							end do
                        end if
                    end if
                    
                    dg_inter(i,j,k)=dg_inter(i,j,k) + IGM_gradnorm_inter - dsqrt(sum(grad_inter**2))
		        end do
	        end do
        end do
        !$OMP END PARALLEL DO
    end do
    close(10)
    dg_inter=dg_inter/nfps
    
else if (iIGMtype==3) then !Obsolete, result is poor. Calculate averaged (regular/IGM) density gradient first, then calculate aIGM
    allocate(dg_inter(nx,ny,nz),frag_grad(3,nx,ny,nz,nIGMfrag))
    open(10,file=filename,status="old")
    write(*,*) "Calculating grid data of averaged density gradient of each fragment..."
    frag_grad=0
    do ifps=1,ifpsend
	    call readxyztrj(10)
	    if (ifps<ifpsstart) cycle
        call showprog(ifps,nfps)
        do ifrag=1,nIGMfrag
            !$OMP PARALLEL DO SHARED(frag_grad) PRIVATE(i,j,k,tmpx,tmpy,tmpz,gradtmp) schedule(dynamic) NUM_THREADS(nthreads)
            do k=1,nz
	            do j=1,ny
		            do i=1,nx
						call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                        !gradtmp is gradient vector of ifrag
                        call IGMgrad_promol(tmpx,tmpy,tmpz,IGMfrag(ifrag,1:IGMfragsize(ifrag)),gradtmp(:),rnouse)
                        frag_grad(:,i,j,k,ifrag)=frag_grad(:,i,j,k,ifrag)+gradtmp(:)
		            end do
	            end do
            end do
            !$OMP END PARALLEL DO
        end do
    end do
    close(10)
    frag_grad=frag_grad/nfps

    write(*,*) "Calculating grid data of delta-g_inter..."
    do k=1,nz
	    do j=1,ny
		    do i=1,nx
                grad_inter=0
                IGM_gradnorm_inter=0
                do ifrag=1,nIGMfrag
                    grad_inter=grad_inter+frag_grad(:,i,j,k,ifrag)
	                IGM_gradnorm_inter=IGM_gradnorm_inter+dsqrt(sum(frag_grad(:,i,j,k,ifrag)**2))
                end do
                dg_inter(i,j,k)=IGM_gradnorm_inter-dsqrt(sum(grad_inter**2))
	        end do
        end do
    end do    
end if

call walltime(iwalltime2)
write(*,"(' Calculation totally took up wall clock time',i10,' s')") iwalltime2-iwalltime1

ymin=0D0
ymax=maxval(dg_inter)
xmin=-0.5D0
xmax=0.5D0

do while (.true.)
	write(*,*)
	write(*,"(' -3 Change range of Y-axis of scatter graph, current:',f11.6,' to',f11.6)") ymin,ymax
	write(*,"(' -2 Change range of X-axis of scatter graph, current:',f11.6,' to',f11.6)") xmin,xmax
	write(*,*) "-1 Draw scatter map between averaged delta-g_inter and sign(lambda2)*rho"
	write(*,*) "0 Exit"
	write(*,*) "1 Save the scatter map to image file"
	write(*,*) "2 Output data of the scatter map to output.txt"
	write(*,"(a)") " 3 Output averaged delta-g_inter and sign(lambda2)*rho to avgdg_inter.cub and avgsl2r.cub in current folder, respectively"
	write(*,"(a)") " 4 Output averaged RDG to avgRDG.cub in current folder"
	write(*,"(a)") " 5 Compute thermal fluctuation index (TFI) and export to thermflu.cub in current folder"
    !I found the effect of mapping TFI onto dg_inter isosurface is poor (two sides of isosurface show very different color), so hidden these options
	!if (iIGMtype==1) write(*,*) "6 Compute TFI(aIGM) and export to TFI_aIGM.cub in current folder"
	!if (iIGMtype==-1) write(*,*) "6 Compute TFI(amIGM) and export to TFI_amIGM.cub in current folder"
	!write(*,"(a)") " 7 Evaluate contribution of atomic pairs and atoms to interfragment interaction (atom and atomic pair delta-g indices as well as IBSIW index)"
	read(*,*) isel
    
	if (isel==-3) then
		write(*,*) "Input lower limit and upper limit of Y axis  e.g. 0,1.5"
		read(*,*) ymin,ymax
    else if (isel==-2) then
		write(*,*) "Input lower limit and upper limit of X axis  e.g. 0,1.5"
		read(*,*) xmin,xmax
        
    else if (isel==0) then
        return
    
    else if (isel==1.or.isel==-1) then
		if (.not.allocated(scatterx)) allocate(scatterx(nx*ny*nz),scattery(nx*ny*nz))
		write(*,*) "Drawing graph, please wait..."
		ii=1
		do k=1,nz
			do j=1,ny
				do i=1,nx
					scatterx(ii)=avgsl2r(i,j,k)
                    scattery(ii)=dg_inter(i,j,k)
					ii=ii+1
				end do
			end do
		end do
		if (isel==-1) then
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"Averaged $sign({\lambda}_2)\rho$ (a.u.)","Averaged ${\delta}g^{inter}$ (a.u.)")
		else
			isavepic=1
			call drawscatter(scatterx,scattery,nx*ny*nz,xmin,xmax,ymin,ymax,1,"Averaged $sign({\lambda}_2)\rho$ (a.u.)","Averaged ${\delta}g^{inter}$ (a.u.)")
			isavepic=0
			write(*,"(a)") " Figure has been saved to "//trim(graphformat)//" file with ""dislin"" prefix in current directory"
		end if
        
    else if (isel==2) then
		open(10,file="output.txt",status="replace")
		write(*,*) "Outputting output.txt..."
		write(10,"(2E16.8)") ((( avgsl2r(i,j,k),dg_inter(i,j,k),k=1,nz),j=1,ny),i=1,nx)
		close(10)
		write(*,"(a)") " Done! Columns 1 and 2 correspond to averaged sign(lambda2)rho and delta-g_inter, respectively"
        
    else if (isel==3) then
        write(*,*) "Exporting averaged delta-g_inter to avgdg_inter.cub..."
        open(10,file="avgdg_inter.cub",status="replace")
		call outcube(dg_inter,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
        write(*,*) "Done!"
        write(*,*) "Exporting averaged sign(lambda2)rho to avgsl2r.cub..."
        open(10,file="avgsl2r.cub",status="replace")
		call outcube(avgsl2r,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
        write(*,*) "Done!"
        
    else if (isel==4) then
        write(*,*) "Exporting averaged RDG to avgRDG.cub..."
        open(10,file="avgRDG.cub",status="replace")
		call outcube(avgRDG,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
        write(*,*) "Done!"
        
    else if (isel==5) then
        call calcexport_TFI(avgdens,ifpsstart,ifpsend)
        
    else if (isel==6) then
		if (iIGMtype==1) write(*,*) "Calculating grid data of TFI(aIGM)..."
		if (iIGMtype==-1) write(*,*) "Calculating grid data of TFI(amIGM)..."
		allocate(TFI_IGM(nx,ny,nz))
        TFI_IGM=0
    
		open(10,file=filename,status="old")
		do ifps=1,ifpsend
			call readxyztrj(10)
			if (ifps<ifpsstart) cycle
			call showprog(ifps,nfps)
			!$OMP PARALLEL DO SHARED(TFI_IGM) PRIVATE(i,j,k,ifrag,gradtmp,grad_inter,IGM_gradnorm_inter,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads)
			do k=1,nz
				do j=1,ny
					do i=1,nx
						call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
						grad_inter=0
						IGM_gradnorm_inter=0
						do ifrag=1,nIGMfrag
							if (iIGMtype==-1) then !amIGM
								call IGMgrad_Hirshpromol(tmpx,tmpy,tmpz,IGMfrag(ifrag,1:IGMfragsize(ifrag)),gradtmp(:),rnouse) !Supports PBC
							else if (iIGMtype==1) then !aIGM 
								call IGMgrad_promol(tmpx,tmpy,tmpz,IGMfrag(ifrag,1:IGMfragsize(ifrag)),gradtmp(:),rnouse) !Supports PBC
							end if
							grad_inter(:)=grad_inter(:)+gradtmp(:)
							IGM_gradnorm_inter=IGM_gradnorm_inter+dsqrt(sum(gradtmp**2))
						end do
						TFI_IGM(i,j,k)=TFI_IGM(i,j,k) + ( IGM_gradnorm_inter-dsqrt(sum(grad_inter**2)) - dg_inter(i,j,k) )**2
					end do
				end do
			end do
			!$OMP END PARALLEL DO
		end do
		close(10)
        
        do k=1,nz
			do j=1,ny
				do i=1,nx
					if (dg_inter(i,j,k)/=0) then
						TFI_IGM(i,j,k)=dsqrt(TFI_IGM(i,j,k)/nfps)/dg_inter(i,j,k)
					else
						TFI_IGM(i,j,k)=0D0
					end if
				end do
			end do
		end do
        if (iIGMtype==1) then
			write(*,*) "Exporting TFI(aIGM) to TFI-aIGM.cub..."
			open(10,file="TFI-aIGM.cub",status="replace")
        else if (iIGMtype==-1) then
			write(*,*) "Exporting TFI(amIGM) to TFI-amIGM.cub..."
			open(10,file="TFI-amIGM.cub",status="replace")
        end if
		call outcube(TFI_IGM,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
        write(*,*) "Done!"
		deallocate(TFI_IGM)
    end if
end do

end subroutine





!!-------------------------------------
!!------ Calculate vdW potential ------
!!-------------------------------------
!For simplicity, only UFF is employed in current code (If using AMBER99 & GAFF, then atom types must be set first)
subroutine vdwpotential
use defvar
use GUI
implicit real*8 (a-h,o-z)
real*8 parmA(ncenter),parmB(ncenter),UFF_A(103),UFF_B(103),tvec(3)
real*8,allocatable :: repulgrid(:,:,:),dispgrid(:,:,:),vdwgrid(:,:,:)
character outcubfile*200,c80tmp*80

write(*,"(/,a)") " !!! If this method is employed in your work, please cite this paper along with Multiwfn original paper:"
write(*,"(a)") " Tian Lu, Qinxue Chen, van der Waals Potential: An Important Complement to Molecular Electrostatic &
&Potential in Studying Intermolecular Interactions. J. Mol. Model., 26, 315 (2020) DOI: 10.1007/s00894-020-04577-0"
write(*,"(a,/)") " Review: Tian Lu, Visualization Analysis of Covalent and Noncovalent Interactions in Real Space, Angew. Chem. Int. Ed., 137, e202504895 (2025) DOI: 10.1002/anie.202504895"

if (ivdwprobe==0) then
    write(*,*) "Input name of probe atom, e.g. Ar"
    read(*,*) c80tmp
    call elename2idx(c80tmp,ivdwprobe)
end if

write(*,*) "Parameters of UFF forcefield are used in this module"
write(*,"(' Element of probe atom: ',a)") ind2name(ivdwprobe)

!call setvdWparm(1,FFtype,parmA,parmB,istatus)
call defineUFFparm(UFF_A,UFF_B)
do iatm=1,ncenter
    parmA(iatm)=UFF_A(a(iatm)%index)
    parmB(iatm)=UFF_B(a(iatm)%index)
end do
parmAj=UFF_A(ivdwprobe)
parmBj=UFF_B(ivdwprobe)
write(*,"(' UFF atomic well depth:',f10.3,' kcal/mol')") parmAj
write(*,"(' UFF atomic radius:    ',f10.3,' Angstrom')") parmBj/2

aug3D=8 !Use 8 Bohr as extension distance because vdW potential largely spread
call setgrid(1,igridsel)
allocate(repulgrid(nx,ny,nz),dispgrid(nx,ny,nz),vdwgrid(nx,ny,nz))
    
write(*,*) "Calculating, please wait..."
ifinish=0
call showprog(0,nz)
!$OMP PARALLEL DO SHARED(ifinish,repulgrid,dispgrid,vdwgrid) PRIVATE(i,j,k,tmpx,tmpy,tmpz,tmprepul,tmpdisp,&
!$OMP iatm,dist,Dij,Xij,ic,jc,kc,icell,jcell,kcell,atmx,atmy,atmz,tvec) schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
            tmprepul=0
            tmpdisp=0
            if (ifPBC==0) then
                do iatm=1,ncenter
                    dist=dsqrt( (a(iatm)%x-tmpx)**2 + (a(iatm)%y-tmpy)**2 + (a(iatm)%z-tmpz)**2 )*b2a
                    if (dist>25) cycle
                    if (dist==0) dist=1E-10 !Avoid singular when grid point happens at nucleus
				    Dij=dsqrt(parmA(iatm)*parmAj) !Well depth
				    Xij=dsqrt(parmB(iatm)*parmBj) !vdW distance
				    tmprepul=tmprepul+Dij*(Xij/dist)**12 !Repulsion
				    tmpdisp=tmpdisp-2*Dij*(Xij/dist)**6 !Dispersion
                end do
            else
                call getpointcell(tmpx,tmpy,tmpz,ic,jc,kc)
                do icell=ic-PBCnx,ic+PBCnx
                    do jcell=jc-PBCny,jc+PBCny
                        do kcell=kc-PBCnz,kc+PBCnz
                            call tvec_PBC(icell,jcell,kcell,tvec)
                            do iatm=1,ncenter
                                atmx=a(iatm)%x+tvec(1)
                                atmy=a(iatm)%y+tvec(2)
                                atmz=a(iatm)%z+tvec(3)
                                dist=dsqrt( (atmx-tmpx)**2 + (atmy-tmpy)**2 + (atmz-tmpz)**2 )*b2a
                                if (dist>25) cycle
                                if (dist==0) dist=1E-10 !Avoid singular when grid point happens at nucleus
				                Dij=dsqrt(parmA(iatm)*parmAj) !Well depth
				                Xij=dsqrt(parmB(iatm)*parmBj) !vdW distance
				                tmprepul=tmprepul+Dij*(Xij/dist)**12 !Repulsion
				                tmpdisp=tmpdisp-2*Dij*(Xij/dist)**6 !Dispersion
                            end do
                        end do
                    end do
                end do
            end if
            repulgrid(i,j,k)=tmprepul
            dispgrid(i,j,k)=tmpdisp
            vdWgrid(i,j,k)=tmprepul+tmpdisp
		end do
	end do
	!$OMP critical
    ifinish=ifinish+1
    call showprog(ifinish,nz)
	!$OMP end critical
end do
!$OMP END PARALLEL DO

write(*,*) "Note: The unit of the grid data is in kcal/mol"
sur_value=1D0
do while(.true.)
	write(*,*)
	write(*,*) "0 Return"
	write(*,*) "1 Show isosurface graph of repulsion potential"
	write(*,*) "2 Show isosurface graph of dispersion potential"
	write(*,*) "3 Show isosurface graph of van der Waals potential"
	write(*,*) "4 Export grid data of repulsion potential as repul.cub in current folder"
	write(*,*) "5 Export grid data of dispersion potential as disp.cub in current folder"
	write(*,*) "6 Export grid data of van der Waals potential as vdW.cub in current folder"
	read(*,*) isel
	if (isel==0) then
		exit
    else if (isel==1.or.isel==2.or.isel==3) then
        if (allocated(cubmat)) deallocate(cubmat)
        allocate(cubmat(nx,ny,nz))
        if (isel==1) cubmat=repulgrid
        if (isel==2) cubmat=dispgrid
        if (isel==3) cubmat=vdwgrid
		call drawisosurgui(1)
        deallocate(cubmat)
	else if (isel==4.or.isel==5.or.isel==6) then
        if (allocated(cubmat)) deallocate(cubmat)
        allocate(cubmat(nx,ny,nz))
        if (isel==4) then
            cubmat=repulgrid
            outcubfile="repul.cub"
        else if (isel==5) then
            cubmat=dispgrid
            outcubfile="disp.cub"
        else if (isel==6) then
            cubmat=vdwgrid
            outcubfile="vdW.cub"
        end if
		open(10,file=outcubfile,status="replace")
		call outcube(cubmat,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
		close(10)
        deallocate(cubmat)
		write(*,"(' Done! Grid data has been exported to ',a,' in current folder')") trim(outcubfile)
    end if
end do
end subroutine