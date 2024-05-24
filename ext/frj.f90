!### Content of this file was contributed by frj, slightly adapted by Tian Lu
    
    
    
!!============================ MBIS ============================!!
!frj addition of new MBIS routines
!
!Wrapper of MBIS module to automatically set radpot and sphpot to proper values
subroutine mbis_wrapper_frj
use defvar
implicit real*8 (a-h,o-z)
nradpotold=radpot
nsphpotold=sphpot
!frj cloned from HI, turn off reducing grid to get 75 x 434 instead of 30 x 170
!if (iautointgrid==1) then
!  	radpot=30
!  	sphpot=170
! 	if (any(a%index>18)) radpot=40
! 	if (any(a%index>36)) radpot=50
! 	if (any(a%index>54)) radpot=60
!end if
!frj accurate grid to test against Gamess with a similar sized grid
!radpot=100
!sphpot=590
call mbis_frj
if (iautointgrid==1) then
	radpot=nradpotold
	sphpot=nsphpotold
end if
end subroutine
!
!
!!--------- Calculate MBIS charge and yield final atomic radial density
!frj : has been tested against grid based integration from Gamess 
subroutine mbis_frj
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
!real*8 molrho(radpot*sphpot),promol(radpot*sphpot),tmpdens(radpot*sphpot),selfdens(radpot*sphpot),molrhoall(ncenter,radpot*sphpot)
real*8 charge(ncenter),lastcharge(ncenter) !Atomic charge of current iter. and last iter.
real*8 beckeweigrid(radpot*sphpot)
!frj new arrays, since Qini and Zini are only implemented up to Ar, there is a maxium of 3 shells, and this is currently hardwired
real*8 shellpop(3,ncenter), zeff(3,ncenter)   !shell populations and shell effective nuclear charges, hardwired to 3 shells maximum 
real*8 shelltmp(3,ncenter), ztmp(3,ncenter)   !new shell populations and shell effective nuclear charges, hardwired to 3 shells maximum 
real*8 wtatm(3,ncenter)                       !shell weights for each atomic shell for each atom
real*8 tmpdens(ncenter,radpot*sphpot)         ! molecular density
real*8 qatom(ncenter), adip(3,ncenter), aquad(6,ncenter), aquadt(6,ncenter)      !atomic charge, dipole and quadrupole, Cartesian and traceless
real*8 shelldip(3,3,ncenter), shellquad(6,3,ncenter)      !shell dipole and (Cartesian) quadrupole, 3 shells hardwired
real*8 mdip(3,3), mquad(6,4), mquadt(6,4)                  ! reconstructed molecular dipole and quadrupole
integer mshell(ncenter)
!frj
character wbsfilename*200,gridfilename*200
integer :: maxcyc=50,ioutmedchg=0,iaddinfo=0
real*8 :: crit=0.0001D0
logical :: ldumpgrids

!Ignore jatm contribution to iatm centered grids if distance between iatm and jatm is larger than 1.5 times of sum of their vdwr
!This can reduce lots of time for large system, the lose of accuracy can be ignored (error is ~0.0001 per atom)
!frj: turn this off for now
!integer :: ignorefar=1
integer :: ignorefar=0
real*8 :: vdwsumcut=2D0
real*8 :: eps=1.0d-14
real*8 :: dencut=1.0d-10

!frj imode is currently unused
!Mode 1 use very low memory but expensive, because most data is computed every iteration
!Mode 2 use large memory but fast, because most data is only computed once at initial stage
!The result of the two modes differ with each other marginally, probably because in mode 1 radial density is related to max(npthigh,nptlow), which is not involved in mode 2
!In principle, result of mode 2 is slightly better
integer :: imode=2
!frj dump raw grid and mbis weights ?
ldumpgrids=.false.
!ldumpgrids=.true.

ntotpot=radpot*sphpot

do while(.true.)
    write(*,*) "     =============== Iterative MBIS (frj) ==============="
    if (iaddinfo==0) write(*,*) "-4 Switch if outputting more information, current: No"
    if (iaddinfo==1) write(*,*) "-4 Switch if outputting more information, current: Yes"
	!if (ignorefar==1) write(*,"(a,f6.3)") " -3 Switch if speeding up calculation using distance cutoff, current: Yes, ratio factor is",vdwsumcut
	!if (ignorefar==0) write(*,*) "-3 Switch if speeding up calculation using distance cutoff, current: No"
	!if (imode==1) write(*,*) "-2 Switch algorithm, current: Slow & low memory requirement"
	!if (imode==2) write(*,*) "-2 Switch algorithm, current: Fast & large memory requirement"
	write(*,*) "1 Start calculation!"
	write(*,"(a,i4)") " 2 Set the maximum number of iterations, current:",maxcyc
	write(*,"(a,f10.6)") " 3 Set convergence criterion of atomic charges, current:",crit
	read(*,*) isel
	if (isel==-4) then
        if (iaddinfo==0) then
            iaddinfo=1
        else
            iaddinfo=0
        end if
	else if (isel==-3) then
        if (ignorefar==1) then
            ignorefar=0
        else
            ignorefar=1
            write(*,*) "Input ratio factor of cutoff, e.g. 2.5"
            write(*,*) "Note: The higher the value, the more accurate the result and the more robust &
            the calculation will be, however the computational cost will be correspondingly higher. The default value is 2.0"
            read(*,*) vdwsumcut
        end if
	else if (isel==-2) then
		if (imode==1) then
			imode=2
		else
			imode=1
			crit=0.001 !mode 1 is more time-consuming, use loose criterion
		end if
	else if (isel==-1) then
		if (ioutmedchg==1) then
			ioutmedchg=0
		else
			ioutmedchg=1
		end if
	else if (isel==0) then
		return
	else if (isel==1) then
		exit
	else if (isel==2) then
		write(*,*) "Input maximum number of iterations, e.g. 30"
		read(*,*) maxcyc
	else if (isel==3) then
		write(*,*) "Input convergence criterion of atomic charges, e.g. 0.001"
		read(*,*) crit
	end if
end do

!====== Start calculation ======!
call walltime(iwalltime1)

!Generate single center integration grid
call gen1cintgrid(gridatmorg,iradcut)
write(*,"(' Radial grids:',i4,'  Angular grids:',i5,'  Total:',i7,'  After pruning:',i7)") radpot,sphpot,radpot*sphpot,radpot*sphpot-iradcut*sphpot

!frj precalculate the molecular density at each grid point
! frj: dump grid info to fort.201
if (ldumpgrids) then
        call path2filename(firstfilename,gridfilename)
        open(201,file=trim(gridfilename)//".grid",status="replace")
        write(201,*)' ### Grid info: ###'
        write(201,*)' iatm, i_grid_point, Grid_x, Grid_y, Grid_z, grid weight, density'
        itmp=ntotpot-iradcut*sphpot
        write(201,*)ncenter,ncenter*itmp
endif
write(*,*) "Calculating molecular density in grid points..."
tmpdens=0.0d0
!frj these loops should be interchanged, keep for now to avoid trouble with the dump-info

if (ldumpgrids) then
    do iatm=1,ncenter
	        gridatm%value=gridatmorg%value   !Weight in this grid point
       	    gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
       	    gridatm%y=gridatmorg%y+a(iatm)%y
       	    gridatm%z=gridatmorg%z+a(iatm)%z
            call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	        do ipt=1+iradcut*sphpot,ntotpot
    !       ipt=1,ntotpot
    !               the molecular density contribution at this point
		        dtmp = fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
                tmpdens(iatm,ipt) = dtmp*gridatm(ipt)%value*beckeweigrid(ipt)
                write(201,"(2i8,3f15.8,2d25.15)")iatm,ipt,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,gridatm(ipt)%value*beckeweigrid(ipt),dtmp
            end do
    end do
else
    !$OMP PARALLEL DO SHARED(tmpdens) PRIVATE(iatm,gridatm,beckeweigrid,dtmp) schedule(dynamic) NUM_THREADS(nthreads)
    do iatm=1,ncenter
	    gridatm%value=gridatmorg%value   !Weight in this grid point
       	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
       	gridatm%y=gridatmorg%y+a(iatm)%y
       	gridatm%z=gridatmorg%z+a(iatm)%z
        call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	    do ipt=1+iradcut*sphpot,ntotpot
		    dtmp = fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
            tmpdens(iatm,ipt) = dtmp*gridatm(ipt)%value*beckeweigrid(ipt)
        end do
    end do
    !$OMP END PARALLEL DO
end if
! do not close 201, as we add converged MBIS parameters later
call walltime(iwalltime2)
write(*,"(' Density calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1
write(*,*) "Done calculating molecular density in grid points..."
!
! frj tmpdens summed over ipt should now be Becke charges..., and indeed it is ...
if (iaddinfo==1) then
qatom = 0.0d0
adip  = 0.0d0
aquad = 0.0d0
do iatm=1,ncenter
        gridatm%x=gridatmorg%x + a(iatm)%x
        gridatm%y=gridatmorg%y + a(iatm)%y
        gridatm%z=gridatmorg%z + a(iatm)%z
        dtmp = 0.0d0
	do ipt=1+iradcut*sphpot,ntotpot
                dx = gridatm(ipt)%x - a(iatm)%x
                dy = gridatm(ipt)%y - a(iatm)%y
                dz = gridatm(ipt)%z - a(iatm)%z
		dtmp = dtmp + tmpdens(iatm,ipt)
		adip(1,iatm) = adip(1,iatm) - tmpdens(iatm,ipt)*dx
		adip(2,iatm) = adip(2,iatm) - tmpdens(iatm,ipt)*dy
		adip(3,iatm) = adip(3,iatm) - tmpdens(iatm,ipt)*dz
		aquad(1,iatm) = aquad(1,iatm) - tmpdens(iatm,ipt)*dx*dx
		aquad(2,iatm) = aquad(2,iatm) - tmpdens(iatm,ipt)*dx*dy
		aquad(3,iatm) = aquad(3,iatm) - tmpdens(iatm,ipt)*dx*dz
		aquad(4,iatm) = aquad(4,iatm) - tmpdens(iatm,ipt)*dy*dy
		aquad(5,iatm) = aquad(5,iatm) - tmpdens(iatm,ipt)*dy*dz
		aquad(6,iatm) = aquad(6,iatm) - tmpdens(iatm,ipt)*dz*dz
        end do
        qatom(iatm) = a(iatm)%charge - dtmp
!        write(*,"('Atom =',i5,' Becke charge and dipole =',f12.5,5x,3f12.5)")iatm,dtmp,dtmpx,dtmpy,dtmpz
end do
         write(*,*)' '
         write(*,*)'  ##### Becke results: ##### '
         write(*,*)' '
         write(*,*)' Atomic charges, un-normalized'
         do iatm=1,ncenter
                 write(*,"(i5,3f15.8)")iatm,qatom(iatm)
         end do
         write(*,*)' Atomic dipoles, in order x, y, z'
         do iatm=1,ncenter
                 write(*,"(i5,3f15.8)")iatm,(adip(j,iatm),j=1,3)
         end do
         write(*,*)' Atomic quadrupoles, Cartesian form, in order xx, xy, xz, yy, yz, zz'
         do iatm=1,ncenter
                 write(*,"(i5,6f15.8)")iatm,(aquad(j,iatm),j=1,6)
         end do
!         convert quadrupole to traceless form
         aquadt = aquad
         write(*,*)' Atomic quadrupoles, Traceless form, in order xx, xy, xz, yy, yz, zz'
         do iatm=1,ncenter
                 qtrace = aquadt(1,iatm) + aquadt(4,iatm) + aquadt(6,iatm)
                 aquadt(1,iatm) = (3.0d0*aquadt(1,iatm) - qtrace) / 2.0d0
                 aquadt(4,iatm) = (3.0d0*aquadt(4,iatm) - qtrace) / 2.0d0
                 aquadt(6,iatm) = (3.0d0*aquadt(6,iatm) - qtrace) / 2.0d0
                 write(*,"(i5,6f15.8)")iatm,(aquadt(j,iatm),j=1,6)
         end do
!
!  construct the molecular multipoles from the atomic ones
         qtot  = 0.d00
         mdip  = 0.d00
         mquad = 0.d00
         do iatm=1,ncenter
                 qtot = qtot + qatom(iatm)
                 mdip(1,1) = mdip(1,1) + qatom(iatm)*a(iatm)%x
                 mdip(2,1) = mdip(2,1) + qatom(iatm)*a(iatm)%y
                 mdip(3,1) = mdip(3,1) + qatom(iatm)*a(iatm)%z
                 mdip(1,2) = mdip(1,2) + adip(1,iatm)
                 mdip(2,2) = mdip(2,2) + adip(2,iatm)
                 mdip(3,2) = mdip(3,2) + adip(3,iatm)
                 mquad(1,1) = mquad(1,1) + qatom(iatm)*a(iatm)%x*a(iatm)%x
                 mquad(2,1) = mquad(2,1) + qatom(iatm)*a(iatm)%x*a(iatm)%y
                 mquad(3,1) = mquad(3,1) + qatom(iatm)*a(iatm)%x*a(iatm)%z
                 mquad(4,1) = mquad(4,1) + qatom(iatm)*a(iatm)%y*a(iatm)%y
                 mquad(5,1) = mquad(5,1) + qatom(iatm)*a(iatm)%y*a(iatm)%z
                 mquad(6,1) = mquad(6,1) + qatom(iatm)*a(iatm)%z*a(iatm)%z
                 mquad(1,2) = mquad(1,2) + adip(1,iatm)*a(iatm)%x
                 mquad(2,2) = mquad(2,2) + adip(1,iatm)*a(iatm)%y
                 mquad(3,2) = mquad(3,2) + adip(1,iatm)*a(iatm)%z
                 mquad(4,2) = mquad(4,2) + adip(2,iatm)*a(iatm)%y
                 mquad(5,2) = mquad(5,2) + adip(2,iatm)*a(iatm)%z
                 mquad(6,2) = mquad(6,2) + adip(3,iatm)*a(iatm)%z
                 mquad(1,2) = mquad(1,2) + adip(1,iatm)*a(iatm)%x
                 mquad(2,2) = mquad(2,2) + adip(2,iatm)*a(iatm)%x
                 mquad(3,2) = mquad(3,2) + adip(3,iatm)*a(iatm)%x
                 mquad(4,2) = mquad(4,2) + adip(2,iatm)*a(iatm)%y
                 mquad(5,2) = mquad(5,2) + adip(3,iatm)*a(iatm)%y
                 mquad(6,2) = mquad(6,2) + adip(3,iatm)*a(iatm)%z
                 mquad(1,3) = mquad(1,3) + aquad(1,iatm)
                 mquad(2,3) = mquad(2,3) + aquad(2,iatm)
                 mquad(3,3) = mquad(3,3) + aquad(3,iatm)
                 mquad(4,3) = mquad(4,3) + aquad(4,iatm)
                 mquad(5,3) = mquad(5,3) + aquad(5,iatm)
                 mquad(6,3) = mquad(6,3) + aquad(6,iatm)
         end do
         do i = 1,3
           do j = 1,2
             mdip(i,3) = mdip(i,3) + mdip(i,j)
           enddo
         enddo
         do i = 1,6
           do j = 1,3
             mquad(i,4) = mquad(i,4) + mquad(i,j)
           enddo
         enddo
         write(*,*)' '
         write(*,"('  Molecular monopole ',f15.8)")qtot
         write(*,*)' Dipole contributions: Atomic rank 0, 1, molecular'
         do i=1,3
                 write(*,"(i5,6f15.8)")i-1,(mdip(j,i),j=1,3)
         end do
         write(*,*)' Cartesian Quadrupole contributions: Atomic rank 0, 1, 2, molecular'
         do i=1,4
                 write(*,"(i5,6f15.8)")i-1,(mquad(j,i),j=1,6)
         end do
!
!        same in traceless form
!        NOTE: if comparing to Gaussian results, Gaussian defines traceless as
!        e.g. xx - trace/3  while Stone defines it as (3xx - trace)/2
!        here the Stone definition is used
         mquadt = mquad
         do j=1,4
                 qtrace = mquad(1,j) + mquad(4,j) + mquad(6,j)
                 mquadt(1,j) = (3.0d0*mquad(1,j) - qtrace) / 2.0d0
                 mquadt(4,j) = (3.0d0*mquad(4,j) - qtrace) / 2.0d0
                 mquadt(6,j) = (3.0d0*mquad(6,j) - qtrace) / 2.0d0
         end do
         write(*,*)' Traceless Quadrupole contributions: Atomic rank 0, 1, 2, molecular'
         do i=1,4
                 write(*,"(i5,6f15.8)")i-1,(mquadt(j,i),j=1,6)
         end do
        call outatmmpl(10,1,qatom(:),adip(:,:),aquad(:,:),aquadt(:,:),mdip(:,:),mquad(:,:),mquadt(:,:) )

        call walltime(iwalltime2x)
         write(*,*)' '
        write(*,"(' Becke integration took up wall clock time',i10,' s')") iwalltime2x-iwalltime2
         write(*,*)' '
         write(*,*)'  ##### end Becke results: ##### '
         write(*,*)' '
         write(*,*)'  ##### Starting MBIS ##### '
end if !iaddinfo

! estimate initial shell values for Zeff and population, only up to Ar at present....
! Znuc = 0 indicate a ghost center, initial as Zeff=1 and with a tiny population
mshell=0
do iatm=1,ncenter
     znuc = a(iatm)%charge
     if (znuc.eq.0.0d0) then
       zeff(1,iatm) = 1.0d0
       shellpop(1,iatm)=1.0d-3
     else
       zeff(1,iatm) = 2.0d0*znuc
     endif
     if (znuc.gt.0.0d0 .and. znuc.le.2.0d0) then
       mshell(iatm) = 1
       shellpop(1,iatm)=znuc
     endif
     if (znuc.gt.2.0d0 .and. znuc.le.10.0d0) then
       mshell(iatm) = 2
       zeff(2,iatm) = 2.0d0
       shellpop(1,iatm)=2.0d0
       shellpop(2,iatm)=znuc-2.0d0
     endif
     if (znuc.gt.10.0d0 .and. znuc.le.18.0d0) then
       mshell(iatm) = 3
       zeff(2,iatm) = 2.0d0*sqrt(znuc)
       zeff(3,iatm) = 2.0d0
       shellpop(1,iatm)=2.0d0
       shellpop(2,iatm)=8.0d0
       shellpop(3,iatm)=znuc-10.0d0
     endif
     if (znuc.gt.18.0d0) then
       write(*,*)' Sorry, no MBIS for atoms beyond Ar yet....'
       exit
     endif
enddo
!write(*,*)'Zeff initial'
!do iatm=1,ncenter
!  write(*,"(2i5,5f12.4)")iatm,mshell(iatm),(zeff(j,iatm),j=1,mshell(iatm))
!enddo
!write(*,*)'Qshell initial'
!do iatm=1,ncenter
!  write(*,"(2i5,5f12.4)")iatm,mshell(iatm),(shellpop(j,iatm),j=1,mshell(iatm))
!enddo
!write(*,*)' '
write(*,*)
write(*,*) "Performing MBIS iterations to refine atomic spaces..."
lastcharge=0.0d0
!Cycle each atom to calculate their charges
do icyc=1,maxcyc
	if (ioutmedchg==1) write(*,*)
	if (icyc==1) then
		write(*,"(' Cycle',i5)") icyc
	else
		write(*,"(' Cycle',i5,'   Maximum change:',f10.6)") icyc,varmax
	end if
	
!       calculate shell weights for each atomic grid point with the current shellpop and zeff
!       calculate the new shell population
        shelltmp = 0.0d0
        do iatm=1,ncenter
               	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
               	gridatm%y=gridatmorg%y+a(iatm)%y
               	gridatm%z=gridatmorg%z+a(iatm)%z
        	do ipt=1+iradcut*sphpot,ntotpot
!        	do ipt=1,ntotpot
                        wtatm=0.0d0
                        wtot = 0.0d0
                        do jatm=1,ncenter
                                dx = gridatm(ipt)%x - a(jatm)%x
                                dy = gridatm(ipt)%y - a(jatm)%y
                                dz = gridatm(ipt)%z - a(jatm)%z
                                dis = dx*dx + dy*dy + dz*dz
                                dis = sqrt(dis)
                                do kshell=1,mshell(jatm)
                                        znuc = zeff(kshell,jatm)
                                        znorm = (znuc**3)/(8.0d0*pi)
                                        qshell = shellpop(kshell,jatm)
!                                       note that we do not need to multiply with the grid weight, as wtatm will only be used as relative values
                                        tmp = qshell*znorm*exp(-znuc*dis)
                                        if (tmp.lt.dencut) tmp = 0.0d0
                                        wtatm(kshell,jatm) = wtatm(kshell,jatm) + tmp
                                        wtot = wtot + wtatm(kshell,jatm)
                                end do
                        end do
                        tmpden = tmpdens(iatm,ipt)
                        if (wtot.gt.0.0d0 .and. tmpden.gt.eps) then
                                do jatm=1,ncenter
                                        do kshell=1,mshell(jatm)
                                                shelltmp(kshell,jatm) = shelltmp(kshell,jatm) + wtatm(kshell,jatm)*tmpden/wtot
                                        end do
                                end do
                        endif
                end do
        end do
!       condense to atoms and possibly print current values
!        write(*,*)' '
!        write(*,*)' Qshell, current, previous, difference'
         do iatm=1,ncenter
                 electmp = 0.0d0
                 do kshell=1,mshell(iatm)
                         electmp = electmp + shelltmp(kshell,iatm)
                         tmp  = shelltmp(kshell,iatm) - shellpop(kshell,iatm)
!                         write(*,"(2i5,3f15.8)")iatm,kshell,shelltmp(kshell,iatm),shellpop(kshell,iatm),tmp
                 end do
                 charge(iatm) = a(iatm)%charge - electmp
         end do
!         write(*,*)' Qatom , current, previous, difference'
!         do iatm=1,ncenter
!                 tmp = charge(iatm) - lastcharge(iatm)
!                 write(*,"(i5,3f15.8)")iatm,charge(iatm),lastcharge(iatm),tmp
!         end do
!         do iatm=1,ncenter
! 		if (ioutmedchg==1) write(*,"(' Charge of atom',i5,'(',a2,')',': ',f12.6,'  Delta:',f12.6)") &
! 		jatm,a(jatm)%name,charge(iatm),charge(iatm)-lastcharge(iatm)
!         end do
!        write(*,*)' '

!       check for convergence
	varmax=maxval(abs(charge-lastcharge))
	if (varmax<crit) then
		write(*,"(a,f10.6)") " All atomic charges have converged to criterion of",crit
		write(*,"(' Sum of all charges:',f14.8)") sum(charge)
		call normalize_atmchg(charge(:))
        call printatmchg(charge(:))
		exit
	else
		if (icyc==maxcyc) then
			write(*,"(/,' Convergence failed within',i4,' cycles!')") maxcyc
			exit
		end if
	end if

!       same procedure for updating the effective nuclear charges with the new shell populations
!       calculate shell weights for each atomic grid point with the current shellpop and zeff
!       calculate the new effective nuclear charge
        shellpop = shelltmp
	lastcharge=charge
        ztmp = 0.0d0
        do iatm=1,ncenter
               	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
               	gridatm%y=gridatmorg%y+a(iatm)%y
               	gridatm%z=gridatmorg%z+a(iatm)%z
        	do ipt=1+iradcut*sphpot,ntotpot
!        	do ipt=1,ntotpot
                        wtatm=0.0d0
                        wtot = 0.0d0
                        do jatm=1,ncenter
                                dx = gridatm(ipt)%x - a(jatm)%x
                                dy = gridatm(ipt)%y - a(jatm)%y
                                dz = gridatm(ipt)%z - a(jatm)%z
                                dis = dx*dx + dy*dy + dz*dz
                                dis = sqrt(dis)
                                do kshell=1,mshell(jatm)
                                        znuc = zeff(kshell,jatm)
                                        znorm = (znuc**3)/(8.0d0*pi)
                                        qshell = shellpop(kshell,jatm)
                                        tmp = qshell*znorm*exp(-znuc*dis)
                                        if (tmp.lt.dencut) tmp = 0.0d0
                                        wtatm(kshell,jatm) = wtatm(kshell,jatm) + tmp
                                        wtot = wtot + wtatm(kshell,jatm)
                                end do
                        end do
                        tmpden = tmpdens(iatm,ipt)
                        if (wtot.gt.0.0d0 .and. tmpden.gt.eps) then
                                do jatm=1,ncenter
                                        dx = gridatm(ipt)%x - a(jatm)%x
                                        dy = gridatm(ipt)%y - a(jatm)%y
                                        dz = gridatm(ipt)%z - a(jatm)%z
                                        dis = dx*dx + dy*dy + dz*dz
                                        dis = sqrt(dis)
                                        do kshell=1,mshell(jatm)
                                                ztmp(kshell,jatm) = ztmp(kshell,jatm) + wtatm(kshell,jatm)*dis*tmpden/wtot
                                        end do
                                end do
                        endif
                end do
        end do
!       convert to actual Zeff and possibly print current values
!        write(*,*)' Zshell, current, previous, difference'
        do iatm=1,ncenter
                do kshell=1,mshell(iatm)
                        ztmp(kshell,iatm) = (3.0d0*shellpop(kshell,iatm))/ztmp(kshell,iatm)
                        tmp  = ztmp(kshell,iatm) - zeff(kshell,iatm)
!                        write(*,"(2i5,3f15.8)")iatm,kshell,ztmp(kshell,iatm),zeff(kshell,iatm),tmp
                end do
        end do
!        write(*,*)' '
!       update zeff and return for another cycle
        zeff = ztmp
end do

call walltime(iwalltime3)
if (iaddinfo==1) then
    write(*,*)
    write(*,"(' MBIS iterations took up wall clock time',i10,' s')") iwalltime3-iwalltime2x
end if
write(*,"(' Calculation took up wall clock time',i10,' s')") iwalltime3-iwalltime1
call outatmchg(10,charge(:))
!frj: add converged MBIS parameters to fort.201, and close
if (ldumpgrids) then
        write(201,*)' MBIS parameters, iatm, kshell, Axyz, Zeff, Npop'
        write(201,'(20i5)')(mshell(i),i=1,ncenter)
        do iatm=1,ncenter
                do kshell=1,mshell(iatm)
                        write(201,"(2i8,3f15.8,2d25.15)")iatm,kshell,a(iatm)%x,a(iatm)%y,a(iatm)%z,zeff(kshell,iatm),shellpop(kshell,iatm)
                end do
        end do
        write(*,*)' ### Grid info and MBIS parameters dumped to file.grid ###'
        close(201)
!frj end dump
endif

!
if (ldumpgrids) then
        call path2filename(firstfilename,wbsfilename)
        open(200,file=trim(wbsfilename)//".wbs",status="replace")
        write(200,*)' ### Grid info: ###', ntotpot-iradcut*sphpot
        write(200,*)' iatm, i_grid_point, jatm, jshell, Grid_x, Grid_y, Grid_z, MBIS weight'
endif
!
if (iaddinfo==1) then
        write(*,*)
        write(*,*)' Calculating MBIS multipole moments up to rank 2'
        shelltmp = 0.0d0
        shelldip = 0.0d0
        shellquad = 0.0d0
        do iatm=1,ncenter
               	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
               	gridatm%y=gridatmorg%y+a(iatm)%y
               	gridatm%z=gridatmorg%z+a(iatm)%z
        	do ipt=1+iradcut*sphpot,ntotpot
!        	do ipt=1,ntotpot
                        wtatm=0.0d0
                        wtot = 0.0d0
                        do jatm=1,ncenter
                                dx = gridatm(ipt)%x - a(jatm)%x
                                dy = gridatm(ipt)%y - a(jatm)%y
                                dz = gridatm(ipt)%z - a(jatm)%z
                                dis = dx*dx + dy*dy + dz*dz
                                dis = sqrt(dis)
                                do kshell=1,mshell(jatm)
                                        znuc = zeff(kshell,jatm)
                                        znorm = (znuc**3)/(8.0d0*pi)
                                        qshell = shellpop(kshell,jatm)
!                                       note that we do not need to multiply with the grid weight, as wtatm will only be used as relative values
                                        tmp = qshell*znorm*exp(-znuc*dis)
                                        if (tmp.lt.dencut) tmp = 0.0d0
                                        wtatm(kshell,jatm) = wtatm(kshell,jatm) + tmp
                                        wtot = wtot + wtatm(kshell,jatm)
                                end do
                        end do
                        tmpden = tmpdens(iatm,ipt)
                         if (wtot.gt.0.0d0 .and. tmpden.gt.eps) then
                                do jatm=1,ncenter
                                        dx = gridatm(ipt)%x - a(jatm)%x
                                        dy = gridatm(ipt)%y - a(jatm)%y
                                        dz = gridatm(ipt)%z - a(jatm)%z
                                        do kshell=1,mshell(jatm)
if (ldumpgrids) then
!                                               for dumping grid and MBIS weight points
                                                wtmp = wtatm(kshell,jatm)/wtot
                                                write(200,"(4i8,3f15.8,d25.15)")iatm,ipt,jatm,kshell,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,wtmp
!                                               end of hack for dumping grid
endif
                                                wtmp = wtatm(kshell,jatm)*tmpden/wtot
                                                shelltmp(kshell,jatm) = shelltmp(kshell,jatm) + wtmp
                                                shelldip(1,kshell,jatm) = shelldip(1,kshell,jatm) + wtmp*dx
                                                shelldip(2,kshell,jatm) = shelldip(2,kshell,jatm) + wtmp*dy
                                                shelldip(3,kshell,jatm) = shelldip(3,kshell,jatm) + wtmp*dz
                                                shellquad(1,kshell,jatm) = shellquad(1,kshell,jatm) + wtmp*dx*dx
                                                shellquad(2,kshell,jatm) = shellquad(2,kshell,jatm) + wtmp*dx*dy
                                                shellquad(3,kshell,jatm) = shellquad(3,kshell,jatm) + wtmp*dx*dz
                                                shellquad(4,kshell,jatm) = shellquad(4,kshell,jatm) + wtmp*dy*dy
                                                shellquad(5,kshell,jatm) = shellquad(5,kshell,jatm) + wtmp*dy*dz
                                                shellquad(6,kshell,jatm) = shellquad(6,kshell,jatm) + wtmp*dz*dz
                                        end do
                                end do
                        else
if (ldumpgrids) then
!                               another hack for making sure the full grid is dumped
                                do jatm=1,ncenter
                                        dx = gridatm(ipt)%x - a(jatm)%x
                                        dy = gridatm(ipt)%y - a(jatm)%y
                                        dz = gridatm(ipt)%z - a(jatm)%z
                                        do kshell=1,mshell(jatm)
!                                                wtmp = wtatm(kshell,jatm)/wtot
                                                wtmp = 0.0d0
                                                write(200,"(4i8,3f15.8,2d25.15)")iatm,ipt,jatm,kshell,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,wtmp
                                        end do
                                end do
endif

                        endif
                end do
        end do
if (ldumpgrids) then
        write(*,*)' ### Grid info dumped to file.wbs ###'
        close(200)
endif
!         write(*,*)' Qshell, char, dip, quad'
         adip  = 0.0d0
         aquad = 0.0d0
         do iatm=1,ncenter
                 electmp = 0.0d0
                 do kshell=1,mshell(iatm)
                         electmp = electmp + shelltmp(kshell,iatm)
!                         write(*,"(2i5,3f15.8)")iatm,kshell,shelltmp(kshell,iatm)
!                         write(*,"(2i5,3f15.8)")iatm,kshell,(shelldip(j,kshell,iatm),j=1,3)
!                         write(*,"(2i5,6f15.8)")iatm,kshell,(shellquad(j,kshell,iatm),j=1,6)
                         adip(1,iatm) = adip(1,iatm) - shelldip(1,kshell,iatm)
                         adip(2,iatm) = adip(2,iatm) - shelldip(2,kshell,iatm)
                         adip(3,iatm) = adip(3,iatm) - shelldip(3,kshell,iatm)
                         aquad(1,iatm) = aquad(1,iatm) - shellquad(1,kshell,iatm)
                         aquad(2,iatm) = aquad(2,iatm) - shellquad(2,kshell,iatm)
                         aquad(3,iatm) = aquad(3,iatm) - shellquad(3,kshell,iatm)
                         aquad(4,iatm) = aquad(4,iatm) - shellquad(4,kshell,iatm)
                         aquad(5,iatm) = aquad(5,iatm) - shellquad(5,kshell,iatm)
                         aquad(6,iatm) = aquad(6,iatm) - shellquad(6,kshell,iatm)
                 end do
                 charge(iatm) = a(iatm)%charge - electmp
         end do
         write(*,*)' Atomic charges, un-normalized'
         do iatm=1,ncenter
                 write(*,"(i5,3f15.8)")iatm,charge(iatm)
         end do
         write(*,*)' Atomic dipoles, in order x, y, z'
         do iatm=1,ncenter
                 write(*,"(i5,3f15.8)")iatm,(adip(j,iatm),j=1,3)
         end do
         write(*,*)' Atomic quadrupoles, Cartesian form, in order xx, xy, xz, yy, yz, zz'
         do iatm=1,ncenter
                 write(*,"(i5,6f15.8)")iatm,(aquad(j,iatm),j=1,6)
         end do
!         convert quadrupole to traceless form
         aquadt = aquad
         write(*,*)' Atomic quadrupoles, Traceless form, in order xx, xy, xz, yy, yz, zz'
         do iatm=1,ncenter
                 qtrace = aquadt(1,iatm) + aquadt(4,iatm) + aquadt(6,iatm)
                 aquadt(1,iatm) = (3.0d0*aquadt(1,iatm) - qtrace) / 2.0d0
                 aquadt(4,iatm) = (3.0d0*aquadt(4,iatm) - qtrace) / 2.0d0
                 aquadt(6,iatm) = (3.0d0*aquadt(6,iatm) - qtrace) / 2.0d0
                 write(*,"(i5,6f15.8)")iatm,(aquadt(j,iatm),j=1,6)
         end do
!
!  construct the molecular multipoles from the atomic ones
         qtot  = 0.d00
         mdip  = 0.d00
         mquad = 0.d00
         do iatm=1,ncenter
                 qtot = qtot + charge(iatm)
                 mdip(1,1) = mdip(1,1) + charge(iatm)*a(iatm)%x
                 mdip(2,1) = mdip(2,1) + charge(iatm)*a(iatm)%y
                 mdip(3,1) = mdip(3,1) + charge(iatm)*a(iatm)%z
                 mdip(1,2) = mdip(1,2) + adip(1,iatm)
                 mdip(2,2) = mdip(2,2) + adip(2,iatm)
                 mdip(3,2) = mdip(3,2) + adip(3,iatm)
                 mquad(1,1) = mquad(1,1) + charge(iatm)*a(iatm)%x*a(iatm)%x
                 mquad(2,1) = mquad(2,1) + charge(iatm)*a(iatm)%x*a(iatm)%y
                 mquad(3,1) = mquad(3,1) + charge(iatm)*a(iatm)%x*a(iatm)%z
                 mquad(4,1) = mquad(4,1) + charge(iatm)*a(iatm)%y*a(iatm)%y
                 mquad(5,1) = mquad(5,1) + charge(iatm)*a(iatm)%y*a(iatm)%z
                 mquad(6,1) = mquad(6,1) + charge(iatm)*a(iatm)%z*a(iatm)%z
                 mquad(1,2) = mquad(1,2) + adip(1,iatm)*a(iatm)%x
                 mquad(2,2) = mquad(2,2) + adip(1,iatm)*a(iatm)%y
                 mquad(3,2) = mquad(3,2) + adip(1,iatm)*a(iatm)%z
                 mquad(4,2) = mquad(4,2) + adip(2,iatm)*a(iatm)%y
                 mquad(5,2) = mquad(5,2) + adip(2,iatm)*a(iatm)%z
                 mquad(6,2) = mquad(6,2) + adip(3,iatm)*a(iatm)%z
                 mquad(1,2) = mquad(1,2) + adip(1,iatm)*a(iatm)%x
                 mquad(2,2) = mquad(2,2) + adip(2,iatm)*a(iatm)%x
                 mquad(3,2) = mquad(3,2) + adip(3,iatm)*a(iatm)%x
                 mquad(4,2) = mquad(4,2) + adip(2,iatm)*a(iatm)%y
                 mquad(5,2) = mquad(5,2) + adip(3,iatm)*a(iatm)%y
                 mquad(6,2) = mquad(6,2) + adip(3,iatm)*a(iatm)%z
                 mquad(1,3) = mquad(1,3) + aquad(1,iatm)
                 mquad(2,3) = mquad(2,3) + aquad(2,iatm)
                 mquad(3,3) = mquad(3,3) + aquad(3,iatm)
                 mquad(4,3) = mquad(4,3) + aquad(4,iatm)
                 mquad(5,3) = mquad(5,3) + aquad(5,iatm)
                 mquad(6,3) = mquad(6,3) + aquad(6,iatm)
         end do
         do i = 1,3
           do j = 1,2
             mdip(i,3) = mdip(i,3) + mdip(i,j)
           enddo
         enddo
         do i = 1,6
           do j = 1,3
             mquad(i,4) = mquad(i,4) + mquad(i,j)
           enddo
         enddo
         write(*,*)' '
         write(*,"('  Molecular monopole ',f15.8)")qtot
         write(*,*)' Dipole contributions: Atomic rank 0, 1, molecular'
         do i=1,3
                 write(*,"(i5,6f15.8)")i-1,(mdip(j,i),j=1,3)
         end do
         write(*,*)' Cartesian Quadrupole contributions: Atomic rank 0, 1, 2, molecular'
         do i=1,4
                 write(*,"(i5,6f15.8)")i-1,(mquad(j,i),j=1,6)
         end do
!
!        same in traceless form
!        NOTE: if comparing to Gaussian results, Gaussian defines traceless as
!        e.g. xx - trace/3  while Stone defines it as (3xx - trace)/2
!        here the Stone definition is used
         mquadt = mquad
         do j=1,4
                 qtrace = mquad(1,j) + mquad(4,j) + mquad(6,j)
                 mquadt(1,j) = (3.0d0*mquad(1,j) - qtrace) / 2.0d0
                 mquadt(4,j) = (3.0d0*mquad(4,j) - qtrace) / 2.0d0
                 mquadt(6,j) = (3.0d0*mquad(6,j) - qtrace) / 2.0d0
         end do
         write(*,*)' Traceless Quadrupole contributions: Atomic rank 0, 1, 2, molecular'
         do i=1,4
                 write(*,"(i5,6f15.8)")i-1,(mquadt(j,i),j=1,6)
         end do
!
call outatmmpl(10,2,charge(:),adip(:,:),aquad(:,:),aquadt(:,:),mdip(:,:),mquad(:,:),mquadt(:,:) )
end if !iaddinfo
end subroutine

!---- Output atomic multipole to a .mpl file
!ifileid is the file id that can be used in this routine
subroutine outatmmpl(ifileid,ilabel,charge,adip,aquad,aquadt,mdip,mquad,mquadt)
use defvar
use util
integer ifileid,ilabel,i
real*8 dtmp
real*8 charge(ncenter)
real*8 adip(3,ncenter), aquad(6,ncenter),aquadt(6,ncenter)
real*8 mdip(3,ncenter), mquad(6,ncenter),mquadt(6,ncenter)
character selectyn,chgfilename*200
call path2filename(firstfilename,chgfilename)
if (ilabel.eq.1) write(*,"(a)") " Output atomic multipoles to "//trim(chgfilename)//".becke_mpl in current folder? (y/n)"
if (ilabel.eq.2) write(*,"(a)") " Output atomic multipoles to "//trim(chgfilename)//".mbis_mpl in current folder? (y/n)"
read(*,*) selectyn
if (selectyn=="y".or.selectyn=="Y") then
	if (ilabel.eq.1) open(ifileid,file=trim(chgfilename)//".becke_mpl",status="replace")
	if (ilabel.eq.2) open(ifileid,file=trim(chgfilename)//".mbis_mpl",status="replace")
        write(ifileid,*)' Atomic charges, un-normalized'
	do i=1,ncenter
		write(ifileid,"(a,f15.8)") a(i)%name,charge(i)
	end do
        write(ifileid,*)' Atomic dipoles, in order x, y, z'
	do i=1,ncenter
		write(ifileid,"(a,3f15.8)") a(i)%name,(adip(j,i),j=1,3)
	end do
        write(ifileid,*)' Atomic quadrupoles, Cartesian form, in order xx, xy, xz, yy, yz, zz'
	do i=1,ncenter
		write(ifileid,"(a,6f15.8)") a(i)%name,(aquad(j,i),j=1,6)
	end do
        write(ifileid,*)' Atomic quadrupoles, Traceless form, in order xx, xy, xz, yy, yz, zz'
	do i=1,ncenter
		write(ifileid,"(a,6f15.8)") a(i)%name,(aquadt(j,i),j=1,6)
	end do
        write(ifileid,*)' '
        write(ifileid,*)' Atomic to molecular condensed quantities'
        do i=1,ncenter
                dtmp = dtmp + charge(i)
        enddo
	write(ifileid,"(a,f15.8)")' Molecular charge',dtmp
        write(ifileid,*)' Dipole contributions: Atomic rank 0, 1, molecular'
        do i=1,3
                write(ifileid,"(i5,6f15.8)")i-1,(mdip(j,i),j=1,3)
        end do
        write(ifileid,*)' Cartesian Quadrupole contributions: Atomic rank 0, 1, 2, molecular'
        do i=1,4
                write(ifileid,"(i5,6f15.8)")i-1,(mquad(j,i),j=1,6)
        end do
        write(ifileid,*)' Traceless Quadrupole contributions: Atomic rank 0, 1, 2, molecular'
        do i=1,4
                write(ifileid,"(i5,6f15.8)")i-1,(mquadt(j,i),j=1,6)
        end do
	close(ifileid)
	if (ilabel.eq.1)write(*,"(a)") " Atomic multipoles have been saved to "//trim(chgfilename)//".becke_mpl in current folder"
	if (ilabel.eq.2)write(*,"(a)") " Atomic multipoles have been saved to "//trim(chgfilename)//".mbis_mpl in current folder"
	write(*,"(a)") " First (unormalized) charges, then dipoles (x,y,z), then Cartesian quadrupoles (xx,xy,xz,yy,yz,zz), then Traceless quadrupoles, units are au"
	write(*,"(a)") " Then the atom to molecule condensed quantities"
end if
end subroutine

