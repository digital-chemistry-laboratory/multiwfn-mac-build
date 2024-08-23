!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================================================================
!**************************************************************************
!  The codes in this file are written specific for Shubin Liu's project
!**************************************************************************
!==========================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fuzzySBL
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
integer,parameter :: nfunc=8,nquant=7 !The number of real space function, the number of quantities to be calculated
real*8 Pvec(ncenter),atmBeckewei(radpot*sphpot),atmHirshwei(radpot*sphpot,ncenter)
real*8 atmcontri(ncenter,0:nfunc,0:nquant)
real*8 promol(radpot*sphpot),atomdens(radpot*sphpot,ncenter),selfdens(radpot*sphpot)
real*8 promolgrad(radpot*sphpot,3),atomgrad(radpot*sphpot,3)
real*8 funcval(radpot*sphpot,0:nfunc),funcgrdn(radpot*sphpot,0:nfunc) !Store function value and gradient norm, respectively
real*8 funcref(radpot*sphpot,0:nfunc) !Store function value of promolecular state (reference state)
real*8 covr_becke(0:nelesupp) !Covalent radii used for Becke partition
type(content) gridatm(radpot*sphpot)
real*8 potx(sphpot),poty(sphpot),potz(sphpot),potw(sphpot)
real*8 arrtmp(nfunc),arrtmp2(nfunc),gradrho(3)
character(len=40) functionname(0:nfunc),quantityname(0:nquant)

if (ifiletype/=2.and.ifiletype/=3) then
	write(*,*) "Note: Using .wfn or .wfx file will make calculation much faster"
end if

nbeckeiter=3
expcutoff=1
covr_becke=covr_tianlu

functionname(0)=" rho"
functionname(1)=" rho/rho0"
functionname(2)=" |der_rho|/rho^(4/3)"
functionname(3)=" der2rho/rho^(5/3)"
functionname(4)=" (tau-t_w)/t_TF"
functionname(5)=" Xi part of SEDD"
functionname(6)=" Theta part of DORI"
functionname(7)=" Spin density"
functionname(8)=" (tau-t_w)/t_w"
quantityname(0)=" Function itself" !Special
quantityname(1)=" Shannon entropy"
quantityname(2)=" Fisher information"
quantityname(3)=" Onicescu information energy of order 2"
quantityname(4)=" Onicescu information energy of order 3"
quantityname(5)=" Information gain"
quantityname(6)=" Relative Renyi entropy of orders 2"
quantityname(7)=" Relative Renyi entropy of orders 3"

write(*,*) "Select partition method:"
write(*,"(a)") " 1 Becke partition"
write(*,"(a)") " 2 Hirshfeld partition"
read(*,*) ipartition

invfunc=0
write(*,*) "If replacing real space functions with their reciprocals?  0=No  1=Yes"
read(*,*) invfunc 

atmcontri=0D0
checkacc=0
call setpromol
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
write(*,*) "Please wait..."
write(*,*)
call walltime(nwalltime1)
call Lebedevgen(sphpot,potx,poty,potz,potw)

! call valaryyLSB(1D0,1D0,1D0,arrtmp,rho,rhogrdn)
! write(*,*) arrtmp(4)
! pause

do iatm=1,ncenter !! Cycle each atom
	write(*,"(/,' Calculating:',i6,'   /',i6)") iatm,ncenter

	!! Prepare grid points on current center
	iradcut=0 !Before where the radial points will be cut
	parm=1D0
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
	gridatm%x=gridatm%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatm%y+a(iatm)%y
	gridatm%z=gridatm%z+a(iatm)%z
	
	
	!! Calculate function values of reference state, including promolecular density and density of present atom, function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
	!index=2: |der_rho|/rho^(4/3)
	!index=3: der2rho/rho^(5/3)
	!index=4: (tau-t_w)/t_TF
	!index=5: Part of SEDD
	!index=6: Part of DORI
	!index=7: Spin density
	!index=8: (tau-t_w)/t_w
	funcref=0D0
	promol=0D0
	promolgrad=0D0
	do jatm=1,ncenter_org
		call dealloall(0)
		call readwfn(custommapname(jatm),1)
		!$OMP parallel do shared(atomdens,atomgrad,funcref) private(ipt,rnowx,rnowy,rnowz,arrtmp,rho,gradrho) num_threads(nthreads)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			rnowx=gridatm(ipt)%x
			rnowy=gridatm(ipt)%y
			rnowz=gridatm(ipt)%z
			call valaryyLSB(rnowx,rnowy,rnowz,arrtmp,rho,gradrho)
			atomdens(ipt,jatm)=rho
			atomgrad(ipt,:)=gradrho
			funcref(ipt,0)=funcref(ipt,0)+rho
			funcref(ipt,1)=funcref(ipt,1)+1D0 !i.e. rho/rho0, the value is unity if present system is a single atom
			funcref(ipt,2:nfunc)=funcref(ipt,2:nfunc)+arrtmp(2:nfunc)
		end do
		!$OMP end parallel do
		promol=promol+atomdens(:,jatm)
		promolgrad=promolgrad+atomgrad(:,:)
		if (jatm==iatm) selfdens=atomdens(:,jatm)
	end do
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule)
	
    
    !! Calculate function values and gradient norm for present molecule, then store them in funcval/funcgrdn. Function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
	!index=2: |der_rho|/rho^(4/3)
	!index=3: der2rho/rho^(5/3)
	!index=4: (tau-t_w)/t_TF
	!index=5: Part of SEDD
	!index=6: Part of DORI
	!index=7: Spin density
	!index=8: (tau-t_w)/t_w
	!$OMP parallel do shared(funcval,funcgrdn) private(ipt,arrtmp,arrtmp2,rho,gradrho,idir) num_threads(nthreads)
	do ipt=1+iradcut*sphpot,radpot*sphpot
		call valgradarrLSB(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z,arrtmp,arrtmp2,rho,gradrho)
		funcval(ipt,0)=rho
		funcgrdn(ipt,0)=dsqrt(sum(gradrho**2))
		!Note: see word document for explicit expression of |der(rho/rho0)|
		funcval(ipt,1)=rho/promol(ipt)
		funcgrdn(ipt,1)=0D0
		do idir=1,3
			funcgrdn(ipt,1)=funcgrdn(ipt,1)+( gradrho(idir)/promol(ipt)-rho/promol(ipt)**2*promolgrad(ipt,idir) )**2
		end do
		funcgrdn(ipt,1)=dsqrt(funcgrdn(ipt,1))
		funcval(ipt,2:nfunc)=arrtmp(2:nfunc)
		funcgrdn(ipt,2:nfunc)=arrtmp2(2:nfunc)
	end do
	!$OMP end parallel do
    
	
	!! Calculate atomic Becke weight for present atom at all of its points
	!$OMP parallel do shared(atmpartwei) private(i) num_threads(nthreads) schedule(dynamic)
	do i=1+iradcut*sphpot,radpot*sphpot
        call Beckeatmwei(iatm,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,atmBeckewei(i),covr_becke,nbeckeiter)
	end do
	!$OMP end parallel do
	
	
	!Calculate Hirshfeld weights of all atoms
	if (ipartition==2) then
		do jatm=1,ncenter
			do i=1+iradcut*sphpot,radpot*sphpot
				if (promol(i)/=0D0) then
					atmHirshwei(i,jatm)=atomdens(i,jatm)/promol(i)
				else
					atmHirshwei(i,jatm)=0D0
				end if
			end do
		end do
	end if
	
	
    !! Calculate all quantities contributed by this atom
    if (ipartition==1) then !Atomic grid integration for Becke partition
		do ipt=1+iradcut*sphpot,radpot*sphpot
			weitot=atmBeckewei(ipt)*gridatm(ipt)%value
			do ifunc=0,nfunc
				if (invfunc==0) then
					fv=funcval(ipt,ifunc)
					fgrdn=funcgrdn(ipt,ifunc)
					fref=funcref(ipt,ifunc)
				else if (invfunc==1) then
					fv=1/funcval(ipt,ifunc)
					fgrdn=1/funcgrdn(ipt,ifunc)
					fref=1/funcref(ipt,ifunc)
				end if
				!Value of function itself
				atmcontri(iatm,ifunc,0)=atmcontri(iatm,ifunc,0)+weitot*fv
				!Shannon entropy
				atmcontri(iatm,ifunc,1)=atmcontri(iatm,ifunc,1)+weitot*( -fv*log(fv) )
				!Fisher information
				atmcontri(iatm,ifunc,2)=atmcontri(iatm,ifunc,2)+weitot*( fgrdn**2 /fv )
				!Onicescu information energy of order 2
				atmcontri(iatm,ifunc,3)=atmcontri(iatm,ifunc,3)+weitot*( fv**2 )
				!Onicescu information energy of order 3
				atmcontri(iatm,ifunc,4)=atmcontri(iatm,ifunc,4)+weitot*( fv**3 )
				!Information gain
				atmcontri(iatm,ifunc,5)=atmcontri(iatm,ifunc,5)+weitot*( fv*log(fv/fref) )
				!Relative Renyi entropy of orders 2
				atmcontri(iatm,ifunc,6)=atmcontri(iatm,ifunc,6)+weitot*( fv**2/fref )
				!Relative Renyi entropy of orders 3
				atmcontri(iatm,ifunc,7)=atmcontri(iatm,ifunc,7)+weitot*( fv**3/fref**2 )
			end do
		end do
	    
		write(*,"(' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
		do ifunc=0,nfunc
			write(*,"(/,' Function:',a)") trim(functionname(ifunc))
			do iquant=0,nquant !0 corresponds to function itself
				write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
			end do
		end do
		
	else if (ipartition==2) then !Molecular grid integration for Hirshfeld partition
		do jatm=1,ncenter
			do ipt=1+iradcut*sphpot,radpot*sphpot
				weitot=atmBeckewei(ipt)*atmHirshwei(ipt,jatm)*gridatm(ipt)%value
				do ifunc=0,nfunc
					if (invfunc==0) then
						fv=funcval(ipt,ifunc)
						fgrdn=funcgrdn(ipt,ifunc)
						fref=funcref(ipt,ifunc)
					else if (invfunc==1) then
						fv=1/funcval(ipt,ifunc)
						fgrdn=1/funcgrdn(ipt,ifunc)
						fref=1/funcref(ipt,ifunc)
					end if
					!Value of function itself
					atmcontri(jatm,ifunc,0)=atmcontri(jatm,ifunc,0)+weitot*fv
					!Shannon entropy
					atmcontri(jatm,ifunc,1)=atmcontri(jatm,ifunc,1)+weitot*( -fv*log(fv) )
					!Fisher information
					atmcontri(jatm,ifunc,2)=atmcontri(jatm,ifunc,2)+weitot*( fgrdn**2 /fv )
					!Onicescu information energy of order 2
					atmcontri(jatm,ifunc,3)=atmcontri(jatm,ifunc,3)+weitot*( fv**2 )
					!Onicescu information energy of order 3
					atmcontri(jatm,ifunc,4)=atmcontri(jatm,ifunc,4)+weitot*( fv**3 )
					!Information gain
					atmcontri(jatm,ifunc,5)=atmcontri(jatm,ifunc,5)+weitot*( fv*log(fv/fref) )
					!Relative Renyi entropy of orders 2
					atmcontri(jatm,ifunc,6)=atmcontri(jatm,ifunc,6)+weitot*( fv**2/fref )
					!Relative Renyi entropy of orders 3
					atmcontri(jatm,ifunc,7)=atmcontri(jatm,ifunc,7)+weitot*( fv**3/fref**2 )
				end do
			end do
		end do
	end if
	
end do !End cycling atoms

!Output atomic contribution obtained from molecular grid integration
if (ipartition==2) then
	do iatm=1,ncenter
		write(*,"(/,' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
		do ifunc=0,nfunc
			write(*,"(/,' Function:',a)") trim(functionname(ifunc))
			do iquant=0,nquant !0 corresponds to function itself
				write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
			end do
		end do
	end do
end if

write(*,*)
write(*,*) "========================================="
write(*,*) "             Overall result"
write(*,*) "========================================="
do ifunc=0,nfunc
	write(*,"(/,' Function:',a)") trim(functionname(ifunc))
	do iquant=0,nquant
		write(*,"(a,':',1E16.8)") quantityname(iquant),sum(atmcontri(:,ifunc,iquant))
	end do
end do

call walltime(nwalltime2)
write(*,"(/,' Calculation took up',i8,' seconds wall clock time')") nwalltime2-nwalltime1

end subroutine



!Get value for all the functions used in LSB's DFRT 2.0 project
!Slot of "valarr" array is shown below. rho and gradient of rho are individually returned
!slot=2: |der_rho|/rho^(4/3)
!slot=3: der2rho/rho^(5/3)
!slot=4: (tau-t_w)/t_TF
!slot=5: Xi part of SEDD
!slot=6: Theta part of DORI
!slot=7: Spin density
!slot=8: (tau-t_w)/t_w
subroutine valaryyLSB(x,y,z,valarr,rho,gradrho)
use defvar
use functions
implicit real*8 (a-h,o-z)
integer,parameter :: nfunc=8
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),gradrho(3),hess(3,3),valarr(nfunc)
real*8 gradrhoa(3),gradrhob(3),MOoccnow
real*8 :: Fc=2.871234000D0 ! Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 ! Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194

call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
rho=0D0
gradrho=0D0
do i=1,nmo
	rho=rho+MOocc(i)*wfnval(i)**2
	gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do

gradrho=2*gradrho
rhogrdn=dsqrt(sum(gradrho**2))
dersqr=sum(gradrho**2)
hess(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
hess(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
hess(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
hess(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
hess(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
hess(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
hess(2,1)=hess(1,2)
hess(3,2)=hess(2,3)
hess(3,1)=hess(1,3)
der2rho=hess(1,1)+hess(2,2)+hess(3,3)

!|der_rho|/rho^(4/3)
valarr(2)=rhogrdn/rho**(4D0/3D0)

!der2rho/rho^(5/3)
valarr(3)=der2rho/rho**(5D0/3D0)

!4: (tau-t_w)/t_TF
!8: (tau-t_w)/t_w
D=0D0
rhoa=0D0
rhob=0D0
gradrhoa=0D0
gradrhob=0D0
!Calculate actual kinetic energy (D), Thomas-Fermi kinetic energy (Dh), Weizsacker kinetic energy (t_w)
if (wfntype==0.or.wfntype==3) then !spin-unpolarized case
	do i=1,nmo
		D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Calculate actual kinetic term
	end do
	D=D/2D0
	Dh=Fc*rho**(5D0/3D0) !Thomas-Fermi uniform electron gas kinetic energy
	if (iKEDsel/=0) D=KED(x,y,z,iKEDsel) !Special case proposed by LSB, use other KED instead of exact KED
    t_w=0
	if (rho/=0D0) t_w=sum(gradrho(:)**2)/rho/8D0
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !spin-polarized case
	do i=1,nmo
		MOoccnow=MOocc(i)
		if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
		if (MOtype(i)==1.or.MOtype(i)==0) then
			rhoa=rhoa+MOoccnow*wfnval(i)**2
			gradrhoa(:)=gradrhoa(:)+2D0*MOoccnow*wfnval(i)*wfnderv(:,i)
		end if
		if (MOtype(i)==2.or.MOtype(i)==0) then
			rhob=rhob+MOoccnow*wfnval(i)**2
			gradrhob(:)=gradrhob(:)+2D0*MOoccnow*wfnval(i)*wfnderv(:,i)
		end if
		D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Calculate actual kinetic term
	end do
	D=D/2D0
	Dh=Fc_pol*(rhoa**(5D0/3D0)+rhob**(5D0/3D0))
	if (iKEDsel/=0) D=KED(x,y,z,iKEDsel) !Special case proposed by LSB, use other KED instead of exact KED
    t_w_a=0
	if (rhoa/=0D0) t_w_a=sum(gradrhoa(:)**2)/rhoa/8
    t_w_b=0
	if (rhob/=0D0) t_w_b=sum(gradrhob(:)**2)/rhob/8
    t_w=t_w_a+t_w_b
end if
valarr(4)=(D-t_w)/Dh
valarr(8)=(D-t_w)/t_w
	
!Xi part of SEDD
tmp1_1=rho*(gradrho(1)*hess(1,1)+gradrho(2)*hess(1,2)+gradrho(3)*hess(1,3))
tmp1_2=gradrho(1)*dersqr
tmp2_1=rho*(gradrho(1)*hess(1,2)+gradrho(2)*hess(2,2)+gradrho(3)*hess(2,3))
tmp2_2=gradrho(2)*dersqr
tmp3_1=rho*(gradrho(1)*hess(1,3)+gradrho(2)*hess(2,3)+gradrho(3)*hess(3,3))
tmp3_2=gradrho(3)*dersqr
valarr(5)=4/rho**8*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )

!Theta part of DORI
valarr(6)=4/dersqr**3*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )

!Spin density
rhoa=0D0
rhob=0D0
do i=1,nmo
	if (MOtype(i)==1) then
		rhoa=rhoa+MOocc(i)*wfnval(i)**2
	else if (MOtype(i)==2) then
		rhob=rhob+MOocc(i)*wfnval(i)**2
	else if (MOtype(i)==0) then
		rhoa=rhoa+MOocc(i)/2D0*wfnval(i)**2
		rhob=rhob+MOocc(i)/2D0*wfnval(i)**2
	end if
end do
valarr(7)=rhoa-rhob

end subroutine


!Get value and norm of numerical gradient for all the functions used in LSB's DFRT 2.0 project
!The definition of slot is identical to valaryyLSB
subroutine valgradarrLSB(x,y,z,valarr,grdnarr,rho,gradrho)
use defvar
implicit real*8 (a-h,o-z)
integer,parameter :: nfunc=8
real*8 x,y,z,rho,gradrho(3),valarr(nfunc),grdnarr(nfunc),tmparr(3)
real*8 xadd(nfunc),xmin(nfunc),yadd(nfunc),ymin(nfunc),zadd(nfunc),zmin(nfunc),grdxarr(nfunc),grdyarr(nfunc),grdzarr(nfunc)
diff=5D-4
denom=2D0*diff
call valaryyLSB(x,y,z,valarr,rho,gradrho)
call valaryyLSB(x+diff,y,z,xadd,tmp,tmparr)
call valaryyLSB(x-diff,y,z,xmin,tmp,tmparr)
call valaryyLSB(x,y+diff,z,yadd,tmp,tmparr)
call valaryyLSB(x,y-diff,z,ymin,tmp,tmparr)
call valaryyLSB(x,y,z+diff,zadd,tmp,tmparr)
call valaryyLSB(x,y,z-diff,zmin,tmp,tmparr)
grdxarr=(xadd-xmin)/denom
grdyarr=(yadd-ymin)/denom
grdzarr=(zadd-zmin)/denom
do ifunc=1,nfunc
	grdnarr(ifunc)=dsqrt(grdxarr(ifunc)**2+grdyarr(ifunc)**2+grdzarr(ifunc)**2)
end do
end subroutine



!***** This code is fully useless now, because integrating a function in fuzzy analysis module has already used molecular grid by default
!!------ Integrate real space function in Hirshfeld space with molecular grid (i.e. the grid is the same as integating over the whole space)
! The grid used in function 1 of fuzzy space analysis is only the grid centered at the atom to be studied, this lead to inaccurate integration result
! when the integrand varies fast at the tail region of the Hirshfeld atom (commonly close to the other nuclei). In this case we must use the molecular grid.
! Because incorporate molecular grid integration into "intatomspace" routine will break the structure of the routine, I decide to write this new routine
! dedicated to this purpose. This routine is mainly used in Shubin's study.
subroutine intHirsh_molgrid
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
real*8 funcval(radpot*sphpot),beckeweigrid(radpot*sphpot),atmdens(radpot*sphpot,ncenter),atmintval(ncenter) !Integration value of each atom
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)

write(*,*) "Select the real space function"
call selfunc_interface(1,ifunc)

call setpromol
call gen1cintgrid(gridatmorg,iradcut)

!funcval: real space function at all grid of current center
!atmdens: Free-state density of every atom at all grid of current center
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot

call walltime(iwalltime1)
atmintval=0
do iatm=1,ncenter !Cycle each atom
	write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	
	!Calculate real space function value
	!$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
	do i=1+iradcut*sphpot,radpot*sphpot
		funcval(i)=calcfuncall(ifunc,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
	end do
	!$OMP end parallel do
	
	!Calculate Becke weight
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	
	!Calculate atom densities for evaluating Hirshfeld weight later
	do jatm=1,ncenter_org
		call dealloall(0)
		call readwfn(custommapname(jatm),1)
		!$OMP parallel do shared(atmdens) private(ipt) num_threads(nthreads)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			atmdens(ipt,jatm)=fdens(gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
		end do
		!$OMP end parallel do
	end do
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve to the first loaded file(whole molecule) to calc real rho again
	
	!Generate Hirshfeld weight and integrate the function
	do i=1+iradcut*sphpot,radpot*sphpot
		promol=sum(atmdens(i,:))
		do jatm=1,ncenter
			if (promol/=0) Hirshwei=atmdens(i,jatm)/promol !Hirshfeld weight of jatm at i point
			atmintval(jatm)=atmintval(jatm)+funcval(i)*Hirshwei*beckeweigrid(i)*gridatmorg(i)%value
		end do
	end do
	
end do

call walltime(iwalltime2)
write(*,"(' Calculation took up wall clock time',i10,'s',/)") iwalltime2-iwalltime1

do iatm=1,ncenter
	write(*,"(' Atom',i6,'(',a2,'):',f20.8)") iatm,a(iatm)%name,atmintval(iatm)
end do
write(*,"(' Total:',f20.8,/)") sum(atmintval(:)) 
end subroutine







!!------- Integrate AIM basins using mixed atomic-center and uniform grids for DFRT 2.0 project
subroutine integratebasinmix_LSB
use defvar
use util
use functions
use basinintmod
use topo
implicit real*8 (a-h,o-z)
real*8 trustrad(numrealatt),grad(3),hess(3,3),k1(3),k2(3),k3(3),k4(3),xarr(nx),yarr(ny),zarr(nz)
real*8,allocatable :: potx(:),poty(:),potz(:),potw(:)
type(content),allocatable :: gridatt(:) !Record x,y,z,weight of grids in trust radius
integer att2atm(numrealatt) !The attractor corresponds to which atom. If =0, means this is a NNA
integer walltime1,walltime2,radpotAIM,sphpotAIM
integer,parameter :: nfunc=8,nquant=7 !The number of real space function, the number of quantities to be calculated
real*8 arrtmp(nfunc),arrtmp2(nfunc),gradrho(3)
character(len=40) functionname(0:nfunc),quantityname(0:nquant)
real*8 atmcontri(ncenter,0:nfunc,0:nquant)
real*8,allocatable :: promol(:),promolgrad(:,:),atomdens(:),atomgrad(:,:),funcval(:,:),funcgrdn(:,:),funcref(:,:)
integer,allocatable :: grdidx(:),grdidy(:),grdidz(:) !Record x,y,z index of all the grids belonging to an attractor

itype=2 !Exact refinement
radpotAIM=200
nbeckeiter=8
call setpromol
expcutoff=1
functionname(0)=" rho"
functionname(1)=" rho/rho0"
functionname(2)=" |der_rho|/rho^(4/3)"
functionname(3)=" der2rho/rho^(5/3)"
functionname(4)=" (tau-t_w)/t_TF"
functionname(5)=" Xi part of SEDD"
functionname(6)=" Theta part of DORI"
functionname(7)=" Spin density"
functionname(8)=" (tau-t_w)/t_w"
quantityname(0)=" Function itself" !Special
quantityname(1)=" Shannon entropy"
quantityname(2)=" Fisher information"
quantityname(3)=" Onicescu information energy of order 2"
quantityname(4)=" Onicescu information energy of order 3"
quantityname(5)=" Information gain"
quantityname(6)=" Relative Renyi entropy of orders 2"
quantityname(7)=" Relative Renyi entropy of orders 3"

invfunc=0
write(*,*) "If replacing real space functions with their reciprocals?  0=No  1=Yes"
read(*,*) invfunc 

call walltime(walltime1)

numcp=0
att2atm=0
atmcontri=0
!Determine trust radius and then integrate in the trust sphere
write(*,*) "Integrating in trust sphere..."

do iatt=1,numrealatt !Cycle each attractors

	!Find one-to-one correspondence between NCP and atom, then the NCP position is defined as nuclear position
	do iatm=1,ncenter
		disttest=dsqrt( (realattxyz(1,iatt)-a(iatm)%x)**2+(realattxyz(2,iatt)-a(iatm)%y)**2+(realattxyz(3,iatt)-a(iatm)%z)**2 )
		if (disttest<0.3D0) then
			att2atm(iatt)=iatm
			write(*,"(' Attractor',i6,' corresponds to atom',i6,' (',a,')')") iatt,iatm,a(iatm)%name
			numcpold=numcp
			!Refine the crude position of attractor by exact newton method
			call findcp(a(iatm)%x,a(iatm)%y,a(iatm)%z,1)
			if (numcp==numcpold) then
				write(*,*) "Note: Unable to locate exact CP position! Use nuclear position"
				numcp=numcp+1
				CPpos(1,numcp)=a(iatm)%x
				CPpos(2,numcp)=a(iatm)%y
				CPpos(3,numcp)=a(iatm)%z
			end if
			exit
		end if
	end do
	if (att2atm(iatt)==0) then
		write(*,"(a,i6,a)") " Warning: Unable to determine the attractor",iatt," belongs to which atom!"
		write(*,*) "Non-nuclear attractor is not supported."
		write(*,*) "Press ENTER button to exit"
		read(*,*)
		return
	end if
	
	!Determine trust radius and set integration points and weight
	parm=1
	isettrustrad=0
	nintgrid=0 !Then number of integration grids within trust radius
	if (allocated(gridatt)) deallocate(gridatt) !Used to record information of grids in trust sphere of this attractor
	allocate(gridatt(radpotAIM*500))
	do ish=1,radpotAIM !Cycle each radial shell. Radius distance is from near to far
		if (isettrustrad==1) exit !The trust radius has been finally determined in last shell cycle
		!Becke, namely the second-kind Gauss-Chebyshev
		itmp=radpotAIM+1-ish !Invert ish to make radr from near to far
		radx=cos(itmp*pi/(radpotAIM+1D0))
		radr=(1+radx)/(1-radx)*parm
		radw=2*pi/(radpotAIM+1)*parm**3 *(1+radx)**2.5D0/(1-radx)**3.5D0 *4*pi
		!Set proper Lebedev grid according to shell radius
		radtmp=covr(a(iatm)%index)
		if (radr<0.2D0*radtmp) then
			sphpotAIM=26
		else if (radr<0.5D0*radtmp) then
			sphpotAIM=74
		else if (radr<0.8D0*radtmp) then
			sphpotAIM=146
		else
			sphpotAIM=194
		end if
		if (allocated(potx)) deallocate(potx,poty,potz,potw)
		allocate(potx(sphpotAIM),poty(sphpotAIM),potz(sphpotAIM),potw(sphpotAIM))
		call Lebedevgen(sphpotAIM,potx,poty,potz,potw)
		!Combine radial point and weights with angular part, and make them centered at current attractor
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%x=radr*potx+CPpos(1,numcp)
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%y=radr*poty+CPpos(2,numcp)
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%z=radr*potz+CPpos(3,numcp)
		gridatt( nintgrid+1:nintgrid+sphpotAIM )%value=radw*potw
		!Find trust radius for present attractor
		angmax=0
		radrinit=0.15D0
		if (a(iatm)%index>2) radrinit=0.5D0
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
					isettrustrad=1
					exit
				end if
			end do
			if (isettrustrad==0) trustrad(iatt)=radr !Passed this shell and temporarily set the radius as trust radius. Continue to enlarge the trust radius, until reached angmax>45 degree
		end if
		nintgrid=nintgrid+sphpotAIM
	end do
	if (isettrustrad==0) trustrad(iatt)=20
	write(*,"(' The trust radius of attractor',i6,' is',f10.3,' Bohr',/)") iatt,trustrad(iatt)
	
	allocate(promol(nintgrid),promolgrad(nintgrid,3),atomdens(nintgrid),atomgrad(nintgrid,3),funcval(nintgrid,0:nfunc),funcgrdn(nintgrid,0:nfunc),funcref(nintgrid,0:nfunc))
	
	!! Calculate function value of reference state, including promolecular density and density of present atom, function index:
    !index=0: rho
    !index=1: rho/rho0, where rho0 is promolecular density
	!index=2: |der_rho|/rho^(4/3)
	!index=3: der2rho/rho^(5/3)
	!index=4: (tau-t_w)/t_TF
	!index=5: Part of SEDD
	!index=6: Part of DORI
	!index=7: Spin density
	!index=8: (tau-t_w)/t_w
	funcref=0D0
	promol=0D0
	promolgrad=0D0
	do jatm=1,ncenter_org
		call dealloall(0)
		call readwfn(custommapname(jatm),1)
		!$OMP parallel do shared(atomdens,atomgrad,funcref) private(ipt,arrtmp,rho,gradrho) num_threads(nthreads)
		do ipt=1,nintgrid
			call valaryyLSB(gridatt(ipt)%x,gridatt(ipt)%y,gridatt(ipt)%z,arrtmp,rho,gradrho)
			atomdens(ipt)=rho
			atomgrad(ipt,:)=gradrho
			funcref(ipt,0)=funcref(ipt,0)+rho
			funcref(ipt,1)=funcref(ipt,1)+1D0 !i.e. rho/rho0, the value is unity if present system is a single atom
			funcref(ipt,2:nfunc)=funcref(ipt,2:nfunc)+arrtmp(2:nfunc)
		end do
		!$OMP end parallel do
		promol=promol+atomdens(:)
		promolgrad=promolgrad+atomgrad(:,:)
	end do
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule)
	
    !! Calculate function value and gradient norm for present molecule, then store them in funcval/funcgrdn
	!$OMP parallel do shared(funcval,funcgrdn) private(ipt,arrtmp,arrtmp2,rho,gradrho,idir) num_threads(nthreads)
	do ipt=1,nintgrid
		call valgradarrLSB(gridatt(ipt)%x,gridatt(ipt)%y,gridatt(ipt)%z,arrtmp,arrtmp2,rho,gradrho)
		funcval(ipt,0)=rho
		funcgrdn(ipt,0)=dsqrt(sum(gradrho**2))
		!Note: see word document for explicit expression of |der(rho/rho0)|
		funcval(ipt,1)=rho/promol(ipt)
		funcgrdn(ipt,1)=0D0
		do idir=1,3
			funcgrdn(ipt,1)=funcgrdn(ipt,1)+( gradrho(idir)/promol(ipt)-rho/promol(ipt)**2*promolgrad(ipt,idir) )**2
		end do
		funcgrdn(ipt,1)=dsqrt(funcgrdn(ipt,1))
		funcval(ipt,2:nfunc)=arrtmp(2:nfunc)
		funcgrdn(ipt,2:nfunc)=arrtmp2(2:nfunc)
	end do
	!$OMP end parallel do
	
	!! Integrate the region inside trust radius
	do ipt=1,nintgrid
		rx=gridatt(ipt)%x-CPpos(1,numcp) !The relative distance between current point to corresponding attractor
		ry=gridatt(ipt)%y-CPpos(2,numcp)
		rz=gridatt(ipt)%z-CPpos(3,numcp)
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
		
		weitot=switchwei*gridatt(ipt)%value
		do ifunc=0,nfunc
			if (invfunc==0) then
				fv=funcval(ipt,ifunc)
				fgrdn=funcgrdn(ipt,ifunc)
				fref=funcref(ipt,ifunc)
			else if (invfunc==1) then
				fv=1/funcval(ipt,ifunc)
				fgrdn=1/funcgrdn(ipt,ifunc)
				fref=1/funcref(ipt,ifunc)
			end if
			!Value of function itself
			atmcontri(iatm,ifunc,0)=atmcontri(iatm,ifunc,0)+weitot*fv
			!Shannon entropy
			atmcontri(iatm,ifunc,1)=atmcontri(iatm,ifunc,1)+weitot*( -fv*log(fv) )
			!Fisher information
			atmcontri(iatm,ifunc,2)=atmcontri(iatm,ifunc,2)+weitot*( fgrdn**2 /fv )
			!Onicescu information energy of order 2
			atmcontri(iatm,ifunc,3)=atmcontri(iatm,ifunc,3)+weitot*( fv**2 )
			!Onicescu information energy of order 3
			atmcontri(iatm,ifunc,4)=atmcontri(iatm,ifunc,4)+weitot*( fv**3 )
			!Information gain
			atmcontri(iatm,ifunc,5)=atmcontri(iatm,ifunc,5)+weitot*( fv*log(fv/fref) )
			!Relative Renyi entropy of orders 2
			atmcontri(iatm,ifunc,6)=atmcontri(iatm,ifunc,6)+weitot*( fv**2/fref )
			!Relative Renyi entropy of orders 3
			atmcontri(iatm,ifunc,7)=atmcontri(iatm,ifunc,7)+weitot*( fv**3/fref**2 )
		end do
	end do
	
! 	write(*,"(' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
! 	do ifunc=0,nfunc
! 		write(*,"(/,' Function:',a)") trim(functionname(ifunc))
! 		do iquant=0,nquant !0 corresponds to function itself
! 			write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
! 		end do
! 	end do
	
	deallocate(promol,promolgrad,atomdens,atomgrad,funcval,funcgrdn,funcref)
end do !End cycle attractors

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

!---------------------------------------------!
!--------- Integrating uniform grids ---------!
!---------------------------------------------!

!--------- Exact refinement boundary grids
write(*,*) "Refining boundary grids..."
if (itype==2.or.itype==3) then
	if (itype==3) then !Calculate grid data of gradient of electron density used to linear interpolation to obtain the value at any point
		write(*,*)
		call gengradmat
	end if
	nrk4lim=100
	nrk4gradswitch=40
	hsizeinit=0.25D0
	ifinish=0
	!$OMP PARALLEL do private(ix,iy,iz,iatt,rnowx,rnowy,rnowz,rx,ry,rz,tmpval,tmpval2,tmpval3,&
	!$OMP rnowxtmp,rnowytmp,rnowztmp,orgxref,orgyref,orgzref,dxref,dyref,dzref,ixref,iyref,izref,nrefine,ndiv,&
	!$OMP k1,k2,k3,k4,dens,denshold,grad,hess,iattref,xtmp,ytmp,ztmp,irk4,hsize,ixtest,iytest,iztest,tmpdist) shared(ifinish) NUM_THREADS(nthreads) schedule(DYNAMIC)
	do iz=2,nz-1
		rnowz=zarr(iz)
		do iy=2,ny-1
			rnowy=yarr(iy)
			do ix=2,nx-1
				rnowx=xarr(ix)
				if (.not.interbasgrid(ix,iy,iz)) cycle
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
									do ixtest=2,nx-1
										tmpdist=abs(xtmp-xarr(ixtest))
										if (tmpdist<dx/2D0) exit
									end do
									do iytest=2,ny-1
										tmpdist=abs(ytmp-yarr(iytest))
										if (tmpdist<dy/2D0) exit
									end do
									do iztest=2,nz-1
										tmpdist=abs(ztmp-zarr(iztest))
										if (tmpdist<dz/2D0) exit
									end do
									iattref=gridbas(ixtest,iytest,iztest)
									do imove=1,26
										if ( gridbas(ixtest+vec26x(imove),iytest+vec26y(imove),iztest+vec26z(imove))/=iattref ) exit
									end do
									if (imove==27) exit !Successfully passed neighbour test
								end do cycrk4
								if (irk4==nrk4lim+1) then !Didn't enter trust radius or didn't approach a grid who and whose neighbour have the same attribution
! 									write(*,*) "Warning: Exceeded the step limit of steepest ascent process!"
									iattref=gridbas(ix,iy,iz) !Use its original attribution
								end if
							end if
							gridbas(ix,iy,iz)=iattref !Update attribution of boundary grids
						end do !End refine grid
					end do
				end do
				
			end do !End cycle ix grid
		end do
		ifinish=ifinish+1
	end do
	!$OMP end PARALLEL do
	call detectinterbasgrd(6)
	write(*,*) "Basin boundary has been updated"
end if

!!---------- Calculate part of quantities contributed from uniform grid, and gain final result
write(*,*)
write(*,*) "Calculating information at uniform grid..."
do iatt=1,numrealatt !Cycle each attractors
	write(*,"(' Processing basin',i5,' /',i5)") iatt,numrealatt
	iatm=att2atm(iatt)
	nintgrid=count(gridbas(2:nx-1,2:ny-1,2:nz-1)==iatt)
	allocate(promol(nintgrid),promolgrad(nintgrid,3),atomdens(nintgrid),atomgrad(nintgrid,3),funcval(nintgrid,0:nfunc),funcgrdn(nintgrid,0:nfunc),funcref(nintgrid,0:nfunc))
	allocate(grdidx(nintgrid),grdidy(nintgrid),grdidz(nintgrid))
	itmp=0
	do iz=2,nz-1
		do iy=2,ny-1
			do ix=2,nx-1
				if (gridbas(ix,iy,iz)==iatt) then !Find grids attributing to this attractor
					itmp=itmp+1
					grdidx(itmp)=ix
					grdidy(itmp)=iy
					grdidz(itmp)=iz
				end if
			end do
		end do
	end do
	
	!Calculate information for promolecular state
	funcref=0D0
	promol=0D0
	promolgrad=0D0
	do jatm=1,ncenter_org
		call dealloall(0)
		call readwfn(custommapname(jatm),1)
		!$OMP parallel do shared(atomdens,atomgrad,funcref) private(ipt,ptx,pty,ptz,arrtmp,rho,gradrho) num_threads(nthreads)
		do ipt=1,nintgrid
			ptx=xarr(grdidx(ipt))
			pty=yarr(grdidy(ipt))
			ptz=zarr(grdidz(ipt))
			call valaryyLSB(ptx,pty,ptz,arrtmp,rho,gradrho)
			atomdens(ipt)=rho
			atomgrad(ipt,:)=gradrho
			funcref(ipt,0)=funcref(ipt,0)+rho
			funcref(ipt,1)=funcref(ipt,1)+1D0 !i.e. rho/rho0, the value is unity if present system is a single atom
			funcref(ipt,2:nfunc)=funcref(ipt,2:nfunc)+arrtmp(2:nfunc)
		end do
		!$OMP end parallel do
		promol=promol+atomdens(:)
		promolgrad=promolgrad+atomgrad(:,:)
	end do
	call dealloall(0)
	call readinfile(firstfilename,1) !Retrieve the firstly loaded file(whole molecule)
	
	!Calculate information for present actual state
	!$OMP PARALLEL do shared(funcval,funcgrdn) private(ipt,ptx,pty,ptz,arrtmp,arrtmp2,rho,gradrho,idir) NUM_THREADS(nthreads) schedule(DYNAMIC)
	do ipt=1,nintgrid
		ptx=xarr(grdidx(ipt))
		pty=yarr(grdidy(ipt))
		ptz=zarr(grdidz(ipt))
		call valgradarrLSB(ptx,pty,ptz,arrtmp,arrtmp2,rho,gradrho)
		funcval(ipt,0)=rho
		funcgrdn(ipt,0)=dsqrt(sum(gradrho**2))
		!Note: See Word document for explicit expression of |der(rho/rho0)|
		funcval(ipt,1)=rho/promol(ipt)
		funcgrdn(ipt,1)=0D0
		do idir=1,3
			funcgrdn(ipt,1)=funcgrdn(ipt,1)+( gradrho(idir)/promol(ipt)-rho/promol(ipt)**2*promolgrad(ipt,idir) )**2
		end do
		funcgrdn(ipt,1)=dsqrt(funcgrdn(ipt,1))
		funcval(ipt,2:nfunc)=arrtmp(2:nfunc)
		funcgrdn(ipt,2:nfunc)=arrtmp2(2:nfunc)
	end do
	!$OMP end PARALLEL do
	
	!Accumulate uniform grid contribution to quantities
	do ipt=1,nintgrid
		ptx=xarr(grdidx(ipt))
		pty=yarr(grdidy(ipt))
		ptz=zarr(grdidz(ipt))
		!Calculate switching function at current grid
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
		
		weitot=switchwei*dvol
		do ifunc=0,nfunc
			if (invfunc==0) then
				fv=funcval(ipt,ifunc)
				fgrdn=funcgrdn(ipt,ifunc)
				fref=funcref(ipt,ifunc)
			else if (invfunc==1) then
				fv=1/funcval(ipt,ifunc)
				fgrdn=1/funcgrdn(ipt,ifunc)
				fref=1/funcref(ipt,ifunc)
			end if
			!Value of function itself
			atmcontri(iatm,ifunc,0)=atmcontri(iatm,ifunc,0)+weitot*fv
			!Shannon entropy
			atmcontri(iatm,ifunc,1)=atmcontri(iatm,ifunc,1)+weitot*( -fv*log(fv) )
			!Fisher information
			atmcontri(iatm,ifunc,2)=atmcontri(iatm,ifunc,2)+weitot*( fgrdn**2 /fv )
			!Onicescu information energy of order 2
			atmcontri(iatm,ifunc,3)=atmcontri(iatm,ifunc,3)+weitot*( fv**2 )
			!Onicescu information energy of order 3
			atmcontri(iatm,ifunc,4)=atmcontri(iatm,ifunc,4)+weitot*( fv**3 )
			!Information gain
			atmcontri(iatm,ifunc,5)=atmcontri(iatm,ifunc,5)+weitot*( fv*log(fv/fref) )
			!Relative Renyi entropy of orders 2
			atmcontri(iatm,ifunc,6)=atmcontri(iatm,ifunc,6)+weitot*( fv**2/fref )
			!Relative Renyi entropy of orders 3
			atmcontri(iatm,ifunc,7)=atmcontri(iatm,ifunc,7)+weitot*( fv**3/fref**2 )
		end do
	end do
	
	deallocate(promol,promolgrad,atomdens,atomgrad,funcval,funcgrdn,funcref,grdidx,grdidy,grdidz)
end do !End cycling attractor

do iatm=1,ncenter
	write(*,"(/,' ========== Contribution of atom',i4,a,':')") iatm,a(iatm)%name
	do ifunc=0,nfunc
		write(*,"(/,' Function:',a)") trim(functionname(ifunc))
		do iquant=0,nquant !0 corresponds to function itself
			write(*,"(a,':',1E16.8)") quantityname(iquant),atmcontri(iatm,ifunc,iquant)
		end do
	end do
end do
write(*,*)
write(*,*) "========================================="
write(*,*) "             Overall result"
write(*,*) "========================================="
do ifunc=0,nfunc
	write(*,"(/,' Function:',a)") trim(functionname(ifunc))
	do iquant=0,nquant
		write(*,"(a,':',1E16.8)") quantityname(iquant),sum(atmcontri(:,ifunc,iquant))
	end do
end do

call walltime(walltime2)
write(*,"(' Integrating basins took up wall clock time',i10,' s')") walltime2-walltime1

end subroutine







!!!------------ Integrate various forms of kinetic energy density (KED) over the whole space
subroutine intKED
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 intval(nKEDmax),funcval(radpot*sphpot,nKEDmax),beckeweigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)
intval=0
do iatm=1,ncenter
	write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	!$OMP parallel do shared(funcval) private(i) num_threads(nthreads)
	do i=1+iradcut*sphpot,radpot*sphpot
		call KEDall(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,funcval(i,:))
	end do
	!$OMP end parallel do
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	do i=1+iradcut*sphpot,radpot*sphpot
		intval=intval+funcval(i,:)*gridatmorg(i)%value*beckeweigrid(i)
	end do
end do

write(*,"(' KED',i3,'(Lagrangian):  ',f20.10)") 2,intval(2)
write(*,"(' KED',i3,'(Thomas-Fermi):',f20.10)") 3,intval(3)
write(*,"(' KED',i3,'(Weizsacker):  ',f20.10)") 4,intval(4)
write(*,"(' KED',i3,'(GEA2):        ',f20.10)") 5,intval(5)
write(*,"(' KED',i3,'(TF5W):        ',f20.10)") 6,intval(6)
write(*,"(' KED',i3,'(TFvW):        ',f20.10)") 7,intval(7)
write(*,"(' KED',i3,'(TF9W):        ',f20.10)") 8,intval(8)
write(*,"(' KED',i3,'(TF-N):        ',f20.10)") 9,intval(9)
write(*,"(' KED',i3,'(Pearson):     ',f20.10)") 10,intval(10)
write(*,"(' KED',i3,'(DK Pade):     ',f20.10)") 11,intval(11)
write(*,"(' KED',i3,'(LLP):         ',f20.10)") 12,intval(12)
write(*,"(' KED',i3,'(OL1):         ',f20.10)") 13,intval(13)
write(*,"(' KED',i3,'(OL2):         ',f20.10)") 14,intval(14)
write(*,"(' KED',i3,'(Thakkar):     ',f20.10)") 15,intval(15)
write(*,"(' KED',i3,'(Becke86A):    ',f20.10)") 16,intval(16)
write(*,"(' KED',i3,'(Becke86B):    ',f20.10)") 17,intval(17)
write(*,"(' KED',i3,'(DK87):        ',f20.10)") 18,intval(18)
write(*,"(' KED',i3,'(PW86):        ',f20.10)") 19,intval(19)
write(*,"(' KED',i3,'(PW91):        ',f20.10)") 20,intval(20)
write(*,"(' KED',i3,'(LG94):        ',f20.10)") 21,intval(21)
write(*,"(' KED',i3,'(ABSP):        ',f20.10)") 22,intval(22)
write(*,"(' KED',i3,'(GR):          ',f20.10)") 23,intval(23)
write(*,"(' KED',i3,'(GEA4):        ',f20.10)") 24,intval(24)
write(*,*) "Press ENTER button to continue"
read(*,*)
end subroutine



!!------- Calculate some quantities involved in Shubin's project in a plane
!itype=1: Calculate the sum of atomic relative Shannon entropy (namely total relative Shannon entropy). sum{rho(A)*ln[rho(A)/rho0(A)]}
!itype=2: Calculate the sum of x=[rhoA-rho0A]/rhoA
!itype=3: Calculate the difference between total relative Shannon entropy and deformation density
!itype=4: Calculate 2nd relative Onicescu information sum{[rho(A)]^2/rho0(A)}
!itype=5: Calculate 3rd relative Onicescu information sum{[rho(A)]^3/[rho0(A)]^2}/2
subroutine genentroplane(itype)
use defvar
use functions
implicit real*8 (a-h,o-z)
integer itype
real*8 planeprodens(ngridnum1,ngridnum2),planedens(ngridnum1,ngridnum2)

if (allocated(planemat)) deallocate(planemat)
allocate(planemat(ngridnum1,ngridnum2))
planeprodens=0D0
planemat=0D0

!Calculate molecular density in the plane and store it to planedens
!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz) shared(planedens) schedule(dynamic) NUM_THREADS(nthreads)
do i=1,ngridnum1
	do j=1,ngridnum2
        call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
		planedens(i,j)=fdens(rnowx,rnowy,rnowz)
	end do
end do
!$OMP END PARALLEL DO

!Calculate promolecular density in the plane and store it to planeprodens
do jatm=1,ncenter_org
	call dealloall(0)
	call readwfn(custommapname(jatm),1)
	!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz) shared(planeprodens) schedule(dynamic) NUM_THREADS(nthreads)
	do i=1,ngridnum1
		do j=1,ngridnum2
            call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
			planeprodens(i,j)=planeprodens(i,j)+fdens(rnowx,rnowy,rnowz)
		end do
	end do
	!$OMP END PARALLEL DO
end do

!Calculate Hirshfeld weight, relative Shannon entropy and x=[rhoA-rho0A]/rhoA for each atom in the plane and accumulate them to planemat
do jatm=1,ncenter_org !Cycle each atom, calculate its contribution in the plane
	call dealloall(0)
	call readwfn(custommapname(jatm),1)
	!$OMP PARALLEL DO private(i,j,rnowx,rnowy,rnowz,rho0A,rhoA,tmpval) shared(planemat) schedule(dynamic) NUM_THREADS(nthreads)
	do i=1,ngridnum1
		do j=1,ngridnum2
            call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
			rho0A=fdens(rnowx,rnowy,rnowz)
			rhoA=planedens(i,j)*rho0A/planeprodens(i,j)
			if (itype==1.or.itype==3) then
				tmpval=rhoA*log(rhoA/rho0A) !Relative Shannon entropy
			else if (itype==2) then
				tmpval=(rhoA-rho0A)/rhoA !x=[rhoA-rho0A]/rhoA
			else if (itype==4) then
				tmpval=rhoA**2/rho0A
			else if (itype==5) then
				tmpval=(rhoA**3/rho0A**2)/2D0
            end if
			planemat(i,j)=planemat(i,j)+tmpval
		end do
	end do
	!$OMP END PARALLEL DO
end do
call dealloall(0)
call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule)
if (itype==3) planemat=planemat-(planedens-planeprodens) !Diff between total relative Shannon entropy and deformation density
end subroutine




!!------- Calculate some quantities involved in Shubin's project as grid data
!Definition is the same as subroutine genentroplane
subroutine genentrocub(itype)
use defvar
use functions
implicit real*8 (a-h,o-z)
integer itype
real*8 cubprodens(nx,ny,nz),cubdens(nx,ny,nz)

call setpromol

if (allocated(cubmat)) deallocate(cubmat)
allocate(cubmat(nx,ny,nz))
cubprodens=0D0
cubmat=0D0

!Calculate molecular density grid data and store it to cubdens
write(*,*) "Calculating real electron density grid data..."
!$OMP PARALLEL DO SHARED(cubdens) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
do k=1,nz
	do j=1,ny
		do i=1,nx
			call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
			cubdens(i,j,k)=fdens(tmpx,tmpy,tmpz)
		end do
	end do
end do
!$OMP END PARALLEL DO

write(*,*) "Calculating promolecular density grid data..."
do jatm=1,ncenter_org
	call dealloall(0)
	call readwfn(custommapname(jatm),1)
	!$OMP PARALLEL DO SHARED(cubprodens) PRIVATE(i,j,k,tmpx,tmpy,tmpz) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                cubprodens(i,j,k)=cubprodens(i,j,k)+fdens(tmpx,tmpy,tmpz)
			end do
		end do
	end do
	!$OMP END PARALLEL DO
    call showprog(jatm,ncenter_org)
end do

write(*,*) "Calculating promolecular grid data of information-theoretic quantities..."
do jatm=1,ncenter_org !Cycle each atom, calculate its contribution
	call dealloall(0)
	call readwfn(custommapname(jatm),1)
	!$OMP PARALLEL DO SHARED(cubmat) PRIVATE(i,j,k,tmpx,tmpy,tmpz,rho0A,rhoA,tmpval) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                rho0A=fdens(tmpx,tmpy,tmpz)
				rhoA=cubdens(i,j,k)*rho0A/cubprodens(i,j,k)
				if (itype==1.or.itype==3) then
					tmpval=rhoA*log(rhoA/rho0A) !Relative Shannon entropy
				else if (itype==2) then
					tmpval=(rhoA-rho0A)/rhoA !x=[rhoA-rho0A]/rhoA
				else if (itype==4) then
					tmpval=rhoA**2/rho0A
				else if (itype==5) then
					tmpval=(rhoA**3/rho0A**2)/2D0
				end if
				cubmat(i,j,k)=cubmat(i,j,k)+tmpval
			end do
		end do
	end do
	!$OMP END PARALLEL DO
    call showprog(jatm,ncenter_org)
end do
call dealloall(0)
call readinfile(firstfilename,1) !Retrieve the first loaded file(whole molecule)
if (itype==3) cubmat=cubmat-(cubdens-cubprodens) !Diff between total relative Shannon entropy and deformation density
end subroutine




!!---------- Calculate g1, g2, g3 terms along a line, all of them rely on promolecular density calculated based on atom .wfn files
subroutine g1g2g3line(orgx1D,orgy1D,orgz1D,transx,transy,transz)
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 orgx1D,orgy1D,orgz1D,transx,transy,transz
real*8 rho(num1Dpoints),derrho(3,num1Dpoints),hessrho(3,3,num1Dpoints),rho0(num1Dpoints),derrho0(3,num1Dpoints),hessrho0(3,3,num1Dpoints)
real*8 ptxyz(3,num1Dpoints),tmparr(3),tmpmat(3,3)

do ipt=1,num1Dpoints
	curvex(ipt)=ipt*dsqrt(transx**2+transy**2+transz**2)
	ptxyz(1,ipt)=orgx1D+(ipt-1)*transx
	ptxyz(2,ipt)=orgy1D+(ipt-1)*transy
	ptxyz(3,ipt)=orgz1D+(ipt-1)*transz
end do

!Calculate molecular density
!$OMP parallel do shared(rho,derrho,hessrho) private(ipt) num_threads(nthreads)
do ipt=1,num1Dpoints
    call calchessmat_dens(2,ptxyz(1,ipt),ptxyz(2,ipt),ptxyz(3,ipt),rho(ipt),derrho(:,ipt),hessrho(:,:,ipt))
end do
!$OMP end parallel do

call setpromol
!Calculate promolecular density
rho0=0
derrho0=0
hessrho0=0
do ipro=1,ncustommap
    filename=custommapname(ipro)
	call dealloall(0)
	write(*,"(' Loading: ',a)") trim(filename)
	call readinfile(filename,1)
    !$OMP parallel do shared(rho0,derrho0,hessrho0) private(ipt,tmprho,tmparr,tmpmat) num_threads(nthreads)
    do ipt=1,num1Dpoints
        call calchessmat_dens(2,ptxyz(1,ipt),ptxyz(2,ipt),ptxyz(3,ipt),tmprho,tmparr,tmpmat)
        rho0(ipt)=rho0(ipt)+tmprho
        derrho0(:,ipt)=derrho0(:,ipt)+tmparr(:)
        hessrho0(:,:,ipt)=hessrho0(:,:,ipt)+tmpmat(:,:)
    end do
    !$OMP end parallel do
end do

call dealloall(0)
write(*,"(' Reloading: ',a)") trim(firstfilename)
call readinfile(firstfilename,1)

!Calculate function values
do ipt=1,num1Dpoints
    rholapl=hessrho(1,1,ipt)+hessrho(2,2,ipt)+hessrho(3,3,ipt)
    rholapl0=hessrho0(1,1,ipt)+hessrho0(2,2,ipt)+hessrho0(3,3,ipt)
    if (iuserfunc==57) then
        curvey(ipt)=rholapl*log(rho(ipt)/rho0(ipt))
    else if (iuserfunc==58) then
        curvey(ipt)=rho(ipt)*(rholapl/rho(ipt)-rholapl0/rho0(ipt))
    else if (iuserfunc==59) then
        xtmp=derrho(1,ipt)/rho(ipt)-derrho0(1,ipt)/rho0(ipt)
        ytmp=derrho(2,ipt)/rho(ipt)-derrho0(2,ipt)/rho0(ipt)
        ztmp=derrho(3,ipt)/rho(ipt)-derrho0(3,ipt)/rho0(ipt)
        curvey(ipt)=rho(ipt)*(xtmp**2+ytmp**2+ztmp**2)
    end if
end do
end subroutine



!!---------- Calculate g1, g2, g3 terms in the plane, all of them rely on promolecular density calculated based on atom .wfn files
subroutine g1g2g3plane
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 rho(ngridnum1,ngridnum2),derrho(3,ngridnum1,ngridnum2),hessrho(3,3,ngridnum1,ngridnum2)
real*8 rho0(ngridnum1,ngridnum2),derrho0(3,ngridnum1,ngridnum2),hessrho0(3,3,ngridnum1,ngridnum2)
real*8 tmparr(3),tmpmat(3,3)

!Calculate molecular density
!$OMP parallel do shared(rho,derrho,hessrho) private(i,j,rnowx,rnowy,rnowz) num_threads(nthreads)
do i=1,ngridnum1
	do j=1,ngridnum2
        call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
        call calchessmat_dens(2,rnowx,rnowy,rnowz,rho(i,j),derrho(:,i,j),hessrho(:,:,i,j))
    end do
end do
!$OMP end parallel do

call setpromol
!Calculate promolecular density
rho0=0
derrho0=0
hessrho0=0
do ipro=1,ncustommap
    filename=custommapname(ipro)
	call dealloall(0)
	write(*,"(' Loading: ',a)") trim(filename)
	call readinfile(filename,1)
    !$OMP parallel do shared(rho,derrho,hessrho) private(i,j,rnowx,rnowy,rnowz,tmprho,tmparr,tmpmat) num_threads(nthreads)
    do i=1,ngridnum1
	    do j=1,ngridnum2
            call get2Dgridxyz(i,j,rnowx,rnowy,rnowz)
            call calchessmat_dens(2,rnowx,rnowy,rnowz,tmprho,tmparr,tmpmat)
            rho0(i,j)=rho0(i,j)+tmprho
            derrho0(:,i,j)=derrho0(:,i,j)+tmparr(:)
            hessrho0(:,:,i,j)=hessrho0(:,:,i,j)+tmpmat(:,:)
        end do
    end do
    !$OMP end parallel do
end do

call dealloall(0)
write(*,"(' Reloading: ',a)") trim(firstfilename)
call readinfile(firstfilename,1)

!Calculate function values
do i=1,ngridnum1
	do j=1,ngridnum2
        rholapl=hessrho(1,1,i,j)+hessrho(2,2,i,j)+hessrho(3,3,i,j)
        rholapl0=hessrho0(1,1,i,j)+hessrho0(2,2,i,j)+hessrho0(3,3,i,j)
        if (iuserfunc==57) then
            planemat(i,j)=rholapl*log(rho(i,j)/rho0(i,j))
        else if (iuserfunc==58) then
            planemat(i,j)=rho(i,j)*(rholapl/rho(i,j)-rholapl0/rho0(i,j))
        else if (iuserfunc==59) then
            xtmp=derrho(1,i,j)/rho(i,j)-derrho0(1,i,j)/rho0(i,j)
            ytmp=derrho(2,i,j)/rho(i,j)-derrho0(2,i,j)/rho0(i,j)
            ztmp=derrho(3,i,j)/rho(i,j)-derrho0(3,i,j)/rho0(i,j)
            planemat(i,j)=rho(i,j)*(xtmp**2+ytmp**2+ztmp**2)
        end if
    end do
end do
end subroutine



!!---------- Calculate g1, g2, g3 terms as grid data, all of them rely on promolecular density calculated based on atom .wfn files
subroutine g1g2g3grid
use defvar
use functions
implicit real*8 (a-h,o-z)
real*8 rho(nx,ny,nz),derrho(3,nx,ny,nz),rholapl(nx,ny,nz)
real*8 rho0(nx,ny,nz),derrho0(3,nx,ny,nz),rholapl0(nx,ny,nz)
real*8 tmparr(3),tmpmat(3,3)

write(*,*) "Calculating electron density and derivatives for actual molecule..."
!$OMP PARALLEL DO SHARED(rho,derrho,rholapl) PRIVATE(i,j,k,tmpx,tmpy,tmpz,tmpmat) schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
            call calchessmat_dens(2,tmpx,tmpy,tmpz,rho(i,j,k),derrho(:,i,j,k),tmpmat)
            rholapl(i,j,k)=tmpmat(1,1)+tmpmat(2,2)+tmpmat(3,3)
        end do
    end do
end do
!$OMP end parallel do

write(*,*) "Calculating electron density and derivatives for promolecule..."
call setpromol
!Calculate promolecular density
rho0=0
derrho0=0
rholapl0=0
do ipro=1,ncustommap
    filename=custommapname(ipro)
	call dealloall(0)
	write(*,"(' Loading: ',a)") trim(filename)
	call readinfile(filename,1)
    !$OMP parallel do shared(rho0,derrho0,rholapl0) private(i,j,k,tmpx,tmpy,tmpz,tmprho,tmparr,tmpmat) num_threads(nthreads)
    do k=1,nz
	    do j=1,ny
		    do i=1,nx
                call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
                call calchessmat_dens(2,tmpx,tmpy,tmpz,tmprho,tmparr,tmpmat)
                rho0(i,j,k)=rho0(i,j,k)+tmprho
                derrho0(:,i,j,k)=derrho0(:,i,j,k)+tmparr(:)
                rholapl0(i,j,k)=rholapl0(i,j,k)+(tmpmat(1,1)+tmpmat(2,2)+tmpmat(3,3))
            end do
        end do
    end do
    !$OMP end parallel do
end do

call dealloall(0)
write(*,"(' Reloading: ',a)") trim(firstfilename)
call readinfile(firstfilename,1)

write(*,*) "Calculating final function values..."
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,tmpx,tmpy,tmpz)
            if (iuserfunc==57) then
                cubmat(i,j,k)=rholapl(i,j,k)*log(rho(i,j,k)/rho0(i,j,k))
            else if (iuserfunc==58) then
                cubmat(i,j,k)=rho(i,j,k)*(rholapl(i,j,k)/rho(i,j,k)-rholapl0(i,j,k)/rho0(i,j,k))
            else if (iuserfunc==59) then
                xtmp=derrho(1,i,j,k)/rho(i,j,k)-derrho0(1,i,j,k)/rho0(i,j,k)
                ytmp=derrho(2,i,j,k)/rho(i,j,k)-derrho0(2,i,j,k)/rho0(i,j,k)
                ztmp=derrho(3,i,j,k)/rho(i,j,k)-derrho0(3,i,j,k)/rho0(i,j,k)
                cubmat(i,j,k)=rho(i,j,k)*(xtmp**2+ytmp**2+ztmp**2)
            end if
        end do
    end do
end do
end subroutine





!------ Obtain information quantities for density difference between two wavefunctions. Adapted from intdiff
subroutine info_rhodiff
use defvar
use functions
use util
implicit real*8 (a-h,o-z)
real*8 intval,intvalold,intvalneg,intvaloldneg,funcval1(radpot*sphpot),funcval2(radpot*sphpot),beckeweigrid(radpot*sphpot)
type(content) gridatm(radpot*sphpot),gridatmorg(radpot*sphpot)
character c80tmp*80,filename2*200,sep
character(len=2) :: statname(-4:4)=(/ "-4","-3","-2","-1","_0","+1","+2","+3","+4" /)
real*8 atmraddens_Np1(ncenter,200),atmraddens_N(ncenter,200),atmraddens_Nn1(ncenter,200) !Radial density of each atom in N+1, N, N-1 states
integer atmradnpt_Np1(ncenter),atmradnpt_N(ncenter),atmradnpt_Nn1(ncenter) !Number of radial density points of each atom in N+1, N, N-1 states

write(*,"(a)") " The first wavefunction file is that you loaded after booting up Multiwfn, now input the path of the second wavefunction file, e.g. C:\yuri.wfn"
read(*,"(a)") filename2

!write(*,"(a)") " Select the type of calculation. The electron density in original formula is replaced with ""f=rho(wfn1)-rho(wfn2)"""
!write(*,*) "1 Fukui Shannon entropy"
!write(*,*) "2 Fukui relative Shannon entropy"
!write(*,*) "2 Fisher information density"
!write(*,*) "3 Second Fisher information density"
!read(*,*) itype
itype=1

!Prepare atomic .rad files for involved states
if (itype==2) then
	write(*,*) "How to calculate f0 in relative Shannon entropy?"
	write(*,*) "1 f0 is calculated as rho0(N+1) - rho0(N)"
	write(*,*) "2 f0 is calculated as rho0(N) - rho0(N-1)"
    read(*,*) if0
	if (if0==1) then
		ilow=-1
        ihigh=0
    else
		ilow=0
        ihigh=1
    end if
	sep='/'
	if (isys==1) sep='\'
	do iatm=1,ncenter
		iele=a_org(iatm)%index
		do istat=ilow,ihigh
			c80tmp="atmrad_spec"//sep//trim(a_org(iatm)%name)//statname(istat)
			inquire(file=trim(c80tmp)//".rad",exist=alive)
			if (alive) cycle
			inquire(file=trim(c80tmp)//".wfn",exist=alive)
			if (.not.alive) then
				write(*,"(' Error: ',a,' was not found!')") trim(c80tmp)//".wfn"
                write(*,*) "Press ENTER button to return"
                read(*,*)
				return
			end if
			write(*,"(' Converting ',a,' to ',a)") trim(c80tmp)//".wfn",trim(c80tmp)//".rad"
			call atmwfn2atmrad(trim(c80tmp)//".wfn",trim(c80tmp)//".rad")
		end do
	end do
    
    !Load radial density and points for each atom of N state
	do iatm=1,ncenter_org
		c80tmp="atmrad_spec"//sep//trim(a_org(iatm)%name)//"_0.rad"
		open(10,file=c80tmp,status="old")
		read(10,*) atmradnpt_N(iatm)
		do ipt=1,atmradnpt_N(iatm)
			read(10,*) rnouse,atmraddens_N(iatm,ipt)
		end do
		close(10)
	end do
    do iatm=1,ncenter_org
		if (if0==1) then !Load radial density and points for each atom of N+1 state
			c80tmp="atmrad_spec"//sep//trim(a_org(iatm)%name)//"-1.rad"
			open(10,file=c80tmp,status="old")
			read(10,*) atmradnpt_Np1(iatm)
			do ipt=1,atmradnpt_Np1(iatm)
				read(10,*) rnouse,atmraddens_Np1(iatm,ipt)
			end do
		else if (if0==2) then !Load radial density and points for each atom of N-1 state
			c80tmp="atmrad_spec"//sep//trim(a_org(iatm)%name)//"+1.rad"
			open(10,file=c80tmp,status="old")
			read(10,*) atmradnpt_Nn1(iatm)
			do ipt=1,atmradnpt_Nn1(iatm)
				read(10,*) rnouse,atmraddens_Nn1(iatm,ipt)
			end do
		end if
        close(10)
    end do
	call dealloall(0)
	call readinfile(firstfilename,1)
end if

!Start calculation
open(11,file="integrate.txt",status="replace")
write(*,"(/,' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)
call walltime(iwalltime1)
intval=0
intvalold=0
intvalneg=0
intvaloldneg=0
izero=0
ineg=0
do iatm=1,ncenter
	write(*,"(/,' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
	gridatm%x=gridatmorg%x+a(iatm)%x !Move quadrature point to actual position in molecule
	gridatm%y=gridatmorg%y+a(iatm)%y
	gridatm%z=gridatmorg%z+a(iatm)%z
	
	!Calculate data for wfn1
	!$OMP parallel do shared(funcval1) private(i) num_threads(nthreads) schedule(DYNAMIC)
	do i=1+iradcut*sphpot,radpot*sphpot
		if (itype==1) funcval1(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
	end do
	!$OMP end parallel do
	
	!Calculate data for wfn2
	call dealloall(0)
	call readinfile(filename2,1)
	!$OMP parallel do shared(funcval2) private(i) num_threads(nthreads) schedule(DYNAMIC)
	do i=1+iradcut*sphpot,radpot*sphpot
		if (itype==1) funcval2(i)=fdens(gridatm(i)%x,gridatm(i)%y,gridatm(i)%z)
	end do
	!$OMP end parallel do
    
	!Recover to wfn1
	call dealloall(0)
	call readinfile(firstfilename,1)
	
	call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
	do i=1+iradcut*sphpot,radpot*sphpot
		rhodiff=funcval1(i)-funcval2(i)
        accum=0
        if (rhodiff==0) then
			izero=izero+1
            inquire(file="zero_rhodiff.txt",number=ifilezero)
            if (ifilezero==-1) open(12,file="zero_rhodiff.txt",status="replace")
            write(12,"(i7,3f12.5,4(1PE16.8))") i,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,accum,rhodiff,gridatmorg(i)%value,beckeweigrid(i) 
        else if (rhodiff<0) then !Negative delta rho
			ineg=ineg+1
			if (itype==1) tmpval=-abs(rhodiff)*log(abs(rhodiff))
			accum=tmpval*gridatmorg(i)%value*beckeweigrid(i)
			intvalneg=intvalneg+accum
            inquire(file="neg_rhodiff.txt",number=ifileneg)
            if (ifileneg==-1) open(13,file="neg_rhodiff.txt",status="replace")
            write(13,"(i7,3f12.5,4(1PE16.8))") i,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,accum,rhodiff,gridatmorg(i)%value,beckeweigrid(i) 
        else !Positive delta rho
			if (itype==1) tmpval=-rhodiff*log(rhodiff)
			accum=tmpval*gridatmorg(i)%value*beckeweigrid(i)
			intval=intval+accum
			write(11,"(i7,3f12.5,4(1PE16.8))") i,gridatm(i)%x,gridatm(i)%y,gridatm(i)%z,accum,rhodiff,gridatmorg(i)%value,beckeweigrid(i)
		end if
	end do
	write(*,"(' Accumulated (+ delta_rho):',f20.10,'  This cen.:',f20.10)") intval,intval-intvalold
	write(*,"(' Accumulated (- delta_rho):',f20.10,'  This cen.:',f20.10)") intvalneg,intvalneg-intvaloldneg
	intvalold=intval
	intvaloldneg=intvalneg
end do

call walltime(iwalltime2)
write(*,"(/,' Calculation took up wall clock time',i10,'s',/)") iwalltime2-iwalltime1
write(*,"(' Final result (based on positive delta_rho):       ',f24.12)") intval
if (ineg>0) write(*,"(' Final result (based on abs of negative delta_rho):',f24.12)") intvalneg
write(*,*)
write(*,*) "integrate.txt has been exported to current folder"
write(*,*) "Column 1: Index of integration points"
write(*,*) "Columns 2~4: X, Y, Z of integration points in Bohr"
write(*,*) "Column 5: Contribution of this point to integral"
write(*,*) "Column 6: Density difference"
write(*,*) "Column 7: Lebedev integration weighting"
write(*,*) "Column 8: Atom weighting function"
write(*,*)
close(11)
if (izero>0) then
	close(12)
	write(*,"(a,i9,a)") " Note:",izero," points have zero density difference, these points were skipped &
	&during calculation. Please check zero_rhodiff.txt in current folder for information of these points"
end if
if (ineg>0) then
	close(13)
	write(*,"(a,i9,a)") " Note:",ineg," points have negative density difference, their absolute values were considered &
	&during integration. Please check neg_rhodiff.txt in current folder for information of these points"
end if
end subroutine



!------ Calculate Fukui Shannon density as |¦¤¦Ñ|*ln(|¦¤¦Ñ|) and store to cubmat. cubmat must contain ¦¤¦Ñ before entering this subroutine
subroutine info_rhodiff_grid
use defvar
implicit real*8 (a-h,o-z)
do k=1,nz
	do j=1,ny
		do i=1,nx
			absrho=abs(cubmat(i,j,k))
			cubmat(i,j,k)=-absrho*log(absrho)
		end do
	end do
end do
write(*,*) "Done!"
end subroutine





!!--------- The Energetic Information Project
!This routine use atomic grid to perform integration over whole space with Becke partition. Because Becke partition is sharp enough, molecular grid is not needed
subroutine energy_info_project
use defvar
use util
use functions
implicit real*8 (a-h,o-z)
character refsysname*200
real*8 intval,funcval(radpot*sphpot),beckeweigrid(radpot*sphpot)
type(content) gridatmorg(radpot*sphpot),gridatm(radpot*sphpot)
real*8 eval(radpot*sphpot),egrad(3,radpot*sphpot),elapl(radpot*sphpot)
real*8 e0val(radpot*sphpot),e0grad(3,radpot*sphpot)
real*8,allocatable :: evalcub(:,:,:),egradcub(:,:,:,:),elaplcub(:,:,:),e0valcub(:,:,:),e0gradcub(:,:,:,:),funcvalcub(:,:,:)
real*8 vectmp(3),hess(3,3)

do while(.true.) !Interface loop
write(*,*)
call menutitle("Energetic Information Project",10,1)
write(*,*) "0 Return"
write(*,*) "Select type of energetic information"
write(*,*) "1 Shannon entropy"
write(*,*) "2 Fisher information"
write(*,*) "3 Alternative Fisher information"
write(*,*) "4 Relative Shannon entropy"
write(*,*) "5 Relative Fisher information"
write(*,*) "6 Alternative relative Fisher information"
read(*,*) ieneinfo
if (ieneinfo==0) return

ider=0 !Evaluate value
if (ieneinfo==2.or.ieneinfo==5) ider=1 !Evaluate value and gradient
if (ieneinfo==3.or.ieneinfo==6) ider=2 !Evaluate value, gradient and Hessian

if (ieneinfo==4.or.ieneinfo==5.or.ieneinfo==6) then
	write(*,*) "Choose type of reference wavefunction"
	write(*,*) "1 Promolecule"
	write(*,*) "2 Another system"
	read(*,*) ireftype
	if (ireftype==1) then
		call generate_promolwfn(0) !Generate promolecular wavefunction and store to CO_pmol, MOocc_pmol, etc.
	else
		write(*,*) "Input path of wavefunction file of another system, e.g. /sob/test.wfn"
		do while(.true.)
			read(*,"(a)") refsysname
			inquire(file=refsysname,exist=alive)
			if (alive) exit
			write(*,*) "Cannot find the file, input again!"
		end do
	end if
end if

write(*,*)
write(*,*) "Select type of energy density"
write(*,"(a)") " 1 Total energy E (electronic energy density)"
write(*,"(a)") " 2 Total kinetic energy Ts (Hamiltonian kinetic energy density)"
write(*,"(a)") " -2 Total kinetic energy Ts (Lagrangian kinetic energy density)"
write(*,"(a)") " 3 Electrostatic energy Ee (negative ESP multiplied by electron density)"
write(*,"(a)") " 4 Exchange-correlation energy Exc (Integrand of Exc. The form is determined by ""iDFTxcsel"" in settings.ini)"
write(*,"(a)") " 5 Weizsacker kinetic energy Tw (closed-shell form)"
write(*,"(a)") " 6 Pauli energy Tp=Ts-Tw (the form of Ts is determined by ""iKEDsel"" in settings.ini)"
write(*,"(a)") " 7 Fermionic quantum energy Eq=Tp+Exc (evaluating Tp based on Hamiltonian kinetic energy density)"
write(*,"(a)") " -7 Fermionic quantum energy Eq=Tp+Exc (evaluating Tp based on Lagrangian kinetic energy density)"
read(*,*) ienedens

if (ienedens==2) then
	ifunc=6
else if (ienedens==-2) then
	ifunc=7
else
	ifunc=100
	if (ienedens==1) iuserfunc=11
	if (ienedens==3) iuserfunc=68
	if (ienedens==4) iuserfunc=1000
	if (ienedens==5) iuserfunc=5
	if (ienedens==6) iuserfunc=114
	if (ienedens==7) iuserfunc=69
	if (ienedens==-7) iuserfunc=-69
end if

write(*,*)
write(*,*) "Use origin form or scaled (normalized) form of the energy density?"
write(*,*) "1 Origin form"
write(*,*) "2 Scaled form"
read(*,*) iform

write(*,*)
write(*,*) "Choose task"
write(*,*) "1 Obtain integral over whole space"
write(*,*) "2 Calculate grid data and export to ITA.cub in current folder"
read(*,*) itask
if (itask==2) then
	call setgrid(0,igridsel)
    allocate(evalcub(nx,ny,nz),egradcub(3,nx,ny,nz),elaplcub(nx,ny,nz),e0valcub(nx,ny,nz),e0gradcub(3,nx,ny,nz),funcvalcub(nx,ny,nz))
end if

call walltime(iwalltime1)
write(*,"(' Radial points:',i5,'    Angular points:',i5,'   Total:',i10,' per center')") radpot,sphpot,radpot*sphpot
call gen1cintgrid(gridatmorg,iradcut)

!Calculate integral of energy density over whole space, so that we can use scaled (normalized) electron density later
if (iform==1) then
	eint=1
    e0int=1
else if (iform==2) then
	write(*,*) "Calculating integral of energy density over whole space..."
	eint=0 !Interal of present molecule to be evaluated
    e0int=0 !Integral of reference state to be evaluated
	do iatm=1,ncenter
		write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
		gridatm%x=gridatmorg%x+a(iatm)%x
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
    
		if (ienedens==3) call doinitlibreta(2)
    
		!Calculate energy density for present system
		!$OMP parallel do shared(eval) private(ipt) num_threads(nthreads)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			eval(ipt)=calcfuncall(ifunc,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
		end do
		!$OMP end parallel do
    
		!Calculate energy density for reference (promolecule)
		if (ieneinfo==4.or.ieneinfo==5.or.ieneinfo==6) then
			if (ireftype==1) then
				deallocate(MOocc,MOtype,MOene,CO)
				allocate(MOocc(nmo_pmol),MOene(nmo_pmol),MOtype(nmo_pmol),CO(nmo_pmol,nprims))
				nmo=nmo_pmol
				MOocc=MOocc_pmol
				MOene=MOene_pmol
				MOtype=MOtype_pmol
				CO=CO_pmol
			else if (ireftype==2) then
				call dealloall(0)
				call readinfile(refsysname,1)
			end if
			if (ienedens==3) call doinitlibreta(2)
			!$OMP parallel do shared(e0val) private(ipt) num_threads(nthreads)
			do ipt=1+iradcut*sphpot,radpot*sphpot
				e0val(ipt)=calcfuncall(ifunc,gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z)
			end do
			!$OMP end parallel do
			call dealloall(0)
			call readinfile(firstfilename,1) !Retrieve the first loaded file
		end if
		
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			eint=eint+eval(ipt)*gridatmorg(ipt)%value*beckeweigrid(ipt)
			e0int=e0int+e0val(ipt)*gridatmorg(ipt)%value*beckeweigrid(ipt)
		end do
	end do
    write(*,"(' Integral of present system:',1PE20.10)") eint
    if (ieneinfo==4.or.ieneinfo==5.or.ieneinfo==6) write(*,"(' Integral of reference system:',1PE20.10)") e0int
end if

!Start formal calculation
if (iform==2) then
	write(*,*)
	write(*,*) "Start formal calculation"
end if

if (itask==1) then !Obtain integral over whole space
	intval=0
	do iatm=1,ncenter
		write(*,"(' Processing center',i6,'(',a2,')   /',i6)") iatm,a(iatm)%name,ncenter
		gridatm%x=gridatmorg%x+a(iatm)%x
		gridatm%y=gridatmorg%y+a(iatm)%y
		gridatm%z=gridatmorg%z+a(iatm)%z
    
		if (ienedens==3) call doinitlibreta(2)
    
		!Calculate value (eval), gradient (egrad) and Laplacian (elapl) of energy density for present system
		!$OMP parallel do shared(eval,egrad,elapl) private(ipt,x,y,z,hess) num_threads(nthreads)
		do ipt=1+iradcut*sphpot,radpot*sphpot
			x=gridatm(ipt)%x
			y=gridatm(ipt)%y
			z=gridatm(ipt)%z
			if (ider==0) then !Only need value
				eval(ipt)=calcfuncall(ifunc,x,y,z)
			else !Need 1st or 1st+2nd derivative
				call gencalchessmat(ider,ifunc,x,y,z,eval(ipt),egrad(:,ipt),hess(:,:),1)
				elapl(ipt)=hess(1,1)+hess(2,2)+hess(3,3)
			end if
		end do
		!$OMP end parallel do
		if (iform==2) then
			eval=eval/eint
			egrad=egrad/eint
			elapl=elapl/eint
		end if
    
		!Calculate value, gradient and Laplacian of energy density for reference (promolecule), e0
		if (ieneinfo==4.or.ieneinfo==5.or.ieneinfo==6) then
			if (ireftype==1) then
				deallocate(MOocc,MOtype,MOene,CO)
				allocate(MOocc(nmo_pmol),MOene(nmo_pmol),MOtype(nmo_pmol),CO(nmo_pmol,nprims))
				nmo=nmo_pmol
				MOocc=MOocc_pmol
				MOene=MOene_pmol
				MOtype=MOtype_pmol
				CO=CO_pmol
			else if (ireftype==2) then
				call dealloall(0)
				call readinfile(refsysname,1)
			end if
			if (ienedens==3) call doinitlibreta(2)
			!$OMP parallel do shared(e0val,e0grad) private(ipt,x,y,z,hess) num_threads(nthreads)
			do ipt=1+iradcut*sphpot,radpot*sphpot
				x=gridatm(ipt)%x
				y=gridatm(ipt)%y
				z=gridatm(ipt)%z
				if (ider==0) then !Only need value
					e0val(ipt)=calcfuncall(ifunc,x,y,z)
				else !Need 1st derivative
					call gencalchessmat(1,ifunc,x,y,z,e0val(ipt),e0grad(:,ipt),hess(:,:),1)
				end if
			end do
			!$OMP end parallel do
			call dealloall(0)
			call readinfile(firstfilename,1) !Retrieve the first loaded file
		end if
		if (iform==2) then
			e0val=e0val/e0int
			e0grad=e0grad/e0int
		end if
    
		!Calculate function value at every point
		do ipt=1+iradcut*sphpot,radpot*sphpot
			if (ieneinfo==1) then
				funcval(ipt)=-eval(ipt)*log(eval(ipt))
			else if (ieneinfo==2) then
				funcval(ipt)=sum(egrad(:,ipt)**2)/eval(ipt)
			else if (ieneinfo==3) then
				funcval(ipt)=-elapl(ipt)*log(eval(ipt))
			else if (ieneinfo==4) then
				funcval(ipt)=eval(ipt)*log(eval(ipt)/e0val(ipt))
			else if (ieneinfo==5) then
				vectmp(:)=egrad(:,ipt)/eval(ipt)-e0grad(:,ipt)/e0val(ipt)
				funcval(ipt)=eval(ipt)*sum(vectmp(:)**2)
			else if (ieneinfo==6) then
				funcval(ipt)=elapl(ipt)*log(eval(ipt)/e0val(ipt))
			end if
			!write(15,*) ipt,funcval(ipt),eval(ipt),e0val(ipt),gridatm(ipt)%x,gridatm(ipt)%y,gridatm(ipt)%z
		end do
    
		call gen1cbeckewei(iatm,iradcut,gridatm,beckeweigrid,covr_tianlu,3)
		valthis=0
		do ipt=1+iradcut*sphpot,radpot*sphpot
			valthis=valthis+funcval(ipt)*gridatmorg(ipt)%value*beckeweigrid(ipt)
		end do
		intval=intval+valthis
		write(*,"(' Contribution of center',i6,':',1PE20.10)") iatm,valthis
	end do
	write(*,*) "Note: The center contributions correspond to Becke partition"

	write(*,"(/,' Result is',1PE20.10)") intval


else if (itask==2) then !Calculate grid data and export to ITA.cub in current folder
	if (ienedens==3) call doinitlibreta(2)
    write(*,"(a)") " Calculating grid data of energy density and its derivatives for actual system..."
	ifinish=0;ishowprog=1
	ntmp=floor(ny*nz/100D0)
	!$OMP PARALLEL DO SHARED(evalcub,egradcub,elaplcub,ifinish,ishowprog) PRIVATE(i,j,k,x,y,z,hess) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
	do k=1,nz
		do j=1,ny
			do i=1,nx
				call getgridxyz(i,j,k,x,y,z)
				if (ider==0) then !Only need value
					evalcub(i,j,k)=calcfuncall(ifunc,x,y,z)
				else !Need 1st or 1st+2nd derivative
					call gencalchessmat(ider,ifunc,x,y,z,evalcub(i,j,k),egradcub(:,i,j,k),hess(:,:),1)
					elaplcub(i,j,k)=hess(1,1)+hess(2,2)+hess(3,3)
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
    
	if (iform==2) then
		evalcub=evalcub/eint
		egradcub=egradcub/eint
		elaplcub=elaplcub/eint
	end if
    
    write(*,"(a)") " Calculating grid data of energy density and its derivatives for reference (promolecule) system..."
	if (ieneinfo==4.or.ieneinfo==5.or.ieneinfo==6) then
		if (ireftype==1) then
			deallocate(MOocc,MOtype,MOene,CO)
			allocate(MOocc(nmo_pmol),MOene(nmo_pmol),MOtype(nmo_pmol),CO(nmo_pmol,nprims))
			nmo=nmo_pmol
			MOocc=MOocc_pmol
			MOene=MOene_pmol
			MOtype=MOtype_pmol
			CO=CO_pmol
		else if (ireftype==2) then
			call dealloall(0)
			call readinfile(refsysname,1)
		end if
		if (ienedens==3) call doinitlibreta(2)
		ifinish=0;ishowprog=1
		!$OMP PARALLEL DO SHARED(e0valcub,e0gradcub,ifinish,ishowprog) PRIVATE(i,j,k,x,y,z,hess) schedule(dynamic) NUM_THREADS(nthreads) collapse(2)
		do k=1,nz
			do j=1,ny
				do i=1,nx
					call getgridxyz(i,j,k,x,y,z)
					if (ider==0) then !Only need value
						e0valcub(i,j,k)=calcfuncall(ifunc,x,y,z)
					else !Need 1st derivative
						call gencalchessmat(ider,ifunc,x,y,z,e0valcub(i,j,k),e0gradcub(:,i,j,k),hess(:,:),1)
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
		call dealloall(0)
		call readinfile(firstfilename,1) !Retrieve the first loaded file
	end if
	if (iform==2) then
		e0valcub=e0valcub/e0int
		e0gradcub=e0gradcub/e0int
	end if
    
    write(*,"(a)") " Calculating grid data of function value..."
	do k=1,nz
		do j=1,ny
			do i=1,nx
				if (ieneinfo==1) then
					funcvalcub(i,j,k)=-evalcub(i,j,k)*log(evalcub(i,j,k))
				else if (ieneinfo==2) then
					funcvalcub(i,j,k)=sum(egradcub(:,i,j,k)**2)/evalcub(i,j,k)
				else if (ieneinfo==3) then
					funcvalcub(i,j,k)=-elaplcub(i,j,k)*log(evalcub(i,j,k))
				else if (ieneinfo==4) then
					funcvalcub(i,j,k)=evalcub(i,j,k)*log(evalcub(i,j,k)/e0valcub(i,j,k))
				else if (ieneinfo==5) then
					vectmp(:)=egradcub(:,i,j,k)/evalcub(i,j,k)-e0gradcub(:,i,j,k)/e0valcub(i,j,k)
					funcvalcub(i,j,k)=evalcub(i,j,k)*sum(vectmp(:)**2)
				else if (ieneinfo==6) then
					funcvalcub(i,j,k)=elaplcub(i,j,k)*log(evalcub(i,j,k)/e0valcub(i,j,k))
				end if
			end do
		end do
	end do
    
    write(*,*) "Exporting ITA.cub in current folder"
	open(10,file="ITA.cub",status="replace")
	call outcube(funcvalcub,nx,ny,nz,orgx,orgy,orgz,gridv1,gridv2,gridv3,10)
    close(10)
    write(*,*) "Done!"
	if (itask==2) deallocate(evalcub,egradcub,elaplcub,e0valcub,e0gradcub,funcvalcub)
end if

call walltime(iwalltime2)
write(*,"(/,' Calculation took up wall clock time',i10,' s')") iwalltime2-iwalltime1


end do !End of interface loop

end subroutine





