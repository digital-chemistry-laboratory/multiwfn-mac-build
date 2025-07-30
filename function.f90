! ======================================
! =========== Module function ==========
! ======================================

module functions
use defvar
use libreta
implicit real*8 (a-h,o-z)

contains


!!-------- Calculate any supported real space function at a given point
real*8 function calcfuncall(ifunc,x,y,z)
integer ifunc
real*8 x,y,z
if (ifunc==1) then
	if (allocated(b)) then
		calcfuncall=fdens(x,y,z)
	else
		calcfuncall=calcprodens(x,y,z,0) !Use promolecular density
	end if
else if (ifunc==2) then
	calcfuncall=fgrad(x,y,z,'t')
else if (ifunc==3) then
	calcfuncall=flapl(x,y,z,'t')
else if (ifunc==4) then
	calcfuncall=fmo(x,y,z,iorbsel)
else if (ifunc==5) then
	calcfuncall=fspindens(x,y,z,'s')
else if (ifunc==6) then
	calcfuncall=Hamkin(x,y,z,0)
else if (ifunc==7) then
	calcfuncall=lagkin(x,y,z,0)
else if (ifunc==8) then
	calcfuncall=nucesp(x,y,z)
else if (ifunc==9) then
	calcfuncall=ELF_LOL(x,y,z,"ELF")
else if (ifunc==10) then
	calcfuncall=ELF_LOL(x,y,z,"LOL")
else if (ifunc==11) then
	calcfuncall=infoentro(1,x,y,z)
else if (ifunc==12) then
	calcfuncall=totesp(x,y,z)
else if (ifunc==13) then
	calcfuncall=fgrad(x,y,z,'r')
else if (ifunc==14) then
	calcfuncall=RDGprodens(x,y,z)
else if (ifunc==15) then
	calcfuncall=signlambda2rho(x,y,z)
else if (ifunc==16) then
	calcfuncall=signlambda2rho_prodens(x,y,z)
else if (ifunc==17) then
	calcfuncall=pairfunc(refx,refy,refz,x,y,z)
else if (ifunc==18) then
	calcfuncall=avglocion(x,y,z)
else if (ifunc==19) then
	calcfuncall=srcfunc(x,y,z,srcfuncmode)
else if (ifunc==20) then
	calcfuncall=edr(x,y,z) 
else if (ifunc==21) then
	calcfuncall=edrdmax(x,y,z)
else if (ifunc==22) then
	calcfuncall=delta_g_promol(x,y,z)
else if (ifunc==23) then
	calcfuncall=delta_g_Hirsh(x,y,z)
else if (ifunc==24) then
    calcfuncall=IRIfunc(x,y,z)
else if (ifunc==25) then
    calcfuncall=vdwpotfunc(x,y,z,1)
else if (ifunc==44) then
    calcfuncall=forbdens(x,y,z,iorbsel)
else if (ifunc==100) then
	calcfuncall=userfunc(x,y,z)
end if
end function




!!--------- User-defined function, the content is needed to be filled by users or selected by "iuserfunc" in settings.ini
real*8 function userfunc(x,y,z)
real*8 x,y,z,vec(3),mat(3,3)
userfunc=1D0 !Default value

select case(iuserfunc)
case (-10) !For test purpose
	userfunc=fmo(x,y,z,13)*fmo(x,y,z,20)
case (-3) !The function value evaluated by cubic spline interpolation from cubmat
    userfunc=splineintp3d(x,y,z,1)
case (-2) !Promolecular density
    userfunc=calcprodens(x,y,z,0)
case (-1) !The function value evaluated by trilinear interpolation from cubmat
    userfunc=linintp3d(x,y,z,1)
case (1) !Alpha density
    userfunc=fspindens(x,y,z,'a')
case (2) !Beta density
    userfunc=fspindens(x,y,z,'b')
case (3) !Integrand of electronic spatial extent <r^2>
    userfunc=(x*x+y*y+z*z)*fdens(x,y,z)
case (4) !Weizsacker potential
    userfunc=weizpot(x,y,z)
case (5) !Integrand of Weizsacker functional
    userfunc=KED(x,y,z,4) !Equivalent to weizsacker(x,y,z)
case (6) !Radial distribution function (assume that density is sphericalized)
    userfunc=4*pi*fdens(x,y,z)*(x*x+y*y+z*z)
case (7) !Local Temperature (Kelvin), PNAS,81,8028
	tmp=fdens(x,y,z)
	if (tmp>uservar) then
	    userfunc=2D0/3D0*lagkin(x,y,z,0)/tmp
    else
		userfunc=0
    end if
case (8) !Average local electrostatic potential, useful in exhibiting atomic shell structure
    userfunc=totesp(x,y,z)/fdens(x,y,z)
case (9) !Shape function
    userfunc=fdens(x,y,z)/nelec
case (10) !Potential energy density, also known as virial field
    userfunc=-Hamkin(x,y,z,0)-lagkin(x,y,z,0)
case (11) !Energy density
    userfunc=-Hamkin(x,y,z,0)
case (-11) !Scaled energy density
    userfunc=-Hamkin(x,y,z,0)*(virialratio-1)
case (12) !Local nuclear attraction potential energy
    userfunc=-nucesp(x,y,z)*fdens(x,y,z)
case (13) !This quantity at bond critical point is useful to discriminate covalent bonding and closed-shell interaction
    userfunc=lagkin(x,y,z,0)/fdens(x,y,z)
case (14) !Electrostatic potential from electrons
    userfunc=eleesp(x,y,z)
case (15) !Bond metallicity
    userfunc=fdens(x,y,z)/flapl(x,y,z,'t')
case (16) !Dimensionless bond metallicity
    userfunc=36*(3*pi*pi)**(2D0/3D0)/5D0*fdens(x,y,z)**(5D0/3D0)/flapl(x,y,z,'t')
case (17) !Energy density per electron
    userfunc=-Hamkin(x,y,z,0)/fdens(x,y,z)
case (18) !Region of Slow Electrons (RoSE), defined in Chem. Phys. Lett., 582, 144 (2013)
	rho=fdens(x,y,z)
	if (wfntype==0.or.wfntype==3) then !Closed-shell cases
		Dh=2.871234000D0*rho**(5D0/3D0)
	else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open shell cases
		rhospin=fspindens(x,y,z,'s') !rhospin=rhoa-rhob, rho=rhoa+rhob
		rhoa=(rhospin+rho)/2D0
		rhob=(rho-rhospin)/2D0
		Dh=4.557799872D0*(rhoa**(5D0/3D0)+rhob**(5D0/3D0)) !kinetic energy of HEG
	end if
	G=Lagkin(x,y,z,0)
	userfunc=(Dh-G)/(Dh+G)
case (19) !SEDD
    userfunc=SEDD(x,y,z)
case (20) !DORI
    userfunc=DORI(x,y,z)
case (21) !Integrand of X component of electric dipole moment
    userfunc=-x*fdens(x,y,z)
case (22) !Integrand of Y component of electric dipole moment
    userfunc=-y*fdens(x,y,z)
case (23) !Integrand of Z component of electric dipole moment
    userfunc=-z*fdens(x,y,z)
case (24) !Approximate form of DFT linear response kernel for closed-shell
    userfunc=linrespkernel(x,y,z)
case (25) !Magnitude of fluctuation of the electronic momentum
    userfunc=fgrad(x,y,z,'t')/fdens(x,y,z)/2D0
case (26) !Thomas-Fermi kinetic energy density
    userfunc=KED(x,y,z,3)
case (27) !Local electron affinity
    userfunc=loceleaff(x,y,z)
case (-27) !Local electron attachment energy
    userfunc=loceleatt(x,y,z)
case (28) !Local Mulliken electronegativity
    userfunc=(avglocion(x,y,z)+loceleaff(x,y,z))/2
case (29) !Local hardness
    userfunc=(avglocion(x,y,z)-loceleaff(x,y,z))/2
case (30) !Ellipticity of electron density
    userfunc=densellip(x,y,z,1)
case (31) !eta index, Angew. Chem. Int. Ed., 53, 2766-2770 (2014)
    userfunc=densellip(x,y,z,2)
case (32) !Modified eta index
    userfunc=densellip(x,y,z,2)-1
case (33) !PAEM, potential acting on one electron in a molecule, defined by Zhongzhi Yang
    userfunc=PAEM(x,y,z,1)
case (34) !The same as 33, but using DFT XC potential directly rather than evaluating the XC potential based on pair density
    userfunc=PAEM(x,y,z,2)
case (35) !|V(r)|/G(r)
	tmpval=lagkin(x,y,z,0)
	userfunc=abs(-Hamkin(x,y,z,0)-tmpval)/tmpval
case (36) !On-top pair density, i.e. r1=r2 case of pair density. paircorrtype affects result
	pairfunctypeold=pairfunctype
	pairfunctype=12
	userfunc=pairfunc(x,y,z,x,y,z)
	pairfunctype=pairfunctypeold
case (37) !SCI
    userfunc=ELF_LOL(x,y,z,"SCI")
case (38) !The angle between the second eigenvector of rho and the plane defined by option 4 of main function 1000
    userfunc=Ang_rhoeigvec_ple(x,y,z,2)
case (39) !ESP without contribution of nuclues "iskipnuc"
    userfunc=totespskip(x,y,z,iskipnuc)
case (40) !Steric energy density
    userfunc=weizsacker(x,y,z)
case (41) !Steric potential
    userfunc=stericpot(x,y,z)
case (42) !Steric charge
    userfunc=stericcharge(x,y,z)
case (43) !The magnitude of steric force
    userfunc=stericforce(x,y,z)
case (44) !Steric potential with damping function to a given constant value
    userfunc=stericpot_damp(x,y,z)
case (45) !Steric force based on damped potential
    userfunc=stericforce_damp(x,y,z)
case (46) !Steric force directly damped to zero
    userfunc=stericforce_directdamp(x,y,z)
case (47) !Steric charge directly damped to zero
    userfunc=stericcharge_directdamp(x,y,z)
case (49) !Relative Shannon entropy, also called information gain
    userfunc=relShannon(x,y,z)
case (50) !Shannon entropy density, see JCP,126,191107 for example
    userfunc=infoentro(2,x,y,z)
case (51) !Fisher information density, see JCP,126,191107 for example
    userfunc=Fisherinfo(1,x,y,z)
case (52) !Second Fisher information density, see JCP,126,191107 for derivation
    userfunc=Fisherinfo(2,x,y,z)
case (53) !Ghosh entropy density with G(r) as kinetic energy density, PNAS,81,8028
    userfunc=Ghoshentro(x,y,z,1)
case (54) !Ghosh entropy density with G(r)-der2rho/8 as kinetic energy density, exactly corresponds to Eq.22 in PNAS,81,8028
    userfunc=Ghoshentro(x,y,z,2)
case (55) !Integrand of quadratic form of Renyi entropy
    userfunc=fdens(x,y,z)**2
case (56) !Integrand of cubic form of Renyi entropy
    userfunc=fdens(x,y,z)**3
case (60) !Pauli potential, Comp. Theor. Chem., 1006, 92-99
    userfunc=paulipot(x,y,z)
case (61) !Magnitude of Pauli force
    userfunc=pauliforce(x,y,z)
case (62) !Pauli charge
    userfunc=paulicharge(x,y,z)
case (63) !Quantum potential
    userfunc=quantumpot(x,y,z)
case (64) !The magnitude of quantum force
    userfunc=quantumforce(x,y,z)
case (65) !Quantum charge
    userfunc=quantumcharge(x,y,z)
case (66) !The magnitude of electrostatic force
    userfunc=elestatforce(x,y,z)
case (67) !Electrostatic charge
    userfunc=elestatcharge(x,y,z)
case (68) !Energy density of electronic part of electrostatic term (Ee) of SBL's energy decomposition
    userfunc=-fdens(x,y,z)*totesp(x,y,z)
case (69) !Energy density of quantum term (Eq) of SBL's energy decomposition using Hamiltonian kinetic energy
    userfunc=Hamkin(x,y,z,0)-weizsacker(x,y,z)+DFTxcfunc(x,y,z)
case (-69) !Energy density of quantum term (Eq) of SBL's energy decomposition using Lagrangian kinetic energy
    userfunc=Lagkin(x,y,z,0)-weizsacker(x,y,z)+DFTxcfunc(x,y,z)
case (70) !Phase-space-defined Fisher information density
    userfunc=4.5D0*fdens(x,y,z)**2/lagkin(x,y,z,0)
case (71) !X component of electron linear momentum density in 3D representation
    userfunc=elemomdens(x,y,z,1)
case (72) !Y component of electron linear momentum density in 3D representation
    userfunc=elemomdens(x,y,z,2)
case (73) !Z component of electron linear momentum density in 3D representation
    userfunc=elemomdens(x,y,z,3)
case (74) !Magnitude of electron linear momentum density in 3D representation
    userfunc=elemomdens(x,y,z,0)
case (75) !X component of magnetic dipole moment density
    userfunc=magmomdens(x,y,z,1)
case (76) !Y component of magnetic dipole moment density
    userfunc=magmomdens(x,y,z,2)
case (77) !Z component of magnetic dipole moment density
    userfunc=magmomdens(x,y,z,3)
case (78) !Magnitude of magnetic dipole moment density
    userfunc=magmomdens(x,y,z,0)
case (79) !Gradient norm of energy density
    userfunc=energydens_grdn(x,y,z)
case (80) !Laplacian of energy density
    userfunc=energydens_lapl(x,y,z)
case (81) !X component of Hamiltonian kinetic energy density
    userfunc=hamkin(x,y,z,1)
case (82) !Y component of Hamiltonian kinetic energy density
    userfunc=hamkin(x,y,z,2)
case (83) !Z component of Hamiltonian kinetic energy density
    userfunc=hamkin(x,y,z,3)
case (84) !X component of Lagrangian kinetic energy density
    userfunc=Lagkin(x,y,z,1)
case (85) !Y component of Lagrangian kinetic energy density
    userfunc=Lagkin(x,y,z,2)
case (86) !Z component of Lagrangian kinetic energy density
    userfunc=Lagkin(x,y,z,3)
case (87) !Local total electron correlation function
    userfunc=localcorr(x,y,z,1)
case (88) !Local dynamic electron correlation function
    userfunc=localcorr(x,y,z,2)
case (89) !Local nondynamic electron correlation function
    userfunc=localcorr(x,y,z,3)
case (90) !Fractional occupation number weighted electron density (FOD)
    userfunc=FODfunc(x,y,z)
case (91) !Delta-g_inter(Hirsh) defined in IGMH model
    userfunc=delta_g_inter_Hirsh(x,y,z)
case (92) !vdW potential
    userfunc=vdwpotfunc(x,y,z,1)
case (93) !Repulsion potential
    userfunc=vdwpotfunc(x,y,z,2)
case (94) !Disperison potential
    userfunc=vdwpotfunc(x,y,z,3)
case (95) !Orbital-weighted f+ Fukui function
    userfunc=orbwei_Fukui(1,x,y,z)
case (96) !Orbital-weighted f- Fukui function
    userfunc=orbwei_Fukui(2,x,y,z)
case (97) !Orbital-weighted f0 Fukui function
    userfunc=orbwei_Fukui(3,x,y,z)
case (98) !Orbital-weighted dual descriptor
    userfunc=orbwei_Fukui(4,x,y,z)
case (99) !Interaction region indicator (IRI)
    userfunc=IRIfunc(x,y,z)
case (100) !Disequilibrium (also known as semi-similarity)
    userfunc=fdens(x,y,z)**2
case (101) !Positive part of ESP
	userfunc=totesp(x,y,z)
	if (userfunc<0D0) userfunc=0D0
case (102) !Negative part of ESP
	userfunc=totesp(x,y,z)
	if (userfunc>0D0) userfunc=0D0
case (103) !Magnitude of electric field
	call gencalchessmat(1,12,x,y,z,value,vec,mat) !Get gradient of ESP
	userfunc=dsqrt(sum(vec**2))
case (110) !Total energy density of the energy components defined by SBL (steric + electrostatic + quantum)
    !That is, userfunc(40) + userfunc(68) + userfunc(69)
    !userfunc = weizsacker(x,y,z) - fdens(x,y,z)*totesp(x,y,z) + (Hamkin(x,y,z,0)-weizsacker(x,y,z)+DFTxcfunc(x,y,z)) !Original expression
    userfunc = -fdens(x,y,z)*totesp(x,y,z) + Hamkin(x,y,z,0) + DFTxcfunc(x,y,z)
case (111) !Total potential of the energy components defined by SBL (steric + electrostatic + quantum)
    userfunc = stericpot(x,y,z) - totesp(x,y,z) + quantumpot(x,y,z)
case (112) !Magnitude of total force of the energy components defined by SBL (steric + electrostatic + quantum)
    !That is, vector sum of userfunc(43), userfunc(66), userfunc(64)
    !userfunc= stericforce(x,y,z) + elestatforce(x,y,z) + quantumforce(x,y,z) !This is sum of magnitude, not what we want
    userfunc = SBLallforce(x,y,z)
case (113) !Total charge of the energy components defined by SBL (steric + electrostatic + quantum)
    !That is, userfunc(42) + userfunc(67) + userfunc(65)
    userfunc = stericcharge(x,y,z) + elestatcharge(x,y,z) + quantumcharge(x,y,z)
case (114) !Pauli kinetic energy density
    userfunc = KED(x,y,z,iKEDsel) - weizsacker(x,y,z)
case (115) !Stiffness
    userfunc = densellip(x,y,z,3)
case (116) !Stress tensor stiffness
    userfunc = stress_stiffness(x,y,z)
case (117) !Stress tensor polarizability
    userfunc = 1D0/stress_stiffness(x,y,z)
case (118) !Stress tensor ellipticity
    userfunc = stress_ellipticity(x,y,z)
case (200) !Random number of [0,1£©
	call RANDOM_NUMBER(userfunc)
case (802:807)
    userfunc=funcvalLSB(iuserfunc-800,x,y,z)
case (812:817)
    userfunc=1/funcvalLSB(iuserfunc-810,x,y,z)
case (819) !Ultrastrong interaction (USI)
    userfunc=flapl(x,y,z,'t')/fdens(x,y,z)**(5D0/3D0)
case (820) !Bonding and noncovalent interaction (BNI): (tau-t_w)/t_w
    tmp=weizsacker(x,y,z)
    userfunc=(lagkin(x,y,z,0)-tmp)/tmp
case (821) !rho^(4/3)/|der_rho|, for shubin
    userfunc=fdens(x,y,z)**(4D0/3D0) / fgrad(x,y,z,'t')
case (900) !X coordinate
    userfunc=x
case (901) !Y coordinate
    userfunc=y
case (902) !Z coordinate
    userfunc=z
case (910) !Hirshfeld atomic weighting function
	call Hirshatmwei(nint(uservar),x,y,z,userfunc)
case (911) !Becke atomic weighting function using covr_tianlu
	call Beckeatmwei(nint(uservar),x,y,z,userfunc,covr_tianlu,3)
case (912) !Becke atomic weighting function using CSD covalent radii
	call Beckeatmwei(nint(uservar),x,y,z,userfunc,covr,3)
case (913) !LT1 atomic weighting function
	call TLatmwei(nint(uservar),x,y,z,userfunc,1)
case (914) !LT1 atomic weighting function
	call TLatmwei(nint(uservar),x,y,z,userfunc,2)
case (999) !Local Hartree-Fock exchange energy
	userfunc=locHFexc(x,y,z)
case (1000) !Various kinds of DFT exchange-correlation functions
    userfunc=DFTxcfunc(x,y,z)
case (1100) !Various kinds of DFT exchange-correlation potentials for closed-shell
    userfunc=DFTxcpot(x,y,z,0)
case (1101) !Various kinds of DFT exchange-correlation potentials of alpha electrons for open-shell
    userfunc=DFTxcpot(x,y,z,1)
case (1102) !Various kinds of DFT exchange-correlation potentials of beta electrons for open-shell
    userfunc=DFTxcpot(x,y,z,2)
case (1200) !Various kinds of electronic kinetic energy density (KED)
    userfunc=KED(x,y,z,iKEDsel)
case (1201) !Get difference between KED selected by iKEDsel and Weizsacker KED
    userfunc=KEDdiff(x,y,z,iKEDsel,1)
case (1202) !Get difference between KED selected by iKEDsel and Lagrangian KED
    userfunc=KEDdiff(x,y,z,iKEDsel,2)
case (1203) !Get absolute difference between KED selected by iKEDsel and Lagrangian KED
    userfunc=KEDdiff(x,y,z,iKEDsel,3)
case (1204) !Local Temperature(Kelvin), PNAS,81,8028, but using user-defined KED by "iKEDsel"
	tmp=fdens(x,y,z)
	if (tmp>uservar) then
	    userfunc=2D0/3D0*KED(x,y,z,iKEDsel)/tmp
    else
		userfunc=0
    end if
case (1210) !Potential of KED
    userfunc=KEDpot(x,y,z)
case (1303) !\xi^\alpha Fractional integrals/derivatives, close to Lagkin, fdens routines
    userfunc=fracderiv(x,y,z)
end select

if (iuserfunc>10000) then
	userfunc=calcfuncall(iuserfunc-10000,x,y,z)
end if

!Below are other examples
! userfunc=hamkin(x,y,z,3)-0.5D0*(hamkin(x,y,z,1)+hamkin(x,y,z,2)) !Anisotropy of Hamiltonian kinetic energy in Z, namely K_Z-0.5*(K_X+K_Y)
! userfunc=-x*y*fdens(x,y,z) !Integrand of XY component of electric quadrupole moment
! userfunc=-x*y*z*fdens(x,y,z)*au2debye*b2a**2 !Integrand of XYZ component of electric octapole moment in Debye-Ang**2
end function




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate wavefunction value of a range of orbitals and their derivatives at a given point, up to third-order
!! istart and iend is the range of the orbitals will be calculated, to calculate all orbitals, use 1,nmo
!! runtype=1: value  =2: value+dx/y/z  =3: value+dxx/yy/zz (diagonal of hess)  =4: value+dx/y/z+Hessian  
!!        =5: value+dx/y/z+hess+3-order derivative tensor 
subroutine orbderv(runtype,istart,iend,x,y,z,wfnval,grad,hess,tens3)
real*8 x,y,z,wfnval(nmo)
real*8,optional :: grad(3,nmo),hess(3,3,nmo),tens3(3,3,3,nmo)
integer runtype,istart,iend

if (ifPBC>0) then !Consider PBC
    call orbderv_PBC(runtype,istart,iend,x,y,z,wfnval,grad,hess,tens3)
    return
end if

wfnval=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0
lastcen=-1 !Arbitrary value

!If the center/exp of current GTF is the same as previous, then we do not need to recalculate them
if (nprims_uniq==0) then !Unique GTF information is not available
	do j=1,nprims
		jcen=b(j)%center
		if (jcen/=lastcen) then
			sftx=x-a(jcen)%x
			sfty=y-a(jcen)%y
			sftz=z-a(jcen)%z
			sftx2=sftx*sftx
			sfty2=sfty*sfty
			sftz2=sftz*sftz
			rr=sftx2+sfty2+sftz2
		end if
        
		ep=b(j)%exp
		tmpval=-ep*rr
		lastcen=jcen
		if (tmpval>expcutoff.or.expcutoff>0) then
			expterm=exp(tmpval)
		else
			cycle
		end if
	
		!Calculate value for current GTF
		jtype=b(j)%type
		ix=type2ix(jtype)
		iy=type2iy(jtype)
		iz=type2iz(jtype)
		if (jtype==1) then !Some functype use manually optimized formula for cutting down computational time
		GTFval=expterm
		else if (jtype==2) then
		GTFval=sftx*expterm
		else if (jtype==3) then
		GTFval=sfty*expterm
		else if (jtype==4) then
		GTFval=sftz*expterm
		else if (jtype==5) then
		GTFval=sftx2*expterm
		else if (jtype==6) then
		GTFval=sfty2*expterm
		else if (jtype==7) then
		GTFval=sftz2*expterm
		else if (jtype==8) then
		GTFval=sftx*sfty*expterm
		else if (jtype==9) then
		GTFval=sftx*sftz*expterm
		else if (jtype==10) then
		GTFval=sfty*sftz*expterm
		else !If above conditions are not satisfied (angular moment higher than f), the function will be calculated explicitly
		GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
		end if
		!Calculate orbital wavefunction value. This is bottle neck of cost of this function
		do imo=istart,iend
			wfnval(imo)=wfnval(imo)+CO(imo,j)*GTFval
		end do
	
		if (runtype>=2) then
			!Calculate 1-order derivative for current GTF
			tx=0D0
			ty=0D0
			tz=0D0
			if (ix/=0) tx=ix*sftx**(ix-1)
			if (iy/=0) ty=iy*sfty**(iy-1)
			if (iz/=0) tz=iz*sftz**(iz-1)
			GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
			GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
			GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
			!Calculate 1-order derivative for orbitals
			do imo=istart,iend
				grad(1,imo)=grad(1,imo)+CO(imo,j)*GTFdx
				grad(2,imo)=grad(2,imo)+CO(imo,j)*GTFdy
				grad(3,imo)=grad(3,imo)+CO(imo,j)*GTFdz
			end do

			if (runtype>=3) then
				!Calculate 2-order derivative for current GTF
				txx=0D0
				tyy=0D0
				tzz=0D0
				if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
				if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
				if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
				GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
				GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
				GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
				ttx=tx-2*ep*sftx**(ix+1)
				tty=ty-2*ep*sfty**(iy+1)
				ttz=tz-2*ep*sftz**(iz+1)
				GTFdxy=sftz**iz *expterm*ttx*tty
				GTFdyz=sftx**ix *expterm*tty*ttz
				GTFdxz=sfty**iy *expterm*ttx*ttz
				!Calculate diagonal Hessian elements for orbitals
				do imo=istart,iend
					hess(1,1,imo)=hess(1,1,imo)+CO(imo,j)*GTFdxx !dxx
					hess(2,2,imo)=hess(2,2,imo)+CO(imo,j)*GTFdyy !dyy
					hess(3,3,imo)=hess(3,3,imo)+CO(imo,j)*GTFdzz !dzz
				end do
				if (runtype>=4) then !Also process nondiagonal elements
					do imo=istart,iend
						hess(1,2,imo)=hess(1,2,imo)+CO(imo,j)*GTFdxy !dxy
						hess(2,3,imo)=hess(2,3,imo)+CO(imo,j)*GTFdyz !dyz
						hess(1,3,imo)=hess(1,3,imo)+CO(imo,j)*GTFdxz !dxz
					end do
					hess(2,1,:)=hess(1,2,:)
					hess(3,2,:)=hess(2,3,:)
					hess(3,1,:)=hess(1,3,:)
				end if
			
				if (runtype>=5) then
					!Calculate 3-order derivative for current GTF
					ep2=ep*2D0
					ep4=ep*4D0
					epep4=ep2*ep2
					epep8=epep4*2D0
					!dxyz
					a1=0D0
					b1=0D0
					c1=0D0
					if (ix>=1) a1=ix*sftx**(ix-1)
					if (iy>=1) b1=iy*sfty**(iy-1)
					if (iz>=1) c1=iz*sftz**(iz-1)
					a2=-ep2*sftx**(ix+1)
					b2=-ep2*sfty**(iy+1)
					c2=-ep2*sftz**(iz+1)
					GTFdxyz=(a1+a2)*(b1+b2)*(c1+c2)*expterm
					!dxyy,dxxy,dxxz,dxzz,dyzz,dyyz
					atmp=0D0
					btmp=0D0
					ctmp=0D0
					if (ix>=2) atmp=ix*(ix-1)*sftx**(ix-2)
					if (iy>=2) btmp=iy*(iy-1)*sfty**(iy-2)
					if (iz>=2) ctmp=iz*(iz-1)*sftz**(iz-2)
					GTFdxyy=(a1+a2)*sftz**iz *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy)
					GTFdxxy=(b1+b2)*sftz**iz *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dyxx
					GTFdxxz=(c1+c2)*sfty**iy *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dzxx
					GTFdxzz=(a1+a2)*sfty**iy *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
					GTFdyzz=(b1+b2)*sftx**ix *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
					GTFdyyz=(c1+c2)*sftx**ix *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !=dzyy
					!dxxx,dyyy,dzzz
					aatmp1=0D0
					bbtmp1=0D0
					cctmp1=0D0
					if (ix>=1) aatmp1=ep2*ix*sftx**(ix-1)
					if (iy>=1) bbtmp1=ep2*iy*sfty**(iy-1)
					if (iz>=1) cctmp1=ep2*iz*sftz**(iz-1)
					aatmp2=0D0
					bbtmp2=0D0
					cctmp2=0D0
					if (ix>=2) aatmp2=ep2*ix*(ix-1)*sftx**(ix-1)
					if (iy>=2) bbtmp2=ep2*iy*(iy-1)*sfty**(iy-1)
					if (iz>=2) cctmp2=ep2*iz*(iz-1)*sftz**(iz-1)
					aatmp3=0D0
					bbtmp3=0D0
					cctmp3=0D0
					if (ix>=3) aatmp3=ix*(ix-1)*(ix-2)*sftx**(ix-3)
					if (iy>=3) bbtmp3=iy*(iy-1)*(iy-2)*sfty**(iy-3)
					if (iz>=3) cctmp3=iz*(iz-1)*(iz-2)*sftz**(iz-3)
					GTFdxxx=sfty**iy*sftz**iz*expterm*( (-2*ix+ep2*sftx2-1)*(-epep4*sftx**(ix+1) + aatmp1) - aatmp2 + epep8*sftx**(ix+1) + aatmp3 )
					GTFdyyy=sftx**ix*sftz**iz*expterm*( (-2*iy+ep2*sfty2-1)*(-epep4*sfty**(iy+1) + bbtmp1) - bbtmp2 + epep8*sfty**(iy+1) + bbtmp3 )
					GTFdzzz=sfty**iy*sftx**ix*expterm*( (-2*iz+ep2*sftz2-1)*(-epep4*sftz**(iz+1) + cctmp1) - cctmp2 + epep8*sftz**(iz+1) + cctmp3 )
				
					!Calculate 3-order derivative tensor for orbital wavefunction
					do imo=istart,iend
						tens3(1,1,1,imo)=tens3(1,1,1,imo)+CO(imo,j)*GTFdxxx !dxxx
						tens3(2,2,2,imo)=tens3(2,2,2,imo)+CO(imo,j)*GTFdyyy !dyyy
						tens3(3,3,3,imo)=tens3(3,3,3,imo)+CO(imo,j)*GTFdzzz !dzzz
						tens3(1,2,2,imo)=tens3(1,2,2,imo)+CO(imo,j)*GTFdxyy !dxyy
						tens3(1,1,2,imo)=tens3(1,1,2,imo)+CO(imo,j)*GTFdxxy !dxxy
						tens3(1,1,3,imo)=tens3(1,1,3,imo)+CO(imo,j)*GTFdxxz !dxxz
						tens3(1,3,3,imo)=tens3(1,3,3,imo)+CO(imo,j)*GTFdxzz !dxzz
						tens3(2,3,3,imo)=tens3(2,3,3,imo)+CO(imo,j)*GTFdyzz !dyzz
						tens3(2,2,3,imo)=tens3(2,2,3,imo)+CO(imo,j)*GTFdyyz !dyyz
						tens3(1,2,3,imo)=tens3(1,2,3,imo)+CO(imo,j)*GTFdxyz !dxyz
					end do
					tens3(1,2,1,:)=tens3(1,1,2,:) !dxyx=dxxy
					tens3(1,3,1,:)=tens3(1,1,3,:) !dxzx=dxxz
					tens3(1,3,2,:)=tens3(1,2,3,:) !dxzy=dxyz
					tens3(2,1,1,:)=tens3(1,1,2,:) !dyxx=dxxy
					tens3(2,1,2,:)=tens3(1,2,2,:) !dyxy=dxyy
					tens3(2,1,3,:)=tens3(1,2,3,:) !dyxz=dxyz
					tens3(2,2,1,:)=tens3(1,2,2,:) !dyyx=dxyy
					tens3(2,3,1,:)=tens3(1,2,3,:) !dyzx=dxyz
					tens3(2,3,2,:)=tens3(2,2,3,:) !dyzy=dyyz
					tens3(3,1,1,:)=tens3(1,1,3,:) !dzxx=dxxz
					tens3(3,1,2,:)=tens3(1,2,3,:) !dzxy=dxyz
					tens3(3,1,3,:)=tens3(1,3,3,:) !dzxz=dxzz
					tens3(3,2,1,:)=tens3(1,2,3,:) !dzyx=dxyz
					tens3(3,2,2,:)=tens3(2,2,3,:) !dzyy=dyyz
					tens3(3,2,3,:)=tens3(2,3,3,:) !dzyz=dyzz
					tens3(3,3,1,:)=tens3(1,3,3,:) !dzzx=dxzz
					tens3(3,3,2,:)=tens3(2,3,3,:) !dzzy=dyzz
				end if !end runtype>=5
			
			end if !end runtype>=3
		end if !end runtype>=2
	end do

else !Unique GTF information has been constructed by gen_GTFuniq
	do j=1,nprims_uniq
		jcen=b_uniq(j)%center
		if (jcen/=lastcen) then
			sftx=x-a(jcen)%x
			sfty=y-a(jcen)%y
			sftz=z-a(jcen)%z
			sftx2=sftx*sftx
			sfty2=sfty*sfty
			sftz2=sftz*sftz
			rr=sftx2+sfty2+sftz2
		end if
        
		ep=b_uniq(j)%exp
		tmpval=-ep*rr
		lastcen=jcen
		if (tmpval>expcutoff.or.expcutoff>0) then
			expterm=exp(tmpval)
		else
			cycle
		end if
	
		!Calculate value for current GTF
		jtype=b_uniq(j)%type
		ix=type2ix(jtype)
		iy=type2iy(jtype)
		iz=type2iz(jtype)
		if (jtype==1) then !Some functype use manually optimized formula for cutting down computational time
		GTFval=expterm
		else if (jtype==2) then
		GTFval=sftx*expterm
		else if (jtype==3) then
		GTFval=sfty*expterm
		else if (jtype==4) then
		GTFval=sftz*expterm
		else if (jtype==5) then
		GTFval=sftx2*expterm
		else if (jtype==6) then
		GTFval=sfty2*expterm
		else if (jtype==7) then
		GTFval=sftz2*expterm
		else if (jtype==8) then
		GTFval=sftx*sfty*expterm
		else if (jtype==9) then
		GTFval=sftx*sftz*expterm
		else if (jtype==10) then
		GTFval=sfty*sftz*expterm
		else !If above conditions are not satisfied (angular moment higher than f), the function will be calculated explicitly
		GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
		end if
		!Calculate orbital wavefunction value. This is bottle neck of cost of this function
		do imo=istart,iend
			wfnval(imo)=wfnval(imo)+CO_uniq(imo,j)*GTFval
		end do
	
		if (runtype>=2) then
			!Calculate 1-order derivative for current GTF
			tx=0D0
			ty=0D0
			tz=0D0
			if (ix/=0) tx=ix*sftx**(ix-1)
			if (iy/=0) ty=iy*sfty**(iy-1)
			if (iz/=0) tz=iz*sftz**(iz-1)
			GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
			GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
			GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
			!Calculate 1-order derivative for orbitals
			do imo=istart,iend
				grad(1,imo)=grad(1,imo)+CO_uniq(imo,j)*GTFdx
				grad(2,imo)=grad(2,imo)+CO_uniq(imo,j)*GTFdy
				grad(3,imo)=grad(3,imo)+CO_uniq(imo,j)*GTFdz
			end do

			if (runtype>=3) then
				!Calculate 2-order derivative for current GTF
				txx=0D0
				tyy=0D0
				tzz=0D0
				if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
				if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
				if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
				GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
				GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
				GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
				ttx=tx-2*ep*sftx**(ix+1)
				tty=ty-2*ep*sfty**(iy+1)
				ttz=tz-2*ep*sftz**(iz+1)
				GTFdxy=sftz**iz *expterm*ttx*tty
				GTFdyz=sftx**ix *expterm*tty*ttz
				GTFdxz=sfty**iy *expterm*ttx*ttz
				!Calculate diagonal Hessian elements for orbitals
				do imo=istart,iend
					hess(1,1,imo)=hess(1,1,imo)+CO_uniq(imo,j)*GTFdxx !dxx
					hess(2,2,imo)=hess(2,2,imo)+CO_uniq(imo,j)*GTFdyy !dyy
					hess(3,3,imo)=hess(3,3,imo)+CO_uniq(imo,j)*GTFdzz !dzz
				end do
				if (runtype>=4) then !Also process nondiagonal elements
					do imo=istart,iend
						hess(1,2,imo)=hess(1,2,imo)+CO_uniq(imo,j)*GTFdxy !dxy
						hess(2,3,imo)=hess(2,3,imo)+CO_uniq(imo,j)*GTFdyz !dyz
						hess(1,3,imo)=hess(1,3,imo)+CO_uniq(imo,j)*GTFdxz !dxz
					end do
					hess(2,1,:)=hess(1,2,:)
					hess(3,2,:)=hess(2,3,:)
					hess(3,1,:)=hess(1,3,:)
				end if
			
				if (runtype>=5) then
					!Calculate 3-order derivative for current GTF
					ep2=ep*2D0
					ep4=ep*4D0
					epep4=ep2*ep2
					epep8=epep4*2D0
					!dxyz
					a1=0D0
					b1=0D0
					c1=0D0
					if (ix>=1) a1=ix*sftx**(ix-1)
					if (iy>=1) b1=iy*sfty**(iy-1)
					if (iz>=1) c1=iz*sftz**(iz-1)
					a2=-ep2*sftx**(ix+1)
					b2=-ep2*sfty**(iy+1)
					c2=-ep2*sftz**(iz+1)
					GTFdxyz=(a1+a2)*(b1+b2)*(c1+c2)*expterm
					!dxyy,dxxy,dxxz,dxzz,dyzz,dyyz
					atmp=0D0
					btmp=0D0
					ctmp=0D0
					if (ix>=2) atmp=ix*(ix-1)*sftx**(ix-2)
					if (iy>=2) btmp=iy*(iy-1)*sfty**(iy-2)
					if (iz>=2) ctmp=iz*(iz-1)*sftz**(iz-2)
					GTFdxyy=(a1+a2)*sftz**iz *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy)
					GTFdxxy=(b1+b2)*sftz**iz *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dyxx
					GTFdxxz=(c1+c2)*sfty**iy *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dzxx
					GTFdxzz=(a1+a2)*sfty**iy *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
					GTFdyzz=(b1+b2)*sftx**ix *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
					GTFdyyz=(c1+c2)*sftx**ix *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !=dzyy
					!dxxx,dyyy,dzzz
					aatmp1=0D0
					bbtmp1=0D0
					cctmp1=0D0
					if (ix>=1) aatmp1=ep2*ix*sftx**(ix-1)
					if (iy>=1) bbtmp1=ep2*iy*sfty**(iy-1)
					if (iz>=1) cctmp1=ep2*iz*sftz**(iz-1)
					aatmp2=0D0
					bbtmp2=0D0
					cctmp2=0D0
					if (ix>=2) aatmp2=ep2*ix*(ix-1)*sftx**(ix-1)
					if (iy>=2) bbtmp2=ep2*iy*(iy-1)*sfty**(iy-1)
					if (iz>=2) cctmp2=ep2*iz*(iz-1)*sftz**(iz-1)
					aatmp3=0D0
					bbtmp3=0D0
					cctmp3=0D0
					if (ix>=3) aatmp3=ix*(ix-1)*(ix-2)*sftx**(ix-3)
					if (iy>=3) bbtmp3=iy*(iy-1)*(iy-2)*sfty**(iy-3)
					if (iz>=3) cctmp3=iz*(iz-1)*(iz-2)*sftz**(iz-3)
					GTFdxxx=sfty**iy*sftz**iz*expterm*( (-2*ix+ep2*sftx2-1)*(-epep4*sftx**(ix+1) + aatmp1) - aatmp2 + epep8*sftx**(ix+1) + aatmp3 )
					GTFdyyy=sftx**ix*sftz**iz*expterm*( (-2*iy+ep2*sfty2-1)*(-epep4*sfty**(iy+1) + bbtmp1) - bbtmp2 + epep8*sfty**(iy+1) + bbtmp3 )
					GTFdzzz=sfty**iy*sftx**ix*expterm*( (-2*iz+ep2*sftz2-1)*(-epep4*sftz**(iz+1) + cctmp1) - cctmp2 + epep8*sftz**(iz+1) + cctmp3 )
				
					!Calculate 3-order derivative tensor for orbital wavefunction
					do imo=istart,iend
						tens3(1,1,1,imo)=tens3(1,1,1,imo)+CO_uniq(imo,j)*GTFdxxx !dxxx
						tens3(2,2,2,imo)=tens3(2,2,2,imo)+CO_uniq(imo,j)*GTFdyyy !dyyy
						tens3(3,3,3,imo)=tens3(3,3,3,imo)+CO_uniq(imo,j)*GTFdzzz !dzzz
						tens3(1,2,2,imo)=tens3(1,2,2,imo)+CO_uniq(imo,j)*GTFdxyy !dxyy
						tens3(1,1,2,imo)=tens3(1,1,2,imo)+CO_uniq(imo,j)*GTFdxxy !dxxy
						tens3(1,1,3,imo)=tens3(1,1,3,imo)+CO_uniq(imo,j)*GTFdxxz !dxxz
						tens3(1,3,3,imo)=tens3(1,3,3,imo)+CO_uniq(imo,j)*GTFdxzz !dxzz
						tens3(2,3,3,imo)=tens3(2,3,3,imo)+CO_uniq(imo,j)*GTFdyzz !dyzz
						tens3(2,2,3,imo)=tens3(2,2,3,imo)+CO_uniq(imo,j)*GTFdyyz !dyyz
						tens3(1,2,3,imo)=tens3(1,2,3,imo)+CO_uniq(imo,j)*GTFdxyz !dxyz
					end do
					tens3(1,2,1,:)=tens3(1,1,2,:) !dxyx=dxxy
					tens3(1,3,1,:)=tens3(1,1,3,:) !dxzx=dxxz
					tens3(1,3,2,:)=tens3(1,2,3,:) !dxzy=dxyz
					tens3(2,1,1,:)=tens3(1,1,2,:) !dyxx=dxxy
					tens3(2,1,2,:)=tens3(1,2,2,:) !dyxy=dxyy
					tens3(2,1,3,:)=tens3(1,2,3,:) !dyxz=dxyz
					tens3(2,2,1,:)=tens3(1,2,2,:) !dyyx=dxyy
					tens3(2,3,1,:)=tens3(1,2,3,:) !dyzx=dxyz
					tens3(2,3,2,:)=tens3(2,2,3,:) !dyzy=dyyz
					tens3(3,1,1,:)=tens3(1,1,3,:) !dzxx=dxxz
					tens3(3,1,2,:)=tens3(1,2,3,:) !dzxy=dxyz
					tens3(3,1,3,:)=tens3(1,3,3,:) !dzxz=dxzz
					tens3(3,2,1,:)=tens3(1,2,3,:) !dzyx=dxyz
					tens3(3,2,2,:)=tens3(2,2,3,:) !dzyy=dyyz
					tens3(3,2,3,:)=tens3(2,3,3,:) !dzyz=dyzz
					tens3(3,3,1,:)=tens3(1,3,3,:) !dzzx=dxzz
					tens3(3,3,2,:)=tens3(2,3,3,:) !dzzy=dyzz
				end if !end runtype>=5
			
			end if !end runtype>=3
		end if !end runtype>=2
	end do
end if
end subroutine




!!--------- The same as orbderv, but explicitly consider PBC
!! Calculate wavefunction value of a range of orbitals and their derivatives at a given point, up to third-order
!! istart and iend is the range of the orbitals will be calculated, to calculate all orbitals, use 1,nmo
!! runtype=1: value  =2: value+dx/y/z  =3: value+dxx/yy/zz(diagonal of hess)  =4: value+dx/y/z+Hessian  
!!        =5: value+dx/y/z+hess+3-order derivative tensor
!! k-point is taken into account if kp1crd,kp2crd,kp3crd has been properly defined (default to 0,0,0), in this case, &
!! itype_in=1: Returns real part  =2: Returns imaginary part (e.g. 4.3i returns 4.3). If not explicitly defined, it is decided by global variable "iorbcomplex"
subroutine orbderv_PBC(runtype,istart,iend,x,y,z,wfnval,grad,hess,tens3,itype_in)
real*8 x,y,z,wfnval(nmo),tvec(3)
real*8,optional :: grad(3,nmo),hess(3,3,nmo),tens3(3,3,3,nmo)
integer runtype,istart,iend
integer,optional :: itype_in

wfnval=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0
if (present(itype_in)) then
	itype=itype_in
else
	itype=iorbcomplex
end if

iquick=0 !Check for this position if quick code can be used. If this position is out of region covered by reduced grid, the slow conventional code must be used
if (allocated(neighGTF)) then
	ix_red=floor((x-orgx_neigh)/spcred)
	iy_red=floor((y-orgy_neigh)/spcred)
	iz_red=floor((z-orgz_neigh)/spcred)
    if (ix_red>=0.and.ix_red<=size(neighnGTF,1)-1 .and. iy_red>=0.and.iy_red<=size(neighnGTF,2)-1 .and. iz_red>=0.and.iz_red<=size(neighnGTF,3)-1) iquick=1
end if

if (iquick==1) then !Utilizing neighbouring GTF list at reduced grids to significantly reduce number of candidate GTFs. Unique GTF cannot be utilized in this case

nloopGTF=neighnGTF(ix_red,iy_red,iz_red)
do iloopGTF=1,nloopGTF !Only loop the neighbouring GTFs corresponding to the reduced grid that the current position belongs to
	j=neighGTF(iloopGTF,ix_red,iy_red,iz_red)
    icell=neighGTFcell(1,iloopGTF,ix_red,iy_red,iz_red) !Obtaining the cell index that the neighbouring GTF attributed to
    jcell=neighGTFcell(2,iloopGTF,ix_red,iy_red,iz_red)
    kcell=neighGTFcell(3,iloopGTF,ix_red,iy_red,iz_red)
    call tvec_PBC(icell,jcell,kcell,tvec)
    if (itype==1) then !Real part of exponent term of wavefunction
		fac=cos(2*pi*(icell*kp1crd+jcell*kp2crd+kcell*kp3crd))
    else !Imaginary part of exponent term of wavefunction
		fac=sin(2*pi*(icell*kp1crd+jcell*kp2crd+kcell*kp3crd))
    end if
	jcen=b(j)%center
	sftx=x-(a(jcen)%x+tvec(1))
	sfty=y-(a(jcen)%y+tvec(2))
	sftz=z-(a(jcen)%z+tvec(3))
    sftx2=sftx*sftx
    sfty2=sfty*sfty
    sftz2=sftz*sftz
	rr=sftx2 + sfty2 + sftz2
	ep=b(j)%exp
	tmpval=-ep*rr
	if (tmpval>expcutoff_PBC.or.expcutoff_PBC>0) then
		expterm=exp(tmpval)*fac !!!! Exponent term of orbital wavefunction is merged into GTF exponent term to make it take effect properly without modifying any other codes
    else
		cycle
    end if
	
	!Calculate value for current GTF
	jtype=b(j)%type
	ix=type2ix(jtype)
	iy=type2iy(jtype)
	iz=type2iz(jtype)
	if (jtype==1) then
	GTFval=expterm
	else if (jtype==2) then
	GTFval=sftx*expterm
	else if (jtype==3) then
	GTFval=sfty*expterm
	else if (jtype==4) then
	GTFval=sftz*expterm
	else if (jtype==5) then
	GTFval=sftx2*expterm
	else if (jtype==6) then
	GTFval=sfty2*expterm
	else if (jtype==7) then
	GTFval=sftz2*expterm
	else if (jtype==8) then
	GTFval=sftx*sfty*expterm
	else if (jtype==9) then
	GTFval=sftx*sftz*expterm
	else if (jtype==10) then
	GTFval=sfty*sftz*expterm
	else
	GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
	end if
	!Calculate orbital wavefunction value
	do imo=istart,iend
		wfnval(imo)=wfnval(imo)+CO(imo,j)*GTFval
	end do
                
	if (runtype>=2) then
		!Calculate 1-order derivative for current GTF
		tx=0D0
		ty=0D0
		tz=0D0
		if (ix/=0) tx=ix*sftx**(ix-1)
		if (iy/=0) ty=iy*sfty**(iy-1)
		if (iz/=0) tz=iz*sftz**(iz-1)
		GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
		GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
		GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
		!Calculate 1-order derivative for orbitals
		do imo=istart,iend
			grad(1,imo)=grad(1,imo)+CO(imo,j)*GTFdx
			grad(2,imo)=grad(2,imo)+CO(imo,j)*GTFdy
			grad(3,imo)=grad(3,imo)+CO(imo,j)*GTFdz
		end do

		if (runtype>=3) then
			!Calculate 2-order derivative for current GTF
			txx=0D0
			tyy=0D0
			tzz=0D0
			if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
			if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
			if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
			GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
			GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
			GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
			ttx=tx-2*ep*sftx**(ix+1)
			tty=ty-2*ep*sfty**(iy+1)
			ttz=tz-2*ep*sftz**(iz+1)
			GTFdxy=sftz**iz *expterm*ttx*tty
			GTFdyz=sftx**ix *expterm*tty*ttz
			GTFdxz=sfty**iy *expterm*ttx*ttz
			!Calculate diagonal Hessian elements for orbitals
			do imo=istart,iend
				hess(1,1,imo)=hess(1,1,imo)+CO(imo,j)*GTFdxx !dxx
				hess(2,2,imo)=hess(2,2,imo)+CO(imo,j)*GTFdyy !dyy
				hess(3,3,imo)=hess(3,3,imo)+CO(imo,j)*GTFdzz !dzz
			end do
			if (runtype>=4) then !Also process nondiagonal elements
				do imo=istart,iend
					hess(1,2,imo)=hess(1,2,imo)+CO(imo,j)*GTFdxy !dxy
					hess(2,3,imo)=hess(2,3,imo)+CO(imo,j)*GTFdyz !dyz
					hess(1,3,imo)=hess(1,3,imo)+CO(imo,j)*GTFdxz !dxz
				end do
				hess(2,1,:)=hess(1,2,:)
				hess(3,2,:)=hess(2,3,:)
				hess(3,1,:)=hess(1,3,:)
			end if
			
			if (runtype>=5) then
				!Calculate 3-order derivative for current GTF
				ep2=ep*2D0
				ep4=ep*4D0
				epep4=ep2*ep2
				epep8=epep4*2D0
				!dxyz
				a1=0D0
				b1=0D0
				c1=0D0
				if (ix>=1) a1=ix*sftx**(ix-1)
				if (iy>=1) b1=iy*sfty**(iy-1)
				if (iz>=1) c1=iz*sftz**(iz-1)
				a2=-ep2*sftx**(ix+1)
				b2=-ep2*sfty**(iy+1)
				c2=-ep2*sftz**(iz+1)
				GTFdxyz=(a1+a2)*(b1+b2)*(c1+c2)*expterm
				!dxyy,dxxy,dxxz,dxzz,dyzz,dyyz
				atmp=0D0
				btmp=0D0
				ctmp=0D0
				if (ix>=2) atmp=ix*(ix-1)*sftx**(ix-2)
				if (iy>=2) btmp=iy*(iy-1)*sfty**(iy-2)
				if (iz>=2) ctmp=iz*(iz-1)*sftz**(iz-2)
				GTFdxyy=(a1+a2)*sftz**iz *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy)
				GTFdxxy=(b1+b2)*sftz**iz *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dyxx
				GTFdxxz=(c1+c2)*sfty**iy *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dzxx
				GTFdxzz=(a1+a2)*sfty**iy *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
				GTFdyzz=(b1+b2)*sftx**ix *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
				GTFdyyz=(c1+c2)*sftx**ix *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !=dzyy
				!dxxx,dyyy,dzzz
				aatmp1=0D0
				bbtmp1=0D0
				cctmp1=0D0
				if (ix>=1) aatmp1=ep2*ix*sftx**(ix-1)
				if (iy>=1) bbtmp1=ep2*iy*sfty**(iy-1)
				if (iz>=1) cctmp1=ep2*iz*sftz**(iz-1)
				aatmp2=0D0
				bbtmp2=0D0
				cctmp2=0D0
				if (ix>=2) aatmp2=ep2*ix*(ix-1)*sftx**(ix-1)
				if (iy>=2) bbtmp2=ep2*iy*(iy-1)*sfty**(iy-1)
				if (iz>=2) cctmp2=ep2*iz*(iz-1)*sftz**(iz-1)
				aatmp3=0D0
				bbtmp3=0D0
				cctmp3=0D0
				if (ix>=3) aatmp3=ix*(ix-1)*(ix-2)*sftx**(ix-3)
				if (iy>=3) bbtmp3=iy*(iy-1)*(iy-2)*sfty**(iy-3)
				if (iz>=3) cctmp3=iz*(iz-1)*(iz-2)*sftz**(iz-3)
				GTFdxxx=sfty**iy*sftz**iz*expterm*( (-2*ix+ep2*sftx2-1)*(-epep4*sftx**(ix+1) + aatmp1) - aatmp2 + epep8*sftx**(ix+1) + aatmp3 )
				GTFdyyy=sftx**ix*sftz**iz*expterm*( (-2*iy+ep2*sfty2-1)*(-epep4*sfty**(iy+1) + bbtmp1) - bbtmp2 + epep8*sfty**(iy+1) + bbtmp3 )
				GTFdzzz=sfty**iy*sftx**ix*expterm*( (-2*iz+ep2*sftz2-1)*(-epep4*sftz**(iz+1) + cctmp1) - cctmp2 + epep8*sftz**(iz+1) + cctmp3 )
				
				!Calculate 3-order derivative tensor for orbital wavefunction
				do imo=istart,iend
					tens3(1,1,1,imo)=tens3(1,1,1,imo)+CO(imo,j)*GTFdxxx !dxxx
					tens3(2,2,2,imo)=tens3(2,2,2,imo)+CO(imo,j)*GTFdyyy !dyyy
					tens3(3,3,3,imo)=tens3(3,3,3,imo)+CO(imo,j)*GTFdzzz !dzzz
					tens3(1,2,2,imo)=tens3(1,2,2,imo)+CO(imo,j)*GTFdxyy !dxyy
					tens3(1,1,2,imo)=tens3(1,1,2,imo)+CO(imo,j)*GTFdxxy !dxxy
					tens3(1,1,3,imo)=tens3(1,1,3,imo)+CO(imo,j)*GTFdxxz !dxxz
					tens3(1,3,3,imo)=tens3(1,3,3,imo)+CO(imo,j)*GTFdxzz !dxzz
					tens3(2,3,3,imo)=tens3(2,3,3,imo)+CO(imo,j)*GTFdyzz !dyzz
					tens3(2,2,3,imo)=tens3(2,2,3,imo)+CO(imo,j)*GTFdyyz !dyyz
					tens3(1,2,3,imo)=tens3(1,2,3,imo)+CO(imo,j)*GTFdxyz !dxyz
				end do
				tens3(1,2,1,:)=tens3(1,1,2,:) !dxyx=dxxy
				tens3(1,3,1,:)=tens3(1,1,3,:) !dxzx=dxxz
				tens3(1,3,2,:)=tens3(1,2,3,:) !dxzy=dxyz
				tens3(2,1,1,:)=tens3(1,1,2,:) !dyxx=dxxy
				tens3(2,1,2,:)=tens3(1,2,2,:) !dyxy=dxyy
				tens3(2,1,3,:)=tens3(1,2,3,:) !dyxz=dxyz
				tens3(2,2,1,:)=tens3(1,2,2,:) !dyyx=dxyy
				tens3(2,3,1,:)=tens3(1,2,3,:) !dyzx=dxyz
				tens3(2,3,2,:)=tens3(2,2,3,:) !dyzy=dyyz
				tens3(3,1,1,:)=tens3(1,1,3,:) !dzxx=dxxz
				tens3(3,1,2,:)=tens3(1,2,3,:) !dzxy=dxyz
				tens3(3,1,3,:)=tens3(1,3,3,:) !dzxz=dxzz
				tens3(3,2,1,:)=tens3(1,2,3,:) !dzyx=dxyz
				tens3(3,2,2,:)=tens3(2,2,3,:) !dzyy=dyyz
				tens3(3,2,3,:)=tens3(2,3,3,:) !dzyz=dyzz
				tens3(3,3,1,:)=tens3(1,3,3,:) !dzzx=dxzz
				tens3(3,3,2,:)=tens3(2,3,3,:) !dzzy=dyzz
			end if !end runtype>=5
			
		end if !end runtype>=3
	end if !end runtype>=2
end do

else !Using very slow method: Looping all candidate GTFs in all possible cells. Unique GTF can be utilized

call getpointcell(x,y,z,ic,jc,kc)
do icell=ic-PBCnx,ic+PBCnx
    do jcell=jc-PBCny,jc+PBCny
        do kcell=kc-PBCnz,kc+PBCnz
			if (itype==1) then !Real part of exponent term of wavefunction
				fac=cos(2*pi*(icell*kp1crd+jcell*kp2crd+kcell*kp3crd))
			else !Imaginary part of exponent term of wavefunction
				fac=sin(2*pi*(icell*kp1crd+jcell*kp2crd+kcell*kp3crd))
			end if
			lastcen=-1 !Arbitrary value
            call tvec_PBC(icell,jcell,kcell,tvec)
            xmove=tvec(1)
            ymove=tvec(2)
            zmove=tvec(3)
            
            if (nprims_uniq==0) then !Unique GTF information is not available
				do j=1,nprims
					jcen=b(j)%center
					if (jcen/=lastcen) then
						sftx=x-(a(jcen)%x+xmove)
						sfty=y-(a(jcen)%y+ymove)
						sftz=z-(a(jcen)%z+zmove)
	            		sftx2=sftx*sftx
	            		sfty2=sfty*sfty
	            		sftz2=sftz*sftz
	            		rr=sftx2+sfty2+sftz2
					end if
					ep=b(j)%exp
					tmpval=-ep*rr
					lastcen=jcen
					if (tmpval>expcutoff_PBC.or.expcutoff_PBC>0) then
                        expterm=exp(tmpval)*fac !!!! Exponent term of orbital wavefunction is merged into GTF exponent term to make it take effect properly without modifying any other codes
					else
						cycle
					end if
	
					!Calculate value for current GTF
					jtype=b(j)%type
					ix=type2ix(jtype)
					iy=type2iy(jtype)
					iz=type2iz(jtype)
					if (jtype==1) then
					GTFval=expterm
					else if (jtype==2) then
					GTFval=sftx*expterm
					else if (jtype==3) then
					GTFval=sfty*expterm
					else if (jtype==4) then
					GTFval=sftz*expterm
					else if (jtype==5) then
					GTFval=sftx2*expterm
					else if (jtype==6) then
					GTFval=sfty2*expterm
					else if (jtype==7) then
					GTFval=sftz2*expterm
					else if (jtype==8) then
					GTFval=sftx*sfty*expterm
					else if (jtype==9) then
					GTFval=sftx*sftz*expterm
					else if (jtype==10) then
					GTFval=sfty*sftz*expterm
					else
					GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
					end if
					!Calculate orbital wavefunction value
					do imo=istart,iend
						wfnval(imo)=wfnval(imo)+CO(imo,j)*GTFval
					end do
                
					if (runtype>=2) then
						!Calculate 1-order derivative for current GTF
						tx=0D0
						ty=0D0
						tz=0D0
						if (ix/=0) tx=ix*sftx**(ix-1)
						if (iy/=0) ty=iy*sfty**(iy-1)
						if (iz/=0) tz=iz*sftz**(iz-1)
						GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
						GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
						GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
						!Calculate 1-order derivative for orbitals
						do imo=istart,iend
							grad(1,imo)=grad(1,imo)+CO(imo,j)*GTFdx
							grad(2,imo)=grad(2,imo)+CO(imo,j)*GTFdy
							grad(3,imo)=grad(3,imo)+CO(imo,j)*GTFdz
						end do

						if (runtype>=3) then
							!Calculate 2-order derivative for current GTF
							txx=0D0
							tyy=0D0
							tzz=0D0
							if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
							if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
							if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
							GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
							GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
							GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
							ttx=tx-2*ep*sftx**(ix+1)
							tty=ty-2*ep*sfty**(iy+1)
							ttz=tz-2*ep*sftz**(iz+1)
							GTFdxy=sftz**iz *expterm*ttx*tty
							GTFdyz=sftx**ix *expterm*tty*ttz
							GTFdxz=sfty**iy *expterm*ttx*ttz
							!Calculate diagonal Hessian elements for orbitals
							do imo=istart,iend
								hess(1,1,imo)=hess(1,1,imo)+CO(imo,j)*GTFdxx !dxx
								hess(2,2,imo)=hess(2,2,imo)+CO(imo,j)*GTFdyy !dyy
								hess(3,3,imo)=hess(3,3,imo)+CO(imo,j)*GTFdzz !dzz
							end do
							if (runtype>=4) then !Also process nondiagonal elements
								do imo=istart,iend
									hess(1,2,imo)=hess(1,2,imo)+CO(imo,j)*GTFdxy !dxy
									hess(2,3,imo)=hess(2,3,imo)+CO(imo,j)*GTFdyz !dyz
									hess(1,3,imo)=hess(1,3,imo)+CO(imo,j)*GTFdxz !dxz
								end do
								hess(2,1,:)=hess(1,2,:)
								hess(3,2,:)=hess(2,3,:)
								hess(3,1,:)=hess(1,3,:)
							end if
			
							if (runtype>=5) then
								!Calculate 3-order derivative for current GTF
								ep2=ep*2D0
								ep4=ep*4D0
								epep4=ep2*ep2
								epep8=epep4*2D0
								!dxyz
								a1=0D0
								b1=0D0
								c1=0D0
								if (ix>=1) a1=ix*sftx**(ix-1)
								if (iy>=1) b1=iy*sfty**(iy-1)
								if (iz>=1) c1=iz*sftz**(iz-1)
								a2=-ep2*sftx**(ix+1)
								b2=-ep2*sfty**(iy+1)
								c2=-ep2*sftz**(iz+1)
								GTFdxyz=(a1+a2)*(b1+b2)*(c1+c2)*expterm
								!dxyy,dxxy,dxxz,dxzz,dyzz,dyyz
								atmp=0D0
								btmp=0D0
								ctmp=0D0
								if (ix>=2) atmp=ix*(ix-1)*sftx**(ix-2)
								if (iy>=2) btmp=iy*(iy-1)*sfty**(iy-2)
								if (iz>=2) ctmp=iz*(iz-1)*sftz**(iz-2)
								GTFdxyy=(a1+a2)*sftz**iz *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy)
								GTFdxxy=(b1+b2)*sftz**iz *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dyxx
								GTFdxxz=(c1+c2)*sfty**iy *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dzxx
								GTFdxzz=(a1+a2)*sfty**iy *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
								GTFdyzz=(b1+b2)*sftx**ix *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
								GTFdyyz=(c1+c2)*sftx**ix *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !=dzyy
								!dxxx,dyyy,dzzz
								aatmp1=0D0
								bbtmp1=0D0
								cctmp1=0D0
								if (ix>=1) aatmp1=ep2*ix*sftx**(ix-1)
								if (iy>=1) bbtmp1=ep2*iy*sfty**(iy-1)
								if (iz>=1) cctmp1=ep2*iz*sftz**(iz-1)
								aatmp2=0D0
								bbtmp2=0D0
								cctmp2=0D0
								if (ix>=2) aatmp2=ep2*ix*(ix-1)*sftx**(ix-1)
								if (iy>=2) bbtmp2=ep2*iy*(iy-1)*sfty**(iy-1)
								if (iz>=2) cctmp2=ep2*iz*(iz-1)*sftz**(iz-1)
								aatmp3=0D0
								bbtmp3=0D0
								cctmp3=0D0
								if (ix>=3) aatmp3=ix*(ix-1)*(ix-2)*sftx**(ix-3)
								if (iy>=3) bbtmp3=iy*(iy-1)*(iy-2)*sfty**(iy-3)
								if (iz>=3) cctmp3=iz*(iz-1)*(iz-2)*sftz**(iz-3)
								GTFdxxx=sfty**iy*sftz**iz*expterm*( (-2*ix+ep2*sftx2-1)*(-epep4*sftx**(ix+1) + aatmp1) - aatmp2 + epep8*sftx**(ix+1) + aatmp3 )
								GTFdyyy=sftx**ix*sftz**iz*expterm*( (-2*iy+ep2*sfty2-1)*(-epep4*sfty**(iy+1) + bbtmp1) - bbtmp2 + epep8*sfty**(iy+1) + bbtmp3 )
								GTFdzzz=sfty**iy*sftx**ix*expterm*( (-2*iz+ep2*sftz2-1)*(-epep4*sftz**(iz+1) + cctmp1) - cctmp2 + epep8*sftz**(iz+1) + cctmp3 )
				
								!Calculate 3-order derivative tensor for orbital wavefunction
								do imo=istart,iend
									tens3(1,1,1,imo)=tens3(1,1,1,imo)+CO(imo,j)*GTFdxxx !dxxx
									tens3(2,2,2,imo)=tens3(2,2,2,imo)+CO(imo,j)*GTFdyyy !dyyy
									tens3(3,3,3,imo)=tens3(3,3,3,imo)+CO(imo,j)*GTFdzzz !dzzz
									tens3(1,2,2,imo)=tens3(1,2,2,imo)+CO(imo,j)*GTFdxyy !dxyy
									tens3(1,1,2,imo)=tens3(1,1,2,imo)+CO(imo,j)*GTFdxxy !dxxy
									tens3(1,1,3,imo)=tens3(1,1,3,imo)+CO(imo,j)*GTFdxxz !dxxz
									tens3(1,3,3,imo)=tens3(1,3,3,imo)+CO(imo,j)*GTFdxzz !dxzz
									tens3(2,3,3,imo)=tens3(2,3,3,imo)+CO(imo,j)*GTFdyzz !dyzz
									tens3(2,2,3,imo)=tens3(2,2,3,imo)+CO(imo,j)*GTFdyyz !dyyz
									tens3(1,2,3,imo)=tens3(1,2,3,imo)+CO(imo,j)*GTFdxyz !dxyz
								end do
								tens3(1,2,1,:)=tens3(1,1,2,:) !dxyx=dxxy
								tens3(1,3,1,:)=tens3(1,1,3,:) !dxzx=dxxz
								tens3(1,3,2,:)=tens3(1,2,3,:) !dxzy=dxyz
								tens3(2,1,1,:)=tens3(1,1,2,:) !dyxx=dxxy
								tens3(2,1,2,:)=tens3(1,2,2,:) !dyxy=dxyy
								tens3(2,1,3,:)=tens3(1,2,3,:) !dyxz=dxyz
								tens3(2,2,1,:)=tens3(1,2,2,:) !dyyx=dxyy
								tens3(2,3,1,:)=tens3(1,2,3,:) !dyzx=dxyz
								tens3(2,3,2,:)=tens3(2,2,3,:) !dyzy=dyyz
								tens3(3,1,1,:)=tens3(1,1,3,:) !dzxx=dxxz
								tens3(3,1,2,:)=tens3(1,2,3,:) !dzxy=dxyz
								tens3(3,1,3,:)=tens3(1,3,3,:) !dzxz=dxzz
								tens3(3,2,1,:)=tens3(1,2,3,:) !dzyx=dxyz
								tens3(3,2,2,:)=tens3(2,2,3,:) !dzyy=dyyz
								tens3(3,2,3,:)=tens3(2,3,3,:) !dzyz=dyzz
								tens3(3,3,1,:)=tens3(1,3,3,:) !dzzx=dxzz
								tens3(3,3,2,:)=tens3(2,3,3,:) !dzzy=dyzz
							end if !end runtype>=5
			
						end if !end runtype>=3
					end if !end runtype>=2
				end do
                
            else !Unique GTF information has been constructed by gen_GTFuniq
				do j=1,nprims_uniq
					jcen=b_uniq(j)%center
					if (jcen/=lastcen) then
						sftx=x-(a(jcen)%x+xmove)
						sfty=y-(a(jcen)%y+ymove)
						sftz=z-(a(jcen)%z+zmove)
	            		sftx2=sftx*sftx
	            		sfty2=sfty*sfty
	            		sftz2=sftz*sftz
	            		rr=sftx2+sfty2+sftz2
					end if
					ep=b_uniq(j)%exp
					tmpval=-ep*rr
					lastcen=jcen
					if (tmpval>expcutoff_PBC.or.expcutoff_PBC>0) then
                        expterm=exp(tmpval)
					else
						cycle
					end if
	
					!Calculate value for current GTF
					jtype=b_uniq(j)%type
					ix=type2ix(jtype)
					iy=type2iy(jtype)
					iz=type2iz(jtype)
					if (jtype==1) then
					GTFval=expterm
					else if (jtype==2) then
					GTFval=sftx*expterm
					else if (jtype==3) then
					GTFval=sfty*expterm
					else if (jtype==4) then
					GTFval=sftz*expterm
					else if (jtype==5) then
					GTFval=sftx2*expterm
					else if (jtype==6) then
					GTFval=sfty2*expterm
					else if (jtype==7) then
					GTFval=sftz2*expterm
					else if (jtype==8) then
					GTFval=sftx*sfty*expterm
					else if (jtype==9) then
					GTFval=sftx*sftz*expterm
					else if (jtype==10) then
					GTFval=sfty*sftz*expterm
					else
					GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
					end if
					!Calculate orbital wavefunction value
					do imo=istart,iend
						wfnval(imo)=wfnval(imo)+CO_uniq(imo,j)*GTFval
					end do
                
					if (runtype>=2) then
						!Calculate 1-order derivative for current GTF
						tx=0D0
						ty=0D0
						tz=0D0
						if (ix/=0) tx=ix*sftx**(ix-1)
						if (iy/=0) ty=iy*sfty**(iy-1)
						if (iz/=0) tz=iz*sftz**(iz-1)
						GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
						GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
						GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
						!Calculate 1-order derivative for orbitals
						do imo=istart,iend
							grad(1,imo)=grad(1,imo)+CO_uniq(imo,j)*GTFdx
							grad(2,imo)=grad(2,imo)+CO_uniq(imo,j)*GTFdy
							grad(3,imo)=grad(3,imo)+CO_uniq(imo,j)*GTFdz
						end do

						if (runtype>=3) then
							!Calculate 2-order derivative for current GTF
							txx=0D0
							tyy=0D0
							tzz=0D0
							if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
							if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
							if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
							GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
							GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
							GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
							ttx=tx-2*ep*sftx**(ix+1)
							tty=ty-2*ep*sfty**(iy+1)
							ttz=tz-2*ep*sftz**(iz+1)
							GTFdxy=sftz**iz *expterm*ttx*tty
							GTFdyz=sftx**ix *expterm*tty*ttz
							GTFdxz=sfty**iy *expterm*ttx*ttz
							!Calculate diagonal Hessian elements for orbitals
							do imo=istart,iend
								hess(1,1,imo)=hess(1,1,imo)+CO_uniq(imo,j)*GTFdxx !dxx
								hess(2,2,imo)=hess(2,2,imo)+CO_uniq(imo,j)*GTFdyy !dyy
								hess(3,3,imo)=hess(3,3,imo)+CO_uniq(imo,j)*GTFdzz !dzz
							end do
							if (runtype>=4) then !Also process nondiagonal elements
								do imo=istart,iend
									hess(1,2,imo)=hess(1,2,imo)+CO_uniq(imo,j)*GTFdxy !dxy
									hess(2,3,imo)=hess(2,3,imo)+CO_uniq(imo,j)*GTFdyz !dyz
									hess(1,3,imo)=hess(1,3,imo)+CO_uniq(imo,j)*GTFdxz !dxz
								end do
								hess(2,1,:)=hess(1,2,:)
								hess(3,2,:)=hess(2,3,:)
								hess(3,1,:)=hess(1,3,:)
							end if
			
							if (runtype>=5) then
								!Calculate 3-order derivative for current GTF
								ep2=ep*2D0
								ep4=ep*4D0
								epep4=ep2*ep2
								epep8=epep4*2D0
								!dxyz
								a1=0D0
								b1=0D0
								c1=0D0
								if (ix>=1) a1=ix*sftx**(ix-1)
								if (iy>=1) b1=iy*sfty**(iy-1)
								if (iz>=1) c1=iz*sftz**(iz-1)
								a2=-ep2*sftx**(ix+1)
								b2=-ep2*sfty**(iy+1)
								c2=-ep2*sftz**(iz+1)
								GTFdxyz=(a1+a2)*(b1+b2)*(c1+c2)*expterm
								!dxyy,dxxy,dxxz,dxzz,dyzz,dyyz
								atmp=0D0
								btmp=0D0
								ctmp=0D0
								if (ix>=2) atmp=ix*(ix-1)*sftx**(ix-2)
								if (iy>=2) btmp=iy*(iy-1)*sfty**(iy-2)
								if (iz>=2) ctmp=iz*(iz-1)*sftz**(iz-2)
								GTFdxyy=(a1+a2)*sftz**iz *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy)
								GTFdxxy=(b1+b2)*sftz**iz *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dyxx
								GTFdxxz=(c1+c2)*sfty**iy *expterm*(-ep4*ix*sftx**ix+epep4*sftx**(ix+2)+atmp-ep2*sftx**ix) !=dzxx
								GTFdxzz=(a1+a2)*sfty**iy *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
								GTFdyzz=(b1+b2)*sftx**ix *expterm*(-ep4*iz*sftz**iz+epep4*sftz**(iz+2)+ctmp-ep2*sftz**iz)
								GTFdyyz=(c1+c2)*sftx**ix *expterm*(-ep4*iy*sfty**iy+epep4*sfty**(iy+2)+btmp-ep2*sfty**iy) !=dzyy
								!dxxx,dyyy,dzzz
								aatmp1=0D0
								bbtmp1=0D0
								cctmp1=0D0
								if (ix>=1) aatmp1=ep2*ix*sftx**(ix-1)
								if (iy>=1) bbtmp1=ep2*iy*sfty**(iy-1)
								if (iz>=1) cctmp1=ep2*iz*sftz**(iz-1)
								aatmp2=0D0
								bbtmp2=0D0
								cctmp2=0D0
								if (ix>=2) aatmp2=ep2*ix*(ix-1)*sftx**(ix-1)
								if (iy>=2) bbtmp2=ep2*iy*(iy-1)*sfty**(iy-1)
								if (iz>=2) cctmp2=ep2*iz*(iz-1)*sftz**(iz-1)
								aatmp3=0D0
								bbtmp3=0D0
								cctmp3=0D0
								if (ix>=3) aatmp3=ix*(ix-1)*(ix-2)*sftx**(ix-3)
								if (iy>=3) bbtmp3=iy*(iy-1)*(iy-2)*sfty**(iy-3)
								if (iz>=3) cctmp3=iz*(iz-1)*(iz-2)*sftz**(iz-3)
								GTFdxxx=sfty**iy*sftz**iz*expterm*( (-2*ix+ep2*sftx2-1)*(-epep4*sftx**(ix+1) + aatmp1) - aatmp2 + epep8*sftx**(ix+1) + aatmp3 )
								GTFdyyy=sftx**ix*sftz**iz*expterm*( (-2*iy+ep2*sfty2-1)*(-epep4*sfty**(iy+1) + bbtmp1) - bbtmp2 + epep8*sfty**(iy+1) + bbtmp3 )
								GTFdzzz=sfty**iy*sftx**ix*expterm*( (-2*iz+ep2*sftz2-1)*(-epep4*sftz**(iz+1) + cctmp1) - cctmp2 + epep8*sftz**(iz+1) + cctmp3 )
				
								!Calculate 3-order derivative tensor for orbital wavefunction
								do imo=istart,iend
									tens3(1,1,1,imo)=tens3(1,1,1,imo)+CO_uniq(imo,j)*GTFdxxx !dxxx
									tens3(2,2,2,imo)=tens3(2,2,2,imo)+CO_uniq(imo,j)*GTFdyyy !dyyy
									tens3(3,3,3,imo)=tens3(3,3,3,imo)+CO_uniq(imo,j)*GTFdzzz !dzzz
									tens3(1,2,2,imo)=tens3(1,2,2,imo)+CO_uniq(imo,j)*GTFdxyy !dxyy
									tens3(1,1,2,imo)=tens3(1,1,2,imo)+CO_uniq(imo,j)*GTFdxxy !dxxy
									tens3(1,1,3,imo)=tens3(1,1,3,imo)+CO_uniq(imo,j)*GTFdxxz !dxxz
									tens3(1,3,3,imo)=tens3(1,3,3,imo)+CO_uniq(imo,j)*GTFdxzz !dxzz
									tens3(2,3,3,imo)=tens3(2,3,3,imo)+CO_uniq(imo,j)*GTFdyzz !dyzz
									tens3(2,2,3,imo)=tens3(2,2,3,imo)+CO_uniq(imo,j)*GTFdyyz !dyyz
									tens3(1,2,3,imo)=tens3(1,2,3,imo)+CO_uniq(imo,j)*GTFdxyz !dxyz
								end do
								tens3(1,2,1,:)=tens3(1,1,2,:) !dxyx=dxxy
								tens3(1,3,1,:)=tens3(1,1,3,:) !dxzx=dxxz
								tens3(1,3,2,:)=tens3(1,2,3,:) !dxzy=dxyz
								tens3(2,1,1,:)=tens3(1,1,2,:) !dyxx=dxxy
								tens3(2,1,2,:)=tens3(1,2,2,:) !dyxy=dxyy
								tens3(2,1,3,:)=tens3(1,2,3,:) !dyxz=dxyz
								tens3(2,2,1,:)=tens3(1,2,2,:) !dyyx=dxyy
								tens3(2,3,1,:)=tens3(1,2,3,:) !dyzx=dxyz
								tens3(2,3,2,:)=tens3(2,2,3,:) !dyzy=dyyz
								tens3(3,1,1,:)=tens3(1,1,3,:) !dzxx=dxxz
								tens3(3,1,2,:)=tens3(1,2,3,:) !dzxy=dxyz
								tens3(3,1,3,:)=tens3(1,3,3,:) !dzxz=dxzz
								tens3(3,2,1,:)=tens3(1,2,3,:) !dzyx=dxyz
								tens3(3,2,2,:)=tens3(2,2,3,:) !dzyy=dyyz
								tens3(3,2,3,:)=tens3(2,3,3,:) !dzyz=dyzz
								tens3(3,3,1,:)=tens3(1,3,3,:) !dzzx=dxzz
								tens3(3,3,2,:)=tens3(2,3,3,:) !dzzy=dyzz
							end if !end runtype>=5
			
						end if !end runtype>=3
					end if !end runtype>=2
				end do
            end if
            
        end do
	end do
end do

end if

end subroutine




!!--------- The same as subroutine orbderv, but calculate for promolecular wavefunction (CO_pmol ...), up to Hessian
!! istart and iend is the range of the orbitals will be calculated, to calculate all orbitals, use 1,nmo_pmol
!! runtype=1: value  =2: value+dx/y/z  =3: value+dxx/yy/zz(diagonal of hess)  =4: value+dx/y/z+Hessian
subroutine orbderv_pmol(runtype,istart,iend,x,y,z,wfnval,grad,hess)
real*8 x,y,z,wfnval(nmo_pmol)
real*8,optional :: grad(3,nmo_pmol),hess(3,3,nmo_pmol)
integer runtype,istart,iend

wfnval=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
lastcen=-1 !Arbitrary value

!If the center/exp of current GTF is the same as previous, then we do not need to recalculate them
do j=1,nprims
	ix=type2ix(b(j)%type)
	iy=type2iy(b(j)%type)
	iz=type2iz(b(j)%type)
	ep=b(j)%exp
	
	if (b(j)%center/=lastcen) then
		sftx=x-a(b(j)%center)%x
		sfty=y-a(b(j)%center)%y
		sftz=z-a(b(j)%center)%z
		sftx2=sftx*sftx
		sfty2=sfty*sfty
		sftz2=sftz*sftz
		rr=sftx2+sfty2+sftz2
	end if
	if (expcutoff>0.or.-ep*rr>expcutoff) then
		expterm=exp(-ep*rr)
	else
		expterm=0D0
	end if
	lastcen=b(j)%center
	if (expterm==0D0) cycle
	
	!Calculate value for current GTF
	GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
	!Calculate orbital wavefunction value
	do imo=istart,iend
		wfnval(imo)=wfnval(imo)+CO_pmol(imo,j)*GTFval
	end do
	
	if (runtype>=2) then
		!Calculate 1-order derivative for current GTF
		tx=0D0
		ty=0D0
		tz=0D0
		if (ix/=0) tx=ix*sftx**(ix-1)
		if (iy/=0) ty=iy*sfty**(iy-1)
		if (iz/=0) tz=iz*sftz**(iz-1)
		GTFdx=sfty**iy *sftz**iz *expterm*(tx-2*ep*sftx**(ix+1))
		GTFdy=sftx**ix *sftz**iz *expterm*(ty-2*ep*sfty**(iy+1))
		GTFdz=sftx**ix *sfty**iy *expterm*(tz-2*ep*sftz**(iz+1))
		!Calculate 1-order derivative for orbitals
		do imo=istart,iend
			grad(1,imo)=grad(1,imo)+CO_pmol(imo,j)*GTFdx
			grad(2,imo)=grad(2,imo)+CO_pmol(imo,j)*GTFdy
			grad(3,imo)=grad(3,imo)+CO_pmol(imo,j)*GTFdz
		end do

		if (runtype>=3) then
			!Calculate 2-order derivative for current GTF
			txx=0D0
			tyy=0D0
			tzz=0D0
			if (ix>=2) txx=ix*(ix-1)*sftx**(ix-2)
			if (iy>=2) tyy=iy*(iy-1)*sfty**(iy-2)
			if (iz>=2) tzz=iz*(iz-1)*sftz**(iz-2)
			GTFdxx=sfty**iy *sftz**iz *expterm*( txx + 2*ep*sftx**ix*(-2*ix+2*ep*sftx2-1) )
			GTFdyy=sftx**ix *sftz**iz *expterm*( tyy + 2*ep*sfty**iy*(-2*iy+2*ep*sfty2-1) )
			GTFdzz=sftx**ix *sfty**iy *expterm*( tzz + 2*ep*sftz**iz*(-2*iz+2*ep*sftz2-1) )
			ttx=tx-2*ep*sftx**(ix+1)
			tty=ty-2*ep*sfty**(iy+1)
			ttz=tz-2*ep*sftz**(iz+1)
			GTFdxy=sftz**iz *expterm*ttx*tty
			GTFdyz=sftx**ix *expterm*tty*ttz
			GTFdxz=sfty**iy *expterm*ttx*ttz
			!Calculate diagonal Hessian elements for orbitals
			do imo=istart,iend
				hess(1,1,imo)=hess(1,1,imo)+CO_pmol(imo,j)*GTFdxx !dxx
				hess(2,2,imo)=hess(2,2,imo)+CO_pmol(imo,j)*GTFdyy !dyy
				hess(3,3,imo)=hess(3,3,imo)+CO_pmol(imo,j)*GTFdzz !dzz
			end do
			if (runtype>=4) then !Also process nondiagonal elements
				do imo=istart,iend
					hess(1,2,imo)=hess(1,2,imo)+CO_pmol(imo,j)*GTFdxy !dxy
					hess(2,3,imo)=hess(2,3,imo)+CO_pmol(imo,j)*GTFdyz !dyz
					hess(1,3,imo)=hess(1,3,imo)+CO_pmol(imo,j)*GTFdxz !dxz
				end do
				hess(2,1,:)=hess(1,2,:)
				hess(3,2,:)=hess(2,3,:)
				hess(3,1,:)=hess(1,3,:)
			end if
		end if !end runtype>=3
	end if !end runtype>=2
end do
end subroutine




!!-------- This subroutine was contributed by i.s.ger
! istart  - first index of processing orbital
! iend    - last  index of processing orbital
! x, y, z - point's coordinates
! a0      - fractional integrals when p = 0
! a1      - fractional derivatives when p = 1
! promol  - for using promolecular wavefunction
subroutine orbfracderiv(istart, iend, x, y, z, a0, a1, promol)
  use mod_2F2, only: GTO_fractional_integral
  implicit none
  ! in variables
  integer, intent(in)  :: istart, iend
  real(8), intent(in)  :: x, y, z
  logical, intent(in), optional :: promol
  ! out variables
  real(8), intent(out), optional :: a0(nmo)
  real(8), intent(out), optional :: a1(3, nmo)
  ! internal variables
  integer :: last_center         ! index of atom in previous step
  integer :: ipmax               ! ipmax is a number of GTOs borned by differentiation (now, 1 is max)
  ! pointer's variables
  logical :: promol_
  integer :: nprims_tmp
  type(primtype), pointer :: b_p(:)
  real(8), pointer :: CO_p(:,:)
  !
  real(8) :: tvec(3)
  integer :: icell, jcell, kcell, ic, jc, kc
  !
  promol_ = .false.
  if (present(promol)) promol_ = promol
  !
  if (present(a0)) then
    a0 = 0._8
    ipmax = 0
  else if (present(a1)) then
    a1 = 0._8
    ipmax = 1
  end if
  !
  last_center = -1
  tvec = 0._8

  !calculate orbderv for promolecular wavefunction
  if (promol_) then
    nprims_tmp = nprims
    b_p => b
    CO_p => CO_pmol
  !If the center/exp of current GTF is the same as previous, then we do not need to recalculate them
  else if (nprims_uniq == 0) then !Unique GTF information is not available
    nprims_tmp = nprims
    b_p => b
    CO_p => CO
  else !Unique GTF information has been constructed by gen_GTFuniq
    nprims_tmp = nprims_uniq
    b_p => b_uniq
    CO_p => CO_uniq
  end if

  if (ifPBC > 0) then !Consider PBC
    call getpointcell(x,y,z,ic,jc,kc)
    do icell=ic-PBCnx,ic+PBCnx
      do jcell=jc-PBCny,jc+PBCny
        do kcell=kc-PBCnz,kc+PBCnz
          last_center=-1 !Arbitrary value
          call tvec_PBC(icell,jcell,kcell,tvec)
          call orbfracderiv_eval(b_p, CO_p, nprims_tmp, tvec)
        end do
      end do
    end do
  else
    call orbfracderiv_eval(b_p, CO_p, nprims_tmp, tvec)
  end if

contains
  subroutine orbfracderiv_eval(b_t, CO_t, nprims_t, tvec)
    type(primtype), intent(in) :: b_t(:)
    real(8), intent(in) :: CO_t(:,:)
    integer, intent(in) :: nprims_t
    real(8), intent(in) :: tvec(3)
    ! internal variables
    integer :: j, imo, ip
    integer :: ix, iy, iz, it      ! decomposition of azimutal quantum number of shell to x,y,z and their sum
    real(8) :: ep                  ! zeta of GTF
    real(8) :: sftx1, sfty1, sftz1 ! X^1, Y^1, Z^1  from center of atom to point
    real(8) :: sftx2, sfty2, sftz2 ! X^2, Y^2, Z^2  from center of atom to point
    real(8) :: rr                  ! R^2 from center of atom to point
    real(8) :: GTFval              ! fractional integral of GTF
    real(8) :: GTFdx, GTFdy, GTFdz ! fractional derivatives of GTF
    real(8) :: t(0:1)              ! 2F2 values (1 is a maximum of ipmax)
    !
    real(8) :: xmove, ymove, zmove

    xmove = tvec(1)
    ymove = tvec(2)
    zmove = tvec(3)

    ! main cycle over shells
    do j = 1, nprims_t
      ix = type2ix(b_t(j)%type)
      iy = type2iy(b_t(j)%type)
      iz = type2iz(b_t(j)%type)
      it = ix + iy + iz
      ep = b_t(j)%exp

      if (b_t(j)%center /= last_center) then
        last_center = b_t(j)%center
        sftx1 = x - (a(b_t(j)%center)%x + xmove)
        sfty1 = y - (a(b_t(j)%center)%y + ymove)
        sftz1 = z - (a(b_t(j)%center)%z + zmove)
        sftx2 = sftx1 * sftx1
        sfty2 = sfty1 * sfty1
        sftz2 = sftz1 * sftz1
        rr    = sftx2 + sfty2 + sftz2
      end if

! commented out since cutoff should depend on p-alpha
! Fox example, when p-alpha = 0, expcutoff can be -40 and results are fine
!          but when p-alpha = 2, expcutoff should be more than 20000 for avoiding numerical problems
!      if (expcutoff>0._8 .or. -ep*rr>-20000._8) then
      if (expcutoff>0._8 .or. .true.) then ! all orbitals are important
      else
        cycle ! just skip it
      end if

      do ip = 0, ipmax
        t(ip) = GTO_fractional_integral(it + ip * 2, -ep*rr)
      end do

      if (present(a0)) then
        ! Calculate fractional integral
        GTFval = sftx1**ix * sfty1**iy * sftz1**iz * t(0)
        ! Calculate fractional integral for orbitals
        do imo=istart,iend
          a0(imo) = a0(imo) + CO_t(imo,j) * GTFval
        end do
      else if (present(a1)) then
        ! Calculate fractional derivatives for current GTF
        ! x-direction derivative
        GTFdx = differentiate_1st(sftx1, ix, ep, t) * sfty1**iy * sftz1**iz
        ! y-direction derivative
        GTFdy = differentiate_1st(sfty1, iy, ep, t) * sftx1**ix * sftz1**iz
        ! z-direction derivative
        GTFdz = differentiate_1st(sftz1, iz, ep, t) * sftx1**ix * sfty1**iy

        !Calculate fractional derivatives for orbitals
        do imo=istart,iend
          a1(1,imo) = a1(1,imo) + CO_t(imo,j) * GTFdx
          a1(2,imo) = a1(2,imo) + CO_t(imo,j) * GTFdy
          a1(3,imo) = a1(3,imo) + CO_t(imo,j) * GTFdz
        end do
      end if

    end do
  end subroutine orbfracderiv_eval
  real(8) function differentiate_1st(u, iu, ep, t)
    real(8), intent(in) :: u, ep, t(0:1)
    integer, intent(in) :: iu
    ! internal variables
    real(8) :: du
    real(8) :: dum1, dup1
    ! initially, we define derivatives over p, then, apply 2F2 function
    ! d^p /d u^p GTO(u, [v, w])
    ! first part of derivative
    ! by definition, we need first evaluate d^p GTO / d u^p, then integrate
    ! integrate u^(iu-1) part
    dum1 = 0._8
    if (iu >= 1) dum1 = iu * u**(iu - 1) * t(0)
    ! integrate u^(iu+1) part
    dup1 = - 2 * ep * u**(iu + 1) * t(1)
    du = dum1 + dup1
    differentiate_1st = du
  end function differentiate_1st
end subroutine orbfracderiv




!!------------- Calculate \xi(r). Contributed by i.s.ger
! the definition is close to local kinetic energy
! if p = 0, compare with fdens routine (without EDF)
! if p = 1, compare with Lagkin routine when idir = 0
real(8) function fracderiv(x, y, z)
  use mod_2F2, only : p
  implicit none
  real(8), intent(in) :: x, y, z
  real(8) :: a0(nmo)
  real(8) :: a1(3, nmo)
  integer :: imo
  fracderiv = 0._8
  select case (p)
    case(0)
      call orbfracderiv(1, nmo, x, y, z, a0=a0)
      do imo = 1, nmo
        fracderiv = fracderiv + MOocc(imo) * a0(imo)**2
      end do
    case(1)
      call orbfracderiv(1, nmo, x, y, z, a1=a1)
      do imo = 1, nmo
        fracderiv = fracderiv + MOocc(imo) * sum(a1(:, imo)**2)
       end do
    case default
      error stop "p can be only 0 or 1"
  end select
  fracderiv = fracderiv / 2._8
end function




!!----------- Calculate values of all GTFs at a point. Adapted from subroutine "orbderv"
subroutine calcGTFval(x,y,z,GTFvalarr)
implicit real*8 (a-h,o-z)
real*8 x,y,z,GTFvalarr(nprims)

lastcen=-1 !Arbitrary value
GTFvalarr=0D0

do j=1,nprims
	jcen=b(j)%center
	if (jcen/=lastcen) then
		sftx=x-a(jcen)%x
		sfty=y-a(jcen)%y
		sftz=z-a(jcen)%z
		sftx2=sftx*sftx
		sfty2=sfty*sfty
		sftz2=sftz*sftz
		rr=sftx2+sfty2+sftz2
	end if
        
	ep=b(j)%exp
	tmpval=-ep*rr
	lastcen=jcen
	if (tmpval>expcutoff.or.expcutoff>0) then
		expterm=exp(tmpval)
	else
		cycle
	end if
	
	!Calculate value for current GTF
	jtype=b(j)%type
	ix=type2ix(jtype)
	iy=type2iy(jtype)
	iz=type2iz(jtype)
	if (jtype==1) then !Some functype use manually optimized formula for cutting down computational time
	GTFval=expterm
	else if (jtype==2) then
	GTFval=sftx*expterm
	else if (jtype==3) then
	GTFval=sfty*expterm
	else if (jtype==4) then
	GTFval=sftz*expterm
	else if (jtype==5) then
	GTFval=sftx2*expterm
	else if (jtype==6) then
	GTFval=sfty2*expterm
	else if (jtype==7) then
	GTFval=sftz2*expterm
	else if (jtype==8) then
	GTFval=sftx*sfty*expterm
	else if (jtype==9) then
	GTFval=sftx*sftz*expterm
	else if (jtype==10) then
	GTFval=sfty*sftz*expterm
	else !If above conditions are not satisfied (angular moment higher than f), the function will be calculated explicitly
	GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
    end if
	GTFvalarr(j)=GTFval
end do
end subroutine




!!----------- Calculate value of all basis functions at a point
!The global COtr matrix must have been allocated and filled by using COtr=transpose(CO)
subroutine calcbasval(x,y,z,basval)
implicit real*8 (a-h,o-z)
real*8 x,y,z,basval(nbasis),GTFvalarr(nprims)

call calcGTFval(x,y,z,GTFvalarr)

if (isphergau==1) then !For each basis function, only loops GTFs in the same shell for reducing cost
    ibas=0
    iGTF=0
    do ish=1,nshell
        nshbas=shtype2nbas(shtype(ish))
        nshbasCar=shtype2nbas(shtypeCar(ish))
        nshGTF=nshbasCar*shcon(ish)
        is=iGTF+1
        ie=iGTF+nshGTF
        do jbas=ibas+1,ibas+nshbas
            basval(jbas)=sum( GTFvalarr(is:ie)*COtr(is:ie,jbas) )
        end do
        ibas=ibas+nshbas
        iGTF=iGTF+nshGTF
    end do
else !All basis functions are Cartesian, below code is faster than the above general code
    do ibas=1,nbasisCar
        is=primstart(ibas)
        ie=primend(ibas)
        basval(ibas)=sum( GTFvalarr(is:ie)*COtr(is:ie,ibas) )
    end do
end if
end subroutine




!!----------- Calculate contribution from EDFs (recorded in wfx file) to density and corresponding derivatives (up to third-order)
!Only S-type GTFs are supported
! In wfx files, GTFs are used to expand core density
! runtype=1: Only calculate rho, =2: rho+dx/dy/dz =3: rho+dx/dy/dz+dxx/dyy/dzz
!        =4: rho+dx/dy/dz+full Hessian =5: rho+dx/dy/dz+full Hessian+tens3
subroutine EDFrho(runtype,x,y,z,value,grad,hess,tens3)
integer runtype
real*8 x,y,z,value
real*8,optional :: grad(3),hess(3,3),tens3(3,3,3)
value=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0

if (ifPBC>0) then !Consider PBC
    call EDFrho_PBC(runtype,x,y,z,value,grad,hess,tens3)
    return
end if

do i=1,nEDFprims
	sftx=x-a(b_EDF(i)%center)%x
	sfty=y-a(b_EDF(i)%center)%y
	sftz=z-a(b_EDF(i)%center)%z
	sftx2=sftx*sftx
	sfty2=sfty*sfty
	sftz2=sftz*sftz
	rr=sftx2+sfty2+sftz2
	ep=b_EDF(i)%exp
	expterm=exp(-ep*rr)
	value=value+CO_EDF(i)*expterm
	if (runtype>=2) then
		tmp=2*CO_EDF(i)*expterm*ep
		grad(1)=grad(1)-tmp*sftx
		grad(2)=grad(2)-tmp*sfty
		grad(3)=grad(3)-tmp*sftz
		if (runtype>=3) then
			hess(1,1)=hess(1,1)+tmp*(2*ep*sftx2-1)
			hess(2,2)=hess(2,2)+tmp*(2*ep*sfty2-1)
			hess(3,3)=hess(3,3)+tmp*(2*ep*sftz2-1)
			if (runtype>=4) then
				epep4=ep*ep*4
				tmp2=CO_EDF(i)*epep4*expterm
				hess(1,2)=hess(1,2)+tmp2*sftx*sfty
				hess(1,3)=hess(1,3)+tmp2*sftx*sftz
				hess(2,3)=hess(2,3)+tmp2*sfty*sftz
				hess(2,1)=hess(1,2)
				hess(3,1)=hess(1,3)
				hess(3,2)=hess(2,3)
				if (runtype>=5) then
					tmp3=CO_EDF(i)*epep4*expterm
					tens3(1,1,1)=tens3(1,1,1)+tmp3*sftx*(3-2*ep*sftx2)
					tens3(2,2,2)=tens3(2,2,2)+tmp3*sfty*(3-2*ep*sfty2)
					tens3(3,3,3)=tens3(3,3,3)+tmp3*sftz*(3-2*ep*sftz2)
					tens3(1,2,2)=tens3(1,2,2)+tmp3*sftx*(1-2*ep*sfty2)
					tens3(1,1,2)=tens3(1,1,2)+tmp3*sfty*(1-2*ep*sftx2)
					tens3(1,1,3)=tens3(1,1,3)+tmp3*sftz*(1-2*ep*sftx2)
					tens3(1,3,3)=tens3(1,3,3)+tmp3*sftx*(1-2*ep*sftz2)
					tens3(2,3,3)=tens3(2,3,3)+tmp3*sfty*(1-2*ep*sftz2)
					tens3(2,2,3)=tens3(2,2,3)+tmp3*sftz*(1-2*ep*sfty2)
					tens3(1,2,3)=tens3(1,2,3)-CO_EDF(i)*8*ep**3*sftx*sfty*sftz*expterm
					tens3(1,2,1)=tens3(1,1,2) !dxyx=dxxy
					tens3(1,3,1)=tens3(1,1,3) !dxzx=dxxz
					tens3(1,3,2)=tens3(1,2,3) !dxzy=dxyz
					tens3(2,1,1)=tens3(1,1,2) !dyxx=dxxy
					tens3(2,1,2)=tens3(1,2,2) !dyxy=dxyy
					tens3(2,1,3)=tens3(1,2,3) !dyxz=dxyz
					tens3(2,2,1)=tens3(1,2,2) !dyyx=dxyy
					tens3(2,3,1)=tens3(1,2,3) !dyzx=dxyz
					tens3(2,3,2)=tens3(2,2,3) !dyzy=dyyz
					tens3(3,1,1)=tens3(1,1,3) !dzxx=dxxz
					tens3(3,1,2)=tens3(1,2,3) !dzxy=dxyz
					tens3(3,1,3)=tens3(1,3,3) !dzxz=dxzz
					tens3(3,2,1)=tens3(1,2,3) !dzyx=dxyz
					tens3(3,2,2)=tens3(2,2,3) !dzyy=dyyz
					tens3(3,2,3)=tens3(2,3,3) !dzyz=dyzz
					tens3(3,3,1)=tens3(1,3,3) !dzzx=dxzz
					tens3(3,3,2)=tens3(2,3,3) !dzzy=dyzz
				end if
			end if
		end if
	end if
end do
end subroutine




!!----------- The same as EDFrho, but for PBC case
subroutine EDFrho_PBC(runtype,x,y,z,value,grad,hess,tens3)
integer runtype
real*8 x,y,z,value,tvec(3)
real*8,optional :: grad(3),hess(3,3),tens3(3,3,3)
value=0D0
if (present(grad)) grad=0D0
if (present(hess)) hess=0D0
if (present(tens3)) tens3=0D0

call getpointcell(x,y,z,ic,jc,kc)
do icell=ic-PBCnx,ic+PBCnx
    do jcell=jc-PBCny,jc+PBCny
        do kcell=kc-PBCnz,kc+PBCnz
            call tvec_PBC(icell,jcell,kcell,tvec)
            xmove=tvec(1)
            ymove=tvec(2)
            zmove=tvec(3)
            do i=1,nEDFprims
	            sftx=x-(a(b_EDF(i)%center)%x+xmove)
	            sfty=y-(a(b_EDF(i)%center)%y+ymove)
	            sftz=z-(a(b_EDF(i)%center)%z+zmove)
	            sftx2=sftx*sftx
	            sfty2=sfty*sfty
	            sftz2=sftz*sftz
	            rr=sftx2+sfty2+sftz2
	            ep=b_EDF(i)%exp
                if (-ep*rr>expcutoff_PBC.or.expcutoff_PBC>0) then
					expterm=exp(-ep*rr)
	            else
                    cycle
	            end if
                
	            value=value+CO_EDF(i)*expterm
                
	            if (runtype>=2) then
		            tmp=2*CO_EDF(i)*expterm*ep
		            grad(1)=grad(1)-tmp*sftx
		            grad(2)=grad(2)-tmp*sfty
		            grad(3)=grad(3)-tmp*sftz
		            if (runtype>=3) then
			            hess(1,1)=hess(1,1)+tmp*(2*ep*sftx2-1)
			            hess(2,2)=hess(2,2)+tmp*(2*ep*sfty2-1)
			            hess(3,3)=hess(3,3)+tmp*(2*ep*sftz2-1)
			            if (runtype>=4) then
				            epep4=ep*ep*4
				            tmp2=CO_EDF(i)*epep4*expterm
				            hess(1,2)=hess(1,2)+tmp2*sftx*sfty
				            hess(1,3)=hess(1,3)+tmp2*sftx*sftz
				            hess(2,3)=hess(2,3)+tmp2*sfty*sftz
				            hess(2,1)=hess(1,2)
				            hess(3,1)=hess(1,3)
				            hess(3,2)=hess(2,3)
				            if (runtype>=5) then
					            tmp3=CO_EDF(i)*epep4*expterm
					            tens3(1,1,1)=tens3(1,1,1)+tmp3*sftx*(3-2*ep*sftx2)
					            tens3(2,2,2)=tens3(2,2,2)+tmp3*sfty*(3-2*ep*sfty2)
					            tens3(3,3,3)=tens3(3,3,3)+tmp3*sftz*(3-2*ep*sftz2)
					            tens3(1,2,2)=tens3(1,2,2)+tmp3*sftx*(1-2*ep*sfty2)
					            tens3(1,1,2)=tens3(1,1,2)+tmp3*sfty*(1-2*ep*sftx2)
					            tens3(1,1,3)=tens3(1,1,3)+tmp3*sftz*(1-2*ep*sftx2)
					            tens3(1,3,3)=tens3(1,3,3)+tmp3*sftx*(1-2*ep*sftz2)
					            tens3(2,3,3)=tens3(2,3,3)+tmp3*sfty*(1-2*ep*sftz2)
					            tens3(2,2,3)=tens3(2,2,3)+tmp3*sftz*(1-2*ep*sfty2)
					            tens3(1,2,3)=tens3(1,2,3)-CO_EDF(i)*8*ep**3*sftx*sfty*sftz*expterm
					            tens3(1,2,1)=tens3(1,1,2) !dxyx=dxxy
					            tens3(1,3,1)=tens3(1,1,3) !dxzx=dxxz
					            tens3(1,3,2)=tens3(1,2,3) !dxzy=dxyz
					            tens3(2,1,1)=tens3(1,1,2) !dyxx=dxxy
					            tens3(2,1,2)=tens3(1,2,2) !dyxy=dxyy
					            tens3(2,1,3)=tens3(1,2,3) !dyxz=dxyz
					            tens3(2,2,1)=tens3(1,2,2) !dyyx=dxyy
					            tens3(2,3,1)=tens3(1,2,3) !dyzx=dxyz
					            tens3(2,3,2)=tens3(2,2,3) !dyzy=dyyz
					            tens3(3,1,1)=tens3(1,1,3) !dzxx=dxxz
					            tens3(3,1,2)=tens3(1,2,3) !dzxy=dxyz
					            tens3(3,1,3)=tens3(1,3,3) !dzxz=dxzz
					            tens3(3,2,1)=tens3(1,2,3) !dzyx=dxyz
					            tens3(3,2,2)=tens3(2,2,3) !dzyy=dyyz
					            tens3(3,2,3)=tens3(2,3,3) !dzyz=dyzz
					            tens3(3,3,1)=tens3(1,3,3) !dzzx=dxzz
					            tens3(3,3,2)=tens3(2,3,3) !dzzy=dyzz
				            end if
			            end if
		            end if
	            end if
            end do
        end do
    end do
end do
end subroutine




!!-------- A general routine used to calculate value, gradient and Hessian matrix for all real space functions
! itype=1 Only calculate value and grad
! itype=2 Calculate value, gradient and Hessian
!If idiagonly_in is presented and is 1, then when full numerical Hessian is calculated, only diagonal part will be evaluated to save cost
subroutine gencalchessmat(itype,ifunc,x,y,z,value,grad,hess,idiagonly_in)
integer ifunc,itype
real*8 x,y,z,value,grad(3),hess(3,3),tens3(3,3,3)
real*8 gradaddx(3),gradminx(3),gradaddy(3),gradminy(3),gradaddz(3),gradminz(3)
character selELFLOL*3
integer,optional :: idiagonly_in

diff=8D-4
denom=2D0*diff
if (ifunc==9) selELFLOL="ELF"
if (ifunc==10) selELFLOL="LOL"

idiagonly=0
if (present(idiagonly_in)) then
	if (idiagonly_in==1) idiagonly=1
end if

!For functions whose both analytic gradient and Hessian are available, evaluate them and then return
!If comment one of them, then the gradient and Hessian of corresponding function will be calculated numerically
if (ifunc==1) then
	call calchessmat_dens(itype,x,y,z,value,grad,hess)
	return
else if (ifunc==2) then
	call calchessmat_rhograd(itype,x,y,z,value,grad,hess)
	return
else if (ifunc==4) then
	call calchessmat_orb(itype,iorbsel,x,y,z,value,grad,hess)
	return
else if (ifunc==11) then !Local information entropy
    call calchessmat_Shannon(1,x,y,z,value,grad,hess)
	return
else if (ifunc==13) then !RDG
    call calchessmat_IRI_RDG(itype,2,x,y,z,value,grad,hess)
	return
else if (ifunc==24) then !IRI
    call calchessmat_IRI_RDG(itype,1,x,y,z,value,grad,hess)
	return
else if (ifunc==25) then !van der Waals potential
    call calchessmat_vdWpot(itype,x,y,z,value,grad,hess)
	return
else if (ifunc==100) then
	if (iuserfunc==49) then !Relative Shannon entropy density
		call calchessmat_relShannon(x,y,z,value,grad,hess)
		return
	else if (iuserfunc==50) then !Shannon entropy density
		call calchessmat_Shannon(2,x,y,z,value,grad,hess)
		return
	else if (iuserfunc==51) then !Fisher information density
		call calchessmat_Fisherinfo(x,y,z,value,grad,hess)
		return
	else if (iuserfunc==52) then !Second Fisher information density
		call calchessmat_second_Fisherinfo(x,y,z,value,grad,hess)
		return
	else if (iuserfunc==99) then !IRI
		call calchessmat_IRI_RDG(itype,1,x,y,z,value,grad,hess)
		return
    end if
end if
	
!Calculate gradient and meantime value
!For other functions aside from above ones, analytic Hessian or even gradient has not been realized
if (ifunc==3) then
	call calchessmat_lapl(1,x,y,z,value,grad,hess)
else if (ifunc==9.or.ifunc==10) then
	if (ELFLOL_type==0) then !Analytic gradient for original ELF/LOL
		call calchessmat_ELF_LOL(1,x,y,z,value,grad,hess,selELFLOL)
	else !Numerical gradient for other definitions of ELF/LOL
		value=ELF_LOL(x,y,z,selELFLOL)
		xadd=ELF_LOL(x+diff,y,z,selELFLOL)
		xmin=ELF_LOL(x-diff,y,z,selELFLOL)
		yadd=ELF_LOL(x,y+diff,z,selELFLOL)
		ymin=ELF_LOL(x,y-diff,z,selELFLOL)
		zadd=ELF_LOL(x,y,z+diff,selELFLOL)
		zmin=ELF_LOL(x,y,z-diff,selELFLOL)
		grad(1)=(xadd-xmin)/denom
		grad(2)=(yadd-ymin)/denom
		grad(3)=(zadd-zmin)/denom
	end if
else !User general interface to calculate value for all other functions
	value=calcfuncall(ifunc,x,y,z)
	xadd=calcfuncall(ifunc,x+diff,y,z)
	xmin=calcfuncall(ifunc,x-diff,y,z)
	yadd=calcfuncall(ifunc,x,y+diff,z)
	ymin=calcfuncall(ifunc,x,y-diff,z)
	zadd=calcfuncall(ifunc,x,y,z+diff)
	zmin=calcfuncall(ifunc,x,y,z-diff)
	grad(1)=(xadd-xmin)/denom
	grad(2)=(yadd-ymin)/denom
	grad(3)=(zadd-zmin)/denom
end if

!Calculate Hessian semi (namely based on analyic gradient) or pure (based on function value) numerically
if (itype==2) then
	if (ifunc==3) then !Use semi-analytic for Hessian of Laplacian
		call calchessmat_lapl(1,x+diff,y,z,tmpval,gradaddx,hess)
		call calchessmat_lapl(1,x-diff,y,z,tmpval,gradminx,hess)
		call calchessmat_lapl(1,x,y+diff,z,tmpval,gradaddy,hess)
		call calchessmat_lapl(1,x,y-diff,z,tmpval,gradminy,hess)
		call calchessmat_lapl(1,x,y,z+diff,tmpval,gradaddz,hess)
		call calchessmat_lapl(1,x,y,z-diff,tmpval,gradminz,hess)
		hess(1,1)=(gradaddx(1)-gradminx(1))/denom
		hess(2,2)=(gradaddy(2)-gradminy(2))/denom
		hess(3,3)=(gradaddz(3)-gradminz(3))/denom
		hess(1,2)=(gradaddy(1)-gradminy(1))/denom
		hess(2,3)=(gradaddz(2)-gradminz(2))/denom
		hess(1,3)=(gradaddz(1)-gradminz(1))/denom
		hess(2,1)=hess(1,2)
		hess(3,2)=hess(2,3)
		hess(3,1)=hess(1,3)
		return !Do not do below procedures to generate pure numerical Hessian
	else if (ifunc==9.or.ifunc==10) then !ELF, LOL
		if (ELFLOL_type==0) then !Use semi-analytic for Hessian of original ELF/LOL
			call calchessmat_ELF_LOL(1,x+diff,y,z,tmpval,gradaddx,hess,selELFLOL)
			call calchessmat_ELF_LOL(1,x-diff,y,z,tmpval,gradminx,hess,selELFLOL)
			call calchessmat_ELF_LOL(1,x,y+diff,z,tmpval,gradaddy,hess,selELFLOL)
			call calchessmat_ELF_LOL(1,x,y-diff,z,tmpval,gradminy,hess,selELFLOL)
			call calchessmat_ELF_LOL(1,x,y,z+diff,tmpval,gradaddz,hess,selELFLOL)
			call calchessmat_ELF_LOL(1,x,y,z-diff,tmpval,gradminz,hess,selELFLOL)
			hess(1,1)=(gradaddx(1)-gradminx(1))/denom
			hess(2,2)=(gradaddy(2)-gradminy(2))/denom
			hess(3,3)=(gradaddz(3)-gradminz(3))/denom
			hess(1,2)=(gradaddy(1)-gradminy(1))/denom
			hess(2,3)=(gradaddz(2)-gradminz(2))/denom
			hess(1,3)=(gradaddz(1)-gradminz(1))/denom
			hess(2,1)=hess(1,2)
			hess(3,2)=hess(2,3)
			hess(3,1)=hess(1,3)
			return
		else !For other definitions of ELF/LOL, use pure numerical Hessian
			xaddxadd=ELF_LOL(x+2D0*diff,y,z,selELFLOL)
			xminxmin=ELF_LOL(x-2D0*diff,y,z,selELFLOL)
			yaddyadd=ELF_LOL(x,y+2D0*diff,z,selELFLOL)
			yminymin=ELF_LOL(x,y-2D0*diff,z,selELFLOL)
			zaddzadd=ELF_LOL(x,y,z+2D0*diff,selELFLOL)
			zminzmin=ELF_LOL(x,y,z-2D0*diff,selELFLOL)
			xaddyadd=ELF_LOL(x+diff,y+diff,z,selELFLOL)
			xminyadd=ELF_LOL(x-diff,y+diff,z,selELFLOL)
			xaddymin=ELF_LOL(x+diff,y-diff,z,selELFLOL)
			xminymin=ELF_LOL(x-diff,y-diff,z,selELFLOL)
			yaddzadd=ELF_LOL(x,y+diff,z+diff,selELFLOL)
			yminzadd=ELF_LOL(x,y-diff,z+diff,selELFLOL)
			yaddzmin=ELF_LOL(x,y+diff,z-diff,selELFLOL)
			yminzmin=ELF_LOL(x,y-diff,z-diff,selELFLOL)
			xaddzadd=ELF_LOL(x+diff,y,z+diff,selELFLOL)
			xminzadd=ELF_LOL(x-diff,y,z+diff,selELFLOL)
			xaddzmin=ELF_LOL(x+diff,y,z-diff,selELFLOL)
			xminzmin=ELF_LOL(x-diff,y,z-diff,selELFLOL)
		end if
	else !User general interface to calculate for all other functions
		xaddxadd=calcfuncall(ifunc,x+2D0*diff,y,z)
		xminxmin=calcfuncall(ifunc,x-2D0*diff,y,z)
		yaddyadd=calcfuncall(ifunc,x,y+2D0*diff,z)
		yminymin=calcfuncall(ifunc,x,y-2D0*diff,z)
		zaddzadd=calcfuncall(ifunc,x,y,z+2D0*diff)
		zminzmin=calcfuncall(ifunc,x,y,z-2D0*diff)
        if (idiagonly==0) then
			xaddyadd=calcfuncall(ifunc,x+diff,y+diff,z)
			xminyadd=calcfuncall(ifunc,x-diff,y+diff,z)
			xaddymin=calcfuncall(ifunc,x+diff,y-diff,z)
			xminymin=calcfuncall(ifunc,x-diff,y-diff,z)
			yaddzadd=calcfuncall(ifunc,x,y+diff,z+diff)
			yminzadd=calcfuncall(ifunc,x,y-diff,z+diff)
			yaddzmin=calcfuncall(ifunc,x,y+diff,z-diff)
			yminzmin=calcfuncall(ifunc,x,y-diff,z-diff)
			xaddzadd=calcfuncall(ifunc,x+diff,y,z+diff)
			xminzadd=calcfuncall(ifunc,x-diff,y,z+diff)
			xaddzmin=calcfuncall(ifunc,x+diff,y,z-diff)
			xminzmin=calcfuncall(ifunc,x-diff,y,z-diff)
        end if
	end if 
	!Collect above temporary data to evaluate pure numerical Hessian
	gradx_yadd=(xaddyadd-xminyadd)/denom
	gradx_ymin=(xaddymin-xminymin)/denom
	grady_zadd=(yaddzadd-yminzadd)/denom
	grady_zmin=(yaddzmin-yminzmin)/denom
	gradx_zadd=(xaddzadd-xminzadd)/denom
	gradx_zmin=(xaddzmin-xminzmin)/denom
	hess(1,1)=(xaddxadd-2*value+xminxmin)/denom**2
	hess(2,2)=(yaddyadd-2*value+yminymin)/denom**2
	hess(3,3)=(zaddzadd-2*value+zminzmin)/denom**2
	hess(1,2)=(gradx_yadd-gradx_ymin)/denom
	hess(2,3)=(grady_zadd-grady_zmin)/denom
	hess(1,3)=(gradx_zadd-gradx_zmin)/denom
	hess(2,1)=hess(1,2)
	hess(3,2)=hess(2,3)
	hess(3,1)=hess(1,3)
end if
end subroutine




!!----------- Calculate orbital wavefunction value (fmo=function for outputting orbital (not necessarily MO))
real*8 function fmo(x,y,z,id)
real*8 x,y,z,orbval(nmo)
integer id
call orbderv(1,id,id,x,y,z,orbval)
fmo=orbval(id)
end function




!!----------- Calculate probability density of an orbital wavefunction value
real*8 function forbdens(x,y,z,id)
real*8 x,y,z,orbval(nmo),orbval_imag(nmo)
integer id
if (ifPBC==0) then
	call orbderv(1,id,id,x,y,z,orbval)
	forbdens=orbval(id)**2
else
	call orbderv_PBC(1,id,id,x,y,z,orbval,itype_in=1)
	call orbderv_PBC(1,id,id,x,y,z,orbval_imag,itype_in=2)
	forbdens=orbval(id)**2+orbval_imag(id)**2
end if
end function




!!----- Calculate orbital wavefunction Hessian matrix, store to hess, output value and gradient vector at same time
!itype=1 Only calculate value and gradient, not Hessian
!itype=2 Calculate value, gradient and Hessian
subroutine calchessmat_orb(itype,id,x,y,z,value,grad,hess)
integer id,itype
real*8 x,y,z,grad(3),hess(3,3),value,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
if (itype==1) then
	call orbderv(2,id,id,x,y,z,wfnval,wfnderv,wfnhess)
else if (itype==2) then
	call orbderv(4,id,id,x,y,z,wfnval,wfnderv,wfnhess)
end if
value=wfnval(id)
grad=wfnderv(:,id)
hess=wfnhess(:,:,id)
end subroutine



!-------- Calculate electron density
real*8 function fdens(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
fdens=0D0
do i=1,nmo
	fdens=fdens+MOocc(i)*wfnval(i)**2
end do
!Add contribution of Electron density function
if (nEDFprims/=0) then
	call EDFrho(1,x,y,z,EDFdens)
	fdens=fdens+EDFdens
end if
end function



!!---- Calculate density using Lagrangian interpolation at (x,y,z) for iatm based on "atmraddens", which is loaded atomic radial density or that generated by e.g. Hirshfeld-I or MBIS procedure
real*8 function fdens_rad(iatm,x,y,z)
use defvar
use util
integer iatm,npt
real*8 x,y,z,r,rnouse
npt=atmradnpt(iatm)
r=dsqrt((a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2)
call lagintpol(atmradpos(1:npt),atmraddens(1:npt,iatm),npt,r,fdens_rad,rnouse,rnouse,1)
end function




!!----------------- Calculate spin or Alpha or Beta electron density
!itype='s' output spin density, ='a' output alpha density, ='b' output beta density
real*8 function fspindens(x,y,z,itype)
real*8 x,y,z,wfnval(nmo)
character itype
call orbderv(1,1,nmo,x,y,z,wfnval)
adens=0D0
bdens=0D0
do i=1,nmo
	if (MOtype(i)==1) then
		adens=adens+MOocc(i)*wfnval(i)**2
	else if (MOtype(i)==2) then
		bdens=bdens+MOocc(i)*wfnval(i)**2
	else if (MOtype(i)==0) then
		adens=adens+MOocc(i)/2D0*wfnval(i)**2
		bdens=bdens+MOocc(i)/2D0*wfnval(i)**2
	end if
end do
if (itype=='s') then
	fspindens=adens-bdens
	if (ipolarpara==1) fspindens=fspindens/(adens+bdens)
else if (itype=='a') then
	fspindens=adens
else if (itype=='b') then
	fspindens=bdens
end if
end function




!!--------------- Calculate gradient norm of rho and RDG (reduced density gradient)
!label=x/y/z output 1-order derivation of x/y/z, =t get norm, =r get RDG
real*8 function fgrad(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),gradrho(3),EDFgrad(3),sumgrad2
character label
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
rho=0D0
gradrho=0D0
do i=1,nmo
	rho=rho+MOocc(i)*wfnval(i)**2
	gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do
gradrho=2*gradrho
!Add in contribution of electron density function (EDF)
if (nEDFprims/=0) then
	call EDFrho(2,x,y,z,EDFdens,EDFgrad)
	rho=rho+EDFdens
	gradrho=gradrho+EDFgrad
end if
if (label=='x') then
	fgrad=gradrho(1)
else if (label=='y') then
	fgrad=gradrho(2)
else if (label=='z') then
	fgrad=gradrho(3)
else if (label=='t') then
	fgrad=dsqrt( sum(gradrho(:)**2) )
else if (label=='r') then
	sumgrad2=sum(gradrho(:)**2)
	if (RDG_maxrho/=0D0.and.rho>=RDG_maxrho) then
		fgrad=100D0
	else if (sumgrad2==0D0.or.rho==0D0) then !This occurs at distant region when exponent cutoff is used, the actual value should be very large. In order to avoid denominator become zero, we set it artifically to a big value
		fgrad=999D0
	else
		fgrad=0.161620459673995D0*dsqrt(sumgrad2)/rho**(4D0/3D0) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
	end if
end if
end function


!!--- Simultaneously generate electron density, gradient norm for alpha and beta electrons, as well as dot product between grada and gradb
!Mainly used to evaluate DFT functional. EDF is not taken into account
!adens/bdens/tdens means the density of alpha/beta/total density, similar for *grad
subroutine gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
real*8 x,y,z,adens,bdens,agrad,bgrad,abgrad,wfnval(nmo),wfnderv(3,nmo),gradrhoa(3),gradrhob(3),gradrhot(3),tmparr(3)
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
adens=0D0
bdens=0D0
gradrhoa=0D0
gradrhob=0D0
do i=1,nmo
	if (MOtype(i)==1) then
		adens=adens+MOocc(i)*wfnval(i)**2
		gradrhoa(:)=gradrhoa(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
	else if (MOtype(i)==2) then
		bdens=bdens+MOocc(i)*wfnval(i)**2
		gradrhob(:)=gradrhob(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
	else if (MOtype(i)==0) then
		tmpval=MOocc(i)/2D0*wfnval(i)**2
		adens=adens+tmpval
		bdens=bdens+tmpval
		tmparr(:)=MOocc(i)/2D0*wfnval(i)*wfnderv(:,i)
		gradrhoa(:)=gradrhoa(:)+tmparr(:)
		gradrhob(:)=gradrhob(:)+tmparr(:)
	end if
end do
tdens=adens+bdens
gradrhoa=gradrhoa*2
gradrhob=gradrhob*2
gradrhot=gradrhoa+gradrhob
agrad=dsqrt(sum(gradrhoa**2))
bgrad=dsqrt(sum(gradrhob**2))
tgrad=dsqrt(sum(gradrhot**2))
abgrad=sum(gradrhoa*gradrhob)
end subroutine


!!--- Simultaneously generate electron density, gradient vector and Laplacian for alpha and beta electrons
!Mainly used to evaluate XC potential of open-shell case. EDF is not taken into account
subroutine gendens_gradvec_lapl_ab(x,y,z,rhoa,rhob,gradrhoa,gradrhob,sigaa,sigbb,sigab,laplrhoa,laplrhob)
real*8 x,y,z,rhoa,rhob,gradrhoa(3),gradrhob(3),sigaa,sigbb,sigab,laplrhoa,laplrhob,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
real*8 gradtmp(3),lapltmp
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
rhoa=0D0
rhob=0D0
gradrhoa=0D0
gradrhob=0D0
laplrhoa=0D0
laplrhob=0D0
do i=1,nmo
	rhotmp=MOocc(i)*wfnval(i)**2
    gradtmp(:)=MOocc(i)*wfnval(i)*wfnderv(:,i)
	xtmp=wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i)
	ytmp=wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i)
	ztmp=wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i)
    lapltmp=MOocc(i)*(xtmp+ytmp+ztmp)
	if (MOtype(i)==1) then
		rhoa=rhoa+rhotmp
		gradrhoa(:)=gradrhoa(:)+gradtmp(:)
        laplrhoa=laplrhoa+lapltmp
	else if (MOtype(i)==2) then
		rhob=rhob+rhotmp
		gradrhob(:)=gradrhob(:)+gradtmp(:)
        laplrhob=laplrhob+lapltmp
	else if (MOtype(i)==0) then
		rhoa=rhoa+rhotmp/2
		rhob=rhob+rhotmp/2
		gradrhoa(:)=gradrhoa(:)+gradtmp(:)/2
		gradrhob(:)=gradrhob(:)+gradtmp(:)/2
        laplrhoa=laplrhoa+lapltmp/2
        laplrhob=laplrhob+lapltmp/2
	end if
end do
gradrhoa=gradrhoa*2
gradrhob=gradrhob*2
laplrhoa=laplrhoa*2
laplrhob=laplrhob*2
sigaa=sum(gradrhoa(:)**2)
sigbb=sum(gradrhob(:)**2)
sigab=sum(gradrhoa(:)*gradrhob(:))
end subroutine




!!----------------- Calculate Laplacian of electron density
!label=x/y/z output 2-order derivative of electron density respect to xx/yy/zz; &
!=t get their summing; =s get der2rho/rho^(5/3), which is used LSB's project
real*8 function flapl(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),laplx,laply,laplz,EDFgrad(3),EDFhess(3,3)
character label
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
laplx=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
laply=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
laplz=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
!Add in contribution of electron density function, assume EDFs are S type
if (nEDFprims/=0) then
	call EDFrho(3,x,y,z,EDFdens,EDFgrad,EDFhess)
	laplx=laplx+EDFhess(1,1)
	laply=laply+EDFhess(2,2)
	laplz=laplz+EDFhess(3,3)
end if
if (label=='t') then
	flapl=laplx+laply+laplz
	flapl=flapl*laplfac !laplfac is an external variable
else if (label=='x') then
	flapl=laplx
else if (label=='y') then
	flapl=laply
else if (label=='z') then
	flapl=laplz
else if (label=='s') then
	dens=sum(MOocc(1:nmo)*wfnval(1:nmo)**2)
	if (nEDFprims/=0) dens=dens+EDFdens
	flapl=(laplx+laply+laplz)/dens**(5D0/3D0)
end if
end function




!!------ Calculate electron density, its gradient and Hessian matrix
!itype=1 Only calculate value and grad, not Hessian
!itype=2 Calculate value, gradient and Hessian
subroutine calchessmat_dens(itype,x,y,z,elerho,elegrad,elehess)
real*8 x,y,z,elerho,elegrad(3),elehess(3,3),wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),EDFgrad(3),EDFhess(3,3)
integer itype
call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
elerho=sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnval(1:nmo) )
do itmp=1,3
	elegrad(itmp)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(itmp,1:nmo) )
end do
if (itype==2) then
	elehess(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
	elehess(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
	elehess(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
	elehess(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
	elehess(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
	elehess(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
	elehess(2,1)=elehess(1,2)
	elehess(3,2)=elehess(2,3)
	elehess(3,1)=elehess(1,3)
end if

!Add in contribution of electron density function, assume EDFs are S type
if (nEDFprims/=0) then
	call EDFrho(4,x,y,z,EDFdens,EDFgrad,EDFhess)
	elerho=elerho+EDFdens
	elegrad=elegrad+EDFgrad
	elehess=elehess+EDFhess
end if
end subroutine




!!------ Same as calchessmat_dens, but for promolecular wavefunction stored in memory
!itype=1 Only calculate value and grad, not Hessian
!itype=2 Calculate value, gradient and Hessian
subroutine calchessmat_dens_promol(itype,x,y,z,elerho,elegrad,elehess)
real*8 x,y,z,elerho,elegrad(3),elehess(3,3),wfnval(nmo_pmol),wfnderv(3,nmo_pmol),wfnhess(3,3,nmo_pmol)
integer itype
call orbderv_pmol(4,1,nmo_pmol,x,y,z,wfnval,wfnderv,wfnhess)
elerho=sum( MOocc_pmol(1:nmo_pmol)*wfnval(1:nmo_pmol)*wfnval(1:nmo_pmol) )
do itmp=1,3
	elegrad(itmp)=2*sum( MOocc_pmol(1:nmo_pmol)*wfnval(1:nmo_pmol)*wfnderv(itmp,1:nmo_pmol) )
end do
if (itype==2) then
	elehess(1,1)=2*sum( MOocc_pmol(1:nmo_pmol)*( wfnderv(1,1:nmo_pmol)**2 + wfnval(1:nmo_pmol)*wfnhess(1,1,1:nmo_pmol) ) )
	elehess(2,2)=2*sum( MOocc_pmol(1:nmo_pmol)*( wfnderv(2,1:nmo_pmol)**2 + wfnval(1:nmo_pmol)*wfnhess(2,2,1:nmo_pmol) ) )
	elehess(3,3)=2*sum( MOocc_pmol(1:nmo_pmol)*( wfnderv(3,1:nmo_pmol)**2 + wfnval(1:nmo_pmol)*wfnhess(3,3,1:nmo_pmol) ) )
	elehess(1,2)=2*sum( MOocc_pmol(1:nmo_pmol)*( wfnderv(1,1:nmo_pmol)*wfnderv(2,1:nmo_pmol)+wfnhess(1,2,1:nmo_pmol)*wfnval(1:nmo_pmol) ) )
	elehess(2,3)=2*sum( MOocc_pmol(1:nmo_pmol)*( wfnderv(2,1:nmo_pmol)*wfnderv(3,1:nmo_pmol)+wfnhess(2,3,1:nmo_pmol)*wfnval(1:nmo_pmol) ) )
	elehess(1,3)=2*sum( MOocc_pmol(1:nmo_pmol)*( wfnderv(1,1:nmo_pmol)*wfnderv(3,1:nmo_pmol)+wfnhess(1,3,1:nmo_pmol)*wfnval(1:nmo_pmol) ) )
	elehess(2,1)=elehess(1,2)
	elehess(3,2)=elehess(2,3)
	elehess(3,1)=elehess(1,3)
end if
end subroutine




!!----- Calculate electron density, its gradient, Hessian and 3rd derivative tensor. EDF is taken into account
!This is the most general and elegant code, can also fully replace calchessmat_dens
subroutine rho_tensor(x,y,z,value,grad,hess,tens3)
implicit real*8 (a-h,o-z)
real*8 x,y,z,value,grad(3),hess(3,3),tens3(3,3,3),wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo),EDFgrad(3),EDFhess(3,3),EDFtens3(3,3,3)

call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)
!rho
value=sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnval(1:nmo) )
!gradient of rho
do i=1,3
	grad(i)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(i,1:nmo) )
end do
!Hessian of rho
hess=0
do i=1,3
	do j=i,3
		tmpval=0
		do imo=1,nmo
			if (MOocc(i)==0) cycle
			tmpval=tmpval+MOocc(imo)*(wfnderv(i,imo)*wfnderv(j,imo)+wfnval(imo)*wfnhess(i,j,imo))
        end do
        hess(i,j)=2*tmpval
        hess(j,i)=hess(i,j)
	end do
end do

!3rd derivative tensor of rho
tens3=0D0
do i=1,3
	do j=i,3
		do k=j,3
			tmpval=0
			do imo=1,nmo
				if (MOocc(i)==0) cycle
                tmpval=tmpval+MOocc(imo)*( wfnderv(j,imo)*wfnhess(i,k,imo) + wfnderv(i,imo)*wfnhess(j,k,imo) + &
                wfnderv(k,imo)*wfnhess(i,j,imo) + wfnval(imo)*wfntens3(i,j,k,imo) )
			end do
            tens3(i,j,k)=2*tmpval
            tens3(i,k,j)=tens3(i,j,k)
            tens3(k,j,i)=tens3(i,j,k)
            tens3(j,i,k)=tens3(i,j,k)
            tens3(k,i,j)=tens3(i,j,k)
            tens3(j,k,i)=tens3(i,j,k)
        end do
    end do
end do

!Add in contribution of electron density function, assume EDFs are S type
if (nEDFprims/=0) then
	call EDFrho(5,x,y,z,EDFdens,EDFgrad,EDFhess,EDFtens3)
	rho=rho+EDFdens
	grad=grad+EDFgrad
	hess=hess+EDFhess
    tens3=tens3+EDFtens3
end if
end subroutine




!!--------- Calculate gradient and Hessian of gradient norm of electron density
!itype=1 Only calculate value and grad, =2 Calculate value, gradient and Hessian
subroutine calchessmat_rhograd(itype,x,y,z,value,grad,hess)
implicit real*8 (a-h,o-z)
integer itype
real*8 x,y,z,value,grad(3),hess(3,3),rhograd(3),rhohess(3,3),rhotens3(3,3,3)

call rho_tensor(x,y,z,rho,rhograd,rhohess,rhotens3)
rhograd2=sum(rhograd**2)
gnorm=dsqrt(rhograd2)
value=gnorm

grad=0
do idir=1,3
	do jdir=1,3
		grad(idir)=grad(idir)+rhograd(jdir)*rhohess(idir,jdir)
    end do
end do
grad=grad/gnorm

if (itype==2) then
	hess=0
	do i=1,3
		do j=1,3
			tmp=0
			tmp2=0
			do k=1,3
				tmp=tmp+rhograd(k)*rhohess(k,i)
				tmp2=tmp2+rhohess(k,i)*rhohess(k,j)+rhograd(k)*rhotens3(k,i,j)
            end do
			hess(i,j)=-1D0/rhograd2*grad(j)*tmp+tmp2/gnorm
		end do
	end do
end if
end subroutine



!!------------- Calculate Laplacian of electron density, its gradient and Hessian matrix
!itype=1 calculate value, gradient
!itype=2 calculate value, gradient and Hessian (Not available)
subroutine calchessmat_lapl(itype,x,y,z,value,grad,hess)
use util
real*8 x,y,z,value,grad(3),hess(3,3)
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo),rhotens3(3,3,3)
real*8 EDFgrad(3),EDFhess(3,3),EDFtens3(3,3,3)
integer itype
!Numerically verify 3-order derivative of orbital wavefunction
! diff=1D-5
! call orbderv(5,1,nmo,x+diff,y,z,wfnval,wfnderv,wfnhess,wfntens3)
! t1=wfnhess(3,3,1)
! call orbderv(5,1,nmo,x-diff,y,z,wfnval,wfnderv,wfnhess,wfntens3)
! t2=wfnhess(3,3,1)
! write(*,*) (t1-t2)/(2*diff)

call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)
rhotens3=0D0
dxx=0D0
dyy=0D0
dzz=0D0
do i=1,nmo
	dxx=dxx+MOocc(i)*(wfnderv(1,i)**2+wfnval(i)*wfnhess(1,1,i))
	dyy=dyy+MOocc(i)*(wfnderv(2,i)**2+wfnval(i)*wfnhess(2,2,i))
	dzz=dzz+MOocc(i)*(wfnderv(3,i)**2+wfnval(i)*wfnhess(3,3,i))
	rhotens3(1,1,1)=rhotens3(1,1,1)+MOocc(i)*( 3*wfnderv(1,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,1,i) )
	rhotens3(2,2,2)=rhotens3(2,2,2)+MOocc(i)*( 3*wfnderv(2,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,2,i) )
	rhotens3(3,3,3)=rhotens3(3,3,3)+MOocc(i)*( 3*wfnderv(3,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,3,i) )
	rhotens3(1,1,2)=rhotens3(1,1,2)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,2,i)+wfnderv(2,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,2,i) )
	rhotens3(1,1,3)=rhotens3(1,1,3)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,3,i)+wfnderv(3,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,3,i) )
	rhotens3(2,2,3)=rhotens3(2,2,3)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,3,i)+wfnderv(3,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,3,i) )
	rhotens3(1,2,2)=rhotens3(1,2,2)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,1,i)+wfnderv(1,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,1,i) ) !=2,2,1 exchange 1<->2 from (1,1,2) to derive this
	rhotens3(1,3,3)=rhotens3(1,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,1,i)+wfnderv(1,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,1,i) )
	rhotens3(2,3,3)=rhotens3(2,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,2,i)+wfnderv(2,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,2,i) )
end do
dxx=2D0*dxx
dyy=2D0*dyy
dzz=2D0*dzz
value=laplfac*(dxx+dyy+dzz) !laplfac is defined in settings.ini
rhotens3=rhotens3*2D0*laplfac
grad(1)=rhotens3(1,1,1)+rhotens3(1,2,2)+rhotens3(1,3,3)
grad(2)=rhotens3(1,1,2)+rhotens3(2,2,2)+rhotens3(2,3,3)
grad(3)=rhotens3(1,1,3)+rhotens3(2,2,3)+rhotens3(3,3,3)

! diff=1D-5
!Check of flapldx,dy,dz is correct!
! write(*,*) grad(:)
! difflapldx=(flapl(x+diff,y,z,'t')-flapl(x-diff,y,z,'t'))/(2*diff)
! difflapldy=(flapl(x,y+diff,z,'t')-flapl(x,y-diff,z,'t'))/(2*diff)
! difflapldz=(flapl(x,y,z+diff,'t')-flapl(x,y,z-diff,'t'))/(2*diff)
! write(*,*) difflapldx,difflapldy,difflapldz

!Check deviation between analytic and numerical solution
! write(*,*) rhotens3(1,1,1),rhotens3(1,2,2),rhotens3(1,3,3)
! diffrhodxxx=(flapl(x+diff,y,z,'x')-flapl(x-diff,y,z,'x'))/(2*diff)
! diffrhodxyy=(flapl(x+diff,y,z,'y')-flapl(x-diff,y,z,'y'))/(2*diff)
! diffrhodxzz=(flapl(x+diff,y,z,'z')-flapl(x-diff,y,z,'z'))/(2*diff)
! write(*,*) diffrhodxxx,diffrhodxyy,diffrhodxzz
! write(*,*)

!Check diagonal term with finite difference
! write(*,*) rhotens3(1,1,1),rhotens3(2,2,2),rhotens3(3,3,3)
! diffrhodxxx=(flapl(x+diff,y,z,'x')-flapl(x-diff,y,z,'x'))/(2*diff)
! diffrhodyyy=(flapl(x,y+diff,z,'y')-flapl(x,y-diff,z,'y'))/(2*diff)
! diffrhodzzz=(flapl(x,y,z+diff,'z')-flapl(x,y,z-diff,'z'))/(2*diff)
! write(*,*) diffrhodxxx,diffrhodyyy,diffrhodzzz

if (nEDFprims/=0) then
	call EDFrho(5,x,y,z,EDFdens,EDFgrad,EDFhess,EDFtens3)
	grad(1)=grad(1)+EDFtens3(1,1,1)+EDFtens3(1,2,2)+EDFtens3(1,3,3)
	grad(2)=grad(2)+EDFtens3(1,1,2)+EDFtens3(2,2,2)+EDFtens3(2,3,3)
	grad(3)=grad(3)+EDFtens3(1,1,3)+EDFtens3(2,2,3)+EDFtens3(3,3,3)
end if

!Don't consider Laplacian currently
if (itype==2) then !Calculate Hessian of laplacian
end if
end subroutine




!!----------------- Calculate Lagrangian kinetic G(r)
!idir=0/1/2/3 means total/x/y/z kinetic energy
real*8 function Lagkin(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo)
integer idir
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
lagkin=0D0
if (idir==0) then
	do imo=1,nmo
		lagkin=lagkin+MOocc(imo)*sum(wfnderv(:,imo)**2)
	end do
else
	do imo=1,nmo
		lagkin=lagkin+MOocc(imo)*wfnderv(idir,imo)**2
	end do
end if
lagkin=lagkin/2D0
end function




!!------------- Calculate Hamiltonian kinetic K(r)
!idir=0/1/2/3 means total/X/Y/Z kinetic energy
real*8 function Hamkin(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
integer idir
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
if (idir==0) then
	hamx=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
	hamy=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
	hamz=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
	Hamkin=hamx+hamy+hamz
else if (idir==1) then
	Hamkin=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
else if (idir==2) then
	Hamkin=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
else if (idir==3) then
	Hamkin=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
end if
Hamkin=-Hamkin/2D0
end function




!!------- Calculate analytic gradient vector of energy density (i.e. negative of K(r))
subroutine energydens_grad(x,y,z,grad)
real*8 x,y,z,grad(3),wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens(3,3,3,nmo)
call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens)
grad=0
do imo=1,nmo
	wfnlapl=wfnhess(1,1,imo)+wfnhess(2,2,imo)+wfnhess(3,3,imo)
	gx=wfnderv(1,imo)*wfnlapl+wfnval(imo)*(wfntens(1,1,1,imo)+wfntens(1,2,2,imo)+wfntens(1,3,3,imo))
	gy=wfnderv(2,imo)*wfnlapl+wfnval(imo)*(wfntens(1,1,2,imo)+wfntens(2,2,2,imo)+wfntens(2,3,3,imo))
	gz=wfnderv(3,imo)*wfnlapl+wfnval(imo)*(wfntens(1,1,3,imo)+wfntens(2,2,3,imo)+wfntens(3,3,3,imo))
	grad(1)=grad(1)+gx*MOocc(imo)
	grad(2)=grad(2)+gy*MOocc(imo)
	grad(3)=grad(3)+gz*MOocc(imo)
end do
grad=grad/2
end subroutine
!!------- Calculate gradient norm of energy density (i.e. negative of K(r))
real*8 function energydens_grdn(x,y,z)
real*8 x,y,z,grad(3)
call energydens_grad(x,y,z,grad)
energydens_grdn=dsqrt(sum(grad**2))
end function
!!------- Calculate Laplacian of energy density (i.e. negative of K(r))
real*8 function energydens_lapl(x,y,z)
real*8 x,y,z,gradaddx(3),gradminx(3),gradaddy(3),gradminy(3),gradaddz(3),gradminz(3)
diff=1D-5
denom=2*diff
call energydens_grad(x+diff,y,z,gradaddx)
call energydens_grad(x-diff,y,z,gradminx)
call energydens_grad(x,y+diff,z,gradaddy)
call energydens_grad(x,y-diff,z,gradminy)
call energydens_grad(x,y,z+diff,gradaddz)
call energydens_grad(x,y,z-diff,gradminz)
xlapl=(gradaddx(1)-gradminx(1))/denom
ylapl=(gradaddy(2)-gradminy(2))/denom
zlapl=(gradaddz(3)-gradminz(3))/denom
energydens_lapl=xlapl+ylapl+zlapl
end function




!!--------- Calculate electron linear momentum density in 3D representation
!idir=0/1/2/3 means magnitude/x/y/z component
real*8 function elemomdens(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),comp(0:3)
integer idir
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
comp=0
do imo=1,nmo
	comp(1:3)=comp(1:3)-MOocc(imo)*wfnval(imo)*wfnderv(1:3,imo)
end do
comp(0)=sum(comp(1:3)**2)
elemomdens=comp(idir)
end function




!!--------- Calculate magnetic dipole moment density
!idir=0/1/2/3 means magnitude/x/y/z component
real*8 function magmomdens(x,y,z,idir)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),comp(0:3)
integer idir
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
comp=0
do imo=1,nmo
	tmpx=-wfnval(imo)*(y*wfnderv(3,imo)-z*wfnderv(2,imo))
	tmpy=-wfnval(imo)*(z*wfnderv(1,imo)-x*wfnderv(3,imo))
	tmpz=-wfnval(imo)*(x*wfnderv(2,imo)-y*wfnderv(1,imo))
	comp(1)=comp(1)+MOocc(imo)*tmpx
	comp(2)=comp(2)+MOocc(imo)*tmpy
	comp(3)=comp(3)+MOocc(imo)*tmpz
end do
comp(0)=sum(comp(1:3)**2)
magmomdens=comp(idir)
end function




!!----- vdW potential (itype=1), repulsion potential (itype=2) and dispersion potential (itype=3)
!UFF parameters are used. The unit is kcal/mol
real*8 function vdwpotfunc(x,y,z,itype)
integer itype
real*8 x,y,z,parmA(ncenter),parmB(ncenter),UFF_A(103),UFF_B(103),tvec(3)

call defineUFFparm(UFF_A,UFF_B)
do iatm=1,ncenter
    parmA(iatm)=UFF_A(a(iatm)%index)
    parmB(iatm)=UFF_B(a(iatm)%index)
end do
parmAj=UFF_A(ivdwprobe)
parmBj=UFF_B(ivdwprobe)

repul=0
disp=0
if (ifPBC==0) then
    do iatm=1,ncenter
        dist=dsqrt( (a(iatm)%x-x)**2 + (a(iatm)%y-y)**2 + (a(iatm)%z-z)**2 )*b2a
        if (dist>25) cycle
        if (dist==0) dist=1E-10 !Avoid singular when grid point happens at nucleus
	    Dij=dsqrt(parmA(iatm)*parmAj) !Well depth
	    Xij=dsqrt(parmB(iatm)*parmBj) !vdW distance
	    repul=repul+Dij*(Xij/dist)**12 !Repulsion
	    disp=disp-2*Dij*(Xij/dist)**6 !Dispersion
	    !disp=disp-2*Dij*(Xij/dist)**3 !Another dispersion definition, this can result in better effect on coloring dispersion potential on vdW surface
    end do
else
    call getpointcell(x,y,z,ic,jc,kc)
    do icell=ic-PBCnx,ic+PBCnx
        do jcell=jc-PBCny,jc+PBCny
            do kcell=kc-PBCnz,kc+PBCnz
                call tvec_PBC(icell,jcell,kcell,tvec)
                do iatm=1,ncenter
					atmx=a(iatm)%x+tvec(1)
					atmy=a(iatm)%y+tvec(2)
					atmz=a(iatm)%z+tvec(3)
					dist=dsqrt( (atmx-x)**2 + (atmy-y)**2 + (atmz-z)**2 )*b2a
					if (dist>25) cycle
					if (dist==0) dist=1E-10 !Avoid singular when grid point happens at nucleus
					Dij=dsqrt(parmA(iatm)*parmAj) !Well depth
					Xij=dsqrt(parmB(iatm)*parmBj) !vdW distance
					repul=repul+Dij*(Xij/dist)**12 !Repulsion
					disp=disp-2*Dij*(Xij/dist)**6 !Dispersion
                end do
            end do
        end do
    end do
end if
if (itype==1) vdwpotfunc=repul+disp
if (itype==2) vdwpotfunc=repul
if (itype==3) vdwpotfunc=disp
end function




!!------ Calculate van der Waals potential, its gradient and Hessian matrix
!itype=1 Only calculate value and gradient, not Hessian
!itype=2 Calculate value, gradient and Hessian
subroutine calchessmat_vdWpot(itype,x,y,z,value,grad,hess)
integer itype
real*8 x,y,z,value,grad(3),hess(3,3),distvec(3),distgrad(3),disthess(3,3),xyzA(3),rvec(3),tmpgrad(3),tmphess(3,3),tvec(3)
real*8 parmA(ncenter),parmB(ncenter),UFF_A(103),UFF_B(103)

call defineUFFparm(UFF_A,UFF_B)
do iatm=1,ncenter
    parmA(iatm)=UFF_A(a(iatm)%index)
    parmB(iatm)=UFF_B(a(iatm)%index)
end do
parmAj=UFF_A(ivdwprobe)
parmBj=UFF_B(ivdwprobe)

value=0
grad=0
hess=0
rvec(1)=x;rvec(2)=y;rvec(3)=z

ilow=1;ihigh=1
jlow=1;jhigh=1
klow=1;khigh=1
if (ifPBC>0) then
	call getpointcell(x,y,z,ic,jc,kc)
    ilow=ic-PBCnx
    ihigh=ic+PBCnx
    jlow=jc-PBCny
    jhigh=jc+PBCny
    klow=kc-PBCnz
    khigh=kc+PBCnz
end if
do icell=ilow,ihigh
	do jcell=jlow,jhigh
		do kcell=klow,khigh
			tvec=0
            if (ifPBC>0) call tvec_PBC(icell,jcell,kcell,tvec)
			do iatm=1,ncenter
				!Get (R_A - r) vector and norm
				xyzA(1)=a(iatm)%x+tvec(1)
				xyzA(2)=a(iatm)%y+tvec(2)
				xyzA(3)=a(iatm)%z+tvec(3)
				distvec(:)=(xyzA(:)-rvec(:))*b2a
				dist=dsqrt(sum(distvec**2))
    
				if (dist>25) cycle
				if (dist==0) dist=1E-10 !Avoid singular when grid point happens at nucleus
				eps=dsqrt(parmA(iatm)*parmAj) !Well depth (epsilon)
				R0=dsqrt(parmB(iatm)*parmBj) !vdW distance (R0)
				R0_6=R0**6
				R0_12=R0_6**2
    
				tmpval = R0_12*dist**(-12) - 2*R0_6*dist**(-6)
				value = value + eps*tmpval
    
				!Get gradient of |R_A - r|
				distgrad(:)=-distvec(:)/dist
				!Calculate gradient of vdW potential
				tmp1= R0_6*dist**(-7) - R0_12*dist**(-13)
				grad(:) = grad(:) + eps*distgrad(:)*tmp1
    
				if (itype==2) then
					!Get Hessian of |R_A - r|
					dist3=dist**3
					do i=1,3
						do j=1,3
							disthess(i,j)=-distvec(i)*distvec(j)/dist3
							if (i==j) disthess(i,j)=disthess(i,j)+1D0/dist
						end do
					end do
					!Calculate Hessian of vdW potential
					tmp2=-7*R0_6*dist**(-8) + 13*R0_12*dist**(-14)
					do i=1,3
						do j=1,3
							tmphess(i,j) = disthess(i,j)*tmp1 + distgrad(i)*distgrad(j)*tmp2
						end do
					end do
					hess(:,:) = hess(:,:) + eps*tmphess(:,:)
				end if
			end do
		end do
	end do
end do

!Gradient and Hessian obtained above are derivatives w.r.t. Angstrom, convert to w.r.t. Bohr
grad=grad*12*b2a
hess=hess*12*b2a*b2a
end subroutine




!!------- Calculate sign(lambda2)*rho, this is a wrapper used to convert subroutine to function form
real*8 function signlambda2rho(x,y,z)
real*8 x,y,z,sl2r,RDG,rho
call signlambda2rho_RDG(x,y,z,sl2r,RDG)
signlambda2rho=sl2r
end function
!!------- Calculate sign(lambda2)*rho and RDG at the same time. rho and grad can also be returned if present
subroutine signlambda2rho_RDG(x,y,z,sl2r,RDG,rho,grad)
use util
real*8 x,y,z,elerho,sl2r,RDG
real*8 eigvecmat(3,3),eigval(3),elehess(3,3),elegrad(3) !Hessian of electron density
real*8,optional :: rho,grad(3)

call calchessmat_dens(2,x,y,z,elerho,elegrad,elehess)
call diagmat(elehess,eigvecmat,eigval,100,1D-10)
call sort(eigval)
if (eigval(2)==0D0) then !When eigval(2)==0D0, eigval(2)/abs(eigval(2)) can't be calculated, elerho generally will be zero, so sign is not important
	sl2r=elerho
else
	sl2r=elerho*eigval(2)/abs(eigval(2))
end if
sumgrad2=sum(elegrad(:)**2)
if (RDG_maxrho/=0D0.and.elerho>=RDG_maxrho) then
	RDG=100D0
!This occurs at distant region when exponent cutoff is used, the actual value should be very large. In order to avoid denominator become zero, we set it artifically to a big value
else if (sumgrad2==0D0.or.elerho==0D0) then
	RDG=999D0
else
	RDG=0.161620459673995D0*dsqrt(sumgrad2)/elerho**(4D0/3D0) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
end if
if (present(rho)) rho=elerho
if (present(grad)) grad=elegrad
end subroutine




!!------- Calculate ELF or LOL or Strong Covalent Interaction (SCI) and variants
!label="ELF" or "LOL" or "SCI"
real*8 function ELF_LOL(x,y,z,label)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo)
real*8 D,Dh,gradrho(3),gradrhoa(3),gradrhob(3),rho,rhoa,rhob,rhospin,MOoccnow
real*8 :: Fc=2.871234000D0 ! Thomas-Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 ! Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194
character label*3

!Calculate Tsirelson version of ELF and LOL, which are only dependent on electron density
!Since rho, nebla-rho, nebla^2-rho support EDF, these functions also support EDF
if (ELFLOL_type==1) then
	rho=fdens(x,y,z)
	if (wfntype==0.or.wfntype==3) then !Closed-shell cases
		Dh=Fc*rho**(5D0/3D0) !Thomas-Fermi kinetic energy
	else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Open shell cases
		rhospin=fspindens(x,y,z,'s') !rhospin=rhoa-rhob, rho=rhoa+rhob
		rhoa=(rhospin+rho)/2D0
		rhob=(rho-rhospin)/2D0
		Dh=Fc_pol*(rhoa**(5D0/3D0)+rhob**(5D0/3D0)) !Thomas-Fermi kinetic energy
	end if
	if (label=="ELF") then !The ELF defined by Tsirelson, CPL, 351, 142
		!Restrictly speaking, the kinetic energy expansion should be replace by polarized form for open-shell
		D=Dh-(1/9D0)*fgrad(x,y,z,'t')**2/rho+(1/6D0)*flapl(x,y,z,'t')
		ELF_LOL=1/(1+(D/Dh)**2)
	else if (label=="LOL") then !The LOL defined by Tsirelson, Acta Cryst. B58, 780-785 (2002)
		D=Dh+(1/72D0)*fgrad(x,y,z,'t')**2/rho+(1/6D0)*flapl(x,y,z,'t')
		t=Dh/D
		ELF_LOL=1D0/(1D0/t+1)
	end if
	return
end if

!Calculate ELF, LOL etc. based on wavefunction
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
D=0D0
rho=0D0
rhoa=0D0
rhob=0D0
gradrho=0D0
gradrhoa=0D0
gradrhob=0D0
if (label=="ELF".or.label=="SCI") then !----- Calculate ELF or SCI based on wavefunction
	if (wfntype==0.or.wfntype==3) then !Closed-shell case
		do i=1,nmo
			rho=rho+MOocc(i)*wfnval(i)**2
			gradrho(:)=gradrho(:)+2D0*MOocc(i)*wfnval(i)*wfnderv(:,i)
			D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Actual kinetic energy density
		end do		
		D=D/2D0
		if (iKEDsel/=0) D=KED(x,y,z,iKEDsel) !Special case proposed by LSB, use other KED instead of exact KED
		if (rho/=0D0) D=D-sum(gradrho(:)**2)/rho/8D0 !Pauli kinetic energy density
		Dh=Fc*rho**(5D0/3D0) !Thomas-Fermi uniform electron gas kinetic energy density
	else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Spin-polarized case
		do i=1,nmo
			MOoccnow=MOocc(i)
			if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !Double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
			if (MOtype(i)==1.or.MOtype(i)==0) then
				rhoa=rhoa+MOoccnow*wfnval(i)**2
				gradrhoa(:)=gradrhoa(:)+2D0*MOoccnow*wfnval(i)*wfnderv(:,i)
			end if
			if (MOtype(i)==2.or.MOtype(i)==0) then
				rhob=rhob+MOoccnow*wfnval(i)**2
				gradrhob(:)=gradrhob(:)+2D0*MOoccnow*wfnval(i)*wfnderv(:,i)
			end if
			D=D+MOocc(i)*(sum(wfnderv(:,i)**2)) !Actual kinetic energy density
		end do
		D=D/2D0
		if (iKEDsel/=0) D=KED(x,y,z,iKEDsel) !Special case proposed by LSB, use other KED instead of exact KED
		if (rhoa/=0D0) D=D-sum(gradrhoa(:)**2)/rhoa/8 !Pauli kinetic energy density
		if (rhob/=0D0) D=D-sum(gradrhob(:)**2)/rhob/8
		Dh=Fc_pol*(rhoa**(5D0/3D0)+rhob**(5D0/3D0)) !Thomas-Fermi uniform electron gas kinetic energy density
	end if
	if (label=="ELF") then
		if (ELFLOL_type==0) then !Conventional ELF
			if (ELF_addminimal==1) D=D+1D-5 !Add 1D-5 to avoid D become zero, leading to unity in infinite
			ELF_LOL=1/(1+(D/Dh)**2)
		else if (ELFLOL_type==2) then !New ELF formalism defined by Tian Lu
			if (ELF_addminimal==1) D=D+1D-5 !add 1D-5 to avoid D become zero, leading to unity in infinite
			ELF_LOL=1/(1+(D/Dh))
		else if (ELFLOL_type==3) then !Only get D/D0 term for special usage
			ELF_LOL=D/Dh
		end if
	else if (label=="SCI") then
		if (ELF_addminimal==1) D=D+1D-5
		ELF_LOL=1/(D/Dh)
	end if

else if (label=="LOL") then !----- Calculate LOL based on wavefunction
	t=0D0
	if (wfntype==0.or.wfntype==3) then !Closed-shell case
		do i=1,nmo !Store actual kinetic energy density to t first
			rho=rho+MOocc(i)*wfnval(i)**2
			t=t+MOocc(i)*(sum(wfnderv(:,i)**2))
		end do
		t=t/2D0
		Dh=Fc*rho**(5D0/3D0)
	else if (wfntype==1.or.wfntype==2.or.wfntype==4) then !Spin-polarized case
		do i=1,nmo !Store actual kinetic energy to t first
			MOoccnow=MOocc(i)
			if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !Doubly occupied, present when wfntype==2 (ROHF), alpha and beta get half part
			if (MOtype(i)==1.or.MOtype(i)==0) rhoa=rhoa+MOoccnow*wfnval(i)**2
			if (MOtype(i)==2.or.MOtype(i)==0) rhob=rhob+MOoccnow*wfnval(i)**2
			t=t+MOocc(i)*(sum(wfnderv(:,i)**2))
		end do
		t=t/2D0
		Dh=Fc_pol*(rhoa**(5D0/3D0)+rhob**(5D0/3D0))
	end if
	!--------- A new definition of LOL, however the value range is not as good as LOL
	! ELF_LOL=t-Dh
	! if (ELF_LOL>0) then
	! 	ELF_LOL=1D0/(1D0+1D0/ELF_LOL)
	! else if (ELF_LOL<0) then
	! 	ELF_LOL=-1D0/(1D0+1D0/abs(ELF_LOL))
	! end if
	!-------------
	!If there is very long distance between molecule and current point, t (above) is zero,
	!and t=Dh/t is also zero (because rho converges faster), but can't be calculate directly, so simply skip
	if (t/=0D0) t=Dh/t
	if (ELFLOL_type==0) ELF_LOL=1D0/(1D0/t+1) !namely t/(1+t). This is default case
	if (ELFLOL_type==2) ELF_LOL=1D0/((1D0/t)**2+1) !New form defined by Tian Lu
end if
end function




!!------ Calculate ELF/LOL, its gradient and Hessian matrix (not supported yet)
!funsel="ELF" or "LOL"
!itype=1 Only calculate value and gradient, not Hessian
!itype=2 Calculate value, gradient and Hessian (Not available)
subroutine calchessmat_ELF_LOL(itype,x,y,z,value,grad,hess,funsel)
use util
integer itype
real*8 x,y,z,value,grad(3),hess(3,3),MOoccnow
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo)
real*8 rho,gradrho(3),hessrho(3,3),Dh,gradDh(3),hessDh(3,3),Ts,gradTs(3),hessTs(3,3),Wei,gradWei(3),hessWei(3,3)
real*8 rhoa,rhob,gradrhoa(3),gradrhob(3),hessrhoa(3,3),hessrhob(3,3),Dha,Dhb,gradDha(3),gradDhb(3),Weia,Weib,gradWeia(3),gradWeib(3)
real*8 :: Fc=2.871234000D0 !Fermi constant = (3/10)*(3*Pi^2)**(2/3) = 2.871234, 1/2.871234=0.34828
real*8 :: Fc_pol=4.557799872D0 !Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3) = 4.5578, 1/4.5578=0.2194
real*8 :: corrELF=1D-5
character funsel*3

if (itype==1) call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess) !Get Hessian of GTF, needn't 3-order tensor
if (itype==2) call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)

!Spin-unpolarized case
if (wfntype==0.or.wfntype==3) then
	rho=0D0
	gradrho=0D0
	Ts=0D0
	gradTs=0D0
	do i=1,nmo
		rho=rho+MOocc(i)*wfnval(i)**2
		gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
		Ts=Ts+MOocc(i)*(sum(wfnderv(:,i)**2))
	end do
	gradrho=2*gradrho
	Ts=Ts/2D0
	Dh=Fc*rho**(5D0/3D0)
	gradDh(:)=5D0/3D0*Fc*rho**(2D0/3D0)*gradrho(:)
	do i=1,nmo
		gradTs(1)=gradTs(1)+MOocc(i)*(wfnderv(1,i)*wfnhess(1,1,i)+wfnderv(2,i)*wfnhess(1,2,i)+wfnderv(3,i)*wfnhess(1,3,i))
		gradTs(2)=gradTs(2)+MOocc(i)*(wfnderv(2,i)*wfnhess(2,2,i)+wfnderv(1,i)*wfnhess(2,1,i)+wfnderv(3,i)*wfnhess(2,3,i))
		gradTs(3)=gradTs(3)+MOocc(i)*(wfnderv(3,i)*wfnhess(3,3,i)+wfnderv(1,i)*wfnhess(3,1,i)+wfnderv(2,i)*wfnhess(3,2,i))
	end do
	if (funsel=="ELF") then
		!Calculate Hessian for rho
		hessrho(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
		hessrho(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
		hessrho(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
		hessrho(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
		hessrho(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
		hessrho(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
		hessrho(2,1)=hessrho(1,2)
		hessrho(3,1)=hessrho(1,3)
		hessrho(3,2)=hessrho(2,3)
		!Calculate Weizsacker functional and its derivatives
		Wei=sum(gradrho(:)**2)/8D0/rho
		D=Ts-Wei+corrELF
		chi=D/Dh
		value=1D0/(1D0+chi**2)
		gradWei(1)= 0.25D0/rho*( gradrho(1)*hessrho(1,1)+gradrho(2)*hessrho(1,2)+gradrho(3)*hessrho(1,3) ) - wei/rho*gradrho(1)
		gradWei(2)= 0.25D0/rho*( gradrho(2)*hessrho(2,2)+gradrho(1)*hessrho(2,1)+gradrho(3)*hessrho(2,3) ) - wei/rho*gradrho(2)
		gradWei(3)= 0.25D0/rho*( gradrho(3)*hessrho(3,3)+gradrho(2)*hessrho(3,2)+gradrho(1)*hessrho(3,1) ) - wei/rho*gradrho(3)
		chidx=(gradTs(1)-gradWei(1))/Dh - gradDh(1)/Dh**2 *D
		chidy=(gradTs(2)-gradWei(2))/Dh - gradDh(2)/Dh**2 *D
		chidz=(gradTs(3)-gradWei(3))/Dh - gradDh(3)/Dh**2 *D
		grad(1)=-2D0*chi/(1+chi**2)**2 * chidx
		grad(2)=-2D0*chi/(1+chi**2)**2 * chidy
		grad(3)=-2D0*chi/(1+chi**2)**2 * chidz
	else if (funsel=="LOL") then
		value=1D0/(1D0+Ts/Dh)
		tmp=-1D0/Dh/(1D0+Ts/Dh)**2
		grad(1)=tmp*(gradTs(1)-Ts/Dh*gradDh(1))
		grad(2)=tmp*(gradTs(2)-Ts/Dh*gradDh(2))
		grad(3)=tmp*(gradTs(3)-Ts/Dh*gradDh(3))
	end if
!Spin-polarized case
else if (wfntype==1.or.wfntype==2.or.wfntype==4) then
	rhoa=0D0
	rhob=0D0
	gradrhoa=0D0
	gradrhob=0D0
	Ts=0D0
	gradTs=0D0
	do i=1,nmo
		MOoccnow=MOocc(i)
		if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0
		if (MOtype(i)==1.or.MOtype(i)==0) then
			rhoa=rhoa+MOoccnow*wfnval(i)**2
			gradrhoa(:)=gradrhoa(:)+MOoccnow*wfnval(i)*wfnderv(:,i)
		end if
		if (MOtype(i)==2.or.MOtype(i)==0) then
			rhob=rhob+MOoccnow*wfnval(i)**2
			gradrhob(:)=gradrhob(:)+MOoccnow*wfnval(i)*wfnderv(:,i)			
		end if
		Ts=Ts+MOocc(i)*(sum(wfnderv(:,i)**2))
	end do
	gradrhoa=2*gradrhoa
	gradrhob=2*gradrhob
	Ts=Ts/2D0
	Dha=Fc_pol*rhoa**(5D0/3D0)
	Dhb=Fc_pol*rhob**(5D0/3D0)
	Dh=Dha+Dhb
	gradDha(:)=5D0/3D0*Fc_pol*rhoa**(2D0/3D0)*gradrhoa(:)
	gradDhb(:)=5D0/3D0*Fc_pol*rhob**(2D0/3D0)*gradrhob(:)
	gradDh=gradDha+gradDhb
	do i=1,nmo
		gradTs(1)=gradTs(1)+MOocc(i)*(wfnderv(1,i)*wfnhess(1,1,i)+wfnderv(2,i)*wfnhess(1,2,i)+wfnderv(3,i)*wfnhess(1,3,i))
		gradTs(2)=gradTs(2)+MOocc(i)*(wfnderv(2,i)*wfnhess(2,2,i)+wfnderv(1,i)*wfnhess(2,1,i)+wfnderv(3,i)*wfnhess(2,3,i))
		gradTs(3)=gradTs(3)+MOocc(i)*(wfnderv(3,i)*wfnhess(3,3,i)+wfnderv(1,i)*wfnhess(3,1,i)+wfnderv(2,i)*wfnhess(3,2,i))
	end do
	
	if (funsel=="ELF") then
! 		Calculate Hessian for rho
		hessrhoa=0D0
		hessrhob=0D0
		do i=1,nmo
			MOoccnow=MOocc(i)
			if (MOtype(i)==0) MOoccnow=MOocc(i)/2D0 !double occupied, present when wfntype==2 (ROHF), alpha and beta get half part
			if (MOtype(i)==1.or.MOtype(i)==0) then
				hessrhoa(1,1)=hessrhoa(1,1)+MOoccnow*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
				hessrhoa(2,2)=hessrhoa(2,2)+MOoccnow*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
				hessrhoa(3,3)=hessrhoa(3,3)+MOoccnow*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
				hessrhoa(1,2)=hessrhoa(1,2)+MOoccnow*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
				hessrhoa(2,3)=hessrhoa(2,3)+MOoccnow*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
				hessrhoa(1,3)=hessrhoa(1,3)+MOoccnow*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
			end if
			if (MOtype(i)==2.or.MOtype(i)==0) then
				hessrhob(1,1)=hessrhob(1,1)+MOoccnow*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
				hessrhob(2,2)=hessrhob(2,2)+MOoccnow*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
				hessrhob(3,3)=hessrhob(3,3)+MOoccnow*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
				hessrhob(1,2)=hessrhob(1,2)+MOoccnow*( wfnderv(1,i)*wfnderv(2,i)+wfnhess(1,2,i)*wfnval(i) )
				hessrhob(2,3)=hessrhob(2,3)+MOoccnow*( wfnderv(2,i)*wfnderv(3,i)+wfnhess(2,3,i)*wfnval(i) )
				hessrhob(1,3)=hessrhob(1,3)+MOoccnow*( wfnderv(1,i)*wfnderv(3,i)+wfnhess(1,3,i)*wfnval(i) )
			end if
		end do
		hessrhoa=hessrhoa*2
		hessrhob=hessrhob*2
		hessrhoa(2,1)=hessrhoa(1,2)
		hessrhoa(3,1)=hessrhoa(1,3)
		hessrhoa(3,2)=hessrhoa(2,3)
		hessrhob(2,1)=hessrhob(1,2)
		hessrhob(3,1)=hessrhob(1,3)
		hessrhob(3,2)=hessrhob(2,3)
! 		!Calculate Weizsacker functional and its derivatives
		Weia=sum(gradrhoa(:)**2)/8D0/rhoa
		Weib=sum(gradrhob(:)**2)/8D0/rhob
		Wei=Weia+Weib
		D=Ts-Wei+corrELF
		chi=D/Dh
		value=1D0/(1D0+chi**2)
		gradWeia(1)= 0.25D0/rhoa*( gradrhoa(1)*hessrhoa(1,1)+gradrhoa(2)*hessrhoa(1,2)+gradrhoa(3)*hessrhoa(1,3) ) - weia/rhoa*gradrhoa(1)
		gradWeia(2)= 0.25D0/rhoa*( gradrhoa(2)*hessrhoa(2,2)+gradrhoa(1)*hessrhoa(2,1)+gradrhoa(3)*hessrhoa(2,3) ) - weia/rhoa*gradrhoa(2)
		gradWeia(3)= 0.25D0/rhoa*( gradrhoa(3)*hessrhoa(3,3)+gradrhoa(2)*hessrhoa(3,2)+gradrhoa(1)*hessrhoa(3,1) ) - weia/rhoa*gradrhoa(3)
		gradWeib(1)= 0.25D0/rhob*( gradrhob(1)*hessrhob(1,1)+gradrhob(2)*hessrhob(1,2)+gradrhob(3)*hessrhob(1,3) ) - weib/rhob*gradrhob(1)
		gradWeib(2)= 0.25D0/rhob*( gradrhob(2)*hessrhob(2,2)+gradrhob(1)*hessrhob(2,1)+gradrhob(3)*hessrhob(2,3) ) - weib/rhob*gradrhob(2)
		gradWeib(3)= 0.25D0/rhob*( gradrhob(3)*hessrhob(3,3)+gradrhob(2)*hessrhob(3,2)+gradrhob(1)*hessrhob(3,1) ) - weib/rhob*gradrhob(3)
		gradWei=gradWeia+gradWeib
		chidx=(gradTs(1)-gradWei(1))/Dh - gradDh(1)/Dh**2 *D
		chidy=(gradTs(2)-gradWei(2))/Dh - gradDh(2)/Dh**2 *D
		chidz=(gradTs(3)-gradWei(3))/Dh - gradDh(3)/Dh**2 *D
		grad(1)=-2D0*chi/(1+chi**2)**2 * chidx
		grad(2)=-2D0*chi/(1+chi**2)**2 * chidy
		grad(3)=-2D0*chi/(1+chi**2)**2 * chidz
	else if (funsel=="LOL") then
		value=1D0/(1D0+Ts/Dh)
		tmp=-1D0/Dh/(1D0+Ts/Dh)**2
		grad(1)=tmp*(gradTs(1)-Ts/Dh*gradDh(1))
		grad(2)=tmp*(gradTs(2)-Ts/Dh*gradDh(2))
		grad(3)=tmp*(gradTs(3)-Ts/Dh*gradDh(3))
	end if
end if
end subroutine



!!-------- Calculate average local ionization energy (ALIE)
real*8 function avglocion(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
avglocion=0D0
rho=0D0
do i=1,nmo
	avglocion=avglocion+abs(MOene(i))*MOocc(i)*wfnval(i)**2
	rho=rho+MOocc(i)*wfnval(i)**2 !Calculate rho
end do
if (rho==0D0) then
	avglocion=0D0 !Avoid at distant region rho becomes zero when exponent cutoff is used
else
	avglocion=avglocion/rho
end if
end function

!!-------- Calculate average local ionization energy and meantime decompose it to occupied orbital contributions
subroutine avglociondecomp(ifileid,x,y,z)
real*8 x,y,z,wfnval(nmo)
integer ifileid
character orbtype*2
call orbderv(1,1,nmo,x,y,z,wfnval)
totALIE=0D0
rho=0D0
do i=1,nmo
	totALIE=totALIE+abs(MOene(i))*MOocc(i)*wfnval(i)**2
	rho=rho+MOocc(i)*wfnval(i)**2 !Calculate rho
end do
if (rho==0D0) then
	totALIE=0D0 !Avoid at distant region rho becomes zero when exponent cutoff is used
else
	totALIE=totALIE/rho
end if
write(ifileid,"(' Average local ionization energy:',f16.10,' a.u.')") totALIE
write(ifileid,*) "Contribution of each orbital to average local ionization energy (a.u.):"
do i=1,nmo
	if (MOtype(i)==0) orbtype="AB"
	if (MOtype(i)==1) orbtype="A "
	if (MOtype(i)==2) orbtype="B "
	write(ifileid,"(' Orbital',i6,'  Ene:',f12.6,'  Occ:',f5.2,'  Type:',a,'  Contribution:',f12.6)") i,MOene(i),MOocc(i),orbtype,abs(MOene(i))*MOocc(i)*wfnval(i)**2/rho
end do
end subroutine




!!-------- Calculate local electron affinity
!Since virtual orbitals are involved, such as .fch/.molden/.gms must be used
real*8 function loceleaff(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
loceleaff=0D0
rho=0D0
do i=1,nmo
	if (MOocc(i)==0) then !Only cycle unoccupied orbitals
		loceleaff=loceleaff-MOene(i)*wfnval(i)**2 !Don't need to multiply "occupation number", because ROHF is not allowed, so all orbitals have the same type
		rho=rho+wfnval(i)**2 !Calculate rho
	end if
end do
if (rho==0D0) then
	loceleaff=0D0 !Avoid at distant region rho become zero when exponent cutoff is used
else
	loceleaff=loceleaff/rho
end if
end function




!!-------- Calculate Local electron attachment energy
!Since virtual orbitals are involved, such as .fch/.molden/.gms must be used
real*8 function loceleatt(x,y,z)
real*8 x,y,z,wfnval(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
loceleatt=0D0
rho=0D0
do i=1,nmo
	tmp=wfnval(i)**2
    rho=rho+MOocc(i)*tmp
	if (MOocc(i)==0.and.MOene(i)<0) then !Only cycle unoccupied orbitals with negative energy
		loceleatt=loceleatt+MOene(i)*tmp
	end if
end do
if (wfntype==0) loceleatt=loceleatt*2
if (rho==0D0) then
	loceleatt=0D0 !Avoid at distant region rho become zero when exponent cutoff is used
else
	loceleatt=loceleatt/rho
end if
end function




!!------ Approximate form of DFT linear response kernel for closed-shell, X(r1,r2), see Eq.3 of PCCP,14,3960. Only applied to the case wfntype=0
! r1 is taken as reference point and determined by refx,refy,refz. x,y,z in the argument is the coordinate of r2
real*8 function linrespkernel(x,y,z)
real*8 x,y,z,orbvalr1(nmo),orbvalr2(nmo)
call orbderv(1,1,nmo,refx,refy,refz,orbvalr1)
call orbderv(1,1,nmo,x,y,z,orbvalr2)
linrespkernel=0D0
do imo=1,nmo !Cycle occupied MOs
	if (nint(MOocc(imo))==2D0) then
		do jmo=idxHOMO+1,nmo !Cycle unoccupied MOs
			if (nint(MOocc(jmo))==0D0) linrespkernel=linrespkernel+orbvalr1(imo)*orbvalr1(jmo)*orbvalr2(jmo)*orbvalr2(imo)/(MOene(imo)-MOene(jmo))
		end do
	end if
end do
linrespkernel=linrespkernel*4
end function

!!------------- Calculate Exchange-correlation density, correlation hole and correlation factor
!rfx,rfy,rfz is reference point (commonly use refx,refy,refz in defvar module), namely r1
!x,y,z in the argument is the coordinate of r2
!Calculate which function is controlled by "pairfunctype" in settings.ini, correlation type is determined by "paircorrtype" in settings.ini
real*8 function pairfunc(rfx,rfy,rfz,x,y,z)
real*8 rfx,rfy,rfz,x,y,z,orbvalr1(nmo),orbvalr2(nmo)
call orbderv(1,1,nmo,rfx,rfy,rfz,orbvalr1)
call orbderv(1,1,nmo,x,y,z,orbvalr2)
!Calculate alpha and beta density at r1 and r2
adensr1=0D0
bdensr1=0D0
adensr2=0D0
bdensr2=0D0
do i=1,nmo
	if (MOtype(i)==0) then
		adensr1=adensr1+MOocc(i)/2*orbvalr1(i)**2
		adensr2=adensr2+MOocc(i)/2*orbvalr2(i)**2
		bdensr1=bdensr1+MOocc(i)/2*orbvalr1(i)**2
		bdensr2=bdensr2+MOocc(i)/2*orbvalr2(i)**2
	else if (MOtype(i)==1) then
		adensr1=adensr1+MOocc(i)*orbvalr1(i)**2
		adensr2=adensr2+MOocc(i)*orbvalr2(i)**2
	else if (MOtype(i)==2) then
		bdensr1=bdensr1+MOocc(i)*orbvalr1(i)**2
		bdensr2=bdensr2+MOocc(i)*orbvalr2(i)**2
	end if
end do
totdensr1=adensr1+bdensr1
totdensr2=adensr2+bdensr2

ntime=1
if (pairfunctype==12) ntime=2 !Will need both alpha and beta information (aXCdens and bXCdens), so process twice
do itime=1,ntime
	!Calculate exchange-correlation density first, and then calculate correlation hole and correlation factor
	!For RHF/ROHF, we calculate them as if present system is open-shell
	if (pairfunctype==1.or.pairfunctype==4.or.pairfunctype==7.or.pairfunctype==10.or.(pairfunctype==12.and.itime==1)) then !Set start and end index of alpha orbitals
		!Cycle alpha orbitals first to obtain aXCdens
		istart=1
		if (wfntype==0.or.wfntype==2.or.wfntype==3) then !RHF,ROHF,R-post-HF
			iend=nmo
		else if (wfntype==1.or.wfntype==4) then !UHF, U-post-HF
			do iend=nmo,1,-1
				if (MOtype(iend)==1) exit
			end do
		end if
	else if (pairfunctype==2.or.pairfunctype==5.or.pairfunctype==8.or.pairfunctype==11.or.(pairfunctype==12.and.itime==2)) then !Set start and end index of beta orbitals
		if (wfntype==0.or.wfntype==3) then !RHF,R-post-HF
			istart=1
			iend=nmo
		else if (wfntype==2) then !ROHF
			istart=1
			do iend=1,nmo
				if (MOtype(iend)==1) exit
			end do
			iend=iend-1
		else if (wfntype==1.or.wfntype==4) then !UHF, U-post-HF
			do istart=1,nmo
				if (MOtype(istart)==2) exit
			end do
			iend=nmo
			if (nint(nbelec)==0) iend=0 !less than istart, so below cycle will be skipped
		end if
	end if

	XCtmp=0D0 !Really X+C density
	Xtmp=0D0 !Only X density
	Ctmp=0D0 !Only C density
	do i=istart,iend
		occi=MOocc(i)
		if (MOtype(i)==0) occi=occi/2D0 !Split closed-shell orbital to spin orbital
		do j=istart,iend
			occj=MOocc(j)
			if (MOtype(j)==0) occj=occj/2D0
			tmpmul=orbvalr1(i)*orbvalr2(j)*orbvalr1(j)*orbvalr2(i)
			XCtmp=XCtmp-dsqrt(occi*occj)*tmpmul
			Xtmp=Xtmp-occi*occj*tmpmul
			Ctmp=Ctmp+(occi*occj-dsqrt(occi*occj))*tmpmul
		end do
	end do
		
	if (pairfunctype==1.or.pairfunctype==4.or.pairfunctype==7.or.pairfunctype==10.or.(pairfunctype==12.and.itime==1)) then
		if (paircorrtype==1) aXCdens=Xtmp
		if (paircorrtype==2) aXCdens=Ctmp
		if (paircorrtype==3) aXCdens=XCtmp
		acorrhole=aXCdens/adensr1
		acorrfac=acorrhole/adensr2
		if (pairfunctype==1) pairfunc=acorrhole
		if (pairfunctype==4) pairfunc=acorrfac
		if (pairfunctype==7) pairfunc=aXCdens
		if (pairfunctype==10) pairfunc=adensr1*adensr2+aXCdens
	else if (pairfunctype==2.or.pairfunctype==5.or.pairfunctype==8.or.pairfunctype==11.or.(pairfunctype==12.and.itime==2)) then
		if (paircorrtype==1) bXCdens=Xtmp
		if (paircorrtype==2) bXCdens=Ctmp
		if (paircorrtype==3) bXCdens=XCtmp
		bcorrhole=bXCdens/bdensr1
		bcorrfac=bcorrhole/bdensr2
		if (pairfunctype==2) pairfunc=bcorrhole
		if (pairfunctype==5) pairfunc=bcorrfac
		if (pairfunctype==8) pairfunc=bXCdens
		if (pairfunctype==11) pairfunc=bdensr1*bdensr2+bXCdens
	end if
end do
if (pairfunctype==12) pairfunc=adensr1*(adensr2+bdensr2)+aXCdens +bdensr1*(adensr2+bdensr2)+bXCdens
end function




!!------------- Calculate source function
real*8 function srcfunc(x,y,z,imode) !Default imode=1
real*8 x,y,z,denomin
integer imode
denomin=4*pi*dsqrt((x-refx)**2+(y-refy)**2+(z-refz)**2)
if (denomin==0D0) denomin=0.001D0
if (imode==1) srcfunc=-flapl(x,y,z,'t')/denomin !Used to study effect of laplacian everywhere on specific point
if (imode==2) srcfunc=-flapl(refx,refy,refz,'t')/denomin !Used to study effect of laplacian at specific point on everywhere
end function




!!------------- Calculate RDG with promolecular approximation
real*8 function RDGprodens(x,y,z)
real*8 x,y,z,elerho,elegrad(3)
call calchessmat_prodens(x,y,z,elerho,elegrad)
elegradnorm=dsqrt(sum(elegrad**2))
if ((RDGprodens_maxrho/=0D0.and.elerho>=RDGprodens_maxrho)) then
	RDGprodens=100D0
else if (elegradnorm==0D0.or.elerho==0D0) then
	RDGprodens=999D0
else
	RDGprodens=0.161620459673995D0*elegradnorm/elerho**(4D0/3D0)
end if
end function
!!------ Calculate sign(lambda2(r))*rho(r) with promolecular approximation
! This is a wrapper to convert subroutine to function form
real*8 function signlambda2rho_prodens(x,y,z)
real*8 x,y,z,sl2r,RDG
call signlambda2rho_RDG_prodens(x,y,z,sl2r,RDG)
signlambda2rho_prodens=sl2r
end function
!!------- Calculate sign(lambda2)*rho and RDG at the same time with promolecular approximation
subroutine signlambda2rho_RDG_prodens(x,y,z,sl2r,RDG)
use util
real*8 x,y,z,elerho,RDG,sl2r,sumgrad2
real*8 eigvecmat(3,3),eigval(3),elehess(3,3),elegrad(3)
call calchessmat_prodens(x,y,z,elerho,elegrad,elehess)
call diagmat(elehess,eigvecmat,eigval,100,1D-6)
call sort(eigval)
if (eigval(2)/=0D0) then
	sl2r=elerho*eigval(2)/abs(eigval(2)) !At nuclei of single atom system, hessian returned may be zero matrix
else
	sl2r=-elerho !Around nuclei, eigval(2)/abs(eigval(2)) always be negative
end if
! sl2r=elerho !Only obtain promolecular density
sumgrad2=sum(elegrad(:)**2)
if ((RDGprodens_maxrho/=0D0.and.elerho>=RDGprodens_maxrho).or.elerho==0D0) then
	RDG=100D0
else if (sumgrad2==0D0.or.elerho==0D0) then
	RDG=999D0
else
	RDG=0.161620459673995D0*dsqrt(sumgrad2)/elerho**(4D0/3D0)
end if
end subroutine




!!----- Calculate electron density, its gradient and Hessian matrix at x,y,z with promolecular approximation
!Currently only employed by RDG analysis!!!!!!!!!!!
!The YWT fitted atomic density used in this subroutine (for element <=Ar) in fact is poor
!Electron density and its gradient are always calculated, Hessian will be calculated when "elehess" is present
!Notice that global array "fragment" must be properly defined! Only the atoms in fragment will be taken into account
subroutine calchessmat_prodens(xin,yin,zin,elerho,elegrad,elehess)
use util
real*8 elerho,xin,yin,zin
real*8,optional :: elegrad(3),elehess(3,3)
real*8 rhoarr(200),tvec(3)
elerho=0D0
derx=0D0
dery=0D0
derz=0D0
dxx=0D0
dyy=0D0
dzz=0D0
dxy=0D0
dyz=0D0
dxz=0D0
idohess=0
if (present(elehess)) idohess=1

call getpointcell(xin,yin,zin,ic,jc,kc)
do icell=ic-PBCnx,ic+PBCnx
    do jcell=jc-PBCny,jc+PBCny
        do kcell=kc-PBCnz,kc+PBCnz
            call tvec_PBC(icell,jcell,kcell,tvec)
            do i=1,nfragatm
	            iatm=fragatm(i)
	            iele=a(iatm)%index
	            !rx=a(iatm)%x+tvec(1)-xin !Wrong code, older than 2022-Sep-18
	            !ry=a(iatm)%y+tvec(2)-yin
	            !rz=a(iatm)%z+tvec(3)-zin
	            rx=xin-tvec(1)-a(iatm)%x !Relative x
	            ry=yin-tvec(2)-a(iatm)%y
	            rz=zin-tvec(3)-a(iatm)%z
	            rx2=rx*rx
	            ry2=ry*ry
	            rz2=rz*rz
	            r2=rx2+ry2+rz2
	            r=dsqrt(r2)
	            if (iele<=18) then !H~Ar, use Weitao Yang's fitted parameters as original RDG paper
	                if (atomdenscut==1) then !Tight cutoff, for CHNO corresponding to cutoff at rho=0.00001
		                if (iele==1.and.r2>25D0) then !H, 6.63^2=43.9569. But this seems to be unnecessarily large, so I use 5^2=25
			                cycle
		                else if (iele==6.and.r2>58.6756D0) then !C, 7.66^2=58.6756
			                cycle
		                else if (iele==7.and.r2>43.917129D0) then !N, 6.627^2=43.917129
			                cycle
		                else if (iele==8.and.r2>34.9281D0) then !O, 5.91^2=34.9281
			                cycle
		                else if (r2>(2.5D0*vdwr(iele))**2) then !Other cases, larger than 2.5 times of its vdw radius will be skipped
			                cycle
		                end if
	                else if (atomdenscut==2) then !Medium cutoff, the result may be not as accurate as atomdenscut==1, but much more cheaper
		                if (r2>(2.2D0*vdwr(iele))**2) cycle
	                else if (atomdenscut==3) then !Loose cutoff, the most inaccurate
		                if (r2>(1.8D0*vdwr(iele))**2) cycle
	                else if (atomdenscut==4) then !Foolish cutoff, you need to know what you are doing
		                if (r2>(1.5D0*vdwr(iele))**2) cycle
	                end if
		            r2_1d5=r2**1.5D0
		            do iSTO=1,3
			            if (YWTatomcoeff(iele,iSTO)==0D0) cycle
			            expterm=YWTatomexp(iele,iSTO)
			            term=YWTatomcoeff(iele,iSTO)*dexp(-r/expterm)
			            elerho=elerho+term
			            if (r==0D0) cycle !Derivative of STO at nuclei is pointless
			            tmp=term/expterm/r
			            derx=derx-tmp*rx !Calculating gradient doesn't cost detectable time, so always calculate it
			            dery=dery-tmp*ry
			            derz=derz-tmp*rz
			            if (idohess==1) then
				            tmp1=1/r2_1d5/expterm
				            tmp2=1/r2/(expterm*expterm)
				            dxx=dxx+term*(tmp1*rx2-1/r/expterm+tmp2*rx2)
				            dyy=dyy+term*(tmp1*ry2-1/r/expterm+tmp2*ry2)
				            dzz=dzz+term*(tmp1*rz2-1/r/expterm+tmp2*rz2)
				            tmp=term*(tmp1+tmp2)
				            dxy=dxy+rx*ry*tmp
				            dyz=dyz+ry*rz*tmp
				            dxz=dxz+rx*rz*tmp
			            end if
		            end do
	            else !Heavier than Ar
                    if (r>atmrhocut(iele)) cycle
		            call genatmraddens(iele,rhoarr,npt) !Extract spherically averaged radial density of corresponding element at specific grids
		            if (idohess==0) then
						call lagintpol(atmradpos(1:npt),rhoarr(1:npt),npt,r,term,der1r,der2r,2)
		            else if (idohess==1) then
						call lagintpol(atmradpos(1:npt),rhoarr(1:npt),npt,r,term,der1r,der2r,3)
                    end if
		            elerho=elerho+term
		            der1rdr=der1r/r
		            derx=derx+der1rdr*rx
		            dery=dery+der1rdr*ry
		            derz=derz+der1rdr*rz
		            if (idohess==1) then
			            tmpval=(der2r-der1rdr)/r2
			            dxx=dxx+der1rdr+tmpval*rx2
			            dyy=dyy+der1rdr+tmpval*ry2
			            dzz=dzz+der1rdr+tmpval*rz2
			            dxy=dxy+tmpval*rx*ry
			            dyz=dyz+tmpval*ry*rz
			            dxz=dxz+tmpval*rx*rz
		            end if
	            end if
            end do
        end do
    end do
end do
if (present(elegrad)) then
	elegrad(1)=derx
	elegrad(2)=dery
	elegrad(3)=derz
end if
if (idohess==1) then
	elehess(1,1)=dxx
	elehess(2,2)=dyy
	elehess(3,3)=dzz
	elehess(1,2)=dxy
	elehess(2,3)=dyz
	elehess(1,3)=dxz
	elehess(2,1)=dxy
	elehess(3,2)=dyz
	elehess(3,1)=dxz
end if
end subroutine




!!------ Return radial electron density of an element at r
!itype defines how the atomic densities will be evaluated
!itype=-2: Fitted by no more than 10 GTFs
!itype=-1: Fitted by a few STOs
!itype=0: Lagrangian interpolation based on built-in atomic radial density
!itype=18: Use STO fitted atomic density for element <=18, quality is quite poor, not normalized to expected electron number, and thus highly deprecated! (From SI of RDG original paper)
!  Accuracy: 0>-2>-1>>18   Cost: -2>=0>-1=18
!  0 is best default choice because most accurate and not expensive. However its definition is truncated at finite distance, so if density at more distant region is needed, use -1(accurate) or -2 (cheaper)
real*8 function eleraddens(iele,r,itype)
integer iele,itype
real*8 r,rhoarr(200),atomcoeff(10),atomexp(10)

eleraddens=0
if (iele==0) return !Bq atom
if (itype==0) then !Interpolation by Lagrangian interpolation based on built-in grid data
    if (r>atmrhocut(iele)) return !r is longer than the maximum distance that rho is prebuilt, unable to calculate
	call genatmraddens(iele,rhoarr,npt) !Extract spherically averaged radial density of corresponding element at specific grids
	call lagintpol(atmradpos(1:npt),rhoarr(1:npt),npt,r,eleraddens,der1r,der2r,1)
else if (itype==-1) then !STOs fitted by Tian Lu
	call genatmraddens_STOfitparm(iele,nSTO,atomcoeff,atomexp)
    do iSTO=1,nSTO
		eleraddens=eleraddens+atomcoeff(iSTO)*exp(-r*atomexp(iSTO))
	end do
else if (itype==-2) then !GTFs fitted by Tian Lu
	call genatmraddens_GTFfitparm(iele,nGTF,atomcoeff,atomexp)
    do iGTF=1,nGTF
		eleraddens=eleraddens+atomcoeff(iGTF)*exp(-r**2*atomexp(iGTF))
	end do
else if (itype==18) then !H~Ar, STO fitted by YWT group
	if (iele>18) return !Cannot evaluate
	do iSTO=1,3
		if (YWTatomcoeff(iele,iSTO)==0D0) cycle
		eleraddens=eleraddens+YWTatomcoeff(iele,iSTO)*exp(-r/YWTatomexp(iele,iSTO))
	end do
end if
end function


!!------ Calculate isolated atomic density. PBC is supported. Meaning of itype is identical to "eleraddens"
real*8 function calcatmdens(iatm,x,y,z,itype)
real*8 x,y,z,tvec(3)
integer iatm,itype
iele=a(iatm)%index
calcatmdens=0
if (ifPBC==0) then
    r=dsqrt( (a(iatm)%x-x)**2 + (a(iatm)%y-y)**2 + (a(iatm)%z-z)**2 )
    calcatmdens=calcatmdens+eleraddens(iele,r,itype)
else !Periodic case
    call getpointcell(x,y,z,ic,jc,kc)
    do icell=ic-PBCnx,ic+PBCnx
        do jcell=jc-PBCny,jc+PBCny
            do kcell=kc-PBCnz,kc+PBCnz
                call tvec_PBC(icell,jcell,kcell,tvec)
                atmx=a(iatm)%x+tvec(1)
                atmy=a(iatm)%y+tvec(2)
                atmz=a(iatm)%z+tvec(3)
                r=dsqrt( (atmx-x)**2 + (atmy-y)**2 + (atmz-z)**2 )
				calcatmdens=calcatmdens+eleraddens(iele,r,itype)
            end do
        end do
    end do
end if
end function


!!------ Calculate promolecular density. PBC is supported. Meaning of itype is identical to "eleraddens"
real*8 function calcprodens(x,y,z,itype)
real*8 x,y,z,r,tvec(3)
integer itype
calcprodens=0
if (ifPBC==0) then
	do i=1,nfragatm
		iatm=fragatm(i) !Global variable
		iele=a(iatm)%index
		r=dsqrt( (a(iatm)%x-x)**2 + (a(iatm)%y-y)**2 + (a(iatm)%z-z)**2 )
		calcprodens=calcprodens+eleraddens(iele,r,itype)
	end do
else
    call getpointcell(x,y,z,ic,jc,kc)
    do icell=ic-PBCnx,ic+PBCnx
        do jcell=jc-PBCny,jc+PBCny
            do kcell=kc-PBCnz,kc+PBCnz
                call tvec_PBC(icell,jcell,kcell,tvec)
				do i=1,nfragatm
					iatm=fragatm(i) !Global variable
					iele=a(iatm)%index
					atmx=a(iatm)%x+tvec(1)
					atmy=a(iatm)%y+tvec(2)
					atmz=a(iatm)%z+tvec(3)
					r=dsqrt( (atmx-x)**2 + (atmy-y)**2 + (atmz-z)**2 )
					calcprodens=calcprodens+eleraddens(iele,r,itype)
                end do
            end do
        end do
    end do
end if
end function




!!------ Calculate delta-g function based on promolecular approximation
real*8 function delta_g_promol(x,y,z)
real*8 x,y,z,grad(3),IGM_gradnorm
call IGMgrad_promol(x,y,z,fragatm,grad,IGM_gradnorm)
delta_g_promol=IGM_gradnorm-dsqrt(sum(grad**2))
end function

!!----- Calculate gradient or IGM gradient for a fragment based on promolecular density
!grad(1:3): Returned usual gradient vector; IGM_grad: Returned IGM type of gradient norm
!Only the atoms in "atmlist" will be taken into account
subroutine IGMgrad_promol(x,y,z,atmlist,grad,IGM_gradnorm)
real*8 x,y,z,atmgrad(3),grad(3),IGM_gradnorm
integer atmlist(:)
grad=0D0
IGM_gradnorm=0D0
do i=1,size(atmlist)
    iatm=atmlist(i)
    call proatmgrad(2,iatm,x,y,z,atmrho,atmgrad)
	grad=grad+atmgrad
	IGM_gradnorm=IGM_gradnorm+dsqrt(sum(atmgrad**2))
end do
end subroutine



!!------ Calculate density and gradient of atom "iatm" in free state fully using built-in density
!itype=1: Lagrangian interpolation density is used
!itype=2: Tian Lu's STO fitted density is used. Accuracy is slightly poorer than 1 but much faster
subroutine proatmgrad(itype,iatm,x,y,z,rho,grad)
real*8 rhoarr(200),rho,grad(3),tvec(3),atomcoeff(10),atomexp(10)
integer itype

iele=a(iatm)%index
rho=0
grad=0
!Two codes for efficiency consideration
if (ifPBC==0) then !Isolated case
	rx=x-a(iatm)%x
	ry=y-a(iatm)%y
	rz=z-a(iatm)%z
	rx2=rx*rx
	ry2=ry*ry
	rz2=rz*rz
	r2=rx2+ry2+rz2
	if (r2<atmrhocutsqr(iele)) then
		r=dsqrt(r2)
        if (itype==1) then
			call genatmraddens(iele,rhoarr,npt) !Extract spherically averaged radial density of corresponding element at specific grids
			call lagintpol(atmradpos(1:npt),rhoarr(1:npt),npt,r,rhotmp,der1r,der2r,2)
			rho=rho+rhotmp
            if (r/=0) then
				der1rdr=der1r/r
				grad(1)=grad(1)+der1rdr*rx
				grad(2)=grad(2)+der1rdr*ry
				grad(3)=grad(3)+der1rdr*rz
            end if
        else
            call genatmraddens_STOfitparm(iele,nSTO,atomcoeff,atomexp)
			do iSTO=1,nSTO
				term=atomcoeff(iSTO)*dexp(-r*atomexp(iSTO))
				rho=rho+term
				if (r/=0) then
					tmp=term*atomexp(iSTO)/r
					grad(1)=grad(1)-tmp*rx
					grad(2)=grad(2)-tmp*ry
					grad(3)=grad(3)-tmp*rz
                end if
            end do
        end if
	end if
else !Periodic case
	call getpointcell(x,y,z,ic,jc,kc)
	do icell=ic-PBCnx,ic+PBCnx
		do jcell=jc-PBCny,jc+PBCny
			do kcell=kc-PBCnz,kc+PBCnz
				call tvec_PBC(icell,jcell,kcell,tvec)
				!rx=a(iatm)%x+tvec(1)-x !Wrong code, older than 2022-Sep-18
				!ry=a(iatm)%y+tvec(2)-y
				!rz=a(iatm)%z+tvec(3)-z
				rx=x-tvec(1)-a(iatm)%x
				ry=y-tvec(2)-a(iatm)%y
				rz=z-tvec(3)-a(iatm)%z
				rx2=rx*rx
				ry2=ry*ry
				rz2=rz*rz
				r2=rx2+ry2+rz2
				if (r2>atmrhocutsqr(iele)) cycle
				r=dsqrt(r2)
				if (itype==1) then
					call genatmraddens(iele,rhoarr,npt) !Extract spherically averaged radial density of corresponding element at specific grids
					call lagintpol(atmradpos(1:npt),rhoarr(1:npt),npt,r,rhotmp,der1r,der2r,2)
					rho=rho+rhotmp
					if (r/=0) then
						der1rdr=der1r/r
						grad(1)=grad(1)+der1rdr*rx
						grad(2)=grad(2)+der1rdr*ry
						grad(3)=grad(3)+der1rdr*rz
					end if
                else
					call genatmraddens_STOfitparm(iele,nSTO,atomcoeff,atomexp)
					do iSTO=1,nSTO
						term=atomcoeff(iSTO)*dexp(-r*atomexp(iSTO))
						rho=rho+term
						if (r/=0) then
							tmp=term*atomexp(iSTO)/r
							grad(1)=grad(1)-tmp*rx
							grad(2)=grad(2)-tmp*ry
							grad(3)=grad(3)-tmp*rz
						end if
					end do
                end if
			end do
		end do
	end do
end if

end subroutine




!!------ Calculate delta-g function defined in IGM based on Hirshfeld partition (IGMH)
!Note that by default fragatm contains the whole system
real*8 function delta_g_Hirsh(x,y,z)
real*8 x,y,z,grad(3),IGM_gradnorm
call IGMgrad_Hirsh(x,y,z,fragatm,grad,IGM_gradnorm)
delta_g_Hirsh=IGM_gradnorm-dsqrt(sum(grad**2))
end function

!!------ Calculate delta-g_inter function defined in IGM based on Hirshfeld partition (IGMH)
real*8 function delta_g_inter_Hirsh(x,y,z)
real*8 x,y,z,grad1(3),IGM_gradnorm1,grad2(3),IGM_gradnorm2,grad(3),IGM_gradnorm
call IGMgrad_Hirsh(x,y,z,frag1,grad1,IGM_gradnorm1)
call IGMgrad_Hirsh(x,y,z,frag2,grad2,IGM_gradnorm2)
grad=grad1+grad2
IGM_gradnorm=dsqrt(sum(grad1**2))+dsqrt(sum(grad2**2))
delta_g_inter_Hirsh=IGM_gradnorm-dsqrt(sum(grad**2))
end function




!!----- Calculate gradient and IGM gradient of a fragment using IGMH definition, namely using Hirshfeld partition of real density
!The fragment is defined by "atmlist" array
!realrhoin and realgradin are real density and gradient of current system. If any of them is not passed in, they will be directly calculated. &
! This design is used to avoid recalculate these expensive data if they are already available
subroutine IGMgrad_Hirsh(x,y,z,atmlist,grad,IGM_gradnorm,realrhoin,realgradin)
implicit real*8 (a-h,o-z)
real*8 x,y,z
real*8 grad(3),IGM_gradnorm
integer atmlist(:)
real*8 realrho,realgrad(3),hess(3,3)
real*8 prorho,prograd(3)
real*8 atmprorho(ncenter),atmprograd(3,ncenter)
real*8 atmgrad(3),Hirshwei(ncenter)
real*8 rhoarr(200)
real*8,optional :: realrhoin,realgradin(3)

!Obtain real density and gradient of the system
if (present(realrhoin).and.present(realgradin)) then
    realrho=realrhoin
    realgrad=realgradin
else
    call calchessmat_dens(1,x,y,z,realrho,realgrad,hess)
end if

!Obtain density and gradient of all atoms in their isolated states
do iatm=1,ncenter
    call proatmgrad(1,iatm,x,y,z,atmprorho(iatm),atmprograd(:,iatm))
end do

!Calculate promolecular density and gradient of current system
prorho=sum(atmprorho(:))
do idir=1,3
    prograd(idir)=sum(atmprograd(idir,:))
end do

!Calculate atomic Hirshfeld weights
do iatm=1,ncenter
    if (prorho==0) then
        Hirshwei(iatm)=0
    else
        Hirshwei(iatm)=atmprorho(iatm)/prorho
    end if
end do

!Calculate density, usual and IGM type of density gradient of the fragment
rho=0D0
grad=0D0
IGM_gradnorm=0D0
do i=1,size(atmlist)
	iatm=atmlist(i)
    
    !Atomic density partitioned by Hirshfeld
    atmrho=Hirshwei(iatm)*realrho
    !Atomic gradient partitioned by Hirshfeld
    do idir=1,3
        t1=Hirshwei(iatm)*realgrad(idir)
        if (prorho==0) then
            t2=0
            t3=0
        else
            t2=realrho/prorho*atmprograd(idir,iatm)
            t3=-realrho*atmprorho(iatm)/prorho**2 * prograd(idir)
        end if
        !atmgrad(idir)=t1+t2+t3 !Mathematically correct free-state atomic gradient but poor effect for IGMH. Older than 2022-Sep-18
        atmgrad(idir)=t1-t2-t3 !IGMH-type "special free-state atomic gradient
    end do
    
    !rho=rho+atmrho !Useless
	grad=grad+atmgrad(:)
	IGM_gradnorm=IGM_gradnorm+dsqrt(sum(atmgrad(:)**2))
end do
end subroutine




!!----- Calculate gradient and IGM gradient of a fragment using mIGM definition, namely using Hirshfeld partition of promolecular density
!The fragment is defined by "atmlist" array
subroutine IGMgrad_Hirshpromol(x,y,z,atmlist,grad,IGM_gradnorm)
implicit real*8 (a-h,o-z)
real*8 x,y,z
integer atmlist(:)
real*8 grad(3),IGM_gradnorm
real*8 prorho,prograd(3)
real*8 atmprorho(ncenter),atmprograd(3,ncenter),atmgrad(3)

!Obtain density and gradient of all atoms in their isolated states
do iatm=1,ncenter
    call proatmgrad(2,iatm,x,y,z,atmprorho(iatm),atmprograd(:,iatm))
end do

!Calculate promolecular density
prorho=sum(atmprorho(:))

grad=0D0
IGM_gradnorm=0D0
if (prorho/=0) then
	!Calculate gradient of promolecular density
	do idir=1,3
		prograd(idir)=sum(atmprograd(idir,:))
	end do
    !Calculate gradient and IGM gradient of this fragment
	do idx=1,size(atmlist)
		iatm=atmlist(idx)
		atmgrad(:)=2*atmprorho(iatm)/prorho*prograd(:) - atmprograd(:,iatm)
		grad=grad+atmgrad
		IGM_gradnorm=IGM_gradnorm+dsqrt(sum(atmgrad**2))
	end do
end if
end subroutine




!!--------------- Calculate Shannon information entropy function
!itype=1 rho/N*ln(rho/N), this is normal definition
!itype=2 rho*ln(rho), this is Shannon information density, see J. Chem. Phys., 126, 191107
real*8 function infoentro(itype,x,y,z)
real*8 x,y,z,rho
integer itype
if (nelec==0D0) then
	infoentro=0D0
else
	rho=fdens(x,y,z)
	if (itype==1) rho=fdens(x,y,z)/nelec
	if (rho<=1D-100) then
		infoentro=0D0
	else
		infoentro=-rho*log(rho)
	end if
end if
end function




!!--------- Calculate gradient and Hessian of local information entropy or Shannon entropy density
!itype=1: -rho/N*ln(rho/N)  local information entropy
!itype=2: -rho*ln(rho)  Shannon entropy density
subroutine calchessmat_Shannon(itype,x,y,z,value,grad,hess)
implicit real*8 (a-h,o-z)
integer itype
real*8 x,y,z,value,grad(3),hess(3,3),rhograd(3),rhohess(3,3)

!Calculate -rho*ln(rho) and its derivatives
call calchessmat_dens(2,x,y,z,rho,rhograd,rhohess)
value=-rho*log(rho)
tmp=log(rho)+1
grad(:)=-rhograd(:)*tmp
do i=1,3
	do j=1,3
		hess(i,j)=-rhohess(i,j)*tmp-rhograd(i)*rhograd(j)/rho
    end do
end do

!Calculate -rho/N*ln(rho/N) and its derivatives based on data above
if (itype==1) then
	grad(:)=grad(:)/nelec + log(nelec)/nelec*rhograd(:)
    hess(:,:)=hess(:,:)/nelec + log(nelec)/nelec*rhohess(:,:)
end if

end subroutine




!!--------------- Calculate total ESP
real*8 function totesp(x,y,z)
real*8 x,y,z
totesp=eleesp(x,y,z)+nucesp(x,y,z)
end function


!!--------------- Calculate total ESP, but skip the nuclei marked by variable "iskipnuc"
real*8 function totespskip(x,y,z,iskip)
real*8 x,y,z
integer iskip
totespskip=0
do i=1,ncenter
	if (i==iskip) cycle
	totespskip=totespskip+a(i)%charge/dsqrt((x-a(i)%x)**2+(y-a(i)%y)**2+(z-a(i)%z)**2)
end do
totespskip=totespskip+eleesp(x,y,z)
end function


!!---------------- Calculate ESP from nuclear or atomic charges
!At nuclear positions, this function returns 1000 instead of infinity to avoid numerical problems
real*8 function nucesp(x,y,z)
nucesp=0D0
do i=1,nfragatm
	dist2mpx=(x-a(fragatm(i))%x)**2
	dist2mpy=(y-a(fragatm(i))%y)**2
	dist2mpz=(z-a(fragatm(i))%z)**2
	dist2=dist2mpx+dist2mpy+dist2mpz
	if (dist2==0D0) then
		 nucesp=1D3
		 return
	end if
	nucesp=nucesp+a(fragatm(i))%charge/dsqrt(dist2)
end do
end function


!!--------------- Calculate ESP due to electrons
!Note that the negative sign of electron has already been taken into account
!If user requests to use the old ESP code, or LIBRETA has not been (e.g. forgotten to be) initialized for present wavefunction, then use the old one
real*8 function eleesp(Cx,Cy,Cz)
real*8 Cx,Cy,Cz
if (iESPcode==1.or.if_initlibreta==0) then
    eleesp=eleesp1(Cx,Cy,Cz)
else
	if (iESPcode==2) then
	    eleesp=eleesp2(Cx,Cy,Cz)
    else if (iESPcode==3) then
		eleesp=eleesp2_slow(Cx,Cy,Cz)
    end if
end if
end function

!!------------ Slow code of calculating ESP from electrons
real*8 function eleesp1(Cx,Cy,Cz)
implicit none
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp
real*8 Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fn(0:10) !Enough for h-type GTF, 5+5=10
real*8 twoepsqPC,tl,tm,tn,espprivate
integer nu,imo,iprim,jprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi

eleesp1=0D0
!$OMP parallel do private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi, &
!$OMP nu,imo,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alri,Amsj,Antk,Fn,&
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,&
!$OMP term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval,espprivate) shared(eleesp1) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
    espprivate=0D0
	icen=b(iprim)%center
	Aexp=b(iprim)%exp
	Ax=a(icen)%x
	Ay=a(icen)%y
	Az=a(icen)%z
	Aix=type2ix(b(iprim)%type)
	Aiy=type2iy(b(iprim)%type)
	Aiz=type2iz(b(iprim)%type)
	sumAi=Aix+Aiy+Aiz
	do jprim=iprim,nprims
		jcen=b(jprim)%center
		Bexp=b(jprim)%exp
		Bix=type2ix(b(jprim)%type)
		Biy=type2iy(b(jprim)%type)
		Biz=type2iz(b(jprim)%type)
		Bx=a(jcen)%x
		By=a(jcen)%y
		Bz=a(jcen)%z
		sumBi=Bix+Biy+Biz
		ep=Aexp+Bexp
		Px=(Ax*Aexp+Bx*Bexp)/ep
		Py=(Ay*Aexp+By*Bexp)/ep
		Pz=(Az*Aexp+Bz*Bexp)/ep
		PAx=Px-Ax
		PAy=Py-Ay
		PAz=Pz-Az
		PBx=Px-Bx
		PBy=Py-By
		PBz=Pz-Bz
		sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
		PCx=Px-Cx
		PCy=Py-Cy
		PCz=Pz-Cz
		sqPC=PCx*PCx+PCy*PCy+PCz*PCz

		tmpval=-Aexp*Bexp*sqAB/ep
		prefac=2*pi/ep*dexp(tmpval)
		
		expngaPC=dexp(-ep*sqPC)
		maxFn=sumAi+sumBi
		Fn(maxFn)=Fmch(maxFn,ep*sqPC,expngaPC)
		nu=maxFn
		twoepsqPC=2*ep*sqPC
		do while (nu>0)
			Fn(nu-1)=(expngaPC+twoepsqPC*Fn(nu))/(2*nu-1) !Cook book p280
			nu=nu-1
		end do

		tmpnuml=0
		do l=0,Aix+Bix
			tl=1D0
			if (mod(l,2)==1) tl=-1D0
			fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
			do r=0,l/2D0
				do i=0,(l-2*r)/2D0
					tmpnuml=tmpnuml+1
					Alri(tmpnuml)=getAfac(l,r,i,PCx,ep,fjtmp)
					maplri(tmpnuml)=l-2*r-i
				end do
			end do
		end do

		tmpnumm=0
		do m=0,Aiy+Biy
			tm=1D0
			if (mod(m,2)==1) tm=-1D0
			fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
			do s=0,m/2D0
				do j=0,(m-2*s)/2D0
					tmpnumm=tmpnumm+1
					Amsj(tmpnumm)=getAfac(m,s,j,PCy,ep,fjtmp)
					mapmsj(tmpnumm)=m-2*s-j
				end do
			end do
		end do

		tmpnumn=0
		do n=0,Aiz+Biz
			tn=1D0
			if (mod(n,2)==1) tn=-1D0
			fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
			do t=0,n/2D0
				do k=0,(n-2*t)/2D0
					tmpnumn=tmpnumn+1
					Antk(tmpnumn)=getAfac(n,t,k,PCz,ep,fjtmp)
					mapntk(tmpnumn)=n-2*t-k
				end do
			end do
		end do

		term=0D0
		!Now calc "term"=<psi(iprim)|1/r_Z|psi(jprim)>
		do l=1,tmpnuml
			do m=1,tmpnumm
				do n=1,tmpnumn
					term=term+Alri(l)*Amsj(m)*Antk(n)*Fn(maplri(l)+mapmsj(m)+mapntk(n))
				end do
			end do
		end do

		if (iprim/=jprim) term=2.0*term
		term=term*prefac
		addesp=0D0
		do imo=1,nmo
			addesp=addesp+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)
		end do
		espprivate=espprivate+addesp*term
	end do !end j primitive
	!$OMP critical
	eleesp1=eleesp1+espprivate
	!$OMP end critical
end do !end i primitive
!$OMP end parallel do
eleesp1=-eleesp1
end function



!!------------ Calculate ESP in a plane. In due time, cubegen will be employed instead of internal code to evaluate ESP
subroutine planeesp
use util
implicit real*8 (a-h,o-z)
character c200tmp*200,c400tmp*400,filename_tmp*200

!Check if it is possible to use cubegen to calculate ESP plane data
alive=.false.
if (cubegenpath/=" ".and.ifiletype==1) then
	inquire(file=cubegenpath,exist=alive)
	if (.not.alive.and.cubegenpath/="none") then
		write(*,"(a)") " Note: Albeit current file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been defined, &
		&the cubegen cannot be found, therefore electrostatic potential will still be calculated using internal code of Multiwfn"
	end if
end if
if (alive.and.ifiletype==1) then !Use cubegen to calculate ESP
	write(*,"(a)") " Since the input file type is fch/fchk/chk and ""cubegenpath"" parameter in settings.ini has been properly defined, &
	&now Multiwfn directly invokes cubegen to calculate electrostatic potential"
	
	!Generate cubegen input file
	open(10,file="cubegenpt.txt",status="replace")
	do ipt=1,ngridnum1
		do jpt=1,ngridnum2
            call get2Dgridxyz(ipt,jpt,rnowx,rnowy,rnowz)
			write(10,"(3f16.8)") rnowx*b2a,rnowy*b2a,rnowz*b2a
		end do
	end do
	close(10)
	
	ncubegenthreads=1 !Parallel implementation of cubegen before G16 is buggy, so test here
	if (index(cubegenpath,"G16")/=0.or.index(cubegenpath,"g16")/=0) ncubegenthreads=nthreads
	
	!if input file is .chk, convert it to .fch before invoking cubegen
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
	do ipt=1,ngridnum1
		do jpt=1,ngridnum2
			read(10,*) rnouse,rnouse,rnouse,planemat(ipt,jpt)
		end do
	end do
	close(10)
	
	!Delete intermediate files
    call delfile("cubegenpt.txt ESPresult.cub nouseout")
	return

else
	!Calculate ESP of electron contribution
    if (iESPcode==1) then !Old slow code, but optimized specifically for plane grid
        call planeeleesp
    else if (iESPcode==2.or.iESPcode==3) then !Based on libreta
        nESPthreads=nthreads
        call doinitlibreta(1)
        if (isys==1.and.nESPthreads>12) nESPthreads=12
        write(*,*)
	    ifinish=0;ishowprog=1
        ntmp=floor(ngridnum1*ngridnum1/100D0)
        !$OMP PARALLEL DO SHARED(planemat,ifinish,ishowprog) PRIVATE(ii,jj,Cx,Cy,Cz) schedule(dynamic) NUM_THREADS(nESPthreads) collapse(2)
	    do ii=0,ngridnum1-1
		    do jj=0,ngridnum2-1
			    Cx=orgx2D+ii*v1x+jj*v2x
			    Cy=orgy2D+ii*v1y+jj*v2y
			    Cz=orgz2D+ii*v1z+jj*v2z
			    planemat(ii+1,jj+1)=eleesp(Cx,Cy,Cz)
				if (ntmp/=0) then
					!$OMP CRITICAL
					ifinish=ifinish+1
					ishowprog=mod(ifinish,ntmp)
					if (ishowprog==0) call showprog(floor(100D0*ifinish/(ngridnum1*ngridnum1)),100)
        			!$OMP END CRITICAL
                end if
		    end do
	    end do
        !$OMP END PARALLEL DO
        if (ishowprog/=0) call showprog(100,100)
    end if
    
	!Combine ESP of nuclear contribution into plane map
	do ii=0,ngridnum1-1
		do jj=0,ngridnum2-1
			Cx=orgx2D+ii*v1x+jj*v2x
			Cy=orgy2D+ii*v1y+jj*v2y
			Cz=orgz2D+ii*v1z+jj*v2z
			planemat(ii+1,jj+1)=planemat(ii+1,jj+1)+nucesp(Cx,Cy,Cz)
		end do
	end do
end if
end subroutine



!!----------- Calculate ESP contributed by electron in a plane and save to planemap, using old and slow code (but optimized for calculate plane)
!We do not allocate Alrivec,Amsjvec,Antkvec dynamically since if we do this, this routine will crash in win7-64bit system.
!The reason may be that dynamical arrays are not fully compatiable with private property in OpenMP
subroutine planeeleesp
use util
implicit real*8 (a-h,o-z)
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,Cxold,Cyold,Czold,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp
real*8 Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fnmat(0:ngridnum1-1,0:ngridnum2-1),Fnvec(0:10) !Enough for h-type GTF, 5+5=10
real*8,allocatable :: Alrivec(:,:),Amsjvec(:,:),Antkvec(:,:)
real*8 twoepsqPC,tl,tm,tn,pleprivate(ngridnum1,ngridnum2) !Store plane contribution of GTFs in each thread, then sum up
integer nu,imo,iprim,jprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,planetype,numx,numy,numz

maxnumgrid=max(ngridnum1,ngridnum2)
allocate(Alrivec(narrmax,maxnumgrid),Amsjvec(narrmax,maxnumgrid),Antkvec(narrmax,maxnumgrid))

planemat=0
ifinish=0
!$OMP PARALLEL DO private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,planetype,numx,numy,numz,&
!$OMP nu,imo,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alrivec,Amsjvec,Antkvec,Alri,Amsj,Antk,Fnmat,Fnvec,&
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,&
!$OMP Cx,Cy,Cz,Cxold,Cyold,Czold,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval,pleprivate) &
!$OMP shared(planemat,ifinish) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
	pleprivate=0D0
	icen=b(iprim)%center
	Aexp=b(iprim)%exp
	Ax=a(icen)%x
	Ay=a(icen)%y
	Az=a(icen)%z
	Aix=type2ix(b(iprim)%type)
	Aiy=type2iy(b(iprim)%type)
	Aiz=type2iz(b(iprim)%type)
	sumAi=Aix+Aiy+Aiz
	do jprim=iprim,nprims
		jcen=b(jprim)%center
		Bexp=b(jprim)%exp
		Bix=type2ix(b(jprim)%type)
		Biy=type2iy(b(jprim)%type)
		Biz=type2iz(b(jprim)%type)
		Bx=a(jcen)%x
		By=a(jcen)%y
		Bz=a(jcen)%z
		ep=Aexp+Bexp
		Px=(Ax*Aexp+Bx*Bexp)/ep
		Py=(Ay*Aexp+By*Bexp)/ep
		Pz=(Az*Aexp+Bz*Bexp)/ep
		PAx=Px-Ax
		PAy=Py-Ay
		PAz=Pz-Az
		PBx=Px-Bx
		PBy=Py-By
		PBz=Pz-Bz
		sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
		tmpval=-Aexp*Bexp*sqAB/ep
		prefac=2*pi/ep*dexp(tmpval)
			
		Cxold=999.99912D0 !An arbitrary number
		Cyold=999.99912D0
		Czold=999.99912D0
		sumBi=Bix+Biy+Biz
		maxFn=sumAi+sumBi

		!! Start cycle grid point
		do ii=0,ngridnum1-1
			do jj=0,ngridnum2-1
				Cx=orgx2D+ii*v1x+jj*v2x
				Cy=orgy2D+ii*v1y+jj*v2y
				Cz=orgz2D+ii*v1z+jj*v2z
				PCx=Px-Cx
				PCy=Py-Cy
				PCz=Pz-Cz
				sqPC=PCx*PCx+PCy*PCy+PCz*PCz
				twoepsqPC=2*ep*sqPC
				term=0D0
                expngaPC=dexp(-ep*sqPC)
				Fnmat(ii,jj)=Fmch(maxFn,ep*sqPC,expngaPC)
				Fnvec(maxFn)=Fnmat(ii,jj)
				do nu=maxFn,1,-1
					Fnvec(nu-1)=(expngaPC+twoepsqPC*Fnvec(nu))/(2*nu-1) !cook book p280
				end do

				if (Cx/=Cxold) then
					tmpnuml=0
					do l=0,Aix+Bix
						tl=1D0
						if (mod(l,2)==1) tl=-1D0
						fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
						do r=0,l/2D0
							do i=0,(l-2*r)/2D0
								tmpnuml=tmpnuml+1
								Alri(tmpnuml)=getAfac(l,r,i,PCx,ep,fjtmp)
								maplri(tmpnuml)=l-2*r-i
							end do
						end do
					end do
					Cxold=Cx
				end if
				if (Cy/=Cyold) then
					tmpnumm=0
					do m=0,Aiy+Biy
						tm=1D0
						if (mod(m,2)==1) tm=-1D0
						fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
						do s=0,m/2D0
							do j=0,(m-2*s)/2D0
								tmpnumm=tmpnumm+1
								Amsj(tmpnumm)=getAfac(m,s,j,PCy,ep,fjtmp)
								mapmsj(tmpnumm)=m-2*s-j
							end do
						end do
					end do
					Cyold=Cy
				end if
				if (Cz/=Czold) then
					tmpnumn=0
					do n=0,Aiz+Biz
						tn=1D0
						if (mod(n,2)==1) tn=-1D0
						fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
						do t=0,n/2D0
							do k=0,(n-2*t)/2D0
								tmpnumn=tmpnumn+1
								Antk(tmpnumn)=getAfac(n,t,k,PCz,ep,fjtmp)
								mapntk(tmpnumn)=n-2*t-k
							end do
						end do
					end do
					Czold=Cz
				end if

				!Now calc "term"=<psi(iprim)|1/r_Z|psi(jprim)>
				do l=1,tmpnuml
					do m=1,tmpnumm
						do n=1,tmpnumn
							term=term+Alri(l)*Amsj(m)*Antk(n)*Fnvec(maplri(l)+mapmsj(m)+mapntk(n))
						end do
					end do
				end do

				if (iprim/=jprim) term=2.0*term
				term=term*prefac
				addesp=0D0
				do imo=1,nmo
					addesp=addesp+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)
				end do
				pleprivate(ii+1,jj+1)=pleprivate(ii+1,jj+1)-addesp*term
			end do !end jj cycle
		end do !end ii cycle

	end do !end j primitive
	ifinish=ifinish+1
	call showprog(ifinish,nprims)
	!$OMP CRITICAL
	planemat=planemat+pleprivate
	!$OMP END CRITICAL
end do !end i primitive
!$OMP END PARALLEL DO
end subroutine




!!----------- Calculate grid data of total ESP
!This subroutine uses old and slow code
subroutine cubesp
use util
implicit none
integer,parameter :: narrmax=396 !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
real*8 Cx,Cy,Cz,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval
real*8 sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp
real*8 Alrivec(narrmax,nx),Amsjvec(narrmax,ny),Antkvec(narrmax,nz),Fnmat(0:nx-1,0:ny-1,0:nz-1),Fnvec(0:10) !Enough for h-type GTF, 5+5=10
real*8 twoepsqPC,tl,tm,tn,time_begin,time_end,cubprivate(nx,ny,nz)
integer nu,imo,iprim,jprim,maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn
integer Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,kk,ifinish

ifinish=0
call showprog(0,nprims)
!$OMP PARALLEL DO private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi,ii,jj,kk, &
!$OMP nu,imo,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn, &
!$OMP twoepsqPC,tl,tm,tn,Alrivec,Amsjvec,Antkvec,Fnmat,Fnvec, &
!$OMP sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,addesp,cubprivate, &
!$OMP Cx,Cy,Cz,term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,tmpval,prefac) &
!$OMP shared(cubmat,ifinish) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
    cubprivate=0D0
	icen=b(iprim)%center
	Aexp=b(iprim)%exp
	Ax=a(icen)%x
	Ay=a(icen)%y
	Az=a(icen)%z
	Aix=type2ix(b(iprim)%type)
	Aiy=type2iy(b(iprim)%type)
	Aiz=type2iz(b(iprim)%type)
	sumAi=Aix+Aiy+Aiz
	do jprim=iprim,nprims
		jcen=b(jprim)%center
		Bexp=b(jprim)%exp
		Bix=type2ix(b(jprim)%type)
		Biy=type2iy(b(jprim)%type)
		Biz=type2iz(b(jprim)%type)
		Bx=a(jcen)%x
		By=a(jcen)%y
		Bz=a(jcen)%z
		ep=Aexp+Bexp
		Px=(Ax*Aexp+Bx*Bexp)/ep
		Py=(Ay*Aexp+By*Bexp)/ep
		Pz=(Az*Aexp+Bz*Bexp)/ep
		PAx=Px-Ax
		PAy=Py-Ay
		PAz=Pz-Az
		PBx=Px-Bx
		PBy=Py-By
		PBz=Pz-Bz
		sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
		tmpval=-Aexp*Bexp*sqAB/ep
		
		prefac=2*pi/ep*dexp(tmpval)
		sumBi=Bix+Biy+Biz
		maxFn=sumAi+sumBi

		do ii=1,nx
			Cx=orgx+(ii-1)*dx
			PCx=Px-Cx
			tmpnuml=0
			do l=0,Aix+Bix
				tl=1D0
				if (mod(l,2)==1) tl=-1D0
				fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
				do r=0,l/2D0
					do i=0,(l-2*r)/2D0
						tmpnuml=tmpnuml+1
						Alrivec(tmpnuml,ii)=getAfac(l,r,i,PCx,ep,fjtmp)
						maplri(tmpnuml)=l-2*r-i
					end do
				end do
			end do
		end do
		do ii=1,ny
			Cy=orgy+(ii-1)*dy
			PCy=Py-Cy
			tmpnumm=0
			do m=0,Aiy+Biy
				tm=1D0
				if (mod(m,2)==1) tm=-1D0
				fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
				do s=0,m/2D0
					do j=0,(m-2*s)/2D0
						tmpnumm=tmpnumm+1
						Amsjvec(tmpnumm,ii)=getAfac(m,s,j,PCy,ep,fjtmp)
						mapmsj(tmpnumm)=m-2*s-j
					end do
				end do
			end do
		end do
		do ii=1,nz
			Cz=orgz+(ii-1)*dz
			PCz=Pz-Cz
			tmpnumn=0
			do n=0,Aiz+Biz
				tn=1D0
				if (mod(n,2)==1) tn=-1D0
				fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
				do t=0,n/2D0
					do k=0,(n-2*t)/2D0
						tmpnumn=tmpnumn+1
						Antkvec(tmpnumn,ii)=getAfac(n,t,k,PCz,ep,fjtmp)
						mapntk(tmpnumn)=n-2*t-k
					end do
				end do
			end do
		end do

		!! Start cycle grid point
		do kk=0,nz-1
			do jj=0,ny-1
				do ii=0,nx-1
					Cx=orgx+ii*dx
					Cy=orgy+jj*dy
					Cz=orgz+kk*dz
					PCx=Px-Cx
					PCy=Py-Cy
					PCz=Pz-Cz
					sqPC=PCx*PCx+PCy*PCy+PCz*PCz
					twoepsqPC=2*ep*sqPC
					term=0D0
                    expngaPC=dexp(-ep*sqPC)
					Fnmat(ii,jj,kk)=Fmch(maxFn,ep*sqPC,expngaPC)
                    Fnvec(maxFn)=Fnmat(ii,jj,kk)
					do nu=maxFn,1,-1
						Fnvec(nu-1)=(expngaPC+twoepsqPC*Fnvec(nu))/(2*nu-1) !Cook book p280
					end do
					do l=1,tmpnuml
						do m=1,tmpnumm
							do n=1,tmpnumn
								term=term+Alrivec(l,ii+1)*Amsjvec(m,jj+1)*Antkvec(n,kk+1)*Fnvec(maplri(l)+mapmsj(m)+mapntk(n))
							end do
						end do
					end do

					if (iprim/=jprim) term=2D0*term
					term=term*prefac
					addesp=0D0
					do imo=1,nmo
						addesp=addesp+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)
					end do
					cubprivate(ii+1,jj+1,kk+1)=cubprivate(ii+1,jj+1,kk+1)-addesp*term
				end do !end ii cycle 
			end do !end jj cycle
		end do !end kk cycle
		
	end do !end j primitive
	ifinish=ifinish+1
	call showprog(ifinish,nprims)
    !$OMP CRITICAL
	cubmat=cubmat+cubprivate
    !$OMP END CRITICAL
end do !end i primitive
!$OMP END PARALLEL DO
if (ifinish<nprims) call showprog(nprims,nprims)

!Combine nuclear contribution into it
!$OMP PARALLEL DO SHARED(cubmat) PRIVATE(i,j,k,Cx,Cy,Cz) schedule(dynamic) NUM_THREADS(nthreads)
do k=1,nz
	do j=1,ny
		do i=1,nx
            call getgridxyz(i,j,k,Cx,Cy,Cz)
			cubmat(i,j,k)=cubmat(i,j,k)+nucesp(Cx,Cy,Cz)
		end do
	end do
end do
!$OMP END PARALLEL DO
end subroutine




!!---------- Generate nuclear attraction potential integral matrix between all GTFs at a given point
!Used by PAEM
subroutine genGTFattmat(x,y,z,GTFattmat)
integer,parameter :: narrmax=396
real*8 x,y,z,GTFattmat(nprims,nprims),Alri(narrmax),Amsj(narrmax),Antk(narrmax),Fn(0:10) !Enough for h-type GTF, 5+5=10. Alri(0:10,0:5,0:5)-->11*6*6=396
integer maxFn,maplri(narrmax),mapmsj(narrmax),mapntk(narrmax),tmpnuml,tmpnumm,tmpnumn,Aix,Aiy,Aiz,Bix,Biy,Biz,r,s,t,sumAi,sumBi
!$OMP parallel do private(Aix,Aiy,Aiz,Bix,Biy,Biz,l,r,i,m,s,j,n,t,k,icen,jcen,sumAi,sumBi, &
!$OMP nu,iprim,jprim,maxFn,maplri,mapmsj,mapntk,tmpnuml,tmpnumm,tmpnumn,&
!$OMP twoepsqPC,tl,tm,tn,Alri,Amsj,Antk,Fn,sqPC,sqAB,expngaPC,PAx,PAy,PAz,PBx,PBy,PBz,PCx,PCy,PCz,fjtmp,&
!$OMP term,ep,Ax,Ay,Az,Bx,By,Bz,Aexp,Bexp,Px,Py,Pz,prefac,tmpval) shared(GTFattmat) schedule(dynamic) NUM_THREADS(nthreads)
do iprim=1,nprims
	icen=b(iprim)%center
	Aexp=b(iprim)%exp
	Ax=a(icen)%x
	Ay=a(icen)%y
	Az=a(icen)%z
	Aix=type2ix(b(iprim)%type)
	Aiy=type2iy(b(iprim)%type)
	Aiz=type2iz(b(iprim)%type)
	sumAi=Aix+Aiy+Aiz
	do jprim=iprim,nprims
		jcen=b(jprim)%center
		Bexp=b(jprim)%exp
		Bix=type2ix(b(jprim)%type)
		Biy=type2iy(b(jprim)%type)
		Biz=type2iz(b(jprim)%type)
		Bx=a(jcen)%x
		By=a(jcen)%y
		Bz=a(jcen)%z
		sumBi=Bix+Biy+Biz
		ep=Aexp+Bexp
		Px=(Ax*Aexp+Bx*Bexp)/ep
		Py=(Ay*Aexp+By*Bexp)/ep
		Pz=(Az*Aexp+Bz*Bexp)/ep
		PAx=Px-Ax
		PAy=Py-Ay
		PAz=Pz-Az
		PBx=Px-Bx
		PBy=Py-By
		PBz=Pz-Bz
		sqAB=(Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2
		PCx=Px-x
		PCy=Py-y
		PCz=Pz-z
		sqPC=PCx*PCx+PCy*PCy+PCz*PCz
		tmpval=-Aexp*Bexp*sqAB/ep
		prefac=2*pi/ep*dexp(tmpval)
        expngaPC=dexp(-ep*sqPC)
		maxFn=sumAi+sumBi
		Fn(maxFn)=Fmch(maxFn,ep*sqPC,expngaPC)
		nu=maxFn
		twoepsqPC=2*ep*sqPC
		do while (nu>0)
			Fn(nu-1)=(expngaPC+twoepsqPC*Fn(nu))/(2*nu-1) !cook book p280
			nu=nu-1
		end do
		
		tmpnuml=0
		do l=0,Aix+Bix
			tl=1D0
			if (mod(l,2)==1) tl=-1D0
			fjtmp=fj(l,Aix,Bix,PAx,PBx)*tl*fact(l)
			do r=0,l/2D0
				do i=0,(l-2*r)/2D0
					tmpnuml=tmpnuml+1
					Alri(tmpnuml)=getAfac(l,r,i,PCx,ep,fjtmp)
					maplri(tmpnuml)=l-2*r-i
				end do
			end do
		end do
		tmpnumm=0
		do m=0,Aiy+Biy
			tm=1D0
			if (mod(m,2)==1) tm=-1D0
			fjtmp=fj(m,Aiy,Biy,PAy,PBy)*tm*fact(m)
			do s=0,m/2D0
				do j=0,(m-2*s)/2D0
					tmpnumm=tmpnumm+1
					Amsj(tmpnumm)=getAfac(m,s,j,PCy,ep,fjtmp)
					mapmsj(tmpnumm)=m-2*s-j
				end do
			end do
		end do
		tmpnumn=0
		do n=0,Aiz+Biz
			tn=1D0
			if (mod(n,2)==1) tn=-1D0
			fjtmp=fj(n,Aiz,Biz,PAz,PBz)*tn*fact(n)
			do t=0,n/2D0
				do k=0,(n-2*t)/2D0
					tmpnumn=tmpnumn+1
					Antk(tmpnumn)=getAfac(n,t,k,PCz,ep,fjtmp)
					mapntk(tmpnumn)=n-2*t-k
				end do
			end do
		end do

		term=0D0
		!Now calc "term"=<psi(iprim)|1/r12|psi(jprim)>, r1 is inputted x,y,z, r2 is the integration variable
		do l=1,tmpnuml
			do m=1,tmpnumm
				do n=1,tmpnumn
					term=term+Alri(l)*Amsj(m)*Antk(n)*Fn(maplri(l)+mapmsj(m)+mapntk(n))
				end do
			end do
		end do
		term=term*prefac
		GTFattmat(iprim,jprim)=term
		GTFattmat(jprim,iprim)=term
	end do !end j primitive
end do !end i primitive
!$OMP end parallel do
end subroutine


!!-------------------- Calculate A-factor at Cook book p245
real*8 function getAfac(l,r,i,PC,gamma,fjtmp)
integer l,r,i,ti
real*8 gamma,PC,comp,PCterm,fjtmp
ti=1D0
if (mod(i,2)==1) ti=-1D0 !faster than ti=(-1)**i
PCterm=1D0
if (l-2*r-2*i/=0) PCterm=PC**(l-2*r-2*i)
comp=ti*PCterm*(0.25D0/gamma)**(r+i) / ( fact(r)*fact(i)*fact(l-2*r-2*i) )
getAfac=fjtmp*comp
! comp=(-1)**i*fact(l)*PC**(l-2*r-2*i)*(1/(4*gamma))**(r+i) / ( fact(r)*fact(i)*fact(l-2*r-2*i) )
! getAfac=(-1)**l * fj(l,l1,l2,PA,PB)*comp
end function      


!!---------------- Calculate fj at Cook book p237
real*8 function fj(j,l,m,aa,bb)
real*8 aa,bb,pre,at,bt
integer j,l,m,k,imin,imax
imax=min(j,l)
imin=max(0,j-m)
fj=0D0
do k=imin,imax
	pre=fact(l)/fact(l-k)/fact(k) * fact(m)/fact(m-j+k)/fact(j-k)
	at=1D0
	bt=1D0
	if (l-k/=0) at=aa**(l-k)  !This determine helps to improve efficient
	if (m+k-j/=0) bt=bb**(m+k-j)
	fj=fj+pre*at*bt
end do
end function


!!------- Calculate Boys function int('t^(2*m)*exp(-x*t^2)','t',0,1), see Cook book p281 for detail
!expnx is input parameter, value should be exp(-x), because calculation of the value is time-consuming and in
!other place this value also need be calculate, so not recalculate in this subroutine
real*8 function Fmch(m,x,expnx)
IMPLICIT none
integer m,i
real*8 x,expnx,a,b,term,partsum,APPROX,xd,FIMULT,NOTRMS,eps,fiprop
eps=1.0D-8  !Convergence precision
Fmch=0D0
if (x<=10) then
	if (expnx==0D0) RETURN
	A=m+0.5D0
	term=1D0/A
	partsum=term
	DO I=2,50
		A=A+1D0
		term=term*X/A
		partsum=partsum+term
		if ( term/partsum < eps) THEN
		   Fmch = 0.5D0*partsum*expnx
		   RETURN
		END IF
	END DO
	write(*,*) "Error: Fmch did not converge"
else !x is big, use suitable method for solve this situation
	A=M
	B=A+0.5D0
	A=A-0.5D0
	XD=1D0/X
	APPROX=0.88622692D0*(dsqrt(XD)*XD**m)
	DO I=1,m
		B=B-1D0
		APPROX=APPROX*B
	END DO
	FIMULT=0.5D0*expnx*XD
	partsum=0.D0
	IF (FIMULT==0D0) THEN
		Fmch=APPROX-FIMULT*partsum
		return
	ELSE
		FIPROP=FIMULT/APPROX
		term=1D0
		partsum=term
		NOTRMS=X
		NOTRMS=NOTRMS+M
		DO I=2,NOTRMS
		   term=term*A*XD
		   partsum=partsum+term
		   IF (dabs(term*FIPROP/partsum)<eps)  THEN
			  Fmch=APPROX-FIMULT*partsum
			  RETURN
		   END IF
		   A=A-1D0
		END DO
		write(*,*) "Error: Fmch did not converge"
	END IF
end if
end function




!!------------- Generate Becke weight function, only used by Sobereva
!If inp2=0, then return atomic space weight of atom inp1, else return overlap weight of atom inp1 and inp2
real*8 function beckewei(x,y,z,inp1,inp2)
real*8 x,y,z
integer inp1,inp2
! integer :: fraglist(13)=(/ 24,12,23,4,11,3,22,1,10,2,18,20,21 /)
real*8 Pvec(ncenter)

!Calculate Becke weight of all atoms (Pvec) at current point
call BeckePvec(x,y,z,Pvec,covr_tianlu,3)
if (inp2==0) then
	beckewei=Pvec(inp1)
else
	beckewei=Pvec(inp1)*Pvec(inp2)
end if
! beckewei=sum(Pvec(fraglist)) !Get fragment Becke weight
end function




!!-------- Calculate electron density ellipticity (itype=1), eta index (itype=2), stiffness (itype=3)
real*8 function densellip(x,y,z,itype)
use util
integer itype
real*8 x,y,z,dens,grad(3),hess(3,3),eigval(3),eigvecmat(3,3)
call calchessmat_dens(2,x,y,z,dens,grad,hess)
call diagmat(hess,eigvecmat,eigval,300,1D-12)
call sort(eigval)
if (itype==1) then
	densellip=eigval(1)/eigval(2)-1
else if (itype==2) then
	densellip=abs(eigval(1))/eigval(3)
else if (itype==3) then
	densellip=abs(eigval(2))/eigval(3)
end if
end function




!!-------------- Calculate single exponential decay detector (SEDD)
!Originally proposed in ChemPhysChem, 13, 3462 (2012), current implementation is based on the new definition in DORI paper (10.1021/ct500490b)
real*8 function SEDD(x,y,z)
real*8 x,y,z,rho,grad(3),hess(3,3)
call calchessmat_dens(2,x,y,z,rho,grad,hess)
dersqr=sum(grad**2)
tmp1_1=rho*(grad(1)*hess(1,1)+grad(2)*hess(1,2)+grad(3)*hess(1,3))
tmp1_2=grad(1)*dersqr
tmp2_1=rho*(grad(1)*hess(1,2)+grad(2)*hess(2,2)+grad(3)*hess(2,3))
tmp2_2=grad(2)*dersqr
tmp3_1=rho*(grad(1)*hess(1,3)+grad(2)*hess(2,3)+grad(3)*hess(3,3))
tmp3_2=grad(3)*dersqr
eps=4/rho**8*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )
SEDD=dlog(1+eps)
end function




!!-------------- Calculate density Overlap Regions Indicator (DORI)
real*8 function DORI(x,y,z)
real*8 x,y,z,rho,grad(3),hess(3,3)
call calchessmat_dens(2,x,y,z,rho,grad,hess)
dersqr=sum(grad**2)
tmp1_1=rho*(grad(1)*hess(1,1)+grad(2)*hess(1,2)+grad(3)*hess(1,3))
tmp1_2=grad(1)*dersqr
tmp2_1=rho*(grad(1)*hess(1,2)+grad(2)*hess(2,2)+grad(3)*hess(2,3))
tmp2_2=grad(2)*dersqr
tmp3_1=rho*(grad(1)*hess(1,3)+grad(2)*hess(2,3)+grad(3)*hess(3,3))
tmp3_2=grad(3)*dersqr
theta=4/dersqr**3*( (tmp1_1-tmp1_2)**2 + (tmp2_1-tmp2_2)**2 + (tmp3_1-tmp3_2)**2 )
DORI=theta/(1+theta)
end function




!!------------- Calculate interaction region indicator (IRI)
real*8 function IRIfunc(x,y,z)
real*8 x,y,z,rho,grad(3),hess(3,3)
call calchessmat_dens(1,x,y,z,rho,grad,hess)
gradnorm=dsqrt(sum(grad**2))
tmp=uservar
if (uservar==0) tmp=1.1D0
IRIfunc=gradnorm/rho**tmp
if (rho==0D0) IRIfunc=5 !Arbitary large value
if (IRI_rhocut/=0.and.(gradnorm==0D0.or.rho<=IRI_rhocut)) IRIfunc=5
end function




!!--------- Calculate gradient and Hessian of IRI or RDG or others with general form: c*|der_rho|/rho^a
!itype=1 Only calculate value and grad, =2 Calculate value, gradient and Hessian
!ifunc=1: IRI, =2: RDG
subroutine calchessmat_IRI_RDG(itype,ifunc,x,y,z,value,grad,hess)
implicit real*8 (a-h,o-z)
integer itype,ifunc
real*8 x,y,z,value,grad(3),hess(3,3),rhograd(3),rhohess(3,3),rhotens3(3,3,3),gnormgrad(3),gnormhess(3,3)

call rho_tensor(x,y,z,rho,rhograd,rhohess,rhotens3)
rhograd2=sum(rhograd**2)
gnorm=dsqrt(rhograd2)

if (ifunc==1) then !IRI
	cfac=1
    afac=1.1D0
else if (ifunc==2) then !RDG
	cfac=0.161620459673995D0 !=1/(2*(3*pi**2)**(1/3))
    afac=4D0/3D0
end if
value=cfac*gnorm/rho**afac

!Calculate gradient of intermediate quantity |der_rho|
gnormgrad=0
do idir=1,3
	do jdir=1,3
		gnormgrad(idir)=gnormgrad(idir)+rhograd(jdir)*rhohess(idir,jdir)
    end do
end do
gnormgrad=gnormgrad/gnorm

!Calculate gradient
tmp1=-afac*gnorm/rho**(afac+1)
do i=1,3
    grad(i)=tmp1*rhograd(i)+gnormgrad(i)/rho**afac
end do
grad=grad*cfac

if (itype==2) then
	!Calculate Hessian of intermediate quantity |der_rho|
	do i=1,3
		do j=1,3
			tmp=0
			tmp2=0
			do k=1,3
				tmp=tmp+rhograd(k)*rhohess(k,i)
				tmp2=tmp2+rhohess(k,i)*rhohess(k,j)+rhograd(k)*rhotens3(k,i,j)
            end do
			gnormhess(i,j)=-1D0/rhograd2*gnormgrad(j)*tmp+tmp2/gnorm
		end do
	end do

	!Calculate Hessian
	do i=1,3
		do j=1,3
			tmp1=-(afac+1)/rho**(afac+2)*rhograd(j)*gnorm*rhograd(i) + 1D0/rho**(afac+1)*gnormgrad(j)*rhograd(i) + gnorm/rho**(afac+1)*rhohess(i,j)
			tmp2=-afac/rho**(afac+1)*rhograd(j)*gnormgrad(i)
			tmp3=gnormhess(i,j)/rho**afac
			hess(i,j) = -afac*tmp1 + tmp2 + tmp3
		end do
	end do
	hess=hess*cfac
end if
end subroutine




!!------------- Simultaneously calculate interaction region indicator and sign(lambda2)rho to reduce cost
subroutine IRI_s2lr(x,y,z,sl2rval,IRIval)
use util
real*8 x,y,z,IRIval,sl2rval,rho,grad(3),hess(3,3)
real*8 eigvecmat(3,3),eigval(3)

call calchessmat_dens(2,x,y,z,rho,grad,hess)
gradnorm=dsqrt(sum(grad**2))
tmp=uservar
if (uservar==0) tmp=1.1D0

IRIval=gradnorm/rho**tmp
if (IRI_rhocut/=0.and.(gradnorm==0D0.or.rho<=IRI_rhocut)) IRIval=5

call diagmat(hess,eigvecmat,eigval,100,1D-10)
call sort(eigval)
if (eigval(2)==0D0) then !When eigval(2)==0, eigval(2)/abs(eigval(2)) can't be calculated, rho generally will be zero, so sign is not important
	sl2rval=rho
else
	sl2rval=rho*eigval(2)/abs(eigval(2))
end if
end subroutine




!!----- Integrand of LSDA exchange functional
real*8 function xLSDA(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
xLSDA=-3D0/2D0*(3D0/4D0/pi)**(1D0/3D0)*(adens**(4D0/3D0)+bdens**(4D0/3D0))
end function

!!----- Integrand of Becke88 exchange functional
real*8 function xBecke88(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
adens4d3=adens**(4D0/3D0)
bdens4d3=bdens**(4D0/3D0)
slatercoeff=-3D0/2D0*(3D0/4D0/pi)**(1D0/3D0)
slaterxa=slatercoeff*adens4d3
slaterxb=slatercoeff*bdens4d3
slaterx=slaterxa+slaterxb
redagrad=agrad/adens4d3 !Reduced density gradient
redbgrad=bgrad/bdens4d3
arshredagrad=log(redagrad+dsqrt(redagrad**2+1))
Beckexa=adens4d3*redagrad**2/(1+6*0.0042D0*redagrad*arshredagrad)
arshredbgrad=log(redbgrad+dsqrt(redbgrad**2+1))
Beckexb=bdens4d3*redbgrad**2/(1+6*0.0042D0*redbgrad*arshredbgrad)
Beckex=-0.0042D0*(Beckexa+Beckexb)
xBecke88=slaterx+Beckex
end function

!!------ Integrand of LYP corelation functional
real*8 function cLYP(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
parma=0.04918D0
parmb=0.132D0
parmc=0.2533D0
parmd=0.349D0
densn1d3=tdens**(-1D0/3D0)
parmw=exp(-parmc*densn1d3) / (1+parmd*densn1d3) * tdens**(-11D0/3D0)
parmdel=parmc*densn1d3+parmd*densn1d3/(1+parmd*densn1d3)
parmCf=3D0/10D0*(3*pi*pi)**(2D0/3D0)
tmp1=-parma*4D0/(1+parmd*densn1d3)*adens*bdens/tdens
tmp2a=2**(11D0/3D0)*parmCf*(adens**(8D0/3D0)+bdens**(8D0/3D0))
tmp2b=(47D0/18D0-7D0/18D0*parmdel)*tgrad**2
tmp2c=-(2.5D0-parmdel/18D0)*(agrad**2+bgrad**2)
tmp2d=-(parmdel-11D0)/9D0*(adens/tdens*agrad**2+bdens/tdens*bgrad**2)
tmp2=adens*bdens*(tmp2a+tmp2b+tmp2c+tmp2d)
tmp3=-2D0/3D0*tdens**2*tgrad**2+(2D0/3D0*tdens**2-adens**2)*bgrad**2+(2D0/3D0*tdens**2-bdens**2)*agrad**2
cLYP=tmp1-parma*parmb*parmw*(tmp2+tmp3)
end function




!!------ Integrand of Weizsacker kinetic functional (steric energy)
! DO NOT consider EDF, because the kinetic energy outputted by quantum chemistry programs is always only for explicitly represented electrons!
real*8 function weizsacker(x,y,z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),gradrho(3) !,EDFgrad(3)
call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
rho=0D0
gradrho=0D0
do i=1,nmo
	rho=rho+MOocc(i)*wfnval(i)**2
	gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do
gradrho=2*gradrho
! if (nEDFprims/=0) then
! 	call EDFrho(2,x,y,z,EDFdens,EDFgrad)
! 	rho=rho+EDFdens
! 	gradrho=gradrho+EDFgrad
! end if
if (rho<1D-30) then
	weizsacker=0
else
	weizsacker=sum(gradrho(:)**2)/8/rho
end if
end function
!!----- Weizsacker potential, essentially identical to stericpot(x,y,z)
real*8 function weizpot(x,y,z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),gradrho(3),laplx,laply,laplz,lapltot
rho=0D0
gradrho=0D0
laplx=0D0
laply=0D0
laplz=0D0
call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
do i=1,nmo
	rho=rho+MOocc(i)*wfnval(i)**2
	gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
	laplx=laplx+MOocc(i)*( wfnderv(1,i)**2 + wfnval(i)*wfnhess(1,1,i) )
	laply=laply+MOocc(i)*( wfnderv(2,i)**2 + wfnval(i)*wfnhess(2,2,i) )
	laplz=laplz+MOocc(i)*( wfnderv(3,i)**2 + wfnval(i)*wfnhess(3,3,i) )
end do
gradrho=2*gradrho
lapltot=2*(laplx+laply+laplz)
weizpot=sum(gradrho(:)**2)/8D0/rho**2-lapltot/4D0/rho
end function




!!------ Steric potential, whose negative value is one-electron potential, essentially identical to weizpot(x,y,z)
real*8 function stericpot(x,y,z)
real*8 x,y,z,gradrho(3),hessrho(3,3),lapltot
call calchessmat_dens(2,x,y,z,rho,gradrho,hessrho)
lapltot=hessrho(1,1)+hessrho(2,2)+hessrho(3,3)
if (rho<steric_potcutrho) then
	stericpot=steric_potcons
	return
end if
stericpot=sum(gradrho**2)/(rho+steric_addminimal)**2/8D0-lapltot/(rho+steric_addminimal)/4D0
end function

!!------ Calculate analytic first-order derivative of steric potential
subroutine stericderv(x,y,z,derv)
real*8 x,y,z,derv(3)
real*8 eleval,elegrad(3),elehess(3,3),laplval,laplgrad(3),rhotens3(3,3,3)
real*8 wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),wfntens3(3,3,3,nmo)
real*8 EDFval,EDFgrad(3),EDFhess(3,3),EDFtens3(3,3,3)
call orbderv(5,1,nmo,x,y,z,wfnval,wfnderv,wfnhess,wfntens3)
eleval=sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnval(1:nmo) )
elegrad(1)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(1,1:nmo) )
elegrad(2)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(2,1:nmo) )
elegrad(3)=2*sum( MOocc(1:nmo)*wfnval(1:nmo)*wfnderv(3,1:nmo) )
elehess(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
elehess(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
elehess(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
elehess(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
elehess(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
elehess(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
elehess(2,1)=elehess(1,2)
elehess(3,2)=elehess(2,3)
elehess(3,1)=elehess(1,3)
laplval=elehess(1,1)+elehess(2,2)+elehess(3,3)
rhotens3=0D0
do i=1,nmo
	rhotens3(1,1,1)=rhotens3(1,1,1)+MOocc(i)*( 3*wfnderv(1,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,1,i) )
	rhotens3(2,2,2)=rhotens3(2,2,2)+MOocc(i)*( 3*wfnderv(2,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,2,i) )
	rhotens3(3,3,3)=rhotens3(3,3,3)+MOocc(i)*( 3*wfnderv(3,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,3,i) )
	rhotens3(1,1,2)=rhotens3(1,1,2)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,2,i)+wfnderv(2,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,2,i) )
	rhotens3(1,1,3)=rhotens3(1,1,3)+MOocc(i)*( 2*wfnderv(1,i)*wfnhess(1,3,i)+wfnderv(3,i)*wfnhess(1,1,i)+wfnval(i)*wfntens3(1,1,3,i) )
	rhotens3(2,2,3)=rhotens3(2,2,3)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,3,i)+wfnderv(3,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,3,i) )
	rhotens3(1,2,2)=rhotens3(1,2,2)+MOocc(i)*( 2*wfnderv(2,i)*wfnhess(2,1,i)+wfnderv(1,i)*wfnhess(2,2,i)+wfnval(i)*wfntens3(2,2,1,i) )
	rhotens3(1,3,3)=rhotens3(1,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,1,i)+wfnderv(1,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,1,i) )
	rhotens3(2,3,3)=rhotens3(2,3,3)+MOocc(i)*( 2*wfnderv(3,i)*wfnhess(3,2,i)+wfnderv(2,i)*wfnhess(3,3,i)+wfnval(i)*wfntens3(3,3,2,i) )
end do
rhotens3=rhotens3*2D0
laplgrad(1)=rhotens3(1,1,1)+rhotens3(1,2,2)+rhotens3(1,3,3)
laplgrad(2)=rhotens3(1,1,2)+rhotens3(2,2,2)+rhotens3(2,3,3)
laplgrad(3)=rhotens3(1,1,3)+rhotens3(2,2,3)+rhotens3(3,3,3)
if (nEDFprims/=0) then
	call EDFrho(5,x,y,z,EDFval,EDFgrad,EDFhess,EDFtens3)
	eleval=eleval+EDFval
	elegrad=elegrad+EDFgrad
	elehess=elehess+EDFhess	
	laplgrad(1)=laplgrad(1)+EDFtens3(1,1,1)+EDFtens3(1,2,2)+EDFtens3(1,3,3)
	laplgrad(2)=laplgrad(2)+EDFtens3(1,1,2)+EDFtens3(2,2,2)+EDFtens3(2,3,3)
	laplgrad(3)=laplgrad(3)+EDFtens3(1,1,3)+EDFtens3(2,2,3)+EDFtens3(3,3,3)
end if
! Above codes can be simplified as below two lines, but will consume additional 1/3 time
! call calchessmat_dens(2,x,y,z,eleval,elegrad,elehess)
! call calchessmat_lapl(1,x,y,z,laplval,laplgrad,laplhess)

eleval=eleval+steric_addminimal
elenorm2=sum(elegrad**2) !Square of norm of electron density gradient
do i=1,3 !x,y,z
	tmp1=(elegrad(1)*elehess(1,i)+elegrad(2)*elehess(2,i)+elegrad(3)*elehess(3,i)) / eleval**2 -elenorm2/eleval**3*elegrad(i)
	tmp2=-laplgrad(i)/eleval+laplval/eleval**2*elegrad(i)
	derv(i)=(tmp1+tmp2)/4D0
end do
end subroutine

!!---- Magnitude of steric force
real*8 function stericforce(x,y,z)
real*8 x,y,z,derv(3)
call stericderv(x,y,z,derv)
stericforce=dsqrt(sum(derv**2))
end function



!memo: Shubin reported that steric potential/force is quite sensitive to steric_addminimal, &
!so 2016-Sep-30 I devised another solution to solve the diverse behavior of steric potential/force via Becke's damping function
!!---- Steric potential with damping function to zero
real*8 function stericpot_damp(x,y,z)
real*8 x,y,z,gradrho(3),hessrho(3,3),lapltot
call calchessmat_dens(2,x,y,z,rho,gradrho,hessrho)
lapltot=hessrho(1,1)+hessrho(2,2)+hessrho(3,3)
stericpotorg=sum(gradrho**2)/rho**2/8D0-lapltot/rho/4D0
weiwidth=2 !e.g. weiwidth=2 and steric_potcutrho=-13 means the Becke damping function of [1,0] corresponds to 1D-11~1D-15
tmps=-(dlog(rho)-steric_potcutrho)/weiwidth
if (tmps<-1) then
	consorg=1
else if (tmps>1) then
	consorg=0
else
	do iter=1,2
		tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
	end do
	consorg=0.5D0*(1-tmps) !The weight to switch to constant value steric_potcons
end if
stericpot_damp=stericpotorg*consorg+steric_potcons*(1-consorg)
end function
!!---- Steric force based on damped steric potential
real*8 function stericforce_damp(x,y,z)
real*8 x,y,z,derv(3)
diffstep=2D-5
derv(1)=(stericpot_damp(x+diffstep,y,z)-stericpot_damp(x-diffstep,y,z))/(2*diffstep)
derv(2)=(stericpot_damp(x,y+diffstep,z)-stericpot_damp(x,y-diffstep,z))/(2*diffstep)
derv(3)=(stericpot_damp(x,y,z+diffstep)-stericpot_damp(x,y,z-diffstep))/(2*diffstep)
stericforce_damp=dsqrt(sum(derv**2))
end function
!!---- Steric force directly damped to zero, rather than based on damped steric potential
real*8 function stericforce_directdamp(x,y,z)
real*8 x,y,z
weiwidth=2
tmps=-(dlog(fdens(x,y,z))-steric_potcutrho)/weiwidth
if (tmps<-1) then
	consorg=1
else if (tmps>1) then
	consorg=0
else
	do iter=1,2
		tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
	end do
	consorg=0.5D0*(1-tmps)
end if
steric_addminimalold=steric_addminimal
steric_addminimal=0
stericforce_directdamp=stericforce(x,y,z)*consorg
steric_addminimal=steric_addminimalold
end function



!!------- Steric charge, =lapl(steric potential)/(-4*pi)
!Based on analytic first-order derivative, using finite difference to obtain d2v/dx2, d2v/dy2 and d2v/dz2
real*8 function stericcharge(x,y,z)
real*8 x,y,z,derv1add(3),derv1min(3)
if (fdens(x,y,z)<steric_potcutrho) then
	stericcharge=0D0
	return
end if
diffstep=2D-5
call stericderv(x+diffstep,y,z,derv1add)
call stericderv(x-diffstep,y,z,derv1min)
derv2x=(derv1add(1)-derv1min(1))/(2*diffstep) !d2v/dx2
call stericderv(x,y+diffstep,z,derv1add)
call stericderv(x,y-diffstep,z,derv1min)
derv2y=(derv1add(2)-derv1min(2))/(2*diffstep) !d2v/dy2
call stericderv(x,y,z+diffstep,derv1add)
call stericderv(x,y,z-diffstep,derv1min)
derv2z=(derv1add(3)-derv1min(3))/(2*diffstep) !d2v/dz2
stericcharge=-(derv2x+derv2y+derv2z)/4D0/pi
end function

!!------- Steric charge directly damped to zero
real*8 function stericcharge_directdamp(x,y,z)
real*8 x,y,z
weiwidth=2
tmps=-(dlog(fdens(x,y,z))-steric_potcutrho)/weiwidth
if (tmps<-1) then
	consorg=1
else if (tmps>1) then
	consorg=0
else
	do iter=1,2
		tmps=1.5D0*(tmps)-0.5D0*(tmps)**3
	end do
	consorg=0.5D0*(1-tmps)
end if
steric_addminimalold=steric_addminimal
steric_addminimal=0
stericcharge_directdamp=stericcharge(x,y,z)*consorg
steric_addminimal=steric_addminimalold
end function



!!------ Calculate Fisher information density
!itype=1 Normal definition
!itype=2 Second Fisher information density, Eq.5 of JCP,126,191107
real*8 function Fisherinfo(itype,x,y,z)
real*8 x,y,z,eleval,elegrad(3),elehess(3,3)
integer itype
Fisherinfo=0
eleval=fdens(x,y,z)
if (eleval<=1D-30) return
if (itype==1) Fisherinfo=fgrad(x,y,z,'t')**2/eleval
if (itype==2) Fisherinfo=-flapl(x,y,z,'t')*log(eleval)
end function




!!--------- Calculate gradient and Hessian of Fisher information density
subroutine calchessmat_Fisherinfo(x,y,z,value,grad,hess)
implicit real*8 (a-h,o-z)
real*8 x,y,z,value,grad(3),hess(3,3),rhograd(3),rhohess(3,3),rhotens3(3,3,3)

call rho_tensor(x,y,z,rho,rhograd,rhohess,rhotens3)

rhograd2=sum(rhograd**2)
value=rhograd2/rho

tmp1=-rhograd2/rho**2
do i=1,3 !dx, dy, dz
	tmp2=0
	do idir=1,3 !Loop x,y,z
		tmp2=tmp2+rhograd(idir)*rhohess(idir,i)
    end do
    grad(i)=tmp1*rhograd(i)+2/rho*tmp2
end do

tmp1=2/rho**3*rhograd2
do i=1,3 !dx, dy, dz
	do j=1,3 !dx, dy, dz
		tmp2=0
		tmp3=0
		tmp4=0
		do idir=1,3 !Loop x,y,z
			tmp2=tmp2+rhograd(idir)*rhohess(idir,j)
			tmp3=tmp3+rhograd(idir)*rhohess(idir,i)
			tmp4=tmp4+rhohess(idir,i)*rhohess(idir,j)+rhograd(idir)*rhotens3(idir,i,j)
		end do
        hess(i,j) = tmp1*rhograd(i)*rhograd(j) - 2/rho**2*rhograd(i)*tmp2 - rhograd2/rho**2*rhohess(i,j) - 2/rho**2*rhograd(j)*tmp3 + 2/rho*tmp4
    end do
end do
end subroutine




!!--------- Calculate gradient (analytically) and Hessian (semi-numerical) of second Fisher information density
subroutine calchessmat_second_Fisherinfo(x,y,z,value,grad,hess)
implicit real*8 (a-h,o-z)
real*8 x,y,z,value,grad(3),hess(3,3)
real*8 gradaddx(3),gradminx(3),gradaddy(3),gradminy(3),gradaddz(3),gradminz(3)
call second_Fisherinfo_grad(x,y,z,value,grad)
diff=8D-4
denom=2D0*diff
call second_Fisherinfo_grad(x+diff,y,z,tmpval,gradaddx)
call second_Fisherinfo_grad(x-diff,y,z,tmpval,gradminx)
call second_Fisherinfo_grad(x,y+diff,z,tmpval,gradaddy)
call second_Fisherinfo_grad(x,y-diff,z,tmpval,gradminy)
call second_Fisherinfo_grad(x,y,z+diff,tmpval,gradaddz)
call second_Fisherinfo_grad(x,y,z-diff,tmpval,gradminz)
hess(1,1)=(gradaddx(1)-gradminx(1))/denom
hess(2,2)=(gradaddy(2)-gradminy(2))/denom
hess(3,3)=(gradaddz(3)-gradminz(3))/denom
hess(1,2)=(gradaddy(1)-gradminy(1))/denom
hess(2,3)=(gradaddz(2)-gradminz(2))/denom
hess(1,3)=(gradaddz(1)-gradminz(1))/denom
hess(2,1)=hess(1,2)
hess(3,2)=hess(2,3)
hess(3,1)=hess(1,3)
end subroutine
!!--------- Analytically calculate gradient of second Fisher information density
subroutine second_Fisherinfo_grad(x,y,z,value,grad)
implicit real*8 (a-h,o-z)
real*8 x,y,z,value,grad(3),rhograd(3),rhohess(3,3),rhotens3(3,3,3)
call rho_tensor(x,y,z,rho,rhograd,rhohess,rhotens3)
rholapl=rhohess(1,1)+rhohess(2,2)+rhohess(3,3)
value=-rholapl*log(rho)
do i=1,3
	grad(i)=-(rhotens3(1,1,i)+rhotens3(2,2,i)+rhotens3(3,3,i))*log(rho)-rholapl/rho*rhograd(i)
end do
end subroutine




!!------ Ghosh entropy density
!If itype==1, G(r) will be used as kinetic energy density
!If itype==2, G(r)-der2rho/8 will be used instead, which is the kinetic energy density exactly corresponding to Eq. 22 of PNAS, 81, 8028.
real*8 function Ghoshentro(x,y,z,itype)
integer itype
real*8 kintot,x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo)
if (itype==1) call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
if (itype==2) call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess) !If K(r) is used, use this
rho=0D0
do i=1,nmo
	rho=rho+MOocc(i)*wfnval(i)**2
end do
ck=2.871234D0
TFkin=ck*rho**(5D0/3D0)
kintot=0D0
!   If we use Hamiltonian kinetic density
! hamx=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
! hamy=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
! hamz=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
! kintot=-(hamx+hamy+hamz)/2
!   If we use Lagrangian kinetic density G(r)
do i=1,nmo
	kintot=kintot+MOocc(i)*sum(wfnderv(:,i)**2)
end do
kintot=kintot/2D0
if (itype==2) then
	xlapl=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
	ylapl=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
	zlapl=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
	kintot=kintot-(xlapl+ylapl+zlapl)/8
end if
rlambda=5D0/3D0+log(4D0*pi*ck/3D0)
if (kintot<0) then
	rlogterm=0
else
	rlogterm=log(kintot/TFkin)
end if
Ghoshentro=1.5D0*rho*(rlambda+rlogterm)
end function




!!------ Relative Shannon entropy density, also known as information gain density
!Must call generate_promolwfn prior to using this function to make CO_pmol, MOocc_pmol, etc. available
real*8 function relShannon(x,y,z)
real*8 x,y,z,wfnval_pmol(nmo_pmol)
!The following commented code borrow regular code to compute
!Calculate rho_0
!deallocate(MOocc,CO)
!allocate(MOocc(nmo_pmol),CO(nmo_pmol,nprims))
!nmo=nmo_pmol
!MOocc=MOocc_pmol
!CO=CO_pmol
!rho_0=fdens(x,y,z)
!!Calculate rho
!deallocate(MOocc,CO)
!nmo=nmo_org
!allocate(MOocc(nmo),CO(nmo,nprims))
!MOocc=MOocc_org
!CO=CO_org

rho=fdens(x,y,z)

call orbderv_pmol(1,1,nmo_pmol,x,y,z,wfnval_pmol)
rho_0=sum( MOocc_pmol(1:nmo_pmol)*wfnval_pmol(1:nmo_pmol)**2 )

relShannon=rho*log(rho/rho_0)
end function




!!--------- Calculate gradient and Hessian of relative Shannon entropy density (information gain)
!Must call generate_promolwfn prior to using this function to make CO_pmol, MOocc_pmol, etc. available
subroutine calchessmat_relShannon(x,y,z,value,grad,hess)
real*8 x,y,z,value,grad(3),hess(3,3)
real*8 rhograd(3),rhohess(3,3),prhograd(3),prhohess(3,3)
call calchessmat_dens(2,x,y,z,rho,rhograd,rhohess)
call calchessmat_dens_promol(2,x,y,z,prho,prhograd,prhohess)
!Calculate gradient
do i=1,3
    grad(i)=rhograd(i)*(log(rho/prho)+1) - rho/prho*prhograd(i)
end do
!Calculate Hessian
do i=1,3
    do j=1,3
        t1= (log(rho/prho)+1)*rhohess(i,j)
        t2= -rho/prho*prhohess(i,j)
        t3= 1/rho*rhograd(i)*rhograd(j)
        t4= rho/prho**2 * prhograd(i)*prhograd(j)
        t5= -1/prho*(rhograd(i)*prhograd(j)+rhograd(j)*prhograd(i))
        hess(i,j)=t1+t2+t3+t4+t5
    end do
end do
end subroutine




!!-------- Pauli potential
!Only suitable for closed-shell cases
real*8 function Paulipot(x,y,z)
real*8 x,y,z
if (ispecial==0) then !Using Eq. 17 of Comput. Theor. Chem., 1006, 92 (2013) and assume miu=0
    paulipot=totesp(x,y,z) -DFTxcpot(x,y,z,0) -weizpot(x,y,z) !Note that the sign of ESP in shubin's CTC paper is inversed
else if (ispecial==1) then !Strict definition
    paulipot=KEDpot(x,y,z)-weizpot(x,y,z)
end if
end function

!!------ The magnitude of Pauli force, namely the gradient norm of negative Pauli potential
real*8 function pauliforce(x,y,z)
real*8 x,y,z
diff=2D-5
forcex=-(paulipot(x+diff,y,z)-paulipot(x-diff,y,z))/(2*diff)
forcey=-(paulipot(x,y+diff,z)-paulipot(x,y-diff,z))/(2*diff)
forcez=-(paulipot(x,y,z+diff)-paulipot(x,y,z-diff))/(2*diff)
pauliforce=dsqrt(forcex**2+forcey**2+forcez**2)
end function

!!------ Pauli charge, =lapl(Pauli potential)/(-4*pi)
real*8 function paulicharge(x,y,z)
real*8 x,y,z
diff=2D-4 !Should not be smaller, otherwise some dirty points will occur
value=paulipot(x,y,z)
valuexaddadd=paulipot(x+2*diff,y,z)
valuexminmin=paulipot(x-2*diff,y,z)
valueyaddadd=paulipot(x,y+2*diff,z)
valueyminmin=paulipot(x,y-2*diff,z)
valuezaddadd=paulipot(x,y,z+2*diff)
valuezminmin=paulipot(x,y,z-2*diff)
xcomp=(valuexaddadd-2*value+valuexminmin)/(2*diff)**2
ycomp=(valueyaddadd-2*value+valueyminmin)/(2*diff)**2
zcomp=(valuezaddadd-2*value+valuezminmin)/(2*diff)**2
paulicharge=(xcomp+ycomp+zcomp)/(-4*pi)
end function

!!------ Quantum potential
real*8 function quantumpot(x,y,z)
real*8 x,y,z
if (ispecial==0) then !Based on Pauli potential by assuming miu=0 (Eq. 17 of Comput. Theor. Chem., 1006, 92 (2013))
    quantumpot=totesp(x,y,z)-weizpot(x,y,z)
else if (ispecial==1) then !Based on strict definition of Pauli potential
    quantumpot=paulipot(x,y,z)+DFTxcpot(x,y,z,0)
end if
end function

!!------ The magnitude of quantum force, namely the gradient norm of quantum potential
real*8 function quantumforce(x,y,z)
real*8 x,y,z
diff=2D-5
forcex=-(quantumpot(x+diff,y,z)-quantumpot(x-diff,y,z))/(2*diff)
forcey=-(quantumpot(x,y+diff,z)-quantumpot(x,y-diff,z))/(2*diff)
forcez=-(quantumpot(x,y,z+diff)-quantumpot(x,y,z-diff))/(2*diff)
quantumforce=dsqrt(forcex**2+forcey**2+forcez**2)
end function

!!------ Quantum charge
real*8 function quantumcharge(x,y,z)
real*8 x,y,z
diff=2D-4 !Should not be smaller, otherwise some dirty points will occur
value=quantumpot(x,y,z)
valuexaddadd=quantumpot(x+2*diff,y,z)
valuexminmin=quantumpot(x-2*diff,y,z)
valueyaddadd=quantumpot(x,y+2*diff,z)
valueyminmin=quantumpot(x,y-2*diff,z)
valuezaddadd=quantumpot(x,y,z+2*diff)
valuezminmin=quantumpot(x,y,z-2*diff)
xcomp=(valuexaddadd-2*value+valuexminmin)/(2*diff)**2
ycomp=(valueyaddadd-2*value+valueyminmin)/(2*diff)**2
zcomp=(valuezaddadd-2*value+valuezminmin)/(2*diff)**2
quantumcharge=(xcomp+ycomp+zcomp)/(-4*pi)
end function

!!------ The magnitude of electrostatic force, namely the gradient norm of electrostatic potential
real*8 function elestatforce(x,y,z)
real*8 x,y,z
diff=2D-5
forcex=-(-totesp(x+diff,y,z)+totesp(x-diff,y,z))/(2*diff)
forcey=-(-totesp(x,y+diff,z)+totesp(x,y-diff,z))/(2*diff)
forcez=-(-totesp(x,y,z+diff)+totesp(x,y,z-diff))/(2*diff)
elestatforce=dsqrt(forcex**2+forcey**2+forcez**2)
end function

!!------ Electrostatic charge
real*8 function elestatcharge(x,y,z)
real*8 x,y,z
diff=2D-4 !Should not be smaller, otherwise some dirty points will be presented
value=-totesp(x,y,z)
valuexaddadd=-totesp(x+2*diff,y,z)
valuexminmin=-totesp(x-2*diff,y,z)
valueyaddadd=-totesp(x,y+2*diff,z)
valueyminmin=-totesp(x,y-2*diff,z)
valuezaddadd=-totesp(x,y,z+2*diff)
valuezminmin=-totesp(x,y,z-2*diff)
xcomp=(valuexaddadd-2*value+valuexminmin)/(2*diff)**2
ycomp=(valueyaddadd-2*value+valueyminmin)/(2*diff)**2
zcomp=(valuezaddadd-2*value+valuezminmin)/(2*diff)**2
elestatcharge=(xcomp+ycomp+zcomp)/(-4*pi)
end function




!!-------- Calculate magnitude of vector sum of steric force, electrostatic force and quantum force
real*8 function SBLallforce(x,y,z)
real*8 x,y,z,sterderv(3)
call stericderv(x,y,z,sterderv)

diff=2D-5
eforcex=-(-totesp(x+diff,y,z)+totesp(x-diff,y,z))/(2*diff)
eforcey=-(-totesp(x,y+diff,z)+totesp(x,y-diff,z))/(2*diff)
eforcez=-(-totesp(x,y,z+diff)+totesp(x,y,z-diff))/(2*diff)

qforcex=-(quantumpot(x+diff,y,z)-quantumpot(x-diff,y,z))/(2*diff)
qforcey=-(quantumpot(x,y+diff,z)-quantumpot(x,y-diff,z))/(2*diff)
qforcez=-(quantumpot(x,y,z+diff)-quantumpot(x,y,z-diff))/(2*diff)

fx=-sterderv(1) + eforcex + qforcex
fy=-sterderv(2) + eforcey + qforcey
fz=-sterderv(3) + eforcez + qforcez

SBLallforce=dsqrt(fx**2+fy**2+fz**2)
end function




!!------- Use trilinear interpolation to obtain value at a given point by using cubmat. Both orthogonal and nonorthogonal grids are supported. PBC can be considered
!itype==1: interpolate from cubmat, =2: from cubmattmp
!Ref.: https://en.wikipedia.org/wiki/Trilinear_interpolation
!Trilinear interpolation is equivalent to linear interpolation between two bilinear interpolations
!In this case, the cell defining fractional coordinate corresponds to the box containing all grids. Origin of fractional coordinate is (orgx,orgy,orgz)
real*8 function linintp3d(x,y,z,itype)
real*8 x,y,z
integer itype
real*8 Cart(3),fract(3),tmpvec(3)
real*8 i1_val,i2_val,i3_val,i1_low,i2_low,i3_low,i1_high,i2_high,i3_high,d1,d2,d3

linintp3d=0
if (itype==1.and.(.not.allocated(cubmat))) return
if (itype==2.and.(.not.allocated(cubmattmp))) return

!Get fractional coordinate of present position
Cart(1)=x;Cart(2)=y;Cart(3)=z
if (ifPBC==3) call move_to_cell(Cart,Cart)
call Cart2fract_grid(Cart,fract)
i1_val=fract(1);i2_val=fract(2);i3_val=fract(3)
d1=1D0/nx !Grid spacing of fractional coordinate in each direction
d2=1D0/ny
d3=1D0/nz
do i1=1,nx
	i1_low=(i1-1)*d1
	i1_high=i1_low+d1
	if (i1_val>=i1_low.and.i1_val<i1_high) exit
end do
do i2=1,ny
	i2_low=(i2-1)*d2
	i2_high=i2_low+d2
	if (i2_val>=i2_low.and.i2_val<i2_high) exit
end do
do i3=1,nz
	i3_low=(i3-1)*d3
	i3_high=i3_low+d3
	if (i3_val>=i3_low.and.i3_val<i3_high) exit
end do

if (ifPBC==0) then
	if (i1<nx.and.i2<ny.and.i3<nz) then
		!Perform two bilinear interpolations first, then further linear interpolation
		if (itype==1) then
			val12_low= ( cubmat(i1,i2,i3    )*(i1_high-i1_val)*(i2_high-i2_val) + cubmat(i1+1,i2,i3    )*(i1_val-i1_low)*(i2_high-i2_val) + &
						 cubmat(i1,i2+1,i3  )*(i1_high-i1_val)*(i2_val-i2_low ) + cubmat(i1+1,i2+1,i3  )*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
			val12_high=( cubmat(i1,i2,i3+1  )*(i1_high-i1_val)*(i2_high-i2_val) + cubmat(i1+1,i2,i3+1  )*(i1_val-i1_low)*(i2_high-i2_val) + &
						 cubmat(i1,i2+1,i3+1)*(i1_high-i1_val)*(i2_val-i2_low ) + cubmat(i1+1,i2+1,i3+1)*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
		else
			val12_low= ( cubmattmp(i1,i2,i3    )*(i1_high-i1_val)*(i2_high-i2_val) + cubmattmp(i1+1,i2,i3    )*(i1_val-i1_low)*(i2_high-i2_val) + &
						 cubmattmp(i1,i2+1,i3  )*(i1_high-i1_val)*(i2_val-i2_low ) + cubmattmp(i1+1,i2+1,i3  )*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
			val12_high=( cubmattmp(i1,i2,i3+1  )*(i1_high-i1_val)*(i2_high-i2_val) + cubmattmp(i1+1,i2,i3+1  )*(i1_val-i1_low)*(i2_high-i2_val) + &
						 cubmattmp(i1,i2+1,i3+1)*(i1_high-i1_val)*(i2_val-i2_low ) + cubmattmp(i1+1,i2+1,i3+1)*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
		end if
		linintp3d=val12_low+(i3_val-i3_low)*(val12_high-val12_low)/d3
	else !Out of grid data range
		linintp3d=0D0
	end if
else !PBC
	i1p1=i1+1
    if (i1==nx) i1p1=1
	i2p1=i2+1
    if (i2==ny) i2p1=1
	i3p1=i3+1
    if (i3==nz) i3p1=1
	if (itype==1) then
		val12_low= ( cubmat(i1,i2,i3    )*(i1_high-i1_val)*(i2_high-i2_val) + cubmat(i1p1,i2,i3    )*(i1_val-i1_low)*(i2_high-i2_val) + &
						cubmat(i1,i2p1,i3  )*(i1_high-i1_val)*(i2_val-i2_low ) + cubmat(i1p1,i2p1,i3  )*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
		val12_high=( cubmat(i1,i2,i3p1  )*(i1_high-i1_val)*(i2_high-i2_val) + cubmat(i1p1,i2,i3p1  )*(i1_val-i1_low)*(i2_high-i2_val) + &
						cubmat(i1,i2p1,i3p1)*(i1_high-i1_val)*(i2_val-i2_low ) + cubmat(i1p1,i2p1,i3p1)*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
	else
		val12_low= ( cubmattmp(i1,i2,i3    )*(i1_high-i1_val)*(i2_high-i2_val) + cubmattmp(i1p1,i2,i3    )*(i1_val-i1_low)*(i2_high-i2_val) + &
						cubmattmp(i1,i2p1,i3  )*(i1_high-i1_val)*(i2_val-i2_low ) + cubmattmp(i1p1,i2p1,i3  )*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
		val12_high=( cubmattmp(i1,i2,i3p1  )*(i1_high-i1_val)*(i2_high-i2_val) + cubmattmp(i1p1,i2,i3p1  )*(i1_val-i1_low)*(i2_high-i2_val) + &
						cubmattmp(i1,i2p1,i3p1)*(i1_high-i1_val)*(i2_val-i2_low ) + cubmattmp(i1p1,i2p1,i3p1)*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
	end if
	linintp3d=val12_low+(i3_val-i3_low)*(val12_high-val12_low)/d3
end if
end function




!!------- Trilinear interpolation of 3D-vector field by using cubmatvec. Both orthogonal and nonorthogonal grids are supported. PBC can be considered
!See comment in subroutine linintp3d for more information
subroutine linintp3dvec(x,y,z,vecintp)
real*8 x,y,z,vecintp(3),val12_low(3),val12_high(3)
real*8 Cart(3),fract(3),tmpvec(3)
real*8 i1_val,i2_val,i3_val,i1_low,i2_low,i3_low,i1_high,i2_high,i3_high,d1,d2,d3

Cart(1)=x;Cart(2)=y;Cart(3)=z
if (ifPBC==3) call move_to_cell(Cart,Cart)
call Cart2fract_grid(Cart,fract)
i1_val=fract(1);i2_val=fract(2);i3_val=fract(3)
d1=1D0/nx !Grid spacing of fractional coordinate in each direction
d2=1D0/ny
d3=1D0/nz
do i1=1,nx
	i1_low=(i1-1)*d1
	i1_high=i1_low+d1
	if (i1_val>=i1_low.and.i1_val<i1_high) exit
end do
do i2=1,ny
	i2_low=(i2-1)*d2
	i2_high=i2_low+d2
	if (i2_val>=i2_low.and.i2_val<i2_high) exit
end do
do i3=1,nz
	i3_low=(i3-1)*d3
	i3_high=i3_low+d3
	if (i3_val>=i3_low.and.i3_val<i3_high) exit
end do

if (ifPBC==0) then
	if (i1<nx.and.i2<ny.and.i3<nz) then
		val12_low(:)= ( cubmatvec(:,i1,i2,i3    )*(i1_high-i1_val)*(i2_high-i2_val) + cubmatvec(:,i1+1,i2,i3    )*(i1_val-i1_low)*(i2_high-i2_val) + &
						cubmatvec(:,i1,i2+1,i3  )*(i1_high-i1_val)*(i2_val-i2_low ) + cubmatvec(:,i1+1,i2+1,i3  )*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
		val12_high(:)=( cubmatvec(:,i1,i2,i3+1  )*(i1_high-i1_val)*(i2_high-i2_val) + cubmatvec(:,i1+1,i2,i3+1  )*(i1_val-i1_low)*(i2_high-i2_val) + &
						cubmatvec(:,i1,i2+1,i3+1)*(i1_high-i1_val)*(i2_val-i2_low ) + cubmatvec(:,i1+1,i2+1,i3+1)*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
		vecintp(:)=val12_low(:)+(i3_val-i3_low)*(val12_high(:)-val12_low(:))/d3
	else !Out of grid data range
		vecintp=0D0
	end if
else !PBC
	i1p1=i1+1
    if (i1==nx) i1p1=1
	i2p1=i2+1
    if (i2==ny) i2p1=1
	i3p1=i3+1
    if (i3==nz) i3p1=1
	val12_low(:)= ( cubmatvec(:,i1,i2,i3    )*(i1_high-i1_val)*(i2_high-i2_val) + cubmatvec(:,i1p1,i2,i3    )*(i1_val-i1_low)*(i2_high-i2_val) + &
					cubmatvec(:,i1,i2p1,i3  )*(i1_high-i1_val)*(i2_val-i2_low ) + cubmatvec(:,i1p1,i2p1,i3  )*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
	val12_high(:)=( cubmatvec(:,i1,i2,i3p1  )*(i1_high-i1_val)*(i2_high-i2_val) + cubmatvec(:,i1p1,i2,i3p1  )*(i1_val-i1_low)*(i2_high-i2_val) + &
					cubmatvec(:,i1,i2p1,i3p1)*(i1_high-i1_val)*(i2_val-i2_low ) + cubmatvec(:,i1p1,i2p1,i3p1)*(i1_val-i1_low)*(i2_val-i2_low ) ) /(d1*d2)
	vecintp(:)=val12_low(:)+(i3_val-i3_low)*(val12_high(:)-val12_low(:))/d3
end if
end subroutine




!!-------- Use cubic spline interpolation to obtain value at a given point by using cubmat. Both orthogonal and nonorthogonal grids are supported. PBC can be considered
!itype==1: interpolate from cubmat, =2: from cubmattmp
!If all grid points in cubmat are used in the interpolation, the cost will be extremely high if very large number of points are to be calculated, &
!therefore I decide only takes a few grids around the present position for the interpolation, and this idea works well and the cost is significantly lowered
!The spline interpolation code comes from https://github.com/jacobwilliams/bspline-fortran, the Bspline.f90 is slightly adapted, &
!See comment of subroutine db3ink and db3val to understand the use of the routines
real*8 function splineintp3D(x,y,z,itype)
use bspline_sub_module
real*8 x,y,z
real*8 Cart(3),fract(3),tmpvec(3)
real*8 i1_val,i2_val,i3_val,i1_low,i2_low,i3_low,i1_high,i2_high,i3_high,d1,d2,d3
integer itype
integer,parameter :: norder=4 !4 corresponds to cubic spline interpolation. order = polynomial degree + 1, order from 3 to 6 are supported
!According to my experience, norder=4~6 has negligible difference. Even norder=3 is quite similar to norder=4. The norder seems doesn't affect cost
integer,parameter :: k1=norder,k2=norder,k3=norder
integer :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1 !Must be set to 1 before call as mentioned in the test file, though I don't know why
integer :: iknot=0  !Automatically determine the knots
integer,parameter :: next=4 !The number of extended grids on both sides. Test showed that 4 leads to converge result, increasing it will not detectably affect result
integer,parameter :: nlocgrd=2+2*next
real*8 :: arr1(nlocgrd),arr2(nlocgrd),arr3(nlocgrd)
real*8 t1(nlocgrd+k1),t2(nlocgrd+k2),t3(nlocgrd+k3) !Not used because I let knot automatically determined, but these arrays must be defined
real*8 cubloc(nlocgrd,nlocgrd,nlocgrd),spline3Dcoeff(nlocgrd,nlocgrd,nlocgrd)

splineintp3D=0
if (itype==1.and.(.not.allocated(cubmat))) return
if (itype==2.and.(.not.allocated(cubmattmp))) return

Cart(1)=x;Cart(2)=y;Cart(3)=z
if (ifPBC==3) call move_to_cell(Cart,Cart)
call Cart2fract_grid(Cart,fract)
i1_val=fract(1);i2_val=fract(2);i3_val=fract(3)
d1=1D0/nx !Grid spacing of fractional coordinate in each direction
d2=1D0/ny
d3=1D0/nz
do i1=1,nx
	i1_low=(i1-1)*d1
	i1_high=i1_low+d1
	if (i1_val>=i1_low.and.i1_val<i1_high) exit
end do
do i2=1,ny
	i2_low=(i2-1)*d2
	i2_high=i2_low+d2
	if (i2_val>=i2_low.and.i2_val<i2_high) exit
end do
do i3=1,nz
	i3_low=(i3-1)*d3
	i3_high=i3_low+d3
	if (i3_val>=i3_low.and.i3_val<i3_high) exit
end do

if (ifPBC==0) then
	if (i1+1>nx-next.or.i2+1>ny-next.or.i3+1>nz-next .or. i1<next+1.or.i2<next+1.or.i3<next+1) then !Out of grid data range
		splineintp3D=0D0
		return
	end if
end if

!Determine the index range of the grids around present point, these grids will be actually used in interpolation
ilow1=i1-next
ihigh1=i1+1+next
ilow2=i2-next
ihigh2=i2+1+next
ilow3=i3-next
ihigh3=i3+1+next
do i1=ilow1,ihigh1
    arr1(i1-ilow1+1)=(i1-1)*d1
end do
do i2=ilow2,ihigh2
    arr2(i2-ilow2+1)=(i2-1)*d2
end do
do i3=ilow3,ihigh3
    arr3(i3-ilow3+1)=(i3-1)*d3
end do
do i1=ilow1,ihigh1
    do i2=ilow2,ihigh2
        do i3=ilow3,ihigh3
			ii1=i1
            ii2=i2
            ii3=i3
			if (ifPBC/=0) call PBCgrididx(ii1,ii2,ii3)
			if (itype==1) cubloc(i1-ilow1+1,i2-ilow2+1,i3-ilow3+1)=cubmat(ii1,ii2,ii3)
			if (itype==2) cubloc(i1-ilow1+1,i2-ilow2+1,i3-ilow3+1)=cubmattmp(ii1,ii2,ii3)
        end do
    end do
end do

!Initialize
call db3ink(arr1,nlocgrd,arr2,nlocgrd,arr3,nlocgrd,cubloc,k1,k2,k3,iknot,t1,t2,t3,spline3Dcoeff,iflag)
if (iflag/=0) write(*,*) "Error when initializing B-spline!"

!Evaluate value
id1=0;id2=0;id3=0 !Do not evaluate derivative
call db3val(i1_val,i2_val,i3_val,id1,id2,id3,t1,t2,t3,nlocgrd,nlocgrd,nlocgrd,k1,k2,k3,spline3Dcoeff,splineintp3D,iflag,inbvx,inbvy,inbvz,iloy,iloz)
end function




!!-------- Calculate various kinds of integrand of DFT exchange-correlation functionals
!The routines are provided by DFT repository (ftp://ftp.dl.ac.uk/qcg/dft_library/index.html)
!The global variable "iDFTxcsel" is used to select the XC functional, see manual
!Note that the inner core density represented by EDF field is not taken into account
real*8 function DFTxcfunc(x,y,z)
real*8 x,y,z
if (wfntype==0.or.wfntype==3) then !Closed-shell
	DFTxcfunc=DFTxcfunc_close(x,y,z)
else !Open-shell
	DFTxcfunc=DFTxcfunc_open(x,y,z)
end if
end function
!---- Calculate various kinds of DFT exchange-correlation potentials, see the comment of DFTxcfunc
!ispin=0: Closed-shell form, =1: Alpha spin, =2: Beta spin
real*8 function DFTxcpot(x,y,z,ispin)
real*8 x,y,z
integer ispin
if ((wfntype==0.or.wfntype==3).and.ispin/=0) then !Closed-shell system but calculate alpha or beta spin
	write(*,"(a)") " Error: This is a closed-shell system, but you request to calculate potential of alpha or beta spin, &
    this is not meaningless. Please use the 1100th user-defined function instead"
    write(*,*) "Press ENTER button to exit program"
    read(*,*)
    stop
end if
if ((wfntype==1.or.wfntype==2.or.wfntype==4).and.ispin==0) then !Open-shell system but calculate all-electron potential
	write(*,"(a)") " Error: This is an open-shell system, but you request to calculate potential of all electrons, &
    this is not meaningless. Please use the 1101th or 1102th user-defined function instead"
    write(*,*) "Press ENTER button to exit program"
    read(*,*)
    stop
end if
	
if (ispin==0) then !Closed-shell
	DFTxcpot=DFTxcpot_close(x,y,z)
else !Open-shell
	DFTxcpot=DFTxcpot_open(x,y,z,ispin)
end if
end function




!!----- Closed-shell form of DFTxcfunc routine
real*8 function DFTxcfunc_close(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
call getXCdata_close(0,tdens,tgrad**2,DFTxcfunc_close,rnouse,rnouse,rnouse,rnouse,rnouse)
end function

!!----- Closed-shell form of DFTxcpot routine
real*8 function DFTxcpot_close(x,y,z)
real*8 x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),gradrho(3),hessrho(3,3),tmparr(3,1),tmpval(1,1),lapltot
rho=0D0
gradrho=0D0
call orbderv(4,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
do i=1,nmo
	rho=rho+MOocc(i)*wfnval(i)**2
	gradrho(:)=gradrho(:)+MOocc(i)*wfnval(i)*wfnderv(:,i)
end do
gradrho=2*gradrho
hessrho(1,1)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)**2 + wfnval(1:nmo)*wfnhess(1,1,1:nmo) ) )
hessrho(2,2)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)**2 + wfnval(1:nmo)*wfnhess(2,2,1:nmo) ) )
hessrho(3,3)=2*sum( MOocc(1:nmo)*( wfnderv(3,1:nmo)**2 + wfnval(1:nmo)*wfnhess(3,3,1:nmo) ) )
hessrho(1,2)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(2,1:nmo)+wfnhess(1,2,1:nmo)*wfnval(1:nmo) ) )
hessrho(2,3)=2*sum( MOocc(1:nmo)*( wfnderv(2,1:nmo)*wfnderv(3,1:nmo)+wfnhess(2,3,1:nmo)*wfnval(1:nmo) ) )
hessrho(1,3)=2*sum( MOocc(1:nmo)*( wfnderv(1,1:nmo)*wfnderv(3,1:nmo)+wfnhess(1,3,1:nmo)*wfnval(1:nmo) ) )
hessrho(2,1)=hessrho(1,2)
hessrho(3,2)=hessrho(2,3)
hessrho(3,1)=hessrho(1,3)
lapltot=hessrho(1,1)+hessrho(2,2)+hessrho(3,3)
sigma=sum(gradrho(:)**2)
tmparr(:,1)=gradrho
tmpval=matmul(2*matmul(transpose(tmparr),hessrho),tmparr) !dot product between grad(sigma) and grad(rho)
call getXCdata_close(2,rho,sigma,value,d1rho,d1sig,d2rho,d2rhosig,d2sig)
DFTxcpot_close=d1rho-2*(d2rhosig*sigma+d2sig*tmpval(1,1)+d1sig*lapltot)
end function




!!---- For closed-shell cases. Input rho and gradrho^2, return the value and its derivative of selected XC (or X/C only) functional
!The global variable "iDFTxcsel" is used to select the XC functional, see manual
!ixcderv=0: only get value, =1: also get d1rho and d1sig, =2: also get d2rho, d2rhosig and d2sig
!rho, sigma: inputted rho and gradrho^2
!value: outputted integrand of functional
!d1rho, d1sig: 1st derivative of functional w.r.t. rho and sigma, respectively
!d2rho, d2sig: 2nd derivative of functional w.r.t. rho and sigma, respectively
!d2rhosig: 1st derv w.r.t. rho and 1st derv w.r.t. sigma
subroutine getXCdata_close(ixcderv,rho,sigma,value,d1rho,d1sig,d2rho,d2rhosig,d2sig)
integer ixcderv
real*8 rho,sigma,value,d1rho,d1sig,d2rho,d2rhosig,d2sig,rhoa1(1),sigmaaa1(1)
real*8 XCzk(1),Xzk(1),Czk(1),XCvrhoa(1),Xvrhoa(1),Cvrhoa(1),XCvsigmaaa(1),Xvsigmaaa(1),Cvsigmaaa(1)
real*8 XCv2rhoa2(1),Xv2rhoa2(1),Cv2rhoa2(1),XCv2rhoasigmaaa(1),Xv2rhoasigmaaa(1),Cv2rhoasigmaaa(1)
real*8 XCv2sigmaaa2(1),Xv2sigmaaa2(1),Cv2sigmaaa2(1)
rhoa1(1)=rho
sigmaaa1(1)=sigma
!X part
if (iDFTxcsel==0.or.iDFTxcsel==80) then
	call rks_x_lda(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
else if (iDFTxcsel==1.or.iDFTxcsel==81.or.iDFTxcsel==82.or.iDFTxcsel==83) then
	call rks_x_b88(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
else if (iDFTxcsel==2.or.iDFTxcsel==84) then
	call rks_x_pbe(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
else if (iDFTxcsel==3.or.iDFTxcsel==85) then
	call rks_x_pw91(ixcderv,1,rhoa1,sigmaaa1,Xzk,Xvrhoa,Xvsigmaaa,Xv2rhoa2,Xv2rhoasigmaaa,Xv2sigmaaa2)
end if
!C part
if (iDFTxcsel==30.or.iDFTxcsel==80) then
	call rks_c_vwn5(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==31.or.iDFTxcsel==81) then
	call rks_c_p86(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==32.or.iDFTxcsel==82) then
	call rks_c_lyp(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==33.or.iDFTxcsel==83.or.iDFTxcsel==85) then
	call rks_c_pw91(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
else if (iDFTxcsel==34.or.iDFTxcsel==84) then
	call rks_c_pbe(ixcderv,1,rhoa1,sigmaaa1,Czk,Cvrhoa,Cvsigmaaa,Cv2rhoa2,Cv2rhoasigmaaa,Cv2sigmaaa2)
end if
!Whole XC
if (iDFTxcsel==70) then
	call rks_xc_b97(ixcderv,1,rhoa1,sigmaaa1,XCzk,XCvrhoa,XCvsigmaaa,XCv2rhoa2,XCv2rhoasigmaaa,XCv2sigmaaa2)
else if (iDFTxcsel==71) then
	call rks_xc_hcth407(ixcderv,1,rhoa1,sigmaaa1,XCzk,XCvrhoa,XCvsigmaaa,XCv2rhoa2,XCv2rhoasigmaaa,XCv2sigmaaa2)
else if (iDFTxcsel>=80.and.iDFTxcsel<99) then
	XCzk=Xzk+Czk
	XCvrhoa=Xvrhoa+Cvrhoa
	XCvsigmaaa=Xvsigmaaa+Cvsigmaaa
	XCv2rhoa2=Xv2rhoa2+Cv2rhoa2
	XCv2rhoasigmaaa=Xv2rhoasigmaaa+Cv2rhoasigmaaa
	XCv2sigmaaa2=Xv2sigmaaa2+Cv2sigmaaa2
end if
!Note that the problem of derivative of the DFT repository is revised by dividing a factor, similarly hereinafter
if (iDFTxcsel<30) then
	value=Xzk(1)
	d1rho=Xvrhoa(1)
	d1sig=Xvsigmaaa(1)/4D0
	d2rho=Xv2rhoa2(1)/2D0
	d2rhosig=Xv2rhoasigmaaa(1)/4D0
	d2sig=Xv2sigmaaa2(1)/16D0
else if (iDFTxcsel<70) then
	value=Czk(1)
	d1rho=Cvrhoa(1)
	d1sig=Cvsigmaaa(1)/4D0
	d2rho=Cv2rhoa2(1)/2D0
	d2rhosig=Cv2rhoasigmaaa(1)/4D0
	d2sig=Cv2sigmaaa2(1)/16D0
else if (iDFTxcsel<100) then
	value=XCzk(1)
	d1rho=XCvrhoa(1)
	d1sig=XCvsigmaaa(1)/4D0
	d2rho=XCv2rhoa2(1)/2D0
	d2rhosig=XCv2rhoasigmaaa(1)/4D0
	d2sig=XCv2sigmaaa2(1)/16D0
end if
end subroutine



!!------ Open-shell form of DFTxcfunc routine
real*8 function DFTxcfunc_open(x,y,z)
real*8 x,y,z
call gendensgradab(x,y,z,adens,bdens,tdens,agrad,bgrad,tgrad,abgrad)
call getXCdata_open(0,adens,bdens,agrad**2,bgrad**2,abgrad,DFTxcfunc_open,&
sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb,sb)
end function


!!------ Open-shell form of DFTxcpot routine
!For simplicity, this function uses numerical way to evaluate some intermediate derivatives
!ispin=1: XC potential of alpha spin, =2: beta spin
real*8 function DFTxcpot_open(x,y,z,ispin)
integer ispin
real*8 x,y,z
real*8 rhoa,rhob,gradrhoa(3),gradrhob(3),laplrhoa,laplrhob
real*8 tmpvec1(3),tmpvec2(3),rhoa_tmp,rhob_tmp,gradrhoa_tmp(3),gradrhob_tmp(3)

diff=1E-5
call gendens_gradvec_lapl_ab(x,y,z,rhoa,rhob,gradrhoa,gradrhob,sigaa,sigbb,sigab,laplrhoa,laplrhob)
call getXCdata_open(1,rhoa,rhob,sigaa,sigbb,sigab,s,d1rhoa,d1rhob,d1sigaa,d1sigbb,d1sigab,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s) !Useless data assign to s
!write(*,*) sigaa,sigbb,sigab
!write(*,*) gradrhoa
!write(*,*) gradrhob
!write(*,*) gradrhoa+gradrhob
if (ispin==1) then !Alpha
	!tmpvec1(i)= Derivative of d1sigaa in direction i
	!tmpvec2(i)= Derivative of d1sigab in direction i
	!dx
	call gendens_gradvec_lapl_ab(x+diff,y,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_xadd,d1sigbb_tmp,d1sigab_xadd,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	call gendens_gradvec_lapl_ab(x-diff,y,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_xmin,d1sigbb_tmp,d1sigab_xmin,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	tmpvec1(1)=(d1sigaa_xadd-d1sigaa_xmin)/(2*diff)
	tmpvec2(1)=(d1sigab_xadd-d1sigab_xmin)/(2*diff)
	!dy
	call gendens_gradvec_lapl_ab(x,y+diff,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_yadd,d1sigbb_tmp,d1sigab_yadd,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	call gendens_gradvec_lapl_ab(x,y-diff,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_ymin,d1sigbb_tmp,d1sigab_ymin,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	tmpvec1(2)=(d1sigaa_yadd-d1sigaa_ymin)/(2*diff)
	tmpvec2(2)=(d1sigab_yadd-d1sigab_ymin)/(2*diff)
	!dz
	call gendens_gradvec_lapl_ab(x,y,z+diff,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_zadd,d1sigbb_tmp,d1sigab_zadd,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	call gendens_gradvec_lapl_ab(x,y,z-diff,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_zmin,d1sigbb_tmp,d1sigab_zmin,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	tmpvec1(3)=(d1sigaa_zadd-d1sigaa_zmin)/(2*diff)
	tmpvec2(3)=(d1sigab_zadd-d1sigab_zmin)/(2*diff)

	DFTxcpot_open = d1rhoa - 2*(sum(tmpvec1(:)*gradrhoa(:)) + d1sigaa*laplrhoa) - (sum(tmpvec2(:)*gradrhob(:)) + d1sigab*laplrhob )
else if (ispin==2) then !Beta
	!tmpvec1(i)= Derivative of d1sigbb in direction i
	!tmpvec2(i)= Derivative of d1sigab in direction i
	!dx
	call gendens_gradvec_lapl_ab(x+diff,y,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_tmp,d1sigbb_xadd,d1sigab_xadd,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	call gendens_gradvec_lapl_ab(x-diff,y,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_tmp,d1sigbb_xmin,d1sigab_xmin,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	tmpvec1(1)=(d1sigbb_xadd-d1sigbb_xmin)/(2*diff)
	tmpvec2(1)=(d1sigab_xadd-d1sigab_xmin)/(2*diff)
	!dy
	call gendens_gradvec_lapl_ab(x,y+diff,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_tmp,d1sigbb_yadd,d1sigab_yadd,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	call gendens_gradvec_lapl_ab(x,y-diff,z,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_tmp,d1sigbb_ymin,d1sigab_ymin,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	tmpvec1(2)=(d1sigbb_yadd-d1sigbb_ymin)/(2*diff)
	tmpvec2(2)=(d1sigab_yadd-d1sigab_ymin)/(2*diff)
	!dz
	call gendens_gradvec_lapl_ab(x,y,z+diff,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_tmp,d1sigbb_zadd,d1sigab_zadd,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	call gendens_gradvec_lapl_ab(x,y,z-diff,rhoa_tmp,rhob_tmp,gradrhoa_tmp,gradrhob_tmp,sigaa,sigbb,sigab,rnouse,rnouse)
	call getXCdata_open(1,rhoa_tmp,rhob_tmp,sigaa,sigbb,sigab,s,s,s,d1sigaa_tmp,d1sigbb_zmin,d1sigab_zmin,s,s,s,s,s,s,s,s,s,s,s,s,s,s,s)
	tmpvec1(3)=(d1sigbb_zadd-d1sigbb_zmin)/(2*diff)
	tmpvec2(3)=(d1sigab_zadd-d1sigab_zmin)/(2*diff)

	DFTxcpot_open = d1rhob - 2*(sum(tmpvec1(:)*gradrhob(:)) + d1sigbb*laplrhob) - (sum(tmpvec2(:)*gradrhoa(:)) + d1sigab*laplrhoa )
end if
end function


!!---- For open-shell cases. Input rho and gradrho^2, return the value and its derivative of selected XC (or X/C only) functional
!The global variable "iDFTxcsel" is used to select the XC functional, see manual
!ixcderv=0: only get value, =1: also get 1st derv., =2: also get 2nd derv.
!Input quantities:
!rhoa=rho_alpha   rhob=rho_beta
!sigaa=|gradrho_alpha|^2  sigbb=|gradrho_beta|^2
!sigab=Dot product between gradrho_alpha and gradrho_beta
subroutine getXCdata_open(ixcderv,rhoa,rhob,sigaa,sigbb,sigab,value,d1rhoa,d1rhob,d1sigaa,d1sigbb,d1sigab,d2rhoaa,d2rhobb,d2rhoab,&
d2rhoasigaa,d2rhoasigab,d2rhoasigbb,d2rhobsigbb,d2rhobsigab,d2rhobsigaa,d2sigaaaa,d2sigaaab,d2sigaabb,d2sigabab,d2sigabbb,d2sigbbbb)
integer ixcderv
!Input arguments
real*8 rhoa,rhob,sigaa,sigbb,sigab,value,d1rhoa,d1rhob,d1sigaa,d1sigbb,d1sigab,d2rhoaa,d2rhobb,d2rhoab,&
d2rhoasigaa,d2rhoasigab,d2rhoasigbb,d2rhobsigbb,d2rhobsigab,d2rhobsigaa,d2sigaaaa,d2sigaaab,d2sigaabb,d2sigabab,d2sigabbb,d2sigbbbb
!Inputted information
real*8 rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1)
!Returned information
real*8 Xzk(1),Xvrhoa(1),Xvrhob(1),Xvsigmaaa(1),Xvsigmabb(1),Xvsigmaab(1),Xv2rhoa2(1),Xv2rhob2(1),Xv2rhoab(1)&
,Xv2rhoasigmaaa(1),Xv2rhoasigmaab(1),Xv2rhoasigmabb(1),Xv2rhobsigmabb(1),Xv2rhobsigmaab(1),Xv2rhobsigmaaa(1)&
,Xv2sigmaaa2(1),Xv2sigmaaaab(1),Xv2sigmaaabb(1),Xv2sigmaab2(1),Xv2sigmaabbb(1),Xv2sigmabb2(1)
real*8 Czk(1),Cvrhoa(1),Cvrhob(1),Cvsigmaaa(1),Cvsigmabb(1),Cvsigmaab(1),Cv2rhoa2(1),Cv2rhob2(1),Cv2rhoab(1)&
,Cv2rhoasigmaaa(1),Cv2rhoasigmaab(1),Cv2rhoasigmabb(1),Cv2rhobsigmabb(1),Cv2rhobsigmaab(1),Cv2rhobsigmaaa(1)&
,Cv2sigmaaa2(1),Cv2sigmaaaab(1),Cv2sigmaaabb(1),Cv2sigmaab2(1),Cv2sigmaabbb(1),Cv2sigmabb2(1)
real*8 XCzk(1),XCvrhoa(1),XCvrhob(1),XCvsigmaaa(1),XCvsigmabb(1),XCvsigmaab(1),XCv2rhoa2(1),XCv2rhob2(1),XCv2rhoab(1)&
,XCv2rhoasigmaaa(1),XCv2rhoasigmaab(1),XCv2rhoasigmabb(1),XCv2rhobsigmabb(1),XCv2rhobsigmaab(1),XCv2rhobsigmaaa(1)&
,XCv2sigmaaa2(1),XCv2sigmaaaab(1),XCv2sigmaaabb(1),XCv2sigmaab2(1),XCv2sigmaabbb(1),XCv2sigmabb2(1)

rhoa1(1)=rhoa
rhob1(1)=rhob
sigmaaa1(1)=sigaa
sigmabb1(1)=sigbb
sigmaab1(1)=sigab
!X part
if (iDFTxcsel==0.or.iDFTxcsel==80) then
	call uks_x_lda(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
	Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
	Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
else if (iDFTxcsel==1.or.iDFTxcsel==81.or.iDFTxcsel==82.or.iDFTxcsel==83) then
	call uks_x_b88(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
	Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
	Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
else if (iDFTxcsel==2.or.iDFTxcsel==84) then
	call uks_x_pbe(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
	Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
	Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
else if (iDFTxcsel==3.or.iDFTxcsel==85) then
	call uks_x_pw91(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Xzk,Xvrhoa,Xvrhob,Xvsigmaaa,Xvsigmabb,Xvsigmaab,Xv2rhoa2,Xv2rhob2,Xv2rhoab,&
	Xv2rhoasigmaaa,Xv2rhoasigmaab,Xv2rhoasigmabb,Xv2rhobsigmabb,Xv2rhobsigmaab,Xv2rhobsigmaaa,&
	Xv2sigmaaa2,Xv2sigmaaaab,Xv2sigmaaabb,Xv2sigmaab2,Xv2sigmaabbb,Xv2sigmabb2)
end if
!C part
if (iDFTxcsel==30.or.iDFTxcsel==80) then
	call uks_c_vwn5(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
	Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
	Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==31.or.iDFTxcsel==81) then
	call uks_c_p86(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
	Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
	Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==32.or.iDFTxcsel==82) then
	call uks_c_lyp(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
	Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
	Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==33.or.iDFTxcsel==83.or.iDFTxcsel==85) then
	call uks_c_pw91(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
	Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
	Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
else if (iDFTxcsel==34.or.iDFTxcsel==84) then
	call uks_c_pbe(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,Czk,Cvrhoa,Cvrhob,Cvsigmaaa,Cvsigmabb,Cvsigmaab,Cv2rhoa2,Cv2rhob2,Cv2rhoab,&
	Cv2rhoasigmaaa,Cv2rhoasigmaab,Cv2rhoasigmabb,Cv2rhobsigmabb,Cv2rhobsigmaab,Cv2rhobsigmaaa,&
	Cv2sigmaaa2,Cv2sigmaaaab,Cv2sigmaaabb,Cv2sigmaab2,Cv2sigmaabbb,Cv2sigmabb2)
end if
!Whole XC
if (iDFTxcsel==70) then
	call uks_xc_b97(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,XCzk,XCvrhoa,XCvrhob,XCvsigmaaa,XCvsigmabb,XCvsigmaab,XCv2rhoa2,XCv2rhob2,XCv2rhoab,&
	XCv2rhoasigmaaa,XCv2rhoasigmaab,XCv2rhoasigmabb,XCv2rhobsigmabb,XCv2rhobsigmaab,XCv2rhobsigmaaa,&
	XCv2sigmaaa2,XCv2sigmaaaab,XCv2sigmaaabb,XCv2sigmaab2,XCv2sigmaabbb,XCv2sigmabb2)
else if (iDFTxcsel==71) then
	call uks_xc_hcth407(ixcderv,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,XCzk,XCvrhoa,XCvrhob,XCvsigmaaa,XCvsigmabb,XCvsigmaab,XCv2rhoa2,XCv2rhob2,XCv2rhoab,&
	XCv2rhoasigmaaa,XCv2rhoasigmaab,XCv2rhoasigmabb,XCv2rhobsigmabb,XCv2rhobsigmaab,XCv2rhobsigmaaa,&
	XCv2sigmaaa2,XCv2sigmaaaab,XCv2sigmaaabb,XCv2sigmaab2,XCv2sigmaabbb,XCv2sigmabb2)
else if (iDFTxcsel>=80.and.iDFTxcsel<99) then
	XCzk=Xzk+Czk
	XCvrhoa=Xvrhoa+Cvrhoa
	XCvrhob=Xvrhob+Cvrhob
	XCvsigmaaa=Xvsigmaaa+Cvsigmaaa
	XCvsigmabb=Xvsigmabb+Cvsigmabb
	XCvsigmaab=Xvsigmaab+Cvsigmaab
	XCv2rhoa2=Xv2rhoa2+Cv2rhoa2
	XCv2rhob2=Xv2rhob2+Cv2rhob2
	XCv2rhoab=Xv2rhoab+Cv2rhoab
	XCv2rhoasigmaaa=Xv2rhoasigmaaa+Cv2rhoasigmaaa
	XCv2rhoasigmaab=Xv2rhoasigmaab+Cv2rhoasigmaab
	XCv2rhoasigmabb=Xv2rhoasigmabb+Cv2rhoasigmabb
	XCv2rhobsigmabb=Xv2rhobsigmabb+Cv2rhobsigmabb
	XCv2rhobsigmaab=Xv2rhobsigmaab+Cv2rhobsigmaab
	XCv2rhobsigmaaa=Xv2rhobsigmaaa+Cv2rhobsigmaaa
	XCv2sigmaaa2= Xv2sigmaaa2+Cv2sigmaaa2
	XCv2sigmaaaab=Xv2sigmaaaab+Cv2sigmaaaab
	XCv2sigmaaabb=Xv2sigmaaabb+Cv2sigmaaabb
	XCv2sigmaab2= Xv2sigmaab2+Cv2sigmaab2
	XCv2sigmaabbb=Xv2sigmaabbb+Cv2sigmaabbb
	XCv2sigmabb2= Xv2sigmabb2+Cv2sigmabb2
end if

if (iDFTxcsel<30) then
	value=Xzk(1)
	d1rhoa=Xvrhoa(1)
	d1rhob=Xvrhob(1)
	d1sigaa=Xvsigmaaa(1)
	d1sigbb=Xvsigmabb(1)
	d1sigab=Xvsigmaab(1)
	d2rhoaa=Xv2rhoa2(1)
	d2rhobb=Xv2rhob2(1)
	d2rhoab=Xv2rhoab(1)
	d2rhoasigaa=Xv2rhoasigmaaa(1)
	d2rhoasigab=Xv2rhoasigmaab(1)
	d2rhoasigbb=Xv2rhoasigmabb(1)
	d2rhobsigbb=Xv2rhobsigmabb(1)
	d2rhobsigab=Xv2rhobsigmaab(1)
	d2rhobsigaa=Xv2rhobsigmaaa(1)
	d2sigaaaa=Xv2sigmaaa2(1)
	d2sigaaab=Xv2sigmaaaab(1)
	d2sigaabb=Xv2sigmaaabb(1)
	d2sigabab=Xv2sigmaab2(1)
	d2sigabbb=Xv2sigmaabbb(1)
	d2sigbbbb=Xv2sigmabb2(1)
else if (iDFTxcsel<70) then
	value=Czk(1)
	d1rhoa=Cvrhoa(1)
	d1rhob=Cvrhob(1)
	d1sigaa=Cvsigmaaa(1)
	d1sigbb=Cvsigmabb(1)
	d1sigab=Cvsigmaab(1)
	d2rhoaa=Cv2rhoa2(1)
	d2rhobb=Cv2rhob2(1)
	d2rhoab=Cv2rhoab(1)
	d2rhoasigaa=Cv2rhoasigmaaa(1)
	d2rhoasigab=Cv2rhoasigmaab(1)
	d2rhoasigbb=Cv2rhoasigmabb(1)
	d2rhobsigbb=Cv2rhobsigmabb(1)
	d2rhobsigab=Cv2rhobsigmaab(1)
	d2rhobsigaa=Cv2rhobsigmaaa(1)
	d2sigaaaa=Cv2sigmaaa2(1)
	d2sigaaab=Cv2sigmaaaab(1)
	d2sigaabb=Cv2sigmaaabb(1)
	d2sigabab=Cv2sigmaab2(1)
	d2sigabbb=Cv2sigmaabbb(1)
	d2sigbbbb=Cv2sigmabb2(1)
else if (iDFTxcsel<100) then
	value=XCzk(1)
	d1rhoa=XCvrhoa(1)
	d1rhob=XCvrhob(1)
	d1sigaa=XCvsigmaaa(1)
	d1sigbb=XCvsigmabb(1)
	d1sigab=XCvsigmaab(1)
	d2rhoaa=XCv2rhoa2(1)
	d2rhobb=XCv2rhob2(1)
	d2rhoab=XCv2rhoab(1)
	d2rhoasigaa=XCv2rhoasigmaaa(1)
	d2rhoasigab=XCv2rhoasigmaab(1)
	d2rhoasigbb=XCv2rhoasigmabb(1)
	d2rhobsigbb=XCv2rhobsigmabb(1)
	d2rhobsigab=XCv2rhobsigmaab(1)
	d2rhobsigaa=XCv2rhobsigmaaa(1)
	d2sigaaaa=XCv2sigmaaa2(1)
	d2sigaaab=XCv2sigmaaaab(1)
	d2sigaabb=XCv2sigmaaabb(1)
	d2sigabab=XCv2sigmaab2(1)
	d2sigabbb=XCv2sigmaabbb(1)
	d2sigbbbb=XCv2sigmabb2(1)
end if
end subroutine




!!---- The distance in Angstrom from a point (x,y,z) to the nearest atom in the array
real*8 function surfana_di(x,y,z,nlen,atmlist)
real*8 x,y,z
integer nlen,atmlist(nlen)
dist2min=1D100
do iatm=1,ncenter
	if (any(atmlist==iatm)) then !The atom is in the list
		dist2=(a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2
		if (dist2<dist2min) dist2min=dist2
	end if
end do
surfana_di=dsqrt(dist2min)*b2a
end function

!!---- The distance in Angstrom from a point (x,y,z) to the nearest atom not in the array
real*8 function surfana_de(x,y,z,nlen,atmlist)
real*8 x,y,z
integer nlen,atmlist(nlen)
dist2min=1D100
do iatm=1,ncenter
	if (all(atmlist/=iatm)) then !The atom is not in the list
		dist2=(a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2
		if (dist2<dist2min) dist2min=dist2
	end if
end do
surfana_de=dsqrt(dist2min)*b2a
end function

!!---- Normalized contact distance (d_norm) in Angstrom. defined in terms of de, di and the vdW radii of the atoms
real*8 function surfana_norm(x,y,z,nlen,atmlist)
real*8 x,y,z
integer nlen,atmlist(nlen)
dist2minin=1D100 !The nearest distance to atoms inside
dist2minext=1D100 !The nearest distance to atoms outside
iminin=0
iminext=0
do iatm=1,ncenter
	dist2=(a(iatm)%x-x)**2+(a(iatm)%y-y)**2+(a(iatm)%z-z)**2
	if (any(atmlist==iatm)) then !Atoms inside
		if (dist2<dist2minin) then
			dist2minin=dist2
			iminin=iatm
		end if
	else !Atoms outside
		if (dist2<dist2minext) then
			dist2minext=dist2
			iminext=iatm
		end if
	end if
end do
if (iminin==0.or.iminext==0) then !In this case there must be a bug in constructing the surface
    surfana_norm=0
    return
end if
di=dsqrt(dist2minin)
de=dsqrt(dist2minext)
rvdwin=vdwr(a(iminin)%index)
rvdwext=vdwr(a(iminext)%index)
surfana_norm=( (di-rvdwin)/rvdwin+(de-rvdwext)/rvdwext )*b2a
end function




!!-------- PAEM, potential acting on one electron in a molecule, defined by Zhongzhi Yang in JCC,35,965(2014)
!If itype=1, evaluate the XC potential based on pair density with Muller approximation, in this case GTFint will be generated first and used throughout
!If itype=2, it will be equivalent to DFT XC potential, in this case libreta will be used to compute ESP
real*8 function PAEM(x,y,z,itype)
integer itype
real*8 x,y,z,wfnval(nmo),GTFint(nprims,nprims)
if (itype==1) then !Based on Muller approximation form
	call genGTFattmat(x,y,z,GTFint) !GTFint will be used in the following parts
	!Evaluate electron contribution to ESP
	rhopot=0
	do imo=1,nmo
		do iprim=1,nprims
			do jprim=1,nprims
				rhopot=rhopot+MOocc(imo)*CO(imo,iprim)*CO(imo,jprim)*GTFint(iprim,jprim)
			end do
		end do
	end do
	!Evaluate XC potential
	call orbderv(1,1,nmo,x,y,z,wfnval)
	rho=sum(MOocc(1:nmo)*wfnval(1:nmo)**2)
	xcpot=0
	if (wfntype==0.or.wfntype==3) then !Closed-shell
		do imo=1,nmo
			if (MOocc(imo)==0D0) cycle
			do jmo=1,nmo
				if (MOocc(jmo)==0D0) cycle
				tmpval=dsqrt(MOocc(imo)*MOocc(jmo))*wfnval(imo)*wfnval(jmo)
				do iprim=1,nprims
					do jprim=1,nprims
						xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
					end do
				end do
			end do
		end do
	else if (wfntype==1.or.wfntype==4) then !Unrestricted open-shell
		do ialphaend=nmo,1,-1 !Find the ending index of alpha MO
			if (MOtype(ialphaend)==1) exit
		end do
		do imo=1,ialphaend !Alpha part
			do jmo=1,ialphaend
				tmpval=dsqrt(MOocc(imo)*MOocc(jmo))*wfnval(imo)*wfnval(jmo)
				do iprim=1,nprims
					do jprim=1,nprims
						xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
					end do
				end do
			end do
		end do
		do imo=ialphaend+1,nmo !Beta part
			do jmo=ialphaend+1,nmo
				tmpval=dsqrt(MOocc(imo)*MOocc(jmo))*wfnval(imo)*wfnval(jmo)
				do iprim=1,nprims
					do jprim=1,nprims
						xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
					end do
				end do
			end do
		end do
	else if (wfntype==2) then !Restricted open-shell
		do imo=1,nmo !Alpha part
			if (MOocc(imo)==0) cycle
			do jmo=1,nmo
				if (MOocc(jmo)==0) cycle
				tmpval=wfnval(imo)*wfnval(jmo) !Every occupied ROHF MOs contributes one alpha electron
				do iprim=1,nprims
					do jprim=1,nprims
						xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
					end do
				end do
			end do
		end do
		do imo=1,nmo !Beta part
			if (MOocc(imo)/=2D0) cycle
			do jmo=1,nmo
				if (MOocc(jmo)/=2D0) cycle
				tmpval=wfnval(imo)*wfnval(jmo) !Every doubly occupied ROHF MOs contributes one beta electron
				do iprim=1,nprims
					do jprim=1,nprims
						xcpot=xcpot-tmpval*CO(imo,iprim)*CO(jmo,jprim)*GTFint(iprim,jprim)
					end do
				end do
			end do
		end do
	end if
	xcpot=xcpot/rho
	PAEM=-nucesp(x,y,z)+rhopot+xcpot
else if (itype==2) then !Directly using DFT XC potential
	PAEM=-totesp(x,y,z)+DFTxcpot(x,y,z,0)
end if

end function




!!------- Angle between the eigenvectors of rho and the plane defined by option 4 of main function 1000
! The plane is represented by global variables pleA,pleB,pleC,pleD
! ivec=1/2/3 means the eigenvector corresponding to the first/second/third highest eigenvalue is calculated
real*8 function Ang_rhoeigvec_ple(x,y,z,ivec)
use util
integer ivec
real*8 x,y,z,eigvecmat(3,3),eigval(3),eigvaltmp(3),elehess(3,3),elegrad(3)
call calchessmat_dens(2,x,y,z,elerho,elegrad,elehess)
call diagmat(elehess,eigvecmat,eigval,100,1D-10)
eigvaltmp=eigval
call sort(eigvaltmp) !From small to large
call invarr(eigvaltmp) !1/2/3=large to small
do i=1,3
	if (eigval(i)==eigvaltmp(ivec)) exit
end do
Ang_rhoeigvec_ple=vecang(eigvecmat(1,i),eigvecmat(2,i),eigvecmat(3,i),pleA,pleB,pleC)
end function


!!------ Local electron correlation function (DOI: 10.1021/acs.jctc.7b00293)
!itype=1: Local total electron correlation function
!itype=2: Local dynamic electron correlation function
!itype=3: Local nondynamic electron correlation function
real*8 function localcorr(x,y,z,itype)
integer itype
real*8 x,y,z,wfnval(nmo),occ(nmo)
call orbderv(1,1,nmo,x,y,z,wfnval)
localcorr=0D0
if (wfntype==3) then
	occ=MOocc/2
	where(occ>1) occ=1 !Remove unphysical larger than unity occupation number
	where(occ<0) occ=0 !Remove unphysical negative occupation number
	if (itype==1) then
		do i=1,nmo
			localcorr=localcorr+ dsqrt(occ(i)*(1-occ(i)))*wfnval(i)**2
		end do
		localcorr=localcorr/4
	else if (itype==2) then
		do i=1,nmo
			localcorr=localcorr+ ( dsqrt(occ(i)*(1-occ(i))) - 2*occ(i)*(1-occ(i)) ) *wfnval(i)**2
		end do
		localcorr=localcorr/4
	else if (itype==3) then
		do i=1,nmo
			localcorr=localcorr+ occ(i)*(1-occ(i))*wfnval(i)**2
		end do
		localcorr=localcorr/2
	end if
	localcorr=localcorr*2
else if (wfntype==4) then
	occ=MOocc
	where(occ>1) occ=1
	where(occ<0) occ=0
	if (itype==1) then
		do i=1,nmo
			localcorr=localcorr+ dsqrt(occ(i)*(1-occ(i)))*wfnval(i)**2
		end do
		localcorr=localcorr/4
	else if (itype==2) then
		do i=1,nmo
			localcorr=localcorr+ ( dsqrt(occ(i)*(1-occ(i))) - 2*occ(i)*(1-occ(i)) ) *wfnval(i)**2
		end do
		localcorr=localcorr/4
	else if (itype==3) then
		do i=1,nmo
			localcorr=localcorr+ occ(i)*(1-occ(i))*wfnval(i)**2
		end do
		localcorr=localcorr/2
	end if
end if
end function




!!------ For visually examine functions used in DFRT2.0 project
real*8 function funcvalLSB(itype,x,y,z)
integer itype
integer,parameter :: nfunc=8
real*8 x,y,z,valarr(nfunc),rho,gradrho(3)
iexpcutoffold=expcutoff
expcutoff=1
call valaryyLSB(x,y,z,valarr,rho,gradrho)
funcvalLSB=valarr(itype)
expcutoff=iexpcutoffold
end function




!!------ Orbital-weighted Fukui function or dual descriptor
!itype=1: f+   =2: f-   =3: f0   =4: Dual descriptor
!The delta parameter is defined by global variable "orbwei_delta"
real*8 function orbwei_Fukui(itype,x,y,z)
integer itype
real*8 x,y,z,wfnval(nmo),expterm(nmo)

idxHOMO=nint(naelec) !It is assumed that the occupation has not been altered
idxLUMO=idxHOMO+1
chempot=(MOene(idxHOMO)+MOene(idxLUMO))/2

cutthres=0.001D0 !Orbitals with weight less than 0.1% will be ignored
if (itype==1.or.itype==3.or.itype==4) then !f+
    do imo=idxLUMO,nmo
        expterm(imo)=exp( -((chempot-MOene(imo))/orbwei_delta)**2 )
        if (expterm(imo)<cutthres) then
            icut=imo-1
            exit
        end if
    end do
    denomin=sum(expterm(idxLUMO:icut))
    call orbderv(1,idxLUMO,icut,x,y,z,wfnval)
    orbwei_Fpos=0
    do imo=idxLUMO,icut
        wei=expterm(imo)/denomin
	    orbwei_Fpos=orbwei_Fpos+wei*wfnval(imo)**2
    end do
end if
if (itype==2.or.itype==3.or.itype==4) then !f-
    do imo=idxHOMO,1,-1
        expterm(imo)=exp( -((chempot-MOene(imo))/orbwei_delta)**2 )
        if (expterm(imo)<cutthres) then
            icut=imo+1
            exit
        end if
    end do
    denomin=sum(expterm(icut:idxHOMO))
    call orbderv(1,icut,idxHOMO,x,y,z,wfnval)
    orbwei_Fneg=0
    do imo=icut,idxHOMO
        wei=expterm(imo)/denomin
	    orbwei_Fneg=orbwei_Fneg+wei*wfnval(imo)**2
    end do
end if

if (itype==1) then !f+
    orbwei_Fukui=orbwei_Fpos
else if (itype==2) then !f-
    orbwei_Fukui=orbwei_Fneg
else if (itype==3) then !f0
    orbwei_Fukui=(orbwei_Fpos+orbwei_Fneg)/2
else if (itype==4) then !DD
    orbwei_Fukui=orbwei_Fpos-orbwei_Fneg
end if
end function




!!----- For three-point numerical fit to evaluate D(r)
!Contributed by Arshad Mehmood
subroutine three_point_interpolation(n,x,y,xmax,ymax)
integer, intent(in) :: n
real*8, intent(in) :: x(max_edr_exponents),y(max_edr_exponents)
real*8, intent(out) :: xmax,ymax
integer i , imax
real*8  x1,x2,x3, y1,y2,y3, a,b
100 format ('XXX ',3F9.5)
ymax = -1D0
imax = -1 
do i=1,n
  if(y(i) .gt. ymax) then 
     ymax = y(i)
     imax = i
     endif 
  end do 
if(imax<1 .or. imax>n) then
	write(*,*) "Error: Bad imax"
	call EXIT()
end if
if(imax .eq. 1 .or. imax.eq.n) then
   xmax = x(imax)**(-0.5d0)
   return 
endif
x1 = x(imax-1)**(-0.5d0)
x2 = x(imax  )**(-0.5d0)
x3 = x(imax+1)**(-0.5d0)
y1 = y(imax-1)
y2 = y(imax  )
y3 = y(imax+1)
a = ( (y3-y2)/(x3-x2) -(y2-y1)/(x2-x1) )/(x3-x1)
b = ( (y3-y2)/(x3-x2)*(x2-x1) + (y2-y1)/(x2-x1)*(x3-x2) )&
      /(x3-x1)
xmax = x2 - b/(2d0*a)
ymax = y2 - b**(2d0)/(4d0*a)
end subroutine 

!!----- Function to calculate EDR(r,d)
!Contributed by Arshad Mehmood
!J. Chem. Phys. 141, 144104(2014), J. Chem. Theory Comput. 12, 79(2016), Angew. Chem. Int. Ed. 56, 6878(2017)
real*8 function edr(x,y,z)
real*8 :: ed(max_edr_exponents),edrval(max_edr_exponents)
nedr=1
ed(1)=dedr**(-2D0)
call EDRcal(1,x,y,z,nedr,ed,edrval,edrdmaxval)
edr=edrval(1)
end function

!!----- Function to calculate D(r)
!Contributed by Arshad Mehmood
!J. Chem. Theory Comput. 12, 3185(2016), Phys. Chem. Chem. Phys. 17, 18305(2015)
real*8 function edrdmax(x,y,z)	
real*8 :: ed(max_edr_exponents),edrdmaxval,edrval(max_edr_exponents)
real*8 :: edrexponent
integer iedr
edrexponent = edrastart  
do iedr=1,nedr
   ed(iedr)=edrexponent
   edrexponent=edrexponent/edrainc
end do
call EDRcal(2,x,y,z,nedr,ed,edrval,edrdmaxval)
edrdmax=edrdmaxval
end function

!!----- Working routine used to evaluate EDR(r;d) and D(r)
!Contributed by Arshad Mehmood
subroutine EDRcal(runtype,x,y,z,nedr,ed,edrval,edrdmaxval) 
real*8, intent(in) :: x,y,z,ed(max_edr_exponents)
integer, intent(in) :: nedr
real*8, intent(out):: edrdmaxval,edrval(max_edr_exponents) 
real*8 :: rho,dmaxdummy  
real*8 :: psi(nmo),AMUVal(max_edr_exponents),Bint(nmo,max_edr_exponents) 
real*8 :: xamu(3,max_edr_exponents),amu0(max_edr_exponents)  
integer :: j,ixyz,i,iedr,runtype
edrval = 0D0
if (runtype==2) then
    edrdmaxval = 0D0
end if

rho=fdens(x,y,z)

if(rho.gt.1D-10) then 
	psi = 0d0
	Bint=0d0 
	do j=1,nprims    
		ix=type2ix(b(j)%type)
		iy=type2iy(b(j)%type)
		iz=type2iz(b(j)%type)
		ep=b(j)%exp	
		sftx=x-a(b(j)%center)%x
		sfty=y-a(b(j)%center)%y
		sftz=z-a(b(j)%center)%z
		sftx2=sftx*sftx
		sfty2=sfty*sfty
		sftz2=sftz*sftz
		rr=sftx2+sfty2+sftz2	 
		expterm=0.0
		amu0 = 0d0
		if (expcutoff>0.or.-ep*rr>expcutoff) then
			expterm=exp(-ep*rr) 
			do iedr=1,nedr
               amu0(iedr)=(2d0*ed(iedr)/pi)**(3d0/4d0) *(pi/(ep+ed(iedr)))**(3d0/2d0)&
               * exp(-ep*ed(iedr)/(ep+ed(iedr))*rr)
			end do
		end if 
 
		if (expterm==0D0) cycle
		 
		GTFval=sftx**ix *sfty**iy *sftz**iz *expterm
		do iedr=1,nedr
			do ixyz=1,3
				ival=ix
				sftval=sftx
				if(ixyz.eq.2) then
					ival=iy
					sftval=sfty
				else if(ixyz.eq.3) then
					ival=iz
					sftval=sftz
				end if
				If(ival.eq.0) then
					xamu(ixyz,iedr)=1d0 
				else if(ival.eq.1) then 
					xamu(ixyz,iedr)=sftval*ed(iedr)/(ed(iedr)+ep)
				else If(ival.eq.2) then
					xamu(ixyz,iedr)=(sftval*ed(iedr)/(ed(iedr)+ep))**2d0 + 1d0/(2d0*(ed(iedr)+ep))
				else If(ival.eq.3) then
					xamu(ixyz,iedr)=(sftval*ed(iedr)/(ed(iedr)+ep))**3d0 + sftval*3d0*ed(iedr)/(2d0*(ed(iedr)+ep)**2d0)
				else If(ival.eq.4) then
					xamu(ixyz,iedr)=sftval**4d0* (ed(iedr)/(ed(iedr)+ep))**4d0 + sftval**2d0 *3d0*ed(iedr)**2d0/(ed(iedr)+ep)**3d0 &
						+ 3d0/(4d0*(ed(iedr)+ep)**2d0)
				else If(ival.eq.5) then
					xamu(ixyz,iedr)=sftval**5d0*  ed(iedr)**5d0/(ed(iedr)+ep)**5d0 + sftval**3d0* 5d0*ed(iedr)**3d0/(ed(iedr)+ep)**4d0 &
						+ sftval *15d0*ed(iedr)/(4d0*(ed(iedr)+ep)**3d0)
				else 
					write(*,*) "Angular momentum out of range"
					Call EXIT()
				end if
			end do
		end do
		do iedr=1,nedr
			AMUVal(iedr)=amu0(iedr)*xamu(1,iedr)*xamu(2,iedr)*xamu(3,iedr)
		end do
		
		do i=1,nmo
			if (nint(MOocc(i)).GE.1D0) then	
				psi(i)=psi(i)+CO(i,j)*GTFval
			end if
		end do

		do iedr=1,nedr
			do i=1,nmo
				if (nint(MOocc(i)).GE.1D0) then	
					Bint(i,iedr)=Bint(i,iedr)+CO(i,j)*AMUVal(iedr)
				end if
		    end do
		end do
	end do  

	edrval = 0d0
	do i=1,nmo
		do iedr=1,nedr
		    edrval(iedr)=edrval(iedr)+psi(i)*Bint(i,iedr)
		end do
    end do
	if (runtype==2) then
		call three_point_interpolation(nedr,ed,edrval,edmax,dmaxdummy)
		edrdmaxval=edmax
	else if (runtype==1) then
    	do iedr=1,nedr
			edrval(iedr)=edrval(iedr)*rho**(-0.5D0)
		end do
	else	
		write(*,*) "EDRcal runtype out of range"
		call EXIT()
	end if
end if  
end subroutine




!----------- Local Hartree-Fock exchange energy
!genPprim must have been called so that density matrix is available
real*8 function locHFexc(x,y,z)
implicit real*8 (a-h,o-z)
real*8 x,y,z,Vprim(nprims,nprims),GTFarr(nprims),Earr(nprims)

call calcGTFval(x,y,z,GTFarr)
call getVmatprim(x,y,z,Vprim)

if (wfntype==0.or.wfntype==3) then !Closed-shell
	do i=1,nprims
		Earr(i)=sum(Ptot_prim(:,i)*GTFarr(:))
	end do
	locHFexc=0
	do ib=1,nprims
		do id=1,nprims
			locHFexc=locHFexc+Earr(ib)*Earr(id)*Vprim(ib,id)
		end do
	end do
	locHFexc=-locHFexc/4D0
else !Open-shell
	locHFexc=0
    !Alpha part
	do i=1,nprims
		Earr(i)=sum(Palpha_prim(:,i)*GTFarr(:))
	end do
	do ib=1,nprims
		do id=1,nprims
			locHFexc=locHFexc+Earr(ib)*Earr(id)*Vprim(ib,id)
		end do
	end do
    !Beta part
	do i=1,nprims
		Earr(i)=sum(Pbeta_prim(:,i)*GTFarr(:))
	end do
	do ib=1,nprims
		do id=1,nprims
			locHFexc=locHFexc+Earr(ib)*Earr(id)*Vprim(ib,id)
		end do
	end do
	locHFexc=-locHFexc/2D0
end if
end function




!!---------- Calculate kinetic energy density (KED)
!The form of KED is determined by iform, the index definition is the same as iKEDsel, which is described in the manual
real*8 function KED(x,y,z,iform)
real*8 KEDarr(nKEDmax),x,y,z
integer iform
call KEDall(x,y,z,KEDarr)
KED=KEDarr(iform)
end function

!!---------- Get difference between the KED selected by iform and specific KED
!idiff=1: selected KED - Weisacker KED, invoked by iuserfunc=1201
!idiff=2: selected KED - Lagrangian KED, invoked by iuserfunc=1202
!idiff=3: |selected KED - Lagrangian KED|, invoked by iuserfunc=1203
real*8 function KEDdiff(x,y,z,iform,idiff)
real*8 KEDarr(nKEDmax),x,y,z
integer iform,idiff
call KEDall(x,y,z,KEDarr)
selKED=KEDarr(iform)
if (idiff==1) KEDdiff=selKED-KEDarr(4)
if (idiff==2) KEDdiff=selKED-KEDarr(2)
if (idiff==3) KEDdiff=abs(selKED-KEDarr(2))
end function

!!---------- Return all kinds of kinetic energy density (KED) as an array
!Most KEDs have general form, and need rho and reduced density gradient (rg)
subroutine KEDall(x,y,z,KEDarr)
implicit real*8 (a-h,o-z)
real*8 KEDarr(nKEDmax),x,y,z,wfnval(nmo),wfnderv(3,nmo),wfnhess(3,3,nmo),EDFgrad(3),EDFhess(3,3)
real*8 KED_TF(0:2),rho(0:2),lapltot(0:2),rg(0:2),rgs(0:2) !0/1/2 = total electron/alpha electron/beta electron
real*8 gradrho(3,0:2),laplrho(3,0:2) !First index 1/2/3=x/y/z, last index 0/1/2 = total electron/alpha electron/beta electron
real*8 :: C_TF=4.557799872D0 !Fermi constant for spin polarized = (3/10)*(6*Pi^2)**(2/3)
real*8 :: C_TF_cls=2.871234D0 !Fermi constant for closed-shell = (3/10)*(3*Pi^2)**(2/3)

ider=1 !Only need 1st derivative of GTFs
if (iKEDsel==1.or.iKEDsel==24.or.uservar/=0) ider=2 !Need second order derivative of GTFs

if (ider==2) then
	call orbderv(3,1,nmo,x,y,z,wfnval,wfnderv,wfnhess)
	if (nEDFprims/=0) call EDFrho(3,x,y,z,EDFdens,EDFgrad,EDFhess)
else
	call orbderv(2,1,nmo,x,y,z,wfnval,wfnderv)
	if (nEDFprims/=0) call EDFrho(2,x,y,z,EDFdens,EDFgrad)
end if

!Get density and its gradient vector
rho=0
gradrho=0
do imo=1,nmo
	if (MOocc(imo)==0) cycle
	if (MOtype(imo)==0) then
		occtmp=MOocc(imo)/2
		rho(1)=rho(1)+occtmp*wfnval(imo)**2
		rho(2)=rho(2)+occtmp*wfnval(imo)**2
		gradrho(:,1)=gradrho(:,1)+occtmp*wfnval(imo)*wfnderv(:,imo)
		gradrho(:,2)=gradrho(:,2)+occtmp*wfnval(imo)*wfnderv(:,imo)
	else if (MOtype(imo)==1) then
		rho(1)=rho(1)+MOocc(imo)*wfnval(imo)**2
		gradrho(:,1)=gradrho(:,1)+MOocc(imo)*wfnval(imo)*wfnderv(:,imo)
	else
		rho(2)=rho(2)+MOocc(imo)*wfnval(imo)**2
		gradrho(:,2)=gradrho(:,2)+MOocc(imo)*wfnval(imo)*wfnderv(:,imo)
	end if
end do
gradrho=2*gradrho
if (nEDFprims/=0) then
	rho(1)=rho(1)+EDFdens/2
	rho(2)=rho(2)+EDFdens/2
	gradrho(:,1)=gradrho(:,1)+EDFgrad/2
	gradrho(:,2)=gradrho(:,2)+EDFgrad/2
end if
rho(0)=rho(1)+rho(2)
gradrho(:,0)=gradrho(:,1)+gradrho(:,2)

!Get Laplacian
if (ider==2) then
	laplrho=0
	do imo=1,nmo
		if (MOocc(imo)==0) cycle
		if (MOtype(imo)==0) then
			do idir=1,3
				tmp=MOocc(imo)/2*( wfnderv(idir,imo)**2 + wfnval(imo)*wfnhess(idir,idir,imo) )
				laplrho(idir,1)=laplrho(idir,1)+tmp
				laplrho(idir,2)=laplrho(idir,2)+tmp
			end do
		else if (MOtype(imo)==1) then
			do idir=1,3
				laplrho(idir,1)=laplrho(idir,1)+MOocc(imo)*( wfnderv(idir,imo)**2 + wfnval(imo)*wfnhess(idir,idir,imo) )
			end do
		else
			do idir=1,3
				laplrho(idir,2)=laplrho(idir,2)+MOocc(imo)*( wfnderv(idir,imo)**2 + wfnval(imo)*wfnhess(idir,idir,imo) )
			end do
		end if
	end do
	laplrho=2*laplrho
	if (nEDFprims/=0) then
		do idir=1,3
			laplrho(idir,1)=laplrho(idir,1)+EDFhess(idir,idir)/2
			laplrho(idir,2)=laplrho(idir,2)+EDFhess(idir,idir)/2
		end do
	end if
	laplrho(:,0)=laplrho(:,1)+laplrho(:,2)
	do i=0,2 !laptot(0/1/2)=laplacian of rho of total/alpha/beta electrons
		lapltot(i)=sum(laplrho(:,i))
	end do
end if

!Get reduced density gradient (with/without prefactor)
if (rho(0)==0D0) then
	rg=0
	rgs=0
else !rg(1),(2),(3) correspond to RDG evaluated based on total, alpha and beta density
	do i=0,2
		rg(i)=dsqrt(sum(gradrho(:,i)**2))/rho(i)**(4D0/3D0)
	end do
	rgs(:)=0.161620459673995D0*rg(:) !0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
end if

KEDarr=0

!Hamltionian KED
if (ider==2) then
	hamx=sum( MOocc(1:nmo)*wfnhess(1,1,1:nmo)*wfnval(1:nmo) )
	hamy=sum( MOocc(1:nmo)*wfnhess(2,2,1:nmo)*wfnval(1:nmo) )
	hamz=sum( MOocc(1:nmo)*wfnhess(3,3,1:nmo)*wfnval(1:nmo) )
	KEDarr(1)=-(hamx+hamy+hamz)/2
end if

!Lagrangian KED
do imo=1,nmo
	KEDarr(2)=KEDarr(2)+MOocc(imo)*sum(wfnderv(:,imo)**2)
end do
KEDarr(2)=KEDarr(2)/2

!TF
do i=1,2
	KED_TF(i)=C_TF*rho(i)**(5D0/3D0)
end do
KED_TF(0)=KED_TF(1)+KED_TF(2)
KEDarr(3)=KED_TF(0)

!Weizsacker
do i=1,2
	KEDarr(4)=KEDarr(4)+sum(gradrho(:,i)**2)/8/rho(i)
end do

!Most below KEDs are written as KED_TF * F_enhancement general form

!GEA2
do i=1,2
	Fenh=1+1D0/(72*C_TF)*rg(i)**2
	KEDarr(5)=KEDarr(5)+ KED_TF(i)*Fenh
end do

!TF5W
do i=1,2
	Fenh=1+1D0/(40*C_TF)*rg(i)**2
	KEDarr(6)=KEDarr(6)+ KED_TF(i)*Fenh
end do

!TFvW
do i=1,2
	Fenh=1+1D0/(8*C_TF)*rg(i)**2
	KEDarr(7)=KEDarr(7)+ KED_TF(i)*Fenh
end do

!TF9W
do i=1,2
	Fenh=1+1.067D0/(72*C_TF)*rg(i)**2
	KEDarr(8)=KEDarr(8)+ KED_TF(i)*Fenh
end do

!TF-N
do i=1,2
	Fenh=1+0.313D0/nelec**(1D0/3D0)-0.187/nelec**(2D0/3D0)
	KEDarr(9)=KEDarr(9)+ KED_TF(i)*Fenh
end do

!Pearson
do i=1,2
	rgs2=rgs(i)/2**(1D0/3D0) !rgs2=[rg/(2*(3*pi**2)**(1/3))] /2**(1D0/3D0) = rg/(2*(6*pi**2)**(1/3))
	t2=sum(gradrho(:,i)**2)/rho(i)/72
	KEDarr(10)=KEDarr(10)+ KED_TF(i)+t2/(1+rgs2**6)
end do
!PS: Closed form
! t2=sum(gradrho(:,0)**2)/rho(0)/72
! KEDarr(10)=KED_TF(0)+t2/(1+rgs(0)**6)

!DK Pade
a1=0.95D0;a2=14.28111D0;a3=-19.57962D0;b1=-0.05D0;b2=9.99802D0;b3=2.96085D0
do i=1,2
	xt=rg(i)**2/(72*C_TF)
	tmp1= 9*b3*xt**4 + a3*xt**3 + a2*xt**2 + a1*xt + 1
	tmp2= b3*xt**3 + b2*xt**2 + b1*xt + 1
	Fenh=tmp1/tmp2
	KEDarr(11)=KEDarr(11)+ KED_TF(i)*Fenh
end do

!LLP
do i=1,2
	Fenh=1+0.0044188D0*rg(i)**2/(1+0.0253D0*rg(i)*asinh(rg(i)))
	KEDarr(12)=KEDarr(12)+ KED_TF(i)*Fenh
end do

!OL1
do i=1,2
	Fenh=1 + 1D0/(72*C_TF)*rg(i)**2 + 0.00187D0*rg(i)
	KEDarr(13)=KEDarr(13)+ KED_TF(i)*Fenh
end do

!OL2
do i=1,2
	Fenh=1 + 1D0/(72*C_TF)*rg(i)**2 + 0.0245D0*rg(i)/(1+2**(5D0/3D0)*rg(i))
	KEDarr(14)=KEDarr(14)+ KED_TF(i)*Fenh
end do

!Thak
do i=1,2
	Fenh=1 + 0.0055D0*rg(i)**2/(1+0.0253D0*rg(i)*asinh(rg(i))) - 0.072D0*rg(i)/(1+2**(5D0/3D0)*rg(i))
	KEDarr(15)=KEDarr(15)+ KED_TF(i)*Fenh
end do

!B86A
do i=1,2
	Fenh=1 + 0.0039D0*rg(i)**2/(1+0.004D0*rg(i)**2)
	KEDarr(16)=KEDarr(16)+ KED_TF(i)*Fenh
end do

!B86B
do i=1,2
	Fenh=1 + 0.00403D0*rg(i)**2/(1+0.007D0*rg(i)**2)**(4D0/5D0)
	KEDarr(17)=KEDarr(17)+ KED_TF(i)*Fenh
end do

!DK87. The result differs from Table 1 of LG94 paper, I don't know why. But my current implementation must be correct
do i=1,2
	Fenh=1 + 0.00132327D0*rg(i)**2 * (1D0+0.861504D0*rg(i))/(1D0+0.044286D0*rg(i)**2) !The coefficient is taken from DK87 original paper
! 	Fenh=1 + 7D0/(324D0*(18D0*pi**4D0)**(1D0/3D0))* rg2**2 * (1D0+0.861504D0*rg2)/(1D0+0.044286D0*rg2**2) !Introduced in LG94 paper, seemingly incorrect
	KEDarr(18)=KEDarr(18)+ KED_TF(i)*Fenh
end do

!PW86
do i=1,2
	rgs2=rgs(i)/2**(1D0/3D0)
	Fenh=(1 + 1.296D0*rgs2**2 + 14*rgs2**4 + 0.2D0*rgs2**6)**(1D0/15D0)
	KEDarr(19)=KEDarr(19)+ KED_TF(i)*Fenh
end do
!PS: Equivalent closed form
! Fenh=(1 + 1.296D0*rgs(0)**2 + 14*rgs(0)**4 + 0.2D0*rgs(0)**6)**(1D0/15D0)
! KEDarr(19)=KED_TF(0)*Fenh

!PW91
a1=0.19645D0;a2=0.2743D0;a3=0.1508D0;a4=0.004D0;btmp=7.7956D0
do i=1,2
	rgs2=rgs(i)/2**(1D0/3D0)
	tmp1=1+a1*rgs2*asinh(btmp*rgs2)+(a2-a3*exp(-100*rgs2**2))*rgs2**2
	tmp2=1+a1*rgs2*asinh(btmp*rgs2)+a4*rgs2**4
	Fenh=tmp1/tmp2
	KEDarr(20)=KEDarr(20)+ KED_TF(i)*Fenh
end do
!PS: Equivalent closed form
! tmp1=1+a1*rgs(0)*asinh(btmp*rgs(0))+(a2-a3*exp(-100*rgs(0)**2))*rgs(0)**2
! tmp2=1+a1*rgs(0)*asinh(btmp*rgs(0))+a4*rgs(0)**4
! Fenh=tmp1/tmp2
! KEDarr(20)=KED_TF(0)*Fenh

!LG94. I cannot reproduce the data in Table I original paper, I really do not known why. My implementation should be correct
btmp=0.024974D0
a2=(1E-8+0.1234D0)/btmp; a4=29.790D0; a6=22.417D0; a8=12.119D0; a10=1570.1D0; a12=55.944D0
do i=1,2
	rgs2=rgs(i)/2**(1D0/3D0)
	tmp1=1+a2*rgs2**2+a4*rgs2**4+a6*rgs2*6+a8*rgs2**8+a10*rgs2**10+a12*rgs2**12
	tmp2=1+1E-8*rgs2**2
	Fenh=tmp1**btmp/tmp2
	KEDarr(21)=KEDarr(21)+ KED_TF(i)*Fenh
end do
!PS: Closed form
! tmp1=1+a2*rgs(0)**2+a4*rgs(0)**4+a6*rgs(0)*6+a8*rgs(0)**8+a10*rgs(0)**10+a12*rgs(0)**12
! tmp2=1+1E-8*rgs(0)**2
! Fenh=tmp1**btmp/tmp2
! KEDarr(21)=KED_TF(0)*Fenh

!ABSP
do i=1,2
	Fenh=rg(i)**2/(8*C_TF) + 1-1.412D0/nelec**(1D0/3D0)
	KEDarr(22)=KEDarr(22)+ KED_TF(i)*Fenh
end do

!GR
do i=1,2
	Fenh=rg(i)**2/(8*C_TF) + (1-2D0/nelec) * (1 - 1.303D0/nelec**(1D0/3D0) + 0.029D0/nelec**(2D0/3D0))
	KEDarr(23)=KEDarr(23)+ KED_TF(i)*Fenh
end do

!GEA4. see LG94 KED paper Eq.3 for closed-shell form and DK Pade paper for open-shell form (used here)
if (ider==2) then
	do i=1,2
		gr2=sum(gradrho(:,i)**2)
		GEA2=KED_TF(i)+gr2/72/rho(i)
		pretmp=(6*pi**2)**(-2D0/3D0)/540 * rho(i)**(1D0/3D0)
		corr4= pretmp*( (lapltot(i)/rho(i))**2 -9D0/8D0*lapltot(i)/rho(i)*gr2/rho(i)**2 + gr2**2/rho(i)**4/3D0 )
		KEDarr(24)=KEDarr(24)+GEA2+corr4
	end do
!Closed form
! 	gr2=sum(gradrho(:,0)**2)
! 	GEA2=KED_TF(0)+gr2/72/rho(0)
! 	pretmp=(3*pi**2)**(-2D0/3D0)/540 * rho(0)**(1D0/3D0)
!  	corr4= pretmp*( (lapltot(0)/rho(0))**2 -9D0/8D0*lapltot(0)/rho(0)*gr2/rho(0)**2 + gr2**2/rho(0)**4/3D0 )
! 	KEDarr(24)=KEDarr(24)+GEA2+corr4
end if

!Add laplacian term into all KEDs. When uservar=6, the KED of iKEDsel=5 (i.e. GEA2) just corresponds to the KED employed in Tsirelson type of ELF/LOL
if (uservar/=0) then
	KEDarr=KEDarr+lapltot(0)/uservar !Open and closed shell forms are equivalent
end if

end subroutine




!!---------- Potential of kinetic energy functional
!Using closed-shell form
real*8 function KEDpot(x,y,z)
implicit real*8 (a-h,o-z)
real*8 x,y,z,gradrho(3),hessrho(3,3)
real*8 :: C_TF=2.871234D0 !Fermi constant for closed-shell = (3/10)*(3*Pi^2)**(2/3)

call calchessmat_dens(2,x,y,z,rho,gradrho,hessrho)
gradnorm=dsqrt(sum(gradrho(:)**2))
grad2=gradnorm**2
rholapl=hessrho(1,1)+hessrho(2,2)+hessrho(3,3)
rg=gradnorm/rho**(4D0/3D0) !Reduced density gradient
rgs=0.161620459673995D0*rg !RDG with factor. 0.161620459673995D0=1/(2*(3*pi**2)**(1/3))
TFKED=C_TF*rho**(5D0/3D0) !Thomas-Fermi KED
TFpot=5D0/3D0*C_TF*rho**(2D0/3D0) !Potential of Thomas-Fermi KED

if (iKEDsel==3) then !Thomas-Fermi
    KEDpot=TFpot
else if (iKEDsel==5) then !Second order GEA
    F_GEA2=1+grad2/rho**(8D0/3D0)/72D0/C_TF
    tmp1=TFpot*F_GEA2+TFKED*(-grad2/rho**(11D0/3D0)/27D0/C_TF)
    tmp2=-(rholapl/rho-grad2/rho**2)/36D0
    KEDpot=tmp1+tmp2
    !KEDpot=TFpot+grad2/rho**2/72D0-rholapl/rho/36D0 !equivalent to above
else if (iKEDsel==7) then !Thomas-Fermi + Weizsacker
    F_TFW=1+grad2/rho**(8D0/3D0)/8D0/C_TF
    tmp1=TFpot*F_TFW+TFKED*(-grad2/rho**(11D0/3D0)/3D0/C_TF)
    tmp2=-(rholapl/rho-grad2/rho**2)/4D0
    KEDpot=tmp1+tmp2
end if
end function



!!--------- Fractional occupation number weighted electron density (FOD) based on fractionally occupied orbitals
real*8 function FODfunc(x,y,z)
implicit real*8 (a-h,o-z)
real*8 x,y,z,wfnval(nmo)
FODfunc=0
call orbderv(1,1,nmo,x,y,z,wfnval)
if (wfntype==3) then !Closed shell
	idxHOMO=nint(nelec/2)
	do imo=1,idxHOMO
		FODfunc=FODfunc+(2-MOocc(imo))*wfnval(imo)**2
    end do
	do imo=idxHOMO+1,nmo
		FODfunc=FODfunc+MOocc(imo)*wfnval(imo)**2
    end do
else if (wfntype==4) then !Open-shell
	do itmp=nmo,1,-1 !Find the last alpha MO
		if (MOtype(itmp)==1) exit
	end do
    !Alpha part
	do imo=1,nint(naelec)
		FODfunc=FODfunc+(1-MOocc(imo))*wfnval(imo)**2
    end do
	do imo=nint(naelec)+1,itmp
		FODfunc=FODfunc+MOocc(imo)*wfnval(imo)**2
    end do
    !Beta part
	do imo=itmp+1,itmp+nint(nbelec)
		FODfunc=FODfunc+(1-MOocc(imo))*wfnval(imo)**2
    end do
	do imo=itmp+nint(nbelec)+1,nmo
		FODfunc=FODfunc+MOocc(imo)*wfnval(imo)**2
    end do
end if
end function



!!--------- Calculate stress tensor
subroutine stress_tensor(x,y,z,mat)
real*8 x,y,z,mat(3,3),wfnval(nmo),grad(3,nmo),hess(3,3,nmo)

call orbderv(4,1,nmo,x,y,z,wfnval,grad,hess)
mat=0
do imo=1,nmo
    if (MOocc(imo)==0) cycle
    do i=1,3
        do j=i,3
            !tmpval = grad(i,imo)*grad(j,imo) + grad(j,imo)*grad(i,imo) - hess(i,j,imo)*wfnval(imo) - wfnval(imo)*hess(i,j,imo)
            tmpval = 2*grad(i,imo)*grad(j,imo) - 2*hess(i,j,imo)*wfnval(imo) !This is valid because only real wavefunction is processed
            mat(i,j)=mat(i,j)+MOocc(imo)*tmpval
        end do
    end do
end do
mat(2,1)=mat(1,2)
mat(3,1)=mat(1,3)
mat(3,2)=mat(2,3)
mat=-mat/4
end subroutine



!!--------- Calculate stress tensor stiffness
real*8 function stress_stiffness(x,y,z)
real*8 x,y,z,mat(3,3),eigval(3),eigvecmat(3,3)
call stress_tensor(x,y,z,mat)
call diagsymat(mat,eigvecmat,eigval,idiagok)
call sort(eigval) !Sort eigenvalues from low to high
stress_stiffness=abs(eigval(1))/abs(eigval(3))
end function



!!-------- Calculate stress tensor ellipticity. See http://sobereva.com/wfnbbs/viewtopic.php?pid=4667
real*8 function stress_ellipticity(x,y,z)
real*8 x,y,z,mat(3,3),eigval(3),eigvecmat(3,3)
call stress_tensor(x,y,z,mat)
call diagsymat(mat,eigvecmat,eigval,idiagok)
call sort(eigval) !Sort eigenvalues from low to high
stress_ellipticity=abs(eigval(1))/abs(eigval(2))-1
end function







!END of function module
!END of function module
!END of function module
!END of function module
end module
!END of function module
!END of function module
!END of function module
!END of function module




!! ----------- Show the list of all supported real space functions
!When user will select a function, should use "selfunc_interface" instead
subroutine funclist
use defvar
write(*,*) "            ----------- Available real space functions -----------"
if (allocated(b)) then
	write(*,*) "1 Electron density (rho)     2 Gradient norm of rho     3 Laplacian of rho"
	write(*,*) "4 Value of orbital wavefunction         44 Orbital probability density"
	if (ipolarpara==0) write(*,*) "5 Electron spin density"
	if (ipolarpara==1) write(*,*) "5 Spin polarization parameter function"
	write(*,*) "6 Hamiltonian kinetic energy density K(r)"
	write(*,*) "7 Lagrangian kinetic energy density G(r)"
	if (ifiletype==4) then
		write(*,*) "8 Electrostatic potential from atomic charges"
	else
		write(*,*) "8 Electrostatic potential from nuclear charges"
	end if
	if (ELFLOL_type==0) write(*,*) "9 Electron localization function (ELF)"
	if (ELFLOL_type==1) write(*,*) "9 Electron localization function (ELF) defined by Tsirelson" 
	if (ELFLOL_type==2) write(*,*) "9 Electron localization function (ELF) defined by Tian Lu" 
	if (ELFLOL_type==0) write(*,*) "10 Localized orbital locator (LOL)"
	if (ELFLOL_type==1) write(*,*) "10 Localized orbital locator (LOL) defined by Tsirelson" 
	if (ELFLOL_type==2) write(*,*) "10 Localized orbital locator (LOL) defined by Tian Lu" 
	write(*,*) "11 Local information entropy"
	write(*,*) "12 Total electrostatic potential (ESP)"
	write(*,*) "13 Reduced density gradient (RDG)       14 RDG with promolecular approximation"
	write(*,*) "15 Sign(lambda2)*rho      16 Sign(lambda2)*rho with promolecular approximation"
	!Fermi hole function only available to single-determinant wavefunction
	if (pairfunctype==1) write(*,"(a,3f10.5)") " 17 Correlation hole for alpha, ref. point:",refx,refy,refz
	if (pairfunctype==2) write(*,"(a,3f10.5)") " 17 Correlation hole for beta, ref. point:",refx,refy,refz
	if (pairfunctype==4) write(*,"(a,3f10.5)") " 17 Correlation factor for alpha, ref. point:",refx,refy,refz
	if (pairfunctype==5) write(*,"(a,3f10.5)") " 17 Correlation factor for beta, ref. point:",refx,refy,refz
	if (pairfunctype==7) write(*,"(a,3f10.5)") " 17 Exc.-corr. density for alpha, ref. point:",refx,refy,refz
	if (pairfunctype==8) write(*,"(a,3f10.5)") " 17 Exc.-corr. density for beta, ref. point:",refx,refy,refz
	if (pairfunctype==10) write(*,"(a,3f10.5)") " 17 Pair density for alpha, ref. point:",refx,refy,refz
	if (pairfunctype==11) write(*,"(a,3f10.5)") " 17 Pair density for beta, ref. point:",refx,refy,refz
	if (pairfunctype==12) write(*,"(a,3f10.5)") " 17 Pair density for all electrons, ref. point:",refx,refy,refz
	write(*,*) "18 Average local ionization energy (ALIE)"
	write(*,"(a,i2,a,3f10.5)") " 19 Source function, mode:",srcfuncmode,", ref. point:",refx,refy,refz
	write(*,*) "20 Electron delocal. range func. EDR(r;d)  21 Orbital overlap dist. func. D(r)"
	write(*,*) "22 Delta-g (promolecular approximation)    23 Delta-g (Hirshfeld partition)"
	write(*,"(a,a,a)") " 24 Interaction region indicator (IRI)    25 van der Waals potential (probe=",ind2name(ivdwprobe),')'
	write(*,"(a,i5,a)") " 100 User-defined function (iuserfunc=",iuserfunc,"), see Section 2.7 of manual"
else !No wavefunction information is available
	write(*,*) "1 Promolecular electron density "
	if (ifiletype==4) then
		write(*,*) "8 Electrostatic potential from atomic charges"
	else
		write(*,*) "8 Electrostatic potential from nuclear charges"
	end if
	write(*,*) "14 Reduced density gradient (RDG) with promolecular approximation"
	write(*,*) "16 Sign(lambda2)*rho with promolecular approximation"
	write(*,*) "22 Delta-g (promolecular approximation)"
	write(*,"(a,a,a)") " 25 van der Waals potential (probe=",ind2name(ivdwprobe),')'
	write(*,"(a,i5,a)") " 100 User-defined function (iuserfunc=",iuserfunc,")  See Section 2.7 of manual"
end if
end subroutine



!!------ Standard interface for selecting real space function
!Note that iorbsel is a global variable
!if itype=1, this routine will be used for normal case
!if itype=2, this routine will be used for plotting plane map, where two orbitals can be selected
subroutine selfunc_interface(itype,ifunc)
use defvar
integer ifunc,edrmaxpara,wrtnumedr
character c80tmp*80
real*8 wrtstart

call funclist
read(*,*) ifunc

if (ifunc==4.or.ifunc==44) then
	write(*,"(a,i6)") " Input orbital index, e.g. 28"
	if (allocated(CObasb)) write(*,"(a)") " Note: Positive index and negative index correspond to alpha orbital and beta orbital, respectively"
	if (wfntype==0.or.wfntype==1) call orblabsel_prompt
	if (ifunc==4.and.itype==2) write(*,"(a)") " If you want to plot contour map for two orbitals in the same map, input two indices, e.g. 8,10"
    read(*,"(a)") c80tmp
	if (index(c80tmp,',')/=0) then !Inputted two orbitals
		read(c80tmp,*) iorbsel,iorbsel2
    else
		if (index(c80tmp,'h')==0.and.index(c80tmp,'l')==0) then !Input orbital index
			read(c80tmp,*) iorbsel
			if (iorbsel<0.and.allocated(CObasb)) iorbsel=abs(iorbsel)+nbasis
			if (iorbsel<=0.or.iorbsel>nmo) then
				write(*,"(' Error: Orbital index should be in the range of 1 to',i6)") nmo
				write(*,*) "Orbital 1 is selected instead"
				iorbsel=1
			end if
		else !Input orbital label
			call orblabsel(c80tmp,iorbsel)
			if (iorbsel==0) then
				write(*,*) "Error: The orbital label you inputted is wrong! Please double check"
				write(*,*) "Orbital 1 is selected instead"
				iorbsel=1
			else
				write(*,"(a,i7)") " The orbital you selected is",iorbsel
			end if
		end if
    end if
else if (ifunc==20) then !Read length scale to evaluate EDR(r;d)
	write(*,*) "The EDR(r;d) computing code was contributed by Arshad Mehmood"
	write(*,"(a,/)") " References: J. Chem. Phys., 141, 144104 (2014); J. Chem. Theory Comput., 12, 79 (2016); Angew. Chem. Int. Ed., 56, 6878 (2017)"
	write(*,*) "Input length scale d (Bohr), e.g. 0.85"
	read(*,*) dedr
else if (ifunc==21) then 
	write(*,*) "The D(r) computing code was contributed by Arshad Mehmood"
	write(*,"(a,/)") " References: J. Chem. Theory Comput., 12, 3185 (2016); Phys. Chem. Chem. Phys., 17, 18305 (2015)"
	write(*,*) "1 Manually input total number, start and increment in EDR exponents"
	write(*,*) "2 Use default values, i.e. 20,2.50,1.50"
	read(*,*) edrmaxpara
	if (edrmaxpara==1) then  
		write(*,*) "Please input in order: exponents start increment, e.g. 20 2.5 1.5"
		write(*,*) "Note: Max. allowed exponents are 50 and min. allowed increment is 1.01"
		read(*,*) nedr,edrastart,edrainc
		if (nedr<1) then
			write(*,*) "Error: Bad Number of EDR exponents. Should be between 1 to 50"
			write(*,*) "Press ENTER button to exit"
			read(*,*)
			stop
		else if (nedr>50) then
			write(*,*) "Error: Bad Number of EDR exponents. Should be between 1 to 50"
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
		edrastart=2.5D0
		edrainc=1.5D0
	end if
	write(*,*) "The following EDR exponents will be used in calculation:"
	wrtstart=edrastart
	do wrtnumedr=1,nedr
		wrtexpo(wrtnumedr)=wrtstart
		wrtstart=wrtstart/edrainc
		write(*,"(E13.5)") wrtexpo(wrtnumedr) 
	end do
	write(*,*)
end if
end subroutine
