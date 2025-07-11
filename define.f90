module deftype
type atomtype
character name*2 !name of atom
integer index !The index in periodic table, if ECP was used, charge will smaller than this value
real*8 x,y,z,charge !Coordinate (Bohr) and (atomic or effective nuclear) charge of atoms. According to .fch file, -Bq atom has charge of 0
!resid and resname are filled when loading pdb/pqr/gro, and in this case "iresinfo" will be 1. &
!In other cases, iresinfo=0, resid will be set to 1 and resname will be space, subroutine readinfile do these
integer resid
character resname*4
end type

type primtype
integer center,type !The number of nuclei that the basis function centered on and its function type
real*8 exp !Exponent
end type

type content !Type for grid data points
real*8 x,y,z,value
end type

end module

!============ Store globally shared information
module defvar
use deftype
integer, private :: i
real*8,parameter :: pi=3.141592653589793D0
real*8,parameter :: b2a=0.529177210903D0 !Bohr to Angstrom, 2018 CODATA value
real*8,parameter :: au2kcal=627.51D0,au2KJ=2625.5D0,au2eV=27.2113838D0,au2cm=219474.6363D0,cal2J=4.184D0
real*8,parameter :: masse=9.10938215D-31,lightc=2.99792458D8,au2debye=2.5417462D0 !masse/chge: Mass/charge of an electron
real*8,parameter :: planckc=6.62606896D-34,h_bar=1.054571628D-34,amu2kg=1.66053878D-27
real*8,parameter :: boltzc=1.3806488D-23,boltzcau=3.1668114D-6,boltzceV=8.6173324D-5 !in J/K, in Hartree/K and in eV/K, respectively
real*8,parameter :: avogacst=6.02214179D23,eV2nm=1239.842D0,au2nm=45.563D0 !eV and nm can be interconverted by dividing eV2nm
integer,parameter :: nelesupp=150 !The number of elements supported, ghost(index=0) is not taken into account
integer,parameter :: maxneigh=30 !Maximum number of neighbours in neighbouring list
real*8 ctrval(1000) !Value of contour lines

!Store important calculated data
real*8,allocatable :: curvex(:),curvey(:),curveytmp(:) !For line plot
real*8,allocatable :: planemat(:,:),planemattmp(:,:) !planemattmp is mainly used to draw contour line of a function on contour map of another function (e.g. vdw surface on ESP contour map)
real*8,allocatable :: cubmat(:,:,:) !cubmat, store density/laplacian...3D-matrix
real*8,allocatable :: cubmatvec(:,:,:,:) !Used to store vector field
real*8,allocatable :: cubmattmp(:,:,:),cubmattmp2(:,:,:) !For temporarily exchanging and storing grid data, with same size as cubmat
real*8,allocatable :: rhocub(:,:,:) !Specifically store electron density grid data, with same size as cubmat
real*8,allocatable :: gradd1(:,:),gradd2(:,:) !Gradient in direction1/2 for gradient line plot
real*8,allocatable :: distmat(:,:) !Distance matrix, in Bohr
character filename*200,firstfilename*200 !firstfilename is the filename loaded when Multiwfn boots up or in option -11 of main menu
character extctrsetting*80 !cmdarg is the parameter of booting multiwfn
character,allocatable :: custommapname(:)*200,customop(:) !Custom operation for custom map/cube file
logical alive
integer :: ifragcontri=0 !=1 means fragment has been defined by users
integer :: nfragatm=0,nfragatm_org=0 !Number of elements of fragatm and fragatm_org
integer,allocatable :: fragatm(:),fragatm_org(:) !Store the index of atoms in fragment, used in various routines like setpromol. has no relationship with frag1/frag2. fragatm_org is used to backup fragatm during custom operation
integer,allocatable :: frag1(:),frag2(:) !These two fragments are only used for bond order analysis/composition analysis etc., store index of basis functions or atoms. Their size just fit their content
integer nfrag1,nfrag2 !Current size of frag1 and frag2, they should identical to size(frag1) and size(frag2)
integer :: ncustommap=0,imodwfn=0 !if 1, means occupation number or orbital type or basis function information has been modified
integer :: iorbsel=1 !Which orbital is selected, and its value will be calculated by fmo and calchessmat_mo
integer :: iorbsel2=0 !Which orbital will be plotted together with iorbsel in plane map
integer :: iorbvis=0 !The index of the orbital selected in main function 0
integer :: if_initlibreta=0 !If LIBRETA has been initialized for present wavefunction by running "call initlibreta"

integer,parameter :: ncolormax=16,ngoodcolor=15
character(len=10) :: colorname(ncolormax)=(/ character(len=10) :: "Red","Green","Blue","White","Black","Gray","Cyan","Yellow","Orange","Magenta","Crimson","Dark green","Purple","Brown","Dark blue","Pink" /) !Color name involved setcolor/selcolor routine
integer :: goodcolor(ngoodcolor)=(/5,1,3,12,9,10,14,13,11,15,16,2,7,8,6/) !Color list suitable for plotting lines, good colors appear prior to bad ones. Black,Red,Blue,Dark green,Orange,Magenta,Brown,Purple,Crimson,Dark blue,Pink,Green,Cyan,Yellow,Gray
!The name for superheavy atoms are consistent with Stuttgart PP website: http://www.tc.uni-koeln.de/PP/clickpse.en.html
character(len=2) :: ind2name(0:nelesupp)=(/ "Bq","H ","He", &   !Bq(number 0) is ghost atom. Bq is recorded in .fch, but X is not recorded
"Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
"Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og","Un","Ux",("??",i=121,nelesupp) /) !104~all. Name is in line with NIST periodic table. Such as Uun is replaced by Un
character(len=2) :: ind2name_up(0:nelesupp)=(/ "BQ","H ","HE", & !Same as ind2name, but all characters are upper case, to cater to .pdb file
"LI","BE","B ","C ","N ","O ","F ","NE", & !3~10
"NA","MG","AL","SI","P ","S ","CL","AR", & !11~18
"K ","CA","SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR", & !19~36
"RB","SR","Y ","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I ","XE", & !37~54
"CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU", & !55~71
"HF","TA","W ","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN", & !72~86
"FR","RA","AC","TH","PA","U ","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR", & !87~103
"RF","DB","SG","BH","HS","MT","DS","RG","CN","UT","FL","UP","LV","US","UO","UN","UX",("??",i=121,nelesupp) /) !104~all
!Bondi vdW radii, from J.Phys.Chem.,1964,68(3),441-451, unit is Angstrom, will be convert to Bohr when multiwfn start
real*8 :: vdwr(0:nelesupp)=(/ 0.3D0,1.2D0,1.4D0,& !Ghost,H,He
1.82D0,1.77D0,1.74D0,1.7D0,1.55D0,1.52D0,1.47D0,1.54D0,& !Li~Ne
2.27D0,1.73D0,1.73D0,2.1D0,1.8D0,1.8D0,1.75D0,1.88D0,& !Na~Ar
(2D0,i=19,27),1.63D0,1.4D0,1.39D0,1.87D0,2D0,1.85D0,1.9D0,1.85D0,2.02D0,& !Ni~Kr(28~36)
(2D0,i=37,45),1.63D0,1.72D0,1.58D0,1.93D0,2.17D0,2D0,2.06D0,1.98D0,2.16D0,& !Pd~Xe(46~54)
(2D0,i=55,77),1.72D0,1.66D0,1.55D0,1.96D0,2.02D0,(2D0,i=83,nelesupp) /) !Pt~Pb(78~82)
!##No use currently!## Modified Bondi vdW radii, but for all main group (except for H,He), use IVA radius in corresponding row. Specifically used to molecular surface decomposition
real*8 :: vdwr_tianlu(0:nelesupp)=(/ 0.4D0,1.7D0,1.7D0,& !Ghost,H,He   H and Ne are augmented to carbon radius
(1.7D0,i=3,10),& !Li~Ne
(2.1D0,i=11,18),& !Na~Ar
1.87D0,1.87D0,  (2D0,i=21,27),1.63D0,1.40D0,1.39D0,  (1.87D0,i=31,36),& !K ,Ca,  Ni~Zn(21~30),  Ga~Kr(31,37)
1.93D0,1.93D0,  (2D0,i=39,45),1.63D0,1.72D0,1.58D0,  (1.93D0,i=49,54),& !Rb,Sr,  Y ~Cd(39~48),  In~Xe(49~54)
1.96D0,1.96D0,  (2D0,i=57,77),1.72D0,1.66D0,1.55D0,  (1.96D0,i=81,86),& !Cs,Ba,  La~Hg(57~80),  Tl~Rn(81~86)
(2D0,i=87,nelesupp) /) !Rn~Mt(87~109,~all)
!UFF vdW radii, from x of Table 1 of JACS, 114, 10024 (1992). The values are half of x
real*8 :: vdwr_UFF(0:nelesupp)=(/ 0.4D0,1.443D0,1.181D0,& !Ghost,H,He
1.2255D0,1.3725D0,2.0415D0,1.9255D0,1.83D0,1.75D0,1.682D0,1.6215D0,& !Li~Ne
1.4915D0,1.5105D0,2.2495D0,2.1475D0,2.0735D0,2.0175D0,1.9735D0,1.934D0,& !Na~Ar
1.906D0,1.6995D0,1.6475D0,1.5875D0,1.572D0,1.5115D0,1.4805D0,1.456D0,1.436D0,&
1.417D0,1.7475D0,1.3815D0,2.1915D0,2.14D0,2.115D0,2.1025D0,2.0945D0,2.0705D0,& !K~Kr(19~36)
2.057D0,1.8205D0,1.6725D0,1.562D0,1.5825D0,1.526D0,1.499D0,1.4815D0,1.4645D0,&
1.4495D0,1.574D0,1.424D0,2.2315D0,2.196D0,2.21D0,2.235D0,2.25D0,2.202D0,& !Rb~Xe(37~54)
2.2585D0,1.8515D0,1.761D0,1.778D0,1.803D0,1.7875D0,1.7735D0,1.76D0,1.7465D0,1.684D0,1.7255D0,1.714D0,1.7045D0,1.6955D0,1.687D0,1.6775D0,1.82D0,& !Cs~Lu(55~71)
1.5705D0,1.585D0,1.5345D0,1.477D0,1.56D0,1.42D0,1.377D0,1.6465D0,1.3525D0,2.1735D0,2.1485D0,2.185D0,2.3545D0,2.375D0,2.3825D0,& !Hf~Rn(72~86)
2.45D0,1.8385D0,1.739D0,1.698D0,1.712D0,1.6975D0,1.712D0,1.712D0,1.6905D0,1.663D0,1.6695D0,1.6565D0,1.6495D0,1.643D0,1.637D0,1.624D0,1.618D0,& !Fr~Lr(87~103)
(2D0,i=104,nelesupp) /) !Rf~(104~all)

!Covalent radius, from "Dalton Trans., 2008, 2832-2838", unit is Angstrom, will be convert to Bohr when Multiwfn boots up
real*8 :: covr(0:nelesupp)=(/ 0.1D0,0.31D0,0.28D0,& !Ghost,H,He(1~2)
1.28D0,0.96D0,0.84D0,0.76D0,0.71D0,0.66D0,0.57D0,0.58D0,& !Li~Ne(3~10)     here C is sp3
1.66D0,1.41D0,1.21D0,1.11D0,1.07D0,1.05D0,1.02D0,1.06D0,& !Na~Ar(11~18)
2.03D0,1.76D0,1.70D0,1.60D0,1.53D0,1.39D0,1.39D0,1.32D0,1.26D0,& !K~Co(19~27)  here MnD0,FeD0,Co is low-spinD0, high spin is 1.61D0,1.52D0,1.50
1.24D0,1.32D0,1.22D0,1.22D0,1.20D0,1.19D0,1.20D0,1.20D0,1.16D0,& !Ni~Kr(28~36)
2.20D0,1.95D0,1.90D0,1.75D0,1.64D0,1.54D0,1.47D0,1.46D0,1.42D0,& !Rb~Rh(37~45)
1.39D0,1.45D0,1.44D0,1.42D0,1.39D0,1.39D0,1.38D0,1.39D0,1.40D0,& !Pd~Xe(46~54)
2.44D0,2.15D0,2.07D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,& !Cs~Tb(55~65)
1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,1.75D0,1.70D0,1.62D0,1.51D0,1.44D0,1.41D0,& !Dy~Ir(66~77)
1.36D0,1.36D0,1.32D0,1.45D0,1.46D0,1.48D0,1.40D0,1.50D0,1.50D0,2.60D0,2.21D0,& !Pt~Ra(78~88)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,i=97,nelesupp) /) !Ac~Cm(89~96),~all
!(Covalent) radius proposed by Suresh, from J. Phys. Chem. A 2001, 105, 5940-5944. For missing values (including all noble gases and very heavy elements), the ones in covr array are used
!Unit is Angstrom, will be convert to Bohr when Multiwfn boots up
real*8 :: covr_Suresh(0:nelesupp)=(/ 0.1D0,0.327D0,0.28D0,& !Ghost,H,He(1~2)
1.219D0,0.911D0,0.793D0,0.766D0,0.699D0,0.658D0,0.633D0,0.58D0,& !Li~Ne(3~10)
1.545D0,1.333D0,1.199D0,1.123D0,1.11D0,1.071D0,1.039D0,1.06D0,& !Na~Ar(11~18)
1.978D0,1.745D0,1.337D0,1.274D0,1.236D0,1.128D0,1.18D0,1.091D0,1.089D0,& !K~Co(19~27)
1.077D0,1.146D0,1.187D0,1.199D0,1.179D0,1.209D0,1.201D0,1.201D0,1.16D0,& !Ni~Kr(28~36)
2.217D0,1.928D0,1.482D0,1.377D0,1.353D0,1.24D0,1.287D0,1.212D0,1.229D0,& !Rb~Rh(37~45)
1.24D0,1.362D0,1.429D0,1.385D0,1.38D0,1.421D0,1.4D0,1.397D0,1.40D0,& !Pd~Xe(46~54)
2.442D0,2.149D0,1.653D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,& !Cs~Tb(55~65)
1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,1.364D0,1.346D0,1.256D0,1.258D0,1.222D0,1.227D0,& !Dy~Ir(66~77)
1.227D0,1.273D0,1.465D0,1.531D0,1.434D0,1.496D0,1.40D0,1.50D0,1.50D0,2.60D0,2.21D0,& !Pt~Ra(78~88)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,i=97,nelesupp) /) !Ac~Cm(89~96),~all
!Covalent radius, from Pyykko "Chem. Eur. J.,15,186 (2009)", unit is in Angstrom, will be convert to Bohr when Multiwfn boots up
real*8 :: covr_pyy(0:nelesupp)=(/ 0.1D0,0.32D0,0.46D0,& !Ghost,H,He(1~2)
1.33D0,1.02D0,0.85D0,0.75D0,0.71D0,0.63D0,0.64D0,0.67D0,& !Li~Ne(3~10)
1.55D0,1.39D0,1.26D0,1.16D0,1.11D0,1.03D0,0.99D0,0.96D0,& !Na~Ar(11~18)
1.96D0,1.71D0,1.48D0,1.36D0,1.34D0,1.22D0,1.19D0,1.16D0,1.11D0,& !K~Co(19~27)
1.10D0,1.12D0,1.18D0,1.24D0,1.21D0,1.21D0,1.16D0,1.14D0,1.17D0,& !Ni~Kr(28~36)
2.10D0,1.85D0,1.63D0,1.54D0,1.47D0,1.38D0,1.28D0,1.25D0,1.25D0,& !Rb~Rh(37~45)
1.20D0,1.28D0,1.36D0,1.42D0,1.40D0,1.40D0,1.36D0,1.33D0,1.31D0,& !Pd~Xe(46~54)
2.32D0,1.96D0,1.80D0,1.63D0,1.76D0,1.74D0,1.73D0,1.72D0,1.68D0,1.69D0,1.68D0,& !Cs~Tb(55~65)
1.67D0,1.66D0,1.65D0,1.64D0,1.70D0,1.62D0,1.52D0,1.46D0,1.37D0,1.31D0,1.29D0,1.22D0,& !Dy~Ir(66~77)
1.23D0,1.24D0,1.34D0,1.44D0,1.44D0,1.51D0,1.45D0,1.47D0,1.42D0,2.23D0,2.01D0,& !Pt~Ra(78~88)
1.86D0,1.75D0,1.69D0,1.70D0,1.71D0,1.72D0,1.66D0,1.66D0,1.68D0,1.68D0,1.65D0,1.67D0,1.73D0,1.76D0,1.61D0,& !Ac~Lr(89~103)
1.57D0,1.49D0,1.43D0,1.41D0,1.34D0,1.29D0,1.28D0,1.21D0,1.22D0,1.36D0,1.43D0,1.62D0,1.75D0,1.65D0,1.57D0,(1.5D0,i=119,nelesupp)  /) !Rf~118(104~118),~all
real*8 :: covr_tianlu(0:nelesupp)=(/ 0.1D0,0.31D0,0.28D0,& !H,Ne(1~2) !Based on CSD radii, but for all main group (except for H,He), use IVA radius in corresponding row
(0.76D0,i=3,10),& !Li~Ne(3~10)
(1.11D0,i=11,18),1.2D0,1.2D0,& !Na~Ar(11~18),K,Ca
1.70D0,1.60D0,1.53D0,1.39D0,1.39D0,1.32D0,1.26D0,1.24D0,1.32D0,1.22D0,& !Sc~Zn(21~30)  here MnD0,FeD0,Co is low-spinD0, high spin is 1.61D0,1.52D0,1.50
(1.2D0,i=31,36),1.42D0,1.42D0,& !Ga~Kr(31~36),Rb,Sr
1.90D0,1.75D0,1.64D0,1.54D0,1.47D0,1.46D0,1.42D0,1.39D0,1.45D0,1.44D0,& !Y~Cd(39~48)
(1.39D0,i=49,54),1.46D0,1.46D0,& !In~Xe(49~54),Cs,Ba
2.07D0,2.04D0,2.03D0,2.01D0,1.99D0,1.98D0,1.98D0,1.96D0,1.94D0,1.92D0,1.92D0,1.89D0,1.90D0,1.87D0,1.87D0,& !La~Lu(57~71)
1.75D0,1.70D0,1.62D0,1.51D0,1.44D0,1.41D0,1.36D0,1.36D0,1.32D0,& !Hf~Hg(72~80)
(1.46D0,i=81,86),1.46D0,1.46D0,&!Tl~Rn(81~86),Fr(still 1.46),Ra(still 1.46)
2.15D0,2.06D0,2.00D0,1.96D0,1.90D0,1.87D0,1.80D0,1.69D0,(1.5D0,i=97,nelesupp) /) !Ac~Cm(89~96),~all
!Radii proposed in Chem. Phys. Lett., 480 (2009) 127-131, the unit is Bohr!
real*8 :: radii_Hugo(0:nelesupp)=(/ 0.10D0,1.00D0,0.74D0,& !Ghost,H,Ne(1~2)
1.59D0,1.21D0,1.28D0,1.10D0,0.97D0,1.00D0,0.88D0,0.79D0,& !Li~Ne(3~10)
1.63D0,1.33D0,1.51D0,1.29D0,1.14D0,1.15D0,1.02D0,0.93D0,& !Na~Ar(11~18)
1.77D0,1.49D0,1.44D0,1.41D0,1.42D0,1.42D0,1.35D0,1.31D0,1.31D0,1.33D0,1.33D0,1.20D0,1.51D0,1.31D0,1.18D0,1.18D0,1.07D0,0.99D0,& !K~Kr
1.80D0,1.55D0,1.48D0,1.43D0,1.42D0,1.38D0,1.37D0,1.36D0,1.35D0,1.28D0,1.34D0,1.23D0,1.53D0,1.36D0,1.26D0,1.23D0,1.14D0,1.06D0,1.87D0,1.62D0,& !Rb~Xe,Cs,Ba
1.56D0,1.57D0,1.58D0,1.57D0,1.56D0,1.55D0,1.55D0,1.49D0,1.52D0,1.51D0,1.50D0,1.49D0,1.48D0,1.47D0,1.58D0,& !La~Lu
1.41D0,1.34D0,1.31D0,1.32D0,1.27D0,1.23D0,1.23D0,1.21D0,1.14D0,1.49D0,1.35D0,1.37D0,1.27D0,1.21D0,1.12D0,1.83D0,1.16D0,& !Hf~Rn,Fr,Ra
1.62D0,1.47D0,1.52D0,1.48D0,1.47D0,1.50D0,1.51D0,1.51D0,1.48D0,1.47D0,1.46D0,1.45D0,1.44D0,1.43D0,1.67D0,1.51D0,(1.5D0,i=105,nelesupp) /) !Ac~Rf,~all

real*8 :: YWTatomcoeff(18,3)=reshape((/ & !Coef. of fitting B3LYP/6-31G* density by Weitao Yang group for the first three rows, see supporting info. of JACS,132,6498. I found the fitting quality is poor!
0.2815D0,2.437D0,11.84D0,31.34D0,67.82D0,120.2D0,190.9D0,289.5D0,406.3D0,561.3D0,760.8D0,1016.0D0,1319.0D0,1658.0D0,2042.0D0,2501.0D0,3024.0D0,3625.0D0, &
0D0,0D0,0.06332D0,0.3694D0,0.8527D0,1.172D0,2.247D0,2.879D0,3.049D0,6.984D0,22.42D0,37.17D0,57.95D0,87.16D0,115.7D0,158.0D0,205.5D0,260.0D0, &
0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0.06358D0,0.3331D0,0.8878D0,0.7888D0,1.465D0,2.17D0,3.369D0,5.211D0 /),(/18,3/))
real*8 :: YWTatomexp(18,3)=reshape((/ & !Corresponding exponent of YWTatom, the value set to 1.0 doesn't have any meaning, only for avoiding divide zero. Note that this is not exponent of STO, but its inverse!
0.5288D0,0.3379D0,0.1912D0,0.139D0,0.1059D0,0.0884D0,0.0767D0,0.0669D0,0.0608D0,0.0549D0,0.0496D0,0.0449D0,0.0411D0,0.0382D0,0.0358D0,0.0335D0,0.0315D0,0.0296D0, &
1.0D0,1.0D0,0.9992D0,0.6945D0,0.53D0,0.548D0,0.4532D0,0.3974D0,0.3994D0,0.3447D0,0.2511D0,0.215D0,0.1874D0,0.1654D0,0.1509D0,0.1369D0,0.1259D0,0.1168D0, &
1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0236D0,0.7753D0,0.5962D0,0.6995D0,0.5851D0,0.5149D0,0.4974D0,0.4412D0 /),(/18,3/))
!The largest distance (Bohr) of non-zero radial density of every element in atmraddens.f90, used to truncate calculation for saving time
real*8 :: atmrhocut(0:nelesupp)=(/ 14D0,7.80D0,5.89D0,14.37D0,11.30D0,10.10D0,8.20D0,7.08D0,6.75D0,6.16D0,5.64D0,15.33D0,12.70D0,12.70D0,10.68D0,9.57D0,9.08D0,8.20D0,7.80D0,16.38D0,&
14.37D0,14.37D0,14.37D0,13.50D0,13.50D0,13.50D0,12.70D0,12.70D0,11.97D0,11.97D0,11.30D0,13.50D0,11.97D0,10.68D0,10.68D0,10.10D0,9.08D0,16.38D0,&
15.33D0,14.37D0,14.37D0,13.50D0,11.97D0,14.37D0,12.70D0,12.70D0,9.57D0,12.70D0,11.97D0,14.37D0,12.70D0,11.97D0,11.30D0,10.68D0,10.10D0,17.53D0,&
16.38D0,15.33D0,15.33D0,15.33D0,15.33D0,15.33D0,15.33D0,15.33D0,14.37D0,14.37D0,14.37D0,14.37D0,14.37D0,14.37D0,14.37D0,14.37D0,11.97D0,11.97D0,&
11.97D0,11.97D0,12.70D0,11.97D0,11.30D0,11.30D0,11.30D0,14.37D0,12.70D0,11.30D0,11.30D0,10.68D0,10.10D0,17.53D0,15.33D0,14.37D0,14.37D0,14.37D0,&
15.33D0,15.33D0,14.37D0,14.37D0,13.50D0,14.37D0,14.37D0,14.37D0,14.37D0,14.37D0,14.37D0,12.70D0,(14D0,i=104,nelesupp) /)
real*8 atmrhocutsqr(0:nelesupp) !Square of atmrhocut, generate when Multiwfn boots up
!The distance (Bohr) of cutting radial density at 1E-5 a.u. for every element in atmraddens.f90, used to truncate calculation for saving time
real*8 :: atmrhocut_1En5(0:nelesupp)=(/ 14D0,&
5.317D0, 4.167D0, 8.622D0, 7.267D0, 6.707D0, 5.977D0, 5.351D0, 5.009D0, 4.645D0, 4.304D0, 8.863D0, 8.071D0, 8.235D0, 7.506D0, 6.847D0, &
6.497D0, 6.115D0, 5.753D0, 9.881D0, 9.214D0, 9.163D0, 8.920D0, 8.703D0, 8.216D0, 8.383D0, 8.142D0, 7.959D0, 7.780D0, 7.737D0, 7.527D0, &
8.243D0, 7.731D0, 7.221D0, 6.971D0, 6.636D0, 6.307D0,10.191D0, 9.608D0, 9.477D0, 9.237D0, 8.446D0, 8.140D0, 8.728D0, 8.059D0, 8.020D0, &
6.519D0, 7.959D0, 7.821D0, 8.536D0, 8.071D0, 7.677D0, 7.458D0, 7.147D0, 6.815D0,10.659D0,10.204D0, 9.644D0, 9.559D0, 9.827D0, 9.750D0, &
9.676D0, 9.609D0, 9.545D0, 9.204D0, 9.386D0, 9.314D0, 9.250D0, 9.188D0, 9.129D0, 9.070D0, 8.850D0, 8.639D0, 8.508D0, 8.355D0, 8.282D0, &
8.047D0, 7.790D0, 7.375D0, 7.272D0, 7.330D0, 8.674D0, 8.209D0, 7.741D0, 7.524D0, 7.227D0, 6.939D0,10.524D0,10.104D0, 9.687D0, 9.247D0, &
9.399D0, 9.701D0, 9.618D0, 9.387D0, 9.327D0, 9.032D0, 9.163D0, 9.099D0, 9.033D0, 8.966D0, 8.906D0, 8.850D0, 8.619D0,(14D0,i=104,nelesupp) /)
real*8 atmrhocutsqr_1En5(0:nelesupp) !Square of atmrhocut_1En5, generate when Multiwfn boots up

!Standard atomic weight. Computed from abundance and isotope masses, the data was obtained from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all
!For radioactive elements, the mass corresponds to longest-living isotope
!Unit is Da, or simply be regarded as relative atomic mass (dimensionless)
real*8 :: atmwei(0:nelesupp)=(/0D0,&
  1.0079407541D0,  4.0026019321D0,  6.9400366029D0,  9.0121830650D0, 10.8110280464D0,&
 12.0107358967D0, 14.0067032114D0, 15.9994049243D0, 18.9984031627D0, 20.1800463805D0,&
 22.9897692820D0, 24.3050516198D0, 26.9815385300D0, 28.0854987057D0, 30.9737619984D0,&
 32.0647874061D0, 35.4529375826D0, 39.9477985636D0, 39.0983009101D0, 40.0780225110D0,&
 44.9559082800D0, 47.8667449627D0, 50.9414650374D0, 51.9961317554D0, 54.9380439100D0,&
 55.8451444339D0, 58.9331942900D0, 58.6933471099D0, 63.5460399458D0, 65.3777825295D0,&
 69.7230660726D0, 72.6275501647D0, 74.9215945700D0, 78.9593885570D0, 79.9035277805D0,&
 83.7979999953D0, 85.4676635956D0, 87.6166444696D0, 88.9058403000D0, 91.2236415971D0,&
 92.9063730000D0, 95.9597885412D0, 97.9072124000D0,101.0649401392D0,102.9054980000D0,&
106.4153275073D0,107.8681496346D0,112.4115578183D0,114.8180866294D0,118.7101125930D0,&
121.7597836735D0,127.6031264847D0,126.9044719000D0,131.2927614478D0,132.9054519610D0,&
137.3268916286D0,138.9054688737D0,140.1157307379D0,140.9076576000D0,144.2415960318D0,&
144.9127559000D0,150.3663557119D0,151.9643781264D0,157.2521306469D0,158.9253547000D0,&
162.4994728194D0,164.9303288000D0,167.2590826497D0,168.9342179000D0,173.0541501663D0,&
174.9668149579D0,178.4849787234D0,180.9478756362D0,183.8417775505D0,186.2067045456D0,&
190.2248596282D0,192.2160516521D0,195.0844568649D0,196.9665687900D0,200.5991670346D0,&
204.3834128394D0,207.2169080630D0,208.9803991000D0,208.9824308000D0,209.9871479000D0,&
222.0175782000D0,223.0197360000D0,226.0254103000D0,227.0277523000D0,232.0380558000D0,&
231.0358842000D0,238.0289104617D0,237.0481736000D0,244.0642053000D0,243.0613813000D0,&
247.0703541000D0,247.0703073000D0,251.0795886000D0,252.0829800000D0,257.0951061000D0,&
258.0984315000D0,259.1010300000D0,266.1198300000D0,267.1217900000D0,268.1256700000D0,&
271.1339300000D0,270.1333600000D0,269.1337500000D0,278.1563100000D0,281.1645100000D0,&
282.1691200000D0,285.1771200000D0,286.1822100000D0,289.1904200000D0,289.1936300000D0,&
293.2044900000D0,294.2104600000D0,294.2139200000D0,(0D0,i=119,nelesupp)/)
 
!Series of Lebedev-Laikov routines
integer :: Lebelist(32)=(/ 6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810 /)
integer :: fact(0:10)=(/ 1,1,2,6,24,120,720,5040,40320,362880,3628800 /) ! Store factorials from 0 to 10 
integer :: isphergau=0 !By default, all basis functions are Cartesian type, =1 means spherical (but some of them can still be Cartesian type)
character(len=5) :: GTFtype2name(-32:56)=(/ & !Definition of such as G-4, H+5 can be found in http://sobereva.com/97
"H 0  ","H+1  ","H-1  ","H+2  ","H-2  ","H+3  ","H-3  ","H+4  ","H-4  ","H+5  ","H-5  ", & !-32:-22
"G 0  ","G+1  ","G-1  ","G+2  ","G-2  ","G+3  ","G-3  ","G+4  ","G-4  ", & !-21:-13
"F 0  ","F+1  ","F-1  ","F+2  ","F-2  ","F+3  ","F-3  ","D 0  ","D+1  ","D-1  ","D+2  ","D-2  ", & !-12:-6,-5:-1
"     ","S    ","X    ","Y    ","Z    ","XX   ","YY   ","ZZ   ","XY   ","XZ   ","YZ   ", & !0~10
"XXX  ","YYY  ","ZZZ  ","XXY  ","XXZ  ","YYZ  ","XYY  ","XZZ  ","YZZ  ","XYZ  ", & !f 11~20
"ZZZZ ","YZZZ ","YYZZ ","YYYZ ","YYYY ","XZZZ ","XYZZ ","XYYZ ","XYYY ","XXZZ ","XXYZ ","XXYY ","XXXZ ","XXXY ","XXXX ", & !g 21~35
"ZZZZZ","YZZZZ","YYZZZ","YYYZZ","YYYYZ","YYYYY","XZZZZ","XYZZZ","XYYZZ","XYYYZ","XYYYY","XXZZZ","XXYZZ","XXYYZ","XXYYY","XXXZZ","XXXYZ","XXXYY","XXXXZ","XXXXY","XXXXX" /) !h 36~56
!I shell (not supported yet): ZZZZZZ, YZZZZZ, YYZZZZ, YYYZZZ, YYYYZZ, YYYYYZ, YYYYYY, XZZZZZ, XYZZZZ, XYYZZZ, XYYYZZ, XYYYYZ, XYYYYY, &
!XXZZZZ, XXYZZZ, XXYYZZ, XXYYYZ, XXYYYY, XXXZZZ, XXXYZZ, XXXYYZ, XXXYYY, XXXXZZ, XXXXYZ, XXXXYY, XXXXXZ, XXXXXY, XXXXXX
character(len=5) :: type2ang(56)=(/ &
"S    ","P    ","P    ","P    ","D    ","D    ","D    ","D    ","D    ","D    ", & !0~10
"F    ","F    ","F    ","F    ","F    ","F    ","F    ","F    ","F    ","F    ", & !f 11~20
"G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ","G    ", & !g 21~35
"H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    ","H    " /) !h 36~56
!The order of spherical basis functions are the same for all kinds of files, however, the order of Cartesian ones may be different in different cases
!Here s,p,d sequences are identical to .wfn, .wfx, .fch, .molden  !Note: Sequence in .fch = sequence in Gaussian (the printed basis functions via pop=full)
!Here f sequence is identical to .wfn, .wfx, but not identical to .fch and .molden
!Here g sequence is identical to .fch, .wfn does not support higher than f function, not identical to .wfx and .molden
!here h sequence is identical to .wfx and .fch, .molden doesn't support h
!Notice: The .wfn produced by G09 B.01 and later supports g and h, the definition is identical to here, and thus can be normally loaded
!Overall, spd: Multiwfn=wfn=wfx=fch=molden   f: Multiwfn=wfn=wfx!=(fch=molden)   g: Multiwfn=fch!=(wfx=molden=Molden2AIM)   h: Multiwfn=wfx=fch
integer :: type2ix(56)=(/ 0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,0,1,1,0,1, 0,0,0,0,0,1,1,1,1,2,2,2,3,3,4, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5 /)
integer :: type2iy(56)=(/ 0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1, 0,1,2,3,4,0,1,2,3,0,1,2,0,1,0, 0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0 /)
integer :: type2iz(56)=(/ 0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,1,0,2,2,1, 4,3,2,1,0,3,2,1,0,2,1,0,1,0,0, 5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0 /)
!Negative value means the shell use spherical Gaussian function. -1=SP (also known as L in GAMESS), and impossible be used in Multiwfn (when detect it, split it as S and P)
character :: shtype2name(-5:5)=(/ "H","G","F","D","L","S","P","D","F","G","H" /)
!Convert shell type to the number of basis functions in the shell: 0=s,1=p,-1=sp,2=6d,-2=5d,3=10f,-3=7f,4=15g,-4=9g,5=21h,-5=11h
integer :: shtype2nbas(-5:5)=(/ 11,9,7,5,4,1,3,6,10,15,21 /) 

!-------- Variables for wavefunction information (_org means using for backuping the first loaded molecule)
integer :: ibasmode=0 !0/1 = GTO/STO is used in current wavefunction
integer :: nmo=0,nprims=0,ncenter=0,ncenter_org=0,nmo_org=0,nprims_org=0 !Number of orbitals, primitive functions, nuclei
integer :: idxHOMO=0,idxHOMOb=0 !shtype of total/alpha and beta-HOMO. Can be determined by subroutine getHOMOidx
integer :: ifiletype=0 !Plain text=0, fch/fchk=1, wfn=2, wfx=3, chg/pqr=4, pdb/xyz=5, NBO .31=6, cube/VASP grid=7, grd/dx/vti=8, molden=9, gms=10, MDL mol/sdf=11, gjf or ORCA inp or mop =12, mol2=13, mwfn=14, gro=15, cp2k input=16, cif=17, POSCAR=18, QE input=19
integer :: wfntype=0 !0/1/2= R/U/RO single determinant wavefunction, 3/4=R/U multiconfiguration wavefunction
real*8 :: totenergy=0,virialratio=2,nelec=0,naelec=0,nbelec=0
integer :: loadmulti=-99,loadcharge=-99 !Spin multiplicity and net charge, loaded directly from input file (e.g. from .gjf, Gaussian .out, ORCA .inp, title line of .xyz), only utilized in rare cases. -99 means unloaded
real*8 :: kp1crd=0,kp2crd=0,kp3crd=0 !k-point fractional coordinate in three directions in reciprocal space for present wavefunction
!-------- Variables for nuclei & GTF & Orbitals. Note: Row and column of CO(:,:) correspond to orbital and GTF, respectively, in contrary to convention
type(atomtype),allocatable :: a(:),a_org(:),a_tmp(:) !a_tmp is only used in local temporary operation, should be destoried immediatedly after using
integer,allocatable :: a_tmp_idx(:) !Record actual index of every atom in a_tmp
type(primtype),allocatable,target :: b(:)
type(primtype),allocatable :: b_org(:),b_tmp(:)
real*8,allocatable :: MOocc(:),MOocc_org(:),MOene(:),MOene_org(:) !Occupation number & energy of orbital
integer,allocatable :: MOtype(:) !The type of orbitals, (alpha&beta)=0/alpha=1/beta=2, not read from .wfn directly
character(len=10) :: orbtypename(0:2)=(/ character(len=10) :: "Alpha&Beta","Alpha","Beta" /)
character(len=3) :: orbtypename_short(0:2)=(/ character(len=10) :: "A+B"," A "," B " /)
character(len=4),allocatable :: MOsym(:) !The symmetry of orbitals, meaningful when .mwfn/molden/gms is used
real*8,allocatable, target :: CO(:,:) !Coefficient matrix of primitive basis functions, including both normalization and contraction coefficients
real*8,allocatable :: CO_org(:,:),CO_tmp(:,:)
real*8,allocatable :: COtr(:,:) !Transposed CO matrix, which is used in some routines for faster calculation than using CO. Must be deallocated after using
!Unique GTFs (the GTFs with identical center, type and exponent are combined together and leave only one). Can be activated after running gen_GTFuniq
integer :: nprims_uniq=0 !0 means uninitialized
type(primtype),allocatable,target :: b_uniq(:) !b of unique GTFs
real*8,allocatable,target :: CO_uniq(:,:) !CO of b_uniq. The coefficients of duplicated GTFs are summed together
!-------- Describe core electron density in EDF section. If nEDFprims=0, then EDF will not be considered in calculation
type(primtype),allocatable :: b_EDF(:)
real*8,allocatable :: CO_EDF(:)
integer :: nEDFprims=0,nEDFelec=0 !Number of EDFs, number of electrons represented by EDFs
integer,allocatable :: nEDFelecatm(:) !Number of electrons represented by EDF of each atom
!-------- Promolecular wavefunction. Share same "a" and "b"
real*8,allocatable,target :: CO_pmol(:,:)
real*8,allocatable :: MOocc_pmol(:),MOene_pmol(:)
integer,allocatable :: MOtype_pmol(:)
integer :: nmo_pmol=0
!-------- Variables when basis functions are basis rather than primitive function as basis
integer :: nbasis=0,nbasisCar=0,nshell=0,nprimshell=0 !The number of basis (actual/Cartesian), basis shell and primitive shell. SP shell is counted as S and P shell separately
integer :: nindbasis=0 !Number of independent basis functions. Can directly load from fch and mwfn. For other case like molden, nindbasis keeps 0, means undetermined
integer,allocatable :: shtype(:),shtypeCar(:),shcon(:),shcen(:) !Type, contraction degree and attributed center of a basis shell
real*8,allocatable :: primshexp(:),primshcoeff(:) !Exponent and contraction coefficient of a primitive shell
   !Note: Sizes of basshell, bascen, etc. may be larger than the number of actual basis functions, because they are initially allocated for Cartesian basis functions, which may have larger number
integer,allocatable :: basshell(:) !The ith element is the shell index that the ith basis attributed to
integer,allocatable :: bascen(:),bastype(:) !Center/type of basis, definition is the same as GTF
integer,allocatable :: basstart(:),basend(:) !The ith element means the basis from where to where is attributed to the ith atom. If any of them is 0, that means this atom does not have basis function
integer,allocatable :: primstart(:),primend(:) !The ith element means the GTF from where to where is attributed to the ith Cartesian basis function (which may be yielded during reading spherical wavefunction)
real*8,allocatable :: primconnorm(:) !element i means the contract. coeff. * normalization coeff. of GTF i, can be used for e.g. constructing basis integral from GTF integral
real*8,allocatable :: Sbas(:,:),Sbas_org(:,:),Sbas_PBC(:,:) !Overlap matrix, its backup and special form for PBC calculation
real*8,allocatable,target :: Dbas(:,:,:) !Electric dipole moment integral matrix, the first index 1,2,3=X,Y,Z, the last two indices are basis index
real*8,allocatable,target :: Quadbas(:,:,:) !Quadrupole moment integral matrix, the first index 1~6=XX,YY,ZZ,XY,YZ,XZ
real*8,allocatable,target :: Octobas(:,:,:) !Octopole moment integral matrix, the first index 1~10=XXX,YYY,ZZZ,YZZ,XZZ,XXZ,YYZ,XXY,XYY,XYZ
real*8,allocatable,target :: Hexdebas(:,:,:) !Hexadecapole moment integral matrix, the first index 1~15=XXXX,YYYY,ZZZZ,XXXY,XXXZ,YYYX,YYYZ,ZZZX,ZZZY,XXYY,XXZZ,YYZZ,XXYZ,YYXZ,ZZXY (same as Gaussian)
real*8,allocatable,target :: DorbA(:,:,:),DorbB(:,:,:) !Electric dipole moment integral matrix between orbitals, the first index 1,2,3=X,Y,Z. A/B corresponds to (closed or Alpha)/Beta orbitals
real*8,allocatable,target :: Tbas(:,:) !Kinetic energy integral matrix
real*8,allocatable,target :: Vbas(:,:) !Nuclear attraction potential integral matrix
real*8,allocatable,target :: Velbas(:,:,:) !Velocity integral matrix, the first index 1,2,3=X,Y,Z
real*8,allocatable,target :: Magbas(:,:,:) !Magnetic integral matrix, the first index 1,2,3=X,Y,Z
real*8,allocatable,target :: MagorbA(:,:,:),MagorbB(:,:,:) !Magnetic dipole moment integral matrix between orbitals, the first index 1,2,3=X,Y,Z. A/B corresponds to (closed or Alpha)/Beta orbitals
!Coefficient matrix for alpha/beta orbital, CObasa(i,j) means the coefficient of ith basis in the jth orbital, differ to CO(:,:)
real*8,allocatable,target :: CObasa(:,:),CObasb(:,:) !wfntype==0,2,3 only allocate CObasa(nbasis,nmo), ==1,4 also allocate CObasb, dimension of both CObasa and CObasb would be (nbasis,nbasis)
real*8,allocatable,target :: CObasa_org(:,:),CObasb_org(:,:) !Used to temporarily backup in certain subroutines, e.g. Wiberg, ETS-NOCV. Do not immediately backup after file loading like a_org, b_org, CO_org, etc.
real*8,allocatable,target :: Ptot(:,:),Palpha(:,:),Pbeta(:,:) !Density matrix of total/alpha/beta, for wfntype==0.or.wfntype==3, only Ptot is filled, for others, all of Ptot,Palpha and Pbeta are filled
!Some matrics based on GTF
real*8,allocatable :: Ptot_prim(:,:),Palpha_prim(:,:),Pbeta_prim(:,:) !Density matrix of total/alpha/beta based on GTF
real*8,allocatable :: Dprim(:,:,:) !Dipole moment integral matrix based on GTF, the first index 1,2,3=X,Y,Z, the last two indices are basis index
real*8,allocatable :: Quadprim(:,:,:) !Quadrupole moment integral matrix based on GTF, the first index 1~6=XX,YY,ZZ,XY,YZ,XZ
real*8,allocatable :: Octoprim(:,:,:) !Octopole moment integral matrix based on GTF, the first index 1~10=XXX,YYY,ZZZ,YZZ,XZZ,XXZ,YYZ,XXY,XYY,XYZ
real*8,allocatable :: Hexdeprim(:,:,:) !Hexadecapole moment integral matrix based on GTF, the first index 1~15=XXXX,YYYY,ZZZZ,XXXY,XXXZ,YYYX,YYYZ,ZZZX,ZZZY,XXYY,XXZZ,YYZZ,XXYZ,YYXZ,ZZXY (same as Gaussian)
!Back up some information prior to calling delvirorb, so that can be restored via delvirorb_back
integer :: nmo_back,ifdelvirorb=0 !If "delvirorb" has been called, then ifdelvirorb=1, otherwise 0
real*8,allocatable :: MOene_back(:),MOocc_back(:)
integer,allocatable :: MOtype_back(:)
real*8,allocatable :: CO_back(:,:)
!Other property
integer :: iresinfo=0 !=0/1: Residue information is unavailable/available in this file

!-------- PBC information
!Translation vectors of the cell in Bohr. The vector of nonperiodic direction(s) is automatically set to a norm vector perpendicular to periodic directions in subroutine "init_PBC"
real*8 :: cellv1(3)=0,cellv2(3)=0,cellv3(3)=0
!The ones with _bk suffix are used to backup PBC information, the ones with _org are original system information
real*8 :: cellv1_org(3),cellv2_org(3),cellv3_org(3),cellv1_bk(3),cellv2_bk(3),cellv3_bk(3)
integer :: ifPBC=0,ifPBC_org,ifPBC_bk !Dimension of periodicity. 0=Isolated system, 1/2/3/=one/two/three dimensions
integer :: PBCnx,PBCny,PBCnz,ifdoPBCx,ifdoPBCy,ifdoPBCz !PBC setting actually used in calculation, they will be specified by init_PBC
integer :: PBCnx_in=1,PBCny_in=1,PBCnz_in=1,ifdoPBCx_in=1,ifdoPBCy_in=1,ifdoPBCz_in=1 !PBC setting read from settings.ini
integer :: iorbcomplex=1 !=1/2: Calculate real/imaginary part of orbital wavefunction in subroutine orbderv_PBC
!Neighbouring information of GTFs (including images) w.r.t. reduced grids
integer,allocatable :: neighGTF(:,:,:,:),neighnGTF(:,:,:) !neighGTF(1:neighnGTF(i,j,k),i,j,k) contains neighbouring GTF indices at reduced grid (i,j,k). Grid index range in X is 0 to size(neighnGTF,1)-1
integer,allocatable :: neighGTFcell(:,:,:,:,:) !neighGTF(1/2/3,m,i,j,k) is cell index in direction 1/2/3 of the m GTF neighbouring to reduced grid (i,j,k)
real*8 spcred !Spacing of reduced grids
real*8 orgx_neigh,orgy_neigh,orgz_neigh !Starting position of reduced grids

!-------- Connectivity matrix
!Loaded from .mol/mol2 using readmol/readmol2 or readmolconn (from mol), value is formal bond order; can also be guessed via genconnmat, value is 1/0 (connected, not connected)
!Special: ar (aromatic) in mol2 is load as 4, am (nitrogen in piptide bond) in mol2 is read as 1, "un = unknown", "nc = not connected" and "du = dummy" are read as 0
integer*2,allocatable :: connmat(:,:) !Diagonal terms are always zero

!-------- Energy related arrays and matrices
real*8,allocatable,target :: FmatA(:,:),FmatB(:,:) !Fock matrix of total/alpha and beta spin
!-------- Trajectory
integer :: nframetraj=0 !The number of frames in the trajectory
real*8,allocatable :: traj(:,:,:) !traj(1/2/3,a,i) corresponds to x/y/z of the ath atom in frame i
!-------- Points loaded from external file
integer :: numextpt=0
real*8,allocatable :: extpt(:,:),extpttmp(:) !extpt(i,1:4) corresponds to X/Y/Z/value of point i, length unit is Bohr. extpttmp only records function value
!-------- Atomic radial densities, may be loaded from .rad file, or generated by e.g. Hirshfeld-I and MBIS procedures
integer,allocatable :: atmradnpt(:) !How many radial points that each atom has
real*8 atmradpos(200) !Position of radial points, shared by all atoms, since it is generated by the same rule
real*8,allocatable :: atmraddens(:,:) !(j,iatm) corresponds to density value at j radial point of iatm atom
real*8,allocatable :: elemraddens(:,:) !(j,iele) corresponds to density value at j radial point of element iele
integer,allocatable :: elemraddens_npt(:) !Number of non-zero elements in elemraddens for each element
!---------- Used for passing data of spectrum plotting, as well as used by DOS
real*8,allocatable :: datax(:),str(:),FWHM(:) !Transition energy, strength and FWHM loaded from only one file
integer iROAtype !Type of ROA spectrum
integer iUVdir !Type of Directions of UV-Vis
real*8 UVdirvec(3) !Direction of oscillating electric field when iUVdir=7
integer :: istrtype=0 !Strength type. Avoid selecting the type every time when loading multiple files


!!!!!!!!!!!!!!!!!!!!!! Parameter !!!!!!!!!!!!!!!!!!!!!!
!For passing Dislin main parent GUI and draw widget identifier
integer idissetlight1,idissetlight2,idissetlight3,idissetlight4,idissetlight5,idissetlightall0,idissetlightall1,idissetangle,idissetzoom,idissetplaneXVU,idissetplaneYVU
integer idisgraph,idiszoomin,idiszoomout,idisisosurscl,idisscrval,idisshowbothsign,idisshowisosur,idisshowmol,idisisosursec,iorbseltext,iorbtxt,iorblis,idisshowdatarange
integer idisorbinfo2,idisorbinfo3,idissetVANG3D
integer idisshowatmlab,idisshowaxis,idisbondradius,idislabelsize,idisbondcrit,idisatmsize,idisshowpathlab !In draw mol GUI
integer idisshowattlab,idisdrawinternalbasin,idisattsize !Draw basin GUI
integer idisshow3n3,idisshow3n1,idisshow3p1,idisshow3p3,idisshowCPlab,idisshowpath,idisshowbassurf
integer idisshowlocminlab,idisshowlocmaxlab,idisshowlocminpos,idisshowlocmaxpos !For molecular surface analysis
!For setting isosurface style, colors
integer idisisosur1style,idisisosur1solid,idisisosur1mesh,idisisosur1point,idisisosur1solidmesh,idisisosur1tpr,idisisosur1opa
integer idisisosur2style,idisisosur2solid,idisisosur2mesh,idisisosur2point,idisisosur2solidmesh,idisisosur2tpr,idisisosur2opa
integer idisisosurallstyle,idisisosurallsolid,idisisosurallmesh,idisisosurallpoint,idisisosurallsolidmesh,idisisosuralltpr
!For setting box
integer idisboxsizeX,idisboxsizeY,idisboxsizeZ,idisboxposX,idisboxposY,idisboxposZ,idisboxspc,idisnpt
!=1: Show system and orbitals  =2: Show relief plane map  =3: Show isosurface  =4: Show system and CPs &
!=5: Show system and surface extrema  =6: Show basin or domain space  =7: Set box in GUI  =8: Minimum GUI
integer GUI_mode 

!Plotting external parameter, can be set in settings.ini
character(len=4) :: graphformat="png ",graphformatname(9)=(/ "png ","gif ","tiff","bmp ","ps  ","eps ","pdf ","wmf ","svg " /)
integer :: graph1Dwidth=1280,graph1Dheight=800,graph2Dwidth=1280,graph2Dheight=1200,graph3Dwidth=1400,graph3Dheight=1400
integer :: itickreverse=0,iticks=2,symbolsize=8,ilenunit1D=1,ilenunit2D=1,iatmlabtype=1,iatmlabtype3D=3,iplaneextdata=0,itransparent=0
integer :: numdigx=2,numdigy=2,numdigz=3,numdiglinex=3,numdigliney=3,numdigctr=3
real*8 :: planestpx=1.5D0,planestpy=1.5D0,planestpz=0.1D0
integer :: fillcoloritpx=5,fillcoloritpy=3,pleatmlabsize=50
real*8 :: disshowlabel=0.5D0,disshowCP=0.5D0,disshowpath=0.5D0
real*8 :: bondclrR=0.1D0,bondclrG=1D0,bondclrB=0.1D0,atmlabclrR=0D0,atmlabclrG=0D0,atmlabclrB=0D0
real*8 :: CPlabclrR=0D0,CPlabclrG=0D0,CPlabclrB=1D0 !CP label color
real*8 :: CP3n3RGB(3)=(/0.72D0,0D0,0.72D0/),CP3n1RGB(3)=(/1D0,0.5D0,0D0/),CP3p1RGB(3)=(/1D0,1D0,0D0/),CP3p3RGB(3)=(/0D0,1D0,0D0/)
real*8 :: CP3n3RGB_2D(3)=(/0.72D0,0D0,0.72D0/),CP3n1RGB_2D(3)=(/0D0,0D0,1D0/),CP3p1RGB_2D(3)=(/1D0,0.5D0,0D0/),CP3p3RGB_2D(3)=(/0D0,1D0,0D0/)
real*8 :: IBSclrR=0.1D0,IBSclrG=1D0,IBSclrB=0.1D0 !Interbasin surface color
real*8 :: atm3Dclr(0:nelesupp,3) !Colors of the atom spheres shown in 3D plots, set in "loadsetting" routine
character(len=200) :: ttfontfile=" "

!Plotting Internal parameter
integer :: imodlayout=0,plotwinsize3D=90,ishowhydrogen=1,ishoworbsel_prt=1,iplotwfndens=1
integer :: idrawbasinidx=0 !Draw which basin
integer :: idrawinternalbasin=0 !=1 Draw internal part of the basin, =0 Only draw boundary grids
integer :: ifixorbsign=0 !if 1, during generating orbital isosurface by drawmolgui, most part will always be positive (namely if sum(cubmat)<0 or sum(cubmattmp)<0, the data sign will be inverted)
integer :: iatom_on_plane=0,iatom_on_plane_far=0,ibond_on_plane=0,plesel,IGRAD_ARROW=0,ILABEL_ON_CONTOUR,ncontour
integer :: ictrlabsize=20,ivdwctrlabsize=0,iwidthvdwctr=10,iwidthposctr=3,iwidthnegctr=3,iwidthgradline=1,iclrindbndlab=14,plane_axistextsize=60,plane_axisnamesize=50,curve_axistextsize=45,curve_axisnamesize=45,bondthick2D=10
integer :: iclrindctrpos=5,iclrindctrneg=5,ivdwclrindctr=3,iclrindgradline=6,vdwctrstyle(2)=(/1,0/),ctrposstyle(2)=(/1,0/),ctrnegstyle(2)=(/10,15/)
integer :: isavepic=0,icurve_vertlinex=0,iclrindatmlab=1,imarkrefpos=0,ilog10y=0,iclrcurve=1
integer :: inucespplot=0,idrawmol=1,idrawisosur=0,isosursec=0,idrawtype=1,idrawcontour=1
integer :: iinvgradvec=0,icolorvecfield=0,vecclrind=30,idrawplanevdwctr=0,iplaneoutall=0,icurvethick=5,iclrtrans=0,ifillctrline=0,ishowclrfill_bar=0
integer :: ORIGIN_3D_X=0,ORIGIN_3D_Y=0
integer,allocatable :: highlightatomlist(:)
character :: stream_intmethod*5="RK2",curveYname*50=" "
character :: clrtransname(0:19)*50=(/ character(len=50) :: "Rainbow & white/black for out-of-limit data","Rainbow","Reversed rainbow","Rainbow starting from white","Spectrum","Reversed Spectrum",&
"Grey","Reversed Grey","Blue-White-Red","Red-White-Blue","Blue-Green-Red","Red-Green-Blue","White-Dark red","Black-Orange-Yellow","White-Dark green","Black-Green","White-Dark blue","Black-Blue-Cyan","Viridis","Yellow-Orange-Black" /)
real*8 :: drawisosurgui_SWGSTP=0.01D0,drawisosurgui_lowlim=-5,drawisosurgui_highlim=5 !Control isovalue bar setting in drawisosurgui
real*8 :: surcolorzmin,surcolorzmax !fillctr is the contour value will be draw on fillcolor map
real*8 :: curve_vertlinex=0D0,curvexyratio=0.618D0 !Gold partition
real*8 :: gradplotstep=0.002D0,gradplotdis=0.01D0,gradplottest=0.2D0,cutgradvec=0.3D0
real*8 :: clrRcub1same=0.3D0,clrGcub1same=0.75D0,clrBcub1same=0.3D0,clrRcub1oppo=0.3D0,clrGcub1oppo=0.45D0,clrBcub1oppo=0.9D0 !Color for isosurface 1 with solid style
real*8 :: clrRcub2same=0.4D0,clrGcub2same=0.5D0,clrBcub2same=0D0,clrRcub2oppo=0.35D0,clrGcub2oppo=0.1D0,clrBcub2oppo=0.9D0 !Color for isosurface 2 with solid style
real*8 :: clrRcub1samemeshpt=0.3D0,clrGcub1samemeshpt=0.75D0,clrBcub1samemeshpt=0.3D0,clrRcub1oppomeshpt=0.3D0,clrGcub1oppomeshpt=0.45D0,clrBcub1oppomeshpt=0.9D0 !Color for isosurface 1 with mesn style
real*8 :: clrRcub2samemeshpt=0.4D0,clrGcub2samemeshpt=0.5D0,clrBcub2samemeshpt=0D0,clrRcub2oppomeshpt=0.35D0,clrGcub2oppomeshpt=0.1D0,clrBcub2oppomeshpt=0.9D0 !Color for isosurface 2 with mesh style
real*8 :: opacitycub1=0.7D0,opacitycub2=0.7D0 !Opacity for isosurface 1 and 2 with transparent style
!About topology information on plane
integer :: imark3n3=1,imark3n1=1,imark3p1=1,imark3p3=1,imarkpath=1,sizemarkcp=30,sizemarkpath=10,sizemark3n1path=7,idrawintbasple=0,isurfstyle=2
integer :: iclrpath=14 !Color of topology paths (Brown)
integer :: iclr3n1path=15 !Color of interbasin paths (Dark blue)
integer,allocatable :: boldlinelist(:)
character(len=3) :: drawsurmesh="ON "
!Parameters for plotting extrema of a function on a contour line of another function
integer :: iextrema_on_contour=0,ifunc_2Dextrema
real*8 ctrval_2Dextrema
!Parameter for drawing molecular structure or 3D map
integer :: ienablelight1=1,ienablelight2=1,ienablelight3=1,ienablelight4=0,ienablelight5=0 !If enable lighting 1~5
integer :: ishowatmlab=1,ishowCPlab=0,ishowpathlab=0,ishowaxis=1,isosurshowboth=1,ishowdatarange=0,ishowcell=0,ishowboundaryatom=1,ishowboundarytopo=1,idrawpath=1,idrawbassurf=1,ishowattlab=0,ishowatt=0
integer :: isosur1style=1,isosur2style=1 !isosurface style,1/2/3/4/5=solid,mesh,points,solid+mesh,transparent
integer :: ishowlocminlab=0,ishowlocmaxlab=0,ishowlocminpos=0,ishowlocmaxpos=0 !For molecular surface analysis
integer :: ishow3n3=0,ishow3n1=0,ishow3p1=0,ishow3p3=0
real*8 :: bondcrit=1.15D0,textheigh=38D0,ratioatmsphere=1D0,ratioCPsphere=1D0,bondradius=0.2D0,attsphsize=0.1D0
real*8 :: XVU=150D0,YVU=30D0,ZVU=6D0,XFAC=1D0,VANG3DANG=28D0,XFOC=0D0,YFOC=0D0,ZFOC=0D0,camrotang=0D0 !3D view angle, they are all default value of DISLIN
integer :: iorthoview=0 !=0: Perspective, =1: Orthographic projection in 3D view
!Parameter for drawing domain defined by isosurfaces as grids
integer :: idrawdomainidx=0,idrawdomain=0,ndomain=0
real*8,allocatable :: gridxyz(:,:) !XYZ coordinate of grid that statisfied criterion
integer,allocatable :: domainsize(:) !The number of grids contained in each domain
integer,allocatable :: domaingrid(:,:) !The grid indices contained in each domain
!For passing ploting parameter from GUI routine to their call-back routine
!sur_value: The value of isosurface will be plotted by drawmol routine when idrawisosur=1 as well as drawisosurgui. sur_value_orb is specific for orbital isovalue
real*8 :: dp_init1,dp_end1,dp_init2,dp_end2,dp_init3,dp_end3,sur_value=0.05D0,sur_value_orb=0.05D0

!!! Other external parameter !!!
integer :: iautointgrid=1,radpot=75,sphpot=434 !sphpot=230/302/434/590/770, low is 50*434, high is 100*590
integer :: ispecial=0 !=0: Normal, =1 specific for Chunying Rong, =2 for Shubin's 2nd project
#ifdef _WIN32
integer :: isys=1 !Windows
#else
integer :: isys=2 !Linux/MacOS
#endif
integer :: igenP=1,iwfntmptype=1,iESPcode=2,outmedinfo=0,iaddprefix=0,intmolcust=0,isilent=0,idelvirorb=1
integer :: ifchprog=1,iloadascart=0,iloadGaugeom=1,iloadORCAgeom=1,maxloadexc=0,iprintLMOorder=0,iMCBOtype=0,ibasinlocmin=0
integer :: iuserfunc=0,iDFTxcsel=84,iKEDsel=0,ispheratm=1,ishowchgtrans=0,uESEinp=0,SpherIVgroup=0,MCvolmethod=2,readEDF=1,isupplyEDF=2,ishowptESP=1,imolsurparmode=1,nPGmaxatm=200
integer :: NICSnptlim=8000,iCDAcomp=1,ESPrhonlay=1
real*8 :: bndordthres=0.05D0,compthres=0.5D0,compthresCDA=1D0,expcutoff=-40D0,expcutoff_PBC=-20D0,ESPrhoiso=0D0
integer :: nthreads=4
integer*8 :: ompstacksize=200000000
character :: lastfile*200="",gaupath*200="",cubegenpath*200="",formchkpath*200="",orcapath*200="",orca_2mklpath*200="",dftd3path*200="",cubegendenstype*80="SCF"
!! About function calculation, external or internal parameters
integer :: RDG_addminimal=1,ELF_addminimal=1,num1Dpoints=3000,atomdenscut=1,nprevorbgrid=120000,paircorrtype=3,pairfunctype=1,srcfuncmode=1
integer :: ELFLOL_type=0,ipolarpara=0,iALIEdecomp=0,iskipnuc=0,ivdwprobe=6
integer :: nKEDmax=24
real*8 :: laplfac=1D0,uservar=0,uservar2=0,orbwei_delta=0.1D0,amIGMvdwscl=2D0
real*8 :: RDG_maxrho=0.05D0,RDGprodens_maxrho=0.1D0,IRI_rhocut=5D-5,aug1D=1.5D0,aug2D=4.5D0,aug3D=6D0,radcut=10D0,cfgcrossthres=0.01D0
real*8 :: refx=0D0,refy=0D0,refz=0D0
real*8 :: pleA=0D0,pleB=0D0,pleC=0D0,pleD=0D0 !!ABCD of the plane defined by main function 1000, used for special aims
real*8 :: globaltmp=0 !A variable can be used anywhere and can be set by option 5 of main function 1000, for debugging purpose avoiding re-compile code
!! About line/plane/grid calculation, inner parameter
!For 3D grid data. If the grid is not rectangle, only gridvec can fully define translation vectors
real*8 :: orgx,orgy,orgz,endx,endy,endz !Origin, end point and translation length in X/Y/Z. dx=0 means the box was not defined before
integer :: nx=80,ny=80,nz=80 !The number of grids in three directions (never necessarily in X,Y,Z!)
real*8 :: boxlenX,boxlenY,boxlenZ,boxcenX,boxcenY,boxcenZ !For temporary exchange data for setting box in GUI
real*8 :: gridv1(3),gridv2(3),gridv3(3) !1/2/3th translation vector of grid. dx,dy,dz corresponds to gridv1(1), gridv2(2), gridv3(3)
real*8 :: dx=0,dy,dz !Translation length in X/Y/Z. dx=0 means the box was not defined before. Using dx,dy,dz should be avoided in the future, always use gridv1/2/3 instead!
!For 2D plane map
real*8 :: v1x,v1y,v2x,v2y,v1z,v2z,a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,d1,d2 !Translation vector 1 and 2, three point in self-defined plane for projecting label, d1,d2=Length of v1,v2
real*8 :: orgx2D,orgy2D,orgz2D !X, Y, Z coordinate of origin of the plane map in molecular Cartesian space
integer :: ngridnum1=100,ngridnum2=100 !The number of points in two directions
!Specific for Shubin's project
real*8 :: steric_addminimal=1D-4,steric_potcutrho=0D0,steric_potcons=0D0
!Other
integer :: ifirstMultiwfn=1 !If 0, means we reload file via option -11 in main menu and don't need to do some initializations again

!Used for EDR(r;d) and D(r)
integer,parameter :: max_edr_exponents=50 !Maximum EDR exponents used to calculate EDR(r;d) and D(r)
real*8 :: dedr,edrastart,edrainc !Length scale to define EDR(r;d), start and increment in exponents to evaluate D(r) 
real*8 :: wrtexpo(max_edr_exponents) !For users write the number of EDR exponents that will be used in calculation
integer nedr !No of EDR exponents used to evaluate D(r)

end module


!-------- Module for topology analysis
module topo
integer :: ifunctopo=1 !Index of currently selected real space function for topology analysis
!CP related: (the ones with _tmp are used to temporarily show CPs including those at cell boundary)
integer,parameter :: maxnumcp=100000 !Maximum number of CPs
integer :: numcp=0 !Number of located CPs
integer CP_tmp_idx(maxnumcp) !Index of temporary CPs
real*8 CPpos(3,maxnumcp),CPpos_tmp(3,maxnumcp) !XYZ of CPs
integer :: CPtype(maxnumcp)=0,CPtype_tmp(maxnumcp) !Type of CPs. 0=unknown 1=(3,-3) 2=(3,-1) 3=(3,+1) 4=(3,+3)
character :: CPtyp2lab(0:4)*6=(/ "  ??  ","(3,-3)","(3,-1)","(3,+1)","(3,+3)" /)
real*8 :: CPstepscale=1D0,gradconv=1D-6,dispconv=1D-7,minicpdis=0.03D0,vdwsumcrit=1.5D0,singularcrit=5D-22,CPsearchlow=0D0,CPsearchhigh=0D0,topotrustrad=0D0
integer :: topomaxcyc=120,ishowsearchlevel=0,itopomethod=1
integer :: lab_oneCP=0 !If only allowing labelling one CP. =0: Labelling all, =x: Only allow show label for CP x
!Path related: (the ones with _tmp are used to temporarily show paths including those at cell boundary)
integer,parameter :: maxnumpath=10000 !Maximum number of paths
integer,parameter :: maxpathpt=1500 !Maximum number of points in each path
integer :: numpath=0 !Number of generated paths
integer path_tmp_idx(maxnumpath) !Index of temporary paths
integer pathnumpt(maxnumcp),pathnumpt_tmp(maxnumcp) !How many actual points in each path
real*8 :: topopath(3,maxpathpt,maxnumpath),topopath_tmp(3,maxpathpt,maxnumpath) !Coordinate of topology paths. {x,y,z}, index of points in each path, index of paths
real*8 :: discritpathfin=0.05D0,pathstepsize=0.03D0
integer :: maxpathpttry=451 !The actual upper limit of point during path generation
integer :: npathtry=30
!Basin surface related:
integer :: nsurfpt=100,nsurfpathpercp=40,numbassurf=0 !Number of points in each surface path, number of paths in each interbasin surface, total number of basin surfaces
integer :: cp2surf(100000)=0 !Convert total index of (3,-1) to surface index, if zero, means no surface corresponds to this CP
real*8 :: surfpathstpsiz=0.02D0 !Step size in interbasin surface path
real*8,allocatable :: bassurpath(:,:,:,:) !Store interbasin paths. {x,y,z}, indices of points, index of path, index of (3,-1) CP
!Interbasin path on plane map (Note: Two direction paths are counted as one path)
integer :: nple3n1path=0 !Already generated in-plane path from (3,-1)
integer :: n3n1plept=200 !Number of points in each direction of path 
integer :: cp2ple3n1path(10000)=0 !Convert total index of (3,-1) on the given plane to interbasin path index, if zero, means no path corresponds to this CP
real*8 :: ple3n1pathstpsiz=0.02D0
real*8,allocatable :: ple3n1path(:,:,:,:) !Store path derived form (3,-1) on given plane. {x,y,z}, indices of points, direction path (1 or 2), index of (3,-1) CP on plane
end module


!--------- Module for surface analysis
module surfvertex
use defvar
type triangtype
	integer idx(3) !Consists of which three surface vertices
	real*8 area
	real*8 value !Mapped function value at geometry center
end type
type surfcor2vtxtype
!if k=surfcor2vtxtype(i,q)%athcor means the two corner with surface corner index of i and k, interpolated to the surface vertex with index of %itpvtx. this information is stored in slot q
	integer athcor !Another corner
	integer itpvtx !interpolated to which surface vertex index
end type
type(surfcor2vtxtype),allocatable :: surfcor2vtx(:,:)
integer,allocatable :: surcor2vtxpos(:) !Will add new interpolation relationship to which slot of surfcor2vtx

type(triangtype),allocatable :: surtriang(:) !Record center of generated surface triangle
integer nsurtri !Temporary accummlated index of triangles in generating process
type(content),allocatable :: survtx(:) !Record x,y,z coordinate of surface vertex, with interpolated function value
integer,allocatable :: abs2suridx(:,:,:) ! Convert absolute indices of corners to surface vertex indices
integer,allocatable :: vtxconn(:,:) !(i,j)=k means the two surface vertices with index of i and k are connected, j is storage slot
integer,allocatable :: vtxconnpos(:) !Records current slot range of vtxconn
real*8 surfisoval,tetravol0,tetravol1,tetravol2,tetravol3 !volume of interpolated tetradrons of type 1,2,3
integer nsurvtx,nsurlocmin,nsurlocmax
integer surlocmaxidx(10000),surlocminidx(10000) !Store indices of local minimum and maximum points. If =0, means this slot is empty or has been discarded
integer :: nbisec=3 !Do how many times bisection before linear interpolation
integer :: ifuncintp=1 !Use which real space function to do bisection interpolation
integer,allocatable :: elimvtx(:),elimtri(:)
end module


!---------- Module for basin analysis
module basinintmod
integer ixlow,ixup,iylow,iyup,izlow,izup !Lower and upper limits of looping grids in basin analysis. For isolated system, 2:(nx/ny/nz)-1, for periodic system, 1:nx/ny/nz
integer vec26x(26),vec26y(26),vec26z(26)
real*8 len26(26)
real*8 :: valcritclus=0.005D0
integer :: ifuncbasin=0 !Which real space function is calculated as grid data to partition the basin
integer :: mergeattdist=5
real*8 :: basinsphsize=0 !Size of spheres for showing basins
integer :: ishowbasinmethod=1 !=1 Show entire basin in GUI, =2: Only show rho>0.001 region
integer,allocatable :: gridbas(:,:,:) !Each grid belongs to which basin(attractor). -2=Boundary grids, -1=Traveled to boundary grid, 0=Unassigned, x=basin index
integer :: numatt=0 !The number of pristine attractors after near-grid method
integer :: numrealatt=0 !The number of actual(real) attractors (the ones left after clustering)
integer,allocatable :: attgrid(:,:) !attgrid(1/2/3,i)=(ix/iy/iz) of the grid that the ith pristine attractor attributes to
real*8,allocatable :: attval(:) !Value of pristine attractors
real*8,allocatable :: attxyz(:,:) !attxyz(1:3,numatt). XYZ coordinate of pristine attractors 
integer,allocatable :: attconv(:) !Attractor conversion list. If attconv(i)=j, means pristine attractor i is belong to actual attractor j. -1 and 0 is also included
integer,allocatable :: nrealatthas(:) !nrealatthas(i)=m means actual attractor i has m pristine attractors
integer,allocatable :: realatttable(:,:) !realatttable(i,j)=k means the jth member of the ith actual attractor is pristine attractor k
real*8,allocatable :: realattval(:),realattxyz(:,:) !Value and XYZ coordinate of actual attractors. For those having multiple pristine attractors, these arrays record average value
logical,allocatable :: interbasgrid(:,:,:) !.true. means this is a boundary grid, else it is an internal grid
logical,allocatable :: grdposneg(:,:,:) !.true. means the value at this grid is positive, .false. means negative. Used in "basinana" and "generatebasin"
real*8,allocatable :: BOM(:,:,:),BOMb(:,:,:) !Basin overlap matrix of total/alpha and BOM of beta. The first two sizes are number of occupied orbitals, the last is numrealatt
real*8,allocatable :: corerhogrid(:,:,:) !Grid data of core electron density, has same grid setting as cubmat
integer :: ifcorerho=0 !=1 means corerhogrid has been generated and will be used in due time
end module


!-------- NAO related arrays and matrices
!For NAOset, NAOocc, NAOene, the second index is spin, 0/1/2=total/alpha/beta. For other NAO related arrays, they are independent of spin
module NAOmod
integer iopshNAO !0: Closed shell, 1: Open shell (total, alpha and beta are respectively analyzed). This variable is set during loadNAOinfo
integer ncenter_NAO !The number of centers in involved in NAO analysis, usually equals to ncenter
character(len=2),allocatable :: atmname_NAO(:) !Name of centers in involved in NAO analysis, e.g. C, H, O, usually equals to a%name
!NAO information
integer numNAO !The number of NAOs
integer,allocatable :: NAOinit(:),NAOend(:) !size of ncenter_NAO. Initial and ending indices of NAOs of atoms
integer,allocatable :: NAOcen(:) !size of numNAO. Center index of NAOs belong to
character,allocatable :: NAOcenname(:)*2 !size of numNAO. The element name loaded from NPA output, e.g. C, H, O
character,allocatable :: NAOset(:,:)*3  !(numNAO,0:2). Set of NAOs attributed to. i.e. Cor, Val, Ryd (Case sensitive!)
character,allocatable :: NAOshname(:)*3  !size of numNAO. Shell name that NAOs attributed to, e.g. 2s, 3p, 3d
character,allocatable :: NAOtype(:)*7 !size of numNAO. Type name of NAOs. e.g. px, dx2y2
real*8,allocatable :: NAOocc(:,:),NAOene(:,:) !size of numNAO. Occupation and energy (a.u.) of NAOs
!Shell information for NAOs
integer numNAOsh !Number of NAO shells
integer,allocatable :: bassh_NAO(:) !NAO basis index to shell index
integer,allocatable :: shcen_NAO(:) !NAO shell index to atom index
character,allocatable :: shname_NAO(:)*3 !NAO shell index to shell name (e.g. 2s, 3p, 3d)
character,allocatable :: shset_NAO(:,:)*3 !NAO shell index to shell set, i.e. Cor, Val, Ryd (Case sensitive!). 0/1/2=total/alpha/beta
!In the following matrices, "nbasis" and NBsUse are the number of basis functions before and after linear dependency elimination, respectively
!Note that AONAO does not distinguish spin, because NAO orbitals are always generated using total density matrix
real*8,allocatable :: DMNAO(:,:),DMNAOa(:,:),DMNAOb(:,:) !size of (numNAO,numNAO) (I found in rare case is (nbasis,nbasis)). Density matrix of total electrons, alpha electrons and beta electrons
real*8,allocatable :: NAOMO(:,:) !size of (numNAO,NBsUse). (i,r) is coeff. of NAO i in MO r. If numNAO<nbasis, the gap is filled by blank. For open shell, this records alpha part.
real*8,allocatable :: NAOMOb(:,:) !NAOMO for beta part
real*8,allocatable :: AONAO(:,:) !size of (nbasis,numNAO)
end module